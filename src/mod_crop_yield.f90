module mod_crop_yield
use mod_constants, only: dp
use mod_grid, only: grid_i
use mod_common, only: balance1_matrices, balance2_matrices
use mod_meteo, only: meteo_mat
use mod_crop_phenology, only: crop_pars_matrices, crop_matrices
implicit none

type yield_3d_t
    character(len=255) :: fn = ''
    real(dp), dimension(:,:,:), allocatable :: mat !(nrows,ncols,n_crops_in_year)
end type yield_3d_t

type yield_4d_t
    character(len=255) :: fn = ''
    real(dp), dimension(:,:,:,:), allocatable :: mat !(nrows,ncols,n_stages,n_crops_in_year)
end type yield_4d_t

type yield_t
    type(yield_3d_t) :: biomass_pot     ! Potential biomass [t/ha]
    type(yield_3d_t) :: yield_pot       ! Potential yield [t/ha]
    type(yield_3d_t) :: yield_act       ! Actual yield [t/ha]
    type(yield_3d_t) :: T_pot_ratio_sum ! Daily transpiration ratio (t_pot/et0) accumulator; used to estimate biomass_pot from water productivity.
    type(yield_3d_t) :: HS_sum          ! Daily Heat stress accumulator. Cool days add 1, stress days add <1. [-]
    type(yield_3d_t) :: f_HS            ! Mean heat-stress yield factor over the sensitive window, calculated as HS_sum/window_days [-]
    type(yield_3d_t) :: f_WS_tot        ! Water stress yield factor calculated from the whole crop cycle [-]
    type(yield_3d_t) :: f_WS_stage      ! Water stress yield factor calculated stage-by-stage [-]
    type(yield_4d_t) :: T_act_sum       ! Cumulative transpiration in each stage [mm]
    type(yield_4d_t) :: T_pot_sum       ! Cumulative potential transpiration in each stage [mm]
    type(yield_4d_t) :: days_in_stage   ! Number of days spent in each stage !TODO: should be an integer
end type yield_t

integer, parameter :: n_stages = 4

contains

subroutine initialize_yield(yld, domain, n_crops)
    type(yield_t), intent(inout) :: yld
    integer, dimension(:,:), intent(in) :: domain
    integer, intent(in) :: n_crops
    integer :: nrows, ncols

    nrows = size(domain, 1)
    ncols = size(domain, 2)
    allocate(yld%biomass_pot%mat    (nrows,ncols,n_crops), source=0.0_dp)
    allocate(yld%yield_pot%mat      (nrows,ncols,n_crops), source=0.0_dp)
    allocate(yld%yield_act%mat      (nrows,ncols,n_crops), source=0.0_dp)
    allocate(yld%f_WS_tot%mat       (nrows,ncols,n_crops), source=0.0_dp)
    allocate(yld%f_WS_stage%mat     (nrows,ncols,n_crops), source=0.0_dp)
    allocate(yld%f_HS%mat           (nrows,ncols,n_crops), source=0.0_dp)
    allocate(yld%HS_sum%mat         (nrows,ncols,n_crops), source=0.0_dp)
    allocate(yld%T_pot_ratio_sum%mat(nrows,ncols,n_crops), source=0.0_dp)
    allocate(yld%T_act_sum%mat      (nrows,ncols,n_stages,n_crops), source=0.0_dp)
    allocate(yld%T_pot_sum%mat      (nrows,ncols,n_stages,n_crops), source=0.0_dp)
    allocate(yld%days_in_stage%mat  (nrows,ncols,n_stages,n_crops), source=0.0_dp)
end subroutine initialize_yield

subroutine destroy_yield(yld)
    type(yield_t), intent(inout) :: yld

    if (allocated(yld%biomass_pot%mat))     deallocate(yld%biomass_pot%mat)
    if (allocated(yld%yield_pot%mat))       deallocate(yld%yield_pot%mat)
    if (allocated(yld%yield_act%mat))       deallocate(yld%yield_act%mat)
    if (allocated(yld%f_WS_tot%mat))        deallocate(yld%f_WS_tot%mat)
    if (allocated(yld%f_WS_stage%mat))      deallocate(yld%f_WS_stage%mat)
    if (allocated(yld%f_HS%mat))            deallocate(yld%f_HS%mat)
    if (allocated(yld%HS_sum%mat))          deallocate(yld%HS_sum%mat)
    if (allocated(yld%T_pot_ratio_sum%mat)) deallocate(yld%T_pot_ratio_sum%mat)
    if (allocated(yld%T_act_sum%mat))       deallocate(yld%T_act_sum%mat)
    if (allocated(yld%T_pot_sum%mat))       deallocate(yld%T_pot_sum%mat)
    if (allocated(yld%days_in_stage%mat))   deallocate(yld%days_in_stage%mat)
end subroutine destroy_yield

subroutine accumulate_daily_yield(yld, pheno, crop_map, meteo, wat_bal1, wat_bal2, domain, doy)
    type(yield_t), intent(inout) :: yld
    type(crop_pars_matrices), intent(inout) :: pheno
    type(crop_matrices), intent(in) :: crop_map
    type(meteo_mat), intent(in) :: meteo
    type(balance1_matrices), intent(in) :: wat_bal1
    type(balance2_matrices), intent(in) :: wat_bal2
    type(grid_i), intent(in) :: domain
    integer, intent(in) :: doy
    integer :: i, j, n, stage
    real(dp) :: h_transp_act, h_transp_pot

    ! TODO: add crop biomass from the previous year for winter cereals.

    do j = 1, size(domain%mat, 2)
        do i = 1, size(domain%mat, 1)
            if (domain%mat(i,j) /= domain%header%nan) then

                n = pheno%n_crop_in_year(i,j)
                stage = pheno%pheno_idx(i,j)
                h_transp_act = wat_bal1%h_transp_act(i,j) + wat_bal2%h_transp_act(i,j)
                h_transp_pot = wat_bal1%h_transp_pot(i,j) + wat_bal2%h_transp_pot(i,j)

                ! Accumulate thermal stress
                if (doy >= crop_map%TSP_low(i,j,n) .and. doy < crop_map%TSP_high(i,j,n)) then
                    if (meteo%T_ave(i,j) < pheno%T_crit(i,j)) then
                        yld%HS_sum%mat(i,j,n) = yld%HS_sum%mat(i,j,n) + 1
                    else if (meteo%T_ave(i,j) >= pheno%T_crit(i,j) .and. meteo%T_ave(i,j) < pheno%T_lim(i,j)) then
                        yld%HS_sum%mat(i,j,n) = yld%HS_sum%mat(i,j,n) + 1 -                                                   &
                                              & (meteo%T_ave(i,j) - pheno%T_crit(i,j)) / (pheno%T_lim(i,j) - pheno%T_crit(i,j))
                    end if
                end if

                ! Infer the four yield-development stages from the daily Kcb curve
                if (pheno%k_cb_low(i,j) == 0) then !%PS%, todo: replace this check with an explicit flag ("is_annual")
                    if (pheno%k_cb(i,j) == pheno%k_cb_low(i,j)) then
                        stage = 0
                    else if (pheno%k_cb(i,j) <= pheno%k_cb_mid(i,j) .and. (stage == 0 .or. stage == 1)) then
                        stage = 1 
                    else if (pheno%k_cb(i,j) < pheno%k_cb_high(i,j) .and. (stage == 1 .or. stage == 2)) then
                        stage = 2
                    else if (pheno%k_cb(i,j) == pheno%k_cb_high(i,j)) then
                        stage = 3
                    else
                        stage = 4
                    end if
                else
                    ! Permanent and pluriannual crops begin at stage one during vernalization and after harvest
                    if (pheno%k_cb(i,j) == pheno%k_cb_low(i,j)) then
                        stage = 1
                    else if (pheno%k_cb(i,j) < pheno%k_cb_high(i,j) .and. (stage == 1 .or. stage == 2)) then
                        stage = 2
                    else if (pheno%k_cb(i,j) == pheno%k_cb_high(i,j)) then
                        stage = 3
                    else
                        stage = 4
                    end if
                end if
                pheno%pheno_idx(i,j) = stage

                ! Accumulate water stress
                if (stage > 0) then
                    yld%T_act_sum%mat(i,j,stage,n) = yld%T_act_sum%mat(i,j,stage,n) + h_transp_act
                    yld%T_pot_sum%mat(i,j,stage,n) = yld%T_pot_sum%mat(i,j,stage,n) + h_transp_pot
                    yld%days_in_stage%mat(i,j,stage,n) = yld%days_in_stage%mat(i,j,stage,n) + 1
                end if

                ! Update transpiration ratio
                if (meteo%et0(i,j) > 0) then
                    yld%T_pot_ratio_sum%mat(i,j,n) = yld%T_pot_ratio_sum%mat(i,j,n) + (h_transp_pot) / meteo%et0(i,j)
                end if
            end if
        end do
    end do
end subroutine accumulate_daily_yield

subroutine calculate_annual_yield(yld, crop_map, domain)
    type(yield_t), intent(inout) :: yld
    type(crop_matrices), intent(in) :: crop_map
    type(grid_i), intent(in) :: domain
    integer :: i, j, z, heat_sensitive_window, stage
    real(dp) :: t_deficit, stage_duration_frac, ky, total_days

    do j = 1, size(domain%mat, 2)
        do i = 1, size(domain%mat, 1)
            do z = 1, size(crop_map%TSP_high, 3)
                if (domain%mat(i,j) /= domain%header%nan) then

                    ! A declared slot can be absent from the daily crop series when one crop overwrites another; skip it
                    if (crop_map%ii0(i,j,z) == 0 .and. crop_map%iie(i,j,z) == 0) cycle

                    ! Calculate potential yield
                    yld%biomass_pot%mat(i,j,z) = crop_map%wp_adj(i,j,z) * yld%T_pot_ratio_sum%mat(i,j,z)
                    yld%yield_pot%mat(i,j,z) = yld%biomass_pot%mat(i,j,z) * crop_map%HI(i,j,z)

                    ! Heat stress !%PS%, todo: due to how TSP_high and _low are calculated, this is not correct for crops crossing the year boundary
                    heat_sensitive_window = crop_map%TSP_high(i,j,z) - crop_map%TSP_low(i,j,z)
                    if ((heat_sensitive_window) /= 0._dp) then
                        yld%f_HS%mat(i,j,z) = yld%HS_sum%mat(i,j,z) / (heat_sensitive_window)
                    else
                        yld%f_HS%mat(i,j,z) = 1._dp !%PS%: no stess window; crop doesn't enter heat stress --> f_HS = 1
                    end if

                    ! "Total" water stress
                    yld%f_WS_tot%mat(i,j,z) = 1 - crop_map%Ky_tot(i,j,z)                                          * &
                                            & (1 - sum(yld%T_act_sum%mat(i,j,:,z)) / sum(yld%T_pot_sum%mat(i,j,:,z)))
                    yld%f_WS_tot%mat(i,j,z) = max(0._dp, yld%f_WS_tot%mat(i,j,z))

                    ! Stage-specific water stress
                    total_days = sum(yld%days_in_stage%mat(i,j,:,z)) !%PS%: note that total_days does not include days in stage 0
                    yld%f_WS_stage%mat(i,j,z) = 1._dp
                    do stage = 1, n_stages
                        ky = crop_map%Ky_pheno(i,j,z,stage)
                        t_deficit = 1 - yld%T_act_sum%mat(i,j,stage,z) / yld%T_pot_sum%mat(i,j,stage,z)
                        stage_duration_frac = (yld%days_in_stage%mat(i,j,stage,z) / total_days)
                        yld%f_WS_stage%mat(i,j,z) = yld%f_WS_stage%mat(i,j,z) * (1 - ky * t_deficit)**stage_duration_frac
                    end do

                    ! Actual yield depends on both heat and water stress
                    yld%yield_act%mat(i,j,z) = yld%yield_pot%mat(i,j,z) * &
                                             & min(yld%f_WS_tot%mat(i,j,z), yld%f_WS_stage%mat(i,j,z)) * yld%f_HS%mat(i,j,z)
                end if
            end do
        end do
    end do
end subroutine calculate_annual_yield

end module mod_crop_yield
