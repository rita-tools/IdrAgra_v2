module mod_crop_yield
use mod_constants, only: dp
use mod_grid, only: grid_i
use mod_common, only: balance1_matrices, balance2_matrices
use mod_meteo, only: meteo_mat
use mod_crop_phenology, only: crop_pars_matrices, crop_matrices
use cli_save_outputs, only: yield_t
implicit none

contains

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
                        yld%f_HS_sum%mat(i,j,n) = yld%f_HS_sum%mat(i,j,n) + 1
                    else if (meteo%T_ave(i,j) >= pheno%T_crit(i,j) .and. meteo%T_ave(i,j) < pheno%T_lim(i,j)) then
                        yld%f_HS_sum%mat(i,j,n) = yld%f_HS_sum%mat(i,j,n) + 1 -                                                 &
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
                    yld%transp_ratio_sum%mat(i,j,n) = yld%transp_ratio_sum%mat(i,j,n) + (h_transp_pot) / meteo%et0(i,j)
                end if
            end if
        end do
    end do
end subroutine accumulate_daily_yield

subroutine calculate_annual_yield(yld, crop_map, domain)
    type(yield_t), intent(inout) :: yld
    type(crop_matrices), intent(in) :: crop_map
    type(grid_i), intent(in) :: domain
    integer :: i, j, z

    do j = 1, size(domain%mat, 2)
        do i = 1, size(domain%mat, 1)
            do z = 1, size(crop_map%TSP_high, 3)
                if (domain%mat(i,j) /= domain%header%nan) then

                    ! A declared slot can be absent from the daily crop series when one crop overwrites another
                    if (crop_map%ii0(i,j,z) == 0 .and. crop_map%iie(i,j,z) == 0) cycle

                    if ((crop_map%TSP_high(i,j,z) - crop_map%TSP_low(i,j,z)) /= 0.0_dp) then
                        yld%f_HS%mat(i,j,z) = yld%f_HS_sum%mat(i,j,z) / &
                          & (crop_map%TSP_high(i,j,z) - crop_map%TSP_low(i,j,z))
                    else
                        yld%f_HS%mat(i,j,z) = real(domain%header%nan)
                    end if

                    ! Calculate potential yield
                    yld%biomass_pot%mat(i,j,z) = crop_map%wp_adj(i,j,z) * yld%transp_ratio_sum%mat(i,j,z)
                    yld%yield_pot%mat(i,j,z) = yld%biomass_pot%mat(i,j,z) * crop_map%HI(i,j,z)

                    ! Calculate production reduction due to water stress
                    yld%f_WS%mat(i,j,z) = 1 - crop_map%Ky_tot(i,j,z) * &
                      & (1 - sum(yld%T_act_sum%mat(i,j,:,z)) / &
                      & sum(yld%T_pot_sum%mat(i,j,:,z)))

                    if (yld%f_WS%mat(i,j,z) < 0) yld%f_WS%mat(i,j,z) = 0

                    yld%f_WS_stage%mat(i,j,z) = (1 - crop_map%Ky_pheno(i,j,z,1) * &
                      & (1 - yld%T_act_sum%mat(i,j,1,z) / yld%T_pot_sum%mat(i,j,1,z))) &
                      & ** (yld%days_in_stage%mat(i,j,1,z) / sum(yld%days_in_stage%mat(i,j,:,z)))

                    yld%f_WS_stage%mat(i,j,z) = (1 - crop_map%Ky_pheno(i,j,z,2) * &
                      & (1 - yld%T_act_sum%mat(i,j,2,z) / yld%T_pot_sum%mat(i,j,2,z))) &
                      & ** (yld%days_in_stage%mat(i,j,2,z) / sum(yld%days_in_stage%mat(i,j,:,z))) * &
                      & yld%f_WS_stage%mat(i,j,z)

                    yld%f_WS_stage%mat(i,j,z) = (1 - crop_map%Ky_pheno(i,j,z,3) * &
                      & (1 - yld%T_act_sum%mat(i,j,3,z) / yld%T_pot_sum%mat(i,j,3,z))) &
                      & ** (yld%days_in_stage%mat(i,j,3,z) / sum(yld%days_in_stage%mat(i,j,:,z))) * &
                      & yld%f_WS_stage%mat(i,j,z)

                    yld%f_WS_stage%mat(i,j,z) = (1 - crop_map%Ky_pheno(i,j,z,4) * &
                      & (1 - yld%T_act_sum%mat(i,j,4,z) / yld%T_pot_sum%mat(i,j,4,z))) &
                      & ** (yld%days_in_stage%mat(i,j,4,z) / sum(yld%days_in_stage%mat(i,j,:,z))) * &
                      & yld%f_WS_stage%mat(i,j,z)

                    yld%yield_act%mat(i,j,z) = yld%yield_pot%mat(i,j,z) * &
                      & min(yld%f_WS%mat(i,j,z), yld%f_WS_stage%mat(i,j,z)) * &
                      & yld%f_HS%mat(i,j,z)
                end if
            end do
        end do
    end do
end subroutine calculate_annual_yield

end module mod_crop_yield
