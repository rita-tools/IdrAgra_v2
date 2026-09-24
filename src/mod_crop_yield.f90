module mod_crop_yield
use mod_constants, only: dp
use mod_grid, only: grid_i
use mod_common, only: balance1_matrices, balance2_matrices
use mod_meteo, only: meteo_mat
use mod_crop_phenology, only: crop_pars_matrices, crop_matrices
use cli_save_outputs, only: yield_t
implicit none

contains

subroutine accumulate_daily_yield(yield, pheno, crop_map, meteo, wat_bal1, wat_bal2, domain, doy)
    type(yield_t), intent(inout) :: yield
    type(crop_pars_matrices), intent(inout) :: pheno
    type(crop_matrices), intent(in) :: crop_map
    type(meteo_mat), intent(in) :: meteo
    type(balance1_matrices), intent(in) :: wat_bal1
    type(balance2_matrices), intent(in) :: wat_bal2
    type(grid_i), intent(in) :: domain
    integer, intent(in) :: doy
    integer :: i, j

    ! TODO: add crop biomass from the previous year for winter cereals.

    do j = 1, size(domain%mat, 2)
        do i = 1, size(domain%mat, 1)
            if (domain%mat(i,j) /= domain%header%nan) then

                ! Accumulate thermal stress
                if (doy >= crop_map%TSP_low(i,j,pheno%n_crop_in_year(i,j)) .and. &
                  & doy < crop_map%TSP_high(i,j,pheno%n_crop_in_year(i,j))) then
                    if (meteo%T_ave(i,j) < pheno%T_crit(i,j)) then
                        yield%f_HS_sum%mat(i,j,pheno%n_crop_in_year(i,j)) = &
                          & yield%f_HS_sum%mat(i,j,pheno%n_crop_in_year(i,j)) + 1
                    else if (meteo%T_ave(i,j) >= pheno%T_crit(i,j) .and. &
                           & meteo%T_ave(i,j) < pheno%T_lim(i,j)) then
                        yield%f_HS_sum%mat(i,j,pheno%n_crop_in_year(i,j)) = &
                          & yield%f_HS_sum%mat(i,j,pheno%n_crop_in_year(i,j)) + 1 - &
                          & (meteo%T_ave(i,j) - pheno%T_crit(i,j)) / &
                          & (pheno%T_lim(i,j) - pheno%T_crit(i,j))
                    end if
                end if

                ! Infer the four yield-development stages from the daily Kcb curve
                if (pheno%k_cb_low(i,j) == 0) then
                    if (pheno%k_cb(i,j) == pheno%k_cb_low(i,j)) then
                        pheno%pheno_idx(i,j) = 0
                    else if (pheno%k_cb(i,j) <= pheno%k_cb_mid(i,j) .and. &
                           & (pheno%pheno_idx(i,j) == 0 .or. pheno%pheno_idx(i,j) == 1)) then
                        pheno%pheno_idx(i,j) = 1
                    else if (pheno%k_cb(i,j) < pheno%k_cb_high(i,j) .and. &
                           & (pheno%pheno_idx(i,j) == 1 .or. pheno%pheno_idx(i,j) == 2)) then
                        pheno%pheno_idx(i,j) = 2
                    else if (pheno%k_cb(i,j) == pheno%k_cb_high(i,j)) then
                        pheno%pheno_idx(i,j) = 3
                    else
                        pheno%pheno_idx(i,j) = 4
                    end if
                else
                    ! Permanent and pluriannual crops begin at stage one during vernalization and after harvest
                    if (pheno%k_cb(i,j) == pheno%k_cb_low(i,j)) then
                        pheno%pheno_idx(i,j) = 1
                    else if (pheno%k_cb(i,j) < pheno%k_cb_high(i,j) .and. &
                           & (pheno%pheno_idx(i,j) == 1 .or. pheno%pheno_idx(i,j) == 2)) then
                        pheno%pheno_idx(i,j) = 2
                    else if (pheno%k_cb(i,j) == pheno%k_cb_high(i,j)) then
                        pheno%pheno_idx(i,j) = 3
                    else
                        pheno%pheno_idx(i,j) = 4
                    end if
                end if

                ! Accumulate water stress
                if (pheno%pheno_idx(i,j) > 0) then
                    yield%T_act_sum%mat(i,j,pheno%pheno_idx(i,j),pheno%n_crop_in_year(i,j)) = &
                      & yield%T_act_sum%mat(i,j,pheno%pheno_idx(i,j),pheno%n_crop_in_year(i,j)) + &
                      & wat_bal1%h_transp_act(i,j) + wat_bal2%h_transp_act(i,j)
                    yield%T_pot_sum%mat(i,j,pheno%pheno_idx(i,j),pheno%n_crop_in_year(i,j)) = &
                      & yield%T_pot_sum%mat(i,j,pheno%pheno_idx(i,j),pheno%n_crop_in_year(i,j)) + &
                      & wat_bal1%h_transp_pot(i,j) + wat_bal2%h_transp_pot(i,j)
                    yield%dev_stage%mat(i,j,pheno%pheno_idx(i,j),pheno%n_crop_in_year(i,j)) = &
                      & yield%dev_stage%mat(i,j,pheno%pheno_idx(i,j),pheno%n_crop_in_year(i,j)) + 1
                end if

                ! Update transpiration ratio
                if (meteo%et0(i,j) > 0) then
                    yield%transp_ratio_sum%mat(i,j,pheno%n_crop_in_year(i,j)) = &
                      & yield%transp_ratio_sum%mat(i,j,pheno%n_crop_in_year(i,j)) + &
                      & (wat_bal1%h_transp_pot(i,j) + wat_bal2%h_transp_pot(i,j)) / meteo%et0(i,j)
                end if
            end if
        end do
    end do
end subroutine accumulate_daily_yield

subroutine calculate_annual_yield(yield, crop_map, domain)
    type(yield_t), intent(inout) :: yield
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
                        yield%f_HS%mat(i,j,z) = yield%f_HS_sum%mat(i,j,z) / &
                          & (crop_map%TSP_high(i,j,z) - crop_map%TSP_low(i,j,z))
                    else
                        yield%f_HS%mat(i,j,z) = real(domain%header%nan)
                    end if

                    ! Calculate potential yield
                    yield%biomass_pot%mat(i,j,z) = crop_map%wp_adj(i,j,z) * yield%transp_ratio_sum%mat(i,j,z)
                    yield%yield_pot%mat(i,j,z) = yield%biomass_pot%mat(i,j,z) * crop_map%HI(i,j,z)

                    ! Calculate production reduction due to water stress
                    yield%f_WS%mat(i,j,z) = 1 - crop_map%Ky_tot(i,j,z) * &
                      & (1 - sum(yield%T_act_sum%mat(i,j,:,z)) / &
                      & sum(yield%T_pot_sum%mat(i,j,:,z)))

                    if (yield%f_WS%mat(i,j,z) < 0) yield%f_WS%mat(i,j,z) = 0

                    yield%f_WS_stage%mat(i,j,z) = (1 - crop_map%Ky_pheno(i,j,z,1) * &
                      & (1 - yield%T_act_sum%mat(i,j,1,z) / yield%T_pot_sum%mat(i,j,1,z))) &
                      & ** (yield%dev_stage%mat(i,j,1,z) / sum(yield%dev_stage%mat(i,j,:,z)))

                    yield%f_WS_stage%mat(i,j,z) = (1 - crop_map%Ky_pheno(i,j,z,2) * &
                      & (1 - yield%T_act_sum%mat(i,j,2,z) / yield%T_pot_sum%mat(i,j,2,z))) &
                      & ** (yield%dev_stage%mat(i,j,2,z) / sum(yield%dev_stage%mat(i,j,:,z))) * &
                      & yield%f_WS_stage%mat(i,j,z)

                    yield%f_WS_stage%mat(i,j,z) = (1 - crop_map%Ky_pheno(i,j,z,3) * &
                      & (1 - yield%T_act_sum%mat(i,j,3,z) / yield%T_pot_sum%mat(i,j,3,z))) &
                      & ** (yield%dev_stage%mat(i,j,3,z) / sum(yield%dev_stage%mat(i,j,:,z))) * &
                      & yield%f_WS_stage%mat(i,j,z)

                    yield%f_WS_stage%mat(i,j,z) = (1 - crop_map%Ky_pheno(i,j,z,4) * &
                      & (1 - yield%T_act_sum%mat(i,j,4,z) / yield%T_pot_sum%mat(i,j,4,z))) &
                      & ** (yield%dev_stage%mat(i,j,4,z) / sum(yield%dev_stage%mat(i,j,:,z))) * &
                      & yield%f_WS_stage%mat(i,j,z)

                    yield%yield_act%mat(i,j,z) = yield%yield_pot%mat(i,j,z) * &
                      & min(yield%f_WS%mat(i,j,z), yield%f_WS_stage%mat(i,j,z)) * &
                      & yield%f_HS%mat(i,j,z)
                end if
            end do
        end do
    end do
end subroutine calculate_annual_yield

end module mod_crop_yield
