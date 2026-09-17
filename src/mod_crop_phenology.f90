module mod_crop_phenology
use mod_grid
use mod_utility, only: round_2darray
implicit none
! Daily spatial crop properties consumed by the existing water balance.
type crop_pars_matrices
    ! store the crop parameters for each calculation cells
    ! Biological state is stored separately in mod_daily_phenology.
    real(dp),dimension(:,:),pointer::k_cb
    real(dp),dimension(:,:),pointer::h
    real(dp),dimension(:,:),pointer::d_r
    real(dp),dimension(:,:),pointer::lai
    integer,dimension(:,:),pointer::cn_day
    real(dp),dimension(:,:),pointer::f_c
    integer,dimension(:,:),pointer::irrigation_class
    integer,dimension(:,:),pointer::cn_class
    real(dp),dimension(:,:),pointer::p
    real(dp),dimension(:,:),pointer::a
    real(dp),dimension(:,:),pointer::d_t_max
    real(dp),dimension(:,:),pointer::RF_t_max
    real(dp),dimension(:,:),pointer::RF_e
    real(dp),dimension(:,:),pointer::RF_t
    real(dp),dimension(:,:),pointer::T_lim
    real(dp),dimension(:,:),pointer::T_crit
    real(dp),dimension(:,:),pointer::HI
    real(dp),dimension(:,:),pointer::Ky_tot
    real(dp),dimension(:,:,:),pointer::Ky_pheno
    real(dp),dimension(:,:),pointer::k_cb_low
    real(dp),dimension(:,:),pointer::k_cb_mid
    real(dp),dimension(:,:),pointer::k_cb_high
    real(dp),dimension(:,:),pointer::wp_adj
    real(dp),dimension(:,:),pointer::p_day
    real(dp),dimension(:,:),pointer::k_cb_old           ! k_cb of previous day
    integer,dimension(:,:),pointer::n_crop_in_year
    integer,dimension(:,:),pointer::pheno_idx           ! phenological stage index
    real(dp),dimension(:,:),pointer::r_stress           ! plant resistance to (water) stress
end type crop_pars_matrices
contains
subroutine calculate_RF_t(d_t, crop_par_mat,domain)
    ! calculate root fraction in both evaporative and transpirative layer
    real(dp), dimension(:,:), intent(in)::d_t
    type(grid_i),intent(in)::domain
    type(crop_pars_matrices),intent(inout)::crop_par_mat

    ! populate pheno%RF_t & pheno%RF_e
    where (domain%mat /= domain%header%nan .and. &
           crop_par_mat%d_t_max > 0.0D0 .and. &
           crop_par_mat%RF_t_max >= 0.0D0 .and. crop_par_mat%RF_t_max <= 1.0D0)
        crop_par_mat%RF_t = d_t/crop_par_mat%d_t_max
        where (crop_par_mat%RF_t_max*crop_par_mat%RF_t + &
               (1.0D0-crop_par_mat%RF_t_max) > tiny(1.0D0))
            crop_par_mat%RF_t = crop_par_mat%RF_t_max*crop_par_mat%RF_t / &
                (crop_par_mat%RF_t_max*crop_par_mat%RF_t + (1.0D0-crop_par_mat%RF_t_max))
        elsewhere
            crop_par_mat%RF_t = 0.0D0
        end where
    elsewhere
        crop_par_mat%RF_t = 0.0D0
    end where

    crop_par_mat%RF_t = round_2darray(crop_par_mat%RF_t,6)
    crop_par_mat%RF_e = 1.0D0-crop_par_mat%RF_t

    where (domain%mat == domain%header%nan)
        crop_par_mat%RF_t = dble(domain%header%nan)
        crop_par_mat%RF_e = dble(domain%header%nan)
    end where

end subroutine calculate_RF_t


end module mod_crop_phenology
