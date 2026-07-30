program test_phenology_events
use mod_constants, only: dp
use mod_grid, only: grid_i
use mod_crop_phenology, only: crop_pheno_info, crop_matrices, crop_pars_matrices, populate_crop_cell
use mod_phenology_events, only: detect_reference_regrowth, make_cut_irrigation_halt_mask
implicit none

integer, parameter :: n_days = 20
type(crop_pheno_info), dimension(1) :: info_pheno
type(crop_matrices) :: crop_map
type(crop_pars_matrices) :: crop_pars
type(grid_i) :: soil_use
integer, dimension(1,1) :: irandom, ws_idx
logical, dimension(1,1) :: halt_mask
integer :: day, first_cut, regrowth_start, regrowth_end

call allocate_test_data()

! Normal calendars retain the random shift: local day 5 maps to reference day 2.
crop_map%suppress_irandom = .false.
call populate_crop_cell(1,1,crop_pars,info_pheno,irandom,5,ws_idx,soil_use,1,n_days,crop_map)
call assert_reference_day(2, 'normal irandom shift')

! Either an external sowing or harvest event suppresses irandom for the cycle.
crop_map%external_calendar = .true.
crop_map%suppress_irandom = .true.
call populate_crop_cell(1,1,crop_pars,info_pheno,irandom,2,ws_idx,soil_use,1,n_days,crop_map)
call assert_reference_day(2, 'external calendar suppresses irandom')

! The reference first cut is day 8 and its reusable regrowth interval is 8:11.
! External cuts are imposed on days 10 and 15.
crop_map%n_external_cuts = 2
crop_map%external_cut_doy(1,1,1,1:2) = [10,15]
crop_map%ref_first_cut = 8
crop_map%ref_regrowth_start = 8
crop_map%ref_regrowth_end = 11

call populate_crop_cell(1,1,crop_pars,info_pheno,irandom,9,ws_idx,soil_use,1,n_days,crop_map)
call assert_reference_day(7, 'reference cut is delayed until first external cut')

call populate_crop_cell(1,1,crop_pars,info_pheno,irandom,10,ws_idx,soil_use,1,n_days,crop_map)
call assert_reference_day(8, 'external cut resets to reference regrowth start')

call populate_crop_cell(1,1,crop_pars,info_pheno,irandom,14,ws_idx,soil_use,1,n_days,crop_map)
call assert_reference_day(11, 'late next cut holds mature regrowth day')

call populate_crop_cell(1,1,crop_pars,info_pheno,irandom,15,ws_idx,soil_use,1,n_days,crop_map)
call assert_reference_day(8, 'second external cut repeats the same reset')

! The heuristic must identify two mature-canopy drops and expose the complete
! interval beginning with the first post-cut day.
info_pheno(1)%k_cb%tab(:,1) = 0.20D0
info_pheno(1)%k_cb%tab(4:5,1) = 1.00D0
info_pheno(1)%k_cb%tab(9:10,1) = 1.00D0
call detect_reference_regrowth(info_pheno(1),1,1,n_days,first_cut,regrowth_start,regrowth_end)
if (first_cut /= 6 .or. regrowth_start /= 6 .or. regrowth_end /= 10) then
    print *, 'FAILED reference regrowth detection; first/start/end: ', &
        & first_cut, regrowth_start, regrowth_end
    error stop 1
end if

! Irrigation halts include ordinary Kcb-inferred cuts.
crop_map%n_external_cuts = 0
call make_cut_irrigation_halt_mask(crop_map,info_pheno,ws_idx,soil_use,irandom,4,n_days,2,halt_mask)
if (.not. halt_mask(1,1)) error stop 'FAILED inferred-cut irrigation halt'
call make_cut_irrigation_halt_mask(crop_map,info_pheno,ws_idx,soil_use,irandom,3,n_days,2,halt_mask)
if (halt_mask(1,1)) error stop 'FAILED inferred-cut halt extent'

! External cut dates supersede inferred dates for that crop cycle.
crop_map%n_external_cuts = 2
call make_cut_irrigation_halt_mask(crop_map,info_pheno,ws_idx,soil_use,irandom,8,n_days,2,halt_mask)
if (.not. halt_mask(1,1)) error stop 'FAILED external-cut irrigation halt'
call make_cut_irrigation_halt_mask(crop_map,info_pheno,ws_idx,soil_use,irandom,6,n_days,0,halt_mask)
if (halt_mask(1,1)) error stop 'FAILED disabled irrigation halt'

print *, 'phenology event regression tests passed'

contains

subroutine assert_reference_day(expected, label)
    integer, intent(in) :: expected
    character(len=*), intent(in) :: label
    integer :: actual
    actual = nint(100.0D0*crop_pars%k_cb(1,1))
    if (actual /= expected) then
        print *, 'FAILED: ', trim(label), '; expected/actual reference day: ', expected, actual
        error stop 1
    end if
    actual = nint(crop_pars%h(1,1))
    if (actual /= expected) then
        print *, 'FAILED full-vector rewind: ', trim(label), '; expected/actual H day: ', expected, actual
        error stop 1
    end if
end subroutine assert_reference_day

subroutine allocate_test_data()
    allocate(soil_use%mat(1,1))
    soil_use%header%nan = -9999
    soil_use%mat = 1
    irandom = 3
    ws_idx = 1

    allocate(info_pheno(1)%cycle_crop_slot(1,1))
    allocate(info_pheno(1)%crop_slot%tab(n_days,1))
    allocate(info_pheno(1)%k_cb%tab(n_days,1))
    allocate(info_pheno(1)%h%tab(n_days,1))
    allocate(info_pheno(1)%z_r%tab(n_days,1))
    allocate(info_pheno(1)%lai%tab(n_days,1))
    allocate(info_pheno(1)%cn_day%tab(n_days,1))
    allocate(info_pheno(1)%f_c%tab(n_days,1))
    allocate(info_pheno(1)%r_stress%tab(n_days,1))
    info_pheno(1)%cycle_crop_slot = 1
    info_pheno(1)%crop_slot%tab = 0
    info_pheno(1)%crop_slot%tab(2:18,1) = 1
    do day=1,n_days
        info_pheno(1)%k_cb%tab(day,1) = dble(day)/100.0D0
        info_pheno(1)%h%tab(day,1) = dble(day)
        info_pheno(1)%z_r%tab(day,1) = dble(day)
        info_pheno(1)%lai%tab(day,1) = dble(day)
        info_pheno(1)%cn_day%tab(day,1) = day
        info_pheno(1)%f_c%tab(day,1) = dble(day)
        info_pheno(1)%r_stress%tab(day,1) = dble(day)
    end do

    allocate(info_pheno(1)%irrigation_class(1,1), info_pheno(1)%cn_class(1,1))
    allocate(info_pheno(1)%p_raw_const(1,1), info_pheno(1)%a(1,1))
    allocate(info_pheno(1)%d_r_max(1,1), info_pheno(1)%max_RF_t(1,1))
    allocate(info_pheno(1)%T_lim(1,1), info_pheno(1)%T_crit(1,1))
    allocate(info_pheno(1)%HI(1,1), info_pheno(1)%Ky_tot(1,1))
    allocate(info_pheno(1)%Ky_pheno(1,1,4), info_pheno(1)%wp_adj(1,1,1))
    allocate(info_pheno(1)%kcb_phases%low(1,1), info_pheno(1)%kcb_phases%mid(1,1))
    allocate(info_pheno(1)%kcb_phases%high(1,1))
    info_pheno(1)%irrigation_class = 1
    info_pheno(1)%cn_class = 1
    info_pheno(1)%p_raw_const = 0.5D0
    info_pheno(1)%a = 1.0D0
    info_pheno(1)%d_r_max = 1.0D0
    info_pheno(1)%max_RF_t = 1.0D0
    info_pheno(1)%T_lim = 40.0D0
    info_pheno(1)%T_crit = 35.0D0
    info_pheno(1)%HI = 0.5D0
    info_pheno(1)%Ky_tot = 1.0D0
    info_pheno(1)%Ky_pheno = 1.0D0
    info_pheno(1)%wp_adj = 1.0D0
    info_pheno(1)%kcb_phases%low = 0.02D0
    info_pheno(1)%kcb_phases%mid = 0.10D0
    info_pheno(1)%kcb_phases%high = 0.18D0

    allocate(crop_map%ii0(1,1,1), crop_map%iie(1,1,1))
    allocate(crop_map%ii0_ref(1,1,1), crop_map%iie_ref(1,1,1))
    allocate(crop_map%k_cb_min(1,1,1), crop_map%k_cb_max(1,1,1))
    allocate(crop_map%external_calendar(1,1,1))
    allocate(crop_map%suppress_irandom(1,1))
    allocate(crop_map%n_external_cuts(1,1,1), crop_map%external_cut_doy(1,1,1,2))
    allocate(crop_map%ref_first_cut(1,1,1), crop_map%ref_regrowth_start(1,1,1))
    allocate(crop_map%ref_regrowth_end(1,1,1))
    crop_map%ii0 = 2
    crop_map%iie = 18
    crop_map%ii0_ref = 2
    crop_map%iie_ref = 18
    crop_map%k_cb_min = 0.20D0
    crop_map%k_cb_max = 1.00D0
    crop_map%external_calendar = .false.
    crop_map%suppress_irandom = .false.
    crop_map%n_external_cuts = 0
    crop_map%external_cut_doy = 0
    crop_map%ref_first_cut = 0
    crop_map%ref_regrowth_start = 0
    crop_map%ref_regrowth_end = 0

    allocate(crop_pars%k_cb(1,1), crop_pars%h(1,1), crop_pars%d_r(1,1))
    allocate(crop_pars%lai(1,1), crop_pars%cn_day(1,1), crop_pars%f_c(1,1))
    allocate(crop_pars%r_stress(1,1), crop_pars%irrigation_class(1,1))
    allocate(crop_pars%cn_class(1,1), crop_pars%p(1,1), crop_pars%a(1,1))
    allocate(crop_pars%d_t_max(1,1), crop_pars%RF_t_max(1,1))
    allocate(crop_pars%T_lim(1,1), crop_pars%T_crit(1,1), crop_pars%HI(1,1))
    allocate(crop_pars%Ky_tot(1,1), crop_pars%Ky_pheno(1,1,4))
    allocate(crop_pars%k_cb_low(1,1), crop_pars%k_cb_mid(1,1), crop_pars%k_cb_high(1,1))
    allocate(crop_pars%wp_adj(1,1), crop_pars%n_crop_in_year(1,1))
    allocate(crop_pars%pheno_idx(1,1))
    crop_pars%n_crop_in_year = 1
    crop_pars%pheno_idx = 1
end subroutine allocate_test_data

end program test_phenology_events
