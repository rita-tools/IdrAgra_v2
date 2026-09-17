program test_crop_kernel
use mod_constants, only: dp
use mod_utility, only: calc_doy
use cli_crop_parameters, only: crop_definition, interpolate_crop
use mod_daily_phenology, only: thermal_units, advance_thermal, window_dates
implicit none
type(crop_definition) :: c
real(dp) :: g,v
integer :: a,b
call near(thermal_units(20.0_dp,20.0_dp,10.0_dp,30.0_dp),10.0_dp,'constant temperature')
call near(thermal_units(5.0_dp,0.0_dp,10.0_dp,30.0_dp),0.0_dp,'below base')
call near(thermal_units(40.0_dp,35.0_dp,10.0_dp,30.0_dp),20.0_dp,'above cutoff')
call near(thermal_units(20.0_dp,0.0_dp,10.0_dp,30.0_dp),10.0_dp/acos(-1.0_dp),'sine threshold crossing')
allocate(c%gdd(2),c%values(2,8))
c%gdd=[100.0_dp,300.0_dp]; c%values=0; c%values(:,1)=[0.3_dp,0.9_dp]
c%values(:,5)=[1.0_dp,2.0_dp]
call near(interpolate_crop(c,200.0_dp,1),0.6_dp,'thermal interpolation')
call near(interpolate_crop(c,50.0_dp,1),0.3_dp,'lower endpoint')
call near(interpolate_crop(c,400.0_dp,1),0.9_dp,'upper endpoint')
call near(interpolate_crop(c,200.0_dp,5),1.0_dp,'discrete CN')
c%tbase=10; c%tcut=30; c%photo=1; c%dl_if=8; c%dl_ins=16
g=0; v=0
call advance_thermal(c,20.0_dp,20.0_dp,20.0_dp,12.0_dp,v,g)
call near(g,5.0_dp,'photoperiod accumulation')
c%photo=0; c%vern=.true.; c%tvmin=10; c%tvmax=25; c%vstart=0; c%vend=10; c%vfmin=0
g=0; v=0
call advance_thermal(c,20.0_dp,20.0_dp,20.0_dp,12.0_dp,v,g)
call near(v,1.0_dp,'vernalization state'); call near(g,1.0_dp,'vernalization factor')
c%sow_min=360; c%sow_delay=20
call window_dates(c,calc_doy(5,1,2022),0,a,b)
if(a/=calc_doy(26,12,2021).or.b/=calc_doy(15,1,2022)) error stop 'cross-year sowing window'
c%sow_min=60; c%sow_delay=0
call window_dates(c,calc_doy(1,1,2024),0,a,b)
if(a/=calc_doy(29,2,2024)) error stop 'leap-year window'
print *, 'Crop kernel checks passed'
contains
subroutine near(actual,expected,name)
    real(dp), intent(in) :: actual,expected
    character(len=*), intent(in) :: name
    if(abs(actual-expected)>1.e-10_dp) then
        print *, name,actual,expected
        error stop 1
    end if
end subroutine
end program
