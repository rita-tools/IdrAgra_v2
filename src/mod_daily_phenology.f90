! Daily array-based phenology. GDD/vernalization/productivity equations adapted
! from cropcoeff, E. A. Chiaradia, GPL-2.0-or-later. No generated crop tables.
module mod_daily_phenology
use mod_constants, only: dp
use mod_utility, only: calc_doy, calc_date
use mod_parameters, only: simulation
use mod_grid, only: grid_i
use mod_meteo, only: meteo_info, meteo_mat
use mod_crop_phenology, only: crop_pars_matrices
use cli_crop_parameters
implicit none
private
public :: begin_crop_run, advance_crops, finish_crop_day, end_crop_run, thermal_units, advance_thermal
public :: crop_state, window_dates
real(dp), parameter :: pi=acos(-1.0_dp)
type crop_state_arrays
    integer, allocatable :: lu(:,:), slot(:,:), occurrence(:,:), cut(:,:), next_slot(:,:)
    integer, allocatable :: sow(:,:), harvest_limit(:,:), last_harvest(:,:), next_sow(:,:), deadline(:,:)
    integer, allocatable :: offset(:,:), pending_offset(:,:), days(:,:,:), hs_days(:,:)
    real(dp), allocatable :: gdd(:,:), vern(:,:), corr(:,:,:), actual(:,:,:), potential(:,:,:), biomass(:,:), hs(:,:)
    logical, allocatable :: terminate(:,:), carried(:,:)
end type
type(crop_state_arrays), save :: crop_state
! Weather history is shared by stations, never duplicated per cell.
real(dp), allocatable, save :: weather(:,:,:) ! day, station, Tmax/Tmin/RHmin/wind
integer, save :: weather_start, weather_end, event_unit, harvest_unit, last_date=0, run_seed=1
contains
subroutine begin_crop_run(sim,stations,domain,warmup)
    type(simulation), intent(in) :: sim
    type(meteo_info), intent(in) :: stations(:)
    type(grid_i), intent(in) :: domain
    logical, intent(in) :: warmup
    integer :: nx,ny,s,u,k,ios,n,delta
    real(dp) :: row(7)
    character(len=512) :: path,prefix
    if(.not.allocated(weather)) then
        weather_start=stations(1)%start%doy; weather_end=stations(1)%finish%doy
        n=weather_end-weather_start+1
        allocate(weather(n,size(stations),4))
        do s=1,size(stations)
            if(stations(s)%start%doy/=weather_start.or.stations(s)%finish%doy/=weather_end) &
                call crop_fail('station dates differ')
            inquire(unit=stations(s)%unit,name=path)
            open(newunit=u,file=trim(path),status='old',action='read',iostat=ios)
            if(ios/=0) call crop_fail('opening weather cache '//trim(path))
            do k=1,4
                read(u,*,iostat=ios)
                if(ios/=0) call crop_fail('weather header '//trim(path))
            end do
            do k=1,n
                read(u,*,iostat=ios) row
                if(ios/=0) call crop_fail('weather record '//trim(path))
                if(row(1)<row(2)) call crop_fail('Tmax below Tmin in '//trim(path))
                weather(k,s,:)=[row(1),row(2),row(5),row(6)]
            end do
            close(u)
        end do
    end if
    nx=size(domain%mat,1); ny=size(domain%mat,2)
    if(.not.allocated(crop_state%slot)) then
        allocate(crop_state%lu(nx,ny),crop_state%slot(nx,ny),crop_state%occurrence(nx,ny),crop_state%cut(nx,ny), &
            crop_state%next_slot(nx,ny),crop_state%sow(nx,ny),crop_state%harvest_limit(nx,ny), &
            crop_state%last_harvest(nx,ny),crop_state%next_sow(nx,ny),crop_state%deadline(nx,ny),crop_state%offset(nx,ny),crop_state%pending_offset(nx,ny), &
            crop_state%days(nx,ny,4),crop_state%hs_days(nx,ny),crop_state%gdd(nx,ny),crop_state%vern(nx,ny), &
            crop_state%corr(nx,ny,2),crop_state%actual(nx,ny,4),crop_state%potential(nx,ny,4), &
            crop_state%biomass(nx,ny),crop_state%hs(nx,ny),crop_state%terminate(nx,ny),crop_state%carried(nx,ny))
        crop_state%lu=0; crop_state%slot=0; crop_state%occurrence=0; crop_state%cut=0; crop_state%next_slot=0
        crop_state%sow=0; crop_state%harvest_limit=0; crop_state%last_harvest=-huge(1)/2
        crop_state%next_sow=0; crop_state%deadline=0; crop_state%offset=0; crop_state%pending_offset=0
        crop_state%gdd=0; crop_state%vern=0; crop_state%corr=0; crop_state%actual=0; crop_state%potential=0
        crop_state%biomass=0; crop_state%hs=0; crop_state%hs_days=0; crop_state%days=0
        crop_state%terminate=.false.; crop_state%carried=.false.
        run_seed=sim%rand_seed
        if(run_seed==-999) run_seed=1
        if(.not.sim%repeatable) call system_clock(count=run_seed)
    else
        ! Preserve biological age when the existing soil warm-up rewinds weather.
        ! This handoff deliberately does not redesign the simulation loop.
        delta=sim%start_simulation%doy-last_date-1
        where(crop_state%slot>0)
            crop_state%sow=crop_state%sow+delta
            crop_state%harvest_limit=crop_state%harvest_limit+delta
            crop_state%carried=.true.
        end where
        where(crop_state%next_slot>0)
            crop_state%next_sow=crop_state%next_sow+delta
            crop_state%deadline=crop_state%deadline+delta
        end where
        crop_state%last_harvest=crop_state%last_harvest+delta
    end if
    prefix=trim(sim%path)
    if(warmup) prefix=trim(prefix)//'warmup_'
    open(newunit=event_unit,file=trim(prefix)//'crop_events.csv',status='replace',iostat=ios)
    if(ios/=0) call crop_fail('cannot write crop event log')
    write(event_unit,'(a)') 'date;i;j;landuse;slot;occurrence;cut;event;gdd;sowing_offset'
    open(newunit=harvest_unit,file=trim(prefix)//'crop_harvests.csv',status='replace',iostat=ios)
    if(ios/=0) call crop_fail('cannot write crop harvest log')
    write(harvest_unit,'(a)') 'date;i;j;landuse;slot;occurrence;cut;reason;sowing_date;gdd;biomass_potential_t_ha;'// &
        'yield_potential_t_ha;water_factor;stage_water_factor;heat_factor;yield_actual_t_ha;contains_warmup'
end subroutine

! Offset shifts the sowing window, never the already established crop's GDD curve.
integer function sowing_offset(sim,i,j,lu,slot,year,map) result(offset)
    type(simulation), intent(in) :: sim
    integer, intent(in) :: i,j,lu,slot,year,map(:,:)
    integer(kind=8) :: h
    if(sim%f_irandom) then
        offset=map(i,j)
        if(abs(offset)>365) call crop_fail("irandom value must be in -365..365 on every active cell")
    else
        ! Stable integer hash, independent of traversal order/thread scheduling.
        h=modulo(int(run_seed,8)+104729_8*i+13007_8*j+8191_8*lu+131_8*slot+year,2147483646_8)+1
        h=modulo(48271_8*h,2147483647_8)
        offset=int(modulo(h,2_8*sim%sowing_range+1))
        if(sim%rand_symmetry) offset=offset-sim%sowing_range
    end if
end function

subroutine window_dates(c,earliest,offset,start,deadline)
    type(crop_definition), intent(in) :: c
    integer, intent(in) :: earliest,offset
    integer, intent(out) :: start,deadline
    integer :: d,m,y,yy,a,b
    call calc_date(earliest,d,m,y)
    start=huge(1); deadline=huge(1)
    ! Include previous year's window when it extends across January.
    do yy=y-1,y+2
        a=calc_doy(1,1,yy)+c%sow_min-1+offset
        b=a+c%sow_delay
        if(b<earliest) cycle
        if(b<deadline) then
            start=a; deadline=b
        end if
    end do
end subroutine

subroutine schedule_crop(sim,i,j,slot,earliest,map)
    type(simulation), intent(in) :: sim
    integer, intent(in) :: i,j,slot,earliest,map(:,:)
    integer :: d,m,y,off,lu,a,b
    lu=crop_state%lu(i,j)
    call calc_date(earliest,d,m,y)
    off=sowing_offset(sim,i,j,lu,slot,y,map)
    call window_dates(rotations(lu)%crops(slot),earliest,off,a,b)
    crop_state%next_slot(i,j)=slot; crop_state%next_sow(i,j)=a; crop_state%deadline(i,j)=b
    crop_state%pending_offset(i,j)=off
end subroutine

subroutine select_first_crop(sim,i,j,date,lu,map)
    type(simulation), intent(in) :: sim
    integer, intent(in) :: i,j,date,lu,map(:,:)
    integer :: slot,best,a,b,deadline,d,m,y,off
    crop_state%lu(i,j)=lu
    best=0; deadline=huge(1)
    call calc_date(date,d,m,y)
    do slot=1,size(rotations(lu)%crops)
        if(maxval(rotations(lu)%crops(slot)%gdd)<=0) cycle
        off=sowing_offset(sim,i,j,lu,slot,y,map)
        call window_dates(rotations(lu)%crops(slot),date,off,a,b)
        if(b<deadline) then
            deadline=b; best=slot
        end if
    end do
    crop_state%next_slot(i,j)=0
    if(best>0) call schedule_crop(sim,i,j,best,date,map)
end subroutine

subroutine advance_crops(sim,date,meteo,domain,landuse,map,indices,weights,ze,pheno)
    type(simulation), intent(in) :: sim
    integer, intent(in) :: date,map(:,:),indices(:,:,:)
    type(meteo_mat), intent(in) :: meteo
    type(grid_i), intent(in) :: domain,landuse
    real(dp), intent(in) :: weights(:,:,:),ze
    type(crop_pars_matrices), intent(inout) :: pheno
    integer :: i,j,lu,slot,next,n,d,m,y,limit,earliest,off
    real(dp) :: smooth,g
    logical :: forced
    do j=1,size(domain%mat,2)
        do i=1,size(domain%mat,1)
            if(domain%mat(i,j)==domain%header%nan) then
                if(crop_state%slot(i,j)>0) call harvest_crop(i,j,date,'domain_change')
                call clear_crop_cell(pheno,i,j)
                crop_state%slot(i,j)=0; crop_state%next_slot(i,j)=0; crop_state%lu(i,j)=0
                cycle
            end if
            lu=landuse%mat(i,j)
            if(lu<1.or.lu>size(rotations)) call crop_fail('land use outside crop database')
            if(crop_state%lu(i,j)/=lu) then
                if(crop_state%slot(i,j)>0) then
                    call harvest_crop(i,j,date,'landuse_change')
                    crop_state%last_harvest(i,j)=date-1
                end if
                crop_state%slot(i,j)=0
                call select_first_crop(sim,i,j,date,lu,map)
            end if
            smooth=smoothed_temperature(sim,i,j,date,indices,weights)
            slot=crop_state%slot(i,j)
            if(slot==0) then
                next=crop_state%next_slot(i,j)
                if(next>0) then
                    associate(c=>rotations(lu)%crops(next))
                    earliest=max(crop_state%next_sow(i,j),crop_state%last_harvest(i,j)+c%gap+1)
                    if(date>=earliest.and.(smooth>=c%tsow.or.date>=crop_state%deadline(i,j))) then
                        slot=next; crop_state%slot(i,j)=slot
                        crop_state%occurrence(i,j)=crop_state%occurrence(i,j)+1
                        crop_state%cut(i,j)=1; crop_state%gdd(i,j)=0; crop_state%vern(i,j)=0
                        crop_state%sow(i,j)=date; crop_state%carried(i,j)=.false.
                        crop_state%offset(i,j)=crop_state%pending_offset(i,j)
                        forced=smooth<c%tsow.or.date>crop_state%deadline(i,j)
                        call calc_date(date,d,m,y)
                        limit=calc_doy(1,1,y)+c%harvest_max-1
                        if(limit<date) limit=calc_doy(1,1,y+1)+c%harvest_max-1
                        crop_state%harvest_limit(i,j)=limit
                        n=size(rotations(lu)%crops); next=modulo(slot,n)+1
                        ! Zero-GDD entries represent bare soil, not living crop occurrences.
                        do d=1,n
                            if(maxval(rotations(lu)%crops(next)%gdd)>0) exit
                            next=modulo(next,n)+1
                        end do
                        earliest=date+1
                        if(next==slot) earliest=crop_state%deadline(i,j)+1
                        call schedule_crop(sim,i,j,next,earliest,map)
                        call project_correction(sim,c,i,j,date,meteo%lat(i,j),indices,weights)
                        if(forced) then
                            call event(i,j,date,'forced_sowing')
                        else
                            call event(i,j,date,'sowing')
                        end if
                    end if
                    end associate
                end if
            end if
            call clear_crop_cell(pheno,i,j)
            if(slot==0) cycle
            associate(c=>rotations(lu)%crops(slot))
            if(crop_state%cut(i,j)>1.and.crop_state%gdd(i,j)==0) &
                call project_correction(sim,c,i,j,date,meteo%lat(i,j),indices,weights)
            ! Preserve cropcoef's zero-development sowing day.
            if(date>crop_state%sow(i,j)) call advance_thermal(c,meteo%T_max(i,j),meteo%T_min(i,j), &
                smooth,daylight(date,meteo%lat(i,j)),crop_state%vern(i,j),crop_state%gdd(i,j))
            g=min(crop_state%gdd(i,j),maxval(c%gdd))
            call fill_crop_cell(c,pheno,i,j,g,ze,sim%crop_co2)
            crop_state%terminate(i,j)=date>=crop_state%harvest_limit(i,j)
            next=crop_state%next_slot(i,j)
            if(next>0) crop_state%terminate(i,j)=crop_state%terminate(i,j).or. &
                date>=crop_state%deadline(i,j)-rotations(lu)%crops(next)%gap-1
            if(crop_state%gdd(i,j)>=maxval(c%gdd).and.crop_state%cut(i,j)>=c%cuts) crop_state%terminate(i,j)=.true.
            end associate
        end do
    end do
    last_date=date
end subroutine

subroutine clear_crop_cell(p,i,j)
    type(crop_pars_matrices), intent(inout) :: p
    integer, intent(in) :: i,j
    p%k_cb(i,j)=0; p%lai(i,j)=0; p%h(i,j)=0; p%d_r(i,j)=0; p%cn_day(i,j)=0; p%f_c(i,j)=0
    p%irrigation_class(i,j)=0; p%cn_class(i,j)=1; p%p(i,j)=0; p%a(i,j)=0; p%d_t_max(i,j)=0
    p%RF_t_max(i,j)=0; p%T_lim(i,j)=0; p%T_crit(i,j)=0; p%HI(i,j)=0; p%Ky_tot(i,j)=0
    p%Ky_pheno(i,j,:)=0; p%k_cb_low(i,j)=0; p%k_cb_mid(i,j)=0; p%k_cb_high(i,j)=0
    p%wp_adj(i,j)=0; p%n_crop_in_year(i,j)=0; p%pheno_idx(i,j)=0; p%r_stress(i,j)=0
end subroutine

subroutine fill_crop_cell(c,p,i,j,g,ze,co2)
    type(crop_definition), intent(in) :: c
    type(crop_pars_matrices), intent(inout) :: p
    integer, intent(in) :: i,j
    real(dp), intent(in) :: g,ze,co2
    real(dp) :: factor,correction
    p%k_cb(i,j)=interpolate_crop(c,g,1)
    correction=crop_state%corr(i,j,1)
    if(g<c%stage_gdd(1)) then
        correction=0
    else if(g<c%stage_gdd(2)) then
        factor=(g-c%stage_gdd(1))/max(tiny(1.0_dp),c%stage_gdd(2)-c%stage_gdd(1))
        correction=correction*factor
    else if(g>c%stage_gdd(3)) then
        factor=(g-c%stage_gdd(3))/max(tiny(1.0_dp),maxval(c%gdd)-c%stage_gdd(3))
        correction=(1-factor)*correction+factor*crop_state%corr(i,j,2)
    end if
    p%k_cb(i,j)=max(0.0_dp,p%k_cb(i,j)+correction)
    p%lai(i,j)=interpolate_crop(c,g,2); p%h(i,j)=interpolate_crop(c,g,3)
    p%d_r(i,j)=interpolate_crop(c,g,4); p%cn_day(i,j)=nint(interpolate_crop(c,g,5))
    p%f_c(i,j)=interpolate_crop(c,g,6); p%r_stress(i,j)=interpolate_crop(c,g,7)
    p%irrigation_class(i,j)=c%irrigated; p%cn_class(i,j)=c%cn_class; p%p(i,j)=c%praw
    p%a(i,j)=c%interception; p%d_t_max(i,j)=max(0.0_dp,maxval(c%values(:,4))-ze); p%RF_t_max(i,j)=c%rft
    p%T_lim(i,j)=c%tlim; p%T_crit(i,j)=c%tcrit; p%HI(i,j)=c%hi; p%Ky_tot(i,j)=c%kyt
    p%Ky_pheno(i,j,:)=c%ky; p%wp_adj(i,j)=adjusted_wp(c,co2)
    p%n_crop_in_year(i,j)=crop_state%slot(i,j)
    p%pheno_idx(i,j)=1
    if(g>=c%stage_gdd(1)) p%pheno_idx(i,j)=2
    if(g>=c%stage_gdd(2)) p%pheno_idx(i,j)=3
    if(g>=c%stage_gdd(3)) p%pheno_idx(i,j)=4
    p%k_cb_low(i,j)=minval(c%values(:,1)); p%k_cb_mid(i,j)=interpolate_crop(c,c%stage_gdd(1),1)
    p%k_cb_high(i,j)=maxval(c%values(:,1))+crop_state%corr(i,j,1)
end subroutine

subroutine finish_crop_day(date,domain,p,meteo,actual,potential)
    integer, intent(in) :: date
    type(grid_i), intent(in) :: domain
    type(crop_pars_matrices), intent(in) :: p
    type(meteo_mat), intent(in) :: meteo
    real(dp), intent(in) :: actual(:,:),potential(:,:)
    integer :: i,j,lu,slot,stage
    real(dp) :: fraction
    do j=1,size(domain%mat,2)
        do i=1,size(domain%mat,1)
            if(domain%mat(i,j)==domain%header%nan) cycle
            slot=crop_state%slot(i,j)
            if(slot==0) cycle
            lu=crop_state%lu(i,j); stage=p%pheno_idx(i,j)
            associate(c=>rotations(lu)%crops(slot))
            crop_state%actual(i,j,stage)=crop_state%actual(i,j,stage)+actual(i,j)
            crop_state%potential(i,j,stage)=crop_state%potential(i,j,stage)+potential(i,j)
            crop_state%days(i,j,stage)=crop_state%days(i,j,stage)+1
            if(meteo%et0(i,j)>0) crop_state%biomass(i,j)=crop_state%biomass(i,j)+p%wp_adj(i,j)*potential(i,j)/meteo%et0(i,j)
            fraction=crop_state%gdd(i,j)/maxval(c%gdd)
            ! Explicit first-pass choice: thermal progress replaces unknown final calendar duration.
            if(fraction>=0.45_dp.and.fraction<0.75_dp) then
                crop_state%hs_days(i,j)=crop_state%hs_days(i,j)+1
                crop_state%hs(i,j)=crop_state%hs(i,j)+max(0.0_dp,min(1.0_dp, &
                    (c%tlim-meteo%T_ave(i,j))/(c%tlim-c%tcrit)))
            end if
            if(crop_state%terminate(i,j)) then
                if(crop_state%gdd(i,j)>=maxval(c%gdd)) then
                    call harvest_crop(i,j,date,'maturity')
                else if(date>=crop_state%harvest_limit(i,j)) then
                    call harvest_crop(i,j,date,'harvest_deadline')
                else
                    call harvest_crop(i,j,date,'rotation_deadline')
                end if
                crop_state%slot(i,j)=0; crop_state%last_harvest(i,j)=date
            else if(crop_state%gdd(i,j)>=maxval(c%gdd)) then
                call harvest_crop(i,j,date,'cut')
                crop_state%cut(i,j)=crop_state%cut(i,j)+1; crop_state%gdd(i,j)=0
            end if
            end associate
        end do
    end do
end subroutine

subroutine harvest_crop(i,j,date,reason)
    integer, intent(in) :: i,j,date
    character(len=*), intent(in) :: reason
    real(dp) :: fw,fs,fh,yp,ya,total,a,p
    integer :: stage,lu,slot
    lu=crop_state%lu(i,j); slot=crop_state%slot(i,j)
    associate(c=>rotations(lu)%crops(slot))
    a=sum(crop_state%actual(i,j,:)); p=sum(crop_state%potential(i,j,:))
    fw=1
    if(p>0) fw=max(0.0_dp,min(1.0_dp,1-c%kyt*(1-a/p)))
    fs=1; total=sum(crop_state%days(i,j,:))
    do stage=1,4
        p=crop_state%potential(i,j,stage)
        if(p<=0.or.total<=0) cycle
        fs=fs*max(0.0_dp,min(1.0_dp,1-c%ky(stage)*(1-crop_state%actual(i,j,stage)/p)))** &
            (crop_state%days(i,j,stage)/total)
    end do
    fh=1
    if(crop_state%hs_days(i,j)>0) fh=crop_state%hs(i,j)/crop_state%hs_days(i,j)
    yp=crop_state%biomass(i,j)*c%hi; ya=yp*min(fw,fs)*fh
    write(harvest_unit,'(*(g0,:,";"))') date_text(date),i,j,lu,slot,crop_state%occurrence(i,j),crop_state%cut(i,j), &
        trim(reason),date_text(crop_state%sow(i,j)),crop_state%gdd(i,j),crop_state%biomass(i,j),yp,fw,fs,fh,ya,crop_state%carried(i,j)
    call event(i,j,date,reason)
    crop_state%actual(i,j,:)=0; crop_state%potential(i,j,:)=0; crop_state%days(i,j,:)=0
    crop_state%biomass(i,j)=0; crop_state%hs(i,j)=0; crop_state%hs_days(i,j)=0
    end associate
end subroutine

subroutine event(i,j,date,name)
    integer, intent(in) :: i,j,date
    character(len=*), intent(in) :: name
    write(event_unit,'(*(g0,:,";"))') date_text(date),i,j,crop_state%lu(i,j),crop_state%slot(i,j), &
        crop_state%occurrence(i,j),crop_state%cut(i,j),trim(name),crop_state%gdd(i,j),crop_state%offset(i,j)
end subroutine
function date_text(date) result(text)
    integer, intent(in) :: date
    character(len=10) :: text
    integer :: d,m,y
    call calc_date(date,d,m,y)
    write(text,'(i4.4,"-",i2.2,"-",i2.2)') y,m,d
end function
subroutine end_crop_run()
    integer :: i,j
    ! Living crops are explicitly reported as unfinished; never harvested merely for output.
    do j=1,size(crop_state%slot,2)
        do i=1,size(crop_state%slot,1)
            if(crop_state%slot(i,j)>0) call event(i,j,last_date,'ongoing_at_run_end')
        end do
    end do
    close(event_unit); close(harvest_unit)
end subroutine

real(dp) function local_weather(sim,i,j,date,column,indices,weights) result(v)
    type(simulation), intent(in) :: sim
    integer, intent(in) :: i,j,date,column,indices(:,:,:)
    real(dp), intent(in) :: weights(:,:,:)
    integer :: index_day,k,d,m,y,lookup,first_year,last_year
    logical :: interpolate
    lookup=date
    if(date<weather_start.or.date>weather_end) then
        ! Boundary lookahead uses the corresponding date of the nearest available year.
        call calc_date(date,d,m,y)
        call calc_date(weather_start,index_day,k,first_year)
        call calc_date(weather_end,index_day,k,last_year)
        y=max(first_year,min(last_year,y)); lookup=calc_doy(d,m,y)
    end if
    index_day=max(1,min(size(weather,1),lookup-weather_start+1))
    interpolate=sim%interpolate_temp
    if(column==3) interpolate=sim%interpolate_hum
    if(column==4) interpolate=sim%interpolate_wind
    v=weather(index_day,indices(i,j,1),column)
    if(.not.interpolate) return
    v=0
    do k=1,size(indices,3)
        v=v+weights(i,j,k)*weather(index_day,indices(i,j,k),column)
    end do
end function
real(dp) function smoothed_temperature(sim,i,j,date,indices,weights) result(t)
    type(simulation), intent(in) :: sim
    integer, intent(in) :: i,j,date,indices(:,:,:)
    real(dp), intent(in) :: weights(:,:,:)
    integer :: k,n
    n=sim%crop_temperature_window; t=0
    do k=-n,n
        t=t+0.5_dp*(local_weather(sim,i,j,date+k,1,indices,weights)+local_weather(sim,i,j,date+k,2,indices,weights))
    end do
    t=t/(2*n+1)
end function
real(dp) function daylight(date,latitude) result(hours)
    integer, intent(in) :: date
    real(dp), intent(in) :: latitude
    integer :: d,m,y,doy
    real(dp) :: delta,argument
    call calc_date(date,d,m,y); doy=date-calc_doy(1,1,y)+1
    delta=0.409_dp*sin(2*pi*doy/365-1.39_dp)
    argument=-tan(latitude*pi/180)*tan(delta)
    hours=24/pi*acos(max(-1.0_dp,min(1.0_dp,argument)))
end function
pure real(dp) function thermal_units(tmax,tmin,base,cutoff) result(gdd)
    real(dp), intent(in) :: tmax,tmin,base,cutoff
    gdd=sine_above(tmax,tmin,base)-sine_above(tmax,tmin,cutoff)
end function
pure real(dp) function sine_above(tmax,tmin,threshold) result(g)
    real(dp), intent(in) :: tmax,tmin,threshold
    real(dp) :: average,amplitude,theta
    average=0.5_dp*(tmax+tmin); amplitude=0.5_dp*(tmax-tmin)
    if(tmin>=threshold) then
        g=average-threshold
    else if(tmax<=threshold) then
        g=0
    else
        theta=asin(max(-1.0_dp,min(1.0_dp,(threshold-average)/amplitude)))
        g=((average-threshold)*(pi/2-theta)+amplitude*cos(theta))/pi
    end if
end function
pure subroutine advance_thermal(c,tmax,tmin,smooth,dl,vern,gdd)
    type(crop_definition), intent(in) :: c
    real(dp), intent(in) :: tmax,tmin,smooth,dl
    real(dp), intent(inout) :: vern,gdd
    real(dp) :: vf,pf,effect
    vf=1; pf=1
    if(c%vern) then
        effect=max(0.0_dp,min(1.0_dp,min((smooth-c%tvmin+c%vslope)/c%vslope, &
            (c%tvmax+c%vslope-smooth)/c%vslope)))
        vern=vern+effect
        ! Retain CropCoef's existing response, including VF=1 outside [Vstart,Vend].
        if(vern>=c%vstart.and.vern<=c%vend) vf=c%vfmin+(1-c%vfmin)*(vern-c%vstart)/(c%vend-c%vstart)
    end if
    if(c%photo==1) pf=max(0.0_dp,min(1.0_dp,(dl-c%dl_if)/(c%dl_ins-c%dl_if)))
    if(c%photo==2) pf=max(0.0_dp,min(1.0_dp,(c%dl_if-dl)/(c%dl_if-c%dl_ins)))
    gdd=gdd+thermal_units(tmax,tmin,c%tbase,c%tcut)*min(vf,pf)
end subroutine

subroutine project_correction(sim,c,i,j,date,latitude,indices,weights)
    type(simulation), intent(in) :: sim
    type(crop_definition), intent(in) :: c
    integer, intent(in) :: i,j,date,indices(:,:,:)
    real(dp), intent(in) :: latitude,weights(:,:,:)
    real(dp) :: g,v,totals(3,2),h,rh,wind,tmax,tmin
    integer :: day,phase,counts(2),limit,next
    crop_state%corr(i,j,:)=0
    if(.not.c%adjust_kcb) return
    totals=0; counts=0; g=crop_state%gdd(i,j); v=crop_state%vern(i,j)
    next=crop_state%next_slot(i,j)
    limit=crop_state%harvest_limit(i,j)
    if(next>0) limit=min(limit,crop_state%deadline(i,j)-rotations(crop_state%lu(i,j))%crops(next)%gap-1)
    do day=date+1,limit
        tmax=local_weather(sim,i,j,day,1,indices,weights); tmin=local_weather(sim,i,j,day,2,indices,weights)
        call advance_thermal(c,tmax,tmin,smoothed_temperature(sim,i,j,day,indices,weights),daylight(day,latitude),v,g)
        phase=0
        if(g>=c%stage_gdd(2).and.g<c%stage_gdd(3)) phase=1
        if(g>=c%stage_gdd(3)) phase=2
        if(phase>0) then
            totals(:,phase)=totals(:,phase)+[local_weather(sim,i,j,day,3,indices,weights), &
                local_weather(sim,i,j,day,4,indices,weights),interpolate_crop(c,g,3)]
            counts(phase)=counts(phase)+1
        end if
        if(g>=maxval(c%gdd)) exit
    end do
    do phase=1,2
        if(counts(phase)==0) cycle ! truncated crops may never reach these stages
        if(phase==1.and.maxval(c%values(:,1))<=0.45_dp) cycle
        if(phase==2.and.c%values(size(c%gdd),1)<=0.45_dp) cycle
        rh=max(20.0_dp,min(80.0_dp,totals(1,phase)/counts(phase)))
        wind=max(1.0_dp,min(6.0_dp,totals(2,phase)/counts(phase)))
        h=max(0.1_dp,min(10.0_dp,totals(3,phase)/counts(phase)))
        crop_state%corr(i,j,phase)=(0.04_dp*(wind-2)-0.004_dp*(rh-45))*(h/3)**0.3_dp
    end do
end subroutine
end module
