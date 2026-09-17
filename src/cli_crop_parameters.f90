! Raw CropCoef-compatible inputs for standalone IdrAgra.
! Physiological definitions/formulas originate in cropcoeff (E. A. Chiaradia), GPL-2.0-or-later.
module cli_crop_parameters
use mod_constants, only: dp
use mod_utility, only: lower_case
use mod_parameters, only: simulation
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
implicit none
private
public :: crop_definition, rotation_definition, rotations, init_crop_database, interpolate_crop, adjusted_wp
public :: crop_fail, words, clean_line
integer, parameter :: max_points=512, max_words=64
real(dp), parameter :: missing=-9999.0_dp
type crop_definition
    character(len=255) :: name=''
    integer :: sow_min=1, sow_delay=0, harvest_max=365, cuts=1, gap=0, photo=0, cn_class=1, irrigated=0
    logical :: vern=.false., adjust_kcb=.true.
    real(dp) :: tsow=0, tbase=0, tcut=30, tvmin=0, tvmax=10, vslope=7, vstart=10, vend=50, vfmin=0
    real(dp) :: dl_if=8, dl_ins=20, wp=0, fsink=0, tcrit=35, tlim=45, hi=0, kyt=1, ky(4)=1
    real(dp) :: praw=0.5_dp, interception=0.5_dp, rft=1
    real(dp), allocatable :: gdd(:), values(:,:) ! Kcb, LAI, height, roots, CN, fc, resistance, Ky
    real(dp) :: stage_gdd(3)=0 ! development start, mid-season start, late-season start
end type
type rotation_definition
    type(crop_definition), allocatable :: crops(:)
end type
type(rotation_definition), allocatable, save :: rotations(:)
contains
subroutine crop_fail(message)
    character(len=*), intent(in) :: message
    print *, 'Crop input error: ',trim(message)
    error stop 1
end subroutine
function clean_line(raw) result(line)
    character(len=*), intent(in) :: raw
    character(len=len(raw)) :: line
    integer :: k
    line=raw
    k=index(line,'#')
    if(k>0) line=line(:k-1)
    do k=1,len_trim(line)
        if(line(k:k)==achar(9)) line(k:k)=' '
    end do
    line=adjustl(line)
end function
subroutine words(line, token, n)
    character(len=*), intent(in) :: line
    character(len=*), intent(out) :: token(:)
    integer, intent(out) :: n
    integer :: i,j
    n=0; i=1; token=''
    do while(i<=len_trim(line))
        if(line(i:i)==' ') then
            i=i+1; cycle
        end if
        j=i
        do while(j<=len_trim(line))
            if(line(j:j)==' ') exit
            j=j+1
        end do
        n=n+1
        if(n>size(token)) call crop_fail('too many fields: '//trim(line))
        if(j-i>len(token)) call crop_fail('field too long: '//trim(line))
        token(n)=line(i:j-1); i=j
    end do
end subroutine
subroutine init_crop_database(sim)
    type(simulation), intent(inout) :: sim
    character(len=4096) :: line, raw
    character(len=255) :: tok(max_words), path
    integer :: u,ios,n,id,k,nc
    if(.not.allocated(rotations)) then
        allocate(rotations(sim%n_lus))
        open(newunit=u,file=trim(sim%rotation_file),status='old',action='read',iostat=ios)
        if(ios/=0) call crop_fail('cannot open '//trim(sim%rotation_file))
        do
            read(u,'(a)',iostat=ios) raw
            if(ios<0) exit
            if(ios/=0) call crop_fail('reading '//trim(sim%rotation_file))
            line=clean_line(raw)
            if(len_trim(line)==0) cycle
            call words(line,tok,n)
            call lower_case(tok(1))
            if(tok(1)=='endtable') exit
            if(tok(1)=='cr_id') cycle
            read(tok(1),*,iostat=ios) id
            if(ios/=0) call crop_fail('expected land-use ID: '//trim(line))
            if(id<1.or.id>sim%n_lus) call crop_fail('land-use ID outside SoilUsesNum')
            if(allocated(rotations(id)%crops)) call crop_fail('duplicate land-use ID')
            nc=count(tok(2:n)/='*')
            allocate(rotations(id)%crops(nc)); nc=0
            do k=2,n
                if(tok(k)=='*') cycle
                nc=nc+1
                path=trim(sim%crop_path)//'/'//trim(tok(k))
                call read_crop(path,rotations(id)%crops(nc))
            end do
        end do
        close(u)
        do id=1,sim%n_lus
            if(.not.allocated(rotations(id)%crops)) call crop_fail('rotation database must define every land-use ID')
        end do
    end if
    sim%n_crops=1
    do id=1,size(rotations)
        sim%n_crops=max(sim%n_crops,size(rotations(id)%crops))
    end do
    if(associated(sim%res_canopy)) deallocate(sim%res_canopy)
    allocate(sim%res_canopy(sim%meteo_years))
    sim%res_canopy=70.0_dp
    if(sim%crop_co2>0) then
        if(sim%crop_co2>=1155) call crop_fail('CropCO2 must be below 1155 ppm for canopy-resistance formula')
        sim%res_canopy=100.0_dp/(1.4_dp-0.4_dp*sim%crop_co2/330.0_dp)/(0.5_dp*24*0.12_dp)
    end if
end subroutine
subroutine read_crop(path,c)
    character(len=*), intent(in) :: path
    type(crop_definition), intent(out) :: c
    character(len=4096) :: raw,line,key
    character(len=255) :: tok(max_words), headers(max_words)
    real(dp) :: table(max_points,9),v
    integer :: u,ios,n,ncol,row,k,p,q,lo,hi,flags(5)
    logical :: in_table
    c%name=path; table=missing; row=0; in_table=.false.; ncol=0; flags=0
    open(newunit=u,file=trim(path),status='old',action='read',iostat=ios)
    if(ios/=0) call crop_fail('cannot open '//trim(path))
    do
        read(u,'(a)',iostat=ios) raw
        if(ios<0) exit
        if(ios/=0) call crop_fail('reading '//trim(path))
        line=clean_line(raw)
        if(len_trim(line)==0) cycle
        key=line; call lower_case(key)
        if(key=='endtable') exit
        if(in_table) then
            call words(line,tok,n)
            if(n/=ncol) call crop_fail('wrong table column count in '//trim(path))
            row=row+1
            if(row>max_points) call crop_fail('too many crop curve points')
            do k=1,n
                if(tok(k)=='*') cycle
                read(tok(k),*,iostat=ios) table(row,k)
                if(ios/=0) call crop_fail('invalid table value in '//trim(path))
                if(.not.ieee_is_finite(table(row,k))) call crop_fail('non-finite crop value')
            end do
            cycle
        end if
        p=index(line,'=')
        if(p==0) then
            call words(key,headers,ncol)
            if(ncol<5.or.ncol>9) call crop_fail('expected GDD table in '//trim(path))
            if(any(headers(1:5)/=[character(len=255)::'gdd','kcb','lai','hc','sr'])) &
                call crop_fail('expected GDD Kcb LAI Hc Sr columns')
            if(ncol>=6) then
                if(headers(6)/='cn') call crop_fail('expected CN column')
            end if
            if(ncol>=7) then
                if(headers(7)/='fc') call crop_fail('expected fc column')
            end if
            if(ncol>=8) then
                if(headers(8)/='r_stress') call crop_fail('expected r_stress column')
            end if
            if(ncol>=9) then
                if(headers(9)/='ky') call crop_fail('expected Ky column')
            end if
            in_table=.true.; cycle
        end if
        key=trim(line(:p-1)); call lower_case(key)
        read(line(p+1:),*,iostat=ios) v
        if(ios/=0) call crop_fail('expected numeric parameter '//trim(key)//' in '//trim(path))
        if(.not.ieee_is_finite(v)) call crop_fail('non-finite parameter '//trim(key))
        select case(trim(key))
        case('sowingdate_min'); c%sow_min=nint(v)
        case('sowingdelay_max'); c%sow_delay=nint(v)
        case('harvestdate_max'); c%harvest_max=nint(v)
        case('harvnum_max'); c%cuts=nint(v)
        case('cropsoverlap'); c%gap=nint(v)
        case('ph_r'); c%photo=nint(v)
        case('cl_cn'); c%cn_class=nint(v)
        case('irrigation'); c%irrigated=nint(v)
        case('tsowing'); c%tsow=v
        case('tdaybase'); c%tbase=v
        case('tcutoff'); c%tcut=v
        case('tv_min'); c%tvmin=v
        case('tv_max'); c%tvmax=v
        case('vslope'); c%vslope=v
        case('vstart'); c%vstart=v
        case('vend'); c%vend=v
        case('vfmin'); c%vfmin=v
        case('daylength_if'); c%dl_if=v
        case('daylength_ins'); c%dl_ins=v
        case('wp'); c%wp=v
        case('fsink'); c%fsink=v
        case('tcrit_hs'); c%tcrit=v
        case('tlim_hs'); c%tlim=v
        case('hi'); c%hi=v
        case('kyt'); c%kyt=v
        case('ky1'); c%ky(1)=v
        case('ky2'); c%ky(2)=v
        case('ky3'); c%ky(3)=v
        case('ky4'); c%ky(4)=v
        case('praw'); c%praw=v
        case('ainterception'); c%interception=v
        case('rft'); c%rft=v
        case('vern'); c%vern=v/=0
        case('adj_flag'); c%adjust_kcb=v/=0
        case('ke','kt') ! obsolete root-fraction parameters
        case default
            call crop_fail('unknown parameter '//trim(key)//' in '//trim(path))
        end select
        select case(trim(key))
        case('sowingdate_min'); flags(1)=1
        case('harvestdate_max'); flags(2)=1
        case('tdaybase'); flags(3)=1
        case('tcutoff'); flags(4)=1
        case('sowingdelay_max'); flags(5)=1
        end select
    end do
    close(u)
    if(row<1.or.any(flags==0)) call crop_fail('missing required parameters/table in '//trim(path))
    if(any(table(1:row,1)<0)) call crop_fail('missing/negative GDD')
    if(any(table(2:row,1)<=table(1:row-1,1))) call crop_fail('GDD must be strictly increasing')
    if(c%tcut<=c%tbase.or.c%sow_min<1.or.c%sow_min>366.or.c%harvest_max<1.or.c%harvest_max>366 &
        .or.c%sow_delay<0.or.c%sow_delay>365.or.c%gap<0.or.c%gap>365.or.c%cuts<1) call crop_fail('invalid crop limits in '//trim(path))
    if(c%vern.and.(c%vend<=c%vstart.or.c%vslope<=0.or.c%tvmax<c%tvmin)) call crop_fail('invalid vernalization limits')
    if(c%photo<0.or.c%photo>2) call crop_fail('invalid photoperiod class')
    if(c%photo==1.and.c%dl_ins<=c%dl_if) call crop_fail('invalid long-day thresholds')
    if(c%photo==2.and.c%dl_if<=c%dl_ins) call crop_fail('invalid short-day thresholds')
    if(c%cn_class<1.or.c%cn_class>7.or.c%rft<0.or.c%rft>1.or.c%tlim<=c%tcrit) call crop_fail('invalid crop water/yield parameters')
    allocate(c%gdd(row),c%values(row,8))
    c%gdd=table(1:row,1); c%values=table(1:row,2:9)
    ! Interpolate missing knots in thermal time, preserving negative fc as the computed-cover sentinel.
    do k=1,8
        if(k==5) cycle ! discrete CN, handled separately
        if(all(c%values(:,k)==missing)) then
            select case(k)
            case(1:4); call crop_fail('missing growth curve in '//trim(path))
            case(6); c%values(:,k)=-1
            case(7); c%values(:,k)=0
            case(8); c%values(:,k)=c%kyt
            end select
        end if
        do p=1,row
            if(c%values(p,k)/=missing) cycle
            lo=p-1; hi=p+1
            do while(lo>=1)
                if(c%values(lo,k)/=missing) exit
                lo=lo-1
            end do
            do while(hi<=row)
                if(c%values(hi,k)/=missing) exit
                hi=hi+1
            end do
            if(lo<1) then
                c%values(p,k)=c%values(hi,k)
            else if(hi>row) then
                c%values(p,k)=c%values(lo,k)
            else
                c%values(p,k)=c%values(lo,k)+(c%values(hi,k)-c%values(lo,k))* &
                    (c%gdd(p)-c%gdd(lo))/(c%gdd(hi)-c%gdd(lo))
            end if
        end do
    end do
    do p=1,row
        if(c%values(p,5)==missing) then
            c%values(p,5)=1
            if(c%values(p,1)==0) c%values(p,5)=0
            if(c%values(p,1)>=0.45_dp) c%values(p,5)=2
        end if
    end do
    if(any(c%values(:,1:4)<0)) call crop_fail('negative growth properties')
    if(any(c%values(:,5)<0).or.any(c%values(:,5)>2)) call crop_fail('CN stage must be 0, 1 or 2')
    ! Infer stage boundaries once from the uncorrected Kcb curve (never from today's weather-adjusted Kcb).
    p=maxloc(c%values(:,1),dim=1); q=p
    do k=p,row
        if(c%values(k,1)==c%values(p,1)) q=k
    end do
    lo=1
    do k=2,p
        if(c%values(k,1)>0.and.c%values(k,1)==c%values(k-1,1)) lo=k
    end do
    c%stage_gdd=[c%gdd(lo),c%gdd(p),c%gdd(q)]
end subroutine
pure real(dp) function interpolate_crop(c,gdd,column) result(v)
    type(crop_definition), intent(in) :: c
    real(dp), intent(in) :: gdd
    integer, intent(in) :: column
    integer :: k,n
    n=size(c%gdd)
    v=c%values(1,column)
    if(gdd<=c%gdd(1)) return
    do k=2,n
        if(gdd<c%gdd(k)) then
            v=c%values(k-1,column)
            if(column/=5) v=v+(c%values(k,column)-v)*(gdd-c%gdd(k-1))/(c%gdd(k)-c%gdd(k-1))
            return
        end if
    end do
    v=c%values(n,column)
end function
pure real(dp) function adjusted_wp(c,co2) result(wp)
    type(crop_definition), intent(in) :: c
    real(dp), intent(in) :: co2
    real(dp) :: ftype,w,fco2
    wp=c%wp
    if(co2<=0) return
    ftype=max(0.0_dp,min(1.0_dp,(40-c%wp*100)/20))
    w=max(0.0_dp,min(1.0_dp,(co2-369.41_dp)/(550-369.41_dp)))
    fco2=(co2/369.41_dp)/(1+(co2-369.41_dp)*((1-w)*0.000138_dp+ &
        w*(c%fsink*0.000138_dp+(1-c%fsink)*0.001165_dp)))
    wp=(1+ftype*(fco2-1))*c%wp
end function
end module
