module cli_crop_parameters!
    use mod_constants, only: sp, dp
    use mod_utility, only: seek_un, lower_case, string_to_integers, string_to_reals, split_string, count_element
    use mod_parameters, only: simulation
    use mod_meteo, only: meteo_info
    use mod_crop_phenology
    use mod_system
    
    implicit none!

    interface read_crop_pars!
        module procedure read_crop_pars_r, read_crop_pars_i
    end interface!
    
    interface spread_col
        module procedure spread_col_i, spread_col_r
    end interface    
    
    contains!

    subroutine open_daily_crop_par_file(file_unit, file_name, error_flag)!
        ! open day-dependent crop parameters file
        character(len=*), intent(in) :: file_name
        integer, intent(out) :: file_unit
        integer, intent(out) :: error_flag
        integer :: ios ! check opening errors
        error_flag = 0!
        call seek_un( error_flag, file_unit) ! find free file unit
        open(unit=file_unit, file=trim(file_name), status='old', action="read", iostat=ios )!
        if (ios /= 0 ) then
            print *, "Cannot open file ", trim(file_name), ". The specified file does not exist. &
                & Execution will be aborted..."
            stop
        end if
        read( file_unit,*)!
    end subroutine open_daily_crop_par_file!
    
    subroutine init_crop_par_from_file(file_name, n_crop, n_crop_alt, string_elements, n_crops_by_year, error_flag)!
        ! init static crop parameters from parameter file
        implicit none
        character(len=*), intent(in) :: file_name!
        integer, intent(in) :: n_crop
        integer, intent(inout) :: n_crop_alt
        integer, intent(out) :: string_elements
        integer, dimension(n_crop), intent(out) :: n_crops_by_year
        integer, intent(out) :: error_flag!
        integer :: free_unit
        integer :: ios
        integer :: line, p
        character(len=n_crop*20) :: buffer, label !EAC: use mcrop_max x 20
        character(len=10), dimension(:), allocatable :: dummy, dummy_clean
        !!
        error_flag = 0!
        call seek_un(error_flag, free_unit) ! find for free unit
        open(unit=free_unit, file=trim(file_name), status='old', action="read", iostat=ios )!
        if (ios /= 0 ) then
            print *, "Cannot open file ", trim(file_name), ". The specified file does not exist. &
                & Execution will be aborted..."
            stop
        end if
        
        do while (ios == 0)
            read (free_unit, '(A)', iostat=ios) buffer
            if (ios == 0) then
                line = line + 1
                buffer = trim(buffer)
                call lower_case(buffer)
                p = scan(buffer, achar(9))  ! find the first tab ---> tab=achar(9)!
                label = buffer(1:p-1)
                buffer = buffer(p+1:)

                select case (label)
                    case ('var')
                        ! read the header line and divide elements
                        allocate(dummy(n_crop*2))    ! TODO: to add 3 crops, edit 2 to 3
                        call split_string(buffer, achar(9), dummy, string_elements)
                        allocate(dummy_clean(string_elements))
                        dummy_clean = dummy(1:string_elements)
                        deallocate(dummy)
                        n_crops_by_year = 0
                        ! find and count duplicates
                        call count_element(dummy_clean,n_crops_by_year)
                        deallocate(dummy_clean)
                        n_crop_alt = maxval(n_crops_by_year)
                   case default
                end select
            end if
        end do
        close (free_unit)
    end subroutine init_crop_par_from_file!
    
    subroutine read_water_prod_file(file_name, string_elements, n_crops_by_year, unit_param, error_flag)!
        ! read water productivity related parameters
        implicit none
        character(len=*), intent(in) :: file_name!
        integer, intent(in) :: string_elements
        integer, dimension(:), intent(in) :: n_crops_by_year
        real(dp), dimension(:,:,:), intent(inout) :: unit_param
        integer, intent(out) :: error_flag!
        integer :: free_unit
        integer :: ios
        integer :: line, p, i
        character(len=string_elements*20) :: buffer, label  !EAC:  use string_elements x 20
        !!
        error_flag = 0!
        call seek_un(error_flag, free_unit)
        open(unit=free_unit, file=trim(file_name), status='old', action="read", iostat=ios )!
        if (ios /= 0 ) then
            print *, "Cannot open file ", trim(file_name), ". The specified file does not exist. &
                & Execution will be aborted..."
            stop
        end if
        
        do while (ios == 0)
            read (free_unit, '(A)', iostat=ios) buffer
            if (ios == 0) then
                line = line + 1
                buffer = trim(buffer)
                call lower_case(buffer)
                p = scan(buffer, achar(9))  ! find the first tab ---> tab=achar(9)!
                label = buffer(1:p-1)
                buffer = buffer(p+1:)
                
                select case (label)
                    case ('year') ! header line
                        i = 1
                    case default
                        call spread_col(buffer, achar(9), string_elements, n_crops_by_year, unit_param(:,:,i))
                        i = i + 1
                end select
            end if
        end do
        close (free_unit)
    end subroutine read_water_prod_file!
    
    subroutine read_canopy_resistance_file(file_name, unit_param, string_elements, error_flag)!
        ! read canopy resistance parameters
        character(len=*), intent(in) :: file_name!
        real(dp), dimension(:), intent(inout) :: unit_param
        integer, intent(in) :: string_elements
        integer, intent(out) :: error_flag!
        integer :: free_unit
        integer :: ios
        integer :: line, p, i
        character(len=string_elements*20) :: buffer, label  ! EAC: string_elements x 20
        !!
        error_flag = 0!
        i = 1
        call seek_un(error_flag, free_unit)
        
        open(unit=free_unit, file=trim(file_name), status='old', action="read", iostat=ios )!
        if (ios /= 0 ) then
            print *, "Cannot open file ", trim(file_name), ". The specified file does not exist. &
                & Execution will be aborted..."
            stop
        end if
        
        read (free_unit, '(A)', iostat=ios)     ! skyp the first line
        do while (ios == 0)
            read (free_unit, '(A)', iostat=ios) buffer
            if (ios == 0) then
                line = line + 1
                buffer = trim(buffer)
                call lower_case(buffer)
                p = scan(buffer, achar(32)) ! find the first space ---> " "=achar(32)! TODO: use always tab?
                label = buffer(1:p-1)
                buffer = buffer(p+1:)
                read (buffer, *) unit_param(i)
                i = i + 1
            end if
        end do
        
        close (free_unit)
    end subroutine read_canopy_resistance_file!
    
    subroutine spread_col_i(string_in, sep, string_el, string_space, string_out)
        implicit none
        character(len=*), intent(in) :: string_in
        character(len=*), intent(in) :: sep
        integer, intent(in) :: string_el
        integer, dimension(:), intent(in) :: string_space
        integer, dimension(:,:), intent(inout) :: string_out   
        integer, dimension(:), allocatable :: dummy
        integer :: i
        
        allocate(dummy(string_el))
        dummy = string_to_integers(string_in(1:len_trim(string_in)-1), sep)
        do i=1, size(string_space)
           string_out(i,1:string_space(i)) = &
                & dummy(sum(string_space(1:i-1))+1:sum(string_space(1:i)))
        end do
    end subroutine spread_col_i
    
    subroutine spread_col_r(string_in, sep, string_el, string_space, string_out)
        implicit none
        character(len=*), intent(in) :: string_in
        character(len=*), intent(in) :: sep
        integer, intent(in) :: string_el
        integer, dimension(:), intent(in) :: string_space
        real(dp), dimension(:,:), intent(inout) :: string_out   
        real(dp), dimension(:), allocatable :: dummy
        integer :: i
        
        allocate(dummy(string_el))
        dummy = string_to_reals(string_in(1:len_trim(string_in)-1), sep)
                       
        do i=1, size(string_space)
           string_out(i,1:string_space(i)) = &
                & dummy(sum(string_space(1:i-1))+1:sum(string_space(1:i)))
        end do
    end subroutine spread_col_r
                                    
    subroutine read_crop_par_file(file_name, string_elements, unit_param, error_flag)!
        ! read static crop parameters file
        ! TODO: merge with init_crop_par_from_file ?
        implicit none
        character(len=*), intent(in) :: file_name!
        integer, intent(in) :: string_elements
        type(crop_pheno_info) :: unit_param
        integer, intent(out) :: error_flag!
        integer :: free_unit
        integer :: ios
        integer :: line, p
        character(len=string_elements*20) :: buffer, label !EAC: use string_elements x 20
        !!
        error_flag = 0!
        call seek_un(error_flag, free_unit)
        open(unit=free_unit, file=trim(file_name), status='old', action="read", iostat=ios )!
        if (ios /= 0 ) then
            print *, "Cannot open file ", trim(file_name), ". The specified file does not exist. &
                & Execution will be aborted..."
            stop
        end if
        
        do while (ios == 0)
            read (free_unit, '(A)', iostat=ios) buffer
            if (ios == 0) then
                line = line + 1
                buffer = trim(buffer)
                call lower_case(buffer)
                p = scan(buffer, achar(9))  ! find the first tab ---> tab=achar(9)!
                label = buffer(1:p-1)
                buffer = buffer(p+1:)

                select case (label)
                    case ('var') ! already initialized
                    case ('irrig')
                        call spread_col(buffer, achar(9), string_elements, unit_param%n_crops_by_year, unit_param%irrigation_class)
                    case ('cnclass')
                        call spread_col(buffer, achar(9), string_elements, unit_param%n_crops_by_year, unit_param%cn_class)
                    case('praw')
                         call spread_col(buffer, achar(9), string_elements, unit_param%n_crops_by_year, unit_param%p_raw_const)
                    case('aint')
                        call spread_col(buffer, achar(9), string_elements, unit_param%n_crops_by_year, unit_param%a)
                    case('tlim')
                        call spread_col(buffer, achar(9), string_elements, unit_param%n_crops_by_year, unit_param%T_lim)
                    case('tcrit')
                        call spread_col(buffer, achar(9), string_elements, unit_param%n_crops_by_year, unit_param%T_crit)
                    case('hi')
                        call spread_col(buffer, achar(9), string_elements, unit_param%n_crops_by_year, unit_param%HI)
                    case('kyt')
                        call spread_col(buffer, achar(9), string_elements, unit_param%n_crops_by_year, unit_param%Ky_tot)
                    case('ky1')
                        call spread_col(buffer, achar(9), string_elements, unit_param%n_crops_by_year, unit_param%Ky_pheno(:,:,1))
                    case('ky2')
                        call spread_col(buffer, achar(9), string_elements, unit_param%n_crops_by_year, unit_param%Ky_pheno(:,:,2))
                    case('ky3')
                        call spread_col(buffer, achar(9), string_elements, unit_param%n_crops_by_year, unit_param%Ky_pheno(:,:,3))
                    case('ky4')
                        call spread_col(buffer, achar(9), string_elements, unit_param%n_crops_by_year, unit_param%Ky_pheno(:,:,4))
                    case('rft')
                         call spread_col(buffer, achar(9), string_elements, unit_param%n_crops_by_year, unit_param%max_RF_t)
                    case default
                        print *, 'Skipping invalid or obsolete label <',trim(label),'> at line', line, &
                            & ' of file: ', file_name
                end select
            end if
        end do
        close (free_unit)
    end subroutine read_crop_par_file!

    subroutine read_crop_pars_r(file_pars,n_days,n_crop)
        ! read crop parameters table from specified unit (real values)
        integer,intent(in) :: n_days                    ! number of days (i.e. 365 o 366)
        integer,intent(in)::n_crop                      ! number of crops
        type(file_phenology_r), intent(inout) :: file_pars!
        integer :: i!
        !!
        allocate(file_pars%tab(n_days,n_crop))!
        do i=1,size(file_pars%tab,1)!
            read(file_pars%unit, *) file_pars%tab(i,:)!
        end do!
        !
    end subroutine read_crop_pars_r!
    !
    subroutine read_crop_pars_i(file_pars,n_days,n_crop)!
    ! read crop parameters table from specified unit (int values)
        integer,intent(in) :: n_days                ! number of days (i.e. 365 o 366)
        integer,intent(in)::n_crop                  ! number of crops
        type(file_phenology_i), intent(inout) :: file_pars!
        integer :: i!
        !!
        allocate(file_pars%tab(n_days,n_crop))!
        do i=1,size(file_pars%tab,1)!
            read(file_pars%unit, *) file_pars%tab(i,:)!
        end do!
        !
    end subroutine read_crop_pars_i!

    subroutine init_crop_phenology_pars(sim,info_pheno,info_meteo, debug)!
        ! init crop parameters and file references for daily parameters
        type(simulation),intent(inout)::sim!
        type(meteo_info),dimension(:),intent(in)::info_meteo!
        logical,intent(in)::debug

        type(crop_pheno_info),dimension(:),allocatable::info_pheno!
        character(len=255)::dir,froot,dir_name!
        integer::i,errorflag,string_elements!
        integer, dimension(sim%n_lus) :: n_crops_by_year
        real(dp), parameter :: nan = -9999.0D0
        integer, parameter :: phases = 4
        !!
        dir= trim(sim%pheno_path)
        froot = sim%pheno_root
        !!
        allocate(info_pheno(sim%n_voronoi)) ! init to the number of weather stations
        
        allocate(sim%res_canopy(sim%meteo_years))!
        sim%res_canopy=0!
        ! init from crop parameters file (actually produced by cropcoeff)
        dir_name = info_meteo(1)%filename(1:(index(trim(info_meteo(1)%filename),"."))-1)  ! directory has the same name as the weather station dataset
        call init_crop_par_from_file(trim(dir)//trim(froot)//trim(dir_name)//delimiter//"CropParam.dat", &
            & sim%n_lus, sim%n_crops, string_elements, n_crops_by_year, ErrorFlag)
            
        call read_canopy_resistance_file(trim(dir)//delimiter//'CanopyRes.dat', sim%res_canopy, string_elements, ErrorFlag)
        
        do i=1,size(info_pheno)!
            dir_name = info_meteo(i)%filename(1:(index(trim(info_meteo(i)%filename),"."))-1)
            call open_daily_crop_par_file(info_pheno(i)%k_cb%unit,trim(dir)//trim(froot)//trim(dir_name)//delimiter//"Kcb.dat",errorflag)!
            call open_daily_crop_par_file(info_pheno(i)%h%unit,trim(dir)//trim(froot)//trim(dir_name)//delimiter//"H.dat",errorflag)!         
            call open_daily_crop_par_file(info_pheno(i)%d_r%unit,trim(dir)//trim(froot)//trim(dir_name)//delimiter//"Sr.dat",errorflag)!
            call open_daily_crop_par_file(info_pheno(i)%lai%unit,trim(dir)//trim(froot)//trim(dir_name)//delimiter//"LAI.dat",errorflag)!
            call open_daily_crop_par_file(info_pheno(i)%cn_day%unit,trim(dir)//trim(froot)//trim(dir_name)//delimiter//"CNvalue.dat",errorflag)!
            call open_daily_crop_par_file(info_pheno(i)%f_c%unit,trim(dir)//trim(froot)//trim(dir_name)//delimiter//"fc.dat",errorflag)
            ! EDIT: add support for seasonal p_raw
            call open_daily_crop_par_file(info_pheno(i)%r_stress%unit,trim(dir)//trim(froot)//trim(dir_name)//delimiter//"r_stress.dat",errorflag)
            !
            ! TODO - add tabulated ky
            
            ! init crop parameters
            allocate(info_pheno(i)%n_crops_by_year    (sim%n_lus))
            allocate(info_pheno(i)%irrigation_class (sim%n_lus, sim%n_crops))
            allocate(info_pheno(i)%cn_class             (sim%n_lus, sim%n_crops))
            allocate(info_pheno(i)%p_raw_const     (sim%n_lus, sim%n_crops))
            allocate(info_pheno(i)%a              (sim%n_lus, sim%n_crops))
            allocate(info_pheno(i)%max_d_r          (sim%n_lus, sim%n_crops))
            allocate(info_pheno(i)%max_RF_t         (sim%n_lus, sim%n_crops))
            allocate(info_pheno(i)%T_lim           (sim%n_lus, sim%n_crops))
            allocate(info_pheno(i)%T_crit          (sim%n_lus, sim%n_crops))
            allocate(info_pheno(i)%HI             (sim%n_lus, sim%n_crops))
            allocate(info_pheno(i)%Ky_tot            (sim%n_lus, sim%n_crops))
            allocate(info_pheno(i)%Ky_pheno            (sim%n_lus, sim%n_crops, phases))
            allocate(info_pheno(i)%wp_adj          (sim%n_lus, sim%n_crops, sim%meteo_years))
            allocate(info_pheno(i)%kcb_phases%low (sim%n_lus, sim%n_crops))
            allocate(info_pheno(i)%kcb_phases%high(sim%n_lus, sim%n_crops))
            allocate(info_pheno(i)%kcb_phases%mid (sim%n_lus, sim%n_crops))
            allocate(info_pheno(i)%ii0            (sim%n_lus, sim%n_crops))
            allocate(info_pheno(i)%iie            (sim%n_lus, sim%n_crops))
            allocate(info_pheno(i)%iid            (sim%n_lus, sim%n_crops))
            info_pheno(i)%n_crops_by_year     = n_crops_by_year
            info_pheno(i)%irrigation_class  = int(nan)
            info_pheno(i)%cn_class              = int(nan)
            info_pheno(i)%p_raw_const               = nan
            info_pheno(i)%a               = nan
            info_pheno(i)%max_d_r           = nan
            info_pheno(i)%max_RF_t          = nan
            info_pheno(i)%T_lim            = nan
            info_pheno(i)%T_crit           = nan
            info_pheno(i)%HI              = nan
            info_pheno(i)%Ky_tot             = nan
            info_pheno(i)%Ky_pheno             = nan
            info_pheno(i)%wp_adj           = nan
            info_pheno(i)%kcb_phases%low  = nan
            info_pheno(i)%kcb_phases%high = nan
            info_pheno(i)%kcb_phases%mid  = nan
            info_pheno(i)%ii0             = 0!
            info_pheno(i)%iie             = 0!
            info_pheno(i)%iid             = 0!
            call read_water_prod_file(trim(dir)//trim(froot)//trim(dir_name)//delimiter//"WPadj.dat", &
                & string_elements, n_crops_by_year, info_pheno(i)%wp_adj, ErrorFlag)
            call read_crop_par_file(trim(dir)//trim(froot)//trim(dir_name)//delimiter//"CropParam.dat", &
                & string_elements, info_pheno(i), ErrorFlag)
            
            ! TODO
            ! EAC: overwrite p values if exits
            ! call open_daily_crop_par_file(info_pheno(i)%p%unit,trim(dir)//trim(froot)//trim(fname)//"\praw.dat",errorflag)
            
        end do!
        if (debug .eqv. .true.) then
            print *,'===== DEBUG: crop parameters initialize ====='
            print *,'path to phenophase: ',  dir
            print *,'root of subfolder: ',  froot
            print *,'# of phenophase: ',  size(info_meteo)
            print *,'===== END DEBUG ====='
        end if
    end subroutine init_crop_phenology_pars!

    subroutine read_all_crop_pars(n_days,n_crop,info_pheno)!
        ! read all daily crop parameters by the number of days over the opened files
        implicit none!
        integer,intent(in)::n_days!
        integer,intent(in)::n_crop!
        type(crop_pheno_info),dimension(:),intent(inout)::info_pheno!
        integer:: i, doy!
        integer:: crop_idx ! current crop index
        integer:: m, cs, ii_idx, ij_idx
        logical, dimension(n_days)::dummy_array
        !!
        do i=1,size(info_pheno)!
            call read_crop_pars(info_pheno(i)%k_cb,n_days,n_crop)!
            call read_crop_pars(info_pheno(i)%h,n_days,n_crop)!
            call read_crop_pars(info_pheno(i)%d_r,n_days,n_crop)!
            call read_crop_pars(info_pheno(i)%lai,n_days,n_crop)!
            call read_crop_pars(info_pheno(i)%cn_day,n_days,n_crop)!
            call read_crop_pars(info_pheno(i)%f_c,n_days,n_crop)! %RR%
            call read_crop_pars(info_pheno(i)%r_stress,n_days,n_crop)! %EAC%
        
            ! find minimum and maximum k_cb value to define stages
            info_pheno(i)%kcb_phases%low(:,1) = minval(info_pheno(i)%k_cb%tab, dim=1) ! low not to be re-calculated
            info_pheno(i)%kcb_phases%high(:,1) = maxval(info_pheno(i)%k_cb%tab, dim=1)
            
            do crop_idx = 1, size(info_pheno(i)%k_cb%tab, 2)
                if (info_pheno(i)%n_crops_by_year(crop_idx) > 1) then   ! if there are more crops, populate also the other row
                    ij_idx = 1 ! init the index to the beginning of the serie
                    do cs = 1, info_pheno(i)%n_crops_by_year(crop_idx)
                        ! low
                        info_pheno(i)%kcb_phases%low(crop_idx,cs) = info_pheno(i)%kcb_phases%low(crop_idx,1) ! low not to be recalculated
                        ! definition of stages
                        do ii_idx = ij_idx, size(info_pheno(i)%k_cb%tab, 1)
                            if (info_pheno(i)%k_cb%tab(ii_idx,crop_idx) > 0) exit ! find the first value > 0 (start crop)
                        end do
                        do ij_idx = ii_idx, size(info_pheno(i)%k_cb%tab, 1)
                            if (info_pheno(i)%k_cb%tab(ij_idx,crop_idx) == 0) exit ! find the following first value == 0 (end crop)
                        end do
                        ! high
                        dummy_array = .false.
                        ! limit to the maximum number of days in the year
                        if (ij_idx > n_days) ij_idx = n_days
                        dummy_array(ii_idx:ij_idx) = .true.
                        info_pheno(i)%kcb_phases%high(crop_idx,cs) = &
                            & maxval(info_pheno(i)%k_cb%tab(:,crop_idx), dim=1, mask=dummy_array) ! select the greatest between ii_idx and ij_idx
                        ! mid
                        do m = ii_idx+1, ij_idx
                            if (info_pheno(i)%k_cb%tab(m,crop_idx) /= info_pheno(i)%kcb_phases%low(crop_idx,1) &
                                & .and. info_pheno(i)%k_cb%tab(m,crop_idx) /= info_pheno(i)%kcb_phases%high(crop_idx,cs) &
                                & .and. info_pheno(i)%k_cb%tab(m,crop_idx) == info_pheno(i)%k_cb%tab(m-1,crop_idx)) then
                                info_pheno(i)%kcb_phases%mid(crop_idx,cs) = info_pheno(i)%k_cb%tab(m,crop_idx) ! find the mid value between ii_idx and ij_idx
                                exit
                            end if
                        end do
                    end do
                else
                    do m = 2, size(info_pheno(i)%k_cb%tab, 1)
                        if (info_pheno(i)%k_cb%tab(m,crop_idx) /= info_pheno(i)%kcb_phases%low(crop_idx,1) &
                            & .and. info_pheno(i)%k_cb%tab(m,crop_idx) /= info_pheno(i)%kcb_phases%high(crop_idx,1) &
                            & .and. info_pheno(i)%k_cb%tab(m,crop_idx) == info_pheno(i)%k_cb%tab(m-1,crop_idx)) then
                            info_pheno(i)%kcb_phases%mid(crop_idx,1) = info_pheno(i)%k_cb%tab(m,crop_idx)
                            exit
                        end if
                    end do
                end if
            end do
        end do!
        
        do i=1,size(info_pheno)
            ! Calculate max rooting depth for annual crops
            info_pheno(i)%max_d_r(:,1) = maxval(info_pheno(i)%d_r%tab, dim=1)
            
            do crop_idx = 1, size(info_pheno(i)%d_r%tab, 2)
                if (info_pheno(i)%n_crops_by_year(crop_idx) > 1) then
                    ij_idx = 1 ! index initialization
                    do cs = 1, info_pheno(i)%n_crops_by_year(crop_idx)
                        do ii_idx = ij_idx, size(info_pheno(i)%d_r%tab, 1)
                            if (info_pheno(i)%d_r%tab(ii_idx,crop_idx) > 0) exit
                        end do
                        do ij_idx = ii_idx, size(info_pheno(i)%d_r%tab, 1)
                            if (info_pheno(i)%d_r%tab(ij_idx,crop_idx) == 0) exit
                        end do
                        dummy_array = .false.
                        if (ij_idx > n_days) ij_idx = n_days
                        dummy_array(ii_idx:ij_idx) = .true.
                        info_pheno(i)%max_d_r(crop_idx,cs) = maxval(info_pheno(i)%d_r%tab(:,crop_idx), dim=1, mask=dummy_array)
                    end do
                end if
            end do
        end do
        
        do i=1,size(info_pheno)  ! loop over phenological dataset
            do crop_idx=1,size(info_pheno(i)%k_cb%tab,2)  ! loop over crops
                if (info_pheno(i)%n_crops_by_year(crop_idx) == 1) then   ! permanent or annual single crops  (also winter wheat)
                    if (info_pheno(i)%kcb_phases%low(crop_idx,1) == 0) then                ! not permanent crops
                        ! calculate the length of the crop cycle
                        info_pheno(i)%iid(crop_idx,1)= &
                            & count(info_pheno(i)%k_cb%tab(:,crop_idx)>info_pheno(i)%kcb_phases%low(crop_idx,1)) + 1
                        
                        ! calculate the emergence day of the crop
                        do doy=1,size(info_pheno(i)%k_cb%tab,1)        ! loop over days 
                            if(info_pheno(i)%k_cb%tab(doy,crop_idx) >  info_pheno(i)%kcb_phases%low(crop_idx,1)) then!
                                info_pheno(i)%ii0(crop_idx,1)=doy       ! ii0: crop emergence doy
                                ! ii0 = 0 in case of winter crop
                                exit
                            end if
                        end do
                        
                        ! calculate the harvest day
                        do doy=n_days-1,1,-1 ! loop over days starting from the last
                            if(info_pheno(i)%k_cb%tab(doy,crop_idx) >  info_pheno(i)%kcb_phases%low(crop_idx,1)) then!
                                info_pheno(i)%iie(crop_idx,1)=doy       ! iee: crop harvest doy
                                exit
                            end if
                        end do
                        
                        if (info_pheno(i)%ii0(crop_idx,1)==1) then                         ! winter crop (ii0==1, TODO: ?)
                            do doy=1,size(info_pheno(i)%k_cb%tab,1)    ! loop over days
                                if(info_pheno(i)%k_cb%tab(doy,crop_idx) ==  0.)then!
                                    info_pheno(i)%iie(crop_idx,1)=doy   ! iie: harvest doy = 1st day with k_cb = 0 after the crop period
                                    info_pheno(i)%ii0(crop_idx,1)= &
                                        & n_days - info_pheno(i)%iid(crop_idx,1) + info_pheno(i)%iie(crop_idx,1) + 1 ! update ii0
                                    exit
                                end if
                            end do
                        end if
                    else ! permanent crops                                                                
                        ! calculate the emergence day
                        do doy=1,size(info_pheno(i)%k_cb%tab,1)    ! loop over days
                            if (info_pheno(i)%k_cb%tab(doy,crop_idx) > info_pheno(i)%kcb_phases%low(crop_idx,1)) then
                                info_pheno(i)%ii0(crop_idx,1)=doy   ! ii0: day of emergence or the end of vernalization
                                exit
                            end if
                        end do
                        
                        ! calculate harvest day
                        do doy=n_days,1,-1 ! loop over days starting from the last
                            if (info_pheno(i)%k_cb%tab(doy,crop_idx) > info_pheno(i)%kcb_phases%low(crop_idx,1)) then
                                info_pheno(i)%iie(crop_idx,1)=doy   ! iie: day of harvest/ beginning of vernalization
                                exit
                            end if
                        end do
                        
                        ! calculate the crop period
                        info_pheno(i)%iid(crop_idx,1) = info_pheno(i)%iie(crop_idx,1) - info_pheno(i)%ii0(crop_idx,1) + 1
                    end if
                else   ! alternate crops
                    ii_idx = 1 ! init the day index at the beginning of the year
                    do cs = 1, info_pheno(i)%n_crops_by_year(crop_idx)
                        ! update the day index for the alternative crop
                        if (cs > 1) then
                            ii_idx = info_pheno(i)%iie(crop_idx,cs-1) + 1 
                        end if
                        ! calculate emergence day
                        do doy = ii_idx,size(info_pheno(i)%k_cb%tab,1)    ! loop over the days
                            if(info_pheno(i)%k_cb%tab(doy,crop_idx) > info_pheno(i)%kcb_phases%low(crop_idx,cs)) then   ! use cs as reference
                                info_pheno(i)%ii0(crop_idx,cs)=doy
                                exit
                            end if
                        end do
                        
                        ! calculate harvest day
                        do doy = ii_idx, size(info_pheno(i)%k_cb%tab,1)
                            if (info_pheno(i)%k_cb%tab(doy,crop_idx) == info_pheno(i)%kcb_phases%high(crop_idx,cs)) then ! find the first "high" k_cb
                                ii_idx = doy
                                exit
                            end if
                        end do
                        do doy = ii_idx,size(info_pheno(i)%k_cb%tab,1)    ! loop over the following days
                            if(info_pheno(i)%k_cb%tab(doy,crop_idx) == info_pheno(i)%kcb_phases%low(crop_idx,cs)) then   ! the first "low" after the first "high" is the harvest date
                                info_pheno(i)%iie(crop_idx,cs)=doy-1
                                exit
                            end if
                        end do
                        
                        ! calculate the duration of the crop cycle
                        info_pheno(i)%iid(crop_idx,cs)=info_pheno(i)%iie(crop_idx,cs) - info_pheno(i)%ii0(crop_idx,cs) + 1
                    end do
                    
                    ! adjust the duration of the crop cycle and the emergence day of the first crop in case ii0==1
                    if (info_pheno(i)%ii0(crop_idx,1) == 1) then
                        ii_idx = count(info_pheno(i)%k_cb%tab(:,crop_idx)>info_pheno(i)%kcb_phases%low(crop_idx,1)) ! sum the duration of the crop cycles
                        ii_idx = ii_idx - sum(info_pheno(i)%iid(crop_idx,:)) + info_pheno(i)%iid(crop_idx,1)        ! duration of the first crop cycle
                        info_pheno(i)%ii0(crop_idx,1) = n_days - ii_idx + info_pheno(i)%iid(crop_idx,1) + 3         ! emergence of the first crop (consider also n_days and xii_idx)
                        info_pheno(i)%iid(crop_idx,1) = ii_idx - 1
                    end if
                end if
            end do
        end do!
    end subroutine read_all_crop_pars!
    !
    subroutine destroy_infofeno_tab(info_pheno)!
    ! dellaocate all crop phenological time series
        type(crop_pheno_info),dimension(:),intent(inout)::info_pheno!
        integer::i!
        !!
        do i=1,size(info_pheno)!
            if(associated(info_pheno(i)%k_cb%tab)) deallocate(info_pheno(i)%k_cb%tab)!
            if(associated(info_pheno(i)%h%tab)) deallocate(info_pheno(i)%h%tab)!
            if(associated(info_pheno(i)%d_r%tab)) deallocate(info_pheno(i)%d_r%tab)!
            if(associated(info_pheno(i)%lai%tab)) deallocate(info_pheno(i)%lai%tab)!
            if(associated(info_pheno(i)%cn_day%tab)) deallocate(info_pheno(i)%cn_day%tab)!
            if(associated(info_pheno(i)%f_c%tab)) deallocate(info_pheno(i)%f_c%tab)!
            if(associated(info_pheno(i)%r_stress%tab)) deallocate(info_pheno(i)%r_stress%tab)!
        end do!
        
    end subroutine destroy_infofeno_tab!
    !
    subroutine close_pheno_file(info_pheno)!
        ! close all phenological opened files
        type(crop_pheno_info),dimension(:),allocatable,intent(inout)::info_pheno!
        integer::i!
        !!
        do i=1,size(info_pheno)!
            close(info_pheno(i)%k_cb%unit)!
            close(info_pheno(i)%h%unit)!
            close(info_pheno(i)%d_r%unit)!
            close(info_pheno(i)%lai%unit)!
            close(info_pheno(i)%cn_day%unit)!
            close(info_pheno(i)%f_c%unit)!
            close(info_pheno(i)%r_stress%unit)!
        end do!
        deallocate(info_pheno)!
    end subroutine close_pheno_file!
    !
    subroutine check_pheno_parameters(info_pheno,info_meteo)!
    ! check if phenological parameters match weather station data
        implicit none!
        type(crop_pheno_info),dimension(:),intent(in)::info_pheno!
        type(meteo_info),dimension(:),intent(in)::info_meteo!
        integer::i!
        !
        do i=1,size(info_pheno)!
            call check_crop_parameters(info_pheno(i),info_meteo(i)%filename)!
        end do!
    end subroutine check_pheno_parameters!

    subroutine check_crop_parameters(info_pheno,weather_station)!
    ! check if phenological parameters match weather station data
    ! if k_cb is null than all the other parameters must be null
        implicit none
        type(crop_pheno_info),intent(in)::info_pheno!
        character(len=*),intent(in)::weather_station!
        integer::error_flag!
        integer::d,k!
        !!
        error_flag = 0!
        !!
        write(*,*)'INPUT CHECK - SOIL USE PARAMETERS - STATION:',trim(weather_station)!
        do k=1, size(info_pheno%k_cb%tab,2) ! loop over crops
            do d=1,size(info_pheno%k_cb%tab,1)    ! loop over data
                if(info_pheno%k_cb%tab(d,k).gt.0)then!
                    if(info_pheno%d_r%tab(d,k).eq.0.) then!
                        write(*,*)' Warning: Sr null. Day: ',d,' Soil use class: ', k !
                        error_flag = -1!
                    end if!
                    if(info_pheno%h%tab(d,k).eq.0.) then!
                        write(*,*)' Warning: H null. Day: ',d,'  Soil use class:', k!
                        error_flag = -1!
                    end if!
                    if(info_pheno%LAI%tab(d,k).eq.0.) then!
                        write(*,*)' Warning: LAI null. Day: ',d,'  Soil use class:', k!
                        error_flag = -1!
                    end if!
                end if!
            end do!
        end do!
        write(*,*)'END INPUT CHECK - SOIL USE PARAMETERS - STATION:',trim(weather_station)!
        !
    end subroutine check_crop_parameters!
    !
end module cli_crop_parameters!
