module cli_crop_parameters!
    use mod_constants, only: sp, dp
    use mod_utility, only: seek_un, lower_case, string_to_integers, string_to_reals, split_string, count_element
    use mod_parameters, only: simulation,parameters
    use mod_meteo, only: meteo_info
    use mod_crop_phenology
    use mod_system
    
    implicit none!

    logical, dimension(:), allocatable, save :: missing_crop_slot_warned ! %PS%

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
                & The CropCoef version you used to generate inputs might be outdated. Execution will be aborted..."
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
                        print *, 'Skipping invalid or obsolete label <',trim(label),'> at line', line, ' of file: ', file_name
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

    subroutine init_crop_phenology_pars(sim,info_pheno,info_meteo, verbose)!
        ! init crop parameters and file references for daily parameters
        type(simulation),intent(inout)::sim!
        type(meteo_info),dimension(:),intent(in)::info_meteo!
        logical,intent(in)::verbose

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
            call open_daily_crop_par_file(info_pheno(i)%z_r%unit,trim(dir)//trim(froot)//trim(dir_name)//delimiter//"Sr.dat",errorflag)!
            call open_daily_crop_par_file(info_pheno(i)%lai%unit,trim(dir)//trim(froot)//trim(dir_name)//delimiter//"LAI.dat",errorflag)!
            call open_daily_crop_par_file(info_pheno(i)%cn_day%unit,trim(dir)//trim(froot)//trim(dir_name)//delimiter//"CNvalue.dat",errorflag)!
            call open_daily_crop_par_file(info_pheno(i)%f_c%unit,trim(dir)//trim(froot)//trim(dir_name)//delimiter//"fc.dat",errorflag)
            ! EDIT: add support for seasonal p_raw
            call open_daily_crop_par_file(info_pheno(i)%r_stress%unit,trim(dir)//trim(froot)//trim(dir_name)//delimiter//"r_stress.dat",errorflag)
            call open_daily_crop_par_file(info_pheno(i)%crop_id%unit,trim(dir)//trim(froot)//trim(dir_name)//delimiter//"CropId.dat",errorflag)
            !
            ! TODO - add tabulated ky
            
            ! init crop parameters
            allocate(info_pheno(i)%n_crops_by_year    (sim%n_lus))
            allocate(info_pheno(i)%irrigation_class (sim%n_lus, sim%n_crops))
            allocate(info_pheno(i)%cn_class             (sim%n_lus, sim%n_crops))
            allocate(info_pheno(i)%p_raw_const     (sim%n_lus, sim%n_crops))
            allocate(info_pheno(i)%a              (sim%n_lus, sim%n_crops))
            allocate(info_pheno(i)%d_r_max          (sim%n_lus, sim%n_crops))
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
            allocate(info_pheno(i)%cycle_crop_slot(sim%n_lus, sim%n_crops))
            info_pheno(i)%n_crops_by_year     = n_crops_by_year
            info_pheno(i)%irrigation_class  = int(nan)
            info_pheno(i)%cn_class              = int(nan)
            info_pheno(i)%p_raw_const               = nan
            info_pheno(i)%a               = nan
            info_pheno(i)%d_r_max           = nan
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
            info_pheno(i)%cycle_crop_slot = 0
            call read_water_prod_file(trim(dir)//trim(froot)//trim(dir_name)//delimiter//"WPadj.dat", &
                & string_elements, n_crops_by_year, info_pheno(i)%wp_adj, ErrorFlag)
            call read_crop_par_file(trim(dir)//trim(froot)//trim(dir_name)//delimiter//"CropParam.dat", &
                & string_elements, info_pheno(i), ErrorFlag)
            
            ! TODO
            ! EAC: overwrite p values if exits
            ! call open_daily_crop_par_file(info_pheno(i)%p%unit,trim(dir)//trim(froot)//trim(fname)//"\praw.dat",errorflag)
            
        end do!
        if (verbose .eqv. .true.) then
            print *,'===== DEBUG: crop parameters initialize ====='
            print *,'path to phenophase: ',  dir
            print *,'root of subfolder: ',  froot
            print *,'# of phenophase: ',  size(info_meteo)
            print *,'===== END DEBUG ====='
        end if
    end subroutine init_crop_phenology_pars!

    subroutine read_all_crop_pars(n_days,n_crop,info_pheno,pars)!
        ! %PS% Read one calendar year and derive crop cycles exclusively from CropId.
        implicit none!
        integer,intent(in)::n_days,n_crop
        type(parameters),intent(in)::pars
        type(crop_pheno_info),dimension(:),intent(inout)::info_pheno
        integer::i

        do i=1,size(info_pheno)
            call read_crop_pars(info_pheno(i)%k_cb,n_days,n_crop)
            call read_crop_pars(info_pheno(i)%h,n_days,n_crop)
            call read_crop_pars(info_pheno(i)%z_r,n_days,n_crop)
            call read_crop_pars(info_pheno(i)%lai,n_days,n_crop)
            call read_crop_pars(info_pheno(i)%cn_day,n_days,n_crop)
            call read_crop_pars(info_pheno(i)%f_c,n_days,n_crop)
            call read_crop_pars(info_pheno(i)%r_stress,n_days,n_crop)
            call read_crop_pars(info_pheno(i)%crop_id,n_days,n_crop)
            call derive_crop_cycles(info_pheno(i), n_days, pars%depth%ze_fix)
        end do
    end subroutine read_all_crop_pars

    subroutine derive_crop_cycles(pheno, n_days, ze_fix)
        type(crop_pheno_info), intent(inout) :: pheno
        integer, intent(in) :: n_days
        real(dp), intent(in) :: ze_fix
        integer :: lu, slot, cycle_idx, day, start_day, end_day, crop_slot, first_end, last_start
        integer :: n_slots, n_present_slots
        logical, dimension(n_days) :: crop_mask
        real(dp) :: low_value, high_value, mid_value

        pheno%ii0 = 0
        pheno%iie = 0
        pheno%iid = 0
        pheno%cycle_crop_slot = 0

        if (.not. allocated(missing_crop_slot_warned)) then
            allocate(missing_crop_slot_warned(size(pheno%crop_id%tab,2)))
            missing_crop_slot_warned = .false.
        end if

        do lu=1,size(pheno%crop_id%tab,2)
            n_slots = pheno%n_crops_by_year(lu)
            if (n_slots < 1) cycle

            do day=1,n_days
                crop_slot = pheno%crop_id%tab(day,lu)
                if (crop_slot < 0 .or. crop_slot > n_slots) then
                    print *, 'Invalid crop slot ', crop_slot, ' at day ', day, ', land-use class ', lu
                    print *, 'Expected a value between 0 and ', n_slots
                    print *, 'Execution will be aborted...'
                    stop
                end if
            end do

            ! Derive crop-specific daily-series parameters by rotation slot.
            n_present_slots = 0
            do slot=1,n_slots
                crop_mask = pheno%crop_id%tab(:,lu) == slot
                if (.not. any(crop_mask)) then
                    if (.not. missing_crop_slot_warned(lu)) then
                        print *, 'Warning: Crop slot ', slot, ' never occurs in land-use class ', lu, &
                                &'. CropCoef likely overwrote an entire crop due to overlapping sow/harvest dates.'
                        missing_crop_slot_warned(lu) = .true.
                    end if
                    cycle
                end if
                n_present_slots = n_present_slots + 1
                ! %PS% Preserve the existing annual/permanent convention: annual
                ! land uses have a zero Kcb outside crop-in-field periods.
                low_value = minval(pheno%k_cb%tab(:,lu))
                high_value = maxval(pheno%k_cb%tab(:,lu), mask=crop_mask)
                mid_value = high_value
                do day=2,n_days
                    if (crop_mask(day) .and. crop_mask(day-1)) then
                        if (pheno%k_cb%tab(day,lu) == pheno%k_cb%tab(day-1,lu) .and. &
                            & pheno%k_cb%tab(day,lu) > low_value .and. pheno%k_cb%tab(day,lu) < high_value) then
                            mid_value = pheno%k_cb%tab(day,lu)
                            exit
                        end if
                    end if
                end do
                pheno%kcb_phases%low(lu,slot) = low_value
                pheno%kcb_phases%high(lu,slot) = high_value
                pheno%kcb_phases%mid(lu,slot) = mid_value
                pheno%d_r_max(lu,slot) = maxval(pheno%z_r%tab(:,lu), mask=crop_mask) - ze_fix
            end do

            cycle_idx = 0
            first_end = 0
            last_start = n_days + 1

            if (pheno%crop_id%tab(1,lu) > 0 .and. &
                & pheno%crop_id%tab(1,lu) == pheno%crop_id%tab(n_days,lu)) then
                crop_slot = pheno%crop_id%tab(1,lu)
                first_end = 1
                do while (first_end < n_days .and. pheno%crop_id%tab(first_end+1,lu) == crop_slot)
                    first_end = first_end + 1
                end do
                cycle_idx = 1
                if (first_end == n_days) then
                    last_start = n_days + 1
                    pheno%ii0(lu,cycle_idx) = 1
                    pheno%iie(lu,cycle_idx) = n_days
                    pheno%iid(lu,cycle_idx) = n_days
                else
                    last_start = n_days
                    do while (last_start > 1 .and. pheno%crop_id%tab(last_start-1,lu) == crop_slot)
                        last_start = last_start - 1
                    end do
                    pheno%ii0(lu,cycle_idx) = last_start
                    pheno%iie(lu,cycle_idx) = first_end
                    pheno%iid(lu,cycle_idx) = n_days-last_start+1+first_end
                end if
                pheno%cycle_crop_slot(lu,cycle_idx) = crop_slot
            end if

            day = first_end + 1
            do while (day <= min(n_days,last_start-1))
                if (pheno%crop_id%tab(day,lu) == 0) then
                    day = day + 1
                    cycle
                end if
                crop_slot = pheno%crop_id%tab(day,lu)
                start_day = day
                do while (day <= min(n_days,last_start-1) .and. pheno%crop_id%tab(day,lu) == crop_slot)
                    day = day + 1
                end do
                end_day = day - 1
                cycle_idx = cycle_idx + 1
                if (cycle_idx > size(pheno%ii0,2)) then
                    print *, 'Too many crop cycles in land-use class ', lu
                    stop
                end if
                pheno%ii0(lu,cycle_idx) = start_day
                pheno%iie(lu,cycle_idx) = end_day
                pheno%iid(lu,cycle_idx) = end_day-start_day+1
                pheno%cycle_crop_slot(lu,cycle_idx) = crop_slot
            end do

            if (cycle_idx /= n_present_slots) then
                print *, 'CropId.dat defines ', cycle_idx, ' crop cycles for land-use class ', lu, &
                    & ', but ', n_present_slots, ' declared slots occur in the daily series.'
                print *, 'Execution will be aborted...'
                stop
            end if
        end do
    end subroutine derive_crop_cycles

    subroutine destroy_infofeno_tab(info_pheno)!
    ! dellaocate all crop phenological time series
        type(crop_pheno_info),dimension(:),intent(inout)::info_pheno!
        integer::i!
        !!
        do i=1,size(info_pheno)!
            if(associated(info_pheno(i)%k_cb%tab)) deallocate(info_pheno(i)%k_cb%tab)!
            if(associated(info_pheno(i)%h%tab)) deallocate(info_pheno(i)%h%tab)!
            if(associated(info_pheno(i)%z_r%tab)) deallocate(info_pheno(i)%z_r%tab)!
            if(associated(info_pheno(i)%lai%tab)) deallocate(info_pheno(i)%lai%tab)!
            if(associated(info_pheno(i)%cn_day%tab)) deallocate(info_pheno(i)%cn_day%tab)!
            if(associated(info_pheno(i)%f_c%tab)) deallocate(info_pheno(i)%f_c%tab)!
            if(associated(info_pheno(i)%r_stress%tab)) deallocate(info_pheno(i)%r_stress%tab)!
            if(associated(info_pheno(i)%crop_id%tab)) deallocate(info_pheno(i)%crop_id%tab)! %PS%
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
            close(info_pheno(i)%z_r%unit)!
            close(info_pheno(i)%lai%unit)!
            close(info_pheno(i)%cn_day%unit)!
            close(info_pheno(i)%f_c%unit)!
            close(info_pheno(i)%r_stress%unit)!
            close(info_pheno(i)%crop_id%unit)
            if(associated(info_pheno(i)%cycle_crop_slot)) deallocate(info_pheno(i)%cycle_crop_slot)
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
                    if(info_pheno%z_r%tab(d,k).eq.0.) then!
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
