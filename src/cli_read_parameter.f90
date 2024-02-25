module cli_read_parameter!
    use mod_parameters
    use mod_common
    use mod_grid
    use mod_irrigation, only: calc_perc_booster_pars
    use cli_watsources
    use mod_TDx_index

    use cli_save_outputs
    implicit none!

    contains

    subroutine read_all_parameters(file_xml, xml, xml_dtx, ErrorFlag, debug)
        implicit none
        character(len=*), intent(in) :: file_xml
        logical, intent(in) :: debug
        integer, intent(out) :: ErrorFlag
        type(TDx_index), intent(inout) :: xml_dtx
        type(parameters), intent(inout) :: xml
        ! Input related variables used in parsing loop
        integer :: ios          ! State variable (0 = ok)
        integer :: line, tablestart
        logical :: file_exists

        ErrorFlag = 0
        line = 0; tablestart = 0
        ios = 0

        call make_default(xml, xml_dtx)

        call read_sim_parameters(file_xml, xml, xml_dtx, ErrorFlag,debug)
        
        inquire(file="cells.txt", exist=xml%sim%f_out_cells) ! update "output_cells"
        
        ! fw - evaporation !%AB% used for not irrigated cells or without irrigation events
        xml%irr%f_w = cost_fwEva
        ! allocate precipitation over the 24 hours
        allocate(xml%f_eff_rain(24))
        ! %AB% used only for evaluate effective precipitation
        xml%f_eff_rain = cost_f_eff_rain
        
        ! allocate ET0 over the 24 hours
        allocate(xml%fet0(24))
        xml%fet0 = cost_fet0
        xml%ms_i%f_exists=.false. ; xml%ms_i%f_allocate_timeserie=.false. ;!
        xml%ms_ii%f_exists=.false. ; xml%ms_ii%f_allocate_timeserie=.false. ;!
        xml%intreu%f_exists=.false. ; xml%intreu%f_allocate_timeserie=.false. ;!
        xml%cr%f_exists=.false. ; xml%cr%f_allocate_timeserie=.false. ;!
        
        !!! %AB% check for time series of water releases
        inquire(file=trim(xml%sim%watsour_path)//trim(xml%sim%mon_sources_i_div_fn), exist=file_exists)    ! file_exists will be TRUE if the file exists
        if (file_exists .eqv. .true.) then
            xml%ms_i%f_exists=.true.!
            xml%ms_i%f_allocate_timeserie=.true.!
        end if!
        
        inquire(file=trim(xml%sim%watsour_path)//trim(xml%sim%mon_sources_ii_div_fn), exist=file_exists)    ! file_exists will be TRUE if the file exists
        if (file_exists .eqv. .true.) then
            xml%ms_ii%f_exists=.true.!
            xml%ms_ii%f_allocate_timeserie=.true.!
        end if!
        inquire(file=trim(xml%sim%watsour_path)//trim(xml%sim%int_reuse_div_fn), exist=file_exists)    ! file_exists will be TRUE if the file exists
        if (file_exists .eqv. .true.) then
            xml%intreu%f_exists=.true.!
            xml%intreu%f_allocate_timeserie=.true.!
        end if!
        
        inquire(file=trim(xml%sim%watsour_path)//trim(xml%sim%cr_sources_list_fn), exist=file_exists)    ! file_exists will be TRUE if the file exists
        if (file_exists .eqv. .true.) then
            xml%cr%f_exists=.true.!
            ! Collective unmonitored sources have no time series of discharges
        end if!
        
        call read_all_irr_methods(xml, ErrorFlag, debug)

    end subroutine read_all_parameters


    subroutine read_sim_parameters(file_xml, xml, xml_dtx, ErrorFlag,debug)
        ! read settings for the simulation
        implicit none
        character(len=*), intent(in) :: file_xml
        integer, intent(out) :: ErrorFlag
        logical, intent(in) :: debug
        type(TDx_index), intent(inout) :: xml_dtx
        type(parameters), intent(INout) :: xml
        integer :: unit_txt
        ! Input related variables used in parsing loop
        character(len=300) :: comment,buffer, label
        integer :: p
        integer :: ios          ! state variable (0 = ok)
        integer :: line
        logical :: dir_exists
        character :: delimiter
        integer :: i, j, k
        integer :: actcroplen
        integer,parameter :: nanvalue=-9999       ! general NaN value
        integer,dimension(:),allocatable :: dummy
        character(len=300) :: dir_ic, dir_fc

        ErrorFlag = 0
        line = 0
        ios = 0
        actcroplen = 0
        dir_ic = xml%sim%initial_condition
        dir_fc = xml%sim%final_condition
        allocate(xml%sim%clock(3))
        
        ! get a free number to open the file and open it
        call seek_un(errorflag,unit_txt)!
        open(unit_txt,file=trim(file_xml), status="old", action="read", IOSTAT = ios)
        
        if (ios /= 0 ) then
            print *, 'Cannot open file ', trim(file_xml), '. The specified file does not exist.'
            stop
        end if
        
        do while (ios == 0)
            read(unit_txt, '(A)', iostat=ios) buffer
            if (ios == 0) then
                line = line + 1
                ! remove white spaces before and after string sequences
                buffer = trim(buffer)
                ! check if there are comment and get only left side of the row
                p = scan(buffer, '#')
                if (p /= 0) then
                    comment = buffer(p+1:)
                    buffer = buffer(1:p-1)
                end if
                ! if buffer is different from blank string, parse it
                if (buffer /= '') then
                    call lower_case(buffer)
                    p = scan(buffer, '=')
                    label = buffer(1:p-1)
                    buffer = buffer(p+1:)
                    select case (label)
                        case ('outputpath')
                            read(buffer, *, iostat=ios) xml%sim%path
                            ! check if the dir exists or make a new dir
                            inquire(file=trim(xml%sim%path), exist=dir_exists)   ! dir_exists will be TRUE if the directory exists
                            if (dir_exists .eqv. .true.) then
                                print *,'The directory ', trim(xml%sim%path), ' already exists and will be updated'
                                print *, " <enter> to continue "
                                read *
                            else
                                ! TODO: intrinsic 'system' not included in std2008
                                call get_environment_variable('DELIMITER',delimiter)
                                call system('mkdir '//delimiter//trim(xml%sim%path))
                            end if
                        case ('inputpath') ! path to spazialized input files !
                            read(buffer, *, iostat=ios) xml%sim%input_path
                        case ('meteopath') ! path to meteo files !
                            read(buffer, *, iostat=ios) xml%sim%meteo_path
                        case ('meteofilename') ! name of the file with the list of weather station !
                            read(buffer, *, iostat=ios) xml%sim%ws_list_fn
                        case ('phenopath') ! path to the phenophase files !
                            read(buffer, *, iostat=ios) xml%sim%pheno_path
                        case ('phenofileroot') ! the sub string to use as root of phenophase subfolder !
                            read(buffer, *, iostat=ios) xml%sim%pheno_root
                        case ('irrmethpath') ! path to meteo files !
                            read(buffer, *, iostat=ios) xml%sim%irr_met_path
                        case ('irrmethfilename') ! name of the file with the list of weather station !
                            read(buffer, *, iostat=ios) xml%sim%irr_met_list_fn
                        case ('watsourpath') ! path to water sources subfolder
                            read(buffer, *, iostat=ios) xml%sim%watsour_path
                        case ('mode') ! running mode: 0 - no irrigation; 1 - consumption; 2 - water needs to field capacity; 3 - water needs with fixed volume
                            read(buffer, *, iostat=ios) xml%sim%mode
                            ! if no irrigation set startIrrSeason and endIrrSeason to -1
                            if (xml%sim%mode == 0) then
                                xml%sim%start_irr_season = -1
                                xml%sim%end_irr_season = -1
                            end if
                        case ('initialthetaflag') ! T: read init moisture from external file
                            read(buffer, *, iostat=ios) xml%sim%f_init_wc
                        case ('finalthetaflag') ! T: if final soil moisture condition is written on files, F if not
                            read(buffer, *, iostat=ios) xml%sim%f_theta_out
                        case ('startsimulation')    ! date of start simulation
                            read (buffer, *, iostat=ios)
                            call split_date(buffer, xml%sim%start_simulation)
                            xml%sim%start_simulation%doy = calc_doy(xml%sim%start_simulation%day, xml%sim%start_simulation%month, &
                                & xml%sim%start_simulation%year)
                        case ('endsimulation')    ! date of end simulation
                            read (buffer, *, iostat=ios)
                            call split_date(buffer, xml%sim%end_simulation)
                            xml%sim%end_simulation%doy = calc_doy(xml%sim%end_simulation%day, xml%sim%end_simulation%month, &
                                & xml%sim%end_simulation%year)
                        case ('initialconditionpath') ! path to initial condition input files !
                            read(buffer, *, iostat=ios) xml%sim%initial_condition
                            inquire(FILE=trim(xml%sim%initial_condition), EXIST=dir_exists)   ! dir_exists will be TRUE if the directory exists
                            if (dir_exists .eqv. .false.) then
                                print *,'The directory ', trim(xml%sim%initial_condition),&
                                ' does not exist. No initial condition can be uploaded.'
                                print *, 'Execution will be aborted...'
                                stop
                            end if
                        case ('initialcondition') ! initial condition input filenames (root)
                            read(buffer, *, iostat=ios) xml%sim%thetaI_0_fn
                            xml%sim%thetaI_0_fn = trim(xml%sim%thetaI_0_fn)//'I'
                            read(buffer, *, iostat=ios) xml%sim%thetaII_0_fn
                            xml%sim%thetaII_0_fn = trim(xml%sim%thetaII_0_fn)//'II'
                        case ('finalconditionpath') ! path to final condition input files !
                            read(buffer, *, iostat=ios) xml%sim%final_condition
                            inquire(FILE=trim(xml%sim%final_condition), EXIST=dir_exists)   ! dir_exists will be TRUE if the directory exists
                            if (dir_exists .eqv. .true.) then
                                print *,'The directory ', trim(xml%sim%final_condition), ' already exists and will be updated'
                            else
                                ! TODO: intrinsic 'system' non inclusa nello standard std2008
                                call get_environment_variable('DELIMITER',delimiter)
                                call system('mkdir '//delimiter//trim(xml%sim%final_condition))
                            end if
                        case ('finalcondition') ! final condition input filenames (root)
                            read(buffer, *, iostat=ios) xml%sim%thetaI_end_fn
                            xml%sim%thetaI_end_fn = trim(xml%sim%thetaI_end_fn)//'I'
                            read(buffer, *, iostat=ios) xml%sim%thetaII_end_fn
                            xml%sim%thetaII_end_fn = trim(xml%sim%thetaII_end_fn)//'II'
                        case ('capillaryflag') ! T: calculate capillary rise
                            read(buffer, *, iostat=ios) xml%sim%f_cap_rise
                        case('soilusevarflag') ! T: land use update yearly
                            read(buffer, *, iostat=ios) xml%sim%f_soiluse
                        case ('meteostatweightnum') ! number of weather stations
                            read(buffer, *, iostat=ios)xml%sim%n_ws
                        case ('meteostattotnum') ! number of areas of weather stations TODO: check
                            read(buffer, *, iostat=ios)xml%sim%n_voronoi
                        case ('monthlyflag') ! monthly outputs flag !
                            select case (trim(adjustl(buffer)))
                                case ('monthly', 'month', 'm', 't')
                                    xml%sim%step_out = 0
                                case ('weekly', 'week', 'w')
                                    xml%sim%step_out = 1
                                case ('periodic', 'f')
                                    xml%sim%step_out = 2
                                case default
                                    xml%sim%step_out = 0
                                    print *, "Monthly output is printed"
                                    read *
                            end select
                        case ('weekday')
                            select case (trim(adjustl(buffer)))
                                case ('monday', 'mon', '1')
                                    xml%sim%weekday = 1
                                case ('tuesday', 'tue', '2')
                                    xml%sim%weekday = 2
                                case ('wednesday', 'wed', '3')
                                    xml%sim%weekday = 3
                                case ('thursday', 'thu', '4')
                                    xml%sim%weekday = 4
                                case ('friday', 'fri', '5')
                                    xml%sim%weekday = 5
                                case ('saturday', 'sat', '6')
                                    xml%sim%weekday = 6
                                case ('sunday', 'sun', '7')
                                    xml%sim%weekday = 7                                                                       
                                case default
                                    xml%sim%weekday = 1
                                    print *, "Weekly output is printed each Monday"
                                    read *
                            end select
                        case ('startdate')
                            read(buffer, *, iostat=ios)xml%sim%clock(1)
                        case ('enddate')
                            read(buffer, *, iostat=ios)xml%sim%clock(2)
                        case ('deltadate')
                            read(buffer, *, iostat=ios)xml%sim%clock(3)
                        case ('soilusesnum') ! number of land uses to be processed
                            read(buffer, *, iostat=ios)xml%sim%n_lus
                        !!!! %AB% TODO: to be checked !!!!!!!!!!!!!!!!
                        case ('simulatedsoiluses') ! list of simulated soil uses
                            read (buffer, *, iostat=ios)
                            allocate(dummy(xml%sim%n_lus))
                            dummy = nanvalue
                            dummy = string_to_integers(buffer(2:len_trim(buffer)), " ")
                        ! case ('randsowdaysseed')
                        !           read(buffer, *, iostat=ios) xml%sim%rand_seed
                        case ('randsowdayssym')
                            select case (trim(adjustl(buffer)))
                                case ('symmetric', 's', '1', 't', 'true')
                                    xml%sim%rand_symmetry = .true.
                                case ('asymmetric', 'as', '0', 'f', 'false')
                                    xml%sim%rand_symmetry = .false.
                                case default
                                    xml%sim%rand_symmetry = .true.
                                    print *, "Emergence date is distributed symmetrically"
                                    read *
                            end select
                        case ('randsowdayswind') ! range of sowinf!
                            read(buffer, *, iostat=ios) xml%sim%sowing_range
                        ! parameter for percolation booster
                        case ('01q_eva')  ! 10° percentile of the first layer
                            read(buffer, *, iostat=ios)xml%sim%quantiles(1,1)
                        case ('09q_eva')  ! 90° percentile of the first layer
                            read(buffer, *, iostat=ios)xml%sim%quantiles(1,2)
                        case ('01q_trasp')  ! 10° percentile of the second layer
                            read(buffer, *, iostat=ios)xml%sim%quantiles(2,1)
                        case ('09q_trasp')  ! 90° percentile of the second layer
                            read(buffer, *, iostat=ios)xml%sim%quantiles(2,2)
                        case ('startirrseason') ! start doy of irrigation season
                            read(buffer, *, iostat=ios)xml%sim%start_irr_season
                            if (xml%sim%mode == 0) xml%sim%start_irr_season = -1
                        case ('endirrseason') ! end doy of irrigation season
                            read(buffer, *, iostat=ios)xml%sim%end_irr_season
                            if (xml%sim%mode == 0) xml%sim%end_irr_season = -1
                        CASE ('lim_prec') ! Minimum meaningfull precipitation
                            READ(buffer, *, IOSTAT=ios)xml%sim%h_prec_lim
                        case ('zevap') ! Ze_fix: thickness of the first layer
                            read(buffer, *, iostat=ios)xml%depth%ze_fix
                        case ('zroot') ! Zr_fix: default thickness of the second layer
                            read(buffer, *, iostat=ios)xml%depth%zr_fix
                        case ('lambdacn') ! lambda_CN: coefficent for calculate the initial abstraction
                            read(buffer, *, iostat=ios)xml%sim%lambda_cn                       
                        case ('dtxmode') ! DTx mode calculation
                            read(buffer, *, iostat=ios)xml_dtx%mode
                        case ('dtxnumxs') ! Number of indicators to calculate
                            read(buffer, *, iostat=ios)xml_dtx%temp%n_ind
                        case ('dtx_x') ! number of days to use tocalculate the dtx indicators
                            allocate(xml_dtx%temp%x(xml_dtx%temp%n_ind))
                            read(buffer, *, iostat=ios) xml_dtx%temp%x
                        case ('dtxdeltadate') ! validity of the report
                            read(buffer, *, iostat=ios)xml_dtx%temp%td
                        case ('dtxdelaydays') ! Shift di report delivering
                            read(buffer, *, iostat=ios)xml_dtx%temp%delay
                        case ('dtxmincard') ! Number of values to obtain meaningfull statistics
                            read(buffer, *, iostat=ios)xml_dtx%n
                        case ('forecast_day') ! Number days to consider for cumulative precipitation
                            read(buffer, *, iostat=ios)xml_dtx%forecast_day
                        case default ! all other cases ... !
                            print *, 'Skipping invalid or obsolete label <',trim(label),'> at line', line, &
                                & ' of file: ', trim(file_xml)
                    end select
                end if
            end if
        end do

        if (xml%sim%step_out == 1) then
            ! take into account the last week of the previos year and the first of the following year
            allocate(xml%sim%intervals(54))
        else if (xml%sim%step_out == 2) then
            allocate(xml%sim%intervals(int(((xml%sim%clock(2)-xml%sim%clock(1))/xml%sim%clock(3)))))
            xml%sim%intervals=xml%sim%clock(3)       
        end if

        if (xml%sim%initial_condition == dir_ic .and. xml%sim%input_path /= dir_ic) then
            xml%sim%initial_condition = xml%sim%input_path
        end if
        if (xml%sim%final_condition == dir_fc .and. xml%sim%path /= dir_fc) then
            xml%sim%final_condition = xml%sim%path
        end if
    
        close(unit_txt)!
        !
        if (allocated(dummy)) then
            actcroplen = size(dummy)
        else
            actcroplen = xml%sim%n_lus
            allocate(dummy(xml%sim%n_lus))
            dummy(1:xml%sim%n_lus) = (/(i, i=1, xml%sim%n_lus, 1)/)
        end if
        !
        allocate(xml%sim%lu_list(actcroplen))
        k = 1
        do j = 1, actcroplen
            do i = k, xml%sim%n_lus
                if (dummy(i) /= nanvalue) then
                    xml%sim%lu_list(j) = dummy(i)
                    k = i + 1
                    exit
                end if
            end do
        end do
        deallocate(dummy)
        !
        allocate(dummy(xml%sim%n_lus))
        dummy(1:xml%sim%n_lus) = (/(i, i=1, xml%sim%n_lus, 1)/)
        do i = 1, actcroplen
            where (dummy == xml%sim%lu_list(i)) dummy = nanvalue
        end do
        !
        allocate(xml%sim%no_lu_list(xml%sim%n_lus - actcroplen))
        k = 1
        
        do j = 1, size(xml%sim%no_lu_list)
            do i = k, xml%sim%n_lus
                if (dummy(i) /= nanvalue) then
                    xml%sim%no_lu_list(j) = dummy(i)
                    k = i + 1
                    exit
                end if
            end do
        end do

        if (debug .eqv. .true.) then
            print *, "Simulated crop IDs:", xml%sim%lu_list
            print *, "Not simulated crop IDs:", xml%sim%no_lu_list
        end if
    end subroutine read_sim_parameters

    ! TODO: %EAC% at the moment the id of the irrigation method define its order in the list
    subroutine read_irr_method(irr_method_fn, met, debug)!
        implicit none!
        character(len=*),intent(in) :: irr_method_fn!
        type(par_method),dimension(:),intent(inout) :: met!
        logical, intent(in), optional :: debug
        ! Input related variables used in parsing loop
        character(len=300) :: comment,buffer, label
        integer :: c
        integer :: p, ErrorFlag
        character(len=300)::i
        integer :: unit_txt
        integer :: line
        integer :: ios      ! state variable (0 = ok)
        integer :: num, e
        line = 0
        ios = 0
        e = 1
        call seek_un(ErrorFlag,unit_txt)!
        open(unit_txt,file=trim(irr_method_fn),action="read", IOSTAT = ios)!
        if (ios /= 0 ) then
            print *, "Cannot open file ", trim(irr_method_fn), ". The specified file does not exist."
            print *, 'Execution will be aborted...'
            stop
        end if
        do while (ios == 0)
            read (unit_txt, '(A)', iostat = ios) buffer
            if (ios == 0) then
                line = line +1
                buffer = trim(buffer)
                c = scan(buffer, '#')
                if (c /= 0) then
                    comment = buffer(c+1:)
                    buffer = buffer(1:c-1)
                end if
                if (buffer /= '') then
                    call lower_case(buffer)
                    c = scan(buffer, '=')
                    label = buffer(1:c-1)
                    buffer = buffer(c+1:)
                    
                    select case (label)
                        case ('id')
                            read (buffer, *, iostat = ios) p
                            if (p > size(met, 1)) then
                                print *, 'The irrigation method Id, ', p, ', in the file "', trim(irr_method_fn), &
                                    & '", is above upper bound of the number of irrigation methods: ', size (met, 1)
                                print *, 'Execution will be aborted...'
                                stop
                            end if
                            allocate(met(p)%freq(24))
                            read (buffer, *, iostat = ios) met(p)%id  ! id of the irrigation method
                            met(p)%f_wet = cost_fwEva ! set to the default value
                        case ('qwat','qadaq')
                            read (buffer, *, iostat = ios) met(p)%h_irr ! irrigation depth [mm] !
                        case ('fw')
                            read (buffer, *, iostat = ios) met(p)%f_wet    ! soil wetted fraction [-][0-1]
                        case ('k_stress','irr_th_mn')
                            read (buffer, *, iostat = ios) met(p)%irr_th_ms ! fraction of RAW that activate irrigation for monitored sources
                        case ('k_stresswells', 'irr_th_unm')
                            read (buffer, *, iostat = ios) met(p)%irr_th_unm ! fraction of RAW that activate irrigation for unmonitored sources
                        case ('min_a')
                            read (buffer, *, iostat = ios) met(p)%a_min ! min value of the "a" parameter for percolation model
                        case ('max_a')
                            read (buffer, *, iostat = ios)  met(p)%a_max ! max value of the "a" parameter for percolation model
                        case ('min_b')
                            read (buffer, *, iostat = ios) met(p)%b_min ! min value of the "b" parameter for percolation model
                        case ('max_b')
                            read (buffer, *, iostat = ios)  met(p)%b_max ! max value of the "b" parameter for percolation model
                        case ('interceptionflag') ! Flag interception
                            read (buffer, *, IOSTAT = ios) i
                            If ((i == 'true') .or. (i == 't'))  then
                                met(p)%f_interception = 1
                            else
                                met(p)%f_interception = 0
                            end if
                        ! irrigation losses model parameters
                        case ('a')
                            read (buffer, *, iostat = ios)  met(p)%a_loss ! a!
                        case ('b')
                            read (buffer, *, iostat = ios)  met(p)%b_loss ! b!
                        case ('c')
                            read (buffer, *, iostat = ios)  met(p)%c_loss ! c!
                        case ('irr_starts')
                            read (buffer, *, iostat = ios)  met(p)%irr_starts ! doy when irrigation season starts
                        case ('irr_ends')
                            read (buffer, *, iostat = ios)  met(p)%irr_ends ! doy when irrigation season starts
                        case default
                            read (label, *, iostat = e) num
                            if ((e==0) .and. (num>=1) .and. (num<=24)) then
                                read (buffer, *, iostat = ios) met(p)%freq(num)
                            else
                                print *, 'Skipping invalid label <',trim(label),'> at line', line, ' of file: ', trim(irr_method_fn)
                                print *, 'Execution will be aborted...'
                                stop
                            end if
                    end select
                end if
            end if
        end do
        close(unit_txt)!

        if (debug .eqv. .true.) then
            print *,'====== DEBUG: irrigation method: ', irr_method_fn, ' ====='
            print *,'Id = ',  p
            print *,'Qwat = ',  met(p)%h_irr
            print *,'fw = ', met(p)%f_wet
            print *,'Irrigation threshold monitored sources = ',  met(p)%irr_th_ms
            Print *,'Irrigation threshold unmonitored sources = ',  met(p)%irr_th_unm
            print *,'min_a = ', met(p)%a_min
            print *,'max_a = ',  met(p)%a_max
            print *,'min_b = ', met(p)%b_min
            print *,'max_b = ',  met(p)%b_max
            If (met(p)%f_interception == 1)  then
                print *, 'InterceptionFlag = TRUE'
            else
                print *, 'InterceptionFlag = FALSE'
            end if
            print *,'a = ',  met(p)%a_loss
            print *,'b = ', met(p)%b_loss
            print *,'c = ', met(p)%c_loss
            print *,'freqirr = ', met(p)%freq
            print *,'===== END DEBUG ====='
        end if
    end subroutine read_irr_method!

    subroutine read_all_irr_methods(xml, ErrorFlag, debug)
        implicit none
        logical, intent(in) :: debug
        integer, intent(out) :: ErrorFlag
        type(parameters), intent(INout) :: xml
        integer :: free_unit
        integer :: i
        character(len=200) :: irr_method_fn
        ! Input related variables used in parsing loop
        character(len=300) :: comment,buffer, label
        integer :: p
        integer :: ios          ! state variable (0 = ok)
        integer :: line, tablestart

        ErrorFlag = 0
        line = 0; tablestart = 0
        ios = 0
        call seek_un(errorflag,free_unit)!
        print*,'In read_irr_parameters, irr_met_path:',xml%sim%irr_met_path
        print*,'In read_irr_parameters, irr_met_filename:',xml%sim%irr_met_list_fn
        
        open(free_unit,file=trim(xml%sim%irr_met_path)//trim(xml%sim%irr_met_list_fn),action="read")!
        line = 0; tablestart = 0
        ios = 0
        do while (ios == 0)
            read(free_unit, '(A)', iostat=ios) buffer
            if (ios == 0) then
                line = line + 1
                buffer = trim(buffer)
                p = scan(buffer, '#')
                if (p /= 0) then
                    comment = buffer(p+1:)
                    buffer = buffer(1:p-1)
                end if
                if (buffer /= '') then
                    call lower_case(buffer)
                    p = scan(buffer, '=')
                    label = buffer(1:p-1)
                    buffer = buffer(p+1:)
                    select case (label)
                        case ('irrmethnum')
                             ! number of irrigation method considered in the simulation
                            read(buffer, *, iostat=ios) xml%sim%n_irr_meth
                            allocate (xml%irr%met(xml%sim%n_irr_meth))!
                            ! set defualt values for irrigation methods from general setup
                            do i = 1, xml%sim%n_irr_meth
                                xml%irr%met(i)%irr_starts = xml%sim%start_irr_season
                                xml%irr%met(i)%irr_ends = xml%sim%end_irr_season
                            end do
                        case ('list')
                            tablestart = line
                            ! List of the file with the irrigation methods
                            do i = 1, xml%sim%n_irr_meth
                                read(free_unit, *, iostat=ios) irr_method_fn  
                                line = line + 1
                                call read_irr_method(trim(xml%sim%irr_met_path)//irr_method_fn, xml%irr%met, debug)!
                            end do
                        case ('endlist')
                            if ((line - tablestart - 1) > xml%sim%n_irr_meth) then
                                stop 'Irrigation methods number (IrrMethNum) is higher than irrigation methods listed'
                            else if ((line - tablestart - 1) < xml%sim%n_irr_meth) then
                                stop 'Irrigation methods number (IrrMethNum) is lower than irrigation methods listed'
                            end if
                        case default ! all other cases ... !
                            print *, 'Skipping invalid label <',trim(label),'> at line ', line, ' of file: ', &
                                & trim(xml%sim%irr_met_list_fn)
                            print *, 'Execution will be aborted...'
                            stop
                    end select
                end if
            end if
        end do
        close (free_unit)
        !
    end subroutine read_all_irr_methods

    
    subroutine read_spatial_info(info_spat, extent, sim, tab_CN2, tab_CN3, theta2_rice, met)!
        ! read all spatially distributed parameters and check parameters
        ! init tab_CN TODO: separate the function
        
        type(bound),intent(in)::extent!
        type(simulation),intent(inout)::sim!
        type(par_method),dimension(:),intent(inout)::met!
        type(spatial_info),intent(out)::info_spat!
        real(dp),dimension(:,:,:),intent(out):: tab_CN2, tab_CN3
        integer::ios!
        type(soil2_rice)::theta2_rice       !
        integer :: line
        !!
        line = 0
        ios = 0

        call read_grid_files(info_spat, extent, sim)
        call check_grid(info_spat, sim)
        !
        select case (sim%mode)
            case (0)
                print *, 'Irrigation distribution parameters are not set'
            case (1:4)
                call read_irr_grid(info_spat, extent, sim, met)
                call check_irr_grid(info_spat, sim)
            case default
                stop 'Simulation mode has been incorrectly set'
        end select
        !
        call read_rice_parameters(sim, theta2_rice)
        call init_cn_table(tab_CN2, tab_CN3)
        !
    end subroutine read_spatial_info!
    !
    subroutine read_grid_files(info_spat, extent, sim)!
        !Read *.asc files content into info_spat variable
        implicit none!
        type(bound),intent(in)::extent!
        type(simulation),intent(inout)::sim!
        type(spatial_info),intent(out)::info_spat!
        integer::k!
        character(len=30)::k_str
        character(len=300)::dir, dir_ic
        character(len=30)::start_year
        !!
        dir = sim%input_path
        dir_ic = sim%initial_condition

        call read_grid(trim(dir)//trim(sim%domain_fn)//'.asc',info_spat%domain,sim,extent)
        ! Water content parameters
        call read_grid(trim(dir)//trim(sim%thetaI_SAT_fn)//'.asc',info_spat%theta(1)%sat,sim,extent);
        call read_grid(trim(dir)//trim(sim%thetaII_SAT_fn)//'.asc',info_spat%theta(2)%sat,sim,extent);
        call read_grid(trim(dir)//trim(sim%thetaI_FC_fn)//'.asc',info_spat%theta(1)%fc,sim,extent)
        call read_grid(trim(dir)//trim(sim%thetaII_FC_fn)//'.asc',info_spat%theta(2)%fc,sim,extent)
        call read_grid(trim(dir)//trim(sim%thetaI_WP_fn)//'.asc',info_spat%theta(1)%wp,sim,extent)
        call read_grid(trim(dir)//trim(sim%thetaII_WP_fn)//'.asc',info_spat%theta(2)%wp,sim,extent)
        call read_grid(trim(dir)//trim(sim%thetaI_r_fn)//'.asc',info_spat%theta(1)%r,sim,extent)
        call read_grid(trim(dir)//trim(sim%thetaII_r_fn)//'.asc',info_spat%theta(2)%r,sim,extent)
        ! Brooks and Corey's equation parameters
        call read_grid(trim(dir)//trim(sim%ksat_I_fn)//'.asc',info_spat%k_sat(1),sim,extent)
        call read_grid(trim(dir)//trim(sim%ksat_II_fn)//'.asc',info_spat%k_sat(2),sim,extent)
        call read_grid(trim(dir)//trim(sim%n_I_fn)//'.asc',info_spat%fact_n(1),sim,extent)
        call read_grid(trim(dir)//trim(sim%n_II_fn)//'.asc',info_spat%fact_n(2),sim,extent)
        ! CN method parameters: slope, hydrologic condition (dren) and hydrologic soil group (hydr_gr)
        call read_grid(trim(dir)//trim(sim%slope_fn)//'.asc',info_spat%slope,sim,extent)
        info_spat%slope%mat = info_spat%slope%mat/100  ! Slopes are converted to percentage: x/100 = 0.x
        call read_grid(trim(dir)//trim(sim%dren_fn)//'.asc',info_spat%drainage,sim,extent)
        call read_grid(trim(dir)//trim(sim%hydr_group_fn)//'.asc',info_spat%hydr_gr,sim,extent)
        !
        ! read initial moisture condition
        if(sim%f_init_wc .eqv. .true.)then  ! from external files
            call read_grid(trim(dir_ic)//trim(sim%thetaI_0_fn)//'.asc',info_spat%theta(1)%old,sim,extent)
            call read_grid(trim(dir_ic)//trim(sim%thetaII_0_fn)//'.asc',info_spat%theta(2)%old,sim,extent)
        else                            ! from the first year of simulation warm-up
            call read_grid(trim(dir)//trim(sim%thetaI_FC_fn)//'.asc',info_spat%theta(1)%old,sim,extent)
            call read_grid(trim(dir)//trim(sim%thetaII_FC_fn)//'.asc',info_spat%theta(2)%old,sim,extent)
        end if!
        !
        ! read initial land use map
        if (sim%f_soiluse .eqv. .false.)then  ! Static land use map
            call read_grid(trim(dir)//trim(sim%soiluse_fn)//'.asc',info_spat%soil_use_id,sim,extent)
            call overlay_domain (info_spat%soil_use_id,info_spat%domain)   ! Set to NaN non simulated soil uses
        else                                ! Dynamic land use map
            write (start_year, *) sim%start_year
            call read_grid(trim(dir)//trim(sim%soiluse_fn)//'_'//trim(adjustl(start_year))//'.asc', &
                & info_spat%soil_use_id,sim,extent)
            print *, "First year soil use has been read"
        end if
        !
        ! read parameters for capillary rise and water table depth
        if(sim%f_cap_rise .eqv. .true.)then!
            call read_grid(trim(dir)//trim(sim%wat_table_fn)//'.asc',info_spat%wat_tab,sim,extent)
            call read_grid(trim(dir)//trim(sim%ParRisCap_a3_fn)//'.asc',info_spat%a3,sim,extent)
            call read_grid(trim(dir)//trim(sim%ParRisCap_a4_fn)//'.asc',info_spat%a4,sim,extent)
            call read_grid(trim(dir)//trim(sim%ParRisCap_b1_fn)//'.asc',info_spat%b1,sim,extent)
            call read_grid(trim(dir)//trim(sim%ParRisCap_b2_fn)//'.asc',info_spat%b2,sim,extent)
            call read_grid(trim(dir)//trim(sim%ParRisCap_b3_fn)//'.asc',info_spat%b3,sim,extent)
            call read_grid(trim(dir)//trim(sim%ParRisCap_b4_fn)//'.asc',info_spat%b4,sim,extent)
        else!
            info_spat%wat_tab = info_spat%domain!
            info_spat%wat_tab%mat=100!
            info_spat%a3 = info_spat%domain
            info_spat%a4 = info_spat%domain
            info_spat%b1 = info_spat%domain
            info_spat%b2 = info_spat%domain
            info_spat%b3 = info_spat%domain
            info_spat%b4 = info_spat%domain
        end if!
        !
        ! read weather station weights
        allocate(info_spat%weight_ws(sim%n_ws))!
        do k=1,size(info_spat%weight_ws)!
            write(k_str,*) k!
            call read_grid(trim(dir)//trim(sim%meteoweight_fn)//'_'//trim(adjustl(k_str))//'.asc',&
                & info_spat%weight_ws(k),sim, extent)!
        end do!
        
        ! read shape area, if the file exists
        inquire (file=trim(dir)//trim(sim%shapearea_fn)//'.asc', exist=sim%f_shapearea)
        if(sim%f_shapearea .eqv. .true.) then
            call read_grid(trim(dir)//trim(sim%shapearea_fn)//'.asc', info_spat%cell_area,sim,extent)
        end if
        
    end subroutine read_grid_files
    !
    subroutine read_irr_grid(info_spat, extent, sim, met)!
        ! distributed parameters for irrigation
        type(bound),intent(in)::extent!
        type(simulation),intent(inout)::sim!
        type(par_method),dimension(:),intent(inout)::met!
        type(spatial_info),intent(out)::info_spat!
        character(len=300)::dir
        character(len=30)::start_year
        !!
        dir = sim%input_path

        select case (sim%mode)
            case (1)
                ! read parameter for irrigation practice
                call read_grid(trim(dir)//trim(sim%irr_units_fn)//'.asc',info_spat%irr_unit_id,sim,extent)
                if (sim%f_soiluse .eqv. .false.) then
                    call read_grid(trim(dir)//trim(sim%id_irr_meth_fn)//'.asc',info_spat%irr_meth_id,sim,extent)
                else
                    write(start_year,*) sim%start_year
                    call read_grid (trim(dir)//trim(sim%id_irr_meth_fn)//'_'//trim(adjustl(start_year))//'.asc',&
                        & info_spat%irr_meth_id,sim,extent)
                end if
                call read_grid(trim(dir)//trim(sim%eff_net_fn)//'.asc',info_spat%eff_net,sim,extent)
            case (2)
                if (sim%f_soiluse .eqv. .false.) then
                    call read_grid(trim(dir)//trim(sim%id_irr_meth_fn)//'.asc',info_spat%irr_meth_id,sim,extent)
                    call read_grid(trim(dir)//trim(sim%eff_irr_fn)//'.asc',info_spat%eff_met,sim,extent)
                else
                    write(start_year,*) sim%start_year
                    call read_grid (trim(dir)//trim(sim%id_irr_meth_fn)//'_'//trim(adjustl(start_year))//'.asc',&
                            & info_spat%irr_meth_id,sim,extent)
                    call read_grid (trim(dir)//trim(sim%eff_irr_fn)//'_'//trim(adjustl(start_year))//'.asc',&
                            & info_spat%eff_met,sim,extent)
                end if
                !call read_matrices(trim(dir)//trim(sim%eff_rete_fn)//'.asc',info_spat%eff_rete,sim,confini) - RR
            case (3)
                if (sim%f_soiluse .eqv. .false.) then
                    call read_grid(trim(dir)//trim(sim%id_irr_meth_fn)//'.asc',info_spat%irr_meth_id,sim,extent)
                else
                    write(start_year,*) sim%start_year
                    call read_grid (trim(dir)//trim(sim%id_irr_meth_fn)//'_'//trim(adjustl(start_year))//'.asc',&
                        & info_spat%irr_meth_id,sim,extent)
                end if
                !call read_matrices(trim(dir)//trim(sim%eff_rete_fn)//'.asc',info_spat%eff_rete,sim,confini) - RR
            case (4)
                call read_grid(trim(dir)//trim(sim%irr_units_fn)//'.asc',info_spat%irr_unit_id,sim,extent)
                call overlay_domain(info_spat%irr_unit_id,info_spat%domain)
                if (sim%f_soiluse .eqv. .false.) then
                    call read_grid(trim(dir)//trim(sim%id_irr_meth_fn)//'.asc',info_spat%irr_meth_id,sim,extent)
                    call read_grid(trim(dir)//trim(sim%eff_irr_fn)//'.asc',info_spat%eff_met,sim,extent)
                else
                    write(start_year,*) sim%start_year
                    call read_grid (trim(dir)//trim(sim%id_irr_meth_fn)//'_'//trim(adjustl(start_year))//'.asc',&
                        & info_spat%irr_meth_id,sim,extent)
                    call read_grid (trim(dir)//trim(sim%eff_irr_fn)//'_'//trim(adjustl(start_year))//'.asc',&
                        & info_spat%eff_met,sim,extent)
                end if
                call read_grid(trim(dir)//trim(sim%eff_net_fn)//'.asc',info_spat%eff_net,sim,extent)
            case default
        end select
        print *, "Irrigation distribution parameters have been read"
        ! calculate the parameter for the percolative booster TODO: move away from the reading routines
        call calc_perc_booster_pars(info_spat,met,sim%quantiles)!
        !
    end subroutine read_irr_grid!
    !
    subroutine check_grid(info_spat, sim)!
        ! check grid extension respect to the domain
        ! Note: only water table depth is not verified
        ! so, it cannot be called inside the reading routine
        ! TODO: reorganize the function
        type(simulation),intent(inout)::sim!
        type(spatial_info),intent(out)::info_spat!
        integer::k!
        character(len=30)::k_str!
        !!
        call overlay_domain(info_spat%theta(1)%sat,info_spat%domain)
        call overlay_domain(info_spat%theta(2)%sat,info_spat%domain)
        call overlay_domain(info_spat%theta(1)%fc,info_spat%domain)
        call overlay_domain(info_spat%theta(2)%fc,info_spat%domain)
        call overlay_domain(info_spat%theta(1)%wp,info_spat%domain)
        call overlay_domain(info_spat%theta(2)%wp,info_spat%domain)
        call overlay_domain(info_spat%theta(1)%r,info_spat%domain)
        call overlay_domain(info_spat%theta(2)%r,info_spat%domain)
        call overlay_domain(info_spat%k_sat(1),info_spat%domain)
        call overlay_domain(info_spat%k_sat(2),info_spat%domain)
        call overlay_domain(info_spat%fact_n(1),info_spat%domain)
        call overlay_domain(info_spat%fact_n(2),info_spat%domain)
        call overlay_domain(info_spat%slope,info_spat%domain)
        call overlay_domain(info_spat%drainage,info_spat%domain)
        call overlay_domain(info_spat%hydr_gr,info_spat%domain)
        !
        if(sim%f_init_wc .eqv. .true.)then  ! Initial moisture condition defined externally
            where(info_spat%domain%mat/=info_spat%domain%header%nan .and. &
                    & info_spat%theta(1)%old%mat==info_spat%theta(1)%old%header%nan)!
                info_spat%theta(1)%old%mat = info_spat%theta(1)%fc%mat ! use theta_fc to fill the gaps 
            end where
            where(info_spat%domain%mat/=info_spat%domain%header%nan .and. &
                    & info_spat%theta(2)%old%mat==info_spat%theta(2)%old%header%nan)!
                info_spat%theta(2)%old%mat = info_spat%theta(2)%fc%mat!
            end where!
        else                                ! Initial condition calculated from the first warm-up period 
            call overlay_domain(info_spat%theta(1)%old,info_spat%domain)
            call overlay_domain(info_spat%theta(2)%old,info_spat%domain)
        end if!
        !
        ! check the parameters of the capillary rise model
        if(sim%f_cap_rise .eqv. .true.)then!
            call overlay_domain(info_spat%wat_tab,info_spat%domain)
            if(minval(info_spat%wat_tab%mat,info_spat%wat_tab%mat/=info_spat%wat_tab%header%nan)<0.)then!
                print*,"A null value has been set where groundwater levels are negative"!
                where(info_spat%wat_tab%mat<0. .and. info_spat%wat_tab%mat/=info_spat%wat_tab%header%nan)info_spat%wat_tab%mat=0.!
            end if!
            call overlay_domain(info_spat%a3,info_spat%domain)
            call overlay_domain(info_spat%a4,info_spat%domain)
            call overlay_domain(info_spat%b1,info_spat%domain)
            call overlay_domain(info_spat%b2,info_spat%domain)
            call overlay_domain(info_spat%b3,info_spat%domain)
            call overlay_domain(info_spat%b4,info_spat%domain)
            ! Set the minimum value of water table depth
            ! TODO: check what happen if the depth is not set
            ! TODO: check the minimum value of water table
            ! TODO: check when root goes into water table or cut roots!
            ! %EAC%: moved in the watre balance module
            !where(info_spat%wattab%mat<1.d0 .and. info_spat%wattab%mat/=info_spat%wattab%intest%nan)info_spat%wattab%mat=1.d0!
        end if!
        !
        ! weather station weights maps 
        do k=1,size(info_spat%weight_ws)!
            write(k_str,*) k!
            call overlay_domain(info_spat%weight_ws(k),info_spat%domain)
        end do!
        !
        info_spat%backup_domain = info_spat%domain
        info_spat%backup_domain%mat = info_spat%domain%mat
        !      
        ! land use
        call overlay_domain(info_spat%soil_use_id,info_spat%domain)
        if(minval(info_spat%soil_use_id%mat,info_spat%soil_use_id%mat/=info_spat%soil_use_id%header%nan) < 1 &
                .or. maxval(info_spat%soil_use_id%mat) > sim%n_lus)then
            print *,"Soil use maps have soil uses not defined in crop database"
            print *,"Verify the maximum allowed crop uses (SoilUsesNum) and soil maps"
            print *, 'Execution will be aborted...'
            stop!
        end if
    end subroutine check_grid!
    !
    subroutine check_irr_grid(info_spat, sim)!
        ! check irrigation parameters grids
        type(simulation),intent(inout)::sim!
        type(spatial_info),intent(out)::info_spat!
        integer::ios!
        character(len=300)::dir
        integer :: line
        !!
        dir = sim%input_path
        line = 0
        ios = 0
        select case (sim%mode)
            case (1)
                call overlay_domain(info_spat%irr_unit_id,info_spat%domain)
                ! fix not defined irrigation methods
                call set_default_par(info_spat%irr_meth_id, info_spat%domain, 1)
                call set_default_par(info_spat%eff_net,info_spat%domain, 1.0D0)
            case (2)
                call set_default_par(info_spat%irr_meth_id, info_spat%domain, 1)
                !call check_mat_irrigation(info_spat%eff_rete,info_spat%domain, 1.0D0) - %RR%
                call set_default_par(info_spat%eff_met, info_spat%domain, 1.0D0)
            case (3)
                call set_default_par(info_spat%irr_meth_id, info_spat%domain, 1)
                !call check_mat_irrigation(info_spat%eff_rete,info_spat%domain, 1.0D0) - %RR%
            case (4)
                call overlay_domain(info_spat%irr_unit_id,info_spat%domain)
                call set_default_par(info_spat%irr_meth_id, info_spat%domain, 1)
                call set_default_par(info_spat%eff_met, info_spat%domain, 1.0D0)
                call set_default_par(info_spat%eff_net,info_spat%domain, 1.0D0)
            case default
        end select
        call set_default_par(info_spat%a_perc(1), info_spat%domain, 1.0D0)
        call set_default_par(info_spat%a_perc(2), info_spat%domain, 1.0D0)
        call set_default_par(info_spat%b_perc(1), info_spat%domain, 1.0D0)
        call set_default_par(info_spat%b_perc(2), info_spat%domain, 1.0D0)
    end subroutine check_irr_grid!

    subroutine read_rice_parameters(sim, theta2_rice)!
        ! read the parameters specific for rice paddy
        implicit none!
        type(simulation),intent(inout)::sim!
        integer::errorflag,ios!
        type(soil2_rice)::theta2_rice       !
        character(len=300)::dir
        character(len=300) :: comment,buffer, label
        integer :: p
        integer :: line
        !!
        dir = sim%input_path
        line = 0
        ios = 0

        call seek_un(errorflag,theta2_rice%unit_soil_rice)!
        open(theta2_rice%unit_soil_rice,file=trim(dir)//trim(sim%soil_prop_x_rice_fn),action='read')!
        do while (ios == 0)
            read (theta2_rice%unit_soil_rice, '(A)', iostat = ios) buffer
            if (ios == 0) then
                line = line +1
                buffer = trim(buffer)
                p = scan(buffer, '#')
                if (p /= 0) then
                    comment = buffer(p+1:)
                    buffer = buffer(1:p-1)
                end if
                if (buffer /= '') then
                    call lower_case(buffer)
                    p = scan(buffer, '=')
                    label = buffer(1:p-1)
                    buffer = buffer(p+1:)
                    select case (label)
                        case ('ksat_ii');   read (buffer, *, iostat = ios) theta2_rice%k_sat_2
                        case ('n_ii');      read (buffer, *, iostat = ios) theta2_rice%n_2
                        case ('tetaii_fc'); read (buffer, *, iostat = ios) theta2_rice%theta2_FC
                        case ('tetaii_r');  read (buffer, *, iostat = ios) theta2_rice%theta2_R
                        case ('tetaii_sat');read (buffer, *, iostat = ios) theta2_rice%theta2_SAT
                        case ('tetaii_wp'); read (buffer, *, iostat = ios) theta2_rice%theta2_WP
                        case default
                            print *, 'Skipping invalid label <',trim(label),'> at line', line, ' of file: ', &
                                & trim(sim%soil_prop_x_rice_fn), '. Execution will be aborted...'
                            stop
                    end select
                end if
            end if
        end do
        close(theta2_rice%unit_soil_rice)!
    end subroutine read_rice_parameters
    !
    subroutine init_cn_table(tab_CN2, tab_CN3)!
        ! init the CN table
        ! TODO: replace with the initialization from external file
        implicit none!
        real(dp),dimension(:,:,:),intent(out):: tab_CN2, tab_CN3
        !
        tab_CN2 = tabCN
        where(tab_CN2 > 0) ! calculate tab_CN3
            tab_CN3=tab_CN2*exp(0.00673*(100.-tab_CN2))!
        end where!
        !!
    end subroutine init_cn_table!

    subroutine write_init_grids(info_spat,mode,path,sim)!
        ! save the grid data for the selected area
        implicit none!
        type(spatial_info),intent(in)::info_spat!
        integer, intent(in)::mode
        character(len=*),intent(in)::path!
        type(simulation),intent(in)::sim
        !!
        integer::errorflag,k!
        character(len=30)::meteo_stringa!
        !!
        call write_grid(trim(path)//'out_'//trim(sim%domain_fn)//'.asc',info_spat%domain,errorflag)!
        call write_grid(trim(path)//'out_'//trim(sim%thetaI_FC_fn)//'.asc',info_spat%theta(1)%fc,errorflag)!
        call write_grid(trim(path)//'out_'//trim(sim%thetaII_FC_fn)//'.asc',info_spat%theta(2)%fc,errorflag)!
        call write_grid(trim(path)//'out_'//trim(sim%thetaI_WP_fn)//'.asc',info_spat%theta(1)%wp,errorflag)!
        call write_grid(trim(path)//'out_'//trim(sim%thetaII_WP_fn)//'.asc',info_spat%theta(2)%wp,errorflag)!
        call write_grid(trim(path)//'out_'//trim(sim%thetaI_r_fn)//'.asc',info_spat%theta(1)%r,errorflag)!
        call write_grid(trim(path)//'out_'//trim(sim%thetaII_r_fn)//'.asc',info_spat%theta(2)%r,errorflag)!
        call write_grid(trim(path)//'out_'//trim(sim%thetaI_SAT_fn)//'.asc',info_spat%theta(1)%sat,errorflag)!
        call write_grid(trim(path)//'out_'//trim(sim%thetaII_SAT_fn)//'.asc',info_spat%theta(2)%sat,errorflag)!
        call write_grid(trim(path)//'out_'//trim(sim%slope_fn)//'.asc',info_spat%slope,errorflag)!
        call write_grid(trim(path)//'out_'//trim(sim%dren_fn)//'.asc',info_spat%drainage,errorflag)!
        call write_grid(trim(path)//'out_'//trim(sim%hydr_group_fn)//'.asc',info_spat%hydr_gr,errorflag)!
        call write_grid(trim(path)//'out_'//trim(sim%ksat_I_fn)//'.asc',info_spat%k_sat(1),errorflag)!
        call write_grid(trim(path)//'out_'//trim(sim%ksat_II_fn)//'.asc',info_spat%k_sat(2),errorflag)!
        call write_grid(trim(path)//'out_'//trim(sim%n_I_fn)//'.asc',info_spat%fact_n(1),errorflag)!
        call write_grid(trim(path)//'out_'//trim(sim%n_II_fn)//'.asc',info_spat%fact_n(2),errorflag)!
        call write_grid(trim(path)//'out_'//trim(sim%soiluse_fn)//'.asc',info_spat%soil_use_id,errorflag)!
        call write_grid(trim(path)//'out_'//trim(sim%wat_table_fn)//'.asc',info_spat%wat_tab,errorflag)!
        call write_grid(trim(path)//'out_'//trim(sim%ParRisCap_a3_fn)//'.asc',info_spat%a3,errorflag)!
        call write_grid(trim(path)//'out_'//trim(sim%ParRisCap_a4_fn)//'.asc',info_spat%a4,errorflag)!
        call write_grid(trim(path)//'out_'//trim(sim%ParRisCap_b1_fn)//'.asc',info_spat%b1,errorflag)!
        call write_grid(trim(path)//'out_'//trim(sim%ParRisCap_b2_fn)//'.asc',info_spat%b2,errorflag)!
        call write_grid(trim(path)//'out_'//trim(sim%ParRisCap_b3_fn)//'.asc',info_spat%b3,errorflag)!
        call write_grid(trim(path)//'out_'//trim(sim%ParRisCap_b4_fn)//'.asc',info_spat%b4,errorflag)!
        !!
        select case (mode)
            case (1)
                call write_grid(trim(path)//'out_'//trim(sim%irr_units_fn)//'.asc',info_spat%irr_unit_id,errorflag)!
                call write_grid(trim(path)//'out_'//trim(sim%eff_net_fn)//'.asc',info_spat%eff_net,errorflag)!
                call write_grid(trim(path)//'out_'//trim(sim%id_irr_meth_fn)//'.asc',info_spat%irr_meth_id,errorflag)!
                call write_grid(trim(path)//"out_am_percI.asc",info_spat%a_perc(1),errorflag)!
                call write_grid(trim(path)//"out_am_percII.asc",info_spat%a_perc(2),errorflag)!
                call write_grid(trim(path)//"out_bm_percI.asc",info_spat%b_perc(1),errorflag)!
                call write_grid(trim(path)//"out_bm_percII.asc",info_spat%b_perc(2),errorflag)!
            case (2)
                call write_grid(trim(path)//'out_'//trim(sim%eff_irr_fn)//'.asc',info_spat%eff_met,errorflag)!
                !call write_matrices(trim(path)//'out_'//trim(sim%eff_rete_fn)//'.asc',info_spat%eff_rete,errorflag)! - RR
                call write_grid(trim(path)//'out_'//trim(sim%id_irr_meth_fn)//'.asc',info_spat%irr_meth_id,errorflag)!
                call write_grid(trim(path)//"out_am_percI.asc",info_spat%a_perc(1),errorflag)!
                call write_grid(trim(path)//"out_am_percII.asc",info_spat%a_perc(2),errorflag)!
                call write_grid(trim(path)//"out_bm_percI.asc",info_spat%b_perc(1),errorflag)!
                call write_grid(trim(path)//"out_bm_percII.asc",info_spat%b_perc(2),errorflag)!
            case (3)
                !call write_matrices(trim(path)//'out_'//trim(sim%eff_rete_fn)//'.asc',info_spat%eff_rete,errorflag)! - RR
                call write_grid(trim(path)//'out_'//trim(sim%id_irr_meth_fn)//'.asc',info_spat%irr_meth_id,errorflag)!
                call write_grid(trim(path)//"out_am_percI.asc",info_spat%a_perc(1),errorflag)!
                call write_grid(trim(path)//"out_am_percII.asc",info_spat%a_perc(2),errorflag)!
                call write_grid(trim(path)//"out_bm_percI.asc",info_spat%b_perc(1),errorflag)!
                call write_grid(trim(path)//"out_bm_percII.asc",info_spat%b_perc(2),errorflag)!
            case (4)
                call write_grid(trim(path)//'out_'//trim(sim%irr_units_fn)//'.asc',info_spat%irr_unit_id,errorflag)!
                call write_grid(trim(path)//'out_'//trim(sim%eff_irr_fn)//'.asc',info_spat%eff_met,errorflag)!
                call write_grid(trim(path)//'out_'//trim(sim%eff_net_fn)//'.asc',info_spat%eff_net,errorflag)!
                call write_grid(trim(path)//'out_'//trim(sim%id_irr_meth_fn)//'.asc',info_spat%irr_meth_id,errorflag)!
                call write_grid(trim(path)//"out_am_percI.asc",info_spat%a_perc(1),errorflag)!
                call write_grid(trim(path)//"out_am_percII.asc",info_spat%a_perc(2),errorflag)!
                call write_grid(trim(path)//"out_bm_percI.asc",info_spat%b_perc(1),errorflag)!
                call write_grid(trim(path)//"out_bm_percII.asc",info_spat%b_perc(2),errorflag)!
            case default
        end select

        !meteo weigths
        do k=1,size(info_spat%weight_ws)!
            write(meteo_stringa,*) k!
            call write_grid(trim(path)//trim(adjustl('out_meteo_'//trim(adjustl(meteo_stringa))//'.asc')), &!
                & info_spat%weight_ws(k),errorflag)!
        end do!
        !!
    end subroutine write_init_grids!
    
    subroutine init_irrigation_units(domain_map,irr_units_map,eff_net,irr_units_tbl,wat_src_tbl,pars,h_met,debug)!
        !allocazione e inizializzazione della variabile IU!
        implicit none!
        type(grid_i),intent(in)::domain_map,irr_units_map!
        type(grid_r),intent(in)::eff_net,h_met!
        type(parameters),intent(in)::pars!
        type(irr_units_table),dimension(:),allocatable,intent(out)::irr_units_tbl
        type(water_sources_table),dimension(:),intent(inout)::wat_src_tbl!
        logical,intent(in)::debug
        !!
        integer::i,free_unit,error_flag,ios!
        integer:: strlen
        character(len=999) :: str                            ! File string
        character(len=1),parameter :: delimiter = achar(9)   ! Delimiter: horizontal tab
        integer::cols, n_irr_units                                 ! File columns (cols) and rows (n_bac = irrigation district number)
        !!
        strlen = 999
        ! Read "<irrdistr_list>.txt"
        call seek_un(error_flag,free_unit)!
        open(free_unit,file=trim(pars%sim%watsour_path)//trim(pars%sim%irrdistr_list_fn),action="read",iostat=ios)!
        if(ios/=0)then!
            print *, "Error opening file ", trim(pars%sim%irrdistr_list_fn)," connected to unit ", free_unit, &
                & " iostat=", ios, '. Execution will be aborted...'
            stop
        end if!
        
        read(free_unit,'(a)') str       ! Read file header (considering also delimiters and spaces: it needs formatting)
        strlen = len(trim(str))
        
        cols = 1                        ! Count the number of delimiters in the first line
        do i = 1, strlen
            if (str(i:i) == delimiter) then
                cols = cols+1
            end if
        end do
        if (cols /= 3) then
            print *, 'Error in file formatting: ', trim(pars%sim%irrdistr_list_fn), ' has ', cols, &
                & ' columns instead of 3. Execution will be aborted...'
            stop
            !!! TODO: %AB%: add wat_shift reading
        end if

        ! Count the number of rows
        n_irr_units = 0
        do
            read (free_unit,*,iostat=ios) str
            if (ios/=0) exit
            n_irr_units = n_irr_units + 1
        end do

        allocate(irr_units_tbl(n_irr_units))          ! Allocate variable with irrigation district number
        
        rewind(free_unit)                ! Rewind of file to store its content
        read(free_unit,*) str            ! Skip header line
        do i = 1, n_irr_units
!~             read(free_unit,*,iostat=ios)IU(i)%id,IU(i)%f_explore,IU(i)%f_un_priv,IU(i)%wat_shift!
            read(free_unit,*,iostat=ios)irr_units_tbl(i)%id,irr_units_tbl(i)%f_explore,irr_units_tbl(i)%f_un_priv!   
            irr_units_tbl(i)%wat_shift = 1      ! set wat_shift to 1 by default
            
            ! Calculate the number of irrigable cells in a day by counting the cells belonging to each irrigation units
            ! The value is yearly updated because landuse can change
            ! The value if rounded to the higher integer
            irr_units_tbl(i)%n_cells=count(irr_units_map%mat==irr_units_tbl(i)%id .and. domain_map%mat/=domain_map%header%nan)
            
            ! Calculate the maximum number of cells processed in a day
            irr_units_tbl(i)%n_day= ceiling (irr_units_tbl(i)%n_cells / irr_units_tbl(i)%wat_shift)
            ! set n_day max to n_cells
            if (irr_units_tbl(i)%n_day > irr_units_tbl(i)%n_cells) irr_units_tbl(i)%n_day = irr_units_tbl(i)%n_cells
            
            ! calculate di average network efficiency: SUM(net_eff)/ n_cells    
            irr_units_tbl(i)%int_distr_eff=sum(eff_net%mat, irr_units_map%mat==irr_units_tbl(i)%id .and. domain_map%mat/=domain_map%header%nan)/irr_units_tbl(i)%n_cells
            ! f2003 compatibility error
            irr_units_tbl(i)%int_distr_eff=merge(0.d0,irr_units_tbl(i)%int_distr_eff,irr_units_tbl(i)%int_distr_eff/=irr_units_tbl(i)%int_distr_eff)
            ! calculate the average irrigation height
            irr_units_tbl(i)%h_irr_mean=sum(h_met%mat,irr_units_map%mat==irr_units_tbl(i)%id .and. domain_map%mat/=domain_map%header%nan)/irr_units_tbl(i)%n_cells
            if(ios/=0)then
                print *, 'Error reading file ', trim(pars%sim%irrdistr_list_fn), ' in line ', n_irr_units, '. Execution will be aborted...'
                stop
            end if!
        end do!
        close(free_unit)

        ! init the reference to the water sources table
        do i=1,size(wat_src_tbl)!
            wat_src_tbl(i)%irr_unit_idx = get_value_index(irr_units_tbl%id,wat_src_tbl(i)%id_irr_unit)!
        end do
        
        ! init the discharges to zeros
        irr_units_tbl%q_day=0.;irr_units_tbl%q_nom=0.!
        irr_units_tbl%q_act_fld(1)=0.;irr_units_tbl%q_act_fld(2)=0.;irr_units_tbl%q_act_fld(3)=0.;irr_units_tbl%q_act_fld(4)=0.
        irr_units_tbl%q_pot_fld(1)=0.;irr_units_tbl%q_pot_fld(2)=0.;irr_units_tbl%q_pot_fld(3)=0.;irr_units_tbl%q_pot_fld(4)=0.
        irr_units_tbl%n_irrigated_cells=0!
        irr_units_tbl%q_un_priv=0!

        if(debug .eqv. .true.)then!
            call seek_un(error_flag,free_unit)!
            open(free_unit,file=(trim(pars%sim%path)//"out_"//trim(pars%sim%watsources_fn)),action="write")!
            write(free_unit,*)'distr_id; source_code; source_type; flow_ratio; distr_column; watsour_column'!
            do i=1,size(wat_src_tbl)!
                write(free_unit,*)wat_src_tbl(i)%id_irr_unit,'; ',wat_src_tbl(i)%id_wat_src,'; ',wat_src_tbl(i)%type_id, &!
                    & '; ',wat_src_tbl(i)%duty_frc,'; ',wat_src_tbl(i)%irr_unit_idx,'; ',wat_src_tbl(i)%wat_src_idx!
            end do!
            close(free_unit)!
            
            !%AB%: add duration
            call init_cell_output_file(free_unit,trim(pars%sim%path)//"irrdistricts_log.csv",&!
                & 'distr_id; cells_no; wat_shift; cell_day; conv_eff; havmet; privatewells')!
            do i=1,n_irr_units!
                write(free_unit,*)irr_units_tbl(i)%id,'; ',irr_units_tbl(i)%n_cells,&
                    & '; ',irr_units_tbl(i)%wat_shift,'; ',irr_units_tbl(i)%n_day,'; ',&!
                    & '; ',irr_units_tbl(i)%int_distr_eff,'; ',&!
                    & irr_units_tbl(i)%h_irr_mean,'; ',irr_units_tbl(i)%f_un_priv!
            end do
            close(free_unit)!
        end if
        
    end subroutine init_irrigation_units!

end module