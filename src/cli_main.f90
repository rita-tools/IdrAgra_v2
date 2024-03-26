program main!
    use mod_utility, only: lower_case, print_execution_time
    use mod_parameters!, only: read_parameters, make_default, print_parameters  ! variables and methods to read file *.xml paramters
    use mod_grid, only: write_grid, min_domain, bound      ! variables and methods to handle spatial input
    use mod_meteo!, only: meteo_series_length, read_meteo_parameters,
        ! close_meteo_file, read_meteo_data
    use mod_TDx_index                                     ! variables and methods to handle meteorological data

    use cli_crop_parameters, only: init_crop_phenology_pars, close_pheno_file, crop_pheno_info ! variables and methods to handle crop parameters
    use cli_watsources                                                          ! variables and methods to handle irrigation units water supply
    use cli_simulation_manager                                                  ! simulation manager: control IO and daily cycle
    use cli_read_parameter

    use mod_system

    implicit none!

    type(parameters) :: xml                                     ! stores parameters of input file
    type(TDx_index) :: xml_TDx                                  ! stores parameters of input file for TDX indices
    type(bound)::boundaries                                     ! stores domain extension (rectangle)
    type(spatial_info)::info_spat                               ! stores spatial inputs
    type(water_sources_table),dimension(:),allocatable::watsour ! stores parameters of "<watsources>.txt" file
    real(dp), dimension(9,3,4):: tab_CN2, tab_CN3               ! CN implementation
    type(source_info)::info_sources                             ! stores water sources series
    type(meteo_info),dimension(:),allocatable::info_meteo       ! stores meteorological series data
    type(crop_pheno_info),dimension(:),allocatable::info_pheno  ! stores phenological parameters series
    type(soil2_rice)::theta2_rice                               ! stores soil parameters data for paddy rice fields
    !!
    integer,dimension(8)::t_start,t_stop!
    integer :: errorflag!
    character(len = 255) :: filename = 'idragra_parameters.txt'
    character(len=50) :: arg
    integer :: i
    ! default seettings
    logical :: debug = .false.
    logical :: summary = .false.
    logical :: showpreview = .false.

    character :: delimiter

    
    ! get options
    do i = 1, iargc()
        call getarg(i, arg)
        if (arg(1:1)=='-') then
            call lower_case(arg)
            select case (arg)
                case ('-verbose', '-v')         ! set debug mode
                    debug = .true.
                case ('-summary','-s')          ! set if only irrigation is printed
                    summary = .true.
                case ('-help', '-h')            ! show help string
                    call show_help()
                case ('-filename', '-f')        ! set filename
                    call getarg(i+1, arg)
                    filename = arg
                case ('-preview', '-p')         ! show all parameters after reading configuration file
                    showpreview = .true.
                case ('-default', '-d')         ! print all default value
                    call make_default(xml, xml_TDx)
                    print *, '=== DEFAULT VALUES ==='
                    call print_parameters(xml,xml_TDx)
                    print *, '=== END DEFAULT ==='
                    stop
                case default                    ! all other cases ... !
                    print *, 'Not supported option <',trim(arg),'>'
                    call show_help()
            end select
        end if
    end do

	! TODO: check if there are other strings without '-', as it can be related to a skipped -f option

    ! set up environment options
    call get_environment_variable('DELIMITER',delimiter)
    call set_command('mkdir',delimiter)
    
    call print_header()
    ! Memorization of time in which simulation starts
    call date_and_time(values=t_start)

    ! Reads simulation input parameters
    call read_all_parameters(filename, xml, xml_TDx, ErrorFlag, debug)!

    ! check if output dir already exists
    call check_out_path(xml)

    if (showpreview .eqv. .true.) then
        print *, '=== PREVIEW ==='
        call print_parameters(xml, xml_TDx)
        print *, '=== END PREVIEW ==='
    end if
    print*,'Simulation parameters have been set'

    ! Evaluates simulation length [days] by counting meteorological files' lines
    call meteo_series_length(xml%sim, debug)!

    ! Defines simulation domain boundaries
    call min_domain(trim(xml%sim%input_path)//trim(xml%sim%domain_fn)//'.asc',boundaries,xml%sim)!
    print*,"Simulation domain boundaries have been set"!

    ! Checks if soil initial condition is an input or it is generated by running first year of simulation
    if(xml%sim%f_init_wc .eqv. .false.)then ! Generates soil initial condition

        ! Initializes info_spat and reads spatial files (*.asc); initializes tab_CN*
        call read_spatial_info(info_spat,boundaries,xml%sim,tab_CN2,tab_CN3,theta2_rice,xml%irr%met)!
        print*,"Variable info_spat has been initialized"!
        if (debug .eqv. .true.) then
            ! Print spatial data matrices in files out_*.asc
            call write_init_grids(info_spat,xml%sim%mode,xml%sim%path,xml%sim)!
            print*,"Input maps have been printed"!
        end if

        print *, '=== INITIAL CONDITION ==='
        ! Initializes info_meteo matrices with meteorological data
        call read_meteo_parameters(xml%sim,info_meteo,debug)!
        print*, 'Variable "info_meteo" has been initialized'

        ! Initializes info_pheno matrices (by associating file units to files)
        call init_crop_phenology_pars(xml%sim,info_pheno,info_meteo, debug)!
        print*, 'Variable "info_pheno" has been initialized'

        ! Initializes watsources and info_sources matrices
        if (xml%sim%mode == 1) then                 ! USE mode
            call init_water_sources_duty(xml,watsour,info_sources,info_meteo)!
            print*, 'Variables "watsour" and "info_sources" have been initialized'
        end if
       
        ! Soil-crop water balance algorithm
        call simulation_manager(xml, xml_TDx, info_spat, watsour, info_sources, info_meteo, &!
            & info_pheno, tab_CN2, tab_CN3, theta2_rice, 1, boundaries, debug, summary)!

        ! Prints initial condition values
        call write_grid(trim(xml%sim%path)//'IC_thetaI.asc',info_spat%theta(1)%old,errorflag)
        call write_grid(trim(xml%sim%path)//'IC_thetaII.asc',info_spat%theta(2)%old,errorflag)

        ! Closes input files
        if (xml%sim%mode == 1) call close_water_sources_dudy(info_sources,xml)        ! USE mode
        call close_meteo_file(info_meteo)!
        call close_pheno_file(info_pheno)!
        print *, '=== INITIAL CONDITION SET ==='
    end if!
    
    ! Initializes info_spat and reads spatial files (*.asc); initializes tab_CN*
    call read_spatial_info(info_spat,boundaries,xml%sim,tab_CN2,tab_CN3,theta2_rice,xml%irr%met)!
    print*,"Variable info_spat has been initialized"!
    if (debug .eqv. .true.) then
        ! Print spatial data matrices in files out_*.asc
        call write_init_grids(info_spat,xml%sim%mode,xml%sim%path,xml%sim)!
        print*,"Input maps have been printed"!
    end if
    
    ! %EAC%: add filter to limit water table to evaporative layer
    where(info_spat%wat_tab%mat<xml%depth%ze_fix .and. &
        info_spat%wat_tab%mat/=info_spat%wat_tab%header%nan)
        info_spat%wat_tab%mat=xml%depth%ze_fix
    end where

    ! Initializes info_meteo matrices with meteorological data
    call read_meteo_parameters(xml%sim,info_meteo,debug)!
    print*, 'Variable "info_meteo" has been initialized'

    ! Initializes info_pheno matrices (by associating file units to files)
    call init_crop_phenology_pars(xml%sim,info_pheno,info_meteo, debug)!
    print*, 'Variable "info_pheno" has been initialized'

    ! Initializes watsources and info_sources matrices
    if (xml%sim%mode == 1) then         ! %AB% USE mode
        call init_water_sources_duty(xml,watsour,info_sources, info_meteo)!
        print*, 'Variables "watsources" and "info_sources" have been initialized'
    end if

    ! Soil-crop water balance algorithm
    print*,"=== SIMULATION ==="!
    call simulation_manager(xml,xml_TDx,info_spat,watsour,info_sources,info_meteo,info_pheno, tab_CN2, tab_CN3, &!
        & theta2_rice,xml%sim%sim_years,boundaries, debug, summary)!

    ! Closes input files
    if (xml%sim%mode == 1) call close_water_sources_dudy(info_sources,xml)        ! USE mode
    call close_meteo_file(info_meteo)!
    call close_pheno_file(info_pheno)!

    ! Memorization of time in which simulation ends
    call date_and_time(values=t_stop)!

    ! Prints simulation computing time
    call print_execution_time(t_start, t_stop)!
    !!
end program main!

subroutine show_help()
    print *, 'List of options:'
    print *, '-h, -help: return some helpful information'
    print *, '-d, -default: print default input values'
    print *, '-p, -preview: show all parameters after reading configuration files'
    print *, '-t, -teta, -theta: print soil water content at the end of each year'
    print *, '-v, -verbose: print all outputs'
    print *, '-s, -summary: print only irrigation - both annual and periodical - outputs'
    print *, '-f, -filename <parameters> : use file <parameters> to pass settings'
    print *, '====================='
    ! stop execution because of a possible error !
    stop
end subroutine show_help

subroutine print_header()
    ! use makefile macros to set version and compilation timestamp
    print *, '============== IDRAGRA - IDRologia AGRAria ================'
    print *, 'contact: claudio.gandolfi@unimi.it'
    print *, 'code version: ',GIT_VERSION
    print *, 'compiled on: ',COMP_DATE
    print *, 'source available on https://github.com/rita-tools/IdrAgra_v2'
    print *, '============================================================'
    print *, ''
end subroutine print_header

subroutine check_out_path(xml)
    use mod_parameters
    type(parameters), intent(in) :: xml
    logical :: dir_exists
    character :: delimiter
    ! check if the dir exists or make a new dir
   

end subroutine

! TODO: replace with direct initialization of the variable
subroutine make_default(xml, xml_dtx)
    use mod_parameters
    use mod_TDx_index
    implicit none
    type(parameters), intent(inout) :: xml
    type(TDx_index), intent(inout) :: xml_dtx
    
    xml%sim%path =                   '.\\sim_results\\'
    xml%sim%input_path =             '.\\spatial_data\\'
    xml%sim%initial_condition =      xml%sim%input_path
    xml%sim%final_condition =        xml%sim%path
    xml%sim%domain_fn =                 'domain'
    xml%sim%thetaI_FC_fn =            'thetaI_FC'
    xml%sim%thetaII_FC_fn =           'thetaII_FC'
    xml%sim%thetaI_WP_fn =            'thetaI_WP'
    xml%sim%thetaII_WP_fn =           'thetaII_WP'
    xml%sim%thetaI_r_fn =             'thetaI_r'
    xml%sim%thetaII_r_fn =            'thetaII_r'
    xml%sim%thetaI_SAT_fn =           'thetaI_SAT'
    xml%sim%thetaII_SAT_fn =          'thetaII_SAT'
    xml%sim%slope_fn =            'slope'
    xml%sim%dren_fn =                'hydr_cond'
    xml%sim%hydr_group_fn =             'hydr_group'
    xml%sim%ksat_I_fn =              'ksat_I'
    xml%sim%ksat_II_fn =             'ksat_II'
    xml%sim%n_I_fn =                 'n_I'
    xml%sim%n_II_fn =                'n_II'
    xml%sim%thetaI_0_fn =             'IC_thetaI'
    xml%sim%thetaII_0_fn =            'IC_thetaII'
    xml%sim%thetaI_end_fn =           'FC_thetaI'
    xml%sim%thetaII_end_fn =          'FC_thetaII'
    xml%sim%soil_prop_x_rice_fn =   'rice_soilparam.txt'
    xml%sim%eff_irr_fn =          'appl_eff'
    xml%sim%eff_net_fn =            'conv_eff'
    xml%sim%irr_units_fn =           'irr_units'
    xml%sim%id_irr_meth_fn =       'irr_meth'
    xml%sim%wat_table_fn =            'waterdepth'
    xml%sim%ParRisCap_a3_fn =        'CapRisePar_a3'
    xml%sim%ParRisCap_a4_fn =        'CapRisePar_a4'
    xml%sim%ParRisCap_b1_fn =        'CapRisePar_b1'
    xml%sim%ParRisCap_b2_fn =        'CapRisePar_b2'
    xml%sim%ParRisCap_b3_fn =        'CapRisePar_b3'
    xml%sim%ParRisCap_b4_fn =        'CapRisePar_b4'
    xml%sim%soiluse_fn =             'soiluse'
    xml%sim%meteoweight_fn =         'meteo'
    xml%sim%shapearea_fn =           'shapearea'

    xml%sim%meteo_path =             '.\\meteo_data\\'
    xml%sim%ws_list_fn =         'weather_stations.dat'
    xml%sim%pheno_path =             '.\\crop_series\\'
    xml%sim%pheno_root =             'pheno_'
    xml%sim%irr_met_path =           '.\\irrmeth_data\\'
    xml%sim%irr_met_list_fn =       'irrmethods.txt'
    xml%sim%watsour_path =           '.\\watsour_data\\'
    xml%sim%watsources_fn =    'watsources.txt'
    xml%sim%mon_sources_i_div_fn =    'monit_sources_i.txt'
    xml%sim%mon_sources_ii_div_fn =   'monit_sources_ii.txt'
    xml%sim%int_reuse_div_fn =          'int_reuse.txt'
    xml%sim%cr_sources_list_fn =        'cr_sources.txt'
    xml%sim%irrdistr_list_fn =          'irr_districts.txt'
    xml%sim%sched_irr_fn =     'scheduled_irrigation.txt'
    xml%sim%step_out =                  0
    xml%sim%mode =                   2
    xml%sim%f_soiluse =                .false.
    xml%sim%h_prec_lim =              5.0
    
    ! Initializes start_simulation and end_simulation (to be overwritten)
    xml%sim%start_simulation%day   = 29
    xml%sim%start_simulation%month = 02
    xml%sim%start_simulation%year  = 1600
    xml%sim%start_simulation%doy = calc_doy(xml%sim%start_simulation%day, xml%sim%start_simulation%month, &
        & xml%sim%start_simulation%year)
    xml%sim%end_simulation%doy = xml%sim%start_simulation%doy

    allocate(xml%sim%clock(3))
    xml%sim%clock(1) =               10
    xml%sim%clock(2) =              100
    xml%sim%clock(3) =               30
    xml%sim%weekday =                 1
    xml%sim%f_init_wc =                 .false.
    xml%sim%rand_seed =             -999 ! Blank
    xml%sim%rand_symmetry =         .true.
    xml%sim%sowing_range =                0
    xml%sim%imax =                    1
    xml%sim%jmax =                    1
    xml%sim%cell_size =                  250
    xml%sim%x0 =                      1.0
    xml%sim%y0 =                      1.0
    xml%sim%n_lus =               1
    xml%sim%n_crops =               1
    xml%sim%f_cap_rise =             .false.
    xml%sim%n_ws =                1
    xml%sim%n_voronoi =              1
    xml%sim%start_irr_season =       91
    xml%sim%end_irr_season =        304
    xml%sim%f_out_cells =          .false.
    xml%sim%lambda_cn =               0.2
    xml%sim%f_shapearea =        .false.
    xml%depth%ze_fix =                0.10
    xml%depth%zr_fix =                0.90
    xml%n_dotaz =                     0
    xml%ms_i%n_withdrawals =                    0
    xml%ms_ii%n_withdrawals =                   0
    xml%intreu%n_withdrawals =                  0
    xml%cr%n_withdrawals =                      0
    
    allocate(xml%sim%quantiles(2,2))
    xml%sim%quantiles = reshape((/0.575118, 0.472116, 8.026400, 7.706101/),(/2,2/))

    xml_dtx%mode =                    'none'
    xml_dtx%n =                       0
end subroutine make_default


subroutine print_parameters(xml,xml_dtx)
    use mod_parameters, only: parameters
    use mod_TDx_index, only: TDx_index
    implicit none
    type(parameters), intent(in) :: xml
    type(TDx_index), intent(in) :: xml_dtx
    
    print *,'OutputPath = ',  xml%sim%path
    print *,'MeteoPath = ', xml%sim%meteo_path
    print *,'MeteoFileName = ', xml%sim%ws_list_fn
    print *,'PhenoPath = ', xml%sim%pheno_path
    print *,'PhenoFileRoot = ', xml%sim%pheno_root
    print *,'IrrMethPath = ', xml%sim%irr_met_path
    print *,'IrrMethFileName = ', xml%sim%irr_met_list_fn
    print *,'WatSourpath = ', xml%sim%watsour_path
    print *,'InitialCondition = ', xml%sim%initial_condition
    print *,'FinalCondition = ', xml%sim%final_condition
    ! spatial inputs
    print *, 'domFileName = ', xml%sim%domain_fn
    print *, 'tetaI_FCFileName = ', xml%sim%thetaI_FC_fn
    print *, 'tetaII_FCFileName = ', xml%sim%thetaII_FC_fn
    print *, 'tetaI_WPFileName = ', xml%sim%thetaI_WP_fn
    print *, 'tetaII_WPFileName = ', xml%sim%thetaII_WP_fn
    print *, 'tetaI_rFileName = ', xml%sim%thetaI_r_fn
    print *, 'tetaII_rFileName = ', xml%sim%thetaII_r_fn
    print *, 'tetaI_SATFileName = ', xml%sim%thetaI_SAT_fn
    print *, 'tetaII_SATFileName = ', xml%sim%thetaII_SAT_fn
    print *, 'pendenzaFileName = ', xml%sim%slope_fn
    print *, 'drenFileName = ', xml%sim%dren_fn
    print *, 'hydr_grFileName = ', xml%sim%hydr_group_fn
    print *, 'ksat_IFileName = ', xml%sim%ksat_I_fn
    print *, 'ksat_IIFileName = ', xml%sim%ksat_II_fn
    print *, 'n_IFileName = ', xml%sim%n_I_fn
    print *, 'n_IIFileName = ', xml%sim%n_II_fn
    print *, 'tetaI_0FileName = ', xml%sim%thetaI_0_fn
    print *, 'tetaII_0FileName = ', xml%sim%thetaII_0_fn
    print *, 'pedologica_x_risoFileName = ', xml%sim%soil_prop_x_rice_fn
    print *, 'eff_metodoFileName = ', xml%sim%eff_irr_fn
    print *, 'eff_reteFileName = ', xml%sim%eff_net_fn
    print *, 'irr_unitsFileName = ', xml%sim%irr_units_fn
    print *, 'codice_metodoFileName = ', xml%sim%id_irr_meth_fn
    print *, 'soggiacenzaFileName = ', xml%sim%wat_table_fn
    print *, 'ParRisCap_a3FileName = ', xml%sim%ParRisCap_a3_fn
    print *, 'ParRisCap_a4FileName = ', xml%sim%ParRisCap_a4_fn
    print *, 'ParRisCap_b1FileName = ', xml%sim%ParRisCap_b1_fn
    print *, 'ParRisCap_b2FileName = ', xml%sim%ParRisCap_b2_fn
    print *, 'ParRisCap_b3FileName = ', xml%sim%ParRisCap_b3_fn
    print *, 'ParRisCap_b4FileName = ', xml%sim%ParRisCap_b4_fn
    print *, 'ShapeAreaFlag = ', xml%sim%f_shapearea
    print *, 'ShapeAreaFileName = ', xml%sim%shapearea_fn

    print *,'mode = ',  xml%sim%mode
    print *,'mflag = ',  xml%sim%step_out
    print *,'time = ', xml%sim%clock
    print *,'InitialThetaFlag = ',  xml%sim%f_init_wc
    print *,'FinalThetaFlag = ',  xml%sim%f_theta_out
    print *,'RadnSowDaysWind = ',  xml%sim%sowing_range
    print *,'SoilUsesNum = ', xml%sim%n_lus
    print *,'CapillaryFlag = ', xml%sim%f_cap_rise
    print *,'SoilUseVarFlag = ', xml%sim%f_soiluse
    print *,'MeteoStatWeightNum = ', xml%sim%n_ws
    print *,'MeteoStatTotNum = ', xml%sim%n_voronoi
    print *,'forecast_day = ', xml%sim%forecast_day
    print *,'zEvap = ', xml%depth%ze_fix
    print *,'zRoot = ', xml%depth%zr_fix
    print *,'LambdaCN = ', xml%sim%lambda_cn
    print *,'StartIrrSeason = ', xml%sim%start_irr_season
    print *,'EndIrrSeason = ', xml%sim%end_irr_season
    print *, 'DTxMode =',xml_dtx%mode
    print *, 'DTxNumXs=',xml_dtx%temp%n_ind
    print *, 'DTx_X=',xml_dtx%temp%x
    print *, 'DTxDeltaDate=',xml_dtx%temp%td
    print *, 'DTxDelayDays=',xml_dtx%temp%delay
    print *, 'DTxMinCard=',xml_dtx%n
end subroutine print_parameters
