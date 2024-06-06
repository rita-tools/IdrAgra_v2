module mod_parameters!
    use mod_utility, only: date, lower_case, seek_un, string_to_integers, calc_doy, split_date
    use mod_constants
    implicit none!

    type simulation!
        character(len=200) :: path                      ! path of output!
        character(len=200) :: meteo_path                ! path of meteo tables
        character(len=200) :: ws_list_fn                ! filename (with complete path) of weather stations list
        character(len=200) :: pheno_path                ! path where phenophases files are stored
        character(len=200) :: pheno_root                ! a common string to use to get folder
        character(len=200) :: irr_met_path              ! path to irrigation methods folder
        character(len=200) :: irr_met_list_fn           ! filename of a list of irrigation methods
        character(len=200) :: watsour_path              ! path to water sources folder
        character(len=200) :: watsources_fn             ! filename of water sources distribution
        character(len=200) :: mon_sources_i_div_fn      ! filename of monitored sources (i) daily diversions
        character(len=200) :: mon_sources_ii_div_fn     ! filename of monitored sources (ii) daily diversions
        character(len=200) :: int_reuse_div_fn          ! filename of internal reuse daily diversions
        character(len=200) :: cr_sources_list_fn        ! filename of collective runtime sources list
        character(len=200) :: irrdistr_list_fn          ! filename of public irrigation districts' list
        character(len=200) :: sched_irr_fn              ! filename of scheduled irrigation file
        character(len=200) :: input_path                ! path to input file
        character(len=200) :: initial_condition         ! path to initial condition input file
        character(len=200) :: final_condition           ! path to final condition input file
        character(len=255) :: domain_fn                 ! domain file name name
        character(len=255) :: thetaI_FC_fn
        character(len=255) :: thetaII_FC_fn
        character(len=255) :: thetaI_WP_fn
        character(len=255) :: thetaII_WP_fn
        character(len=255) :: thetaI_r_fn
        character(len=255) :: thetaII_r_fn
        character(len=255) :: thetaI_SAT_fn
        character(len=255) :: thetaII_SAT_fn
        character(len=255) :: slope_fn
        character(len=255) :: dren_fn
        character(len=255) :: hydr_group_fn
        character(len=255) :: ksat_I_fn
        character(len=255) :: ksat_II_fn
        character(len=255) :: n_I_fn
        character(len=255) :: n_II_fn
        character(len=255) :: thetaI_0_fn
        character(len=255) :: thetaII_0_fn
        character(len=255) :: thetaI_end_fn
        character(len=255) :: thetaII_end_fn
        character(len=255) :: soil_prop_x_rice_fn
        character(len=255) :: eff_irr_fn
        character(len=255) :: eff_net_fn
        character(len=255) :: irr_units_fn
        character(len=255) :: id_irr_meth_fn
        character(len=255) :: wat_table_fn
        character(len=255) :: ParRisCap_a3_fn
        character(len=255) :: ParRisCap_a4_fn
        character(len=255) :: ParRisCap_b1_fn
        character(len=255) :: ParRisCap_b2_fn
        character(len=255) :: ParRisCap_b3_fn
        character(len=255) :: ParRisCap_b4_fn
        character(len=255) :: soiluse_fn
        character(len=255) :: meteoweight_fn
        character(len=255) :: shapearea_fn

        integer :: step_out                         ! monthly output = 0, weekly output = 1, user defined = 2
        integer :: mode                             ! type of simulation:
                                                    ! 0 = without irrigation, 1 = need mode at field capacity, 2 = need mode at fixed volume, 3 = Use mode, 4 = scheduled irrigation
        logical :: f_soiluse                        ! land use change between years (true) otherwise false
        integer :: weekday                          ! day of the week of the weekly outputs (1 = Monday, 2 = Tuesday, ..., 7 = Sunday)
        integer,dimension(:),pointer::clock
        integer,dimension(:),pointer::intervals
        integer :: start_year                       ! first year of simulation
        integer :: sim_years                        ! number of years for the simulation
        integer :: meteo_years                      ! number of available weather time series
        integer,dimension(:),pointer :: year_step   ! number of days for each simulation years
        logical :: f_init_wc                        ! if true, use user provided values otherwise calculate from first year run
        logical :: f_theta_out                      ! id true, write final soil moisture condition
        integer :: rand_seed                        ! seed to initialize the random number generator
        logical :: rand_symmetry                    ! states if the randomization of emergence date is symmetric [-n,+n] or asymmetric [0, 2n]
        integer :: sowing_range                     ! range of variability of the sowing day [days]
        integer :: imax                             ! maximum number of rows in the raster map
        integer :: jmax                             ! maximum number of columns in the raster map
        real(dp) :: cell_size                       ! lateral dimension of the square cell of the raster map
        real(dp) :: x0                              ! x coordinate of the reference corner
        real(dp) :: y0                              ! y coordinate of the reference corner
        integer :: n_voronoi                        ! number of voronoi polygons
        integer :: n_lus                            ! total number of land uses
        integer :: n_crops                          ! max number of crops per year
        integer,dimension(:),pointer :: lu_list     ! list of ids if the simulated land uses
        integer,dimension(:),pointer :: no_lu_list  ! list of ids if the NOT simulated land uses
        logical :: f_cap_rise                       ! if true, capillary rise is calculated
        logical :: f_shapearea                      ! if true, use shapes area (for vectorialization)
        integer :: n_ws                             ! maximum number of weather station
        integer :: start_irr_season                 ! start irrigation season [doy]
        integer :: end_irr_season                   ! end irrigation season [doy]
        integer :: n_irr_meth                       ! number of irrigation methods
        logical :: f_out_cells                      ! if true, print outputs file for control points (aka cells)
        real(dp),dimension(:),pointer :: res_canopy ! plant resistance
        real(dp),dimension(:,:),pointer :: quantiles! extremes of variability of the k_sat (percolation model): 10th e 90th for 1st layer, 10th e 90th for 2nd layer
        real(dp) :: lambda_cn                       ! lambda parameters for curve number [-]
        REAL(dp) :: h_prec_lim = 5.0                ! minimum meaningful precipitation [mm]
        
        type(date) :: start_simulation              ! start simulation date
        type(date) :: end_simulation                ! end simulation date
        integer :: forecast_day = 5                 ! number of days to use to cumulate precipitation
    end type simulation!

    type layer_depth!
        real(dp) :: ze_fix                          ! depth of the evaporative (1st from top) layer [m]
        real(dp) :: zr_fix                          ! depth of the transpirative (2nd from top) layer [m] (used when crop is missing) [m]
    end type layer_depth!

    type par_method!
        integer :: id                               ! id of the irrigation method, linked to the irrigation input map
        real(dp) :: h_irr                           ! irrigation depth [mm]!
        real(dp) :: irr_th_ms                        ! threshold coefficient for the activation of irrigation from monitored source of need mode
        real(dp) :: irr_th_unm                       ! threshold coefficient for the activation of irrigation from unmonitored source of need mode
        real(dp) :: f_wet                           ! fraction of wetted soil [-]
        real(dp) :: a_min = 0.0D0                   ! min value of the "a" parameter for percolation model
        real(dp) :: a_max = 0.0D0                   ! max value of the "a" parameter for percolation model              
        real(dp) :: b_min = 0.0D0                   ! min value of the "b" parameter for percolation model
        real(dp) :: b_max = 0.0D0                   ! max value of the "b" parameter for percolation model
        
        ! irrigation losses estimation model: losses=a_loss+b_loss x wind_vel+c_loss x T_mean
        real(dp) :: a_loss
        real(dp) :: b_loss
        real(dp) :: c_loss

        integer :: f_interception                   ! if 1 (true), irrigation is above the canopy
        real(dp),dimension(:),pointer::freq!

        integer :: irr_starts = 1                    ! irrigation period can start from
        integer :: irr_ends  = 366                   ! irrigation period will end to

    end type par_method!

    type par_irr_event
        real(dp) :: f_w                                 ! wetted fraction                             
        type(par_method),dimension(:),pointer::met      ! irrigation method parameter
        real(dp) :: h_irr                               ! water application depth [mm]!
        real(dp) :: alpha_raw                           ! fraction of RAW at which the irrigation occurs
        ! for USE mode only
        real(dp) :: alpha_raw_uc                        ! fraction of RAW at which the irrigation occurs for unmonitored collective water sources
        real(dp) :: expl_fact                           ! exploration factor: multiply the number of cells potentially irrigable
    end type par_irr_event!

    type water_source!
        integer :: n_withdrawals                    ! number of withdrawals for each water sources
        logical :: f_exists                         ! if true, allocate the variable
        logical :: f_allocate_timeserie             ! if true, allocate the array
    end type water_source!

    type uc_act_rule                                      ! activation rule for unmonitored collective water sources
        integer :: id_well_group                          ! id of the group of wells
        real(dp) :: min_act_trsld                         ! minimum activation threshold
        real(dp), dimension(3) :: flow_rate               ! activation flow rate
        real(dp), dimension(3) :: act_trsld               ! activation threshold
    end type uc_act_rule!

    ! TODO: NOT USED
    !type pond
    !    real(dp)::slope_min
    !    real(dp)::q_max
    !    real(dp)::slope_max
    !    real(dp)::q_min
    !end type pond
    
    type parameters!
        type(simulation) :: sim!
        integer :: n_bac!
        type(layer_depth) :: depth!
        type(par_irr_event) :: irr!
        integer :: n_dotaz!
        type(water_source) :: ms_i!
        type(water_source) :: ms_ii!
        type(water_source) :: intreu!
        type(water_source) :: cr!
        type(uc_act_rule),dimension(:),pointer :: uc_act_rules!
        real(dp),dimension(:),pointer::f_eff_rain!
        real(dp),dimension(:),pointer::fet0!
    end type parameters!

    type water_sources_table
        integer :: id_irr_unit              ! id of the irrigation unit served by the water source
        character(len=10) :: id_wat_src     ! id of the water source which the water duty refers to
        integer :: type_id                  ! type of the water source, (1) monitored sources i, (2) monitored sources ii, (3) internal reuse, (4) collective unmonitored sources
        real :: duty_frc                    ! Rate of the water duty
        integer :: irr_unit_idx             ! column index in the irrigation units variable
        integer :: wat_src_idx              ! column index in the discharge file
    end type water_sources_table

    type monitored_sources_table
        ! store monitored sources parameters and daily data    
        real(dp),dimension(:),pointer::q_nom                ! list of nominal discharge [m3/s]
        character(len=10),dimension(:),pointer::wat_src_id  ! list of water sources id (name of max 10 characters)
        real(dp),dimension(:,:),pointer::q_daily            ! matrix of daily discharges (days x num water sources) [m3/s]
        integer::unit                                       ! unit associated to the input file
        type(date)::start                                   ! start date of the time serie
        type(date)::finish                                  ! end date of the time serie
    end type monitored_sources_table!

    type unmonitored_sources_table
        ! store unmonitored collective sources parameters    
        real(dp),dimension(:),pointer::q_max                ! list of maximum discharge value [m3/s]!
        real(dp),dimension(:),pointer::q_nom                ! list of nominal discharge value [m3/s]!
        character(len=10),dimension(:),pointer::wat_src_id  ! list of water sources id (name of max 10 characters)
        integer,dimension(:),pointer::unit                  ! unit associated to the input file
    end type unmonitored_sources_table!

    type source_info!
    ! collect all the different water sources
        type(monitored_sources_table)::mn_src_tbl1              ! monitored sources (i)
        type(monitored_sources_table)::mn_src_tbl2              ! monitored sources (ii)
        type(monitored_sources_table)::int_reuse_tbl            ! internal reuse
        type(unmonitored_sources_table)::unm_src_tbl            ! unmonitored water sources
    end type source_info!

    type irr_units_table!
        ! store all the data referred to th irrigation units
        integer::id                             ! irrigation unit id
        integer::n_cells                        ! number of cells belonging to the irrigation units
!        integer::n_irr                         ! number of irrigable cells, depending on phenology
        !integer::conta                         ! NOT USED cell counter with specific id
        integer::last_cell_id                   ! id of the latest irrigated cell
        real(dp)::int_distr_eff                 ! internal water distribution efficiency
        real(dp)::h_irr_mean                    ! average irrigation height in the irrigation unit
        real(dp),dimension(4)::q_act_fld        ! discharge available at cell
        real(dp),dimension(4)::q_pot_fld        ! potential discharge available at cell
        real(dp)::q_nom                         ! potential discharge
        real(dp)::q_day                         ! actual daily discharge
        real(dp)::q_trashed                     ! discharge lost at the end of the rotation or the end of the irrigation season
        real(dp)::q_surplus                     ! discharge available the following day
        real(dp)::q_un_priv                     ! discharge extracted from unmonitored private sources (e.g. farm wells)
        integer::n_irrigated_cells              ! number of cell irrigated in the current day
        integer::f_un_priv                      ! flag for activate irrigation from unmonitored private sources (e.g. farm wells)
        real(dp)::n_irrigable_cells             ! number of cell potentially irrigable calculated as v_nom/SUM(v_irr_cells)
        !real(dp)::a_raw                        ! NOT USED: activation threshold for monitored water sources
        !real(dp)::a_raw_priv                   ! NOT USED: activation threshold for private water sources
        real(dp)::f_explore                     ! multiplicative factor to set the number of cells potentially irrigable in a day
        !!! %AB% <------- new variable for the rotation
        real::wat_shift                     ! rotation length (in days or fraction of day)
        integer::n_day                      ! number of cells irrigable in a day (related to the rotation)!
    end type irr_units_table!

    type scheduled_irrigation
        integer::irr_unit_id               ! id of the irrigation unit
        integer::year                      ! year of irrigation
        integer::doy                       ! day of the year [1-->365/6]
        real(dp)::h_irr                    ! water depth in each irrigation
    end type
    
end module mod_parameters!
