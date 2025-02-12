module cli_simulation_manager!
! Contains subroutines:
! soilcropwater_balance
! write_daily_output        <- writes daily output for selected cells (files *.csv)
! write_monthly_output      <- writes monthly/periodic files (files *.asc)
! deallocate_all            <- deallocates balance matrices
! iniz_bil1                 <- initializes evaporative layer matrices
! iniz_bil2                 <- initializes transpirative layer matrices
! eq_bil1                   <- updates evaporative layer matrices (operator "=")
! eq_bil2                   <- updates transpirative layer matrices (operator "=")
! iniz_pheno                <- initializes phenological matrices
! iniz_meteo                <- initializes meteorological matrices
! matrici_meteo_pesate      <- spatializes meteorological data and calculates ET0
! eq_extensive              <- initializes hourly variables' matrices
! interventi_fabbisogni_fix <- calculates irrigation - NEED mode, fixed water supply to each cell
! interventi_fabbisogni_fc  <- calculates irrigation - NEED mode, water supply to reach field capacity
! scheduled_irrigation_duty <- calculates irrigation - scheduled irrigation
! b1_no_iter_eva            <- calculates of variables used in evaporation estimation
! inizializza_aggregati     <- initialization of crop-soil water balance variables
! x_wat                     <- calculates soil content thresholds as a function of root zone depth
! irr_xRiso                 <- calculates irrigation for paddy fieleds
! dim_bil1                  <- allocates and deallocates bil1 variable
! dim_bil2                  <- allocates and deallocates bil2 variable
! dim_bil_hour              <- allocates and deallocates bil_hour variable
! dim_meteo                 <- allocates and deallocates meteo variable
! dim_pheno                 <- allocates and deallocates pheno variable
! dim_wat                   <- allocates and deallocates wat variable
! interventi_consumi        <- calculates irrigation in USE simulation (mode==1)
! Contains functions:
! calc_interception              <- interception calculation
! Pioggia_Efficace          <- effective rainfall calculation
!
    use mod_utility, only: sp, dp, get_value_index, make_numbered_name, get_uniform_sample, days_x_month, calc_date, day_of_week
    use mod_parameters
    use mod_grid, only: read_grid, write_grid, print_mat_as_grid, overlay_domain, &
                        & bound, id_to_par
    use mod_evapotranspiration, only: ET_reference, calculateDLH
    use mod_meteo, only: meteo_info, meteo_mat, read_meteo_data
    use mod_runoff
    use mod_crop_soil_water
    use mod_crop_phenology, only: crop_pheno_info, crop_matrices
    use mod_TDx_index
    use mod_constants, only: tmax_time, tmin_time, pi, cost_fwEva
    use mod_common, only: wat_matrix, soil2_rice, hourly, unit_file_scratch
    use mod_irrigation
    
    use cli_watsources
    use cli_crop_parameters, only: read_all_crop_pars, destroy_infofeno_tab, check_pheno_parameters, k_cb_matrices
    use cli_save_outputs
    use cli_read_parameter
    
    implicit none!

    ! Assignment statement overloading interface
    ! Help to initialize and update derived types
    interface assignment(=) !eq_extensive
        module procedure eq_extensive!
    end interface!
    interface assignment(=) !eq_bil
        module procedure eq_wat_bal1,eq_wat_bal2,init_wat_bal1,init_wat_bal2!
    end interface!
    interface assignment(=) !iniz_meteo
        module procedure init_meteo!
    end interface!
    interface assignment(=) !iniz_pheno
        module procedure init_pheno!
    end interface!
    !
    contains!
    !
    subroutine simulation_manager(pars,pars_TDx,info_spat,wat_src_tbl,info_sources, info_meteo, info_pheno, tab_CN2, tab_CN3,&
        & theta2_rice, sim_years,boundaries, debug, summary) !
        implicit none!

        type(parameters),intent(in)::pars!
        type(TDx_index),intent(in)::pars_TDx!
        real(dp),dimension(:,:,:),intent(in):: tab_CN2, tab_CN3!
        integer,intent(in)::sim_years!
        type(bound),intent(in)::boundaries!
        logical,intent(in)::debug
        logical,intent(in)::summary
        type(soil2_rice),intent(in)::theta2_rice
        type(spatial_info),intent(inout)::info_spat!
        type(water_sources_table),dimension(:),intent(inout)::wat_src_tbl!
        type(source_info),intent(inout)::info_sources!
        type(meteo_info),dimension(:),intent(inout)::info_meteo!
        type(crop_pheno_info),dimension(:),intent(inout)::info_pheno!
        !!
        type(balance1_matrices)::wat_bal1,wat_bal1_old!
        type(balance2_matrices)::wat_bal2,wat_bal2_old!
        type(crop_pars_matrices)::pheno
        type(meteo_mat)::meteo!
        real(dp),dimension(info_spat%domain%header%imax,info_spat%domain%header%jmax)::out_cn_day!
        type(output_CN),dimension(info_spat%domain%header%imax,info_spat%domain%header%jmax)::out_cn!
        type(hourly)::wat_bal_hour!
        type(wat_matrix)::wat!
        type(output_table_list)::out_tbl_list!
        type(step_map)::stp_map!
        type(step_debug_map)::deb_map!
        type(annual_map)::yr_map
        type(annual_debug_map)::yr_deb_map
        type(yield_t)::yield
        type(irr_units_table),dimension(:),allocatable::irr_units      ! Allocated in mod_watsources
        type(scheduled_irrigation),dimension(:),allocatable::irr_sch ! Allocated in 'open_scheduled_irrigation' function
        type(crop_matrices)::crop_map
        !!
        integer:: unit_crop
        integer::i,j,k,y,current_year,doy,hour,z,w                           ! for cycles
        integer,dimension(info_spat%domain%header%imax,info_spat%domain%header%jmax)::irandom
        integer,dimension(info_spat%domain%header%imax,info_spat%domain%header%jmax)::dir_phenofases
        integer,dimension(info_spat%domain%header%imax,info_spat%domain%header%jmax,size(info_spat%weight_ws))::dir_meteo
        real(dp),dimension(info_spat%domain%header%imax,info_spat%domain%header%jmax,size(info_spat%weight_ws))::meteo_weight
        character(len=255),dimension(info_spat%domain%header%imax,info_spat%domain%header%jmax)::code_pmeteo
        real(dp),dimension(info_spat%domain%header%imax,info_spat%domain%header%jmax,pars%sim%n_irr_meth)::h_irr ! z depends on number of irrigation methods 
        real(dp),dimension(info_spat%domain%header%imax,info_spat%domain%header%jmax)::priv_irr!
        real(dp),dimension(info_spat%domain%header%imax,info_spat%domain%header%jmax)::coll_irr!
        integer,dimension(12)::days_in_yr!
        integer::error_flag!
        integer::xx,yy ! Test cells coordinates
        integer,dimension(info_spat%domain%header%imax,info_spat%domain%header%jmax)::iter1,iter2
        logical,parameter::out_asc=.true.
        logical,parameter::out_yearly=.true.
        character(len=33)::str
        !!
        real(dp),dimension(info_spat%domain%header%imax,info_spat%domain%header%jmax)::irr_loss ! Irrigation application losses
        real(dp),dimension(info_spat%domain%header%imax,info_spat%domain%header%jmax)::alpha_ms_map, alpha_unm_map
        real(dp),dimension(info_spat%domain%header%imax,info_spat%domain%header%jmax)::fw_irr
        real(dp),dimension(info_spat%domain%header%imax,info_spat%domain%header%jmax)::fw_day, fw_old! fw daily updated
        real(dp),dimension(info_spat%domain%header%imax,info_spat%domain%header%jmax)::fc ! cover fraction - %RR%
        real(dp),dimension(info_spat%domain%header%imax,info_spat%domain%header%jmax)::a_loss, b_loss, c_loss, f_interception ! application losses model
        real(dp),dimension(info_spat%domain%header%imax,info_spat%domain%header%jmax)::h_irr_sum, h_bypass
        !! TDx
        real(dp),dimension(info_spat%domain%header%imax,info_spat%domain%header%jmax,pars_TDx%temp%n_ind)::tot_deficit      ! TDx sum
        integer,dimension(2)::unit_deficit
        integer::n_day,n_week!
        real(dp),dimension(info_spat%domain%header%imax,info_spat%domain%header%jmax)::TD!
        type(unit_file_scratch),dimension(:),allocatable::unit_Dxi!
        integer::cont_td    ! cycles
        character(len=33)::str_td!
        character(len=255)::str_delete, upfilename
        character(len=55)::s_year,s_years,s_doy
        logical :: file_exists
        !! percolation model
        integer,dimension(info_spat%domain%header%imax,info_spat%domain%header%jmax)::day_from_irr ! days past from the latest irrigation event
        real(dp),dimension(info_spat%domain%header%imax,info_spat%domain%header%jmax,2)::esp_perc  ! exponent of the percolation model
        integer:: shift_days
        integer::jul_day, day_cal, month_cal, year_cal
        integer::tmax_d, tmin_d, time
        ! to account for skipped years
        integer:: s_meteostat, s_step
        real(dp):: value
        !CHARACTER(LEN=20)::fc_name
        TYPE(grid_r)::pheno_grd
        !
        real(dp)::lat_sum = 0.0D0
        integer::lat_num = 0
        real(dp)::lat_mean
        real(dp)::DLH
        
        pheno_grd = info_spat%domain
        
        ! init irrigation time limits
        info_spat%irr_starts = info_spat%domain
        info_spat%irr_ends = info_spat%domain
        
        ! init maximum pond
        info_spat%h_maxpond=info_spat%slope
        info_spat%h_maxpond%mat = 10000.0D0
                
        ! spread irrigation variables
        if (pars%sim%mode>0) then
            info_spat%irr_ends%mat=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%irr_ends)
            info_spat%irr_starts%mat=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%irr_starts)
            info_spat%h_maxpond%mat=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%h_maxpond)
        end if

        ! Variables allocation
        call allocate_all (stp_map, yr_map, deb_map, yr_deb_map, wat_bal1, wat_bal1_old, wat_bal2, wat_bal2_old, wat_bal_hour, meteo, wat, pheno, &
            & info_spat%domain%header%imax,info_spat%domain%header%jmax, info_spat%domain%mat, debug)
        !!
        ! dir_phenofases: for each cell the appropriate meteorological station is selected
        ! (by linking to its progressive number in meteorological stations list)
        dir_phenofases(:,:) = int(info_spat%weight_ws(1)%mat(:,:))
        do j=1,size(info_spat%domain%mat,2)
            do i=1,size(info_spat%domain%mat,1)
                if(info_spat%domain%mat(i,j)/=info_spat%domain%header%nan) then
                    code_pmeteo(i,j) = make_numbered_name(dir_phenofases(i,j),".dat")
                    dir_phenofases(i,j) = get_value_index(info_meteo%filename, code_pmeteo(i,j))
                end if
            end do
        end do
        !!
        ! dir_meteo: for each cell, the appropriate meteorological stations are selected
        ! (by changing meteorological station ID to its progressive number in meteorological stations list)
        do k=1,size(dir_meteo,3)
            dir_meteo(:,:,k) = int(info_spat%weight_ws(k)%mat(:,:))
            do j=1,size(info_spat%domain%mat,2)
                do i=1,size(info_spat%domain%mat,1)
                    if(info_spat%domain%mat(i,j)/=info_spat%domain%header%nan) then
                        code_pmeteo(i,j) = make_numbered_name(dir_meteo(i,j,k),".dat")
                        dir_meteo(i,j,k) = get_value_index(info_meteo%filename,code_pmeteo(i,j))
                        if (dir_meteo(i,j,k)==0) then
                            print*,"The weather station", int(info_spat%weight_ws(k)%mat(i,j)), &
                                & "cannot be found in the &
                                & meteorological stations list. Execution will be aborted..."
                            stop
                        end if
                        meteo_weight(i,j,k) = info_spat%weight_ws(k)%mat(i,j) - int(info_spat%weight_ws(k)%mat(i,j))!
                    end if
                end do
            end do
        end do
       
        ! Creates scratch files for TDx calculation
        scratch_td: if(trim(pars_TDx%mode)/="none")then
            unit_deficit=[(maxval(pars%sim%year_step(:))/pars_TDx%temp%td)+1,4]
            allocate(unit_Dxi(pars_TDx%temp%n_ind))
            do cont_td=1,pars_TDx%temp%n_ind
                write(str_td,*)pars_TDx%temp%x(cont_td)
                str_td="td"//trim(adjustl(str_td))                          ! Creates TDx strings (TD10, TD30, etc.) for scratch files
                call init_TDx(unit_deficit,info_spat%domain,pars%sim%path,trim(str_td))
                allocate(unit_Dxi(cont_td)%dxi(pars_TDx%temp%x(cont_td)))    ! Allocation of scratch units
            end do
        end if scratch_td
        !!
        n_day=0   ! Daily count initialization
        h_irr = 0.
        !!
        ! TODO Explore if is possible to have outputs on a specific day of the week (utility:day_of_week)
        !!
        fw_day = cost_fwEva ! Initialization of fw, set to 1
        fw_old = fw_day
        select case (pars%sim%mode)
            case (0)
                info_spat%irr_meth_id=info_spat%domain
                info_spat%irr_meth_id%mat=1                      ! Methods matrix set to 1
                info_spat%h_meth=info_spat%domain
                info_spat%h_meth%mat=0                           ! Qwat == 0 for each irrigation method
                !
                alpha_ms_map=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%irr_th_ms)  ! Spreads irrigation threshold for each irrigation method
                alpha_unm_map=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%irr_th_unm)! Spreads irrigation threshold for each irrigation method
                !
                f_interception=info_spat%domain%mat
                f_interception=1                                                   ! Flag interception set to 1
            case (1)
                info_spat%h_meth=info_spat%domain
                info_spat%h_meth%mat=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%h_irr)  ! Spreads irrigation height for each irrigation method
                call init_irrigation_units(info_spat%domain,info_spat%irr_unit_id,info_spat%eff_net,irr_units,wat_src_tbl,&
                    &pars,info_spat%h_meth, debug)!
                alpha_ms_map=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%irr_th_ms)  ! Spreads irrigation threshold for each irrigation method
                alpha_unm_map=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%irr_th_unm)! Spreads irrigation threshold for each irrigation method
                fw_irr=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%f_wet)     ! Spreads wetted fraction for each irrigation method
                a_loss=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%a_loss)    ! Spreads irrigation application loss pars for each irrigation method
                b_loss=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%b_loss)
                c_loss=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%c_loss)
                f_interception=id_to_par(info_spat%irr_meth_id,dble(pars%irr%met(:)%f_interception))
            case (2)
                alpha_ms_map=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%irr_th_ms)
                alpha_unm_map=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%irr_th_unm)
                fw_irr=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%f_wet)
                f_interception=id_to_par(info_spat%irr_meth_id,dble(pars%irr%met(:)%f_interception))
            case (3)
                info_spat%h_meth=info_spat%domain
                info_spat%h_meth%mat=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%h_irr)
                alpha_ms_map=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%irr_th_ms)
                alpha_unm_map=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%irr_th_unm)
                fw_irr=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%f_wet)
                a_loss=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%a_loss)
                b_loss=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%b_loss)
                c_loss=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%c_loss)
                f_interception=id_to_par(info_spat%irr_meth_id,dble(pars%irr%met(:)%f_interception))
            case (4)
                info_spat%h_meth=info_spat%domain
                info_spat%h_meth%mat=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%h_irr)
                alpha_ms_map=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%irr_th_ms)
                alpha_unm_map=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%irr_th_unm)
                fw_irr=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%f_wet)
                a_loss=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%a_loss)
                b_loss=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%b_loss)
                c_loss=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%c_loss)
                f_interception=id_to_par(info_spat%irr_meth_id,dble(pars%irr%met(:)%f_interception))          
            case default
        end select
        !!
        ! Yearly simulation cycle
        year_cycle: do y=1,sim_years!
            if (pars%sim%start_simulation%doy >= info_meteo(1)%start%doy + sum(pars%sim%year_step(1:y))) then
                ! Skip all year if not simulated
                do s_meteostat = 1, pars%sim%n_voronoi   ! read in call xmeteo
                    do s_step = 1, pars%sim%year_step(y)
                        read (info_meteo(s_meteostat)%unit, *) value
                    end do
                end do!
                cycle
            else if (pars%sim%start_simulation%doy > &
                & info_meteo(1)%start%doy + sum(pars%sim%year_step(1:(y-1)))) then
                ! Eventually skip some days
                do s_meteostat = 1, pars%sim%n_voronoi
                    do s_step = 1, &
                        & (pars%sim%start_simulation%doy - (info_meteo(1)%start%doy + sum(pars%sim%year_step(1:(y-1)))))
                        read (info_meteo(s_meteostat)%unit, *) value
                    end do
                end do!
            end if
            
            n_week=0!
            
            ! TODO: check if effectively required
            current_year = pars%sim%start_year+y-1
            if (info_meteo(1)%start%month > 1) then
                
                write (s_year,*) current_year
                write (s_years,*) current_year+1
                write (s_years,*) adjustr(trim(s_year))//'-'//adjustl(trim(s_years))
            else
                write (s_year,*) current_year
                write (s_years,*) current_year
            end if
            
            ! Yearly soil use: change soil use map each year at the beginning
            if (pars%sim%f_soiluse .eqv. .true.) then
                info_spat%domain = info_spat%backup_domain
                info_spat%domain%mat = info_spat%backup_domain%mat
                call read_grid (trim(pars%sim%input_path)//trim(pars%sim%soiluse_fn)//'_'//trim(adjustl(s_years))//".asc",&
                    & info_spat%soil_use_id,pars%sim,boundaries)
                if (minval(info_spat%soil_use_id%mat,info_spat%soil_use_id%mat/=info_spat%soil_use_id%header%nan) < 1 &
                    & .or. maxval(info_spat%soil_use_id%mat) > pars%sim%n_lus) then
                    print *,"Soil use maps have soil uses not defined in crop database"
                    print *,"Verify the maximum allowed crop uses (SoilUsesNum) and soil maps"
                    print *, 'Execution will be aborted...'
                    stop!
                end if
                !
                do w=1,size(pars%sim%no_lu_list)
                    where (info_spat%soil_use_id%mat == pars%sim%no_lu_list(w)) info_spat%soil_use_id%mat = info_spat%soil_use_id%header%nan
                end do
                
                ! Soil uses that aren't simulated are removed from domain
                call overlay_domain (info_spat%soil_use_id,info_spat%domain) 
                
                ! update irrigation related parameters map as they can change during the simulation period
                ! TODO: add check if they exist
                ! TODO: make a function for that
                select case(pars%sim%mode)
                    case(1)
                        call read_grid (trim(pars%sim%input_path)// &
                            & trim(pars%sim%id_irr_meth_fn)//'_'//trim(adjustl(s_years))//".asc", &
                            & info_spat%irr_meth_id,pars%sim,boundaries)
                        call set_default_par (info_spat%irr_meth_id, info_spat%domain,1)
                        info_spat%h_meth=info_spat%domain
                        info_spat%h_meth%mat=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%h_irr) ! Spreads h_irr for each irrigation method
                        call init_irrigation_units(info_spat%domain,info_spat%irr_unit_id,info_spat%eff_net,irr_units,wat_src_tbl,&
                            &pars,info_spat%h_meth, debug)!
                            ! Spatial distribution of irrigation districts' (soil water content thresholds for irrigation application)
                        alpha_ms_map=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%irr_th_ms)    ! Spreads activation threshold for each irrigation method
                        alpha_unm_map=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%irr_th_unm)  ! Spreads activation threshold for each irrigation method
                        fw_irr=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%f_wet)             ! Spreads fw for each irrigation method
                        a_loss=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%a_loss)            ! Spreads irrigation application loss for each irrigation method
                        b_loss=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%b_loss)
                        c_loss=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%c_loss)
                        f_interception=id_to_par(info_spat%irr_meth_id,dble(pars%irr%met(:)%f_interception))
                        ! %EAC%: update percolation parameters
                        call calc_perc_booster_pars(info_spat,pars%irr%met,pars%sim%quantiles)
                        
                        if (debug .eqv. .true.) then
                            call write_grid(trim(pars%sim%path)//&
                                & 'out_'//trim(pars%sim%soiluse_fn)//'_'//trim(adjustl(s_years))//'.asc', &
                                & info_spat%soil_use_id, error_flag)
                            call write_grid(trim(pars%sim%path)//&
                                & 'out'//trim(pars%sim%id_irr_meth_fn)//'_'//trim(adjustl(s_years))//'.asc', &
                                & info_spat%irr_meth_id, error_flag)
                        end if
                    case(2)
                        call read_grid (trim(pars%sim%input_path)// &
                            & trim(pars%sim%id_irr_meth_fn)//'_'//trim(adjustl(s_years))//".asc", &
                            & info_spat%irr_meth_id,pars%sim,boundaries)
                        call set_default_par (info_spat%irr_meth_id,info_spat%domain,1)
                        call read_grid (trim(pars%sim%input_path)// &
                            & trim(pars%sim%eff_irr_fn)//'_'//trim(adjustl(s_years))//".asc",&
                            & info_spat%eff_met,pars%sim,boundaries)
                        call set_default_par (info_spat%eff_met,info_spat%domain,1.0D0)
                        alpha_ms_map=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%irr_th_ms)
                        alpha_unm_map=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%irr_th_unm)
                        fw_irr=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%f_wet)
                        f_interception=id_to_par(info_spat%irr_meth_id,dble(pars%irr%met(:)%f_interception))
                        ! %EAC%: update percolation parameters
                        call calc_perc_booster_pars(info_spat,pars%irr%met,pars%sim%quantiles)
                        
                        if (debug .eqv. .true.) then
                            call write_grid(trim(pars%sim%path)//&
                                & 'out_'//trim(pars%sim%soiluse_fn)//'_'//trim(adjustl(s_years))//'.asc', &
                                & info_spat%soil_use_id, error_flag)
                            call write_grid(trim(pars%sim%path)//&
                                & 'out'//trim(pars%sim%id_irr_meth_fn)//'_'//trim(adjustl(s_years))//'.asc', &
                                & info_spat%irr_meth_id, error_flag)
                            call write_grid(trim(pars%sim%path)//&
                                & 'out'//trim(pars%sim%eff_irr_fn)//'_'//trim(adjustl(s_years))//'.asc', &
                                & info_spat%eff_met, error_flag)
                        end if
                    case(3)
                        call read_grid (trim(pars%sim%input_path)// &
                            & trim(pars%sim%id_irr_meth_fn)//'_'//trim(adjustl(s_years))//".asc", &
                            & info_spat%irr_meth_id,pars%sim,boundaries)
                        call set_default_par (info_spat%irr_meth_id,info_spat%domain,1)
                        info_spat%h_meth=info_spat%domain
                        info_spat%h_meth%mat=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%h_irr)
                        alpha_ms_map=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%irr_th_ms)
                        alpha_unm_map=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%irr_th_unm)
                        fw_irr=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%f_wet)
                        a_loss=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%a_loss)
                        b_loss=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%b_loss)
                        c_loss=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%c_loss)
                        f_interception=id_to_par(info_spat%irr_meth_id,dble(pars%irr%met(:)%f_interception))
                        ! EAC: update percolation parameters
                        call calc_perc_booster_pars(info_spat,pars%irr%met,pars%sim%quantiles)
                        
                        if (debug .eqv. .true.) then
                            call write_grid(trim(pars%sim%path)//&
                                & 'out_'//trim(pars%sim%soiluse_fn)//'_'//trim(adjustl(s_years))//'.asc', &
                                & info_spat%soil_use_id, error_flag)
                            call write_grid(trim(pars%sim%path)//&
                                & 'out'//trim(pars%sim%id_irr_meth_fn)//'_'//trim(adjustl(s_years))//'.asc', &
                                & info_spat%irr_meth_id, error_flag)
                        end if
                    case(4)
                        call read_grid (trim(pars%sim%input_path)// &
                            & trim(pars%sim%id_irr_meth_fn)//'_'//trim(adjustl(s_years))//".asc", &
                            & info_spat%irr_meth_id,pars%sim,boundaries)
                        call set_default_par (info_spat%irr_meth_id,info_spat%domain,1)
                        call read_grid (trim(pars%sim%input_path)// &
                            & trim(pars%sim%eff_irr_fn)//'_'//trim(adjustl(s_years))//".asc",&
                            & info_spat%eff_met,pars%sim,boundaries)
                        call set_default_par (info_spat%eff_met,info_spat%domain,1.0D0)
                        info_spat%h_meth=info_spat%domain
                        info_spat%h_meth%mat=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%h_irr)
                        alpha_ms_map=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%irr_th_ms)
                        alpha_unm_map=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%irr_th_unm)
                        fw_irr=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%f_wet)
                        a_loss=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%a_loss)
                        b_loss=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%b_loss)
                        c_loss=id_to_par(info_spat%irr_meth_id,pars%irr%met(:)%c_loss)
                        f_interception=id_to_par(info_spat%irr_meth_id,dble(pars%irr%met(:)%f_interception))  
                        ! %EAC%: update percolation parameters
                        call calc_perc_booster_pars(info_spat,pars%irr%met,pars%sim%quantiles)
                        
                        if (debug .eqv. .true.) then
                            call write_grid(trim(pars%sim%path)//&
                                & 'out_'//trim(pars%sim%soiluse_fn)//'_'//trim(adjustl(s_years))//'.asc', &
                                & info_spat%soil_use_id, error_flag)
                            call write_grid(trim(pars%sim%path)//&
                                & 'out'//trim(pars%sim%id_irr_meth_fn)//'_'//trim(adjustl(s_years))//'.asc', &
                                & info_spat%irr_meth_id, error_flag)
                            call write_grid(trim(pars%sim%path)//&
                                & 'out'//trim(pars%sim%eff_irr_fn)//'_'//trim(adjustl(s_years))//'.asc', &
                                & info_spat%eff_met, error_flag)
                        end if
                    case default
                end select
            end if ! end change soil use map condition

            ! Read all phenological tables and allocation of info_pheno%prm%tab(:,:)!
            call read_all_crop_pars(pars%sim%year_step(y),pars%sim%n_lus,info_pheno)!
            if (debug .eqv. .true.) then
                call check_pheno_parameters(info_pheno,info_meteo)!
                call seek_un(error_flag,unit_crop)!
                call init_cell_output_file(unit_crop,trim(pars%sim%path)//'Kcb_levels.csv',&
                    &'MeteoStat; SoilUse; nCrop; low; mid; high')
                do i=1, size(info_pheno)
                    do j = 1, size(info_pheno(i)%k_cb%tab, 2)
                        do z = 1, info_pheno(i)%n_crops_by_year(j)
                            write(unit_crop,*) trim(info_meteo(i)%filename(1:(index(trim(info_meteo(i)%filename),"."))-1)), &
                                & '; ', j, '; ', z, '; ', &
                                & info_pheno(i)%kcb_phases%low(j,z), '; ', &
                                & info_pheno(i)%kcb_phases%mid(j,z), '; ', info_pheno(i)%kcb_phases%high(j,z)
                        end do
                    end do
                end do
                close(unit_crop)!
                call init_cell_output_file(unit_crop,trim(pars%sim%path)//trim(adjustl(s_years))//'_PhenoLengths.csv',&
                    &'MeteoStat; SoilUse; nCrop; ii0; iie; iid')
                do i=1,size(info_pheno)!
                    do j=1,size(info_pheno(i)%ii0,1)
                        do z=1,size(info_pheno(i)%ii0,2)
                            if (info_pheno(i)%ii0(j,z)>0) then
                                write(unit_crop,*)trim(info_meteo(i)%filename(1:(index(trim(info_meteo(i)%filename),"."))-1)), &
                                    & '; ',j,'; ', z, '; ',info_pheno(i)%ii0(j,z),'; ',&
                                    & info_pheno(i)%iie(j,z),'; ',info_pheno(i)%iid(j,z)
                            end if
                        end do
                    end do
                end do!
                close(unit_crop)!
            end if

            ! Randomization and spatial distribution of crop emergence
            call get_uniform_sample(irandom,pars%sim%sowing_range,pars%sim%rand_symmetry)!
            call allocate_crop_map (crop_map,info_spat%domain%mat,pars%sim%n_crops,info_spat%domain%header%nan)
            ! make_random_emergence calculates reference data to estimate crop emergence date
            ! which will be calculated in populate_crop_pars_matrices
            call make_random_emergence(info_pheno,meteo_weight,dir_meteo,info_spat%domain,info_spat%soil_use_id%mat, &
                & crop_map, irandom, pars%sim%year_step(y))!
            
            if (debug .eqv. .true.) then
                ! write debug files of reference data for crop randomization
                call print_mat_as_grid(trim(pars%sim%path)//trim(adjustl(s_years))//"_irandom.asc", &
                    & info_spat%domain%header,irandom,error_flag)!
                do i=1,size(crop_map%ii0,3)
                    write(str,*)i
                    call print_mat_as_grid(trim(pars%sim%path)//trim(adjustl(s_years))//"_ii0_" & 
                        & //trim(adjustl(str))//".asc",info_spat%domain%header,crop_map%ii0(:,:,i), &
                        & error_flag)!
                    call print_mat_as_grid(trim(pars%sim%path)//trim(adjustl(s_years))//"_iie_" & 
                        & //trim(adjustl(str))//".asc",info_spat%domain%header,crop_map%iie(:,:,i), &
                        & error_flag)!
                    call print_mat_as_grid(trim(pars%sim%path)//trim(adjustl(s_years))//"_dij_" &
                        & //trim(adjustl(str))//".asc",info_spat%domain%header,crop_map%dij(:,:,i), &
                        & error_flag)!
                end do
            end if
            ! Writing crop parameters in cult matrix
            call populate_crop_yield_matrices(info_pheno,dir_phenofases,info_spat%domain,info_spat%soil_use_id%mat,crop_map,y)
            
            ! Inizialization of kcb_low, cult_switch and phase_switch
            pheno%k_cb_low = info_spat%domain%header%nan
            pheno%n_crop_in_year = 1
            pheno%pheno_idx = 1

            ! Allocation of yield variables
            call init_yearly_yield_output(yield, info_spat%domain%mat, size(info_pheno(1)%ii0,2))
            
            select case(pars%sim%mode)
                case (1)                                                ! USE mode
                    ! Read water sources and dynamic allocation of info_sources%deriv%qt(:,:)!
                    call read_water_sources(pars%sim%year_step(y),pars,info_sources)!
                    call nom_water_supply(trim(pars%sim%watsour_path)//trim(pars%sim%watsources_fn), &
                        & irr_units, info_sources, wat_src_tbl, pars%sim%f_shapearea, info_spat%domain%header%cellsize, &
                        & info_spat%cell_area%mat, info_spat%irr_unit_id%mat, debug)
                case (4)                                                ! CALENDAR mode
                    ! Calculating water supply on the basis of irrigation application calendar
                    call open_scheduled_irrigation(trim(pars%sim%watsour_path)//trim(pars%sim%sched_irr_fn), irr_sch, debug)!
                case default
            end select
            
            ! Weekly simulation settings
            if (pars%sim%step_out == 1) then
                if (y == 1) then
                    info_meteo(1)%start%weekday = &
                        & day_of_week(info_meteo(1)%start%day, info_meteo(1)%start%month, info_meteo(1)%start%year)
                    pars%sim%clock(1) = pars%sim%weekday - info_meteo(1)%start%weekday + 3
                else
                    pars%sim%clock(1) = 8 - (mod((pars%sim%clock(2)-pars%sim%clock(1))/pars%sim%clock(3),7) - pars%sim%clock(1)) ! 8 = 7 + 1 To switch to the next day
                end if
                if (pars%sim%clock(1) > 7) pars%sim%clock(1) = pars%sim%clock(1) - 7
                pars%sim%clock(2) = pars%sim%year_step(y)
                pars%sim%clock(3) = 7
                pars%sim%intervals = pars%sim%clock(3)
                pars%sim%intervals(1) = pars%sim%clock(1) - 1 ! Adjusting of first memorized interval
                pars%sim%intervals(54) = mod((pars%sim%clock(2)-pars%sim%clock(1))/pars%sim%clock(3),7)-pars%sim%clock(1)
                if (pars%sim%intervals(54) < 0) then ! Skipping last interval if needed
                    pars%sim%intervals(53) = pars%sim%intervals(53) + pars%sim%intervals(54)
                    pars%sim%intervals(54) = 0
                end if
                if (pars%sim%intervals(1) == 0) then ! Skipping first interval if void
                    pars%sim%intervals = cshift(pars%sim%intervals,1)
                end if
            end if
            
            ! Output files *.csv inizialization 
            if (pars%sim%mode == 1) then
                call init_cell_output_by_year(out_tbl_list, pars%sim%path, s_years, info_meteo%filename, &
                    & pars%sim%mode, pars%sim%f_out_cells, debug, &
                    & irr_units%id, pars%cr%n_withdrawals, info_sources%unm_src_tbl%wat_src_id)!
            else 
                call init_cell_output_by_year(out_tbl_list,pars%sim%path,s_years,info_meteo%filename, &
                    & pars%sim%mode, pars%sim%f_out_cells, debug)!
            end if
            if(pars%sim%f_out_cells .eqv. .true.)then!
                call write_cell_info(info_spat, out_tbl_list%cell_info, pars%sim%mode, pars%sim%f_cap_rise)
            end if!
            if (out_yearly .eqv. .true.) then
                call init_yearly_output_file(yr_map,pars%sim%path,s_years)
                ! TODO: Verify when output_yield_iniz needs to be activated
                call init_yield_output_file(yield,pars%sim%path,s_years)
                if (debug .eqv. .true.) call init_debug_yearly_output_file(yr_deb_map,pars%sim%path,s_years)
            end if

            ! Calendar initializing to take into account start of year
            if (info_meteo(1)%start%month > 2) then
                call days_x_month(days_in_yr, pars%sim%start_year+y)
            else
                call days_x_month(days_in_yr,pars%sim%start_year+y-1)!
            end if
            shift_days = calc_doy(info_meteo(1)%start%day, info_meteo(1)%start%month, pars%sim%start_year+y-1) - calc_doy(1, 1, pars%sim%start_year+y-1)
            days_in_yr = cshift(days_in_yr, info_meteo(1)%start%month-1)
            
            ! Daily simulation cycle
            day_cycle: do doy=1, min(pars%sim%year_step(y), pars%sim%year_step(y) - &
                    & (pars%sim%start_simulation%doy - (info_meteo(1)%start%doy + sum(pars%sim%year_step(1:(y-1))))))
                n_day=n_day+1        ! Counter updating
                coll_irr=0           ! Irrigation matrix for collective water sources
                priv_irr=0           ! Irrigation matrix for private water sources
                h_irr_sum = 0   ! Irrigation matrix
                h_bypass = 0         ! Field irrigation losses matrix
                iter1 = 0
                iter2 = 0
                h_irr = 0.
                write (s_doy,*) doy

                print*,'Simulation day', achar(9), trim(adjustl(s_doy)), achar(9), 'year', achar(9), trim(adjustl(s_year))

              
                ! Updating daily data matrix for water table depth
                if(pars%sim%f_cap_rise .eqv. .true.)then!
                    upfilename = trim(pars%sim%input_path)//trim(pars%sim%wat_table_fn)//"_"//trim(adjustl(s_year))//"_"&
                        & //trim(adjustl(s_doy))//".asc" ! "
                    inquire(file=trim(upfilename), exist=file_exists)   ! file_exists will be TRUE if the file exists
                    if (file_exists .eqv. .true.) then
                        print *,'Water table data are updated: ', upfilename
                        call read_grid(trim(upfilename), info_spat%wat_tab,pars%sim,boundaries)
                        ! fix water table depth to ze
                        where((info_spat%wat_tab%mat<pars%depth%ze_fix) .and. &
                            (info_spat%wat_tab%mat/=info_spat%wat_tab%header%nan))
                            info_spat%wat_tab%mat=pars%depth%ze_fix
                        end where
                    end if
                end if

                ! TODO Variabile inizialization only if in the first day of the interval
                ! TODO As in output_asc_month_iniz
                                
                ! Monthly output *.asc file inizialization
                if(pars%sim%step_out == 0)then!
                    call init_step_output_file(stp_map,pars%sim%path,s_years,doy,days_in_yr,0, 'month')!
                    if (debug .eqv. .true.) then
                        call init_step_debug_output_file(deb_map, pars%sim%path, s_years, doy, days_in_yr, 0, 'month')
                    end if
                else if (pars%sim%step_out == 1) then!
                    ! In output_asc_month_iniz, uses 0 as first day (as in monthly routine)
                    call init_step_output_file(stp_map,pars%sim%path,s_years,doy,pars%sim%intervals,0, 'week')
                    if (debug .eqv. .true.) then
                        call init_step_debug_output_file(deb_map, pars%sim%path, s_years, doy, pars%sim%intervals, 0, 'week')
                    end if
                else!
                    ! In output_asc_month_iniz, uses (StartDate - 1) as first day
                    call init_step_output_file(stp_map,pars%sim%path,s_years,doy,pars%sim%intervals,pars%sim%clock(1)-1, 'step')
                    if (debug .eqv. .true.) then
                        call init_step_debug_output_file(deb_map, pars%sim%path, s_years, doy, pars%sim%intervals, &
                        &   pars%sim%clock(1)-1, 'step')
                    end if
                end if!
                ! Phenological parameters spatialization
                ! Updating of pheno%kcb_old to the last day value - pheno%cult_switch is not updated
                pheno%k_cb_old = pheno%k_cb
                if (y==pars%sim%start_simulation%year - info_meteo(1)%start%year + 1 &
                    & .and. (pars%sim%start_simulation%day > 1 .or. pars%sim%start_simulation%month > 1)) then
                    call populate_crop_pars_matrices(pheno,info_pheno, irandom, &
                        & doy + pars%sim%start_simulation%doy - calc_doy(1, 1, pars%sim%start_simulation%year), &
                        & dir_phenofases,info_spat%domain,info_spat%soil_use_id, &
                        & y, pars%sim%year_step(y), crop_map)!
                else
                    call populate_crop_pars_matrices(pheno,info_pheno,irandom,doy,dir_phenofases,info_spat%domain,info_spat%soil_use_id, &
                    & y, pars%sim%year_step(y), crop_map)!
                end if
                
                ! Output fc in debug %RR%
                !if (debug .eqv. .true.) then
                !    pheno_grd%mat = pheno%fc
                !    call write_matrices(trim(xml%sim%path)//'fc_'//trim(adjustl(s_years))//'_'//&
                !                                    trim(adjustl(s_gg))//'.asc', pheno_grd, errorflag)
                !end if
                !!
                    
                ! Creating cell parameters output on first day of simulation
                if (doy == 1 .and. (pars%sim%f_out_cells .eqv. .true.))then!
                    call write_cell_prod(out_tbl_list%prod_info, crop_map, irandom)
                end if
                
                ! Inizialization on first day of simulation
                first_day: if(doy==1 .and. y==pars%sim%start_simulation%year - info_meteo(1)%start%year + 1) then!
                    day_from_irr=-9999!
                    esp_perc=1.    !
                    wat_bal1 = 0.0D0; wat_bal1_old = 0.0D0!
                    wat_bal2 = 0.0D0; wat_bal2_old = 0.0D0!
                    wat_bal1%d_e = pars%depth%ze_fix!

                    where(info_spat%domain%mat /= info_spat%domain%header%nan)!
                        ! Layer depths inizialization - to calculate h_soil and t_soil
                        ! Comparing root zone to evaporative layer depth
                        where(pheno%d_r > pars%depth%ze_fix)
                            wat_bal2%d_t =  pheno%d_r - pars%depth%ze_fix
                        else where
                            wat_bal2%d_t = pars%depth%zr_fix  ! If Sr < Ze, Zr_fix = 1 - Ze
                        end where!
                        ! Soil water content inizialization [mm]
                        wat_bal1%h_soil = info_spat%theta(1)%old%mat*1000.*pars%depth%ze_fix!
                        wat_bal2%h_soil = info_spat%theta(2)%old%mat*1000.*wat_bal2_old%d_t        !
                        ! Soil water content inizialization [m3/m3]
                        wat_bal1%t_soil = info_spat%theta(1)%old%mat!
                        wat_bal2%t_soil = info_spat%theta(2)%old%mat!
                    end where!
                end if first_day!

                ! Saving bil* values of (gg-1)-th iteration on bil*_old variables
                wat_bal1_old = wat_bal1!
                wat_bal2_old = wat_bal2!
                
                ! %EAC%: reset h_perc each day
                wat_bal1%h_perc = 0.0D0
                wat_bal2%h_perc = 0.0D0
                
                ! %EAC%: limit root depth to water table interface
                where ((info_spat%domain%mat /= info_spat%domain%header%nan) &
                            .and. (info_spat%wat_tab%mat<pheno%d_r))
                    pheno%d_r = info_spat%wat_tab%mat
                END where

                where(info_spat%domain%mat /= info_spat%domain%header%nan)!
                    ! Layer depths update as a function of d_r (phenological parameter - root depth)
                    where(pheno%d_r > pars%depth%ze_fix)!
                        wat_bal2%d_t = pheno%d_r - wat_bal1%d_e
                    else where!
                        wat_bal2%d_t = pars%depth%zr_fix!
                    end where!
                    ! Water table depth inizialization under root zone
                    wat_bal2%depth_under_rz = info_spat%wat_tab%mat - pheno%d_r!
                    ! Soil water content update
                    where(wat_bal2%d_t == wat_bal2_old%d_t)!
                        wat_bal2_old%h_soil = wat_bal2_old%t_soil*1000*wat_bal2_old%d_t!
                    else where(wat_bal2%d_t > wat_bal2_old%d_t)!
                        ! Paddy field correction (only during growth season)
                        where(pheno%cn_class==7 .and. pheno%k_cb>0)
                            wat_bal2_old%h_soil = wat_bal2_old%t_soil*1000*wat_bal2_old%d_t + &
                                                    theta2_rice%theta2_FC*1000*(wat_bal2%d_t-wat_bal2_old%d_t)!
                        else where
                            wat_bal2_old%h_soil = wat_bal2_old%t_soil*1000*wat_bal2_old%d_t + &
                                                    info_spat%theta(2)%fc%mat*1000*(wat_bal2%d_t-wat_bal2_old%d_t)!
                        end where
                    else where!
                        wat_bal2_old%h_soil = wat_bal2_old%t_soil*1000*wat_bal2_old%d_t - &
                                                wat_bal2_old%t_soil*1000*(wat_bal2_old%d_t-wat_bal2%d_t)!
                    end where!

                    pheno%p_day = pheno%p + 0.04*(5.-(wat_bal1_old%h_eva_pot + wat_bal1_old%h_transp_pot + wat_bal2_old%h_transp_pot))!
                    ! pheno%pday amendment if pday values are not in its range [0.1 ; 0.8]
                    where (pheno%p_day <0.1) pheno%p_day=0.1
                    where (pheno%p_day >0.8) pheno%p_day=0.8 ! TODO: larger values will be permitted in order to consider stress irrigation
                end where!
                ! 
                call calculate_RF_t(wat_bal2%d_t, pars%depth%ze_fix, pheno, info_spat%domain)
                !!
                ! Soil water thresholds update (wat variable)
                call x_wat(info_spat%domain,info_spat%theta,wat_bal1%d_e,wat_bal2%d_t,wat,theta2_rice,pheno%cn_class,pheno%k_cb)!
                ! Irrigation application thresholds update
                where(info_spat%domain%mat /= info_spat%domain%header%nan)!
!~                     bil2%RAWbig =   wat%layer(2)%fc - (wat%layer(2)%fc-wat%layer(2)%wp)*pheno%pday*mkraw
!~                     bil2%RAW =      wat%layer(2)%fc - (wat%layer(2)%fc-wat%layer(2)%wp)*pheno%pday
!~                     bil2%RAWinf =   wat%layer(2)%fc - (wat%layer(2)%fc-wat%layer(2)%wp)*((pheno%pday+1)/2)
!~                     bil2%RAWaz =    wat%layer(2)%fc - (wat%layer(2)%fc-wat%layer(2)%wp)*pheno%pday*mkpraw
                    ! wat_bal2%h_raw_sup =   wat%layer(1)%h_fc + wat%layer(2)%h_fc - &
                    !     & (wat%layer(1)%h_fc - wat%layer(1)%h_wp + wat%layer(2)%h_fc - wat%layer(2)%h_wp)*pheno%p_day*(alpha_ms_map+pheno%r_stress)

                    ! wat_bal2%h_raw    =   wat%layer(1)%h_fc + wat%layer(2)%h_fc - &
                    !     & (wat%layer(1)%h_fc - wat%layer(1)%h_wp + wat%layer(2)%h_fc - wat%layer(2)%h_wp)*pheno%p_day

                    ! wat_bal2%h_raw_inf    =   wat%layer(1)%h_fc + wat%layer(2)%h_fc - &
                    !     & (wat%layer(1)%h_fc - wat%layer(1)%h_wp + wat%layer(2)%h_fc - wat%layer(2)%h_wp)*((pheno%p_day+1)/2)
                        
                    ! wat_bal2%h_raw_priv    =   wat%layer(1)%h_fc + wat%layer(2)%h_fc - &
                    !     & (wat%layer(1)%h_fc - wat%layer(1)%h_wp + wat%layer(2)%h_fc - wat%layer(2)%h_wp)*pheno%p_day*(alpha_unm_map+pheno%r_stress)
                
                
                    wat_bal2%h_raw_sup =  (wat%layer(1)%h_fc - (wat%layer(1)%h_fc-wat%layer(1)%h_wp)*pheno%p_day*(alpha_ms_map+pheno%r_stress))*pheno%RF_e + &
                                          (wat%layer(2)%h_fc - (wat%layer(2)%h_fc-wat%layer(2)%h_wp)*pheno%p_day*(alpha_ms_map+pheno%r_stress))*pheno%RF_t 
                    
                    wat_bal2%h_raw     =  (wat%layer(1)%h_fc - (wat%layer(1)%h_fc-wat%layer(1)%h_wp)*pheno%p_day)*pheno%RF_e + &
                                          (wat%layer(2)%h_fc - (wat%layer(2)%h_fc-wat%layer(2)%h_wp)*pheno%p_day)*pheno%RF_t 

                    wat_bal2%h_raw_inf =  (wat%layer(1)%h_fc - (wat%layer(1)%h_fc-wat%layer(1)%h_wp)*((pheno%p_day+1)/2))*pheno%RF_e + &
                                          (wat%layer(2)%h_fc - (wat%layer(2)%h_fc-wat%layer(2)%h_wp)*((pheno%p_day+1)/2))*pheno%RF_t 
                        
                    wat_bal2%h_raw_priv    =   (wat%layer(1)%h_fc - (wat%layer(1)%h_fc-wat%layer(1)%h_wp)*pheno%p_day*(alpha_unm_map+pheno%r_stress))*pheno%RF_e + &
                                               (wat%layer(2)%h_fc - (wat%layer(2)%h_fc-wat%layer(2)%h_wp)*pheno%p_day*(alpha_unm_map+pheno%r_stress))*pheno%RF_t 

                
                end where
                !!
                ! read weather daily data and calculate ET0 for each weather stations 
                call read_meteo_data(info_meteo,doy,pars%sim%res_canopy(y), pars%sim%forecast_day)!
                
                ! spread weather data to the entire domain
                if (y==pars%sim%start_simulation%year - info_meteo(1)%start%year + 1 &
                    & .and. (pars%sim%start_simulation%day > 1 .or. pars%sim%start_simulation%month > 1)) then
                    call create_meteo_matrices(info_meteo,dir_meteo,meteo_weight,meteo,info_spat%domain,&
                        & doy + pars%sim%start_simulation%doy - calc_doy(1, 1, pars%sim%start_simulation%year),&
                        & pars%sim%res_canopy(y))
                else
                    call create_meteo_matrices(info_meteo,dir_meteo,meteo_weight,meteo,info_spat%domain,doy,pars%sim%res_canopy(y))
                end if

                ! calculate average latitude
                if (doy == 1) then 
                    forall (i=1:size(meteo%lat,1), j=1:size(meteo%lat,2), meteo%lat(i,j)/=nan_r)!
                        lat_sum = lat_sum + meteo%lat(i,j)
                        lat_num = lat_num + 1
                    end forall!
                    lat_mean = lat_sum/lat_num
                end if

                ! calculate day length
                call calculateDLH(doy, lat_mean, DLH)
                ! calculate radiation distribution along day and update params
                pars%fet0 = pdf_normal(cost_hrs, 12.5D0, DLH/5)
                
                ! calculate the effective precipitation for the rice (only precipitation is considered)
                wat_bal1%h_interc=calc_interception(meteo%p,pheno)
                wat_bal1%h_eff_rain = net_precipitation(meteo%p,wat_bal1%h_interc)
                
                ! calculate the temperature stress factor
                ! TODO: move to day of the year
                jul_day = calc_doy(info_meteo(1)%start%day, info_meteo(1)%start%month, info_meteo(1)%start%year) + y + doy
                call calc_date(jul_day, day_cal, month_cal, year_cal)
                tmax_d = tmax_time(month_cal)
                tmin_d = tmin_time(month_cal)
                do time = 8, 14                  
                    meteo%T_ave = (meteo%T_max + meteo%T_min)/2 - &
                        & (meteo%T_max - meteo%T_min)/2 * cos(pi*(time - tmin_d)/(tmax_d - tmin_d)) + meteo%T_ave
                end do
                meteo%T_ave = meteo%T_ave / (14-8+1)

                ! TODO: implement separated subroutine for each simulation mode
                !if(doy>=pars%sim%start_irr_season .and. doy<=pars%sim%end_irr_season)then!
                
                ! define irrigation height base on irrigation period and specific condiction
                select case (pars%sim%mode)
                    case (1)                        ! use mode
                        ! %AB%: init the the cumulative value at the beginning of the season
                        !if (doy==pars%sim%start_irr_season) irr_units(:)%q_surplus = 0
                        ! %EAC%: as the irrigation season can change with the irrigation methods,
                        ! q_surplus is updated at the beginning of the year
                        ! TODO: manage condition when irrigation season is in winter
                        if (doy==1) irr_units(:)%q_surplus = 0
                        ! calculate the daily water duty for each irrigation units, considering the water distribution efficiency 
                        call calc_daily_duty(doy, irr_units, info_sources, wat_src_tbl, info_spat%irr_unit_id, &
                            & info_spat%domain, pars, pheno%irrigation_class, &
                            & (wat_bal1_old%h_soil + wat_bal2_old%h_soil), (wat_bal1_old%h_transp_pot + wat_bal2_old%h_transp_pot), &
                            & wat_bal2%h_raw, info_spat%theta(2)%fc%mat, wat_bal2_old%d_t)!
                        ! %EAC%: save irrigation units results
                        call save_irr_unit_debug_data(doy,out_tbl_list,irr_units)

                        call irrigation_use(info_spat%domain, info_spat%irr_unit_id, pheno%irrigation_class,&
                            & info_spat%irr_meth_id, irr_units, &
                            & (wat_bal1_old%h_transp_pot+wat_bal2_old%h_transp_pot), (wat_bal1_old%h_soil + wat_bal2_old%h_soil), &
                            & wat_bal2%h_raw_sup, wat_bal2%h_raw_inf, wat_bal2%h_raw, wat_bal2%h_raw_priv, &
                            & h_irr, doy, priv_irr, coll_irr, day_from_irr, esp_perc, &
                            & info_spat%a_perc, info_spat%b_perc, pars%sim%f_shapearea, info_spat%cell_area%mat, &
                            & info_spat%h_meth%mat, info_spat%irr_starts%mat, info_spat%irr_ends%mat)!pars%sim%end_irr_season)!
                        ! %EAC%: save irrigation units results
                        call save_irr_unit_data(doy,out_tbl_list,irr_units)

                        ! update irrigation losses
                        call calc_irrigation_losses(a_loss, b_loss, c_loss, meteo%Wind_vel, 0.5*(meteo%T_max+meteo%T_min),irr_loss)
                        ! calculate net irrigation
                        do z=1, pars%sim%n_irr_meth
                            h_irr(:,:,z) = h_irr(:,:,z) * (1.0-irr_loss/100.0)
                        end do
                        ! %AB% init the cumulative value
                        ! %EAC% not sure that q_surplus must be initialized to zero at the beginning and the end of the irrigation period
                        !if (doy==pars%sim%end_irr_season) irr_units(:)%q_surplus = 0
                    case (2)                        ! need mode at field capacity
                        call irrigation_need_fc(info_spat, h_irr, wat_bal2, wat_bal2_old, wat_bal1_old, pheno, &
                            & wat_bal1%h_eff_rain, theta2_rice%k_sat_2, day_from_irr, esp_perc)
                        ! if outside the irrigation period, set irrigation height to zero
                        do z=1, pars%sim%n_irr_meth
                            where(doy<info_spat%irr_starts%mat .or. doy>info_spat%irr_ends%mat) h_irr(:,:,z) = 0.
                        end do
                        irr_loss = 0. ! not consider irrigation losses
                    case (3)                        ! need mode with fixed volume
                        call irrigation_need_fix(info_spat, h_irr, wat_bal2, wat_bal2_old, wat_bal1_old, pheno, &
                            & wat_bal1%h_eff_rain, theta2_rice%k_sat_2, day_from_irr, esp_perc)
                        ! update irrigation losses
                        call calc_irrigation_losses(a_loss, b_loss, c_loss, meteo%Wind_vel, 0.5*(meteo%T_max+meteo%T_min),irr_loss)
                        ! if outside the irrigation period, set irrigation height to zero
                        ! and calculate net irrigation
                        do z=1, pars%sim%n_irr_meth
                            where(doy<info_spat%irr_starts%mat .or. doy>info_spat%irr_ends%mat) h_irr(:,:,z) = 0.
                            h_irr(:,:,z) = h_irr(:,:,z) * (1.0-irr_loss/100.0)
                        end do
                    case (4)                        ! scheduled irrigation mode
                        call irrigation_scheduled(info_spat, doy, current_year, irr_sch, pheno, &
                            & h_irr, day_from_irr, esp_perc, debug, wat_bal1_old, wat_bal2, wat_bal2_old, &
                            & a_loss, b_loss, c_loss, meteo%Wind_vel, 0.5*(meteo%T_max+meteo%T_min),irr_loss,&
                            wat_bal1%h_eff_rain, theta2_rice%k_sat_2)
                        ! if outside the irrigation period, set irrigation height to zero
                        ! and calculate net irrigation
                        do z=1, pars%sim%n_irr_meth
                            where(doy<info_spat%irr_starts%mat .or. doy>info_spat%irr_ends%mat) h_irr(:,:,z) = 0.
                            h_irr(:,:,z) = h_irr(:,:,z) * (1.0-irr_loss/100.0)
                        end do
                    case default
                end select
                
                h_irr_sum = sum(h_irr,dim=3)
                
                ! irrigation losses due to the irrigation method
                h_bypass = h_irr_sum * irr_loss / (100 - irr_loss)
                where (h_irr_sum/=0) yr_map%n_irr_events%mat = yr_map%n_irr_events%mat +1
                
                !end if ! end check irrigation season

                ! calculate intercetion according to the Von Hoyningen-Huene and Braden model 
                ! consider both precipitation and above canopy irrigation
                wat_bal1%h_interc = calc_interception(meteo%p+h_irr_sum * f_interception, pheno)
                wat_bal1%h_eff_rain = net_precipitation(meteo%p+h_irr_sum * f_interception, wat_bal1%h_interc)

                ! calculate the CN value
                call CN_table(tab_CN2, tab_CN3,info_spat%drainage,pheno%cn_class,pheno%cn_day, &!
                    & out_cn_day,info_spat%domain,info_spat%hydr_gr, &!
                    & info_spat%slope%mat, out_cn,info_spat%theta, &
                    & wat_bal1%t_soil,wat_bal2%t_soil)

                ! init water balance variables
                call init_water_balance_variables(wat_bal1,wat_bal2)
                
                ! set soil water content of current day to the previous day
                wat_bal_hour%inten%h_soil1 = wat_bal1_old%h_soil!
                wat_bal_hour%inten%h_soil2 = wat_bal2_old%h_soil!
                
                ! the gross available water is the sum of precipitation and above canopy irrigation
                wat_bal1%h_gross_av_water = meteo%p + h_irr_sum*f_interception

                ! the net available precipitation is the net precipitation + ponding
                wat_bal1%h_net_av_water = wat_bal1%h_eff_rain ! + wat_bal1_old%h_pond %CG% 2024-03-29 removed
                !wat_bal1%h_net_av_water = wat_bal1%h_eff_rain  + wat_bal1_old%h_pond

                
                ! calculate the runoff with the CN model
                call CN_runoff(wat_bal1%h_gross_av_water, wat_bal1%h_net_av_water, &
                    & h_irr_sum*(1-f_interception), info_spat%domain, pheno, &
                    & wat_bal1%h_runoff, out_cn_day, pars%sim%lambda_cn)
                
                ! HOURLY LOOP OF THE SIMULATION
                            
                
                hr_loop: do hour = 1,24
                    ! init the variables to zero (except for f_eff_rain, h_net_av_water)
                    wat_bal_hour = 0.0D0
                    ! calculate the precipitation and infiltration hourly distribution
                    where (info_spat%domain%mat /= info_spat%domain%header%nan)
                        wat_bal_hour%esten%h_eff_rain = wat_bal1%h_eff_rain*pars%f_eff_rain(hour)
                        wat_bal_hour%esten%h_inf  = wat_bal1%h_net_av_water*pars%f_eff_rain(hour)
                    end where
                    ! calculate parameters for evaporation that remain constant during the day
                    if(hour==1) then
                        call b1_no_iter_eva(pheno,meteo, pars%sim%h_prec_lim, wat, fw_day, fw_irr, pars%irr%f_w, fc, &
                                            h_irr_sum, f_interception, info_spat%domain, wat_bal1_old, fw_old)
                        wat_bal_hour%inten%h_pond0 = 0.
                        !%CG% 2024-03-29 add total ponding to the first hour
                        wat_bal_hour%esten%h_inf = wat_bal_hour%esten%h_inf + wat_bal1_old%h_pond
                    end if
                    
                    ! compute water balance for each cells
                    do j=1, size(info_spat%domain%mat,2)
                        do i=1,size(info_spat%domain%mat,1)
                            if(info_spat%domain%mat(i,j)/=info_spat%domain%header%nan)then!
                                ! %RR%: add k_r
                                ! water balance for the evaporative layer
                                call water_balance_evap_lay(h_irr(i,j, info_spat%irr_meth_id%mat(i,j)) * &
                                    & pars%irr%met(info_spat%irr_meth_id%mat(i,j))%freq(hour)  &
                                    & * (1-pars%irr%met(info_spat%irr_meth_id%mat(i,j))%f_interception), &
                                    & wat_bal_hour%inten%h_soil1(i,j), wat_bal_hour%esten%h_inf(i,j), &
                                    & wat_bal_hour%esten%h_eva(i,j), wat_bal_hour%esten%h_eva_pot(i,j), &
                                    & wat_bal_hour%esten%h_perc1(i,j), wat_bal_hour%inten%h_pond0(i,j), wat_bal_hour%esten%h_pond(i,j), &
                                    & wat_bal_hour%esten%h_transp_act1(i,j), wat_bal_hour%esten%h_transp_pot1(i,j), &
                                    & pheno%k_cb(i,j), pheno%p_day(i,j),  meteo%et0(i,j)*pars%fet0(hour), wat_bal_hour%esten%k_e(i,j), &
                                    & wat_bal_hour%esten%k_r(i,j), wat_bal_hour%inten%k_s_dry(i,j), wat_bal_hour%inten%k_s_sat(i,j),&
                                    & wat_bal_hour%inten%k_s(i,j), wat%kc_max(i,j), wat%few(i,j), pheno%RF_e(i,j), &
                                    & wat%wat1_rew(i,j), wat%layer(1)%h_sat(i,j), wat%layer(1)%h_fc(i,j), &
                                    & wat%layer(1)%h_wp(i,j), wat%layer(1)%h_r(i,j), &
                                    & info_spat%k_sat(1)%mat(i,j), info_spat%fact_n(1)%mat(i,j), &
                                    & wat_bal_hour%n_iter1(i,j), esp_perc(i,j,1), wat_bal_hour%n_max1(i,j))
                                    
                                ! water balance for the transpirative layer
                                call water_balance_transp_lay(wat_bal_hour%inten%h_soil2(i,j), wat_bal_hour%esten%h_transp_act2(i,j), &!
                                    & wat_bal_hour%esten%h_transp_pot2(i,j), wat_bal_hour%esten%h_perc2(i,j), &!
                                    & wat_bal_hour%esten%h_perc1(i,j), &
                                    & wat_bal_hour%inten%k_s_dry(i,j),wat_bal_hour%inten%k_s_sat(i,j), wat_bal_hour%inten%k_s(i,j), &!
                                    & wat_bal_hour%esten%h_eva_pot(i,j), wat_bal_hour%esten%h_caprise(i,j), &!
                                    & wat_bal_hour%esten%h_rise(i,j), &!
                                    & pheno%d_r(i,j), wat_bal2%d_t(i,j), pheno%RF_t(i,j), &!
                                    & pheno%k_cb(i,j), pheno%p_day(i,j), pheno%cn_class(i,j), &
                                    & meteo%et0(i,j)*pars%fet0(hour), wat%layer(2)%h_sat(i,j), &
                                    & wat%layer(2)%h_fc(i,j), wat%layer(2)%h_wp(i,j), wat%layer(2)%h_r(i,j), &!
                                    & info_spat%k_sat(2)%mat(i,j), info_spat%fact_n(2)%mat(i,j), &!
                                    & info_spat%a3%mat(i,j), info_spat%a4%mat(i,j), &!
                                    & info_spat%b1%mat(i,j), info_spat%b2%mat(i,j), &!
                                    & info_spat%b3%mat(i,j), info_spat%b4%mat(i,j), &!
                                    & wat_bal2%depth_under_rz(i,j), wat_bal_hour%n_iter2(i,j), &!
                                    & esp_perc(i,j,2),pars%sim%f_cap_rise, wat_bal_hour%n_max2(i,j))!
                            end if
                            
                        end do
                    end do
                    
                    ! TODO: %AB% move to subroutine
                    ! update the soil water content of the evaporative layer according to the rise from the transpirative layer
                    where (wat_bal_hour%esten%h_rise > 0)
                        wat_bal_hour%inten%h_soil1 = wat_bal_hour%inten%h_soil1 + wat_bal_hour%esten%h_rise
                        wat_bal_hour%esten%h_pond = merge (wat_bal_hour%esten%h_pond + (wat_bal_hour%inten%h_soil1 - wat%layer(1)%h_sat), &
                            & wat_bal_hour%esten%h_pond, wat_bal_hour%inten%h_soil1 > wat%layer(1)%h_sat)
                        wat_bal_hour%inten%h_soil1 = merge (wat%layer(1)%h_sat, wat_bal_hour%inten%h_soil1, &
                            & wat_bal_hour%inten%h_soil1 > wat%layer(1)%h_sat)
                    end where

                    ! wat_bal_hour%inten%h_soil1 = wat_bal_hour%inten%h_soil1 + wat_bal_hour%esten%h_rise
                    ! wat_bal_hour%esten%h_pond = merge (wat_bal_hour%inten%h_soil1 - wat%layer(1)%h_sat, &
                    !          & 0.0D0, wat_bal_hour%inten%h_soil1 > wat%layer(1)%h_sat)
                    ! wat_bal_hour%inten%h_soil1 = merge (wat%layer(1)%h_sat, wat_bal_hour%inten%h_soil1, &
                    !          & wat_bal_hour%inten%h_soil1 > wat%layer(1)%h_sat)

                    ! update water balance variables
                    wat_bal1%h_soil = wat_bal_hour%inten%h_soil1!
                    wat_bal1%h_inf = wat_bal1%h_inf + wat_bal_hour%esten%h_inf!
                    wat_bal1%h_eva = wat_bal1%h_eva + wat_bal_hour%esten%h_eva!
                    wat_bal1%h_eva_pot = wat_bal1%h_eva_pot + wat_bal_hour%esten%h_eva_pot!
                    wat_bal1%h_transp_act = wat_bal1%h_transp_act + wat_bal_hour%esten%h_transp_act1!
                    wat_bal1%h_transp_pot = wat_bal1%h_transp_pot + wat_bal_hour%esten%h_transp_pot1!
                    wat_bal1%h_perc = wat_bal1%h_perc + wat_bal_hour%esten%h_perc1!
                    wat_bal2%h_soil = wat_bal_hour%inten%h_soil2!
                    wat_bal2%h_transp_act = wat_bal2%h_transp_act + wat_bal_hour%esten%h_transp_act2!
                    wat_bal2%h_transp_pot = wat_bal2%h_transp_pot + wat_bal_hour%esten%h_transp_pot2!
                    wat_bal2%h_perc = wat_bal2%h_perc + wat_bal_hour%esten%h_perc2!
                    wat_bal2%k_s = wat_bal_hour%inten%k_s!
                    wat_bal2%h_caprise = wat_bal2%h_caprise + wat_bal_hour%esten%h_caprise!
                    wat_bal2%h_rise = wat_bal2%h_rise + wat_bal_hour%esten%h_rise!
                    wat_bal_hour%inten%h_pond0 = wat_bal_hour%esten%h_pond
                  

                    ! update the number of iterations
                    iter1 = merge(iter1,wat_bal_hour%n_iter1,iter1>wat_bal_hour%n_iter1)
                    iter2 = merge(iter2,wat_bal_hour%n_iter2,iter2>wat_bal_hour%n_iter2)
                    
                    if (debug .eqv. .true.) then
                        ! print the number of iteration for each control cells
                        if(pars%sim%f_out_cells .eqv. .true.)then!
                            do i=1,size(out_tbl_list%cell_conv)!
                                xx=out_tbl_list%cell_conv(i)%coord%row!
                                yy=out_tbl_list%cell_conv(i)%coord%col!
                                write(out_tbl_list%cell_conv(i)%file%unit,'(i3,a1,i2,a1,(2(i2,a1,i4,a1)))')doy,'; ',hour,'; ',&
                                    & wat_bal_hour%n_max1(xx,yy), '; ', wat_bal_hour%n_iter1(xx,yy), '; ',&
                                    & wat_bal_hour%n_max2(xx,yy), '; ',wat_bal_hour%n_iter2(xx,yy)!
                            end do!
                        end if
                    end if
                end do hr_loop
                
                ! update the soil water content (dimensionless)
                if(doy == pars%sim%year_step(y)) then
                    where(info_spat%domain%mat /= info_spat%domain%header%nan)!
                        info_spat%theta(1)%old%mat = wat_bal1%h_soil/(1000.*wat_bal1%d_e)!
                        info_spat%theta(2)%old%mat = wat_bal2%h_soil/(1000.*wat_bal2%d_t)!
                    end where!
                end if!

                where(info_spat%domain%mat /= info_spat%domain%header%nan)!
                        wat_bal1%t_soil = wat_bal1%h_soil/(1000.*wat_bal1%d_e)!
                        wat_bal2%t_soil = wat_bal2%h_soil/(1000.*wat_bal2%d_t)!
                end where!

                ! update the ponding variable for each day
                wat_bal1%h_runoff = wat_bal1%h_runoff+ max(wat_bal_hour%esten%h_pond-info_spat%h_maxpond%mat,0.0D0)
                wat_bal1%h_pond = min(wat_bal_hour%esten%h_pond,info_spat%h_maxpond%mat)
                !wat_bal1%h_pond = wat_bal_hour%esten%h_pond

                ! Calculate the crop production
                if(out_yearly .eqv. .true.)then
                    ! TODO:
                    ! add the crop biomass from the previuos year for winter cereals (need variables to store previous year)
                
                    ! Update the parameters for the calculation of the thermal stress
                    do j=1,size(info_spat%domain%mat,2)
                        do i=1,size(info_spat%domain%mat,1)
                            if(info_spat%domain%mat(i,j) /= info_spat%domain%header%nan) then
                                if (doy >= crop_map%TSP_low(i,j,pheno%n_crop_in_year(i,j)) .and. &
                                    & doy < crop_map%TSP_high(i,j,pheno%n_crop_in_year(i,j))) then
                                    if (meteo%T_ave(i,j) < pheno%T_crit(i,j)) then
                                        yield%f_HS_sum%mat(i,j,pheno%n_crop_in_year(i,j)) = &
                                            & yield%f_HS_sum%mat(i,j,pheno%n_crop_in_year(i,j)) + 1
                                    else if (meteo%T_ave(i,j) >= pheno%T_crit(i,j) .and. meteo%T_ave(i,j) < pheno%T_lim(i,j)) then
                                        yield%f_HS_sum%mat(i,j,pheno%n_crop_in_year(i,j)) = &
                                            & yield%f_HS_sum%mat (i,j,pheno%n_crop_in_year(i,j)) + 1 - &
                                            & (meteo%T_ave(i,j) - pheno%T_crit(i,j))/ (pheno%T_lim(i,j) - pheno%T_crit(i,j))
                                    end if
                                end if
                            
                                ! Calculate the period of growing
                                if (pheno%k_cb_low(i,j) == 0) then ! annual crop 
                                    if (pheno%k_cb(i,j) == pheno%k_cb_low(i,j)) then
                                        pheno%pheno_idx(i,j) = 0
                                    ! initial stage
                                    else if (pheno%k_cb(i,j) <= pheno%k_cb_mid(i,j) .and. &
                                        & (pheno%pheno_idx(i,j)==0 .or. pheno%pheno_idx(i,j)==1)) then
                                        pheno%pheno_idx(i,j) = 1
                                    ! growing stage
                                    else if (pheno%k_cb(i,j) < pheno%k_cb_high(i,j) .and. &
                                        & (pheno%pheno_idx(i,j)==1 .or. pheno%pheno_idx(i,j)==2)) then
                                        pheno%pheno_idx(i,j) = 2
                                    ! maturity stage
                                    else if (pheno%k_cb(i,j) == pheno%k_cb_high(i,j)) then
                                        pheno%pheno_idx(i,j) = 3
                                    ! senescence stage
                                    else
                                        pheno%pheno_idx(i,j) = 4
                                    end if
                                else  ! permanent, pluriannual crop
                                    ! Vernalization or after the harvest
                                    if (pheno%k_cb(i,j) == pheno%k_cb_low(i,j)) then
                                        pheno%pheno_idx(i,j) = 1
                                    ! growing stage
                                    else if (pheno%k_cb(i,j) < pheno%k_cb_high(i,j) .and. &
                                        & (pheno%pheno_idx(i,j)==1 .or. pheno%pheno_idx(i,j)==2)) then
                                        pheno%pheno_idx(i,j) = 2
                                    ! maturity stage
                                    else if (pheno%k_cb(i,j) == pheno%k_cb_high(i,j)) then
                                        pheno%pheno_idx(i,j) = 3
                                    ! senescence stage
                                    else
                                        pheno%pheno_idx(i,j) = 4
                                    end if
                                end if
                                
                                ! update the parameters for the calculation of the water stress
                                if (pheno%pheno_idx(i,j) > 0) then
                                    yield%T_act_sum%mat(i,j, pheno%pheno_idx(i,j), pheno%n_crop_in_year(i,j)) = &
                                        & yield%T_act_sum%mat(i,j, pheno%pheno_idx(i,j), pheno%n_crop_in_year(i,j)) + &
                                        & wat_bal1%h_transp_act(i,j) + wat_bal2%h_transp_act(i,j)
                                    yield%T_pot_sum%mat(i,j, pheno%pheno_idx(i,j), pheno%n_crop_in_year(i,j)) = &
                                        & yield%T_pot_sum%mat(i,j, pheno%pheno_idx(i,j), pheno%n_crop_in_year(i,j)) + &
                                        & wat_bal1%h_transp_pot(i,j) + wat_bal2%h_transp_pot(i,j)
                                    yield%dev_stage%mat(i,j, pheno%pheno_idx(i,j), pheno%n_crop_in_year(i,j)) = &
                                        & yield%dev_stage%mat(i,j,pheno%pheno_idx(i,j), pheno%n_crop_in_year(i,j)) + 1
                                end if
                                
                                ! Update transpiration ratio
                                if ( meteo%et0(i,j)>0) then
                                    yield%transp_ratio_sum%mat(i,j,pheno%n_crop_in_year(i,j)) = &
                                        &  yield%transp_ratio_sum%mat(i,j,pheno%n_crop_in_year(i,j)) + &
                                        & (wat_bal1%h_transp_pot(i,j) + wat_bal2%h_transp_pot(i,j)) / meteo%et0(i,j)
                                end if    
                            end if
                        end do
                    end do
                end if
                
                ! calculate transpiration deficit index
                if (trim(pars_TDx%mode)/="none") then
                    call sum_TD(wat_bal1%h_transp_act+wat_bal2%h_transp_act,wat_bal1%h_transp_pot+wat_bal2%h_transp_pot,pheno%k_cb,n_day,y,TD)
                    do cont_td=1,pars_TDx%temp%n_ind     ! update TD values
                        call calc_TDx(info_spat%domain, n_day, y, tot_deficit(:,:,cont_td), pheno%k_cb,pars%sim%year_step(y), sim_years, &
                            & pars_TDx%temp%x(cont_td), TD,unit_Dxi(cont_td)%dxi)
                    end do
                    ! TODO: check if the following is necessary in "report" mode
                    if(mod((doy - pars_TDx%temp%delay),pars_TDx%temp%td)==0)then
                        n_week = n_week +1
                        do cont_td=1,pars_TDx%temp%n_ind ! calculate the sum over the integration period
                            write(str_td,*)pars_TDx%temp%x(cont_td)
                            str_td="td"//trim(adjustl(str_td))
                            ! Create the temporary files that calculate TDx for each period
                            call update_TDx_DB(tot_deficit(:,:,cont_td), trim(str_td), n_week, info_spat%domain,pars%sim%path)
                        end do
                    end if
                end if

                call write_daily_output (doy, meteo, info_meteo, pheno, h_irr_sum, wat_bal1, wat_bal2, wat_bal2_old, &
                    & info_spat, pars, wat, wat_bal_hour, fw_day, fw_old, esp_perc, out_cn, out_cn_day, h_bypass, coll_irr, priv_irr, out_tbl_list, &
                    & pars%sim%mode,pars%sim%f_out_cells,debug) !! %RR% fw_old

                ! save output files by step
                if(out_asc .eqv. .true.)then!
                    if (pars%sim%step_out == 0) then
                        call write_outputs_by_step (doy, meteo, h_irr_sum, wat_bal1, wat_bal2, &
                            & info_spat, coll_irr, priv_irr, stp_map, deb_map, h_bypass, days_in_yr, 0, debug,summary)        ! uscite mensili
                    else if (pars%sim%step_out == 1) then
                        call write_outputs_by_step (doy, meteo, h_irr_sum, wat_bal1, wat_bal2, &
                            & info_spat, coll_irr, priv_irr, stp_map, deb_map, h_bypass, pars%sim%intervals, 0, debug,summary)    ! uscite settimanali
                    else
                        call write_outputs_by_step (doy, meteo, h_irr_sum, wat_bal1, wat_bal2, &
                            & info_spat, coll_irr, priv_irr, stp_map, deb_map, h_bypass, pars%sim%intervals, pars%sim%clock(1)-1, &
                            & debug,summary)    ! scheduled outputs
                    end if
                end if
                
                ! save to file the output bu year
                if (out_yearly .eqv. .true.) then
                    yr_map%rain%mat = yr_map%rain%mat + meteo%p
                    yr_map%runoff%mat = yr_map%runoff%mat + wat_bal1%h_runoff
                    yr_map%net_flux_gw%mat = yr_map%net_flux_gw%mat + wat_bal2%h_perc - wat_bal2%h_caprise
                    
                    ! VERY IMPORTANT EDIT
                    ! %EAC%: water release for irrigation should be already controlled by the presence of crop in field
                    ! otherwise is an error
                    yr_map%irr%mat = yr_map%irr%mat + h_irr_sum + h_bypass
                    yr_map%irr_loss%mat = yr_map%irr_loss%mat + h_bypass
                    
                    where(pheno%k_cb/=0) 
                        yr_map%rain_crop_season%mat = yr_map%rain_crop_season%mat + meteo%p
                        yr_map%eva_pot_crop_season%mat = yr_map%eva_pot_crop_season%mat + wat_bal1%h_eva_pot
                        yr_map%eva_act_crop_season%mat = yr_map%eva_act_crop_season%mat + wat_bal1%h_eva
                        yr_map%transp_act%mat = yr_map%transp_act%mat + wat_bal1%h_transp_act + wat_bal2%h_transp_act
                        yr_map%transp_pot%mat = yr_map%transp_pot%mat + wat_bal1%h_transp_pot + wat_bal2%h_transp_pot
                    end where
                    
                    if (debug .eqv. .true.) then
                        yr_deb_map%eva_act_tot%mat = yr_deb_map%eva_act_tot%mat + wat_bal1%h_eva
                        yr_deb_map%iter1%mat = merge(yr_deb_map%iter1%mat,dble(iter1),iter1<yr_deb_map%iter1%mat)!
                        yr_deb_map%iter2%mat = merge(yr_deb_map%iter2%mat,dble(iter2),iter2<yr_deb_map%iter2%mat)!
                    end if
                end if

            end do day_cycle!
            
            ! Calculate productivity
            do j=1,size(info_spat%domain%mat,2)
                do i=1,size(info_spat%domain%mat,1)
                    do z=1, size(crop_map%TSP_high,3)
                        if(info_spat%domain%mat(i,j) /= info_spat%domain%header%nan) then
                            
                            yield%f_HS%mat(i,j,z) = yield%f_HS_sum%mat(i,j,z) / &
                                & (crop_map%TSP_high(i,j,z) - crop_map%TSP_low(i,j,z))
                            
                            yield%biomass_pot%mat(i,j,z) = &
                                & crop_map%wp_adj(i,j,z) * yield%transp_ratio_sum%mat(i,j,z)
                            
                            yield%yield_pot%mat(i,j,z) = &
                                & yield%biomass_pot%mat(i,j,z) * crop_map%HI(i,j,z)
                            
                            ! Calculate the reduction of the production from the water stress
                            yield%f_WS%mat(i,j,z) = 1 - pheno%Ky_tot(i,j) * &
                                & (1- sum(yield%T_act_sum%mat(i,j,:,z)) / &
                                & sum(yield%T_pot_sum%mat(i,j,:,z)))

                            if (yield%f_WS%mat(i,j,z) < 0) yield%f_WS%mat(i,j,z) = 0    ! limit to zero
                            
                            yield%f_WS_stage%mat(i,j,z) = (1 - pheno%Ky_pheno(i,j,1) * &
                                & (1 - (yield%T_act_sum%mat(i,j,1,z) / &
                                & yield%T_pot_sum%mat(i,j,1,z)))) &
                                & ** (yield%dev_stage%mat(i,j,1,z) &
                                & / sum(yield%dev_stage%mat(i,j,:,z)))
                            
                            yield%f_WS_stage%mat(i,j,z) = (1 - pheno%Ky_pheno(i,j,2) * &
                                & (1 - (yield%T_act_sum%mat(i,j,2,z) / &
                                & yield%T_pot_sum%mat(i,j,2,z)))) &
                                & ** (yield%dev_stage%mat(i,j,2,z) &
                                & / sum(yield%dev_stage%mat(i,j,:,z))) * &
                                & yield%f_WS_stage%mat(i,j,z)
                            
                            yield%f_WS_stage%mat(i,j,z) = (1 - pheno%Ky_pheno(i,j,3) * &
                                & (1 - (yield%T_act_sum%mat(i,j,3,z) / &
                                & yield%T_pot_sum%mat(i,j,3,z)))) &
                                & ** (yield%dev_stage%mat(i,j,3,z) &
                                & / sum(yield%dev_stage%mat(i,j,:,z))) * &
                                & yield%f_WS_stage%mat(i,j,z)
                            
                            yield%f_WS_stage%mat(i,j,z) = (1 - pheno%Ky_pheno(i,j,4) * &
                                & (1 - (yield%T_act_sum%mat(i,j,4,z) / &
                                & yield%T_pot_sum%mat(i,j,4,z)))) &
                                & ** (yield%dev_stage%mat(i,j,4,z) &
                                & / sum(yield%dev_stage%mat(i,j,:,z))) * &
                                & yield%f_WS_stage%mat(i,j,z)
                            
                            yield%yield_act%mat(i,j,z) = &
                                & yield%yield_pot%mat(i,j,z) * &
                                & min(yield%f_WS%mat(i,j,z), &
                                & yield%f_WS_stage%mat(i,j,z)) * &
                                & yield%f_HS%mat(i,j,z)
                        end if
                    end do
                end do
            end do
            
            if (out_yearly .eqv. .true.) then
                ! Calculate the annual efficiency for the use of the water inputs (rain and irrigation)
                yr_map%total_eff%mat = (yr_map%eva_act_crop_season%mat + yr_map%transp_act%mat) &
                    & / (yr_map%rain_crop_season%mat + yr_map%irr%mat)
                
                yr_map%h_irr_mean%mat = yr_map%irr%mat /  yr_map%n_irr_events%mat
                
                where (yr_map%n_irr_events%mat == 0)
                    yr_map%h_irr_mean%mat = info_spat%domain%header%nan
                end where
                
                if (summary .eqv. .false.) then
                    call save_yearly_data(yr_map,info_spat%domain,debug)
                else
                    call save_annual_irrigation_data(yr_map,info_spat%domain)
                end if
                
                call save_yield_data(yield,info_spat%domain)
               
                if (debug .eqv. .true.) then
                    yr_deb_map%rain_eff%mat = (yr_map%eva_act_crop_season%mat + yr_map%transp_act%mat) / yr_map%rain_crop_season%mat
                    call save_annual_debug_data(yr_deb_map, info_spat%domain)
                    call save_yield_debug_data(yield, info_spat%domain)
                end if
            end if

            ! close the csv files for cell outputs
            call close_cell_output_by_year(out_tbl_list,pars%sim%mode,pars%sim%f_out_cells,debug)
            ! destroy annual variables
            call destroy_infofeno_tab(info_pheno)
            call destroy_crop(crop_map)
            if (pars%sim%mode ==1) call destroy_water_sources_duty(info_sources)
            call destroy_yield_output(yield)
        end do year_cycle!
        
        ! Save output for the following year
        if (pars%sim%f_theta_out .eqv. .true.) then
            info_spat%theta(1)%old%mat = wat_bal1%h_soil/(1000.*wat_bal1%d_e)!
            info_spat%theta(2)%old%mat = wat_bal2%h_soil/(1000.*wat_bal2%d_t)!
            call write_grid(trim(pars%sim%final_condition)//trim(pars%sim%thetaI_end_fn)//'.asc',info_spat%theta(1)%old,error_flag)
            call write_grid(trim(pars%sim%final_condition)//trim(pars%sim%thetaII_end_fn)//'.asc',info_spat%theta(2)%old,error_flag)
        end if
        
        ! calculate DTx statistics
        if(trim(pars_TDx%mode)=="analysis")then
            print*,"TDx index statistics are calculated"
            do cont_td=1,pars_TDx%temp%n_ind
                write(str_td,*)pars_TDx%temp%x(cont_td)
                str_td="td"//trim(adjustl(str_td))
                call save_TDx_statistics(unit_deficit,info_spat%domain,pars%sim%path,pars_TDx%n,trim(str_td))
            end do
            str_delete="del "//trim(pars%sim%path)//"*.tmp"
            call system(trim(str_delete))
        else if(trim(pars_TDx%mode)=="application")then
            print*, "TDx index is calculated"
            do cont_td=1,pars_TDx%temp%n_ind
                write(str_td,*)pars_TDx%temp%x(cont_td)
                str_td="td"//trim(adjustl(str_td))
                do n_week = 1,(maxval(pars%sim%year_step(:))/pars_TDx%temp%td)+1
                    call make_TDx_report(info_spat%domain,pars%sim%path,n_week,trim(str_td))
                    print *, "TDx index ", trim(adjustl(str_td)), " has been calculated for period ", n_week
                end do
            end do
            str_delete="del "//trim(pars%sim%path)//"*.tmp"
            call system(trim(str_delete))
        else
            print*,"TDx index has not been calculated"
        end if
        !
        call destroy_all(stp_map,yr_map,deb_map,yr_deb_map,wat_bal1,wat_bal1_old,wat_bal2,wat_bal2_old,wat_bal_hour,meteo,wat,pheno, &
            & pars%sim%imax,pars%sim%jmax, debug)
        !!
    end subroutine simulation_manager!

    subroutine write_daily_output (doy, meteo, info_meteo, pheno, h_irr_sum, wat_bal1, wat_bal2, wat_bal2_old, &
        & info_spat, pars, wat, wat_bal_hour, fw, fw_old, esp_perc, out_cn, out_cn_day, h_bypass, coll_irr, priv_irr, out_tbl, mode,cells,debug)
        ! write daily output for each control cells
        implicit none
        integer, intent(in):: doy
        type(meteo_mat), intent(in):: meteo
        type(meteo_info), dimension(:), intent(in):: info_meteo
        type(crop_pars_matrices), intent(in):: pheno
        real(dp), dimension(:,:), intent(in):: h_irr_sum
        type(balance1_matrices), intent(in):: wat_bal1
        type(balance2_matrices), intent(in):: wat_bal2, wat_bal2_old
        type(spatial_info), intent(in):: info_spat
        type(parameters), intent(in)::pars
        real(dp), dimension(:,:), intent(in)::fw, fw_old ! %RR% test
        real(dp), dimension(:,:,:), intent(in):: esp_perc
        type(output_cn),dimension(:,:),intent(in)::out_cn
        real(dp),dimension(:,:),intent(in)::out_cn_day
        real(dp), dimension(:,:), intent(in):: h_bypass
        real(dp), dimension(:,:), intent(in):: coll_irr, priv_irr
        type(output_table_list), intent(in):: out_tbl
        integer, intent(in)::mode
        logical, intent(in)::cells,debug
        
        type(wat_matrix)::wat
        type(hourly)::wat_bal_hour
        integer::i
        integer::xx,yy

        if (cells .eqv. .true.) then
            do i=1,size(out_tbl%sample_cells)!
                xx=out_tbl%sample_cells(i)%coord%row!
                yy=out_tbl%sample_cells(i)%coord%col!
                select case (mode)
                    case (0, 2:4)
                        write(out_tbl%sample_cells(i)%file%unit,*)doy, &!
                            & ';', meteo%p(xx,yy), ';', meteo%T_max(xx,yy), &!
                            & ';', meteo%T_min(xx,yy), ';', meteo%et0(xx,yy), &!
                            & ';', pheno%k_cb(xx,yy), ';',pheno%lai(xx,yy), ';',pheno%p_day(xx,yy), &!
                            & ';', h_irr_sum(xx,yy),&!
                            & ';', wat_bal1%h_eff_rain(xx,yy), ';', wat_bal1%h_gross_av_water(xx,yy), ';',wat_bal1%h_net_av_water(xx,yy), &!
                            & ';', wat_bal1%h_interc(xx,yy), ';',wat_bal1%h_runoff(xx,yy), ';',wat_bal1%h_inf(xx,yy), &!
                            & ';', wat_bal1%h_eva_pot(xx,yy), ';',wat_bal1%h_eva(xx,yy), &!
                            & ';', wat_bal1%h_transp_pot(xx,yy), ';',wat_bal1%h_transp_act(xx,yy), &! 
                            & ';', wat_bal1%h_perc(xx,yy), ';',wat_bal1%h_soil(xx,yy), &!
                            & ';', wat_bal1%h_pond(xx,yy), ';', wat_bal2%h_rise(xx,yy), &!
                            & ';', wat_bal2%h_transp_pot(xx,yy), ';',wat_bal2%h_transp_act(xx,yy), ';',wat_bal2%k_s(xx,yy), &!
                            & ';', wat_bal2%d_t(xx,yy), ';',wat_bal2%depth_under_rz(xx,yy), ';',wat_bal2%h_caprise(xx,yy), &!
                            & ';', wat_bal2%h_perc(xx,yy), ';',wat_bal2%h_soil(xx,yy),';',wat_bal2_old%h_soil(xx,yy), &!
                            & ';', wat_bal2%h_raw_sup(xx,yy), ';',wat_bal2%h_raw_inf(xx,yy), &!
                            & ';', info_spat%wat_tab%mat(xx,yy), &!
                            & ';', 0,';',0, & ! add dummy variables to maintain file structure
                            & ';', esp_perc(xx,yy,1),';',esp_perc(xx,yy,2),';',h_bypass(xx,yy)
                    case (1)
                        write(out_tbl%sample_cells(i)%file%unit,*)doy, &!
                            & ';', meteo%p(xx,yy), ';', meteo%T_max(xx,yy), &!
                            & ';', meteo%T_min(xx,yy), ';', meteo%et0(xx,yy), &!
                            & ';', pheno%k_cb(xx,yy), ';',pheno%lai(xx,yy), ';',pheno%p_day(xx,yy), &!
                            & ';', h_irr_sum(xx,yy),&!
                            & ';', wat_bal1%h_eff_rain(xx,yy), ';', wat_bal1%h_gross_av_water(xx,yy), ';',wat_bal1%h_net_av_water(xx,yy), &!
                            & ';', wat_bal1%h_interc(xx,yy), ';',wat_bal1%h_runoff(xx,yy), ';',wat_bal1%h_inf(xx,yy), &!
                            & ';', wat_bal1%h_eva_pot(xx,yy), ';',wat_bal1%h_eva(xx,yy), &!
                            & ';', wat_bal1%h_transp_pot(xx,yy), ';',wat_bal1%h_transp_act(xx,yy), &!
                            & ';', wat_bal1%h_perc(xx,yy), ';',wat_bal1%h_soil(xx,yy), &!
                            & ';', wat_bal1%h_pond(xx,yy), ';', wat_bal2%h_rise(xx,yy), &!
                            & ';', wat_bal2%h_transp_pot(xx,yy), ';',wat_bal2%h_transp_act(xx,yy), ';',wat_bal2%k_s(xx,yy), &!
                            & ';', wat_bal2%d_t(xx,yy), ';',wat_bal2%depth_under_rz(xx,yy), ';',wat_bal2%h_caprise(xx,yy), &!
                            & ';', wat_bal2%h_perc(xx,yy), ';',wat_bal2%h_soil(xx,yy),';',wat_bal2_old%h_soil(xx,yy), &!
                            & ';', wat_bal2%h_raw_sup(xx,yy), ';',wat_bal2%h_raw_inf(xx,yy), &!
                            & ';', info_spat%wat_tab%mat(xx,yy), &!
                            & ';', coll_irr(xx,yy),';',priv_irr(xx,yy),';',esp_perc(xx,yy,1),&
                            & ';', esp_perc(xx,yy,2),';',h_bypass(xx,yy)
                    case default
                end select
            end do!
        end if
        !
        if (debug .eqv. .true.) then
            write(out_tbl%et0_ws%unit,*) doy,'; ',(info_meteo(i)%et0,'; ',i=1,size(info_meteo))
            if (cells .eqv. .true.) then
                do i=1,size(out_tbl%cell_eva)!
                    xx=out_tbl%cell_eva(i)%coord%row!
                    yy=out_tbl%cell_eva(i)%coord%col!
                    if (h_irr_sum(xx,yy)/=0) then
                         write(out_tbl%cell_eva(i)%file%unit,*)doy, &!
                            & ';', meteo%Wind_vel(xx,yy), ';', meteo%RH_min(xx,yy), &!
                            & ';', meteo%et0(xx,yy), &!
                            & ';', pheno%k_cb(xx,yy), ';', pheno%h(xx,yy), &!
                            & ';', pars%irr%met(info_spat%irr_meth_id%mat(xx,yy))%f_wet, ';', wat%few(xx,yy), &!
                            & ';', pheno%f_c(xx,yy), &!
                            & ';', wat%kc_max(xx,yy), &!
                            & ';', wat%wat1_rew(xx,yy), ';', wat%layer(1)%h_wp(xx,yy), ';', wat_bal_hour%inten%h_soil1(xx,yy), &!
                            & ';', wat_bal_hour%esten%k_e(xx,yy), &!
                            & ';', wat_bal1%h_eva_pot(xx,yy), ';', wat_bal1%h_eva(xx,yy), ';', wat_bal_hour%esten%k_r(xx,yy), ';', fw_old(xx,yy) ! RR test
                    else
                         write(out_tbl%cell_eva(i)%file%unit,*)doy, &!
                            & ';', meteo%Wind_vel(xx,yy), ';', meteo%RH_min(xx,yy), &!
                            & ';', meteo%et0(xx,yy), &!
                            & ';', pheno%k_cb(xx,yy), ';', pheno%h(xx,yy), &!
                            & ';', fw(xx,yy), ';', wat%few(xx,yy), &!
                            & ';', pheno%f_c(xx,yy), &!
                            & ';', wat%kc_max(xx,yy), &!
                            & ';', wat%wat1_rew(xx,yy), ';', wat%layer(1)%h_wp(xx,yy), ';', wat_bal_hour%inten%h_soil1(xx,yy), &!
                            & ';', wat_bal_hour%esten%k_e(xx,yy), &!
                            & ';', wat_bal1%h_eva_pot(xx,yy), ';', wat_bal1%h_eva(xx,yy), ';', wat_bal_hour%esten%k_r(xx,yy), ';', fw_old(xx,yy) !%RR% test
                    end if
                end do!
                do i=1, size(out_tbl%cell_cn)!
                    xx=out_tbl%cell_cn(i)%coord%row!
                    yy=out_tbl%cell_cn(i)%coord%col!
                    write(out_tbl%cell_cn(i)%file%unit,*)doy, &!
                        & ';', info_spat%theta(1)%wp%mat(xx,yy)+info_spat%theta(2)%wp%mat(xx,yy), &
                        & ';', info_spat%theta(1)%fc%mat(xx,yy)+info_spat%theta(2)%fc%mat(xx,yy), &
                        & ';', info_spat%theta(1)%sat%mat(xx,yy)+info_spat%theta(2)%sat%mat(xx,yy), &
                        & ';', info_spat%theta(1)%wp%mat(xx,yy)+info_spat%theta(2)%wp%mat(xx,yy)+ &
                        & (2.0/3.0)*(info_spat%theta(1)%fc%mat(xx,yy)+info_spat%theta(2)%fc%mat(xx,yy)- &
                        & info_spat%theta(1)%wp%mat(xx,yy)-info_spat%theta(2)%wp%mat(xx,yy)), &
                        & ';', wat_bal1%t_soil(xx,yy)+wat_bal2%t_soil(xx,yy), &
                        & ';', out_cn(xx,yy)%tab_cn2, ';', out_cn(xx,yy)%tab_cn3, &
                        & ';', info_spat%slope%mat(xx,yy), &!
                        & ';', pheno%cn_class(xx,yy),';', pheno%cn_day(xx,yy), &!
                        & ';', out_cn(xx,yy)%tab_cn2_baresoil, ';', out_cn(xx,yy)%tab_cn2_slope, &!
                        & ';', out_cn(xx,yy)%cn1_day, ';', out_cn(xx,yy)%cn2_day, &!
                        & ';', out_cn(xx,yy)%cn3_day, ';', out_cn(xx,yy)%cn_day, &!
                        & ';', pars%sim%lambda_cn, &!
                        & ';', 25.4*((1000./out_cn_day(xx,yy))-10.), &!
                        & ';', pars%sim%lambda_cn*25.4*((1000./out_cn_day(xx,yy))-10.), &!
                        & ';', wat_bal1%h_gross_av_water(xx,yy), ';', wat_bal1%h_net_av_water(xx,yy), &!
                        & ';', ((wat_bal1%h_gross_av_water(xx,yy)-pars%sim%lambda_cn*25.4*((1000./out_cn_day(xx,yy))-10.))**2.)/ &
                        & (wat_bal1%h_gross_av_water(xx,yy)+0.8*25.4*((1000./out_cn_day(xx,yy))-10.)), &
                        & ';', wat_bal1%h_runoff(xx,yy)
                end do
            end if
        end if
    end subroutine write_daily_output

    subroutine write_outputs_by_step (doy, meteo, irrigation_sum, bil1, bil2, &
        & info_spat, coll_irr, priv_irr, asc, deb_asc, hbypass, intervals, clock_time, debug,summary)
        ! writes periodic (monthly/weekly/custom) output in *.asc files
        implicit none
        integer, intent(in):: doy
        type(meteo_mat), intent(in):: meteo
        real(dp), dimension(:,:), intent(in):: irrigation_sum
        real(dp), dimension(:,:), intent(in):: hbypass
        type(balance1_matrices), intent(in):: bil1
        type(balance2_matrices), intent(in):: bil2
        type(spatial_info), intent(in):: info_spat
        real(dp), dimension(:,:), intent(in):: coll_irr, priv_irr
        type(step_map), intent(in):: asc
        type(step_debug_map), intent(in):: deb_asc
        integer, dimension(:), intent(in):: intervals
        integer, intent(in):: clock_time
        logical, intent(in)::debug
        logical, intent(in)::summary

        asc%runoff%mat = asc%runoff%mat + bil1%h_runoff
        asc%rain%mat = asc%rain%mat + meteo%p!
        asc%transp_act%mat = asc%transp_act%mat + bil1%h_transp_act + bil2%h_transp_act!
        asc%transp_pot%mat = asc%transp_pot%mat + bil1%h_transp_pot + bil2%h_transp_pot!
        asc%irr%mat = asc%irr%mat + irrigation_sum + hbypass
        asc%irr_loss%mat = asc%irr_loss%mat + hbypass
        asc%cap_rise%mat = asc%cap_rise%mat + bil2%h_caprise!
        asc%irr_nm_priv%mat = asc%irr_nm_priv%mat + priv_irr!
        asc%irr_nm_col%mat = asc%irr_nm_col%mat + coll_irr!
        asc%deep_perc%mat = asc%deep_perc%mat + bil2%h_perc - bil2%h_caprise 
        asc%et_pot%mat = asc%et_pot%mat + bil1%h_eva_pot + bil1%h_transp_pot + bil2%h_transp_pot
        asc%et_act%mat = asc%et_act%mat + bil1%h_eva + bil1%h_transp_act + bil2%h_transp_act
        if (debug .eqv. .true.) THEN
            deb_asc%eva_act%mat = deb_asc%eva_act%mat + bil1%h_eva!
            deb_asc%eff_rain%mat = deb_asc%eff_rain%mat + bil1%h_eff_rain!
            deb_asc%perc1%mat = deb_asc%perc1%mat + bil1%h_perc!
            deb_asc%perc2%mat = deb_asc%perc2%mat + bil2%h_perc!
            deb_asc%h_soil1%mat = bil1%h_soil
            deb_asc%h_soil2%mat = bil2%h_soil
        end if
        
        if (summary .eqv. .false.) then
            call save_step_data(asc,doy,info_spat%domain,intervals,clock_time)
        else
            call save_step_irrigation(asc,doy,info_spat%domain,intervals,clock_time)
        end if
        
        if (debug .eqv. .true.) then
            call save_debug_step_data(deb_asc,doy,info_spat%domain,intervals,clock_time)
        end if

    end subroutine write_outputs_by_step

    subroutine allocate_all (asc,yasc,deb_asc,deb_yasc,bil1,bil1_old,bil2,bil2_old,bil_hour,meteo,wat,pheno,&
        & imax, jmax, domain, debug)
        integer,dimension(:,:),intent(in)::domain
        logical,intent(in)::debug
        type(step_map)::asc
        type(step_debug_map), optional:: deb_asc
        type(annual_map)::yasc
        type(annual_debug_map), optional:: deb_yasc
        type(balance1_matrices)::bil1, bil1_old!
        type(balance2_matrices)::bil2, bil2_old!
        type(hourly)::bil_hour
        type(meteo_mat)::meteo
        type(wat_matrix)::wat
        type(crop_pars_matrices)::pheno
        integer,intent(in)::imax!
        integer,intent(in)::jmax!
        logical::allocazione

        allocazione = .true.
        !allocazione della variabile per gli output *.asc!
        call init_step_output(asc,domain)
        call init_yearly_output(yasc,domain)
        if (debug .eqv. .true.) then
            call init_step_debug_output(deb_asc,domain)
            call init_yearly_debug_output(deb_yasc,domain)
        end if

        !alloca/dealloca le matrici in base al flag TRUE/FALSE!
        call init_wat_bal1_matrices(bil1,imax, jmax, allocazione)!
        call init_wat_bal1_matrices(bil1_old,imax, jmax, allocazione)!
        call init_wat_bal2_matrices(bil2,imax, jmax, allocazione)!
        call init_wat_bal2_matrices(bil2_old,imax, jmax, allocazione)!
        call init_wat_bal_hour(bil_hour,imax, jmax, allocazione)!
        call init_meteo_matrices(meteo,imax, jmax, allocazione)!
        call init_wat_matrices(wat,imax, jmax, allocazione)!
        call init_pheno_matrices(pheno,imax, jmax, allocazione)!
    end subroutine allocate_all

    subroutine destroy_all(asc,yasc,deb_asc,deb_yasc,bil1,bil1_old,bil2,bil2_old,bil_hour,meteo,wat,pheno,imax, jmax, debug)
        implicit none
        type(step_map)::asc!
        type(step_debug_map),optional::deb_asc!
        type(annual_map)::yasc!
        type(annual_debug_map),optional::deb_yasc!
        type(balance1_matrices)::bil1, bil1_old!
        type(balance2_matrices)::bil2, bil2_old!
        type(hourly)::bil_hour
        type(meteo_mat)::meteo
        type(wat_matrix)::wat
        type(crop_pars_matrices)::pheno
        integer,intent(in)::imax!
        integer,intent(in)::jmax!
        logical,intent(in)::debug
        logical::f_allocate

        f_allocate = .false.
        ! destroy reference to output files
        call destroy_step_output(asc)!
        call destroy_annual_output(yasc)
        if (debug .eqv. .true.) then
            call destroy_step_debug_output(deb_asc)
            call destroy_annual_debug_output(deb_yasc)
        end if
        ! allocate/destrot the matrix base on f_allocate flag (true/false)
        call init_wat_bal1_matrices(bil1,imax, jmax, f_allocate)!
        call init_wat_bal1_matrices(bil1_old,imax, jmax, f_allocate)!
        call init_wat_bal2_matrices(bil2,imax, jmax, f_allocate)!
        call init_wat_bal2_matrices(bil2_old,imax, jmax, f_allocate)!
        call init_wat_bal_hour(bil_hour,imax, jmax, f_allocate)!
        call init_meteo_matrices(meteo,imax, jmax, f_allocate)!
        call init_wat_matrices(wat,imax, jmax, f_allocate)!
        call init_pheno_matrices(pheno,imax, jmax, f_allocate)!
    end subroutine destroy_all

    subroutine allocate_crop_map(crop_mat,domain,mcrop_alt,a)
        ! init crop matrix
        implicit none
        type(crop_matrices)::crop_mat
        integer,dimension(:,:),intent(in)::domain
        integer,intent(in)::mcrop_alt
        integer,intent(in)::a
        
        allocate(crop_mat%ii0(size(domain,1),size(domain,2),mcrop_alt))
        allocate(crop_mat%iie(size(domain,1),size(domain,2),mcrop_alt))
        allocate(crop_mat%iid(size(domain,1),size(domain,2),mcrop_alt))
        allocate(crop_mat%ii0_ref(size(domain,1),size(domain,2),mcrop_alt))
        allocate(crop_mat%iie_ref(size(domain,1),size(domain,2),mcrop_alt))
        allocate(crop_mat%iid_ref(size(domain,1),size(domain,2),mcrop_alt))
        allocate(crop_mat%dij(size(domain,1),size(domain,2),mcrop_alt))
        allocate(crop_mat%TSP_high(size(domain,1),size(domain,2),mcrop_alt))
        allocate(crop_mat%TSP_low(size(domain,1),size(domain,2),mcrop_alt))
        allocate(crop_mat%wp_adj(size(domain,1),size(domain,2),mcrop_alt))
        allocate(crop_mat%HI(size(domain,1),size(domain,2),mcrop_alt))
        allocate(crop_mat%Ky_tot(size(domain,1),size(domain,2),mcrop_alt))
        allocate(crop_mat%Ky_pheno(size(domain,1),size(domain,2),mcrop_alt,4))
        allocate(crop_mat%T_crit(size(domain,1),size(domain,2),mcrop_alt))
        allocate(crop_mat%T_lim(size(domain,1),size(domain,2),mcrop_alt))
        allocate(crop_mat%k_cb_min(size(domain,1),size(domain,2),mcrop_alt))
        allocate(crop_mat%k_cb_mid(size(domain,1),size(domain,2),mcrop_alt))
        allocate(crop_mat%k_cb_max(size(domain,1),size(domain,2),mcrop_alt))
        crop_mat%ii0 = a
        crop_mat%iie = a
        crop_mat%iid = 0.d0
        crop_mat%ii0_ref = a
        crop_mat%iie_ref = a
        crop_mat%iid_ref = a
        crop_mat%dij = a
        crop_mat%TSP_high = a
        crop_mat%TSP_low = a
        crop_mat%wp_adj = a
        crop_mat%HI = a
        crop_mat%Ky_tot = a
        crop_mat%Ky_pheno = a
        crop_mat%T_crit = a
        crop_mat%T_lim = a
        crop_mat%k_cb_min = a
        crop_mat%k_cb_mid = a
        crop_mat%k_cb_max = a
    end subroutine allocate_crop_map

    subroutine destroy_crop(crop_map)
        implicit none!
        type(crop_matrices),intent(inout)::crop_map
        
        deallocate(crop_map%ii0)
        deallocate(crop_map%iie)
        deallocate(crop_map%iid)
        deallocate(crop_map%ii0_ref)
        deallocate(crop_map%iie_ref)
        deallocate(crop_map%iid_ref)
        deallocate(crop_map%dij)
        deallocate(crop_map%TSP_high)
        deallocate(crop_map%TSP_low)
    end subroutine destroy_crop

    subroutine init_wat_bal1(bil,a)
        ! overloading of the assignment operator "="
        implicit none!
        real(dp),intent(in)::a!
        type(balance1_matrices),intent(out)::bil!
        !!
        bil%d_e = a!
        bil%h_eva = a!
        bil%h_eva_pot = a!
        bil%h_transp_act = a!
        bil%h_transp_pot = a!
        bil%h_soil = a!
        bil%t_soil = a!
        bil%h_interc = a!
        bil%h_perc = a!
        bil%h_inf = a!
        bil%h_eff_rain = a!
        bil%h_net_av_water = a!
        bil%h_runoff = a!
        bil%h_pond = a
    end subroutine init_wat_bal1!
    !
    subroutine init_wat_bal2(bil,a)!
        ! overloading of the assignment operator "="
        implicit none!
        real(dp),intent(in)::a!
        type(balance2_matrices),intent(out)::bil!
        !!
        bil%d_t = a!
        bil%h_soil = a!
        bil%t_soil = a!
        bil%h_transp_act = a!
        bil%h_transp_pot = a!
        bil%h_perc = a!
        bil%h_raw_sup = a!
        bil%h_raw = a!
        bil%h_raw_inf = a!
        bil%h_raw_priv = a!
        bil%k_s = a!
        bil%depth_under_rz = a!
        bil%h_caprise = a!
        bil%h_rise = a!
    end subroutine init_wat_bal2!
    !
    subroutine eq_wat_bal1(bil_out,bil_in)!
        ! overloading of the assignment operator "="
        implicit none!
        type(balance1_matrices),intent(in)::bil_in!
        type(balance1_matrices),intent(out)::bil_out!
        bil_out%d_e = bil_in%d_e!
        bil_out%h_eva = bil_in%h_eva!
        bil_out%h_eva_pot = bil_in%h_eva_pot!
        bil_out%h_transp_act = bil_in%h_transp_act!
        bil_out%h_transp_pot = bil_in%h_transp_pot!
        bil_out%h_soil = bil_in%h_soil!
        bil_out%h_interc = bil_in%h_interc!
        bil_out%t_soil = bil_in%t_soil!
        bil_out%h_perc = bil_in%h_perc!
        bil_out%h_inf = bil_in%h_inf!
        bil_out%h_eff_rain = bil_in%h_eff_rain!
        bil_out%h_net_av_water = bil_in%h_net_av_water!
        bil_out%h_runoff = bil_in%h_runoff!
        bil_out%h_pond=bil_in%h_pond
    end subroutine eq_wat_bal1!
    !
    subroutine eq_wat_bal2(bil_out,bil_in)!
        implicit none!
        ! overloading of the assignment operator "="
        type(balance2_matrices),intent(in)::bil_in!
        type(balance2_matrices),intent(out)::bil_out!
        !!
        bil_out%d_t = bil_in%d_t!
        bil_out%h_soil = bil_in%h_soil!
        bil_out%t_soil = bil_in%t_soil!
        bil_out%h_transp_act = bil_in%h_transp_act!
        bil_out%h_transp_pot = bil_in%h_transp_pot!
        bil_out%h_perc = bil_in%h_perc!
        bil_out%h_raw_sup = bil_in%h_raw_sup!
        bil_out%h_raw = bil_in%h_raw!
        bil_out%h_raw_inf = bil_in%h_raw_inf!
        bil_out%h_raw_priv = bil_in%h_raw_priv!
        bil_out%k_s = bil_in%k_s!
        bil_out%depth_under_rz = bil_in%depth_under_rz!
        bil_out%h_caprise = bil_in%h_caprise!
        bil_out%h_rise = bil_in%h_rise!
    end subroutine eq_wat_bal2!
    !
    subroutine init_pheno(pheno,a)!
        ! TODO: probably to be updated
        implicit none!
        type(crop_pars_matrices),intent(inout)::pheno!
        real(dp),intent(in)::a!
        !!
        pheno%k_cb_old = pheno%k_cb   ! Store kcb of the current day in kcb_old
        ! pheno%kcb_low = pheno%kcb_low       

        pheno%k_cb = a!
        pheno%h = a!
        pheno%d_r = a!
        pheno%lai = a!
        pheno%cn_day = int(a)!

        pheno%irrigation_class = int(a)!
        pheno%cn_class = int(a)!
        pheno%p = a!
        pheno%a = a!
        pheno%max_RF_t = a!
        pheno%RF_e = a!
        pheno%RF_t = a!
        pheno%T_lim = a!
        pheno%T_crit = a!
        pheno%HI = a!
        pheno%Ky_tot = a!
        pheno%Ky_pheno = a!
        pheno%k_cb_mid = a!
        pheno%k_cb_high = a!
        pheno%wp_adj = a!
        pheno%p_day = a!
        pheno%r_stress = a!
    end subroutine init_pheno!
    !
    subroutine init_meteo(meteo,a)!
        implicit none!
        type(meteo_mat),intent(out)::meteo!
        real(dp),intent(in)::a!
        !!
        meteo%T_max = a!
        meteo%T_min = a!
        meteo%P = a!
        meteo%P_cum = a!
        meteo%RH_max = a!
        meteo%RH_min = a!
        meteo%Wind_vel = a!
        meteo%Rad_sol = a!
        meteo%lat = a!
        meteo%alt = a!
        meteo%et0 = a!
        meteo%T_ave = a!
    end subroutine init_meteo!
    !
    subroutine create_meteo_matrices(info_meteo, dir_meteo, meteo_weight, meteo, domain, doy, res_canopy)!
        ! distribute weather variables to the domain according to the weights of each weather stations
        implicit none!
        type(meteo_info),dimension(:),intent(in)::info_meteo!
        integer,dimension(:,:,:),intent(in)::dir_meteo!
        type(meteo_mat),intent(inout)::meteo!
        real(dp),dimension(:,:,:),intent(in)::meteo_weight!
        type(grid_i),intent(in)::domain!
        integer,intent(in)::doy!
        real(dp),intent(in)::res_canopy
        integer::i,j,k!
        !!
        !inizializzazione della variabile meteo!
        meteo = 0.0D0!
        do k=1,size(meteo_weight,3)!
            forall(i=1:size(domain%mat,1),j=1:size(domain%mat,2),domain%mat(i,j)/=domain%header%nan)!
                            meteo%T_max(i,j)    = info_meteo(dir_meteo(i,j,k))%T_max*meteo_weight(i,j,k) + meteo%T_max(i,j)!
                            meteo%T_min(i,j)    = info_meteo(dir_meteo(i,j,k))%T_min*meteo_weight(i,j,k) + meteo%T_min(i,j)!
                            meteo%P(i,j)        = info_meteo(dir_meteo(i,j,k))%P*meteo_weight(i,j,k) + meteo%P(i,j)!
                            meteo%P_cum(i,j)    = info_meteo(dir_meteo(i,j,k))%P_cum*meteo_weight(i,j,k) + meteo%P_cum(i,j)!
                            meteo%RH_max(i,j)   = info_meteo(dir_meteo(i,j,k))%RH_max*meteo_weight(i,j,k) + meteo%RH_max(i,j)!
                            meteo%RH_min(i,j)   = info_meteo(dir_meteo(i,j,k))%RH_min*meteo_weight(i,j,k) + meteo%RH_min(i,j)!
                            meteo%Wind_vel(i,j) = info_meteo(dir_meteo(i,j,k))%wind_vel*meteo_weight(i,j,k) + meteo%Wind_vel(i,j)!
                            meteo%Rad_sol(i,j)  = info_meteo(dir_meteo(i,j,k))%sol_rad*meteo_weight(i,j,k) + meteo%Rad_sol(i,j)!
                            meteo%lat(i,j)      = info_meteo(dir_meteo(i,j,k))%lat_deg*meteo_weight(i,j,k) + meteo%lat(i,j)!
                            meteo%alt(i,j)      = info_meteo(dir_meteo(i,j,k))%alt_m*meteo_weight(i,j,k) + meteo%alt(i,j)!
            end forall!
        end do!
        ! calculate ET0 from distributed parameters
        meteo%et0 = ET_reference(meteo%T_max, meteo%T_min, meteo%RH_max, meteo%RH_min, meteo%Wind_vel, meteo%Rad_sol,&
                                       meteo%lat, meteo%alt, res_canopy, doy, domain%header%imax, domain%header%jmax)!
    end subroutine create_meteo_matrices!
    
    function calc_interception(p,pheno)!
        ! Calculate the interception according to Von Hoyningen-Hune (1983) and Braden (1985)
        implicit none!
        real(dp),dimension(:,:),intent(in)::p       ! precipitation and any other above canopy irrigation [mm]
        type(crop_pars_matrices),intent(in)::pheno!
        real(dp),dimension(size(p,1),size(p,2))::f_c  !cover fraction [-]!
        real(dp),dimension(size(p,1),size(p,2))::calc_interception!
        !!
        calc_interception = 0.!
        where(p>0.)
            where (pheno%lai>0.)
                ! limit f_c to 1 in order to not have interception greater than precipitation
                f_c=min(1.0d0, pheno%lai/3.) ! %EAC%: fix "Different type kinds" error
                calc_interception = pheno%a * pheno%lai * (1.-(1./(1.+((f_c * p)/(pheno%a * pheno%lai)))))!
            end where!
        end where
    end function calc_interception!
    !
    function net_precipitation(h_gross_precip, h_interception)!
        ! calculate the effective precipitation of each day
        implicit none!
        real(dp),dimension(:,:),intent(in)::h_gross_precip                       ! precipitation + irrigation above canopy [mm]
        real(dp),dimension(size(h_gross_precip,1),size(h_gross_precip,2))::h_interception       ! interception [mm] !TODO: allocated externally?
        real(dp),dimension(size(h_gross_precip,1),size(h_gross_precip,2))::net_precipitation!
        !!
        net_precipitation = 0.
        where(h_gross_precip>0.)
            net_precipitation=h_gross_precip-h_interception!
            where ((net_precipitation<0.) .and. (net_precipitation>-1E-05))
                ! delete rounding error (except for big negative values that are errors!!!)
                net_precipitation=0.
            end where
        end where!
    end function net_precipitation!
    !
    subroutine eq_extensive(a,b)!
        ! overloading of the assignment operator for the variable of type hourly
        implicit none!
        real(dp),intent(in)::b!
        type(hourly),intent(out)::a!
        !!
        a%esten%k_e = b!
        a%esten%k_r = b!
        a%esten%h_eva = b!
        a%esten%h_eva_pot = b!
        a%esten%h_inf = b!
        a%esten%h_perc1 = b!
        a%esten%h_perc2 = b!
        a%esten%h_pond = b!
        a%esten%h_transp_act1 = b!
        a%esten%h_transp_pot1 = b!
        a%esten%h_transp_act2 = b!
        a%esten%h_transp_pot2 = b!
        a%esten%h_caprise = b!
        a%esten%h_rise = b!
        a%n_iter1 = int(b)
        a%n_iter2 = int(b)
        a%n_max1 = 1!
        a%n_max2 = 1!
    end subroutine eq_extensive!
    !
    !
    subroutine b1_no_iter_eva(pheno, meteo, h_rain_lim, wat, fw_day, fw_irr, fw_rain, fc, &
                            & h_irr_sum, f_interception, domain, balance1_mat, fw_old) ! %RR% add fw_old for testing
        ! calculate the elements of the evaporative model that change daily
        implicit none!
        type(grid_i),intent(in)::domain!
        TYPE(balance1_matrices),INTENT(IN):: balance1_mat
        type(crop_pars_matrices),intent(in)::pheno!
        type(meteo_mat),intent(in)::meteo!
        type(wat_matrix),intent(inout)::wat!
        real(dp),dimension(:,:),intent(in)::h_irr_sum
        real(dp),dimension(:,:),intent(in)::f_interception !
        real(dp),dimension(:,:),intent(in)::fw_irr   ! fraction of soil wetted during irrigation event
        real(dp),dimension(:,:),intent(inout)::fw_day, fw_old ! fw daily updated - %RR% test
        real(dp),intent(in):: fw_rain                ! fraction of soil wetted during a precipitation event
        real(dp),parameter::kc_min=0.175             ! minimum value of Kc of bare soil (it changes between 0.15 a 0.20)[-]
        real(dp),dimension(size(domain%mat,1),size(domain%mat,2)):: Kc_max1, Kc_max2, fc!
        real(dp),dimension(size(domain%mat,1),size(domain%mat,2)):: V_eva, TEW, h_irr_net, h_rain_net, fw_day_tmp !, fw_old test
        real(dp),dimension(size(domain%mat,1),size(domain%mat,2)):: zero_mat ! to set minimum value of V_eva - %RR%
        real(dp),intent(in)::h_rain_lim !minimum meaningfull value of precipitation

        zero_mat = 0.
        fw_old = fw_day +0.
        fw_day_tmp = fw_day + 0.

        where(domain%mat/=domain%header%nan)!
            ! init fw
            ! fw = fw_irr      if the surface is wetted by irrigation
            ! fw = (fw_irr*net_irrigation + fwCost*net_rain)/(net_irr+net_rain)           if the surface is wetted by irrigation and rain, fw is the weigthed average of
            !                                                                                                                                           fw_irr and fw_rain (=1) based on infiltration depths
            ! fw = 1                if the surface is wetted by significant rain (i.e. >= 5mm)
            ! if fw_old > fw_day --> weigthed average between fw_old and fw_day based on the water content of the 1st layer - Gandolfi 29/5

            h_irr_net = h_irr_sum - calc_interception(h_irr_sum*f_interception, pheno) !
            h_rain_net = meteo%P - calc_interception(meteo%p, pheno)

            V_eva = max(zero_mat, (balance1_mat%h_soil - wat%layer(1)%h_wp ))
            TEW = wat%layer(1)%h_sat - wat%layer(1)%h_wp

            WHERE(h_irr_sum/=0 .and. meteo%P>= h_rain_lim) !!
                fw_day = (fw_irr * h_irr_net+ fw_rain * h_rain_net)/(h_irr_net + h_rain_net)
            ELSE WHERE(h_irr_sum/=0)
                fw_day = fw_irr
            eLSE WHERE (meteo%P>= h_rain_lim) !
                fw_day = fw_rain
            ELSE WHERE
                fw_day = fw_old
            END WHERE

            fw_day_tmp = fw_day

            WHERE(fw_old>fw_day_tmp) !!
                fw_day = (fw_old * V_eva + fw_day_tmp * min(V_eva + h_irr_net, TEW))/ &
                &                         (V_eva + min(V_eva + h_irr_net, TEW))
            END WHERE

            ! FAO56 eq. 23
            ! Kc_max is the upper limit for evaporation and transpiration
            ! linked to the available energy (it ranges from 1.05 to 1.30) [-]!
            ! after the precipitation or irrigation event
            Kc_max1 = 1.2+(0.04*(meteo%Wind_vel-2) - 0.004*(meteo%RH_min-45))*((pheno%h/3)**0.3)
            Kc_max2 = pheno%k_cb + 0.05 ! always greater than k_cb also in case of complete cover
            wat%kc_max = max(Kc_max1,Kc_max2)
            
            ! FAO56 eq. 76
            ! few  =  exposed and wetted soil fraction [-] !
            ! the evaporation depends from the weetted surface and it is greater with higher fraction of bare soil
            where(pheno%f_c<0.)!
                ! fc calculated with equation 76 [FAO56 p.149] unless it isn't an input %RR%
                WHERE(pheno%k_cb<=kc_min)
                    fc=0.!
                ELSE WHERE
                    fc = ((pheno%k_cb-kc_min)/(wat%kc_max-kc_min))**(1+0.5*pheno%h)
                END WHERE
            else where
                fc = pheno%f_c
            END WHERE

            pheno%f_c = fc

            ! FAO56 eq. 75
            wat%few=min(1-fc,fw_day)!
            !(3) Kr!
            wat%wat1_rew=wat%layer(1)%h_fc-wat%layer(1)%rew ! water content at REW
        end where!

    end subroutine b1_no_iter_eva!
    !
    subroutine init_water_balance_variables(wat_bal1,wat_bal2)!
        implicit none!
        type(balance1_matrices),intent(out)::wat_bal1!
        type(balance2_matrices),intent(out)::wat_bal2!
        !!
        wat_bal1%h_soil = 0.!
        wat_bal1%h_inf = 0.!
        wat_bal1%h_eva = 0.!
        wat_bal1%h_eva_pot = 0.!
        wat_bal1%h_transp_act = 0.!
        wat_bal1%h_transp_pot = 0.!
        wat_bal1%h_perc = 0.!
        wat_bal1%h_runoff = 0.!

        wat_bal2%h_soil = 0.!
        wat_bal2%h_transp_act = 0.!
        wat_bal2%h_transp_pot = 0.!
        wat_bal2%h_perc = 0.!
        wat_bal2%k_s = 0.!
        wat_bal2%h_caprise = 0.!
        wat_bal2%h_rise = 0.!

    end subroutine init_water_balance_variables!
    !
    subroutine x_wat(domain,teta,ze,zr,wat,theta2_rice,k_CN,Kcb)
        ! update water matrix that change with only the thikness of the soil layer
        implicit none!
        type(grid_i),intent(in)::domain!
        type(moisture),dimension(:),intent(in)::teta!
        real(dp),dimension(:,:),intent(in)::ze,zr!
        type(wat_matrix),intent(out)::wat!
        type(soil2_rice),intent(in)::theta2_rice
        integer,dimension(:,:),intent(in)::k_CN
        real(dp),dimension(:,:),intent(in)::Kcb
        !!
        where(domain%mat/=domain%header%nan)!
            wat%layer(1)%h_wp  = 1000*teta(1)%wp%mat*ze             ! water soil content at WP [mm]
            wat%layer(1)%h_fc  = 1000*teta(1)%fc%mat*ze             ! water soil content at FC [mm]
            wat%layer(1)%h_sat = 1000*teta(1)%sat%mat*ze            ! water soil content at saturation [mm] %AB%
            wat%layer(1)%h_r   = 1000*teta(1)%r%mat*ze              ! water soil content at residual humidity [mm]
            wat%layer(1)%rew = ((teta(1)%fc%mat - 0.5*teta(1)%wp%mat)*0.4)*1000*ze!
            where (k_CN ==7 .and. Kcb>0) ! %AB% rice in crop season
                wat%layer(2)%h_wp=1000*theta2_rice%theta2_WP*zr     ! water soil content at WP [mm]
                wat%layer(2)%h_fc=1000*theta2_rice%theta2_FC*zr     ! water soil content at FC [mm]
                wat%layer(2)%h_sat=1000*theta2_rice%theta2_SAT*zr   ! water soil content at saturation [mm] %AB%
                wat%layer(2)%h_r=1000*theta2_rice%theta2_R*zr       ! water soil content at residual humidity [mm]
            else where
                wat%layer(2)%h_wp=1000*teta(2)%wp%mat*zr            ! water soil content at WP [mm]
                wat%layer(2)%h_fc=1000*teta(2)%fc%mat*zr            ! water soil content at FC [mm]
                wat%layer(2)%h_sat=1000*teta(2)%sat%mat*zr          ! water soil content at saturation [mm] %AB%
                wat%layer(2)%h_r=1000*teta(2)%r%mat*zr              ! water soil content at residual humidity [mm]
            end where
        end where!
    end subroutine x_wat!
    !
    subroutine init_wat_bal1_matrices(wat_bal1,imax,jmax,f_allocate)!
        ! init/destroy water balance variable
        implicit none!
        type(balance1_matrices),intent(inout)::wat_bal1!
        integer,intent(in)::imax!
        integer,intent(in)::jmax!
        logical,intent(in)::f_allocate
        integer::checkstat!
        character (len=*),parameter:: error_message = "wat_ball has been wrongly allocated"
        !!
        if(f_allocate)then!
            allocate(wat_bal1%d_e(imax,jmax),stat=checkstat)           ; if(checkstat/=0)print*,error_message!
            allocate(wat_bal1%h_eva(imax,jmax),stat=checkstat)          ; if(checkstat/=0)print*,error_message!
            allocate(wat_bal1%h_eva_pot(imax,jmax),stat=checkstat)      ; if(checkstat/=0)print*,error_message!
            allocate(wat_bal1%h_transp_act(imax,jmax),stat=checkstat)    ; if(checkstat/=0)print*,error_message!
            allocate(wat_bal1%h_transp_pot(imax,jmax),stat=checkstat)    ; if(checkstat/=0)print*,error_message!
            allocate(wat_bal1%h_soil(imax,jmax),stat=checkstat)          ; if(checkstat/=0)print*,error_message!
            allocate(wat_bal1%t_soil(imax,jmax),stat=checkstat)         ; if(checkstat/=0)print*,error_message!
            allocate(wat_bal1%h_interc(imax,jmax),stat=checkstat)       ; if(checkstat/=0)print*,error_message!
            allocate(wat_bal1%h_perc(imax,jmax),stat=checkstat)         ; if(checkstat/=0)print*,error_message!
            allocate(wat_bal1%h_inf(imax,jmax),stat=checkstat)          ; if(checkstat/=0)print*,error_message!
            allocate(wat_bal1%h_eff_rain(imax,jmax),stat=checkstat)     ; if(checkstat/=0)print*,error_message!
            allocate(wat_bal1%h_net_av_water(imax,jmax),stat=checkstat)  ; if(checkstat/=0)print*,error_message!
            allocate(wat_bal1%h_gross_av_water(imax,jmax),stat=checkstat); if(checkstat/=0)print*,error_message!
            allocate(wat_bal1%h_runoff(imax,jmax),stat=checkstat)         ; if(checkstat/=0)print*,error_message!
            allocate(wat_bal1%h_pond(imax,jmax),stat=checkstat)         ; if(checkstat/=0)print*,error_message!
            allocate(wat_bal1%k_s_dry(imax,jmax),stat=checkstat)           ; if(checkstat/=0)print*,error_message!
            allocate(wat_bal1%k_s_sat(imax,jmax),stat=checkstat)           ; if(checkstat/=0)print*,error_message!
            allocate(wat_bal1%k_s(imax,jmax),stat=checkstat)           ; if(checkstat/=0)print*,error_message!
            !
        else!
            deallocate(wat_bal1%d_e) !
            deallocate(wat_bal1%h_eva) !
            deallocate(wat_bal1%h_eva_pot) !
            deallocate(wat_bal1%h_transp_act) !
            deallocate(wat_bal1%h_transp_pot) !
            deallocate(wat_bal1%h_soil) !
            deallocate(wat_bal1%t_soil) !
            deallocate(wat_bal1%h_interc) !
            deallocate(wat_bal1%h_perc) !
            deallocate(wat_bal1%h_inf) !
            deallocate(wat_bal1%h_eff_rain) !
            deallocate(wat_bal1%h_net_av_water) !
            deallocate(wat_bal1%h_gross_av_water) !
            deallocate(wat_bal1%h_runoff)!
            deallocate(wat_bal1%h_pond)!
            deallocate(wat_bal1%k_s_dry)
            deallocate(wat_bal1%k_s_sat)
            deallocate(wat_bal1%k_s)
        end if!
    end subroutine init_wat_bal1_matrices!
    !
    subroutine init_wat_bal2_matrices(wat_bal2,imax,jmax,f_allocate)!
        ! init/destroy water balance variable
        implicit none!
        type(balance2_matrices),intent(inout)::wat_bal2!
        integer,intent(in)::imax!
        integer,intent(in)::jmax!
        logical,intent(in)::f_allocate!
        !!
        integer::checkstat!
        character (len=*),parameter:: error_message = "wat_bal2 has been wrongly allocated"
        !!
        if(f_allocate)then!
            allocate(wat_bal2%d_t(imax,jmax),stat=checkstat)           ; if(checkstat/=0)print*,error_message!
            allocate(wat_bal2%h_soil(imax,jmax),stat=checkstat)          ; if(checkstat/=0)print*,error_message!
            allocate(wat_bal2%t_soil(imax,jmax),stat=checkstat)         ; if(checkstat/=0)print*,error_message!
            allocate(wat_bal2%h_transp_act(imax,jmax),stat=checkstat)    ; if(checkstat/=0)print*,error_message!
            allocate(wat_bal2%h_transp_pot(imax,jmax),stat=checkstat)    ; if(checkstat/=0)print*,error_message!
            allocate(wat_bal2%h_perc(imax,jmax),stat=checkstat)         ; if(checkstat/=0)print*,error_message!
            allocate(wat_bal2%h_raw_sup(imax,jmax),stat=checkstat)       ; if(checkstat/=0)print*,error_message!
            allocate(wat_bal2%h_raw(imax,jmax),stat=checkstat)          ; if(checkstat/=0)print*,error_message!
            allocate(wat_bal2%h_raw_inf(imax,jmax),stat=checkstat)       ; if(checkstat/=0)print*,error_message!
            allocate(wat_bal2%h_raw_priv(imax,jmax),stat=checkstat)        ; if(checkstat/=0)print*,error_message!
            allocate(wat_bal2%k_s_dry(imax,jmax),stat=checkstat)           ; if(checkstat/=0)print*,error_message!
            allocate(wat_bal2%k_s_sat(imax,jmax),stat=checkstat)           ; if(checkstat/=0)print*,error_message!
            allocate(wat_bal2%k_s(imax,jmax),stat=checkstat)           ; if(checkstat/=0)print*,error_message!
            allocate(wat_bal2%depth_under_rz(imax,jmax),stat=checkstat); if(checkstat/=0)print*,error_message!
            allocate(wat_bal2%h_caprise(imax,jmax),stat=checkstat)      ; if(checkstat/=0)print*,error_message!
            allocate(wat_bal2%h_rise(imax,jmax),stat=checkstat)         ; if(checkstat/=0)print*,error_message!
        else!
            deallocate(wat_bal2%d_t) !
            deallocate(wat_bal2%h_soil) !
            deallocate(wat_bal2%t_soil) !
            deallocate(wat_bal2%h_transp_act) !
            deallocate(wat_bal2%h_transp_pot) !
            deallocate(wat_bal2%h_perc) !
            deallocate(wat_bal2%h_raw_sup) !
            deallocate(wat_bal2%h_raw) !
            deallocate(wat_bal2%h_raw_inf) !
            deallocate(wat_bal2%h_raw_priv) !
            deallocate(wat_bal2%k_s_dry) !
            deallocate(wat_bal2%k_s_sat) !
            deallocate(wat_bal2%k_s) !
            deallocate(wat_bal2%depth_under_rz) !
            deallocate(wat_bal2%h_caprise) !
            deallocate(wat_bal2%h_rise) !
        end if!
    end subroutine init_wat_bal2_matrices!
    !
    subroutine init_wat_bal_hour(wat_bal1_hour,imax,jmax,f_allocate)!
        ! init/destroy wat_bal_hourly
        implicit none!
        type(hourly),intent(inout)::wat_bal1_hour!
        integer,intent(in)::imax!
        integer,intent(in)::jmax!
        logical,intent(in)::f_allocate!
        integer::checkstat!
        character (len=*),parameter:: error_message = "wat_bal1_hour has been wrongly allocated"

        if(f_allocate .eqv. .true.)then!
            allocate(wat_bal1_hour%esten%k_e           (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,error_message
            allocate(wat_bal1_hour%esten%k_r           (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,error_message ! %RR% added
            allocate(wat_bal1_hour%esten%h_eva         (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,error_message
            allocate(wat_bal1_hour%esten%h_eva_pot     (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,error_message
            allocate(wat_bal1_hour%esten%h_inf         (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,error_message
            allocate(wat_bal1_hour%esten%h_perc1       (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,error_message
            allocate(wat_bal1_hour%esten%h_perc2      (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,error_message
            allocate(wat_bal1_hour%esten%h_pond        (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,error_message
            allocate(wat_bal1_hour%esten%h_transp_act1  (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,error_message
            allocate(wat_bal1_hour%esten%h_transp_pot1  (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,error_message
            allocate(wat_bal1_hour%esten%h_transp_act2 (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,error_message
            allocate(wat_bal1_hour%esten%h_transp_pot2 (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,error_message
            allocate(wat_bal1_hour%esten%h_caprise     (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,error_message
            allocate(wat_bal1_hour%esten%h_rise        (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,error_message
            allocate(wat_bal1_hour%esten%h_net_av_water (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,error_message
            allocate(wat_bal1_hour%esten%h_gross_av_water (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,error_message
            allocate(wat_bal1_hour%esten%h_eff_rain    (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,error_message
            allocate(wat_bal1_hour%inten%h_soil1        (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,error_message
            allocate(wat_bal1_hour%inten%h_soil2        (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,error_message
            allocate(wat_bal1_hour%inten%k_s_dry        (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,error_message
            allocate(wat_bal1_hour%inten%k_s_sat        (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,error_message
            allocate(wat_bal1_hour%inten%k_s          (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,error_message
            allocate(wat_bal1_hour%inten%h_pond0       (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,error_message
            allocate(wat_bal1_hour%n_iter1            (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,error_message
            allocate(wat_bal1_hour%n_iter2            (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,error_message
            allocate(wat_bal1_hour%n_max1              (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,error_message
            allocate(wat_bal1_hour%n_max2              (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,error_message
        else!
            deallocate(wat_bal1_hour%esten%k_e      )!
            deallocate(wat_bal1_hour%esten%k_r      )! %RR% added
            deallocate(wat_bal1_hour%esten%h_eva      )!
            deallocate(wat_bal1_hour%esten%h_eva_pot      )!
            deallocate(wat_bal1_hour%esten%h_inf      )!
            deallocate(wat_bal1_hour%esten%h_perc1    )!
            deallocate(wat_bal1_hour%esten%h_perc2   )!
            deallocate(wat_bal1_hour%esten%h_pond   )!
            deallocate(wat_bal1_hour%esten%h_transp_act1    ) !
            deallocate(wat_bal1_hour%esten%h_transp_pot1  )!
            deallocate(wat_bal1_hour%esten%h_transp_act2    ) !
            deallocate(wat_bal1_hour%esten%h_transp_pot2  )!
            deallocate(wat_bal1_hour%esten%h_caprise  )!
            deallocate(wat_bal1_hour%esten%h_rise  )!
            deallocate(wat_bal1_hour%esten%h_net_av_water )!
            deallocate(wat_bal1_hour%esten%h_gross_av_water)!
            deallocate(wat_bal1_hour%esten%h_eff_rain     )!
            deallocate(wat_bal1_hour%inten%h_soil1)     !
            deallocate(wat_bal1_hour%inten%h_soil2)     !
            deallocate(wat_bal1_hour%inten%k_s_dry  )     !
            deallocate(wat_bal1_hour%inten%k_s_sat  )     !
            deallocate(wat_bal1_hour%inten%k_s  )     !
            deallocate(wat_bal1_hour%inten%h_pond0   )!
            deallocate(wat_bal1_hour%n_iter1)         !
            deallocate(wat_bal1_hour%n_iter2)         !
            deallocate(wat_bal1_hour%n_max1)         !
            deallocate(wat_bal1_hour%n_max2)         !
        end if!
    end subroutine init_wat_bal_hour!
    !
    subroutine init_meteo_matrices(meteo,imax,jmax,f_allocate)!
        ! init/destroy weather map
        implicit none!
        type(meteo_mat),intent(inout)::meteo!
        integer,intent(in)::imax!
        integer,intent(in)::jmax!
        logical,intent(in)::f_allocate!
        integer::checkstat!
        character (len=*),parameter:: errormessage = "meteo has been wrongly allocated"
        !!
        if(f_allocate)then!
            allocate(meteo%T_max    (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,errormessage!
            allocate(meteo%T_min    (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,errormessage!
            allocate(meteo%P       (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,errormessage!
            allocate(meteo%P_cum    (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,errormessage!
            allocate(meteo%RH_max    (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,errormessage!
            allocate(meteo%RH_min    (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,errormessage!
            allocate(meteo%Wind_vel(imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,errormessage!
            allocate(meteo%Rad_sol      (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,errormessage!
            allocate(meteo%lat     (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,errormessage!
            allocate(meteo%alt     (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,errormessage!
            allocate(meteo%et0     (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,errormessage!
            allocate(meteo%T_ave    (imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,errormessage!
        else!
            deallocate(meteo%T_max    )!
            deallocate(meteo%T_min    )!
            deallocate(meteo%P       )!
            deallocate(meteo%P_cum   )!
            deallocate(meteo%RH_max    )!
            deallocate(meteo%RH_min    )!
            deallocate(meteo%Wind_vel)!
            deallocate(meteo%Rad_sol      )!
            deallocate(meteo%lat     )!
            deallocate(meteo%alt     )!
            deallocate(meteo%et0     )!
            deallocate(meteo%T_ave    )!
        end if!
    end subroutine init_meteo_matrices!
    !
    subroutine init_pheno_matrices(pheno,imax,jmax,f_allocate)!
        ! init/destroy pheno
        implicit none!
        type(crop_pars_matrices),intent(inout)::pheno!
        integer,intent(in)::imax!
        integer,intent(in)::jmax!
        logical,intent(in)::f_allocate!
        integer,parameter::phases=4
        integer::checkstat!
        character (len=*),parameter:: errormessage = "pheno has been wrongly allocated"
        !!
        if(f_allocate)then!
            allocate(pheno%k_cb_old        (imax,jmax),stat=checkstat)        ; if(checkstat/=0)print*,errormessage!
            allocate(pheno%k_cb            (imax,jmax),stat=checkstat)        ; if(checkstat/=0)print*,errormessage!
            allocate(pheno%h              (imax,jmax),stat=checkstat)        ; if(checkstat/=0)print*,errormessage!
            allocate(pheno%d_r             (imax,jmax),stat=checkstat)        ; if(checkstat/=0)print*,errormessage!
            allocate(pheno%lai            (imax,jmax),stat=checkstat)        ; if(checkstat/=0)print*,errormessage!
            allocate(pheno%cn_day         (imax,jmax),stat=checkstat)        ; if(checkstat/=0)print*,errormessage!
             allocate(pheno%f_c            (imax,jmax),stat=checkstat)        ; if(checkstat/=0)print*,errormessage!
            allocate(pheno%irrigation_class (imax,jmax),stat=checkstat)        ; if(checkstat/=0)print*,errormessage!
            allocate(pheno%cn_class             (imax,jmax),stat=checkstat)        ; if(checkstat/=0)print*,errormessage!
            allocate(pheno%p              (imax,jmax),stat=checkstat)        ; if(checkstat/=0)print*,errormessage!
            allocate(pheno%a              (imax,jmax),stat=checkstat)        ; if(checkstat/=0)print*,errormessage!
            allocate(pheno%max_d_r          (imax,jmax),stat=checkstat)        ; if(checkstat/=0)print*,errormessage!
            allocate(pheno%max_RF_t         (imax,jmax),stat=checkstat)        ; if(checkstat/=0)print*,errormessage!
            allocate(pheno%RF_e            (imax,jmax),stat=checkstat)        ; if(checkstat/=0)print*,errormessage!
            allocate(pheno%RF_t            (imax,jmax),stat=checkstat)        ; if(checkstat/=0)print*,errormessage!
            allocate(pheno%T_lim           (imax,jmax),stat=checkstat)        ; if(checkstat/=0)print*,errormessage!
            allocate(pheno%T_crit          (imax,jmax),stat=checkstat)        ; if(checkstat/=0)print*,errormessage!
            allocate(pheno%HI             (imax,jmax),stat=checkstat)        ; if(checkstat/=0)print*,errormessage!
            allocate(pheno%Ky_tot            (imax,jmax),stat=checkstat)        ; if(checkstat/=0)print*,errormessage!
            allocate(pheno%k_cb_low        (imax,jmax),stat=checkstat)        ; if(checkstat/=0)print*,errormessage!
            allocate(pheno%k_cb_mid        (imax,jmax),stat=checkstat)        ; if(checkstat/=0)print*,errormessage!
            allocate(pheno%k_cb_high       (imax,jmax),stat=checkstat)        ; if(checkstat/=0)print*,errormessage!
            allocate(pheno%wp_adj          (imax,jmax),stat=checkstat)        ; if(checkstat/=0)print*,errormessage!
            allocate(pheno%p_day           (imax,jmax),stat=checkstat)        ; if(checkstat/=0)print*,errormessage!
            allocate(pheno%n_crop_in_year    (imax,jmax),stat=checkstat)        ; if(checkstat/=0)print*,errormessage!
            allocate(pheno%pheno_idx   (imax,jmax),stat=checkstat)        ; if(checkstat/=0)print*,errormessage!
            allocate(pheno%Ky_pheno            (imax,jmax,phases),stat=checkstat) ; if(checkstat/=0) print*,errormessage! 3d matrix
            allocate(pheno%r_stress           (imax,jmax),stat=checkstat)        ; if(checkstat/=0)print*,errormessage!
            
        else!
            deallocate(pheno%k_cb_old        )!
            deallocate(pheno%k_cb            )!
            deallocate(pheno%h              )!
            deallocate(pheno%d_r             )!
            deallocate(pheno%lai            )!
            deallocate(pheno%cn_day         )!
            deallocate(pheno%f_c            )!
            deallocate(pheno%irrigation_class )  !
            deallocate(pheno%cn_class             )!
            deallocate(pheno%p              )!
            deallocate(pheno%a              )!
            deallocate(pheno%max_d_r          )!
            deallocate(pheno%max_RF_t         )!
            deallocate(pheno%RF_e            )!
            deallocate(pheno%RF_t            )!
            deallocate(pheno%T_lim           )!
            deallocate(pheno%T_crit          )!
            deallocate(pheno%HI             )!
            deallocate(pheno%Ky_tot            )!
            deallocate(pheno%Ky_pheno            )!
            deallocate(pheno%k_cb_low        )!
            deallocate(pheno%k_cb_mid        )!
            deallocate(pheno%k_cb_high       )!
            deallocate(pheno%wp_adj          )!
            deallocate(pheno%p_day           )!
            deallocate(pheno%n_crop_in_year    )!
            deallocate(pheno%pheno_idx   )!
            deallocate(pheno%r_stress   )!
        end if!
    end subroutine init_pheno_matrices
    !
    subroutine init_wat_matrices(wat,imax,jmax,allocazione)!
        ! init/destroy water matrices
        implicit none!
        type(wat_matrix),intent(inout)::wat!
        integer,intent(in)::imax!
        integer,intent(in)::jmax!
        logical,intent(in)::allocazione!
        integer::checkstat!
        character (len=50):: errormessage = "wat has been wrongly allocated"
        !!
        if(allocazione)then!
            allocate(wat%layer(1)%h_wp(imax,jmax),stat=checkstat)   ; if(checkstat/=0)print*,errormessage
            allocate(wat%layer(1)%h_fc(imax,jmax),stat=checkstat)   ; if(checkstat/=0)print*,errormessage
            allocate(wat%layer(1)%h_sat(imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,errormessage
            allocate(wat%layer(1)%h_r(imax,jmax),stat=checkstat)    ; if(checkstat/=0)print*,errormessage
            allocate(wat%layer(1)%rew(imax,jmax),stat=checkstat)  ; if(checkstat/=0)print*,errormessage
            allocate(wat%layer(2)%h_wp(imax,jmax),stat=checkstat)   ; if(checkstat/=0)print*,errormessage
            allocate(wat%layer(2)%h_fc(imax,jmax),stat=checkstat)   ; if(checkstat/=0)print*,errormessage
            allocate(wat%layer(2)%h_sat(imax,jmax),stat=checkstat) ; if(checkstat/=0)print*,errormessage
            allocate(wat%layer(2)%h_r(imax,jmax),stat=checkstat)    ; if(checkstat/=0)print*,errormessage
            allocate(wat%theta2_rice(1)%h_wp(imax,jmax),stat=checkstat)    ; if(checkstat/=0)print*,errormessage
            allocate(wat%theta2_rice(1)%h_fc(imax,jmax),stat=checkstat)    ; if(checkstat/=0)print*,errormessage
            allocate(wat%theta2_rice(1)%h_sat(imax,jmax),stat=checkstat)  ; if(checkstat/=0)print*,errormessage
            allocate(wat%theta2_rice(1)%h_r(imax,jmax),stat=checkstat)     ; if(checkstat/=0)print*,errormessage
            allocate(wat%wat1_rew(imax,jmax),stat=checkstat)       ; if(checkstat/=0)print*,errormessage
            allocate(wat%few(imax,jmax),stat=checkstat)            ; if(checkstat/=0)print*,errormessage
            allocate(wat%kc_max(imax,jmax),stat=checkstat)          ; if(checkstat/=0)print*,errormessage
        else!
            deallocate(wat%layer(1)%h_wp)!
            deallocate(wat%layer(1)%h_fc)!
            deallocate(wat%layer(1)%h_sat)!
            deallocate(wat%layer(1)%h_r)!
            deallocate(wat%layer(1)%rew)!
            deallocate(wat%layer(2)%h_wp)!
            deallocate(wat%layer(2)%h_fc)!
            deallocate(wat%layer(2)%h_sat)!
            deallocate(wat%layer(2)%h_r)!
            deallocate(wat%theta2_rice(1)%h_wp)!
            deallocate(wat%theta2_rice(1)%h_fc)!
            deallocate(wat%theta2_rice(1)%h_sat)!
            deallocate(wat%theta2_rice(1)%h_r)!
            deallocate(wat%wat1_rew)!
            deallocate(wat%few)!
            deallocate(wat%kc_max)!
        end if!
    end subroutine init_wat_matrices!

end module cli_simulation_manager!
