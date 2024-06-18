module mod_meteo!
    use mod_constants, only: sp, dp
    use mod_utility, only: date, lower_case, days_x_month, calc_doy, seek_un, split_date
    use mod_parameters, only: simulation, par_method
    use mod_evapotranspiration, only: ET_reference
    implicit none!
    ! store weather station data
    type meteo_info!
        integer::unit                   ! unit associated to the input file
        character(len=255)::filename    ! name of the input file
        real(dp)::x_m                   ! longitude
        real(dp)::y_m                   ! latitude
        type(date)::start               ! first day of the time serie
        type(date)::finish              ! last day of the time serie
        real(dp)::alt_m                 ! altitude [m]
        real(dp)::lat_deg               ! latitude [degree]
        real(dp)::T_max                 ! maximum daily temperature [°C]!
        real(dp)::T_min                 ! minimum daily temperature [°C]!
        real(dp)::P                     ! daily precipitation [mm]!
        real(dp)::P_cum                 ! daily cumulative precipitation [mm]!
        real(dp)::RH_max                ! maximum volumetric air moisture [%]!
        real(dp)::RH_min                ! minimum volumetric air moisture [%]!
        real(dp)::wind_vel              ! wind velocity at 2 m height [m/s]!
        real(dp)::sol_rad               ! daily solar radiation [MJ m-2 day-1]!
        real(dp)::et0                   ! reference et0!
    end type meteo_info!

    ! It contains the weather parameters spatially distributed
    ! The distribution is made over the domain with the weights defined in info_spat%weigth_ws(:)
    type meteo_mat
        real(dp),dimension(:,:),pointer::T_ave      ! Daily average temperature [°C]
        real(dp),dimension(:,:),pointer::T_max      ! maximum daily temperature [°C]
        real(dp),dimension(:,:),pointer::T_min      ! minimum daily temperature [°C]
        real(dp),dimension(:,:),pointer::P          ! precipitation [mm]
        real(dp),dimension(:,:),pointer::P_cum      ! cumulate precipitation [mm]
        real(dp),dimension(:,:),pointer::RH_max     ! maximum daily air moisture [%]
        real(dp),dimension(:,:),pointer::RH_min     ! minimum daily air moisture [%]
        real(dp),dimension(:,:),pointer::Wind_vel   ! wind velocity at 2 m height [m/s]
        real(dp),dimension(:,:),pointer::Rad_sol    ! global solar radiation [MJ m^-2 d^-1]
        real(dp),dimension(:,:),pointer::lat        ! Latitude [m]
        real(dp),dimension(:,:),pointer::alt        ! Altitude [m a.s.l.]
        real(dp),dimension(:,:),pointer::et0        ! evapotranspiration [mm]
    end type meteo_mat

    !
    contains!
    !
    subroutine meteo_series_length(sim, debug)!
        ! calculate the length of the weather time series [days]
        ! check if time series have the same length
        implicit none!
        type(simulation),intent(INout)::sim
        logical,optional,intent(in)::debug
        type(meteo_info),dimension(:),allocatable::info_meteo!
        integer::k!
        real(dp)::value!
        integer::count,gg_count,gg_in_yy, gg_diff
        integer,dimension(12)::calendario!

        call read_meteo_parameters(sim,info_meteo,debug)!
        
        ! Set simulation dates, if not already set
        if(sim%start_simulation%doy==calc_doy(29,02,1600)) sim%start_simulation = info_meteo(1)%start
        if(sim%end_simulation%doy==calc_doy(29,02,1600))   sim%end_simulation = info_meteo(1)%finish
        
        ! Verify that meteorological series are coherent
        if(any(info_meteo(:)%start%doy/=info_meteo(1)%start%doy))then
            stop 'Meteorological series start in different days. Execution will be aborted...'
        end if
        if(any(info_meteo(:)%finish%doy/=info_meteo(1)%finish%doy))then
            stop 'Meteorological series end in different days. Execution will be aborted...'
        end if

        ! Verify that meteorological series are coherent with simulation dates
        if(any(info_meteo(:)%start%doy>sim%start_simulation%doy)) then
            stop 'Meteorological series are not coherent with simulation dates. Execution will be aborted...'
        end if       
        if(any(info_meteo(:)%finish%doy<sim%end_simulation%doy)) then
            stop 'Meteorological series are not coherent with simulation dates. Execution will be aborted...'
        end if       
        
        ! Check length of time series
        do k = 1,sim%n_voronoi   ! TODO: n_voronoi or n_ws?
            gg_diff = info_meteo(k)%finish%doy - info_meteo(k)%start%doy + 1
            count=0!
            do while (.true.)
                read (info_meteo(k)%unit, *, end=999) value
                count=count+1!
            end do
            999 continue
            if (count /= gg_diff) then
                print *,"Meteorological series have incoherent lengths"
                print *,"The range between dates is ", gg_diff, " while data series is of ", count, " elements"
                print *, 'Execution will be aborted...'
                stop!
            end if!

            if(k==1)then!
                gg_count=count!
            else
                if(count/=gg_count)then!
                    stop 'Meteorological series have different lengths. Execution will be aborted...'
                end if
            end if
        end do

        ! calculate the number of years in the time series
        sim%meteo_years=0; gg_in_yy=0!
        sim%start_year = info_meteo(1)%start%year ! first year
        do!
            sim%meteo_years=sim%meteo_years+1!
            call days_x_month(calendario,sim%start_year+sim%meteo_years-1)!
            gg_in_yy=gg_in_yy+sum(calendario)!
            if(gg_in_yy>=gg_count) exit!
        end do!
        sim%sim_years=sim%meteo_years
        !
        ! Updating gg_count to take into account the end of simulation date
        if(sim%end_simulation%doy<info_meteo(1)%finish%doy) then
            gg_count = sim%end_simulation%doy - info_meteo(1)%start%doy + 1
            ! calculate the number of years to be simulated
            sim%sim_years=0; gg_in_yy=0!
            sim%start_year = info_meteo(1)%start%year ! first year
            do!
                sim%sim_years=sim%sim_years+1!
                call days_x_month(calendario,sim%start_year+sim%sim_years-1)!
                gg_in_yy=gg_in_yy+sum(calendario)!
                if(gg_in_yy>=gg_count) exit!
            end do!
        end if
        !
        allocate(sim%year_step(sim%sim_years))!
        
        sim%year_step=0!

        ! calculate the number of days for each simulation years
        do k=1,sim%sim_years!
            if(k/=sim%sim_years)then!
                if (info_meteo(1)%start%month > 2) then
                    call days_x_month(calendario,sim%start_year+k)!
                else
                    call days_x_month(calendario,sim%start_year+k-1)!
                end if
                sim%year_step(k)=sum(calendario)!
            else ! the last year fo the dataset could be uncompleted 
                sim%year_step(k)=gg_count-sum(sim%year_step)!
            end if!
        end do!
        call close_meteo_file(info_meteo)!
        !
        !debug!
        if(debug .eqv. .true.)then!
            print *,'===== DEBUG: meteo series length ====='
            print*,"Total number of days: ",gg_count!
            print*,"First simulation year: ",sim%start_year!
            print*,"Simulation length: ",sim%sim_years, " years"
            print*,"Number of days for each year:"!
            do k=1,sim%sim_years!
                print*," ",sim%year_step(k),"<---",sim%start_year+k-1!
            end do!
            print *,'===== END DEBUG ====='
        end if!
    end subroutine meteo_series_length!
    !
    subroutine read_meteo_parameters(sim,info_meteo,debug)!
        ! read the number  of weather station in the file "weather_stations.dat" and init the array
        ! for each station, it will be read: the name of the station, the latitude, the altitude, the starting day of the time series
        implicit none!
        type(simulation),intent(inout)::sim!
        type(meteo_info),dimension(:),allocatable,intent(inout)::info_meteo!
        logical, intent(in)::debug
        integer::errorflag!
        !!
        integer :: free_unit, i!
        integer :: ios
        character(len=255) :: filemeteo_name!
        integer::m_topoieta   ! check information in simulation
        character(len=255)::dir!
        character(len=300) :: comment,buffer, label
        integer :: p
        integer :: line, tablestart
        character(len=300) :: date_string
        character(leN=300) :: date_start, date_end

        ios = 0
        line = 0; tablestart = 0
        dir=trim(sim%meteo_path)
        errorflag=0!
        call seek_un( ErrorFlag, free_unit) ! search for free unit
        filemeteo_name=trim(sim%ws_list_fn)

        open( unit = free_unit, file = filemeteo_name, status = 'old', action="read", iostat = ios )!
        if (ios /= 0 ) then
            print *, 'Cannot open file ', filemeteo_name, '. The specified file does not exist. Execution will be aborted...'
            stop
        end if

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
                        case ('statnum')
                            read(buffer, *, iostat=ios) m_topoieta ! number of weather station
                            if (m_topoieta /= sim%n_voronoi) then
                                print *, "Meteorological stations number (MeteoStatTotNum) in simulation parameter &
                                    & file is not equal to the number (StatNum) in ", filemeteo_name, " database"
                                print *, 'Execution will be aborted...'
                                stop
                            End if
                            allocate(info_meteo(sim%n_voronoi))    ! allocate enough memory to store weather stations information
                        case ('table')
                            tablestart = line
                            read (free_unit, *) ! skip the header of the table
                            do i=1,size(info_meteo)!
                                read(free_unit,*)info_meteo(i)%filename, info_meteo(i)%x_m, info_meteo(i)%y_m!
                                line = line + 1
                                call seek_un(errorflag,info_meteo(i)%unit)!
                                open( unit=info_meteo(i)%unit, file=trim(dir)//trim(info_meteo(i)%filename), status='old',&
                                    & action="read", iostat=ios )!
                                if (ios /= 0 ) then
                                    print *, 'Cannot open file ', trim(info_meteo(i)%filename), '. The specified file does &
                                        & not exist. Execution will be aborted...'
                                    stop
                                end if
                                read(info_meteo(i)%unit,*)  ! skip the header with  Id and name of the weather station
                                read(info_meteo(i)%unit,*)  info_meteo(i)%lat_deg, info_meteo(i)%alt_m ! read latitude and altitude 
                                ! Read the first and the last days, separated by " -> " : e.g. "01/01/1993 -> 31/12/2014" TODO: reverse the order aaaa/mm/dd
                                read(info_meteo(i)%unit, '(a300)') date_string
                                call split_date(date_string, date_start, date_end)
                                call split_date(date_start, info_meteo(i)%start)
                                call split_date(date_end, info_meteo(i)%finish)
                                info_meteo(i)%start%doy = &
                                    & calc_doy(info_meteo(i)%start%day,info_meteo(i)%start%month,info_meteo(i)%start%year)
                                info_meteo(i)%finish%doy = &
                                    & calc_doy(info_meteo(i)%finish%day,info_meteo(i)%finish%month,info_meteo(i)%finish%year)
                                read(info_meteo(i)%unit,*)  ! skip line
                            end do
                        case ('endtable')
                            if ((line - tablestart - 1) > m_topoieta) then
                                stop 'Meteorological stations number (StatNum) is higher than meteorological &
                                    & station listed in the table. Execution will be aborted...'
                            else if ((line - tablestart - 1) < m_topoieta) then
                                stop 'Meteorological stations number (StatNum) is lower than meteorological &
                                    & station listed in the table. Execution will be aborted...'
                            end if
                        case default ! all other cases ... !
                            print *, 'Skipping invalid label <',trim(label),'> at line', line, ' of file: ', trim(filemeteo_name)
                            print *, 'Execution will be aborted...'
                            stop
                    end select
                end if
            end if
        end do
        close(free_unit)!

        if (debug .eqv. .true.) then
            print *,'===== DEBUG: weather station data====='
            print *,'Weather station file: ',  filemeteo_name
            print *,'Meteo files path: ',  dir
            print *,'# Weather station in list file: ',  size(info_meteo)
            print *,'name lon_m lat_m lat_g alt_m dd/mm/yy'
            do i=1,size(info_meteo)
                print *, trim(info_meteo(i)%filename),' ' , info_meteo(i)%x_m,' ', info_meteo(i)%y_m, ' ', &
                    info_meteo(i)%lat_deg, ' ', info_meteo(i)%alt_m, ' ',  &
                    info_meteo(i)%start%day, '/', info_meteo(i)%start%month, '/', info_meteo(i)%start%year
            end do
            print *,'===== END DEBUG ====='
        end if

    end subroutine read_meteo_parameters!
    !
    subroutine close_meteo_file(info_meteo)!
        ! close all files with weather data and set memory free 
        implicit none!
        type(meteo_info),dimension(:),allocatable,intent(inout)::info_meteo!
        !
        integer::i!
        !
        do i=1,size(info_meteo)!
            close(info_meteo(i)%unit)!
        end do!
        deallocate(info_meteo)!
    end subroutine close_meteo_file!
    !
    subroutine read_meteo_data(info_meteo,current_doy,res_surf, forecast_day)
        ! read meteo data and update ET0
        implicit none!
        integer,intent(in)::current_doy!
        integer,intent(in)::forecast_day!
        type(meteo_info),dimension(:),intent(inout)::info_meteo!
        real(dp),intent(in)::res_surf ! surface resistance
        real(dp) :: dummy(3) = 0.0D0
        !!
        integer::i, ii
        integer::IOstatus
        !!
        do i=1,size(info_meteo)!
            read(info_meteo(i)%unit,*,iostat=IOstatus) &
                info_meteo(i)%T_max,info_meteo(i)%T_min,info_meteo(i)%P,info_meteo(i)%RH_max, &!
                info_meteo(i)%RH_min,info_meteo(i)%wind_vel,info_meteo(i)%sol_rad!
            ! init and populate cumulate precipitation
            info_meteo(i)%P_cum = info_meteo(i)%P 
            do ii=1,forecast_day
                read(info_meteo(i)%unit,*,iostat=IOstatus) dummy
                info_meteo(i)%P_cum = info_meteo(i)%P_cum + dummy(3)
            end do
            ! go back
            do ii=1,forecast_day
                backspace(info_meteo(i)%unit,iostat=IOstatus)
            end do
            
        end do!

        ! calculate ET0!
        do i=1,size(info_meteo)!
            info_meteo(i)%et0=ET_reference(info_meteo(i)%T_max,info_meteo(i)%T_min,info_meteo(i)%RH_max, &!
                                            info_meteo(i)%RH_min,info_meteo(i)%wind_vel,info_meteo(i)%sol_rad,&
                                            info_meteo(i)%lat_deg,info_meteo(i)%alt_m,res_surf,current_doy)!
        end do!
    end subroutine read_meteo_data!
    !
        !
end module mod_meteo!
