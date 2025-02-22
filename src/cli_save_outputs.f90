! TODO: some fixes
module cli_save_outputs!
    use mod_constants, only: sp, dp
    use mod_utility, only: seek_un, lower_case
    use mod_grid, only: grid_i, grid_r, print_mat_as_grid
    use mod_common, only: spatial_info
    use mod_parameters
    use cli_watsources
    use mod_crop_phenology
    implicit none!

    type output!
        integer::unit               ! unit linked to the file
        character(len=255)::fn      ! name of the file
    end type output!

    type coordinate!
        integer::row                 ! row cell coordinate
        integer::col                 ! col cell coordinate
    end type coordinate!

    type cell_output!
        type(output)::file          !
        integer::num                ! cell id (only for reading)
        type(coordinate)::coord     !
    end type cell_output!

    type output_table_list!
        type(cell_output),dimension(:),pointer::sample_cells    ! attribute of each selected cells
        type(output)::et0_ws                                    ! ET0 weather stations
        type(output)::q_irr_units                               ! daily discharge for each irrigation units [m3/s]
        type(output)::q_surplus                                 ! daily discharge not used for irrigation for each irrigation units [m3/s]
        type(output)::q_irr                                     ! daily discharge used for irrigation for each irrigation units [m3/s]
        type(output)::q_rem                                     ! daily discharge transferred to the following day for each irrigation units [m3/s]
                                                                ! only volume not used for irrigation and less than the water depth required for a single cell are considered
        type(output)::q_un_priv                                 ! daily discharge estimated for private non-monitored sources [m3/s]
        type(output)::n_inv_cells                               ! number of investigated cells in one day
        type(output)::q_un_col                                  ! daily discharge estimated for collective non-monitored sources [m3/s]
        type(cell_output),dimension(:),pointer::cell_info       ! pedo-hydrological parameters of the sampled cell
        type(cell_output),dimension(:),pointer::cell_conv       ! model resolution convergence
        type(cell_output),dimension(:),pointer::prod_info       ! productivity and crop parameters
        type(cell_output),dimension(:),pointer::cell_eva        ! evaporation process outputs
        type(cell_output),dimension(:),pointer::cell_cn         ! CN & runoff process outputs
    end type output_table_list!

    type out_2d_mat
        character(len=255)::fn                                  ! name of the file
        real(dp),dimension(:,:),pointer::mat                    ! single data matrix
    end type out_2d_mat!

    type out_3d_mat
        character(len=255)::fn                                  ! name of the file
        real(dp),dimension(:,:,:),pointer::mat                  ! data matrix per layer 
    end type out_3d_mat!
    
    type out_4d_mat
        character(len=255)::fn                                  ! name of the file
        real(dp),dimension(:,:,:,:),pointer::mat                ! data matrix per layer and condiction
    end type out_4d_mat!
    
    ! TODO: check meaning
    type yield_t
        type(out_3d_mat):: biomass_pot                   ! potential biomass
        type(out_3d_mat):: yield_pot                     ! potential yield [t/ha]
        type(out_3d_mat):: yield_act                     ! actual yield [t/ha]
        type(out_3d_mat):: f_WS                          ! water-stress yield reduction factor - overall [-]
        type(out_3d_mat):: f_WS_stage                    ! water-stress yield reduction factor - stages [-]
        type(out_3d_mat):: f_HS                          ! heat-stress yield reduction factor [-]
        type(out_3d_mat):: f_HS_sum                      ! somma dei valori per il calcolo del fattore di riduzione correlato allo stress
        type(out_3d_mat):: transp_ratio_sum              ! somma del rapporto tra traspirazione potenziale ed evapotraspirazione di riferimento
        type(out_4d_mat):: T_act_sum                     ! somma della traspirazione effettiva per ciascuna fase del kcb
        type(out_4d_mat):: T_pot_sum                     ! somma della traspirazione potenziale per ciascuna fase del kcb
        type(out_4d_mat):: dev_stage                     ! durata delle fasi del kcb
    end type yield_t

    type step_map!
        ! store the spatial distribution of several parameters aggregated by month
        type(out_2d_mat)::rain          ! rain [mm]
        type(out_2d_mat)::transp_act    ! actual transpiration [mm]
        type(out_2d_mat)::transp_pot    ! potential transpiration [mm]
        type(out_2d_mat)::irr           ! irrigation [mm]
        type(out_2d_mat)::irr_loss      ! irrigation loss [mm]
        type(out_2d_mat)::cap_rise      ! capillary rise [mm]
        type(out_2d_mat)::irr_nm_priv   ! irrigation from unmonitored private sources [mm]
        type(out_2d_mat)::irr_nm_col    ! irrigation from unmonitored collective sources [mm]
        type(out_2d_mat)::deep_perc     ! deep percolation to ground water [mm]
        type(out_2d_mat)::runoff        ! runoff [mm]
        type(out_2d_mat)::et_pot        ! potential evapotranspiration [mm]
        type(out_2d_mat)::et_act        ! actual evapotranspiration [mm]
    end type step_map!

    type step_debug_map
        ! store the spatial distribution of several parameters aggregated by month for DEBUG
        type(out_2d_mat)::eva_act       ! actual evaporation [mm]
        type(out_2d_mat)::perc1         ! percolation from the first layer [mm]
        type(out_2d_mat)::eff_rain      ! effective rain [mm]
        type(out_2d_mat)::perc2         ! percolation from the second layer [mm]!
        type(out_2d_mat)::h_soil1       ! soil water content of the first layer [mm]
        type(out_2d_mat)::h_soil2       ! soil water content of the second layer [mm]
    end type step_debug_map!

    type annual_map
        ! store the spatial distribution of several parameters aggregated by year
        type(out_2d_mat)::rain                  ! rain [mm]
        type(out_2d_mat)::rain_crop_season      ! rain during the crop season [mm]
        type(out_2d_mat)::irr                   ! irrigation [mm]
        type(out_2d_mat)::irr_loss              ! irrigation losses [mm]
        type(out_2d_mat)::eva_pot_crop_season   ! potential evaporation during the crop season [mm]
        type(out_2d_mat)::eva_act_crop_season   ! actual evaporation during the crop season [mm]
        type(out_2d_mat)::transp_act            ! actual transpiration [mm]
        type(out_2d_mat)::transp_pot            ! potential transpiration [mm]
        type(out_2d_mat)::runoff                ! total runoff [mm]
        type(out_2d_mat)::net_flux_gw           ! net flux to ground water [mm]
        type(out_2d_mat)::total_eff             ! efficiency of irrigation and rain [-]
        type(out_2d_mat)::n_irr_events          ! number of irrigation events [-]
        type(out_2d_mat)::h_irr_mean            ! average irrigation depth [mm]
    end type annual_map!

    type annual_debug_map
        ! store the spatial distribution of several parameters aggregated by year for DEBUG
        type(out_2d_mat)::eva_act_tot       ! total actual evaporation [mm]
        type(out_2d_mat)::rain_eff          ! rain efficiency
        type(out_2d_mat)::iter1             ! maximum number of iteration of the first layer [-]
        type(out_2d_mat)::iter2             ! maximum number of iteration of the second layer [-]
    end type annual_debug_map!

    
    interface assignment(=)!
        module procedure assign_step_map!
        module procedure assign_step_debug_map!
        module procedure assign_annual_map
        module procedure assign_annual_debug_map
        module procedure assign_yield_map
    end interface!
    !
    contains!
    
    ! TODO: make a general form
    subroutine init_cell_output_file(free_unit,file_name,table_header,table_subheader)!
        ! create the output file and print the header
        implicit none!
        character(len=*),intent(in)::file_name!
        character(len=*),intent(in)::table_header                               ! column names
        character(len=*),dimension(:),optional,intent(in)::table_subheader      ! additional column names
        integer,intent(out)::free_unit!
        integer::errorflag,ios,i!
        !
        errorflag=0!
        ios=0!
        !
        call seek_un(errorflag,free_unit)    ! in utility
        open(free_unit,file=file_name,action="write",iostat=ios)!
        if(ios/=0)then!
            print*, "Error opening file '",file_name,"' connected to unit ",free_unit," iostat=",ios, &
                & ". Execution will be aborted..."!
            stop
        end if!
        if(.not. present(table_subheader))then!
            write(free_unit,'(1a)')table_header
        else!
            write(free_unit,*)(trim(table_subheader(i)),'; ',i=1,size(table_subheader))!
        end if!
    end subroutine init_cell_output_file!
    
    subroutine init_cell_output_by_year(out_tbl_list,path,yr,id_ws_list, mode, f_cell_exists, debug, id_irr_unit_list, n_nm_pub, id_nm_pub_list)!
        ! read the list of output cells and
        ! prepare the output file for the selected year
        character(len=*),intent(in)::path                                   ! output path
        character(len=*),intent(in)::yr                                     ! year
        character(len=*),dimension(:),intent(in)::id_ws_list                ! weather stations
        integer,dimension(:),intent(in),optional::id_irr_unit_list          ! irrigation units
        integer, intent(in),optional:: n_nm_pub                             ! number of non monitored public sources
        character(len=*),dimension(:),intent(in),optional::id_nm_pub_list   ! id list of non monitored public sources
        type(output_table_list),intent(inout)::out_tbl_list                 ! list of the output table to be printed
        integer, intent(in)::mode                                           ! simulation mode
        logical,intent(in)::f_cell_exists,debug
        character(len=100),dimension(:),allocatable::id_irr_unit_str!
        character(len=100),dimension(:),allocatable::id_ws_str
        character(len=100),dimension(:),allocatable::w_str!
        integer::freeunit,errorflag,ios!
        integer::n_cells,i!
        character(len=55)::row_str,col_str!
        character(len=255)::title!
        character(len=300) :: comment, buffer, label
        integer :: line, table_start, p

        ! TODO: move to a separate subroutine
        line = 0
        if (f_cell_exists .eqv. .true.) then
            ! coordinates of the sample cells (deallocated at the end of each year)
            call seek_un(errorflag,freeunit)
            open(freeunit,file='cells.txt',action="read",status="old",iostat=ios)!
            if(ios/=0)then!
                print*," Error opening file 'cells.txt' connected to unit ",freeunit, &
                    & " iostat=",ios, ". Execution will be aborted..."
                stop
            end if!
            do while (ios == 0)
                read (freeunit, '(A)', iostat = ios) buffer
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
                            case ('ncells')
                                read (buffer, *, iostat = ios) n_cells
                                allocate(out_tbl_list%sample_cells(n_cells),out_tbl_list%cell_info(n_cells),out_tbl_list%prod_info(n_cells))
                                if (debug .eqv. .true.) then
                                    allocate(out_tbl_list%cell_conv(n_cells))!
                                    allocate(out_tbl_list%cell_eva(n_cells))!
                                    allocate(out_tbl_list%cell_cn(n_cells))!
                                end if
                            case ('table')
                                table_start = line
                                read (freeunit, *) ! skip the table header
                                do i=1, n_cells
                                    !~ read (freeunit, *) xls%celle(i)%num,xls%celle(i)%coor%cy,xls%celle(i)%coor%cx!
                                    read (freeunit, *) out_tbl_list%sample_cells(i)%num,out_tbl_list%sample_cells(i)%coord%row,out_tbl_list%sample_cells(i)%coord%col!
                                    line = line + 1
                                end do
                            case ('endtable')
                                if ((line - table_start - 1) > n_cells) then
                                    print *, 'The number of cells is higher than cell listed. Execution will be aborted...'
                                    stop
                                else if ((line - table_start - 1) < n_cells) then
                                    print *, 'The number of cells is lower than cell listed. Execution will be aborted...'
                                    stop
                                end if    
                            case default
                                print *, "Skipping invalid label <",trim(label),"> at line", line, " of file: 'cells.txt'. &
                                    & Execution will be aborted..."
                                stop
                        end select
                    end if
                end if
            end do
            close(freeunit)

            out_tbl_list%cell_info%coord%row=out_tbl_list%sample_cells%coord%row!
            out_tbl_list%cell_info%coord%col=out_tbl_list%sample_cells%coord%col!
            out_tbl_list%prod_info%coord%row=out_tbl_list%sample_cells%coord%row!
            out_tbl_list%prod_info%coord%col=out_tbl_list%sample_cells%coord%col!
            if (debug .eqv. .true.) then
                out_tbl_list%cell_conv%coord%row=out_tbl_list%sample_cells%coord%row!
                out_tbl_list%cell_conv%coord%col=out_tbl_list%sample_cells%coord%col!
                out_tbl_list%cell_eva%coord%row=out_tbl_list%sample_cells%coord%row!
                out_tbl_list%cell_eva%coord%col=out_tbl_list%sample_cells%coord%col!
                out_tbl_list%cell_cn%coord%row=out_tbl_list%sample_cells%coord%row!
                out_tbl_list%cell_cn%coord%col=out_tbl_list%sample_cells%coord%col!
            end if
            !
            ! save file header and connect to the control cell
            do i=1,size(out_tbl_list%sample_cells)!
                write(row_str,*)out_tbl_list%sample_cells(i)%coord%row!
                write(col_str,*)out_tbl_list%sample_cells(i)%coord%col!
                out_tbl_list%sample_cells(i)%file%fn = trim(adjustl(yr))//'_cell_'//trim(adjustl(row_str)) &!
                    & //'_'//trim(adjustl(col_str))//'.csv'!
                call init_cell_output_file(out_tbl_list%sample_cells(i)%file%unit,trim(path)//trim(adjustl(out_tbl_list%sample_cells(i)%file%fn)), &!
                    & 'Giulian_day'//'; '//&!
                    & 'rain_mm'//'; '//'Tmax'//'; '//&!
                    & 'Tmin'//'; '//'et0'//'; '// &!
                    & 'kcb'//'; '//'lai'//'; '//'pday'//'; '//&!
                    & 'irrig_mm'//'; '// &!
                    & 'eff_rain_mm'//'; '//'GrossAv_water_mm'//'; '//'NetAv_water_mm'//'; '// &!
                    & 'interception_mm'//'; '//'runoff_mm'//'; '//'infiltration_mm'//'; '// &!
                    & 'eva_pot_mm'//'; '//'eva_mm'//'; '// &!
                    & 'trasp_pot1_mm'//'; '//'trasp_act1_mm'//'; '// &!
                    & 'perc1_mm'//'; '//'theta1_mm'//'; '// &!
                    & 'ponding_mm'//'; '// 'rise_mm'//'; '// 'trasp_pot2_mm'//'; '// &!
                    & 'trasp_act2_mm'//'; '//'ks'//'; '// &!
                    & 'thickness_II_m'//'; '//'wat_table_depth_under_root_m'//'; '//'capflux_mm'//'; '// &!
                    & 'perc2_mm'//'; '//'theta2_mm'//'; '//'theta_old_mm'//'; '//&!
                    & 'rawbig'//'; '//'rawinf'//'; '// &!
                    & 'wat_table_depth_m'//'; '//&!
                    & 'distr_irr_mm'//'; '//'priv_well_irr_mm'//'; '//&!
                    & 'espperc1'//'; '//'espperc2'//'; '//'irr_loss_mm')
            end do!
                
            ! cell info
            do i=1,size(out_tbl_list%cell_info)!
                write(row_str,*)out_tbl_list%cell_info(i)%coord%row!
                write(col_str,*)out_tbl_list%cell_info(i)%coord%col!
                out_tbl_list%cell_info(i)%file%fn =trim(adjustl(yr))//'_cellinfo_'//trim(adjustl(row_str)) &!
                    & //'_'//trim(adjustl(col_str))//'.csv'!
                call init_cell_output_file(out_tbl_list%cell_info(i)%file%unit,trim(path)//trim(adjustl(out_tbl_list%cell_info(i)%file%fn)), 'input files')!
            end do!

            do i=1,size(out_tbl_list%prod_info)!
            ! cell parameter
                write(row_str,*)out_tbl_list%prod_info(i)%coord%row!
                write(col_str,*)out_tbl_list%prod_info(i)%coord%col!
                out_tbl_list%prod_info(i)%file%fn =trim(adjustl(yr))//'_cellparameters_'//trim(adjustl(row_str)) &!
                    & //'_'//trim(adjustl(col_str))//'.csv'!
                call init_cell_output_file(out_tbl_list%prod_info(i)%file%unit,trim(path)//trim(adjustl(out_tbl_list%prod_info(i)%file%fn)), 'input files')!
            end do!
            
            if (debug .eqv. .true.) then
                ! convergence log
                do i=1,size(out_tbl_list%cell_conv)!
                    write(row_str,*)out_tbl_list%cell_conv(i)%coord%row!
                    write(col_str,*)out_tbl_list%cell_conv(i)%coord%col!
                    out_tbl_list%cell_conv(i)%file%fn =trim(adjustl(yr))//'_convergence_'//trim(adjustl(row_str)) &!
                        & //'_'//trim(adjustl(col_str))//'.csv'!
                    call init_cell_output_file(out_tbl_list%cell_conv(i)%file%unit,trim(path)//trim(adjustl(out_tbl_list%cell_conv(i)%file%fn)), &!
                        & 'date; hour; mmax1; nIter1; mmax2; nIter2')!
                end do!
                ! evaporation terms
                do i=1,size(out_tbl_list%cell_eva)!
                    write(row_str,*)out_tbl_list%cell_eva(i)%coord%row!
                    write(col_str,*)out_tbl_list%cell_eva(i)%coord%col!
                    out_tbl_list%cell_eva(i)%file%fn = trim(adjustl(yr))//'_cellevaporation_'//trim(adjustl(row_str)) &!
                        & //'_'//trim(adjustl(col_str))//'.csv'!
                    call init_cell_output_file(out_tbl_list%cell_eva(i)%file%unit,trim(path)//trim(adjustl(out_tbl_list%cell_eva(i)%file%fn)), &!
                        & 'Giulian_day'//'; '//&!
                        & 'u2'//'; '//'RHmin'//'; '//&                                  ! Meteorological parameters
                        & 'ET0'//'; '//&
                        & 'kcb'//'; '//'h'//'; '// &                                    ! Phenological parameters
                        & 'fw'//'; '//'few'//'; '//'fc'//';'//'Kcmax'//'; '//&          ! Wetted soil fraction, cover fraction & Kcmax
                        & 'wat_REW_mm'//';'//'theta1_wp_mm'//';'//'theta1_mm'//'; '//&  ! 1st layer water content
                        & 'Ke'//'; '//&                                                 ! Evaporation coefficient
                        & 'eva_pot_mm'//'; '//'eva_mm'//';'//'kr'//';'//'fw_old')       ! Evaporation ! %RR% test kr fw_old
                end do!
                ! CN & runoff terms
                do i=1, size(out_tbl_list%cell_cn)!
                    write(row_str,*)out_tbl_list%cell_cn(i)%coord%row!
                    write(col_str,*)out_tbl_list%cell_cn(i)%coord%col!
                    out_tbl_list%cell_cn(i)%file%fn = trim(adjustl(yr))//'_cellrunoff_'//trim(adjustl(row_str)) &!
                        & //'_'//trim(adjustl(col_str))//'.csv'!
                    call init_cell_output_file(out_tbl_list%cell_cn(i)%file%unit,trim(path)//trim(adjustl(out_tbl_list%cell_cn(i)%file%fn)), &!
                        & 'Giulian_day'//'; '//&!
                        & 'teta_wp'//'; '//'teta_fc'//'; '//&                           ! Water content parameters
                        & 'teta_sat'//'; '//'teta_fcwp'//'; '//'teta_end_mm'//'; '//&
                        & 'CN2_tab'//'; '//'CN3_tab'//'; '//&                           ! CN tabulated values
                        & 'slope'//'; '//&                                              ! Slope
                        & 'pheno_type'//'; '//'pheno_phase'//'; '//&                    ! Crop parameters
                        & 'cn_baresoil'//'; '//'cn_slope'//'; '//&                      ! CN parameters
                        & 'cn1_pheno'//'; '//'cn2_pheno'//'; '//&                       ! CN parameters for each moisture condition
                        & 'cn3_pheno'//'; '//'cn_moist'//'; '//&
                        & 'lambda'//'; '//&                                             ! lambda parameters
                        & 'S'//'; '//'Ia'//'; '//&                                      ! Initial abstractions
                        & 'GrossAv'//'; '//'NetAv'//'; '//&                             ! Water availability
                        & 'GrossRunoff'//'; '//'NetRunoff')                             ! Runoff
                end do
            end if
        end if
        
        ! TODO: add aggregation also for need mode 
        if (mode == 1) then
            allocate(id_irr_unit_str(size(id_irr_unit_list)))
            ! daily series of irrigation water supply to each irrigation unit from StartIrrSeason to EndIrrSeason - Qirrunit = Qtot_fonte * proporzione_dafonte
            out_tbl_list%q_irr_units%fn =trim(adjustl(yr))//'_Qirrunits.csv'!
        
            title=""!
            do i=1,size(id_irr_unit_list)!
                write(id_irr_unit_str(i),*)id_irr_unit_list(i)!
                id_irr_unit_str(i)="SubDistr_"//trim(adjustl(id_irr_unit_str(i)))!
            end do!
            allocate(w_str(size(id_irr_unit_str)+2))!
            w_str(1)="DoY"!
            w_str(2)="Source"!
            w_str(3:size(w_str))=id_irr_unit_str(:)!
            call init_cell_output_file(out_tbl_list%q_irr_units%unit,trim(path)//trim(adjustl(out_tbl_list%q_irr_units%fn)),trim(title),w_str)!
            deallocate(w_str)!

            ! daily series of irrigation water supply from collective runtime sources
            if (n_nm_pub >=1) then
                out_tbl_list%q_un_col%fn =trim(adjustl(yr))//'_Qcrs.csv'!
                title=""!
                do i=1,size(id_irr_unit_list)!
                    write(id_irr_unit_str(i),*)id_irr_unit_list(i)!
                    id_irr_unit_str(i)="SubDistr_"//trim(adjustl(id_irr_unit_str(i)))!
                end do!
                allocate(w_str(size(id_irr_unit_str)+1)) !allocate(w_str(n_nm_pub+1))!
                w_str(1)="DoY"!
                w_str(2:size(w_str))=id_irr_unit_str(:) !w_str(2:size(w_str))=id_nm_pub_list(:)!
                call init_cell_output_file(out_tbl_list%q_un_col%unit,trim(path)//trim(adjustl(out_tbl_list%q_un_col%fn)),trim(title),w_str)!
                deallocate(w_str)!
            end if

            ! daily series of water used in each irrigation units
            out_tbl_list%q_irr%fn =trim(adjustl(yr))//'_Qirr.csv'!
            title=""!
            do i=1,size(id_irr_unit_list)!
                write(id_irr_unit_str(i),*)id_irr_unit_list(i)!
                id_irr_unit_str(i)="SubDistr_"//trim(adjustl(id_irr_unit_str(i)))!
            end do!
            allocate(w_str(size(id_irr_unit_str)+1))!
            w_str(1)="DoY"!
            w_str(2:size(w_str))=id_irr_unit_str(:)!
            call init_cell_output_file(out_tbl_list%q_irr%unit,trim(path)//trim(adjustl(out_tbl_list%q_irr%fn)),trim(title),w_str)!
            deallocate(w_str)!
            
            ! daily series of water not used for irrigation and released by each irrigation units
            out_tbl_list%q_surplus%fn =trim(adjustl(yr))//'_Qsurplus.csv'!
            title=""!
            do i=1,size(id_irr_unit_list)!
                write(id_irr_unit_str(i),*)id_irr_unit_list(i)!
                id_irr_unit_str(i)="SubDistr_"//trim(adjustl(id_irr_unit_str(i)))!
            end do!
            allocate(w_str(size(id_irr_unit_str)+1))!
            w_str(1)="DoY"!
            w_str(2:size(w_str))=id_irr_unit_str(:)!
            call init_cell_output_file(out_tbl_list%q_surplus%unit,trim(path)//trim(adjustl(out_tbl_list%q_surplus%fn)),trim(title),w_str)!
            deallocate(w_str)!
            
            ! daily series of water not used for irrigation and available the following day
            out_tbl_list%q_rem%fn =trim(adjustl(yr))//'_Qrem.csv'!
            title=""!
            do i=1,size(id_irr_unit_list)!
                write(id_irr_unit_str(i),*)id_irr_unit_list(i)!
                id_irr_unit_str(i)="SubDistr_"//trim(adjustl(id_irr_unit_str(i)))!
            end do!
            allocate(w_str(size(id_irr_unit_str)+1))!
            w_str(1)="DoY"!
            w_str(2:size(w_str))=id_irr_unit_str(:)!
            call init_cell_output_file(out_tbl_list%q_rem%unit,trim(path)//trim(adjustl(out_tbl_list%q_rem%fn)),trim(title),w_str)!
            deallocate(w_str)!        
            
            ! daily series of water withdrawled by private sources for irrigation
            out_tbl_list%q_un_priv%fn =trim(adjustl(yr))//'_Qprivate.csv'!
            title=""!
            do i=1,size(id_irr_unit_list)!
                write(id_irr_unit_str(i),*)id_irr_unit_list(i)!
                id_irr_unit_str(i)="SubDistr_"//trim(adjustl(id_irr_unit_str(i)))!
            end do!
            allocate(w_str(size(id_irr_unit_str)+1))!
            w_str(1)="DoY"!
            w_str(2:size(w_str))=id_irr_unit_str(:)!
            call init_cell_output_file(out_tbl_list%q_un_priv%unit,trim(path)//trim(adjustl(out_tbl_list%q_un_priv%fn)),trim(title),w_str)!
            deallocate(w_str)!		
	
            ! water shift
            out_tbl_list%n_inv_cells%fn =trim(adjustl(yr))//'_Watshift.csv'!
            title=""!
            do i=1,size(id_irr_unit_list)!
                write(id_irr_unit_str(i),*)id_irr_unit_list(i)!
                id_irr_unit_str(i)="SubDistr_"//trim(adjustl(id_irr_unit_str(i)))!
            end do!
            allocate(w_str(size(id_irr_unit_str)+1))!
            w_str(1)="DoY"!
            w_str(2:size(w_str))=id_irr_unit_str(:)!
            call init_cell_output_file(out_tbl_list%n_inv_cells%unit,trim(path)//trim(adjustl(out_tbl_list%n_inv_cells%fn)),trim(title),w_str)!
            deallocate(w_str)!
        end if
        
        ! TODO: not referred to cells, move away
        if (debug .eqv. .true.) then
            allocate(id_ws_str(size(id_ws_list)))
            ! evapotranspiration for each weather station
            out_tbl_list%et0_ws%fn = trim(adjustl(yr))//'_et0_stations.csv'!
            title=""!
            do i=1,size(id_ws_list)!
                id_ws_str(i)=trim(adjustl(id_ws_list(i)))!
            end do!
            allocate(w_str(size(id_ws_str)+1))!
            w_str(1)="DoY"!
            w_str(2:size(w_str))=id_ws_str(:)!
            call init_cell_output_file(out_tbl_list%et0_ws%unit,&
                               & trim(path)//trim(adjustl(out_tbl_list%et0_ws%fn)),&
                               & trim(title),w_str)!
            deallocate(w_str)!
        end if
        !
    end subroutine init_cell_output_by_year!
    !
    subroutine close_cell_output_by_year(out_tbl_list,mode,f_cell_exists,debug)!
        ! close all opened file
        integer,intent(in)::mode
        logical,intent(in)::f_cell_exists,debug
        type(output_table_list),intent(inout)::out_tbl_list
        integer::i!
        !!
        if (f_cell_exists .eqv. .true.) then
            do i=1,size(out_tbl_list%sample_cells)!
                close(out_tbl_list%sample_cells(i)%file%unit)!
                close(out_tbl_list%cell_info(i)%file%unit)!
                close(out_tbl_list%prod_info(i)%file%unit)!
                if (debug .eqv. .true.) then
                    close(out_tbl_list%cell_conv(i)%file%unit)!
                    close(out_tbl_list%cell_eva(i)%file%unit)!
                    close(out_tbl_list%cell_cn(i)%file%unit)!
                end if
            end do!
            deallocate(out_tbl_list%sample_cells,out_tbl_list%cell_info,out_tbl_list%prod_info)
            if (debug .eqv. .true.) then
                deallocate(out_tbl_list%cell_conv)   !to re-allocate them the year after!
                deallocate(out_tbl_list%cell_eva)
                deallocate(out_tbl_list%cell_cn)
            end if
        end if
        if (mode == 1 ) then
            close(out_tbl_list%q_irr_units%unit)!
            close(out_tbl_list%q_un_col%unit)!
            close(out_tbl_list%q_surplus%unit)!
            close(out_tbl_list%q_irr%unit)!
            close(out_tbl_list%q_rem%unit)!
            close(out_tbl_list%q_un_priv%unit)!
            close(out_tbl_list%n_inv_cells%unit)!
        end if
        if (debug .eqv. .true.) then
            close(out_tbl_list%et0_ws%unit)!
        end if
    end subroutine close_cell_output_by_year!
    !
    subroutine write_cell_info(info_spat,cell_info,mode,f_caprise)!
        ! write parameters of each cells
        type(spatial_info),intent(in)::info_spat!
        type(cell_output),dimension(:),intent(in)::cell_info!
        integer, intent(in)::mode
        logical, intent(in)::f_caprise           ! if true, capillary rise is activated
        integer::x,y!
        integer::i!
        integer::k
        !!
        do i=1,size(cell_info)!
            x=cell_info(i)%coord%row!
            y=cell_info(i)%coord%col!
            if (x > size(info_spat%domain%mat,1)) then
                print *, "Row index of cell ", x, ", ",y, " is higher than domain's number of rows ", size(info_spat%domain%mat,1)
                print *, 'Execution will be aborted...'
                stop
            end if
            if (y > size(info_spat%domain%mat,2)) then
                print *, "Column index of cell ", x, ", ", y, " is higher than domain's number of column ", &
                    & size(info_spat%domain%mat,2)
                print *, 'Execution will be aborted...'
                stop
            end if
            write(cell_info(i)%file%unit,*)'Domain; ', info_spat%domain%mat(x,y)!
            write(cell_info(i)%file%unit,*)'Slope; ', info_spat%slope%mat(x,y)!
            write(cell_info(i)%file%unit,*)'Soil use; ', info_spat%soil_use_id%mat(x,y)!
            select case (mode)
                case (0,2:3)
                    write(cell_info(i)%file%unit,*)'Irrigation subdistrict; null'
                case (1,4)
                    write(cell_info(i)%file%unit,*)'Irrigation subdistrict; ', info_spat%irr_unit_id%mat(x,y)!
                case default
            end select
            select case (mode)
                case (0)
                    write(cell_info(i)%file%unit,*)'Irrigation method; null'
                case (1:4)
                    write(cell_info(i)%file%unit,*)'Irrigation method; ', info_spat%irr_meth_id%mat(x,y)!
                case default
            end select
            select case (mode)
                case (0,1,3)
                    write(cell_info(i)%file%unit,*)'Method efficiency; null'
                case (2,4)
                    write(cell_info(i)%file%unit,*)'Method efficiency; ', info_spat%eff_met%mat(x,y)!
                case default
            end select
            select case (mode)
                case (0, 2, 3)
                    write(cell_info(i)%file%unit,*)'Conveyance efficiency; null' ! %RR%
                case (1, 4)
                    write(cell_info(i)%file%unit,*)'Conveyance efficiency; ', info_spat%eff_net%mat(x,y)! %RR%
                case default
            end select
            write(cell_info(i)%file%unit,*)'Hydrologic soil group; ', info_spat%hydr_gr%mat(x,y)!
            write(cell_info(i)%file%unit,*)'Hydrologic condition; ', info_spat%drainage%mat(x,y)!
            write(cell_info(i)%file%unit,*)'ThetaI_sat; ', info_spat%theta(1)%sat%mat(x,y)!
            write(cell_info(i)%file%unit,*)'ThetaI_fc; ', info_spat%theta(1)%FC%mat(x,y)!
            write(cell_info(i)%file%unit,*)'ThetaI_wp; ', info_spat%theta(1)%WP%mat(x,y)!
            write(cell_info(i)%file%unit,*)'ThetaI_r; ', info_spat%theta(1)%r%mat(x,y)!
            write(cell_info(i)%file%unit,*)'ThetaII_sat; ', info_spat%theta(2)%sat%mat(x,y)!
            write(cell_info(i)%file%unit,*)'ThetaII_fc; ', info_spat%theta(2)%FC%mat(x,y)!
            write(cell_info(i)%file%unit,*)'ThetaII_wp; ', info_spat%theta(2)%WP%mat(x,y)!
            write(cell_info(i)%file%unit,*)'ThetaII_r; ', info_spat%theta(2)%r%mat(x,y)!
            write(cell_info(i)%file%unit,*)'ksat_I; ', info_spat%k_sat(1)%mat(x,y)!
            write(cell_info(i)%file%unit,*)'ksat_II; ', info_spat%k_sat(2)%mat(x,y)!
            write(cell_info(i)%file%unit,*)'expn_I; ', info_spat%fact_n(1)%mat(x,y)!
            write(cell_info(i)%file%unit,*)'expn_II; ', info_spat%fact_n(2)%mat(x,y)!
            if (f_caprise .eqv. .true.) then
                write(cell_info(i)%file%unit,*)'CapFluxParam_a3; ', info_spat%a3%mat(x,y)!
                write(cell_info(i)%file%unit,*)'CapFluxParam_a4; ', info_spat%a4%mat(x,y)!
                write(cell_info(i)%file%unit,*)'CapFluxParam_b1; ', info_spat%b1%mat(x,y)!
                write(cell_info(i)%file%unit,*)'CapFluxParam_b2; ', info_spat%b2%mat(x,y)!
                write(cell_info(i)%file%unit,*)'CapFluxParam_b3; ', info_spat%b3%mat(x,y)!
                write(cell_info(i)%file%unit,*)'CapFluxParam_b4; ', info_spat%b4%mat(x,y)!
                write(cell_info(i)%file%unit,*)'Water depth; ', info_spat%wat_tab%mat(x,y)!
            else
                write(cell_info(i)%file%unit,*)'CapFluxParam_a3; null'
                write(cell_info(i)%file%unit,*)'CapFluxParam_a4; null'
                write(cell_info(i)%file%unit,*)'CapFluxParam_b1; null'
                write(cell_info(i)%file%unit,*)'CapFluxParam_b2; null'
                write(cell_info(i)%file%unit,*)'CapFluxParam_b3; null'
                write(cell_info(i)%file%unit,*)'CapFluxParam_b4; null'
                write(cell_info(i)%file%unit,*)'Water depth; null'
            end if
            !!
            select case (mode)
                case (1:4)
                    write(cell_info(i)%file%unit,*)'am_perc1; ', info_spat%a_perc(1)%mat(x,y)!
                    write(cell_info(i)%file%unit,*)'am_perc2; ', info_spat%a_perc(2)%mat(x,y)!
                    write(cell_info(i)%file%unit,*)'bm_perc1; ', info_spat%b_perc(1)%mat(x,y)!
                    write(cell_info(i)%file%unit,*)'bm_perc2; ', info_spat%b_perc(2)%mat(x,y)!
                case default
            end select
            !
            do k=1,size(info_spat%weight_ws)!
                write(cell_info(i)%file%unit,*) 'Meteorological Station nr', k ,'; ',info_spat%weight_ws(k)%mat(x,y)
            end do
        end do!
    end subroutine write_cell_info!
    
    subroutine write_cell_prod(info_prod, crop, irandom)!
        ! write parameters related to productivity of each cells
        type(cell_output),dimension(:),intent(in)::info_prod!
        type(crop_matrices), intent(in)::crop
        integer, dimension(:,:), intent(in)::irandom
        integer::x,y,z!
        integer::zmax
        integer::i!
        integer::nan=-9999.
        !!
        do i=1,size(info_prod)!
            x=info_prod(i)%coord%row!
            y=info_prod(i)%coord%col!
            zmax = 1
            do z=1,size(crop%wp_adj,3)
                if (crop%wp_adj(x,y,z) /= nan) then
                    zmax = z
                end if
            end do
            write(info_prod(i)%file%unit,*)'WPadj; ', (crop%wp_adj(x,y,z), '; ', z=1,zmax)
            write(info_prod(i)%file%unit,*)'HI; ', (crop%HI(x,y,z), '; ', z=1,zmax)
            write(info_prod(i)%file%unit,*)'KyT; ', (crop%Ky_tot(x,y,z), '; ', z=1,zmax)
            write(info_prod(i)%file%unit,*)'Ky1; ', (crop%Ky_pheno(x,y,z,1), '; ', z=1,zmax)
            write(info_prod(i)%file%unit,*)'Ky2; ', (crop%Ky_pheno(x,y,z,2), '; ', z=1,zmax)
            write(info_prod(i)%file%unit,*)'Ky3; ', (crop%Ky_pheno(x,y,z,3), '; ', z=1,zmax)
            write(info_prod(i)%file%unit,*)'Ky4; ', (crop%Ky_pheno(x,y,z,4), '; ', z=1,zmax)
            write(info_prod(i)%file%unit,*)'Tcrit; ', (crop%T_crit(x,y,z), '; ', z=1,zmax)
            write(info_prod(i)%file%unit,*)'Tlim; ', (crop%T_lim(x,y,z), '; ', z=1,zmax)
            write(info_prod(i)%file%unit,*)'Kcbmin; ', (crop%k_cb_min(x,y,z), '; ', z=1,zmax)
            write(info_prod(i)%file%unit,*)'Kcbini; ', (crop%k_cb_mid(x,y,z), '; ', z=1,zmax)
            write(info_prod(i)%file%unit,*)'Kcbmax; ', (crop%k_cb_max(x,y,z), '; ', z=1,zmax)
            write(info_prod(i)%file%unit,*)'ii0; ', (crop%ii0(x,y,z), '; ', z=1,zmax)
            write(info_prod(i)%file%unit,*)'iie; ', (crop%iie(x,y,z), '; ', z=1,zmax)
            write(info_prod(i)%file%unit,*)'iid; ', (crop%iid(x,y,z), '; ', z=1,zmax)
            write(info_prod(i)%file%unit,*)'dij; ', (crop%dij(x,y,z), '; ', z=1,zmax)
            write(info_prod(i)%file%unit,*)'irandom; ', irandom(x,y)!
        end do!
    end subroutine write_cell_prod!
    
    subroutine init_step_output(a_step_map,domain)!
        type(step_map)::a_step_map!
        integer,dimension(:,:),intent(in)::domain!
        !!
        allocate(a_step_map%rain%mat(size(domain,1),size(domain,2)))!
        allocate(a_step_map%transp_act%mat(size(domain,1),size(domain,2)))!
        allocate(a_step_map%transp_pot%mat(size(domain,1),size(domain,2)))!
        allocate(a_step_map%irr%mat(size(domain,1),size(domain,2)))!
        allocate(a_step_map%irr_loss%mat(size(domain,1),size(domain,2)))!
        allocate(a_step_map%cap_rise%mat(size(domain,1),size(domain,2)))!
        allocate(a_step_map%irr_nm_priv%mat(size(domain,1),size(domain,2)))!
        allocate(a_step_map%irr_nm_col%mat(size(domain,1),size(domain,2)))!
        allocate(a_step_map%deep_perc%mat(size(domain,1),size(domain,2)))!
        allocate(a_step_map%runoff%mat(size(domain,1),size(domain,2)))
        allocate(a_step_map%et_pot%mat(size(domain,1),size(domain,2)))
        allocate(a_step_map%et_act%mat(size(domain,1),size(domain,2)))
    end subroutine init_step_output!
    !
    subroutine init_step_debug_output(step_dbg_map,domain)!
        type(step_debug_map)::step_dbg_map!
        integer,dimension(:,:),intent(in)::domain!
        !!
        allocate(step_dbg_map%eva_act%mat(size(domain,1),size(domain,2)))!
        allocate(step_dbg_map%eff_rain%mat(size(domain,1),size(domain,2)))!
        allocate(step_dbg_map%perc1%mat(size(domain,1),size(domain,2)))!
        allocate(step_dbg_map%perc2%mat(size(domain,1),size(domain,2)))!
        allocate(step_dbg_map%h_soil1%mat(SIZE(domain,1),size(domain,2)))
        allocate(step_dbg_map%h_soil2%mat(SIZE(domain,1),size(domain,2)))
    end subroutine init_step_debug_output!
    !
    subroutine init_yearly_output(yr_map,domain)!
        type(annual_map)::yr_map!
        integer,dimension(:,:),intent(in)::domain!
        !integer,parameter::n_kcb_stages=4 ! not used
        
        allocate(yr_map%rain%mat(size(domain,1),size(domain,2)))!
        allocate(yr_map%rain_crop_season%mat(size(domain,1),size(domain,2)))!
        allocate(yr_map%irr%mat(size(domain,1),size(domain,2)))!
        allocate(yr_map%irr_loss%mat(size(domain,1),size(domain,2)))!
        allocate(yr_map%eva_pot_crop_season%mat(size(domain,1),size(domain,2)))!
        allocate(yr_map%eva_act_crop_season%mat(size(domain,1),size(domain,2)))!
        allocate(yr_map%transp_act%mat(size(domain,1),size(domain,2)))!
        allocate(yr_map%transp_pot%mat(size(domain,1),size(domain,2)))!
        allocate(yr_map%runoff%mat(size(domain,1),size(domain,2)))!
        allocate(yr_map%net_flux_gw%mat(size(domain,1),size(domain,2)))!
        allocate(yr_map%total_eff%mat(size(domain,1),size(domain,2)))!
        allocate(yr_map%n_irr_events%mat(size(domain,1),size(domain,2)))!
        allocate(yr_map%h_irr_mean%mat(size(domain,1),size(domain,2)))!
    end subroutine init_yearly_output!
    
    subroutine init_yearly_yield_output(yield_map,domain,cs)!
        type(yield_t)::yield_map!
        integer,dimension(:,:),intent(in)::domain!
        integer,intent(in)::cs
        integer,parameter::fasi_kcb=4
        !!
        allocate(yield_map%biomass_pot%mat(size(domain,1),size(domain,2),cs))!
        allocate(yield_map%yield_pot%mat(size(domain,1),size(domain,2),cs))!
        allocate(yield_map%yield_act%mat(size(domain,1),size(domain,2),cs))!
        allocate(yield_map%f_WS%mat(size(domain,1),size(domain,2),cs))!
        allocate(yield_map%f_WS_stage%mat(size(domain,1),size(domain,2),cs))!
        allocate(yield_map%f_HS%mat(size(domain,1),size(domain,2),cs))!
        allocate(yield_map%f_HS_sum%mat(size(domain,1),size(domain,2),cs))!
        allocate(yield_map%transp_ratio_sum%mat(size(domain,1),size(domain,2),cs))!
        allocate(yield_map%T_act_sum%mat(size(domain,1),size(domain,2),fasi_kcb, cs))!
        allocate(yield_map%T_pot_sum%mat(size(domain,1),size(domain,2),fasi_kcb, cs))!
        allocate(yield_map%dev_stage%mat(size(domain,1),size(domain,2),fasi_kcb, cs))!
    end subroutine init_yearly_yield_output!
    !
    subroutine init_yearly_debug_output(dbg_yr_map,domain)!
        type(annual_debug_map)::dbg_yr_map!
        integer,dimension(:,:),intent(in)::domain!
        !!
        allocate(dbg_yr_map%eva_act_tot%mat(size(domain,1),size(domain,2)))!
        allocate(dbg_yr_map%rain_eff%mat(size(domain,1),size(domain,2)))!
        allocate(dbg_yr_map%iter1%mat(size(domain,1),size(domain,2)))!
        allocate(dbg_yr_map%iter2%mat(size(domain,1),size(domain,2)))!
    end subroutine init_yearly_debug_output!
    !
    subroutine destroy_step_output(a_step_map)!
        type(step_map)::a_step_map!
        !!
        deallocate(a_step_map%rain%mat)!
        deallocate(a_step_map%transp_act%mat)!
        deallocate(a_step_map%transp_pot%mat)!
        deallocate(a_step_map%irr%mat)!
        deallocate(a_step_map%irr_loss%mat)!
        deallocate(a_step_map%cap_rise%mat)!
        deallocate(a_step_map%irr_nm_priv%mat)!
        deallocate(a_step_map%irr_nm_col%mat)!
        deallocate(a_step_map%deep_perc%mat)!
        deallocate(a_step_map%runoff%mat)
        deallocate(a_step_map%et_pot%mat)
        deallocate(a_step_map%et_act%mat)
    end subroutine destroy_step_output!
    !
    subroutine destroy_step_debug_output(a_step_dbg_map)!
        type(step_debug_map)::a_step_dbg_map!
        !!
        deallocate(a_step_dbg_map%eva_act%mat)!
        deallocate(a_step_dbg_map%perc1%mat)!
        deallocate(a_step_dbg_map%eff_rain%mat)!
        deallocate(a_step_dbg_map%perc2%mat)!
        deallocate(a_step_dbg_map%h_soil1%mat)
        deallocate(a_step_dbg_map%h_soil2%mat)
    end subroutine destroy_step_debug_output!
    !
    subroutine destroy_annual_output(yr_map)!
        type(annual_map)::yr_map!
        deallocate(yr_map%rain%mat)!
        deallocate(yr_map%rain_crop_season%mat)!
        deallocate(yr_map%irr%mat)!
        deallocate(yr_map%irr_loss%mat)!
        deallocate(yr_map%eva_pot_crop_season%mat)!
        deallocate(yr_map%eva_act_crop_season%mat)!
        deallocate(yr_map%transp_act%mat)!
        deallocate(yr_map%transp_pot%mat)!
        deallocate(yr_map%runoff%mat)!
        deallocate(yr_map%net_flux_gw%mat)!
        deallocate(yr_map%total_eff%mat)!
        deallocate(yr_map%n_irr_events%mat)!
        deallocate(yr_map%h_irr_mean%mat)!
    end subroutine destroy_annual_output!
    
    subroutine destroy_yield_output(yld_map)!
        implicit none!
        type(yield_t)::yld_map!
        !!
        deallocate(yld_map%biomass_pot%mat)!
        deallocate(yld_map%yield_pot%mat)!
        deallocate(yld_map%yield_act%mat)!
        deallocate(yld_map%f_WS%mat)!
        deallocate(yld_map%f_WS_stage%mat)!
        deallocate(yld_map%f_HS%mat)!
        deallocate(yld_map%f_HS_sum%mat)!
        deallocate(yld_map%transp_ratio_sum%mat)!
        deallocate(yld_map%T_act_sum%mat)!
        deallocate(yld_map%T_pot_sum%mat)!
        deallocate(yld_map%dev_stage%mat)!
    end subroutine destroy_yield_output!
    !
    subroutine destroy_annual_debug_output(dbg_yr_map)!
        implicit none!
        type(annual_debug_map)::dbg_yr_map!
        !!
        deallocate(dbg_yr_map%eva_act_tot%mat)!
        deallocate(dbg_yr_map%rain_eff%mat)!
        deallocate(dbg_yr_map%iter1%mat)!
        deallocate(dbg_yr_map%iter2%mat)!
    end subroutine destroy_annual_debug_output!
    
    ! TODO: update name of the output file
    subroutine init_step_output_file(a_step_map,path,yr,doy,calendar,init_total, string)!
        character(len=*),intent(in)::path!
        integer,intent(in)::doy!
        character(len=*),intent(in)::yr
        integer,dimension(:),intent(in)::calendar!
        integer,intent(in)::init_total
        character(len=*),intent(in)::string
        type(step_map),intent(out)::a_step_map!
        integer::i,total!
        character(len=33)::year,step!
        !!
        total=init_total!

        do i=1,size(calendar)!
            total = total+calendar(i)!
            if(doy==total-calendar(i)+1)then   ! day of the month/period
                write(step,*)i!
                step=trim(adjustl(string))//trim(adjustl(step))//'_'!
                year=trim(adjustl(yr))//'_'!
                a_step_map%rain%fn=trim(adjustl(trim(path)//trim(adjustl(year))//trim(adjustl(step))//'prec.asc'))!
                a_step_map%transp_act%fn=trim(adjustl(trim(path)//trim(adjustl(year))//trim(adjustl(step))//'trasp_act.asc'))!
                a_step_map%transp_pot%fn=trim(adjustl(trim(path)//trim(adjustl(year))//trim(adjustl(step))//'trasp_pot.asc'))!
                a_step_map%irr%fn=trim(adjustl(trim(path)//trim(adjustl(year))//trim(adjustl(step))//'irr.asc'))!
                a_step_map%irr_loss%fn=trim(adjustl(trim(path)//trim(adjustl(year))//trim(adjustl(step))//'irr_loss.asc'))!
                a_step_map%cap_rise%fn=trim(adjustl(trim(path)//trim(adjustl(year))//trim(adjustl(step))//'caprise.asc'))!
                a_step_map%irr_nm_priv%fn=trim(adjustl(trim(path)//trim(adjustl(year))//trim(adjustl(step))//'irr_privw.asc'))!
                a_step_map%irr_nm_col%fn=trim(adjustl(trim(path)//trim(adjustl(year))//trim(adjustl(step))//'irr_distr.asc'))!
                a_step_map%deep_perc%fn=trim(adjustl(trim(path)//trim(adjustl(year))//trim(adjustl(step))//'flux2.asc'))!
                a_step_map%runoff%fn=trim(adjustl(trim(path)//trim(adjustl(year))//trim(adjustl(step))//'runoff.asc'))!
                a_step_map%et_pot%fn=trim(adjustl(trim(path)//trim(adjustl(year))//trim(adjustl(step))//'et_pot.asc'))!
                a_step_map%et_act%fn=trim(adjustl(trim(path)//trim(adjustl(year))//trim(adjustl(step))//'et_act.asc'))!
                a_step_map = 0.0D0    ! init to zero
                exit!
            else if (doy < total-calendar(i)+1) then ! exit from the cycle
                exit
            else
                cycle!
            end if!
        end do!
    end subroutine init_step_output_file!
    
    
    subroutine init_step_debug_output_file(a_dbg_map,path,yr,doy,calendar,init_total,string)!
        character(len=*),intent(in)::path!
        integer,intent(in)::doy!
        character(len=*),intent(in)::yr
        integer,dimension(:),intent(in)::calendar!
        integer,intent(in)::init_total
        character(len=*),intent(in)::string
        type(step_debug_map),intent(out)::a_dbg_map!
        integer::i,total!
        character(len=33)::year,step!
        
        total=init_total!

        do i=1,size(calendar)!
            total = total+calendar(i)!
            if(doy==total-calendar(i)+1)then   ! first day of the step
                write(step,*)i!
                step=trim(adjustl(string))//trim(adjustl(step))//'_'!
                year=trim(adjustl(yr))//'_'!
                a_dbg_map%eva_act%fn=trim(adjustl(trim(path)//trim(adjustl(year))//trim(adjustl(step))//'eva.asc'))!
                a_dbg_map%eff_rain%fn=trim(adjustl(trim(path)//trim(adjustl(year))//trim(adjustl(step))//'prec_eff.asc'))!
                a_dbg_map%perc1%fn=trim(adjustl(trim(path)//trim(adjustl(year))//trim(adjustl(step))//'perc1.asc'))!
                a_dbg_map%perc2%fn=trim(adjustl(trim(path)//trim(adjustl(year))//trim(adjustl(step))//'perc2.asc'))!
                a_dbg_map%h_soil1%fn=trim(adjustl(trim(path)//trim(adjustl(year))//trim(adjustl(step))//'theta1.asc'))!
                a_dbg_map%h_soil2%fn=trim(adjustl(trim(path)//trim(adjustl(year))//trim(adjustl(step))//'theta2.asc'))!
                a_dbg_map = 0.0D0 ! set to zero after
                exit!
            else if (doy < total-calendar(i)+1) then ! exit if lower
                exit
            else!
                cycle!
            end if!
        end do!
    end subroutine init_step_debug_output_file!
    !
    subroutine init_yearly_output_file(a_yr_map,path,year)
        character(len=*),intent(in)::path!
        character(len=*),intent(in)::year
        type(annual_map),intent(out)::a_yr_map!
        character(len=33)::year_str
        !!
        year_str=trim(adjustl(year))//'_'!
        a_yr_map%rain%fn=trim(adjustl(trim(path)//trim(adjustl(year_str))//'prec_tot.asc'))
        a_yr_map%rain_crop_season%fn=trim(adjustl(trim(path)//trim(adjustl(year_str))//'prec_agr.asc'))
        a_yr_map%irr%fn=trim(adjustl(trim(path)//trim(adjustl(year_str))//'irr_tot.asc'))
        a_yr_map%irr_loss%fn=trim(adjustl(trim(path)//trim(adjustl(year_str))//'irr_loss.asc'))
        a_yr_map%eva_act_crop_season%fn=trim(adjustl(trim(path)//trim(adjustl(year_str))//'eva_act_agr.asc'))
        a_yr_map%eva_pot_crop_season%fn=trim(adjustl(trim(path)//trim(adjustl(year_str))//'eva_pot_agr.asc'))
        a_yr_map%transp_act%fn=trim(adjustl(trim(path)//trim(adjustl(year_str))//'trasp_act_tot.asc'))
        a_yr_map%transp_pot%fn=trim(adjustl(trim(path)//trim(adjustl(year_str))//'trasp_pot_tot.asc'))
        a_yr_map%runoff%fn=trim(adjustl(trim(path)//trim(adjustl(year_str))//'run_tot.asc'))
        a_yr_map%net_flux_gw%fn=trim(adjustl(trim(path)//trim(adjustl(year_str))//'net_flux_gw.asc'))
        a_yr_map%total_eff%fn=trim(adjustl(trim(path)//trim(adjustl(year_str))//'eff_tot.asc'))
        a_yr_map%n_irr_events%fn=trim(adjustl(trim(path)//trim(adjustl(year_str))//'irr_nr.asc'))
        a_yr_map%h_irr_mean%fn=trim(adjustl(trim(path)//trim(adjustl(year_str))//'irr_mean.asc'))
        a_yr_map = 0.0D0    ! init to zero
    end subroutine init_yearly_output_file!
    !
    subroutine init_yield_output_file(yield,path,year)
        character(len=*),intent(in)::path!
        character(len=*),intent(in)::year
        type(yield_t),intent(out)::yield!
        character(len=33)::year_str
        !!
        year_str=trim(adjustl(year))//'_'!
        yield%biomass_pot%fn=trim(adjustl(trim(path)//trim(adjustl(year_str))//'biomass_pot'))
        yield%yield_pot%fn=trim(adjustl(trim(path)//trim(adjustl(year_str))//'yield_pot'))
        yield%yield_act%fn=trim(adjustl(trim(path)//trim(adjustl(year_str))//'yield_act'))
        yield%T_act_sum%fn=trim(adjustl(trim(path)//trim(adjustl(year_str))//'T_act_sum'))
        yield%T_pot_sum%fn=trim(adjustl(trim(path)//trim(adjustl(year_str))//'T_pot_sum'))
        yield%f_WS_stage%fn=trim(adjustl(trim(path)//trim(adjustl(year_str))//'fcCS'))
        yield%f_WS%fn=trim(adjustl(trim(path)//trim(adjustl(year_str))//'fcT'))
        yield%f_HS%fn=trim(adjustl(trim(path)//trim(adjustl(year_str))//'fHS'))
        yield%f_HS_sum%fn=trim(adjustl(trim(path)//trim(adjustl(year_str))//'fHS_sum'))
        yield = 0.0D0    ! set to zero after
    end subroutine init_yield_output_file!
    !
    subroutine init_debug_yearly_output_file(a_yr_dbg_map,path,year)
        character(len=*),intent(in)::path!
        character(len=*),intent(in)::year
        type(annual_debug_map),intent(out)::a_yr_dbg_map!
        character(len=33)::year_str
        
        year_str=trim(adjustl(year))//'_'!
        a_yr_dbg_map%eva_act_tot%fn=trim(adjustl(trim(path)//trim(adjustl(year_str))//'eva_tot.asc'))
        a_yr_dbg_map%rain_eff%fn=trim(adjustl(trim(path)//trim(adjustl(year_str))//'eff_prec_tot.asc'))
        a_yr_dbg_map%iter1%fn=trim(adjustl(trim(path)//trim(adjustl(year_str))//'iter1.asc'))
        a_yr_dbg_map%iter2%fn=trim(adjustl(trim(path)//trim(adjustl(year_str))//'iter2.asc'))
        a_yr_dbg_map = 0.0D0    ! set to zero
    end subroutine init_debug_yearly_output_file!
    !
    subroutine save_step_data(a_step_map,doy,domain,calendar, init_total)
        ! save the results aggregated by step if the current day (doy) is the last day of the step
        type(step_map),intent(in)::a_step_map!
        integer,intent(in)::doy!
        type(grid_i),intent(in)::domain!
        integer,dimension(:),intent(in)::calendar!
        integer,intent(in)::init_total
        integer::i,total!
        integer::errorflag!
        !!
        errorflag=0!
        total=init_total!
        !!
        do i=1,size(calendar)!
            total=total+calendar(i)!
            if(doy==total)then    ! last day of the step
                where(domain%mat==domain%header%nan)!
                    a_step_map%rain%mat=real(domain%header%nan)!
                    a_step_map%transp_act%mat=real(domain%header%nan)!
                    a_step_map%transp_pot%mat=real(domain%header%nan)!
                    a_step_map%irr%mat=real(domain%header%nan)!
                    a_step_map%irr_loss%mat=real(domain%header%nan)!
                    a_step_map%cap_rise%mat=real(domain%header%nan)!
                    a_step_map%irr_nm_priv%mat=real(domain%header%nan)!
                    a_step_map%irr_nm_col%mat=real(domain%header%nan)!
                    a_step_map%deep_perc%mat=real(domain%header%nan)!
                    a_step_map%runoff%mat=real(domain%header%nan)
                    a_step_map%et_pot%mat=real(domain%header%nan)
                    a_step_map%et_act%mat=real(domain%header%nan)
                end where!
                call print_mat_as_grid(trim(a_step_map%rain%fn),domain%header,a_step_map%rain%mat,errorflag)!
                print*,"print file: ", trim(trim(a_step_map%rain%fn))!
                call print_mat_as_grid(trim(a_step_map%transp_act%fn),domain%header,a_step_map%transp_act%mat,errorflag)!
                print*,"print file: ", trim(trim(a_step_map%transp_act%fn))!
                call print_mat_as_grid(trim(a_step_map%transp_pot%fn),domain%header,a_step_map%transp_pot%mat,errorflag)!
                print*,"print file: ", trim(trim(a_step_map%transp_pot%fn))!
                call print_mat_as_grid(trim(a_step_map%irr%fn),domain%header,a_step_map%irr%mat,errorflag)!
                print*,"print file: ", trim(trim(a_step_map%irr%fn))!
                call print_mat_as_grid(trim(a_step_map%irr_loss%fn),domain%header,a_step_map%irr_loss%mat,errorflag)!
                print*,"print file: ", trim(trim(a_step_map%irr_loss%fn))!
                call print_mat_as_grid(trim(a_step_map%cap_rise%fn),domain%header,a_step_map%cap_rise%mat,errorflag)!
                print*,"print file: ", trim(trim(a_step_map%cap_rise%fn))!
                call print_mat_as_grid(trim(a_step_map%irr_nm_priv%fn),domain%header,a_step_map%irr_nm_priv%mat,errorflag)!
                print*,"print file: ", trim(trim(a_step_map%irr_nm_priv%fn))!
                call print_mat_as_grid(trim(a_step_map%irr_nm_col%fn),domain%header,a_step_map%irr_nm_col%mat,errorflag)!
                print*,"print file: ", trim(trim(a_step_map%irr_nm_col%fn))!
                call print_mat_as_grid(trim(a_step_map%deep_perc%fn),domain%header,a_step_map%deep_perc%mat,errorflag)!
                print*,"print file: ", trim(trim(a_step_map%deep_perc%fn))!
                call print_mat_as_grid(trim(a_step_map%runoff%fn),domain%header,a_step_map%runoff%mat,errorflag)!
                print*,"print file: ", trim(trim(a_step_map%runoff%fn))!
                call print_mat_as_grid(trim(a_step_map%et_pot%fn),domain%header,a_step_map%et_pot%mat,errorflag)!
                print*,"print file: ", trim(trim(a_step_map%et_pot%fn))!
                call print_mat_as_grid(trim(a_step_map%et_act%fn),domain%header,a_step_map%et_act%mat,errorflag)!
                print*,"print file: ", trim(trim(a_step_map%et_act%fn))!
                exit!
            else!
                cycle!
            end if!
        end do!
    end subroutine save_step_data!
    !
    subroutine save_step_irrigation(a_step_map,doy,domain,calendar, init_total)
        ! save only irrigation map
        implicit none!
        type(step_map),intent(in)::a_step_map!
        integer,intent(in)::doy!
        type(grid_i),intent(in)::domain!
        integer,dimension(:),intent(in)::calendar!
        integer,intent(in)::init_total
        integer::i,total!
        integer::errorflag!
        !!
        errorflag=0!
        total=init_total!
        !!
        do i=1,size(calendar)!
            total=total+calendar(i)!
            if(doy==total)then    ! last day of the step
                where(domain%mat==domain%header%nan)!
                    a_step_map%irr%mat=real(domain%header%nan)!
                end where!
                call print_mat_as_grid(trim(a_step_map%irr%fn),domain%header,a_step_map%irr%mat,errorflag)!
                print*,"print file: ", trim(trim(a_step_map%irr%fn))!
            else!
                cycle!
            end if!
        end do!
    end subroutine save_step_irrigation!    
    !
    subroutine save_debug_step_data(a_debug_asc,doy,domain,calendar, total_init)
        ! save the results for debug aggregated by step if the current day (doy) is the last day of the step
        type(step_debug_map),intent(in)::a_debug_asc!
        integer,intent(in)::doy!
        type(grid_i),intent(in)::domain!
        integer,dimension(:),intent(in)::calendar!
        integer,intent(in)::total_init
        integer::i,total!
        integer::errorflag!
        !!
        errorflag=0!
        total=total_init!
        !!
        do i=1,size(calendar)!
            total=total+calendar(i)!
            if(doy==total)then    ! last day of the step
                where(domain%mat==domain%header%nan)!
                    a_debug_asc%eva_act%mat=real(domain%header%nan)!
                    a_debug_asc%eff_rain%mat=real(domain%header%nan)!
                    a_debug_asc%perc1%mat=real(domain%header%nan)!
                    a_debug_asc%perc2%mat=real(domain%header%nan)!
                    a_debug_asc%h_soil1%mat=REAL(domain%header%nan)
                    a_debug_asc%h_soil2%mat=REAL(domain%header%nan)
                end where
                call print_mat_as_grid(trim(a_debug_asc%eva_act%fn),domain%header,a_debug_asc%eva_act%mat,errorflag)!
                print*,"print file: ", trim(trim(a_debug_asc%eva_act%fn))!
                call print_mat_as_grid(trim(a_debug_asc%eff_rain%fn),domain%header,a_debug_asc%eff_rain%mat,errorflag)!
                print*,"print file: ", trim(trim(a_debug_asc%eff_rain%fn))!
                call print_mat_as_grid(trim(a_debug_asc%perc1%fn),domain%header,a_debug_asc%perc1%mat,errorflag)!
                print*,"print file: ", trim(trim(a_debug_asc%perc1%fn))!
                call print_mat_as_grid(trim(a_debug_asc%perc2%fn),domain%header,a_debug_asc%perc2%mat,errorflag)!
                print*,"print file: ", trim(trim(a_debug_asc%perc2%fn))!
                call print_mat_as_grid(trim(a_debug_asc%h_soil1%fn),domain%header,a_debug_asc%h_soil1%mat,errorflag)!
                print*,"print file: ", trim(trim(a_debug_asc%h_soil1%fn))!
                call print_mat_as_grid(trim(a_debug_asc%h_soil2%fn),domain%header,a_debug_asc%h_soil2%mat,errorflag)!
                print*,"print file: ", trim(trim(a_debug_asc%h_soil2%fn))!
                exit!
            else!
                cycle!
            end if!
        end do!
    end subroutine save_debug_step_data!
    !
    subroutine save_yearly_data(yr_map,domain,debug)!
        ! save annual outputs (water balance variable and efficiency)
        implicit none!
        type(annual_map),intent(in)::yr_map!
        type(grid_i),intent(in)::domain!
        logical,intent(in)::debug
        !!
        integer::errorflag!
        !!
        errorflag=0!
        where(domain%mat==domain%header%nan)!
            yr_map%rain%mat=real(domain%header%nan)!
            yr_map%rain_crop_season%mat=real(domain%header%nan)
            yr_map%irr%mat=real(domain%header%nan)!
            yr_map%irr_loss%mat=real(domain%header%nan)!
            yr_map%eva_act_crop_season%mat=real(domain%header%nan)!
            yr_map%eva_pot_crop_season%mat=real(domain%header%nan)!
            yr_map%transp_act%mat=real(domain%header%nan)!
            yr_map%transp_pot%mat=real(domain%header%nan)!
            yr_map%runoff%mat=real(domain%header%nan)!
            yr_map%net_flux_gw%mat=real(domain%header%nan)!
            yr_map%total_eff%mat=real(domain%header%nan)!
            yr_map%n_irr_events%mat=real(domain%header%nan)
            yr_map%h_irr_mean%mat=real(domain%header%nan)
        end where!
        !
        call print_mat_as_grid(trim(yr_map%rain%fn),domain%header,yr_map%rain%mat,errorflag)!
        print*,"print file: ", trim(trim(yr_map%rain%fn))!
        call print_mat_as_grid(trim(yr_map%irr%fn),domain%header,yr_map%irr%mat,errorflag)!
        print*,"print file: ", trim(trim(yr_map%irr%fn))!
        call print_mat_as_grid(trim(yr_map%irr_loss%fn),domain%header,yr_map%irr_loss%mat,errorflag)!
        print*,"print file: ", trim(trim(yr_map%irr_loss%fn))!
        call print_mat_as_grid(trim(yr_map%eva_act_crop_season%fn),domain%header,yr_map%eva_act_crop_season%mat,errorflag)!
        print*,"print file: ", trim(trim(yr_map%eva_act_crop_season%fn))!
        call print_mat_as_grid(trim(yr_map%eva_pot_crop_season%fn),domain%header,yr_map%eva_pot_crop_season%mat,errorflag)!
        print*,"print file: ", trim(trim(yr_map%eva_pot_crop_season%fn))!
        call print_mat_as_grid(trim(yr_map%transp_act%fn),domain%header,yr_map%transp_act%mat,errorflag)!
        print*,"print file: ", trim(trim(yr_map%transp_act%fn))!
        call print_mat_as_grid(trim(yr_map%transp_pot%fn),domain%header,yr_map%transp_pot%mat,errorflag)!
        print*,"print file: ", trim(trim(yr_map%transp_pot%fn))!
        call print_mat_as_grid(trim(yr_map%runoff%fn),domain%header,yr_map%runoff%mat,errorflag)!
        print*,"print file: ", trim(trim(yr_map%runoff%fn))!
        call print_mat_as_grid(trim(yr_map%net_flux_gw%fn),domain%header,yr_map%net_flux_gw%mat,errorflag)!
        print*,"print file: ", trim(trim(yr_map%net_flux_gw%fn))!
        call print_mat_as_grid(trim(yr_map%total_eff%fn),domain%header,yr_map%total_eff%mat,errorflag)!
        print*,"print file: ", trim(trim(yr_map%total_eff%fn))!
        call print_mat_as_grid(trim(yr_map%n_irr_events%fn),domain%header,yr_map%n_irr_events%mat,errorflag)!
        print*,"print file: ", trim(trim(yr_map%n_irr_events%fn))!
        call print_mat_as_grid(trim(yr_map%h_irr_mean%fn),domain%header,yr_map%h_irr_mean%mat,errorflag)!
        print*,"print file: ", trim(trim(yr_map%h_irr_mean%fn))!
        !
        if (debug .eqv. .true.) then
            call print_mat_as_grid(trim(yr_map%rain_crop_season%fn),domain%header,yr_map%rain_crop_season%mat,errorflag)!
            print*,"print file: ", trim(trim(yr_map%rain_crop_season%fn))!
        end if
    end subroutine save_yearly_data!
    !
    subroutine save_yield_data(yield,domain)!
        type(yield_t),intent(in)::yield!
        type(grid_i),intent(in)::domain!
        integer::j!
        integer::errorflag!
        character(len=55)::strj!
        !!
        errorflag=0!
        
        do j=1, size(yield%yield_pot%mat,3)
            yield%yield_pot%mat(:,:,j) = merge(yield%yield_pot%mat(:,:,j),dble(domain%header%nan),domain%mat/=domain%header%nan)
            yield%yield_pot%mat(:,:,j) = merge(yield%yield_pot%mat(:,:,j),dble(domain%header%nan),yield%yield_pot%mat(:,:,j)/=0)
            yield%yield_act%mat(:,:,j) = merge(yield%yield_act%mat(:,:,j),dble(domain%header%nan),domain%mat/=domain%header%nan)
            yield%yield_act%mat(:,:,j) = merge(yield%yield_act%mat(:,:,j),dble(domain%header%nan),yield%yield_act%mat(:,:,j)/=0)
            yield%biomass_pot%mat(:,:,j) = merge(yield%biomass_pot%mat(:,:,j), &
                & dble(domain%header%nan),domain%mat/=domain%header%nan)
            yield%biomass_pot%mat(:,:,j) = merge(yield%biomass_pot%mat(:,:,j), &
                & dble(domain%header%nan),yield%biomass_pot%mat(:,:,j)/=0)
        end do
        
        do j=1, size(yield%yield_pot%mat,3)
            write(strj,*)j
            call print_mat_as_grid(trim(trim(yield%yield_pot%fn)//"_"//trim(adjustl(strj))//".asc"), &
                & domain%header,yield%yield_pot%mat(:,:,j),errorflag)!
            print*,"print file: ", (trim(trim(yield%yield_pot%fn)//"_"//trim(adjustl(strj))//".asc"))!
            call print_mat_as_grid(trim(trim(yield%yield_act%fn)//"_"//trim(adjustl(strj))//".asc"), &
                & domain%header,yield%yield_act%mat(:,:,j),errorflag)!
            print*,"print file: ", (trim(trim(yield%yield_act%fn)//"_"//trim(adjustl(strj))//".asc"))!
            call print_mat_as_grid(trim(trim(yield%biomass_pot%fn)//"_"//trim(adjustl(strj))//".asc"), &
                & domain%header,yield%biomass_pot%mat(:,:,j),errorflag)!
            print*,"print file: ", (trim(trim(yield%biomass_pot%fn))//"_"//trim(adjustl(strj))//".asc")!
        end do
    end subroutine save_yield_data!
    !
    subroutine save_annual_irrigation_data(yr_map,domain)!
        type(annual_map),intent(in)::yr_map!
        type(grid_i),intent(in)::domain!
        !!
        integer::errorflag!
        !!
        errorflag=0!
        where(domain%mat==domain%header%nan)!
            yr_map%irr%mat=real(domain%header%nan)!
        end where!
        call print_mat_as_grid(trim(yr_map%irr%fn),domain%header,yr_map%irr%mat,errorflag)!
        print*,"print file: ", trim(trim(yr_map%irr%fn))!
    end subroutine save_annual_irrigation_data!    
    !
    subroutine save_annual_debug_data(yr_dbg_map, domain)!
        type(annual_debug_map),intent(in)::yr_dbg_map!
        type(grid_i),intent(in)::domain!
        integer::errorflag!
        !!
        errorflag=0!
        where(domain%mat==domain%header%nan)!
            yr_dbg_map%eva_act_tot%mat=real(domain%header%nan)!
            yr_dbg_map%rain_eff%mat=real(domain%header%nan)!
            yr_dbg_map%iter1%mat=domain%header%nan
            yr_dbg_map%iter2%mat=domain%header%nan
        end where!
        call print_mat_as_grid(trim(yr_dbg_map%eva_act_tot%fn),domain%header,yr_dbg_map%eva_act_tot%mat,errorflag)!
        print*,"print file: ", trim(trim(yr_dbg_map%eva_act_tot%fn))!
        call print_mat_as_grid(trim(yr_dbg_map%rain_eff%fn),domain%header,yr_dbg_map%rain_eff%mat,errorflag)!
        print*,"print file: ", trim(trim(yr_dbg_map%rain_eff%fn))!
        call print_mat_as_grid(trim(yr_dbg_map%iter1%fn),domain%header,yr_dbg_map%iter1%mat,errorflag)!
        print*,"print file: ", trim(trim(yr_dbg_map%iter1%fn))!
        call print_mat_as_grid(trim(yr_dbg_map%iter2%fn),domain%header,yr_dbg_map%iter2%mat,errorflag)!
        print*,"print file: ", trim(trim(yr_dbg_map%iter2%fn))!
    end subroutine save_annual_debug_data!
    !
    subroutine save_yield_debug_data(yield,domain)!
        implicit none!
        type(yield_t),intent(in)::yield!
        type(grid_i),intent(in)::domain!
        integer::errorflag!
        integer::i,j
        character(len=55)::stri,strj
        !!
        errorflag=0!
        
        do i =1, size(yield%T_act_sum%mat,3)
            write(stri,*)i
            do j=1, size(yield%T_act_sum%mat,4)
                write(strj,*)j
                call print_mat_as_grid &
                    & (trim(trim(yield%T_act_sum%fn)//"_"//trim(adjustl(stri))//"_"//trim(adjustl(strj))//".asc"), &
                    & domain%header,yield%T_act_sum%mat(:,:,i,j),errorflag)!
                print*,"print file: ", trim(trim(yield%T_act_sum%fn)//"_"//trim(adjustl(stri))//"_"//trim(adjustl(strj))//".asc")!
                call print_mat_as_grid &
                    & (trim(trim(yield%T_pot_sum%fn)//"_"//trim(adjustl(stri))//"_"//trim(adjustl(strj))//".asc"), &
                    & domain%header,yield%T_pot_sum%mat(:,:,i,j),errorflag)!
                print*,"print file: ", trim(trim(yield%T_pot_sum%fn)//"_"//trim(adjustl(stri))//"_"//trim(adjustl(strj))//".asc")!
            end do
        end do
        do i=1, size(yield%f_WS%mat,3)
            write(stri,*)i
            call print_mat_as_grid(trim(trim(yield%f_WS%fn)//"_"//trim(adjustl(stri))//".asc"), &
                & domain%header,yield%f_WS%mat(:,:,i),errorflag)!
            print*,"print file: ", (trim(trim(yield%f_WS%fn)//"_"//trim(adjustl(stri))//".asc"))!
            call print_mat_as_grid(trim(trim(yield%f_WS_stage%fn)//"_"//trim(adjustl(stri))//".asc"), &
                & domain%header,yield%f_WS_stage%mat(:,:,i),errorflag)!
            print*,"print file: ", (trim(trim(yield%f_WS_stage%fn)//"_"//trim(adjustl(stri))//".asc"))!
            call print_mat_as_grid(trim(trim(yield%f_HS%fn)//"_"//trim(adjustl(stri))//".asc"), &
                & domain%header,yield%f_HS%mat(:,:,i),errorflag)!
            print*,"print file: ", (trim(trim(yield%f_HS%fn))//"_"//trim(adjustl(stri))//".asc")!
        end do
    end subroutine save_yield_debug_data!
    !
    subroutine assign_step_map(a_step_map,a)!
        type(step_map),intent(out)::a_step_map!
        real(dp),intent(in)::a!
        
        a_step_map%runoff%mat = a!
        a_step_map%rain%mat = a!
        a_step_map%transp_act%mat = a!
        a_step_map%transp_pot%mat = a!
        a_step_map%irr%mat = a!
        a_step_map%irr_loss%mat = a!
        a_step_map%cap_rise%mat = a!
        a_step_map%irr_nm_priv%mat = a!
        a_step_map%irr_nm_col%mat = a!
        a_step_map%deep_perc%mat = a!
        a_step_map%runoff%mat = a!
        a_step_map%et_pot%mat = a!
        a_step_map%et_act%mat = a!
    end subroutine assign_step_map!
    !
    subroutine assign_step_debug_map(a_debug_map,a)!
        type(step_debug_map),intent(out)::a_debug_map!
        real(dp),intent(in)::a!
        !!
        a_debug_map%eva_act%mat = a!
        a_debug_map%perc1%mat = a!
        a_debug_map%eff_rain%mat = a!
        a_debug_map%perc2%mat = a!
        a_debug_map%h_soil1%mat = a!
        a_debug_map%h_soil2%mat = a!
    end subroutine assign_step_debug_map!
    !
    subroutine assign_annual_map(yr_map,a)!
        type(annual_map),intent(out)::yr_map!
        real(dp),intent(in)::a!
        !!
        yr_map%rain%mat = a!
        yr_map%rain_crop_season%mat = a!
        yr_map%irr%mat = a!
        yr_map%irr_loss%mat = a!
        yr_map%eva_pot_crop_season%mat = a!
        yr_map%eva_act_crop_season%mat = a!
        yr_map%transp_act%mat = a!
        yr_map%transp_pot%mat = a!
        yr_map%runoff%mat = a!
        yr_map%net_flux_gw%mat = a!
        yr_map%total_eff%mat = a!
        yr_map%n_irr_events%mat = a!
        yr_map%h_irr_mean%mat = a!
    end subroutine assign_annual_map!
    !
    subroutine assign_yield_map(yield_map,a)!
        type(yield_t),intent(out)::yield_map!
        real(dp),intent(in)::a!
        !!
        yield_map%biomass_pot%mat = a!
        yield_map%yield_pot%mat = a!
        yield_map%yield_act%mat = a!
        yield_map%f_WS%mat = a!
        yield_map%f_WS_stage%mat = a!
        yield_map%f_HS%mat = a!
        yield_map%f_HS_sum%mat = a!
        yield_map%transp_ratio_sum%mat = a!
        yield_map%T_act_sum%mat = a!
        yield_map%T_pot_sum%mat = a!
        yield_map%dev_stage%mat = a!
    end subroutine assign_yield_map!
    !
    subroutine assign_annual_debug_map(yr_debug_map,a)!
        type(annual_debug_map),intent(out)::yr_debug_map!
        real(dp),intent(in)::a!
        !!
        yr_debug_map%eva_act_tot%mat = a!
        yr_debug_map%rain_eff%mat = a!
        yr_debug_map%iter1%mat = a
        yr_debug_map%iter2%mat = a
    end subroutine assign_annual_debug_map!

    subroutine save_irr_unit_data(doy,out_tables,irr_units)
        integer, intent(in):: doy
        type(output_table_list), intent(in)::out_tables
        type(irr_units_table),dimension(:),allocatable, intent(in)::irr_units
        integer:: i
        ! %AB% mettere un check per la stampa
        write(out_tables%q_surplus%unit,*)doy,'; ',(irr_units(i)%q_surplus,'; ',i=1,size(irr_units))
        write(out_tables%q_irr%unit,*)doy,'; ',(irr_units(i)%q_day,'; ',i=1,size(irr_units))
        write(out_tables%q_rem%unit,*)doy,'; ',(irr_units(i)%q_rem,'; ',i=1,size(irr_units))
        write(out_tables%q_un_priv%unit,*)doy,'; ',(irr_units(i)%q_un_priv,'; ',i=1,size(irr_units))
        write(out_tables%n_inv_cells%unit,*)doy,'; ',(irr_units(i)%n_irrigated_cells,'; ',i=1,size(irr_units))
        write(out_tables%q_un_col%unit,*)doy,'; ',(irr_units(i)%q_un_coll,'; ',i=1,size(irr_units))
        ! q_act_fld(4)
    end subroutine save_irr_unit_data

    subroutine save_irr_unit_debug_data(doy,out_tables,irr_units)
        integer, intent(in):: doy
        type(output_table_list), intent(in)::out_tables
        type(irr_units_table),dimension(:),allocatable, intent(in)::irr_units
        integer:: i

        ! print discharges fro each irrigation units
        write(out_tables%q_irr_units%unit,*)doy,'; ','Monitored sources (i)','; ',(irr_units(i)%q_act_fld(1),'; ',i=1,size(irr_units))!
        write(out_tables%q_irr_units%unit,*)doy,'; ','Monitored sources (ii)','; ',(irr_units(i)%q_act_fld(2),'; ',i=1,size(irr_units))!
        write(out_tables%q_irr_units%unit,*)doy,'; ','Internal reuse','; ',(irr_units(i)%q_act_fld(3),'; ',i=1,size(irr_units))!
        write(out_tables%q_irr_units%unit,*)doy,'; ','Collective runtime sources','; ',(irr_units(i)%q_act_fld(4),'; ',i=1,size(irr_units))!
        
    end subroutine save_irr_unit_debug_data

end module cli_save_outputs!
