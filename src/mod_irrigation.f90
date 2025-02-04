module mod_irrigation
    use mod_constants, only: dp
    use mod_utility, only: get_value_index
    use mod_grid
    use mod_common, only: balance1_matrices,balance2_matrices, spatial_info
    use mod_crop_phenology, only: crop_pars_matrices
    use mod_parameters, only: parameters, scheduled_irrigation, irr_units_table, source_info, water_sources_table
    
    implicit none

    contains

    subroutine calc_irrigation_losses(a_loss, b_loss, c_loss, wind_vel, temp_ave,losses)!
        ! calculate the water losses during the irrigation event
        ! according to the formula
        ! losses = a + (b * wind_vel) + (c * temp_ave)
        ! where a, b and c are coefficient
        real(dp),dimension(:,:),intent(in)::a_loss,b_loss,c_loss
        real(dp),dimension(:,:),intent(in)::wind_vel
        real(dp),dimension(:,:),intent(in)::temp_ave
        real(dp),dimension(:,:),intent(inout)::losses
        
        losses = a_loss + b_loss * wind_vel + c_loss * temp_ave
        losses = merge(losses,0.0D0,losses>0.0D0)
    end subroutine calc_irrigation_losses
    
    subroutine calc_perc_booster_pars(info_spat, met, quantiles)!
        ! calculate the parameter for the percolation booster
        type(spatial_info),intent(inout)::info_spat!
        type(par_method),dimension(:),intent(in)::met!
        real(dp),dimension(info_spat%domain%header%imax,info_spat%domain%header%jmax)::a_min,a_max
        real(dp),dimension(info_spat%domain%header%imax,info_spat%domain%header%jmax)::b_min,b_max
        real(dp),dimension(:,:),intent(in)::quantiles
        integer::i!
        !!
        a_min=0.d0;a_max=0.d0!
        b_min=0.d0;b_max=0.d0!
        a_min=id_to_par(info_spat%irr_meth_id,met(:)%a_min)!
        a_max=id_to_par(info_spat%irr_meth_id,met(:)%a_max)!
        b_min=id_to_par(info_spat%irr_meth_id,met(:)%b_min)!
        b_max=id_to_par(info_spat%irr_meth_id,met(:)%b_max)!
        !!
        do i=1,2!
            info_spat%a_perc(i) = info_spat%domain!
            info_spat%a_perc(i)%mat = 0.d0!
            info_spat%b_perc(i) = info_spat%domain!
            info_spat%b_perc(i)%mat = 0.d0!

            where(info_spat%k_sat(i)%mat<=quantiles(i,1))
                info_spat%a_perc(i)%mat=a_max
                info_spat%b_perc(i)%mat=b_max
            else where(info_spat%k_sat(i)%mat>=quantiles(i,2))
                info_spat%a_perc(i)%mat=a_min
                info_spat%b_perc(i)%mat=b_min
            else where
                info_spat%a_perc(i)%mat = &
                    & info_spat%k_sat(i)%mat*((a_max-a_min)/(quantiles(i,1)-quantiles(i,2))) + &
                    & a_min - quantiles(i,2) * (a_max-a_min)/(quantiles(i,1)-quantiles(i,2))!
                info_spat%b_perc(i)%mat = &
                    & info_spat%k_sat(i)%mat*((b_max-b_min)/(quantiles(i,1)-quantiles(i,2))) + &
                    & b_min - quantiles(i,2) * (b_max-b_min)/(quantiles(i,1)-quantiles(i,2))!
            end where!
            
            call set_default_par(info_spat%a_perc(i),info_spat%domain,1.0D0)
            call set_default_par(info_spat%b_perc(i),info_spat%domain,1.0D0)
        end do!
        
        ! a_perc and b_perc are NaN for not irrigated cells
        
    end subroutine calc_perc_booster_pars!


    subroutine update_adj_perco_parameters(info_spat, matrice_irr, day_from_irr, adj_perc_par)
        ! calculate adjustment parameter of the percolation model (eq 1.44)!
        ! from the last irrigation event
        type(spatial_info),intent(in)::info_spat
        real(dp),dimension(:,:),intent(in)::matrice_irr
        integer,dimension(:,:),intent(inout)::day_from_irr
        real(dp),dimension(:,:,:),intent(inout)::adj_perc_par
        
        where(matrice_irr==0)!
            where(day_from_irr==-9999)!
                day_from_irr=-9999!
            else where!
                day_from_irr = day_from_irr +1!
            end where!
        else where!
            day_from_irr=0!
        end where!

        ! %EAC%: limit day_from_irr whene day_from_irr= 20 to prevent calculation error, TODO: to be fixed with a better solution
        ! 20 should be safe for all OS and CPU
        where(day_from_irr>20)
            day_from_irr=20
        end where

        ! espeperc --> H_(x,irr)=(1+a_x e^(-t⋅b_x ))H_x
        adj_perc_par(:,:,1)=merge(1.0D0,1+info_spat%a_perc(1)%mat*exp(-1.0*day_from_irr*info_spat%b_perc(1)%mat),day_from_irr==-9999)!
        adj_perc_par(:,:,2)=merge(1.0D0,1+info_spat%a_perc(2)%mat*exp(-1.0*day_from_irr*info_spat%b_perc(2)%mat),day_from_irr==-9999)!
    
    end subroutine


    subroutine irrigation_need_fix(info_spat, h_irr, bil2, bil2_old, bil1_old, pheno, &!
        & eff_rain, xrice_ksat, day_from_irr, adj_perc_par)!
        ! calculate irrigation needs at field capacity and fixed volume (defined by irrigation methods)
        implicit none!
        real(dp),dimension(:,:,:),intent(out)::h_irr
        type(spatial_info),intent(in)::info_spat
        type(balance1_matrices),intent(in)::bil1_old
        type(balance2_matrices),intent(in)::bil2,bil2_old
        type(crop_pars_matrices),intent(in)::pheno
        integer,dimension(:,:),intent(inout)::day_from_irr
        real(dp),dimension(:,:,:),intent(inout)::adj_perc_par ! parameter that adjust percolation
        real(dp),dimension(:,:),intent(in)::eff_rain
        real(dp),intent(in)::xrice_ksat     ! ksat of the transpirative layer for rice
        real(dp),dimension(size(info_spat%domain%mat,1),size(info_spat%domain%mat,2))::h_irr_temp!
        integer::i,j!
        ! init
        h_irr = 0.!
        ! the irrigation event is estimated only during the irrigation season
        h_irr_temp = 0.!

        where(pheno%irrigation_class==1 .and. pheno%k_cb>0)
            !!! %RR% %CG% %EAC% %AB% (feb-23) network efficiency no more considered in need mode [info_spat%eff_rete%mat]
            where(pheno%cn_class==7)                   ! irrigation matrix for rice
                !h_irr_temp = (bil1_old%h_eva + bil2_old%h_transp_pot)!
                ! %CG% add water to fill soil, ponding and ET
                h_irr_temp = ((info_spat%h_meth%mat-  bil1_old%h_pond)+ & ! fill the ponding layer
                                (info_spat%theta(1)%sat%mat*bil1_old%d_e*1000-bil1_old%h_soil)+ &!
                                (info_spat%theta(2)%sat%mat*bil2%d_t*1000-bil2_old%h_soil) + & ! fill the soil layer to saturation
                                (bil1_old%h_eva + bil2_old%h_transp_pot))/ & ! fill ET
                                1.0 ! don't consider efficiency
            ! TODO: %EAC%: small edit to manage submerged condition
            else where ((info_spat%h_maxpond%mat>10.0D0) .and. (bil1_old%h_pond<(0.9*info_spat%h_maxpond%mat))) ! crop submerged, only if h_maxpond> 10 mm
                h_irr_temp = info_spat%h_meth%mat
            else where ((bil1_old%h_soil + bil2_old%h_soil) < bil2%h_raw_sup) ! all other crops
                ! TODO: %AB% why only 2nd layer RAWbig?
                ! else where (bil2_old%tmm<bil2%RAWbig)
                h_irr_temp = info_spat%h_meth%mat
            end where!
        end where

        call irrigate_rice(h_irr_temp,pheno,eff_rain,xrice_ksat)!

        ! split the irrigation matrix for each irrigation method 
        ! n is the number of irrigation methods
        forall(i=1:size(info_spat%domain%mat,1),j=1:size(info_spat%domain%mat,2),&!
            & info_spat%irr_meth_id%mat(i,j)/=info_spat%irr_meth_id%header%nan) &!
            & h_irr(i,j,info_spat%irr_meth_id%mat(i,j)) = h_irr_temp(i,j)!
        
        call update_adj_perco_parameters(info_spat, h_irr_temp, day_from_irr, adj_perc_par)
        
    end subroutine irrigation_need_fix!

    subroutine irrigation_need_fc(info_spat,h_irr,bil2,bil2_old,bil1_old,pheno,&!
        & eff_rain,xrice_ksat,day_from_irr,adj_perc_par)!
        ! calculate irrigation needs at field capacity
        implicit none!
        real(dp),dimension(:,:,:),intent(out)::h_irr!
        type(spatial_info),intent(in)::info_spat!
        type(balance1_matrices),intent(in)::bil1_old!
        type(balance2_matrices),intent(in)::bil2,bil2_old!
        type(crop_pars_matrices),intent(in)::pheno!
        integer,dimension(:,:),intent(inout)::day_from_irr!
        real(dp),dimension(:,:,:),intent(inout)::adj_perc_par!
        real(dp),dimension(:,:),intent(in)::eff_rain
        integer::i,j
        real(dp),intent(in)::xrice_ksat     ! ksat of the transpirative layer for rice
        real(dp),dimension(size(info_spat%domain%mat,1),size(info_spat%domain%mat,2))::h_irr_temp!
        
        ! init
        h_irr = 0.
        h_irr_temp = 0.
        !!! %RR% %CG% %EAC% %AB% (feb-23) network efficiency no more considered in need mode [info_spat%eff_rete%mat]
        where (pheno%irrigation_class==1 .and. pheno%k_cb>0)
            where (pheno%cn_class==7)                      ! %AB% irrigation matrix for rice
                h_irr_temp = (bil1_old%h_eva + bil2_old%h_transp_pot)/(info_spat%eff_met%mat)

                ! TODO: info_spat%h_meth%mat is not initialized in FC mode
                !h_irr_temp = ((info_spat%h_meth%mat-bil1_old%h_pond)+ & ! fill the ponding layer
                !             (info_spat%theta(1)%sat%mat*bil1_old%d_e*1000-bil1_old%h_soil)+ &!
                !             (info_spat%theta(2)%sat%mat*bil2%d_t*1000-bil2_old%h_soil) + & ! fill the soil layer to saturation
                !             (bil1_old%h_eva + bil2_old%h_transp_pot))/ & ! fill ET
                !             (info_spat%eff_met%mat) ! consider also efficiency

!~             else where (bil2_old%tmm<bil2%RAWbig)       !all other crops
            else where ((bil1_old%h_soil + bil2_old%h_soil) < bil2%h_raw_sup)       !all other crops
                ! TODO: %AB% why only 2nd layer RAWbig?
                where(bil1_old%h_soil <= info_spat%theta(1)%FC%mat*bil1_old%d_e*1000)!
                    h_irr_temp = ((info_spat%theta(1)%FC%mat*bil1_old%d_e*1000-bil1_old%h_soil)+ &!
                        & (info_spat%theta(2)%FC%mat*bil2%d_t*1000-bil2_old%h_soil))/(info_spat%eff_met%mat)
                else where!
                    h_irr_temp = (info_spat%theta(2)%FC%mat*bil2%d_t*1000-bil2_old%h_soil)/&
                        & (info_spat%eff_met%mat)
                end where!
            end where!
        end where

        call irrigate_rice(h_irr_temp,pheno,eff_rain,xrice_ksat)!

        ! split the irrigation matrix for each irrigation method 
        ! n is the number of irrigation methods
        forall(i=1:size(info_spat%domain%mat,1),j=1:size(info_spat%domain%mat,2),&!
            & info_spat%irr_meth_id%mat(i,j)/=info_spat%irr_meth_id%header%nan) &!
            & h_irr(i,j,info_spat%irr_meth_id%mat(i,j)) = h_irr_temp(i,j)!
        
        call update_adj_perco_parameters(info_spat, h_irr_temp, day_from_irr, adj_perc_par)
        
    end subroutine irrigation_need_fc!

    subroutine irrigation_scheduled(info_spat, doy_cur, year_cur, sch_irr, pheno, h_irr, &
        & day_from_irr, adj_perc_par, debug, bil1_old, bil2, bil2_old, &
        & a_loss, b_loss, c_loss, wind_vel, temp_ave,losses,eff_rain,xrice_ksat)
        ! spread irrigation height when scheduled
        ! TODO: need for testing
        ! TODO: include losses calculation inside
        
        type(spatial_info),intent(in)::info_spat!
        integer,intent(in)::doy_cur ! check the day as day of the year
        integer,intent(in)::year_cur ! check the year
        type(scheduled_irrigation),dimension(:),intent(in)::sch_irr!
        type(crop_pars_matrices),intent(in)::pheno!
        real(dp),dimension(:,:,:),intent(out)::h_irr!
        integer,dimension(:,:),intent(inout)::day_from_irr!
        real(dp),dimension(:,:,:),intent(inout)::adj_perc_par!
        logical:: debug
        type(balance1_matrices),intent(inout):: bil1_old!
        type(balance2_matrices),intent(inout)::bil2, bil2_old!
        
        real(dp),dimension(:,:),intent(in)::a_loss,b_loss,c_loss
        real(dp),dimension(:,:),intent(in)::wind_vel
        real(dp),dimension(:,:),intent(in)::temp_ave
        real(dp),dimension(:,:),intent(inout)::losses
        real(dp),dimension(:,:),intent(in)::eff_rain
        real(dp),intent(in)::xrice_ksat     ! ksat of the transpirative layer for rice
        
        real(dp),dimension(size(info_spat%domain%mat,1),size(info_spat%domain%mat,2))::h_irr_temp!
        
        integer::i,j

        !init
        h_irr = 0.!
        h_irr_temp = 0.!
        losses = 0.

        ! loop throw scheduled irrigation list and update 'irrigation'
        do i=1,size(sch_irr)!
            if (sch_irr(i)%year == year_cur .and. sch_irr(i)%doy == doy_cur) then
                if (debug .eqv. .true.) then
                    print *,'current year = ',  year_cur, '; current day = ', doy_cur, &
                        & '; irr.unit.id= ', sch_irr(i)%irr_unit_id, '; water depth = ', sch_irr(i)%h_irr
                end if
                if (sch_irr(i)%h_irr > 0) then
                    ! add water quantity  as scheduled 
                    where (info_spat%irr_unit_id%mat == sch_irr(i)%irr_unit_id .and. pheno%irrigation_class==1)
                        h_irr_temp = sch_irr(i)%h_irr ! %EAC%: like the old version
                    end where

                    call calc_irrigation_losses(a_loss, b_loss, c_loss, wind_vel, temp_ave,losses)
                    
                else if (sch_irr(i)%h_irr == -1) then
                    ! add water with fixed volume
                    where (info_spat%irr_unit_id%mat == sch_irr(i)%irr_unit_id .and. pheno%irrigation_class==1)
                        where(pheno%cn_class==7)                   ! irrigation matrix for rice
                            !h_irr_temp = (bil1_old%h_eva+bil2_old%h_transp_pot)!
                            h_irr_temp = ((info_spat%h_meth%mat-bil1_old%h_pond)+ & ! fill the ponding layer
                                        (info_spat%theta(1)%sat%mat*bil1_old%d_e*1000-bil1_old%h_soil)+ &!
                                        (info_spat%theta(2)%sat%mat*bil2%d_t*1000-bil2_old%h_soil) + & ! fill the soil layer to saturation
                                        (bil1_old%h_eva + bil2_old%h_transp_pot))/ & ! fill ET
                                        1.0 ! don't consider efficiency
                        else where ((bil1_old%h_soil + bil2_old%h_soil) < bil2%h_raw_sup) ! all other crops
                            ! TODO: %AB% why only 2nd layer RAWbig?
                            ! else where (bil2_old%tmm<bil2%RAWbig)
                            h_irr_temp = info_spat%h_meth%mat
                        end where
                    end where
                    
                    call calc_irrigation_losses(a_loss, b_loss, c_loss, wind_vel, temp_ave,losses)

                else if (sch_irr(i)%h_irr == -10) then
                    ! add water to field capacity
                    where (info_spat%irr_unit_id%mat == sch_irr(i)%irr_unit_id .and. pheno%irrigation_class==1)
                        where (pheno%cn_class==7)                      ! %AB% irrigation matrix for rice
                            !h_irr_temp = (bil1_old%h_eva+bil2_old%h_transp_pot)/(info_spat%eff_met%mat)
                            h_irr_temp = ((info_spat%h_meth%mat-bil1_old%h_pond)+ & ! fill the ponding layer
                                        (info_spat%theta(1)%sat%mat*bil1_old%d_e*1000-bil1_old%h_soil)+ &!
                                        (info_spat%theta(2)%sat%mat*bil2%d_t*1000-bil2_old%h_soil) + & ! fill the soil layer to saturation
                                        (bil1_old%h_eva + bil2_old%h_transp_pot))/ & ! fill ET
                                        (info_spat%eff_met%mat) ! consider also efficiency
                            ! else where (bil2_old%tmm<bil2%RAWbig)       !all other crops
                        else where ((bil1_old%h_soil + bil2_old%h_soil) < bil2%h_raw_sup)       !all other crops
                            ! TODO: %AB% why only 2nd layer RAWbig?
                            where(bil1_old%h_soil <= info_spat%theta(1)%FC%mat*bil1_old%d_e*1000)!
                                h_irr_temp = ((info_spat%theta(1)%FC%mat*bil1_old%d_e*1000-bil1_old%h_soil)+ &!
                                    & (info_spat%theta(2)%FC%mat*bil2%d_t*1000-bil2_old%h_soil))/(info_spat%eff_met%mat)
                            else where!
                                h_irr_temp = (info_spat%theta(2)%FC%mat*bil2%d_t*1000-bil2_old%h_soil)/&
                                    & (info_spat%eff_met%mat)
                            end where!
                        end where!
                    end where

                    losses = 0.

                    ! OLD: force soil water content to field capacity
                    ! set the soil water content [mm]!
                    !bil1%tmm = info_spat%theta(1)%fc%mat*1000*ze_fix!
                    !bil2%tmm = info_spat%theta(2)%fc%mat*1000*bil2%zr        !
                    ! set the soil moisture [m3/m3]!
                    !bil1%tvol = info_spat%theta(1)%fc%mat!
                    !bil2%tvol = info_spat%theta(2)%fc%mat!
                ! else do nothing
                end if
            end if
        end do!

        call irrigate_rice(h_irr_temp,pheno,eff_rain,xrice_ksat)!

        ! split the irrigation matrix for each irrigation method 
        ! n is the number of irrigation methods
        forall(i=1:size(info_spat%domain%mat,1),j=1:size(info_spat%domain%mat,2),&
            & info_spat%irr_meth_id%mat(i,j)/=info_spat%irr_meth_id%header%nan) &
            & h_irr(i,j,info_spat%irr_meth_id%mat(i,j)) = h_irr_temp(i,j)!

        !print*,'max h_irr',maxval(h_irr)
        
        call update_adj_perco_parameters(info_spat, h_irr_temp, day_from_irr, adj_perc_par)
        
    end subroutine irrigation_scheduled
    
    subroutine irrigate_rice(h_irr,pheno,eff_rain,k_sat)!
        ! calculate irrigation event for paddy fields
        real(dp),dimension(:,:),intent(inout)::h_irr    ! TODO: externally set to equal the potential crop ET,
                                                        !       the soil water content to saturation and the ponding height
        type(crop_pars_matrices),intent(in)::pheno      ! 
        real(dp),dimension(:,:),intent(in)::eff_rain    ! effective rain
        real(dp),intent(in)::k_sat!
        real(dp)::h_irr_min!
        !!
        ! init the minimum irrigation height (equal to infiltration)
        h_irr_min = 10.*24.*k_sat ! k_sat is in mm
        where(pheno%irrigation_class==1 .and. pheno%cn_class==7)  
            h_irr = max(0.0D0,h_irr_min + h_irr - eff_rain)  ! %CG%: irrigation compensate ET and Percolation, minus effective rain 
            !where(h_irr<=h_irr_min) h_irr = h_irr_min       ! irrigation height at least equal to the minimum irrigation height
            !where(eff_rain>=h_irr) h_irr = 0.               ! if effective precipitation > irrigation height -> zero irrigation height
            ! adjust start and end of the season
            where(pheno%k_cb<pheno%k_cb_old) h_irr = 0       ! %CG%: if last crop period don't irrigate
        end where!
    end subroutine irrigate_rice!

    subroutine irrigation_use(domain, irr_units_map, irr_class, method, &!
        & irr_units, transp_pot, h_old, h_raw_coll,h_raw_half,h_raw,h_raw_priv,h_irr,doy,priv_irr,coll_irr, &!
        & day_from_irr,esp_perc,am_perc,bm_perc, f_shape_area, cell_area, h_met,irr_starts,irr_ends)!
        ! calculate irrigation heights at use mode
        type(grid_i),intent(in)::domain, irr_units_map, method
        integer,dimension(:,:),intent(in)::irr_class!
        type(irr_units_table),dimension(:),intent(inout)::irr_units!
        real(dp),dimension(:,:),intent(in)::h_old                      ! water content of 2nd layer of previous day [mm]
        real(dp),dimension(:,:),intent(in)::transp_pot                 ! potential transpiration of previous day [mm]
        real(dp),dimension(:,:),intent(in)::h_raw, h_raw_half          ! soil water content at RAW and half RAW [mm]
        real(dp),dimension(:,:),intent(in)::h_raw_coll                 ! irrigation threshold for collective water sources [mm]
        real(dp),dimension(:,:),intent(in)::h_raw_priv                 ! irrigation threshold for private water sources [mm]
        real(dp),dimension(:,:,:),intent(out)::h_irr
        integer,intent(in)::doy
        real(dp),dimension(:,:),intent(inout)::priv_irr,coll_irr!
        integer,dimension(:,:),intent(inout)::day_from_irr!
        real(dp),dimension(:,:,:),intent(inout)::esp_perc!
        type(grid_r),dimension(:),intent(in)::am_perc!
        type(grid_r),dimension(:),intent(in)::bm_perc!
        real(dp),dimension(:,:),intent(in)::cell_area
        real(dp),dimension(:,:),intent(in)::h_met!
        integer,dimension(:,:),intent(in)::irr_starts
        integer,dimension(:,:),intent(in)::irr_ends
        logical,intent(in)::f_shape_area
        
        integer::i,j,k,shift,n,p!
        integer,parameter::sec_to_day=24*60*60
        real(dp),dimension(domain%header%imax,domain%header%jmax)::v_irr_cell   ! irrigation volumes [m^3]
        logical,dimension(domain%header%imax,domain%header%jmax)::irr_mask      ! a mask to get all the irrigable cells
        integer,dimension(:,:),allocatable::i_mat,j_mat,id_cell                 ! for movement inside the matrix
        integer,dimension(size(irr_units))::ind!
        integer,dimension(:),allocatable::vi,vj,vid
        integer,dimension(:),allocatable::vcells ! sign the cells already processed
        ! shifted copies of the already defined variable (see above)
        real(dp),dimension(:),allocatable::s_h_old, s_h_raw_coll, s_h_raw, s_h_raw_half, &
                                            & s_h_raw_priv, s_h_transp_pot, s_transp_ratio, &
                                            & s_def_day, s_v_irr_cell, s_h_met
        real(dp)::n_day_to_deficit ! expected number of days to deficit
        integer::dist
        real(dp)::Q_tot_act,Q_tot_pot,q_cell,q_mean!
        
        ! rice cells are not considered in a different way as theta is already corrected
        
        ! TODO: move the initialization of the cell_area outside in order to overcome control
        if (f_shape_area .eqv. .false.) then
            v_irr_cell=1.e-3*h_met*(domain%header%cellsize**2) ! TODO: also NAN are calculated
        else
            v_irr_cell=1.e-3*h_met*cell_area
        end if
        
        h_irr = 0.
        
        ! init the matrix for coordinates and id of the calculation cells
        ! TODO: move the initialization of the cell_area outside in order to overcome control
        ind=0!
        if(.not.(allocated(i_mat))) allocate(i_mat(domain%header%imax,domain%header%jmax))!
        if(.not.(allocated(j_mat))) allocate(j_mat(domain%header%imax,domain%header%jmax))!
        if(.not.(allocated(id_cell))) allocate(id_cell(domain%header%imax,domain%header%jmax))!
        i_mat(:,1)=[(i,i=1,domain%header%imax)]!
        i_mat=spread(i_mat(:,1),dim=2,ncopies=domain%header%jmax)!
        j_mat(1,:)=[(j,j=1,domain%header%jmax)]!
        j_mat=spread(j_mat(1,:),dim=1,ncopies=domain%header%imax)!
        !
        id_cell=0!
        !
        do j=1,domain%header%jmax!
            do i=1,domain%header%imax!
                do k=1,size(irr_units)!
                    if(irr_units_map%mat(i,j)==irr_units(k)%id .and. domain%mat(i,j)/=domain%header%nan)then
                        ind(k)=ind(k)+1!
                        id_cell(i,j)=ind(k)!
                    end if!
                end do!
            end do!
        end do!

        irr_units_loop: do k=1,size(irr_units)!
            ! select the cells that belong to the k-irr_units and are irrigable
            ! the mask includes also the rice cells
            !irr_mask=(irr_units_map%mat==irr_units(k)%id .and. irr_class==1) ! OLD CONTROL
            ! %EAC% add control for the local irrigation season
            irr_mask=(irr_units_map%mat==irr_units(k)%id &
                     .and. irr_class==1 &
                     .and. irr_starts<=doy &
                     .and. irr_ends>=doy)

            n=count(irr_mask)!
            
            if (n==0) then
                !%AB% skyp the cycle if count(mask) == 0, no irrigable cells
                irr_units(k)%q_trashed = irr_units(k)%q_day
                irr_units(k)%q_day = 0
            else 
                ! %AB% in case of irrigable cells
                ! %AB% vi(n) e vj(n) are the coordinates of the irrigable cells
                allocate(vi(n));vi=pack(i_mat,irr_mask)!
                allocate(vj(n));vj=pack(j_mat,irr_mask)!
                allocate(vid(n));vid=pack(id_cell,irr_mask)!
                ! the following are all the required parameters
                allocate(s_h_old(n));s_h_old=pack(h_old,irr_mask)!
                allocate(s_h_raw_coll(n));s_h_raw_coll=pack(h_raw_coll,irr_mask)!
                allocate(s_h_raw(n));s_h_raw=pack(h_raw,irr_mask)!
                allocate(s_h_raw_half(n));s_h_raw_half=pack(h_raw_half,irr_mask)!
                allocate(s_h_raw_priv(n));s_h_raw_priv=pack(h_raw_priv,irr_mask)!
                allocate(s_h_transp_pot(n));s_h_transp_pot=pack(transp_pot,irr_mask)!
                allocate(vcells(n));vcells=0
                allocate(s_v_irr_cell(n));s_v_irr_cell=pack(v_irr_cell,irr_mask)!
                allocate(s_h_met(n));s_h_met=pack(h_met,irr_mask)!
                
                ! find the index of the latest irrigated cell (the day before)
                ! get_value_index returns shift=0 by default
                if (any(vid==irr_units(k)%last_cell_id)) then
                    shift=get_value_index(vid,irr_units(k)%last_cell_id)!
                end if!
                
                ! shift the vectors
                vi=cshift(vi,shift)!
                vj=cshift(vj,shift)!
                vid=cshift(vid,shift)    !
                s_h_old=cshift(s_h_old,shift)!
                s_h_raw_coll=cshift(s_h_raw_coll,shift)!
                s_h_raw=cshift(s_h_raw,shift)!
                s_h_raw_half=cshift(s_h_raw_half,shift)!
                s_h_raw_priv=cshift(s_h_raw_priv,shift)    !
                s_h_transp_pot=cshift(s_h_transp_pot,shift)!
                s_v_irr_cell=cshift(s_v_irr_cell,shift)!
                s_h_met=cshift(s_h_met,shift)!
                
                Q_tot_act=0. ! total discharge actually delivered for irrigation
                Q_tot_pot=0. ! total discharge potentially delivered for irrigation
                
                ! TODO: discharge at cell not in irrigation unit
                ! average net discharge required by the irrigation unit
                ! it doesn't consider the field efficiency
                ! it is estimated from the water depth from the irrigation method
                ! TODO: not used q_mean!
                q_mean = sum(s_v_irr_cell)/(n*sec_to_day) 
                
                ! TODOs:
                ! %AB% check irrigation when k_cb > 0
                ! %AB% check condition when the number of irrigable cells is very
                ! low respect the total number of cell in the irrigation unit
                
                cell_loop: do p=1,n !%AB% loop in irrigable cells where mask = 1
                    n_day_to_deficit=9999. !default (>0) x transp_pot=0!
                    
                    if(s_h_transp_pot(p)/=0.) n_day_to_deficit = (s_h_old(p)-s_h_raw_half(p))/s_h_transp_pot(p)
                    
                    ! discharge assigned to p-cell (not considering the field efficiency)
                    q_cell = s_v_irr_cell(p)/sec_to_day   
                    
                    ! check if the water is enough for irrigation                    
                    ! %AB% if there is enough water to irrigate the current cell
                    ! %AB% and the delivered can irrigate the current cell
                    if((Q_tot_pot+q_cell) < (irr_units(k)%q_day * irr_units(k)%f_explore) .and. & 
                        & (Q_tot_act+q_cell) <= irr_units(k)%q_day) then                           
                
                        ! %AB% record the index of the current cell.
                        ! It could be the latest irrigated in the current day and the first in the following 
                        Q_tot_pot = Q_tot_pot + q_cell
                        irr_units(k)%last_cell_id = vid(p)  
                        vcells(p) = vid(p)!
                    
                        ! check if the cell is irrigable according to the soil water content less than the RAW big threshold
                        ! %AB% alternative consider the number of days to have deficit 
                        if(s_h_old(p)<=s_h_raw_coll(p)) then !.or. dgg <= (n/IU(k)%n_irrigable_cells))then
                            ! %AB% apply the irrigation depth of the method
                            ! %AB% and sum the discharge actually delivered to the cell
                            coll_irr(vi(p),vj(p)) = s_h_met(p)
                            Q_tot_act = Q_tot_act + q_cell
                        end if!
                    else
                        exit cell_loop
                    end if
                end do cell_loop
                !
                ! %AB% update the irrigation period
                if (n /=0 .and. p<=n .and. Q_tot_act /= 0) then
                    irr_units(k)%n_irrigated_cells = p
                else if (p>n) then ! all the cells are verified
                    irr_units(k)%n_irrigated_cells = p-1
                else ! irrigation units with no irrigable cells
                    irr_units(k)%n_irrigated_cells = 0
                end if
            
                ! mass conservation: how much water remains
                irr_units(k)%q_trashed = irr_units(k)%q_day - Q_tot_act
                
                ! %AB% store the amount of water surplus for the next day if:
                ! 1 - the irrigation unit has irrigable cells
                ! 2 - not all the cells are evaluated for irrigation requirements
                ! 3 - current day is in the irrigation season
                !if (n/=0 .and. irr_units(k)%n_irrigated_cells /= n .and. doy < end_irr_season) then
                if (n/=0 .and. irr_units(k)%n_irrigated_cells /= n) then
                    irr_units(k)%q_surplus = irr_units(k)%q_trashed
                    irr_units(k)%q_trashed = 0.  ! %AB% re-init the surplus
                else
                    irr_units(k)%q_surplus = 0.
                end if

                !%AB% store the water actually used
                irr_units(k)%q_day=Q_tot_act
                                
                ! UNMONITORED PRIVATE ONLY
                ! %AB% restart n 
                n = count(irr_mask)!
                
                if (irr_units(k)%f_un_priv==1) then!
                    allocate(s_transp_ratio(n)); s_transp_ratio = s_h_transp_pot/sum(s_h_transp_pot)!
                    allocate(s_def_day(n)); s_def_day = 0.!
                    where(s_h_transp_pot/=0.) s_def_day=(s_h_old-s_h_raw)/s_h_transp_pot!
                    ! %AB% skip at least one cell that could be irrigated by monitored sources
                    ! dist = int(sum(vtmm-vRaw)/sum(vTrasp_pot))*IU(k)%n_irrigable_cells ! OLD version
                    dist = 1 + nint((sum(s_transp_ratio*s_def_day)/sum(s_transp_ratio))*irr_units(k)%n_irrigable_cells) 
                    
                    
                    if(dist>n) dist=n!
                    ! TODO: check shift is the same as before
                    ! vector shift
                    vi=cshift(vi,shift)!
                    vj=cshift(vj,shift)!
                    vid=cshift(vid,shift)    !
                    s_h_old=cshift(s_h_old,shift)!
                    s_h_raw_coll=cshift(s_h_raw_coll,shift)!
                    s_h_raw=cshift(s_h_raw,shift)!
                    s_h_raw_half=cshift(s_h_raw_half,shift)!
                    s_h_raw_priv=cshift(s_h_raw_priv,shift)    !
                    s_h_transp_pot=cshift(s_h_transp_pot,shift)!
                    s_v_irr_cell=cshift(s_v_irr_cell,shift)!
                    s_h_met=cshift(s_h_met,shift)
                    
                    do p=1,n!
                        if((p>=dist).and.(.not.(any(vcells==vid(p)))))then!
                            if(s_h_old(p) <= s_h_raw_priv(p)) then!
                                ! %AB% irrigation volume is equal to the irrigation depth of the method
                                priv_irr(vi(p),vj(p))=s_h_met(p)                    
                                irr_units(k)%q_un_priv=irr_units(k)%q_un_priv+s_v_irr_cell(p)/sec_to_day
                            end if!
                        end if!
                    end do!
                    
                    deallocate(s_transp_ratio)!
                    deallocate(s_def_day)!
                end if
                
                ! free memory
                deallocate(vi)
                deallocate(vj)
                deallocate(vid)
                deallocate(s_h_old)
                deallocate(s_h_raw_coll)
                deallocate(s_h_raw)
                deallocate(s_h_raw_half)
                deallocate(s_h_raw_priv)
                deallocate(s_h_transp_pot)
                deallocate(vcells)
                deallocate(s_v_irr_cell)
                deallocate(s_h_met)
            end if
        end do irr_units_loop ! end irrigation units loop
        !!
        ! check coll_irr(i,j)/=priv_irr(i,j) e record irrigated volumes!
        ! TODO: split irrigation sources
        do j=1,domain%header%jmax!
            do i=1,domain%header%imax!
                if(domain%mat(i,j)/=domain%header%nan)then!
                    if(.not.(priv_irr(i,j)/=0. .and. coll_irr(i,j)/=0.))then!
                        h_irr(i,j,method%mat(i,j)) = priv_irr(i,j)+coll_irr(i,j)!
                    else!
                        print*,"Some cells are irrigated both from wells and diversions"    !
                    end if!
                end if!
            end do!
        end do!
        !!
        
        ! TODO: replace with update_perco_parameters subroutine
        ! calculate the exponential param from the number of days from the latest irrigation event day
        where(priv_irr+coll_irr==0)!
            where(day_from_irr==-9999)!
                day_from_irr=-9999!
            else where!
                day_from_irr=day_from_irr+1!
            end where!
        else where!
            day_from_irr=0!
        end where!

        esp_perc(:,:,1)=merge(1.0D0,1+am_perc(1)%mat*exp(-1.0*day_from_irr*bm_perc(1)%mat),day_from_irr==-9999)!
        esp_perc(:,:,2)=merge(1.0D0,1+am_perc(2)%mat*exp(-1.0*day_from_irr*bm_perc(2)%mat),day_from_irr==-9999)!
        
    end subroutine irrigation_use!

    subroutine calc_daily_duty(cur_doy,irr_units,sources_info,wat_sources,irr_units_map,domain_map,par,&!
        & irrigation_class,h_soil,h_transp_pot,raw,h_fc,zr)!
        ! estimate the water volume for irrigation available in each irrigation units 
        integer,intent(in)::cur_doy ! current day of simulation
        type(irr_units_table), dimension(:),intent(inout)::irr_units!
        type(source_info),intent(in)::sources_info!
        type(water_sources_table),dimension(:),intent(in)::wat_sources!
        type(grid_i),intent(in)::irr_units_map,domain_map!
        type(parameters),intent(in)::par!
        integer,dimension(:,:),intent(in)::irrigation_class!
        real(dp),dimension(:,:),intent(in)::h_soil, h_transp_pot, raw,h_fc, zr!
        !!
        integer,dimension(:,:),allocatable::cells_un_coll       ! map of the cells irrigated by unmonitored collective water sources
        real(dp)::frac_rel_un_coll                              ! fraction of water released respect to the maximum water available
        integer:: n_cells_irr_un_coll                           ! number of cells irrigated by unmonitored collective water sources
        integer:: n_cells_un_coll                               ! number of cells irrigable by unmonitored collective water sources
        real(dp),dimension(par%cr%n_withdrawals)::q_un_coll!
        integer::i,j,k!
        !!
        ! estimate the daily duty for each irrigation units
        irr_units%q_act_fld(1)=0.;irr_units%q_act_fld(2)=0.;irr_units%q_act_fld(3)=0.;irr_units%q_act_fld(4)=0.;!
        irr_units%q_day=0.!
        irr_units%q_trashed=0.
        irr_units%n_irrigated_cells=0.
        irr_units%q_un_priv=0.
        frac_rel_un_coll = 0.
        q_un_coll=0.!
        !
        if(.not.(allocated(cells_un_coll))) allocate(cells_un_coll(domain_map%header%imax, domain_map%header%jmax))
        cells_un_coll = 0
        if (par%cr%f_exists .eqv. .true.) then
            do i=1, size(wat_sources)
                ! %AB% irrigation_class is updated every days, so this cycle must be updated every day
                if (wat_sources(i)%type_id == 4) then
                    where (irr_units_map%mat == wat_sources(i)%id_irr_unit &
                        &.and. domain_map%mat /= domain_map%header%nan &
                        & .and. irrigation_class == 1) &
                        & cells_un_coll = wat_sources(i)%wat_src_idx
                end if
            end do
        end if
        !
        do i=1,size(wat_sources)!
            j=wat_sources(i)%irr_unit_idx;
            if (j==0) then
                print*,'WARNING: water source ',wat_sources(i)%id_wat_src,' returns zero index'
                cycle
            end if
            k=wat_sources(i)%wat_src_idx!
            select case(wat_sources(i)%type_id)!
                case(1)!
                    irr_units(j)%q_act_fld(1) = wat_sources(i)%duty_frc * sources_info%mn_src_tbl1%q_daily(cur_doy,k) * irr_units(j)%int_distr_eff + irr_units(j)%q_act_fld(1)
                case(2)!
                    irr_units(j)%q_act_fld(2) = wat_sources(i)%duty_frc * sources_info%mn_src_tbl2%q_daily(cur_doy,k) * irr_units(j)%int_distr_eff + irr_units(j)%q_act_fld(2)
                case(3)!
                    irr_units(j)%q_act_fld(3) = wat_sources(i)%duty_frc * sources_info%int_reuse_tbl%q_daily(cur_doy,k) * irr_units(j)%int_distr_eff + irr_units(j)%q_act_fld(3)
                case(4)!
                    frac_rel_un_coll = deliverable_ratio_unm_coll(cells_un_coll,h_soil,h_transp_pot,raw,h_fc,zr,irrigation_class,par, &
                                        & domain_map%header%imax, domain_map%header%jmax, k)!
                    n_cells_irr_un_coll = count(cells_un_coll==k .and. irrigation_class==1)!
                    n_cells_un_coll = count(cells_un_coll==k)
                    if(n_cells_un_coll==0) cycle
                    irr_units(j)%q_act_fld(4)= &
                        & frac_rel_un_coll*wat_sources(i)%duty_frc * sources_info%unm_src_tbl%q_max(k) * irr_units(j)%int_distr_eff * (real(n_cells_irr_un_coll)/real(n_cells_un_coll)) + &
                        & irr_units(j)%q_act_fld(4)!
                    
                    q_un_coll(k) = &
                        & frac_rel_un_coll*wat_sources(i)%duty_frc * sources_info%unm_src_tbl%q_max(k)*(real(n_cells_irr_un_coll)/real(n_cells_un_coll))+&
                        & q_un_coll(k)   ! only for printing %AB% gross used water
                case default!
            end select!
        end do!

        ! compute the water duty of each irrigation unit
        irr_units%q_day=irr_units%q_act_fld(1)+irr_units%q_act_fld(2)+irr_units%q_act_fld(3)+irr_units%q_act_fld(4)+irr_units%q_surplus ! %AB% consider also the water surplus of the previous day
    end subroutine calc_daily_duty!

    function deliverable_ratio_unm_coll(cells_un_coll, h_soil, h_transp_pot, h_raw, theta_fc, &
                                & d_r, irr_class, pars, imax, jmax, wat_src_idx)!
        ! calculate the actual deliverable ratio from unmonitored water sources (e.g. collective wells plan)
        ! deliverable ratio is calculated from average water deficit of the served cells
        integer,intent(in)::imax,jmax!
        integer,dimension(:,:),intent(in)::cells_un_coll     ! map of cells irrigated by unmonitored collective sources
        real(dp),dimension(:,:),intent(in)::h_soil, h_transp_pot, h_raw, theta_fc, d_r
        integer,dimension(:,:),intent(in)::irr_class!
        type(parameters),intent(in)::pars
        integer,intent(in)::wat_src_idx

        real(dp)::transp_un_coll
        real(dp)::deficit_un_coll
        real(dp)::raw_un_coll
        real(dp),dimension(imax,jmax)::deficit_slp, raw_slp
        logical::deficit
        
        real(dp)::deliverable_ratio_unm_coll
        !!
        deliverable_ratio_unm_coll=0
        deficit_slp= (h_soil - h_raw) * h_transp_pot!
        raw_slp=(theta_fc*1000*d_r-h_raw)*h_transp_pot!
        ! %AB% for rice, wat%layer(2)%fc=1000*xriso%tetaII_FC*zr
        
        transp_un_coll=sum(h_transp_pot, cells_un_coll == wat_src_idx .and. irr_class==1)!

        if(transp_un_coll/=0.)then
            deficit_un_coll = sum(deficit_slp, cells_un_coll==wat_src_idx .and. irr_class==1)/transp_un_coll
            raw_un_coll=sum(raw_slp,cells_un_coll==wat_src_idx .and. irr_class==1)/transp_un_coll!
        else!
            deficit_un_coll=0.
            raw_un_coll=0.
        end if!
        
        deficit = deficit_un_coll < (pars%uc_act_rules(wat_src_idx)%min_act_trsld*raw_un_coll)!
        
        if(deficit .eqv. .true.)then!
            if(deficit_un_coll <= pars%uc_act_rules(wat_src_idx)%flow_rate(3)*raw_un_coll) then!
                deliverable_ratio_unm_coll = pars%uc_act_rules(wat_src_idx)%act_trsld(3)!
                return!
            else if(deficit_un_coll <= pars%uc_act_rules(wat_src_idx)%flow_rate(2)*raw_un_coll)then!
                deliverable_ratio_unm_coll = pars%uc_act_rules(wat_src_idx)%act_trsld(2)!
                return!
            else if(deficit_un_coll <= pars%uc_act_rules(wat_src_idx)%flow_rate(1)*raw_un_coll)then!
                deliverable_ratio_unm_coll = pars%uc_act_rules(wat_src_idx)%act_trsld(1)!
                return!
            end if!
        end if!
       !
    end function deliverable_ratio_unm_coll!

end module