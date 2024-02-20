module mod_crop_phenology
    use mod_grid
    implicit none

    type file_phenology_r!
        integer::unit                           ! unit associated to the file
        real(dp),dimension(:,:),pointer::tab    ! daily crop parameters for each land uses
    end type file_phenology_r!
    !
    type file_phenology_i!
        integer::unit                           ! unit associated to the file
        integer,dimension(:,:),pointer::tab     ! daily crop parameters for each land uses
    end type file_phenology_i!
    
    type k_cb_matrices
        real(dp), dimension(:,:), pointer::low  ! minimum kcb
        real(dp), dimension(:,:), pointer::high ! maximum kcb
        real(dp), dimension(:,:), pointer::mid  ! kcb between the 1st and 2nd stage
    end type k_cb_matrices

    type crop_pheno_info!
        ! include all crop parameters for each calculation cell
        integer,dimension(:,:),pointer::ii0                 ! start growing day of the crop [doy]
        integer,dimension(:,:),pointer::iid                 ! duration of the growing period [day]
        integer,dimension(:,:),pointer::iie                 ! harvest day [doy]
        integer,dimension(:),pointer::n_crops_by_year       ! number of crop in a year
        type(file_phenology_r)::k_cb                        ! base crop coefficient
        type(file_phenology_r)::h                           ! crop height
        type(file_phenology_r)::d_r                         ! root depth
        type(file_phenology_r)::lai                         ! leaf area index
        type(file_phenology_r)::f_c                         ! cover fraction
        type(file_phenology_r)::p_raw                       ! EDIT: daily variable readly available water factor [-]
        type(file_phenology_i)::cn_day                      ! %AB% soil moisture adjusted cn value 
        type(k_cb_matrices)::kcb_phases                     ! k_cb change points during phenology
        real(dp),dimension(:,:),pointer::p_raw_const        ! readly available water factor [-]
        real(dp),dimension(:,:),pointer::a                  ! interception coefficient according to the Von Hoyningen-Hune & Braden method [-]
        real(dp),dimension(:,:),pointer::max_d_r            ! maximum root depth [m]
        real(dp),dimension(:,:),pointer::max_RF_t           ! maximum fractio of active roots in transpirative layer [-]
        real(dp),dimension(:,:),pointer::T_lim              ! limit temperature threshold - termic stress [°C]
        real(dp),dimension(:,:),pointer::T_crit             ! critical temperature threshold - termic stress [°C]
        real(dp),dimension(:,:),pointer::HI                 ! harvest index
        real(dp),dimension(:,:),pointer::Ky_tot             ! water stress coefficient - overall
        real(dp),dimension(:,:,:),pointer::Ky_pheno         ! water stress coefficient - phases
        integer,dimension(:,:),pointer::irrigation_class    ! 1 = irrigated, 0 = not irrgated
        integer,dimension(:,:),pointer::cn_class            ! CN class
        real(dp),dimension(:,:,:),pointer::wp_adj            ! normalized biomass water productivity 
    end type crop_pheno_info   !es: valday(:)%kcb%unit; valday(:)%kcb%tab(:,:)!

    type crop_pars_matrices!
        ! store the crop parameters for each calculation cells
        ! see crop_pheno_info for details
        real(dp),dimension(:,:),pointer::k_cb
        real(dp),dimension(:,:),pointer::h
        real(dp),dimension(:,:),pointer::d_r
        real(dp),dimension(:,:),pointer::lai
        integer,dimension(:,:),pointer::cn_day
        real(dp),dimension(:,:),pointer::f_c
        integer,dimension(:,:),pointer::irrigation_class
        integer,dimension(:,:),pointer::cn_class
        real(dp),dimension(:,:),pointer::p
        real(dp),dimension(:,:),pointer::a
        real(dp),dimension(:,:),pointer::max_d_r
        real(dp),dimension(:,:),pointer::max_RF_t
        real(dp),dimension(:,:),pointer::RF_e
        real(dp),dimension(:,:),pointer::RF_t
        real(dp),dimension(:,:),pointer::T_lim
        real(dp),dimension(:,:),pointer::T_crit
        real(dp),dimension(:,:),pointer::HI
        real(dp),dimension(:,:),pointer::Ky_tot
        real(dp),dimension(:,:,:),pointer::Ky_pheno
        real(dp),dimension(:,:),pointer::k_cb_low
        real(dp),dimension(:,:),pointer::k_cb_mid
        real(dp),dimension(:,:),pointer::k_cb_high
        real(dp),dimension(:,:),pointer::wp_adj
        real(dp),dimension(:,:),pointer::p_day
        real(dp),dimension(:,:),pointer::k_cb_old           ! k_cb of previous day
        integer,dimension(:,:),pointer::n_crops_by_year
        integer,dimension(:,:),pointer::pheno_idx           ! phenological stage index
    end type crop_pars_matrices!

    type crop_matrices
    ! crop_pheno_info for details    
        integer,dimension(:,:,:),pointer::ii0
        integer,dimension(:,:,:),pointer::iie
        real(dp),dimension(:,:,:),pointer::iid
        integer,dimension(:,:,:),pointer::ii0_ref           ! harvest date in meteorological reference series
        integer,dimension(:,:,:),pointer::iie_ref           ! harvest date in meteorological reference series
        real(dp),dimension(:,:,:),pointer::iid_ref          ! crop cycle length in meteorological reference series
        real(dp),dimension(:,:,:),pointer::dij              ! coefficient of crop cycle expansion in relation to meteorological reference series
        integer,dimension(:,:,:),pointer::TSP_high
        integer,dimension(:,:,:),pointer::TSP_low
        real(dp),dimension(:,:,:),pointer::wp_adj
        real(dp),dimension(:,:,:),pointer::HI
        real(dp),dimension(:,:,:),pointer::Ky_tot
        real(dp),dimension(:,:,:,:),pointer::Ky_pheno
        real(dp),dimension(:,:,:),pointer::T_crit
        real(dp),dimension(:,:,:),pointer::T_lim
        real(dp),dimension(:,:,:),pointer::k_cb_min
        real(dp),dimension(:,:,:),pointer::k_cb_mid
        real(dp),dimension(:,:,:),pointer::k_cb_max
    end type crop_matrices

    contains

    subroutine make_random_emergence(info_pheno,meteo_weight,dir_meteo,domain,soiluse,crop_mat,irandom,year_length)!
        !randomization of emergence date from phenological series at meteorological stations

        implicit none!
        real(dp),dimension(:,:,:),intent(in)::meteo_weight!
        integer,dimension(:,:,:),intent(in)::dir_meteo!
        type(grid_i),intent(in)::domain!
        integer,dimension(:,:),intent(in)::soiluse!
        type(crop_pheno_info),dimension(:),intent(in)::info_pheno!
        type(crop_matrices),intent(out)::crop_mat
        integer,dimension(:,:),intent(in)::irandom
        integer,intent(in)::year_length
        
        real(dp),dimension(size(crop_mat%ii0,1),size(crop_mat%ii0,2),size(crop_mat%ii0,3))::ii0_r
        real(dp),dimension(size(crop_mat%iie,1),size(crop_mat%iie,2),size(crop_mat%iie,3))::iie_r
        integer::i,j,k,z
        
        ! initialization
        ii0_r = 0.d0
        iie_r = 0.d0
        
        ! spatial distribution of emergence date and phenological cycle duration by using meteorological weights
        do k=1,size(meteo_weight,3)!
            do j=1,size(domain%mat,2)!
                do i=1,size(domain%mat,1)!
                    if(domain%mat(i,j)/=domain%header%nan)then!
                        ii0_r(i,j,:)   =info_pheno(dir_meteo(i,j,k))%ii0(soiluse(i,j),:)*meteo_weight(i,j,k) + ii0_r(i,j,:)
                        iie_r(i,j,:)   =info_pheno(dir_meteo(i,j,k))%iie(soiluse(i,j),:)*meteo_weight(i,j,k) + iie_r(i,j,:)
                        crop_mat%iid(i,j,:)=info_pheno(dir_meteo(i,j,k))%iid(soiluse(i,j),:)*meteo_weight(i,j,k) + crop_mat%iid(i,j,:)  ! crop cycle length (real)
                    end if!
                end do!                
            end do!
        end do!

        do z=1,size(crop_mat%ii0,3)
            crop_mat%ii0(:,:,z)=merge(nint(ii0_r(:,:,z)),domain%header%nan,domain%mat/=domain%header%nan)   ! emergence date
            crop_mat%iie(:,:,z)=merge(nint(iie_r(:,:,z)),domain%header%nan,domain%mat/=domain%header%nan)   ! harvest date
        end do
        
        do j=1,size(domain%mat,2)!
            do i=1,size(domain%mat,1)!
                if(domain%mat(i,j)/=domain%header%nan)then!
                    do z=1,size(crop_mat%ii0,3)
                        crop_mat%ii0_ref(i,j,:) = info_pheno(dir_meteo(i,j,1))%ii0(soiluse(i,j),:)   ! emergence date in meteorological reference series
                        crop_mat%iie_ref(i,j,:) = info_pheno(dir_meteo(i,j,1))%iie(soiluse(i,j),:)   ! harvest date in meteorological reference series
                        crop_mat%iid_ref(i,j,:) = info_pheno(dir_meteo(i,j,1))%iid(soiluse(i,j),:)   ! crop cycle length in meteorological reference series
                        crop_mat%dij(i,j,:)     = crop_mat%iid_ref(i,j,:)/crop_mat%iid(i,j,:)        ! coefficient of crop cycle expansion in relation to meteorological reference series
                    end do
                end if
            end do
        end do
        
        do z=1, size(crop_mat%ii0,3)
            crop_mat%ii0_ref(:,:,z) = merge(crop_mat%ii0_ref(:,:,z),domain%header%nan,domain%mat/=domain%header%nan)
            crop_mat%ii0_ref(:,:,z) = merge(crop_mat%ii0_ref(:,:,z),domain%header%nan,crop_mat%ii0_ref(:,:,z)/=0)
            crop_mat%iie_ref(:,:,z) = merge(crop_mat%iie_ref(:,:,z),domain%header%nan,domain%mat/=domain%header%nan)
            crop_mat%iie_ref(:,:,z) = merge(crop_mat%iie_ref(:,:,z),domain%header%nan,crop_mat%iie_ref(:,:,z)/=0)
            crop_mat%iid_ref(:,:,z) = merge(crop_mat%iid_ref(:,:,z),dble(domain%header%nan),domain%mat/=domain%header%nan)
            crop_mat%iid_ref(:,:,z) = merge(crop_mat%iid_ref(:,:,z),dble(domain%header%nan),crop_mat%iid_ref(:,:,z)/=0)
            crop_mat%dij(:,:,z)     = merge(crop_mat%dij(:,:,z),dble(domain%header%nan),domain%mat/=domain%header%nan)
            crop_mat%dij(:,:,z)     = merge(crop_mat%dij(:,:,z),dble(domain%header%nan),crop_mat%iid(:,:,z)/=0)
        end do
        
        do z=1,size(crop_mat%ii0,3)
            crop_mat%TSP_high(:,:,z) = nint(0.75 * crop_mat%iid(:,:,z)) + crop_mat%ii0(:,:,z) + irandom
            crop_mat%TSP_low(:,:,z)  = nint(0.45 * crop_mat%iid(:,:,z)) + crop_mat%ii0(:,:,z) + irandom
        end do
          
        ! values adjustment for double-years crops
        crop_mat%TSP_high = merge(crop_mat%TSP_high-year_length, crop_mat%TSP_high, crop_mat%TSP_high>year_length)
        crop_mat%TSP_low  = merge(crop_mat%TSP_low-year_length,  crop_mat%TSP_low,  crop_mat%TSP_low>year_length)
        
    end subroutine make_random_emergence!

    subroutine populate_crop_yield_matrices(info_pheno,dir_phenofases,domain,soil_use,crop_mat,year)
        ! populate of crop matrices with crop yield parameters
        implicit none!
        type(crop_pheno_info),dimension(:),intent(in)::info_pheno
        integer,dimension(:,:),intent(in)::dir_phenofases!
        type(grid_i),intent(in)::domain
        integer,dimension(:,:),intent(in)::soil_use
        type(crop_matrices),intent(inout)::crop_mat
        integer,intent(in)::year
        integer::i,j
        
        do j=1,size(domain%mat,2)!
            do i=1,size(domain%mat,1)!
                if(domain%mat(i,j)/=domain%header%nan)then!
                    crop_mat%wp_adj(i,j,:) = info_pheno(dir_phenofases(i,j))%wp_adj(soil_use(i,j),:,year)
                    crop_mat%HI(i,j,:) = info_pheno(dir_phenofases(i,j))%HI(soil_use(i,j),:)
                    crop_mat%Ky_tot(i,j,:) = info_pheno(dir_phenofases(i,j))%Ky_tot(soil_use(i,j),:)
                    crop_mat%Ky_pheno(i,j,:,:) = info_pheno(dir_phenofases(i,j))%Ky_pheno(soil_use(i,j),:,:)
                    crop_mat%T_crit(i,j,:) = info_pheno(dir_phenofases(i,j))%T_crit(soil_use(i,j),:)
                    crop_mat%T_lim(i,j,:) = info_pheno(dir_phenofases(i,j))%T_lim(soil_use(i,j),:)
                    crop_mat%k_cb_min(i,j,:) = info_pheno(dir_phenofases(i,j))%kcb_phases%low(soil_use(i,j),:)
                    crop_mat%k_cb_mid(i,j,:) = info_pheno(dir_phenofases(i,j))%kcb_phases%mid(soil_use(i,j),:)
                    crop_mat%k_cb_max(i,j,:) = info_pheno(dir_phenofases(i,j))%kcb_phases%high(soil_use(i,j),:)
                end if
            end do
        end do

    end subroutine populate_crop_yield_matrices

    subroutine populate_crop_pars_matrices(crop_pars_mat,info_pheno,irandom,doy,ws_idx,domain,soil_use,y, year_length, crop_mat)!
        ! populate crop parameters matrices from weather stations time series
        integer,intent(in)::doy,year_length,y!
        type(grid_i),intent(in)::domain,soil_use!
        integer,dimension(:,:),intent(in)::ws_idx                   !index in the list of weather stations
        type(crop_pheno_info),dimension(:),intent(in)::info_pheno
        type(crop_pars_matrices),intent(inout)::crop_pars_mat
        integer,dimension(:,:),intent(inout)::irandom               ! pseudorandom parameter that shifts crop cycle
        type(crop_matrices),intent(in)::crop_mat

        integer::i,j    
        integer::doy_s ! shifted day of the year 
        
        do j=1,size(domain%mat,2)!
            do i=1,size(domain%mat,1)!
                scans_domain: if(domain%mat(i,j) /= domain%header%nan)then!
                    if (crop_mat%ii0(i,j,crop_pars_mat%n_crops_by_year(i,j)) == 0 .and. crop_mat%iie(i,j,crop_pars_mat%n_crops_by_year(i,j)) == 0 ) then  ! no crop
                        doy_s = doy - irandom(i,j)
                    else if (crop_mat%ii0(i,j,crop_pars_mat%n_crops_by_year(i,j)) < crop_mat%iie(i,j,crop_pars_mat%n_crops_by_year(i,j))) then            ! annuals or perennials
                        ! emergence date is shifted as ii0(i,j,cs)-irandom(i,j,cs)
                        ! nint((gg-ii0(i,j))*dij(i,j)) contracts/expands the series
                        ! randomization of emergence date (ii0/irandom) and factor of dilatation (dij) are used to calculate gg1
                        doy_s = crop_mat%ii0_ref(i,j,crop_pars_mat%n_crops_by_year(i,j)) - irandom(i,j) + &
                            & nint((doy-crop_mat%ii0(i,j,crop_pars_mat%n_crops_by_year(i,j))) * crop_mat%dij(i,j,crop_pars_mat%n_crops_by_year(i,j)))
                    else                                                                                                ! biennals
                        if (doy < crop_mat%iie(i,j,crop_pars_mat%n_crops_by_year(i,j)) + &
                            & irandom(i,j)) then
                            ! from 1/1 to harvest date, only the contraction/expansion of crop cycle is taken into account
                            ! in the first part of the year, the limits are 1/1 and harvest date (only contraction/expansion)
                            doy_s = doy * (crop_mat%iie_ref(i,j,crop_pars_mat%n_crops_by_year(i,j))) / &
								& crop_mat%iie(i,j,crop_pars_mat%n_crops_by_year(i,j)) - irandom(i,j)
                        else if (doy >= crop_mat%ii0(i,j,crop_pars_mat%n_crops_by_year(i,j)) + &
                            & irandom(i,j)) then
                            
                            ! from emergence to 31/12, both parameters are taken into account: the randomization of emergence date and the contraction/expansion of crop cycle
                            ! in the second part of the year, the limits are the randomized emergence date and 31/12 (day 365/366)
                            doy_s = (doy - crop_mat%ii0(i,j,crop_pars_mat%n_crops_by_year(i,j))) * &
								& (year_length - crop_mat%ii0_ref(i,j,crop_pars_mat%n_crops_by_year(i,j))) / &
                                & (year_length - crop_mat%ii0(i,j,crop_pars_mat%n_crops_by_year(i,j))) + &
                                & crop_mat%ii0_ref(i,j,crop_pars_mat%n_crops_by_year(i,j)) - irandom (i,j)
                        else
                            doy_s = crop_mat%iie_ref(i,j,crop_pars_mat%n_crops_by_year(i,j)) +1    ! it points to a null Kcb
                        end if
                    end if
                    
                    ! adjust pointing if gg1 doesn't belong to [1,year_length] (it takes into account irandom and not agricultural soil uses)
                    if (doy_s < 1) doy_s = 1
                    if (doy_s > year_length) doy_s = year_length

                    ! TODO: parameters overloading can be moved to a specific subroutine
                    
                    ! conveniently updates phenological data from its series - update occurs only if Kcb varies
                    crop_pars_mat%k_cb(i,j)=info_pheno(ws_idx(i,j))%k_cb%tab(doy_s,soil_use%mat(i,j))
                    
                    if (crop_pars_mat%k_cb(i,j) /= crop_pars_mat%k_cb_low(i,j) .or. crop_pars_mat%k_cb_old(i,j) /= crop_pars_mat%k_cb_low(i,j)) then
                        if (crop_pars_mat%k_cb_old(i,j) > crop_pars_mat%k_cb_low(i,j) .and. crop_pars_mat%k_cb(i,j) == crop_pars_mat%k_cb_low(i,j) &
                            & .and. info_pheno(ws_idx(i,j))%n_crops_by_year(soil_use%mat(i,j))>1) then
                            if (crop_pars_mat%n_crops_by_year(i,j) < info_pheno(ws_idx(i,j))%n_crops_by_year(soil_use%mat(i,j))) then       ! cult_switch cycle
                                crop_pars_mat%n_crops_by_year(i, j) = crop_pars_mat%n_crops_by_year(i, j) + 1
                                crop_pars_mat%pheno_idx(i,j) = 1
                            else if (crop_pars_mat%n_crops_by_year(i,j) == info_pheno(ws_idx(i,j))%n_crops_by_year(soil_use%mat(i,j))) then ! if cult_switch cycle ends, switch is set back to 1
                                crop_pars_mat%n_crops_by_year(i, j) = 1
                                crop_pars_mat%pheno_idx(i,j) = 1
                            end if
                        end if
                        crop_pars_mat%h(i,j)=info_pheno(ws_idx(i,j))%h%tab(doy_s,soil_use%mat(i,j))
                        crop_pars_mat%d_r(i,j)=info_pheno(ws_idx(i,j))%d_r%tab(doy_s,soil_use%mat(i,j))
                        crop_pars_mat%lai(i,j)=info_pheno(ws_idx(i,j))%lai%tab(doy_s,soil_use%mat(i,j))
                        crop_pars_mat%cn_day(i,j)=info_pheno(ws_idx(i,j))%cn_day%tab(doy_s,soil_use%mat(i,j))
                        crop_pars_mat%f_c(i,j)=info_pheno(ws_idx(i,j))%f_c%tab(doy_s,soil_use%mat(i,j))

                       
                        crop_pars_mat%irrigation_class(i,j)=&
                            info_pheno(ws_idx(i,j))%irrigation_class(soil_use%mat(i,j),crop_pars_mat%n_crops_by_year(i,j))
                        crop_pars_mat%cn_class(i,j)=info_pheno(ws_idx(i,j))%cn_class(soil_use%mat(i,j),crop_pars_mat%n_crops_by_year(i,j))
                        ! TEST: replace here with variable p value
                        !crop_pars_mat%p(i,j)=info_pheno(ws_idx(i,j))%p_raw%tab(doy_s,soil_use%mat(i,j))
                        crop_pars_mat%p(i,j)=info_pheno(ws_idx(i,j))%p_raw_const(soil_use%mat(i,j),crop_pars_mat%n_crops_by_year(i,j))

                        crop_pars_mat%a(i,j)=info_pheno(ws_idx(i,j))%a(soil_use%mat(i,j),crop_pars_mat%n_crops_by_year(i,j))
                        crop_pars_mat%max_d_r(i,j)=info_pheno(ws_idx(i,j))%max_d_r(soil_use%mat(i,j),crop_pars_mat%n_crops_by_year(i,j))
                        crop_pars_mat%max_RF_t(i,j)=info_pheno(ws_idx(i,j))%max_RF_t(soil_use%mat(i,j),crop_pars_mat%n_crops_by_year(i,j))
                        crop_pars_mat%T_lim(i,j)=info_pheno(ws_idx(i,j))%T_lim(soil_use%mat(i,j),crop_pars_mat%n_crops_by_year(i,j))
                        crop_pars_mat%T_crit(i,j)=info_pheno(ws_idx(i,j))%T_crit(soil_use%mat(i,j),crop_pars_mat%n_crops_by_year(i,j))
                        crop_pars_mat%HI(i,j)=info_pheno(ws_idx(i,j))%HI(soil_use%mat(i,j),crop_pars_mat%n_crops_by_year(i,j))
                        crop_pars_mat%Ky_tot(i,j)=info_pheno(ws_idx(i,j))%Ky_tot(soil_use%mat(i,j),crop_pars_mat%n_crops_by_year(i,j))
                        crop_pars_mat%Ky_pheno(i,j,:)=info_pheno(ws_idx(i,j))%Ky_pheno(soil_use%mat(i,j),crop_pars_mat%n_crops_by_year(i,j),:)
                        crop_pars_mat%k_cb_low(i,j)=info_pheno(ws_idx(i,j))%kcb_phases%low(soil_use%mat(i,j),crop_pars_mat%n_crops_by_year(i,j))
                        crop_pars_mat%k_cb_mid(i,j)=info_pheno(ws_idx(i,j))%kcb_phases%mid(soil_use%mat(i,j),crop_pars_mat%n_crops_by_year(i,j))
                        crop_pars_mat%k_cb_high(i,j)=info_pheno(ws_idx(i,j))%kcb_phases%high(soil_use%mat(i,j),crop_pars_mat%n_crops_by_year(i,j))
                        crop_pars_mat%wp_adj(i,j) = info_pheno(ws_idx(i,j))%wp_adj(soil_use%mat(i,j),crop_pars_mat%n_crops_by_year(i,j),y)
                    end if
                end if scans_domain
            end do!
        end do!
    end subroutine populate_crop_pars_matrices!
    
    subroutine calculate_RF_t(d_r, d_e_fix, crop_par_mat,domain)!
        ! calculate root fraction in both evaporative and transpirative layer 
        real(dp), dimension(:,:), intent(in)::d_r
        real(dp), intent(in)::d_e_fix
        type(grid_i),intent(in)::domain!
        type(crop_pars_matrices),intent(inout)::crop_par_mat
        
        ! populate pheno%RF_t & pheno%RF_e
        crop_par_mat%RF_t = merge(crop_par_mat%max_RF_t * (d_r / (crop_par_mat%max_d_r-d_e_fix)), crop_par_mat%RF_t, crop_par_mat%max_RF_t /= domain%header%nan)
        crop_par_mat%RF_e = merge(1- crop_par_mat%RF_t, crop_par_mat%RF_e, crop_par_mat%max_RF_t /= domain%header%nan)    

    end subroutine calculate_RF_t!
    
end module mod_crop_phenology
