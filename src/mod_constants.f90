module mod_constants
    implicit none
    
    ! Symbolic names for kind types of single- and double-precision reals
    integer, parameter :: sp = kind(1.0)
    integer, parameter :: dp = kind(1.0D0)

    ! set default NAN values
    integer, parameter :: nan_i = -9999
    real(dp), parameter :: nan_r = -9999.0
        
    real(dp), parameter :: cost_fwEva = 1 ! default fwEvaporation
    real(dp),dimension(24),parameter :: cost_f_eff_rain = &
        &   (/1.0/24.0,&
        &     1.0/24.0,&
        &     1.0/24.0,&
        &     1.0/24.0,&
        &     1.0/24.0,&
        &     1.0/24.0,&
        &     1.0/24.0,&
        &     1.0/24.0,&
        &     1.0/24.0,&
        &     1.0/24.0,&
        &     1.0/24.0,&
        &     1.0/24.0,&
        &     1.0/24.0,&
        &     1.0/24.0,&
        &     1.0/24.0,&
        &     1.0/24.0,&
        &     1.0/24.0,&
        &     1.0/24.0,&
        &     1.0/24.0,&
        &     1.0/24.0,&
        &     1.0/24.0,&
        &     1.0/24.0,&
        &     1.0/24.0,&
        &     1.0/24.0/)
    real(dp),dimension(24),parameter :: cost_fet0 =  &
        &   (/0.001784677,&
        &     0.001640849,&
        &     0.002308392,&
        &     0.002801005,&
        &     0.002505317,&
        &     0.001244509,&
        &     0.00546749,&
        &     0.018271901,&
        &     0.04046744,&
        &     0.074224044,&
        &     0.108566773,&
        &     0.140097667,&
        &     0.166536044,&
        &     0.166524601,&
        &     0.122840875,&
        &     0.072098911,&
        &     0.038270891,&
        &     0.015545331,&
        &     0.006179273,&
        &     0.004320957,&
        &     0.002402308,&
        &     0.002143947,&
        &     0.001850866,&
        &     0.001905931/)

    integer, dimension(9,3,4), parameter:: tabCN = reshape( & ! CN2 table
        & (/74, 65, 61, 58, 30, 32, 99, 77, 25, &    ! # (i,1,1)
        &    0, 72, 65, 66,  0, 43,  0,  0, 36, &    ! # (i,2,1)
        &    0,  0,  0,  0,  0, 57,  0,  0, 45, &    ! # (i,3,1)
        &   83, 75, 73, 72, 58, 58, 99, 86, 55, &    ! # (i,1,2)
        &    0, 81, 76, 77,  0, 65,  0,  0, 60, &    ! # (i,2,2)
        &    0,  0,  0,  0,  0, 73,  0,  0, 66, &    ! # (i,3,2)
        &   88, 82, 81, 81, 71, 72, 99, 91, 70, &    ! # (i,1,3)
        &    0, 88, 84, 85,  0, 76,  0,  0, 73, &    ! # (i,2,3)
        &    0,  0,  0,  0,  0, 82,  0,  0, 77, &    ! # (i,3,3)
        &   90, 86, 84, 85, 78, 79, 99, 94, 77, &    ! # (i,1,4)
        &    0, 91, 88, 89,  0, 82,  0,  0, 79, &    ! # (i,2,4)
        &    0,  0,  0,  0,  0, 86,  0,  0, 83/), &  ! # (i,3,4)
        & (/9,3,4/))                                 ! final form
        ! columns (1° number):  Cover type:        
                ! 1 - Crop Residue Cover; 2 - Row crops; 3 - Small grain
                ! 4 - Close-seeded or broadcast legumes or rotation meadow; 5 - Meadow; 6 - Woods – grass combination (orchard or tree farm)
                ! 7 - Rice; 8 - Fallow - Bare soil; 9 - Woods
        ! rows (2° number):     Hydrological condition
        ! levels (3° number):   Hydrological group
        
    integer, dimension(12), parameter:: tmax_time = &
        & (/14, 14, 14, 15, 15, 16, 15, 15, 15, 14, 14, 14/)    ! Tmax time for the months from Jan to Dec

    integer, dimension(12), parameter:: tmin_time = &
        & (/6, 6, 5, 5, 4, 4, 4, 5, 5, 6, 6, 7/)                ! Tmin time for the months from Jan to Dec
        
    real(dp), parameter :: pi= acos(-1.)
    
    real(dp), parameter :: res_canopy_std = 70.                 ! Standard canopy resistance (alfalfa) for ET0 calculation
end module mod_constants
