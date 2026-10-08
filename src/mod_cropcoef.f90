module mod_cropcoef

use mod_constants, only: dp, pi, nan_r
use mod_crop_phenology, only: crop_definition, crop_rotation, crop_pars_matrices
use mod_evapotranspiration, only: calculateDLH
use mod_grid, only: grid_i
use mod_meteo, only: meteo_mat
use mod_meteo, only: meteo_info
use mod_parameters, only: simulation
use mod_date, only: date, days_between_dates, advance_calendar_date

implicit none
private
public :: advance_crops_daily, adjust_water_productivity, compute_canopy_resistance
public :: crop_weather_cache, load_crop_weather_cache

! Station data are cached once so a newly sown crop can look ahead to the
! complete weather interval used by CropCoef to correct a Kcb plateau.
type crop_weather_cache
    type(date) :: first
    real(dp), allocatable :: tmax(:,:), tmin(:,:), rhmin(:,:), wind(:,:)
end type crop_weather_cache

contains

subroutine load_crop_weather_cache(cache, stations, sim)
    type(crop_weather_cache), intent(inout) :: cache
    type(meteo_info), intent(in) :: stations(:)
    type(simulation), intent(in) :: sim
    integer :: station, day, unit, ios, n_days, header_line
    real(dp) :: rain, rhmax, radiation
    character(len=300) :: header

    cache%first = stations(1)%start
    n_days = days_between_dates(stations(1)%start, stations(1)%finish) + 1
    allocate(cache%tmax(n_days,size(stations)), cache%tmin(n_days,size(stations)), &
             & cache%rhmin(n_days,size(stations)), cache%wind(n_days,size(stations)))

    do station = 1, size(stations)
        open(newunit=unit, file=trim(sim%meteo_path)//trim(stations(station)%filename), &
             & status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, 'Cannot reopen meteorological series for Kcb correction: ', trim(stations(station)%filename)
            stop
        end if
        do header_line = 1, 4
            read(unit, '(A)', iostat=ios) header
            if (ios /= 0) exit
        end do
        do day = 1, n_days
            if (ios /= 0) exit
            read(unit, *, iostat=ios) cache%tmax(day,station), cache%tmin(day,station), &
                & rain, rhmax, cache%rhmin(day,station), cache%wind(day,station), radiation
        end do
        close(unit)
        if (ios /= 0) then
            print *, 'Incomplete meteorological series for Kcb correction: ', trim(stations(station)%filename)
            stop
        end if
    end do
end subroutine load_crop_weather_cache

subroutine project_corrected_kcb(corrected, crop, cache, stations, weights, sim, sowing_date, latitude)
    ! CropCoef corrects declining control points and whole Kcb plateaus with
    ! mean weather from the interval between the corresponding GDD points.
    ! Project the crop once at sowing, then interpolate its corrected points
    ! during the daily water balance. The projection uses this cell's weather.
    real(dp), intent(out) :: corrected(:)
    type(crop_definition), intent(in) :: crop
    type(crop_weather_cache), intent(in) :: cache
    integer, intent(in) :: stations(:)
    real(dp), intent(in) :: weights(:)
    type(simulation), intent(in) :: sim
    type(date), intent(in) :: sowing_date
    real(dp), intent(in) :: latitude

    real(dp), allocatable :: rh(:), wind(:), height(:)
    integer, allocatable :: crossing(:)
    type(date) :: projected_date
    real(dp) :: tmax, tmin, gdd, vern_days, daylength, temperature, correction
    integer :: first_day, n_days, day, point, plateau_end, s, e

    corrected = crop%k_cb
    if (.not. allocated(cache%tmax)) return
    first_day = days_between_dates(cache%first, sowing_date) + 1
    if (first_day < 1 .or. first_day > size(cache%tmax,1)) return

    n_days = min(size(cache%tmax,1) - first_day + 1, 3*366)
    allocate(rh(n_days), wind(n_days), height(n_days), crossing(size(crop%gdd)))
    crossing = 0
    projected_date = sowing_date
    gdd = 0._dp
    vern_days = 0._dp

    do day = 1, n_days
        call cell_projection_weather(cache, first_day + day - 1, stations, weights, sim, &
                                     & tmax, tmin, rh(day), wind(day))
        temperature = (tmax + tmin)/2._dp
        call calculateDLH(projected_date%doy, latitude, daylength)
        gdd = gdd + calculate_daily_gdd(tmax, tmin, crop%base_temp, crop%cutoff_temp) * &
            & min(vernalization_factor(crop, temperature, vern_days), photoperiod_factor(crop, daylength))
        height(day) = curve_at_gdd(crop%gdd, crop%height, gdd)
        do point = 1, size(crossing)
            if (crossing(point) == 0 .and. gdd >= crop%gdd(point)) crossing(point) = day
        end do
        if (gdd >= maxval(crop%gdd)) exit
        call advance_calendar_date(projected_date)
    end do

    point = 1
    do while (point < size(corrected))
        if (crop%k_cb(point) == crop%k_cb(point+1) .and. crop%k_cb(point) > 0.45_dp) then
            plateau_end = point + 1
            do while (plateau_end < size(corrected))
                if (crop%k_cb(plateau_end+1) /= crop%k_cb(point)) exit
                plateau_end = plateau_end + 1
            end do
            s = crossing(point)
            e = crossing(plateau_end)
            if (s > 0 .and. e >= s) then
                correction = k_cb_weather_correction(sum(rh(s:e))/real(e-s+1,dp), &
                    & sum(wind(s:e))/real(e-s+1,dp), sum(height(s:e))/real(e-s+1,dp))
                corrected(point:plateau_end) = crop%k_cb(point:plateau_end) + correction
            end if
            point = plateau_end
        else
            point = point + 1
        end if
    end do

    do point = 2, size(corrected)
        if (crop%k_cb(point) >= crop%k_cb(point-1)) cycle
        if (crop%k_cb(point) <= 0.45_dp .or. crop%k_cb(point-1) <= 0.45_dp) cycle
        s = crossing(point-1)
        e = crossing(point)
        if (s == 0 .or. e < s) cycle
        correction = k_cb_weather_correction(sum(rh(s:e))/real(e-s+1,dp), &
            & sum(wind(s:e))/real(e-s+1,dp), sum(height(s:e))/real(e-s+1,dp))
        corrected(point) = crop%k_cb(point) + correction
    end do
end subroutine project_corrected_kcb

subroutine cell_projection_weather(cache, day, stations, weights, sim, tmax, tmin, rhmin, wind)
    type(crop_weather_cache), intent(in) :: cache
    integer, intent(in) :: day, stations(:)
    real(dp), intent(in) :: weights(:)
    type(simulation), intent(in) :: sim
    real(dp), intent(out) :: tmax, tmin, rhmin, wind
    integer :: k, station

    tmax = cache%tmax(day,stations(1))
    tmin = cache%tmin(day,stations(1))
    rhmin = cache%rhmin(day,stations(1))
    wind = cache%wind(day,stations(1))
    if (sim%interpolate_temp) then
        tmax = 0._dp
        tmin = 0._dp
        do k = 1, size(stations)
            station = stations(k)
            tmax = tmax + cache%tmax(day,station)*weights(k)
            tmin = tmin + cache%tmin(day,station)*weights(k)
        end do
    end if
    if (sim%interpolate_hum) then
        rhmin = 0._dp
        do k = 1, size(stations)
            rhmin = rhmin + cache%rhmin(day,stations(k))*weights(k)
        end do
    end if
    if (sim%interpolate_wind) then
        wind = 0._dp
        do k = 1, size(stations)
            wind = wind + cache%wind(day,stations(k))*weights(k)
        end do
    end if
end subroutine cell_projection_weather

real(dp) function k_cb_weather_correction(rhmin, wind, height) result(correction)
    ! Same clamped FAO-56 correction used by CropCoef's calcKcbCorrFact().
    real(dp), intent(in) :: rhmin, wind, height

    correction = (0.04_dp*(min(6._dp,max(1._dp,wind)) - 2._dp) - &
        & 0.004_dp*(min(80._dp,max(20._dp,rhmin)) - 45._dp)) * &
        & (min(10._dp,max(0.1_dp,height))/3._dp)**0.3_dp
end function k_cb_weather_correction

! CropCoef's CO2 correction of normalized biomass water productivity.
real(dp) function adjust_water_productivity(orig_wp, co2, sink_strength) result(adjusted_wp)
    real(dp), intent(in) :: orig_wp, co2, sink_strength
    real(dp) :: crop_type, weight, co2_factor

    crop_type = min(1._dp, max(0._dp, (40._dp - orig_wp*100._dp)/20._dp))
    weight = min(1._dp, max(0._dp, 1._dp - (550._dp - co2)/(550._dp - 369.41_dp)))
    co2_factor = (co2/369.41_dp) / (1._dp + (co2 - 369.41_dp) * &
        & ((1._dp - weight)*0.000138_dp + weight*(sink_strength*0.000138_dp + &
        & (1._dp - sink_strength)*0.001165_dp)))
    adjusted_wp = (1._dp + crop_type*(co2_factor - 1._dp))*orig_wp
end function adjust_water_productivity

! CropCoef's annual canopy-resistance correction (EPIC/FAO56).
real(dp) function compute_canopy_resistance(co2) result(resistance)
    real(dp), intent(in) :: co2
    real(dp) :: leaf_conductance

    leaf_conductance = 0.01_dp * (1.4_dp - 0.4_dp*(co2/330._dp))
    resistance = (1._dp/leaf_conductance)/(0.5_dp*24._dp*0.12_dp)
end function compute_canopy_resistance

subroutine advance_crops_daily(pheno, definitions, rotations, domain, soil_use, irandom, meteo, &
                               & cache, station_idx, station_weight, sim, current)
    type(crop_pars_matrices), intent(inout) :: pheno
    type(crop_definition), intent(in) :: definitions(:)
    type(crop_rotation), intent(inout) :: rotations(:)
    type(grid_i), intent(in) :: domain, soil_use
    integer, intent(in) :: irandom(:,:)
    type(meteo_mat), intent(in) :: meteo
    type(crop_weather_cache), intent(in) :: cache
    integer, intent(in) :: station_idx(:,:,:)
    real(dp), intent(in) :: station_weight(:,:,:)
    type(simulation), intent(in) :: sim
    type(date), intent(in) :: current

    integer :: i, j, land_use, position, crop_id
    real(dp) :: daily_gdd, daylength, temperature
    integer :: year, doy

    year = current%year
    doy = current%doy

    pheno%k_cb_old = pheno%k_cb

    do j = 1, size(domain%mat, 2)
        do i = 1, size(domain%mat, 1)
            if (domain%mat(i,j) == domain%header%nan) cycle ! Skip inactive cells

            land_use = soil_use%mat(i,j)
            if (land_use < 1 .or. land_use > size(rotations)) then ! Todo: move this check upstream?
                call set_bare_soil_today(pheno, i, j)
                cycle
            end if

            ! Skip landuses that don't have any crop
            if (.not. allocated(rotations(land_use)%crop_ids)) then
                call set_bare_soil_today(pheno, i, j)
                cycle
            end if

            ! Skip landuses that only include crop "0" (baresoil)
            if (size(rotations(land_use)%crop_ids) == 0) then
                call set_bare_soil_today(pheno, i, j)
                cycle
            end if

            position = pheno%rotation_position(i,j)
            ! A new annual land-use map can replace a two-crop rotation with
            ! a shorter one. Its old position is then outside the new rotation.
            if (position > size(rotations(land_use)%crop_ids)) then
                position = 1
                pheno%rotation_position(i,j) = 1
                pheno%harvest_pending(i,j) = .false.
                pheno%crop_id(i,j) = 0
                pheno%gdd(i,j) = 0._dp
                pheno%vernalization_days(i,j) = 0._dp
                pheno%bare_soil_days_left(i,j) = 0
            end if

            ! Harvest the crop (becoming baresoil) if due
            if (pheno%harvest_pending(i,j)) then
                pheno%harvest_pending(i,j) = .false.
                pheno%crop_id(i,j) = 0
                pheno%gdd(i,j) = 0._dp
                pheno%vernalization_days(i,j) = 0._dp
                pheno%cuts_completed(i,j) = 0
                pheno%pheno_idx(i,j) = 0

                ! Move to next crop in the rotation
                position = merge(1, position + 1, position == size(rotations(land_use)%crop_ids))
                pheno%rotation_position(i,j) = position
                crop_id = rotations(land_use)%crop_ids(position)
                if (crop_id >= 1 .and. crop_id <= size(definitions)) then
                    pheno%bare_soil_days_left(i,j) = max(0, definitions(crop_id)%crop_overlap_days)
                end if
            end if

            ! Soil is bare: wait crop_overlap_days before considering sowing
            if (pheno%crop_id(i,j) == 0) then
                call set_bare_soil_today(pheno, i, j)
                if (pheno%bare_soil_days_left(i,j) > 0) then
                    pheno%bare_soil_days_left(i,j) = pheno%bare_soil_days_left(i,j) - 1
                    cycle
                end if

                call skip_missed_sowing_windows(position, rotations(land_use), definitions, doy, irandom(i,j), land_use, year)
                pheno%rotation_position(i,j) = position
            end if

            crop_id = rotations(land_use)%crop_ids(position)
            if (crop_id < 1 .or. crop_id > size(definitions)) then
                call set_bare_soil_today(pheno, i, j)
                cycle
            end if
            if (.not. allocated(definitions(crop_id)%gdd)) then
                call set_bare_soil_today(pheno, i, j)
                cycle
            end if
            if (size(definitions(crop_id)%gdd) == 0) then
                call set_bare_soil_today(pheno, i, j)
                cycle
            end if

            temperature = (meteo%T_max(i,j) + meteo%T_min(i,j))/2._dp

            if (pheno%crop_id(i,j) == 0) then
                if (.not. should_sow_today(definitions(crop_id), doy, irandom(i,j), temperature)) cycle

                ! Sow a new crop: initialize phenology and GDD accumulation
                pheno%crop_id(i,j) = crop_id
                pheno%sowing_year(i,j) = year
                pheno%sowing_doy(i,j) = doy
                pheno%gdd(i,j) = 0._dp
                pheno%vernalization_days(i,j) = 0._dp
                pheno%cuts_completed(i,j) = 0
                pheno%pheno_idx(i,j) = 0
                if (definitions(crop_id)%adjust_k_cb) then
                    call project_corrected_kcb(pheno%corrected_k_cb(i,j,1:size(definitions(crop_id)%k_cb)), &
                        & definitions(crop_id), cache, station_idx(i,j,:), station_weight(i,j,:), &
                        & sim, current, meteo%lat(i,j))
                else
                    pheno%corrected_k_cb(i,j,1:size(definitions(crop_id)%k_cb)) = definitions(crop_id)%k_cb
                end if
            end if

            ! Crop is in field: calculate GDD accumulation
            call calculateDLH(doy, meteo%lat(i,j), daylength)
            daily_gdd = calculate_daily_gdd(meteo%T_max(i,j), meteo%T_min(i,j), &
                                           & definitions(crop_id)%base_temp, definitions(crop_id)%cutoff_temp)
            daily_gdd = daily_gdd * min(vernalization_factor(definitions(crop_id), temperature, pheno%vernalization_days(i,j)), &
                                      & photoperiod_factor(definitions(crop_id), daylength))
            pheno%gdd(i,j) = pheno%gdd(i,j) + daily_gdd

            ! Use the new GDD value to estimate crop state
            call populate_crop_today(pheno, definitions(crop_id), i, j)

            ! If harvest GDD has been reached, mark the crop; it will be harvested tomorrow
            if (pheno%gdd(i,j) >= real(max(1, definitions(crop_id)%max_cuts), dp) * maxval(definitions(crop_id)%gdd)) then
                pheno%harvest_pending(i,j) = .true.
            end if
        end do
    end do
end subroutine advance_crops_daily

subroutine skip_missed_sowing_windows(position, rotation, definitions, doy, random_shift, land_use, year)
    integer, intent(inout) :: position
    type(crop_rotation), intent(inout) :: rotation
    type(crop_definition), intent(in) :: definitions(:)
    integer, intent(in) :: doy, random_shift, land_use, year
    integer :: considered, crop_id, first_day, last_day

    do considered = 1, size(rotation%crop_ids)
        crop_id = rotation%crop_ids(position)
        if (crop_id < 1 .or. crop_id > size(definitions)) return

        call sowing_window(definitions(crop_id), random_shift, first_day, last_day)
        if (doy <= last_day) return

        if (size(rotation%crop_ids) > 1 .and. .not. rotation%missed_sowing_warned) then
            print *, 'Warning: land use ', land_use, ' missed the sowing window for ', &
                     & trim(definitions(crop_id)%crop_name), ' in ', year, ' (day ', doy, ').'
            print *, 'The rotation advances; further missed-sowing warnings for this land use are suppressed.'
            rotation%missed_sowing_warned = .true.
        end if

        position = merge(1, position + 1, position == size(rotation%crop_ids))
    end do
end subroutine skip_missed_sowing_windows

subroutine sowing_window(crop, random_shift, first_day, last_day)
    type(crop_definition), intent(in) :: crop
    integer, intent(in) :: random_shift
    integer, intent(out) :: first_day, last_day

    first_day = max(1, crop%sowing_doy_min + random_shift)
    last_day = first_day + max(0, crop%sowing_delay_max)
end subroutine sowing_window

logical function should_sow_today(crop, doy, random_shift, temperature)
    type(crop_definition), intent(in) :: crop
    integer, intent(in) :: doy, random_shift
    real(dp), intent(in) :: temperature
    integer :: first_day, last_day

    call sowing_window(crop, random_shift, first_day, last_day)

    should_sow_today = doy >= first_day .and. doy <= last_day .and. &
                       & (temperature > crop%sowing_temp .or. doy == last_day)
end function should_sow_today

real(dp) function calculate_daily_gdd(tmax, tmin, base, cutoff) result(gdd)
    ! Single-day form of CropCoef's sine-wave calculateGDD() (Borghi, eq. i-85).
    real(dp), intent(in) :: tmax, tmin, base, cutoff
    real(dp) :: average, amplitude, lower_area, upper_area, theta, phi

    gdd = 0._dp
    if (tmax <= base .or. cutoff <= base) return

    average = (tmax + tmin)/2._dp
    amplitude = (tmax - tmin)/2._dp
    lower_area = 0._dp
    upper_area = 0._dp

    if (tmin >= cutoff) then
        gdd = cutoff - base
        return
    else if (tmin >= base .and. tmax <= cutoff) then
        lower_area = average - base
    else if (tmin < base .and. tmax <= cutoff) then
        theta = asin(max(-1._dp, min(1._dp, (base - average)/amplitude)))
        lower_area = ((average - base)*(pi/2._dp - theta) + amplitude*cos(theta))/pi
    else if (tmin >= base .and. tmax > cutoff) then
        phi = asin(max(-1._dp, min(1._dp, (cutoff - average)/amplitude)))
        lower_area = average - base
        upper_area = ((average - cutoff)*(pi/2._dp - phi) + amplitude*cos(phi))/pi
    else
        theta = asin(max(-1._dp, min(1._dp, (base - average)/amplitude)))
        phi = asin(max(-1._dp, min(1._dp, (cutoff - average)/amplitude)))
        lower_area = ((average - base)*(pi/2._dp - theta) + amplitude*cos(theta))/pi
        upper_area = ((average - cutoff)*(pi/2._dp - phi) + amplitude*cos(phi))/pi
    end if

    gdd = max(0._dp, lower_area - upper_area)
end function calculate_daily_gdd

real(dp) function vernalization_factor(crop, temperature, accumulated_days) result(factor)
    ! CropCoef's vernalization contribution and factor, evaluated incrementally.
    type(crop_definition), intent(in) :: crop
    real(dp), intent(in) :: temperature
    real(dp), intent(inout) :: accumulated_days
    real(dp) :: contribution

    factor = 1._dp
    if (.not. crop%requires_vernalization) return

    contribution = 0._dp
    if (crop%vern_curve_slope > 0._dp) then
        if (temperature >= crop%vern_temp_min - crop%vern_curve_slope .and. &
            temperature < crop%vern_temp_min) then
            contribution = 1._dp - (crop%vern_temp_min - temperature)/crop%vern_curve_slope
        else if (temperature >= crop%vern_temp_max .and. &
                 temperature < crop%vern_temp_max + crop%vern_curve_slope) then
            contribution = 1._dp - (temperature - crop%vern_temp_max)/crop%vern_curve_slope
        end if
    end if
    if (temperature >= crop%vern_temp_min .and. temperature < crop%vern_temp_max) contribution = 1._dp
    accumulated_days = accumulated_days + contribution

    if (accumulated_days < crop%vern_days_start .or. accumulated_days > crop%vern_days_end) return
    if (crop%vern_days_end <= crop%vern_days_start) return

    factor = crop%vern_fact_min + (1._dp - crop%vern_fact_min) * &
             & (accumulated_days - crop%vern_days_start) / (crop%vern_days_end - crop%vern_days_start)
end function vernalization_factor

real(dp) function photoperiod_factor(crop, daylength) result(factor)
    ! Long-day and short-day responses from CropCoef's photoperiod().
    type(crop_definition), intent(in) :: crop
    real(dp), intent(in) :: daylength

    factor = 1._dp
    select case (crop%photoperiod_response)
    case (1)
        if (daylength < crop%daylength_if) then
            factor = 0._dp
        else if (daylength <= crop%daylength_ins .and. crop%daylength_ins > crop%daylength_if) then
            factor = (daylength - crop%daylength_if)/(crop%daylength_ins - crop%daylength_if)
        end if
    case (2)
        if (daylength > crop%daylength_if) then
            factor = 0._dp
        else if (daylength >= crop%daylength_ins .and. crop%daylength_if > crop%daylength_ins) then
            factor = (crop%daylength_if - daylength)/(crop%daylength_if - crop%daylength_ins)
        end if
    end select
end function photoperiod_factor

real(dp) function curve_at_gdd(gdd_points, values, accumulated_gdd) result(value)
    ! Interpolate today's parameter directly from accumulated thermal time.
    real(dp), intent(in) :: gdd_points(:), values(:), accumulated_gdd
    integer :: point

    value = nan_r
    if (size(gdd_points) == 0) return
    if (accumulated_gdd <= gdd_points(1)) then
        value = values(1)
        return
    end if
    do point = 2, size(gdd_points)
        if (accumulated_gdd <= gdd_points(point)) then
            if (gdd_points(point) > gdd_points(point-1)) then
                value = values(point-1) + (values(point) - values(point-1)) * &
                      & (accumulated_gdd - gdd_points(point-1))/(gdd_points(point) - gdd_points(point-1))
            else
                value = values(point)
            end if
            return
        end if
    end do
    value = values(size(values))
end function curve_at_gdd

subroutine populate_crop_today(pheno, crop, i, j)
    type(crop_pars_matrices), intent(inout) :: pheno
    type(crop_definition), intent(in) :: crop
    integer, intent(in) :: i, j
    integer :: point
    real(dp) :: current_gdd

    current_gdd = pheno%gdd(i,j)
    ! Multi-cut crops also need a reset of the within-cut curve at each cut.
    pheno%k_cb(i,j) = curve_at_gdd(crop%gdd, pheno%corrected_k_cb(i,j,1:size(crop%k_cb)), current_gdd)
    pheno%lai(i,j) = curve_at_gdd(crop%gdd, crop%lai, current_gdd)
    pheno%h(i,j) = curve_at_gdd(crop%gdd, crop%height, current_gdd)
    pheno%d_r(i,j) = curve_at_gdd(crop%gdd, crop%root_depth, current_gdd)
    pheno%f_c(i,j) = curve_at_gdd(crop%gdd, crop%cover_fraction, current_gdd)
    pheno%r_stress(i,j) = curve_at_gdd(crop%gdd, crop%stress_resistance, current_gdd)
    pheno%cn_day(i,j) = nint(curve_at_gdd(crop%gdd, crop%cn_value, current_gdd))
    pheno%cn_class(i,j) = crop%cn_class
    pheno%irrigation_class(i,j) = merge(1, 0, crop%is_irrigated)
    pheno%p(i,j) = crop%raw_fraction
    pheno%a(i,j) = crop%interception_coef
    pheno%d_t_max(i,j) = maxval(crop%root_depth)
    pheno%RF_t_max(i,j) = crop%maximum_transpirative_root_fraction
    pheno%T_crit(i,j) = crop%heat_stress_temp_crit
    pheno%T_lim(i,j) = crop%heat_stress_temp_lim
    pheno%k_cb_low(i,j) = minval(pheno%corrected_k_cb(i,j,1:size(crop%k_cb)))
    pheno%k_cb_high(i,j) = maxval(pheno%corrected_k_cb(i,j,1:size(crop%k_cb)))
    pheno%k_cb_mid(i,j) = (pheno%k_cb_low(i,j) + pheno%k_cb_high(i,j))/2._dp
    do point = 2, size(crop%k_cb)
        if (crop%k_cb(point) == crop%k_cb(point-1) .and. &
            & crop%k_cb(point) > minval(crop%k_cb) .and. crop%k_cb(point) < maxval(crop%k_cb)) then
            pheno%k_cb_mid(i,j) = pheno%corrected_k_cb(i,j,point)
            exit
        end if
    end do
    pheno%pheno_idx(i,j) = max(1, pheno%pheno_idx(i,j))
end subroutine populate_crop_today

subroutine set_bare_soil_today(pheno, i, j)
    type(crop_pars_matrices), intent(inout) :: pheno
    integer, intent(in) :: i, j

    pheno%crop_id(i,j) = 0
    pheno%pheno_idx(i,j) = 0
    pheno%k_cb(i,j) = 0._dp
    pheno%h(i,j) = 0._dp
    pheno%d_r(i,j) = 0._dp
    pheno%lai(i,j) = 0._dp
    pheno%cn_day(i,j) = 0
    pheno%f_c(i,j) = 0._dp
    pheno%r_stress(i,j) = 0._dp
    pheno%irrigation_class(i,j) = 0
    pheno%cn_class(i,j) = 1
    pheno%p(i,j) = 0._dp
    pheno%a(i,j) = 0._dp
    pheno%d_t_max(i,j) = 0._dp
    pheno%RF_t_max(i,j) = 0._dp
    pheno%T_lim(i,j) = 0._dp
    pheno%T_crit(i,j) = 0._dp
    pheno%k_cb_low(i,j) = 0._dp
    pheno%k_cb_mid(i,j) = 0._dp
    pheno%k_cb_high(i,j) = 0._dp
end subroutine set_bare_soil_today

end module mod_cropcoef
