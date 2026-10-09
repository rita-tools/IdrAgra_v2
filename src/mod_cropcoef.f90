module mod_cropcoef

use mod_constants, only: dp, pi, nan_r
use mod_crop_phenology, only: crop_definition, crop_rotation, crop_pars_matrices
use mod_evapotranspiration, only: calculateDLH
use mod_grid, only: grid_i
use mod_meteo, only: meteo_mat
use mod_meteo, only: meteo_info
use mod_parameters, only: simulation
use mod_date, only: date, days_between_dates, days_in_year, month_lengths, calendar_date_from_doy

implicit none
private
public :: advance_crops_daily, adjust_water_productivity, compute_canopy_resistance, prepare_crop_stage_boundaries, &
        & compute_kcb_corrections

contains

subroutine compute_kcb_corrections(kcb_corr_fact, crops, stations, sim)
    real(dp), allocatable :: kcb_corr_fact(:,:)
    type(crop_definition), intent(in) :: crops(:)
    type(meteo_info), intent(in) :: stations(:)
    type(simulation), intent(in) :: sim

    real(dp), allocatable :: rhmin(:), wind(:)
    real(dp) :: tmax, tmin, rain, rhmax, radiation, mean_rhmin, mean_wind
    integer :: crop_idx, station_idx, n_days, day_idx, n_years, partial_days, unit, ios, header_line
    character(len=300) :: header

    allocate(kcb_corr_fact(size(crops),size(stations)), source=0._dp)

    do crop_idx = 1, size(crops)
        if (.not. needs_kcb_adjustment(crops(crop_idx))) cycle
        if ((crops(crop_idx)%kcb_corr_start_month == 0) .neqv. (crops(crop_idx)%kcb_corr_end_month == 0)) then
            print *, 'Both KcbCorrectionStart and KcbCorrectionEnd are required in ', &
                & trim(crops(crop_idx)%parameter_file), '.'
            print *, 'Execution will be aborted...'
            stop
        end if
        if (crops(crop_idx)%kcb_corr_start_month == 0) then
            print *, 'Warning: ', trim(crops(crop_idx)%parameter_file), ' has no Kcb correction interval.'
            print *, 'Using 70-95% of the estimated sowing-to-harvest period for climate averages.'
        end if
    end do

    do station_idx = 1, size(stations)
        n_days = days_between_dates(stations(station_idx)%start, stations(station_idx)%finish) + 1
        allocate(rhmin(n_days), wind(n_days))
        open(newunit=unit, file=trim(sim%meteo_path)//trim(stations(station_idx)%filename), &
             & status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, 'Cannot reopen meteorological series for Kcb correction: ', trim(stations(station_idx)%filename)
            print *, 'Execution will be aborted...'
            stop
        end if
        do header_line = 1, 4
            read(unit, '(A)', iostat=ios) header
            if (ios /= 0) exit
        end do
        do day_idx = 1, n_days
            if (ios /= 0) exit
            read(unit, *, iostat=ios) tmax, tmin, rain, rhmax, rhmin(day_idx), wind(day_idx), radiation ! Only store rhmin and wind
        end do
        close(unit)
        if (ios /= 0) then
            print *, 'Incomplete meteorological series for Kcb correction: ', trim(stations(station_idx)%filename)
            print *, 'Execution will be aborted...'
            stop
        end if

        do crop_idx = 1, size(crops)
            if (.not. needs_kcb_adjustment(crops(crop_idx))) cycle

            ! Get average RHmin and wind speed over the crop's Kcb correction interval
            call mean_correction_weather(crops(crop_idx), stations(station_idx), rhmin, wind, &
                                       & mean_rhmin, mean_wind, n_years, partial_days         )
            if (n_years == 0) then
                if (partial_days > 0) then
                    print *, 'Warning: no complete Kcb correction interval for ', trim(crops(crop_idx)%crop_name), &
                        & ' at weather station ', trim(stations(station_idx)%filename), '.'
                    print *, 'Used ', partial_days, ' available day(s) within incomplete interval(s).'
                else
                    print *, 'Warning: no Kcb correction interval overlaps weather station ', &
                        & trim(stations(station_idx)%filename), ' for ', trim(crops(crop_idx)%crop_name), '.'
                    print *, 'Used the entire available weather series instead.'
                end if
            end if

            ! Calculate the correction factor
            !%PS%, note: FAO56 prescribe using the mean plant height during mid-late season; here we use max height for simplicity
            kcb_corr_fact(crop_idx,station_idx) = k_cb_weather_correction(mean_rhmin, mean_wind, maxval(crops(crop_idx)%height))
        end do
        deallocate(rhmin, wind)
    end do
end subroutine compute_kcb_corrections

function needs_kcb_adjustment(crop) result(is_eligible)
    type(crop_definition), intent(in) :: crop
    logical :: is_eligible

    is_eligible = .false.
    if (.not. crop%adjust_k_cb) return
    if (.not. allocated(crop%k_cb)) return
    if (size(crop%k_cb) == 0) return
    is_eligible = crop%mid_start_gdd < huge(0._dp) .and. crop%mid_end_gdd > crop%mid_start_gdd .and. maxval(crop%k_cb) > 0.45_dp
end function needs_kcb_adjustment

! Calculate mean RHmin and wind speed over the crop's Kcb correction interval
subroutine mean_correction_weather(crop, station, rhmin, wind, mean_rhmin, mean_wind, n_years, partial_days)
    type(crop_definition), intent(in) :: crop
    type(meteo_info), intent(in) :: station
    real(dp), intent(in) :: rhmin(:), wind(:)
    real(dp), intent(out) :: mean_rhmin, mean_wind
    integer, intent(out) :: n_years, partial_days

    type(date) :: first_day, last_day
    real(dp) :: partial_rhmin, partial_wind
    integer :: year, first_idx, last_idx, n_days, start_key, end_key, days_in_month(12)
    integer :: sowing_idx, harvest_idx, harvest_year, estimated_sowing_doy, growing_days

    mean_rhmin = 0._dp
    mean_wind = 0._dp
    n_years = 0
    partial_rhmin = 0._dp
    partial_wind = 0._dp
    partial_days = 0
    start_key = 100*crop%kcb_corr_start_month + crop%kcb_corr_start_day
    end_key = 100*crop%kcb_corr_end_month + crop%kcb_corr_end_day

    do year = station%start%year - 1, station%finish%year
        if (crop%kcb_corr_start_month == 0) then
            ! Crop file does not specify a Kcb correction interval; use 70-95% of the estimated sowing-to-harvest period
            estimated_sowing_doy = max(1, crop%sowing_doy_min + max(0, crop%sowing_delay_max)/2)
            sowing_idx = days_between_dates(station%start, calendar_date_from_doy(year, 1)) + estimated_sowing_doy
            harvest_year = year + merge(1, 0, crop%harvest_doy_max <= estimated_sowing_doy)
            harvest_idx = days_between_dates(station%start, calendar_date_from_doy(harvest_year, 1)) + &
                & min(days_in_year(harvest_year), max(1, crop%harvest_doy_max))
            growing_days = harvest_idx - sowing_idx
            if (growing_days <= 0) cycle
            first_idx = sowing_idx + nint(0.70_dp*real(growing_days,dp))
            last_idx = sowing_idx + nint(0.95_dp*real(growing_days,dp))
        else
            first_day%day = crop%kcb_corr_start_day
            first_day%month = crop%kcb_corr_start_month
            first_day%year = year
            last_day%day = crop%kcb_corr_end_day
            last_day%month = crop%kcb_corr_end_month
            last_day%year = year + merge(1, 0, start_key > end_key)
            days_in_month = month_lengths(first_day%year)
            if (first_day%day > days_in_month(first_day%month)) cycle
            days_in_month = month_lengths(last_day%year)
            if (last_day%day > days_in_month(last_day%month)) cycle
            first_idx = days_between_dates(station%start, first_day) + 1
            last_idx = days_between_dates(station%start, last_day) + 1
        end if
        if (last_idx < first_idx) cycle
        if (first_idx < 1 .or. last_idx > size(rhmin)) then
            first_idx = max(1, first_idx)
            last_idx = min(size(rhmin), last_idx)
            if (last_idx >= first_idx) then
                partial_rhmin = partial_rhmin + sum(rhmin(first_idx:last_idx))
                partial_wind = partial_wind + sum(wind(first_idx:last_idx))
                partial_days = partial_days + last_idx - first_idx + 1
            end if
            cycle
        end if
        n_days = last_idx - first_idx + 1
        mean_rhmin = mean_rhmin + sum(rhmin(first_idx:last_idx))/real(n_days,dp)
        mean_wind = mean_wind + sum(wind(first_idx:last_idx))/real(n_days,dp)
        n_years = n_years + 1
    end do
    if (n_years > 0) then
        mean_rhmin = mean_rhmin/real(n_years,dp)
        mean_wind = mean_wind/real(n_years,dp)
    else if (partial_days > 0) then
        mean_rhmin = partial_rhmin/real(partial_days,dp)
        mean_wind = partial_wind/real(partial_days,dp)
    else
        mean_rhmin = sum(rhmin)/real(size(rhmin),dp)
        mean_wind = sum(wind)/real(size(wind),dp)
    end if
end subroutine mean_correction_weather

! Same clamped FAO-56 correction used by CropCoef's calcKcbCorrFact().
real(dp) function k_cb_weather_correction(rhmin, wind, height) result(correction)
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
                               & kcb_corr_fact, station_idx, current)
    type(crop_pars_matrices), intent(inout) :: pheno
    type(crop_definition), intent(in) :: definitions(:)
    type(crop_rotation), intent(inout) :: rotations(:)
    type(grid_i), intent(in) :: domain, soil_use
    integer, intent(in) :: irandom(:,:)
    type(meteo_mat), intent(in) :: meteo
    real(dp), intent(in) :: kcb_corr_fact(:,:)
    integer, intent(in) :: station_idx(:,:,:)
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
                pheno%pheno_stage(i,j) = 0

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
                pheno%pheno_stage(i,j) = 0
                pheno%is_real_crop(i,j) = maxval(definitions(crop_id)%gdd) > 0._dp
            end if

            ! Crop is in field: calculate GDD accumulation
            call calculateDLH(doy, meteo%lat(i,j), daylength)
            daily_gdd = calculate_daily_gdd(meteo%T_max(i,j), meteo%T_min(i,j), &
                                           & definitions(crop_id)%base_temp, definitions(crop_id)%cutoff_temp)
            daily_gdd = daily_gdd * min(vernalization_factor(definitions(crop_id), temperature, pheno%vernalization_days(i,j)), &
                                      & photoperiod_factor(definitions(crop_id), daylength))
            pheno%gdd(i,j) = pheno%gdd(i,j) + daily_gdd

            ! Use the new GDD value to estimate crop state
            call populate_crop_today(pheno, definitions(crop_id), kcb_corr_fact(crop_id,station_idx(i,j,1)), i, j)

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

subroutine populate_crop_today(pheno, crop, kcb_correction, i, j)
    type(crop_pars_matrices), intent(inout) :: pheno
    type(crop_definition), intent(in) :: crop
    real(dp), intent(in) :: kcb_correction
    integer, intent(in) :: i, j
    real(dp) :: current_gdd

    current_gdd = pheno%gdd(i,j)
    ! Multi-cut crops also need a reset of the within-cut curve at each cut.
    pheno%k_cb(i,j) = curve_at_gdd(crop%gdd, crop%k_cb, current_gdd)
    if (current_gdd >= crop%mid_start_gdd .and. pheno%k_cb(i,j) > 0.45_dp) then
        pheno%k_cb(i,j) = pheno%k_cb(i,j) + kcb_correction
    end if
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
    pheno%pheno_stage(i,j) = crop_growth_stage(crop, current_gdd)
end subroutine populate_crop_today

! Calculate GDD thresholds for crop growth stages once per crop and store them for later use
subroutine prepare_crop_stage_boundaries(crop)
    type(crop_definition), intent(inout) :: crop

    integer :: initial_start, initial_end, mid_start, mid_end, point
    real(dp) :: k_cb_mid

    crop%emergence_gdd = huge(0._dp)
    crop%initial_end_gdd = huge(0._dp)
    crop%mid_start_gdd = huge(0._dp)
    crop%mid_end_gdd = huge(0._dp)

    initial_start = 0
    do point = 1, size(crop%k_cb)
        if (crop%k_cb(point) > 0._dp) then
            initial_start = point
            exit
        end if
    end do
    if (initial_start == 0) return

    crop%emergence_gdd = crop%gdd(initial_start)
    initial_end = initial_start
    do point = initial_start + 1, size(crop%k_cb)
        if (crop%k_cb(point) /= crop%k_cb(initial_start)) exit
        initial_end = point
    end do
    if (initial_end == size(crop%k_cb)) return ! Constant positive Kcb stays in stage 1.
    crop%initial_end_gdd = crop%gdd(initial_end)

    k_cb_mid = maxval(crop%k_cb)
    mid_start = 0
    do point = initial_end + 1, size(crop%k_cb)
        if (crop%k_cb(point) == k_cb_mid) then
            mid_start = point
            exit
        end if
    end do
    if (mid_start == 0) return ! The peak lies within the initial segment.

    crop%mid_start_gdd = crop%gdd(mid_start)
    mid_end = mid_start
    do point = mid_start + 1, size(crop%k_cb)
        if (crop%k_cb(point) /= k_cb_mid) exit
        mid_end = point
    end do
    crop%mid_end_gdd = crop%gdd(mid_end)
end subroutine prepare_crop_stage_boundaries

! 0 = pre-emergence, 1 = initial, 2 = development, 3 = mid-season, 4 = late-season.
! Stages 1 and 3 usually correspond to Kcb plateaus - see FAO56
integer function crop_growth_stage(crop, accumulated_gdd) result(stage)
    type(crop_definition), intent(in) :: crop
    real(dp), intent(in) :: accumulated_gdd

    stage = 0
    if (accumulated_gdd < crop%emergence_gdd) return

    stage = 1
    if (accumulated_gdd <= crop%initial_end_gdd) return

    if (crop%mid_start_gdd == huge(0._dp)) then
        stage = 4
        return
    end if

    stage = 2
    if (accumulated_gdd < crop%mid_start_gdd) return

    stage = 3
    if (accumulated_gdd <= crop%mid_end_gdd) return

    stage = 4
end function crop_growth_stage

subroutine set_bare_soil_today(pheno, i, j)
    type(crop_pars_matrices), intent(inout) :: pheno
    integer, intent(in) :: i, j

    pheno%crop_id(i,j) = 0
    pheno%pheno_stage(i,j) = 0
    pheno%is_real_crop(i,j) = .false.
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
end subroutine set_bare_soil_today

end module mod_cropcoef
