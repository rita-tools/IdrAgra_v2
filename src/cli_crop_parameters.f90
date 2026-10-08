module cli_crop_parameters
use mod_constants, only: dp, nan_r
use mod_utility, only: clean_input_line, join_path, lower_case, string_to_integers, string_to_reals, split_string, &
                     & count_elements, replace_str
use mod_parameters, only: simulation
use mod_meteo, only: meteo_info
use mod_crop_phenology
use mod_system
implicit none


logical, dimension(:), allocatable, save :: missing_crop_slot_warned

interface read_crop_pars
    module procedure read_crop_pars_r, read_crop_pars_i
end interface

interface init_daily_crop_par_file
    module procedure init_daily_crop_par_file_r, init_daily_crop_par_file_i
end interface

interface spread_col
    module procedure spread_col_i, spread_col_r
end interface

contains

! Read all of the crop .dat files and the crop rotation file, and store the information in crop_definitions and crop_rotations.
subroutine import_crop_definitions(sim, crop_definitions, crop_rotations)
    type(simulation), intent(in) :: sim
    type(crop_definition), dimension(:), allocatable, intent(out) :: crop_definitions
    type(crop_rotation), dimension(:), allocatable, intent(out) :: crop_rotations

    character(len=255), dimension(:), allocatable :: crop_files
    character(len=255), dimension(:,:), allocatable :: rotation_files
    integer, dimension(:), allocatable :: crop_ids
    character(len=1000) :: line
    character(len=255) :: first_crop, second_crop
    character(len=500) :: rotations_file
    integer :: unit, ios, land_use_id, expected_land_use_id
    integer :: rotation_idx, crop_idx, file_idx, n_rotations, n_crop_files
    integer :: explicit_id, maximum_crop_id

    rotations_file = join_path(sim%crop_inputs_path, sim%soil_uses_fn, delimiter)
    n_rotations = count_crop_rotations(rotations_file)

    allocate(crop_rotations(n_rotations))
    allocate(rotation_files(n_rotations, 2))
    allocate(crop_files(2*n_rotations), crop_ids(2*n_rotations))
    rotation_files = ''
    crop_files = ''
    crop_ids = 0
    n_crop_files = 0

    open(newunit=unit, file=trim(rotations_file), status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print *, 'Cannot open crop rotation file ', trim(rotations_file), '.'
        print *, 'Execution will be aborted...'
        stop
    end if

    rotation_idx = 0
    expected_land_use_id = 1
    do
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        line = clean_input_line(line)
        if (len_trim(line) == 0) cycle

        first_crop = ''
        second_crop = ''
        read(line, *, iostat=ios) land_use_id, first_crop, second_crop
        if (ios /= 0) cycle

        if (land_use_id /= expected_land_use_id) then
            print *, 'Crop rotation IDs must be consecutive integers starting from 1.'
            print *, 'Expected/found: ', expected_land_use_id, land_use_id, ' in ', trim(rotations_file)
            print *, 'Execution will be aborted...'
            stop
        end if

        rotation_idx = rotation_idx + 1
        crop_rotations(rotation_idx)%land_use_id = land_use_id
        rotation_files(rotation_idx, 1) = trim(first_crop)
        if (trim(second_crop) /= '*') rotation_files(rotation_idx, 2) = trim(second_crop)

        do crop_idx = 1, 2
            if (len_trim(rotation_files(rotation_idx, crop_idx)) == 0) cycle
            file_idx = find_crop_file(crop_files, n_crop_files, rotation_files(rotation_idx, crop_idx))
            if (file_idx == 0) then
                n_crop_files = n_crop_files + 1
                crop_files(n_crop_files) = rotation_files(rotation_idx, crop_idx)
            end if
        end do
        expected_land_use_id = expected_land_use_id + 1
    end do
    close(unit)

    ! Numeric filename prefixes are stable crop IDs. Unprefixed files receive
    ! deterministic IDs after the largest explicit one.
    maximum_crop_id = 0
    do file_idx = 1, n_crop_files
        explicit_id = crop_id_from_filename(crop_files(file_idx))
        if (explicit_id <= 0) cycle
        crop_ids(file_idx) = explicit_id
        maximum_crop_id = max(maximum_crop_id, explicit_id)
    end do

    do file_idx = 1, n_crop_files
        if (crop_ids(file_idx) > 0) cycle
        maximum_crop_id = maximum_crop_id + 1
        crop_ids(file_idx) = maximum_crop_id
    end do

    allocate(crop_definitions(maximum_crop_id))
    do file_idx = 1, n_crop_files
        crop_idx = crop_ids(file_idx)
        crop_definitions(crop_idx)%crop_id = crop_idx
        crop_definitions(crop_idx)%parameter_file = trim(crop_files(file_idx))
        call read_crop_file(join_path(sim%crop_parameters_path, crop_files(file_idx), delimiter), &
                          & crop_definitions(crop_idx))
    end do

    do rotation_idx = 1, n_rotations
        crop_idx = 1
        if (len_trim(rotation_files(rotation_idx, 2)) > 0) crop_idx = 2
        allocate(crop_rotations(rotation_idx)%crop_ids(crop_idx))
        do file_idx = 1, crop_idx
            crop_rotations(rotation_idx)%crop_ids(file_idx) = &
                & crop_ids(find_crop_file(crop_files, n_crop_files, rotation_files(rotation_idx, file_idx)))
        end do
    end do
end subroutine import_crop_definitions

integer function count_crop_rotations(file_name)
    character(len=*), intent(in) :: file_name
    character(len=1000) :: line
    character(len=255) :: first_crop, second_crop
    integer :: unit, ios, land_use_id

    count_crop_rotations = 0
    open(newunit=unit, file=trim(file_name), status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print *, 'Cannot open crop rotation file ', trim(file_name), '.'
        print *, 'Execution will be aborted...'
        stop
    end if

    do
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        line = clean_input_line(line)
        read(line, *, iostat=ios) land_use_id, first_crop, second_crop
        if (ios == 0) count_crop_rotations = count_crop_rotations + 1
    end do
    close(unit)
end function count_crop_rotations

integer function find_crop_file(crop_files, n_crop_files, file_name)
    character(len=255), dimension(:), intent(in) :: crop_files
    integer, intent(in) :: n_crop_files
    character(len=*), intent(in) :: file_name
    character(len=255) :: stored_name, requested_name
    integer :: idx

    requested_name = trim(file_name)
    call lower_case(requested_name)
    find_crop_file = 0
    do idx = 1, n_crop_files
        stored_name = trim(crop_files(idx))
        call lower_case(stored_name)
        if (trim(stored_name) == trim(requested_name)) then
            find_crop_file = idx
            return
        end if
    end do
end function find_crop_file

integer function crop_id_from_filename(file_name)
    character(len=*), intent(in) :: file_name
    character(len=255) :: base_name
    integer :: separator, underscore, ios

    separator = max(scan(trim(file_name), '/', back=.true.), scan(trim(file_name), achar(92), back=.true.))
    base_name = file_name(separator + 1:)
    underscore = scan(trim(base_name), '_')
    crop_id_from_filename = 0
    if (underscore <= 1) return
    read(base_name(:underscore - 1), *, iostat=ios) crop_id_from_filename
    if (ios /= 0 .or. crop_id_from_filename <= 0) crop_id_from_filename = 0
end function crop_id_from_filename

subroutine read_crop_file(file_name, crop)
    character(len=*), intent(in) :: file_name
    type(crop_definition), intent(inout) :: crop
    character(len=2000) :: line, value, lower_line
    character(len=100) :: label
    real(dp), dimension(9) :: row
    integer :: unit, ios, separator, binary_flag
    logical :: reading_table

    crop%crop_name = crop_name_from_filename(file_name)
    crop%parameter_file = trim(file_name)
    allocate(crop%gdd(0), crop%k_cb(0), crop%lai(0), crop%height(0), crop%root_depth(0), &
           & crop%yield_response(0), crop%cn_value(0), crop%cover_fraction(0), crop%stress_resistance(0))

    open(newunit=unit, file=trim(file_name), status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print *, 'Cannot open crop parameter file ', trim(file_name), '.'
        print *, 'Execution will be aborted...'
        stop
    end if

    reading_table = .false.
    do
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        line = clean_input_line(line)
        if (len_trim(line) == 0) cycle

        lower_line = trim(line)
        call lower_case(lower_line)
        if (trim(lower_line) == 'endtable') then
            reading_table = .false.
            cycle
        end if
        if (index(adjustl(lower_line), 'gdd') == 1) then
            reading_table = .true.
            cycle
        end if

        if (reading_table) then
            row = nan_r
            value = replace_str(trim(line), '*', '-9999')
            read(value, *, iostat=ios) row
            if (row(1) == nan_r) then
                print *, 'Invalid GDD curve row in ', trim(file_name), ': ', trim(line)
                print *, 'Execution will be aborted...'
                stop
            end if
            call append_crop_curve_row(crop, row)
            cycle
        end if

        separator = scan(line, '=')
        if (separator == 0) cycle
        label = adjustl(trim(line(:separator - 1)))
        value = adjustl(trim(line(separator + 1:)))
        call lower_case(label)
        select case (trim(label))
            case ('sowingdate_min'); read(value, *, iostat=ios) crop%sowing_doy_min
            case ('sowingdelay_max'); read(value, *, iostat=ios) crop%sowing_delay_max
            case ('harvestdate_max'); read(value, *, iostat=ios) crop%harvest_doy_max
            case ('harvnum_max'); read(value, *, iostat=ios) crop%max_cuts
            case ('cropsoverlap'); read(value, *, iostat=ios) crop%crop_overlap_days
            case ('tsowing'); read(value, *, iostat=ios) crop%sowing_temp
            case ('tdaybase'); read(value, *, iostat=ios) crop%base_temp
            case ('tcutoff'); read(value, *, iostat=ios) crop%cutoff_temp
            case ('vern'); read(value, *, iostat=ios) binary_flag; if (ios == 0) crop%requires_vernalization = binary_flag /= 0
            case ('tv_min'); read(value, *, iostat=ios) crop%vern_temp_min
            case ('tv_max'); read(value, *, iostat=ios) crop%vern_temp_max
            case ('vfmin'); read(value, *, iostat=ios) crop%vern_fact_min
            case ('vstart'); read(value, *, iostat=ios) crop%vern_days_start
            case ('vend'); read(value, *, iostat=ios) crop%vern_days_end
            case ('vslope'); read(value, *, iostat=ios) crop%vern_curve_slope
            case ('ph_r'); read(value, *, iostat=ios) crop%photoperiod_response
            case ('daylength_if'); read(value, *, iostat=ios) crop%daylength_if
            case ('daylength_ins'); read(value, *, iostat=ios) crop%daylength_ins
            case ('wp'); read(value, *, iostat=ios) crop%water_productivity
            case ('fsink'); read(value, *, iostat=ios) crop%sink_strength
            case ('tcrit_hs'); read(value, *, iostat=ios) crop%heat_stress_temp_crit
            case ('tlim_hs'); read(value, *, iostat=ios) crop%heat_stress_temp_lim
            case ('hi'); read(value, *, iostat=ios) crop%harvest_index
            case ('kyt'); read(value, *, iostat=ios) crop%yield_response_total
            case ('ky1'); read(value, *, iostat=ios) crop%yield_response_stage(1)
            case ('ky2'); read(value, *, iostat=ios) crop%yield_response_stage(2)
            case ('ky3'); read(value, *, iostat=ios) crop%yield_response_stage(3)
            case ('ky4'); read(value, *, iostat=ios) crop%yield_response_stage(4)
            case ('praw'); read(value, *, iostat=ios) crop%raw_fraction
            case ('ainterception'); read(value, *, iostat=ios) crop%interception_coef
            case ('cl_cn'); read(value, *, iostat=ios) crop%cn_class
            case ('irrigation'); read(value, *, iostat=ios) binary_flag; if (ios == 0) crop%is_irrigated = binary_flag /= 0
            case ('adj_flag'); read(value, *, iostat=ios) binary_flag; if (ios == 0) crop%adjust_k_cb = binary_flag /= 0
            case ('ke'); read(value, *, iostat=ios) crop%evaporative_layer_root_fraction
            case ('kt'); read(value, *, iostat=ios) crop%transpirative_layer_root_fraction
            case ('rft'); read(value, *, iostat=ios) crop%maximum_transpirative_root_fraction
            case default
                print *, 'Skipping invalid or obsolete label <',trim(label),'> in ', trim(file_name)
                cycle
        end select

        if (ios /= 0) then
            print *, 'Invalid value for ', trim(label), ' in ', trim(file_name), '.'
            print *, 'Execution will be aborted...'
            stop
        end if
    end do
    close(unit)

    call prepare_crop_definition_curves(crop)
end subroutine read_crop_file

subroutine append_crop_curve_row(crop, row)
    type(crop_definition), intent(inout) :: crop
    real(dp), dimension(9), intent(in) :: row

    crop%gdd = [crop%gdd, row(1)]
    crop%k_cb = [crop%k_cb, row(2)]
    crop%lai = [crop%lai, row(3)]
    crop%height = [crop%height, row(4)]
    crop%root_depth = [crop%root_depth, row(5)]
    crop%cn_value = [crop%cn_value, row(6)]
    crop%cover_fraction = [crop%cover_fraction, row(7)]
    crop%stress_resistance = [crop%stress_resistance, row(8)]
    crop%yield_response = [crop%yield_response, row(9)]
end subroutine append_crop_curve_row

subroutine prepare_crop_definition_curves(crop)
    ! Match CropCoef's preprocessing of missing curve values so the numerical
    ! module can consume complete parameter curves directly.
    type(crop_definition), intent(inout) :: crop
    integer :: idx

    if (size(crop%gdd) == 0) return

    call fill_missing_curve(crop%k_cb, crop%gdd)
    call fill_missing_curve(crop%lai, crop%gdd)
    call fill_missing_curve(crop%height, crop%gdd)
    call fill_missing_curve(crop%root_depth, crop%gdd)
    call fill_missing_curve(crop%yield_response, crop%gdd)
    call fill_missing_curve(crop%cover_fraction, crop%gdd)

    if (all(crop%stress_resistance == nan_r)) then
        crop%stress_resistance = 0.0_dp
    else
        call fill_missing_curve(crop%stress_resistance)
    end if

    do idx = 1, size(crop%cn_value)
        if (crop%cn_value(idx) /= nan_r) cycle
        crop%cn_value(idx) = 1.0_dp
        if (crop%k_cb(idx) == 0.0_dp) crop%cn_value(idx) = 0.0_dp
        if (crop%k_cb(idx) >= 0.45_dp) crop%cn_value(idx) = 2.0_dp
    end do

    crop%yield_response(1) = crop%yield_response_stage(1)
    crop%yield_response(size(crop%yield_response)) = crop%yield_response_stage(4)
    do idx = 1, size(crop%yield_response)
        if (crop%yield_response(idx) /= nan_r) cycle
        if (crop%k_cb(idx) < 0.45_dp) then
            crop%yield_response(idx) = crop%yield_response_stage(2)
        else
            crop%yield_response(idx) = crop%yield_response_stage(3)
        end if
    end do
end subroutine prepare_crop_definition_curves

! Analogous to CropCoef's fillMissingL()
subroutine fill_missing_curve(values, abscissa)
    real(dp), dimension(:), intent(inout) :: values
    real(dp), dimension(:), optional, intent(in) :: abscissa
    real(dp), dimension(size(values)) :: x
    integer, dimension(size(values)) :: known
    integer :: idx, pair_idx, start_idx, end_idx, n_known
    real(dp) :: slope

    if (present(abscissa)) then
        x = abscissa
    else
        x = [(real(idx, dp), idx=1,size(values))]
    end if

    known = 0
    n_known = 0
    do idx = 1, size(values)
        if (values(idx) == nan_r) cycle
        n_known = n_known + 1
        known(n_known) = idx
    end do
    do pair_idx = 1, n_known - 1
        start_idx = known(pair_idx)
        end_idx = known(pair_idx + 1)
        if (end_idx <= start_idx + 1) cycle
        if (x(end_idx) == x(start_idx)) then
            slope = 0._dp
        else
            slope = (values(end_idx) - values(start_idx))/(x(end_idx) - x(start_idx))
        end if
        do idx = start_idx + 1, end_idx - 1
            values(idx) = values(start_idx) + slope*(x(idx) - x(start_idx))
        end do
    end do
end subroutine fill_missing_curve

function crop_name_from_filename(file_name) result(crop_name)
    character(len=*), intent(in) :: file_name
    character(len=255) :: crop_name
    integer :: separator, extension

    separator = max(scan(trim(file_name), '/', back=.true.), scan(trim(file_name), achar(92), back=.true.))
    crop_name = file_name(separator + 1:)
    extension = scan(trim(crop_name), '.', back=.true.)
    if (extension > 1) crop_name = crop_name(:extension - 1)
end function crop_name_from_filename

subroutine init_daily_crop_par_file_r(file_pars, file_name)
    ! Store and validate a real-valued daily crop parameter file.
    character(len=*), intent(in) :: file_name
    type(file_phenology_r), intent(out) :: file_pars

    call init_daily_crop_par_file_common(file_pars%filename, file_pars%next_pos, file_name)
end subroutine init_daily_crop_par_file_r

subroutine init_daily_crop_par_file_i(file_pars, file_name)
    ! Store and validate an integer-valued daily crop parameter file.
    character(len=*), intent(in) :: file_name
    type(file_phenology_i), intent(out) :: file_pars

    call init_daily_crop_par_file_common(file_pars%filename, file_pars%next_pos, file_name)
end subroutine init_daily_crop_par_file_i

subroutine init_daily_crop_par_file_common(stored_name, next_pos, file_name)
    ! Validate the file, skip its header, remember the first data position, and close it.
    character(len=*), intent(out) :: stored_name
    integer, intent(out) :: next_pos
    character(len=*), intent(in) :: file_name
    integer :: unit, ios
    character(len=500) :: io_message

    stored_name = trim(file_name)

    open(newunit=unit, file=trim(stored_name), status='old', action='read', &
       & access='stream', form='formatted', iostat=ios, iomsg=io_message    )
    if (ios /= 0) then
        print *, 'Cannot open phenology file ', trim(stored_name), ': ', trim(io_message)
        print *, 'Execution will be aborted...'
        stop
    end if

    read(unit, '(A)', iostat=ios, iomsg=io_message)
    if (ios /= 0) then
        print *, 'Cannot read the header of phenology file ', trim(stored_name), ': ', trim(io_message)
        print *, 'Execution will be aborted...'
        stop
    end if

    ! Save the position in the file as next_pos (read_crop_pars() will start reading data from here)
    inquire(unit=unit, pos=next_pos)
    close(unit)
end subroutine init_daily_crop_par_file_common

subroutine init_crop_par_from_file(file_name, n_crop, n_crop_alt, string_elements, n_crops_by_year)
    ! init static crop parameters from parameter file
    character(len=*), intent(in) :: file_name
    integer, intent(in) :: n_crop
    integer, intent(out) :: n_crop_alt
    integer, intent(out) :: string_elements
    integer, dimension(n_crop), intent(out) :: n_crops_by_year
    integer :: free_unit, ios, p
    integer, dimension(:), allocatable :: crop_counts
    character(len=n_crop*20) :: buffer, label !EAC: use mcrop_max x 20
    character(len=10), dimension(:), allocatable :: dummy, dummy_clean

    open(newunit=free_unit, file=trim(file_name), status='old', action="read", iostat=ios)
    if (ios /= 0 ) then
        print *, "Cannot open file ", trim(file_name), ". The specified file does not exist. Execution will be aborted..."
        stop
    end if

    do while (ios == 0)
        read (free_unit, '(A)', iostat=ios) buffer
        if (ios == 0) then
            call lower_case(buffer)
            p = scan(buffer, achar(9))  ! find the first tab ---> tab=achar(9)
            label = buffer(1:p-1)
            buffer = buffer(p+1:)

            select case (label)
                case ('var')
                    ! read the header line and divide elements
                    allocate(dummy(n_crop*2))    ! TODO: to add 3 crops, edit 2 to 3
                    call split_string(buffer, achar(9), dummy, string_elements)
                    allocate(dummy_clean(string_elements))
                    dummy_clean = dummy(1:string_elements)
                    deallocate(dummy)
                    call count_elements(dummy_clean, crop_counts)
                    deallocate(dummy_clean)
                    if (size(crop_counts) > n_crop) then
                        print *, "Invalid crop parameters in ", trim(file_name), "."
                        print *, "The header contains ", size(crop_counts), " unique crop IDs, but SoilUsesNum is ", n_crop, "."
                        print *, "Set SoilUsesNum to at least the number of unique crop IDs in idragra_parameters.txt."
                        stop
                    end if
                    n_crops_by_year = 0
                    n_crops_by_year(:size(crop_counts)) = crop_counts
                    deallocate(crop_counts)
                    n_crop_alt = maxval(n_crops_by_year)
               case default
            end select
        end if
    end do
    close (free_unit)
end subroutine init_crop_par_from_file

subroutine read_water_prod_file(file_name, string_elements, n_crops_by_year, &
                                unit_param, sim_end_year, weath_start_year)
    ! read water productivity related parameters
    character(len=*), intent(in) :: file_name
    integer, intent(in) :: string_elements
    integer, dimension(:), intent(in) :: n_crops_by_year
    real(dp), dimension(:,:,:), intent(inout) :: unit_param
    integer, intent(in) :: sim_end_year, weath_start_year
    integer :: free_unit, ios, line, p, year
    character(len=string_elements*20) :: buffer, label  !EAC:  use string_elements x 20

    line = 0
    open(newunit=free_unit, file=trim(file_name), status='old', action="read", iostat=ios)
    if (ios /= 0 ) then
        print *, "Cannot open file ", trim(file_name), ". The specified file does not exist. &
            & Execution will be aborted..."
        stop
    end if

    do while (ios == 0)
        read (free_unit, '(A)', iostat=ios) buffer
        if (ios == 0) then
            line = line + 1
            buffer = trim(buffer)
            call lower_case(buffer)
            p = scan(buffer, achar(9))  ! find the first tab ---> tab=achar(9)
            label = buffer(1:p-1)
            buffer = buffer(p+1:)

            select case (label)
                case ('year') ! header line
                case default
                    read (label, *, iostat=ios) year
                    if (ios /= 0) then
                        print *, 'Invalid WPadj.dat entry at line ', line, ': ', trim(label), '. Execution will be aborted...'
                        stop
                    end if
                    ! Store values from the first year of weather data to the last year of simulation (weather data might start earlier than simulation)
                    if (year >= weath_start_year .and. year <= sim_end_year) then
                        call spread_col(buffer,achar(9),string_elements,n_crops_by_year,unit_param(:,:,year-weath_start_year+1))
                    end if
            end select
        end if
    end do
    close (free_unit)
end subroutine read_water_prod_file

subroutine read_canopy_resistance_file(file_name, unit_param, sim_end_year, weath_start_year)
    ! read canopy resistance parameters
    character(len=*), intent(in) :: file_name
    real(dp), dimension(:), intent(inout) :: unit_param
    integer, intent(in) :: sim_end_year, weath_start_year
    integer :: free_unit, ios, line, year
    real(dp) :: resistance
    character(len=255) :: full_line

    open(newunit=free_unit, file=trim(file_name), status='old', action="read", iostat=ios)
    if (ios /= 0 ) then
        print *, "Cannot open file ", trim(file_name), ". The specified file does not exist. Execution will be aborted..."
        stop
    end if

    read (free_unit, '(A)', iostat=ios) ! skip the first line
    line = 1
    do while (ios == 0)
        read (free_unit, '(A)', iostat=ios) full_line
        if (ios == 0) then
            line = line + 1
            read (full_line, *, iostat=ios) year, resistance
            if (ios /= 0) then
                print *, 'Invalid CanopyRes.dat entry at line ', line, ': ', trim(full_line), '. Execution will be aborted...'
                stop
            end if
            ! Store values from the first year of weather data to the last year of simulation (weather data might start earlier than simulation)
            if (year >= weath_start_year .and. year <= sim_end_year) then
                unit_param(year - weath_start_year + 1) = resistance
            end if
        end if
    end do

    close (free_unit)
end subroutine read_canopy_resistance_file

subroutine spread_col_i(string_in, sep, string_el, string_space, string_out)
    character(len=*), intent(in) :: string_in
    character(len=*), intent(in) :: sep
    integer, intent(in) :: string_el
    integer, dimension(:), intent(in) :: string_space
    integer, dimension(:,:), intent(inout) :: string_out
    integer, dimension(:), allocatable :: dummy
    integer :: i

    allocate(dummy(string_el))
    dummy = string_to_integers(string_in(1:len_trim(string_in)-1), sep)
    do i=1, size(string_space)
       string_out(i,1:string_space(i)) = &
            & dummy(sum(string_space(1:i-1))+1:sum(string_space(1:i)))
    end do
end subroutine spread_col_i

subroutine spread_col_r(string_in, sep, string_el, string_space, string_out)
    character(len=*), intent(in) :: string_in
    character(len=*), intent(in) :: sep
    integer, intent(in) :: string_el
    integer, dimension(:), intent(in) :: string_space
    real(dp), dimension(:,:), intent(inout) :: string_out
    real(dp), dimension(:), allocatable :: dummy
    integer :: i

    allocate(dummy(string_el))
    dummy = string_to_reals(string_in(1:len_trim(string_in)-1), sep)

    do i=1, size(string_space)
       string_out(i,1:string_space(i)) = &
            & dummy(sum(string_space(1:i-1))+1:sum(string_space(1:i)))
    end do
end subroutine spread_col_r

subroutine read_crop_par_file(file_name, string_elements, ze_fix, unit_param)
    ! read static crop parameters file
    ! TODO: merge with init_crop_par_from_file ?
    character(len=*), intent(in) :: file_name
    integer, intent(in) :: string_elements
    real(dp), intent(in) :: ze_fix
    type(crop_pheno_info) :: unit_param
    integer :: free_unit
    integer :: ios
    integer :: line, p
    character(len=string_elements*20) :: buffer, label !EAC: use string_elements x 20

    open(newunit=free_unit, file=trim(file_name), status='old', action="read", iostat=ios)
    if (ios /= 0 ) then
        print *, "Cannot open file ", trim(file_name), ". The specified file does not exist. &
            & The CropCoef version you used to generate inputs might be outdated. Execution will be aborted..."
        stop
    end if

    do while (ios == 0)
        read (free_unit, '(A)', iostat=ios) buffer
        if (ios == 0) then
            line = line + 1
            buffer = trim(buffer)
            call lower_case(buffer)
            p = scan(buffer, achar(9))  ! find the first tab ---> tab=achar(9)
            label = buffer(1:p-1)
            buffer = buffer(p+1:)

            select case (label)
                case ('var') ! already initialized
                case ('irrig')
                    call spread_col(buffer, achar(9), string_elements, unit_param%n_crops_by_year, unit_param%irrigation_class)
                case ('cnclass')
                    call spread_col(buffer, achar(9), string_elements, unit_param%n_crops_by_year, unit_param%cn_class)
                case('praw')
                     call spread_col(buffer, achar(9), string_elements, unit_param%n_crops_by_year, unit_param%p_raw_const)
                case('aint')
                    call spread_col(buffer, achar(9), string_elements, unit_param%n_crops_by_year, unit_param%a)
                case('tlim')
                    call spread_col(buffer, achar(9), string_elements, unit_param%n_crops_by_year, unit_param%T_lim)
                case('tcrit')
                    call spread_col(buffer, achar(9), string_elements, unit_param%n_crops_by_year, unit_param%T_crit)
                case('hi')
                    call spread_col(buffer, achar(9), string_elements, unit_param%n_crops_by_year, unit_param%HI)
                case('kyt')
                    call spread_col(buffer, achar(9), string_elements, unit_param%n_crops_by_year, unit_param%Ky_tot)
                case('ky1')
                    call spread_col(buffer, achar(9), string_elements, unit_param%n_crops_by_year, unit_param%Ky_pheno(:,:,1))
                case('ky2')
                    call spread_col(buffer, achar(9), string_elements, unit_param%n_crops_by_year, unit_param%Ky_pheno(:,:,2))
                case('ky3')
                    call spread_col(buffer, achar(9), string_elements, unit_param%n_crops_by_year, unit_param%Ky_pheno(:,:,3))
                case('ky4')
                    call spread_col(buffer, achar(9), string_elements, unit_param%n_crops_by_year, unit_param%Ky_pheno(:,:,4))
                case('rft')
                     call spread_col(buffer, achar(9), string_elements, unit_param%n_crops_by_year, unit_param%max_RF_t)
                case('maxsr')
                    call spread_col(buffer, achar(9), string_elements, unit_param%n_crops_by_year, unit_param%d_r_max)
                    unit_param%d_r_max = unit_param%d_r_max - ze_fix
                case default
                    print *, 'Skipping invalid or obsolete label <',trim(label),'> at line', line, ' of file: ', file_name
            end select
        end if
    end do
    close (free_unit)
end subroutine read_crop_par_file

subroutine read_crop_pars_r(file_pars,n_days,n_crop)
    ! read crop phenology series from the specified file (real values)
    integer,intent(in) :: n_days                    ! number of days (i.e. 365 o 366)
    integer,intent(in)::n_crop                      ! number of crops
    type(file_phenology_r), intent(inout) :: file_pars
    integer :: i, unit, ios
    character(len=500) :: io_message

    allocate(file_pars%tab(n_days,n_crop))

    ! Open the required file
    open(newunit=unit, file=trim(file_pars%filename), status='old', action='read', &
       & access='stream', form='formatted', iostat=ios, iomsg=io_message           )
    if (ios /= 0) then
        print *, 'Cannot open phenology file ', trim(file_pars%filename), ': ', trim(io_message)
        print *, 'Execution will be aborted...'
        stop
    end if

    ! Read one year of data starting from next_pos, i.e. the place where last year's reading terminated
    do i=1,size(file_pars%tab,1)
        if (i == 1) then
            read(unit, *, pos=file_pars%next_pos, iostat=ios, iomsg=io_message) file_pars%tab(i,:)
        else
            read(unit, *, iostat=ios, iomsg=io_message) file_pars%tab(i,:)
        end if
        if (ios /= 0) then
            print *, 'Cannot read day ', i, ' from phenology file ', trim(file_pars%filename), ': ', trim(io_message)
            print *, 'Execution will be aborted...'
            stop
        end if
    end do

    ! Save the current position of the cursor before closing the file 
    inquire(unit=unit, pos=file_pars%next_pos)
    close(unit)

end subroutine read_crop_pars_r

subroutine read_crop_pars_i(file_pars,n_days,n_crop)
    ! read crop phenology series from the specified file (int values)
    integer,intent(in) :: n_days                ! number of days (i.e. 365 o 366)
    integer,intent(in)::n_crop                  ! number of crops
    type(file_phenology_i), intent(inout) :: file_pars
    integer :: i, unit, ios
    character(len=500) :: io_message

    allocate(file_pars%tab(n_days,n_crop))

    ! Open the required file
    open(newunit=unit, file=trim(file_pars%filename), status='old', action='read', &
       & access='stream', form='formatted', iostat=ios, iomsg=io_message           )
    if (ios /= 0) then
        print *, 'Cannot open phenology file ', trim(file_pars%filename), ': ', trim(io_message)
        print *, 'Execution will be aborted...'
        stop
    end if

    ! Read one year of data starting from next_pos, i.e. the place where last year's reading terminated
    do i=1,size(file_pars%tab,1)
        if (i == 1) then
            read(unit, *, pos=file_pars%next_pos, iostat=ios, iomsg=io_message) file_pars%tab(i,:)
        else
            read(unit, *, iostat=ios, iomsg=io_message) file_pars%tab(i,:)
        end if
        if (ios /= 0) then
            print *, 'Cannot read day ', i, ' from phenology file ', trim(file_pars%filename), ': ', trim(io_message)
            print *, 'Execution will be aborted...'
            stop
        end if
    end do

    ! Save the current position of the cursor before closing the file 
    inquire(unit=unit, pos=file_pars%next_pos)
    close(unit)

end subroutine read_crop_pars_i

subroutine init_crop_phenology_pars(sim, info_pheno, info_meteo, ze_fix, verbose, last_year)
    ! init crop parameters and file references for daily parameters
    type(simulation),intent(inout)::sim
    type(meteo_info),dimension(:),intent(in)::info_meteo
    real(dp),intent(in) :: ze_fix
    logical,intent(in)::verbose
    integer, optional, intent(in) :: last_year

    type(crop_pheno_info),dimension(:),allocatable::info_pheno
    character(len=255)::dir,froot,dir_name
    integer :: i, string_elements, parameter_end_year
    integer, dimension(sim%n_lus) :: n_crops_by_year
    real(dp), parameter :: nan = -9999.0D0
    integer, parameter :: phases = 4

    dir= trim(sim%pheno_path)
    froot = sim%pheno_root
    parameter_end_year = sim%end%year
    if (present(last_year)) parameter_end_year = max(parameter_end_year, last_year)

    allocate(info_pheno(sim%n_weather_stations)) ! init to the number of weather stations

    allocate(sim%res_canopy(sim%meteo_years))
    sim%res_canopy=0
    ! init from crop parameters file (actually produced by cropcoeff)
    dir_name = info_meteo(1)%filename(1:(index(trim(info_meteo(1)%filename),"."))-1)  ! directory has the same name as the weather station dataset
    call init_crop_par_from_file(trim(dir)//trim(froot)//trim(dir_name)//delimiter//"CropParam.dat", &
        & sim%n_lus, sim%n_crops, string_elements, n_crops_by_year)

    call read_canopy_resistance_file(trim(dir)//delimiter//'CanopyRes.dat', sim%res_canopy, parameter_end_year, sim%start_year)

    do i=1,size(info_pheno)
        dir_name = info_meteo(i)%filename(1:(index(trim(info_meteo(i)%filename),"."))-1)
        call init_daily_crop_par_file(info_pheno(i)%k_cb,      trim(dir)//trim(froot)//trim(dir_name)//delimiter//"Kcb.dat")
        call init_daily_crop_par_file(info_pheno(i)%h,         trim(dir)//trim(froot)//trim(dir_name)//delimiter//"H.dat")
        call init_daily_crop_par_file(info_pheno(i)%z_r,       trim(dir)//trim(froot)//trim(dir_name)//delimiter//"Sr.dat")
        call init_daily_crop_par_file(info_pheno(i)%lai,       trim(dir)//trim(froot)//trim(dir_name)//delimiter//"LAI.dat")
        call init_daily_crop_par_file(info_pheno(i)%cn_day,    trim(dir)//trim(froot)//trim(dir_name)//delimiter//"CNvalue.dat")
        call init_daily_crop_par_file(info_pheno(i)%f_c,       trim(dir)//trim(froot)//trim(dir_name)//delimiter//"fc.dat")
        ! EDIT: add support for seasonal p_raw
        call init_daily_crop_par_file(info_pheno(i)%r_stress,  trim(dir)//trim(froot)//trim(dir_name)//delimiter//"r_stress.dat")
        call init_daily_crop_par_file(info_pheno(i)%crop_slot, trim(dir)//trim(froot)//trim(dir_name)//delimiter//"CropId.dat")

        ! TODO - add tabulated ky

        ! init crop parameters
        allocate(info_pheno(i)%n_crops_by_year    (sim%n_lus))
        allocate(info_pheno(i)%irrigation_class (sim%n_lus, sim%n_crops))
        allocate(info_pheno(i)%cn_class             (sim%n_lus, sim%n_crops))
        allocate(info_pheno(i)%p_raw_const     (sim%n_lus, sim%n_crops))
        allocate(info_pheno(i)%a              (sim%n_lus, sim%n_crops))
        allocate(info_pheno(i)%d_r_max          (sim%n_lus, sim%n_crops))
        allocate(info_pheno(i)%max_RF_t         (sim%n_lus, sim%n_crops))
        allocate(info_pheno(i)%T_lim           (sim%n_lus, sim%n_crops))
        allocate(info_pheno(i)%T_crit          (sim%n_lus, sim%n_crops))
        allocate(info_pheno(i)%HI             (sim%n_lus, sim%n_crops))
        allocate(info_pheno(i)%Ky_tot            (sim%n_lus, sim%n_crops))
        allocate(info_pheno(i)%Ky_pheno            (sim%n_lus, sim%n_crops, phases))
        allocate(info_pheno(i)%wp_adj          (sim%n_lus, sim%n_crops, sim%meteo_years))
        allocate(info_pheno(i)%kcb_phases%low (sim%n_lus, sim%n_crops))
        allocate(info_pheno(i)%kcb_phases%high(sim%n_lus, sim%n_crops))
        allocate(info_pheno(i)%kcb_phases%mid (sim%n_lus, sim%n_crops))
        allocate(info_pheno(i)%ii0            (sim%n_lus, sim%n_crops))
        allocate(info_pheno(i)%iie            (sim%n_lus, sim%n_crops))
        allocate(info_pheno(i)%iid            (sim%n_lus, sim%n_crops))
        allocate(info_pheno(i)%cycle_crop_slot(sim%n_lus, sim%n_crops))
        info_pheno(i)%n_crops_by_year     = n_crops_by_year
        info_pheno(i)%irrigation_class  = int(nan)
        info_pheno(i)%cn_class              = int(nan)
        info_pheno(i)%p_raw_const               = nan
        info_pheno(i)%a               = nan
        info_pheno(i)%d_r_max           = nan
        info_pheno(i)%max_RF_t          = nan
        info_pheno(i)%T_lim            = nan
        info_pheno(i)%T_crit           = nan
        info_pheno(i)%HI              = nan
        info_pheno(i)%Ky_tot             = nan
        info_pheno(i)%Ky_pheno             = nan
        info_pheno(i)%wp_adj           = nan
        info_pheno(i)%kcb_phases%low  = nan
        info_pheno(i)%kcb_phases%high = nan
        info_pheno(i)%kcb_phases%mid  = nan
        info_pheno(i)%ii0             = 0
        info_pheno(i)%iie             = 0
        info_pheno(i)%iid             = 0
        info_pheno(i)%cycle_crop_slot = 0
        call read_water_prod_file(trim(dir)//trim(froot)//trim(dir_name)//delimiter//"WPadj.dat",  &
                                & string_elements, n_crops_by_year, info_pheno(i)%wp_adj,          &
                                & parameter_end_year, sim%start_year)
        call read_crop_par_file(trim(dir)//trim(froot)//trim(dir_name)//delimiter//"CropParam.dat", &
                              & string_elements, ze_fix, info_pheno(i))

        ! TODO
        ! EAC: overwrite p values if exits
        ! call open_daily_crop_par_file(info_pheno(i)%p%unit,trim(dir)//trim(froot)//trim(fname)//"\praw.dat")

    end do
    if (verbose .eqv. .true.) then
        print *,'===== DEBUG: crop parameters initialize ====='
        print *,'path to phenophase: ',  dir
        print *,'root of subfolder: ',  froot
        print *,'# of phenophase: ',  size(info_meteo)
        print *,'===== END DEBUG ====='
    end if
end subroutine init_crop_phenology_pars

subroutine read_all_crop_pars(n_days, n_crop, info_pheno)
    integer,intent(in)::n_days,n_crop
    type(crop_pheno_info),dimension(:),intent(inout)::info_pheno
    integer::i

    do i=1,size(info_pheno) ! i.e. the number of weather stations
        call read_crop_pars(info_pheno(i)%k_cb,n_days,n_crop)
        call read_crop_pars(info_pheno(i)%h,n_days,n_crop)
        call read_crop_pars(info_pheno(i)%z_r,n_days,n_crop)
        call read_crop_pars(info_pheno(i)%lai,n_days,n_crop)
        call read_crop_pars(info_pheno(i)%cn_day,n_days,n_crop)
        call read_crop_pars(info_pheno(i)%f_c,n_days,n_crop)
        call read_crop_pars(info_pheno(i)%r_stress,n_days,n_crop)
        call read_crop_pars(info_pheno(i)%crop_slot,n_days,n_crop)
        call derive_crop_cycles(info_pheno(i), n_days)
    end do
end subroutine read_all_crop_pars

subroutine derive_crop_cycles(pheno, n_days)
    type(crop_pheno_info), intent(inout) :: pheno
    integer, intent(in) :: n_days
    integer :: lu, slot, cycle_idx, day, start_day, end_day, crop_slot, first_end, last_start
    integer :: n_slots, n_present_slots
    logical, dimension(n_days) :: crop_mask
    real(dp) :: low_value, high_value, mid_value

    pheno%ii0 = 0
    pheno%iie = 0
    pheno%iid = 0
    pheno%cycle_crop_slot = 0

    if (.not. allocated(missing_crop_slot_warned)) then
        allocate(missing_crop_slot_warned(size(pheno%crop_slot%tab,2)))
        missing_crop_slot_warned = .false.
    end if

    do lu=1,size(pheno%crop_slot%tab,2)
        n_slots = pheno%n_crops_by_year(lu)
        if (n_slots < 1) cycle

        ! CropId.dat uses 0 for bare soil and 1:n_slots for the declared rotation slots
        do day=1,n_days
            crop_slot = pheno%crop_slot%tab(day,lu)
            if (crop_slot < 0 .or. crop_slot > n_slots) then
                print *, 'Invalid crop slot ', crop_slot, ' at day ', day, ', land-use class ', lu
                print *, 'Expected a value between 0 and ', n_slots
                print *, 'Execution will be aborted...'
                stop
            end if
        end do

        ! Derive crop-specific daily time-series by rotation slot
        n_present_slots = 0
        do slot=1,n_slots
            crop_mask = pheno%crop_slot%tab(:,lu) == slot
            if (.not. any(crop_mask)) then
                if (.not. missing_crop_slot_warned(lu)) then
                    print *, 'Warning: Crop slot ', slot, ' never occurs in land-use class ', lu, &
                            &'. CropCoef likely overwrote an entire crop due to overlapping sow/harvest dates.'
                    missing_crop_slot_warned(lu) = .true.
                end if
                cycle
            end if
            n_present_slots = n_present_slots + 1
            low_value = minval(pheno%k_cb%tab(:,lu))
            high_value = maxval(pheno%k_cb%tab(:,lu), mask=crop_mask)
            mid_value = high_value
            do day=2,n_days
                if (crop_mask(day) .and. crop_mask(day-1)) then
                    if (pheno%k_cb%tab(day,lu) == pheno%k_cb%tab(day-1,lu) .and. &
                        & pheno%k_cb%tab(day,lu) > low_value .and. pheno%k_cb%tab(day,lu) < high_value) then
                        mid_value = pheno%k_cb%tab(day,lu)
                        exit
                    end if
                end if
            end do
            pheno%kcb_phases%low(lu,slot) = low_value
            pheno%kcb_phases%high(lu,slot) = high_value
            pheno%kcb_phases%mid(lu,slot) = mid_value
        end do

        cycle_idx = 0
        first_end = 0
        last_start = n_days + 1

        ! Matching nonzero slots at both year ends are treated as one crop cycle crossing New Year.
        if (pheno%crop_slot%tab(1,lu) > 0 .and. &
            & pheno%crop_slot%tab(1,lu) == pheno%crop_slot%tab(n_days,lu)) then
            crop_slot = pheno%crop_slot%tab(1,lu)
            first_end = 1
            do while (first_end < n_days .and. pheno%crop_slot%tab(first_end+1,lu) == crop_slot)
                first_end = first_end + 1
            end do
            cycle_idx = 1
            if (first_end == n_days) then
                last_start = n_days + 1
                pheno%ii0(lu,cycle_idx) = 1
                pheno%iie(lu,cycle_idx) = n_days
                pheno%iid(lu,cycle_idx) = n_days
            else
                last_start = n_days
                do while (last_start > 1 .and. pheno%crop_slot%tab(last_start-1,lu) == crop_slot)
                    last_start = last_start - 1
                end do
                pheno%ii0(lu,cycle_idx) = last_start
                pheno%iie(lu,cycle_idx) = first_end
                pheno%iid(lu,cycle_idx) = n_days-last_start+1+first_end
            end if
            pheno%cycle_crop_slot(lu,cycle_idx) = crop_slot
        end if

        ! Store the remaining contiguous crop periods as cycles in calendar order
        day = first_end + 1
        do while (day <= min(n_days,last_start-1))
            if (pheno%crop_slot%tab(day,lu) == 0) then
                day = day + 1
                cycle
            end if
            crop_slot = pheno%crop_slot%tab(day,lu)
            start_day = day
            do while (day <= min(n_days,last_start-1) .and. pheno%crop_slot%tab(day,lu) == crop_slot)
                day = day + 1
            end do
            end_day = day - 1
            cycle_idx = cycle_idx + 1
            if (cycle_idx > size(pheno%ii0,2)) then
                print *, 'Too many crop cycles in land-use class ', lu
                stop
            end if
            pheno%ii0(lu,cycle_idx) = start_day
            pheno%iie(lu,cycle_idx) = end_day
            pheno%iid(lu,cycle_idx) = end_day-start_day+1
            pheno%cycle_crop_slot(lu,cycle_idx) = crop_slot
        end do

        if (cycle_idx /= n_present_slots) then
            print *, 'CropId.dat defines ', cycle_idx, ' crop cycles for land-use class ', lu, &
                & ', but ', n_present_slots, ' declared slots occur in the daily series.'
            print *, 'Execution will be aborted...'
            stop
        end if
    end do
end subroutine derive_crop_cycles

subroutine destroy_infofeno_tab(info_pheno)
! dellaocate all crop phenological time series
    type(crop_pheno_info),dimension(:),intent(inout)::info_pheno
    integer::i

    do i=1,size(info_pheno)
        if(associated(info_pheno(i)%k_cb%tab)) deallocate(info_pheno(i)%k_cb%tab)
        if(associated(info_pheno(i)%h%tab)) deallocate(info_pheno(i)%h%tab)
        if(associated(info_pheno(i)%z_r%tab)) deallocate(info_pheno(i)%z_r%tab)
        if(associated(info_pheno(i)%lai%tab)) deallocate(info_pheno(i)%lai%tab)
        if(associated(info_pheno(i)%cn_day%tab)) deallocate(info_pheno(i)%cn_day%tab)
        if(associated(info_pheno(i)%f_c%tab)) deallocate(info_pheno(i)%f_c%tab)
        if(associated(info_pheno(i)%r_stress%tab)) deallocate(info_pheno(i)%r_stress%tab)
        if(associated(info_pheno(i)%crop_slot%tab)) deallocate(info_pheno(i)%crop_slot%tab)
    end do

end subroutine destroy_infofeno_tab

subroutine destroy_info_pheno(info_pheno)
    ! close all phenological opened files
    type(crop_pheno_info),dimension(:),allocatable,intent(inout)::info_pheno
    integer::i

    !%PS%: now pheno files are closed after each year, here we only need to deallocate info_pheno and its components
    !todo: because some components of info_pheno are pointers, deallocating it does not automatically free all of the memory up.
    !      Consider using allocatable instead of pointers or extending this subroutine to properly deallocate all pointers.

    do i=1,size(info_pheno)
        if(associated(info_pheno(i)%cycle_crop_slot)) deallocate(info_pheno(i)%cycle_crop_slot)
    end do
    deallocate(info_pheno)
end subroutine destroy_info_pheno

subroutine check_pheno_parameters(info_pheno,info_meteo)
! check if phenological parameters match weather station data
    type(crop_pheno_info),dimension(:),intent(in)::info_pheno
    type(meteo_info),dimension(:),intent(in)::info_meteo
    integer::i

    do i=1,size(info_pheno)
        call check_crop_parameters(info_pheno(i),info_meteo(i)%filename)
    end do
end subroutine check_pheno_parameters

subroutine check_crop_parameters(info_pheno,weather_station)
! check if phenological parameters match weather station data
! if k_cb is null than all the other parameters must be null
    type(crop_pheno_info),intent(in)::info_pheno
    character(len=*),intent(in)::weather_station
    integer::d,k

    do k=1, size(info_pheno%k_cb%tab,2) ! loop over crops
        do d=1,size(info_pheno%k_cb%tab,1)    ! loop over data
            if(info_pheno%k_cb%tab(d,k).gt.0)then
                if(info_pheno%z_r%tab(d,k).eq.0.) then
                    write(*,*) 'Warning: Sr null. Station: ', trim(weather_station), ' Day: ', d, ' Soil use class: ', k
                end if
                if(info_pheno%h%tab(d,k).eq.0.) then
                    write(*,*) 'Warning: H null. Station: ', trim(weather_station), ' Day: ', d, ' Soil use class: ', k
                end if
                if(info_pheno%LAI%tab(d,k).eq.0.) then
                    write(*,*) 'Warning: LAI null. Station: ', trim(weather_station), ' Day: ', d, ' Soil use class: ', k
                end if
            end if
        end do
    end do
end subroutine check_crop_parameters

end module cli_crop_parameters
