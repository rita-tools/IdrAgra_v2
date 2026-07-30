module mod_phenology_events
use mod_constants, only: dp
use mod_utility, only: lower_case, days_in_calendar_year
use mod_grid, only: grid_i
use mod_crop_phenology, only: crop_pheno_info, crop_matrices
implicit none

integer, parameter :: event_sowing = 1
integer, parameter :: event_cut = 2
integer, parameter :: event_harvest = 3

type phenology_event_list
    integer :: n = 0
    integer, dimension(:), allocatable :: row
    integer, dimension(:), allocatable :: col
    integer, dimension(:), allocatable :: year
    integer, dimension(:), allocatable :: doy
    integer, dimension(:), allocatable :: kind
    integer, dimension(:), allocatable :: crop_slot
end type phenology_event_list

contains

! Import events from their input file
subroutine read_phenology_events(file_name, events, n_rows, n_cols, n_slots)
    character(len=*), intent(in) :: file_name
    type(phenology_event_list), intent(inout) :: events
    integer, intent(in) :: n_rows, n_cols, n_slots

    integer :: unit_id, ios, line_number, n_records, record_idx
    integer :: row, col, year, doy, crop_slot, kind
    character(len=500) :: buffer
    character(len=32) :: event_name
    logical :: file_exists

    call destroy_phenology_events(events)
    if (len_trim(file_name) == 0) return
    inquire(file=trim(file_name), exist=file_exists)
    if (.not. file_exists) return

    open(newunit=unit_id, file=trim(file_name), status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print *, 'Cannot open phenology events file ', trim(file_name), '; iostat=', ios
        print *, 'Execution will be aborted...'
        stop
    end if

    n_records = 0
    line_number = 0
    do
        read(unit_id, '(A)', iostat=ios) buffer
        if (ios < 0) exit
        line_number = line_number + 1
        if (ios > 0) call event_file_error(file_name, line_number, 'cannot read line')
        call clean_event_line(buffer)
        if (len_trim(buffer) == 0) cycle
        read(buffer, *, iostat=ios) row, col, year, doy, event_name, crop_slot
        if (ios /= 0) then
            call lower_case(buffer)
            if (index(buffer, 'row' ) > 0 .and. index(buffer, 'col') > 0 .and. &
              & index(buffer, 'year') > 0 .and. index(buffer, 'doy') > 0) cycle
            call event_file_error(file_name, line_number, 'expected: row col year doy event crop_slot')
        end if
        n_records = n_records + 1
    end do

    rewind(unit_id)
    if (n_records > 0) then
        allocate(events%row(n_records), events%col(n_records), events%year(n_records))
        allocate(events%doy(n_records), events%kind(n_records), events%crop_slot(n_records))
    end if

    record_idx = 0
    line_number = 0
    do
        read(unit_id, '(A)', iostat=ios) buffer
        if (ios < 0) exit
        line_number = line_number + 1
        call clean_event_line(buffer)
        if (len_trim(buffer) == 0) cycle
        read(buffer, *, iostat=ios) row, col, year, doy, event_name, crop_slot
        if (ios /= 0) cycle ! header was validated during the first pass
        call lower_case(event_name)
        select case (trim(event_name))
            case ('sowing', 'sow')
                kind = event_sowing
            case ('cut', 'cutting')
                kind = event_cut
            case ('harvest', 'harvesting')
                kind = event_harvest
            case default
                call event_file_error(file_name, line_number, 'Accepted events are "sowing", "cut", or "harvest"')
        end select

        if (row < 1 .or. row > n_rows .or. col < 1 .or. col > n_cols) then
            call event_file_error(file_name, line_number, 'target row or column is outside the grid')
        end if
        if (crop_slot < 1 .or. crop_slot > n_slots) then
            call event_file_error(file_name, line_number, 'crop_slot is outside 1:n_crops')
        end if
        if (year < 1 .or. doy < 1 .or. doy > days_in_calendar_year(year)) then
            call event_file_error(file_name, line_number, 'invalid year/day-of-year')
        end if

        record_idx = record_idx + 1
        events%row(record_idx) = row
        events%col(record_idx) = col
        events%year(record_idx) = year
        events%doy(record_idx) = doy
        events%kind(record_idx) = kind
        events%crop_slot(record_idx) = crop_slot
    end do
    close(unit_id)
    events%n = record_idx

    call validate_duplicate_events(events, file_name)
    print *, 'External phenology events loaded: ', events%n, ' from ', trim(file_name)
end subroutine read_phenology_events

! Modify CropCoef-prescribed phenology according to cell-level sow/harvest dates specified by the user for this year
! For "cut" events, allocate a cut schedule with pre-calculated useful variables
subroutine apply_phenology_events(events, current_year, year_length, domain, soil_use, ws_idx, irandom, info_pheno, crop_mat)
    type(phenology_event_list), intent(in) :: events
    integer, intent(in) :: current_year, year_length
    type(grid_i), intent(in) :: domain, soil_use
    integer, dimension(:,:), intent(in) :: ws_idx, irandom
    type(crop_pheno_info), dimension(:), intent(in) :: info_pheno
    type(crop_matrices), intent(inout) :: crop_mat

    integer :: e, i, j, cycle_idx, slot, lu, ws, max_cuts, cut_position, effective_shift
    integer :: old_start, old_end, local_length, ref_length
    integer :: first_cut, regrowth_start, regrowth_end
    integer :: duration_warning_count, displacement_warning_count, ignored_cut_count
    logical :: has_year_events

    crop_mat%external_calendar = .false.
    crop_mat%suppress_irandom = .false.
    crop_mat%n_external_cuts = 0
    if (associated(crop_mat%external_cut_doy)) deallocate(crop_mat%external_cut_doy)
    crop_mat%ref_first_cut = 0
    crop_mat%ref_regrowth_start = 0
    crop_mat%ref_regrowth_end = 0

    ! Early return for years without events
    has_year_events = .false.
    do e=1, events%n
        if (events%year(e) == current_year) has_year_events = .true.
    end do
    if (.not. has_year_events) return

    ! Look for sowing and harvest events; if any, suppress irandom for the target cell.
    do e=1, events%n
        if (events%year(e) /= current_year) cycle
        if (events%doy(e) > year_length) then
            print *, 'Phenology event DOY exceeds the current simulation-year length: ', &
                & events%row(e), events%col(e), events%year(e), events%doy(e), year_length
            print *, 'Execution will be aborted...'
            stop
        end if
        if (events%kind(e) == event_cut) cycle
        call resolve_event_cycle(events, e, domain, soil_use, ws_idx, info_pheno, cycle_idx)
        i = events%row(e)
        j = events%col(e)
        crop_mat%external_calendar(i,j,cycle_idx) = .true.
        crop_mat%suppress_irandom(i,j) = .true.
        if (events%kind(e) == event_sowing) crop_mat%ii0(i,j,cycle_idx) = events%doy(e)
        if (events%kind(e) == event_harvest) crop_mat%iie(i,j,cycle_idx) = events%doy(e)
    end do

    ! Modify crop_mat to account for the new 
    duration_warning_count = 0
    displacement_warning_count = 0
    do j=1,size(domain%mat,2)
        do i=1,size(domain%mat,1)
            if (domain%mat(i,j) == domain%header%nan) cycle
            do cycle_idx=1,size(crop_mat%ii0,3)

                ! Recalculate thermal-stress windows using the new in-field period
                if (crop_mat%suppress_irandom(i,j) .and. crop_mat%ii0(i,j,cycle_idx) > 0 .and. &
                    & crop_mat%iie(i,j,cycle_idx) > 0) then
                    local_length = cycle_length(crop_mat%ii0(i,j,cycle_idx), crop_mat%iie(i,j,cycle_idx), year_length)
                    crop_mat%TSP_low(i,j,cycle_idx) = crop_mat%ii0(i,j,cycle_idx) + nint(0.45D0*local_length)
                    crop_mat%TSP_high(i,j,cycle_idx) = crop_mat%ii0(i,j,cycle_idx) + nint(0.75D0*local_length)
                    if (crop_mat%TSP_low(i,j,cycle_idx) > year_length) &
                        & crop_mat%TSP_low(i,j,cycle_idx) = crop_mat%TSP_low(i,j,cycle_idx)-year_length
                    if (crop_mat%TSP_high(i,j,cycle_idx) > year_length) &
                        & crop_mat%TSP_high(i,j,cycle_idx) = crop_mat%TSP_high(i,j,cycle_idx)-year_length
                end if

                ! Warn the user if the new cycle is now too short/long or shifted too much
                if (.not. crop_mat%external_calendar(i,j,cycle_idx)) cycle
                old_start = crop_mat%ii0_ref(i,j,cycle_idx)
                old_end = crop_mat%iie_ref(i,j,cycle_idx)
                local_length = cycle_length(crop_mat%ii0(i,j,cycle_idx), crop_mat%iie(i,j,cycle_idx), year_length)
                ref_length = nint(crop_mat%iid_ref(i,j,cycle_idx))
                if (local_length < 2) then
                    print *, 'Invalid external sowing/harvest interval at row/col/cycle ', i, j, cycle_idx
                    print *, 'Sowing/harvest DOY: ', crop_mat%ii0(i,j,cycle_idx), crop_mat%iie(i,j,cycle_idx)
                    print *, 'Execution will be aborted...'
                    stop
                end if
                crop_mat%iid(i,j,cycle_idx) = dble(local_length)
                if (ref_length > 0) then
                    crop_mat%dij(i,j,cycle_idx) = dble(ref_length)/dble(local_length)
                    if (dble(local_length)/dble(ref_length) < 0.70D0 .or. &
                        & dble(local_length)/dble(ref_length) > 1.30D0) then
                        duration_warning_count = duration_warning_count + 1
                        if (duration_warning_count <= 10) then
                            print *, 'Warning: external crop duration differs by more than 30% at row/col/cycle ', &
                                   & i, j, cycle_idx, '; reference/external days: ', ref_length, local_length
                        end if
                    end if
                end if
                if (calendar_midpoint_distance(crop_mat%ii0(i,j,cycle_idx), crop_mat%iie(i,j,cycle_idx), &
                    & old_start, old_end, year_length) > 30) then
                    displacement_warning_count = displacement_warning_count + 1
                    if (displacement_warning_count <= 10) then
                        print *, 'Warning: external crop-cycle midpoint moves by more than 30 days at row/col/cycle ', &
                            & i, j, cycle_idx
                    end if
                end if
            end do
        end do
    end do
    if (duration_warning_count > 10) print *, 'Additional duration warnings omitted: ', duration_warning_count-10
    if (displacement_warning_count > 10) &
        & print *, 'Additional calendar-displacement warnings omitted: ', displacement_warning_count-10
    call validate_cycle_overlaps(domain, soil_use, ws_idx, info_pheno, crop_mat, year_length)

    ! Allocate an indexed per-cell cut schedule for the current year.
    max_cuts = 0
    do e=1,events%n
        if (events%year(e) /= current_year .or. events%kind(e) /= event_cut) cycle
        call resolve_event_cycle(events, e, domain, soil_use, ws_idx, info_pheno, cycle_idx)
        i = events%row(e)
        j = events%col(e)
        crop_mat%n_external_cuts(i,j,cycle_idx) = crop_mat%n_external_cuts(i,j,cycle_idx)+1
        max_cuts = max(max_cuts, crop_mat%n_external_cuts(i,j,cycle_idx))
    end do
    if (max_cuts == 0) return
    allocate(crop_mat%external_cut_doy(size(domain%mat,1),size(domain%mat,2), size(crop_mat%ii0,3),max_cuts))
    crop_mat%external_cut_doy = 0
    crop_mat%n_external_cuts = 0
    do e=1,events%n
        if (events%year(e) /= current_year .or. events%kind(e) /= event_cut) cycle
        call resolve_event_cycle(events, e, domain, soil_use, ws_idx, info_pheno, cycle_idx)
        i = events%row(e)
        j = events%col(e)
        cut_position = crop_mat%n_external_cuts(i,j,cycle_idx)+1
        crop_mat%external_cut_doy(i,j,cycle_idx,cut_position) = events%doy(e)
        crop_mat%n_external_cuts(i,j,cycle_idx) = cut_position
    end do

    ignored_cut_count = 0
    do j=1,size(domain%mat,2)
        do i=1,size(domain%mat,1)
            if (domain%mat(i,j) == domain%header%nan) cycle
            lu = soil_use%mat(i,j)
            ws = ws_idx(i,j)
            do cycle_idx=1,size(crop_mat%ii0,3)
                if (crop_mat%n_external_cuts(i,j,cycle_idx) == 0) cycle
                call sort_integer_prefix(crop_mat%external_cut_doy(i,j,cycle_idx,:), crop_mat%n_external_cuts(i,j,cycle_idx))
                effective_shift = irandom(i,j)
                if (crop_mat%suppress_irandom(i,j)) effective_shift = 0
                call discard_out_of_cycle_cuts(crop_mat%external_cut_doy(i,j,cycle_idx,:), &
                    & crop_mat%n_external_cuts(i,j,cycle_idx), &
                    & max(1,min(year_length,crop_mat%ii0(i,j,cycle_idx)+effective_shift)), &
                    & max(1,min(year_length,crop_mat%iie(i,j,cycle_idx)+effective_shift)), ignored_cut_count)
                if (crop_mat%n_external_cuts(i,j,cycle_idx) == 0) cycle
                slot = info_pheno(ws)%cycle_crop_slot(lu,cycle_idx)
                call detect_reference_regrowth(info_pheno(ws), lu, slot, year_length, first_cut, regrowth_start, regrowth_end)
                if (first_cut == 0) then
                    print *, 'Warning: external cuts ignored at row/col/crop_slot ', i, j, slot, &
                        & ' because the reference Kcb series has no complete post-cut interval.'
                    ignored_cut_count = ignored_cut_count + crop_mat%n_external_cuts(i,j,cycle_idx)
                    crop_mat%n_external_cuts(i,j,cycle_idx) = 0
                    cycle
                end if
                crop_mat%ref_first_cut(i,j,cycle_idx) = first_cut
                crop_mat%ref_regrowth_start(i,j,cycle_idx) = regrowth_start
                crop_mat%ref_regrowth_end(i,j,cycle_idx) = regrowth_end
            end do
        end do
    end do
    if (ignored_cut_count > 0) print *, 'External cut events ignored for this year: ', ignored_cut_count
end subroutine apply_phenology_events

subroutine destroy_phenology_events(events)
    type(phenology_event_list), intent(inout) :: events
    if (allocated(events%row))       deallocate(events%row)
    if (allocated(events%col))       deallocate(events%col)
    if (allocated(events%year))      deallocate(events%year)
    if (allocated(events%doy))       deallocate(events%doy)
    if (allocated(events%kind))      deallocate(events%kind)
    if (allocated(events%crop_slot)) deallocate(events%crop_slot)
    events%n = 0
end subroutine destroy_phenology_events

subroutine make_cut_irrigation_halt_mask(crop_mat, info_pheno, ws_idx, soil_use, irandom, doy, &
                                       & year_length, halt_days, halt_mask)
    type(crop_matrices), intent(in) :: crop_mat
    type(crop_pheno_info), dimension(:), intent(in) :: info_pheno
    integer, dimension(:,:), intent(in) :: ws_idx, irandom
    type(grid_i), intent(in) :: soil_use
    integer, intent(in) :: doy, year_length, halt_days
    logical, dimension(:,:), intent(out) :: halt_mask
    integer :: i, j, cycle_idx, cut_idx, candidate_doy
    integer :: previous_cycle, current_cycle, previous_ref_doy, current_ref_doy
    integer :: lu, ws, slot, reference_doy

    halt_mask = .false.
    if (halt_days <= 0) return

    do j=1,size(halt_mask,2)
        do i=1,size(halt_mask,1)
            if (soil_use%mat(i,j) == soil_use%header%nan) cycle
            if (ws_idx(i,j) < 1 .or. ws_idx(i,j) > size(info_pheno)) cycle
            ! Externally prescribed cuts supersede the ordinary Kcb-inferred
            ! cuts for their crop cycle.
            crop_cycles: do cycle_idx=1,size(crop_mat%n_external_cuts,3)
                do cut_idx=1,crop_mat%n_external_cuts(i,j,cycle_idx)
                    if (abs(doy-crop_mat%external_cut_doy(i,j,cycle_idx,cut_idx)) <= halt_days) then
                        halt_mask(i,j) = .true.
                        exit crop_cycles
                    end if
                end do
            end do crop_cycles
            if (halt_mask(i,j)) cycle

            lu = soil_use%mat(i,j)
            ws = ws_idx(i,j)
            do candidate_doy=max(2,doy-halt_days),min(year_length,doy+halt_days)
                call map_local_to_reference_day(crop_mat, info_pheno(ws), i, j, lu, irandom(i,j), &
                    & candidate_doy-1, year_length, previous_cycle, previous_ref_doy)
                call map_local_to_reference_day(crop_mat, info_pheno(ws), i, j, lu, irandom(i,j), &
                    & candidate_doy, year_length, current_cycle, current_ref_doy)
                if (current_cycle == 0 .or. current_cycle /= previous_cycle) cycle
                if (crop_mat%n_external_cuts(i,j,current_cycle) > 0) cycle
                if (current_ref_doy <= previous_ref_doy) cycle
                slot = info_pheno(ws)%cycle_crop_slot(lu,current_cycle)
                do reference_doy=previous_ref_doy+1,current_ref_doy
                    if (is_reference_cut(info_pheno(ws),lu,slot,reference_doy, &
                        & crop_mat%k_cb_min(i,j,current_cycle),crop_mat%k_cb_max(i,j,current_cycle))) then
                        halt_mask(i,j) = .true.
                        exit
                    end if
                end do
                if (halt_mask(i,j)) exit
            end do
        end do
    end do
end subroutine make_cut_irrigation_halt_mask

subroutine map_local_to_reference_day(crop_mat, pheno, i, j, lu, irandom, doy, year_length, &
                                    & active_cycle, reference_doy)
    type(crop_matrices), intent(in) :: crop_mat
    type(crop_pheno_info), intent(in) :: pheno
    integer, intent(in) :: i, j, lu, irandom, doy, year_length
    integer, intent(out) :: active_cycle, reference_doy
    integer :: cycle_idx, effective_shift, ii0, iie, ref_ii0, ref_iie
    real(dp) :: autumn_scale

    effective_shift = irandom
    if (crop_mat%suppress_irandom(i,j)) effective_shift = 0
    active_cycle = 0
    reference_doy = 0
    do cycle_idx=1,size(crop_mat%ii0,3)
        if (pheno%cycle_crop_slot(lu,cycle_idx) <= 0) cycle
        ii0 = max(1,min(year_length,crop_mat%ii0(i,j,cycle_idx)+effective_shift))
        iie = max(1,min(year_length,crop_mat%iie(i,j,cycle_idx)+effective_shift))
        if (crop_mat%ii0(i,j,cycle_idx) <= crop_mat%iie(i,j,cycle_idx)) then
            if (doy >= ii0 .and. doy <= iie) active_cycle = cycle_idx
        else
            if (doy <= iie .or. doy >= ii0) active_cycle = cycle_idx
        end if
        if (active_cycle > 0) exit
    end do
    if (active_cycle == 0) return

    ii0 = max(1,min(year_length,crop_mat%ii0(i,j,active_cycle)+effective_shift))
    iie = max(1,min(year_length,crop_mat%iie(i,j,active_cycle)+effective_shift))
    ref_ii0 = crop_mat%ii0_ref(i,j,active_cycle)
    ref_iie = crop_mat%iie_ref(i,j,active_cycle)
    if (crop_mat%ii0(i,j,active_cycle) <= crop_mat%iie(i,j,active_cycle)) then
        if (iie > ii0) then
            reference_doy = ref_ii0+nint(dble(doy-ii0)*dble(ref_iie-ref_ii0)/dble(iie-ii0))
        else
            reference_doy = ref_ii0
        end if
    else if (doy <= iie) then
        reference_doy = nint(dble(doy)*ref_iie/dble(max(1,iie)))
    else
        autumn_scale = dble(year_length-ref_ii0)/dble(max(1,year_length-ii0))
        reference_doy = ref_ii0+nint(dble(doy-ii0)*autumn_scale)
    end if
    reference_doy = max(1,min(year_length,reference_doy))
end subroutine map_local_to_reference_day

logical function is_reference_cut(pheno, lu, slot, day, low_kcb, high_kcb)
    type(crop_pheno_info), intent(in) :: pheno
    integer, intent(in) :: lu, slot, day
    real(dp), intent(in) :: low_kcb, high_kcb
    real(dp) :: amplitude

    is_reference_cut = .false.
    if (day < 2 .or. day > size(pheno%k_cb%tab,1)) return
    if (pheno%crop_slot%tab(day,lu) /= slot .or. pheno%crop_slot%tab(day-1,lu) /= slot) return
    amplitude = high_kcb-low_kcb
    if (amplitude <= tiny(1.0D0)) return
    if (pheno%k_cb%tab(day-1,lu) < low_kcb+0.80D0*amplitude) return
    is_reference_cut = pheno%k_cb%tab(day-1,lu)-pheno%k_cb%tab(day,lu) > 0.20D0*amplitude
end function is_reference_cut

subroutine clean_event_line(buffer)
    character(len=*), intent(inout) :: buffer
    integer :: p
    buffer = adjustl(buffer)
    p = scan(buffer, '#')
    if (p > 0) buffer = buffer(:p-1)
    buffer = trim(buffer)
end subroutine clean_event_line

subroutine event_file_error(file_name, line_number, message)
    character(len=*), intent(in) :: file_name, message
    integer, intent(in) :: line_number
    print *, 'Invalid phenology event at line ', line_number, ' of ', trim(file_name)
    print *, trim(message)
    print *, 'Execution will be aborted...'
    stop
end subroutine event_file_error

subroutine validate_duplicate_events(events, file_name)
    type(phenology_event_list), intent(in) :: events
    character(len=*), intent(in) :: file_name
    integer :: a, b
    do a=1,events%n
        do b=a+1,events%n
            if (events%row(a) /= events%row(b) .or. events%col(a) /= events%col(b) .or. &
              & events%year(a) /= events%year(b) .or. events%crop_slot(a) /= events%crop_slot(b)) cycle
            if (events%kind(a) == events%kind(b) .and. (events%doy(a) == events%doy(b) .or. events%kind(a) /= event_cut)) then
                print *, 'Duplicate/conflicting external phenology events in ', trim(file_name)
                print *, 'Target/year/crop_slot: ', events%row(a), events%col(a), events%year(a), events%crop_slot(a)
                print *, 'Execution will be aborted...'
                stop
            end if
        end do
    end do
end subroutine validate_duplicate_events

subroutine resolve_event_cycle(events, e, domain, soil_use, ws_idx, info_pheno, cycle_idx)
    type(phenology_event_list), intent(in) :: events
    integer, intent(in) :: e
    type(grid_i), intent(in) :: domain, soil_use
    integer, dimension(:,:), intent(in) :: ws_idx
    type(crop_pheno_info), dimension(:), intent(in) :: info_pheno
    integer, intent(out) :: cycle_idx
    integer :: i, j, lu, ws, z

    i = events%row(e)
    j = events%col(e)
    if (domain%mat(i,j) == domain%header%nan) then
        print *, 'Phenology event targets an inactive cell at row/col/year: ', i, j, events%year(e)
        print *, 'Execution will be aborted...'
        stop
    end if
    lu = soil_use%mat(i,j)
    ws = ws_idx(i,j)
    cycle_idx = 0
    do z=1,size(info_pheno(ws)%cycle_crop_slot,2)
        if (info_pheno(ws)%cycle_crop_slot(lu,z) == events%crop_slot(e)) then
            cycle_idx = z
            exit
        end if
    end do
    if (cycle_idx == 0) then
        print *, 'Phenology event crop_slot is not present at row/col/year: ', &
            & i, j, events%year(e), '; crop_slot=', events%crop_slot(e)
        print *, 'Execution will be aborted...'
        stop
    end if
end subroutine resolve_event_cycle

subroutine validate_cycle_overlaps(domain, soil_use, ws_idx, info_pheno, crop_mat, year_length)
    type(grid_i), intent(in) :: domain, soil_use
    integer, dimension(:,:), intent(in) :: ws_idx
    type(crop_pheno_info), dimension(:), intent(in) :: info_pheno
    type(crop_matrices), intent(in) :: crop_mat
    integer, intent(in) :: year_length
    logical, dimension(year_length) :: occupied
    integer :: i, j, z, day, lu, ws

    do j=1,size(domain%mat,2)
        do i=1,size(domain%mat,1)
            if (domain%mat(i,j) == domain%header%nan) cycle
            occupied = .false.
            lu = soil_use%mat(i,j)
            ws = ws_idx(i,j)
            do z=1,size(crop_mat%ii0,3)
                if (info_pheno(ws)%cycle_crop_slot(lu,z) <= 0) cycle
                do day=1,year_length
                    if (.not. day_in_cycle(day,crop_mat%ii0(i,j,z),crop_mat%iie(i,j,z))) cycle
                    if (occupied(day)) then
                        print *, 'External phenology calendar creates overlapping crop cycles at row/col/day ', i, j, day
                        print *, 'Execution will be aborted...'
                        stop
                    end if
                    occupied(day) = .true.
                end do
            end do
        end do
    end do
end subroutine validate_cycle_overlaps

subroutine detect_reference_regrowth(pheno, lu, slot, n_days, first_cut, regrowth_start, regrowth_end)
    type(crop_pheno_info), intent(in) :: pheno
    integer, intent(in) :: lu, slot, n_days
    integer, intent(out) :: first_cut, regrowth_start, regrowth_end
    integer :: day, second_cut
    real(dp) :: low_kcb, high_kcb, amplitude, drop_threshold, mature_threshold
    logical, dimension(n_days) :: slot_mask

    first_cut = 0
    second_cut = 0
    regrowth_start = 0
    regrowth_end = 0
    slot_mask = pheno%crop_slot%tab(:,lu) == slot
    if (.not. any(slot_mask)) return
    low_kcb = minval(pheno%k_cb%tab(:,lu), mask=slot_mask)
    high_kcb = maxval(pheno%k_cb%tab(:,lu), mask=slot_mask)
    amplitude = high_kcb-low_kcb
    if (amplitude <= tiny(1.0D0)) return
    drop_threshold = 0.20D0*amplitude
    mature_threshold = low_kcb+0.80D0*amplitude

    do day=2,n_days
        if (.not. slot_mask(day) .or. .not. slot_mask(day-1)) cycle
        if (pheno%k_cb%tab(day-1,lu) < mature_threshold) cycle
        if (pheno%k_cb%tab(day-1,lu)-pheno%k_cb%tab(day,lu) <= drop_threshold) cycle
        if (first_cut == 0) then
            first_cut = day
        else
            second_cut = day
            exit
        end if
    end do
    if (second_cut == 0) then
        first_cut = 0
        return
    end if
    regrowth_start = first_cut
    regrowth_end = second_cut-1
end subroutine detect_reference_regrowth

subroutine discard_out_of_cycle_cuts(cuts, n_cuts, sowing, harvest, ignored_count)
    integer, dimension(:), intent(inout) :: cuts
    integer, intent(inout) :: n_cuts, ignored_count
    integer, intent(in) :: sowing, harvest
    integer :: source_idx, target_idx
    target_idx = 0
    do source_idx=1,n_cuts
        if (day_in_cycle(cuts(source_idx),sowing,harvest)) then
            target_idx = target_idx+1
            cuts(target_idx) = cuts(source_idx)
        else
            print *, 'Warning: cut DOY ', cuts(source_idx), ' is outside the active sowing/harvest interval and will be ignored.'
            ignored_count = ignored_count+1
        end if
    end do
    n_cuts = target_idx
end subroutine discard_out_of_cycle_cuts

subroutine sort_integer_prefix(values, n)
    integer, dimension(:), intent(inout) :: values
    integer, intent(in) :: n
    integer :: i, j, value
    do i=2,n
        value = values(i)
        j = i-1
        do while (j >= 1)
            if (values(j) <= value) exit
            values(j+1) = values(j)
            j = j-1
        end do
        values(j+1) = value
    end do
end subroutine sort_integer_prefix

pure logical function day_in_cycle(day, sowing, harvest)
    integer, intent(in) :: day, sowing, harvest
    if (sowing <= harvest) then
        day_in_cycle = day >= sowing .and. day <= harvest
    else
        day_in_cycle = day >= sowing .or. day <= harvest
    end if
end function day_in_cycle

pure integer function cycle_length(sowing, harvest, year_length)
    integer, intent(in) :: sowing, harvest, year_length
    if (sowing <= harvest) then
        cycle_length = harvest-sowing+1
    else
        cycle_length = year_length-sowing+1+harvest
    end if
end function cycle_length

pure integer function calendar_midpoint_distance(start_a, end_a, start_b, end_b, year_length)
    integer, intent(in) :: start_a, end_a, start_b, end_b, year_length
    integer :: midpoint_a, midpoint_b, direct_distance
    midpoint_a = modulo(start_a-1+(cycle_length(start_a,end_a,year_length)-1)/2,year_length)+1
    midpoint_b = modulo(start_b-1+(cycle_length(start_b,end_b,year_length)-1)/2,year_length)+1
    direct_distance = abs(midpoint_a-midpoint_b)
    calendar_midpoint_distance = min(direct_distance,year_length-direct_distance)
end function calendar_midpoint_distance

end module mod_phenology_events
