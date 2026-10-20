module mod_phenology_events
use mod_constants, only: dp
use mod_crop_phenology, only: crop_definition, crop_pars_matrices
use mod_date, only: date, days_in_year
use mod_grid, only: grid_i
use mod_utility, only: lower_case
implicit none
private

public :: phenology_event_list, daily_phenology_events, read_phenology_events, events_for_cell, check_bare_soil_events,       &
        & check_standing_crop_events, wait_for_observed_sowing, count_observed_cuts, crop_event_policy, &
        & record_observed_cut, apply_observed_harvest

integer, parameter :: event_sowing = 1, event_cut = 2, event_harvest = 3

type phenology_event_list
    integer :: n = 0
    integer, allocatable :: row(:), col(:), year(:), doy(:), kind(:)
    integer, allocatable :: first_event(:,:), next_event(:)
end type phenology_event_list

type daily_phenology_events
    logical :: sow = .false.
    logical :: cut = .false.
    logical :: harvest = .false.
end type daily_phenology_events

contains

! Read cell-level dates; the active rotation determines which crop each event concerns.
subroutine read_phenology_events(file_name, events, n_rows, n_cols)
    character(len=*), intent(in) :: file_name
    type(phenology_event_list), intent(out) :: events
    integer, intent(in) :: n_rows, n_cols
    character(len=500) :: line
    character(len=32) :: event_name
    character(len=32) :: extra_field
    integer :: unit_id, ios, extra_ios, line_number, n_records, record, row, col, year, doy, kind, e

    if (len_trim(file_name) == 0) return
    open(newunit=unit_id, file=trim(file_name), status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print *, 'Cannot open phenology events file: ', trim(file_name)
        stop
    end if

    n_records = 0
    do
        read(unit_id, '(A)', iostat=ios) line
        if (ios < 0) exit
        if (ios > 0) exit
        call clean_event_line(line)
        if (len_trim(line) == 0) cycle
        read(line, *, iostat=ios) row, col, year, doy, event_name
        if (ios /= 0 .and. index(line, 'row') == 1) cycle
        n_records = n_records + 1
    end do
    rewind(unit_id)
    allocate(events%row(n_records), events%col(n_records), events%year(n_records), events%doy(n_records), &
             events%kind(n_records))
    allocate(events%first_event(n_rows,n_cols), events%next_event(n_records))
    events%first_event = 0

    line_number = 0
    record = 0
    do
        read(unit_id, '(A)', iostat=ios) line
        if (ios < 0) exit
        line_number = line_number + 1
        if (ios > 0) call event_file_error(file_name, line_number, 'cannot read line')
        call clean_event_line(line)
        if (len_trim(line) == 0) cycle
        read(line, *, iostat=ios) row, col, year, doy, event_name
        if (ios /= 0 .and. index(line, 'row') == 1) cycle
        if (ios /= 0) call event_file_error(file_name, line_number, 'expected row col year doy event')
        read(line, *, iostat=extra_ios) row, col, year, doy, event_name, extra_field
        if (extra_ios == 0) call event_file_error(file_name, line_number, 'unexpected field after event')
        call lower_case(event_name)
        select case (trim(event_name))
            case ('sowing', 'sow')
                kind = event_sowing
            case ('cut', 'cutting')
                kind = event_cut
            case ('harvest', 'harvesting')
                kind = event_harvest
            case default
                call event_file_error(file_name, line_number, 'unknown event kind')
        end select
        if (row < 1 .or. row > n_rows .or. col < 1 .or. col > n_cols) &
            call event_file_error(file_name, line_number, 'cell is outside the grid')
        if (year < 1) call event_file_error(file_name, line_number, 'invalid year')
        if (doy < 1 .or. doy > days_in_year(year)) &
            call event_file_error(file_name, line_number, 'invalid day of year')
        e = events%first_event(row,col)
        do while (e > 0)
            if (events%year(e) == year .and. events%doy(e) == doy .and. events%kind(e) == kind) &
                call event_file_error(file_name, line_number, 'duplicate event for this cell and day')
            e = events%next_event(e)
        end do
        record = record + 1
        events%row(record) = row
        events%col(record) = col
        events%year(record) = year
        events%doy(record) = doy
        events%kind(record) = kind
        events%next_event(record) = events%first_event(row,col)
        events%first_event(row,col) = record
    end do
    close(unit_id)
    events%n = record
    print *, 'External phenology events loaded: ', events%n
end subroutine read_phenology_events

subroutine clean_event_line(line)
    character(len=*), intent(inout) :: line
    integer :: comment

    comment = scan(line, '#!')
    if (comment > 0) line(comment:) = ' '
    line = adjustl(line)
end subroutine clean_event_line

subroutine event_file_error(file_name, line_number, message)
    character(len=*), intent(in) :: file_name, message
    integer, intent(in) :: line_number

    print *, trim(file_name), ':', line_number, ': ', trim(message)
    print *, 'Execution will be aborted...'
    stop
end subroutine event_file_error

integer function event_on_day(events, row, col, year, doy, kind) result(event_idx)
    type(phenology_event_list), intent(in) :: events
    integer, intent(in) :: row, col, year, doy, kind
    integer :: e

    event_idx = 0
    if (events%n == 0) return
    e = events%first_event(row,col)
    do while (e > 0)
        if (events%year(e) == year .and. events%doy(e) == doy .and. events%kind(e) == kind) then
            event_idx = e
            return
        end if
        e = events%next_event(e)
    end do
end function event_on_day

logical function later_event_exists(events, row, col, year, doy, kind, latest_year, before_next_sowing) result(found)
    type(phenology_event_list), intent(in) :: events
    integer, intent(in) :: row, col, year, doy, kind, latest_year
    logical, intent(in) :: before_next_sowing
    integer :: e, current_day, event_day, next_sowing_day, next_target_day

    found = .false.
    if (events%n == 0) return
    current_day = year*367 + doy
    next_sowing_day = huge(0)
    next_target_day = huge(0)
    e = events%first_event(row,col)
    do while (e > 0)
        event_day = events%year(e)*367 + events%doy(e)
        if (event_day > current_day .and. events%year(e) <= latest_year) then
            if (events%kind(e) == kind) next_target_day = min(next_target_day, event_day)
            if (before_next_sowing .and. events%kind(e) == event_sowing) &
                next_sowing_day = min(next_sowing_day, event_day)
        end if
        e = events%next_event(e)
    end do
    found = next_target_day < next_sowing_day
end function later_event_exists

subroutine events_for_cell(events, i, j, current, today)
    type(phenology_event_list), intent(in) :: events
    integer, intent(in) :: i, j
    type(date), intent(in) :: current
    type(daily_phenology_events), intent(out) :: today
    integer :: e

    today = daily_phenology_events()
    if (events%n == 0) return
    e = event_on_day(events, i, j, current%year, current%doy, event_sowing)
    today%sow = e > 0
    e = event_on_day(events, i, j, current%year, current%doy, event_cut)
    today%cut = e > 0
    e = event_on_day(events, i, j, current%year, current%doy, event_harvest)
    today%harvest = e > 0
end subroutine events_for_cell

subroutine check_bare_soil_events(today, i, j, current)
    type(daily_phenology_events), intent(in) :: today
    integer, intent(in) :: i, j
    type(date), intent(in) :: current

    if (today%cut .or. today%harvest) &
        call invalid_crop_event(i, j, current, 'cut or harvest targets bare soil')
end subroutine check_bare_soil_events

subroutine check_standing_crop_events(today, i, j, current)
    type(daily_phenology_events), intent(in) :: today
    integer, intent(in) :: i, j
    type(date), intent(in) :: current

    if (today%sow) call invalid_crop_event(i, j, current, 'sowing overlaps a standing crop')
end subroutine check_standing_crop_events

logical function wait_for_observed_sowing(events, i, j, current) result(wait)
    type(phenology_event_list), intent(in) :: events
    integer, intent(in) :: i, j
    type(date), intent(in) :: current

    wait = .false.
    if (events%n == 0) return
    wait = later_event_exists(events, i, j, current%year, current%doy, event_sowing, current%year, .false.)
end function wait_for_observed_sowing

! Select one cut method for the occurrence at sowing. A later sowing or harvest bounds its observed cuts.
integer function count_observed_cuts(events, i, j, current, crop) result(n_cuts)
    type(phenology_event_list), intent(in) :: events
    integer, intent(in) :: i, j
    type(date), intent(in) :: current
    type(crop_definition), intent(in) :: crop
    integer :: e, event_day, sowing_day, next_sowing_day, harvest_day, latest_year

    n_cuts = 0
    if (events%n == 0) return

    sowing_day = current%year*367 + current%doy
    latest_year = current%year
    if (crop%harvest_doy_max <= current%doy) latest_year = latest_year + 1
    next_sowing_day = huge(0)
    harvest_day = huge(0)

    e = events%first_event(i,j)
    do while (e > 0)
        event_day = events%year(e)*367 + events%doy(e)
        if (event_day > sowing_day .and. events%year(e) <= latest_year) then
            if (events%kind(e) == event_sowing) next_sowing_day = min(next_sowing_day, event_day)
            if (events%kind(e) == event_harvest) harvest_day = min(harvest_day, event_day)
        end if
        e = events%next_event(e)
    end do

    e = events%first_event(i,j)
    do while (e > 0)
        event_day = events%year(e)*367 + events%doy(e)
        if (events%kind(e) == event_cut .and. event_day >= sowing_day .and. &
            & event_day < next_sowing_day .and. event_day < harvest_day .and. events%year(e) <= latest_year) &
            n_cuts = n_cuts + 1
        e = events%next_event(e)
    end do

    if (n_cuts >= crop%max_cuts) call invalid_crop_event(i, j, current, 'too many observed cuts for this crop')
end function count_observed_cuts

subroutine crop_event_policy(events, pheno, crop, i, j, current, wait_harvest)
    type(phenology_event_list), intent(in) :: events
    type(crop_pars_matrices), intent(in) :: pheno
    type(crop_definition), intent(in) :: crop
    integer, intent(in) :: i, j
    type(date), intent(in) :: current
    logical, intent(out) :: wait_harvest
    integer :: latest_year

    wait_harvest = .false.
    if (events%n == 0) return

    latest_year = pheno%sowing_year(i,j)
    if (crop%harvest_doy_max <= pheno%sowing_doy(i,j)) latest_year = latest_year + 1
    wait_harvest = later_event_exists(events, i, j, current%year, current%doy, event_harvest, latest_year, .true.)
end subroutine crop_event_policy

subroutine record_observed_cut(pheno, crop, i, j, current, projected_days)
    type(crop_pars_matrices), intent(inout) :: pheno
    type(crop_definition), intent(in) :: crop
    integer, intent(in) :: i, j
    type(date), intent(in) :: current
    integer, intent(in) :: projected_days(:)
    integer :: n_points, first_point, last_point

    if (crop%max_cuts <= 1 .or. pheno%cuts_completed(i,j) >= crop%max_cuts - 1) &
        call invalid_crop_event(i, j, current, 'crop has no remaining cut')

    pheno%cuts_completed(i,j) = pheno%cuts_completed(i,j) + 1
    pheno%cut_gdd_offset(i,j) = pheno%gdd(i,j) - real(pheno%cuts_completed(i,j), dp)*maxval(crop%gdd)
    n_points = size(crop%gdd)
    first_point = pheno%cuts_completed(i,j)*n_points + 1
    last_point = crop%max_cuts*n_points
    pheno%curve_day(i,j,first_point:last_point) = projected_days(1:last_point-first_point+1)
    where (pheno%curve_day(i,j,first_point:last_point) /= huge(0)) &
        pheno%curve_day(i,j,first_point:last_point) = &
            & pheno%curve_day(i,j,first_point:last_point) + pheno%crop_age_days(i,j)
end subroutine record_observed_cut

subroutine apply_observed_harvest(pheno, i, j)
    type(crop_pars_matrices), intent(inout) :: pheno
    integer, intent(in) :: i, j

    pheno%harvest_pending(i,j) = .true.
end subroutine apply_observed_harvest

subroutine invalid_crop_event(i, j, current, message)
    integer, intent(in) :: i, j
    type(date), intent(in) :: current
    character(len=*), intent(in) :: message

    print *, 'Phenology event at row/col/year/doy ', i, j, current%year, current%doy, ': ', trim(message)
    print *, 'Execution will be aborted...'
    stop
end subroutine invalid_crop_event

end module mod_phenology_events
