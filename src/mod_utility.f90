module mod_utility!

    use mod_constants, only: sp, dp, pi
    
    ! type that stores date
    type date!
        integer::day!
        integer::month!
        integer::year!
        integer::doy
        integer::weekday
    end type date!
    
    private get_value_index_i, get_value_index_c

    interface get_value_index
        module procedure get_value_index_i, get_value_index_c
    end interface

    interface split_date
        module procedure split_date_c, split_date_d
    end interface
        
    contains

    ! Count the frequency of unique values from vec
    subroutine count_element(vec,vec_el)
        implicit none
        character(len=*),dimension(:),intent(in) :: vec
        integer,dimension(:),intent(inout) :: vec_el
        
        integer :: i
        logical,dimension(size(vec)) :: mask
        integer,dimension(size(vec)) :: count_vec
 
        ! Select the unique elements
        mask(1)=.true.
        do i=size(vec),2,-1
            mask(i)= .not.(any(vec(:i-1)==vec(i)))
        end do

        ! Count the number of occurrences of elements
        do i=1,size(vec)
            count_vec(i) = count(vec(i)==vec)
        end do
        
        vec_el = pack(count_vec, mask)
    end subroutine count_element

    subroutine split_string(string, delimiter, substrings, substring_count)
        ! split a string into two sides of a delimiter token
        implicit none
        character(len=*), intent(in) :: string
        character, intent(in) :: delimiter
        character(len=*), intent(out) :: substrings(*)
        integer, intent(out) :: substring_count
        
        integer :: start_position, end_position
        
        start_position = 1
        substring_count = 0
        
        do
            end_position = index(string(start_position:), delimiter)
            substring_count = substring_count +1
            if (end_position == 0) then
                substrings(substring_count) = string(start_position:)
                exit
            else
                substrings(substring_count) = string(start_position : start_position + end_position - 2)
                start_position = start_position + end_position
            end if
        end do
    end subroutine split_string

    subroutine lower_case(word)
        ! convert a word to lower case
        character (len=*) , intent(inout) :: word
        integer :: i,ic,nlen
        nlen = len(word)
        if (nlen == 0) stop "Zero characters string"
        do i=1,nlen
            ic = ichar(word(i:i))
            if (ic >= 65 .and. ic <= 90) word(i:i) = char(ic+32)
        end do
    end subroutine lower_case

    subroutine seek_un(ErrorFlag, free_unit)
        ! Find the first free unit for reading and writing operation
        
        integer :: ErrorFlag    ! check for any error
        integer :: free_unit    ! number of the first available free unit from 11 to 99999
        logical :: OP           ! check if the unit is used
        
        ErrorFlag=0 !
        ! find for the first free unit
        do free_unit=11,99999!
            inquire( UNIT=free_unit, OPENED=OP)!
            if( free_unit .eq. 99999 ) ErrorFlag=1  ! no unit is available
            if(.not. OP) exit                       ! the unit is not used
        end do
    end subroutine seek_un
    
    pure function get_value_index_i( val_list, a_value) result( val_idx)
        ! Find the position of an integer value in an array of integers
        integer, dimension(:), intent(in) :: val_list    ! a list of values
        integer, intent(in) :: a_value                   ! the value to be found
        integer :: val_idx                               ! the position in the array
        integer :: i                                     ! counter
        !
        val_idx = 0 ! default value, value not found
        do i=1, size(val_list)
            if(val_list(i) == a_value) then
                val_idx = i
                return
            end if
        end do
    end function get_value_index_i
    
    pure function get_value_index_c( val_list, a_value) result( val_idx)!
        ! Find the position of a character value in an array of characters
        character(len=*), dimension(:), intent(in) :: val_list  ! a list of values
        character(len=*), intent(in) :: a_value                 ! the value to be found        
        integer :: val_idx                                      ! the position in the array
        integer :: i                                            ! counter
        !
        val_idx = 0 ! Default value, value not found
        do i=1, size(val_list)
            if(val_list(i) == a_value) then
                val_idx = i
                return
            end if
        end do!
    end function get_value_index_c

    subroutine calc_time_diff(t_start, t_stop, t_delta)
        ! Calculate the difference between two times
        ! TODO: difference in days not implemented
        !
        integer, dimension(8), intent(in) :: t_start, t_stop!
        integer, dimension(8), intent(out) :: t_delta
        integer :: milliseconds, seconds, minutes, hours, days!
        
        t_delta = 0 ! init to zero

        milliseconds = t_stop(8) - t_start(8)!
        seconds = t_stop(7) - t_start(7)!
        minutes = t_stop(6) - t_start(6)!
        hours = t_stop(5) - t_start(5)!
        days = t_stop(3) - t_start(3)!
        if( milliseconds < 0) then!
            milliseconds = milliseconds + 1000!
            seconds = seconds - 1!
        end if!
        if( seconds < 0 ) then!
            seconds = seconds + 60!
            minutes = minutes - 1!
        end if!
        if( minutes < 0 ) then!
            minutes = minutes + 60!
            hours = hours - 1!
        end if!
        if( hours < 0  ) then!
            hours = hours + 24!
            days = days - 1!
            !per ora nulla!
        end if!
        t_delta(3) = days
        t_delta(5) = days*24+hours ! add hours of completed days
        t_delta(6) = minutes
        t_delta(7) = seconds
        t_delta(8) = milliseconds
    end subroutine calc_time_diff
   
    subroutine print_execution_time(t_start, t_stop)!
        ! print the difference between time
        ! TODO: mode to cli
        integer, dimension(8), intent(in) :: t_start, t_stop
        integer, dimension(8) :: t_delta
        
        call calc_time_diff(t_start, t_stop, t_delta)
        print *, " ===> Simulation duration: ", t_delta(5), "h ", t_delta(6), "' ", t_delta(7), ' " '!
    end subroutine print_execution_time!
    
    pure function make_numbered_name(n,ext) result(num_name)!
        ! Dalla stringa nome del file (espresso in numero), estensione, genera un nome di file completo
        integer, intent(in)::n              ! number of file
        character(len=4),intent(in)::ext    ! file extention
        character(len=30)::num_name         ! complete numbered name
        write(num_name,*)n!
        num_name=trim(adjustl(num_name))//trim(adjustl(ext))!
    end function make_numbered_name!
    
    subroutine get_uniform_sample(irandom, amplitude, rand_symmetry,repeatable)!
        ! Get a matrix of random number from a matrix of integer values
        ! values are included between [-amplitude, +amplitude]
        implicit none!
        integer,dimension(:,:),intent(out)::irandom!
        integer,intent(in)::amplitude!
        logical,intent(in)::rand_symmetry
        logical, intent(in) :: repeatable
        real(dp),dimension(size(irandom,1),size(irandom,2))::rrandom!

        ! note about "random_init(repeatable ,image_distinct)"
        ! repeatable : is true, use the same initialization values
        ! image_distinct - mostly for coarray parallel programs, if true, each image has its random setup
        call random_init(repeatable, .false.)
        call random_number(rrandom)             ! generate a range [0-1]!

        if (rand_symmetry .eqv. .true.) then
            ! TODO: check
            irandom=int(amplitude*(2*rrandom-1)) ! transform R[0,1] in iR [-n,+n] with n = amplitude (e.g. amplitude = 4, iR = 0,-9)
        else
            irandom=int(rrandom*(2.d0*amplitude+1)) ! transform R[0,1] in iR [0,2n+1] with n = amplitude (e.g. amplitude = 4, iR = 0,-9)
        end if
        
    end subroutine get_uniform_sample

    function calc_doy(idd,imm,iyyy)
        ! Calculate the day of the year from day, month, year
        ! source: Numerical recipes in FORTRAN 90
        implicit none
        integer, intent(in) :: imm, idd, iyyy
        integer :: calc_doy
        integer, parameter :: igreg=15+31*(10+12*1582) ! Gregorian Calendar adopted Oct. 15, 1582.
        integer :: ja,jm,jy
        jy=iyyy
        if (jy == 0) stop 'calc_doy: there is no year zero'
        if (jy < 0) jy = jy+1
        if (imm > 2) then
            jm = imm+1
        else
            jy = jy-1
            jm = imm+13
        end if
        calc_doy = floor(365.25 * jy) + floor(30.6001 * jm) + idd + 1720995
        if (idd + 31 * (imm+12*iyyy) >= igreg) then !Test whether to change to Gregorian Calendar.
            ja=floor(0.01 * jy)
            calc_doy = calc_doy + 2 - ja + floor(0.25 * ja)
        end if
    end function calc_doy
    !
    subroutine calc_date(julian_day,idd, imm,iyyy)
        ! Calculate the date from the julian date
        ! source: Numerical recipes in FORTRAN 90
        implicit none
        integer, intent(in) :: julian_day
        integer, intent(out) :: imm, idd, iyyy
        INTEGER :: ja,jalpha,jb,jc,jd,je
        INTEGER, PARAMETER :: igreg=2299161
        if (julian_day >= igreg) then
            jalpha=int(((julian_day-1867216)-0.25)/36524.25)
            ja=julian_day+1+jalpha-int(0.25*jalpha)
        else if (julian_day < 0) then
            ja=julian_day+36525*(1-julian_day/36525)
        else
            ja=julian_day
        end if
        jb=ja+1524
        jc=int(6680.0+((jb-2439870)-122.1)/365.25)
        jd=365*jc+int(0.25*jc)
        je=int((jb-jd)/30.6001)
        idd=jb-jd-int(30.6001*je)
        imm=je-1
        if (imm > 12) imm = imm-12
        iyyy=jc-4715
        if (imm > 2) iyyy=iyyy-1
        if (iyyy <= 0) iyyy=iyyy-1
        if (julian_day < 0) iyyy=iyyy-100*(1-julian_day/36525)
    end subroutine calc_date
    !
    function day_of_week(idd, imm, iyy)
        ! Calculate the day of the week [0- Sat, 1-Sun, ..., 6-Fri]
        ! source: Rosetta Code
        implicit none
        integer, intent(in) :: idd, imm, iyy
        integer :: day_of_week, j, k, mm, yy

        if (imm < 2) then
            mm = imm + 12
            yy = iyy - 1
        end if
        j = yy/100           ! first two digits of the year
        k = mod (yy, 100)    ! last two digits of the year
        day_of_week = mod (idd + (mm+1)*26/10 + k + k/4 + j/4 + 5*j, 7)
    end function day_of_week
    
    subroutine days_x_month(calendar,year)!
        ! Get the number of days for each months in the provided year
        integer,dimension(:),intent(out)::calendar!
        integer,intent(in)::year!
        !
        if(mod(year,400)==0 .or. (mod(year,4)==0 .and. (.not.(mod(year,100)==0)))) then  ! leap year
            calendar=(/31,29,31,30,31,30,31,31,30,31,30,31/)!
        else    ! other
            calendar=(/31,28,31,30,31,30,31,31,30,31,30,31/)  !
        end if!
    end subroutine days_x_month!
    
    subroutine split_date_c(instring, date1, date2)
        ! return the string dates from a string where are separated by delimeter
        character(len=*), intent(in) :: instring
        character(len=300), intent(out):: date1, date2
        character(len=300) :: string
        character(Len=2), parameter:: delimiter = '->'
        integer :: index

        string = TRIM(instring)
        index = SCAN(string, delimiter)
        if (index == 0) then
            print *, 'Input files do not list end date or the delimiter "->" is not used'
            print *, 'Execution will be aborted...'
            stop
        end if
        date1 = string(1:index-1)
        index = SCAN(string, delimiter, .TRUE.)
        date2 = string(index+1:)
    END SUBROUTINE split_date_c
    !
    subroutine split_date_d(instring, outdate)
        ! return the date from a string where are separated by delimiter in the order day/month/year
        ! TODO: general format
        character(len=*), intent(in) :: instring
        type(date), intent(out):: outdate
        character(Len=1), parameter:: delimiter = '/'
        character(len=300) :: string
        character(len=4) :: date_num
        integer :: index
        integer :: ErrorFlag

        ErrorFlag = 0
        string = trim(instring)
        index = scan(string, delimiter)
        if (index == 0) ErrorFlag = -1
        date_num = adjustl(string(1:index-1)) ! select the day
        read(date_num, '(i2)') outdate%day
        string = string(index+1:)
        index = scan(string, delimiter)
        if (index == 0) ErrorFlag = -1
        date_num = string(1:index-1) ! select the month
        read(date_num, '(i2)') outdate%month
        date_num = string(index+1:) ! select the year
        read(date_num, '(i4)') outdate%year
        ! TODO: add control to check date validity
        if (ErrorFlag == -1) then
            print *, 'Input files are not correctly formatted'
            print *, 'Right date format is dd/mm/yyyy'
            print *, 'Execution will be aborted...'
            stop
        end if
    end subroutine split_date_d
    
    function string_to_integers(str, sep) result(a)
        ! return a sequence of integers from a string
        integer, allocatable :: a(:)
        character(*) :: str
        character :: sep
        integer :: i, n_sep

        n_sep = 0
        
        do i = 1, len(str)
          if (str(i:i)==sep) then
            n_sep = n_sep + 1
            str(i:i) = ','
           end if
        end do
        allocate(a(n_sep+1))
        read(str,*) a
    end function string_to_integers
    
    function string_to_reals(str, sep) result(a)
        ! return a sequence of reals from a string
        real(dp), allocatable :: a(:)
        character(*) :: str
        character :: sep
        integer :: i, n_sep

        n_sep = 0
        do i = 1, len(str)
          if (str(i:i)==sep) then
            n_sep = n_sep + 1
            str(i:i) = ','
           end if
        end do
        
        allocate(a(n_sep+1))
        read(str,*) a
    end function string_to_reals

    pure recursive function replace_str(string,search,substitute) result(modifiedString)
        ! https://stackoverflow.com/questions/58938347/how-do-i-replace-a-character-in-the-string-with-another-charater-in-fortran
        implicit none
        character(len=*), intent(in)  :: string, search, substitute
        character(len=:), allocatable :: modifiedString
        integer                       :: i, stringLen, searchLen
        stringLen = len(string)
        searchLen = len(search)
        if (stringLen==0 .or. searchLen==0) then
            modifiedString = ""
            return
        elseif (stringLen<searchLen) then
            modifiedString = string
            return
        end if
        i = 1
        do
            if (string(i:i+searchLen-1)==search) then
                modifiedString = string(1:i-1) // substitute // replace_str(string(i+searchLen:stringLen),search,substitute)
                exit
            end if
            if (i+searchLen>stringLen) then
                modifiedString = string
                exit
            end if
            i = i + 1
            cycle
        end do
    end function replace_str

    function round(val, n)
        implicit none
        real(dp) :: val, round
        integer :: n
        round = anint(val*10.0**n)/10.0**n
    end function round

    pure elemental function pdf_normal(x, x_mean, x_std) result(pdf)
        implicit none
        real(dp), intent(in):: x, x_mean, x_std
        real(dp):: pdf
        
        pdf = (1/((2*pi*x_std**2)**0.5))*exp(-((x-x_mean)**2)/(2*x_std**2))
    end function
!
end module mod_utility!
