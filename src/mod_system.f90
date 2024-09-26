module mod_system
    implicit none
! don't indent macro call
#if WIN == 1
    ! default is unix
    character(len = 10), private :: mkdir_cmd = 'mkdir'
    character :: delimiter = '\\'
#elif
    ! default is unix
    character(len = 10), private :: mkdir_cmd = 'mkdir'
    character :: delimiter = '/'
#endif

    
    contains

    subroutine make_dir(path)
        character(*), intent(in) :: path
        !print*,trim(mkdir_cmd)//' .'//delimiter//trim(path)
        call system(trim(mkdir_cmd)//' .'//delimiter//trim(path))
    end subroutine

    subroutine print_os_settings()
        print*, "** OS Settings **"
        print*,'make dir command: ',mkdir_cmd
        print*,'folder separator: ', delimiter
        print*, '================='
    end subroutine
    
end module 