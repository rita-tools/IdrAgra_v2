module mod_system
    implicit none

    ! default is unix
    character(len = 10), private :: mkdir_cmd = 'mkdir'
    character :: delimiter = '/'
    
    contains

    subroutine make_dir(path)
        character(200), intent(in) :: path
! don't indent macro!!!
#if WIN == 1
        mkdir_cmd = 'mkdir'
        delimiter = '\\'
#endif
        !print*,'mkdir: ',mkdir_cmd,' delimiter: ',delimiter,' path: ',trim(path)
        call system(trim(mkdir_cmd)//' .'//delimiter//trim(path))
    end subroutine
    
end module 