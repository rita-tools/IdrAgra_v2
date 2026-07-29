module mod_system
implicit none

character(len=10), private :: mkdir_cmd = 'mkdir'

! don't indent macro call
#if WIN == 1
!%PS%: used to set delimiter = '\\', but it was being truncated to just '\' because unspecified character length is 1. achar(92) is equivalent to, but safer than, '\'.
#define PATH_DELIMITER achar(92)

#else
! default is unix
#define PATH_DELIMITER '/'
#endif

character(len=1) :: delimiter = PATH_DELIMITER ! %PS%: Variable is assigned only once so that this doesn't look like an error to VS Code/Modern Fortran

#undef PATH_DELIMITER

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