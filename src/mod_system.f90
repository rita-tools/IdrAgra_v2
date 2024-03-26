module mod_system
    implicit none

    character(len = 10) :: mkdir_cmd_priv
    character :: delimiter_priv

    contains

    subroutine set_command(mkdir_os,delimiter_os)
        character(len = *), intent(in) :: mkdir_os
        character, intent(in) :: delimiter_os
        mkdir_cmd_priv = mkdir_os
        delimiter_priv = delimiter_os
    
    end subroutine

    subroutine make_dir(path)
        character(200), intent(in) :: path
        call system(trim(mkdir_cmd_priv)//' '//delimiter_priv//trim(path))
    end subroutine
    
end module 