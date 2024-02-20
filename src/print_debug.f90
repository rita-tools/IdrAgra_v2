subroutine print_hello(val)
    !character(len=20), intent(in) :: name
    integer, intent(in) :: val
    print*,'hello',val
end subroutine print_hello

subroutine print_2darray_real(mat)
    real(8),dimension(:,:),intent(in)::mat
    integer i,j
    do j=1,size(mat,2)!
        do i=1,size(mat,1)!
            print*,mat(i,j)
        end do!
    end do!
end subroutine

subroutine print_2darray_int(mat)
    integer,dimension(:,:),intent(in)::mat
    integer i,j
    do j=1,size(mat,2)!
        do i=1,size(mat,1)!
            print*,mat(i,j)
        end do!
    end do!
end subroutine


subroutine test_exp()
    integer i
    do i =1, 1000
        print*,exp(-1.*i)
    end do!
end subroutine
