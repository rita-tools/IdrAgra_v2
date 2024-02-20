module mod_TDx_index
    ! TODO: split UI functions from model
    use mod_constants, only: sp, dp, nan_r
    use mod_utility, only: seek_un
    use mod_grid, only: grid_i, grid_r, print_mat_as_grid

    ! type for DTx calculation
    type clock!
        integer::n_ind!
        integer,dimension(:),pointer::x!
        integer::td!
        integer::delay           ! time shift of report generation
    end type clock!

    type TDx_index!
        character(len=255)::mode!
        type(clock)::temp!
        integer::n!
    end type TDx_index!

    
    contains

    subroutine calc_TDx(domain,gg,year,DxiTOT,kcb,gg_max,max_year,x,TD,unit_Dxi)
        implicit none
        type(grid_i),intent(in)::domain
        integer,intent(in)::gg,year
        real(dp),dimension(:,:),intent(in)::kcb
        integer,intent(in)::gg_max,max_year
        real(dp),dimension(size(domain%mat,1),size(domain%mat,2)),intent(out)::DxiTOT
        integer,intent(in)::x
        integer,dimension(:),intent(inout)::unit_Dxi
        real(dp),dimension(:,:),intent(in)::TD              ! transpiration deficit
        real(dp),dimension(size(domain%mat,1),size(domain%mat,2))::Dxi_tempN,Dxi_temp
        integer::index_Dxi
        integer::i,errorflag,ios
        real(dp),parameter::nan=-9999.

        Dxi_temp = 0.
        if(gg==1 .and. year==1)then
            do i=1,size(unit_Dxi)
                call seek_un(errorflag,unit_Dxi(i))
                open ( unit = unit_Dxi(i), status = 'scratch', form = 'unformatted', &
                    & action = 'readwrite', iostat = ios)
                if(ios/=0)then
                    print *, "Error opening file scratch connected to unit unit_Dxi(",i,")=", unit_Dxi(i), &
                        & " iostat=", ios, ". Execution will be aborted..."
                    stop
                end if
                write(unit_Dxi(i),iostat=ios)Dxi_temp

                if(ios/=0)then
                    print *, "Error writing unit_Dxi(",i,"). Execution will be aborted..."
                    stop
                end if
                rewind(unit_Dxi(i),iostat=ios)
                if(ios/=0)then
                    print *, "Error rewinding unit_Dxi(",i,"). Execution will be aborted..."
                    stop
                end if
            end do
        end if
        !selection of unit_Dxi layer
        index_Dxi =mod(gg,x)
        if(index_Dxi==0) index_Dxi=x
        !updating of Dxi values in scratch file
        read(unit_Dxi(index_Dxi),iostat=ios)Dxi_tempN
        if(ios/=0)then
            print *, "Error reading unit_Dxi(",index_Dxi,"). . Execution will be aborted..."
            stop
        end if
        rewind(unit_Dxi(index_Dxi),iostat=ios)
        if(ios/=0)then
            print *, "Error rewinding unit_Dxi(",index_Dxi,") after reading the file. Execution will be aborted..."
            stop
        end if
        if(gg>=x)then
            where(kcb/=0.)                  ! from emergence to harvesting
                where(Dxi_tempN/=nan)
                    DxiTOT = (TD-Dxi_tempN)
                else where
                    DxiTOT = nan
                end where
            else where
                DxiTOT = nan
            end where
        else
            DxiTOT = nan
        end if

        !print Dxi values in scratch file
        write(unit_Dxi(index_Dxi),iostat=ios)TD
        if(ios/=0)then
            print *, "Error writing unit_Dxi(",index_Dxi,"). Execution will be aborted..."
            stop
        end if
        rewind(unit_Dxi(index_Dxi),iostat=ios)
        if(ios/=0)then
            print *, "Error rewinding unit_Dxi(",index_Dxi,") after writing the file. Execution will be aborted..."
            stop
        end if
        !deallocates Dxi files on last day of simulation
        if(year==max_year .and. gg==gg_max)then
            do i=1,size(unit_Dxi)
                close(unit_Dxi(i),iostat=ios)
                if(ios/=0)then
                    print*,"Error closing unit_Dxi(",i,"). Execution will be aborted..."
                    stop
                end if
            end do
        end if
    end subroutine calc_TDx

    subroutine init_TDx(unit_deficit,domain,path,TDx)
    ! initializes scratch files to NaN
        implicit none
        integer,dimension(:),intent(in)::unit_deficit
        type(grid_i),intent(in)::domain
        character(len=*),intent(in)::path
        real(dp),dimension(size(domain%mat,1),size(domain%mat,2))::temp
        integer,dimension(size(domain%mat,1),size(domain%mat,2))::num
        integer::i,j,errorflag,ios
        character(len=*),intent(in)::TDx
        integer::free_unit
        character(len=55)::str_y, str7, nfile

        temp=real(domain%header%nan)
        num=0
        do j=1,unit_deficit(2)
            write(str_y,*)j
            do i=1,unit_deficit(1)
                write(str7,*)i
                call seek_un(errorflag,free_unit)
                nfile=TDx//'_'//trim(adjustl(str7))//'_'//trim(adjustl(str_y))//'.tmp'
                open(unit=free_unit,file=trim(path)//trim(nfile),form='unformatted',action='readwrite',iostat=ios)
                if(ios/=0)then
                    print *, "Error opening file tmp connected to unit =", free_unit, " iostat=", ios, &
                        & " nfile=",trim(nfile), ". Execution will be aborted..."
                    stop
                end if
                if(j==1 .or. j==2)then
                    write(free_unit,iostat=ios)temp
                else
                    write(free_unit,iostat=ios)num
                end if
                if(ios/=0)then
                    print *, "Error writing file ",trim(nfile)," in subroutine iniz_TDx; iostat=", ios, &
                        & ". Execution will be aborted..."
                    stop
                end if
                close(free_unit,iostat=ios)
                if(ios/=0)then
                    print*,"Error closing file ",trim(nfile), "iostat=", ios, ". Execution will be aborted..."
                    stop
                end if
            end do
        end do
    end subroutine init_TDx

    subroutine update_TDx_DB(deficit,TDx,time_step,domain,path)
    ! database updating
        implicit none
        real(dp),dimension(:,:),intent(in)::deficit
        integer,intent(in)::time_step
        type(grid_i),intent(in)::domain
        real(dp),dimension(size(domain%mat,1),size(domain%mat,2),2)::temp_vect
        real(dp),dimension(size(domain%mat,1),size(domain%mat,2))::temp
        real(dp),dimension(size(domain%mat,1),size(domain%mat,2))::TD_sum
        real(dp),dimension(size(domain%mat,1),size(domain%mat,2))::deficit_nozeros
        integer,dimension(size(domain%mat,1),size(domain%mat,2))::count_num         ! values count
        integer,dimension(size(domain%mat,1),size(domain%mat,2))::count_zeros       ! zeros count
        character(len=*),intent(in)::TDx
        character(len=*),intent(in)::path
        integer::ios
        integer::free_unit,errorflag
        character(len=55)::str_y,str7,nfile
        real(dp),parameter::nan=-9999.

        write(str7,*)time_step
        ! sum(x) over cycle
        temp_vect(:,:,1)=deficit(:,:)       ! populates x,j,1 of temp_vect with tot_deficit ( tot_deficit(:,:,cont_dt))
        write(str_y,*)1
        nfile=TDx//'_'//trim(adjustl(str7))//'_'//trim(adjustl(str_y))//'.tmp'
        call seek_un(free_unit,errorflag)
        open(unit=free_unit,file=trim(path)//trim(nfile),form='unformatted',action='readwrite',iostat=ios, status='old')
        if(ios/=0)then
            print *, "Error opening file tmp connected to unit =", free_unit, " iostat=", ios, &
                & " nfile=",trim(nfile), ". Execution will be aborted..."
            stop
        end if
        read(free_unit,iostat=ios)temp
        if(ios/=0)then
            print *, "Error reading in file ",trim(nfile)," in subroutine TDx_DB; iostat=",  ios, &
                & ". Execution will be aborted..."
            stop
        end if
        rewind(free_unit,iostat=ios)
        if(ios/=0)then
            print *, "Error rewinding file ",trim(nfile)," in subroutine TDx_DB; iostat=", ios, &
                & ". Execution will be aborted..."
            stop
        end if
        temp_vect(:,:,2)=temp(:,:)       ! populates x,j,2 of temp_vect with TDx_"str7"_1.tmp

        TD_sum=nan
        TD_sum(:,:)=sum(temp_vect(:,:,:),dim=3,mask=(temp_vect(:,:,:)/=nan)) ! sums tot_deficit to TDx_"str7"_1.tmp and populates TD_sum
        ! at cycle end, TD_sum contains TDX sum for each selected period over different years

        write(free_unit,iostat=ios)TD_sum                                ! the result is copied to TDx_"str7"_1.tmp
        if(ios/=0)then
            print *, "Error writing file ",trim(nfile)," in subroutine TDx_DB; iostat=", ios, &
                & ". Execution will be aborted..."
            stop
        end if
        close(free_unit)

        ! sum(ln(x))
        deficit_nozeros=merge(deficit,nan,deficit>0)                ! zeros are not taken into account, as ln(0) is not defined
        temp_vect(:,:,1)=merge(log(deficit_nozeros),nan,deficit_nozeros/=nan)
        write(str_y,*)2
        nfile=TDx//'_'//trim(adjustl(str7))//'_'//trim(adjustl(str_y))//'.tmp'
        call seek_un(free_unit,errorflag)
        open(unit=free_unit,file=trim(path)//trim(nfile),form='unformatted',action='readwrite',iostat=ios, status='old')
        if(ios/=0)then
            print *, "Error opening file tmp connected to unit =", free_unit, " iostat=", ios, " nfile=",trim(nfile), &
                & ". Execution will be aborted..."
            stop
        end if
        read(free_unit,iostat=ios)temp
        if(ios/=0)then
            print *, "Error reading file ",trim(nfile)," in subroutine TDx_DB; iostat=", ios, &
                & ". Execution will be aborted..."
            stop
        end if
        rewind(free_unit,iostat=ios)
        if(ios/=0)then
            print *, "Error rewinding file ",trim(nfile)," in subroutine TDx_DB; iostat=", ios, &
                & ". Execution will be aborted..."
            stop
        end if
        temp_vect(:,:,2)=temp(:,:)
        TD_sum=nan
        TD_sum(:,:)=sum(temp_vect(:,:,:),dim=3,mask=(temp_vect(:,:,:)/=nan)) ! sums ln(tot_deficit) to TDx_"str7"_2.tmp and populates TD_sum
        ! at cycle end, TD_sum contains ln(TDX) sum for each selected period over different years

        write(free_unit,iostat=ios)TD_sum                                ! the result is copied to TDx_"str7"_2.tmp
        if(ios/=0)then
            print *, "Error writing file ",trim(nfile),"in subroutine TDx_DB; iostat=", ios, &
                & ". Execution will be aborted..."
            stop
        end if
        close(free_unit)

        ! count values
        count_num=0
        write(str_y,*)3
        nfile=TDx//'_'//trim(adjustl(str7))//'_'//trim(adjustl(str_y))//'.tmp'
        call seek_un(free_unit,errorflag)
        open(unit=free_unit,file=trim(path)//trim(nfile),form='unformatted',action='readwrite',iostat=ios, status='old')
        if(ios/=0)then
            print *, "Error opening file tmp connected to unit =", free_unit, " iostat=", ios, " nfile=",trim(nfile), &
                & ". Execution will be aborted..."
            stop
        end if
        read(free_unit,iostat=ios)count_num
        if(ios/=0)then
            print *, "Error reading file ",trim(nfile)," in subroutine TDx_DB; iostat=", ios, &
                & ". Execution will be aborted..."
            stop
        end if
        rewind(free_unit,iostat=ios)
        if(ios/=0)then
            print *, "Error rewinding file ",trim(nfile)," in subroutine TDx_DB; iostat=", ios, &
                & ". Execution will be aborted..."
            stop
        end if
        
        ! it takes into account values (different from NaN) of tot_deficit in TDx_"str7"_3.tmp
        ! at cycle end, count_num contains the count of values over which statistics are performed
        where(deficit/=nan)count_num=count_num+1 
        

        write(free_unit,iostat=ios)count_num
        if(ios/=0)then
            print *, "Error writing file ",trim(nfile)," in subroutine TDx_DB; iostat=", ios, &
                & ". Execution will be aborted..."
            stop
        end if
        close(free_unit)

        ! count zeros
        count_zeros=0
        write(str_y,*)4
        nfile=TDx//'_'//trim(adjustl(str7))//'_'//trim(adjustl(str_y))//'.tmp'
        call seek_un(free_unit,errorflag)
        open(unit=free_unit,file=trim(path)//trim(nfile),form='unformatted',action='readwrite',iostat=ios, status='old')
        if(ios/=0)then
            print *, "Error opening file tmp connected to unit =", free_unit, " iostat=", ios, " nfile=",trim(nfile), &
                & ". Execution will be aborted..."
            stop
        end if
        read(free_unit,iostat=ios)count_zeros
        if(ios/=0)then
            print *, "Error reading file ",trim(nfile)," in subroutine TDx_DB; iostat=", ios, &
                & ". Execution will be aborted..."
            stop
        end if
        rewind(free_unit,iostat=ios)
        if(ios/=0)then
            print *, "Error rewinding file ",trim(nfile)," in subroutine TDx_DB; iostat=", ios, &
                & ". Execution will be aborted..."
            stop
        end if

        ! it takes into account zeros of tot_deficit in TDx_"str7"_4.tmp
        ! at cycle end, count_zeros contains the count of zeros over which statistics are performed
        where(deficit==0.)count_zeros=count_zeros+1                         

        write(free_unit,iostat=ios)count_zeros
        if(ios/=0)then
            print *, "Error writing file ",trim(nfile)," xnum_zeros in subroutine TDx_DB; iostat=", ios, &
                & ". Execution will be aborted..."
            stop
        end if
        close(free_unit)

    end subroutine update_TDx_DB

    subroutine save_TDx_statistics(unit_deficit,domain,path,threshold_num,TDx)
        implicit none
        integer,dimension(:),intent(in)::unit_deficit
        type(grid_i),intent(in)::domain
        character(len=*),intent(in)::path
        integer,intent(in)::threshold_num
        real(dp),dimension(size(domain%mat,1),size(domain%mat,2))::TDx_mean,A_gamma,alpha_hat,beta_hat
        integer,dimension(size(domain%mat,1),size(domain%mat,2))::den
        real(dp),dimension(size(domain%mat,1),size(domain%mat,2),2)::TDx_ln_values
        integer,dimension(size(domain%mat,1),size(domain%mat,2),2)::TDx_ln_num
        real(dp),dimension(size(domain%mat,1),size(domain%mat,2))::temp

        integer,dimension(size(domain%mat,1),size(domain%mat,2))::num_nozeros
        real(dp),dimension(size(domain%mat,1),size(domain%mat,2))::zero_prob

        integer::time_step,y,ios
        character(len=255)::filename
        integer::free_unit,errorflag
        real(dp),parameter::nan=-9999.
        character(len=*),intent(in)::TDx
        character(len=55)::str_y,str7,nfile

        time_step_cycle: do time_step=1,unit_deficit(1)

            ! variables initialization
            num_nozeros = int(nan)
            zero_prob = nan
            TDx_mean = nan
            A_gamma = nan
            alpha_hat = nan
            beta_hat = nan

            write(str7,*)time_step
            do y=1,unit_deficit(2)
                write(str_y,*)y
                nfile=TDx//'_'//trim(adjustl(str7))//'_'//trim(adjustl(str_y))//'.tmp'
                call seek_un(free_unit,errorflag)
                open(unit=free_unit,file=trim(path)//trim(nfile),form='unformatted',action='readwrite',iostat=ios, status='old')
                if(ios/=0)then
                    print *, "Error opening file tmp connected to unit =", free_unit, " iostat=", ios, " nfile=",trim(nfile), &
                        & ". Execution will be aborted..."
                    stop
                end if
                if(y==1 .or. y==2)then      ! TD_sum_x and TD_sum_lnx
                    read(free_unit,iostat=ios) temp
                    TDx_ln_values(:,:,y)=temp(:,:)
                else                        ! values and zeros count
                    read(free_unit,iostat=ios) den
                    TDx_ln_num(:,:,y-2)=den(:,:)
                end if
                if(ios/=0)then
                    print *, "Error reading file ",trim(nfile)," in subroutine TDx_statistics; iostat=", ios, &
                        & ". Execution will be aborted..."
                    stop
                end if
                close(free_unit,iostat=ios)
                if(ios/=0)then
                    print*,"Error closing file ",trim(nfile), "iostat=", ios, ". Execution will be aborted..."
                    stop
                end if
            end do

            ! TDx_ln_values contains statistics (sum and logsum)
            ! TDx_ln_num contains days count

            ! evaluation of statistical significance
            where(TDx_ln_num(:,:,1)<threshold_num)TDx_ln_num(:,:,1)=int(nan)    ! comparison between values count and threshold
            where(TDx_ln_num(:,:,1)-TDx_ln_num(:,:,2)<threshold_num) &
                & TDx_ln_num(:,:,1)=int(nan)                                    ! comparison between values (different from zero) count and threshold
            where((TDx_ln_num(:,:,1)/=0) .and. (TDx_ln_num(:,:,1)/=int(nan)))   ! if valid values > threshold, parameters are calculated
                ! nozeros count
                num_nozeros (:,:) = TDx_ln_num(:,:,1)-TDx_ln_num(:,:,2)
                ! probability q of 0 value
                zero_prob (:,:) = float(TDx_ln_num(:,:,2)) / float(TDx_ln_num(:,:,1))   ! float converts integer to real
                ! mean
                TDx_mean(:,:)=TDx_ln_values(:,:,1)/num_nozeros(:,:)
                ! alfa of gamma function
                A_gamma (:,:) = log(TDx_mean(:,:))-TDx_ln_values(:,:,2)/num_nozeros(:,:)
                ! alpha_hat of gamma function
                alpha_hat(:,:) = 1 / (4.*A_gamma(:,:)*(1+sqrt(1.+4./3.*A_gamma(:,:))))
                ! beta_hat of gamma function
                beta_hat(:,:) = TDx_mean(:,:)/alpha_hat(:,:)
            end where

            ! output creation
            filename=trim(adjustl(trim(path)//TDx//'_alpha_'//trim(adjustl(str7))//'.asc'))
            call print_mat_as_grid(filename,domain%header,alpha_hat,errorflag)
            filename=trim(adjustl(trim(path)//TDx//'_beta_'//trim(adjustl(str7))//'.asc'))
            call print_mat_as_grid(filename,domain%header,beta_hat,errorflag)
            filename=trim(adjustl(trim(path)//TDx//'_zero_prob_'//trim(adjustl(str7))//'.asc'))
            call print_mat_as_grid(filename,domain%header,zero_prob,errorflag)
        end do time_step_cycle

    end subroutine save_TDx_statistics

    subroutine make_TDx_report(domain,path,n_week,TDx)
        implicit none
        type(grid_i),intent(in)::domain
        character(len=*),intent(in)::path
        integer,intent(in)::n_week
        integer::ios
        character(len=255)::filename
        integer::free_unit,errorflag
        real(dp),parameter::nan=-9999.
        real(dp),dimension(size(domain%mat,1),size(domain%mat,2))::mat
        character(len=*),intent(in)::TDx
        character(len=55)::str_y,str7,nfile

        ! TODO verify if different years can be written in the same run
        where (domain%mat/=domain%header%nan) mat =nan
        write(str7,*)n_week
        write(str_y,*)1
        nfile=TDx//'_'//trim(adjustl(str7))//'_'//trim(adjustl(str_y))//'.tmp'
        call seek_un(free_unit,errorflag)
        open(unit=free_unit,file=trim(path)//trim(nfile),form='unformatted',action='readwrite',iostat=ios, status='old')
        if(ios/=0)then
            print *, "Error opening file tmp connected to unit =", free_unit, " iostat=", ios, " nfile=",trim(nfile), &
            & ". Execution will be aborted..."
            stop
        end if
        read(free_unit,iostat=ios) mat(:,:)
        if(ios/=0)then
            print *, "Error reading file ",trim(nfile)," in subroutine TDx_report; iostat=", ios, &
                & ". Execution will be aborted..."
            stop
        end if
        close(free_unit,iostat=ios)
        if(ios/=0)then
            print*,"Error closing file ",trim(nfile), "iostat=", ios, ". Execution will be aborted..."
            stop
        end if
        where(domain%mat==domain%header%nan) mat=nan
        ! output creation
        filename=trim(adjustl(trim(path)//TDx//'_mm_'//trim(adjustl(str7))//'.asc'))
        call print_mat_as_grid(filename,domain%header,mat,errorflag)
    end subroutine make_TDx_report

    subroutine sum_TD(transp_act,transp_pot,k_cb,doy,year,TD)
        implicit none
        real(dp),dimension(:,:),intent(in)::transp_act,transp_pot
        integer,intent(in)::doy,year
        real(dp),dimension(:,:),intent(in)::k_cb
        real(dp),dimension(:,:),intent(inout)::TD
        real(dp),parameter::nan = -nan_r

        ! TD is NaN outside the growing period (k_cb /= 0)
        if(doy==1 .and. year==1) TD=0.
        where(k_cb/=0.)
            where(TD==nan)TD=0.
            TD = TD+(transp_pot-transp_act)
        else where
            TD=nan
        end where
        
    end subroutine sum_TD
end module mod_TDx_index
