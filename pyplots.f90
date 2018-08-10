Module pyplots
  use Kinds
  implicit none

  INTERFACE plot  
     MODULE PROCEDURE plot,plot_notitle,stan_plot,stan_plot_notitle,plot_legend,stan_plot_two,stan_plot_two_legend
  END INTERFACE plot

contains
    
  !---------------------------------------------------------------------!
  !                                                                     !
  !                       Multi-Dimensional Plots                       !
  !                                                                     !
  !---------------------------------------------------------------------!

  subroutine plot(x,xaxis,yaxis)
    use Kinds
    implicit none
    real(DP),dimension(:,:),intent(in) :: x
    character(len=*),intent(in)        :: xaxis,yaxis
    character(len=9)                   :: fmt
    integer                            :: ii,j,N

    N = size(x,dim=2)
    
    open(100,file="titles.dat",action="write", &
         & status="replace",form="formatted")        
    open(101,file="output.dat",action="write", &
         & status="replace",form="formatted")

    write(100,'(a15)') xaxis
    write(100,'(a15)') yaxis
    write(100,'(a3)') "nil"
    write(100,'(a3)') "nil"
    write(fmt,'(a1,i1,a7)') '(', N, 'es20.9)'

    do ii=1,size(x,dim=1)
       write(101,fmt) x(ii,:)
    end do

    close(100)
    close(101)

    call system("./pyplots.py")

  end subroutine plot

  subroutine plot_legend(x,xaxis,yaxis,legend)
    use Kinds
    implicit none
    character(len=*),intent(in)              :: xaxis,yaxis
    character(len=*),dimension(:),intent(in) :: legend
    real(DP),dimension(:,:),intent(in)       :: x
    character(len=9)                         :: fmt
    integer                                  :: ii,j,N

    N = size(x,dim=2)
    
    open(100,file="titles.dat",action="write", &
         & status="replace",form="formatted")        
    open(101,file="output.dat",action="write", &
         & status="replace",form="formatted")

    write(100,'(a15)') xaxis
    write(100,'(a15)') yaxis
    write(fmt,'(a1,i1,a7)') '(', N, 'es20.9)'

    do ii=1,size(legend)
       write(100,'(a12)') legend(ii)
    end do

    do ii=1,size(x,dim=1)
       write(101,fmt) x(ii,:)
    end do

    close(100)
    close(101)

    call system("./pyplots.py")

  end subroutine plot_legend

  subroutine plot_notitle(x)
    use Kinds
    implicit none
    real(DP),dimension(:,:),intent(in) :: x
    character(len=9)                   :: fmt
    integer                            :: ii,j,N

    N = size(x,dim=2)
    
    open(100,file="titles.dat",action="write", &
         & status="replace",form="formatted")        
    open(101,file="output.dat",action="write", &
         & status="replace",form="formatted")

    write(100,'(a15)') " "
    write(100,'(a15)') " "
    write(100,'(a3)') "nil"
    write(100,'(a3)') "nil"
    write(fmt,'(a1,i1,a7)') '(', N, 'es20.9)'

    do ii=1,size(x,dim=1)
       write(101,fmt) x(ii,:)
    end do

    close(100)
    close(101)

    call system("./pyplots.py")

  end subroutine plot_notitle

  !---------------------------------------------------------------------!
  !                                                                     !
  !                            Standard Plots                           !
  !                                                                     !
  !---------------------------------------------------------------------!

  subroutine stan_plot(x,y,xaxis,yaxis)
    use Kinds
    implicit none
    real(DP),dimension(:),intent(in)   :: x,y
    character(len=*),intent(in)        :: xaxis,yaxis
    integer                            :: ii
    
    open(100,file="titles.dat",action="write", &
         & status="replace",form="formatted")        
    open(101,file="output.dat",action="write", &
         & status="replace",form="formatted")

    write(100,'(a15)') xaxis
    write(100,'(a15)') yaxis
    write(100,'(a3)') "nil"
    write(100,'(a3)') "nil"
    do ii=1,size(x)
       write(101,'(2es20.9)') x(ii), y(ii)
    end do

    close(100)
    close(101)

    call system("./pyplots.py")

  end subroutine stan_plot

  subroutine stan_plot_notitle(x,y)
    use Kinds
    implicit none
    real(DP),dimension(:),intent(in)   :: x,y
    integer                            :: ii
    
    open(100,file="titles.dat",action="write", &
         & status="replace",form="formatted")        
    open(101,file="output.dat",action="write", &
         & status="replace",form="formatted")

    write(100,'(a15)') " "
    write(100,'(a15)') " "
    write(100,'(a3)') "nil"
    write(100,'(a3)') "nil"
    do ii=1,size(x)
       write(101,'(2es20.9)') x(ii), y(ii)
    end do

    close(100)
    close(101)

    call system("./pyplots.py")

  end subroutine stan_plot_notitle

  subroutine stan_plot_two(x,y,z,w)
    use Kinds
    implicit none
    real(DP),dimension(:),intent(in)   :: x,y,z,w
    integer                            :: ii
    
    open(100,file="titles.dat",action="write", &
         & status="replace",form="formatted")        
    open(101,file="output.dat",action="write", &
         & status="replace",form="formatted")

    write(100,'(a15)') " "
    write(100,'(a15)') " "
    write(100,'(a3)') "nil"
    write(100,'(a3)') "nil"
    do ii=1,size(x)
       write(101,'(4es20.9)') x(ii), y(ii), z(ii), w(ii)
    end do

    close(100)
    close(101)

    call system("./pyplots.py")

  end subroutine stan_plot_two

  subroutine stan_plot_two_legend(x,y,z,w,len1,len2)
    use Kinds
    implicit none
    real(DP),dimension(:),intent(in)   :: x,y,z,w
    character(len=*),intent(in)        :: len1,len2
    integer                            :: ii
    
    open(100,file="titles.dat",action="write", &
         & status="replace",form="formatted")        
    open(101,file="output.dat",action="write", &
         & status="replace",form="formatted")

    write(100,'(a15)') " "
    write(100,'(a15)') " "
    write(100,'(a9)') len1
    write(100,'(a9)') len2
    do ii=1,size(x)
       write(101,'(4es20.9)') x(ii), y(ii), z(ii), w(ii)
    end do

    close(100)
    close(101)

    call system("./pyplots.py")

  end subroutine stan_plot_two_legend

  
end Module pyplots
