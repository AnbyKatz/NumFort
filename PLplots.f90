module PLplots
  use kinds
  use plplot
  implicit none

  interface plot
     module procedure plot,plotmany
  end interface plot

contains

  !---------------------------------------------------------------------!
  !                                                                     !
  !                             PLplot Code                             !
  !                                                                     !
  !---------------------------------------------------------------------!

  ! See documentation for further description of the plotting algorithim

  subroutine plot(x,y,xlabel,ylabel,title)
    use kinds
    use plplot
    implicit none

    real(DP),dimension(:),intent(in) :: x,y
    character(len=*),optional,intent(in) :: xlabel,ylabel,title

    real(DP) :: ymin,ymax,xmin,xmax
    character(len=23) :: xaxis,yaxis,name

    xaxis = ""
    yaxis = ""
    name  = ""

    if(present(xlabel)) xaxis = xlabel
    if(present(ylabel)) yaxis = ylabel
    if(present(title)) name = title

    ymin = minval(y)
    ymin = ymin-0.2_DP*abs(ymin)
    ymax = maxval(y)
    ymax = ymax+0.2_DP*abs(ymax)
    xmin = minval(x)
    xmax = maxval(x)

    call plsdev("xwin")
    call plinit
    call plenv(xmin,xmax,ymin,ymax,0,0)
    call pllab(xaxis,yaxis,name)
    call plcol0(3)
    call plline(x,y)
    call plend

  end subroutine plot

  !---------------------------------------------------------------------!
  !                                                                     !
  !                             Scatter Plot                            !
  !                                                                     !
  !---------------------------------------------------------------------!

  subroutine scatterplot(x,y,style,xlabel,ylabel,title)
    use kinds
    use plplot
    implicit none

    real(DP),dimension(:),intent(in) :: x,y
    character(len=*),intent(in)      :: style
    character(len=*),optional,intent(in) :: xlabel,ylabel,title

    real(DP) :: ymin,ymax,xmin,xmax
    character(len=23) :: xaxis,yaxis,name

    xaxis = ""
    yaxis = ""
    name  = ""

    if(present(xlabel)) xaxis = xlabel
    if(present(ylabel)) yaxis = ylabel
    if(present(title)) name = title

    ymin = minval(y)
    ymin = ymin-0.2_DP*abs(ymin)
    ymax = maxval(y)
    ymax = ymax+0.2_DP*abs(ymax)
    xmin = minval(x)
    xmax = maxval(x)

    call plsdev("xwin")
    call plinit
    call plenv(xmin,xmax,ymin,ymax,0,0)
    call pllab(xaxis,yaxis,name)
    call plcol0(3)
    call plstring(x,y,style)
    call plend

  end subroutine scatterplot

  !---------------------------------------------------------------------!
  !                                                                     !
  !                              MultiPlot                              !
  !                                                                     !
  !---------------------------------------------------------------------!

  subroutine plotMany(data,xlabel,ylabel,title)
    use kinds
    use plplot
    implicit none

    real(DP),dimension(:,:),intent(in) :: data
    character(len=*),optional,intent(in) :: xlabel,ylabel,title


    real(DP) :: ymin,ymax,xmin,xmax,newval
    integer :: i,N
    character(len=23) :: xaxis,yaxis,name

    xaxis = ""
    yaxis = ""
    name  = ""

    if(present(xlabel)) xaxis = xlabel
    if(present(ylabel)) yaxis = ylabel
    if(present(title)) name = title

    N = size(data(1,:))
    xmin = data(1,1)
    xmax = data(1,1)
    ymin = data(1,2)
    ymax = data(1,2)

    do i = 1,N,2
       newval = minval(data(:,i+1))
       if ( newval < ymin ) ymin = newval
       newval = maxval(data(:,i+1))
       if ( newval > ymax ) ymax = newval
       newval = minval(data(:,i))
       if ( newval < xmin ) xmin = newval
       newval = maxval(data(:,i))
       if ( newval > xmax ) xmax = newval
    end do

    ymin = ymin-0.2_DP*abs(ymin)
    ymax = ymax+0.2_DP*abs(ymax)

    call plsdev("xwin")
    call plinit
    call plenv(xmin,xmax,ymin,ymax,0,0)
    call pllab(xaxis,yaxis,name)
    do i = 1,N/2
       call plcol0(i+2)
       call plline(data(:,2*i-1),data(:,2*i))
    end do
    call plend

  end subroutine plotmany

  !---------------------------------------------------------------------!
  !                                                                     !
  !                             Surface Plot                            !
  !                                                                     !
  !---------------------------------------------------------------------!

  subroutine surf(X,Y,Z,xlabel,ylabel,zlabel,title)
    use kinds
    use plplot
    implicit none

    real(DP),dimension(:),intent(in)   :: x,y
    real(DP),dimension(:,:),intent(in) :: Z
    character(len=*),optional,intent(in)        :: xlabel,ylabel,zlabel,title

    integer                :: opt=3,nxsub,nysub,nzsub
    real(DP)               :: xmin,xmax,ymin,ymax,zmin,&
         & zmax,xtick,ytick,ztick
    real(DP)               :: altitude    =  20.0_DP
    real(DP)               :: azimuth     = 298.0_DP
    real(DP)               :: azimuthStep =   0.2_DP
    real(DP)               :: basex  = 1.0_DP
    real(DP)               :: basey  = 1.0_DP
    real(DP)               :: height = 1.0_DP
    real(kind=pl_test_flt) :: alt = 33.0_pl_test_flt
    real(kind=pl_test_flt) :: az  = 24.0_pl_test_flt

    character(len=23) :: xaxis,yaxis,zaxis,name

    xaxis = ""
    yaxis = ""
    zaxis = ""
    name  = ""

    if(present(xlabel)) xaxis = xlabel
    if(present(ylabel)) yaxis = ylabel
    if(present(zlabel)) zaxis = zlabel
    if(present(title)) name = title

    xtick=0
    ytick=0
    ztick=0
    nxsub=0
    nysub=0
    nzsub=0

    xmin = minval(x)
    xmax = maxval(x)
    ymin = minval(y)
    ymax = maxval(y)
    zmin = minval(z)
    zmax = maxval(z)

    call plsdev("xwin")
    call plinit()
    call pladv(0)
    call plvpor(0.0_pl_test_flt,1.0_pl_test_flt,0.0_pl_test_flt,0.9_pl_test_flt)
    call plwind(-1.0_pl_test_flt,1.0_pl_test_flt,-1.0_pl_test_flt,1.5_pl_test_flt)
    call plw3d(1.0_pl_test_flt, 1.0_pl_test_flt, 1.2_pl_test_flt, xmin, &
         xmax, ymin, ymax, zmin, zmax, alt,az)
    call plbox3('bnstu',xaxis, 0.0_pl_test_flt, 0, &
         'bnstu',yaxis, 0.0_pl_test_flt, 0, &
         'bcdmnstuv',zaxis, 0.0_pl_test_flt, 0)
    call plmtex('t', 1.0_pl_test_flt, 0.5_pl_test_flt, 0.5_pl_test_flt,name)
    call plmesh(x, y, z, ior(opt, MAG_COLOR))
    call plend()

  end subroutine surf

  !---------------------------------------------------------------------!
  !                                                                     !
  !                           3D scatter plot                           !
  !                                                                     !
  !---------------------------------------------------------------------!

  subroutine scatter3D(X,Y,Z,xlabel,ylabel,zlabel,title)
    use kinds
    use plplot
    implicit none

    real(DP),dimension(:),intent(in) :: X,Y,Z
    character(len=*),optional,intent(in)        :: xlabel,ylabel,zlabel,title

    integer                :: opt=3,nxsub,nysub,nzsub
    real(DP)               :: xmin,xmax,ymin,ymax,zmin,&
         & zmax,xtick,ytick,ztick
    real(DP)               :: altitude    =  20.0_DP
    real(DP)               :: azimuth     = 298.0_DP
    real(DP)               :: azimuthStep =   0.2_DP
    real(DP)               :: basex  = 1.0_DP
    real(DP)               :: basey  = 1.0_DP
    real(DP)               :: height = 1.0_DP
    real(kind=pl_test_flt) :: alt = 33.0_pl_test_flt
    real(kind=pl_test_flt) :: az  = 24.0_pl_test_flt

    character(len=23) :: xaxis,yaxis,zaxis,name

    xaxis = ""
    yaxis = ""
    zaxis = ""
    name  = ""

    if(present(xlabel)) xaxis = xlabel
    if(present(ylabel)) yaxis = ylabel
    if(present(zlabel)) zaxis = zlabel
    if(present(title)) name = title

    xtick=0
    ytick=0
    ztick=0
    nxsub=0
    nysub=0
    nzsub=0

    xmin = minval(x)
    xmax = maxval(x)
    ymin = minval(y)
    ymax = maxval(y)
    zmin = minval(z)
    zmax = maxval(z)

    call plsdev("xwin")
    call plinit()
    call pladv(0)
    call plvpor(0.0_pl_test_flt,1.0_pl_test_flt,0.0_pl_test_flt,0.9_pl_test_flt)
    call plwind(-1.0_pl_test_flt,1.0_pl_test_flt,-1.0_pl_test_flt,1.5_pl_test_flt)
    call plw3d(1.0_pl_test_flt, 1.0_pl_test_flt, 1.2_pl_test_flt, xmin, &
         xmax, ymin, ymax, zmin, zmax, alt,az)
    call plbox3('bnstu',xaxis, 0.0_pl_test_flt, 0, &
         'bnstu',yaxis, 0.0_pl_test_flt, 0, &
         'bcdmnstuv',zaxis, 0.0_pl_test_flt, 0)
    call plmtex('t', 1.0_pl_test_flt, 0.5_pl_test_flt, 0.5_pl_test_flt,name)
    call plstring3(x,y,z,".")
    call plend()

  end subroutine scatter3D

end module PLplots
