!*****************************************************************************!
!                               NUMFORT                                       !
!              Authors: Anthony Kalaitzis, Curtis Abell,                      !
!                 Derek Leinweber and John Burkardt                           !             
!                                                                             !
! Numerical library for fortran. Use the following module titles              !
! to Move Around easily.                                                      !
!                                                                             !
! 1. QuadPack                                                                 !
! 2. NumFort                                                                  !
!                                                                             !
! Below contains a list of the function names inside numFort for easier       !
! navigation purposes. Jump to the bottom for numFort.                        !
!                                                                             !
! - factorial                                                                 !
! - bessel                                                                    !
! - deriv                                                                     !
! - Trace                                                                     !
! - inv                                                                       !
! - meshgrid                                                                  !
! - FFT                                                                       !
! - splinefit                                                                 !
! - splineval                                                                 !
! - polyfit                                                                   !
! - polyval                                                                   !
! - polyInt                                                                   !
! - minimize                                                                  !
! - linspace                                                                  !
! - PrintTime                                                                 !
! - guessZero                                                                 !
! - Newton1D                                                                  !
! - Euler                                                                     !
! - rk4                                                                       !
! - trapZ                                                                     !
! - integral                                                                  !
! - integralPV                                                                !
! - writeData                                                                 !
! - LoopTimer                                                                 !
! - linspace                                                                  !
! - pyplot                                                                    !
!                                                                             !
!*****************************************************************************!



!==============================================================================
!##############################################################################
!==============================================================================


module Quadpack
  use kinds
  implicit none

  !******************************************************************************
  !
  ! 1. Introduction
  !
  !    Quadpack is a fortran subroutine package for the numerical
  !    computation of definite one-dimensional integrals. It originated
  !    from a joint project of R. Piessens and E. de Doncker (Appl.
  !    Math. and Progr. Div.- k.u.leuven, Belgium), C. Ueberhuber (Inst.
  !    Fuer Math.- techn.u.wien, Austria), and D. Kahaner (Nation. Bur.
  !    of Standards- Washington D.C., U.S.A.).
  !
  ! 2. Survey
  !
  !    - qags : Is an integrator based on globally adaptive interval
  !             subdivision in connection with extrapolation (de doncker,
  !             1978) by the epsilon algorithm (wynn, 1956).
  !
  !    - qagp : Serves the same purposes as qags, but also allows
  !             for eventual user-supplied information, i.e. the
  !             abscissae of internal singularities, discontinuities
  !             and other difficulties of the integrand function.
  !             The algorithm is a modification of that in qags.
  !
  !    - qagi : Handles integration over infinite intervals. The
  !             infinite range is mapped onto a finite interval and
  !             then the same strategy as in qags is applied.
  !
  !    - qawo : Is a routine for the integration of cos(omega*x)*f(x)
  !             or sin(omega*x)*f(x) over a finite interval (a,b).
  !             Omega is specified by the user
  !             the rule evaluation component is based on the
  !             modified clenshaw-curtis technique.
  !             An adaptive subdivision scheme is used connected with
  !             an extrapolation procedure, which is a modification
  !             of that in qags and provides the possibility to deal
  !             even with singularities in f.
  !
  !    - qawf : Calculates the fourier cosine or fourier sine
  !             transform of f(x), for user-supplied interval (a,
  !             infinity), omega, and f. The procedure of qawo is
  !             used on successive finite intervals, and convergence
  !             acceleration by means of the epsilon algorithm (wynn,
  !             1956) is applied to the series of the integral
  !             contributions.
  !
  !    - qaws : Integrates w(x)*f(x) over (a,b) with a < b finite,
  !             and   w(x) = ((x-a)**alfa)*((b-x)**beta)*v(x)
  !             where v(x) = 1 or log(x-a) or log(b-x)
  !                            or log(x-a)*log(b-x)
  !             and   alfa > (-1), beta > (-1).
  !             The user specifies a, b, alfa, beta and the type of
  !             the function v.
  !             A globally adaptive subdivision strategy is applied,
  !             with modified clenshaw-curtis integration on the
  !             subintervals which contain a or b.
  !
  !    - qawc : Computes the cauchy principal value of f(x)/(x-c)
  !             over a finite interval (a,b) and for user-determined c.
  !             The strategy is globally adaptive, and modified
  !             clenshaw-curtis integration is used on the subranges
  !             which contain the point x = c.
  !
  !  Each of the routines above also has a "more detailed" version
  !    with a name ending in e, as qage.  These provide more
  !    information and control than the easier versions.
  !
  !
  !  The preceeding routines are all automatic.  That is, the user
  !      inputs his problem and an error tolerance.  The routine
  !      attempts to perform the integration to within the requested
  !      absolute or relative error.
  !  There are, in addition, a number of non-automatic integrators.
  !      These are most useful when the problem is such that the
  !      user knows that a fixed rule will provide the accuracy
  !      required.  Typically they return an error estimate but make
  !      no attempt to satisfy any particular input error request.
  !
  !      qk15
  !      qk21
  !      qk31
  !      qk41
  !      qk51
  !      qk61
  !           estimate the integral on [a,b] using 15, 21,..., 61
  !           point rule and return an error estimate.
  !      qk15i 15 point rule for (semi)infinite interval.
  !      qk15w 15 point rule for special singular weight functions.
  !      qc25c 25 point rule for cauchy principal values
  !      qc25o 25 point rule for sin/cos integrand.
  !      qmomo integrates k-th degree chebychev polynomial times
  !            function with various explicit singularities.
  !
  ! 3. Guidelines for the use of quadpack
  !
  !    Here it is not our purpose to investigate the question when
  !    automatic quadrature should be used. We shall rather attempt
  !    to help the user who already made the decision to use quadpack,
  !    with selecting an appropriate routine or a combination of
  !    several routines for handling their problem.
  !
  !    For both quadrature over finite and over infinite intervals,
  !    one of the first questions to be answered by the user is
  !    related to the amount of computer time they want to spend,
  !    versus their -own- time which would be needed, for example, for
  !    manual subdivision of the interval or other analytic
  !    manipulations.
  !
  !    (1) The user may not care about computer time, or not be
  !        willing to do any analysis of the problem. Especially when
  !        only one or a few integrals must be calculated, this attitude
  !        can be perfectly reasonable. In this case it is clear that
  !        either the most sophisticated of the routines for finite
  !        intervals, qags, must be used, or its analogue for infinite
  !        intervals, qagi. These routines are able to cope with
  !        rather difficult, even with improper integrals.
  !        This way of proceeding may be expensive. But the integrator
  !        is supposed to give you an answer in return, with additional
  !        information in the case of a failure, through its error
  !        estimate and flag. Yet it must be stressed that the programs
  !        cannot be totally reliable.
  !
  !    (2) The user may want to examine the integrand function.
  !        If bad local difficulties occur, such as a discontinuity, a
  !        singularity, derivative singularity or high peak at one or
  !        more points within the interval, the first advice is to
  !        split up the interval at these points. The integrand must
  !        then be examinated over each of the subintervals separately,
  !        so that a suitable integrator can be selected for each of
  !        them. If this yields problems involving relative accuracies
  !        to be imposed on -finite- subintervals, one can make use of
  !        qagp, which must be provided with the positions of the local
  !        difficulties. However, if strong singularities are present
  !        and a high accuracy is requested, application of qags on the
  !        subintervals may yield a better result.
  !
  !        For quadrature over finite intervals we thus dispose of qags
  !        and
  !        - qng for well-behaved integrands,
  !        - qag for functions with an oscillating behavior of a non
  !          specific type,
  !        - qawo for functions, eventually singular, containing a
  !          factor cos(omega*x) or sin(omega*x) where omega is known,
  !        - qaws for integrands with algebraico-logarithmic end point
  !          singularities of known type,
  !        - qawc for cauchy principal values.
  !
  !        remark
  !
  !        On return, the work arrays in the argument lists of the
  !        adaptive integrators contain information about the interval
  !        subdivision process and hence about the integrand behavior:
  !        the end points of the subintervals, the local integral
  !        contributions and error estimates, and eventually other
  !        characteristics. For this reason, and because of its simple
  !        globally adaptive nature, the routine qag in particular is
  !        well-suited for integrand examination. Difficult spots can
  !        be located by investigating the error estimates on the
  !        subintervals.
  !
  !        For infinite intervals we provide only one general-purpose
  !        routine, qagi. It is based on the qags algorithm applied
  !        after a transformation of the original interval into (0,1).
  !        Yet it may eventuate that another type of transformation is
  !        more appropriate, or one might prefer to break up the
  !        original interval and use qagi only on the infinite part
  !        and so on. These kinds of actions suggest a combined use of
  !        different quadpack integrators. Note that, when the only
  !        difficulty is an integrand singularity at the finite
  !        integration limit, it will in general not be necessary to
  !        break up the interval, as qagi deals with several types of
  !        singularity at the boundary point of the integration range.
  !        It also handles slowly convergent improper integrals, on
  !        the condition that the integrand does not oscillate over
  !        the entire infinite interval. If it does we would advise
  !        to sum succeeding positive and negative contributions to
  !        the integral -e.g. integrate between the zeros- with one
  !        or more of the finite-range integrators, and apply
  !        convergence acceleration eventually by means of quadpack
  !        subroutine qelg which implements the epsilon algorithm.
  !        Such quadrature problems include the fourier transform as
  !        a special case. Yet for the latter we have an automatic
  !        integrator available, qawf.
  !
  ! 4. Summary
  !        qags :: general purpose integration over a finite range
  !        qagi :: general purpose integration over semi-infinite or
  !                infinite ranges
  !        qagp :: general purpose integration with singularities,
  !                discontinuities and other difficulties
  !
  !   If computational speed is an issue try replacing qags with
  !        qag         :: for integrands with oscillatory behaviour
  !                       of unknown type
  !                       (however in general qags will is probably better)
  !        qawc        :: for Cauchy Principle value integration of the
  !                       type f(x)/(x-c)
  !        qawf        :: for Fourier sin or cos integrals
  !        qawo        :: if integrand has sin or cosine factor
  !        qaws        :: for integrands with algebraico-logarithmic
  !                       endpoint singularities of a known type
  !        qng or quad :: for very well behaved integrands


  private
  public :: qag, qagi, qagp, qags, qawc, qawf, qawo, qaws,&
       & qk15, qk21, qk31, qk41, qk51, qk61, qng, quad

contains
  !==============================================================================
  !##############################################################################
  !==============================================================================
  !
  subroutine qag(f,a,b,epsabs,epsrel,key,result,abserr,neval,ifail)

    ! QAG approximates an integral over a finite interval.
    !
    ! Discussion:
    !
    !   The routine calculates an approximation RESULT to a definite integral
    !     I = integral of F over (A,B),
    !   hopefully satisfying
    !     || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
    !
    !   QAG is a simple globally adaptive integrator using the strategy of
    !   Aind (Piessens, 1973).  It is possible to choose between 6 pairs of
    !   Gauss-Kronrod quadrature formulae for the rule evaluation component.
    !   The pairs of high degree of precision are suitable for handling
    !   integration difficulties due to a strongly oscillating integrand.
    !
    ! Reference:
    !
    !   R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
    !   QUADPACK, a Subroutine Package for Automatic Integration,
    !   Springer Verlag, 1983
    !
    ! Parameters: SEE qage
    !
    ! Local parameters:
    !
    !   LIMIT is the maximum number of subintervals allowed in
    !   the subdivision process of QAGE.

    !  REAL(dp), external    :: f
    interface
       function f(x)
         use kinds
         implicit none
         real(dp), intent(in) :: x
         real(dp)             :: f
       end function f
    end interface
    integer,  intent(in)  :: key
    real(dp), intent(in)  :: a,b,epsabs,epsrel
    integer,  intent(out) :: neval,ifail
    real(dp), intent(out) :: result,abserr

    ! Declare local parameters
    integer, parameter :: limit = 500
    integer :: iord(limit),last
    real(dp) :: alist(limit),blist(limit),elist(limit),rlist(limit)

    call qage (f,a,b,epsabs,epsrel,key,limit,result,abserr,neval, &
         ifail,alist,blist,rlist,elist,iord,last)

  end subroutine qag

  !===============================================================================
  !###############################################################################
  !===============================================================================
  subroutine qage(f,a,b,epsabs,epsrel,key, limit,result,abserr,neval, &
       ifail,alist,blist,rlist,elist,iord,last)

    !******************************************************************************
    !
    !! QAGE estimates a definite integral.
    !
    !
    ! Discussion:
    !
    !   The routine calculates an approximation RESULT to a definite integral
    !     I = integral of F over (A,B),
    !   hopefully satisfying
    !     || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
    !
    ! Reference:
    !
    !   R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
    !   QUADPACK, a SUBROUTINE Package for Automatic Integration,
    !   Springer Verlag, 1983
    !
    ! Parameters
    !         ON ENTRY
    !            f      - double precision
    !                     function subprogram defining the integrand
    !                     function f(x). the actual name for f needs to be
    !                     declared e x t e r n a l in the driver program.
    !
    !            a      - double precision
    !                     lower limit of integration
    !
    !            b      - double precision
    !                     upper limit of integration
    !
    !            epsabs - double precision
    !                     absolute accuracy requested
    !            epsrel - double precision
    !                     relative accuracy requested
    !                     if  epsabs.le.0
    !                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
    !                     the routine will end with ifail = 6.
    !
    !            key    - integer
    !                     key for choice of local integration rule
    !                     a gauss-kronrod pair is used with
    !                          7 - 15 points if key.lt.2,
    !                         10 - 21 points if key = 2,
    !                         15 - 31 points if key = 3,
    !                         20 - 41 points if key = 4,
    !                         25 - 51 points if key = 5,
    !                         30 - 61 points if key.gt.5.
    !
    !            limit  - integer
    !                     gives an upperbound on the number of subintervals
    !                     in the partition of (a,b), limit.ge.1.
    !
    !         ON RETURN
    !            result - double precision
    !                     approximation to the integral
    !
    !            abserr - double precision
    !                     estimate of the modulus of the absolute error,
    !                     which should equal or exceed abs(i-result)
    !
    !            neval  - integer
    !                     number of integrand evaluations
    !
    !          ifail    - integer
    !                   ifail = 0 normal and reliable termination of the
    !                             routine. it is assumed that the requested
    !                             accuracy has been achieved.
    !                   ifail.gt.0 abnormal termination of the routine
    !                             the estimates for result and error are
    !                             less reliable. it is assumed that the
    !                             requested accuracy has not been achieved.
    !            error messages
    !                   ifail = 1 maximum number of subdivisions allowed
    !                             has been achieved. one can allow more
    !                             subdivisions by increasing the value
    !                             of limit.
    !                             however, if this yields no improvement it
    !                             is rather advised to analyze the integrand
    !                             in order to determine the integration
    !                             difficulties. if the position of a local
    !                             difficulty can be determined(e.g.
    !                             singularity, discontinuity within the
    !                             interval) one will probably gain from
    !                             splitting up the interval at this point
    !                             and calling the integrator on the
    !                             subranges. if possible, an appropriate
    !                             special-purpose integrator should be used
    !                             which is designed for handling the type of
    !                             difficulty involved.
    !                         = 2 the occurrence of roundoff error is
    !                             detected, which prevents the requested
    !                             tolerance from being achieved.
    !                         = 3 extremely bad integrand behaviour occurs
    !                             at some points of the integration
    !                             interval.
    !                         = 6 the input is invalid, because
    !                             (epsabs.le.0 and
    !                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
    !                             result, abserr, neval, last, rlist(1) ,
    !                             elist(1) and iord(1) are set to zero.
    !                             alist(1) and blist(1) are set to a and b
    !                             respectively.
    !
    !            alist   - double precision
    !                      vector of dimension at least limit, the first
    !                       last  elements of which are the left
    !                      end points of the subintervals in the partition
    !                      of the given integration range (a,b)
    !
    !            blist   - double precision
    !                      vector of dimension at least limit, the first
    !                       last  elements of which are the right
    !                      end points of the subintervals in the partition
    !                      of the given integration range (a,b)
    !
    !            rlist   - double precision
    !                      vector of dimension at least limit, the first
    !                       last  elements of which are the
    !                      integral approximations on the subintervals
    !
    !            elist   - double precision
    !                      vector of dimension at least limit, the first
    !                       last  elements of which are the moduli of the
    !                      absolute error estimates on the subintervals
    !
    !            iord    - integer
    !                      vector of dimension at least limit, the first k
    !                      elements of which are pointers to the
    !                      error estimates over the subintervals,
    !                      such that elist(iord(1)), ...,
    !                      elist(iord(k)) form a decreasing sequence,
    !                      with k = last if last.le.(limit/2+2), and
    !                      k = limit+1-last otherwise
    !
    !            last    - integer
    !                      number of subintervals actually produced in the
    !                      subdivision process

    !  REAL(dp), external    :: f
    interface
       function f(x)
         use kinds
         implicit none
         real(dp), intent(in) :: x
         real(dp)             :: f
       end function f
    end interface
    integer,  intent(in)  :: key,limit
    real(dp), intent(in)  :: a,b,epsabs,epsrel
    integer,  intent(out) :: neval,ifail,iord(limit),last
    real(dp), intent(out) :: result,abserr, &
         alist(limit),blist(limit),elist(limit),rlist(limit)

    ! Declare local variables
    real(dp) :: area,area1,area12,area2,a1,a2,b1,b2,c,defabs,defab1,defab2,errbnd, &
         errmax,error1,error2,erro12,errsum,resabs
    integer  :: iroff1,iroff2,keyf,maxerr,nrmax

    ! Test on validity of parameters.
    ifail    = 0
    neval    = 0
    last     = 0
    result   = 0.0_dp
    abserr   = 0.0_dp
    alist(1) = a
    blist(1) = b
    rlist(1) = 0.0_dp
    elist(1) = 0.0_dp
    iord(1)  = 0

    if ( epsabs < 0.0_dp .and. epsrel < 0.0_dp ) then
       ifail = 6
       return
    endif

    ! First approximation to the integral.
    keyf = key
    keyf = max ( keyf, 1 )
    keyf = min ( keyf, 6 )

    c = keyf
    neval = 0

    if ( keyf == 1 ) then
       call qk15 (f,a,b,result,abserr,defabs,resabs)
    elseif ( keyf == 2 ) then
       call qk21 (f,a,b,result,abserr,defabs,resabs)
    elseif ( keyf == 3 ) then
       call qk31 (f,a,b,result,abserr,defabs,resabs)
    elseif ( keyf == 4 ) then
       call qk41 (f,a,b,result,abserr,defabs,resabs)
    elseif ( keyf == 5 ) then
       call qk51 (f,a,b,result,abserr,defabs,resabs)
    elseif ( keyf == 6 ) then
       call qk61 (f,a,b,result,abserr,defabs,resabs)
    endif

    last = 1
    rlist(1) = result
    elist(1) = abserr
    iord(1) = 1

    ! Test on accuracy.
    errbnd = max ( epsabs, epsrel * abs ( result ) )

    if ( abserr <= 5.0e+01_dp * epsilon ( defabs ) * defabs .and. &
         abserr > errbnd ) then
       ifail = 2
    endif

    if ( limit == 1 ) ifail = 1

    if ( ifail /= 0 .or. &
         ( abserr <= errbnd .and. abserr /= resabs ) .or. &
         abserr == 0.0_dp ) then

       if ( keyf /= 1 ) then
          neval = (10*keyf+1) * (2*neval+1)
       else
          neval = 30 * neval + 15
       endif
       return
    endif

    ! Initialization.
    errmax = abserr
    maxerr = 1
    area   = result
    errsum = abserr
    nrmax  = 1
    iroff1 = 0
    iroff2 = 0

    do last = 2, limit

       ! Bisect the subinterval with the largest error estimate.
       a1 = alist(maxerr)
       b1 = 0.5_dp * ( alist(maxerr) + blist(maxerr) )
       a2 = b1
       b2 = blist(maxerr)

       if     ( keyf == 1 ) then
          call qk15 (f,a1,b1,area1,error1,resabs,defab1)
       elseif ( keyf == 2 ) then
          call qk21 (f,a1,b1,area1,error1,resabs,defab1)
       elseif ( keyf == 3 ) then
          call qk31 (f,a1,b1,area1,error1,resabs,defab1)
       elseif ( keyf == 4 ) then
          call qk41 (f,a1,b1,area1,error1,resabs,defab1)
       elseif ( keyf == 5 ) then
          call qk51 (f,a1,b1,area1,error1,resabs,defab1)
       elseif ( keyf == 6 ) then
          call qk61 (f,a1,b1,area1,error1,resabs,defab1)
       endif

       if     ( keyf == 1 ) then
          call qk15 (f,a2,b2,area2,error2,resabs,defab2)
       elseif ( keyf == 2 ) then
          call qk21 (f,a2,b2,area2,error2,resabs,defab2)
       elseif ( keyf == 3 ) then
          call qk31 (f,a2,b2,area2,error2,resabs,defab2)
       elseif ( keyf == 4 ) then
          call qk41 (f,a2,b2,area2,error2,resabs,defab2)
       elseif ( keyf == 5 ) then
          call qk51 (f,a2,b2,area2,error2,resabs,defab2)
       elseif ( keyf == 6 ) then
          call qk61 (f,a2,b2,area2,error2,resabs,defab2)
       endif

       ! Improve previous approximations to integral and error and
       ! test for accuracy.
       neval  = neval + 1
       area12 = area1 + area2
       erro12 = error1 + error2
       errsum = errsum + erro12 - errmax
       area = area + area12 - rlist(maxerr)

       if ( defab1 /= error1 .and. defab2 /= error2 ) then

          if ( abs ( rlist(maxerr) - area12 ) <= 1.0e-05 * abs ( area12 ) &
               .and. erro12 >= 9.9e-01_dp * errmax ) iroff1 = iroff1 + 1

          if ( last > 10 .and. erro12 > errmax ) iroff2 = iroff2 + 1
       endif

       rlist(maxerr) = area1
       rlist(last) = area2
       errbnd = max ( epsabs, epsrel * abs ( area ) )

       ! Test for roundoff error and eventually set error flag.
       if ( errsum > errbnd ) then

          if ( iroff1 >= 6 .or. iroff2 >= 20 ) ifail = 2

          ! Set error flag in the case that the number of subintervals
          ! equals limit.
          if ( last == limit ) ifail = 1

          ! Set error flag in the case of bad integrand behavior
          ! at a point of the integration range.
          if ( max ( abs ( a1 ), abs ( b2 ) ) <= ( 1.0_dp + c * 1.0e+03_dp * &
               epsilon ( a1 ) ) * ( abs ( a2 ) + 1.0e+04 * tiny ( a2 ) ) ) ifail = 3
       endif

       ! Append the newly-created intervals to the list.
       if ( error2 <= error1 ) then
          alist(last)   = a2
          blist(maxerr) = b1
          blist(last)   = b2
          elist(maxerr) = error1
          elist(last)   = error2
       else
          alist(maxerr) = a2
          alist(last)   = a1
          blist(last)   = b1
          rlist(maxerr) = area2
          rlist(last)   = area1
          elist(maxerr) = error2
          elist(last)   = error1
       endif

       ! Call QSORT to maintain the descending ordering
       ! in the list of error estimates and select the subinterval
       ! with the largest error estimate (to be bisected next).
       call qsort ( limit, last, maxerr, errmax, elist, iord, nrmax )

       if ( ifail /= 0 .or. errsum <= errbnd ) exit
    enddo

    ! Compute final result.
    result = sum ( rlist(1:last) )
    abserr = errsum

    if ( keyf /= 1 ) then
       neval = ( 10 * keyf + 1 ) * ( 2 * neval + 1 )
    else
       neval = 30 * neval + 15
    endif

  end subroutine qage

  !===============================================================================
  !###############################################################################
  !===============================================================================
  subroutine qagi(f,bound,inf,epsabs,epsrel,result,abserr,neval,ifail)

    !******************************************************************************
    !
    ! QAGI estimates an integral over a semi-infinite or infinite interval.
    !
    !
    ! Discussion:
    !
    !   The routine calculates an approximation RESULT to a definite integral
    !     I = integral of F over (A, +Infinity),
    !   or
    !     I = integral of F over (-Infinity,A)
    !   or
    !     I = integral of F over (-Infinity,+Infinity),
    !   hopefully satisfying
    !     || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
    !
    ! Reference:
    !
    !   R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
    !   QUADPACK, a Subroutine Package for Automatic Integration,
    !   Springer Verlag, 1983
    !
    ! Parameters:
    !
    !   Input, external real F, the name of the function routine, of the form
    !     function f ( x )
    !     real f
    !     real x
    !   which evaluates the integrand function.
    !
    !   Input, real BOUND, the value of the finite endpoint of the integration
    !   range, if any, that is, if INF is 1 or -1.
    !
    !   Input, integer INF, indicates the type of integration range.
    !   1:  (  BOUND,    +Infinity),
    !   -1: ( -Infinity,  BOUND),
    !   2:  ( -Infinity, +Infinity).
    !
    !   Input, real EPSABS, EPSREL, the absolute and relative accuracy requested.
    !
    !   Output, real RESULT, the estimated value of the integral.
    !
    !   Output, real ABSERR, an estimate of || I - RESULT ||.
    !
    !   Output, integer NEVAL, the number of times the integral was evaluated.
    !
    !   Output, integer IFAIL, error indicator.
    !   0, normal and reliable termination of the routine.  It is assumed that
    !     the requested accuracy has been achieved.
    !   > 0,  abnormal termination of the routine.  The estimates for result
    !     and error are less reliable.  It is assumed that the requested
    !     accuracy has not been achieved.
    !   1, maximum number of subdivisions allowed has been achieved.  One can
    !     allow more subdivisions by increasing the data value of LIMIT in QAGI
    !     (and taking the according dimension adjustments into account).
    !     However, if this yields no improvement it is advised to analyze the
    !     integrand in order to determine the integration difficulties.  If the
    !     position of a local difficulty can be determined (e.g. singularity,
    !     discontinuity within the interval) one will probably gain from
    !     splitting up the interval at this point and calling the integrator
    !     on the subranges.  If possible, an appropriate special-purpose
    !     integrator should be used, which is designed for handling the type
    !     of difficulty involved.
    !   2, the occurrence of roundoff error is detected, which prevents the
    !     requested tolerance from being achieved.  The error may be
    !     under-estimated.
    !   3, extremely bad integrand behavior occurs at some points of the
    !     integration interval.
    !   4, the algorithm does not converge.  Roundoff error is detected in the
    !     extrapolation table.  It is assumed that the requested tolerance
    !     cannot be achieved, and that the returned result is the best which
    !     can be obtained.
    !   5, the integral is probably divergent, or slowly convergent.  It must
    !     be noted that divergence can occur with any other value of IFAIL.
    !   6, the input is invalid, because INF /= 1 and INF /= -1 and INF /= 2, or
    !     epsabs < 0 and epsrel < 0.  result, abserr, neval are set to zero.
    !
    ! Local parameters:
    !
    !           the dimension of rlist2 is determined by the value of
    !           limexp in QEXTR.
    !
    !          alist     - list of left end points of all subintervals
    !                      considered up to now
    !          blist     - list of right end points of all subintervals
    !                      considered up to now
    !          rlist(i)  - approximation to the integral over
    !                      (alist(i),blist(i))
    !          rlist2    - array of dimension at least (limexp+2),
    !                      containing the part of the epsilon table
    !                      which is still needed for further computations
    !          elist(i)  - error estimate applying to rlist(i)
    !          maxerr    - pointer to the interval with largest error
    !                      estimate
    !          errmax    - elist(maxerr)
    !          erlast    - error on the interval currently subdivided
    !                      (before that subdivision has taken place)
    !          area      - sum of the integrals over the subintervals
    !          errsum    - sum of the errors over the subintervals
    !          errbnd    - requested accuracy max(epsabs,epsrel*
    !                      abs(result))
    !          *****1    - variable for the left subinterval
    !          *****2    - variable for the right subinterval
    !          last      - index for subdivision
    !          nres      - number of calls to the extrapolation routine
    !          numrl2    - number of elements currently in rlist2. if an
    !                      appropriate approximation to the compounded
    !                      integral has been obtained, it is put in
    !                      rlist2(numrl2) after numrl2 has been increased
    !                      by one.
    !          small     - length of the smallest interval considered up
    !                      to now, multiplied by 1.5
    !          erlarg    - sum of the errors over the intervals larger
    !                      than the smallest interval considered up to now
    !          extrap    - logical variable denoting that the routine
    !                      is attempting to perform extrapolation. i.e.
    !                      before subdividing the smallest interval we
    !                      try to decrease the value of erlarg.
    !          noext     - logical variable denoting that extrapolation
    !                      is no longer allowed (true-value)

    !  REAL(dp), external    :: f
    interface
       function f(x)
         use kinds
         implicit none
         real(dp), intent(in) :: x
         real(dp)             :: f
       end function f
    end interface
    integer,  intent(in)  :: inf
    real(dp), intent(in)  :: bound,epsabs,epsrel
    integer,  intent(out) :: neval,ifail
    real(dp), intent(out) :: result,abserr

    ! Declare local varibles
    integer, parameter :: limit = 2000
    logical  :: extrap,noext
    integer  :: id,ierro,iord(limit),iroff1,iroff2,iroff3,jupbnd,k,ksgn,ktmin,last, &
         maxerr,nres,nrmax,numrl2
    real(dp) :: abseps,alist(limit),area,area1,area12,area2,a1,a2,blist(limit),boun, &
         b1,b2,correc,defabs,defab1,defab2,dres,elist(limit),erlarg,erlast,   &
         errbnd,errmax,error1,error2,erro12,errsum,ertest,resabs,reseps,res3la(3), &
         rlist(limit),rlist2(52),small

    ! Test on validity of parameters.
    ifail    = 0
    neval    = 0
    last     = 0
    result   = 0.0_dp
    abserr   = 0.0_dp
    alist(1) = 0.0_dp
    blist(1) = 1.0_dp
    rlist(1) = 0.0_dp
    elist(1) = 0.0_dp
    iord(1)  = 0

    if ( epsabs < 0.0_dp .and. epsrel < 0.0_dp ) then
       ifail = 6
       return
    endif

    ! First approximation to the integral.

    ! Determine the interval to be mapped onto (0,1).
    ! If INF = 2 the integral is computed as i = i1+i2, where
    ! i1 = integral of f over (-infinity,0),
    ! i2 = integral of f over (0,+infinity).
    if ( inf == 2 ) then
       boun = 0.0_dp
    else
       boun = bound
    endif

    call qk15i(f,boun,inf,0.0_dp,1.0_dp,result,abserr,defabs,resabs)

    ! Test on accuracy.
    last     = 1
    rlist(1) = result
    elist(1) = abserr
    iord(1)  = 1
    dres     = abs ( result )
    errbnd   = max ( epsabs, epsrel * dres )

    if ( abserr <= 100.0_dp * epsilon ( defabs ) * defabs .and. &
         abserr > errbnd ) ifail = 2

    if ( limit == 1 ) ifail = 1

    if ( ifail /= 0 .or. (abserr <= errbnd .and. abserr /= resabs ) .or. &
         abserr == 0.0_dp ) go to 130

    ! Initialization
    rlist2(1) = result
    errmax    = abserr
    maxerr    = 1
    area      = result
    errsum    = abserr
    abserr    = huge ( abserr )
    nrmax     = 1
    nres      = 0
    ktmin     = 0
    numrl2    = 2
    extrap    = .false.
    noext     = .false.
    ierro     = 0
    iroff1    = 0
    iroff2    = 0
    iroff3    = 0

    if ( dres >= ( 1.0_dp - 5.0e+01_dp * epsilon ( defabs ) ) * defabs ) then
       ksgn = 1
    else
       ksgn = -1
    endif

    do last = 2, limit

       ! Bisect the subinterval with nrmax-th largest error estimate.
       a1 = alist(maxerr)
       b1 = 5.0e-01_dp * ( alist(maxerr) + blist(maxerr) )
       a2 = b1
       b2 = blist(maxerr)
       erlast = errmax
       call qk15i (f,boun,inf,a1,b1,area1,error1,resabs,defab1)
       call qk15i (f,boun,inf,a2,b2,area2,error2,resabs,defab2)

       ! Improve previous approximations to integral and error
       ! and test for accuracy.
       area12 = area1 + area2
       erro12 = error1 + error2
       errsum = errsum + erro12 - errmax
       area   = area + area12 - rlist(maxerr)

       if ( defab1 /= error1 .and. defab2 /= error2 ) then

          if ( abs ( rlist(maxerr) - area12 ) <= 1.0e-05 * abs ( area12 ) &
               .and. erro12 >= 9.9e-01_dp * errmax ) then

             if ( extrap ) iroff2 = iroff2 + 1

             if ( .not. extrap ) iroff1 = iroff1 + 1
          endif

          if ( last > 10 .and. erro12 > errmax ) iroff3 = iroff3 + 1
       endif

       rlist(maxerr) = area1
       rlist(last) = area2
       errbnd = max ( epsabs, epsrel * abs ( area ) )

       ! Test for roundoff error and eventually set error flag.
       if ( iroff1 + iroff2 >= 10 .or. iroff3 >= 20 ) ifail = 2

       if ( iroff2 >= 5 ) ierro = 3

       ! Set error flag in the case that the number of subintervals equals LIMIT.
       if ( last == limit ) ifail = 1

       ! Set error flag in the case of bad integrand behavior
       ! at some points of the integration range.
       if ( max ( abs(a1), abs(b2) ) <= (1.0_dp + 1.0e+03_dp * epsilon ( a1 ) ) * &
            ( abs(a2) + 1.0e+03_dp * tiny ( a2 ) )) ifail = 4

       ! Append the newly-created intervals to the list.
       if ( error2 <= error1 ) then
          alist(last)   = a2
          blist(maxerr) = b1
          blist(last)   = b2
          elist(maxerr) = error1
          elist(last)   = error2
       else
          alist(maxerr) = a2
          alist(last)   = a1
          blist(last)   = b1
          rlist(maxerr) = area2
          rlist(last)   = area1
          elist(maxerr) = error2
          elist(last)   = error1
       endif

       ! Call QSORT to maintain the descending ordering
       ! in the list of error estimates and select the subinterval
       ! with NRMAX-th largest error estimate (to be bisected next).
       call qsort(limit,last,maxerr,errmax,elist,iord,nrmax)

       if ( errsum <= errbnd ) go to 115

       if ( ifail /= 0 ) exit

       if ( last == 2 ) then
          small = 3.75e-01_dp
          erlarg = errsum
          ertest = errbnd
          rlist2(2) = area
          cycle
       endif

       if ( noext ) cycle

       erlarg = erlarg-erlast

       if ( abs(b1-a1) > small ) erlarg = erlarg+erro12

       ! Test whether the interval to be bisected next is the
       ! smallest interval.
       if ( .not. extrap ) then
          if ( abs(blist(maxerr)-alist(maxerr)) > small ) cycle
          extrap = .true.
          nrmax = 2
       endif

       if ( ierro == 3 .or. erlarg <= ertest ) go to 60

       ! The smallest interval has the largest error.
       ! before bisecting decrease the sum of the errors over the
       ! larger intervals (erlarg) and perform extrapolation.
       id     = nrmax
       jupbnd = last

       if ( last > (2+limit/2) ) then
          jupbnd = limit + 3 - last
       endif

       do k = id, jupbnd
          maxerr = iord(nrmax)
          errmax = elist(maxerr)
          if ( abs ( blist(maxerr) - alist(maxerr) ) > small ) go to 90
          nrmax = nrmax + 1
       enddo

       ! Extrapolate.
60     continue

       numrl2 = numrl2 + 1
       rlist2(numrl2) = area
       call qextr(numrl2,rlist2,reseps,abseps,res3la,nres)
       ktmin = ktmin+1

       if ( ktmin > 5.and.abserr < 1.0e-03_dp*errsum ) ifail = 5

       if ( abseps < abserr ) then
          ktmin  = 0
          abserr = abseps
          result = reseps
          correc = erlarg
          ertest = max ( epsabs, epsrel*abs(reseps) )
          if ( abserr <= ertest ) exit
       endif

       ! Prepare bisection of the smallest interval.
       if ( numrl2 == 1 ) noext = .true.

       if ( ifail == 5 ) exit

       maxerr = iord(1)
       errmax = elist(maxerr)
       nrmax  = 1
       extrap = .false.
       small  = small*5.0e-01_dp
       erlarg = errsum
90     continue
    enddo

    ! Set final result and error estimate.
    if ( abserr == huge ( abserr ) ) go to 115

    if ( (ifail+ierro) == 0 ) go to 110

    if ( ierro == 3 ) abserr = abserr+correc

    if ( ifail == 0 ) ifail = 3

    if ( result /= 0.0_dp .and. area /= 0.0_dp) go to 105
    if ( abserr > errsum) go to 115
    if ( area == 0.0_dp)  go to 130

    go to 110

105 continue
    if ( abserr / abs(result) > errsum / abs(area) ) go to 115

    ! Test on divergence
110 continue

    if ( ksgn == (-1) .and. &
         max ( abs(result), abs(area) ) <=  defabs * 1.0e-02_dp) go to 130

    if ( 1.0e-02_dp > (result/area) .or. &
         (result/area) > 1.0e+02_dp .or. &
         errsum > abs(area)) ifail = 6

    go to 130

    ! Compute global integral sum.
115 continue

    result = sum ( rlist(1:last) )
    abserr = errsum
130 continue

    neval = 30*last-15
    if ( inf == 2 ) neval = 2*neval

    if ( ifail > 2 ) ifail = ifail - 1

  end subroutine qagi

  !===============================================================================
  !###############################################################################
  !===============================================================================
  subroutine qagp(f,a,b,npts2,points,epsabs,epsrel,result,abserr,neval,ifail)

    !******************************************************************************
    !
    !! QAGP computes a definite integral.
    !
    !
    ! Discussion:
    !
    !   The routine calculates an approximation RESULT to a definite integral
    !     I = integral of F over (A,B),
    !   hopefully satisfying
    !     || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
    !
    !   Interior break points of the integration interval,
    !   where local difficulties of the integrand may occur, such as
    !   singularities or discontinuities, are provided by the user.
    !
    ! Reference:
    !
    !   R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
    !   QUADPACK, a Subroutine Package for Automatic Integration,
    !   Springer Verlag, 1983
    !
    ! Parameters:
    !
    !   Input, external real F, the name of the function routine, of the form
    !     function f ( x )
    !     real f
    !     real x
    !   which evaluates the integrand function.
    !
    !   Input, real A, B, the limits of integration.
    !
    !   Input, integer NPTS2, the number of user-supplied break points within
    !   the integration range, plus 2.  NPTS2 must be at least 2.
    !
    !   Input/output, real POINTS(NPTS2), contains the user provided interior
    !   breakpoints in entries 1 through NPTS2-2.  If these points are not
    !   in ascending order on input, they will be sorted.
    !
    !   Input, real EPSABS, EPSREL, the absolute and relative accuracy requested.
    !
    !   Output, real RESULT, the estimated value of the integral.
    !
    !   Output, real ABSERR, an estimate of || I - RESULT ||.
    !
    !   Output, integer NEVAL, the number of times the integral was evaluated.
    !
    !         ifail    - integer
    !                  ifail = 0 normal and reliable termination of the
    !                            routine. it is assumed that the requested
    !                            accuracy has been achieved.
    !                  ifail > 0 abnormal termination of the routine.
    !                            the estimates for integral and error are
    !                            less reliable. it is assumed that the
    !                            requested accuracy has not been achieved.
    !                  ifail = 1 maximum number of subdivisions allowed
    !                            has been achieved. one can allow more
    !                            subdivisions by increasing the data value
    !                            of limit in qagp(and taking the according
    !                            dimension adjustments into account).
    !                            however, if this yields no improvement
    !                            it is advised to analyze the integrand
    !                            in order to determine the integration
    !                            difficulties. if the position of a local
    !                            difficulty can be determined (i.e.
    !                            singularity, discontinuity within the
    !                            interval), it should be supplied to the
    !                            routine as an element of the vector
    !                            points. if necessary, an appropriate
    !                            special-purpose integrator must be used,
    !                            which is designed for handling the type
    !                            of difficulty involved.
    !                        = 2 the occurrence of roundoff error is
    !                            detected, which prevents the requested
    !                            tolerance from being achieved.
    !                            the error may be under-estimated.
    !                        = 3 extremely bad integrand behavior occurs
    !                            at some points of the integration
    !                            interval.
    !                        = 4 the algorithm does not converge. roundoff
    !                            error is detected in the extrapolation
    !                            table. it is presumed that the requested
    !                            tolerance cannot be achieved, and that
    !                            the returned result is the best which
    !                            can be obtained.
    !                        = 5 the integral is probably divergent, or
    !                            slowly convergent. it must be noted that
    !                            divergence can occur with any other value
    !                            of ifail > 0.
    !                        = 6 the input is invalid because
    !                            npts2 < 2 or
    !                            break points are specified outside
    !                            the integration range or
    !                            epsabs < 0 and epsrel < 0,
    !                            or limit < npts2.
    !                            result, abserr, neval are set to zero.
    !
    ! Local parameters:
    !
    !           the dimension of rlist2 is determined by the value of
    !           limexp in QEXTR (rlist2 should be of dimension
    !           (limexp+2) at least).
    !
    !          alist     - list of left end points of all subintervals
    !                      considered up to now
    !          blist     - list of right end points of all subintervals
    !                      considered up to now
    !          rlist(i)  - approximation to the integral over
    !                      (alist(i),blist(i))
    !          rlist2    - array of dimension at least limexp+2
    !                      containing the part of the epsilon table which
    !                      is still needed for further computations
    !          elist(i)  - error estimate applying to rlist(i)
    !          maxerr    - pointer to the interval with largest error
    !                      estimate
    !          errmax    - elist(maxerr)
    !          erlast    - error on the interval currently subdivided
    !                      (before that subdivision has taken place)
    !          area      - sum of the integrals over the subintervals
    !          errsum    - sum of the errors over the subintervals
    !          errbnd    - requested accuracy max(epsabs,epsrel*
    !                      abs(result))
    !          *****1    - variable for the left subinterval
    !          *****2    - variable for the right subinterval
    !          last      - index for subdivision
    !          nres      - number of calls to the extrapolation routine
    !          numrl2    - number of elements in rlist2. if an appropriate
    !                      approximation to the compounded integral has
    !                      obtained, it is put in rlist2(numrl2) after
    !                      numrl2 has been increased by one.
    !          erlarg    - sum of the errors over the intervals larger
    !                      than the smallest interval considered up to now
    !          extrap    - logical variable denoting that the routine
    !                      is attempting to perform extrapolation. i.e.
    !                      before subdividing the smallest interval we
    !                      try to decrease the value of erlarg.
    !          noext     - logical variable denoting that extrapolation is
    !                      no longer allowed (true-value)

    !  REAL(dp), external    :: f
    interface
       function f(x)
         use kinds
         implicit none
         real(dp), intent(in) :: x
         real(dp)             :: f
       end function f
    end interface
    integer,  intent(in)  :: npts2
    real(dp), intent(in)  :: a,b,epsabs,epsrel,points(npts2)
    integer,  intent(out) :: neval,ifail
    real(dp), intent(out) :: result,abserr

    ! Declare local variables
    integer, parameter :: limit = 500

    logical  :: extrap,noext
    integer  :: i,id,ierro,ind1,ind2,iord(limit),ip1,iroff1,iroff2,iroff3,j,jlow, &
         jupbnd,k,ksgn,ktmin,last,levcur,level(limit),levmax,maxerr,ndin(40), &
         nint,npts,nres,nrmax,numrl2
    real(dp) :: abseps,alist(limit),area,area1,area12,area2,a1,a2,blist(limit),b1,b2, &
         correc,defabs,defab1,defab2,dres,elist(limit),erlarg,erlast,errbnd,   &
         errmax,error1,erro12,error2,errsum,ertest,pts(40),resa,resabs,reseps, &
         res3la(3),rlist(limit),rlist2(52),sign

    ! Test on validity of parameters.
    ifail    = 0
    neval    = 0
    last     = 0
    result   = 0.0_dp
    abserr   = 0.0_dp
    alist(1) = a
    blist(1) = b
    rlist(1) = 0.0_dp
    elist(1) = 0.0_dp
    iord(1)  = 0
    level(1) = 0
    npts = npts2-2

    if ( npts2 < 2 ) then
       ifail = 6
       return
    elseif ( limit <= npts .or. (epsabs < 0.0_dp.and. &
         epsrel < 0.0_dp) ) then
       ifail = 6
       return
    endif

    ! If any break points are provided, sort them into an
    ! ascending sequence.
    if ( a > b ) then
       sign = -1.0_dp
    else
       sign = +1.0_dp
    endif

    pts(1) = min ( a,b)

    do i = 1, npts
       pts(i+1) = points(i)
    enddo

    pts(npts+2) = max ( a,b)
    nint = npts+1
    a1   = pts(1)

    if ( npts /= 0 ) then
       do i = 1, nint
          ip1 = i+1
          do j = ip1, nint+1
             if ( pts(i) > pts(j) ) then
                call r_swap ( pts(i), pts(j) )
             endif
          enddo
       enddo
       if ( pts(1) /= min ( a, b ) .or. pts(nint+1) /= max ( a,b) ) then
          ifail = 6
          return
       endif
    endif

    ! Compute first integral and error approximations.
    resabs = 0.0_dp

    do i = 1, nint
       b1 = pts(i+1)
       call qk21 ( f, a1, b1, area1, error1, defabs, resa )
       abserr = abserr+error1
       result = result+area1
       ndin(i) = 0

       if ( error1 == resa .and. error1 /= 0.0_dp ) ndin(i) = 1

       resabs = resabs + defabs
       level(i) = 0
       elist(i) = error1
       alist(i) = a1
       blist(i) = b1
       rlist(i) = area1
       iord(i)  = i
       a1 = b1
    enddo

    errsum = 0.0_dp

    do i = 1, nint
       if ( ndin(i) == 1 ) elist(i) = abserr
       errsum = errsum + elist(i)
    enddo

    ! Test on accuracy.
    last   = nint
    neval  = 21 * nint
    dres   = abs ( result )
    errbnd = max ( epsabs, epsrel * dres )

    if ( abserr <= 1.0e+02_dp * epsilon ( resabs ) * resabs .and. &
         abserr > errbnd ) ifail = 2

    if ( nint /= 1 ) then
       do i = 1, npts
          jlow = i+1
          ind1 = iord(i)
          do j = jlow, nint
             ind2 = iord(j)
             if ( elist(ind1) <= elist(ind2) ) then
                ind1 = ind2
                k = j
             endif
          enddo
          if ( ind1 /= iord(i) ) then
             iord(k) = iord(i)
             iord(i) = ind1
          endif
       enddo
       if ( limit < npts2 ) then
          ifail = 1
       endif
    endif

    if ( ifail /= 0 .or. abserr <= errbnd ) return

    ! Initialization
    rlist2(1) = result
    maxerr = iord(1)
    errmax = elist(maxerr)
    area   = result
    nrmax  = 1
    nres   = 0
    numrl2 = 1
    ktmin  = 0
    extrap = .false.
    noext  = .false.
    erlarg = errsum
    ertest = errbnd
    levmax = 1
    iroff1 = 0
    iroff2 = 0
    iroff3 = 0
    ierro  = 0
    abserr = huge ( abserr )

    if ( dres >= ( 1.0_dp - 0.5_dp * epsilon ( resabs ) ) * resabs ) then
       ksgn = 1
    else
       ksgn = -1
    endif

    do last = npts2, limit

       ! Bisect the subinterval with the nrmax-th largest error estimate.
       levcur = level(maxerr)+1
       a1 = alist(maxerr)
       b1 = 0.5_dp * ( alist(maxerr) + blist(maxerr) )
       a2 = b1
       b2 = blist(maxerr)
       erlast = errmax
       call qk21(f,a1,b1,area1,error1,resa,defab1)
       call qk21(f,a2,b2,area2,error2,resa,defab2)

       ! Improve previous approximations to integral and error
       ! and test for accuracy.
       neval = neval+42
       area12 = area1+area2
       erro12 = error1+error2
       errsum = errsum+erro12-errmax
       area = area+area12-rlist(maxerr)

       if ( defab1 /= error1 .and. defab2 /= error2 ) then

          if ( abs(rlist(maxerr)-area12) <= 1.0e-05*abs(area12) .and. &
               erro12 >= 9.9e-01_dp*errmax ) then

             if ( extrap ) then
                iroff2 = iroff2+1
             else
                iroff1 = iroff1+1
             endif
          endif

          if ( last > 10 .and. erro12 > errmax ) iroff3 = iroff3 + 1
       endif

       level(maxerr) = levcur
       level(last)   = levcur
       rlist(maxerr) = area1
       rlist(last)   = area2
       errbnd = max ( epsabs, epsrel * abs ( area ) )

       ! Test for roundoff error and eventually set error flag.
       if ( iroff1 + iroff2 >= 10 .or. iroff3 >= 20 ) ifail = 2

       if ( iroff2 >= 5 ) ierro = 3

       ! Set error flag in the case that the number of subintervals
       ! equals limit.
       if ( last == limit ) ifail = 1

       ! Set error flag in the case of bad integrand behavior
       ! at a point of the integration range
       if ( max ( abs(a1),abs(b2)) <= (1.0_dp+1.0e+03_dp* epsilon ( a1 ) )* &
            ( abs(a2) + 1.0e+03_dp * tiny ( a2 ) ) ) ifail = 4

       ! Append the newly-created intervals to the list.
       if ( error2 <= error1 ) then
          alist(last)   = a2
          blist(maxerr) = b1
          blist(last)   = b2
          elist(maxerr) = error1
          elist(last)   = error2
       else
          alist(maxerr) = a2
          alist(last)   = a1
          blist(last)   = b1
          rlist(maxerr) = area2
          rlist(last)   = area1
          elist(maxerr) = error2
          elist(last)   = error1
       endif

       ! Call QSORT to maintain the descending ordering
       ! in the list of error estimates and select the subinterval
       ! with nrmax-th largest error estimate (to be bisected next).
       call qsort(limit,last,maxerr,errmax,elist,iord,nrmax)

       if ( errsum <= errbnd ) go to 190

       if ( ifail /= 0 ) exit

       if ( noext ) cycle

       erlarg = erlarg - erlast

       if ( levcur+1 <= levmax ) erlarg = erlarg + erro12

       ! Test whether the interval to be bisected next is the
       ! smallest interval.
       if ( .not. extrap ) then
          if ( level(maxerr)+1 <= levmax ) cycle
          extrap = .true.
          nrmax = 2
       endif

       ! The smallest interval has the largest error.
       ! Before bisecting decrease the sum of the errors over the
       ! larger intervals (erlarg) and perform extrapolation.
       if ( ierro /= 3 .and. erlarg > ertest ) then
          id = nrmax
          jupbnd = last
          if ( last > (2+limit/2) ) jupbnd = limit+3-last
          do k = id, jupbnd
             maxerr = iord(nrmax)
             errmax = elist(maxerr)
             if ( level(maxerr)+1 <= levmax ) go to 160
             nrmax = nrmax+1
          enddo
       endif

       ! Perform extrapolation.
       numrl2 = numrl2+1
       rlist2(numrl2) = area
       if ( numrl2 <= 2 ) go to 155
       call qextr ( numrl2, rlist2, reseps, abseps, res3la, nres )
       ktmin = ktmin+1

       if ( ktmin > 5 .and. abserr < 1.0e-03_dp*errsum ) ifail = 5

       if ( abseps < abserr ) then
          ktmin  = 0
          abserr = abseps
          result = reseps
          correc = erlarg
          ertest = max ( epsabs,epsrel*abs(reseps))
          if ( abserr < ertest ) exit
       endif

       ! Prepare bisection of the smallest interval.
       if ( numrl2 == 1 ) noext = .true.

       if ( ifail >= 5 ) exit

155    continue
       maxerr = iord(1)
       errmax = elist(maxerr)
       nrmax  = 1
       extrap = .false.
       levmax = levmax+1
       erlarg = errsum
160    continue
    enddo

    ! Set the final result.
    if ( abserr == huge ( abserr ) ) go to 190
    if ( (ifail+ierro) == 0 ) go to 180

    if ( ierro == 3 ) abserr = abserr+correc

    if ( ifail == 0 ) ifail = 3

    if ( result /= 0.0_dp.and.area /= 0.0_dp ) go to 175
    if ( abserr > errsum ) go to 190
    if ( area == 0.0_dp ) go to 210
    go to 180

175 continue

    if ( abserr/abs(result) > errsum/abs(area) ) go to 190

    ! Test on divergence.
180 continue

    if ( ksgn == (-1) .and. max ( abs(result),abs(area)) <=  &
         resabs*1.0e-02_dp ) go to 210

    if ( 1.0e-02_dp > (result/area) .or. (result/area) > 1.0e+02_dp .or. &
         errsum > abs(area) ) ifail = 6
    go to 210

    ! Compute global integral sum.
190 continue
    result = sum ( rlist(1:last) )
    abserr = errsum

210 continue
    if ( ifail > 2 ) ifail = ifail - 1

    result = result * sign

  end subroutine qagp

  !===============================================================================
  !###############################################################################
  !===============================================================================
  subroutine qags(f,a,b,epsabs,epsrel,result,abserr,neval,ifail)
    !
    !******************************************************************************
    !
    !! QAGS estimates the integral of a function.
    !
    !
    ! Discussion:
    !
    !   The routine calculates an approximation RESULT to a definite integral
    !     I = integral of F over (A,B),
    !   hopefully satisfying
    !     || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
    !
    ! Reference:
    !
    !   R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
    !   QUADPACK, a Subroutine Package for Automatic Integration,
    !   Springer Verlag, 1983
    !
    ! Parameters:
    !
    !   Input, external real F, the name of the function routine, of the form
    !     function f ( x )
    !     real f
    !     real x
    !   which evaluates the integrand function.
    !
    !   Input, real A, B, the limits of integration.
    !
    !   Input, real EPSABS, EPSREL, the absolute and relative accuracy requested.
    !
    !   Output, real RESULT, the estimated value of the integral.
    !
    !   Output, real ABSERR, an estimate of || I - RESULT ||.
    !
    !   Output, integer NEVAL, the number of times the integral was evaluated.
    !
    !   Output, integer IFAIL, error flag.
    !                  ifail = 0 normal and reliable termination of the
    !                            routine. it is assumed that the requested
    !                            accuracy has been achieved.
    !                  ifail > 0 abnormal termination of the routine
    !                            the estimates for integral and error are
    !                            less reliable. it is assumed that the
    !                            requested accuracy has not been achieved.
    !                        = 1 maximum number of subdivisions allowed
    !                            has been achieved. one can allow more sub-
    !                            divisions by increasing the data value of
    !                            limit in qags (and taking the according
    !                            dimension adjustments into account).
    !                            however, if this yields no improvement
    !                            it is advised to analyze the integrand
    !                            in order to determine the integration
    !                            difficulties. if the position of a
    !                            local difficulty can be determined (e.g.
    !                            singularity, discontinuity within the
    !                            interval) one will probably gain from
    !                            splitting up the interval at this point
    !                            and calling the integrator on the sub-
    !                            ranges. if possible, an appropriate
    !                            special-purpose integrator should be used,
    !                            which is designed for handling the type
    !                            of difficulty involved.
    !                        = 2 the occurrence of roundoff error is detec-
    !                            ted, which prevents the requested
    !                            tolerance from being achieved.
    !                            the error may be under-estimated.
    !                        = 3 extremely bad integrand behavior occurs
    !                            at some  points of the integration
    !                            interval.
    !                        = 4 the algorithm does not converge. roundoff
    !                            error is detected in the extrapolation
    !                            table. it is presumed that the requested
    !                            tolerance cannot be achieved, and that the
    !                            returned result is the best which can be
    !                            obtained.
    !                        = 5 the integral is probably divergent, or
    !                            slowly convergent. it must be noted that
    !                            divergence can occur with any other value
    !                            of ifail.
    !                        = 6 the input is invalid, because
    !                            epsabs < 0 and epsrel < 0,
    !                            result, abserr and neval are set to zero.
    !
    ! Local Parameters:
    !
    !          alist     - list of left end points of all subintervals
    !                      considered up to now
    !          blist     - list of right end points of all subintervals
    !                      considered up to now
    !          rlist(i)  - approximation to the integral over
    !                      (alist(i),blist(i))
    !          rlist2    - array of dimension at least limexp+2 containing
    !                      the part of the epsilon table which is still
    !                      needed for further computations
    !          elist(i)  - error estimate applying to rlist(i)
    !          maxerr    - pointer to the interval with largest error
    !                      estimate
    !          errmax    - elist(maxerr)
    !          erlast    - error on the interval currently subdivided
    !                      (before that subdivision has taken place)
    !          area      - sum of the integrals over the subintervals
    !          errsum    - sum of the errors over the subintervals
    !          errbnd    - requested accuracy max(epsabs,epsrel*
    !                      abs(result))
    !          *****1    - variable for the left interval
    !          *****2    - variable for the right interval
    !          last      - index for subdivision
    !          nres      - number of calls to the extrapolation routine
    !          numrl2    - number of elements currently in rlist2. if an
    !                      appropriate approximation to the compounded
    !                      integral has been obtained it is put in
    !                      rlist2(numrl2) after numrl2 has been increased
    !                      by one.
    !          small     - length of the smallest interval considered
    !                      up to now, multiplied by 1.5
    !          erlarg    - sum of the errors over the intervals larger
    !                      than the smallest interval considered up to now
    !          extrap    - logical variable denoting that the routine is
    !                      attempting to perform extrapolation i.e. before
    !                      subdividing the smallest interval we try to
    !                      decrease the value of erlarg.
    !          noext     - logical variable denoting that extrapolation
    !                      is no longer allowed (true value)

    !  REAL(dp), external    :: f
    interface
       function f(x)
         use kinds
         implicit none
         real(dp), intent(in) :: x
         real(dp)             :: f
       end function f
    end interface
    real(dp), intent(in)  :: a,b,epsabs,epsrel
    integer,  intent(out) :: neval,ifail
    real(dp), intent(out) :: result,abserr

    ! Declare local variables
    integer, parameter :: limit = 500
    logical  :: extrap,noext
    integer  :: id,ierro,iord(limit),iroff1,iroff2,iroff3,jupbnd,k,ksgn,ktmin, &
         last,maxerr,nres,nrmax,numrl2
    real(dp) :: abseps,alist(limit),area,area1,area12,area2,a1,a2,blist(limit),b1,b2, &
         correc,defabs,defab1,defab2,dres,elist(limit),erlarg,erlast,errbnd,   &
         errmax,error1,error2,erro12,errsum,ertest,resabs,reseps,res3la(3),    &
         rlist(limit),rlist2(52),small

    ! The dimension of rlist2 is determined by the value of
    ! limexp in QEXTR (rlist2 should be of dimension (limexp+2) at least).

    ! Test on validity of parameters.
    ifail    = 0
    neval    = 0
    last     = 0
    result   = 0.0_dp
    abserr   = 0.0_dp
    alist(1) = a
    blist(1) = b
    rlist(1) = 0.0_dp
    elist(1) = 0.0_dp

    if ( epsabs < 0.0_dp .and. epsrel < 0.0_dp ) then
       ifail = 6
       return
    endif

    ! First approximation to the integral.
    ierro = 0
    call qk21(f,a,b,result,abserr,defabs,resabs)

    ! Test on accuracy.
    dres = abs ( result )
    errbnd = max ( epsabs, epsrel * dres )
    last = 1
    rlist(1) = result
    elist(1) = abserr
    iord(1) = 1

    if ( abserr <= 1.0e+02_dp * epsilon ( defabs ) * defabs .and. &
         abserr > errbnd ) ifail = 2

    if ( limit == 1 ) ifail = 1

    if ( ifail /= 0 .or. (abserr <= errbnd .and. abserr /= resabs ) .or. &
         abserr == 0.0_dp ) go to 140

    ! Initialization.
    rlist2(1) = result
    errmax = abserr
    maxerr = 1
    area   = result
    errsum = abserr
    abserr = huge ( abserr )
    nrmax  = 1
    nres   = 0
    numrl2 = 2
    ktmin  = 0
    extrap = .false.
    noext  = .false.
    iroff1 = 0
    iroff2 = 0
    iroff3 = 0

    if ( dres >= (1.0_dp-5.0e+01_dp* epsilon ( defabs ) )*defabs ) then
       ksgn = 1
    else
       ksgn = -1
    endif

    do last = 2, limit

       ! Bisect the subinterval with the nrmax-th largest error estimate.
       a1 = alist(maxerr)
       b1 = 5.0e-01_dp*(alist(maxerr)+blist(maxerr))
       a2 = b1
       b2 = blist(maxerr)
       erlast = errmax
       call qk21(f,a1,b1,area1,error1,resabs,defab1)
       call qk21(f,a2,b2,area2,error2,resabs,defab2)

       ! Improve previous approximations to integral and error
       ! and test for accuracy.
       area12 = area1+area2
       erro12 = error1+error2
       errsum = errsum+erro12-errmax
       area = area+area12-rlist(maxerr)

       if ( defab1 == error1 .or. defab2 == error2 ) go to 15

       if ( abs ( rlist(maxerr) - area12) > 1.0e-05 * abs(area12) &
            .or. erro12 < 9.9e-01_dp * errmax ) go to 10

       if ( extrap ) then
          iroff2 = iroff2+1
       else
          iroff1 = iroff1+1
       endif

10     continue

       if ( last > 10.and.erro12 > errmax ) iroff3 = iroff3+1

15     continue

       rlist(maxerr) = area1
       rlist(last) = area2
       errbnd = max ( epsabs,epsrel*abs(area))

       ! Test for roundoff error and eventually set error flag.
       if ( iroff1+iroff2 >= 10 .or. iroff3 >= 20 ) ifail = 2

       if ( iroff2 >= 5 ) ierro = 3

       ! Set error flag in the case that the number of subintervals
       ! equals limit.
       if ( last == limit ) ifail = 1

       ! Set error flag in the case of bad integrand behavior
       ! at a point of the integration range.
       if ( max ( abs(a1),abs(b2)) <= (1.0_dp+1.0e+03_dp* epsilon ( a1 ) )* &
            (abs(a2)+1.0e+03_dp* tiny ( a2 ) ) ) ifail = 4

       ! Append the newly-created intervals to the list.
       if ( error2 <= error1 ) then
          alist(last) = a2
          blist(maxerr) = b1
          blist(last) = b2
          elist(maxerr) = error1
          elist(last) = error2
       else
          alist(maxerr) = a2
          alist(last) = a1
          blist(last) = b1
          rlist(maxerr) = area2
          rlist(last) = area1
          elist(maxerr) = error2
          elist(last) = error1
       endif

       ! Call QSORT to maintain the descending ordering
       ! in the list of error estimates and select the subinterval
       ! with nrmax-th largest error estimate (to be bisected next).
       call qsort(limit,last,maxerr,errmax,elist,iord,nrmax)

       if ( errsum <= errbnd ) go to 115

       if ( ifail /= 0 ) exit

       if ( last == 2 ) go to 80
       if ( noext ) go to 90

       erlarg = erlarg-erlast

       if ( abs(b1-a1) > small ) erlarg = erlarg+erro12

       if ( extrap ) go to 40

       ! Test whether the interval to be bisected next is the
       ! smallest interval.
       if ( abs(blist(maxerr)-alist(maxerr)) > small ) go to 90
       extrap = .true.
       nrmax = 2

40     continue

       ! The smallest interval has the largest error.
       ! Before bisecting decrease the sum of the errors over the
       ! larger intervals (erlarg) and perform extrapolation.
       if ( ierro /= 3 .and. erlarg > ertest ) then
          id = nrmax
          jupbnd = last
          if ( last > (2+limit/2) ) then
             jupbnd = limit+3-last
          endif
          do k = id, jupbnd
             maxerr = iord(nrmax)
             errmax = elist(maxerr)
             if ( abs(blist(maxerr)-alist(maxerr)) > small ) go to 90
             nrmax = nrmax+1
          enddo
       endif

       ! Perform extrapolation.
60     continue

       numrl2 = numrl2+1
       rlist2(numrl2) = area
       call qextr(numrl2,rlist2,reseps,abseps,res3la,nres)
       ktmin = ktmin+1

       if ( ktmin > 5 .and. abserr < 1.0e-03_dp * errsum ) ifail = 5

       if ( abseps < abserr ) then
          ktmin  = 0
          abserr = abseps
          result = reseps
          correc = erlarg
          ertest = max ( epsabs,epsrel*abs(reseps))
          if ( abserr <= ertest ) exit
       endif

       ! Prepare bisection of the smallest interval.
       if ( numrl2 == 1 ) noext = .true.

       if ( ifail == 5 ) exit

       maxerr = iord(1)
       errmax = elist(maxerr)
       nrmax  = 1
       extrap = .false.
       small  = small*5.0e-01_dp
       erlarg = errsum
       go to 90

80     continue

       small  = abs(b-a)*3.75e-01_dp
       erlarg = errsum
       ertest = errbnd
       rlist2(2) = area

90     continue
    enddo

    ! Set final result and error estimate.
    if ( abserr == huge ( abserr ) ) go to 115
    if ( ifail+ierro == 0 ) go to 110

    if ( ierro == 3 ) abserr = abserr+correc

    if ( ifail == 0 ) ifail = 3
    if ( result /= 0.0_dp.and.area /= 0.0_dp ) go to 105
    if ( abserr > errsum ) go to 115
    if ( area == 0.0_dp ) go to 130
    go to 110

105 continue

    if ( abserr/abs(result) > errsum/abs(area) ) go to 115

    ! Test on divergence.
110 continue

    if ( ksgn == (-1).and.max ( abs(result),abs(area)) <=  &
         defabs*1.0e-02_dp ) go to 130

    if ( 1.0e-02_dp > (result/area) .or. (result/area) > 1.0e+02_dp &
         .or. errsum > abs(area) ) ifail = 6

    go to 130

    ! Compute global integral sum.
115 continue

    result = sum ( rlist(1:last) )
    abserr = errsum

130 continue

    if ( ifail > 2 ) ifail = ifail-1

140 continue

    neval = 42*last-21

  end subroutine qags

  !===============================================================================
  !###############################################################################
  !===============================================================================
  subroutine qawc(f,a,b,c,epsabs,epsrel,result,abserr,neval,ifail)
    !
    !******************************************************************************
    !
    !! QAWC computes a Cauchy principal value.
    !
    !
    ! Discussion:
    !
    !   The routine calculates an approximation RESULT to a Cauchy principal
    !   value
    !     I = integral of F*W over (A,B),
    !   with
    !     W(X) = 1 / (X-C),
    !   with C distinct from A and B, hopefully satisfying
    !     || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
    !
    ! Reference:
    !
    !   R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
    !   QUADPACK, a Subroutine Package for Automatic Integration,
    !   Springer Verlag, 1983
    !
    ! Parameters:
    !
    !   Input, external real F, the name of the function routine, of the form
    !     function f ( x )
    !     real f
    !     real x
    !   which evaluates the integrand function.
    !
    !   Input, real A, B, the limits of integration.
    !
    !   Input, real C, a parameter in the weight function, which must
    !   not be equal to A or B.
    !
    !   Input, real EPSABS, EPSREL, the absolute and relative accuracy requested.
    !
    !   Output, real RESULT, the estimated value of the integral.
    !
    !   Output, real ABSERR, an estimate of || I - RESULT ||.
    !
    !   Output, integer NEVAL, the number of times the integral was evaluated.
    !
    !           ifail    - integer
    !                    ifail = 0 normal and reliable termination of the
    !                            routine. it is assumed that the requested
    !                            accuracy has been achieved.
    !                    ifail > 0 abnormal termination of the routine
    !                            the estimates for integral and error are
    !                            less reliable. it is assumed that the
    !                            requested accuracy has not been achieved.
    !                    ifail = 1 maximum number of subdivisions allowed
    !                            has been achieved. one can allow more sub-
    !                            divisions by increasing the data value of
    !                            limit in qawc (and taking the according
    !                            dimension adjustments into account).
    !                            however, if this yields no improvement it
    !                            is advised to analyze the integrand in
    !                            order to determine the integration
    !                            difficulties. if the position of a local
    !                            difficulty can be determined (e.g.
    !                            singularity, discontinuity within the
    !                            interval one will probably gain from
    !                            splitting up the interval at this point
    !                            and calling appropriate integrators on the
    !                            subranges.
    !                        = 2 the occurrence of roundoff error is detec-
    !                            ted, which prevents the requested
    !                            tolerance from being achieved.
    !                        = 3 extremely bad integrand behavior occurs
    !                            at some points of the integration
    !                            interval.
    !                        = 6 the input is invalid, because
    !                            c = a or c = b or
    !                            epsabs < 0 and epsrel < 0,
    !                            result, abserr, neval are set to zero.
    !
    ! Local parameters:
    !
    !   LIMIT is the maximum number of subintervals allowed in the
    !   subdivision process of qawce. take care that limit >= 1.

    !  REAL(dp), external    :: f
    interface
       function f(x)
         use kinds
         implicit none
         real(dp), intent(in) :: x
         real(dp)             :: f
       end function f
    end interface
    real(dp), intent(in)  :: a,b,epsabs,epsrel
    integer,  intent(out) :: neval,ifail
    real(dp), intent(out) :: result,abserr

    ! Declare local variables
    integer, parameter :: limit = 500
    integer  :: iord(limit),last
    real(dp) :: alist(limit),blist(limit),elist(limit),c,rlist(limit)

    call qawce(f,a,b,c,epsabs,epsrel,limit,result,abserr,neval,ifail, &
         alist,blist,rlist,elist,iord,last)

  end subroutine qawc

  !===============================================================================
  !###############################################################################
  !===============================================================================
  subroutine qawce(f,a,b,c,epsabs,epsrel,limit,result,abserr,neval, &
       ifail,alist,blist,rlist,elist,iord,last)
    !
    !******************************************************************************
    !
    !! QAWCE computes a Cauchy principal value.
    !
    !
    ! Discussion:
    !
    !   The routine calculates an approximation RESULT to a Cauchy principal
    !   value
    !     I = integral of F*W over (A,B),
    !   with
    !     W(X) = 1 / ( X - C ),
    !   with C distinct from A and B, hopefully satisfying
    !     | I - RESULT | <= max ( EPSABS, EPSREL * |I| ).
    !
    ! Reference:
    !
    !   R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
    !   QUADPACK, a Subroutine Package for Automatic Integration,
    !   Springer Verlag, 1983
    !
    ! Parameters:
    !
    !   Input, external real F, the name of the function routine, of the form
    !     function f ( x )
    !     real f
    !     real x
    !   which evaluates the integrand function.
    !
    !   Input, real A, B, the limits of integration.
    !
    !   Input, real C, a parameter in the weight function, which cannot be
    !   equal to A or B.
    !
    !   Input, real EPSABS, EPSREL, the absolute and relative accuracy requested.
    !
    !   Input, integer LIMIT, the upper bound on the number of subintervals that
    !   will be used in the partition of [A,B].  LIMIT is typically 500.
    !
    !   Output, real RESULT, the estimated value of the integral.
    !
    !   Output, real ABSERR, an estimate of || I - RESULT ||.
    !
    !   Output, integer NEVAL, the number of times the integral was evaluated.
    !
    !         ifail    - integer
    !                  ifail = 0 normal and reliable termination of the
    !                            routine. it is assumed that the requested
    !                            accuracy has been achieved.
    !                  ifail > 0 abnormal termination of the routine
    !                            the estimates for integral and error are
    !                            less reliable. it is assumed that the
    !                            requested accuracy has not been achieved.
    !                  ifail = 1 maximum number of subdivisions allowed
    !                            has been achieved. one can allow more sub-
    !                            divisions by increasing the value of
    !                            limit. however, if this yields no
    !                            improvement it is advised to analyze the
    !                            integrand, in order to determine the
    !                            integration difficulties.  if the position
    !                            of a local difficulty can be determined
    !                            (e.g. singularity, discontinuity within
    !                            the interval) one will probably gain
    !                            from splitting up the interval at this
    !                            point and calling appropriate integrators
    !                            on the subranges.
    !                        = 2 the occurrence of roundoff error is detec-
    !                            ted, which prevents the requested
    !                            tolerance from being achieved.
    !                        = 3 extremely bad integrand behavior occurs
    !                            at some interior points of the integration
    !                            interval.
    !                        = 6 the input is invalid, because
    !                            c = a or c = b or
    !                            epsabs < 0 and epsrel < 0,
    !                            or limit < 1.
    !                            result, abserr, neval, rlist(1), elist(1),
    !                            iord(1) and last are set to zero.
    !                            alist(1) and blist(1) are set to a and b
    !                            respectively.
    !
    !   Workspace, real ALIST(LIMIT), BLIST(LIMIT), contains in entries 1
    !   through LAST the left and right ends of the partition subintervals.
    !
    !   Workspace, real RLIST(LIMIT), contains in entries 1 through LAST
    !   the integral approximations on the subintervals.
    !
    !   Workspace, real ELIST(LIMIT), contains in entries 1 through LAST
    !   the absolute error estimates on the subintervals.
    !
    !           iord    - integer
    !                     vector of dimension at least limit, the first k
    !                     elements of which are pointers to the error
    !                     estimates over the subintervals, so that
    !                     elist(iord(1)), ...,  elist(iord(k)) with
    !                     k = last if last <= (limit/2+2), and
    !                     k = limit+1-last otherwise, form a decreasing
    !                     sequence.
    !
    !           last    - integer
    !                     number of subintervals actually produced in
    !                     the subdivision process
    !
    ! Local parameters:
    !
    !          alist     - list of left end points of all subintervals
    !                      considered up to now
    !          blist     - list of right end points of all subintervals
    !                      considered up to now
    !          rlist(i)  - approximation to the integral over
    !                      (alist(i),blist(i))
    !          elist(i)  - error estimate applying to rlist(i)
    !          maxerr    - pointer to the interval with largest error
    !                      estimate
    !          errmax    - elist(maxerr)
    !          area      - sum of the integrals over the subintervals
    !          errsum    - sum of the errors over the subintervals
    !          errbnd    - requested accuracy max(epsabs,epsrel*
    !                      abs(result))
    !          *****1    - variable for the left subinterval
    !          *****2    - variable for the right subinterval
    !          last      - index for subdivision

    !  REAL(dp), external    :: f
    interface
       function f(x)
         use kinds
         implicit none
         real(dp), intent(in) :: x
         real(dp)             :: f
       end function f
    end interface
    integer,  intent(in)  :: limit
    real(dp), intent(in)  :: a,b,epsabs,epsrel
    integer,  intent(out) :: neval,ifail,iord(limit),last
    real(dp), intent(out) :: result,abserr, &
         alist(limit),blist(limit),elist(limit),rlist(limit)
    ! Declare local variables
    integer  :: iroff1,iroff2,krule,maxerr,nev,nrmax
    real(dp) :: aa,area,area1,area12,area2,a1,a2,bb,b1,b2,c,errbnd,errmax, &
         error1,error2,erro12,errsum

    ! Test on validity of parameters.
    ifail    = 0
    neval    = 0
    last     = 0
    alist(1) = a
    blist(1) = b
    rlist(1) = 0.0_dp
    elist(1) = 0.0_dp
    iord(1)  = 0
    result   = 0.0_dp
    abserr   = 0.0_dp

    if ( c == a  ) then
       ifail = 6
       return
    elseif ( c == b ) then
       ifail = 6
       return
    elseif ( epsabs < 0.0_dp .and. epsrel < 0.0_dp ) then
       ifail = 6
       return
    endif

    ! First approximation to the integral.
    if ( a <= b ) then
       aa = a
       bb = b
    else
       aa = b
       bb = a
    endif

    krule = 1
    call qc25c(f,aa,bb,c,result,abserr,krule,neval)
    last = 1
    rlist(1) = result
    elist(1) = abserr
    iord(1)  = 1
    alist(1) = a
    blist(1) = b

    ! Test on accuracy.
    errbnd = max ( epsabs, epsrel*abs(result) )

    if ( limit == 1 ) then
       ifail = 1
       go to 70
    endif

    if ( abserr < min ( 1.0e-02_dp*abs(result),errbnd)  ) go to 70

    ! Initialization
    alist(1) = aa
    blist(1) = bb
    rlist(1) = result
    errmax   = abserr
    maxerr   = 1
    area     = result
    errsum   = abserr
    nrmax    = 1
    iroff1   = 0
    iroff2   = 0

    do last = 2, limit

       ! Bisect the subinterval with nrmax-th largest error estimate.
       a1 = alist(maxerr)
       b1 = 5.0e-01_dp*(alist(maxerr)+blist(maxerr))
       b2 = blist(maxerr)
       if ( c <= b1 .and. c > a1 ) b1 = 5.0e-01_dp*(c+b2)
       if ( c > b1 .and. c < b2 ) b1 = 5.0e-01_dp*(a1+c)
       a2 = b1
       krule = 2

       call qc25c (f,a1,b1,c,area1,error1,krule,nev)
       neval = neval+nev

       call qc25c (f,a2,b2,c,area2,error2,krule,nev)
       neval = neval+nev

       ! Improve previous approximations to integral and error
       ! and test for accuracy.
       area12 = area1+area2
       erro12 = error1+error2
       errsum = errsum+erro12-errmax
       area = area+area12-rlist(maxerr)

       if ( abs(rlist(maxerr)-area12) < 1.0e-05*abs(area12) &
            .and.erro12 >= 9.9e-01_dp*errmax .and. krule == 0 ) &
            iroff1 = iroff1+1

       if ( last > 10.and.erro12 > errmax .and. krule == 0 ) iroff2 = iroff2+1

       rlist(maxerr) = area1
       rlist(last)   = area2
       errbnd = max ( epsabs,epsrel*abs(area))

       if ( errsum > errbnd ) then

          ! Test for roundoff error and eventually set error flag.
          if ( iroff1 >= 6 .and. iroff2 > 20 ) ifail = 2

          ! Set error flag in the case that number of interval
          ! bisections exceeds limit.
          if ( last == limit ) ifail = 1

          ! Set error flag in the case of bad integrand behavior at
          ! a point of the integration range.
          if ( max ( abs(a1), abs(b2) ) <= ( 1.0_dp + 1.0e+03_dp * epsilon ( a1 ) ) &
               *(abs(a2)+1.0e+03_dp* tiny ( a2 ) )) ifail = 3
       endif

       ! Append the newly-created intervals to the list.
       if ( error2 <= error1 ) then
          alist(last)   = a2
          blist(maxerr) = b1
          blist(last)   = b2
          elist(maxerr) = error1
          elist(last)   = error2
       else
          alist(maxerr) = a2
          alist(last)   = a1
          blist(last)   = b1
          rlist(maxerr) = area2
          rlist(last)   = area1
          elist(maxerr) = error2
          elist(last)   = error1
       endif

       ! Call QSORT to maintain the descending ordering
       ! in the list of error estimates and select the subinterval
       ! with NRMAX-th largest error estimate (to be bisected next).
       call qsort(limit,last,maxerr,errmax,elist,iord,nrmax)

       if ( ifail /= 0 .or. errsum <= errbnd ) exit
    enddo

    ! Compute final result.
    result = sum ( rlist(1:last) )
    abserr = errsum

70  continue

    if ( aa == b ) result = -result

  end subroutine qawce

  !===============================================================================
  !###############################################################################
  !===============================================================================
  subroutine qawf(f,a,omega,integr,epsabs,result,abserr,neval,ifail)
    !
    !******************************************************************************
    !
    !! QAWF computes Fourier integrals over the interval [ A, +Infinity ).
    !
    !
    ! Discussion:
    !
    !   The routine calculates an approximation RESULT to a definite integral
    !
    !     I = integral of F*COS(OMEGA*X)
    !   or
    !     I = integral of F*SIN(OMEGA*X)
    !
    !   over the interval [A,+Infinity), hopefully satisfying
    !
    !     || I - RESULT || <= EPSABS.
    !
    !   If OMEGA = 0 and INTEGR = 1, the integral is calculated by means
    !   of QAGI, and IFAIL has the meaning as described in the comments of QAGI.
    !
    ! Reference:
    !
    !   R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
    !   QUADPACK, a Subroutine Package for Automatic Integration,
    !   Springer Verlag, 1983
    !
    ! Parameters:
    !
    !   Input, external real F, the name of the function routine, of the form
    !     function f ( x )
    !     real f
    !     real x
    !   which evaluates the integrand function.
    !
    !   Input, real A, the lower limit of integration.
    !
    !   Input, real OMEGA, the parameter in the weight function.
    !
    !   Input, integer INTEGR, indicates which weight functions is used
    !   = 1, w(x) = cos(omega*x)
    !   = 2, w(x) = sin(omega*x)
    !
    !   Input, real EPSABS, the absolute accuracy requested.
    !
    !   Output, real RESULT, the estimated value of the integral.
    !
    !   Output, real ABSERR, an estimate of || I - RESULT ||.
    !
    !   Output, integer NEVAL, the number of times the integral was evaluated.
    !
    !         ifail    - integer
    !                  ifail = 0 normal and reliable termination of the
    !                            routine. it is assumed that the
    !                            requested accuracy has been achieved.
    !                  ifail > 0 abnormal termination of the routine.
    !                            the estimates for integral and error are
    !                            less reliable. it is assumed that the
    !                            requested accuracy has not been achieved.
    !                   if omega /= 0
    !                  ifail = 6 the input is invalid because
    !                            (integr /= 1 and integr /= 2) or
    !                             epsabs <= 0
    !                             result, abserr, neval, lst are set to
    !                             zero.
    !                        = 7 abnormal termination of the computation
    !                            of one or more subintegrals
    !                        = 8 maximum number of cycles allowed
    !                            has been achieved, i.e. of subintervals
    !                            (a+(k-1)c,a+kc) where
    !                            c = (2*int(abs(omega))+1)*pi/abs(omega),
    !                            for k = 1, 2, ...
    !                        = 9 the extrapolation table constructed for
    !                            convergence acceleration of the series
    !                            formed by the integral contributions
    !                            over the cycles, does not converge to
    !                            within the requested accuracy.
    !
    ! Local parameters:
    !
    !   Integer LIMLST, gives an upper bound on the number of cycles, LIMLST >= 3.
    !   if limlst < 3, the routine will end with ifail = 6.
    !
    !   Integer MAXP1, an upper bound on the number of Chebyshev moments which
    !   can be stored, i.e. for the intervals of lengths abs(b-a)*2**(-l),
    !   l = 0,1, ..., maxp1-2, maxp1 >= 1.  if maxp1 < 1, the routine will end
    !   with ifail = 6.

    !  REAL(dp), external    :: f
    interface
       function f(x)
         use kinds
         implicit none
         real(dp), intent(in) :: x
         real(dp)             :: f
       end function f
    end interface
    integer,  intent(in)  :: integr
    real(dp), intent(in)  :: a,omega,epsabs
    integer,  intent(out) :: neval,ifail
    real(dp), intent(out) :: result,abserr

    ! Declare local variables
    integer, parameter :: limit = 500,limlst = 50,maxp1 = 21
    integer  :: iord(limit),ierlst(limlst),last,lst,nnlog(limit)
    real(dp) :: alist(limit),blist(limit),chebmo(maxp1,25),elist(limit), &
         erlst(limlst),rlist(limit),rslst(limlst)

    ifail  = 6
    neval  = 0
    last   = 0
    result = 0.0_dp
    abserr = 0.0_dp

    if ( limlst < 3 .or. maxp1 < 1 ) return

    call qawfe(f,a,omega,integr,epsabs,limlst,limit,maxp1,result, &
         abserr,neval,ifail,rslst,erlst,ierlst,lst,alist,blist,rlist, &
         elist,iord,nnlog,chebmo)

  end subroutine qawf

  !===============================================================================
  !###############################################################################
  !===============================================================================
  subroutine qawfe(f,a,omega,integr,epsabs,limlst,limit,maxp1, &
       result,abserr,neval,ifail,rslst,erlst,ierlst,lst,alist,blist, &
       rlist,elist,iord,nnlog,chebmo )
    !
    !******************************************************************************
    !
    !! QAWFE computes Fourier integrals.
    !
    !
    ! Discussion:
    !
    !   The routine calculates an approximation RESULT to a definite integral
    !     I = integral of F*COS(OMEGA*X) or F*SIN(OMEGA*X) over (A,+Infinity),
    !   hopefully satisfying
    !     || I - RESULT || <= EPSABS.
    !
    ! Reference:
    !
    !   R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
    !   QUADPACK, a Subroutine Package for Automatic Integration,
    !   Springer Verlag, 1983
    !
    ! Parameters:
    !
    !   Input, external real F, the name of the function routine, of the form
    !     function f ( x )
    !     real f
    !     real x
    !   which evaluates the integrand function.
    !
    !   Input, real A, the lower limit of integration.
    !
    !   Input, real OMEGA, the parameter in the weight function.
    !
    !   Input, integer INTEGR, indicates which weight function is used
    !   = 1      w(x) = cos(omega*x)
    !   = 2      w(x) = sin(omega*x)
    !
    !   Input, real EPSABS, the absolute accuracy requested.
    !
    !   Input, integer LIMLST, an upper bound on the number of cycles.
    !   LIMLST must be at least 1.  In fact, if LIMLST < 3, the routine
    !   will end with IFAIL= 6.
    !
    !           limit  - integer
    !                    gives an upper bound on the number of
    !                    subintervals allowed in the partition of
    !                    each cycle, limit >= 1.
    !
    !           maxp1  - integer
    !                    gives an upper bound on the number of
    !                    Chebyshev moments which can be stored, i.e.
    !                    for the intervals of lengths abs(b-a)*2**(-l),
    !                    l=0,1, ..., maxp1-2, maxp1 >= 1
    !
    !   Output, real RESULT, the estimated value of the integral.
    !
    !   Output, real ABSERR, an estimate of || I - RESULT ||.
    !
    !   Output, integer NEVAL, the number of times the integral was evaluated.
    !
    !       ifail    - ifail = 0 normal and reliable termination of
    !                            the routine. it is assumed that the
    !                            requested accuracy has been achieved.
    !                  ifail > 0 abnormal termination of the routine
    !                            the estimates for integral and error
    !                            are less reliable. it is assumed that
    !                            the requested accuracy has not been
    !                            achieved.
    !                   if omega /= 0
    !                  ifail = 6 the input is invalid because
    !                            (integr /= 1 and integr /= 2) or
    !                             epsabs <= 0 or limlst < 3.
    !                             result, abserr, neval, lst are set
    !                             to zero.
    !                        = 7 bad integrand behavior occurs within
    !                            one or more of the cycles. location
    !                            and type of the difficulty involved
    !                            can be determined from the vector ierlst.
    !                            here lst is the number of cycles actually
    !                            needed (see below).
    !                            ierlst(k) = 1 the maximum number of
    !                                          subdivisions (= limit)
    !                                          has been achieved on the
    !                                          k th cycle.
    !                                      = 2 occurence of roundoff
    !                                          error is detected and
    !                                          prevents the tolerance
    !                                          imposed on the k th cycle
    !                                          from being acheived.
    !                                      = 3 extremely bad integrand
    !                                          behavior occurs at some
    !                                          points of the k th cycle.
    !                                      = 4 the integration procedure
    !                                          over the k th cycle does
    !                                          not converge (to within the
    !                                          required accuracy) due to
    !                                          roundoff in the
    !                                          extrapolation procedure
    !                                          invoked on this cycle. it
    !                                          is assumed that the result
    !                                          on this interval is the
    !                                          best which can be obtained.
    !                                      = 5 the integral over the k th
    !                                          cycle is probably divergent
    !                                          or slowly convergent. it
    !                                          must be noted that
    !                                          divergence can occur with
    !                                          any other value of
    !                                          ierlst(k).
    !                        = 8 maximum number of  cycles  allowed
    !                            has been achieved, i.e. of subintervals
    !                            (a+(k-1)c,a+kc) where
    !                            c = (2*int(abs(omega))+1)*pi/abs(omega),
    !                            for k = 1, 2, ..., lst.
    !                            one can allow more cycles by increasing
    !                            the value of limlst (and taking the
    !                            according dimension adjustments into
    !                            account).
    !                            examine the array iwork which contains
    !                            the error flags over the cycles, in order
    !                            to eventual look for local integration
    !                            difficulties.
    !                            if the position of a local difficulty can
    !                            be determined (e.g. singularity,
    !                            discontinuity within the interval)
    !                            one will probably gain from splitting
    !                            up the interval at this point and
    !                            calling appopriate integrators on the
    !                            subranges.
    !                        = 9 the extrapolation table constructed for
    !                            convergence acceleration of the series
    !                            formed by the integral contributions
    !                            over the cycles, does not converge to
    !                            within the required accuracy.
    !                            as in the case of ifail = 8, it is advised
    !                            to examine the array iwork which contains
    !                            the error flags on the cycles.
    !                   if omega = 0 and integr = 1,
    !                   the integral is calculated by means of qagi
    !                   and ifail = ierlst(1) (with meaning as described
    !                   for ierlst(k), k = 1).
    !
    !           rslst  - real
    !                    vector of dimension at least limlst
    !                    rslst(k) contains the integral contribution
    !                    over the interval (a+(k-1)c,a+kc) where
    !                    c = (2*int(abs(omega))+1)*pi/abs(omega),
    !                    k = 1, 2, ..., lst.
    !                    note that, if omega = 0, rslst(1) contains
    !                    the value of the integral over (a,infinity).
    !
    !           erlst  - real
    !                    vector of dimension at least limlst
    !                    erlst(k) contains the error estimate
    !                    corresponding with rslst(k).
    !
    !           ierlst - integer
    !                    vector of dimension at least limlst
    !                    ierlst(k) contains the error flag corresponding
    !                    with rslst(k). for the meaning of the local error
    !                    flags see description of output parameter ifail.
    !
    !           lst    - integer
    !                    number of subintervals needed for the integration
    !                    if omega = 0 then lst is set to 1.
    !
    !           alist, blist, rlist, elist - real
    !                    vector of dimension at least limit,
    !
    !           iord, nnlog - integer
    !                    vector of dimension at least limit, providing
    !                    space for the quantities needed in the
    !                    subdivision process of each cycle
    !
    !           chebmo - real
    !                    array of dimension at least (maxp1,25),
    !                    providing space for the Chebyshev moments
    !                    needed within the cycles
    !
    ! Local parameters:
    !
    !          c1, c2    - end points of subinterval (of length
    !                      cycle)
    !          cycle     - (2*int(abs(omega))+1)*pi/abs(omega)
    !          psum      - vector of dimension at least (limexp+2)
    !                      (see routine qextr)
    !                      psum contains the part of the epsilon table
    !                      which is still needed for further computations.
    !                      each element of psum is a partial sum of
    !                      the series which should sum to the value of
    !                      the integral.
    !          errsum    - sum of error estimates over the
    !                      subintervals, calculated cumulatively
    !          epsa      - absolute tolerance requested over current
    !                      subinterval
    !          chebmo    - array containing the modified Chebyshev
    !                      moments (see also routine qc25o)

    !  REAL(dp), external    :: f
    interface
       function f(x)
         use kinds
         implicit none
         real(dp), intent(in) :: x
         real(dp)             :: f
       end function f
    end interface
    integer,  intent(in)  :: integr,limlst,limit,maxp1
    real(dp), intent(in)  :: a,omega,epsabs
    integer,  intent(out) :: neval,ifail,iord(limit), &
         ierlst(limit),lst,nnlog(limit)
    real(dp), intent(out) :: result,abserr, rslst(limlst),erlst(limit), &
         alist(limit),blist(limit),elist(limit),rlist(limit),chebmo(maxp1,25)

    ! Declare local variables
    real(dp), parameter :: pi = 3.14159265358979323846264338327950288419716939937510_dp
    real(dp), parameter :: p = 0.9_dp
    integer  :: ktmin,l,ll,momcom,nev,nres,numrl2
    real(dp) :: abseps,correc,cycle,c1,c2,dl,dla,drl,ep,eps,epsa,errsum,fact, &
         p1,psum(52),reseps,res3la(3)

    ! The dimension of  psum  is determined by the value of
    ! limexp in QEXTR (psum must be
    ! of dimension (limexp+2) at least).

    ! Test on validity of parameters.
    result = 0.0_dp
    abserr = 0.0_dp
    neval  = 0
    lst    = 0
    ifail  = 0

    if ( (integr /= 1 .and. integr /= 2 ) .or. &
         epsabs <= 0.0_dp .or. limlst < 3 ) then
       ifail = 6
       return
    endif

    if ( omega == 0.0_dp ) then
       if ( integr == 1 ) then
          call qagi (f,0.0_dp,1,epsabs,0.0_dp,result,abserr,neval,ifail)
       else
          result = 0.0_dp
          abserr = 0.0_dp
          neval = 0
          ifail = 0
       endif
       rslst(1) = result
       erlst(1) = abserr
       ierlst(1) = ifail
       lst = 1
       return
    endif

    ! Initializations.
    l      = abs ( omega )
    dl     = 2 * l + 1
    cycle  = dl * pi / abs ( omega )
    ifail  = 0
    ktmin  = 0
    neval  = 0
    numrl2 = 0
    nres   = 0
    c1     = a
    c2     = cycle+a
    p1     = 1.0_dp-p
    eps    = epsabs

    if ( epsabs > tiny ( epsabs ) / p1 ) eps = epsabs * p1

    ep = eps
    fact = 1.0_dp
    correc = 0.0_dp
    abserr = 0.0_dp
    errsum = 0.0_dp

    do lst = 1, limlst

       ! Integrate over current subinterval.
       dla = lst
       epsa = eps*fact
       call qfour(f,c1,c2,omega,integr,epsa,0.0_dp,limit,lst,maxp1, &
            rslst(lst),erlst(lst),nev,ierlst(lst),alist,blist,rlist,elist, &
            iord,nnlog,momcom,chebmo )
       neval = neval + nev
       fact = fact * p
       errsum = errsum + erlst(lst)
       drl = 5.0e+01_dp * abs(rslst(lst))

       ! Test on accuracy with partial sum.
       if ((errsum+drl) <= epsabs.and.lst >= 6) go to 80

       correc = max ( correc,erlst(lst))

       if ( ierlst(lst) /= 0 ) then
          eps = max ( ep,correc*p1)
          ifail = 7
       endif

       if ( ifail == 7 .and. (errsum+drl) <= correc*1.0e+01_dp.and. lst > 5) go to 80

       numrl2 = numrl2+1

       if ( lst <= 1 ) then
          psum(1) = rslst(1)
          go to 40
       endif

       psum(numrl2) = psum(ll)+rslst(lst)

       if ( lst == 2 ) go to 40

       ! Test on maximum number of subintervals
       if ( lst == limlst ) ifail = 8

       ! Perform new extrapolation
       call qextr(numrl2,psum,reseps,abseps,res3la,nres)

       ! Test whether extrapolated result is influenced by roundoff
       ktmin = ktmin+1

       if ( ktmin >= 15 .and. abserr <= 1.0e-03_dp * (errsum+drl) ) ifail = 9

       if ( abseps <= abserr .or. lst == 3 ) then
          abserr = abseps
          result = reseps
          ktmin = 0

          ! If IFAIL is not 0, check whether direct result (partial
          ! sum) or extrapolated result yields the best integral
          ! approximation
          if ( ( abserr + 1.0e+01_dp * correc ) <= epsabs ) exit
          if ( abserr <= epsabs .and. 1.0e+01_dp * correc >= epsabs ) exit
       endif

       if ( ifail /= 0 .and. ifail /= 7 ) exit

40     continue

       ll = numrl2
       c1 = c2
       c2 = c2+cycle
    enddo

    ! Set final result and error estimate.
60  continue

    abserr = abserr + 1.0e+01_dp * correc

    if ( ifail == 0 ) return

    if ( result /= 0.0_dp .and. psum(numrl2) /= 0.0_dp) go to 70

    if ( abserr > errsum ) go to 80

    if ( psum(numrl2) == 0.0_dp ) return

70  continue

    if ( abserr / abs(result) <= (errsum+drl)/abs(psum(numrl2)) ) then

       if ( ifail >= 1 .and. ifail /= 7 ) abserr = abserr + drl
       return
    endif

80  continue

    result = psum(numrl2)
    abserr = errsum + drl

  end subroutine qawfe

  !===============================================================================
  !###############################################################################
  !===============================================================================
  subroutine qawo(f,a,b,omega,integr,epsabs,epsrel,result,abserr,neval,ifail)

    !******************************************************************************
    !
    !! QAWO computes the integrals of oscillatory integrands.
    !
    !
    ! Discussion:
    !
    !   The routine calculates an approximation RESULT to a given
    !   definite integral
    !     I = Integral ( A <= X <= B ) F(X) * cos ( OMEGA * X ) dx
    !   or
    !     I = Integral ( A <= X <= B ) F(X) * sin ( OMEGA * X ) dx
    !   hopefully satisfying following claim for accuracy
    !     | I - RESULT | <= max ( epsabs, epsrel * |I| ).
    !
    ! Reference:
    !
    !   R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
    !   QUADPACK, a Subroutine Package for Automatic Integration,
    !   Springer Verlag, 1983
    !
    ! Parameters:
    !
    !   Input, external real F, the name of the function routine, of the form
    !     function f ( x )
    !     real f
    !     real x
    !   which evaluates the integrand function.
    !
    !   Input, real A, B, the limits of integration.
    !
    !   Input, real OMEGA, the parameter in the weight function.
    !
    !   Input, integer INTEGR, specifies the weight function:
    !   1, W(X) = cos ( OMEGA * X )
    !   2, W(X) = sin ( OMEGA * X )
    !
    !   Input, real EPSABS, EPSREL, the absolute and relative accuracy requested.
    !
    !   Output, real RESULT, the estimated value of the integral.
    !
    !   Output, real ABSERR, an estimate of || I - RESULT ||.
    !
    !   Output, integer NEVAL, the number of times the integral was evaluated.
    !
    !         ifail    - integer
    !                  ifail = 0 normal and reliable termination of the
    !                            routine. it is assumed that the
    !                            requested accuracy has been achieved.
    !                - ifail > 0 abnormal termination of the routine.
    !                            the estimates for integral and error are
    !                            less reliable. it is assumed that the
    !                            requested accuracy has not been achieved.
    !                  ifail = 1 maximum number of subdivisions allowed
    !                            (= leniw/2) has been achieved. one can
    !                            allow more subdivisions by increasing the
    !                            value of leniw (and taking the according
    !                            dimension adjustments into account).
    !                            however, if this yields no improvement it
    !                            is advised to analyze the integrand in
    !                            order to determine the integration
    !                            difficulties. if the position of a local
    !                            difficulty can be determined (e.g.
    !                            singularity, discontinuity within the
    !                            interval) one will probably gain from
    !                            splitting up the interval at this point
    !                            and calling the integrator on the
    !                            subranges. if possible, an appropriate
    !                            special-purpose integrator should
    !                            be used which is designed for handling
    !                            the type of difficulty involved.
    !                        = 2 the occurrence of roundoff error is
    !                            detected, which prevents the requested
    !                            tolerance from being achieved.
    !                            the error may be under-estimated.
    !                        = 3 extremely bad integrand behavior occurs
    !                            at some interior points of the integration
    !                            interval.
    !                        = 4 the algorithm does not converge. roundoff
    !                            error is detected in the extrapolation
    !                            table. it is presumed that the requested
    !                            tolerance cannot be achieved due to
    !                            roundoff in the extrapolation table,
    !                            and that the returned result is the best
    !                            which can be obtained.
    !                        = 5 the integral is probably divergent, or
    !                            slowly convergent. it must be noted that
    !                            divergence can occur with any other value
    !                            of ifail.
    !                        = 6 the input is invalid, because
    !                            epsabs < 0 and epsrel < 0,
    !                            result, abserr, neval are set to zero.
    !
    ! Local parameters:
    ! integer :: ktmin,l,ll,momcom,nev,nres,numrl2
    !   limit is the maximum number of subintervals allowed in the
    !   subdivision process of QFOUR. take care that limit >= 1.
    !
    !   maxp1 gives an upper bound on the number of Chebyshev moments
    !   which can be stored, i.e. for the intervals of lengths
    !   abs(b-a)*2**(-l), l = 0, 1, ... , maxp1-2. take care that
    !   maxp1 >= 1.

    !  REAL(dp), external    :: f
    interface
       function f(x)
         use kinds
         implicit none
         real(dp), intent(in) :: x
         real(dp)             :: f
       end function f
    end interface
    integer,  intent(in)  :: integr
    real(dp), intent(in)  :: a,b,omega,epsabs,epsrel
    integer,  intent(out) :: neval,ifail
    real(dp), intent(out) :: result,abserr

    ! Declare local variables
    integer, parameter :: limit = 500,maxp1 = 21
    integer  :: iord(limit),momcom,nnlog(limit)
    real(dp) :: alist(limit),blist(limit),chebmo(maxp1,25),elist(limit),rlist(limit)

    call qfour(f,a,b,omega,integr,epsabs,epsrel,limit,1,maxp1,result, &
         abserr,neval,ifail,alist,blist,rlist,elist,iord,nnlog,momcom,chebmo)

  end subroutine qawo

  !===============================================================================
  !###############################################################################
  !===============================================================================
  subroutine qaws(f,a,b,alfa,beta,integr,epsabs,epsrel,result,abserr,neval,ifail)
    !
    !******************************************************************************
    !
    !! QAWS estimates integrals with algebraico-logarithmic endpoint singularities.
    !
    !
    ! Discussion:
    !
    !   This routine calculates an approximation RESULT to a given
    !   definite integral
    !     I = integral of f*w over (a,b)
    !   where w shows a singular behavior at the end points, see parameter
    !   integr, hopefully satisfying following claim for accuracy
    !     abs(i-result) <= max(epsabs,epsrel*abs(i)).
    !
    ! Reference:
    !
    !   R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
    !   QUADPACK, a Subroutine Package for Automatic Integration,
    !   Springer Verlag, 1983
    !
    ! Parameters:
    !
    !   Input, external real F, the name of the function routine, of the form
    !     function f ( x )
    !     real f
    !     real x
    !   which evaluates the integrand function.
    !
    !   Input, real A, B, the limits of integration.
    !
    !   Input, real ALFA, BETA, parameters used in the weight function.
    !   ALFA and BETA should be greater than -1.
    !
    !   Input, integer INTEGR, indicates which weight function is to be used
    !   = 1  (x-a)**alfa*(b-x)**beta
    !   = 2  (x-a)**alfa*(b-x)**beta*log(x-a)
    !   = 3  (x-a)**alfa*(b-x)**beta*log(b-x)
    !   = 4  (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)
    !
    !   Input, real EPSABS, EPSREL, the absolute and relative accuracy requested.
    !
    !   Output, real RESULT, the estimated value of the integral.
    !
    !   Output, real ABSERR, an estimate of || I - RESULT ||.
    !
    !   Output, integer NEVAL, the number of times the integral was evaluated.
    !
    !         ifail    - integer
    !                  ifail = 0 normal and reliable termination of the
    !                            routine. it is assumed that the requested
    !                            accuracy has been achieved.
    !                  ifail > 0 abnormal termination of the routine
    !                            the estimates for the integral and error
    !                            are less reliable. it is assumed that the
    !                            requested accuracy has not been achieved.
    !                  ifail = 1 maximum number of subdivisions allowed
    !                            has been achieved. one can allow more
    !                            subdivisions by increasing the data value
    !                            of limit in qaws (and taking the according
    !                            dimension adjustments into account).
    !                            however, if this yields no improvement it
    !                            is advised to analyze the integrand, in
    !                            order to determine the integration
    !                            difficulties which prevent the requested
    !                            tolerance from being achieved. in case of
    !                            a jump discontinuity or a local
    !                            singularity of algebraico-logarithmic type
    !                            at one or more interior points of the
    !                            integration range, one should proceed by
    !                            splitting up the interval at these points
    !                            and calling the integrator on the
    !                            subranges.
    !                        = 2 the occurrence of roundoff error is
    !                            detected, which prevents the requested
    !                            tolerance from being achieved.
    !                        = 3 extremely bad integrand behavior occurs
    !                            at some points of the integration
    !                            interval.
    !                        = 6 the input is invalid, because
    !                            b <= a or alfa <= (-1) or beta <= (-1) or
    !                            integr < 1 or integr > 4 or
    !                            epsabs < 0 and epsrel < 0,
    !                            result, abserr, neval are set to zero.
    !
    ! Local parameters:
    !
    !   LIMIT is the maximum number of subintervals allowed in the
    !   subdivision process of qawse. take care that limit >= 2.

    !  REAL(dp), external    :: f
    interface
       function f(x)
         use kinds
         implicit none
         real(dp), intent(in) :: x
         real(dp)             :: f
       end function f
    end interface
    integer,  intent(in)  :: integr
    real(dp), intent(in)  :: a,b,alfa,beta,epsabs,epsrel
    integer,  intent(out) :: neval,ifail
    real(dp), intent(out) :: result,abserr

    ! Declare local variables
    integer, parameter :: limit = 500
    integer  :: iord(limit),last
    real(dp) :: alist(limit),blist(limit),elist(limit),rlist(limit)

    call qawse(f,a,b,alfa,beta,integr,epsabs,epsrel,limit,result, &
         abserr,neval,ifail,alist,blist,rlist,elist,iord,last)

  end subroutine qaws

  !===============================================================================
  !###############################################################################
  !===============================================================================
  subroutine qawse(f,a,b,alfa,beta,integr,epsabs,epsrel,limit, &
       result,abserr,neval,ifail,alist,blist,rlist,elist,iord,last)
    !
    !******************************************************************************
    !
    !! QAWSE estimates integrals with algebraico-logarithmic endpoint singularities.
    !
    !
    ! Discussion:
    !
    !   This routine calculates an approximation RESULT to an integral
    !     I = integral of F(X) * W(X) over (a,b),
    !   where W(X) shows a singular behavior at the endpoints, hopefully
    !   satisfying:
    !     | I - RESULT | <= max ( epsabs, epsrel * |I| ).
    !
    ! Reference:
    !
    !   R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
    !   QUADPACK, a Subroutine Package for Automatic Integration,
    !   Springer Verlag, 1983
    !
    ! Parameters:
    !
    !   Input, external real F, the name of the function routine, of the form
    !     function f ( x )
    !     real f
    !     real x
    !   which evaluates the integrand function.
    !
    !   Input, real A, B, the limits of integration.
    !
    !   Input, real ALFA, BETA, parameters used in the weight function.
    !   ALFA and BETA should be greater than -1.
    !
    !   Input, integer INTEGR, indicates which weight function is used:
    !   = 1  (x-a)**alfa*(b-x)**beta
    !   = 2  (x-a)**alfa*(b-x)**beta*log(x-a)
    !   = 3  (x-a)**alfa*(b-x)**beta*log(b-x)
    !   = 4  (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)
    !
    !   Input, real EPSABS, EPSREL, the absolute and relative accuracy requested.
    !
    !   Input, integer LIMIT, an upper bound on the number of subintervals
    !   in the partition of (A,B), LIMIT >= 2.  If LIMIT < 2, the routine
    !    will end with IFAIL = 6.
    !
    !   Output, real RESULT, the estimated value of the integral.
    !
    !   Output, real ABSERR, an estimate of || I - RESULT ||.
    !
    !   Output, integer NEVAL, the number of times the integral was evaluated.
    !
    !         ifail    - integer
    !                  ifail = 0 normal and reliable termination of the
    !                            routine. it is assumed that the requested
    !                            accuracy has been achieved.
    !                  ifail > 0 abnormal termination of the routine
    !                            the estimates for the integral and error
    !                            are less reliable. it is assumed that the
    !                            requested accuracy has not been achieved.
    !                        = 1 maximum number of subdivisions allowed
    !                            has been achieved. one can allow more
    !                            subdivisions by increasing the value of
    !                            limit. however, if this yields no
    !                            improvement it is advised to analyze the
    !                            integrand, in order to determine the
    !                            integration difficulties which prevent
    !                            the requested tolerance from being
    !                            achieved. in case of a jump discontinuity
    !                            or a local singularity of algebraico-
    !                            logarithmic type at one or more interior
    !                            points of the integration range, one
    !                            should proceed by splitting up the
    !                            interval at these points and calling the
    !                            integrator on the subranges.
    !                        = 2 the occurrence of roundoff error is
    !                            detected, which prevents the requested
    !                            tolerance from being achieved.
    !                        = 3 extremely bad integrand behavior occurs
    !                            at some points of the integration
    !                            interval.
    !                        = 6 the input is invalid, because
    !                            b <= a or alfa <= (-1) or beta <= (-1) or
    !                            integr < 1 or integr > 4, or
    !                            epsabs < 0 and epsrel < 0,
    !                            or limit < 2.
    !                            result, abserr, neval, rlist(1), elist(1),
    !                            iord(1) and last are set to zero.
    !                            alist(1) and blist(1) are set to a and b
    !                            respectively.
    !
    !   Workspace, real ALIST(LIMIT), BLIST(LIMIT), contains in entries 1
    !   through LAST the left and right ends of the partition subintervals.
    !
    !   Workspace, real RLIST(LIMIT), contains in entries 1 through LAST
    !   the integral approximations on the subintervals.
    !
    !   Workspace, real ELIST(LIMIT), contains in entries 1 through LAST
    !   the absolute error estimates on the subintervals.
    !
    !           iord   - integer
    !                    vector of dimension at least limit, the first k
    !                    elements of which are pointers to the error
    !                    estimates over the subintervals, so that
    !                    elist(iord(1)), ..., elist(iord(k)) with k = last
    !                    if last <= (limit/2+2), and k = limit+1-last
    !                    otherwise, form a decreasing sequence.
    !
    !   Output, integer LAST, the number of subintervals actually produced in
    !   the subdivision process.
    !
    ! Local parameters:
    !
    !          alist     - list of left end points of all subintervals
    !                      considered up to now
    !          blist     - list of right end points of all subintervals
    !                      considered up to now
    !          rlist(i)  - approximation to the integral over
    !                      (alist(i),blist(i))
    !          elist(i)  - error estimate applying to rlist(i)
    !          maxerr    - pointer to the interval with largest error
    !                      estimate
    !          errmax    - elist(maxerr)
    !          area      - sum of the integrals over the subintervals
    !          errsum    - sum of the errors over the subintervals
    !          errbnd    - requested accuracy max(epsabs,epsrel*
    !                      abs(result))
    !          *****1    - variable for the left subinterval
    !          *****2    - variable for the right subinterval
    !          last      - index for subdivision
    !

    !  REAL(dp), external    :: f
    interface
       function f(x)
         use kinds
         implicit none
         real(dp), intent(in) :: x
         real(dp)             :: f
       end function f
    end interface
    integer,  intent(in)  :: integr,limit
    real(dp), intent(in)  :: a,b,alfa,beta,epsabs,epsrel
    integer,  intent(out) :: neval,ifail,iord(limit)
    real(dp), intent(out) :: result,abserr, &
         alist(limit),blist(limit),elist(limit),rlist(limit)

    ! Declare local variables
    integer  :: iroff1,iroff2,last,maxerr,nev,nrmax
    real(dp) :: area,area1,area12,area2,a1,a2,b1,b2,centre,errbnd,errmax,error1, &
         erro12,error2,errsum,resas1,resas2,rg(25),rh(25),ri(25),rj(25)

    ! Test on validity of parameters.
    ifail    = 0
    neval    = 0
    last     = 0
    rlist(1) = 0.0_dp
    elist(1) = 0.0_dp
    iord(1)  = 0
    result   = 0.0_dp
    abserr   = 0.0_dp

    if ( b <= a .or. &
         (epsabs < 0.0_dp .and. epsrel < 0.0_dp)  .or. &
         alfa <= (-1.0_dp) .or. beta <= (-1.0_dp) .or. &
         integr < 1 .or. integr > 4 .or. limit < 2) then
       ifail = 6
       return
    endif

    ! Compute the modified Chebyshev moments.
    call qmomo(alfa,beta,ri,rj,rg,rh,integr)

    ! Integrate over the intervals (a,(a+b)/2) and ((a+b)/2,b).
    centre = 5.0e-01_dp * ( b + a )

    call qc25s(f,a,b,a,centre,alfa,beta,ri,rj,rg,rh,area1,error1,resas1,integr,nev)
    neval = nev

    call qc25s(f,a,b,centre,b,alfa,beta,ri,rj,rg,rh,area2,error2,resas2,integr,nev)

    last   = 2
    neval  = neval+nev
    result = area1+area2
    abserr = error1+error2

    ! Test on accuracy.
    errbnd = max ( epsabs, epsrel * abs ( result ) )

    ! Initialization.
    if ( error2 <= error1 ) then
       alist(1) = a
       alist(2) = centre
       blist(1) = centre
       blist(2) = b
       rlist(1) = area1
       rlist(2) = area2
       elist(1) = error1
       elist(2) = error2
    else
       alist(1) = centre
       alist(2) = a
       blist(1) = b
       blist(2) = centre
       rlist(1) = area2
       rlist(2) = area1
       elist(1) = error2
       elist(2) = error1
    endif

    iord(1) = 1
    iord(2) = 2

    if ( limit == 2 ) then
       ifail = 1
       return
    endif

    if ( abserr <= errbnd ) return

    errmax = elist(1)
    maxerr = 1
    nrmax  = 1
    area   = result
    errsum = abserr
    iroff1 = 0
    iroff2 = 0

    do last = 3, limit

       ! Bisect the subinterval with largest error estimate.
       a1 = alist(maxerr)
       b1 = 5.0e-01_dp * (alist(maxerr)+blist(maxerr))
       a2 = b1
       b2 = blist(maxerr)

       call qc25s(f,a,b,a1,b1,alfa,beta,ri,rj,rg,rh,area1,error1,resas1,integr,nev)
       neval = neval + nev
       call qc25s(f,a,b,a2,b2,alfa,beta,ri,rj,rg,rh,area2,error2,resas2,integr,nev)
       neval = neval + nev

       ! Improve previous approximations integral and error and
       ! test for accuracy.
       area12 = area1+area2
       erro12 = error1+error2
       errsum = errsum+erro12-errmax
       area = area+area12-rlist(maxerr)

       ! Test for roundoff error.
       if ( a /= a1 .and. b /= b2 ) then
          if ( resas1 /= error1 .and. resas2 /= error2 ) then

             if ( abs(rlist(maxerr)-area12) < 1.0e-05*abs(area12) &
                  .and.erro12 >= 9.9e-01_dp*errmax) iroff1 = iroff1+1
             if ( last > 10.and.erro12 > errmax ) iroff2 = iroff2+1

          endif
       endif

       rlist(maxerr) = area1
       rlist(last)   = area2

       ! Test on accuracy.
       errbnd = max ( epsabs, epsrel * abs ( area ) )

       if ( errsum > errbnd ) then

          ! Set error flag in the case that the number of interval
          ! bisections exceeds limit.
          if ( last == limit ) ifail = 1

          ! Set error flag in the case of roundoff error.
          if ( iroff1 >= 6 .or. iroff2 >= 20 ) ifail = 2

          ! Set error flag in the case of bad integrand behavior
          ! at interior points of integration range.
          if ( max ( abs(a1),abs(b2)) <= (1.0_dp+1.0e+03_dp* epsilon ( a1 ) )* &
               (abs(a2)+1.0e+03_dp* tiny ( a2) )) ifail = 3
       endif

       ! Append the newly-created intervals to the list.
       if ( error2 <= error1 ) then
          alist(last) = a2
          blist(maxerr) = b1
          blist(last) = b2
          elist(maxerr) = error1
          elist(last) = error2
       else
          alist(maxerr) = a2
          alist(last) = a1
          blist(last) = b1
          rlist(maxerr) = area2
          rlist(last) = area1
          elist(maxerr) = error2
          elist(last) = error1
       endif

! Call QSORT to maintain the descending ordering
! in the list of error estimates and select the subinterval
! with largest error estimate (to be bisected next).
    call qsort(limit,last,maxerr,errmax,elist,iord,nrmax)

    if (ifail /= 0 .or. errsum <= errbnd ) exit
 enddo

! Compute final result.
 result = sum ( rlist(1:last) )
 abserr = errsum

end subroutine qawse

!===============================================================================
!###############################################################################
!===============================================================================
subroutine qc25c(f,a,b,c,result,abserr,krul,neval)
!
!******************************************************************************
!
!! QC25C returns integration rules for Cauchy Principal Value integrals.
!
!
! Discussion:
!
!   This routine estimates
!     I = integral of F(X) * W(X) over (a,b)
!   with error estimate, where
!     w(x) = 1/(x-c)
!
! Reference:
!
!   R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
!   QUADPACK, a Subroutine Package for Automatic Integration,
!   Springer Verlag, 1983
!
! Parameters:
!
!   Input, external real F, the name of the function routine, of the form
!     function f ( x )
!     real f
!     real x
!   which evaluates the integrand function.
!
!   Input, real A, B, the limits of integration.
!
!   Input, real C, the parameter in the weight function.
!
!   Output, real RESULT, the estimated value of the integral.
!   RESULT is computed by using a generalized Clenshaw-Curtis method if
!   C lies within ten percent of the integration interval.  In the
!   other case the 15-point Kronrod rule obtained by optimal addition
!   of abscissae to the 7-point Gauss rule, is applied.
!
!   Output, real ABSERR, an estimate of || I - RESULT ||.
!
!          krul   - integer
!                   key which is decreased by 1 if the 15-point
!                   Gauss-Kronrod scheme has been used
!
!   Output, integer NEVAL, the number of times the integral was evaluated.
!
! Local parameters:
!
!          fval   - value of the function f at the points
!                   cos(k*pi/24),  k = 0, ..., 24
!          cheb12 - Chebyshev series expansion coefficients, for the
!                   function f, of degree 12
!          cheb24 - Chebyshev series expansion coefficients, for the
!                   function f, of degree 24
!          res12  - approximation to the integral corresponding to the
!                   use of cheb12
!          res24  - approximation to the integral corresponding to the
!                   use of cheb24
!          qwgtc  - external function subprogram defining the weight
!                   function
!          hlgth  - half-length of the interval
!          centr  - mid point of the interval

!  REAL(dp), external    :: f
  interface
    function f(x)
      use kinds
      implicit none
      real(dp), intent(in) :: x
      real(dp)             :: f
    end function f
  end interface
  real(dp), intent(in)  :: a,b,c
  integer,  intent(out) :: neval,krul
  real(dp), intent(out) :: result,abserr

! Declare local variables
  integer  :: i,isym,k,kp
  real(dp) :: ak22,amom0,amom1,amom2,cc,centr,cheb12(13),cheb24(25),fval(25),hlgth, &
              p2,p3,p4,resabs,resasc,res12,res24,u

  real(dp), parameter, dimension ( 11 ) :: x = (/ &
    9.914448613738104e-01_dp, 9.659258262890683e-01_dp, &
    9.238795325112868e-01_dp, 8.660254037844386e-01_dp, &
    7.933533402912352e-01_dp, 7.071067811865475e-01_dp, &
    6.087614290087206e-01_dp, 5.000000000000000e-01_dp, &
    3.826834323650898e-01_dp, 2.588190451025208e-01_dp, &
    1.305261922200516e-01_dp /)

! Check the position of C.
  cc = ( 2.0_dp * c - b - a ) / ( b - a )

! Apply the 15-point Gauss-Kronrod scheme.
  if ( abs ( cc ) >= 1.1_dp ) then
     krul = krul - 1
     call qk15w (f,qwgtc,c,p2,p3,p4,kp,a,b,result,abserr,resabs,resasc)
     neval = 15
     if ( resasc == abserr ) krul = krul+1
     return
  endif

! Use the generalized Clenshaw-Curtis method.
  hlgth = 5.0e-01_dp * ( b - a )
  centr = 5.0e-01_dp * ( b + a )
  neval = 25
  fval(1) = 5.0e-01_dp * f(hlgth+centr)
  fval(13) = f(centr)
  fval(25) = 5.0e-01_dp * f(centr-hlgth)

  do i = 2, 12
     u = hlgth * x(i-1)
     isym = 26 - i
     fval(i) = f(u+centr)
     fval(isym) = f(centr-u)
  enddo

! Compute the Chebyshev series expansion.
  call qcheb(x,fval,cheb12,cheb24)

! The modified Chebyshev moments are computed by forward
! recursion, using AMOM0 and AMOM1 as starting values.
  amom0 = log ( abs ( ( 1.0_dp - cc ) / ( 1.0_dp + cc ) ) )
  amom1 = 2.0_dp + cc * amom0
  res12 = cheb12(1) * amom0 + cheb12(2) * amom1
  res24 = cheb24(1) * amom0 + cheb24(2) * amom1

  do k = 3, 13
     amom2 = 2.0_dp * cc * amom1 - amom0
     ak22 = ( k - 2 ) * ( k - 2 )
     if ( ( k / 2 ) * 2 == k ) amom2 = amom2 - 4.0_dp / ( ak22 - 1.0_dp )
     res12 = res12 + cheb12(k) * amom2
     res24 = res24 + cheb24(k) * amom2
     amom0 = amom1
     amom1 = amom2
  enddo

  do k = 14, 25
     amom2 = 2.0_dp * cc * amom1 - amom0
     ak22 = ( k - 2 ) * ( k - 2 )
     if ( ( k / 2 ) * 2 == k ) amom2 = amom2 - 4.0_dp / ( ak22 - 1.0_dp )
     res24 = res24 + cheb24(k) * amom2
     amom0 = amom1
     amom1 = amom2
  enddo

  result = res24
  abserr = abs ( res24 - res12 )

end subroutine qc25c

!===============================================================================
!###############################################################################
!===============================================================================
subroutine qc25o(f,a,b,omega,integr,nrmom,maxp1,ksave,result, &
                  abserr,neval,resabs,resasc,momcom,chebmo )
!
!******************************************************************************
!
!! QC25O returns integration rules for integrands with a COS or SIN factor.
!
!
! Discussion:
!
!   This routine estimates the integral
!     I = integral of f(x) * w(x) over (a,b)
!   where
!     w(x) = cos(omega*x)
!   or
!     w(x) = sin(omega*x),
!   and estimates
!     J = integral ( A <= X <= B ) |F(X)| dx.
!
!   For small values of OMEGA or small intervals (a,b) the 15-point
!   Gauss-Kronrod rule is used.  In all other cases a generalized
!   Clenshaw-Curtis method is used, that is, a truncated Chebyshev
!   expansion of the function F is computed on (a,b), so that the
!   integrand can be written as a sum of terms of the form W(X)*T(K,X),
!   where T(K,X) is the Chebyshev polynomial of degree K.  The Chebyshev
!   moments are computed with use of a linear recurrence relation.
!
! Reference:
!
!   R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
!   QUADPACK, a Subroutine Package for Automatic Integration,
!   Springer Verlag, 1983
!
! Parameters:
!
!   Input, external real F, the name of the function routine, of the form
!     function f ( x )
!     real f
!     real x
!   which evaluates the integrand function.
!
!   Input, real A, B, the limits of integration.
!
!   Input, real OMEGA, the parameter in the weight function.
!
!   Input, integer INTEGR, indicates which weight function is to be used
!   = 1, w(x) = cos(omega*x)
!   = 2, w(x) = sin(omega*x)
!
!   ?, integer NRMOM, the length of interval (a,b) is equal to the length
!   of the original integration interval divided by
!   2**nrmom (we suppose that the routine is used in an
!   adaptive integration process, otherwise set
!   nrmom = 0).  nrmom must be zero at the first call.
!
!          maxp1  - integer
!                   gives an upper bound on the number of Chebyshev
!                   moments which can be stored, i.e. for the intervals
!                   of lengths abs(bb-aa)*2**(-l), l = 0,1,2, ...,
!                   maxp1-2.
!
!          ksave  - integer
!                   key which is one when the moments for the
!                   current interval have been computed
!
!   Output, real RESULT, the estimated value of the integral.
!
!          abserr - real
!                   estimate of the modulus of the absolute
!                   error, which should equal or exceed abs(i-result)
!
!   Output, integer NEVAL, the number of times the integral was evaluated.
!
!   Output, real RESABS, approximation to the integral J.
!
!   Output, real RESASC, approximation to the integral of abs(F-I/(B-A)).
!
!        on entry and return
!          momcom - integer
!                   for each interval length we need to compute
!                   the Chebyshev moments. momcom counts the number
!                   of intervals for which these moments have already
!                   been computed. if nrmom < momcom or ksave = 1,
!                   the Chebyshev moments for the interval (a,b)
!                   have already been computed and stored, otherwise
!                   we compute them and we increase momcom.
!
!          chebmo - real
!                   array of dimension at least (maxp1,25) containing
!                   the modified Chebyshev moments for the first momcom
!                   interval lengths
!
! Local parameters:
!
!   maxp1 gives an upper bound
!          on the number of Chebyshev moments which can be
!          computed, i.e. for the interval (bb-aa), ...,
!          (bb-aa)/2**(maxp1-2).
!          should this number be altered, the first dimension of
!          chebmo needs to be adapted.
!
!   x contains the values cos(k*pi/24)
!          k = 1, ...,11, to be used for the Chebyshev expansion of f
!
!          centr  - mid point of the integration interval
!          hlgth  - half length of the integration interval
!          fval   - value of the function f at the points
!                   (b-a)*0.5*cos(k*pi/12) + (b+a)*0.5
!                   k = 0, ...,24
!          cheb12 - coefficients of the Chebyshev series expansion
!                   of degree 12, for the function f, in the
!                   interval (a,b)
!          cheb24 - coefficients of the Chebyshev series expansion
!                   of degree 24, for the function f, in the
!                   interval (a,b)
!          resc12 - approximation to the integral of
!                   cos(0.5*(b-a)*omega*x)*f(0.5*(b-a)*x+0.5*(b+a))
!                   over (-1,+1), using the Chebyshev series
!                   expansion of degree 12
!          resc24 - approximation to the same integral, using the
!                   Chebyshev series expansion of degree 24
!          ress12 - the analogue of resc12 for the sine
!          ress24 - the analogue of resc24 for the sine
!

!  REAL(dp), external    :: f
  interface
    function f(x)
      use kinds
      implicit none
      real(dp), intent(in) :: x
      real(dp)             :: f
    end function f
  end interface
  integer,  intent(in)  :: integr,nrmom,maxp1,ksave
  real(dp), intent(in)  :: a,b,omega
  integer,  intent(out) :: neval,momcom
  real(dp), intent(out) :: result,abserr, resabs, resasc,chebmo(maxp1,25)

! Declare local variables
  integer, parameter :: nmac = 28
  integer  :: i,isym,j,k,m,noeq1,noequ
  real(dp) :: ac,an,an2,as,asap,ass,centr,cheb12(13),cheb24(25),conc,cons,cospar, &
              d(28),d1(28),d2(28),d3(28),estc,ests,fval(25),hlgth,parint,par2,par22, &
              p2,p3,p4,resc12,resc24,ress12,ress24,sinpar,v(28)

  real(dp), dimension ( 11 ) :: x = (/ &
    9.914448613738104e-01_dp,     9.659258262890683e-01_dp, &
    9.238795325112868e-01_dp,     8.660254037844386e-01_dp, &
    7.933533402912352e-01_dp,     7.071067811865475e-01_dp, &
    6.087614290087206e-01_dp,     5.000000000000000e-01_dp, &
    3.826834323650898e-01_dp,     2.588190451025208e-01_dp, &
    1.305261922200516e-01_dp /)

  centr  = 5.0e-01_dp*(b+a)
  hlgth  = 5.0e-01_dp*(b-a)
  parint = omega * hlgth

! Compute the integral using the 15-point Gauss-Kronrod
! formula if the value of the parameter in the integrand
! is small or if the length of the integration interval
! is less than (bb-aa)/2**(maxp1-2), where (aa,bb) is the
! original integration interval.
  if ( abs ( parint ) <= 2.0_dp ) then
     call qk15w(f,qwgto,omega,p2,p3,p4,integr,a,b,result,abserr,resabs,resasc)
     neval = 15
     return
  endif

! Compute the integral using the generalized clenshaw-curtis method.
  conc   = hlgth * cos(centr*omega)
  cons   = hlgth * sin(centr*omega)
  resasc = huge ( resasc )
  neval  = 25

! Check whether the Chebyshev moments for this interval
! have already been computed.
  if ( nrmom < momcom .or. ksave == 1 ) go to 140

! Compute a new set of Chebyshev moments.
  m      = momcom+1
  par2   = parint*parint
  par22  = par2+2.0_dp
  sinpar = sin(parint)
  cospar = cos(parint)

! Compute the Chebyshev moments with respect to cosine.
  v(1) = 2.0_dp*sinpar/parint
  v(2) = (8.0_dp*cospar+(par2+par2-8.0_dp)*sinpar/ parint)/par2
  v(3) = (3.2e+01_dp*(par2-1.2e+01_dp)*cospar+(2.0_dp* &
       ((par2-8.0e+01_dp)*par2+1.92e+02_dp)*sinpar)/ &
       parint)/(par2*par2)
  ac = 8.0_dp*cospar
  as = 2.4e+01_dp*parint*sinpar

  if ( abs ( parint ) > 2.4e+01_dp ) go to 70

! Compute the Chebyshev moments as the solutions of a boundary value
! problem with one initial value (v(3)) and one end value computed
! using an asymptotic formula.
  noequ = nmac-3
  noeq1 = noequ-1
  an    = 6.0_dp

  do k = 1, noeq1
     an2    = an*an
     d(k)   = -2.0_dp*(an2-4.0_dp)*(par22-an2-an2)
     d2(k)  = (an-1.0_dp)*(an-2.0_dp)*par2
     d1(k)  = (an+3.0_dp)*(an+4.0_dp)*par2
     v(k+3) = as-(an2-4.0_dp)*ac
     an     = an+2.0_dp
  enddo

  an2 = an*an
  d(noequ) = -2.0_dp*(an2-4.0_dp)*(par22-an2-an2)
  v(noequ+3) = as-(an2-4.0_dp)*ac
  v(4) = v(4)-5.6e+01_dp*par2*v(3)
  ass = parint*sinpar
  asap = (((((2.10e+02_dp*par2-1.0_dp)*cospar-(1.05e+02_dp*par2 &
       -6.3e+01_dp)*ass)/an2-(1.0_dp-1.5e+01_dp*par2)*cospar &
       +1.5e+01_dp*ass)/an2-cospar+3.0_dp*ass)/an2-cospar)/an2
  v(noequ+3) = v(noequ+3)-2.0_dp*asap*par2*(an-1.0_dp)* &
       (an-2.0_dp)

! Solve the tridiagonal system by means of Gaussian
! elimination with partial pivoting.
  d3(1:noequ) = 0.0_dp
  d2(noequ)   = 0.0_dp

  do i = 1, noeq1
     if ( abs(d1(i)) > abs(d(i)) ) then
        an = d1(i)
        d1(i) = d(i)
        d(i) = an
        an = d2(i)
        d2(i) = d(i+1)
        d(i+1) = an
        d3(i) = d2(i+1)
        d2(i+1) = 0.0_dp
        an = v(i+4)
        v(i+4) = v(i+3)
        v(i+3) = an
     endif
     d(i+1) = d(i+1)-d2(i)*d1(i)/d(i)
     d2(i+1) = d2(i+1)-d3(i)*d1(i)/d(i)
     v(i+4) = v(i+4)-v(i+3)*d1(i)/d(i)
  enddo

  v(noequ+3) = v(noequ+3)/d(noequ)
  v(noequ+2) = (v(noequ+2)-d2(noeq1)*v(noequ+3))/d(noeq1)

  do i = 2, noeq1
     k = noequ-i
     v(k+3) = (v(k+3)-d3(k)*v(k+5)-d2(k)*v(k+4))/d(k)
  enddo
  go to 90

! Compute the Chebyshev moments by means of forward recursion
70 continue
  an = 4.0_dp

  do i = 4, 13
     an2 = an*an
     v(i) = ((an2-4.0_dp)*(2.0_dp*(par22-an2-an2)*v(i-1)-ac) &
          +as-par2*(an+1.0_dp)*(an+2.0_dp)*v(i-2))/ &
          (par2*(an-1.0_dp)*(an-2.0_dp))
     an = an+2.0_dp
  enddo
90 continue

  do j = 1, 13
     chebmo(m,2*j-1) = v(j)
  enddo

! Compute the Chebyshev moments with respect to sine.
  v(1) = 2.0_dp*(sinpar-parint*cospar)/par2
  v(2) = (1.8e+01_dp-4.8e+01_dp/par2)*sinpar/par2 &
       +(-2.0_dp+4.8e+01_dp/par2)*cospar/parint
  ac = -2.4e+01_dp*parint*cospar
  as = -8.0_dp*sinpar
  chebmo(m,2) = v(1)
  chebmo(m,4) = v(2)

  if ( abs(parint) <= 2.4e+01_dp ) then

     do k = 3, 12
        an = k
        chebmo(m,2*k) = -sinpar/(an*(2.0_dp*an-2.0_dp)) &
             -2.5e-01_dp*parint*(v(k+1)/an-v(k)/(an-1.0_dp))
     enddo

! Compute the Chebyshev moments by means of forward recursion.
  else
     an = 3.0_dp
     do i = 3, 12
        an2 = an*an
        v(i) = ((an2-4.0_dp)*(2.0_dp*(par22-an2-an2)*v(i-1)+as) &
             +ac-par2*(an+1.0_dp)*(an+2.0_dp)*v(i-2)) &
             /(par2*(an-1.0_dp)*(an-2.0_dp))
        an = an+2.0_dp
        chebmo(m,2*i) = v(i)
    enddo
  endif

140 continue

  if ( nrmom < momcom ) m = nrmom + 1

  if ( momcom < maxp1 - 1 .and. nrmom >= momcom ) momcom = momcom + 1

! Compute the coefficients of the Chebyshev expansions
! of degrees 12 and 24 of the function F.
  fval(1)  = 5.0e-01_dp*f(centr+hlgth)
  fval(13) = f(centr)
  fval(25) = 5.0e-01_dp*f(centr-hlgth)

  do i = 2, 12
     isym       = 26-i
     fval(i)    = f(hlgth*x(i-1)+centr)
     fval(isym) = f(centr-hlgth*x(i-1))
  enddo

  call qcheb(x,fval,cheb12,cheb24)

! Compute the integral and error estimates.
  resc12 = cheb12(13) * chebmo(m,13)
  ress12 = 0.0_dp
  estc = abs ( cheb24(25)*chebmo(m,25))+abs((cheb12(13)- &
       cheb24(13))*chebmo(m,13) )
  ests = 0.0_dp
  k = 11

  do j = 1, 6
     resc12 = resc12+cheb12(k)*chebmo(m,k)
     ress12 = ress12+cheb12(k+1)*chebmo(m,k+1)
     estc = estc+abs((cheb12(k)-cheb24(k))*chebmo(m,k))
     ests = ests+abs((cheb12(k+1)-cheb24(k+1))*chebmo(m,k+1))
     k = k-2
  enddo

  resc24 = cheb24(25)*chebmo(m,25)
  ress24 = 0.0_dp
  resabs = abs(cheb24(25))
  k = 23

  do j = 1, 12
     resc24 = resc24+cheb24(k)*chebmo(m,k)
     ress24 = ress24+cheb24(k+1)*chebmo(m,k+1)
     resabs = resabs+abs(cheb24(k))+abs(cheb24(k+1))
     if ( j <= 5 ) then
        estc = estc+abs(cheb24(k)*chebmo(m,k))
        ests = ests+abs(cheb24(k+1)*chebmo(m,k+1))
     endif
     k = k-2
  enddo

  resabs = resabs * abs ( hlgth )

  if ( integr == 1 ) then
     result = conc * resc24-cons*ress24
     abserr = abs ( conc * estc ) + abs ( cons * ests )
  else
     result = conc*ress24+cons*resc24
     abserr = abs(conc*ests)+abs(cons*estc)
  endif

end subroutine qc25o

!===============================================================================
!###############################################################################
!===============================================================================
subroutine qc25s(f,a,b,bl,br,alfa,beta,ri,rj,rg,rh,result,abserr,resasc,integr,neval)
!
!******************************************************************************
!
!! QC25S returns rules for algebraico-logarithmic end point singularities.
!
!
! Discussion:
!
!   This routine computes
!     i = integral of F(X) * W(X) over (bl,br),
!   with error estimate, where the weight function W(X) has a singular
!   behavior of algebraico-logarithmic type at the points
!   a and/or b.
!
!   The interval (bl,br) is a subinterval of (a,b).
!
! Reference:
!
!   R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
!   QUADPACK, a Subroutine Package for Automatic Integration,
!   Springer Verlag, 1983
!
! Parameters:
!
!   Input, external real F, the name of the function routine, of the form
!     function f ( x )
!     real f
!     real x
!   which evaluates the integrand function.
!
!   Input, real A, B, the limits of integration.
!
!   Input, real BL, BR, the lower and upper limits of integration.
!   A <= BL < BR <= B.
!
!   Input, real ALFA, BETA, parameters in the weight function.
!
!   Input, real RI(25), RJ(25), RG(25), RH(25), modified Chebyshev moments
!   for the application of the generalized Clenshaw-Curtis method,
!   computed in QMOMO.
!
!   Output, real RESULT, the estimated value of the integral, computed by
!   using a generalized clenshaw-curtis method if b1 = a or br = b.
!   In all other cases the 15-point Kronrod rule is applied, obtained by
!   optimal addition of abscissae to the 7-point Gauss rule.
!
!   Output, real ABSERR, an estimate of || I - RESULT ||.
!
!   Output, real RESASC, approximation to the integral of abs(F*W-I/(B-A)).
!
!   Input, integer INTEGR,  determines the weight function
!   1, w(x) = (x-a)**alfa*(b-x)**beta
!   2, w(x) = (x-a)**alfa*(b-x)**beta*log(x-a)
!   3, w(x) = (x-a)**alfa*(b-x)**beta*log(b-x)
!   4, w(x) = (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)
!
!   Output, integer NEVAL, the number of times the integral was evaluated.
!
! Local Parameters:
!
!          fval   - value of the function f at the points
!                   (br-bl)*0.5*cos(k*pi/24)+(br+bl)*0.5
!                   k = 0, ..., 24
!          cheb12 - coefficients of the Chebyshev series expansion
!                   of degree 12, for the function f, in the interval
!                   (bl,br)
!          cheb24 - coefficients of the Chebyshev series expansion
!                   of degree 24, for the function f, in the interval
!                   (bl,br)
!          res12  - approximation to the integral obtained from cheb12
!          res24  - approximation to the integral obtained from cheb24
!          qwgts  - external function subprogram defining the four
!                   possible weight functions
!          hlgth  - half-length of the interval (bl,br)
!          centr  - mid point of the interval (bl,br)
!
!          the vector x contains the values cos(k*pi/24)
!          k = 1, ..., 11, to be used for the computation of the
!          Chebyshev series expansion of f.

!  REAL(dp), external    :: f
  interface
    function f(x)
      use kinds
      implicit none
      real(dp), intent(in) :: x
      real(dp)             :: f
    end function f
  end interface
  integer,  intent(in)  :: integr
  real(dp), intent(in)  :: a,b,bl,br,alfa,beta,ri(25),rj(25),rg(25),rh(25)
  integer,  intent(out) :: neval
  real(dp), intent(out) :: result,abserr,resasc

! Declare local variables
  integer  :: i,isym
  real(dp) :: centr,cheb12(13),cheb24(25),dc,factor,fix,fval(25),hlgth,resabs, &
              res12,res24,u

  real(dp), dimension ( 11 ) :: x = (/ &
    9.914448613738104e-01_dp,     9.659258262890683e-01_dp, &
    9.238795325112868e-01_dp,     8.660254037844386e-01_dp, &
    7.933533402912352e-01_dp,     7.071067811865475e-01_dp, &
    6.087614290087206e-01_dp,     5.000000000000000e-01_dp, &
    3.826834323650898e-01_dp,     2.588190451025208e-01_dp, &
    1.305261922200516e-01_dp /)

  neval = 25

  if ( bl == a .and. (alfa /= 0.0_dp .or. integr == 2 .or. integr == 4)) &
       go to 10

  if ( br == b .and. (beta /= 0.0_dp .or. integr == 3 .or. integr == 4)) &
       go to 140

! If a > bl and b < br, apply the 15-point Gauss-Kronrod scheme.
  call qk15w(f,qwgts,a,b,alfa,beta,integr,bl,br,result,abserr,resabs,resasc)

  neval = 15
  return

! This part of the program is executed only if a = bl.
!
! Compute the Chebyshev series expansion of the function
! f1 = (0.5*(b+b-br-a)-0.5*(br-a)*x)**beta*f(0.5*(br-a)*x+0.5*(br+a))
10 continue

  hlgth = 5.0e-01_dp*(br-bl)
  centr = 5.0e-01_dp*(br+bl)
  fix = b-centr
  fval(1) = 5.0e-01_dp*f(hlgth+centr)*(fix-hlgth)**beta
  fval(13) = f(centr)*(fix**beta)
  fval(25) = 5.0e-01_dp*f(centr-hlgth)*(fix+hlgth)**beta

  do i = 2, 12
     u = hlgth*x(i-1)
     isym = 26-i
     fval(i) = f(u+centr)*(fix-u)**beta
     fval(isym) = f(centr-u)*(fix+u)**beta
  enddo

  factor = hlgth**(alfa+1.0_dp)
  result = 0.0_dp
  abserr = 0.0_dp
  res12  = 0.0_dp
  res24  = 0.0_dp

  if ( integr > 2 ) go to 70

  call qcheb(x,fval,cheb12,cheb24)

! integr = 1  (or 2)
  do i = 1, 13
     res12 = res12+cheb12(i)*ri(i)
     res24 = res24+cheb24(i)*ri(i)
  enddo

  do i = 14, 25
     res24 = res24 + cheb24(i) * ri(i)
  enddo

  if ( integr == 1 ) go to 130

! integr = 2
  dc     = log ( br - bl )
  result = res24 * dc
  abserr = abs((res24-res12)*dc)
  res12  = 0.0_dp
  res24  = 0.0_dp

  do i = 1, 13
     res12 = res12+cheb12(i)*rg(i)
     res24 = res24+cheb24(i)*rg(i)
  enddo

  do i = 14, 25
     res24 = res24+cheb24(i)*rg(i)
  enddo

  go to 130

! Compute the Chebyshev series expansion of the function
! F4 = f1*log(0.5*(b+b-br-a)-0.5*(br-a)*x)
70 continue

  fval(1)  = fval(1) * log ( fix - hlgth )
  fval(13) = fval(13) * log ( fix )
  fval(25) = fval(25) * log ( fix + hlgth )

  do i = 2, 12
     u = hlgth*x(i-1)
     isym = 26-i
     fval(i) = fval(i) * log ( fix - u )
     fval(isym) = fval(isym) * log ( fix + u )
  enddo

  call qcheb(x,fval,cheb12,cheb24)

! integr = 3  (or 4)
  do i = 1, 13
     res12 = res12+cheb12(i)*ri(i)
     res24 = res24+cheb24(i)*ri(i)
  enddo

  do i = 14, 25
     res24 = res24+cheb24(i)*ri(i)
  enddo

  if ( integr == 3 ) go to 130

! integr = 4
  dc     = log ( br - bl )
  result = res24*dc
  abserr = abs((res24-res12)*dc)
  res12  = 0.0_dp
  res24  = 0.0_dp

  do i = 1, 13
     res12 = res12+cheb12(i)*rg(i)
     res24 = res24+cheb24(i)*rg(i)
  enddo

  do i = 14, 25
     res24 = res24+cheb24(i)*rg(i)
  enddo
130 continue

  result = (result+res24)*factor
  abserr = (abserr+abs(res24-res12))*factor
  return

! This part of the program is executed only if b = br.
!
! Compute the Chebyshev series expansion of the function
! f2 = (0.5*(b+bl-a-a)+0.5*(b-bl)*x)**alfa*f(0.5*(b-bl)*x+0.5*(b+bl))
140 continue

  hlgth = 5.0e-01_dp*(br-bl)
  centr = 5.0e-01_dp*(br+bl)
  fix = centr-a
  fval(1) = 5.0e-01_dp*f(hlgth+centr)*(fix+hlgth)**alfa
  fval(13) = f(centr)*(fix**alfa)
  fval(25) = 5.0e-01_dp*f(centr-hlgth)*(fix-hlgth)**alfa

  do i = 2, 12
     u = hlgth*x(i-1)
     isym = 26-i
     fval(i) = f(u+centr)*(fix+u)**alfa
     fval(isym) = f(centr-u)*(fix-u)**alfa
  enddo

  factor = hlgth**(beta+1.0_dp)
  result = 0.0_dp
  abserr = 0.0_dp
  res12  = 0.0_dp
  res24  = 0.0_dp

  if ( integr == 2 .or. integr == 4 ) go to 200

! integr = 1  (or 3)
  call qcheb(x,fval,cheb12,cheb24)

  do i = 1, 13
     res12 = res12+cheb12(i)*rj(i)
     res24 = res24+cheb24(i)*rj(i)
  enddo

  do i = 14, 25
     res24 = res24+cheb24(i)*rj(i)
  enddo

  if ( integr == 1 ) go to 260

! integr = 3
  dc = log ( br - bl )
  result = res24*dc
  abserr = abs((res24-res12)*dc)
  res12  = 0.0_dp
  res24  = 0.0_dp

  do i = 1, 13
     res12 = res12+cheb12(i)*rh(i)
     res24 = res24+cheb24(i)*rh(i)
  enddo

  do i = 14, 25
     res24 = res24+cheb24(i)*rh(i)
  enddo

  go to 260

! Compute the Chebyshev series expansion of the function
! f3 = f2*log(0.5*(b-bl)*x+0.5*(b+bl-a-a))
200 continue

  fval(1)  = fval(1) * log ( hlgth + fix )
  fval(13) = fval(13) * log ( fix )
  fval(25) = fval(25) * log ( fix - hlgth )

  do i = 2, 12
     u = hlgth*x(i-1)
     isym = 26-i
     fval(i) = fval(i) * log(u+fix)
     fval(isym) = fval(isym) * log(fix-u)
  enddo

  call qcheb(x,fval,cheb12,cheb24)

! integr = 2  (or 4)
  do i = 1, 13
     res12 = res12+cheb12(i)*rj(i)
     res24 = res24+cheb24(i)*rj(i)
  enddo

  do i = 14, 25
     res24 = res24+cheb24(i)*rj(i)
  enddo

  if ( integr == 2 ) go to 260

  dc     = log(br-bl)
  result = res24*dc
  abserr = abs((res24-res12)*dc)
  res12  = 0.0_dp
  res24  = 0.0_dp

! integr = 4
  do i = 1, 13
     res12 = res12+cheb12(i)*rh(i)
     res24 = res24+cheb24(i)*rh(i)
  enddo

  do i = 14, 25
     res24 = res24+cheb24(i)*rh(i)
  enddo

260 continue

  result = (result+res24)*factor
  abserr = (abserr+abs(res24-res12))*factor

end subroutine qc25s

!===============================================================================
!###############################################################################
!===============================================================================
subroutine qcheb(x,fval,cheb12,cheb24)
!
!******************************************************************************
!
!! QCHEB computes the Chebyshev series expansion.
!
!
! Discussion:
!
!   This routine computes the Chebyshev series expansion
!   of degrees 12 and 24 of a function using a fast Fourier transform method
!
!     f(x) = sum(k=1, ...,13) (cheb12(k)*t(k-1,x)),
!     f(x) = sum(k=1, ...,25) (cheb24(k)*t(k-1,x)),
!
!   where T(K,X) is the Chebyshev polynomial of degree K.
!
! Reference:
!
!   R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
!   QUADPACK, a Subroutine Package for Automatic Integration,
!   Springer Verlag, 1983
!
! Parameters:
!
!   Input, real X(11), contains the values of COS(K*PI/24), for K = 1 to 11.
!
!   Input/output, real FVAL(25), the function values at the points
!   (b+a+(b-a)*cos(k*pi/24))/2, k = 0, ...,24, where (a,b) is the
!   approximation interval.  FVAL(1) and FVAL(25) are divided by two
!   These values are destroyed at output.
!
!         on return
!          cheb12 - real
!                   vector of dimension 13 containing the Chebyshev
!                   coefficients for degree 12
!
!          cheb24 - real
!                   vector of dimension 25 containing the Chebyshev
!                   coefficients for degree 24

  real(dp), intent(in)  :: x(11)
  real(dp), intent(out) :: fval(25),cheb12(13),cheb24(25)

! Declare local variables
  integer  :: i,j
  real(dp) :: alam,alam1,alam2,part1,part2,part3,v(12)


  do i = 1, 12
     j = 26-i
     v(i) = fval(i)-fval(j)
     fval(i) = fval(i)+fval(j)
  enddo

  alam1 = v(1)-v(9)
  alam2 = x(6)*(v(3)-v(7)-v(11))
  cheb12(4) = alam1+alam2
  cheb12(10) = alam1-alam2
  alam1 = v(2)-v(8)-v(10)
  alam2 = v(4)-v(6)-v(12)
  alam = x(3)*alam1+x(9)*alam2
  cheb24(4) = cheb12(4)+alam
  cheb24(22) = cheb12(4)-alam
  alam = x(9)*alam1-x(3)*alam2
  cheb24(10) = cheb12(10)+alam
  cheb24(16) = cheb12(10)-alam
  part1 = x(4)*v(5)
  part2 = x(8)*v(9)
  part3 = x(6)*v(7)
  alam1 = v(1)+part1+part2
  alam2 = x(2)*v(3)+part3+x(10)*v(11)
  cheb12(2) = alam1+alam2
  cheb12(12) = alam1-alam2
  alam = x(1)*v(2)+x(3)*v(4)+x(5)*v(6)+x(7)*v(8) &
       +x(9)*v(10)+x(11)*v(12)
  cheb24(2) = cheb12(2)+alam
  cheb24(24) = cheb12(2)-alam
  alam = x(11)*v(2)-x(9)*v(4)+x(7)*v(6)-x(5)*v(8) &
       +x(3)*v(10)-x(1)*v(12)
  cheb24(12) = cheb12(12)+alam
  cheb24(14) = cheb12(12)-alam
  alam1 = v(1)-part1+part2
  alam2 = x(10)*v(3)-part3+x(2)*v(11)
  cheb12(6) = alam1+alam2
  cheb12(8) = alam1-alam2
  alam = x(5)*v(2)-x(9)*v(4)-x(1)*v(6) &
       -x(11)*v(8)+x(3)*v(10)+x(7)*v(12)
  cheb24(6) = cheb12(6)+alam
  cheb24(20) = cheb12(6)-alam
  alam = x(7)*v(2)-x(3)*v(4)-x(11)*v(6)+x(1)*v(8) &
       -x(9)*v(10)-x(5)*v(12)
  cheb24(8) = cheb12(8)+alam
  cheb24(18) = cheb12(8)-alam

  do i = 1, 6
     j = 14-i
     v(i) = fval(i)-fval(j)
     fval(i) = fval(i)+fval(j)
  enddo

  alam1 = v(1)+x(8)*v(5)
  alam2 = x(4)*v(3)
  cheb12(3) = alam1+alam2
  cheb12(11) = alam1-alam2
  cheb12(7) = v(1)-v(5)
  alam = x(2)*v(2)+x(6)*v(4)+x(10)*v(6)
  cheb24(3) = cheb12(3)+alam
  cheb24(23) = cheb12(3)-alam
  alam = x(6)*(v(2)-v(4)-v(6))
  cheb24(7) = cheb12(7)+alam
  cheb24(19) = cheb12(7)-alam
  alam = x(10)*v(2)-x(6)*v(4)+x(2)*v(6)
  cheb24(11) = cheb12(11)+alam
  cheb24(15) = cheb12(11)-alam

  do i = 1, 3
     j = 8-i
     v(i) = fval(i)-fval(j)
     fval(i) = fval(i)+fval(j)
  enddo

  cheb12(5) = v(1)+x(8)*v(3)
  cheb12(9) = fval(1)-x(8)*fval(3)
  alam = x(4)*v(2)
  cheb24(5) = cheb12(5)+alam
  cheb24(21) = cheb12(5)-alam
  alam = x(8)*fval(2)-fval(4)
  cheb24(9) = cheb12(9)+alam
  cheb24(17) = cheb12(9)-alam
  cheb12(1) = fval(1)+fval(3)
  alam = fval(2)+fval(4)
  cheb24(1) = cheb12(1)+alam
  cheb24(25) = cheb12(1)-alam
  cheb12(13) = v(1)-v(3)
  cheb24(13) = cheb12(13)
  alam = 1.0_dp/6.0_dp

  do i = 2, 12
     cheb12(i) = cheb12(i)*alam
  enddo

  alam = 5.0e-01_dp*alam
  cheb12(1) = cheb12(1)*alam
  cheb12(13) = cheb12(13)*alam

  do i = 2, 24
     cheb24(i) = cheb24(i)*alam
  enddo

  cheb24(1)  = 0.5_dp * alam*cheb24(1)
  cheb24(25) = 0.5_dp * alam*cheb24(25)

end subroutine qcheb

!===============================================================================
!###############################################################################
!===============================================================================
subroutine qextr(n,epstab,result,abserr,res3la,nres )
!
!******************************************************************************
!
!! QEXTR carries out the Epsilon extrapolation algorithm.
!
!
! Discussion:
!
!   The routine determines the limit of a given sequence of approximations,
!   by means of the epsilon algorithm of P. Wynn.  An estimate of the
!   absolute error is also given.  The condensed epsilon table is computed.
!   Only those elements needed for the computation of the next diagonal
!   are preserved.
!
! Reference:
!
!   R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
!   QUADPACK, a Subroutine Package for Automatic Integration,
!   Springer Verlag, 1983
!
! Parameters:
!
!   Input, integer N, indicates the entry of EPSTAB which contains
!   the new element in the first column of the epsilon table.
!
!   Input/output, real EPSTAB(52), the two lower diagonals of the triangular
!   epsilon table.  The elements are numbered starting at the right-hand
!   corner of the triangle.
!
!   Output, real RESULT, the estimated value of the integral.
!
!   Output, real ABSERR, estimate of the absolute error computed from
!   RESULT and the 3 previous results.
!
!   ?, real RES3LA(3), the last 3 results.
!
!   Input/output, integer NRES, the number of calls to the routine.  This
!   should be zero on the first call, and is automatically updated
!   before return.
!
! Local Parameters:
!
!          e0     - the 4 elements on which the
!          e1       computation of a new element in
!          e2       the epsilon table is based
!          e3                 e0
!                       e3    e1    new
!                             e2
!          newelm - number of elements to be computed in the new
!                   diagonal
!          error  - error = abs(e1-e0)+abs(e2-e1)+abs(new-e2)
!          result - the element in the new diagonal with least value
!                   of error
!          limexp is the maximum number of elements the epsilon table
!          can contain. if this number is reached, the upper diagonal
!          of the epsilon table is deleted.

  integer,  intent(inout) :: n
  real(dp), intent(inout) :: epstab(52)
  integer,  intent(out)   :: nres
  real(dp), intent(out)   :: result,abserr,res3la(3)

! Declare local variables
  integer  :: i,ib,ib2,ie,indx,k1,k2,k3,limexp,newelm,num
  real(dp) :: delta1,delta2,delta3,epsinf,error,err1,err2,err3,e0,e1,e1abs,e2, &
              e3,res,ss,tol1,tol2,tol3

  nres   = nres+1
  abserr = huge ( abserr )
  result = epstab(n)

  if ( n < 3 ) go to 100
  limexp = 50
  epstab(n+2) = epstab(n)
  newelm = (n-1)/2
  epstab(n) = huge ( epstab(n) )
  num = n
  k1 = n

  do i = 1, newelm
     k2 = k1-1
     k3 = k1-2
     res = epstab(k1+2)
     e0 = epstab(k3)
     e1 = epstab(k2)
     e2 = res
     e1abs = abs(e1)
     delta2 = e2-e1
     err2 = abs(delta2)
     tol2 = max ( abs(e2),e1abs)* epsilon ( e2 )
     delta3 = e1-e0
     err3 = abs(delta3)
     tol3 = max ( e1abs,abs(e0))* epsilon ( e0 )

! If e0, e1 and e2 are equal to within machine accuracy, convergence
! is assumed.
     if ( err2 <= tol2 .and. err3 <= tol3 ) then
        result = res
        abserr = err2+err3
        go to 100
     endif

     e3 = epstab(k1)
     epstab(k1) = e1
     delta1 = e1-e3
     err1 = abs(delta1)
     tol1 = max ( e1abs,abs(e3))* epsilon ( e3 )

! If two elements are very close to each other, omit a part
! of the table by adjusting the value of N.
     if ( err1 <= tol1 .or. err2 <= tol2 .or. err3 <= tol3 ) go to 20

     ss = 1.0_dp/delta1+1.0_dp/delta2-1.0_dp/delta3
     epsinf = abs ( ss*e1 )

! Test to detect irregular behavior in the table, and
! eventually omit a part of the table adjusting the value of N.
     if ( epsinf > 1.0e-04_dp ) go to 30

20   continue
     n = i+i-1
     exit

! Compute a new element and eventually adjust the value of RESULT.
30   continue

     res = e1+1.0_dp/ss
     epstab(k1) = res
     k1 = k1-2
     error = err2+abs(res-e2)+err3

     if ( error <= abserr ) then
        abserr = error
        result = res
     endif
  enddo

! Shift the table.
  if ( n == limexp ) then
     n = 2*(limexp/2)-1
  endif

  if ( (num/2)*2 == num ) then
     ib = 2
  else
     ib = 1
  endif

  ie = newelm+1

  do i = 1, ie
     ib2 = ib+2
     epstab(ib) = epstab(ib2)
     ib = ib2
  enddo

  if ( num /= n ) then
     indx = num-n+1
     do i = 1, n
        epstab(i)= epstab(indx)
        indx = indx+1
     enddo
  endif

  if ( nres < 4 ) then
     res3la(nres) = result
     abserr = huge ( abserr )
  else
     abserr = abs(result-res3la(3))+abs(result-res3la(2)) &
          +abs(result-res3la(1))
     res3la(1) = res3la(2)
     res3la(2) = res3la(3)
     res3la(3) = result
  endif

100 continue

  abserr = max ( abserr,0.5_dp* epsilon ( result ) *abs(result))

end subroutine qextr

!===============================================================================
!###############################################################################
!===============================================================================
subroutine qfour(f,a,b,omega,integr,epsabs,epsrel,limit,icall,maxp1,result, &
                  abserr,neval,ifail,alist,blist,rlist,elist,iord,nnlog,momcom,chebmo)
!
!******************************************************************************
!
!! QFOUR estimates the integrals of oscillatory functions.
!
!
! Discussion:
!
!   This routine calculates an approximation RESULT to a definite integral
!     I = integral of F(X) * COS(OMEGA*X)
!   or
!     I = integral of F(X) * SIN(OMEGA*X)
!   over (A,B), hopefully satisfying:
!     | I - RESULT | <= max ( epsabs, epsrel * |I| ) ).
!
!   QFOUR is called by QAWO and QAWF.  It can also be called directly in
!   a user-written program.  In the latter case it is possible for the
!   user to determine the first dimension of array CHEBMO(MAXP1,25).
!   See also parameter description of MAXP1.  Additionally see
!   parameter description of ICALL for eventually re-using
!   Chebyshev moments computed during former call on subinterval
!   of equal length abs(B-A).
!
! Reference:
!
!   R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
!   QUADPACK, a Subroutine Package for Automatic Integration,
!   Springer Verlag, 1983
!
! Parameters:
!
!   Input, external real F, the name of the function routine, of the form
!     function f ( x )
!     real f
!     real x
!   which evaluates the integrand function.
!
!   Input, real A, B, the limits of integration.
!
!   Input, real OMEGA, the multiplier of X in the weight function.
!
!   Input, integer INTEGR, indicates the weight functions to be used.
!   = 1, w(x) = cos(omega*x)
!   = 2, w(x) = sin(omega*x)
!
!   Input, real EPSABS, EPSREL, the absolute and relative accuracy requested.
!
!   Input, integer LIMIT, the maximum number of subintervals of [A,B]
!   that can be generated.
!
!           icall  - integer
!                    if qfour is to be used only once, ICALL must
!                    be set to 1.  assume that during this call, the
!                    Chebyshev moments (for clenshaw-curtis integration
!                    of degree 24) have been computed for intervals of
!                    lenghts (abs(b-a))*2**(-l), l=0,1,2,...momcom-1.
!                    the Chebyshev moments already computed can be
!                    re-used in subsequent calls, if qfour must be
!                    called twice or more times on intervals of the
!                    same length abs(b-a). from the second call on, one
!                    has to put then ICALL > 1.
!                    if ICALL < 1, the routine will end with ifail = 6.
!
!           maxp1  - integer
!                    gives an upper bound on the number of
!                    Chebyshev moments which can be stored, i.e.
!                    for the intervals of lenghts abs(b-a)*2**(-l),
!                    l=0,1, ..., maxp1-2, maxp1 >= 1.
!                    if maxp1 < 1, the routine will end with ifail = 6.
!                    increasing (decreasing) the value of maxp1
!                    decreases (increases) the computational time but
!                    increases (decreases) the required memory space.
!
!   Output, real RESULT, the estimated value of the integral.
!
!   Output, real ABSERR, an estimate of || I - RESULT ||.
!
!   Output, integer NEVAL, the number of times the integral was evaluated.
!
!         ifail    - integer
!                  ifail = 0 normal and reliable termination of the
!                            routine. it is assumed that the
!                            requested accuracy has been achieved.
!                - ifail > 0 abnormal termination of the routine.
!                            the estimates for integral and error are
!                            less reliable. it is assumed that the
!                            requested accuracy has not been achieved.
!                  ifail = 1 maximum number of subdivisions allowed
!                            has been achieved. one can allow more
!                            subdivisions by increasing the value of
!                            limit (and taking according dimension
!                            adjustments into account). however, if
!                            this yields no improvement it is advised
!                            to analyze the integrand, in order to
!                            determine the integration difficulties.
!                            if the position of a local difficulty can
!                            be determined (e.g. singularity,
!                            discontinuity within the interval) one
!                            will probably gain from splitting up the
!                            interval at this point and calling the
!                            integrator on the subranges. if possible,
!                            an appropriate special-purpose integrator
!                            should be used which is designed for
!                            handling the type of difficulty involved.
!                        = 2 the occurrence of roundoff error is
!                            detected, which prevents the requested
!                            tolerance from being achieved.
!                            the error may be under-estimated.
!                        = 3 extremely bad integrand behavior occurs
!                            at some points of the integration
!                            interval.
!                        = 4 the algorithm does not converge. roundoff
!                            error is detected in the extrapolation
!                            table. it is presumed that the requested
!                            tolerance cannot be achieved due to
!                            roundoff in the extrapolation table, and
!                            that the returned result is the best which
!                            can be obtained.
!                        = 5 the integral is probably divergent, or
!                            slowly convergent. it must be noted that
!                            divergence can occur with any other value
!                            of ifail > 0.
!                        = 6 the input is invalid, because
!                            epsabs < 0 and epsrel < 0,
!                            or (integr /= 1 and integr /= 2) or
!                            ICALL < 1 or maxp1 < 1.
!                            result, abserr, neval, last, rlist(1),
!                            elist(1), iord(1) and nnlog(1) are set to
!                            zero. alist(1) and blist(1) are set to a
!                            and b respectively.
!
!   Workspace, real ALIST(LIMIT), BLIST(LIMIT), contains in entries 1
!   through LAST the left and right ends of the partition subintervals.
!
!   Workspace, real RLIST(LIMIT), contains in entries 1 through LAST
!   the integral approximations on the subintervals.
!
!   Workspace, real ELIST(LIMIT), contains in entries 1 through LAST
!   the absolute error estimates on the subintervals.
!
!           iord   - integer
!                    vector of dimension at least limit, the first k
!                    elements of which are pointers to the error
!                    estimates over the subintervals, such that
!                    elist(iord(1)), ..., elist(iord(k)), form
!                    a decreasing sequence, with k = last
!                    if last <= (limit/2+2), and
!                    k = limit+1-last otherwise.
!
!           nnlog  - integer
!                    vector of dimension at least limit, indicating the
!                    subdivision levels of the subintervals, i.e.
!                    iwork(i) = l means that the subinterval numbered
!                    i is of length abs(b-a)*2**(1-l)
!
!        on entry and return
!           momcom - integer
!                    indicating that the Chebyshev moments have been
!                    computed for intervals of lengths
!                    (abs(b-a))*2**(-l), l=0,1,2, ..., momcom-1,
!                    momcom < maxp1
!
!           chebmo - real
!                    array of dimension (maxp1,25) containing the
!                    Chebyshev moments
!
! Local Parameters:
!
!          alist     - list of left end points of all subintervals
!                      considered up to now
!          blist     - list of right end points of all subintervals
!                      considered up to now
!          rlist(i)  - approximation to the integral over
!                      (alist(i),blist(i))
!          rlist2    - array of dimension at least limexp+2 containing
!                      the part of the epsilon table which is still
!                      needed for further computations
!          elist(i)  - error estimate applying to rlist(i)
!          maxerr    - pointer to the interval with largest error
!                      estimate
!          errmax    - elist(maxerr)
!          erlast    - error on the interval currently subdivided
!          area      - sum of the integrals over the subintervals
!          errsum    - sum of the errors over the subintervals
!          errbnd    - requested accuracy max(epsabs,epsrel*
!                      abs(result))
!          *****1    - variable for the left subinterval
!          *****2    - variable for the right subinterval
!          last      - index for subdivision
!          nres      - number of calls to the extrapolation routine
!          numrl2    - number of elements in rlist2. if an appropriate
!                      approximation to the compounded integral has
!                      been obtained it is put in rlist2(numrl2) after
!                      numrl2 has been increased by one
!          small     - length of the smallest interval considered
!                      up to now, multiplied by 1.5
!          erlarg    - sum of the errors over the intervals larger
!                      than the smallest interval considered up to now
!          extrap    - logical variable denoting that the routine is
!                      attempting to perform extrapolation, i.e. before
!                      subdividing the smallest interval we try to
!                      decrease the value of erlarg
!          noext     - logical variable denoting that extrapolation
!                      is no longer allowed (true value)

!  REAL(dp), external    :: f
  interface
    function f(x)
      use kinds
      implicit none
      real(dp), intent(in) :: x
      real(dp)             :: f
    end function f
  end interface
  integer,  intent(in)  :: integr,limit,maxp1,icall
  real(dp), intent(in)  :: a,b,omega,epsabs,epsrel
  integer,  intent(out) :: neval,ifail,momcom,iord(limit),nnlog(limit)
  real(dp), intent(out) :: result,abserr,chebmo(maxp1,25), &
                           alist(limit),blist(limit),elist(limit),rlist(limit)

! Declare local variables
  logical  :: extall,extrap,noext
  integer  :: id,ierro,iroff1,iroff2,iroff3,jupbnd,k,ksgn,ktmin,last,maxerr,nev, &
              nres,nrmax,nrmom,numrl2
  real(dp) :: abseps,area,area1,area12,area2,a1,a2,b1,b2,correc,defab1,defab2,defabs, &
              domega,dres,erlarg,erlast,errbnd,errmax,error1,erro12,error2,errsum,ertest, &
              resabs,reseps,res3la(3),rlist2(52),small,width

! The dimension of rlist2 is determined by  the value of
! limexp in QEXTR (rlist2 should be of dimension
! (limexp+2) at least).

! Test on validity of parameters.
  ifail    = 0
  neval    = 0
  last     = 0
  result   = 0.0_dp
  abserr   = 0.0_dp
  alist(1) = a
  blist(1) = b
  rlist(1) = 0.0_dp
  elist(1) = 0.0_dp
  iord(1)  = 0
  nnlog(1) = 0

  if ( (integr /= 1.and.integr /= 2) .or. (epsabs < 0.0_dp.and. &
       epsrel < 0.0_dp) .or. icall < 1 .or. maxp1 < 1 ) then
     ifail = 6
     return
  endif

! First approximation to the integral.
  domega = abs ( omega )
  nrmom  = 0

  if ( icall <= 1 ) momcom = 0

  call qc25o(f,a,b,domega,integr,nrmom,maxp1,0,result,abserr, &
              neval,defabs,resabs,momcom,chebmo)

! Test on accuracy.
  dres     = abs(result)
  errbnd   = max ( epsabs,epsrel*dres)
  rlist(1) = result
  elist(1) = abserr
  iord(1)  = 1
  if ( abserr <= 1.0e+02_dp* epsilon ( defabs ) *defabs .and. &
       abserr > errbnd ) ifail = 2

  if ( limit == 1 ) ifail = 1

  if ( ifail /= 0 .or. abserr <= errbnd ) go to 200

! Initializations
  errmax = abserr
  maxerr = 1
  area   = result
  errsum = abserr
  abserr = huge ( abserr )
  nrmax  = 1
  extrap = .false.
  noext  = .false.
  ierro  = 0
  iroff1 = 0
  iroff2 = 0
  iroff3 = 0
  ktmin  = 0
  small  = abs(b-a)*7.5e-01_dp
  nres   = 0
  numrl2 = 0
  extall = .false.

  if ( 5.0e-01_dp*abs(b-a)*domega <= 2.0_dp) then
     numrl2 = 1
     extall = .true.
     rlist2(1) = result
  endif

  if ( 2.5e-01_dp * abs(b-a) * domega <= 2.0_dp ) extall = .true.

  if ( dres >= (1.0_dp-5.0e+01_dp* epsilon ( defabs ) )*defabs ) then
     ksgn = 1
  else
     ksgn = -1
  endif

! main do-loop
  do 140 last = 2, limit

! Bisect the subinterval with the nrmax-th largest error estimate.
     nrmom = nnlog(maxerr)+1
     a1 = alist(maxerr)
     b1 = 5.0e-01_dp*(alist(maxerr)+blist(maxerr))
     a2 = b1
     b2 = blist(maxerr)
     erlast = errmax

     call qc25o(f,a1,b1,domega,integr,nrmom,maxp1,0,area1,&
                 error1,nev,resabs,defab1,momcom,chebmo)
     neval = neval+nev
     call qc25o(f,a2,b2,domega,integr,nrmom,maxp1,1,area2,&
                 error2,nev,resabs,defab2,momcom,chebmo)
     neval = neval+nev

! Improve previous approximations to integral and error and
! test for accuracy.
     area12 = area1+area2
     erro12 = error1+error2
     errsum = errsum+erro12-errmax
     area   = area+area12-rlist(maxerr)
     if ( defab1 == error1 .or. defab2 == error2 ) go to 25
     if ( abs(rlist(maxerr)-area12) > 1.0e-05*abs(area12) &
          .or. erro12 < 9.9e-01_dp*errmax ) go to 20
     if ( extrap ) iroff2 = iroff2+1

     if ( .not.extrap ) iroff1 = iroff1+1

20   continue

     if ( last > 10.and.erro12 > errmax ) iroff3 = iroff3+1

25  continue

     rlist(maxerr) = area1
     rlist(last)   = area2
     nnlog(maxerr) = nrmom
     nnlog(last)   = nrmom
     errbnd = max ( epsabs,epsrel*abs(area))

! Test for roundoff error and eventually set error flag
     if ( iroff1+iroff2 >= 10 .or. iroff3 >= 20 ) ifail = 2

     if ( iroff2 >= 5) ierro = 3

! Set error flag in the case that the number of subintervals
! equals limit.
     if ( last == limit ) ifail = 1

! Set error flag in the case of bad integrand behavior at
! a point of the integration range.
    if ( max ( abs(a1),abs(b2)) <= (1.0_dp+1.0e+03_dp* epsilon ( a1 ) ) &
         *(abs(a2)+1.0e+03_dp* tiny ( a2 ) )) ifail = 4

! Append the newly-created intervals to the list.
    if ( error2 <= error1 ) then
       alist(last)   = a2
       blist(maxerr) = b1
       blist(last)   = b2
       elist(maxerr) = error1
       elist(last)   = error2
    else
       alist(maxerr) = a2
       alist(last)   = a1
       blist(last)   = b1
       rlist(maxerr) = area2
       rlist(last)   = area1
       elist(maxerr) = error2
       elist(last)   = error1
    endif

! Call QSORT to maintain the descending ordering
! in the list of error estimates and select the subinterval
! with nrmax-th largest error estimate (to be bisected next).
40  continue

    call qsort(limit,last,maxerr,errmax,elist,iord,nrmax)

    if ( errsum <= errbnd ) go to 170

    if ( ifail /= 0 ) go to 150
    if ( last == 2 .and. extall ) go to 120
    if ( noext ) go to 140
    if ( .not. extall ) go to 50
    erlarg = erlarg-erlast
    if ( abs(b1-a1) > small ) erlarg = erlarg+erro12
    if ( extrap ) go to 70

! Test whether the interval to be bisected next is the
! smallest interval.
50  continue

    width = abs(blist(maxerr)-alist(maxerr))
    if ( width > small ) go to 140
    if ( extall ) go to 60

! Test whether we can start with the extrapolation procedure
! (we do this if we integrate over the next interval with
! use of a Gauss-Kronrod rule - see QC25O).
    small = small*5.0e-01_dp
    if ( 2.5e-01_dp*width*domega > 2.0_dp ) go to 140
    extall = .true.
    go to 130

60  continue

    extrap = .true.
    nrmax = 2

70  continue

    if ( ierro == 3 .or. erlarg <= ertest ) go to 90

! The smallest interval has the largest error.
! Before bisecting decrease the sum of the errors over the
! larger intervals (ERLARG) and perform extrapolation.
    jupbnd = last

    if ( last > (limit/2+2) ) jupbnd = limit+3-last

    id = nrmax

    do k = id, jupbnd
       maxerr = iord(nrmax)
       errmax = elist(maxerr)
       if ( abs(blist(maxerr)-alist(maxerr)) > small ) go to 140
       nrmax = nrmax+1
    enddo

! Perform extrapolation.
90  continue

    numrl2 = numrl2+1
    rlist2(numrl2) = area
    if ( numrl2 < 3 ) go to 110
    call qextr ( numrl2, rlist2, reseps, abseps, res3la, nres )
    ktmin = ktmin+1

    if ( ktmin > 5.and.abserr < 1.0e-03_dp*errsum ) ifail = 5

    if ( abseps >= abserr ) go to 100
    ktmin = 0
    abserr = abseps
    result = reseps
    correc = erlarg
    ertest = max ( epsabs, epsrel*abs(reseps))
    if ( abserr <= ertest ) go to 150

! Prepare bisection of the smallest interval.
100 continue

    if ( numrl2 == 1 ) noext = .true.
    if ( ifail == 5 ) go to 150

110 continue

    maxerr = iord(1)
    errmax = elist(maxerr)
    nrmax  = 1
    extrap = .false.
    small  = small*5.0e-01_dp
    erlarg = errsum
    go to 140

120 continue

    small  = small * 5.0e-01_dp
    numrl2 = numrl2 + 1
    rlist2(numrl2) = area

130 continue

    ertest = errbnd
    erlarg = errsum

140 continue

! set the final result.
150 continue

  if ( abserr == huge ( abserr ) .or. nres == 0 ) go to 170
  if ( ifail+ierro == 0 ) go to 165
  if ( ierro == 3 ) abserr = abserr+correc
  if ( ifail == 0 ) ifail = 3
  if ( result /= 0.0_dp.and.area /= 0.0_dp ) go to 160
  if ( abserr > errsum ) go to 170
  if ( area == 0.0_dp ) go to 190
  go to 165

160 continue

  if ( abserr/abs(result) > errsum/abs(area) ) go to 170

! Test on divergence.
  165 continue

  if ( ksgn == (-1) .and. max ( abs(result),abs(area)) <=  &
       defabs*1.0e-02_dp ) go to 190

  if ( 1.0e-02_dp > (result/area) .or. (result/area) > 1.0e+02_dp &
       .or. errsum >= abs(area) ) ifail = 6

  go to 190

! Compute global integral sum.
170 continue

  result = sum ( rlist(1:last) )
  abserr = errsum

190 continue

  if (ifail > 2) ifail=ifail-1

200 continue

  if ( integr == 2 .and. omega < 0.0_dp ) result = -result

end subroutine qfour

!===============================================================================
!###############################################################################
!===============================================================================
subroutine qk15(f,a,b,result,abserr,resabs,resasc)
!
!******************************************************************************
!
!! QK15 carries out a 15 point Gauss-Kronrod quadrature rule.
!
!
! Discussion:
!
!   This routine approximates
!     I = integral ( A <= X <= B ) F(X) dx
!   with an error estimate, and
!     J = integral ( A <= X <= B ) | F(X) | dx
!
! Reference:
!
!   R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
!   QUADPACK, a Subroutine Package for Automatic Integration,
!   Springer Verlag, 1983
!
! Parameters:
!
!   Input, external real F, the name of the function routine, of the form
!     function f ( x )
!     real f
!     real x
!   which evaluates the integrand function.
!
!   Input, real A, B, the limits of integration.
!
!   Output, real RESULT, the estimated value of the integral.
!   RESULT is computed by applying the 15-point Kronrod rule (RESK)
!   obtained by optimal addition of abscissae to the 7-point Gauss rule
!   (RESG).
!
!   Output, real ABSERR, an estimate of | I - RESULT |.
!
!   Output, real RESABS, approximation to the integral of the absolute
!   value of F.
!
!   Output, real RESASC, approximation to the integral | F-I/(B-A) |
!   over [A,B].
!
! Local Parameters:
!
!          the abscissae and weights are given for the interval (-1,1).
!          because of symmetry only the positive abscissae and their
!          corresponding weights are given.
!
!          xgk    - abscissae of the 15-point Kronrod rule
!                   xgk(2), xgk(4), ...  abscissae of the 7-point
!                   Gauss rule
!                   xgk(1), xgk(3), ...  abscissae which are optimally
!                   added to the 7-point Gauss rule
!
!          wgk    - weights of the 15-point Kronrod rule
!
!          wg     - weights of the 7-point Gauss rule
!
!          centr  - mid point of the interval
!          hlgth  - half-length of the interval
!          absc   - abscissa
!          fval*  - function value
!          resg   - result of the 7-point Gauss formula
!          resk   - result of the 15-point Kronrod formula
!          reskh  - approximation to the mean value of f over (a,b),
!                   i.e. to i/(b-a)

!  REAL(dp), external    :: f
  interface
    function f(x)
      use kinds
      implicit none
      real(dp), intent(in) :: x
      real(dp)             :: f
    end function f
  end interface
  real(dp), intent(in)  :: a,b
  real(dp), intent(out) :: result,abserr,resabs,resasc

! Declare local variables
  integer  :: j,jtw,jtwm1
  real(dp) :: absc,centr,dhlgth,fc,fsum,fval1,fval2,fv1(7),fv2(7),hlgth, &
              resg,resk,reskh,wg(4),wgk(8),xgk(8)

  xgk = (/ 9.914553711208126e-01_dp,   9.491079123427585e-01_dp, &
           8.648644233597691e-01_dp,   7.415311855993944e-01_dp, &
           5.860872354676911e-01_dp,   4.058451513773972e-01_dp, &
           2.077849550078985e-01_dp,   0.0_dp                  /)
  wgk = (/ 2.293532201052922e-02_dp,   6.309209262997855e-02_dp, &
           1.047900103222502e-01_dp,   1.406532597155259e-01_dp, &
           1.690047266392679e-01_dp,   1.903505780647854e-01_dp, &
           2.044329400752989e-01_dp,   2.094821410847278e-01_dp    /)
  wg  = (/ 1.294849661688697e-01_dp,   2.797053914892767e-01_dp, &
           3.818300505051189e-01_dp,   4.179591836734694e-01_dp    /)

  centr  = 5.0e-01_dp*(a+b)
  hlgth  = 5.0e-01_dp*(b-a)
  dhlgth = abs(hlgth)

! Compute the 15-point Kronrod approximation to the integral,
! and estimate the absolute error.
  fc = f(centr)
  resg = fc*wg(4)
  resk = fc*wgk(8)
  resabs = abs(resk)

  do j = 1, 3
     jtw = j*2
     absc = hlgth*xgk(jtw)
     fval1 = f(centr-absc)
     fval2 = f(centr+absc)
     fv1(jtw) = fval1
     fv2(jtw) = fval2
     fsum = fval1+fval2
     resg = resg+wg(j)*fsum
     resk = resk+wgk(jtw)*fsum
     resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
  enddo

  do j = 1, 4
     jtwm1 = j*2-1
     absc = hlgth*xgk(jtwm1)
     fval1 = f(centr-absc)
     fval2 = f(centr+absc)
     fv1(jtwm1) = fval1
     fv2(jtwm1) = fval2
     fsum = fval1+fval2
     resk = resk+wgk(jtwm1)*fsum
     resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
  enddo

  reskh = resk * 5.0e-01_dp
  resasc = wgk(8)*abs(fc-reskh)

  do j = 1, 7
     resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
  enddo

  result = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr = abs((resk-resg)*hlgth)

  if ( resasc /= 0.0_dp.and.abserr /= 0.0_dp ) &
     abserr = resasc*min ( 1.0_dp,(2.0e+02_dp*abserr/resasc)**1.5_dp)

  if ( resabs > tiny ( resabs ) / (5.0e+01_dp* epsilon ( resabs ) ) ) &
     abserr = max (( epsilon ( resabs ) *5.0e+01_dp)*resabs,abserr)

end subroutine qk15

!===============================================================================
!###############################################################################
!===============================================================================
subroutine qk15i(f,boun,inf,a,b,result,abserr,resabs,resasc)
!
!******************************************************************************
!
!! QK15I applies a 15 point Gauss-Kronrod quadrature on an infinite interval.
!
!
! Discussion:
!
!   The original infinite integration range is mapped onto the interval
!   (0,1) and (a,b) is a part of (0,1).  The routine then computes:
!
!   i = integral of transformed integrand over (a,b),
!   j = integral of abs(transformed integrand) over (a,b).
!
! Reference:
!
!   R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
!   QUADPACK, a Subroutine Package for Automatic Integration,
!   Springer Verlag, 1983
!
! Parameters:
!
!   Input, external real F, the name of the function routine, of the form
!     function f ( x )
!     real f
!     real x
!   which evaluates the integrand function.
!
!   Input, real BOUN, the finite bound of the original integration range,
!   or zero if INF is 2.
!
!             inf    - integer
!                      if inf = -1, the original interval is
!                                  (-infinity,BOUN),
!                      if inf = +1, the original interval is
!                                  (BOUN,+infinity),
!                      if inf = +2, the original interval is
!                                  (-infinity,+infinity) and
!                      The integral is computed as the sum of two
!                      integrals, one over (-infinity,0) and one
!                      over (0,+infinity).
!
!   Input, real A, B, the limits of integration, over a subrange of [0,1].
!
!   Output, real RESULT, the estimated value of the integral.
!   RESULT is computed by applying the 15-point Kronrod rule (RESK) obtained
!   by optimal addition of abscissae to the 7-point Gauss rule (RESG).
!
!   Output, real ABSERR, an estimate of | I - RESULT |.
!
!   Output, real RESABS, approximation to the integral of the absolute
!   value of F.
!
!   Output, real RESASC, approximation to the integral of the
!   transformated integrand | F-I/(B-A) | over [A,B].
!
! Local Parameters:
!
!          centr  - mid point of the interval
!          hlgth  - half-length of the interval
!          absc*  - abscissa
!          tabsc* - transformed abscissa
!          fval*  - function value
!          resg   - result of the 7-point Gauss formula
!          resk   - result of the 15-point Kronrod formula
!          reskh  - approximation to the mean value of the transformed
!                   integrand over (a,b), i.e. to i/(b-a)

!  REAL(dp), external    :: f
  interface
    function f(x)
      use kinds
      implicit none
      real(dp), intent(in) :: x
      real(dp)             :: f
    end function f
  end interface
  integer,  intent(in)  :: inf
  real(dp), intent(in)  :: a,b,boun
  real(dp), intent(out) :: result,abserr,resabs,resasc

! Declare local variables
  integer  :: j
  real(dp) :: absc,absc1,absc2,centr,dinf,fc,fsum,fval1,fval2,fv1(7),fv2(7),hlgth, &
              resg,resk,reskh,tabsc1,tabsc2,wg(8),wgk(8),xgk(8)

! the abscissae and weights are supplied for the interval
! (-1,1).  because of symmetry only the positive abscissae and
! their corresponding weights are given.
!
!          xgk    - abscissae of the 15-point Kronrod rule
!                   xgk(2), xgk(4), ... abscissae of the 7-point Gauss
!                   rule
!                   xgk(1), xgk(3), ...  abscissae which are optimally
!                   added to the 7-point Gauss rule
!
!          wgk    - weights of the 15-point Kronrod rule
!
!          wg     - weights of the 7-point Gauss rule, corresponding
!                   to the abscissae xgk(2), xgk(4), ...
!                   wg(1), wg(3), ... are set to zero.

  xgk = (/ 9.914553711208126e-01_dp,     9.491079123427585e-01_dp, &
           8.648644233597691e-01_dp,     7.415311855993944e-01_dp, &
           5.860872354676911e-01_dp,     4.058451513773972e-01_dp, &
           2.077849550078985e-01_dp,     0.0000000000000000_dp   /)

  wgk = (/ 2.293532201052922e-02_dp,     6.309209262997855e-02_dp, &
           1.047900103222502e-01_dp,     1.406532597155259e-01_dp, &
           1.690047266392679e-01_dp,     1.903505780647854e-01_dp, &
           2.044329400752989e-01_dp,     2.094821410847278e-01_dp    /)

  wg  = (/ 0.0000000000000000_dp,     1.294849661688697e-01_dp, &
           0.0000000000000000_dp,     2.797053914892767e-01_dp, &
           0.0000000000000000_dp,     3.818300505051189e-01_dp, &
           0.0000000000000000_dp,     4.179591836734694e-01_dp    /)

  dinf = min ( 1, inf )

  centr  = 5.0e-01_dp*(a+b)
  hlgth  = 5.0e-01_dp*(b-a)
  tabsc1 = boun+dinf*(1.0_dp-centr)/centr
  fval1  = f(tabsc1)
  if ( inf == 2 ) fval1 = fval1+f(-tabsc1)
  fc = (fval1/centr)/centr

! Compute the 15-point Kronrod approximation to the integral,
! and estimate the error.
  resg = wg(8)*fc
  resk = wgk(8)*fc
  resabs = abs(resk)

  do j = 1, 7
     absc   = hlgth*xgk(j)
     absc1  = centr-absc
     absc2  = centr+absc
     tabsc1 = boun+dinf*(1.0_dp-absc1)/absc1
     tabsc2 = boun+dinf*(1.0_dp-absc2)/absc2
     fval1  = f(tabsc1)
     fval2  = f(tabsc2)

    if ( inf == 2 ) then
       fval1 = fval1+f(-tabsc1)
       fval2 = fval2+f(-tabsc2)
    endif

    fval1  = (fval1/absc1)/absc1
    fval2  = (fval2/absc2)/absc2
    fv1(j) = fval1
    fv2(j) = fval2
    fsum   = fval1+fval2
    resg   = resg+wg(j)*fsum
    resk   = resk+wgk(j)*fsum
    resabs = resabs+wgk(j)*(abs(fval1)+abs(fval2))
  enddo

  reskh  = resk * 5.0e-01_dp
  resasc = wgk(8) * abs(fc-reskh)

  do j = 1, 7
     resasc = resasc + wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
  enddo

  result = resk * hlgth
  resasc = resasc * hlgth
  resabs = resabs * hlgth
  abserr = abs ( ( resk - resg ) * hlgth )

  if ( resasc /= 0.0_dp.and.abserr /= 0.0_dp) &
     abserr = resasc* min ( 1.0_dp,(2.0e+02_dp*abserr/resasc)**1.5_dp)

  if ( resabs > tiny ( resabs ) / ( 5.0e+01_dp * epsilon ( resabs ) ) ) &
    abserr = max (( epsilon ( resabs ) *5.0e+01_dp)*resabs,abserr)

end subroutine qk15i

!===============================================================================
!###############################################################################
!===============================================================================
subroutine qk15w(f,w,p1,p2,p3,p4,kp,a,b,result,abserr,resabs,resasc)
!
!******************************************************************************
!
!! QK15W applies a 15 point Gauss-Kronrod rule for a weighted integrand.
!
!
! Discussion:
!
!   This routine approximates
!     i = integral of f*w over (a,b),
!   with error estimate, and
!     j = integral of abs(f*w) over (a,b)
!
! Reference:
!
!   R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
!   QUADPACK, a Subroutine Package for Automatic Integration,
!   Springer Verlag, 1983
!
! Parameters:
!
!   Input, external real F, the name of the function routine, of the form
!     function f ( x )
!     real f
!     real x
!   which evaluates the integrand function.
!
!             w      - real
!                      function subprogram defining the integrand
!                      weight function w(x). the actual name for w
!                      needs to be declared e x t e r n a l in the
!                      calling program.
!
!   ?, real P1, P2, P3, P4, parameters in the weight function
!
!             kp     - integer
!                      key for indicating the type of weight function
!
!   Input, real A, B, the limits of integration.
!
!   Output, real RESULT, the estimated value of the integral.
!   RESULT is computed by applying the 15-point Kronrod rule (RESK) obtained by
!   optimal addition of abscissae to the 7-point Gauss rule (RESG).
!
!   Output, real ABSERR, an estimate of | I - RESULT |.
!
!   Output, real RESABS, approximation to the integral of the absolute
!   value of F.
!
!   Output, real RESASC, approximation to the integral | F-I/(B-A) |
!   over [A,B].
!
! Local Parameters:
!
!          centr  - mid point of the interval
!          hlgth  - half-length of the interval
!          absc*  - abscissa
!          fval*  - function value
!          resg   - result of the 7-point Gauss formula
!          resk   - result of the 15-point Kronrod formula
!          reskh  - approximation to the mean value of f*w over (a,b),
!                   i.e. to i/(b-a)

!  REAL(dp), external    :: f,w
  interface
    function f(x)
      use kinds
      implicit none
      real(dp)             :: f
      real(dp), intent(in) :: x
    end function f
    function w(x,c,p2,p3,p4,kp)
      use kinds
      implicit none
      real(dp)             :: w
      real(dp), intent(in) :: x,c,p2,p3,p4
      integer,  intent(in) :: kp
    end function w
  end interface
  integer,  intent(in)  :: kp
  real(dp), intent(in)  :: a,b, p1, p2, p3, p4
  real(dp), intent(out) :: result,abserr,resabs,resasc

! Declare local variables
  integer  :: j,jtw,jtwm1
  real(dp) :: absc,absc1,absc2,centr,dhlgth,fc,fsum,fval1,fval2,fv1(7),fv2(7), &
              hlgth,resg,resk,reskh,wgk(8),xgk(8),wg(4)

! the abscissae and weights are given for the interval (-1,1).
! because of symmetry only the positive abscissae and their
! corresponding weights are given.
!
!          xgk    - abscissae of the 15-point Gauss-Kronrod rule
!                   xgk(2), xgk(4), ... abscissae of the 7-point Gauss
!                   rule
!                   xgk(1), xgk(3), ... abscissae which are optimally
!                   added to the 7-point Gauss rule
!
!          wgk    - weights of the 15-point Gauss-Kronrod rule
!
!          wg     - weights of the 7-point Gauss rule

  wg  = (/ 1.294849661688697e-01_dp,     2.797053914892767e-01_dp, &
           3.818300505051889e-01_dp,     4.179591836734694e-01_dp   /)

  xgk = (/ 9.914553711208126e-01_dp,     9.491079123427585e-01_dp, &
           8.648644233597691e-01_dp,     7.415311855993944e-01_dp, &
           5.860872354676911e-01_dp,     4.058451513773972e-01_dp, &
           2.077849550789850e-01_dp,     0.000000000000000_dp    /)

  wgk = (/ 2.293532201052922e-02_dp,     6.309209262997855e-02_dp, &
           1.047900103222502e-01_dp,     1.406532597155259e-01_dp, &
           1.690047266392679e-01_dp,     1.903505780647854e-01_dp, &
           2.044329400752989e-01_dp,     2.094821410847278e-01_dp    /)

  centr = 5.0e-01_dp*(a+b)
  hlgth = 5.0e-01_dp*(b-a)
  dhlgth = abs(hlgth)

! Compute the 15-point Kronrod approximation to the integral,
! and estimate the error.
  fc = f(centr)*w(centr,p1,p2,p3,p4,kp)
  resg = wg(4)*fc
  resk = wgk(8)*fc
  resabs = abs(resk)

  do j = 1, 3
     jtw = j*2
     absc = hlgth*xgk(jtw)
     absc1 = centr-absc
     absc2 = centr+absc
     fval1 = f(absc1)*w(absc1,p1,p2,p3,p4,kp)
     fval2 = f(absc2)*w(absc2,p1,p2,p3,p4,kp)
     fv1(jtw) = fval1
     fv2(jtw) = fval2
     fsum = fval1+fval2
     resg = resg+wg(j)*fsum
     resk = resk+wgk(jtw)*fsum
     resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
  enddo

  do j = 1, 4
     jtwm1 = j*2-1
     absc = hlgth*xgk(jtwm1)
     absc1 = centr-absc
     absc2 = centr+absc
     fval1 = f(absc1)*w(absc1,p1,p2,p3,p4,kp)
     fval2 = f(absc2)*w(absc2,p1,p2,p3,p4,kp)
     fv1(jtwm1) = fval1
     fv2(jtwm1) = fval2
     fsum = fval1+fval2
     resk = resk+wgk(jtwm1)*fsum
     resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
  enddo

  reskh = resk*5.0e-01_dp
  resasc = wgk(8)*abs(fc-reskh)

  do j = 1, 7
     resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
  enddo

  result = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr = abs((resk-resg)*hlgth)

  if ( resasc /= 0.0_dp.and.abserr /= 0.0_dp) &
     abserr = resasc*min ( 1.0_dp,(2.0e+02_dp*abserr/resasc)**1.5_dp)

  if ( resabs > tiny ( resabs ) /(5.0e+01_dp* epsilon ( resabs ) ) ) &
     abserr = max ( ( epsilon ( resabs ) * 5.0e+01_dp)*resabs,abserr)

end subroutine qk15w

!===============================================================================
!###############################################################################
!===============================================================================
subroutine qk21(f,a,b,result,abserr,resabs,resasc)

!******************************************************************************
!
!! QK21 carries out a 21 point Gauss-Kronrod quadrature rule.
!
!
! Discussion:
!
!   This routine approximates
!     I = integral ( A <= X <= B ) F(X) dx
!   with an error estimate, and
!     J = integral ( A <= X <= B ) | F(X) | dx
!
! Reference:
!
!   R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
!   QUADPACK, a Subroutine Package for Automatic Integration,
!   Springer Verlag, 1983
!
! Parameters:
!
!   Input, external real F, the name of the function routine, of the form
!     function f ( x )
!     real f
!     real x
!   which evaluates the integrand function.
!
!   Input, real A, B, the limits of integration.
!
!   Output, real RESULT, the estimated value of the integral.
!                      result is computed by applying the 21-point
!                      Kronrod rule (resk) obtained by optimal addition
!                      of abscissae to the 10-point Gauss rule (resg).
!
!   Output, real ABSERR, an estimate of | I - RESULT |.
!
!   Output, real RESABS, approximation to the integral of the absolute
!   value of F.
!
!   Output, real RESASC, approximation to the integral | F-I/(B-A) |
!   over [A,B].

!  REAL(dp), external    :: f
  interface
    function f(x)
      use kinds
      implicit none
      real(dp), intent(in) :: x
      real(dp)             :: f
    end function f
  end interface
  real(dp), intent(in)  :: a,b
  real(dp), intent(out) :: result,abserr,resabs,resasc

! Declare local variables
  integer  :: j,jtw,jtwm1
  real(dp) :: absc,centr,dhlgth,fc,fsum,fval1,fval2,fv1(10),fv2(10),hlgth, &
              resg,resk,reskh,wg(5),wgk(11),xgk(11)

!          the abscissae and weights are given for the interval (-1,1).
!          because of symmetry only the positive abscissae and their
!          corresponding weights are given.
!
!          xgk    - abscissae of the 21-point Kronrod rule
!                   xgk(2), xgk(4), ...  abscissae of the 10-point
!                   Gauss rule
!                   xgk(1), xgk(3), ...  abscissae which are optimally
!                   added to the 10-point Gauss rule
!
!          wgk    - weights of the 21-point Kronrod rule
!
!          wg     - weights of the 10-point Gauss rule

  xgk = (/ 9.956571630258081e-01_dp,     9.739065285171717e-01_dp, &
           9.301574913557082e-01_dp,     8.650633666889845e-01_dp, &
           7.808177265864169e-01_dp,     6.794095682990244e-01_dp, &
           5.627571346686047e-01_dp,     4.333953941292472e-01_dp, &
           2.943928627014602e-01_dp,     1.488743389816312e-01_dp, &
           0.000000000000000_dp                               /)

  wgk = (/ 1.169463886737187e-02_dp,     3.255816230796473e-02_dp, &
           5.475589657435200e-02_dp,     7.503967481091995e-02_dp, &
           9.312545458369761e-02_dp,     1.093871588022976e-01_dp, &
           1.234919762620659e-01_dp,     1.347092173114733e-01_dp, &
           1.427759385770601e-01_dp,     1.477391049013385e-01_dp, &
           1.494455540029169e-01_dp                               /)

  wg  = (/ 6.667134430868814e-02_dp,     1.494513491505806e-01_dp, &
           2.190863625159820e-01_dp,     2.692667193099964e-01_dp, &
           2.955242247147529e-01_dp                               /)

!          list of major variables
!
!          centr  - mid point of the interval
!          hlgth  - half-length of the interval
!          absc   - abscissa
!          fval*  - function value
!          resg   - result of the 10-point Gauss formula
!          resk   - result of the 21-point Kronrod formula
!          reskh  - approximation to the mean value of f over (a,b),
!                   i.e. to i/(b-a)

  centr  = 5.0e-01_dp*(a+b)
  hlgth  = 5.0e-01_dp*(b-a)
  dhlgth = abs(hlgth)

! Compute the 21-point Kronrod approximation to the
! integral, and estimate the absolute error.
  resg = 0.0_dp
  fc = f(centr)
  resk = wgk(11)*fc
  resabs = abs(resk)

  do j = 1, 5
     jtw = 2*j
     absc = hlgth*xgk(jtw)
     fval1 = f(centr-absc)
     fval2 = f(centr+absc)
     fv1(jtw) = fval1
     fv2(jtw) = fval2
     fsum = fval1+fval2
     resg = resg+wg(j)*fsum
     resk = resk+wgk(jtw)*fsum
     resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
  enddo

  do j = 1, 5
     jtwm1 = 2*j-1
     absc = hlgth*xgk(jtwm1)
     fval1 = f(centr-absc)
     fval2 = f(centr+absc)
     fv1(jtwm1) = fval1
     fv2(jtwm1) = fval2
     fsum = fval1+fval2
     resk = resk+wgk(jtwm1)*fsum
     resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
  enddo

  reskh = resk*5.0e-01_dp
  resasc = wgk(11)*abs(fc-reskh)

  do j = 1, 10
     resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
  enddo

  result = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr = abs((resk-resg)*hlgth)

  if ( resasc /= 0.0_dp.and.abserr /= 0.0_dp) &
     abserr = resasc*min ( 1.0_dp,(2.0e+02_dp*abserr/resasc)**1.5_dp)

  if ( resabs > tiny ( resabs ) /(5.0e+01_dp* epsilon ( resabs ) )) &
     abserr = max (( epsilon ( resabs ) *5.0e+01_dp)*resabs,abserr)

end subroutine qk21

!===============================================================================
!###############################################################################
!===============================================================================
subroutine qk31(f,a,b,result,abserr,resabs,resasc)
!
!******************************************************************************
!
!! QK31 carries out a 31 point Gauss-Kronrod quadrature rule.
!
!
! Discussion:
!
!   This routine approximates
!     I = integral ( A <= X <= B ) F(X) dx
!   with an error estimate, and
!     J = integral ( A <= X <= B ) | F(X) | dx
!
! Reference:
!
!   R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
!   QUADPACK, a Subroutine Package for Automatic Integration,
!   Springer Verlag, 1983
!
! Parameters:
!
!   Input, external real F, the name of the function routine, of the form
!     function f ( x )
!     real f
!     real x
!   which evaluates the integrand function.
!
!   Input, real A, B, the limits of integration.
!
!   Output, real RESULT, the estimated value of the integral.
!                      result is computed by applying the 31-point
!                      Gauss-Kronrod rule (resk), obtained by optimal
!                      addition of abscissae to the 15-point Gauss
!                      rule (resg).
!
!   Output, real ABSERR, an estimate of | I - RESULT |.
!
!   Output, real RESABS, approximation to the integral of the absolute
!   value of F.
!
!   Output, real RESASC, approximation to the integral | F-I/(B-A) |
!   over [A,B].

!  REAL(dp), external    :: f
  interface
    function f(x)
      use kinds
      implicit none
      real(dp), intent(in) :: x
      real(dp)             :: f
    end function f
  end interface
  real(dp), intent(in)  :: a,b
  real(dp), intent(out) :: result,abserr,resabs,resasc

! Declare local variables
  integer  :: j,jtw,jtwm1
  real(dp) :: absc,centr,dhlgth,fc,fsum,fval1,fval2,fv1(15),fv2(15),hlgth, &
              resg,resk,reskh,wg(8),wgk(16),xgk(16)

!          the abscissae and weights are given for the interval (-1,1).
!          because of symmetry only the positive abscissae and their
!          corresponding weights are given.
!
!          xgk    - abscissae of the 31-point Kronrod rule
!                   xgk(2), xgk(4), ...  abscissae of the 15-point
!                   Gauss rule
!                   xgk(1), xgk(3), ...  abscissae which are optimally
!                   added to the 15-point Gauss rule
!
!          wgk    - weights of the 31-point Kronrod rule
!
!          wg     - weights of the 15-point Gauss rule

  xgk = (/ 9.980022986933971e-01_dp,   9.879925180204854e-01_dp, &
           9.677390756791391e-01_dp,   9.372733924007059e-01_dp, &
           8.972645323440819e-01_dp,   8.482065834104272e-01_dp, &
           7.904185014424659e-01_dp,   7.244177313601700e-01_dp, &
           6.509967412974170e-01_dp,   5.709721726085388e-01_dp, &
           4.850818636402397e-01_dp,   3.941513470775634e-01_dp, &
           2.991800071531688e-01_dp,   2.011940939974345e-01_dp, &
           1.011420669187175e-01_dp,   0.0_dp                  /)

  wgk = (/ 5.377479872923349e-03_dp,   1.500794732931612e-02_dp, &
           2.546084732671532e-02_dp,   3.534636079137585e-02_dp, &
           4.458975132476488e-02_dp,   5.348152469092809e-02_dp, &
           6.200956780067064e-02_dp,   6.985412131872826e-02_dp, &
           7.684968075772038e-02_dp,   8.308050282313302e-02_dp, &
           8.856444305621177e-02_dp,   9.312659817082532e-02_dp, &
           9.664272698362368e-02_dp,   9.917359872179196e-02_dp, &
           1.007698455238756e-01_dp,   1.013300070147915e-01_dp    /)

  wg  = (/ 3.075324199611727e-02_dp,   7.036604748810812e-02_dp, &
           1.071592204671719e-01_dp,   1.395706779261543e-01_dp, &
           1.662692058169939e-01_dp,   1.861610000155622e-01_dp, &
           1.984314853271116e-01_dp,   2.025782419255613e-01_dp    /)

!          list of major variables
!
!          centr  - mid point of the interval
!          hlgth  - half-length of the interval
!          absc   - abscissa
!          fval*  - function value
!          resg   - result of the 15-point Gauss formula
!          resk   - result of the 31-point Kronrod formula
!          reskh  - approximation to the mean value of f over (a,b),
!                   i.e. to i/(b-a)

  centr  = 5.0e-01_dp*(a+b)
  hlgth  = 5.0e-01_dp*(b-a)
  dhlgth = abs(hlgth)

! Compute the 31-point Kronrod approximation to the integral,
! and estimate the absolute error.
  fc = f(centr)
  resg = wg(8)*fc
  resk = wgk(16)*fc
  resabs = abs(resk)

  do j = 1, 7
     jtw = j*2
     absc = hlgth*xgk(jtw)
     fval1 = f(centr-absc)
     fval2 = f(centr+absc)
     fv1(jtw) = fval1
     fv2(jtw) = fval2
     fsum = fval1+fval2
     resg = resg+wg(j)*fsum
     resk = resk+wgk(jtw)*fsum
     resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
  enddo

  do j = 1, 8
     jtwm1 = j*2-1
     absc = hlgth*xgk(jtwm1)
     fval1 = f(centr-absc)
     fval2 = f(centr+absc)
     fv1(jtwm1) = fval1
     fv2(jtwm1) = fval2
     fsum = fval1+fval2
     resk = resk+wgk(jtwm1)*fsum
     resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
  enddo

  reskh = resk*5.0e-01_dp
  resasc = wgk(16)*abs(fc-reskh)

  do j = 1, 15
     resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
 enddo

  result = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr = abs((resk-resg)*hlgth)

  if ( resasc /= 0.0_dp.and.abserr /= 0.0_dp) &
       abserr = resasc*min ( 1.0_dp,(2.0e+02_dp*abserr/resasc)**1.5_dp)

  if ( resabs > tiny ( resabs ) /(5.0e+01_dp* epsilon ( resabs ) )) &
     abserr = max (( epsilon ( resabs ) *5.0e+01_dp)*resabs,abserr)

end subroutine qk31

!===============================================================================
!###############################################################################
!===============================================================================
subroutine qk41(f,a,b,result,abserr,resabs,resasc)
!
!******************************************************************************
!
!! QK41 carries out a 41 point Gauss-Kronrod quadrature rule.
!
!
! Discussion:
!
!   This routine approximates
!     I = integral ( A <= X <= B ) F(X) dx
!   with an error estimate, and
!     J = integral ( A <= X <= B ) | F(X) | dx
!
! Reference:
!
!   R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
!   QUADPACK, a Subroutine Package for Automatic Integration,
!   Springer Verlag, 1983
!
! Parameters:
!
!   Input, external real F, the name of the function routine, of the form
!     function f ( x )
!     real f
!     real x
!   which evaluates the integrand function.
!
!   Input, real A, B, the limits of integration.
!
!   Output, real RESULT, the estimated value of the integral.
!                      result is computed by applying the 41-point
!                      Gauss-Kronrod rule (resk) obtained by optimal
!                      addition of abscissae to the 20-point Gauss
!                      rule (resg).
!
!   Output, real ABSERR, an estimate of | I - RESULT |.
!
!   Output, real RESABS, approximation to the integral of the absolute
!   value of F.
!
!   Output, real RESASC, approximation to the integral | F-I/(B-A) |
!   over [A,B].
!
! Local Parameters:
!
!          centr  - mid point of the interval
!          hlgth  - half-length of the interval
!          absc   - abscissa
!          fval*  - function value
!          resg   - result of the 20-point Gauss formula
!          resk   - result of the 41-point Kronrod formula
!          reskh  - approximation to mean value of f over (a,b), i.e.
!                   to i/(b-a)

!  REAL(dp), external    :: f
  interface
    function f(x)
      use kinds
      implicit none
      real(dp), intent(in) :: x
      real(dp)             :: f
    end function f
  end interface
  real(dp), intent(in)  :: a,b
  real(dp), intent(out) :: result,abserr,resabs,resasc

! Declare local variables
  integer  :: j,jtw,jtwm1
  real(dp) :: absc,centr,dhlgth,fc,fsum,fval1,fval2,fv1(20),fv2(20),hlgth, &
              resg,resk,reskh,wg(10),wgk(21),xgk(21)
!
!          the abscissae and weights are given for the interval (-1,1).
!          because of symmetry only the positive abscissae and their
!          corresponding weights are given.
!
!          xgk    - abscissae of the 41-point Gauss-Kronrod rule
!                   xgk(2), xgk(4), ...  abscissae of the 20-point
!                   Gauss rule
!                   xgk(1), xgk(3), ...  abscissae which are optimally
!                   added to the 20-point Gauss rule
!
!          wgk    - weights of the 41-point Gauss-Kronrod rule
!
!          wg     - weights of the 20-point Gauss rule

  xgk = (/ 9.988590315882777e-01_dp,   9.931285991850949e-01_dp, &
           9.815078774502503e-01_dp,   9.639719272779138e-01_dp, &
           9.408226338317548e-01_dp,   9.122344282513259e-01_dp, &
           8.782768112522820e-01_dp,   8.391169718222188e-01_dp, &
           7.950414288375512e-01_dp,   7.463319064601508e-01_dp, &
           6.932376563347514e-01_dp,   6.360536807265150e-01_dp, &
           5.751404468197103e-01_dp,   5.108670019508271e-01_dp, &
           4.435931752387251e-01_dp,   3.737060887154196e-01_dp, &
           3.016278681149130e-01_dp,   2.277858511416451e-01_dp, &
           1.526054652409227e-01_dp,   7.652652113349733e-02_dp, &
           0.0_dp                                           /)

  wgk = (/ 3.073583718520532e-03_dp,   8.600269855642942e-03_dp, &
           1.462616925697125e-02_dp,   2.038837346126652e-02_dp, &
           2.588213360495116e-02_dp,   3.128730677703280e-02_dp, &
           3.660016975820080e-02_dp,   4.166887332797369e-02_dp, &
           4.643482186749767e-02_dp,   5.094457392372869e-02_dp, &
           5.519510534828599e-02_dp,   5.911140088063957e-02_dp, &
           6.265323755478117e-02_dp,   6.583459713361842e-02_dp, &
           6.864867292852162e-02_dp,   7.105442355344407e-02_dp, &
           7.303069033278667e-02_dp,   7.458287540049919e-02_dp, &
           7.570449768455667e-02_dp,   7.637786767208074e-02_dp, &
           7.660071191799966e-02_dp                             /)

  wg  = (/ 1.761400713915212e-02_dp,   4.060142980038694e-02_dp, &
           6.267204833410906e-02_dp,   8.327674157670475e-02_dp, &
           1.019301198172404e-01_dp,   1.181945319615184e-01_dp, &
           1.316886384491766e-01_dp,   1.420961093183821e-01_dp, &
           1.491729864726037e-01_dp,   1.527533871307259e-01_dp    /)

  centr  = 5.0e-01_dp*(a+b)
  hlgth  = 5.0e-01_dp*(b-a)
  dhlgth = abs(hlgth)

! Compute 41-point Gauss-Kronrod approximation to the
! the integral, and estimate the absolute error.
  resg = 0.0_dp
  fc = f(centr)
  resk = wgk(21)*fc
  resabs = abs(resk)

  do j = 1, 10
     jtw = j*2
     absc = hlgth*xgk(jtw)
     fval1 = f(centr-absc)
     fval2 = f(centr+absc)
     fv1(jtw) = fval1
     fv2(jtw) = fval2
     fsum = fval1+fval2
     resg = resg+wg(j)*fsum
     resk = resk+wgk(jtw)*fsum
     resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
  enddo

  do j = 1, 10
     jtwm1 = j*2-1
     absc = hlgth*xgk(jtwm1)
     fval1 = f(centr-absc)
     fval2 = f(centr+absc)
     fv1(jtwm1) = fval1
     fv2(jtwm1) = fval2
     fsum = fval1+fval2
     resk = resk+wgk(jtwm1)*fsum
     resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
  enddo

  reskh = resk*5.0e-01_dp
  resasc = wgk(21)*abs(fc-reskh)

  do j = 1, 20
     resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
  enddo

  result = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr = abs((resk-resg)*hlgth)

  if ( resasc /= 0.0_dp.and.abserr /= 0.0_dp) &
       abserr = resasc*min ( 1.0_dp,(2.0e+02_dp*abserr/resasc)**1.5_dp)

  if ( resabs > tiny ( resabs ) /(5.0e+01_dp* epsilon ( resabs ) )) &
     abserr = max (( epsilon ( resabs ) *5.0e+01_dp)*resabs,abserr)

end subroutine qk41

!===============================================================================
!###############################################################################
!===============================================================================
subroutine qk51(f,a,b,result,abserr,resabs,resasc)
!
!******************************************************************************
!
!! QK51 carries out a 51 point Gauss-Kronrod quadrature rule.
!
!
! Discussion:
!
!   This routine approximates
!     I = integral ( A <= X <= B ) F(X) dx
!   with an error estimate, and
!     J = integral ( A <= X <= B ) | F(X) | dx
!
! Reference:
!
!   R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
!   QUADPACK, a Subroutine Package for Automatic Integration,
!   Springer Verlag, 1983
!
! Parameters:
!
!   Input, external real F, the name of the function routine, of the form
!     function f ( x )
!     real f
!     real x
!   which evaluates the integrand function.
!
!   Input, real A, B, the limits of integration.
!
!   Output, real RESULT, the estimated value of the integral.
!                      result is computed by applying the 51-point
!                      Kronrod rule (resk) obtained by optimal addition
!                      of abscissae to the 25-point Gauss rule (resg).
!
!   Output, real ABSERR, an estimate of | I - RESULT |.
!
!   Output, real RESABS, approximation to the integral of the absolute
!   value of F.
!
!   Output, real RESASC, approximation to the integral | F-I/(B-A) |
!   over [A,B].
!
! Local Parameters:
!
!          centr  - mid point of the interval
!          hlgth  - half-length of the interval
!          absc   - abscissa
!          fval*  - function value
!          resg   - result of the 25-point Gauss formula
!          resk   - result of the 51-point Kronrod formula
!          reskh  - approximation to the mean value of f over (a,b),
!                   i.e. to i/(b-a)

!  REAL(dp), external    :: f
  interface
    function f(x)
      use kinds
      implicit none
      real(dp), intent(in) :: x
      real(dp)             :: f
    end function f
  end interface
  real(dp), intent(in)  :: a,b
  real(dp), intent(out) :: result,abserr,resabs,resasc

! Declare local variables
  integer  :: j,jtw,jtwm1
  real(dp) :: absc,centr,dhlgth,fc,fsum,fval1,fval2,fv1(25),fv2(25),hlgth, &
              resg,resk,reskh,wg(13),wgk(26),xgk(26)

!          the abscissae and weights are given for the interval (-1,1).
!          because of symmetry only the positive abscissae and their
!          corresponding weights are given.
!
!          xgk    - abscissae of the 51-point Kronrod rule
!                   xgk(2), xgk(4), ...  abscissae of the 25-point
!                   Gauss rule
!                   xgk(1), xgk(3), ...  abscissae which are optimally
!                   added to the 25-point Gauss rule
!
!          wgk    - weights of the 51-point Kronrod rule
!
!          wg     - weights of the 25-point Gauss rule

  xgk = (/ 9.992621049926098e-01_dp,   9.955569697904981e-01_dp, &
           9.880357945340772e-01_dp,   9.766639214595175e-01_dp, &
           9.616149864258425e-01_dp,   9.429745712289743e-01_dp, &
           9.207471152817016e-01_dp,   8.949919978782754e-01_dp, &
           8.658470652932756e-01_dp,   8.334426287608340e-01_dp, &
           7.978737979985001e-01_dp,   7.592592630373576e-01_dp, &
           7.177664068130844e-01_dp,   6.735663684734684e-01_dp, &
           6.268100990103174e-01_dp,   5.776629302412230e-01_dp, &
           5.263252843347192e-01_dp,   4.730027314457150e-01_dp, &
           4.178853821930377e-01_dp,   3.611723058093878e-01_dp, &
           3.030895389311078e-01_dp,   2.438668837209884e-01_dp, &
           1.837189394210489e-01_dp,   1.228646926107104e-01_dp, &
           6.154448300568508e-02_dp,   0.0_dp                  /)

  wgk = (/ 1.987383892330316e-03_dp,   5.561932135356714e-03_dp, &
           9.473973386174152e-03_dp,   1.323622919557167e-02_dp, &
           1.684781770912830e-02_dp,   2.043537114588284e-02_dp, &
           2.400994560695322e-02_dp,   2.747531758785174e-02_dp, &
           3.079230016738749e-02_dp,   3.400213027432934e-02_dp, &
           3.711627148341554e-02_dp,   4.008382550403238e-02_dp, &
           4.287284502017005e-02_dp,   4.550291304992179e-02_dp, &
           4.798253713883671e-02_dp,   5.027767908071567e-02_dp, &
           5.236288580640748e-02_dp,   5.425112988854549e-02_dp, &
           5.595081122041232e-02_dp,   5.743711636156783e-02_dp, &
           5.868968002239421e-02_dp,   5.972034032417406e-02_dp, &
           6.053945537604586e-02_dp,   6.112850971705305e-02_dp, &
           6.147118987142532e-02_dp,   6.158081806783294e-02_dp    /)

  wg  = (/ 1.139379850102629e-02_dp,   2.635498661503214e-02_dp, &
           4.093915670130631e-02_dp,   5.490469597583519e-02_dp, &
           6.803833381235692e-02_dp,   8.014070033500102e-02_dp, &
           9.102826198296365e-02_dp,   1.005359490670506e-01_dp, &
           1.085196244742637e-01_dp,   1.148582591457116e-01_dp, &
           1.194557635357848e-01_dp,   1.222424429903100e-01_dp, &
           1.231760537267155e-01_dp                             /)

  centr  = 5.0e-01_dp*(a+b)
  hlgth  = 5.0e-01_dp*(b-a)
  dhlgth = abs(hlgth)

! Compute the 51-point Kronrod approximation to the integral,
! and estimate the absolute error.
  fc = f(centr)
  resg = wg(13)*fc
  resk = wgk(26)*fc
  resabs = abs(resk)

  do j = 1, 12
     jtw = j*2
     absc = hlgth*xgk(jtw)
     fval1 = f(centr-absc)
     fval2 = f(centr+absc)
     fv1(jtw) = fval1
     fv2(jtw) = fval2
     fsum = fval1+fval2
     resg = resg+wg(j)*fsum
     resk = resk+wgk(jtw)*fsum
     resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
  enddo

  do j = 1, 13
     jtwm1 = j*2-1
     absc = hlgth*xgk(jtwm1)
     fval1 = f(centr-absc)
     fval2 = f(centr+absc)
     fv1(jtwm1) = fval1
     fv2(jtwm1) = fval2
     fsum = fval1+fval2
     resk = resk+wgk(jtwm1)*fsum
     resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
  enddo

  reskh = resk*5.0e-01_dp
  resasc = wgk(26)*abs(fc-reskh)

  do j = 1, 25
     resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
  enddo

  result = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr = abs((resk-resg)*hlgth)

  if ( resasc /= 0.0_dp.and.abserr /= 0.0_dp) &
     abserr = resasc*min ( 1.0_dp,(2.0e+02_dp*abserr/resasc)**1.5_dp)

  if ( resabs > tiny ( resabs ) / (5.0e+01_dp* epsilon ( resabs ) ) ) &
     abserr = max (( epsilon ( resabs ) *5.0e+01_dp)*resabs,abserr)

end subroutine qk51

!===============================================================================
!###############################################################################
!===============================================================================
subroutine qk61(f,a,b,result,abserr,resabs,resasc)
!
!******************************************************************************
!
!! QK61 carries out a 61 point Gauss-Kronrod quadrature rule.
!
!
! Discussion:
!
!   This routine approximates
!     I = integral ( A <= X <= B ) F(X) dx
!   with an error estimate, and
!     J = integral ( A <= X <= B ) | F(X) | dx
!
! Reference:
!
!   R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
!   QUADPACK, a Subroutine Package for Automatic Integration,
!   Springer Verlag, 1983
!
! Parameters:
!
!   Input, external real F, the name of the function routine, of the form
!     function f ( x )
!     real f
!     real x
!   which evaluates the integrand function.
!
!   Input, real A, B, the limits of integration.
!
!   Output, real RESULT, the estimated value of the integral.
!                   result is computed by applying the 61-point
!                   Kronrod rule (resk) obtained by optimal addition of
!                   abscissae to the 30-point Gauss rule (resg).
!
!   Output, real ABSERR, an estimate of | I - RESULT |.
!
!   Output, real RESABS, approximation to the integral of the absolute
!   value of F.
!
!   Output, real RESASC, approximation to the integral | F-I/(B-A) |
!   over [A,B].
!
! Local Parameters:
!
!          centr  - mid point of the interval
!          hlgth  - half-length of the interval
!          absc   - abscissa
!          fval*  - function value
!          resg   - result of the 30-point Gauss rule
!          resk   - result of the 61-point Kronrod rule
!          reskh  - approximation to the mean value of f
!                   over (a,b), i.e. to i/(b-a)

!  REAL(dp), external    :: f
  interface
    function f(x)
      use kinds
      implicit none
      real(dp), intent(in) :: x
      real(dp)             :: f
    end function f
  end interface
  real(dp), intent(in)  :: a,b
  real(dp), intent(out) :: result,abserr,resabs,resasc

! Declare local variables
  integer  :: j,jtw,jtwm1
  real(dp) :: absc,centr,dhlgth,fc,fsum,fval1,fval2,fv1(30),fv2(30),hlgth, &
              resg,resk,reskh,wg(15),wgk(31),xgk(31)

!          the abscissae and weights are given for the
!          interval (-1,1). because of symmetry only the positive
!          abscissae and their corresponding weights are given.
!
!          xgk   - abscissae of the 61-point Kronrod rule
!                  xgk(2), xgk(4)  ... abscissae of the 30-point
!                  Gauss rule
!                  xgk(1), xgk(3)  ... optimally added abscissae
!                  to the 30-point Gauss rule
!
!          wgk   - weights of the 61-point Kronrod rule
!
!          wg    - weigths of the 30-point Gauss rule

  xgk = (/ 9.994844100504906e-01_dp,     9.968934840746495e-01_dp, &
           9.916309968704046e-01_dp,     9.836681232797472e-01_dp, &
           9.731163225011263e-01_dp,     9.600218649683075e-01_dp, &
           9.443744447485600e-01_dp,     9.262000474292743e-01_dp, &
           9.055733076999078e-01_dp,     8.825605357920527e-01_dp, &
           8.572052335460611e-01_dp,     8.295657623827684e-01_dp, &
           7.997278358218391e-01_dp,     7.677774321048262e-01_dp, &
           7.337900624532268e-01_dp,     6.978504947933158e-01_dp, &
           6.600610641266270e-01_dp,     6.205261829892429e-01_dp, &
           5.793452358263617e-01_dp,     5.366241481420199e-01_dp, &
           4.924804678617786e-01_dp,     4.470337695380892e-01_dp, &
           4.004012548303944e-01_dp,     3.527047255308781e-01_dp, &
           3.040732022736251e-01_dp,     2.546369261678898e-01_dp, &
           2.045251166823099e-01_dp,     1.538699136085835e-01_dp, &
           1.028069379667370e-01_dp,     5.147184255531770e-02_dp, &
           0.0_dp                                             /)

  wgk = (/ 1.389013698677008e-03_dp,     3.890461127099884e-03_dp, &
           6.630703915931292e-03_dp,     9.273279659517763e-03_dp, &
           1.182301525349634e-02_dp,     1.436972950704580e-02_dp, &
           1.692088918905327e-02_dp,     1.941414119394238e-02_dp, &
           2.182803582160919e-02_dp,     2.419116207808060e-02_dp, &
           2.650995488233310e-02_dp,     2.875404876504129e-02_dp, &
           3.090725756238776e-02_dp,     3.298144705748373e-02_dp, &
           3.497933802806002e-02_dp,     3.688236465182123e-02_dp, &
           3.867894562472759e-02_dp,     4.037453895153596e-02_dp, &
           4.196981021516425e-02_dp,     4.345253970135607e-02_dp, &
           4.481480013316266e-02_dp,     4.605923827100699e-02_dp, &
           4.718554656929915e-02_dp,     4.818586175708713e-02_dp, &
           4.905543455502978e-02_dp,     4.979568342707421e-02_dp, &
           5.040592140278235e-02_dp,     5.088179589874961e-02_dp, &
           5.122154784925877e-02_dp,     5.142612853745903e-02_dp, &
           5.149472942945157e-02_dp                               /)

  wg  = (/ 7.968192496166606e-03_dp,     1.846646831109096e-02_dp, &
           2.878470788332337e-02_dp,     3.879919256962705e-02_dp, &
           4.840267283059405e-02_dp,     5.749315621761907e-02_dp, &
           6.597422988218050e-02_dp,     7.375597473770521e-02_dp, &
           8.075589522942022e-02_dp,     8.689978720108298e-02_dp, &
           9.212252223778613e-02_dp,     9.636873717464426e-02_dp, &
           9.959342058679527e-02_dp,     1.017623897484055e-01_dp, &
           1.028526528935588e-01_dp                               /)

  centr = 5.0e-01_dp*(b+a)
  hlgth = 5.0e-01_dp*(b-a)
  dhlgth = abs(hlgth)

! Compute the 61-point Kronrod approximation to the integral,
! and estimate the absolute error.
  resg = 0.0_dp
  fc = f(centr)
  resk = wgk(31)*fc
  resabs = abs(resk)

  do j = 1, 15
     jtw = j*2
     absc = hlgth*xgk(jtw)
     fval1 = f(centr-absc)
     fval2 = f(centr+absc)
     fv1(jtw) = fval1
     fv2(jtw) = fval2
     fsum = fval1+fval2
     resg = resg+wg(j)*fsum
     resk = resk+wgk(jtw)*fsum
     resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
  enddo

  do j = 1, 15
     jtwm1 = j*2-1
     absc = hlgth*xgk(jtwm1)
     fval1 = f(centr-absc)
     fval2 = f(centr+absc)
     fv1(jtwm1) = fval1
     fv2(jtwm1) = fval2
     fsum = fval1+fval2
     resk = resk+wgk(jtwm1)*fsum
     resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
  enddo

  reskh  = resk * 5.0e-01_dp
  resasc = wgk(31)*abs(fc-reskh)

  do j = 1, 30
     resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
  enddo

  result = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr = abs((resk-resg)*hlgth)

  if ( resasc /= 0.0_dp .and. abserr /= 0.0_dp) &
     abserr = resasc*min ( 1.0_dp,(2.0e+02_dp*abserr/resasc)**1.5_dp)

  if ( resabs > tiny ( resabs ) / (5.0e+01_dp* epsilon ( resabs ) )) &
     abserr = max ( ( epsilon ( resabs ) *5.0e+01_dp)*resabs, abserr )

end subroutine qk61

!===============================================================================
!###############################################################################
!===============================================================================
subroutine qmomo(alfa,beta,ri,rj,rg,rh,integr)
!
!******************************************************************************
!
!! QMOMO computes modified Chebyshev moments.
!
!
! Discussion:
!
!   This routine computes modified Chebyshev moments.
!   The K-th modified Chebyshev moment is defined as the
!   integral over (-1,1) of W(X)*T(K,X), where T(K,X) is the
!   Chebyshev polynomial of degree K.
!
! Reference:
!
!   R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
!   QUADPACK, a Subroutine Package for Automatic Integration,
!   Springer Verlag, 1983
!
! Parameters:
!
!   Input, real ALFA, a parameter in the weight function w(x), ALFA > -1.
!
!   Input, real BETA, a parameter in the weight function w(x), BETA > -1.
!
!          ri     - real
!                   vector of dimension 25
!                   ri(k) is the integral over (-1,1) of
!                   (1+x)**alfa*t(k-1,x), k = 1, ..., 25.
!
!          rj     - real
!                   vector of dimension 25
!                   rj(k) is the integral over (-1,1) of
!                   (1-x)**beta*t(k-1,x), k = 1, ..., 25.
!
!          rg     - real
!                   vector of dimension 25
!                   rg(k) is the integral over (-1,1) of
!                   (1+x)**alfa*log((1+x)/2)*t(k-1,x), k = 1, ...,25.
!
!          rh     - real
!                   vector of dimension 25
!                   rh(k) is the integral over (-1,1) of
!                   (1-x)**beta*log((1-x)/2)*t(k-1,x), k = 1, ..., 25.
!
!          integr - integer
!                   input parameter indicating the modified moments
!                   to be computed
!                   integr = 1 compute ri, rj
!                          = 2 compute ri, rj, rg
!                          = 3 compute ri, rj, rh
!                          = 4 compute ri, rj, rg, rh

  integer,  intent(in)  :: integr
  real(dp), intent(in)  :: alfa,beta
  real(dp), intent(out) :: rg(25),rh(25),ri(25),rj(25)

! Declare local variables
  integer  :: i,im1
  real(dp) :: alfp1,alfp2,an,anm1,betp1,betp2,ralf,rbet

  alfp1 = alfa+1.0_dp
  betp1 = beta+1.0_dp
  alfp2 = alfa+2.0_dp
  betp2 = beta+2.0_dp
  ralf  = 2.0_dp**alfp1
  rbet  = 2.0_dp**betp1

! Compute RI, RJ using a forward recurrence relation.
  ri(1) = ralf/alfp1
  rj(1) = rbet/betp1
  ri(2) = ri(1)*alfa/alfp2
  rj(2) = rj(1)*beta/betp2
  an    = 2.0_dp
  anm1  = 1.0_dp

  do i = 3, 25
     ri(i) = -(ralf+an*(an-alfp2)*ri(i-1))/(anm1*(an+alfp1))
     rj(i) = -(rbet+an*(an-betp2)*rj(i-1))/(anm1*(an+betp1))
     anm1 = an
     an = an+1.0_dp
  enddo

  if ( integr == 1 ) go to 70
  if ( integr == 3 ) go to 40

! Compute RG using a forward recurrence relation.
  rg(1) = -ri(1)/alfp1
  rg(2) = -(ralf+ralf)/(alfp2*alfp2)-rg(1)
  an    = 2.0_dp
  anm1  = 1.0_dp
  im1   = 2

  do i = 3, 25
     rg(i) = -(an*(an-alfp2)*rg(im1)-an*ri(im1)+anm1*ri(i))/ &
          (anm1*(an+alfp1))
     anm1 = an
     an   = an+1.0_dp
     im1  = i
  enddo

  if ( integr == 2 ) go to 70

! Compute RH using a forward recurrence relation.
40 continue

  rh(1) = -rj(1) / betp1
  rh(2) = -(rbet+rbet)/(betp2*betp2)-rh(1)
  an    = 2.0_dp
  anm1  = 1.0_dp
  im1   = 2

  do i = 3, 25
     rh(i) = -(an*(an-betp2)*rh(im1)-an*rj(im1)+ &
          anm1*rj(i))/(anm1*(an+betp1))
     anm1 = an
     an = an+1.0_dp
     im1 = i
  enddo

  do i = 2, 25, 2
     rh(i) = -rh(i)
  enddo

   70 continue

  do i = 2, 25, 2
     rj(i) = -rj(i)
  enddo

end subroutine qmomo

!===============================================================================
!###############################################################################
!===============================================================================
subroutine qng(f,a,b,epsabs,epsrel,result,abserr,neval,ifail)
!
!******************************************************************************
!
!! QNG estimates an integral, using non-adaptive integration.
!
!
! Discussion:
!
!   The routine calculates an approximation RESULT to a definite integral
!     I = integral of F over (A,B),
!   hopefully satisfying
!     || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
!
!   The routine is a simple non-adaptive automatic integrator, based on
!   a sequence of rules with increasing degree of algebraic
!   precision (Patterson, 1968).
!
! Reference:
!
!   R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
!   QUADPACK, a Subroutine Package for Automatic Integration,
!   Springer Verlag, 1983
!
! Parameters:
!
!   Input, external real F, the name of the function routine, of the form
!     function f ( x )
!     real f
!     real x
!   which evaluates the integrand function.
!
!   Input, real A, B, the limits of integration.
!
!   Input, real EPSABS, EPSREL, the absolute and relative accuracy requested.
!
!   Output, real RESULT, the estimated value of the integral.
!   RESULT is obtained by applying the 21-point Gauss-Kronrod rule (RES21)
!   obtained  by optimal addition of abscissae to the 10-point Gauss rule
!   (RES10), or by applying the 43-point rule (RES43) obtained by optimal
!   addition of abscissae to the 21-point Gauss-Kronrod rule, or by
!   applying the 87-point rule (RES87) obtained by optimal addition of
!   abscissae to the 43-point rule.
!
!   Output, real ABSERR, an estimate of || I - RESULT ||.
!
!   Output, integer NEVAL, the number of times the integral was evaluated.
!
!      ifail    - ifail = 0 normal and reliable termination of the
!                           routine. it is assumed that the requested
!                           accuracy has been achieved.
!                 ifail > 0 abnormal termination of the routine. it is
!                           assumed that the requested accuracy has
!                           not been achieved.
!                 ifail = 1 the maximum number of steps has been
!                           executed. the integral is probably too
!                           difficult to be calculated by qng.
!                       = 6 the input is invalid, because
!                           epsabs < 0 and epsrel < 0,
!                           result, abserr and neval are set to zero.
!
! Local Parameters:
!
!          centr  - mid point of the integration interval
!          hlgth  - half-length of the integration interval
!          fcentr - function value at mid point
!          absc   - abscissa
!          fval   - function value
!          savfun - array of function values which have already
!                   been computed
!          res10  - 10-point Gauss result
!          res21  - 21-point Kronrod result
!          res43  - 43-point result
!          res87  - 87-point result
!          resabs - approximation to the integral of abs(f)
!          resasc - approximation to the integral of abs(f-i/(b-a))

!  REAL(dp), external    :: f
  interface
    function f(x)
      use kinds
      implicit none
      real(dp), intent(in) :: x
      real(dp)             :: f
    end function f
  end interface
  real(dp), intent(in)  :: a,b,epsabs,epsrel
  integer,  intent(out) :: neval,ifail
  real(dp), intent(out) :: result,abserr

! Declare local variables
  integer  :: ipx,k,l
  real(dp) :: absc,centr,dhlgth,fcentr,fval,fval1,fval2,fv1(5),fv2(5),fv3(5),fv4(5), &
              hlgth,res10,res21,res43,res87,resabs,resasc,reskh,savfun(21),w10(5),   &
              w21a(5),w21b(6),w43a(10),w43b(12),w87a(21),w87b(23),x1(5),x2(5),x3(11),x4(22)

!          the following data statements contain the abscissae
!          and weights of the integration rules used.
!
!          x1      abscissae common to the 10-, 21-, 43- and 87-point
!                  rule
!          x2      abscissae common to the 21-, 43- and 87-point rule
!          x3      abscissae common to the 43- and 87-point rule
!          x4      abscissae of the 87-point rule
!          w10     weights of the 10-point formula
!          w21a    weights of the 21-point formula for abscissae x1
!          w21b    weights of the 21-point formula for abscissae x2
!          w43a    weights of the 43-point formula for absissae x1, x3
!          w43b    weights of the 43-point formula for abscissae x3
!          w87a    weights of the 87-point formula for abscissae x1,
!                  x2 and x3
!          w87b    weights of the 87-point formula for abscissae x4

  x1   = (/ 9.739065285171717e-01_dp,     8.650633666889845e-01_dp, &
            6.794095682990244e-01_dp,     4.333953941292472e-01_dp, &
            1.488743389816312e-01_dp                               /)

  x2   = (/ 9.956571630258081e-01_dp,     9.301574913557082e-01_dp, &
            7.808177265864169e-01_dp,     5.627571346686047e-01_dp, &
            2.943928627014602e-01_dp                               /)

  x3   = (/ 9.993333609019321e-01_dp,     9.874334029080889e-01_dp, &
            9.548079348142663e-01_dp,     9.001486957483283e-01_dp, &
            8.251983149831142e-01_dp,     7.321483889893050e-01_dp, &
            6.228479705377252e-01_dp,     4.994795740710565e-01_dp, &
            3.649016613465808e-01_dp,     2.222549197766013e-01_dp, &
            7.465061746138332e-02_dp                               /)

  x4   = (/ 9.999029772627292e-01_dp,     9.979898959866787e-01_dp, &
            9.921754978606872e-01_dp,     9.813581635727128e-01_dp, &
            9.650576238583846e-01_dp,     9.431676131336706e-01_dp, &
            9.158064146855072e-01_dp,     8.832216577713165e-01_dp, &
            8.457107484624157e-01_dp,     8.035576580352310e-01_dp, &
            7.570057306854956e-01_dp,     7.062732097873218e-01_dp, &
            6.515894665011779e-01_dp,     5.932233740579611e-01_dp, &
            5.314936059708319e-01_dp,     4.667636230420228e-01_dp, &
            3.994248478592188e-01_dp,     3.298748771061883e-01_dp, &
            2.585035592021616e-01_dp,     1.856953965683467e-01_dp, &
            1.118422131799075e-01_dp,     3.735212339461987e-02_dp    /)

  w10  = (/ 6.667134430868814e-02_dp,     1.494513491505806e-01_dp, &
            2.190863625159820e-01_dp,     2.692667193099964e-01_dp, &
            2.955242247147529e-01_dp                               /)

  w21a = (/ 3.255816230796473e-02_dp,     7.503967481091995e-02_dp, &
            1.093871588022976e-01_dp,     1.347092173114733e-01_dp, &
            1.477391049013385e-01_dp                               /)

  w21b = (/ 1.169463886737187e-02_dp,     5.475589657435200e-02_dp, &
            9.312545458369761e-02_dp,     1.234919762620659e-01_dp, &
            1.427759385770601e-01_dp,     1.494455540029169e-01_dp    /)

  w43a = (/ 1.629673428966656e-02_dp,     3.752287612086950e-02_dp, &
            5.469490205825544e-02_dp,     6.735541460947809e-02_dp, &
            7.387019963239395e-02_dp,     5.768556059769796e-03_dp, &
            2.737189059324884e-02_dp,     4.656082691042883e-02_dp, &
            6.174499520144256e-02_dp,     7.138726726869340e-02_dp    /)

  w43b = (/ 1.844477640212414e-03_dp,     1.079868958589165e-02_dp, &
            2.189536386779543e-02_dp,     3.259746397534569e-02_dp, &
            4.216313793519181e-02_dp,     5.074193960018458e-02_dp, &
            5.837939554261925e-02_dp,     6.474640495144589e-02_dp, &
            6.956619791235648e-02_dp,     7.282444147183321e-02_dp, &
            7.450775101417512e-02_dp,     7.472214751740301e-02_dp    /)

  w87a = (/ 8.148377384149173e-03_dp,     1.876143820156282e-02_dp, &
            2.734745105005229e-02_dp,     3.367770731163793e-02_dp, &
            3.693509982042791e-02_dp,     2.884872430211531e-03_dp, &
            1.368594602271270e-02_dp,     2.328041350288831e-02_dp, &
            3.087249761171336e-02_dp,     3.569363363941877e-02_dp, &
            9.152833452022414e-04_dp,     5.399280219300471e-03_dp, &
            1.094767960111893e-02_dp,     1.629873169678734e-02_dp, &
            2.108156888920384e-02_dp,     2.537096976925383e-02_dp, &
            2.918969775647575e-02_dp,     3.237320246720279e-02_dp, &
            3.478309895036514e-02_dp,     3.641222073135179e-02_dp, &
            3.725387550304771e-02_dp                               /)

  w87b = (/ 2.741455637620724e-04_dp,     1.807124155057943e-03_dp, &
            4.096869282759165e-03_dp,     6.758290051847379e-03_dp, &
            9.549957672201647e-03_dp,     1.232944765224485e-02_dp, &
            1.501044734638895e-02_dp,     1.754896798624319e-02_dp, &
            1.993803778644089e-02_dp,     2.219493596101229e-02_dp, &
            2.433914712600081e-02_dp,     2.637450541483921e-02_dp, &
            2.828691078877120e-02_dp,     3.005258112809270e-02_dp, &
            3.164675137143993e-02_dp,     3.305041341997850e-02_dp, &
            3.425509970422606e-02_dp,     3.526241266015668e-02_dp, &
            3.607698962288870e-02_dp,     3.669860449845609e-02_dp, &
            3.712054926983258e-02_dp,     3.733422875193504e-02_dp, &
            3.736107376267902e-02_dp                               /)

! Test on validity of parameters.
  result = 0.0_dp
  abserr = 0.0_dp
  neval  = 0

  if ( epsabs < 0.0_dp .and. epsrel < 0.0_dp ) then
     ifail = 6
     return
  endif

  hlgth  = 5.0e-01_dp*(b-a)
  dhlgth = abs(hlgth)
  centr  = 5.0e-01_dp*(b+a)
  fcentr = f(centr)
  neval  = 21
  ifail  = 1

! Compute the integral using the 10- and 21-point formula.
  do l = 1, 3
    if ( l == 1 ) then
       res10  = 0.0_dp
       res21  = w21b(6) * fcentr
       resabs = w21b(6) * abs(fcentr)
       do k = 1, 5
          absc = hlgth*x1(k)
          fval1 = f(centr+absc)
          fval2 = f(centr-absc)
          fval = fval1+fval2
          res10 = res10+w10(k)*fval
          res21 = res21+w21a(k)*fval
          resabs = resabs+w21a(k)*(abs(fval1)+abs(fval2))
          savfun(k) = fval
          fv1(k) = fval1
          fv2(k) = fval2
       enddo

       ipx = 5

       do k = 1, 5
          ipx = ipx+1
          absc = hlgth*x2(k)
          fval1 = f(centr+absc)
          fval2 = f(centr-absc)
          fval = fval1 + fval2
          res21 = res21 + w21b(k) * fval
          resabs = resabs + w21b(k) * (abs(fval1)+abs(fval2))
          savfun(ipx) = fval
          fv3(k) = fval1
          fv4(k) = fval2
       enddo

! Test for convergence.
       result = res21*hlgth
       resabs = resabs*dhlgth
       reskh = 5.0e-01_dp*res21
       resasc = w21b(6)*abs(fcentr-reskh)

       do k = 1, 5
          resasc = resasc+w21a(k)*(abs(fv1(k)-reskh)+abs(fv2(k)-reskh)) &
               +w21b(k)*(abs(fv3(k)-reskh)+abs(fv4(k)-reskh))
       enddo

       abserr = abs((res21-res10)*hlgth)
       resasc = resasc*dhlgth

! Compute the integral using the 43-point formula.
    elseif ( l == 2 ) then

       res43 = w43b(12)*fcentr
       neval = 43

       do k = 1, 10
          res43 = res43+savfun(k) * w43a(k)
       enddo

      do k = 1, 11
         ipx = ipx+1
         absc = hlgth*x3(k)
         fval = f(absc+centr)+f(centr-absc)
         res43 = res43+fval*w43b(k)
         savfun(ipx) = fval
      enddo

! Test for convergence.
      result = res43 * hlgth
      abserr = abs((res43-res21)*hlgth)

! Compute the integral using the 87-point formula.
   elseif ( l == 3 ) then
      res87 = w87b(23) * fcentr
      neval = 87

      do k = 1, 21
         res87 = res87 + savfun(k) * w87a(k)
      enddo

      do k = 1, 22
         absc = hlgth * x4(k)
         res87 = res87+w87b(k)*(f(absc+centr)+f(centr-absc))
      enddo

      result = res87 * hlgth
      abserr = abs ( ( res87-res43) * hlgth )
    endif

    if ( resasc /= 0.0_dp.and.abserr /= 0.0_dp ) &
       abserr = resasc * min ( 1.0_dp,(2.0e+02_dp*abserr/resasc)**1.5_dp)

    if ( resabs > tiny ( resabs ) / ( 5.0e+01_dp * epsilon ( resabs ) ) ) &
       abserr = max (( epsilon ( resabs ) *5.0e+01_dp) * resabs, abserr )

    if ( abserr <= max ( epsabs, epsrel*abs(result))) ifail = 0

    if ( ifail == 0 ) exit
  enddo

end subroutine qng

!===============================================================================
!###############################################################################
!===============================================================================
      subroutine quad(f,a,b,epsabs,answer,neval,ifail)

!     This quadrature program uses formulae due to T. N. L. Patterson,
!     Mathematics of computation, Volume 22, 1968, pages 847-856, as
!     modified by F. T. Krogh and W. V. Snyder, ACM Transactions on
!     Mathematical Software 17, 4 (December 1991) pp 457-461.  It is a
!     functional replacement for Algorithm 468, T. N. L. Patterson,
!     Communications of the ACM 16, 11 (November 1973) 694-699.

!     *****     Formal Arguments     ***********************************

! INPUT:
!  f       A function subprogram that evaluates the integrand at a given
!          abscissa.  F is invoked F(X).  F must be mentioned in an
!          E X T E R N A L statement in the calling program.
!  a, b    Lower and upper limits of integration, respectively.
!  epsabs  Relative accuracy required.  When the relative difference of
!          two successive formulae does not exceed EPSABS the last formula
!          computed is taken as the result.
! OUTPUT:
!  answer  This is the k-th element of the array RESULT and is the
!          approximation to the integral.
!  neval   Reports the number of integrand evaluations.
!  ifail   On exit normally IFAIL=0.  However if convergence to the
!          accuracy requested is not achieved IFAIL=1 on exit.

!      REAL(dp), external    :: F
      interface
        function f(x)
          use kinds
          implicit none
          real(dp), intent(in) :: x
          real(dp)             :: f
        end function f
      end interface
      real(dp), intent(in)  :: A, B, epsabs
      real(dp), intent(out) :: answer
      integer,  intent(out) :: neval, ifail

!     *****     Local Variables     ************************************

! ACUM    is the accumulating estimate of the integral.
! DELTA   is B - A.
! DIFF    is 0.5 * DELTA.
! FH, FL  contain subscripts to indicate where to store integrand
!         samples in WORK.
! FNCVAL  is an integrand sample, or a sum of two symmetrically placed
!         integrand samples.
! IP      is a subscript used to index P.
! J       is a subscript and loop inductor.
! J1, J2  are bounds of the indexes at which integrand samples are
!         stored in WORK.
! JH, JL  are bounds for the loop that accumulates new function values
!         into the integral estimate.
! KH, KL  are bounds for indexes at which integrand samples are
!         retrieved from WORK in order to begin applying a quadrature
!         formula.
! KK      is a loop inductor and subscript.
! KX      is a list of bounds of subscripts into KH and KL.  KH is
!         indexed by K.  The set of indexes from which to retrieve
!         integrand samples is given by the set of bounds in KL and KH
!         indexed by KX(K-1)+1 to KX(K) inclusively.
! P       contains the coefficients necessary to the quadrature
!         formulae.  Their organization is described below in the
!         section on DATA statements.
! PACUM   is the previous value of ACUM.
! X       is the distance of the abscissa from the boundary of the
!         region .
! WORK    is space in which to store function values to be used in the
!         next formula.

      integer  :: fh(2:8), fl(2:8), kh(11), kk, kl(11), kx(1:8)
      integer  :: k, ip, j, j1, j2, jh, jl
      real(dp) :: acum, delta, diff, fncval, pacum, work(17), x, p(305)
      real(dp) :: result(8)

!     *****     Data Statements     ************************************
      FL = (/2, 3, 5, 9,12,14, 1/)
      FH = (/2, 4, 8,16,17,17, 0/)
      KL = (/1, 1, 1, 1, 1, 3, 5, 9, 5, 9,12/)
      KH = (/1, 2, 4, 8,16, 3, 6,17, 5, 9,17/)
      KX = (/0, 1, 2, 3, 4, 5, 8, 11/)

!     In the comments below, F(K,I) refers to the function value
!     computed for the I'th node of the K'th formula.  The abscissae and
!     weights are stored in order according to the distance from the
!     boundary of the region, not from the center.  Since we store
!     1 - |abscissa|, the first "node" coefficient for each formula is
!     the smallest.

!     Corrections, nodes and weights for the 3-point formula.

!     Correction for F(1,1).
       P(  1) = -0.11111111111111111111e+00_dp
!     Node and weight for F(2,1).
       P(  2) = +0.22540333075851662296e+00_dp
       P(  3) = +0.55555555555555555556e+00_dp

!     Corrections, nodes and weights for the 7-point formula.

!     Corrections for F(1,1) and F(2,1).
       P(  4) = +0.64720942140296979100e-02_dp
       P(  5) = -0.92896879094443370500e-02_dp
!     Nodes and weights for F(3,1-2)
       P(  6) = +0.39508731291979716579e-01_dp
       P(  7) = +0.10465622602646726519e+00_dp
       P(  8) = +0.56575625065319744200e+00_dp
       P(  9) = +0.40139741477596222291e+00_dp

!     Corrections, nodes and weights for the 15-point formula.

!     Corrections for F(1,1), F(2,1), F(3,1-2).
       P( 10) = +0.52230468969616220000e-04_dp
       P( 11) = +0.17121030961750000000e-03_dp
       P( 12) = -0.72483001615389289800e-03_dp
       P( 13) = -0.70178010992090420000e-04_dp
!     Nodes and weights for F(4,1-4).
       P( 14) = +0.61680367872449777899e-02_dp
       P( 15) = +0.17001719629940260339e-01_dp
       P( 16) = +0.11154076712774300110e+00_dp
       P( 17) = +0.92927195315124537686e-01_dp
       P( 18) = +0.37889705326277359705e+00_dp
       P( 19) = +0.17151190913639138079e+00_dp
       P( 20) = +0.77661331357103311837e+00_dp
       P( 21) = +0.21915685840158749640e+00_dp

!     Corrections, nodes and weights for the 31-point formula.

!     Corrections for F(1,1), F(2,1), F(3,1-2), F(4,1-4).
       P( 22) = +0.68216653479200000000e-08_dp
       P( 23) = +0.12667409859336000000e-06_dp
       P( 24) = +0.59565976367837165000e-05_dp
       P( 25) = +0.13923301068260000000e-07_dp
       P( 26) = -0.66294075649023920000e-04_dp
       P( 27) = -0.70439580428230200000e-06_dp
       P( 28) = -0.34518205339241000000e-07_dp
       P( 29) = -0.81448691099600000000e-08_dp
!     Nodes and weights for F(5,1-8).
       P( 30) = +0.90187503233240234038e-03_dp
       P( 31) = +0.25447807915618744154e-02_dp
       P( 32) = +0.18468850446259893130e-01_dp
       P( 33) = +0.16446049854387810934e-01_dp
       P( 34) = +0.70345142570259943330e-01_dp
       P( 35) = +0.35957103307129322097e-01_dp
       P( 36) = +0.16327406183113126449e+00_dp
       P( 37) = +0.56979509494123357412e-01_dp
       P( 38) = +0.29750379350847292139e+00_dp
       P( 39) = +0.76879620499003531043e-01_dp
       P( 40) = +0.46868025635562437602e+00_dp
       P( 41) = +0.93627109981264473617e-01_dp
       P( 42) = +0.66886460674202316691e+00_dp
       P( 43) = +0.10566989358023480974e+00_dp
       P( 44) = +0.88751105686681337425e+00_dp
       P( 45) = +0.11195687302095345688e+00_dp

!     Corrections, nodes and weights for the 63-point formula.

!     Corrections for F(1,1), F(2,1), F(3,1-2), F(4,1-4), F(5,1-8).
       P( 46) = +0.37158300000000000000e-15_dp
       P( 47) = +0.21237877000000000000e-12_dp
       P( 48) = +0.10522629388435000000e-08_dp
       P( 49) = +0.17480290000000000000e-14_dp
       P( 50) = +0.34757189830171600000e-06_dp
       P( 51) = +0.90312761725000000000e-11_dp
       P( 52) = +0.12558916000000000000e-13_dp
       P( 53) = +0.54591000000000000000e-15_dp
       P( 54) = -0.72338395508691963000e-05_dp
       P( 55) = -0.16969957975797700000e-07_dp
       P( 56) = -0.85436390715500000000e-10_dp
       P( 57) = -0.12281300930000000000e-11_dp
       P( 58) = -0.46233482500000000000e-13_dp
       P( 59) = -0.42244055000000000000e-14_dp
       P( 60) = -0.88501000000000000000e-15_dp
       P( 61) = -0.40904000000000000000e-15_dp
!     Nodes and weights for F(6,1-16).
       P( 62) = +0.12711187964238806027e-03_dp
       P( 63) = +0.36322148184553065969e-03_dp
       P( 64) = +0.27937406277780409196e-02_dp
       P( 65) = +0.25790497946856882724e-02_dp
       P( 66) = +0.11315242452570520059e-01_dp
       P( 67) = +0.61155068221172463397e-02_dp
       P( 68) = +0.27817125251418203419e-01_dp
       P( 69) = +0.10498246909621321898e-01_dp
       P( 70) = +0.53657141626597094849e-01_dp
       P( 71) = +0.15406750466559497802e-01_dp
       P( 72) = +0.89628843042995707499e-01_dp
       P( 73) = +0.20594233915912711149e-01_dp
       P( 74) = +0.13609206180630952284e+00_dp
       P( 75) = +0.25869679327214746911e-01_dp
       P( 76) = +0.19305946804978238813e+00_dp
       P( 77) = +0.31073551111687964880e-01_dp
       P( 78) = +0.26024395564730524132e+00_dp
       P( 79) = +0.36064432780782572640e-01_dp
       P( 80) = +0.33709033997521940454e+00_dp
       P( 81) = +0.40715510116944318934e-01_dp
       P( 82) = +0.42280428994795418516e+00_dp
       P( 83) = +0.44914531653632197414e-01_dp
       P( 84) = +0.51638197305415897244e+00_dp
       P( 85) = +0.48564330406673198716e-01_dp
       P( 86) = +0.61664067580126965307e+00_dp
       P( 87) = +0.51583253952048458777e-01_dp
       P( 88) = +0.72225017797817568492e+00_dp
       P( 89) = +0.53905499335266063927e-01_dp
       P( 90) = +0.83176474844779253501e+00_dp
       P( 91) = +0.55481404356559363988e-01_dp
       P( 92) = +0.94365568695340721002e+00_dp
       P( 93) = +0.56277699831254301273e-01_dp

!     Corrections, nodes and weights for the 127-point formula.

!     Corrections for F(3,1), F(4,1-2), F(5,1-3), F(6,1-6).
       P( 94) = +0.10410980000000000000e-15_dp
       P( 95) = +0.24947205459800000000e-10_dp
       P( 96) = +0.55000000000000000000e-20_dp
       P( 97) = +0.29041247599538500000e-07_dp
       P( 98) = +0.36728212600000000000e-13_dp
       P( 99) = +0.55680000000000000000e-18_dp
       P(100) = -0.87117647737697202500e-06_dp
       P(101) = -0.81473242674410000000e-09_dp
       P(102) = -0.88309203370000000000e-12_dp
       P(103) = -0.18018239000000000000e-14_dp
       P(104) = -0.70528000000000000000e-17_dp
       P(105) = -0.50600000000000000000e-19_dp
!     Nodes and weights for F(7,1-32).
       P(106) = +0.17569645108401419961e-04_dp
       P(107) = +0.50536095207862517625e-04_dp
       P(108) = +0.40120032808931675009e-03_dp
       P(109) = +0.37774664632698466027e-03_dp
       P(110) = +0.16833646815926074696e-02_dp
       P(111) = +0.93836984854238150079e-03_dp
       P(112) = +0.42758953015928114900e-02_dp
       P(113) = +0.16811428654214699063e-02_dp
       P(114) = +0.85042788218938676006e-02_dp
       P(115) = +0.25687649437940203731e-02_dp
       P(116) = +0.14628500401479628890e-01_dp
       P(117) = +0.35728927835172996494e-02_dp
       P(118) = +0.22858485360294285840e-01_dp
       P(119) = +0.46710503721143217474e-02_dp
       P(120) = +0.33362148441583432910e-01_dp
       P(121) = +0.58434498758356395076e-02_dp
       P(122) = +0.46269993574238863589e-01_dp
       P(123) = +0.70724899954335554680e-02_dp
       P(124) = +0.61679602220407116350e-01_dp
       P(125) = +0.83428387539681577056e-02_dp
       P(126) = +0.79659974529987579270e-01_dp
       P(127) = +0.96411777297025366953e-02_dp
       P(128) = +0.10025510022305996335e+00_dp
       P(129) = +0.10955733387837901648e-01_dp
       P(130) = +0.12348658551529473026e+00_dp
       P(131) = +0.12275830560082770087e-01_dp
       P(132) = +0.14935550523164972024e+00_dp
       P(133) = +0.13591571009765546790e-01_dp
       P(134) = +0.17784374563501959262e+00_dp
       P(135) = +0.14893641664815182035e-01_dp
       P(136) = +0.20891506620015163857e+00_dp
       P(137) = +0.16173218729577719942e-01_dp
       P(138) = +0.24251603361948636206e+00_dp
       P(139) = +0.17421930159464173747e-01_dp
       P(140) = +0.27857691462990108452e+00_dp
       P(141) = +0.18631848256138790186e-01_dp
       P(142) = +0.31701256890892077191e+00_dp
       P(143) = +0.19795495048097499488e-01_dp
       P(144) = +0.35772335749024048622e+00_dp
       P(145) = +0.20905851445812023852e-01_dp
       P(146) = +0.40059606975775710702e+00_dp
       P(147) = +0.21956366305317824939e-01_dp
       P(148) = +0.44550486736806745112e+00_dp
       P(149) = +0.22940964229387748761e-01_dp
       P(150) = +0.49231224246628339785e+00_dp
       P(151) = +0.23854052106038540080e-01_dp
       P(152) = +0.54086998801016766712e+00_dp
       P(153) = +0.24690524744487676909e-01_dp
       P(154) = +0.59102017877011132759e+00_dp
       P(155) = +0.25445769965464765813e-01_dp
       P(156) = +0.64259616216846784762e+00_dp
       P(157) = +0.26115673376706097680e-01_dp
       P(158) = +0.69542355844328595666e+00_dp
       P(159) = +0.26696622927450359906e-01_dp
       P(160) = +0.74932126969651682339e+00_dp
       P(161) = +0.27185513229624791819e-01_dp
       P(162) = +0.80410249728889984607e+00_dp
       P(163) = +0.27579749566481873035e-01_dp
       P(164) = +0.85957576684743982540e+00_dp
       P(165) = +0.27877251476613701609e-01_dp
       P(166) = +0.91554595991628911629e+00_dp
       P(167) = +0.28076455793817246607e-01_dp
       P(168) = +0.97181535105025430566e+00_dp
       P(169) = +0.28176319033016602131e-01_dp

!     Corrections, nodes and weights for the 255-point formula.

!     Corrections for F(4,1), F(5,1), F(6,1-2), F(7,1-4).
       P(170) = +0.33260000000000000000e-18_dp
       P(171) = +0.11409477047800000000e-11_dp
       P(172) = +0.29524360569703510000e-08_dp
       P(173) = +0.51608328000000000000e-15_dp
       P(174) = -0.11017721965059732300e-06_dp
       P(175) = -0.58656987416475000000e-10_dp
       P(176) = -0.23340340645000000000e-13_dp
       P(177) = -0.12489500000000000000e-16_dp
!     Nodes and weights for F(8,1-64).
       P(178) = +0.24036202515353807630e-05_dp
       P(179) = +0.69379364324108267170e-05_dp
       P(180) = +0.56003792945624240417e-04_dp
       P(181) = +0.53275293669780613125e-04_dp
       P(182) = +0.23950907556795267013e-03_dp
       P(183) = +0.13575491094922871973e-03_dp
       P(184) = +0.61966197497641806982e-03_dp
       P(185) = +0.24921240048299729402e-03_dp
       P(186) = +0.12543855319048853002e-02_dp
       P(187) = +0.38974528447328229322e-03_dp
       P(188) = +0.21946455040427254399e-02_dp
       P(189) = +0.55429531493037471492e-03_dp
       P(190) = +0.34858540851097261500e-02_dp
       P(191) = +0.74028280424450333046e-03_dp
       P(192) = +0.51684971993789994803e-02_dp
       P(193) = +0.94536151685852538246e-03_dp
       P(194) = +0.72786557172113846706e-02_dp
       P(195) = +0.11674841174299594077e-02_dp
       P(196) = +0.98486295992298408193e-02_dp
       P(197) = +0.14049079956551446427e-02_dp
       P(198) = +0.12907472045965932809e-01_dp
       P(199) = +0.16561127281544526052e-02_dp
       P(200) = +0.16481342421367271240e-01_dp
       P(201) = +0.19197129710138724125e-02_dp
       P(202) = +0.20593718329137316189e-01_dp
       P(203) = +0.21944069253638388388e-02_dp
       P(204) = +0.25265540247597332240e-01_dp
       P(205) = +0.24789582266575679307e-02_dp
       P(206) = +0.30515340497540768229e-01_dp
       P(207) = +0.27721957645934509940e-02_dp
       P(208) = +0.36359378430187867480e-01_dp
       P(209) = +0.30730184347025783234e-02_dp
       P(210) = +0.42811783890139037259e-01_dp
       P(211) = +0.33803979910869203823e-02_dp
       P(212) = +0.49884702478705123440e-01_dp
       P(213) = +0.36933779170256508183e-02_dp
       P(214) = +0.57588434808916940190e-01_dp
       P(215) = +0.40110687240750233989e-02_dp
       P(216) = +0.65931563842274211999e-01_dp
       P(217) = +0.43326409680929828545e-02_dp
       P(218) = +0.74921067092924347640e-01_dp
       P(219) = +0.46573172997568547773e-02_dp
       P(220) = +0.84562412844234959360e-01_dp
       P(221) = +0.49843645647655386012e-02_dp
       P(222) = +0.94859641186738404810e-01_dp
       P(223) = +0.53130866051870565663e-02_dp
       P(224) = +0.10581543166444097714e+00_dp
       P(225) = +0.56428181013844441585e-02_dp
       P(226) = +0.11743115975265809315e+00_dp
       P(227) = +0.59729195655081658049e-02_dp
       P(228) = +0.12970694445188609414e+00_dp
       P(229) = +0.63027734490857587172e-02_dp
       P(230) = +0.14264168911376784347e+00_dp
       P(231) = +0.66317812429018878941e-02_dp
       P(232) = +0.15623311732729139895e+00_dp
       P(233) = +0.69593614093904229394e-02_dp
       P(234) = +0.17047780536259859981e+00_dp
       P(235) = +0.72849479805538070639e-02_dp
       P(236) = +0.18537121234486258656e+00_dp
       P(237) = +0.76079896657190565832e-02_dp
       P(238) = +0.20090770903915859819e+00_dp
       P(239) = +0.79279493342948491103e-02_dp
       P(240) = +0.21708060588171698360e+00_dp
       P(241) = +0.82443037630328680306e-02_dp
       P(242) = +0.23388218069623990928e+00_dp
       P(243) = +0.85565435613076896192e-02_dp
       P(244) = +0.25130370638306339718e+00_dp
       P(245) = +0.88641732094824942641e-02_dp
       P(246) = +0.26933547875781873867e+00_dp
       P(247) = +0.91667111635607884067e-02_dp
       P(248) = +0.28796684463774796540e+00_dp
       P(249) = +0.94636899938300652943e-02_dp
       P(250) = +0.30718623022088529711e+00_dp
       P(251) = +0.97546565363174114611e-02_dp
       P(252) = +0.32698116976958152079e+00_dp
       P(253) = +0.10039172044056840798e-01_dp
       P(254) = +0.34733833458998250389e+00_dp
       P(255) = +0.10316812330947621682e-01_dp
       P(256) = +0.36824356228880576959e+00_dp
       P(257) = +0.10587167904885197931e-01_dp
       P(258) = +0.38968188628481359983e+00_dp
       P(259) = +0.10849844089337314099e-01_dp
       P(260) = +0.41163756555233745857e+00_dp
       P(261) = +0.11104461134006926537e-01_dp
       P(262) = +0.43409411457634557737e+00_dp
       P(263) = +0.11350654315980596602e-01_dp
       P(264) = +0.45703433350168850951e+00_dp
       P(265) = +0.11588074033043952568e-01_dp
       P(266) = +0.48044033846254297801e+00_dp
       P(267) = +0.11816385890830235763e-01_dp
       P(268) = +0.50429359208123853983e+00_dp
       P(269) = +0.12035270785279562630e-01_dp
       P(270) = +0.52857493412834112307e+00_dp
       P(271) = +0.12244424981611985899e-01_dp
       P(272) = +0.55326461233797152625e+00_dp
       P(273) = +0.12443560190714035263e-01_dp
       P(274) = +0.57834231337383669993e+00_dp
       P(275) = +0.12632403643542078765e-01_dp
       P(276) = +0.60378719394238406082e+00_dp
       P(277) = +0.12810698163877361967e-01_dp
       P(278) = +0.62957791204992176986e+00_dp
       P(279) = +0.12978202239537399286e-01_dp
       P(280) = +0.65569265840056197721e+00_dp
       P(281) = +0.13134690091960152836e-01_dp
       P(282) = +0.68210918793152331682e+00_dp
       P(283) = +0.13279951743930530650e-01_dp
       P(284) = +0.70880485148175331803e+00_dp
       P(285) = +0.13413793085110098513e-01_dp
       P(286) = +0.73575662758907323806e+00_dp
       P(287) = +0.13536035934956213614e-01_dp
       P(288) = +0.76294115441017027278e+00_dp
       P(289) = +0.13646518102571291428e-01_dp
       P(290) = +0.79033476175681880523e+00_dp
       P(291) = +0.13745093443001896632e-01_dp
       P(292) = +0.81791350324074780175e+00_dp
       P(293) = +0.13831631909506428676e-01_dp
       P(294) = +0.84565318851862189130e+00_dp
       P(295) = +0.13906019601325461264e-01_dp
       P(296) = +0.87352941562769803314e+00_dp
       P(297) = +0.13968158806516938516e-01_dp
       P(298) = +0.90151760340188079791e+00_dp
       P(299) = +0.14017968039456608810e-01_dp
       P(300) = +0.92959302395714482093e+00_dp
       P(301) = +0.14055382072649964277e-01_dp
       P(302) = +0.95773083523463639678e+00_dp
       P(303) = +0.14080351962553661325e-01_dp
       P(304) = +0.98590611358921753738e+00_dp
       P(305) = +0.14092845069160408355e-01_dp

!     *****     Executable Statements     ******************************
      ifail = 0
      delta  = b-a
      diff   = 0.5_dp*delta
      ip = 1
      jh = 0

!     Apply 1-point Gauss formula (Midpoint rule).

      fncval=f(a+diff)
!     Don't write "0.5*(b+a)" above if the radix of arithmetic isn't 2.
      neval=1
      work(1)=fncval
      acum=fncval*delta
      result(1)=acum

      do k = 2,8

!       Go on to the next formula.
         pacum=acum
         acum=0.0

!       Compute contribution to current estimate due to function
!       values used in previous formulae.
         do kk = kx(k-1)+1, kx(k)
            do j = kl(kk), kh(kk)
               acum=acum+real(p(ip)*work(j),dp)
               ip=ip+1
            enddo
         enddo

!       Compute contribution from new function values.
         jl=jh+1
         jh=jl+jl-1
         j1=fl(k)
         j2=fh(k)
         do j = jl, jh
            x=p(ip)*diff
            fncval=f(a+x)+f(b-x)
            neval=neval+2
            acum=acum+real(p(ip+1)*fncval,dp)
            if (j1.le.j2) then
               work(j1)=fncval
               j1=j1+1
            end if
            ip=ip+2
         enddo
         acum=real(diff,dp)*acum+0.5_dp*pacum
         result(k)=acum
         if (abs(result(k)-result(k-1)).le.abs(epsabs*result(k))) then
            answer = result(k)
            return
         endif
      enddo
      ifail=1
      k=8
      answer = result(k)
      end subroutine quad

!===============================================================================
!###############################################################################
!===============================================================================
subroutine qsort(limit,last,maxerr,ermax,elist,iord,nrmax)
!
!******************************************************************************
!
!! QSORT maintains the order of a list of local error estimates.
!
!
! Discussion:
!
!   This routine maintains the descending ordering in the list of the
!   local error estimates resulting from the interval subdivision process.
!   At each call two error estimates are inserted using the sequential
!   search top-down for the largest error estimate and bottom-up for the
!   smallest error estimate.
!
! Reference:
!
!   R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
!   QUADPACK, a Subroutine Package for Automatic Integration,
!   Springer Verlag, 1983
!
! Parameters:
!
!   Input, integer LIMIT, the maximum number of error estimates the list can
!   contain.
!
!   Input, integer LAST, the current number of error estimates.
!
!   Input/output, integer MAXERR, the index in the list of the NRMAX-th
!   largest error.
!
!   Output, real ERMAX, the NRMAX-th largest error = ELIST(MAXERR).
!
!   Input, real ELIST(LIMIT), contains the error estimates.
!
!   Input/output, integer IORD(LAST).  The first K elements contain
!   pointers to the error estimates such that ELIST(IORD(1)) through
!   ELIST(IORD(K)) form a decreasing sequence, with
!     K = LAST
!   if
!     LAST <= (LIMIT/2+2),
!   and otherwise
!     K = LIMIT+1-LAST.
!
!   Input/output, integer NRMAX.

  integer,  intent(in)    :: limit,last
  real(dp), intent(in)    :: elist(limit)
  integer,  intent(inout) :: maxerr,nrmax,iord(last)
  real(dp), intent(out)   :: ermax

! Declare local variables
  integer  :: i,ibeg,isucc,j,jbnd,jupbn,k
  real(dp) :: errmax,errmin

! Check whether the list contains more than two error estimates.
  if ( last <= 2 ) then
     iord(1) = 1
     iord(2) = 2
     go to 90
  endif

! This part of the routine is only executed if, due to a
! difficult integrand, subdivision increased the error
! estimate. in the normal case the insert procedure should
! start after the nrmax-th largest error estimate.
  errmax = elist(maxerr)

  do i = 1, nrmax-1
    isucc = iord(nrmax-1)
    if ( errmax <= elist(isucc) ) exit
    iord(nrmax) = isucc
    nrmax = nrmax-1
  enddo

! Compute the number of elements in the list to be maintained
! in descending order.  This number depends on the number of
! subdivisions still allowed.
  jupbn = last

  if ( last > (limit/2+2) ) jupbn = limit+3-last

  errmin = elist(last)

! Insert errmax by traversing the list top-down, starting
! comparison from the element elist(iord(nrmax+1)).
  jbnd = jupbn-1
  ibeg = nrmax+1

  do i = ibeg, jbnd
     isucc = iord(i)
     if ( errmax >= elist(isucc) ) go to 60
     iord(i-1) = isucc
  enddo

  iord(jbnd) = maxerr
  iord(jupbn) = last
  go to 90

! Insert errmin by traversing the list bottom-up.
60 continue

  iord(i-1) = maxerr
  k = jbnd

  do j = i, jbnd
     isucc = iord(k)
     if ( errmin < elist(isucc) ) go to 80
     iord(k+1) = isucc
     k = k-1
  enddo

  iord(i) = last
  go to 90

80 continue

  iord(k+1) = last

! Set maxerr and ermax.
90 continue

  maxerr = iord(nrmax)
  ermax  = elist(maxerr)

end subroutine qsort

!===============================================================================
!###############################################################################
!===============================================================================
function qwgtc(x,c,p2,p3,p4,kp)
!
!******************************************************************************
!
!! QWGTC defines the weight function used by QC25C.
!
!
! Discussion:
!
!   The weight function has the form 1 / ( X - C ).
!
! Reference:
!
!   R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
!   QUADPACK, a Subroutine Package for Automatic Integration,
!   Springer Verlag, 1983
!
! Parameters:
!
!   Input, real X, the point at which the weight function is evaluated.
!
!   Input, real C, the location of the singularity.
!
!   Input, real P2, P3, P4, parameters that are not used.
!
!   Input, integer KP, a parameter that is not used.
!
!   Output, real QWGTC, the value of the weight function at X.

  real(dp) :: qwgtc
  real(dp), intent(in) :: c,p2,p3,p4,x
  integer,  intent(in) :: kp

  qwgtc = 1.0_dp / ( x - c )

end function qwgtc

!===============================================================================
!###############################################################################
!===============================================================================
function qwgto(x,omega,p2,p3,p4,integr)
!
!******************************************************************************
!
!! QWGTO defines the weight functions used by QC25O.
!
!
! Reference:
!
!   R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
!   QUADPACK, a Subroutine Package for Automatic Integration,
!   Springer Verlag, 1983
!
! Parameters:
!
!   Input, real X, the point at which the weight function is evaluated.
!
!   Input, real OMEGA, the factor multiplying X.
!
!   Input, real P2, P3, P4, parameters that are not used.
!
!   Input, integer INTEGR, specifies which weight function is used:
!   1. W(X) = cos ( OMEGA * X )
!   2, W(X) = sin ( OMEGA * X )
!
!   Output, real QWGTO, the value of the weight function at X.

  real(dp) :: qwgto
  integer,  intent(in) :: integr
  real(dp), intent(in) :: omega,p2,p3,p4,x

  if ( integr == 1 ) then
     qwgto = cos ( omega * x )
  elseif ( integr == 2 ) then
     qwgto = sin ( omega * x )
  endif

end function  qwgto

!===============================================================================
!###############################################################################
!===============================================================================
function qwgts(x,a,b,alfa,beta,integr)
!
!******************************************************************************
!
!! QWGTS defines the weight functions used by QC25S.
!
!
! Reference:
!
!   R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
!   QUADPACK, a Subroutine Package for Automatic Integration,
!   Springer Verlag, 1983
!
! Parameters:
!
!   Input, real X, the point at which the weight function is evaluated.
!
!   Input, real A, B, the endpoints of the integration interval.
!
!   Input, real ALFA, BETA, exponents that occur in the weight function.
!
!   Input, integer INTEGR, specifies which weight function is used:
!   1. W(X) = (X-A)**ALFA * (B-X)**BETA
!   2, W(X) = (X-A)**ALFA * (B-X)**BETA * log (X-A)
!   3, W(X) = (X-A)**ALFA * (B-X)**BETA * log (B-X)
!   4, W(X) = (X-A)**ALFA * (B-X)**BETA * log (X-A) * log(B-X)
!
!   Output, real QWGTS, the value of the weight function at X.

  real(dp) :: qwgts
  integer,  intent(in) :: integr
  real(dp), intent(in) :: a,b,alfa,beta,x

  if ( integr == 1 ) then
     qwgts = ( x - a )**alfa * ( b - x )**beta
  elseif ( integr == 2 ) then
     qwgts = ( x - a )**alfa * ( b - x )**beta * log ( x - a )
  elseif ( integr == 3 ) then
     qwgts = ( x - a )**alfa * ( b - x )**beta * log ( b - x )
  elseif ( integr == 4 ) then
     qwgts = ( x - a )**alfa * ( b - x )**beta * log ( x - a ) * log ( b - x )
  endif

end function qwgts

!===============================================================================
!###############################################################################
!===============================================================================

subroutine r_swap(x,y)
!
!*******************************************************************************
!
!! R_SWAP swaps two real values.
!
!
! Modified:
!
!   01 May 2000
!
! Author:
!
!   John Burkardt
!
! Parameters:
!
!   Input/output, real X, Y.  On output, the values of X and
!                             Y have been interchanged.

  real(dp), intent(inout) :: x, y
  real(dp) :: z

  z = x
  x = y
  y = z

end subroutine r_swap

end module Quadpack

!===============================================================================
!###############################################################################
!===============================================================================


!---------------------------------------------------------------------!
!                                                                     !
!                      NumFort numerical library                      !
!                                                                     !
!---------------------------------------------------------------------!

module numFort
  use kinds
  use quadpack
  use LAPACK95
  use minFun
  implicit none

  interface linspace
     module procedure linspace,linspaceReal,linspaceInt
  end interface linspace
  interface integral
     module procedure integral, integralToInfty, integralOf, integralBreakPts
  end interface integral
  interface rk4
     module procedure rk4,rk4N
  end interface rk4
  interface GuessZero
     module procedure GuessZero,GuessZeroF
  end interface GuessZero
  interface Trace
     module procedure TraceDP,TraceSP,TraceComplexDP,TraceComplexSP
  end interface Trace
  interface inv
     module procedure invMatSP, invMatDP, invMatComplexSP, invMatComplexDP
  end interface inv
  interface pyplot
     module procedure pyplotN,pyplotXY,pyplotXYZW
  end interface pyplot
  interface Euler
     module procedure euler,eulerND
  end interface Euler
  interface writeData
     module procedure writeDataXY,writeDataXYZW,writeDataN
  end interface writeData
  interface trapZ
     module procedure trapZ,trapZZ
  end interface trapZ

contains

  !---------------------------------------------------------------------!
  !                                                                     !
  !                            FFT Algorithim                           !
  !                                                                     !
  !---------------------------------------------------------------------!

  ! In place Cooley-Tukey FFT
  recursive subroutine fft(x)
    complex(DP), dimension(:), intent(inout)  :: x
    complex(DP)                               :: t

    complex(DP), dimension(:), allocatable    :: even, odd
    integer :: N
    integer :: i

    N=size(x)

    if(N .le. 1) return

    allocate(odd((N+1)/2))
    allocate(even(N/2))

    ! divide
    odd  = x(1:N:2)
    even = x(2:N:2)

    ! conquer
    call fft(odd)
    call fft(even)

    ! combine
    do i=1,N/2
       t=exp(cmplx(0.0_DP,-2.0_DP*pi*real(i-1,DP)/real(N,DP),DP))*even(i)
       x(i)     = odd(i) + t
       x(i+N/2) = odd(i) - t
    end do

    deallocate(odd)
    deallocate(even)

  end subroutine fft

  !---------------------------------------------------------------------!
  !                                                                     !
  !                         Write data to a file                        !
  !                                                                     !
  !---------------------------------------------------------------------!

  subroutine writeDataXY(x,y,title)
    use Kinds
    implicit none
    real(DP),dimension(:),intent(in)     :: x,y
    character(len=*),intent(in),optional :: title
    integer :: ii,N

    N = size(x)
    if (present(title)) then
       open(110,file=title,action="write", &
            & status="replace",form="formatted")
    else
       open(110,file="data.dat",action="write", &
            & status="replace",form="formatted")
    end if
    do ii = 1,N
       write(110,'(2e20.10)') x(ii), y(ii)
    end do

    close(110)

  end subroutine writeDataXY

  subroutine writeDataXYZW(x,y,z,w,title)
    use Kinds
    implicit none
    real(DP),dimension(:),intent(in)     :: x,y,z,w
    character(len=*),intent(in),optional :: title
    integer :: ii,N

    N = size(x)
    if (present(title)) then
       open(110,file=title,action="write", &
            & status="replace",form="formatted")
    else
       open(110,file="data.dat",action="write", &
            & status="replace",form="formatted")
    end if
    do ii = 1,N
       write(110,'(4e20.10)') x(ii), y(ii), z(ii), w(ii)
    end do

    close(110)

  end subroutine writeDataXYZW

  subroutine writeDataN(x,title)
    use Kinds
    implicit none
    real(DP),dimension(:,:),intent(in)     :: x
    character(len=*),intent(in),optional   :: title
    character(len=14)                      :: fmt
    integer :: ii,N

    N = size(x,dim=2)
    if ( N < 10) then
       write(fmt,'(a1,i1,a8)') '(', N, 'es20.10)'
    else
       write(fmt,'(a1,i2,a8)') '(', N, 'es20.10)'
    end if

    if (present(title)) then
       open(110,file=title,action="write", &
            & status="replace",form="formatted")
    else
       open(110,file="data.dat",action="write", &
            & status="replace",form="formatted")
    end if
    do ii = 1,size(x,dim=1)
       write(110,fmt) x(ii,:)
    end do

    close(110)

  end subroutine writeDataN

  !---------------------------------------------------------------------!
  !                                                                     !
  !             evaluate Bessel function at x of order 0<n<3            !
  !                                                                     !
  !---------------------------------------------------------------------!

  function bessel(x,n) result(val)
    use kinds
    implicit none

    real(DP),intent(in) :: x
    integer ,intent(in) :: n
    real(DP)            :: val

    if( n > 5 .or. n < 0 ) then
       write(*,'(a)') "Error: n must be an integer between 0 and 5"
    else

       select case(n)
       case(0)
          val = sin(x)/x
       case(1)
          val = sin(x)/x**2-cos(x)/x
       case(2)
          val = (3/x**2-1)*sin(x)/x -3*cos(x)/x**2
       case(3)
          val = (15/x**3-6/x)*sin(x)/x-(15/x**2-1)*sin(x)/x
       case(4)
          val = (sin(x)+4*cos(x)/x-12*sin(x)/x**2-24*cos(x)/x**3+24*sin(x)/x**4)/x
       case(5)
          val = (cos(x)-5*sin(x)/x-20*cos(x)/x**2+60*sin(x)/x**3+120*cos(x)/x**4-120*sin(x)/x**5)/x
       end select
    end if

  end function bessel

  !---------------------------------------------------------------------!
  !                                                                     !
  !                   Calculate the trace of a Matrix                   !
  !                                                                     !
  !---------------------------------------------------------------------!

  function TraceDP(M)
    use kinds
    implicit none

    real(DP),dimension(:,:),intent(in) :: M
    real(DP)                           :: TraceDP

    integer :: rows,cols,i

    rows = size(M(:,1))
    cols = size(M(1,:))
    TraceDP = 0
    if ( rows .ne. cols ) then
       write(*,'(a)') "Error, not a square matrix"
    else
       do i = 1,rows
          TraceDP = TraceDP + M(i,i)
       end do
    end if

  end function TraceDP

  function TraceSP(M)
    use kinds
    implicit none

    real(SP),dimension(:,:),intent(in) :: M
    real(SP)                           :: TraceSP

    integer :: rows,cols,i

    rows = size(M(:,1))
    cols = size(M(1,:))
    TraceSP = 0
    if ( rows .ne. cols ) then
       write(*,'(a)') "Error, not a square matrix"
    else
       do i = 1,rows
          TraceSP = TraceSP + M(i,i)
       end do
    end if

  end function TraceSP

  function TraceComplexDP(M)
    use kinds
    implicit none

    complex(DP),dimension(:,:),intent(in) :: M
    complex(DP)                           :: TraceComplexDP

    integer :: rows,cols,i

    rows = size(M(:,1))
    cols = size(M(1,:))
    TraceComplexDP = 0.0_DP
    if ( rows .ne. cols ) then
       write(*,'(a)') "Error, not a square matrix"
    else
       do i = 1,rows
          TraceComplexDP = TraceComplexDP + M(i,i)
       end do
    end if

  end function TraceComplexDP

  function TraceComplexSP(M)
    use kinds
    implicit none

    complex(SP),dimension(:,:),intent(in) :: M
    complex(SP)                           :: TraceComplexSP

    integer :: rows,cols,i

    rows = size(M(:,1))
    cols = size(M(1,:))
    TraceComplexSP = 0.0_DP
    if ( rows .ne. cols ) then
       write(*,'(a)') "Error, not a square matrix"
    else
       do i = 1,rows
          TraceComplexSP = TraceComplexSP + M(i,i)
       end do
    end if

  end function TraceComplexSP

  !---------------------------------------------------------------------!
  !                                                                     !
  !                  Calculate the inverse of a matrix                  !
  !                                                                     !
  !---------------------------------------------------------------------!

  function invMatSP(A) result(B)
    real(SP), intent(in) :: A(:,:)
    real(SP)             :: B(size(A,1),size(A,2))
    ! Local variables
    integer              :: n

    n = size(A,1)
    select case (n)
    case (1)
       if (A(1,1).eq.0.0) then
          write(*,*) 'Matrix is numerically singular!'
       end if
       B(:,:) = 1.0_SP / A(:,:)
    case (2)
       B(:,:) = invMatSP2(A)
    case (3)
       B(:,:) = invMatSP3(A)
    case (4)
       B(:,:) = invMatSP4(A)
    case (5:)
       B(:,:) = invMatSPN(A)
    end select

  end function invMatSP

  function invMatDP(A) result(B)
    real(DP), intent(in) :: A(:,:)
    real(DP)             :: B(size(A,1),size(A,2))
    ! Local variables
    integer              :: n

    n = size(A,1)
    select case (n)
    case (1)
       if (A(1,1).eq.0.0_DP) then
          write(*,*) 'Matrix is numerically singular!'
       end if
       B(:,:) = 1.0_DP / A(:,:)
    case (2)
       B(:,:) = invMatDP2(A)
    case (3)
       B(:,:) = invMatDP3(A)
    case (4)
       B(:,:) = invMatDP4(A)
    case (5:)
       B(:,:) = invMatDPN(A)
    end select

  end function invMatDP

  function invMatComplexSP(A) result(B)
    complex(SP), intent(in) :: A(:,:)
    complex(SP)             :: B(size(A,1),size(A,2))
    ! Local variables
    integer              :: n

    n = size(A,1)
    select case (n)
    case (1)
       if (A(1,1).eq.cmplx(0.0,0.0)) then
          write(*,*) 'Matrix is numerically singular!'
       end if
       B(:,:) = 1.0 / A(:,:)
    case (2)
       B(:,:) = invMatCmplxSP2(A)
    case (3)
       B(:,:) = invMatCmplxSP3(A)
    case (4)
       B(:,:) = invMatCmplxSP4(A)
    case (5:)
       B(:,:) = invMatCmplxSPN(A)
    end select

  end function invMatComplexSP

  function invMatComplexDP(A) result(B)
    complex(DP), intent(in) :: A(:,:)
    complex(DP)             :: B(size(A,1),size(A,2))
    ! Local variables
    integer              :: n

    n = size(A,1)
    select case (n)
    case (1)
       if (A(1,1).eq.cmplx(0.0_DP,0.0_DP)) then
          write(*,*) 'Matrix is numerically singular!'
       end if
       B(:,:) = 1.0_DP / A(:,:)
    case (2)
       B(:,:) = invMatCmplxDP2(A)
    case (3)
       B(:,:) = invMatCmplxDP3(A)
    case (4)
       B(:,:) = invMatCmplxDP4(A)
    case (5:)
       B(:,:) = invMatCmplxDPN(A)
    end select

  end function invMatComplexDP

  function invMatDP2(A) result(B)
    ! Performs a direct calculation of the inverse of a 2×2 matrix.
    real(dp), intent(in) :: A(2,2)   ! Matrix
    real(dp)             :: B(2,2)   ! Inverse matrix
    real(dp)             :: det

    ! Calculate the inverse determinant of the matrix
    det = A(1,1)*A(2,2) - A(1,2)*A(2,1)

    if (det.eq.0.0_DP) then
       write(*,*) 'Matrix is numerically singular!'
    end if

    ! Calculate the inverse of the matrix
    B(1,1) =  A(2,2)
    B(2,1) = -A(2,1)
    B(1,2) = -A(1,2)
    B(2,2) =  A(1,1)
    B(:,:) =  B(:,:)/det
  end function invMatDP2

  function invMatDP3(A) result(B)
    ! Performs a direct calculation of the inverse of a 3×3 matrix.
    real(dp), intent(in) :: A(3,3)   ! Matrix
    real(dp)             :: B(3,3)   ! Inverse matrix
    real(dp)             :: det

    ! Calculate the inverse determinant of the matrix
    det = (A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
         - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
         + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

    if (det.eq.0.0_DP) then
       write(*,*) 'Matrix is numerically singular!'
    end if

    ! Calculate the inverse of the matrix
    B(1,1) =  (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -(A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) =  (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -(A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) =  (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -(A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) =  (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -(A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) =  (A(1,1)*A(2,2) - A(1,2)*A(2,1))
    B(:,:) =  B(:,:)/det
  end function invMatDP3

  function invMatDP4(A) result(B)
    ! Performs a direct calculation of the inverse of a 4×4 matrix.
    real(dp), intent(in) :: A(4,4)   ! Matrix
    real(dp)             :: B(4,4)   ! Inverse matrix
    real(dp)             :: det

    ! Calculate the inverse determinant of the matrix
    det = &
         1/(A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
         - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
         + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
         - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))

    if (det.eq.0.0_DP) then
       write(*,*) 'Matrix is numerically singular!'
    end if

    ! Calculate the inverse of the matrix
    B(1,1) = (A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
    B(2,1) = (A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
    B(3,1) = (A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(4,1) = (A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(1,2) = (A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
    B(2,2) = (A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
    B(3,2) = (A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(4,2) = (A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(1,3) = (A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
    B(2,3) = (A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
    B(3,3) = (A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
    B(4,3) = (A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
    B(1,4) = (A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
    B(2,4) = (A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
    B(3,4) = (A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
    B(4,4) = (A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
    B(:,:) = B(:,:)*det
  end function invMatDP4

  function invMatDPN(A) result(Ainv)

    real(dp), dimension(:,:), intent(in) :: A
    real(dp), dimension(size(A,1),size(A,2)) :: Ainv

    real(dp), dimension(size(A,1)) :: work  ! work array for LAPACK
    integer, dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: n, info

    ! External procedures defined in LAPACK
    external DGETRF
    external DGETRI

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call dgetrf(n, n, Ainv, n, ipiv, info)

    if (info /= 0) then
       stop 'Matrix is numerically singular!'
    end if

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call dgetri(n, Ainv, n, ipiv, work, n, info)

    if (info /= 0) then
       stop 'Matrix inversion failed!'
    end if
  end function invMatDPN

  function invMatSP2(A) result(B)
    ! Performs a direct calculation of the inverse of a 2×2 matrix.
    real(SP), intent(in) :: A(2,2)   ! Matrix
    real(SP)             :: B(2,2)   ! Inverse matrix
    real(SP)             :: det

    ! Calculate the inverse determinant of the matrix
    det = A(1,1)*A(2,2) - A(1,2)*A(2,1)

    if (det.eq.0.0) then
       write(*,*) 'Matrix is numerically singular!'
    end if

    ! Calculate the inverse of the matrix
    B(1,1) =  A(2,2)
    B(2,1) = -A(2,1)
    B(1,2) = -A(1,2)
    B(2,2) =  A(1,1)
    B(:,:) =  B(:,:)/det
  end function invMatSP2

  function invMatSP3(A) result(B)
    ! Performs a direct calculation of the inverse of a 3×3 matrix.
    real(SP), intent(in) :: A(3,3)   ! Matrix
    real(SP)             :: B(3,3)   ! Inverse matrix
    real(SP)             :: det

    ! Calculate the inverse determinant of the matrix
    det = (A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
         - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
         + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

    if (det.eq.0.0) then
       write(*,*) 'Matrix is numerically singular!'
    end if

    ! Calculate the inverse of the matrix
    B(1,1) =  (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -(A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) =  (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -(A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) =  (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -(A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) =  (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -(A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) =  (A(1,1)*A(2,2) - A(1,2)*A(2,1))
    B(:,:) =  B(:,:)/det
  end function invMatSP3

  function invMatSP4(A) result(B)
    ! Performs a direct calculation of the inverse of a 4×4 matrix.
    real(SP), intent(in) :: A(4,4)   ! Matrix
    real(SP)             :: B(4,4)   ! Inverse matrix
    real(SP)             :: det

    ! Calculate the inverse determinant of the matrix
    det = &
         1/(A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
         - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
         + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
         - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))

    if (det.eq.0.0) then
       write(*,*) 'Matrix is numerically singular!'
    end if

    ! Calculate the inverse of the matrix
    B(1,1) = (A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
    B(2,1) = (A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
    B(3,1) = (A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(4,1) = (A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(1,2) = (A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
    B(2,2) = (A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
    B(3,2) = (A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(4,2) = (A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(1,3) = (A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
    B(2,3) = (A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
    B(3,3) = (A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
    B(4,3) = (A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
    B(1,4) = (A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
    B(2,4) = (A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
    B(3,4) = (A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
    B(4,4) = (A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
    B(:,:) = B(:,:)*det
  end function invMatSP4

  function invMatSPN(A) result(Ainv)

    real(SP), dimension(:,:), intent(in) :: A
    real(SP), dimension(size(A,1),size(A,2)) :: Ainv

    real(SP), dimension(size(A,1)) :: work  ! work array for LAPACK
    integer, dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: n, info

    ! External procedures defined in LAPACK
    external DGETRF
    external DGETRI

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call dgetrf(n, n, Ainv, n, ipiv, info)

    if (info /= 0) then
       stop 'Matrix is numerically singular!'
    end if

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call dgetri(n, Ainv, n, ipiv, work, n, info)

    if (info /= 0) then
       stop 'Matrix inversion failed!'
    end if
  end function invMatSPN

  function invMatCmplxSP2(A) result(B)
    ! Performs a direct calculation of the inverse of a 2×2 matrix.
    complex(SP), intent(in) :: A(2,2)   ! Matrix
    complex(SP)             :: B(2,2)   ! Inverse matrix
    complex(SP)             :: det

    ! Calculate the inverse determinant of the matrix
    det = (A(1,1)*A(2,2) - A(1,2)*A(2,1))

    if (det.eq.cmplx(0.0,0.0)) then
       write(*,*) 'Matrix is numerically singular!'
    end if

    ! Calculate the inverse of the matrix
    B(1,1) =  A(2,2)
    B(2,1) = -A(2,1)
    B(1,2) = -A(1,2)
    B(2,2) =  A(1,1)
    B(:,:) =  B(:,:)/det
  end function invMatCmplxSP2

  function invMatCmplxSP3(A) result(B)
    ! Performs a direct calculation of the inverse of a 3×3 matrix.
    complex(SP), intent(in) :: A(3,3)   ! Matrix
    complex(SP)             :: B(3,3)   ! Inverse matrix
    complex(SP)             :: det

    ! Calculate the inverse determinant of the matrix
    det = (A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
         - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
         + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

    if (det.eq.cmplx(0.0,0.0)) then
       write(*,*) 'Matrix is numerically singular!'
    end if

    ! Calculate the inverse of the matrix
    B(1,1) =  (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -(A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) =  (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -(A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) =  (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -(A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) =  (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -(A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) =  (A(1,1)*A(2,2) - A(1,2)*A(2,1))
    B(:,:) =  B(:,:)/det
  end function invMatCmplxSP3

  function invMatCmplxSP4(A) result(B)
    ! Performs a direct calculation of the inverse of a 4×4 matrix.
    complex(SP), intent(in) :: A(4,4)   ! Matrix
    complex(SP)             :: B(4,4)   ! Inverse matrix
    complex(SP)             :: det

    ! Calculate the inverse determinant of the matrix
    det = &
         (A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
         - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
         + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
         - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))

    if (det.eq.cmplx(0.0,0.0)) then
       write(*,*) 'Matrix is numerically singular!'
    end if

    ! Calculate the inverse of the matrix
    B(1,1) = (A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
    B(2,1) = (A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
    B(3,1) = (A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(4,1) = (A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(1,2) = (A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
    B(2,2) = (A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
    B(3,2) = (A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(4,2) = (A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(1,3) = (A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
    B(2,3) = (A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
    B(3,3) = (A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
    B(4,3) = (A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
    B(1,4) = (A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
    B(2,4) = (A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
    B(3,4) = (A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
    B(4,4) = (A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
    B(:,:) = B(:,:)/det
  end function invMatCmplxSP4

  function invMatCmplxSPN(A) result(Ainv)
    complex(SP), dimension(:,:), intent(in) :: A
    complex(SP), dimension(size(A,1),size(A,2)) :: Ainv

    complex(SP), dimension(size(A,1)) :: work  ! work array for LAPACK
    integer, dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: n, info

    ! External procedures defined in LAPACK
    external DGETRF
    external DGETRI

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call zgetrf(n, n, Ainv, n, ipiv, info)

    if (info .ne. 0) then
       stop 'Matrix is numerically singular!'
    end if

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call zgetri(n, Ainv, n, ipiv, work, n, info)

    if (info .ne. 0) then
       stop 'Matrix inversion failed!'
    end if
  end function invMatCmplxSPN

  function invMatCmplxDP2(A) result(B)
    ! Performs a direct calculation of the inverse of a 2×2 matrix.
    complex(DP), intent(in) :: A(2,2)   ! Matrix
    complex(DP)             :: B(2,2)   ! Inverse matrix
    complex(DP)             :: det

    ! Calculate the inverse determinant of the matrix
    det = (A(1,1)*A(2,2) - A(1,2)*A(2,1))

    if (det.eq.cmplx(0.0_DP,0.0_DP)) then
       write(*,*) 'Matrix is numerically singular!'
    end if

    ! Calculate the inverse of the matrix
    B(1,1) =  A(2,2)
    B(2,1) = -A(2,1)
    B(1,2) = -A(1,2)
    B(2,2) =  A(1,1)
    B(:,:) =  B(:,:)/det
  end function invMatCmplxDP2

  function invMatCmplxDP3(A) result(B)
    ! Performs a direct calculation of the inverse of a 3×3 matrix.
    complex(DP), intent(in) :: A(3,3)   ! Matrix
    complex(DP)             :: B(3,3)   ! Inverse matrix
    complex(DP)             :: det

    ! Calculate the inverse determinant of the matrix
    det = (A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
         - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
         + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

    if (det.eq.cmplx(0.0_DP,0.0_DP)) then
       write(*,*) 'Matrix is numerically singular!'
    end if

    ! Calculate the inverse of the matrix
    B(1,1) =  (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -(A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) =  (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -(A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) =  (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -(A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) =  (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -(A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) =  (A(1,1)*A(2,2) - A(1,2)*A(2,1))
    B(:,:) =  B(:,:)/det
  end function invMatCmplxDP3

  function invMatCmplxDP4(A) result(B)
    ! Performs a direct calculation of the inverse of a 4×4 matrix.
    complex(DP), intent(in) :: A(4,4)   ! Matrix
    complex(DP)             :: B(4,4)   ! Inverse matrix
    complex(DP)             :: det

    ! Calculate the inverse determinant of the matrix
    det = &
         (A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
         - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
         + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
         - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))

    if (det.eq.cmplx(0.0_DP,0.0_DP)) then
       write(*,*) 'Matrix is numerically singular!'
    end if

    ! Calculate the inverse of the matrix
    B(1,1) = (A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
    B(2,1) = (A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
    B(3,1) = (A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(4,1) = (A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(1,2) = (A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
    B(2,2) = (A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
    B(3,2) = (A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(4,2) = (A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(1,3) = (A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
    B(2,3) = (A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
    B(3,3) = (A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
    B(4,3) = (A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
    B(1,4) = (A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
    B(2,4) = (A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
    B(3,4) = (A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
    B(4,4) = (A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
    B(:,:) = B(:,:)/det
  end function invMatCmplxDP4

  function invMatCmplxDPN(A) result(Ainv)
    complex(DP), dimension(:,:), intent(in) :: A
    complex(DP), dimension(size(A,1),size(A,2)) :: Ainv

    complex(DP), dimension(size(A,1)) :: work  ! work array for LAPACK
    integer, dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: n, info

    ! External procedures defined in LAPACK
    external DGETRF
    external DGETRI

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call zgetrf(n, n, Ainv, n, ipiv, info)

    if (info .ne. 0) then
       stop 'Matrix is numerically singular!'
    end if

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call zgetri(n, Ainv, n, ipiv, work, n, info)

    if (info .ne. 0) then
       stop 'Matrix inversion failed!'
    end if
  end function invMatCmplxDPN


  !---------------------------------------------------------------------!
  !                                                                     !
  !                      Useful Mathamatical functions                  !
  !                                                                     !
  !---------------------------------------------------------------------!

  function factorial(n)
    implicit none

    integer,intent(in) :: n
    integer            :: factorial

    integer            :: i

    factorial = 1

    do i = 1,n
       factorial = factorial*i
    end do

  end function factorial

  !---------------------------------------------------------------------!
  !                                                                     !
  !                               MeshGrid                              !
  !                                                                     !
  !---------------------------------------------------------------------!
  ! creates a unique lattice of points given 2 vectors x,y of length    !
  ! N and M. Can be used for 3d plots. See Matlab meshgrid              !
  ! documentation for more details                                      !
  !---------------------------------------------------------------------!

  subroutine meshgrid(x,y,XX,YY)
    use kinds
    implicit none

    real(DP),dimension(:),intent(in)    :: x
    real(DP),dimension(:),intent(in)    :: y
    real(DP),dimension(:,:),intent(out) :: XX,YY

    integer :: i,j,rows,cols

    rows = size(x)
    cols = size(y)

    do j = 1,cols
       XX(j,1:rows) = x(1:rows)
    end do
    do j = 1,rows
       YY(1:cols,j) = y(1:cols)
    end do

  end subroutine meshgrid

  !---------------------------------------------------------------------!
  !                                                                     !
  !                           Cubic Spline Fit                          !
  !                                                                     !
  !---------------------------------------------------------------------!
  ! Fits a cubic spline to input data. First version fits and outputs   !
  ! values at specified interpolated points. Next version returns       !
  ! the coefficients of the cubic spline fit which can then be passed   !
  ! to the next function to evaluate the spline at a certain x.         !
  !                                                                     !
  ! Uses Lapack for necessary linear algebra                            !
  !---------------------------------------------------------------------!

  function splinefit(x,y) result(c)
    use kinds
    use LAPACK95
    implicit none

    real(DP),dimension(:),intent(in)    :: x,y
    real(DP),dimension(size(x))         :: c

    integer                             :: N
    real(DP),dimension(:,:),allocatable :: A,B
    integer ,dimension(:)  ,allocatable :: ipiv
    integer                             :: i,j

    N = size(x)
    allocate(A(N,N),B(N,1),ipiv(N))

    B(:,1) = y(:)
    A(:,1) = 1.0_DP
    A(:,2) = x(:)
    do i = 1,N
       do j = 3,N
          A(i,j) = abs(x(i)-x(j-1))
       end do
    end do
    call getrf(A,ipiv)
    call getrs(A,ipiv,B)
    c(:) = B(:,1)

    deallocate(A,B,ipiv)

  end function splinefit

  function splineVal(c,xj,x)
    use kinds
    implicit none

    real(DP),dimension(:),intent(in) :: c,xj
    real(DP),intent(in)              :: x
    real(DP)                         :: splineval

    splineval = c(1)+c(2)*x+sum(c(3:size(c))*abs(x-xj(2:size(c)-1)))

  end function splineVal

  !---------------------------------------------------------------------!
  !                                                                     !
  !                               PolyFit                               !
  !                                                                     !
  !---------------------------------------------------------------------!
  ! N dimensional polynomial fitting algorithim which outputs the       !
  ! coefficitents of the fit. Can then be passed into polycal to        !
  ! evaluate the fit at a certain point.                                !
  !                                                                     !
  ! Uses Lapack for necessary linear algebra                            !
  !---------------------------------------------------------------------!

  function polyVal(c,x) result(val)
    use kinds
    implicit none

    real(DP),dimension(:),intent(in) :: c
    real(DP)             ,intent(in) :: x

    real(DP)                           :: val
    integer                            :: i,powers(size(c))

    powers = linspace(0,size(c)-1,size(c))
    val = sum(c*x**powers)

  end function polyVal


  function polyfit(vx, vy, d)
    use kinds
    use Lapack95
    implicit none
    integer, intent(in)                   :: d
    real(dp), dimension(d+1)              :: polyfit
    real(dp), dimension(:), intent(in)    :: vx, vy

    real(dp), dimension(:,:), allocatable :: X
    real(dp), dimension(:,:), allocatable :: XT
    real(dp), dimension(:,:), allocatable :: XTX

    integer :: i, j

    integer     :: n, lda, lwork
    integer :: info
    integer, dimension(:), allocatable :: ipiv
    real(dp), dimension(:), allocatable :: work

    n = d+1
    lda = n
    lwork = n

    allocate(ipiv(n))
    allocate(work(lwork))
    allocate(XT(n, size(vx)))
    allocate(X(size(vx), n))
    allocate(XTX(n, n))

    ! prepare the matrix
    do i = 0, d
       do j = 1, size(vx)
          X(j, i+1) = vx(j)**i
       end do
    end do

    XT  = transpose(X)
    XTX = matmul(XT, X)

    ! calls to LAPACK subs DGETRF and DGETRI
    call DGETRF(n, n, XTX, lda, ipiv, info)
    if ( info /= 0 ) then
       print *, "problem"
       return
    end if
    call DGETRI(n, XTX, lda, ipiv, work, lwork, info)
    if ( info /= 0 ) then
       print *, "problem"
       return
    end if

    polyfit = matmul( matmul(XTX, XT), vy)

    deallocate(ipiv)
    deallocate(work)
    deallocate(X)
    deallocate(XT)
    deallocate(XTX)

  end function polyfit

  !---------------------------------------------------------------------!
  !                                                                     !
  !                         Polynomial Integral                         !
  !                                                                     !
  !---------------------------------------------------------------------!

  function polyInt(c,a,b) result(val)
    use kinds
    implicit none

    real(DP), dimension(:),intent(in) :: c
    real(DP)              ,intent(in) :: a,b
    real(DP)                          :: val

    real(DP) :: aval,bval
    integer  :: ii

    bval = 0
    aval = 0

    do ii = 1,size(c)
       bval = bval+c(ii)*(b**ii)/ii
       aval = aval+c(ii)*(a**ii)/ii
    end do

    val = bval - aval

  end function polyInt

  !---------------------------------------------------------------------!
  !                                                                     !
  !                            Eulers method                            !
  !                                                                     !
  !---------------------------------------------------------------------!
  ! Eulers method for solving 1 or N coupled DE's                       !
  !---------------------------------------------------------------------!

  function Euler(f,h,t0,y0)
    use kinds
    implicit none

    real(DP),intent(in) :: h,t0,y0
    interface
       function f(t,y)
         use kinds
         implicit none
         real(DP)             :: f
         real(DP), intent(in) :: t,y
       end function f
    end interface
    real(DP) :: Euler

    Euler = h*f(t0,y0)+y0

  end function Euler

  function eulerND(f,h,t0,y0)
    use kinds
    implicit none

    real(DP),intent(in)              :: t0,h
    real(DP),dimension(:),intent(in) :: y0
    interface
       function f(t0,y0,nEq)
         use kinds
         implicit none
         integer,intent(in)    :: nEq
         real(DP),intent(in)   :: t0,y0(nEq)
         real(DP),dimension(nEq) :: f
       end function f
    end interface
    real(DP),dimension(size(y0)) :: eulerND
    integer :: nEq

    nEq = size(y0)
    eulerND = h*f(t0,y0,nEq)+y0

  end function eulerND

  !---------------------------------------------------------------------!
  !                                                                     !
  !                       Runge Kutta 4th Order                         !
  !                                                                     !
  !---------------------------------------------------------------------!
  ! Runge Kutta calculation for input function(s). For a single step    !
  ! or a range of specified values. Algorithims for a single DE or      !
  ! N coupled DE's.                                                     !
  !---------------------------------------------------------------------!

  function rk4(f,h,t0,y0)
    use kinds
    implicit none

    real(DP),intent(in) :: h,t0,y0
    interface
       function f(t,y)
         use kinds
         implicit none
         real(DP)             :: f
         real(DP), intent(in) :: t,y
       end function f
    end interface
    real(DP) :: rk4

    real(DP) :: k1,k2,k3,k4

    k1 = f(t0,y0)
    k2 = f(t0+0.5_dp*h,y0+0.5_dp*h*k1)
    k3 = f(t0+0.5_dp*h,y0+0.5_dp*h*k2)
    k4 = f(t0+h,y0+k3*h)
    rk4 = y0 + (1/6.0_dp)*(k1+2*k2+2*k3+k4)*h

  end function rk4

  function rk4N(f,h,t0,y0)
    use kinds
    implicit none

    real(DP),intent(in)              :: t0,h
    real(DP),dimension(:),intent(in) :: y0
    interface
       function f(t0,y0,nEq)
         use kinds
         implicit none
         integer,intent(in)    :: nEq
         real(DP),intent(in)   :: t0,y0(nEq)
         real(DP),dimension(nEq) :: f
       end function f
    end interface
    real(DP),dimension(size(y0)) :: rk4N

    real(DP),dimension(size(y0)) :: k1,k2,k3,k4
    integer :: nEq

    nEq = size(y0)
    k1 = f(t0,y0,nEq)
    k2 = f(t0 + 0.5_DP*h,y0 + 0.5_DP*h*k1,nEq)
    k3 = f(t0 + 0.5_DP*h,y0 + 0.5_DP*h*k2,nEq)
    k4 = f(t0 + h,y0 + h*k3,nEq)

    rk4N = y0 + (1/6.0_DP)*(k1 + 2*k2 + 2*k3 + k4)*h

  end function rk4N

  !---------------------------------------------------------------------!
  !                                                                     !
  !                               Guess Zero                            !
  !                                                                     !
  !---------------------------------------------------------------------!
  ! Returns the points where a function changes sign. Can be an input   !
  ! function of an input set of values.                                 !
  !---------------------------------------------------------------------!

  function guessZero(fvals)
    use Kinds
    implicit none

    real(DP),dimension(:),intent(in) :: fvals
    integer                          :: GuessZero

    integer  :: j,sign,newsign

    sign = int(fvals(1)/abs(fvals(1)))
    do j = 2,size(fvals)
       newsign = int(fvals(j)/abs(fvals(j)))
       if ( newsign .ne. sign) exit
    end do
    if ( j .eq. size(fvals) ) then
       write(*,'(a)') "ERROR did not see any sign change/ zero"
       write(*,'(a)') "maybe try again with more points"
    else
       GuessZero = j
    end if

  end function GuessZero

  function guessZeroF(f,a,b)
    use kinds
    implicit none

    real(DP),intent(in) :: a,b
    real(DP)            :: GuessZeroF
    interface
       function f(x)
         use kinds
         implicit none
         real(DP),intent(in) :: x
         real(DP)            :: f
       end function f
    end interface

    real(DP) :: h,fval
    integer :: i,sign,newsign,N

    h = a*1e-1
    N = int((b-a)/h)
    fval = f(a)
    sign = int(fval/abs(fval))

    do i = 1,N
       fval = f(a+i*h)
       newsign = int(fval/abs(fval))
       if ( newsign .ne. sign ) exit
    end do
    if ( i .eq. N ) then
       write(*,*) "Error Max itterations reached"
    else
       GuessZeroF = a+i*h
    end if

  end function GuessZeroF

  !---------------------------------------------------------------------!
  !                                                                     !
  !                         Newtons Method                              !
  !                                                                     !
  !---------------------------------------------------------------------!
  ! Newtons method for an input function of 1 variable                  !
  !                                                                     !
  !                 CODE WRITTEN BY CURTIS ABELL                        !
  !---------------------------------------------------------------------!

  function newton1D(fn,guess)
    use Kinds
    implicit none
    real(DP) :: newton1D
    interface
       function fn(x)
         use kinds
         implicit none
         real(DP)             :: fn
         real(DP), intent(in) :: x
       end function fn
    end interface
    real(DP) :: guess
    ! Local Variables
    real(DP) :: newt_tol = 1e-6
    real(DP) :: h ! = 1e-6_DP
    real(DP) :: fx, fxfh,fxbh, dfdx
    integer :: attempt
    integer :: attempt_limit = 20
    logical :: verbose = .false.

    if (verbose) write(*,*) 'Newton-Raphson Method'
    h = guess * 1.0e-6_DP
    newton1D = guess
    attempt = 1

    do
       fx = fn(newton1D)
       if (abs(fx) <= newt_tol) then
          if (verbose) then
             write(*,'(a,i3,a)') 'Success after', attempt, ' attempts.'
             write(*,*) 'Zero  = ', newton1D
          end if
          exit
       else if (attempt >= attempt_limit) then
          write(*,*) 'Failed, change initial guess or increase attempt limit.'
          exit
       end if

       fxfh = fn(newton1D + h/2)
       fxbh = fn(newton1D - h/2)
       dfdx = (fxfh-fxbh) / h
       newton1D = newton1D - fx / dfdx

       attempt = attempt + 1
    end do

  end function newton1D

  !---------------------------------------------------------------------!
  !                                                                     !
  !                               Linspace                              !
  !                                                                     !
  !---------------------------------------------------------------------!
  ! Creates a linear space of points. See matlab documentation for      !
  ! further details.                                                    !
  !---------------------------------------------------------------------!

  function linspace(start,finish,N)
    use kinds
    implicit none

    integer                :: N
    real(DP), intent(in)   :: start, finish
    real(DP), dimension(N) :: linspace

    real(DP)               :: int
    integer                :: i

    linspace(1) = start
    linspace(N) = finish
    int  = (real(finish)-start)/(N-1)

    do i=2,N-1
       linspace(i)=linspace(i-1)+int
    end do

  end function linspace

  function linspaceReal(start,finish,N)
    use kinds
    implicit none

    integer                :: N
    real, intent(in)       :: start, finish
    real, dimension(N)     :: linspaceReal

    real(DP)               :: int
    integer                :: i

    linspaceReal(1) = start
    linspaceReal(N) = finish
    int  = (real(finish)-start)/(N-1)

    do i=2,N-1
       linspaceReal(i)=linspaceReal(i-1)+int
    end do

  end function linspaceReal

  function linspaceInt(start,finish,N)
    use kinds
    implicit none 

    integer                :: N
    integer, intent(in)    :: start, finish
    integer, dimension(N)  :: linspaceInt

    real(DP)               :: int
    integer                :: i

    linspaceInt(1) = start
    linspaceInt(N) = finish
    int  = (real(finish)-start)/(N-1)

    do i=2,N-1
       linspaceInt(i)=linspaceInt(i-1)+int
    end do

  end function linspaceInt

  !---------------------------------------------------------------------!
  !                                                                     !
  !                         Iteration loop timer                        !
  !                                                                     !
  !---------------------------------------------------------------------!

  subroutine LoopTimer(i,L,t0,tF)
    use kinds
    implicit none

    real(DP), intent(in)  :: t0
    real(DP), intent(out) :: tF
    integer , intent(in)  :: i,L

    character(len=50)  :: bar,barNew
    real(DP):: delta
    integer :: ii,timeTotal

    call cpu_time(tF)

    bar = '='
    barNew = bar
    delta = tF-t0

    timeTotal = int(delta*L-delta*i)

    do ii = 1,((real(i)/L)*100)/2-1
       barNew = trim(barNew)//trim(bar)
    end do

    barNew = trim(barNew)//'>'

    write(*,'(a13,i3,a2,a50,a2,f5.1,a11,f5.1,a3)',advance='no') &
         & '| Iteration: ',i,' |',barNew,'| ',(real(i)/L)*100,'% | It Len ',delta,&
         & 's |'

    if ( timeTotal/60 .le. 9 .and. mod(timeTotal,60) .le. 9)   write(*,'(i2,a2,i1,a2,i1)') (timeTotal/60)/60,':0',(timeTotal/60),':0',mod(timeTotal,60)
    if ( timeTotal/60 .le. 9 .and. mod(timeTotal,60) .ge. 10)  write(*,'(i2,a2,i1,a1,i2)') (timeTotal/60)/60,':0',(timeTotal/60),':',mod(timeTotal,60)
    if ( timeTotal/60 .ge. 10 .and. mod(timeTotal,60) .le. 9)  write(*,'(i2,a1,i2,a2,i1)') (timeTotal/60)/60,':',(timeTotal/60),':0',mod(timeTotal,60)
    if ( timeTotal/60 .ge. 10 .and. mod(timeTotal,60) .ge. 10) write(*,'(i2,a1,i2,a1,i2)') (timeTotal/60)/60,':',(timeTotal/60),':',mod(timeTotal,60)

  end subroutine LoopTimer

  !---------------------------------------------------------------------!
  !                                                                     !
  !                          Finite difference                          !
  !                                                                     !
  !---------------------------------------------------------------------!
  ! Numerically calculates the derivative with finite difference        !
  !---------------------------------------------------------------------!

  function deriv(f,x0)
    use Kinds
    implicit none
    real(DP)                 :: deriv
    interface
       function f(x)
         use kinds
         implicit none
         real(DP)            :: f
         real(DP),intent(in) :: x
       end function f
    end interface
    real(DP),intent(in)      :: x0
    real(DP)                 :: h,fxfh,fxbh

    h     = x0*1e-6_DP
    fxfh  = f(x0+h/2)
    fxbh  = f(x0-h/2)
    deriv = (fxfh-fxbh)/h

  end function deriv


  !---------------------------------------------------------------------!
  !                                                                     !
  !                       Trapezodial Integration                       !
  !                                                                     !
  !---------------------------------------------------------------------!

  function trapZ(f,a,b,relErr) result(val)
    use kinds
    implicit none

    real(DP), intent(in) :: a,b,relErr
    interface
       function f(x)
         use kinds
         implicit none
         real(DP),intent(in) :: x
         real(DP)            :: f
       end function f
    end interface
    real(DP) :: val

    real(DP) :: width,f1,f2
    integer  :: N,i

    val = 0
    width = relErr
    N = int((b-a)/width)

    do i = 2,N+1
       f1 = f(a+width*(i-1))
       f2 = f(a+width*(i-2))
       val = val+(f1+f2)*width/2
    end do

  end function trapZ


  function trapZZ(x,y) result(val)
    use kinds
    implicit none

    real(DP), dimension(:), intent(in) :: x,y
    real(DP) :: val

    real(DP) :: f1,f2,width
    integer  :: N,i

    val = 0
    N = size(x)
    width = abs(x(1)-x(2))

    do i = 2,N
       f1 = y(i)
       f2 = y(i-1)
       val = val+(f1+f2)*width/2
    end do

  end function trapZZ

  !---------------------------------------------------------------------!
  !                                                                     !
  !                     General Minimization Wrapper                    !
  !                                                                     !
  !---------------------------------------------------------------------!

  subroutine Minimize(TestFun, n, x, e, maxStepErrorScaleFactor, minimum)
    use kinds
    implicit none
    interface
       subroutine TestFun(n, x, f)
         USE kinds 
         implicit none
         integer               , intent(in)  :: n
         real(DP), dimension(n), intent(in)  :: x
         real(DP)              , intent(out) :: f
       end subroutine TestFun
    end interface
    integer,                intent(in)    :: n
    real(DP), dimension(n), intent(inout) :: x, e
    real(DP),               intent(in)    :: maxStepErrorScaleFactor
    real(DP),               intent(out)   :: minimum

    ! Local variables to interface with minf
    ! Written by Derek Leinweber
    !
    ! These variables control how minf works
    !
    ! icon=2 helps to avoid local minima.  Highly recommended.
    ! See minf for details on the other parameters
    !
    integer, parameter           :: icon=2
    integer                      :: iprint=2, maxit=100000, iterCount
    integer                      :: nwork
    real(DP), dimension(n*(n+3)) :: w

    nwork = n*(n+3)

    call minf( TestFun, x, e, n, minimum, maxStepErrorScaleFactor, &
         &     iprint, icon, maxit, iterCount, w, nwork )

  end Subroutine Minimize

  !---------------------------------------------------------------------!
  !                                                                     !
  !                      QuadPack Integral wrapper                      !
  !                                                                     !
  !---------------------------------------------------------------------!
  ! Wrappers for the Quadpack integral routine shown above.             !
  !---------------------------------------------------------------------!

  function integral(f, a, b, absErr, relErr)
    implicit none
    real(DP)                  :: integral
    interface                                   ! Interfaces: Sections 5.11, 5.18
       function f(x)
         use kinds
         implicit none
         real(DP)             :: f
         real(DP), intent(in) :: x
       end function f
    end interface
    real(DP), intent(in)      :: a, b, absErr, relErr

    !  Local variables
    !
    real(DP)                  :: integResult, bound=0.0_DP
    integer                   :: inf
    integer  :: neval, ifail
    real(DP) :: errEstimate

    !  Determine if the limits include infinity and call qagi if nessary
    !
    if (a == -Infty) then
       if (b == Infty) then
          inf=2
          call qagi(f, bound, inf, absErr, relErr, integResult, errEstimate, neval, ifail)
          if ( ifail /= 0 ) then
             write(*,*) 'Warning from qagi: the error code is ', ifail
          end if
       else
          inf = -1
          bound = b
          call qagi(f, bound, inf, absErr, relErr, integResult, errEstimate, neval, ifail)
          if ( ifail /= 0 ) then
             write(*,*) 'Warning from qagi: the error code is ', ifail
          end if
       end if
    else
       if (b == Infty) then
          inf = 1
          bound = a
          call qagi(f, bound, inf, absErr, relErr, integResult, errEstimate, neval, ifail)
          if ( ifail /= 0 ) then
             write(*,*) 'Warning from qagi: the error code is ', ifail
          end if
       else
          call qags(f, a, b, absErr, relErr, integResult, errEstimate, neval, ifail)
          if ( ifail /= 0 ) then
             write(*,*) 'Warning from qags: the error code is ', ifail
          end if
       end if
    end if
    integral = integResult

  end function integral


  function integralToInfty(f, bound, absErr, relErr)
    implicit none
    real(DP)                  :: integralToInfty
    interface
       function f(x)
         use kinds
         implicit none
         real(DP)             :: f
         real(DP), intent(in) :: x
       end function f
    end interface
    real(DP), intent(in)      :: bound, absErr, relErr

    !  Local variables
    !
    real(DP)                  :: integResult
    integer, parameter        :: inf=1
    integer  :: neval, ifail
    real(DP) :: errEstimate

    call qagi(f, bound, inf, absErr, relErr, integResult, errEstimate, neval, ifail)
    if ( ifail /= 0 ) then
       write(*,*) 'Warning from qagi: the error code is ', ifail
    end if
    integralToInfty = integResult

  end function integralToInfty


  function integralOf(f, absErr, relErr)
    implicit none
    real(DP)                  :: integralOf
    interface
       function f(x)
         implicit none
         integer,  parameter  :: DP = kind(1.0d0)
         real(DP)             :: f
         real(DP), intent(in) :: x
       end function f
    end interface
    real(DP), intent(in)      :: absErr, relErr

    !  Local variables
    !
    real(DP)                  :: integResult, bound=0.0d0
    integer, parameter        :: inf=2
    integer  :: neval, ifail
    real(DP) :: errEstimate

    call qagi(f, bound, inf, absErr, relErr, integResult, errEstimate, neval, ifail)
    if ( ifail /= 0 ) then
       write(*,*) 'Warning from qagi: the error code is ', ifail
    end if
    integralOf = integResult

  end function integralOf


  function integralBreakPts(f, a, b, absErr, relErr, nBreakPts, BreakPts)
    implicit none
    real(DP)                  :: integralBreakPts
    interface
       function f(x)
         use kinds
         implicit none
         real(DP)             :: f
         real(DP), intent(in) :: x
       end function f
    end interface
    real(DP), intent(in)      :: a, b, absErr, relErr
    integer,  intent(in)      :: nBreakPts
    real(DP), intent(in), dimension(nBreakPts) :: BreakPts

    !  Local variables
    !
    real(DP)                         :: integResult
    real(DP), dimension(nBreakPts+2) :: BreakPtsP2        ! Automatic array.  Similar to allocatable arrays.
    integer  :: neval, ifail
    real(DP) :: errEstimate

    BreakPtsP2(1:nBreakPts) = BreakPts(1:nBreakPts)       ! Array section limits are required here.

    call qagp(f, a, b, nBreakPts+2, BreakPtsP2, absErr, relErr, integResult, errEstimate, neval, ifail)
    if ( ifail /= 0 ) then
       write(*,*) 'Warning from qagp: the error code is ', ifail
    end if
    integralBreakPts = integResult

  end function integralBreakPts

  !---------------------------------------------------------------------!
  !                                                                     !
  !                      Principle Value Integrator                     !
  !                                                                     !
  !---------------------------------------------------------------------!
  ! Evaluates the Cauchy-Principle value integral for f(x)/(x-c)        !
  !---------------------------------------------------------------------!

  function integralPV(f, c, a, b, absErr, relErr)
    use kinds
    implicit none
    real(DP)                  :: integralPV
    interface
       function f(x)
         use kinds
         implicit none
         real(DP)             :: f
         real(DP), intent(in) :: x
       end function f
    end interface
    real(DP), intent(in)      :: c, a, b, absErr, relErr

    ! ---------------------------Local Variables--------------------------
    real(DP)                  :: res, errRequest
    integer                   :: inf ! -1 for -infty, +1 for infty
    ! Used to split the integral if a or b are at plus/minux infty
    !    into semi-infinite integrals and a PV integral. This is done
    !    by moving the upper and/or lower bounds of the integral
    !    away from the pole and splitting the integral into 2 (3).
    real(DP)                  :: intBoundaryHigh, intBoundaryLow
    integer  :: neval, ifail
    real(DP) :: errEstimate


    integralPV = 0.0_DP
    intBoundaryLow = a
    intBoundaryHigh = b

    ! Integral from -Infty to lower bound
    if (a.eq.-Infty) then
       inf = -1
       if (abs(c).lt.2.0_DP) then
          intBoundaryLow = -4.0_DP
       else
          intBoundaryLow = c / 2.0_DP
       end if
       call qagi(f, intBoundaryLow, inf, absErr, relErr, res &
            & , errEstimate, neval, ifail)
       integralPV = integralPV + res
       if ( ifail.ne.0 ) then
          write(*,*) ' Warning from qagi from -Infty: the error code is ', ifail
       end if
    end if

    ! Integral from upper bound to Infty
    if (b.eq.Infty) then
       inf = 1
       if (abs(c).lt.2.0_DP) then
          intBoundaryHigh = 4.0_DP
       else
          intBoundaryHigh = c * 2.0_DP
       end if
       call qagi(fNonSingular, intBoundaryHigh, inf, absErr  &
            & , relErr, res, errEstimate, neval, ifail)
       integralPV = integralPV + res
       if ( ifail.ne.0 ) then
          write(*,*) ' Warning from qagi to Infty: the error code is ', ifail
       end if
    end if

    ! Principle-value integral
    call qawc(f, intBoundaryLow, intBoundaryHigh, c, absErr, relErr &
         & , res, errEstimate, neval, ifail)
    if ( ifail.ne.0 ) then
       write(*,*) ' Warning from qawc: the error code is ', ifail
    end if
    integralPV = integralPV + res

  contains
    ! Function to pass to the semi-infinite integrals
    !    - now includes the (x-c) on the denominator
    function fNonSingular(x)
      use kinds
      implicit none
      real(DP) :: fNonSingular
      real(DP), intent(in) :: x
      fNonSingular = f(x)/(x - c)
    end function fNonSingular

  end function integralPV

  !---------------------------------------------------------------------!
  !                                                                     !
  !                       Multi-Dimensional Pyplots                     !
  !                                                                     !
  !---------------------------------------------------------------------!

  subroutine pyplotN(x,title,xaxis,yaxis,legend)
    use Kinds
    implicit none
    character(len=*),intent(in),optional              :: xaxis,yaxis,title
    character(len=*),dimension(:),intent(in),optional :: legend
    real(DP),dimension(:,:),intent(in)                :: x

    character(len=10),dimension(size(x(1,:))/2)       :: ld
    character(len=14)                                 :: fmt,xlabel,ylabel,name
    integer                                           :: ii,j,N,L

    N = size(x,dim=2)
    L = size(ld)

    open(100,file="titles.dat",action="write", &
         & status="replace",form="formatted")
    open(101,file="output.dat",action="write", &
         & status="replace",form="formatted")

    name   = ""
    xlabel = ""
    ylabel = ""
    ld     = ""

    if(present(title)) name   = title
    if(present(xaxis)) xlabel = xaxis
    if(present(yaxis)) ylabel = yaxis
    if(present(legend))ld(1:L)= legend(1:L)

    write(100,'(a20)') name
    write(100,'(a15)') xlabel
    write(100,'(a15)') ylabel
    write(fmt,'(a1,i1,a7)') '(', N, 'es20.9)'

    do ii=1,L
       write(100,'(a12)') ld(ii)
    end do

    do ii=1,size(x,dim=1)
       write(101,fmt) x(ii,:)
    end do

    close(100)
    close(101)

  end subroutine pyplotN

  !---------------------------------------------------------------------!
  !                                                                     !
  !                          Standard Pyplots                           !
  !                                                                     !
  !---------------------------------------------------------------------!

  subroutine pyplotXY(x,y,title,xaxis,yaxis)
    use Kinds
    implicit none
    character(len=*),intent(in),optional :: xaxis,yaxis,title
    real(DP),dimension(:),intent(in)     :: x,y

    character(len=14)                    :: xlabel,ylabel,name,ld
    integer                              :: ii

    open(100,file="titles.dat",action="write", &
         & status="replace",form="formatted")
    open(101,file="output.dat",action="write", &
         & status="replace",form="formatted")

    name   = ""
    xlabel = ""
    ylabel = ""
    ld     = "empty"

    if(present(title)) name   = title
    if(present(xaxis)) xlabel = xaxis
    if(present(yaxis)) ylabel = yaxis

    write(100,'(a20)') name
    write(100,'(a15)') xlabel
    write(100,'(a15)') ylabel
    write(100,'(a7)') ld

    do ii=1,size(x)
       write(101,'(2es20.9)') x(ii), y(ii)
    end do

    close(100)
    close(101)

  end subroutine pyplotXY

  subroutine pyplotXYZW(x,y,z,w,title,xaxis,yaxis,legend)
    use Kinds
    implicit none
    character(len=*),dimension(2),intent(in),optional :: legend
    character(len=*),intent(in),optional :: xaxis,yaxis,title
    real(DP),dimension(:),intent(in)     :: x,y,z,w

    character(len=14)                    :: xlabel,ylabel,name
    character(len=14),dimension(2)       :: ld
    integer                              :: ii

    open(100,file="titles.dat",action="write", &
         & status="replace",form="formatted")
    open(101,file="output.dat",action="write", &
         & status="replace",form="formatted")

    name   = ""
    xlabel = ""
    ylabel = ""
    ld     = "empty"

    if(present(title)) name   = title
    if(present(xaxis)) xlabel = xaxis
    if(present(yaxis)) ylabel = yaxis
    if(present(yaxis)) ld     = legend

    write(100,'(a20)') name
    write(100,'(a15)') xlabel
    write(100,'(a15)') ylabel
    write(100,'(a7)') ld
    write(100,'(a7)') ld    

    do ii=1,size(x)
       write(101,'(4es20.9)') x(ii), y(ii), z(ii), w(ii)
    end do

    close(100)
    close(101)

  end subroutine pyplotXYZW

end module numFort

!==============================================================================
!##############################################################################
!==============================================================================
