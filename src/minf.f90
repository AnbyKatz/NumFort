MODULE MinFun
  implicit none

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                           !!
!!     SUBROUTINE minf                                                       !!
!!                                                                           !!
!!     Purpose:                                                              !!
!!     minf finds the minimum value of a function by iterative               !!
!!     variation of the function parameters.                                 !!
!!     The parameters need not be independent.                               !!
!!     The method is an updated version of a routine by m.powell             !!
!!     published in the computer journal(1964).                              !!
!!                                                                           !!
!!     Description Of Parameters:                                            !!
!!     n        number of parameters in optimization                         !!
!!     x        array of parameters.                                         !!
!!              on entry to minf x(i) must be set to an approximation of     !!
!!              i-th parameter.                                              !!
!!     e        array containing the absolute accuracies to which the        !!
!!              parameters are required.                                     !!
!!              the magnitudes of the e(i) are assumed to be roughly the     !!
!!              same as the x(i).                                            !!
!!     f        function value returned by funct.                            !!
!!     w        auxilliary work array of dimension nwork.                    !!
!!              nwork must be specified to be at least n*(n+3) in calling    !!
!!              program.                                                     !!
!!     icon     controls the ultimate convergence criteria                   !!
!!              =1      convergence will be assumed when an iteration        !!
!!                      changes each variable by less than 10 percent        !!
!!              =2      a point such as for icon=1 is found and then         !!
!!                      displaced by ten times the required accuracy         !!
!!                      in each parameter.  minimization is continued        !!
!!                      from the new point until a 10 percent accuracy       !!
!!                      is again achieved.                                   !!
!!     iprint   controls printing                                            !!
!!              =0      there will be no printing.                           !!
!!              =1      the variables and function will be printed after     !!
!!                      every search along a line (approximately every       !!
!!                      other function value).                               !!
!!              =2      the variables and function value will be printed     !!
!!                      after every iteration (n+1 searches along a line)    !!
!!     maxit    maximum number of iterations.                                !!
!!     iterc    the number of iterations actually completed.                 !!
!!     escale   limits the maximum change in the parameters.                 !!
!!              x(i) will not be changed by more than escale*e(i) at a       !!
!!              single step.                                                 !!
!!     functions and subroutines required:                                   !!
!!              funct(n,x,f)                                                 !!
!!              external subroutine to caculate the function value for       !!
!!              a given set of the parameters x(i).                          !!
!!                                                                           !!
!!                                                                           !!
!!     NOTE: Code updated to remove obsolete Fortran features by             !!
!!           Stewart V. Wright  8 March 2001                                 !!
!!           Refinements and modularization by                               !!
!!           Derek B. Leinweber 23 May 2007                                  !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE minf( funct, x, e, n, f, escale, iprint, icon, maxit, iterc, &
                 & w, nwork)
    USE kinds
    IMPLICIT NONE
    INTEGER , INTENT(IN)                      :: n
    REAL(DP), INTENT(INOUT), DIMENSION(n)     :: x
    REAL(DP), INTENT(IN)   , DIMENSION(n)     :: e
    REAL(DP), INTENT(OUT)                     :: f
    INTEGER , INTENT(IN)                      :: nwork
    REAL(DP), INTENT(OUT)  , DIMENSION(nwork) :: w
    INTEGER , INTENT(IN)                      :: icon, iprint, maxit
    INTEGER , INTENT(OUT)                     :: iterc
    REAL(DP), INTENT(IN)                      :: escale
    INTERFACE
       SUBROUTINE funct( n , x , f )
         USE kinds
         IMPLICIT none
         INTEGER               , INTENT(IN)  :: n
         REAL(DP), DIMENSION(n), INTENT(IN)  :: x
         REAL(DP)              , INTENT(OUT) :: f
       END SUBROUTINE funct
    END INTERFACE
    REAL(DP) :: a, aaa, b, d, da, dacc, db, dc, dd, ddmag, ddmax, di
    REAL(DP) :: dl, dmag, dmax, fa, fb, fc, fhold, fi, fkeep, fp, fprev
    REAL(DP) :: fstore, scer, sumval
    INTEGER  :: i, idirn, iline, ind, inn, is, isgrad, itone, ixp, j
    INTEGER  :: jil, jj, jjj, k, nfcc

    ddmag = 0.1_dp  * escale
    scer  = 0.05_dp / escale
    jj    = n * n + n
    jjj   = jj + n
    k     = n + 1
    nfcc  = 1
    ind   = 1
    inn   = 1

    do i = 1 , n
       do j = 1 , n
          w(k) = 0.0_dp
          if ( i == j ) then
             w(k) =  abs(e(i))
             w(i) = escale
          end if
          k = k + 1
       end do
    end do

    iterc  = 1
    isgrad = 2
    CALL funct( n , x , f )
    fkeep  = abs( 2.0_dp * f )
5   itone  = 1
    fp     = f
    sumval = 0.0_dp
    ixp    = jj

    do i = 1 , n
       ixp    = ixp + 1
       w(ixp) = x(i)
    end do

    idirn = n + 1
    iline = 1

7   dmax = w( iline )
    dacc = dmax * scer
    dmag = min( ddmag ,  0.1_dp * dmax )
    dmag = max( dmag  , 20.0_dp * dacc )
    ddmax = 10.0_dp * dmag

    if ( itone > 2 ) GOTO 71
    dl = 0.0_dp
    d = dmag
    fprev = f
    is = 5
    fa = f
    da = dl
8   dd = d - dl
    dl = d
58  k = idirn

    do i = 1 , n
       x(i) = x(i) + dd * w(k)
       k    = k + 1
    end do

    fstore = f
    call funct( n , x , f )
    nfcc = nfcc + 1

    if      ( is == 1 ) then
       GOTO 10
    else if ( is == 2 ) then
       GOTO 11
    else if ( is == 3 ) then
       GOTO 12
    else if ( is == 4 ) then
       GOTO 13
    else if ( is == 5 ) then
       GOTO 14
    else if ( is == 6 ) then
       GOTO 96
    end if

14  if ( f < fa ) then
       GOTO 15
    else if ( f == fa ) then
       GOTO 16
    else
       GOTO 24
    end if

16  if( abs(d) .le. dmax ) then
       d = d + d
       GOTO 8
    end if
    write(* , 19)
    GOTO 20
15  fb = f
    db = d
    GOTO 21
24  fb = fa
    db = da
    fa = f
    da = d
21  if( isgrad .le. 1 ) GOTO 83
23  d = db + db - da
    is = 1
    GOTO 8
83  d = 0.5_dp * ( da + db - (fa-fb)/(da-db) )
    is = 4
    if( ((da-d)*(d-db)) .ge. 0.0 ) GOTO 8
25  is = 1
    if( abs(d-db) .le. ddmax ) GOTO 8
26  d = db + sign( ddmax , db-da )
    is = 1
    ddmax = ddmax + ddmax
    ddmag = ddmag + ddmag
    if( ddmax > dmax ) then
       ddmax = dmax
    end if
    GOTO 8

13  if( f .ge. fa ) GOTO 23
28  fc = fb
    dc = db
29  fb = f
    db = d
    GOTO 30
12  if( f .le. fb ) GOTO 28
    fa = f
    da = d
    GOTO 30
11  if( f .ge. fb ) GOTO 10
    fa = fb
    da = db
    GOTO 29
71  dl = 1.0_dp
    ddmax = 5.0_dp
    fa = fp
    da = -1.0_dp
    fb = fhold
    db = 0.0_dp
    d = 1.0_dp
10  fc = f
    dc = d
30  a = (db-dc) * (fa-fc)
    b = (dc-da) * (fb-fc)
    if( ((a+b)*(da-dc)) > 0.0 ) GOTO 34
    fa = fb
    da = db
    fb = fc
    db = dc
    GOTO 26
34  d = 0.5_dp * ( a * (db+dc) + b * (da+dc) ) / (a+b)
    di = db
    fi = fb
    if( fb .le. fc ) GOTO 44
    di = dc
    fi = fc
44  if( itone .le. 2) GOTO 86
    itone = 2
    GOTO 45
86  if( abs(d-di)       .le. dacc               ) GOTO 41
    if( abs(d-di)       .le. (0.03_dp * abs(d)) ) GOTO 41
    if( abs(f-fstore)   .le. (0.01_dp * abs(f)) ) GOTO 41
45  if( ((da-dc)*(dc-d)) <   0.0                ) GOTO 47
    fa = fb
    da = db
    fb = fc
    db = dc
    GOTO 25
47  is = 2
    if( ((db-d)*(d-dc)) .ge. 0 ) GOTO 8
    is = 3
    GOTO 8
41  f = fi
    d = di - dl
    dd = sqrt( (dc-db) * (dc-da) * (da-db) / (a+b) )

    do i = 1 , n
       x(i)     = x(i) + d * w(idirn)
       w(idirn) =       dd * w(idirn)
       idirn    = idirn + 1
    end do

    w(iline) = w(iline) / dd
    iline    = iline + 1

    if( iprint /= 1) GOTO 51

50  write(* , 52) iterc , nfcc , f , (x(i) , i=1,n)

    if( iprint               > 1   ) GOTO 53
51  if( itone                > 1   ) GOTO 38
    if( (fprev - f - sumval) < 0.0 ) GOTO 94

    sumval = fprev - f
    jil    = iline

94  if( idirn .le. jj ) GOTO 7
    if( ind    >   1  ) GOTO 72

92  fhold = f
    is  = 6
    ixp = jj

    do i = 1 , n
       ixp    = ixp  + 1
       w(ixp) = x(i) - w(ixp)
    end do

    dd = 1.0_dp
    GOTO 58

96  if( ind  >   1 ) GOTO 87
    if( fp  .le. f ) GOTO 37

    d = 2.0_dp * (fp + f - 2.0_dp * fhold) / (fp - f)**2
    if( (d * (fp - fhold - sumval)**2 - sumval) .ge. 0.0 ) GOTO 37
87  j = jil * n + 1
    if( j > jj ) GOTO 61

    do i = j , jj
       k    = i - n
       w(k) = w(i)
    end do

    do i = jil , n
       w(i-1) = w(i)
    end do

61  idirn = idirn - n
    itone = 3
    k     = idirn
    ixp   = jj
    aaa   = 0.0_dp

    do i = 1 , n
       ixp  = ixp + 1
       w(k) = w(ixp)
       if( ( aaa - abs( w(k)/e(i) ) ) < 0.0 ) then
          aaa = abs( w(k)/e(i) )
       end if
       k = k + 1
    end do

    ddmag = 1.0_dp
    w(n)  = escale / aaa
    iline = n
    GOTO 7

37  ixp = jj
    aaa = 0.0_dp
    f   = fhold

    do i = 1 , n
       ixp  = ixp + 1
       x(i) = x(i) - w(ixp)
       if( (aaa * abs( e(i) ) - abs( w(ixp) ) ) < 0.0 ) then
          aaa = abs( w(ixp)/e(i) )
       end if
    end do

    GOTO 72

38  aaa = aaa * (1.0_dp + di)
    if( ind     >   1      ) GOTO 106
72  if( iprint .ge. 2      ) GOTO  50
53  if( ind     >   1      ) GOTO  88
    if( aaa     >   0.1_dp ) GOTO  76
    if( icon   .le. 1      ) GOTO  20
    ind = 2
    if( inn     >   1      ) GOTO 101

    inn = 2
    k = jjj

    do i = 1 , n
       k    = k + 1
       w(k) = x(i)
       x(i) = x(i) + 10.0_dp * e(i)
    end do

    fkeep = f
    CALL funct(n,x,f)
    nfcc  = nfcc+1
    ddmag = 0.0_dp
    GOTO 108

76  if( f < fp ) GOTO 35

78  write(*,80)
    GOTO 20

88  ind    = 1
35  ddmag  = 0.4_dp * sqrt( fp - f )
    isgrad  = 1
108 iterc = iterc + 1

    if( iterc .le. maxit ) GOTO 5

    if( f .le. fkeep ) GOTO 20

    f = fkeep

    do i = 1 , n
       jjj  = jjj + 1
       x(i) = w(jjj)
    end do

    GOTO 20

101 jil = 1
    fp = fkeep
    if(      f  < fkeep ) then
       GOTO 105
    else if( f == fkeep ) then
       GOTO 78
    end if

    jil = 2
    fp = f
    f  = fkeep
105 ixp = jj

    do i = 1 , n
       ixp    = ixp + 1
       k      = ixp + n
       if( jil > 1 ) then
          w(ixp) = x(i)
          x(i)   = w(k)
       else
          w(ixp) = w(k)
       end if
    end do

    jil = 2
    GOTO 92
106 if( aaa > 0.1_dp ) GOTO 107
20  RETURN
107 inn = 1
    GOTO 35

    !! The FORMAT statements
19  FORMAT( 5x , 'minf maximum change does not alter the function value.' )
52  FORMAT(/1x , 'iteration' , i5 , i15 , ' function values' ,   &
         &        6x , 'f =' , es21.14/(3es24.14) )
80  FORMAT( 5x , 'minf accuracy is limited by errors in function value.' )

  END SUBROUTINE minf

END MODULE MinFun
