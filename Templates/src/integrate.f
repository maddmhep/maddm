c Duplicates (sets the array element to 10 for now)
      subroutine Duplicates(alist, length)
         double precision alist(length)
         integer i, j

         i = 1
         do while (alist(i).le.1.d0)
           j = i+1
           do while (j.le.length)
                if (alist(i).eq.alist(j)) then
                 alist(j) = 10.d0
                endif
                j = j+1
           enddo
           i = i+1
          enddo
       end subroutine Duplicates
c Bubblesort

        SUBROUTINE Bubble_Sort(a, length)
          double precision a(length)
          double precision :: temp
          INTEGER :: i, j
          LOGICAL :: swapped

          DO j = length-1, 1, -1
            swapped = .FALSE.
            DO i = 1, j
              IF (a(i) > a(i+1)) THEN
                temp = a(i)
                a(i) = a(i+1)
                a(i+1) = temp
                swapped = .TRUE.
              END IF
            END DO
            IF (.NOT. swapped) EXIT
          END DO
        END SUBROUTINE Bubble_Sort

c  coordinate transformation
c      function eps_trans(beta1, beta0, gamma_beta)
c         implicit none
c         include 'maddm.inc'
c         include 'coupl.inc'

c         double precision beta1, beta0, gamma_beta
c
c        eps_trans = atan((beta1*beta1 - beta0*beta0)/(beta0*gamma_beta))
c         return
c      end

c  beta
c      function beta_fun(eps, beta0, gamma_beta)
c         implicit none
c         include 'maddm.inc'
c         include 'coupl.inc'
c
c         double precision eps, beta0, gamma_beta
c         beta_fun = sqrt(beta0**2+tan(eps)*beta0*gamma_beta)
c      return
c      end

c  jacobian
c      function jacobian(beta1, beta0, gamma_beta)
c         implicit none
c         include 'maddm.inc'
c         include 'coupl.inc'
c         double precision beta1, beta0, gamma_beta
c         jacobian = 2*beta1/(beta0*(1.d0+( (beta1**2 - beta0**2)**2/ (beta0**2*gamma_beta**2) )*gamma_beta*beta0  ))
c         return
c      end

c-------------------------------------------------------------------------c
c set up integration grid for relic density. First initialize the whole fix
c length array (about 2000 element, see maddm.inc for more details) to 10.
c then go through from the first one and set up the part from a to b.
c in all cases b is maximally 1. Then add resonance positions after 1.
c Then sort the array. The simpson routine will know to use only the part
c from 0 to 1 omitting all the remaining 10s.
c-------------------------------------------------------------------------c
      subroutine set_up_grid(a, b)
c-------------------------------------------------------------------------c

      implicit none
      include 'maddm.inc'
      include 'coupl.inc'

      double precision  a,b, func, additional_pt, exponent
      double precision left, right, whole, end_pt, start_pt, width

      integer ii, jj, kk, grid_pos
      integer nres_points_local,max_width_factor
      

c    Set up the integration grid
c    Initialize to 10*b, so that we can identify the relevant points
c    by going only to b in the integration grid and neglecting the rest.
          grid_pos = 1
          grid=10.d0*b

c         write(*,*) 'width: ', beta_res_width(1)
          call find_resonances()
c    Make sure that the array is large enough to store the grid
         if (grid_npts.lt.ngrid_init) then
             write(*,*) 'Error: grid array not large enough!'
             call exit(1)
         endif

         do ii=1, ngrid_init+1
             grid(grid_pos) = (b-a)/ngrid_init*(ii-1)
             grid_pos=grid_pos+1
         enddo

c add the resonance positions first
       do ii=1, nres
            if (beta_res(ii).ge.a.and.beta_res(ii).lt.b) then
                grid(grid_pos) = beta_res(ii)
                grid_pos = grid_pos+1
            endif
       enddo
       
       nres_points_local = 20
       max_width_factor = 2 ! allow to go as far as e^N-1 times the width for the extra point

c then add more points around the resonance
       do ii=1, nres
          if (beta_res(ii).ge.a.and.beta_res(ii).lt.b.and.beta_res(ii).ge.0.d0) then
              do jj=1, nres_points_local
c                pts_to_add_adaptive = ceiling(real(pts_to_add_adaptive / 2))
c                do kk = 1, pts_to_add_adaptive
                    additional_pt = beta_res(ii) + beta_res_width(ii)*(DEXP(1d0*jj**max_width_factor/nres_points_local)-1)
c                    additional_pt = beta_res(ii) + (exp(real(beta_res_width(ii))/10.d0*exp(real(jj)))-1.d0)
                    if (additional_pt.lt.b) then
                        grid(grid_pos) = additional_pt
                        grid_pos=grid_pos+1
                        if (grid_npts.lt.grid_pos) then
                            write(*,*) 'Error: grid array not large enough!'
                            call exit(1)
                        endif
                    endif

                    additional_pt = beta_res(ii) - beta_res_width(ii)*(DEXP(1d0*jj**max_width_factor/nres_points_local)-1)
                    if (additional_pt.gt.a) then
                        grid(grid_pos) = additional_pt
                        grid_pos=grid_pos+1
                        if (grid_npts.lt.grid_pos) then
                            write(*,*) 'Error: grid array not large enough!'
                            call exit(1)
                        endif
                    endif
c                enddo
              enddo
           endif
       enddo



c remove duplicates
c sort the grid
       call Duplicates(grid, grid_npts)
       call  Bubble_Sort(grid, grid_npts)

c       write(*,*) 'grid:', grid_pos
c       do ii=1, grid_pos
c           write(*,*) grid(ii)
c       enddo

       end subroutine set_up_grid
c------------------------------------------------------------

c-------------------------------------------------------------------------c
c 1-D function integration using the Simpson's method - OLD!!!
c-------------------------------------------------------------------------c
c      function simpson_taacs(func, a,b, grid, grid_npts)

c         implicit none
c         include 'maddm.inc'
c         include 'coupl.inc'

c         external func


c         include 'resonances.inc'


c       simpson_taacs = simpson(func,0.d0, 1.d0, grid, grid_npts)
c       return

c       end function simpson_taacs

c-------------------------------------------------------------------------c
c 1-D function integration - SIMPSON'S RULE, integrates from 0 to 1
c-------------------------------------------------------------------------c
      function simpson(func,  a, b, grid_simp, gr_npts)

      implicit none
      include 'maddm.inc'
      include 'coupl.inc'

      double precision func, exponent
      double precision  whole, end_pt, start_pt,a, b
      external func
      integer gr_npts, ii
      double precision grid_simp(gr_npts)
      double precision tmp1, tmp2, tmp3
      simpson = 0.d0
      ii = 1

      if (a.lt.grid_simp(1).or.b.gt.grid_simp(gr_npts-1) ) then
           write(*,*) 'Error: Simpson integration bounds outside the grid.'
           return
      endif

      do while (grid_simp(ii+1).le.b)
           start_pt = grid_simp(ii)
           if (start_pt.lt.a) then
               ii= ii+1
               continue
           endif
           end_pt =  grid_simp(ii+1)
c            write(101,*) start_pt, func(start_pt)
c           write(*,*) 'start, end', start_pt, end_pt
c           write(*,*) 'func: ', func(0.d0), func(0.00001d0)

           tmp1 = (func(start_pt) + 4.d0*func(0.5d0*
     .                                  (start_pt+end_pt)) + func(end_pt) )
           if (tmp1.lt.1e-20*simpson)then
              exit
           else if (tmp1.gt.1d-300) then
           simpson = simpson+ (end_pt - start_pt)/6.d0*tmp1
           else 
              exit
           endif
c            simpson = simpson+ (end_pt - start_pt)/8.d0*(max(0.d0, func(start_pt))
c     .             + 3.d0*max(0.d0,func(0.33d0*(2.d0*start_pt+end_pt)))
c     .             + 3.d0*max(0.d0,func(0.33d0*(start_pt+2.d0*end_pt)))  + max(0.d0,func(end_pt)))
           ii = ii+1

      enddo
      return

      end function simpson

c-------------------------------------------------------------------------c
c 1-D function integration
c-------------------------------------------------------------------------c
      subroutine romberg(func, a, b, integral, eps, iter)
c-------------------------------------------------------------------------c
c                                                                         c
c  Returns the integral of the function func from a to b.  Integration    c
c  is performed by Romberg's method of order 2K, where, e.g. K=2 is       c
c  Simpson's rule.  Based on numerical recipe's qromb.                    c
c                                                                         c
c  Parameters:                                                            c
c  func - the user supplied double precision function                     c
c  a - lower bound of the integral                                        c
c  b - upper bound of the integral                                        c
c  integral - result returned by romberg                                  c
c  eps - the fractional accuracy desired as determined by the             c
c        extrapolation error estimate                                     c
c  iter - the minimum number of iterations (to improve accuracy)          c
c                                                                         c
c  jmax - limits the total number of steps                                c
c  k - the number of points used in the extrapolation.                    c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none

c input parameters
      double precision func, a, b, integral, eps
      integer iter
      external func

c parameters for this subroutine only
      integer jmax, jmaxp, k, km
      parameter (jmax=20, jmaxp=jmax+1, k=5, km=k-1)

c uses polint, trapzd
      integer j

c these store the successive trapezoidal approximations and their
c relative stepsizes
      double precision dss, h(jmaxp), s(jmaxp)

      h(1) = 1.d0
      do j=1,jmax
        call trapzd(func,a,b,s(j),j)
        if ((j.ge.k).and.(j.ge.iter)) then
          call polint(h(j-km),s(j-km),k,0.d0,integral,dss)
          if (dabs(dss).le.eps*dabs(integral)) return
        endif
        s(j+1)=s(j)

c This is a key step: The factor is 0.25 even though the stepsize decrease is only 0.5
c This makes the extrapolation a polynomial in h^2.
        h(j+1)=0.25d0*h(j)
      enddo

c      pause 'too many steps in romberg'
      end



c-------------------------------------------------------------------------c
      subroutine trapzd(func, a, b, s, n)
c-------------------------------------------------------------------------c
c                                                                         c
c  This routine computes the nth stage of refinement of an extended       c
c  trapezoidal rule.  func is input as the name of the function to be     c
c  integrated between limits a and b, also input.  When called with n=1,  c
c  the routine returns as s the crudent estimate of the integral.         c
c  Subsequent calls with n=2,3,... will improve the accuracy of s by      c
c  adding 2^(n-2) additional interior points.  s should not be modified   c
c  between sequential calls.                                              c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none

      integer n
      double precision a, b, s, func
      external func
      integer it, j
      double precision del, summ, tnm, x

      if (n.eq.1) then
        s = 0.5d0*(b-a)*(func(a)+func(b))
      else
        it = 2**(n-2)
        tnm = it

c This is the spacing of the points to be added
        del = (b-a)/tnm
        x = a+0.5d0*del
        summ = 0.d0
        do j=1,it
          summ = summ+func(x)
          x = x+del
        enddo

c This replaces s by its refined value
        s = 0.5d0*(s+(b-a)*summ/tnm)
      endif
      return
      end



c-------------------------------------------------------------------------c
      subroutine polint(xa, ya, n, x, y, dy)
c-------------------------------------------------------------------------c
c                                                                         c
c  Given arrays xa and ya, each of length n, and given a value x, this    c
c  routine returns a value y, and an error estimate dy.  If P(x) is the   c
c  polynomial of degree N=1 such that P(xa_i) = ya_i, i=1...n then the    c
c  returned value y = P(x).                                               c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none

      integer n, nmax
      double precision dy, x, y, xa(n), ya(n)
      parameter (nmax=10)
      integer i, m, ns
      double precision den, dif, dift, ho, hp, w, c(nmax), d(nmax)

      ns = 1
      dif = dabs(x-xa(1))

c Here we find the index ns of the closest table entry
      do i=1,n
        dift = dabs(x-xa(i))
        if (dift.lt.dif) then
          ns = i
          dif = dift
        endif

c and initialize the tableau of c's and d's
        c(i) = ya(i)
        d(i) = ya(i)
      enddo

c This is the initial approximation to y
      y = ya(ns)
      ns = ns-1

c For each column of the tableau,
c we loop over the current c's and d's and update them
      do m=1,n-1
        do i=1,n-m
          ho = xa(i)-x
          hp = xa(i+m)-x
          w = c(i+1)-d(i)
          den = ho-hp

c This error can occur only if two input xa's are (to within roundoff) identical.
c          if (den .eq. 0.d0) pause 'failure in polint'
          den = w/den

c Here the c's and d's are updated
          d(i) = hp*den
          c(i) = ho*den
        enddo
        if(2*ns .lt. n-m) then
          dy = c(ns+1)
        else
          dy = d(ns)
          ns = ns-1
        endif
        y = y+dy
      enddo
      return
      end



c-------------------------------------------------------------------------c
c Ordinary Differential Equation integrator
c-------------------------------------------------------------------------c
      subroutine odeint(ystart, nvar, x1, x2, eps, h1, hmin, nok, nbad, derivs, rkqs)
c-------------------------------------------------------------------------c
c                                                                         c
c  Runge-Kutta driver with adaptive stepsize control.  Integrate the      c
c  starting values ystart(1:nvar) from x1 to x2 with accuracy eps,        c
c  storing intermediate results in the common block /path/.  h1 should    c
c  be set as a guessed first stepsize, hmin as the minimum allowed        c
c  stepsize (can be zero). On output nok and nbad are the number of       c
c  good and bad (but retried and fixed) steps taken, and ystart is        c
c  replaced by values at the end of the integration interval.  derivs is  c
c  the user-supplied subroutine to be used. /path/ contains its own info  c
c  about how often an intermediate value is stored.                       c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      integer nbad, nok, nvar, kmaxx, maxstp, nmax
      double precision eps, h1, hmin, x1, x2, ystart(nvar), tiny
      external derivs, rkqs
      parameter (maxstp=100000, nmax=50, kmaxx=200, tiny=1.0d-30)
      integer i, kmax, kount, nstp
      double precision dxsav, h, hdid, hnext, x, xsav, dydx(nmax), xp(kmaxx), y(nmax),
     &      yp(nmax,kmaxx), yscal(nmax)
      common /path/ kmax, kount, dxsav, xp, yp
      x = x1
      h = dsign(h1,x2-x1)
      nok = 0
      nbad = 0
      kount = 0
      do i=1,nvar
        y(i)=ystart(i)
      enddo
      if (kmax .gt. 0) xsav = x-2.d0*dxsav
      do nstp=1,maxstp
        call derivs(x,y,dydx)
        do i=1,nvar
          yscal(i)=dabs(y(i))+dabs(h*dydx(i))+TINY
        enddo
        if (kmax .gt. 0) then
          if (dabs(x-xsav) .gt. dabs(dxsav)) then
            if(kount .lt. kmax-1) then
              kount = kount+1
              xp(kount) = x
              do i=1,nvar
                yp(i,kount) = y(i)
              enddo
              xsav = x
            endif
          endif
        endif
        if((x+h-x2)*(x+h-x1) .gt. 0.0d0) h = x2-x
        call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
        if(hdid.eq.h) then
          nok = nok+1
        else
          nbad = nbad+1
        endif
        if((x-x2)*(x2-x1) .ge. 0.0d0) then
          do i=1,nvar
            ystart(i) = y(i)
          enddo
          if(kmax.ne.0) then
            kount = kount+1
            xp(kount) = x
            do i = 1,nvar
              yp(i,kount) = y(i)
            enddo
          endif
          return
        endif
        if(dabs(hnext) .lt. hmin) then
           write(*,*) 'too small step'
           stop 1
        endif 
        h = hnext
        enddo
c      pause 'too many steps in odeint'
      return
      end



c-------------------------------------------------------------------------c
      subroutine rkck(y, k1, n, x, h, yout, yerr, derivs)
c-------------------------------------------------------------------------c
c                                                                         c
c  This subroutine calculates a step in solving an ODE using the          c
c  Runga-Kutta-Cash-Karp method.                                          c
c                                                                         c
c  derivs = external subroutine for calculating derivatives               c
c  n = number of ODEs being solved                                        c
c  y = inital point y in the step                                         c
c  x = initial point x in the step                                        c
c  h = stepsize                                                           c
c  yout = output of the new y-value                                       c
c  yerr = estimated error of the step                                     c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      include 'maddm.inc'

c input parameters
      integer n
      double precision x,y(n),h,yout(n),yerr(n)
      external derivs

c parameters for this subroutine only
      integer i
      double precision k1(n),k2(n),k3(n),k4(n),k5(n),k6(n),ytemp(n)
      double precision c2,c3,c4,c5,c6
      double precision a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65
      double precision b1,b2,b3,b4,b5,b6
      double precision b1s,b2s,b3s,b4s,b5s,b6s
      parameter (c2=1.d0/5.d0,c3=3.d0/10.d0,c4=3.d0/5.d0,c5=1.d0,c6=7.d0/8.d0,a21=1.d0/5.d0,
     . a31=3.d0/40.d0,a32=9.d0/40.d0,a41=3.d0/10.d0,a42=-9.d0/10.d0,a43=6.d0/5.d0,
     . a51=-11.d0/54.d0,a52=5.d0/2.d0,a53=-70.d0/27.d0,a54=-35.d0/27.d0,
     . a61=1631.d0/55296.d0,a62=175.d0/512.d0,a63=575.d0/13824.d0,a64=44275.d0/110592.d0,
     . a65=253.d0/4096.d0,b1=37.d0/378.d0,b2=0.d0,b3=250.d0/621.d0,b4=125.d0/594.d0,
     . b5=0.d0,b6=512.d0/1771.d0,b1s=2825.d0/27648.d0,b2s=0.d0,b3s=18575.d0/48384.d0,
     . b4s=13525.d0/55296.d0,b5s=277.d0/14336.d0,b6s=1.d0/4.d0)

c many setup calculations are required
c      call derivs(x,y,k1)
      do i=1,n
        ytemp(i) = y(i) + h*a21*k1(i)
      enddo
      call derivs(x+h*c2,ytemp,k2)
      do i=1,n
        ytemp(i) = y(i) + h*(a31*k1(i) + a32*k2(i))
      enddo
      call derivs(x+h*c3,ytemp,k3)
      do i=1,n
        ytemp(i) = y(i) + h*(a41*k1(i) + a42*k2(i) + a43*k3(i))
      enddo
      call derivs(x+h*c4,ytemp,k4)
      do i=1,n
        ytemp(i) = y(i) + h*(a51*k1(i) + a52*k2(i) + a53*k3(i) + a54*k4(i))
      enddo
      call derivs(x+h*c5,ytemp,k5)
      do i=1,n
        ytemp(i) = y(i) + h*(a61*k1(i) + a62*k2(i) + a63*k3(i) + a64*k4(i) + a65*k5(i))
      enddo
      call derivs(x+h*c6,ytemp,k6)

c time to calculate the result and error
      do i=1,n
C         if ((found_xf(i)).and.((x+h).lt.xf(i))) then
C           g_star_S = get_gstar(mdm(1)/(x+h))
C           yout(i) = get_Y_eq(i,(x+h),g_star_S)
C           yerr(i) = 0.d0
C         else
          yout(i) = y(i) + h*(b1s*k1(i) + b2s*k2(i) + b3s*k3(i) + b4s*k4(i) + b5s*k5(i) + b6s*k6(i))
          yerr(i) = h*((b1-b1s)*k1(i) + (b2-b2s)*k2(i) + (b3-b3s)*k3(i) + (b4-b4s)*k4(i) + (b5-b5s)*k5(i) + (b6-b6s)*k6(i))
C         endif
      enddo

      return
      end



c-------------------------------------------------------------------------c
      subroutine rkqs(y, dydx, n, x, htry, eps, yscal, hdid, hnext, derivs)
c-------------------------------------------------------------------------c
c                                                                         c
c  Fifth-order Runge-Kutta step with monitoring of local truncation       c
c  error to ensure accuracy and adjust stepsize. Input are the dependent  c
c  variable y vector(1:n) and its derivative dydx(1:n) at the starting    c
c  value of the independent variable x.  Also input are the stepsize to   c
c  be attempted htry, the required accuracy eps, and the vector           c
c  yscal(1:n) against which the error is scaled.  On output, y and x are  c
c  replaced by their new values, hdid is the stepsize that was actually   c
c  accomplished, and hnext is the estimated next stepsize.  derivs is     c
c  the user supplied subroutine that computes the right-hand side         c
c  derivatives.                                                           c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none

      integer n,nmax
      double precision eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
      external derivs

c     Maximum number of equations
      parameter (nmax=50)
      integer i
      double precision errmax,h,htemp,xnew,yerr(nmax),ytemp(nmax),safety,pgrow,
     &  pshrink,errcon
      parameter (safety=0.9d0,pgrow=-.2d0,pshrink=-.25d0,errcon=1.89d-4)
      h = htry
    1 call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)
      errmax=0.0d0
      do i=1,n
        errmax = max(errmax,dabs(yerr(i)/yscal(i)))
      enddo
      errmax = errmax/eps
      if(errmax .gt. 1.0d0) then
        htemp = safety*h*(errmax**pshrink)
        h = dsign(max(dabs(htemp), 0.1d0*dabs(h)), h)
        xnew = x+h
c        if (xnew.eq.x) pause 'stepsize underflow in rkqs'
        goto 1
      else
        if(errmax.gt.ERRCON) then
          hnext = safety*h*(errmax**pgrow)
        else
          hnext = 5.0d0*h
        endif
        hdid = h
        x = x+h
        do i=1,n
          y(i) = ytemp(i)
        enddo
        return
      endif
      end



c-------------------------------------------------------------------------c
c Monte Carlo Integration
c-------------------------------------------------------------------------c
      subroutine vegas(region, ndim, fxn, init, ncall, itmx,
     . nprn, tgral, sd, chi2a)
c-------------------------------------------------------------------------c
c                                                                         c
c  Monte Carlo integration of external function fxn                       c
c  region - 2*ndim vector of the boundaries of integration                c
c  ndim - number of dimensions of the integration                         c
c  fxn - external function to integrate over                              c
c  init - new grid (0), prev grid w/o (1) or w/ (2) results               c
c  ncall - number of times to call the function                           c
c  itmx - number of iterations (5)                                        c
c  nprn - flag for controlling output: nprn < 0 for nothing               c
c  tgral - best estimate of the integral                                  c
c  sd - standard deviation                                                c
c  chi2a - chi squared per degree of freedom                              c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none

      integer init, itmx, ncall, ndim, nprn, ndmx, mxdim
      double precision tgral, chi2a, sd, region(2*ndim), fxn, alph, tiny
      parameter (alph = 1.5d0, ndmx = 50, mxdim = 10, tiny = 1.d-30)
      external fxn

c uses fxn, random, rebin
      integer i,idum,it,j,k,mds,nd,ndo,ng,npg,ia(mxdim),kg(mxdim)
      double precision calls,dv2g,dxg,f,f2,f2b,fb,rc,ti,tsi,
     . wgt,xjac,xn,xnd,xo,d(ndmx,mxdim),di(ndmx,mxdim),dt(mxdim),
     . dx(mxdim),r(ndmx),x(mxdim),xi(ndmx,mxdim),xin(ndmx),random
      double precision schi,si,swgt

c Means for random number initialization
      common /ranno/ idum
      save

c Normal entry. Enter here on cold start
      if (init.le.0) then

c Change mds = 0 to disable stratified sampling (i.e. use importance
c sampling only)
        mds = 1
        ndo = 1
        do j=1,ndim
          xi(1,j) = 1.d0
        enddo
      endif

c Enter here to inherit the grid from a previous call, but not its answers
      if (init.le.1) then
        si = 0.d0
        swgt = 0.d0
        schi = 0.d0
      endif

c Enter here to inherit the previous grid and its answers
      if (init.le.2) then
        nd = ndmx
        ng = 1

c Set up for stratification
        if (mds.ne.0) then
          ng = int((dble(ncall)/2.d0+0.25d0)**(1.d0/dble(ndim)))
          mds = 1
          if ((2*ng-ndmx).ge.0) then
            mds = -1
            npg = ng/ndmx + 1
            nd = ng/npg
            ng = npg*nd
          endif
        endif
        k = ng**ndim
        npg = max(ncall/k,2)
        calls = dble(npg)*dble(k)
        dxg = 1.d0/dble(ng)
        dv2g = (calls*dxg**ndim)**2/dble(npg)/dble(npg)/(dble(npg)-1.d0)
        xnd = dble(nd)
        dxg = dxg*xnd
        xjac = 1.d0/calls
        do j=1,ndim
          dx(j) = region(j+ndim)-region(j)
          xjac = xjac*dx(j)
        enddo

c Do binning if necessary
        if (nd.ne.ndo) then
          do i=1,max(nd,ndo)
            r(i) = 1.d0
          enddo
          do j=1,ndim
            call rebin(ndo/xnd,nd,r,xin,xi(1,j))
          enddo
          ndo = nd
        endif
        if (nprn.ge.0) write(*,200) ndim,calls,it,itmx,nprn,ALPH,mds,nd,
     *(j,region(j),j,region(j+ndim),j=1,ndim)
      endif

c Main iteration loop
      do it=1,itmx
        ti = 0.d0
        tsi = 0.d0
        do j=1,ndim
          kg(j) = 1
          do i=1,nd
            d(i,j) = 0.d0
            di(i,j) = 0.d0
          enddo
        enddo
10      continue
        fb = 0.d0
        f2b = 0.d0
        do k=1,npg
          wgt = xjac
          do j=1,ndim
            xn = (kg(j)-random(idum))*dxg+1.d0
            ia(j) = max(min(int(xn),NDMX),1)
            if (ia(j).gt.1) then
              xo = xi(ia(j),j)-xi(ia(j)-1,j)
              rc = xi(ia(j)-1,j)+(xn-ia(j))*xo
            else
              xo = xi(ia(j),j)
              rc = (xn-ia(j))*xo
            endif
            x(j) = region(j)+rc*dx(j)
            wgt = wgt*xo*xnd
          enddo
          f = wgt*fxn(x,wgt)
          f2 = f*f
          fb = fb+f
          f2b = f2b+f2
          do j=1,ndim
            di(ia(j),j) = di(ia(j),j)+f
            if (mds.ge.0) d(ia(j),j)=d(ia(j),j)+f2
          enddo
        enddo
        f2b = dsqrt(f2b*npg)
        f2b = (f2b-fb)*(f2b+fb)
        if (f2b.le.0.d0) f2b=tiny
        ti = ti+fb
        tsi = tsi+f2b

c Use stratified sampling
        if (mds.lt.0) then
          do j=1,ndim
            d(ia(j),j) = d(ia(j),j)+f2b
          enddo
        endif
        do k=ndim,1,-1
          kg(k) = mod(kg(k),ng)+1
          if (kg(k).ne.1) goto 10
        enddo

c Compute final results with this iteration
        tsi = tsi*dv2g
        wgt = 1.d0/tsi
        si = si+dble(wgt)*dble(ti)
        schi = schi+dble(wgt)*dble(ti)**2
        swgt = swgt+dble(wgt)
        tgral = si/swgt
        chi2a = max((schi-si*tgral)/(it-.99d0),0.d0)
        sd=dsqrt(1.d0/swgt)
        tsi=dsqrt(tsi)
        if (nprn.ge.0) then
          write(*,201) it,ti,tsi,tgral,sd,chi2a
          if (nprn.ne.0) then
            do j=1,ndim
              write(*,202) j,(xi(i,j),di(i,j),i=1+nprn/2,nd,nprn)
            enddo
          endif
        endif

c Refine the grid
        do j=1,ndim
          xo = d(1,j)
          xn = d(2,j)
          d(1,j) = (xo+xn)/2.d0
          dt(j) = d(1,j)
          do i=2,nd-1
            rc = xo+xn
            xo = xn
            xn = d(i+1,j)
            d(i,j) = (rc+xn)/3.d0
            dt(j) = dt(j)+d(i,j)
          enddo
          d(nd,j) = (xo+xn)/2.d0
          dt(j) = dt(j)+d(nd,j)
        enddo
        do j=1,ndim
          rc = 0.d0
          do i=1,nd
            if (d(i,j).lt.tiny) d(i,j)=tiny
            r(i) = ((1.d0-d(i,j)/dt(j))/(dlog(dt(j))-dlog(d(i,j))))**ALPH
            rc = rc+r(i)
          enddo
          call rebin(rc/xnd,nd,r,xin,xi(1,j))
        enddo
      enddo

200   format(/' input parameters for vegas:  ndim=',i3,'  ncall=',
     *f8.0/28x,'  it=',i5,'  itmx=',i5/28x,'  nprn=',i3,'  alph=',
     *f5.2/28x,'  mds=',i3,'   nd=',i4/(30x,'xl(',i2,')= ',g11.4,' xu(',
     *i2,')= ',g11.4))
201   format(/' iteration no.',I3,': ','integral =',g14.7,'+/- ',g9.2/
     *' all iterations:   integral =',g14.7,'+/- ',g9.2,
     *' chi**2/it''n =',g9.2)
202   format(/' data for axis ',I2/'    x       delta i       ',
     *'   x       delta i       ','    x       delta i       ',/(1x,
     *f7.5,1x,g11.4,5x,f7.5,1x,g11.4,5x,f7.5,1x,g11.4))

      return
      end



c-------------------------------------------------------------------------c
      function random(idum)
c-------------------------------------------------------------------------c
c                                                                         c
c  Randomly generates a double precision number between 0 and 1.  First   c
c  set idum to a negative number to initialize the seed and afterward     c
c  just call random(idum)                                                 c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none

      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      double precision random,am,eps,rnmx
      parameter (im1=2147483563,im2=2147483399,am=1.d0/im1,imm1=im1-1,
     . ia1=40014, ia2=40692, iq1=53668, iq2=52774, ir1=12211, ir2=3791,
     . ntab=32, ndiv=1+imm1/ntab, eps=1.2d-7, rnmx=1.d0-eps)
      integer idum2, j, k, iv(ntab), iy
      save iv, iy, idum2
      data idum2/123456789/, iv/ntab*0/, iy/0/

c initialize
      if (idum.le.0) then

c prevents idum = 0
       idum = max(-idum,1)
       idum2 = idum

c Load the shuffle table (after 8 warm-ups)
       do j=ntab+8,1,-1
        k = idum/iq1
        idum = ia1*(idum-k*iq1)-k*ir1
        if (idum.lt.0) idum = idum + im1
        if (j.le.ntab) iv(j) = idum
       enddo
       iy=iv(1)
      endif

c Start here when not initializing
      k = idum/iq1

c Compute idum=mod(ia1*idum,im1) without overflows by Schrage's method
      idum = ia1*(idum-k*iq1)-k*ir1
      if (idum.lt.0) idum = idum + im1
      k = idum/iq2

c Compute idum2=mod(ia2*idum2,im2) likewise
      idum2 = ia2*(idum2-k*iq2)-k*ir2
      if(idum2.lt.0) idum2 = idum2 + im2

c Will be in the range 1:ntab
      j = 1+iy/ndiv

c Here idum is shuffled, idum and idum2 are combined to generate output
      iy = iv(j)-idum2
      iv(j) = idum
      if(iy.lt.1) iy = iy + imm1

c The endpoints are not expected in random number generators
      random = min(am*iy,rnmx)

      return
      end



c-------------------------------------------------------------------------c
      subroutine rebin(rc,nd,r,xin,xi)
c-------------------------------------------------------------------------c
c                                                                         c
c  Utility routine used by vegas to rebin a vector of densities xi into   c
c  new bins defined by a vector r                                         c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none

      integer nd
      double precision rc, r(*), xi(*), xin(*)
      integer i, k
      double precision dr, xn, xo
      k = 0
      xo = 0.d0
      dr = 0.d0
      do i=1,nd-1
1      if (rc.gt.dr) then
        k = k+1
        dr = dr+r(k)
        goto 1
       endif
       if (k.gt.1) xo = xi(k-1)
       xn = xi(k)
       dr = dr-rc
       xin(i) = xn-(xn-xo)*dr/r(k)
      enddo
      do i=1,nd-1
       xi(i) = xin(i)
      enddo
      xi(nd) = 1.d0
      return
      end
