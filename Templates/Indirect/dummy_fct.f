      subroutine get_dummy_x1(sjac, X1, R, pbeam1, pbeam2, stot, shat)

      implicit none
      include 'maxparticles.inc'
      include 'run.inc'
      include 'nexternal.inc'
c      include 'genps.inc'
      double precision sjac ! jacobian. should be updated not reinit
      double precision X1   ! bjorken X. output
      double precision R    ! random value after grid transfrormation. between 0 and 1
      double precision pbeam1(0:3) ! momentum of the first beam (input and/or output)
      double precision pbeam2(0:3) ! momentum of the second beam (input and/or output)
      double precision stot        ! total energy  (input and /or output)
      double precision shat        ! output

c     global variable to set (or not)
      double precision cm_rap
      logical set_cm_rap
      common/to_cm_rap/set_cm_rap,cm_rap
c     mass of initial particles
      double precision pmass(nexternal)
      common/to_mass/  pmass
      double precision m

      double precision pi
      data pi /3.1415926535897932d0/
      double precision S, T
      double precision thres
      data thres /1e-7/

      double precision ymin      
      double precision besselk
c     parameter for newton method
      integer i
      integer sqrts ! temporary value
      double precision f, df
c     The integral to compute is 
c     d\sqrt(S) \sart(s)^2 K1(\sqrt{s}/T)
c     to generate events we consider the approximate
c     d\sqrt(S) \sart(s)^2 exp[-\sqrt(S)/T]
c     This integral is know analytically:
c     -T* exp(-x/T) *(T^2 + (x+T)^2)
c     
c     Therefore we will apply the following change of variable:
c     y = -T exp(-x/T) *(T^2 + (x+T)^2))
c    The associated jacobian being:
c     dy = x^2 exp[-x/T] dx
c
c     Therefore the integral to compute is 
c
c     dy K1(x(y)/T)*exp(x(y)/T)
c
c     the bound of integration are
c     ymin = -T exp(-2m/T) (T^2 + (2m+T)^2)) 
c     ymax = 0 
c
c
c     In order to map this integral between 0 and 1, we then define
c     the following change of variable
c     R = (y - ymin)/(-ymin) = -y/ymin +1
c     dR = dy/(-ymin) = x^2 exp[-x/T]/ymin dx
      m =pmass(1)
      if (lpp(1).eq.9)then
         T = ebeam(1)
      else
         T = ebeam(2)
      endif
      if (pmass(1).ne.pmass(2))then
         stop 1
      endif
      ymin = -T * DEXP(-1*(2*m)/T) *(T**2+(2*m+T)**2)
c
c     Solve x=sqrts for the following equation
c     R = (y - ymin)/(-ymin) = (-T exp(-x/T) *(T^2 + (x+T)^2))-ymin)/(-ymin)
c
      sqrts = T
      do i=1,1000
         f = (-T*DEXP(-sqrts/T) *(T**2 + (sqrts+T)**2)-ymin)/(-ymin)-X1
         if (f.lt.thres) then
            exit
         endif
         df = sqrts**2 *DEXP(-sqrts/T)/(-ymin)
         sqrts=sqrts-f/df
      enddo
      if (f.gt.thres) then
         sqrts = 0d0
         S = 0d0
         sjac = -1d0
         return
      endif
c
c      setting some global variable/output variable
c
       stot =sqrts**2
       shat = stot
       cm_rap = 0d0
       set_cm_rap=.true.
       pbeam1(0)=dsqrt(stot)/2
       pbeam1(3)=sqrt(max(pbeam1(0)**2-m**2, 0d0))
       pbeam2(0)=pbeam1(0)
       pbeam2(3)=-sqrt(max(pbeam2(0)**2-m**2, 0d0))
       ebeam(1) = pbeam1(0)
       ebeam(2) = pbeam2(0)
       x1 = 1d0
c
c      Setting the Jacobian AND the associate PDF
c
c      normally we should do:
c      sjac = sjac * dabs(ymin /S * DEXP[sqrts/T])   
c      But we will multiply here by the PDF!
c      2*sqrts**2*K1(sqrts/T) / (4 pi (m)**2 *T K2(m/T))
c      so the S depedence drops
      sjac = sjac * dabs(2* ymin) * besselk(1, sqrts/T, 2) ! the 2 means that it is convoluted with exp(sqrts/T)
      sjac = sjac / (4*pi*m**2*T*besselk(2, m/T,1))        ! the 1 means no convolution 

      return 
      end

