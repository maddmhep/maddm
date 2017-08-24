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
      double precision S, ve,v
      data ve /0d0/
      save ve
      double precision thres
      data thres /1e-7/

      double precision ymin      
      double precision besselk
c     parameter for newton method
      integer i
      double precision sqrts ! temporary value
      double precision f, df
c
c     need to integrate 
c     dv \sqrt{2/pi} v^3/ve^3 EXP(-v^2/ve^2) * sigma(ebeam1=ebeam2=\sqrt(M/(1-(v/2)^2))
c
c     The functions x^3 EXP(-x^2/a^2) has the following primitive:
c     -1/2 a^2 EXP[-x^2/a^2] * (a^2+x^2)
c
c     therefore we use the following change of variable (
c     y = -1/2 ve^2 EXP[-v^2/ve^2] * (ve^2+v^2)
c     which has the following derivative
c     dy = v^3 EXP[-v^2/ve^2] dv
c
c     v=0 -> ymin= -1/2 ve^4
c     v=1 -> ymax= 0  # technically not 0, but to avoid numerical accuracy setting it to zero and discard wrong events
c     
c     to pass to a number between 0 and 1, we then redefine
c      
cc     R = (y - ymin)/(-ymin) = -y/ymin +1
cc     dR = dy/(-ymin) = 2 v^3/ve^4 EXP[-v^2/ve^2] dv 
c
c     The integral to evaluate is then
c
c     \int_0^1 dR  v2/(\sqrt(2pi) * sigma(v(R))
c
c 
      m =pmass(1)
      if (ve.eq.0d0)then
         if (lpp(1).eq.9)then
            ve = ebeam(1)
         else
            ve = ebeam(2)
         endif
         if (pmass(1).ne.pmass(2))then
            stop 1
         endif
      endif
      ymin = -1/2*ve**4
c
c     Solve v for the following equation (newton method)
c     R = (y - ymin)/(-ymin) = (1/2 ve^2 EXP[-v^2/ve^2] * (ve^2+v^2)/ymin +1) = (-1/ve^2 EXP[-v^2/ve^2] * (ve^2+v^2) +1)
c
      v = ve ! initial guess
      do i=1,1000
         f = (-1/ve**2 * DEXP(-v**2/ve**2) * (ve**2+v**2) +1) - R
         if (abs(f).lt.thres) then
c            write(*,*) 'v is ', v, 'precision',f
            exit
         endif
         df = 2d0 * v**3/ve**4 *DEXP(-v**2/ve**2)
         v=v-f/df
      enddo
      if (abs(f).gt.thres.or.v.le.0d0.or.v.ge.1d0) then
         v = 0
         S = 0d0
         sjac = -1
         return
      endif
      cm_rap =0d0
      set_cm_rap=.true.
       pbeam1(0)=dsqrt( m**2 / (1-(v/2d0)**2))
       pbeam1(1) = 0d0
       pbeam1(2) = 0d0
       pbeam1(3)=sqrt(max(pbeam1(0)**2-m**2, 0d0))
       pbeam2(0)=pbeam1(0)
       pbeam2(1) = 0d0
       pbeam2(2) = 0d0
       pbeam2(3)=-sqrt(max(pbeam2(0)**2-m**2, 0d0))
       ebeam(1) = pbeam1(0)
       ebeam(2) = pbeam2(0)      
       sqrts = 2* pbeam1(0)
       stot=sqrts**2
       shat = stot
c       write(*,*) 'sqrts', sqrts
c
c      Setting the Jacobian AND the associate PDF
c
c      normally we should do:
c      sjac = sjac /( 2 v^3/ve^4 EXP[-v^2/ve^2])
c      But we will multiply here by the PDF!
c      \sqrt{2/pi} v^3/ve^3 EXP(-v^2/ve^2)      
c      so many term cancel out
       sjac = 1d0/DSQRT(2d0*pi)*ve    
       return 
      end

c      ==============================
c      OLD CODE FOR RELATIVISTIC CASE
c      ==============================
c 
c        has numerical issue with EXP(-4000)
c
cc     The integral to compute is 
cc     d\sqrt(S) \sart(s)^2 K1(\sqrt{s}/T)
cc     to generate events we consider the approximate
cc     d\sqrt(S) \sart(s)^2 exp[-\sqrt(S)/T]
cc     This integral is know analytically:
cc     -T* exp(-x/T) *(T^2 + (x+T)^2)
cc     
cc     Therefore we will apply the following change of variable:
cc     y = -T exp(-x/T) *(T^2 + (x+T)^2))
cc    The associated jacobian being:
cc     dy = x^2 exp[-x/T] dx
cc
cc     Therefore the integral to compute is 
cc
cc     dy K1(x(y)/T)*exp(x(y)/T)
cc
cc     the bound of integration are
cc     ymin = -T exp(-2m/T) (T^2 + (2m+T)^2)) 
cc     ymax = 0 
cc
cc
cc     In order to map this integral between 0 and 1, we then define
cc     the following change of variable
cc     R = (y - ymin)/(-ymin) = -y/ymin +1
cc     dR = dy/(-ymin) = x^2 exp[-x/T]/ymin dx
c      m =pmass(1)
c      if (lpp(1).eq.9)then
c         T = ebeam(1)
c      else
c         T = ebeam(2)
c      endif
c      if (pmass(1).ne.pmass(2))then
c         stop 1
c      endif
c      ymin = -T * DEXP(-1*(2*m)/T) *(T**2+(2*m+T)**2)
cc
cc     Solve x=sqrts for the following equation
cc     R = (y - ymin)/(-ymin) = (-T exp(-x/T) *(T^2 + (x+T)^2))-ymin)/(-ymin)
cc
c      sqrts = T
c      do i=1,1000
c         f = (-T*DEXP(-sqrts/T) *(T**2 + (sqrts+T)**2)-ymin)/(-ymin)-X1
c         if (f.lt.thres) then
c            exit
c         endif
c         df = sqrts**2 *DEXP(-sqrts/T)/(-ymin)
c         sqrts=sqrts-f/df
c      enddo
c      if (f.gt.thres) then
c         sqrts = 0d0
c         S = 0d0
c         sjac = -1d0
c         return
c      endif
cc
cc      setting some global variable/output variable
cc
c       stot =sqrts**2
c       shat = stot
c       cm_rap = 0d0
c       set_cm_rap=.true.
c       pbeam1(0)=dsqrt(stot)/2
c       pbeam1(3)=sqrt(max(pbeam1(0)**2-m**2, 0d0))
c       pbeam2(0)=pbeam1(0)
c       pbeam2(3)=-sqrt(max(pbeam2(0)**2-m**2, 0d0))
c       ebeam(1) = pbeam1(0)
c       ebeam(2) = pbeam2(0)
c       x1 = 1d0
cc
cc      Setting the Jacobian AND the associate PDF
cc
cc      normally we should do:
cc      sjac = sjac * dabs(ymin /S * DEXP[sqrts/T])   
cc      But we will multiply here by the PDF!
cc      2*sqrts**2*K1(sqrts/T) / (4 pi (m)**2 *T K2(m/T))
cc      so the S depedence drops
c      sjac = sjac * dabs(2* ymin) * besselk(1, sqrts/T, 2) ! the 2 means that it is convoluted with exp(sqrts/T)
c      sjac = sjac / (4*pi*m**2*T*besselk(2, m/T,1))        ! the 1 means no convolution 
c
c
c     cccccccc TEST WITH V^5 instead of V^3 # newton method fails
ccc
ccc     need to integrate 
ccc     dv \sqrt{2/pi} v^3/ve^3 EXP(-v^2/ve^2) * sigma(ebeam1=ebeam2=\sqrt(M/(1-(v/2)^2))
ccc     We assume that sigma is here proportionel to v^2
ccc
ccc     The functions x^5 EXP(-x^2/a^2) has the following primitive:
ccc     -1/2 a^2 EXP[-x^2/a^2] * (2a^4+2a^2x^2+x^4)
ccc
ccc     therefore we use the following change of variable (
ccc     y = -1/2 ve^2 EXP[-v^2/ve^2] * (2ve^4+ 2ve^2v^2+v^4)
ccc     which has the following derivative
ccc     dy = v^5 EXP[-v^2/ve^2] dv
ccc
ccc     v=0 -> ymin= -1 ve^6
ccc     v=1 -> ymax= 0  # technically not 0, but to avoid numerical accuracy setting it to zero and discard wrong events
ccc     
ccc     to pass to a number between 0 and 1, we then redefine
ccc      
cccc     R = (y - ymin)/(-ymin) = -y/ymin +1
cccc     dR = dy/(-ymin) = v^5/ve^6 EXP[-v^2/ve^2] dv
ccc
ccc     The integral to evaluate is then
ccc
ccc     \int_0^1 dR  1/(\sqrt(2pi)) * ve^3/v(R)^2 * sigma(v(R))
ccc
cc      m =pmass(1)
cc      if (ve.eq.0d0)then
cc         if (lpp(1).eq.9)then
cc            ve = ebeam(1)
cc         else
cc            ve = ebeam(2)
cc         endif
cc         if (pmass(1).ne.pmass(2))then
cc            stop 1
cc         endif
cc      endif
cc      ymin = -ve**6
ccc
ccc     Solve v for the following equation (newton method)
ccc     R = (y - ymin)/(-ymin) = (1/2 ve^2 EXP[-v^2/ve^2] * (ve^2+v^2)/ymin +1) = (-1/ve^2 EXP[-v^2/ve^2] * (ve^2+v^2) +1)
ccc
cc      v = ve ! initial guess
cc      do i=1,1000
cc         f = 1d0/2d0 * DEXP(-v**2/ve**2) * (2*ve**4+ 2*ve**2*v**2+v**4)/ve**4+1-R
cc         if (abs(f).lt.thres) then
ccc            write(*,*) 'v is ', v, 'precision',f
cc            exit
cc         endif
cc         df = v**5/ve**6 * DEXP(-v**2/ve**2)
cc         v=v-f/df
c      enddo
cc      if (abs(f).gt.thres.or.v.le.0d0.or.v.ge.1d0) then
c         v = 0
c         S = 0d0
c         sjac = -1
c         return
c      endif
c      cm_rap =0d0
c      set_cm_rap=.true.
c       pbeam1(0)=dsqrt( m**2 / (1-(v/2d0)**2))
c       pbeam1(3)=sqrt(max(pbeam1(0)**2-m**2, 0d0))
c       pbeam2(0)=pbeam1(0)
c       pbeam2(3)=-sqrt(max(pbeam2(0)**2-m**2, 0d0))
c       ebeam(1) = pbeam1(0)
c       ebeam(2) = pbeam2(0)      
c       sqrts = 2* pbeam1(0)
c       stot=sqrts**2
c       shat = stot
c       write(*,*) 'sqrts', sqrts
c
c      Setting the Jacobian AND the associate PDF
c      so many term cancel out
c       sjac =  1d0/(dsqrt(2d0*pi)) * ve**3/v**2    
c       return 
c       end
