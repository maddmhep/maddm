      subroutine find_resonances()

      implicit none
      include 'maddm.inc'
      include 'coupl.inc'

c parameters for the odeintegrator
      integer nres_above_threshold, ii

      include 'resonances.inc'

c convert the resonance locations into velocities
c and find the number of resonances above threshold
        nres_above_threshold = 0
        if (nres.gt.0) then
            do ii=1, nres
c               write(*,*) 'Found a resonance with mass: ', resonances(ii)
               if (2.d0*mdm(1).le.resonances(ii)) then
c                 write(*,*) 'Above threshold resonance found!'
                  nres_above_threshold = nres_above_threshold + 1
                  beta_res(ii) = dsqrt(1.d0 - 4.d0*mdm(1)**2 / (resonances(ii)**2))
                  beta_res_width(ii) = dsqrt(1.d0 - 4.d0*mdm(1)**2 / (resonances(ii)
     .                +resonance_widths(ii))**2) - beta_res(ii)
c                  write(*,*) beta_res(ii), beta_res_width(ii)
                else
                    beta_res(ii) = -1.d0
                    beta_res_width(ii) = -1.d0
                endif
            enddo
         endif

       end

        function relic_density(relic_can)
c-----------------------------------------------------------
c	This is the driver for the relic density calculation
c	relic_canon is a logical flag which determines whether
c	the calculation should proceed via coupled differential
c	equations, or via a summed approximation (as in case of
c	co-annihilations), also referred to as 'can'.
c   ndmparticles is the number of DM states and is relevant
c   only for relic_canonical=.false. scenario.
c   x_f is a global variable representing the freezeout x = m/T
c   sigmav_xf is also a global variable representing the vel. ave.
c   cross section at x_f.
c-----------------------------------------------------------
c         calculate the Wijs for all the annihilation channels
          implicit none
          include 'maddm.inc'
          include 'coupl.inc'

		  logical relic_can

		  include 'resonances.inc'

          call find_resonances()
		  call calculate_Wij_ann()

c         calculate the relic abundance of the DM candidate using either the canonical or coupled method
		  if (relic_can) then
			nvar = 1
			relic_density = relic_canon(x_f)
			sigmav_xf = taacs_canon(x_f)
c		    write(*,*) Oh2
c           In order to call odeint_check make sure dxsav, kmax, and kount are set in relicdensity.f
c           call odeint_check(0)

		  else
			nvar = ndmparticles
			relic_density =  relic_coupled()
			x_f = -1d0
			sigmav_xf = -1d0
c			write(*,*) Oh2

	      endif
         return
         end

