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

		  logical relic_can

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

