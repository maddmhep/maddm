c-------------------------------------------------------------------------c
      function relic_canon(xd_start)
c-------------------------------------------------------------------------c
c                                                                         c
c  This is the canonical calculation of the relic abundance of DM with    c
c  coannihilation particles.  If there are no coupled DM density          c
c  evolution equations then this should give the same result as           c
c  relic_canon.                                                           c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      include 'maddm.inc'
      include 'coupl.inc'

c parameters for the odeintegrator
      integer nok, nbad, nres_above_threshold, ii
      double precision Y0, Y1, x_1, x_2, h1, hmin, min_mass, Bfact, x_1_approx, xd_start
      external derivative_canon, rkqs

      include 'resonances.inc'

      if (xd_start.le.1d0) then
			    write(*,*) 'Freezeout occurs too early! Relic density too big to comprehend. Seeting Omegah^2 = -1.0 \n'
			    relic_canon=-1.0
			    return
      endif

c setting up the boundaries and set sizes for the ODE integrator
      nvar = 1    !may not necessarily be true in the future
      x_1 = xd_start
      x_2 = x_end
      h1 = 0.1d0
      hmin = 0.0d0
      Y0 = 0.d0
      Y1 = 1.d-20

c 	      We start at x_i = x_start and step back the integration until the results of two concecutive
c    	  claculations return a result that has a relative difference of epsilon.
		  do while (dabs((Y1-Y0)/Y1).gt.eps_ode)

c   	  	  Allows us to compare this iteration to the previous one
		  	  Y0 = Y1


c    		  Initial value of the relativistic degrees of freedom and number density per entropy density, Y
			  g_star_S = get_gstar(mdm(1)/x_1)
C             g_star_S = lineint(mdm(1)/x_1,gstarS_temperature,gstarS,100)

c    		  Initial thermal equilibrium of the number density per comoving volume
c   		  Y_1 = Y_eq = 45/(4*pi^4) g_1/g_*S x^2 K_2(x), with x=m/T
			  Y1 = 0.d0
			  do x1=1,ndmparticles
			      Y1 = Y1 + get_Y_eq(mdm(x1),mdm(1)/x_1,dof_total(x1),g_star_S)
			  enddo

c             passes the temperature so that other functions can use it
			  x_pass = x_1

C c           if we wish to save the ODE integration and print it to a file or the screen we need to set
C c           the step size to save and the maximum number of saved points
C             dxsav = 0.5d0
C             kmax = 200
C             kount = 0

c             calls the ODE integrator to get the final number density per entropy density.
			  call odeint(Y1,nvar,x_1,x_2,eps_ode,h1,hmin,nok,nbad,derivative_canon,rkqs)

C c         For testing purposes we can print out the saved steps of each iteration to the screen
C           call odeint_check(0)
C
C c         Also, for testing purposes we can print out the initial x_i and the resulting Y_inf
c            write(*,*) x_1, Y1

c             steps back the start of the ODE integration until we find the correct starting point below freezeout

c           If using the freezeout approximation
            if (xd_approx) then
                if (xd_start.le.10.d0) then
                    xd_start =  10.d0
                endif
            endif

			x_1 = x_1 - dx_step
           if (x_1.le.1d0) then
			    write(*,*) 'Freezeout occurs too early! Relic density too big to comprehend. Seeting Omegah^2 = -1.0 \n'
			    relic_canon=-1.0
			    return
			endif

		  enddo


c find the minimum DM particle mass (to calculate the overall relic density)
      min_mass = 1.d20
      do x1=1,ndmparticles
        min_mass = min(min_mass,mdm(x1))
      enddo

c finish the calculation of omega h^2
      relic_canon = Y1*2889.2d0*min_mass/(1.05d-5)
      x_f = x_1
c      if (print_out) write(*,*) 'Relic density from the canonical density evolution equation ',relic_canon

      return
      end



c-------------------------------------------------------------------------c
      subroutine derivative_canon(x_i,Y_i,dYdx)
c-------------------------------------------------------------------------c
c                                                                         c
c  This subroutine defines the derivative of Y with respect to x at any   c
c  point x along the ODE integrator.                                      c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      include 'maddm.inc'
      include 'coupl.inc'

c input parameters
      double precision x_i, Y_i, dYdx

c parameters in this subroutine only
      double precision s, Y_eq, taacs_value

      include 'resonances.inc'

c finding the number of relativistic degrees of freedom given temperature x_i
      g_star_S = get_gstar(mdm(1)/x_i)
      g_star = g_star_S

C c old g_star and g_star_S interpolation
C c for testing purposes only
C       g_star = lineint(T,gstar_temperature,gstar,100)
C       g_star_S = lineint(T,gstarS_temperature,gstarS,100)

c Thermal equilibrium of the number density per comoving volume
c Y_eq = 45/(4*pi^4) g_1/g_*S x^2 K_2(x), with x=m/T
      Y_eq = 0.d0
      do x1=1,ndmparticles
        Y_eq = Y_eq + get_Y_eq(mdm(x1),mdm(1)/x_i,dof_total(x1),g_star_S)
      enddo

c thermally averaged annihilation cross section
      taacs_value = taacs_canon(x_i)

c putting everything together, including the factor of two which is necessary to
c properly count the number of interactions
      dYdx = -dsqrt(pi/(45.d0*g_star))*g_star_S*mpl*mdm(1)/x_i**2
     .      *(2.d0*taacs_value)*(Y_i-Y_eq)*(Y_i+Y_eq)

C c this is for testing purposes only
C       write(*,fmt='(4(ES14.8,1x))') x_i,Y_i,Y_eq,dYdx

      return
      end

      double precision function cosine(x)
         implicit none
         double precision x, cos
         cosine = cos(x)
         return
      end

c-------------------------------------------------------------------------c
      function taacs_canon(x)
c-------------------------------------------------------------------------c
c                                                                         c
c  This function calculates the thermally averaged annihilation cross     c
c  section.                                                               c
c                                                                         c
c-------------------------------------------------------------------------c
      include 'maddm.inc'
      include 'coupl.inc'
c input parameters
      double precision x,  cosine

c external function to be integrated over
      external taacs_integrand_canon, cosine

      include 'resonances.inc'


c passing x so that taacs_integrand can use it
      x_pass = x

c This is used if we are doing a 1-d integration of taacs
c      call romberg(taacs_integrand_canon, 0.d0, 1.d0, taacs_canon, eps_taacs, iter_taacs)

c call simpson's rule (assumes grid is set up)
        taacs_canon =simpson(taacs_integrand_canon, 0.d0, 1.d0, grid, grid_npts)
c       write(*,*) 'test: ', testvar

C c Return a constant value for taacs (for testing purposes ONLY)
C       taacs_canon = 3.0d-9

      return
      end


c-------------------------------------------------------------------------c
      function taacs_integrand_canon(beta)
c-------------------------------------------------------------------------c
c                                                                         c
c  This is the integrand of the thermally averaged annihilation cross     c
c  section.  It is integrated over the velocity parameter beta.           c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      include 'maddm.inc'
      include 'coupl.inc'

c input parameters
      double precision beta

c parameters used in this function only
      double precision T, s, Wij_value, Wij_sum, K2_sum


      include 'resonances.inc'

c making sure that we don't hit any beta = 1 singularities
      if ((1.d0-beta**2).eq.0.d0) then
        taacs_integrand_canon = 0.d0
        return
      endif

c setting the temperature and the center of mass energy squared
      T = mdm(1)/x_pass
      s = 4.d0*mdm(1)**2/(1.d0-beta**2)

c the K2_sum appears in the denominator of taacs
c \sum_{i}^{n} g_{i} m_{i}^{2} K_{2}(m_{i}/T)/K_{2}(m_{1}/T)
c the bessel functions are expanded for large values of m/T
      K2_sum = 0.d0
      do x1=1, ndmparticles
        K2_sum = K2_sum + dof_total(x1)*mdm(x1)**1.5d0*dsqrt(mdm(1))*dexp((mdm(1)-mdm(x1))/T)
     .      *(1.d0 + 15.d0*T/(8.d0*mdm(x1)) + 105.d0*T**2/(128.d0*mdm(x1)**2) - 315.d0*T**3/(1024.d0*mdm(x1)**3))
     .      /(1.d0 + 15.d0*T/(8.d0*mdm(1)) + 105.d0*T**2/(128.d0*mdm(1)**2) - 315.d0*T**3/(1024.d0*mdm(1)**3))
      enddo

c To get the effective thermally averaged annihilation cross section the individual Wij's
c need to be weighted by a momentum factor for each \chi_{i}, \chi_{j} pair
      Wij_sum = 0.d0
      do x1=1, ndmparticles
        do x2=x1, ndmparticles
          dmi = x1
          dmj = x2
          if (s.ge.(mdm(x1)+mdm(x2))**2) then
            Wij_value = get_Wij_value(beta,x1,x2,0,1)
            Wij_sum = Wij_sum + max(0.d0,dsqrt(lambda(s,mdm(x1)**2,mdm(x2)**2)/s)*Wij_value)
          endif
        enddo
      enddo

c integrand for the thermally averaged annihilation cross section
      if (2.d0*x_pass*(1.d0-1.d0/(1.d0-beta**2)**(0.5d0)).lt.-695)then
         taacs_integrand_canon =mdm(1)**2*beta/(1.d0-beta**2)**2 * Wij_sum / (2.d0*T*K2_sum**2)
     .    * 0d0
     .    * (1.d0 + 3.d0*T/(8.d0*dsqrt(s)) - 15.d0*T**2/(128.d0*s) + 105.d0*T**3/(1024.d0*s**1.5d0))
     .    / (1.d0 + 15.d0*T/(8.d0*mdm(1)) + 105.d0*T**2/(128.d0*mdm(1)**2) - 315.d0*T**3/(1024.d0*mdm(1)**3))**2
      else
      taacs_integrand_canon = mdm(1)**2*beta/(1.d0-beta**2)**2 * Wij_sum / (2.d0*T*K2_sum**2)
     .    * dsqrt(2.d0/(pi*dsqrt(s)*T))*mdm(1)*dexp(2.d0*x_pass*(1.d0-1.d0/(1.d0-beta**2)**(0.5d0)))
     .    * (1.d0 + 3.d0*T/(8.d0*dsqrt(s)) - 15.d0*T**2/(128.d0*s) + 105.d0*T**3/(1024.d0*s**1.5d0))
     .    / (1.d0 + 15.d0*T/(8.d0*mdm(1)) + 105.d0*T**2/(128.d0*mdm(1)**2) - 315.d0*T**3/(1024.d0*mdm(1)**3))**2
      endif
      return
      end
