c-------------------------------------------------------------------------c
      function relic_coupled()
c-------------------------------------------------------------------------c
c                                                                         c
c  This subroutine calculates the relic abundance given any model that    c
c  contains a dark matter candidate.  The algorithm basic is as follows.  c
c                                                                         c
c  - Set up the model parameters (i.e. reading in the param card)         c
c  - Calculate the Wij's with respect to s                                c
c  - Calculate the thermally averaged annihilation cross section wrt m/T  c
c  - integrate the chemical rate equation to find the relic abundance     c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      include 'maddm.inc'
      include 'coupl.inc'

c parameters for the odeintegrator
      integer nok, nbad
      double precision Y0(maxdms), Y1(maxdms), x_1, x_2, h1, hmin
      external derivative, rkqs

c parameters used in this subroutine only
      double precision g_star_S_check, Y_eq_check
      logical ODEs_solved

c calculate the Wij's for the DM -> DM and DM/SM scattering processes
      call calculate_Wij_other(2)
      call calculate_Wij_other(3)

c calculate the taacs for all process groups
      if (calc_taacs_ann_array) call calculate_taacs_array(1)
      if (calc_taacs_dm2dm_array) call calculate_taacs_array(2)
      if (calc_taacs_scattering_array) call calculate_taacs_array(3)

c output to show overall progress in relic density calculation
      if (print_out) write(*,*) 'Intergrating ODE for relic abundance'

c setting up the boundaries and set sizes for the ODE integrator
      nvar = ndmparticles
      x_1 = 31.d0
      x_2 = 1000.d0
      h1 = 0.1d0
      hmin = 0.0d0
      do x1=1, ndmparticles
        Y1(x1) = 1.d-20
      enddo
      ODEs_solved = .false.

c We start at x_i = 30 and step back the integration until the results of two concecutive
c calculations return a result that has a relative difference of epsilon for each ODE.
      do while (.not. solved_ODE(1))

c steps back the start of the ODE integration until we find the correct starting point below freezeout
        x_1 = x_1 - 1.d0

c Allows us to compare this calculation to the previous one
        do x1=1, ndmparticles
          Y0(x1) = Y1(x1)
        enddo

c get the number of relativistic degrees of freedom
        g_star_S = get_gstar(mdm(1)/x_1)
C         g_star_S = lineint(mdm(1)/x_1,gstarS_temperature,gstarS,100) 

c Initial thermal equilibrium of the number density per comoving volume
c Y_i = Y_eq = 45/(4*pi^4) g_i/g_*S x^2 K_2(x), with x=m_i/T
        do x1=1,ndmparticles
          Y1(x1) = get_Y_eq(mdm(x1),mdm(1)/x_1,dof_total(x1),g_star_S)
        enddo

c passese the temperature so that other functions can use it
        x_pass = x_1

C c if we wish to save the ODE integration and print it to a file or the screen we need to set
C c the step size to save and the maximum number of saved points
C         dxsav = 0.5d0
C         kmax = 200
C         kount = 0

c calls the ODE integrator to get the final number density per entropy density.
        call odeint(Y1,nvar,x_1,x_2,eps_ode,h1,hmin,nok,nbad,derivative,rkqs)

C c For testing purposes we can print out the saved steps of each iteration to the screen
C         call odeint_check(1)
C 
C c Also, for testing purposes we can print out the initial x_i and the resulting Y_inf
C         write(*,*) x_1, (Y1(x1), x1=1,ndmparticles)

C c print out the saved ODE information
C         do x1=1,ndmparticles
C           write(*,*) 'particle ',x1
C           do x2=1,kount
C             write(*,*) xp(x2),yp(x1,x2)
C           enddo
C         enddo

C c Using the saved information of the differential equation we see if the solution starts to match
C c thermal equilibrium so we can just use thermal equilibrium for the solution below this x value
C         g_star_S_check = get_gstar(mdm(1)/xp(2))
C         do x1=1, ndmparticles
C           if (.not. found_xf(x1)) then
C             Y_eq_check = get_Y_eq(mdm(x1),mdm(1)/xp(2),dof_total(x1),g_star_S_check)
C             if ((dabs(yp(x1,2)-Y_eq_check)/Y_eq_check).le.0.0001d0) then
C               xf(x1) = x_1-1
C               found_xf(x1) = .true.
C               if (print_out) write(*,*) 'Found xf for differential equation ', x1,' at x=',xf(x1)
C             endif
C           endif
C         enddo

c Checks to see if the results of this ODE calculation are within eps_ode of the previous calculation
c If so then set flag that the ODE was solved and the xf where the solution occured
        do x1=1, ndmparticles
          if (dabs((Y1(x1)-Y0(x1))/Y1(x1)) .le. eps_ode) then
            solved_ODE(x1) = .true.
          endif
        enddo

c Check to see if all of the differential equations have been solved
        ODEs_solved = .true.
        do x1=1, ndmparticles
          if (.not. solved_ODE(x1)) then
            ODEs_solved = .false.
          endif
        enddo

      enddo

c The relic density is summed over the individual results for each particle's differential equation
      relic_coupled = 0.d0
      do x1=1, ndmparticles
        relic_coupled = relic_coupled + Y1(x1)*2889.2d0*mdm(x1)/(1.05d-5)
      enddo

      if (print_out) write(*,*) 'Relic density from the density evolution equation ',relic_coupled

      return
      end



c-------------------------------------------------------------------------c
      subroutine calculate_Wij_other(group)
c-------------------------------------------------------------------------c
c                                                                         c
c  This routine calculates the Wij's for all of the DM -> DM processes.   c
c  Using the same betas that were used to calculate the normal Wij's to   c
c  SM particles and the process_dm2dm flag we get all the necessary info  c
c  to properly use these processes in the boltzman equations.             c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      include 'maddm.inc'
      include 'coupl.inc'

c input parameters
      integer group

c parameters used in this subroutine only
      integer group_pass
      common/group_pass/ group_pass

c external integrand for Wij integral
      external Wij_integrand_other

c set the group_pass variable to pass it into the integrand
      group_pass = group

c Calculate the Wij's for DM -> DM processes
      if (group.eq.2) then
        if (print_out) write(*,*) 'Calculating Wijs for DM -> DM processes'

c loop over all DM particle and process combinations
        do x1=1,ndmparticles
          do x2=x1,ndmparticles
            do x3=1,dm2dm_nprocesses(x1,x2)

c loop over all the betas and calculate the Wijs
              do x4=1,nWij

c setting variables to pass to Wij_integreand_other
                beta_pass = betas(x4)
                dmi = x1
                dmj = x2
                p_k = x3

c integrate the 2 -> 2 process over x = (cos(theta)+1)/2
                call romberg(Wij_integrand_other, 0.d0, 1.d0, Wij_dm2dm(x1,x2,x3,x4), eps_wij, iter_wij)

              enddo
            enddo
          enddo
        enddo

        if (print_out) write(*,*) 'Done!'
        Wij_dm2dm_calc = .true.

c Calculate the Wij's for the DM/SM scattering processes        
      else if (group.eq.3) then
        if (print_out) write(*,*) 'Calculating Wijs for DM/SM scattering processes'

c loop over all DM particle and process combinations
        do x1=1,ndmparticles
          do x2=1,ndmparticles
            do x3=1,scattering_nprocesses(x1,x2)

c loop over all the betas and calculate the Wijs
              do x4=1,nWij

c setting variables to pass to Wij_integrand_other
                beta_pass = betas(x4)
                dmi = x1
                dmj = x2
                p_k = x3

c integrate the 2 -> 2 process over x = (cos(theta)+1)/2
                call romberg(Wij_integrand_other, 0.d0, 1.d0, Wij_scattering(x1,x2,x3,x4), eps_wij, iter_wij)

              enddo
            enddo
          enddo
        enddo
        
        if (print_out) write(*,*) 'Done!'
      endif

      return
      end



c-------------------------------------------------------------------------c
      function Wij_integrand_other(x)
c-------------------------------------------------------------------------c
c                                                                         c
c  This function is the integrand called by the get_Wij function to       c
c  calculate the matrix elements squared integrated over the final state  c
c  phase space.                                                           c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      include 'maddm.inc'
      include 'coupl.inc'

c input parameters
      double precision x

c parameters used in this subroutine only
      double precision beta, s, msq, p_init(0:3), p_ext(0:3,4), ps_point(2), wgt

c group pass common block
      integer group_pass
      common/group_pass/ group_pass

c get the final state masses using the pmass(i)_(j)_(k).inc files
      if (group_pass.eq.2) then
        call getfsmasses_dm2dm(pmass,dmi,dmj,p_k)
      else if (group_pass.eq.3) then
        call getfsmasses_scattering(pmass,dmi,dmj,p_k)
      endif

c calculating the center of mass energy squared given the passed variable beta_pass
      beta = beta_pass
      s = (pmass(1)+pmass(2))**2/(1.d0-beta**2)

c ps_point determines cos(theta) and phi given two random numbers from 0 to 1
c 2 to 2 processes are azimuthally symmetric so we set phi=0
      ps_point(1) = x
      ps_point(2) = 0.d0

c center of mass 4-momentum
      p_init(0) = dsqrt(s)
      p_init(1) = 0.d0
      p_init(2) = 0.d0
      p_init(3) = 0.d0

c setting up the external particle momenta for the matrix element      
      p_ext(0,1) = (s+pmass(1)**2-pmass(2)**2)/(2.d0*dsqrt(s))
      p_ext(0,2) = (s+pmass(2)**2-pmass(1)**2)/(2.d0*dsqrt(s))
      p_ext(1,1) = 0.d0
      p_ext(1,2) = 0.d0
      p_ext(2,1) = 0.d0
      p_ext(2,2) = 0.d0
      p_ext(3,1) = dsqrt(lambda(s,pmass(1)**2,pmass(2)**2))/(2.d0*dsqrt(s))
      p_ext(3,2) = -dsqrt(lambda(s,pmass(1)**2,pmass(2)**2))/(2.d0*dsqrt(s))

c initialize the integrand
      Wij_integrand_other = 0.d0

c check to see if the center of mass energy is enough to produce the final states
      if (s.ge.(pmass(3)+pmass(4))**2) then

c two body phase space subroutine
        call decay2(p_init,pmass(3),pmass(4),pf1,pf2,wgt,ps_point,1)

        do x1=0,3
          p_ext(x1,3) = pf1(x1)
          p_ext(x1,4) = pf2(x1)
        enddo

c calculate the matrix element      
        if (group_pass.eq.2) then
          msq = smatrix_dm2dm(p_ext,dmi,dmj,p_k)
          if (dm2dm_process_iden_init(dmi,dmj,p_k)) then
            msq = msq/2.d0
          endif
        else if (group_pass.eq.3) then
          msq = smatrix_scattering(p_ext,dmi,dmj,p_k)
        endif

c MadGraph automatically averages over the initial state degrees of freedom.  Our definition
c of Wij requires inital states to be summed, hence factors of dof_dm and dof_SM
        if (group_pass.eq.2) then
          Wij_integrand_other = dof_dm(dmi)*dof_dm(dmj)*msq*wgt
        else if (group_pass.eq.3) then
          Wij_integrand_other = dof_dm(dmi)*dof_SM(dmi,dmj,p_k)*msq*wgt
        endif

      endif

      return
      end



c-------------------------------------------------------------------------c
      subroutine calculate_taacs_array(group)
c-------------------------------------------------------------------------c
c                                                                         c
c  This subroutine calculates the thermally averaged annihilation cross   c
c  for all the DM -> DM processes over a range of x.  When integrating    c
c  the ODEs these values are then interpolated.                           c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      include 'maddm.inc'
      include 'coupl.inc'

c input parameter
      integer group

c Calculate the taacs array for the annihilation processes
      if (group.eq.1) then
        if (print_out) write(*,*) 'Calculating the taacs arrays for annihilation processes'

c loop over all the combinations of initial state DM particles
        do x1=1,ndmparticles
          do x2=x1,ndmparticles

c loop over all the values of x_taacs
            do x3=1,100
              taacs_ann(x1,x2,x3) = taacs_coupled(x_taacs(x3),x1,x2,0,group)
              if (x1.ne.x2) taacs_ann(x2,x1,x3) = taacs_ann(x1,x2,x3)
            enddo
          enddo
        enddo

        if (print_out) write(*,*) 'Done!'
        taacs_ann_calc = .true.

c Calculate the taacs array for the DM -> DM processes
      else if (group.eq.2) then
        if (print_out) write(*,*) 'Calculating the taacs arrays for DM -> DM processes'

c loop over all the combinations of initial state DM particles and processes
        do x1=1,ndmparticles
          do x2=x1,ndmparticles
            do x3=1,dm2dm_nprocesses(x1,x2)

c loop over all the values of x_taacs
              do x4=1, 100
                taacs_dm2dm(x1,x2,x3,x4) = taacs_coupled(x_taacs(x4),x1,x2,x3,group)
              enddo
            enddo
          enddo
        enddo
        
        if (print_out) write(*,*) 'Done!'
        taacs_dm2dm_calc = .true.

c Calculate the taacs array for the DM/SM scattering processes
      else if (group.eq.3) then
        if (print_out) write(*,*) 'Calculating the taacs arrays for the DM/SM scattering processes'
        
c loop over all the combinations of initial state DM particles
        do x1=1,ndmparticles
          do x2=1,ndmparticles
            do x3=1,scattering_nprocesses(x1,x2)

c loop over all the values of x_taacs
              do x4=1,100
                taacs_scattering(x1,x2,x3,x4) = taacs_coupled(x_taacs(x4),x1,x2,x3,group)
              enddo
            enddo
          enddo
        enddo
        
        if (print_out) write(*,*) 'Done!'
        taacs_scattering_calc = .true.
        
      endif

      return
      end



c-------------------------------------------------------------------------c
      function taacs_coupled(x,i,j,k,group)
c-------------------------------------------------------------------------c
c                                                                         c
c  This function calculates the thermally averaged annihilation cross     c
c  section. to SM particles given temperature x and initial states        c
c  i and j.                                                               c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      include 'maddm.inc'

c input parameters
      double precision x
      integer i, j, k, group

c holding variable just in case we're at some point x that can't be integrated easily
      double precision taacs_ann_hold(maxdms,maxdms), taacs_dm2dm_hold(maxdms,maxdms,max_dm2dm)
      double precision taacs_scattering_hold(maxdms,maxdms,max_dm2dm)
      common/taacs_hold/ taacs_ann_hold, taacs_dm2dm_hold, taacs_scattering_hold

c group passing variable so that the integrand knows what group of processes we are calculating
      integer group_pass
      common/group_pass/ group_pass

c external function to be integrated over
      external taacs_integrand

c setting the input variables to be passed to taacs_integrand
      dmi = i
      dmj = j
      p_k = k
      x_pass = x
      group_pass = group

c This is used if we are doing a 1-d integration of taacs
      call romberg(taacs_integrand, 0.d0, 1.d0, taacs_coupled, eps_taacs, iter_taacs)

c if the integration worked (i.e. didn't return a nan) then set the taacs variable
      if (taacs_coupled.eq.taacs_coupled) then
        if (group.eq.1) then
          taacs_ann_hold(i,j) = taacs_coupled
        else if (group.eq.2) then
          taacs_dm2dm_hold(i,j,k) = taacs_coupled
        else if (group.eq.3) then
          taacs_scattering_hold(i,j,k) = taacs_coupled
        endif
      else
        if (group.eq.1) then
          taacs_coupled = taacs_ann_hold(i,j)
        else if (group.eq.2) then
          taacs_coupled = taacs_dm2dm_hold(i,j,k)
        else if (group.eq.3) then
          taacs_coupled = taacs_scattering_hold(i,j,k)
        endif
      endif

C c Return a constant value for taacs (for testing purposes ONLY)
C       taacs = 1.0d-10

      return
      end



c-------------------------------------------------------------------------c
      function taacs_integrand(beta)
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
      double precision T, s, Wij_value, K2_prod

c group passing variable so that the integrand knows what group of processes we are calculating
      integer group_pass
      common/group_pass/ group_pass

c If beta is zero or if we hit any beta = 1 singularities then automatically return zero
      taacs_integrand = 0.d0
      if (((1.d0-beta**2).eq.0.d0).or.(beta.eq.0)) return

c setting the temperature and the center of mass energy squared
      T = mdm(1)/x_pass

c taacs_integrand for the annihilation processes
      if ((group_pass.eq.1).or.(group_pass.eq.2)) then

c the K2_product appears in the denominator of taacs
c K2_prod = \prod(i=dmi,dmj) g_{i} m_{i}^{2} K_{2}(m_{i}/T)/K_{2}(m_{1}/T)
c the bessel functions are expanded for large values of m/T

c dmi
        K2_prod = dof_total(dmi)*mdm(dmi)**1.5d0*dsqrt(mdm(1))*dexp((mdm(1)-mdm(dmi))/T)
     .        * (1.d0 + 15.d0*T/(8.d0*mdm(dmi)) + 105.d0*T**2/(128.d0*mdm(dmi)**2) - 315.d0*T**3/(1024.d0*mdm(dmi)**3))
     .        / (1.d0 + 15.d0*T/(8.d0*mdm(1)) + 105.d0*T**2/(128.d0*mdm(1)**2) - 315.d0*T**3/(1024.d0*mdm(1)**3))

c dmj
        K2_prod = K2_prod * dof_total(dmj)*mdm(dmj)**1.5d0*dsqrt(mdm(1))*dexp((mdm(1)-mdm(dmj))/T)
     .        * (1.d0 + 15.d0*T/(8.d0*mdm(dmj)) + 105.d0*T**2/(128.d0*mdm(dmj)**2) - 315.d0*T**3/(1024.d0*mdm(dmj)**3))
     .        / (1.d0 + 15.d0*T/(8.d0*mdm(1)) + 105.d0*T**2/(128.d0*mdm(1)**2) - 315.d0*T**3/(1024.d0*mdm(1)**3))

        s = (mdm(dmi)+mdm(dmj))**2/(1.d0-beta**2)

c We interpolate the Wij value given beta
        Wij_value = get_Wij_value(beta,dmi,dmj,p_k,group_pass)

c integrand for the thermally averaged annihilation cross section
        taacs_integrand = (mdm(dmi)+mdm(dmj))**2*beta/(1.d0-beta**2)**2 * Wij_value / (8.d0*T*K2_prod)
     .      * dsqrt(lambda(s,mdm(dmi)**2,mdm(dmj)**2)/s)
     .      * dsqrt(2.d0/(pi*dsqrt(s)*T))*mdm(1)*dexp(2.d0*x_pass-dsqrt(s)/T)
     .      * (1.d0 + 3.d0*T/(8.d0*dsqrt(s)) - 15.d0*T**2/(128.d0*s) + 105.d0*T**3/(1024.d0*s**1.5d0))
     .      / (1.d0 + 15.d0*T/(8.d0*mdm(1)) + 105.d0*T**2/(128.d0*mdm(1)**2) - 315.d0*T**3/(1024.d0*mdm(1)**3))**2
        
      else if (group_pass.eq.3) then

c get the external masses to properly calculate the center of mass energy
        call getfsmasses_scattering(pmass, dmi, dmj, p_k)

c the K2_product appears in the denominator of taacs
c K2_prod = \prod(i=dmi,dmj) g_{i} m_{i}^{2} K_{2}(m_{i}/T)/K_{2}(m_{1}/T)
c the bessel functions are expanded for large values of m/T

c dmi (pmass(1))
        K2_prod = dof_total(dmi)*mdm(dmi)**1.5d0*dsqrt(mdm(1))*dexp((mdm(1)-mdm(dmi))/T)
     .        * (1.d0 + 15.d0*T/(8.d0*mdm(dmi)) + 105.d0*T**2/(128.d0*mdm(dmi)**2) - 315.d0*T**3/(1024.d0*mdm(dmi)**3))
     .        / (1.d0 + 15.d0*T/(8.d0*mdm(1)) + 105.d0*T**2/(128.d0*mdm(1)**2) - 315.d0*T**3/(1024.d0*mdm(1)**3))

c pmass(2)
c if pmass(2) is much larger than the temperature then we can use the same trick of normalization as the DM iniital states
        if (pmass(2).gt.5.d0*T) then
          K2_prod = K2_prod * dof_SM_total(dmi,dmj,p_k)*pmass(2)**1.5d0*dsqrt(mdm(1))*dexp((mdm(1)-pmass(2))/T)
     .          * (1.d0 + 15.d0*T/(8.d0*pmass(2)) + 105.d0*T**2/(128.d0*pmass(2)**2) - 315.d0*T**3/(1024.d0*pmass(2)**3))
     .          / (1.d0 + 15.d0*T/(8.d0*mdm(1)) + 105.d0*T**2/(128.d0*mdm(1)**2) - 315.d0*T**3/(1024.d0*mdm(1)**3))

c if pmass(2) is around the order of the temperature then we numerically calculate the bessel function
        else if ((pmass(2).le.5.d0*T).and.(pmass(2).gt.T/5.d0)) then
          K2_prod = K2_prod * dof_SM_total(dmi,dmj,p_k)*pmass(2)**2*bessk(2,pmass(2)/T)

c if pmass(2) is zero (or if its smaller than the temperature) then we have a nice analytic expression
        else if (pmass(2).le.T/5.d0) then
          K2_prod = K2_prod * 2.d0*dof_SM_total(dmi,dmj,p_k)*T**2
        endif

        s = (pmass(1)+pmass(2))**2/(1.d0-beta**2)

c check to see if the center of mass energy is large enough to create the final state particles
        if (s.gt.(pmass(3)+pmass(4))**2) then

c We interpolate the Wij value given beta
          Wij_value = get_Wij_value(beta,dmi,dmj,p_k,group_pass)
          
c integrand for the thermally averaged annihilation cross section
          taacs_integrand = (pmass(1)+pmass(2))**2*beta/(1.d0-beta**2)**2 * Wij_value / (8.d0*T*K2_prod)
     .        * dsqrt(lambda(s,pmass(1)**2,pmass(2)**2)/s)

c depending on the mass of the SM scattering particle we may or may not have factored in an K_2(m_{1}/T)
          if (pmass(2).gt.5.d0*T) then
            taacs_integrand = taacs_integrand * dsqrt(2.d0/(pi*dsqrt(s)*T))*mdm(1)*dexp(2.d0*x_pass-dsqrt(s)/T)
     .            * (1.d0 + 3.d0*T/(8.d0*dsqrt(s)) - 15.d0*T**2/(128.d0*s) + 105.d0*T**3/(1024.d0*s**1.5d0))
     .            / (1.d0 + 15.d0*T/(8.d0*mdm(1)) + 105.d0*T**2/(128.d0*mdm(1)**2) - 315.d0*T**3/(1024.d0*mdm(1)**3))**2
          else
            taacs_integrand = taacs_integrand * dsqrt(mdm(1)/dsqrt(s))*dexp(x_pass-dsqrt(s)/T)
     .            * (1.d0 + 3.d0*T/(8.d0*dsqrt(s)) - 15.d0*T**2/(128.d0*s) + 105.d0*T**3/(1024.d0*s**1.5d0))
     .            / (1.d0 + 15.d0*T/(8.d0*mdm(1)) + 105.d0*T**2/(128.d0*mdm(1)**2) - 315.d0*T**3/(1024.d0*mdm(1)**3))
          endif
        endif
      endif

      return
      end



c-------------------------------------------------------------------------c
      subroutine derivative(x_i,Y_i,dYdx)
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
      double precision x_i, Y_i(maxdms), dYdx(maxdms)

c parameters in this subroutine only
      double precision T, Y_eq(maxdms), taacs_value, rate_factor, Y_eq_SM
      double precision taacs_value2, Y_eq_SM2
      double precision dYdx_ann, dYdx_dm2dm, dYdx_scattering

c finding the number of relativistic degrees of freedom given temperature x_i
      T = mdm(1)/x_i
      g_star_S = get_gstar(T)
      g_star = g_star_S

C c old g_star and g_star_S interpolation
C c for testing purposes only
C       g_star = lineint(T,gstar_temperature,gstar,100)
C       g_star_S = lineint(T,gstarS_temperature,gstarS,100)

c Thermal equilibrium of the number density per comoving volume
c Y_eq = 45/(4*pi^4) g_1/g_*S x^2 K_2(x), with x=m/T
      do x1=1,ndmparticles
        Y_eq(x1) = get_Y_eq(mdm(x1),T,dof_total(x1),g_star_S)
      enddo

      do x1=1, ndmparticles

c initialize all the parts of the derivative
        dYdx_ann = 0.d0
        dYdx_dm2dm = 0.d0
        dYdx_scattering = 0.d0

c loop over all the possible particles that annihilate with x1
        do x2=1, ndmparticles

c thermally averaged annihilation cross section x1, x2 -> SM particles
          if (calc_taacs_ann_array) then
            taacs_value = get_taacs_value(x_i,x1,x2,0,1)
          else
            taacs_value = taacs_coupled(x_i,x1,x2,0,1)
          endif

c add the contribution of the annihilation channel
c If x1 = x2 then you have a factor of two to properly count the interaction rate
          if (x1.eq.x2) then
            dYdx_ann = dYdx_ann - dsqrt(pi/(45.d0*g_star))*g_star_S*mpl*mdm(1)/x_i**2
     .          *(2.d0*taacs_value)*(Y_i(x1)*Y_i(x2)-Y_eq(x1)*Y_eq(x2))
          else
            dYdx_ann = dYdx_ann - dsqrt(pi/(45.d0*g_star))*g_star_S*mpl*mdm(1)/x_i**2
     .          *(taacs_value)*(Y_i(x1)*Y_i(x2)-Y_eq(x1)*Y_eq(x2))
          endif
        enddo


c loop over all the DM -> DM processes and see if they need to be included
        do x2=1,ndmparticles
          do x3=x2,ndmparticles
            do x4=1,dm2dm_nprocesses(x2,x3)
              
c If the DM particle doesn't appear in the initial or final state then it doesn't contribute
              if ((x1.eq.x2).or.(x1.eq.x3).or.
     .              (x1.eq.dm2dm_fs(x2,x3,x4,1)).or.(x1.eq.dm2dm_fs(x2,x3,x4,2))) then

c Properly count the numer of interaction rates
                rate_factor = 0.d0
                if (x1.eq.x2) rate_factor = rate_factor - 1.d0
                if (x1.eq.x3) rate_factor = rate_factor - 1.d0
                if (x1.eq.dm2dm_fs(x2,x3,x4,1)) rate_factor = rate_factor + 1.d0
                if (x1.eq.dm2dm_fs(x2,x3,x4,2)) rate_factor = rate_factor + 1.d0

                if (calc_taacs_dm2dm_array) then
                  taacs_value = get_taacs_value(x_i,x2,x3,x4,2)
                else 
                  taacs_value = taacs_coupled(x_i,x2,x3,x4,2)
                endif

c the rate_factor contains the relative minus sign so we simply add the contribution
                dYdx_dm2dm = dYdx_dm2dm + dsqrt(pi/(45.d0*g_star))*g_star_S*mpl*mdm(1)/x_i**2
     .                  *(rate_factor*taacs_value)*Y_i(x2)*Y_i(x3)
              endif
            enddo
          enddo
        enddo


c loop over all the DM/SM scattering processes and see if they need to be included
        do x2=1,ndmparticles
          if (x1.eq.x2) then
            do x3=1,ndmparticles

c processes that reduce the number of x1 particles
              do x4=1,scattering_nprocesses(x2,x3)

                if (calc_taacs_scattering_array) then
                  taacs_value = get_taacs_value(x_i,x2,x3,x4,3)
                  taacs_value2 = get_taacs_value(x_i,x3,x2,x4,3)
                else 
                  taacs_value = taacs_coupled(x_i,x2,x3,x4,3)
                  taacs_value2 = taacs_coupled(x_i,x3,x2,x4,3)
                endif

                call getfsmasses_scattering(pmass,x2,x3,x4)
                Y_eq_SM = get_Y_eq(pmass(2),T,dof_SM_total(x2,x3,x4),g_star_S)
                call getfsmasses_scattering(pmass,x3,x2,x4)
                Y_eq_SM2 = get_Y_eq(pmass(2),T,dof_SM_total(x3,x2,x4),g_star_S)

c the rate_factor contains the relative minus sign so we simply add the contribution
                dYdx_scattering = -dsqrt(pi/(45.d0*g_star))*g_star_S*mpl*mdm(1)/x_i**2
     .                  *(taacs_value*Y_i(x2)*Y_eq_SM-taacs_value2*Y_i(x3)*Y_eq_SM2)
              enddo
            enddo
          endif
        enddo

c if there are scattering diagrams then it should display a warning to use the canonical method
        if ((dYdx_scattering.ne.0.d0).and.(.not.dm_sm_scattering_warn)) then
          write(*,*) 'There are non-zero DM/SM scattering diagrams!!'
          write(*,*) 'Use the canonical relic abundance calculation'
          write(*,*) 'Set relic_canonical to true in include/maddm_card.inc'
          dm_sm_scattering_warn = .true.
        endif
        
c Add up all the contributions to the derivative for x1
        dYdx(x1) = dYdx_ann + dYdx_dm2dm + dYdx_scattering

      enddo

      return
      end
