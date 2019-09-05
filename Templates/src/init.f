c-------------------------------------------------------------------------c
      subroutine init_relic()
c-------------------------------------------------------------------------c
c                                                                         c
c  This subroutine initializes some of the parameters needed to do the    c
c  relic abundance calculation.                                           c
c                                                                         c
c  Initialized variables:                                                 c
c  - the masses and number of degrees of freedom for the DM particles     c
c  - arrays used to calculate the total number of relativistic degrees    c
c       of freedom                                                        c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      include 'maddm.inc'
      include 'coupl.inc'

c parameters used in this subroutine only
      character*30 filename

      if (print_out) write(*,*) 'Initializing numerical session'

c Setting flags to tell if certain numerical quantities have been calculated
      Wij_ann_calc = .false.
      Wij_dm2dm_calc = .false.
      Wij_scattering_calc = .false.
      taacs_ann_calc = .false.
      taacs_dm2dm_calc = .false.
      taacs_scattering_calc = .false.
      dm_sm_scattering_warn = .false.

c Reading in the model info
      filename = 'include/model_info.txt'
      open(unit=42,file=filename)
      read(42,*) nparticles
      do x1=1, nparticles
        read(42,*) masses(x1), dofs(x1), isfermion(x1)
      enddo
      close(42)

c Setting up the DM masses and degrees of freedom
      include 'dm_info.inc'
      do x1=1,ndmparticles
        dof_total(x1) = dof_dm(x1)
        if (dm_names(x1).ne.dm_antinames(x1)) then
          dof_total(x1) = 2.d0*dof_dm(x1)
        endif
      enddo

      include 'process_names.inc'

c initialize the arrays for the relativistic degrees of freedom
      filename = 'include/fermion.txt'
      call readinfile(filename,401,gstar_x,gstar_fermion)
      filename = 'include/boson.txt'
      call readinfile(filename,401,gstar_x,gstar_boson)

C c old initialization of digitized g_star and g_star_S curves
C c for testing purposes only!
C       filename = 'include/g_star.txt'
C       call readinfile(filename,100,gstar_temperature,gstar)
C       filename = 'include/g_star_S.txt'
C       call readinfile(filename,100,gstarS_temperature,gstarS)

c set up the array of x values to calcuate taacs
      do x1=0,99
        if (x1.le.90) then
          x_taacs(x1+1) = 10.d0 + dble(x1)
        else
          x_taacs(x1+1) = 100.d0 + 100.d0*dble(x1-90)
        endif
      enddo

c We initialize the logical variable that keeps track of which ODE's have been solved
c and at what value of xf the solution appears.
      do x1=1, ndmparticles
        xf(x1) = 1.d0
        found_xf(x1) = .false.
        solved_ODE(x1) = .false.
      enddo

      return
      end



c-------------------------------------------------------------------------c
      function get_Y_eq(mass,temp,dof,rel_dofs)
c-------------------------------------------------------------------------c
c                                                                         c
c  This function returns the equilibrium denstiy of DM particle i given   c
c  the temperature x and the number of relativistic degrees of freedom.   c
c                                                                         c
c  Since x is defined by the mass of DM particle 1 over the temperature   c
c  this is where the ratio of dm_mass(i)/dm_mass(1) come from.            c
c                                                                         c
c  Generally before this is run one calls get_gstar(temp) and then Y_eq   c
c  can be calculated for all DM particles.                                c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      include 'maddm.inc'

c input parameters
      double precision mass, temp, dof, rel_dofs

c Y_eq = 45/(4*pi^4) g_i/g_*S x^2 K_2(x), with x=m_i/T
c K_2(x) is expanded for large values of x ( > 10)
      if (mass.gt.5.d0*temp) then
        if (mass/temp.gt.706)then
           get_Y_eq = 0d0
        else
        get_Y_eq = 45.d0/(4.d0*pi**4)*(dof/rel_dofs)*(mass/temp)**(1.5d0)
     .        * dsqrt(pi/2.d0)*dexp(-mass/temp) * (1.d0 + 15.d0/(8.d0*(mass/temp))
     .        + 105.d0/(128.d0*(mass/temp)**2) - 315.d0/(1024.d0*(mass/temp)**3))
      endif
      else if (mass.gt.temp/5.d0) then
        get_Y_eq = 45.d0/(4.d0*pi**4)*(dof/rel_dofs)*(mass/temp)**2
     .        * bessk(2,mass/temp)
      else
        get_Y_eq = 45.d0/(2.d0*pi**4)*(dof/rel_dofs)
      endif

c check to see if the result is a nan.  If so then set Y_eq to 0
      if (.not.(get_Y_eq.eq.get_Y_eq)) get_Y_eq = 0.d0

      return
      end



c-------------------------------------------------------------------------c
      function get_gstar(T)
c-------------------------------------------------------------------------c
c                                                                         c
c  This function returns the number of relativisitic degrees of freedom   c
c  that is calculated for each particle in the model given the            c
c  temperature T using the gstar_boson and gstar_fermion arrays.          c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      include 'maddm.inc'

c input parameters
      double precision T

c parameters used in this subroutine only
      double precision x
      integer gstar_index

c initialzing get_gstar
      get_gstar = 0.d0

c loop over all the particles in the model
      do x1=1, nparticles

c if the mass is effectively zero then the particle is fully relativisitic
        if ((masses(x1)/T).lt.1.0d-2) then
          get_gstar = get_gstar + dofs(x1)*(7.d0/8.d0)**isfermion(x1)

c for masses \approx T then we interpolate gstar_fermion and gstar_boson to
c get the appropriate number of degrees of freedom
        else if ((masses(x1)/T).lt.1.0d2) then
          x = masses(x1)/T
          gstar_index = getindex(x,gstar_x,401)
          if (isfermion(x1).eq.1) then
            get_gstar = get_gstar + dofs(x1)*lineint(x,gstar_x,gstar_fermion,401)
          else if (isfermion(x1).eq.0) then
            get_gstar = get_gstar + dofs(x1)*lineint(x,gstar_x,gstar_boson,401)
          endif
        endif

c any mass that is too large does not contribute to the total entropy density
      enddo 

      return
      end



c-------------------------------------------------------------------------c
      subroutine calculate_Wij_ann()
c-------------------------------------------------------------------------c
c                                                                         c
c  This function calculates for each pair of DM particles \chi_i, \chi_j  c
c  the sum of all the matrix elements squared, integrated over phase      c
c  space over a range of center of mass energies determined by beta.      c
c                                                                         c
c  The setpsize for beta varies depending on the change in consecutive    c
c  values of Wij.  This helps to pick up narrow width ressonances.        c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      include 'maddm.inc'
      include 'coupl.inc'

c parameters used in this subroutine only
      double precision Wij_hold, Wij_total, beta_step
      logical flag1, flag2, flag3

      integer i, ii
      double precision next_beta(nres), tmp_beta
      integer res_step(nres)
      integer nres_points_local, max_width_factor

c output to the screen to show progress
      if (print_out) then
        write(*,*) 'Calculating all the Wijs'
        flag1 = .true.
        flag2 = .true.
        flag3 = .true.
      endif

      nres_points_local = 100
      max_width_factor = 6 ! this is e^N-1 times the width for the further point

      do ii=1, nres
         res_step(ii) = - nres_points_local
         tmp_beta = -1d0
         if (beta_res(ii).lt.0d0) then
            res_step(ii) = nres_points_local +1
            cycle
         endif
         do while (tmp_beta.lt.0d0)
            tmp_beta = beta_res(ii) - beta_res_width(ii)*(DEXP(1d0*res_step(ii)*max_width_factor/nres_points_local)-1)
            res_step(ii) = res_step(ii) + 1
         enddo
c         write(*,*) 'resonances',ii,beta_res(ii), beta_res_width(ii), res_step(ii)
      enddo


c setting up step sizes
      beta_step = 1e-3
      beta_step_min =  1e-6
      beta_step_max = 1e-2
c starting the calculation of the Wij's
      nWij = 1
      betas(nWij) = 0.001d0
      Wij_hold = get_Wij_ann(betas(nWij))

c continue to do beta steps until we are near beta = 1
      do while (betas(nWij).le.0.99d0-beta_step)
        nWij = nWij+1
        tmp_beta = betas(nWij-1) + beta_step
        do ii=1, nres
           if (res_step(ii).le.0)then
              next_beta(ii) = beta_res(ii) - beta_res_width(ii)*(DEXP(1d0*res_step(ii)*max_width_factor/nres_points_local)-1)
           elseif(res_step(ii).le.nres_points_local)then
               next_beta(ii) = beta_res(ii) + beta_res_width(ii)*(DEXP(1d0*res_step(ii)*max_width_factor/nres_points_local)-1)
           else
              next_beta(ii) = 1d0
            endif
            tmp_beta = min(tmp_beta, next_beta(ii),1d0)
         enddo
        betas(nWij) = tmp_beta
        Wij_total = get_Wij_ann(betas(nWij))


c if the new Wij is 100% larger than the previous Wij, then reduce the stepsize
        if (((dabs(Wij_total-Wij_hold)/Wij_hold).gt.1.d0).and.(beta_step.gt.beta_step_min)) then
c           write(*,*) 'reduce', betas(nWij), beta_step, ((dabs(Wij_total-Wij_hold)/Wij_hold)), Wij_total, Wij_hold

           if (((dabs(Wij_total-Wij_hold)/Wij_hold).gt.2.d0))then
c              write(*,*) 'rewind',  betas(nWij-1)
              nWij = nWij-1
              beta_step = beta_step/2.d0
             cycle
           endif

          beta_step = beta_step/2.d0
          if (beta_step .lt. beta_step_min) then
            beta_step = beta_step_min
          endif
        endif

c if the new Wij is within 10% of the previous Wij, then increase the stepsize
        if (((dabs(Wij_total-Wij_hold)/Wij_hold).lt.0.1d0).and.(beta_step.lt.beta_step_max)) then
          beta_step = beta_step*2d0
          if (beta_step .gt. beta_step_max) then
            beta_step = beta_step_max
          endif
c          write(*,*) 'increase', betas(nWij), beta_step
        endif

         do ii=1, nres
            if (next_beta(ii).eq.tmp_beta) then
               res_step(ii) = res_step(ii)+1
               beta_step = max(beta_step,1e-3)
c               if (beta_step.eq.1e-3) then
c                  write(*,*) 'reset'
c               endif
            endif
         enddo


c Hold the Wij_value so that if the Wij integral is too small we have a stored value ot use
        Wij_hold = Wij_total

c Simple way to update the user on the progress of calculating the Wij's
        if (print_out) then
          if ((betas(nWij).gt.0.25d0).and.(flag1)) then
            write(*,*) '25% done'
            flag1 = .false.
          else if ((betas(nWij).gt.0.50d0).and.(flag2)) then
            write(*,*) '50% done'
            flag2 = .false.
          else if ((betas(nWij).gt.0.75d0).and.(flag3)) then
            write(*,*) '75% done'
            flag3 = .false.
          endif
        endif

      enddo

c Calculate the final Wij at beta = 0.999
      nWij = nWij+1
      betas(nWij) = 0.999d0
      Wij_total = get_Wij_ann(betas(nWij))

      if (print_out)  write(*,*) 'Done!', nWij
      Wij_ann_calc = .true.

      return
      end



c-------------------------------------------------------------------------c
      function get_Wij_ann(beta)
c-------------------------------------------------------------------------c
c                                                                         c
c  Given the input velocity beta, this function integrates the matrix     c
c  elements generated by MadGraph over the solid angle for all            c
c  combinations of \chi_i, \chi_j to all final states.                    c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      include 'maddm.inc'

c input parameters
      double precision beta

c external function used in this subroutine
      external Wij_ann_integrand

c set beta_pass to pass beta to Wij_integrand
      beta_pass = beta
c loop over all the dm particles
      do x1=1, ndmparticles
        do x2=x1, ndmparticles
          dmi = x1
          dmj = x2

c integrate the 2 -> 2 process over x = (cos(theta)+1)/2.d0
          call romberg(Wij_ann_integrand, 0.d0, 1.d0, Wij_ann(x1,x2,nWij), eps_wij, iter_wij)

c the Wijs for (DMi, DMj) are the same for (DMj, DMi)
          if (x1.ne.x2) Wij_ann(x2,x1,nWij) = Wij_ann(x1,x2,nWij)
        enddo
      enddo

c the return of this function will be the sum of all the Wij's 
c so the next beta step can be determined
      get_Wij_ann = 0.d0
      do x1=1, ndmparticles
        do x2=x1, ndmparticles
          get_Wij_ann = get_Wij_ann + Wij_ann(x1,x2,nWij)
        enddo
      enddo

      return
      end

c-------------------------------------------------------------------------c
      function get_Wij_ann_nogrid(beta)
c-------------------------------------------------------------------------c
c                                                                         c
c  Given the input velocity beta, this function integrates the matrix     c
c  elements generated by MadGraph over the solid angle for all            c
c  combinations of \chi_i, \chi_j to all final states.                    c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      include 'maddm.inc'

c input parameters
      double precision beta
      double precision Wij_ann_nogrid(maxdms,maxdms)
c external function used in this subroutine
      external Wij_ann_integrand

c set beta_pass to pass beta to Wij_integrand
      beta_pass = beta

      get_Wij_ann_nogrid =  0d0
c loop over all the dm particles
      do x1=1, ndmparticles
        do x2=x1, ndmparticles
          dmi = x1
          dmj = x2

c integrate the 2 -> 2 process over x = (cos(theta)+1)/2.d0
          call romberg(Wij_ann_integrand, 0.d0, 1.d0, Wij_ann_nogrid(x1,x2), eps_wij, iter_wij)

c the Wijs for (DMi, DMj) are the same for (DMj, DMi)
          if (x1.ne.x2) Wij_ann_nogrid(x2,x1) = Wij_ann_nogrid(x1,x2)
        enddo
      enddo

c the return of this function will be the sum of all the Wij's 
c so the next beta step can be determined
      get_Wij_ann = 0.d0
      do x1=1, ndmparticles
        do x2=x1, ndmparticles
          get_Wij_ann_nogrid = get_Wij_ann_nogrid + Wij_ann_nogrid(x1,x2)
        enddo
      enddo

      return
      end



c-------------------------------------------------------------------------c
      function Wij_ann_integrand(x)
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
      integer p_index

c calculating the center of mass energy squared given the passed variable beta_pass
      beta = beta_pass
      s = (mdm(dmi)+mdm(dmj))**2/(1.d0-beta**2)

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
      p_ext(0,1) = (s+mdm(dmi)**2-mdm(dmj)**2)/(2.d0*dsqrt(s))
      p_ext(0,2) = (s+mdm(dmj)**2-mdm(dmi)**2)/(2.d0*dsqrt(s))
      p_ext(1,1) = 0.d0
      p_ext(1,2) = 0.d0
      p_ext(2,1) = 0.d0
      p_ext(2,2) = 0.d0
      p_ext(3,1) = dsqrt(lambda(s,mdm(dmi)**2,mdm(dmj)**2))/(2.d0*dsqrt(s))
      p_ext(3,2) = -dsqrt(lambda(s,mdm(dmi)**2,mdm(dmj)**2))/(2.d0*dsqrt(s))

c initializing the Wij_integrand
      Wij_ann_integrand = 0.d0

      do x1=1, ann_nprocesses(dmi,dmj)        

c gets the final state masses using the pmass(i)_(j)_(k).inc files
        call getfsmasses_ann(pmass,dmi,dmj,x1)

c check to see if the center of mass energy is enough to produce the final states
        if (s.ge.(pmass(3)+pmass(4))**2) then

c two body phase space subroutine
          call decay2(p_init,pmass(3),pmass(4),pf1,pf2,wgt,ps_point,1)

c          write(*,*) 'diagram: ', x1, beta

          do x2=0,3
            p_ext(x2,3) = pf1(x2)
            p_ext(x2,4) = pf2(x2)
          enddo

c calculate the matrix element      
          msq = smatrix_ann(p_ext,dmi,dmj,x1)

c If the process being evaluated has identical initial state particles then we 
c preemptively divide by two because eventually we will be integrating over the
c initial state phase space which double counts the number of interactions
          if (ann_process_iden_init(dmi,dmj,x1)) then
            msq = msq/2.d0
          endif

c MadGraph automatically averages over the initial state degrees of freedom.  Our definition
c of Wij requires inital states to be summed, hence factors of dof_dm
          Wij_ann_integrand = dof_dm(dmi)*dof_dm(dmj)*msq*wgt + Wij_ann_integrand

        endif     
      enddo

      return
      end



c-------------------------------------------------------------------------c
      subroutine get_process_index(i, j, k, group, p_index)
c-------------------------------------------------------------------------c
c                                                                         c
c  This subroutine gives the process index giving the initial state       c
c  particle ids (i and j), the process number k, and the process group.   c
c  The output is given as p_index.                                        c
c                                                                         c
c  group = 1 for annihilation processes.                                  c
c  group = 2 for dm2dm processes.                                         c
c  group = 3 for DM/SM scattering processes.                              c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      include 'maddm.inc'

c input variables
      integer i, j, k, group, p_index

c initialize the index
      p_index = 0

c finding p_index is different depending on which process group we want
c annihilation processes
      if (group.eq.1) then

c loop over all combinations that happen before the desired initial combinations and sum up
c the number of processes
        do x1=1, ndmparticles
          do x2=x1, ndmparticles
            if ((x1.lt.i).or.(x2.lt.j)) then
              p_index = p_index + ann_nprocesses(x1,x2)
            endif
            if ((x1.eq.i).and.(x2.eq.j)) then
              p_index = p_index + k
              return
            endif
          enddo
        enddo

c DM -> DM processes
      else if (group.eq.2) then
        p_index = p_index + ann_num_processes

c loop over all combinations that happen before the desired initial combinations and sum up
c the number of processes
        do x1=1, ndmparticles
          do x2=x1, ndmparticles
            if ((x1.lt.i).or.(x2.lt.j)) then
              p_index = p_index + dm2dm_nprocesses(x1,x2)
            endif
            if ((x1.eq.i).and.(x2.eq.j)) then
              p_index = p_index + k
              return
            endif
          enddo
        enddo

c DM/SM scattering processes
      else if (group.eq.3) then
        p_index = p_index + ann_num_processes + dm2dm_num_processes

c loop over all combinations that happen before the desired initial combinations and sum up
c the number of processes
        do x1=1, ndmparticles
          do x2=1, ndmparticles
            if ((x1.lt.i).or.(x2.lt.j)) then
              p_index = p_index + scattering_nprocesses(x1,x2)
            endif
            if ((x1.eq.i).and.(x2.eq.j)) then
              p_index = p_index + k
              return
            endif
          enddo
        enddo
      endif

c we shouldn't get here but just in case
      write(*,*) 'Something went wrong in get_process_index'
      read(*,*)

      return
      end
