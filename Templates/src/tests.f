c-------------------------------------------------------------------------c
      subroutine dof_check(Tmin, Tmax, nsteps, logscale)
c-------------------------------------------------------------------------c
c                                                                         c
c  This subroutine prints out the number of relativistic degrees of       c
c  freedom over a range of temperatures defined by Tmin and Tmax with     c
c  nsteps number of steps.  logscale chooses whether the scale will be    c
c  logarithmic or linear. (logscale = 1 for log, logscale = 0 for linear) c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      include 'maddm.inc'

c input parameters
      double precision Tmin, Tmax
      integer nsteps, logscale

c parameters used in this subroutine only
      double precision temp

c open the output file
      open(unit=42, file='output/dof_check.txt')

c relativistic degrees of freedom check output
      write(42,fmt='(A25)') '#   temp           g_star'

      do x1=1, nsteps
        if (logscale.eq.0) then
          temp = (Tmax-Tmin)*dble(x1-1)/dble(nsteps-1) + Tmin
        else if (logscale.eq.1) then
          temp = 10.d0**(dlog(Tmax/Tmin)/dlog(10.d0)*dble(x1-1)/dble(nsteps-1) + dlog(Tmin)/dlog(10.d0))
        else
          write(*,*) 'Invalid value of logscale'
          write(*,*) 'logscale = 1 for log scale'
          write(*,*) 'logscale = 0 for linear scale'
        endif
        write(42,fmt='(2(ES16.8))') temp, get_gstar(temp)
      enddo

      close(42)

      return
      end



c-------------------------------------------------------------------------c
      subroutine Wij_check()
c-------------------------------------------------------------------------c
c                                                                         c
c  This subroutine prints out the Wij values over a range of betas.       c
c  This should be called after calculate_Wij is called.                   c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      include 'maddm.inc'

c parameters used in this subroutine only
      character(50) output_file, pname
      integer proc_index, pname_index

c Wij check for the annihilation processes
      do x1=1,nvar
        do x2=x1,nvar

c open up the output file
          output_file = 'output/wij_check_'//dm_names(x1)(1:dm_index(x1))//'_'//dm_names(x2)(1:dm_index(x2))//'.txt'
          open(unit=42, file=output_file)          

c Wij check output
          write(42,fmt='(A18,I4)') '# number of betas:',nWij
          write(42,fmt='(A22)') '#  beta            Wij'

          do x3=1, nWij
            write(42,fmt='(2(ES16.8))') betas(x3), Wij_ann(x1,x2,x3)
          enddo

          close(42)

        enddo
      enddo

c Wij check for the DM -> DM processes
      do x1=1,ndmparticles
        do x2=x1,ndmparticles
          do x3=1,dm2dm_nprocesses(x1,x2)
            
c open up the output file
            call get_process_index(x1,x2,x3,2,proc_index) 
            pname = process_names(proc_index)
            pname_index = index(pname,' ')-1
            output_file = 'output/wij_check_'//pname(1:pname_index)//'.txt'
            open(unit=42, file=output_file)

c Wij check output
            write(42,fmt='(A18,I4)') '# number of betas:',nWij
            write(42,fmt='(A22)') '#  beta            Wij'

            do x4=1, nWij
              write(42,fmt='(2(ES16.8))') betas(x4), Wij_dm2dm(x1,x2,x3,x4)
            enddo

            close(42)
          enddo
        enddo
      enddo

C c We aren't using the DM/SM scattering processes in MadDM so we can comment this out
C c Wij check for the DM/SM processes
C       do x1=1,ndmparticles
C         do x2=1,ndmparticles
C           do x3=1,scattering_nprocesses(x1,x2)
C             
C c open up the output file
C             call get_process_index(x1,x2,x3,3,proc_index) 
C             pname = process_names(proc_index)
C             pname_index = index(pname,' ')-1
C             output_file = 'output/wij_check_'//pname(1:pname_index)//'.txt'
C             open(unit=42, file=output_file)
C 
C c Wij check output
C             write(42,fmt='(A18,I4)') '# number of betas:',nWij
C             write(42,fmt='(A22)') '#  beta            Wij'
C 
C             do x4=1, nWij
C               write(42,fmt='(2(ES16.8))') betas(x4), Wij_scattering(x1,x2,x3,x4)
C             enddo
C 
C             close(42)
C           enddo
C         enddo
C       enddo

      return
      end



c-------------------------------------------------------------------------c
      subroutine taacs_check(xmin, xmax, nsteps, logscale)
c-------------------------------------------------------------------------c
c                                                                         c
c  This subroutine prints out the thermally averaged annihilation cross   c
c  section over a range of x values defined by xmin and xmax.  nsteps     c
c  is the number of steps that are calculated and logscale chooses        c
c  whether the scale will be logrithmic or linear.  (logscale = 1 for     c
c  log, logscale = 0 for linear)                                          c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      include 'maddm.inc'

c input parameters
      double precision xmin, xmax
      integer nsteps, logscale

c parameters used in this subroutine only
      double precision x
      character(50) output_file, pname
      integer proc_index, pname_index, i, j, k

c taacs check for the annihilation processes
      do x1=1,nvar
        do x2=x1,nvar

c open up the output file          
          output_file = 'output/taacs_check_'//dm_names(x1)(1:dm_index(x1))//'_'//dm_names(x2)(1:dm_index(x2))//'.txt'
          open(unit=42, file=output_file)

c thermally averaged annihilation cross section check
          write(42,fmt='(A23)') '#  x              taacs'
          do x3=1, nsteps
            if (logscale.eq.0) then
              x = (xmax-xmin)*dble(x3-1)/dble(nsteps-1) + xmin
            else if (logscale.eq.1) then
              x = 10.d0**(dlog(xmax/xmin)/dlog(10.d0)*dble(x3-1)/dble(nsteps-1) + dlog(xmin)/dlog(10.d0))
            else
              write(*,*) 'Invalid value of logscale'
              write(*,*) 'logscale = 1 for log scale'
              write(*,*) 'logscale = 0 for linear scale'
            endif
            if (relic_canonical) then
              write(42,fmt='(2(ES16.8))') x, taacs_canon(x)
            else
              write(42,fmt='(2(ES16.8))') x, taacs_coupled(x,x1,x2,0,1)
            endif
          enddo

          close(42)
        enddo
      enddo

c The rest of these are only used in the coupled scenario
      if (.not.relic_canonical) then
        
c taacs check for the DM -> DM processes
        do x1=1,nvar
          do x2=x1,nvar
            do x3=1,dm2dm_nprocesses(x1,x2)
            
c open up the output file
              call get_process_index(x1,x2,x3,2,proc_index) 
              pname = process_names(proc_index)
              pname_index = index(pname,' ')-1
              output_file = 'output/taacs_check_'//pname(1:pname_index)//'.txt'
              open(unit=42, file=output_file)        

c thermally averaged annihilation cross section check
              write(42,fmt='(A23)') '#  x              taacs'

              do x4=1, nsteps
                if (logscale.eq.0) then
                  x = (xmax-xmin)*dble(x4-1)/dble(nsteps-1) + xmin
                else if (logscale.eq.1) then
                  x = 10.d0**(dlog(xmax/xmin)/dlog(10.d0)*dble(x4-1)/dble(nsteps-1) + dlog(xmin)/dlog(10.d0))
                else
                  write(*,*) 'Invalid value of logscale'
                  write(*,*) 'logscale = 1 for log scale'
                  write(*,*) 'logscale = 0 for linear scale'
                endif

c get the process information to get the initial states and call the taacs_dm2dm integrator
                write(42,fmt='(2(ES16.8))') x, taacs_coupled(x,x1,x2,x3,2)
              enddo

              close(42)
            enddo
          enddo
        enddo

C c taacs check for the DM/SM scattering processes
C         do x1=1,nvar
C           do x2=1,nvar
C             do x3=1,scattering_nprocesses(x1,x2)
C             
C c open up the output file
C               call get_process_index(x1,x2,x3,3,proc_index) 
C               pname = process_names(proc_index)
C               pname_index = index(pname,' ')-1
C               output_file = 'output/taacs_check_'//pname(1:pname_index)//'.txt'
C               open(unit=42, file=output_file)        
C 
C c thermally averaged annihilation cross section check
C               write(42,fmt='(A23)') '#  x              taacs'
C 
C               do x4=1, nsteps
C                 if (logscale.eq.0) then
C                   x = (xmax-xmin)*dble(x4-1)/dble(nsteps-1) + xmin
C                 else if (logscale.eq.1) then
C                   x = 10.d0**(dlog(xmax/xmin)/dlog(10.d0)*dble(x4-1)/dble(nsteps-1) + dlog(xmin)/dlog(10.d0))
C                 else
C                   write(*,*) 'Invalid value of logscale'
C                   write(*,*) 'logscale = 1 for log scale'
C                   write(*,*) 'logscale = 0 for linear scale'
C                 endif
C 
C c get the process information to get the initial states and call the taacs_dm2dm integrator
C                 write(42,fmt='(2(ES16.8))') x, taacs_coupled(x,x1,x2,x3,3)
C               enddo
C 
C               close(42)
C             enddo
C           enddo
C         enddo
      endif

      return
      end



c-------------------------------------------------------------------------c
      subroutine odeint_check(canon_or_coupled)
c-------------------------------------------------------------------------c
c                                                                         c
c  During the integration of the chemistry equation we can print out the  c
c  results from the ODE integrator to a file to actually see when         c
c  freezeout occurs.                                                      c
c                                                                         c
c  canon_or_coupled = 0, prints output for the sum of all DM particles    c
c  canon_or_coupled = 1, prints output for each individual DM particle    c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      include 'maddm.inc'

c input parameters
      integer canon_or_coupled

c parameters used in this subroutine only
      double precision Y_eq
      character*50 output_file

c checks to make sure that the ODE information has been saved
      if (dxsav.eq.0.d0) then
        write(*,*) '---------------------------------------------------------'
        write(*,*) 'odeint_check:'
        write(*,*) '  ODE information not saved!!!'
        write(*,*) '  dxsav, kount, and kmax need to be set in relicdensity.f'
        write(*,*) '---------------------------------------------------------'
        return
      endif

c output of the new relic density
      if (canon_or_coupled.eq.1) then
        do x1=1,nvar

c output of the ode integration
          output_file = 'output/odeint_coupled'//dm_names(x1)(1:dm_index(x1))//'.txt'
          open(unit=42, file=output_file)

          write(42,fmt='(A21)') '#  x                Y'
          do x2=1, kount
            write(42,fmt='(2(ES16.8))') xp(x2), yp(x1,x2)
          enddo
          close(42)

c output file for Y_eq to compare with the ode integration
          output_file = 'output/odeint_coupled_Yeq_'//dm_names(x1)(1:dm_index(x1))//'.txt'
          open(unit=42, file=output_file)
          write(42,fmt='(A23)') '#  x                Yeq'

c loop over x=1 to 100
          do x2=1, 100

c for each point we need to calculate Y_eq
            g_star_S = get_gstar(mdm(1)/dble(x2))
            Y_eq = get_Y_eq(mdm(x1),mdm(1)/x2,dof_total(x1),g_star_S)

c write the results to the output file          
            write(42,fmt='(2(ES16.8))') dble(x2), Y_eq
          enddo

          close(42)
        enddo

c output of the canonical relic density
      else if (canon_or_coupled.eq.0) then

c output of the ode integration
        output_file = 'output/odeint_canon.txt'
        open(unit=42, file=output_file)

        write(42,fmt='(A21)') '#  x                Y'
        do x1=1, kount
          write(42,fmt='(2(ES16.8))') xp(x1), yp(1,x1)
        enddo
        close(42)

c output file for the Y_eq to compare with the ode integration
        output_file = 'output/odeint_canon_Yeq.txt'
        open(unit=42, file=output_file)
        write(42,fmt='(A23)') '#  x                Yeq'

c loop over x=1 to 100
        do x1=1, 100

c for each point we need to calculate the sum of Y_eq for all particles
          g_star_S = get_gstar(mdm(1)/dble(x1))
          Y_eq = 0.d0
          do x2=1,nvar 
            Y_eq = Y_eq + get_Y_eq(mdm(x2),mdm(1)/x1,dof_total(x2),g_star_S)
          enddo

c write the results to the output file          
          write(42,fmt='(2(ES16.8))') dble(x1), Y_eq
        enddo

        close(42)
      endif          

      return
      end



c-------------------------------------------------------------------------c
      subroutine cross_check_scan_all(x_min, x_max, nsteps, p_or_E, logscale)
c-------------------------------------------------------------------------c
c                                                                         c
c  This subroutine calculates the cross section for all of the            c
c  annihilation processes generated in madgraph over a range of inital    c
c  state momenta and prints everything to a file for each process.        c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      include 'maddm.inc'
      include 'coupl.inc'

c input parameters
      double precision x_min, x_max
      integer nsteps, p_or_E, logscale

c loops over all initial state combinations.
      do x1=1,ndmparticles
        do x2=x1,ndmparticles

c loop over the individual annihilation processes
          do x3=1,ann_nprocesses(x1,x2)
            call cross_check_scan(x1,x2,x3,1,x_min,x_max,nsteps,p_or_E,logscale)
          enddo

c loop over the DM -> DM processes
          do x3=1,dm2dm_nprocesses(x1,x2)
            call cross_check_scan(x1,x2,x3,2,x_min,x_max,nsteps,p_or_E,logscale)
          enddo
        enddo

C c loop over the DM/SM scattering processes
C         do x2=1,ndmparticles
C           do x3=1,scattering_nprocesses(x1,x2)
C             call cross_check_scan(x1,x2,x3,3,x_min,x_max,nsteps,p_or_E,logscale)
C           enddo
C         enddo
      enddo

      return
      end



c-------------------------------------------------------------------------c
      subroutine cross_check_scan(i, j, k, group, x_min, x_max, nsteps, p_or_E, logscale)
c-------------------------------------------------------------------------c
c                                                                         c
c  This subroutine calculates the annihilation cross section for each     c
c  individual process for a given range of initial momenta used for       c
c  both initial state particles.  This output can then be easily          c
c  compared to other calculations (e.g. from CalcHEP)                     c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      include 'maddm.inc'
      include 'coupl.inc'

c input parameters
      double precision x_min, x_max
      integer i, j, k, group, nsteps, p_or_E, logscale

c parameters used in this subroutine only
      double precision xinit, roots, cross_section, get_cross_section
      integer process_index, pname_index
      character*50 pname, output_file

c opening the appropriate output file
      call get_process_index(i,j,k,group,process_index)
      pname = process_names(process_index)
      pname_index = index(pname,' ') - 1
      output_file = 'output/cross_check_'//pname(1:pname_index)//'.txt'
      open(unit=42, file=output_file)

c header for the output file
      write(42,fmt='(A14,I3,A15,A10)') '# Process no.:',process_index,' Process name: ',pname(1:pname_index)
      if (p_or_E.eq.1) then
        write(42,fmt='(A35)') '# pinit (GeV)    cross section (pb)'
      else if (p_or_E.eq.0) then
        write(42,fmt='(A35)') '# Einit (GeV)    cross section (pb)'
      endif

      do x1=1, nsteps
        if (logscale.eq.0) then
          xinit = (x_max-x_min)*dble(x1-1)/dble(nsteps-1) + x_min
        else if (logscale.eq.1) then
          xinit = 10.d0**(dlog(x_max/x_min)/dlog(10.d0)*dble(x1-1)/dble(nsteps-1) + dlog(x_min)/dlog(10.d0))
        else
          write(*,*) 'Invalid value of logscale'
          write(*,*) 'logscale = 1 for log scale'
          write(*,*) 'logscale = 0 for linear scale'
          return
        endif

c get the final state masses using the pmass(i)_(j)_(k).inc files
        if (group.eq.1) then
          call getfsmasses_ann(pmass,i,j,k)
        else if (group.eq.2) then
          call getfsmasses_dm2dm(pmass,i,j,k)
C         else if (group.eq.3) then
C           call getfsmasses_scattering(pmass,i,j,k)
        endif

        if (p_or_E.eq.1) then
          roots = dsqrt(pmass(1)**2 + xinit**2) + dsqrt(pmass(2)**2 + xinit**2)
        else if (p_or_E.eq.0) then
          roots = dsqrt(4.d0*xinit**2 - (dsqrt(xinit**2-pmass(1)**2)-dsqrt(xinit**2-pmass(2)**2))**2)
        else
          write(*,*) 'Invalid value of p_or_E'
          write(*,*) 'p_or_E = 1 for initial momentum'
          write(*,*) 'p_or_E = 0 for initial energy'
          return
        endif

        cross_section = get_cross_section(roots,i,j,k,group)
        write(42,fmt='(2(ES16.8))') xinit, cross_section
      enddo

      close(42)
      return
      end

c-------------------------------------------------------------------------c
      subroutine cross_check_all_process(x1init, x2init, p_or_E)
c-------------------------------------------------------------------------c
c                                                                         c
c  This subroutine calculates the annihilation cross section for each     c
c  individual process for a given the initial momentum or energy for both c
c  initial state particles x1init and x2init, which is chosen by p_or_E.  c
c  This output can then be easily compared to other calculations          c
c  (e.g. from CalcHEP and MadGraph)                                       c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      include 'maddm.inc'
      include 'coupl.inc'

c input parameters
      double precision x1init, x2init
      integer p_or_E

c parameters used in this subroutine only
      double precision roots, cross_section, get_cross_section, total_cross
      character*30 pname
      integer process_index, pname_index

c initialize the total cross section
      total_cross = 0.d0

c header for the output to the screen
      if (p_or_E.eq.1) then
        write(*,*) '# initial momentum for particle 1 = ',x1init,' GeV'
        write(*,*) '# initial momentum for particle 2 = ',x2init,' GeV'
      else if (p_or_E.eq.0) then
        write(*,*) '# initial energy for particle 1 = ',x1init,' GeV'
        write(*,*) '# initial energy for particle 2 = ',x2init,' GeV'
      else
        write(*,*) 'Invalid value of p_or_E'
        write(*,*) 'p_or_E = 1 for initial momentum'
        write(*,*) 'p_or_E = 0 for initial energy'
        return
      endif
      write(*,*) '# DM1 DM2 process no.    process name    cross section (pb)'


c loop over all the pairs of DM initial states
        do x1=1,ndmparticles
         do x2=x1,ndmparticles
c loop over all the annihilation diagrams for each DM particle pair
       do x3=1,ann_nprocesses(x1,x2)


c gets the final state masses using the pmass(i)_(j)_(k).inc files

             call getfsmasses_ann(pmass,x1,x2,x3)


            call get_process_index(x1,x2,x3,1,process_index)

             if (p_or_E.eq.1) then
                  roots = dsqrt((dsqrt(pmass(1)**2 + x1init**2) + dsqrt(pmass(2)**2 + x2init**2))**2 - (x1init-x2init)**2)
             else if (p_or_E.eq.0) then
                  roots = dsqrt((x1init+x2init)**2 - (dsqrt(x1init**2-pmass(1)**2) - dsqrt(x2init**2-pmass(2)**2))**2)
             endif

            pname = process_names(process_index)
            pname_index = index(pname,' ') - 1
            cross_section = get_cross_section(roots,x1,x2,x3,1)

            write(*,fmt='(I5,I4,I7,2X,A20,4X,ES16.8)') x1,x2,x3,pname(1:pname_index), cross_section

c add up all the individual cross sections
            total_cross = total_cross + cross_section
          enddo

c loop over all the DM -> DM diagrams for each DM particle pair
          do x3=1,dm2dm_nprocesses(x1,x2)
            call get_process_index(x1,x2,x3,2,process_index)
            pname = process_names(process_index)
            pname_index = index(pname,' ') - 1

             if (p_or_E.eq.1) then
                  roots = dsqrt((dsqrt(pmass(1)**2 + x1init**2) + dsqrt(pmass(2)**2 + x2init**2))**2 - (x1init-x2init)**2)
             else if (p_or_E.eq.0) then
                  roots = dsqrt((x1init+x2init)**2 - (dsqrt(x1init**2-pmass(1)**2) - dsqrt(x2init**2-pmass(2)**2))**2)
             endif

            cross_section = get_cross_section(roots,x1,x2,x3,2)

            write(*,fmt='(I5,I4,I7,2X,A20,4X,ES16.8)') x1,x2,x3,pname(1:pname_index), cross_section

c add up all the individual cross sections
            total_cross = total_cross + cross_section
          enddo          
        enddo
      enddo

      write(*,fmt='(A31,ES16.8)') ' # Total cross section (in pb): ',total_cross

      return
      end



c-------------------------------------------------------------------------c
      function cross_check_process(i, j, k, group, x1init, x2init, p_or_E)
c-------------------------------------------------------------------------c
c                                                                         c
c  This subroutine calculates the annihilation cross section for each     c
c  individual process for a given the initial momentum pinit used for     c
c  both initial state particles.  This output can then be easily          c
c  compared to other calculations (e.g. from CalcHEP)                     c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      include 'maddm.inc'
      include 'coupl.inc'

c input parameters
      integer i, j, k, group, p_or_E
      double precision x1init, x2init

c parameters used in this subroutine only
      double precision roots, cross_section, get_cross_section

c gets the final state masses using the pmass(i)_(j)_(k).inc files
      if (group.eq.1) then
        call getfsmasses_ann(pmass,i,j,k)
      else if (group.eq.2) then
        call getfsmasses_dm2dm(pmass,i,j,k)
C       else if (group.eq.3) then
C         call getfsmasses_scattering(pmass,i,j,k)
      endif

      if (p_or_E.eq.1) then
        roots = dsqrt((dsqrt(pmass(1)**2 + x1init**2) + dsqrt(pmass(2)**2 + x2init**2))**2 - (x1init-x2init)**2)
      else if (p_or_E.eq.0) then
        roots = dsqrt((x1init+x2init)**2 - (dsqrt(x1init**2-pmass(1)**2) - dsqrt(x2init**2-pmass(2)**2))**2)
      else
        write(*,*) 'Invalid value of p_or_E'
        write(*,*) 'p_or_E = 1 for initial momentum'
        write(*,*) 'p_or_E = 0 for initial energy'
      cross_check_process = 0d0
        return
      endif

      cross_check_process = get_cross_section(roots,i,j,k,group)

      return
      end



c-------------------------------------------------------------------------c
      function get_cross_section(roots, i, j, k, group)
c-------------------------------------------------------------------------c
c                                                                         c
c  This test function returns the annihilation cross section of for the   c
c  process dm_i + dm_j -> process_k with a center of mass energy equal    c
c  to roots.                                                              c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      include 'maddm.inc'
      include 'coupl.inc'

c input parameters
      double precision roots
      integer i, j, k, group

c parameters to pass to crossijk
      double precision roots_pass
      common/cross_pass1/ roots_pass
      integer i_pass, j_pass, k_pass, group_pass
      common/cross_pass2/ i_pass, j_pass, k_pass, group_pass
      external crossijk

c function initiation
      double precision get_cross_section

c setting up the input varibales to pass to crossijk
      roots_pass = roots
      i_pass = i
      j_pass = j
      k_pass = k
      group_pass = group

c integrate the 2 -> 2 process over x = (cos(theta)+1)/2.d0
      call romberg(crossijk,0.d0,1.d0,get_cross_section,eps_wij,iter_wij)      

      return
      end



c-------------------------------------------------------------------------c
      function crossijk(x)
c-------------------------------------------------------------------------c
c                                                                         c
c  This function is the integrand called by the get_cross_section         c
c  function to calculate the cross section of a given process.            c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      include 'maddm.inc'
      include 'coupl.inc'

c input parameters
      double precision x

c parameters used in this subroutine only
      double precision roots_pass
      common/cross_pass1/ roots_pass
      integer i_pass, j_pass, k_pass, group_pass
      common/cross_pass2/ i_pass, j_pass, k_pass, group_pass

c parameters for the phase space point and matrix element
      double precision msq, p_init(0:3), p_ext(0:3,4)
      double precision ps_point(2), wgt, mf1, mf2

c initialize the function
      crossijk = 0.d0

c two body phase space point (cos(theat),phi)
      ps_point(1) = x
      ps_point(2) = 0.d0


c gets the final state masses using the pmass(i)_(j)_(k).inc files
      if (group_pass.eq.1) then
        call getfsmasses_ann(pmass,i_pass,j_pass,k_pass)
      else if (group_pass.eq.2) then
        call getfsmasses_dm2dm(pmass,i_pass,j_pass,k_pass)
C       else if (group_pass.eq.3) then
C         call getfsmasses_scattering(pmass,i_pass,j_pass,k_pass)
      endif

c center of mass 4-momentum
      p_init(0) = roots_pass
      p_init(1) = 0.d0
      p_init(2) = 0.d0
      p_init(3) = 0.d0

c setting up the external particle momenta for the matrix element
      p_ext(0,1) = (roots_pass**2+pmass(1)**2-pmass(2)**2)/(2.d0*roots_pass)
      p_ext(0,2) = (roots_pass**2+pmass(1)**2-pmass(2)**2)/(2.d0*roots_pass)
      p_ext(1,1) = 0.d0
      p_ext(1,2) = 0.d0
      p_ext(2,1) = 0.d0
      p_ext(2,2) = 0.d0
      p_ext(3,1) = dsqrt(lambda(roots_pass**2,pmass(1)**2,pmass(2)**2))/(2.d0*roots_pass)
      p_ext(3,2) = -dsqrt(lambda(roots_pass**2,pmass(1)**2,pmass(2)**2))/(2.d0*roots_pass)


c gets the final particle 4-momenta       
      mf1 = pmass(3)
      mf2 = pmass(4)

      if (roots_pass.ge.(mf1+mf2)) then

c two body phase space subroutine
        call decay2(p_init,mf1,mf2,pf1,pf2,wgt,ps_point,1)

c put the final state 4-momenta into p_ext to feed into MG5 matrix elements
        do x1=0,3
          p_ext(x1,3) = pf1(x1)
          p_ext(x1,4) = pf2(x1)
        enddo

c calculate the matrix element
        if (group_pass.eq.1) then      
          msq = smatrix_ann(p_ext,i_pass,j_pass,k_pass)
        else if (group_pass.eq.2) then
          msq = smatrix_dm2dm(p_ext,i_pass,j_pass,k_pass)
        else
           msq = 0d0
           write(*,*) 'undefined msq'
           stop 1
C         else if (group_pass.eq.3) then
C           msq = smatrix_scattering(p_ext,i_pass,j_pass,k_pass)
        endif

        crossijk = 1.d0/(2.d0*dsqrt(lambda(roots_pass**2,pmass(1)**2,pmass(2)**2)))*msq*wgt*3.894d8
      endif     

      return
      end



c-------------------------------------------------------------------------c
      subroutine matrix_element_check_all(x1init, x2init, p_or_E)
c-------------------------------------------------------------------------c
c                                                                         c
c  This subroutine outputs the value of the matrix element squared over   c
c  a range of cos(theta) for all of the processes with a given initial    c
c  state energy or momentum.                                              c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      include 'maddm.inc'
      include 'coupl.inc'

c input parameters
      double precision x1init, x2init
      integer p_or_E

c loop over all the combinations of initial state particles
      do x1=1,ndmparticles
        do x2=x1,ndmparticles
          
c loop over all the annihilation processes 
          do x3=1,ann_nprocesses(x1,x2)
            call matrix_element_check_process(x1, x2, x3, 1, x1init, x2init, p_or_E)
          enddo

c loop over all the DM -> DM processes
          do x3=1,dm2dm_nprocesses(x1,x2)
            call matrix_element_check_process(x1, x2, x3, 2, x1init, x2init, p_or_E)            
          enddo
        enddo
      enddo

      return
      end



c-------------------------------------------------------------------------c
      subroutine matrix_element_check_process(i, j, k, group, x1init, x2init, p_or_E)
c-------------------------------------------------------------------------c
c                                                                         c
c  This subroutine outputs the matrix element for given initial state     c
c  energy or momentum for a specific process chosen by the initial        c
c  states i and j, and the individual process k over a range of           c
c  cos(thetas) for the final state particle orientation.                  c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      include 'maddm.inc'
      include 'coupl.inc'

c input parameters
      integer i, j, k, group, p_or_E
      double precision x1init, x2init

c parameters for the phase space point and matrix element
      double precision roots, msq, p_init(0:3), p_ext(0:3,4)
      double precision ps_point(2), wgt, mf1, mf2
      character*50 pname, output_file
      integer process_index, pname_index


c gets the final state masses using the pmass(i)_(j)_(k).inc files
      if (group.eq.1) then
        call getfsmasses_ann(pmass,i,j,k)
      else if (group.eq.2) then
        call getfsmasses_dm2dm(pmass,i,j,k)
C       else if (group.eq.3) then
C         call getfsmasses_scattering(pmass,i,j,k)
      endif

c calculate the center of mass energy given the input value for p_or_E
      if (p_or_E.eq.1) then
        roots = dsqrt((dsqrt(pmass(1)**2 + x1init**2) + dsqrt(pmass(2)**2 + x2init**2))**2 - (x1init-x2init)**2)
      else if (p_or_E.eq.0) then
        roots = dsqrt((x1init+x2init)**2 - (dsqrt(x1init**2-pmass(1)**2) - dsqrt(x2init**2-pmass(2)**2))**2)
      else
        write(*,*) 'Invalid value of p_or_E'
        write(*,*) 'p_or_E = 1 for initial momentum'
        write(*,*) 'p_or_E = 0 for initial energy'
        return
      endif

c check to see if the center of mass energy is valid
      if ((p_or_E.eq.0).and.((x1init.lt.pmass(1)).or.(x2init.lt.pmass(2)))) then
        write(*,*) 'Invalid initial state energy'
        return
      endif

c opening the appropriate output file
      call get_process_index(i,j,k,group,process_index)
      pname = process_names(process_index)
      pname_index = index(pname,' ') - 1
      output_file = 'output/matrix_check_'//pname(1:pname_index)//'.txt'
      open(unit=42, file=output_file)

c center of mass 4-momentum
      p_init(0) = roots
      p_init(1) = 0.d0
      p_init(2) = 0.d0
      p_init(3) = 0.d0

c all the matrix elements are evaluated at \phi = 0 
      ps_point(2) = 0.d0

c setting up the external particle momenta for the matrix element
      p_ext(0,1) = (roots**2+pmass(1)**2-pmass(2)**2)/(2.d0*roots)
      p_ext(0,2) = (roots**2+pmass(2)**2-pmass(1)**2)/(2.d0*roots)
      p_ext(1,1) = 0.d0
      p_ext(1,2) = 0.d0
      p_ext(2,1) = 0.d0
      p_ext(2,2) = 0.d0
      p_ext(3,1) = dsqrt(lambda(roots**2,pmass(1)**2,pmass(2)**2))/(2.d0*roots)
      p_ext(3,2) = -dsqrt(lambda(roots**2,pmass(1)**2,pmass(2)**2))/(2.d0*roots)


c gets the final particle 4-momenta       
      mf1 = pmass(3)
      mf2 = pmass(4)

c check if the center of mass energy is large enough for the final state particles
      if (roots.ge.(mf1+mf2)) then

c loop over a bunch of phase space points for cos(theta)
        do x1=0, 200
          ps_point(1) = dble(x1)/200.d0

c two body phase space subroutine
          call decay2(p_init,mf1,mf2,pf1,pf2,wgt,ps_point,1)

c put the final state 4-momenta into p_ext to feed into MG5 matrix elements
          do x2=0,3
            p_ext(x2,3) = pf1(x2)
            p_ext(x2,4) = pf2(x2)
          enddo

c calculate the matrix element
          if (group.eq.1) then
            msq = smatrix_ann(p_ext,i,j,k)
          else if (group.eq.2) then
            msq = smatrix_dm2dm(p_ext,i,j,k)
C           else if (group.eq.3) then
C             msq = smatrix_dm2dm(p_ext,i,j,k)
          endif
        enddo    
      endif

      close(42)

      return
      end



c-------------------------------------------------------------------------c
      subroutine bessel_check(xmin, xmax, nsteps, logscale)
c-------------------------------------------------------------------------c
c                                                                         c
c  This subroutine checks the accuracy of the bessel function             c
c  approximations used in the relic abundance calculation over a range    c
c  of x values defined by xmin and xmax.  We compare the approximations   c
c  to the bessel functions from numerical recipes and it shows where we   c
c  can use the approximations.                                            c
c                                                                         c
c  nsteps is the number of steps in the scan logscale chooses whether     c
c  the scale will be logarithmic or linear. (logscale = 1 for log,        c
c  logscale = 0 for linear)                                               c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      include 'maddm.inc'

c input parameters
      double precision xmin, xmax
      integer nsteps, logscale

c parameters used in this subroutine only
      double precision x, K1_value, K1_approx, K1_diff, K2_value, K2_approx, K2_diff

c open the output file
      open(unit=42, file='output/bessel_text.txt')

c relativistic degrees of freedom check output
      write(42,fmt='(A15,6(A18))') '#             x','K_1(x)','K(1) approx','% diff','K_2(x)',' K_2(x) approx','% diff'

      do x1=1, nsteps
        if (logscale.eq.0) then
          x = (xmax-xmin)*dble(x1-1)/dble(nsteps-1) + xmin
        else if (logscale.eq.1) then
          x = 10.d0**(dlog(xmax/xmin)/dlog(10.d0)*dble(x1-1)/dble(nsteps-1) + dlog(xmin)/dlog(10.d0))
        else
          write(*,*) 'Invalid value of logscale'
          write(*,*) 'logscale = 1 for log scale'
          write(*,*) 'logscale = 0 for linear scale'
        endif

c K_1(x)        
        K1_value = bessk1(x)
        K1_approx = dsqrt(pi/(2.d0*x)) * dexp(-x)
     .        *(1.d0 + 3.d0/(8.d0*x) - 15.d0/(128.d0*x**2) + 315.d0/(3072.d0*x**3))
        K1_diff = dabs(K1_value-K1_approx)/K1_value

c K_2(x)
        K2_value = bessk(2,x)
        K2_approx = dsqrt(pi/(2.d0*x)) * dexp(-x)
     .        *(1.d0 + 15.d0/(8.d0*x) + 105.d0/(128.d0*x**2) - 945.d0/(3072.d0*x**3))
        K2_diff = dabs(K2_value-K2_approx)/K2_value

c output to the file
        write(42,fmt='(7(ES18.8))') x, K1_value, K1_approx, K1_diff, K2_value, K2_approx, K2_diff
      enddo

      close(42)

      return
      end

