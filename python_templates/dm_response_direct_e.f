c-------------------------------------------------------------------------c
      function dm_response()
c-------------------------------------------------------------------------c
c                                                                         c
c     Compute the dark matter response for the DM-electron scattering,    c
c                                                                         c
c     Please cite: arXiv:1912.08204v2 [hep-ph] 19 Aug 2020                c
c                                                                         c
c-------------------------------------------------------------------------c
            
            integer dm_spin
            double precision get_dm_response

            include 'maddm.inc'
            include 'coupl.inc'
            include 'process_names.inc'
            include 'maddm_card.inc'
            include 'dm_info.inc'

c           Initialize the parameters
            dm_spin = dof_dm(1)     ! spin of the dm particle


c --------------------------------------------------------------------------------
c Select scalar, fermion or vector DM particle and compute the dark matter
c response function
c --------------------------------------------------------------------------------
            select case(dm_spin)
                  case(1) ! will be implmented in future work
                  write (*,*) "scalar DM-e scattering still not implemented!"
                  case(2)
                  dm_response = get_dm_response()
                  case(3) ! will be implmented in future work
                  write (*,*) "vector DM-e scattering still not implemented!"
            end select



      end function
            



            
c-------------------------------------------------------------------------c
      function get_dm_response()
c-------------------------------------------------------------------------c
c
c Compute the non relativistic coefficients and use them to get the       c
c dark matter response function.                                          c
c                                                                         c
c-------------------------------------------------------------------------c
            implicit none
            
            integer j, k
            integer j_dm_e, j_antidm_e
            integer j_eff_SI_dm_e, j_tot_SI_dm_e, j_eff_SI_dm_p, j_tot_SI_dm_p
            integer j_eff_SD_dm_e, j_tot_SD_dm_e, j_eff_SD_dm_p, j_tot_SD_dm_p
            double precision p_ext(0:3,4)
            double precision M_e, M_dm
            double precision Minterf_SI_e, Minterf_SI_p, Minterf_SD_e, Minterf_SD_p
            double precision alpha_SI_even, alpha_SI_odd, alpha_SD_even, alpha_SD_odd, c_1, c_4
            double precision get_dm_response
            character(2) siorsd
            integer len_process_name
            character(50) process_name , dm_e_process_name, dm_p_process_name

            include 'maddm.inc'
            include 'coupl.inc'
            include 'process_names.inc'
            include 'maddm_card.inc'

            %(electron_mass)s       ! electron mass M_e (GeV)
            M_dm = mdm(1)           ! dark matter mass (GeV)

            dm_e_process_name = trim(DM_NAMES(1)) // 'em_' // trim(DM_NAMES(1)) // 'em'
            dm_p_process_name = trim(DM_NAMES(1)) // 'ep_' // trim(DM_NAMES(1)) // 'ep'
c            write (*,*) 'dm_e_process_name = ' , dm_e_process_name
c            write (*,*) 'dm_p_process_name = ' , dm_p_process_name
c            write(*,*) 'dm info : ' , dof_dm(1)
c            write(*,*) 'number of processes: ', dd_num_processes

c           Four momenta of the particles involved in the scattering
            p_ext(0,1) = M_dm
            p_ext(1,1) = 0.d0
            p_ext(2,1) = 0.d0
            p_ext(3,1) = 0.d0

            p_ext(0,2) = M_e
            p_ext(1,2) = 0.d0
            p_ext(2,2) = 0.d0
            p_ext(3,2) = 0.d0
                  
            p_ext(0,3) = M_dm
            p_ext(1,3) = 0.d0
            p_ext(2,3) = 0.d0
            p_ext(3,3) = 0.d0

            p_ext(0,4) = M_e
            p_ext(1,4) = 0.d0
            p_ext(2,4) = 0.d0
            p_ext(3,4) = 0.d0

c--------------------------------------------------------------------------------
c     Find the indices of the dm e- scatt and dm e+ scatt for the full process,
c     the effective and the full + effective
c--------------------------------------------------------------------------------
            do j=1,(dd_num_processes) 

c                 j is the index that identifies the matrix element of the FULL Lagrangian
                  if (DD_PROCESS_NAMES(j).eq.dm_e_process_name) then
                        j_dm_e = j
                  endif
            
                  if (DD_PROCESS_NAMES(j).eq.dm_p_process_name) then
                        j_antidm_e = j
                  endif
            enddo
            

            do k=1, 2*dd_num_processes

c                 j_eff is the index of the matrix_element in the EFF case for the
c                 SI/SD case, and dark matter - electorn/proton process
c                 JTOT same as JEFF but for the EFF+FULL

                  process_name = 'EFT_SI_' // dm_e_process_name
                  if (DD_EFF_PROCESS_NAMES(k).eq.process_name) then
                        j_eff_SI_dm_e = k
                  endif

                  process_name = 'EFT_SI_' // dm_p_process_name
                  if (DD_EFF_PROCESS_NAMES(k).eq.process_name) then
                        j_eff_SI_dm_p = k
                  endif

                  process_name = 'TOT_SI_' // dm_e_process_name
                  if (DD_TOT_PROCESS_NAMES(k).eq.process_name) then
                        j_tot_SI_dm_e = k
                  endif

                  process_name = 'TOT_SI_' // dm_p_process_name
                  if (DD_TOT_PROCESS_NAMES(k).eq.process_name) then
                        j_tot_SI_dm_p = k
                  endif

                  process_name = 'EFT_SD_' // dm_e_process_name
                  if (DD_EFF_PROCESS_NAMES(k).eq.process_name) then
                        j_eff_SD_dm_e = k
                  endif

                  process_name = 'EFT_SD_' // dm_p_process_name
                  if (DD_EFF_PROCESS_NAMES(k).eq.process_name) then
                        j_eff_SD_dm_p = k
                  endif

                  process_name = 'TOT_SD_' // dm_e_process_name
                  if (DD_TOT_PROCESS_NAMES(k).eq.process_name) then
                        j_tot_SD_dm_e = k
                  endif

                  process_name = 'TOT_SD_' // dm_p_process_name
                  if (DD_TOT_PROCESS_NAMES(k).eq.process_name) then
                        j_tot_SD_dm_p = k
                  endif
            enddo

c            write (*,*) "j_dm_e        :" , j_dm_e
c            write (*,*) "j_eff_SI_dm_e :" , j_eff_SI_dm_e
c            write (*,*) "j_tot_SI_dm_e :" , j_tot_SI_dm_e

c            write (*,*) "j_antidm_e    :" , j_antidm_e
c            write (*,*) "j_eff_SI_dm_p :" , j_eff_SI_dm_p
c            write (*,*) "j_tot_SI_dm_p :" , j_tot_SI_dm_p

c            write (*,*) "j_eff_SD_dm_e :" , j_eff_SD_dm_e
c            write (*,*) "j_tot_SD_dm_e :" , j_tot_SD_dm_e
c            write (*,*) "j_eff_SD_dm_p :" , j_eff_SD_dm_p
c            write (*,*) "j_tot_SD_dm_p :" , j_tot_SD_dm_p
            

c--------------------------------------------------------------------------------
c Compute the interference terms and use them to get the even and odd effective 
c coefficients in the SI and SD case
c--------------------------------------------------------------------------------

            Minterf_SI_e = smatrix_dd_tot(p_ext, 1,1,j_tot_SI_dm_e) - max(0d0, smatrix_dd(p_ext,1,1,j_dm_e)) 
     &              - smatrix_dd_eff(p_ext,1,1,j_eff_SI_dm_e)
            Minterf_SI_e = 0.5d0*Minterf_SI_e
            Minterf_SI_e = Minterf_SI_e / smatrix_dd_eff(p_ext,1,1,j_eff_SI_dm_e)

            Minterf_SI_p = smatrix_dd_tot(p_ext, 1,1,j_tot_SI_dm_p) - max(0d0, smatrix_dd(p_ext,1,1,j_antidm_e)) 
     &                   - smatrix_dd_eff(p_ext,1,1,j_eff_SI_dm_p)
            Minterf_SI_p = 0.5d0*Minterf_SI_p
            Minterf_SI_p = Minterf_SI_p / smatrix_dd_eff(p_ext,1,1,j_eff_SI_dm_p)

            Minterf_SD_e = smatrix_dd_tot(p_ext, 1,1,j_tot_SD_dm_e) - max(0d0, smatrix_dd(p_ext,1,1,j_dm_e)) 
     &                   - smatrix_dd_eff(p_ext,1,1,j_eff_SD_dm_e)
            Minterf_SD_e = 0.5d0*Minterf_SD_e
            Minterf_SD_e = Minterf_SD_e / smatrix_dd_eff(p_ext,1,1,j_eff_SD_dm_e)
            
            Minterf_SD_p = smatrix_dd_tot(p_ext, 1,1,j_tot_SD_dm_p) - max(0d0, smatrix_dd(p_ext,1,1,j_antidm_e)) 
     &                   - smatrix_dd_eff(p_ext,1,1,j_eff_SD_dm_p)
            Minterf_SD_p = 0.5d0*Minterf_SD_p
            Minterf_SD_p = Minterf_SD_p / smatrix_dd_eff(p_ext,1,1,j_eff_SD_dm_p)

            
            alpha_SI_even = 0.5d0 * (Minterf_SI_e + Minterf_SI_p)
            alpha_SI_odd  = 0.5d0 * (Minterf_SI_e - Minterf_SI_p)

            alpha_SD_even = 0.5d0 * (Minterf_SD_e + Minterf_SD_p)
            alpha_SD_odd  = 0.5d0 * (Minterf_SD_e - Minterf_SD_p)

c--------------------------------------------------------------------------------
c Compute the non relativistic coefficients and use them to get the
c dark matter response function.
c--------------------------------------------------------------------------------
            c_1 = 4*mdm(1)*M_e*(alpha_SI_odd + alpha_SI_even)
            c_4 = 16*mdm(1)*M_e*(2*alpha_SD_odd - alpha_SD_even)
            get_dm_response = (c_1)**2 + (3.d0/16)*(c_4)**2
            
            if(ISNAN(get_dm_response)) then
                  get_dm_response = -1
            endif

c            write(*,*) '---------------------------------------------------------------'
c            write(*,*) '                    dm_response_direct_e.f                     '
c            write(*,*) '---------------------------------------------------------------'
c          write(*,*) 'Minterf_SI_e   : ',Minterf_SI_e
c          write(*,*) 'Minterf_SI_p   : ',Minterf_SI_p
c          write(*,*) 'Minterf_SD_e   : ',Minterf_SD_e
c           write(*,*) 'Minterf_SD_p   : ',Minterf_SD_p
c            write(*,*) 'alpha_SI_even  : ', alpha_SI_even
c            write(*,*) 'alpha_SI_odd   : ', alpha_SI_odd
c            write(*,*) 'alpha_SD_even  : ', alpha_SD_even
c            write(*,*) 'alpha_SD_odd   : ', alpha_SD_odd
c            write(*,*) 'c_1            : ', c_1
c            write(*,*) 'c_4            : ', c_4
c            write(*,*) 'dm_response    : ', get_dm_response
c            write(*,*) '---------------------------------------------------------------'
      
      end function
