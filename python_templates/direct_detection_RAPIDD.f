      subroutine dmquark_alpha(quark, qalphas)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                            c
c compute DM-nucleon cross-section                           c
c if proton=1 do the sigma_p computation otherwise sigma_n   c
c                                                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            implicit none
      
            integer quark, quarks, dm_spin
            integer j
            !double precision get_dm_alpha
            logical do_get_dm_response
            double precision result_value  ! Placeholder for actual return value
            character(len=1), dimension(6) :: quark_names
            character(50) dm_q_process_name, dm_qx_process_name
            character(50) initial_process_name, initial_antiprocess_name
            double precision qalphas(4)  ! Return array aSIe, aSIo, aSDe, aSDo
            double precision aSIe, aSIo, aSDe, aSDo



            include 'maddm.inc'
            include 'coupl.inc'
            include 'process_names.inc'
            include 'maddm_card.inc'
            include 'dm_info.inc'
c           Initialize the parameters
            dm_spin = dof_dm(1)     ! spin of the dm particle

            ! write(*,*) "dm spin:", dm_spin
            data quark_names /"d", "u", "s", "c", "b", "t"/
            ! write(*,*) "DM NAME:", trim(DM_NAMES(1))
            ! write(*,*) "quark number call 1:", quark
            ! write(*,*) "quark name:", trim(quark_names(quark))


            ! dm_q_process_name = trim(DM_NAMES(1)) // trim(quark_names(quark))
            initial_process_name = trim(DM_NAMES(1)) // trim(quark_names(quark)) 
            initial_antiprocess_name = trim(DM_NAMES(1)) // trim(quark_names(quark)) // "x" 
            call replace_tilde(initial_process_name)
            call replace_tilde(initial_antiprocess_name)
            dm_q_process_name = trim(initial_process_name) // "_" // trim(initial_process_name) 
            ! dm_q_process_name = trim(DM_NAMES(1)) // trim(quark_names(quark)) &
            !       & // "_" // trim(DM_NAMES(1)) // trim(quark_names(quark))

            dm_qx_process_name = trim(initial_antiprocess_name) // "_" // trim(initial_antiprocess_name) 

            !write (*, *) "process name :", dm_q_process_name
            !write (*, *) "process name :", dm_qx_process_name


            do_get_dm_response = .false.
c           Check if there are DM-q/DM-qx processes
            do j=1,(dd_num_processes)
c                 j is the index that identifies the matrix element of the FULL Lagrangian
                  if ((DD_PROCESS_NAMES(j).eq.dm_q_process_name).or.(DD_PROCESS_NAMES(j).eq.dm_qx_process_name)) then
                        do_get_dm_response = .true.
                  endif
            enddo

            !write (*, *) "is there a process:" , do_get_dm_response


            ! Placeholder: Set function return value
            ! result_value = 0.0d0  ! Temporary return value (replace with actual logic)
            ! dmquark_alpha = result_value  ! Assign to function name

c --------------------------------------------------------------------------------
c Select scalar, fermion or vector DM particle and compute the dark matter
c response function
c --------------------------------------------------------------------------------
            !dmquark_alpha = -1.0d0
            if (do_get_dm_response) then
                  select case(dm_spin)
                        case(1)
                        call get_dmquark_alpha(dm_spin, quark,aSIe, aSIo, aSDe, aSDo)
                        case(2)
                        call get_dmquark_alpha(dm_spin, quark, aSIe, aSIo, aSDe, aSDo)
                        case(3) ! will be implmented in future work
                        write (*,*) "vector DM-e scattering still not implemented!"
                  end select
            endif

            qalphas(1) = aSIe
            qalphas(2) = aSIo
            qalphas(3) = aSDe
            qalphas(4) = aSDo

            !write (*, *) "alpha SI even :", aSIe
            !write (*, *) "alpha SI even :", qalphas(1)

            !write (*, *) "get dm :", dmquark_alpha


      end subroutine dmquark_alpha


c-------------------------------------------------------------------------c
      subroutine get_dmquark_alpha(dm_spin, quark, alpha_SI_even, alpha_SI_odd, alpha_SD_even, alpha_SD_odd)
c-------------------------------------------------------------------------c
c
c Compute the non relativistic coefficients and use them to get the       c
c dark matter response function.                                          c
c                                                                         c
c-------------------------------------------------------------------------c
            implicit none
                        
            integer j, k, dm_spin, quark
            integer j_dm_q, j_dm_qx
            integer j_eff_SI_dm_q, j_tot_SI_dm_q, j_eff_SI_dm_qx, j_tot_SI_dm_qx
            integer j_eff_SD_dm_q, j_tot_SD_dm_q, j_eff_SD_dm_qx, j_tot_SD_dm_qx
            double precision p_ext(0:3,4)
            double precision M(6)
            character(len=1), dimension(6) :: quark_names
            character(50) dm_q_process_name, dm_qx_process_name
            character(50) initial_process_name, initial_antiprocess_name
            character(50) process_name
            double precision Minterf_SI_q, Minterf_SI_qx, Minterf_SD_q, Minterf_SD_qx
            double precision alpha_SI_even, alpha_SI_odd, alpha_SD_even, alpha_SD_odd, c_1, c_4

            include 'maddm.inc'
            include 'coupl.inc'
            include 'process_names.inc'
            include 'maddm_card.inc'
            include 'dm_info.inc'

c initialize the quark masses array (u,d,s,c,b,t)
c--------------------------------------------------------------------------------

            %(quark_masses)s

            !get_dmquark_alpha = 0.0d0


            !write (*,*) "here "

            ! write (*,*) "dm mass", mdm(1)
            !write (*,*) "quark number:", quark
            ! write (*,*) "q mass", M(quark)

            data quark_names /"d", "u", "s", "c", "b", "t"/


            dm_q_process_name = trim(DM_NAMES(1)) // trim(quark_names(quark))
            initial_process_name = trim(DM_NAMES(1)) // trim(quark_names(quark)) 
            initial_antiprocess_name = trim(DM_NAMES(1)) // trim(quark_names(quark)) // "x" 
            call replace_tilde(initial_process_name)
            call replace_tilde(initial_antiprocess_name)
            dm_q_process_name = trim(initial_process_name) // "_" // trim(initial_process_name) 
            dm_qx_process_name = trim(initial_antiprocess_name) // "_" // trim(initial_antiprocess_name) 

            ! write (*, *) "process name :", dm_q_process_name
            ! write (*, *) "process name :", dm_qx_process_name

c HERE DEFINE THE FOUR MOMENTA - RIGHT NOW JUST ANY NUMBERS
            p_ext(0,1) = mdm(1)
            p_ext(1,1) = 0.d0
            p_ext(2,1) = 0.d0
            p_ext(3,1) = 0.d0
            
            p_ext(0,2) = M(quark)
            p_ext(1,2) = 0.d0
            p_ext(2,2) = 0.d0
            p_ext(3,2) = 0.d0
                  
            p_ext(0,3) = mdm(1)
            p_ext(1,3) = 0.d0
            p_ext(2,3) = 0.d0
            p_ext(3,3) = 0.d0

            p_ext(0,4) = M(quark)
            p_ext(1,4) = 0.d0
            p_ext(2,4) = 0.d0
            p_ext(3,4) = 0.d0          

c--------------------------------------------------------------------------------
c     Find the indices of the dm e- scatt and dm e+ scatt for the full process,
c     the effective and the full + effective
c--------------------------------------------------------------------------------
            do j=1,(dd_num_processes) 
c    j is the index that identifies the matrix element of the FULL Lagrangian
                        if (DD_PROCESS_NAMES(j).eq.dm_q_process_name) then
                              j_dm_q = j
                        endif
                  
                        if (DD_PROCESS_NAMES(j).eq.dm_qx_process_name) then
                              j_dm_qx = j
                        endif
                  enddo
                  
      
                  do k=1, 2*dd_num_processes
      
c                 j_eff is the index of the matrix_element in the EFF case for the
c                 SI/SD case, and dark matter - electorn/proton process.
c                 JTOT same as JEFF but for the EFF+FULL
      
                        process_name = 'EFT_SI_' // dm_q_process_name
                        if (DD_EFF_PROCESS_NAMES(k).eq.process_name) then
                              j_eff_SI_dm_q = k
                        endif
      
                        process_name = 'EFT_SI_' // dm_qx_process_name
                        if (DD_EFF_PROCESS_NAMES(k).eq.process_name) then
                              j_eff_SI_dm_qx = k
                        endif
      
                        process_name = 'TOT_SI_' // dm_q_process_name
                        if (DD_TOT_PROCESS_NAMES(k).eq.process_name) then
                              j_tot_SI_dm_q = k
                        endif
      
                        process_name = 'TOT_SI_' // dm_qx_process_name
                        if (DD_TOT_PROCESS_NAMES(k).eq.process_name) then
                              j_tot_SI_dm_qx = k
                        endif
      
                        if(dm_spin.ne.1) then
                              process_name = 'EFT_SD_' // dm_q_process_name
                              if (DD_EFF_PROCESS_NAMES(k).eq.process_name) then
                                    j_eff_SD_dm_q = k
                              endif
      
                              process_name = 'EFT_SD_' // dm_qx_process_name
                              if (DD_EFF_PROCESS_NAMES(k).eq.process_name) then
                                    j_eff_SD_dm_qx = k
                              endif
      
                              process_name = 'TOT_SD_' // dm_q_process_name
                              if (DD_TOT_PROCESS_NAMES(k).eq.process_name) then
                                    j_tot_SD_dm_q = k
                              endif
      
                              process_name = 'TOT_SD_' // dm_qx_process_name
                              if (DD_TOT_PROCESS_NAMES(k).eq.process_name) then
                                    j_tot_SD_dm_qx = k
                              endif
                        endif
                  enddo

            ! write (*,*) "j_dm_q        :" , j_dm_q
            ! write (*,*) "j_eff_SI_dm_q :" , j_eff_SI_dm_q
            ! write (*,*) "j_tot_SI_dm_q :" , j_tot_SI_dm_q

            ! write (*,*) "j_dm_qx        :" , j_dm_qx
            ! write (*,*) "j_eff_SI_dm_qx :" , j_eff_SI_dm_qx
            ! write (*,*) "j_tot_SI_dm_qx :" , j_tot_SI_dm_qx

            ! write (*,*) "j_eff_SD_dm_q :" , j_eff_SD_dm_q
            ! write (*,*) "j_tot_SD_dm_q :" , j_tot_SD_dm_q
            ! write (*,*) "j_eff_SD_dm_qx :" , j_eff_SD_dm_qx
            ! write (*,*) "j_tot_SD_dm_qx :" , j_tot_SD_dm_qx
 

            


            Minterf_SI_q = smatrix_dd_tot(p_ext, 1,1,j_tot_SI_dm_q) - max(0d0, smatrix_dd(p_ext,1,1,j_dm_q)) 
     &              - smatrix_dd_eff(p_ext,1,1,j_eff_SI_dm_q)
            Minterf_SI_q = 0.5d0*Minterf_SI_q
            Minterf_SI_q = Minterf_SI_q / smatrix_dd_eff(p_ext,1,1,j_eff_SI_dm_q)

            Minterf_SI_qx = smatrix_dd_tot(p_ext, 1,1,j_tot_SI_dm_qx) - max(0d0, smatrix_dd(p_ext,1,1,j_dm_qx)) 
     &                   - smatrix_dd_eff(p_ext,1,1,j_eff_SI_dm_qx)
            Minterf_SI_qx = 0.5d0*Minterf_SI_qx
            Minterf_SI_qx = Minterf_SI_qx / smatrix_dd_eff(p_ext,1,1,j_eff_SI_dm_qx)

            alpha_SI_even = 0.5d0 * (Minterf_SI_q + Minterf_SI_qx)
            alpha_SI_odd  = 0.5d0 * (Minterf_SI_q - Minterf_SI_qx)

            if(dm_spin.ne.1) then
                  Minterf_SD_q = smatrix_dd_tot(p_ext, 1,1,j_tot_SD_dm_q) - max(0d0, smatrix_dd(p_ext,1,1,j_dm_q)) 
     &                         - smatrix_dd_eff(p_ext,1,1,j_eff_SD_dm_q)
                  Minterf_SD_q = 0.5d0*Minterf_SD_q
                  Minterf_SD_q = Minterf_SD_q / smatrix_dd_eff(p_ext,1,1,j_eff_SD_dm_q)
                  
                  Minterf_SD_qx = smatrix_dd_tot(p_ext, 1,1,j_tot_SD_dm_qx) - max(0d0, smatrix_dd(p_ext,1,1,j_dm_qx)) 
     &                         - smatrix_dd_eff(p_ext,1,1,j_eff_SD_dm_qx)
                  Minterf_SD_qx = 0.5d0*Minterf_SD_qx
                  Minterf_SD_qx = Minterf_SD_qx / smatrix_dd_eff(p_ext,1,1,j_eff_SD_dm_qx)

                  alpha_SD_even = 0.5d0 * (Minterf_SD_q + Minterf_SD_qx)
                  alpha_SD_odd  = 0.5d0 * (Minterf_SD_q - Minterf_SD_qx)
            endif

            
            ! write (*, *) "alpha SI even :", alpha_SI_even
            ! write (*, *) "alpha SI odd:", alpha_SI_odd
            ! write (*, *) "alpha SD even :", alpha_SD_even
            ! write (*, *) "alpha SD odd:", alpha_SD_odd
      end subroutine get_dmquark_alpha
      
      

c -------------------------------------------------------------
c replace the tilde in process name
c -------------------------------------------------------------

      subroutine replace_tilde(input_string)
            implicit none

            character(len=*), intent(inout) :: input_string
            integer :: i, length

            length = len_trim(input_string)

            do i = 1, length
                  if (input_string(i:i) == '~') then
                        input_string(i:i) = 'x'
                  end if
            end do

            end subroutine replace_tilde
