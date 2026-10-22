      subroutine dmquark_alpha(quark, qalphas)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                            c
c compute the DM-quark alphas (aSIe, aSIo, aSDe, aSDo)       c
c quark = 1..6 for d, u, s, c, b, t                          c
c                                                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            implicit none

            integer quark
            double precision qalphas(4)  ! Return array aSIe, aSIo, aSDe, aSDo
            double precision m_quark

            call dm_alpha('q', quark, qalphas, m_quark)

      end subroutine dmquark_alpha



      subroutine dmelectron_alpha(ealphas, m_e)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                            c
c compute the DM-electron alphas (aSIe, aSIo, aSDe, aSDo)    c
c and return the electron mass used (GeV)                    c
c                                                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            implicit none

            double precision ealphas(4)  ! Return array aSIe, aSIo, aSDe, aSDo
            double precision m_e

            call dm_alpha('e', 0, ealphas, m_e)

      end subroutine dmelectron_alpha



c-------------------------------------------------------------------------c
      subroutine dm_alpha(target_type, idx, alphas, m_target)
c-------------------------------------------------------------------------c
c                                                                         c
c Switch between the DM-quark (target_type='q', idx = quark 1..6) and     c
c the DM-electron (target_type='e', idx unused) alphas.                   c
c If the process does not exist all the alphas are set to -1.             c
c For scalar DM (dm_spin=1) the SD alphas are set to 0.                   c
c                                                                         c
c-------------------------------------------------------------------------c
            implicit none

            character(1) target_type
            integer idx, dm_spin, j
            logical do_get_alpha
            double precision alphas(4)  ! Return array aSIe, aSIo, aSDe, aSDo
            double precision m_target
            double precision aSIe, aSIo, aSDe, aSDo
            double precision M(6), M_e
            character(len=1), dimension(6) :: quark_names
            character(50) dm_t_process_name, dm_tx_process_name
            character(50) initial_process_name, initial_antiprocess_name

            include 'maddm.inc'
            include 'coupl.inc'
            include 'process_names.inc'
            include 'maddm_card.inc'
            include 'dm_info.inc'

            data quark_names /"d", "u", "s", "c", "b", "t"/

c initialize the quark masses array (d,u,s,c,b,t) and the electron mass
c--------------------------------------------------------------------------------

            %(quark_masses)s
            %(electron_mass)s

            dm_spin = dof_dm(1)     ! spin of the dm particle

c           Process names: DM-q / DM-qx for quarks (with ~ replaced by x),
c           DM-e- / DM-e+ for electrons
            if (target_type.eq.'e') then
                  dm_t_process_name = trim(DM_NAMES(1)) // 'em_' // trim(DM_NAMES(1)) // 'em'
                  dm_tx_process_name = trim(DM_NAMES(1)) // 'ep_' // trim(DM_NAMES(1)) // 'ep'
                  m_target = M_e
            else
                  initial_process_name = trim(DM_NAMES(1)) // trim(quark_names(idx))
                  initial_antiprocess_name = trim(DM_NAMES(1)) // trim(quark_names(idx)) // "x"
                  call replace_tilde(initial_process_name)
                  call replace_tilde(initial_antiprocess_name)
                  dm_t_process_name = trim(initial_process_name) // "_" // trim(initial_process_name)
                  dm_tx_process_name = trim(initial_antiprocess_name) // "_" // trim(initial_antiprocess_name)
                  m_target = M(idx)
            endif

c           Missing process -> alphas = -1
            do j=1,4
                  alphas(j) = -1d0
            enddo

            do_get_alpha = .false.
c           Check if there are DM-target/DM-antitarget processes
            do j=1,(dd_num_processes)
c                 j is the index that identifies the matrix element of the FULL Lagrangian
                  if ((DD_PROCESS_NAMES(j).eq.dm_t_process_name).or.
     &                (DD_PROCESS_NAMES(j).eq.dm_tx_process_name)) then
                        do_get_alpha = .true.
                  endif
            enddo

            if (.not.do_get_alpha) return

c --------------------------------------------------------------------------------
c Select scalar, fermion or vector DM particle and compute the alphas
c --------------------------------------------------------------------------------
            select case(dm_spin)
                  case(1,2)
                  aSDe = 0d0
                  aSDo = 0d0
                  call get_dm_target_alpha(dm_spin, dm_t_process_name, dm_tx_process_name,
     &                                     m_target, aSIe, aSIo, aSDe, aSDo)
                  case(3) ! will be implmented in future work
                  write (*,*) "vector DM-target scattering still not implemented!"
                  return
                  case default
                  return
            end select

            if (ISNAN(aSIe).or.ISNAN(aSIo).or.ISNAN(aSDe).or.ISNAN(aSDo)) return

            alphas(1) = aSIe
            alphas(2) = aSIo
            alphas(3) = aSDe
            alphas(4) = aSDo

      end subroutine dm_alpha


c-------------------------------------------------------------------------c
      subroutine get_dm_target_alpha(dm_spin, proc_name, antiproc_name, m_target,
     &                               alpha_SI_even, alpha_SI_odd, alpha_SD_even, alpha_SD_odd)
c-------------------------------------------------------------------------c
c                                                                         c
c Compute the even and odd effective coefficients (alphas) from the       c
c interference between the full and the effective matrix elements, for   c
c the DM-target (proc_name) and DM-antitarget (antiproc_name) processes   c
c at rest. The SD alphas are computed only if dm_spin.ne.1.               c
c                                                                         c
c-------------------------------------------------------------------------c
            implicit none

            integer j, k, dm_spin
            integer j_dm_t, j_dm_tx
            integer j_eff_SI_dm_t, j_tot_SI_dm_t, j_eff_SI_dm_tx, j_tot_SI_dm_tx
            integer j_eff_SD_dm_t, j_tot_SD_dm_t, j_eff_SD_dm_tx, j_tot_SD_dm_tx
            double precision p_ext(0:3,4)
            double precision m_target
            character(50) proc_name, antiproc_name
            character(50) process_name
            double precision Minterf_SI_t, Minterf_SI_tx, Minterf_SD_t, Minterf_SD_tx
            double precision alpha_SI_even, alpha_SI_odd, alpha_SD_even, alpha_SD_odd

            include 'maddm.inc'
            include 'coupl.inc'
            include 'process_names.inc'
            include 'maddm_card.inc'
            include 'dm_info.inc'

c Four momenta of the particles involved in the scattering (at rest)
            p_ext(0,1) = mdm(1)
            p_ext(1,1) = 0.d0
            p_ext(2,1) = 0.d0
            p_ext(3,1) = 0.d0

            p_ext(0,2) = m_target
            p_ext(1,2) = 0.d0
            p_ext(2,2) = 0.d0
            p_ext(3,2) = 0.d0

            p_ext(0,3) = mdm(1)
            p_ext(1,3) = 0.d0
            p_ext(2,3) = 0.d0
            p_ext(3,3) = 0.d0

            p_ext(0,4) = m_target
            p_ext(1,4) = 0.d0
            p_ext(2,4) = 0.d0
            p_ext(3,4) = 0.d0

c--------------------------------------------------------------------------------
c     Find the indices of the dm-target and dm-antitarget scatt for the full process,
c     the effective and the full + effective
c--------------------------------------------------------------------------------
            do j=1,(dd_num_processes)
c    j is the index that identifies the matrix element of the FULL Lagrangian
                  if (DD_PROCESS_NAMES(j).eq.proc_name) then
                        j_dm_t = j
                  endif

                  if (DD_PROCESS_NAMES(j).eq.antiproc_name) then
                        j_dm_tx = j
                  endif
            enddo

            do k=1, 2*dd_num_processes

c                 j_eff is the index of the matrix_element in the EFF case for the
c                 SI/SD case, and dark matter - target/antitarget process.
c                 JTOT same as JEFF but for the EFF+FULL

                  process_name = 'EFT_SI_' // proc_name
                  if (DD_EFF_PROCESS_NAMES(k).eq.process_name) then
                        j_eff_SI_dm_t = k
                  endif

                  process_name = 'EFT_SI_' // antiproc_name
                  if (DD_EFF_PROCESS_NAMES(k).eq.process_name) then
                        j_eff_SI_dm_tx = k
                  endif

                  process_name = 'TOT_SI_' // proc_name
                  if (DD_TOT_PROCESS_NAMES(k).eq.process_name) then
                        j_tot_SI_dm_t = k
                  endif

                  process_name = 'TOT_SI_' // antiproc_name
                  if (DD_TOT_PROCESS_NAMES(k).eq.process_name) then
                        j_tot_SI_dm_tx = k
                  endif

                  if(dm_spin.ne.1) then
                        process_name = 'EFT_SD_' // proc_name
                        if (DD_EFF_PROCESS_NAMES(k).eq.process_name) then
                              j_eff_SD_dm_t = k
                        endif

                        process_name = 'EFT_SD_' // antiproc_name
                        if (DD_EFF_PROCESS_NAMES(k).eq.process_name) then
                              j_eff_SD_dm_tx = k
                        endif

                        process_name = 'TOT_SD_' // proc_name
                        if (DD_TOT_PROCESS_NAMES(k).eq.process_name) then
                              j_tot_SD_dm_t = k
                        endif

                        process_name = 'TOT_SD_' // antiproc_name
                        if (DD_TOT_PROCESS_NAMES(k).eq.process_name) then
                              j_tot_SD_dm_tx = k
                        endif
                  endif
            enddo

c--------------------------------------------------------------------------------
c           Compute the interference terms and use them to get the even and odd effective
c           coefficients in the SI and SD case
c--------------------------------------------------------------------------------
            Minterf_SI_t = smatrix_dd_tot(p_ext, 1,1,j_tot_SI_dm_t) - max(0d0, smatrix_dd(p_ext,1,1,j_dm_t))
     &              - smatrix_dd_eff(p_ext,1,1,j_eff_SI_dm_t)
            Minterf_SI_t = 0.5d0*Minterf_SI_t
            Minterf_SI_t = Minterf_SI_t / smatrix_dd_eff(p_ext,1,1,j_eff_SI_dm_t)

            Minterf_SI_tx = smatrix_dd_tot(p_ext, 1,1,j_tot_SI_dm_tx) - max(0d0, smatrix_dd(p_ext,1,1,j_dm_tx))
     &                   - smatrix_dd_eff(p_ext,1,1,j_eff_SI_dm_tx)
            Minterf_SI_tx = 0.5d0*Minterf_SI_tx
            Minterf_SI_tx = Minterf_SI_tx / smatrix_dd_eff(p_ext,1,1,j_eff_SI_dm_tx)

            alpha_SI_even = 0.5d0 * (Minterf_SI_t + Minterf_SI_tx)
            alpha_SI_odd  = 0.5d0 * (Minterf_SI_t - Minterf_SI_tx)

            if(dm_spin.ne.1) then
                  Minterf_SD_t = smatrix_dd_tot(p_ext, 1,1,j_tot_SD_dm_t) - max(0d0, smatrix_dd(p_ext,1,1,j_dm_t))
     &                         - smatrix_dd_eff(p_ext,1,1,j_eff_SD_dm_t)
                  Minterf_SD_t = 0.5d0*Minterf_SD_t
                  Minterf_SD_t = Minterf_SD_t / smatrix_dd_eff(p_ext,1,1,j_eff_SD_dm_t)

                  Minterf_SD_tx = smatrix_dd_tot(p_ext, 1,1,j_tot_SD_dm_tx) - max(0d0, smatrix_dd(p_ext,1,1,j_dm_tx))
     &                         - smatrix_dd_eff(p_ext,1,1,j_eff_SD_dm_tx)
                  Minterf_SD_tx = 0.5d0*Minterf_SD_tx
                  Minterf_SD_tx = Minterf_SD_tx / smatrix_dd_eff(p_ext,1,1,j_eff_SD_dm_tx)

                  alpha_SD_even = 0.5d0 * (Minterf_SD_t + Minterf_SD_tx)
                  alpha_SD_odd  = 0.5d0 * (Minterf_SD_t - Minterf_SD_tx)
            endif

      end subroutine get_dm_target_alpha



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
