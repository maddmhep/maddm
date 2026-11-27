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
            logical do_get_alpha, found
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
                  call replace_tilde(dm_t_process_name)
                  call replace_tilde(dm_tx_process_name)
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
     &                                     m_target, aSIe, aSIo, aSDe, aSDo, found)
                  if (.not.found) return
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
     &                               alpha_SI_even, alpha_SI_odd, alpha_SD_even, alpha_SD_odd, found)
c-------------------------------------------------------------------------c
c                                                                         c
c Compute the even and odd effective coefficients (alphas) from the       c
c interference between the full and the effective matrix elements, for   c
c the DM-target (proc_name) and DM-antitarget (antiproc_name) processes   c
c at rest. The SD alphas are computed only if dm_spin.ne.1.               c
c found is .false. (and the alphas are not computed) if one of the full,  c
c effective or full+effective processes is missing.                       c
c                                                                         c
c-------------------------------------------------------------------------c
            use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
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
            logical found
            character(6) op_SI, op_SD
            double precision interference_ratio

            include 'maddm.inc'
            include 'coupl.inc'
            include 'process_names.inc'
            include 'maddm_card.inc'
            include 'dm_info.inc'

            found = .false.
            j_dm_t = 0
            j_dm_tx = 0
            j_eff_SI_dm_t = 0
            j_eff_SI_dm_tx = 0
            j_tot_SI_dm_t = 0
            j_tot_SI_dm_tx = 0
            j_eff_SD_dm_t = 0
            j_eff_SD_dm_tx = 0
            j_tot_SD_dm_t = 0
            j_tot_SD_dm_tx = 0

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

c           j_eff / j_tot are the indices of the matrix element in the EFF and EFF+FULL
c           cases for the SI/SD dark matter - target/antitarget processes. Each array is
c           scanned up to its own size (nb_me_dd_eff, nb_me_dd_tot in maddm.inc).
            do k=1, size(DD_EFF_PROCESS_NAMES)
                  process_name = 'EFT_SI_' // proc_name
                  if (DD_EFF_PROCESS_NAMES(k).eq.process_name) j_eff_SI_dm_t = k

                  process_name = 'EFT_SI_' // antiproc_name
                  if (DD_EFF_PROCESS_NAMES(k).eq.process_name) j_eff_SI_dm_tx = k

                  if(dm_spin.ne.1) then
                        process_name = 'EFT_SD_' // proc_name
                        if (DD_EFF_PROCESS_NAMES(k).eq.process_name) j_eff_SD_dm_t = k

                        process_name = 'EFT_SD_' // antiproc_name
                        if (DD_EFF_PROCESS_NAMES(k).eq.process_name) j_eff_SD_dm_tx = k
                  endif
            enddo

            do k=1, size(DD_TOT_PROCESS_NAMES)
                  process_name = 'TOT_SI_' // proc_name
                  if (DD_TOT_PROCESS_NAMES(k).eq.process_name) j_tot_SI_dm_t = k

                  process_name = 'TOT_SI_' // antiproc_name
                  if (DD_TOT_PROCESS_NAMES(k).eq.process_name) j_tot_SI_dm_tx = k

                  if(dm_spin.ne.1) then
                        process_name = 'TOT_SD_' // proc_name
                        if (DD_TOT_PROCESS_NAMES(k).eq.process_name) j_tot_SD_dm_t = k

                        process_name = 'TOT_SD_' // antiproc_name
                        if (DD_TOT_PROCESS_NAMES(k).eq.process_name) j_tot_SD_dm_tx = k
                  endif
            enddo

c           If one of the full, effective or full+effective processes is missing
c           (e.g. no effective operator for this DM type) the alphas cannot be computed
            if ((j_dm_t*j_dm_tx*j_eff_SI_dm_t*j_eff_SI_dm_tx*j_tot_SI_dm_t*j_tot_SI_dm_tx).eq.0) return
            if ((dm_spin.ne.1).and.((j_eff_SD_dm_t*j_eff_SD_dm_tx*j_tot_SD_dm_t*j_tot_SD_dm_tx).eq.0)) return

c--------------------------------------------------------------------------------
c           Compute the interference terms and use them to get the even and odd effective
c           coefficients in the SI and SD case
c--------------------------------------------------------------------------------
            if (dm_spin.eq.1) then
                  op_SI = 'SIEFFS'
            else
                  op_SI = 'SIEFFF'
            endif
            op_SD = 'SDEFFF'

            Minterf_SI_t = interference_ratio(p_ext, j_dm_t, j_eff_SI_dm_t, j_tot_SI_dm_t, op_SI)
            Minterf_SI_tx = interference_ratio(p_ext, j_dm_tx, j_eff_SI_dm_tx, j_tot_SI_dm_tx, op_SI)

            alpha_SI_even = 0.5d0 * (Minterf_SI_t + Minterf_SI_tx)
            alpha_SI_odd  = 0.5d0 * (Minterf_SI_t - Minterf_SI_tx)

            if(dm_spin.ne.1) then
                  Minterf_SD_t = interference_ratio(p_ext, j_dm_t, j_eff_SD_dm_t, j_tot_SD_dm_t, op_SD)
                  Minterf_SD_tx = interference_ratio(p_ext, j_dm_tx, j_eff_SD_dm_tx, j_tot_SD_dm_tx, op_SD)

                  alpha_SD_even = 0.5d0 * (Minterf_SD_t + Minterf_SD_tx)
                  alpha_SD_odd  = 0.5d0 * (Minterf_SD_t - Minterf_SD_tx)
            endif

c           found stays .false. if the alphas are not finite (e.g. smatrix_dd_eff = 0)
            if (.not.all(ieee_is_finite((/alpha_SI_even, alpha_SI_odd, alpha_SD_even, alpha_SD_odd/)))) return
            found = .true.

      end subroutine get_dm_target_alpha


c-------------------------------------------------------------------------c
      double precision function interference_ratio(p_ext, j_full, j_eff, j_tot, op)
c-------------------------------------------------------------------------c
c                                                                         c
c Re(M_full M_eff*) / |M_eff|^2 for the effective operator op (coupling   c
c order SIEFFS, SIEFFF or SDEFFF) at unit effective coupling, i.e. the    c
c effective coefficient (alpha) of the model for that operator.           c
c                                                                         c
c It is extracted from |M_tot|^2 - |M_full|^2 - |M_eff|^2. The effective  c
c couplings are first rescaled by c = 2**kexp (couplings of order k by    c
c c**(k/2), so that M_eff -> c M_eff), with c chosen so that |M_eff| is   c
c of the order of |M_full|: otherwise, for a model coefficient much       c
c smaller than the unit effective coupling (e.g. a feeble or kinetically  c
c mixed mediator), the interference is lost in the double precision       c
c cancellation and alpha comes out as zero. kexp is even, so the          c
c rescaling and its inverse are exact powers of two.                      c
c                                                                         c
c-------------------------------------------------------------------------c
            implicit none

            double precision p_ext(0:3,4)
            integer j_full, j_eff, j_tot
            character(*) op
            double precision M_full, M_eff1, M_eff, M_tot, M_interf
            integer kexp
            logical scaled, restored

            include 'maddm.inc'

            M_full = max(0d0, smatrix_dd(p_ext,1,1,j_full))
            M_eff1 = smatrix_dd_eff(p_ext,1,1,j_eff)

            kexp = 0
            if ((M_full.gt.0d0).and.(M_eff1.gt.0d0)) then
c                 |M_eff(c)/M_full| ~ 1 for c = sqrt(M_full/M_eff1), rounded to 2**kexp
                  kexp = 2*nint(log(M_full/M_eff1)/(4d0*log(2d0)))
                  kexp = max(-1000, min(1000, kexp))
            endif

            scaled = .false.
            if (kexp.ne.0) call scale_eff_couplings(op, kexp, scaled)
            if (.not.scaled) kexp = 0

            if (scaled) then
                  M_eff = smatrix_dd_eff(p_ext,1,1,j_eff)
            else
                  M_eff = M_eff1
            endif
            M_tot = smatrix_dd_tot(p_ext,1,1,j_tot)

            if (scaled) call scale_eff_couplings(op, -kexp, restored)

c           An interference below the rounding error of the subtraction is zero (e.g. the SD
c           operator for a vector mediator): it would otherwise give a ~1e-32 relative noise
            M_interf = M_tot - M_full - M_eff
            if (abs(M_interf).le.1d3*epsilon(1d0)*(abs(M_tot) + M_full + abs(M_eff))) M_interf = 0d0

            interference_ratio = 0.5d0*M_interf / M_eff * 2d0**kexp

      end function interference_ratio


c-------------------------------------------------------------------------c
      subroutine scale_eff_couplings(op, kexp, scaled)
c-------------------------------------------------------------------------c
c                                                                         c
c Multiply the couplings of the effective operator op by 2**(kexp*k/2),   c
c k being the coupling order (written by MGoutput.py from the model       c
c interactions). scaled is .false. if op has no couplings that can be     c
c rescaled alone (missing, or shared with other interactions).            c
c                                                                         c
c-------------------------------------------------------------------------c
            implicit none

            character(*) op
            integer kexp
            logical scaled

            include 'coupl.inc'

            scaled = .false.
%(eff_coupling_scaling)s

      end subroutine scale_eff_couplings



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
