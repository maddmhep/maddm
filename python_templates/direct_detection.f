      function sigma_nucleon(proton,SI)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                            c
c compute DM-nucleon cross-section                           c
c if proton=1 do the sigma_p computation otherwise sigma_n   c
c                                                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      
      integer proton, spin, even, SI
      double precision mr,MN,M(6)
c      double precision ff_N
      complex*16 ff_N, ff_N_odd
      complex*16 NucleonFF
      include 'maddm.inc'
      include 'coupl.inc'
      include 'process_names.inc'
      include 'maddm_card.inc'
      include 'dm_info.inc'

c initialize the quark masses array (u,d,s,c,b,t)
c--------------------------------------------------------------------------------

      %(quark_masses)s

      if (proton .eq. 1) then
         MN=m_neutron
      else
         MN=m_proton
      endif
        
      mr=(MN*mdm(1))/(MN+mdm(1))
      
      spin = dof_dm(1) !get the spin of the dm particle

      if (proton .eq. 1) then
      
      if (SI.eq.1) then
      ff_N=0
      ff_N_odd=0
c--------------------------------------------------------------------------------
c              SPIN INDEPENDENT
c--------------------------------------------------------------------------------
      
c compute nucleon form factor for scalar and vector interactions (so SI interactions)
c--------------------------------------------------------------------------------
      write (*,*) '#################################################################'
      write (*,*) '#                               SI                              #'
      write (*,*) '#################################################################'
      select case(spin) !select scalar, fermion or vector DM particle
        case(1)
          even=1
          ff_N=NucleonFF(M,proton,SI,even) !projection over EVEN operator (scalar interaction)
          ff_N=ff_N*MN
          even=0
          ff_N_odd=NucleonFF(M,proton,SI,even) !projection over ODD operator (vector interaction)
          ff_N=ff_N+(ff_N_odd)
        case(2)
          even=1
          ff_N=NucleonFF(M,proton,SI,even) !projection over EVEN operator (scalar interaction)
          ff_N=(ff_N*MN)
c          even=0
c          ff_N_odd=NucleonFF(M,proton,SI,even) !projection over ODD operator (vector interaction)
c          ff_N=ff_N+(ff_N_odd)
        case(3)
          even=1
          ff_N=NucleonFF(M,proton,SI,even) !projection over EVEN operator (scalar interaction)
          ff_N=ff_N*MN
          even=0
          ff_N_odd=NucleonFF(M,proton,SI,even) !projection over ODD operator (vector interaction)
          ff_N=ff_N+(ff_N_odd)
      end select
      
c compute SI nucleon-DM cross-section
c--------------------------------------------------------------------------------
      if (dm_names(1) .NE. dm_antinames(1)) then !distinguish REAL/COMPLEX field
            ff_N=ff_N/2.0
      endif
      sigma_nucleon=(4.d0/(pi))*(mr**2)*(ff_N*dconjg(ff_N))

      else
c--------------------------------------------------------------------------------
c              SPIN DEPENDENT
c--------------------------------------------------------------------------------
      write (*,*) '#################################################################'
      write (*,*) '#                               SD                              #'
      write (*,*) '#################################################################'
      ff_N=0
      ff_N_odd=0
      select case(spin) !select fermion or vector DM particle
        case(2)
          even=1
          ff_N=NucleonFF(M,proton,SI,even) !projection over EVEN operator (axial-vector interaction)
c          even=0
c          ff_N_odd=NucleonFF(M,proton,SI,even) !projection over ODD operator (tensor interaction)
c          ff_N=ff_N+ff_N_odd
        case(3)
          even=1
          ff_N=NucleonFF(M,proton,SI,even) !projection over EVEN operator (axial-vector interaction)
          even=0
          ff_N_odd=NucleonFF(M,proton,SI,even) !projection over ODD operator (tensor interaction)
          ff_N=(ff_N+ff_N_odd)*2.0
      end select
      if (dm_names(1) .NE. dm_antinames(1)) then !distinguish REAL/COMPLEX field
            ff_N=ff_N/2.0
      endif
      sigma_nucleon=(12.d0/(pi))*(mr**2)*(ff_N*dconjg(ff_N))
      endif

      endif
      
      end function
      








      
      
      
      
      function NucleonFF(M,proton,SI,even)
      
      implicit none
      
      integer i, j, jeff,jtot, k ,proton, even, SI
      integer  j_dm_e, j_eff_dm_e, j_tot_dm_e, j_antidm_e, j_eff_antidm_e, j_tot_antidm_e
      double precision p_ext(0:3,4)
      double precision fNt(6) 
      double precision M(6), M_e
      double precision Minterf_e, Minterf_qb, Minterf, Minterf_dm, Minterf_antidm
      double precision alpha_SI_even, alpha_SI_odd, alpha_SD_even, alpha_SD_odd, c_1, c_4
      complex*16 NucleonFF
      character(2) siorsd
      integer len_process_name
      character(50) process_name , dm_e_process_name, antidm_e_process_name

      include 'maddm.inc'
      include 'coupl.inc'
      include 'process_names.inc'
      include 'maddm_card.inc'
      
      NucleonFF=0.0d0
      M_e = 0.0005109989
      
c--------------------------------------------------------------------------------
c compute the form factors for scalar and vector interactions
c define reduced mass (mr) of the DM-nucleon system
c--------------------------------------------------------------------------------
      if (proton.eq.1) then
        if (SI .eq. 1) then
          select case(even)
    	    case(1) !scalar interaction DONE !
    	      fNt(1)=SPd/M(1)
              fNt(2)=SPu/M(2)
              fNt(3)=SPs/M(3)
              fNt(4)=(2.d0/27.d0)*SPg/M(4)
              fNt(5)=(2.d0/27.d0)*SPg/M(5)
              fNt(6)=(2.d0/27.d0)*SPg/M(6)
            case(0) !vector interaction DONE !
              fNt(1)=VPd
              fNt(2)=VPu
              fNt(3)=0.0d0
              fNt(4)=0.0d0
              fNt(5)=0.0d0
              fNt(6)=0.0d0
          end select
        else
          select case(even)
    	    case(1) !Axial-Vector interaction
    	      fNt(1)=AVPd
              fNt(2)=AVPu
              fNt(3)=AVPs
              fNt(4)=0.0d0
              fNt(5)=0.0d0
              fNt(6)=0.0d0
            case(0) !sigma(mu nu) interaction
              fNt(1)=SigPd
              fNt(2)=SigPu
              fNt(3)=SigPs
              fNt(4)=0.0d0
              fNt(5)=0.0d0
              fNt(6)=0.0d0
          end select
        endif
      else
        if (SI .eq. 1) then
          select case(even)
            case(1) !scalar interaction DONE !
              fNt(1)=SNd/M(1)
              fNt(2)=SNu/M(2)
              fNt(3)=SNs/M(3)
              fNt(4)=(2.0/27)*SNg/M(4)
              fNt(5)=(2.0/27)*SNg/M(5)
              fNt(6)=(2.0/27)*SNg/M(6)
            case(0) !vector interaction DONE !
              fNt(1)=VNd
              fNt(2)=VNu
              fNt(3)=0.0d0
              fNt(4)=0.0d0
              fNt(5)=0.0d0
              fNt(6)=0.0d0
          end select
        else
          select case(even)
    	    case(1) !Axial-Vector interaction
    	      fNt(1)=AVNd
              fNt(2)=AVNu
              fNt(3)=AVNs
              fNt(4)=0.0d0
              fNt(5)=0.0d0
              fNt(6)=0.0d0
            case(0) !sigma(mu nu) interaction
              fNt(1)=SigNd
              fNt(2)=SigNu
              fNt(3)=SigNs
              fNt(4)=0.0d0
              fNt(5)=0.0d0
              fNt(6)=0.0d0
          end select
        endif
      endif

      dm_e_process_name = trim(DM_NAMES(1)) // 'em_' // trim(DM_NAMES(1)) // 'em'
      antidm_e_process_name = trim(DM_NAMES(1)) // 'xem_' // trim(DM_NAMES(1)) // 'xem'
      write (*,*) 'dm_e_process_name = ' , dm_e_process_name
      write (*,*) 'antidm_e_process_name = ' , antidm_e_process_name
      
c      write(*,*) 'dm info : ' , dof_dm(1)
      
c      write(*,*) 'number of processes: ', dd_num_processes
c--------------------------------------------------------------------------------
c loop over the number of processes 
c--------------------------------------------------------------------------------
      do j=1,(dd_num_processes) 
C
C        J is the index for the matrix element in the FULL LAGRANGIAN
C        I is the PDG code of the quark associated to that ME
C        JEFF is the index of the matrix_element in the effective case
C             for the SI/SD case with the same quark associated
C        JTOT same as JEFF but for the EFF+FULL

      i = dd_process_ids(j)
C     Determine JEFF
      if (SI.eq.1)then
         do k=1, 2*dd_num_processes
          siorsd = dd_eff_process_names(k)(5:6) ! return 'SI' or 'SD'
          len_process_name = len_trim(dd_eff_process_names(k))
          process_name = dd_eff_process_names(k)(8:len_process_name)
          if (siorsd.eq.'SI'.and.process_name.eq.DD_PROCESS_NAMES(j))then
             jeff = k
             exit ! exit the loop
          endif
         enddo
      else
         do k=1, 2*dd_num_processes
            siorsd = dd_eff_process_names(k)(5:6) ! return 'SI' or 'SD'
            len_process_name = len_trim(dd_eff_process_names(k))
            process_name = dd_eff_process_names(k)(8:len_process_name)
            if (siorsd.eq.'SD'.and.process_name.eq.DD_PROCESS_NAMES(j))then
               jeff = k
               exit ! exit the loop
            endif
         enddo
      endif
C     SAME for JTOT. FIRST Shortcut test JTOT=JEFF
      if (dd_tot_process_ids(jeff).eq.i.and.dd_tot_process_names(jeff)(5:6).eq.
     &                                       dd_eff_process_names(jeff)(5:6)) then
         jtot = jeff
      else
C     LONG Determination of  JTOT
         if (SI.eq.1)then
            do k=1, 2*dd_num_processes
               siorsd = dd_tot_process_names(k)(5:6) ! return 'SI' or 'SD'
               len_process_name = len_trim(dd_tot_process_names(k))
               process_name = dd_tot_process_names(k)(8:len_process_name)
               if (siorsd.eq.'SI'.and.process_name.eq.DD_PROCESS_NAMES(j))then
                  jtot = k
                  exit          ! exit the loop prevent to loop over undefined part of the array
               endif
            enddo
         else
            do k=1, 2*dd_num_processes
               siorsd = dd_tot_process_names(k)(5:6) ! return 'SI' or 'SD'
               len_process_name = len_trim(dd_tot_process_names(k))
               process_name = dd_tot_process_names(k)(8:len_process_name)
               if (siorsd.eq.'SD'.and.process_name.eq.DD_PROCESS_NAMES(j))then
                  jtot = k
                  exit          ! exit the loop
               endif
            enddo
         endif
      endif





c HERE DEFINE THE FOUR MOMENTA - RIGHT NOW JUST ANY NUMBERS
      p_ext(0,1) = mdm(1)
      p_ext(1,1) = 0.d0
      p_ext(2,1) = 0.d0
      p_ext(3,1) = 0.d0
	 
      p_ext(0,2) = M_e
      p_ext(1,2) = 0.d0
      p_ext(2,2) = 0.d0
      p_ext(3,2) = 0.d0
		  
      p_ext(0,3) = mdm(1)
      p_ext(1,3) = 0.d0
      p_ext(2,3) = 0.d0
      p_ext(3,3) = 0.d0

      p_ext(0,4) = M_e
      p_ext(1,4) = 0.d0
      p_ext(2,4) = 0.d0
      p_ext(3,4) = 0.d0
      	

      Minterf_e = smatrix_dd_tot(p_ext, 1,1,jtot) - max(0d0, smatrix_dd(p_ext,1,1,j)) 
     &           - smatrix_dd_eff(p_ext,1,1,jeff)
      Minterf_e = 0.5d0*Minterf_e
      Minterf_e = Minterf_e / smatrix_dd_eff(p_ext,1,1,jeff)

cc    Add computation of Minterf_p (positron) or implement the same thing for
cc    dm and antidm scattering, in the end is the same thing

      write(*,*) '-----------------------------------------------------------------'
      write(*,*) 'process     :    ', DD_PROCESS_NAMES(j) , 'j =   ' , j
      write(*,*) 'EFF process :    ', DD_EFF_PROCESS_NAMES(jeff) , 'jeff =' , jeff
      write(*,*) 'TOT process :    ', DD_TOT_PROCESS_NAMES(jtot) , 'jtot =' , jtot
      write(*,*) 'full + eff  : ', smatrix_dd_tot(p_ext,1,1,jtot)
      write(*,*) 'full        : ', smatrix_dd(p_ext,1,1,j)
      write(*,*) 'eff         : ', smatrix_dd_eff(p_ext,1,1,jeff)
      write(*,*) 'Minterf_e   : ', Minterf_e
      write(*,*) '-----------------------------------------------------------------'

cc    Find the indices for the dm e- scatt and anti-dm e- scatt in order
cc    to compute alpha_even and alpha_odd

      if (DD_PROCESS_NAMES(j).eq.dm_e_process_name) then
      j_dm_e = j
      endif

      if (DD_PROCESS_NAMES(j).eq.antidm_e_process_name) then
            j_antidm_e = j
      endif

      do k=1, 2*dd_num_processes
            if (SI.eq.1) then

                  process_name = 'EFT_SI_' // dm_e_process_name
                  if (DD_EFF_PROCESS_NAMES(k).eq.process_name) then
                        j_eff_dm_e = k
                  endif

                  process_name = 'EFT_SI_' // antidm_e_process_name
                  if (DD_EFF_PROCESS_NAMES(k).eq.process_name) then
                        j_eff_antidm_e = k
                  endif

                  process_name = 'TOT_SI_' // dm_e_process_name
                  if (DD_TOT_PROCESS_NAMES(k).eq.process_name) then
                        j_tot_dm_e = k
                  endif

                  process_name = 'TOT_SI_' // antidm_e_process_name
                  if (DD_TOT_PROCESS_NAMES(k).eq.process_name) then
                        j_tot_antidm_e = k
                  endif

            else
                  process_name = 'EFT_SD_' // dm_e_process_name
                  if (DD_EFF_PROCESS_NAMES(k).eq.process_name) then
                        j_eff_dm_e = k
                  endif

                  process_name = 'EFT_SD_' // antidm_e_process_name
                  if (DD_EFF_PROCESS_NAMES(k).eq.process_name) then
                        j_eff_antidm_e = k
                  endif

                  process_name = 'TOT_SD_' // dm_e_process_name
                  if (DD_TOT_PROCESS_NAMES(k).eq.process_name) then
                        j_tot_dm_e = k
                  endif

                  process_name = 'TOT_SD_' // antidm_e_process_name
                  if (DD_TOT_PROCESS_NAMES(k).eq.process_name) then
                        j_tot_antidm_e = k
                  endif
            endif
      enddo


c compute ff_p or ff_n
c      if (even .eq. 1) then 
c        ! symmetric so do not have to distinguish quark from anti-quark
c        Minterf=0.5d0*Minterf_q
c        write(*,*) "Even contribution : ", Minterf
c      else
c        ! antisymmetric so we have to distinguish quark from anti-quark
c         if (i.ge.0) then
c            Minterf=0.5d0*Minterf_q
c         else
c            Minterf=-0.5d0*Minterf_q
c         endif
c	write(*,*) "Odd contribution : ", Minterf
c      endif
cc      write(*,*) 'M_interf: ', Minterf
c      NucleonFF=NucleonFF+(fNt(ABS(i))*Minterf)
c      write(*,*) 'fN: ', NucleonFF

      enddo

      write (*,*) "j_dm_e =" , j_dm_e
      write (*,*) "j_eff_dm_e =" , j_eff_dm_e
      write (*,*) "j_tot_dm_e =" , j_tot_dm_e
      write (*,*) "j_antidm_e =" , j_antidm_e
      write (*,*) "j_eff_antidm_e =" , j_eff_antidm_e
      write (*,*) "j_tot_antidm_e =" , j_tot_antidm_e

      Minterf_dm = smatrix_dd_tot(p_ext, 1,1,j_tot_dm_e) - max(0d0, smatrix_dd(p_ext,1,1,j_dm_e)) 
     &           - smatrix_dd_eff(p_ext,1,1,j_eff_dm_e)
      Minterf_dm = 0.5d0*Minterf_dm
      Minterf_dm = Minterf_dm / smatrix_dd_eff(p_ext,1,1,j_eff_dm_e)

      Minterf_antidm = smatrix_dd_tot(p_ext, 1,1,j_tot_antidm_e) - max(0d0, smatrix_dd(p_ext,1,1,j_antidm_e)) 
     &           - smatrix_dd_eff(p_ext,1,1,j_eff_antidm_e)
      Minterf_antidm = 0.5d0*Minterf_antidm
      Minterf_antidm = Minterf_antidm / smatrix_dd_eff(p_ext,1,1,j_eff_antidm_e)

      if (SI.eq.1) then
            alpha_SI_even = 0.5d0 * (Minterf_dm + Minterf_antidm)
            alpha_SI_odd  = 0.5d0 * (Minterf_dm - Minterf_antidm)
            c_1 = 1d0 / (4*mdm(1)*M_e) * (alpha_SI_odd + alpha_SI_even)
      else
            alpha_SD_even = 0.5d0 * (Minterf_dm + Minterf_antidm)
            alpha_SD_odd  = 0.5d0 * (Minterf_dm - Minterf_antidm)
            c_4 = 1d0 / (32*mdm(1)*M_e) * (alpha_SD_odd - 2*alpha_SD_even)
      endif

c     If I want to mantain the structure of this subroutine, in the end I have to output the c_1 and
c     c_4 coeff, depending on SI=1 or SI=0.
c     Then outside the routine we compute the dm response function: R_1DM = c_1**2 + (3d0/16d0)*c_4**2 

      write(*,*) 'alpha_SI_even  : ', alpha_SI_even
      write(*,*) 'alpha_SI_odd   : ', alpha_SI_odd
      write(*,*) 'alpha_SD_even  : ', alpha_SD_even
      write(*,*) 'alpha_SD_odd   : ', alpha_SD_odd
      write(*,*) 'c_1            : ', c_1
      write(*,*) 'c_4            : ', c_4
      
      end function
      
