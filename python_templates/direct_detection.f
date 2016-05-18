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

c initialize the quark masses array (u,d,s,c,b,t)
c--------------------------------------------------------------------------------

      %(quark_masses)s

      if (proton .eq. 1) then
         MN=0.9383d0
      else
         MN=0.9396d0
      endif
        
      mr=(MN*mdm(1))/(MN+mdm(1))
      
      spin = dof_dm(1) !get the spin of the dm particle
      
      if (SI.eq.1) then
      ff_N=0
      ff_N_odd=0
c--------------------------------------------------------------------------------
c              SPIN INDEPENDENT
c--------------------------------------------------------------------------------
      
c compute nucleon form factor for scalar and vector interactions (so SI interactions)
c--------------------------------------------------------------------------------
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
          even=0
          ff_N_odd=NucleonFF(M,proton,SI,even) !projection over ODD operator (vector interaction)
          ff_N=ff_N+(ff_N_odd)
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
      ff_N=0
      ff_N_odd=0
      select case(spin) !select fermion or vector DM particle
        case(2)
          even=1
          ff_N=NucleonFF(M,proton,SI,even) !projection over EVEN operator (axial-vector interaction)
          even=0
          ff_N_odd=NucleonFF(M,proton,SI,even) !projection over ODD operator (tensor interaction)
          ff_N=ff_N+ff_N_odd
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
      
      end function
      








      
      
      
      
      function NucleonFF(M,proton,SI,even)
      
      implicit none
      
      integer i, j ,proton, even, SI
      double precision p_ext(0:3,4)
      double precision fNt(6) 
      double precision M(6)
      double precision Minterf_q, Minterf_qb, Minterf
      complex*16 NucleonFF
      
      include 'maddm.inc'
      include 'coupl.inc'
      include 'process_names.inc'
      include 'maddm_card.inc'
      
      NucleonFF=0.0d0
      
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
      
c      write(*,*) 'dm info : ' , dof_dm(1)
      
c      write(*,*) 'number of processes: ', dd_num_processes
c--------------------------------------------------------------------------------
c loop over the number of processes (min=1 , max=6)
c--------------------------------------------------------------------------------
      do j=1, (dd_num_processes/2) !divide by 2 since dd_num_processes from 1 -> 12
      
c HERE BE CAREFUL. THERE ARE ONLY 2 MATRIX ELEMENTS FOR FULL LAGRANGIAN
c BUT 6 FOR EFF and EFF+FULL. THE INDEX i IS USED TO MATCH THE FULL LAGRANGIAN ME
C WITH THE EFF+FULL AND EFF. WE WILl FIGURE OUT A BETTER WAY TO DO THIS LATER

      i = dd_process_ids(j)
      
c      write(*,*) 'process id: ', i
c HERE DEFINE THE FOUR MOMENTA - RIGHT NOW JUST ANY NUMBERS
      p_ext(0,1) = mdm(1)
      p_ext(1,1) = 0.d0
      p_ext(2,1) = 0.d0
      p_ext(3,1) = 0.d0
	 
      p_ext(0,2) = M(i)
      p_ext(1,2) = 0.d0
      p_ext(2,2) = 0.d0
      p_ext(3,2) = 0.d0
		  
      p_ext(0,3) = mdm(1)
      p_ext(1,3) = 0.d0
      p_ext(2,3) = 0.d0
      p_ext(3,3) = 0.d0

      p_ext(0,4) = M(i)
      p_ext(1,4) = 0.d0
      p_ext(2,4) = 0.d0
      p_ext(3,4) = 0.d0
      	
c compute ff_p or ff_n
c      write(*,*) '--------------------------------'
c      write(*,*) 'full + eff: ',smatrix_dd_tot(p_ext,1,1,i)
c      write(*,*) 'full      : ',smatrix_dd(p_ext,1,1,j)
c      write(*,*) 'eff       : ',smatrix_dd_eff(p_ext,1,1,i)
c      write(*,*) '--------------------------------'
      
      if (SI .eq. 1) then
c compute \lambda_e + \lambda_o      
        Minterf_q = smatrix_dd_tot(p_ext, 1,1,i) - max(0d0, smatrix_dd(p_ext,1,1,j)) - smatrix_dd_eff(p_ext,1,1,i)
        Minterf_q = 0.5d0*Minterf_q
        Minterf_q = Minterf_q / smatrix_dd_eff(p_ext,1,1,i)
c compute \lambda_e - \lambda_o   ! i+6 and j+6 to get anti-quark  
        Minterf_qb=smatrix_dd_tot(p_ext,1,1,i+6)-max(0d0,smatrix_dd(p_ext,1,1,j+(dd_num_processes/2)))
     $    -smatrix_dd_eff(p_ext,1,1,i+6)
        Minterf_qb = 0.5d0*Minterf_qb
        Minterf_qb = Minterf_qb / smatrix_dd_eff(p_ext,1,1,i+6)
c solve the system to get \lambda_e or \lambda_o  
      else
c compute \lambda_e + \lambda_o      
        Minterf_q = smatrix_dd_tot(p_ext, 1,1,i+12) - max(0d0, smatrix_dd(p_ext,1,1,j)) - smatrix_dd_eff(p_ext,1,1,i+12)
        Minterf_q = 0.5d0*Minterf_q
        Minterf_q = Minterf_q / smatrix_dd_eff(p_ext,1,1,i+12)
c compute \lambda_e - \lambda_o   ! i+18 to get anti-quark  
        Minterf_qb=smatrix_dd_tot(p_ext,1,1,i+18)-max(0d0,smatrix_dd(p_ext,1,1,j+(dd_num_processes/2)))
     $    -smatrix_dd_eff(p_ext,1,1,i+18)
        Minterf_qb = 0.5d0*Minterf_qb
        Minterf_qb = Minterf_qb / smatrix_dd_eff(p_ext,1,1,i+18)
c solve the system to get \lambda_e or \lambda_o  
      endif      
      if (even .eq. 1) then 
        Minterf=0.5d0*(Minterf_q+Minterf_qb)
c        write(*,*) "Even contribution : ", Minterf
      else
        Minterf=0.5d0*(Minterf_q-Minterf_qb)
c	write(*,*) "Odd contribution : ", Minterf
      endif
c      write(*,*) 'M_interf: ', Minterf
      NucleonFF=NucleonFF+(fNt(dd_process_ids(j))*Minterf)
c      write(*,*) 'fN: ', NucleonFF
      enddo
      
      end function
      
