c-------------------------------------------------------------------------c
      module class_atomic_response_nl
c-------------------------------------------------------------------------c
c
c This module contains the name of the shell and the corresponding
c response function, that can be read using the approprieate procedure.
c
c-------------------------------------------------------------------------c
            implicit none 

            type, public :: atomic_response_nl
                  integer shell_num, shell_tot_num
                  character(2) shell_name
                  real(kind=10) atomic_response_matrix(100,100)
                  contains
                  procedure :: read_atomic_response
            end type atomic_response_nl

      contains

            subroutine read_atomic_response(this,shell_filename,atom_resp_path)

                  Use, intrinsic :: iso_fortran_env, Only : iostat_end

                  character(20) shell_filename
                  character(100) atom_resp_path
                  integer j,k,error
                  class(atomic_response_nl), intent(inout)::this

                  open(2, FILE = trim(atom_resp_path)//trim(shell_filename), action="read")
c                            write(*,*) 'Opened ', trim(atom_resp_path) // trim(shell_filename)
                  this%%shell_name = shell_filename(4:5)
c                  Write(*, *) 'Printing shell : ', this(i)%%shell_name
                  Do j=1,100
                        Read(2, *, iostat = error) (this%%atomic_response_matrix(j,k), k=1,100)
                        Select Case(error)
                        Case(0)
c                              do k=1,100
c                                          Write(*, *) '(', j , ',' , k , ')' , this%%atomic_response_matrix(j,k)
c                              enddo
c                                    Write(*,*) ''
                        Case(iostat_end)
                              Exit
                        Case Default
                              Write(*, *) 'Error in reading file'
                              Stop
                        End Select
                  enddo
                  close(2)

            end subroutine read_atomic_response

      end module class_atomic_response_nl





c-----------------------------------------------------------------------------c
      subroutine get_shell_filenames(atom_resp_path,target,shell_filenames,num_shell)
c-----------------------------------------------------------------------------c
c
c Read all names of the files that contains the tabulated atomic response
c function of the chosen target. There is one file for each shell.
c
c-----------------------------------------------------------------------------c
            Use, intrinsic :: iso_fortran_env, Only : iostat_end

            character(100) atom_resp_path
            character(2) target
            character(20) filename, shell_filenames(50)
            integer error, num_shell

            num_shell = 0

            call system('ls ' // trim(atom_resp_path) // ' > ' // trim(atom_resp_path) // 'fileNames.txt')
            open(1, FILE = trim(atom_resp_path)//'fileNames.txt', action="read")
                  Do
                        Read(1, *, iostat = error) filename
                        Select Case(error)
                        Case(0)
                              if(filename(1:3).eq.(trim(target)//'_')) then
                                    num_shell = num_shell + 1
                                    shell_filenames(num_shell) = filename
                              endif
                        Case(iostat_end)
                              Exit
                        Case Default
                              Write(*, *) 'Error in reading file'
                              Stop
                        End Select
                  enddo
            close(1)

      end subroutine get_shell_filenames





c-------------------------------------------------------------------------c
      function sigma_nucleon(proton,SI)
c-------------------------------------------------------------------------c
c                                                                         c
c     Compute the squared ionization amplitude for the DM-electron        c
c     scattering, considering the atomic response function of the target  c
c     atom.                                                               c
c     (remember to rename the function as 'ioniz_amplitude')              c
c                                                                         c
c     Please cite: arXiv:1912.08204v2 [hep-ph] 19 Aug 2020                c
c                                                                         c
c-------------------------------------------------------------------------c

            use class_atomic_response_nl
            use, intrinsic :: iso_fortran_env, Only : iostat_end
            
            integer proton, SI, dm_spin
            double precision M_e
            double precision dm_res, dm_response
            real(kind=10) ioniz_amplitude(100,100)

            integer i, j, k, tot_num_shell
            character(20) shell_filenames(50)
            character(2) target
            character(100) maddm_path, atom_resp_path
            type(atomic_response_nl) atom_res_nl(50)


            include 'maddm.inc'
            include 'coupl.inc'
            include 'process_names.inc'
            include 'maddm_card.inc'
            include 'dm_info.inc'

c           Initialize the parameters
            %(electron_mass)s       ! electron mass M_e
            dm_spin = dof_dm(1)     ! spin of the dm particle
            %(maddm_path)s          ! MadDM path
            atom_resp_path = trim(maddm_path) // '/Atomic_responses/'
            do j=1,100
                  do k=1,100
                              ioniz_amplitude(j,k) = 0
                  enddo
            enddo


            if (proton.eq.1.and.SI.eq.1) then        ! to avoid that the function is called four times

c --------------------------------------------------------------------------------
c Read the atomic response functions
c --------------------------------------------------------------------------------
            target = 'Xe'

c           Read all the names of the files that contains the tabulated atomic response function
c           of the chosen target. There is one file for each shell.
            call get_shell_filenames(atom_resp_path,target,shell_filenames,tot_num_shell)

c           Read the atomic response function of each file.
            do i=1,(tot_num_shell)
                  atom_res_nl(i)%%shell_num = i
                  atom_res_nl(i)%%shell_tot_num = tot_num_shell
                  call atom_res_nl(i)%%read_atomic_response(shell_filenames(i),atom_resp_path)
            enddo

c --------------------------------------------------------------------------------
c Select scalar, fermion or vector DM particle and compute the dark matter
c response function
c --------------------------------------------------------------------------------
            select case(dm_spin)
                  case(1) ! will be implmented in future work
                  write (*,*) "scalar DM-e scattering still not implemented!"
                  case(2)
                  dm_res = dm_response(M_e)
                  case(3) ! will be implmented in future work
                  write (*,*) "vector DM-e scattering still not implemented!"
            end select

c --------------------------------------------------------------------------------
c Multiply the dark matter response function with the atomic one, and sum over
c all electronic shells
c --------------------------------------------------------------------------------
            do i=1,(tot_num_shell)
                  do j=1,100
                        do k=1,100
                              ioniz_amplitude(j,k) = ioniz_amplitude(j,k) + atom_res_nl(i)%%atomic_response_matrix(j,k)
                        enddo
                  enddo
            enddo

            do j=1,100
                  do k=1,100
                        ioniz_amplitude(j,k) = dm_res * ioniz_amplitude(j,k)
                  enddo
            enddo

c            write(*,*) 'DM response: ', dm_res
c            write(*,*) 'Ionization amplitude:'
c            do j=1,5
c                  do k=1,5
c                        write(*,*) ioniz_amplitude(j,k)
c                  enddo
c                  write(*,*) ''
c            enddo

            endif
            
      end function
            



            
c-------------------------------------------------------------------------c
      function dm_response(M_e)
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
            double precision M_e
            double precision Minterf_SI_e, Minterf_SI_p, Minterf_SD_e, Minterf_SD_p
            double precision alpha_SI_even, alpha_SI_odd, alpha_SD_even, alpha_SD_odd, c_1, c_4
            double precision dm_response
            character(2) siorsd
            integer len_process_name
            character(50) process_name , dm_e_process_name, dm_p_process_name

            include 'maddm.inc'
            include 'coupl.inc'
            include 'process_names.inc'
            include 'maddm_card.inc'

            dm_e_process_name = trim(DM_NAMES(1)) // 'em_' // trim(DM_NAMES(1)) // 'em'
            dm_p_process_name = trim(DM_NAMES(1)) // 'ep_' // trim(DM_NAMES(1)) // 'ep'
c            write (*,*) 'dm_e_process_name = ' , dm_e_process_name
c            write (*,*) 'dm_p_process_name = ' , dm_p_process_name
c            write(*,*) 'dm info : ' , dof_dm(1)
c            write(*,*) 'number of processes: ', dd_num_processes

c           Four momenta of the particles involved in the scattering
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
c Compute the non relativistic coefficients and use them to compute the
c dark matter response function.
c--------------------------------------------------------------------------------
            c_1 = 4*mdm(1)*M_e*(alpha_SI_odd + alpha_SI_even)
            c_4 = 16*mdm(1)*M_e*(2*alpha_SD_odd - alpha_SD_even)
            dm_response = (c_1)**2 + (3.d0/16)*(c_4)**2

            write(*,*) '-------------------------------------------'
c          write(*,*) 'Minterf_SI_e   : ',Minterf_SI_e
c          write(*,*) 'Minterf_SI_p   : ',Minterf_SI_p
c          write(*,*) 'Minterf_SD_e   : ',Minterf_SD_e
c           write(*,*) 'Minterf_SD_p   : ',Minterf_SD_p
            write(*,*) 'alpha_SI_even  : ', alpha_SI_even
            write(*,*) 'alpha_SI_odd   : ', alpha_SI_odd
            write(*,*) 'alpha_SD_even  : ', alpha_SD_even
            write(*,*) 'alpha_SD_odd   : ', alpha_SD_odd
c          write(*,*) 'c_1            : ', c_1
c          write(*,*) 'c_4            : ', c_4
c          write(*,*) 'dm_response    : ', dm_response
            write(*,*) '-------------------------------------------'
      
      end function
