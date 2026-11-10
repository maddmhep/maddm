c     Inclusive Galactic annihilation rates averaged over the shared SHM.
      subroutine write_maddm_shm_rates(want_cont,want_line)
      implicit none
      include 'maddm.inc'
      logical want_cont, want_line, has_rule
      integer shm_unit, shm_io, shm_n, shm_k, shm_j
      double precision shm_v(128), shm_w(128), shm_p, shm_sum
      character(len=*), parameter :: shm_path = 'Cards/mw_shm_quadrature.dat'

      if (.not.(want_cont.or.want_line)) return
      inquire(file=shm_path,exist=has_rule)
      if (.not.has_rule) return
      open(newunit=shm_unit,file=shm_path,status='old',iostat=shm_io)
      if (shm_io.ne.0) stop 1
      read(shm_unit,*,iostat=shm_io) shm_n
      if (shm_io.ne.0.or.shm_n.lt.8.or.shm_n.gt.128) stop 1
      do shm_j=1,shm_n
         read(shm_unit,*,iostat=shm_io) shm_v(shm_j),shm_w(shm_j)
         if (shm_io.ne.0) stop 1
      enddo
      close(shm_unit)

      do shm_k=1,ANN_NUM_PROCESSES
         if (.not.want_cont.and..not.IS_LINE_PROCESS(shm_k)) cycle
         shm_sum = 0.d0
         do shm_j=1,shm_n
c           MadDM uses the nonrelativistic pair speed as the Moller flux
c           speed (2 times each incoming speed in the centre-of-mass frame).
            shm_p = mdm(1)*shm_v(shm_j)/sqrt(4.d0-shm_v(shm_j)**2)
            shm_sum = shm_sum + shm_w(shm_j)*shm_v(shm_j)*
     &         cross_check_process(1,1,shm_k,1,shm_p,shm_p,1)
         enddo
         write(33,fmt='(A,A,A,ES24.16,A)') 'mw_shm_sv:',
     &      trim(PROCESS_NAMES(shm_k)), ' ', shm_sum*pbtocm3, ' cm^3/s'
      enddo
      end
