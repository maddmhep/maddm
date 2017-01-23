c the subroutine provides a sigmav at average halo velocity
       function sigmav_ID(channel)

            implicit none
            include 'maddm.inc'

            double precision p1, p2, sigma_ID, vdm_natural
            integer i, j,channel

            vdm_natural = vmp/299792.d0 !vmp is the halo velocity in km/s
c    HERE AT SOME POINT CHECK THAT i =1, n and j=1,n are not double counting 1 <-> 2 etc.
            sigma_ID=0.d0
            do i=1, ndmparticles
                    do j =1, ndmparticles
c                        write(*,*) i, j, channel
                        call getfsmasses_ann(pmass,i,j,channel)
c                        write(*,*) pmass(1)
                        p1 = pmass(1)*vdm_natural
                        p2 = pmass(2)*vdm_natural
                        sigma_ID = sigma_ID + cross_check_process(i, j, channel, 1, p1,p2, 1)
                    enddo
            enddo
            sigmav_ID = sigma_ID*vdm_natural
       return
       end
