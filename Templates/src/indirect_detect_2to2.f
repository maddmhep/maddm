c-------------------------------------------------
c the functioneprovides a sigma*v at velocity v (not thermally averaged!!!)
c-------------------------------------------------
       function sigmav_ID(channel, v)

            implicit none
            include 'maddm.inc'

            double precision e1, e2, sigma_ID, v
            integer i, j,channel

            if (v.le.0.d0) then
                sigmav_ID = 0.d0
                return
            endif
c            vdm_natural = vmp/299792.d0 !vmp is the halo velocity in km/s
c    HERE AT SOME POINT CHECK THAT i =1, n and j=1,n are not double counting 1 <-> 2 etc.
            sigma_ID=0.d0
            do i=1, ndmparticles
                    do j =1, ndmparticles
c                        write(*,*) i, j, channel
                        call getfsmasses_ann(pmass,i,j,channel)
c                       write(*,*) pmass(1)
                        e1 = dsqrt(pmass(1)*pmass(1)/(1.d0-v*v))
                        e2 = e1
                        sigma_ID = sigma_ID + cross_check_process(i, j, channel, 1, e1,e2, 0)
                    enddo
            enddo
            sigmav_ID = sigma_ID*v
       return
       end

c-----------------------------------------------------------------------
       subroutine set_up_grid_ID(vave)

            implicit none
            include 'maddm.inc'
            double precision width, additional_pt, vave
            integer channel, kk, jj, grid_pos

            width = vave ! width of the maxwellian = most probable vel.
            grid_pos = grid_npts_ID
            grid_ID(grid_pos)=width
            grid_pos=grid_pos-1
            do jj=1, nres_points
                    additional_pt = width +5.d0*width/nres_points*jj
                    if (additional_pt.lt.1.d0) then
                        grid_ID(grid_pos) = additional_pt
                        grid_pos=grid_pos-1
                        if (grid_npts_ID.lt.grid_pos.or.grid_pos.lt.1) then
                            write(*,*) 'Error: grid array not large enough!'
                            call exit(1)
                        endif
                    endif

                    additional_pt = width - 5.d0*width/nres_points*jj
                    if (additional_pt.gt.0.d0) then
                        grid_ID(grid_pos) = additional_pt
                        grid_pos=grid_pos-1
                        if (grid_npts_ID.lt.grid_pos.or.grid_pos.lt.1) then
                            write(*,*) 'Error: grid array not large enough!'
                            call exit(1)
                        endif
                    endif
               enddo


            call Duplicates(grid_ID, grid_npts_ID)
            call  Bubble_Sort(grid_ID, grid_npts_ID)

       end

c---------------------------------------------------------------------
c the function provides a thermally averaged cross section at vave
c---------------------------------------------------------------------
       function taacs_ID(channel, vave)

            implicit none
            include 'maddm.inc'
            external integrand
            double precision vave
            integer channel

c           Here set the global variables used by integrand bexfore you integrate it
            vave_glob = vave
            channel_glob=channel
c------------------------------------>>>>>>>>
            taacs_ID = simpson(integrand, 0.d0, min(10.d0*vave_glob,1.d0), grid_ID, grid_npts_ID)

            return
       end function taacs_ID

c---------------------------------------------------------------------
c   Non relativistic approximation to the velocity dist.
c-------------------------------------------------------------------
       function dm_vel_dist_ID(v, vave)

            implicit none
            include 'maddm.inc'
            double precision v, vave

            dm_vel_dist_ID = dsqrt(2.d0/pi)*(1.d0/(vave*vave))**(1.5)*v*v*dexp(-v*v/(2.d0*vave*vave))
            return
       end function dm_vel_dist_ID


c---------------------------------------------------------------------
c  integrand needs to be a function only of v, so vave_glob and channel_glob
c   will be the global variables which it uses to know which channel and which
c   average velocity to use
c----------------------------------------------------------------------
       function integrand(v)
            implicit none
            include 'maddm.inc'
            double precision v

            integrand = sigmav_ID(channel_glob, v/2.d0) *dm_vel_dist_ID(v, vave_glob)
            if (integrand.lt.1d-50) integrand=0.d0
            return
       end function integrand