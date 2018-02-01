c-------------------------------------------------
c      computation of sigmav at present time
c-------------------------------------------------

      function taacs_ID(channel, vave)
      implicit none
      include 'maddm.inc'
      double precision vave
      double precision channel

      taacs_ID = cross_check_process(1, 1, channel, 1,mdm(1)*vave,mdm(1)*vave , 1)
c      print*, taacs_ID

      return
      end function taacs_ID


c$$$
c-------------------------------------------------
c$$$c the function provides a sigma*v at RELATIVE velocity v (not thermally averaged!!!)
c$$$c-------------------------------------------------
c$$$       function sigmav_ID(channel, v)
c$$$
c$$$            implicit none
c$$$            include 'maddm.inc'
c$$$
c$$$            double precision e1, e2, sigma_ID, v
c$$$            integer i, j,channel
c$$$
c$$$            if (v.le.0.d0) then
c$$$                sigmav_ID = 0.d0
c$$$                return
c$$$            endif
c$$$c            vdm_natural = vmp/299792.d0 !vmp is the halo velocity in km/s
c$$$c    HERE AT SOME POINT CHECK THAT i =1, n and j=1,n are not double counting 1 <-> 2 etc.
c$$$            sigma_ID=0.d0
c$$$            do i=1, ndmparticles
c$$$                    do j =1, ndmparticles
c$$$c                        write(*,*) i, j, channel
c$$$                        call getfsmasses_ann(pmass,i,j,channel)
c$$$c                       write(*,*) pmass(1)
c$$$                        e1 = dsqrt(pmass(1)*pmass(1)/(1.d0-v*v/4.d0)) ! divide by 4 for vcm
c$$$                        e2 = e1
c$$$                        sigma_ID = sigma_ID + cross_check_process(i, j, channel, 1, e1,e2, 0)
c$$$                    enddo
c$$$            enddo
c$$$            sigmav_ID = sigma_ID*v
c$$$       return
c$$$       end
c$$$
c$$$c-----------------------------------------------------------------------
c$$$       subroutine set_up_grid_ID(vave)
c$$$
c$$$            implicit none
c$$$            include 'maddm.inc'
c$$$            double precision width, additional_pt, vave
c$$$            integer channel, kk, jj, grid_pos
c$$$
c$$$c take as the initial grid the one from relic density
c$$$c which has all the resonance locations already in
c$$$            grid_ID = grid
c$$$
c$$$            width = vave ! width of the maxwellian = most probable vel.
c$$$
c$$$c navigate the grid position to the location of v =1.0
c$$$            grid_pos = 1
c$$$            do while (grid_ID(grid_pos).le.1.d0)
c$$$                grid_pos = grid_pos + 1
c$$$            enddo
c$$$
c$$$            grid_ID(grid_pos)=width
c$$$
c$$$            do jj=1, nres_points
c$$$                    additional_pt = width +5.d0*width/nres_points*jj
c$$$                    if (additional_pt.lt.1.d0) then
c$$$                        grid_ID(grid_pos) = additional_pt
c$$$                        grid_pos=grid_pos+1
c$$$                        if (grid_npts.lt.grid_pos.or.grid_pos.lt.1) then
c$$$                            write(*,*) 'Error 1 (ID grid): grid array not large enough!'
c$$$                            call exit(1)
c$$$                        endif
c$$$                    endif
c$$$
c$$$                    additional_pt = width - 5.d0*width/nres_points*jj
c$$$                    if (additional_pt.gt.0.d0) then
c$$$                        grid_ID(grid_pos) = additional_pt
c$$$                        grid_pos=grid_pos+1
c$$$                        if (grid_npts.lt.grid_pos.or.grid_pos.lt.1) then
c$$$                            write(*,*) 'Error 2 (ID grid): grid array not large enough!'
c$$$                            call exit(1)
c$$$                        endif
c$$$                    endif
c$$$               enddo
c$$$
c$$$
c$$$            call Duplicates(grid_ID, grid_npts)
c$$$            call  Bubble_Sort(grid_ID, grid_npts)
c$$$
c$$$       end
c$$$
c$$$c---------------------------------------------------------------------
c$$$c the function provides a thermally averaged cross section at vave
c$$$c---------------------------------------------------------------------
c$$$       function taacs_ID(channel, vave)
c$$$
c$$$            implicit none
c$$$            include 'maddm.inc'
c$$$            external integrand
c$$$            double precision vave
c$$$            integer channel
c$$$
c$$$c           Here set the global variables used by integrand bexfore you integrate it
c$$$            vave_glob = vave
c$$$            channel_glob=channel
c$$$
c$$$            taacs_ID = simpson(integrand, 0.d0, min(20.d0*vave,1.d0), grid_ID, grid_npts)
c$$$
c$$$            return
c$$$       end function taacs_ID
c$$$
c$$$c---------------------------------------------------------------------
c$$$c   Non relativistic approximation to the velocity dist (relative velocity).
c$$$c-------------------------------------------------------------------
c$$$       function dm_vel_dist_ID(v, vave)
c$$$
c$$$            implicit none
c$$$            include 'maddm.inc'
c$$$            double precision v, vave
c$$$
c$$$            dm_vel_dist_ID = dsqrt(2.d0/pi)*(1.d0/(vave*vave))**(1.5)*v*v*dexp(-v*v/(2.d0*vave*vave))
c$$$            return
c$$$       end function dm_vel_dist_ID
c$$$
c$$$
c$$$c---------------------------------------------------------------------
c$$$c  integrand needs to be a function only of v, so vave_glob and channel_glob
c$$$c   will be the global variables which it uses to know which channel and which
c$$$c   average velocity to use
c$$$c----------------------------------------------------------------------
c$$$       function integrand(v)
c$$$            implicit none
c$$$            include 'maddm.inc'
c$$$            double precision v
c$$$
c$$$            integrand = sigmav_ID(channel_glob, v) *dm_vel_dist_ID(v, vave_glob) ! v is vrelative
c$$$c            print*,sigmav_ID(channel_glob, v), dm_vel_dist_ID(v, vave_glob), v, channel_glob, vave_glob
c$$$            if (integrand.lt.1d-50) then
c$$$                integrand=0.d0
c$$$            endif
c$$$            return
c$$$       end function integrand
