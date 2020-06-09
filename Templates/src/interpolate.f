c-------------------------------------------------------------------------c
      function get_Wij_value(beta,i,j,k,group)
c-------------------------------------------------------------------------c
c                                                                         c
c  This function finds the linear interpolation of the Wij values         c
c  calculated for any group of diagrams given an input velocity beta and  c
c  the desired group of diagrams.                                         c
c                                                                         c
c  group = 1 for annihilation diagrams.                                   c
c  group = 2 for DM -> DM diagrams.                                       c
c  group = 3 for DM/SM scattering diagrams.                               c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      include 'maddm.inc'

c input parameters
      double precision beta
      integer i, j, k, group

c parameters used to find the index
      integer index, lowindex, highindex, mid

      get_Wij_value = 0d0
c If the xvalue is outside the given array then set the index to the boundary
      if (beta .lt. betas(1)) then
        index = 1
      else if (beta .gt. betas(nWij-1)) then
        index = nWij-1
      endif

c Binary search
      lowindex = 1
      highindex = nWij

      do while ((highindex-lowindex).gt.1)
        mid = (highindex+lowindex)/2
        if (betas(mid) .le. beta) lowindex = mid
        if (betas(mid) .gt. beta) highindex = mid
      enddo

      index = lowindex

c interpolate the array with a line
      if (group.eq.1) then
c         write(*,*) 'bypass interpolation'
c        get_Wij_value = get_Wij_ann_nogrid(beta)
        get_Wij_value = (Wij_ann(i,j,index+1)-Wij_ann(i,j,index))/(betas(index+1)-betas(index))
     .        *(beta-betas(index))+Wij_ann(i,j,index)
      else if (group.eq.2) then
        get_Wij_value = (Wij_dm2dm(i,j,k,index+1)-Wij_dm2dm(i,j,k,index))/(betas(index+1)-betas(index))
     .        *(beta-betas(index))+Wij_dm2dm(i,j,k,index)
      else if (group.eq.3) then
        get_Wij_value = (Wij_scattering(i,j,k,index+1)-Wij_scattering(i,j,k,index))/(betas(index+1)-betas(index))
     .        *(beta-betas(index))+Wij_scattering(i,j,k,index)
      endif
    
      return
      end



c-------------------------------------------------------------------------c
      function get_taacs_value(x,i,j,k,group)
c-------------------------------------------------------------------------c
c                                                                         c
c  Given an input value of x and the group of processes this function     c
c  returns the interpolated value of the thermally averaged annihilation  c
c  cross section for that particular group of processes.                  c
c                                                                         c
c  group = 1 for annihilation processes.                                  c
c  group = 2 for DM -> DM processes (these require initial state ids      c
c  i and j as well as process k).                                         c
c  group = 3 for DM/SM scattering processes.                              c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      include 'maddm.inc'

c input parameters
      double precision x
      integer i, j, k, group

c parameters used to find the index
      integer x_index, lowindex, highindex, mid
      get_taacs_value = 0d0
c If the xvalue is outside the given array then set the index to the boundary
      if (x .lt. x_taacs(1)) then
        x_index = 1
      else if (x .gt. x_taacs(99)) then
        x_index = 99
      else

c Binary search
        lowindex = 1
        highindex = 100

        do while ((highindex-lowindex).gt.1)
          mid = (highindex+lowindex)/2
          if (x_taacs(mid) .le. x) lowindex = mid
          if (x_taacs(mid) .gt. x) highindex = mid
        enddo

        x_index = lowindex
      endif

c interpolate the array with a line
      if (group.eq.1) then
        get_taacs_value = (taacs_ann(i,j,x_index+1)-taacs_ann(i,j,x_index))
     .        /(x_taacs(x_index+1)-x_taacs(x_index))*(x-x_taacs(x_index))+taacs_ann(i,j,x_index)
      else if (group.eq.2) then
        get_taacs_value = (taacs_dm2dm(i,j,k,x_index+1)-taacs_dm2dm(i,j,k,x_index))
     .        /(x_taacs(x_index+1)-x_taacs(x_index))*(x-x_taacs(x_index))+taacs_dm2dm(i,j,k,x_index)
      else if (group.eq.3) then
        get_taacs_value = (taacs_scattering(i,j,k,x_index+1)-taacs_scattering(i,j,k,x_index))
     .        /(x_taacs(x_index+1)-x_taacs(x_index))*(x-x_taacs(x_index))+taacs_scattering(i,j,k,x_index)
      endif

      return
      end



c-------------------------------------------------------------------------c
      function getindex(xvalue,xarray,n)
c-------------------------------------------------------------------------c
c                                                                         c
c  This function returns the index of the array, i, such that xarray(i)   c
c  < xvalue < xarray(i+1).  If xvalue is outside of the array then        c
c  getindex returns the boundary.                                         c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      include 'maddm.inc'

c input parameters
      integer n
      double precision xvalue, xarray(n)

c parameters used in this subroutine only      
      integer lowindex, highindex, mid
      double precision low, high
      
c If the xvalue is outside the given array then set the index to the boundary
      if (xvalue .lt. xarray(1)) then
        getindex = 1
        return
      else if (xvalue .gt. xarray(n)) then
        getindex = n-1
        return
      endif
      
c Binary search
      lowindex = 1
      highindex = n

      do while ((highindex-lowindex).gt.1)
        mid = (highindex+lowindex)/2
        if (xarray(mid).le.xvalue) lowindex = mid
        if (xarray(mid).gt.xvalue) highindex = mid
      enddo

c Now we have our index
      getindex = lowindex

      return
      end
    
    

c-------------------------------------------------------------------------c
      function lineint(xvalue,xarray,yarray,n)
c-------------------------------------------------------------------------c
c                                                                         c
c  This function finds the linear interpolation yvalue given an input     c
c  value xvalue and arrays xarray and yarray (both arrays are length n).  c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      include 'maddm.inc'
      
c input parameters
      integer n
      double precision xvalue, xarray(n), yarray(n)

c paramters used in this subroutine only
      integer index
      
c find the index of the array we need    
      index = getindex(xvalue, xarray, n)

c interpolate the array with a line      
      lineint = (yarray(index+1)-yarray(index))/(xarray(index+1)-xarray(index))
     .  *(xvalue-xarray(index))+yarray(index)

      return
      end
      
      
      
c-------------------------------------------------------------------------c
      subroutine readinfile(filename,n,xarray,yarray)
c-------------------------------------------------------------------------c
c                                                                         c
c  This subroutine is designed to make it very easy to read in a file     c
c  consisting of two columns of double precision data.  In conjunction    c
c  with getvalue we can interpolate these curves at any point.            c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none

      integer n,i
      double precision xarray(n),yarray(n) 
      character*30 filename

      open(unit=10,file=filename,status='old')
      
      do i=1,n
       read(10,*) xarray(i),yarray(i)
      enddo

      close(10)

      return
      end



c-------------------------------------------------------------------------c
c bessel functions (taken from numerical recipes)
c-------------------------------------------------------------------------c
      function bessk(n,x)
c-------------------------------------------------------------------------c
c                                                                         c
c  Returns the modified Bessel function K_n(x) for positive x and n >= 2  c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      
      integer n
      double precision bessk, x

c uses bessk0, bessk1
      integer i
      double precision bk, bkm, bkp, tox, bessk0, bessk1

c      if (n.lt.2) pause 'bad argument n in bessk'
      tox = 2.d0/x
      bkm = bessk0(x)
      bk = bessk1(x)
      do i=1,n-1
        bkp = bkm + i*tox*bk
        bkm = bk
        bk = bkp
      enddo
      bessk = bk

      return
      end



c-------------------------------------------------------------------------c
      function bessk1(x)
c-------------------------------------------------------------------------c
c                                                                         c
c  Returns the modified Bessel function K_{1}(x) for positive real x.     c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      
      double precision bessk1, x

c uses bessi1
      double precision bessi1
      double precision p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,y
      save p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
      data p1,p2,p3,p4,p5,p6,p7/1.0d0,0.15443144d0,-0.67278579d0,
     . -0.18156897d0,-0.1919402d-1,-0.110404d-2,-0.4686d-4/
      data q1,q2,q3,q4,q5,q6,q7/1.25331414d0,0.23498619d0,-0.3655620d-1,
     . 0.1504268d-1,-0.780353d-2,0.325614d-2,-0.68245d-3/

      if (x.le.2.d0) then
       y=x*x/4.d0
       bessk1=(dlog(x/2.d0)*bessi1(x))+(1.d0/x)*(p1+y*(p2+y*(p3+y*(p4+y*
     .  (p5+y*(p6+y*p7))))))
      else
       y=2.d0/x
       bessk1=(dexp(-x)/dsqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*
     .  q7))))))
      endif

      return
      end



c-------------------------------------------------------------------------c
      function bessk0(x)
c-------------------------------------------------------------------------c
c                                                                         c
c  Returns the modified Bessel function K_{0}(x) for positive real x.     c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none

      double precision bessk0,x
      
c uses bessi0
      double precision bessi0
      double precision p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,y
      save p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
      data p1,p2,p3,p4,p5,p6,p7/-0.57721566d0,0.42278420d0,0.23069756d0,
     . 0.3488590d-1,0.262698d-2,0.10750d-3,0.74d-5/
      data q1,q2,q3,q4,q5,q6,q7/1.25331414d0,-0.7832358d-1,0.2189568d-1,
     . -0.1062446d-1,0.587872d-2,-0.251540d-2,0.53208d-3/

      if (x.le.2.d0) then
        y=x*x/4.d0
        bessk0=(-dlog(x/2.d0)*bessi0(x))+(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*
     . (p6+y*p7))))))
      else
        y=(2.d0/x)
        bessk0=(dexp(-x)/dsqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*
     . q7))))))
      endif

      return
      end



c-------------------------------------------------------------------------c
      function bessi1(x)
c-------------------------------------------------------------------------c
c                                                                         c
c  Retuns the modified Bessel function I_{1}(x) for real x.               c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      
      double precision bessi1, x
      double precision ax
      double precision p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
      save p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
      data p1,p2,p3,p4,p5,p6,p7/0.5d0,0.87890594d0,0.51498869d0,
     *0.15084934d0,0.2658733d-1,0.301532d-2,0.32411d-3/
      data q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,-0.3988024d-1,
     *-0.362018d-2,0.163801d-2,-0.1031555d-1,0.2282967d-1,-0.2895312d-1,
     *0.1787654d-1,-0.420059d-2/

      if (dabs(x).lt.3.75d0) then
        y=(x/3.75d0)**2
        bessi1=x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
      else
        ax=dabs(x)
        y=3.75d0/ax
        bessi1=(dexp(ax)/dsqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*
     *(q7+y*(q8+y*q9))))))))
        if (x.lt.0.d0) bessi1=-bessi1
      endif

      return
      end



c-------------------------------------------------------------------------c
      function bessi0(x)
c-------------------------------------------------------------------------c
c                                                                         c
c  Returns the modified Bessel function I_{0}(x) for real x.              c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      
      double precision bessi0, x
      double precision ax
      double precision p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
      save p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
      data p1,p2,p3,p4,p5,p6,p7/1.0d0,3.5156229d0,3.0899424d0,
     . 1.2067492d0,0.2659732d0,0.360768d-1,0.45813d-2/
      data q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,0.1328592d-1,
     . 0.225319d-2,-0.157565d-2,0.916281d-2,-0.2057706d-1,0.2635537d-1,
     . -0.1647633d-1,0.392377d-2/

      if (dabs(x).lt.3.75d0) then
        y=(x/3.75d0)**2
        bessi0=p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
      else
        ax=dabs(x)
        y=3.75d0/ax
        bessi0=(dexp(ax)/dsqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*
     . (q7+y*(q8+y*q9))))))))
      endif

      return
      end
