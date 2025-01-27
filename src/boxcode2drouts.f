c
c     common boxcode tasks
c

      subroutine bc2d_time(t)
      implicit none
      real *8 :: t
c$    real *8 :: omp_get_wtime       

      call cpu_time(t)
c$    t = omp_get_wtime()

      return
      end
      
      subroutine bc2d_dformcoefs(nboxes,ichildbox,npts,npoly,pt2poly,
     1     fval,fcoef)
      implicit real *8 (a-h,o-z)
      integer :: ichildbox(4,*)
      real *8 :: fval(npts,*),fcoef(npoly,*),pt2poly(npoly,npts)

      do i = 1,nboxes
         if (ichildbox(1,i) .le. 0) then
            call bc2d_dmatvec(npoly,npts,pt2poly,fval(1,i),fcoef(1,i))
         endif
      enddo

      return
      end

c
c     low-level matrix/vector operations wrappers
c      
      
      subroutine bc2d_dscalzmatdvec(m,n,scal,a,x,y)
      implicit none
      real *8 scal, x(n)
      complex *16 a(m,n), y(m)
      integer m,n,i

      do i = 1,m
         y(i)=0
      enddo

      call bc2d_dscalzmatdvec_add(m,n,scal,a,x,y)

      return
      end

      subroutine bc2d_dscalzmatdvec_add(m,n,scal,a,x,y)
      implicit none
      real *8 scal, x(n)
      complex *16 a(m,n), y(m)
      integer m,n
c     local
      complex *16 xj
      integer i,j

      do j = 1,n
         xj = x(j)*scal
         do i = 1,m
            y(i) = y(i)+a(i,j)*xj
         enddo
      enddo

      return
      end

      subroutine bc2d_dscalzmatzvec_dout_add(m,n,scal,a,x,y)
      implicit none
      real *8 y(m), scal
      complex *16 a(m,n), x(n)
      integer m,n
c     local
      complex *16 xj
      integer i,j

      do j = 1,n
         xj = x(j)*scal
         do i = 1,m
            y(i) = y(i)+a(i,j)*xj
         enddo
      enddo

      return
      end

      subroutine bc2d_dscalzmatvec(m,n,scal,a,x,y)
      implicit none
      real *8 scal
      complex *16 a(m,n), x(n), y(m)
      integer m,n
c     local
      complex *16 xj
      integer i,j

      do i = 1,m
         y(i) = 0
      enddo

      call bc2d_dscalzmatvec_add(m,n,scal,a,x,y)

      return
      end

      subroutine bc2d_dscalzmatvec_add(m,n,scal,a,x,y)
      implicit none
      real *8 scal
      complex *16 a(m,n), x(n), y(m)
      integer m,n
c     local
      complex *16 xj
      integer i,j

      do j = 1,n
         xj = x(j)*scal
         do i = 1,m
            y(i) = y(i)+a(i,j)*xj
         enddo
      enddo

      return
      end

      subroutine bc2d_zmatzvec_dout_add(m,n,a,x,y)
      implicit none
      real *8 y(m)
      complex *16 a(m,n), x(n)
      integer m,n
c     local
      complex *16 xj
      integer i,j

      do j = 1,n
         xj = x(j)
         do i = 1,m
            y(i) = y(i)+a(i,j)*xj
         enddo
      enddo

      return
      end
      
      
      subroutine bc2d_dmatvec(m,n,a,x,y)
      implicit none
      real *8 x(n), a(m,n), y(m)
      integer m,n
c     local
      integer i
      
      do i = 1,m
         y(i) = 0
      enddo

      call bc2d_dmatvec_add(m,n,a,x,y)   
      
      return
      end
      
      subroutine bc2d_dmatvec_add(m,n,a,x,y)
      implicit none
      real *8 x(n),a(m,n), y(m)
      integer m,n
c     local
      integer i, j
      real *8 xj

      do j = 1,n
         xj = x(j)
         do i = 1,m
            y(i) = y(i)+a(i,j)*xj
         enddo
      enddo

      return
      end
      
      subroutine bc2d_zmatdvec(m,n,a,x,y)
      implicit none
      real *8 x(n)
      complex *16 a(m,n), y(m)
      integer m,n
c     local
      integer i
      
      do i = 1,m
         y(i) = 0
      enddo

      call bc2d_zmatdvec_add(m,n,a,x,y)   
      
      return
      end
      
      subroutine bc2d_zmatdvec_add(m,n,a,x,y)
      implicit none
      real *8 x(n)
      complex *16 a(m,n), y(m)
      integer m,n
c     local
      integer i, j
      real *8 xj

      do j = 1,n
         xj = x(j)
         do i = 1,m
            y(i) = y(i)+a(i,j)*xj
         enddo
      enddo

      return
      end
      
      subroutine bc2d_dscalmatvec_add(m,n,scal,a,x,y)
c
c     y <- y + scal*a*x
c      
      implicit none
      real *8 a(m,n), x(n), y(m), scal
      real *8 xj
      integer m,n,i,j
      
      do j = 1,n
         xj = x(j)*scal
         do i = 1,m
            y(i) = y(i)+a(i,j)*xj
         enddo
      enddo

      return
      end

      subroutine bc2d_dvecaddscal(n,alpha,y)
      implicit none
      real *8 y(n), alpha
      integer n, i

      do i = 1,n
         y(i) = y(i)+alpha
      enddo

      return
      end

      subroutine bc2d_zvecaddvec(n,x,y)
      implicit none
      complex *16 :: x(n), y(n)
      integer :: n, i

      do i = 1,n
         y(i) = y(i) + x(i)
      enddo
      
      return
      end

      subroutine bc2d_realaddvec(n,x,y)
      implicit none
      real *8 :: x(2,n), y(n)
      integer n,i

      do i = 1,n
         y(i) = y(i) + x(1,i)
      enddo
      return
      end

      subroutine bc2d_zdiagmatvec_add(n,d,x,y)
      implicit none
      complex *16 :: d(n), x(n), y(n)
      integer n,i

      do i = 1,n
         y(i) = y(i) + d(i)*x(i)
      enddo
      return
      end
