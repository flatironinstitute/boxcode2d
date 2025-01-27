      program test_pbox2drouts
c
c     set a  tree and evaluate the convolution by the box code
c     test against adaptive integration
c      
      implicit real *8 (a-h,o-z)

      
      real *8 :: zll(2), blength, timeinfo(10)
      real *8, allocatable :: pts(:,:,:), fvals(:,:)
      real *8, allocatable :: pot(:,:), grad(:,:,:), hess(:,:,:)
      
      integer, allocatable :: levelbox(:), icolbox(:), irowbox(:)
      integer, allocatable :: iparentbox(:), nblevel(:),icolleagbox(:,:)
      integer, allocatable :: ichildbox(:,:), iboxlev(:)
      integer, allocatable :: istartlev(:), itree(:)

      real *8 :: gradtmp(2)
      
      character :: pttype, polytype

      integer :: ipass(20)

      ntest = 0
      do i = 1,20
         ipass(i) = 0
      enddo
      
      iseed = 1234
      dburn = hkrand(iseed)
      
      call prini(6,13)

      polytype = 'T'
      pttype = 'F'
      norder = 2

      eps = 1d-8
      
      zll(1) = -0.23d0
      zll(2) = -0.17d0
      zll(1) = 0
      zll(2) = 0
      blength = 0.3d0
      blength = 1

      maxboxes = 400000
      maxlevel = 50
      
      allocate(levelbox(maxboxes),icolbox(maxboxes),
     1     irowbox(maxboxes),icolleagbox(9,maxboxes),
     1     iparentbox(maxboxes),ichildbox(4,maxboxes),
     2     nblevel(0:maxlevel),
     2     iboxlev(maxboxes),istartlev(0:maxlevel))      

c     and solve a poisson problem
      
      imode = 3 
      call lrt2d_set(levelbox,icolbox,irowbox,nboxes,nlev,
     1     ichildbox,iparentbox,nblevel,istartlev,iboxlev,imode)
c      nlev = 6
c      call lrt2d_uni(levelbox,icolbox,irowbox,nboxes,nlev, 
c     1     ichildbox, iparentbox, nblevel, istartlev, iboxlev)

      call prinf('nlev *',nlev,1)
      call prinf('nboxes *',nboxes,1)      
      
      litree = 10*nboxes+5*(nlev+1)+100
      allocate(itree(litree))
      ifnbr = 0
      ifcolleag = 0
      ier = 0
      call lrt2d_pack(itree,litree,levelbox,icolbox,irowbox,
     1     nboxes,nlev,ichildbox,iparentbox,nblevel,istartlev,iboxlev,
     2     ifcolleag,icolleagbox,ifnbr,neighbors,nnbrs,ier)
      call prinf('ier, after tree pack*',ier,1)

      npts = norder**2
      allocate(pts(2,npts,nboxes),fvals(npts,nboxes))
      call compositegrid2d_getpts0(norder, pttype,
     1     levelbox, icolbox, irowbox, nboxes, zll, blength,
     2     pts, npts, ier)

      fvals = 0
      call fillfvals(fvals,pts,npts,nboxes,ichildbox)

      ifpot = 1
      ifgrad = 1
      ifhess = 0
      allocate(pot(npts,nboxes),grad(2,npts,nboxes),hess(3,npts,nboxes))

      nleaf = 0
      do i = 1,nboxes
         if (ichildbox(1,i) .gt. 0) nleaf = nleaf+1
      enddo
      
      call prinf('npts *',npts,1)
      call prinf('nleaf *',nleaf,1)
      call prinf('npts total *',npts*nleaf,1)

      call bc2d_time(t0)
      call poissbox2d(eps,litree,itree,blength,
     1     nd,nboxes,npts,norder,pttype,polytype,fvals,
     2     ifpot,pot,ifgrad,grad,ifhess,hess,ier,timeinfo)
      call bc2d_time(t1)

      call prinf('ier *',ier,1)
      call prin2('timeinfo *',timeinfo,5)

      call prin2('points per second *',nleaf*npts/(t1-t0),1)

c     test

      epsabs = 1d-12
      epsrel = 1d-12
      maxpts = 100000
      pi = 4*atan(1d0)

      derrmax = 0
      derrmaxgrad = 0
      
      do i = 1,nboxes
         if (ichildbox(1,i) .le. 0) then
            write(*,*) 'box ', i, nboxes            
            do j = 1,npts
               call getconv(pts(1,j,i),zll,blength,epsabs,epsrel,
     1              maxpts,val,ier)
               call getconvgrad(pts(1,j,i),zll,blength,epsabs,epsrel,
     1              maxpts,gradtmp,ier)

               derr = abs(val-pot(j,i))/max(abs(val),1d0)
               derrgrad1 = abs(grad(1,j,i)-gradtmp(1))/
     1              max(abs(gradtmp(1)),1d0)
               derrgrad2 = abs(grad(2,j,i)-gradtmp(2))/
     1              max(abs(gradtmp(2)),1d0)
               write(*,*) derr, derrgrad1, derrgrad2
c               write(*,*) derrgrad1, derrgrad2
               derrmax = max(derr,derrmax)
               derrmaxgrad = max(derrgrad1,derrmaxgrad)
               derrmaxgrad = max(derrgrad2,derrmaxgrad)
            enddo
c     exit
         endif
      enddo

      write(*,*) 'err max (pot, grad)'
      write(*,*) derrmax, derrmaxgrad

      ntest = ntest+1
      if (derrmax .lt. eps) ipass(ntest) = 1
      
      npass = 0
      do i = 1,ntest
         npass = npass + ipass(i)
      enddo
      
      open(unit=33,file='print_testres.txt',access='append')
      write(33,'(a,i2,a,i2,a)') 'Successfully completed ',npass,
     1   ' out of ',ntest,' tests in pbox2drouts testing suite'
      close(33)
      
      stop
      end


      subroutine fillfvals(fvals,pts,npts,nboxes,ichildbox)
      implicit real *8 (a-h,o-z)
      real *8 fvals(npts,nboxes), pts(2,npts,nboxes)
      integer ichildbox(4,nboxes)
      integer nboxes

      do i = 1,nboxes
         if (ichildbox(1,i) .le. 0) then
            do j = 1,npts
               call fun1(pts(1,j,i),fvals(j,i))
            enddo
         endif
      enddo
      return
      end


      subroutine fun1(z,f)
      implicit real *8 (a-h,o-z)
      real *8 :: z(2),f


      x = z(1)
      y = z(2)

      f = 4*x**2 + 3*y**3 + x - y + 0.2d0
      f = 4*x - y + 0.2d0

c     zero out except for box 6
c      if (x .ge. 0.25 .or. y .ge. 0.25 .or. x .le. 0
c     1     .or. y .le. 0) f = 0


c     zero out except for box 3
c      if (x .ge. 1 .or. y .ge. 0.5 .or. x .le. 0.5
c     1     .or. y .le. 0) f = 0
      return
      end

      subroutine getconvgrad(x0,zll,blength,epsabs,epsrel,maxpts,grad,
     1     ier)
      implicit real *8 (a-h,o-z)
      real *8 :: x0(2), zll(2)
c     local
      real *8 :: a(2), b(2), grad(2)
      real *8, allocatable :: work(:)
      integer restar
      external loggrad

      nw = maxpts*10 + 100
      allocate(work(nw))

      ndim = 2
      numfun = 2
      key = 0
      minpts = min(16,maxpts)

      a(1) = zll(1)
      a(2) = zll(2)
      b(1) = zll(1)+blength
      b(2) = zll(2)+blength

      restar = 0

      call loggradinit(x0)
      
      call dcuhre(ndim,numfun,a,b,minpts,maxpts,loggrad,epsabs,
     1     epsrel,key,nw,restar,grad,abserr,neval,ier,
     2     work)

c      call prinf('dcuhre, ier *',ier,1)
c      call prin2('abserr *',abserr,1)
c      call prinf('neval *',neval,1)
      

      return
      end

      subroutine loggrad(ndim,x,numfun,funvls)
      implicit real *8 (a-h,o-z)
      real *8 :: funvls(*), x(2), x07(2)
      real *8 :: x0(2), ff, r2, diff(2)
      
      save x0

      pi2 = 8*atan(1d0)
      
      call fun1(x,ff)

      diff(1) = x0(1)-x(1)
      diff(2) = x0(2)-x(2)
      r2 = diff(1)**2 + diff(2)**2
      funvls(1) = -ff*diff(1)/pi2/r2
      funvls(2) = -ff*diff(2)/pi2/r2

      return
      
      entry loggradinit(x07)
      x0(1) = x07(1)
      x0(2) = x07(2)

      return
      end

      subroutine getconv(x0,zll,blength,epsabs,epsrel,maxpts,val,ier)
      implicit real *8 (a-h,o-z)
      real *8 :: x0(2), zll(2)
c     local
      real *8 :: a(2), b(2)
      real *8, allocatable :: work(:)
      integer restar
      external logfun

      nw = maxpts*10 + 100
      allocate(work(nw))

      ndim = 2
      numfun = 1
      key = 0
      minpts = min(16,maxpts)

      a(1) = zll(1)
      a(2) = zll(2)
      b(1) = zll(1)+blength
      b(2) = zll(2)+blength

      restar = 0

      call logfuninit(x0)
      
      call dcuhre(ndim,numfun,a,b,minpts,maxpts,logfun,epsabs,
     1     epsrel,key,nw,restar,val,abserr,neval,ier,
     2     work)

c      call prinf('dcuhre, ier *',ier,1)
c      call prin2('abserr *',abserr,1)
c      call prinf('neval *',neval,1)
      

      return
      end

      subroutine logfun(ndim,x,numfun,funvls)
      implicit real *8 (a-h,o-z)
      real *8 :: funvls(*), x(2), x07(2)
      real *8 :: x0(2), ff, r2
      
      save x0

      pi4 = 16*atan(1d0)
      
      call fun1(x,ff)

      r2 = (x0(1)-x(1))**2 + (x0(2)-x(2))**2
      funvls(1) = -ff*log(r2)/pi4

      return
      
      entry logfuninit(x07)
      x0(1) = x07(1)
      x0(2) = x07(2)

      return
      end

