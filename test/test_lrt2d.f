      program test_lrt2d
c
c
c     unit tests for level-restricted tree
c     
c      
      implicit real *8 (a-h,o-z)

      
      real *8, allocatable :: pts(:,:), pts_sort(:,:)
      real *8 :: zll(2), blength
      
      integer, allocatable :: levelbox(:), icolbox(:), irowbox(:)
      integer, allocatable :: iparentbox(:), nblevel(:),icolleagbox(:,:)
      integer, allocatable :: ichildbox(:,:), iboxlev(:)
      integer, allocatable :: istartlev(:), iptsladder(:,:)
      integer, allocatable :: isort(:), itemparray(:)

      integer ierrs(10)

      integer ipass(20)

      ntest = 0
      do i = 1,20
         ipass(i) = 0
      enddo
      
      iseed = 1234
      dburn = hkrand(iseed)
      
      call prini(6,13)

      zll(1) = 0d0
      zll(2) = 0d0
      blength = 1d0

      npts = 10000

c     set up points
      
      allocate(pts(2,npts),pts_sort(2,npts),isort(npts))

      do i = 1,npts
         pts(1,i) = zll(1) + blength*sqrt(hkrand(0))
         pts(2,i) = zll(2) + blength*sqrt(hkrand(0))
      enddo
      
      nstart = 20*npts
      maxiter = 5
      maxnodes = 40

      nlev = -1
      nboxes = -1
      ier = 0
      
      call lrt2d_mktpts_query(pts,pts_sort,isort,npts,
     1     maxnodes,zll,blength,nstart,maxiter,nboxes,nlev,ier)

      call prinf('nboxes *',nboxes,1)
      call prinf('nlev *',nlev,1)

      if (ier .ne. 0) then
         write(*,*) 'query failed, abort test'
         stop
      endif
      
      allocate(levelbox(nboxes),icolbox(nboxes),irowbox(nboxes),
     1     iparentbox(nboxes),ichildbox(4,nboxes),nblevel(0:nlev),
     2     iboxlev(nboxes),istartlev(0:nlev),itemparray(nboxes),
     3     iptsladder(2,nboxes),icolleagbox(9,nboxes))

      maxboxes = nboxes
      maxlevel = nlev
      
      call lrt2d_mktpts(levelbox, icolbox, irowbox, nboxes, nlev,
     1     iparentbox, ichildbox, nblevel, iboxlev, istartlev,
     2     maxboxes, itemparray, maxlevel, pts, pts_sort, isort,
     3     iptsladder, npts, maxnodes, zll, blength, ier)

      call prinf('nblevel *',nblevel,nlev+1)

      call test_pts_sort(levelbox,icolbox,irowbox,nboxes,nlev,
     1     iparentbox, ichildbox, nblevel, iboxlev, istartlev,
     2     pts, pts_sort, isort, iptsladder, npts, maxnodes,
     3     zll, blength, ier)

      call prinf('ier *', ier,1)

      ntest = ntest + 1
      if (ier .eq. 0) ipass(ntest) = 1

      call lrt2d_testtree(levelbox,iparentbox,ichildbox,icolbox, 
     1     irowbox,nboxes,nlev,nblevel,iboxlev,istartlev,
     2     ierrs)

      call lrt2d_mkcolls(icolbox,irowbox,icolleagbox,nboxes,nlev,
     2      iparentbox,ichildbox,nblevel,iboxlev,istartlev)
      

      call prinf('ierrs *',ierrs,4)

      ntest = ntest + 1
      if (ierrs(1) .eq. 0 .and. ierrs(2) .eq. 0 .and. ierrs(3) .eq. 0
     1     .and. ierrs(4) .eq. 0) ipass(ntest) = 1

c      iw=1
c      call lrt2d_plotleaves(iw,zll,blength,levelbox,ichildbox,
c     1     icolbox,irowbox,nboxes)

c      ic = 366
c      ibox = iparentbox(ic)
c      iw = 2
c      call lrt2d_plotlists(iw,ibox,zll,blength,levelbox,icolleagbox,
c     1     ichildbox,icolbox,irowbox,nboxes)
      

      if (maxboxes .lt. 10000 .or. maxlevel .lt. 20) then
         maxboxes = 10000
         maxlevel = 20
         
         deallocate(levelbox,icolbox,irowbox,iparentbox,ichildbox,
     2        nblevel,iboxlev,istartlev)
         
         allocate(levelbox(maxboxes),icolbox(maxboxes),
     1        irowbox(maxboxes),
     1        iparentbox(maxboxes),ichildbox(4,maxboxes),
     2        nblevel(0:maxlevel),
     2        iboxlev(maxboxes),istartlev(0:maxlevel))      
      endif

      do imode = 2,4
         write(*,*) 'testing ', imode, ' level adaptive tree'
         call lrt2d_set(levelbox,icolbox,irowbox,nboxes,nlev,
     1        ichildbox,iparentbox,nblevel,istartlev,iboxlev,imode)
      
         call lrt2d_testtree(levelbox,iparentbox,ichildbox,icolbox, 
     1        irowbox,nboxes,nlev,nblevel,iboxlev,istartlev,
     2        ierrs)

         ntest = ntest + 1
         if (ierrs(1) .eq. 0 .and. ierrs(2) .eq. 0 .and. ierrs(3) .eq. 0
     1        .and. ierrs(4) .eq. 0) ipass(ntest) = 1
         
         call prinf('ierrs *',ierrs,4)
      enddo


      nlev = 5
      write(*,*) 'testing ',nlev, ' level uniform tree'
      call lrt2d_uni(levelbox,icolbox,irowbox,nboxes,nlev, 
     1     ichildbox, iparentbox, nblevel, istartlev, iboxlev)

      call lrt2d_testtree(levelbox,iparentbox,ichildbox,icolbox, 
     1     irowbox,nboxes,nlev,nblevel,iboxlev,istartlev,
     2     ierrs)
      call prinf('ierrs *',ierrs,4)

      ntest = ntest + 1
      if (ierrs(1) .eq. 0 .and. ierrs(2) .eq. 0 .and. ierrs(3) .eq. 0
     1     .and. ierrs(4) .eq. 0) ipass(ntest) = 1
      

      npass = 0
      do i = 1,ntest
         npass = npass + ipass(i)
      enddo
      
      open(unit=33,file='print_testres.txt',access='append')
      write(33,'(a,i2,a,i2,a)') 'Successfully completed ',npass,
     1   ' out of ',ntest,' tests in lrt2d testing suite'
      close(33)
      
      
      stop
      end
      

      subroutine test_pts_sort(levelbox,icolbox,irowbox,nboxes,nlev,
     1     iparentbox, ichildbox, nblevel, iboxlev, istartlev,
     2     pts, pts_sort, isort, iptsladder, npts, maxnodes,
     3     zll, blength, ier)
c
c     check tree based on sorting points into boxes.
c
c     checks:
c     (1) sorting is correct
c     (2) all points accounted for
c     (3) points are in the boxes it says they are
c     (4) maxnodes restriction actually met
c      
      implicit real *8 (a-h,o-z)
      integer :: levelbox(*), icolbox(*), irowbox(*), nboxes
      integer :: iparentbox(*), ichildbox(4,*), nblevel(0:*)
      integer :: iboxlev(*), istartlev(0:*), isort(*)
      integer :: iptsladder(2,*)
      real *8 :: pts(2,*), pts_sort(2,*), zll(2), blength
      integer nlev, npts, maxnodes, ier
c
      integer, allocatable :: itemp(:)

      tol = 1d-15
      ier = 0

      allocate(itemp(npts))
      
c     (1) check sorting

      do i = 1,npts
         if (pts(1,isort(i)) .ne. pts_sort(1,i) .or.
     1        pts(2,isort(i)) .ne. pts_sort(2,i)) then
            ier = 1
            return
         endif
      enddo

      do i = 1,npts
         itemp(i) = 0
      enddo

      do i = 1,npts
         itemp(isort(i)) = itemp(isort(i)) + 1
      enddo

      do i = 1,npts
         if (itemp(i) .ne. 1) then
            ier = 2
            return
         endif
      enddo

c     (2) all points accounted for?

      do i = 1,npts
         itemp(i) = 0
      enddo
      
      do i = 1,nboxes
         istart = iptsladder(1,i)
         iend = iptsladder(2,i)
         if (ichildbox(1,i) .le. 0) then
            do j = istart,iend
               itemp(j) = itemp(j)+1
            enddo
         endif
      enddo

      do i = 1,npts
         if (itemp(i) .ne. 1) then
            ier = 3
            return
         endif
      enddo

c     (3) check points landed in right boxes

      do i = 1,nboxes
         irow = irowbox(i)
         icol = icolbox(i)
         lev = levelbox(i)
         scal = 2.0d0**lev
         xleft = zll(1) + blength*(icol-1)/scal - blength*tol
         xright = zll(1) + blength*(icol)/scal + blength*tol
         ybot = zll(2) + blength*(irow-1)/scal - blength*tol
         ytop = zll(2) + blength*(irow)/scal + blength*tol
         istart = iptsladder(1,i)
         iend = iptsladder(2,i)
         do j = istart,iend
            xx = pts_sort(1,j)
            yy = pts_sort(2,j)
            if (xx .lt. xleft .or. xx .gt. xright .or.
     1           yy .lt. ybot .or. yy .gt. ytop) then
               ier = 4
               return
            endif
         enddo
      enddo
      
c     (4) maxnodes restriction

      do i = 1,nboxes
         istart = iptsladder(1,i)
         iend = iptsladder(2,i)
         if (ichildbox(1,i) .le. 0) then
            if (iend - istart + 1 .gt. maxnodes) then
               ier = 5
               return
            endif
         endif
      enddo
      
      return
      end
