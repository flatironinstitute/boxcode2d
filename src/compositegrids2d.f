      subroutine compositegrid2d_getpts0(norder, pttype,
     1     levelbox, icolbox, irowbox, nboxes, zll, blength,
     2     pts, npts, ier)
c****************************************************************
c     get equispaced points of the requested order,
c     scaled to each box in the tree hierarchy
C
C     INPUT:
C
c     NORDER - integer. the order of grid points to return
c     PTTYPE - character. point type
c          PTTYPE .eq. 'F' -> full tensor product Legendre 
c                             grid (npts = norder**2 points)
C     LEVELBOX is an array determining the level of each box
C     NBOXES is the total number of boxes
C     ICOLBOX denotes the column of each box
C     IROWBOX denotes the row of each box
c     ZLL - real *8(2) array. lower left corner of root box
c     BLENGTH - real *8. size of root box
c     NPTS - second dimension of pts array. 
C
C     OUTPUT:
C
c     PTS - real *8(2,NPTS,NBOXES) array. the
c     points of the requested order and type,
c     scaled to each box.
c     IER - integer. error flag
c             IER .eq. 1 -> npts not consistent with requested
c                     norder and pttype 
C             IER .eq. 2 -> unavailable pttype or norder
C*****************************************************************
      implicit real *8 (a-h,o-z)
c-----Global variables
      integer norder, npts, nboxes
      integer icolbox(*), irowbox(*), levelbox(*)
      real *8 pts(2,npts,*), zll(2), blength
      character pttype 
c-----Local variables
      real *8, allocatable :: x(:,:)
      character polytypetmp

      polytypetmp = 'T'
c

      if (pttype .eq. 'F' .or. pttype .eq. 'f') then
         npts2 = norder**2
         if (npts .ne. npts2) then
            ier = 1
            return
         endif
         itype=0
         allocate(x(2,npts))
         call legetens_exps_2d(itype,norder,polytypetmp,x,
     1        u,1,v,1,w)
      else
         ier = 2
         return
      endif

c$OMP PARALLEL DO DEFAULT(shared) PRIVATE(i,xlength,scal)
c$OMP$ PRIVATE(xshift,yshift,j)      
      do i = 1,nboxes
         xlength = blength/(2d0**levelbox(i))
         scal = xlength/2
         xshift  =  zll(1)+dble(icolbox(i) - 1)*xlength
         yshift  =  zll(2)+dble(irowbox(i) - 1)*xlength
         do j = 1,npts
            pts(1,j,i) = xshift + (x(1,j)+1)*scal
            pts(2,j,i) = yshift + (x(2,j)+1)*scal
         enddo
      enddo
c$OMP END PARALLEL DO
      
      return
      end

