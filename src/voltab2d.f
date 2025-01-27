c
c
c     routines for building volume tables in 2d from reference
c     tables, assuming certain symmetries of the target points
c
c     this file takes the following conventions
c
c     a tensor grid of points is traversed with x on the inner
c     loop and y on the loop above that
c     e.g. for a 2x2 grid, we have the order:
c            (x1,y1), (x2,y1), (x1,y2), (x2,y2)
c      
c     tensor product polynomials are numbered analagously
c        e.g. for full degree at order 3, we would have the order
c           1, x, x^2, y, xy, x^2y, y^2, xy^2, x^2 y^2, 
c        e.g. for total degree, we would have the order
c           1, x, x^2, y, xy, y^2
c


      subroutine bc2d_dreftab2full(norder,pttype,polytype,npt,npoly,
     1     tabcol,tabstob,tabbtos,gradtabcol,gradtabstob,
     2     gradtabbtos)      
c     
c     uses basic permutation info to extract the requested
c     table from the reference set
c
c     input
c
c     ndeg - integer, degree of polynomial basis
c     type - character, basis type. 'T' total degree, 'F' full degree      
c     npt - integer, number of points in grid
c     npol - integer, number of polynomials in basis
c
c     output
c      
c     tabout - real *8 (npt,npol), the desired table after the
c                necessary permutations and transformations
      implicit none
      integer norder, npt, npoly
      real *8 :: tabcol(npt,npoly,*), tabstob(npt,npoly,*),
     1     tabbtos(npt,npoly,*)
      real *8 :: gradtabcol(2,npt,npoly,*), gradtabstob(2,npt,npoly,*),
     1     gradtabbtos(2,npt,npoly,*)
      character pttype, polytype
c     local
      integer n, j, k, ic
      integer, allocatable :: ipperm(:), icperm(:), icsign(:)
      integer :: iref(12), idimp(2,12), iflip(2,12)
      integer :: igperm(2), igsign(2)
      
      allocate(ipperm(npt),icperm(npoly),icsign(npoly))

c     expand colleague tables
c     WARNING: the order of the outer loop matters
      
      call loadsyms2dc(iref,idimp,iflip)
      do ic = 9,3,-1
         call buildperm2d(idimp(1,ic),iflip(1,ic),norder,
     1        pttype,polytype,ipperm,icperm,icsign,igperm,igsign)

         do j = 1,npoly
            do k = 1,npt
               tabcol(ipperm(k),j,ic) =
     1              icsign(j)*tabcol(k,icperm(j),iref(ic))
               gradtabcol(igperm(1),ipperm(k),j,ic) =
     1              igsign(1)*icsign(j)*
     2              gradtabcol(1,k,icperm(j),iref(ic))
               gradtabcol(igperm(2),ipperm(k),j,ic) =
     1              igsign(2)*icsign(j)*
     2              gradtabcol(2,k,icperm(j),iref(ic))
            enddo
         enddo
      enddo

c     expand stob and btos tables
      
      call loadsyms2dstob(iref,idimp,iflip)
      do ic = 3,12
         call buildperm2d(idimp(1,ic),iflip(1,ic),norder,
     1        pttype,polytype,ipperm,icperm,icsign,igperm,igsign)

         do j = 1,npoly
            do k = 1,npt
               tabstob(ipperm(k),j,ic) =
     1              icsign(j)*tabstob(k,icperm(j),iref(ic))
               gradtabstob(igperm(1),ipperm(k),j,ic) =
     1              igsign(1)*icsign(j)*
     2              gradtabstob(1,k,icperm(j),iref(ic))
               gradtabstob(igperm(2),ipperm(k),j,ic) =
     1              igsign(2)*icsign(j)*
     2              gradtabstob(2,k,icperm(j),iref(ic))
            enddo
         enddo
      enddo

      call loadsyms2dbtos(iref,idimp,iflip)
      do ic = 3,12
         call buildperm2d(idimp(1,ic),iflip(1,ic),norder,
     1        pttype,polytype,ipperm,icperm,icsign,igperm,igsign)

         do j = 1,npoly
            do k = 1,npt
               tabbtos(ipperm(k),j,ic) =
     1              icsign(j)*tabbtos(k,icperm(j),iref(ic))
               gradtabbtos(igperm(1),ipperm(k),j,ic) =
     1              igsign(1)*icsign(j)*
     2              gradtabbtos(1,k,icperm(j),iref(ic))
               gradtabbtos(igperm(2),ipperm(k),j,ic) =
     1              igsign(2)*icsign(j)*
     2              gradtabbtos(2,k,icperm(j),iref(ic))
            enddo
         enddo
      enddo

      return
      end


      
      subroutine buildperm2d(idimp,iflip,nq,pttype,polytype,
     1     ipperm,icperm,icsign,igperm,igsign)
c
c
c     combines dimension permutation and flip information into
c     index permutation (for target points) and coefficient
c     permutation (for source densities) information. also returns
c     overall sign flips per coefficient of a source density
c            
      implicit real*8 (a-h,o-z)
      dimension idimp(2), iflip(2), ipperm(*), icperm(*)
      dimension icsign(*), igperm(2), igsign(2)
      dimension indmap(nq,2), iskip(2), ilist(nq,nq)
      dimension itemp(2)
      character polytype, pttype

      iskip(1) = 1
      iskip(2) = nq

      igperm(1) = idimp(1)
      igperm(2) = idimp(2)

      igsign(1) = iflip(idimp(1))
      igsign(2) = iflip(idimp(2))
      
      do ii = 1,2
         do jj = 1,nq
            if (iflip(ii) .eq. -1) then
               indmap(jj,ii) = nq-jj+1
            else
               indmap(jj,ii) = jj
            endif
         enddo
      enddo

c     only point permutation implemented is full tensor grid
      
      do iy = 1,nq
         jy = indmap(iy,idimp(2))
         do ix = 1,nq
            jx = indmap(ix,idimp(1))
            ind1 = (iy-1)*nq + ix
            ind2 =  (jy-1)*iskip(idimp(2)) +
     1           (jx-1)*iskip(idimp(1)) + 1
            ipperm(ind1) = ind2
         enddo
      enddo

      if (polytype .eq. 't' .or. polytype .eq. 'T') then
c     total degree order
         ii = 0
         do iy = 1,nq
            do ix = 1,nq+1-iy
               ii = ii+1
               ilist(ix,iy) = ii
            enddo
         enddo

         ii = 0
         do iy = 1,nq
            do ix = 1,nq+1-iy
               ii = ii+1
               itemp(1) = ix
               itemp(2) = iy
               jx = itemp(idimp(1))
               jy = itemp(idimp(2))
               jj = ilist(jx,jy)
               
               icperm(ii) = jj
               icsign(ii) = (iflip(idimp(1)))**(jx-1)*
     2              (iflip(idimp(2)))**(jy-1)
            enddo
         enddo
         
      else if (polytype .eq. 'f' .or. polytype .eq. 'F') then
c     full tensor product
         ii = 0
         do iy = 1,nq
            do ix = 1,nq
               ii = ii+1
               ilist(ix,iy) = ii
            enddo
         enddo

         ii = 0
         do iy = 1,nq
            do ix = 1,nq
               ii = ii+1
               itemp(1) = ix
               itemp(2) = iy
               jx = itemp(idimp(1))
               jy = itemp(idimp(2))
               jj = ilist(jx,jy)
               
               icperm(ii) = jj
               icsign(ii) = (iflip(idimp(1)))**(jx-1)*
     2              (iflip(idimp(2)))**(jy-1)
            enddo
         enddo
         
      endif
      return
      end
      
c
      subroutine alltargs2d_grid(grid0,ngrid,bs,xyc,wc,wbtos,wstob)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     generate reference points for all targets in colleague,
c     big to small, and small to big interactions
c
      
      implicit real *8 (a-h,o-z)
      dimension wc(2,ngrid,9), wbtos(2,ngrid,12), wstob(2,ngrid,12)
      dimension xyc(2)
      dimension grid0(2,ngrid)
      real *8, allocatable :: grid(:,:)

      allocate(grid(2,ngrid))

      do i=1,ngrid
        grid(1,i) = (grid0(1,i) + 1)*bs/2
        grid(2,i) = (grid0(2,i) + 1)*bs/2
      enddo

c     lowest corner of cube (each coordinate is smallest)
      
      xc = xyc(1)
      yc = xyc(2)

c     get corresponding meshgrid

c     colleagues are straightforward

      ind = 0
      do iy = -1,1
         do ix = -1,1
            ind = ind+1
            do i = 1,ngrid
               wc(1,i,ind) = xc + ix*bs + grid(1,i)
               wc(2,i,ind) = yc + iy*bs + grid(2,i)
            enddo
         enddo
      enddo


c     stob and btos are analogous but a
c     little more complicated

      bsh = bs/2

      ind = 0
      do iy = -1,2
         if (iy .eq. -1 .or. iy .eq. 2) then
c     all x positions possible
            do ix = -1,2
               ind = ind+1
               do i = 1,ngrid
                  wbtos(1,i,ind) = xc + ix*bsh + grid(1,i)/2
                  wbtos(2,i,ind) = yc + iy*bsh + grid(2,i)/2
                  
                  wstob(1,i,ind) = xc + (ix-1)*bs + grid(1,i)*2
                  wstob(2,i,ind) = yc + (iy-1)*bs + grid(2,i)*2
               enddo
            enddo
         else
c     only outer x positions possible
            do ix = -1,2,3
               ind = ind+1
               do i = 1,ngrid
                  wbtos(1,i,ind) = xc + ix*bsh + grid(1,i)/2
                  wbtos(2,i,ind) = yc + iy*bsh + grid(2,i)/2

                  wstob(1,i,ind) = xc + (ix-1)*bs + grid(1,i)*2
                  wstob(2,i,ind) = yc + (iy-1)*bs + grid(2,i)*2
               enddo
            enddo
         endif
      enddo
         
      return
      end
c

      subroutine tensrefpts2d_grid(grid0,ngrid,bs,xyc,wc,wbtos,wstob)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     generate reference points for the limited subset of
c     points in colleague, big-to-small, and small-to-big
c     interactions which can be used to obtain the other
c     interactions by symmetries, where points on unit box
c     are specified
c
      
      implicit real *8 (a-h,o-z)
      dimension grid0(2,ngrid)
      dimension wc(2,ngrid,3), wbtos(2,ngrid,2), wstob(2,ngrid,2)
      dimension xshift(3), yshift(3), xyc(2)
      allocatable :: grid(:,:)

c     lowest corner of cube (each coordinate is smallest)
      
      xc = xyc(1)
      yc = xyc(2)

c     get corresponding meshgrid

      allocate(grid(2,ngrid))
      do i = 1,ngrid
         grid(1,i) = (grid0(1,i)+1)*bs/2
         grid(2,i) = (grid0(2,i)+1)*bs/2
      enddo

      xshift(1) = -1
      xshift(2) = 0
      xshift(3) = 0
      yshift(1) = -1
      yshift(2) = -1
      yshift(3) = 0

      bsh = bs/2
      
      do ii = 1,3
         do i = 1,ngrid
            wc(1,i,ii) = xc + xshift(ii)*bs + grid(1,i)
            wc(2,i,ii) = yc + yshift(ii)*bs + grid(2,i)

            if (ii .lt. 3) then
               wbtos(1,i,ii) = xc + xshift(ii)*bsh + grid(1,i)/2
               wbtos(2,i,ii) = yc + yshift(ii)*bsh + grid(2,i)/2

               wstob(1,i,ii) = xc + (xshift(ii)-1)*bs + grid(1,i)*2
               wstob(2,i,ii) = yc + (yshift(ii)-1)*bs + grid(2,i)*2
            endif
         enddo
      enddo

      return
      end

