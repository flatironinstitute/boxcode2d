C
c
c     TODO
c     - add (nd,) dimensioning
C*****
******


************************************************************
      subroutine poissbox2d(eps,litree,itree,blength,
     1     nd,nboxes,npts,norder,pttype,polytype,fval,
     2     ifpot,pot,ifgrad,grad,ifhess,hess,ier,timeinfo)
c***********************************************************************
c
c     Computes the volume potential:
c
c     pot(x) = -1/(2*pi) int_B log|x-y| f(y) dV(y)
c
c     Where B is a box of side length blength and f is the polynomial
c     interpolant of f on a composite grid. Specifically B is given as
c     the union of the leaf boxes of a level-restricted tree
c     and f is specified by its values on a scaled grid of points on
c     each leaf box. The potential, pot, is evaluated that the same
c     points where f is specified.
c
c     input
c
c     eps - real *8, requested precision (req 1d-3 .le. eps .le. 1d-15)
c     litree - integer, length of itree array
c     itree - integer(litree), packed array storage of 2d level
c        restricted tree.
c     blength - real *8, side length of root box in the tree
c     nd - integer, UNUSED FOR NOW
c     nboxes - integer, number of boxes in tree
c     npts - integer, number of points per leaf box
c     norder - integer, order of polynomial basis (one more than degree,
c         typically)
c     pttype - character, type of points to use on each leaf box
c       pttype (IGNORED FOR NOW, ALWAYS 'F'):
c         'F' -> full tensor product norder x norder Legendre grid
c     polytype - character, polynomial basis to use is the span
c       of { x^m y^n where deg(x^m y^n) < norder }. degree function 
c       depends on polytype (IGNORED FOR NOW, ALWAYS 'T'):
c         'T' -> total degree, i.e. deg(x^m y^n) = m+n
c     fval - real *8(npts,nboxes) array of point values of function f
c        fval(j,i) - if the ith box is a leaf box, this is the value of f
c        at the jth grid point scaled to box i. If the ith box is not
c        a leaf, this value is ignored.
c     ifpot - integer, flag, ifpot = 1 -> compute potential 
c     ifgrad - integer, flag, ifgrad = 1 -> compute gradient of potential
c     ifhess - integer, flag, ifhess = 1 -> DOES NOTHING FOR NOW
c     
c     output
c
c     pot - real *8(npts,nboxes) potential (if requested) at the same
c        points where f is specified 
c     grad - real *8(2,npts,nboxes) potential (if requested) at the same
c        points where f is specified 
c     ier - error flag
c        ier = 1, failed to determine nterms for fmm. bad value of eps?
c        ier = 2, npts incompatible with pttype and norder 
c        ier = 3, pttype not available 
c        ier = 4, nboxes not compatible with itree structure 
c        ier = 5, tables not found for requested order, pttype and polytype
c     timeinfo - real *8(5)
c        timeinfo(1) - seconds spent on precomputes 
c        timeinfo(2) - seconds spent on upward pass (mp2mp)
c        timeinfo(3) - seconds spent on sideways work (mp2loc)
c        timeinfo(4) - seconds spent on downward pass (loc2loc)
c        timeinfo(5) - seconds spent on direct
c      
C***********************************************************************
      implicit real *8 (a-h,o-z)
C-----Global variables
      integer litree, itree(litree), nd, nboxes, npts, ipttype
      integer norder
      real *8 blength, fval(npts,nboxes), pot(npts,*)
      real *8 grad(2,npts,*), hess(3,npts,*), timeinfo(*)
      character :: pttype, polytype
c-----Local variables
      integer, allocatable :: icolleagbox(:,:), localonoff(:)
      real *8, allocatable :: xnodes(:), wnodes(:), ratio(:),
     1     choose(:,:), comp(:,:)
      complex *16, allocatable :: poly2mpmat(:,:), zs(:,:,:,:),
     1     poly2pw_btosfar(:,:,:)

      call bc2d_time(time0)

c     nd not implemented 
      nd = 1
      
      ier=0 
      
c      pttype = 'F'
c      polytype = 'T'
      ndeg = norder-1
      call legetens_npol_2d(ndeg,polytype,npoly)

      ier1 = 0
      call l2dterms(eps,nterms,ier1)
      if (ier1 .ne. 0) then
         ier = 1
         return
      endif

      if (pttype .eq. 'F') then
         if (npts .ne. norder**2) then
            ier = 2
            return
         endif
      else
         ier = 3
         return
      endif
      
      call lwtsexp2b_terms(eps,nnodes)

      call lrt2d_unpack1(itree,nboxes1,nlev,ifcolleag1,ifnbr1,
     1     ilevelbox,iicolbox,iirowbox,iichildbox,iiparentbox,
     2     iiboxlev,inblevel,iistartlev,iicolleagbox,ineighbors,
     3     innbrs)

      if (nboxes .ne. nboxes1) then
         ier = 4
         return
      endif

      allocate(icolleagbox(9,nboxes))
      call lrt2d_mkcolls(itree(iicolbox),itree(iirowbox),icolleagbox,
     1     nboxes,nlev,itree(iiparentbox),itree(iichildbox),
     2     itree(inblevel),itree(iiboxlev),itree(iistartlev))


cccccccccccccccccccccccccccccc
c      
c     precomputes
c
cccccccccccccccccccccccccccccc

      call bc2d_time(time00)

      allocate(poly2mpmat(0:nterms,npoly))
      nleg = norder + 14
      call pbox2d_genpoly2mpmat(poly2mpmat,nterms,npoly,polytype,ndeg,
     1     nleg)

c     planewave and translation precomputations
c     
      allocate(xnodes(nnodes),wnodes(nnodes),choose(2*nterms,2*nterms),
     1     ratio(nnodes),comp(0:nterms,nnodes))
      call lwtsexp2b(nnodes,xnodes,wnodes,errnodes)
      call l2dpw_precomp(comp,nterms,nnodes,xnodes,wnodes,ratio,choose,
     1     xsum)

      allocate(poly2pw_btosfar(0:nnodes,4,npoly))
      call pbox2d_genpoly2pwmat_btosfar(norder,npoly,polytype,
     1     nnodes,wnodes,xnodes,poly2pw_btosfar)

c     precompute matrix from coefficients to multipoles
      
c     precompute planewave shifts
      allocate(zs(0:nnodes,-3:3,-3:3,0:4))
      call l2dpw_mkshifts(xnodes,nnodes,zs)

      allocate(localonoff(nboxes))
      do i = 1,nboxes
         localonoff(i) = 1
      enddo

      call bc2d_time(time11)

      timeinfo(1) = time11-time00

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     main work
c


      call poissbox2d0(nboxes,itree(ilevelbox),itree(iicolbox),
     1     itree(iirowbox),itree(iichildbox),itree(iiparentbox),
     2     itree(iiboxlev),nlev,itree(inblevel),itree(iistartlev),
     2     icolleagbox,blength,localonoff,npts,norder,npoly,
     3     pttype,polytype,fval,nterms,
     3     poly2mpmat,nnodes,xnodes,wnodes,choose,comp,xsum,ratio,zs,
     4     poly2pw_btosfar,ifpot,pot,ifgrad,grad,ifhess,
     5     hess,timeinfo,ier)

      
      return
      end

      
      subroutine poissbox2d0(nboxes,levelbox,icolbox,irowbox,
     1     ichildbox,iparentbox,iboxlev,nlev,nblevel,istartlev,
     2     icolleagbox,blength,localonoff,npts,norder,npoly,pttype,
     3     polytype,fval,nterms,poly2mpmat,nnodes,xnodes,wnodes,choose,
     4     comp,xsum,ratio,zs,poly2pw_btosfar,ifpot,pot,ifgrad,grad,
     5     ifhess,hess,timeinfo,ier)
      implicit real *8 (a-h,o-z)
      complex *16 :: poly2mpmat(0:nterms,*)
      integer levelbox(*), icolbox(*), irowbox(*), ichildbox(4,*)
      integer iparentbox(*), iboxlev(*), icolleagbox(9,*), nblevel(0:*)
      integer istartlev(0:*), localonoff(*)
      real *8 :: blength, wnodes(*), comp(0:nterms,nnodes),xsum,ratio(*)
      real *8 :: fval(npts,nboxes), choose(2*nterms,*), xnodes(*)
      real *8 :: pot(npts,*), grad(2,npts,*), hess(3,npts,*),timeinfo(*)
      character :: pttype, polytype
      complex *16 :: zs(0:nnodes,-3:3,-3:3,0:4),
     1     poly2pw_btosfar(0:nnodes,4,*)
c     local
      real *8, allocatable, dimension(:,:,:) :: tabcol, tabstob, tabbtos
      real *8, allocatable, dimension(:,:,:,:) :: gradtabcol,
     1     gradtabstob, gradtabbtos
      complex *16, allocatable :: mpole(:,:), expnesw(:,:,:),
     1     expneswbig(:,:,:), locexp(:,:), tmpexp(:),tmpvals(:)
      complex *16, allocatable :: mapslocpot(:,:), mapslocgrad(:,:),
     1     mapslochess(:,:)
      complex *16, allocatable :: mapstobpot(:,:,:), mapstobgrad(:,:,:),
     1     mapstobhess(:,:,:)
      integer, allocatable :: iflagnesw(:,:)
      real *8, allocatable :: pt2polymat(:,:), pts(:,:), wts(:)
      real *8, allocatable :: coeffs(:,:), xlengths(:), dlogxlengths(:)

c     load and expand tables

      allocate(tabcol(npts,npoly,9),tabstob(npts,npoly,12),
     1     tabbtos(npts,npoly,12))
      allocate(xlengths(0:(nlev+1)),dlogxlengths(0:(nlev+1)))
      allocate(gradtabcol(2,npts,npoly,9),gradtabstob(2,npts,npoly,12),
     1     gradtabbtos(2,npts,npoly,12))
      
      call pbox2dreftabs(norder,pttype,polytype,npts,npoly,
     1     tabcol,tabstob,tabbtos,gradtabcol,gradtabstob,
     2     gradtabbtos,ier1)

      if (ier1 .ne. 0) then
         ier = 5
         return
      endif

      call bc2d_dreftab2full(norder,pttype,polytype,npts,npoly,
     1     tabcol,tabstob,tabbtos,gradtabcol,gradtabstob,
     2     gradtabbtos)      

      
c     form coefficients
      
      allocate(pt2polymat(npoly,npts),pts(2,npts),wts(npts))
      if (pttype .eq. 'F') then
         kind = 3
         call legetens_exps_2d(kind,norder,polytype,pts,
     1        pt2polymat,npoly,dummy,1,wts)
      endif

      allocate(coeffs(npoly,nboxes))
      call bc2d_dformcoefs(nboxes,ichildbox,npts,npoly,pt2polymat,
     1     fval,coeffs)

c     length at each level

      xlengths(0) = blength
      do i = 1,(nlev+1)
         xlengths(i) = xlengths(i-1)/2
      enddo
      do i = 1,(nlev+1)
         dlogxlengths(i) = log(xlengths(i))
      enddo

c     local evaluation tables

      allocate(mapslocpot(npts,0:nterms),mapslocgrad(npts,0:nterms),
     1     mapslochess(npts,0:nterms))

      call pbox2d_locevaltabs(npts,pts,nterms,mapslocpot,
     1     mapslocgrad,mapslochess,ier1)

c     precompute eval exponentials for stobfar
      allocate(mapstobpot(npts,0:nnodes,4),mapstobgrad(npts,0:nnodes,4),
     1     mapstobhess(npts,0:nnodes,4))
      call pbox2d_stobfartabs(npts,pts,nnodes,xnodes,
     1     mapstobpot,mapstobgrad,mapstobhess)
      
      
c     initialize expansions
      
      allocate(mpole(0:nterms,nboxes),expnesw(0:nnodes,4,nboxes),
     1     expneswbig(0:nnodes,4,nboxes),iflagnesw(4,nboxes),
     2     locexp(0:nterms,nboxes),tmpexp(0:nterms),tmpvals(npts))
      
C$OMP PARALLEL DO PRIVATE(i,j) IF(nboxes*nterms .gt. 100)
C$OMP& SCHEDULE(static)
      do i = 1, nboxes
         do j = 0, nterms
            mpole(j,i) = 0 
            locexp(j,i) = 0 
         enddo
         do j = 0, nnodes
            expnesw(j,1,i) = 0
            expnesw(j,2,i) = 0
            expnesw(j,3,i) = 0
            expnesw(j,4,i) = 0
            expneswbig(j,1,i) = 0
            expneswbig(j,2,i) = 0
            expneswbig(j,3,i) = 0
            expneswbig(j,4,i) = 0
         enddo
         iflagnesw(1,i) = 0
         iflagnesw(2,i) = 0
         iflagnesw(3,i) = 0
         iflagnesw(4,i) = 0
         if (ifpot .eq. 1) then
            do j = 1,npts
               pot(j,i) = 0
            enddo
         endif
         if (ifgrad .eq. 1) then
            do j = 1,npts
               grad(1,j,i) = 0
               grad(2,j,i) = 0
            enddo
         endif
         if (ifhess .eq. 1) then
            do j = 1,npts
               hess(1,j,i) = 0
               hess(2,j,i) = 0
               hess(3,j,i) = 0
            enddo
         endif
      enddo
C$OMP END PARALLEL DO


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     upward pass (formmp or density-to-outgoing)
      
      ntermsp1 = nterms+1

      call bc2d_time(t00)

      do i = nlev, 2, -1
         istart = istartlev(i)
         iend = istart + nblevel(i) - 1
         xlength = xlengths(i)
         scal = (xlength/2)**2
C$OMP PARALLEL DO PRIVATE(ii,j,ic1,ic2,ic3,ic4)
C$OMP& PRIVATE(b)
C$OMP& IF(iend-istart .gt. 100)
C$OMP& SCHEDULE(static,1)
         do ii = istart,iend
            j = iboxlev(ii)
            if(ichildbox(1,j) .le. 0)then
               call bc2d_dscalzmatdvec(ntermsp1,npoly,scal,poly2mpmat,
     1              coeffs(1,j),mpole(0,j))
            else
               ic1 = ichildbox(1,j)
               ic2 = ichildbox(2,j)
               ic3 = ichildbox(3,j)
               ic4 = ichildbox(4,j)
               call l2dchildpar(mpole(0,j),mpole(0,ic1),mpole(0,ic2),
     1              mpole(0,ic3),mpole(0,ic4),nterms,choose)
            endif
         enddo
C$OMP END PARALLEL DO
      enddo

      call bc2d_time(t11)
      timeinfo(2) = t11-t00

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     planewave translations

      
c     this goes over boxes so that none of the interactions
c     within a loop affects the same piece of data

      call bc2d_time(t00)

      do iiii = 0,2
      do jjjj = 0,2
C$OMP PARALLEL DO
C$OMP& PRIVATE(ii,j,i)
C$OMP& PRIVATE(jrow,jcol,xlength,dlogxlength)
C$OMP& IF(nboxes .gt. 1000)
C$OMP& SCHEDULE(dynamic,10)
         do ii = 1,nboxes
            j = iboxlev(ii)

            jrow = irowbox(j)
            jcol = icolbox(j)

            i = levelbox(j)

            if ( mod(jrow,3) .eq. iiii .and. mod(jcol,3) .eq. jjjj) then

               xlength = xlengths(i)
               dlogxlength = dlogxlengths(i+1)
               if(ichildbox(1,j) .gt. 0) then
c     Box has children: 

c     ship plane waves for same level and small-to-big far.
c     expnesw and expneswbig updated for every box in these
c     interactions. iflagnesw updated if a small-to-big far
c     interaction came into that box

                  call l2dpw_processall(j,icolleagbox,ichildbox,
     1                 icolbox,irowbox,localonoff,nterms,mpole,
     2                 nnodes,wnodes,comp,xsum,ratio,zs,dlogxlength,
     3                 expnesw,expneswbig,iflagnesw)
                  
               elseif (ichildbox(1,j) .le. 0) then

c     Box is childless:
c     (1) Do neighbor direct interactions with this box as target
c     (2) Send any big to small far interactions for children of
c        colleagues

                  call pbox2d_btosfarall(j,icolleagbox,ichildbox,
     1                 nnodes,npoly,poly2pw_btosfar,zs,xsum,xlength,
     2                 coeffs,expnesw)
                  
               endif
               
            endif
         enddo
C$OMP END PARALLEL DO
      enddo
      enddo

      call bc2d_time(t11)
      timeinfo(3) = t11-t00

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     downward pass (exp-to-loc and loc-loc, or exp-to-incoming
c                          and incoming-to-incoming)
c
      
      call bc2d_time(t00)
      do i = 0,nlev
c     
         istart = istartlev(i)
         iend = istart + nblevel(i) - 1
         xlength = xlengths(i)
         scal1 = xlength*2
         scal2 = xlength

C$OMP PARALLEL DO PRIVATE(jj,tmpexp,ii,ic1,ic2,ic3,ic4,j) 
C$OMP& IF(iend-istart .gt. 100) 
C$OMP& SCHEDULE(static,1)
         do jj = istart,iend
            j = iboxlev(jj)
            if(ichildbox(1,j).gt.0 .and. localonoff(j).eq.1)then
               ic1 = ichildbox(1,j)
               ic2 = ichildbox(2,j)
               ic3 = ichildbox(3,j)
               ic4 = ichildbox(4,j)

c     split parent expansion

               call l2dparchild(locexp(0,j),locexp(0,ic1),
     1              locexp(0,ic2),locexp(0,ic3),
     2              locexp(0,ic4),nterms, choose)

c     convert incoming planewaves to local (on children)
               
               do ii = 1, 4
                  ic1 = ichildbox(ii,j)
                  call l2dpw_exp4local(tmpexp,nterms,nnodes,
     1                 expnesw(0,1,ic1),expnesw(0,2,ic1),
     2                 expnesw(0,3,ic1),expnesw(0,4,ic1),comp)
                  call bc2d_zvecaddvec(ntermsp1,tmpexp,locexp(0,ic1))
               end do
            endif
         enddo
C$OMP END PARALLEL DO
      enddo

      call bc2d_time(t11)
      timeinfo(4) = t11-t00

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     LOCAL EVAL
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      

      call bc2d_time(t00)

      nnodesp1 = nnodes+1
      
C$OMP PARALLEL DO PRIVATE(ii,i,j,xlength,tmpvals,scal0,scal1,scal2,id)
C$OMP& IF(nboxes .gt. 1000)
C$OMP& SCHEDULE(static,1)
      do ii = 1,nboxes
         i = iboxlev(ii)
         if (ichildbox(1,i) .lt. 0 .and. localonoff(i) .eq. 1) then
            j = levelbox(i)
            xlength = xlengths(j)
            scal0 = 1
            scal1 = 2/xlength
            scal2 = scal1*scal1
            
            if (ifpot .eq. 1) then
               call bc2d_dscalzmatvec(npts,ntermsp1,scal0,
     1              mapslocpot,locexp(0,i),tmpvals)
               call bc2d_realaddvec(npts,tmpvals,pot(1,i))
               do id = 1,4
                  if (iflagnesw(id,i) .eq. 1) then
                     call bc2d_dscalzmatvec(npts,nnodesp1,scal0,
     1                    mapstobpot(1,0,id),expneswbig(0,id,i),tmpvals)
                     call bc2d_realaddvec(npts,tmpvals,pot(1,i))
                  endif
               enddo
            endif
            if (ifgrad .eq. 1) then
               call bc2d_dscalzmatvec(npts,ntermsp1,scal1,
     1              mapslocgrad,locexp(0,i),tmpvals)
               call pbox2d_convertzdertograd_add(npts,tmpvals,
     1              grad(1,1,i))
               do id = 1,4
                  if (iflagnesw(id,i) .eq. 1) then
                     call bc2d_dscalzmatvec(npts,nnodesp1,scal1,
     1                    mapstobgrad(1,0,id),expneswbig(0,id,i),
     2                    tmpvals)
                     call pbox2d_convertzdertograd_add(npts,tmpvals,
     1                    grad(1,1,i))
                  endif
               enddo
            endif
            if (ifhess .eq. 1) then
               call bc2d_dscalzmatvec(npts,ntermsp1,scal2,
     1              mapslochess,locexp(0,i),tmpvals)
               call pbox2d_convertzder2tohess_add(npts,tmpvals,
     1              hess(1,1,i))
            endif

            call pbox2d_direct(i,ichildbox,iparentbox,icolleagbox,npoly,
     1           coeffs,
     1           npts,tabcol,gradtabcol,hesstabcol,tabstob,gradtabstob,
     2           hesstabstob,tabbtos,gradtabbtos,hesstabbtos,xlength,
     1           ifpot,pot,ifgrad,grad,ifhess,hess)
            
         endif
      enddo
C$OMP END PARALLEL DO

      call bc2d_time(t11)
      timeinfo(5) = t11-t00

      return
      end
