c
c     The routines in this file are support routines for the
c     poisson box fmm. They are not intended to be user
c     callable. 
c      


      subroutine pbox2d_genpoly2mpmat(poly2mpmat,nterms,
     1     npoly,type,ndeg,nleg)
c
c     Generate the (nterms+1) x npoly matrix of integrals
c
c     int_{ [-1,1] x [-1,1] } ((x +iy)/2)^l p_k(x,y) dx dy 
c
c     Where p_k is the kth function in the requested
c     bivariate basis. The ordering convention is as in
c     legetens.f.
c      
c     This matrix maps coefficients in the polynomial basis
c     to a multipole expansion.
c
c     input:
c
c     nterms - order of mp expansion
c     npoly - total number of polynomials in basis
c     type - character. 'T' -> total degree polynomials,
c                       'F' -> full tensor product polynomials
c                       (see legetens.f for details)
c     ndeg - the degree of the polynomial basis
c     nleg - order of integral rule to use (recommend: norder+14)
c     tleg - real *8(nleg). precomputed 1D legendre nodes
c     wleg - real *8(nleg). precomputed 1D legendre weights
c      
      implicit real *8 (a-h,o-z)
      complex * 16 :: poly2mpmat(0:nterms,npoly)
      real *8, allocatable :: wleg(:), tleg(:)
      character :: type

      real *8, allocatable :: pols(:)
      complex *16, allocatable :: zmp(:)
      complex *16 :: im, zmp1, zmult
      real *8 :: pt(2), pi2
      data im / (0.0d0,1.0d0) /

      pi2 = 8*atan(1.0d0)
      
      allocate(pols(npoly),zmp(0:nterms),wleg(nleg),tleg(nleg))

      itype = 1
      call legeexps(itype,nleg,tleg,u,v,wleg)
      
      do k = 1,npoly
         do l = 0,nterms
            poly2mpmat(l,k) = 0
         enddo
      enddo
      
      do j = 1,nleg
         tj = tleg(j)
         wj = wleg(j)
         do i = 1,nleg
            ti = tleg(i)
            wi = wleg(i)
            wij = wi*wj
            zmp1 = 1
            zmult = (ti+im*tj)/2
            zmp(0) = -zmp1/pi2
            do k = 1,nterms
               zmp1 = zmp1*zmult
               zmp(k) = zmp1/k/pi2
            enddo
            pt(1) = ti
            pt(2) = tj
            call legetens_pols_2d(pt,ndeg,type,pols)
            do k = 1,npoly
               pk = pols(k)*wij
               do l = 0,nterms
                  poly2mpmat(l,k) = poly2mpmat(l,k)+pk*zmp(l)
               enddo
            enddo
         enddo
      enddo
                  
      return
      end

      subroutine pbox2d_stobfar_shifts(npts,pts,nnodes,
     1     xnodes,shifts)
      implicit real *8 (a-h,o-z)
      real *8 pts(2,npts), xnodes(0:nnodes)
      complex *16 shifts(npts,0:nnodes,4)
c     local
      complex *16 tempw, tempe, tempn, temps, imag
      data imag / (0.0d0,1.0d0) /

      do jj = 1,npts
         tempw = 2.0d0*(pts(1,jj) + imag*pts(2,jj))
         tempe = -tempshiftwest
         tempn = imag*tempshiftwest
         temps = -tempshiftnorth
         do ii = 1, nnodes
            shifts(jj,ii,2) = -cdexp(xnodes(ii) * tempe)
            shifts(jj,ii,4) = -cdexp(xnodes(ii) * tempw)
            shifts(jj,ii,1) = -cdexp(xnodes(ii) * tempn)
            shifts(jj,ii,3) = -cdexp(xnodes(ii) * temps)
         enddo
      enddo

      return
      end
      
      subroutine pbox2d_direct(ibox,ichildbox,iparentbox,icolleagbox,
     1     npoly,coeffs,npts,tabpot,tabgrad,tabhess,tabstobpot,
     2     tabstobgrad,tabstobhess,tabbtospot,tabbtosgrad,tabbtoshess,
     3     blength,ifpot,pot,ifgrad,grad,ifhess,hess)
cc
c     this routine does all incoming neighbor work for the box ibox
c     using tables
c
c     WARNING: the shift for the logarithmic potential is
c     computed assuming that the polynomial basis is the Legendre
c     basis AND that the first polynomial is the constant.
c
c     WARNING: b
c      
      implicit real *8 (a-h,o-z)
      integer ibox, ichildbox(4,*), icolleagbox(9,*), ifpgh,
     1     iparentbox(*)
      real *8 coeffs(npoly,*),tabpot(npts,npoly,9),
     1     tabgrad(2,npts,npoly,9),tabhess(3,npts,npoly,9),
     1     tabstobpot(npts,npoly,12), tabbtospot(npts,npoly,12),
     1     tabstobgrad(2,npts,npoly,12),tabstobhess(3,npts,npoly,12),
     1     tabbtosgrad(2,npts,npoly,12),tabbtoshess(3,npts,npoly,12)      
      real *8 pot(npts,*), grad(2,npts,*), hess(3,npts,*)
      real *8 blength
c     local
      integer nb, iout, ic1,ic2,ic3,ic4, iclose(2), nclose, nbtemp(2)
      integer nc, ic, nbt, npts2, ipar

      pi2 = 8*atan(1d0)

c
c     colleague work
c

      h = blength/2
      h2 = h*h
      h2dlogh = -4*log(h)*h2/pi2

      scalpot = h2
      scalgrad = h2/(blength/2)
      scalhess = 1
      
      do nb = 1,9
         iout = icolleagbox(nb,ibox)
         if (iout .lt. 0) cycle
         if (ichildbox(1,iout) .le. 0) then
c     childless colleague, do colleague local interaction
            nbt = 10-nb
            if (ifpot .eq. 1) then
               call bc2d_dscalmatvec_add(npts,npoly,scalpot,
     1              tabpot(1,1,nbt),coeffs(1,iout),pot(1,ibox))
               shiftpot = h2dlogh*coeffs(1,iout)
               call bc2d_dvecaddscal(npts,shiftpot,pot(1,ibox))
            endif
            if (ifgrad .eq. 1) then
               npts2 = npts*2
               call bc2d_dscalmatvec_add(npts2,npoly,scalgrad,
     1              tabgrad(1,1,1,nbt),
     1              coeffs(1,iout),grad(1,1,ibox))
            endif
            if (ifhess .eq. 1) then
               npts2 = npts*3
               call bc2d_dscalmatvec_add(npts2,npoly,scalhess,
     1              tabhess(1,1,1,nbt),coeffs(1,iout),hess(1,1,ibox))
            endif
         endif
      enddo

c
c     small to big work
c     

      h = blength/2/2
      h2 = h*h
      h2dlogh = -4*log(h)*h2/pi2

      scalpot = h2
      scalgrad = h
      scalhess = 1
      
      
      do nb = 1,9
         iout = icolleagbox(nb,ibox)
         if (iout .lt. 0) cycle
         if (ichildbox(1,iout) .gt. 0) then
            ic1 = ichildbox(1,iout)
            ic2 = ichildbox(2,iout)
            ic3 = ichildbox(3,iout)
            ic4 = ichildbox(4,iout)
            
            if(nb .eq. 1)then
c     lower left corner: one box not well sep, 3 are.
               nclose = 1
               iclose(1) = ic4
               nbtemp(1) = 12
            elseif(nb .eq. 2) then
C     immediately below, two boxes well separated
               nclose = 2
               iclose(1) = ic3
               iclose(2) = ic4
               nbtemp(1) = 11
               nbtemp(2) = 10
            elseif(nb .eq. 3)then
c     lower right corner, one box not well sep., 3 are.
               nclose = 1
               iclose(1) = ic3
               nbtemp(1) = 9
            elseif(nb .eq. 4)then
C     immediate left, two boxes well separated
               nclose = 2
               iclose(1) = ic2
               iclose(2) = ic4
               nbtemp(1) = 8
               nbtemp(2) = 6
            elseif(nb .eq. 5) then
               cycle
            elseif(nb .eq. 6)then
C     immediate right, two boxes well separated
               nclose = 2
               iclose(1) = ic1
               iclose(2) = ic3
               nbtemp(1) = 7
               nbtemp(2) = 5
            elseif(nb .eq. 7)then
c     upper left corner, one box not well sep., 3 are.
               nclose = 1
               iclose(1) = ic2
               nbtemp(1) = 4
            elseif(nb .eq. 8)then
C     immediately above, two boxes well separated
               nclose = 2
               iclose(1) = ic1
               iclose(2) = ic2
               nbtemp(1) = 3
               nbtemp(2) = 2
            elseif(nb .eq. 9)then
c     upper right corner, one box not well sep., 3 are
               nclose = 1
               iclose(1) = ic1
               nbtemp(1) = 1
            endif
            do nc = 1,nclose
               ic = iclose(nc)
               nbt = nbtemp(nc)
               if (ifpot .eq. 1) then
                  call bc2d_dscalmatvec_add(npts,npoly,scalpot,
     1                 tabstobpot(1,1,nbt),coeffs(1,ic),pot(1,ibox))
                  shiftpot = h2dlogh*coeffs(1,ic)                  
                  call bc2d_dvecaddscal(npts,shiftpot,pot(1,ibox))
               endif
               if (ifgrad .eq. 1) then
                  npts2 = npts*2
                  call bc2d_dscalmatvec_add(npts2,npoly,scalgrad,
     1                 tabstobgrad(1,1,1,nbt),
     1                 coeffs(1,ic),grad(1,1,ibox))
               endif
               if (ifhess .eq. 1) then
                  npts2 = npts*3
                  call bc2d_dscalmatvec_add(npts2,npoly,scalhess,
     1                 tabstobhess(1,1,1,nbt),coeffs(1,ic),
     2                 hess(1,1,ibox))
               endif
            enddo
         endif
      enddo

c
c     big to small work
c      

      h = blength
      h2 = h*h
      h2dlogh = -4*log(h)*h2/pi2

      scalpot = h2
      scalgrad = h
      scalhess = 1

      ipar = iparentbox(ibox)
      
      if (ipar .le. 0) return
      
      icself = -1
      do i = 1,4
         if (ichildbox(i,ipar) .eq. ibox) then
            icself = i
            exit
         endif
      enddo
      
c     loop over parent's colleagues
c     check for childless colleagues of parent
c     that neighbor this box
      
      do nb = 1,9
         iout = icolleagbox(nb,ipar)
         if (iout .lt. 0) cycle
         if (ichildbox(1,iout) .le. 0) then
            if(nb .eq. 1)then
c     lower left corner
               if (icself .eq. 1) then
                  nbt = 12
               else
                  cycle
               endif
            elseif(nb .eq. 2)then
c     immediately below
               if (icself .eq. 1) then
                  nbt = 10
               else if (icself .eq. 2) then
                  nbt = 11
               else
                  cycle
               endif
            elseif(nb .eq. 3) then
c     lower right corner
               if (icself .eq. 2) then
                  nbt = 9
               else
                  cycle
               endif
            elseif(nb .eq. 4)then
c     immediate left
               if (icself .eq. 1) then
                  nbt = 6
               else if (icself .eq. 3) then
                  nbt = 8
               else
                  cycle
               endif
            elseif(nb .eq. 5)then
c     parent, no interaction
               cycle
            elseif(nb .eq. 6)then
c     immediate right
               if (icself .eq. 2) then
                  nbt = 5
               else if (icself .eq. 4) then
                  nbt = 7
               else
                  cycle
               endif
            elseif(nb .eq. 7)then
c     top left
               if (icself .eq. 3) then
                  nbt = 4
               else
                  cycle
               endif
            elseif(nb .eq. 8)then
c     immediately above
               if (icself .eq. 3) then
                  nbt = 2
               else if (icself .eq. 4) then
                  nbt = 3
               else
                  cycle
               endif
            elseif(nb .eq. 9)then
c     top right
               if (icself .eq. 4) then
                  nbt = 1
               else
                  cycle
               endif
            endif
            if (ifpot .eq. 1) then
               call bc2d_dscalmatvec_add(npts,npoly,scalpot,
     1              tabbtospot(1,1,nbt),coeffs(1,iout),pot(1,ibox))
               shiftpot = h2dlogh*coeffs(1,iout)               
               call bc2d_dvecaddscal(npts,shiftpot,pot(1,ibox))
            endif
            if (ifgrad .eq. 1) then
               npts2 = npts*2
               call bc2d_dscalmatvec_add(npts2,npoly,scalgrad,
     1              tabbtosgrad(1,1,1,nbt),
     1              coeffs(1,iout),grad(1,1,ibox))
            endif
            if (ifhess .eq. 1) then
               npts2 = npts*3
               call bc2d_dscalmatvec_add(npts2,npoly,scalhess,
     1              tabbtoshess(1,1,1,nbt),coeffs(1,iout),
     2              hess(1,1,ibox))
            endif
         endif
      enddo
      
      return
      end
         
      subroutine pbox2d_genpoly2pwmat(norder,npoly,polytype,nnodes,
     1     wnodes,xnodes,poly2pw)
c
c     generate the matrix that directly forms the planewave expansions
c     on a box
c      
      implicit real *8 (a-h,o-z)
      real *8 :: xnodes(*), wnodes(*)
      complex *16 :: poly2pw(0:nnodes,4,npoly)
      character :: polytype
      real *8, allocatable :: t(:), w(:), polyvals(:,:), polyintsx(:,:)
      complex *16, allocatable :: polyintsy(:,:)
      complex *16 ztmp, im
      data im /(0d0,1d0)/

      pi = 4*atan(1d0)
      pi2 = 2*pi
      
      ndeg = norder - 1

      xmax = 0
      do j = 1,nnodes
         xmax = max(xmax,xnodes(j))
      enddo

c     estimate the number of legendre points to use based on largest
c     planewave node. (legendre nodes pi/2 less efficient than
c     equispaced). max freq is xmax/2
c     lam = 4pi/xmax. nwave = 2/lam = xmax/2 pi
c     4 pts per wave -> 2 xmax/pi points. convert to lege -> xmax points
      
      nleg = min(xmax+norder,300d0)
      allocate(t(nleg),w(nleg),polyvals(norder,nleg))
      itype = 1
      call legeexps(itype,nleg,t,u,v,w)

      do j = 1,nleg
         call legepols(t(j),ndeg,polyvals(1,j))
      enddo
      
c     use separability to evaluate integrals 
      
      allocate(polyintsx(norder,nnodes),polyintsy(norder,nnodes))

      do i = 1,nnodes
         do j = 1,norder
            polyintsx(j,i)=0
            polyintsy(j,i)=0
         enddo
      enddo

      do i = 1,nnodes
         do k = 1,nleg
            dtmp = exp(xnodes(i)*t(k)/2)*wnodes(i)*w(k)/xnodes(i)/pi2
            ztmp = exp(im*xnodes(i)*t(k)/2)*w(k)
            do j = 1,norder
               polyintsx(j,i)=polyintsx(j,i)+dtmp*polyvals(j,k)
               polyintsy(j,i)=polyintsy(j,i)+ztmp*polyvals(j,k)
            enddo
         enddo
      enddo

      
c     combine 1D integrals

      do i = 1,npoly
         do j = 0,nnodes
            poly2pw(j,1,i) = 0
            poly2pw(j,2,i) = 0
            poly2pw(j,3,i) = 0
            poly2pw(j,4,i) = 0
         enddo
      enddo
      
      if (polytype .eq. 'T') then
         do inode = 1,nnodes
            ipol = 0
            do i=1,norder
               do j=1,norder+1-i
                  ipol = ipol + 1
c     north 
                  poly2pw(inode,1,ipol) =
     1                 polyintsx(i,inode)*polyintsy(j,inode)*(-1)**(j-1)
c     east
                  poly2pw(inode,2,ipol) =
     1                 polyintsx(j,inode)*polyintsy(i,inode)
c     south
                  poly2pw(inode,3,ipol) =
     1                 polyintsx(i,inode)*polyintsy(j,inode)*(-1)**(i-1)
c     west
                  poly2pw(inode,4,ipol) =
     1                 polyintsx(j,inode)*polyintsy(i,inode)*(-1)**(i+j)
               enddo
            enddo
         enddo
      endif
               
      if (polytype .eq. 'F') then
         do inode = 1,nnodes
            ipol = 0
            do i=1,norder
               do j=1,norder
                  ipol = ipol + 1
c     north 
                  poly2pw(inode,1,ipol) =
     1                 polyintsx(i,inode)*polyintsy(j,inode)*(-1)**(j-1)
c     east
                  poly2pw(inode,2,ipol) =
     1                 polyintsx(j,inode)*polyintsy(i,inode)
c     south
                  poly2pw(inode,3,ipol) =
     1                 polyintsx(i,inode)*polyintsy(j,inode)*(-1)**(i-1)
c     west
                  poly2pw(inode,4,ipol) =
     1                 polyintsx(j,inode)*polyintsy(i,inode)*(-1)**(i+j)
               enddo
            enddo
         enddo
      endif
               
      return
      end
      
      subroutine pbox2d_genpoly2pwmat_btosfar(norder,npoly,polytype,
     1     nnodes,wnodes,xnodes,poly2pw)
c
c     generate the matrix that maps the density on a box to
c     planewave expansions on 1/2 the scale of the box (if it had children,
c     the size of those children). planewaves are recentered at
c     what would be the first child box.
c      
      implicit real *8 (a-h,o-z)
      real *8 :: xnodes(*), wnodes(*)
      complex *16 :: poly2pw(0:nnodes,4,npoly)
      character :: polytype
      real *8, allocatable :: t(:), w(:), polyvals(:,:), polyintsx(:,:)
      complex *16, allocatable :: polyintsy(:,:)
      complex *16 ztmp, im, zn, ze, zs, zw
      data im /(0d0,1d0)/
      
      pi = 4*atan(1d0)
      pi2 = 2*pi
      
      ndeg = norder - 1

      xmax = 0
      do j = 1,nnodes
         xmax = max(xmax,xnodes(j))
      enddo

c     estimate the number of legendre points to use based on largest
c     planewave node. (legendre nodes pi/2 less efficient than
c     equispaced). max freq is xmax
c     lam = 2pi/xmax. nwave = 2/lam = xmax/pi
c     4 pts per wave -> 4 xmax/pi points. convert to lege -> 2 xmax points
      
      nleg = min(2*xmax+norder,300d0)
      allocate(t(nleg),w(nleg),polyvals(norder,nleg))
      itype = 1
      call legeexps(itype,nleg,t,u,v,w)

      do j = 1,nleg
         call legepols(t(j),ndeg,polyvals(1,j))
      enddo
      
c     use separability to evaluate integrals 
      
      allocate(polyintsx(norder,nnodes),polyintsy(norder,nnodes))

      do i = 1,nnodes
         do j = 1,norder
            polyintsx(j,i)=0
            polyintsy(j,i)=0
         enddo
      enddo

      do i = 1,nnodes
         do k = 1,nleg
            dtmp = exp(xnodes(i)*t(k))*wnodes(i)*w(k)/xnodes(i)/pi2
            ztmp = exp(im*xnodes(i)*t(k))*w(k)
            do j = 1,norder
               polyintsx(j,i)=polyintsx(j,i)+dtmp*polyvals(j,k)
               polyintsy(j,i)=polyintsy(j,i)+ztmp*polyvals(j,k)
            enddo
         enddo
      enddo

      
c     combine 1D integrals

      do i = 1,npoly
         do j = 0,nnodes
            poly2pw(j,1,i) = 0
            poly2pw(j,2,i) = 0
            poly2pw(j,3,i) = 0
            poly2pw(j,4,i) = 0
         enddo
      enddo
      
      if (polytype .eq. 'T') then
         do inode = 1,nnodes
c     shift to first child
            zn = exp(xnodes(inode)*(0.5d0-im*0.5d0))
            ze = exp(xnodes(inode)*(0.5d0+im*0.5d0))
            zs = exp(xnodes(inode)*(-0.5d0+im*0.5d0))
            zw = exp(xnodes(inode)*(-0.5d0-im*0.5d0))
            ipol = 0
            do i=1,norder
               do j=1,norder+1-i
                  ipol = ipol + 1
c     north 
                  poly2pw(inode,1,ipol) = zn*
     1                 polyintsx(i,inode)*polyintsy(j,inode)*(-1)**(j-1)
c     east
                  poly2pw(inode,2,ipol) = ze*
     1                 polyintsx(j,inode)*polyintsy(i,inode)
c     south
                  poly2pw(inode,3,ipol) = zs*
     1                 polyintsx(i,inode)*polyintsy(j,inode)*(-1)**(i-1)
c     west
                  poly2pw(inode,4,ipol) = zw*
     1                 polyintsx(j,inode)*polyintsy(i,inode)*(-1)**(i+j)
               enddo
            enddo
         enddo
      endif
      
      if (polytype .eq. 'F') then
         do inode = 1,nnodes
c     shift to first child
            zn = exp(xnodes(inode)*(0.5d0-im*0.5d0))
            ze = exp(xnodes(inode)*(0.5d0+im*0.5d0))
            zs = exp(xnodes(inode)*(-0.5d0+im*0.5d0))
            zw = exp(xnodes(inode)*(-0.5d0-im*0.5d0))
            ipol = 0
            do i=1,norder
               do j=1,norder
                  ipol = ipol + 1
c     north 
                  poly2pw(inode,1,ipol) = zn*
     1                 polyintsx(i,inode)*polyintsy(j,inode)*(-1)**(j-1)
c     east
                  poly2pw(inode,2,ipol) = ze*
     1                 polyintsx(j,inode)*polyintsy(i,inode)
c     south
                  poly2pw(inode,3,ipol) = zs*
     1                 polyintsx(i,inode)*polyintsy(j,inode)*(-1)**(i-1)
c     west
                  poly2pw(inode,4,ipol) = zw*
     1                 polyintsx(j,inode)*polyintsy(i,inode)*(-1)**(i+j)
               enddo
            enddo
         enddo
      endif
               
      return
      end
      
      
      subroutine pbox2d_btosfarall(ibox,icolleagbox,ichildbox,
     1     nnodes,npoly,poly2pw_btosfar,zs,xsum,xlength,coeffs,
     2     expnesw)

      implicit real *8 (a-h,o-z)
      integer ibox, icolleagbox(9,*),ichildbox(4,*)
      real *8 :: coeffs(npoly,*)
      complex *16 :: expnesw(0:nnodes,4,*),
     1     poly2pw_btosfar(0:nnodes,4,npoly)
c     local
      logical ifmadebtos
      complex *16 zmul
      complex *16, allocatable :: extempnesw(:,:)
      complex *16 :: zs(0:nnodes,-3:3,-3:3,0:4)
      integer ispinord(2,4), ispinsign(2,4)
      data ispinord /   2,1, 1,2, 2,1, 1,2 /
      data ispinsign /  1,-1, 1,1, -1,1, -1,-1 /
      integer ifar(3), idir(3), nfar, nb, m, iout, itarg(2,3)
      
      ifmadebtos = .false.

      do nb = 1,9
         iout = icolleagbox(nb,ibox)
         if(iout .le. 0) cycle
         if (ichildbox(1,iout) .gt. 0) then
            ic1 = ichildbox(1,iout)
            ic2 = ichildbox(2,iout)
            ic3 = ichildbox(3,iout)
            ic4 = ichildbox(4,iout)
            
            if (.not. ifmadebtos) then
               m = 4*(nnodes+1)
               allocate(extempnesw(0:nnodes,4))
               scal = (xlength/2)**2
               call bc2d_dscalzmatdvec(m,npoly,scal,poly2pw_btosfar,
     1              coeffs(1,ibox),extempnesw)
               pi2 = 8*atan(1d0)
               dtmp = -4*scal*(xsum+log(xlength/2))*coeffs(1,ibox)/pi2
               extempnesw(0,1) = dtmp
               extempnesw(0,2) = dtmp
               extempnesw(0,3) = dtmp
               extempnesw(0,4) = dtmp
               ifmadebtos = .true.
            endif

            if(nb .eq. 1)then
c     lower left corner: one box not well sep, 3 are.
               nfar = 3
               ifar(1) = ic1
               ifar(2) = ic2
               ifar(3) = ic3
               itarg(1,1) = -2
               itarg(2,1) = -2
               itarg(1,2) = -1
               itarg(2,2) = -2
               itarg(1,3) = -2
               itarg(2,3) = -1
               idir(1) = 3
               idir(2) = 3
               idir(3) = 4
            elseif(nb .eq. 2) then
C     immediately below, two boxes well separated
               nfar = 2
               ifar(1) = ic1
               ifar(2) = ic2
               itarg(1,1) = 0
               itarg(2,1) = -2
               itarg(1,2) = 1
               itarg(2,2) = -2
               idir(1) = 3
               idir(2) = 3
            elseif(nb .eq. 3)then
c     lower right corner, one box not well sep., 3 are.
               nfar = 3
               ifar(1) = ic1  
               ifar(2) = ic2 
               ifar(3) = ic4
               itarg(1,1) = 2
               itarg(2,1) = -2
               itarg(1,2) = 3
               itarg(2,2) = -2
               itarg(1,3) = 3
               itarg(2,3) = -1
               idir(1) = 3
               idir(2) = 3
               idir(3) = 2
            elseif(nb .eq. 4)then
C     immediate left, two boxes well separated
               nfar = 2
               ifar(1) = ic1 
               ifar(2) = ic3
               itarg(1,1) = -2
               itarg(2,1) = 0
               itarg(1,2) = -2
               itarg(2,2) = 1
               idir(1) = 4
               idir(2) = 4
            elseif(nb .eq. 6)then
C     immediate right, two boxes well separated
               nfar = 2
               ifar(1) = ic4
               ifar(2) = ic2
               itarg(1,1) = 3
               itarg(2,1) = 1
               itarg(1,2) = 3
               itarg(2,2) = 0
               idir(1) = 2
               idir(2) = 2
            elseif(nb .eq. 7)then
c     upper left corner, one box not well sep., 3 are.
               nfar = 3
               ifar(1) = ic1
               ifar(2) = ic3
               ifar(3) = ic4
               itarg(1,1) = -2
               itarg(2,1) = 2
               itarg(1,2) = -2
               itarg(2,2) = 3
               itarg(1,3) = -1
               itarg(2,3) = 3
               idir(1) = 4
               idir(2) = 1
               idir(3) = 1
            elseif(nb .eq. 8)then
C     immediately above, two boxes well separated
               nfar = 2
               ifar(1) = ic3
               ifar(2) = ic4
               itarg(1,1) = 0
               itarg(2,1) = 3
               itarg(1,2) = 1
               itarg(2,2) = 3
               idir(1) = 1
               idir(2) = 1
            elseif(nb .eq. 9)then
c     upper right corner, one box not well sep., 3 are
               nfar = 3
               ifar(1) = ic3  
               ifar(2) = ic4
               ifar(3) = ic2
               itarg(1,1) = 2
               itarg(2,1) = 3
               itarg(1,2) = 3
               itarg(2,2) = 3
               itarg(1,3) = 3
               itarg(2,3) = 2
               idir(1) = 1
               idir(2) = 1
               idir(3) = 2
            endif
            
            do nf = 1,nfar
               ic = ifar(nf)
               id = idir(nf)

               io = ispinord(1,id)
               jo = ispinord(2,id)
               is = ispinsign(1,id)
               js = ispinsign(2,id)
               
               ii = is*itarg(io,nf)
               jj = js*itarg(jo,nf)
               
               do j = 0,nnodes
                  zmul = zs(j,ii,jj,0)
                  expnesw(j,id,ic) = expnesw(j,id,ic)
     1                 + extempnesw(j,id)*zmul
               enddo
               
            enddo
            
         endif
         
      enddo
            
      return
      end

      subroutine pbox2d_locevaltabs(npts,pts,nterms,
     1     mapslocpot,mapslocgrad,mapslochess,ier1)
      implicit real *8 (a-h,o-z)
      integer npts,nterms,ier1
      real *8 :: pts(2,npts)
      complex *16 :: mapslocpot(npts,0:nterms),
     1     mapslocgrad(npts,0:nterms), mapslochess(npts,0:nterms)
      
      complex *16 :: z, zj, zjm1, zjm2, im
      data im / (0d0,1d0) /
      
      do i = 1,npts
         z = (pts(1,i) + im*pts(2,i))/2
         zj = 1
         zjm1 = 0
         zjm2 = 0
         do j = 0,nterms
            mapslocpot(i,j) = zj
            mapslocgrad(i,j) = zjm1*j/2d0

            zjm2 = zjm1
            zjm1 = zj
            zj = zj*z
         enddo
      enddo

      return
      end

      subroutine pbox2d_stobfartabs(npts,pts,nnodes,xnodes,
     1     mapstobpot,mapstobgrad,mapstobhess,ier1)
      implicit real *8 (a-h,o-z)
      integer npts,nnodes,ier1
      real *8 :: pts(2,npts), xnodes(*)
      complex *16 :: mapstobpot(npts,0:nnodes,4),
     1     mapstobgrad(npts,0:nnodes,4), mapstobhess(npts,0:nnodes,4)
      
      complex *16 :: z, im, zn, ze, zs, zw
      data im / (0d0,1d0) /
      
      do i = 1,npts
         z = (pts(1,i) + im*pts(2,i))
         zn = im*z
         ze = -z
         zs = -im*z
         zw = z
         
         mapstobpot(i,0,1) = 1
         mapstobpot(i,0,2) = 1
         mapstobpot(i,0,3) = 1
         mapstobpot(i,0,4) = 1

         mapstobgrad(i,0,1) = 0
         mapstobgrad(i,0,2) = 0
         mapstobgrad(i,0,3) = 0
         mapstobgrad(i,0,4) = 0
         
         do j = 1,nnodes
            mapstobpot(i,j,1) = exp(xnodes(j)*zn)
            mapstobpot(i,j,2) = exp(xnodes(j)*ze)
            mapstobpot(i,j,3) = exp(xnodes(j)*zs)
            mapstobpot(i,j,4) = exp(xnodes(j)*zw)

            mapstobgrad(i,j,1) = exp(xnodes(j)*zn)*xnodes(j)*im
            mapstobgrad(i,j,2) = exp(xnodes(j)*ze)*xnodes(j)*(-1)
            mapstobgrad(i,j,3) = exp(xnodes(j)*zs)*xnodes(j)*(-im)
            mapstobgrad(i,j,4) = exp(xnodes(j)*zw)*xnodes(j)*(1)
         enddo
      enddo

      return
      end

      subroutine pbox2d_convertzdertograd_add(n,zder,grad)
      implicit real *8 (a-h,o-z)
      real *8 :: zder(2,n)
      real *8 :: grad(2,n)

      do j = 1,n
         grad(1,j) = grad(1,j) + zder(1,j)
         grad(2,j) = grad(2,j) - zder(2,j)
      enddo 
      return
      end

      subroutine pbox2d_convertzder2tohess_add(n,zder2,hess)
      implicit real *8 (a-h,o-z)
      real *8 :: zder2(2,n)
      real *8 :: hess(3,n)

      do j = 1,n
         hess(1,j) = hess(1,j) + zder2(1,j)
         hess(2,j) = hess(2,j) - zder2(2,j)
         hess(3,j) = hess(3,j) - zder2(1,j)
      enddo
      return
      end
