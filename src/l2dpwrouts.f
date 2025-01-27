C
C     plane wave routines for 2D Laplace fmm code 
C

      subroutine l2dpw_eval(ztarg,xnodes,nnodes,rscale,zbox,bexp,idir,
     3     ifpgh,pot,grad,hess)
      
      implicit real *8 (a-h,o-z)
      real *8 ztarg(2), zbox(2), xnodes(1:nnodes)
      complex *16 :: bexp(0:nnodes)
      complex *16 :: ztmp
      real *8 pot, grad(2), hess(3)
      complex *16 im, xfac, yfac
      data im /(0d0,1d0)/

      xtmp = ztarg(1)-zbox(1)
      ytmp = ztarg(2)-zbox(2)
      if (idir .eq. 1) then
         ztmp = im*(xtmp + im*ytmp)/rscale
         xfac = im/rscale
         yfac = -1/rscale
      endif
      if (idir .eq. 2) then
         ztmp = -(xtmp + im*ytmp)/rscale
         xfac = -1/rscale
         yfac = -im/rscale
      endif
      if (idir .eq. 3) then
         ztmp = -im*(xtmp + im*ytmp)/rscale
         xfac = -im/rscale
         yfac = 1/rscale
      endif
      if (idir .eq. 4) then
         ztmp = (xtmp + im*ytmp)/rscale
         xfac = 1/rscale
         yfac = im/rscale
      endif

      pot = bexp(0)
      do j = 1,nnodes
         pot = pot + exp(ztmp*xnodes(j))*bexp(j)
      enddo

      if (ifpgh .eq. 2 .or. ifpgh .eq. 3) then
         grad(1) = 0
         grad(2) = 0
         do j = 1,nnodes
            grad(1) = grad(1)
     1           + xfac*xnodes(j)*exp(ztmp*xnodes(j))*bexp(j)
            grad(2) = grad(2)
     1           + yfac*xnodes(j)*exp(ztmp*xnodes(j))*bexp(j)
         enddo
      endif
         
      
      return
      end

      
c***********************************************************************
c      subroutine precompute
c***********************************************************************
C     computes several arrays and values that
C     are needed for later use in the code.  
C
C     INPUT:
C 
C     NNODES  order of the plane wave expansions
C     NTERMS  order of the multipole expansions
C     XNODES  nodes in the plane wave expansions
C     WNODES  weights in the plane wave expansions
C
C     OUTPUT:
C 
C     COMP    array that is a combination of factorial
c                and nodal terms
C     RATIO   the ratio of the weights and nodes
C     CHOOSE  array of binomial coefficients 
C     XSUM    precomputed terms needed in exponential expansions
C
C***********************************************************************
      subroutine l2dpw_precomp(comp,nterms,nnodes,xnodes,wnodes,
     2         ratio,choose,xsum)
      implicit none
c-----GLOBAL VARIABLES
      integer  nnodes, nterms
      real *8  choose(2*nterms, 2*nterms)
      real *8  wnodes(*), xnodes(*), ratio(*)
      real *8  comp(0:nterms,nnodes), xsum, x
c-----lOCAL VARIABLES
      integer  i, j
      real *8, allocatable :: fact(:)
c
C     create factorial array:
C
      allocate(fact(0:nterms))
      fact(0) = 1.0d0
      do i = 1, nterms
         fact(i) = fact(i-1)*dble(i)
      end do
c
C     Compute CHOOSE(i,j) to be [(i-1) choose (j-1)].
C
      choose(1,1) = 1.0d0
      do i = 2, 2*nterms
         choose(1,i) = 0.0d0
         choose(i,1) = 1.0d0
      enddo
      do i = 2, 2*nterms
         do j = 2, 2*nterms
            choose(i,j) = choose(i-1,j-1) + choose(i-1,j)
         enddo
      enddo
      do i = 1, nnodes
         ratio(i) = wnodes(i) / xnodes(i)
      end do
C
C     Compute array COMP(i,j) = xnodes(i)^j / j!
C
      do i = 1, nnodes
         x = 1.0d0
         do j = 0, nterms
            comp(j,i) = x / fact(j)
            x = x * xnodes(i)
         enddo
      enddo
C
      xsum = 0.0d0
      do i = 1, nnodes
         xsum = xsum + dexp(-xnodes(i))*wnodes(i)/xnodes(i)
      enddo
      return
      end

      subroutine l2dpw_formbtos(xlength,src,ifcharge,charge,ifdipole,
     1     dipstr,ns,zbox,expnall,expeall,expsall,expwall,xsum,ratio,
     2     xnodes,wnodes,nnodes)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This subroutine forms planewaves for big to small far
c     interactions due to a collection of sources. 
c
c     The planewaves are centered at the 1st ghost child of the
c     given box and scaled for a box the size of the ghost children.
c
c     INPUT:
c
c     xlength - REAL *8, box length of current box (parent of ghosts)
c     src - REAL *8 ARRAY (2,*), array of source locations in box
c     ifcharge, ifdipole - INTEGER, flags for charges and dipoles
c     charge, dipstr - COMPLEX *16, strengths of charges and dipoles
c     ns - INTEGER, number of sources in box
c     zbox - REAL *8 (2), center of box
c     xsum - REAL *8, precomputed value (see PRECOMPUTE)
c     ratio - REAL *8 ARRAY, precomputed values (wnodes(i)/xnodes(i))
c     xnodes - REAL *8 ARRAY, planewave nodes
c     wnodes - REAL *8 ARRAY, planewave weights
c     nnodes - INTEGER, number of nodes for planewaves
c
c     OUTPUT:
c     
c     expnall, expsall, expwall, expeall - north, south, east, and west
c     plane wave expansions centered at the 4th ghost child. Used for
c     big to small far interactions
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      integer ifcharge, ifdipole, ns, nnodes
      real *8 xlength, src(2,*), zbox(2), xsum, ratio(*), xnodes(*)
      real *8 wnodes(*)
      complex *16 charge(*), dipstr(*), expnall(0:1), expsall(0:1)
      complex *16 expeall(0:1), expwall(0:1)
c     local
      integer i, j
      real *8 xlength2, zboxc1(2), dlogx2
      complex *16 expall(0:100), expall1(0:100), expall2(0:100)
      complex *16 eye, zero, one
      complex *16 spin, zshift
      data eye /(0.0d0,1.0d0)/
      data zero /(0.0d0,0.0d0)/
      data one /(1.0d0,0.0d0)/

      xlength2 = 0.5d0*xlength
      dlogx2 = dlog(xlength2)
      zboxc1(1) = zbox(1) - 0.25d0*xlength
      zboxc1(2) = zbox(2) - 0.25d0*xlength

      do j = 0,nnodes
         expnall(j) = zero
         expsall(j) = zero
         expeall(j) = zero
         expwall(j) = zero
      enddo

      do i = 1,ns
         do j = 0,nnodes
            expall1(j) = zero
            expall2(j) = zero
         enddo

c     set up planewave centered at source

         if (ifcharge .eq. 1) then
            expall1(0) = expall1(0) + charge(i)*(xsum+dlogx2)
            do j = 1,nnodes
               expall1(j) = expall1(j) + ratio(j)*charge(i)
            enddo
         endif
         if (ifdipole .eq. 1) then
            do j = 1,nnodes
               expall2(j) = expall2(j) - wnodes(j)*dipstr(i)/xlength2
            enddo
         endif

c     shift planewave

         zshift = zboxc1(1)-src(1,i) + (zboxc1(2)-src(2,i))*eye
         zshift = zshift/xlength2
         
         do j = 0,nnodes
            expall(j) = expall1(j) - eye*expall2(j)
         enddo
         spin = eye
         call l2dpw_expshift(expall,nnodes,xnodes,spin,zshift)
         call l2dpw_addexp(expall,expnall,nnodes)

         do j = 0,nnodes
            expall(j) = expall1(j) + eye*expall2(j)
         enddo
         spin = -eye
         call l2dpw_expshift(expall,nnodes,xnodes,spin,zshift)
         call l2dpw_addexp(expall,expsall,nnodes)

         do j = 0,nnodes
            expall(j) = expall1(j) + expall2(j)
         enddo
         spin = -one
         call l2dpw_expshift(expall,nnodes,xnodes,spin,zshift)
         call l2dpw_addexp(expall,expeall,nnodes)

         do j = 0,nnodes
            expall(j) = expall1(j) - expall2(j)
         enddo
         spin = one
         call l2dpw_expshift(expall,nnodes,xnodes,spin,zshift)
         call l2dpw_addexp(expall,expwall,nnodes)
      enddo

      return
      end
C
C
C
c***********************************************************************
c     subroutine exp4local
c***********************************************************************
C     takes four plane wave expansions (from four
C     different directions) and converts them to a single
C     local expansion.  
C
C     INPUT:
C
C     NTERMS   order of the multipole expansions
C     NNODES   order of the plane wave expansions
C     BETAW    incoming west exponential coefficients
C     BETAS    incoming south exponential coefficients
C     BETAE    incoming east exponential coefficients
C     BETAN    incoming north exponential coefficients
C     COMP     precomputed term involving the ratio of
C                the weights and factorial terms
C
C     OUTPUT:
C
C     BETAHAT  the local coefficients of the parent box
C
C***********************************************************************
      subroutine l2dpw_exp4local(betahat,nterms,nnodes,betan,
     1                     betae,betas,betaw,comp)
      implicit none
c-----Global variables
      integer nnodes, nterms
      real *8 comp(0:nterms,nnodes)
      complex *16 betahat(0:nterms)
      complex *16 betan(0:nnodes),betas(0:nnodes)
      complex *16 betae(0:nnodes),betaw(0:nnodes)
c-----Local variables
      integer i, j
      complex *16 bsum1,bsum2,bsum3,bsum4
      complex *16 bs14,bs3m2,bs1m4,bs32
      complex *16 imag
      data imag/(0.0d0, 1.0d0)/
C
C     local coefficients are computed below using Taylor expansion
C
C     Upon input, the exponential expansions are in the 
C     form sum_(i=0,nnodes)beta(i)*exp(spin*xnodes(i)*z).
C     Upon output, the local expansion is of the form
C     sum(i 0 to p) betahat(i)*z^i.
C
      betahat(0) = 0.0d0
      do i = 1, nnodes
         betahat(0) = betahat(0) + (betaw(i) + betae(i) +
     1                betas(i) + betan(i))
      enddo
      do j = 1, nterms
         betahat(j) = 0.0d0
      enddo
      do i = 1, nnodes
         bsum1 = imag*(betan(i) - betas(i))
         bsum2 = betae(i) + betaw(i)
         bsum3 = betan(i) + betas(i)
         bsum4 = betaw(i) - betae(i)
         bs14  = bsum4 + bsum1
         bs3m2 = bsum3 - bsum2
         bs1m4 = bsum1 - bsum4
         bs32  = bsum3 + bsum2
         do j = 1, nterms-3, 4
           betahat(j)   = betahat(j)   + comp(j,i)*bs14
           betahat(j+1) = betahat(j+1) - comp(j+1,i)*bs3m2
           betahat(j+2) = betahat(j+2) - comp(j+2,i)*bs1m4
           betahat(j+3) = betahat(j+3) + comp(j+3,i)*bs32
         end do
      end do
      betahat(0) = betahat(0) + betaw(0) + betae(0)
     1           + betas(0) + betan(0)
      return
      end
C
C
C
C
c***********************************************************************
C     subroutine expevalall
c***********************************************************************
C     used to evaluate the 'Big'
C     expansions.  These are the expansions that result from the
C     Stobfar interactions in plane wave form.
C
C     INPUT:
C
C     NNODES    plane wave expansion length
C     EXPW, EXPS, EXPE, and EXPN are the exponential expansions that
C     represent the incoming Stobfar interactions.  
C     IFLAGEAST, IFLAGWEST, IFLAGNORTH, and IFLAGSOUTH are flags that
C     indicate whether or not a given expansion needs to be evaluated.
C     TEMPEAST,TEMPNORTH,TEMPWEST,TEMPSOUTH are precomputed arrays 
C     that may be needed on the evaluation.
C
C     IFGRAD flag = 1, then compute gradient of potential
C     IFHESS flag = 1, then compute Hessian of potential
C
C     OUTPUT:
C
C     POT is the *incremented* potential
C     GRAD is the *incremented* gradient
C     HESS is the *incremented* Hessian
C***********************************************************************
      subroutine l2dpw_expevalall(nt,ztarg,
     1     ifpot,pot,ifgrad,grad,ifhess,hess,
     2     xnodes,nnodes,rscale,zbox,expn,expe,exps,expw,
     3     iflagnorth,iflageast,iflagsouth,iflagwest)
      implicit none
c-----Global variables
      integer  nnodes, nt
      integer  iflageast, iflagnorth
      integer  iflagwest, iflagsouth
      integer  ifgrad, ifhess, ifpot
      complex *16  pot(*), grad(2,*), hess(3,*)
      real *8  xnodes(*), ztarg(2,*), zbox(2), rscale
      complex *16  expn(0:nnodes), exps(0:nnodes)
      complex *16  expe(0:nnodes), expw(0:nnodes)
c-----Local variables
      integer  i, j
      real *8 tempr, tempi, h2, xtemp, xtemps, xtarg, ytarg
      complex *16 tempshiftwest, tempshifteast, tempshiftnorth
      complex *16 tempshiftsouth, tempwest, tempeast, tempnorth
      complex *16 tempsouth, eye, ctemp
      data eye / (0.0d0,1.0d0) /
      h2 = rscale

      do j = 1,nt
         xtarg = ztarg(1,j)-zbox(1)
         ytarg = ztarg(2,j)-zbox(2)
         if(iflagwest.eq.1) then
            tempshiftwest = (xtarg+eye*ytarg)/rscale
            pot(j) = pot(j) + expw(0)
            do i = 1, nnodes
               tempwest = -cdexp(xnodes(i)*tempshiftwest)
               ctemp = expw(i)*tempwest
               tempr = dreal(ctemp)
               tempi = dimag(ctemp)
               xtemp = xnodes(i)/h2
               xtemps = xtemp*xtemp
               pot(j) = pot(j) + ctemp
               grad(1,j) = grad(1,j) + tempr*xtemp
               grad(2,j) = grad(2,j) - tempi*xtemp
               hess(1,j) = hess(1,j) + tempr*xtemps
               hess(2,j) = hess(2,j) - tempi*xtemps
               hess(3,j) = hess(3,j) - tempr*xtemps
            enddo
         endif
         if(iflageast.eq.1) then
            tempshifteast = -(xtarg+eye*ytarg)/rscale
            pot(j) = pot(j) + expe(0)
            do i = 1, nnodes
               tempeast = -cdexp(xnodes(i)*tempshifteast)
               ctemp = expe(i)*tempeast
               tempr = dreal(ctemp)
               tempi = dimag(ctemp)
               xtemp = xnodes(i)/h2
               xtemps = xtemp*xtemp
               pot(j) = pot(j) + ctemp
               grad(1,j) = grad(1,j) - tempr*xtemp
               grad(2,j) = grad(2,j) + tempi*xtemp
               hess(1,j) = hess(1,j) + tempr*xtemps
               hess(2,j) = hess(2,j) - tempi*xtemps
               hess(3,j) = hess(3,j) - tempr*xtemps
            enddo
         endif
         if(iflagnorth.eq.1) then
            tempshiftnorth = eye*(xtarg+eye*ytarg)/rscale
            pot(j) = pot(j) + expn(0)
            do i = 1, nnodes
               tempnorth = -cdexp(xnodes(i)*tempshiftnorth)
               ctemp = expn(i)*tempnorth
               tempr = dreal(ctemp)
               tempi = dimag(ctemp)
               xtemp = xnodes(i)/h2
               xtemps = xtemp*xtemp
               pot(j) = pot(j) + ctemp
               grad(1,j) = grad(1,j) - tempi*xtemp
               grad(2,j) = grad(2,j) - tempr*xtemp
               hess(1,j) = hess(1,j) - tempr*xtemps
               hess(2,j) = hess(2,j) + tempi*xtemps
               hess(3,j) = hess(3,j) + tempr*xtemps
            enddo
         endif
         if(iflagsouth .eq. 1)then
            tempshiftsouth = -eye*(xtarg+eye*ytarg)/rscale
            pot(j) = pot(j) + exps(0)
            do i = 1, nnodes
               tempsouth = -cdexp(xnodes(i)*tempshiftsouth)
               ctemp = exps(i)*tempsouth
               tempr = dreal(ctemp)
               tempi = dimag(ctemp)
               xtemp = xnodes(i)/h2
               xtemps = xtemp*xtemp
               pot(j) = pot(j) + ctemp
               grad(1,j) = grad(1,j) + tempi*xtemp
               grad(2,j) = grad(2,j) + tempr*xtemp
               hess(1,j) = hess(1,j) - tempr*xtemps
               hess(2,j) = hess(2,j) + tempi*xtemps
               hess(3,j) = hess(3,j) + tempr*xtemps
            enddo
         endif
      enddo

      return
      end
C
C
C
c***********************************************************************
C     subroutine expbigshiftnesw
c***********************************************************************
C
C     INPUT:
C
C     NNODES    plane wave expansion length
C     EXPW, EXPS, EXPE, and EXPN are the exponential expansions that
C     represent the incoming Stobfar interactions.  
C     IFLAGEAST, IFLAGWEST, IFLAGNORTH, and IFLAGSOUTH are flags that
C     indicate whether or not a given expansion needs to be evaluated.
C
C     IFGRAD flag = 1, then compute gradient of potential
C     IFHESS flag = 1, then compute Hessian of potential
C
C     OUTPUT:
C     
C     betaw - shifted west expansion
C     betas - shifted south expansion
C     betae - shifted east expansion
C     betan - shifted north expansion
c
C***********************************************************************
      subroutine l2dpw_expshiftnesw(ztarg,xnodes,nnodes,rscale,zbox,
     1     expn,expe,exps,expw,betan,betae,betas,betaw,
     2     iflagnorth,iflageast,iflagsouth,iflagwest)
      implicit none
c-----Global variables
      integer  nnodes
      integer  iflageast, iflagnorth
      integer  iflagwest, iflagsouth
      complex *16 betan(0:nnodes), betae(0:nnodes)
      complex *16 betas(0:nnodes), betaw(0:nnodes)
      real *8  xnodes(*), ztarg(2), zbox(2), rscale
      complex *16  expn(0:nnodes), exps(0:nnodes)
      complex *16  expe(0:nnodes), expw(0:nnodes)
c-----Local variables
      integer  i
      real *8 xtarg, ytarg
      complex *16 tempshift, temp, eye
      data eye / (0.0d0,1.0d0) /

      xtarg = (ztarg(1)-zbox(1))
      ytarg = (ztarg(2)-zbox(2))
      if(iflagwest.eq.1) then
         tempshift = (xtarg+eye*ytarg)/rscale
         betaw(0) = expw(0)
         do i = 1, nnodes
            temp = cdexp(xnodes(i)*tempshift)
            betaw(i) = expw(i)*temp
         enddo
      endif
      if(iflageast.eq.1) then
         tempshift = -(xtarg+eye*ytarg)/rscale
         betae(0) = expe(0)
         do i = 1, nnodes
            temp = cdexp(xnodes(i)*tempshift)
            betae(i) = expe(i)*temp
         enddo
      endif
      if(iflagnorth.eq.1) then
         tempshift = eye*(xtarg+eye*ytarg)/rscale
         betan(0) = expn(0)
         do i = 1, nnodes
            temp = cdexp(xnodes(i)*tempshift)
            betan(i) = expn(i)*temp
         enddo
      endif
      if(iflagsouth.eq.1) then
         tempshift = -eye*(xtarg+eye*ytarg)/rscale
         betas(0) = exps(0)
         do i = 1, nnodes
            temp = cdexp(xnodes(i)*tempshift)
            betas(i) = exps(i)*temp
         enddo
      endif

      return
      end
C     
C
C
C
c***********************************************************************
C     subroutine mkshifts2d
C***********************************************************************
C     This subroutine computes the tables of exponentials needed
C     for translating exponential representations of harmonic
C     functions, discretized via Norman's quadratures.
C
C     INPUT:
C
C     XNODES  nodes in the plane wave expansions
C     NNODES  order of the plane wave expansions
C
C     ON OUTPUT:
C
C     ZS(K,N,M)   e^{- \lambda_k (n + im)}.
C     ZSEAST(K,N,M) -> extra scaling to 
C     ZSWEST(K,N,M)
C     ZSNORTH(K,N,M)
C     ZSSOUTH(K,N,M)
C***********************************************************************
      subroutine l2dpw_mkshifts(xnodes,nnodes,zs)
      implicit none
      complex *16  zs(0:nnodes,-3:3,-3:3,0:4)
      complex *16 tempexp
      real *8     xnodes(nnodes)
      integer  nnodes, k, m, n
c
      do m = -3,3
         do n = -3,3
            do k = 1,nnodes
               zs(k,n,m,0)=cdexp(-xnodes(k)*dcmplx(dble(n),dble(m)))
            enddo
         enddo
      enddo
      do m = -3,3
         do n = -3,3
            zs(0,n,m,0)=1.0d0
         enddo
      enddo
      do k = 1,nnodes
         tempexp = cdexp(-xnodes(k)*dcmplx(0.5d0,-0.5d0))
         do n = -3,3
            do m = -3,3
               zs(k,n,m,1)=zs(k,n,m,0) * tempexp
            enddo
         enddo
      enddo
      do n = -3,3
         do m = -3,3
            zs(0,n,m,1)=1.0d0
         enddo
      enddo
      do k = 1,nnodes
         tempexp = cdexp(-xnodes(k)*dcmplx(0.5d0,0.5d0))
         do n = -3,3
            do m = -3,3
               zs(k,n,m,2)=zs(k,n,m,0) * tempexp
            enddo
         enddo
      enddo
      do n = -3,3
         do m = -3,3
            zs(0,n,m,2)=1.0d0
         enddo
      enddo
      do k = 1,nnodes
         tempexp = cdexp(-xnodes(k)*dcmplx(-0.5d0,0.5d0))
         do n = -3,3
            do m = -3,3
               zs(k,n,m,3)=zs(k,n,m,0) * tempexp
            enddo
         enddo
      enddo
      do n = -3,3
         do m = -3,3
            zs(0,n,m,3)=1.0d0
         enddo
      enddo
      do k = 1,nnodes
         tempexp = cdexp(-xnodes(k)*dcmplx(-0.5d0,-0.5d0))
         do n = -3,3
            do m = -3,3
               zs(k,n,m,4)=zs(k,n,m,0) * tempexp
            enddo
         enddo
      enddo
      do n = -3,3
         do m = -3,3
            zs(0,n,m,4)=1.0d0
         enddo
      enddo

      return
      end
c
c
c***********************************************************************
      subroutine l2dpw_mkexp(ibox,nterms,mpole,nnodes,mexnall,mexn34,
     2           mexeall,mexe14,mexee4,mexee1,mexsall,mexs12,
     3           mexwall,mexw23,mexww3,mexww2,zs,comp,wnodes,
     4           scale,xsum,ratio,ichildbox)
c***********************************************************************
C     This subroutine creates the 
c     north (+y), south (-y), east and west plane wave 
C     expansions for a parent box due to all four children. 
C
C     children are ordered
C
C                      | 1  | 2  |
C                      |----|----|
C                      | 4  | 3  |
C
C     Some intelligence is used in the order of summation. 
C
C     INPUT:
C
C     ibox       box number of parent box under consideration.
c     nterms     order of the multipole expansions
c     mpole(0:nterms,nboxes)  array of 
c                             multipole expansions for all boxes
c     nnodes     order of the plane wave expansions
c     zs(0:nnodes,-3:3,-3:3) precomputed array
c              of shifts for the plane wave expansions
c     comp(0:nterms,nnodes) (real *8) precomputed array containing
c                             wnodes(i) xnodes(i)^j / j!
c     wnodes     quadrature weights in plane wave formula
c     scale      box scaling parameter
c     xsum       sum of exp(-xnodes(i)) * wnodes(i) / xnodes(i)
c     ratio       ratio(i) = wnodes(i) / xnodes(i)
c     ichildbox(4,nboxes)     array containing the
c                             children of each box
C     OUTPUT:
C
C     mexnall(0:nnodes) contribution to north
c           expansion for all boxes
c     mexn34(0:nnodes ) contribution to north
c           expansion for boxes 3 and 4
c     etc.
c***********************************************************************
      implicit none
      integer  ibox
      integer  nnodes,ic
      integer  ichild(4)
      integer  ichildbox(4,1)
      integer  nterms, jj
      real *8  xlength
      real *8  comp(0:nterms,nnodes)
      real *8  wnodes(nnodes),scale, xsum, sum1
      real *8  ratio(nnodes)
      complex *16  mpole(0:nterms,1)
      complex *16  zs(0:nnodes,-3:3,-3:3)
      complex *16  expn(0:100),exps(0:100),expe(0:100),expw(0:100)
      complex *16  mexnall(0:nnodes),mexn34(0:nnodes)
      complex *16  mexsall(0:nnodes),mexs12(0:nnodes)
      complex *16  mexeall(0:nnodes),mexe14(0:nnodes)
      complex *16  mexee4(0:nnodes),mexee1(0:nnodes)
      complex *16  mexwall(0:nnodes),mexw23(0:nnodes)
      complex *16  mexww3(0:nnodes),mexww2(0:nnodes)
C***********************************************************************
      xlength = 1.0d0/scale
      sum1 = xsum + dlog(xlength)
      ichild(1) = ichildbox(1,ibox)
      ichild(2) = ichildbox(2,ibox)
      ichild(3) = ichildbox(3,ibox)
      ichild(4) = ichildbox(4,ibox)
c
C     include contributions from child 4
C
      ic = ichild(4)
      call l2dpw_expcoeff(expn,expe,exps,expw,nnodes,mpole(0,ic),
     1              nterms,comp,wnodes,sum1,ratio)
      do jj = 0,nnodes
         mexn34(jj) = expn(jj)
         mexsall(jj) = exps(jj)
      enddo
      do jj = 0,nnodes
         mexee4(jj) = expe(jj)
         mexwall(jj) = expw(jj)
      enddo
c
C     include contributions from child 3
C
      ic = ichild(3)
      call l2dpw_expcoeff(expn,expe,exps,expw,nnodes,mpole(0,ic),
     1              nterms,comp,wnodes,sum1,ratio)
      do jj = 0,nnodes
         expn(jj) = expn(jj)*zs(jj,0,1)
         mexn34(jj) = mexn34(jj) + expn(jj)
         exps(jj) = exps(jj)*zs(jj,0,-1)
         mexsall(jj) = mexsall(jj) + exps(jj)
      enddo

      do jj = 0,nnodes
         mexeall(jj) = expe(jj)*zs(jj,-1,0)
         mexww3(jj) = expw(jj)*zs(jj,1,0)
      enddo
c
C     include contributions from child 1
C
      ic = ichild(1)
      call l2dpw_expcoeff(expn,expe,exps,expw,nnodes,mpole(0,ic),
     1              nterms,comp,wnodes,sum1,ratio)
      do jj = 0,nnodes
         expn(jj) = expn(jj)*zs(jj,-1,0)
         mexnall(jj) = expn(jj) + mexn34(jj)
         exps(jj) = exps(jj)*zs(jj,1,0)
         mexs12(jj) = exps(jj)
      enddo
      do jj = 0,nnodes
         mexee1(jj) = expe(jj)*zs(jj,0,-1)
         mexe14(jj) = mexee4(jj) + mexee1(jj)
         mexeall(jj) = mexeall(jj) + mexe14(jj)
         expw(jj) = expw(jj)*zs(jj,0,1)
         mexwall(jj) = mexwall(jj) + expw(jj)
      enddo
c
C     include contributions from child 2
C
      ic = ichild(2)
      call l2dpw_expcoeff(expn,expe,exps,expw,nnodes,mpole(0,ic),
     1              nterms,comp,wnodes,sum1,ratio)
      do jj = 0,nnodes
         expn(jj) = expn(jj)*zs(jj,-1,1)
         mexnall(jj) = mexnall(jj) + expn(jj)
         exps(jj) = exps(jj)*zs(jj,1,-1)
         mexs12(jj) = mexs12(jj) + exps(jj)
         mexsall(jj) = mexsall(jj) + mexs12(jj)
      enddo
      do jj = 0,nnodes
         expe(jj) = expe(jj)*zs(jj,-1,-1)
         mexeall(jj) = mexeall(jj) + expe(jj)
         mexww2(jj) = expw(jj)*zs(jj,1,1)
         mexw23(jj) = mexww3(jj) + mexww2(jj)
         mexwall(jj) = mexwall(jj) + mexw23(jj)
      enddo
      return
      end
C

      subroutine l2dpw_processall(ibox,icolleagbox,ichildbox,
     1     icolbox,irowbox,localonoff,nterms,mpole,
     2     nnodes,wnodes,comp,xsum,ratio,zs,dlogxlength,
     3     expnesw,expneswbig,iflagnesw)

      implicit none
      integer ibox,icolleagbox(9,*),ichildbox(4,*),icolbox(*)
      integer irowbox(*),localonoff(*),nterms,nnodes
      complex *16 mpole(0:nterms,*),zs(0:nnodes,-3:3,-3:3,0:4)
      real *8 wnodes(*),comp(0:nterms,nnodes),xsum,ratio(*)
      real *8 scale
      complex *16 expnesw(0:nnodes,4,*), expneswbig(0:nnodes,4,*)
      integer iflagnesw(4,*)
c     local variables

c     list storage
      integer  inall(6),nnall,iynall(6)
      integer  isall(6),nsall,iysall(6)
      integer  ieall(4),neall,iyeall(4)
      integer  iwall(4),nwall,iywall(4)
      integer  in12(4),nn12,iy12(4)
      integer  is34(4),ns34,iy34(4)
      integer  iw24(2),nw24,iy24(2)
      integer  ie13(2),ne13,iy13(2)
      integer  iee1(1),nee1,iey1(1)
      integer  iee3(1),nee3,iey3(1)
      integer  iww4(1),nww4,iwy4(1)
      integer  iww2(1),nww2,iwy2(1)
      integer  inbig12(3), isbig34(3)
      integer  iebig13(1), iwbig24(1)
      integer  iebig3(1), iwbig2(1)
      integer  iebig1(1), iwbig4(1)

c     temporary expansions
      integer lexp
      parameter (lexp=60)
      complex *16 extempn(0:lexp),extempe(0:lexp),extemps(0:lexp)
      complex *16 extempw(0:lexp)
      complex *16 exnall(0:lexp), exn12(0:lexp), exsall(0:lexp)
      complex *16 exs34(0:lexp), exeall(0:lexp), exwall(0:lexp)
      complex *16 exe13(0:lexp), exw24(0:lexp)
      complex *16 :: zmul
      
c     indices
      integer iy, i, j, ii, ic1, ic2, ic3, ic4, nnodesp1

      real *8 dlogxlength, sum1

      nnodesp1 = nnodes+1
      
c     get list info

      call lrt2d_mklists(ibox,inall,nnall,iynall,in12,nn12,iy12,
     1    isall,nsall,iysall,is34,ns34,iy34,ieall,neall,
     2    iyeall,ie13,ne13,iy13,iwall,nwall,iywall,
     3    iw24,nw24,iy24,iww2,iwy2,nww2,iww4,iwy4,nww4,
     4    iee1,iey1,nee1,iee3,iey3,nee3,inbig12,isbig34,iebig13,
     6    iwbig24,iebig1,iwbig2,iebig3,iwbig4,icolleagbox,ichildbox,
     8    icolbox,irowbox,iflagnesw)

c     make temporary expansions and agglomerate when appropriate
c     some optimization is done in the ordering of the agglomeration
c     (these are original to the Ethridge code)
c     warning: the temporary arrays are intialized by the first
c     write.
c     immediately send single child interactions (ee1,ee3,etc)

      sum1 = xsum + dlogxlength

      ic1 = ichildbox(1,ibox)
      ic2 = ichildbox(2,ibox)
      ic3 = ichildbox(3,ibox)
      ic4 = ichildbox(4,ibox)      

c     child 1

      call l2dpw_expcoeff(extempn,extempe,extemps,extempw,nnodes,
     1     mpole(0,ic1),nterms,comp,wnodes,sum1,ratio)

      do j = 0,nnodes
         exn12(j) = extempn(j)
         exe13(j) = extempe(j)
         exsall(j) = extemps(j)
         exwall(j) = extempw(j)
      enddo
      
      do i = 1,nee1
         if(localonoff(iee1(i)) .eq. 1)then
            ii = iee1(i)
            iy = iey1(i)
            call bc2d_zdiagmatvec_add(nnodesp1,zs(0,2,iy,0),extempe,
     1           expnesw(0,2,ii))
         endif
      enddo

      ii = iebig1(1)
      if(ii .gt. 0 .and. localonoff(ii) .eq. 1) then
         call bc2d_zdiagmatvec_add(nnodesp1,zs(0,2,-2,2),extempe,
     1        expneswbig(0,2,ii))
      endif

c     child 2

      call l2dpw_expcoeff(extempn,extempe,extemps,extempw,nnodes,
     1     mpole(0,ic2),nterms,comp,wnodes,sum1,ratio)

      do j = 0,nnodes
         extempn(j) = extempn(j)*zs(j,0,1,0)
         exn12(j) = exn12(j) + extempn(j)

         extempe(j) = extempe(j)*zs(j,-1,0,0)
         exeall(j) = extempe(j)
         
         extemps(j) = extemps(j)*zs(j,0,-1,0)
         exsall(j) = exsall(j) + extemps(j)
         
         extempw(j) = extempw(j)*zs(j,1,0,0)
         exw24(j) = extempw(j)
      enddo
      
      do i = 1,nww2
         ii = iww2(i)
         if(localonoff(ii) .eq. 1)then
            iy = -iwy2(i)
            call bc2d_zdiagmatvec_add(nnodesp1,zs(0,1,iy,0),extempw,
     1           expnesw(0,4,ii))
         endif
      enddo

      ii = iwbig2(1)
      if(ii .gt. 0 .and. localonoff(ii) .eq. 1) then
         call bc2d_zdiagmatvec_add(nnodesp1,zs(0,2,2,4),extempw,
     1        expneswbig(0,4,ii))
      endif

c     child 3
c
      call l2dpw_expcoeff(extempn,extempe,extemps,extempw,nnodes,
     1     mpole(0,ic3),nterms,comp,wnodes,sum1,ratio)

      do j = 0,nnodes
         extempn(j) = extempn(j)*zs(j,-1,0,0)
         exnall(j) = exn12(j) + extempn(j)

         extempe(j) = extempe(j)*zs(j,0,-1,0)
         exe13(j) = exe13(j) + extempe(j)
         exeall(j) = exeall(j) + exe13(j)
         
         extemps(j) = extemps(j)*zs(j,1,0,0)
         exs34(j) = extemps(j)

         extempw(j) = extempw(j)*zs(j,0,1,0)
         exwall(j) = exwall(j) + extempw(j)
      enddo
      
      ii = iebig3(1)
      if(ii .gt. 0 .and. localonoff(ii) .eq. 1) then
         call bc2d_zdiagmatvec_add(nnodesp1,zs(0,2,2,2),extempe,
     1        expneswbig(0,2,ii))
      endif         

      do i = 1,nee3
         ii = iee3(i)
         if(localonoff(ii) .eq. 1)then
            iy = iey3(i)
            call bc2d_zdiagmatvec_add(nnodesp1,zs(0,2,iy,0),extempe,
     1           expnesw(0,2,ii))
         endif
      enddo
c
c     child 4

      call l2dpw_expcoeff(extempn,extempe,extemps,extempw,nnodes,
     1     mpole(0,ic4),nterms,comp,wnodes,sum1,ratio)

      do j = 0,nnodes
         extempn(j) = extempn(j)*zs(j,-1,1,0)
         exnall(j) = exnall(j) + extempn(j) 

         extempe(j) = extempe(j)*zs(j,-1,-1,0)
         exeall(j) = exeall(j) + extempe(j)
         
         extemps(j) = extemps(j)*zs(j,1,-1,0)
         exs34(j) = exs34(j) + extemps(j)
         exsall(j) = exsall(j) + exs34(j)

         extempw(j) = extempw(j)*zs(j,1,1,0)
         exw24(j) = exw24(j) + extempw(j)
         exwall(j) = exwall(j) + exw24(j)
      enddo
      
      do i = 1,nww4
         ii = iww4(i)
         if(localonoff(ii) .eq. 1)then
            iy = -iwy4(i)
            call bc2d_zdiagmatvec_add(nnodesp1,zs(0,1,iy,0),extempw,
     1           expnesw(0,4,ii))
         endif
      enddo
c
      ii = iwbig4(1)
      if(ii .gt. 0 .and. localonoff(ii) .eq. 1) then
         call bc2d_zdiagmatvec_add(nnodesp1,zs(0,2,-2,4),extempw,
     1        expneswbig(0,4,ii))
      endif

c     agglomeration and individual child interactions
c     done. send agglomerated expansions

c     take care of north work.

      do i = 1,nnall
         ii = inall(i)
         if(localonoff(ii) .eq. 1)then
            iy = -iynall(i)
            call bc2d_zdiagmatvec_add(nnodesp1,zs(0,3,iy,0),exnall,
     1           expnesw(0,1,ii))
         endif
      enddo
c     
      do i = 1,nn12
         ii = in12(i)
         if(localonoff(ii) .eq. 1)then
            iy = -iy12(i)
            call bc2d_zdiagmatvec_add(nnodesp1,zs(0,2,iy,0),exn12,
     1           expnesw(0,1,ii))
         endif
      enddo
c
      do i = 1, 3
         ii = inbig12(i)
         if (ii .gt. 0 .and. localonoff(ii) .eq. 1) then
            iy = -2*(i-2)
            call bc2d_zdiagmatvec_add(nnodesp1,zs(0,2,iy,1),exn12,
     1           expneswbig(0,1,ii))
         endif
      enddo

c     take care of south work

      do i = 1,nsall
         ii = isall(i)
         if(localonoff(ii) .eq. 1)then
            iy = iysall(i)
            call bc2d_zdiagmatvec_add(nnodesp1,zs(0,2,iy,0),exsall,
     1           expnesw(0,3,ii))
         endif
      enddo
c
      do i = 1,ns34
         ii = is34(i)
         if(localonoff(ii) .eq. 1)then
            iy = iy34(i)
            call bc2d_zdiagmatvec_add(nnodesp1,zs(0,1,iy,0),exs34,
     1           expnesw(0,3,ii))
         endif
      enddo
c
      do i = 1, 3
         ii = isbig34(i)
         if(ii .gt. 0 .and. localonoff(ii) .eq. 1)then
            iy = 2*(i-2)
            call bc2d_zdiagmatvec_add(nnodesp1,zs(0,2,iy,3),exs34,
     1           expneswbig(0,3,ii))
         endif
      enddo

c     take care of east work that wasn't already done
      
      do i = 1,neall
         ii = ieall(i)
         if(localonoff(ii) .eq. 1)then
            iy = iyeall(i)
            call bc2d_zdiagmatvec_add(nnodesp1,zs(0,3,iy,0),exeall,
     1           expnesw(0,2,ii))
         endif
      enddo
c
      do i = 1,ne13
         ii = ie13(i)
         if(localonoff(ii) .eq. 1)then
            iy = iy13(i)
            call bc2d_zdiagmatvec_add(nnodesp1,zs(0,2,iy,0),exe13,
     1           expnesw(0,2,ii))
         endif
      enddo
c
      ii = iebig13(1)
      if(ii .gt. 0 .and. localonoff(ii) .eq. 1) then
         call bc2d_zdiagmatvec_add(nnodesp1,zs(0,2,0,2),exe13,
     1        expneswbig(0,2,ii))
      endif
c
      
c     take care of west work that wasn't already done
c
      do i = 1,nwall
         ii = iwall(i)
         if(localonoff(ii) .eq. 1)then
            iy = -iywall(i)
            call bc2d_zdiagmatvec_add(nnodesp1,zs(0,2,iy,0),exwall,
     1           expnesw(0,4,ii))
         endif
      enddo
c
      do i = 1,nw24
         ii = iw24(i)
         if(localonoff(ii) .eq. 1)then
            iy = -iy24(i)
            call bc2d_zdiagmatvec_add(nnodesp1,zs(0,1,iy,0),exw24,
     1           expnesw(0,4,ii))
         endif
      enddo
c
      ii = iwbig24(1)
      if(ii .gt. 0 .and. localonoff(ii) .eq. 1)then
         call bc2d_zdiagmatvec_add(nnodesp1,zs(0,2,0,4),exw24,
     1        expneswbig(0,4,ii))
      endif
c
      return
      end
C
C
C
c***********************************************************************
      subroutine l2dpw_expshift(beta, nnodes, xnodes, spin, zshift)
c***********************************************************************
C     shifts an exponential expansion to a new center.
C
C     INPUT:
C
C     BETA          array of plane-wave coefficients being shifted
C     NNODES        order of plane-wave expansion
C     XNODES        nodes associated with plane wave expansion
C     SPIN          parameter determines what direction the 
C                   exponential expansions decay in
C     ZSHIFT        distance by which the expansion is to be shifted 
C                   (suitably scaled)
C
C     OUTPUT:
C
C     BETA          overwritten to account for the shift
C***********************************************************************
      implicit none
c-----Global variables
      integer nnodes
      real *8 xnodes(nnodes)
      complex *16 beta(0:nnodes), spin, zshift
c-----Local variables
      integer i
      complex *16 zspin
C 
c     Now adjust the coefficients to account for the shift.
C     Leave the first one unchanged.
C 
      zspin = spin*zshift
      do i = 1, nnodes
         beta(i) = beta(i) * cdexp(zspin*xnodes(i))
      enddo
      return
      end
C
C
C
c***********************************************************************
      subroutine l2dpw_expcoeff(betan,betae,betas,betaw,nnodes,mpole,
     2     nterms,comp,wnodes,sum,ratio)
c***********************************************************************
C     This subroutine generates four exponential expansions given
C     one multipole expansion as input.
C
C     INPUT:
C
C     NNODES      order of plane wave expansions
C     MPOLE       array or multipole coefficients
C     NTERMS      order of multipole expansions
C     COMP        precomputed array involving the ratio of
C                 quadrature weights and factorial terms
C     WNODES      array of exponential weights
C     SUM         precomputed real term
C     RATIO(I)     table with data WNODES(I) / XNODES(I)
C
C     OUTPUT:
C
C     BETAE       east exponential coefficients
C     BETAW       west exponential coefficients
C     BETAN       north exponential coefficients
C     BETAS       south exponential coefficients
C
C***********************************************************************
      implicit none
c-----Global variables
      integer  nnodes, nterms
      real *8  comp(0:nterms,nnodes)
      real *8  wnodes(nnodes)
      real *8  ratio(nnodes)
      real *8  sum, rscale
      complex *16  betae(0:nnodes), betaw(0:nnodes)
      complex *16  betan(0:nnodes), betas(0:nnodes)
      complex *16  mpole(0:nterms)
c-----Local variables
      integer  i, j
      complex *16  imag
      complex *16 sum1, sum2, sum5, sum6
      complex *16 sum7, sum8
      complex *16 temp1
      real *8 rscalej
      data imag/(0.0d0, 1.0d0)/
c
      temp1 = mpole(0)*sum
      betae(0) = temp1
      betaw(0) = temp1
      betas(0) = temp1
      betan(0) = temp1
      do i = 1, nnodes

c     1/z^j contribution
         
         sum1 = 0.0d0
         sum2 = 0.0d0
         sum5 = 0.0d0
         sum6 = 0.0d0
         do j = 1, nterms-3, 4
            sum2 = sum2 + comp(j-1,i) * mpole(j)
            sum1 = sum1 + comp(j,i)   * mpole(j+1)
            sum5 = sum5 + comp(j+1,i) * mpole(j+2)
            sum6 = sum6 + comp(j+2,i) * mpole(j+3)
         enddo
         sum8 = imag*(sum2 - sum5)
         sum7 = sum6 - sum1
         sum1 = sum1 + sum6
         sum2 = sum2 + sum5
         betae(i) = wnodes(i) * (sum1 + sum2)
         betaw(i) = wnodes(i) * (sum1 - sum2)
         betas(i) = wnodes(i) * (sum7 + sum8)
         betan(i) = wnodes(i) * (sum7 - sum8)

c     log contribution
         
         temp1 =  mpole(0)*ratio(i)
         betae(i) = betae(i) - temp1
         betaw(i) = betaw(i) - temp1
         betas(i) = betas(i) - temp1
         betan(i) = betan(i) - temp1
         
      enddo

      return
      end

C
C**********************************************************************
c     subroutine childpar
C**********************************************************************
C     merges the multipole expansions
C     of four child boxes together to form the multipole expansion
C     of the parent.  
C     The algorithm works by passing through the loop mod 4 and
C     rescaling the new coefficients so that they are on the parent's
C     level.
C
C     INPUT:
C 
C     MPOLE1 multipole coefficients for bottom left corner child
C     MPOLE2 multipole coefficients for bottom right corner child
C     MPOLE3 multipole coefficients for upper left corner child
C     MPOLE4 multipole coefficients for upper right corner child
C     NTERMS is the order of the multipole expansions
C     C      array of binomial coefficients 
C
C     OUTPUT:
C 
C     MPOLEPAR multipole coefficients of the parent box
C  
C     To AVOID REPEATED ALLOCS, many expansions are limited in 
C     length to 60. This is sufficient for precisions up to fifteen
C     digits.
C
C**********************************************************************
      subroutine l2dchildpar(mpolepar, mpole1, mpole2, mpole3,
     1                    mpole4, nterms, c)
      implicit none
c-----Global variables
      integer nterms
      real *8 c(2*nterms, *), rscale1, rscale2
      complex *16 mpole1(0:nterms), mpole2(0:nterms)
      complex *16 mpole3(0:nterms), mpole4(0:nterms)
      complex *16 mpolepar(0:nterms)
c-----Local variables
      integer  i,m,k,m1,m2,m3,k1,k2,k3
      real *8  xntemp, xntemp1
      complex *16 atemp(0:100),atemp2(0:100),atemp3(0:100),atemp4(0:100)
      complex *16  b1(100)
      complex *16  z00,z0p(0:100),z0p2(0:100)
      complex *16  cd,cdd,z0,imag
      complex *16  temp1,temp2,temp3,temp4
      data imag/(0.0d0,1.0d0)/
C
      z0 = (-.5d0 + .5d0*imag)
      z00=z0/2d0
      cd=1d0 / z0
      z0 = z00
      cdd=cd
      z0p(0) =1.0D0
      z0p2(0)=1.0D0
      do i=1,nterms
         z0p(i)  = z00
         z0p2(i) = cdd
         cdd = cdd*cd
         z00 = z00*z0
      end do
      do i=1,nterms
         b1(i)= 0.0d0
      end do
C
      
      mpolepar(0) = mpolepar(0) + mpole3(0)
     1      + mpole4(0) + mpole2(0) + mpole1(0)
      do m = 1, nterms
         temp1 = mpole3(m) + mpole2(m)
         temp2 = mpole3(m) - mpole2(m)
         temp3 = mpole4(m) + mpole1(m)
         temp4 = imag*(mpole4(m) - mpole1(m))
         atemp(m) =  (temp1 + temp3)*z0p2(m)
         atemp2(m) = (temp1 - temp3)*z0p2(m)
         atemp3(m) = (temp2 - temp4)*z0p2(m)
         atemp4(m) = (temp2 + temp4)*z0p2(m)
      end do
      xntemp = 2.0d0
      do m=1,nterms-3,4
         m1 = m + 1
         m2 = m + 2
         m3 = m + 3
         do k=1,nterms-3,4
            k1 = k + 1
            k2 = k + 2
            k3 = k + 3
            b1(m)=b1(m) + atemp(k)*c(m,k)
     1       + atemp4(k1)*c(m,k1)
     2       + atemp2(k2)*c(m,k2)
     3       + atemp3(k3)*c(m,k3)
            b1(m1)=b1(m1) + atemp3(k)*c(m1,k)
     1       + atemp(k1)*c(m1,k1)
     2       + atemp4(k2)*c(m1,k2)
     3       + atemp2(k3)*c(m1,k3)
            b1(m2)=b1(m2) + atemp2(k)*c(m2,k)
     1       + atemp3(k1)*c(m2,k1)
     2       + atemp(k2)*c(m2,k2)
     3       + atemp4(k3)*c(m2,k3)
            b1(m3)=b1(m3) + atemp4(k)*c(m3,k)
     1       + atemp2(k1)*c(m3,k1)
     2       + atemp3(k2)*c(m3,k2)
     3       + atemp(k3)*c(m3,k3)
         end do
         temp1 = mpole3(0) - mpole2(0)
         temp2 = imag*(mpole4(0) - mpole1(0))
         temp3 = mpole3(0) + mpole2(0)
         temp4 = mpole4(0) + mpole1(0)
         mpolepar(m) = mpolepar(m) + z0p(m)*(b1(m)
     1       +(temp2 - temp1)/dble(m))
         mpolepar(m1) = mpolepar(m1) + z0p(m1)*(b1(m1)
     1       + (temp4 - temp3)/dble(m1))
         mpolepar(m2) = mpolepar(m2) + z0p(m2)*(b1(m2)
     1       -(temp1 + temp2)/dble(m2))
         mpolepar(m3) = mpolepar(m3) + z0p(m3)*(b1(m3)
     1       -(temp3 + temp4)/dble(m3))
    
      enddo
      return
      end

c**********************************************************************
c     subroutine parentchild
c**********************************************************************
C     shifts the local expansions from a parent box to its 
C     four children.  This is used in the downward pass. 
C
C     INPUT:
C 
C     BETAHATPAR multipole coefficients of the parent box
C     NTERMS     order of the multipole expansions
C     C          array of binomial coefficients 
C
C     OUTPUT:
C
C     BETAHAT1 multipole coefficients for child 1, etc.
c
c     update: 2023/05/24. fix child ordering 
C**********************************************************************
      subroutine l2dparchild(betahatpar,
     1      beta1hat, beta2hat, beta3hat,
     2      beta4hat, nterms, c)
      implicit none
c-----Global variables
      integer nterms, iswitch
      real *8 c(2*nterms,*), rscale1, rscale2
      complex *16 beta1hat(0:nterms), beta2hat(0:nterms)
      complex *16 beta3hat(0:nterms), beta4hat(0:nterms)
      complex *16 betahatpar(0:nterms)
c-----Local variables
      integer i, m, k
      integer k1,k2,k3,k4
      integer m1,m2,m3
      complex *16 temp1,temp2,temp3,temp4,temp5
      complex *16 temp6,temp7,temp8,temp9,temp10
      complex *16 z0,z00,z0p(0:100),z0p2(0:100),anew(0:100)
      complex *16 cd,cdd,imag
      data imag/(0.0d0, 1.0d0)/

C     
C     Initialize array containing
C     powers of the shift vector (-.25 + i*.25)
C     Set Z0P(i) equal to Z0^i and
C     Z0P2(i) equal to 1 / Z0^i:
C
      z0 = 2*(-.25d0 + imag*.25d0)
      z00= z0/2d0
      cdd = 1d0/(z0)

      z0 = z00
      cd = cdd
      
      z0p(0)=1.0d0
      z0p2(0) = 1.0d0
C
      do i=1,nterms
         z0p(i)=z00
         z00= z00*z0
         z0p2(i) = cdd
         cdd = cdd*cd
      enddo
C
c     Initialize all local expansions to 
C     zero and set up the array ANEW:
      do i=0, nterms
         beta1hat(i) = 0.0d0
         beta2hat(i) = 0.0d0
         beta3hat(i) = 0.0d0
         beta4hat(i) = 0.0d0
         anew(i) = betahatpar(i)*z0p(i)
      enddo
C
c     Now go through a mod 4 loop in which the
C     expansions are shifted:
      do k=0, nterms-3, 4
         k1 = k + 1
         k2 = k + 2
         k3 = k + 3
         k4 = k + 4
         do m=0, k
            m1 = m + 1
            temp5 = anew(k)*c(k1,m1)
            temp6 = anew(k1)*c(k2,m1)
            temp7 = anew(k2)*c(k3,m1)
            temp8 = anew(k3)*c(k4,m1)
            temp1 = temp5 + temp7
            temp2 = temp5 - temp7
            temp3 = temp6 + temp8
            temp4 = imag*(temp6 - temp8)
            beta3hat(m) = beta3hat(m) + temp1 + temp3
            beta4hat(m) = beta4hat(m) + temp2 - temp4
            beta2hat(m) = beta2hat(m) + temp1 - temp3
            beta1hat(m) = beta1hat(m) + temp2 + temp4
         enddo
         temp7 = anew(k2)*c(k3,k2)
         temp8 = anew(k3)*c(k4,k2)
         temp9 = imag*(anew(k1) - temp8)
         temp10 = anew(k1) + temp8
         beta3hat(k1) = beta3hat(k1) + temp7 + temp10
         beta4hat(k1) = beta4hat(k1) - temp9 - temp7
         beta2hat(k1) = beta2hat(k1) - temp10 + temp7
         beta1hat(k1) = beta1hat(k1) + temp9  - temp7
         temp8 = anew(k3)*c(k4,k3)
         temp9 = imag*temp8
         beta3hat(k2) = beta3hat(k2) + anew(k2) + temp8
         beta4hat(k2) = beta4hat(k2) - anew(k2) + temp9
         beta2hat(k2) = beta2hat(k2) + anew(k2) - temp8
         beta1hat(k2) = beta1hat(k2) - anew(k2) - temp9
         temp8 = imag*anew(k3)
         beta3hat(k3) = beta3hat(k3) + anew(k3)
         beta4hat(k3) = beta4hat(k3) + temp8
         beta2hat(k3) = beta2hat(k3) - anew(k3)
         beta1hat(k3) = beta1hat(k3) - temp8
      enddo
      do m = 1, nterms-3, 4
         m1 = m + 1
         m2 = m + 2
         m3 = m + 3
         temp1 = z0p2(m)*imag
         temp2 = z0p2(m2)*imag
         beta3hat(m)  = beta3hat(m)*z0p2(m)
         beta3hat(m1) = beta3hat(m1)*z0p2(m1)
         beta3hat(m2) = beta3hat(m2)*z0p2(m2)
         beta3hat(m3) = beta3hat(m3)*z0p2(m3)
         beta4hat(m)  =  beta4hat(m)*temp1
         beta4hat(m1) = -beta4hat(m1)*z0p2(m1)
         beta4hat(m2) = -beta4hat(m2)*temp2
         beta4hat(m3) =  beta4hat(m3)*z0p2(m3)
         beta2hat(m)  = -beta2hat(m)*z0p2(m)
         beta2hat(m1) =  beta2hat(m1)*z0p2(m1)
         beta2hat(m2) = -beta2hat(m2)*z0p2(m2)
         beta2hat(m3) =  beta2hat(m3)*z0p2(m3)
         beta1hat(m)  = -beta1hat(m)*temp1
         beta1hat(m1) = -beta1hat(m1)*z0p2(m1)
         beta1hat(m2) =  beta1hat(m2)*temp2
         beta1hat(m3) =  beta1hat(m3)*z0p2(m3)
      enddo
      return
      end



c**********************************************************************
      subroutine l2dpw_addexp(b,a,nterms)
C**********************************************************************
C     INPUT :     complex arrays B,A of length (0:NTERMS)
C     OUTPUT:     A is over written by (A+B).
C***********************************************************************
      implicit none
      integer i, nterms
      complex *16 b(0:1),a(0:1)
      do i = 0,nterms
         a(i) = a(i) + b(i)
      enddo
      return
      end

c**********************************************************************
      subroutine l2dpw_setexpzero(b,nterms)
C**********************************************************************
C     INPUT :     complex arrays B of length (0:NTERMS)
C     OUTPUT:     B is over written by zeros.
C***********************************************************************
      implicit none
      integer i, nterms
      complex *16 b(0:1)
      complex *16 zero
      data zero /(0.0d0,0.0d0)/
      do i = 0,nterms
         b(i) = zero
      enddo
      return
      end

C
c**********************************************************************
      subroutine l2dpw_rescaleta(b,rscale1,rscale2,nterms)
C**********************************************************************
C     INPUT :     complex array B of length (0:NTERMS)
C     OUTPUT:     B is over written by B(i) = B(i)*(rscale2/rscale1)^i.
C***********************************************************************
      implicit none
c     global
      integer i, nterms
      complex *16 b(0:1)
      real *8 rscale1, rscale2
c     local
      real *8 temp1, temp2, done
      data done / 1.0d0 /
      temp1 = (rscale2/rscale1)
      temp2 = done
      do i = 0,nterms
         b(i) = b(i)*temp2
         temp2 = temp2*temp1
      enddo
      return
      end
C
      
