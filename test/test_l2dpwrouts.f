      program test_l2dpwrouts
c
c
c     unit tests for level-restricted tree
c     
c      
      implicit real *8 (a-h,o-z)

      
      real *8, allocatable :: x(:,:), w(:), u(:,:), v(:,:)
      real *8, allocatable :: tleg(:), wleg(:), fcoefs(:,:), fvals(:,:)
      real *8, allocatable :: pols(:), c(:,:), comp(:,:), ratio(:)
      real *8, allocatable :: xnodes(:), wnodes(:), targs(:,:)
      real *8 :: pt(2), center(2), rpot1, rgrad1(2), rhess1(3),
     1     center2(2), tctr(2)
      complex *16, allocatable :: poly2mpmat(:,:),mpole(:,:),
     1     locexp(:,:),poly2pw(:,:,:),beta_all(:,:),beta_all_btos(:,:),
     2     poly2pw_btos(:,:,:)
      
      
      complex *16, allocatable, dimension(:) :: betan,betae,betas,betaw      
      complex *16, allocatable, dimension(:) :: betani,betaei,
     1     betasi,betawi,mpolepar
      complex *16, allocatable :: zs(:,:,:,:)
      complex *16 :: pot0, grad0(2), hess0(3)
      complex *16 :: pot1, grad1(2), hess1(3)
      complex *16 :: pot2, grad2(2), hess2(3)
      complex *16 :: pot3, grad3(2), hess3(3)
      complex *16 :: im

      real *8 :: gradtrue(2)
      
      data im /(0.0d0,1.0d0)/

      integer :: ipass(20)


      character :: type

      ntest = 0
      do i = 1,20
         ipass(i) = 0
      enddo
      
      iseed = 1234
      dburn = hkrand(iseed)
      
      call prini(6,13)

      ndeg = 15
      norder = ndeg+1
      npt = ndeg+1
      npt2 = npt*npt
      type = 'T'
      call legetens_npol_2d(ndeg,type,npol)

      allocate(x(2,npt2),u(npol,npt2),v(npt2,npol),w(npt2))
      itype = 2
      ldu = npol
      ldv = npt2
      call legetens_exps_2d(itype,npt,type,x,u,ldu,v,ldv,w)

      eps = 1d-12
      call l2dterms(eps,nterms,ier1)
      nterms = nterms+4-mod(nterms,4)
      call prinf('l2dterms, ier1 *',ier1,1)
      call prinf('nterms *',nterms,1)
      call lwtsexp2b_terms(eps,nnodes)
      call prinf('nnodes *',nnodes,1)
      allocate(c(2*nterms,2*nterms))
      allocate(comp(0:nterms,nnodes),ratio(nnodes))
      allocate(xnodes(nnodes),wnodes(nnodes))
      call lwtsexp2b(nnodes,xnodes,wnodes,errnodes)
      call l2dpw_precomp(comp,nterms,nnodes,xnodes,wnodes,ratio,c,xsum)


      allocate(poly2mpmat(0:nterms,npol))
      nleg = npt+14
      itype = 1
      alpha = 0.0d0
      beta = 0.0d0
      allocate(tleg(nleg),wleg(nleg))
      call legeexps(itype,nleg,tleg,alpha,beta,wleg)

      call cpu_time(t0)
      call pbox2d_genpoly2mpmat(poly2mpmat,nterms,npol,type,ndeg,
     1     nleg)
      call cpu_time(t1)

      call prin2('time for poly2mpmat *',t1-t0,1)

      allocate(fvals(npt2,5))
      
      do i = 1,npt2
         call testfun(x(1,i),fvals(i,1))

         pt(1) = x(1,i)/2 - 0.5d0
         pt(2) = x(2,i)/2 - 0.5d0
         call testfun(pt,fvals(i,1+1))
         pt(1) = x(1,i)/2 + 0.5d0
         pt(2) = x(2,i)/2 - 0.5d0
         call testfun(pt,fvals(i,1+2))
         pt(1) = x(1,i)/2 - 0.5d0
         pt(2) = x(2,i)/2 + 0.5d0
         call testfun(pt,fvals(i,1+3))
         pt(1) = x(1,i)/2 + 0.5d0
         pt(2) = x(2,i)/2 + 0.5d0
         call testfun(pt,fvals(i,1+4))
      enddo
      
      allocate(mpole(0:nterms,5),fcoefs(npol,5))
      rscale = 0.5d0
      do i = 1,5
         call bc2d_dmatvec(npol,npt2,u,fvals(1,i),fcoefs(1,i))
         m = nterms+1
         if (i .gt. 1) then
            rscales = rscale**2
         else
            rscales=1
         endif
         call bc2d_dscalzmatdvec(m,npol,rscales,poly2mpmat,
     1        fcoefs(1,i),mpole(0,i))

      enddo
      
c     test multipole expansion on parent and children
      
      pt(1) = 3d0
      pt(2) = 0.01d0

      ifgrad = 1
      ifhess = 1
      rscale2 = 2
      center(1)=0
      center(2)=0
      call l2dmpeval(rscale2,center,mpole,nterms,pt,pot1,ifgrad,grad1,
     1     ifhess,hess1)

      call prin2('pot1 *',pot1,2)

      pot2 = 0
      grad2(1) = 0
      grad2(2) = 0
      rscale = 1d0
      center(1)=-0.5d0
      center(2)=-0.5d0
      call l2dmpeval(rscale,center,mpole(0,1+1),nterms,pt,pot0,
     1     ifgrad,grad0,ifhess,hess0)
      pot2 = pot2+pot0
      grad2(1) = grad2(1)+grad0(1)
      grad2(2) = grad2(2)+grad0(2)

      center(1)=0.5d0
      center(2)=-0.5d0
      call l2dmpeval(rscale,center,mpole(0,1+2),nterms,pt,pot0,
     1     ifgrad,grad0,ifhess,hess0)
      pot2 = pot2+pot0
      grad2(1) = grad2(1)+grad0(1)
      grad2(2) = grad2(2)+grad0(2)
      
      center(1)=-0.5d0
      center(2)=0.5d0
      call l2dmpeval(rscale,center,mpole(0,1+3),nterms,pt,pot0,
     1     ifgrad,grad0,ifhess,hess0)
      pot2 = pot2+pot0
      grad2(1) = grad2(1)+grad0(1)
      grad2(2) = grad2(2)+grad0(2)
      
      center(1)=0.5d0
      center(2)=0.5d0
      call l2dmpeval(rscale,center,mpole(0,1+4),nterms,pt,pot0,
     1     ifgrad,grad0,ifhess,hess0)
      pot2 = pot2+pot0
      grad2(1) = grad2(1)+grad0(1)
      grad2(2) = grad2(2)+grad0(2)

      call prin2('pot1 *',pot1,2)
      call prin2('pot2 *',pot2,2)

      call prin2('grad1 *',grad1,4)
      call prin2('grad2 *',grad2,4)

      call testconv(pt,pottrue)
      call testconvgrad(pt,gradtrue)

      call prin2('pottrue *',pottrue,1)
      call prin2('gradtrue *',gradtrue,2)

      write(*,*) abs(pottrue-real(pot1))
      write(*,*) abs(pottrue-real(pot2))
      write(*,*) abs(gradtrue(1)-real(grad1(1)))
      write(*,*) abs(gradtrue(2)-real(grad1(2)))
      write(*,*) abs(gradtrue(1)-real(grad2(1)))
      write(*,*) abs(gradtrue(2)-real(grad2(2)))

c     test multipole expansion merge

      allocate(mpolepar(0:nterms))
      do i = 0,nterms
         mpolepar(i)=0
      enddo
      rscale1 = 1d0
      rscale2 = 2d0
      call l2dchildpar(mpolepar, mpole(0,2), mpole(0,3), mpole(0,4),
     1     mpole(0,5), nterms, c)
      

      center(1)=0
      center(2)=0
      call l2dmpeval(rscale2,center,mpolepar,nterms,pt,pot3,
     1     ifgrad,grad3,ifhess,hess3)

      write(*,*) abs(pottrue-real(pot3))

c     form planewaves

      
      allocate(betan(0:nnodes),betae(0:nnodes),betas(0:nnodes),
     1     betaw(0:nnodes))

      sum1 = xsum+dlog(rscale2)
      call l2dpw_expcoeff(betan,betae,betas,betaw,nnodes,mpolepar,
     2     nterms,comp,wnodes,sum1,ratio)

      allocate(poly2pw(0:nnodes,4,npol))
      call cpu_time(t00)
      call pbox2d_genpoly2pwmat(norder,npol,type,nnodes,
     1     wnodes,xnodes,poly2pw)
      call cpu_time(t11)

      write(*,*) 'time for generating pw ', t11-t00
      
      allocate(poly2pw_btos(0:nnodes,4,npol))
      call pbox2d_genpoly2pwmat_btosfar(norder,npol,type,nnodes,
     1     wnodes,xnodes,poly2pw_btos)

      allocate(beta_all(0:nnodes,4),beta_all_btos(0:nnodes,4))
      m = 4*(nnodes+1)
      xlength = 2
      scal = (xlength/2)**2
      call bc2d_dscalzmatdvec(m,npol,scal,poly2pw,
     1     fcoefs(1,1),beta_all)
      pi2 = 8*atan(1d0)
      dtmp = -4*scal*(xsum+log(xlength))*fcoefs(1,1)/pi2
      write(*,*) dtmp/betan(0)
      beta_all(0,1) = dtmp
      beta_all(0,2) = dtmp
      beta_all(0,3) = dtmp
      beta_all(0,4) = dtmp

      m = 4*(nnodes+1)
      xlength = 2
      scal = (xlength/2)**2
      call bc2d_dscalzmatdvec(m,npol,scal,poly2pw_btos,
     1     fcoefs(1,1),beta_all_btos)
      pi2 = 8*atan(1d0)
      dtmp = -4*scal*(xsum+log(xlength/2))*fcoefs(1,1)/pi2
      write(*,*) dtmp/betan(0)
      beta_all_btos(0,1) = dtmp
      beta_all_btos(0,2) = dtmp
      beta_all_btos(0,3) = dtmp
      beta_all_btos(0,4) = dtmp

      center(1)=0
      center(2)=0
      ifpgh = 2
      rscale = rscale2

      pt(1)=0.01d0
      pt(2)=3d0
      idir=1
      call l2dpw_eval(pt,xnodes,nnodes,rscale,center,betan,idir,
     1     ifpgh,rpot1,rgrad1,rhess1)
      call testconv(pt,pottrue)
      call testconvgrad(pt,gradtrue)
      write(*,*) 'error in planewave (north)'
      write(*,*) abs(pottrue-rpot1)
      write(*,*) abs(gradtrue(1)-rgrad1(1))
      write(*,*) abs(gradtrue(2)-rgrad1(2))

      ntest = ntest+1
      if (abs(pottrue-rpot1) .lt. 1d-12) ipass(ntest) = 1
      
      idir=1
      call l2dpw_eval(pt,xnodes,nnodes,rscale,center,beta_all(0,1),idir,
     1     ifpgh,rpot1,rgrad1,rhess1)
      call testconv(pt,pottrue)
      call testconvgrad(pt,gradtrue)
      write(*,*) 'error in planewave (formed directly, north)'
      write(*,*) abs(pottrue-rpot1)
      write(*,*) abs(gradtrue(1)-rgrad1(1))
      write(*,*) abs(gradtrue(2)-rgrad1(2))

      ntest = ntest+1
      if (abs(pottrue-rpot1) .lt. 1d-12) ipass(ntest) = 1

      idir=1
      pt(1)=0.01d0
      pt(2)=2d0
      center(1) = -0.5d0
      center(2) = -0.5d0
      rscale = 1
      call l2dpw_eval(pt,xnodes,nnodes,rscale,center,beta_all_btos(0,1),
     1     idir,
     1     ifpgh,rpot1,rgrad1,rhess1)
      call testconv(pt,pottrue)
      call testconvgrad(pt,gradtrue)
      write(*,*) 'error in planewave (btosfar, north)'
      write(*,*) abs(pottrue-rpot1)
      write(*,*) abs(gradtrue(1)-rgrad1(1))
      write(*,*) abs(gradtrue(2)-rgrad1(2))
      
      ntest = ntest+1
      if (abs(pottrue-rpot1) .lt. 1d-12) ipass(ntest) = 1
      
      pt(2)=0.01d0
      pt(1)=3d0
      idir=2
      rscale = 2
      center(1)=0
      center(2)=0
      
      call l2dpw_eval(pt,xnodes,nnodes,rscale,center,betae,idir,
     1     ifpgh,rpot1,rgrad1,rhess1)
      call testconv(pt,pottrue)
      call testconvgrad(pt,gradtrue)
      write(*,*) 'error in planewave (east)'
      write(*,*) abs(pottrue-rpot1)
      write(*,*) abs(gradtrue(1)-rgrad1(1))
      write(*,*) abs(gradtrue(2)-rgrad1(2))

      ntest = ntest+1
      if (abs(pottrue-rpot1) .lt. 1d-12) ipass(ntest) = 1
      
      idir=2
      pt(2)=0.01d0
      pt(1)=2d0
      center(1) = -0.5d0
      center(2) = -0.5d0
      rscale = 1
      call l2dpw_eval(pt,xnodes,nnodes,rscale,center,beta_all_btos(0,2),
     1     idir,
     1     ifpgh,rpot1,rgrad1,rhess1)
      call testconv(pt,pottrue)
      call testconvgrad(pt,gradtrue)
      write(*,*) 'error in planewave (btosfar, east)'
      write(*,*) abs(pottrue-rpot1)
      write(*,*) abs(gradtrue(1)-rgrad1(1))
      write(*,*) abs(gradtrue(2)-rgrad1(2))

      ntest = ntest+1
      if (abs(pottrue-rpot1) .lt. 1d-12) ipass(ntest) = 1
      
      pt(1)=0.01d0
      pt(2)=-3d0
      idir=3
      rscale = 2
      center(1)=0
      center(2)=0
      
      call l2dpw_eval(pt,xnodes,nnodes,rscale,center,betas,idir,
     1     ifpgh,rpot1,rgrad1,rhess1)
      call testconv(pt,pottrue)
      call testconvgrad(pt,gradtrue)
      write(*,*) 'error in planewave (south)'
      write(*,*) abs(pottrue-rpot1)
      write(*,*) abs(gradtrue(1)-rgrad1(1))
      write(*,*) abs(gradtrue(2)-rgrad1(2))

      ntest = ntest+1
      if (abs(pottrue-rpot1) .lt. 1d-12) ipass(ntest) = 1

      idir=3
      pt(1)=0.01d0
      pt(2)=-2d0
      center(1) = -0.5d0
      center(2) = -0.5d0
      rscale = 1
      call l2dpw_eval(pt,xnodes,nnodes,rscale,center,beta_all_btos(0,3),
     1     idir,
     1     ifpgh,rpot1,rgrad1,rhess1)
      call testconv(pt,pottrue)
      call testconvgrad(pt,gradtrue)
      write(*,*) 'error in planewave (btosfar, south)'
      write(*,*) abs(pottrue-rpot1)
      write(*,*) abs(gradtrue(1)-rgrad1(1))
      write(*,*) abs(gradtrue(2)-rgrad1(2))

      ntest = ntest+1
      if (abs(pottrue-rpot1) .lt. 1d-12) ipass(ntest) = 1
      
      
      pt(2)=0.01d0
      pt(1)=-3d0
      idir=4
      center(1) = 0
      center(2) = 0
      rscale = 2
      call l2dpw_eval(pt,xnodes,nnodes,rscale,center,betaw,idir,
     1     ifpgh,rpot1,rgrad1,rhess1)
      call testconv(pt,pottrue)
      call testconvgrad(pt,gradtrue)

      write(*,*) 'error in planewave (west)'
      write(*,*) abs(pottrue-rpot1)
      write(*,*) abs(gradtrue(1)-rgrad1(1))
      write(*,*) abs(gradtrue(2)-rgrad1(2))

      ntest = ntest+1
      if (abs(pottrue-rpot1) .lt. 1d-12) ipass(ntest) = 1
      
      idir=4
      pt(2)=0.01d0
      pt(1)=-2d0
      center(1) = -0.5d0
      center(2) = -0.5d0
      rscale = 1
      call l2dpw_eval(pt,xnodes,nnodes,rscale,center,beta_all_btos(0,4),
     1     idir,
     1     ifpgh,rpot1,rgrad1,rhess1)
      call testconv(pt,pottrue)
      call testconvgrad(pt,gradtrue)
      write(*,*) 'error in planewave (btosfar, west)'
      write(*,*) abs(pottrue-rpot1)
      write(*,*) abs(gradtrue(1)-rgrad1(1))
      write(*,*) abs(gradtrue(2)-rgrad1(2))
      
      ntest = ntest+1
      if (abs(pottrue-rpot1) .lt. 1d-12) ipass(ntest) = 1

c     test shift and pw to local

      allocate(zs(0:nnodes,-3:3,-3:3,0:4))
      call cpu_time(t0)
      call l2dpw_mkshifts(xnodes,nnodes,zs)
      call cpu_time(t1)
      
      call prin2('time for mkshifts *',t1-t0,1)

      allocate(betani(0:nnodes),betaei(0:nnodes),betasi(0:nnodes),
     1     betawi(0:nnodes))

      ix = 2
      iy = 0
      do i = 0,nnodes
         betani(i)=0
         betaei(i)=0
         betasi(i)=0
         betawi(i)=0

         betani(i) = betan(i)*zs(i,ix,iy,0)
      enddo

      center(1)=0
      center(2)=4
      ifpgh = 2
      rscale = rscale2

      pt(1)=0.01d0
      pt(2)=3d0
      idir=1
      call l2dpw_eval(pt,xnodes,nnodes,rscale,center,betani,idir,
     1     ifpgh,rpot1,rgrad1,rhess1)
      call testconv(pt,pottrue)
      call testconvgrad(pt,gradtrue)
      write(*,*) 'error in shifted planewave (north)'
      write(*,*) abs(pottrue-rpot1)
      write(*,*) abs(gradtrue(1)-rgrad1(1))
      write(*,*) abs(gradtrue(2)-rgrad1(2))
      
      ntest = ntest+1
      if (abs(pottrue-rpot1) .lt. 1d-12) ipass(ntest) = 1
      
c     convert to local and eval

      allocate(locexp(0:nterms,5))
      call l2dpw_exp4local(locexp(0,5),nterms,nnodes,betani,
     1     betaei,betasi,betawi,comp)
      
      call l2dtaeval(rscale,center,locexp(0,5),nterms,pt,
     1     pot1,ifgrad,grad1,ifhess,hess1)

      write(*,*)
     1     'error in shifted planewave converted to local (north)'
      write(*,*) abs(pottrue-real(pot1))
      write(*,*) abs(gradtrue(1)-real(grad1(1)))
      write(*,*) abs(gradtrue(2)-real(grad1(2)))

      ntest = ntest+1
      if (abs(pottrue-rpot1) .lt. 1d-12) ipass(ntest) = 1
      
c     send to children, eval

      rscale2 = 2
      rscale1 = 1
      call l2dparchild(locexp(0,5),locexp(0,1),locexp(0,2),locexp(0,3),
     1     locexp(0,4), nterms, c)

c     eval in child 1
      
      center(1)=0-0.5
      center(2)=4-0.5
      center2(1) = 0
      center2(2) = 0
      
      rscale2 = 2
      rscale1 = 1
      
      call l2dtaeval(rscale1,center,locexp(0,1),nterms,pt,
     1     pot1,ifgrad,grad1,ifhess,hess1)      

      write(*,*)
     1     'error in local passed to child'
      write(*,*) abs(pottrue-real(pot1))
      write(*,*) abs(gradtrue(1)-real(grad1(1)))
      write(*,*) abs(gradtrue(2)-real(grad1(2)))

      ntest = ntest+1
      if (abs(pottrue-rpot1) .lt. 1d-12) ipass(ntest) = 1
      

c
c     similar for btosfar
c      
      
      ix = 2
      iy = -1
      do i = 0,nnodes
         betani(i)=0
         betaei(i)=0
         betasi(i)=0
         betawi(i)=0

         betawi(i) = beta_all_btos(i,4)*zs(i,ix,iy,0)
      enddo

      center(1)=-2-0.5d0
      center(2)=1-0.5d0
      ifpgh = 2
      rscale = 1
      pt(1)=-2d0
      pt(2)=0.51d0
      idir=4
      call l2dpw_eval(pt,xnodes,nnodes,rscale,center,betawi,idir,
     1     ifpgh,rpot1,rgrad1,rhess1)
      call testconv(pt,pottrue)
      call testconvgrad(pt,gradtrue)
      write(*,*) 'error in shifted planewave (west)'
      write(*,*) abs(pottrue-rpot1)
      write(*,*) abs(gradtrue(1)-rgrad1(1))
      write(*,*) abs(gradtrue(2)-rgrad1(2))
      
      
      ntest = ntest+1
      if (abs(pottrue-rpot1) .lt. 1d-12) ipass(ntest) = 1
      
c     convert to local and eval

      call l2dpw_exp4local(locexp(0,5),nterms,nnodes,betani,
     1     betaei,betasi,betawi,comp)
      
      call l2dtaeval(rscale,center,locexp(0,5),nterms,pt,
     1     pot1,ifgrad,grad1,ifhess,hess1)

      write(*,*)
     1     'error in shifted planewave converted to local (west)'
      write(*,*) abs(pottrue-real(pot1))
      write(*,*) abs(gradtrue(1)-real(grad1(1)))
      write(*,*) abs(gradtrue(2)-real(grad1(2)))

      ntest = ntest+1
      if (abs(pottrue-rpot1) .lt. 1d-12) ipass(ntest) = 1
      
      ix = 2
      iy = 1
      do i = 0,nnodes
         betani(i)=0
         betaei(i)=0
         betasi(i)=0
         betawi(i)=0

         betasi(i) = beta_all_btos(i,3)*zs(i,ix,iy,0)
      enddo

      center(1)=1-0.5d0
      center(2)=-2-0.5d0
      ifpgh = 2
      rscale = 1
      pt(1)=0.51d0
      pt(2)=-2
      idir=3
      call l2dpw_eval(pt,xnodes,nnodes,rscale,center,betasi,idir,
     1     ifpgh,rpot1,rgrad1,rhess1)
      call testconv(pt,pottrue)
      call testconvgrad(pt,gradtrue)
      write(*,*) 'error in shifted planewave (south)'
      write(*,*) abs(pottrue-rpot1)
      write(*,*) abs(gradtrue(1)-rgrad1(1))
      write(*,*) abs(gradtrue(2)-rgrad1(2))
      
      ntest = ntest+1
      if (abs(pottrue-rpot1) .lt. 1d-12) ipass(ntest) = 1
      
c     convert to local and eval

      call l2dpw_exp4local(locexp(0,5),nterms,nnodes,betani,
     1     betaei,betasi,betawi,comp)
      
      call l2dtaeval(rscale,center,locexp(0,5),nterms,pt,
     1     pot1,ifgrad,grad1,ifhess,hess1)

      write(*,*)
     1     'error in shifted planewave converted to local (south)'
      write(*,*) abs(pottrue-real(pot1))
      write(*,*) abs(gradtrue(1)-real(grad1(1)))
      write(*,*) abs(gradtrue(2)-real(grad1(2)))

      ntest = ntest+1
      if (abs(pottrue-rpot1) .lt. 1d-12) ipass(ntest) = 1

      npass = 0
      do i = 1,ntest
         npass = npass + ipass(i)
      enddo
      
      open(unit=33,file='print_testres.txt',access='append')
      write(33,'(a,i2,a,i2,a)') 'Successfully completed ',npass,
     1   ' out of ',ntest,' tests in l2dpwrouts testing suite'
      close(33)
      
      stop
      end

      subroutine testconv(x,val)
      implicit real *8 (a-h,o-z)
      real *8, allocatable :: tleg(:), wleg(:)
      real *8 pt(2), x(2)
      
      pi4 = 16*atan(1.0d0)
      nleg = 50
      allocate(tleg(nleg),wleg(nleg))
      itype = 1
      alpha = 0.0d0
      beta = 0.0d0
      call legeexps(itype,nleg,tleg,alpha,beta,wleg)

      val = 0
      do j = 1,nleg
         do i = 1,nleg
            wij = wleg(i)*wleg(j)
            pt(1) = tleg(i)
            pt(2) = tleg(j)
            call testfun(pt,f)
            val = val-f*log((pt(1)-x(1))**2+(pt(2)-x(2))**2)/pi4*wij
         enddo
      enddo
      
      return
      end
      
      subroutine testconvgrad(x,grad)
      implicit real *8 (a-h,o-z)
      real *8, allocatable :: tleg(:), wleg(:)
      real *8 pt(2), x(2), grad(2), diff(2)
      
      pi2 = 8*atan(1.0d0)
      nleg = 50
      allocate(tleg(nleg),wleg(nleg))
      itype = 1
      alpha = 0.0d0
      beta = 0.0d0
      call legeexps(itype,nleg,tleg,alpha,beta,wleg)

      grad(1) = 0
      grad(2) = 0
      do j = 1,nleg
         do i = 1,nleg
            wij = wleg(i)*wleg(j)
            pt(1) = tleg(i)
            pt(2) = tleg(j)
            diff(1) = tleg(i)-x(1)
            diff(2) = tleg(j)-x(2)
            call testfun(pt,f)
            grad(1) = grad(1)+wij*f*diff(1)/(diff(1)**2+diff(2)**2)/pi2
            grad(2) = grad(2)+wij*f*diff(2)/(diff(1)**2+diff(2)**2)/pi2
         enddo
      enddo
      
      return
      end
      
      
      subroutine testfun(x,f)
      implicit real *8 (a-h,o-z)
      real *8 :: x(2), f

      xx = x(1)
      yy = x(2)
      f = cos(xx - yy) + (xx-0.2d0)**4 + (yy+0.1d0)**3+0.1d0

      return
      end
