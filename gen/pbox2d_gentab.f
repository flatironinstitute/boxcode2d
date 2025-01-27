cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     author: Travis Askham
c
c     this file is being released under an Apache License (version 2.0
c     See LICENSE in home directory for details 
c     

c
c     subroutines for generating tables for poisson interactions 
c
      
      subroutine pbox2d_genreftab(pttype,polytype,norder,npoly,npt,
     1     tabcol,tabstob,tabbtos,gradtabcol,gradtabstob,gradtabbtos,
     2     maxcls,ier)
      implicit real *8 (a-h,o-z)
      character :: pttype, polytype
      integer :: ndeg, npoly, npt, maxcls
      real *8 :: tabcol(npt,npoly,*), tabstob(npt,npoly,*),
     1     tabbtos(npt,npoly,*)
      real *8 :: gradtabcol(2,npt,npoly,*), gradtabstob(2,npt,npoly,*),
     1     gradtabbtos(2,npt,npoly,*)
c     local
      real *8, allocatable :: pt(:,:), vals(:), grads(:,:)
      real *8 :: ptsc(2), ctrs(2,3)
      
      ier = 0

      ndeg = norder-1
      call legetens_npol_2d(ndeg,polytype,npoly2)
      if (npoly .ne. npoly2) then
         ier = 1
         return
      endif
      
      if (pttype .eq. 'F' .or. pttype .eq. 'f') then
         npt2 = (ndeg+1)**2
         if (npt .ne. npt2) then
            ier = 2
            return
         endif
         allocate(pt(2,npt))
         itype = 0
         ldu = 1
         ldv = 1
         npt0 = ndeg+1
         call legetens_exps_2d(itype,npt0,polytype,pt,u,ldu,v,ldv,wts)
      endif


      call pbox2d_genfun_init0(polytype,ndeg,npoly)
      call pbox2d_genfunx_init0(polytype,ndeg,npoly)
      call pbox2d_genfuny_init0(polytype,ndeg,npoly)
      
      allocate(vals(npoly),grads(2,npoly))

c     colleagues

      ier2 = 0

      
      ctrs(1,1) = -2
      ctrs(2,1) = -2
      ctrs(1,2) = 0
      ctrs(2,2) = -2
      ctrs(1,3) = 0
      ctrs(2,3) = 0
      dlen = 2
      do l = 1,3
         write(*,*) 'computing colleague ref table ', l
c$OMP PARALLEL DO DEFAULT(shared) PRIVATE(ptsc,vals,grads)         
c$OMP1 PRIVATE(j,ier1) REDUCTION(+:ier2)
         do i = 1,npt
            ptsc(1) = ctrs(1,l) + pt(1,i)*dlen/2
            ptsc(2) = ctrs(2,l) + pt(2,i)*dlen/2
            call pbox2d_genreftab_int1(ptsc,polytype,ndeg,npoly,vals,
     1           grads,maxcls,ier1)
            if (ier1 .ne. 0) then
               ier2 = ier2+1
            endif
            do j = 1,npoly
               tabcol(i,j,l) = vals(j)
               gradtabcol(1,i,j,l) = grads(1,j)
               gradtabcol(2,i,j,l) = grads(2,j)
            enddo
         enddo
c$OMP END PARALLEL DO
         if (ier2 .ne. 0) then
            ier = 4
            return
         endif
      enddo
      
c     stob
      
      ctrs(1,1) = -3
      ctrs(2,1) = -3
      ctrs(1,2) = -1
      ctrs(2,2) = -3
      dlen = 4
      do l = 1,2
         write(*,*) 'computing stob ref table ', l         
c$OMP PARALLEL DO DEFAULT(shared) PRIVATE(ptsc,vals,grads)         
c$OMP1 PRIVATE(j,ier1) REDUCTION(+:ier2)
         do i = 1,npt
            ptsc(1) = ctrs(1,l) + pt(1,i)*dlen/2
            ptsc(2) = ctrs(2,l) + pt(2,i)*dlen/2
            call pbox2d_genreftab_int1(ptsc,polytype,ndeg,npoly,vals,
     1           grads,maxcls,ier1)
            if (ier1 .ne. 0) then
               ier2 = ier2+1
            endif
            do j = 1,npoly
               tabstob(i,j,l) = vals(j)
               gradtabstob(1,i,j,l) = grads(1,j)
               gradtabstob(2,i,j,l) = grads(2,j)
            enddo
         enddo
c$OMP END PARALLEL DO
         if (ier2 .ne. 0) then
            ier = 4
            return
         endif
      enddo
      
c     btos
      
      ctrs(1,1) = -1.5d0
      ctrs(2,1) = -1.5d0
      ctrs(1,2) = -0.5d0
      ctrs(2,2) = -1.5d0
      dlen = 1
      do l = 1,2
         write(*,*) 'computing btos ref table ', l         
c$OMP PARALLEL DO DEFAULT(shared) PRIVATE(ptsc,vals,grads)         
c$OMP1 PRIVATE(j,ier1) REDUCTION(+:ier2)
         do i = 1,npt
            ptsc(1) = ctrs(1,l) + pt(1,i)*dlen/2
            ptsc(2) = ctrs(2,l) + pt(2,i)*dlen/2
            call pbox2d_genreftab_int1(ptsc,polytype,ndeg,npoly,vals,
     1           grads,maxcls,ier1)
            if (ier1 .ne. 0) then
               ier2 = ier2+1
            endif
            do j = 1,npoly
               tabbtos(i,j,l) = vals(j)
               gradtabbtos(1,i,j,l) = grads(1,j)
               gradtabbtos(2,i,j,l) = grads(2,j)
            enddo
         enddo
c$OMP END PARALLEL DO
         if (ier2 .ne. 0) then
            ier = 4
            return
         endif
      enddo
      
      return
      end

      subroutine pbox2d_genreftab_int1(pt,polytype,ndeg,npoly,
     1     vals,grads,maxcls,ier1)
      implicit real *8 (a-h,o-z)
      real *8 :: vals(npoly), grads(2,npoly), pt(2)
      character :: polytype
      integer ndeg, npoly
c     local
      external pbox2d_genfun, pbox2d_genfunx, pbox2d_genfuny
      integer key, n, nf, ndim, mincls, maxcls, ifail, neval, nw
      parameter (ndim = 2)
      real *8 a(ndim), b(ndim)
      real *8 absreq, relreq
      real *8, allocatable :: absest(:), finest(:), work(:)

      ier1=0
      
      mincls = 0
      num=65
      maxsub = maxcls/(2*num) + 1
      nw = maxsub*(2*ndim+2*npoly+2) + 17*npoly + 1 + 1000000
      
      allocate(absest(npoly),finest(npoly),work(nw))

      do n = 1,ndim
         a(n) = -1d0
         b(n) =  1d0
      enddo

      key = 0
      absreq = 1d-14
      relreq = 1d-14

      call pbox2d_genfun_init1(pt)
      call pbox2d_genfunx_init1(pt)
      call pbox2d_genfuny_init1(pt)

      call dcuhre(ndim, npoly, a,b,mincls, maxcls, pbox2d_genfun, 
     1     absreq, relreq, key, nw, 0, finest, absest, neval,
     2     ier1, work)
      
      
      do i = 1,npoly
         vals(i) = finest(i)
      enddo
      call dcuhre(ndim, npoly, a, b,mincls, maxcls,pbox2d_genfunx, 
     1     absreq, relreq, key, nw, 0, finest, absest, neval,
     2     ier1, work)
      
      do i = 1,npoly
         grads(1,i) = finest(i)
      enddo
      call dcuhre(ndim, npoly,a,b,mincls, maxcls, pbox2d_genfuny, 
     1     absreq, relreq, key, nw, 0, finest, absest, neval,
     2     ier1, work)
      
      do i = 1,npoly
         grads(2,i) = finest(i)
      enddo
      
      return
      end


      subroutine pbox2d_genfun(ndim, z, nfun, f)
      implicit none
c$    integer omp_get_thread_num            
      integer ndim, nfun, ndeg, npoly, id
      real *8 z(ndim), f(nfun), rx, ry, pi4
      real *8 rr,reps,dgreen
      real *8 xtarg,ytarg
      real *8 :: pt(2,200), pt0(2)
      character :: polytype, polytype0
      integer ndeg0, npoly0, i
      save pt, ndeg, npoly, pi4, polytype

      id = 1
c$    id = omp_get_thread_num()+1   

      rx = z(1) - pt(1,id)
      ry = z(2) - pt(2,id)
      
      rr = rx*rx +ry*ry
      reps = 1.0d-30
      if (rr.gt.reps) then
         dgreen = -dlog(rr)/pi4
         call legetens_pols_2d(z,ndeg,polytype,f)
         do i = 1,npoly
            f(i) = dgreen * f(i)
         enddo
      else
         do i = 1,npoly
            f(i) = 0.0d0
         enddo
      endif

      return
      
      entry pbox2d_genfun_init0(polytype0,ndeg0,npoly0)
      pi4 = 16.0d0*datan(1.0d0)
      ndeg=ndeg0
      npoly=npoly0
      polytype=polytype0
      return

      entry pbox2d_genfun_init1(pt0)
      id = 1
c$    id = omp_get_thread_num()+1
      pt(1,id)=pt0(1)
      pt(2,id)=pt0(2)
      return

      
      end

      subroutine pbox2d_genfunx(ndim, z, nfun, f)
      implicit none
c$    integer omp_get_thread_num            
      integer ndim, nfun, ndeg, npoly, id
      real *8 z(ndim), f(nfun), rx, ry, pi2
      real *8 rr,reps,dgreen
      real *8 xtarg,ytarg
      real *8 :: pt(2,200), pt0(2)
      character :: polytype, polytype0
      integer ndeg0, npoly0, i
      save pt, ndeg, npoly, pi2, polytype

      id = 1
c$    id = omp_get_thread_num()  +1 
      
      rx = z(1) - pt(1,id)
      ry = z(2) - pt(2,id)
      rr = rx*rx +ry*ry
      reps = 1.0d-30
      if (rr.gt.reps) then
         dgreen = rx/(rr*pi2)
         call legetens_pols_2d(z,ndeg,polytype,f)
         do i = 1,npoly
            f(i) = dgreen * f(i)
         enddo
      else
         do i = 1,npoly
            f(i) = 0.0d0
         enddo
      endif

      return

      entry pbox2d_genfunx_init0(polytype0,ndeg0,npoly0)
      pi2 = 8.0d0*datan(1.0d0)
      ndeg=ndeg0
      npoly=npoly0
      polytype=polytype0
      return

      entry pbox2d_genfunx_init1(pt0)
      id = 1
c$    id = omp_get_thread_num()+1
      pt(1,id)=pt0(1)
      pt(2,id)=pt0(2)
      return
      
      return
      end

      subroutine pbox2d_genfuny(ndim, z, nfun, f)
      implicit none
      integer ndim, nfun, ndeg, npoly, id
c$    integer omp_get_thread_num      
      real *8 z(ndim), f(nfun), rx, ry, pi2
      real *8 rr,reps,dgreen
      real *8 xtarg,ytarg
      real *8 :: pt(2,200), pt0(2)
      character :: polytype, polytype0
      integer ndeg0, npoly0, i
      save pt, ndeg, npoly, pi2, polytype

      id = 1
c$    id = omp_get_thread_num() +1     
      
      rx = z(1) - pt(1,id)
      ry = z(2) - pt(2,id)
      rr = rx*rx +ry*ry
      reps = 1.0d-30
      if (rr.gt.reps) then
         dgreen = ry/(rr*pi2)
         call legetens_pols_2d(z,ndeg,polytype,f)
         do i = 1,npoly
            f(i) = dgreen * f(i)
         enddo
      else
         do i = 1,npoly
            f(i) = 0.0d0
         enddo
      endif

      return

      entry pbox2d_genfuny_init0(polytype0,ndeg0,npoly0)
      pi2 = 8.0d0*datan(1.0d0)
      ndeg=ndeg0
      npoly=npoly0
      polytype=polytype0
      return

      entry pbox2d_genfuny_init1(pt0)
      id = 1
c$    id = omp_get_thread_num()+1
      pt(1,id)=pt0(1)
      pt(2,id)=pt0(2)
      return
      
      
      return
      end

