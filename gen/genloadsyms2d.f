cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     author: Travis Askham
c
c     this file is being released under an Apache License (version 2.0
c     See LICENSE in home directory for details 
c     
c     determine symmetries relating all of the table types to
c     the reference tables (for colleague, stob, and btos)
c
c     this is done by brute force
c
c     subroutines are printed to the file loadsyms2d.f
      
      implicit real *8 (a-h,o-z)

      iw = 20
      nq = 6

      open(unit = iw, file='loadsyms2d.f')
      call prini(6,13)

      call testit(iw,nq)

      stop
      end




      subroutine testit(iw,nq)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      implicit real *8 (a-h,o-z)
      dimension wc(2,nq**2,9), wbtos(2,nq**2,12), wstob(2,nq**2,12)
      dimension wcref(2,nq**2,3), wbtosref(2,nq**2,2)
      dimension wstobref(2,nq**2,2)
      dimension xt(4*nq**2), yt(4*nq**2)
      dimension xts(nq**2), yts(nq**2)
      dimension xq(nq), xyc(2)
      dimension dlens(nq**2), dlensref(nq**2)
      dimension alldist(nq**2,nq**2), alldistref(nq**2,nq**2)
      dimension ipperm(nq**2), icperm(nq**2), icsign(nq**2)
      dimension ippermold(nq**2)
      dimension ibtosdimperm(12), ibtosflipopt(12), ibtosref(12)
      dimension istobdimperm(12), istobflipopt(12), istobref(12)
      dimension icdimperm(12), icflipopt(12), icref(12)
      dimension idimperms(2,2)
      data idimperms / 1,2, 2,1 /
      data ndimperms / 2 /
      dimension iflipopts(2,4)      
      data iflipopts / 1,1, -1,1, 1,-1, -1,-1 /
      
      data nflipopts / 4 /

      real *8, allocatable :: grid0(:,:)

      character *32 title
      character cp
      logical iwork

      cp = 't'
      

      nq2 = nq**2

c     tolerance for a match

      tol = 1.0d-14

c     define box size, box corner for self

      bs = 2.0d0
      xyc(1) = -1
      xyc(2) = -1

c     create tensor chebyshev grid

      allocate(grid0(2,nq2))
      itype = 0
      ldu = 1
      ldv = 1
      call legetens_exps_2d(itype,nq,cp,grid0,u,ldu,v,ldv,w)

c     grab reference points
      
      call tensrefpts2d_grid(grid0,nq2,bs,xyc,wcref,wbtosref,wstobref)

c     grab all targets 

      call alltargs2d_grid(grid0,nq2,bs,xyc,wc,wbtos,wstob)      

      iwp = 1
      itype = 2
c      call zpyplot(iwp,wc,nq2*9,itype,'all colleag*')
      iwp = 2
      itype = 2
c      call zpyplot(iwp,wstob,nq2*12,itype,'all stob*')
      iwp = 3
      itype = 2
c      call zpyplot(iwp,wbtos,nq2*12,itype,'all btos*')

c     find permutations which match distances for btos
      
      do i = 1,12
         iworkever = 0
         call dists2d(wbtos(1,1,i),nq2,wcref(1,1,3),nq2,alldist)
         do j = 1,2
            call dists2d(wbtosref(1,1,j),nq2,wcref(1,1,3),nq2,
     1           alldistref)
            do i1 = 1,ndimperms
               do i2 = 1,nflipopts
c     for given permutation of dimensions and flipping of
c     index orders, check if you've got it
                  call buildperm2d(idimperms(1,i1),
     1                 iflipopts(1,i2),nq,cp,ipperm,icperm,icsign)
c                  call prinf('ipperm *',ipperm,nq**2)
                  iwork = .true.
                  do jj = 1,nq2
                     do ii = 1,nq2
                        iwork = iwork .and.
     1                       (abs(alldist(ipperm(ii),ipperm(jj))
     1                       -alldistref(ii,jj)) .lt. tol)
                     enddo
                  enddo

c     save if new
                  if (iwork .and. iworkever .eq. 0) then
                     ibtosref(i) = j
                     ibtosdimperm(i) = i1
                     ibtosflipopt(i) = i2
                     do iii = 1,nq2
                        ippermold(iii) = ipperm(iii)
                     enddo
                  endif
                  
                  if (iwork) iworkever = iworkever + 1
               enddo
            enddo
         enddo
         if (iworkever .eq. 0) then
            write(*,*) 'fail', i
         else
            write(*,*) 'success', i, ibtosref(i), iworkever
         endif
         
      enddo
      
c     find permutations which match distances for stob
      
      do i = 1,12
         iworkever = 0
         call dists2d(wstob(1,1,i),nq2,wcref(1,1,3),nq2,alldist)
         do j = 1,3
            call dists2d(wstobref(1,1,j),nq2,wcref(1,1,3),nq2,
     1           alldistref)
            do i1 = 1,ndimperms
               do i2 = 1,nflipopts
c     for given permutation of dimensions and flipping of
c     index orders, check if you've got it
                  call buildperm2d(idimperms(1,i1),
     1                 iflipopts(1,i2),nq,cp,ipperm,icperm,icsign)
c                  call prinf('ipperm *',ipperm,nq**2)
                  iwork = .true.
                  do jj = 1,nq2
                     do ii = 1,nq2
                        iwork = iwork .and.
     1                       (abs(alldist(ipperm(ii),ipperm(jj))
     1                       -alldistref(ii,jj)) .lt. tol)
                     enddo
                  enddo

c     save if new
                  if (iwork .and. iworkever .eq. 0) then
                     istobref(i) = j
                     istobdimperm(i) = i1
                     istobflipopt(i) = i2
                     do iii = 1,nq2
                        ippermold(iii) = ipperm(iii)
                     enddo
                  endif
                  
                  if (iwork) iworkever = iworkever + 1
               enddo
            enddo
         enddo
         if (iworkever .eq. 0) then
            write(*,*) 'fail', i
         else
            write(*,*) 'success', i, istobref(i), iworkever
         endif
         
      enddo
      
c     find permutations which match distances for colleague
      
      do i = 1,9
         iworkever = 0
         call dists2d(wc(1,1,i),nq2,wcref(1,1,3),nq2,alldist)
         do j = 1,4
            call dists2d(wcref(1,1,j),nq2,wcref(1,1,3),nq2,
     1           alldistref)
            do i1 = 1,ndimperms
               do i2 = 1,nflipopts
c     for given permutation of dimensions and flipping of
c     index orders, check if you've got it
                  call buildperm2d(idimperms(1,i1),
     1                 iflipopts(1,i2),nq,cp,ipperm,icperm,icsign)
c                  call prinf('ipperm *',ipperm,nq**2)
                  iwork = .true.
                  do jj = 1,nq2
                     do ii = 1,nq2
                        iwork = iwork .and.
     1                       (abs(alldist(ipperm(ii),ipperm(jj))
     1                       -alldistref(ii,jj)) .lt. tol)
                     enddo
                  enddo

c     save if new
                  if (iwork .and. iworkever .eq. 0) then
                     icref(i) = j
                     icdimperm(i) = i1
                     icflipopt(i) = i2
                     do iii = 1,nq2 
                        ippermold(iii) = ipperm(iii)
                     enddo
                  endif
                  
                  if (iwork) iworkever = iworkever + 1
               enddo
            enddo
         enddo
         if (iworkever .eq. 0) then
            write(*,*) 'fail', i
         else
            write(*,*) 'success', i, icref(i), iworkever
         endif
         
      enddo

c     print to file

      write(iw,'(a)') 'c'
      write(iw,'(a)') 'c'            
      write(iw,'(a)') 'c       this file was generated automatically'
      write(iw,'(a)') 'c       it contains subroutines which load'
      write(iw,'(a)') 'c       symmetry info for volume code tables'
      write(iw,'(a)') 'c'      
      write(iw,'(a)') 'c'      
      write(iw,'(a)') ''      
      write(iw,'(a)') ''      
      write(iw,'(a)') ''      
      write(iw,'(a)') ''      
      write(iw,*) '      subroutine loadsyms2dc(iref,idimp,iflip)'
      write(iw,*) '      implicit real *8 (a-h,o-z)'
      write(iw,*) ''
      write(iw,*) '      dimension iref(*), idimp(2,*), iflip(2,*)'
      write(iw,*) ''
      do i = 1,9
         i1 = icdimperm(i)
         i2 = icflipopt(i)         
         write(iw,'(a,I2,a,I2)') '      iref(', i, ')  = ', icref(i)
         write(iw,'(a,I2,a,I2)') '      idimp(1,', i, ') = ',
     1        idimperms(1,i1)
         write(iw,'(a,I2,a,I2)') '      idimp(2,', i, ') = ',
     1        idimperms(2,i1)
         write(iw,'(a,I2,a,I2)') '      iflip(1,', i, ') = ',
     1        iflipopts(1,i2)
         write(iw,'(a,I2,a,I2)') '      iflip(2,', i, ') = ',
     1        iflipopts(2,i2)
      enddo
      write(iw,*) ''
      write(iw,*) '      return'
      write(iw,*) '      end'
      write(iw,*) ''      
      write(iw,*) ''
      
      write(iw,*) '      subroutine loadsyms2dbtos(iref,idimp,iflip)'
      write(iw,*) '      implicit real *8 (a-h,o-z)'
      write(iw,*) ''
      write(iw,*) '      dimension iref(*), idimp(2,*), iflip(2,*)'
      write(iw,*) ''
      do i = 1,12
         i1 = ibtosdimperm(i)
         i2 = ibtosflipopt(i)         
         write(iw,'(a,I2,a,I2)') '      iref(', i, ')  = ', ibtosref(i)
         write(iw,'(a,I2,a,I2)') '      idimp(1,', i, ') = ',
     1        idimperms(1,i1)
         write(iw,'(a,I2,a,I2)') '      idimp(2,', i, ') = ',
     1        idimperms(2,i1)
         write(iw,'(a,I2,a,I2)') '      iflip(1,', i, ') = ',
     1        iflipopts(1,i2)
         write(iw,'(a,I2,a,I2)') '      iflip(2,', i, ') = ',
     1        iflipopts(2,i2)
      enddo
      write(iw,*) ''
      write(iw,*) '      return'
      write(iw,*) '      end'
      write(iw,*) ''      
      write(iw,*) ''
      
      write(iw,*) '      subroutine loadsyms2dstob(iref,idimp,iflip)'
      write(iw,*) '      implicit real *8 (a-h,o-z)'
      write(iw,*) ''
      write(iw,*) '      dimension iref(*), idimp(2,*), iflip(2,*)'
      write(iw,*) ''
      do i = 1,12
         i1 = istobdimperm(i)
         i2 = istobflipopt(i)         
         write(iw,'(a,I2,a,I2)') '      iref(', i, ')  = ', istobref(i)
         write(iw,'(a,I2,a,I2)') '      idimp(1,', i, ') = ',
     1        idimperms(1,i1)
         write(iw,'(a,I2,a,I2)') '      idimp(2,', i, ') = ',
     1        idimperms(2,i1)
         write(iw,'(a,I2,a,I2)') '      iflip(1,', i, ') = ',
     1        iflipopts(1,i2)
         write(iw,'(a,I2,a,I2)') '      iflip(2,', i, ') = ',
     1        iflipopts(2,i2)
      enddo
      write(iw,*) ''
      write(iw,*) '      return'
      write(iw,*) '      end'
      write(iw,*) ''      
      write(iw,*) ''
      
      return
      end

