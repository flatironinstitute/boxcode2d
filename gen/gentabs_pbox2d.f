cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     author: Travis Askham
c
c     this file is being released under an Apache License (version 2.0
c     See LICENSE in home directory for details 
c
c     this program runs the table generator and prints
c     the table to file

      
      program gentabs_pbox2d
c
c
c     unit tests for level-restricted tree
c     
c      
      implicit real *8 (a-h,o-z)

      
      real *8, allocatable, dimension(:,:,:) :: tabcol,tabbtos,tabstob
      real *8, allocatable, dimension(:,:,:,:) ::gradtabcol,gradtabbtos,
     1     gradtabstob
      character :: pttype, polytype
      character *40 :: filename, subname
      character *1 :: fname(40), sname(40), ci1, ci2
      equivalence (fname,filename), (sname,subname)
      
      call prini(6,13)


      write(*,*) "MENU -------------------"
      write(*,*) "norder options: any positive integer (rec: norder<=32"
      write(*,*) "point type options (only symmetric sets here): "
c      write(*,*) "- 'T', symmetric node set approx size of 'T' poly set"
      write(*,*) "- 'F', full tensor product grid of legendre nodes"
      write(*,*) "polynomial type options: "
      write(*,*) "- 'T', total degree polys degree norder-1"
      write(*,*) "- 'F', full tensor product polys degree norder-1"

      write(*,*) 'enter norder: '
      read(*,*) norder
      write(*,*) 'enter point type: '
      read(*,*) pttype
      write(*,*) 'enter polynomial type: '
      read(*,*) polytype

      if (pttype .eq. 'F') then
         npt = (norder)**2
      else
         write(*,*) 'selected point type unavailable. no table produced'
         stop
      endif

      ndeg = norder-1
      
      if (polytype .eq. 'F' .or. polytype .eq. 'T') then
         call legetens_npol_2d(ndeg,polytype,npoly)
      else
         write(*,*) 'selected poly type unavailable. no table produced'
         stop
      endif

      i1 = norder/10
      i2 = norder-10*(norder/10)

      write(ci1,'(i1)') i1
      write(ci2,'(i1)') i2
      
      filename = 'pbox2dtab_ptX_polyX_norderXX.f'
      fname(13) = pttype
      fname(19) = polytype
      fname(27) = ci1
      fname(28) = ci2

      subname = 'pbox2dtab_ptX_polyX_norderXX_getref'
      sname(13) = pttype
      sname(19) = polytype
      sname(27) = ci1
      sname(28) = ci2

      write(*,*) 'the reference table subroutine named: '
      write(*,*) '     ',  trim(subname)
      write(*,*) 'will be written to: '
      write(*,*) '     ', filename

      allocate(tabcol(npt,npoly,3),tabstob(npt,npoly,2),
     1     tabbtos(npt,npoly,2))
      allocate(gradtabcol(2,npt,npoly,3),gradtabstob(2,npt,npoly,2),
     1     gradtabbtos(2,npt,npoly,2))


      write(*,*) 'npt ', npt
      write(*,*) 'npoly ', npoly

      maxcls = max(10000000,(10000000*norder/16))

      call cpu_time(t0)
c$    t0 = omp_get_wtime()      
      call pbox2d_genreftab(pttype,polytype,norder,npoly,npt,
     1     tabcol,tabstob,tabbtos,gradtabcol,gradtabstob,gradtabbtos,
     2     maxcls,ier)
      call cpu_time(t1)
c$    t1 = omp_get_wtime()      

      write(*,*) 'after generating tables, ier = ', ier
      write(*,*) 'time for computation ', t1-t0
      write(*,*) 'writing tables to file ...'
      
      call gentabs_writedreftabs(filename,subname,npt,npoly,
     1     tabcol,tabstob,tabbtos,gradtabcol,gradtabstob,gradtabbtos)

      stop
      end
      
      subroutine gentabs_writedreftabs(filename,subname,npt,npoly,
     1     tabcol,tabstob,tabbtos,gradtabcol,gradtabstob,gradtabbtos)
      implicit real *8 (a-h,o-z)
      character *40 :: filename,subname
      real *8 :: tabcol(npt,npoly,3), tabstob(npt,npoly,2),
     1     tabbtos(npt,npoly,2)
      real *8 :: gradtabcol(2,npt,npoly,3),
     1     gradtabstob(2,npt,npoly,2),gradtabbtos(2,npt,npoly,2)
c     local
      character *30 :: colnamelist(2), namelist(4), vname
      
      
      iw = 888
      open(unit=iw, file=filename)

      write(iw,*) '     subroutine ', trim(subname), '(npt,npoly,'
      write(iw,*) '    1 tabcol,tabstob,tabbtos,'
      write(iw,*) '    2 gradtabcol,gradtabstob,gradtabbtos)'
      write(iw,*) '     implicit real *8 (a-h,o-z)'
 2002 format('      dimension ',a,'(',i4,',',i4,',*)')
 2003 format('      dimension ',a,'(',i1,',',i4,',',i4,',*)')

      write(iw,2002) 'tabcol', npt, npoly
      write(iw,2002) 'tabstob', npt, npoly
      write(iw,2002) 'tabbtos', npt, npoly
      write(iw,2003) 'gradtabcol', 2, npt, npoly
      write(iw,2003) 'gradtabstob', 2, npt, npoly
      write(iw,2003) 'gradtabbtos', 2, npt, npoly 

      write(iw,'(a,i4)') '      npoly = ',  npoly
      write(iw,'(a,i4)') '      npt = ',  npt
      
 1003 format('      ',a,'(',i4,',',i4,',',i4,') = ',d25.18)
 1004 format('      ',a,'(',i4,',',i4,',',i4,',',i4,') = ',d25.18)

      colnamelist(1) = 'tabcol'
      colnamelist(2) = 'gradtabcol'
      do l = 1,2
         vname = colnamelist(l)
         do jj = 1,3
            j2 = jj
            do i = 1,npoly
            do j = 1,npt
               if (l .eq. 1) then
                  write(iw,1003) trim(vname), j, i, j2, tabcol(j,i,jj)
               endif
               if (l .eq. 2) then
                  write(iw,1004) trim(vname), 1, j, i, j2,
     1                 gradtabcol(1,j,i,jj)
                  write(iw,1004) trim(vname), 2, j, i, j2,
     1                 gradtabcol(2,j,i,jj)
               endif
            enddo
            enddo
         enddo
      enddo

      namelist(1) = 'tabstob'
      namelist(2) = 'gradtabstob'
      namelist(3) = 'tabbtos'
      namelist(4) = 'gradtabbtos'
      do l = 1,4
         vname = namelist(l)
         do jj = 1,2
            j2 = jj
            do i = 1,npoly
            do j = 1,npt
               if (l .eq. 1) then
                  write(iw,1003) trim(vname), j, i, j2, tabstob(j,i,jj)
               endif
               if (l .eq. 2) then
                  write(iw,1004) trim(vname), 1, j, i, j2,
     1                 gradtabstob(1,j,i,jj)
                  write(iw,1004) trim(vname), 2, j, i, j2,
     1                 gradtabstob(2,j,i,jj)
               endif
               if (l .eq. 3) then
                  write(iw,1003) trim(vname), j, i, j2, tabbtos(j,i,jj)
               endif
               if (l .eq. 4) then
                  write(iw,1004) trim(vname), 1, j, i, j2,
     1                 gradtabbtos(1,j,i,jj)
                  write(iw,1004) trim(vname), 2, j, i, j2,
     1                 gradtabbtos(2,j,i,jj)
               endif
            enddo
            enddo
         enddo
      enddo

      write(iw,*) '      return '
      write(iw,*) '      end '

      close(iw)

      return
      end
