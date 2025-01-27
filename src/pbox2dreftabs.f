
      subroutine pbox2dreftabs(norder,pttype,polytype,npt,npoly,
     1     tabcol,tabstob,tabbtos,gradtabcol,gradtabstob,gradtabbtos,
     2     ier)
c
c     load the requested reference table if it's available
c
c     inputs
c     
c     pttype - character *1, point type
c       pttype .eq. 'F' -> full tensor product grid of legendre pts of
c                          order norder
c       pttype .eq. 'T' -> symmetric grid of cardinality
c                          approximately equal to norder*(norder+1)/2
c     polytype - character *1, polynomial type 
c       polytype .eq. 'F' -> full tensor product set of lege polynomials
c       polytype .eq. 'T' -> lege products of total degree <= norder-1
c     npt - integer. number of points in pt set, first dimension
c              of tabcol, tabstob, etc
c     npoly - integer. number of polynomials, second dimension
c              of tabcol, tabstob, etc
c
c     outputs
c
c     tabcols, tabstob, etc. -> the requested tables
c     ier - integer. error flag
c       ier .eq. 0 -> normal execution
c       ier .eq. 1 -> npts or npoly inconsistent with settings for
c                     norder, pttype, polytype
c       ier .eq. 2 -> no rule found for requested settings 
c

      implicit real *8 (a-h,o-z)
      character :: pttype, polytype
      integer :: norder, npoly, npt
      real *8 :: tabcol(npt,npoly,*), tabstob(npt,npoly,*),
     1     tabbtos(npt,npoly,*)
      real *8 :: gradtabcol(2,npt,npoly,*), gradtabstob(2,npt,npoly,*),
     1     gradtabbtos(2,npt,npoly,*)

      ier = 0

      if (norder .eq. 2 .and. pttype .eq. 'F' .and. polytype .eq. 'T')
     1     then
         call pbox2dtab_ptF_polyT_norder02_getref(npt,npoly,
     1        tabcol,tabstob,tabbtos,
     2        gradtabcol,gradtabstob,gradtabbtos)
         return
      endif

      if (norder .eq. 2 .and. pttype .eq. 'F' .and. polytype .eq. 'F')
     1     then
         call pbox2dtab_ptF_polyF_norder02_getref(npt,npoly,
     1        tabcol,tabstob,tabbtos,
     2        gradtabcol,gradtabstob,gradtabbtos)
         return
      endif

      if (norder .eq. 3 .and. pttype .eq. 'F' .and. polytype .eq. 'T')
     1     then
         call pbox2dtab_ptF_polyT_norder03_getref(npt,npoly,
     1        tabcol,tabstob,tabbtos,
     2        gradtabcol,gradtabstob,gradtabbtos)
         return
      endif

      if (norder .eq. 3 .and. pttype .eq. 'F' .and. polytype .eq. 'F')
     1     then
         call pbox2dtab_ptF_polyF_norder03_getref(npt,npoly,
     1        tabcol,tabstob,tabbtos,
     2        gradtabcol,gradtabstob,gradtabbtos)
         return
      endif

      if (norder .eq. 4 .and. pttype .eq. 'F' .and. polytype .eq. 'T')
     1     then
         call pbox2dtab_ptF_polyT_norder04_getref(npt,npoly,
     1        tabcol,tabstob,tabbtos,
     2        gradtabcol,gradtabstob,gradtabbtos)
         return
      endif

      if (norder .eq. 4 .and. pttype .eq. 'F' .and. polytype .eq. 'F')
     1     then
         call pbox2dtab_ptF_polyF_norder04_getref(npt,npoly,
     1        tabcol,tabstob,tabbtos,
     2        gradtabcol,gradtabstob,gradtabbtos)
         return
      endif

      if (norder .eq. 5 .and. pttype .eq. 'F' .and. polytype .eq. 'T')
     1     then
         call pbox2dtab_ptF_polyT_norder05_getref(npt,npoly,
     1        tabcol,tabstob,tabbtos,
     2        gradtabcol,gradtabstob,gradtabbtos)
         return
      endif

      if (norder .eq. 5 .and. pttype .eq. 'F' .and. polytype .eq. 'F')
     1     then
         call pbox2dtab_ptF_polyF_norder05_getref(npt,npoly,
     1        tabcol,tabstob,tabbtos,
     2        gradtabcol,gradtabstob,gradtabbtos)
         return
      endif

      if (norder .eq. 6 .and. pttype .eq. 'F' .and. polytype .eq. 'T')
     1     then
         call pbox2dtab_ptF_polyT_norder06_getref(npt,npoly,
     1        tabcol,tabstob,tabbtos,
     2        gradtabcol,gradtabstob,gradtabbtos)
         return
      endif

      if (norder .eq. 6 .and. pttype .eq. 'F' .and. polytype .eq. 'F')
     1     then
         call pbox2dtab_ptF_polyF_norder06_getref(npt,npoly,
     1        tabcol,tabstob,tabbtos,
     2        gradtabcol,gradtabstob,gradtabbtos)
         return
      endif

      if (norder .eq. 7 .and. pttype .eq. 'F' .and. polytype .eq. 'T')
     1     then
         call pbox2dtab_ptF_polyT_norder07_getref(npt,npoly,
     1        tabcol,tabstob,tabbtos,
     2        gradtabcol,gradtabstob,gradtabbtos)
         return
      endif

      if (norder .eq. 7 .and. pttype .eq. 'F' .and. polytype .eq. 'F')
     1     then
         call pbox2dtab_ptF_polyF_norder07_getref(npt,npoly,
     1        tabcol,tabstob,tabbtos,
     2        gradtabcol,gradtabstob,gradtabbtos)
         return
      endif

      if (norder .eq. 8 .and. pttype .eq. 'F' .and. polytype .eq. 'T')
     1     then
         call pbox2dtab_ptF_polyT_norder08_getref(npt,npoly,
     1        tabcol,tabstob,tabbtos,
     2        gradtabcol,gradtabstob,gradtabbtos)
         return
      endif

      if (norder .eq. 8 .and. pttype .eq. 'F' .and. polytype .eq. 'F')
     1     then
         call pbox2dtab_ptF_polyF_norder08_getref(npt,npoly,
     1        tabcol,tabstob,tabbtos,
     2        gradtabcol,gradtabstob,gradtabbtos)
         return
      endif

      if (norder .eq. 9 .and. pttype .eq. 'F' .and. polytype .eq. 'T')
     1     then
         call pbox2dtab_ptF_polyT_norder09_getref(npt,npoly,
     1        tabcol,tabstob,tabbtos,
     2        gradtabcol,gradtabstob,gradtabbtos)
         return
      endif

      if (norder .eq. 9 .and. pttype .eq. 'F' .and. polytype .eq. 'F')
     1     then
         call pbox2dtab_ptF_polyF_norder09_getref(npt,npoly,
     1        tabcol,tabstob,tabbtos,
     2        gradtabcol,gradtabstob,gradtabbtos)
         return
      endif

      if (norder .eq. 10 .and. pttype .eq. 'F' .and. polytype .eq. 'T')
     1     then
         call pbox2dtab_ptF_polyT_norder10_getref(npt,npoly,
     1        tabcol,tabstob,tabbtos,
     2        gradtabcol,gradtabstob,gradtabbtos)
         return
      endif

      if (norder .eq. 10 .and. pttype .eq. 'F' .and. polytype .eq. 'F')
     1     then
         call pbox2dtab_ptF_polyF_norder10_getref(npt,npoly,
     1        tabcol,tabstob,tabbtos,
     2        gradtabcol,gradtabstob,gradtabbtos)
         return
      endif

      if (norder .eq. 12 .and. pttype .eq. 'F' .and. polytype .eq. 'T')
     1     then
         call pbox2dtab_ptF_polyT_norder12_getref(npt,npoly,
     1        tabcol,tabstob,tabbtos,
     2        gradtabcol,gradtabstob,gradtabbtos)
         return
      endif

      if (norder .eq. 12 .and. pttype .eq. 'F' .and. polytype .eq. 'F')
     1     then
         call pbox2dtab_ptF_polyF_norder12_getref(npt,npoly,
     1        tabcol,tabstob,tabbtos,
     2        gradtabcol,gradtabstob,gradtabbtos)
         return
      endif

      if (norder .eq. 16 .and. pttype .eq. 'F' .and. polytype .eq. 'T')
     1     then
         call pbox2dtab_ptF_polyT_norder16_getref(npt,npoly,
     1        tabcol,tabstob,tabbtos,
     2        gradtabcol,gradtabstob,gradtabbtos)
         return
      endif

      if (norder .eq. 16 .and. pttype .eq. 'F' .and. polytype .eq. 'F')
     1     then
         call pbox2dtab_ptF_polyF_norder16_getref(npt,npoly,
     1        tabcol,tabstob,tabbtos,
     2        gradtabcol,gradtabstob,gradtabbtos)
         return
      endif

      ier = 2
      return
      end
