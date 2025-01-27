c
c     the following are accessory plotting routines for
c     visualizing level restricted trees, their lists,
c     and composite grids on these trees. these depend
c     on pyplot.f
c
c     Author: Travis Askham <askhamwhat@gmail.com>
c
      subroutine lrt2d_plotleaves(iw,zll,blength,levelbox,ichildbox,
     1     icolbox,irowbox,nboxes)
      implicit real *8 (a-h,o-z)
      real *8 :: zll(2), blength
      integer :: levelbox(*), ichildbox(4,*), icolbox(*), irowbox(*)
      real *8, allocatable :: x(:), y(:), w(:), h(:), rot(:)
      character *5, allocatable :: names(:)
      integer, allocatable :: ileaf(:)

      allocate(ileaf(nboxes))
      
      nleaf = 0
      do i = 1,nboxes
         if (ichildbox(1,i) .le. 0) then
            nleaf = nleaf+1
            ileaf(nleaf) = i
         endif
      enddo

      allocate(x(nleaf),y(nleaf),w(nleaf),h(nleaf),rot(nleaf))
      allocate(names(nleaf))

      do i = 1,nleaf
         ii = ileaf(i)
         l = levelbox(ii)
         wid = blength/(2d0**l)
         x(i) = zll(1)+(icolbox(ii)-1)*wid
         y(i) = zll(2)+(irowbox(ii)-1)*wid
         w(i) = wid
         h(i) = wid
         rot(i) = 0
      enddo

      ifid = 1
      ifnames = 0
      itype = 1
      call pyplotrects(iw,x,y,w,h,rot,ileaf,names,nleaf,ifid,ifnames,
     1     itype,'leaf boxes *')
      
      return
      end

      subroutine lrt2d_plotlists(iw,ibox,zll,blength,levelbox,
     1     icolleagbox,ichildbox,icolbox,irowbox,nboxes)
      implicit real *8 (a-h,o-z)
      real *8 :: zll(2), blength
      integer :: levelbox(*), ichildbox(4,*), icolbox(*), irowbox(*)
      integer :: icolleagbox(9,*)
c     local
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
      integer  ibox
      integer  inbig12(3), isbig34(3)
      integer  iebig13(1), iwbig24(1)
      integer  iebig3(1), iwbig2(1)
      integer  iebig1(1), iwbig4(1)
      real *8, allocatable :: x(:), y(:), w(:), h(:), rot(:)
      character *5, allocatable :: names(:)
      character *5 :: cname
      integer, allocatable :: iflagnesw(:,:), id(:)

      allocate(iflagnesw(4,nboxes))

      call lrt2d_mklists(ibox,inall,nnall,iynall,in12,nn12,iy12,
     1    isall,nsall,iysall,is34,ns34,iy34,ieall,neall,
     2    iyeall,ie13,ne13,iy13,iwall,nwall,iywall,
     3    iw24,nw24,iy24,iww2,iwy2,nww2,iww4,iwy4,nww4,
     4    iee1,iey1,nee1,iee3,iey3,nee3,inbig12,isbig34,iebig13,
     6    iwbig24,iebig1,iwbig2,iebig3,iwbig4,icolleagbox,ichildbox,
     8     icolbox,irowbox,iflagnesw)


      nmax = 200
      allocate(names(nmax),x(nmax),y(nmax),w(nmax),h(nmax),rot(nmax),
     1     id(nmax))
      
      nfound = 0

      do i = 1,4
         nfound = nfound+1
         ii = ichildbox(i,ibox)
         l = levelbox(ii)
         wid = blength/(2d0**l)
         x(nfound) = zll(1)+(icolbox(ii)-1)*wid
         y(nfound) = zll(2)+(irowbox(ii)-1)*wid
         w(nfound) = wid
         h(nfound) = wid
         rot(nfound) = 0
         write(cname,'(A4,I1)') 'chld', i
         names(nfound) = cname
         id(nfound) = ii
      enddo         
      
      do i = 1,nnall
         nfound = nfound+1
         ii = inall(i)
         l = levelbox(ii)
         wid = blength/(2d0**l)
         x(nfound) = zll(1)+(icolbox(ii)-1)*wid
         y(nfound) = zll(2)+(irowbox(ii)-1)*wid
         w(nfound) = wid
         h(nfound) = wid
         rot(nfound) = 0
         names(nfound) = 'north'
         id(nfound) = ii
      enddo

      do i = 1,nsall
         nfound = nfound+1
         ii = isall(i)
         l = levelbox(ii)
         wid = blength/(2d0**l)
         x(nfound) = zll(1)+(icolbox(ii)-1)*wid
         y(nfound) = zll(2)+(irowbox(ii)-1)*wid
         w(nfound) = wid
         h(nfound) = wid
         rot(nfound) = 0
         names(nfound) = 'south'
         id(nfound) = ii
      enddo

      do i = 1,neall
         nfound = nfound+1
         ii = ieall(i)
         l = levelbox(ii)
         wid = blength/(2d0**l)
         x(nfound) = zll(1)+(icolbox(ii)-1)*wid
         y(nfound) = zll(2)+(irowbox(ii)-1)*wid
         w(nfound) = wid
         h(nfound) = wid
         rot(nfound) = 0
         names(nfound) = 'east '
         id(nfound) = ii
      enddo

      do i = 1,nwall
         nfound = nfound+1
         ii = iwall(i)
         l = levelbox(ii)
         wid = blength/(2d0**l)
         x(nfound) = zll(1)+(icolbox(ii)-1)*wid
         y(nfound) = zll(2)+(irowbox(ii)-1)*wid
         w(nfound) = wid
         h(nfound) = wid
         rot(nfound) = 0
         names(nfound) = 'west '
         id(nfound) = ii
      enddo

      do i = 1,nn12
         nfound = nfound+1
         ii = in12(i)
         l = levelbox(ii)
         wid = blength/(2d0**l)
         x(nfound) = zll(1)+(icolbox(ii)-1)*wid
         y(nfound) = zll(2)+(irowbox(ii)-1)*wid
         w(nfound) = wid
         h(nfound) = wid
         rot(nfound) = 0
         names(nfound) = 'n12  '
         id(nfound) = ii
      enddo

      do i = 1,ns34
         nfound = nfound+1
         ii = is34(i)
         l = levelbox(ii)
         wid = blength/(2d0**l)
         x(nfound) = zll(1)+(icolbox(ii)-1)*wid
         y(nfound) = zll(2)+(irowbox(ii)-1)*wid
         w(nfound) = wid
         h(nfound) = wid
         rot(nfound) = 0
         names(nfound) = 's34  '
         id(nfound) = ii
      enddo

      do i = 1,ne13
         nfound = nfound+1
         ii = ie13(i)
         l = levelbox(ii)
         wid = blength/(2d0**l)
         x(nfound) = zll(1)+(icolbox(ii)-1)*wid
         y(nfound) = zll(2)+(irowbox(ii)-1)*wid
         w(nfound) = wid
         h(nfound) = wid
         rot(nfound) = 0
         names(nfound) = 'e13  '
         id(nfound) = ii
      enddo


      do i = 1,nw24
         nfound = nfound+1
         ii = iw24(i)
         l = levelbox(ii)
         wid = blength/(2d0**l)
         x(nfound) = zll(1)+(icolbox(ii)-1)*wid
         y(nfound) = zll(2)+(irowbox(ii)-1)*wid
         w(nfound) = wid
         h(nfound) = wid
         rot(nfound) = 0
         names(nfound) = 'w24  '
         id(nfound) = ii
      enddo

      do i = 1,nee1
         nfound = nfound+1
         ii = iee1(i)
         l = levelbox(ii)
         wid = blength/(2d0**l)
         x(nfound) = zll(1)+(icolbox(ii)-1)*wid
         y(nfound) = zll(2)+(irowbox(ii)-1)*wid
         w(nfound) = wid
         h(nfound) = wid
         rot(nfound) = 0
         names(nfound) = 'ee1  '
         id(nfound) = ii
      enddo

      do i = 1,nee3
         nfound = nfound+1
         ii = iee3(i)
         l = levelbox(ii)
         wid = blength/(2d0**l)
         x(nfound) = zll(1)+(icolbox(ii)-1)*wid
         y(nfound) = zll(2)+(irowbox(ii)-1)*wid
         w(nfound) = wid
         h(nfound) = wid
         rot(nfound) = 0
         names(nfound) = 'ee3  '
         id(nfound) = ii
      enddo

      do i = 1,nww2
         nfound = nfound+1
         ii = iww2(i)
         l = levelbox(ii)
         wid = blength/(2d0**l)
         x(nfound) = zll(1)+(icolbox(ii)-1)*wid
         y(nfound) = zll(2)+(irowbox(ii)-1)*wid
         w(nfound) = wid
         h(nfound) = wid
         rot(nfound) = 0
         names(nfound) = 'ww2  '
         id(nfound) = ii
      enddo

      do i = 1,nww4
         nfound = nfound+1
         ii = iww4(i)
         l = levelbox(ii)
         wid = blength/(2d0**l)
         x(nfound) = zll(1)+(icolbox(ii)-1)*wid
         y(nfound) = zll(2)+(irowbox(ii)-1)*wid
         w(nfound) = wid
         h(nfound) = wid
         rot(nfound) = 0
         names(nfound) = 'ww4  '
         id(nfound) = ii
      enddo

      do i = 1,3
         if (inbig12(i) .gt. 0) then
            nfound = nfound+1
            ii = inbig12(i)
            l = levelbox(ii)
            wid = blength/(2d0**l)
            x(nfound) = zll(1)+(icolbox(ii)-1)*wid
            y(nfound) = zll(2)+(irowbox(ii)-1)*wid
            w(nfound) = wid
            h(nfound) = wid
            rot(nfound) = 0
            names(nfound) = 'nbg12'
            id(nfound) = ii
         endif
      enddo

      do i = 1,3
         if (isbig34(i) .gt. 0) then
            nfound = nfound+1
            ii = isbig34(i)
            l = levelbox(ii)
            wid = blength/(2d0**l)
            x(nfound) = zll(1)+(icolbox(ii)-1)*wid
            y(nfound) = zll(2)+(irowbox(ii)-1)*wid
            w(nfound) = wid
            h(nfound) = wid
            rot(nfound) = 0
            names(nfound) = 'sbg34'
            id(nfound) = ii
         endif
      enddo

      if (iebig13(1) .gt. 0) then
         nfound = nfound+1
         ii = iebig13(1)
         l = levelbox(ii)
         wid = blength/(2d0**l)
         x(nfound) = zll(1)+(icolbox(ii)-1)*wid
         y(nfound) = zll(2)+(irowbox(ii)-1)*wid
         w(nfound) = wid
         h(nfound) = wid
         rot(nfound) = 0
         names(nfound) = 'ebg13'
         id(nfound) = ii
      endif

      if (iwbig24(1) .gt. 0) then
         nfound = nfound+1
         ii = iwbig24(1)
         l = levelbox(ii)
         wid = blength/(2d0**l)
         x(nfound) = zll(1)+(icolbox(ii)-1)*wid
         y(nfound) = zll(2)+(irowbox(ii)-1)*wid
         w(nfound) = wid
         h(nfound) = wid
         rot(nfound) = 0
         names(nfound) = 'wbg24'
         id(nfound) = ii
      endif

      if (iwbig2(1) .gt. 0) then
         nfound = nfound+1
         ii = iwbig2(1)
         l = levelbox(ii)
         wid = blength/(2d0**l)
         x(nfound) = zll(1)+(icolbox(ii)-1)*wid
         y(nfound) = zll(2)+(irowbox(ii)-1)*wid
         w(nfound) = wid
         h(nfound) = wid
         rot(nfound) = 0
         names(nfound) = 'wbig2'
         id(nfound) = ii
      endif

      if (iwbig4(1) .gt. 0) then
         nfound = nfound+1
         ii = iwbig4(1)
         l = levelbox(ii)
         wid = blength/(2d0**l)
         x(nfound) = zll(1)+(icolbox(ii)-1)*wid
         y(nfound) = zll(2)+(irowbox(ii)-1)*wid
         w(nfound) = wid
         h(nfound) = wid
         rot(nfound) = 0
         names(nfound) = 'wbig4'
         id(nfound) = ii
      endif

      if (iebig1(1) .gt. 0) then
         nfound = nfound+1
         ii = iebig1(1)
         l = levelbox(ii)
         wid = blength/(2d0**l)
         x(nfound) = zll(1)+(icolbox(ii)-1)*wid
         y(nfound) = zll(2)+(irowbox(ii)-1)*wid
         w(nfound) = wid
         h(nfound) = wid
         rot(nfound) = 0
         names(nfound) = 'ebig1'
         id(nfound) = ii
      endif

      if (iebig3(1) .gt. 0) then
         nfound = nfound+1
         ii = iebig3(1)
         l = levelbox(ii)
         wid = blength/(2d0**l)
         x(nfound) = zll(1)+(icolbox(ii)-1)*wid
         y(nfound) = zll(2)+(irowbox(ii)-1)*wid
         w(nfound) = wid
         h(nfound) = wid
         rot(nfound) = 0
         names(nfound) = 'ebig3'
         id(nfound) = ii
      endif

      ifid = 0
      ifnames = 1
      itype = 1
      call pyplotrects(iw,x,y,w,h,rot,id,names,nfound,ifid,ifnames,
     1     itype,'leaf boxes *')
      

      return
      end
      subroutine compositegrid2d_plotleafpoints(iw,pts,npts,
     1     levelbox,icolbox,irowbox,ichildbox,nboxes,zll,blength)
      implicit real *8 (a-h,o-z)
      real *8 :: pts(2,npts,*), zll(2)
      integer :: levelbox(*), icolbox(*), irowbox(*)
      integer :: ichildbox(4,*)
c     local
      real *8, allocatable :: x2(:), y2(:)
      real *8, allocatable :: x(:), y(:), w(:), h(:), rot(:)
      character *5 :: names
      integer :: id
      integer, allocatable :: ileaf(:)
      character *30 :: title

      allocate(ileaf(nboxes))
      
      nleaf = 0
      do i = 1,nboxes
         if (ichildbox(1,i) .le. 0) then
            nleaf = nleaf+1
            ileaf(nleaf) = i
         endif
      enddo
      
      allocate(x(nleaf),y(nleaf),w(nleaf),h(nleaf),rot(nleaf))
      do i = 1,nleaf
         ii = ileaf(i)
         l = levelbox(ii)
         wid = blength/(2d0**l)
         x(i) = zll(1)+(icolbox(ii)-1)*wid
         y(i) = zll(2)+(irowbox(ii)-1)*wid
         w(i) = wid
         h(i) = wid
         rot(i) = 0
      enddo

      n2 = nleaf*npts
      allocate(x2(n2),y2(n2))
      
      do i = 1,nleaf
         ii = ileaf(i)
         istart = (i-1)*npts
         do j = 1,npts
            jj = istart + j
            x2(jj) = pts(1,j,ii)
            y2(jj) = pts(2,j,ii)
         enddo
      enddo

      ifid=0
      ifnames=0
      title = 'leaf boxes and leaf nodes *'
      call pyplotrectsandpts(iw,x,y,w,h,rot,id,names,nleaf,
     1     ifid,ifnames,x2,y2,n2,itype,title)
      
      
      return
      end

      


c
c
      subroutine pyplotrects(iw,x,y,w,h,rot,id,names,n,ifid,ifnames,
     1     itype,title)
        implicit real *8 (a-h,o-z)
        real *8 x(*),y(*),w(*),h(*),rot(*)
        integer id(*)
        character *5 names(*)
        character *1 a1,a10,file1(11),file11(9),title(1),file1p(10),
     1      temp(32), file19(11), file110(11)
        character *2 ls1
        character *9 file4
        character *11 file8, file9, file10
        character *32 title2
        character *10 file8p
c
        equivalence (file1,file8), (file11,file4), (temp,title2),
     1        (file1p,file8p), (file19,file9), (file110,file10)

c
c       plots the points determined by x,y using itype style
c
        i1=mod(iw,10)
        i10=(iw-i1)/10
        call int2char2(i1,a1)
        call int2char2(i10,a10)

        file8='plotiw.dat1'
        file1(5)=a10
        file1(6)=a1
        file9='plotiw.dat2'
        file19(5)=a10
        file19(6)=a1
        file10='plotiw.dat3'
        file110(5)=a10
        file110(6)=a1
c
        file4='plotiw.py'
        file11(5)=a10
        file11(6)=a1

        file8p='plotiw.pdf'
        file1p(5)=a10
        file1p(6)=a1
c
c
c       print out the contents of the scripting file, plotiw.py
c
        iun87=877
        open(unit=iun87,file=file4)
c
        call quamesslen3(title,nchar)
c
        title2=''
        do 1200 i=1,nchar
        temp(i)=title(i)
 1200 continue
c
        ls1='k.'
        if (itype .eq. 2) ls1='kx'
        if (itype .eq. 3) ls1='k-'
c
        write(iun87,'(a)') '#!/usr/bin/python'
        write(iun87,'(a)') 'import matplotlib.pyplot as pt'
        write(iun87,'(a)') 'import numpy as np'
        write(iun87,'(a)') 'from matplotlib.patches import Rectangle'
        write(iun87,'(a)')
     1       'from matplotlib.collections import PatchCollection'
        
        write(iun87,'(a)') ''
c
        write(iun87,'(a,a,a)') 'x = np.loadtxt("',file8,'")'
        write(iun87,'(a)') 'x = x.reshape(int(np.size(x)/5),5)'
        write(iun87,'(a)')
     1   'rcts = [Rectangle([x_,y_],w_,h_,r_) for x_,y_,w_,h_,r_ in x]'

        write(iun87,'(a)')
     1 'clctn = PatchCollection(rcts,facecolors="None",edgecolors="k")'

        write(iun87,'(a)') 'ax = pt.gca()'
        write(iun87,'(a)') 'ax.add_collection(clctn)'
        
        if (ifid .eq. 1) then
           write(iun87,'(a,a,a)') 'id = np.loadtxt("',file9,'")'
           write(iun87,'(a)') 'for i in range(len(rcts)):'
           write(iun87,'(a)') '   R = rcts[i]'
           write(iun87,'(a)') '   rx, ry = R.get_xy()'
           write(iun87,'(a)') '   cx = rx + R.get_width()/2'
           write(iun87,'(a)') '   cy = ry + R.get_height()/2'
           write(iun87,'(a)')
     1  '   ax.annotate(str(int(id[i])), (cx, cy), color="k",'
           write(iun87,'(a)')
     1  '    weight="bold",fontsize=6, ha="center", va="center")'
           write(iun87,'(a)') ''
           
        endif   
        
        if (ifnames .eq. 1) then
           write(iun87,'(a,a,a)') 'names = np.loadtxt("',file10,'",'
           write(iun87,'(a)') '      dtype="string")'
           write(iun87,'(a)') 'for i in range(len(rcts)):'
           write(iun87,'(a)') '   R = rcts[i]'
           write(iun87,'(a)') '   rx, ry = R.get_xy()'
           write(iun87,'(a)') '   cx = rx + R.get_width()/2'
           write(iun87,'(a)') '   cy = ry + R.get_height()/2'
           write(iun87,'(a)')
     1  '   ax.annotate(str(names[i]), (cx, cy), color="k",'
           write(iun87,'(a)')
     1  '    weight="bold",fontsize=6, ha="center", va="center")'
           write(iun87,'(a)') ''
           
        endif   
        
c
        write(iun87,'(a,a,a)') 'pt.title("', title2, '")'
        write(iun87, '(a)') 'pt.gca().set_aspect("equal")'
ccc        write(iun87,'(a)') 'pt.show()'
        write(iun87,'(a,a,a)') 'pt.savefig("',file8p,'")'
        close(iun87)

c
c       now print out the data file, plotiw.dat
c
        iun88=888
        open(unit=iun88,file=file8)
        do i=1,n
          write(iun88,*) x(i), y(i), w(i), h(i), rot(i)
        enddo
        close(iun88)

        if (ifid .eq. 1) then
           iun88=888
           open(unit=iun88,file=file9)
           do i=1,n
              write(iun88,*) id(i)
           enddo
           close(iun88)
        endif

        if (ifnames .eq. 1) then
           iun88=888
           open(unit=iun88,file=file10)
           do i=1,n
              write(iun88,*) names(i)
           enddo
           close(iun88)
        endif

c
        return
        end
c
c
      subroutine pyplotrectsandpts(iw,x,y,w,h,rot,id,names,n,
     1     ifid,ifnames,x2,y2,n2,itype,title)
        implicit real *8 (a-h,o-z)
        real *8 x(*),y(*),w(*),h(*),rot(*),x2(*),y2(*)
        integer id(*)
        character *5 names(*)
        character *1 a1,a10,file18(11),file14(9),title(1),file18p(10),
     1      temp(32), file19(11), file110(11), file111(11)
        character *2 ls1
        character *9 file4
        character *11 file8, file9, file10, file11
        character *32 title2
        character *10 file8p
c
        equivalence (file18,file8), (file14,file4), (temp,title2),
     1       (file18p,file8p), (file19,file9), (file110,file10),
     2       (file111,file11)

c
c       plots the points determined by x,y using itype style
c
        i1=mod(iw,10)
        i10=(iw-i1)/10
        call int2char2(i1,a1)
        call int2char2(i10,a10)

        file8='plotiw.dat1'
        file18(5)=a10
        file18(6)=a1
        file9='plotiw.dat2'
        file19(5)=a10
        file19(6)=a1
        file10='plotiw.dat3'
        file110(5)=a10
        file110(6)=a1
        file11='plotiw.dat4'
        file111(5)=a10
        file111(6)=a1
c
        file4='plotiw.py'
        file14(5)=a10
        file14(6)=a1

        file8p='plotiw.pdf'
        file18p(5)=a10
        file18p(6)=a1
c
c
c       print out the contents of the scripting file, plotiw.py
c
        iun87=877
        open(unit=iun87,file=file4)
c
        call quamesslen3(title,nchar)
c
        title2=''
        do 1200 i=1,nchar
        temp(i)=title(i)
 1200 continue
c
        ls1='k.'
        if (itype .eq. 2) ls1='kx'
        if (itype .eq. 3) ls1='k-'
c
        write(iun87,'(a)') '#!/usr/bin/python'
        write(iun87,'(a)') 'import matplotlib.pyplot as pt'
        write(iun87,'(a)') 'import numpy as np'
        write(iun87,'(a)') 'from matplotlib.patches import Rectangle'
        write(iun87,'(a)')
     1       'from matplotlib.collections import PatchCollection'
        
        write(iun87,'(a)') ''
c
        write(iun87,'(a,a,a)') 'x = np.loadtxt("',file8,'")'
        write(iun87,'(a)') 'x = x.reshape(int(np.size(x)/5),5)'
        write(iun87,'(a)')
     1   'rcts = [Rectangle([x_,y_],w_,h_,r_) for x_,y_,w_,h_,r_ in x]'

        write(iun87,'(a)')
     1 'clctn = PatchCollection(rcts,facecolors="None",edgecolors="k")'

        write(iun87,'(a)') 'ax = pt.gca()'
        write(iun87,'(a)') 'ax.add_collection(clctn)'
        
        if (ifid .eq. 1) then
           write(iun87,'(a,a,a)') 'id = np.loadtxt("',file9,'")'
           write(iun87,'(a)') 'for i in range(len(rcts)):'
           write(iun87,'(a)') '   R = rcts[i]'
           write(iun87,'(a)') '   rx, ry = R.get_xy()'
           write(iun87,'(a)') '   cx = rx + R.get_width()/2'
           write(iun87,'(a)') '   cy = ry + R.get_height()/2'
           write(iun87,'(a)')
     1  '   ax.annotate(str(int(id[i])), (cx, cy), color="k",'
           write(iun87,'(a)')
     1  '    weight="bold",fontsize=6, ha="center", va="center")'
           write(iun87,'(a)') ''
           
        endif   
        
        if (ifnames .eq. 1) then
           write(iun87,'(a,a,a)') 'names = np.loadtxt("',file10,'",'
           write(iun87,'(a)') '      dtype="string")'
           write(iun87,'(a)') 'for i in range(len(rcts)):'
           write(iun87,'(a)') '   R = rcts[i]'
           write(iun87,'(a)') '   rx, ry = R.get_xy()'
           write(iun87,'(a)') '   cx = rx + R.get_width()/2'
           write(iun87,'(a)') '   cy = ry + R.get_height()/2'
           write(iun87,'(a)')
     1  '   ax.annotate(str(names[i]), (cx, cy), color="k",'
           write(iun87,'(a)')
     1  '    weight="bold",fontsize=6, ha="center", va="center")'
           write(iun87,'(a)') ''
           
        endif

        write(iun87,'(a,a,a)') 'x2 = np.loadtxt("',file11,'")'
        write(iun87,'(a)') 'x2 = x2.reshape(np.size(x2)/2,2)'
        write(iun87,'(a)') 'pt.plot(x2[:,0],x2[:,1],"kx")'
        
c
        write(iun87,'(a,a,a)') 'pt.title("', title2, '")'
        write(iun87, '(a)') 'pt.gca().set_aspect("equal")'
ccc        write(iun87,'(a)') 'pt.show()'
        write(iun87,'(a,a,a)') 'pt.savefig("',file8p,'")'
        close(iun87)

c
c       now print out the data file, plotiw.dat
c
        iun88=888
        open(unit=iun88,file=file8)
        do i=1,n
           write(iun88,*) x(i), y(i), w(i), h(i), rot(i)
        enddo
        close(iun88)

        if (ifid .eq. 1) then
           iun88=888
           open(unit=iun88,file=file9)
           do i=1,n
              write(iun88,*) id(i)
           enddo
           close(iun88)
        endif

        if (ifnames .eq. 1) then
           iun88=888
           open(unit=iun88,file=file10)
           do i=1,n
              write(iun88,*) names(i)
           enddo
           close(iun88)
        endif

        iun88=888
        open(unit=iun88,file=file11)
        do i=1,n2
          write(iun88,*) x2(i), y2(i)
        enddo
        close(iun88)
        
c
        return
        end
c
c
c
      
      
