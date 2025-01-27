c
c
c       this file was generated automatically
c       it contains subroutines which load
c       symmetry info for volume code tables
c
c




       subroutine loadsyms2dc(iref,idimp,iflip)
       implicit real *8 (a-h,o-z)
 
       dimension iref(*), idimp(2,*), iflip(2,*)
 
      iref( 1)  =  1
      idimp(1, 1) =  1
      idimp(2, 1) =  2
      iflip(1, 1) =  1
      iflip(2, 1) =  1
      iref( 2)  =  2
      idimp(1, 2) =  1
      idimp(2, 2) =  2
      iflip(1, 2) =  1
      iflip(2, 2) =  1
      iref( 3)  =  1
      idimp(1, 3) =  1
      idimp(2, 3) =  2
      iflip(1, 3) = -1
      iflip(2, 3) =  1
      iref( 4)  =  2
      idimp(1, 4) =  2
      idimp(2, 4) =  1
      iflip(1, 4) =  1
      iflip(2, 4) =  1
      iref( 5)  =  3
      idimp(1, 5) =  1
      idimp(2, 5) =  2
      iflip(1, 5) =  1
      iflip(2, 5) =  1
      iref( 6)  =  2
      idimp(1, 6) =  2
      idimp(2, 6) =  1
      iflip(1, 6) = -1
      iflip(2, 6) =  1
      iref( 7)  =  1
      idimp(1, 7) =  1
      idimp(2, 7) =  2
      iflip(1, 7) =  1
      iflip(2, 7) = -1
      iref( 8)  =  2
      idimp(1, 8) =  1
      idimp(2, 8) =  2
      iflip(1, 8) =  1
      iflip(2, 8) = -1
      iref( 9)  =  1
      idimp(1, 9) =  1
      idimp(2, 9) =  2
      iflip(1, 9) = -1
      iflip(2, 9) = -1
 
       return
       end
 
 
       subroutine loadsyms2dbtos(iref,idimp,iflip)
       implicit real *8 (a-h,o-z)
 
       dimension iref(*), idimp(2,*), iflip(2,*)
 
      iref( 1)  =  1
      idimp(1, 1) =  1
      idimp(2, 1) =  2
      iflip(1, 1) =  1
      iflip(2, 1) =  1
      iref( 2)  =  2
      idimp(1, 2) =  1
      idimp(2, 2) =  2
      iflip(1, 2) =  1
      iflip(2, 2) =  1
      iref( 3)  =  2
      idimp(1, 3) =  1
      idimp(2, 3) =  2
      iflip(1, 3) = -1
      iflip(2, 3) =  1
      iref( 4)  =  1
      idimp(1, 4) =  1
      idimp(2, 4) =  2
      iflip(1, 4) = -1
      iflip(2, 4) =  1
      iref( 5)  =  2
      idimp(1, 5) =  2
      idimp(2, 5) =  1
      iflip(1, 5) =  1
      iflip(2, 5) =  1
      iref( 6)  =  2
      idimp(1, 6) =  2
      idimp(2, 6) =  1
      iflip(1, 6) = -1
      iflip(2, 6) =  1
      iref( 7)  =  2
      idimp(1, 7) =  2
      idimp(2, 7) =  1
      iflip(1, 7) =  1
      iflip(2, 7) = -1
      iref( 8)  =  2
      idimp(1, 8) =  2
      idimp(2, 8) =  1
      iflip(1, 8) = -1
      iflip(2, 8) = -1
      iref( 9)  =  1
      idimp(1, 9) =  1
      idimp(2, 9) =  2
      iflip(1, 9) =  1
      iflip(2, 9) = -1
      iref(10)  =  2
      idimp(1,10) =  1
      idimp(2,10) =  2
      iflip(1,10) =  1
      iflip(2,10) = -1
      iref(11)  =  2
      idimp(1,11) =  1
      idimp(2,11) =  2
      iflip(1,11) = -1
      iflip(2,11) = -1
      iref(12)  =  1
      idimp(1,12) =  1
      idimp(2,12) =  2
      iflip(1,12) = -1
      iflip(2,12) = -1
 
       return
       end
 
 
       subroutine loadsyms2dstob(iref,idimp,iflip)
       implicit real *8 (a-h,o-z)
 
       dimension iref(*), idimp(2,*), iflip(2,*)
 
      iref( 1)  =  1
      idimp(1, 1) =  1
      idimp(2, 1) =  2
      iflip(1, 1) =  1
      iflip(2, 1) =  1
      iref( 2)  =  2
      idimp(1, 2) =  1
      idimp(2, 2) =  2
      iflip(1, 2) =  1
      iflip(2, 2) =  1
      iref( 3)  =  2
      idimp(1, 3) =  1
      idimp(2, 3) =  2
      iflip(1, 3) = -1
      iflip(2, 3) =  1
      iref( 4)  =  1
      idimp(1, 4) =  1
      idimp(2, 4) =  2
      iflip(1, 4) = -1
      iflip(2, 4) =  1
      iref( 5)  =  2
      idimp(1, 5) =  2
      idimp(2, 5) =  1
      iflip(1, 5) =  1
      iflip(2, 5) =  1
      iref( 6)  =  2
      idimp(1, 6) =  2
      idimp(2, 6) =  1
      iflip(1, 6) = -1
      iflip(2, 6) =  1
      iref( 7)  =  2
      idimp(1, 7) =  2
      idimp(2, 7) =  1
      iflip(1, 7) =  1
      iflip(2, 7) = -1
      iref( 8)  =  2
      idimp(1, 8) =  2
      idimp(2, 8) =  1
      iflip(1, 8) = -1
      iflip(2, 8) = -1
      iref( 9)  =  1
      idimp(1, 9) =  1
      idimp(2, 9) =  2
      iflip(1, 9) =  1
      iflip(2, 9) = -1
      iref(10)  =  2
      idimp(1,10) =  1
      idimp(2,10) =  2
      iflip(1,10) =  1
      iflip(2,10) = -1
      iref(11)  =  2
      idimp(1,11) =  1
      idimp(2,11) =  2
      iflip(1,11) = -1
      iflip(2,11) = -1
      iref(12)  =  1
      idimp(1,12) =  1
      idimp(2,12) =  2
      iflip(1,12) = -1
      iflip(2,12) = -1
 
       return
       end
 
 
