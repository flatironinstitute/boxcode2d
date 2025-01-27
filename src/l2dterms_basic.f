c
c      
      
      subroutine l2dterms(eps,nterms,ier)
c
c     determine number of multipole terms to meet given
c     precision, using the formula
c
c     ( sqrt(2)/2 )^nterms/ (3/2)^(nterms+1) < eps
c      

      implicit none
      real *8 eps, tmp
      integer nterms, ier

      ier = 0

      tmp = log(1.5d0*eps)/log(sqrt(2d0)/3d0)
      
      nterms = tmp
      if (nterms .lt. tmp) nterms = nterms + 1

      nterms = max(nterms,1)

      return
      end
