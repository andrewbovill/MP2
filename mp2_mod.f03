      module mp2_mod
!
!     This module supports the program mp2.
!
!     -A. J. Bovill, 2023.
!
!
!     USE Connections
!
      use mqc_general
      use mqc_molecule
      use mqc_gaussian
      use mqc_algebra2
      use mqc_algebra
      use iso_fortran_env
      use OMP_LIB
!
!     Variable Declarations
!
      implicit none
      integer(kind=int64),parameter::IOut=6
!
!     Module Procedures
!
      CONTAINS
!
!PROCEDURE integralTransformationN8
      subroutine integralTransformationN8sameSpin(nBasis,nBasisUse,C,aoInts,moInts)
!
!     This subroutine carries out AO-->MO integral transformations for a same-spin case.
!
!     H. P. Hratchian, 2022
!
      implicit none
      integer(kind=int64),intent(in)::nBasis,nBasisUse
      real(kind=real64),dimension(nBasis,nBasisUse)::C
      real(kind=real64),dimension(nBasis,nBasis,nBasis,nBasis)::aoInts
      real(kind=real64),dimension(nBasisUse,nBasisUse,nBasisUse,nBasisUse)::moInts
      integer(kind=int64)::p,q,r,s,mu,nu,lambda,sigma
!
!     Do the transformation with the straightforward O(N^8) algorithm.
!
      do p = 1,nBasisUse
        do q = 1,nBasisUse
          do r = 1,nBasisUse
            do s = 1,nBasisUse
              do mu = 1,nBasis
                do nu = 1,nBasis
                  do lambda = 1,nBasis
                    do sigma = 1,nBasis
                      moInts(p,q,r,s) = moInts(p,q,r,s)  &
                        + C(mu,p)*C(nu,q)*C(lambda,r)*C(sigma,s)*aoInts(mu,nu,lambda,sigma)
                    endDo
                  endDo
                endDo
              endDo
            endDo
          endDo
        endDo
      endDo
!
      return
      end subroutine integralTransformationN8sameSpin
!
!
      end module mp2_mod
