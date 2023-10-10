      module mp2_mod
!
!     This module supports the program mp2.
!
!     Original Written by
!     -H. P. Hratchian, 2022.
!
!     Modified by
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

!PROCEDURE commandLineArgs
      subroutine commandLineArgs(iPrint,nOMP,matrixFilename,doN8,  &
        doSlowN5,doRegularN5,useBLAS,fail)
!
!     This subroutine is used to process the command line arguments.
!
      integer(kind=int64),intent(OUT)::iPrint,nOMP
      character(len=512),intent(OUT)::matrixFilename
      logical,intent(OUT)::doN8,doSlowN5,doRegularN5,useBLAS,fail
!
      integer::nCommands,nMatrixFilenames,i
      character(len=512)::tmpString,lowercase
      logical::getNProc
!
!     Format statements.
!
 9000 format(1x,'Unknown command line switch: ',A,'.')
 9100 format(1x,'More than 1 matrix filename found!')
!
!     Set defaults.
!
      iPrint = 0
      nOMP = 1
      doN8 = .false.
      doSlowN5 = .false.
      doRegularN5 = .false.
      useBLAS = .true.
      fail = .false.
      nMatrixFilenames = 0
      getNProc = .false.
!
!     Determine the number of command line arguments. Then, loop through the
!     list of arguments to set options.
!
      nCommands = command_argument_count()
      do i = 1,nCommands
        call get_command_argument(i,tmpString)
        if(getNProc) then
          read(tmpString,'(I3)') nOMP
          getNProc = .false.
          cycle
        endIf
        if(tmpString(1:1).eq.'-') then
          call String_Change_Case(tmpString(2:),'l',lowercase)
          select case(TRIM(lowercase))
          case('debug')
            iPrint = 1
          case('nproc')
            getNProc = .true.
          case('don8')
            doN8 = .true.
          case('skipn8')
            doN8 = .false.
          case('doslown5')
            doSlowN5 = .true.
          case('skipslown5')
            doSlowN5 = .false.
          case('doregularn5')
            doRegularN5 = .true.
          case('skipregularn5')
            doRegularN5 = .false.
          case('useblas')
            useBLAS = .true.
          case('usematmul')
            useBLAS = .false.
          case default
            fail = .true.
            write(iOut,9000) TRIM(tmpString)
            return
          endSelect
        else
          nMatrixFilenames = nMatrixFilenames + 1
          if(nMatrixFilenames.gt.1) then
            fail = .true.
            write(iOut,9100)
            return
          endIf
          matrixFilename = tmpString
        endIf
      endDo
!
      return
      end subroutine commandLineArgs

!
!PROCEDURE integralTransformationN8
      subroutine integralTransformationN8sameSpin(nBasis,nBasisUse,C,aoInts,moInts)
!
!     This subroutine carries out AO-->MO integral transformations for a same-spin case.
!
!     Written by 
!     H. P. Hratchian, 2022
!
!     Modified by 
!     A. J. Bovill, 2023
!
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
      write(*,*) "Andrew got here"
!
      return
      end subroutine integralTransformationN8sameSpin

!PROCEDURE integralTransformationN5
      subroutine integralTransformationN5sameSpin(nBasis,nBasisUse,C,aoInts,moInts)
!
!     This subroutine carries out AO-->MO integral transformations for a same-spin case.
!
!     Written by 
!     H. P. Hratchian, 2022
!
!     Modified by 
!     A. J. Bovill, 2023
!
      implicit none
      integer(kind=int64),intent(in)::nBasis,nBasisUse
      real(kind=real64),dimension(nBasis,nBasisUse)::C
      real(kind=real64),dimension(nBasis,nBasis,nBasis,nBasis)::aoInts
      real(kind=real64),dimension(nBasisUse,nBasisUse,nBasisUse,nBasisUse)::moInts
      integer(kind=int64)::p,q,r,s,mu,nu,lambda,sigma
!
!     Do the transformation with the Noddy O(N^5) algorithm.
!
!     if(doRegularN5) then
!       call cpu_time(time0)
!       do nu = 1,nBasis
!         do lambda = 1,nBasis
!           do sigma = 1,nBasis
!             do p = 1,nBasisUse
!               tmpReal = float(0)
!               do mu = 1,nBasis
!                 partialInts1(p,nu,lambda,sigma) = partialInts1(p,nu,lambda,sigma)  &
!                   + CAlpha(mu,p)*aoInts(mu,nu,lambda,sigma)
!               endDo
!             endDo
!           endDo
!         endDo
!       endDo
!       call cpu_time(time1)
!       write(iOut,5000) 'Quarter Transformation 1b',time1-time0
!       flush(iOut)
!     endIf
!
!     Allocate(partialInts2(nBasisUse,nBasisUse,nBasis,nBasis))
!     if(doSlowN5) then
!       call cpu_time(time0)
!       do p = 1,nBasisUse
!         do q = 1,nBasisUse
!           do nu = 1,nBasis
!             do lambda = 1,nBasis
!               do sigma = 1,nBasis
!                 partialInts2(p,q,lambda,sigma) = partialInts2(p,q,lambda,sigma)  &
!                   + CAlpha(nu,q)*partialInts1(p,nu,lambda,sigma)
!               endDo
!             endDo
!           endDo
!         endDo
!       endDo
!       call cpu_time(time1)
!       write(iOut,5000) 'Quarter Transformation 2a',time1-time0
!       flush(iOut)
!       partialInts2 = float(0)
!     endIf
!
!     write(*,*) "Andrew check 3"
!     if(doRegularN5) then
!       call cpu_time(time0)
!       do lambda = 1,nBasis
!         do sigma = 1,nBasis
!           do q = 1,nBasisUse
!             do p = 1,nBasisUse
!               do nu = 1,nBasis
!                 partialInts2(p,q,lambda,sigma) = partialInts2(p,q,lambda,sigma)  &
!                   + CAlpha(nu,q)*partialInts1(p,nu,lambda,sigma)
!               endDo
!             endDo
!           endDo
!         endDo
!       endDo
!       call cpu_time(time1)
!       write(iOut,5000) 'Quarter Transformation 2b',time1-time0
!       flush(iOut)
!     endIf
!
!     DeAllocate(partialInts1)
!     Allocate(partialInts1(nBasisUse,nBasisUse,nBasisUse,nBasis))
!     partialInts1 = float(0)
!
!     if(doSlowN5) then
!       call cpu_time(time0)
!       do p = 1,nBasisUse
!         do q = 1,nBasisUse
!           do r = 1,nBasisUse
!             do lambda = 1,nBasis
!               do sigma = 1,nBasis
!                 partialInts1(p,q,r,sigma) = partialInts1(p,q,r,sigma)  &
!                   + CAlpha(lambda,r)*partialInts2(p,q,lambda,sigma)
!               endDo
!             endDo
!           endDo
!         endDo
!       endDo
!       call cpu_time(time1)
!       write(iOut,5000) 'Quarter Transformation 3a',time1-time0
!       flush(iOut)
!       partialInts1 = float(0)
!     endIf
!
!     if(doRegularN5) then
!       call cpu_time(time0)
!       do p = 1,nBasisUse
!         do q = 1,nBasisUse
!               do sigma = 1,nBasis
!           do r = 1,nBasisUse
!             do lambda = 1,nBasis
!                 partialInts1(p,q,r,sigma) = partialInts1(p,q,r,sigma)  &
!                   + CAlpha(lambda,r)*partialInts2(p,q,lambda,sigma)
!               endDo
!             endDo
!           endDo
!         endDo
!       endDo
!       call cpu_time(time1)
!       write(iOut,5000) 'Quarter Transformation 3b',time1-time0
!       flush(iOut)
!     endIf
!
!     DeAllocate(partialInts2)
!     moInts = float(0)
!
!     Error in the loop below
!
!     write(*,*) "Andrew here"
!     if(doRegularN5.or.doSlowN5) then
!     write(*,*) "Andrew here"
!       call cpu_time(time0)
!       do p = 1,nBasisUse
!         do q = 1,nBasisUse
!           do r = 1,nBasisUse
!             do s = 1,nBasisUse
!               do sigma = 1,nBasis
!                 moInts(p,q,r,s) = moInts(p,q,r,s)  &
!                   + CAlpha(sigma,s)*partialInts1(p,q,r,sigma)
!               endDo
!             endDo
!           endDo
!         endDo
!       endDo
!       call cpu_time(time1)
!       write(iOut,5000) 'Quarter Transformation 4',time1-time0
!       flush(iOut)
!       DeAllocate(partialInts1)
 
      end subroutine integralTransformationN5sameSpin

      subroutine E2(nBasis,nBasisUse,nElectronsAlpha,nElectronsBeta, & 
        moInts,moEnergiesAlpha,moEnergiesBeta,E2AA,E2AB) 
      real(8),dimension(:),allocatable::moEnergiesAlpha,moEnergiesBeta
      real(8),dimension(:,:,:,:),allocatable::moInts  
      integer(8),intent(in)::nBasis,nBasisUse,nElectronsAlpha,nElectronsBeta
      integer(8)::i,j,a,b
      real(8):: deltaIJAB, numerator,time0,time1
      real(8),intent(out):: E2AA,E2AB
!
!     Evaluate the E(2) AA, BB, AB, and BA contribution.
!
      write(*,*)
      write(*,*)' Same Spin E2...'
      E2AA = float(0)
      call cpu_time(time0)
      do i = 1,nElectronsAlpha
        do j = 1,nElectronsAlpha
          do a = nElectronsAlpha+1,nBasisUse
            do b = nElectronsAlpha+1,nBasisUse
              deltaIJAB = moEnergiesAlpha(i) + moEnergiesAlpha(j)  &
                - moEnergiesAlpha(a) - moEnergiesAlpha(b)
              numerator = moInts(i,a,j,b) - moInts(i,b,j,a)
              numerator = numerator*numerator
              E2AA = E2AA + numerator/(float(4)*deltaIJAB)
            endDo
          endDo
        endDo
      endDo
      call cpu_time(time1)

      write(*,*)
      write(*,*)' Opposite Spin E2...'
      E2AB = float(0)
      call cpu_time(time0)
      do i = 1,nElectronsAlpha
        do j = 1,nElectronsBeta
          do a = nElectronsAlpha+1,nBasisUse
            do b = nElectronsBeta+1,nBasisUse
              deltaIJAB = moEnergiesAlpha(i) + moEnergiesBeta(j)  &
                - moEnergiesAlpha(a) - moEnergiesBeta(b)
              numerator = moInts(i,a,j,b)*moInts(i,a,j,b)
              E2AB = E2AB + numerator/deltaIJAB
            endDo
          endDo
        endDo
      endDo
      call cpu_time(time1)
      end subroutine E2

      subroutine dpReshape4(N1,N2,N3,N4,arrayIn,r4ArrayOut)
!
      use iso_fortran_env
      implicit none
      integer(kind=int64)::N1,N2,N3,N4
      real(kind=real64),dimension(N1,N2,N3,N4)::arrayIn,r4ArrayOut
!
      r4ArrayOut = arrayIn
!
      return
      end subroutine dpReshape4
!
!
      end module mp2_mod
