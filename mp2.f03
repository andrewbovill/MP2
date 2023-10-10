INCLUDE 'mp2_mod.f03'
      Program mp2
!
!     This program reads AO integrals from a Gaussian matrix file and times
!     AO-to-MO integral transformations.
!
!     Original Written by
!     -H. P. Hratchian, 2022.
!
!     Modified by
!     -A. J. Bovill, 2023.
!
!     USE Connections
!
      use mp2_mod
!
!     Variable Declarations
!
      implicit none
      integer(kind=int64)::nCommands,iPrint=0,nOMP,nAtoms,nBasis,  &
        nBasisUse,nElectrons,nElectronsAlpha,nElectronsBeta
      integer(kind=int64)::mu,nu,lambda,sigma,p,q,r,s,pq,rs,pqrs,iCount,i,j,a,b
      real(kind=real64)::timeStart,timeEnd,time0,time1,tmpReal,  &
        deltaIJAB,numerator,E2AA,E2BB,E2AB,E2BA
      real(kind=real64),dimension(:),allocatable::moEnergiesAlpha,moEnergiesBeta
      real(kind=real64),dimension(:,:),allocatable::CAlpha,CBeta,  &
        tmpMatrix1,tmpMatrix2
      real(kind=real64),dimension(:,:,:,:),allocatable::aoInts,moInts,  &
        partialInts1,partialInts2
      type(MQC_Variable)::ERIs,mqcTmpArray
      character(len=512)::tmpString,matrixFilename
      type(mqc_gaussian_unformatted_matrix_file)::GMatrixFile
      logical::fail=.false.,doN8=.false.,doSlowN5=.true.,  &
        doRegularN5=.true.,useBLAS=.true.
!
!     Andrew debug (all integers passed into MQC subroutines must be "kind=int64")
!
      integer(kind=int64):: newerThanMajor,newerThanMinor,newerThanRevision 
!
!     Format Statements
!
 1000 Format(1x,'Enter Test Program mp2.')
 1010 Format(1x,'Matrix File: ',A,/)
 1020 Format(1x,'Use ',I3,' shared memory processors.')
 1100 Format(1x,'nAtoms    =',I4,6x,'nBasis  =',I4,6x,'nBasisUse=',I4,/,  &
             1x,'nElectrons=',I4,6x,'nElAlpha=',I4,6x,'nElBeta  =',I4)
 1500 Format(/,1x,'Carrying out O(N^8) transformation.')
 1510 Format(/,1x,'Skipping O(N^8) transformation.')
 2000 Format(1x,'<',I3,',',I3,' || ',I3,',',I3,' > ... pq=',I3,'  rs=',I3,'  pqrs=',I3)
 3000 Format(/,1x,'E(2)-SS = ',f15.10,' a.u.',4x,'E(2)-OS = ',f15.10,' a.u.')
 5000 Format(1x,'Time (',A,'): ',f8.1,' s.')
 8999 Format(/,1x,'END OF PROGRAM mp2.')
!
!
      call cpu_time(timeStart)
!
!     Open the Gaussian matrix file and load the number of atomic centers.

      nCommands = command_argument_count()
      if(nCommands.eq.0)  &
        call mqc_error('No command line arguments provided. The input Gaussian matrix file name is required.')
      call get_command_argument(1,tmpString)
      call  commandLineArgs(iPrint,nOMP,matrixFilename,doN8,  &
        doSlowN5,doRegularN5,useBLAS,fail)
      call omp_set_num_threads(nOMP)
      call GMatrixFile%load(matrixFilename)
      write(IOut,1010) TRIM(matrixFilename)
      write(iOut,1020) nOMP
      nAtoms = GMatrixFile%getVal('nAtoms')
      nBasis = Int(GMatrixFile%getVal('nbasis'))
      nBasisUse = Int(GMatrixFile%getVal('nbasisuse'))
      nElectrons = Int(GMatrixFile%getVal('nelectrons'))
      nElectronsAlpha = Int(GMatrixFile%getVal('nAlpha'))
      nElectronsBeta = Int(GMatrixFile%getVal('nBeta'))
      write(IOut,1100) nAtoms,nBasis,nBasisUse,nElectrons,  &
        nElectronsAlpha,nElectronsBeta
!
!     Load the orbital eigenvalues.
!
      call GMatrixFile%getArray('ALPHA ORBITAL ENERGIES',mqcVarOut=mqcTmpArray)
      call mqcTmpArray%print(header='Alpha MO Energies')
      flush(iOut)
      moEnergiesAlpha = mqcTmpArray
      if(GMatrixFile%isUnrestricted()) then
        call mqc_error('UHF/UKS NYI.')
      else
        moEnergiesBeta = moEnergiesAlpha
      endIf
!
!     Load the MO coefficients.
!
      write(*,*)' Hrant -  isUnrestricted: ',GMatrixFile%isUnrestricted()
      call GMatrixFile%getArray('ALPHA MO COEFFICIENTS',mqcVarOut=mqcTmpArray)
      CAlpha = mqcTmpArray
      if(GMatrixFile%isUnrestricted()) then
        call GMatrixFile%getArray('BETA  MO COEFFICIENTS',mqcVarOut=mqcTmpArray)
        CBeta = mqcTmpArray
      else
        CBeta = CAlpha
      endIf
!
!     Read in and report out the (AO) ERIs.
!
      write(*,*) "Andrew check 1"
      call GMatrixFile%getArray('REGULAR 2E INTEGRALS',mqcVarOut=ERIs)
      if(iPrint.ge.2) call ERIs%print(IOut,' ERIs=')
!
!     Allocate space for AO and MO integrals. Then fill the intrinsic array of
!     AO integrals.
!
      allocate(aoInts(nBasis,nBasis,nBasis,nBasis),  &
        moInts(nBasisUse,nBasisUse,nBasisUse,nBasisUse))
      flush(iOut)

!hph+
      write(*,*) "Andrew check 1"
      aoInts = ERIs
      call dpReshape4(nBasis,nBasis,nBasis,nBasis,ERIs%realArray,aoInts)
!hph-
      
      write(*,*) "Andrew check 1"
      flush(iOut)
!
!     Do N^8 AO --> MO transformation.
!     Andrew-- only same spin atm.
!
      call cpu_time(time0)
      call integralTransformationN8sameSpin(nBasis,nBasisUse,CAlpha,aoInts,moInts)
      call cpu_time(time1)
      write(iOut,5000) 'Integral Transformation N8',time1-time0

!     call cpu_time(time0)
!     call integralTransformationN8sameSpin(nBasis,nBasisUse,CAlpha,aoInts,moInts)
!     call cpu_time(time1)

!
!     Evaluate the E(2) AA, BB, AB, and BA contribution.
!

      call cpu_time(time0)
      call E2(nBasis,nBasisUse,nElectronsAlpha,nElectronsBeta, & 
        moInts,moEnergiesAlpha,moEnergiesBeta,E2AA,E2AB) 
      call cpu_time(time1)
      write(iOut,5000) 'E2',time1-time0

      flush(iOut)
!
!
!     Load up AA MO ERIs from the matrixfile and print them out to ensure the
!     explicit AO-->MO transformations above gave the right answers.
!

!
!     ^^^ Andrew add later ^^^ (Ask Hrant,doesnt work) ¯\_(ツ)_/¯
!

      E2BB = E2AA
      E2BA = E2AB
      write(iOut,3000) E2AA,E2AB
      write(*,*)' Total E2 = ',E2AA+E2BB+E2AB
!
  999 Continue
      call cpu_time(timeEnd)
      write(iOut,5000) 'TOTAL JOB TIME',timeEnd-timeStart
      write(iOut,8999)
      end program mp2

