INCLUDE 'mp2_mod.f03'
      Program mp2
!
!     This program reads AO integrals from a Gaussian matrix file and times
!     AO-to-MO integral transformations.
!
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
      call GMatrixFile%getArray('REGULAR 2E INTEGRALS',mqcVarOut=ERIs)
      if(iPrint.ge.2) call ERIs%print(IOut,' ERIs=')
!
!     Allocate space for AO and MO integrals. Then fill the intrinsic array of
!     AO integrals.
!
      allocate(aoInts(nBasis,nBasis,nBasis,nBasis),  &
        moInts(nBasisUse,nBasisUse,nBasisUse,nBasisUse))

      flush(iOut)

!
!          
!

      aoInts = ERIs
      call dpReshape4(nBasis,nBasis,nBasis,nBasis,ERIs%realArray,aoInts)

      if(iPrint.ge.2) call mqc_print_rank4Tensor_array_real(aoInts,IOut,header='Intrinsic AO Integrals')

      flush(iOut)

!
!     Do N^8 AO --> MO transformation.
!
      call integralTransformationN8sameSpin(nBasis,nBasisUse,CAlpha,aoInts,moInts)
      flush(iOut)
      
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
              if(iPrint.ge.1) write(*,*)' num, denom = ',numerator,deltaIJAB
              E2AA = E2AA + numerator/(float(4)*deltaIJAB)
            endDo
          endDo
        endDo
      endDo
      call cpu_time(time1)
      write(iOut,5000) 'Same Spin E2',time1-time0
      E2BB = E2AA
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
              if(iPrint.ge.1) write(*,*)' num, denom = ',numerator,deltaIJAB
              E2AB = E2AB + numerator/deltaIJAB
            endDo
          endDo
        endDo
      endDo
      call cpu_time(time1)
      write(iOut,5000) 'Opposite Spin E2',time1-time0
      E2BA = E2AB

      write(iOut,*) "E2AA",E2AA
      write(iOut,*) "E2BB",E2BB
      write(iOut,*) "E2AB",E2AB
      write(iOut,*) "E2BA",E2BA

      write(*,*)' Total E2 = ',E2AA+E2BB+E2AB
!
  999 Continue
      call cpu_time(timeEnd)
      write(iOut,5000) 'TOTAL JOB TIME',timeEnd-timeStart
      write(iOut,8999)
      end program mp2

      subroutine dpReshape4(N1,N2,N3,N4,arrayIn,r4ArrayOut)
!
      use iso_fortran_env
      implicit none
      integer(kind=int64)::N1,N2,N3,N4
      real(kind=real64),dimension(N1,N2,N3,N4)::arrayIn,r4ArrayOut

      write(*,*) "Andrew within subroutine"
      r4ArrayOut = arrayIn
!
      return
      end subroutine dpReshape4
