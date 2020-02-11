module ramancal
implicit none
integer, parameter :: dp = selected_real_kind(15, 307)
real(kind = dp), parameter :: THzToHartree = 1.0_dp/6579.683920729_dp
real(kind = dp), parameter :: HartreeToEv  = 27.21138386_dp
real(kind = dp), parameter :: eVToHartree  = 1.0_dp/27.21138386_dp
real(kind = dp), parameter :: tpi =2.0_dp*3.14159265358979323846
character(len=10), parameter :: output = 'status.out'
character(len=8), parameter :: input  = 'input.in'
real(kind = dp) :: kT, temperature
character(len = 256) :: phononsInput
real(kind = dp), allocatable :: atomD(:,:), atomM(:), phonQ(:,:), phonF(:), phonD(:,:,:,:), Omegaj(:)
real(kind = dp), allocatable :: genCoord(:), Sj(:), nj(:)
integer :: nOfqPoints, nAtoms, nModes

namelist /ramanInput/ temperature, phononsInput

contains

subroutine readInputs()
implicit none
logical :: file_exists
! Check output file
inquire(file=output, exist=file_exists)
if ( file_exists ) then
  open(unit=12, file=output, status='old')
  close(unit=12, status='delete')
end if
! Read input file parameters
open(13, file=input)
read(13, ramanInput)
close(13)
! Update kT unit from K
kT = temperature*8.6173324e-5_dp*eVToHartree
! Read phonon information
call readPhonons()
end subroutine

subroutine readPhonons()
implicit none
integer :: iAtom, iMode, iq
real(kind = dp) :: dummyD, freqInTHz, tempA(3)
CHARACTER :: dummyC
open(1, file=trim(phononsInput), status="old")
read(1,*) nOfqPoints, nAtoms, nModes
read (1,*)
allocate( atomD(3,nAtoms), atomM(nAtoms) )
atomD = 0.0_dp
atomM = 0.0_dp
!The unit of atomD is in Bohr, converting is needed
do iAtom = 1, nAtoms
  read(1,*) atomD(1,iAtom), atomD(2,iAtom), atomD(3,iAtom), atomM(iAtom)
enddo
read(1,*)
do iAtom = 1, nAtoms
  read(1,*) tempA(1),tempA(2),tempA(3)
  atomD(:,iAtom)=tempA(:)-atomD(:,iAtom)
enddo
! A to Bohr
atomD=atomD*1.889725989
read(1,*)
allocate( phonQ(3,nOfqPoints), phonF(nModes), Omegaj(nModes), phonD(3,nAtoms,nModes,nOfqPoints) )
phonQ = 0.0_dp
phonF = 0.0_dp
phonD = 0.0_dp
do iq = 1, nOfqPoints
  read (1,*) dummyC, dummyC, dummyC, phonQ(1,iq), phonQ(2,iq), phonQ(3,iq), dummyC
  read(1,*)
  do iMode = 1, nModes
    read(1,*) dummyC, dummyC, dummyC, dummyC, freqInTHz
    phonF(iMode) = dble(freqInTHz)*THzToHartree
    Omegaj(iMode) = dble(freqInTHz)*tpi
    do iAtom = 1, nAtoms
      read(1,*) dummyC, phonD(1,iAtom,iMode,iq), dummyD, phonD(2,iAtom,iMode,iq), dummyD, phonD(3,iAtom,iMode,iq), dummyC
    enddo
  enddo
enddo
close(1)
return
end subroutine

subroutine computeSj()
implicit none
integer :: iq, iMode, iAtom
allocate( genCoord(nModes) )
allocate( Sj(nModes) )
allocate( nj(nModes) )
do iq = 1, nOfqPoints
  do iMode = 1, nModes
    genCoord(iMode) = 0.0_dp
    do iAtom = 1, nAtoms
      genCoord(iMode) = genCoord(iMode) + sqrt(1822.88833218_dp*atomM(iAtom))*sum(phonD(:,iAtom,iMode,iq)*atomD(:,iAtom))
    enddo
  enddo
enddo
Sj(:) = 0.5_dp*phonF(:)*genCoord(:)*genCoord(:)
do iMode = 1, nModes
  nj(iMode)=1.0d0/(exp(phonF(iMode)/kT)-1.0d0)
enddo
return
end subroutine

subroutine outputSj()
implicit none
integer :: iq, iMode
real(kind=dp) :: S1, S2
open(14,file='Sj.out')
!~ write(14,'("iMode      Sj          nj          S+          S-")')
write(14,'("iMode      Omegaj")')
write(14, '(i5)')nModes
do iq = 1, nOfqPoints
  do iMode = 1, nModes
!~     S1=(nj(iMode)+1.0d0)*Sj(iMode)
!~     S2=nj(iMode)*Sj(iMode)
    write(14,'(i5,4ES12.4E2)')iMode, Sj(iMode), Omegaj(iMode)
    write(*,'(i5,3ES16.6E2)')iMode, Sj(iMode), Omegaj(iMode)
!~     , nj(iMode), S1, S2
  enddo
enddo
close(14)
end subroutine

subroutine computeInt()
implicit none
end subroutine

subroutine finialize()
implicit none
deallocate(atomD, atomM)
deallocate(phonQ, phonF, phonD)
deallocate(genCoord,Sj,nj)
end subroutine

end module
  
program raman
use ramancal
implicit none
call readInputs()
call computeSj()
call outputSj()
!~ call computeInt()
print *,"job done"
call finialize()
end program
