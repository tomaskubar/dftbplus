!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Chemical potential equilibration
!>
module dftbp_cpecalc

  use dftbp_accuracy
  use dftbp_assert
  use dftbp_lapackroutines
  use dftbp_message
  use dftbp_cpeinp

  implicit none

  private

  public :: TCpeCalc

  !> Calculation of CPE
  type :: TCpeCalc

    !> Number of species
    integer :: nSpecies

    !> Electronegativity per species, size (nSpecies)
    real(dp), allocatable :: electronegativity(:)

    !> Chemical hardness per species, size (nSpecies)
    real(dp), allocatable :: hardness(:)

    !> Covalent radius per species, size (nSpecies)
    real(dp), allocatable :: radius(:)

    !> Total charge of the molecule
    real(dp) :: totalCharge

    !> Number of atoms
    integer :: nAtom
    
    !> Species for each atom, size (nAtom)
    integer, allocatable :: species(:)

    !> Chemical symbol for each species, size (nSpecies)
    character(mc), allocatable :: speciesName(:)

    !> Coordinates of atoms, size (3, nAtom)
    real(dp), allocatable :: coord(:,:)

    !> .true. if the number of atoms, species and coordinates have been set
    logical :: moleculeIsKnown

    !> Atomic charges
    real(dp), allocatable :: charge(:)

  contains

    !> Assign atom number, coordinates and species
    procedure :: init

    !> Assign atom number, coordinates and species
    procedure :: setup

    !> Run the actual calculation
    procedure :: calculate

    !> Evaluate the gamma matrix
    procedure :: calcGamma

  end type TCpeCalc

  contains

  !> Initialize CPE data from CPE input
  subroutine init(this, nType, inp, nAtom, species, speciesName, coord)
    
    !> data type for CPE
    class(TCpeCalc), intent(inout) :: this

    !> number of species
    integer, intent(in) :: nType

    !> data type for CPE input
    type(TCpeInp), intent(in) :: inp

    !> number of atoms
    integer, intent(in), optional :: nAtom

    !> species of atoms
    integer, intent(in), optional :: species(:)

    !> names of species
    character(mc), intent(in), optional :: speciesName(:)

    !> coordinates of atoms
    real(dp), intent(in), optional :: coord(:,:)

    this%nSpecies = nType

    allocate(this%electronegativity(this%nSpecies))
    this%electronegativity = inp%electronegativity

    allocate(this%hardness(this%nSpecies))
    this%hardness = inp%hardness

    allocate(this%radius(this%nSpecies))
    this%radius = inp%radius

    this%totalCharge = inp%totalCharge

    if (present(nAtom)) then
      this%nAtom = nAtom

      allocate(this%species(this%nAtom))
      this%species = species

      allocate(this%speciesName(this%nSpecies))
      this%speciesName = speciesName

      allocate(this%coord(3, this%nAtom))
      this%coord = coord

      allocate(this%charge(this%nAtom+1))

      this%moleculeIsKnown = .true.
    else
      this%moleculeIsKnown = .false.
    end if

  end subroutine init

  !> Set up CPE calculation from information about the system
  subroutine setup(this, nAtom, species, coord)
    
    !> data type for CPE
    class(TCpeCalc), intent(inout) :: this

    !> number of atoms
    integer, intent(in) :: nAtom

    !> species of atoms
    integer, intent(in) :: species(:)

    !> coordinates of atoms
    real(dp), intent(in) :: coord(:,:)

  ! @:ASSERT(size(species) == nAtom)
  ! @:ASSERT(size(coord, dim=1) == 3)
  ! @:ASSERT(size(coord, dim=2) == nAtom)

    if (.not. this%moleculeIsKnown) then
      this%nAtom = nAtom
      
      allocate(this%species(nAtom))
      this%species = species
      
      allocate(this%coord(3, nAtom))
      this%coord = coord

      allocate(this%charge(this%nAtom+1))

      this%moleculeIsKnown = .true.
    end if

  end subroutine setup

  !> Carry out the CPE calculation
  subroutine calculate(this)
    
    !> data type for CPE
    class(TCpeCalc), intent(inout) :: this

    integer :: iSp, iAt, jAt, iDir
    real(dp), allocatable :: gammaMat(:,:), atomElectronegativity(:,:)

    write (*,*) 'speciesNames'
    do iSp=1, this%nSpecies
      write (*,*) iSp, this%speciesName(iSp)
    end do

    do iAt=1, this%nAtom
      write (*,*) this%species(iAt), (this%coord(iDir, iAt), iDir=1, 3) 
    end do

    allocate(gammaMat(this%nAtom+1, this%nAtom+1))
    call this%calcGamma(gammaMat)
    write (*,*) "gamma matrix"
    write (*,'(13F8.4)') gammaMat

    allocate(atomElectronegativity(this%nAtom+1, 1))
    do iAt=1, this%nAtom
      atomElectronegativity(iAt, 1) = - this%electronegativity(this%species(iAt))
    end do
    atomElectronegativity(this%nAtom+1, 1) = this%totalCharge
    write (*,*) "atom electronegativities"
    write (*,'(F8.4)') atomElectronegativity

    call symmatinv(gammaMat)
    do iAt=2, this%nAtom+1
      do jAt=1, iAt-1
        gammaMat(jAt, iAt) = gammaMat(iAt, jAt)
      end do
    end do
    write (*,*) "inverted matrix"
    write (*,'(13F8.4)') gammaMat

    this%charge = matmul(gammaMat, atomElectronegativity(:,1))

  ! call dtrsm('L', 'L', 'N', 'N', this%nAtom+1, 1, 1._dp, gammaMat, this%nAtom+1, atomElectronegativity, this%nAtom+1)

    write (*,*) "resulting charges"
    write (*,'(F8.4)') this%charge

  end subroutine calculate

  subroutine calcGamma(this, gammaMat)

    !> data type for CPE
    class(TCpeCalc), intent(inout) :: this

    !> result -- gamma matrix
    real(dp), intent(out) :: gammaMat(:,:)

    integer :: iAt, jAt
    real(dp) :: radius1, radius2, averageRadius, gammaValue, distance

    do iAt = 1, this%nAtom
      radius1 = this%radius(this%species(iAt))
      gammaMat(iAt, iAt) = this%hardness(this%species(iAt)) + &
                         & 1._dp / sqrt(4._dp * atan(1._dp)) / radius1 
      do jAt = 1, iAt-1
        radius2 = this%radius(this%species(jAt))
        averageRadius = sqrt(radius1**2 + radius2**2)
        distance = sqrt((this%coord(1,iAt) - this%coord(1,jAt))**2 + &
                      & (this%coord(2,iAt) - this%coord(2,jAt))**2 + &
                      & (this%coord(2,iAt) - this%coord(2,jAt))**2)
        gammaValue = erf(distance / sqrt(2._dp) / averageRadius) / distance
        gammaMat(iAt, jAt) = gammaValue
        gammaMat(jAt, iAt) = gammaValue
      end do
    end do

    do iAt = 1, this%nAtom
      gammaMat(this%nAtom+1, iAt) = 1._dp
      gammaMat(iAt, this%nAtom+1) = 1._dp
    end do
    gammaMat(this%nAtom+1, this%nAtom+1) = 0._dp

  end subroutine calcGamma

! function calcGammaValue(hard1, hard2, distance)
!   !> chemical hardness of both atoms
!   real(dp), intent(in) :: hard1, hard2
!   !> distance of the atoms
!   real(dp), intent(in) :: distance
!   !> result
!   real(dp), intent(out) :: calcGammaValue
!
!   calcGammaValue = ...
! end function calcGammaValue

end module dftbp_cpecalc
