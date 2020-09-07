#:include 'common.fypp'

!> Offers everything which is publicly available when dealing with machine learning.

module dftbp_machinelearning
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_constants
  use dftbp_machinelearning_sf
  use dftbp_machinelearning_nn
  implicit none
  public

  !> Types of the model based on machine learning
  type :: TMachineLearningInp

    !> Symmetry functions
    type(TMLSymmetryFunctionsInp), allocatable :: sf

    !> Neural network
    type(TMLNeuralNetInp), allocatable :: nn
  
  end type TMachineLearningInp

  !> Types of the model based on machine learning
  type :: TMachineLearning

    !> The input data
    type(TMachineLearningInp), pointer :: input

    !> Symmetry function related data
    type(TMLSymmetryFunctions) :: sf

    !> Neural network related data
    type(TMLNeuralNet) :: nn

  contains

    procedure :: init
    procedure :: getEnergy
    procedure :: addGradients

  end type TMachineLearning

contains

  !> Initialize with data from input
  subroutine init(this, input, nAt, nSp, species)

    !> instance
    class(TMachineLearning), intent(inout) :: this

    !> the input structure to be linked into here
    type(TMachineLearningInp), intent(in), target :: input

    !> number of atoms
    integer, intent(in) :: nAt

    !> number of species/elements
    integer, intent(in) :: nSp

    !> species of each atom
    integer, intent(in) :: species(:)

    write (*,*) "MACHINE LEARNING INIT"

    @:ASSERT(size(species) == nAt)

    this%input => input

    call this%sf%init(input%sf, nAt, nSp)

    call this%nn%init(input%nn, species)

    ! store this information:
    ! atomic indices, atomic element/species, and symmetry function hyperparameters

  end subroutine init

  !> Get energy contributions from machine learning
  function getEnergy(this, coords, img2CentCell) result(energy)

    !> instance
    class(TMachineLearning), intent(inout) :: this

    !> Current coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Updated mapping to central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Resulting energy contribution
    real(dp) :: energy

    !> Energy per atom -- TODO output that as well?
    real(dp), allocatable :: energyAtom(:)

    integer :: iAt

    write (*,*) "MACHINE LEARNING ENERGY"

    allocate(energyAtom(size(coords, dim=2)))

    ! calculate all of the symmetry functions
    call this%sf%prepare(coords)
    call this%sf%evaluate()

    ! feed those values into the neural net
    call this%nn%evaluate(this%sf%sf, energyAtom)

    ! the unit is kcal/mol throughout the machine learning code
    energyAtom(:) = energyAtom(:) * kcal_mol__Hartree

    energy = sum(energyAtom)

    write (*,'(A,F15.10)') "MACHINE_LEARNING_ENERGY ", energy

!   do iAt = 1, size(energyAtom)
!     write (*,'(F15.10)') energyAtom(iAt)
!   end do

    deallocate(energyAtom)

  end function getEnergy


  !> Gradient contribution from machine learning
  subroutine addGradients(this, derivs, img2CentCell)

    !> instance
    class(TMachineLearning), intent(inout) :: this

    !> Derivatives to add contribution to 
    real(dp), intent(inout) :: derivs(:,:)

    !> Updated mapping to central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Derivatives per atom -- TODO output that as well?
    real(dp), allocatable :: energyDerivsAtom(:,:,:)

    !> Derivative contributions to be added
    real(dp), allocatable :: derivsAdd(:,:)

    integer :: iAt

    write (*,*) "MACHINE LEARNING GRADIENTS"

    allocate(energyDerivsAtom(3, this%nn%nAt, this%nn%nAt))
    allocate(derivsAdd(3, this%nn%nAt))

    ! derivatives of the symmetry functions
    call this%sf%evaluateDerivs()

    ! feed those values into the neural net
    call this%nn%evaluateDerivs(this%sf%dsfdr, energyDerivsAtom)

    ! the unit is kcal/mol throughout the machine learning code
    ! also 1/AA is converted to 1/bohr
    energyDerivsAtom(:,:,:) = energyDerivsAtom(:,:,:) * kcal_mol__Hartree * Bohr__AA

    ! sum contributions from all of the neural nets (one for each atom)
    derivsAdd = sum(energyDerivsAtom, 3)

    derivs = derivs + derivsAdd

    do iAt = 1, size(derivsAdd, dim=2)
      write (*,'(3F15.10)') derivsAdd(:,iAt)
    end do

    deallocate(energyDerivsAtom)
    deallocate(derivsAdd)

  end subroutine addGradients


end module dftbp_machinelearning
