#:include 'common.fypp'

!> Offers everything which is publicly available when dealing with machine learning.

module dftbp_machinelearning
  use dftbp_common_assert
  use dftbp_common_accuracy
  use dftbp_common_constants
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

    procedure :: init => TMachineLearning_init
    procedure :: getEnergy => TMachineLearning_getEnergy
    procedure :: addGradients => TMachineLearning_addGradients

  end type TMachineLearning

contains

  !> Initialize with data from input
  subroutine TMachineLearning_init(this, input)

    !> instance
    class(TMachineLearning), intent(inout) :: this

    !> the input structure to be linked into here
    type(TMachineLearningInp), intent(in), target :: input

    write (*,*) "MACHINE LEARNING INIT"

    @:ASSERT(size(input%sf%species) == input%sf%nAtom)

    this%input => input

    call this%sf%init(input%sf)

    call this%nn%init(input%nn, input%sf%species)

    ! store this information:
    ! atomic indices, atomic element/species, and symmetry function hyperparameters

  end subroutine TMachineLearning_init

  !> Get energy contributions from machine learning
  function TMachineLearning_getEnergy(this, coordsIn, img2CentCell) result(energy)

    !> instance
    class(TMachineLearning), intent(inout) :: this

    !> Current coordinates
    real(dp), intent(in) :: coordsIn(:,:)

    !> Updated mapping to central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Resulting energy contribution
    real(dp) :: energy

    !> Energy per atom -- TODO output that as well?
    real(dp), allocatable :: energyAtom(:)

    integer :: iAt

    ! Coordinates of atoms to be treated with machine learning
    real(dp), allocatable :: coords(:,:)

    ! Reduce the coordinates to only those atoms which are treated with machine learning
    allocate(coords(3, size(this%sf%indAtomsML)))
    do iAt = 1, size(this%sf%indAtomsML)
      coords(:, iAt) = coordsIn(:, this%sf%indAtomsML(iAt))
    end do

    write (*,*) "MACHINE LEARNING ENERGY"

    allocate(energyAtom(size(coords, dim=2)))

    ! calculate all of the symmetry functions
    call this%sf%prepare(coords)
    deallocate(coords)
    call this%sf%evaluate()

!   write (*,*) "  Symmetry functions evaluated:"
!   do iAt = 1, this%sf%nSF
!     write (*,'(A,I4,A,30F10.6)') "    SF ", iAt, ": ", this%sf%sf(iAt,:)
!   end do
!   write (*,*) "  Symmetry functions -- END"

    ! feed those values into the neural net
    call this%nn%evaluate(this%sf%sf, energyAtom)

    ! if the unit is kcal/mol throughout the machine learning code
    if (this%nn%tUnitKcalMol) then
      ! then convert to Hartree
      energyAtom(:) = energyAtom(:) * kcal_mol__Hartree
    end if

    energy = sum(energyAtom)

    ! scale the energy if requested
    if (this%nn%tScaleEnergy) then
      energy = energy * this%nn%scaleFactors(2) + this%nn%scaleFactors(1)
    end if

    write (*,'(A,F15.10)') "MACHINE_LEARNING_ENERGY ", energy

!   do iAt = 1, size(energyAtom)
!     write (*,'(F15.10)') energyAtom(iAt)
!   end do

    deallocate(energyAtom)

  end function TMachineLearning_getEnergy


  !> Gradient contribution from machine learning
  subroutine TMachineLearning_addGradients(this, derivs, img2CentCell)

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

    ! if the unit is kcal/mol throughout the machine learning code
    if (this%nn%tUnitKcalMol) then
      ! then convert gradients to Hartree
      energyDerivsAtom(:,:,:) = energyDerivsAtom(:,:,:) * kcal_mol__Hartree
    end if

    ! if the unit is Angstrom throughout the machine learning code
    if (this%sf%tUnitAngstrom) then
      ! then convert 1/AA to 1/bohr
      energyDerivsAtom(:,:,:) = energyDerivsAtom(:,:,:) * Bohr__AA
    end if

    ! sum contributions from all of the neural nets (one for each atom)
    derivsAdd = sum(energyDerivsAtom, 3)

    ! scale the gradients if requested
    if (this%nn%tScaleEnergy) then
      derivsAdd = derivsAdd * this%nn%scaleFactors(2)
    end if

    ! Add the contributions to the total gradients
    !derivs = derivs + derivsAdd

    ! Add the contributions to the correct atoms
    !   (those which are treated with machine learning):
    do iAt = 1, size(derivsAdd, dim=2)
      derivs(:, this%sf%indAtomsML(iAt)) = derivs(:, this%sf%indAtomsML(iAt)) + derivsAdd(:, iAt)
      write (*,'(3F15.10)') derivsAdd(:,iAt)
    end do

    deallocate(energyDerivsAtom)
    deallocate(derivsAdd)

  end subroutine TMachineLearning_addGradients


end module dftbp_machinelearning
