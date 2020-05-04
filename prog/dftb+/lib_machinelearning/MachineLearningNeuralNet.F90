#:include 'common.fypp'

!> Neural nets

module dftbp_machinelearning_nn
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_machinelearning_sf
  implicit none
  private

  public :: TMLNeuralNet, TMLNeuralNetInp, TMLNeuralNetSpeciesInp, TMLNeuralNetLayerInp, testNeuralNetInp

  type :: TMLNeuralNetLayerInp

    !> Number of neurons
    integer :: nNeuron

    !> Number of neurons in the previous layer
    integer :: nNeuronPrev

    !> Activation function
    !>  = 1: tanh
    !>  = 2: sigmoid
    !>  = 3: linear
    integer :: activation

    !> Array of bias weights
    real(dp), allocatable :: weightBias(:)

    !> Array of weights (index1 -- this layer, index2 -- previous layer)
    real(dp), allocatable :: weight(:,:)

  end type TMLNeuralNetLayerInp

  !> Contains the input data for the neural nets for every element
  type :: TMLNeuralNetSpeciesInp

    !> Number of layers
    integer :: nLayer

    !> Each individual layer
    type(TMLNeuralNetLayerInp), allocatable :: layer(:)

  end type TMLNeuralNetSpeciesInp

  !> Contains the input data for the neural nets for every element
  type :: TMLNeuralNetInp

    !> number of atoms
    integer :: nAt

    !> number of species/elements
    integer :: nSp

    !> number of symmetry functions (input neurons)
    integer :: nSF

    !> Array of individual nets, one for each element
    type(TMLNeuralNetSpeciesInp), allocatable :: species(:)

  end type TMLNeuralNetInp

  type :: TMLNeuralNetLayer

    !> activity of neurons "act(b + sum w val)"
    real(dp), allocatable :: valueAct(:)

    !> arguments to activation function "b + sum w val"
    real(dp), allocatable :: valueArg(:)

    !> derivatives of activity w.r.t. input parameters (-> force)
    !> 3rd dim -- which neuron (nNeuron)
    !> 2nd dim -- derivatives w.r.t coords of which atom (nAt)
    !> 1st dim -- w.r.t which coord, x/y/z (3)
    real(dp), allocatable :: derivAct(:,:,:)

    !> derivatives of activity w.r.t. valueArg
    !> one value for each neuron (dimension nNeuron)
    real(dp), allocatable :: derivArg(:)

  end type TMLNeuralNetLayer

  !> Type for neural net for 1 atom
  type :: TMLNeuralNetAtom

    !> species of that atom
    integer :: iSp

    !> Number of layers
    integer :: nLayer

    !> Each individual layer
    type(TMLNeuralNetLayer), allocatable :: layer(:)

    !> resulting energy
    real(dp) :: energy

    !> resulting derivatives of energy w.r.t. atom coordinates
    real(dp), allocatable :: deriv(:,:)

  contains

    procedure :: init => NeuralNetAtom_init
    procedure :: evaluate => NeuralNetAtom_evaluate
    procedure :: evaluateDerivs => NeuralNetAtom_evaluateDerivs

  end type TMLNeuralNetAtom

  !> The entire neural net
  type :: TMLNeuralNet

    !> number of species/elements
    integer :: nSp

    !> number of atoms
    integer :: nAt

    !> number of symmetry functions (input neurons)
    integer :: nSF

    !> Array of input for the individual nets, one for each element
    type(TMLNeuralNetSpeciesInp), allocatable :: species(:)

    !> Array of individual nets, one for each atom
    type(TMLNeuralNetAtom), allocatable :: atom(:)

  contains
  
    procedure :: init => NeuralNet_init
    procedure :: evaluate => NeuralNet_evaluate
    procedure :: evaluateDerivs => NeuralNet_evaluateDerivs

  end type TMLNeuralNet

contains

  subroutine NeuralNetAtom_evaluate(this, input, sf, energyAtom)

    !> instance
    class(TMLNeuralNetAtom), intent(inout) :: this

    !> input for neural net atom
    type(TMLNeuralNetSpeciesInp), intent(in) :: input

    !> values of symmetry functions / neuron of 1st layer
    real(dp), intent(in) :: sf(:)

    !> output energy
    real(dp), intent(out) :: energyAtom

    integer :: iLayer, iNeuronThis, iNeuronPrev

    @:ASSERT(input%nLayer == this%nLayer)
    @:ASSERT(size(sf) == size(this%layer(1)%valueAct))

    ! set the 1st layer activation to the values of symmetry functions
    this%layer(0)%valueAct = sf

    do iLayer = 1, this%nLayer
      do iNeuronThis = 1, size(this%layer(iLayer)%valueAct)
        this%layer(iLayer)%valueArg(iNeuronThis) = input%layer(iLayer)%weightBias(iNeuronThis)
        do iNeuronPrev = 1, size(this%layer(iLayer-1)%valueAct)
          this%layer(iLayer)%valueArg(iNeuronThis) = this%layer(iLayer)%valueArg(iNeuronThis) + &
              & input%layer(iLayer)%weight(iNeuronThis, iNeuronPrev) * this%layer(iLayer-1)%valueAct(iNeuronPrev)
        end do
        select case (input%layer(iLayer)%activation)
          case (1)
            this%layer(iLayer)%valueAct(iNeuronThis) = tanh(this%layer(iLayer)%valueArg(iNeuronThis))
          case (2)
            this%layer(iLayer)%valueAct(iNeuronThis) = sigmoid(this%layer(iLayer)%valueArg(iNeuronThis))
          case (3)
            ! no explicit linear activation function
            this%layer(iLayer)%valueAct(iNeuronThis) = this%layer(iLayer)%valueArg(iNeuronThis)
        end select
      end do
    end do

    ! there is only 1 neuron in the last layer -- the energy
    energyAtom = this%layer(this%nLayer)%valueAct(1)

  end subroutine NeuralNetAtom_evaluate


  subroutine NeuralNet_evaluate(this, sf, energyAtom)

    !> instance
    class(TMLNeuralNet), intent(inout) :: this

    !> values of symmetry functions
    real(dp), intent(in) :: sf(:,:)

    !> array of energies per atom
    real(dp), intent(out) :: energyAtom(:)

    integer :: iAt

    write (*,*) "  NEURAL NET EVALUATE"

    @:ASSERT(size(sf, dim=1) == this%input%nSF)
    @:ASSERT(size(sf, dim=2) == this%nAt)
    @:ASSERT(size(energyAtom) == this%nAt)

    do iAt = 1, this%nAt
      call this%atom(iAt)%evaluate(this%species(this%atom(iAt)%iSp), sf(:,iAt), energyAtom(iAt))
    end do

  end subroutine NeuralNet_evaluate


  subroutine NeuralNetAtom_evaluateDerivs(this, input, dsfdr, derivAtom)

    !> instance
    class(TMLNeuralNetAtom), intent(inout) :: this

    !> input for neural net atom
    type(TMLNeuralNetSpeciesInp), intent(in) :: input

    !> derivaties of symmetry functions / activities of neurons in 1st layer
    !> 3rd dim -- which symmetry function, iSf (nSF)
    !> 2nd dim -- derivatives w.r.t coords of which atom (nAt)
    !> 1st dim -- w.r.t which coord, x/y/z (3)
    real(dp), intent(in) :: dsfdr(:,:,:)

    !> output -- derivatives
    !> dim 1 -- xyz (3)
    !> dim 2 -- deriv w.r.t. coords of which atom (nAt)
    real(dp), intent(out) :: derivAtom(:,:)

    integer :: iLayer, iNeuronThis, iNeuronPrev

    @:ASSERT(size(derivAtom) == size(this%layer(1)%derivAct))

    ! set the 1st layer activation derivatives to the derivaties of symmetry functions
    this%layer(0)%derivAct = dsfdr

    do iLayer = 1, this%nLayer
      this%layer(iLayer)%derivAct = 0._dp
      do iNeuronThis = 1, size(this%layer(iLayer)%valueAct)
        do iNeuronPrev = 1, size(this%layer(iLayer-1)%valueAct)
          this%layer(iLayer)%derivAct(:,:,iNeuronThis) = this%layer(iLayer)%derivAct(:,:,iNeuronThis) + &
              & input%layer(iLayer)%weight(iNeuronThis, iNeuronPrev) * this%layer(iLayer-1)%derivAct(:,:,iNeuronPrev)
        end do
        select case (input%layer(iLayer)%activation)
          case (1)
            this%layer(iLayer)%derivArg(iNeuronThis) = tanhDeriv(this%layer(iLayer)%valueArg(iNeuronThis))
          case (2)
            this%layer(iLayer)%derivArg(iNeuronThis) = sigmoidDeriv(this%layer(iLayer)%valueArg(iNeuronThis))
          case (3)
            ! no explicit linear activation function
            this%layer(iLayer)%derivArg(iNeuronThis) = 1._dp
        end select
        this%layer(iLayer)%derivAct(:,:,iNeuronThis) = this%layer(iLayer)%derivAct(:,:,iNeuronThis) * &
            & this%layer(iLayer)%derivArg(iNeuronThis)
      end do
    end do

    ! there is just 1 output neuron -- energyAtom,
    ! and we want its derivatives
    derivAtom(:,:) = this%layer(this%nLayer)%derivAct(:,:,1)

  end subroutine NeuralNetAtom_evaluateDerivs


  subroutine NeuralNet_evaluateDerivs(this, dsfdr, derivAtom)

    !> instance
    class(TMLNeuralNet), intent(inout) :: this

    !> derivatives of symmetry functions / activities in 1st layer
    !> 4th dim -- neural net of which atom (nAt)
    !> 3rd dim -- which symmetry function, iSf (nSF)
    !> 2nd dim -- derivatives w.r.t coords of which atom (nAt)
    !> 1st dim -- w.r.t which coord, x/y/z (3)
    real(dp), intent(in) :: dsfdr(:,:,:,:)

    !> array of derivative vectors
    !> dim 1 = xyz (3)
    !> dim 2 = deriv w.r.t. coords of this atom (nAt)
    !> dim 3 = deriv from the neural net of this atom (nAt)
    real(dp), intent(out) :: derivAtom(:,:,:)

    integer :: iAt

    @:ASSERT(size(dsfdr, dim=1) == 3)
    @:ASSERT(size(dsfdr, dim=2) == this%nAt)
    @:ASSERT(size(dsfdr, dim=3) == this%nSF)
    @:ASSERT(size(dsfdr, dim=4) == this%nAt)
    @:ASSERT(size(derivAtom, dim=1) == 3)
    @:ASSERT(size(derivAtom, dim=2) == this%nAt)
    @:ASSERT(size(derivAtom, dim=3) == this%nAt)

    do iAt = 1, this%nAt
      call this%atom(iAt)%evaluateDerivs(this%species(this%atom(iAt)%iSp), dsfdr(:,:,:,iAt), derivAtom(:,:,iAt))
    end do

  end subroutine NeuralNet_evaluateDerivs


  subroutine testNeuralNetInp(nn, sf)

    !> input for the neural network 
    type(TMLNeuralNetInp), intent(in), target :: nn

    !> input for the symmetry functions
    type(TMLSymmetryFunctionsInp), intent(in) :: sf

    type(TMLNeuralNetSpeciesInp), pointer :: net

    integer :: i, iSf, iSp, iLayer
    character(len = 10) :: act_fn

    write (*,*) "testNeuralNetInp"

    write (*,*) "NumberOfSymmetryFunctions is ", sf%nSymmetryFunctions

    write (*,*) "RadialCutoff = ", sf%radialCutoff

    write (*,*) "RadialParameters"
    do iSf = 1, sf%nRadialFunction
      write (*,*) sf%radialParameters(:,iSf)
    end do

    write (*,*) "AngularCutoff = ", sf%angularCutoff

    write (*,*) "AngularParameters"
    do iSf = 1, sf%nAngularFunction
      write (*,*) sf%angularParameters(:,iSf)
    end do

    write (*,*)
    write (*,*) "Neural nets for each species / element"
    write (*,*)

    do iSp = 1, nn%nSp
      net => nn%species(iSp)
      write (*,*) "iSp = ", iSp, ", # of layers = ", net%nLayer

      do iLayer = 1, net%nLayer
        select case(net%layer(iLayer)%activation)
        case (1)
           act_fn = "tanh"
        case (2)
           act_fn = "sigmoid"
        case (3)
           act_fn = "linear"
        end select
        write (*,*) "layer no. ", iLayer, ", # of neurons = ", net%layer(iLayer)%nNeuron, &
            & ", activation fn = ", act_fn

        if (net%layer(iLayer)%nNeuron > 1) then
          write (*,"(A,3F12.7)") "weightBias 1,3,5", ((net%layer(iLayer)%weightBias(i)), i=1,5,2)
          write (*,"(A,3F12.7)") "weight (1,3),(3,5),(7,10)", &
              & net%layer(iLayer)%weight(1,3), &
              & net%layer(iLayer)%weight(3,5), &
              & net%layer(iLayer)%weight(7,10)
        else
          write (*,"(A,F12.7)") "weightBias 1 ", net%layer(iLayer)%weightBias(1)
          write (*,"(A,F12.7)") "weight (1,3) ", net%layer(iLayer)%weight(1,3)
        end if
      end do
    end do

  end subroutine testNeuralNetInp


  subroutine NeuralNetAtom_init(this, iSp, netInput, nAt, nSF)

  !> instance
  class(TMLNeuralNetAtom), intent(inout) :: this

  !> which species/element
  integer, intent(in) :: iSp

  !> number of layers
  type(TMLNeuralNetSpeciesInp), intent(in) :: netInput

  !> number of atoms in the system
  integer, intent(in) :: nAt

  !> number of symmetry functions
  integer, intent(in) :: nSF

  integer :: iLayer

  this%iSp = iSp
  this%nLayer = netInput%nLayer
  allocate(this%layer(0:this%nLayer))

  ! input layer (iLayer == 0) only has these two arrays
  allocate(this%layer(0)%valueAct(nSF))
  allocate(this%layer(0)%derivAct(3, nAt, nSF))

  ! hidden and output layers have the full set of arrays
  do iLayer = 1, this%nLayer
    allocate(this%layer(iLayer)%valueAct(netInput%layer(iLayer)%nNeuron))
    allocate(this%layer(iLayer)%derivAct(3, nAt, netInput%layer(iLayer)%nNeuron))
    allocate(this%layer(iLayer)%valueArg(netInput%layer(iLayer)%nNeuron))
    allocate(this%layer(iLayer)%derivArg(netInput%layer(iLayer)%nNeuron))
  end do
  allocate(this%deriv(3, nAt))

  end subroutine NeuralNetAtom_init


  subroutine NeuralNet_init(this, input, species)

  !> instance
  class(TMLNeuralNet), target, intent(inout) :: this

  !> input structure
  type(TMLNeuralNetInp), intent(in) :: input

  !> species/element of each atom, dimension (nAt)
  integer, intent(in) :: species(:)

  integer :: iAt, iSp, iLayer

  !> shorthand
  type(TMLNeuralNetLayerInp), pointer :: currentLayer

  write (*,*) "  NEURAL NET INIT"

  @:ASSERT(dim(species) == input%nSp)

  this%nAt = input%nAt
  this%nSp = input%nSp
  this%nSF = input%nSF

  write (*,*) "nAt, nSp, nSF", this%nAt, this%nSp, this%nSF 

  ! the input structure "species"
  allocate(this%species(this%nSp))
  do iSp = 1, this%nSp
    this%species(iSp)%nLayer = input%species(iSp)%nLayer
    allocate(this%species(iSp)%layer(this%species(iSp)%nLayer))

    do iLayer = 1, this%species(iSp)%nLayer
      currentLayer => this%species(iSp)%layer(iLayer)

      currentLayer%nNeuron = input%species(iSp)%layer(iLayer)%nNeuron
      currentLayer%nNeuronPrev = input%species(iSp)%layer(iLayer)%nNeuronPrev
      currentLayer%activation = input%species(iSp)%layer(iLayer)%activation

      allocate(currentLayer%weightBias(currentLayer%nNeuron))
      allocate(currentLayer%weight(currentLayer%nNeuron, currentLayer%nNeuronPrev))
      currentLayer%weightBias = input%species(iSp)%layer(iLayer)%weightBias
      currentLayer%weight = input%species(iSp)%layer(iLayer)%weight
    end do
  end do

  ! the structure "atom"
  allocate(this%atom(this%nAt))

  do iAt = 1, this%nAt
    call this%atom(iAt)%init(species(iAt), this%species(species(iAt)), this%nAt, input%nSF)
  end do

  end subroutine NeuralNet_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!            ARITHMETICS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure function sigmoid(x)

  real(dp), intent(in) :: x
  real(dp) :: sigmoid
   
  sigmoid = 1._dp / (1._dp + exp(-x))

  end function sigmoid

  pure function sigmoidDeriv(x)

  real(dp), intent(in) :: x
  real(dp) :: sigmoidDeriv
  real(dp) :: sigmoid
   
  sigmoid = 1._dp / (1._dp + exp(-x))
  sigmoidDeriv = sigmoid * (1._dp - sigmoid)

  end function sigmoidDeriv

! pure function tanh -- is Fortran intrinsic

  pure function tanhDeriv(x)

  real(dp), intent(in) :: x
  real(dp) :: tanhDeriv

  tanhDeriv = 1._dp - tanh(x)**2

  end function tanhDeriv

end module dftbp_machinelearning_nn
