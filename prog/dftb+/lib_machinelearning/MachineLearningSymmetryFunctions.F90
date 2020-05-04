#:include 'common.fypp'

!> Symmetry functions for machine learning

module dftbp_machinelearning_sf
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_constants
  implicit none
  private

  public :: TMLSymmetryFunctionsInp, TMLSymmetryFunctions


  !> Contains the initialisation data for the Slater-Kirkwood module.
  type :: TMLSymmetryFunctionsInp

    !> Number of symmetry functions (for each species/element)
    integer :: nSymmetryFunctions

    !> Atom species by atom, ordered by atomic number (dimension nAt)
    integer, allocatable :: speciesOrder(:)

    !> Radial functions -- cut-off
    real(dp) :: radialCutoff

    !> Number of radial functions in the set
    integer :: nRadialFunction

    !> Radial functions -- parameters Rs and eta (2, nRadialFunction)
    real(dp), allocatable :: radialParameters(:,:)

    !> Angular functions -- cut-off
    real(dp) :: angularCutoff

    !> Number of angular functions in the set
    integer :: nAngularFunction

    !> Angular functions -- parameters eta, zeta and lambda (3, nAngularFunction)
    real(dp), allocatable :: angularParameters(:,:)

    !> Use neighborsearching to evaluate symmetry functions or not
    logical :: tNeighborSearching

  end type TMLSymmetryFunctionsInp


  !> Contains the data for the symmetry function module
  type :: TMLSymmetryFunctions

    !> Number of atoms
    integer :: nAt

    !> Number of species
    integer :: nSp

!   !> Number of species pairs
!   integer :: nSpPair

    !> Use neighborsearching to evaluate symmetry functions or not
    logical :: tNeighborSearching

    !> Number of symmetry functions
    integer :: nSF

    !> Radial functions -- cut-off (COPY FROM INPUT)
    real(dp) :: radialCutoff

    !> Number of radial functions in the set (COPY FROM INPUT)
    integer :: nRadialFunction

    !> Radial functions -- parameters Rs and eta (2, nRadialFunction) (COPY FROM INPUT)
    real(dp), allocatable :: radialParameters(:,:)

    !> Angular functions -- cut-off (COPY FROM INPUT)
    real(dp) :: angularCutoff

    !> Number of angular functions in the set (COPY FROM INPUT)
    integer :: nAngularFunction

    !> Angular functions -- parameters eta, zeta and lambda (3, nAngularFunction) (COPY FROM INPUT)
    real(dp), allocatable :: angularParameters(:,:)

    !> array of symmetry function values
    !> 1st dim -- indiv. symmetry functions (nSymmetryFunctions), 2nd dim -- atoms (nAt)
    real(dp), allocatable :: sf(:,:)

    !> derivatives of symmetry functions w.r.t. atom coordinates
    !> 4th -- atoms (nAt), 3rd -- iSf (nSF), 2nd -- w.r.t which atom (nAt), 1st -- w.r.t x/y/z (3)
    real(dp), allocatable :: dsfdr(:,:,:,:)

    !> Atom species by atom, ordered by atomic number (dimension nAt)
    integer, allocatable :: speciesOrder(:)

    !> Atom coordinates (3 x nAt matrix)
    real(dp), allocatable :: coords(:,:)

    !> Interatomic distances (symmetric nAt x nAt matrix)
    real(dp), allocatable :: distance(:,:)

    !> Neighborlist, 1st dim -- of which atom, 2nd dim -- # of neighbor
    integer, allocatable :: neighborListArr(:,:)

    !> Number of neighbors in the list, 1st dim -- of which atom
    integer, allocatable :: neighborListCount(:)

    !> Neighbor pair list, 1st dim -- of which atom, 2nd dim -- # of neighbor pair, 3rd -- 2 (2 atoms)
    integer, allocatable :: neighborPairArr(:,:,:)

    !> Number of neighbor pairs in the list, 1st dim -- of which atom
    integer, allocatable :: neighborPairCount(:)

  contains

    procedure :: init => SymmetryFunctions_init
    procedure :: calculateDistance => SymmetryFunctions_calculateDistance
    procedure :: getNeighborsReally => SymmetryFunctions_getNeighborsReally
    procedure :: getNeighborsAll => SymmetryFunctions_getNeighborsAll
    procedure :: prepare => SymmetryFunctions_prepare
    procedure :: getSpeciesPair => SymmetryFunctions_getSpeciesPair
    procedure :: evaluate => SymmetryFunctions_evaluate
    procedure :: evaluateDerivs => SymmetryFunctions_evaluateDerivs

  end type TMLSymmetryFunctions


contains


  !> Do this once, at the start of the DFTB+ run
  subroutine SymmetryFunctions_init(this, input, nAt, nSp)

  !> instance
  class(TMLSymmetryFunctions), intent(inout) :: this

  !> input structure
  type(TMLSymmetryFunctionsInp), intent(in) :: input

  !> number of atoms
  integer, intent(in) :: nAt

  !> number of species
  integer, intent(in) :: nSp

  write (*,*) "  SYMMETRY FUNCTIONS INIT"

  this%nAt = nAt
  this%nSp = nSp
! this%nSpPair = nSp * (nSp + 1) / 2
  this%tNeighborSearching = input%tNeighborSearching
  this%nSF = input%nSymmetryFunctions

  this%nRadialFunction = input%nRadialFunction
  this%radialCutoff = input%radialCutoff
  this%nAngularFunction = input%nAngularFunction
  this%angularCutoff = input%angularCutoff

  allocate(this%radialParameters(2, this%nRadialFunction))
  this%radialParameters = input%radialParameters
  allocate(this%angularParameters(3, this%nAngularFunction))
  this%angularParameters = input%angularParameters

  allocate(this%speciesOrder(nAt))
  this%speciesOrder = input%speciesOrder

  allocate(this%coords(3, nAt))
  allocate(this%distance(nAt, nAt))
  allocate(this%neighborListArr(nAt, nAt-1))
  allocate(this%neighborListCount(nAt))
  allocate(this%neighborPairArr(nAt, (nAt-1)*(nAt-2)/2, 2))
  allocate(this%neighborPairCount(nAt))

  allocate(this%sf(this%nSF, this%nAt))
  allocate(this%dsfdr(3, this%nAt, this%nSF, this%nAt))

  end subroutine SymmetryFunctions_init


  !> get interatomic distances
  subroutine SymmetryFunctions_calculateDistance(this)

  !> instance
  class(TMLSymmetryFunctions), intent(inout) :: this

  ! indices to atoms
  integer :: iAt1, iAt2

  do iAt1 = 1, this%nAt
    this%distance(iAt1, iAt1) = 0._dp
    do iAt2 = iAt1+1, this%nAt
      this%distance(iAt1, iAt2) = norm2(this%coords(:,iAt1) - this%coords(:,iAt2))
      this%distance(iAt2, iAt1) = this%distance(iAt1, iAt2)
    end do
  end do

  end subroutine SymmetryFunctions_calculateDistance


  !> create the list of neighbors and neighbor pairs for each atom
  subroutine SymmetryFunctions_getNeighborsReally(this, radialCutoff, angularCutoff)

  !> instance
  class(TMLSymmetryFunctions), intent(inout) :: this

  !> radial cutoff
  real(dp), intent(in) :: radialCutoff

  !> angular cutoff
  real(dp), intent(in) :: angularCutoff

  ! indices to atoms
  integer :: iAt1, iAt2, iAt3

  ! counter in arrays
  integer :: counterNeig, counterPair

  do iAt1 = 1, this%nAt
    counterNeig = 0
    counterPair = 0
    do iAt2 = 1, this%nAt
      ! neighbor list
      if ((iAt1 /= iAt2) .and. (this%distance(iAt1, iAt2) < radialCutoff)) then
        counterNeig = counterNeig + 1
        this%neighborListArr(iAt1, counterNeig) = iAt2
      end if
      ! neighbor pair list
      if ((iAt1 /= iAt2) .and. (this%distance(iAt1, iAt2) < angularCutoff)) then
        do iAt3 = iAt2 + 1, this%nAt
          if ((iAt1 /= iAt3) .and. (this%distance(iAt1, iAt3) < angularCutoff)) then
            counterPair = counterPair + 1
            this%neighborPairArr(iAt1, counterPair, 1) = iAt2
            this%neighborPairArr(iAt1, counterPair, 2) = iAt3
          end if
        end do
      end if
      ! done for iAt2
    end do
    this%neighborListCount(iAt1) = counterNeig
    this%neighborPairCount(iAt1) = counterPair
  end do

! ! test -- print the lists
! do iAt1 = 1, this%nAt
!   write (*,'(A,I3,A)',advance="no") "neighbors to atom ", iAt1, ":"
!   write (*,'(10I3)') (this%neighborListArr(iAt1, counterNeig), counterNeig=1,this%neighborListCount(iAt1))
!   write (*,'(A,I3,A,I3,A)',advance="no") "neighbor pairs to atom ", iAt1, "(there are ", this%neighborPairCount(iAt1), " :"
!   do counterPair = 1, this%neighborPairCount(iAt1)
!     write (*,'(T10,2I3)') this%neighborPairArr(iAt1, counterPair, 1), this%neighborPairArr(iAt1, counterPair, 2)
!   end do
! end do

  end subroutine SymmetryFunctions_getNeighborsReally

  !> create the list of neighbors and neighbor pairs for each atom
  !> this version: rudimentary, contains all atoms, not only neighbors!
  subroutine SymmetryFunctions_getNeighborsAll(this)

  !> instance
  class(TMLSymmetryFunctions), intent(inout) :: this

  ! indices to atoms
  integer :: iAt1, iAt2, iAt3

  ! counter in arrays
  integer :: counter

  ! neighbor list
  do iAt1 = 1, this%nAt
    this%neighborListCount(iAt1) = this%nAt - 1
    counter = 0
    do iAt2 = 1, this%nAt
      if (iAt1 /= iAt2) then
         counter = counter + 1
         this%neighborListArr(iAt1, counter) = iAt2
      end if
    end do
  end do

  ! neighbor pair list
  do iAt1 = 1, this%nAt
    this%neighborPairCount(iAt1) = (this%nAt - 1) * (this%nAt - 2) / 2
    counter = 0
    do iAt2 = 1, this%nAt
      if (iAt1 /= iAt2) then
        do iAt3 = iAt2 + 1, this%nAt
          if (iAt1 /= iAt3) then
            counter = counter + 1
            this%neighborPairArr(iAt1, counter, 1) = iAt2
            this%neighborPairArr(iAt1, counter, 2) = iAt3
          end if
        end do
      end if
    end do
  end do

  end subroutine SymmetryFunctions_getNeighborsAll


  !> Do this at the beginning of each calculation of energy
  subroutine SymmetryFunctions_prepare(this, coords)

  !> instance
  class(TMLSymmetryFunctions), intent(inout) :: this

  !> Current coordinates
  real(dp), intent(in) :: coords(:,:)

  @:ASSERT(all(shape(coords) = shape(this%coords)))

  ! conversion is OK
  this%coords = coords * Bohr__AA
  
  ! Calculate interatomic distances
  call this%calculateDistance()

  ! Create the neighborlist
  if (this%tNeighborSearching) then
    call this%getNeighborsReally(this%radialCutoff, this%angularCutoff)
  else
    call this%getNeighborsAll()
  end if
  
  end subroutine SymmetryFunctions_prepare


  pure function SymmetryFunctions_getSpeciesPair(this, iAt1, iAt2)

  !> instance
  class(TMLSymmetryFunctions), intent(in) :: this

  !> atom indices
  integer, intent(in) :: iAt1, iAt2

  !> output
  integer :: SymmetryFunctions_getSpeciesPair

  ! species of atoms 1 and 2; the smaller and the larger species number
  integer :: iSp1, iSp2, iSpSm, iSpLa

  iSp1 = this%speciesOrder(iAt1)
  iSp2 = this%speciesOrder(iAt2)

  if (iSp1 < iSp2) then
    iSpSm = iSp1
    iSpLa = iSp2
  else
    iSpSm = iSp2
    iSpLa = iSp1
  end if

  ! expression that is unique and always > 0 and <= nSp * (nSp + 1) / 2
  SymmetryFunctions_getSpeciesPair = this%nSp * (iSpSm - 1) - iSpSm * (iSpSm - 1) / 2 + iSpLa

  end function SymmetryFunctions_getSpeciesPair


  subroutine SymmetryFunctions_evaluate(this)

    !> instance
    class(TMLSymmetryFunctions), target, intent(inout) :: this

    integer :: iAt1, iAt2, iAt3, iSp2, iNeig, iSF, speciesPair, iSymmFuncInd
    real(dp) :: R12, R13, R23, Rs, eta, zeta, lambda

    this%sf = 0._dp

    do iAt1 = 1, this%nAt

      ! radial symmetry functions
      do iNeig = 1, this%neighborListCount(iAt1)
        ! index of the neighbor atom
        iAt2 = this%neighborListArr(iAt1, iNeig)
        iSp2 = this%speciesOrder(iAt2)
        R12 = this%distance(iAt1, iAt2)
        do iSF = 1, this%nRadialFunction
          Rs = this%radialParameters(1, iSF)
          eta = this%radialParameters(2, iSF)
          iSymmFuncInd = (iSp2 - 1) * this%nRadialFunction + iSF
          this%sf(iSymmFuncInd, iAt1) = this%sf(iSymmFuncInd, iAt1) + radialFilter(Rs, eta, R12)
        end do
      end do
      
      ! angular symmetry functions
      do iNeig = 1, this%neighborPairCount(iAt1)
        ! index of the neighbor atoms
        iAt2 = this%neighborPairArr(iAt1, iNeig, 1)
        iAt3 = this%neighborPairArr(iAt1, iNeig, 2)
        speciesPair = this%getSpeciesPair(iAt2, iAt3)
        R12 = this%distance(iAt1, iAt2)
        R13 = this%distance(iAt1, iAt3)
        R23 = this%distance(iAt2, iAt3)
        do iSF = 1, this%nAngularFunction
          eta = this%angularParameters(1, iSF)
          zeta = this%angularParameters(2, iSF)
          lambda = this%angularParameters(3, iSF)
          iSymmFuncInd = this%nSp * this%nRadialFunction &
              & + (speciesPair - 1) * this%nAngularFunction + iSF
          this%sf(iSymmFuncInd, iAt1) = this%sf(iSymmFuncInd, iAt1) &
              & + angularFilter(eta, zeta, lambda, R12, R13, R23)
        end do
      end do

    end do
        
  end subroutine SymmetryFunctions_evaluate


  subroutine SymmetryFunctions_evaluateDerivs(this)

    !> instance
    class(TMLSymmetryFunctions), target, intent(inout) :: this

    integer :: iAt0, iAt1, iAt2, iAt3, iSp2, iNeig, iSF, iSpPair23, iSymmFuncInd
    real(dp) :: R12, R13, R23, Rs, eta, zeta, lambda

    real(dp), dimension(3) :: derivAdd, xyz1, xyz2, xyz3
    logical :: tAtom0Is1, tAtom0Is2, tAtom0Is3

    this%dsfdr = 0._dp

    do iAt1 = 1, this%nAt

      xyz1 = this%coords(:,iAt1)

      ! derivatives w.r.t. coordinates of atom iAt0
      do iAt0 = 1, this%nAt

        ! are we differentiating w.r.t. iAt1?
        if (iAt0 == iAt1) then
          tAtom0Is1 = .true.
        else
          tAtom0Is1 = .false.
        end if

        ! radial symmetry functions
        do iNeig = 1, this%neighborListCount(iAt1)
          ! index of the neighbor atom
          iAt2 = this%neighborListArr(iAt1, iNeig)

          ! are we differentiating w.r.t. iAt2?
          if (iAt0 == iAt2) then
            tAtom0Is2 = .true.
          else
            tAtom0Is2 = .false.
          end if

          ! if neither iAt1 nor iAt2 is iAt0, then there is no contribution
          if (.not. (tAtom0Is1 .or. tAtom0Is2)) then
            cycle
          end if

          iSp2 = this%speciesOrder(iAt2)
          R12 = this%distance(iAt1, iAt2)
          xyz2 = this%coords(:,iAt2)

          do iSF = 1, this%nRadialFunction
            Rs = this%radialParameters(1, iSF)
            eta = this%radialParameters(2, iSF)
            iSymmFuncInd = (iSp2 - 1) * this%nRadialFunction + iSF
            derivAdd = radialFilterDeriv(Rs, eta, R12, xyz1, xyz2, tAtom0Is1, tAtom0Is2)
            this%dsfdr(:, iAt0, iSymmFuncInd, iAt1) = this%dsfdr(:, iAt0, iSymmFuncInd, iAt1) + derivAdd
          end do
        end do
        
        ! angular symmetry functions
        do iNeig = 1, this%neighborPairCount(iAt1)
          ! index of the neighbor atoms
          iAt2 = this%neighborPairArr(iAt1, iNeig, 1)
          iAt3 = this%neighborPairArr(iAt1, iNeig, 2)

          ! are we differentiating w.r.t. iAt2 or iAt3?
          if (iAt0 == iAt2) then
            tAtom0Is2 = .true.
          else
            tAtom0Is2 = .false.
          end if
          if (iAt0 == iAt3) then
            tAtom0Is3 = .true.
          else
            tAtom0Is3 = .false.
          end if

          ! if neither iAt1 nor iAt2 nor iAt3 is iAt0, then there is no contribution
          if (.not. (tAtom0Is1 .or. tAtom0Is2 .or. tAtom0Is3)) then
            cycle
          end if

          iSpPair23 = this%getSpeciesPair(iAt2, iAt3)
          R12 = this%distance(iAt1, iAt2)
          R13 = this%distance(iAt1, iAt3)
          R23 = this%distance(iAt2, iAt3)
          xyz2 = this%coords(:,iAt2)
          xyz3 = this%coords(:,iAt3)

          do iSF = 1, this%nAngularFunction
            eta = this%angularParameters(1, iSF)
            zeta = this%angularParameters(2, iSF)
            lambda = this%angularParameters(3, iSF)
            iSymmFuncInd = this%nSp * this%nRadialFunction &
                & + (iSpPair23 - 1) * this%nAngularFunction + iSF
            derivAdd = angularFilterDeriv(eta, zeta, lambda, R12, R13, R23, xyz1, xyz2, xyz3, &
                & tAtom0Is1, tAtom0Is2, tAtom0Is3)
            this%dsfdr(:, iAt0, iSymmFuncInd, iAt1) = this%dsfdr(:, iAt0, iSymmFuncInd, iAt1) + derivAdd
          end do
        end do

      end do ! iAt0

    end do ! iAt
        
  end subroutine SymmetryFunctions_evaluateDerivs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!                   ARITHMETICS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Radial filter for symmetry functions
  pure function radialFilter(Rs, eta, Rij)

    !> radial symmetry function parameters
    real(dp), intent(in) :: Rs, eta

    !> distance values between two given atoms i and j;
    real(dp), intent(in) :: Rij

    !> output
    real(dp) :: radialFilter

    radialFilter = exp(-eta * (Rij - Rs)**2)

  end function radialFilter


  !> Derivative of the radial symmetry function
  pure function radialFilterDeriv(Rs, eta, R12, xyz1, xyz2, tAtom0Is1, tAtom0Is2)

    !> radial symmetry function parameters
    real(dp), intent(in) :: Rs, eta

    !> distance 
    real(dp), intent(in) :: R12

    !> coordinates of atoms 1 and 2
    real(dp), dimension(3), intent(in) :: xyz1, xyz2

    !> are we differentiating w.r.t. atom 1 or 2? (never both)
    logical, intent(in) :: tAtom0Is1, tAtom0Is2

    !> output
    real(dp), dimension(3) :: radialFilterDeriv

    real(dp) :: dG_dR12
  ! real(dp), dimension(3) :: vector

    @:ASSERT(tAtom0Is1 .or. tAtom0Is2)

    dG_dR12 = - 2._dp * exp(- eta * (R12 - Rs)**2) * eta * (R12 - Rs)

    if (tAtom0Is1) then
      radialFilterDeriv = dG_dR12 * (xyz1(:) - xyz2(:)) / R12
    end if

    if (tAtom0Is2) then
      radialFilterDeriv = dG_dR12 * (xyz2(:) - xyz1(:)) / R12
    end if

  end function radialFilterDeriv


  pure function angularFilter(eta, zeta, lambda, Rij, Rik, Rjk)

    !> angular symmetry function parameters
    real(dp), intent(in) :: eta, zeta, lambda

    !> distances among three atoms i, j, k
    real(dp), intent(in) :: Rij, Rik, Rjk

    !> output
    real(dp) :: angularFilter

    real(dp) :: cosAngle, radFilter

    cosAngle = (Rij**2 + Rik**2 - Rjk**2)/(2._dp * Rij * Rik)
    radFilter = exp(-eta * (Rij + Rik + Rjk)**2)
    angularFilter = 2._dp**(1._dp - zeta) * (1._dp + lambda * cosAngle)**zeta * radFilter

  end function angularFilter


  !> Derivative of the angular symmetry function
  pure function angularFilterDeriv(eta, zeta, lambda, R12, R13, R23, xyz1, xyz2, xyz3, &
      & tAtom0Is1, tAtom0Is2, tAtom0Is3)

    !> angular symmetry function parameters
    real(dp), intent(in) :: eta, zeta, lambda

    !> distances
    real(dp), intent(in) :: R12, R13, R23

    !> coordinates of atoms 1, 2 and 3
    real(dp), dimension(3), intent(in) :: xyz1, xyz2, xyz3

    !> are we differentiating w.r.t. atom 1 or 2 or 3? (never two or all)
    logical, intent(in) :: tAtom0Is1, tAtom0Is2, tAtom0Is3

    !> output
    real(dp), dimension(3) :: angularFilterDeriv

    real(dp) :: expZetaMinusOne, cosAngle, onePlusLambdaCosAngle, sumR123, gaussDist
    real(dp) :: dG_dR12, dG_dR13, dG_dR23
  ! real(dp), dimension(3) :: vector12, vector13, vector21, vector23, vector31, vector32

    @:ASSERT(tAtom0Is1 .or. tAtom0Is2 .or. tAtom0Is3)

    cosAngle = (R12**2 + R13**2 - R23**2)/(2._dp * R12 * R13)
    onePlusLambdaCosAngle = 1._dp + lambda * cosAngle
    expZetaMinusOne = (onePlusLambdaCosAngle / 2.)**(zeta - 1._dp)
    sumR123 = R12 + R13 + R23
    gaussDist = exp(- eta * (R12 + R13 + R23)**2)

    if (tAtom0Is1 .or. tAtom0Is2) then
    ! dG_dR12 = - 2._dp**(2._dp - zeta) * gaussDist * eta * sumR123 * &
    !         &   onePlusLambdaCosAngle**zeta &
    !         & + 2._dp**(1._dp - zeta) * gaussDist * lambda * (1._dp / R13 - cosAngle / R12) * &
    !         &   onePlusLambdaCosAngle**(zeta - 1._dp) * zeta
      dG_dR12 = expZetaMinusOne * gaussDist * &
              & ( - 2._dp * eta * sumR123 * onePlusLambdaCosAngle &
              &   + lambda * (1._dp / R13 - cosAngle / R12) * zeta)
    end if

    if (tAtom0Is1 .or. tAtom0Is3) then
    ! dG_dR13 = - 2._dp**(2._dp - zeta) * gaussDist * eta * sumR123 * &
    !         &   onePlusLambdaCosAngle**zeta &
    !         & + 2._dp**(1._dp - zeta) * gaussDist * lambda * (1._dp / R12 - cosAngle / R13) * &
    !         &   onePlusLambdaCosAngle**(zeta - 1._dp) * zeta
      dG_dR13 = expZetaMinusOne * gaussDist * &
              & ( - 2._dp * eta * sumR123 * onePlusLambdaCosAngle &
              &   + lambda * (1._dp / R12 - cosAngle / R13) * zeta)
    end if

    if (tAtom0Is2 .or. tAtom0Is3) then
    ! dG_dR23 = - 2._dp**(2._dp - zeta) * gaussDist * eta * sumR123 * &
    !         &   onePlusLambdaCosAngle**zeta &
    !         & - 2._dp**(1._dp - zeta) * gaussDist * lambda * R23 / R12 / R13 * &
    !         &   onePlusLambdaCosAngle**(zeta - 1._dp) * zeta
      dG_dR23 = - expZetaMinusOne * gaussDist * &
              & ( 2._dp * eta * sumR123 * onePlusLambdaCosAngle &
              & + lambda * R23 / R13 / R12 * zeta)
    end if

    if (tAtom0Is1) then
    ! vector21 = (xyz1(:) - xyz2(:)) / R12
    ! vector31 = (xyz1(:) - xyz3(:)) / R13
      angularFilterDeriv = dG_dR12 * (xyz1(:) - xyz2(:)) / R12 + dG_dR13 * (xyz1(:) - xyz3(:)) / R13
    end if

    if (tAtom0Is2) then
    ! vector12 = (xyz2(:) - xyz1(:)) / R12
    ! vector32 = (xyz2(:) - xyz3(:)) / R23
      angularFilterDeriv = dG_dR12 * (xyz2(:) - xyz1(:)) / R12 + dG_dR23 * (xyz2(:) - xyz3(:)) / R23
    end if

    if (tAtom0Is3) then
    ! vector13 = (xyz3(:) - xyz1(:)) / R13
    ! vector23 = (xyz3(:) - xyz2(:)) / R23
      angularFilterDeriv = dG_dR13 * (xyz3(:) - xyz1(:)) / R13 + dG_dR23 * (xyz3(:) - xyz2(:)) / R23
    end if

  end function angularFilterDeriv


end module dftbp_machinelearning_sf
