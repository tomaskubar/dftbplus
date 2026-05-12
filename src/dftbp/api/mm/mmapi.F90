!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Provides DFTB+ API for MM-type high level access
module dftbp_mmapi
  use, intrinsic :: iso_fortran_env, only : output_unit
  use dftbp_common_accuracy, only : dp
  use dftbp_common_environment, only : TEnvironment, TEnvironment_init
  use dftbp_common_file, only : closeFile, openFile, TFileDescr
  use dftbp_common_globalenv, only : destructGlobalEnv, initGlobalEnv, instanceSafeBuild, withMpi
  use dftbp_common_status, only : TStatus
  use dftbp_dftb_sparse2dense, only : unpackHS
  use dftbp_dftbplus_hsdhelpers, only : doPostParseJobs
  use dftbp_dftbplus_initprogram, only : TDftbPlusMain
  use dftbp_dftbplus_inputdata, only : TInputData
  use dftbp_dftbplus_mainapi, only : checkSpeciesNames, doOneTdStep, finalizeTimeProp,&
      & getAtomicMasses, getCM5Charges, getCutOff, getElStatPotential, getEnergy,&
      & getExtChargeGradients, getGradients, getGrossCharges, getLocalKS, getRefCharges,&
      & getStressTensor, getTdForces, initializeTimeProp, nrOfAtoms, nrOfKPoints, nrOfLocalKS,&
      & nrOfSpin, setExternalCharges, setExternalEfield, setExternalPotential, setGeometry,&
      & setNeighbourList, setQDepExtPotProxy, setRefCharges, setTdCoordsAndVelos,&
      & setTdElectricField, updateDataDependentOnSpeciesOrdering, nrOfOrbitals,&
      & getEigenValues, getEigenVectors, getHamilOverl, getFragmentBasedHamiltonian
  use dftbp_dftbplus_parser, only : parseHsdTree, readHsdFile, rootTag, TParserFlags
  use dftbp_dftbplus_qdepextpotgen, only : TQDepExtPotGen, TQDepExtPotGenWrapper
  use dftbp_dftbplus_qdepextpotproxy, only : TQDepExtPotProxy, TQDepExtPotProxy_init
  use dftbp_extlibs_xmlf90, only : appendChild, createDocumentNode, createElement, destroyNode,&
      & fnode
  use dftbp_fmo, only: TPointersToPhase1, checkInvertPhase
  use dftbp_fmogradient, only: fmoGradients
  use dftbp_io_charmanip, only : newline
  use dftbp_io_hsdutils, only : getChild
  use dftbp_io_message, only : error
  use dftbp_type_linkedlist, only : append, asArray, get, init, len, TListString
  use dftbp_type_typegeometry, only : TGeometry
  implicit none
  private

  public :: TDftbPlus, getDftbPlusBuild, getDftbPlusApi
  public :: TDftbPlus_init, TDftbPlus_destruct
  public :: TDftbPlusAtomList
  public :: TDftbPlusInput, TDftbPlusInput_destruct
  public :: TQDepExtPotGen
  public :: getMaxAngFromSlakoFile, convertAtomTypesToSpecies


  !> List of QM atoms and species for DFTB+ calculation
  type :: TDftbPlusAtomList
    !> Number of atoms
    integer :: nAtom
    !> Linked list of chemical symbols of elements (species names), size=nSpecies
    type(TListString) :: speciesNames
    !> Array of species for each atom, size=nAtom
    integer, allocatable :: species(:)
  contains
    !> Read list of atoms
    procedure :: get => TDftbPlusAtomList_get
    !> Insert the list of atoms into the input data structure
    procedure :: add => TDftbPlusAtomList_addToInpData
  end type TDftbPlusAtomList


  !> Input tree for DFTB+ calculation
  type :: TDftbPlusInput
    !> Tree for HSD format input
    type(fnode), pointer :: hsdTree => null()
  contains
    !> Obtain the root of the tree of input
    procedure :: getRootNode => TDftbPlusInput_getRootNode
    !> Finaliser
    final :: TDftbPlusInput_final
  end type TDftbPlusInput


  !> A DFTB+ calculation
  type :: TDftbPlus
    private
    !> Computational environment
    type(TEnvironment), allocatable :: env
    !> Calculation instance
    type(TDftbPlusMain), allocatable :: main
    !> Number of fragments in the system (for FMO calculations)
    integer :: nSite = 1
    !> Pointers to the data from phase 1 of FMO.
    type(TPointersToPhase1), allocatable :: ptrsPhase1(:)
    !> Has this been initialised and ready to use
    logical :: isInitialised = .false.
  contains
    !> Read input from a file
    procedure :: getInputFromFile => TDftbPlus_getInputFromFile
    !> Get an empty input to populate
    procedure :: getEmptyInput => TDftbPlus_getEmptyInput
    !> Set up a DFTB+ calculator from input tree
    procedure :: setupCalculator => TDftbPlus_setupCalculator
    !> Set/replace the geometry of a calculator
    procedure :: setGeometry => TDftbPlus_setGeometry
    !> Set/replace the neighbour list
    procedure :: setNeighbourList => TDftbPlus_setNeighbourList
    !> Add an external potential to a calculator
    procedure :: setExternalPotential => TDftbPlus_setExternalPotential
    !> Add an external electric field to a calculator
    procedure :: setExternalEfield => TDftbPlus_setExternalEfield
    !> Add external charges to a calculator
    procedure :: setExternalCharges => TDftbPlus_setExternalCharges
    !> Add reactive external charges to a calculator
    procedure :: setQDepExtPotGen => TDftbPlus_setQDepExtPotGen
    !> Obtain the DFTB+ energy
    procedure :: getEnergy => TDftbPlus_getEnergy
    !> Obtain the DFTB+ gradients
    procedure :: getGradients => TDftbPlus_getGradients
    !> Obtain the DFTB+ stress tensor
    procedure :: getStressTensor => TDftbPlus_getStressTensor
    !> Obtain the gradients of the external charges
    procedure :: getExtChargeGradients => TDftbPlus_getExtChargeGradients
    !> Get the gross (Mulliken) DFTB+ charges
    procedure :: getGrossCharges => TDftbPlus_getGrossCharges
    !> Get the CM5 DFTB+ charges
    procedure :: getCM5Charges => TDftbPlus_getCM5Charges
    !> Get the reference charges for neutral DFTB+ atoms
    procedure :: getRefCharges => TDftbPlus_getRefCharges
    !> Set the reference charges for neutral DFTB+ atoms
    procedure :: setRefCharges => TDftbPlus_setRefCharges
    !> Get electrostatic potential at specified points
    procedure :: getElStatPotential => TDftbPlus_getElStatPotential
    !> Return the number of DFTB+ atoms in the system
    procedure :: nrOfAtoms => TDftbPlus_nrOfAtoms
    !> Return the number of spin channels in the system
    procedure :: nrOfSpin => TDftbPlus_nrOfSpin
    !> Return the number of (k-point,spin chanel) pairs in the process group
    procedure :: nrOfLocalKS => TDftbPlus_nrOfLocalKS
    !> get (k-point,spin chanel) pairs in current process group
    procedure :: getLocalKS => TDftbPlus_getLocalKS
    !> Queries weights of k-points
    procedure :: getKWeights => TDftbPlus_getKWeights
    !> Returns size of the basis set
    procedure :: getBasisSize => TDftbPlus_getBasisSize
    !> Whether the system is described with real matrices (complex otherwise)
    procedure :: isHSReal => TDftbPlus_isHSReal
    !> Register callback function to be invoked on each evaluation of the density matrix
    procedure :: registerDMCallback => TDftbPlus_registerDMCallback
    !> Register callback function to be invoked on the first evaluation of the overlap matrix
    procedure :: registerSCallback => TDftbPlus_registerSCallback
    !> Register callback function to be invoked to import the overlap matrix
    procedure :: registerSetSCallback => TDftbPlus_registerSetSCallback
    !> Register callback function to be invoked on the first evaluation of the hamiltonian matrix
    procedure :: registerHCallback => TDftbPlus_registerHCallback
    !> Register callback function to be invoked to import the hamiltonian matrix
    procedure :: registerSetHCallback => TDftbPlus_registerSetHCallback
    !> Return the number of k-points in the DFTB+ calculation (1 if non-repeating)
    procedure :: nrOfKPoints => TDftbPlus_nrOfKPoints
    !> Return the number of spin channels in the DFTB+ calculation (1 if spin free, 2 for z spin
    !> polarised and 4 for non-collinear/spin-orbit)
    procedure :: nrOfSpinChannels => TDftbPlus_nrOfSpinChannels
    !> Check that the list of species names has not changed
    procedure :: checkSpeciesNames => TDftbPlus_checkSpeciesNames
    !> Replace species and redefine all quantities that depend on it
    procedure :: setSpeciesAndDependents => TDftbPlus_setSpeciesAndDependents
    !> Initialise electron and nuclear Ehrenfest dynamics
    procedure :: initializeTimeProp => TDftbPlus_initializeTimeProp
    !> Finalizes electron and nuclear Ehrenfest dynamics
    procedure :: finalizeTimeProp => TDftbPlus_finalizeTimeProp
    !> Do one propagator step for electrons and, if enabled, nuclei
    procedure :: doOneTdStep => TDftbPlus_doOneTdStep
    !> Set electric field for current propagation step of electrons and nuclei
    procedure :: setTdElectricField => TDftbPlus_setTdElectricField
    !> Set electric field for current propagation step of electrons and nuclei
    procedure :: setTdCoordsAndVelos => TDftbPlus_setTdCoordsAndVelos
    !> Set electric field for current propagation step of electrons and nuclei
    procedure :: getTdForces => TDftbPlus_getTdForces
    !> Check instance of DFTB+ is initialised
    procedure, private :: checkInit => TDftbPlus_checkInit
    !> Return the masses for each atom in the system
    procedure :: getAtomicMasses => TDftbPlus_getAtomicMasses
    !> Return the number of basis functions for each atom in the system
    procedure :: getNOrbitalsOnAtoms => TDftbPlus_getNOrbAtoms
    !> Get the maximum cutoff distance
    procedure :: getCutOff => TDftbPlus_getCutOff
    !> Return the number of orbitals
    procedure :: nrOfOrbitals => TDftbPlus_nrOfOrbitals
    !> get the DFTB+ eigenvalues
    procedure :: getEigenValues => TDftbPlus_getEigenValues
    !> get the DFTB+ eigenvectors
    procedure :: getEigenVectors => TDftbPlus_getEigenVectors
    !> get the DFTB+ hamiltonian and overlap matrices
    procedure :: getHamilOverl => TDftbPlus_getHamilOverl
    !> check and possibly invert the phase of frontier orbitals
    procedure :: checkInvertPhase => TDftbPlus_checkInvertPhase
    !> obtain the DFTB+ gradients due to FMO orbital/s
    procedure :: getFmoGradients => TDftbPlus_getFmoGradients
    !> get pointers to phase 1 of the DFTB-FMO calculation
    procedure :: getPointersToPhase1 => TDftbPlus_getPointersToPhase1
    !> init pointers to phase 1 of the DFTB-FMO calculation
    procedure :: initPointersToPhase1 => TDftbPlus_initPointersToPhase1
    !> set pointers to phase 1 of the DFTB-FMO calculation
    procedure :: setPointersToPhase1 => TDftbPlus_setPointersToPhase1
    !> run the phase 2 of the DFTB-FMO calculation, and get the Hamiltonian in FMO basis set
    procedure :: getFragmentBasedHamiltonian => TDftbPlus_getFragmentBasedHamiltonian
    !> Finalizer
    final :: TDftbPlus_final
  end type TDftbPlus


#:if not INSTANCE_SAFE_BUILD

  !> Nr. of existing instances (if build is not instance safe)
  integer :: nInstance_ = 0

#:endif


contains


  !> Return the version string for the current DFTB+ build
  subroutine getDftbPlusBuild(version)

    !> Version string for DFTB+
    character(:), allocatable, intent(out) :: version

    version = '${RELEASE}$'

  end subroutine getDftbPlusBuild


  !> Returns the DFTB+ API version
  subroutine getDftbPlusApi(major, minor, patch, instanceSafe)

    !> Major version number
    integer, intent(out) :: major

    !> Minor version number
    integer, intent(out) :: minor

    !> Patch level for API
    integer, intent(out) :: patch

    !> Whether API is instance safe
    logical, optional, intent(out) :: instanceSafe

    major = ${APIMAJOR}$
    minor = ${APIMINOR}$
    patch = ${APIPATCH}$
    if (present(instanceSafe)) then
      instanceSafe = instanceSafeBuild
    end if

  end subroutine getDftbPlusApi


  !> Finalizer for the DFTB+ input type
  subroutine TDftbPlusInput_final(this)

    !> Instance.
    type(TDftbPlusInput), intent(inout) :: this

    call TDftbPlusInput_destruct(this)

  end subroutine TDftbPlusInput_final


  !> Destructs the DFTB+ input type
  subroutine TDftbPlusInput_destruct(this)

    !> Instance.
    type(TDftbPlusInput), intent(inout) :: this

    if (associated(this%hsdTree)) then
      call destroyNode(this%hsdTree)
      this%hsdTree => null()
    end if

  end subroutine TDftbPlusInput_destruct


  !> Returns the root node of the input, so that it can be further processed
  subroutine TDftbPlusInput_getRootNode(this, root)

    !> Instance.
    class(TDftbPlusInput), intent(in) :: this

    !> Pointer to root node
    type(fnode), pointer, intent(out) :: root

    if (.not. associated(this%hsdTree)) then
      call error("Input has not been created yet!")
    end if
    call getChild(this%hsdTree, rootTag, root)

  end subroutine TDftbPlusInput_getRootNode


  !> Passes the information about the QM region to DFTB+
  subroutine TDftbPlusAtomList_get(instance, nAtom, speciesNames, species)

    !> Input containing the tree representation of the parsed HSD file.
    class(TDftbPlusAtomList), intent(out) :: instance

    !> Number of atoms
    integer, intent(in) :: nAtom

    !> Linked list of chemical symbols of elements (species names), size=nSpecies
    type(TListString), intent(inout) :: speciesNames

    !> Array of species for each atom, size=nAtom
    integer, intent(in) :: species(:)

    integer :: i
    character(3) :: s

    instance%nAtom = nAtom

    call init(instance%speciesNames)
    do i=1,len(speciesNames)
      call get(speciesNames, s, i)
      call append(instance%speciesNames, s)
    end do

    allocate(instance%species(nAtom))
    instance%species(1:nAtom) = species(1:nAtom)

  end subroutine TDftbPlusAtomList_get


  !> Insert the list of atoms into the input data structure
  subroutine TDftbPlusAtomList_addToInpData(instance, inpData)

    !> Input structure of the API
    class(TDftbPlusAtomList), intent(inout) :: instance

    !> Input data structure that will be in turn filled by parsing the HSD tree
    type(TInputData), intent(out), target :: inpData

    type(TGeometry), pointer :: geo

    ! adopted from subroutine readTGeometryGen_help
    geo => inpData%geom

    geo%nAtom = instance%nAtom
    geo%tPeriodic = .false.
    geo%tFracCoord = .false.
    geo%tHelical = .false.

    geo%nSpecies = len(instance%speciesNames)
    allocate(geo%speciesNames(geo%nSpecies))
    call asArray(instance%speciesNames, geo%speciesNames)

    ! Read in sequential and species indices.
    allocate(geo%species(geo%nAtom))
    allocate(geo%coords(3, geo%nAtom))

    geo%species(1:geo%nAtom) = instance%species(1:geo%nAtom)

    if (geo%nSpecies /= maxval(geo%species) .or. minval(geo%species) /= 1) then
      call error("Nr. of species and nr. of specified elements do not match.")
    end if

  end subroutine TDftbPlusAtomList_addToInpData


  !> Initialises a DFTB+ instance
  !!
  !! Note: due to some remaining global variables in the DFTB+ core, only one instance can be
  !! initialised within one process. Therefore, this routine can not be called twice, unless the
  !! TDftbPlus_destruct() has been called in between the inits (or the instance had already been
  !! finalized).  Otherwise the subroutine will stop.
  !!
  subroutine TDftbPlus_init(this, outputUnit, mpiComm, devNull)

    !> Instance
    type(TDftbPlus), intent(out) :: this

    !> Unit where to write the output (note: also error messages are written here)
    integer, intent(in), optional :: outputUnit

    !> MPI-communicator to use
    integer, intent(in), optional :: mpiComm

    !> Unit of the null device (you must open the null device and pass its unit number, if you open
    !> multiple TDftbPlus instances within an MPI-process)
    integer, intent(in), optional :: devNull

    integer :: stdOut

  #:if not INSTANCE_SAFE_BUILD
    if (nInstance_ /= 0) then
      call error("This build does not support multiple DFTB+ instances")
    end if
    nInstance_ = 1
  #:endif

    if (present(mpiComm) .and. .not. withMpi) then
      call error("MPI Communicator supplied to initialise a serial DFTB+ instance")
    end if

    if (present(outputUnit)) then
      stdOut = outputUnit
    else
      stdOut = output_unit
    end if

    call initGlobalEnv(outputUnit=outputUnit, mpiComm=mpiComm, devNull=devNull)
    allocate(this%env)
    allocate(this%main)
    call TEnvironment_init(this%env)
    this%env%tAPICalculation = .true.
    this%isInitialised = .true.

  end subroutine TDftbPlus_init


  !> Finalizer for TDftbPlus.
  subroutine TDftbPlus_final(this)

    !> Instance
    type(TDftbPlus), intent(inout) :: this

    call TDftbPlus_destruct(this)

  end subroutine TDftbPlus_final


  !> Destroys a DFTB+ calculation instance
  subroutine TDftbPlus_destruct(this)

    !> Instance
    type(TDftbPlus), intent(inout) :: this

    if (.not. this%isInitialised) return
    call this%checkInit()

    call this%main%destructProgramVariables()
    call this%env%destruct()
    deallocate(this%main, this%env)
    call destructGlobalEnv()
    this%isInitialised = .false.

    #:if not INSTANCE_SAFE_BUILD
      nInstance_ = 0
    #:endif

  end subroutine TDftbPlus_destruct


  !> Fills up the input by parsing an HSD file
  subroutine TDftbPlus_getInputFromFile(this, fileName, input)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Name of the file to parse
    character(len=*), intent(in) :: fileName

    !> Input containing the tree representation of the parsed HSD file.
    type(TDftbPlusInput), intent(out) :: input

    call this%checkInit()

    call readHsdFile(fileName, input%hsdTree)

  end subroutine TDftbPlus_getInputFromFile


  !> Creates an input with no entries.
  subroutine TDftbPlus_getEmptyInput(this, input)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Instance.
    type(TDftbPlusInput), intent(out) :: input

    type(fnode), pointer :: root, dummy

    call this%checkInit()

    input%hsdTree => createDocumentNode()
    root => createElement(rootTag)
    dummy => appendChild(input%hsdTree, root)

  end subroutine TDftbPlus_getEmptyInput


  !> Sets up the calculator using a given input.
  subroutine TDftbPlus_setupCalculator(this, input, atomList, outFileNameStub)

    !> Instance.
    class(TDftbPlus), target, intent(inout) :: this

    !> Representation of the DFTB+ input.
    type(TDftbPlusInput), intent(inout) :: input

    !> List of atoms and species for the QM region.
    type(TDftbPlusAtomList), intent(inout), optional :: atomList

    !> Output file name stub
    character(len=*), intent(in), optional :: outFileNameStub

    type(TParserFlags) :: parserFlags
    type(TInputData) :: inpData

    call this%checkInit()

    if (present(atomList)) then
      call atomList%add(inpData)
    end if

    call parseHsdTree(input%hsdTree, inpData, parserFlags)
    call doPostParseJobs(input%hsdTree, parserFlags)
    call this%main%initProgramVariables(inpData, this%env, outFileNameStub)

  end subroutine TDftbPlus_setupCalculator


  !> Sets the geometry in the calculator.
  subroutine TDftbPlus_setGeometry(this, coords, latVecs, origin)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Atomic coordinates in Bohr units. Shape: (3, nAtom).
    real(dp), intent(in) :: coords(:,:)

    !> Lattice vectors in Bohr units, stored column-wise. Shape: (3, 3).
    real(dp), intent(in), optional :: latVecs(:,:)

    !> Coordinate origin in Bohr units. Shape: (3).
    real(dp), intent(in), optional :: origin(:)

    call this%checkInit()

    call setGeometry(this%env, this%main, coords, latVecs, origin)

  end subroutine TDftbPlus_setGeometry


  !> Sets the neighbour list and skips the neighbour list creation in DFTB+
  subroutine TDftbPlus_setNeighbourList(this, nNeighbour, iNeighbour, neighDist, cutOff,&
      & coordNeighbours, neighbour2CentCell)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Number of neighbours of an atom in the central cell
    integer, intent(in) :: nNeighbour(:)

    !> References to the neighbour atoms for an atom in the central cell
    integer, intent(in) :: iNeighbour(:,:)

    !> Distances to the neighbour atoms for an atom in the central cell
    real(dp), intent(in) :: neighDist(:,:)

    !> Cutoff distance used for this neighbour list
    real(dp), intent(in) :: cutOff

    !> Coordinates of all neighbours
    real(dp), intent(in) :: coordNeighbours(:,:)

    !> Mapping between neighbour reference and atom index in the central cell
    integer, intent(in) :: neighbour2CentCell(:)

    call this%checkInit()

    call setNeighbourList(this%env, this%main, nNeighbour, iNeighbour, neighDist, cutOff,&
        & coordNeighbours, neighbour2CentCell)

  end subroutine TDftbPlus_setNeighbourList


  !> Sets an external potential.
  subroutine TDftbPlus_setExternalPotential(this, atomPot, potGrad)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Potential acting on each atom. Shape: (nAtom)
    real(dp), intent(in), optional :: atomPot(:)

    !> Gradient of the potential  on each atom. Shape: (3, nAtom)
    real(dp), intent(in), optional :: potGrad(:,:)

    call this%checkInit()

    if (allocated(this%main%solvation)) then
      if (this%main%solvation%isEFieldModified()) then
        call error("External fields currently unsupported for this solvent model")
      end if
    end if

    call setExternalPotential(this%main, atomPot=atomPot, potGrad=potGrad)

  end subroutine TDftbPlus_setExternalPotential


  !> Sets an external electric field.
  subroutine TDftbPlus_setExternalEfield(this, EFieldStr, EfieldVec)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Electric field amplitude
    real(dp), intent(in) :: EFieldStr

    !> Unitary electric field vector. Shape: (3)
    real(dp), intent(in) :: EfieldVec(:)

    call this%checkInit()

    call setExternalEfield(this%main, EFieldStr, EfieldVec)

  end subroutine TDftbPlus_setExternalEfield


  !> Sets external point charges.
  subroutine TDftbPlus_setExternalCharges(this, chargeCoords, chargeQs, blurWidths)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Coordinate of the external charges
    real(dp), intent(in) :: chargeCoords(:,:)

    !> Charges of the external point charges (sign convention: electron is negative)
    real(dp), intent(in) :: chargeQs(:)

    !> Widths of the Gaussian for each charge used for blurring (0.0 = no blurring)
    real(dp), intent(in), optional :: blurWidths(:)

    call this%checkInit()

    if (allocated(this%main%solvation)) then
      if (this%main%solvation%isEFieldModified()) then
        call error("External fields currently unsupported for this solvent model")
      end if
    end if

    call setExternalCharges(this%main, chargeCoords, chargeQs, blurWidths)

  end subroutine TDftbPlus_setExternalCharges


  !> Sets the generator for the population dependant external potential.
  subroutine TDftbPlus_setQDepExtPotGen(this, extPotGen)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Population dependant external potential generator
    class(TQDepExtPotGen), intent(in) :: extPotGen

    type(TQDepExtPotGenWrapper) :: extPotGenWrapper
    type(TQDepExtPotProxy) :: extPotProxy

    call this%checkInit()

    if (allocated(this%main%solvation)) then
      if (this%main%solvation%isEFieldModified()) then
        call error("External fields currently unsupported for this solvent model")
      end if
    end if

    allocate(extPotGenWrapper%instance, source=extPotGen)
    call TQDepExtPotProxy_init(extPotProxy, [extPotGenWrapper])
    call setQDepExtPotProxy(this%main, extPotProxy)

  end subroutine TDftbPlus_setQDepExtPotGen


  !> Return the energy of the current system.
  subroutine TDftbPlus_getEnergy(this, merminEnergy)

    !> Instance.
    class(TDftbPlus), intent(inout) :: this

    !> Mermin free energy.
    real(dp), intent(out) :: merminEnergy

    call this%checkInit()

    call getEnergy(this%env, this%main, merminEnergy)

  end subroutine TDftbPlus_getEnergy


  !> Returns the gradient on the atoms in the system.
  subroutine TDftbPlus_getGradients(this, gradients)

    !> Instance.
    class(TDftbPlus), intent(inout) :: this

    !> Gradients on the atoms.
    real(dp), intent(out) :: gradients(:,:)

    call this%checkInit()

    call getGradients(this%env, this%main, gradients)

  end subroutine TDftbPlus_getGradients


  !> Returns the stress tensor of the periodic system.
  subroutine TDftbPlus_getStressTensor(this, stresstensor)

    !> Instance.
    class(TDftbPlus), intent(inout) :: this

    !> Gradients on the atoms.
    real(dp), intent(out) :: stresstensor(:,:)

    call this%checkInit()

    call getStressTensor(this%env, this%main, stresstensor)

  end subroutine TDftbPlus_getStressTensor


  !> Returns the gradients on the external charges.
  !!
  !! This function may only be called if TDftbPlus_setExternalCharges was called before it
  !!
  subroutine TDftbPlus_getExtChargeGradients(this, gradients)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Gradients acting on the external charges.
    real(dp), intent(out) :: gradients(:,:)

    call this%checkInit()

    call getExtChargeGradients(this%main, gradients)

  end subroutine TDftbPlus_getExtChargeGradients


  !> Returns the gross (Mulliken) charges of each atom
  subroutine TDftbPlus_getGrossCharges(this, atomCharges)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Atomic gross charges.
    real(dp), intent(out) :: atomCharges(:)

    call this%checkInit()

    call getGrossCharges(this%env, this%main, atomCharges)

  end subroutine TDftbPlus_getGrossCharges


  !> Returns the CM5 charges of each atom
  subroutine TDftbPlus_getCM5Charges(this, atomCharges)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Atomic gross charges.
    real(dp), intent(out) :: atomCharges(:)

    call this%checkInit()

    call getCM5Charges(this%env, this%main, atomCharges)

  end subroutine TDftbPlus_getCM5Charges


  !> Get the reference atomic charges for the atoms of the system to be neutral
  subroutine TDftbPlus_getRefCharges(this, z0)

    !> Instance
    class(TDftbPlus), intent(in) :: this

    !> Atomic valence reference charges
    real(dp), intent(out) :: z0(:)

    real(dp), allocatable :: q0(:, :, :)
    integer :: mOrb, nAtom, nSpin

    call this%checkInit()

    if (this%main%uniqHubbU%mHubbU > 1) then
      call error("Reference charge call unsupported for shell resolved models")
    end if
    mOrb = this%main%orb%mOrb
    nAtom = nrOfAtoms(this%main)
    nSpin = this%main%nSpin
    allocate(q0(mOrb, nAtom, nspin))
    call getRefCharges(this%main, q0)
    z0(:) = sum(q0(:,:,1), dim=1)

  end subroutine TDftbPlus_getRefCharges


  !> Set the reference atomic charges for the atoms of the system to be neutral
  subroutine TDftbPlus_setRefCharges(this, z0)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Atomic valence reference charges
    real(dp), intent(in) :: z0(:)

    real(dp), allocatable :: q0(:, :, :)
    integer :: mOrb, nAtom, nSpin

    call this%checkInit()

    if (this%main%uniqHubbU%mHubbU > 1) then
      call error("Reference charge call unsupported for shell resolved models")
    end if
    mOrb = this%main%orb%mOrb
    nAtom = nrOfAtoms(this%main)
    nSpin = this%main%nSpin
    allocate(q0(mOrb, nAtom, nspin), source=0.0_dp)
    q0(1,:,1) = z0
    call setRefCharges(this%env, this%main, q0)

  end subroutine TDftbPlus_setRefCharges


  !> Returns electrostatic potential at specified points
  subroutine TDftbPlus_getElStatPotential(this, pot, locations)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Resulting potentials
    real(dp), intent(out) :: pot(:)

    !> Sites at which to calculate potential
    real(dp), intent(in) :: locations(:,:)

    call this%checkInit()

    call getElStatPotential(this%env, this%main, pot, locations)

  end subroutine TDftbPlus_getElStatPotential


  !> Returns the nr. of atoms in the system.
  function TDftbPlus_nrOfAtoms(this) result(nAtom)

    !> Instance
    class(TDftbPlus), intent(in) :: this

    !> Nr. of atoms
    integer :: nAtom

    call this%checkInit()

    nAtom = nrOfAtoms(this%main)

  end function TDftbPlus_nrOfAtoms


  !> Get (k-point,spin chanel) pairs in current process group
  subroutine TDftbPlus_getLocalKS(this, localKS)

    !> Instance
    class(TDftbPlus), intent(in) :: this

    !> The (K, S) tuples of the local processor group (localKS(1:2,iKS))
    !! Usage: iK = localKS(1, iKS); iS = localKS(2, iKS)
    integer, intent(out) :: localKS(:,:)

    call getLocalKS(this%main, localKS)

  end subroutine TDftbPlus_getLocalKS


  !> Queries weights of k-points
  subroutine TDftbPlus_getKWeights(this, KWeights)

    !> Instance
    class(TDftbPlus), intent(in) :: this

    !> Weights of k-points
    real(dp), intent(out) :: KWeights(:)

    KWeights(:) = this%main%kweight(:)


  end subroutine TDftbPlus_getKWeights


  !> Returns the nr. of spin channels in the system.
  function TDftbPlus_nrOfSpin(this) result(nSpin)

    !> Instance
    class(TDftbPlus), intent(in) :: this

    !> Nr. of spins channels
    integer :: nSpin

    call this%checkInit()

    nSpin = nrOfSpin(this%main)

  end function TDftbPlus_nrOfSpin


  !> Return the number of (k-point,spin chanel) pairs in the process group.
  function TDftbPlus_nrOfLocalKS(this) result(nLocalKS)

    !> Instance
    class(TDftbPlus), intent(in) :: this

    !> Nr. of (k-point,spin chanel) pairs
    integer :: nLocalKS

    call this%checkInit()

    nLocalKS = nrOfLocalKS(this%main)

  end function TDftbPlus_nrOfLocalKS


  !> Returns size of the basis set
  function TDftbPlus_getBasisSize(this) result(basisSize)

    !> Instance
    class(TDftbPlus), intent(in) :: this

    integer :: basisSize

    call this%checkInit()

    basisSize = this%main%denseDesc%fullSize

  end function TDftbPlus_getBasisSize


  !> Whether the system is described with real matrices
  function TDftbPlus_isHSReal(this) result(isHSReal)

    !> Instance
    class(TDftbPlus), intent(in) :: this

    logical isHSReal

    call this%checkInit()

    isHSReal = this%main%tRealHS

  end function TDftbPlus_isHSReal


  !> Register callback function to be invoked on each evaluation of the density matrix
  subroutine TDftbPlus_registerDMCallback(this, callback, aux_ptr)
    use dftbp_dftbplus_apicallback, only : TAPICallback, TDMHSCallbackFunc

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> callback function for DM export
    procedure(TDMHSCallbackFunc) :: callback

    !> pointer to a context object for the DM callback
    class(*), pointer :: aux_ptr

    call this%checkInit()
    call this%main%apicallback%registerDM(callback, aux_ptr)

  end subroutine TDftbPlus_registerDMCallback


  !> Register callback function to be invoked on the first evaluation of the overlap matrix
  subroutine TDftbPlus_registerSCallback(this, callback, aux_ptr)
    use dftbp_dftbplus_apicallback, only : TAPICallback, TDMHSCallbackFunc

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> callback function for S export
    procedure(TDMHSCallbackFunc) :: callback

    !> pointer to a context object for the S callback
    class(*), pointer :: aux_ptr

    call this%checkInit()
    call this%main%apicallback%registerS(callback, aux_ptr)

  end subroutine TDftbPlus_registerSCallback


  !> Register callback function to be invoked to import overlap matrix
  subroutine TDftbPlus_registerSetSCallback(this, callback, aux_ptr)
    use dftbp_dftbplus_apicallback, only : TAPICallback, TSetDMHSCallbackFunc

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> callback function for S import
    procedure(TSetDMHSCallbackFunc) :: callback

    !> pointer to a context object for the S callback
    class(*), pointer :: aux_ptr

    call this%checkInit()
    call this%main%apicallback%registerSetS(callback, aux_ptr)

  end subroutine TDftbPlus_registerSetSCallback


  !> Register callback function to be invoked on the first evaluation of the hamiltonian matrix
  subroutine TDftbPlus_registerHCallback(this, callback, aux_ptr)
    use dftbp_dftbplus_apicallback, only : TAPICallback, TDMHSCallbackFunc

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> callback function for H export
    procedure(TDMHSCallbackFunc) :: callback

    !> pointer to a context object for the H callback
    class(*), pointer :: aux_ptr

    call this%checkInit()
    call this%main%apicallback%registerH(callback, aux_ptr)

  end subroutine TDftbPlus_registerHCallback


  !> Register callback function to be invoked to import hamiltonian matrix
  subroutine TDftbPlus_registerSetHCallback(this, callback, aux_ptr)
    use dftbp_dftbplus_apicallback, only : TAPICallback, TSetDMHSCallbackFunc

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> callback function for H export
    procedure(TSetDMHSCallbackFunc) :: callback

    !> pointer to a context object for the H callback
    class(*), pointer :: aux_ptr

    call this%checkInit()
    call this%main%apicallback%registerSetH(callback, aux_ptr)

  end subroutine TDftbPlus_registerSetHCallback


  !> Returns the nr. of k-points describing the system.
  function TDftbPlus_nrOfKPoints(this) result(nKpts)

    !> Instance
    class(TDftbPlus), intent(in) :: this

    !> Nr. of k-points
    integer :: nKpts

    call this%checkInit()

    nKpts = nrOfKPoints(this%main)

  end function TDftbPlus_nrOfKPoints


  !> Returns the nr. of spin channels
  function TDftbPlus_nrOfSpinChannels(this) result(nSpin)

    !> Instance
    class(TDftbPlus), intent(in) :: this

    !> Nr. of spin channels
    integer :: nSpin

    call this%checkInit()

    nSpin = this%main%nSpin

  end function TDftbPlus_nrOfSpinChannels


  !> Returns the atomic masses for each atom in the system.
  subroutine TDftbPlus_getAtomicMasses(this, mass)

    !> Instance
    class(TDftbPlus), intent(in) :: this

    !> Masses for each species of the system
    real(dp), intent(out) :: mass(:)

    call this%checkInit()

    call getAtomicMasses(this%main, mass)

  end subroutine TDftbPlus_getAtomicMasses


  !> Returns the number of orbitals for each atom in the system
  subroutine TDftbPlus_getNOrbAtoms(this, nOrbs)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Number of basis functions associated with each atom
    integer, intent(out) :: nOrbs(:)

    nOrbs(:) = this%main%orb%nOrbAtom

  end subroutine TDftbPlus_getNOrbAtoms


  !> Gets the cutoff distance used for interactions
  function TDftbPlus_getCutOff(this) result(cutOff)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Cutoff distance
    real(dp) :: cutOff

    call this%checkInit()

    cutOff = getCutOff(this%main)

  end function TDftbPlus_getCutOff


  !> Checks whether the type is already initialized and stops the code if not.
  subroutine TDftbPlus_checkInit(this)

    !> Instance.
    class(TDftbPlus), intent(in) :: this

    if (.not. this%isInitialised) then
      call error("Received uninitialized TDftbPlus instance")
    end if

  end subroutine TDftbPlus_checkInit


  !> Reads out the atomic angular momenta from an SK-file
  !!
  !! NOTE: This only works with handcrafted (non-standard) SK-files, where the nr. of shells
  !!   has been added as 3rd entry to the first line of the homo-nuclear SK-files.
  !!
  function getMaxAngFromSlakoFile(slakoFile) result(maxAng)

    !> Instance.
    character(len=*), intent(in) :: slakoFile

    !> Maximal angular momentum found in the file
    integer :: maxAng

    real(dp) :: dr
    integer :: nGridPoints, nShells
    type(TFileDescr) :: fd

    call openFile(fd, slakoFile, mode="r")
    read(fd%unit, *) dr, nGridPoints, nShells
    call closeFile(fd)
    maxAng = nShells - 1

  end function getMaxAngFromSlakoFile


  !> Converts atom types to species
  subroutine convertAtomTypesToSpecies(typeNumbers, species, speciesNames, typeNames)

    !> Type number of each atom is the system. It can be arbitrary number (e.g. atomic number)
    integer, intent(in) :: typeNumbers(:)

    !> Species index for each atom (1 for the first atom type found, 2 for the second, etc.)
    integer, allocatable, intent(out) :: species(:)

    !> Names of each species, usually X1, X2 unless typeNames have been specified.
    character(len=*), allocatable, intent(out) :: speciesNames(:)

    !> Array of type names, indexed by the type numbers.
    character(len=*), intent(in), optional :: typeNames(:)

    integer, allocatable :: uniqueTypeNumbers(:)
    integer :: nAtom, nSpecies
    integer :: iAt, iSp, curType

    nAtom = size(typeNumbers)

    allocate(uniqueTypeNumbers(nAtom))
    nSpecies = 0
    do iAt = 1, nAtom
      curType = typeNumbers(iAt)
      if (.not. any(uniqueTypeNumbers(1:nSpecies) == curType)) then
        nSpecies = nSpecies + 1
        uniqueTypeNumbers(nSpecies) = curType
      end if
    end do

    allocate(species(nAtom))
    do iSp = 1, nSpecies
      where (typeNumbers == uniqueTypeNumbers(iSp))
        species = iSp
      end where
    end do

    allocate(speciesNames(nSpecies))
    do iSp = 1, nSpecies
      if (present(typeNames)) then
        speciesNames(iSp) = typeNames(uniqueTypeNumbers(iSp))
      else
        write(speciesNames(iSp), "(A,I0)") "X", iSp
      end if
    end do

  end subroutine convertAtomTypesToSpecies


  !> Check whether speciesNames has changed between calls to DFTB+
  subroutine TDftbPlus_checkSpeciesNames(this, inputSpeciesNames)

    !> Instance
    class(TDftbPlus),  intent(inout) :: this

    !> Chemical species labels
    character(len=*), intent(in) :: inputSpeciesNames(:)

    logical :: tSpeciesNameChanged

    call this%checkInit()

    tSpeciesNameChanged = checkSpeciesNames(this%env, this%main, inputSpeciesNames)

    if(tSpeciesNameChanged)then
      call error('speciesNames has changed between calls to DFTB+. This will cause erroneous&
          & results.' // newline // 'Instead call destruct and then fully re-initialize.')
    else
       continue
    endif

  end subroutine TDftbPlus_checkSpeciesNames


  !> Set species and all variables/data dependent on it
  subroutine TDftbPlus_setSpeciesAndDependents(this, inputSpeciesNames, inputSpecies)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Type of each atom (nAllAtom)
    integer, intent(in) :: inputSpecies(:)

    !> Labels of atomic species (nSpecies)
    character(len=*), intent(in) :: inputSpeciesNames(:)

    call this%checkInit()
    call this%checkSpeciesNames(inputSpeciesNames)
    call updateDataDependentOnSpeciesOrdering(this%env, this%main, inputSpecies)

  end subroutine TDftbPlus_setSpeciesAndDependents


  !> Initialise propagators for electron and nuclei dynamics
  subroutine TDftbPlus_initializeTimeProp(this, dt, tdFieldThroughAPI, tdCoordsAndVelosThroughAPI)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Time step
    real(dp), intent(in) :: dt

    !> Field will be provided through the API?
    logical, intent(in) :: tdFieldThroughAPI

    !> Coords and velocities will be provided at each step through the API?
    logical, intent(in) :: tdCoordsAndVelosThroughAPI

    call initializeTimeProp(this%env, this%main, dt, tdFieldThroughAPI, tdCoordsAndVelosThroughAPI)

  end subroutine TDftbPlus_initializeTimeProp


  !> Initialise propagators for electron and nuclei dynamics
  subroutine TDftbPlus_finalizeTimeProp(this)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    call finalizeTimeProp(this%main)

  end subroutine TDftbPlus_finalizeTimeProp


  !> Propagate one time step for electron and nuclei dynamics
  subroutine TDftbPlus_doOneTdStep(this, iStep, dipole, energy, atomNetCharges,&
      & coord, force, occ)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Present step of dynamics
    integer, intent(in) :: iStep

    !> Dipole moment
    real(dp), optional, intent(out) :: dipole(:,:)

    !> Data type for energy components and total
    real(dp), optional, intent(out) :: energy

    !> Negative gross charge
    real(dp), optional, intent(out) :: atomNetCharges(:,:)

    !> Atomic coordinates
    real(dp), optional, intent(out) :: coord(:,:)

    !> Forces (3, nAtom)
    real(dp), optional, intent(out) :: force(:,:)

    !> Molecular orbital projected populations
    real(dp), optional, intent(out) :: occ(:)

    call doOneTdStep(this%env, this%main, iStep, dipole=dipole, energy=energy,&
        & atomNetCharges=atomNetCharges, coordOut = coord, force=force, occ=occ)

  end subroutine TDftbPlus_doOneTdStep


  !> Sets electric field for td propagation
  subroutine TDftbPlus_setTdElectricField(this, field)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Electric field components
    real(dp), intent(in) :: field(3)

    if (allocated(this%main%solvation)) then
      if (this%main%solvation%isEFieldModified()) then
        call error("External fields currently unsupported for this solvent model")
      end if
    end if

    call setTdElectricField(this%main, field)

  end subroutine TDftbPlus_setTdElectricField


  !> Set atomic coordinates and velocities for MD
  subroutine TDftbPlus_setTdCoordsAndVelos(this, coords, velos)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Coordinates
    real(dp), intent(in) :: coords(3, this%main%nAtom)

    !> Velocities
    real(dp), intent(in) :: velos(3, this%main%nAtom)

    call setTdCoordsAndVelos(this%main, coords, velos)

  end subroutine TDftbPlus_setTdCoordsAndVelos


  !> Returns forces from time dependent propagation
  subroutine TDftbPlus_getTdForces(this, forces)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Forces (3, nAtom)
    real(dp), intent(out) :: forces(:,:)

    call getTdForces(this%main, forces)

  end subroutine TDftbPlus_getTdForces


  !> Returns the nr. of atoms in the system.
  function TDftbPlus_nrOfOrbitals(this) result(nOrb)

    !> Instance
    class(TDftbPlus), intent(in) :: this

    !> Nr. of atoms
    integer :: nOrb

    call this%checkInit()

    nOrb = nrOfOrbitals(this%main)

  end function TDftbPlus_nrOfOrbitals


  !> Returns the eigenvalues (of each orbital)
  subroutine TDftbPlus_getEigenValues(this, eigVal)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Eigenvalues.
    real(dp), intent(out) :: eigVal(:)

    call this%checkInit()

    call getEigenValues(this%env, this%main, eigVal)

  end subroutine TDftbPlus_getEigenValues


  !> Returns the eigenvectors (orbitals)
  subroutine TDftbPlus_getEigenVectors(this, eigVec)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Eigenvectors.
    real(dp), intent(out) :: eigVec(:,:)

    call this%checkInit()

    call getEigenVectors(this%env, this%main, eigVec)

  end subroutine TDftbPlus_getEigenVectors


  !> Returns the Hamiltonian and overlap matrices
  subroutine TDftbPlus_getHamilOverl(this, hamil, overl)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Hamiltonian matrix.
    real(dp), intent(out) :: hamil(:,:)

    !> Overlap matrix.
    real(dp), intent(out) :: overl(:,:)

    !> Error status
    type(TStatus) :: errStatus

    call this%checkInit()

    call getHamilOverl(this%env, this%main, hamil, overl, errStatus)

  end subroutine TDftbPlus_getHamilOverl


  !> Check and possibly invert the phase of frontier orbitals
  subroutine TDftbPlus_checkInvertPhase(this, atomIndexSign, nFrontiers, frontiers, firstStep)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> set of three atoms to define the plane of the molecule
    integer, intent(in) :: atomIndexSign(3)

    !> number of frontier orbitals to check the phase
    integer, intent(in) :: nFrontiers

    !> indices of frontier orbitals to check the phase
    integer, intent(in) :: frontiers(:)

    !> is this the first step of the simulation?
    logical, intent(in) :: firstStep

  ! call checkInvertPhase(this%env, this%main, atomIndexSign, nFrontiers, frontiers(1:nFrontiers), firstStep)
    call checkInvertPhase(this%env, atomIndexSign, nFrontiers, frontiers(1:nFrontiers), firstStep,&
        & this%main%denseDesc%nOrb, this%main%denseDesc%iAtomStart, this%main%coord0,&
        & this%main%SSqrReal, this%main%eigVecsReal(:,:,1), this%main%oldEigVecsReal,&
        & this%main%frontierOverlap)

  end subroutine TDftbPlus_checkInvertPhase


  !> Returns the gradient due to FMO orbital/s
  subroutine TDftbPlus_getFmoGradients(this, filling, gradients)

    !> Instance.
    class(TDftbPlus), intent(inout) :: this

    !> Filling of the orbitals (typically, one of the orbitals has 1, and the others have 0)
    real(dp), intent(in) :: filling(:)

    !> Gradients on the atoms.
    real(dp), intent(out) :: gradients(:,:)

    !> Error status
    type(TStatus) :: errStatus

    call this%checkInit()

    call fmoGradients(this%env, this%main%scc, this%main%isExtField, this%main%nonSccDeriv,&
        & this%main%eigVecsReal(:,:,1), this%main%eigen(:,1,1), filling, this%main%qOutput,&
        & this%main%q0, this%main%skHamCont, this%main%skOverCont, this%main%neighbourList, this%main%symNeighbourList,&
        & this%main%nNeighbourSK, this%main%nNeighbourCamSym, this%main%species, this%main%img2CentCell, this%main%orb,&
        & this%main%potential, this%main%coord, this%main%thirdOrd, this%main%qDepExtPot,&
        & this%main%hybridXc, this%main%SSqrReal, this%main%ints, this%main%denseDesc,&
        & this%main%iSparseStart, this%main%tRealHS, this%main%densityMatrix, gradients, errStatus)

  end subroutine TDftbPlus_getFmoGradients


  !> This is called before after finishing the phase 2 of the FMO calculation:
  !> Get the pointers to data from the phase 1 of the FMO DFTB+ calculation
  subroutine TDftbPlus_getPointersToPhase1(this, ptrsPhase1)

    !> Instance.
    class(TDftbPlus), intent(inout), target :: this

    !> Structure to store the data into
    type(TPointersToPhase1), intent(out) :: ptrsPhase1

    ! number of atoms
    ptrsPhase1%nAtom = this%main%nAtom

    ! number of orbitals
    ptrsPhase1%nOrb = this%main%nOrb

    ! sparse matrix, overlap
  ! ptrsPhase1%denseOver = this%main%SSqrReal
  ! ptrsPhase1%denseOver => this%main%SSqrReal
    allocate(ptrsPhase1%denseOver(this%main%nOrb,this%main%nOrb))
    call unpackHS(ptrsPhase1%denseOver, this%main%ints%overlap, this%main%neighbourList%iNeighbour,&
        & this%main%nNeighbourSK, this%main%denseDesc%iAtomStart, this%main%iSparseStart,&
        & this%main%img2CentCell)

    ! sparse matrix, charge-independent Hamiltonian
    allocate(ptrsPhase1%denseH0(this%main%nOrb,this%main%nOrb))
    call unpackHS(ptrsPhase1%denseH0, this%main%H0, this%main%neighbourList%iNeighbour,&
        & this%main%nNeighbourSK, this%main%denseDesc%iAtomStart, this%main%iSparseStart,&
        & this%main%img2CentCell)

    ! dense matrix, eigenvectors
    ptrsPhase1%eigVec = this%main%eigVecsReal(:,:,1)

    ! 1D array, eigenvalues (orbital energies)
    ptrsPhase1%eigVal = this%main%eigen(:,1,1)

    ! charge per atom (mOrb, atom, spin channel); spin channel == 1
  ! ptrsPhase1%qInput(:,:,:) = this%main%qOutput
    ptrsPhase1%qInput => this%main%qOutput

    ! charge per atomic shell (shell, atom, spin channel); spin channel == 1
  ! ptrsPhase1%chargePerShell(:,:,:) = this%main%chargePerShell
    ptrsPhase1%chargePerShell => this%main%chargePerShell

  end subroutine TDftbPlus_getPointersToPhase1


  !> Sets the number of sites for DFTB phase 2 calculation
  subroutine TDftbPlus_initPointersToPhase1(this, nSite)

    !> Instance.
    class(TDftbPlus), intent(inout) :: this

    !> Number of sites / fragments for FMO
    integer :: nSite

    call this%checkInit()

    this%nSite = nSite
    if (allocated(this%ptrsPhase1)) then
      deallocate(this%ptrsPhase1)
    end if
    allocate(this%ptrsPhase1(nSite))

  end subroutine TDftbPlus_initPointersToPhase1


  !> This is called before starting the phase 2 of the FMO calculation:
  !> Set the pointers to data from the phase 1 of the FMO DFTB+ calculation
 !subroutine TDftbPlus_setPointersToPhase1(this, iSite, nAtom, nOrb, nFO, iHOMO, denseOver,&
 !    & denseH0, eigVec, eigVal, qInput, chargePerShell)
  subroutine TDftbPlus_setPointersToPhase1(this, iSite, ptrsPhase1, nFO, iHOMO)

    !> Instance.
    class(TDftbPlus), intent(inout) :: this

    !> Which site? (index into the array of TPointersToPhase1)
    integer, intent(in) :: iSite

    !> Structure to copy the data from
    type(TPointersToPhase1), intent(in) :: ptrsPhase1

  ! ! number of atoms
  ! integer, intent(in) :: nAtom

  ! ! number of orbitals
  ! integer, intent(in) :: nOrb

    ! number of frontier / fragment orbitals to consider
    integer, intent(in) :: nFO

    ! which orbital is HOMO/LUMO?
    integer, intent(in) :: iHOMO

  ! ! sparse matrix, overlap
  ! real(dp), allocatable, intent(in) :: denseOver(:,:)

  ! ! sparse matrix, charge-independent Hamiltonian
  ! real(dp), allocatable, intent(in) :: denseH0(:,:)

  ! ! dense matrix, eigenvectors
  ! real(dp), allocatable, intent(in) :: eigVec(:,:)

  ! ! 1D array, eigenvalues (orbital energies)
  ! real(dp), allocatable, intent(in) :: eigVal(:)

  ! ! charge per atom (mOrb, atom, spin channel); spin channel == 1
  ! real(dp), allocatable, intent(in) :: qInput(:,:,:)

  ! ! charge per atomic shell (shell, atom, spin channel); spin channel == 1
  ! real(dp), allocatable, intent(in) :: chargePerShell(:,:,:)

    call this%checkInit()

    this%ptrsPhase1(iSite) = ptrsPhase1

  ! ! number of atoms
  ! this%ptrsPhase1(iSite)%nAtom = nAtom

  ! ! number of orbitals
  ! this%ptrsPhase1(iSite)%nOrb = nOrb

  ! ! number of frontier / fragment orbitals to consider
  ! this%ptrsPhase1(iSite)%nFO = nFO

  ! ! which orbital is HOMO/LUMO?
  ! this%ptrsPhase1(iSite)%iHOMO = iHOMO

  ! ! sparse matrix, overlap
  ! this%ptrsPhase1(iSite)%denseOver = denseOver

  ! ! sparse matrix, charge-independent Hamiltonian
  ! this%ptrsPhase1(iSite)%denseH0 = denseH0

  ! ! dense matrix, eigenvectors
  ! this%ptrsPhase1(iSite)%eigVec = eigVec

  ! ! 1D array, eigenvalues (orbital energies)
  ! this%ptrsPhase1(iSite)%eigVal = eigVal

  ! ! charge per atom (mOrb, atom, spin channel); spin channel == 1
  ! this%ptrsPhase1(iSite)%qInput = qInput

  ! ! charge per atomic shell (shell, atom, spin channel); spin channel == 1
  ! this%ptrsPhase1(iSite)%chargePerShell = chargePerShell

    this%ptrsPhase1(iSite)%nFO = nFO
    this%ptrsPhase1(iSite)%iHOMO = iHOMO

  end subroutine TDftbPlus_setPointersToPhase1


  !> Returns the Hamiltonian matrix in the FMO basis
  !> This invokes a DFTB phase 2 calculation
  subroutine TDftbPlus_getFragmentBasedHamiltonian(this, TijOrtho)

    !> Instance.
    class(TDftbPlus), intent(inout) :: this

    !> Mermin free energy.
    real(dp), allocatable, intent(out) :: TijOrtho(:,:)

    !> Error status
    type(TStatus) :: errStatus

    call this%checkInit()

    call getFragmentBasedHamiltonian(this%env, this%main, this%nSite, this%ptrsPhase1, TijOrtho,&
        & errStatus)

  end subroutine TDftbPlus_getFragmentBasedHamiltonian

end module dftbp_mmapi
