!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Chemical potential equilibration (CPE)
!>
module dftbp_cpeinp

  use dftbp_accuracy
  use dftbp_machinelearning
! use dftbp_machinelearning_sf
! use dftbp_machinelearning_nn

  implicit none

  private

  public :: TCpeInp

  !> Data type for initial values for CPE
  type :: TCpeInp

    !> Electronegativity is given simply by numerical values
    logical :: tElectronegValues

    !> Electronegativity values per species
    real(dp), allocatable :: electronegativity(:)

    !> Electronegativity is represented by neural nets
    logical :: tElectronegNeuralNet

  ! !> Symmetry function related data
  ! type(TMLSymmetryFunctionsInp), allocatable :: electronegativitySF

  ! !> Neural network related data
  ! type(TMLNeuralNetInp), allocatable :: electronegativityNN

    !> Neural net and symmetry function data for calculation of electronegativities
    type(TMachineLearningInp), allocatable :: electronegativityML

    !> Chemical hardness per species
    real(dp), allocatable :: hardness(:)

    !> Covalent radius per species
    real(dp), allocatable :: radius(:)

    !> Total charge of the molecule
    real(dp) :: totalCharge

  end type TCpeInp

end module dftbp_cpeinp

