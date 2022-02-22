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

  implicit none

  private

  public :: TCpeInp

  !> Data type for initial values for CPE
  type :: TCpeInp

    !> Electronegativity per species
    real(dp), allocatable :: electronegativity(:)

    !> Chemical hardness per species
    real(dp), allocatable :: hardness(:)

    !> Covalent radius per species
    real(dp), allocatable :: radius(:)

    !> Total charge of the molecule
    real(dp) :: totalCharge

  end type TCpeInp

end module dftbp_cpeinp

