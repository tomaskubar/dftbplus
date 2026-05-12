!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Modified copy of densitymatrix.F90
!> Change from the previous version:
!>   the "full" procedures are copied, and not the previously used "sparse" ones
module dftbp_dftb_fmo_densitymatrix
  use dftbp_common_accuracy, only : dp
  use dftbp_math_blasroutines, only : herk

  implicit none

  public :: fmo_density_matrix_real
  public :: fmo_energy_density_matrix_real


contains


  !> Make a regular density matrix for the real wave-function case
  !! adopted from: subroutine fullDensityMatrix_real(dm, eigenvecs, filling)
  subroutine fmo_density_matrix_real(dm, eigenvecs, filling)

    !> The resulting nOrb*nOrb density matrix
    !! with only the elements of interest calculated, instead of the whole dm
    real(dp), intent(out) :: dm(:,:)

    !> The eigenvectors of the system
    real(dp), intent(in) :: eigenvecs(:,:)

    !> The occupation numbers of the orbitals
    real(dp), intent(in) :: filling(:)

    integer :: ii, nLevels
    real(dp), allocatable :: tmpEigen(:,:)

    @:ASSERT(all(shape(eigenvecs) == shape(dm)))
    @:ASSERT(size(eigenvecs,dim=1) == size(eigenvecs,dim=2))
    @:ASSERT(size(eigenvecs,dim=1) == size(filling))

    dm(:,:) = 0.0_dp
    ! do not restrict the procedure to occupied levels
    nLevels = size(filling)

    allocate(tmpEigen(nLevels, nLevels))

    do ii = 1, nLevels
      tmpEigen(:,ii) = sqrt(filling(ii)) * eigenvecs(:,ii)
    end do

    call herk(dm, tmpEigen)

    deallocate(tmpEigen)

  end subroutine fmo_density_matrix_real


  !> Make an energy weighted density matrix for the real wave-function case
  !! adopted from: subroutine fullEnergyDensityMatrix_real(dm, eigenvecs, filling, eigen)
  subroutine fmo_energy_density_matrix_real(dm, eigenvecs, filling, eigen)

    !> The resulting nOrb*nOrb density matrix
    real(dp), intent(out) :: dm(:,:)

    !> The eigenvectors of the system
    real(dp), intent(in) :: eigenvecs(:,:)

    !> The occupation numbers of the orbitals
    real(dp), intent(in) :: filling(:)

    !> Eigenvalues of the system
    real(dp), intent(in) :: eigen(:)

    integer :: ii, nLevels
    real(dp) :: fillProduct(size(filling))
    real(dp), allocatable :: tmpEigen(:,:)

    @:ASSERT(all(shape(eigenvecs) == shape(dm)))
    @:ASSERT(size(eigenvecs,dim=1) == size(eigenvecs,dim=2))
    @:ASSERT(size(eigenvecs,dim=1) == size(filling))
    @:ASSERT(size(eigen) == size(filling))

    dm(:,:) = 0.0_dp
    ! do not restrict the procedure to occupied levels
    nLevels = size(filling)

    fillProduct(1:nLevels) = filling(1:nLevels) * eigen(1:nLevels)

    allocate(tmpEigen(nLevels, nLevels))
    
    do ii = 1, nLevels
      tmpEigen(:,ii) = sqrt(abs(fillProduct(ii))) * eigenvecs(:,ii)
    end do

    call herk(dm, tmpEigen, alpha=sign(1.0_dp, maxval(fillProduct(1:nLevels))))

    deallocate(tmpEigen)

  end subroutine fmo_energy_density_matrix_real


end module dftbp_dftb_fmo_densitymatrix
