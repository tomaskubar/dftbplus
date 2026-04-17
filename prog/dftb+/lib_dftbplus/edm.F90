#:include 'common.fypp'

module dftbp_fmo_edm

  use dftbp_accuracy
  use dftbp_assert
  use dftbp_orbitals, only: TOrbitals
  use dftbp_sorting, only: heap_sort, unique

  public :: fmo_sp_energy_density_matrix_real

  contains

  ! COPIED AND MODIFIED FROM:
  !   densitymatrix.F90, makeDensityMatrix -> sp_energy_density_matrix_real
  !> Make an energy weighted density matrix for the real wave-function case
  subroutine fmo_sp_energy_density_matrix_real(dm, eigenvecs, filling, eigen, iNeighbour,&
      & nNeighbourSK, orb, iAtomStart, img2CentCell)

    !> the resulting nOrb*nOrb density matrix with only the elements of interest
    !> calculated, instead of the whole dm
    real(dp), intent(out) :: dm(:,:)

    !> the eigenvectors of the system
    real(dp), intent(in) :: eigenvecs(:,:)

    !> the occupation numbers of the orbitals
    real(dp), intent(in) :: filling(:)

    !> eigenvalues of the system
    real(dp), intent(in) :: eigen(:)

    !> Neighbour list for each atom (First index from 0!)
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for each atom (incl. itself).
    integer, intent(in) :: nNeighbourSK(:)

    !> Information about the orbitals.
    type(TOrbitals), intent(in) :: orb

    !> Atom offset for the squared Hamiltonian
    integer, intent(in) :: iAtomStart(:)

    !> Atomic mapping indexes.
    integer, intent(in) :: img2CentCell(:)

    integer :: iAt1, iNeigh1, jj, nOrb1, nOrb2, nAtom
    integer :: nLevels, start1, start2
    integer, allocatable :: inCellNeighbour(:,:)
    integer, allocatable :: nInCellNeighbour(:)
    real(dp), allocatable :: tmpEigen(:,:)

    @:ASSERT(all(shape(eigenvecs) == shape(dm)))
    @:ASSERT(size(eigenvecs,dim=1) == size(eigenvecs,dim=2))
    @:ASSERT(size(eigenvecs,dim=1) == size(filling))
    @:ASSERT(size(eigen) == size(filling))

    allocate(inCellNeighbour(0:size(iNeighbour,dim=1),size(iNeighbour,dim=2)))
    allocate(nInCellNeighbour(size(iNeighbour,dim=2)))

    nAtom = size(orb%nOrbAtom)
    dm(:,:) = 0.0_dp

    inCellNeighbour(:,:) = 0
    nInCellNeighbour(:) = 0

    ! do not restrict the procedure to occupied levels
    nLevels = size(filling)
!   if (abs(filling(ii)) >= epsilon(1.0_dp)) then

    allocate(tmpEigen(nLevels,orb%mOrb))
!   !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iAt1) SCHEDULE(RUNTIME)
    do iAt1 = 1, nAtom
      inCellNeighbour(0:nNeighbourSK(iAt1),iAt1) =&
          & img2CentCell(iNeighbour(0:nNeighbourSK(iAt1),iAt1))
      call heap_sort(inCellNeighbour(:nNeighbourSK(iAt1),iAt1))
      nInCellNeighbour(iAt1) = unique(inCellNeighbour(:,iAt1), nNeighbourSK(iAt1)+1) - 1
    end do
!   !$OMP  END PARALLEL DO

    do iAt1 = 1, nAtom
      nOrb1 = orb%nOrbAtom(iAt1)
      start1 = iAtomStart(iAt1)
      tmpEigen(1:nLevels,1:nOrb1) = transpose(eigenvecs(start1:start1+nOrb1-1,1:nLevels))

!     !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(jj) SCHEDULE(RUNTIME)
      do jj = 1, nLevels
        tmpEigen(jj,1:nOrb1) = filling(jj)*eigen(jj)*tmpEigen(jj,1:nOrb1)
      end do
!     !$OMP  END PARALLEL DO
!     !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iNeigh1,start2,nOrb2) SCHEDULE(RUNTIME)
      do iNeigh1 = 0, nInCellNeighbour(iAt1)
        start2 = iAtomStart(inCellNeighbour(iNeigh1,iAt1))
        nOrb2 = orb%nOrbAtom(inCellNeighbour(iNeigh1,iAt1))
        dm(start2:start2+nOrb2-1, start1:start1+nOrb1-1) =&
            & matmul(eigenvecs(start2:start2+nOrb2-1,1:nLevels), tmpEigen(1:nLevels,1:nOrb1))
      end do
!     !$OMP  END PARALLEL DO
    end do

  end subroutine fmo_sp_energy_density_matrix_real

end module dftbp_fmo_edm