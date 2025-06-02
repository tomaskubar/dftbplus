!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> The main routines for DFTB+
module dftbp_fmogradient
  use dftbp_assert
  use dftbp_constants
  use dftbp_globalenv
  use dftbp_environment
  use dftbp_densedescr
  use dftbp_nonscc
  use dftbp_eigenvects
  use dftbp_populations
  use dftbp_forces
  use dftbp_scc
  use dftbp_hamiltonian
! use dftbp_externalcharges
  use dftbp_potentials
  use dftbp_sparse2dense
 !use dftbp_blasroutines
  use dftbp_shift
  use dftbp_mainio
  use dftbp_thirdorder, only : TThirdOrder
  use dftbp_rangeseparated, only : TRangeSepFunc
  use dftbp_slakocont
  use dftbp_lapackroutines
  use dftbp_elstatpot, only : TElStatPotentials
  use dftbp_qdepextpotproxy, only : TQDepExtPotProxy
  use dftbp_initprogram
! use dftbp_main, only : handleCoordinateChange
  use dftbp_sorting
  use dftbp_periodic, only : getSparseDescriptor
  use dftbp_fmo

  implicit none
  private

  public :: fmoGradients, fmoGradientsOffdiag

contains

  !> Calculates the gradients
  subroutine fmoGradients(env, sccCalc, nonSccDeriv, eigVecs, eigen, filling, qOutput, q0,&
      & skHamCont, skOverCont, neighbourList, nNeighbourSK, species, img2CentCell, orb, potential,&
      & coord, thirdOrd, qDepExtPot, rangeSep, SSqrReal, over, denseDesc, iSparseStart,&
      & deltaRhoOutSqr, derivs)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> SCC module internal variables
    type(TScc), allocatable, intent(in) :: sccCalc

    !> method for calculating derivatives of S and H0
    type(TNonSccDiff), intent(in) :: nonSccDeriv

    !> the eigenvectors of the system
    real(dp), intent(in) :: eigVecs(:,:)

    !> the eigenvalues
    real(dp), intent(in) :: eigen(:)

    !> filling of the fragment orbitals -- active orbitals only (1 for a single FMO, or so...)
    real(dp), intent(in) :: filling(:)

    !> electron populations (may be unallocated for non-scc case)
    real(dp), allocatable, intent(in) :: qOutput(:,:,:)

    !> reference atomic charges (may be unallocated for non-scc case)
    real(dp), allocatable, intent(in) :: q0(:,:,:)

    !> non-SCC hamiltonian information
    type(TSlakoCont), intent(in) :: skHamCont

    !> overlap information
    type(TSlakoCont), intent(in) :: skOverCont

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> species of all atoms in the system
    integer, intent(in) :: species(:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !>  potential acting on the system
    type(TPotentials), intent(in) :: potential

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> Is 3rd order SCC being used
    type(TThirdOrder), intent(inout), allocatable :: thirdOrd

    !> Population dependant external potential
    type(TQDepExtPotProxy), intent(inout), allocatable :: qDepExtPot

    !> Data from rangeseparated calculations
    type(TRangeSepFunc), intent(inout), allocatable :: rangeSep

    !> dense overlap matrix, required for rangeSep
    real(dp), intent(inout), allocatable :: SSqrReal(:,:)

    !> sparse overlap matrix, required for rangeSep
    real(dp), intent(in) :: over(:)

    !> Dense matrix descriptor,required for rangeSep
    type(TDenseDescr), intent(in) :: denseDesc

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> Change in density matrix during this SCC step for rangesep
    real(dp), pointer, intent(in) :: deltaRhoOutSqr(:,:,:)

    !> derivatives of energy wrt to atomic positions
    real(dp), intent(out) :: derivs(:,:)

    !> sparse density matrix
    real(dp), allocatable :: rho(:)

    !> sparse energy weighted density matrix
    real(dp), allocatable :: ERho(:)

    !> dense density matrix
    real(dp), allocatable :: rhoSqr(:,:)

    !> dense energy weighted density matrix
    real(dp), allocatable :: ERhoSqr(:,:)

    ! Locals
    real(dp), allocatable :: dQ(:,:,:)
    logical :: tSccCalc
    integer :: sparseSize, nAtom, nOrb


    tSccCalc = allocated(sccCalc)
    nAtom = size(derivs, dim=2)
    nOrb = size(eigVecs, dim=1)
    sparseSize = size(over)

    ! First, we need to create density matrix and energy-weighted density matrix
    !   corresponding to the (active) frontier orbitals

    ! Allocate arrays for the matrices
    allocate(rho(sparseSize))
    allocate(ERho(sparseSize))
    allocate(rhoSqr(nOrb,nOrb))
    allocate(ERhoSqr(nOrb,nOrb))

    ! Create the matrices in square format
    call fmo_sp_density_matrix_real(rhoSqr, eigVecs, filling,&
        & neighbourList%iNeighbour, nNeighbourSK, orb, denseDesc%iAtomStart, img2CentCell)
    call fmo_sp_energy_density_matrix_real(ERhoSqr, eigVecs, filling, eigen,&
        & neighbourList%iNeighbour, nNeighbourSK, orb, denseDesc%iAtomStart, img2CentCell)

  ! write(*,*) "RHOSQR"
  ! write(*,'(66F9.5)') rhoSqr
  ! write(*,*) "RHOSQR END"
  ! write(*,*) "ERHOSQR"
  ! write(*,'(66F9.5)') ERhoSqr
  ! write(*,*) "ERHOSQR END"

    ! Pack them
    rho(:) = 0._dp
    call packHS(rho, rhoSqr, neighbourlist%iNeighbour, nNeighbourSK, orb%mOrb,&
        & denseDesc%iAtomStart, iSparseStart, img2CentCell)
    ERho(:) = 0._dp
    call packHS(ERho, ERhoSqr, neighbourlist%iNeighbour, nNeighbourSK, orb%mOrb,&
        & denseDesc%iAtomStart, iSparseStart, img2CentCell)

    derivs(:,:) = 0.0_dp

    if (.not. (tSccCalc)) then
      ! No external or internal potentials
      call fmo_derivative_nonscc(env, derivs, nonSccDeriv, rho, ERho, skHamCont, skOverCont, coord,&
          & species, neighbourList%iNeighbour, nNeighbourSK, img2CentCell, iSparseStart, orb)
    else
      call fmo_derivative_shift(env, derivs, nonSccDeriv, rho, ERho, skHamCont, skOverCont, coord,&
          & species, neighbourList%iNeighbour, nNeighbourSK, img2CentCell, iSparseStart, orb,&
          & potential%intBlock)

      call sccCalc%addForceDc(env, derivs, species, neighbourList%iNeighbour, img2CentCell)

      if (allocated(thirdOrd)) then
        call thirdOrd%addGradientDc(neighbourList, species, coord, img2CentCell, derivs)
      end if

      if (allocated(qDepExtPot)) then
        allocate(dQ(orb%mShell, nAtom, size(qOutput, dim=3)))
        call getChargePerShell(qOutput, orb, species, dQ, qRef=q0)
        call qDepExtPot%addGradientDc(sum(dQ(:,:,1), dim=1), dQ(:,:,1), derivs)
      end if
    end if

    if (allocated(rangeSep)) then
      call unpackHS(SSqrReal, over, neighbourList%iNeighbour, nNeighbourSK, denseDesc%iAtomStart,&
          & iSparseStart, img2CentCell)
      call rangeSep%addLRGradients(derivs, nonSccDeriv, deltaRhoOutSqr, skOverCont, coord,&
          & species, orb, denseDesc%iAtomStart, SSqrReal, neighbourList%iNeighbour, nNeighbourSK)
    end if

    deallocate(rho)
    deallocate(ERho)
    deallocate(rhoSqr)
    deallocate(ERhoSqr)

  end subroutine fmoGradients


  ! COPIED AND MODIFIED FROM:
  !   forces.F90, derivative_shift -> derivative_nonSCC -> derivativeNonSccEuclidian
  !> The non-SCC electronic force contribution for all atoms from the matrix derivatives and the
  !> density and energy-density matrices
  subroutine fmo_derivative_NonScc(env, deriv, derivator, DM, EDM, skHamCont, skOverCont, coords,&
      & species, iNeighbour, nNeighbourSK, img2CentCell, iPair, orb)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> x,y,z derivatives for each real atom in the system
    real(dp), intent(out) :: deriv(:,:)

    !> Differentiatior for the non-scc components
    class(TNonSccDiff), intent(in) :: derivator

    !> density matrix in packed format
    real(dp), intent(in) :: DM(:)

    !> energy-weighted density matrix in packed format
    real(dp), intent(in) :: EDM(:)

    !> Container for SK Hamiltonian integrals
    type(TSlakoCont), intent(in) :: skHamCont

    !> Container for SK overlap integrals
    type(TSlakoCont), intent(in) :: skOverCont

    !> list of all atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> list of all atomic species
    integer, intent(in) :: species(:)

    !> neighbour list for atoms
    integer, intent(in) :: iNeighbour(0:,:)

    !> number of neighbours of each atom
    integer, intent(in) :: nNeighbourSK(:)

    !> indexing array for periodic image atoms
    integer, intent(in) :: img2CentCell(:)

    !> indexing array for the Hamiltonian
    integer, intent(in) :: iPair(0:,:)

    !> Information about the shells and orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    integer :: iOrig, ii, nAtom, iNeigh, iAtom1, iAtom2, iAtom2f, nOrb1, nOrb2, iAtFirst, iAtLast
    real(dp) :: sqrDMTmp(orb%mOrb,orb%mOrb), sqrEDMTmp(orb%mOrb,orb%mOrb)
    real(dp) :: hPrimeTmp(orb%mOrb,orb%mOrb,3), sPrimeTmp(orb%mOrb,orb%mOrb,3)

    @:ASSERT(size(deriv,dim=1) == 3)

    nAtom = size(orb%nOrbAtom)
    deriv(:,:) = 0.0_dp

    ! remove MPI support
  ! call distributeRangeInChunks(env, 1, nAtom, iAtFirst, iAtLast)
    iAtFirst = 1
    iAtLast = nAtom

    !$OMP PARALLEL DO PRIVATE(nOrb1, iNeigh, iAtom2, iAtom2f, nOrb2, iOrig, sqrDMTmp, sqrEDMTmp,&
    !$OMP& hPrimeTmp, sPrimeTmp,ii) DEFAULT(SHARED) SCHEDULE(RUNTIME) REDUCTION(+:deriv)
    do iAtom1 = iAtFirst, iAtLast
      nOrb1 = orb%nOrbAtom(iAtom1)
      !! loop from 1 as no contribution from the atom itself
      do iNeigh = 1, nNeighbourSK(iAtom1)
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        if (iAtom1 /= iAtom2f) then
          nOrb2 = orb%nOrbAtom(iAtom2f)
          iOrig = iPair(iNeigh,iAtom1)
          sqrDMTmp(:,:) = 0.0_dp
          sqrEDMTmp(:,:) = 0.0_dp
          hPrimeTmp(:,:,:) = 0.0_dp
          sPrimeTmp(:,:,:) = 0.0_dp
          sqrDMTmp(1:nOrb2,1:nOrb1) = reshape(DM(iOrig+1:iOrig+nOrb1*nOrb2), (/nOrb2,nOrb1/))
          sqrEDMTmp(1:nOrb2,1:nOrb1) = reshape(EDM(iOrig+1:iOrig+nOrb1*nOrb2), (/nOrb2,nOrb1/))
          call derivator%getFirstDeriv(hPrimeTmp, skHamCont, coords, species, iAtom1, iAtom2, orb)
          call derivator%getFirstDeriv(sPrimeTmp, skOverCont, coords, species, iAtom1, iAtom2, orb)
          ! note factor of 2 for implicit summation over lower triangle of density matrix:
          do ii = 1, 3
            deriv(ii,iAtom1) = deriv(ii,iAtom1)&
                & + sum(sqrDMTmp(1:nOrb2,1:nOrb1) * 2.0_dp*hPrimeTmp(1:nOrb2,1:nOrb1,ii))&
                & - sum(sqrEDMTmp(1:nOrb2,1:nOrb1) * 2.0_dp*sPrimeTmp(1:nOrb2,1:nOrb1,ii))
          end do
          ! Add contribution to the force from atom 1 onto atom 2f using the symmetry in the
          ! blocks, and note that the skew symmetry in the derivatives is being used
          do ii = 1, 3
            deriv(ii,iAtom2f) = deriv(ii,iAtom2f)&
                & - sum(sqrDMTmp(1:nOrb2,1:nOrb1) * 2.0_dp*hPrimeTmp(1:nOrb2,1:nOrb1,ii))&
                & + sum(sqrEDMTmp(1:nOrb2,1:nOrb1) * 2.0_dp*sPrimeTmp(1:nOrb2,1:nOrb1,ii))
          end do
        end if
      end do
    end do
    !$OMP END PARALLEL DO

    ! remove MPI support
  ! call assembleChunks(env, deriv)

  end subroutine fmo_derivative_NonScc


  ! COPIED AND MODIFIED FROM:
  !   forces.F90, derivative_shift -> derivative_block -> derivative_blockEuclidian
  !> The SCC and spin electronic force contribution for all atoms from the matrix derivatives, self
  !> consistent potential and the density and energy-density matrices
  subroutine fmo_derivative_shift(env, deriv, derivator, DM, EDM, skHamCont, skOverCont,&
      & coords, species, iNeighbour, nNeighbourSK, img2CentCell, iPair, orb, shift)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> x,y,z derivatives for each real atom in the system
    real(dp), intent(out) :: deriv(:,:)

    !> Differentiatior for the non-scc components
    class(TNonSccDiff), intent(in) :: derivator

    !> density matrix in packed format
    real(dp), intent(in) :: DM(:)

    !> energy-weighted density matrix in packed format
    real(dp), intent(in) :: EDM(:)

    !> Container for SK Hamiltonian integrals
    type(TSlakoCont) :: skHamCont

    !> Container for SK overlap integrals
    type(TSlakoCont) :: skOverCont

    !> list of all atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> list of all atomic species
    integer, intent(in) :: species(:)

    !> neighbour list for atoms
    integer, intent(in) :: iNeighbour(0:,:)

    !> number of neighbours of each atom
    integer, intent(in) :: nNeighbourSK(:)

    !> indexing array for periodic image atoms
    integer, intent(in) :: img2CentCell(:)

    !> indexing array for the Hamiltonian
    integer, intent(in) :: iPair(0:,:)

    !> Information about the shells and orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    !> block shift from the potential
    real(dp), intent(in) :: shift(:,:,:,:)

    integer :: iOrig, ii, nAtom, iNeigh, iAtom1, iAtom2, iAtom2f, iSp1, iSp2
    integer :: nOrb1, nOrb2, iAtFirst, iAtLast

    real(dp) :: sqrDMTmp(orb%mOrb,orb%mOrb), sqrEDMTmp(orb%mOrb,orb%mOrb)
    real(dp) :: shiftSprime(orb%mOrb,orb%mOrb)
    real(dp) :: hPrimeTmp(orb%mOrb,orb%mOrb,3), sPrimeTmp(orb%mOrb,orb%mOrb,3)
    real(dp) :: derivTmp(3)

    nAtom = size(orb%nOrbAtom)
    ! nSpin = size(shift,dim=4)
    @:ASSERT(size(shift,dim=4) == 1)
    ! ^ not 2 or 4
    @:ASSERT(size(deriv,dim=1) == 3)
    @:ASSERT(size(deriv,dim=2)==nAtom)
    @:ASSERT(size(DM)==size(EDM))
    @:ASSERT(size(shift,dim=1)==orb%mOrb)
    @:ASSERT(size(shift,dim=2)==orb%mOrb)
    @:ASSERT(size(shift,dim=3)==nAtom)

    deriv(:,:) = 0.0_dp

    ! remove MPI support
  ! call distributeRangeInChunks(env, 1, nAtom, iAtFirst, iAtLast)
    iAtFirst = 1
    iAtLast = nAtom

    !$OMP PARALLEL DO PRIVATE(iAtom1,iSp1,nOrb1,iNeigh,iAtom2,iAtom2f,iSp2,nOrb2,iOrig,sqrDMTmp, &
    !$OMP& sqrEDMTmp,hPrimeTmp,sPrimeTmp,derivTmp,shiftSprime,ii) DEFAULT(SHARED) &
    !$OMP& SCHEDULE(RUNTIME) REDUCTION(+:deriv)
    do iAtom1 = iAtFirst, iAtLast
      iSp1 = species(iAtom1)
      nOrb1 = orb%nOrbSpecies(iSp1)
      do iNeigh = 1, nNeighbourSK(iAtom1)
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        iSp2 = species(iAtom2f)
        if (iAtom1 /= iAtom2f) then
          nOrb2 = orb%nOrbSpecies(iSp2)
          iOrig = iPair(iNeigh,iAtom1) + 1
          sqrDMTmp (1:nOrb2,1:nOrb1) = reshape(DM (iOrig:iOrig+nOrb1*nOrb2-1),(/nOrb2,nOrb1/))
          sqrEDMTmp(1:nOrb2,1:nOrb1) = reshape(EDM(iOrig:iOrig+nOrb1*nOrb2-1),(/nOrb2,nOrb1/))
          call derivator%getFirstDeriv(hPrimeTmp, skHamCont, coords, species, iAtom1, iAtom2, orb)
          call derivator%getFirstDeriv(sPrimeTmp, skOverCont, coords, species, iAtom1, iAtom2,&
              & orb)

          derivTmp(:) = 0.0_dp
          ! note factor of 2 for implicit summation over lower triangle of density matrix:
          do ii = 1, 3
            derivTmp(ii) = 2.0_dp * (&
                & sum(sqrDMTmp(1:nOrb2,1:nOrb1)*hPrimeTmp(1:nOrb2,1:nOrb1,ii))&
                & - sum(sqrEDMTmp(1:nOrb2,1:nOrb1)*sPrimeTmp(1:nOrb2,1:nOrb1,ii)))
          end do

          do ii = 1, 3
            shiftSprime(1:nOrb2,1:nOrb1) = 0.5_dp * (&
                & matmul(sPrimeTmp(1:nOrb2,1:nOrb1,ii), shift(1:nOrb1,1:nOrb1,iAtom1,1) )&
                & + matmul(shift(1:nOrb2,1:nOrb2,iAtom2f,1), sPrimeTmp(1:nOrb2,1:nOrb1,ii)))
            ! again factor of 2 from lower triangle, cf published force expressions for SCC:
            derivTmp(ii) = derivTmp(ii) + 2.0_dp * ( sum(shiftSprime(1:nOrb2,1:nOrb1) *&
                & reshape(DM(iOrig:iOrig+nOrb1*nOrb2-1),(/nOrb2,nOrb1/)) ) )
          end do

          ! forces from atom 1 on atom 2f and 2f onto 1
          deriv(:,iAtom1) = deriv(:,iAtom1) + derivTmp(:)
          deriv(:,iAtom2f) = deriv(:,iAtom2f) - derivTmp(:)

        end if
      enddo
    enddo
    !$OMP END PARALLEL DO

    ! remove MPI support
  ! call assembleChunks(env, deriv)

  end subroutine fmo_derivative_shift


  ! COPIED AND MODIFIED FROM:
  !   densitymatrix.F90, makeDensityMatrix -> sp_density_matrix_real
  !> Make a regular density matrix for the real wave-function case
  subroutine fmo_sp_density_matrix_real(dm, eigenvecs, filling, iNeighbour, nNeighbourSK, orb,&
      & iAtomStart, img2CentCell)

    !> the resulting nOrb*nOrb density matrix with only the elements of interest
    !> calculated, instead of the whole dm
    real(dp), intent(out) :: dm(:,:)

    !> the eigenvectors of the system
    real(dp), intent(in) :: eigenvecs(:,:)

    !> the occupation numbers of the orbitals
    real(dp), intent(in) :: filling(:)

    !> Neighbour list for each atom (First index from 0!)
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for each atom (incl. itself).
    integer, intent(in) :: nNeighbourSK(:)

    !> Information about the orbitals of the atoms
    type(TOrbitals), intent(in) :: orb

    !> Atom offset for the squared Hamiltonian
    integer, intent(in) :: iAtomStart(:)

    !> Atomic mapping indexes.
    integer, intent(in) :: img2CentCell(:)

    integer :: iAt1, iNeigh1, nOrb1, nOrb2, jj
    integer :: nAtom, nLevels, start1, start2
    integer, allocatable :: inCellNeighbour(:,:)
    integer, allocatable :: nInCellNeighbour(:)
    real(dp), allocatable :: tmpEigen(:,:)

    @:ASSERT(all(shape(eigenvecs) == shape(dm)))
    @:ASSERT(size(eigenvecs,dim=1) == size(eigenvecs,dim=2))
    @:ASSERT(size(eigenvecs,dim=1) == size(filling))

    allocate(inCellNeighbour(0:size(iNeighbour,dim=1),size(iNeighbour,dim=2)))
    allocate(nInCellNeighbour(size(iNeighbour,dim=2)))

    nAtom = size(orb%nOrbAtom)
    dm(:,:) = 0.0_dp

    inCellNeighbour(:,:) = 0
    nInCellNeighbour(:) = 0

    ! do not restrict the procedure to occupied levels
    nLevels = size(filling)
    ! if (abs(filling(ii)) >= epsilon(1.0_dp)) then

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
        tmpEigen(jj,1:nOrb1) = filling(jj) * tmpEigen(jj,1:nOrb1)
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

  end subroutine fmo_sp_density_matrix_real


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


  !> Calculates the gradients stemming from interaction between pairs of different fragments
  subroutine fmoGradientsOffdiag(env, sccCalc, nonSccDeriv, nAtom1, eigVecs1, eigen1, filling1, &
      & nAtom2, eigVecs2, eigen2, filling2, qOutput, q0,&
      & skHamCont, skOverCont, neighbourList, nNeighbourSK, species, img2CentCell, orb, potential,&
      & coord, thirdOrd, qDepExtPot, rangeSep, SSqrReal, over, denseDesc, iSparseStart,&
      & deltaRhoOutSqr, derivs)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> SCC module internal variables
    type(TScc), allocatable, intent(in) :: sccCalc

    !> method for calculating derivatives of S and H0
    type(TNonSccDiff), intent(in) :: nonSccDeriv

    !> number of atoms in fragment 1
    integer, intent(in) :: nAtom1

    !> the eigenvectors of fragment 1
    real(dp), intent(in) :: eigVecs1(:,:)

    !> the eigenvalues of fragment 1
    real(dp), intent(in) :: eigen1(:)

    !> filling of the FOs in fragment 1 -- active orbitals only (1 for a single FMO, or so...)
    real(dp), intent(in) :: filling1(:)

    !> number of atoms in fragment 2
    integer, intent(in) :: nAtom2

    !> the eigenvectors of fragment 2
    real(dp), intent(in) :: eigVecs2(:,:)

    !> the eigenvalues of fragment 2
    real(dp), intent(in) :: eigen2(:)

    !> filling of the FOs in fragment 2
    real(dp), intent(in) :: filling2(:)

    !> electron populations (may be unallocated for non-scc case)
    real(dp), allocatable, intent(in) :: qOutput(:,:,:)

    !> reference atomic charges (may be unallocated for non-scc case)
    real(dp), allocatable, intent(in) :: q0(:,:,:)

    !> non-SCC hamiltonian information
    type(TSlakoCont), intent(in) :: skHamCont

    !> overlap information
    type(TSlakoCont), intent(in) :: skOverCont

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> species of all atoms in the system
    integer, intent(in) :: species(:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !>  potential acting on the system
    type(TPotentials), intent(in) :: potential

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> Is 3rd order SCC being used
    type(TThirdOrder), intent(inout), allocatable :: thirdOrd

    !> Population dependant external potential
    type(TQDepExtPotProxy), intent(inout), allocatable :: qDepExtPot

    !> Data from rangeseparated calculations
    type(TRangeSepFunc), intent(inout), allocatable :: rangeSep

    !> dense overlap matrix, required for rangeSep
    real(dp), intent(inout), allocatable :: SSqrReal(:,:)

    !> sparse overlap matrix, required for rangeSep
    real(dp), intent(in) :: over(:)

    !> Dense matrix descriptor,required for rangeSep
    type(TDenseDescr), intent(in) :: denseDesc

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> Change in density matrix during this SCC step for rangesep
    real(dp), pointer, intent(in) :: deltaRhoOutSqr(:,:,:)

    !> derivatives of energy wrt to atomic positions
    real(dp), intent(out) :: derivs(:,:)

    !> sparse density matrix
    real(dp), allocatable :: rho(:)

    !> sparse energy weighted density matrix
    real(dp), allocatable :: ERho(:)

    !> dense density matrix
    real(dp), allocatable :: rhoSqr(:,:)

    !> dense energy weighted density matrix
    real(dp), allocatable :: ERhoSqr(:,:)

    ! Locals
    real(dp), allocatable :: dQ(:,:,:)
    logical :: tSccCalc
    integer :: sparseSize, nAtom, nOrb1, nOrb2, nOrb
    real(dp), allocatable :: eigen(:), filling(:), eigVecs(:,:)



    tSccCalc = allocated(sccCalc)
  ! nAtom = size(derivs, dim=2)
    nAtom = nAtom1 + nAtom2
    nOrb1 = size(eigVecs1, dim=1)
    nOrb2 = size(eigVecs2, dim=1)
    nOrb = nOrb1 + nOrb2
    sparseSize = size(over)

    ! First, we need to create density matrix and energy-weighted density matrix
    !   corresponding to the (active) frontier orbitals

    ! Allocate arrays for the matrices
    allocate(eigen(nOrb))
    allocate(filling(nOrb))
    allocate(rho(sparseSize))
    allocate(ERho(sparseSize))
    allocate(rhoSqr(nOrb,nOrb))
    allocate(ERhoSqr(nOrb,nOrb))
    allocate(eigVecs(nOrb,nOrb))

    eigVecs(1:nOrb1,1:nOrb1) = eigVecs1
    eigVecs(nOrb1+1:nOrb2,nOrb1+1:nOrb2) = eigVecs2

    ! Create the matrices in square format
    call fmo_sp_density_matrix_real(rhoSqr, eigVecs, filling,&
        & neighbourList%iNeighbour, nNeighbourSK, orb, denseDesc%iAtomStart, img2CentCell)
    call fmo_sp_energy_density_matrix_real(ERhoSqr, eigVecs, filling, eigen,&
        & neighbourList%iNeighbour, nNeighbourSK, orb, denseDesc%iAtomStart, img2CentCell)

  ! write(*,*) "RHOSQR"
  ! write(*,'(66F9.5)') rhoSqr
  ! write(*,*) "RHOSQR END"
  ! write(*,*) "ERHOSQR"
  ! write(*,'(66F9.5)') ERhoSqr
  ! write(*,*) "ERHOSQR END"

    ! Pack them
    rho(:) = 0._dp
    call packHS(rho, rhoSqr, neighbourlist%iNeighbour, nNeighbourSK, orb%mOrb,&
        & denseDesc%iAtomStart, iSparseStart, img2CentCell)
    ERho(:) = 0._dp
    call packHS(ERho, ERhoSqr, neighbourlist%iNeighbour, nNeighbourSK, orb%mOrb,&
        & denseDesc%iAtomStart, iSparseStart, img2CentCell)

    derivs(:,:) = 0.0_dp

    if (.not. (tSccCalc)) then
      ! No external or internal potentials
      call fmo_derivative_nonscc(env, derivs, nonSccDeriv, rho, ERho, skHamCont, skOverCont, coord,&
          & species, neighbourList%iNeighbour, nNeighbourSK, img2CentCell, iSparseStart, orb)
    else
      call fmo_derivative_shift(env, derivs, nonSccDeriv, rho, ERho, skHamCont, skOverCont, coord,&
          & species, neighbourList%iNeighbour, nNeighbourSK, img2CentCell, iSparseStart, orb,&
          & potential%intBlock)

      call sccCalc%addForceDc(env, derivs, species, neighbourList%iNeighbour, img2CentCell)

      if (allocated(thirdOrd)) then
        call thirdOrd%addGradientDc(neighbourList, species, coord, img2CentCell, derivs)
      end if

      if (allocated(qDepExtPot)) then
        allocate(dQ(orb%mShell, nAtom, size(qOutput, dim=3)))
        call getChargePerShell(qOutput, orb, species, dQ, qRef=q0)
        call qDepExtPot%addGradientDc(sum(dQ(:,:,1), dim=1), dQ(:,:,1), derivs)
      end if
    end if

    if (allocated(rangeSep)) then
      call unpackHS(SSqrReal, over, neighbourList%iNeighbour, nNeighbourSK, denseDesc%iAtomStart,&
          & iSparseStart, img2CentCell)
      call rangeSep%addLRGradients(derivs, nonSccDeriv, deltaRhoOutSqr, skOverCont, coord,&
          & species, orb, denseDesc%iAtomStart, SSqrReal, neighbourList%iNeighbour, nNeighbourSK)
    end if

    deallocate(rho)
    deallocate(ERho)
    deallocate(rhoSqr)
    deallocate(ERhoSqr)

  end subroutine fmoGradientsOffdiag

! !  USE PARTS OF THIS PROCEDURE TO CREATE AN INTERFACE TO OFFDIAG GRADIENTS!
! !> Run the 2nd phase of a fragment-molecular orbital calculation (FMO)
! subroutine processGeometryPhase2(this, env, nSite, ptrsPhase1, filling, derivs)

!   !> Global variables
!   type(TDftbPlusMain), intent(inout) :: this

!   !> Environment settings
!   type(TEnvironment), intent(inout) :: env

!   !> Number of sites / fragments
!   integer, intent(in) :: nSite

!   !> Array of structures containing pointers to data structures in phase1
!   type(TPointersToPhase1), allocatable, intent(in) :: ptrsPhase1(:)

!   !> filling of the FOs in all fragments
!   real(dp), intent(in) :: filling(:)

!   !> derivatives of energy wrt to atomic positions
!   real(dp), intent(out) :: derivs(:,:)

!   integer, allocatable :: iSiteAtomStart(:), iSiteOrbStart(:)

!   !> dense hamiltonian matrix
!   real(dp), allocatable :: HSqrReal(:,:)

!   !> dense overlap matrix
!   real(dp), allocatable :: SSqrReal(:,:)

!   !> temporary array for charges
!   real(dp), allocatable :: dQ(:,:,:)

!   !> Hamiltonian and overlap matrices in the FMO basis set
!   real(dp), allocatable :: Tij(:,:), Sij(:,:)

!   integer :: iSite, jSite, nAtom, iAtom, iAtomBeg, iAtomEnd, iOrb, jOrb, iOrbBeg, iOrbEnd
!   integer :: nFO, iFO, jFO, iAO, jAO
!   integer, allocatable :: indFO(:)

!   ! number of atoms and orbitals in each site and totals
!   allocate(iSiteAtomStart(nSite+1))
!   allocate(iSiteOrbStart(nSite+1))
!   iSiteAtomStart(1) = 1
!   iSiteOrbStart(1) = 1
!   do iSite = 1, nSite
!     iSiteAtomStart(iSite+1) = iSiteAtomStart(iSite) + ptrsPhase1(iSite)%nAtom
!     iSiteOrbStart(iSite+1) = iSiteOrbStart(iSite) + ptrsPhase1(iSite)%nOrb
!   end do

!   ! total number of fragment orbitals to consider = dimension of resulting Hamiltonian
!   nFO = 0
!   do iSite = 1, nSite
!     nFO = nFO + ptrsPhase1(iSite)%nFO
!   end do
!   ! indFO(i) is the index into Tij & Sij at which the orbitals on fragment (i) start
!   allocate(indFO(iSite))
!   indFO(1) = 1
!   do iSite = 1, nSite-1
!     indFO(iSite+1) = indFO(iSite) + ptrsPhase1(iSite)%nFO
!   end do

!   call env%globalTimer%startTimer(globalTimers%preSccInit)

!   ! CHARGE-INDEPENDENT MATRICES -- CONVERT TO DENSE
!   allocate(HSqrReal(this%nOrb,this%nOrb))
!   allocate(SSqrReal(this%nOrb,this%nOrb))
!   call unpackHS(HSqrReal, this%H0, this%neighbourList%iNeighbour, this%nNeighbourSK,&
!       & this%denseDesc%iAtomStart, this%iSparseStart, this%img2CentCell)
!   call unpackHS(SSqrReal, this%over, this%neighbourList%iNeighbour, this%nNeighbourSK,&
!       & this%denseDesc%iAtomStart, this%iSparseStart, this%img2CentCell)

!   ! CHARGE-INDEPENDENT MATRICES -- INSERT DIAGONAL BLOCKS FROM PHASE 1
!   do iSite = 1, nSite
!     iOrbBeg = iSiteOrbStart(iSite)
!     iOrbEnd = iSiteOrbStart(iSite+1) - 1
!     HSqrReal(iOrbBeg:iOrbEnd,iOrbBeg:iOrbEnd) = ptrsPhase1(iSite)%denseH0(:,:)
!     SSqrReal(iOrbBeg:iOrbEnd,iOrbBeg:iOrbEnd) = ptrsPhase1(iSite)%denseOver(:,:)
!   end do

!   ! CHARGE-INDEPENDENT MATRICES -- CONVERT TO SPARSE
!   this%h0 = 0._dp
!   call packHS(this%h0, HSqrReal, this%neighbourList%iNeighbour, this%nNeighbourSK,&
!       & this%orb%mOrb, this%denseDesc%iAtomStart, this%iSparseStart, this%img2CentCell)
!   this%over = 0._dp
!   call packHS(this%over, SSqrReal, this%neighbourList%iNeighbour, this%nNeighbourSK,&
!       & this%orb%mOrb, this%denseDesc%iAtomStart, this%iSparseStart, this%img2CentCell)

!   ! EXTERNAL POTENTIALS -- UNIFORM ELECTRIC FIELD
!   call resetExternalPotentials(this%refExtPot, this%potential)
!   if (this%tEField) then
!     nAtom = size(this%nNeighbourSK)
!     this%EField(:) = this%EFieldStrength * this%EfieldVector
!     do iAtom = 1, nAtom
!       this%potential%extAtom(iAtom,1) = this%potential%extAtom(iAtom,1)&
!           & + dot_product(this%coord(:, iAtom), this%EField)
!     end do
!     this%potential%extGrad(:,:) = this%potential%extGrad + spread(this%EField, 2, nAtom)
!   else
!     this%EField(:) = 0.0_dp
!   end if
!   call mergeExternalPotentials(this%orb, this%species, this%potential)

!   call env%globalTimer%stopTimer(globalTimers%preSccInit)

!   call env%globalTimer%startTimer(globalTimers%scc)

!   call resetInternalPotentials(.false., this%xi, this%orb, this%species, this%potential)

!   ! GET ATOM CHARGES FROM THE INDIVIDUAL FRAGMENTS / SITES
!   ! qInput(mOrb, nAtom, nSpin); nSpin == 1
!   do iSite = 1, nSite
!     iAtomBeg = iSiteAtomStart(iSite)
!     iAtomEnd = iSiteAtomStart(iSite+1)-1
!     this%qInput(:,iAtomBeg:iAtomEnd,1) = ptrsPhase1(iSite)%qInput(:,:,1)
!   end do

!   if (this%tSccCalc) then

!     ! GET SHELL CHARGES FROM THE INDIVIDUAL FRAGMENTS / SITES
!     call getChargePerShell(this%qInput, this%orb, this%species, this%chargePerShell)

!     call addChargePotentials(env, this%sccCalc, this%qInput, this%q0, this%chargePerShell,&
!         & this%orb, this%species, this%neighbourList, this%img2CentCell, this%spinW,&
!         & this%solvation, this%thirdOrd, this%potential, this%electrostatics, this%tPoisson,&
!         & this%tUpload, this%shiftPerLUp, this%dispersion)

!   end if

!   ! All potentials are added up into intBlock
!   this%potential%intBlock = this%potential%intBlock + this%potential%extBlock

!   ! DO WE NEED THIS?
!   if (allocated(this%qDepExtPot)) then
!     allocate(dQ(this%orb%mShell, this%nAtom, this%nSpin))
!     call getChargePerShell(this%qInput, this%orb, this%species, dQ, qRef=this%q0)
!     call this%qDepExtPot%addPotential(sum(dQ(:,:,1), dim=1), dQ(:,:,1), this%orb,&
!         & this%species, this%potential%intBlock)
!   end if

!   ! SELF-CONSISTENT / CHARGE-DEPENDENT HAMILTONIAN
!   call getSccHamiltonian(this%H0, this%over, this%nNeighbourSK, this%neighbourList,&
!       & this%species, this%orb, this%iSparseStart, this%img2CentCell, this%potential,&
!       & allocated(this%reks), this%ham, this%iHam)
!   call unpackHS(HSqrReal, this%ham(:,1), this%neighbourList%iNeighbour, this%nNeighbourSK,&
!       & this%denseDesc%iAtomStart, this%iSparseStart, this%img2CentCell)

!   call env%globalTimer%stopTimer(globalTimers%scc)

!   call env%globalTimer%startTimer(globalTimers%postSCC)

!   ! fill the other triangle of both matrices
!   do iAO = 2, this%nOrb
!     do jAO = 1, iAO
!       HSqrReal(jAO,iAO) = HSqrReal(iAO,jAO)
!       SSqrReal(jAO,iAO) = SSqrReal(iAO,jAO)
!     end do
!   end do

!   ! CALCULATE THE FMO HAMILTONIAN AND OVERLAP
!   allocate(Tij(nFO,nFO))
!   allocate(Sij(nFO,nFO))
!   Tij = 0._dp
!   Sij = 0._dp

!   do iSite = 1, nSite
!     do jSite = iSite, nSite
!       do iOrb = 1, ptrsPhase1(iSite)%nFO
!         iFO = indFO(iSite) + iOrb - 1
!         do jOrb = 1, ptrsPhase1(iSite)%nFO
!           jFO = indFO(jSite) + jOrb - 1

!           Tij(iFO,jFO) = 0._dp
!           Sij(iFO,jFO) = 0._dp

!           do iAO = 1, ptrsPhase1(iSite)%nOrb ! iSiteOrbStart(iSite), iSiteOrbStart(iSite+1)-1
!             do jAO = 1, ptrsPhase1(jSite)%nOrb ! iSiteOrbStart(jSite), iSiteOrbStart(jSite+1)-1
!               Tij(iFO,jFO) = Tij(iFO,jFO)&
!                 ! & + ptrsPhase1(iSite)%eigVec(iAO,iSiteOrbStart(iSite)+ptrsPhase1(iSite)%iHOMO-iOrb)&
!                   & + ptrsPhase1(iSite)%eigVec(iAO,ptrsPhase1(iSite)%iHOMO+1-iOrb)&
!                   & * HSqrReal(iSiteOrbStart(iSite)+iAO-1,iSiteOrbStart(jSite)+jAO-1)&
!                 ! & * ptrsPhase1(jSite)%eigVec(jAO,iSiteOrbStart(jSite)+ptrsPhase1(jSite)%iHOMO-jOrb)
!                   & * ptrsPhase1(jSite)%eigVec(jAO,ptrsPhase1(jSite)%iHOMO+1-jOrb)
!               Sij(iFO,jFO) = Sij(iFO,jFO)&
!                 ! & + ptrsPhase1(iSite)%eigVec(iAO,iSiteOrbStart(iSite)+ptrsPhase1(iSite)%iHOMO-iOrb)&
!                   & + ptrsPhase1(iSite)%eigVec(iAO,ptrsPhase1(iSite)%iHOMO+1-iOrb)&
!                   & * SSqrReal(iSiteOrbStart(iSite)+iAO-1,iSiteOrbStart(jSite)+jAO-1)&
!                 ! & * ptrsPhase1(jSite)%eigVec(jAO,iSiteOrbStart(jSite)+ptrsPhase1(jSite)%iHOMO-jOrb)
!                   & * ptrsPhase1(jSite)%eigVec(jAO,ptrsPhase1(jSite)%iHOMO+1-jOrb)
!             end do
!           end do
!         end do ! jOrb
!       end do ! iOrb
!     end do ! jSite
!   end do ! iSite

!   ! COPY TO THE OTHER TRIANGLE OF THE MATRIX
!   do iFO = 1, nFO
!     do jFO = iFO+1, nFO
!       Tij(jFO,iFO) = Tij(iFO,jFO)
!       Sij(jFO,iFO) = Sij(iFO,jFO)
!     end do
!   end do

! end subroutine processGeometryPhase2

end module dftbp_fmogradient
