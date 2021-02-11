!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> The main routines for DFTB+
module dftbp_fmo
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
  use dftbp_externalcharges
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

  implicit none
  private

  public :: TPointersToPhase1
  public :: processGeometryPhase2

  ! there will be an array of these structures, one for each fragment/site
  type :: TPointersToPhase1

    ! number of atoms
    integer :: nAtom

    ! number of orbitals
    integer :: nOrb

    ! number of frontier / fragment orbitals to consider
    integer :: nFO

    ! which orbital is HOMO/LUMO?
    integer :: iHOMO

    ! sparse matrix, overlap
    real(dp), allocatable :: denseOver(:,:)

    ! sparse matrix, charge-independent Hamiltonian
    real(dp), allocatable :: denseH0(:,:)

    ! dense matrix, eigenvectors
    real(dp), allocatable :: eigVec(:,:)

    ! 1D array, eigenvalues (orbital energies)
    real(dp), allocatable :: eigVal(:)

    ! charge per atom (mOrb, atom, spin channel); spin channel == 1
    real(dp), allocatable :: qInput(:,:,:)

    ! charge per atomic shell (shell, atom, spin channel); spin channel == 1
    real(dp), allocatable :: chargePerShell(:,:,:)

    ! an entry for gradients?

  end type TPointersToPhase1

contains

  !> Run the 2nd phase of a fragment-molecular orbital calculation (FMO)
  subroutine processGeometryPhase2(this, env, nSite, ptrsPhase1, TijOrtho)

    !> Global variables
    type(TDftbPlusMain), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Number of sites / fragments
    integer, intent(in) :: nSite

    !> Array of structures containing pointers to data structures in phase1
    type(TPointersToPhase1), allocatable, intent(in) :: ptrsPhase1(:)

    !> Orthogonalized FMO Hamiltonian
    real(dp), allocatable, intent(out) :: TijOrtho(:,:)

    integer, allocatable :: iSiteAtomStart(:), iSiteOrbStart(:)

    !> dense hamiltonian matrix
    real(dp), allocatable :: HSqrReal(:,:)

    !> dense overlap matrix
    real(dp), allocatable :: SSqrReal(:,:)

    !> temporary array for charges
    real(dp), allocatable :: dQ(:,:,:)

    !> Hamiltonian and overlap matrices in the FMO basis set
    real(dp), allocatable :: Tij(:,:), Sij(:,:)

    integer :: iSite, jSite, nAtom, iAtom, iAtomBeg, iAtomEnd, iOrb, jOrb, iOrbBeg, iOrbEnd
    integer :: nFO, iFO, jFO, iAO, jAO
    integer, allocatable :: indFO(:)

    ! number of atoms and orbitals in each site and totals
    allocate(iSiteAtomStart(nSite+1))
    allocate(iSiteOrbStart(nSite+1))
    iSiteAtomStart(1) = 1
    iSiteOrbStart(1) = 1
    do iSite = 1, nSite
      iSiteAtomStart(iSite+1) = iSiteAtomStart(iSite) + ptrsPhase1(iSite)%nAtom
      iSiteOrbStart(iSite+1) = iSiteOrbStart(iSite) + ptrsPhase1(iSite)%nOrb
    end do

    ! total number of fragment orbitals to consider = dimension of resulting Hamiltonian
    nFO = 0
    do iSite = 1, nSite
      nFO = nFO + ptrsPhase1(iSite)%nFO
    end do
    ! indFO(i) is the index into Tij & Sij at which the orbitals on fragment (i) start
    allocate(indFO(iSite))
    indFO(1) = 1
    do iSite = 1, nSite-1
      indFO(iSite+1) = indFO(iSite) + ptrsPhase1(iSite)%nFO
    end do

    call env%globalTimer%startTimer(globalTimers%preSccInit)

    call handleCoordinateChange(env, this%coord0, this%latVec, this%invLatVec, this%species0,&
        & this%cutOff, this%orb, this%tPeriodic, this%tHelical, this%sccCalc, this%dispersion,&
        & this%solvation, this%thirdOrd, this%rangeSep, this%reks, this%img2CentCell,&
        & this%iCellVec, this%neighbourList, this%nAllAtom, this%coord0Fold, this%coord,&
        & this%species, this%rCellVec, this%nNeighbourSk, this%nNeighbourRep, this%nNeighbourLC,&
        & this%ham, this%over, this%H0, this%rhoPrim, this%iRhoPrim, this%iHam, this%ERhoPrim,&
        & this%iSparseStart, this%tPoisson, this%cm5Cont)

    ! CHARGE-INDEPENDENT MATRICES -- SPARSE FORM
    call env%globalTimer%startTimer(globalTimers%sparseH0S)
    call buildH0(env, this%H0, this%skHamCont, this%atomEigVal, this%coord, this%nNeighbourSk,&
        & this%neighbourList%iNeighbour, this%species, this%iSparseStart, this%orb)
    call buildS(env, this%over, this%skOverCont, this%coord, this%nNeighbourSk,&
        & this%neighbourList%iNeighbour, this%species, this%iSparseStart, this%orb)
    call env%globalTimer%stopTimer(globalTimers%sparseH0S)

    ! CHARGE-INDEPENDENT MATRICES -- CONVERT TO DENSE
    allocate(HSqrReal(this%nOrb,this%nOrb))
    allocate(SSqrReal(this%nOrb,this%nOrb))
    call unpackHS(HSqrReal, this%H0, this%neighbourList%iNeighbour, this%nNeighbourSK,&
        & this%denseDesc%iAtomStart, this%iSparseStart, this%img2CentCell)
    call unpackHS(SSqrReal, this%over, this%neighbourList%iNeighbour, this%nNeighbourSK,&
        & this%denseDesc%iAtomStart, this%iSparseStart, this%img2CentCell)

    ! CHARGE-INDEPENDENT MATRICES -- INSERT DIAGONAL BLOCKS FROM PHASE 1
    do iSite = 1, nSite
      iOrbBeg = iSiteOrbStart(iSite)
      iOrbEnd = iSiteOrbStart(iSite+1) - 1
      HSqrReal(iOrbBeg:iOrbEnd,iOrbBeg:iOrbEnd) = ptrsPhase1(iSite)%denseH0(:,:)
      SSqrReal(iOrbBeg:iOrbEnd,iOrbBeg:iOrbEnd) = ptrsPhase1(iSite)%denseOver(:,:)
    end do

    ! CHARGE-INDEPENDENT MATRICES -- CONVERT TO SPARSE
    call packHS(this%h0, HSqrReal, this%neighbourList%iNeighbour, this%nNeighbourSK,&
        & this%orb%mOrb, this%denseDesc%iAtomStart, this%iSparseStart, this%img2CentCell)
    call packHS(this%over, SSqrReal, this%neighbourList%iNeighbour, this%nNeighbourSK,&
        & this%orb%mOrb, this%denseDesc%iAtomStart, this%iSparseStart, this%img2CentCell)

    ! EXTERNAL POTENTIALS -- UNIFORM ELECTRIC FIELD
    call resetExternalPotentials(this%refExtPot, this%potential)
    if (this%tEField) then
      nAtom = size(this%nNeighbourSK)
      this%EField(:) = this%EFieldStrength * this%EfieldVector
      do iAtom = 1, nAtom
        this%potential%extAtom(iAtom,1) = this%potential%extAtom(iAtom,1)&
            & + dot_product(this%coord(:, iAtom), this%EField)
      end do
      this%potential%extGrad(:,:) = this%potential%extGrad + spread(this%EField, 2, nAtom)
    else
      this%EField(:) = 0.0_dp
    end if
    call mergeExternalPotentials(this%orb, this%species, this%potential)

    call env%globalTimer%stopTimer(globalTimers%preSccInit)

    call env%globalTimer%startTimer(globalTimers%scc)

!   lpSCC: do iSccIter = 1, this%maxSccIter

    call resetInternalPotentials(.false., this%xi, this%orb, this%species, this%potential)

    ! GET ATOM CHARGES FROM THE INDIVIDUAL FRAGMENTS / SITES
    ! qInput(mOrb, nAtom, nSpin); nSpin == 1
    do iSite = 1, nSite
      iAtomBeg = iSiteAtomStart(iSite)
      iAtomEnd = iSiteAtomStart(iSite+1)-1
      this%qInput(:,iAtomBeg:iAtomEnd,1) = ptrsPhase1(iSite)%qInput(:,:,1)
    end do

    ! GET SHELL CHARGES FROM THE INDIVIDUAL FRAGMENTS / SITES
    call getChargePerShell(this%qInput, this%orb, this%species, this%chargePerShell)

    call addChargePotentials(env, this%sccCalc, this%qInput, this%q0, this%chargePerShell,&
        & this%orb, this%species, this%neighbourList, this%img2CentCell, this%spinW,&
        & this%solvation, this%thirdOrd, this%potential, this%electrostatics, this%tPoisson,&
        & this%tUpload, this%shiftPerLUp, this%dispersion)

    ! All potentials are added up into intBlock
    this%potential%intBlock = this%potential%intBlock + this%potential%extBlock

    ! DO WE NEED THIS?
    if (allocated(this%qDepExtPot)) then
      allocate(dQ(this%orb%mShell, this%nAtom, this%nSpin))
      call getChargePerShell(this%qInput, this%orb, this%species, dQ, qRef=this%q0)
      call this%qDepExtPot%addPotential(sum(dQ(:,:,1), dim=1), dQ(:,:,1), this%orb,&
          & this%species, this%potential%intBlock)
    end if

    ! SELF-CONSISTENT / CHARGE-DEPENDENT HAMILTONIAN
    call getSccHamiltonian(this%H0, this%over, this%nNeighbourSK, this%neighbourList,&
        & this%species, this%orb, this%iSparseStart, this%img2CentCell, this%potential,&
        & allocated(this%reks), this%ham, this%iHam)
    call unpackHS(HSqrReal, this%ham(:,1), this%neighbourList%iNeighbour, this%nNeighbourSK,&
        & this%denseDesc%iAtomStart, this%iSparseStart, this%img2CentCell)

!   end do lpSCC

    call env%globalTimer%stopTimer(globalTimers%scc)

    call env%globalTimer%startTimer(globalTimers%postSCC)

!   if (this%isLinResp) then
!     if (.not. this%isRS_LinResp) then
!       call calculateLinRespExcitations(env, this%linearResponse, this%parallelKS, this%sccCalc,&
!           & this%qOutput, this%q0, this%over, this%eigvecsReal, this%eigen(:,1,:),&
!           & this%filling(:,1,:), this%coord, this%species, this%speciesName, this%orb,&
!           & this%skHamCont, this%skOverCont, autotestTag, this%taggedWriter, this%runId,&
!           & this%neighbourList, this%nNeighbourSK, this%denseDesc, this%iSparseStart,&
!           & this%img2CentCell, this%tWriteAutotest, this%tCasidaForces, this%tLinRespZVect,&
!           & this%tPrintExcitedEigvecs, this%tPrintEigvecsTxt, this%nonSccDeriv,&
!           & this%dftbEnergy(1), this%energiesCasida, this%SSqrReal, this%rhoSqrReal,&
!           & this%excitedDerivs, this%dQAtomEx, this%occNatural)
!     else
!       call calculateLinRespExcitations_RS(env, this%linearResponse, this%parallelKS,&
!           & this%sccCalc, this%qOutput, this%q0, this%over, this%eigvecsReal, this%eigen(:,1,:),&
!           & this%filling(:,1,:), this%coord0, this%species, this%speciesName, this%orb,&
!           & this%skHamCont, this%skOverCont, autotestTag, this%taggedWriter, this%runId,&
!           & this%neighbourList, this%nNeighbourSK, this%denseDesc, this%iSparseStart,&
!           & this%img2CentCell, this%tWriteAutotest, this%tCasidaForces, this%tLinRespZVect,&
!           & this%tPrintExcitedEigvecs, this%tPrintEigvecsTxt, this%nonSccDeriv,&
!           & this%dftbEnergy(1), this%energiesCasida, this%SSqrReal, this%deltaRhoOutSqr,&
!           & this%excitedDerivs, this%dQAtomEx, this%occNatural, this%rangeSep)
!     end if
!   end if

!   if (allocated(this%ppRPA)) then
!     call unpackHS(this%SSqrReal, this%over, this%neighbourList%iNeighbour, this%nNeighbourSK,&
!         & this%denseDesc%iAtomStart, this%iSparseStart, this%img2CentCell)
!     call blockSymmetrizeHS(this%SSqrReal, this%denseDesc%iAtomStart)
!     if (withMpi) then
!       call error("pp-RPA calc. does not work with MPI yet")
!     end if
!     call ppRPAenergies(this%ppRPA, this%denseDesc, this%eigvecsReal, this%eigen(:,1,:),&
!         & this%sccCalc, this%SSqrReal, this%species0, this%nEl(1), this%neighbourList%iNeighbour,&
!         & this%img2CentCell, this%orb, this%tWriteAutotest, autotestTag, this%taggedWriter)
!   end if

!   if (this%tWriteDetailedOut  .and. this%deltaDftb%nDeterminant() == 1) then
!     call writeDetailedOut4(this%fdDetailedOut, this%tSccCalc, tConverged, this%isXlbomd,&
!         & this%isLinResp, this%isGeoOpt, this%tMD, this%tPrintForces, this%tStress,&
!         & this%tPeriodic, this%dftbEnergy(this%deltaDftb%iDeterminant), this%totalStress,&
!         & this%totalLatDeriv, this%derivs, this%chrgForces, this%indMovedAtom, this%cellVol,&
!         & this%intPressure, this%geoOutFile, this%iAtInCentralRegion)
!   end if

    ! CALCULATE THE FMO HAMILTONIAN AND OVERLAP
    allocate(Tij(nFO,nFO))
    allocate(Sij(nFO,nFO))
    Tij = 0._dp
    Sij = 0._dp

    do iSite = 1, nSite
      do jSite = iSite, nSite
        do iOrb = 1, ptrsPhase1(iSite)%nFO
          iFO = indFO(iSite) + iOrb - 1
          do jOrb = 1, ptrsPhase1(iSite)%nFO
            jFO = indFO(jSite) + jOrb - 1

            Tij(iFO,jFO) = 0._dp
            Sij(iFO,jFO) = 0._dp

            do iAO = 1, ptrsPhase1(iSite)%nOrb ! iSiteOrbStart(iSite), iSiteOrbStart(iSite+1)-1
              do jAO = 1, ptrsPhase1(jSite)%nOrb ! iSiteOrbStart(jSite), iSiteOrbStart(jSite+1)-1
                Tij(iFO,jFO) = Tij(iFO,jFO)&
                    & + ptrsPhase1(iSite)%eigVec(iAO,iSiteOrbStart(iSite)+ptrsPhase1(iSite)%iHOMO-iOrb)&
                    & * HSqrReal(iSiteOrbStart(iSite)+iAO-1,iSiteOrbStart(jSite)+jAO-1)&
                    & * ptrsPhase1(jSite)%eigVec(jAO,iSiteOrbStart(jSite)+ptrsPhase1(jSite)%iHOMO-jOrb)
                Sij(iFO,jFO) = Sij(iFO,jFO)&
                    & + ptrsPhase1(iSite)%eigVec(iAO,iSiteOrbStart(iSite)+ptrsPhase1(iSite)%iHOMO-iOrb)&
                    & * SSqrReal(iSiteOrbStart(iSite)+iAO-1,iSiteOrbStart(jSite)+jAO-1)&
                    & * ptrsPhase1(jSite)%eigVec(jAO,iSiteOrbStart(jSite)+ptrsPhase1(jSite)%iHOMO-jOrb)
              end do
            end do
          end do ! jOrb
        end do ! iOrb
      end do ! jSite
    end do ! iSite

    ! COPY TO THE OTHER TRIANGLE OF THE MATRIX
    do iFO = 1, nFO
      do jFO = iFO+1, nFO
        Tij(jFO,iFO) = Tij(iFO,jFO)
        Sij(jFO,iFO) = Sij(iFO,jFO)
      end do
    end do
        
    ! PUT THE EIGENVALUES FROM PHASE 1 TO THE DIAGONAL
    iFO = 0
    do iSite = 1, nSite
      do iOrb = 1, ptrsPhase1(iSite)%nFO
        iFO = iFO + 1
        Tij(iFO,iFO) = ptrsPhase1(iSite)%eigVal(ptrsPhase1(iSite)%iHOMO + 1 - iOrb)
      end do
    end do

    ! ORTHOGONALIZE THE FMO HAMILTONIAN
    allocate(TijOrtho(nFO,nFO))
    call orthogonalizeHamiltonian(nFO, Tij, Sij, TijOrtho)

  end subroutine processGeometryPhase2


  subroutine orthogonalizeHamiltonian(n, tij, sij, tijOrtho)

    !> size of matrices
    integer, intent(in) :: n

    !> non-orthogonal Hamiltonian
    real(dp), allocatable, intent(in) :: tij(:,:)

    !> overlap
    real(dp), allocatable, intent(inout) :: sij(:,:)

    !> resulting, orthogonalized Hamiltonian
    real(dp), allocatable, intent(out) :: tijOrtho(:,:)

    ! for the orthogonalization of FMO Hamiltonian
    character, parameter :: uplo = 'U', jobz = 'V', range = 'A'
    real(dp), parameter :: absTol = 1.e-8_dp
    integer :: i, j, info, nEval, dummyInt
    real(dp) :: dummyReal
    real(dp), allocatable :: eval(:), evec(:,:), work(:), sqrt_eval(:,:)
    integer, allocatable :: isuppz(:), iwork(:)

    ! Cholesky-factorize the overlap matrix
    call dpotrf(uplo, n, sij, n, info);
    ! invert the overlap matrix
    call dpotri(uplo, n, sij, n, info);
    ! fill the other triangle
    do i = 1, n
      do j = 1, i-1
        sij(i,j) = sij(j,i)
      end do
    end do

    ! diagonalize the inverse of overlap
    allocate(eval(n))
    allocate(evec(n,n))
    allocate(isuppz(2*n))
    allocate(work(26*n))
    allocate(iwork(10*n))
    call dsyevr(jobz, range, uplo, n, sij, n, dummyReal, dummyReal, dummyInt, dummyInt,&
        & absTol, nEval, eval, evec, n, isuppz, work, 26*n, iwork, 10*n, info)
    deallocate(isuppz)
    deallocate(work)
    deallocate(iwork)

    ! Sij = evec * sqrt(diag) * evec-transp
    allocate(sqrt_eval(1,n))
    sqrt_eval(1,:) = sqrt(eval(:))
    sij = matmul(evec, matmul(sqrt_eval,transpose(evec)))
    deallocate(sqrt_eval)
    deallocate(evec)
    deallocate(eval)

    ! THamilOr = sij * tij * sij
    tijOrtho = matmul(sij, matmul(tij,sij))

  end subroutine orthogonalizeHamiltonian

end module dftbp_fmo
