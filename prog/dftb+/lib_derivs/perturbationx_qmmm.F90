!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module for linear response derivative calculations -- QM/MM derivatives
module dftbp_perturbxderivs_qmmm
  use dftbp_accuracy
  use dftbp_blasroutines
  use dftbp_commontypes
  use dftbp_constants
  use dftbp_densedescr
  use dftbp_dftbplusu
  use dftbp_environment
  use dftbp_finitethelper
  use dftbp_globalenv
  use dftbp_mainio
  use dftbp_message
  use dftbp_mixer
#:if WITH_MPI
  use dftbp_mpifx
#:endif
  use dftbp_nonscc, only : TNonSccDiff
  use dftbp_onsitecorrection
  use dftbp_orbitalequiv
  use dftbp_periodic
  use dftbp_populations
  use dftbp_potentials
  use dftbp_rangeseparated, only : TRangeSepFunc
  use dftbp_rotateDegenerateOrbs
  use dftbp_scalapackfx
  use dftbp_scc
  use dftbp_shift
  use dftbp_slakocont
  use dftbp_sparse2dense
  use dftbp_spin
  use dftbp_taggedoutput
  use dftbp_thirdorder, only : TThirdOrder

  implicit none

  private
  public :: dPsidxQMMM

  !> Direction labels
  character(len=1), parameter :: direction(3) = ['x','y','z']

contains

  !> Static (frequency independent) perturbation at q=0
  subroutine dPsidxQMMM(env, parallelKS, filling, eigvals, eigVecsReal, qOrb, q0, ham, over, orb,&
      & nAtom, species, neighbourList, nNeighbourSK, denseDesc, iSparseStart, img2CentCell, coord,&
      & sccCalc, maxSccIter, sccTol, nMixElements, nIneqMixElements, iEqOrbitals, tempElec, Ef,&
      & tFixEf, spinW, thirdOrd, dftbU, iEqBlockDftbu, onsMEs, iEqBlockOnSite, rangeSep,&
      & nNeighbourLC, pChrgMixer, taggedWriter, tWriteAutoTest, autoTestTagFile, tWriteTaggedOut,&
      & taggedResultsFile, tWriteDetailedOut, fdDetailedOut, tMulliken)
      
    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Fillings of unperturbed system
    real(dp), intent(in) :: filling(:,:,:)

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(in) :: eigvals(:,:,:)

    !> ground state eigenvectors
    real(dp), intent(in), allocatable :: eigVecsReal(:,:,:)

    !> Electrons in each atomic orbital
    real(dp), intent(in) :: qOrb(:,:,:)

    !> reference charges
    real(dp), intent(in) :: q0(:,:,:)

    !> Sparse Hamiltonian
    real(dp), intent(in) :: ham(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Number of central cell atoms
    integer, intent(in) :: nAtom

    !> chemical species
    integer, intent(in) :: species(:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> SCC module internal variables
    type(TScc), intent(inout), allocatable :: sccCalc

    !> maximal number of SCC iterations
    integer, intent(in) :: maxSccIter

    !> Tolerance for SCC convergence
    real(dp), intent(in) :: sccTol

    !> nr. of elements to go through the mixer - may contain reduced orbitals and also orbital
    !> blocks (if tDFTBU or onsite corrections)
    integer, intent(in) :: nMixElements

    !> nr. of inequivalent charges
    integer, intent(in) :: nIneqMixElements

    !> Equivalence relations between orbitals
    integer, intent(in) :: iEqOrbitals(:,:,:)

    !> onsite matrix elements for shells (elements between s orbitals on the same shell are ignored)
    real(dp), intent(in), allocatable :: onsMEs(:,:,:,:)

    !> Equivalences for onsite block corrections if needed
    integer, intent(in), allocatable :: iEqBlockOnSite(:,:,:,:)

    !> Electron temperature
    real(dp), intent(in) :: tempElec

    !> Fermi level(s)
    real(dp), intent(in) :: Ef(:)

    !> Whether fixed Fermi level(s) should be used. (No charge conservation!)
    logical, intent(in) :: tFixEf

    !> spin constants
    real(dp), intent(in), allocatable :: spinW(:,:,:)

    !> Third order SCC interactions
    type(TThirdOrder), allocatable, intent(inout) :: thirdOrd

    !> Are there orbital potentials present
    type(TDftbU), intent(in), allocatable :: dftbU

    !> equivalence mapping for dual charge blocks
    integer, intent(in), allocatable :: iEqBlockDftbu(:,:,:,:)

    !> Data for range-separated calculation
    type(TRangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> Number of neighbours for each of the atoms for the exchange contributions in the long range
    !> functional
    integer, intent(inout), allocatable :: nNeighbourLC(:)

    !> Charge mixing object
    type(TMixer), intent(inout) :: pChrgMixer

    !> Tagged writer object
    type(TTaggedWriter), intent(inout) :: taggedWriter

    !> should regression test data be written
    logical, intent(in) :: tWriteAutoTest

    !> File name for regression data
    character(*), intent(in) :: autoTestTagFile

    !> should machine readable output data be written
    logical, intent(in) :: tWriteTaggedOut

    !> File name for machine readable results data
    character(*), intent(in) :: taggedResultsFile

    !> should detailed.out be written to
    logical, intent(in) :: tWriteDetailedOut

    !> File id for detailed.out
    integer, intent(in) :: fdDetailedOut

    !> Should Mulliken populations be generated/output
    logical, intent(in) :: tMulliken

    integer :: iS, iK, iKS, iAt, iCart, iLev, iSh, iSp, jAt, jCart

    integer :: nSpin, nKpts, nOrbs, nIndepHam

    ! maximum allowed number of electrons in a single particle state
    real(dp) :: maxFill

    integer, allocatable :: nFilled(:,:), nEmpty(:,:)

    integer :: iSCCIter
    logical :: tStopSCC

    ! matrices for derivatives of terms in hamiltonian and outputs
    real(dp), allocatable :: dHam(:,:)

    ! overlap derivative terms in potential, omega dS + d(delta q gammma) S
    ! QM/MM: omega dS == 0, so we only have d(delta q gammma) S, formerly sOmega(:,:,2)
    real(dp), allocatable :: sOmega(:,:)

    real(dp) :: dRho(size(over),size(ham, dim=2))
    real(dp) :: dqIn(orb%mOrb,nAtom,size(ham, dim=2))
    real(dp) :: dqOut(orb%mOrb, nAtom, size(ham, dim=2), 3, nAtom)
    real(dp) :: dqInpRed(nMixElements), dqOutRed(nMixElements)
    real(dp) :: dqDiffRed(nMixElements), sccErrorQ
    real(dp) :: dqPerShell(orb%mShell,nAtom,size(ham, dim=2))

    real(dp), allocatable :: vAt(:,:), vdgamma(:,:,:)

    real(dp), allocatable :: dqBlockIn(:,:,:,:), SSqrReal(:,:)
    real(dp), allocatable :: dqBlockOut(:,:,:,:)
    real(dp), allocatable :: dummy(:,:,:,:)

    ! derivative of potentials
    type(TPotentials) :: dPotential

    real(dp), allocatable :: shellPot(:,:,:), atomPot(:,:)

    logical :: tSccCalc, tConverged
    logical, allocatable :: tMetallic(:)

    real(dp), allocatable :: dEi(:,:,:,:)
    real(dp), allocatable :: dPsiReal(:,:,:,:)

    integer :: fdResults

    ! used for range separated contributions, note this stays in the up/down representation
    ! throughout if spin polarised
    real(dp), pointer :: dRhoOutSqr(:,:,:), dRhoInSqr(:,:,:)
    real(dp), allocatable, target :: dRhoOut(:), dRhoIn(:)

    real(dp) :: dDipole(3)

    ! QM/MM specific
    ! coordinates and charges of MM atoms
    integer :: nExtCharge
    real(dp), allocatable :: extCoord(:,:), extCharge(:), dgammaQMMM(:,:,:)

    ! obtain the coordinates and charges of MM atoms from SCC structures
    call sccCalc%getExternalCharges(nExtCharge, extCoord, extCharge)
    if (nExtCharge <= 0) then
      write (*,*) "No MM atoms, nothing to do in dPsidxQMMM."
      return
    end if

    if (tFixEf) then
      call error("Perturbation expressions not currently implemented for fixed Fermi energy")
    end if

    nSpin = size(ham, dim=2)
    select case(nSpin)
    case(1,4)
      nIndepHam = 1
    case(2)
      nIndepHam = 2
    end select
    select case(nSpin)
    case(1)
      maxFill = 2.0_dp
    case(2,4)
      maxFill = 1.0_dp
    end select

    allocate(tMetallic(nIndepHam))

    nOrbs = size(filling,dim=1)
    nKpts = size(filling,dim=2)

    allocate(dEi(nOrbs, nAtom, nSpin, 3))

    allocate(dHam(size(ham,dim=1),nSpin))

    tSccCalc = allocated(sccCalc)

    ! terms v S' and v' S
    if (tSccCalc) then
      allocate(sOmega(size(ham,dim=1),nSpin))
      allocate(vAt(nAtom,nSpin))
      allocate(vdgamma(orb%mShell,nAtom,nSpin))
    end if

    if (allocated(rangeSep)) then
      allocate(SSqrReal(nOrbs, nOrbs))
      SSqrReal(:,:) = 0.0_dp
      call unpackHS(SSqrReal, over, neighbourList%iNeighbour, nNeighbourSK, denseDesc%iAtomStart,&
          & iSparseStart, img2CentCell)
      allocate(dRhoOut(nOrbs * nOrbs * nSpin))
      dRhoOutSqr(1:nOrbs, 1:nOrbs, 1:nSpin) => dRhoOut(:nOrbs*nOrbs*nSpin)
      allocate(dRhoIn(nOrbs * nOrbs * nSpin))
      dRhoInSqr(1:nOrbs, 1:nOrbs, 1:nSpin) => dRhoIn(:nOrbs*nOrbs*nSpin)
    else
      dRhoInSqr => null()
      dRhoOutSqr => null()
    end if

    if (allocated(dftbU) .or. allocated(onsMEs)) then
      allocate(dqBlockIn(orb%mOrb,orb%mOrb,nAtom,nSpin))
      allocate(dqBlockOut(orb%mOrb,orb%mOrb,nAtom,nSpin))
    end if

    call init(dPotential, orb, nAtom, nSpin)

    allocate(nFilled(nIndepHam, nKpts))
    allocate(nEmpty(nIndepHam, nKpts))

    if (allocated(spinW) .or. allocated(thirdOrd)) then
      allocate(shellPot(orb%mShell, nAtom, nSpin))
    end if
    if (allocated(thirdOrd)) then
      allocate(atomPot(nAtom, nSpin))
    end if

    nFilled(:,:) = -1
    do iS = 1, nIndepHam
      do iK = 1, nKPts
        do iLev = 1, nOrbs
          if ( filling(iLev,iK,iS) < epsilon(1.0) ) then
            nFilled(iS,iK) = iLev - 1
            exit
          end if
        end do
        if (nFilled(iS, iK) < 0) then
          nFilled(iS, iK) = nOrbs
        end if
      end do
    end do
    nEmpty(:,:) = -1
    do iS = 1, nIndepHam
      do iK = 1, nKpts
        do iLev = 1, nOrbs
          if ( abs( filling(iLev,iK,iS) - maxFill ) > epsilon(1.0)) then
            nEmpty(iS, iK) = iLev
            exit
          end if
        end do
        if (nEmpty(iS, iK) < 0) then
          nEmpty(iS, iK) = 1
        end if
      end do
    end do

    do iS = 1, nIndepHam
      tMetallic(iS) = (.not.all(nFilled(iS,:) == nEmpty(iS,:) -1))
      !write(stdOut,*)'Fractionally filled range'
      !do iK = 1, nKpts
      !  write(stdOut,*) nEmpty(:,iK), ':', nFilled(:,iK)
      !end do
    end do

    dqOut(:,:,:,:,:) = 0.0_dp
    dEi(:,:,:,:) = 0.0_dp

    ! derivatives of QM--MM 1/r w.r.t. coordinates of MM atoms
    allocate(dgammaQMMM(3, nAtom, nExtCharge))
    call calcInvRPrimeQMMM(nAtom, nExtCharge, coord, extCoord, extCharge, dgammaQMMM)

    ! Displaced MM atom to differentiate wrt
    lpAtom: do iAt = 1, nExtCharge

      ! any non-variational QM/MM contribution?

      ! perturbation direction
      lpCart: do iCart = 1, 3

        write (stdOut,*) 'Calculating derivative for displacement along ',&
            & trim(direction(iCart)),' for MM charge number', iAt

        if (tSccCalc) then
          sOmega(:,:) = 0.0_dp
        end if

        if (tSccCalc) then

          vdgamma(:,:,:) = 0.0_dp
          vAt(:,:) = 0.0_dp

          call sccCalc%updateCoords(env, coord, species, neighbourList)
          call sccCalc%updateCharges(env, qOrb, q0, orb, species)
          vAt(:,1) = dgammaQMMM(iCart,:,iAt)
          call total_shift(vdgamma, vAt, orb, species)

        end if

        dqIn(:,:,:) = 0.0_dp
        if (allocated(dftbU) .or. allocated(onsMEs)) then
          dqBlockIn(:,:,:,:) = 0.0_dp
          dqBlockOut(:,:,:,:) = 0.0_dp
        end if

        if (tSccCalc) then
          call reset(pChrgMixer, nMixElements)
          dqInpRed(:) = 0.0_dp
          dqPerShell(:,:,:) = 0.0_dp
          if (allocated(rangeSep)) then
            dRhoIn(:) = 0.0_dp
            dRhoOut(:) = 0.0_dp
          end if
        end if

        if (tSccCalc) then
          write (stdOut, "(1X,A,T12,A)") 'SCC Iter' , 'Error'
        end if

        iSCCIter = 1
        tStopSCC = .false.
        lpSCC: do while (iSCCiter <= maxSccIter)

          dPotential%intAtom(:,:) = 0.0_dp
          dPotential%intShell(:,:,:) = 0.0_dp
          dPotential%intBlock(:,:,:,:) = 0.0_dp

          if (allocated(dftbU) .or. allocated(onsMEs)) then
            dPotential%orbitalBlock(:,:,:,:) = 0.0_dp
          end if

          if (tSccCalc .and. iSCCiter>1) then
            call sccCalc%updateCharges(env, dqIn, orb=orb, species=species)
            call sccCalc%updateShifts(env, orb, species, neighbourList%iNeighbour, img2CentCell)
            call sccCalc%getShiftPerAtom(dPotential%intAtom(:,1))
            call sccCalc%getShiftPerL(dPotential%intShell(:,:,1))

            if (allocated(spinW)) then
              call getChargePerShell(dqIn, orb, species, dqPerShell)
              call getSpinShift(shellPot, dqPerShell, species, orb, spinW)
              dPotential%intShell(:,:,:) = dPotential%intShell + shellPot
            end if

            if (allocated(thirdOrd)) then
              atomPot(:,:) = 0.0_dp
              shellPot(:,:,:) = 0.0_dp
              call thirdOrd%getdShiftdQ(atomPot(:,1), shellPot(:,:,1), species, neighbourList,&
                  & dqIn, img2CentCell, orb)
              dPotential%intAtom(:,1) = dPotential%intAtom(:,1) + atomPot(:,1)
              dPotential%intShell(:,:,1) = dPotential%intShell(:,:,1) + shellPot(:,:,1)
            end if

            if (allocated(dftbU)) then
              ! note the derivatives of both FLL and pSIC are the same (pSIC, i.e. case 2 in module)
              call dftbU%getDftbUShift(dPotential%orbitalBlock, dqBlockIn, species, orb)
            end if
            if (allocated(onsMEs)) then
              ! onsite corrections
              call addOnsShift(dPotential%orbitalBlock, dPotential%iOrbitalBlock,&
                  & dqBlockIn, dummy, onsMEs=onsMEs, species=species, orb=orb)
            end if

          end if

          call total_shift(dPotential%intShell, dPotential%intAtom, orb, species)
          call total_shift(dPotential%intBlock, dPotential%intShell, orb, species)
          if (allocated(dftbU) .or. allocated(onsMEs)) then
            dPotential%intBlock(:,:,:,:) = dPotential%intBlock + dPotential%orbitalBlock
          end if

          if (tSccCalc) then
            sOmega(:,:) = 0.0_dp
            ! add the (Delta q) * d gamma / dx term
            call add_shift(sOmega(:,:), over, nNeighbourSK, neighbourList%iNeighbour, species,&
                & orb, iSparseStart, nAtom, img2CentCell, vdgamma)
            ! and add gamma * d (Delta q) / dx
            call add_shift(sOmega(:,:), over, nNeighbourSK, neighbourList%iNeighbour, species,&
                & orb, iSparseStart, nAtom, img2CentCell, dPotential%intBlock)
          end if

          dHam = 0.0_dp

          if (tSccCalc) then
            dHam(:,:) = dHam + sOmega(:,:)
          end if

          if (nSpin > 1) then
            dHam(:,:) = 2.0_dp * dHam(:,:)
          end if

          call qm2ud(dHam)

          ! evaluate derivative of density matrix
          dRho(:,:) = 0.0_dp
          do iKS = 1, parallelKS%nLocalKS

            iS = parallelKS%localKS(2, iKS)

            if (allocated(dRhoOut)) then
              ! replace with case that will get updated in dRhoReal
              dRhoOutSqr(:,:,iS) = dRhoInSqr(:,:,iS)
            end if

            call dRhoRealQMMM(env, dHam, neighbourList, nNeighbourSK, iSparseStart,&
                & img2CentCell, denseDesc, iKS, parallelKS, nFilled(:,1), nEmpty(:,1),&
                & eigVecsReal, eigVals, Ef, tempElec, orb, dRho(:,iS), iCart, dRhoOutSqr,&
                & rangeSep, over, nNeighbourLC, tMetallic, filling / maxFill, dEi,&
                & dPsiReal, iAt)
          end do

          dRho(:,:) = maxFill * dRho
          if (allocated(dRhoOut)) then
            dRhoOut(:) = maxFill * dRhoOut
          end if
          call ud2qm(dRho)

          dqOut(:, :, :, iCart, iAt) = 0.0_dp
          do iS = 1, nSpin
            call mulliken(dqOut(:, :, iS, iCart, iAt), over, dRho(:,iS), orb,&
                & neighbourList%iNeighbour, nNeighbourSK, img2CentCell, iSparseStart)
            if (allocated(dftbU) .or. allocated(onsMEs)) then
              dqBlockOut(:,:,:,iS) = 0.0_dp
              call mulliken(dqBlockOut(:,:,:,iS), over, dRho(:,iS), orb, neighbourList%iNeighbour,&
                  & nNeighbourSK, img2CentCell, iSparseStart)
            end if
          end do

          if (tSccCalc) then

            if (allocated(rangeSep)) then
              dqDiffRed(:) = dRhoOut - dRhoIn
            else
              dqOutRed = 0.0_dp
              call OrbitalEquiv_reduce(dqOut(:, :, :, iCart, iAt), iEqOrbitals, orb,&
                  & dqOutRed(:nIneqMixElements))
              if (allocated(dftbU)) then
                call appendBlockReduced(dqBlockOut, iEqBlockDFTBU, orb, dqOutRed)
              end if
              if (allocated(onsMEs)) then
                call onsBlock_reduce(dqBlockOut, iEqBlockOnSite, orb, dqOutRed)
              end if
              dqDiffRed(:) = dqOutRed - dqInpRed
            end if
            sccErrorQ = maxval(abs(dqDiffRed))

            write(stdOut,"(1X,I0,T10,E20.12)")iSCCIter, sccErrorQ
            tConverged = (sccErrorQ < sccTol)

            if ((.not. tConverged) .and. iSCCiter /= maxSccIter) then
              if (iSCCIter == 1) then
                if (allocated(rangeSep)) then
                  dRhoIn(:) = dRhoOut
                  call denseMulliken(dRhoInSqr, SSqrReal, denseDesc%iAtomStart, dqIn)
                else
                  dqIn(:,:,:) = dqOut(:, :, :, iCart, iAt)
                  dqInpRed(:) = dqOutRed(:)
                  if (allocated(dftbU) .or. allocated(onsMEs)) then
                    dqBlockIn(:,:,:,:) = dqBlockOut(:,:,:,:)
                  end if
                end if

              else

                if (allocated(rangeSep)) then

                  call mix(pChrgMixer, dRhoIn, dqDiffRed)
                  call denseMulliken(dRhoInSqr, SSqrReal, denseDesc%iAtomStart, dqIn)

                else

                  call mix(pChrgMixer, dqInpRed, dqDiffRed)
                #:if WITH_MPI
                  ! Synchronise charges in order to avoid mixers that store a history drifting apart
                  call mpifx_allreduceip(env%mpi%globalComm, dqInpRed, MPI_SUM)
                  dqInpRed(:) = dqInpRed / env%mpi%globalComm%size
                #:endif

                  call OrbitalEquiv_expand(dqInpRed(:nIneqMixElements), iEqOrbitals, orb, dqIn)

                  if (allocated(dftbU) .or. allocated(onsMEs)) then
                    dqBlockIn(:,:,:,:) = 0.0_dp
                    if (allocated(dftbU)) then
                      call dftbU%expandBlock(dqInpRed, iEqBlockDFTBU, orb, dqBlockIn, species(:nAtom),&
                          & orbEquiv=iEqOrbitals)
                    else
                      call Onsblock_expand(dqInpRed, iEqBlockOnSite, orb, dqBlockIn,&
                          & orbEquiv=iEqOrbitals)
                    end if
                  end if

                end if

              end if

              if (allocated(rangeSep)) then
                call ud2qm(dqIn)
              end if

            end if

          else

            tConverged = .true.

          end if

          if (tConverged) then
            exit lpSCC
          end if

          if (allocated(spinW)) then
            dqPerShell = 0.0_dp
            do jAt = 1, nAtom
              iSp = species(jAt)
              do iSh = 1, orb%nShell(iSp)
                dqPerShell(iSh,jAt,:nSpin) = dqPerShell(iSh,jAt,:nSpin) +&
                    & sum(dqIn(orb%posShell(iSh,iSp): orb%posShell(iSh+1,iSp)-1,jAt,:nSpin),dim=1)
              end do
            end do

          end if

          iSCCIter = iSCCIter +1

        end do lpSCC

      end do lpCart

    end do lpAtom

    write (stdOut, *) 'dEi'
    do iCart = 1, 3
      write (stdOut, *) iCart
      do iS = 1, nSpin
        do iAt = 1, nAtom
          write (stdOut, *) dEi(:, iAt, iS, iCart) ! * Hartree__eV
        end do
      end do
    end do

    if (tMulliken .or. tSccCalc) then
      write (stdOut, *)
      write (stdOut, *) 'Charge derivatives'
      do iAt = 1, nAtom
        write (stdOut,"(A,I0)") '/d Atom_', iAt
        do iS = 1, nSpin
          do jAt = 1, nAtom
            write (stdOut, *) jAt, -sum(dqOut(:, jAt, iS, :, iAt), dim=1)
          end do
          write (stdOut, *)
        end do
      end do
      write (stdOut, *)

      write (stdOut, *) 'Born effective charges'
      ! i.e. derivative of dipole moment wrt to atom positions, or equivalently derivative of forces
      ! wrt to a homogeneous electric field
      do iAt = 1, nAtom
        do iCart = 1, 3
          do jCart = 1, 3
            dDipole(jCart) = -sum(sum(dqOut(:, : , 1, iCart, iAt), dim=1) * coord(jCart, :))
          end do
          dDipole(iCart) = dDipole(iCart) -sum(qOrb(:,iAt,1) - q0(:,iAt,1))
          write (stdOut,*) dDipole
        end do
        write (stdOut, *)
      end do
      write (stdOut, *)

    end if

    if (tWriteAutoTest) then
      open(newunit=fdResults, file=autoTestTagFile, position="append")
      call taggedWriter%write(fdResults, tagLabels%dqdx, dqOut)
      close(fdResults)
    end if
    if (tWriteTaggedOut) then
      open(newunit=fdResults, file=taggedResultsFile, position="append")
      call taggedWriter%write(fdResults, tagLabels%dqdx, dqOut)
      close(fdResults)
    end if

  end subroutine dPsidxQMMM


  !> Calculate the derivative of density matrix from derivative of hamiltonian in static case at
  !> q=0, k=0
  subroutine dRhoRealQMMM(env, dHam, neighbourList, nNeighbourSK, iSparseStart, img2CentCell,&
      & denseDesc, iKS, parallelKS, nFilled, nEmpty, eigVecsReal, eigVals, Ef, tempElec, orb,&
      & dRhoSparse, iCart, dRhoSqr, rangeSep, over, nNeighbourLC, tMetallic, filling,&
      & dEi, dPsi, iAtom)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Derivative of the hamiltonian
    real(dp), intent(in) :: dHam(:,:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Particular spin/k-point
    integer, intent(in) :: iKS

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> ground state eigenvectors
    real(dp), intent(in) :: eigVecsReal(:,:,:)

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(in) :: eigvals(:,:,:)

    !> Fermi level(s)
    real(dp), intent(in) :: Ef(:)

    !> Last (partly) filled level in each spin channel
    integer, intent(in) :: nFilled(:)

    !> First (partly) empty level in each spin channel
    integer, intent(in) :: nEmpty(:)

    !> Electron temperature
    real(dp), intent(in) :: tempElec

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> returning dRhoSparse on exit
    real(dp), intent(out) :: dRhoSparse(:)

    !> Cartesian direction of perturbation
    integer, intent(in) :: iCart

    !> Derivative of rho as a square matrix, if needed
    real(dp), pointer :: dRhoSqr(:,:,:)

    !> Data for range-separated calculation
    type(TRangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> Number of neighbours for each of the atoms for the exchange contributions in the long range
    !> functional
    integer, intent(inout), allocatable :: nNeighbourLC(:)

    !> Is this a metallic system
    logical, intent(in) :: tMetallic(:)

    !> Fillings of unperturbed system
    real(dp), intent(in) :: filling(:,:,:)

    !> Derivative of single particle eigenvalues
    real(dp), intent(out) :: dEi(:,:,:,:)

    !> Optional derivatives of single particle wavefunctions
    real(dp), allocatable, intent(inout) :: dPsi(:,:,:,:)

    !> Atom with which the the derivative is being calculated
    integer, intent(in) :: iAtom

    integer :: ii, iFilled, iEmpty, iS, iK, nOrb
    real(dp) :: workLocal(size(eigVecsReal,dim=1), size(eigVecsReal,dim=2))
    real(dp) :: work3Local(size(eigVecsReal,dim=1), size(eigVecsReal,dim=2))
    real(dp), allocatable :: dRho(:,:)
    type(TDegeneracyTransform) :: transform

    real(dp), allocatable :: dFilling(:)

    iK = parallelKS%localKS(1, iKS)
    iS = parallelKS%localKS(2, iKS)

    if (tMetallic(iS)) then
      allocate(dFilling(size(dEi, dim=1)))
    end if

    call transform%init()

    dEi(:, iAtom, iS, iCart) = 0.0_dp
    if (allocated(dPsi)) then
      dPsi(:, :, iS, iCart) = 0.0_dp
    end if

    workLocal(:,:) = 0.0_dp
    allocate(dRho(size(eigVecsReal,dim=1), size(eigVecsReal,dim=2)))
    dRho(:,:) = 0.0_dp

    ! serial case
    nOrb = size(dRho, dim = 1)

    ! dH matrix in square form

    dRho(:,:) = 0.0_dp
    call unpackHS(dRho, dHam(:,iS), neighbourList%iNeighbour, nNeighbourSK, denseDesc%iAtomStart,&
        & iSparseStart, img2CentCell)

    if (allocated(rangeSep)) then
      call unpackHS(workLocal, over, neighbourList%iNeighbour, nNeighbourSK, denseDesc%iAtomStart,&
          & iSparseStart, img2CentCell)
      call rangeSep%addLRHamiltonian(env, dRhoSqr(:,:,iS), over, neighbourList%iNeighbour,&
          & nNeighbourLC, denseDesc%iAtomStart, iSparseStart, orb, dRho, workLocal)
    end if

    ! form H' |c>
    call symm(workLocal, 'l', dRho, eigVecsReal(:,:,iS))

    ! form |c> H' - e S' <c| -- actually: |c> H' <c|
    workLocal(:,:) = matmul(transpose(eigVecsReal(:,:,iS)), workLocal)

    call transform%generateUnitary(workLocal, eigvals(:,iK,iS))
    call transform%degenerateTransform(workLocal)

    ! diagonal elements of workLocal are now derivatives of eigenvalues
    do ii = 1, nOrb
      dEi(ii, iAtom, iS, iCart) = workLocal(ii, ii)
    end do

    if (tMetallic(iS)) then
      call dEida(dFilling, filling(:,iK,iS), dEi(:,iAtom, iS, iCart), tempElec)
      !write(stdOut,*)'dEf', dEfda(filling(:,iK,iS), dEi(:,iAtom, iS, iCart))
      !write(stdOut,*)dFilling
    end if

    work3Local = eigVecsReal(:,:,iS)
    call transform%applyUnitary(work3Local)

    ! Form actual perturbation U matrix for eigenvectors by weighting the elements
    do iFilled = 1, nFilled(iS)
      do iEmpty = 1, nOrb
        if (iFilled == iEmpty) then ! diagonal elements vanish
          workLocal(iFilled, iFilled) = 0.0_dp
        else
          if (.not. transform%degenerate(iFilled, iEmpty)) then
            workLocal(iEmpty, iFilled) = workLocal(iEmpty, iFilled)&
                & / (eigvals(iFilled, iK, iS) - eigvals(iEmpty, iK, iS))
          else
            workLocal(iEmpty, iFilled) = 0.0_dp
            workLocal(iFilled, iEmpty) = 0.0_dp
          end if
        end if
      end do
    end do

    ! calculate the derivatives of the eigenvectors
    workLocal(:, :nFilled(iS)) = matmul(work3Local, workLocal(:, :nFilled(iS)))

    if (allocated(dPsi)) then
      dPsi(:, :, iS, iCart) = work3Local
    end if

    do iFilled = 1, nOrb
      workLocal(:, iFilled) = workLocal(:, iFilled) * filling(iFilled, iK, iS)
    end do

    ! form the derivative of the density matrix
    dRho(:,:) = matmul(workLocal, transpose(work3Local)) + matmul(work3Local, transpose(workLocal))

    if (tMetallic(iS)) then
      ! extra contribution from change in Fermi level leading to change in occupations
      do iFilled = nEmpty(iS), nFilled(iS)
        workLocal(:, iFilled) = work3Local(:, iFilled) * dFilling(iFilled)
      end do
      dRho(:,:) = dRho + 0.5_dp * matmul(workLocal(:, nEmpty(iS):nFilled(iS)),&
          & transpose(work3Local(:, nEmpty(iS):nFilled(iS))))&
          & + 0.5 * matmul(work3Local(:, nEmpty(iS):nFilled(iS)),&
          & transpose(workLocal(:, nEmpty(iS):nFilled(iS))))
    end if

    dRhoSparse(:) = 0.0_dp
    call packHS(dRhoSparse, dRho, neighbourList%iNeighbour, nNeighbourSK, orb%mOrb,&
        & denseDesc%iAtomStart, iSparseStart, img2CentCell)

    if (associated(dRhoSqr)) then
      dRhoSqr(:,:,iS) = dRho
    end if

    call transform%destroy()

  end subroutine dRhoRealQMMM


  !> Calculates the -1/R**2 derivative of QM--MM distances w.r.t. coords of MM atoms 
  !>   directly multiplied by the MM charges (as that is included in the derivative of the shift)
  !> adopted from subroutine addInvRPrimeClusterAsymm in coulomb.F90
  subroutine calcInvRPrimeQMMM(nAtom0, nAtom1, coord0, coord1, charge1, deriv)

    !> Number of atoms in the first group -- QM atoms
    integer, intent(in) :: nAtom0

    !> Number of atoms in the second group -- MM atoms
    integer, intent(in) :: nAtom1

    !> List of atomic coordinates.
    real(dp), intent(in) :: coord0(:,:)

    !> List of the point charge coordinates
    real(dp), intent(in) :: coord1(:,:)

    !> Charge of the point charges.
    real(dp), intent(in) :: charge1(:)

    !> Derivative, deriv(3, #QM, #MM)
    real(dp), intent(out) :: deriv(:,:,:)

    integer :: iAt0, iAt1
    real(dp) :: dist, vect(3)

    @:ASSERT(size(deriv, dim=1) == 3)
    @:ASSERT(size(deriv, dim=2) == nAtom0)
    @:ASSERT(size(deriv, dim=3) == nAtom1)

    do iAt0 = 1, nAtom0
      do iAt1 = 1, nAtom1
        vect(:) = coord0(:,iAt0) - coord1(:,iAt1)
        dist = sqrt(sum(vect(:)**2))
        deriv(:,iAt0,iAt1) = charge1(iAt1) / (dist**3) * vect(:) ! sign looks OK
      end do
    end do

  end subroutine calcInvRPrimeQMMM


end module dftbp_perturbxderivs_qmmm
