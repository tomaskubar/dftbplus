!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module for linear response derivative calculations using perturbation methods
module dftbp_perturbxderivs_tk
  use dftbp_accuracy
  use dftbp_commontypes
  use dftbp_constants
  use dftbp_densedescr
  use dftbp_environment
  use dftbp_finitethelper
  use dftbp_globalenv
  use dftbp_mainio
  use dftbp_mixer
  use dftbp_nonscc, only : TNonSccDiff
  use dftbp_periodic
  use dftbp_potentials
  use dftbp_scc
  use dftbp_shift
  use dftbp_slakocont
  use dftbp_sparse2dense
  use dftbp_thirdorder, only : TThirdOrder

  implicit none

  private
  public :: dPsidx_TK

  !> Direction labels
  character(len=1), parameter :: direction(3) = ['x','y','z']

contains

  !> Static (frequency independent) perturbation at q=0
  subroutine dPsidx_TK(env, filling, eigVals, eigVecsReal, potential, qOrb, q0, overSparse,&
      & skHamCont, skOverCont, nonSccDeriv, orb, nAtom, species, neighbourList, nNeighbourSK,&
      & denseDesc, iSparseStart, img2CentCell, coord, sccCalc, maxSccIter, sccTol, tempElec,&
      & thirdOrd, pChrgMixer)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Fillings of unperturbed system
    real(dp), intent(in) :: filling(:)

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(in) :: eigVals(:)

    !> ground state eigenvectors
    real(dp), intent(in) :: eigVecsReal(:,:)

    !> Unperturbed potentials
    type(TPotentials), intent(in) :: potential

    !> Electrons in each atomic orbital
    real(dp), intent(in) :: qOrb(:,:,:)

    !> reference charges
    real(dp), intent(in) :: q0(:,:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: overSparse(:)

    !> Container for SK Hamiltonian integrals
    type(TSlakoCont), intent(in) :: skHamCont

    !> Container for SK overlap integrals
    type(TSlakoCont), intent(in) :: skOverCont

    !> method for calculating derivatives of S and H0
    type(TNonSccDiff), intent(in) :: nonSccDeriv

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

    !> Electron temperature
    real(dp), intent(in) :: tempElec

    !> Third order SCC interactions
    type(TThirdOrder), allocatable, intent(inout) :: thirdOrd

    !> Charge mixing object
    type(TMixer), intent(inout) :: pChrgMixer

    integer :: iAt, iCart, jAt, iOrb, jOrb

    integer :: nOrbs

    ! maximum allowed number of electrons in a single particle state
    real(dp) :: maxFill

    integer :: nFilled, nEmpty

    integer :: iSCCIter
    logical :: tStopSCC

    ! matrices for derivatives of terms in hamiltonian and outputs
    real(dp), allocatable :: over(:,:), dOverSparse(:,:), dH0Sparse(:,:)
    real(dp), allocatable :: dHam(:,:), dOver(:,:), dH0(:,:), dHamA(:,:), dOverA(:,:), U(:,:) ! xiA(:,:)
    real(dp), allocatable :: dens(:,:), wMatrix(:,:)

    ! overlap derivative terms in potential, omega dS + d(delta q gammma) S
    real(dp), allocatable :: sOmega(:,:,:)

    real(dp) :: qAtom(nAtom)
  ! real(dp) :: dRho(size(over),1)
    real(dp) :: dqIn(nAtom)
    real(dp) :: dqOut(nAtom, 3, nAtom)
    real(dp) :: dqDiff(nAtom), sccErrorQ

    ! eigenvalue weighted vectors
    real(dp), allocatable :: eCiReal(:, :)

    real(dp) :: shift(nAtom)

  ! real(dp), allocatable :: vAt(:,:), vdgamma(:,:,:)
    real(dp), allocatable :: dGammaComplete(:,:,:), dGamma(:,:)
    real(dp), allocatable :: gamma2(:,:), sumGammaTimesDqDa(:), sumDGammaTimesQ(:)

    logical :: tSccCalc, tThirdOrd, tConverged, tMetallic

  ! real(dp), allocatable :: dEi(:,:,:,:)
  ! real(dp), allocatable :: dPsiReal(:,:,:,:)

  ! integer :: fdResults

    ! non-variational part of charge derivative
    real(dp), allocatable :: dqNonVariational(:)

    ! derivative of shift due to external charges w.r.t. coordinates of an atom
    real(dp), allocatable :: extShiftDerivative(:,:)

  ! real(dp) :: dDipole(3)

    ! possibly missing in the 3rd order calculation?
    integer :: iAt1, iAt2
    real(dp), allocatable :: dGamma3(:,:)

    real(dp), allocatable :: testMatrix(:,:)

    real(dp) :: iElement, jElement
    integer :: mu, nu, muBeg, muEnd, nuBeg, nuEnd ! rho, sigma
    real(dp), allocatable :: dqCurrent(:)
    real(dp), allocatable :: gamma3ab(:,:), gamma3ba(:,:)


    tSccCalc = allocated(sccCalc)
    tThirdOrd = allocated(thirdOrd)

    maxFill = 2.0_dp

    nOrbs = size(filling)

  ! allocate(dEi(nOrbs, nAtom, 1, 3))
    allocate(dqNonVariational(nAtom))
    allocate(extShiftDerivative(3, nAtom))

    allocate(eCiReal(size(eigVecsReal,dim=1), size(eigVecsReal,dim=2)))
    do iOrb = 1, size(eigVecsReal,dim=2)
      eCiReal(:,iOrb) = eigVecsReal(:,iOrb) * eigVals(iOrb)
    end do
  ! write (*,*) "eCiReal"
  ! write (*,'(5F9.5)') eCiReal

    shift(:) = potential%intBlock(1, 1, :, 1)
    qAtom(:) = sum(qOrb(:,:,1) - q0(:,:,1), dim=1)
  ! write (*,*) "shift"
  ! write (*,'(F9.5)') shift
  ! write (*,*) "qAtom"
  ! write (*,'(F9.5)') qAtom

  ! allocate(dHamSparse(size(ham,dim=1))
    allocate(over(nOrbs,nOrbs))
    allocate(dOverSparse(size(overSparse),3))
    allocate(dH0Sparse(size(overSparse),3))
    allocate(dHam(nOrbs,nOrbs))
    allocate(dOver(nOrbs,nOrbs))
    allocate(dH0(nOrbs,nOrbs))
    allocate(dOverA(nOrbs,nOrbs))
    allocate(dHamA(nOrbs,nOrbs))
  ! allocate(xiA(nOrbs,nOrbs))
    allocate(U(nOrbs,nOrbs))
    allocate(dens(nOrbs,nOrbs))
    allocate(wMatrix(nOrbs,nOrbs))

    ! terms v S' and v' S
    allocate(sOmega(nOrbs,nOrbs,3))
  ! allocate(vAt(nAtom,1))
  ! allocate(vdgamma(orb%mShell,nAtom,1))

    allocate(dGamma(nAtom,nAtom))
    allocate(sumGammaTimesDqDa(nAtom))
    allocate(sumDGammaTimesQ(nAtom))
    allocate(dqCurrent(nAtom))

  ! if (tThirdOrd) then
  !   allocate(vAt3(nAtom))
  ! ! allocate(vdgamma3(orb%mShell,nAtom))
  ! end if

    nFilled = -1
    do iOrb = 1, nOrbs
      if (filling(iOrb) < epsilon(1.0)) then
        nFilled = iOrb - 1
        exit
      end if
    end do
    if (nFilled < 0) then
      nFilled = nOrbs
    end if

    nEmpty = -1
    do iOrb = 1, nOrbs
      if (abs(filling(iOrb) - maxFill) > epsilon(1.0)) then
        nEmpty = iOrb
        exit
      end if
    end do
    if (nEmpty < 0) then
      nEmpty = 1
    end if

    ! Are there any fractionally occupied orbitals?
    tMetallic = .not. (nFilled == nEmpty - 1)
  ! write(stdOut,*)'Fractionally filled range'
  ! write(stdOut,*) nEmpty, ':', nFilled

    dqOut(:,:,:) = 0.0_dp
  ! dEi(:,:,:,:) = 0.0_dp

    ! 2nd order gamma matrix
    allocate(gamma2(nAtom,nAtom))
    call sccCalc%getAtomicGammaMatrix(gamma2, neighbourList%iNeighbour, img2CentCell)
    do iAt = 1, nAtom
      do jAt = 1, iAt-1
        gamma2(jAt, iAt) = gamma2(iAt, jAt)
      end do
    end do
  ! write (*,*) "gamma2"
  ! write (*,'(2F9.5)') gamma2

    ! derivatives of 2nd order gamma matrices w.r.t. atom coordinates
    allocate(dGammaComplete(nAtom, nAtom, 3))
    call sccCalc%getGammaDeriv(env, species, neighbourList%iNeighbour, img2CentCell, dGammaComplete)
    do iAt = 1, nAtom
      do jAt = 1, iAt-1
        dGammaComplete(jAt, iAt, :) = - dGammaComplete(iAt, jAt, :)
      end do
    end do
  ! write (*,*) "dGammaComplete -- iCart 1"
  ! write (*,'(2F9.5)') dGammaComplete(:,:,1)
  ! write (*,*) "dGammaComplete -- iCart 2"
  ! write (*,'(2F9.5)') dGammaComplete(:,:,2)
  ! write (*,*) "dGammaComplete -- iCart 3"
  ! write (*,'(2F9.5)') dGammaComplete(:,:,3)

    ! 3rd order gamma matrix
    if (tThirdOrd) then
      allocate(gamma3ab(nAtom,nAtom))
      allocate(gamma3ba(nAtom,nAtom))
      call thirdOrd%getGamma3(species, coord, gamma3ab, gamma3ba)
    ! write (*,*) "gamma3ab"
    ! write (*,'(2F9.5)') gamma3ab
    ! write (*,*) "gamma3ba"
    ! write (*,'(2F9.5)') gamma3ba
    end if

    ! derivatives of 3rd order Gamma matrices w.r.t. atom coordinates
    if (tThirdOrd) then
      allocate(dGamma3(nAtom, nAtom))
    end if

  ! write (*,*) "potential intBlock"
  ! write (*,'(4F8.5)') potential%intBlock

    call unpackHS(over, overSparse, neighbourList%iNeighbour, nNeighbourSK, denseDesc%iAtomStart,&
        & iSparseStart, img2CentCell)
    do mu = 1, nOrbs
      do nu = 1, mu-1
        over(nu,mu) = over(mu,nu)
      end do
    end do
  ! write (*,*) "over"
  ! write (*,'(5F9.5)') over

    ! density matrix
    dens = 0._dp
    do mu = 1, nOrbs
      do nu = 1, nOrbs
        do iOrb = 1, nOrbs
          dens(mu,nu) = dens(mu,nu) + filling(iOrb) * eigVecsReal(mu,iOrb) * eigVecsReal(nu,iOrb)
        end do
      end do
    end do

    ! Displaced atom to differentiate wrt
    lpAtom: do iAt = 1, nAtom

      call nonSccDeriv%getFirstDerivWhole(dOverSparse, env, skOverCont, coord, species, iAt, orb,&
          & nNeighbourSK, neighbourList%iNeighbour, iSparseStart, img2centcell)

      call nonSccDeriv%getFirstDerivWhole(dH0Sparse, env, skHamCont, coord, species, iAt, orb,&
          & nNeighbourSK, neighbourList%iNeighbour, iSparseStart, img2centcell)

      extShiftDerivative = 0.0_dp
      call sccCalc%getExtShiftDerivative(env, extShiftDerivative(:,:), iAt, coord)

      ! perturbation direction
      lpCart: do iCart = 1, 3

      ! write (stdOut,*) 'Calculating derivative for displacement along ',&
      !     & trim(direction(iCart)),' for atom', iAt

        call unpackHS(dOver, dOverSparse(:,iCart), neighbourList%iNeighbour, nNeighbourSK,&
            & denseDesc%iAtomStart, iSparseStart, img2CentCell)
        call unpackHS(dH0, dH0Sparse(:,iCart), neighbourList%iNeighbour, nNeighbourSK,&
            & denseDesc%iAtomStart, iSparseStart, img2CentCell)
        do mu = 1, nOrbs
          do nu = 1, mu-1
            dOver(nu,mu) = dOver(mu,nu)
            dH0(nu,mu) = dH0(mu,nu)
          end do
        end do
      ! write (*,*) "dOver"
      ! write (*,'(5F9.5)') dOver
      ! write (*,*) "dH0"
      ! write (*,'(5F9.5)') dH0

        ! NEW -- "derivative of overlap" in the basis of molecular orbitals
        dOverA = matmul(transpose(eigVecsReal(:,:)), matmul(dOver, eigVecsReal(:,:)))

        ! derivatives of 2nd order Gamma matrices w.r.t. atom coordinates
        dGamma(:,:) = 0.0_dp
        dGamma(:,iAt) = dGammaComplete(:,iAt,iCart)
        dGamma(iAt,:) = dGammaComplete(:,iAt,iCart)
      ! write (*,*) "dGamma2 matrix"
      ! write (*,'(2F9.5)') dGamma

        do jAt = 1, nAtom
          sumDGammaTimesQ(jAt) = dot_product(dGamma(jAt,:), qAtom(:))
        end do

      ! write (*,*) "sumDGammaTimesQ before 3rd order"
      ! write (*,'(2X,F12.6)') sumDGammaTimesQ

        ! derivatives of 3rd order Gamma matrices w.r.t. atom coordinates
        if (tThirdOrd) then
          call thirdOrd%getGamma3Deriv(species, coord, iAt, iCart, dGamma3)
        ! write (*,*) "dGamma3"
        ! write (*,'(2F12.7)') dGamma3

          ! 1/3 sum_C ( dGamma_CA/dx Dq_C^2 + 2 dGamma_AC/dx Dq_A Dq_C )
          do iAt1 = 1, nAtom
            sumDGammaTimesQ(:) = sumDGammaTimesQ(:) &
                & + dGamma3(iAt1, :) * qAtom(iAt1)**2._dp &
                & + 2._dp * dGamma3(:, iAt1) * qAtom(iAt1) * qAtom(:)
          end do

        ! write (*,*) "sumDGammaTimesQ"
        ! write (*,'(2X,F12.6)') sumDGammaTimesQ
        end if

        sOmega(:,:,:) = 0.0_dp

        do iAt1 = 1, nAtom
          do iAt2 = 1, nAtom
            muBeg = denseDesc%iAtomStart(iAt1)
            muEnd = denseDesc%iAtomStart(iAt1 + 1) - 1
            nuBeg = denseDesc%iAtomStart(iAt2)
            nuEnd = denseDesc%iAtomStart(iAt2 + 1) - 1
            do mu = muBeg, muEnd
              do nu = nuBeg, nuEnd
                ! First part, omega dS
                sOmega(mu, nu, 1) = dOver(mu, nu) * 0.5_dp * (shift(iAt1) + shift(iAt2))
                ! Third part, dOmega from QM/MM interactions * S
                sOmega(mu, nu, 3) = over(mu, nu) * 0.5_dp *&
                    & (extShiftDerivative(iCart,iAt1) + extShiftDerivative(iCart,iAt2))
              end do
            end do
          end do
        end do

      ! write (*,*) "extShiftDerivative"
      ! write (*,'(2X,F12.6)') extShiftDerivative(iCart,:)

      ! write (*,*) "overlap"
      ! write (*,'(5F8.4)') over
      ! write (*,*) "sOmega 1"
      ! write (*,'(5F8.4)') sOmega(:,:,1)
      ! write (*,*) "sOmega 3"
      ! write (*,'(5F8.4)') sOmega(:,:,3)

        ! NEW -- W matrix from Nishimoto
        wMatrix(:,:) = - 0.5_dp * matmul(dens, matmul(dOver, dens))

        ! non-variational part of the charge change due to basis derivatives
        do jAt = 1, nAtom
          muBeg = denseDesc%iAtomStart(jAt)
          muEnd = denseDesc%iAtomStart(jAt + 1) - 1
          ! sum sum (dens * dOver)
          dqNonVariational(jAt) = sum(dens(muBeg:muEnd,:) * dOver(muBeg:muEnd,:))&
          ! sum sum (W-matrix * over)
                              & + sum(wMatrix(muBeg:muEnd,:) * over(muBeg:muEnd,:))
        end do

      ! write (*,*) "dqNonVariational TK"
      ! write (*,'(F8.5)') dqNonVariational

        call reset(pChrgMixer, nAtom)
        dqIn(:) = 0.0_dp

      ! write (stdOut, "(1X,A,T12,A)") 'SCC Iter' , 'Error'

        iSCCIter = 1
        tStopSCC = .false.
        lpSCC: do while (iSCCiter <= maxSccIter)

          if (iSCCiter > 1) then

            dqCurrent(:) = dqIn(:) + dqNonVariational(:)
            do jAt = 1, nAtom
              sumGammaTimesDqDa(jAt) = dot_product(gamma2(jAt,:), dqCurrent)
              if (tThirdOrd) then
                sumGammaTimesDqDa(jAt) = sumGammaTimesDqDa(jAt)&
                    & + 2._dp / 3._dp * dot_product(gamma3ab(jAt,:), qAtom(:)) * dqCurrent(jAt)&
                    & + 2._dp / 3._dp * dot_product(gamma3ab(jAt,:), dqCurrent(:)) * qAtom(jAt)&
                    & + 2._dp / 3._dp * dot_product(gamma3ba(jAt,:), dqCurrent(:) * qAtom(:))
              end if
            end do

          else
            ! is a better initialization in the first iteration possible?
            sumGammaTimesDqDa(:) = 0._dp
          end if

          do iAt1 = 1, nAtom
            do iAt2 = 1, nAtom
              muBeg = denseDesc%iAtomStart(iAt1)
              muEnd = denseDesc%iAtomStart(iAt1 + 1) - 1
              nuBeg = denseDesc%iAtomStart(iAt2)
              nuEnd = denseDesc%iAtomStart(iAt2 + 1) - 1
              do mu = muBeg, muEnd
                do nu = nuBeg, nuEnd
                  sOmega(mu, nu, 2) = over(mu, nu) * ( &
                      ! (Delta q) * d gamma / dx
                      & 0.5_dp * (sumDGammaTimesQ(iAt1) + sumDGammaTimesQ(iAt2)) +&
                      ! gamma * d (Delta q) / dx
                      & 0.5_dp * (sumGammaTimesDqDa(iAt1) + sumGammaTimesDqDa(iAt2)) )
                end do
              end do
            end do
          end do
         
        ! write (*,*) "sOmega 2"
        ! write (*,'(5F8.4)') sOmega(:,:,2)

          ! derivative of Hamiltonian in the basis of atomic functions
          dHam(:,:) = dH0 + sOmega(:,:,1) + sOmega(:,:,2) + sOmega(:,:,3)

        ! write (*,*) "dH0"
        ! write (*,'(5F8.4)') testMatrix

        ! write (*,*) "dHam without sOmega 3"
        ! write (*,'(5F8.4)') testMatrix

        ! write (*,*) "dHam with sOmega 3"
        ! write (*,'(5F8.4)') dHam

          ! NEW -- derivative of Hamiltonian in the basis of molecular orbitals
          dHamA = matmul(transpose(eigVecsReal(:,:)), matmul(dHam, eigVecsReal(:,:)))

        ! write (*,*) "dHamA"
        ! write (*,'(5F8.4)') dHamA

        ! ! NEW -- Xi matrix (possible problem of degenerate orbitals
        ! !         -- if necessary, apply the known solution...)
        ! xiA(:,:) = 0.
        ! do iOrb = 1, nOrbs ! 1, nFilled ! occupied
        !   do jOrb = 1, nOrbs ! nEmpty, nOrbs ! virtual
        !     xiA(iOrb, jOrb) = dHamA(iOrb, jOrb) * (filling(jOrb) - filling(iOrb)) &
        !     &                                   / (eigVals(jOrb) - eigVals(iOrb)) &
        !     & - dOverA(iOrb, jOrb) * (filling(jOrb) * eigVals(jOrb) - &
        !     &   filling(iOrb) * eigVals(iOrb)) / (eigVals(jOrb) - eigVals(iOrb))
        !   end do
        ! end do
        ! xiA(:,:) = 1. / maxFill * xiA

        ! write (*,*) "XiA"
        ! write (*,'(5F8.4)') xiA

          ! NEW -- perturbation matrix U
          do iOrb = 1, nOrbs
            U(iOrb,iOrb) = -0.5_dp * dOverA(iOrb,iOrb)
          end do

          do iOrb = 1, nOrbs
            do jOrb = 1, nOrbs
              if (iOrb .ne. jOrb) then
                U(iOrb,jOrb) = (dHamA(iOrb,jOrb) - eigVals(jOrb) * dOverA(iOrb,jOrb)) &
                           & / (eigVals(jOrb) - eigVals(iOrb))
              end if
            end do
          end do

        ! write (*,*) "perturbation matrix U"
        ! write (*,'(5F8.4)') U

          ! NEW -- charge derivatives
          dqOut(:,iCart,iAt) = 0._dp
          do jAt = 1, nAtom
            muBeg = denseDesc%iAtomStart(jAt)
            muEnd = denseDesc%iAtomStart(jAt + 1) - 1
            do mu = muBeg, muEnd
              do nu = 1, nOrbs
                do iOrb = 1, nFilled
                  do jOrb = nEmpty, nOrbs
                    iElement = eigVecsReal(mu,iOrb) * eigVecsReal(nu,jOrb) * over(mu,nu)
                    jElement = eigVecsReal(mu,jOrb) * eigVecsReal(nu,iOrb) * over(mu,nu)
                    dqOut(jAt,iCart,iAt) = dqOut(jAt,iCart,iAt) +&
                        & filling(iOrb) * U(jOrb,iOrb) * (iElement + jElement)
                  end do
                end do
              end do
            end do
          end do

          dqDiff(:) = dqOut(:,iCart,iAt) - dqIn(:)

          sccErrorQ = maxval(abs(dqDiff))

        ! do jAt = 1, nAtom
        !   write (stdOut, '(A,I3,3F11.6)') "iter dQ/da ", jAt,&
        !       & dqOut(jAt, iCart, iAt) + dqNonVariational(jAt)
        ! end do
        ! write(stdOut,"(1X,I0,T10,E20.12)") iSCCIter, sccErrorQ
          tConverged = (sccErrorQ < sccTol)

          if ((.not. tConverged) .and. iSCCiter /= maxSccIter) then
            if (iSCCIter == 1) then
              dqIn(:) = dqOut(:, iCart, iAt)
            else
              call mix(pChrgMixer, dqIn, dqDiff)
            end if
          end if

          if (tConverged) then
            exit lpSCC
          end if

          iSCCIter = iSCCIter +1

        end do lpSCC

        dqOut(:, iCart, iAt) = dqOut(:, iCart, iAt) + dqNonVariational(:)

      end do lpCart

    end do lpAtom

    write (stdOut, *)
    write (stdOut, *) 'Charge derivatives TK'
    do iAt = 1, nAtom
      write (stdOut,"(A,I0)") '/d Atom_', iAt
      do jAt = 1, nAtom
        write (stdOut, '(I3,3F11.6)') jAt, -dqOut(jAt, :, iAt)
      end do
      write (stdOut, *)
    end do
    write (stdOut, *)

  end subroutine dPsidx_TK

end module dftbp_perturbxderivs_tk
