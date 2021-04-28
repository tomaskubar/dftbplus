!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module for linear response derivative calculations using perturbation methods
module dftbp_perturbxderivs_qmmm_tk
  use dftbp_accuracy
  use dftbp_assert
  use dftbp_commontypes
  use dftbp_constants
  use dftbp_coulomb
  use dftbp_densedescr
  use dftbp_finitethelper
  use dftbp_globalenv
  use dftbp_mainio
  use dftbp_mixer
  use dftbp_periodic
  use dftbp_scc
  use dftbp_slakocont
  use dftbp_sparse2dense
  use dftbp_thirdorder, only : TThirdOrder

  implicit none

  private
  public :: dPsidxQMMM_TK

  !> Direction labels
  character(len=1), parameter :: direction(3) = ['x','y','z']

contains

  !> Static (frequency independent) perturbation at q=0
  subroutine dPsidxQMMM_TK(filling, eigVals, eigVecsReal, qOrb, q0, overSparse, nAtom, species,&
      & neighbourList, nNeighbourSK, denseDesc, iSparseStart, img2CentCell, coord, sccCalc,&
      & maxSccIter, sccTol, tempElec, thirdOrd, pChrgMixer)

    !> Fillings of unperturbed system
    real(dp), intent(in) :: filling(:)

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(in) :: eigVals(:)

    !> ground state eigenvectors
    real(dp), intent(in) :: eigVecsReal(:,:)

    !> Electrons in each atomic orbital
    real(dp), intent(in) :: qOrb(:,:,:)

    !> reference charges
    real(dp), intent(in) :: q0(:,:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: overSparse(:)

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

    integer :: iCart, jAt, iOrb, jOrb

    integer :: nOrbs

    ! maximum allowed number of electrons in a single particle state
    real(dp) :: maxFill

    integer :: nFilled, nEmpty

    integer :: iSCCIter
    logical :: tStopSCC

    ! matrices for derivatives of terms in hamiltonian and outputs
    real(dp), allocatable :: over(:,:)
    real(dp), allocatable :: dHam(:,:), dHamA(:,:), U(:,:) !, xiA(:,:)

    ! overlap derivative terms in potential, omega dS + d(delta q gammma) S
    real(dp), allocatable :: sOmega(:,:)

    real(dp) :: qAtom(nAtom)
  ! real(dp) :: dRho(size(over),1)
    real(dp) :: dqIn(nAtom)
    real(dp), allocatable :: dqOut(:,:,:) ! (nAtom, 3, nExtCharge)
    real(dp) :: dqDiff(nAtom), sccErrorQ

    ! eigenvalue weighted vectors
  ! real(dp), allocatable :: eCiReal(:, :)

  ! real(dp), allocatable :: vAt(:,:), vdgamma(:,:,:)
    real(dp), allocatable :: gamma2(:,:), sumGammaTimesDqDa(:), sumDGammaTimesQ(:)

    logical :: tSccCalc, tThirdOrd, tConverged, tMetallic

  ! real(dp), allocatable :: dEi(:,:,:,:)
  ! real(dp), allocatable :: dPsiReal(:,:,:,:)

  ! integer :: fdResults

  ! real(dp) :: dDipole(3)

    integer :: iAt1, iAt2

    real(dp), allocatable :: testMatrix(:,:)

    real(dp) :: iElement, jElement
    integer :: mu, nu, muBeg, muEnd, nuBeg, nuEnd ! rho, sigma
    real(dp), allocatable :: gamma3ab(:,:), gamma3ba(:,:)

    ! QM/MM specific
    ! coordinates and charges of MM atoms
    integer :: iExtChg, nExtCharge
    real(dp), allocatable :: extCoord(:,:), extCharge(:), dGammaQMMM(:,:) !,:)


    ! obtain the coordinates and charges of MM atoms from SCC structures
    call sccCalc%getExternalCharges(nExtCharge, extCoord, extCharge)
    if (nExtCharge <= 0) then
      write (*,*) "No MM atoms, nothing to do in dPsidxQMMM."
      return
    end if
  ! write (*,*) "EXTERNAL CHARGES: NUMBER = ", nExtCharge
  ! do iExtchg=1, nExtcharge
  !   write (*,'(4F10.5)') extCoord(:,iExtchg) / AA__Bohr, extCharge(iExtchg)
  ! end do

    ! allocate the array/s that need to know the number of ext. charges
    allocate(dqOut(nAtom, 3, nExtCharge))

    tSccCalc = allocated(sccCalc)
    tThirdOrd = allocated(thirdOrd)

    maxFill = 2.0_dp

    nOrbs = size(filling)

  ! allocate(dEi(nOrbs, nAtom, 1, 3))

    qAtom(:) = sum(qOrb(:,:,1) - q0(:,:,1), dim=1)
  ! write (*,*) "qAtom"
  ! write (*,'(F9.5)') qAtom

  ! allocate(dHamSparse(size(ham,dim=1))
    allocate(over(nOrbs,nOrbs))
    allocate(dHam(nOrbs,nOrbs))
    allocate(dHamA(nOrbs,nOrbs))
  ! allocate(xiA(nOrbs,nOrbs))
    allocate(U(nOrbs,nOrbs))

    ! terms v S' and v' S -- for QM/MM, only the latter
    allocate(sOmega(nOrbs,nOrbs))

    allocate(sumGammaTimesDqDa(nAtom))
    allocate(sumDGammaTimesQ(nAtom))

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

    ! derivatives of QM--MM 1/r w.r.t. coordinates of MM atoms
  ! allocate(dGammaQMMM(3, nAtom, nExtCharge))
  ! call calcInvRPrimeQMMM(nAtom, nExtCharge, coord, extCoord, extCharge, dGammaQMMM)
    allocate(dGammaQMMM(3, nAtom))

    ! 2nd order gamma matrix
    allocate(gamma2(nAtom,nAtom))
    call sccCalc%getAtomicGammaMatrix(gamma2, neighbourList%iNeighbour, img2CentCell)
    do iAt1 = 1, nAtom
      do iAt2 = 1, iAt1-1
        gamma2(iAt2, iAt1) = gamma2(iAt1, iAt2)
      end do
    end do
  ! write (*,*) "gamma2"
  ! write (*,'(2F9.5)') gamma2

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

    call unpackHS(over, overSparse, neighbourList%iNeighbour, nNeighbourSK, denseDesc%iAtomStart,&
        & iSparseStart, img2CentCell)
    do mu = 1, nOrbs
      do nu = 1, mu-1
        over(nu,mu) = over(mu,nu)
      end do
    end do
  ! write (*,*) "over"
  ! write (*,'(5F9.5)') over

    ! Displaced atom to differentiate wrt
    lpAtom: do iExtChg = 1, nExtCharge

      ! any non-variational QM/MM contribution?

      call calcInvRPrimeAsymm(nAtom, coord, nExtCharge, extCoord, extCharge, iExtChg,&
          & dGammaQMMM, tDerivWrtExtCharges=.true.)

      ! perturbation direction
      lpCart: do iCart = 1, 3

      ! write (stdOut,*) 'Calculating derivative for displacement along ',&
      !     & trim(direction(iCart)),' for MM charge number', iExtChg

      ! sumDGammaTimesQ(:) = dGammaQMMM(iCart, :, iExtchg) ! vAt
        sumDGammaTimesQ(:) = dGammaQMMM(iCart, :) ! vAt

      ! write (*,*) "sumDGammaTimesQ"
      ! write (*,'(2X,F12.6)') sumDGammaTimesQ

        sOmega(:,:) = 0.0_dp

        call reset(pChrgMixer, nAtom)
        dqIn(:) = 0.0_dp

      ! write (stdOut, "(1X,A,T12,A)") 'SCC Iter' , 'Error'

        iSCCIter = 1
        tStopSCC = .false.
        lpSCC: do while (iSCCiter <= maxSccIter)

          if (iSCCiter > 1) then

            do jAt = 1, nAtom
              sumGammaTimesDqDa(jAt) = dot_product(gamma2(jAt,:), dqIn)
              if (tThirdOrd) then
                sumGammaTimesDqDa(jAt) = sumGammaTimesDqDa(jAt)&
                    & + 2._dp / 3._dp * dot_product(gamma3ab(jAt,:), qAtom(:)) * dqIn(jAt)&
                    & + 2._dp / 3._dp * dot_product(gamma3ab(jAt,:), dqIn(:)) * qAtom(jAt)&
                    & + 2._dp / 3._dp * dot_product(gamma3ba(jAt,:), dqIn(:) * qAtom(:))
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
                  sOmega(mu, nu) = over(mu, nu) * ( &
                      ! (Delta q) * d gamma / dx
                      & 0.5_dp * (sumDGammaTimesQ(iAt1) + sumDGammaTimesQ(iAt2)) +&
                      ! gamma * d (Delta q) / dx
                      & 0.5_dp * (sumGammaTimesDqDa(iAt1) + sumGammaTimesDqDa(iAt2)) )
                end do
              end do
            end do
          end do
         
        ! write (*,*) "sOmega 2"
        ! write (*,'(5F8.4)') sOmega(:,:)

          ! derivative of Hamiltonian in the basis of atomic functions
          dHam(:,:) = sOmega(:,:)

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
        !     &                                   / (eigVals(jOrb) - eigVals(iOrb))
        !   end do
        ! end do
        ! xiA(:,:) = 1. / maxFill * xiA

        ! write (*,*) "XiA"
        ! write (*,'(5F8.4)') xiA

          ! NEW -- perturbation matrix U
          do iOrb = 1, nOrbs
            U(iOrb,iOrb) = 0._dp
          end do

          do iOrb = 1, nOrbs
            do jOrb = 1, nOrbs
              if (iOrb .ne. jOrb) then
                U(iOrb,jOrb) = dHamA(iOrb,jOrb) / (eigVals(jOrb) - eigVals(iOrb))
              end if
            end do
          end do

        ! write (*,*) "perturbation matrix U"
        ! write (*,'(5F8.4)') U

          ! NEW -- charge derivatives
          dqOut(:,iCart,iExtChg) = 0._dp
          do jAt = 1, nAtom
            muBeg = denseDesc%iAtomStart(jAt)
            muEnd = denseDesc%iAtomStart(jAt + 1) - 1
            do mu = muBeg, muEnd
              do nu = 1, nOrbs
                do iOrb = 1, nFilled
                  do jOrb = nEmpty, nOrbs
                    iElement = eigVecsReal(mu,iOrb) * eigVecsReal(nu,jOrb) * over(mu,nu)
                    jElement = eigVecsReal(mu,jOrb) * eigVecsReal(nu,iOrb) * over(mu,nu)
                    dqOut(jAt,iCart,iExtChg) = dqOut(jAt,iCart,iExtChg) +&
                        & filling(iOrb) * U(jOrb,iOrb) * (iElement + jElement)
                  end do
                end do
              end do
            end do
          end do

          dqDiff(:) = dqOut(:,iCart,iExtChg) - dqIn(:)

          sccErrorQ = maxval(abs(dqDiff))

        ! do jAt = 1, nAtom
        !   write (stdOut, '(A,I3,3F11.6)') "iter dQ/da ", jAt, dqOut(jAt, iCart, iExtChg)
        ! end do
        ! write(stdOut,"(1X,I0,T10,E20.12)") iSCCIter, sccErrorQ
          tConverged = (sccErrorQ < sccTol)

          if ((.not. tConverged) .and. iSCCiter /= maxSccIter) then
            if (iSCCIter == 1) then
              dqIn(:) = dqOut(:, iCart, iExtChg)
            else
              call mix(pChrgMixer, dqIn, dqDiff)
            end if
          end if

          if (tConverged) then
            exit lpSCC
          end if

          iSCCIter = iSCCIter +1

        end do lpSCC

      end do lpCart

    end do lpAtom

    write (stdOut, *)
    write (stdOut, *) 'Charge derivatives TK'
    do iExtChg = 1, nExtCharge
      write (stdOut,"(A,I0)") '/d MMcharge_', iExtChg
      do jAt = 1, nAtom
        write (stdOut, '(I3,3F11.6)') jAt, -dqOut(jAt, :, iExtChg)
      end do
      write (stdOut, *)
    end do
    write (stdOut, *)

  end subroutine dPsidxQMMM_TK

end module dftbp_perturbxderivs_qmmm_tk
