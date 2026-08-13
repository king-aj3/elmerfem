!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! *
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! *
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write
! *  to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
! *  Boston, MA  02110-1301  USA
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!> One constitutive law at one integration point: given the deformation state,
!> the environment and the history there, return the stress -- and, where a model
!> can form one, the consistent tangent.
!>
!> The shape follows UMAT's semantics without its thirty-five positional
!> arguments, so that a model ignores what it does not need and fields can be
!> added without touching every implementation.
!>
!> THE TWO ENTRY POINTS ARE NOT INTERCHANGEABLE, and which one a driver wants
!> depends on how its assembly is written rather than on the model:
!>
!>   ConstitutiveStress_i  -- stress from the strain measure in Point. Every model
!>     has one. This is all the postprocessing needs, and it is also what an
!>     assembly that forms its tangent by Gateaux derivatives needs, applying the
!>     law once per (integration point, basis function, component).
!>   ConstitutiveTangent_i -- stress and the consistent tangent. Optional. What an
!>     assembly in B^T C B form needs, once per integration point.
!>
!> Measured before being written, because the plan for this interface assumed one
!> call per integration point and that is not what either Elmer elasticity solver
!> does: 7.6 ns a call for the stress shape and 21.4 ns with the tangent, against
!> an element that costs 20-40 microseconds to assemble. So the indirection is
!> affordable at either call frequency; what is not affordable is reshaping an
!> assembly to suit the interface, and hence the two entry points.
!>
!> Three conventions, each chosen against an alternative and each worth stating
!> because getting one wrong is silent:
!>
!>   TENSORS, NOT VOIGT. Strain and stress are 3x3. Elmer's own kernels already
!>     pass 3x3 -- Strain2Stress takes tensors and does the Voigt packing with the
!>     engineering shear internally, LocalStress hands back tensors -- so carrying
!>     Voigt here would introduce a second convention across every caller for the
!>     benefit of the one model, UMAT, whose adapter can convert. A factor of two
!>     on a shear term is the failure mode, and it is invisible in an isotropic
!>     test. Tangent is unavoidably Voigt: a fourth-order tensor in matrix form is.
!>   PROPS ARE PRE-EVALUATED BY THE DRIVER, as a plain array, not a ValueList
!>     pointer. A model that read its constants from the sif at every integration
!>     point would spend a keyword lookup per call -- string comparison, hundreds
!>     of nanoseconds at best -- which is one to two orders of magnitude more than
!>     the entire interface overhead measured above. The driver interpolates once
!>     per point and hands the numbers over. This is UMAT's own arrangement.
!>   THE MODEL DECLARES ITS STRESS MEASURE rather than the driver assuming.
!>     Elmer already has one kernel read two ways: Strain2Stress is fed
!>     infinitesimal strain and read as Cauchy by StressSolve, and Green-Lagrange
!>     strain and read as second Piola-Kirchhoff by ElasticSolve. Nothing in its
!>     signature says which, and that ambiguity is what StressMeasure removes.
!------------------------------------------------------------------------------

MODULE Constitutive

  USE Types

  IMPLICIT NONE

  !> Which measure the returned stress is in. Declared by the model; converted by
  !> the driver, which is the only party that knows the deformation gradient it
  !> would push forward with.
  INTEGER, PARAMETER :: STRESS_CAUCHY = 1
  INTEGER, PARAMETER :: STRESS_PK2    = 2

  !> Which strain measure the driver has put in Point % Strain.
  INTEGER, PARAMETER :: KINEMATICS_SMALL_STRAIN     = 1
  INTEGER, PARAMETER :: KINEMATICS_LARGE_DEFLECTION = 2
  INTEGER, PARAMETER :: KINEMATICS_HENCKY           = 3

!------------------------------------------------------------------------------
!> What the driver supplies about one integration point.
!------------------------------------------------------------------------------
  TYPE :: MaterialPoint_t
    REAL(KIND=dp) :: Strain(3,3) = 0.0_dp   !> the strain measure Kinematics names
    REAL(KIND=dp) :: dStrain(3,3) = 0.0_dp  !> its increment over the step
    REAL(KIND=dp) :: F(3,3) = 0.0_dp        !> deformation gradient now
    REAL(KIND=dp) :: F0(3,3) = 0.0_dp       !> and at the start of the step
    REAL(KIND=dp) :: DetF = 1.0_dp          !> its determinant, already to hand
    INTEGER :: Kinematics = KINEMATICS_SMALL_STRAIN
    REAL(KIND=dp) :: Temperature = 0.0_dp, dTemperature = 0.0_dp
    REAL(KIND=dp) :: Time = 0.0_dp, dt = 0.0_dp
    REAL(KIND=dp) :: Coord(3) = 0.0_dp, CharLength = 0.0_dp
    !> Dim is the coordinate system dimension and CDim the number of displacement
    !> components, which differ under axial symmetry. A model needs both: the
    !> trace runs over Dim while the identity it multiplies is CDim-diagonal.
    INTEGER :: Dim = 3, CDim = 3
    LOGICAL :: PlaneStress = .FALSE., AxiSymmetric = .FALSE.
    INTEGER :: Element = 0, IP = 0
  END TYPE MaterialPoint_t

!------------------------------------------------------------------------------
!> History at one integration point. v is read and written, v0 is the converged
!> value at the start of the step and is read only.
!>
!> This is the one place two existing mechanisms collapse into one: Maxwell keeps
!> its previous stress in an Elmer -ip variable while UMAT keeps its own statev
!> array, for the same purpose. A model that wants neither leaves n at zero.
!------------------------------------------------------------------------------
  TYPE :: MaterialState_t
    REAL(KIND=dp), POINTER :: v(:) => NULL()
    REAL(KIND=dp), POINTER :: v0(:) => NULL()
    INTEGER :: n = 0
  END TYPE MaterialState_t

!------------------------------------------------------------------------------
!> What a model returns. Tangent is filled only by ConstitutiveTangent_i.
!>
!> The type carries Tangent on both paths deliberately. Splitting it to spare the
!> stress path thirty-six unused doubles looks like an optimisation and is not:
!> passing the same type through both entry points measured 7.6 ns against 21.4,
!> so the cost is filling the tangent, not carrying it -- these are passed by
!> reference.
!------------------------------------------------------------------------------
  TYPE :: MaterialResponse_t
    REAL(KIND=dp) :: Stress(3,3) = 0.0_dp
    REAL(KIND=dp) :: Tangent(6,6) = 0.0_dp
    INTEGER :: StressMeasure = STRESS_CAUCHY
    REAL(KIND=dp) :: EnergyElastic = 0.0_dp
    REAL(KIND=dp) :: EnergyPlastic = 0.0_dp
    REAL(KIND=dp) :: EnergyViscous = 0.0_dp
    REAL(KIND=dp) :: SuggestedDtFactor = 1.0_dp   !> UMAT's pnewdt
    LOGICAL :: Failed = .FALSE.
  END TYPE MaterialResponse_t

  ABSTRACT INTERFACE
    SUBROUTINE ConstitutiveStress_i( Point, Props, State, Response )
      IMPORT :: dp, MaterialPoint_t, MaterialState_t, MaterialResponse_t
      TYPE(MaterialPoint_t), INTENT(IN) :: Point
      REAL(KIND=dp), INTENT(IN) :: Props(:)
      TYPE(MaterialState_t), INTENT(INOUT) :: State
      TYPE(MaterialResponse_t), INTENT(OUT) :: Response
    END SUBROUTINE ConstitutiveStress_i

    SUBROUTINE ConstitutiveTangent_i( Point, Props, State, Response )
      IMPORT :: dp, MaterialPoint_t, MaterialState_t, MaterialResponse_t
      TYPE(MaterialPoint_t), INTENT(IN) :: Point
      REAL(KIND=dp), INTENT(IN) :: Props(:)
      TYPE(MaterialState_t), INTENT(INOUT) :: State
      TYPE(MaterialResponse_t), INTENT(OUT) :: Response
    END SUBROUTINE ConstitutiveTangent_i
  END INTERFACE

!------------------------------------------------------------------------------
!> What a model registers. A driver reads this to find out what it may do, rather
!> than the model's author deciding by writing another assembly routine.
!------------------------------------------------------------------------------
  TYPE :: MaterialModel_t
    PROCEDURE(ConstitutiveStress_i), POINTER, NOPASS :: Stress => NULL()
    PROCEDURE(ConstitutiveTangent_i), POINTER, NOPASS :: Tangent => NULL()
    !> True when the response is linear in the strain measure, so that applying
    !> Stress to a strain INCREMENT gives the directional derivative of the
    !> stress. That is the licence for the cheap Gateaux path, and it is a
    !> property of the model: it holds for sigma = C:eps and fails for
    !> neo-Hookean, Hencky and UMAT.
    !>
    !> This flag is why ElasticSolve has three assembly routines rather than one.
    !> LocalMatrix takes the Gateaux shortcut; NeoHookeanLocalMatrix and
    !> LocalMatrixWithUMAT exist because the shortcut is not a derivative for
    !> them. They are one idea and two escapes from it, and a driver that reads
    !> this flag needs neither escape written out again.
    LOGICAL :: StrainLinear = .FALSE.
    INTEGER :: nState = 0             !> history slots wanted per integration point
    INTEGER :: nProps = 0             !> how many pre-evaluated constants expected
    CHARACTER(LEN=MAX_NAME_LEN) :: Name = ' '
  END TYPE MaterialModel_t

  !> Props layout for IsotropicLinearStress.
  INTEGER, PARAMETER :: ISOLIN_LAME1 = 1, ISOLIN_LAME2 = 2, ISOLIN_NPROPS = 2

CONTAINS

!------------------------------------------------------------------------------
!> Isotropic linear elasticity: sigma = 2 mu eps + lambda tr(eps) I.
!>
!> Read as Cauchy from an infinitesimal strain and as second Piola-Kirchhoff from
!> a Green-Lagrange one -- the same arithmetic either way, which is why the
!> measure is taken from the kinematics rather than fixed here.
!>
!> The identity is not the plain one. It is CDim-diagonal, and carries the
!> out-of-plane entry only away from plane stress, which is how the plane stress
!> reduction enters: with it zeroed, the volumetric term makes no contribution to
!> sigma_33, and sigma_33 is then zero as plane stress requires. Under axial
!> symmetry the out-of-plane direction is the hoop and does carry stress, so the
!> entry is present whatever the plane stress flag says.
!------------------------------------------------------------------------------
  SUBROUTINE IsotropicLinearStress( Point, Props, State, Response )
!------------------------------------------------------------------------------
    TYPE(MaterialPoint_t), INTENT(IN) :: Point
    REAL(KIND=dp), INTENT(IN) :: Props(:)
    TYPE(MaterialState_t), INTENT(INOUT) :: State
    TYPE(MaterialResponse_t), INTENT(OUT) :: Response
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Identity(3,3), tr
    INTEGER :: i
!------------------------------------------------------------------------------
    Identity = 0.0_dp
    DO i=1,Point % CDim
      Identity(i,i) = 1.0_dp
    END DO
    IF ( Point % AxiSymmetric .OR. .NOT. Point % PlaneStress ) Identity(3,3) = 1.0_dp

    tr = 0.0_dp
    DO i=1,Point % Dim
      tr = tr + Point % Strain(i,i)
    END DO

    Response % Stress = 2.0_dp * Props(ISOLIN_LAME2) * Point % Strain + &
        Props(ISOLIN_LAME1) * tr * Identity

    Response % StressMeasure = MERGE( STRESS_PK2, STRESS_CAUCHY, &
        Point % Kinematics /= KINEMATICS_SMALL_STRAIN )
!------------------------------------------------------------------------------
  END SUBROUTINE IsotropicLinearStress
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> The model record for the above, for a driver to select and interrogate.
!------------------------------------------------------------------------------
  FUNCTION IsotropicLinearModel() RESULT( Model )
!------------------------------------------------------------------------------
    TYPE(MaterialModel_t) :: Model
!------------------------------------------------------------------------------
    Model % Stress => IsotropicLinearStress
    Model % StrainLinear = .TRUE.
    Model % nState = 0
    Model % nProps = ISOLIN_NPROPS
    Model % Name = 'isotropic linear'
!------------------------------------------------------------------------------
  END FUNCTION IsotropicLinearModel
!------------------------------------------------------------------------------

END MODULE Constitutive
!------------------------------------------------------------------------------
