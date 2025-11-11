MODULE Fourier
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!    MODULE Fourier
!       |_SUBROUTINE PlanFFT
!       |_SUBROUTINE PlanFST
!       |_SUBROUTINE GetFFTDims
!       |_SUBROUTINE GetFFTComplexDims
!       |_INTERFACE FFT
!         |_FUNCTION ForwardFFT_4D (Private)
!         |_FUNCTION BackFFT_4D (Private)
!         |_FUNCTION ForwardFFT_3D (Private)
!         |_FUNCTION BackFFT_3D (Private)
!       |_FUNCTION ForwardFST
!       |_FUNCTION BackFST
!       |_SUBROUTINE CleanFFT
!
! DESCRIPTION:
!   This module interfaces with the Fastest Fourier Transform in the West 
!   (FFTW) public library to provide Fourier transform facilities for our 
!   quantities. Each Fourier transform has to be planned for first (PlanFFT) 
!   then executed as many times as necessary (FFT() ) and finally cleaned up 
!   (free up the memory) using CleanFFT.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!   You are encouraged to consult the FFTW3 manual online at 
!   http://www.fftw.org
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/20/2003  File created.  (VLL)
!   12/15/2003  Changed INTEGER*8 to INTEGER(KIND=8) to make the compiler happy
!               Also reformatted a bunch of stuff, made blurbs (GSH)
!
!------------------------------------------------------------------------------
                             ! << GLOBAL >>
USE CONSTANTS, ONLY : DP , I4B         

IMPLICIT NONE

INCLUDE 'fftw3.inc'

PRIVATE :: &
  ForwardFFT_4D, & ! Use FFT
  BackFFT_4D, &    ! Use FFT
  ForwardFFT_3D, & ! Use FFT
  BackFFT_3D !  , & ! Use FFT
!  Forward_3D_c2cFFT ,& !Use 
!  Back_3D_c2cFFT !Use

REAL(kind=DP), DIMENSION(:,:,:), ALLOCATABLE, SAVE, PRIVATE :: &
  reelRA                        ! This is a permanent but local array to 
                                ! transform data without losing it. 

COMPLEX(kind=DP), DIMENSION(:,:,:), ALLOCATABLE, SAVE, PRIVATE :: &
    cplxRA   &                     ! This one is its complex counterpart. 
 &, Cin_mat  &                     ! RA = array. Get it?
 &, Cout_mat
INTEGER(kind=8), SAVE, PRIVATE :: & 
  planRtoC, planCtoR, &         ! I'm not sure what this is.
  planRtoR_F, planRtoR_B ,&
  planC2C,planC2C_inv

INTEGER(I4B), SAVE :: &
  offset                        ! parity of reelRA's X-size.

! This interface picks the right transform to perform based on the nature of
! the incomming array: if it's real the FFT is done forward, if complex the 
! back transform is done. All the calls in OFDFT should be of this type: 
! FFT(f).
INTERFACE FFT
  MODULE PROCEDURE ForwardFFT_4D
  MODULE PROCEDURE BackFFT_4D
  MODULE PROCEDURE ForwardFFT_3D
  MODULE PROCEDURE BackFFT_3D
END INTERFACE

INTERFACE CFFT
  MODULE PROCEDURE Forward_3D_c2cFFT
  MODULE PROCEDURE Back_3D_c2cFFT
ENDINTERFACE
CONTAINS

SUBROUTINE PlanFFT(dimX,dimY,dimZ)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This is the initialization procedure that first gets the system name as is
!   called as an argument to OFDFT, and turns it into the various input file
!   names.  Then, it calls all the programs necessary to set variables to 
!   default values, then reads the geometry file to get all the variables sets
!   to the correct values.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!   reelRA, cplxRA, offset
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/20/2003  File created.  (VLL)
!
!------------------------------------------------------------------------------
  USE constants
  USE parameters , ONLY : finite_order=>nfd
  IMPLICIT NONE
                        !>> EXTERNAL VARIABLES <<!
   INTEGER(I4B) :: iert
  INTEGER(I4B), INTENT(IN) :: &
    dimX, dimY, dimZ       ! The dimensions of the cell to be FFT'd

                       !>> INTERNAL VARIABLES <<! 
                         !>> INITIALIZATION <<! 
                         !>> FUNCTION BODY <<!
  !xq
  CALL cleanFFT()

  ALLOCATE(reelRA(dimX, dimY, dimZ))
  ALLOCATE(cplxRA(dimX/2 + 1, dimY, dimZ))  

  CALL dfftw_plan_dft_r2c_3d(planRtoC, dimX, dimY, dimZ, reelRA, cplxRA,&
                             FFTW_ESTIMATE)
  CALL dfftw_plan_dft_c2r_3d(planCtoR, dimX, dimY, dimZ, cplxRA, reelRA,&
                             FFTW_ESTIMATE)
  offset = MOD(dimX, 2)
  !xq add c2c
  IF(finite_order==0.AND..FALSE.)THEN !nouse now
     ALLOCATE(Cin_mat(dimX,dimY,dimZ))
     ALLOCATE(Cout_mat(dimX,dimY,dimZ))
     !Forward
     CALL dfftw_plan_dft_3d(planC2C,dimX,dimY,dimZ,Cin_mat,Cout_mat,-1,FFTW_ESTIMATE)
     !Backword
     CALL dfftw_plan_dft_3d(planC2C_inv,dimX,dimY,dimZ,Cin_mat,Cout_mat,+1,FFTW_ESTIMATE)
  ENDIF

END SUBROUTINE PlanFFT


SUBROUTINE PlanFST(dimX,dimY,dimZ)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This is the same as PlanFFT, except for the Fast Sine Transform.
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/20/2003  File created.  (VLL)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                        !>> EXTERNAL VARIABLES <<!
  INTEGER(I4B), INTENT(IN) :: &
    dimX, dimY, dimZ       ! The dimensions of the cell to be FFT'd

                       !>> INTERNAL VARIABLES <<! 
  REAL(kind=DP), DIMENSION(dimX, dimY, dimZ) :: &
    reelRA_RtoR
  
                         !>> INITIALIZATION <<! 
                         !>> FUNCTION BODY <<!

  CALL dfftw_plan_r2r_3d(planRtoR_F, dimX, dimY, dimZ, &
                             reelRA_RtoR, reelRA_RtoR, &
                             FFTW_RODFT10, FFTW_RODFT10, FFTW_RODFT10, &
                             FFTW_ESTIMATE)

  CALL dfftw_plan_r2r_3d(planRtoR_B, dimX, dimY, dimZ, &
                             reelRA_RtoR, reelRA_RtoR, &
                             FFTW_RODFT01, FFTW_RODFT01, FFTW_RODFT01, &
                             FFTW_ESTIMATE)

END SUBROUTINE PlanFST


SUBROUTINE GetFFTDims(dimX,dimY,dimZ)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Gets the dimensions of the FFT (real-space part)
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   4/25/2006  Added (GSH)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                        !>> EXTERNAL VARIABLES <<!
  INTEGER(I4B), INTENT(OUT) :: &
    dimX, dimY, dimZ       ! The dimensions of the cell to be FFT'd

                       !>> INTERNAL VARIABLES <<! 
                         !>> INITIALIZATION <<! 
                         !>> FUNCTION BODY <<!

  dimX = SIZE(reelRA,1)
  dimY = SIZE(reelRA,2)
  dimZ = SIZE(reelRA,3)

END SUBROUTINE GetFFTDims


SUBROUTINE GetFFTComplexDims(dimX,dimY,dimZ)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Gets the dimensions of the FFT (reciprocal space part)
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   4/26/2006 Added (GSH)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                        !>> EXTERNAL VARIABLES <<!
  INTEGER(I4B), INTENT(OUT) :: &
    dimX, dimY, dimZ       ! The dimensions of the cell to be FFT'd

                       !>> INTERNAL VARIABLES <<! 
                         !>> INITIALIZATION <<! 
                         !>> FUNCTION BODY <<!

  dimX = SIZE(cplxRA,1)
  dimY = SIZE(cplxRA,2)
  dimZ = SIZE(cplxRA,3)

END SUBROUTINE GetFFTComplexDims


FUNCTION ForwardFFT_4D(array) RESULT(transform)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function is not called directly from the OFDFT code. Use the FFT 
!   interface instead. It performs the transformation of a real 4-dimensional 
!   array into its complex 4-dimensional transform. The first dimension is 
!   halved.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/20/2003  File created.  (VLL)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                       !>> EXTERNAL VARIABLES <<!
  REAL(kind=DP), DIMENSION(:,:,:,:) :: &
    array             ! The array to transform

  COMPLEX(kind=DP), DIMENSION(SIZE(array,1)/2+1,SIZE(array,2), &
                              SIZE(array,3),SIZE(array,4)) :: &
    transform         ! The answer

                       !>> INTERNAL VARIABLES <<! 
  INTEGER(I4B) :: &
    is                ! dummy counter

  REAL(kind=DP) :: &
    normalize         ! normalization constant?

                         !>> INITIALIZATION <<! 
                         !>> FUNCTION BODY <<!
  DO is=1, SIZE(array,4)
    reelRA = array(1:SIZE(array,1), 1:SIZE(array,2), 1:SIZE(array,3), is)
    CALL dfftw_execute(planRtoC)
    transform(1:SIZE(transform,1), 1:SIZE(transform,2), 1:SIZE(transform,3),&
              is) = cplxRA
  END DO !is

  ! The forward transform needs to be renormalized afterwards.
  normalize = REAL(SIZE(array,4),kind=DP) / REAL(SIZE(array),kind=DP)
  transform = transform * normalize

END FUNCTION ForwardFFT_4D


FUNCTION BackFFT_4D(array) RESULT(transform)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function is not called directly from the OFDFT code, but rather 
!   through the FFT interface. It performs the reverse Fourier transform of 
!   a complex function over the half-box in reciprocal space back to real 
!   space. It acts on 4-dimensional arrays, the fourth dimension being spin.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/20/2003  File created.  (VLL)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                      !>> EXTERNAL VARIABLES <<!
  COMPLEX(kind=DP), DIMENSION(:,:,:,:) :: &
    array             ! The array to be back FFT'd

  ! The x-size of the returned array is computed from the reciprocal-space
  ! size. This is ambiguous, as a size of 2k or 2k+1 in real space will give
  ! k+1 in reciprocal space. We solve the problem by storing the parity.
  REAL(kind=DP), DIMENSION(2*(SIZE(array,1)-1)+offset,SIZE(array,2), &
                              SIZE(array,3), SIZE(array,4)) :: &
    transform         ! The answer
 
                       !>> INTERNAL VARIABLES <<! 
  INTEGER(I4B) :: &
    is                ! dummy counter

                        !>> INITIALIZATION <<!   
                        !>> FUNCTION BODY <<!
  DO is=1, SIZE(array,4)
    cplxRA = array(1:SIZE(array,1), 1:SIZE(array,2), 1:SIZE(array,3), is)
    CALL dfftw_execute(planCtoR)
    transform(1:SIZE(transform,1), 1:SIZE(transform,2), 1:SIZE(transform,3),&
              is) = reelRA
  END DO 

END FUNCTION BackFFT_4D


FUNCTION ForwardFFT_3D(nx,ny,nz,array) RESULT(transform)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function is not called directly from the OFDFT code. Use the FFT 
!   interface instead. It performs the transformation of a real 4-dimensional 
!   array into its complex 4-dimensional transform. The first dimension is 
!   halved.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/20/2003  File created.  (VLL)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                       !>> EXTERNAL VARIABLES <<!
  INTEGER(I4B),INTENT(IN) :: nx,ny,nz
  REAL(kind=DP), DIMENSION(nx,ny,nz) :: &
    array             ! The array to transform

  COMPLEX(kind=DP), DIMENSION(nx/2+1,ny,nz) :: &
    transform         ! The answer

                       !>> INTERNAL VARIABLES <<! 
                         !>> INITIALIZATION <<! 
                         !>> FUNCTION BODY <<!
!-----------------------------------------------------------------------
!INTEGER                     :: it1,it2
  reelRA = array
  CALL dfftw_execute(planRtoC)
  transform = cplxRA

  ! The forward transform needs to be renormalized afterwards.
  transform = transform / REAL(SIZE(array),kind=DP)

END FUNCTION ForwardFFT_3D


FUNCTION BackFFT_3D(nx,ny,nz,array) RESULT(transform)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function is not called directly from the OFDFT code, but rather 
!   through the FFT interface. It performs the reverse Fourier transform of a 
!   complex function over the half-box in reciprocal space back to real 
!   space. It acts on 3-dimensional arrays.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/20/2003  File created.  (VLL)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                      !>> EXTERNAL VARIABLES <<!
  INTEGER(I4B),INTENT(IN) :: nx,ny,nz
  COMPLEX(kind=DP), DIMENSION(nx,ny,nz) :: &
    array             ! The array to be back FFT'd

  ! The x-size of the returned array is computed from the reciprocal-space
  ! size. This is ambiguous, as a size of 2k or 2k+1 in real space will give
  ! k+1 in reciprocal space. We solve the problem by assuming an odd real size
  REAL(kind=DP), DIMENSION(2*(nx-1)+offset,ny, &
                              nz) :: &
    transform         ! The answer
 
                       !>> INTERNAL VARIABLES <<! 
                         !>> INITIALIZATION <<!   
                         !>> FUNCTION BODY <<!

  cplxRA = array
  CALL dfftw_execute(planCtoR)
  transform = reelRA
 
END FUNCTION BackFFT_3D


SUBROUTINE ForwardFST(array)
!------------------------------------------------------------------------------
! DESCRIPTION:
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                       !>> EXTERNAL VARIABLES <<!
  REAL(kind=DP), DIMENSION(:,:,:,:) :: &
    array             ! The array to transform

                       !>> INTERNAL VARIABLES <<! 
                         !>> INITIALIZATION <<! 
                         !>> FUNCTION BODY <<!

  CALL dfftw_execute_r2r(planRtoR_F, array(:,:,:,1), array(:,:,:,1))

  ! The forward transform needs to be renormalized afterwards.
  array = array / REAL(SIZE(array)*8,kind=DP)

END SUBROUTINE ForwardFST


SUBROUTINE BackFST(array)
!------------------------------------------------------------------------------
! DESCRIPTION:
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                       !>> EXTERNAL VARIABLES <<!

  REAL(kind=DP), DIMENSION(:,:,:,:) :: &
    array             ! The array to transform

                       !>> INTERNAL VARIABLES <<! 
                         !>> INITIALIZATION <<! 
                         !>> FUNCTION BODY <<!

  CALL dfftw_execute_r2r(planRtoR_B, array(:,:,:,1), array(:,:,:,1))

END SUBROUTINE BackFST


SUBROUTINE CleanFFT
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine is called at the end of the run to free the memory 
!   associated with the plan.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!   reelRA, cplxRA
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/20/2003  File created.  (VLL)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
   INTEGER(I4B)  :: iert
                      !>> EXTERNAL VARIABLES <<!  
                      !>> INTERNAL VARIABLES <<!   
                       !>> INITIALIZATION <<!   
                        !>> FUNCTION BODY <<!

  CALL dfftw_destroy_plan(planRtoC)
  CALL dfftw_destroy_plan(planCtoR)
  IF (ALLOCATED(reelRA)) DEALLOCATE(reelRA)
  IF (ALLOCATED(cplxRA)) DEALLOCATE(cplxRA)

END SUBROUTINE CleanFFT

!by QiangXu for c2c Fourier transform

FUNCTION Forward_3D_c2cFFT(array) RESULT(transform)
   IMPLICIT NONE
                       !>> EXTERNAL VARIABLES <<!
   COMPLEX(kind=DP), DIMENSION(:,:,:) :: &
    array             ! The array to transform
   COMPLEX(kind=DP), DIMENSION(SIZE(array,1),SIZE(array,2), &
                              SIZE(array,3)) :: &
    transform         ! The answer
   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   Cin_mat=array
   CALL dfftw_execute(planC2C)
   transform=Cout_mat
   !normalize
   transform = transform / REAL(SIZE(array),kind=DP)
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDFUNCTION  Forward_3D_c2cFFT


FUNCTION Back_3D_c2cFFT(array,str) RESULT(transform)
   IMPLICIT NONE
                       !>> EXTERNAL VARIABLES <<!
   COMPLEX(kind=DP), DIMENSION(:,:,:) :: &
    array             ! The array to transform
   CHARACTER :: str
   COMPLEX(kind=DP), DIMENSION(SIZE(array,1),SIZE(array,2), &
                              SIZE(array,3)) :: &
    transform         ! The answer
   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   Cin_mat=array
   CALL dfftw_execute(planC2C_inv)
   transform=Cout_mat
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDFUNCTION Back_3D_c2cFFT

END MODULE Fourier
