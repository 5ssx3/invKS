MODULE parameters
   USE constants
   IMPLICIT NONE
   !---------------------------------------------------------
   !basic parameters
   INTEGER(I4B)  :: ixc=1        !type of XCDF
   INTEGER(I4B)  :: nfd=8 !the order of finite difference
   INTEGER(I4B)  :: nskip  !skip lines
   INTEGER(I4B)  :: nadds=10 , &  !number of eigenstate we need to add
                &   nev=10  ,  & !number of states
                &   Idiag=1 ,  &
                &   CheM=12 ,  &      !the order of chebyshev we use
                &   CheM0=16,  &       !the free diag 
                &   gridn(3)=-1 ,  &  !grid mesh
                &   kgrid(3)=-1 ,  &  !k-grid mesh
                &   nssp=0           !number of simulation steps
   INTEGER(I4B)  :: IGamma=-1
   !        
   REAL(DP)      :: kspacing=0.5   !input k-spacing 
   REAL(DP)      :: kshift(3)=0.d0   !input k-spacing 
   !
   !For Chebyshev
   REAL(DP)      :: wexict=1.0d0
   !
   CHARACTER(30) :: system_name='BEAST WAR' !system name,no use
   !switch for spin
   INTEGER(I4B)  :: nspin=1   !switch for spin
   !smearing
   INTEGER(I4B) :: nsmear=0         !The order for smearing
   REAL(DP)     :: wsmear=0.1d0     !The width of broadening
   !For WY tolerance
   INTEGER(I4B) :: iopt=3
   REAL(DP) :: Chetol=1e-4
   REAL(DP) :: penLambda=1e-5
   REAL(DP) :: RTOL=1e-5
   REAL(DP) :: ETOL=5e-6
   INTEGER(I4B) :: Isym=2 !Isym<0:no 0:time_inv 1:point sym 2:k-point + time 
   !Cheby params
   LOGICAL :: LRROrthNorm=.FALSE.
   LOGICAL :: lfast = .TRUE.
   !---------------------------------------------------------
ENDMODULE parameters
