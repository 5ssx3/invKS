MODULE constants
   !##########################################################
   !* CREATED_TIME  : 2015-05-08 20:04:24
   !* AUTHOR        : Yanchao Wang
   !* CHANGE        : Xuecheng Shao
   !* ADD           : Xuecheng Shao
   !* DESCRIPTION   :              
   !     ------                    
   !* REFERENCES    :              
   !     1. Mohr P J, Taylor B N, Newell D B. CODATA 
   !     Recommended Values of the Fundamental Physical Constants: 2010a)[J]. Journal of
   !     Physical and Chemical Reference Data, 2012, 41(4): 043109.
   !
   !     2. http://physics.nist.gov/cuu/Constants/index.html
   !     3. https://en.wikipedia.org/wiki/Atomic_units
   !* LOG           :              
   !     2015-05-08 :              
   !* LAST MODIFIED : 2015-08-03 10:22:12 AM
   !##########################################################
   IMPLICIT NONE
   INTEGER,  PARAMETER      :: I4B = SELECTED_INT_KIND(9)
   INTEGER,  PARAMETER      :: I2B = SELECTED_INT_KIND(4)
   INTEGER,  PARAMETER      :: I1B = SELECTED_INT_KIND(2)
   INTEGER,  PARAMETER      :: SP = KIND(1.0)
   INTEGER,  PARAMETER      :: DP = KIND(1.0D0)
   INTEGER,  PARAMETER      :: SCP = KIND((1.0,1.0))
   INTEGER,  PARAMETER      :: DCP = KIND((1.0D0,1.0D0))
!   INTEGER,  PARAMETER      :: DCP = KIND((1.0D0))
   INTEGER,  PARAMETER      :: LGT = KIND(.TRUE.)
   !REAL(SP), PARAMETER     :: sqrt2=1.41421356237309504880168872420969807856967_SP
   !REAL(SP), PARAMETER     :: euler=0.5772156649015328606065120900824024310422_DP
   REAL(DP), PARAMETER      :: pi=3.141592653589793238462643383279502884197_DP
   REAL(DP), PARAMETER      :: pio2=1.57079632679489661923132169163975144209858_DP
   REAL(DP), PARAMETER      :: twopi=6.283185307179586476925286766559005768394_DP
   COMPLEX(DP), PARAMETER   :: imag = (0.0_dp, 1.0_dp)
   INTEGER(I4B),  PARAMETER :: inputunit=12
   INTEGER(I4B),  PARAMETER :: errunit=13
   INTEGER(I4B),  PARAMETER :: outputunit=14
   INTEGER(I4B),  PARAMETER :: mdposunit=2123
   REAL(DP), PARAMETER      :: rydberg = 13.6058d0
   REAL(DP), PARAMETER      :: bohr2ang = 0.529177208607388d0
   REAL(DP), PARAMETER      :: ang2bohr = 1.d0/0.529177208607388d0
   REAL(DP), PARAMETER      :: hart2ev = 27.2113834279111d0
   REAL(DP), PARAMETER      :: ev2hart = 1.d0/27.2113834279111d0
   REAL(DP), PARAMETER      :: scf_tol=1.0d-6
   REAL(DP), PARAMETER      :: au2ev_force=hart2ev/bohr2ang
   REAL(DP), PARAMETER      :: au2gpa = 2.9421912d4
   REAL(DP), PARAMETER      :: golden = 0.38196601125010515179541316563436_DP
   REAL(DP), PARAMETER      :: vlight = 299792458  ! m/s
   REAL(DP),PARAMETER       :: xtiny = 1.0D-14
   !-----------------------------------------------------------------------
   !-----------------------------------------------------------------------
   ! THE GOLDEN RATIO
   REAL(DP), PARAMETER :: CONST_ME_AU    = 9.10938356D-31
   ! <M_E>              MASS                    (KG)
   REAL(DP), PARAMETER :: CONST_E_AU     = 1.6021766208D-19
   ! <E>                CHARGE                  (C)
   REAL(DP), PARAMETER :: CONST_EH_AU    = 4.359744650D-18
   ! <E_H>              ENERGY                  (J)
   REAL(DP), PARAMETER :: CONST_LEN_AU   = 0.52917721067
   ! <A_0>              LENGTH                  (Å)
   REAL(DP), PARAMETER :: CONST_HBAR_AU  = 1.054571800D-34
   ! <H_BAR>            ACTION                  (J*S)
   REAL(DP), PARAMETER :: CONST_EP_AU    = 27.21138602
   ! <E_H/E>            ELECTRIC POTENTIAL      (V)
   REAL(DP), PARAMETER :: CONST_F_AU     = 8.23872336D-8
   ! <E_H/A_0>          FORCE                   (N)
   REAL(DP), PARAMETER :: CONST_MT_AU    = 1.992851882D-24
   ! <H_BAR/A_0>        MOMENTUM                (KG*M/S)
   REAL(DP), PARAMETER :: CONST_TIME_AU  = 2.418884326509D-17
   ! <H_BAR/E_H>        TIME                    (S)
   REAL(DP), PARAMETER :: CONST_I_AU     = 6.623618183D-3
   ! <E*E_H/H_BAR>      CURRENT                 (A)
   REAL(DP), PARAMETER :: CONST_TEMP_AU  = 3.1577464D5
   !3<E_H/K_B>          TEMPERATURE             (K)
   REAL(DP), PARAMETER :: CONST_P_AU     = 2.9421912D13
   !3<E_H/A_0^3>        PRESSURE                (N/M^2)
   REAL(DP), PARAMETER :: CONST_V_AU     = 2.1876912633D6
   !3<A_0*E_H/H_BAR>    VELOCITY                (A_0*E_H/H_BAR)
   REAL(DP), PARAMETER :: CONST_KE_AU    = 8.9875517873681D9
   !3<K_E>              COULOMB FORCE CONSTANT  (KG*M^3*/S^2/C^2)
   REAL(DP), PARAMETER :: CONST_MU_AU    = 1.660539040D-27
   ! <M_U>              ATOMIC MASS CONSTANT    (KG)
   REAL(DP), PARAMETER :: CONST_MA_AU    = 1.82288854642446D3
   ! <M_U/M_E>
   REAL(DP), PARAMETER :: CONST_KB_SI    = 1.3806488D-23    !J/K
   ! KB
END MODULE
