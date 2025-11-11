PROGRAM Beast
   !##########################BEAST WAR#############################!
   !*BEAST   : real space Kohn-Sham DFT Software                    !
   !*Author  : Qiang Xu                                             !
   !*E-mail  : xq@calypso.cn                                        !
   !*Date    : 2017/07/01                                           !
   !*Diag H  : ARPACK                                               !
   !################################################################!
   USE constants
   USE read_module
   USE Begin_module, ONLY : Initial_Grid
   USE invKS_module
   IMPLICIT NONE
   !INTEGER(I4B) :: Isp !,it1,it2
   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   WRITE(*,*) '============================Real-space invKSDFT==========================='
   WRITE(*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>Get 3D Veff and Ts<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
   !Read file
   CALL read_file('param.in')
   !initialize grids
   CALL Initial_Grid()
   PRINT*,'[TASK] invKS calculations'
   CALL iks_optim_WY()
   !CALL iks_main()
   STOP 
   WRITE(*,*) '=================================Well Done==================================='
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDPROGRAM Beast
