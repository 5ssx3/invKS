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
#ifdef MPI
   USE smpi_math_module
#endif
   IMPLICIT NONE
   !INTEGER(I4B) :: Isp !,it1,it2
   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#ifdef MPI
   !> Initialize the mpi environment
   CALL smpi_init()
#endif
#ifdef MPI
   IF (parallel%isroot) THEN
#endif
      WRITE(*,*) '============================Real-space invKSDFT==========================='
      WRITE(*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>Get 3D Veff and Ts<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
#ifdef MPI
   ENDIF
#endif
   !Read file
   CALL read_file('param.in')
   !initialize grids
   CALL Initial_Grid()
#ifdef MPI
   IF (parallel%isroot) THEN
#endif
      PRINT*,'[TASK] invKS calculations'
#ifdef MPI
   ENDIF
#endif
   CALL iks_optim_WY()
   !CALL iks_main()
#ifdef MPI
   IF (parallel%isroot) THEN
#endif
      WRITE(*,*) '=================================Well Done==================================='
#ifdef MPI
   ENDIF
   CALL smpi_exit()
#endif
   STOP 
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDPROGRAM Beast
