MODULE Begin_module
   USE constants
#ifdef MPI
   USE smpi_math_module
#endif
   IMPLICIT NONE
CONTAINS
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  SUBROUTINE Initial_Grid()
     USE math ,ONLY : lat2matrix,inv_33,det,dir2car
     USE struct_module , ONLY : lat_mat,lat_para,recip_lat, reclat_para &
                            & , volume, dvol &
                            & , struct
     USE parameters, ONLY : nfd,kspacing
     USE grid_module , ONLY : Build_rgrid,Build_kgrid &
                          &  ,Build_eigen,nr1,nr2,nr3 &
                          & ,ng1,ng2,ng3,nk1,nk2,nk3,nk,gap & 
                          &,FillQTable !,FillRTable &
                          !&,build_ISO_sphere_grid
     USE FOURIER, ONLY : FFT,PlanFFT,CleanFFT
     USE finite_module , ONLY : init_finite
     IMPLICIT NONE
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     !set cell
     CALL lat2matrix(lat_para,lat_mat,2)
     recip_lat(:,:)=2.d0*pi*TRANSPOSE(inv_33(lat_mat))
     CALL lat2matrix(reclat_para,recip_lat,2)
     volume = ABS(det(lat_mat))
     !find point-symmetry
     CALL Find_Symmetry()
     !creat grids
     CALL Build_rgrid()
     !creat k-points
     CALL Build_kgrid()
#ifdef MPI
      CALL smpi_init_2D()
      !debug
      CALL MPI_Barrier(parallel%comm, mpinfo)
      WRITE(6, '(A,I0,A,I0,A,I0)') "[DEBUG] Process ", parallel%myid, &
                                    ": commx_numprocs = ", parallel%commx_numprocs, &
                                    ", commx_myid = ", parallel%commx_myid
      WRITE(6, '(A,I0,A,I0,A,I0)') "[DEBUG] Process ", parallel%myid, &
                                    ": commy_numprocs = ", parallel%commy_numprocs, &
                                    ", commy_myid = ", parallel%commy_myid
      WRITE(6, '(A,I0,A,I0,A,I0)') "[DEBUG] Process ", parallel%myid, &
                                    ": rankx = ", parallel%rankx, ", ranky = ", parallel%ranky
      WRITE(*, '(A,I0,A,I0,I0,I0,I0)') &
         '[DEBUG] Process ', parallel%myid, &
         ': comm/comm2d/commx/commy=', parallel%comm, &
         parallel%comm2d, parallel%commx, parallel%commy
      !PRINT*,'rank2sum=',parallel%comm2d_rank2sum
      CALL MPI_Barrier(parallel%comm, mpinfo)
#endif
     !creat eigen-data
     CALL Build_eigen()
     !>>>Finite difference
     CALL init_finite(gap)
     !init FFT
     CALL PlanFFT(nr1,nr2,nr3)
     !recip grid mesh data
     CALL FillQTable()
#ifdef MPI
   IF (parallel%isroot) THEN
#endif
      WRITE(6,*)'R-space-GRIDS:',nr1,nr2,nr3
      WRITE(6,*)'G-space-GRIDS:',ng1,ng2,ng3
      WRITE(6,*)'K-space-GRIDS:',nk1,nk2,nk3
      IF(kspacing>0.d0)THEN
         WRITE(6,*)'Num of K-used:',nk
      ELSE
         WRITE(6,*)'Only gamma point is used'
      ENDIF
#ifdef MPI
   ENDIF
#endif
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE Initial_Grid
  !------------------inichrg---------------------------
  SUBROUTINE Find_Symmetry()
     USE parameters , ONLY : Isym
     USE math , ONLY : Det
     USE struct_module , ONLY : naty,struct,natom,lat_mat &
                 &, Opsym ,Otrans, nsym, c_i , num_t !, Odet
     USE Cellsym_module , ONLY : GetStructOp
     IMPLICIT NONE
     REAL(DP),PARAMETER :: opoint_tol=0.01 &
                   &,      trans_tol=0.01
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     IF(Isym>0)THEN
         CALL GetStructOp(naty,struct%nati,natom,TRANSPOSE(lat_mat), &
      &       struct%posdir, opoint_tol, trans_tol, &
      &       Opsym ,Otrans, num_t, c_i, nsym)
#ifdef MPI
      IF (parallel%isroot) THEN
#endif
       WRITE(*,*) 'Find point operator number:',nsym
#ifdef MPI
      ENDIF
#endif
     ELSE
       Otrans(:,:)=0._DP
       Opsym(:,:,:)=0._DP
       Opsym(1,1,:)=1._DP
       Opsym(2,2,:)=1._DP
       Opsym(3,3,:)=1._DP
       nsym=1
     ENDIF
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE Find_Symmetry
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDMODULE Begin_module
