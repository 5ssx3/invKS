MODULE Begin_module
   USE constants
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
     !creat eigen-data
     CALL Build_eigen()
     !>>>Finite difference
     CALL init_finite(gap)
     !init FFT
     CALL PlanFFT(nr1,nr2,nr3)
     !recip grid mesh data
     CALL FillQTable()
     WRITE(6,*)'R-space-GRIDS:',nr1,nr2,nr3
     WRITE(6,*)'G-space-GRIDS:',ng1,ng2,ng3
     WRITE(6,*)'K-space-GRIDS:',nk1,nk2,nk3
     IF(kspacing>0.d0)THEN
        WRITE(6,*)'Num of K-used:',nk
     ELSE
        WRITE(6,*)'Only gamma point is used'
     ENDIF
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
       WRITE(*,*) 'Find point operator number:',nsym
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
