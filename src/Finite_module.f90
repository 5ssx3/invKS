MODULE finite_module
   !##########################################################!{{{
   !* CREATED_TIME  : 2017-07-07
   !* AUTHOR        : Qiang Xu
   !* DESCRIPTION   :              
   !     ------
   !* REFERENCES    :              
   !     ------
   !* LOG           :              
   !##########################################################!}}}
   USE constants
   USE grid_module , ONLY : nr1,nr2,nr3,nr
   IMPLICIT NONE
   !REAL(DP),ALLOCATABLE :: Laplcoe(:) !laplace coe
   REAL(DP),ALLOCATABLE :: Lapl(:,:)  !laplace total coe
   !REAL(DP),ALLOCATABLE :: Gradcoe(:) !gradient coe
   REAL(DP),ALLOCATABLE :: Grad(:,:)  !total gradient coe
   REAL(DP) :: tBmat(3,3) ! for non-orthrog grids
   !REAL(DP) :: lap_gap(3)
   INTEGER(I4B) :: lap_add(3)
   !transform full
   INTEGER :: cell_mu(3,3)
   REAL(DP) :: cell_factor(3)
   !Interface
CONTAINS
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !###################################################################!
   !       For complex wave function finite difference operator        !
   !###################################################################!
   !------------------------destroy  finite------------------------------
   SUBROUTINE destroy_finite()
      IMPLICIT NONE
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF(ALLOCATED(Lapl)) DEALLOCATE(Lapl)
      IF(ALLOCATED(Grad)) DEALLOCATE(Grad)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE destroy_finite
   !---------------------------initial finite----------------------------
   SUBROUTINE init_finite(h)
      USE parameters , ONLY : finite_order=>nfd
      USE math , ONLY : finite_factor,finite_factor_new,inv_33,Det
      USE struct_module , ONLY : lat_mat,lat_para
      IMPLICIT NONE
      !IN/OUT
      REAL(DP) :: h(3)  !gridsize
      !LOCAL
      REAL(DP),ALLOCATABLE,DIMENSION(:) ::  &
            &   Laplcoe  &
            & , Gradcoe
      REAL(DP) :: Amat(3,3),Bmat(3,3)
      REAL(DP) :: lap_gap(3),factor(6)
      INTEGER(I4B) :: i
      INTEGER(I4B) :: norder
      REAL(DP) :: err(6)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CALL destroy_finite()
      !exit
      IF(finite_order/=0)THEN
         norder=finite_order
      ELSE
         norder=8
         finite_order=8
      ENDIF
      !
      ALLOCATE(Laplcoe(-norder:norder))
      ALLOCATE(Gradcoe(-norder:norder))
      ALLOCATE(Lapl(-norder:norder,6))
      ALLOCATE(Grad(-norder:norder,3))
      !set data
!      CALL finite_factor(1,norder,Gradcoe)
!      CALL finite_factor(2,norder,Laplcoe)
!print*,'Grad coe'
!print*,Gradcoe
!print*,'Laplcoe'
!print*,Laplcoe
      CALL finite_factor_new(1,norder,Gradcoe)
      CALL finite_factor_new(2,norder,Laplcoe)
!print*,'our'
!print*,'Grad coe'
!print*,Gradcoe
!print*,'Laplcoe'
!print*,Laplcoe


      !transform to cartensien coor
      Amat(:,1)=lat_mat(:,1)/lat_para(1) !SQRT(SUM(lat_mat(:,1)**2))
      Amat(:,2)=lat_mat(:,2)/lat_para(2)!SQRT(SUM(lat_mat(:,2)**2))
      Amat(:,3)=lat_mat(:,3)/lat_para(3)!SQRT(SUM(lat_mat(:,3)**2))
      Bmat=inv_33(Amat)
      tBmat=TRANSPOSE(Bmat)
      !tBmat=Bmat
      !grad
      Grad(:,1)=Gradcoe(:)/h(1)
      Grad(:,2)=Gradcoe(:)/h(2)
      Grad(:,3)=Gradcoe(:)/h(3)
      !Lapl         
      !total coefficient in grid
      !full
      CALL trans_mat_full(lat_mat,factor,cell_mu,lap_gap,err)
      DO i=1,3
         Lapl(:,i)=factor(i)*Laplcoe(:)/h(i)**2
         Lapl(:,i+3)=factor(i+3)*Laplcoe(:)/lap_gap(i)**2
      ENDDO
      !deallocate
      DEALLOCATE(Laplcoe,Gradcoe)
!print*,'gap'
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE init_finite
   !-------------------------trans_mat_full------------------------------
   subroutine trans_mat_full(mat,factor,int_miu,lad_gap,err)!{{{
      !##########################################################
      !* CREATED_TIME  : 2018-08-31 
      !* AUTHOR        : Xuecheng Shao & Qiang Xu
      !* CHANGE        : Xuecheng Shao
      !* ADD           : Xuecheng Shao
      !* DESCRIPTION   :              
      !     ------                    
      !* REFERENCES    :              
      !     1.PHYSICAL REVIEW B 78, 075109 (2008) 
      !       "Real-space pseudopotential method for first principles calculations 
      !       of general periodic and partially periodic systems"
      !* LOG           :              
      !     2015-05-10
      !* LAST MODIFIED : 2015-05-10 08:39:18 PM
      !##########################################################
      USE constants,   ONLY : DP,I4B
      use math
      !USE Lapack_module , ONLY : invmat_real
      USE grid_module ,ONLY : gap
      USE struct_module , ONLY : lat_para
      implicit none
      !IN/OUT
      REAL(DP),INTENT(IN) :: mat(3,3) !lattice matrix
      REAL(DP),INTENT(OUT) :: factor(6),lad_gap(3)
      INTEGER(I4B),INTENT(OUT) :: int_miu(3,3)
      REAL(DP),INTENT(OUT) :: err(6)
      !LOCAL
      INTEGER(I4B)                     :: i,j,k,i1,i2,i3,ucount
      !INTEGER(I4B)                     :: offset(3),bound(3),n_one,idx(27)
      INTEGER(I4B) :: intmiu(3,27),n_one
      INTEGER(I4B)                     :: sum_miu(3) !,int_tmp
      real(DP),dimension(3,3)     :: A,F_mat,inva,M_mat,invM,new_r,miu
      real(DP)                    :: f_vec(3),vec_tmp(3,27),scal_one(27)
      real(DP)                    :: b_vec(3),r(3),lgap(3)
      real(DP)                    :: vol,rmod
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      do i=1,3
         A(:,i)=mat(:,i)/lat_para(i)
      ENDDO
      invA=inv_33(A)
      F_mat=matmul(invA,transpose(invA))
      f_vec(1)=F_mat(1,2)+F_mat(2,1)
      f_vec(2)=F_mat(1,3)+F_mat(3,1)
      f_vec(3)=F_mat(2,3)+F_mat(3,2)
      k=0 !count total number
      ucount=0 !count useful number
      do i3=-1,1
      do i2=-1,1
      do i1=-1,1
         k=k+1
         intmiu(1,k)=i1
         intmiu(2,k)=i2
         intmiu(3,k)=i3
         !delete the cell lattice vector
         if ((abs(i1)+abs(i2)+abs(i3))<=1) then
            vec_tmp(:,k)=0.d0
            scal_one(k)=10000.d0
         else
            !r(1:3)=i3*a(:,3)+i2*a(:,2)+i1*a(:,1)
            ucount=ucount+1
            r(1:3)=gap(1)*i1*A(:,1)+gap(2)*i2*A(:,2)+gap(3)*i3*A(:,3)
            vec_tmp(:,k)=r
            scal_one(k)=SQRT(SUM(r**2))
         endif
      enddo
      enddo
      enddo
      n_one=27    ! n_one=27
!
!print*,'lgap'
!print*,scal_one
!print*,vec_tmp(1,:)
      !sort by length
      CALL realInt_sort(n_one,scal_one,intmiu,vec_tmp)
!      print*,'sort'
!      print*,scal_one
!print*,vec_tmp(1,:)
      !STOP
      !sort the directions by length
      !call sort_id(scal_one,n_one,idx)
      !select the direction
      !changed by Qiang Xu
      !offset(:)=2
      !bound=3
      findit:DO i=1,n_one
         new_r(:,1)=vec_tmp(:,i)
         int_miu(1,:)=intmiu(:,i)
         lad_gap(1)=scal_one(i)
      DO j=1,n_one
         new_r(:,2)=vec_tmp(:,j)
         int_miu(2,:)=intmiu(:,j)
         lad_gap(2)=scal_one(j)
      DO k=1,n_one
         new_r(:,3)=vec_tmp(:,k)
         int_miu(3,:)=intmiu(:,k)
         lad_gap(3)=scal_one(k)
         !get real miu
         miu(:,:)=REAL(int_miu(:,:),DP)
         r = miu(1,1)*A(:,1) + miu(1,2) * A(:,2) + miu(1,3)*A(:,3)
         rmod=SQRT(SUM(r*r))
         IF(rmod<xtiny) CYCLE
         miu(1,:)=miu(1,:)/rmod
         r = miu(2,1)*A(:,1) + miu(2,2) * A(:,2) + miu(2,3)*A(:,3)
         rmod=SQRT(SUM(r*r))
         IF(rmod<xtiny) CYCLE
         miu(2,:)=miu(2,:)/rmod
         r = miu(3,1)*A(:,1) + miu(3,2) * A(:,2) + miu(3,3)*A(:,3)
         rmod=SQRT(SUM(r*r))
         IF(rmod<xtiny) CYCLE
         miu(3,:)=miu(3,:)/rmod
         !get M_mat=2*M_mu
         do i1=1,3
            M_mat(1,i1)=2*miu(i1,1)*miu(i1,2)
            M_mat(2,i1)=2*miu(i1,1)*miu(i1,3)
            M_mat(3,i1)=2*miu(i1,2)*miu(i1,3)
         enddo
         
         vol=ABS(Det(M_mat))
!print*,int_miu(:,:)
!print*,miu(:,:)
!print*,'vol',vol
!pause
         IF(vol<1e-7) CYCLE
         !satisfy
         EXIT findit

      ENDDO
      ENDDO
      ENDDO findit
!print*,miu(1,:)
!print*,miu(2,:)
!print*,miu(3,:)
!print*,vol
!STOP

      !get factors ang gaps

!print*,i,j,k
!!print*,
      !vol=ABS(Det(M_mat))
      IF(vol>xtiny)THEN
         !print *,"M",M_mat
         !call gauss_u(M_mat,f_vec,b_vec)
         invM=inv_33(M_mat)
         b_vec(:)=MATMUL(invM,f_vec)
         DO i=1,3
            IF(ABS(b_vec(i))<5e-7) b_vec(i)=0._DP
         ENDDO
         !print *,"b",b_vec
         factor(1)=F_mat(1,1) - sum(b_vec(:)*(miu(:,1)**2))
         factor(2)=F_mat(2,2) - sum(b_vec(:)*(miu(:,2)**2))
         factor(3)=F_mat(3,3) - sum(b_vec(:)*(miu(:,3)**2))
         factor(4)=b_vec(1)
         factor(5)=b_vec(2)
         factor(6)=b_vec(3)
         !gap
         !DO i=1,3
         !   !r(:)= gap(1)*a(:,1)*int_miu(i,1) + &
         !   !    & gap(2)*a(:,2)*int_miu(i,2) + &
         !   !    & gap(3)*a(:,3)*int_miu(i,3)
         !   r(:)=new_r(:,i)
         !   lad_gap(i)=SQRT(SUM(r*r))
         !ENDDO

      ELSE
         !volum is zero
         print*,'transmat_full : invert is not exit!!!'
         STOP
         !CALL transmat(factor,lad_gap,lad)
      ENDIF
      !err analysis
      err(1)=factor(1)*gap(1)**2
      err(2)=factor(2)*gap(2)**2
      err(3)=factor(3)*gap(3)**2
      err(4)=factor(4)*lad_gap(1)**2
      err(5)=factor(5)*lad_gap(2)**2
      err(6)=factor(6)*lad_gap(3)**2
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   END subroutine trans_mat_full
   !--------------------------PARTING LINE-------------------------------
   !####################################################################!
   !             For real/complex finite difference operator            !
   !*Author : Qiang Xu                                                  !
   !*Date   : 2017/08/28                                                !
   !####################################################################!
   !--------------------------PARTING LINE-------------------------------
   SUBROUTINE real_pbc_nabla1(ifun,derf,mgfun)
      USE parameters , ONLY : norder=>nfd
      USE math, ONLY : pullbackPBC
      IMPLICIT NONE
      !
      REAL(DP),INTENT(IN)  :: ifun(nr1,nr2,nr3)
      REAL(DP),INTENT(OUT) :: derf(3,nr)
      REAL(DP),OPTIONAL :: mgfun(nr)
      !
      INTEGER(I4B) :: i,ish,ip,ix,iy,iz,irx,iry,irz
      LOGICAL :: lmgfun
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      lmgfun=PRESENT(mgfun)
      !
      derf(:,:) = 0.d0
      !df/dx
      ip=0
      DO iz=1,nr3
      DO iy=1,nr2
      DO ix=1,nr1
         ip=ip+1
         DO ish=-norder,norder
            !3D map
            CALL pullbackPBC(ix+ish,nr1,irx) !MOD(MOD((ix+ish-1),nr1)+nr1,nr1)+1
            CALL pullbackPBC(iy+ish,nr2,iry) 
            CALL pullbackPBC(iz+ish,nr3,irz)

            derf(1,ip)=derf(1,ip)+Grad(ish,1)*ifun(irx,iy,iz)
            derf(2,ip)=derf(2,ip)+Grad(ish,2)*ifun(ix,iry,iz)
            derf(3,ip)=derf(3,ip)+Grad(ish,3)*ifun(ix,iy,irz)
         ENDDO
         derf(:,ip)=MATMUL(tBmat(:,:),derf(:,ip))
      ENDDO
      ENDDO
      ENDDO
      !mod of grad
      IF(lmgfun)THEN
         DO ip=1,nr
            mgfun(ip)=SQRT(derf(1,ip)**2+derf(2,ip)**2+derf(3,ip)**2)
         ENDDO
      ENDIF
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE real_pbc_nabla1
   !--------------------------PARTING LINE-------------------------------
   SUBROUTINE real_pbc_nabla2(ifun,ofun)!{{{
      !##########################################################
      !* CREATED_TIME  : 2013-03-12
      !* AUTHOR        : Yanchao Wang
      !* CHANGE        : Xuecheng Shao
      !* ADD           : Xuecheng Shao
      !* DESCRIPTION   :              
      !     ------
      !* REFERENCES    :              
      !     ------
      !* LOG           :              
      !     2015-05-08 :              
      !* LAST MODIFIED : 2015-05-10 07:47:12 PM
      !##########################################################
      USE parameters , ONLY : norder=>nfd
      USE math, ONLY : pullbackPBC
      !
      implicit none
      !IN/OUT
      REAL(DP),INTENT(IN) :: ifun(nr1,nr2,nr3)
      REAL(DP),INTENT(OUT) :: ofun(nr1,nr2,nr3)
      !LOCAL
      !real(dp)             :: coeke(-norder:norder,6)
      !
      INTEGER(I4B) :: i,ish,ix,iy,iz
      INTEGER(I4B) :: irx,iry,irz
      !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>
      !
      DO iz=1,nr3
      DO iy=1,nr2
      DO ix=1,nr1
         ofun(ix,iy,iz) = 0.d0
         DO ish = -norder,norder
              CALL pullbackPBC(ix+ish,nr1,irx) 
              CALL pullbackPBC(iy+ish,nr2,iry) 
              CALL pullbackPBC(iz+ish,nr3,irz)
              ofun(ix,iy,iz)=ofun(ix,iy,iz)+&
                          &  Lapl(ish,1)*(ifun(irx,iy,iz))+ &
                          &  Lapl(ish,2)*(ifun(ix,iry,iz))+ &
                          &  Lapl(ish,3)*(ifun(ix,iy,irz))
             

             IF(ABS(Lapl(ish,4))>xtiny)THEN
                 CALL pullbackPBC(ix+cell_mu(1,1)*ish,nr1,irx) 
                 CALL pullbackPBC(iy+cell_mu(1,2)*ish,nr2,iry) 
                 CALL pullbackPBC(iz+cell_mu(1,3)*ish,nr3,irz)
                 ofun(ix,iy,iz)=ofun(ix,iy,iz) + Lapl(ish,4)*ifun(irx,iry,irz)
             ENDIF

             IF(ABS(Lapl(ish,5))>xtiny)THEN
                 CALL pullbackPBC(ix+cell_mu(2,1)*ish,nr1,irx) 
                 CALL pullbackPBC(iy+cell_mu(2,2)*ish,nr2,iry) 
                 CALL pullbackPBC(iz+cell_mu(2,3)*ish,nr3,irz)
                 ofun(ix,iy,iz)=ofun(ix,iy,iz) + Lapl(ish,5)*ifun(irx,iry,irz)
             ENDIF

             IF(ABS(Lapl(ish,6))>xtiny)THEN
                 CALL pullbackPBC(ix+cell_mu(3,1)*ish,nr1,irx) 
                 CALL pullbackPBC(iy+cell_mu(3,2)*ish,nr2,iry) 
                 CALL pullbackPBC(iz+cell_mu(3,3)*ish,nr3,irz)
                 ofun(ix,iy,iz)=ofun(ix,iy,iz) + Lapl(ish,6)*ifun(irx,iry,irz)
             ENDIF
         ENDDO
      ENDDO
      ENDDO
      ENDDO
      !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>
   ENDSUBROUTINE real_pbc_nabla2!}}}
   !----------------------------PARTING LINE-----------------------------
   SUBROUTINE cmplx_pbc_nabla1(ifun,derf)
      USE parameters , ONLY : norder=>nfd
      USE math, ONLY : pullbackPBC
      IMPLICIT NONE
      !
      COMPLEX(DCP),INTENT(IN)  :: ifun(nr1,nr2,nr3)
      COMPLEX(DCP),INTENT(OUT) :: derf(3,nr)
      !
      INTEGER(I4B) :: i,ish,ix,iy,iz,ip,irx,iry,irz
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !
      derf(:,:) = CMPLX(0.d0,0.d0)
      !df/dx
      ip=0
      DO iz=1,nr3
      DO iy=1,nr2
      DO ix=1,nr1
         ip=ip+1
         DO ish=-norder,norder
            CALL pullbackPBC(ix+ish,nr1,irx) !MOD(MOD((ix+ish-1),nr1)+nr1,nr1)+1
            CALL pullbackPBC(iy+ish,nr2,iry) 
            CALL pullbackPBC(iz+ish,nr3,irz)

            derf(1,ip)=derf(1,ip)+Grad(ish,1)*ifun(irx,iy,iz)
            derf(2,ip)=derf(2,ip)+Grad(ish,2)*ifun(ix,iry,iz)
            derf(3,ip)=derf(3,ip)+Grad(ish,3)*ifun(ix,iy,irz)
         ENDDO
         derf(:,ip)=MATMUL(tBmat(:,:),derf(:,ip))
      ENDDO
      ENDDO
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE cmplx_pbc_nabla1
   !----------------------------PARTING LINE-----------------------------
   SUBROUTINE cmplx_pbc_nabla2(ifun,ofun)!{{{
      !##########################################################
      !* CREATED_TIME  : 2013-03-12
      !* AUTHOR        : Yanchao Wang
      !* CHANGE        : Xuecheng Shao
      !* ADD           : Xuecheng Shao
      !* DESCRIPTION   :              
      !     ------
      !* REFERENCES    :              
      !     ------
      !* LOG           :              
      !     2015-05-08 :              
      !* LAST MODIFIED : 2015-05-10 07:47:12 PM
      !##########################################################
      USE parameters , ONLY : norder=>nfd
      USE math, ONLY : pullbackPBC
      !
      implicit none
      !IN/OUT
      COMPLEX(DP),INTENT(IN) :: ifun(nr1,nr2,nr3)
      COMPLEX(DP),INTENT(OUT) :: ofun(nr1,nr2,nr3)
      !LOCAL
      !real(dp)             :: coeke(-norder:norder,6)
      !
      INTEGER(I4B) :: i,ish,ix,iy,iz
      INTEGER(I4B) :: irx,iry,irz
      !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>
      !
      DO iz=1,nr3
      DO iy=1,nr2
      DO ix=1,nr1
         ofun(ix,iy,iz) = 0.d0
         DO ish = -norder,norder
              CALL pullbackPBC(ix+ish,nr1,irx) 
              CALL pullbackPBC(iy+ish,nr2,iry) 
              CALL pullbackPBC(iz+ish,nr3,irz)
              ofun(ix,iy,iz)=ofun(ix,iy,iz)+&
                          &  Lapl(ish,1)*(ifun(irx,iy,iz))+ &
                          &  Lapl(ish,2)*(ifun(ix,iry,iz))+ &
                          &  Lapl(ish,3)*(ifun(ix,iy,irz))
             

             IF(ABS(Lapl(ish,4))>xtiny)THEN
                 CALL pullbackPBC(ix+cell_mu(1,1)*ish,nr1,irx) 
                 CALL pullbackPBC(iy+cell_mu(1,2)*ish,nr2,iry) 
                 CALL pullbackPBC(iz+cell_mu(1,3)*ish,nr3,irz)
                 ofun(ix,iy,iz)=ofun(ix,iy,iz) + Lapl(ish,4)*ifun(irx,iry,irz)
             ENDIF

             IF(ABS(Lapl(ish,5))>xtiny)THEN
                 CALL pullbackPBC(ix+cell_mu(2,1)*ish,nr1,irx) 
                 CALL pullbackPBC(iy+cell_mu(2,2)*ish,nr2,iry) 
                 CALL pullbackPBC(iz+cell_mu(2,3)*ish,nr3,irz)
                 ofun(ix,iy,iz)=ofun(ix,iy,iz) + Lapl(ish,5)*ifun(irx,iry,irz)
             ENDIF

             IF(ABS(Lapl(ish,6))>xtiny)THEN
                 CALL pullbackPBC(ix+cell_mu(3,1)*ish,nr1,irx) 
                 CALL pullbackPBC(iy+cell_mu(3,2)*ish,nr2,iry) 
                 CALL pullbackPBC(iz+cell_mu(3,3)*ish,nr3,irz)
                 ofun(ix,iy,iz)=ofun(ix,iy,iz) + Lapl(ish,6)*ifun(irx,iry,irz)
             ENDIF
         ENDDO
      ENDDO
      ENDDO
      ENDDO
      !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>
   ENDSUBROUTINE cmplx_pbc_nabla2!}}}
   !----------------------------PARTING LINE-----------------------------
   !####################################################################!
   !             For real/complex FFT nabla operators                   !
   !*Author : Qiang Xu                                                  !
   !*Date   : 2017/08/28                                                !
   !####################################################################!
   !----------------------------PARTING LINE-----------------------------
   SUBROUTINE fft_real_nabla1(ifun,derf)
      USE FOURIER
      USE grid_module , ONLY : ng1,ng2,ng3,grid
      IMPLICIT NONE
      REAL(DP),INTENT(IN) :: ifun(nr1,nr2,nr3)
      REAL(DP),INTENT(OUT) :: derf(nr1,nr2,nr3,3)
      !
      COMPLEX(DP),DIMENSION(ng1,ng2,ng3) :: recf,tmp
      INTEGER(I4B) :: ix,iy,iz,ik,I
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !compute grad
      recf=FFT(nr1,nr2,nr3,ifun)
      ! dir 1
      DO ik=1,3
         I=0
         DO iz=1,ng3
         DO iy=1,ng2
         DO ix=1,ng1
            I=I+1
            tmp(ix,iy,iz)=IMAG*grid%gVec(ik,I)*recf(ix,iy,iz)
         ENDDO
         ENDDO
         ENDDO
         derf(:,:,:,ik)=FFT(ng1,ng2,ng3,tmp)
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE fft_real_nabla1
   !----------------------------PARTING LINE-----------------------------
   SUBROUTINE fft_real_nabla2(ifun,ofun)
      USE FOURIER
      USE grid_module , ONLY : ng1,ng2,ng3,grid
      IMPLICIT NONE
      REAL(DP),INTENT(IN) :: ifun(nr1,nr2,nr3)
      REAL(DP),INTENT(OUT) :: ofun(nr1,nr2,nr3)
      !
      COMPLEX(DP) :: recf(ng1,ng2,ng3)
      INTEGER(I4B) :: ix,iy,iz,I
      !-\nabla^2 infun
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      recf=FFT(nr1,nr2,nr3,ifun)
      I=0
      DO iz=1,ng3
      DO iy=1,ng2
      DO ix=1,ng1
         I=I+1
         recf(ix,iy,iz)=recf(ix,iy,iz)*grid%gVec(4,I)**2
      ENDDO
      ENDDO
      ENDDO
      ofun=-FFT(ng1,ng2,ng3,recf)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<    
   ENDSUBROUTINE fft_real_nabla2
   !----------------------------PARTING LINE-----------------------------
   SUBROUTINE fft_cmplx_nabla1(ifun,ofun)
      IMPLICIT NONE
      COMPLEX(DP),INTENT(IN) :: ifun(nr1,nr2,nr3)
      COMPLEX(DP),INTENT(OUT) :: ofun(nr1,nr2,nr3,3)
      !
      REAL(DP) :: tmp(nr1,nr2,nr3,3,2)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !real part
      CALL fft_real_nabla1(REAL(ifun),tmp(:,:,:,:,1)) 
      !imag part
      CALL fft_real_nabla1(AIMAG(ifun),tmp(:,:,:,:,2))
      !combine
      ofun(:,:,:,:)=CMPLX(tmp(:,:,:,:,1),tmp(:,:,:,:,2))
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE fft_cmplx_nabla1
   !----------------------------PARTING LINE-----------------------------
   SUBROUTINE fft_cmplx_nabla2(ifun,ofun)
      USE grid_module , ONLY : grid,ng1,ng2,ng3
      USE Fourier
      IMPLICIT NONE
      COMPLEX(DP),INTENT(IN) :: ifun(nr1,nr2,nr3)
      COMPLEX(DP),INTENT(OUT) :: ofun(nr1,nr2,nr3)
      !
      REAL(DP) :: tmp(nr1,nr2,nr3,2)
      !INTEGER(I4B) :: ix,iy,iz,I
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !real part
      CALL fft_real_nabla2(REAL(ifun),tmp(:,:,:,1)) 
      !imag part
      CALL fft_real_nabla2(AIMAG(ifun),tmp(:,:,:,2)) 
      !Conbine
      ofun(:,:,:)=CMPLX(tmp(:,:,:,1),tmp(:,:,:,2))
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE fft_cmplx_nabla2
   !----------------------------PARTING LINE-----------------------------
   !----------------------------PARTING LINE-----------------------------
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDMODULE finite_module
