MODULE matvec_module
   USE constants
   interface matvec!{{{
      module procedure cmplx_matvec
      module procedure real_matvec
   end interface!}}}
CONTAINS
!-----------------------for cmplex matver------------------------
  SUBROUTINE cmplx_matvec(Ik,veff1d,p,q,dimen)
     IMPLICIT NONE
     !
     INTEGER(I4B),INTENT(IN) :: Ik,dimen
     REAL(DP),INTENT(IN) :: veff1d(dimen)
     COMPLEX(DCP),INTENT(IN) :: p(dimen)
     COMPLEX(DCP),INTENT(OUT) :: q(dimen)
     !> local
     COMPLEX(DCP)     :: kp(dimen) !kinetic operator phi
     INTEGER(I4B) :: I,Ix,Iy,Iz
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     !nonlocal part
     !CALL cmplx_nlocmatvec(Ik,p,q)
     !calculate the kinetic part(see Bloch Theorm.)
     CALL cmplx_keop(p,Ik,kp)
     !collection
     !q(:)=q(:)+kp(:)+veff1d(:)*p(:)
     q(:)=kp(:)+veff1d(:)*p(:)
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE cmplx_matvec
!-----------------------DIVIDER-LINE--------------------------
  SUBROUTINE real_matvec(veff1d,p,q,dimen)
     USE finite_module, ONLY : real_pbc_nabla2
     use grid_module,only:grid
     IMPLICIT NONE
     !
     INTEGER(I4B),INTENT(IN) :: dimen
     REAL(DP),INTENT(IN) :: veff1d(dimen)
     REAL(DP),INTENT(IN) :: p(dimen)
     REAL(DP),INTENT(OUT) :: q(dimen)
     !> local
     REAL(DP)     :: kp(dimen) !kinetic operator phi
     INTEGER(I4B) :: I,Ix,Iy,Iz
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     !nonlocal part
     !CALL real_nlocmatvec(p,q)
     !calculate the kinetic part(kp=-0.5*nabla2*p)
     CALL real_pbc_nabla2(p,kp)
     kp(:)=-0.5_DP*kp(:)
     !collection
     !q(:)=q(:)+kp(:) +veff1d(:)*p(:)
     q(:)=kp(:) +veff1d(:)*p(:)
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE real_matvec
!!-----------------------DIVIDER-LINE--------------------------
!  SUBROUTINE cmplx_nlocmatvec(Ik,p,q)
!     USE pspot_module , ONLY : psp ,max_nproj
!     USE struct_module , ONLY : naty,natom,struct
!     USE nlpot_module , ONLY : nlpot
!     USE grid_module ,ONLY : dvol
!     IMPLICIT NONE
!     !IN/OUT
!     INTEGER(I4B),INTENT(IN) :: Ik
!     COMPLEX(DCP),INTENT(IN)     :: p(:)
!     COMPLEX(DCP),INTENT(OUT)    :: q(:)
!     !LOCAL
!     INTEGER(I4B) :: Ity,Ia,Ipj,Ip,Id
!     COMPLEX(DCP) :: dots(max_nproj,natom),tmp0
!     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     !
!     q(:)=0.d0
!     IF(max_nproj==0)RETURN
!     !
!     dots(:,:)=0.d0
!     !all type
!     DO Ity=1,naty
!        IF(psp(Ity)%nproj==0) CYCLE
!        !all atom
!        DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
!           !all projectors
!           DO Ipj=1,psp(Ity)%nproj
!              tmp0=0.d0
!              !all points
!              DO Ip=1,nlpot(Ia)%npts
!                 Id=nlpot(Ia)%Id(Ip)
!                 !<proj_lm|p>
!                 tmp0=tmp0+CONJG(nlpot(Ia)%proj_phs(Ip,Ipj,Ik))*p(Id)
!              ENDDO
!              !save dots
!              ! dots(Ipj,Ia)=tmp0*psp(Ity)%D0(Ipj,Ipj)
!              dots(1:psp(Ity)%nproj,Ia)=dots(1:psp(Ity)%nproj,Ia)  &
!                       & + tmp0*psp(Ity)%Dij(:,Ipj)
!           ENDDO
!        ENDDO
!     ENDDO
!     !scale
!     dots(:,:)=dots(:,:)*dvol
!     !Multiply the nonlocal vectors
!     DO Ity=1,naty
!        !all atom
!        DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
!           !all projectors
!           DO Ipj=1,psp(Ity)%nproj
!              tmp0=dots(Ipj,Ia)
!              !all points
!              DO Ip=1,nlpot(Ia)%npts
!                 Id=nlpot(Ia)%Id(Ip)
!                 !q=\SUM_lm{dots_lm*|proj_lm>}
!                 q(Id)=q(Id)+tmp0*nlpot(Ia)%proj_phs(Ip,Ipj,Ik)
!              ENDDO
!           ENDDO
!        ENDDO
!     ENDDO
!     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!  ENDSUBROUTINE cmplx_nlocmatvec
!!-----------------------DIVIDER-LINE--------------------------
!  SUBROUTINE real_nlocmatvec(p,q)
!     USE pspot_module , ONLY : psp ,max_nproj
!     USE struct_module , ONLY : naty,natom,struct
!     USE nlpot_module , ONLY : nlpot
!     USE grid_module ,ONLY : dvol
!     IMPLICIT NONE
!     !IN/OUT
!     REAL(DP),INTENT(IN)     :: p(:)
!     REAL(DP),INTENT(OUT)    :: q(:)
!     !LOCAL
!     INTEGER(I4B) :: Ity,Ia,Ipj,Ip,Id
!     REAL(DP) :: dots(max_nproj,natom),tmp0
!     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     !
!     q(:)=0.d0
!!return
!     IF(max_nproj==0)RETURN
!     !
!     dots(:,:)=0.d0
!     !all type
!     DO Ity=1,naty
!        IF(psp(Ity)%nproj==0) CYCLE
!        !all atom
!        DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
!           !all projectors
!           DO Ipj=1,psp(Ity)%nproj
!              tmp0=0.d0
!              !all points
!              DO Ip=1,nlpot(Ia)%npts
!                 Id=nlpot(Ia)%Id(Ip)
!                 !<proj_lm|p>
!                 tmp0=tmp0+nlpot(Ia)%proj(Ip,Ipj)*p(Id)
!              ENDDO
!              !save dots
!              ! dots(Ipj,Ia)=tmp0*psp(Ity)%D0(Ipj,Ipj)
!              dots(1:psp(Ity)%nproj,Ia)=dots(1:psp(Ity)%nproj,Ia)  &
!                       & + tmp0*psp(Ity)%Dij(:,Ipj)
!           ENDDO
!        ENDDO
!     ENDDO
!     !scale
!     dots(:,:)=dots(:,:)*dvol
!     !Multiply the nonlocal vectors
!     DO Ity=1,naty
!        !all atom
!        DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
!           !all projectors
!           DO Ipj=1,psp(Ity)%nproj
!              tmp0=dots(Ipj,Ia)
!              !all points
!              DO Ip=1,nlpot(Ia)%npts
!                 Id=nlpot(Ia)%Id(Ip)
!                 !q=\SUM_lm{dots_lm*|proj_lm>}
!                 q(Id)=q(Id)+tmp0*nlpot(Ia)%proj(Ip,Ipj)
!              ENDDO
!           ENDDO
!        ENDDO
!     ENDDO
!     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!  ENDSUBROUTINE real_nlocmatvec
!-----------------------DIVIDER-LINE--------------------------
  SUBROUTINE real_matvec_m(mat,p,q,dimen)
     IMPLICIT NONE
     !
     INTEGER(I4B),INTENT(IN) :: dimen
     REAL(DP),INTENT(IN) :: mat(dimen,dimen),p(dimen)
     REAL(DP),INTENT(OUT) :: q(dimen)
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     q(:)=MATMUL(mat,p)
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE real_matvec_m
!-----------------------PARTING-LINE--------------------------
   !####################################################################!
   !             For real/complex Kinetic operator (Tu)                 !
   !*Author : Qiang Xu                                                  !
   !*Date   : 2017/08/28                                                !
   !####################################################################!
   SUBROUTINE cmplx_keop(uk,Ik,ts)
      !#################################################################!
      !ksKE_op\psi=Ts1+Ts2+k^2/2\psi                                    !
      !#################################################################!
      USE struct_module, ONLY : recip_lat
      USE grid_module , ONLY : kpt,nr
      USE Finite_module, ONLY : cmplx_pbc_nabla2,cmplx_pbc_nabla1
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN)  :: Ik
      COMPLEX(DCP),INTENT(IN)  :: uk(nr)
      COMPLEX(DCP),INTENT(OUT) :: ts(nr)
      !LOCAL
      REAL(DP) :: k2,mkpt(3)
      COMPLEX(DCP) :: Gur(3,nr)
      INTEGER(I4B) :: ip
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      k2=SUM(kpt%vcar(:,Ik)*kpt%vcar(:,Ik))/2
      !1)Ts1=-1/2*(\nabla)^2\psi
      CALL cmplx_pbc_nabla2(uk,ts)
      ts=-0.5_dp*ts

      IF(k2>xtiny)THEN
         !2)Ts2=-i\kvec\cdot\nabla 
         mkpt(:) = -IMAG*kpt%vcar(:,Ik)
         CALL cmplx_pbc_nabla1(uk,Gur)
         DO ip=1,nr
            ts(ip)= ts(ip) + SUM(mkpt(:)*Gur(:,ip))
         ENDDO
         !3)Ts3=k^2/2\psi
         ts(:) = ts(:) + k2*uk(:)
      ENDIF
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE cmplx_keop
!-----------------------PARTING-LINE--------------------------
ENDMODULE matvec_module
