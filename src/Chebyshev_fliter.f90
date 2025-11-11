!#############################################################!
!        Chebyshev filter Method for diag(H)                  !
!*Author : Qiang Xu                                           !
!*Date   : 2017/11/28                                         !
!#############################################################!
MODULE chebyshev_module
   USE constants
   USE parameters , ONLY : CheM
   IMPLICIT NONE
CONTAINS
   !##############################################################!
   !*First :            Pseudo Subspace Method                    !
   !##############################################################!
   SUBROUTINE BuildSubspace(nps,nev,veff,eig)
      USE parameters , ONLY : nspin,CheM0,IGamma
      USE grid_module , ONLY : eigen_type,nk
      IMPLICIT NONE
      !INOUT
      INTEGER(I4B),INTENT(IN) :: nps,nev !points and states
      REAL(DP),INTENT(IN) :: veff(nps,nspin)
      TYPE(eigen_type),INTENT(INOUT) :: eig
      !LOCAL
      INTEGER(I4B) :: Is,Ik,Ii,Ip  &
              &, Ity,Ib,nlm,Ia
      REAL(DP) :: randt,randr,radius=1._DP
      INTEGER(I4B) :: Nwfa  &! total orbital-like states we have
               &, Nrand,irand1  !# of random states in subspace
      REAL(DP) :: X0(nps,nev)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      PRINT*,'Building Cheby-Subspace for this task ...'
      Nrand=nev
      irand1=1
      !random states
      IF(Nrand>0)THEN
         !initialize the temp subspace
         CALL random_seed()
         DO Ii=irand1,nev
            DO Ip=1,nps
               CALL random_number(randt)
               !CALL random_number(randr)
               randt = -radius+randt*2.d0*radius
               !randr = -radius+randr*10.d0*radius
               X0(Ip,Ii) = randt !(randt+randr)/2
            ENDDO
         ENDDO

      ENDIF
      !Raleigh-Ritz step
      !Gamma point
      DO Is=1,Nspin !spiner

         DO Ik=1,nk !k-points
            !print*,'Ik',Ik

            IF(Ik/=IGamma)THEN
            !for non-Gamma k-points
               eig%wvf(:,:,Ik,Is)=X0(:,:)
               IF(CheM0>0)THEN
                   !filter
                   CALL cmplx_first_filter(nps,nev,Ik,veff(:,Is), &
                     & eig%wvf(:,:,Ik,Is),eig%val(:,Ik,Is))
               ELSE
                   !RR
                   CALL cmplx_first_RRstep(nps,nev,Ik,veff(:,Is), &
                      &   eig%wvf(:,:,Ik,Is),eig%val(:,Ik,Is))
               ENDIF
            ELSE
            !for Gamma k-points
               eig%wvfG(:,:,Is)=X0(:,:)
               IF(CheM0>0)THEN
                   !filter
                   CALL real_first_filter(nps,nev,veff(:,Is), &
                     & eig%wvfG(:,:,Is),eig%val(:,Ik,Is))
               ELSE
                   !RR
                   CALL real_first_RRstep(nps,nev,veff(:,Is), &
                      &   eig%wvfG(:,:,Is),eig%val(:,Ik,Is))
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE BuildSubspace
   !----------------------------------------------------------
   SUBROUTINE CheF_all(nps,nev,veff,eig)
      !need to be improve for spin
      USE parameters , ONLY : nspin,IGamma
      USE grid_module , ONLY : nk,eigen_type
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN) :: nps,nev
      REAL(DP),INTENT(IN) :: veff(nps,nspin)
      TYPE(eigen_type),INTENT(INOUT) :: eig
      !LOCAL
      INTEGER(I4B) :: Is,Ik
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      DO Is=1,nspin
      DO Ik=1,nk
         IF(Ik==IGamma)THEN
            CALL real_filtering(nps,nev,veff(:,Is),eig%wvfG(:,:,Is),eig%val(:,Ik,Is))
         ELSE
            CALL cmplx_filtering(nps,nev,Ik,veff(:,Is),eig%wvf(:,:,Ik,Is),eig%val(:,Ik,Is))
         ENDIF
      ENDDO
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE CheF_all
   !----------------------------------------------------------
   SUBROUTINE real_first_RRstep(nps,nev,veff,X,D)
      USE Lapack_module , ONLY : lapk_OrthNorm, lapk_Eigen, lapk_MM
      USE parameters , ONLY : CheM0,LRROrthNorm
      IMPLICIT NONE
      !IN/OUT
      INTEGER(I4B),INTENT(IN) :: nps,nev
      REAL(DP),INTENT(IN) :: veff(nps)
      REAL(DP),INTENT(INOUT) :: X(nps,nev) !subspace
      REAL(DP),INTENT(OUT) :: D(:)      !eigenvalue
      !LOCAL
      REAL(DP) :: Xnew(nps,nev)
      !REAL(DP) :: a,b,al
      REAL(DP),DIMENSION(nev,nev) :: Shat,Hhat,Qs
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !Raleigh-Ritz step
      IF(LRROrthNorm)THEN
         !OrthNorm
         CALL lapk_OrthNorm(X)
         !RR
         CALL real_Rayleigh_quotient(nps,nev,veff,X,Hhat)
         !eigen-decomposion
         CALL lapk_Eigen(Hhat,Qs,D)
      ELSE
         !Overlap matrix
         CALL lapk_MM(X,X,'T','N',1.0_dp,0.0_dp,Shat)
         !projected hamiltonian
         CALL real_Rayleigh_quotient(nps,nev,veff,X,Hhat)
         !eigen-decomposion
         CALL lapk_Eigen(nev,Hhat,Shat,Qs,D)
      ENDIF
      !-------------------
      !rotation
      CALL lapk_MM(X,Qs,'N','N',1.0_dp,0.0_dp,Xnew)
      !eigen-value
      X(:,:)=Xnew(:,:)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE real_first_RRstep
   !----------------------------------------------------------
   SUBROUTINE real_init_uplow(nps,k,veff,v,a,b,al)
      !
      USE matvec_module , ONLY : matvec
      USE Lapack_module , ONLY : lapk_Eigen
      IMPLICIT NONE
      !INOUT
      INTEGER(I4B),INTENT(IN) :: nps,k
      !hamiltonian
      REAL(DP),INTENT(IN) :: veff(nps)
      REAL(DP),INTENT(INOUT) :: v(nps)
      REAL(DP),INTENT(OUT) :: a,b,al
      !LOCAL
      REAL(DP) :: v0(nps),f(nps),alpha
      REAL(DP) :: beta, &
              &   fbeta=0.5d0
      REAL(DP) :: T(k,k),evec(k,k),eval(k)
      INTEGER(I4B) :: J,Nj
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      T(:,:)=0.d0
      Nj=MIN(k,10)
      !
      CALL matvec(veff,v,f,nps)
      alpha=DOT_PRODUCT(f,v)
      f=f-alpha*v
      T(1,1)=alpha
      DO J=2,Nj
         !beta=SQRT(DOT_PRODUCT(f,f))
         beta=REAL(DOT_PRODUCT(f,f),DP)
         beta=SQRT(beta)
         v0=v
         v=f/beta
         CALL matvec(veff,v,f,nps)
         f=f-beta*v0
         alpha=DOT_PRODUCT(f,v)
         f=f-alpha*v
         T(J,J-1)=beta
         T(J-1,J)=beta
         T(J,J)=alpha
      ENDDO
      !Rayleigh-Ritz value
      !CALL diagM_real(T,evec,eval)
      !print*,'Tmat',T(:,1)
      !STOP

      CALL lapk_Eigen(T,evec,eval)

      beta=REAL(DOT_PRODUCT(f,f),8)

      a=fbeta*eval(1)+(1.d0-fbeta)*eval(k)
      al=eval(1)
      b=eval(k) + SQRT( beta)*ABS(evec(k,k)) !+1e-10
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE real_init_uplow
   !----------------------------------------------------------
  SUBROUTINE real_first_filter(nps,nst,veff,X,eval)
      USE parameters , ONLY : CF0=>CheM0,LRROrthNorm
      USE Lapack_module , ONLY : lapk_Eigen,lapk_OrthNorm,lapk_MM
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN) :: nps,nst !number of points/states
      REAL(DP),INTENT(IN) :: veff(nps)
      REAL(DP),INTENT(INOUT) :: X(:,:)
      REAL(DP),INTENT(OUT) :: eval(:)
      !LOCAL
      REAL(DP) :: a,b,al,t
      INTEGER(I4B) :: I,Niter=4
      REAL(DP) :: deval,TOL=1e-8
      REAL(DP) :: evald(nst)
      REAL(DP),DIMENSION(nst,nst)  :: Hhat,Qs,Shat
      REAL(DP)  :: Xnew(nps,nst),vec(nps)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !low up bound
      vec(:)=X(:,nst)
      CALL real_init_uplow(nps,7,veff,vec,a,b,al)
      !print*,'a,b,al',a,b,al
      !
      evald(:)=612509.d0
      DO I=1,Niter
         !
         CALL real_chebyshev_filter_scaled(nps,nst,veff,X,CF0,a,b,al)
         IF(LRROrthNorm)THEN
            !OrthNorm
            CALL lapk_OrthNorm(X)
            !xHx
            CALL real_Rayleigh_quotient(nps,nst,veff,X,Hhat)
            !eigen-decomposion
            CALL lapk_Eigen(Hhat,Qs,eval)
         ELSE
            CALL lapk_MM(X,X,'T','N',1._dp,0._dp,Shat)
            !projected hamiltonian
            CALL real_Rayleigh_quotient(nps,nst,veff,X,Hhat)
            !eigen-decomposion
            CALL lapk_Eigen(nst,Hhat,Shat,Qs,eval)
         ENDIF
         deval=SUM(ABS(eval-evald))/nst
         !-----------------
         IF(deval<TOL) EXIT
         !-----------------
         !store old eigenvalue
         evald(:)=eval(:)
         !update the new bound
         a=eval(nst)
         al=eval(1)
      ENDDO

      !rotation
      CALL lapk_MM(X,Qs,'N','N',1._dp,0._dp,Xnew)
      X=Xnew
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE real_first_filter
   !----------------------------------------------------------
   SUBROUTINE cmplx_first_RRstep(nps,nev,Ik,veff,X,D)
      USE Lapack_module , ONLY : lapk_OrthNorm, lapk_Eigen, lapk_MM
      USE parameters , ONLY : CheM0,LRROrthNorm
      IMPLICIT NONE
      !IN/OUT
      INTEGER(I4B),INTENT(IN) :: nps,nev,Ik
      REAL(DP),INTENT(IN) :: veff(nps)
      COMPLEX(DCP),INTENT(INOUT) :: X(nps,nev) !subspace
      REAL(DP),INTENT(OUT) :: D(:)      !eigenvalue
      !LOCAL
      COMPLEX(DP) :: Xnew(nps,nev)
      COMPLEX(DP),DIMENSION(nev,nev) :: Shat,Hhat,Qs
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !Raleigh-Ritz step
      IF(LRROrthNorm)THEN
         !OrthNorm
         CALL lapk_OrthNorm(X)
         !RR
         CALL cmplx_Rayleigh_quotient(nps,nev,Ik,veff,X,Hhat)
         !eigen-decomposion
         CALL lapk_Eigen(Hhat,Qs,D)
      ELSE
         !Overlap matrix
         CALL lapk_MM(X,X,'C','N',cmplx(1._dp,0._dp,DCP),cmplx(0._dp,0._dp,DCP),Shat)
         !projected hamiltonian
         CALL cmplx_Rayleigh_quotient(nps,nev,Ik,veff,X,Hhat)
         !eigen-decomposion
         CALL lapk_Eigen(nev,Hhat,Shat,Qs,D)
      ENDIF
      !-------------------
      !rotation
      CALL lapk_MM(X,Qs,'N','N',cmplx(1._dp,0._dp,DCP),cmplx(0._dp,0._dp,DCP),Xnew)
      !eigen-value
      X(:,:)=Xnew(:,:)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE cmplx_first_RRstep
   !----------------------------------------------------------
   SUBROUTINE cmplx_init_uplow(nps,k,Ik,veff,v,a,b,al)
      !
      USE matvec_module , ONLY : matvec
      USE Lapack_module , ONLY : lapk_Eigen
      IMPLICIT NONE
      !INOUT
      INTEGER(I4B),INTENT(IN) :: nps,k,Ik
      !hamiltonian
      REAL(DP),INTENT(IN) :: veff(nps)
      COMPLEX(DCP),INTENT(INOUT) :: v(nps)
      REAL(DP),INTENT(OUT) :: a,b,al
      !LOCAL
      COMPLEX(DCP) :: v0(nps),f(nps),alpha
      REAL(DP) :: beta, &
              &   fbeta=0.5d0
      COMPLEX(DCP) :: T(k,k),evec(k,k)
      REAL(DP) :: eval(k)
      INTEGER(I4B) :: J,Nj
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      T(:,:)=0.d0
      Nj=MIN(k,10)
      !
      CALL matvec(Ik,veff,v,f,nps)
      alpha=DOT_PRODUCT(f,v)
      f=f-alpha*v
      T(1,1)=alpha
      DO J=2,Nj
         beta=REAL(DOT_PRODUCT(f,f),DP)
         beta=SQRT(beta)
         v0=v
         v=f/beta
         CALL matvec(Ik,veff,v,f,nps)
         f=f-beta*v0
         alpha=DOT_PRODUCT(f,v)
         f=f-alpha*v
         T(J,J-1)=beta
         T(J-1,J)=beta
         T(J,J)=alpha
      ENDDO
      !Rayleigh-Ritz value
      !CALL diagM_real(T,evec,eval)

      !PRINT*,'T',SIZE(T),SUM(abs(T))
      !PRINT*,'V',SIZE(evec)
      !PRINT*,'E',SIZE(eval)
      CALL lapk_Eigen(T,evec,eval)

      beta=REAL(DOT_PRODUCT(f,f),DP)

      a=fbeta*eval(1)+(1.d0-fbeta)*eval(k)
      al=eval(1)
      b=eval(k) + SQRT( beta)*ABS(evec(k,k)) !+1e-10
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE cmplx_init_uplow
   !----------------------------------------------------------
  SUBROUTINE cmplx_first_filter(nps,nst,Ik,veff,X,eval)
      USE parameters , ONLY : CF0=>CheM0,LRROrthNorm
      USE Lapack_module , ONLY : lapk_Eigen,lapk_OrthNorm,lapk_MM
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN) :: nps,nst,Ik !number of points/states
      REAL(DP),INTENT(IN) :: veff(nps)
      COMPLEX(DCP),INTENT(INOUT) :: X(:,:)
      REAL(DP),INTENT(OUT) :: eval(:)
      !LOCAL
      REAL(DP) :: a,b,al,t
      INTEGER(I4B) :: I,Niter=4
      REAL(DP) :: deval,TOL=1e-8
      REAL(DP) :: evald(nst)
      COMPLEX(DCP),DIMENSION(nst,nst)  :: Hhat,Qs,Shat
      COMPLEX(DCP)  :: Xnew(nps,nst),vec(nps)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !low up bound
      vec(:)=X(:,nst)
      CALL cmplx_init_uplow(nps,7,Ik,veff,vec,a,b,al)
      !
      evald(:)=612509.d0
      DO I=1,Niter
         !
         CALL cmplx_chebyshev_filter_scaled(nps,nst,Ik,veff,X,CF0,a,b,al)
         IF(LRROrthNorm)THEN
            !OrthNorm
            CALL lapk_OrthNorm(X)
            !xHx
            CALL cmplx_Rayleigh_quotient(nps,nst,Ik,veff,X,Hhat)
            !eigen-decomposion
            CALL lapk_Eigen(Hhat,Qs,eval)
         ELSE
            CALL lapk_MM(X,X,'T','N',cmplx(1._dp,0._dp,DCP),cmplx(0._dp,0._dp,DCP),Shat)
            !projected hamiltonian
            CALL cmplx_Rayleigh_quotient(nps,nst,Ik,veff,X,Hhat)
            !eigen-decomposion
            CALL lapk_Eigen(nst,Hhat,Shat,Qs,eval)
         ENDIF
         deval=SUM(ABS(eval-evald))/nst
         !-----------------
         IF(deval<TOL) EXIT
         !-----------------
         !store old eigenvalue
         evald(:)=eval(:)
         !update the new bound
         a=eval(nst)
         al=eval(1)
      ENDDO

      !rotation
      CALL lapk_MM(X,Qs,'N','N',cmplx(1._dp,0._dp,DCP),cmplx(0._dp,0._dp,DCP),Xnew)
      X=Xnew
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE cmplx_first_filter
   !##############################################################!
   !*For    :     real Filtering (Standard RR)                    !
   !*Author : Qiang Xu                                            !
   !*Date   : 2018-03-05                                          !
   !##############################################################!
   !---------------------Rayleigh-quotient-------------------------
   SUBROUTINE real_Rayleigh_quotient(nps,nst,veff,x,xhx)
      !xhx=(X,H,X)
      USE Lapack_module , ONLY : lapk_MM
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN) :: nps,nst
      REAL(DP),INTENT(IN) :: veff(nps),x(nps,nst)
      REAL(DP),INTENT(OUT) :: xhx(nst,nst)
      !LOCAL
      REAL(DP) :: hx(nps,nst)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CALL real_HX(nps,nst,veff,x,hx)
      CALL lapk_MM(x,hx,'T','N',1._dp,0._dp,xhx)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE real_Rayleigh_quotient
   !-------------------------------HX------------------------------
   SUBROUTINE real_HX(nps,nst,veff,V,HV)
      USE matvec_module , ONLY : matvec
      IMPLICIT NONE
      !INOUT
      INTEGER(I4B),INTENT(IN) :: nps,nst
      REAL(DP),INTENT(IN)  :: veff(nps) !effective potential
      REAL(DP),INTENT(IN)  :: V(nps,nst)
      REAL(DP),INTENT(OUT) :: HV(nps,nst)
      !LOCAL
      INTEGER(I4B) :: Is
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      DO Is=1,nst,1
         CALL matvec(veff,V(:,Is),HV(:,Is),nps)
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE real_HX
   !-------------------------upper bound---------------------------
   SUBROUTINE real_Estupb(nps,k,veff,vec,b)
      !
      USE matvec_module , ONLY : matvec
      USE Lapack_module , ONLY : lapk_Eigen
      IMPLICIT NONE
      !INOUT
      INTEGER(I4B),INTENT(IN) :: nps,k
      !hamiltonian
      REAL(DP),INTENT(IN) :: veff(nps)
      REAL(DP),INTENT(IN) :: vec(nps) 
      REAL(DP),INTENT(OUT) :: b
      !LOCAL
      REAL(DP) :: v0(nps),v(nps),f(nps),alpha
      REAL(DP) :: beta,eval(k),mz
      REAL(DP) :: T(k,k),evec(k,k)
      INTEGER(I4B) :: J,Nj
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      v=vec !should normalize
      T(:,:)=0.d0
      Nj=MIN(k,10)
      !
      CALL matvec(veff,v,f,nps)
      alpha=DOT_PRODUCT(f,v)
      f=f-alpha*v
      T(1,1)=alpha
      DO J=2,Nj
         beta=REAL(DOT_PRODUCT(f,f),DP)
         beta=SQRT(beta)
         v0=v
         v=f/beta
         CALL matvec(veff,v,f,nps)
         f=f-beta*v0
         alpha=DOT_PRODUCT(f,v)
         f=f-alpha*v
         T(J,J-1)=beta
         T(J-1,J)=beta
         T(J,J)=alpha
      ENDDO
      !NORM2(T)
      !b=Norm_2(T,k) + SQRT(REAL(DOT_PRODUCT(f,f),8))
      CALL lapk_Eigen(T,evec,eval)
      beta=REAL(DOT_PRODUCT(f,f),DP)
      !mz=MAX( ABS(evec(k,k)), ABS(evec(k,k-1)) , ABS(evec(k,k-2)) )
      mz=ABS(evec(k,k))
      b=eval(k) + SQRT(beta)*mz
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE real_Estupb
   !-------------------chebyshev_filter---------------------
   SUBROUTINE real_chebyshev_filter(nps,nst,veff,X,m,a,b)
      !
      IMPLICIT NONE
      !INOUT
      INTEGER(I4B),INTENT(IN) :: nps,nst !number of states
      REAL(DP),INTENT(IN) :: veff(nps) !veff
      REAL(DP),INTENT(INOUT) :: X(nps,nst)
      INTEGER(I4B),INTENT(IN)   :: m ! the m degree Chebyshev polynomial we used
      REAL(DP),INTENT(IN) :: a,b  !interval [a,b] to be filter
      !LOCAL        
      REAL(DP) :: e &   !(b-a)/2
            & ,   c    !(b+a)/2
      !
      REAL(DP),DIMENSION(nps,nst)  :: HV,Y,Ynew
      INTEGER(I4B) :: Ic
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF(m<2)THEN
         WRITE(*,*) 'Chebyshev_filter: the degree m must larger than 1'
      ENDIF 
      
      e = ( b - a ) / 2
      c = ( b + a ) / 2
      
      CALL real_HX(nps,nst,veff,X,HV)
      Y(:,:)=( HV(:,:) - c*X(:,:) ) / e

      !Chebyshev filtering
      DO Ic=2,m
      
         !CALL HY
         CALL real_HX(nps,nst,veff,Y,HV)
         Ynew=2.d0*( HV-c*Y ) / e - X
         !store the Y

         X=Y

         Y=Ynew


      ENDDO
      !out put filted space
      X=Ynew
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE real_chebyshev_filter
   !-------------------chebyshev_filter_scaled---------------------
   SUBROUTINE real_chebyshev_filter_scaled(nps,nst,veff,X,m,a,b,al)
      !
      IMPLICIT NONE
      !INOUT
      INTEGER(I4B),INTENT(IN) :: nps,nst !number of states
      REAL(DP),INTENT(IN) :: veff(nps) !veff
      REAL(DP),INTENT(INOUT) :: X(nps,nst)
      INTEGER(I4B),INTENT(IN)   :: m ! the m degree Chebyshev polynomial we used
      REAL(DP),INTENT(IN) :: a,b,al  !interval [a,b] to be filter
      !LOCAL        
      REAL(DP) :: e &   !(b-a)/2
            & ,   c &   !(b+a)/2
            & , sigma,sigmanew,tau
      !
      REAL(DP),DIMENSION(nps,nst)  :: HV,Y,Ynew
      INTEGER(I4B) :: Ic
      REAL(DP) :: temp1,temp2
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF(m<2)THEN
         WRITE(*,*) 'Chebyshev_filter: the degree m must larger than 1'
      ENDIF 
      
      e = ( b - a ) / 2
      c = ( b + a ) / 2
      sigma=e / (c-al)
      tau=2.d0 / sigma
      CALL real_HX(nps,nst,veff,X,HV)
      temp1=sigma / e
      Y= ( HV - c*X ) * temp1
      DO Ic=2,m
         sigmanew=1.d0 / (tau-sigma)
         CALL real_HX(nps,nst,veff,Y,HV)
         temp1=2.d0*sigmanew/e
         temp2=sigma*sigmanew
         Ynew=( HV - c*Y )*temp1 - temp2*X
         X=Y
         Y=Ynew
         sigma=sigmanew
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE real_chebyshev_filter_scaled
   !-----------------------GRayleigh_Ritz--------------------------
   SUBROUTINE real_GRayleigh_Ritz(nps,nev,veff,X,D)
      USE Lapack_module , ONLY : lapk_Eigen,lapk_MM
      IMPLICIT NONE
      !IN/OUT
      INTEGER(I4B),INTENT(IN) :: nps,nev
      REAL(DP),INTENT(IN) :: veff(nps)
      REAL(DP),INTENT(INOUT) :: X(nps,nev)
      REAL(DP),INTENT(OUT) :: D(:)
      !
      REAL(DP),DIMENSION(nev,nev) :: S_hat,H_hat,Q
      REAL(DP) :: Xnew(nps,nev)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !Calculate the overlap matrix
      CALL lapk_MM(X,X,'T','N',1._dp,0._dp,S_hat)
      !Calculate the project hamiltion
      CALL real_Rayleigh_quotient(nps,nev,veff,X,H_hat)
      !solve the generalized eigenvalue problem
      CALL lapk_Eigen(nev,H_hat,S_hat,Q,D)
      !X=XQ
      CALL lapk_MM(X,Q,'N','N',1._dp,0._dp,Xnew)
      X=Xnew
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE real_GRayleigh_Ritz
   !---------------------Rayleigh-Ritz step------------------------
   SUBROUTINE real_Rayleigh_Ritz(nps,sn,veff,X,D)
      USE Lapack_module ,  ONLY : lapk_OrthNorm,lapk_Eigen,lapk_MM
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN) :: nps,sn
      REAL(DP),INTENT(IN) :: veff(nps)
      REAL(DP),INTENT(INOUT) :: X(nps,sn) 
      REAL(DP),INTENT(OUT) :: D(:)
      !LOCAL
      REAL(DP) :: Xnew(nps,sn)
      REAL(DP),DIMENSION(sn,sn) :: Hhat,Qs
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !F0:OrthNorm
      CALL lapk_OrthNorm(X)
      !F1:Rayleigh_quotient
      CALL real_Rayleigh_quotient(nps,sn,veff,X,Hhat)
      !F2:eigen-decomposition Q,D
      CALL lapk_Eigen(Hhat,Qs,D)
      !F3:X=XQ
      !X=MATMUL( X , Q )
      CALL lapk_MM(X,Qs,'N','N',1.0_dp,0.0_dp,Xnew)
      X=Xnew
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE real_Rayleigh_Ritz
   !-----------------non-OrthNorm Chebyshev_filter ----------------
   SUBROUTINE real_filtering(nps,nev,veff,X,D)
      USE parameters , ONLY : LRROrthNorm
      !To avoid OrthNorm
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN) :: nps,nev
      REAL(DP),INTENT(IN) :: veff(nps)
      REAL(DP),INTENT(INOUT) :: X(nps,nev)
      REAL(DP),INTENT(INOUT) :: D(:)  !rayleigh-ritz value
      !LOCAL
      REAL(DP) :: a,b,al
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      a=MAXVAL(D)
      al=MINVAL(D)
      !up boundary
      CALL real_Estupb(nps,7,veff,X(:,nev),b)
      !filtering (a,b) 
      CALL real_chebyshev_filter_scaled(nps,nev,veff,X,CheM,a,b,al)
      !CALL chebyshev_filter_real(nps,nev,veff,X,CheM,a,b)
      IF(LRROrthNorm)THEN
         !RR (Rayleigh-Ritz Step)
         CALL real_Rayleigh_Ritz(nps,nev,veff,X,D)
      ELSE
         !GRR (Rayleigh-Ritz Step)
         CALL real_GRayleigh_Ritz(nps,nev,veff,X,D)
      ENDIF
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE real_filtering
   !##############################################################!
   !*For    :     cmplx Filtering (Standard RR)                   !
   !*Author : Qiang Xu                                            !
   !*Date   : 2018-03-05                                          !
   !##############################################################!
   !---------------------Rayleigh-quotient-------------------------
   SUBROUTINE cmplx_Rayleigh_quotient(nps,nst,Ik,veff,x,xhx)
      !xhx=(X,H,X)
      USE Lapack_module , ONLY : lapk_MM
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN) :: nps,nst,Ik
      REAL(DP),INTENT(IN) :: veff(nps)
      COMPLEX(DCP),INTENT(IN) :: x(nps,nst)
      COMPLEX(DCP),INTENT(OUT) :: xhx(nst,nst)
      !LOCAL
      COMPLEX(DCP) :: hx(nps,nst)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CALL cmplx_HX(nps,nst,Ik,veff,x,hx)
      CALL lapk_MM(x,hx,'C','N',cmplx(1.0,0.0,DCP),cmplx(0._dp,0._dp,DCP),xhx)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE cmplx_Rayleigh_quotient
   !-------------------------------HX------------------------------
   SUBROUTINE cmplx_HX(nps,nst,Ik,veff,V,HV)
      USE matvec_module , ONLY : matvec
      IMPLICIT NONE
      !INOUT
      INTEGER(I4B),INTENT(IN) :: nps,nst,Ik
      REAL(DP),INTENT(IN)  :: veff(nps) !effective potential
      COMPLEX(DCP),INTENT(IN)  :: V(nps,nst)
      COMPLEX(DCP),INTENT(OUT) :: HV(nps,nst)
      !LOCAL
      INTEGER(I4B) :: Is
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      DO Is=1,nst,1
         CALL matvec(Ik,veff(:),V(:,Is),HV(:,Is),nps)
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE cmplx_HX
   !-------------------------upper bound---------------------------
   SUBROUTINE cmplx_Estupb(nps,k,Ik,veff,vec,b)
      !
      USE matvec_module , ONLY : matvec
      USE Lapack_module , ONLY : lapk_Eigen
      IMPLICIT NONE
      !INOUT
      INTEGER(I4B),INTENT(IN) :: nps,k,Ik
      !hamiltonian
      REAL(DP),INTENT(IN) :: veff(nps)
      COMPLEX(DCP),INTENT(IN) :: vec(nps) 
      REAL(DP),INTENT(OUT) :: b
      !LOCAL
      COMPLEX(DCP) :: v0(nps),v(nps),f(nps),alpha
      REAL(DP) :: beta,eval(k),mz
      COMPLEX(DCP) :: T(k,k),evec(k,k)
      INTEGER(I4B) :: J,Nj
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      v=vec !should normalize
      T(:,:)=0.d0
      Nj=MIN(k,10)
      !
      CALL matvec(Ik,veff,v,f,nps)
      alpha=DOT_PRODUCT(f,v)
      f=f-alpha*v
      T(1,1)=alpha
      DO J=2,Nj
         beta=REAL(DOT_PRODUCT(f,f),DP)
         beta=SQRT(beta)
         v0=v
         v=f/beta
         CALL matvec(Ik,veff,v,f,nps)
         f=f-beta*v0
         alpha=DOT_PRODUCT(f,v)
         f=f-alpha*v
         T(J,J-1)=beta
         T(J-1,J)=beta
         T(J,J)=alpha
      ENDDO
      !NORM2(T)
      !b=Norm_2(T,k) + SQRT(REAL(DOT_PRODUCT(f,f),8))
      CALL lapk_Eigen(T,evec,eval)
      beta=REAL(DOT_PRODUCT(f,f),DP)
      mz=ABS(evec(k,k))
      b=eval(k) + SQRT(beta)*mz
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE cmplx_Estupb
   !-------------------chebyshev_filter_scaled---------------------
   SUBROUTINE cmplx_chebyshev_filter_scaled(nps,nst,Ik,veff,X,m,a,b,al)
      !
      IMPLICIT NONE
      !INOUT
      INTEGER(I4B),INTENT(IN) :: nps,nst,Ik !number of states
      REAL(DP),INTENT(IN) :: veff(nps) !veff
      COMPLEX(DCP),INTENT(INOUT) :: X(nps,nst)
      INTEGER(I4B),INTENT(IN)   :: m ! the m degree Chebyshev polynomial we used
      REAL(DP),INTENT(IN) :: a,b,al  !interval [a,b] to be filter
      !LOCAL        
      REAL(DP) :: e &   !(b-a)/2
            & ,   c &   !(b+a)/2
            & , sigma,sigmanew,tau
      !
      COMPLEX(DP),DIMENSION(nps,nst)  :: HV,Y,Ynew
      INTEGER(I4B) :: Ic
      REAL(DP) :: temp1,temp2
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF(m<2)THEN
         WRITE(*,*) 'Chebyshev_filter: the degree m must larger than 1'
      ENDIF 
      
      e = ( b - a ) / 2
      c = ( b + a ) / 2
      sigma=e / (c-al)
      tau=2.d0 / sigma
      CALL cmplx_HX(nps,nst,Ik,veff,X,HV)
      temp1=sigma / e
      Y= ( HV - c*X ) * temp1
      DO Ic=2,m
         sigmanew=1.d0 / (tau-sigma)
         CALL cmplx_HX(nps,nst,Ik,veff,Y,HV)
         temp1=2.d0*sigmanew/e
         temp2=sigma*sigmanew
         Ynew=( HV - c*Y )*temp1 - temp2*X
         X=Y
         Y=Ynew
         sigma=sigmanew
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE cmplx_chebyshev_filter_scaled
   !-----------------------GRayleigh_Ritz--------------------------
   SUBROUTINE cmplx_GRayleigh_Ritz(nps,nev,Ik,veff,X,D)
      USE Lapack_module , ONLY : lapk_Eigen,lapk_MM
      IMPLICIT NONE
      !IN/OUT
      INTEGER(I4B),INTENT(IN) :: nps,nev,Ik
      REAL(DP),INTENT(IN) :: veff(nps)
      COMPLEX(DCP),INTENT(INOUT) :: X(nps,nev)
      REAL(DP),INTENT(OUT) :: D(:)
      !
      COMPLEX(DCP),DIMENSION(nev,nev) :: S_hat,H_hat,Q
      COMPLEX(DCP) :: Xnew(nps,nev)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !Calculate the overlap matrix
      CALL lapk_MM(X,X,'T','N',cmplx(1._dp,0._dp,DCP),cmplx(0._dp,0._dp,DCP),S_hat)
      !Calculate the project hamiltion
      CALL cmplx_Rayleigh_quotient(nps,nev,Ik,veff,X,H_hat)
      !solve the generalized eigenvalue problem
      CALL lapk_Eigen(nev,H_hat,S_hat,Q,D)
      !X=XQ
      CALL lapk_MM(X,Q,'N','N',cmplx(1._dp,0._dp,DCP),cmplx(0._dp,0._dp,DCP),Xnew)
      X=Xnew
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE cmplx_GRayleigh_Ritz
   !---------------------Rayleigh-Ritz step------------------------
   SUBROUTINE cmplx_Rayleigh_Ritz(nps,sn,Ik,veff,X,D)
      USE Lapack_module ,  ONLY : lapk_OrthNorm,lapk_Eigen,lapk_MM
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN) :: nps,sn,Ik
      REAL(DP),INTENT(IN) :: veff(nps)
      COMPLEX(DCP),INTENT(INOUT) :: X(nps,sn) 
      REAL(DP),INTENT(OUT) :: D(:)
      !LOCAL
      COMPLEX(DCP) :: Xnew(nps,sn)
      COMPLEX(DCP),DIMENSION(sn,sn) :: Hhat,Qs
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !F0:OrthNorm
      CALL lapk_OrthNorm(X)
      !F1:Rayleigh_quotient
      CALL cmplx_Rayleigh_quotient(nps,sn,Ik,veff,X,Hhat)
      !F2:eigen-decomposition Q,D
      CALL lapk_Eigen(Hhat,Qs,D)
      !F3:X=XQ
      !X=MATMUL( X , Q )
      CALL lapk_MM(X,Qs,'N','N',cmplx(1._dp,0._dp,DCP),cmplx(0._dp,0._dp,DCP),Xnew)
      X=Xnew
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE cmplx_Rayleigh_Ritz
   !-----------------non-OrthNorm Chebyshev_filter ----------------
   SUBROUTINE cmplx_filtering(nps,nev,Ik,veff,X,D)
      USE parameters , ONLY : LRROrthNorm
      !To avoid OrthNorm
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN) :: nps,nev,Ik
      REAL(DP),INTENT(IN) :: veff(nps)
      COMPLEX(DCP),INTENT(INOUT) :: X(nps,nev)
      REAL(DP),INTENT(INOUT) :: D(:)  !rayleigh-ritz value
      !LOCAL
      REAL(DP) :: a,b,al
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      a=MAXVAL(D)
      al=MINVAL(D)
      !up boundary
      CALL cmplx_Estupb(nps,7,Ik,veff,X(:,nev),b)
      !filtering (a,b) 
      CALL cmplx_chebyshev_filter_scaled(nps,nev,Ik,veff,X,CheM,a,b,al)
      !CALL chebyshev_filter_real(nps,nev,veff,X,CheM,a,b)
      IF(LRROrthNorm)THEN
         !RR (Rayleigh-Ritz Step)
         CALL cmplx_Rayleigh_Ritz(nps,nev,Ik,veff,X,D)
      ELSE
         !GRR (Rayleigh-Ritz Step)
         CALL cmplx_GRayleigh_Ritz(nps,nev,Ik,veff,X,D)
      ENDIF
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE cmplx_filtering
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDMODULE chebyshev_module
