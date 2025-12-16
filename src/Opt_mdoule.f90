MODULE opt_module
  USE constants
#ifdef MPI
  USE smpi_math_module
#endif
  IMPLICIT NONE
  REAL(DP),ALLOCATABLE :: D(:),GOLD(:),W(:)
  !For LBFGS-B
  INTEGER(I4B) :: myicount=0,maxcount=200
  INTEGER(I4B) :: nline=0,nBFGS=0
  REAL(DP)     :: dsave(29),histW(4)=0.d0
  INTEGER(I4B) :: isave(44)
  CHARACTER(LEN=60) :: task,csave
  LOGICAL :: lsave(4)
  REAL(DP),ALLOCATABLE :: wa(:)
  INTEGER(I4B),ALLOCATABLE :: iwa(:)
  !For CG
  INTEGER(I4B) :: IFLAG,Irest ,ICALL
  LOGICAL :: finish
  !For FIRE
  REAL(DP) :: TIMEmax  &
           &, TIM      &
           &,FIREalpha,fd
  REAL(DP),ALLOCATABLE :: Gd(:)
  REAL(DP),ALLOCATABLE :: Vel(:)
  INTEGER(I4B) :: NEG
  LOGICAL :: lNEG,Lfirst
  INTEGER(I4B) :: It1,It2
CONTAINS
  !---------------LBFGS interface---------------------
  !##################################################!
  !CALL setulb(n,m,x,l,u,nbd,f,g,factr,pgtol,wa,iwa &!
  !  &  task,iprint,csave,lsave,isave,dsave)         !
  !##################################################!
  SUBROUTINE initLBFGS(n,m)
     IMPLICIT NONE
     INTEGER(I4B),INTENT(IN) :: n,m
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     IF(ALLOCATED(wa))  DEALLOCATE(wa)
     IF(ALLOCATED(iwa)) DEALLOCATE(iwa)
     ALLOCATE(wa(2*m*n+5*n+11*m*m+8*m))
     ALLOCATE(iwa(3*n))
     task='LBFGS-START'
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE initLBFGS

  SUBROUTINE LBFGS_optm(Iter,n,m,x,f,g)
     IMPLICIT NONE
     !INOUT
     INTEGER(I4B),INTENT(IN) :: Iter,n,m
     REAL(DP),INTENT(INOUT) :: x(n)
              !x(n)---IN : x(k) , OUT: x(k+1)
     REAL(DP),INTENT(IN) :: f   &  !vaule f(x)
                       & ,  g(n)   !gradient g(x)
     !LOCAL
     REAL(DP),DIMENSION(n) :: l,u !lower and upper bond
     REAL(DP) :: factr=1D7    &  !
              &, pgtol=1.0D-5
     INTEGER(I4B) :: nbd(n) &
                  &, iprint=99
     !
     INTEGER(I4B) :: Ip
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     !initialize
     nbd(:)=0
     l(:)=-1D2
     u(:)= 1D2
     !PRINT*,'LBFGS-TASK:  ',task
     !search
     !CALL setulb( n,m,x,l,u,nbd,f,g,factr,pgtol,   &
     CALL setulb( n,m,x,l,u,nbd,f,g,0._DP,0._DP,   &
           &      wa,iwa,task,iprint,csave,lsave,  &
           &      isave,dsave)
     !
     x(:)=x(:)-SUM(x)/n
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE LBFGS_optm

  !-----------------CG+ interface---------------------
  !##################################################!
  !CALL cgfam(n,x,f,g,D,gold,iprint,eps,W,Iflag,&    !
  !     & irest,method,finish )                      !
  !##################################################!
  !-----------------initialize CG+--------------------
  SUBROUTINE initCG(ndim)
     IMPLICIT NONE
     INTEGER(I4B),INTENT(IN) :: ndim
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     Irest=1
     IFLAG=0
     ICALL=0
     finish=.false.
     IF(ALLOCATED(D))  DEALLOCATE(D)
     IF(ALLOCATED(GOLD))  DEALLOCATE(GOLD)
     IF(ALLOCATED(W))  DEALLOCATE(W)
     ALLOCATE(D(ndim))
     ALLOCATE(GOLD(ndim))
     ALLOCATE(W(ndim))
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE initCG
  SUBROUTINE CGplus_optm(Iter,ndim,x,f,g)
     IMPLICIT NONE
     INTEGER(I4B),INTENT(IN) :: Iter,ndim
     REAL(DP),INTENT(INOUT) :: x(ndim)
     REAL(DP),INTENT(IN) :: f   &  !vaule f(x)
                       & ,  g(ndim)   !gradient g(x)
     !LOCAL
     INTEGER(I4B) :: method=2,iprint(2),I
     REAL(DP) :: TLEV,eps=1.0D-5
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     iprint(1)=1
     iprint(2)=0
2016 CONTINUE
     CALL CGFAM(ndim,x,f,g,D,GOLD,iprint,EPS,W,IFLAG, &
              & Irest,method,finish)
     IF(IFLAG==2)THEN
        TLEV=eps*(1.d0+DABS(f))
        I=0
40      I=I+1
        IF(I>ndim)THEN
           finish=.TRUE.
           GOTO 2016
        ENDIF
        IF(DABS(g(I))>TLEV)THEN
           GOTO 2016
        ELSE
           GOTO 40
        ENDIF
      ENDIF
     x(:)=x(:)-SUM(x)/ndim
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE CGplus_optm
  !-----------------FIRE interface--------------------
  !##################################################!
  !         FIRE                                     !
  !##################################################!
  SUBROUTINE initFIRE(n,TIMES)
     IMPLICIT NONE
     INTEGER(I4B),INTENT(IN) :: n
     REAL(DP),INTENT(IN) :: TIMES
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     !Ldone=.FALSE.
     Lfirst=.TRUE.
     IF(ALLOCATED(Gd)) DEALLOCATE(Gd)
     IF(ALLOCATED(Vel)) DEALLOCATE(Vel)
     ALLOCATE(Gd(n),Vel(n))
     !
     NEG=0
     lNEG=.TRUE.
     TIMEmax=10.d0*TIMES
     TIM=TIMES
     FIREalpha=0.1d0
     Vel(:)=0.d0
     fd=509612.d0
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE initFIRE
  !
    SUBROUTINE FIRE_optm(n,X,f,G)
     IMPLICIT NONE
     !INOUT
     INTEGER(I4B),INTENT(IN) :: n
     REAL(DP),INTENT(INOUT) :: X(n)
     REAL(DP),INTENT(IN) :: G(n),f
     !LOCAL
     REAL(DP) :: nG(n)
     REAL(DP) :: Norm2g,Normg,Normv,Power
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     nG(:)=-G(:)
     Norm2g=DOT_PRODUCT(ng,ng)
     !
     !IF(SQRT(Norm2g)<wtol)THEN
     !   WRITE()
     !ENDIF
     IF(Lfirst)THEN
        Gd(:)=nG(:)
        Lfirst=.FALSE.
     ENDIF
     !---MD integrator---
     !velocity
     Vel(:)=Vel(:)+0.5d0*TIM*Gd(:)
     !move x
     X(:)=X(:)+TIM*Vel(:)
     !velocity
     Vel(:)=Vel(:)+0.5d0*TIM*nG(:)
     !store force
     Gd(:)=nG(:)
     !correct the velocity
     !F1---power
     Power=DOT_PRODUCT(nG,Vel)
     !-------------------------------------------
     !F2---correct
     Normg=SQRT(Norm2g)
     Normv=SQRT(DOT_PRODUCT(Vel,Vel))
     Vel(:)=(1.d0-FIREalpha)*Vel(:)+FIREalpha*ng(:)/Normg*Normv
     !F3-F4---
     IF(Power<0.d0)THEN
        TIM=TIM*0.5
        FIREalpha=0.1d0
        !set 0
        Vel(:)=0.d0
        !neg
        NEG=0
        !lneg
        !lNEG=.TRUE.
     ELSEIF(NEG>5)THEN
        TIM=MIN(TIM*1.1,TIMEmax)
        FIREalpha=FIREalpha*0.99
        !lNEG=.FALSE. 
     ELSE
        NEG=NEG+1
     ENDIF
     !IF(lNEG) NEG=NEG+1
     !Fix center?
     X(:)=X(:)-SUM(X)/n
     !-----------------------------------------------
#ifdef MPI
     IF (parallel%isroot) THEN
#endif
     PRINT*,'***********FIRE OPTIMIATION************'
     WRITE(*,*) 'root-mean-squre(g)',SQRT(Normg/n)
     PRINT*,'nWs changed',f-fd
     PRINT*,'present time step',TIM
     PRINT*,'present power','(',Power,')'
     PRINT*,'***************************************'
#ifdef MPI
     ENDIF
#endif
     fd=f
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE FIRE_optm

     !-----------------------------------------------
!-------------------------PARTING LINE-------------------------
ENDMODULE opt_module
