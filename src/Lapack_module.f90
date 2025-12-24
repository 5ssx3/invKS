MODULE Lapack_module
!###########################################################!
!*For: Lapack                                               !
!*Author:Qiang Xu                                           !
!*Date:2017-7-20  & 2023-11-01 resort                       !
!###########################################################!
   USE constants
   !USE MKL95_PRECISION
   IMPLICIT NONE
   INTEGER(I4B) :: MMAX=1000000 !max dim of matrix
   !EigenSolver
   interface lapk_Eigen!{{{
      module procedure lapk_real_SEigen
      module procedure lapk_real_GEigen
      module procedure lapk_cmplx_SEigen
      module procedure lapk_cmplx_GEigen
   end interface lapk_Eigen!}}}
!OrthNorm
   interface lapk_OrthNorm!{{{
      module procedure lapk_real_OrthNorm
      module procedure lapk_cmplx_OrthNorm
   end interface lapk_OrthNorm!}}}
!MatMat
   interface lapk_MM!{{{
      module procedure lapk_real_matmat
      module procedure lapk_cmplx_matmat
   end interface lapk_MM!}}}
!MatVex
   interface lapk_MV!{{{
      module procedure lapk_real_matvec
      module procedure lapk_cmplx_matvec
   end interface lapk_MV!}}}
!invM
   interface lapk_invM!{{{
      module procedure lapk_real_invmat
      module procedure lapk_cmplx_invmat
   end interface lapk_invM!}}}
!Norm2
   interface lapk_Norm2!{{{
      module procedure lapk_cmplx_Norm2
   end interface lapk_Norm2!}}}
CONTAINS
   !----------------------------------------------------------
   !----------------------------------------------------------
  !#########################################################!
  !  CALL ZHEEV(jobs,uplo,n,a,lda,w,work,lwork,rwork,info)  !
  !#########################################################!
  SUBROUTINE lapk_cmplx_SEigen(mat,evec,eval)
     !
     IMPLICIT NONE
     COMPLEX(DP),INTENT(IN)  :: mat(:,:) !IN:martix
     COMPLEX(DP),INTENT(OUT) :: evec(:,:) !OUT:engivector
     REAL(DP),INTENT(OUT) :: eval(:)     !engivalue
     !
     INTEGER(I4B) :: lda            !see Lapack ZHEEV
     INTEGER(I4B) :: lwork !,lwmax=1000
     COMPLEX(DCP) :: work_tmp(10)
     COMPLEX(DCP),ALLOCATABLE :: work(:)
     REAL(DP),ALLOCATABLE     :: rwork(:)
     INTEGER(I4B) :: info
     INTEGER(I4B) :: dime
     !
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     dime=SIZE(mat,1)
     IF(dime>MMAX)THEN
        WRITE(6,*) 'diagM:Check the dimension of MMAX'
        STOP
     ENDIF
     lda=dime
     !PRINT*,SIZE(mat),SIZE(evec),SIZE(eval)
     !---creat arrays
     ALLOCATE(rwork(3*dime-2))
     !PRINT*,SIZE(rwork)
     !use mat for lapack
     evec(:,:)=mat(:,:)
     !query the optimal workspace
     lwork=-1
     CALL ZHEEV('V','U',dime,evec,lda,eval,work_tmp,lwork,rwork,info)
     !lwork=MIN(lwmax,INT(work(1)))
     lwork=MIN(2*dime,INT(work_tmp(1)))
     ALLOCATE(work(lwork))
     !solve engivalue
     CALL ZHEEV('V','U',dime,evec,lda,eval,work,lwork,rwork,info)
     IF(info/=0)THEN
        print*,'ZHEEV err:info',info
        STOP
     ENDIF
     !---destroy array
     DEALLOCATE(rwork,work)
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE lapk_cmplx_SEigen
  !#########################################################!
  !  CALL DSYEV(jobs,uplo,n,a,lda,w,work,lwork,info)  !
  !#########################################################!
  SUBROUTINE lapk_real_SEigen(mat,evec,eval)
     !
     IMPLICIT NONE
     REAL(DP),INTENT(IN)  :: mat(:,:) !IN:martix
     REAL(DP),INTENT(OUT) :: evec(:,:) !OUT:engivector
     REAL(DP),INTENT(OUT) :: eval(:)     !engivalue
     !
     INTEGER(I4B) :: lda            !see Lapack ZHEEV
     INTEGER(I4B) :: lwork !,lwmax=1000
     REAL(DP) :: wkt(10)
     REAL(DP),ALLOCATABLE :: work(:)
     INTEGER(I4B) :: info
     INTEGER(I4B) :: dime
     !
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     dime=SIZE(mat,1)
     IF(dime>MMAX)THEN
        WRITE(6,*) 'diagM:Check the dimension of MMAX'
        STOP
     ENDIF
     lda=dime
     !use mat for lapack
     evec(:,:)=mat(:,:)

     !query the optimal workspace
     lwork=-1
     CALL DSYEV('V','U',dime,evec,lda,eval,wkt,lwork,info)
     !lwork=MIN(lwmax,INT(work(1)))
     lwork=INT(wkt(1))
     !allocate
     ALLOCATE(work(lwork))

     !solve engivalue
     CALL DSYEV('V','U',dime,evec,lda,eval,work,lwork,info)
     IF(info/=0)THEN
        print*,'DSYEV err:info',info
        STOP
     ENDIF

     !---destroy array
     DEALLOCATE(work)
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE lapk_real_SEigen
  !----------------------------------------------------------
  !###################################################################!
  ! CALL ZHEGV(itype,JOBZ,uplo,N,A,lDA,B,LDB,w,work,lwork,rwork,INFO) !
  ! itype = 1 : Ax=eBx                                                !
  !         2 : ABx=ex                                                !
  !         3 : BAx=ex                                                !
  !###################################################################!
  SUBROUTINE lapk_cmplx_GEigen(dime,matA,matB,evec,eval)
     IMPLICIT NONE
     INTEGER(I4B),INTENT(IN) :: dime
     COMPLEX(DCP),INTENT(IN)  :: matA(:,:),matB(:,:)
     COMPLEX(DCP),INTENT(OUT) :: evec(:,:)
     REAL(DP),INTENT(OUT) :: eval(:)
     !LOCAL
     INTEGER(I4B) :: lda,ldb,info,lwork
     COMPLEX(DCP) :: U(dime,dime)
     COMPLEX(DCP) :: work_tmp(10)
     COMPLEX(DCP), ALLOCATABLE :: work(:)
     REAL(DP)    :: rwork(3*dime-2)
     !
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     IF(dime>MMAX)THEN
        WRITE(6,*) 'diagM:Check the dimension of MMAX'
        STOP
     ENDIF
     lda=dime
     ldb=dime
     !
     evec(:,:)=matA(:,:)
     U(:,:)=matB(:,:)
     !query the optimal workspace
     lwork=-1
     info=0
     CALL ZHEGV(1,'V','U',dime,evec,lda,U,ldb,eval,work_tmp,lwork,rwork,info)
     ! lwork=MIN(2*dime-1,INT(work_tmp(1)))
     lwork=MAX(2*dime-1,INT(work_tmp(1)))
     ALLOCATE(work (lwork))
     !solve the eigen problem
     CALL ZHEGV(1,'V','U',dime,evec,lda,U,ldb,eval,work,lwork,rwork,info)
     DEALLOCATE(work)
     !test
     IF(info/=0)THEN
        print*,'ZHEGV err:info',info
        STOP
     ENDIF
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE lapk_cmplx_GEigen
  !###################################################################!
  ! CALL DSYGV(itype,JOBZ,uplo,N,A,lDA,B,LDB,W,work,lwork,INFO)       !
  ! itype = 1 : Ax=eBx                                                !
  !         2 : ABx=ex                                                !
  !         3 : BAx=ex                                                !
  !###################################################################!
  SUBROUTINE lapk_real_GEigen(dime,matA,matB,evec,eval)
     IMPLICIT NONE
     INTEGER(I4B),INTENT(IN) :: dime
     REAL(DP),INTENT(IN)  :: matA(:,:),matB(:,:)
     REAL(DP),INTENT(OUT) :: evec(:,:)
     REAL(DP),INTENT(OUT) :: eval(:)
     !LOCAL
     INTEGER(I4B) :: lda,ldb,info,lwork
     REAL(DP) :: U(dime,dime)
     REAL(DP) :: wkt(10)
     REAL(DP),ALLOCATABLE :: work(:)
     !
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     IF(dime>MMAX)THEN
        WRITE(6,*) 'diagM:Check the dimension of MMAX'
        STOP
     ENDIF
     lda=dime
     ldb=dime
     !
     evec(:,:)=matA(:,:)
     U(:,:)=matB(:,:)
     !query the optimal workspace
     lwork=-1
     info=0
     CALL DSYGV(1,'V','U',dime,evec,lda,U,ldb,eval,wkt,lwork,info)
     lwork=INT(wkt(1))
     !allocate
     ALLOCATE(work(lwork))
     !solve the eigen problem
     CALL DSYGV(1,'V','U',dime,evec,lda,U,ldb,eval,work,lwork,info)
     !test
     IF(info/=0)THEN
        print*,'DSYGV err:info',info
        STOP
     ENDIF
     !destroy
     DEALLOCATE(work)
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE lapk_real_GEigen
  !#########################################################!
  !      CALL ZGEQRF(M,N,A,LDA,TAU,WORK,LWORK,INFO)         !
  !#########################################################!
   SUBROUTINE lapk_cmplx_OrthNorm(mat)
      !
      IMPLICIT NONE
      COMPLEX(DCP) :: mat(:,:)
      !LOCAL
      INTEGER(I4B) :: m,n,lda,lwork,info
      COMPLEX(DCP) :: tau(SIZE(mat,2))
      COMPLEX(DCP) :: work_tmp(10)
      COMPLEX(DCP),ALLOCATABLE :: work(:)
      !
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      m=SIZE(mat,1)
      n=SIZE(mat,2)
      lda=MAX(1,m)
      !
      !query the optimal workspace
      lwork=-1
      CALL ZGEQRF(m,n,mat,lda,tau,work,lwork,info)
      lwork=MAX(1,INT(work(1)))
      ALLOCATE(work(lwork))
      !QR
      CALL ZGEQRF(m,n,mat,lda,tau,work,lwork,info)
      IF(info/=0)THEN
        print*,'ZGEQRF err:info',info
        STOP
      ENDIF
      !orth-normal
      CALL ZUNGQR(m,n,n,mat,lda,tau,work,lwork,info)
      IF(info/=0)THEN
        print*,'ZUNGQR err:info',info
        STOP
      ENDIF
      !
      DEALLOCATE(work)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE lapk_cmplx_OrthNorm
  !#########################################################!
  !      CALL DGEQRF(M,N,A,LDA,TAU,WORK,LWORK,INFO)         !
  !#########################################################!
   SUBROUTINE lapk_real_OrthNorm(mat)
      !
      IMPLICIT NONE
      REAL(DP) :: mat(:,:)
      !LOCAL
      INTEGER(I4B) :: m,n,lda,lwork,info
      REAL(DP) :: tau(SIZE(mat,2))
      !
      REAL(DP) :: wkt(10)
      REAL(DP),ALLOCATABLE :: work(:)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !
      m=SIZE(mat,1)
      n=SIZE(mat,2)
      lda=MAX(1,m)
      !
      !!!!xlt!IF(m>MMAX .OR. n>NMAX .OR. m<n) THEN
      !!!!xlt!   print*,'m,n',m,n
      !!!!xlt!   WRITE(6,*) 'OrthNorm:Check the dimension of MMAX,NMAX,M,N'
      !!!!xlt!   STOP
      !!!!xlt!ENDIF
      !query the optimal workspace
      lwork=-1
      CALL DGEQRF(m,n,mat,lda,tau,wkt,lwork,info)
      lwork=INT(wkt(1))
      !ALLOCATE
      ALLOCATE(work(lwork))
      !QR
      CALL DGEQRF(m,n,mat,lda,tau,work,lwork,info)
      IF(info/=0)THEN
        print*,'DGEQRF err:info',info
        STOP
      ENDIF
      !orth-normal
      CALL DORGQR(m,n,n,mat,lda,tau,work,lwork,info)
      IF(info/=0)THEN
        print*,'DORGQR err:info',info
        STOP
      ENDIF
      !Destory
      DEALLOCATE(work)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE lapk_real_OrthNorm
  !-------------------------Norm2----------------------------
  FUNCTION lapk_cmplx_Norm2(mat,k)
     !
     IMPLICIT NONE
     COMPLEX(DCP),INTENT(IN) :: mat(:,:)
     INTEGER(I4B),INTENT(IN) :: k
     REAL(DP) :: lapk_cmplx_Norm2
     !
     COMPLEX(DCP) :: MM(k,k),evec(k,k)
     REAL(DP) :: eval(k)
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     MM(:,:)=MATMUL(TRANSPOSE(CONJG(mat)),mat)
     CALL lapk_cmplx_SEigen(MM,evec,eval)
     !
     lapk_cmplx_Norm2=SQRT((eval(k)))
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDFUNCTION lapk_cmplx_Norm2
  !#############################################################!
  ! matrix: C = op(A) * op(B)                                   !
  ! CALL ZGEMM(TRANSA,TRANSB,M,N,K,alpha,A,LDA,B,LDB,beta,c,LDC)!
  ! op*='N','T','C' for op(X)=(X), (X)' and (X*)'               !
  !#############################################################!
  !-------------------------matmat---------------------------
   SUBROUTINE lapk_cmplx_matmat(matA,matB,opA,opB,cab,cc,matC)
      IMPLICIT NONE
      COMPLEX(DCP),INTENT(IN)  :: matA(:,:),matB(:,:)
      COMPLEX(DCP),INTENT(IN) :: cab,cc
      CHARACTER(1),INTENT(IN) :: opA,opB
      COMPLEX(DCP),INTENT(INOUT) :: matC(:,:)
      !
      INTEGER(I4B) :: LDA,LDB,LDC,M,N,K
      COMPLEX(DCP)  :: alpha  &
               &    ,  beta
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !
      alpha=cab
      beta=cc
      !C
      LDC=SIZE(matC,1)
      M=LDC
      N=SIZE(matC,2)
      !A
      IF((opA=='N').OR.(opA=='n'))THEN
         K=SIZE(matA,2)
         LDA=MAX(1,M)
      ELSE
         K=SIZE(matA,1)
         LDA=MAX(1,K)
      ENDIF
      !B
      IF((opB=='N').OR.(opB=='n'))THEN
         LDB=MAX(1,K)
      ELSE
         LDB=MAX(1,N)
      ENDIF
      !
      CALL ZGEMM(opA,opB,M,N,K,alpha,matA,LDA,matB,LDB,beta,matC,LDC)
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE lapk_cmplx_matmat
  !#############################################################!
  ! matrix: C = op(A) * op(B)                                   !
  ! CALL DGEMM(TRANSA,TRANSB,M,N,K,alpha,A,LDA,B,LDB,beta,c,LDC)!
  ! op*='N','T','C' for op(X)=(X), (X)' and (X*)'               !
  !#############################################################!
   SUBROUTINE lapk_real_matmat(matA,matB,opA,opB,cab,cc,matC)
      IMPLICIT NONE
      REAL(DP),INTENT(IN)  :: matA(:,:),matB(:,:)
      REAL(DP),INTENT(IN) :: cab,cc
      REAL(DP),INTENT(OUT) :: matC(:,:)
      CHARACTER(1),INTENT(IN) :: opA,opB
      !
      INTEGER(I4B) :: LDA,LDB,LDC,M,N,K
      REAL(DP)  :: alpha  &
               &    ,  beta 
      !
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      alpha=cab
      beta=cc
      !C
      matC(:,:)=0.d0
      LDC=SIZE(matC,1)
      M=LDC
      N=SIZE(matC,2)
      !A
      IF((opA=='N').OR.(opA=='n'))THEN
         K=SIZE(matA,2)
         LDA=MAX(1,M)
      ELSE
         K=SIZE(matA,1)
         LDA=MAX(1,K)
      ENDIF
      !B
      IF((opB=='N').OR.(opB=='n'))THEN
         LDB=MAX(1,K)
      ELSE
         LDB=MAX(1,N)
      ENDIF
      !
      CALL DGEMM(opA,opB,M,N,K,alpha,matA,LDA,matB,LDB,beta,matC,LDC)
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE lapk_real_matmat
  !-------------------------matvec---------------------------
  SUBROUTINE lapk_cmplx_matvec(dime,matA,matB,opA,opB,cab,cc,matC)
     IMPLICIT NONE
     INTEGER(I4B),INTENT(IN) :: dime
     COMPLEX(DCP),INTENT(IN)  :: matA(dime,dime),matB(dime,1)
     COMPLEX(DCP),INTENT(IN) :: cab,cc
     CHARACTER(1),INTENT(IN) :: opA,opB
     COMPLEX(DCP),INTENT(INOUT) :: matC(dime,1)
     !
     INTEGER(I4B) :: LDA,LDB,LDC,M,N,K
     COMPLEX(DCP)  :: alpha &  !=CMPLX(1._DP,0._DP)  &
             &    ,  beta !=CMPLX(0._DP,0._DP)
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     !
     alpha=cab
     beta=cc
     !C
     LDC=SIZE(matC,1)
     M=LDC
     N=SIZE(matC,2)
     !A
     IF((opA=='N').OR.(opA=='n'))THEN
        K=SIZE(matA,2)
        LDA=MAX(1,M)
     ELSE
        K=SIZE(matA,1)
        LDA=MAX(1,K)
     ENDIF
     !B
     IF((opB=='N').OR.(opB=='n'))THEN
        LDB=MAX(1,K)
     ELSE
        LDB=MAX(1,N)
     ENDIF
     !
     CALL ZGEMM(opA,opB,M,N,K,alpha,matA,LDA,matB,LDB,beta,matC,LDC)
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE lapk_cmplx_matvec
  !-------------------------matvec---------------------------
  SUBROUTINE lapk_real_matvec(dime,matA,vecB,opA,opB,cab,cc,vecC)
     IMPLICIT NONE
     INTEGER(I4B),INTENT(IN) :: dime
     REAL(DP),INTENT(IN)  :: matA(dime,dime),vecB(dime,1) &
              &, cab,cc
     REAL(DP),INTENT(OUT) :: vecC(dime,1)
     CHARACTER(1),INTENT(IN) :: opA,opB
     !
     INTEGER(I4B) :: LDA,LDB,LDC,M,N,K
     REAL(DP)  :: alpha & !=1.d0  &
             &    ,  beta !=0.d0
     !
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     alpha=cab
     beta=cc
     !C
     vecC(:,:)=0.d0
     LDC=SIZE(vecC,1)
     M=LDC
     N=SIZE(vecC,2)
     !A
     IF((opA=='N').OR.(opA=='n'))THEN
        K=SIZE(matA,2)
        LDA=MAX(1,M)
     ELSE
        K=SIZE(matA,1)
        LDA=MAX(1,K)
     ENDIF
     !B
     IF((opB=='N').OR.(opB=='n'))THEN
        LDB=MAX(1,K)
     ELSE
        LDB=MAX(1,N)
     ENDIF
     !
     CALL DGEMM(opA,opB,M,N,K,alpha,matA,LDA,vecB,LDB,beta,vecC,LDC)
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE lapk_real_matvec
  !-------------------------invmat---------------------------
  SUBROUTINE lapk_cmplx_invmat(mat)
     IMPLICIT NONE
     COMPLEX(DCP),INTENT(INOUT) :: mat(:,:)
     !LOCAL
     INTEGER(I4B) :: info,lda,m,n
     INTEGER(I4B) :: IPIV(SIZE(mat,1))
     !
     INTEGER(I4B) :: lwork
     COMPLEX(DCP) :: wkt(10)
     COMPLEX(DCP),ALLOCATABLE :: work(:)
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     !
     n=SIZE(mat,1)
     m=n
     lda=n
     !LU decomposed
     CALL ZGETRF(m,n,mat,lda,IPIV,info)
     IF(info/=0)THEN
        print*,'DGETRF err:info',info
     ENDIF
     !inverse
     !inquire the work matrix
     lwork=-1
     CALL ZGETRI(n,mat,lda,IPIV,wkt,lwork,info)
     lwork=INT(REAL(wkt(1),DP))
     !allocate
     ALLOCATE(work(lwork))
     CALL ZGETRI(n,mat,lda,IPIV,work,lwork,info)
     IF(info/=0)THEN
        print*,'ZGETRF err:info',info
     ENDIF
     !deallocate
     DEALLOCATE(work) 
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE lapk_cmplx_invmat
  !-------------------------invmat---------------------------
  SUBROUTINE lapk_real_invmat(mat)
     IMPLICIT NONE
     REAL(DP),INTENT(INOUT) :: mat(:,:)
     !LOCAL
     INTEGER(I4B) :: info,lda,m,n
     INTEGER(I4B) :: IPIV(SIZE(mat,1))
     !
     INTEGER(I4B) :: lwork
     REAL(DP) :: wkt(10)
     REAL(DP),ALLOCATABLE :: work(:)
     !
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     n=SIZE(mat,1)
     m=n
     lda=n
     !LU decomposed
     CALL DGETRF(m,n,mat,lda,IPIV,info)
     IF(info/=0)THEN
        print*,'DGETRF err:info',info
     ENDIF
     !inverse
     !inquire the work matrix
     lwork=-1
     CALL DGETRI(n,mat,lda,IPIV,wkt,lwork,info)
     lwork=wkt(1)
     !allocate
     ALLOCATE(work(lwork))
     CALL DGETRI(n,mat,lda,IPIV,work,lwork,info)
     IF(info/=0)THEN
        print*,'DGETRF err:info',info
     ENDIF
     !deallocate
     DEALLOCATE(work) 
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE lapk_real_invmat
  !------------------selected eigenvalues--------------------
  SUBROUTINE lapk_real_SEigenX(mat,dime,num,Il,Iu,evec,eval)
     IMPLICIT NONE
     !INOUT
     INTEGER(I4B),INTENT(IN) :: Il & !the lower index of eigen-pair
                            & , Iu,num & !the upper index of eigen-pair
                            & , dime !the dimension of mat
     REAL(DP),INTENT(IN) :: mat(dime,dime) !matrix
     REAL(DP),INTENT(OUT) :: evec(dime,num) &!eigen vector
                         &,  eval(num)   !eigen value
     !LOCAL
     REAL(DP) :: A(dime,dime)
     REAL(DP) :: vl,vu,wkt(10),w(dime) !,z(dime,num)
     INTEGER(I4B) :: lda,lwork,iwork(5*dime),m
     REAL(DP) :: abstol=1e-6
     INTEGER(I4B) :: ifail(dime),info,ldz
     REAL(DP),ALLOCATABLE :: work(:)
     !
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     A(:,:)=mat(:,:)
     lda=dime
     ldz=dime
     m=num
     !query the optimal workspace
     lwork=-1
     CALL DSYEVX('V','I','U',dime,A,lda,vl,vu,Il,Iu,abstol &
         & ,m,w,evec,ldz,wkt,lwork,iwork,ifail,info)
     lwork=INT(wkt(1))
     !work array
     ALLOCATE(work(lwork))
     CALL DSYEVX('V','I','U',dime,A,lda,vl,vu,Il,Iu,abstol &
         & ,m,w,evec,ldz,work,lwork,iwork,ifail,info)
     IF(info/=0)THEN
        print*,'DSYEVX err:info',info
        STOP
     ENDIF
     !destory
     !DEALLOCATE(work)
     !-------------------OUT PUT-----------------------------
     !eigen-values
     eval(1:num)=w(1:num)
     DEALLOCATE(work)
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE lapk_real_SEigenX
  !#########################################################!
  !          CALL DPOTRF(UPLO,N,A,LDA,INFO)                 !
  !#########################################################!
  !-------------------Cholesky factorization-----------------
  SUBROUTINE lapk_real_CholeskyF(mat)
     IMPLICIT NONE
     !IN/OUT
     REAL(DP),INTENT(INOUT) :: mat(:,:)
     !LOCAL
     INTEGER(I4B) :: info,lda,dimen
     !
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     dimen=SIZE(mat,1)
     lda=dimen
     !
     CALL DPOTRF('U',dimen,mat,lda,info)
     !test
     IF(info/=0)THEN
        print*,'Cholesky err:info',info
        STOP
     ENDIF
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE lapk_real_CholeskyF
  !-----------------------inv_mat----------------------------
  !----------------------------------------------------------
   ELEMENTAL FUNCTION UPCASE(c) RESULT(res)
      CHARACTER(1), INTENT(IN) :: c
      CHARACTER(1) :: res
      IF (c >= 'a' .AND. c <= 'z') THEN
         res = ACHAR(IACHAR(c) - 32)
      ELSE
         res = c
      ENDIF
   ENDFUNCTION UPCASE
   !----------------------------------------------------------
ENDMODULE Lapack_module


#ifdef MPI
MODULE ScaLapack_module
   USE mpi
   USE smpi_math_module
   USE parameters, ONLY: nev
   USE constants
   IMPLICIT NONE
   INTEGER(I4B) :: MYROW, MYCOL
   INTEGER(I4B) :: NPROW, NPCOL
   INTEGER(I4B) :: blacs_context
   INTEGER(I4B) :: info, ierr
   LOGICAL,SAVE    :: INIT_CALLED = .false.
   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
   !OrthNorm
   interface scalapk_OrthNorm!{{{
      module procedure scalapk_real_OrthNorm
      module procedure scalapk_cmplx_OrthNorm
   end interface scalapk_OrthNorm!}}}
   !MatMat
   interface scalapk_MM!{{{
      module procedure scalapk_real_matmat
      module procedure scalapk_cmplx_matmat
   end interface scalapk_MM!}}}
CONTAINS
!------------------------------------------------------------------------
   SUBROUTINE Init_scala() 
      IMPLICIT NONE       
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
      CALL Finalize_scala()  !In case called twice   
      blacs_context = parallel%commy           
      NPROW = 1    
      NPCOL = parallel%commy_numprocs
      !CALL BLACS_GET( -1, 0, blacs_context ) 
      CALL BLACS_GRIDINIT( blacs_context, 'Col-major', NPROW, NPCOL )    
      CALL BLACS_GRIDINFO( blacs_context, NPROW, NPCOL, MYROW, MYCOL )
      INIT_CALLED=.true.
   ENDSUBROUTINE Init_scala
!------------------------------------------------------------------------
   SUBROUTINE Finalize_scala() 
      IMPLICIT NONE       
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    
      IF (INIT_CALLED) THEN
         CALL BLACS_GRIDEXIT( blacs_context )    
         INIT_CALLED=.false.    
      ENDIF
   ENDSUBROUTINE Finalize_scala
!------------------------------------------------------------------------
   SUBROUTINE scalapk_cmplx_OrthNorm(mat)
      USE parameters, ONLY:nev
      IMPLICIT NONE
      !inout
      COMPLEX(DCP),INTENT(INOUT)    :: mat(:,:)
      !ScaLapack
      COMPLEX(DCP), ALLOCATABLE     :: TAU(:), WORK(:)
      INTEGER(I4B)                  :: MB, NB, M, N, K, LLD
      INTEGER(I4B)                  :: desc_mat(9)
      INTEGER(I4B)                  :: s_local
      INTEGER(I4B)                  :: LWORK
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF (.NOT. INIT_CALLED) THEN
         PRINT*, 'Error: lapk_cmplx_OrthNorm: Init_scala not called!'
         CALL MPI_ABORT(parallel%commy, 1, ierr)
      ENDIF
      IF (MOD(nev, parallel%commy_numprocs) /= 0) THEN
         IF (parallel%isroot) PRINT*, 'Error: nev=',nev,'Cannot be processed by the ncore=',parallel%commy_numprocs,'divided'
         CALL MPI_ABORT(parallel%commy, 1, ierr)
         RETURN
      ENDIF
      !
      s_local = parallel%nstate_proc
      NPCOL = parallel%commy_numprocs
      NB = s_local
      MB = SIZE(mat,1)
      M = SIZE(mat,1)
      N = nev
      K = MIN(M,N)
      LLD = M
      !
      IF (SIZE(mat, 2) /= s_local) THEN
         PRINT*, 'Rank', parallel%commy_myid, ': mat row err.expected =', s_local, 'actually =', SIZE(mat,2)
         CALL MPI_ABORT(parallel%commy, 99, ierr)
      ENDIF
      CALL DESCINIT(desc_mat, M, N, MB, NB, 0, 0, blacs_context, LLD, info)
      IF (info /= 0) THEN
         PRINT*, 'Rank',parallel%commy_myid,': desc init err,info=', info
         CALL MPI_ABORT(parallel%commy, 7, ierr)
         CALL BLACS_GRIDEXIT(blacs_context)
         RETURN
      ENDIF
      !
      ALLOCATE(TAU(nev), STAT=info)
      IF (info /= 0) THEN
         PRINT*, 'Rank',parallel%commy_myid,': TAU allocated err,info=', info
         CALL MPI_ABORT(parallel%commy, 8, ierr)
         CALL BLACS_GRIDEXIT(blacs_context)
         RETURN
      ENDIF
      ALLOCATE(WORK(1), STAT=info)
      IF (info /= 0) THEN
         PRINT*, 'Rank',parallel%commy_myid,': WORK allocated err,info=', info
         CALL MPI_ABORT(parallel%commy, 9, ierr)
         DEALLOCATE(TAU, STAT=info)
         CALL BLACS_GRIDEXIT(blacs_context)
         RETURN
      ENDIF

      LWORK = -1
      CALL PZGEQRF(M, N, mat, 1, 1, desc_mat, TAU, WORK, LWORK, info)
      IF (info /= 0 .AND. info /= -1000) THEN
         PRINT*, 'Rank',parallel%commy_myid,': PZGEQRF failed to query,info=', info
         CALL MPI_ABORT(parallel%commy, 10, ierr)
         DEALLOCATE(TAU, WORK, STAT=info)
         CALL BLACS_GRIDEXIT(blacs_context)
         RETURN
      ENDIF
      LWORK = MAX(1, INT(WORK(1)))
      DEALLOCATE(WORK)
      ALLOCATE(WORK(LWORK), STAT=info)
      IF (info /= 0) THEN
         PRINT*, 'Rank',parallel%commy_myid,': WORK allocated err,info=', info
         CALL MPI_ABORT(parallel%commy, 11, ierr)
         DEALLOCATE(TAU, STAT=info)
         CALL BLACS_GRIDEXIT(blacs_context)
         RETURN
      ENDIF
      !
      CALL PZGEQRF(M, N, mat, 1, 1, desc_mat, TAU, WORK, LWORK, info)
      IF (info /= 0) THEN
         PRINT*, 'Rank',parallel%commy_myid,': PZGEQRF QR err,info=', info
         CALL MPI_ABORT(parallel%commy, 12, ierr)
         DEALLOCATE(TAU, STAT=info)
         CALL BLACS_GRIDEXIT(blacs_context)
         RETURN
      ENDIF
      !
      CALL PZUNGQR(M, N, K, mat, 1, 1, desc_mat, TAU, WORK, LWORK, info)
      IF (info /= 0) THEN
         PRINT*, 'Rank',parallel%commy_myid,': PZUNGQR err,info=', info
         CALL MPI_ABORT(parallel%commy, 13, ierr)
         DEALLOCATE(TAU, STAT=info)
         CALL BLACS_GRIDEXIT(blacs_context)
         RETURN
      ENDIF
      !
      IF (ALLOCATED(TAU)) DEALLOCATE(TAU, STAT=info)
      IF (ALLOCATED(WORK)) DEALLOCATE(WORK, STAT=info)
   ENDSUBROUTINE scalapk_cmplx_OrthNorm
!------------------------------------------------------------------------
   SUBROUTINE scalapk_real_OrthNorm(mat)
      USE parameters, ONLY:nev
      IMPLICIT NONE
      !inout
      REAL(DP),INTENT(INOUT)        :: mat(:,:)
      !ScaLapack
      REAL(DP), ALLOCATABLE         :: TAU(:), WORK(:)
      INTEGER(I4B)                  :: MB, NB, M, N, K, LLD
      INTEGER(I4B)                  :: desc_mat(9)
      INTEGER(I4B)                  :: s_local
      INTEGER(I4B)                  :: LWORK
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF (.NOT. INIT_CALLED) THEN
         PRINT*, 'Error: lapk_real_OrthNorm: Init_scala not called!'
         CALL MPI_ABORT(parallel%commy, 1, ierr)
      ENDIF
      IF (MOD(nev, parallel%commy_numprocs) /= 0) THEN
         IF (parallel%isroot) PRINT*, 'Error: nev=',nev,'Cannot be processed by the ncore=',parallel%commy_numprocs,'divided'
         CALL MPI_ABORT(parallel%commy, 1, ierr)
         RETURN
      ENDIF
      !
      s_local = parallel%nstate_proc
      NPCOL = parallel%commy_numprocs
      NB = s_local
      MB = SIZE(mat,1)
      M = SIZE(mat,1)
      N = nev
      K = MIN(M,N)
      LLD = M
      !
      IF (SIZE(mat, 2) /= s_local) THEN
         PRINT*, 'Rank', parallel%commy_myid, ': mat row err.expected =', s_local, 'actually =', SIZE(mat,2)
         CALL MPI_ABORT(parallel%commy, 99, ierr)
      ENDIF
      CALL DESCINIT(desc_mat, M, N, MB, NB, 0, 0, blacs_context, LLD, info)
      IF (info /= 0) THEN
         PRINT*, 'Rank',parallel%commy_myid,': desc init err,info=', info
         CALL MPI_ABORT(parallel%commy, 7, ierr)
         CALL BLACS_GRIDEXIT(blacs_context)
         RETURN
      ENDIF
      !
      ALLOCATE(TAU(nev), STAT=info)
      IF (info /= 0) THEN
         PRINT*, 'Rank',parallel%commy_myid,': TAU allocated err,info=', info
         CALL MPI_ABORT(parallel%commy, 8, ierr)
         CALL BLACS_GRIDEXIT(blacs_context)
         RETURN
      ENDIF
      ALLOCATE(WORK(1), STAT=info)
      IF (info /= 0) THEN
         PRINT*, 'Rank',parallel%commy_myid,': WORK allocated err,info=', info
         CALL MPI_ABORT(parallel%commy, 9, ierr)
         DEALLOCATE(TAU, STAT=info)
         CALL BLACS_GRIDEXIT(blacs_context)
         RETURN
      ENDIF

      LWORK = -1
      CALL PDGEQRF(M, N, mat, 1, 1, desc_mat, TAU, WORK, LWORK, info)
      IF (info /= 0 .AND. info /= -1000) THEN
         PRINT*, 'Rank',parallel%commy_myid,': PDGEQRF failed to query,info=', info
         CALL MPI_ABORT(parallel%commy, 10, ierr)
         DEALLOCATE(TAU, WORK, STAT=info)
         CALL BLACS_GRIDEXIT(blacs_context)
         RETURN
      ENDIF
      LWORK = MAX(1, INT(WORK(1)))
      DEALLOCATE(WORK)
      ALLOCATE(WORK(LWORK), STAT=info)
      IF (info /= 0) THEN
         PRINT*, 'Rank',parallel%commy_myid,': WORK allocated err,info=', info
         CALL MPI_ABORT(parallel%commy, 11, ierr)
         DEALLOCATE(TAU, STAT=info)
         CALL BLACS_GRIDEXIT(blacs_context)
         RETURN
      ENDIF
      !
      CALL PDGEQRF(M, N, mat, 1, 1, desc_mat, TAU, WORK, LWORK, info)
      IF (info /= 0) THEN
         PRINT*, 'Rank',parallel%commy_myid,': PDGEQRF QR err,info=', info
         CALL MPI_ABORT(parallel%commy, 12, ierr)
         DEALLOCATE(TAU, STAT=info)
         CALL BLACS_GRIDEXIT(blacs_context)
         RETURN
      ENDIF
      !
      CALL PDORGQR(M, N, K, mat, 1, 1, desc_mat, TAU, WORK, LWORK, info)
      IF (info /= 0) THEN
         PRINT*, 'Rank',parallel%commy_myid,': PDUNGQR err,info=', info
         CALL MPI_ABORT(parallel%commy, 13, ierr)
         DEALLOCATE(TAU, STAT=info)
         CALL BLACS_GRIDEXIT(blacs_context)
         RETURN
      ENDIF
      !
      IF (ALLOCATED(TAU)) DEALLOCATE(TAU, STAT=info)
      IF (ALLOCATED(WORK)) DEALLOCATE(WORK, STAT=info)
   ENDSUBROUTINE scalapk_real_OrthNorm
!-------------------------------------------------------------------------
   SUBROUTINE scalapk_cmplx_matmat(matA,matB,opA,opB,cab,cc,matC,matC_glob)
      IMPLICIT NONE
      !in/out
      COMPLEX(DCP), INTENT(IN)  :: matA(:,:), matB(:,:)
      CHARACTER(1), INTENT(IN)  :: opA, opB
      COMPLEX(DCP), INTENT(IN)  :: cab, cc
      COMPLEX(DCP), INTENT(INOUT) :: matC(:,:)
      COMPLEX(DCP),OPTIONAL,INTENT(OUT) :: matC_glob(:,:)
      !local
      COMPLEX(DCP)   :: alpha, beta
      !ScaLapack
      INTEGER(I4B) :: M_A, N_A, MB_A, NB_A, LLD_A
      INTEGER(I4B) :: M_B, N_B, MB_B, NB_B, LLD_B
      INTEGER(I4B) :: M_C, N_C, MB_C, NB_C, LLD_C
      INTEGER(I4B) :: M, N ,K
      INTEGER(I4B) :: desc_matA(9), desc_matB(9), desc_matC(9)
      CHARACTER(1) :: opA_up, opB_up
      INTEGER(I4B) :: desc_matC_glob(9)
      INTEGER(I4B) :: IA, JA, IB, JB
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      alpha = cab
      beta = cc
      matC(:,:) = CMPLX(0.d0, 0.d0, KIND=DCP)
      opA_up = UPCASE(opA)
      opB_up = UPCASE(opB)
      IF (.NOT. INIT_CALLED) THEN
         PRINT*, 'Error: lapk_cmplx_OrthNorm: Init_scala not called!'
         CALL MPI_ABORT(parallel%commy, 1, ierr)
      ENDIF
      IF (.NOT. (opA_up == 'N' .OR. opA_up == 'T' .OR. opA_up == 'C')) THEN
         WRITE(*, '(A,I0,A,A)') 'commy proc ', parallel%commy_myid, ': Invalid opA=', opA
         info = -1
         GOTO 999
      ENDIF
      IF (.NOT. (opB_up == 'N' .OR. opB_up == 'T' .OR. opB_up == 'C')) THEN
         WRITE(*, '(A,I0,A,A)') 'commy proc ', parallel%commy_myid, ': Invalid opB=', opB
         info = -2
         GOTO 999
      ENDIF
      IF (MOD(nev, NPCOL) /= 0) THEN
         info = -3
         WRITE(*, '(A,I0,A,I0,A,I0)') 'commy proc ', parallel%commy_myid,': nev=',nev,' not divisible by NPCOL=',NPCOL
         GOTO 999
      ENDIF
      !
      SELECT CASE(opA_up)
         CASE('N')
            M = SIZE(matA, 1, KIND=I4B)
         CASE('T','C')
            M = nev
      END SELECT
      SELECT CASE(opB_up)
         CASE('N')
            N = nev
            K = SIZE(matB, 1, KIND=I4B) 
         CASE('T','C')
            N  = SIZE(matB, 1, KIND=I4B)
            K = nev
      END SELECT
      M_A = SIZE(matA, 1, KIND=I4B)
      N_A = nev
      MB_A = M_A
      NB_A = parallel%nstate_proc
      LLD_A = SIZE(matA, 1, KIND=I4B)
      CALL DESCINIT(desc_matA, M_A, N_A, MB_A, NB_A, 0, 0, blacs_context, LLD_A, info)
      IF (info /= 0) THEN 
         WRITE(*,*) 'proc ', parallel%myid,': DESCINIT matA failed, info=',info; GOTO 999 
      ENDIF
      M_B = SIZE(matB, 1, KIND=I4B)
      N_B = nev
      MB_B = M_B
      NB_B = parallel%nstate_proc
      LLD_B = SIZE(matB, 1, KIND=I4B)
      CALL DESCINIT(desc_matB, M_B, N_B, MB_B, NB_B, 0, 0, blacs_context, LLD_B, info)
      IF (info /= 0) THEN
         WRITE(*,*) 'proc ', parallel%myid,': DESCINIT matB failed, info=',info; GOTO 999
      ENDIF
      M_C = M
      N_C = N
      MB_C = M_C
      NB_C = parallel%nstate_proc
      LLD_C = M_C
      IF (SIZE(matC,1,KIND=I4B) /= LLD_C .OR. SIZE(matC,2,KIND=I4B) /= NB_C) THEN
         info = -5
         WRITE(*, '(A,I0,A,I0,A,I0,A,I0,A,I0)') 'commy proc ', parallel%commy_myid, &
            ': matC dim (', SIZE(matC,1), ',', SIZE(matC,2), ') ≠ expected (', LLD_C, ',', NB_C, ')'
         GOTO 999
      ENDIF
      CALL DESCINIT(desc_matC, M_C, N_C, MB_C, NB_C, 0, 0, blacs_context, LLD_C, info)
      IF (info /= 0) THEN 
         WRITE(*,*) 'commy proc ', parallel%commy_myid,': DESCINIT matC failed, info=',info; GOTO 999 
      ENDIF
      !
      CALL PZGEMM(opA_up, opB_up, M, N, K, alpha, &
                  matA, 1_I4B, 1_I4B, desc_matA, &
                  matB, 1_I4B, 1_I4B, desc_matB, &
                  beta, matC, 1_I4B, 1_I4B, desc_matC)
      !
      IF (PRESENT(matC_glob)) THEN
         IF (SIZE(matC_glob,1) /= M_C .OR. SIZE(matC_glob,2) /= N_C) THEN
               info = -1; WRITE(*,*) 'matC_glob dim mismatch: expect (',M_C,',',N_C,') got (',SIZE(matC_glob,1),',',SIZE(matC_glob,2),')'; GOTO 999
         ENDIF
         CALL DESCINIT(desc_matC_glob, M_C, N_C, M_C, N_C, 0, 0, blacs_context, M_C, info)
         IF (info /= 0) THEN
            WRITE(*,*) 'DESCINIT desc_matC_glob failed at commy root, info=',info
            GOTO 999
         ENDIF
         IA = 1; JA = 1
         IB = 1; JB = 1
         CALL PZGEMR2D(M_C, N_C, &
                        matC, IA, JA, desc_matC, &
                        matC_glob, IB, JB, desc_matC_glob, &
                        blacs_context)
         CALL MPI_BCAST(matC_glob, M_C*N_C, MPI_DOUBLE_COMPLEX, &
                        0, parallel%commy, ierr)
      ENDIF
      !
      999 CONTINUE
      IF (info /= 0) THEN
         CALL BLACS_GRIDEXIT(blacs_context)
         WRITE(*, '(A,I0,A,I0)') 'proc ', parallel%myid, ': lapk_cmplx_matmat failed, final info=', info
      ENDIF
   ENDSUBROUTINE scalapk_cmplx_matmat
!------------------------------------------------------------------------
   SUBROUTINE scalapk_real_matmat(matA,matB,opA,opB,cab,cc,matC,matC_glob)
      IMPLICIT NONE
      !inout
      REAL(DP),INTENT(IN)  :: matA(:,:),matB(:,:)
      REAL(DP),INTENT(IN) :: cab,cc
      REAL(DP),INTENT(INOUT) :: matC(:,:)
      CHARACTER(1),INTENT(IN) :: opA,opB
      REAL(DP),OPTIONAL,INTENT(OUT) :: matC_glob(:,:)
      !local
      REAL(DP)       :: alpha, beta
      !ScaLapack
      INTEGER(I4B) :: M_A, N_A, MB_A, NB_A, LLD_A
      INTEGER(I4B) :: M_B, N_B, MB_B, NB_B, LLD_B
      INTEGER(I4B) :: M_C, N_C, MB_C, NB_C, LLD_C
      INTEGER(I4B) :: M, N ,K
      INTEGER(I4B) :: desc_matA(9), desc_matB(9), desc_matC(9)
      CHARACTER(1) :: opA_up, opB_up
      INTEGER(I4B) :: desc_matC_glob(9)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      alpha=cab
      beta=cc
      matC(:,:)=0.d0
      opA_up = UPCASE(opA)
      opB_up = UPCASE(opB)
      IF (.NOT. INIT_CALLED) THEN
         PRINT*, 'Error: lapk_cmplx_OrthNorm: Init_scala not called!'
         CALL MPI_ABORT(parallel%commy, 1, ierr)
      ENDIF
      IF (.NOT. (opA_up == 'N' .OR. opA_up == 'T' .OR. opA_up == 'C')) THEN
         WRITE(*, '(A,I0,A,A)') 'commy proc ', parallel%commy_myid, ': Invalid opA=', opA
         info = -1
         GOTO 999
      ENDIF
      IF (.NOT. (opB_up == 'N' .OR. opB_up == 'T' .OR. opB_up == 'C')) THEN
         WRITE(*, '(A,I0,A,A)') 'commy proc ', parallel%commy_myid, ': Invalid opB=', opB
         info = -2
         GOTO 999
      ENDIF
      !
      IF (MOD(nev, NPCOL) /= 0) THEN
         info = -3
         WRITE(*, '(A,I0,A,I0,A,I0)') 'commy proc ', parallel%commy_myid,': nev=',nev,' not divisible by NPCOL=',NPCOL
         GOTO 999
      ENDIF
      !
      SELECT CASE(opA_up)
         CASE('N')
            M = SIZE(matA, 1, KIND=I4B)
         CASE('T','C')
            M = nev
      END SELECT
      SELECT CASE(opB_up)
         CASE('N')
            N = nev
            K = SIZE(matB, 1, KIND=I4B) 
         CASE('T','C')
            N  = SIZE(matB, 1, KIND=I4B)
            K = nev
      END SELECT
      M_A = SIZE(matA, 1, KIND=I4B)
      N_A = nev
      MB_A = M_A
      NB_A = parallel%nstate_proc
      LLD_A = SIZE(matA, 1, KIND=I4B)
      CALL DESCINIT(desc_matA, M_A, N_A, MB_A, NB_A, 0, 0, blacs_context, LLD_A, info)
      IF (info /= 0) THEN 
         WRITE(*,*) 'proc ', parallel%myid,': DESCINIT matA failed, info=',info; GOTO 999 
      ENDIF
      M_B = SIZE(matB, 1, KIND=I4B)
      N_B = nev
      MB_B = M_B
      NB_B = parallel%nstate_proc
      LLD_B = SIZE(matB, 1, KIND=I4B)
      CALL DESCINIT(desc_matB, M_B, N_B, MB_B, NB_B, 0, 0, blacs_context, LLD_B, info)
      IF (info /= 0) THEN
         WRITE(*,*) 'proc ', parallel%myid,': DESCINIT matB failed, info=',info; GOTO 999
      ENDIF
      M_C = M
      N_C = N
      MB_C = M_C
      NB_C = parallel%nstate_proc
      LLD_C = M_C
      IF (SIZE(matC,1,KIND=I4B) /= LLD_C .OR. SIZE(matC,2,KIND=I4B) /= NB_C) THEN
         info = -5
         WRITE(*, '(A,I0,A,I0,A,I0,A,I0,A,I0)') 'proc ', parallel%myid, &
            ': matC dim (', SIZE(matC,1), ',', SIZE(matC,2), ') ≠ expected (', LLD_C, ',', NB_C, ')'
         GOTO 999
      ENDIF
      CALL DESCINIT(desc_matC, M_C, N_C, MB_C, NB_C, 0, 0, blacs_context, LLD_C, info)
      IF (info /= 0) THEN 
         WRITE(*,*) 'commy proc ', parallel%commy_myid,': DESCINIT matC failed, info=',info; GOTO 999 
      ENDIF
      !
      CALL PDGEMM(opA_up, opB_up, M, N, K, alpha, &
                  matA, 1_I4B, 1_I4B, desc_matA, &
                  matB, 1_I4B, 1_I4B, desc_matB, &
                  beta, matC, 1_I4B, 1_I4B, desc_matC)
      !
      IF (PRESENT(matC_glob)) THEN
         IF (SIZE(matC_glob,1) /= M_C .OR. SIZE(matC_glob,2) /= N_C) THEN
               info = -1; WRITE(*,*) 'matC_glob dim mismatch: expect (',M_C,',',N_C,') got (',SIZE(matC_glob,1),',',SIZE(matC_glob,2),')'; GOTO 999
         ENDIF

         CALL DESCINIT(desc_matC_glob, M_C, N_C, M_C, N_C, 0, 0, blacs_context, M_C, info)
         IF (info /= 0) THEN
            WRITE(*,*) 'DESCINIT desc_matC_glob failed at commy root, info=',info
            GOTO 999
         ENDIF

         CALL PDGEMR2D(M_C, N_C, &
                        matC, 1, 1, desc_matC, &
                        matC_glob, 1, 1, desc_matC_glob, &
                        blacs_context)

         CALL MPI_BCAST(matC_glob, M_C*N_C, MPI_REAL8, &
                        0, parallel%commy, ierr)
      ENDIF
   999 CONTINUE
   IF (info /= 0) THEN
      CALL BLACS_GRIDEXIT(blacs_context)
      WRITE(*, '(A,I0,A,I0)') 'proc ', parallel%myid, ': lapk_real_matmat failed, final info=', info
   ENDIF
   ENDSUBROUTINE scalapk_real_matmat
!----------------------------------------------------------

!----------------------------------------------------------
   ELEMENTAL FUNCTION UPCASE(c) RESULT(res)
      CHARACTER(1), INTENT(IN) :: c
      CHARACTER(1) :: res
      IF (c >= 'a' .AND. c <= 'z') THEN
         res = ACHAR(IACHAR(c) - 32)
      ELSE
         res = c
      ENDIF
   ENDFUNCTION UPCASE
!----------------------------------------------------------
END MODULE ScaLapack_module
#endif