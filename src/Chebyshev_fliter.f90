!#############################################################!
!        Chebyshev filter Method for diag(H)                  !
!*Author : Qiang Xu                                           !
!*Date   : 2017/11/28                                         !
!#############################################################!
MODULE chebyshev_module
   USE constants
   USE parameters , ONLY : CheM, nev
#ifdef MPI
   USE smpi_math_module
   USE ScaLapack_module
#endif
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
      INTEGER(I4B) :: seed_size
      INTEGER(I4B) :: seed(2)
#ifdef MPI
      INTEGER(I4B) :: ierr , i
      INTEGER(I4B) :: Ik_global
      INTEGER(I4B) :: Ik_local
      INTEGER(I4B) :: recvcount
      INTEGER(I4B),ALLOCATABLE :: recevcounts(:),displs(:)
      REAL(DP),ALLOCATABLE :: X0(:,:)
      REAL(DP),ALLOCATABLE :: X0_global(:,:)
#else
      REAL(DP) :: X0(nps,nev)
#endif
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#ifdef MPI 
      IF (ALLOCATED(X0)) DEALLOCATE(X0)
      IF (ALLOCATED(X0_global)) DEALLOCATE(X0_global)
      IF (ALLOCATED(recevcounts)) DEALLOCATE(recevcounts, STAT=ierr)
      IF (ALLOCATED(displs)) DEALLOCATE(displs, STAT=ierr)
      IF (parallel%isroot) THEN
         ALLOCATE(X0_global(nps,nev))
         ALLOCATE(X0(nps,parallel%nstate_proc))
      ELSE
         ALLOCATE(X0(nps,parallel%nstate_proc))
      ENDIF
      IF (parallel%isroot) THEN
#endif
         PRINT*,'Building Cheby-Subspace for this task ...'
#ifdef MPI
      ENDIF
#endif
      Nrand=nev
      irand1=1
      !random states
#ifdef MPI
   IF (parallel%isroot) THEN
         IF(Nrand>0)THEN
         !initialize the temp subspace
         seed = [2147483562, 2147483398]
         CALL random_seed(put=seed)
         print *, "put seed: ", seed
         DO Ii=irand1,nev
            DO Ip=1,nps
               CALL random_number(randt)
               !CALL random_number(randr)
               randt = -radius+randt*2.d0*radius
               !randr = -radius+randr*10.d0*radius
               X0_global(Ip,Ii) = randt !(randt+randr)/2
            ENDDO
         ENDDO
      ENDIF
   ENDIF
#else
      IF(Nrand>0)THEN
         !initialize the temp subspace
         CALL random_seed(size=seed_size)
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
#endif
#ifdef MPI
      ALLOCATE(recevcounts(parallel%commy_numprocs))
      ALLOCATE(displs(parallel%commy_numprocs))
      CALL MPI_ALLGATHER(parallel%nstate_proc,1,MPI_INTEGER4,recevcounts,1,MPI_INTEGER4,parallel%commy,mpinfo)
      recevcounts = nps*recevcounts
      displs(1) = 0
      DO i= 2,parallel%commy_numprocs
         displs(i) = recevcounts(i-1) + displs(i-1)
      ENDDO
      CALL MPI_Barrier(parallel%commy, mpinfo)
      recvcount = nps * parallel%nstate_proc

      IF (parallel%isroot) THEN
         CALL MPI_Scatterv(X0_global,recevcounts,displs,MPI_REAL8,X0,recvcount,MPI_REAL8, 0, parallel%commy, mpinfo)
      ELSE
         CALL MPI_Scatterv(MPI_BOTTOM, [0], [0], MPI_REAL8, X0, recvcount, MPI_REAL8, 0, parallel%commy, mpinfo)
      ENDIF
      CALL MPI_BCAST(X0, parallel%nstate_proc*nps, MPI_REAL8, 0, parallel%commx, mpinfo)
      IF (ALLOCATED(recevcounts)) DEALLOCATE(recevcounts)
      IF (ALLOCATED(displs)) DEALLOCATE(displs)
      IF (parallel%isroot) THEN
         IF (ALLOCATED(X0_global)) DEALLOCATE(X0_global)
      ENDIF
      CALL Init_scala()
      CALL MPI_Barrier(parallel%comm, mpinfo)
#endif
      !Raleigh-Ritz step
      !Gamma point
#ifdef MPI
      DO Is=1,Nspin !spiner

         DO Ik_local=1,parallel%mygrid_range(3) !k-points
            Ik_global=parallel%global_gridrange(1,parallel%commx_myid+1)+Ik_local-1
            PRINT*
            IF(Ik_global /= IGamma)THEN
               !for non-Gamma k-points
               eig%wvf(:,:,Ik_local,Is)=X0(:,:)
               IF(CheM0>0)THEN
                   !filter
                  CALL cmplx_first_filter(nps,parallel%nstate_proc,Ik_global,veff(:,Is), &
                     & eig%wvf(:,:,Ik_local,Is),eig%val(:,Ik_local,Is))
               ELSE
                   !RR
                  CALL cmplx_first_RRstep(nps,parallel%nstate_proc,Ik_global,veff(:,Is), &
                      &   eig%wvf(:,:,Ik_local,Is),eig%val(:,Ik_local,Is))
               ENDIF
            ELSE
            !for Gamma k-points
               eig%wvfG(:,:,Is)=X0(:,:)
               IF(CheM0>0)THEN
                   !filter
                   CALL real_first_filter(nps,parallel%nstate_proc,veff(:,Is), &
                     & eig%wvfG(:,:,Is),eig%val(:,Ik_local,Is))
               ELSE
                   !RR
                   CALL real_first_RRstep(nps,parallel%nstate_proc,veff(:,Is), &
                      &   eig%wvfG(:,:,Is),eig%val(:,Ik_local,Is))
               ENDIF
            ENDIF
         ENDDO
      ENDDO
#else
      DO Is=1,Nspin !spiner

         DO Ik=1,nk !k-points
            !print*,'Ik',Ik
            IF(Ik /= IGamma)THEN
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
#endif
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE BuildSubspace
   !----------------------------------------------------------
   SUBROUTINE CheF_all(nps,nev,veff,eig)
      !need to be improve for spin
      USE parameters , ONLY : nspin,IGamma
      USE grid_module , ONLY : nk,eigen_type
#ifdef MPI
      USE smpi_math_module
#endif
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN) :: nps,nev
      REAL(DP),INTENT(IN) :: veff(nps,nspin)
      TYPE(eigen_type),INTENT(INOUT) :: eig
      !LOCAL
      INTEGER(I4B) :: Is,Ik
#ifdef MPI
      INTEGER(I4B) :: Ik_global
#endif
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! PRINT*,'process:',parallel%myid,'Entering CheF_all'
      DO Is=1,nspin
#ifdef MPI
      DO Ik=1,parallel%mygrid_range(3)
         Ik_global=parallel%global_gridrange(1,parallel%commx_myid+1)+Ik-1
         IF(Ik_global==IGamma)THEN
            CALL real_filtering(nps,parallel%nstate_proc,veff(:,Is),eig%wvfG(:,:,Is),eig%val(:,Ik,Is))
         ELSE
            CALL cmplx_filtering(nps,parallel%nstate_proc,Ik_global,veff(:,Is),eig%wvf(:,:,Ik,Is),eig%val(:,Ik,Is))
         ENDIF
      ENDDO
#else
      DO Ik=1,nk
         IF(Ik==IGamma)THEN
            CALL real_filtering(nps,nev,veff(:,Is),eig%wvfG(:,:,Is),eig%val(:,Ik,Is))
         ELSE
            CALL cmplx_filtering(nps,nev,Ik,veff(:,Is),eig%wvf(:,:,Ik,Is),eig%val(:,Ik,Is))
         ENDIF
      ENDDO
#endif
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE CheF_all
   !----------------------------------------------------------
   SUBROUTINE real_first_RRstep(nps,nst,veff,X,D)
      USE Lapack_module , ONLY : lapk_OrthNorm, lapk_Eigen, lapk_MM
      USE parameters , ONLY : CheM0,LRROrthNorm,nev
      IMPLICIT NONE
      !IN/OUT
      INTEGER(I4B),INTENT(IN) :: nps,nst
      REAL(DP),INTENT(IN) :: veff(nps)
#ifdef MPI
      REAL(DP),INTENT(INOUT) :: X(nps,nst) !subspace
      REAL(DP),INTENT(OUT) :: D(:)      !eigenvalue
      !LOCAL
      REAL(DP) :: Xnew(nps,nst)
      !REAL(DP) :: a,b,al
      REAL(DP),DIMENSION(nev,nev) :: Shat,Hhat,Qs
      REAL(DP) :: Shat_local(nev,nst), Qs_local(nev,nst)
      INTEGER(I4B),ALLOCATABLE :: displs(:), recevcounts(:),i,s_start,s_end
#else
      REAL(DP),INTENT(INOUT) :: X(nps,nev) !subspace
      REAL(DP),INTENT(OUT) :: D(:)      !eigenvalue
      !LOCAL
      REAL(DP) :: Xnew(nps,nev)
      !REAL(DP) :: a,b,al
      REAL(DP),DIMENSION(nev,nev) :: Shat,Hhat,Qs
#endif
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! PRINT*,'process:',parallel%myid,'Entering real_first_RRstep'
      !Raleigh-Ritz step
      IF(LRROrthNorm)THEN
         !OrthNorm
         CALL lapk_OrthNorm(X)
         !RR
         CALL real_Rayleigh_quotient(nps,nst,veff,X,Hhat)
         !eigen-decomposion
         CALL lapk_Eigen(Hhat,Qs,D)
      ELSE
         !Overlap matrix
#ifdef MPI
         CALL lapk_MM(X,X,'T','N',1.0_dp,0.0_dp,Shat_local)
         IF (ALLOCATED(displs)) DEALLOCATE(displs)
         IF (ALLOCATED(recevcounts)) DEALLOCATE(recevcounts)
         ALLOCATE(displs(parallel%commy_numprocs))
         ALLOCATE(recevcounts(parallel%commy_numprocs))
         CALL MPI_ALLGATHER(parallel%nstate_proc,1,MPI_INTEGER4, &
                     recevcounts,1,MPI_INTEGER4,parallel%commy,mpinfo)
         displs(1) = 0
         recevcounts = nev * recevcounts
         DO i = 2,parallel%commy_numprocs
            displs(i) = displs(i-1) +recevcounts(i-1)
         ENDDO
         CALL MPI_ALLGATHERV(Shat_local, nev*nst, MPI_REAL8, &
                             Shat,recevcounts,displs,MPI_REAL8,parallel%commy,mpinfo)
#else
         CALL lapk_MM(X,X,'T','N',1.0_dp,0.0_dp,Shat)
#endif
         !projected hamiltonian
         CALL real_Rayleigh_quotient(nps,nst,veff,X,Hhat)
         !eigen-decomposion
         CALL lapk_Eigen(nev,Hhat,Shat,Qs,D)
      ENDIF
      !-------------------
      !rotation
#ifdef MPI
      s_start = parallel%sub2sum(1,parallel%commy_myid+1)
      s_end = parallel%sub2sum(parallel%nstate_proc,parallel%commy_myid+1)
      Qs_local = Qs(:,s_start:s_end)
      CALL lapk_MM(X,Qs_local,'N','N',1.0_dp,0.0_dp,Xnew)
#else
      CALL lapk_MM(X,Qs,'N','N',1.0_dp,0.0_dp,Xnew)
#endif
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
      ! PRINT*,'process:',parallel%myid,'Entering real_init_uplow'
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
      USE Lapack_module , ONLY : lapk_Eigen, lapk_OrthNorm, lapk_MM
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
      REAL(DP),DIMENSION(nev,nev)  :: Hhat,Qs,Shat
      REAL(DP)  :: Xnew(nps,nst),vec(nps)
#ifdef MPI
      INTEGER(I4B) :: BCAST_ID_local, BCAST_ID
      REAL(DP),DIMENSION(nev,nst) :: Shat_local, Qs_local
      INTEGER(I4B),ALLOCATABLE :: displs(:), recevcounts(:)
      INTEGER(I4B) :: j,s_start,s_end
#endif
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !low up bound
#ifdef MPI
      ! PRINT*,'process:',parallel%myid,'Entering real_first_filter'
      vec(:)=X(:,nst)
      BCAST_ID_local = -1
      BCAST_ID = -1
      IF (parallel%sub2sum(parallel%nstate_proc,parallel%commy_myid+1) == nev) THEN
         BCAST_ID_local = parallel%commy_myid
         CALL real_init_uplow(nps,7,veff,vec,a,b,al)
      ENDIF
      CALL MPI_ALLREDUCE(BCAST_ID_local,BCAST_ID,1,MPI_INTEGER4,MPI_MAX,parallel%commy,mpinfo)
      CALL MPI_BCAST(a,1,MPI_REAL8,BCAST_ID,parallel%commy,mpinfo)
      CALL MPI_BCAST(al,1,MPI_REAL8,BCAST_ID,parallel%commy,mpinfo)
      CALL MPI_BCAST(b,1,MPI_REAL8,BCAST_ID,parallel%commy,mpinfo)
      ! PRINT*,'process',parallel%myid,'compeleted MPI_BCAST in real_first_filter'
#else
      vec(:)=X(:,nst)
      CALL real_init_uplow(nps,7,veff,vec,a,b,al)
#endif
      !print*,'a,b,al',a,b,al
      !
      evald(:)=612509.d0
      DO I=1,Niter
         !
         CALL real_chebyshev_filter_scaled(nps,nst,veff,X,CF0,a,b,al)
         ! PRINT*,'process',parallel%myid,'compeleted real_chebyshev_filter_scaled in real_first_filter'
         IF(LRROrthNorm)THEN
            !OrthNorm
            CALL lapk_OrthNorm(X)
            !xHx
            CALL real_Rayleigh_quotient(nps,nst,veff,X,Hhat)
            !eigen-decomposion
            CALL lapk_Eigen(Hhat,Qs,eval)
         ELSE
            !Overlap matrix
#ifdef MPI
            ! PRINT*,'process',parallel%myid,'started to call scalapk_MM(X,X) in real_first_filter'
            CALL scalapk_MM(X,X,'T','N',1.0_dp,0.0_dp,Shat_local)
            ! PRINT*,'process',parallel%myid,'compeleted scalapk_MM(X,X) in real_first_filter'
            IF (ALLOCATED(displs)) DEALLOCATE(displs)
            IF (ALLOCATED(recevcounts)) DEALLOCATE(recevcounts)
            ALLOCATE(displs(parallel%commy_numprocs))
            ALLOCATE(recevcounts(parallel%commy_numprocs))
            CALL MPI_ALLGATHER(parallel%nstate_proc,1,MPI_INTEGER4, &
                        recevcounts,1,MPI_INTEGER4,parallel%commy,mpinfo)
            displs(1) = 0
            recevcounts = nev * recevcounts
            DO j = 2,parallel%commy_numprocs
               displs(j) = displs(j-1) +recevcounts(j-1)
            ENDDO
            CALL MPI_ALLGATHERV(Shat_local, nev*nst, MPI_REAL8, &
                              Shat,recevcounts,displs,MPI_REAL8,parallel%commy,mpinfo)
            ! PRINT*,'process',parallel%myid,'compeleted MPI_ALLGATHERV in real_first_filter'
#else
            CALL lapk_MM(X,X,'T','N',1.0_dp,0.0_dp,Shat)
#endif
            !projected hamiltonian
            CALL real_Rayleigh_quotient(nps,nst,veff,X,Hhat)
            !eigen-decomposion
#ifdef MPI
            CALL lapk_Eigen(nev,Hhat,Shat,Qs,eval)
            ! PRINT*,'process',parallel%myid,'compeleted lapk_Eigen in real_first_filter'
#else
            CALL lapk_Eigen(nst,Hhat,Shat,Qs,eval)
#endif
         ENDIF
         deval=SUM(ABS(eval-evald))/nst
         !-----------------
         IF(deval<TOL) EXIT
         !-----------------
         !store old eigenvalue
         evald(:)=eval(:)
         !update the new bound
         a=eval(nev)
         al=eval(1)
      ENDDO
      !rotation
#ifdef MPI
      s_start = parallel%sub2sum(1,parallel%commy_myid+1)
      s_end = parallel%sub2sum(parallel%nstate_proc,parallel%commy_myid+1)
      Qs_local = Qs(:,s_start:s_end)
      CALL scalapk_MM(X,Qs_local,'N','N',1.0_dp,0.0_dp,Xnew)
      ! PRINT*,'process',parallel%myid,'compeleted lapk_MM(X,Qs_local) in real_first_filter'
#else
      CALL lapk_MM(X,Qs,'N','N',1.0_dp,0.0_dp,Xnew)
#endif
      X=Xnew
      ! PRINT*,'process:',parallel%myid,'completed X=Xnew in real_first_filter'
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
! #ifdef MPI
!       PRINT*,'process:',parallel%myid,'Entering cmplx_first_RRstep'
! #endif
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
! #ifdef MPI
!       PRINT*,'process:',parallel%myid,'Entering cmplx_init_uplow'
! #endif
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
      USE parameters , ONLY : CF0=>CheM0, LRROrthNorm, nev
      USE Lapack_module , ONLY : lapk_Eigen, lapk_OrthNorm, lapk_MM
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN) :: nps,nst,Ik !number of points/states
      REAL(DP),INTENT(IN) :: veff(nps)
      COMPLEX(DCP),INTENT(INOUT) :: X(:,:)
      REAL(DP),INTENT(OUT) :: eval(:)
      !LOCAL
      REAL(DP) :: a,b,al,t
      INTEGER(I4B) :: I,Niter=4
      REAL(DP) :: deval,TOL=1e-8
      REAL(DP) :: evald(nev)
      COMPLEX(DCP),DIMENSION(nev,nev)  :: Hhat,Qs,Shat
      COMPLEX(DCP)  :: Xnew(nps,nst),vec(nps)
#ifdef MPI
      INTEGER(I4B) :: BCAST_ID_local, BCAST_ID
      COMPLEX(DCP),DIMENSION(nev,nst) :: Shat_local, Qs_local
      INTEGER(I4B),ALLOCATABLE :: displs(:), recevcounts(:)
      INTEGER(I4B) :: j, s_start, s_end
#endif
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !low up bound
#ifdef MPI
      ! PRINT*,'process:',parallel%myid,'Entering cmplx_first_filter'
      ! PRINT*,'process:',parallel%myid,'size(X)=',SIZE(X,1),SIZE(X,2)
      vec(:)=X(:,nst)
      IF (parallel%commy_numprocs > 1) THEN
         BCAST_ID_local = -1
         BCAST_ID = -1
         IF (parallel%sub2sum(parallel%nstate_proc,parallel%commy_myid+1) == nev) THEN
            BCAST_ID_local = parallel%commy_myid
            CALL cmplx_init_uplow(nps,7,Ik,veff,vec,a,b,al)
         ENDIF
         CALL MPI_ALLREDUCE(BCAST_ID_local,BCAST_ID,1,MPI_INTEGER4,MPI_MAX,parallel%commy,mpinfo)
         CALL MPI_BCAST(a,1,MPI_REAL8,BCAST_ID,parallel%commy,mpinfo)
         CALL MPI_BCAST(al,1,MPI_REAL8,BCAST_ID,parallel%commy,mpinfo)
         CALL MPI_BCAST(b,1,MPI_REAL8,BCAST_ID,parallel%commy,mpinfo)
         ! PRINT*,'process',parallel%myid,'completed MPI_BCAST in cmplx_first_filter'
      ELSE
         CALL cmplx_init_uplow(nps,7,Ik,veff,vec,a,b,al)
      ENDIF
#else
      vec(:)=X(:,nst)
      CALL cmplx_init_uplow(nps,7,Ik,veff,vec,a,b,al)
#endif
      !
      evald(:)=612509.d0
      DO I=1,Niter
         !
         CALL cmplx_chebyshev_filter_scaled(nps,nst,Ik,veff,X,CF0,a,b,al)
         ! PRINT*,'process',parallel%myid,'completed cmplx_chebyshev_filter_scaled in cmplx_first_filter'
         IF(LRROrthNorm)THEN
            !OrthNorm
            CALL lapk_OrthNorm(X)
            !xHx
            CALL cmplx_Rayleigh_quotient(nps,nst,Ik,veff,X,Hhat)
            !eigen-decomposion
            CALL lapk_Eigen(Hhat,Qs,eval)
         ELSE
#ifdef MPI
            ! PRINT*,'process',parallel%myid,'started to call scalapk_MM(X,X) in cmplx_first_filter'
            CALL scalapk_MM(X, X,'C','N', cmplx(1._dp,0._dp,DCP), cmplx(0._dp,0._dp,DCP), Shat_local)
            ! PRINT*,'process',parallel%myid,'completed scalapk_MM(X,X) in cmplx_first_filter'
            IF (ALLOCATED(displs)) DEALLOCATE(displs)
            IF (ALLOCATED(recevcounts)) DEALLOCATE(recevcounts)
            ALLOCATE(displs(parallel%commy_numprocs))
            ALLOCATE(recevcounts(parallel%commy_numprocs))
            CALL MPI_ALLGATHER(parallel%nstate_proc,1,MPI_INTEGER4, &
                        recevcounts,1,MPI_INTEGER4,parallel%commy,mpinfo)
            displs(1) = 0
            recevcounts = nev * recevcounts
            DO j = 2,parallel%commy_numprocs
               displs(j) = displs(j-1) +recevcounts(j-1)
            ENDDO
            CALL MPI_ALLGATHERV(Shat_local, nev*nst, MPI_DOUBLE_COMPLEX, &
                              Shat,recevcounts,displs,MPI_DOUBLE_COMPLEX,parallel%commy,mpinfo)
            ! PRINT*,'process',parallel%myid,'compeleted MPI_ALLGATHERV in cmplx_first_filter'
#else
            CALL lapk_MM(X,X,'T','N',cmplx(1._dp,0._dp,DCP),cmplx(0._dp,0._dp,DCP),Shat)
#endif
            !projected hamiltonian
            CALL cmplx_Rayleigh_quotient(nps,nst,Ik,veff,X,Hhat)
            ! PRINT*,'process',parallel%myid,'completed cmplx_Rayleigh_quotient in cmplx_first_filter'
            !eigen-decomposion
#ifdef MPI
            CALL lapk_Eigen(nev,Hhat,Shat,Qs,eval)
            ! PRINT*,'process',parallel%myid,'completed lapk_Eigen in cmplx_first_filter'
#else
            CALL lapk_Eigen(nst,Hhat,Shat,Qs,eval)
#endif
         ENDIF
         deval=SUM(ABS(eval-evald))/nst
         !-----------------
         IF(deval<TOL) EXIT
         !-----------------
         !store old eigenvalue
         evald(:)=eval(:)
         !update the new bound
         a=eval(nev)
         al=eval(1)
      ENDDO
      !rotation
#ifdef MPI
      s_start = parallel%sub2sum(1,parallel%commy_myid+1)
      s_end = parallel%sub2sum(parallel%nstate_proc,parallel%commy_myid+1)
      Qs_local = Qs(:,s_start:s_end)
      CALL scalapk_MM(X,Qs_local,'N','N',cmplx(1._dp,0._dp,DCP),cmplx(0._dp,0._dp,DCP),Xnew)
      ! PRINT*,'process',parallel%myid,'completed scalapk_MM(X,Qs_local) in cmplx_first_filter'
#else
      CALL lapk_MM(X,Qs,'N','N',cmplx(1._dp,0._dp,DCP),cmplx(0._dp,0._dp,DCP),Xnew)
#endif
      X=Xnew
      !
      ! PRINT*,'process:',parallel%myid,'completed X=Xnew cmplx_first_filter'
      !
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
      REAL(DP),INTENT(OUT) :: xhx(nev,nev)
      !LOCAL
      REAL(DP) :: hx(nps,nst)
#ifdef MPI
      INTEGER(I4B) :: q
      INTEGER(I4B),ALLOCATABLE :: displs(:), recevcounts(:)
      REAL(DP) :: xhx_temp(nst,nst), xhx_local(nev,nst),V_q(nps,nst)
#endif
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! PRINT*,'process:',parallel%myid,'Entering real_Rayleigh_quotient'
      CALL real_HX(nps,nst,veff,x,hx)
      ! PRINT*,'process:',parallel%myid,'completed real_HX in real_Rayleigh_quotient'
#ifdef MPI
      xhx_local = 0.0_DP
      xhx_temp = 0.0_DP
         DO q = 0,parallel%commy_numprocs-1
            IF (q == parallel%commy_myid) THEN
               V_q = x
            ENDIF
            call MPI_Bcast(V_q, nps*nst, MPI_REAL8, q, parallel%commy, mpinfo)
            IF (mpinfo /= 0) THEN
               IF (ALLOCATED(displs)) DEALLOCATE(displs)
               IF (ALLOCATED(recevcounts)) DEALLOCATE(recevcounts)
               WRITE(*,*) 'MPI_Bcast failed at q=',q,'process ID=',parallel%myid
               CALL MPI_Abort(parallel%commy, 1, mpinfo)
            ENDIF
            CALL lapk_MM(V_q, hx, 'T', 'N', 1._dp, 0._dp, xhx_temp)
            ! PRINT*,'process:',parallel%myid,'completed lapk_MM(V_q, hx) in real_Rayleigh_quotient'
            xhx_local(q*nst+1:(q+1)*nst, :) = xhx_temp
            ! PRINT*,'process:',parallel%myid,'completed xhx_local(q*nst+1:(q+1)*nst, :) = xhx_temp in real_Rayleigh_quotient, q=',q
         ENDDO
      IF (ALLOCATED(displs)) DEALLOCATE(displs)
      IF (ALLOCATED(recevcounts)) DEALLOCATE(recevcounts)
      ALLOCATE(displs(parallel%commy_numprocs))
      ALLOCATE(recevcounts(parallel%commy_numprocs))
      CALL MPI_ALLGATHER(parallel%nstate_proc,1,MPI_INTEGER4,recevcounts,1,MPI_INTEGER4,parallel%commy,mpinfo)
      recevcounts = nev * recevcounts
      displs(1) = 0
      DO q = 2,parallel%commy_numprocs
         displs(q) = recevcounts(q-1) + displs(q-1)
      ENDDO
      CALL MPI_ALLGATHERV(xhx_local,nst*nev,MPI_REAL8,xhx,recevcounts,displs,MPI_REAL8,parallel%commy,mpinfo)
      ! PRINT*,'process:',parallel%myid,'completed MPI_ALLGATHERV(xhx) in real_Rayleigh_quotient'
      DEALLOCATE(displs, recevcounts)
#else
      CALL lapk_MM(x,hx,'T','N',1._dp,0._dp,xhx)
#endif
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE real_Rayleigh_quotient
   !-------------------------------HX------------------------------
   SUBROUTINE real_HX(nps,nst,veff,V,HV)
      USE matvec_module , ONLY : matvec
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN)    :: nps,nst
      REAL(DP),INTENT(IN)        :: veff(nps), V(nps,nst)
      REAL(DP),INTENT(OUT)       :: HV(nps,nst)
      INTEGER(I4B)               :: Is, ierr
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      DO Is=1,nst,1
         CALL matvec(veff,V(:,Is),HV(:,Is),nps)
      ENDDO
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
   SUBROUTINE real_GRayleigh_Ritz(nps,nst,veff,X,D)
      USE Lapack_module , ONLY : lapk_Eigen,lapk_MM
      IMPLICIT NONE
      !IN/OUT
      INTEGER(I4B),INTENT(IN) :: nps,nst
      REAL(DP),INTENT(IN) :: veff(nps)
      REAL(DP),INTENT(INOUT) :: X(nps,nst)
      REAL(DP),INTENT(OUT) :: D(:)
      !
      REAL(DP),DIMENSION(nev,nev) :: S_hat,H_hat,Q
      REAL(DP) :: Xnew(nps,nst)
#ifdef MPI
      REAL(DP) :: S_hat_local(nev,nst), Q_local(nev,nst)
      INTEGER(I4B) :: j, s_end, s_start
#endif
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !Calculate the overlap matrix
#ifdef MPI
      CALL scalapk_MM(X,X,'T','N',1._dp,0._dp,S_hat_local)
      CALL MPI_ALLGATHER(parallel%nstate_proc,1,MPI_INTEGER4, &
                        recvcounts,1,MPI_INTEGER4,parallel%commy,mpinfo)
      displs(1) = 0
      recvcounts = nev * recvcounts
      DO j = 2,parallel%commy_numprocs
         displs(j) = displs(j-1) +recvcounts(j-1)
      ENDDO
      CALL MPI_ALLGATHERV(S_hat_local, nev*nst, MPI_REAL8, &
                        S_hat,recvcounts,displs,MPI_REAL8,parallel%commy,mpinfo)
#else
      CALL lapk_MM(X,X,'T','N',1._dp,0._dp,S_hat)
#endif
      !Calculate the project hamiltion
      CALL real_Rayleigh_quotient(nps,nst,veff,X,H_hat)
      !solve the generalized eigenvalue problem
      CALL lapk_Eigen(nev,H_hat,S_hat,Q,D)
      !X=XQ
#ifdef MPI
      s_start = parallel%sub2sum(1,parallel%commy_myid+1)
      s_end = parallel%sub2sum(parallel%nstate_proc,parallel%commy_myid+1)
      Q_local = Q(:,s_start:s_end)
      CALL scalapk_MM(X,Q_local,'N','N',1._dp,0._dp,Xnew)
#else
      CALL lapk_MM(X,Q,'N','N',1._dp,0._dp,Xnew)
#endif
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
   SUBROUTINE real_filtering(nps,nst,veff,X,D)
      USE parameters , ONLY : LRROrthNorm, nev
      !To avoid OrthNorm
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN) :: nps,nst
      REAL(DP),INTENT(IN) :: veff(nps)
      REAL(DP),INTENT(INOUT) :: X(nps,nst)
      REAL(DP),INTENT(INOUT) :: D(:)  !rayleigh-ritz value
      !LOCAL
      REAL(DP) :: a,b,al
#ifdef MPI
      INTEGER(I4B) :: BCAST_ID_local, BCAST_ID
#endif
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      a=MAXVAL(D)
      al=MINVAL(D)
      !up boundary
#ifdef MPI
      BCAST_ID_local = -1
      BCAST_ID = -1
      b=0.0_DP
      IF (parallel%sub2sum(parallel%nstate_proc,parallel%commy_myid+1) == nev) THEN
         CALL real_Estupb(nps,7,veff,X(:,nst),b)
         BCAST_ID_local = parallel%commy_myid
      ENDIF
      CALL MPI_ALLREDUCE(BCAST_ID_local,BCAST_ID,1,MPI_INTEGER4,MPI_MAX,parallel%commy,mpinfo)
      IF (mpinfo /= 0) THEN
         PRINT*,'ERROR: process',parallel%myid, 'MPI_ALLREDUCE failed! mpinfo=',mpinfo
         PRINT*,'  parallel%commy=',parallel%commy
         CALL MPI_Abort(parallel%commy, 1, mpinfo)
      ENDIF
      CALL MPI_BCAST(b,1,MPI_REAL8,BCAST_ID,parallel%commy,mpinfo)
#else
      CALL real_Estupb(nps,7,veff,X(:,nst),b)
#endif
      !filtering (a,b) 
      CALL real_chebyshev_filter_scaled(nps,nst,veff,X,CheM,a,b,al)
      !CALL chebyshev_filter_real(nps,nev,veff,X,CheM,a,b)
      IF(LRROrthNorm)THEN
         !RR (Rayleigh-Ritz Step)
         CALL real_Rayleigh_Ritz(nps,nst,veff,X,D)
      ELSE
         !GRR (Rayleigh-Ritz Step)
         CALL real_GRayleigh_Ritz(nps,nst,veff,X,D)
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
      COMPLEX(DCP),INTENT(OUT) :: xhx(nev,nev)
      !LOCAL
      COMPLEX(DCP) :: hx(nps,nst)
#ifdef MPI
      INTEGER(I4B) :: q
      INTEGER(I4B),ALLOCATABLE :: displs(:), recevcounts(:)
      COMPLEX(DCP) :: xhx_temp(nst,nst), xhx_local(nev,nst),V_q(nps,nst)
#endif
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CALL cmplx_HX(nps,nst,Ik,veff,x,hx)
#ifdef MPI
      xhx_local = (0.0_DP, 0.0_DP)
      xhx_temp = (0.0_DP, 0.0_DP)
      DO q = 0,parallel%commy_numprocs-1
         IF (q == parallel%commy_myid) THEN
            V_q = x
         ENDIF
         call MPI_Bcast(V_q, nps*nst, MPI_DOUBLE_COMPLEX, q, parallel%commy, mpinfo)
         IF (mpinfo /= 0) THEN
            IF (ALLOCATED(displs)) DEALLOCATE(displs)
            IF (ALLOCATED(recevcounts)) DEALLOCATE(recevcounts)
            WRITE(*,*) 'MPI_Bcast failed at q=',q,'process ID=',parallel%myid
            CALL MPI_Abort(parallel%commy, 1, mpinfo)
         ENDIF
         CALL lapk_MM(V_q, hx, 'C', 'N', cmplx(1.0,0.0,DCP),cmplx(0._dp,0._dp,DCP), xhx_temp)
         xhx_local(q*nst+1:(q+1)*nst, :) = xhx_temp
      ENDDO
      IF (ALLOCATED(displs)) DEALLOCATE(displs)
      IF (ALLOCATED(recevcounts)) DEALLOCATE(recevcounts)
      ALLOCATE(displs(parallel%commy_numprocs))
      ALLOCATE(recevcounts(parallel%commy_numprocs))
      CALL MPI_ALLGATHER(parallel%nstate_proc,1,MPI_INTEGER4,recevcounts,1,MPI_INTEGER4,parallel%commy,mpinfo)
      recevcounts = nev * recevcounts
      displs(1) = 0
      DO q = 2,parallel%commy_numprocs
         displs(q) = recevcounts(q-1) + displs(q-1)
      ENDDO
      CALL MPI_ALLGATHERV(xhx_local,nst*nev,MPI_DOUBLE_COMPLEX,xhx,recevcounts,displs,MPI_DOUBLE_COMPLEX,parallel%commy,mpinfo)
      DEALLOCATE(displs, recevcounts)
#else
      CALL lapk_MM(x,hx,'C','N',cmplx(1.0,0.0,DCP),cmplx(0._dp,0._dp,DCP),xhx)
#endif
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
   SUBROUTINE cmplx_GRayleigh_Ritz(nps,nst,Ik,veff,X,D)
      USE Lapack_module , ONLY : lapk_Eigen, lapk_MM
      IMPLICIT NONE
      !IN/OUT
      INTEGER(I4B),INTENT(IN) :: nps,nst,Ik
      REAL(DP),INTENT(IN) :: veff(nps)
      COMPLEX(DCP),INTENT(INOUT) :: X(nps,nst)
      REAL(DP),INTENT(OUT) :: D(:)
      !
      COMPLEX(DCP),DIMENSION(nev,nev) :: S_hat,H_hat,Q
      COMPLEX(DCP) :: Xnew(nps,nst)
#ifdef MPI
      COMPLEX(DCP) :: S_hat_local(nev,nst), Q_local(nev,nst)
      INTEGER(I4B) :: j, s_end, s_start
#endif
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !Calculate the overlap matrix
#ifdef MPI
      CALL scalapk_MM(X,X,'C','N',cmplx(1._dp,0._dp,DCP),cmplx(0._dp,0._dp,DCP),S_hat_local)
      CALL MPI_ALLGATHER(parallel%nstate_proc,1,MPI_INTEGER4, &
                        recvcounts,1,MPI_INTEGER4,parallel%commy,mpinfo)
      displs(1) = 0
      recvcounts = nev * recvcounts
      DO j = 2,parallel%commy_numprocs
         displs(j) = recvcounts(j-1) + displs(j-1)
      ENDDO
      CALL MPI_ALLGATHERV(S_hat_local, nev*nst, MPI_DOUBLE_COMPLEX, &
                        S_hat,recvcounts,displs,MPI_DOUBLE_COMPLEX,parallel%commy,mpinfo)
#else
      CALL lapk_MM(X,X,'T','N',cmplx(1._dp,0._dp,DCP),cmplx(0._dp,0._dp,DCP),S_hat)
#endif
      !Calculate the project hamiltion
      CALL cmplx_Rayleigh_quotient(nps,nst,Ik,veff,X,H_hat)
      !solve the generalized eigenvalue problem
      CALL lapk_Eigen(nev,H_hat,S_hat,Q,D)
      !X=XQ
#ifdef MPI
      s_start = parallel%sub2sum(1,parallel%commy_myid+1)
      s_end = parallel%sub2sum(parallel%nstate_proc,parallel%commy_myid+1)
      Q_local = Q(:,s_start:s_end)
      CALL scalapk_MM(X,Q_local,'N','N',cmplx(1._dp,0._dp,DCP),cmplx(0._dp,0._dp,DCP),Xnew)
#else 
      CALL lapk_MM(X,Q,'N','N',cmplx(1._dp,0._dp,DCP),cmplx(0._dp,0._dp,DCP),Xnew)
#endif
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
#ifdef MPI
      COMPLEX(DP) :: Qs_local(nev,sn)
      INTEGER(I4B) :: s_end, s_start
#endif
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !F0:OrthNorm
#ifdef MPI
      CALL scalapk_OrthNorm(X)
#else
      CALL lapk_OrthNorm(X)
#endif
      !F1:Rayleigh_quotient
      CALL cmplx_Rayleigh_quotient(nps,sn,Ik,veff,X,Hhat)
      !F2:eigen-decomposition Q,D
      CALL lapk_Eigen(Hhat,Qs,D)
      !F3:X=XQ
      !X=MATMUL( X , Q )
#ifdef MPI
      s_start = parallel%sub2sum(1,parallel%commy_myid+1)
      s_end = parallel%sub2sum(parallel%nstate_proc,parallel%commy_myid+1)
      Qs_local = Qs(:,s_start:s_end)
      CALL scalapk_MM(X,Qs_local,'N','N',cmplx(1._dp,0._dp,DCP),cmplx(0._dp,0._dp,DCP),Xnew)
#else
      CALL lapk_MM(X,Qs,'N','N',cmplx(1._dp,0._dp,DCP),cmplx(0._dp,0._dp,DCP),Xnew)
#endif
      X=Xnew
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE cmplx_Rayleigh_Ritz
   !-----------------non-OrthNorm Chebyshev_filter ----------------
   SUBROUTINE cmplx_filtering(nps,nst,Ik,veff,X,D)
      USE parameters , ONLY : LRROrthNorm
      !To avoid OrthNorm
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN) :: nps,nst,Ik
      REAL(DP),INTENT(IN) :: veff(nps)
      COMPLEX(DCP),INTENT(INOUT) :: X(nps,nst)
      REAL(DP),INTENT(INOUT) :: D(:)  !rayleigh-ritz value
      !LOCAL
      REAL(DP) :: a,b,al
#ifdef MPI
      INTEGER(I4B) :: BCAST_ID_local, BCAST_ID
#endif
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      a=MAXVAL(D)
      al=MINVAL(D)
      !up boundary
#ifdef MPI
      BCAST_ID_local = -1
      BCAST_ID = -1
      IF (parallel%sub2sum(parallel%nstate_proc,parallel%commy_myid+1) == nev) THEN
         CALL cmplx_Estupb(nps,7,Ik,veff,X(:,nst),b)
         BCAST_ID_local = parallel%commy_myid
      ENDIF
      CALL MPI_ALLREDUCE(BCAST_ID_local,BCAST_ID,1,MPI_INTEGER4,MPI_MAX,parallel%commy,mpinfo)
      CALL MPI_BCAST(b,1,MPI_REAL8,BCAST_ID,parallel%commy,mpinfo)
#else
      CALL cmplx_Estupb(nps,7,Ik,veff,X(:,nst),b)
#endif
      !filtering (a,b) 
      CALL cmplx_chebyshev_filter_scaled(nps,nst,Ik,veff,X,CheM,a,b,al)
      !CALL chebyshev_filter_real(nps,nev,veff,X,CheM,a,b)
      IF(LRROrthNorm)THEN
         !RR (Rayleigh-Ritz Step)
         CALL cmplx_Rayleigh_Ritz(nps,nst,Ik,veff,X,D)
      ELSE
         !GRR (Rayleigh-Ritz Step)
         CALL cmplx_GRayleigh_Ritz(nps,nst,Ik,veff,X,D)
      ENDIF
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE cmplx_filtering
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDMODULE chebyshev_module
