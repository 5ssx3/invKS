MODULE invKS_module
  USE constants
  USE parameters , ONLY : wtol=>etol,gtol=>Rtol &
             & ,iopt,nskip
  USE struct_module, ONLY : lat=>lat_mat 
#ifdef MPI
  USE smpi_math_module
#endif
  IMPLICIT NONE
  CONTAINS
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !##################################################!
  !           SCF for screen potential               !
  !##################################################!
  !-------------------optim_potential-----------------
  SUBROUTINE iks_optim_WY()
     USE parameters , ONLY : nev,nssp,Idiag, Chetol,nspin,sigma=>Wsmear
     USE grid_module , ONLY : nr1,nr2,nr3,nr,nk,kpt,nrs,eigen,grid,eigen_type, destroy_eigen
     !USE scf_module , ONLY : smear_updaterho,kfilter,ksolver
     USE smearing_module , ONLY : smear_init,Fermilevel,smear_updaterho &
               & ,wke,fme,ets
     USE struct_module , ONLY : ne_r,dvol,volume
      USE chebyshev_module, ONLY : BuildSubspace,CheF_all
     USE opt_module
#ifdef MPI
       USE smpi_math_module
#endif
     IMPLICIT NONE
     !IN/OUT
     !LOCAL
     REAL(DP) :: nwsp,penfun,fmed
     INTEGER(I4B) :: Iter,m=10,Ifilter
     !d
     REAL(DP) :: drho,gNorm,dt=0.3_dp
     !xqtest 
     INTEGER(I4B) :: Is
     REAL(DP) :: ngvec(nr)
     TYPE(eigen_type) :: eig
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     !max iterations for WY functional
     SELECT CASE(Iopt)
     CASE(1)
        CALL initCG(nr)
     CASE(2)
        CALL initLBFGS(nr,m)
     CASE default
        CALL initFIRE(nr,dt)
     END SELECT
     !init
     CALL smear_init(nev)
     !first eigepairs
      DO is=1,nspin
         grid%veff(:,is)=grid%vef0(:)
      ENDDO

      CALL BuildSubspace(nr,nev,grid%veff,eigen)

      CALL smear_updaterho(nr,nev,ne_r,eigen,grid%rhoS,grid%rho)

#ifdef MPI
      IF (ALLOCATED(recvcounts)) DEALLOCATE(recvcounts)
      IF (ALLOCATED(displs)) DEALLOCATE(displs)
      ALLOCATE(recvcounts(parallel%commy_numprocs))
      ALLOCATE(displs(parallel%commy_numprocs))
      IF (parallel%isroot) THEN
#endif
      PRINT*,'Total electrons testing:',SUM(grid%rho)*dvol
#ifdef MPI
   ENDIF
#endif
     !PRINT*,eigen%val(:,1,1)
     maxWY : DO Iter=1,nssp
        !WY functional
        CALL iks_Ws(nr,grid%rho,grid%rho0,grid%veff(:,1) &
           & ,eigen%val,wke,ets,nwsp,penfun,ngvec)
        !opt
        SELECT CASE(Iopt)
        CASE(1)
           WRITE(*,*)'================== PERFORM CGplus =================='
           CALL CGplus_optm(Iter,nr,grid%veff(:,1),nwsp,ngvec)
           !PRINT*,ICALL,'Ws=',-nWs,'penalty=',penfun
           !WRITE(*,*)'================= CHECK THE FLAG ================='
           IF(IFLAG<=0.OR.ICALL>=10000)THEN
              STOP
           ELSE IF(IFLAG==1)THEN
              ICALL=ICALL+1
              !PRINT*,'Continue search'
           ENDIF
        CASE(2)
           WRITE(*,*)'================== PERFORM BFGS =================='
           myicount=myicount+1
           CALL LBFGS_optm(Iter,nr,m,grid%veff(:,1),nwsp,ngvec)
           print*,'LBFGS task :',task
           !============================================
           WRITE(*,*)'================= CHECK THE FLAG ================='
           IF(task(1:5)=='NEW_X')THEN
              nline=0
              IF(myicount>maxcount)THEN
                 PRINT*,'Maxcount BFGS is reached'
                 OPEN(1111,FILE='Veff.dat')
                     WRITE(1111,*) 'Ne',ne_r
                     WRITE(1111,*) nr1,nr2,nr3
                     WRITE(1111,*) grid%veff(:,1)*hart2eV
                 CLOSE(1111)
              ENDIF
              !store current Ws
              histW(1)=histW(2)
              histW(2)=histW(3)
              histW(3)=histW(4)
              histW(4)=nwsp
              IF(myicount>5)THEN
                  IF(abs(histW(4)-histW(3))<wtol .AND. &
                  &  abs(histW(3)-histW(2))<wtol .AND. &
                  &  abs(histW(2)-histW(1))<wtol       &
                  !&  gNorm<gtol) THEN
                  &  ) THEN
                     PRINT*,'Ws has reach the accury:',wtol*hart2eV,' eV'
                     OPEN(1111,FILE='Veff.dat')
                         WRITE(1111,*) 'Ne',ne_r
                         WRITE(1111,*) nr1,nr2,nr3
                         WRITE(1111,*) grid%veff(:,1)*hart2eV
                     CLOSE(1111)
                     !difference
                     OPEN(200,FILE='RhoF.dat')
                         WRITE(200,*) nr1,nr2,nr3
                         WRITE(200,*) grid%rho(:)*volume
                     CLOSE(200)

                     IF(gNorm<gtol)THEN
                        EXIT maxWY
                     ELSE
                        task='START'
                     ENDIF
                  ENDIF
              ENDIF
           ELSEIF(task(1:2)=='FG')THEN
              IF(nline>0)THEN
                 PRINT*,'TASK=FG, continue line search'
              ENDIF
              nline=nline+1
           ELSE
              PRINT*,'L-BFGS task warning , RE-START anyway,task=>'
              PRINT*,task
              task='START'
              Chetol=Chetol*0.9d0
           ENDIF
        CASE default
#ifdef MPI
         CALL MPI_BARRIER(parallel%comm,mpinfo)
         IF (parallel%isroot) THEN
#endif
           WRITE(*,*)'================== PERFORM FIRE =================='
#ifdef MPI
         ENDIF
#endif
            CALL FIRE_optm(nr,grid%veff(:,1),nwsp,ngvec)
            ! PRINT*,'process:',parallel%myid,'Iter=',Iter,'grid%Veff(1:5,1)=', grid%veff(1:5,1)
           myicount=myicount+1
           histW(1)=histW(2)
           histW(2)=histW(3)
           histW(3)=histW(4)
           histW(4)=nwsp
           IF(myicount>5)THEN
               IF(abs(histW(4)-histW(3))<wtol .AND. &
               &  abs(histW(3)-histW(2))<wtol .AND. &
               &  abs(histW(2)-histW(1))<wtol       &
               &   )THEN
#ifdef MPI
         IF (parallel%isroot) THEN
#endif
                  PRINT*,'Ws has reach the accury:',wtol*hart2ev,' eV'
                  OPEN(1111,FILE='Veff.dat')
                      WRITE(1111,*) 'Ne',ne_r
                      WRITE(1111,*) nr1,nr2,nr3
                      WRITE(1111,*) grid%veff(:,:)*hart2eV
                  CLOSE(1111)
                  !difference
                  OPEN(200,FILE='RhoF.dat')
                      WRITE(200,*) nr1,nr2,nr3
                      WRITE(200,*) grid%rho(:)*volume
                  CLOSE(200)
#ifdef MPI
         ENDIF
#endif            
                  IF(gNorm<gtol)  EXIT maxWY
               ENDIF
           ENDIF
           !============================================
        END SELECT
        !gnorm
        gNorm=SUM(ABS(ngvec))/nr
        drho=SUM(ABS(grid%rho-grid%rho0))/nr
#ifdef MPI
      IF (parallel%isroot) THEN
#endif
         PRINT*,Iter,'Ws(eV)=',-nwsp*hart2ev,'penalty(eV)=',penfun*hart2ev
         PRINT*,'AverNorm/min/max(density difference) in a.u.'
         PRINT*,drho,MINVAL(grid%rho-grid%rho0),MAXVAL(grid%rho-grid%rho0)
         PRINT*,'AverNorm/min/max(nGvec) in a.u.'
         PRINT*,gnorm,MINVAL(ngvec),MAXVAL(ngvec)
         PRINT*,'Total rho difference should be 0:',SUM(grid%rho-grid%rho0)*dvol
#ifdef MPI
      ENDIF
      CALL MPI_BARRIER(parallel%comm,mpinfo)
#endif
        !Getting next eiegnpairs by CheFSI
        DO Ifilter=1,10
           fmed=fme
            ! PRINT*,'process:',parallel%myid,'Ifilter=',Ifilter
            CALL CheF_all(nr,nev,grid%veff,eigen)
            !Find Fermi Level
#ifdef MPI
            ! PRINT*, 'process',parallel%myid,'Calling Fermilevel in iks_optim_WY'
            CALL Fermilevel(ne_r,nev,parallel%mygrid_range(3),kpt%wk,eigen%val,sigma)
#else
            CALL Fermilevel(ne_r,nev,nk,kpt%wk,eigen%val,sigma)
#endif           
           IF(ABS(fme-fmed)<Chetol) EXIT
        ENDDO
        !filter acucury
#ifdef MPI
   IF (parallel%isroot) THEN
#endif
        PRINT*,MIN(Ifilter,10),'Chebyshev Accury:',ABS(fme-fmed)*hart2ev,'eV'
#ifdef MPI
   ENDIF
#endif

      CALL smear_updaterho(nr,nev,ne_r,eigen,grid%rhoS,grid%rho)

     ENDDO maxWY
     !
#ifdef MPI
      IF (parallel%isroot) THEN
#endif
         PRINT*,'Ws has reach the accury:',wtol*hart2ev,' eV'
         OPEN(1111,FILE='Veff.dat')
            WRITE(1111,*) 'Ne',ne_r
            WRITE(1111,*) nr1,nr2,nr3
            WRITE(1111,*) grid%veff(:,:)*hart2eV
         CLOSE(1111)
         !difference
         OPEN(200,FILE='RhoF.dat')
            WRITE(200,*) nr1,nr2,nr3
            WRITE(200,*) grid%rho(:)*volume
         CLOSE(200)
#ifdef MPI
   ENDIF
   CALL MPI_Barrier(parallel%comm, mpinfo)
   IF (ALLOCATED(parallel%recvcounts)) DEALLOCATE(parallel%recvcounts)
   IF (ALLOCATED(parallel%displs)) DEALLOCATE(parallel%displs)
   IF (ALLOCATED(parallel%global_gridrange)) DEALLOCATE(parallel%global_gridrange)
   IF (ALLOCATED(displs)) DEALLOCATE(displs)
   IF (ALLOCATED(recvcounts)) DEALLOCATE(recvcounts)
   CALL destroy_eigen()
#endif
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE iks_optim_WY
!-------------------------PARTING LINE-------------------------
  SUBROUTINE iks_Ws(nps,rho,rho0,vs,evals,wks,etpy,nwsp,penfun,ngvec)
     USE parameters , ONLY : penLambda
     USE finite_module , ONLY : real_pbc_nabla2,real_pbc_nabla1
     USE grid_module , ONLY : nr
     USE struct_module, ONLY: dvol
     IMPLICIT NONE
     !IN/OUT
     INTEGER(I4B),INTENT(IN) :: nps
     REAL(DP),INTENT(IN) :: rho(nps),rho0(nps),vs(nps)
     REAL(DP),DIMENSION(:,:,:),INTENT(IN) :: &
           &  evals  & !eigen-vaule
           &, wks    !weight
     REAL(DP),INTENT(IN) :: etpy !entropy
     REAL(DP),INTENT(OUT) :: nwsp,penfun !-WY functional=-Ws
     REAL(DP),INTENT(OUT) :: ngvec(nps)
     !LOCAL
     REAL(DP) :: dvs(3,nps),d2vs(nps) &!penalty function
             &, Ws
#ifdef MPI
     REAL(DP) :: Ws_local
#endif
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     !Penalty function Lambda*SUM(|\nabla V|^2)*dvol
      CALL real_pbc_nabla1(vs,dvs)
      CALL real_pbc_nabla2(vs,d2vs)
      penfun= penLambda*SUM(dvs(:,:)*dvs(:,:))*dvol
      !Ws=Es-TS+SUM( Vs(rho-rho0) ) *dvol= Eband-TS-SUM(vs*rho0)*dvol
#ifdef MPI
      Ws_local=SUM(wks*evals)
      CALL MPI_ALLREDUCE(Ws_local,Ws,1,MPI_REAL8,MPI_SUM,parallel%commx,mpinfo)
      Ws=Ws - etpy - SUM(vs*rho0)*dvol
#else
      Ws=SUM(wks*evals) - etpy - SUM(vs*rho0)*dvol
#endif
      !L=-Ws + PF
      nwsp= -Ws + penfun
      !gradient of WY functional with penalty function
      ngvec(:)=  rho0(:) - rho(:)  - 2.d0*penLambda*d2vs(:)


     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE iks_Ws

!-------------------------PARTING LINE-------------------------
ENDMODULE invKS_module
