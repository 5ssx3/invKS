MODULE struct_module
   !!####################################################!!
   !!*************  author: Qiang Xu          ***********!!
   !!*************  date  : 2017-07-11        ***********!!
   !!*************  Save Structure Information        ***!!
   !!####################################################!!
   USE constants
   IMPLICIT NONE
   !-------------------------------------------------------------
   TYPE struct_type
      INTEGER(I4B),ALLOCATABLE  :: nati(:)  !number of atoms for each atom type
      REAL(DP),ALLOCATABLE      :: posdir(:,:)  !direct
      REAL(DP),ALLOCATABLE      :: poscar(:,:) ! car
   END TYPE struct_type
   INTEGER(I4B)              :: natom! total number of atoms
   INTEGER(I4B)              :: naty ! number of atom types 
   INTEGER(I4B),ALLOCATABLE  :: nati ! number of atoms for each types 
   INTEGER(I4B)              :: ne_i    ! total change in integer
   REAL(DP)                  :: ne_r ! real type charge
   REAL(DP)                  :: volume,dvol
   REAL(DP)                  :: lat_mat(3,3)    
   REAL(DP)                  :: lat_para(6)
   REAL(DP)                  :: recip_lat(3,3) !!! the lattice matrix in reciprocal space 
   REAL(DP)                  :: reclat_para(6)
   REAL(DP)                  :: energy(10)      ! energy
   !cell symmetry
   REAL(DP)                  :: Opsym(3,3,48) !operator of symmetry
   REAL(DP)                  :: Otrans(3,48)
   INTEGER(I4B)              :: nsym  !number of symmetry we used
   INTEGER(I4B)              :: num_t  !number of translations we used
   INTEGER(I4B)              :: c_i(8)
   INTEGER(I4B)              :: Odet(48)
   TYPE(struct_type) :: struct
CONTAINS
   !--------------------------------------------------------------
   SUBROUTINE creat_struct(elen)
      IMPLICIT NONE
      !
      INTEGER(I4B),INTENT(IN) :: elen(naty)
      !TYPE(struct_type),INTENT(OUT) :: struct
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      natom=sum(elen(:))
      CALL destroy_struct()
      ALLOCATE(struct%nati(naty))
      ALLOCATE(struct%posdir(3,natom))
      ALLOCATE(struct%poscar(3,natom))
      struct%nati(:)=elen(:)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE creat_struct
   !--------------------------
   SUBROUTINE destroy_struct()
      IMPLICIT NONE
      !
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF(ALLOCATED(struct%nati))      DEALLOCATE(struct%nati)
      IF(ALLOCATED(struct%posdir))       DEALLOCATE(struct%posdir)
      IF(ALLOCATED(struct%poscar))    DEALLOCATE(struct%poscar)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE destroy_struct
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDMODULE struct_module


MODULE grid_module
  !##############################################!
  !*For    : real space grid mesh and relate mesh!
  !*Author : Qiang Xu                            !
  !*Date   : 2017-7-12                           !
  !##############################################!
  USE constants
  IMPLICIT NONE
  !defined
  TYPE grid_type
     !3D real space grid mesh
     REAL(DP),ALLOCATABLE :: rho0(:),vef0(:)
     REAL(DP),ALLOCATABLE :: rhoS(:,:)
     REAL(DP),ALLOCATABLE :: rho(:)
     REAL(DP),ALLOCATABLE :: veff(:,:)
     !recip g table
     REAL(DP),ALLOCATABLE :: gvec(:,:)
     LOGICAL ,ALLOCATABLE :: gMask(:)
  ENDTYPE grid_type
  !k-points data
  TYPE kgrid_type
     REAL(DP),ALLOCATABLE :: vdir(:,:)
     REAL(DP),ALLOCATABLE :: vcar(:,:)
     REAL(DP),ALLOCATABLE :: wk(:)
  ENDTYPE kgrid_type
  !eigen data
  TYPE eigen_type
     !eigen-states
     COMPLEX(DCP),ALLOCATABLE :: wvf(:,:,:,:)
     !eigen-values
     REAL(DP),ALLOCATABLE :: val(:,:,:)
     REAL(DP),ALLOCATABLE :: wvfG(:,:,:)
  ENDTYPE eigen_type
  !------------------Basic data--------------------
  INTEGER(I4B) :: nr1,nr2,nr3,nr,nrs
  INTEGER(I4B) :: ng1,ng2,ng3,ng
  !grid info
  REAL(DP) :: gap(3)=0.2
  !-------------------kgrids-----------------------
  INTEGER(I4B) :: nk1,nk2,nk3,nk
  !-------------------states-----------------------
  REAL(DP)     :: kdispl(3)
  !basic grid mesh
  TYPE(grid_type)  ::  grid
  !band structrue kpoints
  TYPE(kgrid_type) ::  kpt
  !eigen states and value(k-represent)
  TYPE(eigen_type) :: eigen
CONTAINS
  !-------------------PARTING LINE--------------------
  SUBROUTINE read_data()
     USE parameters, ONLY : nskip,nspin,nev,nadds
     USE struct_module, ONLY : volume,dvol,ne_i,ne_r
#ifdef MPI
      USE smpi_math_module, ONLY:parallel
#endif
     IMPLICIT NONE
     INTEGER(I4B)  :: fs1=0,fs2=0
     INTEGER(I4B) :: ip,nx,ny,nz
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#ifdef MPI
      IF (parallel%isroot) THEN
#endif
     PRINT*,'Reading CHGCAR and LOCPOT files ...'
#ifdef MPI
      ENDIF
#endif
     OPEN(110,FILE='CHGCAR',IOSTAT=fs1)
          IF(fs1/=0)THEN
             WRITE(*,*) 'STOP: Could not open',fs1
             STOP
          ENDIF
     OPEN(111,FILE='LOCPOT',IOSTAT=fs2)
          IF(fs2/=0)THEN
             WRITE(*,*) 'STOP: Could not open',fs2
             STOP
          ENDIF
#ifdef MPI
   IF (parallel%isroot) THEN
#endif
        PRINT*,'nskip', nskip
#ifdef MPI
   ENDIF
#endif
        DO ip =1,nskip
           READ(110,*)
           READ(111,*)
        ENDDO
        READ(110,*) nr1,nr2,nr3
        READ(111,*) nx,ny,nz
        IF(nr1/=nx .OR. nr2/=ny .OR. nr3/=nz)THEN
             STOP "Check the grids in CHGCAR and LOCPOT"
        ENDIF
#ifdef MPI
   IF (parallel%isroot) THEN
#endif
        PRINT*,'Grids(nx,ny,nz):',nr1,nr2,nr3
#ifdef MPI
   ENDIF
#endif
        nr=nr1*nr2*nr3
        nrs=nr*nspin
        ALLOCATE(grid%rho0(nr),grid%vef0(nr))
        READ(110,*) grid%rho0(:)
        READ(111,*) grid%vef0(:)
     CLOSE(110)
     CLOSE(111)
     !convert unit to a.u.
     grid%rho0=grid%rho0/volume
     grid%vef0=grid%vef0/hart2ev
     dvol=volume/nr
     !FIND CHARGES
     ne_r=SUM(grid%rho0)*dvol
     ne_i=ANINT(ne_r)
     ne_r=REAL(ne_i,DP)
#ifdef MPI
   IF (parallel%isroot) THEN
#endif
      PRINT*,'Total # of electrons are:', ne_r
#ifdef MPI
   ENDIF
#endif
     nev=ne_i/2+MAX(2,nadds)
#ifdef MPI
   IF (parallel%isroot) THEN
#endif
      PRINT*,'# of states would be calculated', nev
#ifdef MPI
   ENDIF
#endif
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE read_data
  !-------------------PARTING LINE--------------------
  SUBROUTINE Build_rgrid()
     USE parameters , ONLY : nspin
     USE struct_module, ONLY : lat_mat,lat_para &
               &, reclat_para,recip_lat
     IMPLICIT NONE
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     CALL read_data()
     !grid size
     gap(1)=lat_para(1)/nr1
     gap(2)=lat_para(2)/nr2
     gap(3)=lat_para(3)/nr3
     !total mesh points in real space
     nr=nr1*nr2*nr3
     nrs=nr*nspin
     !
     ng1=nr1/2+1
     ng2=nr2
     ng3=nr3
     ng=ng1*ng2*ng3
     !
     CALL destroy_rgrid()
     !3D
     ALLOCATE(grid%rho(nr))
     ALLOCATE(grid%rhoS(nr,nspin))
     ALLOCATE(grid%veff(nr,nspin)) 
     !recip grids
     ALLOCATE(grid%gvec(4,ng))
     ALLOCATE(grid%gMask(ng))
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE Build_rgrid
  !-------------------PARTING LINE--------------------
  SUBROUTINE destroy_rgrid()
     IMPLICIT NONE
     !
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     !3D
     IF(ALLOCATED(grid%rhoS))  DEALLOCATE(grid%rhoS)
     IF(ALLOCATED(grid%rho))  DEALLOCATE(grid%rho)
     IF(ALLOCATED(grid%veff))  DEALLOCATE(grid%veff)
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE destroy_rgrid
  !-------------------PARTING LINE--------------------
  SUBROUTINE Build_kgrid()
     USE parameters , ONLY : kspacing,kgrid,IGamma
     USE struct_module , ONLY : recip_lat,lat_para
#ifdef MPI
      USE smpi_math_module, ONLY : parallel
#endif
     IMPLICIT NONE
     LOGICAL  :: lkmesh
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     !PRINT*,'Building k-points grid ...'
     CALL destroy_kpt()
     lkmesh=((kgrid(1)>0).AND.(kgrid(2)>0).AND.(kgrid(3)>0))
     !case for kspacing work
     IF((kspacing>0.d0).OR.lkmesh)THEN
        !symmetry operator
        CALL Symm_kgrid()
     ELSE
        !only gamma point
        IGamma=1
        nk1=1
        nk2=1
        nk3=1
        nk=1
        !>>>>>>>>>>>>>>>>>>>>>>>>>>
        ALLOCATE(kpt%vdir(3,nk))
        ALLOCATE(kpt%vcar(3,nk))
        ALLOCATE(kpt%wk(nk))
        !<<<<<<<<<<<<<<<<<<<<<<<<<<
        kpt%wk(:)=1
        kpt%vdir(:,:)=0.d0
        kpt%vcar(:,:)=0.d0
     ENDIF
#ifdef MPI
   IF (parallel%isroot) THEN
#endif
     PRINT*,'k-point meshes:',nk1,nk2,nk3
#ifdef MPI
   ENDIF
#endif
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE Build_kgrid
  !-------------------PARTING LINE--------------------
  SUBROUTINE destroy_kpt()
     IMPLICIT NONE
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     IF(ALLOCATED(kpt%wk))    DEALLOCATE(kpt%wk)
     IF(ALLOCATED(kpt%vdir))  DEALLOCATE(kpt%vdir)
     IF(ALLOCATED(kpt%vcar))  DEALLOCATE(kpt%vcar)
     !IF(ALLOCATED(KPT%nksym))  DEALLOCATE(KPT%nksym)
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE destroy_kpt
  !-------------------PARTING LINE--------------------
  SUBROUTINE Build_eigen()
      USE parameters , ONLY : nspin,nev,IGamma
#ifdef MPI
      USE smpi_math_module, ONLY : grid_split, array_split, parallel
#endif
     IMPLICIT NONE
#ifdef MPI
      INTEGER(I4B) :: Ik, Ik_global, mpinfo
      LOGICAL :: has_gamma
#endif
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#ifdef MPI
   IF (.NOT. ALLOCATED(parallel%recvcounts)) ALLOCATE(parallel%recvcounts(parallel%commx_numprocs))
   IF (.NOT. ALLOCATED(parallel%displs)) ALLOCATE(parallel%displs(parallel%commx_numprocs))
   IF (.NOT. ALLOCATED(parallel%global_gridrange)) ALLOCATE(parallel%global_gridrange(3,parallel%commx_numprocs))
   CALL grid_split(nk,parallel%commx_numprocs,parallel%commx,parallel%commx_myid,parallel%mygrid_range,parallel%recvcounts,parallel%displs,parallel%global_gridrange)
      ! Print the grid_split results (for each process)
      
   CALL array_split(nev, parallel%commy_numprocs, parallel%commy, parallel%commy_myid, parallel%nstate_proc, parallel%sub2sum)
      ! Print the results of array_split (for each process)
#endif
     CALL destroy_eigen()
     IF(IGamma<=0)THEN
#ifdef MPI
         ALLOCATE(eigen%wvf(nr,parallel%nstate_proc,parallel%mygrid_range(3),Nspin))
#else
         ALLOCATE(eigen%wvf(nr,nev,nk,Nspin))
#endif
     ELSE !have gamma point
         IF(nk>=2)THEN

#ifdef MPI
         ALLOCATE(eigen%wvf(nr,parallel%nstate_proc,parallel%mygrid_range(3),Nspin))
#else
         ALLOCATE(eigen%wvf(nr,nev,nk-1,Nspin))
#endif
         ENDIF
        !
#ifdef MPI
        has_gamma = .FALSE.
        DO Ik = 1, parallel%mygrid_range(3)
            Ik_global = parallel%displs(parallel%commx_myid+1) + Ik
            IF (Ik_global == IGamma) THEN
                has_gamma = .TRUE.
                EXIT
            ENDIF
        ENDDO
         IF (has_gamma) THEN
            ALLOCATE(eigen%wvfG(nr, parallel%nstate_proc, Nspin))
         ENDIF
#else
         ALLOCATE(eigen%wvf(nr,nev,nk-1,Nspin))
         ALLOCATE(eigen%wvfG(nr,parallel%nstate_proc,Nspin))
#endif
     ENDIF
     !eigen value
#ifdef MPI
      ALLOCATE(eigen%val(nev,parallel%mygrid_range(3),nspin))
      DEALLOCATE(parallel%recvcounts)
      DEALLOCATE(parallel%displs)
#else
      ALLOCATE(eigen%val(nev,nk,nspin))
#endif
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE Build_eigen
  !-------------------PARTING LINE--------------------
  SUBROUTINE destroy_eigen()
     IMPLICIT NONE
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     IF(ALLOCATED(eigen%wvf)) DEALLOCATE(eigen%wvf)
     IF(ALLOCATED(eigen%wvfG)) DEALLOCATE(eigen%wvfG)
     IF(ALLOCATED(eigen%val)) DEALLOCATE(eigen%val)
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE destroy_eigen
  !---------------------------------------------------
  SUBROUTINE Symm_kgrid()
     USE math , ONLY : thr2mat,norm
     USE parameters , ONLY : kspacing,Isym,kgrid,IGamma
     USE struct_module , ONLY : Opsym,nsym,recip_lat,reclat_para
     IMPLICIT NONE
     INTEGER(I4B) :: nkr  !total k-points
     REAL(DP),ALLOCATABLE :: xk(:,:),wk(:),wk0(:),kvec(:,:)
     INTEGER(I4B),ALLOCATABLE :: equiv(:)
     REAL(DP) :: xkr(3),xx,yy,zz,fact
     !
     INTEGER(I4B) :: k1,k2,k3,Ik,Is,I,J,K,Il,offsetk(3)
     INTEGER(I4B) :: nk0
     LOGICAL :: time_reversal,linlist
     REAL(DP),PARAMETER :: eps=1D-5
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     IF(Isym>=0)THEN
        time_reversal=.TRUE.
     ELSE
        time_reversal=.FALSE.
     ENDIF
     !kgrid
     IF((kgrid(1)>0).AND.(kgrid(2)>0).AND.(kgrid(3)>0))THEN
       nk1=kgrid(1)
       nk2=kgrid(2)
       nk3=kgrid(3)
     ELSE
       nk1=MAX(1,CEILING(reclat_para(1)/kspacing))
       nk2=MAX(1,CEILING(reclat_para(2)/kspacing))
       nk3=MAX(1,CEILING(reclat_para(3)/kspacing))
     ENDIF

     nkr=nk1*nk2*nk3
     ALLOCATE(xk(3,nkr),kvec(3,nkr),wk(nkr),wk0(nkr))
     ALLOCATE(equiv(nkr))
     !generate the k-points by Monkhorst and Pack sampling
     Ik=0
     DO k3=1,nk3
     DO k2=1,nk2
     DO k1=1,nk1
        Ik=Ik+1 
        xk(1,Ik)=(2.0_DP*k1-nk1-1)/2.0_DP/nk1
        xk(2,Ik)=(2.0_DP*k2-nk2-1)/2.0_DP/nk2
        xk(3,Ik)=(2.0_DP*k3-nk3-1)/2.0_DP/nk3
     ENDDO
     ENDDO
     ENDDO
     offsetk(1)=1-nk1
     offsetk(2)=1-nk2
     offsetk(3)=1-nk3
     !find the Irreducible Brillouin Zon(IBZ) k-points
     DO Ik=1,nkr
        equiv(Ik)=Ik
     ENDDO
     !
     DO Ik=1,nkr
    !  check if this k-point has already been found equivalent to another
        IF(equiv(Ik)==Ik)THEN
           wk(Ik)=1.0_DP
        !  check if there are equivalent k-point to this in the list
        !  (excepted those previously found to be equivalent to another)
        !  check both k and -k
           DO Is=1,nsym
              xkr(:)=MATMUL(TRANSPOSE(Opsym(:,:,Is)),xk(:,Ik))
              xkr(:)=xkr(:)-NINT(xkr(:))
              !IF(Odet(Ik)==-1)
              xx=xkr(1)*nk1-0.5_DP*offsetk(1)
              yy=xkr(2)*nk2-0.5_DP*offsetk(2)
              zz=xkr(3)*nk3-0.5_DP*offsetk(3)
              linlist= ABS(xx-NINT(xx))<=eps .AND. &
                    &  ABS(yy-NINT(yy))<=eps .AND. &
                    &  ABS(zz-NINT(zz))<=eps
              IF(linlist)THEN
                 I=MOD(NINT(xkr(1)*nk1-0.5_DP*offsetk(1)+2*nk1),nk1)+1
                 J=MOD(NINT(xkr(2)*nk2-0.5_DP*offsetk(2)+2*nk2),nk2)+1
                 K=MOD(NINT(xkr(3)*nk3-0.5_DP*offsetk(3)+2*nk3),nk3)+1
                 !find the index
                 CALL thr2mat(nk1,nk2,nk3,I,J,K,Il)
                 !havn't been used?
                 IF(Il>Ik.AND.equiv(Il)==Il)THEN
                    equiv(Il)=Ik
                    wk(Ik)=wk(Ik)+1._DP
                 ELSE

                    IF(equiv(Il)/=Ik.OR.Il<Ik)THEN
                       WRITE(*,*) 'Symmetry k-points: Some thing wrong 1'
                       STOP
                    ENDIF

                 ENDIF

              ENDIF
              !time reversal
              IF(time_reversal)THEN
                 xx=-xkr(1)*nk1-0.5_DP*offsetk(1)
                 yy=-xkr(2)*nk2-0.5_DP*offsetk(2)
                 zz=-xkr(3)*nk3-0.5_DP*offsetk(3)
                 linlist=ABS(xx-NINT(xx))<=eps .AND. &
                     &   ABS(yy-NINT(yy))<=eps .AND. &
                     &   ABS(zz-NINT(zz))<=eps
                 IF(linlist)THEN
                    I=MOD(NINT(-xkr(1)*nk1-0.5_DP*offsetk(1)+2*nk1),nk1)+1
                    J=MOD(NINT(-xkr(2)*nk2-0.5_DP*offsetk(2)+2*nk2),nk2)+1
                    K=MOD(NINT(-xkr(3)*nk3-0.5_DP*offsetk(3)+2*nk3),nk3)+1
                    CALL thr2mat(nk1,nk2,nk3,I,J,K,Il)
                    IF(Il>Ik.AND.equiv(Il)==Il)THEN
                       equiv(Il)=Ik
                       wk(Ik)=wk(Ik)+1._DP
                    ELSE

                       IF(equiv(Il)/=Ik.OR.Il<Ik)THEN
                          WRITE(*,*) 'Symmetry k-points: Some thing wrong 2'
                          STOP
                       ENDIF

                    ENDIF

                 ENDIF

              ENDIF
           !all symmetry operator
           ENDDO

        ENDIF
     !cycle all k
     ENDDO
     !
     nk0=0
     fact=0._DP
     wk0(:)=0._DP
     kvec(:,:)=0._DP
     DO Ik=1,nkr
        IF(equiv(Ik)==Ik)THEN
          nk0=nk0+1
          IF(nk0>nkr) STOP 'Symmetry k-points: Some thing wrong 3'
          wk0(nk0)=wk(Ik)
          fact=fact+wk0(nk0)
          !store kvec
          kvec(:,nk0)=xk(:,Ik)
        ENDIF
     ENDDO
     !normalize
     wk0(1:nk0)=wk0(1:nk0)/fact
     !build data
     nk=nk0
     !-----------------------
     ALLOCATE(kpt%vdir(3,nk))
     ALLOCATE(kpt%vcar(3,nk))
     ALLOCATE(kpt%wk(nk))
     !ALLOCATE(KPT%nksym(nk))
     !-----------------------
     !print*,'kmesh used',sum(KPT%wk(:)),fact
     kpt%vdir(:,1:nk)=kvec(:,1:nk)
     kpt%wk(1:nk)=wk0(1:nk)
     DO Ik=1,nk
        kpt%vcar(:,Ik)=MATMUL(recip_lat,kvec(:,Ik))
        !print*,KPT%vec(:,Ik),KPT%wk(Ik)
        IF(norm(kpt%vcar(:,Ik))<xtiny)THEN
           IGamma=Ik
        ENDIF
     ENDDO
     !reset
     IF(IGamma>0)THEN
        !noGamma
        I=0
        DO Ik=1,nk
           IF(Ik/=IGamma)THEN
               I=I+1
               kpt%vdir(:,I)=kvec(:,Ik)
               kpt%vcar(:,I)=MATMUL(recip_lat,kvec(:,Ik))
               kpt%wk(I)=wk0(Ik)
           ENDIF
        ENDDO
        !Gamma
        kpt%vdir(:,nk)=kvec(:,IGamma)
        kpt%vcar(:,nk)=MATMUL(recip_lat,kvec(:,IGamma))
        kpt%wk(nk)=wk0(IGamma)
        !reset Gamma
        IGamma=nk
     ENDIF
     !deallocate something
     DEALLOCATE(xk,kvec,wk,wk0,equiv)
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE Symm_kgrid
  !---------------------------------------------------
    SUBROUTINE FillQTable()
     !##########################################################
     !* CREATED_TIME  : 2015-05-08
     !* AUTHOR        : Yuan Zhou
     !* CHANGE        : Xuecheng Shao
     !* ADD           :
     !* DESCRIPTION   :
     !     This routine updates the table of norms that tells energy calculators in
     !     reciprocal space whether to include a given q-point in the sum or leave it
     !     out (based on whether it is within the energyCutoff sphere in recip. sp.)
     !     It should be called every time the cell dimensions are altered.
     !* REFERENCES    :
     !     ------
     !* LOG           :
     !     2015-05-08 :
     !* LAST MODIFIED : 2015-05-08 08:04:24 PM
     !##########################################################
     USE struct_module , ONLY : recip_lat
     USE math , ONLY : Norm
     IMPLICIT NONE
     !
     INTEGER(I4B) :: ix,iy,iz,Ig
     REAL(DP)     :: mVector(3),qPoint(3)
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     Ig=0
     DO iz = 1, ng3
        mVector(3) = iz - 1
        ! Wrap around if index is greater than 1/2 the table size
        IF (mVector(3)>ng3/2) mVector(3)=mVector(3)-ng3
     DO iy = 1, ng2
        mVector(2) = iy - 1
        ! Wrap around if index is greater than 1/2 the table size
        IF (mVector(2) > ng2/2) mVector(2)=mVector(2)-ng2
     DO ix = 1, ng1
        mVector(1) = ix - 1
        !
        Ig=Ig+1
        ! Calculate the qPoint cartesian coordinates given by this mVector.
        !qPoint = MATMUL(mVector,recip_lat)
        qPoint = MATMUL(recip_lat,mVector)
        ! Assign this point in the qTable and qVectors.
        grid%gVec(1:3,Ig) = qPoint(:)
        grid%gVec(4,Ig) = Norm(qPoint)
        IF(((ix >= 2) &
           .OR. (ix == 1 .AND. iy >= 2 .AND. iy <= ng2/2+1) &
           .OR. (ix == 1 .AND. iy == 1 &
           .AND. iz >= 2 .AND. iz <= ng3/2+1))) THEN
           grid%gMask(Ig) = .TRUE.
        ELSE
           grid%gMask(Ig) = .FALSE.
        END IF
     ENDDO
     ENDDO
     ENDDO
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE FillQTable
  !--------------------symm-density-------------------
  SUBROUTINE Symm_density(rho)
     USE struct_module , ONLY : nsym !,Opsym,Otrans,lat_mat
     IMPLICIT NONE
     !IN/OUT
     REAL(DP),INTENT(INOUT) :: rho(nr1,nr2,nr3)
     !LOCAL
     INTEGER(I4B) :: Ix,Iy,Iz,Ip,Isy,Idir
     INTEGER(I4B) :: Ixyz(3),srxyz(3)  &
                  & , sizen(3) , indexs(3,nsym) ! , Nop
     REAL(DP) :: weight,rhot !weight
     LOGICAL  :: lflag(nr1,nr2,nr3),ltmp ! ,lsym(nsym)
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     IF(nsym<=1) RETURN
     !start symmetry the density
     !store
     sizen(1)=nr1
     sizen(2)=nr2
     sizen(3)=nr3
     !initialize
     lflag(:,:,:)=.FALSE.
     !set weight
     weight=1._DP/nsym
     !all real space points
     DO Iz=1,nr3
     DO Iy=1,nr2
     DO Ix=1,nr1
        IF(lflag(Ix,Iy,Iz)) CYCLE
        !set integer mesh coordinates and
        !move to center of cell in real coordinates
        Ixyz(1)=Ix
        Ixyz(2)=Iy
        Ixyz(3)=Iz
        !set rhot=0
        rhot=0._DP
        !cycle all symmetry operators
        DO Isy=1,nsym
           !operate on vector
           CALL Isymm_apply(Isy,sizen,Ixyz,srxyz,ltmp)
           !store index
           indexs(:,Isy)=srxyz(:)
           rhot=rhot+rho(srxyz(1),srxyz(2),srxyz(3))
           IF(.NOT.((rhot<1000.d0).AND.(rhot>0.d0)))THEN
               print*,rhot
               print*,srxyz
               STOP
           ENDIF
!print*,'xyz',Ixyz
!print*,'newxyz',srxyz
!pause
        ENDDO
        !apply weight
        rhot=weight*rhot
        DO Isy=1,nsym
           rho(indexs(1,Isy),indexs(2,Isy),indexs(3,Isy)) &
          &  = rhot
           lflag(indexs(1,Isy),indexs(2,Isy),indexs(3,Isy))=.TRUE.
        ENDDO
!        print*,'origin density',rhoin(Ix,Iy,Iz)
!        print*,'density',rho(Ix,Iy,Iz)
!pause
     ENDDO
     ENDDO
     ENDDO
     !test
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE Symm_density
  !---------------------------------------------------
  SUBROUTINE Isymm_apply(Isy,nsize,vin,vout,lsymm)
     USE math , ONLY : inv_33
     USE struct_module , ONLY : nsym,Opsym,Otrans,lat_mat
     IMPLICIT NONE
     !IN/OUT
     INTEGER(I4B),INTENT(IN)  :: Isy,nsize(3),vin(3)
     INTEGER(I4B),INTENT(OUT) :: vout(3) !out point
     LOGICAL,INTENT(OUT) :: lsymm
     !LOCAL
     REAL(DP) :: vec(3),svec(3) !temp
     REAL(DP) :: ntrans(3),Op(3,3),offset(3) !operator
     LOGICAL :: linlist !list of point
     REAL(DP),PARAMETER :: eps=1D-5 !
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     !set vec(in)
     offset(:)=1._DP
     vec(:)=(vin(:)-offset(:))/nsize(:)
     !point operater
     svec(:)=MATMUL(Opsym(:,:,Isy),vec(:))+Otrans(:,Isy)
     !translation
     !move back
     svec(:)=svec(:)*nsize(:)+offset(:)
     !Does it in list?
     linlist=ABS(svec(1)-NINT(svec(1)))<=eps .AND. &
           & ABS(svec(2)-NINT(svec(2)))<=eps .AND. &
           & ABS(svec(3)-NINT(svec(3)))<=eps
     IF(linlist)THEN
       !find the grid
       lsymm=.TRUE.
       !output new mesh
       vout(:) = MOD( NINT(svec(:)+10*nsize(:)) , nsize(:) )
       IF(vout(1)==0) vout(1)=nr1
       IF(vout(2)==0) vout(2)=nr2
       IF(vout(3)==0) vout(3)=nr3
       !test
       IF((vout(1)<0).OR.(vout(2)<0).OR.(vout(3)<0))THEN
          WRITE(*,*) 'Bug in symmetry density'
          STOP
       ENDIF
     ELSE
       !didn't find the mesh
       lsymm=.FALSE.
       !output origin mesh
       vout(:)=vin(:)
     ENDIF
!print*,'vout',vout
!pause
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE Isymm_apply
  !---------------------------------------------------
  SUBROUTINE sumrhoS(nps,rhoS,rho)
     USE parameters , ONLY : nspin
     IMPLICIT NONE
     !IN/OUT
     INTEGER(I4B)            :: nps
     REAL(DP),INTENT(IN)     :: rhoS(nps,nspin)
     REAL(DP),INTENT(OUT)    :: rho(nps)
     !LOCAL
     INTEGER(I4B) :: Is
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     IF(nspin==1)THEN
        rho(:)=rhoS(:,1)
     ELSEIF(nspin==2)THEN
        rho(:)=rhoS(:,1)+rhoS(:,2)
     ELSE
        STOP 'sumrhoS errors'
     ENDIF
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE sumrhoS
  !---------------------------------------------------
ENDMODULE grid_module
