MODULE smpi_math_module

   USE CONSTANTS
   USE mpi
   IMPLICIT NONE
      INTEGER(I4B) :: mykgrid_range(3)
      INTEGER(I4B), ALLOCATABLE :: recvcounts(:)
      INTEGER(I4B), ALLOCATABLE :: displs(:)
      INTEGER(I4B), ALLOCATABLE :: gridrange_sum(:,:)
   !INCLUDE 'mpif.h'
   type parallel_type  !{{{
      INTEGER(I4B)              :: comm
      INTEGER(I4B)              :: myid
      INTEGER(I4B)              :: numprocs
      INTEGER(I4B)              :: rootid
      LOGICAL                   :: isroot
      !>> number of states needed by per process
      INTEGER(I4B)              :: nstate_proc
      !>> array_map('id in alignment of n_states','process id')
      INTEGER(I4B),ALLOCATABLE  :: sub2sum(:,:)
      INTEGER(I4B)              :: mygrid_range(3)
      INTEGER(I4B),ALLOCATABLE  :: recvcounts(:)
      INTEGER(I4B),ALLOCATABLE  :: displs(:)
      INTEGER(I4B),ALLOCATABLE  :: global_gridrange(:,:)
      INTEGER(I4B),ALLOCATABLE  :: comm2d_rank2sum(:,:)
      INTEGER(I4B)              :: comm2d,commx,commy,rankx,ranky,  &
           & periods(2),reorder,ndims,dims(2), &
           & commx_numprocs, commy_numprocs, commx_myid, commy_myid
      LOGICAL :: remainX(2), remainY(2)
   end type parallel_type

   type (parallel_type)         :: parallel
   !-----------------------------------------------------------------------
   TYPE smpi_root_type
      INTEGER(I4B), allocatable :: natom_group(:, :)
   end type smpi_root_type
   type(smpi_root_type)         :: smpi_root
   !-----------------------------------------------------------------------
   TYPE smpi_comm_type
      INTEGER(I4B), allocatable :: atoms(:)
      INTEGER(I4B), allocatable :: displs(:)
   end type smpi_comm_type
   type( smpi_comm_type )       :: smpi_comm
   !-----------------------------------------------------------------------
   type time_type
      character(len=100)       :: label
      real(dp)                  :: tic
      real(dp)                  :: toc
      real(dp)                  :: total
      real(dp)                  :: sum_total
      INTEGER(I4B)              :: num
   end type time_type
   !-----------------------------------------------------------------------
   type mem_type
      character(len=100)       :: label
      real(dp)                  :: memic
      real(dp)                  :: total
      INTEGER(I4B)              :: num
   end type mem_type
   !-----------------------------------------------------------------------
   type(time_type),ALLOCATABLE  :: timedat(:)
   !type(time_type)  :: timedat(100)
   type(mem_type)               :: memdat(100)
   INTEGER(I4B),private,save    :: ntime=0
   REAL(DP)                     :: rtic,rtoc
   !-----------------------------------------------------------------------
   TYPE grid_diff_map_type
      INTEGER(I4B),allocatable :: nz_map(:,:) !> the up id and down id for per nz
      INTEGER(I4B) :: mycomm_cores(2)         !> number of cores for communcation
      INTEGER(I4B),allocatable :: mycomm_size(:,:)
      INTEGER(I4B),allocatable :: mysend_size(:,:)
      INTEGER(I4B),allocatable :: local_map(:,:)
      INTEGER(I4B)             :: boundary(2,3)
   ENDTYPE grid_diff_map_type
   !-----------------------------------------------------------------------
   TYPE sphere_map
      INTEGER(I4B)             :: Length !> sphere size
      INTEGER(I4B),ALLOCATABLE :: map3d(:,:) !> (4,length)
   ENDTYPE sphere_map
   !-----------------------------------------------------------------------
   type(sphere_map) :: sphere
   type(grid_diff_map_type) :: diff_map
   !-----------------------------------------------------------------------
   INTEGER(I4B)                 :: mpinfo
   INTEGER(I4B)                 :: smpi_status(MPI_STATUS_SIZE)  !}}}
   !------------------------------- DIVIDER LINE --------------------------------
INTERFACE SompSum  !{{{
   MODULE PROCEDURE sum_real_1d
   MODULE PROCEDURE sum_real_2d
   MODULE PROCEDURE sum_real_3d
   MODULE PROCEDURE sum_cplx_1d
   MODULE PROCEDURE sum_cplx_2d
   MODULE PROCEDURE sum_cplx_3d
END INTERFACE  

INTERFACE SmpiSum  
   MODULE PROCEDURE smpi_sum_int_1s
   MODULE PROCEDURE smpi_sum_cplx_1s
   MODULE PROCEDURE smpi_sum_real_1s
   MODULE PROCEDURE smpi_sum_real_1d
   MODULE PROCEDURE smpi_sum_real_2d
   MODULE PROCEDURE smpi_sum_real_3d
END INTERFACE  

! INTERFACE SmpiSumSub 
!   MODULE PROCEDURE smpi_sum_int_1s_sub
!   MODULE PROCEDURE smpi_sum_cplx_1s_sub
   ! MODULE PROCEDURE smpi_sum_real_1s_sub 
!    MODULE PROCEDURE smpi_sum_real_1d_sub
!    MODULE PROCEDURE smpi_sum_real_2d_sub
!    MODULE PROCEDURE smpi_sum_real_3d_sub
! END INTERFACE  

! INTERFACE SmpiSumSlab
!    !MODULE PROCEDURE smpi_sum_int_1s
!    !MODULE PROCEDURE smpi_sum_cplx_1s
!    MODULE PROCEDURE smpi_sum_real_1s_slab
!    MODULE PROCEDURE smpi_sum_real_1d
!    !MODULE PROCEDURE smpi_sum_real_2d
!    !MODULE PROCEDURE smpi_sum_real_3d
! END INTERFACE  

INTERFACE SmpiSumMem
    MODULE PROCEDURE smpi_sum_mem_1d
    MODULE PROCEDURE smpi_sum_mem_2d
    MODULE PROCEDURE smpi_sum_mem_3d
END INTERFACE
       

INTERFACE SmpiReduceSum 
   MODULE PROCEDURE smpi_reduce_sum_real_1d
   MODULE PROCEDURE smpi_reduce_sum_int_1d
   MODULE PROCEDURE smpi_reduce_sum_cplx_1d
   MODULE PROCEDURE smpi_reduce_sum_real_2d
END INTERFACE  

! INTERFACE SmpiReduceSumSub
!    MODULE PROCEDURE smpi_reduce_sum_real_1d_sub
!    MODULE PROCEDURE smpi_reduce_sum_cplx_1d_sub
!    MODULE PROCEDURE smpi_reduce_sum_real_2d_sub
! !   MODULE PROCEDURE smpi_reduce_sum_cplx_1d
! !   MODULE PROCEDURE smpi_reduce_sum_real_2d
! END INTERFACE  

! INTERFACE SmpiReduceSumSlab
!    MODULE PROCEDURE smpi_reduce_sum_real_1d_slab
!    MODULE PROCEDURE smpi_reduce_sum_cplx_1d_slab
! !   MODULE PROCEDURE smpi_reduce_sum_cplx_1d
! !   MODULE PROCEDURE smpi_reduce_sum_real_2d
! END INTERFACE  

INTERFACE SompSumPow
   module procedure sum_pow_int
   module procedure sum_pow_real
end interface

INTERFACE SmpiSumPow
   module procedure smpi_sum_pow_int
   module procedure smpi_sum_pow_real
end interface!}}}
!---------------------------- DIVIDER LINE -----------------------------
CONTAINS
  SUBROUTINE smpi_init()
    IMPLICIT NONE

    CALL MPI_INIT(mpinfo)

    CALL MPI_COMM_RANK(MPI_COMM_WORLD, parallel%myid, mpinfo)
    CALL MPI_Comm_SIZE(MPI_COMM_WORLD, parallel%numprocs, mpinfo)

    parallel%comm  = MPI_COMM_WORLD

    parallel%rootid = 0
    if(parallel%rootid == parallel%myid)then
       parallel%isroot = .true.
    else
       parallel%isroot = .false.
    endif
  END SUBROUTINE smpi_init
  !-----------------------divided line---------------------------
   SUBROUTINE smpi_init_2D()
      USE grid_module , ONLY : nk
      IMPLICIT NONE
      INTEGER(I4B) :: coords(2), ierr
      INTEGER(I4B) :: i,rank_local_x,rank_local_y,rank_global
      !CALL smpi_init()
      parallel%ndims = 2
      parallel%dims = [0, 0]
      call MPI_Dims_create(parallel%numprocs, parallel%ndims, parallel%dims, mpinfo)
      IF (mpinfo /= MPI_SUCCESS) THEN
         IF (parallel%isroot) WRITE(6, '(A)') "[ERROR] MPI_Dims_create failed"
         CALL MPI_Abort(parallel%comm, 103, ierr)
      ENDIF
      parallel%periods = [.false., .false.]
      parallel%reorder = .false.
      CALL MPI_Cart_create(parallel%comm, parallel%ndims, parallel%dims, parallel%periods, parallel%reorder, parallel%comm2d, mpinfo)
      IF (mpinfo /= MPI_SUCCESS) THEN
         IF (parallel%isroot) WRITE(6, '(A)') "[ERROR] MPI_Cart_create failed"
         CALL MPI_Abort(parallel%comm, 103, ierr)
      ENDIF
      CALL MPI_Cart_coords(parallel%comm2d, parallel%myid, parallel%ndims, coords, mpinfo)
      IF (mpinfo /= MPI_SUCCESS) THEN
         IF (parallel%isroot) WRITE(6, '(A)') "[ERROR] MPI_Cart_coords failed"
         CALL MPI_Abort(parallel%comm, 104, ierr)
      ENDIF
      parallel%rankx = coords(1) + 1
      parallel%ranky = coords(2) + 1
      !
      parallel%remainX = [.TRUE., .FALSE.]
      parallel%remainY = [.FALSE., .TRUE.]
      CALL MPI_Cart_sub(parallel%comm2d, parallel%remainX, parallel%commx, mpinfo)
      IF (mpinfo /= MPI_SUCCESS) THEN
         IF (parallel%isroot) WRITE(6, '(A)') "[ERROR] MPI_Cart_sub failed for commx"
         CALL MPI_Abort(parallel%comm, 105, ierr)
      ENDIF
      CALL MPI_Cart_sub(parallel%comm2d, parallel%remainY, parallel%commy, mpinfo)
      IF (mpinfo /= MPI_SUCCESS) THEN
         IF (parallel%isroot) WRITE(6, '(A)') "[ERROR] MPI_Cart_sub failed for commy"
         CALL MPI_Abort(parallel%comm, 106, ierr)
      ENDIF
      !
      CALL MPI_Comm_size(parallel%commx, parallel%commx_numprocs, mpinfo)
      CALL MPI_Comm_size(parallel%commy, parallel%commy_numprocs, mpinfo)
      CALL MPI_Comm_rank(parallel%commx, parallel%commx_myid, mpinfo)
      CALL MPI_Comm_rank(parallel%commy, parallel%commy_myid, mpinfo)
      !
      IF (ALLOCATED(parallel%comm2d_rank2sum)) DEALLOCATE(parallel%comm2d_rank2sum)
      ALLOCATE(parallel%comm2d_rank2sum(parallel%commx_numprocs,parallel%commy_numprocs))
      parallel%comm2d_rank2sum = -2
      rank_global = 1
      DO rank_local_x = 1,parallel%commx_numprocs
         DO rank_local_y = 1,parallel%commy_numprocs
             parallel%comm2d_rank2sum(rank_local_x,rank_local_y)=rank_global
             rank_global = rank_global + 1
         ENDDO
      ENDDO
      !PRINT*,'parallel%comm2d_rank2sum',parallel%comm2d_rank2sum(parallel%rankx,parallel%ranky),'=','process',parallel%myid+1
      IF (parallel%isroot) PRINT*, "2D Topology created: dims=", parallel%dims
   END SUBROUTINE smpi_init_2D
  !-----------------------divided line---------------------------
  !-----------------------divided line---------------------------
  !-----------------------divided line---------------------------

   SUBROUTINE smpi_exit()!{{{
      IMPLICIT NONE
      call MPI_finalize( mpinfo )
      ! STOP
      RETURN
   END SUBROUTINE!}}}

   SUBROUTINE smpi_stop(message)!{{{
      IMPLICIT NONE

      CHARACTER (LEN=*) message

      WRITE (*,*) message

      call MPI_abort(parallel%comm, 1, mpinfo )
      STOP

      RETURN
   END SUBROUTINE!}}}

   SUBROUTINE smpi_stop_info(message)!{{{
      IMPLICIT NONE

      CHARACTER (LEN=*) message

      WRITE (*,*) message, mpinfo

      call MPI_abort(MPI_comm_world , 1, mpinfo )
      STOP

      RETURN
   END SUBROUTINE!}}}
   !---------------------------- DIVIDER LINE -----------------------------
 !  Subroutine  nstates_split(m,np)  !{{{
 !  implicit none 
 !  !integer(i4b) , intent(inout) :: m
 !  integer(i4b) :: np
 !  !-------------------------------------
 !  ! ms = mod(m ,np)
 !  ! if (ms == 0 ) then 
 !  !   m = m /np
 !  ! else 
 !  !   !if ( ms >= (parallel%myid + 1) ) then 
 !  !   !  m = m/parallel%numprocs+1
 !  !   !else 
 !  !   !  m = m / parallel%numprocs
 !  !   !end if 
 !  !   !if (parallel%myid == 0 ) then 
 !  !   if (parallel%coords(1) == 0 ) then 
 !  !     m = m/np + ms 
 !  !   else 
 !  !     m = m / np
 !  !   end if 
 !  ! end if 
 !  
 !  End Subroutine nstates_split  !}}}
   !---------------------------- DIVIDER LINE -----------------------------
 !  Subroutine  nstates_split_2(m,np)  !{{{
 !  implicit none 
 !  !integer(i4b) , intent(inout) :: m
 !  integer(i4b) :: np
 !  !-------------------------------------
 !  ! ms = mod(m ,np)
 !  ! if (ms == 0 ) then 
 !  !   m = m /np
 !  ! else 
 !  !   !if ( ms >= (parallel%myid + 1) ) then 
 !  !   !  m = m/parallel%numprocs+1
 !  !   !else 
 !  !   !  m = m / parallel%numprocs
 !  !   !end if 
 !  !   !if (parallel%myid == 0 ) then 
 !  !   if (parallel%coords(2) == 0 ) then 
 !  !     m = m/np + ms 
 !  !   else 
 !  !     m = m / np
 !  !   end if 
 !  ! end if 
 !  
 !  End Subroutine nstates_split_2  !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   FUNCTION sum_real_1d(amat) result(totals)!{{{
      REAL(DP)                      :: totals
      REAL(DP)                      :: amat(:,:,:)
      INTEGER(I4B)                  :: n1,n2,n3
      INTEGER(I4B)                  :: i,j,k
      n1 = size(amat,1)
      n2 = size(amat,2)
      n3 = size(amat,3)
      totals = 0.d0
      !$OMP PARALLEL DO REDUCTION(+:totals) PRIVATE(k,j,i)
      DO k = 1, n3
         DO j = 1, n2
            DO i = 1, n1
               totals = totals + amat(i,j,k)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END PARALLEL DO
   END FUNCTION !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   FUNCTION sum_cplx_1d(amat) result(totals)!{{{
      complex(DP)                      :: totals
      complex(DP)                      :: amat(:,:,:)
      INTEGER(I4B)                  :: n1,n2,n3
      INTEGER(I4B)                  :: i,j,k
      n1 = size(amat,1)
      n2 = size(amat,2)
      n3 = size(amat,3)
      totals = 0.d0
      !$OMP PARALLEL DO REDUCTION(+:totals) PRIVATE(k,j,i)
      DO k = 1, n3
         DO j = 1, n2
            DO i = 1, n1
               totals = totals + amat(i,j,k)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END PARALLEL DO
   END FUNCTION !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   FUNCTION sum_real_2d(amat, bmat) result(totals)!{{{
      REAL(DP)                      :: totals
      REAL(DP)                      :: amat(:,:,:)
      REAL(DP)                      :: bmat(:,:,:)
      INTEGER(I4B)                  :: n1,n2,n3
      INTEGER(I4B)                  :: i,j,k
      n1 = size(amat,1)
      n2 = size(amat,2)
      n3 = size(amat,3)
      totals = 0.d0
      !$OMP PARALLEL DO REDUCTION(+:totals) PRIVATE(k,j,i)
      DO k = 1, n3
         DO j = 1, n2
            DO i = 1, n1
               totals = totals + amat(i,j,k)*bmat(i,j,k)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END PARALLEL DO
   END FUNCTION !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   FUNCTION sum_cplx_2d(amat, bmat) result(totals)!{{{
      complex(DP)                      :: totals
      complex(DP)                      :: amat(:,:,:)
      complex(DP)                      :: bmat(:,:,:)
      INTEGER(I4B)                  :: n1,n2,n3
      INTEGER(I4B)                  :: i,j,k
      n1 = size(amat,1)
      n2 = size(amat,2)
      n3 = size(amat,3)
      totals = 0.d0
      !$OMP PARALLEL DO REDUCTION(+:totals) PRIVATE(k,j,i)
      DO k = 1, n3
         DO j = 1, n2
            DO i = 1, n1
               totals = totals + amat(i,j,k)*bmat(i,j,k)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END PARALLEL DO
   END FUNCTION !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   FUNCTION sum_real_3d( amat, bmat, cmat) result(totals)!{{{
      REAL(DP)                      :: totals
      REAL(DP)                      :: amat(:,:,:)
      REAL(DP)                      :: bmat(:,:,:),cmat(:,:,:)
      INTEGER(I4B)                  :: n1,n2,n3
      INTEGER(I4B)                  :: i,j,k
      n1 = size(amat,1)
      n2 = size(amat,2)
      n3 = size(amat,3)
      totals = 0.d0
      !$OMP PARALLEL DO REDUCTION(+:totals) PRIVATE(k,j,i)
      DO k = 1, n3
         DO j = 1, n2
            DO i = 1, n1
               totals = totals + amat(i,j,k)*bmat(i,j,k)*cmat(i,j,k)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END PARALLEL DO
   END FUNCTION !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   FUNCTION sum_cplx_3d( amat, bmat, cmat) result(totals)!{{{
      complex(DP)                      :: totals
      complex(DP)                      :: amat(:,:,:)
      complex(DP)                      :: bmat(:,:,:),cmat(:,:,:)
      INTEGER(I4B)                  :: n1,n2,n3
      INTEGER(I4B)                  :: i,j,k
      n1 = size(amat,1)
      n2 = size(amat,2)
      n3 = size(amat,3)
      totals = 0.d0
      !$OMP PARALLEL DO REDUCTION(+:totals) PRIVATE(k,j,i)
      DO k = 1, n3
         DO j = 1, n2
            DO i = 1, n1
               totals = totals + amat(i,j,k)*bmat(i,j,k)*cmat(i,j,k)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END PARALLEL DO
   END FUNCTION !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   FUNCTION sum_pow_real(amat, pow) result(totals)!{{{
      REAL(DP)                      :: totals
      REAL(DP)                      :: amat(:,:,:)
      REAL(DP)                      :: pow
      INTEGER(I4B)                  :: n1,n2,n3
      INTEGER(I4B)                  :: i,j,k
      n1 = size(amat,1)
      n2 = size(amat,2)
      n3 = size(amat,3)
      totals = 0.d0
      !$OMP PARALLEL DO REDUCTION(+:totals) PRIVATE(k,j,i)
      DO k = 1, n3
         DO j = 1, n2
            DO i = 1, n1
               totals = totals + amat(i,j,k)** pow
            ENDDO
         ENDDO
      ENDDO
      !$OMP END PARALLEL DO
   END FUNCTION !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   FUNCTION sum_pow_int(amat, pow)  result(totals)!{{{
      REAL(DP)                      :: totals
      REAL(DP)                      :: amat(:,:,:)
      INTEGER(I4B)                  :: pow
      INTEGER(I4B)                  :: n1,n2,n3
      INTEGER(I4B)                  :: i,j,k
      n1 = size(amat,1)
      n2 = size(amat,2)
      n3 = size(amat,3)
      totals = 0.d0
      !$OMP PARALLEL DO REDUCTION(+:totals) PRIVATE(k,j,i)
      DO k = 1, n3
         DO j = 1, n2
            DO i = 1, n1
               totals = totals + amat(i,j,k)** pow
            ENDDO
         ENDDO
      ENDDO
      !$OMP END PARALLEL DO
   END FUNCTION !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   function smpi_sum_int_1s(x) result(sumx)!{{{

      IMPLICIT NONE
      INTEGER(I4B)  :: x,sumx

      CALL MPI_ALLREDUCE(x, sumx, 1, MPI_INTEGER4, MPI_SUM, parallel%comm, mpinfo)

   end function smpi_sum_int_1s!}}}
   !---------------------------- DIVIDER LINE -----------------------------
   function smpi_sum_cplx_1s(x) result(sumx)!{{{

      IMPLICIT NONE
      complex(DP)                     :: x,sumx

      CALL MPI_ALLREDUCE(x, sumx, 1, MPI_REAL8, MPI_SUM, parallel%comm, mpinfo)

   end function smpi_sum_cplx_1s!}}}
   !---------------------------- DIVIDER LINE -----------------------------
   function smpi_sum_real_1s(x) result(sumx)!{{{

      IMPLICIT NONE
      REAL(DP)                     :: x,sumx

      CALL MPI_ALLREDUCE(x, sumx, 1, MPI_REAL8, MPI_SUM, parallel%comm, mpinfo)

   end function smpi_sum_real_1s!}}}
   !---------------------------- DIVIDER LINE -----------------------------
   ! function smpi_sum_real_1s_sub(x) result(sumx)!{{{

   !    IMPLICIT NONE
   !    REAL(DP)                     :: x,sumx

   !    CALL MPI_ALLREDUCE(x, sumx, 1, MPI_REAL8, MPI_SUM, parallel%subcomm, mpinfo)

   ! end function smpi_sum_real_1s_sub!}}}
   !---------------------------- DIVIDER LINE -----------------------------
   ! function smpi_sum_real_1s_slab(x) result(sumx)!{{{

   !    IMPLICIT NONE
   !    REAL(DP)                     :: x,sumx

   !    CALL MPI_ALLREDUCE(x, sumx, 1, MPI_REAL8, MPI_SUM, parallel%slabcomm, mpinfo)

   ! end function smpi_sum_real_1s_slab!}}}
   !---------------------------- DIVIDER LINE -----------------------------
   FUNCTION smpi_sum_real_1d(amat) result(suma)!{{{
      REAL(DP)                      :: suma, totals
      REAL(DP)                      :: amat(:,:,:)
      totals = sum_real_1d(amat)
      call mpi_allreduce(totals, suma, 1, mpi_real8, mpi_sum, parallel%comm, mpinfo)
   END FUNCTION !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   ! FUNCTION smpi_sum_real_1d_sub(amat) result(suma)!{{{
   !    REAL(DP)                      :: suma, totals
   !    REAL(DP)                      :: amat(:,:,:)
   !    totals = sum_real_1d(amat)
   !    call mpi_allreduce(totals, suma, 1, mpi_real8, mpi_sum, parallel%subcomm, mpinfo)
   ! END FUNCTION !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   FUNCTION smpi_sum_real_2d(amat, bmat) result(suma)!{{{
      REAL(DP)                      :: suma, totals
      REAL(DP)                      :: amat(:,:,:)
      REAL(DP)                      :: bmat(:,:,:)
      totals = sum_real_2d(amat,bmat)
      call mpi_allreduce(totals, suma, 1, mpi_real8, mpi_sum, parallel%comm, mpinfo)
   END FUNCTION !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   ! FUNCTION smpi_sum_real_2d_sub(amat, bmat) result(suma)!{{{
   !    REAL(DP)                      :: suma, totals
   !    REAL(DP)                      :: amat(:,:,:)
   !    REAL(DP)                      :: bmat(:,:,:)
   !    totals = sum_real_2d(amat,bmat)
   !    call mpi_allreduce(totals, suma, 1, mpi_real8, mpi_sum, parallel%subcomm, mpinfo)
   ! END FUNCTION !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   FUNCTION smpi_sum_real_3d( amat, bmat, cmat) result(suma)!{{{
      REAL(DP)                      :: suma, totals
      REAL(DP)                      :: amat(:,:,:)
      REAL(DP)                      :: bmat(:,:,:),cmat(:,:,:)
      totals = sum_real_2d(amat,bmat)
      call mpi_allreduce(totals, suma, 1, mpi_real8, mpi_sum, parallel%comm, mpinfo)
   END FUNCTION !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   ! FUNCTION smpi_sum_real_3d_sub( amat, bmat, cmat) result(suma)!{{{
   !    REAL(DP)                      :: suma, totals
   !    REAL(DP)                      :: amat(:,:,:)
   !    REAL(DP)                      :: bmat(:,:,:),cmat(:,:,:)
   !    totals = sum_real_2d(amat,bmat)
   !    call mpi_allreduce(totals, suma, 1, mpi_real8, mpi_sum, parallel%subcomm, mpinfo)
   ! END FUNCTION !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   FUNCTION smpi_sum_pow_real(amat, pow) result(suma)!{{{
      REAL(DP)                      :: suma, totals
      REAL(DP)                      :: amat(:,:,:)
      REAL(DP)                      :: pow
      totals = sum_pow_real(amat,pow)
      call mpi_allreduce(totals, suma, 1, mpi_real8, mpi_sum, parallel%comm, mpinfo)
   END FUNCTION !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   FUNCTION smpi_sum_pow_int(amat, pow)  result(suma)!{{{
      REAL(DP)                      :: suma, totals
      REAL(DP)                      :: amat(:,:,:)
      INTEGER(I4B)                  :: pow
      totals = sum_pow_int(amat,pow)
      call mpi_allreduce(totals, suma, 1, mpi_real8, mpi_sum, parallel%comm, mpinfo)
   END FUNCTION !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   SUBROUTINE smpi_reduce_sum_int_1d(amat,ramat) !{{{
      INTEGER(I4B)                      :: amat(:)
      INTEGER(I4B),OPTIONAL             :: ramat(:)
      INTEGER(I4B), dimension(size(amat)) :: asuma
      INTEGER(I4B)                 :: na
      na = size(amat)
      call mpi_allreduce(amat, asuma, na, mpi_integer4, mpi_sum, parallel%comm, mpinfo)
      if (present(ramat)) then
         ramat=asuma
      else
         amat=asuma
      endif
   END SUBROUTINE !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   SUBROUTINE smpi_reduce_sum_real_1d(amat,ramat) !{{{
      REAL(DP)                      :: amat(:)
      REAL(DP),OPTIONAL            :: ramat(:)
      REAL(DP), dimension(size(amat)) :: asuma
      INTEGER(I4B)                 :: na
      na = size(amat)
      call mpi_allreduce(amat, asuma, na, mpi_real8, mpi_sum, parallel%comm, mpinfo)
      if (present(ramat)) then
         ramat=asuma
      else
         amat=asuma
      endif
   END SUBROUTINE !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   SUBROUTINE smpi_reduce_sum_real(amat,na,ramat) !{{{
      REAL(DP)                      :: amat(na)
      REAL(DP),OPTIONAL            :: ramat(na)
      REAL(DP), dimension(size(amat)) :: asuma
      INTEGER(I4B)                 :: na
      call mpi_allreduce(amat, asuma, na, mpi_real8, mpi_sum, parallel%comm, mpinfo)
      if (present(ramat)) then
         ramat=asuma
      else
         amat=asuma
      endif
   END SUBROUTINE !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   SUBROUTINE smpi_reduce_sum_cplx_1d(amat,ramat) !{{{
      complex(DP)                      :: amat(:)
      complex(DP),OPTIONAL            :: ramat(:)
      complex(DP), dimension(size(amat)) :: asuma
      INTEGER(I4B)                 :: na
      na = size(amat)
      call mpi_allreduce(amat, asuma, na, MPI_COMPLEX16, mpi_sum, parallel%comm, mpinfo)
      if (present(ramat)) then
         ramat=asuma
      else
         amat=asuma
      endif
   END SUBROUTINE !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   SUBROUTINE smpi_reduce_sum_real_2d(amat,ramat) !{{{
      REAL(DP)                      :: amat(:,:)
      REAL(DP),OPTIONAL            :: ramat(:,:)
      REAL(DP), dimension(size(amat,1),size(amat,2)) :: asuma
      INTEGER(I4B)                 :: nxy
      nxy = size(amat,1) * size(amat,2)
      call mpi_allreduce(amat, asuma, nxy, mpi_real8, mpi_sum, parallel%comm, mpinfo)
      if (present(ramat)) then
         ramat=asuma
      else
         amat=asuma
      endif
   END SUBROUTINE !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   subroutine enlarge_timedat()
     implicit none
     INTEGER(I4B) :: size_timedat
     type(time_type),allocatable   :: temp_timedat(:)

     !> enlarge timedat
     size_timedat=size(timedat)
     if(ntime.ge.size_timedat)then
        allocate(temp_timedat(size_timedat+100))
        temp_timedat(1:size_timedat)=timedat
        deallocate(timedat)
        allocate(timedat(size_timedat+100))
        timedat(1:size_timedat)=temp_timedat(1:size_timedat)
     endif
   endsubroutine enlarge_timedat
   !---------------------------- DIVIDER LINE -----------------------------
   subroutine start_time(inlabel,flag,tic)!{{{
      implicit none
      character(len=*) :: inlabel
      character(len=100) :: label
      LOGICAL  :: flag
      REAL(DP),optional :: tic
      INTEGER(I4B)                 :: i
   

      if( (.not.allocated(timedat)) ) allocate(timedat(100))
      if( (.not.flag) ) return

      label = trim(inlabel)
      do i = 1, ntime
         if (label.EQ.timedat(i)%label) then
            timedat(i)%tic= mpi_wtime()
            if (present(tic)) tic = timedat(i)%tic
            timedat(i)%num=timedat(i)%num+1
            return
         endif
      enddo
      if(ntime.ge.size(timedat)) call enlarge_timedat()
      ntime=ntime+1
      timedat(ntime)%label = label
      timedat(ntime)%tic= mpi_wtime()
      timedat(ntime)%sum_total= 0.d0
      if (present(tic)) tic = timedat(ntime)%tic
   end subroutine !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   subroutine end_time(inlabel,flag,toc)!{{{
      implicit none
      character(len=*) :: inlabel
      character(len=100) :: label
      LOGICAL  :: flag
      REAL(DP),optional :: toc
      INTEGER(I4B)                 :: i

      if(.not.flag)return
      label = trim(inlabel)
      do i = 1, ntime
         if (label==timedat(i)%label) then
            timedat(i)%toc= mpi_wtime()
            if (present(toc)) toc = timedat(i)%toc
            timedat(i)%total = timedat(i)%toc - timedat(i)%tic
            timedat(i)%sum_total = timedat(i)%sum_total + timedat(i)%total
            return
         endif
      enddo
      if (parallel%isroot) then
         WRITE(6,*) "ERROR: Haven't found ", inlabel
      endif
   end subroutine !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   subroutine write_time(inlabel,flag)!{{{
      character(len=*) :: inlabel
      character(len=100) :: label
      LOGICAL            :: flag
      INTEGER(I4B)                 :: i

      if(.not.flag)return
      ! if (.not. parallel%isroot) return
      label = trim(inlabel)
      do i = 1, ntime
         if (label==timedat(i)%label) then
            WRITE(6,'(A30,2X,ES12.4,1X,1A,2X,I10,I3)') inlabel, timedat(i)%total,'s' &
                & , timedat(i)%num,parallel%myid
            return
         endif
      enddo
   end subroutine write_time!}}}
   !---------------------------- DIVIDER LINE -----------------------------
   subroutine write_sum_time(inlabel,flag)!{{{
      character(len=*) :: inlabel
      character(len=100) :: label
      LOGICAL            :: flag
      INTEGER(I4B)                 :: i

      if(.not.flag)return
      ! if (.not. parallel%isroot) return
      label = trim(inlabel)
      do i = 1, ntime
         if (label==timedat(i)%label) then
            WRITE(6,*) inlabel," >> ",timedat(i)%sum_total
            return
         endif
      enddo
    end subroutine write_sum_time!}}}
   !---------------------------- DIVIDER LINE -----------------------------
   subroutine print_time(inlabel,t)!{{{
      character(len=*) :: inlabel
      character(len=100) :: label
      real(dp) :: t
      INTEGER(I4B)                 :: i
      if (.not. parallel%isroot) return
      label = trim(inlabel)
      do i = 1, ntime
         if (label==timedat(i)%label) then

            t = timedat(i)%toc-timedat(i)%tic
            return
         endif
      enddo
   end subroutine print_time!}}}
   !---------------------------- DIVIDER LINE -----------------------------
   FUNCTION smpi_sum_mem_1d(munit,amat) result(summem)!{{{
      REAL(DP)                      :: summem, mem1
      REAL(DP)                      :: amat(:)
      character(8) :: munit
      if (munit == 'G') then
        mem1 = size(amat,1)*8/1024/1024/1024
      else if (munit == 'M') then
        mem1 = size(amat,1)*8/1024/1024
      end if 

      call mpi_allreduce(mem1, summem, 1, mpi_real8, mpi_sum, parallel%comm, mpinfo)

   END FUNCTION !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   FUNCTION smpi_sum_mem_2d(munit,amat) result(summem)!{{{
      REAL(DP)                      :: summem, mem1
      REAL(DP)                      :: amat(:,:)
      character(8) :: munit
      if (munit == 'G' .or. munit == 'g' ) then
        mem1 = size(amat,1)*size(amat,2)*8/1024/1024/1024
      else if (munit == 'M' .or. munit == 'm' ) then
        mem1 = size(amat,1)*size(amat,2)*8/1024/1024
      end if 

      call mpi_allreduce(mem1, summem, 1, mpi_real8, mpi_sum, parallel%comm, mpinfo)

   END FUNCTION !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   FUNCTION smpi_sum_mem_3d(munit,amat) result(summem)!{{{
      REAL(DP)                      :: summem, mem1
      REAL(DP)                      :: amat(:,:,:)
      character(8) :: munit
      if (munit == 'G' .or. munit == 'g' ) then
        mem1 = size(amat,1)*size(amat,2)*size(amat,3)*8/1024/1024/1024
      else if (munit == 'M' .or. munit == 'm' ) then
        mem1 = size(amat,1)*size(amat,2)*size(amat,3)*8/1024/1024
      end if 

      call mpi_allreduce(mem1, summem, 1, mpi_real8, mpi_sum, parallel%comm, mpinfo)

   END FUNCTION smpi_sum_mem_3d !}}}
   !-------------------------------------------------------------------------
   ! Subroutine states_split(nev)
   !   USE parameters , ONLY: BLOCK_MBNB
   !   IMPLICIT NONE
   !   !>> total states to calculate
   !   !INTEGER(I4B)             :: ncore
   !   INTEGER(I4B),INTENT(INOUT)  ::  nev
   !   !>
   !   INTEGER(I4B)             :: average,redund
   !   !>=================================
   !   !>> spili eigen solver
   !   average=nev/(parallel%dims(1)*BLOCK_MBNB)
   !   IF(average==0)THEN
   !      IF(parallel%isroot)WRITE(6,*)"[ ERROR ]: states not enough, more states or less cores"
   !      stop
   !   ENDIF
   !   redund=mod(nev,parallel%dims(1)*BLOCK_MBNB)
   !   IF(parallel%myid*BLOCK_MBNB.lt.redund)THEN
   !      nev=average*BLOCK_MBNB+min(BLOCK_MBNB,redund-parallel%myid*BLOCK_MBNB)
   !   ELSE
   !      nev=average*BLOCK_MBNB
   !   ENDIF
   ! END Subroutine states_split
   !---------------------------------------------------------------------
   ! Subroutine array_split(nev)
   !   USE parameters , ONLY: BLOCK_MBNB
   !   IMPLICIT NONE
   !   !>> total states to calculate
   !   INTEGER(I4B),INTENT(IN)  ::  nev
   !   INTEGER(I4B)             :: average,redund
   !   ! INTEGER(I4B)             :: i,j,counter,n_temp
   !   INTEGER(I4B) :: i,j,n_temp,seq_local,seq_global
   !   !>=================================
   !   !>> spili eigen solver
   !   average=nev/(parallel%dims(1)*BLOCK_MBNB)
   !   IF(average==0)THEN
   !      IF(parallel%isroot)WRITE(6,*)"[ ERROR ]: smaller block or parallel cores"
   !      stop
   !   ENDIF
   !   redund=mod(nev,parallel%dims(1)*BLOCK_MBNB)
   !   IF(parallel%ranky*BLOCK_MBNB.lt.redund)THEN
   !      parallel%nstate_proc=average*BLOCK_MBNB+min(BLOCK_MBNB,redund-parallel%ranky*BLOCK_MBNB)
   !   ELSE
   !      parallel%nstate_proc=average*BLOCK_MBNB
   !   ENDIF 
   !   !>> set array map associated processes
   !   ALLOCATE(parallel%sub2sum((average+1)*BLOCK_MBNB,parallel%dims(1)))
   !   !>> set initial value
   !   parallel%sub2sum=-2
   !   seq_local=0
   !   seq_global=0
   !   setmap:DO j=1,nev,1
   !      seq_local = seq_local + 1
   !      DO i=1,parallel%dims(1),1
   !         seq_global = seq_global + 1
   !         parallel%sub2sum(seq_local,i) = seq_global
   !         if(seq_global == nev)exit setmap
   !      ENDDO
   !   ENDDO setmap

   !   ! counter=0
   !   ! !>> assignment
   !   ! DO j=1,parallel%numprocs,1
   !   !    IF(parallel%myid.eq.j-1)n_temp=parallel%nstate_proc
   !   !    CALL MPI_BCAST(n_temp,1,MPI_INTEGER4,j-1,parallel%comm,mpinfo)
   !   ! DO i=1,(average+1)*BLOCK_MBNB,1
   !   !    IF(i.le.n_temp)THEN
   !   !       counter=counter+1
   !   !       parallel%sub2sum(i,j)=counter
   !   !    ENDIF
   !   ! ENDDO
   !   ! ENDDO
   !   ! print *,"---->",parallel%nstate_proc,parallel%myid
   !   ! IF(parallel%isroot)print *,"sub2sum",parallel%sub2sum
   ! ENDSUBROUTINE array_split
   !---------------------------------------------------------------------
   Subroutine array_split(nev, ncore, comm, id, nstate_proc, sub2sum)
   IMPLICIT NONE
   !in/out
   INTEGER(I4B),INTENT(IN)  ::  nev, ncore, comm, id
   INTEGER(I4B),INTENT(OUT) :: nstate_proc
   INTEGER(I4B), INTENT(OUT), ALLOCATABLE :: sub2sum(:,:)
   !local
   INTEGER(I4B)             :: average,redund
   INTEGER(I4B) :: i, j, n_temp, seq_local, seq_global, max_nstate
   INTEGER(I4B) :: ierr
   !
   IF (nev < 1) THEN
      IF (parallel%isroot) WRITE(6, '(A,I0)') "[ERROR] array_split: nev must be ≥1 (got ", nev, ")"
      CALL MPI_ABORT(parallel%comm, 20, ierr)
   ENDIF
   IF (ncore < 1) THEN
      IF (parallel%isroot) WRITE(6, '(A,I0)') "[ERROR] array_split: ncore must be ≥1 (got ", ncore, ")"
      CALL MPI_ABORT(parallel%comm, 20, ierr)
   ENDIF
   !
   average=nev/ncore
   redund=mod(nev,ncore)
   IF(average==0)THEN
     IF(parallel%isroot)WRITE(6,*)"[ ERROR ]: smaller block or parallel cores"
     CALL MPI_ABORT(parallel%comm, 20, ierr)
   ENDIF
   !
   IF(id.lt.redund)THEN
      nstate_proc=average+1
   ELSE
      nstate_proc=average
   ENDIF 
   max_nstate = MERGE(average+1, average, redund > 0)
   IF (ALLOCATED(sub2sum)) THEN
      DEALLOCATE(sub2sum, STAT=ierr)
      IF (ierr /= 0 .AND. parallel%isroot) THEN
         WRITE(6, '(A,I0)') "[ERROR] array_split: Deallocate sub2sum failed (ierr=", ierr, ")"
         CALL MPI_ABORT(parallel%comm, 21, ierr)
      ENDIF
   ENDIF
   ALLOCATE(sub2sum(max_nstate, ncore), STAT=ierr)
   IF (ierr /= 0) THEN
      IF (parallel%isroot) WRITE(6, '(A,I0,A,I0)') &
         "[ERROR] array_split: Allocate sub2sum failed (size=", max_nstate, "×", ncore, "), ierr=", ierr
      CALL MPI_ABORT(parallel%comm, 22, ierr)
   ENDIF
   sub2sum=-2
   seq_local=0
   seq_global=0
   setmap:DO j=1,ncore,1
            IF (j-1 .lt. redund) THEN
               n_temp=average+1
            ELSE
               n_temp=average
            ENDIF
            DO i= 1,n_temp,1
               seq_local = seq_local + 1
               seq_global = seq_global+1
               sub2sum(seq_local,j) = seq_global
               if(seq_global == nev) exit setmap
            ENDDO
            seq_local=0
      ENDDO setmap

   IF (seq_global /= nev .AND. parallel%isroot) THEN
      WRITE(6, '(A,I0,A,I0)') "[WARNING] array_split: Mapped states (", seq_global, ") ≠ total nev (", nev, ")"
   ENDIF
  ENDSUBROUTINE array_split
   !---------------------------------------------------------------------
   SUBROUTINE grid_split(ngrid,ncore,comm,id,grid_range,recvcounts,displs,gridrange_sum)
     IMPLICIT NONE
     INTEGER(I4B),INTENT(IN)  :: ngrid,ncore,comm,id
     INTEGER(I4B),INTENT(OUT) :: grid_range(3)
     INTEGER(I4B),INTENT(OUT) :: recvcounts(ncore)
     INTEGER(I4B),INTENT(OUT) :: displs(ncore)
     INTEGER(I4B),INTENT(OUT),optional :: gridrange_sum(3,ncore)
     !> local
     INTEGER(I4B) :: i,j
     INTEGER(I4B) :: average,redund
     INTEGER(I4B) :: recvcounts_temp(ncore)
     INTEGER(I4B) :: displs_temp(ncore)
     INTEGER(I4B) :: mpinfo
     !> split the grid into "ncore" cores
     average=ngrid/ncore
     redund=mod(ngrid,ncore)
     if(id.lt.redund)then
        !> id*(l+1) + 1, ..., (id*(l+1)+1) + (l+1) - 1
        grid_range(1)=id*(average+1)+1
        grid_range(2)=(id+1)*(average+1)
        grid_range(3)=average+1
     else
        !> redund*(l+1)+(id-redund)*l+1, ..., (redund*(l+1)+(id-redund)*l+1) + l - 1
        grid_range(1)=redund+id*average+1
        grid_range(2)=redund+(id+1)*average
        grid_range(3)=average
     endif

     !> set the recvcounts and displs
     do i=0,ncore-1,1
        j=i+1
        if(i.lt.redund)then
           !> id*(l+1) + 1, ..., (id*(l+1)+1) + (l+1) - 1
           recvcounts(j)=average+1
           displs(j)=i*(average+1)+1-1
        else
           !> redund*(l+1)+(id-redund)*l+1, ..., (redund*(l+1)+(id-redund)*l+1) + l - 1
           recvcounts(j)=average
           displs(j)=redund+i*average+1-1
        endif
     enddo

     !> set the global_gridrange
     if(present(gridrange_sum))then
        recvcounts_temp=3
        displs_temp=(/(i*3,i=0,ncore-1,1)/)
        CALL MPI_ALLGATHERV(grid_range,3,MPI_INTEGER4,gridrange_sum,recvcounts_temp,&
             displs_temp,MPI_INTEGER4,comm,mpinfo)
     endif

     !> small test
     !print *,"myid",id,"start",grid_range(1),"end",grid_range(2),"size",grid_range(3)
     ! if(parallel%isroot)print*,"displs",displs
     ! if(parallel%isroot)print*,"recvcounts",recvcounts
   ENDSUBROUTINE grid_split
   !-------------------------------------------------------------------
   Subroutine atom_split(natom,mysize,atom_index)
     !USE struct_module , ONLY : natom
     IMPLICIT NONE
     !>> the states to calculate in cores
     INTEGER(I4B)             :: natom
     INTEGER(I4B)             :: atom_index(:)
     INTEGER(I4B)             :: average,redund
     INTEGER(I4B)             :: mysize
     INTEGER(I4B)             :: i,j,counter,n_temp
     !>=================================
     n_temp=0
     !>> splicit atoms
     average=natom/(parallel%numprocs)
     IF(average==0)THEN
        IF(parallel%isroot)WRITE(6,*)"[ ERROR ]: exist zero in the number of atoms of cores"
        ! stop
     ENDIF
     redund=mod(natom,parallel%numprocs)
     IF(parallel%myid.lt.redund)THEN
        mysize=average+1
     ELSE
        mysize=average
     ENDIF
     !>> set array map associated processes
     !ALLOCATE(atom_index(average+1))
     !>> set initial value
     atom_index=-2
     counter=0
     !>> assignment
     DO j=1,parallel%numprocs,1
        IF(parallel%myid.eq.j-1)n_temp=mysize
        CALL MPI_BCAST(n_temp,1,MPI_INTEGER4,j-1,parallel%comm,mpinfo)
        DO i=1,n_temp,1
           counter=counter+1
           IF(parallel%myid.eq.j-1)atom_index(i)=counter
        ENDDO
     ENDDO
     ! print *,"---->",parallel%nstate_proc,parallel%myid
     ! IF(parallel%isroot)print *,"sub2sum",parallel%sub2sum
   ENDSUBROUTINE atom_split

   ! SUBROUTINE grid_sphere_init(n1,n2,n3,norder,grid_pos,orign,radius,sphere,diff_map)

   SUBROUTINE set_wrap_grid_old(myrho,wrap_box)
     IMPLICIT NONE
     REAL(DP),INTENT(In) :: myrho(:)
     REAL(DP),INTENT(OUT) :: wrap_box(diff_map%boundary(1,1):diff_map%boundary(2,1),&
           & diff_map%boundary(1,2):diff_map%boundary(2,2),&
           & diff_map%boundary(1,3):diff_map%boundary(2,3))
     !> local
     INTEGER(I4B) :: i,j,len0
     REAL(DP),allocatable :: rho_send(:)
     REAL(DP),allocatable :: rho_recv(:)
     integer(I4B),save                 :: tag1=100
     integer(I4B)                 :: tag2
     integer(I4B),allocatable     :: rq(:),req(:)
     integer(I4B),allocatable     :: recflag(:)
     INTEGER(I4B) :: status1(MPI_STATUS_SIZE,diff_map%mycomm_cores(2):parallel%myid+1-1)
     INTEGER(I4B) :: status2(MPI_STATUS_SIZE,parallel%myid+2:diff_map%mycomm_cores(1))

     !> data communcation
     tag1=tag1+1
     tag2=200
     if(tag1.gt.tag2) tag1=100
     allocate(rq(diff_map%mycomm_cores(2):diff_map%mycomm_cores(1)))
     do i=diff_map%mycomm_cores(2),diff_map%mycomm_cores(1),1
        j=i-1
        if( parallel%rankx .EQ. j ) cycle
        len0=diff_map%mysend_size(2,i)-diff_map%mysend_size(1,i)+1
        ! allocate(rho_send(len))
        ! rho_send=myrho(diff_map%mysend_size(1,i):diff_map%mysend_size(2,i))
        call MPI_ISEND(myrho(diff_map%mysend_size(1,i)),len0,MPI_REAL8,j,&
             & tag1,parallel%commx,rq(i),mpinfo)
        ! deallocate(rho_send)
     enddo

!     call mpi_barrier(parallel%comm,mpinfo)

     allocate(rho_recv(diff_map%mycomm_size(1,diff_map%mycomm_cores(2)):diff_map%mycomm_size(2,diff_map%mycomm_cores(1))))
     allocate(req(diff_map%mycomm_cores(2):diff_map%mycomm_cores(1)))
     allocate(recflag(diff_map%mycomm_cores(2):diff_map%mycomm_cores(1)))

     do i=diff_map%mycomm_cores(1),diff_map%mycomm_cores(2),-1
        j=i-1
        if( parallel%rankx .EQ. j )cycle
        len0=diff_map%mycomm_size(2,i)-diff_map%mycomm_size(1,i)+1
        call MPI_IRECV(rho_recv(diff_map%mycomm_size(1,i)),len0,MPI_REAL8,j,&
             & tag1,parallel%commx,req(i),mpinfo)
     enddo

!     call mpi_barrier(parallel%comm,mpinfo)
     !> deal with the local part
     wrap_box=0
     rho_recv(parallel%mygrid_range(1):parallel%mygrid_range(2))=myrho
     do i=parallel%mygrid_range(1),parallel%mygrid_range(2),1
        wrap_box(diff_map%local_map(1,i),diff_map%local_map(2,i),diff_map%local_map(3,i))=rho_recv(i)
     enddo

!!     call mpi_barrier(parallel%comm,mpinfo)

     if( (parallel%rankx+1-diff_map%mycomm_cores(2)) .NE. 0)then
        call MPI_WAITALL(parallel%rankx+1-diff_map%mycomm_cores(2),&
             & rq(diff_map%mycomm_cores(2):parallel%rankx+1-1),status1,mpinfo)
        call MPI_WAITALL(parallel%rankx+1-diff_map%mycomm_cores(2),&
             & req(diff_map%mycomm_cores(2):parallel%rankx+1-1),status1,mpinfo)
     endif
!print*,'warp4'
!
!!     call mpi_barrier(parallel%comm,mpinfo)

     if((diff_map%mycomm_cores(1)-parallel%rankx-1) .NE.0)then
        call MPI_WAITALL(diff_map%mycomm_cores(1)-parallel%rankx-1,&
             & rq(parallel%rankx+2:diff_map%mycomm_cores(1)),status2,mpinfo)
        call MPI_WAITALL(diff_map%mycomm_cores(1)-parallel%rankx-1,&
             & req(parallel%rankx+2:diff_map%mycomm_cores(1)),status2,mpinfo)
     endif
!print*,'warp5'

     !call mpi_barrier(parallel%comm,mpinfo)
     ! !> deal with the local part
     ! wrap_box=0
     ! rho_recv(parallel%mygrid_range(1):parallel%mygrid_range(2))=myrho
     ! do i=parallel%mygrid_range(1),parallel%mygrid_range(2),1
     !    wrap_box(diff_map%local_map(1,i),diff_map%local_map(2,i),diff_map%local_map(3,i))=rho_recv(i)
     ! enddo

     !> sphere segment change to wrap box
     do i=diff_map%mycomm_size(1,diff_map%mycomm_cores(2)),parallel%mygrid_range(1)-1,1
        wrap_box(diff_map%local_map(1,i),diff_map%local_map(2,i),diff_map%local_map(3,i))=rho_recv(i)
     enddo
!print*,'warp6'

     do i=parallel%mygrid_range(2)+1,diff_map%mycomm_size(2,diff_map%mycomm_cores(1)),1
        wrap_box(diff_map%local_map(1,i),diff_map%local_map(2,i),diff_map%local_map(3,i))=rho_recv(i)
     enddo
!print*,'warp7'

     call mpi_barrier(parallel%commx,mpinfo)
     !DEALLOCATE
     IF(ALLOCATED(rho_send)) DEALLOCATE(rho_send)
     IF(ALLOCATED(rho_recv)) DEALLOCATE(rho_recv)
     IF(ALLOCATED(rq))  DEALLOCATE(rq)
     IF(ALLOCATED(req)) DEALLOCATE(req)
     IF(ALLOCATED(recflag)) DEALLOCATE(recflag)
   ENDSUBROUTINE set_wrap_grid_old

   SUBROUTINE set_wrap_grid(myrho,wrap_box)
     IMPLICIT NONE
     REAL(DP),INTENT(IN) :: myrho(:)
     REAL(DP),INTENT(OUT) :: wrap_box(diff_map%boundary(1,1):diff_map%boundary(2,1),&
           & diff_map%boundary(1,2):diff_map%boundary(2,2),&
           & diff_map%boundary(1,3):diff_map%boundary(2,3))
     !> local
     INTEGER(I4B) :: i,j !,len0
     !REAL(DP),allocatable :: rho_send(:)
     REAL(DP),allocatable :: rho_recv(:)
     !integer(I4B),save                 :: tag1=100
     !integer(I4B)                 :: tag2
     !integer(I4B),allocatable     :: rq(:),req(:)
     !integer(I4B),allocatable     :: recflag(:)
     !INTEGER(I4B) :: status1(MPI_STATUS_SIZE,diff_map%mycomm_cores(2):parallel%myid+1-1)
     !INTEGER(I4B) :: status2(MPI_STATUS_SIZE,parallel%myid+2:diff_map%mycomm_cores(1))
     !> alltoallv setting
     INTEGER(I4B),dimension(parallel%numprocs) :: scount,sdispls&
          &, rcount,rdispls

     !> data communcation
     scount=0
     sdispls=0
     rcount=0
     rdispls=0
     do i=diff_map%mycomm_cores(2),diff_map%mycomm_cores(1),1
        j=i-1
        if( parallel%rankx .EQ. j )cycle
        scount(i)=diff_map%mysend_size(2,i)-diff_map%mysend_size(1,i)+1
        sdispls(i)=diff_map%mysend_size(1,i)-1
     enddo

     allocate(rho_recv(diff_map%mycomm_size(1,diff_map%mycomm_cores(2)) &
          &   :diff_map%mycomm_size(2,diff_map%mycomm_cores(1))))
!print*,'warp0'
     do i=diff_map%mycomm_cores(1),diff_map%mycomm_cores(2),-1
        j=i-1
        if( parallel%rankx .EQ. j )cycle
        rcount(i)=diff_map%mycomm_size(2,i)-diff_map%mycomm_size(1,i)+1
        rdispls(i)=diff_map%mycomm_size(1,i)-diff_map%mycomm_size(1,diff_map%mycomm_cores(2))
     enddo
!print*,'warp1'

     CALL MPI_ALLTOALLV(myrho,scount&
          & ,sdispls,MPI_REAL8&
          & ,rho_recv,rcount&
          & ,rdispls,MPI_REAL8&
          & ,parallel%commx,mpinfo)
!print*,'warp2'

     !> deal with the local part
     wrap_box=0
     rho_recv(parallel%mygrid_range(1):parallel%mygrid_range(2))=myrho
     do i=parallel%mygrid_range(1),parallel%mygrid_range(2),1
        wrap_box(diff_map%local_map(1,i),diff_map%local_map(2,i),diff_map%local_map(3,i))=rho_recv(i)
     enddo
!print*,'warp3'

     !> sphere segment change to wrap box
     do i=diff_map%mycomm_size(1,diff_map%mycomm_cores(2)),parallel%mygrid_range(1)-1,1
        wrap_box(diff_map%local_map(1,i),diff_map%local_map(2,i),diff_map%local_map(3,i))=rho_recv(i)
     enddo
     
     do i=parallel%mygrid_range(2)+1,diff_map%mycomm_size(2,diff_map%mycomm_cores(1)),1
        wrap_box(diff_map%local_map(1,i),diff_map%local_map(2,i),diff_map%local_map(3,i))=rho_recv(i)
     enddo

     call mpi_barrier(parallel%commx,mpinfo)
     deallocate(rho_recv)
!print*,'warp4'

   ENDSUBROUTINE set_wrap_grid
   !---------------------------- DIVIDER LINE -----------------------------
END MODULE smpi_math_module