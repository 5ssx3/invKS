!Find Cell symmetry
MODULE Cellsym_module
  USE constants , ONLY : DP, i4b
  type struct_spg
       integer(i4b)           :: nsg !!!! the Number of Space Group 
       integer(i4b)           :: ncs !!! the Number of Crystal System 
       integer(i4b)           :: setting ! the choice of coordinate system
       integer(i4b)           :: npg   !!! the Number of Points Group 
       integer(i4b)           :: np_op !!!! number of symmetry operation 
       character              :: C_center
       character(len=2)       :: bravis
       real(DP)               :: op_p(3,3,48) !!!! point symmetry operation
       real(DP)               :: op_t(3,48)   !!!! tanslation symmetry operation
       !!!!!!!!!!!!!!!!!the intrinsic symmetry of the structure!!!!!!!!!!!!!!
       integer(i4b)           :: s_op_p(3,3,48)
       real(DP)               :: s_op_t(3,48)
       real(DP)               :: T_nig2bla(3,3)
       character              :: s_center
  endtype
  type(struct_spg)               :: struct_sg
  !type struct_latpos_type
     !integer(i4b)               :: ntyp
     !integer(i4b),allocatable   :: typ(:)   
     !integer(i4b),allocatable   :: typ_unit(:)   
     !integer(i4b),allocatable   :: typ_prim(:)   
     !integer(i4b)               :: natom
     !real(DP)                   :: lat_old(3,3)
     !real(DP)                   :: lat_matrix(3,3)
     !real(DP)                   :: lat_unit(3,3)
     !real(DP)                   :: lat_prim(3,3)
     !real(DP),allocatable       :: pos_old(:,:)
     !real(DP),allocatable       :: pos(:,:)
     !real(DP),allocatable       :: pos_unit(:,:)
     !real(DP),allocatable       :: pos_prim(:,:)
  !end type 
  !type(struct_sp)               :: struct_sg
  !type(struct_latpos_type)      :: struct_latpos
CONTAINS
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !-------------------------PARTING LINE------------------------
   SUBROUTINE GetStructOp(ntyp, typ, natom, lat_matrix, pos, opoint_tol, trans_tol, &
         &  pgop_niggli,disp,num_t,C_i, num_ssym)
     USE constants, ONLY : dp, i4b
     IMPLICIT NONE
     !INPUT
     INTEGER(i4b),INTENT(IN) :: ntyp & !number of type
                           &,   natom  !total atoms
     INTEGER(i4b),INTENT(IN) :: typ(ntyp) 
     REAL(DP),INTENT(IN) :: lat_matrix(3,3) & !lattice matrix
                         &, pos(3,natom)    & !atomic position
                         &, trans_tol,opoint_tol !tolerance
     !OUT
     INTEGER(I4B),INTENT(OUT) :: num_ssym,num_t,C_i(8) !symmetry number
     REAL(DP),INTENT(OUT) :: pgop_niggli(3,3,48),disp(3,48) !symmetry Operator
     !WAIT
     INTEGER(i4b) :: num_psym,pgop(3,3,48)
     real(DP),dimension(3,3)    :: T_niggli
     REAL(DP)				   :: lat_niggli(3,3)
     real(DP)                   :: pos_niggli(3,natom)
     logical                   :: flag
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     call find_pointoperation(lat_matrix,num_psym,pgop,opoint_tol, flag)
     if (flag==.false.) goto 200
     call find_translation(ntyp, typ, natom, pos,lat_matrix,pgop,num_psym,num_t,C_i,trans_tol)
     pgop_niggli = 0.d0
     200 if (flag==.false.) then
     	struct_sg%np_op=1
     	struct_sg%op_p(:,:,1)=0.d0
     	struct_sg%op_p(1,1,1)=1.d0
     	struct_sg%op_p(2,2,1)=1.d0
     	struct_sg%op_p(3,3,1)=1.d0
     endif
     num_ssym = struct_sg%np_op
     pgop_niggli(:,:,1:num_ssym) = struct_sg%op_p(:,:,1:num_ssym)

     disp(:,1:num_ssym) = struct_sg%op_t(:,1:num_ssym)
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE GetStructOp
   !-------------------find_pointoperation---------------------
   SUBROUTINE find_pointoperation(lat_niggli,num_psym,pgop,opoint_tol,op_flag)
     !
     !USE constants
     implicit none
     ! in
     real(DP)                 :: lat_niggli(3,3)
     REAL(DP)				   :: opoint_tol
     ! out
     integer(i4b)             :: num_psym
     integer(i4b)             :: pgop(3,3,48)
     ! local
     integer(i4b)             :: i,j,k,ix1,iy1,iz1,ix2,iy2,iz2,ix3,iy3,iz3
     integer(i4b)             :: point_op(3,3)
     integer(i4b)             :: num_sym,det
     real(DP),dimension(3,3)  :: metric_orig,metric,lats
     logical                  :: flag,op_flag
     num_sym=0
     op_flag=.true.
     !call lat_mat(lat_niggli,metric_orig)
     metric_orig=matmul(lat_niggli,transpose(lat_niggli))
     do ix1=-1,1
        do ix2=-1,1
           do ix3=-1,1
              if (ix1==0.and.ix2==0.and.ix3==0) cycle
              do iy1=-1,1
                 do iy2=-1,1
                    do iy3=-1,1
                       if (iy1==0.and.iy2==0.and.iy3==0) cycle
                       do iz1=-1,1
                          if (ix1==0.and.iy1==0.and.iz1==0) cycle
                          do iz2=-1,1
                             if (ix2==0.and.iy2==0.and.iz2==0) cycle
                             do iz3=-1,1
                                if (iz1==0.and.iz2==0.and.iz3==0) cycle
                                if (ix3==0.and.iy3==0.and.iz3==0) cycle
                                point_op(1,1)=ix1
                                point_op(2,1)=ix2
                                point_op(3,1)=ix3
                                point_op(1,2)=iy1
                                point_op(2,2)=iy2
                                point_op(3,2)=iy3
                                point_op(1,3)=iz1
                                point_op(2,3)=iz2
                                point_op(3,3)=iz3
                                det=sym_det_i3(point_op)
                                if (det==1 .or. det==-1) then
                                   !call mul_di3(point_op,lat_niggli,lats)
                                   lats=matmul(point_op,lat_niggli)
                                   !call lat_mat(lats,metric)
                                   metric=matmul(lats,transpose(lats))
                                   call identity(metric,metric_orig,opoint_tol,flag)
                                   if (flag) then
                                      num_sym=num_sym+1
                                      if (num_sym > 48) then
                                        !write(*,*) "tolerance too large"
                                         op_flag=.false.
                                         return
                                      endif
                                      pgop(:,:,num_sym)=transpose(point_op(:,:))
   								   !do i=1,3
                                        !print *, pgop(i,:,num_sym)
   								   !enddo
                                   endif
                                endif
                             enddo
                          enddo
                       enddo
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
     num_psym=num_sym
     !print *, 'num_psym',num_psym
   ENDSUBROUTINE find_pointoperation
   !--------------------------identity---------------------------
   subroutine identity(matrix,matrix_orig,tol,flag)
      USE constants,            ONLY : DP,i4b
      implicit none
      real(DP)     :: a(3,3),b(3,3)
      real(DP)     :: angle_tolerance,tol
      real(DP)     :: cos1,cos2,sinx,sin_d,matrix(3,3),matrix_orig(3,3)
      real(DP)     :: length_av,length(3),length_orig(3)
      integer(i4b) :: i,j,l
      logical      :: flag
      flag=.true.
      angle_tolerance=tol
      do i=1,3
         length(i)=sqrt(matrix(i,i))
         length_orig(i)=sqrt(matrix_orig(i,i))
         !print *,"lent", abs(length(i)-length_orig(i))
         if (abs(length(i)-length_orig(i))> tol) then
            flag=.false.
            return
         endif
      enddo
      do i=1,3
         j=i+1
         if (j==4) j=1
           !print *, "ang12",get_angle(matrix,i,j),get_angle(matrix_orig,i,j) 
           !print *, "ang",abs(get_angle(matrix,i,j)-get_angle(matrix_orig,i,j)) 
           if (abs(get_angle(matrix,i,j)-get_angle(matrix_orig,i,j)) > angle_tolerance) then
              flag=.false.
              return
           endif
      enddo
   ENDSUBROUTINE identity
   !-------------------------get_angle---------------------------
   FUNCTION get_angle(matrix,i,j)
      !
      USE constants, ONLY : DP,i4b
      implicit none
      real(dp), parameter :: pi = 3.14159265358979323846_dp
      real(DP)     :: a(3,3),b(3,3)
      real(DP)     :: get_angle,li,lj,matrix(3,3)
      integer(i4b) :: i,j
      li=sqrt(matrix(i,i))
      lj=sqrt(matrix(j,j))
      get_angle=acos(matrix(i,j)/li/lj)/pi
   ENDFUNCTION get_angle
   !----------------------find_translation-----------------------
subroutine find_translation(ntyp, typ, natom, pos,lat_matrix, pgop,num_psym,num_t,C_i,trans_tol)
   !

   implicit none
   ! in
   integer(i4b)          :: num_psym
   INTEGER(I4B)						:: natom
   integer(i4b)          :: pgop(3,3,48)
   integer(i4b)               :: ntyp
   integer(i4b)               :: typ(ntyp)
   real(DP)                   :: lat_matrix(3,3)
   real(DP)                   :: pos(3,natom)
   REAL(DP)							:: trans_tol
   integer(i4b)          :: C_i(8)
   integer(i4b)          :: ele,num_ssym,num_t,id
   ! local
   integer(i4b)          :: i,j,i1,j1,i2,j2,k,m,num
   real(DP)              :: op(3,3),t_all(3,100),tot_dist(100)
   real(DP)              :: origin(3),pos_test(3),trans_test(3)
   real(DP)              :: dis(3),dis_lat(3),t_trans(3),tau(3)
   real(DP)              :: trans_al(3,8)
   real(DP)              :: dis2,dist_sum
   logical               :: flag,t_flag
   character             :: center
   INTEGER(I4B)						:: ia1,ia2
   INTEGER(I4B)						:: eleid(ntyp+1)
   t_flag=.true.
   eleid(1)=1
   do i =1, ntyp
      eleid(i+1)=eleid(i)+typ(i)
   enddo
   flag=.false.
   ele=minloc(typ,1)
   num_ssym=0
   do num=1,num_psym
      op(:,:)=real(pgop(:,:,num),DP)
      !print *,"op",op
      !call mul_31_id3(op,pos(:,1,ele),origin)
      origin=matmul(op,pos(:,eleid(ele)))
      do k=1,3
         origin(k)=origin(k)-floor(origin(k))
      enddo
      !print *,"origin",origin
      num_t=0
      loop1:do i=eleid(ele),eleid(ele+1)-1
         do j=1,3
            trans_test(j)=pos(j,i)-origin(j)
         enddo
         !print *,"pos",pos(:,i,ele)
         !print *,"i",i,ele
         !print *,"trans",trans_test
         dist_sum=0.d0
         do i1=1,ntyp
            do j1=eleid(i1),eleid(i1+1)-1
               !call mul_31_id3(op,pos(:,j1,i1),pos_test)
               pos_test=matmul(op,pos(:,j1))
               pos_test=pos_test+trans_test
               t_flag=.false.
               do j2=eleid(i1),eleid(i1+1)-1
                  dis=pos_test-pos(:,j2)
                  do k=1,3
                     dis(k)=dis(k)-nint(dis(k))
                  enddo
                  !call mul_vm(dis,lat_matrix,dis_lat)
                  !dis_lat=dis
                  dis_lat=matmul(dis,lat_matrix)
                  dis2=dis_lat(1)*dis_lat(1)+dis_lat(2)*dis_lat(2)+dis_lat(3)*dis_lat(3)
                  if (dis2 < trans_tol*trans_tol ) then
                     dist_sum=dist_sum+dis2
                     t_flag=.true.
                     exit
                  endif
               enddo
               if (t_flag==.false.) cycle loop1
            enddo
         enddo
         if (t_flag) then
            num_t=num_t+1
            !print *,"numttt",num_t
            t_all(:,num_t)=trans_test(:)
            tot_dist(num_t)=dist_sum
         endif
      enddo loop1
      !pause
      !if (t_flag) then
      if (num_t>0) then
         !print *,'num',num
         !print *,'nttt',num_t
         num_ssym=num_ssym+1
         id=minloc(tot_dist(1:num_t),1)
         t_trans=t_all(:,id)
         do k=1,3
            t_trans(k)=t_trans(k)-floor(t_trans(k))
            if (t_trans(k)>0.999d0) t_trans(k)=0.d0
         enddo
         struct_sg%s_op_t(:,num_ssym)=t_trans(:)
         struct_sg%s_op_p(:,:,num_ssym)=op(:,:)
         !do i=1,3
         !   print *, op(i,:)
         !enddo
      endif
   enddo 
   struct_sg%np_op=num_ssym
   struct_sg%op_p(:,:,1:num_ssym)=real(struct_sg%s_op_p(:,:,1:num_ssym),kind=DP)
   struct_sg%op_t(:,1:num_ssym)=struct_sg%s_op_t(:,1:num_ssym)
   trans_al=reshape((/&
      & 0.5d0, 0.5d0, 0.5d0, &
      & 0.0d0, 0.5d0, 0.5d0, &
      & 0.5d0, 0.0d0, 0.5d0, &
      & 0.5d0, 0.5d0, 0.0d0, &
      & 2.d0/3, 1.d0/3, 1.d0/3, &
      & 1.d0/3, 2.d0/3, 2.d0/3, &
      & 1.d0/3, 2.d0/3, 1.d0/3, &
      & 2.d0/3, 1.d0/3, 2.d0/3 /),(/3,8/))
   num_t=0
   lnum:do num=1,8
      t_trans=trans_al(:,num)
      do i=1,ntyp
         do j=eleid(i),eleid(i+1)-1
            pos_test=pos(:,j)+t_trans
            t_flag=.false.
            do j1=eleid(i),eleid(i+1)-1
               if (j1==j) cycle
               dis=pos_test-pos(:,j1)
               do k=1,3
                  dis(k)=dis(k)-nint(dis(k))
               enddo
               dis_lat=matmul(dis,lat_matrix)
               dis2=dis_lat(1)*dis_lat(1)+dis_lat(2)*dis_lat(2)+dis_lat(3)*dis_lat(3)
               dis2=sqrt(dis2)
               if (dis2 < trans_tol ) then
                  t_flag=.true.
                  exit
               endif
            enddo
            if (t_flag==.false.) then
               cycle lnum
            endif
         enddo
      enddo
      if (t_flag) then
         num_t=num_t+1
         C_i(num_t)=num
      endif
   enddo lnum
   !print *,"num_t",num_t
   !print *,'trans_tol',trans_tol
   if (num_t==0) then
      center='P'
   else if (num_t==1) then
      if (C_i(1)==1) then
         center='I'
      else if (C_i(1)==2) then
         center='A'
      else if (C_i(1)==3) then
         center='B'
      else if (C_i(1)==4) then
         center='C'
      endif
   else if (num_t==2) then
      center='R'
   else
      center='F'
   endif
   struct_sg%C_center=center
   !print *,"center",center
   !pause
   !!!!!!!!!!!!!!!!!!!
   !open(111,file="oldopandT")
   !do i=1,num_ssym
   !write(111,*) "num=",i
   !do j=1,3
   !write(111,*) struct_sg%s_op_p(j,:,i)
   !enddo
   !write(111,*) struct_sg%op_t(:,i)
   !enddo
   !close(111)
   !!print *,"op_t",struct_sg%op_t(:,1:num_ssym)
   !print *, 'num_ssym',num_ssym
end subroutine 
   !------------------------Math maybe need----------------------
function sym_det_i3(a)

  implicit none

  integer(i4b) :: a(3,3)
  integer(i4b) :: sym_det_I3,det_i3
   det_I3=0
   det_I3=det_I3+ a(1,1) * (a(2,2) * a(3,3) - a(2,3) * a(3,2)) &
                + a(1,2) * (a(2,3) * a(3,1) - a(2,1) * a(3,3)) &
                + a(1,3) * (a(2,1) * a(3,2) - a(2,2) * a(3,1))
   sym_det_I3=det_I3
end
   !-------------------------PARTING LINE------------------------
FUNCTION sym_det(matrix) 


  IMPLICIT NONE

  real(DP),  intent(in)  :: Matrix(3,3)
  real(DP) :: sym_det


  sym_det = Matrix(1,1)*(Matrix(2,2)*Matrix(3,3)-Matrix(2,3)*Matrix(3,2))&
        -Matrix(1,2)*(Matrix(2,1)*Matrix(3,3)-Matrix(2,3)*Matrix(3,1))&
        +Matrix(1,3)*(Matrix(2,1)*Matrix(3,2)-Matrix(3,1)*Matrix(2,2))

END FUNCTION
   !-------------------------PARTING LINE------------------------
subroutine sym_dis_mat(A,B,dist) 
    implicit none
    integer(i4b)    :: i,j,k
    real(DP)        :: A(3,3),B(3,3)
    real(DP)        :: dist
    dist=0.0
    do i=1,3
       do j=1,3
          dist=dist+abs(A(i,j)-B(i,j))
       enddo
    enddo
end subroutine
   !-------------------------PARTING LINE------------------------
   subroutine sym_inv_mat(A,invA)
      implicit none
      integer(i4b)    :: i,j,k,l,is(3),js(3)
      real(DP)        :: A(3,3),invA(3,3),A_old(3,3)
      real(DP)        :: t,d
      A_old=A
      l=1
      do k=1,3
         d=1D-3
         do i=k,3
            do j=k,3
               if (abs(a(i,j)).gt.d) then
                  d=abs(a(i,j))
                  is(k)=i
                  js(k)=j
               endif
            enddo
         enddo
         if ((d+1.0).eq.1.0) then
            l=0
         endif
         do j=1,3
            t=a(k,j)
            a(k,j)=a(is(k),j)
            a(is(k),j)=t
         enddo
         do i=1,3
            t=a(i,k)
            a(i,k)=a(i,js(k))
            a(i,js(k))=t
         enddo
         a(k,k)=1/a(k,k)
         do j=1,3
            if(j.ne.k) then
               a(k,j)=a(k,j)*a(k,k)
            endif
         enddo
         do i=1,3
            if (i.ne.k) then
               do j=1,3
                  if (j.ne.k) then
                     a(i,j)=a(i,j)-a(i,k)*a(k,j)
                  endif
               enddo
            endif
         enddo
         do i=1,3
            if(i.ne.k) then
               a(i,k)=-a(i,k)*a(k,k)
            endif
         enddo
      enddo
      do k=3,1,-1
         do j=1,3
            t=a(k,j)
            a(k,j)=a(js(k),j)
            a(js(k),j)=t
         enddo
         do i=1,3
            t=a(i,k)
            a(i,k)=a(i,is(k))
            a(i,is(k))=t
         enddo
      enddo
      invA=A
      A=A_old
   end subroutine
   !-------------------------PARTING LINE------------------------
subroutine sym_lat2matrix(lat2,matrix,flag)
  !
  !
  implicit none
  !
  integer(i4b)        :: i,j,k,n,flag
  real(kind=DP)   :: lat2(6),matrix(3,3),alfa 
  real(DP) :: ra,rb,rc,cosinea,cosineb,cosinec,anglea,angleb,anglec
  if (flag==1) then
     matrix=0.0
     matrix(1,1) = lat2(1)
     matrix(2,1) = lat2(2)*cos(lat2(6))
     matrix(2,2) = lat2(2)*sin(lat2(6))
     matrix(3,1) = lat2(3)*cos(lat2(5))
     matrix(3,2) = lat2(3)*cos(lat2(4))*sin(lat2(6))-((lat2(3)*cos(lat2(5))-lat2(3)*cos(lat2(4))*cos(lat2(6)))/tan(lat2(6)))
     matrix(3,3) = sqrt(lat2(3)**2 -matrix(3,1)**2 - matrix(3,2)**2)
  else
     lat2=0.0
     ra=sqrt(matrix(1,1)**2+matrix(1,2)**2+matrix(1,3)**2)
     rb=sqrt(matrix(2,1)**2+matrix(2,2)**2+matrix(2,3)**2)
     rc=sqrt(matrix(3,1)**2+matrix(3,2)**2+matrix(3,3)**2)
     cosinea=(matrix(2,1)*matrix(3,1)+matrix(2,2)*matrix(3,2)+matrix(2,3)*matrix(3,3))/rb/rc
     cosineb=(matrix(1,1)*matrix(3,1)+matrix(1,2)*matrix(3,2)+matrix(1,3)*matrix(3,3))/ra/rc
     cosinec=(matrix(1,1)*matrix(2,1)+matrix(1,2)*matrix(2,2)+matrix(1,3)*matrix(2,3))/ra/rb
     anglea=acos(cosinea)
     angleb=acos(cosineb)
     anglec=acos(cosinec)
     lat2(1)=ra
     lat2(2)=rb
     lat2(3)=rc
     lat2(4)=anglea
     lat2(5)=angleb
     lat2(6)=anglec
  endif
end subroutine 
   !-------------------------PARTING LINE------------------------
subroutine det_pg(pg_flag)
   !
   implicit none
   ! out
   logical               :: pg_flag
   ! local
   integer(i4b)          :: ncs ! the number of crystal system
   !    !!!!!!!!!!!!!!!!!!!!!!!!crystal system!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !    ! 1: Triclinic   2: Monolinic    3: Orthorhombic  4: Tetragonal !
   !    ! 5: Trigonal    6: Hexagonal    7: Cubic                       !
   !    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer(i4b)          :: npg  ! the number of point group
   integer(i4b)          :: axis ! the rotation type 
   integer(i4b)          :: np_op! the number of operations
   integer(i4b)          :: op_p(3,3,48)
   integer(i4b)          :: k,tr_w,det_w,nm
   integer(i4b)          :: rot4,rot_4,rot3,rot_3,rot2,rot_2,rot6,rot_6
   logical               :: invE,rot_flag
   npg=1
   ncs=1
   np_op=struct_sg%np_op
   op_p(:,:,1:np_op)=struct_sg%s_op_p(:,:,1:np_op)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!! determine crystal system  !!!!!!!!!!!!
   !!!!     Acta Cryst.(1999) A55, 383-395     !!!!!!
   !!!!                Table IV                !!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   rot6=0;rot_6=0;rot4=0;rot_4=0
   rot3=0;rot_3=0;rot2=0;rot_2=0
   InvE=.false.
   pg_flag=.true.
   ncs=0;npg=0;
   do k=1,np_op
      det_w=sym_det_i3(op_p(:,:,k)) ! the determinant of operation matrix
      tr_w=op_p(1,1,k)+op_p(2,2,k)+op_p(3,3,k) ! trace of operation matrix
      if (tr_w==-3) invE=.true.
      call rot_typ(tr_w,det_w,axis,rot_flag)! determine the rotation part
      if(rot_flag==.false.) then
         pg_flag=.false.
         return
      endif
      if (axis==6) then
         rot6=rot6+1
      else if (axis==-6) then
         rot_6=rot_6+1
      else if(axis==3) then
         rot3=rot3+1
      else if(axis==-3) then
         rot_3=rot_3+1
      else if (axis==4) then
         rot4=rot4+1
      else if (axis==-4) then
         rot_4=rot_4+1
      else if(axis==2) then
         rot2=rot2+1
      else if(axis==-2) then
         rot_2=rot_2+1
      endif
   enddo
   !print *,"rot2",rot2
   !print *,"rot3",rot3
   !print *,"rot4",rot4
   !print *,"rot6",rot6
   !print *,"rot_2",rot_2
   !print *,"rot_3",rot_3
   !print *,"rot_4",rot_4
   !print *,"rot_6",rot_6
   if (rot3>=4) then
      ncs=7
   else if (rot6>=1 .or. rot_6>=1) then
      ncs=6
   else if (rot3>=1 .or. rot_3>=1) then
      ncs=5
   else if (rot4>=1 .or. rot_4>=1) then
      ncs=4
   else if (rot2>=2 .or. rot_2>=2) then
      ncs=3
   else if (rot2>=1 .or. rot_2>=1) then
      ncs=2
   else
      ncs=1
   endif
   if (ncs==0) then
      pg_flag=.false.
      return
   endif
   !!!!!!!!!!!Ensure Point Group!!!!!!!!
   nm=np_op
   select case(ncs)
   case(1)
      if (invE) then
         npg=2
      else
         npg=1
      endif
   case(2)
      if (invE) then
         npg=5
      else if(rot_2==1) then
         npg=4
      else if(rot2==1) then
         npg=3
      endif
   case(3)
      if (invE) then
         npg=8
      else if(rot_2==2) then
         npg=7
      else if(rot2==3) then
         npg=6
      endif
   case(4)
      if (invE) then
         nm=np_op/2
         if (nm==8) then
            npg=15
         else if(nm==4) then
            npg=11
         endif
      else
         if (nm==8) then
            if (rot_4==2) then
               npg=14
            else if(rot_2==4) then
               npg=13
            else if(rot4==2 .and. rot2==5) then
               npg=12
            endif
         else if (nm==4) then
            if (rot_4==2) then
               npg=10
            else if(rot4==2) then
               npg=9
            endif
         endif
      endif
   case(5)
      if (invE) then
         nm=np_op/2
         if (nm==6) then
            npg=20
         else if(nm==3) then
            npg=17
         endif
      else 
         if (nm==6) then
            if (rot_2==3) then
               npg=19
            else if (rot2==3) then
               npg=18
            endif
         else if (nm==3) then
            npg=16
         endif
      endif
   case(6)
      if (invE) then
         nm=np_op/2
         if (nm==12) then
            npg=27
         else if(nm==6) then
            npg=23
         endif
      else
         if (nm==12) then
            if (rot_6==2) then
               npg=26
            else if (rot_2==6) then
               npg=25
            else if (rot6==2 .and. rot2==7) then
               npg=24
            endif
         else if (nm==6) then
            if (rot_6==2) then
               npg=22
            else if (rot6==2) then
               npg=21
            endif
         endif
      endif
   case(7) 
      if (invE) then
         nm=np_op/2
         if (nm==24) then
            npg=32
         else if (nm==12) then
            npg=29
         endif
      else 
         if (nm==24) then
            if (rot_4==6) then
               npg=31
            else if (rot4==6) then
               npg=30
            endif
         else if (nm==12) then
            npg=28
         endif
      endif
   end select
   if (npg==0) then
      pg_flag=.false.
      return
   endif
   struct_sg%npg=npg
   struct_sg%ncs= ncs
   !print *,"npg",npg
   !print *,"ncs",ncs
end subroutine
   !-------------------------PARTING LINE------------------------
subroutine rot_typ(tr_w,det_w,rot_c,rot_flag)
   implicit none
   !in
   integer(i4b)               :: tr_w,det_w
   !out
   logical                    :: rot_flag
   !local
   integer(i4b)               :: rot_c
   rot_flag=.true.
   if (det_w==1) then
      select case(tr_w)
      case(-1)
         rot_c=2
      case(0)
         rot_c=3
      case(1)
         rot_c=4
      case(2)
         rot_c=6
      case(3)
         rot_c=1
      case default
         write(*,*)"ERROR!!!The trace of op is wrong"
         rot_flag=.false.
         return
      end select
   elseif (det_w==-1) then
      select case(tr_w)
      case(-3)
         rot_c=-1
      case(-2)
         rot_c=-6
      case(-1)
         rot_c=-4
      case(0)
         rot_c=-3
      case(1)
         rot_c=-2
      case default
         write(*,*)"ERROR!!!The trace of op is wrong"
         rot_flag=.false.
         return
      end select
   else
      write(*,*)"ERROR!!!The det of op is wrong"
      rot_flag=.false.
      return
   endif
end subroutine
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDMODULE Cellsym_module
