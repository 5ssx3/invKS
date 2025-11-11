Module Math
   !##########################################################!{{{
   !* CREATED_TIME  : 2013-03-12
   !* AUTHOR        : Yanchao Wang
   !* CHANGE        : Xuecheng Shao
   !* ADD           : Qiang Xu 
   !* DESCRIPTION   :              
   !     ------
   !* REFERENCES    :              
   !     ------
   !* LOG           :              
   !     2017-07-24 :              
   !* LAST MODIFIED : 2017-07-24 11:24:35 AM
   !##########################################################!}}}
   use constants
   implicit none
   !INTERFACE gasdev!{{{
   !   MODULE PROCEDURE gasdev_s_sp, gasdev_s_dp, &
   !         gasdev_v_sp, gasdev_v_dp
   !END INTERFACE!}}}
   !---------------------------- DIVIDER LINE -----------------------------
   interface norm!{{{
      module procedure norm_real
      module procedure norm_complex
   end interface!}}}
   !---------------------------- DIVIDER LINE -----------------------------
   interface cross!{{{
      module procedure cross_real
      module procedure cross_complex
   end interface!}}}
   !---------------------------- DIVIDER LINE -----------------------------
   interface sort_id!{{{
      module procedure integer_index
      module procedure real_index
   end interface!}}}
   !---------------------------- DIVIDER LINE -----------------------------
CONTAINS
   SUBROUTINE change_case(instr,str,fun)!{{{
      !##########################################################
      !* CREATED_TIME  : 2015-05-12
      !* AUTHOR        : Xuecheng Shao
      !* CHANGE        : Xuecheng Shao
      !* ADD           : Xuecheng Shao
      !* DESCRIPTION   :              
      !     ------
      !* REFERENCES    :              
      !     ------
      !* LOG           :              
      !     2015-05-08 :              
      !* LAST MODIFIED : 2015-05-12 04:50:12 PM
      !##########################################################
      implicit none
      INTEGER(I4B)              :: fun
      character(len=*)     :: instr
      character(len=*)     :: str
      !
      INTEGER(I4B)              :: ia ! 32
      INTEGER(I4B)              :: i,l_str
      !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>
      l_str=len_trim(instr)
      str = instr
      ia = ichar('a') - ichar('A') ! 32
      DO i=1,l_str
         IF (fun == 1 ) then ! lower to upper
            !if(str(i:i) >= 'a' .and. str(i:i) <= 'z') then
            if( lge(instr(i:i),'a') .and. lle(instr(i:i),'z') ) then
               str(i:i)=char( ichar(instr(i:i)) - ia )
            else
               str(i:i)=instr(i:i)
            ENDIF
         ELSE ! upper to lowe
            !if(str(i:i) >= 'A' .and. str(i:i) <= 'Z') then
            if( lge(instr(i:i),'A') .and. lle(instr(i:i),'Z') ) then
               str(i:i)=char( ichar(instr(i:i)) + ia )
            ELSE
               str(i:i)=instr(i:i)
            ENDIF
         endif
      ENDDO
   END subroutine change_case!}}}
   !---------------------------- DIVIDER LINE -----------------------------
   subroutine find_keywords(str,ch_mark,id_key,id_value)!{{{
      !##########################################################
      !* CREATED_TIME  : 2015-05-12
      !* AUTHOR        : Xuecheng Shao
      !* CHANGE        : Xuecheng Shao
      !* ADD           : Xuecheng Shao
      !* DESCRIPTION   :              
      !     ------
      !* REFERENCES    :              
      !     ------
      !* LOG           :              
      !     2015-05-08 :              
      !* LAST MODIFIED : 2015-05-12 04:50:12 PM
      !##########################################################
      implicit none
      character(len=*)          :: str
      CHARACTER                 :: ch_mark
      INTEGER(I4B)              :: id_key
      INTEGER(I4B)              :: id_value
      INTEGER(I4B)              :: i,l_str
      !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>
      l_str=len_trim(str)
      DO i=1,l_str
         IF ( (str(i:i) == " ") .or. (str(i:i) == ch_mark) ) exit
      ENDDO
      id_key=i
      DO i=id_key,l_str
         IF ( (str(i:i) /= " ") .and. (str(i:i) /= ch_mark) ) exit
      ENDDO
      id_value=i
      id_key=id_key-1
   END subroutine find_keywords!}}}
   !---------------------------- DIVIDER LINE -----------------------------
   subroutine find_nword(str,ch_comma,nword)!{{{
      !##########################################################
      !* CREATED_TIME  : 2015-05-12
      !* AUTHOR        : Xuecheng Shao
      !* CHANGE        : Xuecheng Shao
      !* ADD           : Xuecheng Shao
      !* DESCRIPTION   :              
      !     ------
      !* REFERENCES    :              
      !     ------
      !* LOG           :              
      !     2015-05-08 :              
      !* LAST MODIFIED : 2015-05-12 04:50:12 PM
      !##########################################################
      implicit none
      character(len=*)          :: str
      CHARACTER                 :: ch_comma
      INTEGER(I4B)              :: nword
      INTEGER(I4B)              :: i,l_str
      !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>
      l_str=len_trim(str)
      nword = 1
      DO i = 1, l_str -1
         IF ( ((str(i:i) == " " ).or. (str(i:i) == ch_comma )) .and. &
            (str(i+1:i+1) /= " ") ) then
            nword = nword + 1
         ENDIF 
      ENDDO
   END subroutine !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   function norm_real(a)!{{{
      real(dp), intent(in) :: a(3)
      real(dp)             :: norm_real
      norm_real= sqrt(DOT_PRODUCT(a,a))
   end function!}}}
   !---------------------------- DIVIDER LINE -----------------------------
   function norm_complex(a)!{{{
      complex(dcp), intent(in) :: a(3)
      REAL(dp)             :: norm_complex
      norm_complex= sqrt(REAL(DOT_PRODUCT(a,a),DP))
   end function!}}}
   !---------------------------- DIVIDER LINE -----------------------------
   function cross_complex(a,b)!{{{
      complex(dcp), intent(in) :: a(3),b(3)
      complex(dcp)             :: cross_complex(3)

      cross_complex(1) = a(2)*b(3) - a(3)*b(2)
      cross_complex(2) = a(3)*b(1) - a(1)*b(3)
      cross_complex(3) = a(1)*b(2) - a(2)*b(1)
   end function!}}}
   !---------------------------- DIVIDER LINE -----------------------------
   function cross_real(a,b)!{{{
      real(dp), intent(in) :: a(3),b(3)
      real(dp)             :: cross_real(3)

      cross_real(1) = a(2)*b(3) - a(3)*b(2)
      cross_real(2) = a(3)*b(1) - a(1)*b(3)
      cross_real(3) = a(1)*b(2) - a(2)*b(1)

   end function!}}}
   !---------------------------- DIVIDER LINE -----------------------------
   FUNCTION Det(matrix)!{{{
      real(DP),  intent(in) :: Matrix(3,3)
      real(DP)              :: Det


      Det = Matrix(1,1)*(Matrix(2,2)*Matrix(3,3)-Matrix(2,3)*Matrix(3,2))&
         -Matrix(1,2)*(Matrix(2,1)*Matrix(3,3)-Matrix(2,3)*Matrix(3,1))&
         +Matrix(1,3)*(Matrix(2,1)*Matrix(3,2)-Matrix(3,1)*Matrix(2,2))

   END FUNCTION!}}}
   !---------------------------- DIVIDER LINE -----------------------------
   function inv_33(M)!{{{
   !##########################################################!{{{
   !* CREATED_TIME  : 2015-05-19
   !* AUTHOR        : Xuecheng Shao
   !* CHANGE        : Xuecheng Shao
   !* ADD           : Xuecheng Shao
   !* DESCRIPTION   :              
   !         a b c               1   ei-hf  -(bi-hc)  bf-ce
   !    mat=[d e f]  inv_33(mat)=-- [fg-id  -(cg-ia)  cd-af]
   !         g h i          det(mat) dh-ge  -(ah-gb)  ae-bd
   !* REFERENCES    :              
   !     ------
   !* LOG           :              
   !     2015-05-08 :              
   !* LAST MODIFIED : 2015-05-19 09:15:08 AM
   !##########################################################!}}}
      real(dp),intent(in) :: M(3,3)
      real(dp)            :: Inv(3,3),inv_33(3,3)
      real(dp)            ::  d
      !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>
      d = Det(M)
      inv(1,1)=M(2,2)*M(3,3) - M(2,3)*M(3,2)
      inv(2,1)=M(2,3)*M(3,1) - M(2,1)*M(3,3)
      inv(3,1)=M(2,1)*M(3,2) - M(2,2)*M(3,1)
      inv(1,2)=M(1,3)*M(3,2) - M(1,2)*M(3,3)
      inv(2,2)=M(1,1)*M(3,3) - M(1,3)*M(3,1)
      inv(3,2)=M(1,2)*M(3,1) - M(1,1)*M(3,2)
      inv(1,3)=M(1,2)*M(2,3) - M(1,3)*M(2,2)
      inv(2,3)=M(1,3)*M(2,1) - M(1,1)*M(2,3)
      inv(3,3)=M(1,1)*M(2,2) - M(1,2)*M(2,1)
      inv_33=inv/d
      !<<<<<<<<<<<<<<<<<<<<< End Body  <<<<<<<<<<<<<<<<<<<<<<<<
   end function!}}}
   !---------------------------- DIVIDER LINE -----------------------------
   FUNCTION LindG(eta,lambda,mu)!{{{
   !!!!!!!!!!!!!!!!!!!!!!!This function from PROFESS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REAL(kind=DP), INTENT(IN)      :: &
         eta, &                 ! the point at which the function gets evaluated.
         lambda, &              ! the TF multiplier for compensating.
         mu                     ! the vW multiplier

      REAL(kind=DP)                  :: &
         LindG                  ! The Lindhard G function as described in [1]

      !>> INTERNAL VARIABLES <<!
      REAL(kind=DP)                  :: &
         eta2, &                ! eta ** 2
         invEta2                ! 1/eta**2
      !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>

      IF (eta<0._DP) THEN
         LindG = 0._DP

         ! Limit for small eta
      ELSE IF (eta < 1E-10_DP) THEN
         LindG = 1._DP - lambda + eta**2 * (1._DP / 3._DP - 3._DP * mu)

         ! Around the singularity
      ELSE IF (ABS(eta-1._DP) < 1E-10_DP) THEN
         LindG = 2._DP - lambda - 3._DP * mu + 20._DP * (eta-1._DP)

         ! Taylor expansion for high eta
      ELSE IF (eta > 3.65_DP) THEN ! we determined empircally that 3.65 was a
         ! good crossover point to the taylor expansion
         eta2 = eta**2
         invEta2 = 1._DP/eta2
         LindG = 3._DP*(1._DP-mu)*eta2 &
            -lambda-0.6_DP &
            + invEta2 * (-0.13714285714285712_DP &
            + invEta2 * (-6.39999999999999875E-2_DP &
            + invEta2 * (-3.77825602968460128E-2_DP &
            + invEta2 * (-2.51824061652633074E-2_DP &
            + invEta2 * (-1.80879839616166146E-2_DP &
            + invEta2 * (-1.36715733124818332E-2_DP &
            + invEta2 * (-1.07236045520990083E-2_DP &
            + invEta2 * (-8.65192783339199453E-3_DP &
            + invEta2 * (-7.1372762502456763E-3_DP &
            + invEta2 * (-5.9945117538835746E-3_DP &
            + invEta2 * (-5.10997527675418131E-3_DP &
            + invEta2 * (-4.41060829979912465E-3_DP &
            + invEta2 * (-3.84763737842981233E-3_DP &
            + invEta2 * (-3.38745061493813488E-3_DP &
            + invEta2 * (-3.00624946457977689E-3_DP)))))))))))))))

      ELSE
         LindG = 1._DP / (0.5_DP + 0.25_DP * (1._DP-eta**2) * LOG((1._DP + eta)/&
            abs(1._DP-eta))/eta) - 3._DP * mu * eta**2 - lambda
      END IF
      !<<<<<<<<<<<<<<<<<<<<< End Body  <<<<<<<<<<<<<<<<<<<<<<<<
   END FUNCTION!}}}
   !---------------------------- DIVIDER LINE -----------------------------
   Function int_to_char( int )!{{{
      !-----------------------------------------------------------------------
      !
      INTEGER(I4B), intent(in) :: int
      character (len=6)   :: int_to_char
      if ( int <= -10 ) then
         write( unit = int_to_char , fmt = "(i3)" ) int
      else if(int == 0 ) then
         write(unit = int_to_char,fmt= "(i1)") int
      else if(int< 0 ) then
         write( unit = int_to_char , fmt = "(i2)" ) int
      else if ( int < 10 ) then
         write( unit = int_to_char , fmt = "(i1)" ) int
      else if ( int < 100 ) then
         write( unit = int_to_char , fmt = "(i2)" ) int
      else if ( int < 1000 ) then
         write( unit = int_to_char , fmt = "(i3)" ) int
      else if ( int < 10000 ) then
         write( unit = int_to_char , fmt = "(i4)" ) int
      else
         write( unit = int_to_char , fmt = "(i5)" ) int
      end if
      return
   end function!}}}
   !---------------------------- DIVIDER LINE -----------------------------
   subroutine lat2matrix(lat_para,lat_mat,flag)!{{{
      !##########################################################!{{{
      !* CREATED_TIME  : 2015-05-11 10:58:47
      !* AUTHOR        : Xuecheng Shao
      !* CHANGE        : Xuecheng Shao
      !* ADD           : Xuecheng Shao
      !* DESCRIPTION   :              
      !     ------                    
      !    [ A   B   C  ]
      !      M11 M12 M13  M13 = lA * cos(b) 
      !      0   M22 M23  M23 = ( lB*lC*cos(x) - M13*M12 ) / M22
      !      0   0   M33  
      !* REFERENCES    :              
      !     ------                    
      !* LOG           :              
      !     2015-05-11
      !* LAST MODIFIED : 2015-05-14 10:27:21 AM
      !##########################################################!}}}
      !
      USE constants,     only : DP
      !
      implicit none
      INTEGER(I4B)              :: flag
      real(DP)             :: lat_para(6)
      real(DP)             :: lat_mat(3,3)
      real(DP)             :: angle(3)
      if (flag==1) then
         lat_mat=0.0
         lat_mat(1,1) = lat_para(1)
         lat_mat(1,2) = lat_para(2)*cos(lat_para(6))
         lat_mat(2,2) = lat_para(2)*sin(lat_para(6))
         lat_mat(1,3) = lat_para(3)*cos(lat_para(5))
         lat_mat(2,3) = ( lat_para(2)*lat_para(3)*cos(lat_para(4)) - lat_mat(1,3)*lat_mat(1,2) )/ lat_mat(2,2) 
         lat_mat(3,3) = sqrt(lat_para(3)**2 - lat_mat(1,3)**2 - lat_mat(2,3)**2)
      else
         lat_para(1)=dsqrt(sum(lat_mat(:,1)**2))
         lat_para(2)=dsqrt(sum(lat_mat(:,2)**2))
         lat_para(3)=dsqrt(sum(lat_mat(:,3)**2))
         angle(1)=(dot_product(lat_mat(:,2),lat_mat(:,3))) / (lat_para(2)*lat_para(3))
         angle(2)=(dot_product(lat_mat(:,1),lat_mat(:,3))) / (lat_para(1)*lat_para(3))
         angle(3)=(dot_product(lat_mat(:,1),lat_mat(:,2))) / (lat_para(1)*lat_para(2))
         lat_para(4)=acos(angle(1))
         lat_para(5)=acos(angle(2))
         lat_para(6)=acos(angle(3))
      endif
   end subroutine!}}}
   !---------------------------- DIVIDER LINE -----------------------------
   subroutine one2three(id,n_dens,pos)!{{{
      implicit none
      INTEGER(I4B)                      :: id,n_dens(3),pos(3)
      INTEGER(I4B)                      :: ixx,iyz,iyy,izz
      ixx=mod(id,n_dens(1))
      iyz=id/n_dens(1)
      iyy=mod(iyz,n_dens(2))
      izz=iyz/n_dens(2)
      if (ixx==0) then
         pos(1)=n_dens(1)
      else
         pos(1)=ixx
      endif
      if ((ixx==0) .and. (iyy==0)) then
         pos(2)=n_dens(2)
         pos(3)=izz
      else
         pos(3)=izz+1
         if (ixx==0) then
            pos(2)=iyy
         else
            pos(2)=iyy+1
         endif
      endif
   END subroutine!}}}
   !---------------------------- DIVIDER LINE -----------------------------
   subroutine real_index(array,n,id) !{{{
      USE constants,         only : DP
      implicit none
      INTEGER(I4B)                     :: n
      real(DP)                    :: array(n)
      INTEGER(I4B)                     :: id(n)
      INTEGER(I4B)                     :: i,j,k,l,m
      real(DP),parameter          :: eps=1.d-14
      do i=1,n
         id(i)=i
      enddo
      if (n==1) return
      l=n/2+1
      k=n
      do while(.true.)
         if (l>1) then
            l=l-1
            m=id(l)
         else
            m=id(k)
            id(k)=id(1)
            k=k-1
            if (k==1) then
               id(1)=m
               return
            endif
         endif
         i=l
         j=l+1
         do while(j<=k)
            if (j<k) then
               if (array(id(j))<(array(id(j+1))+eps)) j=j+1
            endif
            if (array(m)<(array(id(j))+eps)) then
               id(i)=id(j)
               i=j
               j=j+j
            else
               j=k+1
            endif
         enddo
         id(i)=m
      enddo
   end subroutine!}}}
   !---------------------------- DIVIDER LINE -----------------------------
   subroutine integer_index(array,n,id) !{{{
      USE constants,    only : DP
      implicit none
      INTEGER(I4B)                     :: n
      INTEGER(I4B)                     :: array(n)
      INTEGER(I4B)                     :: id(n)
      INTEGER(I4B)                     :: i,j,k,l,m
      real(DP),parameter          :: eps=1.d-14
      do i=1,n
         id(i)=i
      enddo
      if (n==1) return
      l=n/2+1
      k=n
      do while(.true.)
         if (l>1) then
            l=l-1
            m=id(l)
         else
            m=id(k)
            id(k)=id(1)
            k=k-1
            if (k==1) then
               id(1)=m
               return
            endif
         endif
         i=l
         j=l+1
         do while(j<=k)
            if (j<k) then
               if (array(id(j))<(array(id(j+1))+eps)) j=j+1
            endif
            if (array(m)<(array(id(j))+eps)) then
               id(i)=id(j)
               i=j
               j=j+j
            else
               j=k+1
            endif
         enddo
         id(i)=m
      enddo
   end subroutine!}}}
   !---------------------------- DIVIDER LINE -----------------------------
   subroutine init_random_seed()!{{{
      !##########################################################
      !* CREATED_TIME  : 2015-07-20 19:38:02
      !* AUTHOR        : Xuecheng Shao
      !* CHANGE        : Xuecheng Shao
      !* Mail          : sxc@calypso.cn
      !* DESCRIPTION   :              
      !     ------                    
      !* REFERENCES    :              
      !   https://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html#RANDOM_005fSEED 
      !* LOG           :              
      !     2015-07-20
      !* LAST MODIFIED : 2015-07-20 08:16:23 PM
      !##########################################################
      use iso_fortran_env, only: int64
      ! In order to use the function GETPID 
      ! ifort need add ifport module, gfortran no need
      ! use ifport
      implicit none
      INTEGER(I4B), allocatable :: seed(:)
      INTEGER(I4B) :: i, n, un, istat, dt(8), pid
      integer(int64) :: t
      call random_seed(size = n)
      allocate(seed(n))
      ! First try if the OS provides a random number generator
      open(newunit=un, file="/dev/urandom", access="stream", &
         form="unformatted", action="read", status="old", iostat=istat)
      if (istat == 0) then
         read(un) seed
         close(un)
      else
         ! Fallback to XOR:ing the current time and pid. The PID is
         ! useful in case one launches multiple instances of the same
         ! program in parallel.
         ! pid = getpid()
         pid = 1024
         call system_clock(t)
         if (t == 0) then
            call date_and_time(values=dt)
            t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
               + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
               + dt(3) * 24_int64 * 60 * 60 * 1000 &
               + dt(5) * 60 * 60 * 1000 &
               + dt(6) * 60 * 1000 + dt(7) * 1000 &
               + dt(8)
         end if
         t = ieor(t, int(pid, kind(t)))
         do i = 1, n
            seed(i) = lcg(t)
         end do
      end if
      call random_seed(put=seed)
   contains
      ! This simple PRNG might not be good enough for real work, but is
      ! sufficient for seeding a better PRNG.
      function lcg(s)
         integer(i4b) :: lcg
         integer(int64) :: s
         if (s == 0) then
            s = 104729
         else
            s = mod(s, 4294967296_int64)
         end if
         s = mod(s * 279470273_int64, 4294967291_int64)
         lcg = int(mod(s, int(huge(0), int64)), kind(0))
      end function lcg
   end subroutine init_random_seed!}}}
   !---------------------------- DIVIDER LINE -----------------------------
   !---------------------------- DIVIDER LINE -----------------------------
   !---------------------------- DIVIDER LINE -----------------------------
   !---------------------------- DIVIDER LINE -----------------------------
   !---------------------------- DIVIDER LINE -----------------------------
   SUBROUTINE atom_mass(atom_name,mass)!{{{
      use constants,  only : DP, CONST_MA_AU
      implicit none
      !--------------------------------------------
      character(len=*)    :: atom_name
      integer             :: i
      real(DP)            :: mass
      !------------------------------------------------
      type :: elements
         integer          ::atom_number 
         character(len=3) :: symbols
         real(DP)         :: mass
      end type elements
      type(elements) ::atom(104)
      !-------------------------------------------
      !print *,"name",atom_name
      do i=1,104
         if(adjustl(atom_name)==atom(i)%symbols) then
            mass=atom(i)%mass*CONST_MA_AU
            !print *,"masssss",mass
            !print *,"name",atom_name
         endif
      enddo
      data atom(1:104) /&
         & elements(  1, "H",  1.0079470000d0),&
         & elements(  2, "He", 4.0026022000d0),&
         & elements(  3, "Li", 6.9412000000d0),&
         & elements(  4, "Be", 9.0121823000d0),&
         & elements(  5, "B",  10.811700000d0),&
         & elements(  6, "C",  12.010780000d0),&
         & elements(  7, "N",  14.006720000d0),&
         & elements(  8, "O",  15.999400000d0),&
         & elements(  9, "F",  18.998403250d0),&
         & elements( 10, "Ne", 20.179760000d0),&
         & elements( 11, "Na", 22.989769282d0),&
         & elements( 12, "Mg", 24.305060000d0),&
         & elements( 13, "Al", 26.981538680d0),&
         & elements( 14, "Si", 28.085530000d0),&
         & elements( 15, "P",  30.973762200d0),&
         & elements( 16, "S",  32.065500000d0),&
         & elements( 17, "Cl", 35.453200000d0),&
         & elements( 18, "Ar", 39.948100000d0),&
         & elements( 19, "K",  39.098310000d0),&
         & elements( 20, "Ca", 40.078400000d0),&
         & elements( 21, "Sc", 44.955912600d0),&
         & elements( 22, "Ti", 47.867100000d0),&
         & elements( 23, "V",  50.941510000d0),&
         & elements( 24, "Cr", 51.996160000d0),&
         & elements( 25, "Mn", 54.938045500d0),&
         & elements( 26, "Fe", 55.845200000d0),&
         & elements( 27, "Co", 58.933195220d0),&
         & elements( 28, "Ni", 58.693442000d0),&
         & elements( 29, "Cu", 63.546300000d0),&
         & elements( 30, "Zn", 65.382000000d0),&
         & elements( 31, "Ga", 69.723500000d0),&
         & elements( 32, "Ge", 72.641000000d0),&
         & elements( 33, "As", 74.921595600d0),&
         & elements( 34, "Se", 78.963000000d0),&
         & elements( 35, "Br", 79.904100000d0),&
         & elements( 36, "Kr", 83.798200000d0),&
         & elements( 37, "Rb", 85.467830000d0),&
         & elements( 38, "Sr", 87.621000000d0),&
         & elements( 39, "Y",  88.905852000d0),&
         & elements( 40, "Zr", 91.224000000d0),&
         & elements( 41, "Nb", 92.906382000d0),&
         & elements( 42, "Mo", 95.962000000d0),&
         & elements( 43, "Tc", 98.000000000d0),&
         & elements( 44, "Ru", 101.07200000d0),&
         & elements( 45, "Rh", 102.90550300d0),&
         & elements( 46, "Pd", 106.42100000d0),&
         & elements( 47, "Ag", 107.86822000d0),&
         & elements( 48, "Cd", 112.41180000d0),&
         & elements( 49, "In", 114.81830000d0),&
         & elements( 50, "Sn", 118.71070000d0),&
         & elements( 51, "Sb", 121.76010000d0),&
         & elements( 52, "Te", 127.60300000d0),&
         & elements( 53, "I",  126.90447300d0),&
         & elements( 54, "Xe", 131.29360000d0),&
         & elements( 55, "Cs", 132.90545192d0),&
         & elements( 56, "Ba", 137.32770000d0),&
         & elements( 57, "La", 138.90547700d0),&
         & elements( 58, "Ce", 140.11610000d0),&
         & elements( 59, "Pr", 140.90765200d0),&
         & elements( 60, "Nd", 144.24230000d0),&
         & elements( 61, "Pm", 145.00000000d0),&
         & elements( 62, "Sm", 150.36200000d0),&
         & elements( 63, "Eu", 151.96410000d0),&
         & elements( 64, "Gd", 157.25300000d0),&
         & elements( 65, "Tb", 158.92535200d0),&
         & elements( 66, "Dy", 162.50010000d0),&
         & elements( 67, "Ho", 164.93032300d0),&
         & elements( 68, "Er", 167.25930000d0),&
         & elements( 69, "Tm", 168.93421200d0),&
         & elements( 70, "Yb", 173.05450000d0),&
         & elements( 71, "Lu", 174.96681000d0),&
         & elements( 72, "Hf", 178.49200000d0),&
         & elements( 73, "Ta", 180.94788200d0),&
         & elements( 74, "W",  183.84100000d0),&
         & elements( 75, "Re", 186.20710000d0),&
         & elements( 76, "Os", 190.23000000d0),&
         & elements( 77, "Ir", 192.21730000d0),&
         & elements( 78, "Pt", 195.08490000d0),&
         & elements( 79, "Au", 196.96656940d0),&
         & elements( 80, "Hg", 200.59200000d0),&
         & elements( 81, "Ti", 204.38332000d0),&
         & elements( 82, "Pb", 207.21000000d0),&
         & elements( 83, "Bi", 208.98040100d0),&
         & elements( 84, "Po", 209.00000000d0),&
         & elements( 85, "At", 210.00000000d0),&
         & elements( 86, "Rn", 222.00000000d0),&
         & elements( 87, "Fr", 223.00000000d0),&
         & elements( 88, "Ra", 226.00000000d0),&
         & elements( 89, "Ac", 227.00000000d0),&
         & elements( 90, "Th", 232.03806200d0),&
         & elements( 91, "Pa", 231.03588200d0),&
         & elements( 92, "U",  238.02891300d0),&
         & elements( 93, "Np", 237.00000000d0),&
         & elements( 94, "Pu", 244.00000000d0),&
         & elements( 95, "Am", 243.00000000d0),&
         & elements( 96, "Cm", 247.00000000d0),&
         & elements( 97, "Bk", 247.00000000d0),&
         & elements( 98, "Cf", 251.00000000d0),&
         & elements( 99, "Es", 252.00000000d0),&
         & elements(100, "Fm", 257.00000000d0),&
         & elements(101, "Md", 258.00000000d0),&
         & elements(102, "No", 259.00000000d0),&
         & elements(103, "Lr", 262.00000000d0),& 
         & elements(104, "Rf", 265.00000000d0)/
      ! & elements(105, "Db",268.000000)
      ! & elements(106, "Sg",271.000000)
      ! & elements(107, "Bh",272.000000)
      ! & elements(108, "Hs",277.000000)
      ! & elements(109, "Mt",276.000000)
      ! & elements(110, "Ds",281.000000)
      ! & elements(111, "Rg",280.000000)
      ! & elements(112, "Cn",285.000000)
      ! & elements(113, "Uut",284.000000)
      ! & elements(114, "Fl",289.000000)
      ! & elements(115, "Uup",288.000000)
      ! & elements(116, "Lv",293.000000)
      ! & elements(117, "Uus",294.000000)
      ! & elements(118, "Uuo",294.000000) /
      !--------------------------------------------------
   END SUBROUTINE !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   subroutine newton_inter(n,x,y,m,tx,ty)!{{{
      !       《数值计算方法》
      !       林成森  科学出版社
      !----------------------------------------------------
      INTEGER(i4b) :: n
      INTEGER :: m
      REAL(DP)  :: tx(m),ty(m)
      REAL(DP)  :: x(0:n-1),y(0:n-1)
      !
      INTEGER(i4b) :: i,j,k
      REAL(DP)  :: qx(n)
      REAL(DP)  :: Q(0:n-1,0:n-1)
      q(:,0)=y(:)
      do i=1,n-1
         do j=1,i
            Q(i,j)=(Q(i,j-1)-Q(i-1,j-1))/(x(i)-x(i-j))
         end do
      end do
      qx(n)=Q(n-1,n-1)
      DO i = 1, m
         do k=n-1,1,-1
            qx(k)=Q(k-1,k-1)+qx(k+1)*(tx(i)-x(k-1))
         end do
         ty(i)=qx(1)
      ENDDO 
   end subroutine !}}}
   !---------------------------- DIVIDER LINE -----------------------------
SUBROUTINE diag(N, inA, W, Q)!{{{
   ! ----------------------------------------------------------------------------
   ! Numerical diagonalization of nxn matrcies
   ! Copyright (C) 2006  Joachim Kopp
   ! ----------------------------------------------------------------------------
   ! This library is free software; you can redistribute it and/or
   ! modify it under the terms of the GNU Lesser General Public
   ! License as published by the Free Software Foundation; either
   ! version 2.1 of the License, or (at your option) any later version.

   ! This library is distributed in the hope that it will be useful,
   ! but WITHOUT ANY WARRANTY; without even the implied warranty of
   ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   ! Lesser General Public License for more details.

   ! You should have received a copy of the GNU Lesser General Public
   ! License along with this library; if not, write to the Free Software
   ! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
   ! ----------------------------------------------------------------------------
   ! ----------------------------------------------------------------------------
   ! Calculates the eigenvalues and normalized eigenvectors of a symmetric nxn
   ! matrix A using the Jacobi algorithm.
   ! The upper triangular part of A is destroyed during the calculation,
   ! the diagonal elements are read but not destroyed, and the lower
   ! triangular elements are not referenced at all.
   ! ----------------------------------------------------------------------------
   ! Parameters:
   !   A: The symmetric input matrix
   !   Q: Storage buffer for eigenvectors
   !   W: Storage buffer for eigenvalues
   ! ----------------------------------------------------------------------------
   !     .. Arguments ..
   implicit none
   !     .. Parameters ..
   INTEGER(I4B)                     :: N
   real(DP) :: inA(N,N)
   real(DP) :: A(N,N)
   real(DP) :: Q(N,N)
   real(DP) :: W(N)
   !     .. Local Variables ..
   real(DP) :: SD, SO
   real(DP) :: S, C, T
   real(DP) :: G, H, Z, THETA
   real(DP) :: THRESH
   INTEGER(I4B)                     :: I, X, Y, R

   !     Initialize Q to the identitity matrix
   !     --- This loop can be omitted if only the eigenvalues are desired ---
   A = inA
   !WRITE(*,*)"AA"
   !WRITE(*,*)A
   DO X = 1, N
      Q(X,X) = 1.0D0
      DO Y = 1, X-1
         Q(X, Y) = 0.0D0
         Q(Y, X) = 0.0D0
      END DO
   END DO

   !     Initialize W to diag(A)
   DO X = 1, N
      W(X) = A(X, X)
   END DO

   !     Calculate SQR(tr(A))
   SD = 0.0D0
   DO X = 1, N
      SD = SD + ABS(W(X))
   END DO
   SD = SD**2

   !     Main iteration loop
   DO I = 1, 50
      !       Test for convergence
      SO = 0.0D0
      DO X = 1, N
         DO Y = X+1, N
            SO = SO + ABS(A(X, Y))
         END DO
      END DO
      IF (SO == 0.0D0) THEN
         RETURN
      END IF

      IF (I < 4) THEN
         THRESH = 0.2D0 * SO / N**2
      ELSE
         THRESH = 0.0D0
      END IF

      !       Do sweep
      DO X = 1, N
         DO Y = X+1, N
            G = 100.0D0 * ( ABS(A(X, Y)) )
            IF ( (I > 4) .AND. (ABS(W(X)) + G == ABS(W(X)) ) &
               .AND. (ABS(W(Y)) + G == ABS(W(Y))) ) THEN
               A(X, Y) = 0.0D0
            ELSEIF (ABS(A(X, Y)) > THRESH) THEN
               !             Calculate Jacobi transformation
               H = W(Y) - W(X)
               IF ( (ABS(H) + G) == ABS(H) ) THEN
                  T = A(X, Y) / H
               ELSE
                  THETA = 0.5D0 * H / A(X, Y)
                  IF (THETA < 0.0D0) THEN
                     T = -1.0D0 / (SQRT(1.0D0 + THETA**2) - THETA)
                  ELSE
                     T = 1.0D0 / (SQRT(1.0D0 + THETA**2) + THETA)
                  END IF
               END IF

               C = 1.0D0 / SQRT( 1.0D0 + T**2 )
               S = T * C
               Z = T * A(X, Y)

               !             Apply Jacobi transformation
               A(X, Y) = 0.0D0
               W(X)    = W(X) - Z
               W(Y)    = W(Y) + Z
               DO R = 1, X-1
                  T       = A(R, X)
                  A(R, X) = C * T - S * A(R, Y)
                  A(R, Y) = S * T + C * A(R, Y)
               END DO
               DO R = X+1, Y-1
                  T       = A(X, R)
                  A(X, R) = C * T - S * A(R, Y)
                  A(R, Y) = S * T + C * A(R, Y)
               END DO
               DO R = Y+1, N
                  T       = A(X, R)
                  A(X, R) = C * T - S * A(Y, R)
                  A(Y, R) = S * T + C * A(Y, R)
               END DO
               !Update eigenvectors
               !--- This loop can be omitted if only the eigenvalues are desired ---
               DO R = 1, N
                  T       = Q(R, X)
                  Q(R, X) = C * T - S * Q(R, Y)
                  Q(R, Y) = S * T + C * Q(R, Y)
               END DO
            ENDIF
         END DO
      END DO
   END DO

   PRINT *, "No convergence."

END SUBROUTINE!}}}
   !-----------------------------------------------------------------------
function boltzmann_distribution(rnull,width)!{{{
   real(DP) :: num1,num2
   real(DP) :: boltzmann_distribution,width,rnull
   REAL(DP),parameter       :: tpi = 8.d0*atan(1.d0)
   !real(q),parameter :: twopi = 6.283185307179586_q 

   call random_number(num1)
   num2=0.d0
   do
      call random_number(num2)
      num2=abs(num2)
      if (num2 .gt. 1e-08) exit
   enddo
   !write(*,*) 'num1,num2',num1,num2
   boltzmann_distribution= cos( tpi*num1 ) * sqrt( 2.d0 *abs(log(num2)) )
   !write(*,*) 'boltzmann_distribution_',boltzmann_distribution,width,rnull
   boltzmann_distribution = width * boltzmann_distribution  +  rnull
   !write(*,*) 'boltzmann_distribution',boltzmann_distribution
END FUNCTION !}}}
!>>>>ADD by Qiang Xu
!#####################################################!
!cry_cyyoordinates to ort_coordinates      
!#####################################################!
SUBROUTINE dir2car(cry_coo,ort_coo,lat)
IMPLICIT NONE
INTEGER :: i,IDEM
REAL(DP),DIMENSION(:,:) :: cry_coo,ort_coo
REAL(DP) :: lat(3,3)
!>>>>>>>>>>>>>>>>>>>>>Main body>>>>>>>>>>>>>>>>>>>>>>>
   IDEM = size(cry_coo,2)
   DO i = 1,IDEM
      ort_coo(:,i) =  matmul(lat,cry_coo(:,i))
   ENDDO
!<<<<<<<<<<<<<<<<<<<<<End body<<<<<<<<<<<<<<<<<<<<<<<<
END SUBROUTINE
!#####################################################!
!ort_coordinates to cry_coordinates
!#####################################################!
!----------------------ort to cry------------------------      
SUBROUTINE car2dir(ort_coo,cry_coo,lat)
IMPLICIT NONE
INTEGER :: i,IDEM
REAL(DP),DIMENSION(:,:) :: cry_coo,ort_coo
REAL(DP):: inv_lattice(3,3),lat(3,3)
!>>>>>>>>>>>>>>>>>>>>>Main body>>>>>>>>>>>>>>>>>>>>>>>>
   IDEM = size(ort_coo,2)
   inv_lattice =inv_33(lat)
   DO i=1,IDEM
      cry_coo(:,i)= matmul(inv_lattice,ort_coo(:,i))
   ENDDO
!<<<<<<<<<<<<<<<<<<<<<End body<<<<<<<<<<<<<<<<<<<<<<<<<
END SUBROUTINE
!----------------------3D-to-matrix--------------------------
SUBROUTINE thr2mat(n1,n2,n3,I,J,K,dimnu)
   IMPLICIT NONE
   INTEGER(I4B),INTENT(IN) :: I,J,K,n1,n2,n3
   INTEGER(I4B),INTENT(OUT) :: dimnu
   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   dimnu=(K-1)*n2*n1+(J-1)*n1+I
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDSUBROUTINE thr2mat
!---------------------matrix-to-3D---------------------------
SUBROUTINE mat2thr(n1,n2,n3,I,Ix,Iy,Iz,offset)
   IMPLICIT NONE
   INTEGER(I4B),INTENT(IN) :: n1,n2,n3,I
   INTEGER(I4B),INTENT(OUT) :: Ix,Iy,Iz
   INTEGER(I4B),INTENT(IN),OPTIONAL :: offset(3)
   !
   INTEGER(I4B) :: timz,timy,delta,deltad
   INTEGER(I4B) :: n21
   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   n21=n2*n1
   !Iz
   timz=I/n21
   !delta=I-timz*n21
   delta=MOD(I,n21)
   IF(delta==0)THEN
      Iz=timz
      Iy=n2
      Ix=n1
   ELSE
      Iz=timz+1
      !------------------------
      timy=delta/n1
      !deltad=delta-timy*n1
      deltad=MOD(delta,n1)
      IF(deltad==0)THEN
         Iy=timy
         Ix=n1
      ELSE
         Iy=timy+1
         Ix=deltad
      ENDIF
   ENDIF
   !
   IF(PRESENT(offset))THEN
      Ix=Ix-offset(1)
      Iy=Iy-offset(2)
      Iz=Iz-offset(3)
   ENDIF
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDSUBROUTINE mat2thr
!------------------------------------------------------------
   SUBROUTINE SOPO(A,LDA,N,B,LDB,M,W)
     !
      IMPLICIT NONE
 
      INTEGER, INTENT(IN) :: & 
         N, M, LDA, LDB 
      REAL(KIND=8), DIMENSION(:,:), INTENT(INOUT) :: & 
         A
      REAL(KIND=8), DIMENSION(:,:), INTENT(INOUT) :: &
         B
      REAL(KIND=8), DIMENSION(:), INTENT(INOUT) :: &
         W
     !C---------------------------------------------------------------
     !C      SOLUTION OF A SET OF LINEAR EQUATIONS WITH A POSITIVE
     !C      DEFINITE REAL SYMMETRIC MATRIX AND WITH MULTIPLE R.H.S.
     !C---------------------------------------------------------------

      INTEGER :: &
         I, J, L, MQ, Ntemp

      REAL(KIND=8) :: &
         DUM
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      Ntemp=LDA+LDB
!C                    ROZKLAD
      DO 310 L=1,N-1
      A(L,L)=SQRT(A(L,L))
      DUM=1.d0/A(L,L)
      DO 311 I=L+1,N
      W(I)=A(I,L)*DUM
      A(I,L)=W(I)
311   CONTINUE
      DO 312 J=L+1,N
      DO 313 I=J,N
      A(I,J)=A(I,J)-W(I)*W(J)
313   CONTINUE
312   CONTINUE
310   CONTINUE
      A(N,N)=SQRT(A(N,N))

!C                   CYKLUS PRES PRAVE STRANY
      DO 320 MQ=1,M
!C                   INVERZE DOLNI TROJUH. MATICE
      DO 321 L=1,N-1
      W(L)=B(L,MQ)/A(L,L)
      DO 322 I=L+1,N
      B(I,MQ)=B(I,MQ)-A(I,L)*W(L)
322   CONTINUE
321   CONTINUE
      W(N)=B(N,MQ)/A(N,N)
!C                   INVERZE HORNI TROJUH. MATICE
      B(N,MQ)=W(N)/A(N,N)
      DO 331 L=N-1,1,-1
      DO 332 I=N,L+1,-1
      W(L)=W(L)-B(I,MQ)*A(I,L)
332   CONTINUE
      B(L,MQ)=W(L)/A(L,L)
331   CONTINUE

320   CONTINUE
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE SOPO
!------------------------------------------------------------
   SUBROUTINE csort_eigen(nev,arr,brr)
   !#################################################################!
   ! For a given set of eigenvalues (for all representations), sort  !
   ! them in ascending order. Use straight insertion, no need for    !
   ! fancy algorithms since the arrays are always short.             !
   ! ISBN 7-03-010217-7                                              !
   !#################################################################!
      !
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN) :: nev
      REAL(DP),INTENT(INOUT) :: arr(:)
      COMPLEX(DCP),INTENT(INOUT) :: brr(:,:)
      !
      REAL(DP) :: tmpa
      COMPLEX(DCP) :: tmpb(SIZE(brr,1))
      INTEGER(I4B) :: I,J
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      DO J=2,nev
!print*,'J',J
         tmpa=arr(J)
         tmpb(:)=brr(:,J)
         DO I=J-1,1,-1
            IF(arr(I)<=tmpa) EXIT
!print*,'I+1,I',I+1,I
            arr(I+1)=arr(I)
            brr(:,I+1)=brr(:,I)
            IF(I==1) EXIT
         ENDDO

!print*,'I+1,I',I+1,I
!pause
         IF(arr(I)<=tmpa)THEN
            arr(I+1)=tmpa
         ELSE
            I=0
            arr(I+1)=tmpa
         ENDIF
         brr(:,I+1)=tmpb(:)
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE csort_eigen
!--------------------------------------------------------------------------
   SUBROUTINE rsort_eigen(nev,arr,brr)
   !#################################################################!
   ! For a given set of eigenvalues (for all representations), sort  !
   ! them in ascending order. Use straight insertion, no need for    !
   ! fancy algorithms since the arrays are always short.             !
   ! ISBN 7-03-010217-7                                              !
   !#################################################################!
      !
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN) :: nev
      REAL(DP),INTENT(INOUT) :: arr(:)
      REAL(DP),INTENT(INOUT) :: brr(:,:)
      !
      REAL(DP) :: tmpa
      REAL(DP) :: tmpb(SIZE(brr,1))
      INTEGER(I4B) :: I,J
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      DO J=2,nev
         tmpa=arr(J)
         tmpb(:)=brr(:,J)
         DO I=J-1,1,-1
            IF(arr(I)<=tmpa) EXIT
            arr(I+1)=arr(I)
            brr(:,I+1)=brr(:,I)
            IF(I==1) EXIT
         ENDDO

         IF(arr(I)<=tmpa)THEN
            arr(I+1)=tmpa
         ELSE
            I=0
            arr(I+1)=tmpa
         ENDIF
         brr(:,I+1)=tmpb(:)
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE rsort_eigen
!--------------------------------------------------------------------------
   SUBROUTINE realInt_sort(nev,arr,brr,crr)
   !#################################################################!
   ! For a given set of eigenvalues (for all representations), sort  !
   ! them in ascending order. Use straight insertion, no need for    !
   ! fancy algorithms since the arrays are always short.             !
   ! ISBN 7-03-010217-7                                              !
   !#################################################################!
      !
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN) :: nev
      REAL(DP),INTENT(INOUT) :: arr(:)
      INTEGER(I4B),INTENT(INOUT) :: brr(:,:)
      REAL(DP),OPTIONAL,INTENT(INOUT) :: crr(3,nev)
      !
      REAL(DP) :: tmpa
      INTEGER(I4B) :: tmpb(SIZE(brr,1))
      REAL(DP) :: tmpc(3)
      INTEGER(I4B) :: I,J
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      DO J=2,nev
         tmpa=arr(J)
         tmpb(:)=brr(:,J)
         IF(PRESENT(crr)) tmpc(:)=crr(:,J)
         DO I=J-1,1,-1
            IF(arr(I)<=tmpa) EXIT
            arr(I+1)=arr(I)
            brr(:,I+1)=brr(:,I)
            IF(PRESENT(crr)) crr(:,I+1)=crr(:,I)
            IF(I==1) EXIT
         ENDDO

         IF(arr(I)<=tmpa)THEN
            arr(I+1)=tmpa
         ELSE
            I=0
            arr(I+1)=tmpa
         ENDIF
         brr(:,I+1)=tmpb(:)
         IF(PRESENT(crr)) crr(:,I+1)=tmpc(:)
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE realInt_sort
!-----------------------------------------------------------------------
   SUBROUTINE sort_eigval(n,arr)
      !
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN) :: n
      REAL(DP),INTENT(INOUT)  :: arr(n)
      !LOCAL
      INTEGER(I4B) :: i,j
      REAL(DP) :: a
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      DO j=2,n

         a=arr(j)

         DO i=j-1,1,-1
            IF(arr(i)<=a) EXIT
            arr(i+1)=arr(i)
         ENDDO
 
         IF(arr(i)<=a)THEN
            arr(i+1)=a
         ELSE
            i=0
            arr(i+1)=a
         ENDIF

      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE sort_eigval
!-------------------point group function operator-----------------------
SUBROUTINE PGFO(Omat,n1,n2,n3,Ix,Iy,Iz,Ex,Ey,Ez)
   !
   IMPLICIT NONE
   !INOUT
   REAL(DP),INTENT(IN) :: Omat(3,3)
   INTEGER(I4B),INTENT(IN)  :: n1,n2,n3  & !number of grids in real space
                &    ,         Ix,Iy,Iz    !The point we like
   INTEGER(I4B),INTENT(OUT) :: Ex,Ey,Ez    !the point to read in origin fun
   !LOCAL
   REAL(DP) :: inv_Omat(3,3)
   REAL(DP) :: rin(3),rout(3)
   INTEGER(I4B) :: temp(3)
   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !real for matrix product
   rin(1)=REAL(Ix,DP)
   rin(2)=REAL(Ix,DP)
   rin(3)=REAL(Ix,DP)
   !inver Omat
   inv_Omat(:,:)=inv_33( Omat(:,:) )
   !operator
   rout=MATMUL(inv_Omat,rin)
   !to grid mesh
   temp(:)=NINT(rout(:))
   !regularlizetion
   Ex=MODULO(temp(1),n1)
   IF(Ex==0) Ex=1

   Ey=MODULO(temp(2),n2)
   IF(Ey==0) Ey=1

   Ez=MODULO(temp(3),n3)
   IF(Ez==0) Ez=1
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDSUBROUTINE PGFO
!-----------------------CubicsplineInterpolation-------------------------
FUNCTION CubicSplineInterp(fun,ddfdx2,xmax,dx,x,Zion)
   !just for uniform grid now , x>0
   !IN/OUT
   REAL(DP),INTENT(IN) :: fun(:)  &
                       &, ddfdx2(:) &
                       &,xmax , dx , x
   REAL(DP),OPTIONAL,INTENT(IN)  :: Zion
   REAL(DP) :: CubicSplineInterp
   !LOCAL
   REAL(DP) :: yval  &
           &, ypp_left &
           &, ypp_right &
           &, y_left  &
           &, y_right &
           &, pos &
           &, dt
   INTEGER(I4B) :: left  &
               &,  right
   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !
   IF(x<0.d0)THEN
      WRITE(*,*) 'CubicSplineInterp: r must > 0'
      STOP
   ENDIF
   !Some special case
   IF(x==0.d0)THEN
      CubicSplineInterp=fun(1)
      RETURN
   ENDIF
   !For r > rmax , fun=0.d0
   IF(x>xmax)THEN
      CubicSplineInterp=0.d0
      RETURN
   ENDIF
   !vlpp
   IF(PRESENT(Zion))THEN
      ! read the following words to see why I do this check here
      ! since NaN is not equal to anyting in FORTRAN
      ! we use the following to see if 4*pi/qNorm is too big to case a NaN
      IF (-4._DP*pi/x**2_DP .NE. -4_DP*pi/x**2_DP ) THEN
         WRITE(*,*)&
            ' There is another case which should be considered: very large bulk. ', &
            '  because larger system will give denser q points and will put q points ', &
            ' very closer to q=0 and this will make -4*pi*Z/q**2 to -infinity',&
            ' If computer find -4*pi*Z/q**2 is too big, it will generate NaN',&
            ' But currently, in our group, we have not meet a system large enough to ',&
            ' cause CPU to generate NaN, if you meet such kind of problem, do something!, code STOP!!'
         STOP
      ENDIF
   ENDIF
   !interpolation:
   !1.Find the left and right point , dt
   pos=x / dx + 1.d0
   !left
   left=FLOOR(pos)
   !right
   right=left+1
   !dt
   dt=pos-REAL(left,DP)
   !2.Evaulate the polynomial
   y_left=fun(left)
   y_right=fun(right)
   ypp_left=ddfdx2(left)
   ypp_right=ddfdx2(right)
   !just do it
   yval = y_left + dt * ( y_right - y_left &
      - ( ypp_right / 6.0_DP + ypp_left / 3.0_DP ) &
      + dt * ( 0.5_DP * ypp_left &
      + dt * ( ( ypp_right - ypp_left ) / 6.0_DP ) ) )
   !output
   IF(PRESENT(Zion))THEN
      CubicSplineInterp=yval - 4.d0*pi*REAL(Zion,DP) / x**2._DP
   ELSE
      CubicSplineInterp=yval
   ENDIF

   RETURN
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDFUNCTION CubicSplineInterp
!------------------------------------------------------------------------
   subroutine finite_factor(fnor,norder,coe)!{{{
      !##########################################################
      !* CREATED_TIME  : 2013-03-12
      !* AUTHOR        : Yanchao Wang
      !* CHANGE        : Xuecheng Shao
      !* ADD           : 
      !* DESCRIPTION   :              
      !     ------
      !* REFERENCES    :              
      !     ------
      !* LOG           :              
      !     2015-05-08 :              
      !* LAST MODIFIED : 2015-07-24 05:31:05 PM
      !##########################################################
      use constants
      !
      implicit none
      INTEGER(I4B), intent(in)   :: norder
      INTEGER(I4B), intent(in)   :: fnor
      real(dp), intent(out) :: coe(-norder:norder)
      INTEGER(I4B)               :: i
      !
      !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>
      if(fnor.eq.1) then
         select case (norder)
         case (1)
            coe(1) =  0.50000000000000D+00
         case (2)
            coe(1:2) = [ &
               0.666666666666667      , -8.333333333333333E-002 ]
         case (3)
            coe(1:3) = [ &
               0.750000000000000      , -0.150000000000000      , &
               1.666666666666667E-002 ]
         case (4)
            coe(1:4) = [ &
               0.800000000000000      , -0.200000000000000      , &
               3.809523809523809E-002 , -3.571428571428571E-003 ]
         case (5)
            coe(1:5) = [ &
               0.833333333333333      , -0.238095238095238      , &
               5.952380952380952E-002 , -9.920634920634920E-003 , &
               7.936507936507937E-004 ]
         case (6)
            coe(1:6) = [ &
               0.857142857142857      , -0.267857142857143      , &
               7.936507936507936E-002 , -1.785714285714286E-002 , &
               2.597402597402598E-003 , -1.803751803751804E-004 ]
         case (7)
            coe(1:7) = [ &
               0.875000000000000      , -0.291666666666667      , &
               9.722222222222224E-002 , -2.651515151515151E-002 , &
               5.303030303030303E-003 , -6.798756798756799E-004 , &
               4.162504162504163E-005 ]
         case (8)
            coe(1:8) = [ &
               0.888888888888889      , -0.311111111111111      , &
               0.113131313131313      , -3.535353535353535E-002 , &
               8.702408702408702E-003 , -1.554001554001554E-003 , &
               1.776001776001776E-004 , -9.712509712509713E-006]
         case (9)
            coe(1:9) = [ &
               0.900000000000000      , -0.327272727272727      , &
               0.127272727272727      , -4.405594405594405E-002 , &
               1.258741258741259E-002 , -2.797202797202797E-003 , &
               4.495504495504496E-004 , -4.627725215960510E-005 , &
               2.285296402943462E-006 ]
         case (10)
            coe(1:10) = [ &
               0.909090909090909      , -0.340909090909091      , &
               0.139860139860140      , -5.244755244755244E-002 , &
               1.678321678321678E-002 , -4.370629370629371E-003 , &
               8.814714697067639E-004 , -1.285479226655697E-004 , &
               1.202787580496559E-005 , -5.412544112234515E-007 ]
         end select
         coe(0) = 0.d0
         do i=1,norder
            coe(-i) =-coe(i)
         enddo
      else if (fnor.eq.2) then

         select case (norder)
         case (1)
            coe(0:1) = [ &
               -2.00000000000000, 1.00000000000000]
         case (2)
            coe(0:2) = [ &
               -2.50000000000000, 1.33333333333333, &
               -8.333333333333333E-002]
         case (3)
            coe(0:3) = [ &
               -2.72222222222222  , 1.50000000000000       , &
               -0.150000000000000 , 1.111111111111111E-002 ]
         case (4)
            coe(0:4) = [&
               -2.84722222222222  , 1.60000000000000       , &
               -0.200000000000000 , 2.539682539682540E-002 , &
               -1.785714285714286E-003 ]
         case (5)
            coe(0:5) = [ &
               -2.92722222222222       , 1.66666666666667       , &
               -0.238095238095238      , 3.968253968253968E-002 , &
               -4.960317460317460E-003 , 3.174603174603175E-004 ]
         case (6)
            coe(0:6) = [ &
               -2.98277777777778       , 1.71428571428571       , &
               -0.267857142857143      , 5.291005291005291E-002 , &
               -8.928571428571428E-003 , 1.038961038961039E-003 , &
               -6.012506012506013E-005 ]
         case (7)
            coe(0:7) = [ &
               -3.02359410430839       , 1.75000000000000       , &
               -0.291666666666667      , 6.481481481481481E-002 , &
               -1.325757575757576E-002 , 2.121212121212121E-003 , &
               -2.266252266252267E-004 , 1.189286903572618E-005 ]
         case (8)
            coe(0:8) = [ &
               -3.05484410430839       , 1.77777777777778       , &
               -0.311111111111111      , 7.542087542087542E-002 , &
               -1.767676767676768E-002 , 3.480963480963481E-003 , &
               -5.180005180005181E-004 , 5.074290788576503E-005 , &
               -2.428127428127428E-006 ]
         case (9)
            coe(0:9) = [ &
               -3.07953546233308       , 1.80000000000000       , &
               -0.327272727272727      , 8.484848484848484E-002 , &
               -2.202797202797203E-002 , 5.034965034965034E-003 , &
               -9.324009324009329E-004 , 1.284429855858427E-004 , &
               -1.156931303990128E-005 , 5.078436450985472E-007 ]
         case (10)
            coe(0:10) = [ &
               -3.09953546233308       , 1.81818181818182       , &
               -0.340909090909091      , 9.324009324009326E-002 , &
               -2.622377622377622E-002 , 6.713286713286712E-003 , &
               -1.456876456876457E-003 , 2.518489913447896E-004 , &
               -3.213698066639244E-005 , 2.672861289992354E-006 , &
               -1.082508822446903E-007 ]
         end select
         do i=1,norder
            coe(-i) = coe(i)
         enddo
      endif
      !<<<<<<<<<<<<<<<<<<<<< End Body  <<<<<<<<<<<<<<<<<<<<<<<<
   end subroutine finite_factor
!----------------------finite_factor-------------------------
   SUBROUTINE finite_factor_new(fnor,norder,coe)
      !##########################################################
      !* CREATED_TIME  : 2018-09-03
      !* AUTHOR        : Qiang Xu
      !* Ref.          : 'Real-space grid representation of
      ! momentum and kinetic energy operators for energy operators
      ! for electronic structure calculations'
      !*               : Physics meaning
      !##########################################################
      use constants
      !
      implicit none
      INTEGER(I4B), intent(in)   :: norder
      INTEGER(I4B), intent(in)   :: fnor
      real(dp), intent(out) :: coe(-norder:norder)
      INTEGER(I4B)               :: i
      REAL(DP) :: Omega(norder)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CALL OmegaMm(Omega)
      IF(fnor==1)THEN
          !Gradient operator
          coe(0)=0._DP
          DO i=1,norder
             coe(i)=1._DP/(2._DP*i*Omega(i))
             coe(-i)=-coe(i)
          ENDDO
      ELSEIF(fnor==2)THEN
          !Laplace
          !coe(i==0)
          coe(0)=0._DP
          DO i=1,norder
             coe(0)=coe(0)-2._DP/(i**2*Omega(i))
          ENDDO
          !coe(i/=0)
          DO i=1,norder
             coe(i)=1._DP/(i**2*Omega(i))
             coe(-i)=coe(i)
          ENDDO
      ENDIF
      !
      CONTAINS
      !Omega function
      SUBROUTINE OmegaMm(Omeg)
         IMPLICIT NONE
         REAL(DP),INTENT(OUT) :: Omeg(norder)
         INTEGER(I4B) :: m,l
         REAL(DP) :: tmp
         !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
         DO m=1,norder
            tmp=1._DP
            DO l=1,norder
               IF(l/=m) tmp=tmp*(1._DP-(REAL(m,DP)/l)**2)
            ENDDO 
            Omeg(m)=tmp
         ENDDO
         !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ENDSUBROUTINE OmegaMm
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE finite_factor_new
!--------------------------dfdr--------------------------------
   SUBROUTINE dfdr(np,h,f,df)
      !just for first order of the don-dimension fun,r : uniform grid
      USE parameters , ONLY : finite_order=>nfd
      IMPLICIT NONE
      !INOUT
      INTEGER(I4B),INTENT(IN) :: np
      REAL(DP),INTENT(IN)  :: f(np)  &   !f
                        & , h  !grid size
      REAL(DP),INTENT(OUT) :: df(np)  !df_dr
      !LOCAL
      REAL(DP) :: coe(-finite_order:finite_order)
      REAL(DP) :: ft(-finite_order:np),tmp
      INTEGER(I4B) :: i,ish
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ft(1:np)=f(1:np)
      !get coe
      CALL finite_factor(1,finite_order,coe)
      coe(:)=coe(:)/h
      !center symmetry of origin set now
      DO i=-finite_order,0
         ft(i)=f(2-i)
      ENDDO
      !finite difference
      DO i=finite_order,np-finite_order
         tmp=0.d0
         DO ish=-finite_order,finite_order,1
            tmp=tmp + coe(ish)*ft(ish+i)
         ENDDO
         df(i) = tmp
      ENDDO
      !first few point
      df(1:finite_order-1)=0.d0
      !Last few points
      df(np-finite_order+1:np)=0.d0
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<     
   ENDSUBROUTINE dfdr
!--------------------interploate value-------------------------
FUNCTION CubicHermiteInterp(fun,dfdx,xmax,h,x)
   !##############################################!
   !y=A0*y0+A1*y1+B0*dy0+B1*dy1                   !
   !defined : dx0=x-x0 and dx1=x-x1 and h=x1-x0   !
   !A0=( 1+2*dx0/h )*(dx1/h)**2                   !
   !A1=( 1-2*dx1/h )*(dx0/h)**2                   !
   !B0=dx0*(dx1/h)**2                             !
   !B1=dx1*(dx0/h)**2                             !
   !##############################################!
   IMPLICIT NONE
   REAL(DP),INTENT(IN) :: fun(:)   &  !in fun
                      &,  dfdx(:)  &  !first derivertive
                      &,  xmax     &  !max x
                      &,  h        &  !grid size
                      &,  x           !in x
   REAL(DP) :: CubicHermiteInterp
   !LOCAL
   REAL(DP) :: pos,y0,y1,dy0,dy1
   REAL(DP) :: dx0,dx1,dx0_h,dx1_h,dx0_h2,dx1_h2
   INTEGER(I4B) :: left,right
   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   IF(x<0.d0)THEN
      WRITE(*,*) 'CubicSplineInterp: r must > 0'
      STOP
   ENDIF
   !Some special case
   IF(x==0.d0)THEN
      CubicHermiteInterp=fun(1)
      RETURN
   ENDIF
   !For r > rmax , fun=0.d0
   IF(x>xmax)THEN
      CubicHermiteInterp=0.d0
      RETURN
   ENDIF
   !interpolation:
   !1.Find the left and right point , dt
   pos=x / h + 1.d0
   !left
   left=FLOOR(pos)
   !right
   right=left+1
   !dxi
   dx0 = x-( left - 1 )*h  !x-x0
   dx1=  dx0-h              !x-x1
   !dxi/h
   dx0_h=dx0/h
   dx1_h=dx1/h
   !(dxi/h)**2
   dx0_h2=dx0_h**2
   dx1_h2=dx1_h**2
   !2.interpolate
   y0=fun(left)
   y1=fun(right)
   dy0=dfdx(left)
   dy1=dfdx(right)
   !
   CubicHermiteInterp= ( y0*(1.d0+2.d0*dx0_h)+dy0*dx0 )*dx1_h2 &
                 & +   ( y1*(1.d0-2.d0*dx1_h)+dy1*dx1 )*dx0_h2
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDFUNCTION CubicHermiteInterp
!--------------------simple interploate------------------------
FUNCTION SimpleInterp(fun,xmax,h,x)
   IMPLICIT NONE
   REAL(DP),INTENT(IN) :: fun(:)   &  !in fun
                      &,  xmax     &  !max x
                      &,  h        &  !grid size
                      &,  x           !in x
   REAL(DP) :: SimpleInterp
   !LOCAL
   INTEGER(I4B) :: left,right
   REAL(DP) :: y0,y1,dx,pos
   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   IF(x<0.d0)THEN
      WRITE(*,*) 'CubicSplineInterp: r must > 0'
      STOP
   ENDIF
   !Some special case
   IF(x==0.d0)THEN
      SimpleInterp=fun(1)
      RETURN
   ENDIF
   !For r > rmax , fun=0.d0
   IF(x>xmax)THEN
      SimpleInterp=0.d0
      RETURN
   ENDIF
   !interpolation:
   !1.Find the left and right point , dt
   pos=x / h + 1.d0
   !left
   left=FLOOR(pos)
   !right
   right=left+1
   !dx
   dx= x - (left-1)*h 
   !
   y0=fun(left)
   y1=fun(right)
   !2.interpolate
   SimpleInterp=y0+dx*(y1-y0)/h
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDFUNCTION SimpleInterp
!-----------------------dYlm-------------------------------
SUBROUTINE r_dYlm(l,m,x,y,z,rmod,f)
   IMPLICIT NONE
   !IN/OUT
   INTEGER(I4B),INTENT(IN) :: l
   INTEGER(I4B),INTENT(IN) :: m
   REAL(DP),INTENT(IN) :: x,y,z,rmod
   REAL(DP),INTENT(INOUT) :: f(0:3)
   !LOCAL
   REAL(DP) :: scal
   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   IF( (rmod<=0.d0) .AND. (l/=0) )THEN
       f(:)=0.d0
       RETURN
   ENDIF
   !
   scal= 1.0_dp /SQRT(4.0_dp*pi) !*(-IMAG)**l
   !Apply Y_lm* (spherical harmonic)
   select case(l)
   case(0)
         select case(m)
         case(0)
                  f(0) = scal*(1)
                  f(1:3)=0
         case default
                  WRITE(*,*) 'Apply Ylm:abs(m)>l'
                  STOP
         end select
   case(1)        
         select case(m)
         case(-1)
                  !f = f*scal*(sqrt(3.0_dp)*x)
                  f(0) = scal*(sqrt(3.0_dp)*x)
                  f(1) = scal*(sqrt(3.0_dp)*(1-x**2)/rmod)
                  f(2) = scal*(sqrt(3.0_dp)*(-x*y/rmod))
                  f(3) = scal*(sqrt(3.0_dp)*(-x*z/rmod))
         case(0)
                  !f = f*scal*(sqrt(3.0_dp)*z)
                  f(0) = scal*(sqrt(3.0_dp)*z)
                  f(1) = scal*(sqrt(3.0_dp)*(-z*x/rmod))
                  f(2) = scal*(sqrt(3.0_dp)*(-z*y/rmod))
                  f(3) = scal*(sqrt(3.0_dp)*(1-z**2)/rmod)
         case(1)
                  !f = f*scal*(-(sqrt(3.0_dp)*y))
                  f(0) = scal*(sqrt(3.0_dp)*y)
                  f(1) = scal*(sqrt(3.0_dp)*(-y*x/rmod))
                  f(2) = scal*(sqrt(3.0_dp)*(1-y**2)/rmod)
                  f(3) = scal*(sqrt(3.0_dp)*(-y*z/rmod))
         case default
                  WRITE(*,*) 'Apply Ylm:abs(m)>l'
                  STOP
         end select
   case(2)
         select case(m)
         case(-2)
                  !f = f*scal*(-(sqrt(15.0_dp)*x*y))
                  f(0) = scal*(-(sqrt(15.0_dp)*x*y))
                  f(1) = scal*(-(sqrt(15.0_dp)*y*(1-2*x**2)/rmod))
                  f(2) = scal*(-(sqrt(15.0_dp)*x*(1-2*y**2)/rmod))
                  f(3) = scal*(-(sqrt(15.0_dp)*(-2*x*y*z)/rmod))
         case(-1)
                  !f = f*scal*(sqrt(15.0_dp)*x*z)
                  f(0) = scal*((sqrt(15.0_dp)*x*z))
                  f(1) = scal*((sqrt(15.0_dp)*z*(1-2*x**2)/rmod))
                  f(2) = scal*((sqrt(15.0_dp)*(-2*x*y*z)/rmod))
                  f(3) = scal*((sqrt(15.0_dp)*x*(1-2*z**2)/rmod))
         case(0)
                  !f = f*scal*((sqrt(5.0_dp)*(-1 + 3*z**2))/2)
                  f(0) = scal*((sqrt(5.0_dp)*(-1 + 3*z**2))/2)
                  f(1) = scal*sqrt(5.0_dp)*3*(-x*z**2/rmod)
                  f(2) = scal*(sqrt(5.0_dp)*(3*(-y*z**2/rmod)))
                  f(3) = scal*((sqrt(5.0_dp)*(3*z*(1-z**2)/rmod)))
         case(1)
                  !f = f*scal*(-(sqrt(15.0_dp)*y*z))
                  f(0) = scal*(-(sqrt(15.0_dp)*y*z))
                  f(1) = scal*(-(sqrt(15.0_dp)*(-2*x*y*z)/rmod))
                  f(2) = scal*(-(sqrt(15.0_dp)*z*(1-2*y**2)/rmod))
                  f(3) = scal*(-(sqrt(15.0_dp)*y*(1-2*z**2)/rmod))
         case(2)
                  !f = f*scal*((sqrt(15.0_dp)*(-x**2 + y**2))/2)
                  f(0) = scal*(sqrt(15.0_dp)*(-x**2 + y**2)/2)
                  f(1) = scal*(sqrt(15.0_dp)*(-x*(1-x**2+y**2))/rmod)
                  f(2) = scal*(sqrt(15.0_dp)*y*(1+x**2-y**2)/rmod)
                  f(3) = scal*(sqrt(15.0_dp)*z*(y**2-x**2)/rmod)
         case default
                  WRITE(*,*) 'Apply Ylm:abs(m)>l'
                  STOP
         end select
   case(3)
         select case(m)
         case(-3)
                  !f = f*scal*((sqrt(17.5_dp)*x*(-1 + 4*y**2 + z**2))/2)
                  f(0) = scal*((sqrt(17.5_dp)*x*(-1 + 4*y**2 + z**2))/2)
                  f(1) = scal*(sqrt(17.5_dp)*((y**2+z**2)/rmod*(-1 + 4*y**2 + z**2)+x*(4*(-2*x*y**2)+(-2*x*z**2))/rmod)/2)
                  f(2) = scal*((sqrt(17.5_dp)*((-x*y/rmod)*(-1 + 4*y**2 + z**2)+x*(4*2*y*(x**2+z**2)+(-2*y*z**2))/rmod))/2)
                  f(3) = scal*((sqrt(17.5_dp)*((-x*z/rmod)*(-1 + 4*y**2 + z**2)+x*(4*2*z*(-2*z*y**2)+(2*z*(x**2+y**2)))/rmod))/2)
         case(-2)
                  !f = f*scal*(-(sqrt(105.0_dp)*x*y*z))
                  f(0) = scal*(-(sqrt(105.0_dp)*x*y*z))
                  f(1) = scal*(-(sqrt(105.0_dp)*(y*z*(y**2+z**2)-2*x**2*y*z))/rmod)
                  f(2) = scal*(-(sqrt(105.0_dp)*(x*z*(x**2+z**2)-2*x*y**2*z))/rmod)
                  f(3) = scal*(-(sqrt(105.0_dp)*(x*y*(x**2+y**2)-2*x*y*z**2))/rmod)
         case(-1)
                  !f = f*scal*((sqrt(10.5_dp)*x*(-1 + 5*z**2))/2)
                  f(0) = scal*((sqrt(10.5_dp)*x*(-1 + 5*z**2))/2)
                  f(1) = scal*((sqrt(10.5_dp)*((y**2+z**2)/rmod*(-1 + 5*z**2)+x*5*(-2*x*z**2)/rmod))/2)
                  f(2) = scal*((sqrt(10.5_dp)*((-x*y/rmod)*(-1 + 5*z**2)+x*5*(-2*y*z**2)/rmod))/2)
                  f(3) = scal*((sqrt(10.5_dp)*((-x*z/rmod)*(-1 + 5*z**2)+x*5*(2*z*(x**2+y**2))/rmod))/2)
         case(0)
                  !f = f*scal*((sqrt(7.0_dp)*z*(-3 + 5*z**2))/2)
                  f(0) = scal*((sqrt(7.0_dp)*z*(-3 + 5*z**2))/2)
                  f(1) = scal*((sqrt(7.0_dp)*((-z*x/rmod)*(-3 + 5*z**2)+z*5*(-2*x*z**2)/rmod))/2)
                  f(2) = scal*((sqrt(7.0_dp)*((-z*y/rmod)*(-3 + 5*z**2)+z*5*(-2*y*z**2)/rmod))/2)
                  f(3) = scal*((sqrt(7.0_dp)*((x**2+y**2)/rmod*(-3 + 5*z**2)+z*5*(2*z*(x**2+y**2))/rmod))/2)
         case(1)
                  !f = f*scal*((sqrt(10.5_dp)*y*(1 - 5*z**2))/2)
                  f(0) = scal*((sqrt(10.5_dp)*y*(1 - 5*z**2))/2)
                  f(1) = scal*((sqrt(10.5_dp)*((-y*x/rmod)*(-3 + 5*z**2)+y*5*(-2*x*z**2)/rmod))/2)
                  f(2) = scal*((sqrt(10.5_dp)*((x**2+z**2)/rmod*(-3 + 5*z**2)+y*5*(-2*y*z**2)/rmod))/2)
                  f(3) = scal*((sqrt(10.5_dp)*((-y*z/rmod)*(-3 + 5*z**2)+y*5*(2*z*(x**2+y**2))/rmod))/2)
         case(2)
                  !f = f*scal*((sqrt(105.0_dp)*(-x**2 + y**2)*z)/2)
                  f(0) = scal*((sqrt(105.0_dp)*(-x**2 + y**2)*z)/2)
                  f(1) = scal*((sqrt(105.0_dp)*((-z*x/rmod)*(-x**2 + y**2)+z*(-2*x*(y**2+z**2)+(-2*x*y**2))/rmod))/2)
                  f(2) = scal*((sqrt(105.0_dp)*((-z*y/rmod)*(-x**2 + y**2)+z*(-(-2*y*x**2)+2*y*(x**2+z**2))/rmod))/2)
                  f(3) = scal*((sqrt(105.0_dp)*((x**2+z**2)/rmod*(-x**2 + y**2)+z*(-(-2*z*x**2)+(-2*z*y**2))/rmod))/2)
         case(3)
                  !f = f*scal*(-(sqrt(17.5_dp)*y*(-3*x**2 + y**2))/2)
                  f(0) = scal*(-(sqrt(17.5_dp)*y*(-3*x**2 + y**2))/2)
                  f(1) = scal*(-(sqrt(17.5_dp)*((-y*x/rmod)*(-3*x**2 + y**2)+y*(-3*2*x*(1-x**2)+(-2*x*y**2))/rmod))/2)
                  f(2) = scal*(-(sqrt(17.5_dp)*((1-y**2)/rmod*(-3*x**2 + y**2)+y*(-3*(-2*y*x**2)+(2*y*(1-y**2)))/rmod))/2)
                  f(3) = scal*(-(sqrt(17.5_dp)*((-y*z/rmod)*(-3*x**2 +  y**2)+y*(-3*(-2*z*x**2)+(-2*z*y**2))/rmod))/2)
         case default
                  WRITE(*,*) 'Apply Ylm:abs(m)>l'
                  STOP
         end select
   case default
         WRITE(*,*) 'Apply Ylm:l>3 not programmed'
         STOP
   end select
   return
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDSUBROUTINE r_dYlm 
!----------------------rcsnTable---------------------------
   SUBROUTINE rcsnTable(n1,n2,n3,dn,h,srmax,rcs,cns,Sphindx,nspt)
      IMPLICIT NONE
      !
      INTEGER(I4B),INTENT(IN) :: n1,n2,n3,dn
      !grid size , sphere rmax
      REAL(DP),INTENT(IN) :: h,srmax
      !(r,cos,sin)
      REAL(DP),DIMENSION(3,-dn+1:n1+dn,-dn+1:n2+dn,-dn+1:n3+dn),INTENT(OUT) :: &
        & rcs
      COMPLEX(DCP),DIMENSION(-dn+1:n1+dn,-dn+1:n2+dn,-dn+1:n3+dn),INTENT(OUT) :: &
        & cns
      LOGICAL,DIMENSION(n1,n2,n3),INTENT(OUT)  :: Sphindx
      INTEGER(I4B),INTENT(OUT) :: nspt
      !LOCAL
      INTEGER(I4B) :: ix,iy,iz,nt
      REAL(DP) :: orig(3),rt,cost,sint,cosp,sinp
      LOGICAL :: lb
      !>>>>>>>>>>>>>>>>>>>>>>>>>START F_DR>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !cubic cell?
      IF((n1/=n2).OR.(n2/=n3).OR.(n1/=n3)) THEN
         WRITE(*,*) 'CGPHI:STOP!!! need cubic cell'
         STOP
      ENDIF
      orig(:)=0.5d0*h*(n1-1)
      nspt=0
      nt=n1+dn
      !print*,'nt',nt,size(Sphindx,1)
      DO iz=-dn+1,nt
      DO iy=-dn+1,nt
      DO ix=-dn+1,nt
         CALL car2spe(ORIG,h,ix,iy,iz,rt,cost,sint,cosp,sinp)
         rcs(1,ix,iy,iz)=rt
         rcs(2,ix,iy,iz)=cost
         rcs(3,ix,iy,iz)=sint
         !boundary?
         lb = (  (ix<1) .OR.(iy<1) .OR.(iz<1)   &
            &.OR.(ix>n1).OR.(iy>n2).OR.(iz>n3) )
         IF(lb)THEN
            !rcsn(4,ix,iy,iz)=CMPLX(cosp,-sinp)
            cns(ix,iy,iz)=CMPLX(cosp,-sinp)
            CYCLE
         ELSE
            !rcsn(4,ix,iy,iz)=CMPLX(cosp,sinp)
            cns(ix,iy,iz)=CMPLX(cosp,sinp)

            IF(rt <= srmax)THEN
               IF((ix==n1).OR.(iy==n2).OR.(iz==n3))THEN
                  PRINT*, 'Spherical reigon exceed boundary'
                  PRINT*, 'ix,iy,iz',ix,iy,iz
                  STOP
               ENDIF
               nspt=nspt+1
               Sphindx(ix,iy,iz) = .TRUE.
            ELSE
               Sphindx(ix,iy,iz) = .FALSE.
            ENDIF

         ENDIF

      ENDDO
      ENDDO
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<< END  F_DR<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE rcsnTable
!----------------------rcsnTable---------------------------
   SUBROUTINE rcsnTable_Atoms(n1,n2,n3,dn,na,h,poscar,srmax,atomR,rVec,rcs,cns,Sphindx,nspt)
      IMPLICIT NONE
      !grids, dgrids, #atom
      INTEGER(I4B),INTENT(IN) :: n1,n2,n3,dn,na
      !grid size , sphere rmax, atomic rcut
      REAL(DP),INTENT(IN) :: h,srmax,atomR
      !atomic position (Cartensian)
      REAL(DP),INTENT(IN) :: poscar(3,na)
      !3D-index
      REAL(DP),DIMENSION(4,-dn+1:n1+dn,-dn+1:n2+dn,-dn+1:n3+dn),INTENT(IN) :: &
        & rVec
      !(r,cos,sin)
      REAL(DP),DIMENSION(3,-dn+1:n1+dn,-dn+1:n2+dn,-dn+1:n3+dn),INTENT(OUT) :: &
        & rcs
      COMPLEX(DCP),DIMENSION(-dn+1:n1+dn,-dn+1:n2+dn,-dn+1:n3+dn),INTENT(OUT) :: &
        & cns
      LOGICAL,DIMENSION(n1,n2,n3),INTENT(OUT)  :: Sphindx
      INTEGER(I4B),INTENT(OUT) :: nspt
      !LOCAL
      INTEGER(I4B) :: ix,iy,iz,nt,Ia
      REAL(DP) :: orig(3),rt,cost,sint,cosp,sinp,ra(3),rat
      LOGICAL :: lb
      !>>>>>>>>>>>>>>>>>>>>>>>>>START F_DR>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !cubic cell?
      IF((n1/=n2).OR.(n2/=n3).OR.(n1/=n3)) THEN
         WRITE(*,*) 'CGPHI:STOP!!! need cubic cell'
         STOP
      ENDIF
      orig(:)=0.5d0*h*(n1-1)
      Sphindx(:,:,:) = .FALSE.
      nspt=0
      nt=n1+dn
      !print*,'nt',nt,size(Sphindx,1)
      DO iz=-dn+1,nt
      DO iy=-dn+1,nt
      DO ix=-dn+1,nt
         CALL car2spe(ORIG,h,ix,iy,iz,rt,cost,sint,cosp,sinp)
         rcs(1,ix,iy,iz)=rt
         rcs(2,ix,iy,iz)=cost
         rcs(3,ix,iy,iz)=sint
         !boundary?
         lb = (  (ix<1) .OR.(iy<1) .OR.(iz<1)   &
            &.OR.(ix>n1).OR.(iy>n2).OR.(iz>n3) )
         IF(lb)THEN
            !rcsn(4,ix,iy,iz)=CMPLX(cosp,-sinp)
            cns(ix,iy,iz)=CMPLX(cosp,-sinp)
            CYCLE
         ELSE
            !rcsn(4,ix,iy,iz)=CMPLX(cosp,sinp)
            cns(ix,iy,iz)=CMPLX(cosp,sinp)

            IF(rt <= srmax)THEN
               !is it boundary?
               IF((ix==n1).OR.(iy==n2).OR.(iz==n3))THEN
                  PRINT*, 'Spherical reigon exceed boundary'
                  PRINT*, 'ix,iy,iz',ix,iy,iz
                  STOP
               ENDIF

               IF(AtomR>0._DP)THEN
               !Atomic test
                  DO Ia=1,na
                     !avoid double counting
                     IF(sphindx(ix,iy,iz)) CYCLE

                     ra(1:3)=rVec(1:3,ix,iy,iz)-poscar(1:3,Ia)
                     rat=SQRT( SUM(ra(:)**2) )
                     IF( rat<AtomR ) THEN

                        nspt=nspt+1
                        Sphindx(ix,iy,iz) = .TRUE.

                     ENDIF

                  ENDDO
               ELSE
               !old version
                  nspt=nspt+1
                  Sphindx(ix,iy,iz) = .TRUE.
               ENDIF

            ENDIF

         ENDIF

      ENDDO
      ENDDO
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<< END  F_DR<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE rcsnTable_Atoms
!------------------------PARTING LINE------------------------------
     SUBROUTINE car2spe(ORIG,h,x,y,z,r,cost,sint,cosp,sinp)
       IMPLICIT NONE
       INTEGER,INTENT(IN) :: x,y,z
       REAL(8),INTENT(IN) :: h,ORIG(3)
       REAL(8),INTENT(OUT) :: r,cost,sint,cosp,sinp
       REAL(8) ::rr,rx,ry,rz
       INTEGER :: l
       !>>>>>>>>>>>>>>>>>>>>>>>start CAR2SPE>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       rx=h*(x-1)-ORIG(1)
       ry=h*(y-1)-ORIG(2)
       rz=h*(z-1)-ORIG(3)
       r=SQRT(rx**2+ry**2+rz**2)
       rr=SQRT(rx**2+ry**2)
       IF(r>0.d0)THEN
          cost=rz/r
          sint=rr/r
       ELSE
          !cost=1.d0
          !sint=0.d0
          WRITE(*,*) 'Multi-pole:the origin on grid'
          STOP
       ENDIF

       IF(rr>0.d0)THEN
         cosp=rx/rr
         sinp=ry/rr
       ELSE
         cosp=1.d0
         sinp=0.d0
         ! WRITE(*,*) 'Multi-pole:the origin on grid'
         ! STOP
       ENDIF
       !<<<<<<<<<<<<<<<<<<<<<<<<<end CAR2SPE<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE car2spe
   !------------------------------------------------------------------------
   SUBROUTINE calclm(Lmax,clm)
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN) :: Lmax
      REAL(DP),INTENT(OUT) ::  clm(0:Lmax,0:Lmax)
      !LOCAL
      INTEGER(I4B) :: l,m
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      clm(:,:)=0._DP
      DO l=0,Lmax
         DO m=0,l
            clm(l,m)=c(l,m)
         ENDDO
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE calclm
   !----------------------------------------------------------------------- 
   FUNCTION c(l,m)
     INTEGER(I4B):: l,m,i
     REAL(DP)    :: c
        if(m==0) then
          c=1.0
          return
        else
          c=1.0
          do i=l-m+1,l+m
             c=i*c
          end do
          c=1.0d0/c
        end if
        return
   ENDFUNCTION c
   !----------------------------------------------------------------------- 
   SUBROUTINE cal_plm(Lmax,x,sx,plm)
      INTEGER(I4B),INTENT(IN) :: Lmax
      REAL(DP),INTENT(IN)     :: x,sx
      REAL(DP),INTENT(OUT)    :: plm(0:Lmax,0:Lmax)
      REAL(DP)    :: pt(0:Lmax,0:Lmax)
      INTEGER(I4B) :: l,m
      REAL(DP)     :: fact,pmm,y

      y=ABS(sx)
      pmm=1.d0
      fact=1.d0
      pt(0,0)=1.d0
      do m=1,Lmax
         pmm=-pmm*fact*y
         pt(m,m)=pmm
         fact=fact+2.d0
      enddo

      !if(l.eq.m) then
      !   p=pmm
      !else
      DO l=0,Lmax-1
         pt(l+1,l)=x*(2*l+1)*pt(l,l)
      ENDDO

      DO m=0,Lmax
         do l=m+2,Lmax
            pt(l,m)=(x*(2*l-1)*pt(l-1,m)-(l+m-1)*pt(l-2,m))/(l-m)
         enddo
      ENDDO
      plm(:,:)=pt(:,:)
   ENDSUBROUTINE cal_plm
   !--------------------------------------------------------
   FUNCTION interp(np,f,r,rnorm,Z)
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN) :: np
      REAL(DP),INTENT(IN) :: f(:),r(:),rnorm
      REAL(DP),OPTIONAL :: Z
      REAL(DP) :: interp
      !LOCAL
      REAL(DP) :: rmax,rmin
      INTEGER(I4B) :: I
      REAL(DP) :: c(3)
      LOGICAL  :: lz
      LOGICAL  :: lfind=.FALSE.
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      lz=PRESENT(Z)
      c(:)=0._DP
      rmax=r(np)
      rmin=r(1)
      IF(rnorm>=rmax)THEN
         IF(lz)THEN
            interp=-Z/rnorm !for local part potential
            RETURN
         ELSE
            interp=0._DP
            RETURN
         ENDIF
      ELSEIF(rnorm<=rmin)THEN
         interp=f(1)
         RETURN
      ELSE
         lfind=.FALSE.
         DO I=1,np-1
            IF((rnorm>=r(I)).AND.(rnorm<r(I+1)))THEN
               lfind=.TRUE.
               EXIT
            ENDIF
         ENDDO
         !Check
         IF(.NOT.lfind)THEN
            WRITE(*,*) 'STOP: Check the interp()',rnorm,rmax
            STOP
         ENDIF
         !interpole it 
         IF(I>=(np-2))THEN
            interp=polynom(0,3,r(np-2:np),f(np-2:np),c,rnorm)  
            RETURN
         ELSEIF(I<=2)THEN
            interp=polynom(0,3,r(1:3),f(1:3),c,rnorm)  
            RETURN
         ELSE
            interp=polynom(0,3,r(I-1:I+1),f(I-1:I+1),c,rnorm)  
            RETURN
         ENDIF
      ENDIF
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDFUNCTION interp
   !--------------------------------------------------------
   FUNCTION interp_dnf(n,np,f,r,rnorm,Z)
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN) :: n,np
      REAL(DP),INTENT(IN) :: f(:),r(:),rnorm
      REAL(DP),OPTIONAL :: Z
      REAL(DP) :: interp_dnf
      !LOCAL
      REAL(DP) :: rmax,rmin
      INTEGER(I4B) :: I
      REAL(DP) :: c(5)
      LOGICAL  :: lz
      LOGICAL  :: lfind=.FALSE.
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      lz=PRESENT(Z)
      IF(n<0) STOP 'Check n in interp_dnf'
      c(:)=0._DP
      rmax=r(np)
      rmin=r(1)
      IF(rnorm>=rmax)THEN
         IF(lz)THEN
            interp_dnf=-Z/rnorm !for local part potential
            IF(n>0)THEN
              DO I=1,n
                 interp_dnf=(-I)*interp_dnf/rnorm
              ENDDO
            ENDIF
            RETURN
         ELSE
            interp_dnf=0._DP
            RETURN
         ENDIF
      ELSEIF(rnorm<rmin)THEN
         interp_dnf=f(1)
         IF(n>0) interp_dnf=0._DP
         RETURN
      ELSE
         lfind=.FALSE.
         DO I=1,np-1
            IF((rnorm>=r(I)).AND.(rnorm<r(I+1)))THEN
               lfind=.TRUE.
               EXIT
            ENDIF
         ENDDO
         !Check
         IF(.NOT.lfind)THEN
            WRITE(*,*) 'STOP: Check the interp_dnf()',rnorm,rmax
            STOP
         ENDIF
         !interpole it 
         !interpole it 
         IF(I>=(np-2))THEN
            interp_dnf=polynom(n,5,r(np-4:np),f(np-4:np),c,rnorm)
            RETURN
         ELSEIF(I<=3)THEN
            interp_dnf=polynom(n,5,r(1:5),f(1:5),c,rnorm)  
            RETURN
         ELSE
            interp_dnf=polynom(n,5,r(I-2:I+2),f(I-2:I+2),c,rnorm)
            RETURN
         ENDIF
      ENDIF
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDFUNCTION interp_dnf
   !--------------------------------------------------------
   FUNCTION polynom(m,np,xa,ya,c,x)
   ! !INPUT/OUTPUT PARAMETERS:
   !   m  : order of derivative (in,integer)
   !   np : number of points to fit (in,integer)
   !   xa : abscissa array (in,real(np))
   !   ya : ordinate array (in,real(np))
   !   c  : work array (out,real(np))
   !   x  : evaluation abscissa (in,real)
   ! !DESCRIPTION:
   !   Fits a polynomial of order $n_p-1$ to a set of $n_p$ points. If $m\ge 0$ the
   !   function returns the $m$th derviative of the polynomial at $x$, while for
   !   $m<0$ the integral of the polynomial from the first point in the array to
   !   $x$ is returned.
   !
   ! !REVISION HISTORY:
   !   Created October 2002 (JKD)
   !EOP
   !BOC
      implicit none
      ! argmuments
      REAL(DP)                        :: polynom
      integer(i4b), intent(in) :: m
      integer(i4b), intent(in) :: np
      real(DP), intent(in) :: xa(np),ya(np)
      real(DP), intent(out) :: c(np)
      real(DP), intent(in) :: x
      ! local variables
      integer(i4b) :: i,j,k
      real(DP) :: x0,x1,x2,x3,y0,y1,y2,y3
      real(DP) :: t0,t1,t2,t3,t4,t5,t6
      real(DP) :: c1,c2,c3,sum
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ! fast evaluations for small np
      select case(np)
      case(1)
        select case(m)
        case(:-1)
          polynom=ya(1)*(x-xa(1))
        case(0)
          polynom=ya(1)
        case default
          polynom=0.d0
        end select
        return
      case(2)
        c1=(ya(2)-ya(1))/(xa(2)-xa(1))
        t1=x-xa(1)
        select case(m)
        case(:-1)
          polynom=t1*(ya(1)+0.5d0*c1*t1)
        case(0)
          polynom=c1*t1+ya(1)
        case(1)
          polynom=c1
        case default
          polynom=0.d0
        end select
        return
      case(3)
        x0=xa(1)
        x1=xa(2)-x0
        x2=xa(3)-x0
        y0=ya(1)
        y1=ya(2)-y0
        y2=ya(3)-y0
        t0=1.d0/(x1*x2*(x2-x1))
        t1=x1*y2
        t2=x2*y1
        c1=x2*t2-x1*t1
        c2=t1-t2
        t1=x-x0
        select case(m)
        case(:-1)
          polynom=t1*(y0+t0*t1*(0.5d0*c1+0.3333333333333333333d0*c2*t1))
        case(0)
          polynom=y0+t0*t1*(c1+c2*t1)
        case(1)
          polynom=t0*(2.d0*c2*t1+c1)
        case(2)
          polynom=t0*2.d0*c2
        case default
          polynom=0.d0
        end select
        return
      case(4)
        x0=xa(1)
        x1=xa(2)-x0
        x2=xa(3)-x0
        x3=xa(4)-x0
        y0=ya(1)
        y1=ya(2)-y0
        y2=ya(3)-y0
        y3=ya(4)-y0
        t0=1.d0/(x1*x2*x3*(x1-x2)*(x1-x3)*(x2-x3))
        t1=x1*x2*y3
        t2=x2*x3*y1
        t3=x3*x1*y2
        c3=t1*(x1-x2)+t2*(x2-x3)+t3*(x3-x1)
        t6=x3**2
        t5=x2**2
        t4=x1**2
        c2=t1*(t5-t4)+t2*(t6-t5)+t3*(t4-t6)
        c1=t1*(x2*t4-x1*t5)+t2*(x3*t5-x2*t6)+t3*(x1*t6-x3*t4)
        t1=x-x0
        select case(m)
        case(:-1)
          polynom=t1*(y0+t0*t1*(0.5d0*c1+t1*(0.3333333333333333333d0*c2 &
           +0.25d0*c3*t1)))
        case(0)
          polynom=y0+t0*t1*(c1+t1*(c2+c3*t1))
        case(1)
          polynom=t0*(c1+t1*(2.d0*c2+3.d0*c3*t1))
        case(2)
          polynom=t0*(6.d0*c3*t1+2.d0*c2)
        case(3)
          polynom=t0*6.d0*c3
        case default
          polynom=0.d0
        end select
        return
      end select
      if (np.le.0) then
        write(*,*)
        write(*,'("Error(polynom): np <= 0 : ",I8)') np
        write(*,*)
        stop
      end if
      if (m.ge.np) then
        polynom=0.d0
        return
      end if
      ! find the polynomial coefficients in divided differences form
      c(:)=ya(:)
      do i=2,np
        do j=np,i,-1
          c(j)=(c(j)-c(j-1))/(xa(j)-xa(j+1-i))
        end do
      end do
      ! special case m=0
      if (m.eq.0) then
        sum=c(1)
        t1=1.d0
        do i=2,np
          t1=t1*(x-xa(i-1))
          sum=sum+c(i)*t1
        end do
        polynom=sum
        return
      end if
      x0=xa(1)
      ! convert to standard form
      do j=1,np-1
        do i=1,np-j
          k=np-i
          c(k)=c(k)+(x0-xa(k-j+1))*c(k+1)
        end do
      end do
      if (m.gt.0) then
      ! take the m th derivative
        do j=1,m
          do i=m+1,np
            c(i)=c(i)*dble(i-j)
          end do
        end do
        t1=c(np)
        t2=x-x0
        do i=np-1,m+1,-1
          t1=t1*t2+c(i)
        end do
        polynom=t1
      else
      ! find the integral
        t1=c(np)/dble(np)
        t2=x-x0
        do i=np-1,1,-1
          t1=t1*t2+c(i)/dble(i)
        end do
        polynom=t1*t2
      end if
      return
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDFUNCTION
   !--------------------------------------------------------
   subroutine fourier_1d(nr,rr,rab,vr,ll,nql,yp,vql,vt)
 
     !use constants
     implicit none
     !     Input/Output variables:
     !     number of grid points
     integer, intent(in) :: nr
     !     coordinates of grid points and differential, dr
     real(dp), intent(in) :: rr(nr), rab(nr)
     !     input function in real space
     real(dp), intent(in) :: vr(nr)
     !     order of Bessel function
     integer, intent(in) :: ll
     !     number of reciprocal space points
     integer, intent(in) :: nql
     !     coordinates of reciprocal space points
     real(dp), intent(in) :: yp(nql)
     !     output function
     real(dp), intent(out) :: vql(nql)
     !     compensating charge
     real(dp), intent(in) :: vt
     !     Work variables:
     !     counters
     integer :: ii, j
     real(dp), dimension(nr+4) :: &
          &     vtmp,     &           ! input function minus compensating function
          &     y                    ! integrand in 1-D Fourier transform
     !
     !     External function:
     !real(dp), external :: besselj
     y = 0._DP
     vtmp = 0._DP
     vql = 0._DP
     do ii = 2, nr
        vtmp(ii) = rr(ii)*(rr(ii)*vr(ii)+vt)
     enddo
     !open(1111,FILE='vtmp.dat')
     !   write(1111,*) vtmp(:)
     !close(1111)
     !g=0
     CALL integ_new(rab,vtmp,vql(1))
     vql(1)=vql(1)*4*pi
     !g>0
     do j = 2, nql
        do ii = 2, nr
           y(ii) = vtmp(ii)*sphbess(ll,yp(j)*rr(ii))
        enddo
        CALL integ_new(rab,y,vql(j))
        !do ii = 1, nr, 4
        !   vql(j) = vql(j) + 7.0*y(ii) + 32.0*y(ii+1) +  &
        !   &    12.0*y(ii+2) + 32.0*y(ii+3) + 7.0*y(ii+4)
        !enddo
        !vql(j) = (2._DP*yp(j)*yp(j)*vql(j)/45.0 - vt)*2._DP/pi
        vql(j) = (yp(j)*yp(j)*vql(j) - vt)*2._DP/pi
        vql(j) = vql(j)*2._DP*pi**2/yp(j)**2
     enddo
     !vql(1)=-2._DP*vt/pi
     !open(1112,FILE='vql.dat')
     !   write(1112,*) vql(:)
     !close(1112)
   end subroutine fourier_1d
   !--------------------------------------------------------
   SUBROUTINE integ_new(rab,y,f)
     IMPLICIT NONE
     REAL(DP),INTENT(IN) :: rab(:),y(:)
     REAL(DP),INTENT(OUT) :: f
     !LOCAL
     INTEGER(I4B) :: nr,i
     REAL(DP) :: yp(SIZE(y))
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     nr=SIZE(rab)
     f=0._DP
     yp=0._DP
     DO i=2,nr
        yp(i)=y(i)*rab(i)
     ENDDO
     !
     DO i=1,nr,4
        f=f+ 7._DP*yp(i) + 32._DP*yp(i+1) + 12._DP*yp(i+2) &
             & + 32._DP*yp(i+3) + 7._DP*yp(i+4)
     ENDDO
     f=f*2._DP/45._DP
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE integ_new
   !--------------------------------------------------------
   FUNCTION sphbess(l,x)
     ! Arguments
     integer, intent(in) :: l
     real(DP), intent(in) :: x
     REAL(DP) :: sphbess
     ! Local variables
     integer :: j
     real(kind=DP), parameter :: third = 1.0_DP / 3.0_DP
     real(kind=DP), parameter :: ftnth = 1.0_DP / 14.0_DP
     real(kind=DP) :: x2, sb0, sb1, by, bym, byp, ux
 
     if (abs(x) > 0.001_DP) then
        sb0 = sin(x)/x
     else
        x2 = 0.5_DP*x*x
        sb0 = 1.0_DP - third*x2*(1.0_DP - 0.1_DP*x2)
     end if
     if (l == 0) then
        sphbess = sb0
     else
        if (abs(x) > 0.001_DP) then
           sb1 = (sb0 - cos(x)) / x
        else
           sb1 = third*x*(1.0_DP - (0.2_DP*x2)*(1.0_DP - ftnth*x2))
        end if
        if (l == 1) then
           sphbess = sb1
        else if (x == 0.0_DP) then
           sphbess = 0.0_DP
        else
           by = sb1
           bym = sb0
           ux = 1.0_DP / x
           do j=1,l-1
              byp = real(2*J+1,DP)*ux*by - bym
              bym = by
              by = byp
           end do
           sphbess = by
        end if
     end if
   end function sphbess
   !--------------------------------------------------------
   SUBROUTINE pullbackPBC(indi,nn,indf)
      !By Qiang XU
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN) :: indi,nn
      INTEGER(I4B),INTENT(OUT) :: indf 
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      indf=MOD(MOD((indi-1),nn)+nn,nn) + 1
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE pullbackPBC
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
end module math

MODULE MathSplines
!
! 		   GNU LESSER GENERAL PUBLIC LICENSE
!                        Version 3, 29 June 2007
!
!  Copyright (C) 2007 Free Software Foundation, Inc. <http://fsf.org/>
!  Everyone is permitted to copy and distribute verbatim copies
!  of this license document, but changing it is not allowed.
!
!
!   This version of the GNU Lesser General Public License incorporates
! the terms and conditions of version 3 of the GNU General Public
! License, supplemented by the additional permissions listed below.
!
!   0. Additional Definitions.
!
!   As used herein, "this License" refers to version 3 of the GNU Lesser
! General Public License, and the "GNU GPL" refers to version 3 of the GNU
! General Public License.
!
!   "The Library" refers to a covered work governed by this License,
! other than an Application or a Combined Work as defined below.
!
!   An "Application" is any work that makes use of an interface provided
! by the Library, but which is not otherwise based on the Library.
! Defining a subclass of a class defined by the Library is deemed a mode
! of using an interface provided by the Library.
!
!   A "Combined Work" is a work produced by combining or linking an
! Application with the Library.  The particular version of the Library
! with which the Combined Work was made is also called the "Linked
! Version".
!
!   The "Minimal Corresponding Source" for a Combined Work means the
! Corresponding Source for the Combined Work, excluding any source code
! for portions of the Combined Work that, considered in isolation, are
! based on the Application, and not on the Linked Version.
!
!   The "Corresponding Application Code" for a Combined Work means the
! object code and/or source code for the Application, including any data
! and utility programs needed for reproducing the Combined Work from the
! Application, but excluding the System Libraries of the Combined Work.
!
!   1. Exception to Section 3 of the GNU GPL.
!
!   You may convey a covered work under sections 3 and 4 of this License
! without being bound by section 3 of the GNU GPL.
!
!   2. Conveying Modified Versions.
!
!   If you modify a copy of the Library, and, in your modifications, a
! facility refers to a function or data to be supplied by an Application
! that uses the facility (other than as an argument passed when the
! facility is invoked), then you may convey a copy of the modified
! version:
!
!    a) under this License, provided that you make a good faith effort to
!    ensure that, in the event an Application does not supply the
!    function or data, the facility still operates, and performs
!    whatever part of its purpose remains meaningful, or
!
!    b) under the GNU GPL, with none of the additional permissions of
!    this License applicable to that copy.
!
!   3. Object Code Incorporating Material from Library Header Files.
!
!   The object code form of an Application may incorporate material from
! a header file that is part of the Library.  You may convey such object
! code under terms of your choice, provided that, if the incorporated
! material is not limited to numerical parameters, data structure
! layouts and accessors, or small macros, inline functions and templates
! (ten or fewer lines in length), you do both of the following:
!
!    a) Give prominent notice with each copy of the object code that the
!    Library is used in it and that the Library and its use are
!    covered by this License.
!
!    b) Accompany the object code with a copy of the GNU GPL and this license
!    document.
!
!   4. Combined Works.
!
!   You may convey a Combined Work under terms of your choice that,
! taken together, effectively do not restrict modification of the
! portions of the Library contained in the Combined Work and reverse
! engineering for debugging such modifications, if you also do each of
! the following:
!
!    a) Give prominent notice with each copy of the Combined Work that
!    the Library is used in it and that the Library and its use are
!    covered by this License.
!
!    b) Accompany the Combined Work with a copy of the GNU GPL and this license
!    document.
!
!    c) For a Combined Work that displays copyright notices during
!    execution, include the copyright notice for the Library among
!    these notices, as well as a reference directing the user to the
!    copies of the GNU GPL and this license document.
!
!    d) Do one of the following:
!
!        0) Convey the Minimal Corresponding Source under the terms of this
!        License, and the Corresponding Application Code in a form
!        suitable for, and under terms that permit, the user to
!        recombine or relink the Application with a modified version of
!        the Linked Version to produce a modified Combined Work, in the
!        manner specified by section 6 of the GNU GPL for conveying
!        Corresponding Source.
!
!        1) Use a suitable shared library mechanism for linking with the
!        Library.  A suitable mechanism is one that (a) uses at run time
!        a copy of the Library already present on the user's computer
!        system, and (b) will operate properly with a modified version
!        of the Library that is interface-compatible with the Linked
!        Version.
!
!    e) Provide Installation Information, but only if you would otherwise
!    be required to provide such information under section 6 of the
!    GNU GPL, and only to the extent that such information is
!    necessary to install and execute a modified version of the
!    Combined Work produced by recombining or relinking the
!    Application with a modified version of the Linked Version. (If
!    you use option 4d0, the Installation Information must accompany
!    the Minimal Corresponding Source and Corresponding Application
!    Code. If you use option 4d1, you must provide the Installation
!    Information in the manner specified by section 6 of the GNU GPL
!    for conveying Corresponding Source.)
!
!   5. Combined Libraries.
!
!   You may place library facilities that are a work based on the
! Library side by side in a single library together with other library
! facilities that are not Applications and are not covered by this
! License, and convey such a combined library under terms of your
! choice, if you do both of the following:
!
!    a) Accompany the combined library with a copy of the same work based
!    on the Library, uncombined with any other library facilities,
!    conveyed under the terms of this License.
!
!    b) Give prominent notice with the combined library that part of it
!    is a work based on the Library, and explaining where to find the
!    accompanying uncombined form of the same work.
!
!   6. Revised Versions of the GNU Lesser General Public License.
!
!   The Free Software Foundation may publish revised and/or new versions
! of the GNU Lesser General Public License from time to time. Such new
! versions will be similar in spirit to the present version, but may
! differ in detail to address new problems or concerns.
!
!   Each version is given a distinguishing version number. If the
! Library as you received it specifies that a certain numbered version
! of the GNU Lesser General Public License "or any later version"
! applies to it, you have the option of following the terms and
! conditions either of that published version or of any later version
! published by the Free Software Foundation. If the Library as you
! received it does not specify a version number of the GNU Lesser
! General Public License, you may choose any version of the GNU Lesser
! General Public License ever published by the Free Software Foundation.
!
!   If the Library as you received it specifies that a proxy can decide
! whether future versions of the GNU Lesser General Public License shall
! apply, that proxy's public statement of acceptance of any version is
! permanent authorization for you to choose that version for the
! Library.
!
! THIS SOFTWARE WAS OBTAINED FROM http://people.sc.fsu.edu/~jburkardt/f_src/spline/spline.f90
! IT HAS BEEN SHORTENED IN ORDER TO BE BETTER USABLE, THE ABOVE LICENSING APPLIES SOLELY
! TO THESE FUNCTIONS!


  USE CONSTANTS, ONLY : DP, I4B

  CONTAINS

subroutine spline_cubic_set ( n, t, y, ibcbeg, ybcbeg, ibcend, ybcend, ypp )
!*******************************************************************************
!
!! SPLINE_CUBIC_SET computes the second derivatives of a cubic spline.
!
!
!  Discussion:
!
!    For data interpolation, the user must call SPLINE_CUBIC_SET to
!    determine the second derivative data, passing in the data to be
!    interpolated, and the desired boundary conditions.
!
!    The data to be interpolated, plus the SPLINE_CUBIC_SET output,
!    defines the spline.  The user may then call SPLINE_CUBIC_VAL to
!    evaluate the spline at any point.
!
!    The cubic spline is a piecewise cubic polynomial.  The intervals
!    are determined by the "knots" or abscissas of the data to be
!    interpolated.  The cubic spline has continous first and second
!    derivatives over the entire interval of interpolation.
!
!    For any point T in the interval T(IVAL), T(IVAL+1), the form of
!    the spline is
!
!      SPL(T) = A(IVAL)
!             + B(IVAL) * ( T - T(IVAL) )
!             + C(IVAL) * ( T - T(IVAL) )**2
!             + D(IVAL) * ( T - T(IVAL) )**3
!
!    If we assume that we know the values Y(*) and YPP(*), which represent
!    the values and second derivatives of the spline at each knot, then
!    the coefficients can be computed as:
!
!      A(IVAL) = Y(IVAL)
!      B(IVAL) = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
!        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
!      C(IVAL) = YPP(IVAL) / 2
!      D(IVAL) = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
!
!    Since the first derivative of the spline is
!
!      SPL'(T) =     B(IVAL)
!              + 2 * C(IVAL) * ( T - T(IVAL) )
!              + 3 * D(IVAL) * ( T - T(IVAL) )**2,
!
!    the requirement that the first derivative be continuous at interior
!    knot I results in a total of N-2 equations, of the form:
!
!      B(IVAL-1) + 2 C(IVAL-1) * (T(IVAL)-T(IVAL-1))
!      + 3 * D(IVAL-1) * (T(IVAL) - T(IVAL-1))**2 = B(IVAL)
!
!    or, setting H(IVAL) = T(IVAL+1) - T(IVAL)
!
!      ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
!      - ( YPP(IVAL) + 2 * YPP(IVAL-1) ) * H(IVAL-1) / 6
!      + YPP(IVAL-1) * H(IVAL-1)
!      + ( YPP(IVAL) - YPP(IVAL-1) ) * H(IVAL-1) / 2
!      =
!      ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
!      - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * H(IVAL) / 6
!
!    or
!
!      YPP(IVAL-1) * H(IVAL-1) + 2 * YPP(IVAL) * ( H(IVAL-1) + H(IVAL) )
!      + YPP(IVAL) * H(IVAL)
!      =
!      6 * ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
!      - 6 * ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
!
!    Boundary conditions must be applied at the first and last knots.
!    The resulting tridiagonal system can be solved for the YPP values.
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, INTEGER(I4B) N, the number of data points; N must be at least 2.
!
!    Input, real T(N), the points where data is specified.
!    The values should be distinct, and increasing.
!
!    Input, real Y(N), the data values to be interpolated.
!
!    Input, INTEGER(I4B) IBCBEG, the left boundary condition flag:
!
!      0: the spline should be a quadratic over the first interval;
!      1: the first derivative at the left endpoint should be YBCBEG;
!      2: the second derivative at the left endpoint should be YBCBEG.
!
!    Input, real YBCBEG, the left boundary value, if needed.
!
!    Input, INTEGER(I4B) IBCEND, the right boundary condition flag:
!
!      0: the spline should be a quadratic over the last interval;
!      1: the first derivative at the right endpoint should be YBCEND;
!      2: the second derivative at the right endpoint should be YBCEND.
!
!    Input, real YBCEND, the right boundary value, if needed.
!
!    Output, real YPP(N), the second derivatives of the cubic spline.
!
  implicit none
!
  INTEGER(I4B) :: n
!
  real(kind=dp) ::  diag(n)
  INTEGER(I4B) :: i
  INTEGER(I4B) :: ibcbeg
  INTEGER(I4B) :: ibcend
  real(kind=dp) ::  sub(2:n)
  real(kind=dp) ::  sup(1:n-1)
  real(kind=dp) ::  t(n)
  real(kind=dp) ::  y(n)
  real(kind=dp) ::  ybcbeg
  real(kind=dp) ::  ybcend
  real(kind=dp) ::  ypp(n)
!
!  Check.
!
  if ( n <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
    write ( *, '(a)' ) '  The number of knots must be at least 2.'
    write ( *, '(a,i6)' ) '  The input value of N = ', n
    stop
  end if

  do i = 1, n-1
    if ( t(i) >= t(i+1) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
      write ( *, '(a)' ) '  The knots must be strictly increasing, but'
      write ( *, '(a,i6,a,g14.6)' ) '  T(',  i,') = ', t(i)
      write ( *, '(a,i6,a,g14.6)' ) '  T(',i+1,') = ', t(i+1)
      stop
    end if
  end do
!
!  Set the first equation.
!
  if ( ibcbeg == 0 ) then
    ypp(1) = 0.0E+00
    diag(1) = 1.0E+00
    sup(1) = -1.0E+00
  else if ( ibcbeg == 1 ) then
    ypp(1) = ( y(2) - y(1) ) / ( t(2) - t(1) ) - ybcbeg
    diag(1) = ( t(2) - t(1) ) / 3.0E+00
    sup(1) = ( t(2) - t(1) ) / 6.0E+00
  else if ( ibcbeg == 2 ) then
    ypp(1) = ybcbeg
    diag(1) = 1.0E+00
    sup(1) = 0.0E+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
    write ( *, '(a)' ) '  The boundary flag IBCBEG must be 0, 1 or 2.'
    write ( *, '(a,i6)' ) '  The input value is IBCBEG = ', ibcbeg
    stop
  end if
!
!  Set the intermediate equations.
!
  do i = 2, n-1
    ypp(i) = ( y(i+1) - y(i) ) / ( t(i+1) - t(i) ) &
           - ( y(i) - y(i-1) ) / ( t(i) - t(i-1) )
    sub(i) = ( t(i) - t(i-1) ) / 6.0E+00
    diag(i) = ( t(i+1) - t(i-1) ) / 3.0E+00
    sup(i) = ( t(i+1) - t(i) ) / 6.0E+00
  end do
!
!  Set the last equation.
!
  if ( ibcend == 0 ) then
    ypp(n) = 0.0E+00
    sub(n) = -1.0E+00
    diag(n) = 1.0E+00
  else if ( ibcend == 1 ) then
    ypp(n) = ybcend - ( y(n) - y(n-1) ) / ( t(n) - t(n-1) )
    sub(n) = ( t(n) - t(n-1) ) / 6.0E+00
    diag(n) = ( t(n) - t(n-1) ) / 3.0E+00
  else if ( ibcend == 2 ) then
    ypp(n) = ybcend
    sub(n) = 0.0E+00
    diag(n) = 1.0E+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
    write ( *, '(a)' ) '  The boundary flag IBCEND must be 0, 1 or 2.'
    write ( *, '(a,i6)' ) '  The input value is IBCEND = ', ibcend
    stop
  end if
!
!  Special case:
!    N = 2, IBCBEG = IBCEND = 0.
!
  if ( ((n == 2) .and. (ibcbeg == 0)) .and. (ibcend == 0) ) then

    ypp(1) = 0.0E+00
    ypp(2) = 0.0E+00
!
!  Solve the linear system.
!
  else

    call s3_fs ( sub, diag, sup, n, ypp, ypp )

  end if

  return
end  SUBROUTINE

subroutine spline_cubic_val ( n, t, y, ypp, tval, yval, ypval, yppval )
!
!*******************************************************************************
!
!! SPLINE_CUBIC_VAL evaluates a cubic spline at a specific point.
!
!
!  Discussion:
!
!    SPLINE_CUBIC_SET must have already been called to define the
!    values of YPP.
!
!    For any point T in the interval T(IVAL), T(IVAL+1), the form of
!    the spline is
!
!      SPL(T) = A
!             + B * ( T - T(IVAL) )
!             + C * ( T - T(IVAL) )**2
!             + D * ( T - T(IVAL) )**3
!
!    Here:
!      A = Y(IVAL)
!      B = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
!        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
!      C = YPP(IVAL) / 2
!      D = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
!
!  Modified:
!
!    20 November 2000
!    30 January 2014 (intents and optimizations)
!
!  Author:
!
!    John Burkardt
!    Johannes M Dieterich
!
!  Parameters:
!
!    Input, INTEGER(I4B) N, the number of data values.
!
!    Input, real T(N), the knot values.
!
!    Input, real Y(N), the data values at the knots.
!
!    Input, real YPP(N), the second derivatives of the spline at the knots.
!
!    Input, real TVAL, a point, typically between T(1) and T(N), at
!    which the spline is to be evalulated.  If TVAL lies outside
!    this range, extrapolation is used.
!
!    Output, real YVAL, YPVAL, YPPVAL, the value of the spline, and
!    its first two derivatives at TVAL.
!
  implicit none
!
  INTEGER(I4B),intent(in):: n
  real(kind=dp),intent(in)::  t(n)
  real(kind=dp),intent(in)::  tval
  real(kind=dp),intent(in)::  y(n)
  real(kind=dp),intent(in)::  ypp(n)
  real(kind=dp),intent(out)::  yppval
  real(kind=dp),intent(out)::  ypval
  real(kind=dp),intent(out)::  yval
  !
  real(kind=dp) :: dt
  real(kind=dp) :: h
  real(kind=dp) :: tmp1, tmp2, tmp3, tmp4
  INTEGER(I4B) :: left
  INTEGER(I4B) :: right
!
!  Determine the interval [T(LEFT), T(RIGHT)] that contains TVAL.
!  Values below T(1) or above T(N) use extrapolation.
!
  call rvec_bracket ( n, t, tval, left, right )
!
!  Evaluate the polynomial.
!
  dt = tval - t(left)
  h = t(right) - t(left)

  tmp1 = ( y(right) - y(left) ) / h
  tmp2 = - ( ypp(right) / 6.0E+00 + ypp(left) / 3.0E+00 ) * h
  tmp3 = ( ypp(right) - ypp(left) )
  tmp4 = tmp3/h

  yval = y(left) &
       + dt * ( ( tmp1 ) &
       + tmp2 &
       + dt * ( 0.5E+00 * ypp(left) &
       + dt * ( tmp3 / ( 6.0E+00 * h ) ) ) )

  ypval = ( tmp1 ) &
       + tmp2 &
       + dt * ( ypp(left) &
       + dt * ( 0.5E+00 * tmp4 ) )

  yppval = ypp(left) + dt * tmp4


  return
end SUBROUTINE

  SUBROUTINE rvec_bracket ( n, x, xval, left, right )
  !  this subroutine is used for spline by spline_cubic_set() and
  !  spline_cubic_val() above.
  !*******************************************************************************
  !
  !! RVEC_BRACKET searches a sorted array for successive brackets of a value.
  !
  !
  !  Discussion:
  !
  !    If the values in the vector are thought of as defining intervals
  !    on the real line, then this routine searches for the interval
  !    nearest to or containing the given value.
  !
  !  Modified:
  !
  !    06 April 1999
  !    30 January 2014 (intents and opts)
  !
  !  Author:
  !
  !    John Burkardt
  !    Johannes M Dieterich
  !
  !  Parameters:
  !
  !    Input, INTEGER(I4B) N, length of input array.
  !
  !    Input, real X(N), an array sorted into ascending order.
  !
  !    Input, real XVAL, a value to be bracketed.
  !
  !    Output, INTEGER(I4B) LEFT, RIGHT, the results of the search.
  !    Either:
  !      XVAL < X(1), when LEFT = 1, RIGHT = 2;
  !      XVAL > X(N), when LEFT = N-1, RIGHT = N;
  !    or
  !      X(LEFT) <= XVAL <= X(RIGHT).
  !
    implicit none
    !
    INTEGER(I4B),intent(in)::n
    INTEGER(I4B),intent(out)::left
    INTEGER(I4B),intent(out)::right
    real(kind=dp),intent(in)::x(n)
    real(kind=dp),intent(in)::xval
    !
    INTEGER(I4B) i1, i2, i3!, i
!goto 20
    !
    ! a faster version, has been tested for Al fcc bulk,
    ! JMD: no idea who was the original author here?!
    !      but it certainly wasn't as fast as it could
    !      have been and now is
    !
    i1 = 1
    i3 = n
    i2 = (i3+i1)/2 !floor(real(i3+i1,kind=dp)*0.5_DP)

  do while (i2/=i1)
    if (xval == x(i2)) then
      left = i2
      right = i2+1
      return
    else if (xval < x(i2)) then
      i3 = i2
    else if (xval > x(i2)) then
      i1 = i2
    end if
    i2 = (i3+i1)/2 !floor(real(i3+i1,kind=dp)*0.5_DP)
  end do

  if (xval == x(i3) ) then
    left = i3
    right = i3 + 1
  else if (xval==x(i1)) then
    left = i1
    right = i1 + 1
  else
    left = i1
    right = i3
  endif
  return

!
! old version
!
!20 continue
!    do i = 2, n - 1
!       if ( xval < x(i) ) then
!           left = i - 1
!           right = i
!           return
!       end if
!    end do
!
!    left = n - 1
!    right = n
!
!    return
  END SUBROUTINE rvec_bracket



  subroutine s3_fs ( a1, a2, a3, n, b, x )
  !
  !*******************************************************************************
  !
  !! S3_FS factors and solves a tridiagonal linear system.
  !
  !
  !  Note:
  !
  !    This algorithm requires that each diagonal entry be nonzero.
  !
  !  Modified:
  !
  !    05 December 1998
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, real A1(2:N), A2(1:N), A3(1:N-1).
  !    On input, the nonzero diagonals of the linear system.
  !    On output, the data in these vectors has been overwritten
  !    by factorization information.
  !
  !    Input, INTEGER(I4B) N, the order of the linear system.
  !
  !    Input/output, real B(N).
  !    On input, B contains the right hand side of the linear system.
  !    On output, B has been overwritten by factorization information.
  !
  !    Output, real X(N), the solution of the linear system.
  !
    implicit none
    !
      INTEGER(I4B) :: n
      !
        real(kind=dp) :: a1(2:n)
        real(kind=dp) :: a2(1:n)
        real(kind=dp) :: a3(1:n-1)
        real(kind=dp) :: b(n)
        INTEGER(I4B) :: i
        real(kind=dp) :: x(n)
        real(kind=dp) :: xmult
        !
                    !  The diagonal entries can't be zero.
                    !
         do i = 1, n
             if ( a2(i) == 0.0E+00 ) then
                 write ( *, '(a)' ) ' '
                 write ( *, '(a)' ) 'S3_FS - Fatal error!'
                 write ( *, '(a,i6,a)' ) '  A2(', i, ') = 0.'
                 return
             end if
         end do

         do i = 2, n-1
             xmult = a1(i) / a2(i-1)
             a2(i) = a2(i) - xmult * a3(i-1)

             b(i) = b(i) - xmult * b(i-1)

         end do


         xmult = a1(n) / a2(n-1)
         a2(n) = a2(n) - xmult * a3(n-1)

         x(n) = ( b(n) - xmult * b(n-1) ) / a2(n)
         do i = n-1, 1, -1
             x(i) = ( b(i) - a3(i) * x(i+1) ) / a2(i)
         end do

         return
END SUBROUTINE
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: polynom
! !INTERFACE:
function polynom(m,np,xa,ya,c,x)
! !INPUT/OUTPUT PARAMETERS:
!   m  : order of derivative (in,integer)
!   np : number of points to fit (in,integer)
!   xa : abscissa array (in,real(np))
!   ya : ordinate array (in,real(np))
!   c  : work array (out,real(np))
!   x  : evaluation abscissa (in,real)
! !DESCRIPTION:
!   Fits a polynomial of order $n_p-1$ to a set of $n_p$ points. If $m\ge 0$ the
!   function returns the $m$th derviative of the polynomial at $x$, while for
!   $m<0$ the integral of the polynomial from the first point in the array to
!   $x$ is returned.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
implicit none
! argmuments
REAL(DP)                        :: polynom
integer(i4b), intent(in) :: m
integer(i4b), intent(in) :: np
real(DP), intent(in) :: xa(np),ya(np)
real(DP), intent(out) :: c(np)
real(DP), intent(in) :: x
! local variables
integer(i4b) :: i,j,k
real(DP) :: x0,x1,x2,x3,y0,y1,y2,y3
real(DP) :: t0,t1,t2,t3,t4,t5,t6
real(DP) :: c1,c2,c3,sum
! fast evaluations for small np
select case(np)
case(1)
  select case(m)
  case(:-1)
    polynom=ya(1)*(x-xa(1))
  case(0)
    polynom=ya(1)
  case default
    polynom=0.d0
  end select
  return
case(2)
  c1=(ya(2)-ya(1))/(xa(2)-xa(1))
  t1=x-xa(1)
  select case(m)
  case(:-1)
    polynom=t1*(ya(1)+0.5d0*c1*t1)
  case(0)
    polynom=c1*t1+ya(1)
  case(1)
    polynom=c1
  case default
    polynom=0.d0
  end select
  return
case(3)
  x0=xa(1)
  x1=xa(2)-x0
  x2=xa(3)-x0
  y0=ya(1)
  y1=ya(2)-y0
  y2=ya(3)-y0
  t0=1.d0/(x1*x2*(x2-x1))
  t1=x1*y2
  t2=x2*y1
  c1=x2*t2-x1*t1
  c2=t1-t2
  t1=x-x0
  select case(m)
  case(:-1)
    polynom=t1*(y0+t0*t1*(0.5d0*c1+0.3333333333333333333d0*c2*t1))
  case(0)
    polynom=y0+t0*t1*(c1+c2*t1)
  case(1)
    polynom=t0*(2.d0*c2*t1+c1)
  case(2)
    polynom=t0*2.d0*c2
  case default
    polynom=0.d0
  end select
  return
case(4)
  x0=xa(1)
  x1=xa(2)-x0
  x2=xa(3)-x0
  x3=xa(4)-x0
  y0=ya(1)
  y1=ya(2)-y0
  y2=ya(3)-y0
  y3=ya(4)-y0
  t0=1.d0/(x1*x2*x3*(x1-x2)*(x1-x3)*(x2-x3))
  t1=x1*x2*y3
  t2=x2*x3*y1
  t3=x3*x1*y2
  c3=t1*(x1-x2)+t2*(x2-x3)+t3*(x3-x1)
  t6=x3**2
  t5=x2**2
  t4=x1**2
  c2=t1*(t5-t4)+t2*(t6-t5)+t3*(t4-t6)
  c1=t1*(x2*t4-x1*t5)+t2*(x3*t5-x2*t6)+t3*(x1*t6-x3*t4)
  t1=x-x0
  select case(m)
  case(:-1)
    polynom=t1*(y0+t0*t1*(0.5d0*c1+t1*(0.3333333333333333333d0*c2 &
     +0.25d0*c3*t1)))
  case(0)
    polynom=y0+t0*t1*(c1+t1*(c2+c3*t1))
  case(1)
    polynom=t0*(c1+t1*(2.d0*c2+3.d0*c3*t1))
  case(2)
    polynom=t0*(6.d0*c3*t1+2.d0*c2)
  case(3)
    polynom=t0*6.d0*c3
  case default
    polynom=0.d0
  end select
  return
end select
if (np.le.0) then
  WRITE(6,*)
  WRITE(6,'("Error(polynom): np <= 0 : ",I8)') np
  WRITE(6,*)
  stop
end if
if (m.ge.np) then
  polynom=0.d0
  return
end if
! find the polynomial coefficients in divided differences form
c(:)=ya(:)
do i=2,np
  do j=np,i,-1
    c(j)=(c(j)-c(j-1))/(xa(j)-xa(j+1-i))
  end do
end do
! special case m=0
if (m.eq.0) then
  sum=c(1)
  t1=1.d0
  do i=2,np
    t1=t1*(x-xa(i-1))
    sum=sum+c(i)*t1
  end do
  polynom=sum
  return
end if
x0=xa(1)
! convert to standard form
do j=1,np-1
  do i=1,np-j
    k=np-i
    c(k)=c(k)+(x0-xa(k-j+1))*c(k+1)
  end do
end do
if (m.gt.0) then
! take the m th derivative
  do j=1,m
    do i=m+1,np
      c(i)=c(i)*dble(i-j)
    end do
  end do
  t1=c(np)
  t2=x-x0
  do i=np-1,m+1,-1
    t1=t1*t2+c(i)
  end do
  polynom=t1
else
! find the integral
  t1=c(np)/dble(np)
  t2=x-x0
  do i=np-1,1,-1
    t1=t1*t2+c(i)/dble(i)
  end do
  polynom=t1*t2
end if
return
end function

end module MathSplines
