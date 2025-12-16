MODULE read_module
   !#########################READ MODULE#########################!
   !* FOR    : read the input data                               !
   !* Author : Qiang Xu                                          !
   !* Date   : 2017/07/04                                        !
   !#############################################################!
   USE constants
#ifdef MPI
   USE smpi_math_module
#endif
   IMPLICIT NONE
#ifdef MPI
   INTEGER(I4B) :: ierr
#endif
CONTAINS
   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !------------------------read input.dat------------------------
   SUBROUTINE read_file(infile)
      USE parameters
      USE math, ONLY:change_case,find_keywords
      USE struct_module, ONLY : naty
      IMPLICIT NONE
      !-----------------------------------------------------------
      CHARACTER(8),INTENT(IN) :: infile
      !-----------------------------------------------------------
      INTEGER(I4B)    :: l_str
      INTEGER(I4B)    :: ios,id_ex,id_pound,id_key,id_value !id of key words
      INTEGER(I4B)    :: i,j,k
      INTEGER(I4B)    :: nele !the number of elements
      INTEGER(I4B)    :: vmajor,vminor,vmicro ! XC version
      !
      CHARACTER       :: ch_mark
      CHARACTER(100)  :: instr
      CHARACTER(100)  :: str
      !
      LOGICAL         :: lexist
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !-------start input.dat-----------
      INQUIRE(FILE=infile,EXIST=lexist)
      !test the input.dat
      IF(.NOT.lexist)THEN
#ifdef MPI
         IF (parallel%isroot) THEN
            WRITE(*,*) '>>>WARNING<<<:input file is not exist.stop!!!'
            CALL MPI_ABORT(parallel%comm,mpinfo,ierr)
         ENDIF
#else
         WRITE(*,*) '>>>WARNING<<<:input file is not exist.stop!!!'
         STOP
#endif
      ENDIF
      !
      OPEN(100,FILE=infile,STATUS='old')
      ch_mark='='
      DO WHILE(.TRUE.)
         READ(100,'(A100)',IOSTAT=ios) instr
         IF( ios /= 0 ) EXIT
         instr=ADJUSTL(instr)
         !change the upper case to lower case
         CALL change_case(instr,str,2)
         !the length of strings
         l_str=LEN_TRIM(str)
         !mark '!' '#'
         id_ex=INDEX(str(1:l_str),'!')
         id_pound=INDEX(str(1:l_str),'#')
         !cancel the strings after '!' or '#'
         IF (id_ex > 0 .AND. id_pound > 0 ) THEN
            l_str=MIN(l_str,id_ex-1,id_pound-1)
         ELSEIF(id_ex > 0 .AND. id_pound == 0 ) THEN
            l_str=MIN(l_str,id_ex-1)
         ELSEIF(id_ex == 0 .AND. id_pound > 0 ) THEN
            l_str=MIN(l_str,id_pound-1)
         ENDIF
         l_str=LEN_TRIM(str(1:l_str))
         IF(l_str<2) CYCLE
         !
         CALL find_keywords(str,ch_mark,id_key,id_value)
         !System name,no use
         IF(str(1:id_key) =='system')THEN
            str=instr
            READ(str(id_value:l_str),*) system_name
#ifdef MPI
      IF (parallel%isroot) THEN
#endif
            WRITE(*,*) 'Task name >>> ',system_name
#ifdef MPI
      ENDIF
#endif
         ENDIF
         !k-grids
         IF(str(1:id_key)=='kspacing')THEN
            READ(str(id_value:l_str),*) kspacing
            kspacing=kspacing*bohr2ang !/(2._DP*pi)
         ENDIF
         !kgrid mesh
         IF(str(1:id_key)=='kgrids')THEN
            READ(str(id_value:l_str),*) kgrid
         ENDIF
         !exchange-correlation type
         IF( str(1:id_key) == 'xc' )THEN
            READ(str(id_value:l_str),*) ixc
         ENDIF
         !# of elments
         IF( str(1:id_key) == 'nele' )THEN
            READ(str(id_value:l_str),*) naty
         ENDIF
         !finite difference order
         !finite difference order
         IF( str(1:id_key) == 'norder' )THEN
            READ(str(id_value:l_str),*) nfd
            nfd=nfd/2
            IF (nfd>10 .or. nfd<=0) THEN
               WRITE(*,*) "STOP!!:the order of finite difference should in [2-20]"
               STOP
            ENDIF
#ifdef MPI
      IF (parallel%isroot) THEN
#endif
            WRITE(*,*) "The order of finite difference is:",nfd*2
#ifdef MPI
      ENDIF
#endif
         ENDIF
         !diag(H)
         IF(str(1:id_key)=='nadds')THEN
            READ(str(id_value:l_str),*) nadds
         ENDIF
         !>>>spin
         IF(str(1:id_key)=='lspin')THEN
            IF(str(id_value:id_value)=='t')THEN
                nspin=2
#ifdef MPI
      IF (parallel%isroot) THEN
#endif
                WRITE(*,*) "Spin polarized"
#ifdef MPI
      ENDIF
#endif
            ELSE
                nspin=1
#ifdef MPI
      IF (parallel%isroot) THEN
#endif
                WRITE(*,*) "Spin unpolarized"
#ifdef MPI
      ENDIF
#endif
            ENDIF
         ENDIF
         !---
         !>>>smearing for metal
         IF( str(1:id_key) == 'nsmear' )THEN
            READ(str(id_value:l_str),*) Nsmear
         ENDIF
         !---
         IF( str(1:id_key) == 'wsmear' )THEN
            READ(str(id_value:l_str),*) Wsmear
            !Wsmear to hartree
            Wsmear=Wsmear/hart2ev
         ENDIF
         !CheM
         IF(str(1:id_key) == 'chem')THEN
            READ(str(id_value:l_str),*) CheM
            IF(Idiag/=0)THEN
#ifdef MPI
      IF (parallel%isroot) THEN
#endif
                WRITE(*,*) "Chebyshev filter is used,order:",CheM
#ifdef MPI
      ENDIF
#endif
            ENDIF
         ENDIF
         !CheM0
         IF(str(1:id_key) == 'chem0')THEN
            READ(str(id_value:l_str),*) CheM0
            IF(Idiag/=0)THEN
#ifdef MPI
      IF (parallel%isroot) THEN
#endif
                WRITE(*,*) "first Chebyshev filter order:",CheM0
#ifdef MPI
      ENDIF
#endif
            ENDIF
         ENDIF
         !rtol
         IF(str(1:id_key)=='rtol')THEN
            READ(str(id_value:l_str),*) rtol
         ENDIF
         !etol
         IF(str(1:id_key)=='etol')THEN
            READ(str(id_value:l_str),*) etol
            etol=etol/hart2ev
         ENDIF
         !max step
         IF( str(1:id_key) == 'nssp' )THEN
            READ(str(id_value:l_str),*) nssp
            nssp=MAX(nssp,0)
#ifdef MPI
      IF (parallel%isroot) THEN
#endif
            WRITE(*,*) 'Max simulate steps is:',nssp
#ifdef MPI
      ENDIF
#endif
         ENDIF
         !=====================optimization=====================
         IF( str(1:id_key) == 'iopt' )THEN
            READ(str(id_value:l_str),*) iopt
         ENDIF
         !Chebyshev tolerance
         IF(str(1:id_key)=='chetol')THEN
            READ(str(id_value:l_str),*) Chetol
            Chetol=Chetol/hart2ev
         ENDIF
         !pen
         IF(str(1:id_key)=='penlam')THEN
            READ(str(id_value:l_str),*) penLambda
         ENDIF
         !isym
         IF(str(1:id_key)=='isym')THEN
             !symmetry
             READ(str(id_value:l_str),*) isym
         ENDIF

         !>>>fast
         IF(str(1:id_key)=='lfast')THEN
            IF(str(id_value:id_value)=='t')THEN
                lfast=.TRUE.
            ELSE
                lfast=.FALSE.
            ENDIF
         ENDIF


      ENDDO

      CLOSE(100)
      !read POSCAR
      CALL read_pos()
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<      
   END SUBROUTINE read_file
   !---------------------------read poscar------------------------
   SUBROUTINE read_pos()
      USE math, ONLY: dir2car,car2dir
      USE parameters, ONLY : nskip
      USE struct_module,ONLY: lat_mat,natom,naty,creat_struct,struct
      IMPLICIT NONE
      !IN/OUT
      !
      INTEGER(I4B)  :: i,j,icont
      INTEGER(I4B)  :: filestatu
      INTEGER(I4B)  :: ele_n(naty)
      REAL(DP)      :: lat_ratio
      REAL(DP)      :: lat_read(3,3)
      CHARACTER(1)  :: chr
      CHARACTER(30) :: ctmp
      LOGICAL       :: ldir
      !
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      OPEN(101,FILE='CHGCAR',IOSTAT=filestatu)

         IF(filestatu/=0)THEN
#ifdef MPI
      IF (parallel%isroot) THEN
            WRITE(*,*) 'STOP:Could not open CHGCAR file'
            CALL MPI_ABORT(parallel%comm,mpinfo,ierr)
      ENDIF
#else
            WRITE(*,*) 'STOP:Could not open CHGCAR file'
            STOP
#endif
         ENDIF

          READ(101,*)          !title,no use
          READ(101,*) lat_ratio
          !read lattice matrix
          DO i=1,3
             READ(101,*) lat_read(:,i)
          ENDDO
          !turn to atomic unit
          lat_mat(:,:)=lat_read(:,:)*lat_ratio*ang2bohr
          !read over poscar
          icont=5
          DO
             READ(101,'(A)') ctmp
             icont=icont+1
             j=INDEX(ctmp,'Dir')
             IF(j/=0)THEN
                ldir=.TRUE.
                EXIT
             ENDIF

             j=INDEX(ctmp,'Car')
             IF(j/=0)THEN
                ldir=.false.
                EXIT
             ENDIF

             READ(ctmp,*) chr
             IF(LGE(chr,'0').AND.LLE(chr,'9'))THEN
                BACKSPACE(101)
                READ(101,*) ele_n(1:naty)
             ENDIF
          ENDDO
          !creat struct
          CALL creat_struct(ele_n)
          !total atom
          nskip=icont+natom
          !read atom positions
          IF (ldir) THEN
             DO i=1,natom
                READ(101,*) struct%posdir(:,i)
             ENDDO
             CALL dir2car(struct%posdir,struct%poscar,lat_mat)
          ELSE
             DO i=1,natom
                READ(101,*) struct%poscar(:,i)
             ENDDO
             struct%poscar=struct%poscar*ang2bohr
             CALL car2dir(struct%poscar,struct%posdir,lat_mat)
          ENDIF
#ifdef MPI
      IF (parallel%isroot) THEN
#endif
          PRINT*,'Skip lines in CHGCAR and LOCPOT',nskip
#ifdef MPI
      ENDIF
#endif
      CLOSE(101)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE read_pos
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
END MODULE read_module   
