MODULE smearing_module
  !###########################################################!
  !*For metal integration in BZ                               !
  !*Author : Qiang Xu                                         !
  !*Date   : 2017/08/18                                       !
  !*Methfessel and A.T.Paxton smearing method                 !
  !*Ref    : Phy.Rev B 40,6(1989)                             !
  !###########################################################!
  USE constants
  USE parameters , ONLY : Nsmear,Wsmear
  IMPLICIT NONE
  REAL(DP),ALLOCATABLE :: wke(:,:,:)
  REAL(DP) :: fme & !Fermi level
            &, ets !entropy
CONTAINS
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  SUBROUTINE smear_init(nev)
     !
     USE parameters , ONLY : NSPIN
     USE grid_module , ONLY : nk
     IMPLICIT NONE
     INTEGER(I4B) :: nev !number of state we get
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     CALL destroy_smear()
     IF(.NOT.ALLOCATED(wke))THEN
        ALLOCATE(wke(nev,nk,nspin))
     ENDIF
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE smear_init
  !-------------------------PARTING LINE----------------------
  SUBROUTINE destroy_smear()
     IMPLICIT NONE
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     IF(ALLOCATED(wke))THEN
        DEALLOCATE(wke)
     ENDIF
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE destroy_smear
  !-------------------------PARTING LINE----------------------
  FUNCTION hp(x,n)
  !C  Returns the value of the Hermite polynomial of degree N
  !C  evaluated at X.
  !C  H_0  (x) = 1
  !C  H_1  (x) = 2x
  !C  ...
  !C  H_n+1(x) = 2 x H_n(x) - 2 n H_n-1(x)
  !C  P. Ordejon, June 2003
  !C Passed variables
     IMPLICIT NONE
     !
     INTEGER(I4B)  :: n
     REAL(DP) :: x,hp

     !C Local variables
     INTEGER(I4B)  :: i
     REAL(DP) :: hm1 
     REAL(DP) :: hm2 
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     IF(n .gt. 1000)THEN
        WRITE(6,*)'Fermid: Order of Hermite polynomial too large'
        STOP
     ENDIF

     hp = 1.0d0
     hm2 = 0.0d0
     hm1 = 1.0d0
     DO i = 1,n

        hp = 2.0_dp * (x * hm1 - (i-1) * hm2)
        hm2 = hm1
        hm1 = hp

     ENDDO
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDFUNCTION hp
  !------------------------PARTING LINE-----------------------
   SUBROUTINE smearSN(y,y0,w,sn)
     !
     IMPLICIT NONE
     !
     REAL(DP),INTENT(IN) :: y,y0,w
     REAL(DP),INTENT(OUT) :: sn
     !
     REAL(DP) :: x,e_x2,An
     INTEGER(I4B) :: I,J
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     x=(y-y0)/w
      !
     IF(Nsmear<=0)THEN
        !Fermi-Dirac-distribution
        IF (x.gt.100.D0) THEN
          sn = 0.D0
        ELSEIF (x.lt.-100.d0) THEN
          sn = 1.d0
        ELSE
          sn = 1.d0 / ( 1.d0 + EXP(x) )
        ENDIF
 
     ELSE
        !M-P smearing
        sn =  0.5d0 * erfc(x)

        IF(Nsmear>1)THEN
           e_x2 = EXP(-x*x)
           IF(e_x2>1e-15)THEN
              An = 1.0d0/SQRT(pi)
              DO i = 1,Nsmear

                 !Get coefficients in Hermite-Gauss expansion
                 An = -An / (I * 4.0d0)

                 !Get contribution to step function at order I
                 J = 2*I -1  
                 sn = sn + An * hp(x,j) * e_x2

              ENDDO
           ENDIF
        ENDIF
     ENDIF
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE smearSN
  !-------------------------PARTING LINE----------------------
 SUBROUTINE Fermilevel(ne,nev,nk,wk,eval,sigma)
    !########################################!
    !*   This code reference siesta-4.1      !
    !*Author: Qiang Xu                       !
    !*Date  : 2017/08/17                     ! 
    !########################################!
    !Find Fermi level and k-point weight
    USE parameters , ONLY : NSPIN
    IMPLICIT NONE  
    !INPUT
    REAL(DP),INTENT(IN) :: ne
    REAL(DP),INTENT(IN) :: eval(:,:,:)
    INTEGER(I4B),INTENT(IN) :: nk , &!total k points
              &    nev    !total number of states
    REAL(DP),INTENT(IN) :: wk(:) & !wk0
                        &,sigma !width of smearing
    !OUTPUT
    ! 
    REAL(DP) :: emin,emax
    REAL(DP) :: sumq,E_mu,fi
    INTEGER(I4B) :: Ik,Ispin,Ioc,iter
    INTEGER(I4B) :: NITER=150 !total iter step
    REAL(DP) :: T,drange,totq,sn
    REAL(DP) :: TOL=1e-10
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    wke(:,:,:)=0.d0
    totq=ne
    !
    emin=MINVAL(eval(:,1,:))
    emax=MAXVAL(eval(:,1,:))
    sumq=0.d0
    DO Ispin=1,NSPIN
       DO Ik=1,nk
          DO Ioc=1,nev
             !4.d0 for time reverse symmtery
             wke(Ioc,Ik,Ispin)=wk(Ik)*2.d0/NSPIN
             sumq=sumq+wke(Ioc,Ik,Ispin)
             emin=MIN(emin,eval(Ioc,Ik,Ispin))
             emax=MAX(emax,eval(Ioc,Ik,Ispin))
          ENDDO
       ENDDO
    ENDDO
    !
    fme=emax
    IF(ABS(sumq-totq)<TOL) RETURN
    !
    IF(sumq<totq)THEN
       WRITE(*,*) 'Fermilevel:The states not enough'
       WRITE(*,*) 'totq,sumq=',totq,sumq
       STOP
    ENDIF
    !
    T =MAX(sigma,1.d-6)
    drange = T*SQRT(-LOG(tol*0.01d0))
    emin = emin - drange
    emax = emax + drange
    DO iter = 1,NITER
        fme = 0.5d0*(emin + emax)
        sumq = 0.0d0
        DO Ispin=1,NSPIN
           DO Ik = 1,nk
              DO Ioc=1,nev
                 CALL smearSN(eval(Ioc,Ik,Ispin),fme,T,sn)
                 !4.d0 for time inverse symmtery
                 wke(Ioc,Ik,Ispin) = wk(Ik)*sn*2.d0/NSPIN
                 sumq = sumq + wke(Ioc,Ik,Ispin)
              ENDDO
           ENDDO
        ENDDO
        !find the fermi level ?
        IF(ABS(sumq-totq)<TOL) EXIT
        !next
        IF (sumq<totq) emin = fme
        IF (sumq>totq) emax = fme

    ENDDO
    
    IF (iter>=150) THEN
       WRITE(*,*) 'Fermileval : search fermi level failed'
       WRITE(*,*) 'totq,sumq=',totq,sumq
       STOP
    ENDIF
!print*,'sumq',sumq
!print*,'wke'
!print*,wke
!print*,'eval'
!print*,eval
!STOP
    !entropy
    ets = 0.0d0
    DO Ispin=1,NSPIN
       DO Ik = 1, nk
          DO Ioc = 1,nev

              fi = (NSPIN/2.d0)*wke(Ioc,Ik,Ispin) / wk(Ik)
              E_mu = (eval(Ioc,Ik,Ispin)-fme) / T

              ets = ets + ( 2.d0 * wk(Ik)/NSPIN ) * enpy(E_mu,fi)

          ENDDO
       ENDDO
    ENDDO
    !E=TS
    ets=T*ets
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE Fermilevel
  !-------------------------enpy------------------------------
  FUNCTION enpy(E_mu,fi)
     USE constants
     IMPLICIT NONE
     !IN/OUT
     REAL(DP),INTENT(IN)  :: E_mu , & !(e-\mu)/T
                      &      fi     !occ
     REAL :: enpy !electronic entropy
     !LOCAL
     REAL(DP),PARAMETER :: rtiny=1e-15
     REAL(DP) :: fo,fe
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     IF(Nsmear<=0)THEN
        fo=MAX(fi,rtiny)
        fe=1.d0-fo
        fe=MAX(fe,rtiny)

        enpy = -1.d0 *( fo*log(fo) + fe*log(fe) )

     ELSE
        
        enpy = whg(E_mu,Nsmear)
        
     ENDIF
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDFUNCTION enpy
  !--------------------whg function------------------------
  FUNCTION whg(x,n)
     USE constants
     IMPLICIT NONE
     !Computes the factors to get the entropy term 
     !for the Methfessel-Paxton smearing with Hermite
     !polynomials of order N
     !P. Ordejon, June '03
     !Passed variables
     INTEGER(I4B),INTENT(IN)  :: n
     REAL(DP),INTENT(IN)     :: x
     REAL(DP) :: whg
     !Local variables
     INTEGER(I4B)   :: i
     REAL(DP)       :: a
     REAL(DP)       :: gauss
     REAL(DP)       :: x2
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     x2 = x**2.0d0
     !Get coefficients

      a = 1.0d0/SQRT(pi)
      do i = 1,n
        a = - a / (DBLE(I) * 4.0d0)
      enddo

      gauss = DEXP(-x2)
      whg = 0.0d0

      if (gauss .gt. 1.0d-20) whg = 0.5d0*a*hp(x,2*n)*gauss

      return
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDFUNCTION whg
  !------------------updaterho---------------------
  SUBROUTINE Updaterho_PBC(nps,nev,eig,wk,rhoS,rho)
     USE parameters , ONLY : nspin,Isym,IGamma
     USE grid_module , ONLY : nk,nr &
                    &, eigen_type,sumrhoS
     USE Grid_symmetry, ONLY : Symm_density_r
     USE struct_module, ONLY : dvol
     IMPLICIT NONE
     !INOUT
     INTEGER(I4B),INTENT(IN) :: nps,nev !num. of points and states
     TYPE(eigen_type),INTENT(IN) :: eig
     REAL(DP),INTENT(IN) :: wk(nev,nk,nspin)
     REAL(DP),INTENT(OUT) :: rhoS(nr,nspin),rho(nr)
     !LOCAL
     INTEGER(I4B) :: lft,rit,Is,Ik,Ii,ip,ix,iy,iz,Ii1
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     DO Is=1,Nspin
        rho(:)=0.d0
        !gamma point
        DO Ik=1,nk

           IF(Ik==IGamma)THEN
              !Gamma point
              DO Ii=1,nev

                 IF(wk(Ii,Ik,Is)<xtiny) CYCLE
                 !calculation
                 DO Ip=1,nps

                    rho(Ip)=rho(Ip)+wk(Ii,Ik,Is)*eig%wvfG(Ip,Ii,Is)**2

                 ENDDO
              ENDDO

           ELSE
              !Non-Gamma-k
              DO Ii=1,nev

                 IF(wk(Ii,Ik,Is)<xtiny) CYCLE

                 DO Ip=1,nr

                    rho(Ip)=rho(Ip)+wk(Ii,Ik,Is)* &
                  &   (REAL(eig%wvf(Ip,Ii,Ik,Is))**2+AIMAG(eig%wvf(Ip,Ii,Ik,Is))**2)

                 ENDDO
              ENDDO

           ENDIF

        ENDDO

        !symmetried density
        IF(Isym>=1)THEN
           CALL Symm_density_r(rho)
        ENDIF

        rhoS(:,Is)=rho(:)

     ENDDO
     rhoS(:,:)=rhoS(:,:)/dvol
     !sum rhoS
     CALL sumrhoS(nps,rhoS,rho)
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE Updaterho_PBC
  !-------------------------PARTING LINE----------------------
  SUBROUTINE smear_updaterho(nps,nev,ne,eig,rhoS,rho)
     USE grid_module , ONLY : nk,eigen_type,kpt,sumrhoS
     USE struct_module, ONLY : dvol
     USE parameters, ONLY: nspin,Wsmear
     !
     IMPLICIT NONE
     !INPUT
     INTEGER(I4B),INTENT(IN) :: nps,nev   !number of points/states
     REAL(DP),INTENT(IN)     :: ne !number of charge
     TYPE(eigen_type) :: eig
     !OUT PUT
     REAL(DP),INTENT(OUT) :: rhoS(nps,nspin),rho(nps)
     !LOCAL
     REAL(DP) :: tele
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     CALL Fermilevel(ne,nev,nk,kpt%wk,eig%val,Wsmear)

     !updating the charge density acrroding wke
     CALL updaterho_PBC(nps,nev,eig,wke,rhoS,rho)
     !Checking
     tele=SUM(rho)*dvol
     IF(ABS(tele-ne)>1e-6)THEN
        WRITE(*,*) 'The Density Updating is WROUNG',tele,ne
        !> rescale
        rhoS=rhoS/tele*ne
        call sumrhoS(nps,rhoS,rho)
     ENDIF
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE smear_updaterho
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDMODULE smearing_module
