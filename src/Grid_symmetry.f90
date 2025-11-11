!Apply Grid symmetry
    MODULE Grid_symmetry
    USE constants
    USE struct_module , ONLY : nsym,num_t,c_i,Otrans,Opsym
    IMPLICIT NONE
    integer(i4b),save  ::  invmap(48)
    CONTAINS
    !==========================================================
    SUBROUTINE symm_density_r(rho)
       USE grid_module , ONLY : nr1,nr2,nr3,ng1,ng2,ng3
       USE Fourier  !, ONLY : FFT
       IMPLICIT NONE
       REAL(DP),INTENT(INOUT) :: rho(nr1,nr2,nr3)
       !LOCAL
       COMPLEX(DP) :: rhog(ng1,ng2,ng3)
       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       IF(nsym<=1) RETURN
       rhog=FFT(nr1,nr2,nr3,rho)
       CALL symmdensityg(rhog)
       rho=FFT(ng1,ng2,ng3,rhog)
       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ENDSUBROUTINE symm_density_r
    
    subroutine symmdensityg(rhog)!{{{
    use Grid_module , only : grid,ng1,ng2,ng3
    implicit none   
    complex(dp), intent(inout)  :: rhog(ng1,ng2,ng3)
    complex(dp),dimension(2*(ng1-1),ng2,ng3)  :: rhogfull, rhogsym
    integer(i4b)  ::  tng1,tng2,tng3,npcell,nrot
    real(dp)  ::  spg_trans(3,48), sym_op(3,3,48) 
    real(dp)  ::  ptrans(3,8),trans_al(3,8)
    integer(i4b)  :: itrans, i 
    integer(i4b),allocatable  ::  ngrid(:)
    !-----------------------------------------------
    tng1 = ng1
    tng2 = ng2
    tng3 = ng3
    !num. of rotation
    nrot = nsym
    !number of prime cell
    npcell = num_t + 1
    !store the operation
    spg_trans(:,:) = otrans(:,:)
    sym_op(:,:,:) = opsym(:,:,:)
    !invers
    call invgrp(sym_op,nsym)
    !full FFT grid
    call fillfullq(tng1,tng2,tng3,rhog,rhogfull)
    !
    trans_al=reshape((/&
      & 0.5d0, 0.5d0, 0.5d0, &
      & 0.0d0, 0.5d0, 0.5d0, &
      & 0.5d0, 0.0d0, 0.5d0, &
      & 0.5d0, 0.5d0, 0.0d0, &
      & 2.d0/3, 1.d0/3, 1.d0/3, &
      & 1.d0/3, 2.d0/3, 2.d0/3, &
      & 1.d0/3, 2.d0/3, 1.d0/3, &
      & 2.d0/3, 1.d0/3, 2.d0/3 /),(/3,8/))

    ptrans(:,1) = 0.d0 
    do itrans = 2,num_t + 1
      ptrans(:,itrans) = trans_al(:,c_i(itrans-1))
    end do 
!if (parallel%isroot) print *,'npcell' , npcell,'c_i',c_i
!if (parallel%isroot) print *,'trnas_al'
!if (parallel%isroot) print *,trans_al(:,c_i(1))
!print *,'ptrans' 
!print *,ptrans
!stop

    if ( num_t == 0) then 
      call rhosymg(2*(tng1-1),tng2,tng3,npcell,nrot,invmap,spg_trans,sym_op,rhogfull,rhogsym)
    else 
      call rhosymg(2*(tng1-1),tng2,tng3,npcell,nrot,invmap,spg_trans,sym_op,rhogfull,rhogsym,ptrans)
    end if 
    !return rhog 
    call fillhalfq(tng1,tng2,tng3,rhogsym,rhog)
    end subroutine symmdensityg!}}}

    subroutine invgrp(s,nrot)!{{{
    implicit none 
    integer(i4b)  :: nrot
    real(dp)  :: s(3,3,48) ,  g(9) 
    integer(i4b)  :: irot , jrot
    real(dp)  :: e(9) = (/1,0,0,0,1,0,0,0,1/)
    
    do irot = 1, nrot
      do jrot = 1, nrot
        !call sgrprd(s(:,:,irot),s(:,:,8),g(:))
        call sgrprd(s(:,:,irot),s(:,:,jrot),g(:))
        if ( sgreql( g(:),e(:) ) ) then 
          invmap(irot) = jrot
          exit
        end if 
      end do 
    end do 

    contains 
      logical function sgreql(mat1,mat2)
      implicit none 
      real(dp)  ::  mat1(9), mat2(9) 
      integer(i4b)  :: i
      sgreql = .false.
      do i = 1,9 
        if ( nint(mat1(i)) /= nint(mat2(i)) ) return 
      end do
      sgreql = .true. 

      return
      end function sgreql
    end subroutine invgrp!}}}

    subroutine sgrprd(s1,s2,s)!{{{
    implicit none
    real(dp)  :: s1(3,3), s2(3,3), s(3,3), help(3,3)
    integer(i4b)  :: sssum
    integer(i4b)  :: i,j,k
    do i = 1,3
      do j = 1,3
        sssum = 0.d0
        do k = 1,3
          sssum = sssum + nint(s1(i,k))*nint(s2(k,j))
        end do 
        help(i,j) = float(sssum)
      end do 
      call sgrcop(help,s)
    end do 
    
    end subroutine !}}}
    
    subroutine sgrcop(s1,s2)!{{{
    implicit none 
    real(dp)  ::  s1(9), s2(9) 
    integer(i4b)  :: i

    do  i = 1,9
      s2(i) = s1(i)
    end do 
    return 
    end subroutine sgrcop!}}}

    subroutine fillfullq(id1,id2,id3,mat1,mat2)!{{{
    implicit none 
    integer(i4b),intent(in)  :: id1,id2,id3
    complex(dp),intent(in)  ::  mat1(:,:,:)
    complex(dp),intent(out) ::  mat2(1:2*(id1-1),id2,id3)
    logical  :: lconjg
    integer(i4b)  :: g1,g2,g3,ig1,ig2,ig3, dim1
    integer(i4b)  :: i,j,k
    integer(i4b)  :: ix,iy,iz
    !--------------------------------------------------------
    mat2 = (0.d0,0.d0)
    dim1 = 2*( id1 - 1 )

    do iz = 1, id3 
      g3= iz - 1
      if ( g3 > id3/2 ) g3 = g3 - id3
      do iy = 1, id2
        g2 = iy - 1
        if ( g2 > id2/2 ) g2 = g2 - id2
        do ix = 1, dim1
          g1 = ix - 1
          if ( g1 > id1 - 1 ) g1 = g1 - dim1
          
          !if ( g1 >= id1 ) cycle

          !if ( ( g1>=0 ) .and. ( g2>=0 ) .and. ( g3>=0)  ) then 
          if ( ( g1>=0 ) ) then 
            !ig1 = g1 + 1
            !ig2 = g2 + 1
            !ig3 = g3 + 1
            !ig1 = mod( g1 + 6*dim1 , 2*id1 ) + 1
            ig1 = mod( g1 + 6*dim1 , dim1 ) + 1
            ig2 = mod( g2 + 6*id2 , id2 ) + 1
            ig3 = mod( g3 + 6*id3 , id3 ) + 1
            lconjg = .false.
          else  
            ig1 = mod( dim1 - g1 , dim1 ) + 1
            ig2 = mod( id2 - g2 , id2 ) + 1
            ig3 = mod( id3 - g3 , id3 ) + 1
            lconjg = .true.
          end if 
           
!if (parallel%isroot) print *,'ix',ix,iy,iz
!if (parallel%isroot) print *,'g1',g1,g2,g3
!if (parallel%isroot) print *,'gx',ig1-1,ig2-1,ig3-1
!if (parallel%isroot) print *,'ig',ig1,ig2,ig3
!if (parallel%isroot) print *,'mat1=',mat1(ig1,ig2,ig3),'wheter conjg',lconjg
!if (parallel%isroot) print *,'---------------------------------------'
          mat2(ix,iy,iz) = mat1(ig1,ig2,ig3)
          if (lconjg) mat2(ix,iy,iz) = conjg( mat2( ix,iy,iz ) )
          
        end do
      end do
    end do 
!open(111,file = 'mat2.dat')
!write(111,*) mat2( id1+1:2*(id1-1) , : , : )
!close(111)
!
!open(111,file = 'mat1.dat')
!write(111,*) mat1(2:id1,:,:)
!close(111)
!
!stop
    
    end subroutine fillfullq!}}}

    subroutine fillhalfq(id1,id2,id3,mat1,mat2)!{{{
    implicit none 
    integer(i4b),intent(in)  :: id1,id2,id3
    complex(dp),intent(in) ::  mat1(1:2*(id1-1),id2,id3)
    complex(dp),intent(out)  ::  mat2(:,:,:)
    logical  :: lconjg
    integer(i4b)  :: g1,g2,g3,ig1,ig2,ig3
    integer(i4b)  :: i,j,k
    integer(i4b)  :: ix,iy,iz
    !--------------------------------------------------------
    do iz = 1, id3
      g3= iz - 1
      if ( g3 > id3/2 ) g3 = g3 - id3
      do iy = 1, id2
        g2 = iy - 1
        if ( g2 > id2/2 ) g2 = g2 - id2
        do ix = 1, id1 
          g1 = ix - 1
          !if ( g3 > id3/2 ) g3 = g3 - id3
          
          !if (g1  <= ( id1 - 1) ) then 
          !  ig1 = g1
          !  ig2 = g2 
          !  ig3 = g3 
          !  lconjg = .false.
          !else 
          !  ig1 = mod( 2*(id1-1)-g1 , 2*(id1-1) ) + 1
          !  ig2 = mod( id2-g2 , id2) + 1 
          !  ig3 = mod( id3-g3 , id3) + 1
          !  lconjg = .true.
          !end if 
           
!if (parallel%isroot) print *,'ig',ig1,ig2,ig3
          mat2(ix,iy,iz) = mat1(ix,iy,iz)
          !mat2(ix,iy,iz) = mat2(ix,iy,iz) * exp( cmplx( 0.d0, -1.d0 * float(g1 + g2 + g3) ) )
          !mat2(ix,iy,iz) = mat2(ix,iy,iz) * cmplx( cos( float(g1 + g2 + g3) ), -1.d0 * sin( float(g1 + g2 + g3) ) ) 
          !if (lconjg) mat2(ix,iy,iz) = conjg( mat2( ix,iy,iz ) )
          
        end do
      end do
    end do 
!stop
    
    end subroutine fillhalfq!}}}

    subroutine rhosymg(id1,id2,id3,npcell,nrot,invmap,spg_trans,sym_op,rho,rhosym,ptrans)!{{{
    implicit none
    integer(i4b),intent(in)  :: id1,id2,id3,npcell,nrot,invmap(48)
    real(dp)  :: spg_trans(3,48) , sym_op(3,3,48)
    complex(dp),intent(in)  :: rho(id1*id2*id3)
    complex(dp),intent(out) :: rhosym(id1*id2*id3)
    real(dp),optional  ::  ptrans(3,8)
    real(dp)  :: done(id1*id2*id3)
    complex(dp)  :: rhosum
    integer(i4b) :: ng,nstar
    real(dp)      :: scal 
    complex(dp)      :: phase(48) , phase1
    real(dp)         :: gphase , cosgphase , singphase
    logical       :: lconjg

    complex(dp)   :: rho1,rho2,rho3
    integer(i4b)  :: g(3) , sg(3) , tindex(4,48)
    integer(i4b)  :: ir,ir1,ir2,ir3,itrans,irot,invirot,igrid,igrid2
    integer(i4b)  :: ig,ig1,ig2,ig3
    !----------------------------------------------------------
    do ir = 1,id1*id2*id3
      done( ir ) = -1.d0
    end do
    
    rhosym = (0.d0,0.d0)
!print *,'id',id1,id2,id3
    do ir3 = 0, id3 - 1
      g(3) = mod(ir3+id3/2-1,id3) - id3/2 + 1
      !g(3) = ir3
      !if ( g(3) > id3/2 ) g(3) = g(3) - id3
      do ir2 = 0,id2 - 1
        g(2) = mod(ir2+id2/2-1,id2) - id2/2 + 1 
        !g(2) = ir2
        !if  ( g(2) > id2/2 ) g(2) = g(2) - id2
        do ir1 = 0, id1-1
!write(*,*) 'invmap'
!write(*,*) invmap
!if (parallel%isroot) write(*,*) '---------------------'
          if ( done( id1*(id2*ir3+ir2) + ir1 + 1 ) > 0.d0 ) cycle
          g(1) = mod( ir1 + id1/2 -1 ,id1) - id1/2 + 1
          !g(1) = ir1 
          !if ( g(1) > id1/2 ) g(1) = g(1) - id1 
!if (parallel%isroot) print *,'ir',ir1,ir2,ir3
!if (parallel%isroot) print *,'g',g(1),g(2),g(3)
          rhosum = (0.d0,0.d0)
          ng = 0
          nstar = 0 
          do irot = 1,nrot
            invirot = invmap(irot) 
            !invirot = irot
            sg(1) = nint( sym_op(1,1,invirot) ) * g(1) + nint( sym_op(2,1,invirot) ) * g(2) + nint( sym_op(3,1,invirot) ) * g(3)
            sg(2) = nint( sym_op(1,2,invirot) ) * g(1) + nint( sym_op(2,2,invirot) ) * g(2) + nint( sym_op(3,2,invirot) ) * g(3)
            sg(3) = nint( sym_op(1,3,invirot) ) * g(1) + nint( sym_op(2,3,invirot) ) * g(2) + nint( sym_op(3,3,invirot) ) * g(3)
            
            !sg(1) = nint( sym_op(1,1,invirot) ) * g(1) + nint( sym_op(1,2,invirot) ) * g(2) + nint( sym_op(1,3,invirot) ) * g(3)
            !sg(2) = nint( sym_op(2,1,invirot) ) * g(1) + nint( sym_op(2,2,invirot) ) * g(2) + nint( sym_op(2,3,invirot) ) * g(3)
            !sg(3) = nint( sym_op(3,1,invirot) ) * g(1) + nint( sym_op(3,2,invirot) ) * g(2) + nint( sym_op(3,3,invirot) ) * g(3)
            
            if ( ( sg(1) < -id1/2 + 1 ) .or. ( sg(1) > id1/2 ) .or. & 
              &  ( sg(2) < -id2/2 + 1 ) .or. ( sg(2) > id2/2 ) .or. & 
              &  ( sg(3) < -id3/2 + 1 ) .or. ( sg(3) > id3/2 ) ) cycle 
           
!print *,sg(1) , sg(2) , sg(3)
            sg(1) = mod( sg(1) + 6*id1, id1)
            sg(2) = mod( sg(2) + 6*id2, id2)
            sg(3) = mod( sg(3) + 6*id3, id3)
            igrid = id1*( id2*sg(3) + sg(2) ) + sg(1) + 1
!print *,sg(1) , sg(2) , sg(3)
            if ( sg(1) <= (id1/2) ) then 
              ig1 = sg(1)
              ig2 = sg(2)
              ig3 = sg(3)
              lconjg = .false.
            else 
              IG1=MOD(id1-SG(1),id1)
              IG2=MOD(id2-SG(2),id2)
              IG3=MOD(id3-SG(3),id3)
              LCONJG=.TRUE.
              !ig1 = mod( sg(1) + 6*id1, id1) 
              !ig2 = mod( sg(2) + 6*id2, id2) 
              !ig3 = mod( sg(3) + 6*id3, id3) 
              !ig1 = id1 - ig1
              !ig2 = id2 - ig2
              !ig3 = id2 - ig3
            end if 

!            IF ( ( SG(3)<= (id3/2) ) .and. ( sg(2) <= (id2/2) ) .and. ( sg(3) <= (id3/2)) ) THEN 
!              IG1=SG(1)
!              IG2=SG(2)
!              IG3=SG(3)
!              LCONJG=.FALSE.
!            ELSE 
!              IG1=MOD(id1-SG(1),id1)
!              IG2=MOD(id2-SG(2),id2)
!              IG3=MOD(id3-SG(3),id3)
!              LCONJG=.TRUE.
!            ENDIF
            IGRID2=ID1*(ID2*IG3+IG2)+IG1+1
!print *,ig1 , ig2 , ig3
!  write(*,*) '-------------------------------'
!if ( (id1*(id2*ir3+ir2) + ir1 + 1==2) .and. (irot == 2) ) stop
            
            !igrid2=id1*(id2*ig3+ig2)+ig1+1
            
            gphase = twopi * ( float( g(1) ) *  spg_trans(1,invirot)   + &
                      &        float( g(2) ) *  spg_trans(2,invirot)   + &
                      &        float( g(3) ) *  spg_trans(3,invirot)  )
            scal = 1.d0 
            cosgphase = cos(gphase) * scal
            singphase = sin(gphase) * scal
            
            do itrans = 2, npcell
              gphase = ( ( spg_trans(1,invirot) + ptrans(1,itrans) ) * float(g(1)) + & 
                     &   ( spg_trans(2,invirot) + ptrans(2,itrans) ) * float(g(2)) + & 
                     &   ( spg_trans(3,invirot) + ptrans(3,itrans) ) * float(g(3)) )
              scal = 1.d0
              cosgphase = cosgphase + (cos(twopi*gphase)*scal)
              singphase = singphase + (sin(twopi*gphase)*scal)
            end do 
            cosgphase = cosgphase / float(npcell )
            singphase = singphase / float(npcell )
            phase1 = cmplx( cosgphase , singphase , kind = dp ) 
            
!if (parallel%isroot) write(*,*) 'phasg=',phase1
!if (parallel%isroot) write(*,*) '------------------------------'
            rho1 = rho(igrid2)
            if (lconjg) rho1 = conjg(rho1)
!if (parallel%isroot) write(*,*) 'symop',invirot
!if (parallel%isroot) write(*,*) nint(sym_op(:,:,invirot))
!if (parallel%isroot) write(*,*) spg_trans(:,irot)
!if (parallel%isroot) write(*,*) 'igrid',igrid2 , lconjg
!if (parallel%isroot) write(*,*) 'rhosum=',rhosum
!if (parallel%isroot) write(*,*) 'rho1=',rho(igrid2)
            rhosum = rhosum + rho1 * conjg(phase1)
!if (parallel%isroot) write(*,*) 'phager=',cosgphase
!if (parallel%isroot) write(*,*) 'phagei=',singphase
!if (parallel%isroot) write(*,*) 'irot=',invirot,'igrid2=',igrid2 
!if (parallel%isroot) write(*,*) 'rhosum=',rhosum
!if (parallel%isroot) write(*,*) '------------------------------'
            
            nstar = nstar + 1 
            
            if ( done(igrid) < 0 ) then 
              ng = ng + 1 
              tindex(1,ng) = sg(1) 
              tindex(2,ng) = sg(2) 
              tindex(3,ng) = sg(3) 
              tindex(4,ng) = 1
              phase(ng) = phase1
            else 
              do ig = 1,ng
                if ( ( sg(1) == tindex(1,ig) ) .and. &
                  &  ( sg(2) == tindex(2,ig) ) .and. &
                  &  ( sg(3) == tindex(3,ig) ) ) goto 110
              
              end do 
        110   continue
              tindex(4,ig) = tindex(4,ig) + 1
              phase(ig) = phase(ig) + phase1
            end if 
    
            done(igrid) = 1.d0
          end do 
!if (parallel%isroot) write(*,*) 'phasig',phase(:)
          
          do ig = 1,ng
            phase(ig) = phase(ig) / float( tindex(4,ig) )
          end do 
    
!if (parallel%isroot) print *,'ir',ir1,ir2,ir3,float(nstar)
!if (parallel%isroot) print *,'rhosum',rhosum
          rhosum = rhosum / float( nstar ) 
          
          do 170 ig = 1,ng
            if ( tindex(1,ig) > ( id1/2 ) )  goto 170
            igrid = id1 * ( id2 * tindex(3,ig) + tindex(2,ig) ) + tindex(1,ig) + 1           
            rhosym(igrid) = rhosum * phase(ig) 
         170 continue
!if (parallel%isroot .and.  ( (id1*(id2*ir3+ir2) + ir1 + 1) == 21270 ) ) then 
!if (parallel%isroot ) then 
!  write(*,*) 'grid ',ir1,ir2,ir3, id1*(id2*ir3+ir2) + ir1 + 1 
!  write(*,*) tindex(1:3,:)
!  write(*,*) 'tindex ng',tindex(4,:)
!  write(*,*) 'nstar = ', nstar
!  write(*,*) 'phase',phase(:)
!  write(*,*) 'rhosum = ', rhosum
!  write(*,*) '-------------------------------'
!end if 
!if (   id1*(id2*ir3+ir2) + ir1 + 1  == 2 ) stop
        end do
      end do
    end do
!if (parallel%isroot) then 
!open(111,file = 'rhogsym.dat')
!do ir3 = 0, id3 - 1 
!do ir2 = 0, id2 - 1
!do ir1 = 0, id1 - 1
!ir = id1*( id2*ir3 + ir2 ) + ir1 + 1
!write(111,*) 'igrid=',ir1,ir2,ir3,ir
!write(111,*) rhosym(ir)
!write(111,*) '----------------------'
!end do 
!end do 
!end do 
!close(111)
!end if 
!call mpi_barrier(parallel%comm,mpinfo)
!
!if (parallel%isroot) then 
!open(111,file = 'rhogsym2.dat')
!write(111,*) rhosym
!close(111)
!end if 
!call mpi_barrier(parallel%comm,mpinfo)
!stop
    
    return
    
    end subroutine rhosymg  !}}}

END MODULE Grid_symmetry

