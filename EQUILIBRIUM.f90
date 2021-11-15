!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!SLAB AT EQUILIBRIUM COUPLED TO TWO EXTARNAL LEADS
!Author G Mazza
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE SLAB_EQU
  use global
  use LANCELOT_simple_double
  
  !##################################################################################
  implicit none
  private
  real(8)                                     :: energy
  real(8),dimension(:),allocatable            :: eZ,eZL,eZR
  complex(8),dimension(:),allocatable         :: eZnn,eZnnL,eZnnR
  complex(8),dimension(:),allocatable         :: hop_plus,hop_minus
  complex(8)                                  :: Delta_1L,Delta_NR

  complex(8),dimension(:),allocatable         :: R
  real(8),dimension(:),allocatable            :: dop
  real(8)                                     :: hop_ene,hubb_ene
  real(8)                                     :: e_layers,e_orth

  type(gz_projector),dimension(:),allocatable :: opt_phi

  real(8) :: beta_equ
  real(8) :: beta_neq

  integer :: Lneq

  public :: GZ_equilibrium
    
CONTAINS

  !+-----------------------------------------------------------------------------+!
  !                         GUTZWILLER MINIMITAZION ALGORITHM                     !
  !     -start from some GZ projectors (not interacting or U eq bulk values)      !
  !     -Diagonalize H_star                                                       !
  !     -Obtain new GZ projectors                                                 !
  !     -Iterate                                                                  !
  !+-----------------------------------------------------------------------------+!

  SUBROUTINE GZ_equilibrium(temp_,Gequ,PHIequ,Lneq_)
    integer :: i,j,ik,iL
    integer,optional :: Lneq_
    real(8) :: Qstep,nki,nki_plus
    real(8),dimension(2) :: temp_
    real(8),dimension(:),allocatable :: Zold,test_norm
    complex(8),dimension(:,:,:),allocatable :: Gequ
    type(gz_projector),dimension(:),allocatable :: PHIequ


    Lneq=L
    if(present(Lneq_)) Lneq=Lneq_

    allocate(R(L),dop(L),Zold(L),test_norm(L))
    allocate(eZ(L),eZL(N),eZR(N))
    allocate(eZnn(L-1),eZnnL(N-1),eZnnR(N-1))
    allocate(hop_plus(L-1),hop_minus(L-1))    
    allocate(opt_phi(L))

    opt_Phi = GZ_init_half(U)
    R = GZ_hop (opt_Phi)
    dop = gz_doping(opt_Phi)

    beta_equ=1.d0/temp_(1)
    beta_neq=1.d0/temp_(2)

    do i=1,n_max


       call slater_step
       call projector_step
       call check_convergence(i)

       energy = hop_ene + hubb_ene/(1.d0*L)
       write(*,'(A,I4,5f18.10)')'step numr',i, energy,Qstep,hubb_ene/(1.d0*L), &
            hop_ene

       if(Qstep.lt.conv_treshold) then
          exit
       end if

    end do

    test_norm = gz_normalization(opt_Phi)
    R = GZ_hop (opt_Phi)
    dop = gz_doping(opt_Phi)



    if(allocated(Gequ)) deallocate(Gequ)
    allocate(Gequ(Nk_tot,L,L))
    call slater_step(0,Gequ)
    eZ=0.d0
    hop_plus=0.d0
    do iL=1,L
       do ik=1,Nk_tot
          nki = 1.d0-Zi*Gequ(ik,iL,iL)
          eZ(iL) = eZ(iL) + nki*wt(ik)*epsik(ik)
          if(iL.lt.L) then
             nki_plus = -Zi*Gequ(ik,iL,iL+1)
             hop_plus(iL) = hop_plus(iL) + 2.d0*nki_plus*wt(ik)*Hslab(iL,iL+1)
          end if
       end do
    end do

    if(allocated(PHIequ)) deallocate(PHIequ)
    allocate(PHIequ(L))
    PHIequ(:)%p0 = opt_phi(:)%p0
    PHIequ(:)%p1 = opt_phi(:)%p1
    PHIequ(:)%p2 = opt_phi(:)%p2


    open(10,file='equ_gz.out',position='append')
    do i=1,L
       if(i.lt.L) then
          write(10,"(I4,20(F18.10))")  &
               i                     &
               ,real(opt_Phi(i)%p0)    &
               ,real(opt_Phi(i)%p2)    &
               ,real(opt_Phi(i)%p1)    &
               , test_norm(i)        &
               ,ABS(R(i))**2        &
               ,eZ(i) &
               ,hop_plus(i)
       else
          write(10,"(I4,20(F18.10))")  &
               i                     &
               ,real(opt_Phi(i)%p0)    &
               ,real(opt_Phi(i)%p2)    &
               ,real(opt_Phi(i)%p1)    &
               , test_norm(i)        &
               ,ABS(R(i))**2        &
               ,eZ(i)   

       end if
    end do
    write(10,*)
    close(10)

    deallocate(R,dop,Zold,test_norm)
    deallocate(eZ,eZL,eZR)
    deallocate(eZnn,eZnnL,eZnnR)
    deallocate(hop_plus,hop_minus)
    deallocate(opt_phi)



  CONTAINS


    SUBROUTINE check_convergence(iter)
      integer, intent(in) :: iter
      integer :: j
      if(iter.eq.1) then
         Zold = ABS(R)**2+0.5
      end if
      Qstep = ABS(SUM(ABS(R)**2-Zold))/(1.d0*L)
      Zold = ABS(R)**2
    END SUBROUTINE check_convergence

  END SUBROUTINE GZ_equilibrium




  !+----------------------------------------+!
  !  Slater Determinant optimization step    !
  !  Given the Gutzeiller projectors build   !
  !  up and diagonalize the real space       !
  !  hopping hamiltonian.                    !
  !  Compute                                 !
  !  - average intra-layer hop. energy eZ    !
  !  - average inter-layer hop. energy eZnn  !
  !  - average hybridization energy          !
  !  - total hopping energy                  !
  !+----------------------------------------+!
  SUBROUTINE slater_step(flag_nz,Gequ)
    integer,intent(in),optional  :: flag_nz
    complex(8),dimension(Nk_tot,L,L),optional :: Gequ
    integer     :: ix,iy,ik,i,j
    integer     :: i_max!highest occupied state index
    type(vec2D) :: k
    real(8)     :: ek,ntot,fk
    complex(8),dimension(:,:),allocatable       :: H
    real(8),dimension(:),allocatable            :: w
    real(8),dimension(:),allocatable            :: nz,nkz
    real(8)     :: R_ave,beta_qp
    character(len=1)        :: jobz
    character(len=1)        :: uplo
    jobz = 'V'
    uplo = 'L'


    allocate(H(L,L),w(L),nz(L),nkz(L))

    nz=0.d0
    ik = 0 !k_index
    call init_to_zero

    if(present(flag_nz)) then
       Gequ = Z0
    end if


    if(.not.present(flag_nz)) then
       do ik=1,Nk_tot
          call slab_k(ik)            
          call accumulate_ksum(ik)  !accumulate k sums of hopping matrix elements
          nz = nz + nkz   
       end do
    else

       if(present(Gequ)) then
          R_ave = sum(R)/dble(L)
          if(abs(R_ave).gt.1.d-4) then
             write(*,*) 'metallic',R_ave
             do ik=1,Nk_tot
                call slab_k(ik)            
                call accumulate_ksum(ik)  
                call get_equilibriumG(ik,Gequ(ik,:,:))
                nz = nz +nkz
             end do
          else
             write(*,*) 'insulating',R_ave
             do ik=1,Nk_tot
                ek = epsik(ik)
                H(:,:) = 0.d0
                do i=1,L
                   H(i,i) =  ek
                   if(i.lt.L) then
                      H(i,i+1) = -t_perp!1.d0
                      H(i+1,i) = -t_perp!1.d0
                   end if
                end do
                if(pbc) then
                   H(1,L) = -t_perp!1.d0
                   H(L,1) = -t_perp!1.d0
                end if
                call matrix_diagonalize(H,w,jobz,uplo)
                call get_equilibriumG(ik,Gequ(ik,:,:))
                nkz(:) = 0.d0
                do i=1,L
                   do j=1,L
                      nkz(j) = nkz(j) + abs(H(j,i))**2*wt(ik)*fermi(w(i),beta_equ)
                   end do
                end do
                nz = nz + nkz
             end do
          end if
       end if
    end if

    if(present(flag_nz)) then
       open(30,file='equ_profile.out')
       ntot =0.d0
       do i=1,L
          write(30,*) i,nz(i)
          ntot = ntot + nz(i)
       end do
       write(30,*)
       write(30,*)
       write(30,*) dble(ntot)/dble(L)
       close(30)
    end if

    call  slab_tot_ene

    deallocate(H,w,nz,nkz)

  CONTAINS

    SUBROUTINE slab_k(ik)
      !type(vec2D), intent(in) :: k
      integer                 :: ik
      character(len=1)        :: jobz
      character(len=1)        :: uplo
      jobz = 'V'
      uplo = 'L'
      call build_hop(ik)
      call  matrix_diagonalize(H,w,jobz,uplo)
    END SUBROUTINE slab_k

    SUBROUTINE ho_state(ef)
      real(8) :: ef
      integer :: i
      i_max=L

      do i=1,L
         if(w(i).gt.ef) then
            i_max = i-1
            !            write(*,*) w(i),w(i-1),ik
            goto 100
         end if
      end do

100   continue
    END SUBROUTINE ho_state



    SUBROUTINE build_hop(ik)
      !type(vec2D), intent(in) :: k
      integer :: ik
      integer :: i,j
      real(8) :: ek
      ek = epsik(ik)
      H(:,:) = 0.d0
      !+-  SLAB  -+!
      do i=1,L

         H(i,i) =  ek*abs(r(i))**2*Hslab(i,i)
         do j=i+1,L
            H(i,j) = conjg(r(i))*r(j)*Hslab(i,j)
            H(j,i) = conjg(r(j))*r(i)*Hslab(j,i)
         end do
      end do
      !
      ! do i=1,L
      !    H(i,i) =  ek*abs(r(i))**2*Hslab(i,i)
      !    if(i.lt.L) then
      !       H(i,i+1) = -1.d0*t_perp*conjg(r(i))*r(i+1)
      !       H(i+1,i) = -1.d0*t_perp*r(i)*conjg(r(i+1))
      !    end if
      ! end do


      if(pbc) then
         H(1,L) = conjg(r(1))*r(L)*Hslab(1,2)
         H(L,1) = conjg(r(L))*r(1)*Hslab(2,1)
      end if




    END SUBROUTINE build_hop

    SUBROUTINE get_equilibriumG(ik,Gk)
      integer,intent(in) :: ik
      complex(8),dimension(L,L) :: Gk
      integer :: i,j,iE
      
      do i=1,L
         do j=1,L
            Gk(i,j) = 0.d0
            do iE = 1,L
               if(i.le.Lneq.and.j.le.Lneq) then
                  Gk(i,j) = Gk(i,j) + Zi*conjg(H(j,iE))*H(i,iE)*fermi(w(iE),beta_neq)
               else
                  Gk(i,j) = Gk(i,j) + Zi*conjg(H(j,iE))*H(i,iE)*fermi(w(iE),beta_equ)
               end if
            end do
            if(i.eq.j) then
               Gk(i,j) = Gk(i,j) - Zi

            end if
         end do
      end do
    END SUBROUTINE get_equilibriumG


    SUBROUTINE accumulate_ksum(ik)
      integer,intent(in) :: ik
      integer :: i,j
      integer :: i_slab
      integer :: iL,iR
      real(8) :: ek,beta_qp
      beta_qp = beta_equ
      nkz(:) = 0.d0
      ek = epsik(ik)
      do i=1,L
         do j=1,L
            ! if(i.le.Lneq.and.j.le.Lneq) then
            !    beta_qp=beta_neq
            ! else
            !    beta_qp=beta_equ
            ! end if
            nkz(j) = nkz(j) + abs(H(j,i))**2*wt(ik)*fermi(w(i),beta_qp)
            eZ(j) = eZ(j) + 2.d0*ek*abs(H(j,i))**2*wt(ik)*fermi(w(i),beta_qp)
            if(j.lt.L) then
               hop_plus(j) = hop_plus(j) + 2.d0*conjg(H(j+1,i))*H(j,i)*wt(ik)*fermi(w(i),beta_qp)
               hop_minus(j) = hop_minus(j) + 2.d0*conjg(H(j,i))*H(j+1,i)*wt(ik)*fermi(w(i),beta_qp)
            end if
         end do
      end do
    END SUBROUTINE accumulate_ksum


    SUBROUTINE slab_tot_ene
      integer  :: i,j

      !slab

      hop_ene = 0.d0

      do i=1,L

         e_layers = e_layers + eZ(i)*abs(r(i))**2

         hop_ene = hop_ene + eZ(i)*abs(r(i))**2

         if (i.lt.L) then
            hop_ene = hop_ene - conjg(r(i+1))*r(i)*hop_plus(i)
         end if

         if(i.gt.1) then
            hop_ene = hop_ene - conjg(r(i-1))*r(i)*hop_minus(i-1)
         end if


      end do

      hop_ene = hop_ene/(1.d0*L)

    END SUBROUTINE slab_tot_ene


    SUBROUTINE init_to_zero
      Delta_1L  = 0.d0
      Delta_NR  = 0.d0
      eZ(:)     = 0.d0
      eZnn(:)   = 0.d0
      eZL(:)    = 0.d0
      eZnnL(:)  = 0.d0
      eZR(:)    = 0.d0
      eZnnR(:)  = 0.d0
      nz(:)     = 0.d0
      hop_plus(:) = 0.d0
      hop_minus(:) = 0.d0

      hop_ene   = 0.d0
      e_layers  = 0.d0
      e_orth    = 0.d0
    END SUBROUTINE init_to_zero

  END SUBROUTINE slater_step

        
  !+-----------------------------------------+!
  !  Gutzwiller projectors optimization step  !
  !  Build the p-h symmetric hamiltonian for  !
  !  the GZ projectors.                       !
  !  Diagonalize and compute the Hubbard      !
  !  energy                                   !
  !+-----------------------------------------+!
  SUBROUTINE projector_step
    complex(8),dimension(:),allocatable :: H0,H1,H2
    complex(8),dimension(:),allocatable :: S0,S1,S2
    real(8) :: S_tot

    allocate(H0(L),H1(L),H2(L))
    allocate(S0(L),S1(L),S2(L))
    
    call build_Hphi_equ
    call diag_Hphi_equ

    deallocate(H0,H1,H2,S0,S1,S2)
    
  CONTAINS
    
    SUBROUTINE build_Hphi_equ
      integer :: i
      
      complex(8) :: p0,p1,p2

      S_tot=0.d0

      do i=1,L
         
         p0 = opt_Phi(i)%p0
         p1 = opt_Phi(i)%p1
         p2 = opt_Phi(i)%p2

         H0(i) = Uz(i)*.5
         H2(i) = Uz(i)*.5
         H1(i) = conjg(r(i))*eZ(i)
         
         if(i.gt.1) then
            H1(i) = H1(i) + conjg(r(i-1))*hop_minus(i-1)*Hslab(i-1,i)
         end if
         
         if(i.lt.L) then
            H1(i) = H1(i) + conjg(r(i+1))*hop_plus(i)*Hslab(i+1,i)
         end if
         H1(i) = H1(i)*sqrt(2.d0)
         
         ! S0(i) = log(4.d0*ABS(eqPhi(i)%p0)**2) + 1.d0
         ! S1(i) = log(2.d0*ABS(eqPhi(i)%p1)**2) + 1.d0
         ! S2(i) = log(4.d0*ABS(eqPhi(i)%p2)**2) + 1.d0
         
         ! S_tot = S_tot + abs(p0)**2*log(abs(p0)**2/4.d0) + &
         !      + abs(p2)**2*log(abs(p2)**2/4.d0) + &
         !      + abs(p1)**2*log(abs(p1)**2/2.d0)


      end do

      
    END SUBROUTINE build_Hphi_equ
    
    SUBROUTINE diag_Hphi_equ
      complex(8),dimension(3,3) :: GZ_H
      real(8),dimension(3)      :: GZ_EIGEN
      integer :: i,j,i_L
      character(len=1) ::  jobz
      character(len=1) ::  uplo
      jobz = 'V'
      uplo = 'L'
      
      do i_L=1,L
         !  build p-h symmetric GZ HAMILTONIAN  !
         GZ_H(1,1) = H0(i_L) 
         GZ_H(1,2) = H1(i_L)
         GZ_H(1,3) = Z0

         GZ_H(2,1) = conjg(H1(i_L))
         GZ_H(2,2) = 0.d0 
         GZ_H(2,3) = H1(i_L)


         GZ_H(3,1) = Z0
         GZ_H(3,2) = conjg(H1(i_L))
         GZ_H(3,3) = H2(i_L) 
         
         call matrix_diagonalize(GZ_H,GZ_eigen,jobz,uplo)
         
         opt_Phi(i_L)%p0 = GZ_H(1,1)
         opt_Phi(i_L)%p1 = GZ_H(2,1)
         opt_Phi(i_L)%p2 = GZ_H(3,1)

      end do
      
      !  compute new wfc renormalization factors  !
      r   = gz_hop(opt_Phi)
      dop = gz_doping(opt_Phi)
      hubb_ene = slab_hubb_ene(opt_Phi,L,Uz)
      
    END SUBROUTINE diag_Hphi_equ

    
  END SUBROUTINE projector_step
      
  
END MODULE SLAB_EQU
     
