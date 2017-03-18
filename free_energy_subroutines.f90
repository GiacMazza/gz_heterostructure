subroutine free_energy_GZ_parameters( x, f, i )
  use global
  !use slab_equ
  implicit none
  !+- routine variables -+!
  real(8), intent( in )           :: x( : )
  real(8), intent( out )          :: f
  integer, intent( in ), optional :: i

  !+- inner variables -+!
  integer                         :: iL,jL
  type(gz_projector),dimension(L) :: phi_gz
  complex(8),dimension(L)         :: Rhop
  real(8),dimension(L)            :: Docc,xgz


  !+- slater determinant optimization -+!
  complex(8),dimension(L,L)       :: H_star
  real(8),dimension(L)            :: w
  !character(len=1)        :: jobz
  !character(len=1)        :: uplo

  type(vec2D) :: k
  integer                         :: ik
  real(8)  :: ek
  integer  :: i_theta,i_gm,i_phi

  real(8) :: theta,gm,phi
  real(8) :: e_hop,e_hubb,entropy
  !real(8),dimension(L) :: s_star
  real(8) :: s_star_slab
  complex(8) :: p0,p1,p2
  real(8)                                     :: n_temp
  real(8)  :: free_entropy
  
  do iL=1,L
     phi_gz(iL)%p0 = x((iL-1)*3 + 1)
     phi_gz(iL)%p2 = x((iL-1)*3 + 2)
     phi_gz(iL)%p1 = x((iL-1)*3 + 3)
  end do


  if ( .not. present( i ) ) then
     !+- FREE ENERGY ESTIMATION WITHIN GZ APPROX -+!


     !+- GZ parameter reconstruction -+!
     Rhop = GZ_hop(phi_gz)
     Docc = GZ_double(phi_gz)
     xgz  = GZ_doping(phi_gz) 
     !+-------------------------------+!


     !+- minimize F_star with respect to layer dependent chemical potentials -+!
     mu_star = 0.d0  
     call fzero_broyden(nstar,mu_star)
     do iL=1,L
        write(10,'(10(F18.10))') dble(iL),mu_star(iL),tmp_dop(iL),abs(Rhop(iL))
     end do
     !+-----------------------------------------------------------------------+!

     !+- trace over uncorrelated distribution -+!
     hop_intra    = 0.d0
     delta_plus   = 0.d0
     delta_minus  = 0.d0
     free_entropy = 0.d0
     do ik=1,Nk_tot !+- k points loop -+!
        !+- build and diagonalize real space hamiltonian -+!
        k = vec_k(ik)
        ek = square_lattice_disp(k)
        H_star=0.d0
        w = 0.d0
        do iL=1,L
           H_star(iL,iL) =  ek*abs(Rhop(iL))**2 + mu_star(iL)
           if(iL.lt.L) then
              H_star(iL,iL+1) = -1.d0*conjg(Rhop(iL))*Rhop(iL+1)
              H_star(iL+1,iL) = -1.d0*Rhop(iL)*conjg(Rhop(iL+1))
           end if
        end do
        if(pbc) then
           H_star(1,L) = -1.d0*conjg(Rhop(1))*Rhop(L)
           H_star(L,1) = -1.d0*conjg(Rhop(L))*Rhop(1)
        end if
        call  matrix_diagonalize(H_star,w,'V','L')
        !+------------------------------------------------+!
        do iL=1,L
           do jL=1,L
              hop_intra(jL) = hop_intra(jL) + 2.d0*ek*abs(H_star(jL,iL))**2*wt(ik)*fermi(w(iL),beta)
              if(jL.lt.L) then
                 delta_plus(jL)  = delta_plus(jL) + 2.d0*conjg(H_star(jL+1,iL))*H_star(jL,iL)*wt(ik)*fermi(w(iL),beta)
                 delta_minus(jL) = delta_minus(jL) + 2.d0*conjg(H_star(jL,iL))*H_star(jL+1,iL)*wt(ik)*fermi(w(iL),beta)
                 !+- free electrons entropy -+!
                 if(jL.eq.1) then
                    n_temp = fermi(w(iL),beta)
                    if(n_temp.gt.1.d-10) then
                       free_entropy = free_entropy - ( n_temp*log(n_temp))*wt(ik)*2.d0 
                    end if
                    if(abs(1.d0-n_temp).gt.1.d-10) then
                       free_entropy = free_entropy - (1.d0-n_temp)*log(1.d0-n_temp )*wt(ik)*2.d0
                    end if
                 end if
                 !+--------------------------+!
              end if
           end do
        end do
     end do !+- end k points loop -+! 



     e_hop  = 0.d0
     e_hubb = 0.d0
     entropy = 0.d0
     do iL=1,L
        !+- hopping energy -+!
        e_hop = e_hop + hop_intra(iL)*abs(Rhop(iL))**2
        if (iL.lt.L) then
           e_hop = e_hop - conjg(Rhop(iL+1))*Rhop(iL)*delta_plus(iL)
        end if
        if(iL.gt.1) then
           e_hop = e_hop - conjg(Rhop(iL-1))*Rhop(iL)*delta_minus(iL-1)
        end if
        !+------------------+!

        !+- hubbard energy -+!
        e_hubb = e_hubb + Uz(iL)*Docc(iL)  + local_field(iL)*(1.d0 - xgz(iL))
        !+------------------+!


        p0 = 0.25d0*(1.d0+tmp_dop(iL))**2.d0
        p1 = 0.5d0*(1.d0-tmp_dop(iL)**2.d0)
        p2 = 0.25d0*(1.d0-tmp_dop(iL))**2.d0

        !+- projectors entropy contribution -+!
        if(abs(phi_gz(iL)%p1).gt.1.d-8) then
           entropy = entropy - abs(phi_gz(iL)%p0)**2*log(abs(phi_gz(iL)%p0)**2/p0)
        end if
        if(abs(phi_gz(iL)%p2).gt.1.d-8) then
           entropy = entropy - abs(phi_gz(iL)%p2)**2*log(abs(phi_gz(iL)%p2)**2/p2)
        end if
        if(abs(phi_gz(iL)%p1).gt.1.d-8) then
           entropy = entropy - abs(phi_gz(iL)%p1)**2*log(abs(phi_gz(iL)%p1)**2/p1)
        end if
        !+-----------------------------------+!
     end do

     !+- free energy -+!
     f = e_hop + e_hubb - Temp*(entropy + free_entropy)
     write(*,*) e_hop
     !+---------------+!

  else

     !+- CONSTRAINTS ON GUTZWILLER PARAMETERS -+!
     if(i.le.L) then
        f = abs(phi_gz(i)%p0)**2 - abs(phi_gz(i)%p2)**2 - tmp_dop(i)        !+----------------------------------------+!
     else
        f = abs(phi_gz(i-L)%p0)**2 + abs(phi_gz(i-L)%p2)**2 + abs(phi_gz(i-L)%p1)**2 - 1.d0
     end if

  end if


contains



  FUNCTION NSTAR(mu)
    USE global
    real(8),dimension(:) :: mu
    real(8),dimension(size(mu)) :: nstar


    real(8) :: fstar_lagrange
    real(8) :: entropy,energy
    real(8) :: n_temp

    hop_intra   = 0.d0
    delta_plus  = 0.d0
    delta_minus = 0.d0
    nstar = 0.d0
    do ik=1,Nk_tot !+- k points loop -+!
       k = vec_k(ik)
       ek = square_lattice_disp(k)
       H_star=0.d0
       w = 0.d0
       do iL=1,L
          H_star(iL,iL) =  ek*abs(Rhop(iL))**2 + mu(iL)
          if(iL.lt.L) then
             H_star(iL,iL+1) = -1.d0*conjg(Rhop(iL))*Rhop(iL+1)
             H_star(iL+1,iL) = -1.d0*Rhop(iL)*conjg(Rhop(iL+1))
          end if
       end do
       if(pbc) then
          H_star(1,L) = -1.d0*conjg(Rhop(1))*Rhop(L)
          H_star(L,1) = -1.d0*conjg(Rhop(L))*Rhop(1)
       end if
       call  matrix_diagonalize(H_star,w,'V','L')
       do iL = 1,L
          do jL=1,L
             nstar(jL) = nstar(jL) + abs(H_star(jL,iL))**2*fermi(w(iL),beta)*wt(ik)*2.d0
             hop_intra(jL) = hop_intra(jL) + 2.d0*ek*abs(H_star(jL,iL))**2*wt(ik)*fermi(w(iL),beta)
             if(jL.lt.L) then
                delta_plus(jL)  = delta_plus(jL) + 2.d0*conjg(H_star(jL+1,iL))*H_star(jL,iL)*wt(ik)*fermi(w(iL),beta)
                delta_minus(jL) = delta_minus(jL) + 2.d0*conjg(H_star(jL,iL))*H_star(jL+1,iL)*wt(ik)*fermi(w(iL),beta)
             end if
          end do
       end do
    end do !+- end k points loop -+! 

    nstar = nstar - (1.d0 - tmp_dop)

    return
  END FUNCTION NSTAR




end subroutine free_energy_GZ_parameters







subroutine fun ( x, f, i )
  !.............................................................................
  use global
  real(8), intent( in )   :: x( : )
  real(8), intent( out )  :: f
  integer, intent( in ), optional   :: i

  if ( .not. present( i ) ) then
     !       the objective function value (user defined)
     !==============================================================================
     f = 100.0d0*(x(2)-x(1)**2)**2 +(1.0d0-x(1))**2                      !
     !==============================================================================
  else
     select case ( i )
     case ( 1 )
        !               the equality constraint value (user defined)
        !==============================================================================
        f = x(1)+3.0d0*x(2)-3.0d0                                   !
        !==============================================================================
     case ( 2 )
        !               the inequality constraint value (user defined)
        !==============================================================================
        f = x(1)**2+x(2)**2-4.0d0                                    !
        !==============================================================================
     end select
  end if
  return
end subroutine fun



