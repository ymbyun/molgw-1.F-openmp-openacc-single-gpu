!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! basic tools to play with the calculated self-energies
!
!=========================================================================
module m_selfenergy_tools
 use m_definitions
 use m_warning
 use m_mpi
 use m_inputparam

 !
 ! frozen core approximation parameters
 integer,protected :: ncore_G
 integer,protected :: nvirtual_G

 !
 ! Range of states to evaluate the self-energy
 integer,protected :: nsemin
 integer,protected :: nsemax

 !
 ! Highest occupied state
 integer,protected :: nhomo_G

 !
 ! Selfenergy evaluated on a frequency grid
 type selfenergy_grid
   integer                 :: nomega
   integer                 :: nomegai
   complex(dp),allocatable :: omega(:)
   complex(dp),allocatable :: omegai(:)
   real(dp),allocatable    :: energy0(:,:)
   complex(dp),allocatable :: sigma(:,:,:)
   complex(dp),allocatable :: sigmai(:,:,:)
 end type


contains


!=========================================================================
subroutine selfenergy_set_state_range(nstate,occupation)
 implicit none
 integer             :: nstate
 real(dp),intent(in) :: occupation(:,:)
!=====
 integer :: istate
!=====

 if( nstate > SIZE( occupation(:,:) , DIM=1 ) ) then
   call die('selfenergy_set_state_range: nstate is too large')
 endif

 ncore_G      = ncoreg
 nvirtual_G   = MIN(nvirtualg,nstate+1)

 if(is_frozencore) then
   if( ncore_G == 0) ncore_G = atoms_core_states()
 endif

 if( ncore_G > 0 ) then
   write(msg,'(a,i4,2x,i4)') 'frozen core approximation in G switched on up to state = ',ncore_G
   call issue_warning(msg)
 endif

 if( nvirtual_G <= nstate ) then
   write(msg,'(a,i4,2x,i4)') 'frozen virtual approximation in G switched on starting with state = ',nvirtual_G
   call issue_warning(msg)
 endif

 ! Find the HOMO index
 nhomo_G = 1
 do istate=1,nstate
   if( .NOT. ANY( occupation(istate,:) < completely_empty ) ) then
     nhomo_G = MAX(nhomo_G,istate)
   endif
 enddo

 nsemin = MAX(ncore_G+1,selfenergy_state_min,1,nhomo_G-selfenergy_state_range)
 nsemax = MIN(nvirtual_G-1,selfenergy_state_max,nstate,nhomo_G+selfenergy_state_range)

 write(stdout,'(a,i4,a,i4)') ' Calculate state range from ',nsemin,' to ',nsemax

end subroutine selfenergy_set_state_range


!=========================================================================
! ymbyun 2018/07/11
! The QP energy Eqp is obtained from the spectral function A(w) while A(w) is being written to a file.
#ifdef ENABLE_YMBYUN
 subroutine write_selfenergy_omega(filename_root,nstate,exchange_m_vxc,energy0,se, energy_qp_A)
#else
 subroutine write_selfenergy_omega(filename_root,nstate,exchange_m_vxc,energy0,se)
#endif
 implicit none

 character(len=*)    :: filename_root
 integer,intent(in)  :: nstate
 real(dp),intent(in) :: exchange_m_vxc(nstate,nspin)
 real(dp),intent(in) :: energy0(nstate,nspin)
 type(selfenergy_grid),intent(in) :: se

! ymbyun 2018/07/11
#ifdef ENABLE_YMBYUN
 real(dp),intent(out),optional    :: energy_qp_A(nstate,nspin)
#endif
!=====
 character(len=3)   :: ctmp
 character(len=256) :: filename
 integer :: selfenergyfile
 integer :: pstate
 integer :: iomega
 real(dp) :: spectral_function_w(nspin)

! ymbyun 2018/07/11
#ifdef ENABLE_YMBYUN
 integer  :: ii
 integer  :: pspin
 real(dp) :: Z_spe(nspin), Z_der(nspin)
 real(dp) :: Eqp(nspin), Eqp_left(nspin), Eqp_right(nspin)
 real(dp) :: A_max(nspin), A_left(nspin), A_right(nspin)
 real(dp) :: dSdw_all(nstate,nspin), dSdw(nspin)
#endif
!=====

 ! Just the master writes
 if( .NOT. is_iomaster ) return

 write(stdout,'(/,1x,a)') 'Write Sigma(omega) on file'

! ymbyun 2018/07/11
#ifdef ENABLE_YMBYUN
 if( se%nomega > 0 .AND. PRESENT(energy_qp_A) ) then

! Too much effort for too little return
! First touch to reduce NUMA effects using memory affinity
! Both (positive) affinity and (negative) overhead occur at the same time.
! A performance test is needed to find which one has a stronger effect on the speedup.
!
! ymbyun 2019/03/11
! Let's keep it simple. This will make GnWn safe when not all eigenvalues are updated.
!
!#ifdef ENABLE_OPENMP_AFFINITY
!!$OMP PARALLEL
!!$OMP DO PRIVATE(pstate)
!   do pstate=nsemin,nsemax
!     energy_qp_A(pstate,:) = energy0(pstate,:)
!     dSdw_all(pstate,:) = 0.0_dp
!   enddo
!!$OMP END DO
!!$OMP END PARALLEL
!#else
!   energy_qp_A(:,:) = energy0(:,:)
!   dSdw_all(:,:) = 0.0_dp
!#endif

   energy_qp_A(:,:) = energy0(:,:)
   dSdw_all(:,:) = 0.0_dp
 endif
#endif

 !
 ! omega is defined with respect to energy0_a
 ! Absolute omega is omega + energy0_a
 !
 ! ymbyun 2018/07/11
 ! NOTE: I'm not sure if it's safe for multiple threads to write multiple files at the same time.
 !
!$OMP PARALLEL
#ifdef ENABLE_YMBYUN
!$OMP DO PRIVATE(ctmp,filename,selfenergyfile,pstate,iomega,spectral_function_w, pspin,Z_spe,Z_der,Eqp,Eqp_left,Eqp_right,A_max,A_left,A_right,dSdw)
#else
!$OMP DO PRIVATE(ctmp,filename,selfenergyfile,pstate,iomega,spectral_function_w)
#endif
 do pstate=nsemin,nsemax
   write(ctmp,'(i3.3)') pstate
   filename = TRIM(filename_root) // '_state' // TRIM(ctmp) // '.dat'
   write(stdout,'(1x,a,a)') 'Writing selfenergy in file: ', TRIM(filename)
   open(newunit=selfenergyfile,file=filename)

! ymbyun 2019/09/26
#ifdef ENABLE_YMBYUN
   write(selfenergyfile,'(a)') '#       omega (eV)          Re SigmaC (eV)     Im SigmaC (eV)    omega - e_gKS - Vxc + SigmaX (eV)     A (eV^-1)     Z(omega)     dSigma/domega'
#else
   write(selfenergyfile,'(a)') '#       omega (eV)          Re SigmaC (eV)     Im SigmaC (eV)    omega - e_gKS - Vxc + SigmaX (eV)     A (eV^-1)'
#endif

#ifdef ENABLE_YMBYUN
   if( se%nomega > 0 .AND. PRESENT(energy_qp_A) ) then
     A_max(:) = 0.0_dp
   endif
#endif

   do iomega=-se%nomega,se%nomega
     spectral_function_w(:) = 1.0_dp / pi * ABS(   &
                                       AIMAG( 1.0_dp   &
                                         / ( se%energy0(pstate,:)+se%omega(iomega) - energy0(pstate,:)  &
                                               - exchange_m_vxc(pstate,:) - se%sigma(iomega,pstate,:) ) ) )

! ymbyun 2019/09/26
! PT2 uses real numbers, not complex numbers. So, we use Z(w) instead of A(w) to address the multi-solution issue from the graphical-solution method.
#ifdef ENABLE_YMBYUN
     if( iomega == -se%nomega ) then
       dSdw(:) = ( REAL(se%sigma(iomega+1,pstate,:),dp) - REAL(se%sigma(iomega-0,pstate,:),dp) ) / ( se%omega(iomega+1) - se%omega(iomega-0) )
     else if( iomega == se%nomega ) then
       dSdw(:) = ( REAL(se%sigma(iomega+0,pstate,:),dp) - REAL(se%sigma(iomega-1,pstate,:),dp) ) / ( se%omega(iomega+0) - se%omega(iomega-1) )
     else
       dSdw(:) = ( REAL(se%sigma(iomega+1,pstate,:),dp) - REAL(se%sigma(iomega-1,pstate,:),dp) ) / ( se%omega(iomega+1) - se%omega(iomega-1) )
     endif
     Z_der(:) = 1.0_dp / ( 1.0_dp - dSdw(:) )

     write(selfenergyfile,'(20(f16.8,2x))') ( se%omega(iomega) + se%energy0(pstate,:) )*Ha_eV,             &
                                            REAL(se%sigma(iomega,pstate,:),dp) * Ha_eV,                    &
                                            AIMAG(se%sigma(iomega,pstate,:)) * Ha_eV,                      &
                                            (REAL(se%omega(iomega),dp)+se%energy0(pstate,:) - energy0(pstate,:) - exchange_m_vxc(pstate,:) )*Ha_eV, &
                                            spectral_function_w(:) / Ha_eV,                                &
                                            Z_der(:),                                                      &
                                            dSdw(:)
#else
     write(selfenergyfile,'(20(f16.8,2x))') ( se%omega(iomega) + se%energy0(pstate,:) )*Ha_eV,             &
                                            REAL(se%sigma(iomega,pstate,:),dp) * Ha_eV,                    &
                                            AIMAG(se%sigma(iomega,pstate,:)) * Ha_eV,                      &
                                            (REAL(se%omega(iomega),dp)+se%energy0(pstate,:) - energy0(pstate,:) - exchange_m_vxc(pstate,:) )*Ha_eV, &
                                            spectral_function_w(:) / Ha_eV
#endif

#ifdef ENABLE_YMBYUN
     if( se%nomega > 0 .AND. PRESENT(energy_qp_A) ) then
       if( iomega /= -se%nomega .AND. iomega /= se%nomega ) then
         ! NOTE: Z_spe and Z_der here are obtained by using the discrete w grid. A fine w grid makes them more accurate.
         Eqp(:)   = REAL(se%omega(iomega),dp) + se%energy0(pstate,:)
         Z_spe(:) = (Eqp(:)+energy0(pstate,:)*(-1.0_dp)) / (exchange_m_vxc(pstate,:)+REAL(se%sigma(0,pstate,:),dp))

         dSdw(:)  = ( REAL(se%sigma(iomega+1,pstate,:),dp) - REAL(se%sigma(iomega-1,pstate,:),dp) ) / ( se%omega(iomega+1) - se%omega(iomega-1) )
         Z_der(:) = 1.0_dp / ( 1.0_dp - dSdw(:) )

         do pspin=1,nspin
           if( spectral_function_w(pspin) > A_max(pspin) ) then
             ! NOTE: I'm not sure which one to use between Z_spe and Z_der. They yield identical results for energy levels near HOMO/LUMO.
             if( Z_spe(pspin) > 0.0_dp .AND. Z_spe(pspin) < 1.0_dp ) then
!             if( Z_der(pspin) > 0.0_dp .AND. Z_der(pspin) < 1.0_dp ) then
               ! First, spectral function A(w)
               A_max(pspin) = spectral_function_w(pspin)

               ! Second, quasi-particle energy Eqp
               A_left(pspin)  = 1.0_dp / pi * ABS(   &
                                       AIMAG( 1.0_dp   &
                                         / ( se%energy0(pstate,pspin)+se%omega(iomega-1) - energy0(pstate,pspin)  &
                                               - exchange_m_vxc(pstate,pspin) - se%sigma(iomega-1,pstate,pspin) ) ) )
               A_right(pspin) = 1.0_dp / pi * ABS(   &
                                       AIMAG( 1.0_dp   &
                                         / ( se%energy0(pstate,pspin)+se%omega(iomega+1) - energy0(pstate,pspin)  &
                                               - exchange_m_vxc(pstate,pspin) - se%sigma(iomega+1,pstate,pspin) ) ) )
               Eqp_left(pspin)  = REAL(se%omega(iomega-1),dp) + se%energy0(pstate,pspin)
               Eqp_right(pspin) = REAL(se%omega(iomega+1),dp) + se%energy0(pstate,pspin)
!               energy_qp_A(pstate,pspin) = Eqp(pspin)
               energy_qp_A(pstate,pspin) = 1.0_dp / ( A_left(pspin) + A_max(pspin) + A_right(pspin) ) &
                                                  * ( A_left(pspin)*Eqp_left(pspin) + A_max(pspin)*Eqp(pspin) + A_right(pspin)*Eqp_right(pspin) )

               ! Thrid, derivative of Sigma(omega)
               dSdw_all(pstate,pspin) = dSdw(pspin)
!               dSdw_all(pstate,pspin) =   ( REAL(se%sigma(iomega+1,pstate,pspin),dp) - REAL(se%sigma(iomega-1,pstate,pspin),dp) ) &
!                                        / ( se%omega(iomega+1) - se%omega(iomega-1) )
             endif
           endif
         enddo
       endif
     endif
#endif
   enddo
   if( se%nomegai > 0 ) then
     do iomega=-se%nomegai,se%nomegai
       write(selfenergyfile,'(20(f16.8,2x))') ( se%omegai(iomega) + se%energy0(pstate,:) )*Ha_eV,     &
                                              REAL(se%sigmai(iomega,pstate,:),dp) * Ha_eV,            &
                                              AIMAG(se%sigmai(iomega,pstate,:)) * Ha_eV,              &
                                              0.0_dp,0.0_dp
     enddo
   endif
   write(selfenergyfile,*)
   close(selfenergyfile)

! ymbyun 2018/07/11
! This happens for energy levels near HOMO/LUMO + large step_sigma.
#ifdef ENABLE_YMBYUN
   if( se%nomega > 0 .AND. PRESENT(energy_qp_A) ) then
     do pspin=1,nspin
       if( A_max(pspin) == 0.0_dp ) then
         write(stdout,'(a,i4,2x,i4)') ' (ymbyun) The spectral function A(w) cannot yield the QP energy Eqp for state,spin: ',pstate,pspin
         call issue_warning('(ymbyun) The QP energy Eqp cannot be obtained from the spectral function A(w)')
       endif
     enddo
   endif
#endif

 enddo
!$OMP END DO
!$OMP END PARALLEL

! ymbyun 2018/07/11
! For (i) the derivative of Sigma(omega) at Eqp and (ii) the high precision
#ifdef ENABLE_YMBYUN
 if( se%nomega > 0 .AND. PRESENT(energy_qp_A) ) then
   write(stdout,*)
   write(stdout,'(a)') ' (ymbyun) QP energies Eqp (eV) obtained from spectral function A(w)'

   if( nspin == 1 ) then
     write(stdout,'(3x,a, 7x,a, 7x,a, 5x,a, 6x,a, 6x,a, 6x,a, 5x,a)') '#','E0','Eqp^spe','Eqp-E0','Z^spe','Z^der','dS/dw','ReSigC'
   else
     write(stdout,'(3x,a,13x,a,17x,a,18x,a,16x,a,16x,a,16x,a,16x,a)') '#','E0','Eqp^spe','Eqp-E0','Z^spe','Z^der','dS/dw','ReSigC'
     write(stdout,'(9x,16(a6,5x))') ('  up  ',' down ',ii=1,7)
   endif

   do pstate=nsemin,nsemax
     write(stdout,'(i4,a,20(1x,f10.4))')  pstate,                                                    &
                                          'y',                                                       & ! y means Young-Moo.
                                          energy0(pstate,:)*Ha_eV,                                   &
                                          energy_qp_A(pstate,:)*Ha_eV,                               &
                                          (energy_qp_A(pstate,:)+energy0(pstate,:)*(-1.0_dp))*Ha_eV, &
                                          (energy_qp_A(pstate,:)+energy0(pstate,:)*(-1.0_dp))/(exchange_m_vxc(pstate,:)+(REAL(se%sigma(0,pstate,:),dp))), &
                                          1.0_dp / ( 1.0_dp - dSdw_all(pstate,:) ),                  &
                                          dSdw_all(pstate,:),                                        &
                                          REAL(se%sigma(0,pstate,:)) * Ha_eV                           ! ymbyun 2018/10/22 Re{SigC} at Eks for the linearization method
   enddo
 endif
#endif

end subroutine write_selfenergy_omega


!=========================================================================
subroutine find_qp_energy_linearization(se,nstate,exchange_m_vxc,energy0,energy_qp_z,zz)
 implicit none

 type(selfenergy_grid),intent(in) :: se
 integer,intent(in)               :: nstate
 real(dp),intent(in)              :: exchange_m_vxc(nstate,nspin),energy0(nstate,nspin)
 real(dp),intent(out)             :: energy_qp_z(nstate,nspin)
 real(dp),intent(out),optional    :: zz(nsemin:nsemax,nspin)
!=====
 integer  :: pstate,pspin
 real(dp) :: zz_p(nspin)
!=====

 ! First, a dummy initialization
! ymbyun 2018/07/12
! First touch to reduce NUMA effects using memory affinity
!
! ymbyun 2019/03/12
! Let's keep it simple.
!
!#ifdef ENABLE_OPENMP_AFFINITY
!!$OMP PARALLEL
!!$OMP DO PRIVATE(pstate)
! do pstate=nsemin,nsemax
!   energy_qp_z(pstate,:) = energy0(pstate,:)
! enddo
!!$OMP END DO
!#else
 energy_qp_z(:,:) = energy0(:,:)
!!$OMP PARALLEL
!#endif

 ! Then overwrite the interesting energy with the calculated GW one
! ymbyun 2018/07/12
!$OMP PARALLEL
!$OMP DO PRIVATE(pstate,pspin,zz_p)
 do pstate=nsemin,nsemax

   if( se%nomega > 0 .AND. PRESENT(zz) ) then
     zz_p(:) = ( REAL(se%sigma(1,pstate,:),dp) - REAL(se%sigma(-1,pstate,:),dp) ) / ( se%omega(1) - se%omega(-1) )
     zz_p(:) = 1.0_dp / ( 1.0_dp - zz_p(:) )
     ! Conrain Z to be in [0:1] to avoid crazy values
     ! ymbyun 2018/07/02
     ! Z falls out of [0:1] when a weak self-energy pole is very close to Eks.
     ! Z out of [0:1] is an indicator for whether it happened or not.
#ifndef ENABLE_YMBYUN
     do pspin=1,nspin
       zz_p(pspin) = MIN( MAX(zz_p(pspin),0.0_dp) , 1.0_dp )
     enddo
#endif

     zz(pstate,:)          = zz_p(:)
     energy_qp_z(pstate,:) = se%energy0(pstate,:)  &
             + zz_p(:) * ( energy0(pstate,:) - se%energy0(pstate,:)  &
                            + REAL(se%sigma(0,pstate,:),dp) + exchange_m_vxc(pstate,:) )
   else
     energy_qp_z(pstate,:) = energy0(pstate,:) + REAL(se%sigma(0,pstate,:),dp) + exchange_m_vxc(pstate,:)
   endif

 enddo
!$OMP END DO
!$OMP END PARALLEL

end subroutine find_qp_energy_linearization


!=========================================================================
subroutine find_qp_energy_graphical(se,nstate,exchange_m_vxc,energy0,energy_qp_g)
 implicit none

 type(selfenergy_grid),intent(in) :: se
 integer,intent(in)   :: nstate
 real(dp),intent(in)  :: exchange_m_vxc(nstate,nspin),energy0(nstate,nspin)
 real(dp),intent(out) :: energy_qp_g(nstate,nspin)
!=====
 integer  :: pstate,pspin
 integer  :: info
 real(dp) :: sigma_xc_m_vxc(-se%nomega:se%nomega)
!=====

 ! First, a dummy initialization
! ymbyun 2018/07/12
! First touch to reduce NUMA effects using memory affinity
!
! ymbyun 2019/03/12
! Let's keep it simple.
!
!#ifdef ENABLE_OPENMP_AFFINITY
!!$OMP PARALLEL
!!$OMP DO PRIVATE(pstate)
! do pstate=nsemin,nsemax
!   energy_qp_g(pstate,:) = 0.0_dp
! enddo
!!$OMP END DO
!#else
 energy_qp_g(:,:) = 0.0_dp
!!$OMP PARALLEL
!#endif

 ! Then overwrite the interesting energy with the calculated GW one
! ymbyun 2018/07/12
!$OMP PARALLEL
!$OMP DO PRIVATE(pstate,pspin,info,sigma_xc_m_vxc)
 do pstate=nsemin,nsemax

   if( MODULO(pstate-nsemin,nproc_world) /= rank_world ) cycle

   do pspin=1,nspin
     sigma_xc_m_vxc(:) = REAL(se%sigma(:,pstate,pspin),dp) + exchange_m_vxc(pstate,pspin)
     !
     ! QP equation:
     ! E_GW = E0 + \omega =  E_gKS + \Sigma_c(E0+\omega) + \Sigma_x - v_xc
     !
     energy_qp_g(pstate,pspin) = find_fixed_point(se%nomega,REAL(se%omega,dp)+se%energy0(pstate,pspin),sigma_xc_m_vxc(:)+energy0(pstate,pspin),info)
     if( info /= 0 ) then
       write(stdout,'(1x,a,i4,2x,i4)') 'Unreliable graphical solution of the QP equation for state,spin: ',pstate,pspin
       call issue_warning('No fixed point found for the QP equation. Try to increase nomega_sigma or step_sigma.')
     endif
   enddo

 enddo
!$OMP END DO
!$OMP END PARALLEL

 call xsum_world(energy_qp_g)

 energy_qp_g(:nsemin-1,:) = energy0(:nsemin-1,:)
 energy_qp_g(nsemax+1:,:) = energy0(nsemax+1:,:)

end subroutine find_qp_energy_graphical


!=========================================================================
function find_fixed_point(nx,xx,fx,info) result(fixed_point)
 implicit none
 integer,intent(in)  :: nx
 real(dp),intent(in) :: xx(-nx:nx)
 real(dp),intent(in) :: fx(-nx:nx)
 integer,intent(out) :: info
 real(dp)            :: fixed_point
!=====
 integer             :: ix,imin1,imin2
 real(dp)            :: rmin
 real(dp)            :: gx(-nx:nx)
 real(dp)            :: gpx
!=====

 !
 ! g(x) contains f(x) - x
 gx(:) = fx(:) - xx(:)

 ! Find the sign change in g(x) which is the closest to ix=0
 ! Search positive x
 imin1 = 1000000
 do ix=0,nx-1
   if( gx(ix) * gx(ix+1) < 0.0_dp ) then
     imin1 = ix
     exit
   endif
 enddo
 ! Search negative x
 imin2 = 1000000
 do ix=0,-nx+1,-1
   if( gx(ix) * gx(ix-1) < 0.0_dp ) then
     imin2 = ix
     exit
   endif
 enddo

 if( imin1 == 1000000 .AND. imin2 == 1000000 ) then

   if( gx(0) > 0.0_dp ) then 
     info =  1
     fixed_point = xx(nx)
   else
     info = -1
     fixed_point = xx(-nx)
   endif

 else
   info = 0
   ! 
   ! If sign change found, interpolate the abscissa between 2 grid points
   if( ABS(imin1) <= ABS(imin2) )  then
     gpx = ( gx(imin1+1) - gx(imin1) ) / ( xx(imin1+1) - xx(imin1) )
     fixed_point = xx(imin1) - gx(imin1) / gpx 
   else
     gpx = ( gx(imin2) - gx(imin2-1) ) / ( xx(imin2) - xx(imin2-1) )
     fixed_point = xx(imin2-1) - gx(imin2-1) / gpx 
   endif
 endif


end function find_fixed_point


!=========================================================================
! ymbyun 2018/07/11
! The QP energy Eqp is obtained by three different methods:
!   1. Linearization
!   2. Graph
!   3. Spectral function A(w)
#ifdef ENABLE_YMBYUN
 subroutine output_qp_energy(calcname,nstate,energy0,exchange_m_vxc,ncomponent,se,energy1,energy2,zz, energy3)
#else
 subroutine output_qp_energy(calcname,nstate,energy0,exchange_m_vxc,ncomponent,se,energy1,energy2,zz)
#endif
 implicit none
 
 character(len=*)             :: calcname
 integer                      :: nstate,ncomponent
 real(dp),intent(in)          :: energy0(nstate,nspin),exchange_m_vxc(nstate,nspin)
 type(selfenergy_grid),intent(in) :: se
 real(dp),intent(in)          :: energy1(nstate,nspin)
 real(dp),intent(in),optional :: energy2(nstate,nspin),zz(nsemin:nsemax,nspin)

! ymbyun 2018/07/11
#ifdef ENABLE_YMBYUN
 real(dp),intent(in),optional :: energy3(nstate,nspin)
#endif
!=====
 integer           :: pstate,ii
 character(len=36) :: sigc_label
!=====

 if( ncomponent > 2 ) call die('output_qp_energy: too many components. Not implemented yet')
 if( ncomponent < 1 ) call die('output_qp_energy: too few components. Something is not correct in the coding.')

 if( ncomponent == 1 ) then
   sigc_label = 'SigC'
 else
   if( nspin == 1 ) then
     sigc_label = 'SigC_1      SigC_2'
   else
     sigc_label = 'SigC_1                  SigC_2'
   endif
 endif

 write(stdout,'(/,1x,a,1x,a)') TRIM(calcname),'eigenvalues (eV)'

! ymbyun 2018/03/13
#ifdef ENABLE_YMBYUN
 if( PRESENT(zz) .AND. PRESENT(energy2) .AND. PRESENT(energy3) ) then
   if(nspin==1) then
     write(stdout,'(3x,a,5x,a,4x,a,5x,a,4x,a,4x,a,4x,a,2x,a,2x,a,2x,a,3x,a)') &
                  '#','E0','SigX-Vxc',TRIM(sigc_label),'Z^lin','Z^gra','Z^spe','Eqp^lin','Eqp^gra','Eqp^spe','Eqp-E0'
   else
     write(stdout,'(3x,a,9x,a,15x,a,12x,a,14x,a,13x,a,13x,a,11x,a,11x,a,11x,a,11x,a)') &
                  '#','E0','SigX-Vxc',TRIM(sigc_label),'Z^lin','Z^gra','Z^spe','Eqp^lin','Eqp^gra','Eqp^spe','Eqp-E0'
     write(stdout,'(9x,20(a4,5x))') (' up ','down',ii=1,9+ncomponent)
   endif

   do pstate=nsemin,nsemax
     write(stdout,'(i4,a,24(1x,f8.2))')   pstate,                               &
                                          'f',                                  & ! f means Fabien.
                                          energy0(pstate,:)*Ha_eV,              &
                                          exchange_m_vxc(pstate,:)*Ha_eV,       &
                                          REAL(se%sigma(0,pstate,:)*Ha_eV,dp),  &
                                          zz(pstate,:),                         &
                                          (energy2(pstate,:)+energy0(pstate,:)*(-1.0_dp))/(exchange_m_vxc(pstate,:)+(REAL(se%sigma(0,pstate,:),dp))),  &
                                          (energy3(pstate,:)+energy0(pstate,:)*(-1.0_dp))/(exchange_m_vxc(pstate,:)+(REAL(se%sigma(0,pstate,:),dp))),  &
                                          energy1(pstate,:)*Ha_eV,              &
                                          energy2(pstate,:)*Ha_eV,              &
                                          energy3(pstate,:)*Ha_eV,              &
                                          (energy3(pstate,:)+energy0(pstate,:)*(-1.0_dp))*Ha_eV
   enddo
 else
   if(nspin==1) then
     write(stdout,'(3x,a, 7x,a, 6x,a, 7x,a, 4x,a, 5x,a)') '#','E0','SigX-Vxc',TRIM(sigc_label),'Eqp^Z=1','Eqp-E0'
   else
     write(stdout,'(3x,a,12x,a,17x,a,16x,a,17x,a,16x,a)') '#','E0','SigX-Vxc',TRIM(sigc_label),'Eqp^Z=1','Eqp-E0'
     write(stdout,'(9x,16(a4,7x))') (' up ','down',ii=1,4+ncomponent)
   endif

   do pstate=nsemin,nsemax
     write(stdout,'(i4,a,20(1x,f10.4))')  pstate,                               &
                                          'f',                                  & ! f means Fabien.
                                          energy0(pstate,:)*Ha_eV,              &
                                          exchange_m_vxc(pstate,:)*Ha_eV,       &
                                          REAL(se%sigma(0,pstate,:)*Ha_eV,dp),  &
                                          energy1(pstate,:)*Ha_eV,              &
                                          (energy1(pstate,:)+energy0(pstate,:)*(-1.0_dp))*Ha_eV
   enddo
 endif
#else
 if( PRESENT(zz) .AND. PRESENT(energy2) ) then
   if(nspin==1) then
     write(stdout,'(3x,a,8x,a,9x,a,7x,a,10x,a,8x,a,5x,a)') '#','E0','SigX-Vxc',TRIM(sigc_label),'Z','E_qp^lin','E_qp^graph'
   else
     write(stdout,'(3x,a,15x,a,22x,a,19x,a,24x,a,21x,a,17x,a)') '#','E0','SigX-Vxc',TRIM(sigc_label),'Z','E_qp^lin','E_qp^graph'
     write(stdout,'(12x,14(a4,9x))') (' up ','down',ii=1,5+ncomponent)
   endif

   do pstate=nsemin,nsemax
     write(stdout,'(i4,1x,20(1x,f12.6))') pstate,                          &
                                          energy0(pstate,:)*Ha_eV,         &
                                          exchange_m_vxc(pstate,:)*Ha_eV,  &
                                          REAL(se%sigma(0,pstate,:)*Ha_eV,dp),  &
                                          zz(pstate,:),                    &
                                          energy1(pstate,:)*Ha_eV,         &
                                          energy2(pstate,:)*Ha_eV
   enddo
 else
   if(nspin==1) then
     write(stdout,'(3x,a,8x,a,9x,a,7x,a,9x,a)') '#','E0','SigX-Vxc',TRIM(sigc_label),'E_qp'
   else
     write(stdout,'(3x,a,15x,a,22x,a,20x,a,22x,a)') '#','E0','SigX-Vxc',TRIM(sigc_label),'E_qp'
     write(stdout,'(12x,10(a4,9x))') (' up ','down',ii=1,3+ncomponent)
   endif

   do pstate=nsemin,nsemax
     write(stdout,'(i4,1x,20(1x,f12.6))') pstate,energy0(pstate,:)*Ha_eV,       &
                                          exchange_m_vxc(pstate,:)*Ha_eV,       &
                                          REAL(se%sigma(0,pstate,:)*Ha_eV,dp),  &
                                          energy1(pstate,:)*Ha_eV
   enddo
 endif
#endif

end subroutine output_qp_energy


!=========================================================================
subroutine init_selfenergy_grid(selfenergy_technique,nstate,energy0,se)
 use m_atoms
 implicit none

 integer,intent(in)                  :: selfenergy_technique,nstate
 real(dp),intent(in)                 :: energy0(nstate,nspin)
 type(selfenergy_grid),intent(inout) :: se
!=====
 integer            :: iomega,pstate
!=====

 se%nomegai = 0
 se%nomega  = 0

 select case(selfenergy_technique)
 
 case(EVSC,static_selfenergy)

   allocate(se%omega(-se%nomega:se%nomega))
   se%omega(0) = 0.0_dp

 case(imaginary_axis)
   !
   ! Set the final sampling points for Sigma
   se%nomega = nomega_sigma/2
   allocate(se%omega(-se%nomega:se%nomega))
   do iomega=-se%nomega,se%nomega
     se%omega(iomega) = step_sigma * iomega
   enddo

   !
   ! Set the calculated sampling points for Sigma
   se%nomegai = nomega_sigma / 2
   allocate(se%omegai(-se%nomegai:se%nomegai))
   do iomega=-se%nomegai,se%nomegai
     se%omegai(iomega) = step_sigma * iomega * im * 3.0_dp
   enddo


 case(one_shot)
   select case(calc_type%selfenergy_approx)
   case(GSIGMA)
     se%nomega = 1

   case default
     se%nomega = nomega_sigma/2
     allocate(se%omega(-se%nomega:se%nomega))
     do iomega=-se%nomega,se%nomega
       se%omega(iomega) = step_sigma * iomega
     enddo

   end select

 end select


 !
 ! Set the central point of the grid
 allocate(se%energy0(nsemin:nsemax,nspin))

 select case(selfenergy_technique)
 case(imaginary_axis)
   ! Find the center of the HOMO-LUMO gap
   forall(pstate=nsemin:nsemax)
     se%energy0(pstate,:) = 0.5_dp * ( energy0(nhomo_G,:) + energy0(nhomo_G+1,:) )
   end forall
 case default
   se%energy0(nsemin:nsemax,:) = energy0(nsemin:nsemax,:)
 end select

 !
 ! Set the central point of the grid
 allocate(se%sigma(-se%nomega:se%nomega,nsemin:nsemax,nspin))
 select case(selfenergy_technique)
 case(imaginary_axis)
   allocate(se%sigmai(-se%nomegai:se%nomegai,nsemin:nsemax,nspin))
 end select


end subroutine init_selfenergy_grid


!=========================================================================
subroutine destroy_selfenergy_grid(se)
 implicit none
 type(selfenergy_grid),intent(inout) :: se
!=====

 se%nomega  = 0
 se%nomegai = 0
 if( ALLOCATED(se%omega) )  deallocate(se%omega)
 if( ALLOCATED(se%omegai) ) deallocate(se%omegai)
 deallocate(se%energy0)
 deallocate(se%sigma)
 if( ALLOCATED(se%sigmai) ) deallocate(se%sigmai)

end subroutine destroy_selfenergy_grid


!=========================================================================
subroutine setup_exchange_m_vxc_diag(basis,nstate,occupation,energy,c_matrix,hamiltonian_fock,exchange_m_vxc_diag)
 use m_inputparam
 use m_basis_set
 use m_dft_grid
 use m_hamiltonian
 use m_hamiltonian_sca
 implicit none

 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: nstate
 real(dp),intent(in)        :: occupation(nstate,nspin)
 real(dp),intent(in)        :: energy(nstate,nspin)
 real(dp),intent(in)        :: c_matrix(basis%nbf,nstate,nspin)
 real(dp),intent(in)        :: hamiltonian_fock(basis%nbf,basis%nbf,nspin)
 real(dp),intent(out)       :: exchange_m_vxc_diag(nstate,nspin)
!=====
 integer,parameter    :: BATCH_SIZE=64
 integer              :: ispin,istate
 real(dp)             :: exc,eexx
 real(dp),allocatable :: occupation_tmp(:,:)
 real(dp),allocatable :: p_matrix_tmp(:,:,:)
 real(dp),allocatable :: hxc_val(:,:,:),hexx_val(:,:,:),hxmxc(:,:,:)
!=====

 !
 ! Testing the core/valence splitting
 !
 if(dft_core > 0) then
   if( alpha_hybrid_lr > 0.001 ) then
     call die('RSH not implemented yet')
   endif
   write(msg,'(a,i4,2x,i4)') 'DFT core-valence interaction switched on up to state = ',dft_core
   call issue_warning(msg)

   allocate(occupation_tmp(nstate,nspin))
   allocate(p_matrix_tmp(basis%nbf,basis%nbf,nspin))
   allocate(hxc_val(basis%nbf,basis%nbf,nspin))
   allocate(hexx_val(basis%nbf,basis%nbf,nspin))
   allocate(hxmxc(basis%nbf,basis%nbf,nspin))

   ! Override the occupation of the core electrons
   occupation_tmp(:,:) = occupation(:,:)
   occupation_tmp(1:dft_core,:) = 0.0_dp

   if( calc_type%is_dft ) then
     call init_dft_grid(basis,grid_level,dft_xc_needs_gradient,.FALSE.,BATCH_SIZE)
     call setup_density_matrix(basis%nbf,nstate,c_matrix,occupation_tmp,p_matrix_tmp)
     call dft_exc_vxc_batch(BATCH_SIZE,basis,nstate,occupation_tmp,c_matrix,hxc_val,exc)
     call destroy_dft_grid()
   endif

   if( .NOT. has_auxil_basis ) then
     call setup_exchange(print_matrix_,basis%nbf,p_matrix_tmp,hexx_val,eexx)
   else
     if( parallel_ham ) then
       call die('setup_exchange_m_vxc_diag: case not implemented')
     else
       call setup_exchange_ri(basis%nbf,nstate,occupation_tmp,c_matrix,p_matrix_tmp,hexx_val,eexx)
     endif
   endif

   hxc_val(:,:,:) = hxc_val(:,:,:) + alpha_hybrid * hexx_val(:,:,:)
   hxmxc(:,:,:) = hexx_val(:,:,:) - hxc_val(:,:,:) 

   deallocate(occupation_tmp,p_matrix_tmp)

   ! ymbyun 2018/07/12
   ! This part is not parallelized because of two reasons:
   !   1. The speedup by parallelization is very small.
   !   2. My OpenMP and OpenMP-MATMUL may conflict with each other.
   !
   ! Calculate the diagonal of the matrix Sigma_x - Vxc
   ! for the forthcoming GW corrections
   exchange_m_vxc_diag(:,:) = 0.0_dp
   do ispin=1,nspin
     do istate=1,nstate

        exchange_m_vxc_diag(istate,ispin) =  DOT_PRODUCT(  c_matrix(:,istate,ispin) , &
                                                MATMUL( hxmxc(:,:,ispin) , c_matrix(:,istate,ispin) ) )
     enddo
   enddo
   deallocate(hxc_val,hexx_val,hxmxc)

 else

   ! ymbyun 2018/07/12
   ! Same as above
   !
   ! Calculate the diagonal of the matrix Sigma_x - Vxc
   ! for the forthcoming GW corrections
   exchange_m_vxc_diag(:,:) = 0.0_dp
   do ispin=1,nspin
     do istate=1,nstate

        exchange_m_vxc_diag(istate,ispin) =  DOT_PRODUCT(  c_matrix(:,istate,ispin) , &
                                                MATMUL( hamiltonian_fock(:,:,ispin) , c_matrix(:,istate,ispin) ) ) &
                                              - energy(istate,ispin)
     enddo
   enddo

 endif


end subroutine setup_exchange_m_vxc_diag


!=========================================================================
subroutine apply_qs_approximation(nbf,nstate,s_matrix,c_matrix,selfenergy)
 implicit none

 integer,intent(in)     :: nbf,nstate
 real(dp),intent(in)    :: s_matrix(nbf,nbf),c_matrix(nbf,nstate,nspin)
 real(dp),intent(inout) :: selfenergy(nbf,nbf,nspin)
!=====
 integer :: ispin
!=====

 !
 ! Kotani's Hermitianization trick
 !
 do ispin=1,nspin
   selfenergy(:,:,ispin) = 0.5_dp * ( selfenergy(:,:,ispin) + TRANSPOSE(selfenergy(:,:,ispin)) )

   ! Transform the matrix elements back to the AO basis
   ! do not forget the overlap matrix S
   ! C^T S C = I
   ! the inverse of C is C^T S
   ! the inverse of C^T is S C
   selfenergy(:,:,ispin) = MATMUL( MATMUL( s_matrix(:,:) , c_matrix(:,nsemin:nsemax,ispin) ) , &
                             MATMUL( selfenergy(nsemin:nsemax,nsemin:nsemax,ispin), &
                               MATMUL( TRANSPOSE(c_matrix(:,nsemin:nsemax,ispin)), s_matrix(:,:) ) ) )

 enddo

end subroutine apply_qs_approximation
 

!=========================================================================
subroutine self_energy_fit(nstate,energy,se)
 use m_lbfgs
 implicit none

 integer,intent(in)                  :: nstate
 real(dp),intent(in)                 :: energy(nstate,nspin)
 type(selfenergy_grid),intent(inout) :: se
!=====
 integer :: pstate,pspin
 integer :: iomega
 integer :: ilbfgs,ipp,ii,info
 real(dp),parameter :: dq=1.0e-4_dp
 integer,parameter :: nlbfgs=100
 integer,parameter :: npp=6
 integer,parameter :: nparam=6
 real(dp)          :: coeff(nparam*npp)
 real(dp)          :: coeffdq(nparam*npp)
 real(dp)          :: dchi2(nparam*npp)
 real(dp)          :: chi2
 type(lbfgs_state) :: lbfgs_plan
!=====

 do pspin=1,nspin
   do pstate=nsemin,nsemax

     !
     ! Optimization: first guess
     do ipp=1,npp
       coeff(1+(ipp-1)*nparam) = 0.5_dp  * ipp
       coeff(2+(ipp-1)*nparam) = 0.5_dp  * ipp
       coeff(3+(ipp-1)*nparam) = 0.01_dp * ipp
       coeff(4+(ipp-1)*nparam) = 0.01_dp * ipp
       coeff(5+(ipp-1)*nparam) = 0.3_dp  / ipp
       coeff(6+(ipp-1)*nparam) = 0.3_dp  / ipp
     enddo
     chi2 = eval_chi2(coeff)

     call lbfgs_init(lbfgs_plan,nparam*npp,5,diag_guess=1.0_dp)
     do ilbfgs=1,nlbfgs
       write(stdout,*) 'chi2=',chi2
       do ipp=1,npp
         do ii=1,nparam
           coeffdq(:) = coeff(:)
           coeffdq(ii+(ipp-1)*nparam) = coeffdq(ii+(ipp-1)*nparam) + dq 
           dchi2(ii+(ipp-1)*nparam) = ( eval_chi2(coeffdq) - eval_chi2(coeff) ) / dq
         enddo
       enddo
       info = lbfgs_execute(lbfgs_plan,coeff,chi2,dchi2)
       chi2 = eval_chi2(coeff)
       if( chi2 < 5.0e-8_dp ) exit
     enddo
     call lbfgs_destroy(lbfgs_plan)
     write(stdout,'(/,1x,a)') '=========================='
     write(stdout,'(1x,a)') '   #           Coeff              Re Pole            Im Pole'
     do ipp=1,npp
       write(stdout,'(1x,i4,4(2x,f18.8))') 2*ipp-1,                   &
                                           coeff(5+(ipp-1)*nparam)**2, &
                                           coeff(1+(ipp-1)*nparam)**2, &
                                          -coeff(3+(ipp-1)*nparam)**2
       write(stdout,'(1x,i4,4(2x,f18.8))') 2*ipp,                     &
                                           coeff(6+(ipp-1)*nparam)**2, &
                                          -coeff(2+(ipp-1)*nparam)**2, &
                                           coeff(4+(ipp-1)*nparam)**2
     enddo
     write(stdout,'(1x,a)') '=========================='


     do iomega=-se%nomegai,se%nomegai
       write(500+pstate,'(6(1x,f18.8))') (se%omegai(iomega) + se%energy0(pstate,pspin))*Ha_eV, &
                                         se%sigmai(iomega,pstate,pspin)*Ha_eV, &
                                         eval_func(coeff, se%omegai(iomega) + se%energy0(pstate,pspin) )*Ha_eV
     enddo

      
     ! Extrapolate to the final desired omega's
     do iomega=-se%nomega,se%nomega
       se%sigma(iomega,pstate,pspin) = eval_func(coeff, se%omega(iomega) + se%energy0(pstate,pspin) )
     enddo

   enddo
 enddo




contains

function eval_func(coeff_in,zz)
 implicit none
 real(dp),intent(in)    :: coeff_in(nparam*npp)
 complex(dp),intent(in) :: zz
 complex(dp)            :: eval_func
!=====
 integer :: ipp
!=====

 eval_func = 0.0_dp
 do ipp=1,npp
   eval_func = eval_func + coeff_in(5+(ipp-1)*nparam)**2  &
                             / ( zz - ( coeff_in(1+(ipp-1)*nparam)**2 - im * coeff_in(3+(ipp-1)*nparam)**2 ) ) &
                         + coeff_in(6+(ipp-1)*nparam)**2  &
                             / ( zz + ( coeff_in(2+(ipp-1)*nparam)**2 - im * coeff_in(4+(ipp-1)*nparam)**2 ) )
 enddo

end function eval_func


function eval_chi2(coeff_in)
 implicit none
 real(dp),intent(in)    :: coeff_in(nparam*npp)
 real(dp)               :: eval_chi2
!=====
 integer  :: iomegai
 real(dp) :: weight
 real(dp) :: norm
!=====

 eval_chi2 = 0.0_dp
 norm      = 0.0_dp
 do iomegai=-se%nomegai,se%nomegai
   weight = 1.0_dp / ABS(1.0_dp+se%omegai(iomegai))**2
   eval_chi2 = eval_chi2         &
                + ABS( se%sigmai(iomegai,pstate,pspin) &
                      - eval_func(coeff_in, se%omegai(iomegai) + se%energy0(pstate,pspin) ) )**2 &
                      * weight
   norm = norm + weight
                              
 enddo 
 eval_chi2 = eval_chi2 / norm


end function eval_chi2


end subroutine self_energy_fit


!=========================================================================
subroutine self_energy_pade(nstate,energy,se)
 use m_tools,only: pade
 implicit none

 integer,intent(in)                  :: nstate
 real(dp),intent(in)                 :: energy(nstate,nspin)
 type(selfenergy_grid),intent(inout) :: se
!=====
 integer  :: pstate,pspin
 integer  :: iomega
 real(dp) :: sign_eta
!=====

 do pspin=1,nspin
   do pstate=nsemin,nsemax
     do iomega=-se%nomega,se%nomega
       sign_eta = -SIGN( 1.0_dp , REAL(se%omega(iomega),dp) )
       se%sigma(iomega,pstate,pspin) = pade( 2*se%nomegai+1, se%omegai(:) + se%energy0(pstate,pspin), se%sigmai(:,pstate,pspin)  , &
                                              se%omega(iomega) + se%energy0(pstate,pspin) + ieta * sign_eta )
     enddo
   enddo
 enddo


end subroutine self_energy_pade


!=========================================================================
end module m_selfenergy_tools
!=========================================================================
