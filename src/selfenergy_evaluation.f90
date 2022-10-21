!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains
! the driver for the different self-energy methods:
! PT2, GW, evGW, COHSEX, GWGamma, etc.
!
! ymbyun 2019/03/20
! QSGW does not use this file.
!
!=========================================================================
subroutine selfenergy_evaluation(basis,auxil_basis,nstate,occupation,energy,c_matrix, &
                                 exchange_m_vxc_diag)
 use m_definitions
 use m_timing
 use m_warning
 use m_memory
 use m_inputparam
 use m_eri
 use m_eri_calculate
 use m_eri_ao_mo
 use m_dft_grid
 use m_scf,only: en
 use m_hamiltonian
 use m_spectral_function
 use m_selfenergy_tools
 use m_virtual_orbital_space
 implicit none

 type(basis_set),intent(in) :: basis
 type(basis_set),intent(in) :: auxil_basis
 integer,intent(in)         :: nstate
 real(dp),intent(in)        :: occupation(nstate,nspin)
 real(dp),intent(inout)     :: energy(nstate,nspin)
 real(dp),intent(inout)     :: c_matrix(basis%nbf,nstate,nspin)
 real(dp),intent(in)        :: exchange_m_vxc_diag(nstate,nspin)
!=====
 type(selfenergy_grid)   :: se,se2,se3
 character(len=36)       :: selfenergy_tag
 integer                 :: reading_status
 integer                 :: ispin
 integer                 :: nstate_small
 type(spectral_function) :: wpol
 real(dp),allocatable    :: matrix_tmp(:,:,:)
 real(dp),allocatable    :: sigc(:,:)
 real(dp)                :: energy_g(nstate,nspin)
 real(dp)                :: energy_w(nstate,nspin)
 real(dp),allocatable    :: zz(:,:)
 real(dp),allocatable    :: energy_qp_new(:,:),energy_qp_z(:,:)
#ifdef ENABLE_YMBYUN
! ymbyun 2018/07/11
 real(dp),allocatable    :: energy_qp_gra(:,:),energy_qp_A(:,:)
! ymbyun 2019/03/11
! integer                 :: istate
! ymbyun 2019/04/24
 logical                 :: has_mixing_ev
 integer                 :: mixing_ev_file
 real(dp)                :: alpha_local
! ymbyun 2019/08/04
 logical                 :: has_0GW
 integer                 :: file_0GW
 integer                 :: use_0,use_G,use_W
! ymbyun 2019/09/22
 logical                 :: has_se_n_range
 integer                 :: file_se_n_range
 integer                 :: se_n_range
#endif
 integer                 :: iomega
#ifdef COHSEX_DEVEL
 integer,parameter       :: BATCH_SIZE=64
 type(calculation_type)  :: calc_type_tmp
 real(dp),allocatable    :: p_matrix(:,:,:)
 integer                 :: istate
 real(dp)                :: exc
#endif
!=====

 write(stdout,'(/,/,1x,a)') '=================================================='
 write(stdout,'(1x,a)')     'Self-energy evaluation starts here'
 write(stdout,'(1x,a,/)')   '=================================================='

 select case(calc_type%selfenergy_approx)
 case(GW,GnW0,GnWn,G0W0_IOMEGA)
   selfenergy_tag='GW'
 case(PT2)
   selfenergy_tag='PT2'
 case(ONE_RING)
   selfenergy_tag='ONE_RING'
 case(SOX)
   selfenergy_tag='SOX'
 case(G0W0SOX0)
   selfenergy_tag='GWSOX'
 case(G0W0Gamma0)
   selfenergy_tag='GWGamma'
 case(COHSEX,COHSEX_DEVEL,TUNED_COHSEX)
   selfenergy_tag='COHSEX'
 case default
   write(stdout,*) 'selfenergy approx not listed:',calc_type%selfenergy_approx
   call die('selfenergy_evaluation: bug')
 end select

 !
 ! Small imaginary part of the poles in the Green's function
 ! output here
 write(msg,'(es9.2)') AIMAG(ieta)
 call issue_warning('small complex number is '//msg)

#ifdef ENABLE_YMBYUN
 ! ymbyun 2019/09/22
 ! This enables G0W0@scCOHSEX.
 ! NOTE: selfenergy_state_range is changed from protected to public, which breaks the rule (see m_inputparam.f90).
 ! ymbyun 2019/10/31
 ! This enables PT2@scCOHSEX, too.
 if( calc_type%scf_name == 'COHSEX' ) then
   inquire(file='manual_selfenergy_state_range',exist=has_se_n_range)
   if(has_se_n_range) then
     open(newunit=file_se_n_range,file='manual_selfenergy_state_range',form='formatted',status='old')
     read(file_se_n_range,*) se_n_range
     if(se_n_range <= 0) then
       call issue_warning('(ymbyun) File manual_selfenergy_state_range found: se_n_range should be positive')
     else
       call issue_warning('(ymbyun) File manual_selfenergy_state_range found: se_n_range will be used for postscf')
       selfenergy_state_range = se_n_range
     endif
     close(file_se_n_range)
   endif
 endif
#endif

 !
 ! Set the range of states on which to evaluate the self-energy
 call selfenergy_set_state_range(nstate,occupation)

 !
 ! If requested,
 ! prepare an optmized virtual subspace based on 
 ! Frozen Natural Orbitals technique
 if( is_virtual_fno ) then
   !
   ! Be aware that the energies and the c_matrix for virtual orbitals are altered after this point
   ! and until they are restored in destroy_fno
   !
   call virtual_fno(basis,nstate,occupation,energy,c_matrix)
 endif
 !
 ! Or alternatively use the small basis technique
 if( has_small_basis ) then
   if( scalapack_nprow == 1 .AND. scalapack_npcol == 1 ) then
     call setup_virtual_smallbasis(basis,nstate,occupation,nsemax,energy,c_matrix,nstate_small)
   else
     call setup_virtual_smallbasis_sca(basis,nstate,occupation,nsemax,energy,c_matrix,nstate_small)
   endif
   !
   ! Set the range again after the change of the virtual space
   ! to nstate
   call selfenergy_set_state_range(nstate_small,occupation)
 else
   nstate_small = nstate
 endif




 !
 ! Choose which one-electron energies to use in G and in W
 !
 if( calc_type%selfenergy_technique == EVSC ) then
   call read_energy_qp(nstate,energy_g,reading_status)
   if(reading_status/=0) then
     call issue_warning('File energy_qp not found: assuming 1st iteration')
     energy_g(:,:) = energy(:,:)
   endif

   ! 
   ! For GnWn, update both the energy in G and in W
   if( calc_type%selfenergy_approx == GnWn ) then
     energy_w(:,:) = energy_g(:,:)
   else
     energy_w(:,:) = energy(:,:)
   endif

 else
   energy_g(:,:) = energy(:,:)
   energy_w(:,:) = energy(:,:)
 endif

#ifdef ENABLE_YMBYUN
!
! ymbyun 2019/04/04
! 0-COHSEX = one_shot:COHSEX
! evCOHSEX = EVSC:COHSEX (Just for explanation. This does not exist in m_inputparam.f90.)
! scCOHSEX = QS:COHSEX
!
! ymbyun 2019/07/22
! Let's keep it simple: evCOHSEX is viewed as a simplified verison of GnWn.
! evG0W0 (0<Z<1) = eigenvectors (scf) + eigenvalues (energy_qp)
!
! ymbyun 2019/10/01
! evPT2 can't use the Z=1 trick to avoid the multi-solution issue.
! Therefore, evPT2 (0<Z<1) QP energies for states distant from HOMO and LUMO are unreliable.
! Therefore, don't use evPT2 unless you know what you're doing.
!
 if( calc_type%selfenergy_technique == one_shot ) then
   if( calc_type%selfenergy_approx == COHSEX ) then
     call read_energy_qp(nstate,energy_g,reading_status)
     if(reading_status==0) then
       call issue_warning('(ymbyun) File energy_qp found: evCOHSEX will be calculated')
       energy_w(:,:) = energy_g(:,:)  ! energy_g is not used in evCOHSEX.
     endif
   else if( calc_type%selfenergy_approx == GW ) then
     call read_energy_qp(nstate,energy_g,reading_status)
     if(reading_status==0) then
       call issue_warning('(ymbyun) File energy_qp found: evG0W0 with 0<Z<1 will be calculated')
       ! ymbyun 2019/08/04
       ! Update energy_0, energy_G, and energy_W independently.
       use_0 = 0
       use_G = 1
       use_W = 1
       inquire(file='manual_energy_qp_0GW',exist=has_0GW)
       if(has_0GW) then
         open(newunit=file_0GW,file='manual_energy_qp_0GW',form='formatted',status='old')
         read(file_0GW,*) use_0,use_G,use_W
         if(use_W == 0) then
           call issue_warning('(ymbyun) File manual_energy_qp_0GW found: energy_qp will be used only for G')
           energy_w(:,:) = energy(:,:)
         else
           call issue_warning('(ymbyun) File manual_energy_qp_0GW found: energy_qp will be used for both G and W')
           energy_w(:,:) = energy_g(:,:)
         endif
         close(file_0GW)
       else
         call issue_warning('(ymbyun) File manual_energy_qp_0GW not found: energy_qp will be used for both G and W')
         energy_w(:,:) = energy_g(:,:)
       endif
     endif
   endif
!
! ymbyun 2019/10/01
! I'm not sure whether energy == energy_g or energy != energy_g for evPT2 (0<Z<1).
! This part is for energy == energy_g.
! ymbyun 2019/10/26
! MP2 is changed to PT2 (see m_inputparam.f90).
!
! else if( calc_type%selfenergy_technique == EVSC ) then
!   if( calc_type%selfenergy_approx == PT2 ) then
!     call read_energy_qp(nstate,energy_g,reading_status)
!     if(reading_status==0) then
!       call issue_warning('(ymbyun) File energy_qp found: evMP2 will be calculated')
!       energy(:,:)   = energy_g(:,:)
!       energy_w(:,:) = energy_g(:,:)  ! energy_w is not used in evMP2.
!     endif
!   endif
 endif
#endif

 ! ymbyun 2019/03/11
 ! NOTE: istate should be declared first.
! write(stdout,'(/,1x,a,/)') '(ymbyun) Checking eigenvalues...'
!
! do ispin=1,nspin
!   do istate=1,nstate
!     select case(nspin)
!     case(1)
!       write(stdout,'(1x,i3,3(1x,f12.5))') istate,energy(istate,:),energy_g(istate,:),energy_w(istate,:)
!     case(2)
!       write(stdout,'(1x,i3,3(2(1x,f12.5)))') istate,energy(istate,:),energy_g(istate,:),energy_w(istate,:)
!     end select
!   enddo
! enddo

#ifdef ENABLE_YMBYUN
 ! ymbyun 2019/08/04
 !   G0W0: energy == energy_g
 ! evG0W0: energy != energy_g
 if( calc_type%selfenergy_technique == one_shot .AND. calc_type%selfenergy_approx == GW ) then
   call init_selfenergy_grid(calc_type%selfenergy_technique,nstate,energy  ,se)
 else
   call init_selfenergy_grid(calc_type%selfenergy_technique,nstate,energy_g,se)
 endif
#else
 call init_selfenergy_grid(calc_type%selfenergy_technique,nstate,energy_g,se)
#endif


 !
 ! selfenergy = GW or COHSEX
 !
 if(     calc_type%selfenergy_approx == GV .OR. calc_type%selfenergy_approx == GSIGMA .OR.  calc_type%selfenergy_approx == LW &
    .OR. calc_type%selfenergy_approx == LW2 &
    .OR. calc_type%selfenergy_approx == G0W0_IOMEGA &
    .OR. calc_type%selfenergy_approx == GW   .OR. calc_type%selfenergy_approx == COHSEX   &
    .OR. calc_type%selfenergy_approx == GnW0 .OR. calc_type%selfenergy_approx == GnWn   ) then


   call init_spectral_function(nstate_small,occupation,nomega_imag,wpol)

   ! Try to read a spectral function file in order to skip the polarizability calculation
   ! Skip the reading if GnWn (=evGW) is requested
   !
   ! ymbyun 2019/04/04
   ! Not only GnWn, but also evCOHSEX.
   ! TODO: one_shot:COHSEX should be distinguished from QS:COHSEX.
#ifdef ENABLE_YMBYUN
   if( calc_type%selfenergy_approx /= GnWn .AND. calc_type%selfenergy_approx /= COHSEX ) then
     call read_spectral_function(wpol,reading_status)
   else
     write(stdout,'(/,1x,a)') 'For GnWn and evCOHSEX calculations, never try to read file SCREENED_COULOMB'
     reading_status = 1
   endif
#else
   if( calc_type%selfenergy_approx /= GnWn ) then
     call read_spectral_function(wpol,reading_status)
   else
     write(stdout,'(/,1x,a)') 'For GnWn calculations, never try to read file SCREENED_COULOMB'
     reading_status = 1
   endif
#endif
   ! If reading has failed, then do the calculation
   if( reading_status /= 0 ) then
     if( calc_type%selfenergy_technique /= imaginary_axis ) then
       call polarizability(basis,nstate,occupation,energy_w,c_matrix,en%rpa,wpol)
     else
       call polarizability_grid_scalapack(basis,nstate,occupation,energy_w,c_matrix,en%rpa,wpol)
     endif
   endif

   en%tot = en%tot + en%rpa
   if( calc_type%is_dft ) en%tot = en%tot - en%xc - en%exx_hyb + en%exx 
   if( ABS(en%rpa) > 1.e-6_dp) then
     write(stdout,'(/,a,f19.10)') ' RPA Total energy (Ha): ',en%tot
   endif


#ifdef HAVE_SCALAPACK
   ! The SCALAPACK implementation only works for plain vanilla GW
   ! TODO: extend it to COHSEX
   if( has_auxil_basis &
     .AND. (calc_type%selfenergy_approx == GW .OR. calc_type%selfenergy_approx == GnW0  &
       .OR. calc_type%selfenergy_approx == GnWn .OR. calc_type%selfenergy_approx == G0W0_IOMEGA) ) then
     if( calc_type%selfenergy_technique /= imaginary_axis ) then
       call gw_selfenergy_scalapack(calc_type%selfenergy_approx,nstate,basis,occupation,energy_g,c_matrix,wpol,se)
     else
       call gw_selfenergy_imag_scalapack(basis,nstate,occupation,energy_g,c_matrix,wpol,se)
       call self_energy_pade(nstate,energy_g,se)
     endif
   else
     call gw_selfenergy(calc_type%selfenergy_approx,nstate,basis,occupation,energy_g,c_matrix,wpol,se,en%gw)
   endif
#else
   if( calc_type%selfenergy_technique /= imaginary_axis ) then
     call gw_selfenergy(calc_type%selfenergy_approx,nstate,basis,occupation,energy_g,c_matrix,wpol,se,en%gw)
   else
     call gw_selfenergy_imag_scalapack(basis,nstate,occupation,energy_g,c_matrix,wpol,se)
     call self_energy_pade(nstate,energy_g,se)
   endif
#endif

   

   if( ABS(en%gw) > 1.0e-5_dp ) then
     write(stdout,'(/,a,f19.10)') ' Galitskii-Migdal Total energy (Ha): ',en%tot - en%rpa + en%gw
   endif

   call destroy_spectral_function(wpol)

   if( has_small_basis ) then
     !
     ! Output the G0W0 results in the small basis first
     allocate(energy_qp_z(nstate,nspin))
     allocate(energy_qp_new(nstate,nspin))
     allocate(zz(nsemin:nsemax,nspin))
     call find_qp_energy_linearization(se,nstate,exchange_m_vxc_diag,energy,energy_qp_z,zz)
     call find_qp_energy_graphical(se,nstate,exchange_m_vxc_diag,energy,energy_qp_new)
     call output_qp_energy('GW small basis',nstate,energy,exchange_m_vxc_diag,1,se,energy_qp_z,energy_qp_new,zz)
     deallocate(zz)
     deallocate(energy_qp_z)
     call output_new_homolumo('GW small basis',nstate,occupation,energy_qp_new,nsemin,nsemax)
     deallocate(energy_qp_new)

     call init_selfenergy_grid(static_selfenergy,nstate,energy,se2)
     call init_selfenergy_grid(static_selfenergy,nstate,energy,se3)
  
     ! Sigma^2 = Sigma^{1-ring}_small
     call onering_selfenergy(ONE_RING,nstate_small,basis,occupation(1:nstate_small,:), &
                             energy_g(1:nstate_small,:),c_matrix(:,1:nstate_small,:),se2,en%mp2)

     ! Reset wavefunctions, eigenvalues and number of virtual orbitals in G
     call destroy_fno(basis,nstate,energy,c_matrix)
     energy_g(:,:) = energy(:,:)
     call selfenergy_set_state_range(nstate,occupation)

     ! Sigma^3 = Sigma^{1-ring}_big
     call onering_selfenergy(ONE_RING,nstate,basis,occupation,energy_g,c_matrix,se3,en%mp2)

     if( print_sigma_ ) then
       call write_selfenergy_omega('selfenergy_GW_small'   ,nstate,exchange_m_vxc_diag,energy_g,se)
       call write_selfenergy_omega('selfenergy_1ring_big'  ,nstate,exchange_m_vxc_diag,energy_g,se3)
       call write_selfenergy_omega('selfenergy_1ring_small',nstate,exchange_m_vxc_diag,energy_g,se2)
     endif

     !
     ! Extrapolated Sigma(omega) = Sigma^{GW}_small(omega) + Sigma^{1-ring}_big(0) - Sigma^{1-ring}_small(0)
     do iomega=-se%nomega,se%nomega
       se%sigma(iomega,:,:) = se%sigma(iomega,:,:) + se3%sigma(0,:,:) - se2%sigma(0,:,:)
     enddo

     call destroy_selfenergy_grid(se2)
     call destroy_selfenergy_grid(se3)

   endif

 endif

 !
 ! GWGamma
 !
 if( calc_type%selfenergy_approx == G0W0GAMMA0 .OR. calc_type%selfenergy_approx == G0W0SOX0 ) then
   call init_spectral_function(nstate,occupation,0,wpol)
   call read_spectral_function(wpol,reading_status)
   ! If reading has failed, then do the calculation
   if( reading_status /= 0 ) then
     call polarizability(basis,nstate,occupation,energy_w,c_matrix,en%rpa,wpol)
   endif

   call gw_selfenergy(GW,nstate,basis,occupation,energy_g,c_matrix,wpol,se,en%gw)

   !
   ! Output the G0W0 results first
   allocate(energy_qp_z(nstate,nspin))
   allocate(energy_qp_new(nstate,nspin))
   allocate(zz(nsemin:nsemax,nspin))
   call find_qp_energy_linearization(se,nstate,exchange_m_vxc_diag,energy,energy_qp_z,zz)
   call find_qp_energy_graphical(se,nstate,exchange_m_vxc_diag,energy,energy_qp_new)
   call output_qp_energy('GW',nstate,energy,exchange_m_vxc_diag,1,se,energy_qp_z,energy_qp_new,zz)
   deallocate(zz)
   deallocate(energy_qp_z)
   call output_new_homolumo('GW',nstate,occupation,energy_qp_new,nsemin,nsemax)
   deallocate(energy_qp_new)


   call gwgamma_selfenergy(nstate,basis,occupation,energy_g,c_matrix,wpol,se)
   call destroy_spectral_function(wpol)
 endif


 !
 ! Selfenergy = PT2
 ! 
 if(   calc_type%selfenergy_approx == PT2       &
  .OR. calc_type%selfenergy_approx == ONE_RING  &
  .OR. calc_type%selfenergy_approx == SOX ) then

   call pt2_selfenergy(calc_type%selfenergy_approx,nstate,basis,occupation,energy_g,c_matrix,se,en%mp2)

   if( ABS( en%mp2 ) > 1.0e-8 ) then
     write(stdout,'(a,2x,f19.10)') ' MP2 Energy       (Ha):',en%mp2
     write(stdout,*)
     en%tot = en%nuc_nuc + en%kin + en%nuc + en%hart + en%exx + en%mp2

     write(stdout,'(a,2x,f19.10)') ' MP2 Total Energy (Ha):',en%tot
     write(stdout,'(a,2x,f19.10)') ' SE+MP2  Total En (Ha):',en%tot+en%se
     write(stdout,*)
   endif

 endif


 !
 ! EXPERIMENTAL COHSEX implementation
 ! final evaluation for perturbative COHSEX
 !
 if( calc_type%selfenergy_approx == COHSEX_DEVEL .OR. calc_type%selfenergy_approx == TUNED_COHSEX ) then

   if( .NOT. has_auxil_basis ) call die('cohsex needs an auxiliary basis')
   call init_spectral_function(nstate,occupation,1,wpol)
   call calculate_eri_3center_eigen(basis%nbf,nstate,c_matrix,ncore_W+1,nhomo_W,nlumo_W,nvirtual_W-1)
   !
   ! Calculate v^{1/2} \chi v^{1/2}
   call static_polarizability(nstate,occupation,energy_w,wpol)

   call destroy_eri_3center_eigen()

   !
   allocate(matrix_tmp(basis%nbf,basis%nbf,nspin))
   allocate(sigc(nstate,nspin))

#ifdef COHSEX_DEVEL
   ! Calculate the DFT potential part
   if( ABS( delta_cohsex ) > 1.0e-6_dp ) then

     allocate(p_matrix(basis%nbf,basis%nbf,nspin))
     call init_dft_grid(basis,grid_level,.TRUE.,.FALSE.,BATCH_SIZE)
     call setup_density_matrix(basis%nbf,nstate,c_matrix,occupation,p_matrix)

     ! Override the DFT XC correlation settings
     calc_type_tmp = calc_type
     call init_dft_type('HJSx',calc_type_tmp)
#ifdef HAVE_LIBXC
     call xc_f90_gga_x_hjs_set_par(calc_type_tmp%xc_func(1),1.0_dp/rcut_mbpt)
#endif
     call dft_exc_vxc_batch(BATCH_SIZE,basis,nstate,occupation,c_matrix,matrix_tmp,exc)
 
     write(stdout,*) '===== SigX SR ======'
     do ispin=1,nspin
       do istate=1,nstate
         sigc(istate,ispin) = DOT_PRODUCT( c_matrix(:,istate,ispin) , &
                                   MATMUL( matrix_tmp(:,:,ispin) , c_matrix(:,istate,ispin ) ) )
         write(stdout,*) istate,ispin,sigc(istate,ispin) * Ha_eV
       enddo
     enddo
     write(stdout,*) '===================='
     sigc(istate,ispin) = sigc(istate,ispin) * delta_cohsex 

     deallocate(p_matrix)
     call destroy_dft_grid()

   else

     sigc(:,:) = 0.0_dp

   endif

#endif

   call cohsex_selfenergy(nstate,basis,occupation,energy_g,exchange_m_vxc_diag, & 
                          c_matrix,wpol,matrix_tmp,sigc,en%gw)


   !
   ! A section under development for the range-separated RPA
   if( calc_type%is_lr_mbpt ) then

     ! 2-center integrals
     call calculate_eri_2center_scalapack(auxil_basis,rcut_mbpt)
     ! 3-center integrals
     call calculate_eri_3center_scalapack(basis,auxil_basis,rcut_mbpt)

     call cohsex_selfenergy_lr(nstate,basis,occupation,energy_g,exchange_m_vxc_diag, &
                               c_matrix,wpol,matrix_tmp,sigc,en%gw)
   endif

   deallocate(matrix_tmp)
   deallocate(sigc)

 endif ! COHSEX
 !
 ! end of EXPERIMENTAL COHSEX implementation
 !


 !
 ! Output the quasiparticle energies, the self-energy etc.
 !

#ifdef ENABLE_YMBYUN
! ymbyun 2018/07/11
! Moved up (see below)
 allocate(energy_qp_new(nstate,nspin))
#endif

 if( print_sigma_ ) then
! ymbyun 2018/07/11
! The QP energy Eqp is obtained from the spectral function A(w) while A(w) is being written to a file.
#ifdef ENABLE_YMBYUN
   ! ymbyun 2019/08/04
   !   G0W0: energy == energy_g
   ! evG0W0: energy != energy_g
   if( calc_type%selfenergy_technique == one_shot .AND. calc_type%selfenergy_approx == GW ) then
     call write_selfenergy_omega('selfenergy_'//TRIM(selfenergy_tag),nstate,exchange_m_vxc_diag,energy  ,se, energy_qp_new)
   else
     ! ymbyun 2019/09/26
     ! At first, I warn users without disabling this part. It it works, I need to do nothing. If not, I have to disable this part.
     ! ymbyun 2019/10/25
     ! Not only GnW0 and GnWn, but also COHSEX
     if( calc_type%selfenergy_approx == GnW0 .OR. calc_type%selfenergy_approx == GnWn .OR. calc_type%selfenergy_approx == COHSEX ) then
       call issue_warning('(ymbyun) Please set print_sigma to ''no'' for GnW0, GnWn, and COHSEX')
     endif

     ! ymbyun 2019/09/26
     ! PT2 uses real numbers, not complex numbers.
     ! ymbyun 2019/09/27
     ! Both one-shot PT2 and evPT2 (see m_inputparam.f90)
     if( calc_type%selfenergy_approx == PT2 ) then
       call write_selfenergy_omega('selfenergy_'//TRIM(selfenergy_tag),nstate,exchange_m_vxc_diag,energy_g,se)
     else
       call write_selfenergy_omega('selfenergy_'//TRIM(selfenergy_tag),nstate,exchange_m_vxc_diag,energy_g,se, energy_qp_new)
     endif
   endif
#else
   call write_selfenergy_omega('selfenergy_'//TRIM(selfenergy_tag),nstate,exchange_m_vxc_diag,energy_g,se)
#endif
 endif

#ifndef ENABLE_YMBYUN
! ymbyun 2018/07/11
! Moved up (see above)
 allocate(energy_qp_new(nstate,nspin))
#endif

 select case(calc_type%selfenergy_approx)
! ymbyun 2019/09/21
! Unlike GnW0, GnWn, and evCOHSEX, evPT2 can't use the Z=1 trick to avoid the multi-soluton issue.
! Don't use evPT2 unless you know what you're doing (see my comment above).
 case(GW,PT2,ONE_RING,SOX,G0W0Gamma0,G0W0SOX0,G0W0_IOMEGA)
! ymbyun 2018/07/11
! The QP energy Eqp is obtained by three different methods:
!   1. Linearization of the non-linear QP equation
!   2. Graph
!   3. Spectral function A(w)
#ifdef ENABLE_YMBYUN
   allocate(energy_qp_gra(nstate,nspin))
   allocate(energy_qp_z(nstate,nspin))
   allocate(zz(nsemin:nsemax,nspin))

   call find_qp_energy_linearization(se,nstate,exchange_m_vxc_diag,energy,energy_qp_z,zz)
   call find_qp_energy_graphical(se,nstate,exchange_m_vxc_diag,energy,energy_qp_gra)

   if( .NOT. print_sigma_ ) then
     energy_qp_new(:,:) = energy_qp_gra(:,:)
   endif

   ! ymbyun 2019/09/26
   ! PT2 uses real numbers, not complex numbers.
   ! ymbyun 2019/09/27
   ! Both one-shot PT2 and evPT2 (see m_inputparam.f90)
   if( calc_type%selfenergy_approx == PT2 ) then
     energy_qp_new(:,:) = energy_qp_gra(:,:)
   endif

   call output_qp_energy(TRIM(selfenergy_tag),nstate,energy,exchange_m_vxc_diag,1,se,energy_qp_z,energy_qp_gra,zz, energy_qp_new)

   deallocate(energy_qp_gra)
   deallocate(energy_qp_z)
   deallocate(zz)
#else
   allocate(energy_qp_z(nstate,nspin))
   allocate(zz(nsemin:nsemax,nspin))

   call find_qp_energy_linearization(se,nstate,exchange_m_vxc_diag,energy,energy_qp_z,zz)
   call find_qp_energy_graphical(se,nstate,exchange_m_vxc_diag,energy,energy_qp_new)
  
   call output_qp_energy(TRIM(selfenergy_tag),nstate,energy,exchange_m_vxc_diag,1,se,energy_qp_z,energy_qp_new,zz)

   deallocate(energy_qp_z)
   deallocate(zz)
#endif
 case(GnWn,GnW0,GV,COHSEX,COHSEX_DEVEL,TUNED_COHSEX)
   call find_qp_energy_linearization(se,nstate,exchange_m_vxc_diag,energy,energy_qp_new)

! ymbyun 2019/04/24
! TODO: one_shot:COHSEX should be distinguished from QS:COHSEX.
#ifdef ENABLE_YMBYUN
   if( calc_type%selfenergy_approx == GnW0 .OR. calc_type%selfenergy_approx == GnWn .OR. calc_type%selfenergy_approx == COHSEX ) then
     inquire(file='manual_mixing_ev',exist=has_mixing_ev)
     if(has_mixing_ev) then
       open(newunit=mixing_ev_file,file='manual_mixing_ev',status='old')
       read(mixing_ev_file,*) alpha_local
       close(mixing_ev_file)
       write(msg,'(a,f12.6)') '(ymbyun) mixing evGW or evCOHSEX eigenvalues with alpha = ',alpha_local
       call issue_warning(msg)

       ! ymbyun 2019/07/22
       ! Let's keep it simple: evCOHSEX is viewed as a simplified verison of GnWn.
!       if( calc_type%selfenergy_approx == COHSEX ) then
!         energy_qp_new(:,:) = alpha_local * energy_qp_new(:,:) + (1.0_dp - alpha_local) * energy_w(:,:)
!       else
         energy_qp_new(:,:) = alpha_local * energy_qp_new(:,:) + (1.0_dp - alpha_local) * energy_g(:,:)
!       endif
     endif
   endif
#endif

   call output_qp_energy(TRIM(selfenergy_tag),nstate,energy,exchange_m_vxc_diag,1,se,energy_qp_new)
 end select

 !
 ! Write the QP energies on disk: ENERGY_QP file
 ! 
 call write_energy_qp(nstate,energy_qp_new)

 !
 ! Output the new HOMO and LUMO energies
 !
 call output_new_homolumo(TRIM(selfenergy_tag),nstate,occupation,energy_qp_new,nsemin,nsemax)



 deallocate(energy_qp_new)



 !
 ! Deallocations
 !
 call destroy_selfenergy_grid(se)


end subroutine selfenergy_evaluation


!=========================================================================
