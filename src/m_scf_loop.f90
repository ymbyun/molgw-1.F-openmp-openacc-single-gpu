!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the main SCF loop for Hartree-Fock or generalized Kohn-Sham
!
!=========================================================================
module m_scf_loop
 use m_definitions
 use m_timing
 use m_warning
 use m_memory
 use m_inputparam

 integer,parameter,private :: BATCH_SIZE = 64

contains


!=========================================================================
subroutine scf_loop(is_restart,& 
                    basis,auxil_basis,&
                    nstate,m_ham,n_ham,m_c,n_c,&
                    s_matrix_sqrt_inv,s_matrix,&
                    hamiltonian_kinetic,hamiltonian_nucleus,&
                    occupation, &
                    energy, &
                    hamiltonian_fock,&
                    c_matrix)
 use m_tools
 use m_atoms
 use m_basis_set
 use m_scf
 use m_eri
 use m_eri_calculate
 use m_eri_ao_mo
 use m_dft_grid
 use m_spectral_function
 use m_hamiltonian
 use m_hamiltonian_sca
 use m_hamiltonian_buffer
 use m_selfenergy_tools
 implicit none

!=====
 logical,intent(in)                 :: is_restart
 type(basis_set),intent(in)         :: basis
 type(basis_set),intent(in)         :: auxil_basis
 integer,intent(in)                 :: nstate,m_ham,n_ham,m_c,n_c
 real(dp),intent(in)                :: s_matrix_sqrt_inv(m_c,n_c)
 real(dp),intent(in)                :: s_matrix(m_ham,n_ham)
 real(dp),intent(in)                :: hamiltonian_kinetic(m_ham,n_ham)
 real(dp),intent(in)                :: hamiltonian_nucleus(m_ham,n_ham)
 real(dp),intent(inout)             :: occupation(nstate,nspin)
 real(dp),intent(out)               :: energy(nstate,nspin)
 real(dp),allocatable,intent(inout) :: hamiltonian_fock(:,:,:)
 real(dp),allocatable,intent(inout) :: c_matrix(:,:,:)
!=====
 type(spectral_function) :: wpol
 logical                 :: is_converged,stopfile_found
 integer                 :: ispin,iscf,istate
 real(dp)                :: energy_tmp
 real(dp),allocatable    :: p_matrix(:,:,:)
 real(dp),allocatable    :: hamiltonian(:,:,:)
 real(dp),allocatable    :: hamiltonian_hartree(:,:)
 real(dp),allocatable    :: hamiltonian_exx(:,:,:)
 real(dp),allocatable    :: hamiltonian_xc(:,:,:)
 real(dp),allocatable    :: matrix_tmp(:,:,:)
 real(dp),allocatable    :: energy_exx(:,:)
 real(dp),allocatable    :: c_matrix_exx(:,:,:)
!=============================


 call start_clock(timing_scf)

 ! Old Fock operator will be updated
 ! Get rid of it!
 call clean_deallocate('Fock operator F',hamiltonian_fock) ! Never distributed

 !
 ! Initialize the SCF mixing procedure
 call init_scf(m_ham,n_ham,m_c,n_c,basis%nbf,nstate)

 !
 ! Allocate the main arrays
 call clean_allocate('Total Hamiltonian H',hamiltonian,m_ham,n_ham,nspin)
 call clean_allocate('Hartree potential Vh',hamiltonian_hartree,m_ham,n_ham)
 call clean_allocate('Exchange operator Sigx',hamiltonian_exx,m_ham,n_ham,nspin)
 call clean_allocate('XC operator Vxc',hamiltonian_xc,m_ham,n_ham,nspin)
 call clean_allocate('Density matrix P',p_matrix,m_ham,n_ham,nspin)


 if( calc_type%is_dft ) then
   !
   ! Setup the grids for the quadrature of DFT potential/energy
   call init_dft_grid(basis,grid_level,dft_xc_needs_gradient,.TRUE.,BATCH_SIZE)
   ! The following is coded but not used... yet!
!   call setup_bf_radius(basis)
 endif

 !
 ! Setup the density matrix: p_matrix
 if( parallel_ham ) then
   call setup_density_matrix_sca(basis%nbf,nstate,m_c,n_c,c_matrix,occupation,m_ham,n_ham,p_matrix)
 else
   call setup_density_matrix(basis%nbf,nstate,c_matrix,occupation,p_matrix)
 endif


 !
 ! start the big scf loop
 !
 do iscf=1,nscf
   write(stdout,'(/,a)') '-------------------------------------------'
   write(stdout,'(a,1x,i4,/)') ' *** SCF cycle No:',iscf


   if( cntxt_ham > 0 ) then
     en%kin  = SUM( hamiltonian_kinetic(:,:) * SUM(p_matrix(:,:,:),DIM=3) )
     en%nuc  = SUM( hamiltonian_nucleus(:,:) * SUM(p_matrix(:,:,:),DIM=3) )
   else
     en%kin  = 0.0_dp
     en%nuc  = 0.0_dp
   endif
   call xsum_trans(en%kin)
   call xsum_trans(en%nuc)

   !
   ! Setup kinetic and nucleus contributions (that are independent of the
   ! density matrix and therefore of spin channel)
   !
   hamiltonian(:,:,1) = hamiltonian_kinetic(:,:) + hamiltonian_nucleus(:,:)
   if(nspin==2) then
     hamiltonian(:,:,nspin) = hamiltonian_kinetic(:,:) + hamiltonian_nucleus(:,:) 
   endif

   !
   ! Hartree contribution to the Hamiltonian
   !
   if( .NOT. has_auxil_basis ) then
     call setup_hartree(print_matrix_,basis%nbf,p_matrix,hamiltonian_hartree,en%hart)
   else
     if( parallel_ham ) then
       if( parallel_buffer ) then
         call setup_hartree_ri_buffer_sca(print_matrix_,basis%nbf,m_ham,n_ham,p_matrix,hamiltonian_hartree,en%hart)
       else
         call setup_hartree_ri_sca(print_matrix_,basis%nbf,m_ham,n_ham,p_matrix,hamiltonian_hartree,en%hart)
       endif
     else
       call setup_hartree_ri(print_matrix_,basis%nbf,p_matrix,hamiltonian_hartree,en%hart)
     endif
   endif
   ! calc_type%is_core is an inefficient way to get the Kinetic+Nucleus Hamiltonian
   if( calc_type%is_core ) then
     hamiltonian_hartree(:,:) = 0.0_dp
     en%hart = 0.0_dp
   endif
   do ispin=1,nspin
     hamiltonian(:,:,ispin) = hamiltonian(:,:,ispin) + hamiltonian_hartree(:,:)
   enddo


   !
   !  XC part of the Hamiltonian
   !
   hamiltonian_xc(:,:,:) = 0.0_dp
   en%exx_hyb = 0.0_dp

   !
   ! DFT XC potential is added here
   ! hamiltonian_xc is used as a temporary matrix
   if( calc_type%is_dft ) then

     if( parallel_ham ) then
       if( parallel_buffer ) then
         call dft_exc_vxc_buffer_sca(basis,nstate,m_c,n_c,m_ham,n_ham,occupation,c_matrix,p_matrix,hamiltonian_xc,en%xc)
       else
         call issue_warning('Exc calculation with SCALAPACK is not coded yet. Just skip it')
         hamiltonian_xc(:,:,:) = 0.0_dp
         en%xc = 0.0_dp
       endif
     else
       call dft_exc_vxc_batch(BATCH_SIZE,basis,nstate,occupation,c_matrix,hamiltonian_xc,en%xc)
     endif

   endif

   !
   ! LR Exchange contribution to the Hamiltonian
   ! Use hamiltonian_exx as a temporary matrix (no need to save it for later use)
   if(calc_type%need_exchange_lr) then

     if( .NOT. has_auxil_basis ) then
       call setup_exchange_longrange(print_matrix_,basis%nbf,p_matrix,hamiltonian_exx,energy_tmp)
     else
       if( parallel_ham ) then
         if( parallel_buffer ) then
           call setup_exchange_longrange_ri_buffer_sca(basis%nbf,nstate,m_c,n_c,m_ham,n_ham,occupation,c_matrix,p_matrix,hamiltonian_exx,energy_tmp)
         else
           call die('Range-separated functionals not implemented with full SCALAPACK yet')
         endif
       else
         call setup_exchange_longrange_ri(basis%nbf,nstate,occupation,c_matrix,p_matrix,hamiltonian_exx,energy_tmp)
       endif
     endif
     ! Rescale with alpha_hybrid_lr for range-separated hybrid functionals
     en%exx_hyb = en%exx_hyb + alpha_hybrid_lr * energy_tmp
     hamiltonian_xc(:,:,:) = hamiltonian_xc(:,:,:) + hamiltonian_exx(:,:,:) * alpha_hybrid_lr
   endif

   !
   ! Exchange contribution to the Hamiltonian
   if( calc_type%need_exchange ) then

     if( .NOT. has_auxil_basis ) then
       call setup_exchange(print_matrix_,basis%nbf,p_matrix,hamiltonian_exx,en%exx)
     else
       if( parallel_ham ) then
         if( parallel_buffer ) then
#ifndef SCASCA
           call setup_exchange_ri_buffer_sca(basis%nbf,nstate,m_c,n_c,m_ham,n_ham,occupation,c_matrix,p_matrix,hamiltonian_exx,en%exx)
#else
           call issue_warning('FBFB devel SCASCA')
           call setup_exchange_ri_sca(basis%nbf,nstate,m_c,n_c,m_ham,n_ham,occupation,c_matrix,p_matrix,hamiltonian_exx,en%exx)
#endif
         else
           call die('Exchange with fully distributed hamiltonian: case not implemented yet')
         endif
       else
         call setup_exchange_ri(basis%nbf,nstate,occupation,c_matrix,p_matrix,hamiltonian_exx,en%exx)
       endif
     endif
     ! Rescale with alpha_hybrid for hybrid functionals
     en%exx_hyb = en%exx_hyb + alpha_hybrid * en%exx
     hamiltonian_xc(:,:,:) = hamiltonian_xc(:,:,:) + hamiltonian_exx(:,:,:) * alpha_hybrid
   endif


   !
   ! QSGW or COHSEX self energy
   if( ( calc_type%selfenergy_approx == GW .OR. calc_type%selfenergy_approx == COHSEX ) &
        .AND. calc_type%selfenergy_technique == QS  &
        .AND. ( iscf > 5 .OR. is_restart ) ) then

     if( parallel_ham ) call die('QSGW not implemented with parallel_ham')

     call init_spectral_function(nstate,occupation,0,wpol)
     call polarizability(basis,nstate,occupation,energy,c_matrix,en%rpa,wpol)

     if( ABS(en%rpa) > 1.e-6_dp) then
       en%tot = en%tot + en%rpa
       write(stdout,'(/,a,f19.10)') ' RPA Total energy (Ha): ',en%tot
     endif

     !
     ! Set the range of states on which to evaluate the self-energy
     call selfenergy_set_state_range(nstate,occupation)

     allocate(matrix_tmp(m_ham,n_ham,nspin))
     call gw_selfenergy_qs(nstate,basis,occupation,energy,c_matrix,s_matrix,wpol,matrix_tmp)

     call dump_out_matrix(print_matrix_,'=== Self-energy ===',basis%nbf,nspin,matrix_tmp)
     call destroy_spectral_function(wpol)
     hamiltonian(:,:,:) = hamiltonian(:,:,:) + matrix_tmp(:,:,:)
     deallocate(matrix_tmp)

   endif

   !
   ! QSPT2
   if( calc_type%selfenergy_approx == PT2 .AND. calc_type%selfenergy_technique == QS .AND. ( iscf > 5 .OR. is_restart ) ) then

     if( parallel_ham ) call die('QSPT2 not implemented with parallel_ham')

     !
     ! Set the range of states on which to evaluate the self-energy
     call selfenergy_set_state_range(nstate,occupation)

     allocate(matrix_tmp(m_ham,n_ham,nspin))
     call pt2_selfenergy_qs(nstate,basis,occupation,energy,c_matrix,s_matrix,matrix_tmp,en%mp2)

     write(stdout,'(a,2x,f19.10)') ' MP2 Energy       (Ha):',en%mp2
     write(stdout,*) 
     en%tot = en%tot + en%mp2
     write(stdout,'(a,2x,f19.10)') ' MP2 Total Energy (Ha):',en%tot

     call dump_out_matrix(print_matrix_,'=== Self-energy ===',basis%nbf,nspin,matrix_tmp)

     hamiltonian(:,:,:) = hamiltonian(:,:,:) + matrix_tmp(:,:,:)
     deallocate(matrix_tmp)

   endif

!FBFB electrostatic potential
!   write(stdout,*) '== FBFB ==',(hamiltonian_nucleus(basis%nbf,basis%nbf)+hamiltonian_hartree(basis%nbf,basis%nbf))*Ha_eV

   !
   ! Add the XC part of the hamiltonian to the total hamiltonian
   hamiltonian(:,:,:) = hamiltonian(:,:,:) + hamiltonian_xc(:,:,:)
   
   ! All the components of the energy have been calculated at this stage
   ! Sum up to get the total energy
   en%tot = en%nuc_nuc + en%kin + en%nuc + en%hart + en%exx_hyb + en%xc

   !
   ! If requested, the level shifting procedure is triggered: 
   ! All the unoccupied states are penalized with an energy =  level_shifting_energy
   if( level_shifting_energy > 1.0e-6_dp ) then
     if( parallel_ham ) call die('level_shifting: not implemented with parallel_ham')
     call level_shifting_up(basis%nbf,nstate,s_matrix,c_matrix,occupation,level_shifting_energy,hamiltonian)
   endif


   ! DIIS or simple mixing on the hamiltonian
   call hamiltonian_prediction(s_matrix,s_matrix_sqrt_inv,p_matrix,hamiltonian)

  
   !
   ! Diagonalize the Hamiltonian H
   ! Generalized eigenvalue problem with overlap matrix S
   ! H \varphi = E S \varphi
   ! save the old eigenvalues
   if( parallel_ham ) then
     call diagonalize_hamiltonian_sca(1,nspin,desc_ham,hamiltonian,desc_c,s_matrix_sqrt_inv,energy,c_matrix)
   else
     ! This subroutine works with or without scalapack
     call diagonalize_hamiltonian_scalapack(nspin,basis%nbf,nstate,hamiltonian,s_matrix_sqrt_inv,energy,c_matrix)
   endif

   !
   ! When level_shifting is used, the unoccupied state energies have to be brought
   ! back to their original value,
   ! So that the "physical" energies are written down
   if( level_shifting_energy > 1.0e-6_dp ) then
     if( parallel_ham ) call die('level_shifting: not implemented with parallel_ham')
     call level_shifting_down(basis%nbf,nstate,s_matrix,c_matrix,occupation,level_shifting_energy,energy,hamiltonian)
   endif
  
   call dump_out_energy('=== Energies ===',nstate,nspin,occupation,energy)

   call output_new_homolumo('gKS',nstate,occupation,energy,1,nstate)


   !
   ! Output the total energy and its components
   write(stdout,*)
   write(stdout,'(a25,1x,f19.10)') 'Nucleus-Nucleus (Ha):',en%nuc_nuc
   write(stdout,'(a25,1x,f19.10)') 'Kinetic Energy  (Ha):',en%kin
   write(stdout,'(a25,1x,f19.10)') 'Nucleus Energy  (Ha):',en%nuc
   write(stdout,'(a25,1x,f19.10)') 'Hartree Energy  (Ha):',en%hart
   if(calc_type%need_exchange) then
     write(stdout,'(a25,1x,f19.10)') 'Exchange Energy (Ha):',en%exx_hyb
   endif
   if( calc_type%is_dft ) then
     write(stdout,'(a25,1x,f19.10)') 'XC Energy       (Ha):',en%xc
   endif
   write(stdout,'(/,a25,1x,f19.10,/)') 'Total Energy    (Ha):',en%tot


   ! If fractional occupancies are allowed, then recalculate the occupations
   if( temperature > 1.0e-8_dp ) then
     call set_occupation(nstate,temperature,electrons,magnetization,energy,occupation)
   endif

   !
   ! Setup the new density matrix: p_matrix
   ! Save the old one for the convergence criterium
   if( parallel_ham ) then
     call setup_density_matrix_sca(basis%nbf,nstate,m_c,n_c,c_matrix,occupation,m_ham,n_ham,p_matrix)
   else
     call setup_density_matrix(basis%nbf,nstate,c_matrix,occupation,p_matrix)
   endif

   is_converged = check_converged(p_matrix)
   inquire(file='STOP',exist=stopfile_found)

   if( is_converged .OR. stopfile_found ) exit

   !
   ! Write down a "small" RESTART file at each step
   ! Skip writing when parallel_ham since all the information is not available on the master proc.
   if( print_restart_ .AND. .NOT. parallel_ham ) then
     call clean_allocate('Fock operator F',hamiltonian_fock,basis%nbf,basis%nbf,nspin)
     call get_fock_operator(hamiltonian,hamiltonian_xc,hamiltonian_exx,hamiltonian_fock)
     call write_restart(SMALL_RESTART,basis,nstate,occupation,c_matrix,energy,hamiltonian_fock)
     call clean_deallocate('Fock operator F',hamiltonian_fock)
   endif

   
 !
 ! end of the big SCF loop
 enddo


 write(stdout,'(/,1x,a)') '=================================================='
 write(stdout,'(1x,a)') 'The SCF loop ends here'
 write(stdout,'(1x,a)') '=================================================='

 !
 ! Cleanly deallocate the integral grid information
 ! and the scf mixing information
 !
 call destroy_scf()
 if( calc_type%is_dft ) call destroy_dft_grid()


 !
 ! Get the exchange operator if not already calculated
 !
 if( .NOT. has_auxil_basis ) then
   if( ABS(en%exx) < 1.0e-6_dp ) call setup_exchange(print_matrix_,basis%nbf,p_matrix,hamiltonian_exx,en%exx)
 else
   if( ABS(en%exx) < 1.0e-6_dp ) then
     if( parallel_ham ) then
       if( parallel_buffer ) then
         call setup_exchange_ri_buffer_sca(basis%nbf,nstate,m_c,n_c,m_ham,n_ham,occupation,c_matrix,p_matrix,hamiltonian_exx,en%exx)
       else
         call die('Exchange with fully distributed hamiltonian: case not implemented yet')
       endif
     else
       call setup_exchange_ri(basis%nbf,nstate,occupation,c_matrix,p_matrix,hamiltonian_exx,en%exx)
     endif
   endif
 endif

 !
 ! Obtain the Fock matrix and store it
 !
 call clean_allocate('Fock operator F',hamiltonian_fock,basis%nbf,basis%nbf,nspin)
 call get_fock_operator(hamiltonian,hamiltonian_xc,hamiltonian_exx,hamiltonian_fock)

 !
 ! Cleanly deallocate the arrays
 !
 call clean_deallocate('Density matrix P',p_matrix)
 call clean_deallocate('Total Hamiltonian H',hamiltonian)
 call clean_deallocate('Hartree potential Vh',hamiltonian_hartree)
 call clean_deallocate('Exchange operator Sigx',hamiltonian_exx)
 call clean_deallocate('XC operator Vxc',hamiltonian_xc)

 write(stdout,'(/,/,a25,1x,f19.10,/)') 'SCF Total Energy (Ha):',en%tot
 write(stdout,'(a25,1x,f19.10)')       '      EXX Energy (Ha):',en%exx
 write(stdout,'(a25,1x,f19.10)')       'Total EXX Energy (Ha):',en%nuc_nuc + en%kin + en%nuc + en%hart + en%exx


 !
 ! Deallocate the buffer here
 if( parallel_ham .AND. parallel_buffer ) call destroy_parallel_buffer()
 !
 ! At this point, all procs get the complete c_matrix
 !
 call form_c_matrix_global(basis%nbf,nstate,c_matrix)


 !
 ! Single excitation term
 !
 call single_excitations(nstate,basis%nbf,energy,occupation,c_matrix,hamiltonian_fock,en%se)
 write(stdout,'(a25,1x,f19.10)') 'Singles correction (Ha):',en%se
 write(stdout,'(a25,1x,f19.10,/)')   'Est. HF Energy (Ha):',en%nuc_nuc + en%kin + en%nuc + en%hart + en%exx + en%se

 !
 ! Evaluate spin contamination
 call evaluate_s2_operator(basis%nbf,nstate,occupation,c_matrix,s_matrix)


 ! A dirty section for the Luttinger-Ward functional
 if(calc_type%selfenergy_approx==LW .OR. calc_type%selfenergy_approx==LW2 .OR. calc_type%selfenergy_approx==GSIGMA) then
   allocate(energy_exx(nstate,nspin))
   allocate(c_matrix_exx(basis%nbf,nstate,nspin))
   call issue_warning('ugly coding here write temp file fort.1000 and fort.1001')

   do ispin=1,nspin
     write(stdout,*) 'Diagonalization H_exx for spin channel',ispin
     call diagonalize_generalized_sym(basis%nbf,&
                                      hamiltonian_fock(:,:,ispin),s_matrix(:,:),&
                                      energy_exx(:,ispin),c_matrix_exx(:,:,ispin))
   enddo
   write(stdout,*) 'FBFB LW sum(      epsilon) + Eii -<vxc> - EH + Ex',en%nuc_nuc + en%kin + en%nuc + en%hart + en%exx
   write(stdout,*) 'FBFB LW sum(tilde epsilon) + Eii - EH - Ex       ',SUM( occupation(:,:)*energy_exx(:,:) ) + en%nuc_nuc - en%hart - en%exx
   open(1000,form='unformatted')
   do ispin=1,nspin
     do istate=1,nstate
       write(1000) c_matrix_exx(:,istate,ispin)
     enddo
   enddo
   close(1000)
   open(1001,form='unformatted')
   write(1001) energy_exx(:,:)
   close(1001)
   deallocate(energy_exx,c_matrix_exx)
 endif



 !
 ! Big RESTART file written if converged
 ! 
 if( is_converged .AND. print_bigrestart_ ) then
   call write_restart(BIG_RESTART,basis,nstate,occupation,c_matrix,energy,hamiltonian_fock)
 else  
   if( print_restart_ ) then
     call write_restart(SMALL_RESTART,basis,nstate,occupation,c_matrix,energy,hamiltonian_fock)
   endif
 endif


 call stop_clock(timing_scf)

end subroutine scf_loop


!=========================================================================
subroutine calculate_hamiltonian_hxc_ri(basis,nstate,m_ham,n_ham,m_c,n_c,occupation,c_matrix,p_matrix,hamiltonian_hxc)
 use m_scalapack
 use m_basis_set
 use m_hamiltonian
 use m_hamiltonian_sca
 use m_hamiltonian_buffer
 implicit none

 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: m_ham,n_ham
 integer,intent(in)         :: nstate
 integer,intent(in)         :: m_c,n_c
 real(dp),intent(in)        :: occupation(nstate,nspin)
 real(dp),intent(in)        :: c_matrix(m_c,n_c,nspin)
 real(dp),intent(in)        :: p_matrix(m_ham,n_ham,nspin)
 real(dp),intent(out)       :: hamiltonian_hxc(m_ham,n_ham,nspin)
!=====
 integer  :: ispin
 real(dp) :: hamiltonian_tmp(m_ham,n_ham)
 real(dp) :: hamiltonian_spin_tmp(m_ham,n_ham,nspin)
 real(dp) :: ehart,exc,eexx,eexx_hyb
!=====



 !
 ! Hartree contribution to the Hamiltonian
 !
 if( parallel_ham ) then
   if( parallel_buffer ) then
     call setup_hartree_ri_buffer_sca(print_matrix_,basis%nbf,m_ham,n_ham,p_matrix,hamiltonian_tmp,ehart)
   else
     call setup_hartree_ri_sca(print_matrix_,basis%nbf,m_ham,n_ham,p_matrix,hamiltonian_tmp,ehart)
   endif
 else
   call setup_hartree_ri(print_matrix_,basis%nbf,p_matrix,hamiltonian_tmp,ehart)
 endif

 do ispin=1,nspin
   hamiltonian_hxc(:,:,ispin) = hamiltonian_tmp(:,:)
 enddo


 !
 !  XC part of the Hamiltonian
 !

 !
 ! DFT XC potential is added here
 ! 
 if( calc_type%is_dft ) then
   hamiltonian_spin_tmp(:,:,:) = 0.0_dp

   if( parallel_ham ) then
     if( parallel_buffer ) then
       call dft_exc_vxc_buffer_sca(basis,nstate,m_c,n_c,m_ham,n_ham,occupation,c_matrix,p_matrix,hamiltonian_spin_tmp,exc)
     else
       call issue_warning('Exc calculation with SCALAPACK is not coded yet. Just skip it')
       hamiltonian_spin_tmp(:,:,:) = 0.0_dp
       exc = 0.0_dp
     endif
   else
     call dft_exc_vxc_batch(BATCH_SIZE,basis,nstate,occupation,c_matrix,hamiltonian_spin_tmp,exc)
   endif

   hamiltonian_hxc(:,:,:) = hamiltonian_hxc(:,:,:) + hamiltonian_spin_tmp(:,:,:) 
 endif


 !
 ! LR Exchange contribution to the Hamiltonian
 !
 if(calc_type%need_exchange_lr) then
   hamiltonian_spin_tmp(:,:,:) = 0.0_dp

   if( parallel_ham ) then
     if( parallel_buffer ) then
       call setup_exchange_longrange_ri_buffer_sca(basis%nbf,nstate,m_c,n_c,m_ham,n_ham,occupation,c_matrix,p_matrix,hamiltonian_spin_tmp,eexx)
     else
       call die('Range-separated functionals not implemented with full SCALAPACK yet')
     endif
   else
     call setup_exchange_longrange_ri(basis%nbf,nstate,occupation,c_matrix,p_matrix,hamiltonian_spin_tmp,eexx)
   endif
   ! Rescale with alpha_hybrid_lr for range-separated hybrid functionals
   eexx_hyb = alpha_hybrid_lr * eexx
   hamiltonian_hxc(:,:,:) = hamiltonian_hxc(:,:,:) + hamiltonian_spin_tmp(:,:,:) * alpha_hybrid_lr
 endif


 !
 ! Exchange contribution to the Hamiltonian
 !
 if( calc_type%need_exchange ) then
   hamiltonian_spin_tmp(:,:,:) = 0.0_dp

   if( parallel_ham ) then
     if( parallel_buffer ) then
       call setup_exchange_ri_buffer_sca(basis%nbf,nstate,m_c,n_c,m_ham,n_ham,occupation,c_matrix,p_matrix,hamiltonian_spin_tmp,eexx)
     else
       call die('Exchange with fully distributed hamiltonian: case not implemented yet')
     endif
   else
     call setup_exchange_ri(basis%nbf,nstate,occupation,c_matrix,p_matrix,hamiltonian_spin_tmp,eexx)
   endif
   ! Rescale with alpha_hybrid for hybrid functionals
   eexx_hyb = eexx_hyb + alpha_hybrid * eexx
   hamiltonian_hxc(:,:,:) = hamiltonian_hxc(:,:,:) + hamiltonian_spin_tmp(:,:,:) * alpha_hybrid
 endif




end subroutine  calculate_hamiltonian_hxc_ri


!=========================================================================
subroutine get_fock_operator(hamiltonian,hamiltonian_xc,hamiltonian_exx,hamiltonian_fock)
 use m_mpi
 use m_scalapack
 implicit none

 real(dp),intent(in)    :: hamiltonian(:,:,:)
 real(dp),intent(in)    :: hamiltonian_xc(:,:,:)
 real(dp),intent(in)    :: hamiltonian_exx(:,:,:)
 real(dp),intent(out)   :: hamiltonian_fock(:,:,:)
!=====
 real(dp),allocatable   :: hfock_local(:,:,:)
!=====

 if( parallel_ham ) then

   call clean_allocate('Local Fock operator F',hfock_local,SIZE(hamiltonian,DIM=1),SIZE(hamiltonian,DIM=2),nspin)
   hfock_local(:,:,:) = hamiltonian(:,:,:) - hamiltonian_xc(:,:,:) + hamiltonian_exx(:,:,:)
   call gather_distributed_copy(desc_ham,hfock_local,hamiltonian_fock)
   call clean_deallocate('Local Fock operator F',hfock_local)

 else
   hamiltonian_fock(:,:,:) = hamiltonian(:,:,:) - hamiltonian_xc(:,:,:) + hamiltonian_exx(:,:,:)

 endif

end subroutine get_fock_operator


!=========================================================================
subroutine form_c_matrix_global(nbf,nstate,c_matrix)
 use m_mpi
 use m_scalapack
 implicit none

 integer,intent(in)                 :: nbf,nstate
 real(dp),allocatable,intent(inout) :: c_matrix(:,:,:)
!=====
 real(dp),allocatable :: c_matrix_local(:,:,:)
!=====

 if( .NOT. parallel_ham ) return    ! Nothing to do

 write(stdout,'(/,1x,a)') 'Form the C matrix on all procs'


 call clean_allocate('Local wfn coeff C',c_matrix_local,SIZE(c_matrix,DIM=1),SIZE(c_matrix,DIM=2),nspin)
! if( cntxt_ham > 0 ) then
 c_matrix_local(:,:,:) = c_matrix(:,:,:)
! endif
 call clean_deallocate('Wavefunctions C',c_matrix)
 call clean_allocate('Wavefunctions C',c_matrix,nbf,nstate,nspin)

 call gather_distributed_copy(desc_c,c_matrix_local,c_matrix)

 call clean_deallocate('Local wfn coeff C',c_matrix_local)


 write(stdout,'(1x,a)') 'C matrix on all procs formed'

end subroutine form_c_matrix_global



!=========================================================================
end module m_scf_loop
!=========================================================================
