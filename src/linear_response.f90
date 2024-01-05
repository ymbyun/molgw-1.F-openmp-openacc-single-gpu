!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains
! the routines to calculate the polarizability within RPA, TDDFT or BSE
! and the corresponding optical spectra
!
!=========================================================================


!=========================================================================
subroutine polarizability(basis,nstate,occupation,energy,c_matrix,rpa_correlation,wpol_out)
 use m_definitions
 use m_timing
 use m_warning
 use m_memory
 use m_inputparam
 use m_mpi
 use m_scalapack
 use m_cart_to_pure
 use m_block_diago
 use m_basis_set
 use m_spectral_function
 use m_eri_ao_mo
 implicit none

 type(basis_set),intent(in)            :: basis
 integer,intent(in)                    :: nstate
 real(dp),intent(in)                   :: occupation(nstate,nspin)
 real(dp),intent(in)                   :: energy(nstate,nspin),c_matrix(basis%nbf,nstate,nspin)
 real(dp),intent(out)                  :: rpa_correlation
 type(spectral_function),intent(inout) :: wpol_out
!=====
 type(spectral_function)   :: wpol_static
 integer                   :: nmat,nexc
 real(dp)                  :: energy_gm
 real(dp)                  :: alpha_local
 real(dp),allocatable      :: amb_diag_rpa(:)
 real(dp),allocatable      :: amb_matrix(:,:),apb_matrix(:,:)
 real(dp),allocatable      :: a_diag(:)
 real(dp),allocatable      :: xpy_matrix(:,:),xmy_matrix(:,:)
 real(dp),allocatable      :: eigenvalue(:)
 real(dp)                  :: energy_qp(nstate,nspin)
 logical                   :: is_tddft
 logical                   :: is_ij
 logical                   :: is_rpa
#ifdef ENABLE_YMBYUN
 ! ymbyun 2018/06/19
 logical                   :: skip_fxc
 ! ymbyun 2019/07/31
 logical                   :: write_mo_integrals
 integer                   :: mo_integrals_file
! integer,parameter         :: iajb_max = 100
! integer                   :: iajb_max_read,iajb_i
 integer                   :: iajb_max,iajb_i
 integer                   :: istate,astate,jstate,bstate
 integer                   :: iaspin,jbspin
 integer                   :: jbmin
 real(dp),allocatable      :: eri_eigenstate_jbmin(:,:,:,:)
 real(dp)                  :: eri_eigen_iajb
 logical                   :: j_is_jbmin
#endif
 logical                   :: has_manual_tdhf
 integer                   :: reading_status
 integer                   :: tdhffile
 integer                   :: m_apb,n_apb,m_x,n_x
! Scalapack variables
 integer                   :: desc_apb(NDEL),desc_x(NDEL)
!=====

 call start_clock(timing_pola)

 write(stdout,'(/,a)') ' Calculating the polarizability'
 if(is_triplet) then
   write(stdout,'(a)') ' Triplet final state'
 else
   write(stdout,'(a)') ' Singlet final state'
 endif

 if( has_auxil_basis ) call calculate_eri_3center_eigen(basis%nbf,nstate,c_matrix,ncore_W+1,nvirtual_W-1,ncore_W+1,nvirtual_W-1)

 ! Set up all the switches to be able to treat
 ! GW, BSE, TDHF, TDDFT (semilocal or hybrid)

 !
 ! Set up flag is_rpa
 inquire(file='manual_rpa',exist=is_rpa)
 if(is_rpa) then
   msg='RPA calculation is enforced'
   call issue_warning(msg)
 endif
 is_rpa   = calc_type%is_gw .OR. is_rpa

 ! ymbyun 2018/06/19
 ! This is quick and dirty.
#ifdef ENABLE_YMBYUN
 inquire(file='manual_rpa_optical_spectrum',exist=skip_fxc)
 if(skip_fxc) then
   ! ymbyun 2018/09/15
   msg='(ymbyun) RPA optical spectrum is enforced. fxc for hybrid functionals is set to zero correctly now.'
   call issue_warning(msg)
 endif
#endif

 ! 
 ! Set up flag is_tddft
 is_tddft = calc_type%is_td .AND. calc_type%is_dft .AND. .NOT. is_rpa

 ! 
 ! Set up exchange content alpha_local
 ! manual_tdhf can override anything
 inquire(file='manual_tdhf',exist=has_manual_tdhf)
 if(has_manual_tdhf) then
   open(newunit=tdhffile,file='manual_tdhf',status='old')
   read(tdhffile,*) alpha_local
   close(tdhffile)
   write(msg,'(a,f12.6,3x,f12.6)') 'calculating the TDHF polarizability with alpha ',alpha_local
   call issue_warning(msg)
 else
   if(is_rpa) then
     alpha_local = 0.0_dp
   else if(is_tddft) then
   ! ymbyun 2018/09/15
#ifdef ENABLE_YMBYUN
     if(skip_fxc) then
       alpha_local = 0.0_dp
     else
       alpha_local = alpha_hybrid
     endif
#else
     alpha_local = alpha_hybrid
#endif
   else ! TDHF or BSE case
     alpha_local = 1.0_dp
   endif
 endif


 call start_clock(timing_build_h2p)

 !
 ! Prepare the QP energies
 !
 if( calc_type%is_bse ) then
   ! Get energy_qp 
   call get_energy_qp(nstate,energy,occupation,energy_qp)
 else
   ! For any other type of calculation, just fill energy_qp array with energy
   energy_qp(:,:) = energy(:,:)
 endif

 ! 
 ! BSE needs the static screening from a previous calculation
 ! It is stored in object wpol_static
 !
 if( calc_type%is_bse ) then
   call init_spectral_function(nstate,occupation,1,wpol_static)
   call read_spectral_function(wpol_static,reading_status)

   ! If a SCREENED_COULOMB file cannot be found,
   ! then recalculate it from scratch
   if( reading_status /= 0 ) then
     wpol_static%nprodbasis = nauxil_3center
     call static_polarizability(nstate,occupation,energy,wpol_static)
   endif

 endif

 !
 ! Prepare the big matrices (A+B) and (A-B)
 ! 
 nmat = wpol_out%npole_reso_apb
 !
 ! The distribution of the two matrices have to be the same for A-B and A+B
 ! This is valid also when SCALAPACK is not used!
 call init_desc('S',nmat,nmat,desc_apb,m_apb,n_apb)
 call clean_allocate('A+B',apb_matrix,m_apb,n_apb)
 call clean_allocate('A-B',amb_matrix,m_apb,n_apb)
 allocate(amb_diag_rpa(nmat))

 ! A diagonal is owned by all procs (= no distribution)
 ! wpol_out%npole_reso_spa are the pole not explictely counted in wpol_out%npole_reso_apb
 allocate(a_diag(wpol_out%npole_reso_spa))

 !
 ! Build the (A+B) and (A-B) matrices in 3 steps
 ! to span all the possible approximations
 ! Only the lower triangle is calculated
 ! the upper part will be filled later by symmetry
 !

 ! Calculate the diagonal separately: it is needed for the single pole approximation
 if( nvirtual_SPA < nvirtual_W .AND. is_rpa ) & 
     call build_a_diag_common(nmat,basis%nbf,nstate,c_matrix,energy_qp,wpol_out,a_diag)

 apb_matrix(:,:) = 0.0_dp
 amb_matrix(:,:) = 0.0_dp
 write(stdout,'(/,a)') ' Build the electron-hole hamiltonian'

 ! ymbyun 2019/07/31
#ifdef ENABLE_YMBYUN
 inquire(file='manual_mo_integrals',exist=write_mo_integrals)
 if(write_mo_integrals) then
   call issue_warning('(ymbyun) Writing MO integrals')

   open(newunit=mo_integrals_file,file='manual_mo_integrals',form='formatted',status='old')

   read(mo_integrals_file,*) iajb_max
   if( iajb_max <= 0 ) then
     call issue_warning('(ymbyun) The number of MO integrals to write should be greater than 0')
   else
     if(.NOT. has_auxil_basis) then
       allocate(eri_eigenstate_jbmin(nstate,nstate,nstate,nspin))
       eri_eigenstate_jbmin(:,:,:,:) = 0.0_dp
     endif

     write(stdout,*)
     do iajb_i=1,iajb_max
       read(mo_integrals_file,*) istate,astate,iaspin,jstate,bstate,jbspin
       if(has_auxil_basis) then
         ! ymbyun 2019/08/01
         ! MPI reduction is needed.
!         eri_eigen_iajb = eri_eigen_ri_paral(istate,astate,iaspin,jstate,bstate,jbspin)
         eri_eigen_iajb = eri_eigen_ri(istate,astate,iaspin,jstate,bstate,jbspin)
       else
         jbmin = MIN(jstate,bstate)
         j_is_jbmin = (jbmin == jstate)
         call calculate_eri_4center_eigen(basis%nbf,nstate,c_matrix,jbmin,jbspin,eri_eigenstate_jbmin)
         if(j_is_jbmin) then ! treating (j,b)
           eri_eigen_iajb = eri_eigenstate_jbmin(bstate,istate,astate,iaspin)
         else                ! treating (b,j)
           eri_eigen_iajb = eri_eigenstate_jbmin(jstate,istate,astate,iaspin)
         endif
       endif
       write(stdout,'(a,6(2x,i4),2x,f16.8)') ' (i a sigma | j b sigma'') ',istate,astate,iaspin,jstate,bstate,jbspin,eri_eigen_iajb
     enddo
     write(stdout,*)

     if(ALLOCATED(eri_eigenstate_jbmin)) deallocate(eri_eigenstate_jbmin)
   endif

   close(mo_integrals_file)
 endif
#endif

 if( has_auxil_basis) then

   !
   ! Step 1
   call build_amb_apb_diag_auxil(nmat,nstate,energy_qp,wpol_out,m_apb,n_apb,amb_matrix,apb_matrix,amb_diag_rpa)

#ifdef HAVE_SCALAPACK
   call build_apb_hartree_auxil_scalapack(desc_apb,wpol_out,m_apb,n_apb,apb_matrix)
#else
   call build_apb_hartree_auxil(desc_apb,wpol_out,m_apb,n_apb,apb_matrix)
#endif

   call get_rpa_correlation(nmat,m_apb,n_apb,amb_matrix,apb_matrix,rpa_correlation)

   !
   ! Step 2
   ! ymbyun 2018/06/19
#ifdef ENABLE_YMBYUN
   if(is_tddft .AND. .NOT. skip_fxc) call build_apb_tddft(nmat,nstate,basis,c_matrix,occupation,wpol_out,m_apb,n_apb,apb_matrix)
#else
   if(is_tddft)                      call build_apb_tddft(nmat,nstate,basis,c_matrix,occupation,wpol_out,m_apb,n_apb,apb_matrix)
#endif

   !
   ! Step 3
   if(alpha_local > 1.0e-6_dp) then
     call build_amb_apb_screened_exchange_auxil(alpha_local,desc_apb,wpol_out,wpol_static,m_apb,n_apb,amb_matrix,apb_matrix)
   endif

   if(calc_type%is_bse) then
     call destroy_spectral_function(wpol_static)
   endif


 else

   !
   ! Step 1
   call build_amb_apb_common(nmat,basis%nbf,nstate,c_matrix,energy_qp,wpol_out,alpha_local, &
                             m_apb,n_apb,amb_matrix,apb_matrix,amb_diag_rpa,rpa_correlation)

   !
   ! Step 2
   ! ymbyun 2018/06/19
#ifdef ENABLE_YMBYUN
   if(is_tddft .AND. .NOT. skip_fxc) call build_apb_tddft(nmat,nstate,basis,c_matrix,occupation,wpol_out,m_apb,n_apb,apb_matrix)
#else
   if(is_tddft)                      call build_apb_tddft(nmat,nstate,basis,c_matrix,occupation,wpol_out,m_apb,n_apb,apb_matrix)
#endif


   !
   ! Step 3
   if(calc_type%is_bse .AND. .NOT. is_rpa) then
     call build_amb_apb_bse(basis%nbf,nstate,wpol_out,wpol_static,m_apb,n_apb,amb_matrix,apb_matrix)
     call destroy_spectral_function(wpol_static)
   endif

 endif



 ! Warning if Tamm-Dancoff flag is on
 if(is_tda) then
   msg='Tamm-Dancoff approximation is switched on'
   call issue_warning(msg)
   ! Tamm-Dancoff approximation consists in setting B matrix to zero
   ! Then A+B = A-B = A
   apb_matrix(:,:) = 0.5_dp * ( apb_matrix(:,:) + amb_matrix(:,:) )
   amb_matrix(:,:) = apb_matrix(:,:) 
 endif
 ! Construction done!
 if(has_auxil_basis) call destroy_eri_3center_eigen()

 call stop_clock(timing_build_h2p)

 if( is_rpa .AND. .NOT. is_tda ) call clean_deallocate('A-B',amb_matrix)
 

 !
 ! Prepare the second dimension of xpy_matrix and xmy_matrix
 nexc = nexcitation
 if( nexc == 0 ) nexc = wpol_out%npole_reso_apb

 allocate(eigenvalue(nmat))

 ! Allocate (X+Y) 
 ! Allocate (X-Y) only if actually needed
 call init_desc('S',nmat,nexc,desc_x,m_x,n_x)

 call clean_allocate('X+Y',xpy_matrix,m_x,n_x)
 if( .NOT. is_rpa .OR. is_tda ) &
   call clean_allocate('X-Y',xmy_matrix,m_x,n_x)

 !
 ! Diago using the 4 block structure and the symmetry of each block
 ! With or Without SCALAPACK
 !
 if( .NOT. is_rpa .OR. is_tda ) then
   if( nexcitation == 0 ) then
     ! The following call works with AND without SCALAPACK
     call diago_4blocks_chol(nmat,desc_apb,m_apb,n_apb,amb_matrix,apb_matrix,eigenvalue,&
                             desc_x,m_x,n_x,xpy_matrix,xmy_matrix)

   else ! Partial diagonalization with Davidson
     ! The following call works with AND without SCALAPACK
     call diago_4blocks_davidson(toldav,nexcitation,nmat,amb_diag_rpa,desc_apb,m_apb,n_apb,amb_matrix,apb_matrix, & 
                                 eigenvalue,desc_x,m_x,n_x,xpy_matrix,xmy_matrix)
   endif
 else
   ! The following call works with AND without SCALAPACK
   call diago_4blocks_rpa_sca(nmat,desc_apb,m_apb,n_apb,amb_diag_rpa,apb_matrix,eigenvalue,&
                              desc_x,m_x,n_x,xpy_matrix)
 endif


 ! Deallocate the non-necessary matrices
 deallocate(amb_diag_rpa)
 write(stdout,*) 'Deallocate (A+B) and possibly (A-B)'
 call clean_deallocate('A+B',apb_matrix)
 !
 ! (A-B) may have been already deallocated earlier in the case of RPA 
 ! Relax: this is indeed tolerated by clean_deallocate
 call clean_deallocate('A-B',amb_matrix)


 !
 ! Second part of the RPA correlation energy: sum over positive eigenvalues
 rpa_correlation = rpa_correlation + 0.50_dp * SUM( ABS(eigenvalue(:)) )
 if(is_rpa) then
   write(stdout,'(/,a)') ' Calculate the RPA energy using the Tamm-Dancoff decomposition'
   write(stdout,'(a)')   ' Eq. (9) from J. Chem. Phys. 132, 234114 (2010)'
   write(stdout,'(/,a,f16.10)') ' RPA correlation energy (Ha): ',rpa_correlation
 endif

 write(stdout,'(/,a,f12.6)') ' Lowest neutral excitation energy (eV):',MINVAL(ABS(eigenvalue(1:nexc)))*Ha_eV

 if( has_auxil_basis ) call calculate_eri_3center_eigen(basis%nbf,nstate,c_matrix,ncore_W+1,nhomo_W,nlumo_W,nvirtual_W-1)

 !
 ! Calculate the optical sprectrum
 ! and the dynamic dipole tensor
 !
 if( calc_type%is_td .OR. calc_type%is_bse ) then
   call optical_spectrum(nstate,basis,occupation,c_matrix,wpol_out,m_x,n_x,xpy_matrix,xmy_matrix,eigenvalue)
!   call stopping_power(nstate,basis,occupation,c_matrix,wpol_out,m_x,n_x,xpy_matrix,eigenvalue)
 endif

 !
 ! Now only the sum ( X + Y ) is needed in fact.
 ! Free ( X - Y ) at once.
 call clean_deallocate('X-Y',xmy_matrix)

 !
 ! Calculate Wp= v * chi * v    if necessary
 ! and then write it down on file
 !
 if( print_w_ .OR. calc_type%is_gw ) then
   if( has_auxil_basis) then
     call chi_to_sqrtvchisqrtv_auxil(basis%nbf,desc_x,m_x,n_x,xpy_matrix,eigenvalue,wpol_out,energy_gm)
     ! This following coding of the Galitskii-Migdal correlation energy is only working with
     ! an auxiliary basis
     if(is_rpa) write(stdout,'(a,f16.10,/)') ' Correlation energy in the Galitskii-Migdal formula (Ha): ',energy_gm
     
     ! Add the single pole approximation for the poles that have been neglected
     ! in the diagonalization
     if( nvirtual_SPA < nvirtual_W .AND. is_rpa ) & 
        call chi_to_sqrtvchisqrtv_auxil_spa(basis%nbf,a_diag,wpol_out)

   else
     call chi_to_vchiv(basis%nbf,nstate,c_matrix,xpy_matrix,eigenvalue,wpol_out)
   endif
  
 
   ! If requested write the spectral function on file
   if( print_w_ ) call write_spectral_function(wpol_out)

 endif

 if( .NOT. calc_type%is_gw ) call destroy_spectral_function(wpol_out)

 write(stdout,*) 'Deallocate eigenvector array'
 call clean_deallocate('X+Y',xpy_matrix)

 if(has_auxil_basis) call destroy_eri_3center_eigen()

 if(ALLOCATED(eigenvalue)) deallocate(eigenvalue)
 if(ALLOCATED(a_diag))     deallocate(a_diag)

 call stop_clock(timing_pola)


end subroutine polarizability


!=========================================================================
subroutine polarizability_onering(basis,nstate,occupation,energy,c_matrix,vchi0v)
 use m_definitions
 use m_timing
 use m_warning
 use m_memory
 use m_inputparam
 use m_mpi
 use m_scalapack
 use m_cart_to_pure
 use m_block_diago
 use m_basis_set
 use m_spectral_function
 use m_eri_ao_mo
 implicit none

 type(basis_set),intent(in)            :: basis
 integer,intent(in)                    :: nstate
 real(dp),intent(in)                   :: occupation(nstate,nspin)
 real(dp),intent(in)                   :: energy(nstate,nspin),c_matrix(basis%nbf,nstate,nspin)
 type(spectral_function),intent(inout) :: vchi0v
!=====
 integer :: t_jb
 integer :: jstate,bstate,jbspin
!=====

 call allocate_spectral_function(nauxil_3center,vchi0v)

 call calculate_eri_3center_eigen(basis%nbf,nstate,c_matrix,ncore_W+1,nhomo_W,nlumo_W,nvirtual_W-1)


 do t_jb=1,vchi0v%npole_reso
   jstate = vchi0v%transition_table_apb(1,t_jb)
   bstate = vchi0v%transition_table_apb(2,t_jb)
   jbspin = vchi0v%transition_table_apb(3,t_jb)

   vchi0v%residue_left(:,t_jb) = eri_3center_eigen(:,jstate,bstate,jbspin) * SQRT(spin_fact)
   vchi0v%pole(t_jb)           = energy(bstate,jbspin) - energy(jstate,jbspin)

 enddo

 call destroy_eri_3center_eigen()

end subroutine polarizability_onering


!=========================================================================
subroutine optical_spectrum(nstate,basis,occupation,c_matrix,chi,m_x,n_x,xpy_matrix,xmy_matrix,eigenvalue)
 use m_definitions
 use m_timing
 use m_warning
 use m_memory
 use m_mpi
 use m_scalapack
 use m_inputparam
 use m_basis_set
 use m_dft_grid
 use m_spectral_function
 use m_hamiltonian_onebody
 implicit none

 integer,intent(in)                 :: nstate,m_x,n_x
 type(basis_set),intent(in)         :: basis
 real(dp),intent(in)                :: occupation(nstate,nspin),c_matrix(basis%nbf,nstate,nspin)
 type(spectral_function),intent(in) :: chi
 real(dp),intent(in)                :: xpy_matrix(m_x,n_x)
 real(dp),intent(in)                :: xmy_matrix(m_x,n_x)
 real(dp),intent(in)                :: eigenvalue(chi%npole_reso_apb)
!=====
 integer                            :: gt
 integer                            :: nexc
 integer                            :: t_ia,t_jb
 integer                            :: t_ia_global,t_jb_global
 integer                            :: istate,astate,iaspin
 integer                            :: mstate,pstate,mpspin
 integer                            :: ibf,jbf
 integer                            :: iomega,idir,jdir
 integer,parameter                  :: nomega=600
 complex(dp)                        :: omega(nomega)
 real(dp)                           :: coeff(2*chi%npole_reso_apb),trace
 real(dp)                           :: dynamical_pol(nomega,3,3),photoabsorp_cross(nomega,3,3)
 real(dp)                           :: static_polarizability(3,3)
 real(dp)                           :: oscillator_strength,trk_sumrule,mean_excitation
 real(dp),allocatable               :: dipole_basis(:,:,:),dipole_state(:,:,:,:)
 real(dp),allocatable               :: dipole_cart(:,:,:)
 real(dp),allocatable               :: residue(:,:)
 integer                            :: dynpolfile
 integer                            :: photocrossfile
 integer                            :: parityi,parityj,reflectioni,reflectionj
 integer,external                   :: wfn_parity
 integer,external                   :: wfn_reflection
 character(len=32)                  :: symsymbol
!=====


 call start_clock(timing_spectrum)
 !
 ! Calculate the spectrum now
 !

 write(stdout,'(/,a)') ' Calculate the optical spectrum'

 gt = get_gaussian_type_tag(basis%gaussian_type)

 nexc = nexcitation
 if( nexc == 0 ) nexc = chi%npole_reso_apb

 !
 ! First precalculate all the needed dipole in the basis set
 !
 call calculate_dipole_basis(basis,dipole_basis)

 !
 ! Get the dipole oscillator strength on states
 allocate(dipole_state(nstate,nstate,nspin,3))

 do idir=1,3
   do mpspin=1,nspin
     dipole_state(:,:,mpspin,idir) = MATMUL( TRANSPOSE( c_matrix(:,:,mpspin) ) , &
                                             MATMUL( dipole_basis(:,:,idir) , c_matrix(:,:,mpspin) ) )
   enddo
 enddo

 deallocate(dipole_basis)


 allocate(residue(3,nexc))

 residue(:,:) = 0.0_dp
 do t_ia=1,m_x
   t_ia_global = rowindex_local_to_global(iprow_sd,nprow_sd,t_ia)
   istate = chi%transition_table_apb(1,t_ia_global)
   astate = chi%transition_table_apb(2,t_ia_global)
   iaspin = chi%transition_table_apb(3,t_ia_global)

   ! Let use (i <-> j) symmetry to halve the loop
   do t_jb=1,n_x
     t_jb_global = colindex_local_to_global(ipcol_sd,npcol_sd,t_jb)

     if( t_jb_global <= nexc) then
       residue(:,t_jb_global) = residue(:,t_jb_global) &
                    + dipole_state(istate,astate,iaspin,:) * xpy_matrix(t_ia,t_jb) * SQRT(spin_fact)
     endif
   enddo

 enddo
 call xsum_world(residue)

 deallocate(dipole_state)



 write(stdout,'(/,a)') ' Excitation energies (eV)     Oscil. strengths   [Symmetry] '  
 trk_sumrule=0.0_dp
 mean_excitation=0.0_dp
 do t_jb_global=1,nexc
   t_jb = colindex_global_to_local('S',t_jb_global)

   if( is_triplet ) then 
     oscillator_strength = 0.0_dp
   else
     oscillator_strength = 2.0_dp/3.0_dp * DOT_PRODUCT(residue(:,t_jb_global),residue(:,t_jb_global)) * eigenvalue(t_jb_global)
   endif
   trk_sumrule = trk_sumrule + oscillator_strength
   mean_excitation = mean_excitation + oscillator_strength * LOG( eigenvalue(t_jb_global) )

   if(t_jb_global<=30) then

     if( is_triplet ) then
       symsymbol='3'
     else
       symsymbol='1'
     endif

     !
     ! Test the parity in case of molecule with inversion symmetry
    
     t_ia_global = 0
     do t_ia=1,m_x
       ! t_jb is zero if the proc is not in charge of this process
       if( t_jb /=0 ) then 
         if( 0.5_dp * ABS( xpy_matrix(t_ia,t_jb) + xmy_matrix(t_ia,t_jb) ) > 0.1_dp ) then
           t_ia_global = rowindex_local_to_global(iprow_sd,nprow_sd,t_ia)
           exit
         endif
       endif
     enddo
     call xmax_world(t_ia_global)
     if( t_ia_global == 0 ) cycle

     istate = chi%transition_table_apb(1,t_ia_global)
     astate = chi%transition_table_apb(2,t_ia_global)
     iaspin = chi%transition_table_apb(3,t_ia_global)
     if(planar) then
       reflectioni = wfn_reflection(nstate,basis,c_matrix,istate,iaspin)
       reflectionj = wfn_reflection(nstate,basis,c_matrix,astate,iaspin)
       select case(reflectioni*reflectionj)
       case( 1)
         symsymbol=TRIM(symsymbol)//'(A1, B2 or Ap )'
       case(-1)
         symsymbol=TRIM(symsymbol)//'(A2, B1 or App)'
       end select
     endif
     if(inversion) then
       parityi = wfn_parity(nstate,basis,c_matrix,istate,iaspin)
       parityj = wfn_parity(nstate,basis,c_matrix,astate,iaspin)
       select case(parityi*parityj)
       case( 1)
         symsymbol=TRIM(symsymbol)//'g'
       case(-1)
         symsymbol=TRIM(symsymbol)//'u'
       end select
     endif

     write(stdout,'(1x,i4.4,a3,2(f18.8,2x),5x,a32)') t_jb_global,' : ', &
                  eigenvalue(t_jb_global)*Ha_eV,oscillator_strength,symsymbol

     !
     ! Output the transition coefficients
     coeff(:) = 0.0_dp
     do t_ia=1,m_x
       t_ia_global = rowindex_local_to_global('S',t_ia)
       istate = chi%transition_table_apb(1,t_ia_global)
       astate = chi%transition_table_apb(2,t_ia_global)
       if( t_jb /= 0 ) then
         ! Resonant
         coeff(                     t_ia_global) = 0.5_dp * ( xpy_matrix(t_ia,t_jb) + xmy_matrix(t_ia,t_jb) ) / SQRT(2.0_dp)
         ! Anti-Resonant
         coeff(chi%npole_reso_apb + t_ia_global) = 0.5_dp * ( xpy_matrix(t_ia,t_jb) - xmy_matrix(t_ia,t_jb) ) / SQRT(2.0_dp)
       endif
     enddo
     call xsum_world(coeff)
     

     do t_ia_global=1,chi%npole_reso_apb
       istate = chi%transition_table_apb(1,t_ia_global)
       astate = chi%transition_table_apb(2,t_ia_global)
       ! Resonant
       if( ABS(coeff(                   t_ia_global)) > 0.1_dp )  &
         write(stdout,'(8x,i4,a,i4,1x,f12.5)') istate,' -> ',astate,coeff(t_ia_global)
       ! Anti-Resonant
       if( ABS(coeff(chi%npole_reso_apb+t_ia_global)) > 0.1_dp )  &
         write(stdout,'(8x,i4,a,i4,1x,f12.5)') istate,' <- ',astate,coeff(chi%npole_reso_apb+t_ia_global)
     enddo

     write(stdout,*)

   endif
 enddo

 !
 ! For some calculation conditions, the rest of the subroutine is irrelevant
 ! So skip it! Skip it!
 if( is_triplet .OR. nexc /= chi%npole_reso_apb ) then
   deallocate(residue)
   return
 endif


 !
 ! Calculate the dynamical dipole polarizability
 ! and the static dipole polarizability
 !
 ! Set the frequency mesh
 omega(1)     =MAX( 0.0_dp      ,MINVAL(ABS(eigenvalue(:)))-10.00/Ha_eV)
 omega(nomega)=MIN(50.0_dp/Ha_eV,MAXVAL(ABS(eigenvalue(:)))+10.00/Ha_eV)
 do iomega=2,nomega-1
   omega(iomega) = omega(1) + ( omega(nomega)-omega(1) ) /REAL(nomega-1,dp) * (iomega-1) 
 enddo
 ! Add the broadening
 omega(:) = omega(:) + ieta 

 dynamical_pol(:,:,:) = 0.0_dp
 static_polarizability(:,:) = 0.0_dp
 do t_ia=1,nexc
   do idir=1,3
     do jdir=1,3
       dynamical_pol(:,idir,jdir) = dynamical_pol(:,idir,jdir) &
                            + residue(idir,t_ia) * residue(jdir,t_ia) &
                              * ( AIMAG( -1.0_dp  / ( omega(:) - eigenvalue(t_ia) ) ) - AIMAG( -1.0_dp  / ( omega(:) + eigenvalue(t_ia) ) ) )
       static_polarizability(idir,jdir) = static_polarizability(idir,jdir) &
                      + 2.0_dp * residue(idir,t_ia) * residue(jdir,t_ia) / eigenvalue(t_ia)
     enddo
   enddo
 enddo
 !
 ! Get the photoabsorption cross section
 do iomega=1,nomega
   photoabsorp_cross(iomega,:,:) = 4.0_dp * pi * REAL(omega(iomega),dp) / c_speedlight * dynamical_pol(iomega,:,:)
 enddo

 write(stdout,'(/,a)')     ' TRK sum rule: the two following numbers should compare well'
 write(stdout,'(a,f12.6)') ' Sum over oscillator strengths: ',trk_sumrule
 write(stdout,'(a,f12.6)') '   Number of valence electrons: ',SUM( occupation(ncore_W+1:,:) )

 write(stdout,'(/,a,f12.6)') ' Mean excitation energy (eV): ',EXP( mean_excitation / trk_sumrule ) * Ha_eV

 write(stdout,'(/,a)') ' Static dipole polarizability:'
 trace = 0.0_dp
 do idir=1,3
   write(stdout,'(3(4x,f12.6))') static_polarizability(idir,:)
   trace = trace + static_polarizability(idir,idir) / 3.0_dp
 enddo
 write(stdout,'(a,f12.6)') ' Static dipole polarizability trace: ',trace

 if( is_iomaster ) then

   open(newunit=dynpolfile,file='dynamical_dipole_polarizability.dat',form='formatted')
   open(newunit=photocrossfile,file='photoabsorption_cross_section.dat',form='formatted')
   write(dynpolfile,'(a)') '#  Imaginary part of dynamical dipole polarizability'
   write(dynpolfile,'(a)') '#  omega (eV)   Average     xx    yx    zx    xy    yy    zy    xz    yz    zz'
   write(photocrossfile,'(a)') '#  Imaginary part of dynamical dipole polarizability'
   write(photocrossfile,'(a)') '#  omega (eV)   Average     xx    yx    zx    xy    yy    zy    xz    yz    zz'
   do iomega=1,nomega
     write(dynpolfile,'(11(e18.8,2x))') REAL(omega(iomega),dp)*Ha_eV,                                      &
                                          (dynamical_pol(iomega,1,1)+dynamical_pol(iomega,2,2)+dynamical_pol(iomega,3,3))/3.0_dp, &
                                          dynamical_pol(iomega,:,:)
     write(photocrossfile,'(11(e18.8,2x))') REAL(omega(iomega),dp)*Ha_eV,                                      &
                                              (photoabsorp_cross(iomega,1,1)+photoabsorp_cross(iomega,2,2)+photoabsorp_cross(iomega,3,3))/3.0_dp, &
                                              photoabsorp_cross(iomega,:,:)
   enddo 

   close(dynpolfile)
   close(photocrossfile)

 endif


 deallocate(residue)

 call stop_clock(timing_spectrum)

end subroutine optical_spectrum


!=========================================================================
subroutine stopping_power(nstate,basis,occupation,c_matrix,chi,m_x,n_x,xpy_matrix,eigenvalue)
 use m_definitions
 use m_timing
 use m_warning
 use m_memory
 use m_mpi
 use m_scalapack
 use m_cart_to_pure
 use m_inputparam
 use m_basis_set
 use m_dft_grid
 use m_spectral_function
 implicit none

 integer,intent(in)                 :: nstate,m_x,n_x
 type(basis_set),intent(in)         :: basis
 real(dp),intent(in)                :: occupation(nstate,nspin),c_matrix(basis%nbf,nstate,nspin)
 type(spectral_function),intent(in) :: chi
 real(dp),intent(in)                :: xpy_matrix(m_x,n_x)
 real(dp),intent(in)                :: eigenvalue(chi%npole_reso_apb)
!=====
 integer                            :: gt
 integer                            :: t_ia,t_jb
 integer                            :: t_ia_global,t_jb_global
 integer                            :: nmat
 integer                            :: istate,astate,iaspin
 integer                            :: mpspin
 integer                            :: ishell,jshell
 integer                            :: ibf1,ibf2,jbf1,jbf2,ibf1_cart,jbf1_cart
 integer                            :: ibf,jbf
 integer                            :: ni,nj,li,lj,ni_cart,nj_cart,i_cart,j_cart
 integer                            :: iomega,idir,jdir
 integer,parameter                  :: nomega=600
 complex(dp)                        :: omega(nomega)
 real(dp)                           :: coeff
 real(dp)                           :: dynamical_pol(nomega),structure_factor(nomega)
 complex(dp)                        :: bethe_sumrule
 complex(dp),allocatable            :: gos_basis(:,:),gos_state(:,:,:)
 complex(dp),allocatable            :: gos_cart(:,:)
 complex(dp),allocatable            :: residue(:)
 real(dp)                           :: qvec(3)
 integer,parameter                  :: nq=0 ! 1000
 integer                            :: iq
 real(dp)                           :: fnq(chi%npole_reso_apb)
 integer,parameter                  :: nv=20
 integer                            :: iv
 real(dp)                           :: stopping(nv)
 real(dp)                           :: vv
!=====


 call start_clock(timing_spectrum)
 !
 ! Calculate the spectrum now
 !

 write(stdout,'(/,a)') ' Calculate the stopping power'
 gt = get_gaussian_type_tag(basis%gaussian_type)

 if (nspin/=1) then
   msg='no nspin/=1 allowed'
   call issue_warning(msg)
   return
 endif

 !
 ! Prepare the precalculated table of coefficients
 call setup_gos_llp()


 !
 ! Calculate the dynamical dipole polarizability
 ! and the static dipole polarizability
 !
 ! Set the frequency mesh
 omega(1)     =0.1_dp ! MAX( 0.0_dp      ,MINVAL(ABS(eigenvalue(:)))-3.00/Ha_eV)
 omega(nomega)=4.0_dp ! MIN(20.0_dp/Ha_eV,MAXVAL(ABS(eigenvalue(:)))+3.00/Ha_eV)
 do iomega=2,nomega-1
   omega(iomega) = omega(1) + ( omega(nomega)-omega(1) ) /REAL(nomega-1,dp) * (iomega-1)
 enddo
 ! Add the broadening
 omega(:) = omega(:) + im * 0.10/Ha_eV
  

 
 do iq=1,nq
   qvec(1) = 0.0_dp
   qvec(2) = 0.0_dp
   qvec(3) = iq*0.01_dp

   !
   ! First precalculate all the needed GOS in the basis set
   !
   allocate(gos_basis(basis%nbf,basis%nbf))

   do jshell=1,basis%nshell
     lj      = basis%shell(jshell)%am
     nj_cart = number_basis_function_am('CART',lj)
     nj      = number_basis_function_am(basis%gaussian_type,lj)
     jbf1    = basis%shell(jshell)%istart
     jbf1_cart = basis%shell(jshell)%istart_cart
     jbf2    = basis%shell(jshell)%iend
  
     do ishell=1,basis%nshell
       li      = basis%shell(ishell)%am
       ni_cart = number_basis_function_am('CART',li)
       ni      = number_basis_function_am(basis%gaussian_type,li)
       ibf1    = basis%shell(ishell)%istart
       ibf1_cart = basis%shell(ishell)%istart_cart
       ibf2    = basis%shell(ishell)%iend



       allocate(gos_cart(ni_cart,nj_cart))
  
       do i_cart=1,ni_cart
         do j_cart=1,nj_cart
           call gos_basis_function(basis%bfc(ibf1_cart+i_cart-1),basis%bfc(jbf1_cart+j_cart-1),qvec,gos_cart(i_cart,j_cart))
         enddo
       enddo
  
       gos_basis(ibf1:ibf2,jbf1:jbf2) = MATMUL( TRANSPOSE( cart_to_pure(li,gt)%matrix(:,:) ) , &
                                                MATMUL(  gos_cart(:,:) , cart_to_pure(lj,gt)%matrix(:,:) ) )
  
       deallocate(gos_cart)
  
     enddo
   enddo
  
   !
   ! Get the gos oscillator strength on states
   allocate(gos_state(basis%nbf,basis%nbf,nspin))
  
   do mpspin=1,nspin
     gos_state(:,:,mpspin) = MATMUL( TRANSPOSE( c_matrix(:,:,mpspin) ) ,  MATMUL( gos_basis(:,:) , c_matrix(:,:,mpspin) ) )
   enddo
   deallocate(gos_basis)
  
  
   nmat=chi%npole_reso_apb
   allocate(residue(chi%npole_reso_apb))
  
   residue(:) = 0.0_dp
   do t_ia=1,m_x
     t_ia_global = rowindex_local_to_global(iprow_sd,nprow_sd,t_ia)
     istate = chi%transition_table_apb(1,t_ia_global)
     astate = chi%transition_table_apb(2,t_ia_global)
     iaspin = chi%transition_table_apb(3,t_ia_global)
  
     ! Let use (i <-> j) symmetry to halve the loop
     do t_jb=1,n_x
       t_jb_global = colindex_local_to_global(ipcol_sd,npcol_sd,t_jb)
  
       residue(t_jb_global) = residue(t_jb_global) &
                    + gos_state(istate,astate,iaspin) * xpy_matrix(t_ia,t_jb) * SQRT(spin_fact)
     enddo
  
   enddo
   call xsum_world(residue)
  
   deallocate(gos_state)

   fnq(:) = 2.0_dp * ABS( residue(:) )**2 * eigenvalue(:) / SUM( qvec(:)**2 )

   write(stdout,*) 'bethe_sumrule',NORM2(qvec(:)),SUM(fnq(:))


  
   dynamical_pol(:) = 0.0_dp
   do t_ia=1,nmat
     dynamical_pol(:) = dynamical_pol(:) &
                       + ABS(residue(t_ia))**2 &
                        * ( AIMAG( -1.0_dp  / ( omega(:) - eigenvalue(t_ia) ) ) - AIMAG( -1.0_dp  / ( omega(:) + eigenvalue(t_ia) ) ) )
   enddo
!   !
!   ! Get the structure factor
!   write(999,*) '# qvec',qvec(:)
!   do iomega=1,nomega
!     structure_factor(iomega) = 4.0_dp * pi * REAL(omega(iomega),dp) / c_speedlight * dynamical_pol(iomega) * SUM( qvec(:)**2 )
!     write(999,*) REAL(omega(iomega),dp)*Ha_eV,structure_factor(iomega)
!   enddo
!   write(999,*)


   deallocate(residue)

!   write(998,*) SUM( qvec(:)**2 ), fnq(6)

!   do iv=1,nv
!     vv = iv * 0.1_dp
!     do t_ia=1,nmat
!       if( NORM2(qvec) < eigenvalue(t_ia) / vv )   &
!          stopping(iv) = stopping(iv) + 1.0_dp / ( pi * vv**2 )  * fnq(t_ia)  * NORM2(qvec)**2
!     enddo
!
!   enddo


 enddo 

! do iv=1,nv
!   vv = iv * 0.1_dp
!   write(997,*) vv,stopping(iv)
! enddo



end subroutine stopping_power


!=========================================================================
subroutine get_energy_qp(nstate,energy,occupation,energy_qp)
 use m_definitions
 use m_timing
 use m_warning
 use m_memory
 use m_inputparam
 use m_mpi
 implicit none

 integer,intent(in)                  :: nstate
 real(dp),intent(in)                 :: energy(nstate,nspin)
 real(dp),intent(in)                 :: occupation(nstate,nspin)
 real(dp),intent(out)                :: energy_qp(nstate,nspin)
!=====
 integer  :: reading_status
 real(dp) :: scissor_energy(nspin)
 integer  :: mspin,mstate
!=====

 ! If the keyword scissor is used in the input file,
 ! then use it and ignore the ENERGY_QP file
 if( ABS(scissor) > 1.0e-5_dp ) then

   call issue_warning('Using a manual scissor to open up the fundamental gap')

   write(stdout,'(a,2(1x,f12.6))') ' Scissor operator with value (eV):',scissor*Ha_eV
   do mspin=1,nspin
     do mstate=1,nstate
       if( occupation(mstate,mspin) > completely_empty/spin_fact ) then
         energy_qp(mstate,mspin) = energy(mstate,mspin)
       else
         energy_qp(mstate,mspin) = energy(mstate,mspin) + scissor
       endif
     enddo
   enddo
   write(stdout,'(/,a)') ' Scissor updated energies'
   do mstate=1,nstate
     write(stdout,'(i5,4(2x,f16.6))') mstate,energy(mstate,:)*Ha_eV,energy_qp(mstate,:)*Ha_eV
   enddo
   write(stdout,*)

 else

   call read_energy_qp(nstate,energy_qp,reading_status)

   select case(reading_status)
   case(0)
     write(stdout,'(a)') ' Reading OK'
   case(1,2)
     write(stdout,'(a,/,a)') ' Something happened during the reading of energy_qp file',' Fill up the QP energies with KS energies'
     energy_qp(:,:) = energy(:,:)
   case default
     call die('reading_status BUG')
   end select

 endif

end subroutine get_energy_qp


!=========================================================================
subroutine chi_to_vchiv(nbf,nstate,c_matrix,xpy_matrix,eigenvalue,wpol)
! ymbyun 2021/12/27
#ifdef _OPENACC
#if _OPENACC >= 201711
 use openacc
#endif
#endif
! ymbyun 2023/08/30
#ifdef ENABLE_TIME
#ifdef _OPENMP
 use omp_lib
#endif
#endif

 use m_definitions
 use m_warning
 use m_basis_set
 use m_eri_ao_mo
 use m_spectral_function
 implicit none
 
 integer,intent(in)                    :: nbf,nstate
 real(dp),intent(in)                   :: c_matrix(nbf,nstate,nspin)
 type(spectral_function),intent(inout) :: wpol
 real(dp),intent(in)                   :: xpy_matrix(wpol%npole_reso_apb,wpol%npole_reso_apb)
 real(dp),intent(in)                   :: eigenvalue(wpol%npole_reso_apb)
!=====
 integer                               :: t_jb,jbspin,mpspin
 integer                               :: mstate,pstate,jstate,bstate,mpstate_spin
 integer                               :: kbstate_min
 integer                               :: kbstate_max
 integer                               :: nmat,nprodbasis
 real(dp)                              :: eri_eigen_klij
 real(dp),allocatable                  :: eri_eigenstate_klmin(:,:,:,:)
! ymbyun 2018/05/23
#if defined(_OPENMP) || defined(_OPENACC)
 integer                               :: spstst2in(nspin,nstate,nstate),temp_in
#endif
! ymbyun 2021/12/27
#ifdef _OPENACC
#if _OPENACC >= 201711
 integer                               :: device_num
 integer(acc_device_kind)              :: device_type
 integer(kind=8)                       :: gpu_memory   ,gpu_free_memory
 integer                               :: gpu_memory_mb,gpu_free_memory_mb
#endif
#endif
! ymbyun 2023/08/30
! NOTE: It seems like omp_get_wtime() is more stable than cpu_time() and system_clock() for large calculations.
! NOTE: Time and clock routines are compiler specific.
#ifdef ENABLE_TIME
 real(dp)                              :: time_all_start,time_all_end
 real(dp)                              :: time_eri_start,time_eri_end
 real(dp)                              :: time_acc_start,time_acc_end
 integer(kind=8)                       :: clock_all_start,clock_all_end,clock_all_rate
 integer(kind=8)                       :: clock_eri_start,clock_eri_end,clock_eri_rate
 integer(kind=8)                       :: clock_acc_start,clock_acc_end,clock_acc_rate
#ifdef _OPENMP
 real(dp)                              :: omp_all_start,omp_all_end
 real(dp)                              :: omp_eri_start,omp_eri_end
 real(dp)                              :: omp_acc_start,omp_acc_end
#endif
#endif
!=====

 call start_clock(timing_vchiv)

! ymbyun 2023/08/30
#ifdef ENABLE_TIME
 call cpu_time(time_all_start)
 call system_clock(count_rate=clock_all_rate)
 call system_clock(count=clock_all_start)
#ifdef _OPENMP
 omp_all_start = omp_get_wtime()
#endif
#endif

 write(stdout,'(/,a)') ' Build W = v * chi * v'
 if(has_auxil_basis) then
   call die('you should not be here')
 endif

 allocate(eri_eigenstate_klmin(nbf,nbf,nbf,nspin))
 ! Set this to zero and then enforce the calculation of the first array of Coulomb integrals
 ! ymbyun 2018/05/30
 ! This array is initialized (i.e. first touched) in calculate_eri_4center_eigen().
 ! ymbyun 2018/11/12
 ! No, it shouldn't be.
 eri_eigenstate_klmin(:,:,:,:) = 0.0_dp

 nprodbasis = index_prodstate(nvirtual_W-1,nvirtual_W-1) * nspin
 call allocate_spectral_function(nprodbasis,wpol)
 wpol%pole(1:wpol%npole_reso_apb) = eigenvalue(:)
 nmat = wpol%npole_reso_apb
 wpol%residue_left(:,:) = 0.0_dp

! ymbyun 2018/05/23
! ymbyun 2021/12/06
#if defined(_OPENMP) || defined(_OPENACC)
 temp_in = 0

 do mpspin=1,nspin
   do pstate=1,nstate
     do mstate=1,pstate
       temp_in = temp_in + 1
       spstst2in(mpspin,pstate,mstate) = temp_in
     enddo
   enddo
 enddo
#endif

!!$acc data copyin(xpy_matrix,eri_eigenstate_klmin,spstst2in)
!!$acc enter data copyin(wpol,wpol%residue_left)
 do t_jb=1,nmat
   jstate = wpol%transition_table_apb(1,t_jb)
   bstate = wpol%transition_table_apb(2,t_jb)
   jbspin = wpol%transition_table_apb(3,t_jb)

   kbstate_min = MIN(jstate,bstate)
   kbstate_max = MAX(jstate,bstate)

   ! ymbyun 2018/05/22
   ! NOTE: calculate_eri_4center_eigen() is OpenMP parallelized without being orphaned.

! ymbyun 2023/08/30
#ifdef ENABLE_TIME
   call cpu_time(time_eri_start)
   call system_clock(count_rate=clock_eri_rate)
   call system_clock(count=clock_eri_start)
#ifdef _OPENMP
   omp_eri_start = omp_get_wtime()
#endif
#endif

   call calculate_eri_4center_eigen(nbf,nstate,c_matrix,kbstate_min,jbspin,eri_eigenstate_klmin)

! ymbyun 2023/08/30
#ifdef ENABLE_TIME
   call cpu_time(time_eri_end)
   if( (time_eri_end - time_eri_start) > 0.000001_dp ) then
     print '("(ymbyun) chi_to_vchiv (time_eri)  : ",f30.6," sec")',time_eri_end-time_eri_start
   endif

   call system_clock(count=clock_eri_end)
   if( ((clock_eri_end-clock_eri_start)/real(clock_eri_rate)) > 0.000001_dp ) then
     print '("(ymbyun) chi_to_vchiv (clock_eri) : ",f30.6," sec")',(clock_eri_end-clock_eri_start)/real(clock_eri_rate)
   endif

#ifdef _OPENMP
   omp_eri_end = omp_get_wtime()
   if( (omp_eri_end - omp_eri_start) > 0.000001_dp ) then
     print '("(ymbyun) chi_to_vchiv (omp_eri)   : ",f30.6," sec")',omp_eri_end-omp_eri_start
   endif
#endif
#endif

   mpstate_spin = 0

! ymbyun 2018/05/22
! NOTE: COLLAPSE is used because nspin is much smaller than # of threads.
! NOTE: mpstate_spin should be used carefully.
! NOTE: COLLAPSE(2) can't be replaced by COLLAPSE(3) because mstate depends on pstate.
! TODO: PARALLEL should be moved from the inner loop to the outer one to reduce the overhead.

! ymbyun 2021/12/06
! OpenACC is added.
! NOTE: COPYIN/UPDATE and ENTER COPYIN/EXIT COPYOUT give the same result for performance.
! NOTE: Manual deep copy is used, because OpenACC 2.6 is used. Another option is true deep copy (-ta=tesla:deepcopy or OpenACC 3.0).
! NOTE: DATA here outperforms DATA at the outermost loop.
! NOTE: CUDA unified memory (-ta=tesla:managed) outperforms manual deep copy.

! ymbyun 2023/08/30
#ifdef ENABLE_TIME
   call cpu_time(time_acc_start)
   call system_clock(count_rate=clock_acc_rate)
   call system_clock(count=clock_acc_start)
#ifdef _OPENMP
   omp_acc_start = omp_get_wtime()
#endif
#endif

#if defined (_OPENMP) && !defined (_OPENACC)
 !$OMP PARALLEL
 !$OMP DO PRIVATE(mpspin,mstate,pstate,eri_eigen_klij,mpstate_spin) COLLAPSE(2)
#endif
#if !defined (_OPENMP) && defined (_OPENACC)
!!$acc data copyin(xpy_matrix,eri_eigenstate_klmin,spstst2in) copyin(wpol,wpol%residue_left)
 !$acc data copyin(xpy_matrix,eri_eigenstate_klmin,spstst2in)
 !$acc enter data copyin(wpol,wpol%residue_left)

#if _OPENACC >= 201711
   if( t_jb == 1 ) then
     device_type        = acc_get_device_type()
     device_num         = acc_get_device_num(device_type)
  
     gpu_memory         = acc_get_property(device_num,device_type,acc_property_memory     )
     gpu_free_memory    = acc_get_property(device_num,device_type,acc_property_free_memory)
  
     gpu_memory_mb      = gpu_memory      / 1024 / 1024
     gpu_free_memory_mb = gpu_free_memory / 1024 / 1024
  
     write(stdout,*) 'OpenACC total/free/used memory (MB):',gpu_memory_mb,gpu_free_memory_mb,(gpu_memory_mb-gpu_free_memory_mb)
  
     if (gpu_free_memory_mb < 1024) then
       call issue_warning('(ymbyun) Less than 1 GB of GPU memory is free now. Maybe this job is using too much GPU memory.')
     endif
   endif
#endif

 !$acc parallel loop independent collapse(2)
#endif

   do mpspin=1,nspin
     do pstate=1,nstate
       do mstate=1,pstate
! ymbyun 2018/05/23
! ymbyun 2021/12/06
#if defined(_OPENMP) || defined(_OPENACC)
         mpstate_spin = spstst2in(mpspin,pstate,mstate)
#else
         mpstate_spin = mpstate_spin + 1
#endif
         eri_eigen_klij = eri_eigenstate_klmin(kbstate_max,mstate,pstate,mpspin)

         ! Use the symmetry ( k l | i j ) to regroup (kl) and (lk) contributions
         ! and the block structure of eigenvector | X  Y |
         !                                        | Y  X |
         wpol%residue_left(mpstate_spin,:) = wpol%residue_left(mpstate_spin,:) &
                              + eri_eigen_klij * xpy_matrix(t_jb,:)
       enddo
     enddo
   enddo
#if !defined (_OPENMP) && defined (_OPENACC)
 !$acc end parallel
!!$acc update self(wpol%residue_left)
 !$acc exit data copyout(wpol%residue_left,wpol)
 !$acc end data
#endif
#if defined (_OPENMP) && !defined (_OPENACC)
 !$OMP END DO
 !$OMP END PARALLEL
#endif

! ymbyun 2023/08/30
#ifdef ENABLE_TIME
   call cpu_time(time_acc_end)
   print '("(ymbyun) chi_to_vchiv (time_acc)  : ",f30.6," sec")',time_acc_end-time_acc_start

   call system_clock(count=clock_acc_end)
   print '("(ymbyun) chi_to_vchiv (clock_acc) : ",f30.6," sec")',(clock_acc_end-clock_acc_start)/real(clock_acc_rate)

#ifdef _OPENMP
   omp_acc_end = omp_get_wtime()
   print '("(ymbyun) chi_to_vchiv (omp_acc)   : ",f30.6," sec")',omp_acc_end-omp_acc_start
#endif
#endif

 enddo
!!$acc exit data copyout(wpol%residue_left,wpol)
!!$acc end data

! ymbyun 2018/08/03
! NOTE: A performance test is needed.
!$OMP PARALLEL
!$OMP WORKSHARE
 wpol%residue_left(:,:) = wpol%residue_left(:,:) * SQRT(spin_fact)
!$OMP END WORKSHARE
!$OMP END PARALLEL

 if(ALLOCATED(eri_eigenstate_klmin)) deallocate(eri_eigenstate_klmin)

 call stop_clock(timing_vchiv)

! ymbyun 2023/08/30
#ifdef ENABLE_TIME
 call cpu_time(time_all_end)
 print '("(ymbyun) chi_to_vchiv (time_all)  : ",f30.6," sec")',time_all_end-time_all_start

 call system_clock(count=clock_all_end)
 print '("(ymbyun) chi_to_vchiv (clock_all) : ",f30.6," sec")',(clock_all_end-clock_all_start)/real(clock_all_rate)

#ifdef _OPENMP
 omp_all_end = omp_get_wtime()
 print '("(ymbyun) chi_to_vchiv (omp_all)   : ",f30.6," sec")',omp_all_end-omp_all_start
#endif
#endif

end subroutine chi_to_vchiv


!=========================================================================
subroutine chi_to_sqrtvchisqrtv_auxil(nbf,desc_x,m_x,n_x,xpy_matrix,eigenvalue,wpol,energy_gm)
 use m_definitions
 use m_warning
 use m_scalapack
 use m_basis_set
 use m_eri_ao_mo
 use m_spectral_function
 implicit none
 
 integer,intent(in)                    :: nbf,m_x,n_x
 integer,intent(in)                    :: desc_x(NDEL)
 real(dp),intent(inout)                :: xpy_matrix(m_x,n_x)
 type(spectral_function),intent(inout) :: wpol
 real(dp),intent(in)                   :: eigenvalue(wpol%npole_reso_apb)
 real(dp),intent(out)                  :: energy_gm
!=====
 integer                               :: t_jb,t_jb_global,jbspin
 integer                               :: nmat
 integer                               :: jstate,bstate
 real(dp),allocatable                  :: eri_3tmp(:,:)
 real(dp),allocatable                  :: eri_3tmp_sd(:,:)
 real(dp),allocatable                  :: xpy_block(:,:)
 real(dp),allocatable                  :: vsqrt_xpy(:,:)
 integer                               :: desc_auxil(NDEL),desc_sd(NDEL)
 integer                               :: mlocal,nlocal
 integer                               :: info
!=====

 call start_clock(timing_vchiv)

 write(stdout,'(/,a)') ' Build v^{1/2} * chi * v^{1/2}'

 call allocate_spectral_function(nauxil_3center,wpol)
 wpol%pole(1:wpol%npole_reso_apb) = eigenvalue(:)

 nmat = wpol%npole_reso_apb

#ifndef HAVE_SCALAPACK

 allocate(eri_3tmp(nauxil_3center,nmat))
 do t_jb=1,nmat
   jstate = wpol%transition_table_apb(1,t_jb)
   bstate = wpol%transition_table_apb(2,t_jb)
   jbspin = wpol%transition_table_apb(3,t_jb)
   eri_3tmp(:,t_jb) = eri_3center_eigen(:,jstate,bstate,jbspin)
 enddo

 ! Use the symmetry ( I | k l ) to regroup (kl) and (lk) contributions
 ! and the block structure of eigenvector | X  Y |
 !                                        | Y  X |
 ! => only needs (X+Y)
 wpol%residue_left(:,:) = MATMUL( eri_3tmp , xpy_matrix(:,:) ) * SQRT(spin_fact)

 energy_gm = 0.5_dp * ( SUM( wpol%residue_left(:,:)**2 ) - spin_fact * SUM( eri_3tmp(:,:)**2 ) )
 !
 ! Since wpol%residue_left and eri_3tmp are distributed, we have to sum up
 call xsum_auxil(energy_gm)

 deallocate(eri_3tmp)

#else 

 call clean_allocate('TMP 3-center integrals',eri_3tmp,nauxil_3center,nmat)
 do t_jb=1,nmat
   jstate = wpol%transition_table_apb(1,t_jb)
   bstate = wpol%transition_table_apb(2,t_jb)
   jbspin = wpol%transition_table_apb(3,t_jb)
   eri_3tmp(:,t_jb) = eri_3center_eigen(:,jstate,bstate,jbspin)
 enddo

 !
 ! Descriptors
 mlocal = NUMROC(nauxil_2center,MBLOCK_AUXIL,iprow_auxil,first_row,nprow_auxil)
 call DESCINIT(desc_auxil,nauxil_2center,nmat,MBLOCK_AUXIL,NBLOCK_AUXIL,first_row,first_col,cntxt_auxil,MAX(1,mlocal),info)

 mlocal = NUMROC(nauxil_2center,block_row,iprow_sd,first_row,nprow_sd)
 nlocal = NUMROC(nmat          ,block_col,ipcol_sd,first_col,npcol_sd)
 call DESCINIT(desc_sd,nauxil_2center,nmat,block_row,block_col,first_row,first_col,cntxt_sd,MAX(1,mlocal),info)
 call clean_allocate('TMP 3-center integrals',eri_3tmp_sd,mlocal,nlocal)

 call PDGEMR2D(nauxil_2center,nmat,eri_3tmp,1,1,desc_auxil, &
                                eri_3tmp_sd,1,1,desc_sd,cntxt_sd)

 call clean_deallocate('TMP 3-center integrals',eri_3tmp)
 
 !
 !   SQRT(spin_fact) * v**1/2 * ( X + Y )
 call clean_allocate('TMP v**1/2 * (X+Y)',vsqrt_xpy,mlocal,nlocal)
 call PDGEMM('N','N',nauxil_2center,nmat,nmat, &
               DSQRT(spin_fact),eri_3tmp_sd,1,1,desc_sd,  &
                                 xpy_matrix,1,1,desc_x,   &
                      0.0_dp,     vsqrt_xpy,1,1,desc_sd)

 call clean_deallocate('TMP 3-center integrals',eri_3tmp_sd)


 call PDGEMR2D(nauxil_2center,nmat,vsqrt_xpy,1,1,desc_sd, &
                            wpol%residue_left,1,1,desc_auxil,cntxt_sd)
 !
 ! Do not forget ortho parallelization direction
 if( nproc_ortho > 1 ) then
   call xbcast_ortho(0,wpol%residue_left)
 endif

 call clean_deallocate('TMP v**1/2 * (X+Y)',vsqrt_xpy)


 energy_gm = 0.0_dp
 do t_jb_global=1,nmat
   jstate = wpol%transition_table_apb(1,t_jb_global)
   bstate = wpol%transition_table_apb(2,t_jb_global)
   jbspin = wpol%transition_table_apb(3,t_jb_global)
   energy_gm = energy_gm - SUM( eri_3center_eigen(:,jstate,bstate,jbspin)**2 ) * spin_fact * 0.5_dp
 enddo

 energy_gm = energy_gm + 0.5_dp * ( SUM( wpol%residue_left(:,:)**2 ) )
 call xsum_auxil(energy_gm)


#endif



 call stop_clock(timing_vchiv)

end subroutine chi_to_sqrtvchisqrtv_auxil


!=========================================================================
subroutine chi_to_sqrtvchisqrtv_auxil_spa(nbf,a_diag,wpol)
 use m_definitions
 use m_warning
 use m_basis_set
 use m_eri_ao_mo
 use m_spectral_function
 implicit none
 
 integer,intent(in)                    :: nbf
 type(spectral_function),intent(inout) :: wpol
 real(dp),intent(in)                   :: a_diag(wpol%npole_reso_spa)
!=====
 integer                               :: t_jb,jbspin
 integer                               :: jstate,bstate
!=====

 call start_clock(timing_vchiv)

 write(stdout,'(/,a)') ' Build v^{1/2} * chi * v^{1/2} part from single pole approximation'

 wpol%pole(wpol%npole_reso_apb+1:wpol%npole_reso) = a_diag(:)

 do t_jb=1,wpol%npole_reso_spa
   jstate = wpol%transition_table_spa(1,t_jb)
   bstate = wpol%transition_table_spa(2,t_jb)
   jbspin = wpol%transition_table_spa(3,t_jb)


   wpol%residue_left(:,wpol%npole_reso_apb+t_jb) = eri_3center_eigen(:,jstate,bstate,jbspin) * SQRT(spin_fact) 

 enddo


 call stop_clock(timing_vchiv)

end subroutine chi_to_sqrtvchisqrtv_auxil_spa


!=========================================================================
