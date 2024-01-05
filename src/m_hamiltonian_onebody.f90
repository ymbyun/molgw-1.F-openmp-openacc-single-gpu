!=========================================================================
! This file is part of MOLGW.
! Authors: Fabien Bruneval, Ivan Maliyov
!
! This module contains
! the methods to evaluate the Kohn-Sham Hamiltonian
! with no distribution of the memory
!
!=========================================================================
module m_hamiltonian_onebody
 use m_definitions
 use m_timing
 use m_mpi
 use m_scalapack
 use m_warning
 use m_memory
 use m_cart_to_pure
 use m_inputparam,only: nspin,spin_fact,scalapack_block_min
 use m_basis_set
 use m_libint_tools





contains


!=========================================================================
subroutine setup_overlap(print_matrix_,basis,s_matrix)
 implicit none
 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: s_matrix(basis%nbf,basis%nbf)
!=====
 integer              :: ishell,jshell
 integer              :: ibf1,ibf2,jbf1,jbf2
 integer              :: ibf_cart,jbf_cart
 integer              :: i_cart,j_cart,ij
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 character(len=100)   :: title
 real(dp),allocatable :: matrix(:,:)
 real(dp)             :: overlap

 real(C_DOUBLE),allocatable :: array_cart(:)
 integer(C_INT)             :: amA,contrdepthA
 real(C_DOUBLE)             :: A(3)
 real(C_DOUBLE),allocatable :: alphaA(:)
 real(C_DOUBLE),allocatable :: cA(:)
 integer(C_INT)             :: amB,contrdepthB
 real(C_DOUBLE)             :: B(3)
 real(C_DOUBLE),allocatable :: alphaB(:)
 real(C_DOUBLE),allocatable :: cB(:)
!=====

 call start_clock(timing_overlap)
 write(stdout,'(/,a)') ' Setup overlap matrix S (LIBINT)'

 do jshell=1,basis%nshell
   lj      = basis%shell(jshell)%am
   nj_cart = number_basis_function_am('CART',lj)
   nj      = number_basis_function_am(basis%gaussian_type,lj)
   jbf1    = basis%shell(jshell)%istart
   jbf2    = basis%shell(jshell)%iend

   call set_libint_shell(basis%shell(jshell),amB,contrdepthB,B,alphaB,cB)

   do ishell=jshell,basis%nshell
     li      = basis%shell(ishell)%am
     ni_cart = number_basis_function_am('CART',li)
     ni      = number_basis_function_am(basis%gaussian_type,li)
     ibf1    = basis%shell(ishell)%istart
     ibf2    = basis%shell(ishell)%iend

     call set_libint_shell(basis%shell(ishell),amA,contrdepthA,A,alphaA,cA)


     allocate(array_cart(ni_cart*nj_cart))

#ifdef HAVE_LIBINT_ONEBODY

     call libint_overlap(amA,contrdepthA,A,alphaA,cA, &
                         amB,contrdepthB,B,alphaB,cB, &
                         array_cart)

     call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart,matrix)
#else
     ij = 0
     do i_cart=1,ni_cart
       do j_cart=1,nj_cart
         ij = ij + 1
         ibf_cart = basis%shell(ishell)%istart_cart + i_cart - 1
         jbf_cart = basis%shell(jshell)%istart_cart + j_cart - 1
         call overlap_basis_function(basis%bfc(ibf_cart),basis%bfc(jbf_cart),overlap)
         array_cart(ij) = overlap
       enddo
     enddo
     call transform_molgw_to_molgw(basis%gaussian_type,li,lj,array_cart,matrix)
#endif

     deallocate(alphaA,cA)


     s_matrix(ibf1:ibf2,jbf1:jbf2) = matrix(:,:)
     s_matrix(jbf1:jbf2,ibf1:ibf2) = TRANSPOSE(matrix(:,:))

     deallocate(array_cart,matrix)

   enddo
   deallocate(alphaB,cB)
 enddo

 title='=== Overlap matrix S (LIBINT) ==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,s_matrix)

 
 call stop_clock(timing_overlap)


end subroutine setup_overlap


!=========================================================================
subroutine setup_overlap_mixedbasis(print_matrix_,basis1,basis2,s_matrix)
 implicit none
 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis1,basis2
 real(dp),intent(out)       :: s_matrix(basis1%nbf,basis2%nbf)
!=====
 integer              :: ishell,jshell
 integer              :: ibf1,ibf2,jbf1,jbf2
 integer              :: ibf_cart,jbf_cart
 integer              :: i_cart,j_cart,ij
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 character(len=100)   :: title
 real(dp),allocatable :: matrix(:,:)
 real(dp)             :: overlap

 real(C_DOUBLE),allocatable :: array_cart(:)
 integer(C_INT)             :: amA,contrdepthA
 real(C_DOUBLE)             :: A(3)
 real(C_DOUBLE),allocatable :: alphaA(:)
 real(C_DOUBLE),allocatable :: cA(:)
 integer(C_INT)             :: amB,contrdepthB
 real(C_DOUBLE)             :: B(3)
 real(C_DOUBLE),allocatable :: alphaB(:)
 real(C_DOUBLE),allocatable :: cB(:)
!=====

 call start_clock(timing_overlap)
 write(stdout,'(/,a)') ' Setup mixed overlap matrix S (LIBINT)'

 if( basis1%gaussian_type /= basis2%gaussian_type ) call die('setup_overlap_mixedbasis_libint: case not implemented')

 do jshell=1,basis2%nshell
   lj      = basis2%shell(jshell)%am
   nj_cart = number_basis_function_am('CART',lj)
   nj      = number_basis_function_am(basis2%gaussian_type,lj)
   jbf1    = basis2%shell(jshell)%istart
   jbf2    = basis2%shell(jshell)%iend

   call set_libint_shell(basis2%shell(jshell),amB,contrdepthB,B,alphaB,cB)

   do ishell=1,basis1%nshell
     li      = basis1%shell(ishell)%am
     ni_cart = number_basis_function_am('CART',li)
     ni      = number_basis_function_am(basis1%gaussian_type,li)
     ibf1    = basis1%shell(ishell)%istart
     ibf2    = basis1%shell(ishell)%iend

     call set_libint_shell(basis1%shell(ishell),amA,contrdepthA,A,alphaA,cA)


     allocate(array_cart(ni_cart*nj_cart))

#ifdef HAVE_LIBINT_ONEBODY

     call libint_overlap(amA,contrdepthA,A,alphaA,cA, &
                         amB,contrdepthB,B,alphaB,cB, &
                         array_cart)

     call transform_libint_to_molgw(basis1%gaussian_type,li,lj,array_cart,matrix)
#else
     ij = 0
     do i_cart=1,ni_cart
       do j_cart=1,nj_cart
         ij = ij + 1
         ibf_cart = basis1%shell(ishell)%istart_cart + i_cart - 1
         jbf_cart = basis2%shell(jshell)%istart_cart + j_cart - 1
         call overlap_basis_function(basis1%bfc(ibf_cart),basis2%bfc(jbf_cart),overlap)
         array_cart(ij) = overlap
       enddo
     enddo
     call transform_molgw_to_molgw(basis1%gaussian_type,li,lj,array_cart,matrix)
#endif

     deallocate(alphaA,cA)


     s_matrix(ibf1:ibf2,jbf1:jbf2) = matrix(:,:)

     deallocate(array_cart,matrix)

   enddo
   deallocate(alphaB,cB)
 enddo

 
 call stop_clock(timing_overlap)


end subroutine setup_overlap_mixedbasis


!=========================================================================
subroutine setup_overlap_grad(print_matrix_,basis,s_matrix_grad)
 implicit none
 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: s_matrix_grad(basis%nbf,basis%nbf,3)
!=====
 integer              :: ishell,jshell
 integer              :: ibf1,ibf2,jbf1,jbf2
 integer              :: i_cart,j_cart
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 character(len=100)   :: title
 real(dp),allocatable :: matrix(:,:)

 real(C_DOUBLE),allocatable        :: array_cart_gradx(:)
 real(C_DOUBLE),allocatable        :: array_cart_grady(:)
 real(C_DOUBLE),allocatable        :: array_cart_gradz(:)
 integer(C_INT)                    :: amA,contrdepthA
 real(C_DOUBLE)                    :: A(3)
 real(C_DOUBLE),allocatable        :: alphaA(:)
 real(C_DOUBLE),allocatable        :: cA(:)
 integer(C_INT)                    :: amB,contrdepthB
 real(C_DOUBLE)                    :: B(3)
 real(C_DOUBLE),allocatable        :: alphaB(:)
 real(C_DOUBLE),allocatable        :: cB(:)
!=====

 call start_clock(timing_overlap)
 write(stdout,'(/,a)') ' Setup gradient of the overlap matrix S (LIBINT)'

 do jshell=1,basis%nshell
   lj      = basis%shell(jshell)%am
   nj_cart = number_basis_function_am('CART',lj)
   nj      = number_basis_function_am(basis%gaussian_type,lj)
   jbf1    = basis%shell(jshell)%istart
   jbf2    = basis%shell(jshell)%iend

   call set_libint_shell(basis%shell(jshell),amB,contrdepthB,B,alphaB,cB)

   do ishell=1,basis%nshell
     li      = basis%shell(ishell)%am
     ni_cart = number_basis_function_am('CART',li)
     ni      = number_basis_function_am(basis%gaussian_type,li)
     ibf1    = basis%shell(ishell)%istart
     ibf2    = basis%shell(ishell)%iend

     call set_libint_shell(basis%shell(ishell),amA,contrdepthA,A,alphaA,cA)

     allocate(array_cart_gradx(ni_cart*nj_cart))
     allocate(array_cart_grady(ni_cart*nj_cart))
     allocate(array_cart_gradz(ni_cart*nj_cart))

#ifdef HAVE_LIBINT_ONEBODY
     call libint_overlap_grad(amA,contrdepthA,A,alphaA,cA, &
                              amB,contrdepthB,B,alphaB,cB, &
                              array_cart_gradx,array_cart_grady,array_cart_gradz)
#else
     call die('overlap gradient not implemented without LIBINT one-body terms')
#endif

     deallocate(alphaA,cA)

     ! X
     call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradx,matrix)
     s_matrix_grad(ibf1:ibf2,jbf1:jbf2,1) = matrix(:,:)

     ! Y
     call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_grady,matrix)
     s_matrix_grad(ibf1:ibf2,jbf1:jbf2,2) = matrix(:,:)

     ! Z
     call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradz,matrix)
     s_matrix_grad(ibf1:ibf2,jbf1:jbf2,3) = matrix(:,:)


     deallocate(array_cart_gradx)
     deallocate(array_cart_grady)
     deallocate(array_cart_gradz)
     deallocate(matrix)

   enddo
   deallocate(alphaB,cB)
 enddo

 title='=== Overlap matrix S (LIBINT) X ==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,s_matrix_grad(:,:,1))
 title='=== Overlap matrix S (LIBINT) Y ==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,s_matrix_grad(:,:,2))
 title='=== Overlap matrix S (LIBINT) Z ==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,s_matrix_grad(:,:,3))

 call stop_clock(timing_overlap)


end subroutine setup_overlap_grad


!=========================================================================
subroutine setup_kinetic(print_matrix_,basis,hamiltonian_kinetic)
 implicit none
 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: hamiltonian_kinetic(basis%nbf,basis%nbf)
!=====
 integer              :: ishell,jshell
 integer              :: ibf1,ibf2,jbf1,jbf2
 integer              :: ibf_cart,jbf_cart
 integer              :: i_cart,j_cart,ij
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 character(len=100)   :: title
 real(dp),allocatable :: matrix(:,:)
 real(dp)             :: kinetic

 real(C_DOUBLE),allocatable :: array_cart(:)
 integer(C_INT)             :: amA,contrdepthA
 real(C_DOUBLE)             :: A(3)
 real(C_DOUBLE),allocatable :: alphaA(:)
 real(C_DOUBLE),allocatable :: cA(:)
 integer(C_INT)             :: amB,contrdepthB
 real(C_DOUBLE)             :: B(3)
 real(C_DOUBLE),allocatable :: alphaB(:)
 real(C_DOUBLE),allocatable :: cB(:)
!=====

 call start_clock(timing_hamiltonian_kin)
 write(stdout,'(/,a)') ' Setup kinetic part of the Hamiltonian (LIBINT)'


 do jshell=1,basis%nshell
   lj      = basis%shell(jshell)%am
   nj_cart = number_basis_function_am('CART',lj)
   nj      = number_basis_function_am(basis%gaussian_type,lj)
   jbf1    = basis%shell(jshell)%istart
   jbf2    = basis%shell(jshell)%iend

   call set_libint_shell(basis%shell(jshell),amB,contrdepthB,B,alphaB,cB)

   do ishell=jshell,basis%nshell
     li      = basis%shell(ishell)%am
     ni_cart = number_basis_function_am('CART',li)
     ni      = number_basis_function_am(basis%gaussian_type,li)
     ibf1    = basis%shell(ishell)%istart
     ibf2    = basis%shell(ishell)%iend

     call set_libint_shell(basis%shell(ishell),amA,contrdepthA,A,alphaA,cA)


     allocate(array_cart(ni_cart*nj_cart))


#ifdef HAVE_LIBINT_ONEBODY
     call libint_kinetic(amA,contrdepthA,A,alphaA,cA, &
                         amB,contrdepthB,B,alphaB,cB, &
                         array_cart)
     call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart,matrix)
#else
     ij = 0
     do i_cart=1,ni_cart
       do j_cart=1,nj_cart
         ij = ij + 1
         ibf_cart = basis%shell(ishell)%istart_cart + i_cart - 1
         jbf_cart = basis%shell(jshell)%istart_cart + j_cart - 1
         call kinetic_basis_function(basis%bfc(ibf_cart),basis%bfc(jbf_cart),kinetic)
         array_cart(ij) = kinetic
       enddo
     enddo
     call transform_molgw_to_molgw(basis%gaussian_type,li,lj,array_cart,matrix)
#endif
     deallocate(alphaA,cA)



     hamiltonian_kinetic(ibf1:ibf2,jbf1:jbf2) = matrix(:,:)
     hamiltonian_kinetic(jbf1:jbf2,ibf1:ibf2) = TRANSPOSE(matrix(:,:))


     deallocate(array_cart,matrix)

   enddo
   deallocate(alphaB,cB)
 enddo

 title='===  Kinetic energy contribution (LIBINT) ==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,hamiltonian_kinetic)

 call stop_clock(timing_hamiltonian_kin)

end subroutine setup_kinetic


!=========================================================================
subroutine setup_kinetic_grad(print_matrix_,basis,hamiltonian_kinetic_grad)
 implicit none
 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: hamiltonian_kinetic_grad(basis%nbf,basis%nbf,3)
!=====
 integer              :: ishell,jshell
 integer              :: ibf1,ibf2,jbf1,jbf2
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 character(len=100)   :: title
 real(dp),allocatable :: matrix(:,:)

 real(C_DOUBLE),allocatable        :: array_cart_gradx(:)
 real(C_DOUBLE),allocatable        :: array_cart_grady(:)
 real(C_DOUBLE),allocatable        :: array_cart_gradz(:)
 integer(C_INT)                    :: amA,contrdepthA
 real(C_DOUBLE)                    :: A(3)
 real(C_DOUBLE),allocatable        :: alphaA(:)
 real(C_DOUBLE),allocatable        :: cA(:)
 integer(C_INT)                    :: amB,contrdepthB
 real(C_DOUBLE)                    :: B(3)
 real(C_DOUBLE),allocatable        :: alphaB(:)
 real(C_DOUBLE),allocatable        :: cB(:)
!=====

 call start_clock(timing_hamiltonian_kin)
 write(stdout,'(/,a)') ' Setup gradient of the kinetic part of the Hamiltonian (LIBINT)'

 do jshell=1,basis%nshell
   lj      = basis%shell(jshell)%am
   nj_cart = number_basis_function_am('CART',lj)
   nj      = number_basis_function_am(basis%gaussian_type,lj)
   jbf1    = basis%shell(jshell)%istart
   jbf2    = basis%shell(jshell)%iend

   call set_libint_shell(basis%shell(jshell),amB,contrdepthB,B,alphaB,cB)

   do ishell=1,basis%nshell
     li      = basis%shell(ishell)%am
     ni_cart = number_basis_function_am('CART',li)
     ni      = number_basis_function_am(basis%gaussian_type,li)
     ibf1    = basis%shell(ishell)%istart
     ibf2    = basis%shell(ishell)%iend

     call set_libint_shell(basis%shell(ishell),amA,contrdepthA,A,alphaA,cA)


     allocate(array_cart_gradx(ni_cart*nj_cart))
     allocate(array_cart_grady(ni_cart*nj_cart))
     allocate(array_cart_gradz(ni_cart*nj_cart))

#ifdef HAVE_LIBINT_ONEBODY
     call libint_kinetic_grad(amA,contrdepthA,A,alphaA,cA, &
                              amB,contrdepthB,B,alphaB,cB, &
                              array_cart_gradx,array_cart_grady,array_cart_gradz)
#else
     call die('kinetic operator gradient not implemented without LIBINT one-body terms')
#endif
     deallocate(alphaA,cA)

     ! X
     call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradx,matrix)
     hamiltonian_kinetic_grad(ibf1:ibf2,jbf1:jbf2,1) = matrix(:,:)

     ! Y
     call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_grady,matrix)
     hamiltonian_kinetic_grad(ibf1:ibf2,jbf1:jbf2,2) = matrix(:,:)

     ! Z
     call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradz,matrix)
     hamiltonian_kinetic_grad(ibf1:ibf2,jbf1:jbf2,3) = matrix(:,:)


     deallocate(array_cart_gradx)
     deallocate(array_cart_grady)
     deallocate(array_cart_gradz)
     deallocate(matrix)

   enddo
   deallocate(alphaB,cB)
 enddo

 title='===  Kinetic energy contribution (LIBINT) X ==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,hamiltonian_kinetic_grad(:,:,1))
 title='===  Kinetic energy contribution (LIBINT) Y ==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,hamiltonian_kinetic_grad(:,:,2))
 title='===  Kinetic energy contribution (LIBINT) Z ==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,hamiltonian_kinetic_grad(:,:,3))

 call stop_clock(timing_hamiltonian_kin)

end subroutine setup_kinetic_grad


!=========================================================================
subroutine setup_nucleus(print_matrix_,basis,hamiltonian_nucleus)
 use m_atoms
 implicit none
 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: hamiltonian_nucleus(basis%nbf,basis%nbf)
!=====
 integer              :: ishell,jshell
 integer              :: ibf1,ibf2,jbf1,jbf2
 integer              :: natom_local
 integer              :: ibf_cart,jbf_cart,ij
 integer              :: i_cart,j_cart
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 integer              :: iatom
 character(len=100)   :: title
 real(dp),allocatable :: matrix(:,:)
 real(dp)             :: nucleus

 real(C_DOUBLE),allocatable        :: array_cart(:)
 real(C_DOUBLE),allocatable        :: array_cart_C(:)
 integer(C_INT)                    :: amA,contrdepthA
 real(C_DOUBLE)                    :: A(3)
 real(C_DOUBLE),allocatable        :: alphaA(:)
 real(C_DOUBLE),allocatable        :: cA(:)
 integer(C_INT)                    :: amB,contrdepthB
 real(C_DOUBLE)                    :: B(3)
 real(C_DOUBLE),allocatable        :: alphaB(:)
 real(C_DOUBLE),allocatable        :: cB(:)
 real(C_DOUBLE)                    :: C(3)
!=====

 call start_clock(timing_hamiltonian_nuc)
 write(stdout,'(/,a)') ' Setup nucleus-electron part of the Hamiltonian (LIBINT)'
 if( nproc_world > 1 ) then
   natom_local=0
   do iatom=1,natom
     if( rank_world /= MODULO(iatom-1,nproc_world) ) cycle
     natom_local = natom_local + 1
   enddo
   write(stdout,'(a)')         '   Parallelizing over atoms'
   write(stdout,'(a,i5,a,i5)') '   this proc treats ',natom_local,' over ',natom
 endif


 do jshell=1,basis%nshell
   lj      = basis%shell(jshell)%am
   nj_cart = number_basis_function_am('CART',lj)
   nj      = number_basis_function_am(basis%gaussian_type,lj)
   jbf1    = basis%shell(jshell)%istart
   jbf2    = basis%shell(jshell)%iend

   call set_libint_shell(basis%shell(jshell),amB,contrdepthB,B,alphaB,cB)

   do ishell=jshell,basis%nshell
     li      = basis%shell(ishell)%am
     ni_cart = number_basis_function_am('CART',li)
     ni      = number_basis_function_am(basis%gaussian_type,li)
     ibf1    = basis%shell(ishell)%istart
     ibf2    = basis%shell(ishell)%iend

     call set_libint_shell(basis%shell(ishell),amA,contrdepthA,A,alphaA,cA)


     allocate(array_cart(ni_cart*nj_cart))
     allocate(array_cart_C(ni_cart*nj_cart))
     array_cart(:) = 0.0_dp

     do iatom=1,natom
       if( rank_world /= MODULO(iatom-1,nproc_world) ) cycle

       C(:) = x(:,iatom)
#ifdef HAVE_LIBINT_ONEBODY
       call libint_elecpot(amA,contrdepthA,A,alphaA,cA, &
                           amB,contrdepthB,B,alphaB,cB, &
                           C,array_cart_C)
       array_cart(:) = array_cart(:) - zatom(iatom) * array_cart_C(:) 
#else
       ij = 0
       do i_cart=1,ni_cart
         do j_cart=1,nj_cart
           ij = ij + 1
           ibf_cart = basis%shell(ishell)%istart_cart + i_cart - 1
           jbf_cart = basis%shell(jshell)%istart_cart + j_cart - 1
           call nucleus_basis_function(basis%bfc(ibf_cart),basis%bfc(jbf_cart),zatom(iatom),x(:,iatom),nucleus)
           array_cart(ij) = array_cart(ij) + nucleus
         enddo
       enddo
#endif

     enddo
     deallocate(alphaA,cA)

#ifdef HAVE_LIBINT_ONEBODY
     call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart,matrix)
#else
     call transform_molgw_to_molgw(basis%gaussian_type,li,lj,array_cart,matrix)
#endif

     hamiltonian_nucleus(ibf1:ibf2,jbf1:jbf2) = matrix(:,:)
     hamiltonian_nucleus(jbf1:jbf2,ibf1:ibf2) = TRANSPOSE(matrix(:,:))


     deallocate(array_cart,array_cart_C,matrix)

   enddo
   deallocate(alphaB,cB)
 enddo

 !
 ! Reduce operation
 call xsum_world(hamiltonian_nucleus)

 title='===  Nucleus potential contribution (LIBINT) ==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,hamiltonian_nucleus)

 call stop_clock(timing_hamiltonian_nuc)

end subroutine setup_nucleus


!=========================================================================
subroutine setup_nucleus_grad(print_matrix_,basis,hamiltonian_nucleus_grad)
 use m_atoms
 implicit none
 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: hamiltonian_nucleus_grad(basis%nbf,basis%nbf,natom+1,3)
!=====
 integer              :: ishell,jshell
 integer              :: ibf1,ibf2,jbf1,jbf2
 integer              :: natom_local
 integer              :: i_cart,j_cart
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 integer              :: iatom
 character(len=100)   :: title
 real(dp)             :: vnucleus_ij
 real(dp),allocatable :: matrixA(:,:)
 real(dp),allocatable :: matrixB(:,:)

 real(C_DOUBLE),allocatable        :: array_cart_gradAx(:)
 real(C_DOUBLE),allocatable        :: array_cart_gradAy(:)
 real(C_DOUBLE),allocatable        :: array_cart_gradAz(:)
 real(C_DOUBLE),allocatable        :: array_cart_gradBx(:)
 real(C_DOUBLE),allocatable        :: array_cart_gradBy(:)
 real(C_DOUBLE),allocatable        :: array_cart_gradBz(:)
 integer(C_INT)                    :: amA,contrdepthA
 real(C_DOUBLE)                    :: A(3)
 real(C_DOUBLE),allocatable        :: alphaA(:)
 real(C_DOUBLE),allocatable        :: cA(:)
 integer(C_INT)                    :: amB,contrdepthB
 real(C_DOUBLE)                    :: B(3)
 real(C_DOUBLE),allocatable        :: alphaB(:)
 real(C_DOUBLE),allocatable        :: cB(:)
 real(C_DOUBLE)                    :: C(3)
!=====

 call start_clock(timing_hamiltonian_nuc)
 write(stdout,'(/,a)') ' Setup nucleus-electron part of the Hamiltonian gradient (LIBINT)'
 if( nproc_world > 1 ) then
   natom_local=0
   do iatom=1,natom
     if( rank_world /= MODULO(iatom-1,nproc_world) ) cycle
     natom_local = natom_local + 1
   enddo
   write(stdout,'(a)')         '   Parallelizing over atoms'
   write(stdout,'(a,i5,a,i5)') '   this proc treats ',natom_local,' over ',natom
 endif

 hamiltonian_nucleus_grad(:,:,:,:) = 0.0_dp

 do jshell=1,basis%nshell
   lj      = basis%shell(jshell)%am
   nj_cart = number_basis_function_am('CART',lj)
   nj      = number_basis_function_am(basis%gaussian_type,lj)
   jbf1    = basis%shell(jshell)%istart
   jbf2    = basis%shell(jshell)%iend

   call set_libint_shell(basis%shell(jshell),amB,contrdepthB,B,alphaB,cB)

   do ishell=1,basis%nshell
     li      = basis%shell(ishell)%am
     ni_cart = number_basis_function_am('CART',li)
     ni      = number_basis_function_am(basis%gaussian_type,li)
     ibf1    = basis%shell(ishell)%istart
     ibf2    = basis%shell(ishell)%iend

     call set_libint_shell(basis%shell(ishell),amA,contrdepthA,A,alphaA,cA)

     allocate(array_cart_gradAx(ni_cart*nj_cart))
     allocate(array_cart_gradAy(ni_cart*nj_cart))
     allocate(array_cart_gradAz(ni_cart*nj_cart))
     allocate(array_cart_gradBx(ni_cart*nj_cart))
     allocate(array_cart_gradBy(ni_cart*nj_cart))
     allocate(array_cart_gradBz(ni_cart*nj_cart))


     do iatom=1,natom
       if( rank_world /= MODULO(iatom-1,nproc_world) ) cycle


       C(:) = x(:,iatom)

#ifdef HAVE_LIBINT_ONEBODY
       call libint_elecpot_grad(amA,contrdepthA,A,alphaA,cA, &
                                amB,contrdepthB,B,alphaB,cB, &
                                C,                           &
                                array_cart_gradAx,array_cart_gradAy,array_cart_gradAz, &
                                array_cart_gradBx,array_cart_gradBy,array_cart_gradBz)
#else
       call die('nuclear potential gradient not implemented without LIBINT one-body terms')
#endif
       array_cart_gradAx(:) = array_cart_gradAx(:) * (-zatom(iatom))
       array_cart_gradAy(:) = array_cart_gradAy(:) * (-zatom(iatom))
       array_cart_gradAz(:) = array_cart_gradAz(:) * (-zatom(iatom))
       array_cart_gradBx(:) = array_cart_gradBx(:) * (-zatom(iatom))
       array_cart_gradBy(:) = array_cart_gradBy(:) * (-zatom(iatom))
       array_cart_gradBz(:) = array_cart_gradBz(:) * (-zatom(iatom))

       ! X
       call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradAx,matrixA)
       call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradBx,matrixB)
       hamiltonian_nucleus_grad(ibf1:ibf2,jbf1:jbf2,iatom  ,1) = -matrixA(:,:) - matrixB(:,:)
       hamiltonian_nucleus_grad(ibf1:ibf2,jbf1:jbf2,natom+1,1) = hamiltonian_nucleus_grad(ibf1:ibf2,jbf1:jbf2,natom+1,1) + matrixA(:,:)

       ! Y
       call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradAy,matrixA)
       call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradBy,matrixB)
       hamiltonian_nucleus_grad(ibf1:ibf2,jbf1:jbf2,iatom  ,2) = -matrixA(:,:) - matrixB(:,:)
       hamiltonian_nucleus_grad(ibf1:ibf2,jbf1:jbf2,natom+1,2) = hamiltonian_nucleus_grad(ibf1:ibf2,jbf1:jbf2,natom+1,2) + matrixA(:,:)

       ! Z
       call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradAz,matrixA)
       call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradBz,matrixB)
       hamiltonian_nucleus_grad(ibf1:ibf2,jbf1:jbf2,iatom  ,3) = -matrixA(:,:) - matrixB(:,:)
       hamiltonian_nucleus_grad(ibf1:ibf2,jbf1:jbf2,natom+1,3) = hamiltonian_nucleus_grad(ibf1:ibf2,jbf1:jbf2,natom+1,3) + matrixA(:,:)

     enddo
     deallocate(alphaA,cA)


     deallocate(array_cart_gradAx)
     deallocate(array_cart_gradAy)
     deallocate(array_cart_gradAz)
     deallocate(array_cart_gradBx)
     deallocate(array_cart_gradBy)
     deallocate(array_cart_gradBz)
     deallocate(matrixA)
     deallocate(matrixB)

   enddo
   deallocate(alphaB,cB)

 enddo

 !
 ! Reduce operation
 call xsum_world(hamiltonian_nucleus_grad)

 title='===  Nucleus potential contribution (LIBINT) C1X==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,hamiltonian_nucleus_grad(:,:,1,1))
 title='===  Nucleus potential contribution (LIBINT) C1Y==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,hamiltonian_nucleus_grad(:,:,1,2))
 title='===  Nucleus potential contribution (LIBINT) C1Z==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,hamiltonian_nucleus_grad(:,:,1,3))

 call stop_clock(timing_hamiltonian_nuc)

end subroutine setup_nucleus_grad


!=========================================================================
subroutine calculate_dipole_basis(basis,dipole_basis)
 implicit none
 type(basis_set),intent(in)         :: basis
 real(dp),allocatable,intent(out)   :: dipole_basis(:,:,:)
!=====
 integer              :: gt
 integer              :: ishell,jshell
 integer              :: ibf1,ibf2,jbf1,jbf2,ibf1_cart,jbf1_cart
 integer              :: li,lj,ni_cart,nj_cart,i_cart,j_cart
 integer              :: idir
 real(dp),allocatable :: dipole_cart(:,:,:)
!=====
 integer :: unitfile,var_i

 gt = get_gaussian_type_tag(basis%gaussian_type)

 allocate(dipole_basis(basis%nbf,basis%nbf,3))


 do jshell=1,basis%nshell
   lj        = basis%shell(jshell)%am
   nj_cart   = number_basis_function_am('CART',lj)
   jbf1      = basis%shell(jshell)%istart
   jbf1_cart = basis%shell(jshell)%istart_cart
   jbf2      = basis%shell(jshell)%iend

   do ishell=1,basis%nshell
     li        = basis%shell(ishell)%am
     ni_cart   = number_basis_function_am('CART',li)
     ibf1      = basis%shell(ishell)%istart
     ibf1_cart = basis%shell(ishell)%istart_cart
     ibf2      = basis%shell(ishell)%iend


     allocate(dipole_cart(ni_cart,nj_cart,3))

     do i_cart=1,ni_cart
       do j_cart=1,nj_cart
         call basis_function_dipole(basis%bfc(ibf1_cart+i_cart-1),basis%bfc(jbf1_cart+j_cart-1),dipole_cart(i_cart,j_cart,:))
       enddo
     enddo

     do idir=1,3
       dipole_basis(ibf1:ibf2,jbf1:jbf2,idir) = MATMUL( TRANSPOSE( cart_to_pure(li,gt)%matrix(:,:) ) , &
             MATMUL(  dipole_cart(:,:,idir) , cart_to_pure(lj,gt)%matrix(:,:) ) )
     enddo

     deallocate(dipole_cart)

   enddo
 enddo


end subroutine calculate_dipole_basis


!=========================================================================
subroutine calculate_quadrupole_basis(basis,quadrupole_basis)
 implicit none
 type(basis_set),intent(in)         :: basis
 real(dp),allocatable,intent(out)   :: quadrupole_basis(:,:,:,:)
!=====
 integer              :: gt
 integer              :: ishell,jshell
 integer              :: ibf1,ibf2,jbf1,jbf2,ibf1_cart,jbf1_cart
 integer              :: li,lj,ni_cart,nj_cart,i_cart,j_cart
 integer              :: idir,jdir
 real(dp),allocatable :: quadrupole_cart(:,:,:,:)
!=====
 integer :: unitfile,var_i

 gt = get_gaussian_type_tag(basis%gaussian_type)

 allocate(quadrupole_basis(basis%nbf,basis%nbf,3,3))


 do jshell=1,basis%nshell
   lj        = basis%shell(jshell)%am
   nj_cart   = number_basis_function_am('CART',lj)
   jbf1      = basis%shell(jshell)%istart
   jbf1_cart = basis%shell(jshell)%istart_cart
   jbf2      = basis%shell(jshell)%iend

   do ishell=1,basis%nshell
     li        = basis%shell(ishell)%am
     ni_cart   = number_basis_function_am('CART',li)
     ibf1      = basis%shell(ishell)%istart
     ibf1_cart = basis%shell(ishell)%istart_cart
     ibf2      = basis%shell(ishell)%iend


     allocate(quadrupole_cart(ni_cart,nj_cart,3,3))

     do i_cart=1,ni_cart
       do j_cart=1,nj_cart
         call basis_function_quadrupole(basis%bfc(ibf1_cart+i_cart-1),basis%bfc(jbf1_cart+j_cart-1),quadrupole_cart(i_cart,j_cart,:,:))
       enddo
     enddo

     do jdir=1,3
       do idir=1,3
         quadrupole_basis(ibf1:ibf2,jbf1:jbf2,idir,jdir) = MATMUL( TRANSPOSE( cart_to_pure(li,gt)%matrix(:,:) ) , &
               MATMUL(  quadrupole_cart(:,:,idir,jdir) , cart_to_pure(lj,gt)%matrix(:,:) ) )
       enddo
     enddo

     deallocate(quadrupole_cart)

   enddo
 enddo


end subroutine calculate_quadrupole_basis


!=========================================================================
subroutine setup_nucleus_ecp(print_matrix_,basis,hamiltonian_nucleus)
 use m_atoms
 use m_dft_grid
 use m_ecp
 implicit none
 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis
 real(dp),intent(inout)     :: hamiltonian_nucleus(basis%nbf,basis%nbf)
!=====
 integer              :: natom_local
 integer              :: ibf,jbf
 integer              :: iatom
 real(dp)             :: vnucleus_ij

 integer              :: iecp
 integer              :: iproj,nproj
 integer              :: mm
 integer              :: igrid
 real(dp)             :: rr(3)
 real(dp)             :: weight
 real(dp)             :: basis_function_r(basis%nbf)
 integer              :: iradial
 integer              :: i1,n1
 real(dp)             :: xtmp,phi,cos_theta
 real(dp)             :: wxa(nradial_ecp),xa(nradial_ecp)
 real(dp)             :: w1(nangular_ecp),x1(nangular_ecp),y1(nangular_ecp),z1(nangular_ecp)
 real(dp),allocatable :: int_fixed_r(:,:)
 real(dp),external    :: real_spherical_harmonics
 integer              :: necp,ie
 character(len=100)   :: title
 logical              :: element_has_ecp
!=====

 ! Check if there are some ECP
 if( nelement_ecp == 0 ) return


 call start_clock(timing_ecp)

 !
 ! Since there will be an allreduce operation in the end, 
 ! anticipate by dividing the input value of Hnucl by the number of procs
 if( nproc_world > 1 ) then
   hamiltonian_nucleus(:,:) = hamiltonian_nucleus(:,:) / nproc_world
 endif

 n1 = nangular_ecp
 select case(nangular_ecp)
 case(6)
   call ld0006(x1,y1,z1,w1,n1)
 case(14)
   call ld0014(x1,y1,z1,w1,n1)
 case(26)
   call ld0026(x1,y1,z1,w1,n1)
 case(38)
   call ld0038(x1,y1,z1,w1,n1)
 case(50)
   call ld0050(x1,y1,z1,w1,n1)
 case(74)
   call ld0074(x1,y1,z1,w1,n1)
 case(86)
   call ld0086(x1,y1,z1,w1,n1)
 case(110)
   call ld0110(x1,y1,z1,w1,n1)
 case(146)
   call ld0146(x1,y1,z1,w1,n1)
 case(170)
   call ld0170(x1,y1,z1,w1,n1)
 case(230)
   call ld0230(x1,y1,z1,w1,n1)
 case(302)
   call ld0302(x1,y1,z1,w1,n1)
 case(434)
   call ld0434(x1,y1,z1,w1,n1)
 case default
   write(stdout,*) 'grid points: ',nangular_ecp
   call die('setup_nucleus_ecp: Lebedev grid is not available')
 end select


 do iradial=1,nradial_ecp
   xtmp = ( iradial - 0.5_dp ) / REAL(nradial_ecp,dp)
   xa(iradial)   = -5.0_dp * log( 1.0_dp - xtmp**3)
   wxa(iradial)  = 3.0_dp * 5.0_dp * xtmp**2 / ( 1.0_dp - xtmp**3 ) / REAL(nradial_ecp,dp)
 enddo


 do iatom=1,natom
   element_has_ecp = .FALSE.
   do ie=1,nelement_ecp
     if( element_ecp(ie) == basis_element(iatom) ) then
       element_has_ecp = .TRUE.
       exit
     endif
   enddo 

   if( .NOT. element_has_ecp ) cycle

   necp = ecp(ie)%necp
     

   nproj = 0
   do iecp=1,necp
     nproj = nproj + number_basis_function_am('PURE',ecp(ie)%lk(iecp))
   enddo
  
  
   allocate(int_fixed_r(basis%nbf,nproj))
   do iradial=1,nradial_ecp
     if( MODULO(iradial-1,nproc_world) /= rank_world ) cycle

     int_fixed_r(:,:) = 0.0_dp
     do i1=1,nangular_ecp
       rr(1) = xa(iradial) * x1(i1) + x(1,iatom)
       rr(2) = xa(iradial) * y1(i1) + x(2,iatom)
       rr(3) = xa(iradial) * z1(i1) + x(3,iatom)
       call calculate_basis_functions_r(basis,rr,basis_function_r)
  
       cos_theta = z1(i1)
       phi       = ATAN2(y1(i1),x1(i1))
  
       iproj = 0
       do iecp=1,necp
         do mm=-ecp(ie)%lk(iecp),ecp(ie)%lk(iecp)
           iproj = iproj + 1
           int_fixed_r(:,iproj) = int_fixed_r(:,iproj) + basis_function_r(:) &
                                     * real_spherical_harmonics(ecp(ie)%lk(iecp),mm,cos_theta,phi) &
                                        * w1(i1) * 4.0_dp * pi  
         enddo
       enddo
     enddo
  
  
     iproj = 0
     do iecp=1,necp
       do mm=-ecp(ie)%lk(iecp),ecp(ie)%lk(iecp)
         iproj = iproj + 1
         do jbf=1,basis%nbf
           do ibf=1,basis%nbf
             hamiltonian_nucleus(ibf,jbf) = hamiltonian_nucleus(ibf,jbf)  &
                 + int_fixed_r(ibf,iproj) * int_fixed_r(jbf,iproj) * wxa(iradial) * xa(iradial)**2  &
                    * ecp(ie)%dk(iecp) * EXP( -ecp(ie)%zetak(iecp) * xa(iradial)**2 ) * xa(iradial)**(ecp(ie)%nk(iecp)-2)
           enddo
         enddo
       enddo
     enddo
  
   enddo
  
   deallocate(int_fixed_r)

 enddo 

 call xsum_world(hamiltonian_nucleus)

 title='=== ECP Nucleus potential contribution ==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,hamiltonian_nucleus)

 call stop_clock(timing_ecp)

end subroutine setup_nucleus_ecp


end module m_hamiltonian_onebody
!=========================================================================
