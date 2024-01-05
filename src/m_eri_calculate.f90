!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the method to calculate the 2-, 3-, and 4-center Coulomb integrals
!
!=========================================================================
module m_eri_calculate
 use m_definitions
 use m_mpi
 use m_scalapack
 use m_memory
 use m_warning
 use m_basis_set
 use m_timing
 use m_cart_to_pure
 use m_inputparam,only: scalapack_block_min
 use m_eri
 use m_libint_tools


 integer,protected :: nauxil_2center     ! size of the 2-center matrix
 integer,protected :: nauxil_2center_lr  ! size of the 2-center LR matrix

 real(dp),private,allocatable :: eri_2center(:,:)
 real(dp),private,allocatable :: eri_2center_lr(:,:)
 integer,private              :: desc_2center(NDEL)


contains


!=========================================================================
subroutine calculate_eri(print_eri_,basis,rcut)
 implicit none
 logical,intent(in)           :: print_eri_
 type(basis_set),intent(in)   :: basis
 real(dp),intent(in)          :: rcut
!=====

 call start_clock(timing_eri_4center)

 ! ymbyun 2018/05/03
 ! From now on, the array index for 4-center Coulomb integrals can be greater than 2^32.
#ifdef ENABLE_YMBYUN
 write(stdout,'(/,a,i16)') ' Number of integrals to be stored: ',nsize
#else
 write(stdout,'(/,a,i12)') ' Number of integrals to be stored: ',nsize
#endif

 ! ymbyun 2018/05/30
 ! This array is initialized (i.e. first touched) in calculate_eri_4center().
 if( rcut < 1.0e-12_dp ) then
! ymbyun 2021/12/06
! Currently, NVIDIA (PGI) compilers 21.9 do not understand clean_allocate(), leading to FATAL ERROR.
#ifdef _OPENACC
   allocate(eri_4center(nsize))
#else
   call clean_allocate('4-center integrals',eri_4center,nsize)
#endif
!   eri_4center(:) = 0.0_dp
 else
   call clean_allocate('4-center LR integrals',eri_4center_lr,nsize)
!   eri_4center_lr(:) = 0.0_dp
 endif

 if( .NOT. read_eri(rcut) ) then
   call calculate_eri_4center(basis,rcut)
   if( print_eri_ ) then
     call dump_out_eri(rcut)
   endif
 endif

 call stop_clock(timing_eri_4center)

end subroutine calculate_eri


!=========================================================================
subroutine calculate_eri_4center(basis,rcut)
 implicit none
 type(basis_set),intent(in)   :: basis
 real(dp),intent(in)          :: rcut
!=====
 logical                      :: is_longrange
 integer                      :: ishell,jshell,kshell,lshell
 integer                      :: ijshellpair,klshellpair
 integer                      :: n1c,n2c,n3c,n4c
 integer                      :: ni,nj,nk,nl
 integer                      :: ami,amj,amk,aml
 integer                      :: ibf,jbf,kbf,lbf
 real(dp),allocatable         :: integrals(:,:,:,:)
!=====
! variables used to call C
 real(C_DOUBLE)               :: rcut_libint
 integer(C_INT)               :: ng1,ng2,ng3,ng4
 integer(C_INT)               :: am1,am2,am3,am4
 real(C_DOUBLE)               :: x01(3),x02(3),x03(3),x04(3)
 real(C_DOUBLE),allocatable   :: coeff1(:),coeff2(:),coeff3(:),coeff4(:)
 real(C_DOUBLE),allocatable   :: alpha1(:),alpha2(:),alpha3(:),alpha4(:)
 real(C_DOUBLE),allocatable   :: int_shell(:)
!=====

 is_longrange = (rcut > 1.0e-12_dp)
 rcut_libint = rcut
 if( .NOT. is_longrange ) then 
   write(stdout,'(/,a)') ' Calculate and store the 4-center Coulomb integrals'
 else
   write(stdout,'(/,a)') ' Calculate and store the 4-center Long-Range-only Coulomb integrals'
 endif


! ymbyun 2018/05/30
! First touch to reduce NUMA effects using memory affinity
#ifdef ENABLE_OPENMP_AFFINITY
!$OMP PARALLEL
!$OMP DO PRIVATE(klshellpair,ijshellpair, kshell,lshell,ishell,jshell, amk,aml,ami,amj, ni,nj,nk,nl, lbf,kbf,jbf,ibf)
 do klshellpair=1,nshellpair
   kshell = index_shellpair(1,klshellpair)
   lshell = index_shellpair(2,klshellpair)

   amk = basis%shell(kshell)%am
   aml = basis%shell(lshell)%am

   do ijshellpair=1,nshellpair
     ishell = index_shellpair(1,ijshellpair)
     jshell = index_shellpair(2,ijshellpair)

     ami = basis%shell(ishell)%am
     amj = basis%shell(jshell)%am
     if( amk+aml < ami+amj ) cycle

     ni = number_basis_function_am( basis%gaussian_type , ami )
     nj = number_basis_function_am( basis%gaussian_type , amj )
     nk = number_basis_function_am( basis%gaussian_type , amk )
     nl = number_basis_function_am( basis%gaussian_type , aml )

     if( .NOT. is_longrange ) then
       do lbf=1,nl
         do kbf=1,nk
           do jbf=1,nj
             do ibf=1,ni
               eri_4center( index_eri(basis%shell(ishell)%istart+ibf-1, &
                                      basis%shell(jshell)%istart+jbf-1, &
                                      basis%shell(kshell)%istart+kbf-1, &
                                      basis%shell(lshell)%istart+lbf-1) ) = 0.0_dp
             enddo
           enddo
         enddo
       enddo
     else
       do lbf=1,nl
         do kbf=1,nk
           do jbf=1,nj
             do ibf=1,ni
               eri_4center_lr( index_eri(basis%shell(ishell)%istart+ibf-1, &
                                         basis%shell(jshell)%istart+jbf-1, &
                                         basis%shell(kshell)%istart+kbf-1, &
                                         basis%shell(lshell)%istart+lbf-1) ) = 0.0_dp
             enddo
           enddo
         enddo
       enddo
     endif

   enddo

 enddo
!$OMP END DO
#else
 if( .NOT. is_longrange ) then
   eri_4center(:) = 0.0_dp
 else
   eri_4center_lr(:) = 0.0_dp
 endif
!$OMP PARALLEL
#endif

! ymbyun 2018/05/21
! NOTE: Worker threads use a very large private variable (i.e. integrals), so OMP_STACKSIZE should be set appropriately at run time.
!$OMP DO PRIVATE(ishell,jshell,kshell,lshell, ijshellpair,klshellpair, n1c,n2c,n3c,n4c, ni,nj,nk,nl, ami,amj,amk,aml, ibf,jbf,kbf,lbf, integrals, ng1,ng2,ng3,ng4, am1,am2,am3,am4, x01,x02,x03,x04, coeff1,coeff2,coeff3,coeff4, alpha1,alpha2,alpha3,alpha4, int_shell)
 do klshellpair=1,nshellpair
   kshell = index_shellpair(1,klshellpair)
   lshell = index_shellpair(2,klshellpair)

   !
   ! The angular momenta are already ordered so that libint is pleased
   ! 1) amk+aml >= ami+amj
   ! 2) amk>=aml
   ! 3) ami>=amj
   amk = basis%shell(kshell)%am
   aml = basis%shell(lshell)%am

   do ijshellpair=1,nshellpair
     ishell = index_shellpair(1,ijshellpair)
     jshell = index_shellpair(2,ijshellpair)

     ami = basis%shell(ishell)%am
     amj = basis%shell(jshell)%am
     if( amk+aml < ami+amj ) cycle

     ni = number_basis_function_am( basis%gaussian_type , ami )
     nj = number_basis_function_am( basis%gaussian_type , amj )
     nk = number_basis_function_am( basis%gaussian_type , amk )
     nl = number_basis_function_am( basis%gaussian_type , aml )

     am1 = basis%shell(ishell)%am
     am2 = basis%shell(jshell)%am
     am3 = basis%shell(kshell)%am
     am4 = basis%shell(lshell)%am
     n1c = number_basis_function_am( 'CART' , ami )
     n2c = number_basis_function_am( 'CART' , amj )
     n3c = number_basis_function_am( 'CART' , amk )
     n4c = number_basis_function_am( 'CART' , aml )
     ng1 = basis%shell(ishell)%ng
     ng2 = basis%shell(jshell)%ng
     ng3 = basis%shell(kshell)%ng
     ng4 = basis%shell(lshell)%ng
     allocate(alpha1(ng1),alpha2(ng2),alpha3(ng3),alpha4(ng4))
     alpha1(:) = basis%shell(ishell)%alpha(:) 
     alpha2(:) = basis%shell(jshell)%alpha(:)
     alpha3(:) = basis%shell(kshell)%alpha(:)
     alpha4(:) = basis%shell(lshell)%alpha(:)
     x01(:) = basis%shell(ishell)%x0(:)
     x02(:) = basis%shell(jshell)%x0(:)
     x03(:) = basis%shell(kshell)%x0(:)
     x04(:) = basis%shell(lshell)%x0(:)
     allocate(coeff1(basis%shell(ishell)%ng))
     allocate(coeff2(basis%shell(jshell)%ng))
     allocate(coeff3(basis%shell(kshell)%ng))
     allocate(coeff4(basis%shell(lshell)%ng))
     coeff1(:)=basis%shell(ishell)%coeff(:)
     coeff2(:)=basis%shell(jshell)%coeff(:)
     coeff3(:)=basis%shell(kshell)%coeff(:)
     coeff4(:)=basis%shell(lshell)%coeff(:)

     allocate(int_shell(n1c*n2c*n3c*n4c))

     call libint_4center(am1,ng1,x01,alpha1,coeff1, &
                         am2,ng2,x02,alpha2,coeff2, &
                         am3,ng3,x03,alpha3,coeff3, &
                         am4,ng4,x04,alpha4,coeff4, &
                         rcut_libint,int_shell)

     call transform_libint_to_molgw(basis%gaussian_type,ami,amj,amk,aml,int_shell,integrals)

     if( .NOT. is_longrange ) then
       do lbf=1,nl
         do kbf=1,nk
           do jbf=1,nj
             do ibf=1,ni
               eri_4center( index_eri(basis%shell(ishell)%istart+ibf-1, &
                                      basis%shell(jshell)%istart+jbf-1, &
                                      basis%shell(kshell)%istart+kbf-1, &
                                      basis%shell(lshell)%istart+lbf-1) ) = integrals(ibf,jbf,kbf,lbf)
             enddo
           enddo
         enddo
       enddo
     else
       do lbf=1,nl
         do kbf=1,nk
           do jbf=1,nj
             do ibf=1,ni
               eri_4center_lr( index_eri(basis%shell(ishell)%istart+ibf-1, &
                                         basis%shell(jshell)%istart+jbf-1, &
                                         basis%shell(kshell)%istart+kbf-1, &
                                         basis%shell(lshell)%istart+lbf-1) ) = integrals(ibf,jbf,kbf,lbf)
             enddo
           enddo
         enddo
       enddo
     endif

     deallocate(integrals)
     deallocate(int_shell)
     deallocate(alpha1,alpha2,alpha3,alpha4)
     deallocate(coeff1)
     deallocate(coeff2)
     deallocate(coeff3)
     deallocate(coeff4)

   enddo
 enddo
!$OMP END DO
!$OMP END PARALLEL

 write(stdout,'(a,/)') ' All ERI have been calculated'

! ymbyun 2021/12/06
!$acc update device(eri_4center)

end subroutine calculate_eri_4center


!=========================================================================
subroutine calculate_eri_4center_shell(basis,rcut,ijshellpair,klshellpair,&
                                       shellABCD)
 implicit none
 type(basis_set),intent(in)       :: basis
 real(dp),intent(in)              :: rcut
 integer,intent(in)               :: ijshellpair,klshellpair
 real(dp),allocatable,intent(out) :: shellABCD(:,:,:,:)
!=====
 logical                      :: is_longrange
 integer                      :: ishell,jshell,kshell,lshell
 integer                      :: n1c,n2c,n3c,n4c
 integer                      :: ni,nj,nk,nl
 integer                      :: ami,amj,amk,aml
 integer                      :: ibf,jbf,kbf,lbf
!=====
! variables used to call C
 real(C_DOUBLE)               :: rcut_libint
 integer(C_INT)               :: ng1,ng2,ng3,ng4
 integer(C_INT)               :: am1,am2,am3,am4
 real(C_DOUBLE)               :: x01(3),x02(3),x03(3),x04(3)
 real(C_DOUBLE),allocatable   :: coeff1(:),coeff2(:),coeff3(:),coeff4(:)
 real(C_DOUBLE),allocatable   :: alpha1(:),alpha2(:),alpha3(:),alpha4(:)
 real(C_DOUBLE),allocatable   :: shell_libint(:)
!=====

 is_longrange = (rcut > 1.0e-12_dp)
 rcut_libint = rcut



 kshell = index_shellpair(1,klshellpair)
 lshell = index_shellpair(2,klshellpair)

 !
 ! The angular momenta are already ordered so that libint is pleased
 ! 1) amk+aml >= ami+amj
 ! 2) amk>=aml
 ! 3) ami>=amj
 amk = basis%shell(kshell)%am
 aml = basis%shell(lshell)%am


 ishell = index_shellpair(1,ijshellpair)
 jshell = index_shellpair(2,ijshellpair)

 ami = basis%shell(ishell)%am
 amj = basis%shell(jshell)%am
 if( amk+aml < ami+amj ) call die('calculate_4center_shell: wrong ordering')

 ni = number_basis_function_am( basis%gaussian_type , ami )
 nj = number_basis_function_am( basis%gaussian_type , amj )
 nk = number_basis_function_am( basis%gaussian_type , amk )
 nl = number_basis_function_am( basis%gaussian_type , aml )


 am1 = basis%shell(ishell)%am
 am2 = basis%shell(jshell)%am
 am3 = basis%shell(kshell)%am
 am4 = basis%shell(lshell)%am
 n1c = number_basis_function_am( 'CART' , ami )
 n2c = number_basis_function_am( 'CART' , amj )
 n3c = number_basis_function_am( 'CART' , amk )
 n4c = number_basis_function_am( 'CART' , aml )
 ng1 = basis%shell(ishell)%ng
 ng2 = basis%shell(jshell)%ng
 ng3 = basis%shell(kshell)%ng
 ng4 = basis%shell(lshell)%ng
 allocate(alpha1(ng1),alpha2(ng2),alpha3(ng3),alpha4(ng4))
 alpha1(:) = basis%shell(ishell)%alpha(:) 
 alpha2(:) = basis%shell(jshell)%alpha(:)
 alpha3(:) = basis%shell(kshell)%alpha(:)
 alpha4(:) = basis%shell(lshell)%alpha(:)
 x01(:) = basis%shell(ishell)%x0(:)
 x02(:) = basis%shell(jshell)%x0(:)
 x03(:) = basis%shell(kshell)%x0(:)
 x04(:) = basis%shell(lshell)%x0(:)
 allocate(coeff1(basis%shell(ishell)%ng))
 allocate(coeff2(basis%shell(jshell)%ng))
 allocate(coeff3(basis%shell(kshell)%ng))
 allocate(coeff4(basis%shell(lshell)%ng))
 coeff1(:)=basis%shell(ishell)%coeff(:)
 coeff2(:)=basis%shell(jshell)%coeff(:)
 coeff3(:)=basis%shell(kshell)%coeff(:)
 coeff4(:)=basis%shell(lshell)%coeff(:)

 allocate(shell_libint(n1c*n2c*n3c*n4c))

 call libint_4center(am1,ng1,x01,alpha1,coeff1, &
                     am2,ng2,x02,alpha2,coeff2, &
                     am3,ng3,x03,alpha3,coeff3, &
                     am4,ng4,x04,alpha4,coeff4, &
                     rcut_libint, &
                     shell_libint)

 call transform_libint_to_molgw(basis%gaussian_type,ami,amj,amk,aml,shell_libint,shellABCD)


 deallocate(shell_libint)
 deallocate(alpha1,alpha2,alpha3,alpha4)
 deallocate(coeff1)
 deallocate(coeff2)
 deallocate(coeff3)
 deallocate(coeff4)



end subroutine calculate_eri_4center_shell


!=========================================================================
subroutine calculate_eri_4center_shell_grad(basis,rcut,ijshellpair,klshellpair,&
                                            shell_gradA,shell_gradB,shell_gradC,shell_gradD)
 implicit none
 type(basis_set),intent(in)   :: basis
 real(dp),intent(in)          :: rcut
 integer,intent(in)           :: ijshellpair,klshellpair
 real(dp),intent(out)         :: shell_gradA(:,:,:,:,:)
 real(dp),intent(out)         :: shell_gradB(:,:,:,:,:)
 real(dp),intent(out)         :: shell_gradC(:,:,:,:,:)
 real(dp),intent(out)         :: shell_gradD(:,:,:,:,:)
!=====
 logical                      :: is_longrange
 integer                      :: ishell,jshell,kshell,lshell
 integer                      :: n1c,n2c,n3c,n4c
 integer                      :: ni,nj,nk,nl
 integer                      :: ami,amj,amk,aml
 integer                      :: ibf,jbf,kbf,lbf
 real(dp),allocatable         :: grad_tmp(:,:,:,:)
!=====
! variables used to call C
 real(C_DOUBLE)               :: rcut_libint
 integer(C_INT)               :: ng1,ng2,ng3,ng4
 integer(C_INT)               :: am1,am2,am3,am4
 real(C_DOUBLE)               :: x01(3),x02(3),x03(3),x04(3)
 real(C_DOUBLE),allocatable   :: coeff1(:),coeff2(:),coeff3(:),coeff4(:)
 real(C_DOUBLE),allocatable   :: alpha1(:),alpha2(:),alpha3(:),alpha4(:)
 real(C_DOUBLE),allocatable   :: gradAx(:)
 real(C_DOUBLE),allocatable   :: gradAy(:)
 real(C_DOUBLE),allocatable   :: gradAz(:)
 real(C_DOUBLE),allocatable   :: gradBx(:)
 real(C_DOUBLE),allocatable   :: gradBy(:)
 real(C_DOUBLE),allocatable   :: gradBz(:)
 real(C_DOUBLE),allocatable   :: gradCx(:)
 real(C_DOUBLE),allocatable   :: gradCy(:)
 real(C_DOUBLE),allocatable   :: gradCz(:)
 real(C_DOUBLE),allocatable   :: gradDx(:)
 real(C_DOUBLE),allocatable   :: gradDy(:)
 real(C_DOUBLE),allocatable   :: gradDz(:)
!=====

 is_longrange = (rcut > 1.0e-12_dp)
 rcut_libint = rcut



 kshell = index_shellpair(1,klshellpair)
 lshell = index_shellpair(2,klshellpair)

 !
 ! The angular momenta are already ordered so that libint is pleased
 ! 1) amk+aml >= ami+amj
 ! 2) amk>=aml
 ! 3) ami>=amj
 amk = basis%shell(kshell)%am
 aml = basis%shell(lshell)%am


 ishell = index_shellpair(1,ijshellpair)
 jshell = index_shellpair(2,ijshellpair)

 ami = basis%shell(ishell)%am
 amj = basis%shell(jshell)%am
 if( amk+aml < ami+amj ) call die('calculate_4center_shell_grad: wrong ordering')

 ni = number_basis_function_am( basis%gaussian_type , ami )
 nj = number_basis_function_am( basis%gaussian_type , amj )
 nk = number_basis_function_am( basis%gaussian_type , amk )
 nl = number_basis_function_am( basis%gaussian_type , aml )


 am1 = basis%shell(ishell)%am
 am2 = basis%shell(jshell)%am
 am3 = basis%shell(kshell)%am
 am4 = basis%shell(lshell)%am
 n1c = number_basis_function_am( 'CART' , ami )
 n2c = number_basis_function_am( 'CART' , amj )
 n3c = number_basis_function_am( 'CART' , amk )
 n4c = number_basis_function_am( 'CART' , aml )
 ng1 = basis%shell(ishell)%ng
 ng2 = basis%shell(jshell)%ng
 ng3 = basis%shell(kshell)%ng
 ng4 = basis%shell(lshell)%ng
 allocate(alpha1(ng1),alpha2(ng2),alpha3(ng3),alpha4(ng4))
 alpha1(:) = basis%shell(ishell)%alpha(:) 
 alpha2(:) = basis%shell(jshell)%alpha(:)
 alpha3(:) = basis%shell(kshell)%alpha(:)
 alpha4(:) = basis%shell(lshell)%alpha(:)
 x01(:) = basis%shell(ishell)%x0(:)
 x02(:) = basis%shell(jshell)%x0(:)
 x03(:) = basis%shell(kshell)%x0(:)
 x04(:) = basis%shell(lshell)%x0(:)
 allocate(coeff1(basis%shell(ishell)%ng))
 allocate(coeff2(basis%shell(jshell)%ng))
 allocate(coeff3(basis%shell(kshell)%ng))
 allocate(coeff4(basis%shell(lshell)%ng))
 coeff1(:)=basis%shell(ishell)%coeff(:)
 coeff2(:)=basis%shell(jshell)%coeff(:)
 coeff3(:)=basis%shell(kshell)%coeff(:)
 coeff4(:)=basis%shell(lshell)%coeff(:)

 allocate(gradAx(n1c*n2c*n3c*n4c))
 allocate(gradAy(n1c*n2c*n3c*n4c))
 allocate(gradAz(n1c*n2c*n3c*n4c))
 allocate(gradBx(n1c*n2c*n3c*n4c))
 allocate(gradBy(n1c*n2c*n3c*n4c))
 allocate(gradBz(n1c*n2c*n3c*n4c))
 allocate(gradCx(n1c*n2c*n3c*n4c))
 allocate(gradCy(n1c*n2c*n3c*n4c))
 allocate(gradCz(n1c*n2c*n3c*n4c))
 allocate(gradDx(n1c*n2c*n3c*n4c))
 allocate(gradDy(n1c*n2c*n3c*n4c))
 allocate(gradDz(n1c*n2c*n3c*n4c))

#ifdef HAVE_LIBINT_ONEBODY
 call libint_4center_grad(am1,ng1,x01,alpha1,coeff1, &
                          am2,ng2,x02,alpha2,coeff2, &
                          am3,ng3,x03,alpha3,coeff3, &
                          am4,ng4,x04,alpha4,coeff4, &
                          rcut_libint, &
                          gradAx,gradAy,gradAz, &
                          gradBx,gradBy,gradBz, &
                          gradCx,gradCy,gradCz, &
                          gradDx,gradDy,gradDz)
#endif

 call transform_libint_to_molgw(basis%gaussian_type,ami,amj,amk,aml,gradAx,grad_tmp)
 shell_gradA(:,:,:,:,1) = grad_tmp(:,:,:,:) 
 call transform_libint_to_molgw(basis%gaussian_type,ami,amj,amk,aml,gradAy,grad_tmp)
 shell_gradA(:,:,:,:,2) = grad_tmp(:,:,:,:) 
 call transform_libint_to_molgw(basis%gaussian_type,ami,amj,amk,aml,gradAz,grad_tmp)
 shell_gradA(:,:,:,:,3) = grad_tmp(:,:,:,:) 
 call transform_libint_to_molgw(basis%gaussian_type,ami,amj,amk,aml,gradBx,grad_tmp)
 shell_gradB(:,:,:,:,1) = grad_tmp(:,:,:,:) 
 call transform_libint_to_molgw(basis%gaussian_type,ami,amj,amk,aml,gradBy,grad_tmp)
 shell_gradB(:,:,:,:,2) = grad_tmp(:,:,:,:) 
 call transform_libint_to_molgw(basis%gaussian_type,ami,amj,amk,aml,gradBz,grad_tmp)
 shell_gradB(:,:,:,:,3) = grad_tmp(:,:,:,:) 
 call transform_libint_to_molgw(basis%gaussian_type,ami,amj,amk,aml,gradCx,grad_tmp)
 shell_gradC(:,:,:,:,1) = grad_tmp(:,:,:,:) 
 call transform_libint_to_molgw(basis%gaussian_type,ami,amj,amk,aml,gradCy,grad_tmp)
 shell_gradC(:,:,:,:,2) = grad_tmp(:,:,:,:) 
 call transform_libint_to_molgw(basis%gaussian_type,ami,amj,amk,aml,gradCz,grad_tmp)
 shell_gradC(:,:,:,:,3) = grad_tmp(:,:,:,:) 
 call transform_libint_to_molgw(basis%gaussian_type,ami,amj,amk,aml,gradDx,grad_tmp)
 shell_gradD(:,:,:,:,1) = grad_tmp(:,:,:,:) 
 call transform_libint_to_molgw(basis%gaussian_type,ami,amj,amk,aml,gradDy,grad_tmp)
 shell_gradD(:,:,:,:,2) = grad_tmp(:,:,:,:) 
 call transform_libint_to_molgw(basis%gaussian_type,ami,amj,amk,aml,gradDz,grad_tmp)
 shell_gradD(:,:,:,:,3) = grad_tmp(:,:,:,:) 


 deallocate(grad_tmp)
 deallocate(gradAx)
 deallocate(gradAy)
 deallocate(gradAz)
 deallocate(gradBx)
 deallocate(gradBy)
 deallocate(gradBz)
 deallocate(gradCx)
 deallocate(gradCy)
 deallocate(gradCz)
 deallocate(gradDx)
 deallocate(gradDy)
 deallocate(gradDz)
 deallocate(alpha1,alpha2,alpha3,alpha4)
 deallocate(coeff1)
 deallocate(coeff2)
 deallocate(coeff3)
 deallocate(coeff4)



end subroutine calculate_eri_4center_shell_grad


!=========================================================================
subroutine calculate_eri_2center_scalapack(auxil_basis,rcut)
 implicit none
 type(basis_set),intent(in)   :: auxil_basis
 real(dp),intent(in)          :: rcut
!=====
 logical                      :: is_longrange
 integer                      :: ishell,kshell
 integer                      :: n1c,n3c
 integer                      :: ni,nk
 integer                      :: ami,amk
 integer                      :: ibf,kbf
 integer                      :: agt
 integer                      :: info
 integer                      :: ibf_auxil,jbf_auxil
 integer                      :: nauxil_neglect,nauxil_kept
 real(dp)                     :: eigval(auxil_basis%nbf)
 real(dp),allocatable         :: integrals(:,:)
 real(dp)                     :: symmetrization_factor
 real(dp),allocatable         :: eri_2center_sqrt(:,:)
 real(dp),allocatable         :: eri_2center_tmp(:,:)
 integer                      :: mlocal,nlocal
 integer                      :: iglobal,jglobal,ilocal,jlocal
 integer                      :: kglobal,klocal
 integer                      :: desc2center(NDEL)
 logical                      :: skip_shell
!=====
! variables used to call C
 real(C_DOUBLE)               :: rcut_libint
 integer(C_INT)               :: am1,am3
 integer(C_INT)               :: ng1,ng3
 real(C_DOUBLE),allocatable   :: alpha1(:),alpha3(:)
 real(C_DOUBLE)               :: x01(3),x03(3)
 real(C_DOUBLE),allocatable   :: coeff1(:),coeff3(:)
 real(C_DOUBLE),allocatable   :: int_shell(:)
!=====

 call start_clock(timing_eri_2center)



 is_longrange = (rcut > 1.0e-12_dp)
 rcut_libint = rcut
 agt = get_gaussian_type_tag(auxil_basis%gaussian_type)

 if( .NOT. is_longrange ) then
#ifdef HAVE_SCALAPACK
   write(stdout,'(a,i4,a,i4)') ' 2-center integrals distributed using a SCALAPACK grid (LIBINT): ',nprow_3center,' x ',npcol_3center
#else
   write(stdout,'(a)') ' 2-center integrals (LIBINT)'
#endif
 else
#ifdef HAVE_SCALAPACK
   write(stdout,'(a,i4,a,i4)') ' 2-center LR integrals distributed using a SCALAPACK grid (LIBINT): ',nprow_3center,' x ',npcol_3center
#else
   write(stdout,'(a)') ' 2-center LR integrals (LIBINT)'
#endif
 endif

 !
 ! Enforce that all proc participate to the distribution
 ! Much easier to code then
 if( cntxt_3center <= 0 ) call die('distribution not permitted')



 ! Set mlocal => auxil_basis%nbf
 ! Set nlocal => auxil_basis%nbf
 mlocal = NUMROC(auxil_basis%nbf,block_row,iprow_3center,first_row,nprow_3center)
 nlocal = NUMROC(auxil_basis%nbf,block_col,ipcol_3center,first_col,npcol_3center)
 call DESCINIT(desc2center,auxil_basis%nbf,auxil_basis%nbf,block_row,block_col,first_row,first_col,cntxt_3center,MAX(1,mlocal),info)

 call clean_allocate('tmp 2-center integrals',eri_2center_tmp,mlocal,nlocal)


 ! Initialization need since we are going to symmetrize the matrix then
 eri_2center_tmp(:,:) = 0.0_dp


 do kshell=1,auxil_basis%nshell
   amk = auxil_basis%shell(kshell)%am
   nk  = number_basis_function_am( auxil_basis%gaussian_type , amk )

   ! Check if this shell is actually needed for the local matrix
   skip_shell = .TRUE.
   do kbf=1,nk
     kglobal = auxil_basis%shell(kshell)%istart + kbf - 1
     skip_shell = skip_shell .AND. .NOT. ( ipcol_3center == INDXG2P(kglobal,block_col,0,first_col,npcol_3center) )
   enddo

   if( skip_shell ) cycle

 
   do ishell=1,auxil_basis%nshell
     ami = auxil_basis%shell(ishell)%am
     ni = number_basis_function_am( auxil_basis%gaussian_type , ami )

     !
     ! Order the angular momenta so that libint is pleased
     ! 1) am3 >= am1
     if( amk < ami ) cycle
     if( amk == ami ) then
       symmetrization_factor = 0.5_dp
     else
       symmetrization_factor = 1.0_dp
     endif

     ! Check if this shell is actually needed for the local matrix
     skip_shell = .TRUE.
     do ibf=1,ni
       iglobal = auxil_basis%shell(ishell)%istart + ibf - 1
       skip_shell = skip_shell .AND. .NOT. ( iprow_3center == INDXG2P(iglobal,block_row,0,first_row,nprow_3center) )
     enddo

     if( skip_shell ) cycle


     am1 = auxil_basis%shell(ishell)%am
     am3 = auxil_basis%shell(kshell)%am
     n1c = number_basis_function_am( 'CART' , ami )
     n3c = number_basis_function_am( 'CART' , amk )
     allocate( int_shell( n1c*n3c ) )
     ng1 = auxil_basis%shell(ishell)%ng
     ng3 = auxil_basis%shell(kshell)%ng
     allocate(alpha1(ng1),alpha3(ng3))
     alpha1(:) = auxil_basis%shell(ishell)%alpha(:) 
     alpha3(:) = auxil_basis%shell(kshell)%alpha(:)
     x01(:) = auxil_basis%shell(ishell)%x0(:)
     x03(:) = auxil_basis%shell(kshell)%x0(:)
     allocate(coeff1(auxil_basis%shell(ishell)%ng))
     allocate(coeff3(auxil_basis%shell(kshell)%ng))
     coeff1(:)=auxil_basis%shell(ishell)%coeff(:) * cart_to_pure_norm(0,agt)%matrix(1,1)
     coeff3(:)=auxil_basis%shell(kshell)%coeff(:) * cart_to_pure_norm(0,agt)%matrix(1,1)


     call libint_2center(am1,ng1,x01,alpha1,coeff1, &
                         am3,ng3,x03,alpha3,coeff3, &
                         rcut_libint,int_shell)


     deallocate(alpha1,alpha3)
     deallocate(coeff1,coeff3)

     call transform_libint_to_molgw(auxil_basis%gaussian_type,ami,amk,int_shell,integrals)

     deallocate(int_shell)

     
     do kbf=1,nk
       kglobal = auxil_basis%shell(kshell)%istart + kbf - 1

       if( ipcol_3center == INDXG2P(kglobal,block_col,0,first_col,npcol_3center) ) then
         klocal = INDXG2L(kglobal,block_col,0,first_col,npcol_3center)
       else
         cycle
       endif

       do ibf=1,ni
         iglobal = auxil_basis%shell(ishell)%istart + ibf - 1

         if( iprow_3center == INDXG2P(iglobal,block_row,0,first_row,nprow_3center) ) then
           ilocal = INDXG2L(iglobal,block_row,0,first_row,nprow_3center)
         else
           cycle
         endif


         eri_2center_tmp(ilocal,klocal) = integrals(ibf,kbf) * symmetrization_factor

       enddo
     enddo
 
     deallocate(integrals)

   enddo   ! ishell
 enddo   ! kshell

 !
 ! Symmetrize and then diagonalize the 2-center integral matrix
 !
#ifdef HAVE_SCALAPACK

 call clean_allocate('2-center integrals sqrt',eri_2center_sqrt,mlocal,nlocal)

 ! B = A
 call PDLACPY('A',auxil_basis%nbf,auxil_basis%nbf,eri_2center_tmp,1,1,desc2center,eri_2center_sqrt,1,1,desc2center)
 ! A = A + B**T
 call PDGEADD('T',auxil_basis%nbf,auxil_basis%nbf,1.0d0,eri_2center_sqrt,1,1,desc2center,1.0d0,eri_2center_tmp,1,1,desc2center)
 ! Diagonalize
 call diagonalize_sca(auxil_basis%nbf,desc2center,eri_2center_tmp,eigval,desc2center,eri_2center_sqrt)
 call clean_deallocate('tmp 2-center integrals',eri_2center_tmp)

#else

 ! Symmetrize
 eri_2center_tmp(:,:) = eri_2center_tmp(:,:) + TRANSPOSE( eri_2center_tmp(:,:) )
 ! Diagonalize
 call diagonalize_scalapack(scalapack_block_min,auxil_basis%nbf,eri_2center_tmp,eigval)
 call move_alloc(eri_2center_tmp,eri_2center_sqrt)

#endif


 !
 ! Skip the too small eigenvalues
 nauxil_kept = COUNT( eigval(:) > TOO_LOW_EIGENVAL )
 nauxil_neglect = auxil_basis%nbf - nauxil_kept
 if( .NOT. is_longrange ) then
   nauxil_2center = nauxil_kept
   ! Prepare the distribution of the 3-center integrals
   ! nauxil_3center variable is now set up
   call distribute_auxil_basis(nauxil_2center)
 else
   nauxil_2center_lr = nauxil_kept
   ! Prepare the distribution of the 3-center integrals
   ! nauxil_3center_lr variable is now set up
   call distribute_auxil_basis_lr(nauxil_2center_lr)
 endif

 !
 ! Now resize the 2-center matrix accordingly
 ! Set mlocal => nauxil_3center
 ! Set nlocal => nauxil_kept < auxil_basis%nbf
 mlocal = NUMROC(auxil_basis%nbf,block_row,iprow_3center,first_row,nprow_3center)
 nlocal = NUMROC(nauxil_kept    ,block_col,ipcol_3center,first_col,npcol_3center)
 call DESCINIT(desc_2center,auxil_basis%nbf,nauxil_kept,block_row,block_col,first_row,first_col,cntxt_3center,MAX(1,mlocal),info)


 if( .NOT. is_longrange ) then
   call clean_allocate('Distributed 2-center integrals',eri_2center,mlocal,nlocal)
 else
   call clean_allocate('Distributed LR 2-center integrals',eri_2center_lr,mlocal,nlocal)
 endif

#ifdef HAVE_SCALAPACK
 call clean_allocate('tmp 2-center integrals',eri_2center_tmp,mlocal,nlocal)
 !
 ! Create a rectangular matrix with only 1 / SQRT( eigval) on a diagonal
 eri_2center_tmp(:,:) = 0.0_dp
 do jlocal=1,nlocal
   jglobal = INDXL2G(jlocal,block_col,ipcol_3center,first_col,npcol_3center)
   do ilocal=1,mlocal
     iglobal = INDXL2G(ilocal,block_row,iprow_3center,first_row,nprow_3center)

     if( iglobal == jglobal + nauxil_neglect ) eri_2center_tmp(ilocal,jlocal) = 1.0_dp / SQRT( eigval(jglobal+nauxil_neglect) )

   enddo
 enddo


 if( .NOT. is_longrange ) then
   call PDGEMM('N','N',auxil_basis%nbf,nauxil_2center,auxil_basis%nbf, &
               1.0_dp,eri_2center_sqrt ,1,1,desc2center,  &
                      eri_2center_tmp,1,1,desc_2center,   &
               0.0_dp,eri_2center    ,1,1,desc_2center)
 else
   call PDGEMM('N','N',auxil_basis%nbf,nauxil_2center_lr,auxil_basis%nbf, &
               1.0_dp,eri_2center_sqrt ,1,1,desc2center,  &
                      eri_2center_tmp,1,1,desc_2center,   &
               0.0_dp,eri_2center_lr ,1,1,desc_2center)
 endif

 call clean_deallocate('tmp 2-center integrals',eri_2center_tmp)

#else

 ilocal = 0
 do jlocal=1,auxil_basis%nbf
   if( eigval(jlocal) < TOO_LOW_EIGENVAL ) cycle
   ilocal = ilocal + 1
   eri_2center_sqrt(:,ilocal) = eri_2center_sqrt(:,jlocal) / SQRT( eigval(jlocal) )
 enddo

 if( .NOT. is_longrange ) then
   do ibf_auxil=1,nauxil_3center
     jbf_auxil = ibf_auxil_g(ibf_auxil)
     eri_2center(:,ibf_auxil) = eri_2center_sqrt(:,jbf_auxil)
   enddo
 else
   do ibf_auxil=1,nauxil_3center_lr
     jbf_auxil = ibf_auxil_g_lr(ibf_auxil)
     eri_2center_lr(:,ibf_auxil) = eri_2center_sqrt(:,jbf_auxil)
   enddo
 endif

#endif

 call clean_deallocate('2-center integrals sqrt',eri_2center_sqrt)



 write(stdout,'(/,1x,a)')      'All 2-center integrals have been calculated, diagonalized and stored'
 write(stdout,'(1x,a,i6)')     'Some have been eliminated due to too large overlap ',nauxil_neglect
 write(stdout,'(1x,a,es16.6)') 'because their eigenvalue was lower than:',TOO_LOW_EIGENVAL


 call stop_clock(timing_eri_2center)

end subroutine calculate_eri_2center_scalapack



!=========================================================================
subroutine calculate_eri_3center_scalapack(basis,auxil_basis,rcut)
 implicit none
 type(basis_set),intent(in)   :: basis
 type(basis_set),intent(in)   :: auxil_basis
 real(dp),intent(in)          :: rcut
!=====
 logical                      :: is_longrange
 integer                      :: agt
 integer                      :: ishell,kshell,lshell
 integer                      :: klshellpair
 integer                      :: n1c,n3c,n4c
 integer                      :: ni,nk,nl
 integer                      :: ami,amk,aml
 integer                      :: ibf,kbf,lbf
 integer                      :: ibf_auxil,jbf_auxil
 integer                      :: info
 real(dp),allocatable         :: integrals(:,:,:)
 real(dp),allocatable         :: eri_3center_tmp(:,:)
 integer                      :: klpair_global
 integer                      :: ibf_auxil_local,ibf_auxil_global
 integer                      :: mlocal,nlocal
 integer                      :: iglobal,jglobal,ilocal,jlocal
 integer                      :: kglobal,klocal
 integer                      :: desc3center(NDEL)
 integer                      :: desc3tmp(NDEL)
 integer                      :: desc3final(NDEL)
 integer                      :: nauxil_kept
 logical                      :: skip_shell
!=====
! variables used to call C
 real(C_DOUBLE)               :: rcut_libint
 integer(C_INT)               :: am1,am3,am4
 integer(C_INT)               :: ng1,ng3,ng4
 real(C_DOUBLE),allocatable   :: alpha1(:),alpha3(:),alpha4(:)
 real(C_DOUBLE)               :: x01(3),x03(3),x04(3)
 real(C_DOUBLE),allocatable   :: coeff1(:),coeff3(:),coeff4(:)
 real(C_DOUBLE),allocatable   :: int_shell(:)
!=====

 call start_clock(timing_eri_3center)

 is_longrange = (rcut > 1.0e-12_dp)
 rcut_libint = rcut
 agt = get_gaussian_type_tag(auxil_basis%gaussian_type)

 if( .NOT. is_longrange ) then
   nauxil_kept = nauxil_2center
 else
   nauxil_kept = nauxil_2center_lr
 endif

 if( .NOT. is_longrange ) then
   write(stdout,'(/,a)')    ' Calculate and store all the 3-center Electron Repulsion Integrals (LIBINT 3center)'
 else
   write(stdout,'(/,a)')    ' Calculate and store all the LR 3-center Electron Repulsion Integrals (LIBINT 3center)'
 endif
#ifdef HAVE_SCALAPACK
 write(stdout,'(a,i4,a,i4)') ' 3-center integrals distributed using a SCALAPACK grid: ',nprow_3center,' x ',npcol_3center
#endif

 !
 ! Enforce that all proc participate to the distribution
 ! Much easier to code then
 if( cntxt_3center <= 0 ) call die('distribution not permitted')


 ! Set mlocal => auxil_basis%nbf
 ! Set nlocal => npair
 mlocal = NUMROC(auxil_basis%nbf,block_row,iprow_3center,first_row,nprow_3center)
! ymbyun 2018/05/21
! From now on, npair is an 8-byte integer.
#ifdef ENABLE_YMBYUN
 nlocal = NUMROC(INT(npair)     ,block_col,ipcol_3center,first_col,npcol_3center)
 call DESCINIT(desc3center,auxil_basis%nbf,INT(npair),block_row,block_col,first_row,first_col,cntxt_3center,MAX(1,mlocal),info)
#else
 nlocal = NUMROC(npair          ,block_col,ipcol_3center,first_col,npcol_3center)
 call DESCINIT(desc3center,auxil_basis%nbf,npair,block_row,block_col,first_row,first_col,cntxt_3center,MAX(1,mlocal),info)
#endif

 !  Allocate the 3-center integral array
 !
 ! 3-CENTER INTEGRALS 
 !
 call clean_allocate('TMP 3-center integrals',eri_3center_tmp,mlocal,nlocal)



 do klshellpair=1,nshellpair
   kshell = index_shellpair(1,klshellpair)
   lshell = index_shellpair(2,klshellpair)

   amk = basis%shell(kshell)%am
   aml = basis%shell(lshell)%am
   nk = number_basis_function_am( basis%gaussian_type , amk )
   nl = number_basis_function_am( basis%gaussian_type , aml )


   ! Check if this shell is actually needed for the local matrix
   skip_shell = .TRUE.
   do lbf=1,nl
     do kbf=1,nk
       klpair_global = index_pair(basis%shell(kshell)%istart+kbf-1,basis%shell(lshell)%istart+lbf-1)

       skip_shell = skip_shell .AND. .NOT. ( ipcol_3center == INDXG2P(klpair_global,block_col,0,first_col,npcol_3center) )
     enddo
   enddo

   if( skip_shell ) cycle


   do ishell=1,auxil_basis%nshell

     ami = auxil_basis%shell(ishell)%am
     ni = number_basis_function_am( auxil_basis%gaussian_type , ami )

     ! Check if this shell is actually needed for the local matrix
     skip_shell = .TRUE.
     do ibf=1,ni
       iglobal = auxil_basis%shell(ishell)%istart + ibf - 1
       skip_shell = skip_shell .AND. .NOT. ( iprow_3center == INDXG2P(iglobal,block_row,0,first_row,nprow_3center) )
     enddo

     if( skip_shell ) cycle


     am1 = ami
     am3 = amk
     am4 = aml
     n1c = number_basis_function_am( 'CART' , ami )
     n3c = number_basis_function_am( 'CART' , amk )
     n4c = number_basis_function_am( 'CART' , aml )
     ng1 = auxil_basis%shell(ishell)%ng
     ng3 = basis%shell(kshell)%ng
     ng4 = basis%shell(lshell)%ng
     allocate(alpha1(ng1),alpha3(ng3),alpha4(ng4))
     allocate(coeff1(ng1),coeff3(ng3),coeff4(ng4))
     alpha1(:) = auxil_basis%shell(ishell)%alpha(:) 
     alpha3(:) = basis%shell(kshell)%alpha(:)
     alpha4(:) = basis%shell(lshell)%alpha(:)
     coeff1(:) = auxil_basis%shell(ishell)%coeff(:) * cart_to_pure_norm(0,agt)%matrix(1,1)
     coeff3(:) = basis%shell(kshell)%coeff(:)
     coeff4(:) = basis%shell(lshell)%coeff(:)
     x01(:) = auxil_basis%shell(ishell)%x0(:)
     x03(:) = basis%shell(kshell)%x0(:)
     x04(:) = basis%shell(lshell)%x0(:)


     allocate( int_shell(n1c*n3c*n4c) )

     call libint_3center(am1,ng1,x01,alpha1,coeff1, &
                         am3,ng3,x03,alpha3,coeff3, &
                         am4,ng4,x04,alpha4,coeff4, &
                         rcut_libint,int_shell)

     call transform_libint_to_molgw(auxil_basis%gaussian_type,ami,basis%gaussian_type,amk,aml,int_shell,integrals)

     
     do lbf=1,nl
       do kbf=1,nk
         klpair_global = index_pair(basis%shell(kshell)%istart+kbf-1,basis%shell(lshell)%istart+lbf-1)
         if( ipcol_3center /= INDXG2P(klpair_global,block_col,0,first_col,npcol_3center) ) cycle
         jlocal = INDXG2L(klpair_global,block_col,0,first_col,npcol_3center)

         do ibf=1,ni
           ibf_auxil_global = auxil_basis%shell(ishell)%istart+ibf-1
           if( iprow_3center /= INDXG2P(ibf_auxil_global,block_row,0,first_row,nprow_3center) ) cycle
           ilocal = INDXG2L(ibf_auxil_global,block_row,0,first_row,nprow_3center)

           eri_3center_tmp(ilocal,jlocal) = integrals(ibf,kbf,lbf)

         enddo
       enddo
     enddo


     deallocate(integrals)
     deallocate(int_shell)
     deallocate(alpha1,alpha3,alpha4)
     deallocate(coeff1,coeff3,coeff4)

   enddo ! ishell

 enddo ! klshellpair


 ! Set mlocal => nauxil_kept = nauxil_2center OR nauxil_2center_lr
 ! Set nlocal => npair
 mlocal = NUMROC(nauxil_kept,block_row,iprow_3center,first_row,nprow_3center)

! ymbyun 2018/05/21
! From now on, npair is an 8-byte integer.
#ifdef ENABLE_YMBYUN
 nlocal = NUMROC(INT(npair) ,block_col,ipcol_3center,first_col,npcol_3center)
 call DESCINIT(desc3tmp,nauxil_kept,INT(npair),block_row,block_col,first_row,first_col,cntxt_3center,MAX(1,mlocal),info)
#else
 nlocal = NUMROC(npair      ,block_col,ipcol_3center,first_col,npcol_3center)
 call DESCINIT(desc3tmp,nauxil_kept,npair,block_row,block_col,first_row,first_col,cntxt_3center,MAX(1,mlocal),info)
#endif


#ifdef HAVE_SCALAPACK
 call clean_allocate('3-center integrals SCALAPACK',eri_3center_sca,mlocal,nlocal)

 if( .NOT. is_longrange ) then
   call PDGEMM('T','N',nauxil_kept,npair,auxil_basis%nbf, &
               1.0_dp,eri_2center    ,1,1,desc_2center,  &
                      eri_3center_tmp,1,1,desc3center,   &
               0.0_dp,eri_3center_sca,1,1,desc3tmp)
 else
   call PDGEMM('T','N',nauxil_kept,npair,auxil_basis%nbf, &
               1.0_dp,eri_2center_lr ,1,1,desc_2center,  &
                      eri_3center_tmp,1,1,desc3center,   &
               0.0_dp,eri_3center_sca,1,1,desc3tmp)
 endif

 write(stdout,'(a,i8,a,i4)') ' Final 3-center integrals distributed using a SCALAPACK grid: ',nprow_auxil,' x ',npcol_auxil

 if( cntxt_auxil > 0 ) then
   mlocal = NUMROC(nauxil_kept,MBLOCK_AUXIL,iprow_auxil,first_row,nprow_auxil)
   nlocal = NUMROC(npair      ,NBLOCK_AUXIL,ipcol_auxil,first_col,npcol_auxil)
 else
   mlocal = -1
   nlocal = -1
 endif
 call xmax_ortho(mlocal)
 call xmax_ortho(nlocal)

 call clean_deallocate('TMP 3-center integrals',eri_3center_tmp)

 !
 ! Final distribution on Full-Range or Long-Range arrays
 !
 call DESCINIT(desc3final,nauxil_kept,npair,MBLOCK_AUXIL,NBLOCK_AUXIL,first_row,first_col,cntxt_auxil,MAX(1,mlocal),info)
 if( .NOT. is_longrange ) then
   call clean_deallocate('Distributed 2-center integrals',eri_2center)
   call clean_allocate('3-center integrals',eri_3center,mlocal,nlocal)
   call PDGEMR2D(nauxil_kept,npair,eri_3center_sca,1,1,desc3tmp,eri_3center,1,1,desc3final,cntxt_3center)
   !
   ! Propagate to the ortho MPI direction
   if( cntxt_auxil <= 0 ) then
     eri_3center(:,:) = 0.0_dp
   endif
   call xsum_ortho(eri_3center)

   ! Development version
#ifndef SCASCA
   call clean_deallocate('3-center integrals SCALAPACK',eri_3center_sca)
#endif

 else
   call clean_deallocate('Distributed LR 2-center integrals',eri_2center_lr)
   call clean_allocate('LR 3-center integrals',eri_3center_lr,mlocal,nlocal)
   call PDGEMR2D(nauxil_kept,npair,eri_3center_sca,1,1,desc3tmp,eri_3center_lr,1,1,desc3final,cntxt_3center)
   !
   ! Propagate to the ortho MPI direction
   if( cntxt_auxil <= 0 ) then
     eri_3center_lr(:,:) = 0.0_dp
   endif
   call xsum_ortho(eri_3center_lr)

  ! Development version
#ifndef SCASCA
   call clean_deallocate('3-center integrals SCALAPACK',eri_3center_sca)
#endif
 endif


#else

 if( .NOT. is_longrange ) then

   call clean_allocate('3-center integrals',eri_3center,mlocal,nlocal)

   call DGEMM('T','N',nauxil_kept,npair,auxil_basis%nbf, &
               1.0_dp,eri_2center,auxil_basis%nbf,eri_3center_tmp,auxil_basis%nbf,0.0_dp,eri_3center,nauxil_kept)

   write(stdout,'(a)') ' Final 3-center integrals not distributed! since there is no SCALAPACK'

   call clean_deallocate('Distributed 2-center integrals',eri_2center)

 else

   call clean_allocate('LR 3-center integrals',eri_3center_lr,mlocal,nlocal)

   call DGEMM('T','N',nauxil_kept,npair,auxil_basis%nbf, &
               1.0_dp,eri_2center_lr,auxil_basis%nbf,eri_3center_tmp,auxil_basis%nbf,0.0_dp,eri_3center_lr,nauxil_kept)

   write(stdout,'(a)') ' Final LR 3-center integrals not distributed! since there is no SCALAPACK'

   call clean_deallocate('Distributed 2-center integrals',eri_2center_lr)

 endif
 call clean_deallocate('TMP 3-center integrals',eri_3center_tmp)

 
#endif



 write(stdout,'(/,1x,a)') 'All 3-center integrals have been calculated and stored'


 call stop_clock(timing_eri_3center)


end subroutine calculate_eri_3center_scalapack


!=========================================================================
subroutine calculate_eri_approximate_hartree(basis,mv,nv,x0_rho,ng_rho,coeff_rho,alpha_rho,vhrho)
 implicit none
 type(basis_set),intent(in)   :: basis
 integer,intent(in)           :: mv,nv
 real(dp),intent(in)          :: x0_rho(3)
 integer,intent(in)           :: ng_rho
 real(dp),intent(in)          :: coeff_rho(ng_rho),alpha_rho(ng_rho)
 real(dp),intent(inout)       :: vhrho(mv,nv)
!=====
 integer                      :: kshell,lshell
 integer                      :: klshellpair
 integer                      :: n3c,n4c
 integer                      :: nk,nl
 integer                      :: amk,aml
 integer                      :: kbf,lbf
 real(dp),allocatable         :: integrals(:,:,:)
 integer                      :: ilocal,jlocal,iglobal,jglobal
!=====
! variables used to call C
 integer(C_INT)               :: am1,am3,am4
 integer(C_INT)               :: ng1,ng3,ng4
 real(C_DOUBLE),allocatable   :: alpha1(:),alpha3(:),alpha4(:)
 real(C_DOUBLE)               :: x01(3),x03(3),x04(3)
 real(C_DOUBLE),allocatable   :: coeff1(:),coeff3(:),coeff4(:)
 real(C_DOUBLE),allocatable   :: int_shell(:)
!=====

 if( parallel_ham .AND. parallel_buffer ) then    
   if( mv /= basis%nbf .OR. nv /= basis%nbf ) call die('calculate_eri_approximate_hartree: wrong dimension for the buffer')
 endif

 ! Nullify vhrho just for safety.
 ! I guess this is useless.
 if( .NOT. ( parallel_ham .AND. parallel_buffer )  ) then    
   vhrho(:,:) = 0.0_dp
 endif

 do klshellpair=1,nshellpair
   kshell = index_shellpair(1,klshellpair)
   lshell = index_shellpair(2,klshellpair)

   amk = basis%shell(kshell)%am
   aml = basis%shell(lshell)%am

   nk = number_basis_function_am( basis%gaussian_type , amk )
   nl = number_basis_function_am( basis%gaussian_type , aml )

   am1 = 0
   am3 = amk
   am4 = aml
   n3c = number_basis_function_am( 'CART' , amk )
   n4c = number_basis_function_am( 'CART' , aml )
   ng1 = ng_rho
   ng3 = basis%shell(kshell)%ng
   ng4 = basis%shell(lshell)%ng
   allocate(alpha1(ng1),alpha3(ng3),alpha4(ng4))
   allocate(coeff1(ng1),coeff3(ng3),coeff4(ng4))
   alpha1(:) = alpha_rho(:)
   alpha3(:) = basis%shell(kshell)%alpha(:)
   alpha4(:) = basis%shell(lshell)%alpha(:)
   coeff1(:) = coeff_rho(:) / 2.0_dp**1.25_dp / pi**0.75_dp * alpha_rho(:)**1.5_dp * cart_to_pure_norm(0,PUREG)%matrix(1,1)
   coeff3(:) = basis%shell(kshell)%coeff(:)
   coeff4(:) = basis%shell(lshell)%coeff(:)
   x01(:) = x0_rho(:)
   x03(:) = basis%shell(kshell)%x0(:)
   x04(:) = basis%shell(lshell)%x0(:)

   allocate( int_shell(n3c*n4c) )

   call libint_3center(am1,ng1,x01,alpha1,coeff1, &
                       am3,ng3,x03,alpha3,coeff3, &
                       am4,ng4,x04,alpha4,coeff4, &
                       0.0_C_DOUBLE,int_shell)

   call transform_libint_to_molgw(basis%gaussian_type,0,basis%gaussian_type,amk,aml,int_shell,integrals)

     

   if( parallel_ham .AND. parallel_buffer ) then    
     do lbf=1,nl
       do kbf=1,nk
         iglobal = basis%shell(kshell)%istart+kbf-1
         jglobal = basis%shell(lshell)%istart+lbf-1
         ilocal = iglobal
         jlocal = jglobal
         if( kshell == lshell ) then ! To avoid double-counting   
           vhrho(ilocal,jlocal) = vhrho(ilocal,jlocal) + integrals(1,kbf,lbf)  * 0.5_dp 
         else
           vhrho(ilocal,jlocal) = vhrho(ilocal,jlocal) + integrals(1,kbf,lbf) 
         endif
       enddo
     enddo

   else

     do lbf=1,nl
       do kbf=1,nk
         iglobal = basis%shell(kshell)%istart+kbf-1
         jglobal = basis%shell(lshell)%istart+lbf-1
         ilocal = rowindex_global_to_local('H',iglobal)
         jlocal = colindex_global_to_local('H',jglobal)
         if( ilocal*jlocal /= 0 ) then
           vhrho(ilocal,jlocal) = integrals(1,kbf,lbf)
         endif
         ! And the symmetric too
         iglobal = basis%shell(lshell)%istart+lbf-1
         jglobal = basis%shell(kshell)%istart+kbf-1
         ilocal = rowindex_global_to_local('H',iglobal)
         jlocal = colindex_global_to_local('H',jglobal)
         if( ilocal*jlocal /= 0 ) then
           vhrho(ilocal,jlocal) = integrals(1,kbf,lbf)
         endif
       enddo
     enddo
   endif


   deallocate(integrals)
   deallocate(int_shell)
   deallocate(alpha1,alpha3,alpha4)
   deallocate(coeff1,coeff3,coeff4)

 enddo


end subroutine calculate_eri_approximate_hartree


!=========================================================================
end module m_eri_calculate
