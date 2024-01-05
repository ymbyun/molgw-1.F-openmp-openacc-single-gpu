!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains
! the perturbation theory to 2nd order evaluation of the self-energy
!
!=========================================================================
subroutine pt2_selfenergy(selfenergy_approx,nstate,basis,occupation,energy,c_matrix,se,emp2)
! ymbyun 2021/12/27
#ifdef _OPENACC
#if _OPENACC >= 201711
 use openacc
#endif
#endif
 use m_definitions
 use m_mpi
 use m_mpi_ortho
 use m_warning
 use m_timing
 use m_basis_set
 use m_eri_ao_mo
 use m_inputparam
 use m_selfenergy_tools
 implicit none

 integer,intent(in)         :: selfenergy_approx,nstate
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: occupation(nstate,nspin),energy(nstate,nspin)
 real(dp),intent(in)        :: c_matrix(basis%nbf,nstate,nspin)
 type(selfenergy_grid),intent(inout) :: se
 real(dp),intent(out)       :: emp2
!=====
 integer               :: pstate,qstate
 real(dp),allocatable  :: selfenergy_ring(:,:,:)
 real(dp),allocatable  :: selfenergy_sox(:,:,:)
 integer               :: iomega
 integer               :: istate,jstate,kstate
 integer               :: pqispin,jkspin
 real(dp)              :: fact_occ1,fact_occ2
 real(dp)              :: fi,fj,fk,ei,ej,ek
 real(dp)              :: omega
 real(dp)              :: fact_real,fact_energy
 real(dp)              :: emp2_sox,emp2_ring
 real(dp),allocatable  :: eri_eigenstate_i(:,:,:,:)
 real(dp)              :: coul_iqjk,coul_ijkq,coul_ipkj
! ymbyun 2021/12/27
#ifdef _OPENACC
#if _OPENACC >= 201711
 integer                  :: device_num
 integer(acc_device_kind) :: device_type
 integer(kind=8)          :: gpu_memory   ,gpu_free_memory
 integer                  :: gpu_memory_mb,gpu_free_memory_mb
#endif
#endif
!=====

 call start_clock(timing_mp2_self)

 emp2_ring = 0.0_dp
 emp2_sox  = 0.0_dp


 write(stdout,'(/,a)') ' Perform the second-order self-energy calculation'
 write(stdout,*) 'with the perturbative approach'

 

 if(has_auxil_basis) then
   call calculate_eri_3center_eigen(basis%nbf,nstate,c_matrix,ncore_G+1,nvirtual_G-1,ncore_G+1,nvirtual_G-1)
 else
   allocate(eri_eigenstate_i(nstate,nstate,nstate,nspin))
 endif



 allocate(selfenergy_ring(-se%nomega:se%nomega,nsemin:nsemax,nspin))
 allocate(selfenergy_sox (-se%nomega:se%nomega,nsemin:nsemax,nspin))


 selfenergy_ring(:,:,:) = 0.0_dp
 selfenergy_sox(:,:,:)  = 0.0_dp

! ymbyun 2021/12/07
! OpenACC is added.
! NOTE: DATA can be either here or at the outermost loop.
!       Unlike chi_to_vchiv() in linear_response.f90, the position of DATA has no effect on the performance.

#if !defined (_OPENMP) && defined (_OPENACC)
 !$acc data copyin(occupation,energy,eri_eigenstate_i,se,se%omega) copy(emp2_ring,emp2_sox,selfenergy_ring,selfenergy_sox)
#endif

#ifdef _OPENACC
#if _OPENACC >= 201711
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
#endif
#endif

 do pqispin=1,nspin
   do istate=ncore_G+1,nvirtual_G-1 !LOOP of the first Green's function
     if( MODULO( istate - (ncore_G+1) , nproc_ortho ) /= rank_ortho ) cycle

     if( .NOT. has_auxil_basis ) then
       call calculate_eri_4center_eigen(basis%nbf,nstate,c_matrix,istate,pqispin,eri_eigenstate_i)
     endif

     fi = occupation(istate,pqispin)
     ei = energy(istate,pqispin)

! ymbyun 2019/09/23
! NOTE: OpenMP 2.0 (2000) supports array reductions for Fortran, while OpenMP 4.5 (2015) supports array reductions for C/C++.
! NOTE: If pstate is moved into the kstate loop, performance increases slightly but readability decreases significantly.

! ymbyun 2020/10/19
! NOTE: The following things are changed to make OpenMP work for KISTI Nurion with Intel Xeon Phi CPUs and Intel compilers 18.0.3.
! NOTE: 1. OMP is moved up from the jkspin loop to the pstate loop.
! NOTE: 2. qstate is added to PRIVATE.
! NOTE: 3. Array reduction is removed.
! NOTE: 4. COLLAPSE is removed.
! TODO: Currently, the OpenMP speedup is minimal when only HOMO and LUMO are calculated.
! TODO: A way to add COLLAPSE should be found to utilize as many cores as possbile.

! ymbyun 2021/12/07
! OpenACC is added.
! NOTE: Reduction for emp2_ring and emp2_sox is implicitly generated.

#if defined (_OPENMP) && !defined (_OPENACC)
 !$OMP PARALLEL
 !$OMP DO PRIVATE(pstate,qstate, jkspin, jstate,fj,ej, kstate,fk,ek,fact_occ1,fact_occ2,coul_ipkj,coul_iqjk,coul_ijkq, iomega,omega,fact_real,fact_energy) REDUCTION(+:emp2_ring,emp2_sox)
#endif
#if !defined (_OPENMP) && defined (_OPENACC)
 !$acc parallel loop independent
#endif

!!REDUCTION(+:emp2_ring,emp2_sox,selfenergy_ring,selfenergy_sox) COLLAPSE(2)

     do pstate=nsemin,nsemax ! external loop ( bra )
       qstate=pstate         ! external loop ( ket )
       !$acc loop seq
       do jkspin=1,nspin
         !$acc loop seq
         do jstate=ncore_G+1,nvirtual_G-1  !LOOP of the second Green's function
           fj = occupation(jstate,jkspin)
           ej = energy(jstate,jkspin)
           !$acc loop seq
           do kstate=ncore_G+1,nvirtual_G-1 !LOOP of the third Green's function
             fk = occupation(kstate,jkspin)
             ek = energy(kstate,jkspin)

             fact_occ1 = (spin_fact-fi) *            fj  * (spin_fact-fk) / spin_fact**3
             fact_occ2 =            fi  * (spin_fact-fj) *            fk  / spin_fact**3

             if( fact_occ1 < completely_empty .AND. fact_occ2 < completely_empty ) cycle

! ymbyun 2021/12/07
! This is an ugly way to avoid accessing eri_eigen_ri().
! Due to this, n2_mp2_ri.in in the test suite fails.

#ifdef _OPENACC
             if( has_auxil_basis ) then
!               call die('(ymbyun) Currently, OpenACC cannot be used together with RI and MPI')
               write(stdout,*) '(ymbyun) Currently, OpenACC cannot be used together with RI and MPI'
               stop
             else
               coul_ipkj = eri_eigenstate_i(pstate,kstate,jstate,jkspin)
               coul_iqjk = eri_eigenstate_i(qstate,jstate,kstate,jkspin)
               if( pqispin == jkspin ) then
                 coul_ijkq = eri_eigenstate_i(jstate,kstate,qstate,pqispin) 
               endif
             endif
#else
             if( has_auxil_basis ) then
               coul_ipkj = eri_eigen_ri(istate,pstate,pqispin,kstate,jstate,jkspin)
               coul_iqjk = eri_eigen_ri(istate,qstate,pqispin,jstate,kstate,jkspin)
               if( pqispin == jkspin ) then
                 coul_ijkq = eri_eigen_ri(istate,jstate,pqispin,kstate,qstate,pqispin) 
               endif
             else
               coul_ipkj = eri_eigenstate_i(pstate,kstate,jstate,jkspin)
               coul_iqjk = eri_eigenstate_i(qstate,jstate,kstate,jkspin)
               if( pqispin == jkspin ) then
                 coul_ijkq = eri_eigenstate_i(jstate,kstate,qstate,pqispin) 
               endif
             endif
#endif

             !$acc loop seq
             do iomega=-se%nomega,se%nomega
               omega = energy(qstate,pqispin) + se%omega(iomega)

               fact_real   = REAL( fact_occ1 / (omega-ei+ej-ek+ieta) + fact_occ2 / (omega-ei+ej-ek-ieta) , dp)
               fact_energy = REAL( fact_occ1 / (energy(pstate,pqispin)-ei+ej-ek+ieta) , dp )

               selfenergy_ring(iomega,pstate,pqispin) = selfenergy_ring(iomega,pstate,pqispin) &
                        + fact_real * coul_ipkj * coul_iqjk * spin_fact

               if(iomega==0 .AND. occupation(pstate,pqispin)>completely_empty) then
                 emp2_ring = emp2_ring + occupation(pstate,pqispin) &
                                       * fact_energy * coul_ipkj * coul_iqjk * spin_fact
               endif

               if( pqispin == jkspin ) then

                 selfenergy_sox(iomega,pstate,pqispin) = selfenergy_sox(iomega,pstate,pqispin) &
                          - fact_real * coul_ipkj * coul_ijkq

                 if(iomega==0 .AND. occupation(pstate,pqispin)>completely_empty) then
                   emp2_sox = emp2_sox - occupation(pstate,pqispin) &
                             * fact_energy * coul_ipkj * coul_ijkq
                 endif

               endif
  
  
             enddo ! iomega

           enddo
         enddo
       enddo
     enddo
#if !defined (_OPENMP) && defined (_OPENACC)
 !$acc end parallel
#endif
#if defined (_OPENMP) && !defined (_OPENACC)
 !$OMP END DO
 !$OMP END PARALLEL
#endif
   enddo 
 enddo ! pqispin

#if !defined (_OPENMP) && defined (_OPENACC)
 !$acc end data
#endif

 call xsum_ortho(selfenergy_ring)
 call xsum_ortho(selfenergy_sox)
 call xsum_ortho(emp2_ring)
 call xsum_ortho(emp2_sox)

 emp2_ring = 0.5_dp * emp2_ring
 emp2_sox  = 0.5_dp * emp2_sox

 if( selfenergy_approx == ONE_RING ) then
   emp2_sox = 0.0_dp
   selfenergy_sox(:,:,:) = 0.0_dp
 endif
 if( selfenergy_approx == SOX ) then
   emp2_ring = 0.0_dp
   selfenergy_ring(:,:,:) = 0.0_dp
 endif

 if( nsemin <= ncore_G+1 .AND. nsemax >= nhomo_G ) then
   emp2 = emp2_ring + emp2_sox
   write(stdout,'(/,a)')       ' MP2 Energy'
   write(stdout,'(a,f14.8)')   ' 2-ring diagram  :',emp2_ring
   write(stdout,'(a,f14.8)')   ' SOX diagram     :',emp2_sox
   write(stdout,'(a,f14.8,/)') ' MP2 correlation :',emp2
 else
   emp2 = 0.0_dp
 endif

! ymbyun 2019/09/25
! NOTE: A performance test is needed.

!$OMP PARALLEL
!$OMP WORKSHARE
 se%sigma(:,:,:) = selfenergy_ring(:,:,:) + selfenergy_sox(:,:,:)
!$OMP END WORKSHARE
!$OMP END PARALLEL

 if( ALLOCATED(eri_eigenstate_i) ) deallocate(eri_eigenstate_i)
 deallocate(selfenergy_ring)
 deallocate(selfenergy_sox)
 if(has_auxil_basis) call destroy_eri_3center_eigen()

 call stop_clock(timing_mp2_self)

end subroutine pt2_selfenergy


!=========================================================================
subroutine onering_selfenergy(selfenergy_approx,nstate,basis,occupation,energy,c_matrix,se,emp2)
 use m_definitions
 use m_mpi
 use m_warning
 use m_basis_set
 use m_eri_ao_mo
 use m_inputparam
 use m_spectral_function
 use m_selfenergy_tools
 implicit none

 integer,intent(in)         :: selfenergy_approx,nstate
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: occupation(nstate,nspin),energy(nstate,nspin)
 real(dp),intent(in)        :: c_matrix(basis%nbf,nstate,nspin)
 type(selfenergy_grid),intent(inout) :: se
 real(dp),intent(out)       :: emp2
!=====
 type(spectral_function) :: vchi0v
 integer                 :: jstate,bstate,jbspin,t_jb
!=====

 call start_clock(timing_mp2_self)

 if( .NOT. has_auxil_basis ) &
   call die('onering_selfenergy: only implemented when an auxiliary basis is available')

 emp2 = 0.0_dp


 write(stdout,'(/,a)') ' Perform the one-ring self-energy calculation'
 write(stdout,*) 'with the perturbative approach'

 call init_spectral_function(nstate,occupation,0,vchi0v)

 call polarizability_onering(basis,nstate,occupation,energy,c_matrix,vchi0v)

#ifdef HAVE_SCALAPACK
 call gw_selfenergy_scalapack(ONE_RING,nstate,basis,occupation,energy,c_matrix,vchi0v,se)
#else
 call gw_selfenergy(ONE_RING,nstate,basis,occupation,energy,c_matrix,vchi0v,se,emp2)
#endif
 
 call destroy_spectral_function(vchi0v)

 call stop_clock(timing_mp2_self)


end subroutine onering_selfenergy


!=========================================================================
subroutine pt2_selfenergy_qs(nstate,basis,occupation,energy,c_matrix,s_matrix,selfenergy,emp2)
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
 use m_mpi
 use m_warning
 use m_basis_set
 use m_eri_ao_mo
 use m_inputparam
 use m_selfenergy_tools
 implicit none

 integer,intent(in)         :: nstate
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: occupation(nstate,nspin),energy(nstate,nspin)
 real(dp),intent(in)        :: c_matrix(basis%nbf,nstate,nspin)
 real(dp),intent(in)        :: s_matrix(basis%nbf,basis%nbf)
 real(dp),intent(out)       :: selfenergy(basis%nbf,basis%nbf,nspin)
 real(dp),intent(out)       :: emp2
!=====
 integer               :: pstate,qstate,qstate2
 real(dp),allocatable  :: selfenergy_ring(:,:,:)
 real(dp),allocatable  :: selfenergy_sox(:,:,:)
 integer               :: istate,jstate,kstate
 integer               :: pqispin,jkspin
 real(dp)              :: fact_occ1,fact_occ2
 real(dp)              :: fi,fj,fk,ei,ej,ek,ep,eq
 real(dp)              :: fact_real,fact_energy
 real(dp)              :: emp2_sox,emp2_ring
 real(dp),allocatable  :: eri_eigenstate_i(:,:,:,:)
 real(dp)              :: coul_iqjk,coul_ijkq,coul_ipkj
! ymbyun 2021/12/27
#ifdef _OPENACC
#if _OPENACC >= 201711
 integer                  :: device_num
 integer(acc_device_kind) :: device_type
 integer(kind=8)          :: gpu_memory   ,gpu_free_memory
 integer                  :: gpu_memory_mb,gpu_free_memory_mb
#endif
#endif
! ymbyun 2023/08/30
! NOTE: It seems like omp_get_wtime() is more stable than cpu_time() and system_clock() for large calculations.
! NOTE: Time and clock routines are compiler specific.
#ifdef ENABLE_TIME
 real(dp)              :: time_all_start,time_all_end
 real(dp)              :: time_eri_start,time_eri_end
 real(dp)              :: time_acc_start,time_acc_end
 integer(kind=8)       :: clock_all_start,clock_all_end,clock_all_rate
 integer(kind=8)       :: clock_eri_start,clock_eri_end,clock_eri_rate
 integer(kind=8)       :: clock_acc_start,clock_acc_end,clock_acc_rate
#ifdef _OPENMP
 real(dp)              :: omp_all_start,omp_all_end
 real(dp)              :: omp_eri_start,omp_eri_end
 real(dp)              :: omp_acc_start,omp_acc_end
#endif
#endif
!=====

 call start_clock(timing_mp2_self)

! ymbyun 2023/08/30
#ifdef ENABLE_TIME
 call cpu_time(time_all_start)
 call system_clock(count_rate=clock_all_rate)
 call system_clock(count=clock_all_start)
#ifdef _OPENMP
 omp_all_start = omp_get_wtime()
#endif
#endif

 emp2_ring = 0.0_dp
 emp2_sox  = 0.0_dp


 write(stdout,'(/,a)') ' Perform the second-order self-energy calculation'
 write(stdout,*) 'with the QP self-consistent approach'

 

 if(has_auxil_basis) then
   call calculate_eri_3center_eigen(basis%nbf,nstate,c_matrix,ncore_G+1,nvirtual_G-1,ncore_G+1,nvirtual_G-1)
 else
   allocate(eri_eigenstate_i(nstate,nstate,nstate,nspin))
 endif



 allocate(selfenergy_ring (nsemin:nsemax,nsemin:nsemax,nspin))
 allocate(selfenergy_sox  (nsemin:nsemax,nsemin:nsemax,nspin))


 selfenergy_ring(:,:,:) = 0.0_dp
 selfenergy_sox(:,:,:)  = 0.0_dp

! ymbyun 2021/12/07
! OpenACC is added.
! NOTE: DATA can be either here or at the outermost loop.
!       Unlike chi_to_vchiv() in linear_response.f90, the position of DATA has no effect on the performance.

#if !defined (_OPENMP) && defined (_OPENACC)
 !$acc data copyin(occupation,energy,eri_eigenstate_i) copy(emp2_ring,emp2_sox,selfenergy_ring,selfenergy_sox)
#endif

#ifdef _OPENACC
#if _OPENACC >= 201711
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
#endif
#endif

 do pqispin=1,nspin
   do istate=ncore_G+1,nvirtual_G-1 !LOOP of the first Green's function
     if( MODULO( istate - (ncore_G+1) , nproc_ortho ) /= rank_ortho ) cycle

     if( .NOT. has_auxil_basis ) then
! ymbyun 2023/08/30
#ifdef ENABLE_TIME
       call cpu_time(time_eri_start)
       call system_clock(count_rate=clock_eri_rate)
       call system_clock(count=clock_eri_start)
#ifdef _OPENMP
       omp_eri_start = omp_get_wtime()
#endif
#endif

       call calculate_eri_4center_eigen(basis%nbf,nstate,c_matrix,istate,pqispin,eri_eigenstate_i)

! ymbyun 2023/08/30
#ifdef ENABLE_TIME
       call cpu_time(time_eri_end)
       print '("(ymbyun) pt2_selfenergy_qs (time_eri)  : ",f30.6," sec")',time_eri_end-time_eri_start

       call system_clock(count=clock_eri_end)
       print '("(ymbyun) pt2_selfenergy_qs (clock_eri) : ",f30.6," sec")',(clock_eri_end-clock_eri_start)/real(clock_eri_rate)

#ifdef _OPENMP
       omp_eri_end = omp_get_wtime()
       print '("(ymbyun) pt2_selfenergy_qs (omp_eri)   : ",f30.6," sec")',omp_eri_end-omp_eri_start
#endif
#endif
     endif

     fi = occupation(istate,pqispin)
     ei = energy(istate,pqispin)

! ymbyun 2019/09/23
! NOTE: OpenMP 2.0 (2000) supports array reductions for Fortran, while OpenMP 4.5 (2015) supports array reductions for C/C++.

! ymbyun 2020/10/19
! NOTE: The following thing is changed to make OpenMP work for KISTI Nurion with Intel Xeon Phi CPUs and Intel compilers 18.0.3.
! NOTE: 1. Array reduction is removed.

! ymbyun 2021/12/07
! OpenACC is added.
! NOTE: Reduction for emp2_ring and emp2_sox is implicitly generated.
! NOTE: COLLAPSE(4) works well for OpenMP, but causes an error for OpenACC.
! NOTE: Maybe, COLLAPSE(3) should be used for OpenMP for safety purposes.

! ymbyun 2021/12/13
! NOTE: For the fair performance comparsion between OpenMP (CPU) and OpenACC (GPU), COLLAPSE(4) in OpenMP is changed to COLLAPSE(3).

! ymbyun 2023/07/25
! NOTE: Due to open-shell systems, COLLAPSE(3) in OpenMP and OpenACC is changed to COLLAPSE(2).
! NOTE: COLLAPSE works differently for OpenMP and OpenACC [i.e. COLLAPSE(4) in OpenMP gives the same result as COLLAPSE(2) in OpenACC].

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
 !$OMP DO PRIVATE(pstate, qstate, jkspin, jstate,fj,ej, kstate,fk,ek,fact_occ1,fact_occ2,coul_ipkj,coul_iqjk,coul_ijkq,ep,eq,fact_real,fact_energy) REDUCTION(+:emp2_ring,emp2_sox) COLLAPSE(2)
#endif
#if !defined (_OPENMP) && defined (_OPENACC)
 !$acc parallel loop independent collapse(2)
#endif

!!REDUCTION(+:emp2_ring,emp2_sox,selfenergy_ring,selfenergy_sox) COLLAPSE(4)

     do pstate=nsemin,nsemax ! external loop ( bra )
       do qstate=nsemin,nsemax   ! external loop ( ket )

         !$acc loop seq
         do jkspin=1,nspin
           !$acc loop seq
           do jstate=ncore_G+1,nvirtual_G-1  !LOOP of the second Green's function
             fj = occupation(jstate,jkspin)
             ej = energy(jstate,jkspin)

             !$acc loop seq
             do kstate=ncore_G+1,nvirtual_G-1 !LOOP of the third Green's function
               fk = occupation(kstate,jkspin)
               ek = energy(kstate,jkspin)
  
               fact_occ1 = (spin_fact-fi) *            fj  * (spin_fact-fk) / spin_fact**3
               fact_occ2 =            fi  * (spin_fact-fj) *            fk  / spin_fact**3
 
               if( fact_occ1 < completely_empty .AND. fact_occ2 < completely_empty ) cycle

! ymbyun 2021/12/07
! This is an ugly way to avoid accessing eri_eigen_ri().
! Due to this, be_qspt2_ri.in in the test suite fails.

#ifdef _OPENACC
               if( has_auxil_basis ) then
!                 call die('(ymbyun) Currently, OpenACC cannot be used together with RI and MPI')
                 write(stdout,*) '(ymbyun) Currently, OpenACC cannot be used together with RI and MPI'
                 stop
               else
                 coul_ipkj = eri_eigenstate_i(pstate,kstate,jstate,jkspin)
                 coul_iqjk = eri_eigenstate_i(qstate,jstate,kstate,jkspin)
                 if( pqispin == jkspin ) then
                   coul_ijkq = eri_eigenstate_i(jstate,kstate,qstate,pqispin) 
                 endif
               endif
#else
               if( has_auxil_basis ) then
                 ! ymbyun 2019/10/01
                 ! NOTE: This part is a bottleneck.
                 ! NOTE: We should minimize the number of calls to MPI_Allreduce at this part using mpi_nproc_ortho.
                 coul_ipkj = eri_eigen_ri(istate,pstate,pqispin,kstate,jstate,jkspin)
                 coul_iqjk = eri_eigen_ri(istate,qstate,pqispin,jstate,kstate,jkspin)
                 if( pqispin == jkspin ) then
                   coul_ijkq = eri_eigen_ri(istate,jstate,pqispin,kstate,qstate,pqispin) 
                 endif
               else
                 coul_ipkj = eri_eigenstate_i(pstate,kstate,jstate,jkspin)
                 coul_iqjk = eri_eigenstate_i(qstate,jstate,kstate,jkspin)
                 if( pqispin == jkspin ) then
                   coul_ijkq = eri_eigenstate_i(jstate,kstate,qstate,pqispin) 
                 endif
               endif
#endif

               ep = energy(pstate,pqispin) 
               eq = energy(qstate,pqispin) 
  
               fact_real   = REAL( fact_occ1 / ( eq - ei + ej - ek + ieta) &
                                 + fact_occ2 / ( eq - ei + ej - ek - ieta) , dp)
               fact_energy = REAL( fact_occ1 / ( ep - ei + ej - ek + ieta) , dp )
  
               selfenergy_ring(pstate,qstate,pqispin) = selfenergy_ring(pstate,qstate,pqispin) &
                        + fact_real * coul_ipkj * coul_iqjk * spin_fact

               if(pstate==qstate .AND. occupation(pstate,pqispin)>completely_empty) then
                 emp2_ring = emp2_ring + occupation(pstate,pqispin) &
                                       * fact_energy * coul_ipkj * coul_iqjk * spin_fact
               endif
 
               if( pqispin == jkspin ) then

                 selfenergy_sox(pstate,qstate,pqispin) = selfenergy_sox(pstate,qstate,pqispin) &
                          - fact_real * coul_ipkj * coul_ijkq

                 if(pstate==qstate .AND. occupation(pstate,pqispin)>completely_empty) then
                   emp2_sox = emp2_sox - occupation(pstate,pqispin) &
                             * fact_energy * coul_ipkj * coul_ijkq
                 endif

               endif
    
    
    
             enddo
           enddo
         enddo
       enddo 
     enddo
#if !defined (_OPENMP) && defined (_OPENACC)
 !$acc end parallel
#endif
#if defined (_OPENMP) && !defined (_OPENACC)
 !$OMP END DO
 !$OMP END PARALLEL
#endif

! ymbyun 2023/08/30
#ifdef ENABLE_TIME
     call cpu_time(time_acc_end)
     print '("(ymbyun) pt2_selfenergy_qs (time_acc)  : ",f30.6," sec")',time_acc_end-time_acc_start

     call system_clock(count=clock_acc_end)
     print '("(ymbyun) pt2_selfenergy_qs (clock_acc) : ",f30.6," sec")',(clock_acc_end-clock_acc_start)/real(clock_acc_rate)

#ifdef _OPENMP
     omp_acc_end = omp_get_wtime()
     print '("(ymbyun) pt2_selfenergy_qs (omp_acc)   : ",f30.6," sec")',omp_acc_end-omp_acc_start
#endif
#endif

   enddo
 enddo ! pqispin
#if !defined (_OPENMP) && defined (_OPENACC)
 !$acc end data
#endif

 ! ymbyun 2019/10/06
 ! NOTE: This part is a bottleneck.
 ! NOTE: Currently, no-RI (OpenMP) is faster than RI (MPI).
 ! NOTE: For intra-node MPI communication, Ethernet (MPICH) is faster than InfiniBand (MVAPICH2) and Aries (Cray MPICH).
 ! TODO: Tests for communication time vs computation time
 ! TODO: Tests for inter-node MPI communication
 ! TODO: Tests for MV2_SHMEM_ALLREDUCE_MSG and MV2_ALLREDUCE_2LEVEL_MSG in MVAPICH2
 ! TODO: Tests for MPI-3 shared memory
 ! TODO: Tests for OpenMPI
 call xsum_ortho(selfenergy_ring)
 call xsum_ortho(selfenergy_sox)
 call xsum_ortho(emp2_ring)
 call xsum_ortho(emp2_sox)

 emp2_ring = 0.5_dp * emp2_ring
 emp2_sox  = 0.5_dp * emp2_sox
 if( nsemin <= ncore_G+1 .AND. nsemax >= nhomo_G ) then
   emp2 = emp2_ring + emp2_sox
   write(stdout,'(/,a)')       ' MP2 Energy'
   write(stdout,'(a,f14.8)')   ' 2-ring diagram  :',emp2_ring
   write(stdout,'(a,f14.8)')   ' SOX diagram     :',emp2_sox
   write(stdout,'(a,f14.8,/)') ' MP2 correlation :',emp2
 else
   emp2 = 0.0_dp
 endif

! ymbyun 2019/09/25
! NOTE: A performance test is needed.

!$OMP PARALLEL
!$OMP WORKSHARE
 selfenergy(:,:,:) = selfenergy_ring(:,:,:) + selfenergy_sox(:,:,:)
!$OMP END WORKSHARE
!$OMP END PARALLEL

 call apply_qs_approximation(basis%nbf,nstate,s_matrix,c_matrix,selfenergy)


 if( ALLOCATED(eri_eigenstate_i) ) deallocate(eri_eigenstate_i)
 deallocate(selfenergy_ring)
 deallocate(selfenergy_sox)
 if(has_auxil_basis) call destroy_eri_3center_eigen()

 call stop_clock(timing_mp2_self)

! ymbyun 2023/08/30
#ifdef ENABLE_TIME
 call cpu_time(time_all_end)
 print '("(ymbyun) pt2_selfenergy_qs (time_all)  : ",f30.6," sec")',time_all_end-time_all_start

 call system_clock(count=clock_all_end)
 print '("(ymbyun) pt2_selfenergy_qs (clock_all) : ",f30.6," sec")',(clock_all_end-clock_all_start)/real(clock_all_rate)

#ifdef _OPENMP
 omp_all_end = omp_get_wtime()
 print '("(ymbyun) pt2_selfenergy_qs (omp_all)   : ",f30.6," sec")',omp_all_end-omp_all_start
#endif
#endif

end subroutine pt2_selfenergy_qs


!=========================================================================
