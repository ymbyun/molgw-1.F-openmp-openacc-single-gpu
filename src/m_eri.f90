!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the method to prepare and store the 2-, 3-, and 4-center Coulomb integrals
!
!=========================================================================
module m_eri
 use m_definitions
 use m_mpi
 use m_memory
 use m_warning
 use m_basis_set
 use m_timing
 use m_cart_to_pure
 use m_libint_tools


 real(dp),parameter,public :: TOO_LOW_EIGENVAL=1.0e-6_dp

 !
 ! max length of a record in the ERI file
 integer,parameter,private :: line_length=1000

 real(dp),private           :: TOL_INT


 real(dp),public,allocatable :: eri_4center(:)
 real(dp),public,allocatable :: eri_4center_lr(:)
 real(dp),public,allocatable :: eri_3center(:,:)
 real(dp),public,allocatable :: eri_3center_lr(:,:)


 logical,protected,allocatable      :: negligible_shellpair(:,:)
! ymbyun 2018/05/03
#ifdef ENABLE_YMBYUN
 integer(kind=8),private  ,allocatable      :: index_pair_1d(:)
#else
 integer,private  ,allocatable     :: index_pair_1d(:)
#endif
 integer,protected,allocatable      :: index_basis(:,:)
 integer,protected,allocatable      :: index_shellpair(:,:)
 integer,protected                  :: nshellpair

 integer,private,allocatable        :: shell_bf(:)


 integer,private   :: nbf_eri            ! local copy of nbf
! ymbyun 2018/05/03
! i4 is too small to handle the 4-center Coulomb integrals using AE 5Z without RI for transition metals.
! Eventually, MOLGW will have to use i8 by default like NWChem.
#ifdef ENABLE_YMBYUN
 integer(kind=8),protected :: nsize      ! size of the eri_4center array
 integer(kind=8),protected :: npair      ! number of independent pairs (i,j) with i<=j 
#else
 integer,protected :: nsize             ! size of the eri_4center array
 integer,protected :: npair             ! number of independent pairs (i,j) with i<=j 
#endif

! ymbyun 2021/11/30
! Scalars
!$acc declare create(npair)
!$acc declare create(nbf_eri)

! Arrays using allocate()
!$acc declare create(negligible_shellpair)
!$acc declare create(shell_bf)

! Arrays using clean_allocate(). Special care is needed.
!$acc declare create(eri_4center)
!$acc declare create(index_pair_1d)

 integer,protected :: nauxil_3center     ! size of the 3-center matrix
                                         ! may differ from the total number of 3-center integrals due to
                                         ! data distribution
 integer,protected :: nauxil_3center_lr  ! size of the 3-center matrix
                                         ! may differ from the total number of 3-center integrals due to
                                         ! data distribution

 real(dp),allocatable,public  :: eri_3center_sca(:,:)

! Parallelization information for the auxiliary basis
 integer,allocatable,protected :: iproc_ibf_auxil(:)
 integer,allocatable,protected :: ibf_auxil_g(:)       ! auxil bf index from local to global
 integer,allocatable,protected :: ibf_auxil_l(:)       ! auxil bf index from global to local

! Parallelization information for the auxiliary basis (LR part)
 integer,allocatable,protected :: iproc_ibf_auxil_lr(:)
 integer,allocatable,protected :: ibf_auxil_g_lr(:)
 integer,allocatable,protected :: ibf_auxil_l_lr(:)


contains


!=========================================================================
subroutine prepare_eri(basis)
 use m_inputparam,only: integral_level
 implicit none
!===== 
 type(basis_set),intent(in) :: basis
!===== 
 logical            :: file_exists
!===== 

 nbf_eri = basis%nbf

! ymbyun 2021/12/06
!$acc update device(nbf_eri)

 select case(integral_level)
 case(low)       ! accuracy not guaranted, just for quick test runs
   TOL_INT = 1.0e-04_dp
 case(medium)    ! 10 meV accuracy on potentials
   TOL_INT = 1.0e-06_dp
 case(high)      !  1 meV accuracy on potentials
   TOL_INT = 1.0e-08_dp
 case(very_high) ! almost perfect potentials
   TOL_INT = 1.0e-10_dp
 case(insane)    ! No screening of any integral
   TOL_INT = 0.0_dp
 case default
   call die('integration quality not recognized')
 end select
 write(stdout,'(/,a,es9.2)') ' Tolerance on integrals set to ',TOL_INT


 if(.NOT.ALLOCATED(negligible_shellpair)) then
   call setup_shell_index(basis)
   allocate(negligible_shellpair(basis%nshell,basis%nshell))
   call identify_negligible_shellpair(basis)
   call setup_shellpair(basis)
   call setup_basispair()
 endif

 nsize = (npair*(npair+1))/2

end subroutine prepare_eri


!=========================================================================
subroutine deallocate_eri_4center()
 implicit none
!=====

 if(ALLOCATED(eri_4center)) then
! ymbyun 2021/12/06
! Currently, NVIDIA (PGI) compilers 21.9 do not understand clean_allocate(), leading to FATAL ERROR.
#ifdef _OPENACC
   deallocate(eri_4center)
#else
   call clean_deallocate('4-center integrals',eri_4center)
#endif
 endif

end subroutine deallocate_eri_4center


!=========================================================================
subroutine deallocate_eri_4center_lr()
 implicit none
!=====

 if(ALLOCATED(eri_4center_lr)) then
   call clean_deallocate('4-center LR integrals',eri_4center_lr)
 endif

end subroutine deallocate_eri_4center_lr


!=========================================================================
subroutine deallocate_eri()
 implicit none

 integer :: ishell
!=====

 if(ALLOCATED(eri_4center)) then
! ymbyun 2021/12/06
! Currently, NVIDIA (PGI) compilers 21.9 do not understand clean_allocate(), leading to FATAL ERROR.
#ifdef _OPENACC
   deallocate(eri_4center)
#else
   call clean_deallocate('4-center integrals',eri_4center)
#endif
 endif
 if(ALLOCATED(eri_4center_lr)) then
   call clean_deallocate('4-center LR integrals',eri_4center_lr)
 endif
 if(ALLOCATED(negligible_shellpair))   deallocate(negligible_shellpair)
 if(ALLOCATED(index_shellpair))        deallocate(index_shellpair)
! ymbyun 2021/12/06
! Currently, NVIDIA (PGI) compilers 21.9 do not understand clean_allocate(), leading to FATAL ERROR.
#ifdef _OPENACC
 if(ALLOCATED(index_pair_1d))          deallocate(index_pair_1d)
#else
 if(ALLOCATED(index_pair_1d)) call clean_deallocate('index pair',index_pair_1d)
#endif
 if(ALLOCATED(index_basis))   call clean_deallocate('index basis',index_basis)

 if(ALLOCATED(shell_bf))              deallocate(shell_bf)


end subroutine deallocate_eri


!=========================================================================
subroutine deallocate_index_pair()
 implicit none

!=====

! ymbyun 2021/12/06
! Currently, NVIDIA (PGI) compilers 21.9 do not understand clean_allocate(), leading to FATAL ERROR.
#ifdef _OPENACC
 if(ALLOCATED(index_pair_1d)) deallocate(index_pair_1d)
#else
 if(ALLOCATED(index_pair_1d)) call clean_deallocate('index pair',index_pair_1d)
#endif

end subroutine deallocate_index_pair


!=========================================================================
function index_eri(ibf,jbf,kbf,lbf)
! ymbyun 2021/11/30
!$acc routine seq

 implicit none

 integer,intent(in) :: ibf,jbf,kbf,lbf
! ymbyun 2018/05/04
#ifdef ENABLE_YMBYUN
 integer(kind=8)    :: index_eri
!=====
 integer(kind=8)    :: klmin,ijmax
 integer(kind=8)    :: index_ij,index_kl
#else
 integer            :: index_eri
!=====
 integer            :: klmin,ijmax
 integer            :: index_ij,index_kl
#endif
!===== 

 index_ij = index_pair(ibf,jbf)
 index_kl = index_pair(kbf,lbf)

 ijmax=MAX(index_ij,index_kl)
 klmin=MIN(index_ij,index_kl)

 index_eri = (klmin-1)*npair - (klmin-1)*(klmin-2)/2 + ijmax-klmin+1

end function index_eri


!=========================================================================
function index_pair(ibf,jbf)
! ymbyun 2021/11/30
!$acc routine seq

 implicit none

 integer,intent(in) :: ibf,jbf
! ymbyun 2018/05/05
#ifdef ENABLE_YMBYUN
 integer(kind=8)    :: index_pair
#else
 integer            :: index_pair
#endif
!=====
 integer            :: ijmin,ijmax
!=====

 if( ibf == jbf ) then
   index_pair = ibf
 else
   ijmax=MAX(ibf,jbf)
   ijmin=MIN(ibf,jbf)

   index_pair = (ijmin-1) * (nbf_eri-1) - (ijmin-1) * (ijmin-2)/2     + ijmax - ijmin + nbf_eri
   index_pair = index_pair_1d(index_pair)
 endif


end function index_pair


!=========================================================================
function eri(ibf,jbf,kbf,lbf)
! ymbyun 2021/11/30
!$acc routine seq

 implicit none
 integer,intent(in) :: ibf,jbf,kbf,lbf
 real(dp)           :: eri
!=====

 if( negligible_basispair(ibf,jbf) .OR. negligible_basispair(kbf,lbf) ) then
   eri = 0.0_dp
 else
   eri = eri_4center(index_eri(ibf,jbf,kbf,lbf))
 endif

end function eri


!=========================================================================
function eri_lr(ibf,jbf,kbf,lbf)
 implicit none
 integer,intent(in) :: ibf,jbf,kbf,lbf
 real(dp)           :: eri_lr
!=====

 if( negligible_basispair(ibf,jbf) .OR. negligible_basispair(kbf,lbf) ) then
   eri_lr = 0.0_dp
 else
   eri_lr = eri_4center_lr(index_eri(ibf,jbf,kbf,lbf))
 endif

end function eri_lr


!=========================================================================
function eri_ri(ibf,jbf,kbf,lbf)
 implicit none
 integer,intent(in) :: ibf,jbf,kbf,lbf
 real(dp)           :: eri_ri
!=====
! ymbyun 2018/05/05
#ifdef ENABLE_YMBYUN
 integer(kind=8)    :: index_ij,index_kl
#else
 integer            :: index_ij,index_kl
#endif
! real(dp)           :: eri_1(1,1)
!=====

 if( negligible_basispair(ibf,jbf) .OR. negligible_basispair(kbf,lbf) ) then
   eri_ri = 0.0_dp
 else
   index_ij = index_pair(ibf,jbf)
   index_kl = index_pair(kbf,lbf)
  
   eri_ri = DOT_PRODUCT( eri_3center(:,index_ij) , eri_3center(:,index_kl) )

   call xsum_auxil(eri_ri)

 endif

end function eri_ri


!=========================================================================
function eri_ri_lr(ibf,jbf,kbf,lbf)
 implicit none
 integer,intent(in) :: ibf,jbf,kbf,lbf
 real(dp)           :: eri_ri_lr
!=====
! ymbyun 2018/05/04
#ifdef ENABLE_YMBYUN
 integer(kind=8)    :: index_ij,index_kl
#else
 integer            :: index_ij,index_kl
#endif
!=====

 if( negligible_basispair(ibf,jbf) .OR. negligible_basispair(kbf,lbf) ) then
   eri_ri_lr = 0.0_dp
 else
   index_ij = index_pair(ibf,jbf)
   index_kl = index_pair(kbf,lbf)

   eri_ri_lr = DOT_PRODUCT( eri_3center_lr(:,index_ij) , eri_3center_lr(:,index_kl) )

   call xsum_auxil(eri_ri_lr)

 endif

end function eri_ri_lr


!=========================================================================
subroutine setup_shell_index(basis)
 implicit none

 type(basis_set),intent(in)   :: basis
!=====
 integer :: ibf
!=====

 allocate(shell_bf(basis%nbf))
 shell_bf(:) = basis%bff(:)%shell_index

! ymbyun 2021/12/06
!$acc update device(shell_bf)

end subroutine setup_shell_index


!=========================================================================
subroutine setup_basispair()
 implicit none
!=====
 integer :: ishell,jshell
 integer :: ibf,jbf,ijbf
!=====

 npair = 0
 do jbf=1,nbf_eri
   do ibf=1,jbf
     if( .NOT. negligible_basispair(ibf,jbf) ) then
       npair = npair + 1
     endif
   enddo
 enddo

! ymbyun 2021/12/06
! Currently, NVIDIA (PGI) compilers 21.9 do not understand clean_allocate(), leading to FATAL ERROR.
#ifdef _OPENACC
 allocate(index_pair_1d((nbf_eri*(nbf_eri+1))/2))
#else
 call clean_allocate('index pair',index_pair_1d,(nbf_eri*(nbf_eri+1))/2)
#endif
 call clean_allocate('index basis',index_basis,2,npair)

 !
 ! Specific ordering where the first nbf pairs contain the diagonal terms ibf==jbf
 !
 npair = 0
 index_pair_1d(:) = 0
 do jbf=1,nbf_eri
   if( negligible_basispair(jbf,jbf) ) then
     call die('setup_negligible_basispair: this should not happen')
   endif
   npair = npair + 1
   index_pair_1d(jbf) = npair
   index_pair_1d(jbf) = npair
   index_basis(1,npair) = jbf
   index_basis(2,npair) = jbf
 enddo

 ijbf = nbf_eri

 do ibf=1,nbf_eri 
   do jbf=ibf+1,nbf_eri  ! Skip the diagonal terms since it is already included 
     ijbf = ijbf + 1
     if( .NOT. negligible_basispair(ibf,jbf) ) then
       npair = npair + 1
       index_pair_1d(ijbf) = npair
       index_basis(1,npair) = ibf
       index_basis(2,npair) = jbf
     endif
   enddo
 enddo

! ymbyun 2021/12/06
!$acc update device(npair)
!$acc update device(index_pair_1d)

end subroutine setup_basispair


!=========================================================================
function negligible_basispair(ibf,jbf)
! ymbyun 2021/11/30
!$acc routine seq

 implicit none

 integer,intent(in) :: ibf,jbf
 logical :: negligible_basispair
!=====
 integer  :: ishell,jshell
!=====

 ishell = shell_bf(ibf)
 jshell = shell_bf(jbf)

 negligible_basispair = negligible_shellpair(ishell,jshell)

end function negligible_basispair


!=========================================================================
!
! Find negligible shell pairs with
! Cauchy-Schwarz inequality: (ij|1/r|kl)**2 <= (ij|1/r|ij) (kl|1/r|(kl) 
!
!=========================================================================
subroutine identify_negligible_shellpair(basis)
 implicit none

 type(basis_set),intent(in)   :: basis
!=====
 integer                      :: info,ip
 integer                      :: ibf,jbf
 integer                      :: n1c,n2c
 integer                      :: ni,nj
 integer                      :: ami,amj
 integer                      :: ishell,jshell
 real(dp),allocatable         :: integrals(:,:,:,:)
 real(dp)                     :: workload(nproc_world)
 integer                      :: shell_proc(basis%nshell)
!=====
! variables used to call C
 integer(C_INT)               :: am1,am2
 integer(C_INT)               :: ng1,ng2
 real(C_DOUBLE),allocatable   :: alpha1(:),alpha2(:)
 real(C_DOUBLE)               :: x01(3),x02(3)
 real(C_DOUBLE),allocatable   :: coeff1(:),coeff2(:)
 real(C_DOUBLE),allocatable   :: int_shell(:)
!=====

 call start_clock(timing_eri_screening)
 write(stdout,'(/,a)')    ' Cauchy-Schwartz screening of the 3- or 4-center integrals'

 !
 ! Load balancing
 workload(:) = 0.0_dp
 do jshell=1,basis%nshell
   amj = basis%shell(jshell)%am
   ip = MINLOC(workload(:),DIM=1)
   !
   ! Cost function was evaluated from a few runs
   workload(ip) = workload(ip) + cost_function_eri(amj)
   shell_proc(jshell) = ip - 1
 enddo


 negligible_shellpair(:,:) = .TRUE.

 do jshell=1,basis%nshell
   !
   ! Workload is distributed here
   if( shell_proc(jshell) /= rank_world ) cycle

   amj = basis%shell(jshell)%am
   nj  = number_basis_function_am( basis%gaussian_type , amj )
   n2c = number_basis_function_am( 'CART' , amj )
   am2 = basis%shell(jshell)%am
   ng2 = basis%shell(jshell)%ng

   do ishell=1,basis%nshell
     ami = basis%shell(ishell)%am
     if( ami < amj ) cycle

     ni = number_basis_function_am( basis%gaussian_type , ami )
     n1c = number_basis_function_am( 'CART' , ami )
     am1 = basis%shell(ishell)%am
     ng1 = basis%shell(ishell)%ng

     allocate(alpha1(ng1),alpha2(ng2))
     allocate(coeff1(ng1),coeff2(ng2))
     alpha1(:) = basis%shell(ishell)%alpha(:)
     alpha2(:) = basis%shell(jshell)%alpha(:)
     x01(:) = basis%shell(ishell)%x0(:)
     x02(:) = basis%shell(jshell)%x0(:)
     coeff1(:) = basis%shell(ishell)%coeff(:)
     coeff2(:) = basis%shell(jshell)%coeff(:)

     allocate( int_shell( n1c*n2c*n1c*n2c ) )


     call libint_4center(am1,ng1,x01,alpha1,coeff1, &
                         am2,ng2,x02,alpha2,coeff2, &
                         am1,ng1,x01,alpha1,coeff1, &
                         am2,ng2,x02,alpha2,coeff2, &
                         0.0_C_DOUBLE,int_shell)

     call transform_libint_to_molgw(basis%gaussian_type,ami,amj,ami,amj,int_shell,integrals)


     do ibf=1,ni
       do jbf=1,nj
         if( ABS( integrals(ibf,jbf,ibf,jbf) ) > TOL_INT**2 ) negligible_shellpair(ishell,jshell) = .FALSE.
       enddo
     enddo

     !
     ! Symmetrize
     negligible_shellpair(jshell,ishell) = negligible_shellpair(ishell,jshell)

     deallocate(integrals)
     deallocate(int_shell)
     deallocate(alpha1,alpha2)
     deallocate(coeff1,coeff2)

   enddo
 enddo

 call xand_world(negligible_shellpair)

 call stop_clock(timing_eri_screening)

! ymbyun 2021/12/06
!$acc update device(negligible_shellpair)

end subroutine identify_negligible_shellpair


!=========================================================================
subroutine setup_shellpair(basis)
 implicit none

 type(basis_set),intent(in) :: basis
!=====
 integer :: ishell,jshell
 integer :: ami,amj
 integer :: ishellpair,jshellpair
!=====

 ishellpair = 0
 jshellpair = 0
 do jshell=1,basis%nshell
   do ishell=1,jshell 
     jshellpair = jshellpair + 1
     ! skip the identified negligible shell pairs
     if( negligible_shellpair(ishell,jshell) ) cycle
     ami = basis%shell(ishell)%am
     amj = basis%shell(jshell)%am
     ishellpair = ishellpair + 1

   enddo
 enddo
 nshellpair = ishellpair
 write(stdout,'(/,1x,a,i8,a,i8)') 'Non negligible shellpairs to be computed',nshellpair,'  over a total of',jshellpair
 write(stdout,'(1x,a,f12.4)')     'Saving (%): ', REAL(jshellpair-nshellpair,dp)/REAL(jshellpair,dp)

 allocate(index_shellpair(2,nshellpair))

 ishellpair = 0
 do jshell=1,basis%nshell
   do ishell=1,jshell 
     ! skip the identified negligible shell pairs
     if( negligible_shellpair(ishell,jshell) ) cycle
     ami = basis%shell(ishell)%am
     amj = basis%shell(jshell)%am
     ishellpair = ishellpair + 1
     ! Reverse if needed the order of the shell so to maximize the angular
     ! momentum of the first shell
     if( ami >= amj ) then
       index_shellpair(1,ishellpair) = ishell
       index_shellpair(2,ishellpair) = jshell
     else
       index_shellpair(1,ishellpair) = jshell
       index_shellpair(2,ishellpair) = ishell
     endif

   enddo
 enddo


end subroutine setup_shellpair


!=================================================================
subroutine destroy_eri_3center()
 implicit none
!=====

 if(ALLOCATED(iproc_ibf_auxil)) then
   deallocate(iproc_ibf_auxil)
 endif
 if(ALLOCATED(ibf_auxil_g)) then
   deallocate(ibf_auxil_g)
 endif
 if(ALLOCATED(ibf_auxil_l)) then
   deallocate(ibf_auxil_l)
 endif
 if(ALLOCATED(eri_3center)) then
   call clean_deallocate('3-center integrals',eri_3center)
 endif

#ifdef SCASCA
 if(ALLOCATED(eri_3center_sca)) then
   call clean_deallocate('3-center integrals SCALAPACK',eri_3center_sca)
 endif
#endif

end subroutine destroy_eri_3center


!=================================================================
subroutine destroy_eri_3center_lr()
 implicit none
!=====

 if(ALLOCATED(iproc_ibf_auxil_lr)) then
   deallocate(iproc_ibf_auxil_lr)
 endif
 if(ALLOCATED(ibf_auxil_g_lr)) then
   deallocate(ibf_auxil_g_lr)
 endif
 if(ALLOCATED(ibf_auxil_l_lr)) then
   deallocate(ibf_auxil_l_lr)
 endif
 if(ALLOCATED(eri_3center_lr)) then
   call clean_deallocate('3-center LR integrals',eri_3center_lr)
 endif

end subroutine destroy_eri_3center_lr


!=========================================================================
subroutine negligible_eri(tol)
 implicit none
 real(dp),intent(in) :: tol
!=====
 integer             :: icount,ibf,jbf,kbf,lbf,jcount
! ymbyun 2018/05/04
#ifdef ENABLE_YMBYUN
 integer(kind=8)     :: ibuffer
#else
 integer             :: ibuffer
#endif
 real(dp)            :: integral_ij(nbf_eri,nbf_eri)
!=====

 icount=0
 do ibuffer=1,nsize
   if( ABS( eri_4center(ibuffer) ) < tol ) icount=icount+1
 enddo

 write(stdout,*) ' number of negligible integrals <',tol
 write(stdout,*) icount, ' / ',nsize,REAL(icount,dp)/REAL(nsize,dp)*100.0_dp,' [%]'


 do ibf=1,nbf_eri
   do jbf=1,nbf_eri
     integral_ij(ibf,jbf) = eri(ibf,jbf,ibf,jbf)
   enddo
 enddo

 write(stdout,*) 'testing Cauchy-Schwarz condition'
 icount=0
 jcount=0
 do ibf=1,nbf_eri
   do jbf=1,nbf_eri
     do kbf=1,nbf_eri
       do lbf=1,nbf_eri
         if( SQRT( integral_ij(ibf,jbf) * integral_ij(kbf,lbf) ) < tol ) icount = icount + 1
         if( ABS( eri(ibf,jbf,kbf,lbf) ) < tol ) jcount = jcount + 1
       enddo
     enddo
   enddo
 enddo
 write(stdout,*) ' number of negligible integrals <',tol
 write(stdout,*) icount, ' / ',nbf_eri**4,REAL(icount,dp)/REAL(nbf_eri,dp)**4*100.0_dp,' [%]'
 write(stdout,*) jcount, ' / ',nbf_eri**4,REAL(jcount,dp)/REAL(nbf_eri,dp)**4*100.0_dp,' [%]'


end subroutine negligible_eri


!=========================================================================
subroutine dump_out_eri(rcut)
 implicit none
 real(dp),intent(in) :: rcut
!=====
 character(len=50) :: filename
! ymbyun 2018/05/04
#ifdef ENABLE_YMBYUN
 integer(kind=8)   :: nline,iline,icurrent
#else
 integer           :: nline,iline,icurrent
#endif
 integer           :: erifile
!=====

 if(rcut < 1.0e-6_dp) then
   filename='molgw_eri.data'
 else
   filename='molgw_eri_lr.data'
 endif
 write(stdout,*) 'Dump out the ERI into file'
 write(stdout,*) 'Size of file [bytes]',REAL(nsize,dp)*dp

 if( is_iomaster ) then
   open(newunit=erifile,file=TRIM(filename),form='unformatted')
   write(erifile) nsize
   write(erifile) rcut

   nline = nsize / line_length + 1
   icurrent=0
   do iline=1,nline
     write(erifile) eri_4center(icurrent+1:MIN(nsize,icurrent+line_length+1))
     icurrent = icurrent + line_length + 1
   enddo

   close(erifile)
 endif

 write(stdout,'(a,/)') ' file written'

end subroutine dump_out_eri


!=========================================================================
logical function read_eri(rcut)
 implicit none
 real(dp),intent(in) :: rcut
!=====
 character(len=50) :: filename
! ymbyun 2018/05/04
#ifdef ENABLE_YMBYUN
 integer(kind=8)   :: nline,iline,icurrent
 integer(kind=8)   :: integer_read
#else
 integer           :: nline,iline,icurrent
 integer           :: integer_read
#endif
 real(dp)          :: real_read
 integer           :: erifile
!=====

 if(rcut < 1.0e-6_dp) then
   filename='molgw_eri.data'
 else
   filename='molgw_eri_lr.data'
 endif
 
 inquire(file=TRIM(filename),exist=read_eri)

 if(read_eri) then

   write(stdout,*) 'Try to read ERI file'
   open(newunit=erifile,file=TRIM(filename),form='unformatted',status='old')
   read(erifile) integer_read
   if(integer_read /= nsize) read_eri=.FALSE.
   read(erifile) real_read
   if(ABS(real_read-rcut) > 1.0e-6_dp) read_eri=.FALSE.

   if(read_eri) then

     nline = nsize / line_length + 1
     icurrent=0
     do iline=1,nline
       read(erifile) eri_4center(icurrent+1:MIN(nsize,icurrent+line_length+1))
       icurrent = icurrent + line_length + 1
     enddo
     write(stdout,'(a,/)') ' ERI file read'

   else
     write(stdout,'(a,/)') ' reading aborted'
   endif

   close(erifile)

 endif


end function read_eri


!=========================================================================
! Rough evaluation of the CPU time to get an ERI as a function of the 
! angular momentum
! 
!=========================================================================
function cost_function_eri(am)
 implicit none
 integer,intent(in)  :: am
 real(dp)            :: cost_function_eri
!=====

 cost_function_eri = am**2 + 4.6_dp 

end function cost_function_eri


!=========================================================================
subroutine distribute_auxil_basis(nbf_auxil_basis)
 use m_scalapack
 implicit none

 integer,intent(in)  :: nbf_auxil_basis
!=====
 integer :: ibf
 integer :: ibf_local
 integer :: iproc
 integer :: ilocal,iglobal
 integer :: nbf_local_iproc(0:nproc_auxil-1)
!=====

 if( parallel_buffer ) then

#if 1
   
   call set_auxil_block_size(nbf_auxil_basis/(nprow_auxil*4))

   do iproc=0,nprow_auxil-1
     nbf_local_iproc(iproc) = NUMROC(nbf_auxil_basis,MBLOCK_AUXIL,iproc,first_row,nprow_auxil)
   enddo

   nauxil_3center = nbf_local_iproc(iprow_auxil)

   allocate(ibf_auxil_g(nauxil_3center))
   do ilocal=1,nauxil_3center
     ibf_auxil_g(ilocal) = INDXL2G(ilocal,MBLOCK_AUXIL,iprow_auxil,first_row,nprow_auxil)
   enddo
   allocate(ibf_auxil_l(nbf_auxil_basis))
   allocate(iproc_ibf_auxil(nbf_auxil_basis))
   do iglobal=1,nbf_auxil_basis
     ibf_auxil_l(iglobal)     = INDXG2L(iglobal,MBLOCK_AUXIL,0,first_row,nprow_auxil)
     iproc_ibf_auxil(iglobal) = INDXG2P(iglobal,MBLOCK_AUXIL,0,first_row,nprow_auxil)
   enddo

#else

  
   allocate(iproc_ibf_auxil(nbf_auxil_basis))
  
   iproc              = nproc_auxil - 1
   nbf_local_iproc(:) = 0
   do ibf=1,nbf_auxil_basis
  
     iproc = MODULO(iproc+1,nproc_auxil)
  
     iproc_ibf_auxil(ibf) = iproc
  
     nbf_local_iproc(iproc) = nbf_local_iproc(iproc) + 1
  
   enddo
  
   nauxil_3center = nbf_local_iproc(rank_auxil)
  
   allocate(ibf_auxil_g(nauxil_3center))
   allocate(ibf_auxil_l(nbf_auxil_basis))
   ibf_auxil_l(:) = 0
   ibf_local = 0
   do ibf=1,nbf_auxil_basis
     if( rank_auxil == iproc_ibf_auxil(ibf) ) then
       ibf_local = ibf_local + 1
       ibf_auxil_g(ibf_local) = ibf
       ibf_auxil_l(ibf)       = ibf_local
     endif
   enddo
#endif
  
 else

   ! Use SCALAPACK routines to distribute the auxiliary basis
   ! Assume a processor grid: nproc_auxil x 1

#if 1
   
   call set_auxil_block_size(nbf_auxil_basis/nprow_auxil/2)

   nauxil_3center = NUMROC(nbf_auxil_basis,MBLOCK_AUXIL,iprow_auxil,first_row,nprow_auxil)
   allocate(ibf_auxil_g(nauxil_3center))
   do ilocal=1,nauxil_3center
     ibf_auxil_g(ilocal) = INDXL2G(ilocal,MBLOCK_AUXIL,iprow_auxil,first_row,nprow_auxil)
   enddo
   allocate(ibf_auxil_l(nbf_auxil_basis))
   allocate(iproc_ibf_auxil(nbf_auxil_basis))
   do iglobal=1,nbf_auxil_basis
     ibf_auxil_l(iglobal)     = INDXG2L(iglobal,MBLOCK_AUXIL,0,first_row,nprow_auxil)
     iproc_ibf_auxil(iglobal) = INDXG2P(iglobal,MBLOCK_AUXIL,0,first_row,nprow_auxil)
   enddo

#else

   allocate(iproc_ibf_auxil(nbf_auxil_basis))
  
   iproc              = nproc_local-1
   nbf_local_iproc(:) = 0
   do ibf=1,nbf_auxil_basis
  
     iproc = MODULO(iproc+1,nproc_local)
  
     iproc_ibf_auxil(ibf) = iproc
  
     nbf_local_iproc(iproc) = nbf_local_iproc(iproc) + 1
  
   enddo
  
   nauxil_3center = nbf_local_iproc(rank_local)
  
   allocate(ibf_auxil_g(nauxil_3center))
   allocate(ibf_auxil_l(nbf_auxil_basis))
   ibf_auxil_l(:) = 0
   ibf_local = 0
   do ibf=1,nbf_auxil_basis
     if( rank_local == iproc_ibf_auxil(ibf) ) then
       ibf_local = ibf_local + 1
       ibf_auxil_g(ibf_local) = ibf
       ibf_auxil_l(ibf)       = ibf_local
     endif
   enddo

#endif
  
 endif

 write(stdout,'(/,a)') ' Distribute auxiliary basis functions among processors'
 write(stdout,'(1x,a,i4,a,i6,a)') 'Max auxiliary basis functions ',MAXVAL(nbf_local_iproc(:)),' for processor ',MAXLOC(nbf_local_iproc,DIM=1)
 write(stdout,'(1x,a,i4,a,i6,a)') 'Min auxiliary basis functions ',MINVAL(nbf_local_iproc(:)),' for processor ',MINLOC(nbf_local_iproc,DIM=1)


end subroutine distribute_auxil_basis


!=========================================================================
subroutine distribute_auxil_basis_lr(nbf_auxil_basis)
 use m_scalapack
 implicit none

 integer,intent(in)  :: nbf_auxil_basis
!=====
 integer :: ibf
 integer :: ibf_local
 integer :: iproc
 integer :: ilocal,iglobal
 integer :: nbf_local_iproc_lr(0:nproc_auxil-1)
!=====

   
#ifdef HAVE_SCALAPACK

 do iproc=0,nprow_auxil-1
   nbf_local_iproc_lr(iproc) = NUMROC(nbf_auxil_basis,MBLOCK_AUXIL,iproc,first_row,nprow_auxil)
 enddo

 nauxil_3center_lr = nbf_local_iproc_lr(iprow_auxil)

 allocate(ibf_auxil_g_lr(nauxil_3center_lr))
 do ilocal=1,nauxil_3center_lr
   ibf_auxil_g_lr(ilocal) = INDXL2G(ilocal,MBLOCK_AUXIL,iprow_auxil,first_row,nprow_auxil)
 enddo
 allocate(ibf_auxil_l_lr(nbf_auxil_basis))
 allocate(iproc_ibf_auxil_lr(nbf_auxil_basis))
 do iglobal=1,nbf_auxil_basis
   ibf_auxil_l_lr(iglobal)     = INDXG2L(iglobal,MBLOCK_AUXIL,0,first_row,nprow_auxil)
   iproc_ibf_auxil_lr(iglobal) = INDXG2P(iglobal,MBLOCK_AUXIL,0,first_row,nprow_auxil)
 enddo

#else

 allocate(iproc_ibf_auxil_lr(nbf_auxil_basis))

 iproc = nproc_auxil - 1
 nbf_local_iproc_lr(:) = 0
 do ibf=1,nbf_auxil_basis

   iproc = MODULO(iproc+1,nproc_auxil)

   iproc_ibf_auxil_lr(ibf) = iproc

   nbf_local_iproc_lr(iproc) = nbf_local_iproc_lr(iproc) + 1

 enddo

 nauxil_3center_lr = nbf_local_iproc_lr(rank_auxil)

 allocate(ibf_auxil_g_lr(nauxil_3center_lr))
 allocate(ibf_auxil_l_lr(nbf_auxil_basis))
 ibf_auxil_l_lr(:) = 0
 ibf_local = 0
 do ibf=1,nbf_auxil_basis
   if( rank_auxil == iproc_ibf_auxil_lr(ibf) ) then
     ibf_local = ibf_local + 1
     ibf_auxil_g_lr(ibf_local) = ibf
     ibf_auxil_l_lr(ibf)       = ibf_local
   endif
 enddo

#endif

 write(stdout,'(/,a)') ' Distribute LR auxiliary basis functions among processors'
 write(stdout,'(1x,a,i4,a,i6,a)')   'Max auxiliary basis functions ',MAXVAL(nbf_local_iproc_lr(:)),' for processor ',MAXLOC(nbf_local_iproc_lr,DIM=1)
 write(stdout,'(1x,a,i4,a,i6,a)')   'Min auxiliary basis functions ',MINVAL(nbf_local_iproc_lr(:)),' for processor ',MINLOC(nbf_local_iproc_lr,DIM=1)


end subroutine distribute_auxil_basis_lr


!=========================================================================
end module m_eri
