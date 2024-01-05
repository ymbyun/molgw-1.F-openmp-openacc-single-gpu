!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This modules contains
! the memory management
!=========================================================================
module m_memory
 use m_definitions
 use m_warning,only: die

 real(dp),private :: total_memory=0.0_dp     ! Total memory occupied 
                                             ! by the big arrays in Mb
 real(dp),private :: peak_memory=0.0_dp      ! Max memory occupied 
                                             ! by the big arrays in Mb

 interface clean_allocate
  module procedure clean_allocate_i1d
#ifdef ENABLE_YMBYUN
  module procedure clean_allocate_i1d_i8_4   ! ymbyun 2018/05/04
#endif
  module procedure clean_allocate_i2d
#ifdef ENABLE_YMBYUN
  module procedure clean_allocate_i2d_i4_48  ! ymbyun 2018/05/03
  module procedure clean_allocate_i2d_i4_84  ! ymbyun 2018/05/04
#endif
  module procedure clean_allocate_s1d
  module procedure clean_allocate_s2d
  module procedure clean_allocate_s3d
  module procedure clean_allocate_s4d
  module procedure clean_allocate_1d
#ifdef ENABLE_YMBYUN
  module procedure clean_allocate_1d_8       ! ymbyun 2018/05/03
#endif
  module procedure clean_allocate_2d
  module procedure clean_allocate_2d_range
  module procedure clean_allocate_3d
  module procedure clean_allocate_3d_range
  module procedure clean_allocate_4d
  module procedure clean_allocate_4d_range
 end interface

 interface clean_deallocate
  module procedure clean_deallocate_i1d
#ifdef ENABLE_YMBYUN
  module procedure clean_deallocate_i1d_i8   ! ymbyun 2018/05/05
#endif
  module procedure clean_deallocate_i2d
  module procedure clean_deallocate_s1d
  module procedure clean_deallocate_s2d
  module procedure clean_deallocate_s3d
  module procedure clean_deallocate_s4d
  module procedure clean_deallocate_1d
  module procedure clean_deallocate_2d
  module procedure clean_deallocate_3d
  module procedure clean_deallocate_4d
 end interface


contains


!=========================================================================
subroutine total_memory_statement()
 implicit none
!=====

 write(stdout,'(/,a,/)') '                 ----------------------------'
 write(stdout,'(a)')     '                 --- Memory usage summary ---'

 write(stdout,'(/,a)') ' Memory that was not deallocated properly'
 if( total_memory < 500._dp ) then
   write(stdout,'(a30,f9.3)') ' Memory (Mb): ',total_memory
 else
   write(stdout,'(a30,f9.3)') ' Memory (Gb): ',total_memory / 1024._dp
 endif

 write(stdout,'(/,a)') ' Maximum memory used during the run'
 if( peak_memory < 500._dp ) then
   write(stdout,'(a30,f9.3)') ' Peak memory (Mb): ',peak_memory
 else
   write(stdout,'(a30,f9.3)') ' Peak memory (Gb): ',peak_memory / 1024._dp
 endif

 write(stdout,*)

end subroutine total_memory_statement


!=========================================================================
subroutine clean_allocate_i1d(array_name,array,n1)
 implicit none

 character(len=*),intent(in)       :: array_name
 integer,allocatable,intent(inout) :: array(:)
 integer,intent(in)                :: n1
!=====
 integer             :: info
 real(dp)            :: mem_mb
!=====

 if( ALLOCATED(array) ) then
   call die('clean_allocate: Cannot allocate. This array is already allocated -> '//TRIM(array_name))
 endif

 mem_mb = REAL(4,dp) * REAL(n1,dp) / 1024._dp**2

 ! The allocation itself
 allocate(array(n1),stat=info)

 if(info/=0) then
   write(stdout,'(a,a)')    ' Failure when allocating ',array_name
   write(stdout,'(a,f9.3)') ' with size (Mb) ',mem_mb
   call die('clean_allocate: Not enough memory. Buy a bigger computer')
 endif


 total_memory = total_memory + mem_mb
 peak_memory = MAX(peak_memory,total_memory)

 call write_memory_allocate(array_name,mem_mb)

end subroutine clean_allocate_i1d

#ifdef ENABLE_YMBYUN
! ymbyun 2018/05/04
!=========================================================================
subroutine clean_allocate_i1d_i8_4(array_name,array,n1)
 implicit none

 character(len=*),intent(in)               :: array_name
! integer,allocatable,intent(inout)        :: array(:)
 integer(kind=8),allocatable,intent(inout) :: array(:)
 integer,intent(in)                        :: n1
!=====
 integer             :: info
 real(dp)            :: mem_mb
!=====

 if( ALLOCATED(array) ) then
   call die('clean_allocate: Cannot allocate. This array is already allocated -> '//TRIM(array_name))
 endif

! ymbyun 2018/05/04
! integer(kind=8) is used instead of integer(kind=4)
! mem_mb = REAL(4,dp) * REAL(n1,dp) / 1024._dp**2
 mem_mb = REAL(8,dp) * REAL(n1,dp) / 1024._dp**2

 ! The allocation itself
 allocate(array(n1),stat=info)

 if(info/=0) then
   write(stdout,'(a,a)')    ' Failure when allocating ',array_name
   write(stdout,'(a,f9.3)') ' with size (Mb) ',mem_mb
   call die('clean_allocate: Not enough memory. Buy a bigger computer')
 endif


 total_memory = total_memory + mem_mb
 peak_memory = MAX(peak_memory,total_memory)

 call write_memory_allocate(array_name,mem_mb)

end subroutine clean_allocate_i1d_i8_4
#endif

!=========================================================================
subroutine clean_allocate_i2d(array_name,array,n1,n2)
 implicit none

 character(len=*),intent(in)       :: array_name
 integer,allocatable,intent(inout) :: array(:,:)
 integer,intent(in)                :: n1,n2
!=====
 integer             :: info
 real(dp)            :: mem_mb
!=====

 if( ALLOCATED(array) ) then
   call die('clean_allocate: Cannot allocate. This array is already allocated -> '//TRIM(array_name))
 endif

 mem_mb = REAL(4,dp) * REAL(n1,dp) * REAL(n2,dp) / 1024._dp**2

 ! The allocation itself
 allocate(array(n1,n2),stat=info)

 if(info/=0) then
   write(stdout,'(a,a)')    ' Failure when allocating ',array_name
   write(stdout,'(a,f9.3)') ' with size (Mb) ',mem_mb
   call die('clean_allocate: Not enough memory. Buy a bigger computer')
 endif


 total_memory = total_memory + mem_mb
 peak_memory = MAX(peak_memory,total_memory)

 call write_memory_allocate(array_name,mem_mb)

end subroutine clean_allocate_i2d

#ifdef ENABLE_YMBYUN
! ymbyun 2018/05/03
!=========================================================================
subroutine clean_allocate_i2d_i4_48(array_name,array,n1,n2)
 implicit none

 character(len=*),intent(in)       :: array_name
 integer,allocatable,intent(inout) :: array(:,:)
! integer,intent(in)               :: n1,n2
 integer,intent(in)                :: n1
 integer(kind=8),intent(in)        :: n2
!=====
 integer             :: info
 real(dp)            :: mem_mb
!=====

 if( ALLOCATED(array) ) then
   call die('clean_allocate: Cannot allocate. This array is already allocated -> '//TRIM(array_name))
 endif

 mem_mb = REAL(4,dp) * REAL(n1,dp) * REAL(n2,dp) / 1024._dp**2

 ! The allocation itself
 allocate(array(n1,n2),stat=info)

 if(info/=0) then
   write(stdout,'(a,a)')    ' Failure when allocating ',array_name
   write(stdout,'(a,f9.3)') ' with size (Mb) ',mem_mb
   call die('clean_allocate: Not enough memory. Buy a bigger computer')
 endif


 total_memory = total_memory + mem_mb
 peak_memory = MAX(peak_memory,total_memory)

 call write_memory_allocate(array_name,mem_mb)

end subroutine clean_allocate_i2d_i4_48

! ymbyun 2018/05/04
!=========================================================================
subroutine clean_allocate_i2d_i4_84(array_name,array,n1,n2)
 implicit none

 character(len=*),intent(in)       :: array_name
 integer,allocatable,intent(inout) :: array(:,:)
! integer,intent(in)               :: n1,n2
 integer(kind=8),intent(in)        :: n1
 integer,intent(in)                :: n2
!=====
 integer             :: info
 real(dp)            :: mem_mb
!=====

 if( ALLOCATED(array) ) then
   call die('clean_allocate: Cannot allocate. This array is already allocated -> '//TRIM(array_name))
 endif

 mem_mb = REAL(4,dp) * REAL(n1,dp) * REAL(n2,dp) / 1024._dp**2

 ! The allocation itself
 allocate(array(n1,n2),stat=info)

 if(info/=0) then
   write(stdout,'(a,a)')    ' Failure when allocating ',array_name
   write(stdout,'(a,f9.3)') ' with size (Mb) ',mem_mb
   call die('clean_allocate: Not enough memory. Buy a bigger computer')
 endif


 total_memory = total_memory + mem_mb
 peak_memory = MAX(peak_memory,total_memory)

 call write_memory_allocate(array_name,mem_mb)

end subroutine clean_allocate_i2d_i4_84
#endif

!=========================================================================
subroutine clean_allocate_1d(array_name,array,n1)
 implicit none

 character(len=*),intent(in)       :: array_name
 real(dp),allocatable,intent(inout) :: array(:)
 integer,intent(in)                 :: n1
!=====
 integer             :: info
 real(dp)            :: mem_mb
!=====

 if( ALLOCATED(array) ) then
   call die('clean_allocate: Cannot allocate. This array is already allocated -> '//TRIM(array_name))
 endif

 mem_mb = REAL(dp,dp) * REAL(n1,dp) / 1024._dp**2

 ! The allocation itself
 allocate(array(n1),stat=info)

 if(info/=0) then
   write(stdout,*) 'failure'
   call die('clean_allocate: Not enough memory. Buy a bigger computer')
 endif


 total_memory = total_memory + mem_mb
 peak_memory = MAX(peak_memory,total_memory)

 call write_memory_allocate(array_name,mem_mb) 

end subroutine clean_allocate_1d

#ifdef ENABLE_YMBYUN
! ymbyun 2018/05/03
!=========================================================================
subroutine clean_allocate_1d_8(array_name,array,n1)
 implicit none

 character(len=*),intent(in)        :: array_name
 real(dp),allocatable,intent(inout) :: array(:)
! integer,intent(in)                :: n1
 integer(kind=8),intent(in)         :: n1
!=====
 integer             :: info
 real(dp)            :: mem_mb
!=====

 if( ALLOCATED(array) ) then
   call die('clean_allocate: Cannot allocate. This array is already allocated -> '//TRIM(array_name))
 endif

 mem_mb = REAL(dp,dp) * REAL(n1,dp) / 1024._dp**2

 ! The allocation itself
 allocate(array(n1),stat=info)

 if(info/=0) then
   write(stdout,*) 'failure'
   call die('clean_allocate: Not enough memory. Buy a bigger computer')
 endif


 total_memory = total_memory + mem_mb
 peak_memory = MAX(peak_memory,total_memory)

 call write_memory_allocate(array_name,mem_mb) 

end subroutine clean_allocate_1d_8
#endif

!=========================================================================
subroutine clean_allocate_2d(array_name,array,n1,n2)
 implicit none

 character(len=*),intent(in)       :: array_name
 real(dp),allocatable,intent(inout) :: array(:,:)
 integer,intent(in)                 :: n1,n2
!=====
 integer             :: info
 real(dp)            :: mem_mb
!=====

 if( ALLOCATED(array) ) then
   call die('clean_allocate: Cannot allocate. This array is already allocated -> '//TRIM(array_name))
 endif

 mem_mb = REAL(dp,dp) * REAL(n1,dp) * REAL(n2,dp) / 1024._dp**2

 ! The allocation itself
 allocate(array(n1,n2),stat=info)

 if(info/=0) then
   write(stdout,*) 'failure'
   call die('clean_allocate: Not enough memory. Buy a bigger computer')
 endif


 total_memory = total_memory + mem_mb
 peak_memory = MAX(peak_memory,total_memory)

 call write_memory_allocate(array_name,mem_mb) 

end subroutine clean_allocate_2d


!=========================================================================
subroutine clean_allocate_2d_range(array_name,array,n1s,n1f,n2s,n2f)
 implicit none

 character(len=*),intent(in)        :: array_name
 real(dp),allocatable,intent(inout) :: array(:,:)
 integer,intent(in)                 :: n1s,n1f,n2s,n2f
!=====
 integer             :: info
 real(dp)            :: mem_mb
!=====

 if( ALLOCATED(array) ) then
   call die('clean_allocate: Cannot allocate. This array is already allocated -> '//TRIM(array_name))
 endif

 mem_mb = REAL(dp,dp) * REAL(n1f-n1s+1,dp) * REAL(n2f-n2s+1,dp) / 1024._dp**2

 ! The allocation itself
 allocate(array(n1s:n1f,n2s:n2f),stat=info)

 if(info/=0) then
   write(stdout,*) 'failure'
   call die('clean_allocate: Not enough memory. Buy a bigger computer')
 endif


 total_memory = total_memory + mem_mb
 peak_memory = MAX(peak_memory,total_memory)

 call write_memory_allocate(array_name,mem_mb) 

end subroutine clean_allocate_2d_range


!=========================================================================
subroutine clean_allocate_3d(array_name,array,n1,n2,n3)
 implicit none

 character(len=*),intent(in)       :: array_name
 real(dp),allocatable,intent(inout) :: array(:,:,:)
 integer,intent(in)                 :: n1,n2,n3
!=====
 integer             :: info
 real(dp)            :: mem_mb
!=====

 if( ALLOCATED(array) ) then
   call die('clean_allocate: Cannot allocate. This array is already allocated -> '//TRIM(array_name))
 endif

 mem_mb = REAL(dp,dp) * REAL(n1,dp) * REAL(n2,dp) * REAL(n3,dp) / 1024._dp**2

 ! The allocation itself
 allocate(array(n1,n2,n3),stat=info)

 if(info/=0) then
   write(stdout,*) 'failure'
   call die('clean_allocate: Not enough memory. Buy a bigger computer')
 endif


 total_memory = total_memory + mem_mb
 peak_memory = MAX(peak_memory,total_memory)

 call write_memory_allocate(array_name,mem_mb) 

end subroutine clean_allocate_3d


!=========================================================================
subroutine clean_allocate_3d_range(array_name,array,n1s,n1f,n2s,n2f,n3s,n3f)
 implicit none

 character(len=*),intent(in)        :: array_name
 real(dp),allocatable,intent(inout) :: array(:,:,:)
 integer,intent(in)                 :: n1s,n1f,n2s,n2f,n3s,n3f
!=====
 integer             :: info
 real(dp)            :: mem_mb
!=====

 if( ALLOCATED(array) ) then
   call die('clean_allocate: Cannot allocate. This array is already allocated -> '//TRIM(array_name))
 endif

 mem_mb = REAL(dp,dp) * REAL(n1f-n1s+1,dp) * REAL(n2f-n2s+1,dp) * REAL(n3f-n3s+1,dp) / 1024._dp**2

 ! The allocation itself
 allocate(array(n1s:n1f,n2s:n2f,n3s:n3f),stat=info)

 if(info/=0) then
   write(stdout,*) 'failure'
   call die('clean_allocate: Not enough memory. Buy a bigger computer')
 endif


 total_memory = total_memory + mem_mb
 peak_memory = MAX(peak_memory,total_memory)

 call write_memory_allocate(array_name,mem_mb) 

end subroutine clean_allocate_3d_range


!=========================================================================
subroutine clean_allocate_4d(array_name,array,n1,n2,n3,n4)
 implicit none

 character(len=*),intent(in)       :: array_name
 real(dp),allocatable,intent(inout) :: array(:,:,:,:)
 integer,intent(in)                 :: n1,n2,n3,n4
!=====
 integer             :: info
 real(dp)            :: mem_mb
!=====

 if( ALLOCATED(array) ) then
   call die('clean_allocate: Cannot allocate. This array is already allocated -> '//TRIM(array_name))
 endif

 mem_mb = REAL(dp,dp) * REAL(n1,dp) * REAL(n2,dp) * REAL(n3,dp) * REAL(n4,dp) / 1024._dp**2

 ! The allocation itself
 allocate(array(n1,n2,n3,n4),stat=info)

 if(info/=0) then
   write(stdout,*) 'failure'
   call die('clean_allocate: Not enough memory. Buy a bigger computer')
 endif


 total_memory = total_memory + mem_mb
 peak_memory = MAX(peak_memory,total_memory)

 call write_memory_allocate(array_name,mem_mb) 

end subroutine clean_allocate_4d


!=========================================================================
subroutine clean_allocate_4d_range(array_name,array,n1s,n1f,n2s,n2f,n3s,n3f,n4s,n4f)
 implicit none

 character(len=*),intent(in)        :: array_name
 real(dp),allocatable,intent(inout) :: array(:,:,:,:)
 integer,intent(in)                 :: n1s,n1f,n2s,n2f,n3s,n3f,n4s,n4f
!=====
 integer             :: info
 real(dp)            :: mem_mb
!=====

 if( ALLOCATED(array) ) then
   call die('clean_allocate: Cannot allocate. This array is already allocated -> '//TRIM(array_name))
 endif

 mem_mb = REAL(dp,dp) * REAL(n1f-n1s+1,dp) * REAL(n2f-n2s+1,dp) * REAL(n3f-n3s+1,dp) * REAL(n4f-n4s+1,dp) / 1024._dp**2

 ! The allocation itself
 allocate(array(n1s:n1f,n2s:n2f,n3s:n3f,n4s:n4f),stat=info)

 if(info/=0) then
   write(stdout,*) 'failure'
   call die('clean_allocate: Not enough memory. Buy a bigger computer')
 endif


 total_memory = total_memory + mem_mb
 peak_memory = MAX(peak_memory,total_memory)

 call write_memory_allocate(array_name,mem_mb) 

end subroutine clean_allocate_4d_range


!=========================================================================
subroutine clean_deallocate_i1d(array_name,array)
 implicit none

 character(len=*),intent(in)               :: array_name
 integer,allocatable,intent(inout)        :: array(:)
!=====
 integer             :: info
 real(dp)            :: mem_mb
 integer             :: n1
!=====

 if( .NOT. ALLOCATED(array) ) return

 n1 = SIZE(array(:),DIM=1)

 mem_mb = REAL(4,dp) * REAL(n1,dp) / 1024._dp**2

 ! The allocation itself
 deallocate(array)

 total_memory = total_memory - mem_mb

 call write_memory_deallocate(array_name,mem_mb)

end subroutine clean_deallocate_i1d

#ifdef ENABLE_YMBYUN
! ymbyun 2018/05/05
!=========================================================================
subroutine clean_deallocate_i1d_i8(array_name,array)
 implicit none

 character(len=*),intent(in)               :: array_name
! integer,allocatable,intent(inout)        :: array(:)
 integer(kind=8),allocatable,intent(inout) :: array(:)
!=====
 integer             :: info
 real(dp)            :: mem_mb
 integer             :: n1
!=====

 if( .NOT. ALLOCATED(array) ) return

 n1 = SIZE(array(:),DIM=1)

! ymbyun 2018/05/05
! integer(kind=8) is used instead of integer(kind=4).
! mem_mb = REAL(4,dp) * REAL(n1,dp) / 1024._dp**2
 mem_mb = REAL(8,dp) * REAL(n1,dp) / 1024._dp**2

 ! The allocation itself
 deallocate(array)

 total_memory = total_memory - mem_mb

 call write_memory_deallocate(array_name,mem_mb)

end subroutine clean_deallocate_i1d_i8
#endif

!=========================================================================
subroutine clean_deallocate_i2d(array_name,array)
 implicit none

 character(len=*),intent(in)       :: array_name
 integer,allocatable,intent(inout) :: array(:,:)
!=====
 integer             :: info
 real(dp)            :: mem_mb
! ymbyun 2018/05/05
#ifdef ENABLE_YMBYUN
 integer(kind=8)     :: n1,n2
#else
 integer             :: n1,n2
#endif
!=====

 if( .NOT. ALLOCATED(array) ) return

! ymbyun 2018/05/05
! I'm not sure if SIZE() can return an 8-byte integer.
 n1 = SIZE(array(:,:),DIM=1)
 n2 = SIZE(array(:,:),DIM=2)

 mem_mb = REAL(4,dp) * REAL(n1,dp) * REAL(n2,dp) / 1024._dp**2

 ! The allocation itself
 deallocate(array)

 total_memory = total_memory - mem_mb

 call write_memory_deallocate(array_name,mem_mb)

end subroutine clean_deallocate_i2d


!=========================================================================
subroutine clean_deallocate_1d(array_name,array)
 implicit none

 character(len=*),intent(in)        :: array_name
 real(dp),allocatable,intent(inout) :: array(:)
!=====
 integer             :: info
 real(dp)            :: mem_mb
! ymbyun 2018/05/03
#ifdef ENABLE_YMBYUN
 integer(kind=8)     :: n1
#else
 integer             :: n1
#endif
!=====

 if( .NOT. ALLOCATED(array) ) return

! ymbyun 2018/05/05
! I'm not sure if SIZE() can return an 8-byte integer.
 n1 = SIZE(array(:))

 mem_mb = REAL(dp,dp) * REAL(n1,dp) / 1024._dp**2

 ! The allocation itself
 deallocate(array)

 total_memory = total_memory - mem_mb

 call write_memory_deallocate(array_name,mem_mb) 

end subroutine clean_deallocate_1d


!=========================================================================
subroutine clean_deallocate_2d(array_name,array)
 implicit none

 character(len=*),intent(in)       :: array_name
 real(dp),allocatable,intent(inout) :: array(:,:)
!=====
 integer             :: info
 real(dp)            :: mem_mb
 integer             :: n1,n2
!=====

 if( .NOT. ALLOCATED(array) ) return

 n1 = SIZE(array(:,:),DIM=1)
 n2 = SIZE(array(:,:),DIM=2)

 mem_mb = REAL(dp,dp) * REAL(n1,dp) * REAL(n2,dp) / 1024._dp**2

 ! The allocation itself
 deallocate(array)

 total_memory = total_memory - mem_mb

 call write_memory_deallocate(array_name,mem_mb) 

end subroutine clean_deallocate_2d


!=========================================================================
subroutine clean_deallocate_3d(array_name,array)
 implicit none

 character(len=*),intent(in)       :: array_name
 real(dp),allocatable,intent(inout) :: array(:,:,:)
!=====
 integer             :: info
 real(dp)            :: mem_mb
 integer             :: n1,n2,n3
!=====

 if( .NOT. ALLOCATED(array) ) return

 n1 = SIZE(array(:,:,:),DIM=1)
 n2 = SIZE(array(:,:,:),DIM=2)
 n3 = SIZE(array(:,:,:),DIM=3)

 mem_mb = REAL(dp,dp) * REAL(n1,dp) *REAL(n2,dp) * REAL(n3,dp) / 1024._dp**2

 ! The allocation itself
 deallocate(array)

 total_memory = total_memory - mem_mb

 call write_memory_deallocate(array_name,mem_mb) 

end subroutine clean_deallocate_3d


!=========================================================================
subroutine clean_deallocate_4d(array_name,array)
 implicit none

 character(len=*),intent(in)       :: array_name
 real(dp),allocatable,intent(inout) :: array(:,:,:,:)
!=====
 integer             :: info
 real(dp)            :: mem_mb
 integer             :: n1,n2,n3,n4
!=====

 if( .NOT. ALLOCATED(array) ) return

 n1 = SIZE(array(:,:,:,:),DIM=1)
 n2 = SIZE(array(:,:,:,:),DIM=2)
 n3 = SIZE(array(:,:,:,:),DIM=3)
 n4 = SIZE(array(:,:,:,:),DIM=4)

 mem_mb = REAL(dp,dp) * REAL(n1,dp) *REAL(n2,dp) * REAL(n3,dp) * REAL(n4,dp)/ 1024._dp**2

 ! The allocation itself
 deallocate(array)

 total_memory = total_memory - mem_mb

 call write_memory_deallocate(array_name,mem_mb) 

end subroutine clean_deallocate_4d


!=========================================================================
subroutine clean_allocate_s1d(array_name,array,n1)
 implicit none

 character(len=*),intent(in)       :: array_name
 real(sp),allocatable,intent(inout) :: array(:)
 integer,intent(in)                 :: n1
!=====
 integer             :: info
 real(dp)            :: mem_mb
!=====

 if( ALLOCATED(array) ) then
   call die('clean_allocate: Cannot allocate. This array is already allocated -> '//TRIM(array_name))
 endif

 mem_mb = REAL(sp,dp) * REAL(n1,dp) / 1024._dp**2

 ! The allocation itself
 allocate(array(n1),stat=info)

 if(info/=0) then
   write(stdout,*) 'failure'
   call die('clean_allocate: Not enough memory. Buy a bigger computer')
 endif


 total_memory = total_memory + mem_mb
 peak_memory = MAX(peak_memory,total_memory)

 call write_memory_allocate(array_name,mem_mb) 

end subroutine clean_allocate_s1d


!=========================================================================
subroutine clean_allocate_s2d(array_name,array,n1,n2)
 implicit none

 character(len=*),intent(in)       :: array_name
 real(sp),allocatable,intent(inout) :: array(:,:)
 integer,intent(in)                 :: n1,n2
!=====
 integer             :: info
 real(dp)            :: mem_mb
!=====

 if( ALLOCATED(array) ) then
   call die('clean_allocate: Cannot allocate. This array is already allocated -> '//TRIM(array_name))
 endif

 mem_mb = REAL(sp,dp) * REAL(n1,dp) * REAL(n2,dp) / 1024._dp**2

 ! The allocation itself
 allocate(array(n1,n2),stat=info)

 if(info/=0) then
   write(stdout,*) 'failure'
   call die('clean_allocate: Not enough memory. Buy a bigger computer')
 endif


 total_memory = total_memory + mem_mb
 peak_memory = MAX(peak_memory,total_memory)

 call write_memory_allocate(array_name,mem_mb) 

end subroutine clean_allocate_s2d


!=========================================================================
subroutine clean_allocate_s3d(array_name,array,n1,n2,n3)
 implicit none

 character(len=*),intent(in)       :: array_name
 real(sp),allocatable,intent(inout) :: array(:,:,:)
 integer,intent(in)                 :: n1,n2,n3
!=====
 integer             :: info
 real(dp)            :: mem_mb
!=====

 if( ALLOCATED(array) ) then
   call die('clean_allocate: Cannot allocate. This array is already allocated -> '//TRIM(array_name))
 endif

 mem_mb = REAL(sp,dp) * REAL(n1,dp) * REAL(n2,dp) * REAL(n3,dp) / 1024._dp**2

 ! The allocation itself
 allocate(array(n1,n2,n3),stat=info)

 if(info/=0) then
   write(stdout,*) 'failure'
   call die('clean_allocate: Not enough memory. Buy a bigger computer')
 endif


 total_memory = total_memory + mem_mb
 peak_memory = MAX(peak_memory,total_memory)

 call write_memory_allocate(array_name,mem_mb) 

end subroutine clean_allocate_s3d


!=========================================================================
subroutine clean_allocate_s4d(array_name,array,n1,n2,n3,n4)
 implicit none

 character(len=*),intent(in)       :: array_name
 real(sp),allocatable,intent(inout) :: array(:,:,:,:)
 integer,intent(in)                 :: n1,n2,n3,n4
!=====
 integer             :: info
 real(dp)            :: mem_mb
!=====

 if( ALLOCATED(array) ) then
   call die('clean_allocate: Cannot allocate. This array is already allocated -> '//TRIM(array_name))
 endif

 mem_mb = REAL(sp,dp) * REAL(n1,dp) * REAL(n2,dp) * REAL(n3,dp) * REAL(n4,dp) / 1024._dp**2

 ! The allocation itself
 allocate(array(n1,n2,n3,n4),stat=info)

 if(info/=0) then
   write(stdout,*) 'failure'
   call die('clean_allocate: Not enough memory. Buy a bigger computer')
 endif


 total_memory = total_memory + mem_mb
 peak_memory = MAX(peak_memory,total_memory)

 call write_memory_allocate(array_name,mem_mb) 

end subroutine clean_allocate_s4d


!=========================================================================
subroutine clean_deallocate_s1d(array_name,array)
 implicit none

 character(len=*),intent(in)       :: array_name
 real(sp),allocatable,intent(inout) :: array(:)
!=====
 integer             :: info
 real(dp)            :: mem_mb
 integer             :: n1
!=====

 if( .NOT. ALLOCATED(array) ) return

 n1 = SIZE(array(:))

 mem_mb = REAL(sp,dp) * REAL(n1,dp) / 1024._dp**2

 ! The allocation itself
 deallocate(array)

 total_memory = total_memory - mem_mb

 call write_memory_deallocate(array_name,mem_mb) 

end subroutine clean_deallocate_s1d


!=========================================================================
subroutine clean_deallocate_s2d(array_name,array)
 implicit none

 character(len=*),intent(in)       :: array_name
 real(sp),allocatable,intent(inout) :: array(:,:)
!=====
 integer             :: info
 real(dp)            :: mem_mb
 integer             :: n1,n2
!=====

 if( .NOT. ALLOCATED(array) ) return

 n1 = SIZE(array(:,:),DIM=1)
 n2 = SIZE(array(:,:),DIM=2)

 mem_mb = REAL(sp,dp) * REAL(n1,dp) * REAL(n2,dp) / 1024._dp**2

 ! The allocation itself
 deallocate(array)

 total_memory = total_memory - mem_mb

 call write_memory_deallocate(array_name,mem_mb) 

end subroutine clean_deallocate_s2d


!=========================================================================
subroutine clean_deallocate_s3d(array_name,array)
 implicit none

 character(len=*),intent(in)       :: array_name
 real(sp),allocatable,intent(inout) :: array(:,:,:)
!=====
 integer             :: info
 real(dp)            :: mem_mb
 integer             :: n1,n2,n3
!=====

 if( .NOT. ALLOCATED(array) ) return

 n1 = SIZE(array(:,:,:),DIM=1)
 n2 = SIZE(array(:,:,:),DIM=2)
 n3 = SIZE(array(:,:,:),DIM=3)

 mem_mb = REAL(sp,dp) * REAL(n1,dp) *REAL(n2,dp) * REAL(n3,dp) / 1024._dp**2

 ! The allocation itself
 deallocate(array)

 total_memory = total_memory - mem_mb

 call write_memory_deallocate(array_name,mem_mb) 

end subroutine clean_deallocate_s3d


!=========================================================================
subroutine clean_deallocate_s4d(array_name,array)
 implicit none

 character(len=*),intent(in)       :: array_name
 real(sp),allocatable,intent(inout) :: array(:,:,:,:)
!=====
 integer             :: info
 real(dp)            :: mem_mb
 integer             :: n1,n2,n3,n4
!=====

 if( .NOT. ALLOCATED(array) ) return

 n1 = SIZE(array(:,:,:,:),DIM=1)
 n2 = SIZE(array(:,:,:,:),DIM=2)
 n3 = SIZE(array(:,:,:,:),DIM=3)
 n4 = SIZE(array(:,:,:,:),DIM=4)

 mem_mb = REAL(sp,dp) * REAL(n1,dp) *REAL(n2,dp) * REAL(n3,dp) * REAL(n4,dp)/ 1024._dp**2

 ! The allocation itself
 deallocate(array)

 total_memory = total_memory - mem_mb

 call write_memory_deallocate(array_name,mem_mb) 

end subroutine clean_deallocate_s4d


!=========================================================================
subroutine write_memory_allocate(array_name,mem_mb)
 implicit none
!=====
 character(len=*),intent(in)       :: array_name
 real(dp),intent(in)               :: mem_mb
!=====

 if( mem_mb < 500.0_dp .AND. total_memory < 500.0_dp ) then

   write(stdout,'(1x,a,a35,a,sp,f8.3,s,a,a,f8.3,a)')   &
         '  Allocate: ',TRIM(array_name),        &
             ',  Mem: ', mem_mb      ,' (Mb)',   &
       ',  Total Mem: ', total_memory,' (Mb)'

 else

   write(stdout,'(1x,a,a35,a,sp,f8.3,s,a,a,f8.3,a)')   &
         '  Allocate: ',TRIM(array_name),        &
             ',  Mem: ', mem_mb/1024.0_dp      ,' (Gb)',   &
       ',  Total Mem: ', total_memory/1024.0_dp,' (Gb)'

 endif


end subroutine write_memory_allocate


!=========================================================================
subroutine write_memory_deallocate(array_name,mem_mb)
 implicit none
!=====
 character(len=*),intent(in)       :: array_name
 real(dp),intent(in)               :: mem_mb
!=====

 if( mem_mb < 500.0_dp .AND. total_memory < 500.0_dp ) then

   write(stdout,'(1x,a,a35,a,sp,f8.3,s,a,a,f8.3,a)')   &
         'Deallocate: ',TRIM(array_name),        &
             ',  Mem: ',-mem_mb      ,' (Mb)',   &
       ',  Total Mem: ', total_memory,' (Mb)'

 else

   write(stdout,'(1x,a,a35,a,sp,f8.3,s,a,a,f8.3,a)')   &
         'Deallocate: ',TRIM(array_name),        &
             ',  Mem: ',-mem_mb/1024.0_dp      ,' (Gb)',   &
       ',  Total Mem: ', total_memory/1024.0_dp,' (Gb)'

 endif


end subroutine write_memory_deallocate


!! ymbyun 2021/12/27
!! NOTE: This is not working. I don't know why yet.
!!=========================================================================
!subroutine write_openacc_memory()
! use openacc
! implicit none
!!=====
! integer                  :: device_num
! integer(acc_device_kind) :: device_type
! integer(kind=8)          :: gpu_memory   , gpu_free_memory
! integer                  :: gpu_memory_mb, gpu_free_memory_mb
!!=====
!
! gpu_memory         = acc_get_property(device_num,device_type,acc_property_memory     )
! gpu_free_memory    = acc_get_property(device_num,device_type,acc_property_free_memory)
!
! gpu_memory_mb      = gpu_memory      / 1024 / 1024
! gpu_free_memory_mb = gpu_free_memory / 1024 / 1024
!
! write(stdout,*) 'OpenACC total/free/used memory (MB):',gpu_memory_mb,gpu_free_memory_mb,(gpu_memory_mb-gpu_free_memory_mb)
!end subroutine write_openacc_memory

end module m_memory
!=========================================================================
