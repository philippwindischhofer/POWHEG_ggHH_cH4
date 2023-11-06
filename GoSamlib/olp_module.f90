module     olp_module
   use, intrinsic :: iso_c_binding;
   implicit none
   private
   public :: OLP_Start, OLP_EvalSubProcess, OLP_Finalize, OLP_Option
   public :: OLP_EvalSubProcess2, OLP_Polvec, OLP_SetParameter, OLP_Info
   public :: OLP_PrintParameter

contains

   subroutine     OLP_Start(contract_file_name,ierr,stage,rndseed) &
   & bind(C,name="olp_start_")
      use, intrinsic :: iso_c_binding
      use p0_part21part21_part25part25_matrix, only: p0_part21part21_part25part25_initgolem => initgolem
      use p0_part21part21_part25part25_config, only: p0_part21part21_part25part25_PSP_rescue => PSP_rescue, &
           & p0_part21part21_part25part25_PSP_verbosity => PSP_verbosity, &
           & p0_part21part21_part25part25_PSP_chk_th1 => PSP_chk_th1, &
           & p0_part21part21_part25part25_PSP_chk_th2 => PSP_chk_th2, &
           & p0_part21part21_part25part25_PSP_chk_th3 => PSP_chk_th3, &
           & p0_part21part21_part25part25_PSP_chk_kfactor => PSP_chk_kfactor
      use pb_part21part21_part25part25_matrix, only: pb_part21part21_part25part25_initgolem => initgolem
      use pb_part21part21_part25part25_config, only: pb_part21part21_part25part25_PSP_rescue => PSP_rescue, &
           & pb_part21part21_part25part25_PSP_verbosity => PSP_verbosity, &
           & pb_part21part21_part25part25_PSP_chk_th1 => PSP_chk_th1, &
           & pb_part21part21_part25part25_PSP_chk_th2 => PSP_chk_th2, &
           & pb_part21part21_part25part25_PSP_chk_th3 => PSP_chk_th3, &
           & pb_part21part21_part25part25_PSP_chk_kfactor => PSP_chk_kfactor
      use p2_part21part21_part25part25part21_matrix, only: p2_part21part21_part25part25part21_initgolem => initgolem
      use p2_part21part21_part25part25part21_config, only: p2_part21part21_part25part25part21_PSP_rescue => PSP_rescue, &
           & p2_part21part21_part25part25part21_PSP_verbosity => PSP_verbosity, &
           & p2_part21part21_part25part25part21_PSP_chk_th1 => PSP_chk_th1, &
           & p2_part21part21_part25part25part21_PSP_chk_th2 => PSP_chk_th2, &
           & p2_part21part21_part25part25part21_PSP_chk_th3 => PSP_chk_th3, &
           & p2_part21part21_part25part25part21_PSP_chk_kfactor => PSP_chk_kfactor
      use p9_part21part21_part25part25part21_matrix, only: p9_part21part21_part25part25part21_initgolem => initgolem
      use p9_part21part21_part25part25part21_config, only: p9_part21part21_part25part25part21_PSP_rescue => PSP_rescue, &
           & p9_part21part21_part25part25part21_PSP_verbosity => PSP_verbosity, &
           & p9_part21part21_part25part25part21_PSP_chk_th1 => PSP_chk_th1, &
           & p9_part21part21_part25part25part21_PSP_chk_th2 => PSP_chk_th2, &
           & p9_part21part21_part25part25part21_PSP_chk_th3 => PSP_chk_th3, &
           & p9_part21part21_part25part25part21_PSP_chk_kfactor => PSP_chk_kfactor
      use p3_part1part21_part25part25part1_matrix, only: p3_part1part21_part25part25part1_initgolem => initgolem
      use p3_part1part21_part25part25part1_config, only: p3_part1part21_part25part25part1_PSP_rescue => PSP_rescue, &
           & p3_part1part21_part25part25part1_PSP_verbosity => PSP_verbosity, &
           & p3_part1part21_part25part25part1_PSP_chk_th1 => PSP_chk_th1, &
           & p3_part1part21_part25part25part1_PSP_chk_th2 => PSP_chk_th2, &
           & p3_part1part21_part25part25part1_PSP_chk_th3 => PSP_chk_th3, &
           & p3_part1part21_part25part25part1_PSP_chk_kfactor => PSP_chk_kfactor
      use p10_part1part21_part25part25part1_matrix, only: p10_part1part21_part25part25part1_initgolem => initgolem
      use p10_part1part21_part25part25part1_config, only: p10_part1part21_part25part25part1_PSP_rescue => PSP_rescue, &
           & p10_part1part21_part25part25part1_PSP_verbosity => PSP_verbosity, &
           & p10_part1part21_part25part25part1_PSP_chk_th1 => PSP_chk_th1, &
           & p10_part1part21_part25part25part1_PSP_chk_th2 => PSP_chk_th2, &
           & p10_part1part21_part25part25part1_PSP_chk_th3 => PSP_chk_th3, &
           & p10_part1part21_part25part25part1_PSP_chk_kfactor => PSP_chk_kfactor
 !     use p4_part21part1_part25part25part1_matrix, only: p4_part21part1_part25part25part1_initgolem => initgolem
 !     use p11_part21part1_part25part25part1_matrix, only: p11_part21part1_part25part25part1_initgolem => initgolem
 !     use p5_partm1part21_part25part25partm1_matrix, only: p5_partm1part21_part25part25partm1_initgolem => initgolem
 !     use p12_partm1part21_part25part25partm1_matrix, only: p12_partm1part21_part25part25partm1_initgolem => initgolem
 !     use p6_part21partm1_part25part25partm1_matrix, only: p6_part21partm1_part25part25partm1_initgolem => initgolem
 !     use p13_part21partm1_part25part25partm1_matrix, only: p13_part21partm1_part25part25partm1_initgolem => initgolem
 !     use p7_part1partm1_part25part25part21_matrix, only: p7_part1partm1_part25part25part21_initgolem => initgolem
 !     use p14_part1partm1_part25part25part21_matrix, only: p14_part1partm1_part25part25part21_initgolem => initgolem
 !     use p8_partm1part1_part25part25part21_matrix, only: p8_partm1part1_part25part25part21_initgolem => initgolem
 !     use p15_partm1part1_part25part25part21_matrix, only: p15_partm1part1_part25part25part21_initgolem => initgolem



      implicit none
      character(kind=c_char,len=1), intent(in) :: contract_file_name
      integer(kind=c_int), intent(out) :: ierr
      integer(kind=c_int), intent(in) :: stage, rndseed
      interface
         function strlen(s) bind(C,name='strlen')
            use, intrinsic :: iso_c_binding
            implicit none
            character(kind=c_char,len=1), intent(in) :: s
            integer(kind=c_int) :: strlen
         end function strlen
      end interface

      integer :: l, ferr
      character(len=128) :: line_buf
      character(len=9) :: kw
      integer :: PSP_verbosity, PSP_chk_th1, PSP_chk_th2, PSP_chk_th3, PSP_chk_kfactor
      logical :: PSP_rescue

      ierr = 1
      l = strlen(contract_file_name)

      open(unit=21, file=contract_file_name(1:l), &
          & status='old', action='read', iostat=ferr)

      if (ferr .ne. 0) then
         write(7,*) "In OLP_Start: ", contract_file_name(1:l), " not found!"
         ierr = -1
      end if

      do while (ferr .eq. 0)
         read(unit=21,fmt='(A128)',iostat=ferr) line_buf
         if (ferr .ne. 0) exit
         line_buf = adjustl(line_buf)
         kw = line_buf(1:9)
         do
            l = scan(kw, "DEFILMO")
            if (l .eq. 0) exit
            kw(l:l) = achar(ichar(kw(l:l)) - ichar('A') + ichar('a'))
         end do
         if (kw .eq. "modelfile") then
            line_buf = adjustl(line_buf(10:128))
            l = scan(line_buf, "|") - 1
            if(l .lt. 1) l = len(line_buf)
            l = len_trim(line_buf(1:l))
            exit
         end if
      end do

      close(unit=21)

      if (ierr .eq. 1) then
         call read_slha_file(line_buf(1:l))
      end if

      ! Uncomment to change rescue system setting on all suprocesses
      ! PSP_rescue = .true.
      ! PSP_verbosity = .false.
      ! PSP_chk_th1 = 8
      ! PSP_chk_th2 = 3
      ! PSP_chk_th3 = 5
      ! PSP_chk_kfactor = 1000
      ! p0_part21part21_part25part25_PSP_rescue = PSP_rescue
      ! p0_part21part21_part25part25_PSP_verbosity =  PSP_verbosity
      ! p0_part21part21_part25part25_PSP_chk_th1 = PSP_chk_th1
      ! p0_part21part21_part25part25_PSP_chk_th2 = PSP_chk_th2
      ! p0_part21part21_part25part25_PSP_chk_th3 = PSP_chk_th3
      ! p0_part21part21_part25part25_PSP_chk_kfactor = PSP_chk_kfactor
      if(stage.lt.0) then
         call p0_part21part21_part25part25_initgolem(.true.)
         call pb_part21part21_part25part25_initgolem(.false.)
         call p2_part21part21_part25part25part21_initgolem(.false.)
         call p9_part21part21_part25part25part21_initgolem(.false.)
         call p3_part1part21_part25part25part1_initgolem(.false.)
         call p10_part1part21_part25part25part1_initgolem(.false.)
!         call p4_part25part1_part25part25part1_initgolem(.false.)
!         call p11_part25part1_part25part25part1_initgolem(.false.)
!         call p5_partm1part21_part25part25partm1_initgolem(.false.)
!         call p12_partm1part21_part25part25partm1_initgolem(.false.)
!         call p6_part21partm1_part25part25partm1_initgolem(.false.)
!         call p13_part21partm1_part25part25partm1_initgolem(.false.)
!         call p7_part1partm1_part25part25part21_initgolem(.false.)
!         call p14_part1partm1_part25part25part21_initgolem(.false.)
      else
         call p0_part21part21_part25part25_initgolem(.true.,stage,rndseed)
         call pb_part21part21_part25part25_initgolem(.false.,stage,rndseed)
         call p2_part21part21_part25part25part21_initgolem(.false.,stage,rndseed)
         call p9_part21part21_part25part25part21_initgolem(.false.,stage,rndseed)
         call p3_part1part21_part25part25part1_initgolem(.false.,stage,rndseed)
         call p10_part1part21_part25part25part1_initgolem(.false.,stage,rndseed)
!         call p5_partm1part21_partm1part25part25_initgolem(.false.,stage,rndseed)
!         call p6_part21partm1_partm1part25part25_initgolem(.false.,stage,rndseed)
!         call p7_part1partm1_part21part25part25_initgolem(.false.,stage,rndseed)
      end if

   end subroutine OLP_Start

   subroutine OLP_Info(olp_name, olp_version, message)&
   & bind(C,name="olp_info_")
   use, intrinsic :: iso_c_binding, only: C_CHAR, C_NULL_CHAR
   use p0_part21part21_part25part25_version, only: gosamversion, gosamrevision

   implicit none
   character(kind=c_char), intent(inout), dimension(20)  :: olp_name
   character(kind=c_char), intent(inout), dimension(20)  :: olp_version
   character(kind=c_char), intent(inout), dimension(255) :: message
   character(len=254) :: msg
   character(len=7)   :: rev
   character(len=1)   :: ver1, ver2

   interface
     subroutine strncpy(dest, src, n) bind(C)
       import
       character(kind=c_char),  intent(out) :: dest(*)
       character(kind=c_char),  intent(in)  :: src(*)
       integer(c_size_t), value, intent(in) :: n
     end subroutine strncpy
   end interface

   write(ver1,'(I1)') gosamversion(1)
   write(ver2,'(I1)') gosamversion(2)
   write(rev,'(Z7)')  gosamrevision

   msg = new_line('C')//'Please cite the following references when using this program:'//&
        &new_line('C')//'Peraro:2014cba, Mastrolia:2012bu, '//&
        &'Mastrolia:2010nb, Guillet:2013msa, Cullen:2011kv, Cullen:2014yla'

   call strncpy(olp_name, c_char_'GoSam-'//ver1//'.'//ver2//C_NULL_CHAR, &
        len(c_char_'GoSam-'//ver1//'.'//ver2//C_NULL_CHAR,kind=c_size_t))

   call strncpy(olp_version, c_char_'git rev-'//trim(adjustl(rev))//C_NULL_CHAR, &
        len(c_char_'git rev-'//trim(adjustl(rev))//C_NULL_CHAR,kind=c_size_t))

   call strncpy(message, c_char_'GoSam:'//trim(adjustl(msg))//C_NULL_CHAR, &
        len(c_char_'GoSam:'//trim(adjustl(msg))//C_NULL_CHAR,kind=c_size_t))

   end subroutine OLP_Info


   subroutine OLP_SetParameter(variable_name, real_part, imag_part, success)&
   & bind(C,name="olp_setparameter_")
      use, intrinsic :: iso_c_binding
      use p0_part21part21_part25part25_model, only: p0_part21part21_part25part25_set_parameter => set_parameter
      use pb_part21part21_part25part25_model, only: pb_part21part21_part25part25_set_parameter => set_parameter
      use p2_part21part21_part25part25part21_model, only: p2_part21part21_part25part25part21_set_parameter => set_parameter
      use p9_part21part21_part25part25part21_model, only: p9_part21part21_part25part25part21_set_parameter => set_parameter
      use p3_part1part21_part25part25part1_model, only: p3_part1part21_part25part25part1_set_parameter => set_parameter
      use p10_part1part21_part25part25part1_model, only: p10_part1part21_part25part25part1_set_parameter => set_parameter
!      use p5_partm1part21_partm1part25part25_model, only: p5_partm1part21_partm1part25part25_set_parameter => set_parameter
!      use p6_part21partm1_partm1part25part25_model, only: p6_part21partm1_partm1part25part25_set_parameter => set_parameter
!      use p7_part1partm1_part21part25part25_model, only: p7_part1partm1_part21part25part25_set_parameter => set_parameter
      implicit none
      character(kind=c_char,len=1), intent(in) :: variable_name
      real(kind=c_double), intent(in) :: real_part, imag_part
      integer(kind=c_int), intent(out) :: success

      interface
         function strlen(s) bind(C,name='strlen')
            use, intrinsic :: iso_c_binding
            implicit none
            character(kind=c_char,len=1), intent(in) :: s
            integer(kind=c_int) :: strlen
         end function strlen
      end interface

      integer :: l;

      l = strlen(variable_name)
      call p0_part21part21_part25part25_set_parameter(variable_name(1:l),real_part,imag_part,success)
      if(success==0) then ! return immediately on error
          return
      end if
      call pb_part21part21_part25part25_set_parameter(variable_name(1:l),real_part,imag_part,success)
      if(success==0) then ! return immediately on error
          return
      end if
      call p2_part21part21_part25part25part21_set_parameter(variable_name(1:l),real_part,imag_part,success)
      if(success==0) then ! return immediately on error
          return
      end if
      call p9_part21part21_part25part25part21_set_parameter(variable_name(1:l),real_part,imag_part,success)
      if(success==0) then ! return immediately on error
          return
      end if
      call p3_part1part21_part25part25part1_set_parameter(variable_name(1:l),real_part,imag_part,success)
      if(success==0) then ! return immediately on error
          return
      end if
      call p10_part1part21_part25part25part1_set_parameter(variable_name(1:l),real_part,imag_part,success)
      if(success==0) then ! return immediately on error
          return
      end if
!      call p5_partm1part21_partm1part25part25_set_parameter(variable_name(1:l),real_part,imag_part,success)
!      if(success==0) then ! return immediately on error
!          return
!      end if
!      call p6_part21partm1_partm1part25part25_set_parameter(variable_name(1:l),real_part,imag_part,success)
!      if(success==0) then ! return immediately on error
!          return
!      end if
!      call p7_part1partm1_part21part25part25_set_parameter(variable_name(1:l),real_part,imag_part,success)
!      if(success==0) then ! return immediately on error
!          return
!      end if
   end subroutine


   ! BLHA1 interface
   subroutine     OLP_EvalSubProcess(label, momenta, mu,  parameters, res, icheck) &
   & bind(C,name="olp_evalsubprocess_")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int), intent(in) :: label
      real(kind=c_double), intent(in) :: mu
      real(kind=c_double), dimension(50), intent(in) :: momenta
      real(kind=c_double), dimension(10), intent(in) :: parameters
      real(kind=c_double), dimension(60), intent(out) :: res
      integer, intent(inout) :: icheck
      integer(kind=c_int) :: succ
      real(kind=c_double) :: alpha_s
      real(kind=c_double), parameter :: one_over_2pi = 0.15915494309189533577d0

      alpha_s = parameters(1)

      select case(label)
      ! Born 1
      case(0)
              call eval0(momenta(1:25), mu, parameters, res, blha1_mode=.true.)
              res(1:3) = alpha_s * one_over_2pi * res(1:3)
      ! Born 2
      case(1)
              call eval1(momenta(1:25), mu, parameters, res, blha1_mode=.true.)
              res(1:3) = alpha_s * one_over_2pi * res(1:3)
      ! gg -> ghh 1
      case(2)
              call eval2(momenta(1:25), mu, parameters, res, icheck, blha1_mode=.true.)
              res(1:3) = alpha_s * one_over_2pi * res(1:3)
      ! gg -> ghh 2
      case(9)
              call eval9(momenta(1:25), mu, parameters, res, icheck, blha1_mode=.true.)
              res(1:3) = alpha_s * one_over_2pi * res(1:3)
      ! dg -> dhh 1
      case(3)
              call eval3(momenta(1:25), mu, parameters, res, icheck, blha1_mode=.true.)
              res(1:3) = alpha_s * one_over_2pi * res(1:3)
      ! dg -> dhh 2
      case(10)
              call eval10(momenta(1:25), mu, parameters, res, icheck, blha1_mode=.true.)
              res(1:3) = alpha_s * one_over_2pi * res(1:3)
      ! gd -> dhh 1
      case(4)
              call eval4(momenta(1:25), mu, parameters, res, icheck, blha1_mode=.true.)
              res(1:3) = alpha_s * one_over_2pi * res(1:3)
      ! gd -> dhh 2
      case(11)
              call eval11(momenta(1:25), mu, parameters, res, icheck, blha1_mode=.true.)
              res(1:3) = alpha_s * one_over_2pi * res(1:3)
      ! db g -> db h h 1
      case(5)
              call eval5(momenta(1:25), mu, parameters, res, icheck, blha1_mode=.true.)
              res(1:3) = alpha_s * one_over_2pi * res(1:3)
      ! db g -> db h h 2
      case(12)
              call eval12(momenta(1:25), mu, parameters, res, icheck, blha1_mode=.true.)
              res(1:3) = alpha_s * one_over_2pi * res(1:3)
      ! g db -> db h h 1
      case(6)
              call eval6(momenta(1:25), mu, parameters, res, icheck, blha1_mode=.true.)
              res(1:3) = alpha_s * one_over_2pi * res(1:3)
      ! g db -> db h h 2
      case(13)
              call eval13(momenta(1:25), mu, parameters, res, icheck, blha1_mode=.true.)
              res(1:3) = alpha_s * one_over_2pi * res(1:3)
      ! d db -> h h g
      case(7)
              call eval7(momenta(1:25), mu, parameters, res, icheck, blha1_mode=.true.)
              res(1:3) = alpha_s * one_over_2pi * res(1:3)
      ! d db -> h h g 2
      case(14)
              call eval14(momenta(1:25), mu, parameters, res, icheck, blha1_mode=.true.)
              res(1:3) = alpha_s * one_over_2pi * res(1:3)
      ! db d -> h h g 1
      case(8)
              call eval8(momenta(1:25), mu, parameters, res, icheck, blha1_mode=.true.)
              res(1:3) = alpha_s * one_over_2pi * res(1:3)
      ! d db -> h h g 2
      case(15)
              call eval15(momenta(1:25), mu, parameters, res, icheck, blha1_mode=.true.)
              res(1:3) = alpha_s * one_over_2pi * res(1:3)
      case default
         res(:) = 0.0d0
      end select

   end subroutine OLP_EvalSubProcess

   ! BLHA2 interface
   subroutine     OLP_EvalSubProcess2(label, momenta, mu, res, acc, icheck) &
   & bind(C,name="olp_evalsubprocess2_")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int), intent(in) :: label
      real(kind=c_double), intent(in) :: mu
      real(kind=c_double), dimension(50), intent(in) :: momenta
      real(kind=c_double), dimension(60), intent(out) :: res
      real(kind=c_double), intent(out) :: acc
      integer icheck

      real(kind=c_double), dimension(10) :: parameters

      real(kind=c_double), parameter :: one_over_2pi = 0.15915494309189533577d0

      select case(label)
      case(0)
              call eval0(momenta(1:25), mu, parameters, res, acc)
      case(1)
              call eval1(momenta(1:25), mu, parameters, res, acc)
      case(2)
              call eval2(momenta(1:25), mu, parameters, res, icheck, acc)
      case(3)
              call eval3(momenta(1:25), mu, parameters, res, icheck, acc)
      case(4)
              call eval4(momenta(1:25), mu, parameters, res, icheck, acc)
      case(5)
              call eval5(momenta(1:25), mu, parameters, res, icheck, acc)
      case(6)
              call eval6(momenta(1:25), mu, parameters, res, icheck, acc)
      case(7)
              call eval7(momenta(1:25), mu, parameters, res, icheck, acc)
      case(8)
              call eval8(momenta(1:25), mu, parameters, res, icheck, acc)
      case(9)
              call eval9(momenta(1:25), mu, parameters, res, icheck, acc=acc)
      case(10)
              call eval10(momenta(1:25), mu, parameters, res, icheck, acc)
      case(11)
              call eval11(momenta(1:25), mu, parameters, res, icheck, acc)
      case(12)
              call eval12(momenta(1:25), mu, parameters, res, icheck, acc)
      case(13)
              call eval13(momenta(1:25), mu, parameters, res, icheck, acc)
      case(14)
              call eval14(momenta(1:25), mu, parameters, res, icheck, acc)
      case(15)
              call eval15(momenta(1:25), mu, parameters, res, icheck, acc)
      case default
         res(:) = 0.0d0
      end select
   end subroutine OLP_EvalSubProcess2

   subroutine     OLP_Finalize() &
   & bind(C,name="olp_finalize_")
      use, intrinsic :: iso_c_binding
      use p0_part21part21_part25part25_matrix, only: p0_part21part21_part25part25_exitgolem => exitgolem
      use pb_part21part21_part25part25_matrix, only: pb_part21part21_part25part25_exitgolem => exitgolem
      use p2_part21part21_part25part25part21_matrix, only: p2_part21part21_part25part25part21_exitgolem => exitgolem
      use p9_part21part21_part25part25part21_matrix, only: p9_part21part21_part25part25part21_exitgolem => exitgolem
      use p3_part1part21_part25part25part1_matrix, only: p3_part1part21_part25part25part1_exitgolem => exitgolem
      use p10_part1part21_part25part25part1_matrix, only: p10_part1part21_part25part25part1_exitgolem => exitgolem
      implicit none
      call p0_part21part21_part25part25_exitgolem(.true.)
      call pb_part21part21_part25part25_exitgolem(.false.)
      call p2_part21part21_part25part25part21_exitgolem(.false.)
      call p9_part21part21_part25part25part21_exitgolem(.false.)
      call p3_part1part21_part25part25part1_exitgolem(.false.)
      call p10_part1part21_part25part25part1_exitgolem(.false.)
   end subroutine OLP_Finalize

   subroutine     OLP_Option(line,stat) &
   & bind(C,name="olp_option_")
      use, intrinsic :: iso_c_binding
      use p0_part21part21_part25part25_model, only: p0_part21part21_part25part25_parseline => parseline
      use pb_part21part21_part25part25_model, only: pb_part21part21_part25part25_parseline => parseline
      use p2_part21part21_part25part25part21_model, only: p2_part21part21_part25part25part21_parseline => parseline
      use p2_part21part21_part25part25part21_model_qp, only: p2_part21part21_part25part25part21_parseline_qp => parseline
      use p9_part21part21_part25part25part21_model, only: p9_part21part21_part25part25part21_parseline => parseline
      use p9_part21part21_part25part25part21_model_qp, only: p9_part21part21_part25part25part21_parseline_qp => parseline
      use p3_part1part21_part25part25part1_model, only: p3_part1part21_part25part25part1_parseline => parseline
      use p3_part1part21_part25part25part1_model_qp, only: p3_part1part21_part25part25part1_parseline_qp => parseline
      use p10_part1part21_part25part25part1_model, only: p10_part1part21_part25part25part1_parseline => parseline
      use p10_part1part21_part25part25part1_model_qp, only: p10_part1part21_part25part25part1_parseline_qp => parseline
!      use p4_part21part1_part25part25part1_model, only: p4_part25part1_part25part25part1_parseline => parseline
!      use p11_part21part1_part25part25part1_model, only: p11_part25part1_part25part25part1_parseline => parseline
!      use p5_partm1part21_part25part25partm1_model, only: p5_partm1part21_part25part25partm1_parseline => parseline
!      use p12_partm1part21_part25part25partm1_model, only: p12_partm1part21_part25part25partm1_parseline => parseline

!      use p5_partm1part21_partm1part25part25_model, only: p5_partm1part21_partm1part25part25_parseline => parseline
!      use p6_part21partm1_partm1part25part25_model, only: p6_part21partm1_partm1part25part25_parseline => parseline
!      use p7_part1partm1_part21part25part25_model, only: p7_part1partm1_part21part25part25_parseline => parseline

      implicit none
      character(kind=c_char,len=1), intent(in) :: line
      integer(kind=c_int), intent(out) :: stat
      integer :: l, ios

      interface
         function strlen(s) bind(C,name='strlen')
            use, intrinsic :: iso_c_binding
            implicit none
            character(kind=c_char,len=1), intent(in) :: s
            integer(kind=c_int) :: strlen
         end function strlen
      end interface

      l = strlen(line)
      call p0_part21part21_part25part25_parseline(line(1:l),ios)
      if (ios .ne. 0) then
         stat = 0
         return
      end if
      call pb_part21part21_part25part25_parseline(line(1:l),ios)
      if (ios .ne. 0) then
         stat = 0
         return
      end if
      call p2_part21part21_part25part25part21_parseline(line(1:l),ios)
      call p2_part21part21_part25part25part21_parseline_qp(line(1:l),ios)
      if (ios .ne. 0) then
         stat = 0
         return
      end if
      call p9_part21part21_part25part25part21_parseline(line(1:l),ios)
      call p9_part21part21_part25part25part21_parseline_qp(line(1:l),ios)
      if (ios .ne. 0) then
         stat = 0
         return
      end if
      call p3_part1part21_part25part25part1_parseline(line(1:l),ios)
      call p3_part1part21_part25part25part1_parseline_qp(line(1:l),ios)
      if (ios .ne. 0) then
         stat = 0
         return
      end if
      call p10_part1part21_part25part25part1_parseline(line(1:l),ios)
      call p10_part1part21_part25part25part1_parseline_qp(line(1:l),ios)
      if (ios .ne. 0) then
         stat = 0
         return
      end if

!      call p5_partm1part21_partm1part25part25_parseline(line(1:l),ios)
!      if (ios .ne. 0) then
!         stat = 0
!         return
!      end if
!      call p6_part21partm1_partm1part25part25_parseline(line(1:l),ios)
!      if (ios .ne. 0) then
!         stat = 0
!         return
!      end if
!      call p7_part1partm1_part21part25part25_parseline(line(1:l),ios)
!      if (ios .ne. 0) then
!         stat = 0
!         return
!      end if
      stat = 1
   end subroutine OLP_Option
   !---#[ init_event_parameters :
   subroutine     init_event_parameters(sp, parameters)
      use, intrinsic :: iso_c_binding
      implicit none
      integer, intent(in) :: sp
      real(kind=c_double), dimension(10), intent(in) :: parameters
      !
      ! User hook for propagating scale dependent parameters to the
      ! model parameters in the subprocesses.
      !
      ! sp specifies the subprocess
      !
   end subroutine init_event_parameters
   !---#] init_event_parameters :

   !---#[ subroutine eval0 :
   subroutine     eval0(momenta, mu, parameters, res, acc, blha1_mode)
      use, intrinsic :: iso_c_binding
      use p0_part21part21_part25part25_config, only: ki, PSP_chk_th3, nlo_prefactors, PSP_check
      use p0_part21part21_part25part25_model, only: parseline
      use p0_part21part21_part25part25_kinematics, only: boost_to_cms
      use p0_part21part21_part25part25_matrix, only: samplitude, OLP_spin_correlated_lo2, OLP_color_correlated
      use p0_part21part21_part25part25_groups, only: ninja_exit

      implicit none
      real(kind=c_double), dimension(20), intent(in) :: momenta
      real(kind=c_double), intent(in) :: mu
      real(kind=c_double), dimension(10), intent(in) :: parameters
      real(kind=c_double), dimension(60), intent(out) :: res

      real(kind=ki), dimension(4,4) :: vecs
      real(kind=ki), dimension(4) :: amp
      real(kind=c_double), optional :: acc
      logical, optional :: blha1_mode
      real(kind=ki) :: zero
      integer :: i, prec, orig_nlo_prefactors
      logical :: ok

      call init_event_parameters(0, parameters)

      if(present(blha1_mode)) then
         if(blha1_mode) then
            ! save nlo_prefactors and restore later
            orig_nlo_prefactors=nlo_prefactors
            nlo_prefactors=0
         end if
     end if

      vecs(:,1) = real(momenta(1::5),ki)
      vecs(:,2) = real(momenta(2::5),ki)
      vecs(:,3) = real(momenta(3::5),ki)
      vecs(:,4) = real(momenta(4::5),ki)

      call boost_to_cms(vecs)

      call samplitude(vecs, real(mu,ki)*real(mu,ki), amp, prec, ok)
      call ninja_exit()
      if (ok) then
         !
      else
         !
      end if
      if(present(acc)) then
         acc=10.0_ki**(-prec) ! point accuracy
      else
         if(prec.lt.PSP_chk_th3 .and. PSP_check) then
            ! Give back a Nan so that point is discarded
            zero = log(1.0_ki)
            amp(2)= 1.0_ki/zero
        end if
        ! Cannot be assigned if present(acc)=F --> commented out!
        ! acc=1E5_ki ! dummy accuracy which is not used
      end if

      
      res(1) = real(amp(4), c_double)
      res(2) = real(amp(3), c_double)
      res(3) = real(amp(2), c_double)
      res(4) = real(amp(1), c_double)

      if(present(blha1_mode)) then
         if(blha1_mode) then
            ! restore nlo_prefactors
            nlo_prefactors = orig_nlo_prefactors
         end if
     end if

   end subroutine eval0
   !---#] subroutine eval0 :

   !---#[ subroutine eval1 :
   subroutine     eval1(momenta, mu, parameters, res, acc, blha1_mode)
      use, intrinsic :: iso_c_binding
      use pb_part21part21_part25part25_config, only: ki, PSP_chk_th3, nlo_prefactors, PSP_check
      use pb_part21part21_part25part25_model, only: parseline
      use pb_part21part21_part25part25_kinematics, only: boost_to_cms
      use pb_part21part21_part25part25_matrix, only: samplitude, OLP_spin_correlated_lo2, OLP_color_correlated
      use pb_part21part21_part25part25_groups, only: ninja_exit

      implicit none
      real(kind=c_double), dimension(20), intent(in) :: momenta
      real(kind=c_double), intent(in) :: mu
      real(kind=c_double), dimension(10), intent(in) :: parameters
      real(kind=c_double), dimension(60), intent(out) :: res

      real(kind=ki), dimension(4,4) :: vecs
      real(kind=ki), dimension(4) :: amp
      real(kind=c_double), optional :: acc
      logical, optional :: blha1_mode
      real(kind=ki) :: zero
      integer :: i, prec, orig_nlo_prefactors
      logical :: ok

      call init_event_parameters(0, parameters)

      if(present(blha1_mode)) then
         if(blha1_mode) then
            ! save nlo_prefactors and restore later
            orig_nlo_prefactors=nlo_prefactors
            nlo_prefactors=0
         end if
     end if

      vecs(:,1) = real(momenta(1::5),ki)
      vecs(:,2) = real(momenta(2::5),ki)
      vecs(:,3) = real(momenta(3::5),ki)
      vecs(:,4) = real(momenta(4::5),ki)

      call boost_to_cms(vecs)

      call samplitude(vecs, real(mu,ki)*real(mu,ki), amp, prec, ok)
      call ninja_exit()
      if (ok) then
         !
      else
         !
      end if
      if(present(acc)) then
         acc=10.0_ki**(-prec) ! point accuracy
      else
         if(prec.lt.PSP_chk_th3 .and. PSP_check) then
            ! Give back a Nan so that point is discarded
            zero = log(1.0_ki)
            amp(2)= 1.0_ki/zero
        end if
        ! Cannot be assigned if present(acc)=F --> commented out!
        ! acc=1E5_ki ! dummy accuracy which is not used
      end if

      
      res(1) = real(amp(4), c_double)
      res(2) = real(amp(3), c_double)
      res(3) = real(amp(2), c_double)
      res(4) = real(amp(1), c_double)

      if(present(blha1_mode)) then
         if(blha1_mode) then
            ! restore nlo_prefactors
            nlo_prefactors = orig_nlo_prefactors
         end if
     end if

   end subroutine eval1
   !---#] subroutine eval1 :

   !---#[ subroutine eval2 :
   subroutine     eval2(momenta, mu, parameters, res, icheck, acc, blha1_mode)
      use, intrinsic :: iso_c_binding
      use p2_part21part21_part25part25part21_config, only: ki, PSP_chk_li3, nlo_prefactors, PSP_check
      use p2_part21part21_part25part25part21_model, only: parseline
      use p2_part21part21_part25part25part21_kinematics, only: boost_to_cms
      use p2_part21part21_part25part25part21_matrix, only: samplitude, OLP_spin_correlated_lo2, OLP_color_correlated
      use p2_part21part21_part25part25part21_groups, only: ninja_exit, quadninja_exit

      implicit none
      real(kind=c_double), dimension(25), intent(in) :: momenta
      real(kind=c_double), intent(in) :: mu
      real(kind=c_double), dimension(10), intent(in) :: parameters
      real(kind=c_double), dimension(60), intent(out) :: res
      integer icheck

      real(kind=ki), dimension(5,4) :: vecs
      real(kind=ki), dimension(4) :: amp
      real(kind=c_double), optional :: acc
      logical, optional :: blha1_mode
      real(kind=ki) :: zero
      integer :: i, prec, orig_nlo_prefactors
      logical :: ok

      call init_event_parameters(0, parameters)

      if(present(blha1_mode)) then
         if(blha1_mode) then
            ! save nlo_prefactors and restore later
            orig_nlo_prefactors=nlo_prefactors
            nlo_prefactors=0
         end if
     end if

      vecs(:,1) = real(momenta(1::5),ki)
      vecs(:,2) = real(momenta(2::5),ki)
      vecs(:,3) = real(momenta(3::5),ki)
      vecs(:,4) = real(momenta(4::5),ki)

      call boost_to_cms(vecs)

      call samplitude(vecs, real(mu,ki)*real(mu,ki), amp, prec, icheck, ok)
      call ninja_exit()
      call quadninja_exit()
      if (ok) then
         !
      else
         !
      end if
      if(present(acc)) then
         acc=10.0_ki**(-prec) ! point accuracy
      else
         if(prec.lt.PSP_chk_li3 .and. PSP_check) then
            ! Give back a Nan so that point is discarded
            zero = log(1.0_ki)
            amp(2)= 1.0_ki/zero
        end if
        ! Cannot be assigned if present(acc)=F --> commented out!
        ! acc=1E5_ki ! dummy accuracy which is not used
      end if


      res(1) = real(amp(4), c_double)
      res(2) = real(amp(3), c_double)
      res(3) = real(amp(2), c_double)
      res(4) = real(amp(1), c_double)

      if(present(blha1_mode)) then
         if(blha1_mode) then
            ! restore nlo_prefactors
            nlo_prefactors = orig_nlo_prefactors
         end if
     end if

   end subroutine eval2
   !---#] subroutine eval2 :
   !---#[ subroutine eval3 :
   subroutine     eval3(momenta, mu, parameters, res, icheck, acc, blha1_mode)
      use, intrinsic :: iso_c_binding
      use p3_part1part21_part25part25part1_config, only: ki, PSP_chk_li3, nlo_prefactors, PSP_check
      use p3_part1part21_part25part25part1_model, only: parseline
      use p3_part1part21_part25part25part1_kinematics, only: boost_to_cms
      use p3_part1part21_part25part25part1_matrix, only: samplitude, OLP_spin_correlated_lo2, OLP_color_correlated
      use p3_part1part21_part25part25part1_groups, only: ninja_exit, quadninja_exit

      implicit none
      real(kind=c_double), dimension(25), intent(in) :: momenta
      real(kind=c_double), intent(in) :: mu
      real(kind=c_double), dimension(10), intent(in) :: parameters
      real(kind=c_double), dimension(60), intent(out) :: res
      integer icheck

      real(kind=ki), dimension(5,4) :: vecs
      real(kind=ki), dimension(4) :: amp
      real(kind=c_double), optional :: acc
      logical, optional :: blha1_mode
      real(kind=ki) :: zero
      integer :: i, prec, orig_nlo_prefactors
      logical :: ok

      call init_event_parameters(0, parameters)

      if(present(blha1_mode)) then
         if(blha1_mode) then
            ! save nlo_prefactors and restore later
            orig_nlo_prefactors=nlo_prefactors
            nlo_prefactors=0
         end if
     end if

      vecs(:,1) = real(momenta(1::5),ki)
      vecs(:,2) = real(momenta(2::5),ki)
      vecs(:,3) = real(momenta(3::5),ki)
      vecs(:,4) = real(momenta(4::5),ki)

      call boost_to_cms(vecs)

      call samplitude(vecs, real(mu,ki)*real(mu,ki), amp, prec, icheck, ok)
      call ninja_exit()
      call quadninja_exit()
      if (ok) then
         !
      else
         !
      end if
      if(present(acc)) then
         acc=10.0_ki**(-prec) ! point accuracy
      else
         if(prec.lt.PSP_chk_li3 .and. PSP_check) then
            ! Give back a Nan so that point is discarded
            zero = log(1.0_ki)
            amp(2)= 1.0_ki/zero
        end if
        ! Cannot be assigned if present(acc)=F --> commented out!
        ! acc=1E5_ki ! dummy accuracy which is not used
      end if


      res(1) = real(amp(4), c_double)
      res(2) = real(amp(3), c_double)
      res(3) = real(amp(2), c_double)
      res(4) = real(amp(1), c_double)

      if(present(blha1_mode)) then
         if(blha1_mode) then
            ! restore nlo_prefactors
            nlo_prefactors = orig_nlo_prefactors
         end if
     end if

   end subroutine eval3
   !---#] subroutine eval3 :
   !---#[ subroutine eval4 :
   subroutine     eval4(momenta, mu, parameters, res, icheck, acc, blha1_mode)
      use, intrinsic :: iso_c_binding
      use p3_part1part21_part25part25part1_config, only: ki, PSP_chk_li3, nlo_prefactors, PSP_check
      use p3_part1part21_part25part25part1_model, only: parseline
      use p3_part1part21_part25part25part1_kinematics, only: boost_to_cms
      use p4_part21part1_part25part25part1_matrix, only: samplitude, OLP_spin_correlated_lo2, OLP_color_correlated
      use p3_part1part21_part25part25part1_groups, only: ninja_exit, quadninja_exit

      implicit none
      real(kind=c_double), dimension(25), intent(in) :: momenta
      real(kind=c_double), intent(in) :: mu
      real(kind=c_double), dimension(10), intent(in) :: parameters
      real(kind=c_double), dimension(60), intent(out) :: res
      integer icheck

      real(kind=ki), dimension(5,4) :: vecs
      real(kind=ki), dimension(4) :: amp
      real(kind=c_double), optional :: acc
      logical, optional :: blha1_mode
      real(kind=ki) :: zero
      integer :: i, prec, orig_nlo_prefactors
      logical :: ok

      call init_event_parameters(2, parameters)

      if(present(blha1_mode)) then
         if(blha1_mode) then
            ! save nlo_prefactors and restore later
            orig_nlo_prefactors=nlo_prefactors
            nlo_prefactors=0
         end if
     end if

      vecs(:,1) = real(momenta(1::5),ki)
      vecs(:,2) = real(momenta(2::5),ki)
      vecs(:,3) = real(momenta(3::5),ki)
      vecs(:,4) = real(momenta(4::5),ki)

      call boost_to_cms(vecs)

      call samplitude(vecs, real(mu,ki)*real(mu,ki), amp, prec, icheck, ok)
      call ninja_exit()
      call quadninja_exit()
      if (ok) then
         !
      else
         !
      end if
      if(present(acc)) then
         acc=10.0_ki**(-prec) ! point accuracy
      else
         if(prec.lt.PSP_chk_li3 .and. PSP_check) then
            ! Give back a Nan so that point is discarded
            zero = log(1.0_ki)
            amp(2)= 1.0_ki/zero
        end if
        ! Cannot be assigned if present(acc)=F --> commented out!
        ! acc=1E5_ki ! dummy accuracy which is not used
      end if


      res(1) = real(amp(4), c_double)
      res(2) = real(amp(3), c_double)
      res(3) = real(amp(2), c_double)
      res(4) = real(amp(1), c_double)

      if(present(blha1_mode)) then
         if(blha1_mode) then
            ! restore nlo_prefactors
            nlo_prefactors = orig_nlo_prefactors
         end if
     end if

   end subroutine eval4
   !---#] subroutine eval4 :
   !---#[ subroutine eval5 :
   subroutine     eval5(momenta, mu, parameters, res, icheck, acc, blha1_mode)
      use, intrinsic :: iso_c_binding
      use p3_part1part21_part25part25part1_config, only: ki, PSP_chk_li3, nlo_prefactors, PSP_check
      use p3_part1part21_part25part25part1_model, only: parseline
      use p3_part1part21_part25part25part1_kinematics, only: boost_to_cms
      use p5_partm1part21_part25part25partm1_matrix, only: samplitude, OLP_spin_correlated_lo2, OLP_color_correlated
      use p3_part1part21_part25part25part1_groups, only: ninja_exit, quadninja_exit

      implicit none
      real(kind=c_double), dimension(25), intent(in) :: momenta
      real(kind=c_double), intent(in) :: mu
      real(kind=c_double), dimension(10), intent(in) :: parameters
      real(kind=c_double), dimension(60), intent(out) :: res
      integer icheck

      real(kind=ki), dimension(5,4) :: vecs
      real(kind=ki), dimension(4) :: amp
      real(kind=c_double), optional :: acc
      logical, optional :: blha1_mode
      real(kind=ki) :: zero
      integer :: i, prec, orig_nlo_prefactors
      logical :: ok

      call init_event_parameters(2, parameters)

      if(present(blha1_mode)) then
         if(blha1_mode) then
            ! save nlo_prefactors and restore later
            orig_nlo_prefactors=nlo_prefactors
            nlo_prefactors=0
         end if
     end if

      vecs(:,1) = real(momenta(1::5),ki)
      vecs(:,2) = real(momenta(2::5),ki)
      vecs(:,3) = real(momenta(3::5),ki)
      vecs(:,4) = real(momenta(4::5),ki)

      call boost_to_cms(vecs)

      call samplitude(vecs, real(mu,ki)*real(mu,ki), amp, prec, icheck, ok)
      call ninja_exit()
      call quadninja_exit()
      if (ok) then
         !
      else
         !
      end if
      if(present(acc)) then
         acc=10.0_ki**(-prec) ! point accuracy
!      else
!         if(prec.lt.PSP_chk_li3 .and. PSP_check) then
!            ! Give back a Nan so that point is discarded
!            zero = log(1.0_ki)
!            amp(2)= 1.0_ki/zero
!        end if
        ! Cannot be assigned if present(acc)=F --> commented out!
        ! acc=1E5_ki ! dummy accuracy which is not used
      end if


      res(1) = real(amp(4), c_double)
      res(2) = real(amp(3), c_double)
      res(3) = real(amp(2), c_double)
      res(4) = real(amp(1), c_double)

      if(present(blha1_mode)) then
         if(blha1_mode) then
            ! restore nlo_prefactors
            nlo_prefactors = orig_nlo_prefactors
         end if
     end if

   end subroutine eval5
   !---#] subroutine eval5 :
   !---#[ subroutine eval6 :
   subroutine     eval6(momenta, mu, parameters, res, icheck, acc, blha1_mode)
      use, intrinsic :: iso_c_binding
      use p3_part1part21_part25part25part1_config, only: ki, PSP_chk_li3, nlo_prefactors, PSP_check
      use p3_part1part21_part25part25part1_model, only: parseline
      use p3_part1part21_part25part25part1_kinematics, only: boost_to_cms
      use p6_part21partm1_part25part25partm1_matrix, only: samplitude, OLP_spin_correlated_lo2, OLP_color_correlated
      use p3_part1part21_part25part25part1_groups, only: ninja_exit, quadninja_exit

      implicit none
      real(kind=c_double), dimension(25), intent(in) :: momenta
      real(kind=c_double), intent(in) :: mu
      real(kind=c_double), dimension(10), intent(in) :: parameters
      real(kind=c_double), dimension(60), intent(out) :: res
      integer icheck

      real(kind=ki), dimension(5,4) :: vecs
      real(kind=ki), dimension(4) :: amp
      real(kind=c_double), optional :: acc
      logical, optional :: blha1_mode
      real(kind=ki) :: zero
      integer :: i, prec, orig_nlo_prefactors
      logical :: ok

      call init_event_parameters(2, parameters)

      if(present(blha1_mode)) then
         if(blha1_mode) then
            ! save nlo_prefactors and restore later
            orig_nlo_prefactors=nlo_prefactors
            nlo_prefactors=0
         end if
     end if

      vecs(:,1) = real(momenta(1::5),ki)
      vecs(:,2) = real(momenta(2::5),ki)
      vecs(:,3) = real(momenta(3::5),ki)
      vecs(:,4) = real(momenta(4::5),ki)

      call boost_to_cms(vecs)

      call samplitude(vecs, real(mu,ki)*real(mu,ki), amp, prec, icheck, ok)
      call ninja_exit()
      call quadninja_exit()
      if (ok) then
         !
      else
         !
      end if
      if(present(acc)) then
         acc=10.0_ki**(-prec) ! point accuracy
      else
         if(prec.lt.PSP_chk_li3 .and. PSP_check) then
            ! Give back a Nan so that point is discarded
            zero = log(1.0_ki)
            amp(2)= 1.0_ki/zero
        end if
        ! Cannot be assigned if present(acc)=F --> commented out!
        ! acc=1E5_ki ! dummy accuracy which is not used
      end if


      res(1) = real(amp(4), c_double)
      res(2) = real(amp(3), c_double)
      res(3) = real(amp(2), c_double)
      res(4) = real(amp(1), c_double)

      if(present(blha1_mode)) then
         if(blha1_mode) then
            ! restore nlo_prefactors
            nlo_prefactors = orig_nlo_prefactors
         end if
     end if

   end subroutine eval6
   !---#] subroutine eval6 :
   !---#[ subroutine eval7 :
   subroutine     eval7(momenta, mu, parameters, res, icheck, acc, blha1_mode)
      use, intrinsic :: iso_c_binding
      use p3_part1part21_part25part25part1_config, only: ki, PSP_chk_li3, nlo_prefactors, PSP_check
      use p3_part1part21_part25part25part1_model, only: parseline
      use p3_part1part21_part25part25part1_kinematics, only: boost_to_cms
      use p7_part1partm1_part25part25part21_matrix, only: samplitude, OLP_spin_correlated_lo2, OLP_color_correlated
      use p3_part1part21_part25part25part1_groups, only: ninja_exit, quadninja_exit

      implicit none
      real(kind=c_double), dimension(25), intent(in) :: momenta
      real(kind=c_double), intent(in) :: mu
      real(kind=c_double), dimension(10), intent(in) :: parameters
      real(kind=c_double), dimension(60), intent(out) :: res
      integer icheck

      real(kind=ki), dimension(5,4) :: vecs
      real(kind=ki), dimension(4) :: amp
      real(kind=c_double), optional :: acc
      logical, optional :: blha1_mode
      real(kind=ki) :: zero
      integer :: i, prec, orig_nlo_prefactors
      logical :: ok

      call init_event_parameters(2, parameters)

      if(present(blha1_mode)) then
         if(blha1_mode) then
            ! save nlo_prefactors and restore later
            orig_nlo_prefactors=nlo_prefactors
            nlo_prefactors=0
         end if
     end if

      vecs(:,1) = real(momenta(1::5),ki)
      vecs(:,2) = real(momenta(2::5),ki)
      vecs(:,3) = real(momenta(3::5),ki)
      vecs(:,4) = real(momenta(4::5),ki)

      call boost_to_cms(vecs)

      call samplitude(vecs, real(mu,ki)*real(mu,ki), amp, prec, icheck, ok)
      call ninja_exit()
      call quadninja_exit()
      if (ok) then
         !
      else
         !
      end if
      if(present(acc)) then
         acc=10.0_ki**(-prec) ! point accuracy
      else
         if(prec.lt.PSP_chk_li3 .and. PSP_check) then
            ! Give back a Nan so that point is discarded
            zero = log(1.0_ki)
            amp(2)= 1.0_ki/zero
        end if
        ! Cannot be assigned if present(acc)=F --> commented out!
        ! acc=1E5_ki ! dummy accuracy which is not used
      end if


      res(1) = real(amp(4), c_double)
      res(2) = real(amp(3), c_double)
      res(3) = real(amp(2), c_double)
      res(4) = real(amp(1), c_double)

      if(present(blha1_mode)) then
         if(blha1_mode) then
            ! restore nlo_prefactors
            nlo_prefactors = orig_nlo_prefactors
         end if
     end if

   end subroutine eval7
   !---#] subroutine eval7 :
   !---#[ subroutine eval8 :
   subroutine     eval8(momenta, mu, parameters, res, icheck, acc, blha1_mode)
      use, intrinsic :: iso_c_binding
      use p3_part1part21_part25part25part1_config, only: ki, PSP_chk_li3, nlo_prefactors, PSP_check
      use p3_part1part21_part25part25part1_model, only: parseline
      use p3_part1part21_part25part25part1_kinematics, only: boost_to_cms
      use p7_part1partm1_part25part25part21_matrix, only: samplitude, OLP_spin_correlated_lo2, OLP_color_correlated
      use p3_part1part21_part25part25part1_groups, only: ninja_exit, quadninja_exit

      implicit none
      real(kind=c_double), dimension(25), intent(in) :: momenta
      real(kind=c_double), intent(in) :: mu
      real(kind=c_double), dimension(10), intent(in) :: parameters
      real(kind=c_double), dimension(60), intent(out) :: res
      integer icheck

      real(kind=ki), dimension(5,4) :: vecs
      real(kind=ki), dimension(4) :: amp
      real(kind=c_double), optional :: acc
      logical, optional :: blha1_mode
      real(kind=ki) :: zero
      integer :: i, prec, orig_nlo_prefactors
      logical :: ok

      call init_event_parameters(2, parameters)

      if(present(blha1_mode)) then
         if(blha1_mode) then
            ! save nlo_prefactors and restore later
            orig_nlo_prefactors=nlo_prefactors
            nlo_prefactors=0
         end if
     end if

      vecs(:,1) = real(momenta(1::5),ki)
      vecs(:,2) = real(momenta(2::5),ki)
      vecs(:,3) = real(momenta(3::5),ki)
      vecs(:,4) = real(momenta(4::5),ki)

      call boost_to_cms(vecs)

      call samplitude(vecs, real(mu,ki)*real(mu,ki), amp, prec, icheck, ok)
      call ninja_exit()
      call quadninja_exit()
      if (ok) then
         !
      else
         !
      end if
      if(present(acc)) then
         acc=10.0_ki**(-prec) ! point accuracy
      else
         if(prec.lt.PSP_chk_li3 .and. PSP_check) then
            ! Give back a Nan so that point is discarded
            zero = log(1.0_ki)
            amp(2)= 1.0_ki/zero
        end if
        ! Cannot be assigned if present(acc)=F --> commented out!
        ! acc=1E5_ki ! dummy accuracy which is not used
      end if


      res(1) = real(amp(4), c_double)
      res(2) = real(amp(3), c_double)
      res(3) = real(amp(2), c_double)
      res(4) = real(amp(1), c_double)

      if(present(blha1_mode)) then
         if(blha1_mode) then
            ! restore nlo_prefactors
            nlo_prefactors = orig_nlo_prefactors
         end if
     end if

   end subroutine eval8
   !---#] subroutine eval8 :
   !---#[ subroutine eval9 :
   subroutine     eval9(momenta, mu, parameters, res, ichecked, acc, blha1_mode)
      use, intrinsic :: iso_c_binding
      use p9_part21part21_part25part25part21_config, only: ki, PSP_chk_li3, nlo_prefactors, PSP_check
      use p9_part21part21_part25part25part21_model, only: parseline
      use p9_part21part21_part25part25part21_kinematics, only: boost_to_cms
      use p9_part21part21_part25part25part21_matrix, only: samplitude, OLP_spin_correlated_lo2, OLP_color_correlated
      use p9_part21part21_part25part25part21_groups, only: ninja_exit, quadninja_exit

      implicit none
      real(kind=c_double), dimension(25), intent(in) :: momenta
      real(kind=c_double), intent(in) :: mu
      real(kind=c_double), dimension(10), intent(in) :: parameters
      real(kind=c_double), dimension(60), intent(out) :: res
      integer, intent(in) :: ichecked

      real(kind=ki), dimension(5,4) :: vecs
      real(kind=ki), dimension(4) :: amp
      real(kind=c_double), optional :: acc
      logical, optional :: blha1_mode
      real(kind=ki) :: zero
      integer :: i, prec, orig_nlo_prefactors
      logical :: ok

      call init_event_parameters(0, parameters)

      if(present(blha1_mode)) then
         if(blha1_mode) then
            ! save nlo_prefactors and restore later
            orig_nlo_prefactors=nlo_prefactors
            nlo_prefactors=0
         end if
     end if

      vecs(:,1) = real(momenta(1::5),ki)
      vecs(:,2) = real(momenta(2::5),ki)
      vecs(:,3) = real(momenta(3::5),ki)
      vecs(:,4) = real(momenta(4::5),ki)

      call boost_to_cms(vecs)

      call samplitude(vecs, real(mu,ki)*real(mu,ki), amp, prec, ichecked, ok)
      call ninja_exit()
      call quadninja_exit()
      if (ok) then
         !
      else
         !
      end if
      if(present(acc)) then
         acc=10.0_ki**(-prec) ! point accuracy
      else
         if(prec.lt.PSP_chk_li3 .and. PSP_check) then
            ! Give back a Nan so that point is discarded
            zero = log(1.0_ki)
            amp(2)= 1.0_ki/zero
        end if
        ! Cannot be assigned if present(acc)=F --> commented out!
        ! acc=1E5_ki ! dummy accuracy which is not used
      end if


      res(1) = real(amp(4), c_double)
      res(2) = real(amp(3), c_double)
      res(3) = real(amp(2), c_double)
      res(4) = real(amp(1), c_double)

      if(present(blha1_mode)) then
         if(blha1_mode) then
            ! restore nlo_prefactors
            nlo_prefactors = orig_nlo_prefactors
         end if
     end if

   end subroutine eval9
   !---#] subroutine eval9 :
   !---#[ subroutine eval10 :
   subroutine     eval10(momenta, mu, parameters, res, ichecked, acc, blha1_mode)
      use, intrinsic :: iso_c_binding
      use p10_part1part21_part25part25part1_config, only: ki, PSP_chk_li3, nlo_prefactors, PSP_check
      use p10_part1part21_part25part25part1_model, only: parseline
      use p10_part1part21_part25part25part1_kinematics, only: boost_to_cms
      use p10_part1part21_part25part25part1_matrix, only: samplitude, OLP_spin_correlated_lo2, OLP_color_correlated
      use p10_part1part21_part25part25part1_groups, only: ninja_exit, quadninja_exit

      implicit none
      real(kind=c_double), dimension(25), intent(in) :: momenta
      real(kind=c_double), intent(in) :: mu
      real(kind=c_double), dimension(10), intent(in) :: parameters
      real(kind=c_double), dimension(60), intent(out) :: res
      integer, intent(in) :: ichecked

      real(kind=ki), dimension(5,4) :: vecs
      real(kind=ki), dimension(4) :: amp
      real(kind=c_double), optional :: acc
      logical, optional :: blha1_mode
      real(kind=ki) :: zero
      integer :: i, prec, orig_nlo_prefactors
      logical :: ok

      call init_event_parameters(0, parameters)

      if(present(blha1_mode)) then
         if(blha1_mode) then
            ! save nlo_prefactors and restore later
            orig_nlo_prefactors=nlo_prefactors
            nlo_prefactors=0
         end if
     end if

      vecs(:,1) = real(momenta(1::5),ki)
      vecs(:,2) = real(momenta(2::5),ki)
      vecs(:,3) = real(momenta(3::5),ki)
      vecs(:,4) = real(momenta(4::5),ki)

      call boost_to_cms(vecs)

      call samplitude(vecs, real(mu,ki)*real(mu,ki), amp, prec, ichecked, ok)
      call ninja_exit()
      call quadninja_exit()
      if (ok) then
         !
      else
         !
      end if
      if(present(acc)) then
         acc=10.0_ki**(-prec) ! point accuracy
      else
         if(prec.lt.PSP_chk_li3 .and. PSP_check) then
            ! Give back a Nan so that point is discarded
            zero = log(1.0_ki)
            amp(2)= 1.0_ki/zero
        end if
        ! Cannot be assigned if present(acc)=F --> commented out!
        ! acc=1E5_ki ! dummy accuracy which is not used
      end if


      res(1) = real(amp(4), c_double)
      res(2) = real(amp(3), c_double)
      res(3) = real(amp(2), c_double)
      res(4) = real(amp(1), c_double)

      if(present(blha1_mode)) then
         if(blha1_mode) then
            ! restore nlo_prefactors
            nlo_prefactors = orig_nlo_prefactors
         end if
     end if

   end subroutine eval10
   !---#] subroutine eval10 :
   !---#[ subroutine eval11 :
   subroutine     eval11(momenta, mu, parameters, res, ichecked, acc, blha1_mode)
      use, intrinsic :: iso_c_binding
      use p10_part1part21_part25part25part1_config, only: ki, PSP_chk_li3, nlo_prefactors, PSP_check
      use p10_part1part21_part25part25part1_model, only: parseline
      use p10_part1part21_part25part25part1_kinematics, only: boost_to_cms
      use p11_part21part1_part25part25part1_matrix, only: samplitude, OLP_spin_correlated_lo2, OLP_color_correlated
      use p10_part1part21_part25part25part1_groups, only: ninja_exit, quadninja_exit

      implicit none
      real(kind=c_double), dimension(25), intent(in) :: momenta
      real(kind=c_double), intent(in) :: mu
      real(kind=c_double), dimension(10), intent(in) :: parameters
      real(kind=c_double), dimension(60), intent(out) :: res
      integer, intent(in) :: ichecked

      real(kind=ki), dimension(5,4) :: vecs
      real(kind=ki), dimension(4) :: amp
      real(kind=c_double), optional :: acc
      logical, optional :: blha1_mode
      real(kind=ki) :: zero
      integer :: i, prec, orig_nlo_prefactors
      logical :: ok

      call init_event_parameters(2, parameters)

      if(present(blha1_mode)) then
         if(blha1_mode) then
            ! save nlo_prefactors and restore later
            orig_nlo_prefactors=nlo_prefactors
            nlo_prefactors=0
         end if
     end if

      vecs(:,1) = real(momenta(1::5),ki)
      vecs(:,2) = real(momenta(2::5),ki)
      vecs(:,3) = real(momenta(3::5),ki)
      vecs(:,4) = real(momenta(4::5),ki)

      call boost_to_cms(vecs)

      call samplitude(vecs, real(mu,ki)*real(mu,ki), amp, prec, ichecked, ok)
      call ninja_exit()
      call quadninja_exit()
      if (ok) then
         !
      else
         !
      end if
      if(present(acc)) then
         acc=10.0_ki**(-prec) ! point accuracy
      else
         if(prec.lt.PSP_chk_li3 .and. PSP_check) then
            ! Give back a Nan so that point is discarded
            zero = log(1.0_ki)
            amp(2)= 1.0_ki/zero
        end if
        ! Cannot be assigned if present(acc)=F --> commented out!
        ! acc=1E5_ki ! dummy accuracy which is not used
      end if


      res(1) = real(amp(4), c_double)
      res(2) = real(amp(3), c_double)
      res(3) = real(amp(2), c_double)
      res(4) = real(amp(1), c_double)

      if(present(blha1_mode)) then
         if(blha1_mode) then
            ! restore nlo_prefactors
            nlo_prefactors = orig_nlo_prefactors
         end if
     end if

   end subroutine eval11
   !---#] subroutine eval11 :
   !---#[ subroutine eval12 :
   subroutine     eval12(momenta, mu, parameters, res, ichecked, acc, blha1_mode)
      use, intrinsic :: iso_c_binding
      use p10_part1part21_part25part25part1_config, only: ki, PSP_chk_li3, nlo_prefactors, PSP_check
      use p10_part1part21_part25part25part1_model, only: parseline
      use p10_part1part21_part25part25part1_kinematics, only: boost_to_cms
      use p12_partm1part21_part25part25partm1_matrix, only: samplitude, OLP_spin_correlated_lo2, OLP_color_correlated
      use p10_part1part21_part25part25part1_groups, only: ninja_exit, quadninja_exit

      implicit none
      real(kind=c_double), dimension(25), intent(in) :: momenta
      real(kind=c_double), intent(in) :: mu
      real(kind=c_double), dimension(10), intent(in) :: parameters
      real(kind=c_double), dimension(60), intent(out) :: res
      integer, intent(in) :: ichecked

      real(kind=ki), dimension(5,4) :: vecs
      real(kind=ki), dimension(4) :: amp
      real(kind=c_double), optional :: acc
      logical, optional :: blha1_mode
      real(kind=ki) :: zero
      integer :: i, prec, orig_nlo_prefactors
      logical :: ok

      call init_event_parameters(2, parameters)

      if(present(blha1_mode)) then
         if(blha1_mode) then
            ! save nlo_prefactors and restore later
            orig_nlo_prefactors=nlo_prefactors
            nlo_prefactors=0
         end if
     end if

      vecs(:,1) = real(momenta(1::5),ki)
      vecs(:,2) = real(momenta(2::5),ki)
      vecs(:,3) = real(momenta(3::5),ki)
      vecs(:,4) = real(momenta(4::5),ki)

      call boost_to_cms(vecs)

      call samplitude(vecs, real(mu,ki)*real(mu,ki), amp, prec, ichecked, ok)
      call ninja_exit()
      call quadninja_exit()
      if (ok) then
         !
      else
         !
      end if
      if(present(acc)) then
         acc=10.0_ki**(-prec) ! point accuracy
      else
         if(prec.lt.PSP_chk_li3 .and. PSP_check) then
            ! Give back a Nan so that point is discarded
            zero = log(1.0_ki)
            amp(2)= 1.0_ki/zero
        end if
        ! Cannot be assigned if present(acc)=F --> commented out!
        ! acc=1E5_ki ! dummy accuracy which is not used
      end if


      res(1) = real(amp(4), c_double)
      res(2) = real(amp(3), c_double)
      res(3) = real(amp(2), c_double)
      res(4) = real(amp(1), c_double)

      if(present(blha1_mode)) then
         if(blha1_mode) then
            ! restore nlo_prefactors
            nlo_prefactors = orig_nlo_prefactors
         end if
     end if

   end subroutine eval12
   !---#] subroutine eval12 :
   !---#[ subroutine eval13 :
   subroutine     eval13(momenta, mu, parameters, res, ichecked, acc, blha1_mode)
      use, intrinsic :: iso_c_binding
      use p10_part1part21_part25part25part1_config, only: ki, PSP_chk_li3, nlo_prefactors, PSP_check
      use p10_part1part21_part25part25part1_model, only: parseline
      use p10_part1part21_part25part25part1_kinematics, only: boost_to_cms
      use p13_part21partm1_part25part25partm1_matrix, only: samplitude, OLP_spin_correlated_lo2, OLP_color_correlated
      use p10_part1part21_part25part25part1_groups, only: ninja_exit, quadninja_exit

      implicit none
      real(kind=c_double), dimension(25), intent(in) :: momenta
      real(kind=c_double), intent(in) :: mu
      real(kind=c_double), dimension(10), intent(in) :: parameters
      real(kind=c_double), dimension(60), intent(out) :: res
      integer, intent(in) :: ichecked

      real(kind=ki), dimension(5,4) :: vecs
      real(kind=ki), dimension(4) :: amp
      real(kind=c_double), optional :: acc
      logical, optional :: blha1_mode
      real(kind=ki) :: zero
      integer :: i, prec, orig_nlo_prefactors
      logical :: ok

      call init_event_parameters(2, parameters)

      if(present(blha1_mode)) then
         if(blha1_mode) then
            ! save nlo_prefactors and restore later
            orig_nlo_prefactors=nlo_prefactors
            nlo_prefactors=0
         end if
     end if

      vecs(:,1) = real(momenta(1::5),ki)
      vecs(:,2) = real(momenta(2::5),ki)
      vecs(:,3) = real(momenta(3::5),ki)
      vecs(:,4) = real(momenta(4::5),ki)

      call boost_to_cms(vecs)

      call samplitude(vecs, real(mu,ki)*real(mu,ki), amp, prec, ichecked, ok)
      call ninja_exit()
      call quadninja_exit()
      if (ok) then
         !
      else
         !
      end if
      if(present(acc)) then
         acc=10.0_ki**(-prec) ! point accuracy
      else
         if(prec.lt.PSP_chk_li3 .and. PSP_check) then
            ! Give back a Nan so that point is discarded
            zero = log(1.0_ki)
            amp(2)= 1.0_ki/zero
        end if
        ! Cannot be assigned if present(acc)=F --> commented out!
        ! acc=1E5_ki ! dummy accuracy which is not used
      end if


      res(1) = real(amp(4), c_double)
      res(2) = real(amp(3), c_double)
      res(3) = real(amp(2), c_double)
      res(4) = real(amp(1), c_double)

      if(present(blha1_mode)) then
         if(blha1_mode) then
            ! restore nlo_prefactors
            nlo_prefactors = orig_nlo_prefactors
         end if
     end if

   end subroutine eval13
   !---#] subroutine eval13 :
   !---#[ subroutine eval14 :
   subroutine     eval14(momenta, mu, parameters, res, ichecked, acc, blha1_mode)
      use, intrinsic :: iso_c_binding
      use p10_part1part21_part25part25part1_config, only: ki, PSP_chk_li3, nlo_prefactors, PSP_check
      use p10_part1part21_part25part25part1_model, only: parseline
      use p10_part1part21_part25part25part1_kinematics, only: boost_to_cms
      use p14_part1partm1_part25part25part21_matrix, only: samplitude, OLP_spin_correlated_lo2, OLP_color_correlated
      use p10_part1part21_part25part25part1_groups, only: ninja_exit, quadninja_exit

      implicit none
      real(kind=c_double), dimension(25), intent(in) :: momenta
      real(kind=c_double), intent(in) :: mu
      real(kind=c_double), dimension(10), intent(in) :: parameters
      real(kind=c_double), dimension(60), intent(out) :: res
      integer, intent(in) :: ichecked

      real(kind=ki), dimension(5,4) :: vecs
      real(kind=ki), dimension(4) :: amp
      real(kind=c_double), optional :: acc
      logical, optional :: blha1_mode
      real(kind=ki) :: zero
      integer :: i, prec, orig_nlo_prefactors
      logical :: ok

      call init_event_parameters(2, parameters)

      if(present(blha1_mode)) then
         if(blha1_mode) then
            ! save nlo_prefactors and restore later
            orig_nlo_prefactors=nlo_prefactors
            nlo_prefactors=0
         end if
     end if

      vecs(:,1) = real(momenta(1::5),ki)
      vecs(:,2) = real(momenta(2::5),ki)
      vecs(:,3) = real(momenta(3::5),ki)
      vecs(:,4) = real(momenta(4::5),ki)

      call boost_to_cms(vecs)

      call samplitude(vecs, real(mu,ki)*real(mu,ki), amp, prec, ichecked, ok)
      call ninja_exit()
      call quadninja_exit()
      if (ok) then
         !
      else
         !
      end if
      if(present(acc)) then
         acc=10.0_ki**(-prec) ! point accuracy
      else
         if(prec.lt.PSP_chk_li3 .and. PSP_check) then
            ! Give back a Nan so that point is discarded
            zero = log(1.0_ki)
            amp(2)= 1.0_ki/zero
        end if
        ! Cannot be assigned if present(acc)=F --> commented out!
        ! acc=1E5_ki ! dummy accuracy which is not used
      end if


      res(1) = real(amp(4), c_double)
      res(2) = real(amp(3), c_double)
      res(3) = real(amp(2), c_double)
      res(4) = real(amp(1), c_double)

      if(present(blha1_mode)) then
         if(blha1_mode) then
            ! restore nlo_prefactors
            nlo_prefactors = orig_nlo_prefactors
         end if
     end if

   end subroutine eval14
   !---#] subroutine eval14 :
   !---#[ subroutine eval15 :
   subroutine     eval15(momenta, mu, parameters, res, ichecked, acc, blha1_mode)
      use, intrinsic :: iso_c_binding
      use p10_part1part21_part25part25part1_config, only: ki, PSP_chk_li3, nlo_prefactors, PSP_check
      use p10_part1part21_part25part25part1_model, only: parseline
      use p10_part1part21_part25part25part1_kinematics, only: boost_to_cms
      use p14_part1partm1_part25part25part21_matrix, only: samplitude, OLP_spin_correlated_lo2, OLP_color_correlated
      use p10_part1part21_part25part25part1_groups, only: ninja_exit, quadninja_exit

      implicit none
      real(kind=c_double), dimension(25), intent(in) :: momenta
      real(kind=c_double), intent(in) :: mu
      real(kind=c_double), dimension(10), intent(in) :: parameters
      real(kind=c_double), dimension(60), intent(out) :: res
      integer, intent(in) :: ichecked

      real(kind=ki), dimension(5,4) :: vecs
      real(kind=ki), dimension(4) :: amp
      real(kind=c_double), optional :: acc
      logical, optional :: blha1_mode
      real(kind=ki) :: zero
      integer :: i, prec, orig_nlo_prefactors
      logical :: ok

      call init_event_parameters(2, parameters)

      if(present(blha1_mode)) then
         if(blha1_mode) then
            ! save nlo_prefactors and restore later
            orig_nlo_prefactors=nlo_prefactors
            nlo_prefactors=0
         end if
     end if

      vecs(:,1) = real(momenta(1::5),ki)
      vecs(:,2) = real(momenta(2::5),ki)
      vecs(:,3) = real(momenta(3::5),ki)
      vecs(:,4) = real(momenta(4::5),ki)

      call boost_to_cms(vecs)

      call samplitude(vecs, real(mu,ki)*real(mu,ki), amp, prec, ichecked, ok)
      call ninja_exit()
      call quadninja_exit()
      if (ok) then
         !
      else
         !
      end if
      if(present(acc)) then
         acc=10.0_ki**(-prec) ! point accuracy
      else
         if(prec.lt.PSP_chk_li3 .and. PSP_check) then
            ! Give back a Nan so that point is discarded
            zero = log(1.0_ki)
            amp(2)= 1.0_ki/zero
        end if
        ! Cannot be assigned if present(acc)=F --> commented out!
        ! acc=1E5_ki ! dummy accuracy which is not used
      end if


      res(1) = real(amp(4), c_double)
      res(2) = real(amp(3), c_double)
      res(3) = real(amp(2), c_double)
      res(4) = real(amp(1), c_double)

      if(present(blha1_mode)) then
         if(blha1_mode) then
            ! restore nlo_prefactors
            nlo_prefactors = orig_nlo_prefactors
         end if
     end if

   end subroutine eval15
   !---#] subroutine eval15 :

   !---#[ OLP Polarization vector:
   subroutine OLP_Polvec(p,q,eps) &
       & bind(C,name="olp_polvec_")
      use, intrinsic :: iso_c_binding
      use p0_part21part21_part25part25_config , only:ki
      use p0_part21part21_part25part25_model
      use p0_part21part21_part25part25_kinematics, only: Spab3, Spaa  
      implicit none
      real(kind=c_double), dimension(0:3), intent(in) :: p,q
      real(kind=c_double), dimension(0:7), intent(out) :: eps
      complex(kind=ki), dimension(4) :: eps_complex
      complex(kind=ki), dimension(0:3) :: Sp

      Sp=Spab3(real(q,ki), real(p,ki))

      eps_complex(:)=Sp(:)/Spaa(real(q,ki),real(p,ki))/sqrt2
      eps(0)=real(eps_complex(1),c_double)
      eps(1)=real(aimag(eps_complex(1)),c_double)
      eps(2)=real(eps_complex(2),c_double)
      eps(3)=real(aimag(eps_complex(2)),c_double)
      eps(4)=real(eps_complex(3),c_double)
      eps(5)=real(aimag(eps_complex(3)),c_double)
      eps(6)=real(eps_complex(4),c_double)
      eps(7)=real(aimag(eps_complex(4)),c_double)

   end subroutine OLP_Polvec
   !---#] OLP Polarization vector:

   !---#[ OLP_PrintParameter
   subroutine OLP_PrintParameter(filename) &
       & bind(C,name="olp_printparameter_")

      use, intrinsic :: iso_c_binding
      use p0_part21part21_part25part25_model, only: p0_part21part21_part25part25_print_parameter => print_parameter
      use pb_part21part21_part25part25_model, only: pb_part21part21_part25part25_print_parameter => print_parameter
      use p2_part21part21_part25part25part21_model, only: p2_part21part21_part25part25part21_print_parameter => print_parameter
      implicit none
      character(kind=c_char,len=1), intent(in) :: filename
      integer :: ierr, l
      logical :: exists

      interface
         function strlen(s) bind(C,name='strlen')
            use, intrinsic :: iso_c_binding
            implicit none
            character(kind=c_char,len=1), intent(in) :: s
            integer(kind=c_int) :: strlen
         end function strlen
      end interface

      l = strlen(filename)

      inquire(file=filename(1:l), exist=exists)
      if (exists) then
         open(unit=27,file=filename(1:l),status="old", position="append", action="write",iostat=ierr)
      else
         open(unit=27,file=filename(1:l),status="new",iostat=ierr)
      end if
      if (ierr .ne. 0) then
         write(7,*) "OLP_PrintParameter: Could not open/create:", filename(1:l), "!"
         ierr = -1
      end if
      write (27, "(A)") "####### Setup of SubProcess p0_part21part21_part25part25 #######"
      call p0_part21part21_part25part25_print_parameter(.false.,27)
      write (27, *)
      write (27, "(A)") "####### Setup of SubProcess pb_part21part21_part25part25 #######"
      call pb_part21part21_part25part25_print_parameter(.true.,27)
      write (27, *)
      write (27, "(A)") "####### Setup of SubProcess p2_part21part21_part25part25part21 #######"
      call p2_part21part21_part25part25part21_print_parameter(.true.,27)
      write (27, *)

      close(27)

   end subroutine OLP_PrintParameter
   !---#] OLP_PrintParameter

   subroutine     read_slha_file(line)
      use p0_part21part21_part25part25_model, only: p0_part21part21_part25part25_read_slha => read_slha
      use pb_part21part21_part25part25_model, only: pb_part21part21_part25part25_read_slha => read_slha
      use p2_part21part21_part25part25part21_model, only: p2_part21part21_part25part25part21_read_slha => read_slha
      use p2_part21part21_part25part25part21_model_qp, only: p2_part21part21_part25part25part21_read_slha_qp => read_slha
      use p9_part21part21_part25part25part21_model, only: p9_part21part21_part25part25part21_read_slha => read_slha
      use p9_part21part21_part25part25part21_model_qp, only: p9_part21part21_part25part25part21_read_slha_qp => read_slha
      use p3_part1part21_part25part25part1_model, only: p3_part1part21_part25part25part1_read_slha => read_slha
      use p3_part1part21_part25part25part1_model_qp, only: p3_part1part21_part25part25part1_read_slha_qp => read_slha
      use p10_part1part21_part25part25part1_model, only: p10_part1part21_part25part25part1_read_slha => read_slha
      use p10_part1part21_part25part25part1_model_qp, only: p10_part1part21_part25part25part1_read_slha_qp => read_slha
      implicit none
      character(len=*), intent(in) :: line
      character(len=512) :: file_name
      integer :: ierr

      call unescape_file_name(line, file_name)
      open(unit=27,file=file_name,status='old',iostat=ierr)
      if(ierr.ne.0) then
         print*, "Could not find SLHA model file"
      else
         call p0_part21part21_part25part25_read_slha(27)
         rewind(unit=27)
         call pb_part21part21_part25part25_read_slha(27)
         rewind(unit=27)
         call p2_part21part21_part25part25part21_read_slha(27)
         call p2_part21part21_part25part25part21_read_slha_qp(27)
         rewind(unit=27)
         call p9_part21part21_part25part25part21_read_slha(27)
         call p9_part21part21_part25part25part21_read_slha_qp(27)
         rewind(unit=27)
         call p3_part1part21_part25part25part1_read_slha(27)
         call p3_part1part21_part25part25part1_read_slha_qp(27)
         rewind(unit=27)
         call p10_part1part21_part25part25part1_read_slha(27)
         call p10_part1part21_part25part25part1_read_slha_qp(27)
         close(27)
      end if
   end subroutine read_slha_file

   subroutine     unescape_file_name(source, dest)
      implicit none
      character(len=*), intent(in) :: source
      character(len=512), intent(out) :: dest
      integer :: is, id, l, hex, hexdigit, hexpos
      character(len=512) :: buf
      logical :: special

      is = scan(source, "|")

      if (is > 1) then
         buf = trim(source(1:is-1))
      else
         buf = trim(source)
      end if

      l = len(buf)
      id = 1
      special = .false.
      hexpos = 0
      if (buf(1:1) .eq. '"') then
         ! double quoted string
         do is = 2, l - 1
            if (special) then
               ! after a backslash or in \x.. escape
               if (hexpos == 1 .or. hexpos == 2) then
                  ! interpret hex digit
                  if ("0" .le. buf(is:is) .and. buf(is:is) .le. "9") then
                     hexdigit = ichar(buf(is:is)) - ichar("0")
                  elseif ("A" .le. buf(is:is) .and. buf(is:is) .le. "F") then
                     hexdigit = ichar(buf(is:is)) - ichar("A") + 10
                  elseif ("a" .le. buf(is:is) .and. buf(is:is) .le. "f") then
                     hexdigit = ichar(buf(is:is)) - ichar("a") + 10
                  else
                     print*, "Invalid hex escape sequence in file name"
                     stop
                  end if

                  if (hexpos == 1) then
                     hex = 16 * hexdigit
                     hexpos = 2
                  else
                     hex = hex + hexdigit
                     hexpos = 0
                     special = .false.
                     dest(id:id) = achar(hex)
                     id = id + 1
                  end if
               elseif (buf(is:is) .eq. "n") then
                  dest(id:id) = achar(10)
                  id = id + 1
                  special = .false.
               elseif (buf(is:is) .eq. "r") then
                  dest(id:id) = achar(13)
                  id = id + 1
                  special = .false.
               elseif (buf(is:is) .eq. "f") then
                  dest(id:id) = achar(12)
                  id = id + 1
                  special = .false.
               elseif (buf(is:is) .eq. "t") then
                  dest(id:id) = achar(9)
                  id = id + 1
                  special = .false.
               elseif (buf(is:is) .eq. "x") then
                  hexpos = 1
               else
                  dest(id:id) = buf(is:is)
                  id = id + 1
                  special = .false.
               end if
            else
               if(buf(is:is) .eq. "\") then
                  special = .true.
               else
                  dest(id:id) = source(is:is)
                  id = id + 1
               end if
            end if
         end do
      elseif (buf(1:1) .eq. '"') then
         ! single quoted string
         do is = 2, l - 1
            if (special) then
               dest(id:id) = buf(is:is)
               id = id + 1
               special = .false.
            elseif (buf(is:is) .eq. "'") then
               special = .true.
            else
               dest(id:id) = buf(is:is)
               id = id + 1
            end if
         end do
      else
         ! assume backslash escaped string
         do is = 1, l
            if (special) then
               dest(id:id) = buf(is:is)
               id = id + 1
               special = .false.
            elseif (buf(is:is) .eq. "\") then
               special = .true.
            else
               dest(id:id) = buf(is:is)
               id = id + 1
            end if
         end do
      end if
   end subroutine unescape_file_name
end module olp_module
