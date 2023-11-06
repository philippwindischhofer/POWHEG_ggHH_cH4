module    p2_part21part21_part25part25part21_amplitudeh0
   use p2_part21part21_part25part25part21_config, only: ki, &
       & reduction_interoperation
   use p2_part21part21_part25part25part21_color, only: numcs
   use p2_part21part21_part25part25part21_groups
   use ninjago_module, only: ki_nin
   use p2_part21part21_part25part25part21_ninjah0
   
   implicit none
   private

   public :: finite_renormalisation, samplitude
contains
!---#[ function finite_renormalisation:
   function     finite_renormalisation(scale2) result(amp)
      use p2_part21part21_part25part25part21_util, only: square
      use p2_part21part21_part25part25part21_color, only: CF, CA
      use p2_part21part21_part25part25part21_kinematics, only: &
      & num_light_quarks, num_gluons
      implicit none
      real(ki), intent(in) :: scale2
      real(ki) :: amp
      amp = 0.0_ki
   end function finite_renormalisation
   !---#] function finite_renormalisation:

   !---#[ function samplitude:
   function     samplitude(scale2,ok,rational2,the_col0,opt_perm)
      use p2_part21part21_part25part25part21_config, only: include_eps_terms, include_eps2_terms, &
      & logfile, debug_nlo_diagrams
      use p2_part21part21_part25part25part21_globalsl1, only: col0,perm, use_perm, epspow
      use p2_part21part21_part25part25part21_globalsh0, &
     & only: init_lo, rat2
      use p2_part21part21_part25part25part21_abbrevd38h0, only: init_abbrevd38 => init_abbrev
      use p2_part21part21_part25part25part21_abbrevd56h0, only: init_abbrevd56 => init_abbrev
      use p2_part21part21_part25part25part21_abbrevd125h0, only: init_abbrevd125 => init_abbrev
      use p2_part21part21_part25part25part21_abbrevd149h0, only: init_abbrevd149 => init_abbrev
      use p2_part21part21_part25part25part21_abbrevd206h0, only: init_abbrevd206 => init_abbrev
      use p2_part21part21_part25part25part21_abbrevd196h0, only: init_abbrevd196 => init_abbrev
      use p2_part21part21_part25part25part21_abbrevd121h0, only: init_abbrevd121 => init_abbrev
      use p2_part21part21_part25part25part21_abbrevd204h0, only: init_abbrevd204 => init_abbrev
      use p2_part21part21_part25part25part21_abbrevd200h0, only: init_abbrevd200 => init_abbrev
      use p2_part21part21_part25part25part21_abbrevd192h0, only: init_abbrevd192 => init_abbrev
      use p2_part21part21_part25part25part21_abbrevd188h0, only: init_abbrevd188 => init_abbrev
      use p2_part21part21_part25part25part21_abbrevd137h0, only: init_abbrevd137 => init_abbrev
      use p2_part21part21_part25part25part21_abbrevd202h0, only: init_abbrevd202 => init_abbrev
      use p2_part21part21_part25part25part21_abbrevd133h0, only: init_abbrevd133 => init_abbrev
      use p2_part21part21_part25part25part21_abbrevd198h0, only: init_abbrevd198 => init_abbrev
      use p2_part21part21_part25part25part21_abbrevd52h0, only: init_abbrevd52 => init_abbrev
      use p2_part21part21_part25part25part21_abbrevd129h0, only: init_abbrevd129 => init_abbrev
      use p2_part21part21_part25part25part21_abbrevd145h0, only: init_abbrevd145 => init_abbrev
      use p2_part21part21_part25part25part21_abbrevd194h0, only: init_abbrevd194 => init_abbrev
      use p2_part21part21_part25part25part21_abbrevd117h0, only: init_abbrevd117 => init_abbrev
      use p2_part21part21_part25part25part21_abbrevd184h0, only: init_abbrevd184 => init_abbrev
      use p2_part21part21_part25part25part21_abbrevd34h0, only: init_abbrevd34 => init_abbrev
      use p2_part21part21_part25part25part21_abbrevd42h0, only: init_abbrevd42 => init_abbrev
      use p2_part21part21_part25part25part21_abbrevd54h0, only: init_abbrevd54 => init_abbrev
      use p2_part21part21_part25part25part21_abbrevd109h0, only: init_abbrevd109 => init_abbrev
      use p2_part21part21_part25part25part21_abbrevd141h0, only: init_abbrevd141 => init_abbrev
      use p2_part21part21_part25part25part21_abbrevd147h0, only: init_abbrevd147 => init_abbrev
      use p2_part21part21_part25part25part21_abbrevd190h0, only: init_abbrevd190 => init_abbrev
      use p2_part21part21_part25part25part21_abbrevd113h0, only: init_abbrevd113 => init_abbrev
      use p2_part21part21_part25part25part21_abbrevd186h0, only: init_abbrevd186 => init_abbrev
      use p2_part21part21_part25part25part21_groups
      implicit none
      real(ki), intent(in) :: scale2
      logical, intent(out) :: ok
      real(ki), intent(out) :: rational2
      integer, intent(in) :: the_col0
      integer, dimension(numcs), intent(in), optional :: opt_perm
      complex(ki), dimension(-2:0) :: samplitude

      complex(ki), dimension(-2:0) :: acc
      complex(ki), dimension(0:2,-2:0) :: samp_part

      logical :: acc_ok

      ok = .true.
      rational2 = 0.0_ki

      samplitude(:) = 0.0_ki
      col0 = the_col0
      if (present(opt_perm)) then
         use_perm = .true.
         perm = opt_perm
      else
         use_perm = .false.
      end if

      rat2 = (0.0_ki, 0.0_ki)
      call init_lo()
        call init_abbrevd38()
        call init_abbrevd56()
        call init_abbrevd125()
        call init_abbrevd149()
        call init_abbrevd206()
        call init_abbrevd196()
        call init_abbrevd121()
        call init_abbrevd204()
        call init_abbrevd200()
        call init_abbrevd192()
        call init_abbrevd188()
        call init_abbrevd137()
        call init_abbrevd202()
        call init_abbrevd133()
        call init_abbrevd198()
        call init_abbrevd52()
        call init_abbrevd129()
        call init_abbrevd145()
        call init_abbrevd194()
        call init_abbrevd117()
        call init_abbrevd184()
        call init_abbrevd34()
        call init_abbrevd42()
        call init_abbrevd54()
        call init_abbrevd109()
        call init_abbrevd141()
        call init_abbrevd147()
        call init_abbrevd190()
        call init_abbrevd113()
        call init_abbrevd186()
      epspow=0
      samplitude(-2) = 0.0_ki
      samplitude(-1) = 0.0_ki
      if(debug_nlo_diagrams) then
         write(logfile,'(A22,G24.16,A6,G24.16,A4)') &
         & "<result name='r2' re='", real(rat2, ki), &
         &                 "' im='", aimag(rat2), "' />"
      end if
      samplitude(0) = rat2
         call evaluate_group0(scale2, acc, acc_ok)
         ok = ok .and. acc_ok
         samplitude(:) = samplitude(:) + acc
         call evaluate_group1(scale2, acc, acc_ok)
         ok = ok .and. acc_ok
         samplitude(:) = samplitude(:) + acc
         call evaluate_group2(scale2, acc, acc_ok)
         ok = ok .and. acc_ok
         samplitude(:) = samplitude(:) + acc
         call evaluate_group3(scale2, acc, acc_ok)
         ok = ok .and. acc_ok
         samplitude(:) = samplitude(:) + acc
         call evaluate_group4(scale2, acc, acc_ok)
         ok = ok .and. acc_ok
         samplitude(:) = samplitude(:) + acc
         call evaluate_group5(scale2, acc, acc_ok)
         ok = ok .and. acc_ok
         samplitude(:) = samplitude(:) + acc
         call evaluate_group6(scale2, acc, acc_ok)
         ok = ok .and. acc_ok
         samplitude(:) = samplitude(:) + acc
         call evaluate_group7(scale2, acc, acc_ok)
         ok = ok .and. acc_ok
         samplitude(:) = samplitude(:) + acc
         call evaluate_group8(scale2, acc, acc_ok)
         ok = ok .and. acc_ok
         samplitude(:) = samplitude(:) + acc
         call evaluate_group9(scale2, acc, acc_ok)
         ok = ok .and. acc_ok
         samplitude(:) = samplitude(:) + acc
         call evaluate_group10(scale2, acc, acc_ok)
         ok = ok .and. acc_ok
         samplitude(:) = samplitude(:) + acc
         call evaluate_group11(scale2, acc, acc_ok)
         ok = ok .and. acc_ok
         samplitude(:) = samplitude(:) + acc
   end function samplitude
   !---#] function samplitude:
!---#[ subroutine evaluate_group0:
subroutine     evaluate_group0(scale2,samplitude,ok)
   use p2_part21part21_part25part25part21_config, only: &
      & logfile, debug_nlo_diagrams
   use p2_part21part21_part25part25part21_globalsl1, only: epspow
   use p2_part21part21_part25part25part21_ninjah0, only: ninja_reduce => ninja_reduce_group0
   implicit none
   real(ki), intent(in) :: scale2
   logical, intent(out) :: ok
   complex(ki), dimension(-2:0), intent(out) :: samplitude
   complex(ki_nin), dimension(-2:0) :: tot
   complex(ki_nin) :: totr

   if(debug_nlo_diagrams) then
      write(logfile,*) "<diagram-group index='0'>"
      write(logfile,*) "<param name='epspow' value='", epspow, "'/>"
   end if
   select case(reduction_interoperation)
   case(2) ! use Ninja only
      call ninja_reduce(real(scale2, ki_nin), tot, totr, ok)
      samplitude(:) = cmplx(real(tot(:), ki_nin), aimag(tot(:)), ki)
   case default
      print*, "Your current choice of reduction_interoperation is", &
            & reduction_interoperation
      print*, "This choice is not valid for your current setup."
      print*, "* This code was generated without support for Samurai."
      print*, "* This code was generated with support for Ninja."
      print*, "* This code was generated without support for Golem95."
   end select

   if(debug_nlo_diagrams) then
      write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
         & "<result kind='nlo-finite' re='", real(samplitude(0), ki), &
         & "' im='", aimag(samplitude(0)), "'/>"
      write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
         & "<result kind='nlo-single' re='", real(samplitude(-1), ki), &
         & "' im='", aimag(samplitude(-1)), "'/>"
      write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
         & "<result kind='nlo-double' re='", real(samplitude(-2), ki), &
         & "' im='", aimag(samplitude(-2)), "'/>"
      if(ok) then
         write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
      else
         write(logfile,'(A29)') "<flag name='ok' status='no'/>"
      end if
      write(logfile,*) "</diagram-group>"
   end if
end subroutine evaluate_group0
!---#] subroutine evaluate_group0:
!---#[ subroutine evaluate_group1:
subroutine     evaluate_group1(scale2,samplitude,ok)
   use p2_part21part21_part25part25part21_config, only: &
      & logfile, debug_nlo_diagrams
   use p2_part21part21_part25part25part21_globalsl1, only: epspow
   use p2_part21part21_part25part25part21_ninjah0, only: ninja_reduce => ninja_reduce_group1
   implicit none
   real(ki), intent(in) :: scale2
   logical, intent(out) :: ok
   complex(ki), dimension(-2:0), intent(out) :: samplitude
   complex(ki_nin), dimension(-2:0) :: tot
   complex(ki_nin) :: totr

   if(debug_nlo_diagrams) then
      write(logfile,*) "<diagram-group index='1'>"
      write(logfile,*) "<param name='epspow' value='", epspow, "'/>"
   end if
   select case(reduction_interoperation)
   case(2) ! use Ninja only
      call ninja_reduce(real(scale2, ki_nin), tot, totr, ok)
      samplitude(:) = cmplx(real(tot(:), ki_nin), aimag(tot(:)), ki)
   case default
      print*, "Your current choice of reduction_interoperation is", &
            & reduction_interoperation
      print*, "This choice is not valid for your current setup."
      print*, "* This code was generated without support for Samurai."
      print*, "* This code was generated with support for Ninja."
      print*, "* This code was generated without support for Golem95."
   end select

   if(debug_nlo_diagrams) then
      write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
         & "<result kind='nlo-finite' re='", real(samplitude(0), ki), &
         & "' im='", aimag(samplitude(0)), "'/>"
      write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
         & "<result kind='nlo-single' re='", real(samplitude(-1), ki), &
         & "' im='", aimag(samplitude(-1)), "'/>"
      write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
         & "<result kind='nlo-double' re='", real(samplitude(-2), ki), &
         & "' im='", aimag(samplitude(-2)), "'/>"
      if(ok) then
         write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
      else
         write(logfile,'(A29)') "<flag name='ok' status='no'/>"
      end if
      write(logfile,*) "</diagram-group>"
   end if
end subroutine evaluate_group1
!---#] subroutine evaluate_group1:
!---#[ subroutine evaluate_group2:
subroutine     evaluate_group2(scale2,samplitude,ok)
   use p2_part21part21_part25part25part21_config, only: &
      & logfile, debug_nlo_diagrams
   use p2_part21part21_part25part25part21_globalsl1, only: epspow
   use p2_part21part21_part25part25part21_ninjah0, only: ninja_reduce => ninja_reduce_group2
   implicit none
   real(ki), intent(in) :: scale2
   logical, intent(out) :: ok
   complex(ki), dimension(-2:0), intent(out) :: samplitude
   complex(ki_nin), dimension(-2:0) :: tot
   complex(ki_nin) :: totr

   if(debug_nlo_diagrams) then
      write(logfile,*) "<diagram-group index='2'>"
      write(logfile,*) "<param name='epspow' value='", epspow, "'/>"
   end if
   select case(reduction_interoperation)
   case(2) ! use Ninja only
      call ninja_reduce(real(scale2, ki_nin), tot, totr, ok)
      samplitude(:) = cmplx(real(tot(:), ki_nin), aimag(tot(:)), ki)
   case default
      print*, "Your current choice of reduction_interoperation is", &
            & reduction_interoperation
      print*, "This choice is not valid for your current setup."
      print*, "* This code was generated without support for Samurai."
      print*, "* This code was generated with support for Ninja."
      print*, "* This code was generated without support for Golem95."
   end select

   if(debug_nlo_diagrams) then
      write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
         & "<result kind='nlo-finite' re='", real(samplitude(0), ki), &
         & "' im='", aimag(samplitude(0)), "'/>"
      write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
         & "<result kind='nlo-single' re='", real(samplitude(-1), ki), &
         & "' im='", aimag(samplitude(-1)), "'/>"
      write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
         & "<result kind='nlo-double' re='", real(samplitude(-2), ki), &
         & "' im='", aimag(samplitude(-2)), "'/>"
      if(ok) then
         write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
      else
         write(logfile,'(A29)') "<flag name='ok' status='no'/>"
      end if
      write(logfile,*) "</diagram-group>"
   end if
end subroutine evaluate_group2
!---#] subroutine evaluate_group2:
!---#[ subroutine evaluate_group3:
subroutine     evaluate_group3(scale2,samplitude,ok)
   use p2_part21part21_part25part25part21_config, only: &
      & logfile, debug_nlo_diagrams
   use p2_part21part21_part25part25part21_globalsl1, only: epspow
   use p2_part21part21_part25part25part21_ninjah0, only: ninja_reduce => ninja_reduce_group3
   implicit none
   real(ki), intent(in) :: scale2
   logical, intent(out) :: ok
   complex(ki), dimension(-2:0), intent(out) :: samplitude
   complex(ki_nin), dimension(-2:0) :: tot
   complex(ki_nin) :: totr

   if(debug_nlo_diagrams) then
      write(logfile,*) "<diagram-group index='3'>"
      write(logfile,*) "<param name='epspow' value='", epspow, "'/>"
   end if
   select case(reduction_interoperation)
   case(2) ! use Ninja only
      call ninja_reduce(real(scale2, ki_nin), tot, totr, ok)
      samplitude(:) = cmplx(real(tot(:), ki_nin), aimag(tot(:)), ki)
   case default
      print*, "Your current choice of reduction_interoperation is", &
            & reduction_interoperation
      print*, "This choice is not valid for your current setup."
      print*, "* This code was generated without support for Samurai."
      print*, "* This code was generated with support for Ninja."
      print*, "* This code was generated without support for Golem95."
   end select

   if(debug_nlo_diagrams) then
      write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
         & "<result kind='nlo-finite' re='", real(samplitude(0), ki), &
         & "' im='", aimag(samplitude(0)), "'/>"
      write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
         & "<result kind='nlo-single' re='", real(samplitude(-1), ki), &
         & "' im='", aimag(samplitude(-1)), "'/>"
      write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
         & "<result kind='nlo-double' re='", real(samplitude(-2), ki), &
         & "' im='", aimag(samplitude(-2)), "'/>"
      if(ok) then
         write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
      else
         write(logfile,'(A29)') "<flag name='ok' status='no'/>"
      end if
      write(logfile,*) "</diagram-group>"
   end if
end subroutine evaluate_group3
!---#] subroutine evaluate_group3:
!---#[ subroutine evaluate_group4:
subroutine     evaluate_group4(scale2,samplitude,ok)
   use p2_part21part21_part25part25part21_config, only: &
      & logfile, debug_nlo_diagrams
   use p2_part21part21_part25part25part21_globalsl1, only: epspow
   use p2_part21part21_part25part25part21_ninjah0, only: ninja_reduce => ninja_reduce_group4
   implicit none
   real(ki), intent(in) :: scale2
   logical, intent(out) :: ok
   complex(ki), dimension(-2:0), intent(out) :: samplitude
   complex(ki_nin), dimension(-2:0) :: tot
   complex(ki_nin) :: totr

   if(debug_nlo_diagrams) then
      write(logfile,*) "<diagram-group index='4'>"
      write(logfile,*) "<param name='epspow' value='", epspow, "'/>"
   end if
   select case(reduction_interoperation)
   case(2) ! use Ninja only
      call ninja_reduce(real(scale2, ki_nin), tot, totr, ok)
      samplitude(:) = cmplx(real(tot(:), ki_nin), aimag(tot(:)), ki)
   case default
      print*, "Your current choice of reduction_interoperation is", &
            & reduction_interoperation
      print*, "This choice is not valid for your current setup."
      print*, "* This code was generated without support for Samurai."
      print*, "* This code was generated with support for Ninja."
      print*, "* This code was generated without support for Golem95."
   end select

   if(debug_nlo_diagrams) then
      write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
         & "<result kind='nlo-finite' re='", real(samplitude(0), ki), &
         & "' im='", aimag(samplitude(0)), "'/>"
      write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
         & "<result kind='nlo-single' re='", real(samplitude(-1), ki), &
         & "' im='", aimag(samplitude(-1)), "'/>"
      write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
         & "<result kind='nlo-double' re='", real(samplitude(-2), ki), &
         & "' im='", aimag(samplitude(-2)), "'/>"
      if(ok) then
         write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
      else
         write(logfile,'(A29)') "<flag name='ok' status='no'/>"
      end if
      write(logfile,*) "</diagram-group>"
   end if
end subroutine evaluate_group4
!---#] subroutine evaluate_group4:
!---#[ subroutine evaluate_group5:
subroutine     evaluate_group5(scale2,samplitude,ok)
   use p2_part21part21_part25part25part21_config, only: &
      & logfile, debug_nlo_diagrams
   use p2_part21part21_part25part25part21_globalsl1, only: epspow
   use p2_part21part21_part25part25part21_ninjah0, only: ninja_reduce => ninja_reduce_group5
   implicit none
   real(ki), intent(in) :: scale2
   logical, intent(out) :: ok
   complex(ki), dimension(-2:0), intent(out) :: samplitude
   complex(ki_nin), dimension(-2:0) :: tot
   complex(ki_nin) :: totr

   if(debug_nlo_diagrams) then
      write(logfile,*) "<diagram-group index='5'>"
      write(logfile,*) "<param name='epspow' value='", epspow, "'/>"
   end if
   select case(reduction_interoperation)
   case(2) ! use Ninja only
      call ninja_reduce(real(scale2, ki_nin), tot, totr, ok)
      samplitude(:) = cmplx(real(tot(:), ki_nin), aimag(tot(:)), ki)
   case default
      print*, "Your current choice of reduction_interoperation is", &
            & reduction_interoperation
      print*, "This choice is not valid for your current setup."
      print*, "* This code was generated without support for Samurai."
      print*, "* This code was generated with support for Ninja."
      print*, "* This code was generated without support for Golem95."
   end select

   if(debug_nlo_diagrams) then
      write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
         & "<result kind='nlo-finite' re='", real(samplitude(0), ki), &
         & "' im='", aimag(samplitude(0)), "'/>"
      write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
         & "<result kind='nlo-single' re='", real(samplitude(-1), ki), &
         & "' im='", aimag(samplitude(-1)), "'/>"
      write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
         & "<result kind='nlo-double' re='", real(samplitude(-2), ki), &
         & "' im='", aimag(samplitude(-2)), "'/>"
      if(ok) then
         write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
      else
         write(logfile,'(A29)') "<flag name='ok' status='no'/>"
      end if
      write(logfile,*) "</diagram-group>"
   end if
end subroutine evaluate_group5
!---#] subroutine evaluate_group5:
!---#[ subroutine evaluate_group6:
subroutine     evaluate_group6(scale2,samplitude,ok)
   use p2_part21part21_part25part25part21_config, only: &
      & logfile, debug_nlo_diagrams
   use p2_part21part21_part25part25part21_globalsl1, only: epspow
   use p2_part21part21_part25part25part21_ninjah0, only: ninja_reduce => ninja_reduce_group6
   implicit none
   real(ki), intent(in) :: scale2
   logical, intent(out) :: ok
   complex(ki), dimension(-2:0), intent(out) :: samplitude
   complex(ki_nin), dimension(-2:0) :: tot
   complex(ki_nin) :: totr

   if(debug_nlo_diagrams) then
      write(logfile,*) "<diagram-group index='6'>"
      write(logfile,*) "<param name='epspow' value='", epspow, "'/>"
   end if
   select case(reduction_interoperation)
   case(2) ! use Ninja only
      call ninja_reduce(real(scale2, ki_nin), tot, totr, ok)
      samplitude(:) = cmplx(real(tot(:), ki_nin), aimag(tot(:)), ki)
   case default
      print*, "Your current choice of reduction_interoperation is", &
            & reduction_interoperation
      print*, "This choice is not valid for your current setup."
      print*, "* This code was generated without support for Samurai."
      print*, "* This code was generated with support for Ninja."
      print*, "* This code was generated without support for Golem95."
   end select

   if(debug_nlo_diagrams) then
      write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
         & "<result kind='nlo-finite' re='", real(samplitude(0), ki), &
         & "' im='", aimag(samplitude(0)), "'/>"
      write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
         & "<result kind='nlo-single' re='", real(samplitude(-1), ki), &
         & "' im='", aimag(samplitude(-1)), "'/>"
      write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
         & "<result kind='nlo-double' re='", real(samplitude(-2), ki), &
         & "' im='", aimag(samplitude(-2)), "'/>"
      if(ok) then
         write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
      else
         write(logfile,'(A29)') "<flag name='ok' status='no'/>"
      end if
      write(logfile,*) "</diagram-group>"
   end if
end subroutine evaluate_group6
!---#] subroutine evaluate_group6:
!---#[ subroutine evaluate_group7:
subroutine     evaluate_group7(scale2,samplitude,ok)
   use p2_part21part21_part25part25part21_config, only: &
      & logfile, debug_nlo_diagrams
   use p2_part21part21_part25part25part21_globalsl1, only: epspow
   use p2_part21part21_part25part25part21_ninjah0, only: ninja_reduce => ninja_reduce_group7
   implicit none
   real(ki), intent(in) :: scale2
   logical, intent(out) :: ok
   complex(ki), dimension(-2:0), intent(out) :: samplitude
   complex(ki_nin), dimension(-2:0) :: tot
   complex(ki_nin) :: totr

   if(debug_nlo_diagrams) then
      write(logfile,*) "<diagram-group index='7'>"
      write(logfile,*) "<param name='epspow' value='", epspow, "'/>"
   end if
   select case(reduction_interoperation)
   case(2) ! use Ninja only
      call ninja_reduce(real(scale2, ki_nin), tot, totr, ok)
      samplitude(:) = cmplx(real(tot(:), ki_nin), aimag(tot(:)), ki)
   case default
      print*, "Your current choice of reduction_interoperation is", &
            & reduction_interoperation
      print*, "This choice is not valid for your current setup."
      print*, "* This code was generated without support for Samurai."
      print*, "* This code was generated with support for Ninja."
      print*, "* This code was generated without support for Golem95."
   end select

   if(debug_nlo_diagrams) then
      write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
         & "<result kind='nlo-finite' re='", real(samplitude(0), ki), &
         & "' im='", aimag(samplitude(0)), "'/>"
      write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
         & "<result kind='nlo-single' re='", real(samplitude(-1), ki), &
         & "' im='", aimag(samplitude(-1)), "'/>"
      write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
         & "<result kind='nlo-double' re='", real(samplitude(-2), ki), &
         & "' im='", aimag(samplitude(-2)), "'/>"
      if(ok) then
         write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
      else
         write(logfile,'(A29)') "<flag name='ok' status='no'/>"
      end if
      write(logfile,*) "</diagram-group>"
   end if
end subroutine evaluate_group7
!---#] subroutine evaluate_group7:
!---#[ subroutine evaluate_group8:
subroutine     evaluate_group8(scale2,samplitude,ok)
   use p2_part21part21_part25part25part21_config, only: &
      & logfile, debug_nlo_diagrams
   use p2_part21part21_part25part25part21_globalsl1, only: epspow
   use p2_part21part21_part25part25part21_ninjah0, only: ninja_reduce => ninja_reduce_group8
   implicit none
   real(ki), intent(in) :: scale2
   logical, intent(out) :: ok
   complex(ki), dimension(-2:0), intent(out) :: samplitude
   complex(ki_nin), dimension(-2:0) :: tot
   complex(ki_nin) :: totr

   if(debug_nlo_diagrams) then
      write(logfile,*) "<diagram-group index='8'>"
      write(logfile,*) "<param name='epspow' value='", epspow, "'/>"
   end if
   select case(reduction_interoperation)
   case(2) ! use Ninja only
      call ninja_reduce(real(scale2, ki_nin), tot, totr, ok)
      samplitude(:) = cmplx(real(tot(:), ki_nin), aimag(tot(:)), ki)
   case default
      print*, "Your current choice of reduction_interoperation is", &
            & reduction_interoperation
      print*, "This choice is not valid for your current setup."
      print*, "* This code was generated without support for Samurai."
      print*, "* This code was generated with support for Ninja."
      print*, "* This code was generated without support for Golem95."
   end select

   if(debug_nlo_diagrams) then
      write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
         & "<result kind='nlo-finite' re='", real(samplitude(0), ki), &
         & "' im='", aimag(samplitude(0)), "'/>"
      write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
         & "<result kind='nlo-single' re='", real(samplitude(-1), ki), &
         & "' im='", aimag(samplitude(-1)), "'/>"
      write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
         & "<result kind='nlo-double' re='", real(samplitude(-2), ki), &
         & "' im='", aimag(samplitude(-2)), "'/>"
      if(ok) then
         write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
      else
         write(logfile,'(A29)') "<flag name='ok' status='no'/>"
      end if
      write(logfile,*) "</diagram-group>"
   end if
end subroutine evaluate_group8
!---#] subroutine evaluate_group8:
!---#[ subroutine evaluate_group9:
subroutine     evaluate_group9(scale2,samplitude,ok)
   use p2_part21part21_part25part25part21_config, only: &
      & logfile, debug_nlo_diagrams
   use p2_part21part21_part25part25part21_globalsl1, only: epspow
   use p2_part21part21_part25part25part21_ninjah0, only: ninja_reduce => ninja_reduce_group9
   implicit none
   real(ki), intent(in) :: scale2
   logical, intent(out) :: ok
   complex(ki), dimension(-2:0), intent(out) :: samplitude
   complex(ki_nin), dimension(-2:0) :: tot
   complex(ki_nin) :: totr

   if(debug_nlo_diagrams) then
      write(logfile,*) "<diagram-group index='9'>"
      write(logfile,*) "<param name='epspow' value='", epspow, "'/>"
   end if
   select case(reduction_interoperation)
   case(2) ! use Ninja only
      call ninja_reduce(real(scale2, ki_nin), tot, totr, ok)
      samplitude(:) = cmplx(real(tot(:), ki_nin), aimag(tot(:)), ki)
   case default
      print*, "Your current choice of reduction_interoperation is", &
            & reduction_interoperation
      print*, "This choice is not valid for your current setup."
      print*, "* This code was generated without support for Samurai."
      print*, "* This code was generated with support for Ninja."
      print*, "* This code was generated without support for Golem95."
   end select

   if(debug_nlo_diagrams) then
      write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
         & "<result kind='nlo-finite' re='", real(samplitude(0), ki), &
         & "' im='", aimag(samplitude(0)), "'/>"
      write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
         & "<result kind='nlo-single' re='", real(samplitude(-1), ki), &
         & "' im='", aimag(samplitude(-1)), "'/>"
      write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
         & "<result kind='nlo-double' re='", real(samplitude(-2), ki), &
         & "' im='", aimag(samplitude(-2)), "'/>"
      if(ok) then
         write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
      else
         write(logfile,'(A29)') "<flag name='ok' status='no'/>"
      end if
      write(logfile,*) "</diagram-group>"
   end if
end subroutine evaluate_group9
!---#] subroutine evaluate_group9:
!---#[ subroutine evaluate_group10:
subroutine     evaluate_group10(scale2,samplitude,ok)
   use p2_part21part21_part25part25part21_config, only: &
      & logfile, debug_nlo_diagrams
   use p2_part21part21_part25part25part21_globalsl1, only: epspow
   use p2_part21part21_part25part25part21_ninjah0, only: ninja_reduce => ninja_reduce_group10
   implicit none
   real(ki), intent(in) :: scale2
   logical, intent(out) :: ok
   complex(ki), dimension(-2:0), intent(out) :: samplitude
   complex(ki_nin), dimension(-2:0) :: tot
   complex(ki_nin) :: totr

   if(debug_nlo_diagrams) then
      write(logfile,*) "<diagram-group index='10'>"
      write(logfile,*) "<param name='epspow' value='", epspow, "'/>"
   end if
   select case(reduction_interoperation)
   case(2) ! use Ninja only
      call ninja_reduce(real(scale2, ki_nin), tot, totr, ok)
      samplitude(:) = cmplx(real(tot(:), ki_nin), aimag(tot(:)), ki)
   case default
      print*, "Your current choice of reduction_interoperation is", &
            & reduction_interoperation
      print*, "This choice is not valid for your current setup."
      print*, "* This code was generated without support for Samurai."
      print*, "* This code was generated with support for Ninja."
      print*, "* This code was generated without support for Golem95."
   end select

   if(debug_nlo_diagrams) then
      write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
         & "<result kind='nlo-finite' re='", real(samplitude(0), ki), &
         & "' im='", aimag(samplitude(0)), "'/>"
      write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
         & "<result kind='nlo-single' re='", real(samplitude(-1), ki), &
         & "' im='", aimag(samplitude(-1)), "'/>"
      write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
         & "<result kind='nlo-double' re='", real(samplitude(-2), ki), &
         & "' im='", aimag(samplitude(-2)), "'/>"
      if(ok) then
         write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
      else
         write(logfile,'(A29)') "<flag name='ok' status='no'/>"
      end if
      write(logfile,*) "</diagram-group>"
   end if
end subroutine evaluate_group10
!---#] subroutine evaluate_group10:
!---#[ subroutine evaluate_group11:
subroutine     evaluate_group11(scale2,samplitude,ok)
   use p2_part21part21_part25part25part21_config, only: &
      & logfile, debug_nlo_diagrams
   use p2_part21part21_part25part25part21_globalsl1, only: epspow
   use p2_part21part21_part25part25part21_ninjah0, only: ninja_reduce => ninja_reduce_group11
   implicit none
   real(ki), intent(in) :: scale2
   logical, intent(out) :: ok
   complex(ki), dimension(-2:0), intent(out) :: samplitude
   complex(ki_nin), dimension(-2:0) :: tot
   complex(ki_nin) :: totr

   if(debug_nlo_diagrams) then
      write(logfile,*) "<diagram-group index='11'>"
      write(logfile,*) "<param name='epspow' value='", epspow, "'/>"
   end if
   select case(reduction_interoperation)
   case(2) ! use Ninja only
      call ninja_reduce(real(scale2, ki_nin), tot, totr, ok)
      samplitude(:) = cmplx(real(tot(:), ki_nin), aimag(tot(:)), ki)
   case default
      print*, "Your current choice of reduction_interoperation is", &
            & reduction_interoperation
      print*, "This choice is not valid for your current setup."
      print*, "* This code was generated without support for Samurai."
      print*, "* This code was generated with support for Ninja."
      print*, "* This code was generated without support for Golem95."
   end select

   if(debug_nlo_diagrams) then
      write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
         & "<result kind='nlo-finite' re='", real(samplitude(0), ki), &
         & "' im='", aimag(samplitude(0)), "'/>"
      write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
         & "<result kind='nlo-single' re='", real(samplitude(-1), ki), &
         & "' im='", aimag(samplitude(-1)), "'/>"
      write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
         & "<result kind='nlo-double' re='", real(samplitude(-2), ki), &
         & "' im='", aimag(samplitude(-2)), "'/>"
      if(ok) then
         write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
      else
         write(logfile,'(A29)') "<flag name='ok' status='no'/>"
      end if
      write(logfile,*) "</diagram-group>"
   end if
end subroutine evaluate_group11
!---#] subroutine evaluate_group11:
end module p2_part21part21_part25part25part21_amplitudeh0
