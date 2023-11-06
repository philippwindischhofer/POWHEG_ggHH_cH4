!
module     p10_part1part21_part25part25part1_ninjah2
   ! This file has been generated for ninja
   use ninjago_module, only: ki_nin
   use p10_part1part21_part25part25part1_config
   implicit none
   private
   public :: ninja_reduce_group0
   public :: ninja_reduce_group1
   public :: ninja_reduce_group2
contains
!---#[ reduce groups with ninja:
!-----#[ subroutine ninja_reduce_group0:
subroutine     ninja_reduce_group0(scale2,tot,totr,ok)
   use iso_c_binding, only: c_ptr, c_loc, c_int
   use ninjago_module
   use p10_part1part21_part25part25part1_kinematics
   use p10_part1part21_part25part25part1_model
   use p10_part1part21_part25part25part1_d84h2l1, only: numerator_diagram84 => &
   &numerator_ninja
   use p10_part1part21_part25part25part1_d84h2l121, only: numerator_tmu_diagram8&
   &4 => numerator_tmu
   use p10_part1part21_part25part25part1_d84h2l131, only: numerator_t3_diagram84 &
   &=> numerator_t3
   use p10_part1part21_part25part25part1_d84h2l132, only: numerator_t2_diagram84 &
   &=> numerator_t2
   use p10_part1part21_part25part25part1_globalsl1, only: epspow

   implicit none
   real(ki_nin), intent(in) :: scale2
   complex(ki_nin), dimension(-2:0), intent(out) :: tot
   complex(ki_nin), intent(out) :: totr
   logical, intent(out) :: ok

   complex(ki_nin), dimension(-2:0) :: acc
   complex(ki_nin) :: accr
   integer(c_int) :: acc_ok

   integer :: istopm, istop0

   integer, parameter :: effective_group_rank = 4

   !-----------#[ invariants for ninja:
   real(ki_nin), dimension(4,4) :: s_mat
   !-----------#] initialize invariants:
   real(ki_nin), dimension(4) :: msq
   real(ki_nin), dimension(4,4) :: Vi

   ok = .true.

   if (ninja_test.eq.1) then
      istopm = 1
      istop0 = 1
   else
      istopm = ninja_istop
      istop0 = max(2,ninja_istop)
   end if
   msq(1) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,1) = real(-k2+k4, ki_nin)
   msq(2) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,2) = real(-k2, ki_nin)
   msq(3) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,3) = real(0, ki_nin)
   msq(4) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,4) = real(-k3, ki_nin)
   !-----------#[ initialize invariants:
   s_mat(1,1) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(1,2) = real(-2.0_ki*mdlMT**2+mdlMh**2, ki_nin)
   s_mat(2,1) = s_mat(1,2)
   s_mat(1,3) = real(-2.0_ki*mdlMT**2+2.0_ki*mdlMh**2-es23+es51-es34, ki_nin)
   s_mat(3,1) = s_mat(1,3)
   s_mat(1,4) = real(-2.0_ki*mdlMT**2+es51, ki_nin)
   s_mat(4,1) = s_mat(1,4)
   s_mat(2,2) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(2,3) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(3,2) = s_mat(2,3)
   s_mat(2,4) = real(-2.0_ki*mdlMT**2+es23, ki_nin)
   s_mat(4,2) = s_mat(2,4)
   s_mat(3,3) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(3,4) = real(-2.0_ki*mdlMT**2+mdlMh**2, ki_nin)
   s_mat(4,3) = s_mat(3,4)
   s_mat(4,4) = real(-2.0_ki*mdlMT**2, ki_nin)
   !-----------#] initialize invariants


      !------#[ sum over reduction of single diagrams:
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='84'>"
         end if
         call ninja_diagram(numerator_diagram84, numerator_t3_diagram84, numera&
         &tor_t2_diagram84, &
          &  numerator_tmu_diagram84, &
          & 4, 4, 4, (/1,2,3,4/), &
          & Vi, msq, s_mat, scale2, istopm, &
          & acc, accr, acc_ok)
            ok = ok .and. (acc_ok.eq.NINJA_SUCCESS)
         if(debug_nlo_diagrams) then
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-finite' re='", -real(acc(0), ki), &
               & "' im='", aimag(acc(0)), "'/>"
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-single' re='", -real(acc(-1), ki), &
               & "' im='", aimag(acc(-1)), "'/>"
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-double' re='", -real(acc(-2), ki), &
               & "' im='", aimag(acc(-2)), "'/>"
            write(logfile,'(A32,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-rational' re='", -real(accr, ki), &
               & "' im='", aimag(accr), "'/>"
            if(ok) then
               write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
            else
               write(logfile,'(A29)') "<flag name='ok' status='no'/>"
            end if
            write(logfile,*) "</diagram>"
         end if
         tot =  + acc
         totr =  + accr
      !------#] sum over reduction of single diagrams:
end subroutine ninja_reduce_group0
!-----#] subroutine ninja_reduce_group0:
!-----#[ subroutine ninja_reduce_group1:
subroutine     ninja_reduce_group1(scale2,tot,totr,ok)
   use iso_c_binding, only: c_ptr, c_loc, c_int
   use ninjago_module
   use p10_part1part21_part25part25part1_kinematics
   use p10_part1part21_part25part25part1_model
   use p10_part1part21_part25part25part1_d11h2l1, only: numerator_diagram11 => &
   &numerator_ninja
   use p10_part1part21_part25part25part1_d11h2l121, only: numerator_tmu_diagram1&
   &1 => numerator_tmu
   use p10_part1part21_part25part25part1_d11h2l131, only: numerator_t3_diagram11 &
   &=> numerator_t3
   use p10_part1part21_part25part25part1_d11h2l132, only: numerator_t2_diagram11 &
   &=> numerator_t2
   use p10_part1part21_part25part25part1_d86h2l1, only: numerator_diagram86 => &
   &numerator_ninja
   use p10_part1part21_part25part25part1_d86h2l121, only: numerator_tmu_diagram8&
   &6 => numerator_tmu
   use p10_part1part21_part25part25part1_d86h2l131, only: numerator_t3_diagram86 &
   &=> numerator_t3
   use p10_part1part21_part25part25part1_d86h2l132, only: numerator_t2_diagram86 &
   &=> numerator_t2
   use p10_part1part21_part25part25part1_globalsl1, only: epspow

   implicit none
   real(ki_nin), intent(in) :: scale2
   complex(ki_nin), dimension(-2:0), intent(out) :: tot
   complex(ki_nin), intent(out) :: totr
   logical, intent(out) :: ok

   complex(ki_nin), dimension(-2:0) :: acc
   complex(ki_nin) :: accr
   integer(c_int) :: acc_ok

   integer :: istopm, istop0

   integer, parameter :: effective_group_rank = 4

   !-----------#[ invariants for ninja:
   real(ki_nin), dimension(4,4) :: s_mat
   !-----------#] initialize invariants:
   real(ki_nin), dimension(4) :: msq
   real(ki_nin), dimension(4,4) :: Vi

   ok = .true.

   if (ninja_test.eq.1) then
      istopm = 1
      istop0 = 1
   else
      istopm = ninja_istop
      istop0 = max(2,ninja_istop)
   end if
   msq(1) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,1) = real(-k3-k4, ki_nin)
   msq(2) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,2) = real(-k3, ki_nin)
   msq(3) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,3) = real(0, ki_nin)
   msq(4) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,4) = real(-k2, ki_nin)
   !-----------#[ initialize invariants:
   s_mat(1,1) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(1,2) = real(-2.0_ki*mdlMT**2+mdlMh**2, ki_nin)
   s_mat(2,1) = s_mat(1,2)
   s_mat(1,3) = real(-2.0_ki*mdlMT**2+es34, ki_nin)
   s_mat(3,1) = s_mat(1,3)
   s_mat(1,4) = real(-2.0_ki*mdlMT**2+es51, ki_nin)
   s_mat(4,1) = s_mat(1,4)
   s_mat(2,2) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(2,3) = real(-2.0_ki*mdlMT**2+mdlMh**2, ki_nin)
   s_mat(3,2) = s_mat(2,3)
   s_mat(2,4) = real(-2.0_ki*mdlMT**2+es23, ki_nin)
   s_mat(4,2) = s_mat(2,4)
   s_mat(3,3) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(3,4) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(4,3) = s_mat(3,4)
   s_mat(4,4) = real(-2.0_ki*mdlMT**2, ki_nin)
   !-----------#] initialize invariants


      !------#[ sum over reduction of single diagrams:
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='11'>"
         end if
         call ninja_diagram(numerator_diagram11, numerator_t3_diagram11, numera&
         &tor_t2_diagram11, &
          &  numerator_tmu_diagram11, &
          & 4, 3, 3, (/1,3,4/), &
          & Vi, msq, s_mat, scale2, istopm, &
          & acc, accr, acc_ok)
            ok = ok .and. (acc_ok.eq.NINJA_SUCCESS)
         if(debug_nlo_diagrams) then
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-finite' re='", -real(acc(0), ki), &
               & "' im='", aimag(acc(0)), "'/>"
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-single' re='", -real(acc(-1), ki), &
               & "' im='", aimag(acc(-1)), "'/>"
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-double' re='", -real(acc(-2), ki), &
               & "' im='", aimag(acc(-2)), "'/>"
            write(logfile,'(A32,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-rational' re='", -real(accr, ki), &
               & "' im='", aimag(accr), "'/>"
            if(ok) then
               write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
            else
               write(logfile,'(A29)') "<flag name='ok' status='no'/>"
            end if
            write(logfile,*) "</diagram>"
         end if
         tot =  + acc
         totr =  + accr
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='86'>"
         end if
         call ninja_diagram(numerator_diagram86, numerator_t3_diagram86, numera&
         &tor_t2_diagram86, &
          &  numerator_tmu_diagram86, &
          & 4, 4, 4, (/1,2,3,4/), &
          & Vi, msq, s_mat, scale2, istopm, &
          & acc, accr, acc_ok)
            ok = ok .and. (acc_ok.eq.NINJA_SUCCESS)
         if(debug_nlo_diagrams) then
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-finite' re='", -real(acc(0), ki), &
               & "' im='", aimag(acc(0)), "'/>"
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-single' re='", -real(acc(-1), ki), &
               & "' im='", aimag(acc(-1)), "'/>"
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-double' re='", -real(acc(-2), ki), &
               & "' im='", aimag(acc(-2)), "'/>"
            write(logfile,'(A32,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-rational' re='", -real(accr, ki), &
               & "' im='", aimag(accr), "'/>"
            if(ok) then
               write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
            else
               write(logfile,'(A29)') "<flag name='ok' status='no'/>"
            end if
            write(logfile,*) "</diagram>"
         end if
         tot = tot  + acc
         totr = totr  + accr
      !------#] sum over reduction of single diagrams:
end subroutine ninja_reduce_group1
!-----#] subroutine ninja_reduce_group1:
!-----#[ subroutine ninja_reduce_group2:
subroutine     ninja_reduce_group2(scale2,tot,totr,ok)
   use iso_c_binding, only: c_ptr, c_loc, c_int
   use ninjago_module
   use p10_part1part21_part25part25part1_kinematics
   use p10_part1part21_part25part25part1_model
   use p10_part1part21_part25part25part1_d88h2l1, only: numerator_diagram88 => &
   &numerator_ninja
   use p10_part1part21_part25part25part1_d88h2l121, only: numerator_tmu_diagram8&
   &8 => numerator_tmu
   use p10_part1part21_part25part25part1_d88h2l131, only: numerator_t3_diagram88 &
   &=> numerator_t3
   use p10_part1part21_part25part25part1_d88h2l132, only: numerator_t2_diagram88 &
   &=> numerator_t2
   use p10_part1part21_part25part25part1_globalsl1, only: epspow

   implicit none
   real(ki_nin), intent(in) :: scale2
   complex(ki_nin), dimension(-2:0), intent(out) :: tot
   complex(ki_nin), intent(out) :: totr
   logical, intent(out) :: ok

   complex(ki_nin), dimension(-2:0) :: acc
   complex(ki_nin) :: accr
   integer(c_int) :: acc_ok

   integer :: istopm, istop0

   integer, parameter :: effective_group_rank = 4

   !-----------#[ invariants for ninja:
   real(ki_nin), dimension(4,4) :: s_mat
   !-----------#] initialize invariants:
   real(ki_nin), dimension(4) :: msq
   real(ki_nin), dimension(4,4) :: Vi

   ok = .true.

   if (ninja_test.eq.1) then
      istopm = 1
      istop0 = 1
   else
      istopm = ninja_istop
      istop0 = max(2,ninja_istop)
   end if
   msq(1) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,1) = real(-k3-k4, ki_nin)
   msq(2) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,2) = real(-k4, ki_nin)
   msq(3) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,3) = real(0, ki_nin)
   msq(4) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,4) = real(-k2, ki_nin)
   !-----------#[ initialize invariants:
   s_mat(1,1) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(1,2) = real(-2.0_ki*mdlMT**2+mdlMh**2, ki_nin)
   s_mat(2,1) = s_mat(1,2)
   s_mat(1,3) = real(-2.0_ki*mdlMT**2+es34, ki_nin)
   s_mat(3,1) = s_mat(1,3)
   s_mat(1,4) = real(-2.0_ki*mdlMT**2+es51, ki_nin)
   s_mat(4,1) = s_mat(1,4)
   s_mat(2,2) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(2,3) = real(-2.0_ki*mdlMT**2+mdlMh**2, ki_nin)
   s_mat(3,2) = s_mat(2,3)
   s_mat(2,4) = real(-2.0_ki*mdlMT**2+2.0_ki*mdlMh**2-es23+es51-es34, ki_nin)
   s_mat(4,2) = s_mat(2,4)
   s_mat(3,3) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(3,4) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(4,3) = s_mat(3,4)
   s_mat(4,4) = real(-2.0_ki*mdlMT**2, ki_nin)
   !-----------#] initialize invariants


      !------#[ sum over reduction of single diagrams:
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='88'>"
         end if
         call ninja_diagram(numerator_diagram88, numerator_t3_diagram88, numera&
         &tor_t2_diagram88, &
          &  numerator_tmu_diagram88, &
          & 4, 4, 4, (/1,2,3,4/), &
          & Vi, msq, s_mat, scale2, istopm, &
          & acc, accr, acc_ok)
            ok = ok .and. (acc_ok.eq.NINJA_SUCCESS)
         if(debug_nlo_diagrams) then
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-finite' re='", -real(acc(0), ki), &
               & "' im='", aimag(acc(0)), "'/>"
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-single' re='", -real(acc(-1), ki), &
               & "' im='", aimag(acc(-1)), "'/>"
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-double' re='", -real(acc(-2), ki), &
               & "' im='", aimag(acc(-2)), "'/>"
            write(logfile,'(A32,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-rational' re='", -real(accr, ki), &
               & "' im='", aimag(accr), "'/>"
            if(ok) then
               write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
            else
               write(logfile,'(A29)') "<flag name='ok' status='no'/>"
            end if
            write(logfile,*) "</diagram>"
         end if
         tot =  + acc
         totr =  + accr
      !------#] sum over reduction of single diagrams:
end subroutine ninja_reduce_group2
!-----#] subroutine ninja_reduce_group2:
!---#] reduce groups with ninja:
end module p10_part1part21_part25part25part1_ninjah2
