!
module     p9_part21part21_part25part25part21_ninjah0
   ! This file has been generated for ninja
   use ninjago_module, only: ki_nin
   use p9_part21part21_part25part25part21_config
   implicit none
   private
   public :: ninja_reduce_group0
   public :: ninja_reduce_group1
   public :: ninja_reduce_group2
   public :: ninja_reduce_group3
   public :: ninja_reduce_group4
   public :: ninja_reduce_group5
   public :: ninja_reduce_group6
   public :: ninja_reduce_group7
   public :: ninja_reduce_group8
   public :: ninja_reduce_group9
   public :: ninja_reduce_group10
   public :: ninja_reduce_group11
contains
!---#[ reduce groups with ninja:
!-----#[ subroutine ninja_reduce_group0:
subroutine     ninja_reduce_group0(scale2,tot,totr,ok)
   use iso_c_binding, only: c_ptr, c_loc, c_int
   use ninjago_module
   use p9_part21part21_part25part25part21_kinematics
   use p9_part21part21_part25part25part21_model
   use p9_part21part21_part25part25part21_d38h0l1, only: numerator_diagram38 =>&
   & numerator_ninja
   use p9_part21part21_part25part25part21_d38h0l121, only: numerator_tmu_diagram3&
   &8 => numerator_tmu
   use p9_part21part21_part25part25part21_d38h0l131, only: numerator_t3_diagram3&
   &8 => numerator_t3
   use p9_part21part21_part25part25part21_d38h0l132, only: numerator_t2_diagram3&
   &8 => numerator_t2
   use p9_part21part21_part25part25part21_d56h0l1, only: numerator_diagram56 =>&
   & numerator_ninja
   use p9_part21part21_part25part25part21_d56h0l121, only: numerator_tmu_diagram5&
   &6 => numerator_tmu
   use p9_part21part21_part25part25part21_d56h0l131, only: numerator_t3_diagram5&
   &6 => numerator_t3
   use p9_part21part21_part25part25part21_d56h0l132, only: numerator_t2_diagram5&
   &6 => numerator_t2
   use p9_part21part21_part25part25part21_d125h0l1, only: numerator_diagram125 &
   &=> numerator_ninja
   use p9_part21part21_part25part25part21_d125h0l121, only: numerator_tmu_diagr&
   &am125 => numerator_tmu
   use p9_part21part21_part25part25part21_d125h0l131, only: numerator_t3_diagram1&
   &25 => numerator_t3
   use p9_part21part21_part25part25part21_d125h0l132, only: numerator_t2_diagram1&
   &25 => numerator_t2
   use p9_part21part21_part25part25part21_d149h0l1, only: numerator_diagram149 &
   &=> numerator_ninja
   use p9_part21part21_part25part25part21_d149h0l121, only: numerator_tmu_diagr&
   &am149 => numerator_tmu
   use p9_part21part21_part25part25part21_d149h0l131, only: numerator_t3_diagram1&
   &49 => numerator_t3
   use p9_part21part21_part25part25part21_d149h0l132, only: numerator_t2_diagram1&
   &49 => numerator_t2
   use p9_part21part21_part25part25part21_d206h0l1, only: numerator_diagram206 &
   &=> numerator_ninja
   use p9_part21part21_part25part25part21_d206h0l121, only: numerator_tmu_diagr&
   &am206 => numerator_tmu
   use p9_part21part21_part25part25part21_d206h0l131, only: numerator_t3_diagram2&
   &06 => numerator_t3
   use p9_part21part21_part25part25part21_d206h0l132, only: numerator_t2_diagram2&
   &06 => numerator_t2
   use p9_part21part21_part25part25part21_globalsl1, only: epspow

   implicit none
   real(ki_nin), intent(in) :: scale2
   complex(ki_nin), dimension(-2:0), intent(out) :: tot
   complex(ki_nin), intent(out) :: totr
   logical, intent(out) :: ok

   complex(ki_nin), dimension(-2:0) :: acc
   complex(ki_nin) :: accr
   integer(c_int) :: acc_ok

   integer :: istopm, istop0

   integer, parameter :: effective_group_rank = 5

   !-----------#[ invariants for ninja:
   real(ki_nin), dimension(5,5) :: s_mat
   !-----------#] initialize invariants:
   real(ki_nin), dimension(5) :: msq
   real(ki_nin), dimension(4,5) :: Vi

   ok = .true.

   if (ninja_test.eq.1) then
      istopm = 1
      istop0 = 1
   else
      istopm = ninja_istop
      istop0 = max(2,ninja_istop)
   end if
   msq(1) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,1) = real(-k2+k5, ki_nin)
   msq(2) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,2) = real(-k2, ki_nin)
   msq(3) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,3) = real(0, ki_nin)
   msq(4) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,4) = real(-k3, ki_nin)
   msq(5) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,5) = real(-k3-k4, ki_nin)
   !-----------#[ initialize invariants:
   s_mat(1,1) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(1,2) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(2,1) = s_mat(1,2)
   s_mat(1,3) = real(-2.0_ki*mdlMT**2-es12+es34-es51, ki_nin)
   s_mat(3,1) = s_mat(1,3)
   s_mat(1,4) = real(-2.0_ki*mdlMT**2+mdlMh**2+es23-es45-es51, ki_nin)
   s_mat(4,1) = s_mat(1,4)
   s_mat(1,5) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(5,1) = s_mat(1,5)
   s_mat(2,2) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(2,3) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(3,2) = s_mat(2,3)
   s_mat(2,4) = real(-2.0_ki*mdlMT**2+es23, ki_nin)
   s_mat(4,2) = s_mat(2,4)
   s_mat(2,5) = real(-2.0_ki*mdlMT**2+es51, ki_nin)
   s_mat(5,2) = s_mat(2,5)
   s_mat(3,3) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(3,4) = real(-2.0_ki*mdlMT**2+mdlMh**2, ki_nin)
   s_mat(4,3) = s_mat(3,4)
   s_mat(3,5) = real(-2.0_ki*mdlMT**2+es34, ki_nin)
   s_mat(5,3) = s_mat(3,5)
   s_mat(4,4) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(4,5) = real(-2.0_ki*mdlMT**2+mdlMh**2, ki_nin)
   s_mat(5,4) = s_mat(4,5)
   s_mat(5,5) = real(-2.0_ki*mdlMT**2, ki_nin)
   !-----------#] initialize invariants


      !------#[ sum over reduction of single diagrams:
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='38'>"
         end if
         call ninja_diagram(numerator_diagram38, numerator_t3_diagram38, numera&
         &tor_t2_diagram38, &
          &  numerator_tmu_diagram38, &
          & 5, 3, 3, (/2,3,5/), &
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
            write(logfile,*) "<diagram index='56'>"
         end if
         call ninja_diagram(numerator_diagram56, numerator_t3_diagram56, numera&
         &tor_t2_diagram56, &
          &  numerator_tmu_diagram56, &
          & 5, 4, 4, (/1,2,3,5/), &
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
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='125'>"
         end if
         call ninja_diagram(numerator_diagram125, numerator_t3_diagram125, nume&
         &rator_t2_diagram125, &
          &  numerator_tmu_diagram125, &
          & 5, 4, 4, (/2,3,4,5/), &
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
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='149'>"
         end if
         call ninja_diagram(numerator_diagram149, numerator_t3_diagram149, nume&
         &rator_t2_diagram149, &
          &  numerator_tmu_diagram149, &
          & 5, 4, 4, (/1,2,3,5/), &
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
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='206'>"
         end if
         call ninja_diagram(numerator_diagram206, numerator_t3_diagram206, nume&
         &rator_t2_diagram206, &
          &  numerator_tmu_diagram206, &
          & 5, 5, 5, (/1,2,3,4,5/), &
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
end subroutine ninja_reduce_group0
!-----#] subroutine ninja_reduce_group0:
!-----#[ subroutine ninja_reduce_group1:
subroutine     ninja_reduce_group1(scale2,tot,totr,ok)
   use iso_c_binding, only: c_ptr, c_loc, c_int
   use ninjago_module
   use p9_part21part21_part25part25part21_kinematics
   use p9_part21part21_part25part25part21_model
   use p9_part21part21_part25part25part21_d196h0l1, only: numerator_diagram196 &
   &=> numerator_ninja
   use p9_part21part21_part25part25part21_d196h0l121, only: numerator_tmu_diagr&
   &am196 => numerator_tmu
   use p9_part21part21_part25part25part21_d196h0l131, only: numerator_t3_diagram1&
   &96 => numerator_t3
   use p9_part21part21_part25part25part21_d196h0l132, only: numerator_t2_diagram1&
   &96 => numerator_t2
   use p9_part21part21_part25part25part21_globalsl1, only: epspow

   implicit none
   real(ki_nin), intent(in) :: scale2
   complex(ki_nin), dimension(-2:0), intent(out) :: tot
   complex(ki_nin), intent(out) :: totr
   logical, intent(out) :: ok

   complex(ki_nin), dimension(-2:0) :: acc
   complex(ki_nin) :: accr
   integer(c_int) :: acc_ok

   integer :: istopm, istop0

   integer, parameter :: effective_group_rank = 5

   !-----------#[ invariants for ninja:
   real(ki_nin), dimension(5,5) :: s_mat
   !-----------#] initialize invariants:
   real(ki_nin), dimension(5) :: msq
   real(ki_nin), dimension(4,5) :: Vi

   ok = .true.

   if (ninja_test.eq.1) then
      istopm = 1
      istop0 = 1
   else
      istopm = ninja_istop
      istop0 = max(2,ninja_istop)
   end if
   msq(1) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,1) = real(-k2+k5+k4, ki_nin)
   msq(2) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,2) = real(-k2+k5, ki_nin)
   msq(3) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,3) = real(-k2, ki_nin)
   msq(4) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,4) = real(0, ki_nin)
   msq(5) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,5) = real(-k3, ki_nin)
   !-----------#[ initialize invariants:
   s_mat(1,1) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(1,2) = real(-2.0_ki*mdlMT**2+mdlMh**2, ki_nin)
   s_mat(2,1) = s_mat(1,2)
   s_mat(1,3) = real(-2.0_ki*mdlMT**2+es45, ki_nin)
   s_mat(3,1) = s_mat(1,3)
   s_mat(1,4) = real(-2.0_ki*mdlMT**2+mdlMh**2+es45-es12-es23, ki_nin)
   s_mat(4,1) = s_mat(1,4)
   s_mat(1,5) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(5,1) = s_mat(1,5)
   s_mat(2,2) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(2,3) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(3,2) = s_mat(2,3)
   s_mat(2,4) = real(-2.0_ki*mdlMT**2-es12+es34-es51, ki_nin)
   s_mat(4,2) = s_mat(2,4)
   s_mat(2,5) = real(-2.0_ki*mdlMT**2+mdlMh**2+es23-es45-es51, ki_nin)
   s_mat(5,2) = s_mat(2,5)
   s_mat(3,3) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(3,4) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(4,3) = s_mat(3,4)
   s_mat(3,5) = real(-2.0_ki*mdlMT**2+es23, ki_nin)
   s_mat(5,3) = s_mat(3,5)
   s_mat(4,4) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(4,5) = real(-2.0_ki*mdlMT**2+mdlMh**2, ki_nin)
   s_mat(5,4) = s_mat(4,5)
   s_mat(5,5) = real(-2.0_ki*mdlMT**2, ki_nin)
   !-----------#] initialize invariants


      !------#[ sum over reduction of single diagrams:
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='196'>"
         end if
         call ninja_diagram(numerator_diagram196, numerator_t3_diagram196, nume&
         &rator_t2_diagram196, &
          &  numerator_tmu_diagram196, &
          & 5, 5, 5, (/1,2,3,4,5/), &
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
end subroutine ninja_reduce_group1
!-----#] subroutine ninja_reduce_group1:
!-----#[ subroutine ninja_reduce_group2:
subroutine     ninja_reduce_group2(scale2,tot,totr,ok)
   use iso_c_binding, only: c_ptr, c_loc, c_int
   use ninjago_module
   use p9_part21part21_part25part25part21_kinematics
   use p9_part21part21_part25part25part21_model
   use p9_part21part21_part25part25part21_d121h0l1, only: numerator_diagram121 &
   &=> numerator_ninja
   use p9_part21part21_part25part25part21_d121h0l121, only: numerator_tmu_diagr&
   &am121 => numerator_tmu
   use p9_part21part21_part25part25part21_d121h0l131, only: numerator_t3_diagram1&
   &21 => numerator_t3
   use p9_part21part21_part25part25part21_d121h0l132, only: numerator_t2_diagram1&
   &21 => numerator_t2
   use p9_part21part21_part25part25part21_d204h0l1, only: numerator_diagram204 &
   &=> numerator_ninja
   use p9_part21part21_part25part25part21_d204h0l121, only: numerator_tmu_diagr&
   &am204 => numerator_tmu
   use p9_part21part21_part25part25part21_d204h0l131, only: numerator_t3_diagram2&
   &04 => numerator_t3
   use p9_part21part21_part25part25part21_d204h0l132, only: numerator_t2_diagram2&
   &04 => numerator_t2
   use p9_part21part21_part25part25part21_globalsl1, only: epspow

   implicit none
   real(ki_nin), intent(in) :: scale2
   complex(ki_nin), dimension(-2:0), intent(out) :: tot
   complex(ki_nin), intent(out) :: totr
   logical, intent(out) :: ok

   complex(ki_nin), dimension(-2:0) :: acc
   complex(ki_nin) :: accr
   integer(c_int) :: acc_ok

   integer :: istopm, istop0

   integer, parameter :: effective_group_rank = 5

   !-----------#[ invariants for ninja:
   real(ki_nin), dimension(5,5) :: s_mat
   !-----------#] initialize invariants:
   real(ki_nin), dimension(5) :: msq
   real(ki_nin), dimension(4,5) :: Vi

   ok = .true.

   if (ninja_test.eq.1) then
      istopm = 1
      istop0 = 1
   else
      istopm = ninja_istop
      istop0 = max(2,ninja_istop)
   end if
   msq(1) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,1) = real(-k3-k5, ki_nin)
   msq(2) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,2) = real(-k3, ki_nin)
   msq(3) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,3) = real(0, ki_nin)
   msq(4) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,4) = real(-k2, ki_nin)
   msq(5) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,5) = real(-k2+k4, ki_nin)
   !-----------#[ initialize invariants:
   s_mat(1,1) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(1,2) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(2,1) = s_mat(1,2)
   s_mat(1,3) = real(-2.0_ki*mdlMT**2+es12+2.0_ki*mdlMh**2-es34-es45, ki_nin)
   s_mat(3,1) = s_mat(1,3)
   s_mat(1,4) = real(-2.0_ki*mdlMT**2+mdlMh**2+es23-es45-es51, ki_nin)
   s_mat(4,1) = s_mat(1,4)
   s_mat(1,5) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(5,1) = s_mat(1,5)
   s_mat(2,2) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(2,3) = real(-2.0_ki*mdlMT**2+mdlMh**2, ki_nin)
   s_mat(3,2) = s_mat(2,3)
   s_mat(2,4) = real(-2.0_ki*mdlMT**2+es23, ki_nin)
   s_mat(4,2) = s_mat(2,4)
   s_mat(2,5) = real(-2.0_ki*mdlMT**2+es51, ki_nin)
   s_mat(5,2) = s_mat(2,5)
   s_mat(3,3) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(3,4) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(4,3) = s_mat(3,4)
   s_mat(3,5) = real(-2.0_ki*mdlMT**2+2.0_ki*mdlMh**2-es23+es51-es34, ki_nin)
   s_mat(5,3) = s_mat(3,5)
   s_mat(4,4) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(4,5) = real(-2.0_ki*mdlMT**2+mdlMh**2, ki_nin)
   s_mat(5,4) = s_mat(4,5)
   s_mat(5,5) = real(-2.0_ki*mdlMT**2, ki_nin)
   !-----------#] initialize invariants


      !------#[ sum over reduction of single diagrams:
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='121'>"
         end if
         call ninja_diagram(numerator_diagram121, numerator_t3_diagram121, nume&
         &rator_t2_diagram121, &
          &  numerator_tmu_diagram121, &
          & 5, 4, 4, (/2,3,4,5/), &
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
            write(logfile,*) "<diagram index='204'>"
         end if
         call ninja_diagram(numerator_diagram204, numerator_t3_diagram204, nume&
         &rator_t2_diagram204, &
          &  numerator_tmu_diagram204, &
          & 5, 5, 5, (/1,2,3,4,5/), &
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
end subroutine ninja_reduce_group2
!-----#] subroutine ninja_reduce_group2:
!-----#[ subroutine ninja_reduce_group3:
subroutine     ninja_reduce_group3(scale2,tot,totr,ok)
   use iso_c_binding, only: c_ptr, c_loc, c_int
   use ninjago_module
   use p9_part21part21_part25part25part21_kinematics
   use p9_part21part21_part25part25part21_model
   use p9_part21part21_part25part25part21_d200h0l1, only: numerator_diagram200 &
   &=> numerator_ninja
   use p9_part21part21_part25part25part21_d200h0l121, only: numerator_tmu_diagr&
   &am200 => numerator_tmu
   use p9_part21part21_part25part25part21_d200h0l131, only: numerator_t3_diagram2&
   &00 => numerator_t3
   use p9_part21part21_part25part25part21_d200h0l132, only: numerator_t2_diagram2&
   &00 => numerator_t2
   use p9_part21part21_part25part25part21_globalsl1, only: epspow

   implicit none
   real(ki_nin), intent(in) :: scale2
   complex(ki_nin), dimension(-2:0), intent(out) :: tot
   complex(ki_nin), intent(out) :: totr
   logical, intent(out) :: ok

   complex(ki_nin), dimension(-2:0) :: acc
   complex(ki_nin) :: accr
   integer(c_int) :: acc_ok

   integer :: istopm, istop0

   integer, parameter :: effective_group_rank = 5

   !-----------#[ invariants for ninja:
   real(ki_nin), dimension(5,5) :: s_mat
   !-----------#] initialize invariants:
   real(ki_nin), dimension(5) :: msq
   real(ki_nin), dimension(4,5) :: Vi

   ok = .true.

   if (ninja_test.eq.1) then
      istopm = 1
      istop0 = 1
   else
      istopm = ninja_istop
      istop0 = max(2,ninja_istop)
   end if
   msq(1) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,1) = real(-k2+k5+k4, ki_nin)
   msq(2) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,2) = real(-k2+k4, ki_nin)
   msq(3) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,3) = real(-k2, ki_nin)
   msq(4) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,4) = real(0, ki_nin)
   msq(5) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,5) = real(-k3, ki_nin)
   !-----------#[ initialize invariants:
   s_mat(1,1) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(1,2) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(2,1) = s_mat(1,2)
   s_mat(1,3) = real(-2.0_ki*mdlMT**2+es45, ki_nin)
   s_mat(3,1) = s_mat(1,3)
   s_mat(1,4) = real(-2.0_ki*mdlMT**2+mdlMh**2+es45-es12-es23, ki_nin)
   s_mat(4,1) = s_mat(1,4)
   s_mat(1,5) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(5,1) = s_mat(1,5)
   s_mat(2,2) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(2,3) = real(-2.0_ki*mdlMT**2+mdlMh**2, ki_nin)
   s_mat(3,2) = s_mat(2,3)
   s_mat(2,4) = real(-2.0_ki*mdlMT**2+2.0_ki*mdlMh**2-es23+es51-es34, ki_nin)
   s_mat(4,2) = s_mat(2,4)
   s_mat(2,5) = real(-2.0_ki*mdlMT**2+es51, ki_nin)
   s_mat(5,2) = s_mat(2,5)
   s_mat(3,3) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(3,4) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(4,3) = s_mat(3,4)
   s_mat(3,5) = real(-2.0_ki*mdlMT**2+es23, ki_nin)
   s_mat(5,3) = s_mat(3,5)
   s_mat(4,4) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(4,5) = real(-2.0_ki*mdlMT**2+mdlMh**2, ki_nin)
   s_mat(5,4) = s_mat(4,5)
   s_mat(5,5) = real(-2.0_ki*mdlMT**2, ki_nin)
   !-----------#] initialize invariants


      !------#[ sum over reduction of single diagrams:
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='200'>"
         end if
         call ninja_diagram(numerator_diagram200, numerator_t3_diagram200, nume&
         &rator_t2_diagram200, &
          &  numerator_tmu_diagram200, &
          & 5, 5, 5, (/1,2,3,4,5/), &
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
end subroutine ninja_reduce_group3
!-----#] subroutine ninja_reduce_group3:
!-----#[ subroutine ninja_reduce_group4:
subroutine     ninja_reduce_group4(scale2,tot,totr,ok)
   use iso_c_binding, only: c_ptr, c_loc, c_int
   use ninjago_module
   use p9_part21part21_part25part25part21_kinematics
   use p9_part21part21_part25part25part21_model
   use p9_part21part21_part25part25part21_d192h0l1, only: numerator_diagram192 &
   &=> numerator_ninja
   use p9_part21part21_part25part25part21_d192h0l121, only: numerator_tmu_diagr&
   &am192 => numerator_tmu
   use p9_part21part21_part25part25part21_d192h0l131, only: numerator_t3_diagram1&
   &92 => numerator_t3
   use p9_part21part21_part25part25part21_d192h0l132, only: numerator_t2_diagram1&
   &92 => numerator_t2
   use p9_part21part21_part25part25part21_globalsl1, only: epspow

   implicit none
   real(ki_nin), intent(in) :: scale2
   complex(ki_nin), dimension(-2:0), intent(out) :: tot
   complex(ki_nin), intent(out) :: totr
   logical, intent(out) :: ok

   complex(ki_nin), dimension(-2:0) :: acc
   complex(ki_nin) :: accr
   integer(c_int) :: acc_ok

   integer :: istopm, istop0

   integer, parameter :: effective_group_rank = 5

   !-----------#[ invariants for ninja:
   real(ki_nin), dimension(5,5) :: s_mat
   !-----------#] initialize invariants:
   real(ki_nin), dimension(5) :: msq
   real(ki_nin), dimension(4,5) :: Vi

   ok = .true.

   if (ninja_test.eq.1) then
      istopm = 1
      istop0 = 1
   else
      istopm = ninja_istop
      istop0 = max(2,ninja_istop)
   end if
   msq(1) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,1) = real(-k3-k5-k4, ki_nin)
   msq(2) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,2) = real(-k3-k4, ki_nin)
   msq(3) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,3) = real(-k3, ki_nin)
   msq(4) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,4) = real(0, ki_nin)
   msq(5) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,5) = real(-k2, ki_nin)
   !-----------#[ initialize invariants:
   s_mat(1,1) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(1,2) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(2,1) = s_mat(1,2)
   s_mat(1,3) = real(-2.0_ki*mdlMT**2+es45, ki_nin)
   s_mat(3,1) = s_mat(1,3)
   s_mat(1,4) = real(-2.0_ki*mdlMT**2+es12, ki_nin)
   s_mat(4,1) = s_mat(1,4)
   s_mat(1,5) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(5,1) = s_mat(1,5)
   s_mat(2,2) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(2,3) = real(-2.0_ki*mdlMT**2+mdlMh**2, ki_nin)
   s_mat(3,2) = s_mat(2,3)
   s_mat(2,4) = real(-2.0_ki*mdlMT**2+es34, ki_nin)
   s_mat(4,2) = s_mat(2,4)
   s_mat(2,5) = real(-2.0_ki*mdlMT**2+es51, ki_nin)
   s_mat(5,2) = s_mat(2,5)
   s_mat(3,3) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(3,4) = real(-2.0_ki*mdlMT**2+mdlMh**2, ki_nin)
   s_mat(4,3) = s_mat(3,4)
   s_mat(3,5) = real(-2.0_ki*mdlMT**2+es23, ki_nin)
   s_mat(5,3) = s_mat(3,5)
   s_mat(4,4) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(4,5) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(5,4) = s_mat(4,5)
   s_mat(5,5) = real(-2.0_ki*mdlMT**2, ki_nin)
   !-----------#] initialize invariants


      !------#[ sum over reduction of single diagrams:
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='192'>"
         end if
         call ninja_diagram(numerator_diagram192, numerator_t3_diagram192, nume&
         &rator_t2_diagram192, &
          &  numerator_tmu_diagram192, &
          & 5, 5, 5, (/1,2,3,4,5/), &
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
end subroutine ninja_reduce_group4
!-----#] subroutine ninja_reduce_group4:
!-----#[ subroutine ninja_reduce_group5:
subroutine     ninja_reduce_group5(scale2,tot,totr,ok)
   use iso_c_binding, only: c_ptr, c_loc, c_int
   use ninjago_module
   use p9_part21part21_part25part25part21_kinematics
   use p9_part21part21_part25part25part21_model
   use p9_part21part21_part25part25part21_d188h0l1, only: numerator_diagram188 &
   &=> numerator_ninja
   use p9_part21part21_part25part25part21_d188h0l121, only: numerator_tmu_diagr&
   &am188 => numerator_tmu
   use p9_part21part21_part25part25part21_d188h0l131, only: numerator_t3_diagram1&
   &88 => numerator_t3
   use p9_part21part21_part25part25part21_d188h0l132, only: numerator_t2_diagram1&
   &88 => numerator_t2
   use p9_part21part21_part25part25part21_globalsl1, only: epspow

   implicit none
   real(ki_nin), intent(in) :: scale2
   complex(ki_nin), dimension(-2:0), intent(out) :: tot
   complex(ki_nin), intent(out) :: totr
   logical, intent(out) :: ok

   complex(ki_nin), dimension(-2:0) :: acc
   complex(ki_nin) :: accr
   integer(c_int) :: acc_ok

   integer :: istopm, istop0

   integer, parameter :: effective_group_rank = 5

   !-----------#[ invariants for ninja:
   real(ki_nin), dimension(5,5) :: s_mat
   !-----------#] initialize invariants:
   real(ki_nin), dimension(5) :: msq
   real(ki_nin), dimension(4,5) :: Vi

   ok = .true.

   if (ninja_test.eq.1) then
      istopm = 1
      istop0 = 1
   else
      istopm = ninja_istop
      istop0 = max(2,ninja_istop)
   end if
   msq(1) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,1) = real(-k3-k5-k4, ki_nin)
   msq(2) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,2) = real(-k3-k5, ki_nin)
   msq(3) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,3) = real(-k3, ki_nin)
   msq(4) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,4) = real(0, ki_nin)
   msq(5) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,5) = real(-k2, ki_nin)
   !-----------#[ initialize invariants:
   s_mat(1,1) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(1,2) = real(-2.0_ki*mdlMT**2+mdlMh**2, ki_nin)
   s_mat(2,1) = s_mat(1,2)
   s_mat(1,3) = real(-2.0_ki*mdlMT**2+es45, ki_nin)
   s_mat(3,1) = s_mat(1,3)
   s_mat(1,4) = real(-2.0_ki*mdlMT**2+es12, ki_nin)
   s_mat(4,1) = s_mat(1,4)
   s_mat(1,5) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(5,1) = s_mat(1,5)
   s_mat(2,2) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(2,3) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(3,2) = s_mat(2,3)
   s_mat(2,4) = real(-2.0_ki*mdlMT**2+es12+2.0_ki*mdlMh**2-es34-es45, ki_nin)
   s_mat(4,2) = s_mat(2,4)
   s_mat(2,5) = real(-2.0_ki*mdlMT**2+mdlMh**2+es23-es45-es51, ki_nin)
   s_mat(5,2) = s_mat(2,5)
   s_mat(3,3) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(3,4) = real(-2.0_ki*mdlMT**2+mdlMh**2, ki_nin)
   s_mat(4,3) = s_mat(3,4)
   s_mat(3,5) = real(-2.0_ki*mdlMT**2+es23, ki_nin)
   s_mat(5,3) = s_mat(3,5)
   s_mat(4,4) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(4,5) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(5,4) = s_mat(4,5)
   s_mat(5,5) = real(-2.0_ki*mdlMT**2, ki_nin)
   !-----------#] initialize invariants


      !------#[ sum over reduction of single diagrams:
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='188'>"
         end if
         call ninja_diagram(numerator_diagram188, numerator_t3_diagram188, nume&
         &rator_t2_diagram188, &
          &  numerator_tmu_diagram188, &
          & 5, 5, 5, (/1,2,3,4,5/), &
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
end subroutine ninja_reduce_group5
!-----#] subroutine ninja_reduce_group5:
!-----#[ subroutine ninja_reduce_group6:
subroutine     ninja_reduce_group6(scale2,tot,totr,ok)
   use iso_c_binding, only: c_ptr, c_loc, c_int
   use ninjago_module
   use p9_part21part21_part25part25part21_kinematics
   use p9_part21part21_part25part25part21_model
   use p9_part21part21_part25part25part21_d137h0l1, only: numerator_diagram137 &
   &=> numerator_ninja
   use p9_part21part21_part25part25part21_d137h0l121, only: numerator_tmu_diagr&
   &am137 => numerator_tmu
   use p9_part21part21_part25part25part21_d137h0l131, only: numerator_t3_diagram1&
   &37 => numerator_t3
   use p9_part21part21_part25part25part21_d137h0l132, only: numerator_t2_diagram1&
   &37 => numerator_t2
   use p9_part21part21_part25part25part21_d202h0l1, only: numerator_diagram202 &
   &=> numerator_ninja
   use p9_part21part21_part25part25part21_d202h0l121, only: numerator_tmu_diagr&
   &am202 => numerator_tmu
   use p9_part21part21_part25part25part21_d202h0l131, only: numerator_t3_diagram2&
   &02 => numerator_t3
   use p9_part21part21_part25part25part21_d202h0l132, only: numerator_t2_diagram2&
   &02 => numerator_t2
   use p9_part21part21_part25part25part21_globalsl1, only: epspow

   implicit none
   real(ki_nin), intent(in) :: scale2
   complex(ki_nin), dimension(-2:0), intent(out) :: tot
   complex(ki_nin), intent(out) :: totr
   logical, intent(out) :: ok

   complex(ki_nin), dimension(-2:0) :: acc
   complex(ki_nin) :: accr
   integer(c_int) :: acc_ok

   integer :: istopm, istop0

   integer, parameter :: effective_group_rank = 5

   !-----------#[ invariants for ninja:
   real(ki_nin), dimension(5,5) :: s_mat
   !-----------#] initialize invariants:
   real(ki_nin), dimension(5) :: msq
   real(ki_nin), dimension(4,5) :: Vi

   ok = .true.

   if (ninja_test.eq.1) then
      istopm = 1
      istop0 = 1
   else
      istopm = ninja_istop
      istop0 = max(2,ninja_istop)
   end if
   msq(1) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,1) = real(-k2+k5, ki_nin)
   msq(2) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,2) = real(-k2, ki_nin)
   msq(3) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,3) = real(0, ki_nin)
   msq(4) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,4) = real(-k4, ki_nin)
   msq(5) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,5) = real(-k3-k4, ki_nin)
   !-----------#[ initialize invariants:
   s_mat(1,1) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(1,2) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(2,1) = s_mat(1,2)
   s_mat(1,3) = real(-2.0_ki*mdlMT**2-es12+es34-es51, ki_nin)
   s_mat(3,1) = s_mat(1,3)
   s_mat(1,4) = real(-2.0_ki*mdlMT**2+mdlMh**2+es45-es12-es23, ki_nin)
   s_mat(4,1) = s_mat(1,4)
   s_mat(1,5) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(5,1) = s_mat(1,5)
   s_mat(2,2) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(2,3) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(3,2) = s_mat(2,3)
   s_mat(2,4) = real(-2.0_ki*mdlMT**2+2.0_ki*mdlMh**2-es23+es51-es34, ki_nin)
   s_mat(4,2) = s_mat(2,4)
   s_mat(2,5) = real(-2.0_ki*mdlMT**2+es51, ki_nin)
   s_mat(5,2) = s_mat(2,5)
   s_mat(3,3) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(3,4) = real(-2.0_ki*mdlMT**2+mdlMh**2, ki_nin)
   s_mat(4,3) = s_mat(3,4)
   s_mat(3,5) = real(-2.0_ki*mdlMT**2+es34, ki_nin)
   s_mat(5,3) = s_mat(3,5)
   s_mat(4,4) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(4,5) = real(-2.0_ki*mdlMT**2+mdlMh**2, ki_nin)
   s_mat(5,4) = s_mat(4,5)
   s_mat(5,5) = real(-2.0_ki*mdlMT**2, ki_nin)
   !-----------#] initialize invariants


      !------#[ sum over reduction of single diagrams:
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='137'>"
         end if
         call ninja_diagram(numerator_diagram137, numerator_t3_diagram137, nume&
         &rator_t2_diagram137, &
          &  numerator_tmu_diagram137, &
          & 5, 4, 4, (/1,3,4,5/), &
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
            write(logfile,*) "<diagram index='202'>"
         end if
         call ninja_diagram(numerator_diagram202, numerator_t3_diagram202, nume&
         &rator_t2_diagram202, &
          &  numerator_tmu_diagram202, &
          & 5, 5, 5, (/1,2,3,4,5/), &
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
end subroutine ninja_reduce_group6
!-----#] subroutine ninja_reduce_group6:
!-----#[ subroutine ninja_reduce_group7:
subroutine     ninja_reduce_group7(scale2,tot,totr,ok)
   use iso_c_binding, only: c_ptr, c_loc, c_int
   use ninjago_module
   use p9_part21part21_part25part25part21_kinematics
   use p9_part21part21_part25part25part21_model
   use p9_part21part21_part25part25part21_d133h0l1, only: numerator_diagram133 &
   &=> numerator_ninja
   use p9_part21part21_part25part25part21_d133h0l121, only: numerator_tmu_diagr&
   &am133 => numerator_tmu
   use p9_part21part21_part25part25part21_d133h0l131, only: numerator_t3_diagram1&
   &33 => numerator_t3
   use p9_part21part21_part25part25part21_d133h0l132, only: numerator_t2_diagram1&
   &33 => numerator_t2
   use p9_part21part21_part25part25part21_d198h0l1, only: numerator_diagram198 &
   &=> numerator_ninja
   use p9_part21part21_part25part25part21_d198h0l121, only: numerator_tmu_diagr&
   &am198 => numerator_tmu
   use p9_part21part21_part25part25part21_d198h0l131, only: numerator_t3_diagram1&
   &98 => numerator_t3
   use p9_part21part21_part25part25part21_d198h0l132, only: numerator_t2_diagram1&
   &98 => numerator_t2
   use p9_part21part21_part25part25part21_globalsl1, only: epspow

   implicit none
   real(ki_nin), intent(in) :: scale2
   complex(ki_nin), dimension(-2:0), intent(out) :: tot
   complex(ki_nin), intent(out) :: totr
   logical, intent(out) :: ok

   complex(ki_nin), dimension(-2:0) :: acc
   complex(ki_nin) :: accr
   integer(c_int) :: acc_ok

   integer :: istopm, istop0

   integer, parameter :: effective_group_rank = 5

   !-----------#[ invariants for ninja:
   real(ki_nin), dimension(5,5) :: s_mat
   !-----------#] initialize invariants:
   real(ki_nin), dimension(5) :: msq
   real(ki_nin), dimension(4,5) :: Vi

   ok = .true.

   if (ninja_test.eq.1) then
      istopm = 1
      istop0 = 1
   else
      istopm = ninja_istop
      istop0 = max(2,ninja_istop)
   end if
   msq(1) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,1) = real(-k4, ki_nin)
   msq(2) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,2) = real(0, ki_nin)
   msq(3) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,3) = real(-k2, ki_nin)
   msq(4) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,4) = real(-k2+k5, ki_nin)
   msq(5) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,5) = real(k3-k2+k5, ki_nin)
   !-----------#[ initialize invariants:
   s_mat(1,1) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(1,2) = real(-2.0_ki*mdlMT**2+mdlMh**2, ki_nin)
   s_mat(2,1) = s_mat(1,2)
   s_mat(1,3) = real(-2.0_ki*mdlMT**2+2.0_ki*mdlMh**2-es23+es51-es34, ki_nin)
   s_mat(3,1) = s_mat(1,3)
   s_mat(1,4) = real(-2.0_ki*mdlMT**2+mdlMh**2+es45-es12-es23, ki_nin)
   s_mat(4,1) = s_mat(1,4)
   s_mat(1,5) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(5,1) = s_mat(1,5)
   s_mat(2,2) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(2,3) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(3,2) = s_mat(2,3)
   s_mat(2,4) = real(-2.0_ki*mdlMT**2-es12+es34-es51, ki_nin)
   s_mat(4,2) = s_mat(2,4)
   s_mat(2,5) = real(-2.0_ki*mdlMT**2+mdlMh**2+es23-es45-es51, ki_nin)
   s_mat(5,2) = s_mat(2,5)
   s_mat(3,3) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(3,4) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(4,3) = s_mat(3,4)
   s_mat(3,5) = real(-2.0_ki*mdlMT**2+es12+2.0_ki*mdlMh**2-es34-es45, ki_nin)
   s_mat(5,3) = s_mat(3,5)
   s_mat(4,4) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(4,5) = real(-2.0_ki*mdlMT**2+mdlMh**2, ki_nin)
   s_mat(5,4) = s_mat(4,5)
   s_mat(5,5) = real(-2.0_ki*mdlMT**2, ki_nin)
   !-----------#] initialize invariants


      !------#[ sum over reduction of single diagrams:
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='133'>"
         end if
         call ninja_diagram(numerator_diagram133, numerator_t3_diagram133, nume&
         &rator_t2_diagram133, &
          &  numerator_tmu_diagram133, &
          & 5, 4, 4, (/1,2,4,5/), &
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
            write(logfile,*) "<diagram index='198'>"
         end if
         call ninja_diagram(numerator_diagram198, numerator_t3_diagram198, nume&
         &rator_t2_diagram198, &
          &  numerator_tmu_diagram198, &
          & 5, 5, 5, (/1,2,3,4,5/), &
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
end subroutine ninja_reduce_group7
!-----#] subroutine ninja_reduce_group7:
!-----#[ subroutine ninja_reduce_group8:
subroutine     ninja_reduce_group8(scale2,tot,totr,ok)
   use iso_c_binding, only: c_ptr, c_loc, c_int
   use ninjago_module
   use p9_part21part21_part25part25part21_kinematics
   use p9_part21part21_part25part25part21_model
   use p9_part21part21_part25part25part21_d52h0l1, only: numerator_diagram52 =>&
   & numerator_ninja
   use p9_part21part21_part25part25part21_d52h0l121, only: numerator_tmu_diagram5&
   &2 => numerator_tmu
   use p9_part21part21_part25part25part21_d52h0l131, only: numerator_t3_diagram5&
   &2 => numerator_t3
   use p9_part21part21_part25part25part21_d52h0l132, only: numerator_t2_diagram5&
   &2 => numerator_t2
   use p9_part21part21_part25part25part21_d129h0l1, only: numerator_diagram129 &
   &=> numerator_ninja
   use p9_part21part21_part25part25part21_d129h0l121, only: numerator_tmu_diagr&
   &am129 => numerator_tmu
   use p9_part21part21_part25part25part21_d129h0l131, only: numerator_t3_diagram1&
   &29 => numerator_t3
   use p9_part21part21_part25part25part21_d129h0l132, only: numerator_t2_diagram1&
   &29 => numerator_t2
   use p9_part21part21_part25part25part21_d145h0l1, only: numerator_diagram145 &
   &=> numerator_ninja
   use p9_part21part21_part25part25part21_d145h0l121, only: numerator_tmu_diagr&
   &am145 => numerator_tmu
   use p9_part21part21_part25part25part21_d145h0l131, only: numerator_t3_diagram1&
   &45 => numerator_t3
   use p9_part21part21_part25part25part21_d145h0l132, only: numerator_t2_diagram1&
   &45 => numerator_t2
   use p9_part21part21_part25part25part21_d194h0l1, only: numerator_diagram194 &
   &=> numerator_ninja
   use p9_part21part21_part25part25part21_d194h0l121, only: numerator_tmu_diagr&
   &am194 => numerator_tmu
   use p9_part21part21_part25part25part21_d194h0l131, only: numerator_t3_diagram1&
   &94 => numerator_t3
   use p9_part21part21_part25part25part21_d194h0l132, only: numerator_t2_diagram1&
   &94 => numerator_t2
   use p9_part21part21_part25part25part21_globalsl1, only: epspow

   implicit none
   real(ki_nin), intent(in) :: scale2
   complex(ki_nin), dimension(-2:0), intent(out) :: tot
   complex(ki_nin), intent(out) :: totr
   logical, intent(out) :: ok

   complex(ki_nin), dimension(-2:0) :: acc
   complex(ki_nin) :: accr
   integer(c_int) :: acc_ok

   integer :: istopm, istop0

   integer, parameter :: effective_group_rank = 5

   !-----------#[ invariants for ninja:
   real(ki_nin), dimension(5,5) :: s_mat
   !-----------#] initialize invariants:
   real(ki_nin), dimension(5) :: msq
   real(ki_nin), dimension(4,5) :: Vi

   ok = .true.

   if (ninja_test.eq.1) then
      istopm = 1
      istop0 = 1
   else
      istopm = ninja_istop
      istop0 = max(2,ninja_istop)
   end if
   msq(1) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,1) = real(-k3-k5-k4, ki_nin)
   msq(2) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,2) = real(-k3-k4, ki_nin)
   msq(3) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,3) = real(-k4, ki_nin)
   msq(4) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,4) = real(0, ki_nin)
   msq(5) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,5) = real(-k2, ki_nin)
   !-----------#[ initialize invariants:
   s_mat(1,1) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(1,2) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(2,1) = s_mat(1,2)
   s_mat(1,3) = real(-2.0_ki*mdlMT**2+es12+2.0_ki*mdlMh**2-es34-es45, ki_nin)
   s_mat(3,1) = s_mat(1,3)
   s_mat(1,4) = real(-2.0_ki*mdlMT**2+es12, ki_nin)
   s_mat(4,1) = s_mat(1,4)
   s_mat(1,5) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(5,1) = s_mat(1,5)
   s_mat(2,2) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(2,3) = real(-2.0_ki*mdlMT**2+mdlMh**2, ki_nin)
   s_mat(3,2) = s_mat(2,3)
   s_mat(2,4) = real(-2.0_ki*mdlMT**2+es34, ki_nin)
   s_mat(4,2) = s_mat(2,4)
   s_mat(2,5) = real(-2.0_ki*mdlMT**2+es51, ki_nin)
   s_mat(5,2) = s_mat(2,5)
   s_mat(3,3) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(3,4) = real(-2.0_ki*mdlMT**2+mdlMh**2, ki_nin)
   s_mat(4,3) = s_mat(3,4)
   s_mat(3,5) = real(-2.0_ki*mdlMT**2+2.0_ki*mdlMh**2-es23+es51-es34, ki_nin)
   s_mat(5,3) = s_mat(3,5)
   s_mat(4,4) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(4,5) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(5,4) = s_mat(4,5)
   s_mat(5,5) = real(-2.0_ki*mdlMT**2, ki_nin)
   !-----------#] initialize invariants


      !------#[ sum over reduction of single diagrams:
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='52'>"
         end if
         call ninja_diagram(numerator_diagram52, numerator_t3_diagram52, numera&
         &tor_t2_diagram52, &
          &  numerator_tmu_diagram52, &
          & 5, 4, 4, (/1,2,4,5/), &
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
            write(logfile,*) "<diagram index='129'>"
         end if
         call ninja_diagram(numerator_diagram129, numerator_t3_diagram129, nume&
         &rator_t2_diagram129, &
          &  numerator_tmu_diagram129, &
          & 5, 4, 4, (/2,3,4,5/), &
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
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='145'>"
         end if
         call ninja_diagram(numerator_diagram145, numerator_t3_diagram145, nume&
         &rator_t2_diagram145, &
          &  numerator_tmu_diagram145, &
          & 5, 4, 4, (/1,2,4,5/), &
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
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='194'>"
         end if
         call ninja_diagram(numerator_diagram194, numerator_t3_diagram194, nume&
         &rator_t2_diagram194, &
          &  numerator_tmu_diagram194, &
          & 5, 5, 5, (/1,2,3,4,5/), &
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
end subroutine ninja_reduce_group8
!-----#] subroutine ninja_reduce_group8:
!-----#[ subroutine ninja_reduce_group9:
subroutine     ninja_reduce_group9(scale2,tot,totr,ok)
   use iso_c_binding, only: c_ptr, c_loc, c_int
   use ninjago_module
   use p9_part21part21_part25part25part21_kinematics
   use p9_part21part21_part25part25part21_model
   use p9_part21part21_part25part25part21_d117h0l1, only: numerator_diagram117 &
   &=> numerator_ninja
   use p9_part21part21_part25part25part21_d117h0l121, only: numerator_tmu_diagr&
   &am117 => numerator_tmu
   use p9_part21part21_part25part25part21_d117h0l131, only: numerator_t3_diagram1&
   &17 => numerator_t3
   use p9_part21part21_part25part25part21_d117h0l132, only: numerator_t2_diagram1&
   &17 => numerator_t2
   use p9_part21part21_part25part25part21_d184h0l1, only: numerator_diagram184 &
   &=> numerator_ninja
   use p9_part21part21_part25part25part21_d184h0l121, only: numerator_tmu_diagr&
   &am184 => numerator_tmu
   use p9_part21part21_part25part25part21_d184h0l131, only: numerator_t3_diagram1&
   &84 => numerator_t3
   use p9_part21part21_part25part25part21_d184h0l132, only: numerator_t2_diagram1&
   &84 => numerator_t2
   use p9_part21part21_part25part25part21_globalsl1, only: epspow

   implicit none
   real(ki_nin), intent(in) :: scale2
   complex(ki_nin), dimension(-2:0), intent(out) :: tot
   complex(ki_nin), intent(out) :: totr
   logical, intent(out) :: ok

   complex(ki_nin), dimension(-2:0) :: acc
   complex(ki_nin) :: accr
   integer(c_int) :: acc_ok

   integer :: istopm, istop0

   integer, parameter :: effective_group_rank = 5

   !-----------#[ invariants for ninja:
   real(ki_nin), dimension(5,5) :: s_mat
   !-----------#] initialize invariants:
   real(ki_nin), dimension(5) :: msq
   real(ki_nin), dimension(4,5) :: Vi

   ok = .true.

   if (ninja_test.eq.1) then
      istopm = 1
      istop0 = 1
   else
      istopm = ninja_istop
      istop0 = max(2,ninja_istop)
   end if
   msq(1) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,1) = real(-k3-k5-k4, ki_nin)
   msq(2) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,2) = real(-k5-k4, ki_nin)
   msq(3) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,3) = real(-k4, ki_nin)
   msq(4) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,4) = real(0, ki_nin)
   msq(5) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,5) = real(-k2, ki_nin)
   !-----------#[ initialize invariants:
   s_mat(1,1) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(1,2) = real(-2.0_ki*mdlMT**2+mdlMh**2, ki_nin)
   s_mat(2,1) = s_mat(1,2)
   s_mat(1,3) = real(-2.0_ki*mdlMT**2+es12+2.0_ki*mdlMh**2-es34-es45, ki_nin)
   s_mat(3,1) = s_mat(1,3)
   s_mat(1,4) = real(-2.0_ki*mdlMT**2+es12, ki_nin)
   s_mat(4,1) = s_mat(1,4)
   s_mat(1,5) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(5,1) = s_mat(1,5)
   s_mat(2,2) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(2,3) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(3,2) = s_mat(2,3)
   s_mat(2,4) = real(-2.0_ki*mdlMT**2+es45, ki_nin)
   s_mat(4,2) = s_mat(2,4)
   s_mat(2,5) = real(-2.0_ki*mdlMT**2+mdlMh**2+es45-es12-es23, ki_nin)
   s_mat(5,2) = s_mat(2,5)
   s_mat(3,3) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(3,4) = real(-2.0_ki*mdlMT**2+mdlMh**2, ki_nin)
   s_mat(4,3) = s_mat(3,4)
   s_mat(3,5) = real(-2.0_ki*mdlMT**2+2.0_ki*mdlMh**2-es23+es51-es34, ki_nin)
   s_mat(5,3) = s_mat(3,5)
   s_mat(4,4) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(4,5) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(5,4) = s_mat(4,5)
   s_mat(5,5) = real(-2.0_ki*mdlMT**2, ki_nin)
   !-----------#] initialize invariants


      !------#[ sum over reduction of single diagrams:
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='117'>"
         end if
         call ninja_diagram(numerator_diagram117, numerator_t3_diagram117, nume&
         &rator_t2_diagram117, &
          &  numerator_tmu_diagram117, &
          & 5, 4, 4, (/1,2,3,4/), &
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
            write(logfile,*) "<diagram index='184'>"
         end if
         call ninja_diagram(numerator_diagram184, numerator_t3_diagram184, nume&
         &rator_t2_diagram184, &
          &  numerator_tmu_diagram184, &
          & 5, 5, 5, (/1,2,3,4,5/), &
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
end subroutine ninja_reduce_group9
!-----#] subroutine ninja_reduce_group9:
!-----#[ subroutine ninja_reduce_group10:
subroutine     ninja_reduce_group10(scale2,tot,totr,ok)
   use iso_c_binding, only: c_ptr, c_loc, c_int
   use ninjago_module
   use p9_part21part21_part25part25part21_kinematics
   use p9_part21part21_part25part25part21_model
   use p9_part21part21_part25part25part21_d34h0l1, only: numerator_diagram34 =>&
   & numerator_ninja
   use p9_part21part21_part25part25part21_d34h0l121, only: numerator_tmu_diagram3&
   &4 => numerator_tmu
   use p9_part21part21_part25part25part21_d34h0l131, only: numerator_t3_diagram3&
   &4 => numerator_t3
   use p9_part21part21_part25part25part21_d34h0l132, only: numerator_t2_diagram3&
   &4 => numerator_t2
   use p9_part21part21_part25part25part21_d42h0l1, only: numerator_diagram42 =>&
   & numerator_ninja
   use p9_part21part21_part25part25part21_d42h0l121, only: numerator_tmu_diagram4&
   &2 => numerator_tmu
   use p9_part21part21_part25part25part21_d42h0l131, only: numerator_t3_diagram4&
   &2 => numerator_t3
   use p9_part21part21_part25part25part21_d42h0l132, only: numerator_t2_diagram4&
   &2 => numerator_t2
   use p9_part21part21_part25part25part21_d54h0l1, only: numerator_diagram54 =>&
   & numerator_ninja
   use p9_part21part21_part25part25part21_d54h0l121, only: numerator_tmu_diagram5&
   &4 => numerator_tmu
   use p9_part21part21_part25part25part21_d54h0l131, only: numerator_t3_diagram5&
   &4 => numerator_t3
   use p9_part21part21_part25part25part21_d54h0l132, only: numerator_t2_diagram5&
   &4 => numerator_t2
   use p9_part21part21_part25part25part21_d109h0l1, only: numerator_diagram109 &
   &=> numerator_ninja
   use p9_part21part21_part25part25part21_d109h0l121, only: numerator_tmu_diagr&
   &am109 => numerator_tmu
   use p9_part21part21_part25part25part21_d109h0l131, only: numerator_t3_diagram1&
   &09 => numerator_t3
   use p9_part21part21_part25part25part21_d109h0l132, only: numerator_t2_diagram1&
   &09 => numerator_t2
   use p9_part21part21_part25part25part21_d141h0l1, only: numerator_diagram141 &
   &=> numerator_ninja
   use p9_part21part21_part25part25part21_d141h0l121, only: numerator_tmu_diagr&
   &am141 => numerator_tmu
   use p9_part21part21_part25part25part21_d141h0l131, only: numerator_t3_diagram1&
   &41 => numerator_t3
   use p9_part21part21_part25part25part21_d141h0l132, only: numerator_t2_diagram1&
   &41 => numerator_t2
   use p9_part21part21_part25part25part21_d147h0l1, only: numerator_diagram147 &
   &=> numerator_ninja
   use p9_part21part21_part25part25part21_d147h0l121, only: numerator_tmu_diagr&
   &am147 => numerator_tmu
   use p9_part21part21_part25part25part21_d147h0l131, only: numerator_t3_diagram1&
   &47 => numerator_t3
   use p9_part21part21_part25part25part21_d147h0l132, only: numerator_t2_diagram1&
   &47 => numerator_t2
   use p9_part21part21_part25part25part21_d190h0l1, only: numerator_diagram190 &
   &=> numerator_ninja
   use p9_part21part21_part25part25part21_d190h0l121, only: numerator_tmu_diagr&
   &am190 => numerator_tmu
   use p9_part21part21_part25part25part21_d190h0l131, only: numerator_t3_diagram1&
   &90 => numerator_t3
   use p9_part21part21_part25part25part21_d190h0l132, only: numerator_t2_diagram1&
   &90 => numerator_t2
   use p9_part21part21_part25part25part21_globalsl1, only: epspow

   implicit none
   real(ki_nin), intent(in) :: scale2
   complex(ki_nin), dimension(-2:0), intent(out) :: tot
   complex(ki_nin), intent(out) :: totr
   logical, intent(out) :: ok

   complex(ki_nin), dimension(-2:0) :: acc
   complex(ki_nin) :: accr
   integer(c_int) :: acc_ok

   integer :: istopm, istop0

   integer, parameter :: effective_group_rank = 5

   !-----------#[ invariants for ninja:
   real(ki_nin), dimension(5,5) :: s_mat
   !-----------#] initialize invariants:
   real(ki_nin), dimension(5) :: msq
   real(ki_nin), dimension(4,5) :: Vi

   ok = .true.

   if (ninja_test.eq.1) then
      istopm = 1
      istop0 = 1
   else
      istopm = ninja_istop
      istop0 = max(2,ninja_istop)
   end if
   msq(1) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,1) = real(-k3-k5-k4, ki_nin)
   msq(2) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,2) = real(-k3-k5, ki_nin)
   msq(3) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,3) = real(-k5, ki_nin)
   msq(4) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,4) = real(0, ki_nin)
   msq(5) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,5) = real(-k2, ki_nin)
   !-----------#[ initialize invariants:
   s_mat(1,1) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(1,2) = real(-2.0_ki*mdlMT**2+mdlMh**2, ki_nin)
   s_mat(2,1) = s_mat(1,2)
   s_mat(1,3) = real(-2.0_ki*mdlMT**2+es34, ki_nin)
   s_mat(3,1) = s_mat(1,3)
   s_mat(1,4) = real(-2.0_ki*mdlMT**2+es12, ki_nin)
   s_mat(4,1) = s_mat(1,4)
   s_mat(1,5) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(5,1) = s_mat(1,5)
   s_mat(2,2) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(2,3) = real(-2.0_ki*mdlMT**2+mdlMh**2, ki_nin)
   s_mat(3,2) = s_mat(2,3)
   s_mat(2,4) = real(-2.0_ki*mdlMT**2+es12+2.0_ki*mdlMh**2-es34-es45, ki_nin)
   s_mat(4,2) = s_mat(2,4)
   s_mat(2,5) = real(-2.0_ki*mdlMT**2+mdlMh**2+es23-es45-es51, ki_nin)
   s_mat(5,2) = s_mat(2,5)
   s_mat(3,3) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(3,4) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(4,3) = s_mat(3,4)
   s_mat(3,5) = real(-2.0_ki*mdlMT**2-es12+es34-es51, ki_nin)
   s_mat(5,3) = s_mat(3,5)
   s_mat(4,4) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(4,5) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(5,4) = s_mat(4,5)
   s_mat(5,5) = real(-2.0_ki*mdlMT**2, ki_nin)
   !-----------#] initialize invariants


      !------#[ sum over reduction of single diagrams:
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='34'>"
         end if
         call ninja_diagram(numerator_diagram34, numerator_t3_diagram34, numera&
         &tor_t2_diagram34, &
          &  numerator_tmu_diagram34, &
          & 5, 3, 3, (/1,3,4/), &
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
            write(logfile,*) "<diagram index='42'>"
         end if
         call ninja_diagram(numerator_diagram42, numerator_t3_diagram42, numera&
         &tor_t2_diagram42, &
          &  numerator_tmu_diagram42, &
          & 5, 3, 3, (/1,3,5/), &
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
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='54'>"
         end if
         call ninja_diagram(numerator_diagram54, numerator_t3_diagram54, numera&
         &tor_t2_diagram54, &
          &  numerator_tmu_diagram54, &
          & 5, 4, 4, (/1,3,4,5/), &
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
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='109'>"
         end if
         call ninja_diagram(numerator_diagram109, numerator_t3_diagram109, nume&
         &rator_t2_diagram109, &
          &  numerator_tmu_diagram109, &
          & 5, 4, 4, (/1,2,3,4/), &
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
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='141'>"
         end if
         call ninja_diagram(numerator_diagram141, numerator_t3_diagram141, nume&
         &rator_t2_diagram141, &
          &  numerator_tmu_diagram141, &
          & 5, 4, 4, (/1,2,3,5/), &
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
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='147'>"
         end if
         call ninja_diagram(numerator_diagram147, numerator_t3_diagram147, nume&
         &rator_t2_diagram147, &
          &  numerator_tmu_diagram147, &
          & 5, 4, 4, (/1,3,4,5/), &
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
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='190'>"
         end if
         call ninja_diagram(numerator_diagram190, numerator_t3_diagram190, nume&
         &rator_t2_diagram190, &
          &  numerator_tmu_diagram190, &
          & 5, 5, 5, (/1,2,3,4,5/), &
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
end subroutine ninja_reduce_group10
!-----#] subroutine ninja_reduce_group10:
!-----#[ subroutine ninja_reduce_group11:
subroutine     ninja_reduce_group11(scale2,tot,totr,ok)
   use iso_c_binding, only: c_ptr, c_loc, c_int
   use ninjago_module
   use p9_part21part21_part25part25part21_kinematics
   use p9_part21part21_part25part25part21_model
   use p9_part21part21_part25part25part21_d113h0l1, only: numerator_diagram113 &
   &=> numerator_ninja
   use p9_part21part21_part25part25part21_d113h0l121, only: numerator_tmu_diagr&
   &am113 => numerator_tmu
   use p9_part21part21_part25part25part21_d113h0l131, only: numerator_t3_diagram1&
   &13 => numerator_t3
   use p9_part21part21_part25part25part21_d113h0l132, only: numerator_t2_diagram1&
   &13 => numerator_t2
   use p9_part21part21_part25part25part21_d186h0l1, only: numerator_diagram186 &
   &=> numerator_ninja
   use p9_part21part21_part25part25part21_d186h0l121, only: numerator_tmu_diagr&
   &am186 => numerator_tmu
   use p9_part21part21_part25part25part21_d186h0l131, only: numerator_t3_diagram1&
   &86 => numerator_t3
   use p9_part21part21_part25part25part21_d186h0l132, only: numerator_t2_diagram1&
   &86 => numerator_t2
   use p9_part21part21_part25part25part21_globalsl1, only: epspow

   implicit none
   real(ki_nin), intent(in) :: scale2
   complex(ki_nin), dimension(-2:0), intent(out) :: tot
   complex(ki_nin), intent(out) :: totr
   logical, intent(out) :: ok

   complex(ki_nin), dimension(-2:0) :: acc
   complex(ki_nin) :: accr
   integer(c_int) :: acc_ok

   integer :: istopm, istop0

   integer, parameter :: effective_group_rank = 5

   !-----------#[ invariants for ninja:
   real(ki_nin), dimension(5,5) :: s_mat
   !-----------#] initialize invariants:
   real(ki_nin), dimension(5) :: msq
   real(ki_nin), dimension(4,5) :: Vi

   ok = .true.

   if (ninja_test.eq.1) then
      istopm = 1
      istop0 = 1
   else
      istopm = ninja_istop
      istop0 = max(2,ninja_istop)
   end if
   msq(1) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,1) = real(-k3-k5-k4, ki_nin)
   msq(2) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,2) = real(-k5-k4, ki_nin)
   msq(3) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,3) = real(-k5, ki_nin)
   msq(4) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,4) = real(0, ki_nin)
   msq(5) = real(mdlMT*mdlMT, ki_nin)
   Vi(:,5) = real(-k2, ki_nin)
   !-----------#[ initialize invariants:
   s_mat(1,1) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(1,2) = real(-2.0_ki*mdlMT**2+mdlMh**2, ki_nin)
   s_mat(2,1) = s_mat(1,2)
   s_mat(1,3) = real(-2.0_ki*mdlMT**2+es34, ki_nin)
   s_mat(3,1) = s_mat(1,3)
   s_mat(1,4) = real(-2.0_ki*mdlMT**2+es12, ki_nin)
   s_mat(4,1) = s_mat(1,4)
   s_mat(1,5) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(5,1) = s_mat(1,5)
   s_mat(2,2) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(2,3) = real(-2.0_ki*mdlMT**2+mdlMh**2, ki_nin)
   s_mat(3,2) = s_mat(2,3)
   s_mat(2,4) = real(-2.0_ki*mdlMT**2+es45, ki_nin)
   s_mat(4,2) = s_mat(2,4)
   s_mat(2,5) = real(-2.0_ki*mdlMT**2+mdlMh**2+es45-es12-es23, ki_nin)
   s_mat(5,2) = s_mat(2,5)
   s_mat(3,3) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(3,4) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(4,3) = s_mat(3,4)
   s_mat(3,5) = real(-2.0_ki*mdlMT**2-es12+es34-es51, ki_nin)
   s_mat(5,3) = s_mat(3,5)
   s_mat(4,4) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(4,5) = real(-2.0_ki*mdlMT**2, ki_nin)
   s_mat(5,4) = s_mat(4,5)
   s_mat(5,5) = real(-2.0_ki*mdlMT**2, ki_nin)
   !-----------#] initialize invariants


      !------#[ sum over reduction of single diagrams:
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='113'>"
         end if
         call ninja_diagram(numerator_diagram113, numerator_t3_diagram113, nume&
         &rator_t2_diagram113, &
          &  numerator_tmu_diagram113, &
          & 5, 4, 4, (/1,2,3,4/), &
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
            write(logfile,*) "<diagram index='186'>"
         end if
         call ninja_diagram(numerator_diagram186, numerator_t3_diagram186, nume&
         &rator_t2_diagram186, &
          &  numerator_tmu_diagram186, &
          & 5, 5, 5, (/1,2,3,4,5/), &
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
end subroutine ninja_reduce_group11
!-----#] subroutine ninja_reduce_group11:
!---#] reduce groups with ninja:
end module p9_part21part21_part25part25part21_ninjah0
