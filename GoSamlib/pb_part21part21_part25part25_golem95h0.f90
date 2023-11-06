module     pb_part21part21_part25part25_golem95h0
   use precision_golem, only: ki_gol => ki
   use pb_part21part21_part25part25_config, only: ki
   implicit none
   private
   interface reconstruct_group
      module procedure reconstruct_group0
      module procedure reconstruct_group1
      module procedure reconstruct_group2
   end interface

   public :: reconstruct_group
contains
!---#[ subroutine reconstruct_group0:
subroutine     reconstruct_group0(coeffs)
   use tens_rec
   use pb_part21part21_part25part25_config
   use pb_part21part21_part25part25_groups, only: tensrec_info_group0
   use pb_part21part21_part25part25_d13h0l1, only: numerator_d13 => numerator_g&
   &olem95
   use pb_part21part21_part25part25_d13h0l1d, only: reconstruct_d13
   implicit none
   type(tensrec_info_group0), intent(out) :: coeffs
   !------#[ Diagram 13:
      if (tens_rec_by_derivatives) then
         call reconstruct_d13(coeffs)
      else
         call reconstruct4(numerator_d13, coeffs%coeffs_13)
      end if
   !------#] Diagram 13:
end subroutine reconstruct_group0
!---#] subroutine reconstruct_group0:
!---#[ subroutine reconstruct_group1:
subroutine     reconstruct_group1(coeffs)
   use tens_rec
   use pb_part21part21_part25part25_config
   use pb_part21part21_part25part25_groups, only: tensrec_info_group1
   use pb_part21part21_part25part25_d2h0l1, only: numerator_d2 => numerator_gol&
   &em95
   use pb_part21part21_part25part25_d2h0l1d, only: reconstruct_d2
   use pb_part21part21_part25part25_d7h0l1, only: numerator_d7 => numerator_gol&
   &em95
   use pb_part21part21_part25part25_d7h0l1d, only: reconstruct_d7
   use pb_part21part21_part25part25_d11h0l1, only: numerator_d11 => numerator_g&
   &olem95
   use pb_part21part21_part25part25_d11h0l1d, only: reconstruct_d11
   implicit none
   type(tensrec_info_group1), intent(out) :: coeffs
   !------#[ Diagram 2:
      if (tens_rec_by_derivatives) then
         call reconstruct_d2(coeffs)
      else
         call reconstruct3(numerator_d2, coeffs%coeffs_2)
      end if
   !------#] Diagram 2:
   !------#[ Diagram 7:
      if (tens_rec_by_derivatives) then
         call reconstruct_d7(coeffs)
      else
         call reconstruct3(numerator_d7, coeffs%coeffs_7)
      end if
   !------#] Diagram 7:
   !------#[ Diagram 11:
      if (tens_rec_by_derivatives) then
         call reconstruct_d11(coeffs)
      else
         call reconstruct4(numerator_d11, coeffs%coeffs_11)
      end if
   !------#] Diagram 11:
end subroutine reconstruct_group1
!---#] subroutine reconstruct_group1:
!---#[ subroutine reconstruct_group2:
subroutine     reconstruct_group2(coeffs)
   use tens_rec
   use pb_part21part21_part25part25_config
   use pb_part21part21_part25part25_groups, only: tensrec_info_group2
   use pb_part21part21_part25part25_d9h0l1, only: numerator_d9 => numerator_gol&
   &em95
   use pb_part21part21_part25part25_d9h0l1d, only: reconstruct_d9
   implicit none
   type(tensrec_info_group2), intent(out) :: coeffs
   !------#[ Diagram 9:
      if (tens_rec_by_derivatives) then
         call reconstruct_d9(coeffs)
      else
         call reconstruct4(numerator_d9, coeffs%coeffs_9)
      end if
   !------#] Diagram 9:
end subroutine reconstruct_group2
!---#] subroutine reconstruct_group2:
end module pb_part21part21_part25part25_golem95h0
