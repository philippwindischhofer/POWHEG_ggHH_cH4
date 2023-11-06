module     p3_part1part21_part25part25part1_groups
   use p3_part1part21_part25part25part1_config, only: ki 
   implicit none
   save

   private
   public :: ninja_exit
   public :: quadninja_exit
contains
!---#[ subroutine ninja_exit:
subroutine ninja_exit()
  use ninjago_module, only: ninja_clear_integral_cache
  implicit none
  !------#[ call ninja_clear_integral_cache():
  call ninja_clear_integral_cache()
  !------#] call ninja_clear_integral_cache():
end subroutine ninja_exit
!---#] subroutine ninja_exit:
!---#[ subroutine quadninja_exit:
subroutine quadninja_exit()
  use quadninjago_module, only: quadninja_clear_integral_cache
  implicit none
  !------#[ call quadninja_clear_integral_cache():
  call quadninja_clear_integral_cache()
  !------#] call quadninja_clear_integral_cache():
end subroutine quadninja_exit
!---#] subroutine quadninja_exit:
end module p3_part1part21_part25part25part1_groups
