module     p0_part21part21_part25part25_globalsh0
   use p0_part21part21_part25part25_config, only: ki
   use p0_part21part21_part25part25_color, only:&
      & c1v => c1

   implicit none
   private
   complex(ki), public :: c1

   public :: init_lo

   complex(ki), public :: rat2
contains

subroutine     init_lo()
   use p0_part21part21_part25part25_globalsl1, only: epspow, col0
   implicit none
   c1 = c1v(col0)
end subroutine init_lo

end module p0_part21part21_part25part25_globalsh0
