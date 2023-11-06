module     p2_part21part21_part25part25part21_globalsh0
   use p2_part21part21_part25part25part21_config, only: ki
   use p2_part21part21_part25part25part21_color, only:&
      & c1v => c1,&
      & c2v => c2

   implicit none
   private
   complex(ki), public :: c1
   complex(ki), public :: c2

   public :: init_lo

   complex(ki), public :: rat2
contains

subroutine     init_lo()
   use p2_part21part21_part25part25part21_globalsl1, only: epspow, col0
   implicit none
   c1 = c1v(col0)
   c2 = c2v(col0)
end subroutine init_lo

end module p2_part21part21_part25part25part21_globalsh0
