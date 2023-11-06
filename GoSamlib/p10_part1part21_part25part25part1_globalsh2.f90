module     p10_part1part21_part25part25part1_globalsh2
   use p10_part1part21_part25part25part1_config, only: ki
   use p10_part1part21_part25part25part1_color, only:&
      & c1v => c1

   implicit none
   private
   complex(ki), public :: c1

   public :: init_lo

   complex(ki), public :: rat2
contains

subroutine     init_lo()
   use p10_part1part21_part25part25part1_globalsl1, only: epspow, ccontract, amp0
   implicit none
   c1 = ccontract(c1v, amp0)
end subroutine init_lo

end module p10_part1part21_part25part25part1_globalsh2
