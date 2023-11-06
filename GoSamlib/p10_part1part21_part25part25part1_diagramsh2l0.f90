module     p10_part1part21_part25part25part1_diagramsh2l0
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p10_part1part &
   ! &21_part25part25part1/helicity2diagramsl0.f90
   ! generator: buildfortranborn.py
   use p10_part1part21_part25part25part1_color, only: numcs
   use p10_part1part21_part25part25part1_config, only: ki
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   complex(ki), dimension(numcs), parameter :: zero_col = 0.0_ki
   public :: amplitude
contains
!---#[ function amplitude:
   function amplitude()
      use p10_part1part21_part25part25part1_model
      use p10_part1part21_part25part25part1_kinematics
      use p10_part1part21_part25part25part1_color
      use p10_part1part21_part25part25part1_config, only: debug_lo_diagrams, &
        & use_sorted_sum
      use p10_part1part21_part25part25part1_accu, only: sorted_sum
      use p10_part1part21_part25part25part1_util, only: inspect_lo_diagram
      implicit none
      complex(ki), dimension(numcs) :: amplitude
      complex(ki), dimension(5) :: abb
!      complex(ki), dimension(2,numcs) :: diagrams
      integer :: i
      amplitude(:) = 0.0_ki
      abb(1)=es51**(-1)
      abb(2)=1.0_ki/(-mdlMh**2+es34)
      abb(3)=mdlGC38*abb(2)*mdlGC43
      abb(4)=mdlGC31*i_
      abb(3)=abb(3)-abb(4)
      abb(4)=-spak1k2*abb(3)*spae2k5*spbk5k2
      abb(3)=abb(3)*spak1e2
      abb(5)=-es34*abb(3)
      abb(4)=abb(5)+abb(4)
      abb(4)=abb(1)*abb(4)
      abb(3)=abb(3)+abb(4)
      abb(3)=1.0_ki/2.0_ki*mdlGC7*spbk5e2*abb(3)
      amplitude=c1*abb(3)
      if (debug_lo_diagrams) then
         write(*,*) "Using Born optimization, debug_lo_diagrams not implemented&
         &."
      end if
!      if (use_sorted_sum) then
!         do i=1,numcs
!            amplitude(i) = sorted_sum(diagrams(i))
!         end do
!      else
!         do i=1,numcs
!            amplitude(i) = sum(diagrams(i))
!         end do
!      end if
   end function     amplitude
!---#] function amplitude:
end module p10_part1part21_part25part25part1_diagramsh2l0
