module     p9_part21part21_part25part25part21_diagramsh0l0_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p9_part21part &
   ! &21_part25part25part21/helicity0diagramsl0_qp.f90
   ! generator: buildfortranborn.py
   use p9_part21part21_part25part25part21_color_qp, only: numcs
   use p9_part21part21_part25part25part21_config, only: ki => ki_qp
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   complex(ki), dimension(numcs), parameter :: zero_col = 0.0_ki
   public :: amplitude
contains
!---#[ function amplitude:
   function amplitude()
      use p9_part21part21_part25part25part21_model_qp
      use p9_part21part21_part25part25part21_kinematics_qp
      use p9_part21part21_part25part25part21_color_qp
      use p9_part21part21_part25part25part21_config, only: debug_lo_diagrams, &
        & use_sorted_sum
      use p9_part21part21_part25part25part21_accu_qp, only: sorted_sum
      use p9_part21part21_part25part25part21_util_qp, only: inspect_lo_diagram
      implicit none
      complex(ki), dimension(numcs) :: amplitude
      complex(ki), dimension(14) :: abb
!      complex(ki), dimension(2,numcs) :: diagrams
      integer :: i
      amplitude(:) = 0.0_ki
      abb(1)=TR**(-1)
      abb(2)=es12**(-1)
      abb(3)=es51**(-1)
      abb(4)=1.0_ki/(-mdlMh**2+es34)
      abb(5)=1.0_ki/(es34-es51-es12)
      abb(6)=abb(2)+abb(3)
      abb(7)=1.0_ki/2.0_ki*es34
      abb(6)=abb(6)*abb(7)
      abb(6)=abb(6)-1.0_ki
      abb(7)=1.0_ki/2.0_ki*abb(2)
      abb(7)=abb(7)*es51
      abb(7)=abb(6)-abb(7)
      abb(8)=spbe2e1*spae1e2
      abb(7)=abb(7)*abb(8)
      abb(9)=spbk5e1*spae1k5
      abb(10)=1.0_ki/2.0_ki*abb(9)
      abb(11)=spbe2k1*spak1e2
      abb(12)=abb(2)*abb(10)*abb(11)
      abb(7)=abb(7)-abb(12)
      abb(12)=spbe5k1*spak1e5
      abb(7)=abb(7)*abb(12)
      abb(13)=es12*abb(5)
      abb(14)=abb(3)*es34
      abb(13)=abb(14)+abb(13)-1.0_ki
      abb(14)=spbe5e2*spae2e5
      abb(10)=abb(10)*abb(14)*abb(13)
      abb(13)=abb(3)*es12
      abb(6)=-abb(6)+1.0_ki/2.0_ki*abb(13)
      abb(11)=abb(11)*spbe5e1*spae1e5
      abb(6)=abb(6)*abb(11)
      abb(6)=abb(7)+abb(10)+abb(6)
      abb(7)=abb(4)*mdlGC43
      abb(10)=mdlGC38*abb(7)*i_
      abb(10)=abb(10)+mdlGC31
      abb(6)=mdlGC6*abb(10)*abb(6)
      abb(9)=abb(14)*abb(9)
      abb(8)=abb(8)*abb(12)
      abb(8)=abb(9)+abb(8)-abb(11)
      abb(9)=i_*mdlGC32
      abb(7)=abb(7)*mdlGC39
      abb(7)=abb(9)-abb(7)
      abb(7)=abb(8)*abb(7)
      abb(6)=abb(6)-1.0_ki/2.0_ki*abb(7)
      abb(7)=1.0_ki/2.0_ki*abb(1)
      abb(6)=abb(6)*abb(7)
      amplitude=c2*abb(6)-c1*abb(6)
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
end module p9_part21part21_part25part25part21_diagramsh0l0_qp
