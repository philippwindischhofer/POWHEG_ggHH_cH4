module     p9_part21part21_part25part25part21_d34h0l1d
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p9_part21part &
   ! &21_part25part25part21/helicity0d34h0l1d.f90
   ! generator: buildfortran_d.py
   use p9_part21part21_part25part25part21_config, only: ki
   use p9_part21part21_part25part25part21_util, only: cond, d => metric_tensor
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   integer, private :: iv0
   integer, private :: iv1
   integer, private :: iv2
   integer, private :: iv3
   real(ki), dimension(4), private :: qshift
   public :: derivative
contains
!---#[ function brack_1:
   pure function brack_1(Q, mu2) result(brack)
      use p9_part21part21_part25part25part21_model
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_color
      use p9_part21part21_part25part25part21_abbrevd34h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(1) :: acd34
      complex(ki) :: brack
      acd34(1)=abb34(4)
      brack=acd34(1)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p9_part21part21_part25part25part21_model
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_color
      use p9_part21part21_part25part25part21_abbrevd34h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(36) :: acd34
      complex(ki) :: brack
      acd34(1)=spvak1k5(iv1)
      acd34(2)=abb34(12)
      acd34(3)=spvak2k5(iv1)
      acd34(4)=abb34(20)
      acd34(5)=spvak5k1(iv1)
      acd34(6)=abb34(16)
      acd34(7)=spvak5k2(iv1)
      acd34(8)=abb34(19)
      acd34(9)=spvak1e5(iv1)
      acd34(10)=abb34(9)
      acd34(11)=spvae5k1(iv1)
      acd34(12)=abb34(13)
      acd34(13)=spvak2e5(iv1)
      acd34(14)=abb34(17)
      acd34(15)=spvae5k2(iv1)
      acd34(16)=abb34(15)
      acd34(17)=spvak5e1(iv1)
      acd34(18)=abb34(11)
      acd34(19)=spvae1k5(iv1)
      acd34(20)=abb34(10)
      acd34(21)=spvae1e5(iv1)
      acd34(22)=abb34(7)
      acd34(23)=spvae5e1(iv1)
      acd34(24)=abb34(6)
      acd34(25)=-acd34(2)*acd34(1)
      acd34(26)=-acd34(4)*acd34(3)
      acd34(27)=-acd34(6)*acd34(5)
      acd34(28)=-acd34(8)*acd34(7)
      acd34(29)=-acd34(10)*acd34(9)
      acd34(30)=-acd34(12)*acd34(11)
      acd34(31)=-acd34(14)*acd34(13)
      acd34(32)=-acd34(16)*acd34(15)
      acd34(33)=-acd34(18)*acd34(17)
      acd34(34)=-acd34(20)*acd34(19)
      acd34(35)=-acd34(22)*acd34(21)
      acd34(36)=-acd34(24)*acd34(23)
      brack=acd34(25)+acd34(26)+acd34(27)+acd34(28)+acd34(29)+acd34(30)+acd34(3&
      &1)+acd34(32)+acd34(33)+acd34(34)+acd34(35)+acd34(36)
   end function brack_2
!---#] function brack_2:
!---#[ function brack_3:
   pure function brack_3(Q, mu2) result(brack)
      use p9_part21part21_part25part25part21_model
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_color
      use p9_part21part21_part25part25part21_abbrevd34h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(20) :: acd34
      complex(ki) :: brack
      acd34(1)=d(iv1,iv2)
      acd34(2)=abb34(5)
      acd34(3)=spvak1e5(iv1)
      acd34(4)=spvae5k1(iv2)
      acd34(5)=abb34(18)
      acd34(6)=spvak1e5(iv2)
      acd34(7)=spvae5k1(iv1)
      acd34(8)=spvak2e5(iv1)
      acd34(9)=spvae5k2(iv2)
      acd34(10)=spvak2e5(iv2)
      acd34(11)=spvae5k2(iv1)
      acd34(12)=spvae1e5(iv1)
      acd34(13)=spvae5e1(iv2)
      acd34(14)=abb34(8)
      acd34(15)=spvae1e5(iv2)
      acd34(16)=spvae5e1(iv1)
      acd34(17)=acd34(4)*acd34(3)
      acd34(18)=acd34(7)*acd34(6)
      acd34(19)=-acd34(9)*acd34(8)
      acd34(20)=-acd34(11)*acd34(10)
      acd34(17)=acd34(20)+acd34(19)+acd34(18)+acd34(17)
      acd34(17)=acd34(5)*acd34(17)
      acd34(18)=acd34(13)*acd34(12)
      acd34(19)=acd34(16)*acd34(15)
      acd34(18)=acd34(19)+acd34(18)
      acd34(18)=acd34(14)*acd34(18)
      acd34(19)=acd34(2)*acd34(1)
      brack=acd34(17)+acd34(18)+2.0_ki*acd34(19)
   end function brack_3
!---#] function brack_3:
!---#[ function brack_4:
   pure function brack_4(Q, mu2) result(brack)
      use p9_part21part21_part25part25part21_model
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_color
      use p9_part21part21_part25part25part21_abbrevd34h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(1) :: acd34
      complex(ki) :: brack
      brack=0.0_ki
   end function brack_4
!---#] function brack_4:
!---#[ function derivative:
   function derivative(mu2,i1,i2,i3) result(numerator)
      use p9_part21part21_part25part25part21_globalsl1, only: epspow
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_abbrevd34h0
      implicit none
      complex(ki), intent(in) :: mu2
      integer, intent(in), optional :: i1
      integer, intent(in), optional :: i2
      integer, intent(in), optional :: i3
      complex(ki) :: numerator
      complex(ki) :: loc
      integer :: t1
      integer :: deg
      complex(ki), dimension(4), parameter :: Q = (/ (0.0_ki,0.0_ki),(0.0_ki,0.&
      &0_ki),(0.0_ki,0.0_ki),(0.0_ki,0.0_ki)/)
      qshift = 0
      numerator = 0.0_ki
      deg = 0
      if(present(i1)) then
          iv1=i1
          deg=1
      else
          iv1=1
      end if
      if(present(i2)) then
          iv2=i2
          deg=2
      else
          iv2=1
      end if
      if(present(i3)) then
          iv3=i3
          deg=3
      else
          iv3=1
      end if
      t1 = 0
      if(deg.eq.0) then
         numerator = cond(epspow.eq.t1,brack_1,Q,mu2)
         return
      end if
      if(deg.eq.1) then
         numerator = cond(epspow.eq.t1,brack_2,Q,mu2)
         return
      end if
      if(deg.eq.2) then
         numerator = cond(epspow.eq.t1,brack_3,Q,mu2)
         return
      end if
      if(deg.eq.3) then
         numerator = cond(epspow.eq.t1,brack_4,Q,mu2)
         return
      end if
   end function derivative
!---#] function derivative:
end module     p9_part21part21_part25part25part21_d34h0l1d
