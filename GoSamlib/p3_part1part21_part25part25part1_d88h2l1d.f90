module     p3_part1part21_part25part25part1_d88h2l1d
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p3_part1part2 &
   ! &1_part25part25part1/helicity2d88h2l1d.f90
   ! generator: buildfortran_d.py
   use p3_part1part21_part25part25part1_config, only: ki
   use p3_part1part21_part25part25part1_util, only: cond, d => metric_tensor
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   integer, private :: iv0
   integer, private :: iv1
   integer, private :: iv2
   integer, private :: iv3
   integer, private :: iv4
   real(ki), dimension(4), private :: qshift
   public :: derivative
contains
!---#[ function brack_1:
   pure function brack_1(Q, mu2) result(brack)
      use p3_part1part21_part25part25part1_model
      use p3_part1part21_part25part25part1_kinematics
      use p3_part1part21_part25part25part1_color
      use p3_part1part21_part25part25part1_abbrevd88h2
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(1) :: acd88
      complex(ki) :: brack
      acd88(1)=abb88(15)
      brack=acd88(1)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p3_part1part21_part25part25part1_model
      use p3_part1part21_part25part25part1_kinematics
      use p3_part1part21_part25part25part1_color
      use p3_part1part21_part25part25part1_abbrevd88h2
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(36) :: acd88
      complex(ki) :: brack
      acd88(1)=spvak1k2(iv1)
      acd88(2)=abb88(33)
      acd88(3)=spvak1l4(iv1)
      acd88(4)=abb88(27)
      acd88(5)=spvak2k5(iv1)
      acd88(6)=abb88(17)
      acd88(7)=spval4k5(iv1)
      acd88(8)=abb88(8)
      acd88(9)=spvak1e2(iv1)
      acd88(10)=abb88(26)
      acd88(11)=spvak2e2(iv1)
      acd88(12)=abb88(12)
      acd88(13)=spvae2k2(iv1)
      acd88(14)=abb88(37)
      acd88(15)=spval3e2(iv1)
      acd88(16)=abb88(36)
      acd88(17)=spvae2l3(iv1)
      acd88(18)=abb88(19)
      acd88(19)=spval4e2(iv1)
      acd88(20)=abb88(11)
      acd88(21)=spvae2l4(iv1)
      acd88(22)=abb88(29)
      acd88(23)=spvae2k5(iv1)
      acd88(24)=abb88(31)
      acd88(25)=-acd88(2)*acd88(1)
      acd88(26)=-acd88(4)*acd88(3)
      acd88(27)=-acd88(6)*acd88(5)
      acd88(28)=-acd88(8)*acd88(7)
      acd88(29)=-acd88(10)*acd88(9)
      acd88(30)=-acd88(12)*acd88(11)
      acd88(31)=-acd88(14)*acd88(13)
      acd88(32)=-acd88(16)*acd88(15)
      acd88(33)=-acd88(18)*acd88(17)
      acd88(34)=-acd88(20)*acd88(19)
      acd88(35)=-acd88(22)*acd88(21)
      acd88(36)=-acd88(24)*acd88(23)
      brack=acd88(25)+acd88(26)+acd88(27)+acd88(28)+acd88(29)+acd88(30)+acd88(3&
      &1)+acd88(32)+acd88(33)+acd88(34)+acd88(35)+acd88(36)
   end function brack_2
!---#] function brack_2:
!---#[ function brack_3:
   pure function brack_3(Q, mu2) result(brack)
      use p3_part1part21_part25part25part1_model
      use p3_part1part21_part25part25part1_kinematics
      use p3_part1part21_part25part25part1_color
      use p3_part1part21_part25part25part1_abbrevd88h2
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(40) :: acd88
      complex(ki) :: brack
      acd88(1)=d(iv1,iv2)
      acd88(2)=abb88(24)
      acd88(3)=spvak1k2(iv1)
      acd88(4)=spvak2e2(iv2)
      acd88(5)=abb88(32)
      acd88(6)=spvak1k2(iv2)
      acd88(7)=spvak2e2(iv1)
      acd88(8)=spvae2k5(iv2)
      acd88(9)=abb88(13)
      acd88(10)=spvae2k5(iv1)
      acd88(11)=spvak1l4(iv1)
      acd88(12)=spval4e2(iv2)
      acd88(13)=abb88(28)
      acd88(14)=spvak1l4(iv2)
      acd88(15)=spval4e2(iv1)
      acd88(16)=abb88(14)
      acd88(17)=spvak2k5(iv1)
      acd88(18)=spvae2k2(iv2)
      acd88(19)=abb88(22)
      acd88(20)=spvak2k5(iv2)
      acd88(21)=spvae2k2(iv1)
      acd88(22)=spvak1e2(iv2)
      acd88(23)=abb88(16)
      acd88(24)=spvak1e2(iv1)
      acd88(25)=spval4k5(iv1)
      acd88(26)=spvae2l4(iv2)
      acd88(27)=abb88(21)
      acd88(28)=spval4k5(iv2)
      acd88(29)=spvae2l4(iv1)
      acd88(30)=abb88(20)
      acd88(31)=abb88(30)
      acd88(32)=acd88(12)*acd88(16)
      acd88(33)=acd88(4)*acd88(9)
      acd88(34)=acd88(22)*acd88(31)
      acd88(32)=acd88(34)+acd88(32)+acd88(33)
      acd88(32)=acd88(10)*acd88(32)
      acd88(33)=acd88(15)*acd88(16)
      acd88(34)=acd88(7)*acd88(9)
      acd88(35)=acd88(24)*acd88(31)
      acd88(33)=acd88(35)+acd88(33)+acd88(34)
      acd88(33)=acd88(8)*acd88(33)
      acd88(34)=acd88(29)*acd88(28)
      acd88(35)=acd88(26)*acd88(25)
      acd88(34)=acd88(35)+acd88(34)
      acd88(34)=acd88(27)*acd88(34)
      acd88(35)=acd88(21)*acd88(20)
      acd88(36)=acd88(18)*acd88(17)
      acd88(35)=acd88(36)+acd88(35)
      acd88(35)=acd88(19)*acd88(35)
      acd88(36)=acd88(15)*acd88(14)
      acd88(37)=acd88(12)*acd88(11)
      acd88(36)=acd88(37)+acd88(36)
      acd88(36)=acd88(13)*acd88(36)
      acd88(37)=acd88(7)*acd88(6)
      acd88(38)=acd88(4)*acd88(3)
      acd88(37)=acd88(38)+acd88(37)
      acd88(37)=acd88(5)*acd88(37)
      acd88(38)=acd88(26)*acd88(30)
      acd88(39)=acd88(18)*acd88(23)
      acd88(38)=acd88(38)+acd88(39)
      acd88(38)=acd88(24)*acd88(38)
      acd88(39)=acd88(29)*acd88(30)
      acd88(40)=acd88(21)*acd88(23)
      acd88(39)=acd88(39)+acd88(40)
      acd88(39)=acd88(22)*acd88(39)
      acd88(40)=acd88(1)*acd88(2)
      brack=acd88(32)+acd88(33)+acd88(34)+acd88(35)+acd88(36)+acd88(37)+acd88(3&
      &8)+acd88(39)+2.0_ki*acd88(40)
   end function brack_3
!---#] function brack_3:
!---#[ function brack_4:
   pure function brack_4(Q, mu2) result(brack)
      use p3_part1part21_part25part25part1_model
      use p3_part1part21_part25part25part1_kinematics
      use p3_part1part21_part25part25part1_color
      use p3_part1part21_part25part25part1_abbrevd88h2
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(45) :: acd88
      complex(ki) :: brack
      acd88(1)=d(iv1,iv2)
      acd88(2)=spvak1e2(iv3)
      acd88(3)=abb88(23)
      acd88(4)=spvae2k5(iv3)
      acd88(5)=abb88(25)
      acd88(6)=d(iv1,iv3)
      acd88(7)=spvak1e2(iv2)
      acd88(8)=spvae2k5(iv2)
      acd88(9)=d(iv2,iv3)
      acd88(10)=spvak1e2(iv1)
      acd88(11)=spvae2k5(iv1)
      acd88(12)=spvak2k5(iv2)
      acd88(13)=spvae2k2(iv3)
      acd88(14)=abb88(18)
      acd88(15)=spvak2k5(iv3)
      acd88(16)=spvae2k2(iv2)
      acd88(17)=spval4k5(iv2)
      acd88(18)=spvae2l4(iv3)
      acd88(19)=abb88(7)
      acd88(20)=spval4k5(iv3)
      acd88(21)=spvae2l4(iv2)
      acd88(22)=spvak2k5(iv1)
      acd88(23)=spvae2k2(iv1)
      acd88(24)=spval4k5(iv1)
      acd88(25)=spvae2l4(iv1)
      acd88(26)=spvak1k2(iv2)
      acd88(27)=spvak2e2(iv3)
      acd88(28)=spvak1k2(iv3)
      acd88(29)=spvak2e2(iv2)
      acd88(30)=spvak1l4(iv2)
      acd88(31)=spval4e2(iv3)
      acd88(32)=spvak1l4(iv3)
      acd88(33)=spval4e2(iv2)
      acd88(34)=spvak1k2(iv1)
      acd88(35)=spvak2e2(iv1)
      acd88(36)=spvak1l4(iv1)
      acd88(37)=spval4e2(iv1)
      acd88(38)=-acd88(32)*acd88(33)
      acd88(39)=-acd88(30)*acd88(31)
      acd88(38)=acd88(38)+acd88(39)
      acd88(38)=acd88(11)*acd88(38)
      acd88(39)=-acd88(20)*acd88(21)
      acd88(40)=-acd88(17)*acd88(18)
      acd88(39)=acd88(39)+acd88(40)
      acd88(39)=acd88(10)*acd88(39)
      acd88(40)=-acd88(32)*acd88(37)
      acd88(41)=-acd88(31)*acd88(36)
      acd88(40)=acd88(40)+acd88(41)
      acd88(40)=acd88(8)*acd88(40)
      acd88(41)=-acd88(20)*acd88(25)
      acd88(42)=-acd88(18)*acd88(24)
      acd88(41)=acd88(41)+acd88(42)
      acd88(41)=acd88(7)*acd88(41)
      acd88(42)=-acd88(33)*acd88(36)
      acd88(43)=-acd88(30)*acd88(37)
      acd88(42)=acd88(42)+acd88(43)
      acd88(42)=acd88(4)*acd88(42)
      acd88(43)=-acd88(21)*acd88(24)
      acd88(44)=-acd88(17)*acd88(25)
      acd88(43)=acd88(43)+acd88(44)
      acd88(43)=acd88(2)*acd88(43)
      acd88(38)=acd88(43)+acd88(42)+acd88(41)+acd88(40)+acd88(38)+acd88(39)
      acd88(38)=acd88(19)*acd88(38)
      acd88(39)=-acd88(28)*acd88(29)
      acd88(40)=-acd88(26)*acd88(27)
      acd88(39)=acd88(39)+acd88(40)
      acd88(39)=acd88(11)*acd88(39)
      acd88(40)=-acd88(15)*acd88(16)
      acd88(41)=-acd88(12)*acd88(13)
      acd88(40)=acd88(40)+acd88(41)
      acd88(40)=acd88(10)*acd88(40)
      acd88(41)=-acd88(28)*acd88(35)
      acd88(42)=-acd88(27)*acd88(34)
      acd88(41)=acd88(41)+acd88(42)
      acd88(41)=acd88(8)*acd88(41)
      acd88(42)=-acd88(15)*acd88(23)
      acd88(43)=-acd88(13)*acd88(22)
      acd88(42)=acd88(42)+acd88(43)
      acd88(42)=acd88(7)*acd88(42)
      acd88(43)=-acd88(29)*acd88(34)
      acd88(44)=-acd88(26)*acd88(35)
      acd88(43)=acd88(43)+acd88(44)
      acd88(43)=acd88(4)*acd88(43)
      acd88(44)=-acd88(16)*acd88(22)
      acd88(45)=-acd88(12)*acd88(23)
      acd88(44)=acd88(44)+acd88(45)
      acd88(44)=acd88(2)*acd88(44)
      acd88(39)=acd88(44)+acd88(43)+acd88(42)+acd88(41)+acd88(39)+acd88(40)
      acd88(39)=acd88(14)*acd88(39)
      acd88(40)=-acd88(11)*acd88(9)
      acd88(41)=-acd88(8)*acd88(6)
      acd88(42)=-acd88(4)*acd88(1)
      acd88(40)=acd88(42)+acd88(41)+acd88(40)
      acd88(40)=acd88(5)*acd88(40)
      acd88(41)=-acd88(10)*acd88(9)
      acd88(42)=-acd88(7)*acd88(6)
      acd88(43)=-acd88(2)*acd88(1)
      acd88(41)=acd88(43)+acd88(42)+acd88(41)
      acd88(41)=acd88(3)*acd88(41)
      acd88(40)=acd88(41)+acd88(40)
      brack=acd88(38)+acd88(39)+2.0_ki*acd88(40)
   end function brack_4
!---#] function brack_4:
!---#[ function brack_5:
   pure function brack_5(Q, mu2) result(brack)
      use p3_part1part21_part25part25part1_model
      use p3_part1part21_part25part25part1_kinematics
      use p3_part1part21_part25part25part1_color
      use p3_part1part21_part25part25part1_abbrevd88h2
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(21) :: acd88
      complex(ki) :: brack
      acd88(1)=d(iv1,iv2)
      acd88(2)=spvak1e2(iv3)
      acd88(3)=spvae2k5(iv4)
      acd88(4)=abb88(9)
      acd88(5)=spvak1e2(iv4)
      acd88(6)=spvae2k5(iv3)
      acd88(7)=d(iv1,iv3)
      acd88(8)=spvak1e2(iv2)
      acd88(9)=spvae2k5(iv2)
      acd88(10)=d(iv1,iv4)
      acd88(11)=d(iv2,iv3)
      acd88(12)=spvak1e2(iv1)
      acd88(13)=spvae2k5(iv1)
      acd88(14)=d(iv2,iv4)
      acd88(15)=d(iv3,iv4)
      acd88(16)=acd88(13)*acd88(11)
      acd88(17)=acd88(9)*acd88(7)
      acd88(18)=acd88(6)*acd88(1)
      acd88(16)=acd88(18)+acd88(16)+acd88(17)
      acd88(16)=acd88(5)*acd88(16)
      acd88(17)=acd88(13)*acd88(14)
      acd88(18)=acd88(9)*acd88(10)
      acd88(19)=acd88(3)*acd88(1)
      acd88(17)=acd88(19)+acd88(17)+acd88(18)
      acd88(17)=acd88(2)*acd88(17)
      acd88(18)=acd88(9)*acd88(12)
      acd88(19)=acd88(8)*acd88(13)
      acd88(18)=acd88(18)+acd88(19)
      acd88(18)=acd88(15)*acd88(18)
      acd88(19)=acd88(12)*acd88(14)
      acd88(20)=acd88(8)*acd88(10)
      acd88(19)=acd88(19)+acd88(20)
      acd88(19)=acd88(6)*acd88(19)
      acd88(20)=acd88(12)*acd88(11)
      acd88(21)=acd88(8)*acd88(7)
      acd88(20)=acd88(20)+acd88(21)
      acd88(20)=acd88(3)*acd88(20)
      acd88(16)=acd88(17)+acd88(20)+acd88(16)+acd88(19)+acd88(18)
      brack=2.0_ki*acd88(16)*acd88(4)
   end function brack_5
!---#] function brack_5:
!---#[ function derivative:
   function derivative(mu2,i1,i2,i3,i4) result(numerator)
      use p3_part1part21_part25part25part1_globalsl1, only: epspow
      use p3_part1part21_part25part25part1_kinematics
      use p3_part1part21_part25part25part1_abbrevd88h2
      implicit none
      complex(ki), intent(in) :: mu2
      integer, intent(in), optional :: i1
      integer, intent(in), optional :: i2
      integer, intent(in), optional :: i3
      integer, intent(in), optional :: i4
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
      if(present(i4)) then
          iv4=i4
          deg=4
      else
          iv4=1
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
      if(deg.eq.4) then
         numerator = cond(epspow.eq.t1,brack_5,Q,mu2)
         return
      end if
   end function derivative
!---#] function derivative:
end module     p3_part1part21_part25part25part1_d88h2l1d
