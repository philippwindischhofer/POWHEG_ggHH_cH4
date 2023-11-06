module     pb_part21part21_part25part25_d11h0l1d
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/pb_part21part &
   ! &21_part25part25/helicity0d11h0l1d.f90
   ! generator: buildfortran_d.py
   use pb_part21part21_part25part25_config, only: ki
   use pb_part21part21_part25part25_util, only: cond, d => metric_tensor
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   integer, private :: iv0
   integer, private :: iv1
   integer, private :: iv2
   integer, private :: iv3
   integer, private :: iv4
   real(ki), dimension(4), private :: qshift
   public :: derivative , reconstruct_d11
contains
!---#[ function brack_1:
   pure function brack_1(Q, mu2) result(brack)
      use pb_part21part21_part25part25_model
      use pb_part21part21_part25part25_kinematics
      use pb_part21part21_part25part25_color
      use pb_part21part21_part25part25_abbrevd11h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(64) :: acd11
      complex(ki) :: brack
      acd11(1)=dotproduct(qshift,qshift)
      acd11(2)=dotproduct(qshift,spvae1e2)
      acd11(3)=dotproduct(qshift,spvae2e1)
      acd11(4)=abb11(24)
      acd11(5)=abb11(23)
      acd11(6)=abb11(22)
      acd11(7)=abb11(7)
      acd11(8)=abb11(9)
      acd11(9)=dotproduct(qshift,spvak1k2)
      acd11(10)=abb11(27)
      acd11(11)=dotproduct(qshift,spvak1e1)
      acd11(12)=dotproduct(qshift,spvae2k1)
      acd11(13)=abb11(4)
      acd11(14)=abb11(36)
      acd11(15)=dotproduct(qshift,spvak2e1)
      acd11(16)=dotproduct(qshift,spvae2k2)
      acd11(17)=abb11(14)
      acd11(18)=abb11(21)
      acd11(19)=dotproduct(qshift,spval4e1)
      acd11(20)=dotproduct(qshift,spvae2l4)
      acd11(21)=abb11(28)
      acd11(22)=abb11(8)
      acd11(23)=abb11(11)
      acd11(24)=dotproduct(qshift,spvak2k1)
      acd11(25)=abb11(25)
      acd11(26)=dotproduct(qshift,spvae1k1)
      acd11(27)=dotproduct(qshift,spvak1e2)
      acd11(28)=abb11(35)
      acd11(29)=dotproduct(qshift,spvae1k2)
      acd11(30)=dotproduct(qshift,spvak2e2)
      acd11(31)=abb11(17)
      acd11(32)=dotproduct(qshift,spvae1l4)
      acd11(33)=dotproduct(qshift,spval4e2)
      acd11(34)=abb11(26)
      acd11(35)=abb11(12)
      acd11(36)=abb11(18)
      acd11(37)=abb11(10)
      acd11(38)=abb11(33)
      acd11(39)=abb11(13)
      acd11(40)=abb11(34)
      acd11(41)=abb11(20)
      acd11(42)=abb11(15)
      acd11(43)=abb11(16)
      acd11(44)=abb11(19)
      acd11(45)=abb11(32)
      acd11(46)=abb11(29)
      acd11(47)=abb11(31)
      acd11(48)=abb11(30)
      acd11(49)=abb11(5)
      acd11(50)=acd11(19)*acd11(20)
      acd11(51)=-acd11(11)*acd11(12)
      acd11(50)=acd11(50)+acd11(51)
      acd11(50)=acd11(13)*acd11(50)
      acd11(51)=acd11(9)*acd11(10)
      acd11(52)=acd11(20)*acd11(22)
      acd11(53)=acd11(19)*acd11(21)
      acd11(54)=acd11(15)*acd11(17)
      acd11(54)=acd11(18)+acd11(54)
      acd11(54)=acd11(16)*acd11(54)
      acd11(55)=acd11(11)*acd11(14)
      acd11(56)=-acd11(1)*acd11(5)
      acd11(57)=acd11(1)*acd11(4)
      acd11(57)=acd11(8)+acd11(57)
      acd11(57)=acd11(3)*acd11(57)
      acd11(50)=acd11(57)+acd11(56)+acd11(50)+acd11(55)+acd11(54)+acd11(53)+acd&
      &11(52)-acd11(23)+acd11(51)
      acd11(50)=acd11(2)*acd11(50)
      acd11(51)=acd11(32)*acd11(33)
      acd11(52)=-acd11(26)*acd11(27)
      acd11(51)=acd11(51)+acd11(52)
      acd11(51)=acd11(13)*acd11(51)
      acd11(52)=acd11(24)*acd11(25)
      acd11(53)=acd11(33)*acd11(35)
      acd11(54)=acd11(32)*acd11(34)
      acd11(55)=acd11(17)*acd11(29)
      acd11(55)=acd11(31)+acd11(55)
      acd11(55)=acd11(30)*acd11(55)
      acd11(56)=acd11(26)*acd11(28)
      acd11(57)=-acd11(1)*acd11(6)
      acd11(51)=acd11(57)+acd11(51)+acd11(56)+acd11(55)+acd11(54)+acd11(53)-acd&
      &11(36)+acd11(52)
      acd11(51)=acd11(3)*acd11(51)
      acd11(52)=-acd11(29)*acd11(43)
      acd11(53)=-acd11(27)*acd11(40)
      acd11(54)=-acd11(15)*acd11(41)
      acd11(55)=-acd11(12)*acd11(38)
      acd11(56)=-acd11(33)*acd11(48)
      acd11(57)=-acd11(32)*acd11(47)
      acd11(58)=-acd11(30)*acd11(44)
      acd11(59)=-acd11(26)*acd11(39)
      acd11(60)=-acd11(20)*acd11(46)
      acd11(61)=-acd11(19)*acd11(45)
      acd11(62)=-acd11(16)*acd11(42)
      acd11(63)=-acd11(11)*acd11(37)
      acd11(64)=acd11(1)*acd11(7)
      brack=acd11(49)+acd11(50)+acd11(51)+acd11(52)+acd11(53)+acd11(54)+acd11(5&
      &5)+acd11(56)+acd11(57)+acd11(58)+acd11(59)+acd11(60)+acd11(61)+acd11(62)&
      &+acd11(63)+acd11(64)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use pb_part21part21_part25part25_model
      use pb_part21part21_part25part25_kinematics
      use pb_part21part21_part25part25_color
      use pb_part21part21_part25part25_abbrevd11h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(82) :: acd11
      complex(ki) :: brack
      acd11(1)=qshift(iv1)
      acd11(2)=dotproduct(qshift,spvae1e2)
      acd11(3)=dotproduct(qshift,spvae2e1)
      acd11(4)=abb11(24)
      acd11(5)=abb11(23)
      acd11(6)=abb11(22)
      acd11(7)=abb11(7)
      acd11(8)=spvae1e2(iv1)
      acd11(9)=dotproduct(qshift,qshift)
      acd11(10)=abb11(9)
      acd11(11)=dotproduct(qshift,spvak1k2)
      acd11(12)=abb11(27)
      acd11(13)=dotproduct(qshift,spvak1e1)
      acd11(14)=dotproduct(qshift,spvae2k1)
      acd11(15)=abb11(4)
      acd11(16)=abb11(36)
      acd11(17)=dotproduct(qshift,spvak2e1)
      acd11(18)=dotproduct(qshift,spvae2k2)
      acd11(19)=abb11(14)
      acd11(20)=abb11(21)
      acd11(21)=dotproduct(qshift,spval4e1)
      acd11(22)=dotproduct(qshift,spvae2l4)
      acd11(23)=abb11(28)
      acd11(24)=abb11(8)
      acd11(25)=abb11(11)
      acd11(26)=spvae2e1(iv1)
      acd11(27)=dotproduct(qshift,spvak2k1)
      acd11(28)=abb11(25)
      acd11(29)=dotproduct(qshift,spvae1k1)
      acd11(30)=dotproduct(qshift,spvak1e2)
      acd11(31)=abb11(35)
      acd11(32)=dotproduct(qshift,spvae1k2)
      acd11(33)=dotproduct(qshift,spvak2e2)
      acd11(34)=abb11(17)
      acd11(35)=dotproduct(qshift,spvae1l4)
      acd11(36)=dotproduct(qshift,spval4e2)
      acd11(37)=abb11(26)
      acd11(38)=abb11(12)
      acd11(39)=abb11(18)
      acd11(40)=spvak1k2(iv1)
      acd11(41)=spvak2k1(iv1)
      acd11(42)=spvak1e1(iv1)
      acd11(43)=abb11(10)
      acd11(44)=spvae2k1(iv1)
      acd11(45)=abb11(33)
      acd11(46)=spvae1k1(iv1)
      acd11(47)=abb11(13)
      acd11(48)=spvak1e2(iv1)
      acd11(49)=abb11(34)
      acd11(50)=spvak2e1(iv1)
      acd11(51)=abb11(20)
      acd11(52)=spvae2k2(iv1)
      acd11(53)=abb11(15)
      acd11(54)=spvae1k2(iv1)
      acd11(55)=abb11(16)
      acd11(56)=spvak2e2(iv1)
      acd11(57)=abb11(19)
      acd11(58)=spval4e1(iv1)
      acd11(59)=abb11(32)
      acd11(60)=spvae2l4(iv1)
      acd11(61)=abb11(29)
      acd11(62)=spvae1l4(iv1)
      acd11(63)=abb11(31)
      acd11(64)=spval4e2(iv1)
      acd11(65)=abb11(30)
      acd11(66)=acd11(46)*acd11(30)
      acd11(67)=-acd11(36)*acd11(62)
      acd11(68)=-acd11(35)*acd11(64)
      acd11(69)=acd11(29)*acd11(48)
      acd11(66)=acd11(69)+acd11(68)+acd11(66)+acd11(67)
      acd11(66)=acd11(15)*acd11(66)
      acd11(67)=-acd11(56)*acd11(32)
      acd11(68)=-acd11(33)*acd11(54)
      acd11(67)=acd11(67)+acd11(68)
      acd11(67)=acd11(19)*acd11(67)
      acd11(68)=-acd11(28)*acd11(41)
      acd11(69)=-acd11(64)*acd11(38)
      acd11(70)=-acd11(62)*acd11(37)
      acd11(71)=-acd11(56)*acd11(34)
      acd11(72)=-acd11(46)*acd11(31)
      acd11(73)=2.0_ki*acd11(1)
      acd11(74)=acd11(6)*acd11(73)
      acd11(75)=acd11(9)*acd11(4)
      acd11(75)=acd11(75)+acd11(10)
      acd11(76)=-acd11(8)*acd11(75)
      acd11(66)=acd11(76)+acd11(66)+acd11(67)+acd11(74)+acd11(72)+acd11(71)+acd&
      &11(70)+acd11(68)+acd11(69)
      acd11(66)=acd11(3)*acd11(66)
      acd11(67)=acd11(42)*acd11(14)
      acd11(68)=-acd11(22)*acd11(58)
      acd11(69)=-acd11(21)*acd11(60)
      acd11(70)=acd11(13)*acd11(44)
      acd11(67)=acd11(70)+acd11(69)+acd11(67)+acd11(68)
      acd11(67)=acd11(15)*acd11(67)
      acd11(68)=-acd11(3)*acd11(4)
      acd11(68)=acd11(68)+acd11(5)
      acd11(68)=acd11(73)*acd11(68)
      acd11(69)=-acd11(52)*acd11(17)
      acd11(70)=-acd11(18)*acd11(50)
      acd11(69)=acd11(69)+acd11(70)
      acd11(69)=acd11(19)*acd11(69)
      acd11(70)=-acd11(12)*acd11(40)
      acd11(71)=-acd11(60)*acd11(24)
      acd11(72)=-acd11(58)*acd11(23)
      acd11(74)=-acd11(52)*acd11(20)
      acd11(76)=-acd11(42)*acd11(16)
      acd11(75)=-acd11(26)*acd11(75)
      acd11(67)=acd11(67)+acd11(75)+acd11(69)+acd11(76)+acd11(74)+acd11(72)+acd&
      &11(70)+acd11(71)+acd11(68)
      acd11(67)=acd11(2)*acd11(67)
      acd11(68)=-acd11(35)*acd11(36)
      acd11(69)=acd11(29)*acd11(30)
      acd11(68)=acd11(68)+acd11(69)
      acd11(68)=acd11(15)*acd11(68)
      acd11(69)=-acd11(19)*acd11(32)
      acd11(69)=acd11(69)-acd11(34)
      acd11(69)=acd11(33)*acd11(69)
      acd11(70)=-acd11(28)*acd11(27)
      acd11(71)=-acd11(36)*acd11(38)
      acd11(72)=-acd11(35)*acd11(37)
      acd11(74)=-acd11(29)*acd11(31)
      acd11(75)=acd11(9)*acd11(6)
      acd11(68)=acd11(68)+acd11(75)+acd11(74)+acd11(72)+acd11(71)+acd11(39)+acd&
      &11(70)+acd11(69)
      acd11(68)=acd11(26)*acd11(68)
      acd11(69)=-acd11(21)*acd11(22)
      acd11(70)=acd11(13)*acd11(14)
      acd11(69)=acd11(69)+acd11(70)
      acd11(69)=acd11(15)*acd11(69)
      acd11(70)=-acd11(19)*acd11(17)
      acd11(70)=acd11(70)-acd11(20)
      acd11(70)=acd11(18)*acd11(70)
      acd11(71)=-acd11(12)*acd11(11)
      acd11(72)=-acd11(22)*acd11(24)
      acd11(74)=-acd11(21)*acd11(23)
      acd11(75)=-acd11(13)*acd11(16)
      acd11(76)=acd11(9)*acd11(5)
      acd11(69)=acd11(69)+acd11(76)+acd11(75)+acd11(74)+acd11(72)+acd11(25)+acd&
      &11(71)+acd11(70)
      acd11(69)=acd11(8)*acd11(69)
      acd11(70)=acd11(54)*acd11(55)
      acd11(71)=acd11(50)*acd11(51)
      acd11(72)=acd11(48)*acd11(49)
      acd11(74)=acd11(44)*acd11(45)
      acd11(75)=acd11(64)*acd11(65)
      acd11(76)=acd11(62)*acd11(63)
      acd11(77)=acd11(60)*acd11(61)
      acd11(78)=acd11(58)*acd11(59)
      acd11(79)=acd11(56)*acd11(57)
      acd11(80)=acd11(52)*acd11(53)
      acd11(81)=acd11(46)*acd11(47)
      acd11(82)=acd11(42)*acd11(43)
      acd11(73)=-acd11(7)*acd11(73)
      brack=acd11(66)+acd11(67)+acd11(68)+acd11(69)+acd11(70)+acd11(71)+acd11(7&
      &2)+acd11(73)+acd11(74)+acd11(75)+acd11(76)+acd11(77)+acd11(78)+acd11(79)&
      &+acd11(80)+acd11(81)+acd11(82)
   end function brack_2
!---#] function brack_2:
!---#[ function brack_3:
   pure function brack_3(Q, mu2) result(brack)
      use pb_part21part21_part25part25_model
      use pb_part21part21_part25part25_kinematics
      use pb_part21part21_part25part25_color
      use pb_part21part21_part25part25_abbrevd11h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(80) :: acd11
      complex(ki) :: brack
      acd11(1)=d(iv1,iv2)
      acd11(2)=dotproduct(qshift,spvae1e2)
      acd11(3)=dotproduct(qshift,spvae2e1)
      acd11(4)=abb11(24)
      acd11(5)=abb11(23)
      acd11(6)=abb11(22)
      acd11(7)=abb11(7)
      acd11(8)=qshift(iv1)
      acd11(9)=spvae1e2(iv2)
      acd11(10)=spvae2e1(iv2)
      acd11(11)=qshift(iv2)
      acd11(12)=spvae1e2(iv1)
      acd11(13)=spvae2e1(iv1)
      acd11(14)=dotproduct(qshift,qshift)
      acd11(15)=abb11(9)
      acd11(16)=spvak1k2(iv2)
      acd11(17)=abb11(27)
      acd11(18)=spvak1e1(iv2)
      acd11(19)=dotproduct(qshift,spvae2k1)
      acd11(20)=abb11(4)
      acd11(21)=abb11(36)
      acd11(22)=spvae2k1(iv2)
      acd11(23)=dotproduct(qshift,spvak1e1)
      acd11(24)=spvak2e1(iv2)
      acd11(25)=dotproduct(qshift,spvae2k2)
      acd11(26)=abb11(14)
      acd11(27)=spvae2k2(iv2)
      acd11(28)=dotproduct(qshift,spvak2e1)
      acd11(29)=abb11(21)
      acd11(30)=spval4e1(iv2)
      acd11(31)=dotproduct(qshift,spvae2l4)
      acd11(32)=abb11(28)
      acd11(33)=spvae2l4(iv2)
      acd11(34)=dotproduct(qshift,spval4e1)
      acd11(35)=abb11(8)
      acd11(36)=spvak1k2(iv1)
      acd11(37)=spvak1e1(iv1)
      acd11(38)=spvae2k1(iv1)
      acd11(39)=spvak2e1(iv1)
      acd11(40)=spvae2k2(iv1)
      acd11(41)=spval4e1(iv1)
      acd11(42)=spvae2l4(iv1)
      acd11(43)=spvak2k1(iv2)
      acd11(44)=abb11(25)
      acd11(45)=spvae1k1(iv2)
      acd11(46)=dotproduct(qshift,spvak1e2)
      acd11(47)=abb11(35)
      acd11(48)=spvak1e2(iv2)
      acd11(49)=dotproduct(qshift,spvae1k1)
      acd11(50)=spvae1k2(iv2)
      acd11(51)=dotproduct(qshift,spvak2e2)
      acd11(52)=spvak2e2(iv2)
      acd11(53)=dotproduct(qshift,spvae1k2)
      acd11(54)=abb11(17)
      acd11(55)=spvae1l4(iv2)
      acd11(56)=dotproduct(qshift,spval4e2)
      acd11(57)=abb11(26)
      acd11(58)=spval4e2(iv2)
      acd11(59)=dotproduct(qshift,spvae1l4)
      acd11(60)=abb11(12)
      acd11(61)=spvak2k1(iv1)
      acd11(62)=spvae1k1(iv1)
      acd11(63)=spvak1e2(iv1)
      acd11(64)=spvae1k2(iv1)
      acd11(65)=spvak2e2(iv1)
      acd11(66)=spvae1l4(iv1)
      acd11(67)=spval4e2(iv1)
      acd11(68)=acd11(51)*acd11(64)
      acd11(69)=acd11(65)*acd11(53)
      acd11(68)=acd11(68)+acd11(69)
      acd11(68)=acd11(26)*acd11(68)
      acd11(69)=acd11(44)*acd11(61)
      acd11(70)=acd11(67)*acd11(60)
      acd11(71)=acd11(66)*acd11(57)
      acd11(72)=acd11(65)*acd11(54)
      acd11(73)=acd11(62)*acd11(47)
      acd11(74)=2.0_ki*acd11(8)
      acd11(75)=-acd11(6)*acd11(74)
      acd11(76)=2.0_ki*acd11(4)
      acd11(77)=acd11(76)*acd11(2)
      acd11(78)=acd11(8)*acd11(77)
      acd11(79)=acd11(4)*acd11(14)
      acd11(79)=acd11(79)+acd11(15)
      acd11(80)=acd11(12)*acd11(79)
      acd11(68)=acd11(80)+acd11(68)+acd11(78)+acd11(75)+acd11(73)+acd11(72)+acd&
      &11(71)+acd11(69)+acd11(70)
      acd11(68)=acd11(10)*acd11(68)
      acd11(69)=acd11(25)*acd11(39)
      acd11(70)=acd11(40)*acd11(28)
      acd11(69)=acd11(69)+acd11(70)
      acd11(69)=acd11(26)*acd11(69)
      acd11(70)=acd11(17)*acd11(36)
      acd11(71)=acd11(42)*acd11(35)
      acd11(72)=acd11(41)*acd11(32)
      acd11(73)=acd11(40)*acd11(29)
      acd11(75)=acd11(37)*acd11(21)
      acd11(74)=-acd11(5)*acd11(74)
      acd11(76)=acd11(76)*acd11(3)
      acd11(78)=acd11(8)*acd11(76)
      acd11(79)=acd11(13)*acd11(79)
      acd11(69)=acd11(79)+acd11(69)+acd11(78)+acd11(74)+acd11(75)+acd11(73)+acd&
      &11(72)+acd11(70)+acd11(71)
      acd11(69)=acd11(9)*acd11(69)
      acd11(70)=-acd11(62)*acd11(48)
      acd11(71)=acd11(58)*acd11(66)
      acd11(72)=acd11(55)*acd11(67)
      acd11(73)=-acd11(45)*acd11(63)
      acd11(70)=acd11(73)+acd11(72)+acd11(70)+acd11(71)
      acd11(70)=acd11(3)*acd11(70)
      acd11(71)=-acd11(37)*acd11(22)
      acd11(72)=acd11(33)*acd11(41)
      acd11(73)=acd11(30)*acd11(42)
      acd11(74)=-acd11(18)*acd11(38)
      acd11(71)=acd11(74)+acd11(73)+acd11(71)+acd11(72)
      acd11(71)=acd11(2)*acd11(71)
      acd11(72)=-acd11(48)*acd11(49)
      acd11(73)=acd11(58)*acd11(59)
      acd11(74)=acd11(55)*acd11(56)
      acd11(75)=-acd11(45)*acd11(46)
      acd11(72)=acd11(75)+acd11(74)+acd11(72)+acd11(73)
      acd11(72)=acd11(13)*acd11(72)
      acd11(73)=-acd11(22)*acd11(23)
      acd11(74)=acd11(33)*acd11(34)
      acd11(75)=acd11(30)*acd11(31)
      acd11(78)=-acd11(18)*acd11(19)
      acd11(73)=acd11(78)+acd11(75)+acd11(73)+acd11(74)
      acd11(73)=acd11(12)*acd11(73)
      acd11(74)=-acd11(49)*acd11(63)
      acd11(75)=acd11(67)*acd11(59)
      acd11(78)=acd11(66)*acd11(56)
      acd11(79)=-acd11(62)*acd11(46)
      acd11(74)=acd11(79)+acd11(78)+acd11(74)+acd11(75)
      acd11(74)=acd11(10)*acd11(74)
      acd11(75)=-acd11(23)*acd11(38)
      acd11(78)=acd11(42)*acd11(34)
      acd11(79)=acd11(41)*acd11(31)
      acd11(80)=-acd11(37)*acd11(19)
      acd11(75)=acd11(80)+acd11(79)+acd11(75)+acd11(78)
      acd11(75)=acd11(9)*acd11(75)
      acd11(70)=acd11(75)+acd11(74)+acd11(73)+acd11(72)+acd11(70)+acd11(71)
      acd11(70)=acd11(20)*acd11(70)
      acd11(71)=acd11(50)*acd11(51)
      acd11(72)=acd11(52)*acd11(53)
      acd11(71)=acd11(71)+acd11(72)
      acd11(71)=acd11(26)*acd11(71)
      acd11(72)=acd11(44)*acd11(43)
      acd11(73)=acd11(58)*acd11(60)
      acd11(74)=acd11(55)*acd11(57)
      acd11(75)=acd11(52)*acd11(54)
      acd11(78)=acd11(45)*acd11(47)
      acd11(79)=2.0_ki*acd11(11)
      acd11(80)=-acd11(6)*acd11(79)
      acd11(77)=acd11(11)*acd11(77)
      acd11(71)=acd11(71)+acd11(77)+acd11(80)+acd11(78)+acd11(75)+acd11(74)+acd&
      &11(72)+acd11(73)
      acd11(71)=acd11(13)*acd11(71)
      acd11(72)=acd11(24)*acd11(25)
      acd11(73)=acd11(27)*acd11(28)
      acd11(72)=acd11(72)+acd11(73)
      acd11(72)=acd11(26)*acd11(72)
      acd11(73)=acd11(17)*acd11(16)
      acd11(74)=acd11(33)*acd11(35)
      acd11(75)=acd11(30)*acd11(32)
      acd11(77)=acd11(27)*acd11(29)
      acd11(78)=acd11(18)*acd11(21)
      acd11(79)=-acd11(5)*acd11(79)
      acd11(76)=acd11(11)*acd11(76)
      acd11(72)=acd11(72)+acd11(76)+acd11(79)+acd11(78)+acd11(77)+acd11(75)+acd&
      &11(73)+acd11(74)
      acd11(72)=acd11(12)*acd11(72)
      acd11(73)=acd11(65)*acd11(50)
      acd11(74)=acd11(52)*acd11(64)
      acd11(73)=acd11(73)+acd11(74)
      acd11(73)=acd11(3)*acd11(73)
      acd11(74)=acd11(40)*acd11(24)
      acd11(75)=acd11(27)*acd11(39)
      acd11(74)=acd11(74)+acd11(75)
      acd11(74)=acd11(2)*acd11(74)
      acd11(73)=acd11(73)+acd11(74)
      acd11(73)=acd11(26)*acd11(73)
      acd11(74)=-acd11(6)*acd11(3)
      acd11(75)=acd11(3)*acd11(4)
      acd11(75)=-acd11(5)+acd11(75)
      acd11(75)=acd11(2)*acd11(75)
      acd11(74)=acd11(75)+acd11(7)+acd11(74)
      acd11(74)=acd11(1)*acd11(74)
      brack=acd11(68)+acd11(69)+acd11(70)+acd11(71)+acd11(72)+acd11(73)+2.0_ki*&
      &acd11(74)
   end function brack_3
!---#] function brack_3:
!---#[ function brack_4:
   pure function brack_4(Q, mu2) result(brack)
      use pb_part21part21_part25part25_model
      use pb_part21part21_part25part25_kinematics
      use pb_part21part21_part25part25_color
      use pb_part21part21_part25part25_abbrevd11h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(65) :: acd11
      complex(ki) :: brack
      acd11(1)=d(iv1,iv2)
      acd11(2)=spvae1e2(iv3)
      acd11(3)=dotproduct(qshift,spvae2e1)
      acd11(4)=abb11(24)
      acd11(5)=abb11(23)
      acd11(6)=spvae2e1(iv3)
      acd11(7)=dotproduct(qshift,spvae1e2)
      acd11(8)=abb11(22)
      acd11(9)=d(iv1,iv3)
      acd11(10)=spvae1e2(iv2)
      acd11(11)=spvae2e1(iv2)
      acd11(12)=d(iv2,iv3)
      acd11(13)=spvae1e2(iv1)
      acd11(14)=spvae2e1(iv1)
      acd11(15)=qshift(iv1)
      acd11(16)=qshift(iv2)
      acd11(17)=qshift(iv3)
      acd11(18)=spvak1e1(iv2)
      acd11(19)=spvae2k1(iv3)
      acd11(20)=abb11(4)
      acd11(21)=spvak1e1(iv3)
      acd11(22)=spvae2k1(iv2)
      acd11(23)=spvak2e1(iv2)
      acd11(24)=spvae2k2(iv3)
      acd11(25)=abb11(14)
      acd11(26)=spvak2e1(iv3)
      acd11(27)=spvae2k2(iv2)
      acd11(28)=spval4e1(iv2)
      acd11(29)=spvae2l4(iv3)
      acd11(30)=spval4e1(iv3)
      acd11(31)=spvae2l4(iv2)
      acd11(32)=spvak1e1(iv1)
      acd11(33)=spvae2k1(iv1)
      acd11(34)=spvak2e1(iv1)
      acd11(35)=spvae2k2(iv1)
      acd11(36)=spval4e1(iv1)
      acd11(37)=spvae2l4(iv1)
      acd11(38)=spvae1k1(iv2)
      acd11(39)=spvak1e2(iv3)
      acd11(40)=spvae1k1(iv3)
      acd11(41)=spvak1e2(iv2)
      acd11(42)=spvae1k2(iv2)
      acd11(43)=spvak2e2(iv3)
      acd11(44)=spvae1k2(iv3)
      acd11(45)=spvak2e2(iv2)
      acd11(46)=spvae1l4(iv2)
      acd11(47)=spval4e2(iv3)
      acd11(48)=spvae1l4(iv3)
      acd11(49)=spval4e2(iv2)
      acd11(50)=spvae1k1(iv1)
      acd11(51)=spvak1e2(iv1)
      acd11(52)=spvae1k2(iv1)
      acd11(53)=spvak2e2(iv1)
      acd11(54)=spvae1l4(iv1)
      acd11(55)=spval4e2(iv1)
      acd11(56)=-acd11(48)*acd11(49)
      acd11(57)=-acd11(46)*acd11(47)
      acd11(58)=acd11(40)*acd11(41)
      acd11(59)=acd11(38)*acd11(39)
      acd11(56)=acd11(59)+acd11(58)+acd11(56)+acd11(57)
      acd11(56)=acd11(14)*acd11(56)
      acd11(57)=-acd11(30)*acd11(31)
      acd11(58)=-acd11(28)*acd11(29)
      acd11(59)=acd11(21)*acd11(22)
      acd11(60)=acd11(18)*acd11(19)
      acd11(57)=acd11(60)+acd11(59)+acd11(57)+acd11(58)
      acd11(57)=acd11(13)*acd11(57)
      acd11(58)=-acd11(48)*acd11(55)
      acd11(59)=-acd11(47)*acd11(54)
      acd11(60)=acd11(40)*acd11(51)
      acd11(61)=acd11(39)*acd11(50)
      acd11(58)=acd11(61)+acd11(60)+acd11(58)+acd11(59)
      acd11(58)=acd11(11)*acd11(58)
      acd11(59)=-acd11(30)*acd11(37)
      acd11(60)=-acd11(29)*acd11(36)
      acd11(61)=acd11(21)*acd11(33)
      acd11(62)=acd11(19)*acd11(32)
      acd11(59)=acd11(62)+acd11(61)+acd11(59)+acd11(60)
      acd11(59)=acd11(10)*acd11(59)
      acd11(60)=-acd11(49)*acd11(54)
      acd11(61)=-acd11(46)*acd11(55)
      acd11(62)=acd11(41)*acd11(50)
      acd11(63)=acd11(38)*acd11(51)
      acd11(60)=acd11(63)+acd11(62)+acd11(60)+acd11(61)
      acd11(60)=acd11(6)*acd11(60)
      acd11(61)=-acd11(31)*acd11(36)
      acd11(62)=-acd11(28)*acd11(37)
      acd11(63)=acd11(22)*acd11(32)
      acd11(64)=acd11(18)*acd11(33)
      acd11(61)=acd11(64)+acd11(63)+acd11(61)+acd11(62)
      acd11(61)=acd11(2)*acd11(61)
      acd11(56)=acd11(61)+acd11(60)+acd11(59)+acd11(58)+acd11(56)+acd11(57)
      acd11(56)=acd11(20)*acd11(56)
      acd11(57)=-acd11(44)*acd11(45)
      acd11(58)=-acd11(42)*acd11(43)
      acd11(57)=acd11(57)+acd11(58)
      acd11(57)=acd11(14)*acd11(57)
      acd11(58)=-acd11(26)*acd11(27)
      acd11(59)=-acd11(23)*acd11(24)
      acd11(58)=acd11(58)+acd11(59)
      acd11(58)=acd11(13)*acd11(58)
      acd11(59)=-acd11(44)*acd11(53)
      acd11(60)=-acd11(43)*acd11(52)
      acd11(59)=acd11(59)+acd11(60)
      acd11(59)=acd11(11)*acd11(59)
      acd11(60)=-acd11(26)*acd11(35)
      acd11(61)=-acd11(24)*acd11(34)
      acd11(60)=acd11(60)+acd11(61)
      acd11(60)=acd11(10)*acd11(60)
      acd11(61)=-acd11(45)*acd11(52)
      acd11(62)=-acd11(42)*acd11(53)
      acd11(61)=acd11(61)+acd11(62)
      acd11(61)=acd11(6)*acd11(61)
      acd11(62)=-acd11(27)*acd11(34)
      acd11(63)=-acd11(23)*acd11(35)
      acd11(62)=acd11(62)+acd11(63)
      acd11(62)=acd11(2)*acd11(62)
      acd11(57)=acd11(62)+acd11(61)+acd11(60)+acd11(59)+acd11(57)+acd11(58)
      acd11(57)=acd11(25)*acd11(57)
      acd11(58)=-acd11(1)*acd11(7)
      acd11(59)=-acd11(13)*acd11(16)
      acd11(60)=-acd11(10)*acd11(15)
      acd11(58)=acd11(60)+acd11(58)+acd11(59)
      acd11(58)=acd11(6)*acd11(58)
      acd11(59)=-acd11(1)*acd11(3)
      acd11(60)=-acd11(14)*acd11(16)
      acd11(61)=-acd11(11)*acd11(15)
      acd11(59)=acd11(61)+acd11(59)+acd11(60)
      acd11(59)=acd11(2)*acd11(59)
      acd11(60)=-acd11(9)*acd11(7)
      acd11(61)=-acd11(13)*acd11(17)
      acd11(60)=acd11(60)+acd11(61)
      acd11(60)=acd11(11)*acd11(60)
      acd11(61)=-acd11(9)*acd11(3)
      acd11(62)=-acd11(14)*acd11(17)
      acd11(61)=acd11(61)+acd11(62)
      acd11(61)=acd11(10)*acd11(61)
      acd11(62)=acd11(14)*acd11(12)
      acd11(63)=-acd11(7)*acd11(62)
      acd11(64)=acd11(13)*acd11(12)
      acd11(65)=-acd11(3)*acd11(64)
      acd11(58)=acd11(59)+acd11(58)+acd11(61)+acd11(60)+acd11(63)+acd11(65)
      acd11(58)=acd11(4)*acd11(58)
      acd11(59)=acd11(11)*acd11(9)
      acd11(60)=acd11(6)*acd11(1)
      acd11(59)=acd11(60)+acd11(59)+acd11(62)
      acd11(59)=acd11(8)*acd11(59)
      acd11(60)=acd11(10)*acd11(9)
      acd11(61)=acd11(2)*acd11(1)
      acd11(60)=acd11(61)+acd11(60)+acd11(64)
      acd11(60)=acd11(5)*acd11(60)
      acd11(58)=acd11(58)+acd11(60)+acd11(59)
      brack=acd11(56)+acd11(57)+2.0_ki*acd11(58)
   end function brack_4
!---#] function brack_4:
!---#[ function brack_5:
   pure function brack_5(Q, mu2) result(brack)
      use pb_part21part21_part25part25_model
      use pb_part21part21_part25part25_kinematics
      use pb_part21part21_part25part25_color
      use pb_part21part21_part25part25_abbrevd11h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(21) :: acd11
      complex(ki) :: brack
      acd11(1)=d(iv1,iv2)
      acd11(2)=spvae1e2(iv3)
      acd11(3)=spvae2e1(iv4)
      acd11(4)=abb11(24)
      acd11(5)=spvae1e2(iv4)
      acd11(6)=spvae2e1(iv3)
      acd11(7)=d(iv1,iv3)
      acd11(8)=spvae1e2(iv2)
      acd11(9)=spvae2e1(iv2)
      acd11(10)=d(iv1,iv4)
      acd11(11)=d(iv2,iv3)
      acd11(12)=spvae1e2(iv1)
      acd11(13)=spvae2e1(iv1)
      acd11(14)=d(iv2,iv4)
      acd11(15)=d(iv3,iv4)
      acd11(16)=acd11(13)*acd11(11)
      acd11(17)=acd11(9)*acd11(7)
      acd11(18)=acd11(6)*acd11(1)
      acd11(16)=acd11(18)+acd11(16)+acd11(17)
      acd11(16)=acd11(5)*acd11(16)
      acd11(17)=acd11(13)*acd11(14)
      acd11(18)=acd11(9)*acd11(10)
      acd11(19)=acd11(3)*acd11(1)
      acd11(17)=acd11(19)+acd11(17)+acd11(18)
      acd11(17)=acd11(2)*acd11(17)
      acd11(18)=acd11(9)*acd11(12)
      acd11(19)=acd11(8)*acd11(13)
      acd11(18)=acd11(18)+acd11(19)
      acd11(18)=acd11(15)*acd11(18)
      acd11(19)=acd11(12)*acd11(14)
      acd11(20)=acd11(8)*acd11(10)
      acd11(19)=acd11(19)+acd11(20)
      acd11(19)=acd11(6)*acd11(19)
      acd11(20)=acd11(12)*acd11(11)
      acd11(21)=acd11(8)*acd11(7)
      acd11(20)=acd11(20)+acd11(21)
      acd11(20)=acd11(3)*acd11(20)
      acd11(16)=acd11(17)+acd11(20)+acd11(16)+acd11(19)+acd11(18)
      brack=2.0_ki*acd11(16)*acd11(4)
   end function brack_5
!---#] function brack_5:
!---#[ function derivative:
   function derivative(mu2,i1,i2,i3,i4) result(numerator)
      use pb_part21part21_part25part25_globalsl1, only: epspow
      use pb_part21part21_part25part25_kinematics
      use pb_part21part21_part25part25_abbrevd11h0
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
      qshift = k2
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
!---#[ subroutine reconstruct_d11:
   subroutine     reconstruct_d11(coeffs)
      use pb_part21part21_part25part25_groups, only: tensrec_info_group1
      implicit none
      complex(ki), parameter :: czip = (0.0_ki, 0.0_ki)
      complex(ki), parameter :: cone = (1.0_ki, 0.0_ki)
      complex(ki), parameter :: ctwo = (2.0_ki, 0.0_ki)
      type(tensrec_info_group1), intent(out) :: coeffs
      ! rank 4 case :
      !---[# reconstruct coeffs%coeffs_11:
      coeffs%coeffs_11%c0 = derivative(czip)
      coeffs%coeffs_11%c1(1,1) = derivative(czip,1)
      coeffs%coeffs_11%c1(1,2) = derivative(czip,1,1)/2.0_ki
      coeffs%coeffs_11%c1(1,3) = derivative(czip,1,1,1)/6.0_ki
      coeffs%coeffs_11%c1(1,4) = derivative(czip,1,1,1,1)/24.0_ki
      coeffs%coeffs_11%c1(2,1) = -derivative(czip,2)
      coeffs%coeffs_11%c1(2,2) = derivative(czip,2,2)/2.0_ki
      coeffs%coeffs_11%c1(2,3) = -derivative(czip,2,2,2)/6.0_ki
      coeffs%coeffs_11%c1(2,4) = derivative(czip,2,2,2,2)/24.0_ki
      coeffs%coeffs_11%c1(3,1) = -derivative(czip,3)
      coeffs%coeffs_11%c1(3,2) = derivative(czip,3,3)/2.0_ki
      coeffs%coeffs_11%c1(3,3) = -derivative(czip,3,3,3)/6.0_ki
      coeffs%coeffs_11%c1(3,4) = derivative(czip,3,3,3,3)/24.0_ki
      coeffs%coeffs_11%c1(4,1) = -derivative(czip,4)
      coeffs%coeffs_11%c1(4,2) = derivative(czip,4,4)/2.0_ki
      coeffs%coeffs_11%c1(4,3) = -derivative(czip,4,4,4)/6.0_ki
      coeffs%coeffs_11%c1(4,4) = derivative(czip,4,4,4,4)/24.0_ki
      coeffs%coeffs_11%c2(1,1) = -derivative(czip,1,2)
      coeffs%coeffs_11%c2(1,2) = derivative(czip,1,2,2)/2.0_ki
      coeffs%coeffs_11%c2(1,3) = -derivative(czip,1,2,2,2)/6.0_ki
      coeffs%coeffs_11%c2(1,4) = -derivative(czip,1,1,2)/2.0_ki
      coeffs%coeffs_11%c2(1,5) = derivative(czip,1,1,2,2)/4.0_ki
      coeffs%coeffs_11%c2(1,6) = -derivative(czip,1,1,1,2)/6.0_ki
      coeffs%coeffs_11%c2(2,1) = -derivative(czip,1,3)
      coeffs%coeffs_11%c2(2,2) = derivative(czip,1,3,3)/2.0_ki
      coeffs%coeffs_11%c2(2,3) = -derivative(czip,1,3,3,3)/6.0_ki
      coeffs%coeffs_11%c2(2,4) = -derivative(czip,1,1,3)/2.0_ki
      coeffs%coeffs_11%c2(2,5) = derivative(czip,1,1,3,3)/4.0_ki
      coeffs%coeffs_11%c2(2,6) = -derivative(czip,1,1,1,3)/6.0_ki
      coeffs%coeffs_11%c2(3,1) = -derivative(czip,1,4)
      coeffs%coeffs_11%c2(3,2) = derivative(czip,1,4,4)/2.0_ki
      coeffs%coeffs_11%c2(3,3) = -derivative(czip,1,4,4,4)/6.0_ki
      coeffs%coeffs_11%c2(3,4) = -derivative(czip,1,1,4)/2.0_ki
      coeffs%coeffs_11%c2(3,5) = derivative(czip,1,1,4,4)/4.0_ki
      coeffs%coeffs_11%c2(3,6) = -derivative(czip,1,1,1,4)/6.0_ki
      coeffs%coeffs_11%c2(4,1) = derivative(czip,2,3)
      coeffs%coeffs_11%c2(4,2) = -derivative(czip,2,3,3)/2.0_ki
      coeffs%coeffs_11%c2(4,3) = derivative(czip,2,3,3,3)/6.0_ki
      coeffs%coeffs_11%c2(4,4) = -derivative(czip,2,2,3)/2.0_ki
      coeffs%coeffs_11%c2(4,5) = derivative(czip,2,2,3,3)/4.0_ki
      coeffs%coeffs_11%c2(4,6) = derivative(czip,2,2,2,3)/6.0_ki
      coeffs%coeffs_11%c2(5,1) = derivative(czip,2,4)
      coeffs%coeffs_11%c2(5,2) = -derivative(czip,2,4,4)/2.0_ki
      coeffs%coeffs_11%c2(5,3) = derivative(czip,2,4,4,4)/6.0_ki
      coeffs%coeffs_11%c2(5,4) = -derivative(czip,2,2,4)/2.0_ki
      coeffs%coeffs_11%c2(5,5) = derivative(czip,2,2,4,4)/4.0_ki
      coeffs%coeffs_11%c2(5,6) = derivative(czip,2,2,2,4)/6.0_ki
      coeffs%coeffs_11%c2(6,1) = derivative(czip,3,4)
      coeffs%coeffs_11%c2(6,2) = -derivative(czip,3,4,4)/2.0_ki
      coeffs%coeffs_11%c2(6,3) = derivative(czip,3,4,4,4)/6.0_ki
      coeffs%coeffs_11%c2(6,4) = -derivative(czip,3,3,4)/2.0_ki
      coeffs%coeffs_11%c2(6,5) = derivative(czip,3,3,4,4)/4.0_ki
      coeffs%coeffs_11%c2(6,6) = derivative(czip,3,3,3,4)/6.0_ki
      coeffs%coeffs_11%c3(1,1) = derivative(czip,1,2,3)
      coeffs%coeffs_11%c3(1,2) = -derivative(czip,1,2,3,3)/2.0_ki
      coeffs%coeffs_11%c3(1,3) = -derivative(czip,1,2,2,3)/2.0_ki
      coeffs%coeffs_11%c3(1,4) = derivative(czip,1,1,2,3)/2.0_ki
      coeffs%coeffs_11%c3(2,1) = derivative(czip,1,2,4)
      coeffs%coeffs_11%c3(2,2) = -derivative(czip,1,2,4,4)/2.0_ki
      coeffs%coeffs_11%c3(2,3) = -derivative(czip,1,2,2,4)/2.0_ki
      coeffs%coeffs_11%c3(2,4) = derivative(czip,1,1,2,4)/2.0_ki
      coeffs%coeffs_11%c3(3,1) = derivative(czip,1,3,4)
      coeffs%coeffs_11%c3(3,2) = -derivative(czip,1,3,4,4)/2.0_ki
      coeffs%coeffs_11%c3(3,3) = -derivative(czip,1,3,3,4)/2.0_ki
      coeffs%coeffs_11%c3(3,4) = derivative(czip,1,1,3,4)/2.0_ki
      coeffs%coeffs_11%c3(4,1) = -derivative(czip,2,3,4)
      coeffs%coeffs_11%c3(4,2) = derivative(czip,2,3,4,4)/2.0_ki
      coeffs%coeffs_11%c3(4,3) = derivative(czip,2,3,3,4)/2.0_ki
      coeffs%coeffs_11%c3(4,4) = derivative(czip,2,2,3,4)/2.0_ki
      coeffs%coeffs_11%c4(1,1) = -derivative(czip,1,2,3,4)
      !---#] reconstruct coeffs%coeffs_11:
   end subroutine reconstruct_d11
!---#] subroutine reconstruct_d11:
end module     pb_part21part21_part25part25_d11h0l1d
