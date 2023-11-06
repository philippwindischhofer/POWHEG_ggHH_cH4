module     p2_part21part21_part25part25part21_d133h0l1d_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d133h0l1d_qp.f90
   ! generator: buildfortran_d.py
   use p2_part21part21_part25part25part21_config, only: ki => ki_qp
   use p2_part21part21_part25part25part21_util_qp, only: cond, d => metric_tens&
   &or
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
      use p2_part21part21_part25part25part21_model_qp
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_color_qp
      use p2_part21part21_part25part25part21_abbrevd133h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(77) :: acd133
      complex(ki) :: brack
      acd133(1)=dotproduct(qshift,qshift)
      acd133(2)=abb133(12)
      acd133(3)=dotproduct(qshift,spvak1k2)
      acd133(4)=abb133(27)
      acd133(5)=dotproduct(qshift,spvak1k5)
      acd133(6)=abb133(52)
      acd133(7)=dotproduct(qshift,spvak2k1)
      acd133(8)=abb133(51)
      acd133(9)=dotproduct(qshift,spvak5k1)
      acd133(10)=abb133(8)
      acd133(11)=dotproduct(qshift,spvak1e1)
      acd133(12)=abb133(42)
      acd133(13)=dotproduct(qshift,spvae1k1)
      acd133(14)=abb133(26)
      acd133(15)=dotproduct(qshift,spvak2e1)
      acd133(16)=abb133(17)
      acd133(17)=dotproduct(qshift,spvae1k2)
      acd133(18)=abb133(15)
      acd133(19)=dotproduct(qshift,spval3e1)
      acd133(20)=abb133(45)
      acd133(21)=dotproduct(qshift,spvae1l3)
      acd133(22)=abb133(40)
      acd133(23)=dotproduct(qshift,spval4e1)
      acd133(24)=abb133(32)
      acd133(25)=dotproduct(qshift,spvae1l4)
      acd133(26)=abb133(28)
      acd133(27)=abb133(11)
      acd133(28)=abb133(30)
      acd133(29)=abb133(22)
      acd133(30)=abb133(39)
      acd133(31)=abb133(55)
      acd133(32)=abb133(53)
      acd133(33)=abb133(7)
      acd133(34)=abb133(48)
      acd133(35)=abb133(20)
      acd133(36)=abb133(50)
      acd133(37)=abb133(49)
      acd133(38)=abb133(47)
      acd133(39)=abb133(46)
      acd133(40)=abb133(44)
      acd133(41)=abb133(37)
      acd133(42)=abb133(34)
      acd133(43)=abb133(29)
      acd133(44)=abb133(10)
      acd133(45)=abb133(24)
      acd133(46)=abb133(23)
      acd133(47)=abb133(14)
      acd133(48)=abb133(36)
      acd133(49)=abb133(41)
      acd133(50)=abb133(18)
      acd133(51)=abb133(35)
      acd133(52)=abb133(43)
      acd133(53)=abb133(33)
      acd133(54)=abb133(38)
      acd133(55)=abb133(31)
      acd133(56)=abb133(25)
      acd133(57)=dotproduct(qshift,spvak5e1)
      acd133(58)=dotproduct(qshift,spvae1k5)
      acd133(59)=abb133(9)
      acd133(60)=abb133(19)
      acd133(61)=abb133(16)
      acd133(62)=abb133(21)
      acd133(63)=-acd133(7)*acd133(8)
      acd133(64)=-acd133(3)*acd133(4)
      acd133(65)=-acd133(25)*acd133(26)
      acd133(66)=-acd133(23)*acd133(24)
      acd133(67)=-acd133(13)*acd133(14)
      acd133(68)=-acd133(11)*acd133(12)
      acd133(69)=-acd133(9)*acd133(10)
      acd133(70)=-acd133(5)*acd133(6)
      acd133(71)=-acd133(21)*acd133(22)
      acd133(72)=-acd133(19)*acd133(20)
      acd133(73)=-acd133(17)*acd133(18)
      acd133(74)=-acd133(15)*acd133(16)
      acd133(75)=acd133(1)*acd133(2)
      acd133(63)=acd133(75)+acd133(74)+acd133(73)+acd133(72)+acd133(71)+acd133(&
      &70)+acd133(69)+acd133(68)+acd133(67)+acd133(66)+acd133(65)+acd133(64)+ac&
      &d133(27)+acd133(63)
      acd133(63)=acd133(1)*acd133(63)
      acd133(64)=acd133(25)*acd133(46)
      acd133(65)=acd133(13)*acd133(41)
      acd133(66)=acd133(9)*acd133(35)
      acd133(67)=acd133(21)*acd133(45)
      acd133(68)=acd133(17)*acd133(44)
      acd133(64)=acd133(68)+acd133(67)+acd133(66)+acd133(65)-acd133(47)+acd133(&
      &64)
      acd133(64)=acd133(15)*acd133(64)
      acd133(65)=acd133(3)*acd133(28)
      acd133(66)=acd133(23)*acd133(53)
      acd133(67)=acd133(11)*acd133(39)
      acd133(68)=acd133(5)*acd133(31)
      acd133(65)=acd133(68)+acd133(67)+acd133(66)-acd133(54)+acd133(65)
      acd133(65)=acd133(21)*acd133(65)
      acd133(66)=acd133(7)*acd133(33)
      acd133(67)=acd133(25)*acd133(51)
      acd133(68)=acd133(13)*acd133(42)
      acd133(69)=acd133(9)*acd133(36)
      acd133(66)=acd133(69)+acd133(68)+acd133(67)-acd133(52)+acd133(66)
      acd133(66)=acd133(19)*acd133(66)
      acd133(67)=acd133(23)*acd133(49)
      acd133(68)=acd133(11)*acd133(38)
      acd133(69)=acd133(5)*acd133(30)
      acd133(70)=acd133(19)*acd133(48)
      acd133(67)=acd133(70)+acd133(69)+acd133(68)-acd133(50)+acd133(67)
      acd133(67)=acd133(17)*acd133(67)
      acd133(68)=-acd133(58)*acd133(61)
      acd133(69)=acd133(58)*acd133(59)
      acd133(69)=-acd133(60)+acd133(69)
      acd133(69)=acd133(57)*acd133(69)
      acd133(70)=-acd133(7)*acd133(34)
      acd133(71)=-acd133(3)*acd133(29)
      acd133(72)=-acd133(25)*acd133(56)
      acd133(73)=-acd133(23)*acd133(55)
      acd133(74)=-acd133(13)*acd133(43)
      acd133(75)=-acd133(11)*acd133(40)
      acd133(76)=-acd133(9)*acd133(37)
      acd133(77)=-acd133(5)*acd133(32)
      brack=acd133(62)+acd133(63)+acd133(64)+acd133(65)+acd133(66)+acd133(67)+a&
      &cd133(68)+acd133(69)+acd133(70)+acd133(71)+acd133(72)+acd133(73)+acd133(&
      &74)+acd133(75)+acd133(76)+acd133(77)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p2_part21part21_part25part25part21_model_qp
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_color_qp
      use p2_part21part21_part25part25part21_abbrevd133h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(96) :: acd133
      complex(ki) :: brack
      acd133(1)=qshift(iv1)
      acd133(2)=dotproduct(qshift,qshift)
      acd133(3)=abb133(12)
      acd133(4)=dotproduct(qshift,spvak1k2)
      acd133(5)=abb133(27)
      acd133(6)=dotproduct(qshift,spvak1k5)
      acd133(7)=abb133(52)
      acd133(8)=dotproduct(qshift,spvak2k1)
      acd133(9)=abb133(51)
      acd133(10)=dotproduct(qshift,spvak5k1)
      acd133(11)=abb133(8)
      acd133(12)=dotproduct(qshift,spvak1e1)
      acd133(13)=abb133(42)
      acd133(14)=dotproduct(qshift,spvae1k1)
      acd133(15)=abb133(26)
      acd133(16)=dotproduct(qshift,spvak2e1)
      acd133(17)=abb133(17)
      acd133(18)=dotproduct(qshift,spvae1k2)
      acd133(19)=abb133(15)
      acd133(20)=dotproduct(qshift,spval3e1)
      acd133(21)=abb133(45)
      acd133(22)=dotproduct(qshift,spvae1l3)
      acd133(23)=abb133(40)
      acd133(24)=dotproduct(qshift,spval4e1)
      acd133(25)=abb133(32)
      acd133(26)=dotproduct(qshift,spvae1l4)
      acd133(27)=abb133(28)
      acd133(28)=abb133(11)
      acd133(29)=spvak1k2(iv1)
      acd133(30)=abb133(30)
      acd133(31)=abb133(22)
      acd133(32)=spvak1k5(iv1)
      acd133(33)=abb133(39)
      acd133(34)=abb133(55)
      acd133(35)=abb133(53)
      acd133(36)=spvak2k1(iv1)
      acd133(37)=abb133(7)
      acd133(38)=abb133(48)
      acd133(39)=spvak5k1(iv1)
      acd133(40)=abb133(20)
      acd133(41)=abb133(50)
      acd133(42)=abb133(49)
      acd133(43)=spvak1e1(iv1)
      acd133(44)=abb133(47)
      acd133(45)=abb133(46)
      acd133(46)=abb133(44)
      acd133(47)=spvae1k1(iv1)
      acd133(48)=abb133(37)
      acd133(49)=abb133(34)
      acd133(50)=abb133(29)
      acd133(51)=spvak2e1(iv1)
      acd133(52)=abb133(10)
      acd133(53)=abb133(24)
      acd133(54)=abb133(23)
      acd133(55)=abb133(14)
      acd133(56)=spvae1k2(iv1)
      acd133(57)=abb133(36)
      acd133(58)=abb133(41)
      acd133(59)=abb133(18)
      acd133(60)=spval3e1(iv1)
      acd133(61)=abb133(35)
      acd133(62)=abb133(43)
      acd133(63)=spvae1l3(iv1)
      acd133(64)=abb133(33)
      acd133(65)=abb133(38)
      acd133(66)=spval4e1(iv1)
      acd133(67)=abb133(31)
      acd133(68)=spvae1l4(iv1)
      acd133(69)=abb133(25)
      acd133(70)=spvak5e1(iv1)
      acd133(71)=dotproduct(qshift,spvae1k5)
      acd133(72)=abb133(9)
      acd133(73)=abb133(19)
      acd133(74)=spvae1k5(iv1)
      acd133(75)=dotproduct(qshift,spvak5e1)
      acd133(76)=abb133(16)
      acd133(77)=acd133(8)*acd133(9)
      acd133(78)=acd133(4)*acd133(5)
      acd133(79)=acd133(26)*acd133(27)
      acd133(80)=acd133(24)*acd133(25)
      acd133(81)=acd133(14)*acd133(15)
      acd133(82)=acd133(12)*acd133(13)
      acd133(83)=acd133(10)*acd133(11)
      acd133(84)=acd133(6)*acd133(7)
      acd133(85)=acd133(22)*acd133(23)
      acd133(86)=acd133(20)*acd133(21)
      acd133(87)=acd133(18)*acd133(19)
      acd133(88)=acd133(16)*acd133(17)
      acd133(89)=acd133(2)*acd133(3)
      acd133(77)=-2.0_ki*acd133(89)+acd133(88)+acd133(87)+acd133(86)+acd133(85)&
      &+acd133(84)+acd133(83)+acd133(82)+acd133(81)+acd133(80)+acd133(79)+acd13&
      &3(78)-acd133(28)+acd133(77)
      acd133(77)=acd133(1)*acd133(77)
      acd133(78)=acd133(36)*acd133(9)
      acd133(79)=acd133(29)*acd133(5)
      acd133(80)=acd133(68)*acd133(27)
      acd133(81)=acd133(66)*acd133(25)
      acd133(82)=acd133(47)*acd133(15)
      acd133(83)=acd133(43)*acd133(13)
      acd133(84)=acd133(39)*acd133(11)
      acd133(85)=acd133(32)*acd133(7)
      acd133(86)=acd133(63)*acd133(23)
      acd133(87)=acd133(60)*acd133(21)
      acd133(88)=acd133(56)*acd133(19)
      acd133(89)=acd133(51)*acd133(17)
      acd133(78)=acd133(89)+acd133(88)+acd133(87)+acd133(86)+acd133(85)+acd133(&
      &84)+acd133(83)+acd133(82)+acd133(81)+acd133(80)+acd133(78)+acd133(79)
      acd133(78)=acd133(2)*acd133(78)
      acd133(79)=-acd133(4)*acd133(30)
      acd133(80)=-acd133(24)*acd133(64)
      acd133(81)=-acd133(12)*acd133(45)
      acd133(82)=-acd133(6)*acd133(34)
      acd133(83)=-acd133(16)*acd133(53)
      acd133(79)=acd133(83)+acd133(82)+acd133(81)+acd133(80)+acd133(65)+acd133(&
      &79)
      acd133(79)=acd133(63)*acd133(79)
      acd133(80)=-acd133(8)*acd133(37)
      acd133(81)=-acd133(26)*acd133(61)
      acd133(82)=-acd133(14)*acd133(49)
      acd133(83)=-acd133(10)*acd133(41)
      acd133(84)=-acd133(18)*acd133(57)
      acd133(80)=acd133(84)+acd133(83)+acd133(82)+acd133(81)+acd133(62)+acd133(&
      &80)
      acd133(80)=acd133(60)*acd133(80)
      acd133(81)=-acd133(24)*acd133(58)
      acd133(82)=-acd133(12)*acd133(44)
      acd133(83)=-acd133(6)*acd133(33)
      acd133(84)=-acd133(20)*acd133(57)
      acd133(85)=-acd133(16)*acd133(52)
      acd133(81)=acd133(85)+acd133(84)+acd133(83)+acd133(82)+acd133(59)+acd133(&
      &81)
      acd133(81)=acd133(56)*acd133(81)
      acd133(82)=-acd133(26)*acd133(54)
      acd133(83)=-acd133(14)*acd133(48)
      acd133(84)=-acd133(10)*acd133(40)
      acd133(85)=-acd133(22)*acd133(53)
      acd133(86)=-acd133(18)*acd133(52)
      acd133(82)=acd133(86)+acd133(85)+acd133(84)+acd133(83)+acd133(55)+acd133(&
      &82)
      acd133(82)=acd133(51)*acd133(82)
      acd133(83)=-acd133(29)*acd133(30)
      acd133(84)=-acd133(66)*acd133(64)
      acd133(85)=-acd133(43)*acd133(45)
      acd133(86)=-acd133(32)*acd133(34)
      acd133(83)=acd133(86)+acd133(85)+acd133(83)+acd133(84)
      acd133(83)=acd133(22)*acd133(83)
      acd133(84)=-acd133(36)*acd133(37)
      acd133(85)=-acd133(68)*acd133(61)
      acd133(86)=-acd133(47)*acd133(49)
      acd133(87)=-acd133(39)*acd133(41)
      acd133(84)=acd133(87)+acd133(86)+acd133(84)+acd133(85)
      acd133(84)=acd133(20)*acd133(84)
      acd133(85)=-acd133(66)*acd133(58)
      acd133(86)=-acd133(43)*acd133(44)
      acd133(87)=-acd133(32)*acd133(33)
      acd133(85)=acd133(87)+acd133(85)+acd133(86)
      acd133(85)=acd133(18)*acd133(85)
      acd133(86)=-acd133(68)*acd133(54)
      acd133(87)=-acd133(47)*acd133(48)
      acd133(88)=-acd133(39)*acd133(40)
      acd133(86)=acd133(88)+acd133(86)+acd133(87)
      acd133(86)=acd133(16)*acd133(86)
      acd133(87)=-acd133(72)*acd133(75)
      acd133(87)=acd133(76)+acd133(87)
      acd133(87)=acd133(74)*acd133(87)
      acd133(88)=-acd133(72)*acd133(71)
      acd133(88)=acd133(73)+acd133(88)
      acd133(88)=acd133(70)*acd133(88)
      acd133(89)=acd133(36)*acd133(38)
      acd133(90)=acd133(29)*acd133(31)
      acd133(91)=acd133(68)*acd133(69)
      acd133(92)=acd133(66)*acd133(67)
      acd133(93)=acd133(47)*acd133(50)
      acd133(94)=acd133(43)*acd133(46)
      acd133(95)=acd133(39)*acd133(42)
      acd133(96)=acd133(32)*acd133(35)
      brack=2.0_ki*acd133(77)+acd133(78)+acd133(79)+acd133(80)+acd133(81)+acd13&
      &3(82)+acd133(83)+acd133(84)+acd133(85)+acd133(86)+acd133(87)+acd133(88)+&
      &acd133(89)+acd133(90)+acd133(91)+acd133(92)+acd133(93)+acd133(94)+acd133&
      &(95)+acd133(96)
   end function brack_2
!---#] function brack_2:
!---#[ function brack_3:
   pure function brack_3(Q, mu2) result(brack)
      use p2_part21part21_part25part25part21_model_qp
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_color_qp
      use p2_part21part21_part25part25part21_abbrevd133h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(90) :: acd133
      complex(ki) :: brack
      acd133(1)=d(iv1,iv2)
      acd133(2)=dotproduct(qshift,qshift)
      acd133(3)=abb133(12)
      acd133(4)=dotproduct(qshift,spvak1k2)
      acd133(5)=abb133(27)
      acd133(6)=dotproduct(qshift,spvak1k5)
      acd133(7)=abb133(52)
      acd133(8)=dotproduct(qshift,spvak2k1)
      acd133(9)=abb133(51)
      acd133(10)=dotproduct(qshift,spvak5k1)
      acd133(11)=abb133(8)
      acd133(12)=dotproduct(qshift,spvak1e1)
      acd133(13)=abb133(42)
      acd133(14)=dotproduct(qshift,spvae1k1)
      acd133(15)=abb133(26)
      acd133(16)=dotproduct(qshift,spvak2e1)
      acd133(17)=abb133(17)
      acd133(18)=dotproduct(qshift,spvae1k2)
      acd133(19)=abb133(15)
      acd133(20)=dotproduct(qshift,spval3e1)
      acd133(21)=abb133(45)
      acd133(22)=dotproduct(qshift,spvae1l3)
      acd133(23)=abb133(40)
      acd133(24)=dotproduct(qshift,spval4e1)
      acd133(25)=abb133(32)
      acd133(26)=dotproduct(qshift,spvae1l4)
      acd133(27)=abb133(28)
      acd133(28)=abb133(11)
      acd133(29)=qshift(iv1)
      acd133(30)=qshift(iv2)
      acd133(31)=spvak1k2(iv2)
      acd133(32)=spvak1k5(iv2)
      acd133(33)=spvak2k1(iv2)
      acd133(34)=spvak5k1(iv2)
      acd133(35)=spvak1e1(iv2)
      acd133(36)=spvae1k1(iv2)
      acd133(37)=spvak2e1(iv2)
      acd133(38)=spvae1k2(iv2)
      acd133(39)=spval3e1(iv2)
      acd133(40)=spvae1l3(iv2)
      acd133(41)=spval4e1(iv2)
      acd133(42)=spvae1l4(iv2)
      acd133(43)=spvak1k2(iv1)
      acd133(44)=spvak1k5(iv1)
      acd133(45)=spvak2k1(iv1)
      acd133(46)=spvak5k1(iv1)
      acd133(47)=spvak1e1(iv1)
      acd133(48)=spvae1k1(iv1)
      acd133(49)=spvak2e1(iv1)
      acd133(50)=spvae1k2(iv1)
      acd133(51)=spval3e1(iv1)
      acd133(52)=spvae1l3(iv1)
      acd133(53)=spval4e1(iv1)
      acd133(54)=spvae1l4(iv1)
      acd133(55)=abb133(30)
      acd133(56)=abb133(39)
      acd133(57)=abb133(55)
      acd133(58)=abb133(7)
      acd133(59)=abb133(20)
      acd133(60)=abb133(50)
      acd133(61)=abb133(47)
      acd133(62)=abb133(46)
      acd133(63)=abb133(37)
      acd133(64)=abb133(34)
      acd133(65)=abb133(10)
      acd133(66)=abb133(24)
      acd133(67)=abb133(23)
      acd133(68)=abb133(36)
      acd133(69)=abb133(41)
      acd133(70)=abb133(35)
      acd133(71)=abb133(33)
      acd133(72)=spvak5e1(iv1)
      acd133(73)=spvae1k5(iv2)
      acd133(74)=abb133(9)
      acd133(75)=spvak5e1(iv2)
      acd133(76)=spvae1k5(iv1)
      acd133(77)=-acd133(27)*acd133(42)
      acd133(78)=-acd133(25)*acd133(41)
      acd133(79)=-acd133(15)*acd133(36)
      acd133(80)=-acd133(13)*acd133(35)
      acd133(81)=-acd133(11)*acd133(34)
      acd133(82)=-acd133(9)*acd133(33)
      acd133(83)=-acd133(7)*acd133(32)
      acd133(84)=-acd133(5)*acd133(31)
      acd133(85)=-acd133(40)*acd133(23)
      acd133(86)=-acd133(39)*acd133(21)
      acd133(87)=-acd133(38)*acd133(19)
      acd133(88)=-acd133(37)*acd133(17)
      acd133(89)=acd133(30)*acd133(3)
      acd133(77)=4.0_ki*acd133(89)+acd133(88)+acd133(87)+acd133(86)+acd133(85)+&
      &acd133(84)+acd133(83)+acd133(82)+acd133(81)+acd133(80)+acd133(79)+acd133&
      &(77)+acd133(78)
      acd133(77)=acd133(29)*acd133(77)
      acd133(78)=acd133(3)*acd133(2)
      acd133(79)=-acd133(27)*acd133(26)
      acd133(80)=-acd133(25)*acd133(24)
      acd133(81)=-acd133(23)*acd133(22)
      acd133(82)=-acd133(21)*acd133(20)
      acd133(83)=-acd133(19)*acd133(18)
      acd133(84)=-acd133(17)*acd133(16)
      acd133(85)=-acd133(15)*acd133(14)
      acd133(86)=-acd133(13)*acd133(12)
      acd133(87)=-acd133(11)*acd133(10)
      acd133(88)=-acd133(9)*acd133(8)
      acd133(89)=-acd133(7)*acd133(6)
      acd133(90)=-acd133(5)*acd133(4)
      acd133(78)=acd133(90)+acd133(89)+acd133(88)+acd133(87)+acd133(86)+acd133(&
      &85)+acd133(84)+acd133(83)+acd133(82)+acd133(81)+acd133(80)+acd133(79)+ac&
      &d133(28)+2.0_ki*acd133(78)
      acd133(78)=acd133(1)*acd133(78)
      acd133(79)=-acd133(27)*acd133(54)
      acd133(80)=-acd133(25)*acd133(53)
      acd133(81)=-acd133(15)*acd133(48)
      acd133(82)=-acd133(13)*acd133(47)
      acd133(83)=-acd133(11)*acd133(46)
      acd133(84)=-acd133(9)*acd133(45)
      acd133(85)=-acd133(7)*acd133(44)
      acd133(86)=-acd133(5)*acd133(43)
      acd133(87)=-acd133(52)*acd133(23)
      acd133(88)=-acd133(51)*acd133(21)
      acd133(89)=-acd133(50)*acd133(19)
      acd133(90)=-acd133(49)*acd133(17)
      acd133(79)=acd133(90)+acd133(89)+acd133(88)+acd133(87)+acd133(86)+acd133(&
      &85)+acd133(84)+acd133(83)+acd133(82)+acd133(81)+acd133(79)+acd133(80)
      acd133(79)=acd133(30)*acd133(79)
      acd133(77)=acd133(78)+acd133(79)+acd133(77)
      acd133(78)=acd133(43)*acd133(55)
      acd133(79)=acd133(53)*acd133(71)
      acd133(80)=acd133(47)*acd133(62)
      acd133(81)=acd133(44)*acd133(57)
      acd133(82)=acd133(49)*acd133(66)
      acd133(78)=acd133(82)+acd133(81)+acd133(80)+acd133(78)+acd133(79)
      acd133(78)=acd133(40)*acd133(78)
      acd133(79)=acd133(45)*acd133(58)
      acd133(80)=acd133(54)*acd133(70)
      acd133(81)=acd133(48)*acd133(64)
      acd133(82)=acd133(46)*acd133(60)
      acd133(83)=acd133(50)*acd133(68)
      acd133(79)=acd133(83)+acd133(82)+acd133(81)+acd133(79)+acd133(80)
      acd133(79)=acd133(39)*acd133(79)
      acd133(80)=acd133(53)*acd133(69)
      acd133(81)=acd133(47)*acd133(61)
      acd133(82)=acd133(44)*acd133(56)
      acd133(83)=acd133(51)*acd133(68)
      acd133(84)=acd133(49)*acd133(65)
      acd133(80)=acd133(84)+acd133(83)+acd133(82)+acd133(80)+acd133(81)
      acd133(80)=acd133(38)*acd133(80)
      acd133(81)=acd133(54)*acd133(67)
      acd133(82)=acd133(48)*acd133(63)
      acd133(83)=acd133(46)*acd133(59)
      acd133(84)=acd133(52)*acd133(66)
      acd133(85)=acd133(50)*acd133(65)
      acd133(81)=acd133(85)+acd133(84)+acd133(83)+acd133(81)+acd133(82)
      acd133(81)=acd133(37)*acd133(81)
      acd133(82)=acd133(31)*acd133(55)
      acd133(83)=acd133(41)*acd133(71)
      acd133(84)=acd133(35)*acd133(62)
      acd133(85)=acd133(32)*acd133(57)
      acd133(82)=acd133(85)+acd133(84)+acd133(82)+acd133(83)
      acd133(82)=acd133(52)*acd133(82)
      acd133(83)=acd133(33)*acd133(58)
      acd133(84)=acd133(42)*acd133(70)
      acd133(85)=acd133(36)*acd133(64)
      acd133(86)=acd133(34)*acd133(60)
      acd133(83)=acd133(86)+acd133(85)+acd133(83)+acd133(84)
      acd133(83)=acd133(51)*acd133(83)
      acd133(84)=acd133(41)*acd133(69)
      acd133(85)=acd133(35)*acd133(61)
      acd133(86)=acd133(32)*acd133(56)
      acd133(84)=acd133(86)+acd133(84)+acd133(85)
      acd133(84)=acd133(50)*acd133(84)
      acd133(85)=acd133(42)*acd133(67)
      acd133(86)=acd133(36)*acd133(63)
      acd133(87)=acd133(34)*acd133(59)
      acd133(85)=acd133(87)+acd133(85)+acd133(86)
      acd133(85)=acd133(49)*acd133(85)
      acd133(86)=acd133(75)*acd133(76)
      acd133(87)=acd133(72)*acd133(73)
      acd133(86)=acd133(86)+acd133(87)
      acd133(86)=acd133(74)*acd133(86)
      brack=2.0_ki*acd133(77)+acd133(78)+acd133(79)+acd133(80)+acd133(81)+acd13&
      &3(82)+acd133(83)+acd133(84)+acd133(85)+acd133(86)
   end function brack_3
!---#] function brack_3:
!---#[ function brack_4:
   pure function brack_4(Q, mu2) result(brack)
      use p2_part21part21_part25part25part21_model_qp
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_color_qp
      use p2_part21part21_part25part25part21_abbrevd133h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(70) :: acd133
      complex(ki) :: brack
      acd133(1)=d(iv1,iv2)
      acd133(2)=qshift(iv3)
      acd133(3)=abb133(12)
      acd133(4)=spvak1k2(iv3)
      acd133(5)=abb133(27)
      acd133(6)=spvak1k5(iv3)
      acd133(7)=abb133(52)
      acd133(8)=spvak2k1(iv3)
      acd133(9)=abb133(51)
      acd133(10)=spvak5k1(iv3)
      acd133(11)=abb133(8)
      acd133(12)=spvak1e1(iv3)
      acd133(13)=abb133(42)
      acd133(14)=spvae1k1(iv3)
      acd133(15)=abb133(26)
      acd133(16)=spvak2e1(iv3)
      acd133(17)=abb133(17)
      acd133(18)=spvae1k2(iv3)
      acd133(19)=abb133(15)
      acd133(20)=spval3e1(iv3)
      acd133(21)=abb133(45)
      acd133(22)=spvae1l3(iv3)
      acd133(23)=abb133(40)
      acd133(24)=spval4e1(iv3)
      acd133(25)=abb133(32)
      acd133(26)=spvae1l4(iv3)
      acd133(27)=abb133(28)
      acd133(28)=d(iv1,iv3)
      acd133(29)=qshift(iv2)
      acd133(30)=spvak1k2(iv2)
      acd133(31)=spvak1k5(iv2)
      acd133(32)=spvak2k1(iv2)
      acd133(33)=spvak5k1(iv2)
      acd133(34)=spvak1e1(iv2)
      acd133(35)=spvae1k1(iv2)
      acd133(36)=spvak2e1(iv2)
      acd133(37)=spvae1k2(iv2)
      acd133(38)=spval3e1(iv2)
      acd133(39)=spvae1l3(iv2)
      acd133(40)=spval4e1(iv2)
      acd133(41)=spvae1l4(iv2)
      acd133(42)=d(iv2,iv3)
      acd133(43)=qshift(iv1)
      acd133(44)=spvak1k2(iv1)
      acd133(45)=spvak1k5(iv1)
      acd133(46)=spvak2k1(iv1)
      acd133(47)=spvak5k1(iv1)
      acd133(48)=spvak1e1(iv1)
      acd133(49)=spvae1k1(iv1)
      acd133(50)=spvak2e1(iv1)
      acd133(51)=spvae1k2(iv1)
      acd133(52)=spval3e1(iv1)
      acd133(53)=spvae1l3(iv1)
      acd133(54)=spval4e1(iv1)
      acd133(55)=spvae1l4(iv1)
      acd133(56)=4.0_ki*acd133(3)
      acd133(57)=-acd133(2)*acd133(56)
      acd133(58)=acd133(4)*acd133(5)
      acd133(59)=acd133(6)*acd133(7)
      acd133(60)=acd133(8)*acd133(9)
      acd133(61)=acd133(10)*acd133(11)
      acd133(62)=acd133(12)*acd133(13)
      acd133(63)=acd133(14)*acd133(15)
      acd133(64)=acd133(16)*acd133(17)
      acd133(65)=acd133(18)*acd133(19)
      acd133(66)=acd133(20)*acd133(21)
      acd133(67)=acd133(22)*acd133(23)
      acd133(68)=acd133(24)*acd133(25)
      acd133(69)=acd133(26)*acd133(27)
      acd133(57)=acd133(69)+acd133(68)+acd133(67)+acd133(66)+acd133(65)+acd133(&
      &64)+acd133(63)+acd133(62)+acd133(61)+acd133(60)+acd133(59)+acd133(57)+ac&
      &d133(58)
      acd133(57)=acd133(1)*acd133(57)
      acd133(58)=-acd133(29)*acd133(56)
      acd133(59)=acd133(30)*acd133(5)
      acd133(60)=acd133(31)*acd133(7)
      acd133(61)=acd133(32)*acd133(9)
      acd133(62)=acd133(33)*acd133(11)
      acd133(63)=acd133(34)*acd133(13)
      acd133(64)=acd133(35)*acd133(15)
      acd133(65)=acd133(36)*acd133(17)
      acd133(66)=acd133(37)*acd133(19)
      acd133(67)=acd133(38)*acd133(21)
      acd133(68)=acd133(39)*acd133(23)
      acd133(69)=acd133(40)*acd133(25)
      acd133(70)=acd133(41)*acd133(27)
      acd133(58)=acd133(70)+acd133(69)+acd133(68)+acd133(67)+acd133(66)+acd133(&
      &65)+acd133(64)+acd133(63)+acd133(62)+acd133(61)+acd133(60)+acd133(59)+ac&
      &d133(58)
      acd133(58)=acd133(28)*acd133(58)
      acd133(56)=-acd133(43)*acd133(56)
      acd133(59)=acd133(44)*acd133(5)
      acd133(60)=acd133(45)*acd133(7)
      acd133(61)=acd133(46)*acd133(9)
      acd133(62)=acd133(47)*acd133(11)
      acd133(63)=acd133(48)*acd133(13)
      acd133(64)=acd133(49)*acd133(15)
      acd133(65)=acd133(50)*acd133(17)
      acd133(66)=acd133(51)*acd133(19)
      acd133(67)=acd133(52)*acd133(21)
      acd133(68)=acd133(53)*acd133(23)
      acd133(69)=acd133(54)*acd133(25)
      acd133(70)=acd133(55)*acd133(27)
      acd133(56)=acd133(70)+acd133(69)+acd133(68)+acd133(67)+acd133(66)+acd133(&
      &65)+acd133(64)+acd133(63)+acd133(62)+acd133(61)+acd133(60)+acd133(59)+ac&
      &d133(56)
      acd133(56)=acd133(42)*acd133(56)
      acd133(56)=acd133(56)+acd133(58)+acd133(57)
      brack=2.0_ki*acd133(56)
   end function brack_4
!---#] function brack_4:
!---#[ function brack_5:
   pure function brack_5(Q, mu2) result(brack)
      use p2_part21part21_part25part25part21_model_qp
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_color_qp
      use p2_part21part21_part25part25part21_abbrevd133h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(10) :: acd133
      complex(ki) :: brack
      acd133(1)=d(iv1,iv2)
      acd133(2)=d(iv3,iv4)
      acd133(3)=abb133(12)
      acd133(4)=d(iv1,iv3)
      acd133(5)=d(iv2,iv4)
      acd133(6)=d(iv1,iv4)
      acd133(7)=d(iv2,iv3)
      acd133(8)=acd133(2)*acd133(1)
      acd133(9)=acd133(5)*acd133(4)
      acd133(10)=acd133(7)*acd133(6)
      acd133(8)=acd133(10)+acd133(8)+acd133(9)
      brack=8.0_ki*acd133(8)*acd133(3)
   end function brack_5
!---#] function brack_5:
!---#[ function derivative:
   function derivative(mu2,i1,i2,i3,i4) result(numerator)
      use p2_part21part21_part25part25part21_globalsl1_qp, only: epspow
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_abbrevd133h0_qp
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
      qshift = -k3+k2-k5
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
end module     p2_part21part21_part25part25part21_d133h0l1d_qp
