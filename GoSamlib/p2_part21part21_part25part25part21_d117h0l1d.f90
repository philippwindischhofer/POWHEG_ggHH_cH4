module     p2_part21part21_part25part25part21_d117h0l1d
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d117h0l1d.f90
   ! generator: buildfortran_d.py
   use p2_part21part21_part25part25part21_config, only: ki
   use p2_part21part21_part25part25part21_util, only: cond, d => metric_tensor
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
      use p2_part21part21_part25part25part21_model
      use p2_part21part21_part25part25part21_kinematics
      use p2_part21part21_part25part25part21_color
      use p2_part21part21_part25part25part21_abbrevd117h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(96) :: acd117
      complex(ki) :: brack
      acd117(1)=dotproduct(qshift,qshift)
      acd117(2)=abb117(35)
      acd117(3)=dotproduct(qshift,spvak1k5)
      acd117(4)=abb117(34)
      acd117(5)=dotproduct(qshift,spvak2k5)
      acd117(6)=abb117(38)
      acd117(7)=dotproduct(qshift,spvak5k1)
      acd117(8)=abb117(41)
      acd117(9)=dotproduct(qshift,spvak5k2)
      acd117(10)=abb117(28)
      acd117(11)=dotproduct(qshift,spvak2e5)
      acd117(12)=abb117(10)
      acd117(13)=dotproduct(qshift,spvae5k2)
      acd117(14)=abb117(15)
      acd117(15)=dotproduct(qshift,spval3e5)
      acd117(16)=abb117(13)
      acd117(17)=dotproduct(qshift,spvae5l3)
      acd117(18)=abb117(60)
      acd117(19)=dotproduct(qshift,spval4e5)
      acd117(20)=abb117(18)
      acd117(21)=dotproduct(qshift,spvae5l4)
      acd117(22)=abb117(44)
      acd117(23)=dotproduct(qshift,spvak5e1)
      acd117(24)=abb117(72)
      acd117(25)=dotproduct(qshift,spvae1k5)
      acd117(26)=abb117(69)
      acd117(27)=dotproduct(qshift,spvak5e5)
      acd117(28)=abb117(66)
      acd117(29)=dotproduct(qshift,spvae5k5)
      acd117(30)=abb117(56)
      acd117(31)=abb117(22)
      acd117(32)=abb117(30)
      acd117(33)=abb117(42)
      acd117(34)=abb117(77)
      acd117(35)=abb117(40)
      acd117(36)=abb117(25)
      acd117(37)=abb117(39)
      acd117(38)=abb117(36)
      acd117(39)=abb117(74)
      acd117(40)=abb117(29)
      acd117(41)=abb117(20)
      acd117(42)=abb117(8)
      acd117(43)=abb117(17)
      acd117(44)=abb117(12)
      acd117(45)=abb117(26)
      acd117(46)=abb117(16)
      acd117(47)=abb117(9)
      acd117(48)=abb117(37)
      acd117(49)=abb117(31)
      acd117(50)=abb117(33)
      acd117(51)=abb117(24)
      acd117(52)=abb117(11)
      acd117(53)=abb117(78)
      acd117(54)=abb117(80)
      acd117(55)=abb117(73)
      acd117(56)=abb117(7)
      acd117(57)=abb117(65)
      acd117(58)=abb117(76)
      acd117(59)=abb117(71)
      acd117(60)=abb117(51)
      acd117(61)=abb117(54)
      acd117(62)=abb117(46)
      acd117(63)=abb117(70)
      acd117(64)=abb117(68)
      acd117(65)=abb117(57)
      acd117(66)=abb117(52)
      acd117(67)=dotproduct(qshift,spvak1e5)
      acd117(68)=dotproduct(qshift,spvae5k1)
      acd117(69)=abb117(19)
      acd117(70)=abb117(27)
      acd117(71)=abb117(23)
      acd117(72)=dotproduct(qshift,spvae1e5)
      acd117(73)=dotproduct(qshift,spvae5e1)
      acd117(74)=abb117(49)
      acd117(75)=abb117(43)
      acd117(76)=abb117(48)
      acd117(77)=abb117(21)
      acd117(78)=-acd117(9)*acd117(10)
      acd117(79)=-acd117(5)*acd117(6)
      acd117(80)=-acd117(29)*acd117(30)
      acd117(81)=-acd117(27)*acd117(28)
      acd117(82)=-acd117(25)*acd117(26)
      acd117(83)=-acd117(23)*acd117(24)
      acd117(84)=-acd117(21)*acd117(22)
      acd117(85)=-acd117(19)*acd117(20)
      acd117(86)=-acd117(7)*acd117(8)
      acd117(87)=-acd117(3)*acd117(4)
      acd117(88)=-acd117(17)*acd117(18)
      acd117(89)=-acd117(15)*acd117(16)
      acd117(90)=-acd117(13)*acd117(14)
      acd117(91)=-acd117(11)*acd117(12)
      acd117(92)=acd117(1)*acd117(2)
      acd117(78)=acd117(92)+acd117(91)+acd117(90)+acd117(89)+acd117(88)+acd117(&
      &87)+acd117(86)+acd117(85)+acd117(84)+acd117(83)+acd117(82)+acd117(81)+ac&
      &d117(80)+acd117(79)+acd117(31)+acd117(78)
      acd117(78)=acd117(1)*acd117(78)
      acd117(79)=acd117(29)*acd117(46)
      acd117(80)=acd117(25)*acd117(45)
      acd117(81)=acd117(21)*acd117(44)
      acd117(82)=acd117(3)*acd117(32)
      acd117(83)=acd117(17)*acd117(43)
      acd117(84)=acd117(13)*acd117(42)
      acd117(79)=acd117(84)+acd117(83)+acd117(82)+acd117(81)+acd117(80)-acd117(&
      &47)+acd117(79)
      acd117(79)=acd117(11)*acd117(79)
      acd117(80)=acd117(9)*acd117(40)
      acd117(81)=acd117(27)*acd117(59)
      acd117(82)=acd117(23)*acd117(58)
      acd117(83)=acd117(19)*acd117(57)
      acd117(84)=acd117(7)*acd117(38)
      acd117(80)=acd117(84)+acd117(83)+acd117(82)+acd117(81)-acd117(60)+acd117(&
      &80)
      acd117(80)=acd117(17)*acd117(80)
      acd117(81)=acd117(5)*acd117(35)
      acd117(82)=acd117(29)*acd117(55)
      acd117(83)=acd117(25)*acd117(54)
      acd117(84)=acd117(21)*acd117(53)
      acd117(85)=acd117(3)*acd117(33)
      acd117(81)=acd117(85)+acd117(84)+acd117(83)+acd117(82)-acd117(56)+acd117(&
      &81)
      acd117(81)=acd117(15)*acd117(81)
      acd117(82)=acd117(27)*acd117(51)
      acd117(83)=acd117(23)*acd117(50)
      acd117(84)=acd117(19)*acd117(49)
      acd117(85)=acd117(7)*acd117(37)
      acd117(86)=acd117(15)*acd117(48)
      acd117(82)=acd117(86)+acd117(85)+acd117(84)+acd117(83)-acd117(52)+acd117(&
      &82)
      acd117(82)=acd117(13)*acd117(82)
      acd117(83)=-acd117(73)*acd117(76)
      acd117(84)=acd117(73)*acd117(74)
      acd117(84)=-acd117(75)+acd117(84)
      acd117(84)=acd117(72)*acd117(84)
      acd117(85)=-acd117(68)*acd117(71)
      acd117(86)=-acd117(68)*acd117(69)
      acd117(86)=-acd117(70)+acd117(86)
      acd117(86)=acd117(67)*acd117(86)
      acd117(87)=-acd117(9)*acd117(41)
      acd117(88)=-acd117(5)*acd117(36)
      acd117(89)=-acd117(29)*acd117(66)
      acd117(90)=-acd117(27)*acd117(65)
      acd117(91)=-acd117(25)*acd117(64)
      acd117(92)=-acd117(23)*acd117(63)
      acd117(93)=-acd117(21)*acd117(62)
      acd117(94)=-acd117(19)*acd117(61)
      acd117(95)=-acd117(7)*acd117(39)
      acd117(96)=-acd117(3)*acd117(34)
      brack=acd117(77)+acd117(78)+acd117(79)+acd117(80)+acd117(81)+acd117(82)+a&
      &cd117(83)+acd117(84)+acd117(85)+acd117(86)+acd117(87)+acd117(88)+acd117(&
      &89)+acd117(90)+acd117(91)+acd117(92)+acd117(93)+acd117(94)+acd117(95)+ac&
      &d117(96)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p2_part21part21_part25part25part21_model
      use p2_part21part21_part25part25part21_kinematics
      use p2_part21part21_part25part25part21_color
      use p2_part21part21_part25part25part21_abbrevd117h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(119) :: acd117
      complex(ki) :: brack
      acd117(1)=qshift(iv1)
      acd117(2)=dotproduct(qshift,qshift)
      acd117(3)=abb117(35)
      acd117(4)=dotproduct(qshift,spvak1k5)
      acd117(5)=abb117(34)
      acd117(6)=dotproduct(qshift,spvak2k5)
      acd117(7)=abb117(38)
      acd117(8)=dotproduct(qshift,spvak5k1)
      acd117(9)=abb117(41)
      acd117(10)=dotproduct(qshift,spvak5k2)
      acd117(11)=abb117(28)
      acd117(12)=dotproduct(qshift,spvak2e5)
      acd117(13)=abb117(10)
      acd117(14)=dotproduct(qshift,spvae5k2)
      acd117(15)=abb117(15)
      acd117(16)=dotproduct(qshift,spval3e5)
      acd117(17)=abb117(13)
      acd117(18)=dotproduct(qshift,spvae5l3)
      acd117(19)=abb117(60)
      acd117(20)=dotproduct(qshift,spval4e5)
      acd117(21)=abb117(18)
      acd117(22)=dotproduct(qshift,spvae5l4)
      acd117(23)=abb117(44)
      acd117(24)=dotproduct(qshift,spvak5e1)
      acd117(25)=abb117(72)
      acd117(26)=dotproduct(qshift,spvae1k5)
      acd117(27)=abb117(69)
      acd117(28)=dotproduct(qshift,spvak5e5)
      acd117(29)=abb117(66)
      acd117(30)=dotproduct(qshift,spvae5k5)
      acd117(31)=abb117(56)
      acd117(32)=abb117(22)
      acd117(33)=spvak1k5(iv1)
      acd117(34)=abb117(30)
      acd117(35)=abb117(42)
      acd117(36)=abb117(77)
      acd117(37)=spvak2k5(iv1)
      acd117(38)=abb117(40)
      acd117(39)=abb117(25)
      acd117(40)=spvak5k1(iv1)
      acd117(41)=abb117(39)
      acd117(42)=abb117(36)
      acd117(43)=abb117(74)
      acd117(44)=spvak5k2(iv1)
      acd117(45)=abb117(29)
      acd117(46)=abb117(20)
      acd117(47)=spvak2e5(iv1)
      acd117(48)=abb117(8)
      acd117(49)=abb117(17)
      acd117(50)=abb117(12)
      acd117(51)=abb117(26)
      acd117(52)=abb117(16)
      acd117(53)=abb117(9)
      acd117(54)=spvae5k2(iv1)
      acd117(55)=abb117(37)
      acd117(56)=abb117(31)
      acd117(57)=abb117(33)
      acd117(58)=abb117(24)
      acd117(59)=abb117(11)
      acd117(60)=spval3e5(iv1)
      acd117(61)=abb117(78)
      acd117(62)=abb117(80)
      acd117(63)=abb117(73)
      acd117(64)=abb117(7)
      acd117(65)=spvae5l3(iv1)
      acd117(66)=abb117(65)
      acd117(67)=abb117(76)
      acd117(68)=abb117(71)
      acd117(69)=abb117(51)
      acd117(70)=spval4e5(iv1)
      acd117(71)=abb117(54)
      acd117(72)=spvae5l4(iv1)
      acd117(73)=abb117(46)
      acd117(74)=spvak5e1(iv1)
      acd117(75)=abb117(70)
      acd117(76)=spvae1k5(iv1)
      acd117(77)=abb117(68)
      acd117(78)=spvak5e5(iv1)
      acd117(79)=abb117(57)
      acd117(80)=spvae5k5(iv1)
      acd117(81)=abb117(52)
      acd117(82)=spvak1e5(iv1)
      acd117(83)=dotproduct(qshift,spvae5k1)
      acd117(84)=abb117(19)
      acd117(85)=abb117(27)
      acd117(86)=spvae5k1(iv1)
      acd117(87)=dotproduct(qshift,spvak1e5)
      acd117(88)=abb117(23)
      acd117(89)=spvae1e5(iv1)
      acd117(90)=dotproduct(qshift,spvae5e1)
      acd117(91)=abb117(49)
      acd117(92)=abb117(43)
      acd117(93)=spvae5e1(iv1)
      acd117(94)=dotproduct(qshift,spvae1e5)
      acd117(95)=abb117(48)
      acd117(96)=acd117(10)*acd117(11)
      acd117(97)=acd117(6)*acd117(7)
      acd117(98)=acd117(30)*acd117(31)
      acd117(99)=acd117(28)*acd117(29)
      acd117(100)=acd117(26)*acd117(27)
      acd117(101)=acd117(24)*acd117(25)
      acd117(102)=acd117(22)*acd117(23)
      acd117(103)=acd117(20)*acd117(21)
      acd117(104)=acd117(8)*acd117(9)
      acd117(105)=acd117(4)*acd117(5)
      acd117(106)=acd117(18)*acd117(19)
      acd117(107)=acd117(16)*acd117(17)
      acd117(108)=acd117(14)*acd117(15)
      acd117(109)=acd117(12)*acd117(13)
      acd117(110)=acd117(2)*acd117(3)
      acd117(96)=-2.0_ki*acd117(110)+acd117(109)+acd117(108)+acd117(107)+acd117&
      &(106)+acd117(105)+acd117(104)+acd117(103)+acd117(102)+acd117(101)+acd117&
      &(100)+acd117(99)+acd117(98)+acd117(97)-acd117(32)+acd117(96)
      acd117(96)=acd117(1)*acd117(96)
      acd117(97)=acd117(44)*acd117(11)
      acd117(98)=acd117(37)*acd117(7)
      acd117(99)=acd117(80)*acd117(31)
      acd117(100)=acd117(78)*acd117(29)
      acd117(101)=acd117(76)*acd117(27)
      acd117(102)=acd117(74)*acd117(25)
      acd117(103)=acd117(72)*acd117(23)
      acd117(104)=acd117(70)*acd117(21)
      acd117(105)=acd117(40)*acd117(9)
      acd117(106)=acd117(33)*acd117(5)
      acd117(107)=acd117(65)*acd117(19)
      acd117(108)=acd117(60)*acd117(17)
      acd117(109)=acd117(54)*acd117(15)
      acd117(110)=acd117(47)*acd117(13)
      acd117(97)=acd117(110)+acd117(109)+acd117(108)+acd117(107)+acd117(106)+ac&
      &d117(105)+acd117(104)+acd117(103)+acd117(102)+acd117(101)+acd117(100)+ac&
      &d117(99)+acd117(97)+acd117(98)
      acd117(97)=acd117(2)*acd117(97)
      acd117(98)=-acd117(10)*acd117(45)
      acd117(99)=-acd117(28)*acd117(68)
      acd117(100)=-acd117(24)*acd117(67)
      acd117(101)=-acd117(20)*acd117(66)
      acd117(102)=-acd117(8)*acd117(42)
      acd117(103)=-acd117(12)*acd117(49)
      acd117(98)=acd117(103)+acd117(102)+acd117(101)+acd117(100)+acd117(99)+acd&
      &117(69)+acd117(98)
      acd117(98)=acd117(65)*acd117(98)
      acd117(99)=-acd117(6)*acd117(38)
      acd117(100)=-acd117(30)*acd117(63)
      acd117(101)=-acd117(26)*acd117(62)
      acd117(102)=-acd117(22)*acd117(61)
      acd117(103)=-acd117(4)*acd117(35)
      acd117(104)=-acd117(14)*acd117(55)
      acd117(99)=acd117(104)+acd117(103)+acd117(102)+acd117(101)+acd117(100)+ac&
      &d117(64)+acd117(99)
      acd117(99)=acd117(60)*acd117(99)
      acd117(100)=-acd117(28)*acd117(58)
      acd117(101)=-acd117(24)*acd117(57)
      acd117(102)=-acd117(20)*acd117(56)
      acd117(103)=-acd117(8)*acd117(41)
      acd117(104)=-acd117(16)*acd117(55)
      acd117(105)=-acd117(12)*acd117(48)
      acd117(100)=acd117(105)+acd117(104)+acd117(103)+acd117(102)+acd117(101)+a&
      &cd117(59)+acd117(100)
      acd117(100)=acd117(54)*acd117(100)
      acd117(101)=-acd117(30)*acd117(52)
      acd117(102)=-acd117(26)*acd117(51)
      acd117(103)=-acd117(22)*acd117(50)
      acd117(104)=-acd117(4)*acd117(34)
      acd117(105)=-acd117(18)*acd117(49)
      acd117(106)=-acd117(14)*acd117(48)
      acd117(101)=acd117(106)+acd117(105)+acd117(104)+acd117(103)+acd117(102)+a&
      &cd117(53)+acd117(101)
      acd117(101)=acd117(47)*acd117(101)
      acd117(102)=-acd117(44)*acd117(45)
      acd117(103)=-acd117(78)*acd117(68)
      acd117(104)=-acd117(74)*acd117(67)
      acd117(105)=-acd117(70)*acd117(66)
      acd117(106)=-acd117(40)*acd117(42)
      acd117(102)=acd117(106)+acd117(105)+acd117(104)+acd117(102)+acd117(103)
      acd117(102)=acd117(18)*acd117(102)
      acd117(103)=-acd117(37)*acd117(38)
      acd117(104)=-acd117(80)*acd117(63)
      acd117(105)=-acd117(76)*acd117(62)
      acd117(106)=-acd117(72)*acd117(61)
      acd117(107)=-acd117(33)*acd117(35)
      acd117(103)=acd117(107)+acd117(106)+acd117(105)+acd117(103)+acd117(104)
      acd117(103)=acd117(16)*acd117(103)
      acd117(104)=-acd117(78)*acd117(58)
      acd117(105)=-acd117(74)*acd117(57)
      acd117(106)=-acd117(70)*acd117(56)
      acd117(107)=-acd117(40)*acd117(41)
      acd117(104)=acd117(107)+acd117(106)+acd117(104)+acd117(105)
      acd117(104)=acd117(14)*acd117(104)
      acd117(105)=-acd117(80)*acd117(52)
      acd117(106)=-acd117(76)*acd117(51)
      acd117(107)=-acd117(72)*acd117(50)
      acd117(108)=-acd117(33)*acd117(34)
      acd117(105)=acd117(108)+acd117(107)+acd117(105)+acd117(106)
      acd117(105)=acd117(12)*acd117(105)
      acd117(106)=-acd117(91)*acd117(94)
      acd117(106)=acd117(95)+acd117(106)
      acd117(106)=acd117(93)*acd117(106)
      acd117(107)=acd117(84)*acd117(87)
      acd117(107)=acd117(107)+acd117(88)
      acd117(107)=acd117(86)*acd117(107)
      acd117(108)=-acd117(91)*acd117(90)
      acd117(108)=acd117(92)+acd117(108)
      acd117(108)=acd117(89)*acd117(108)
      acd117(109)=acd117(84)*acd117(83)
      acd117(109)=acd117(85)+acd117(109)
      acd117(109)=acd117(82)*acd117(109)
      acd117(110)=acd117(44)*acd117(46)
      acd117(111)=acd117(37)*acd117(39)
      acd117(112)=acd117(80)*acd117(81)
      acd117(113)=acd117(78)*acd117(79)
      acd117(114)=acd117(76)*acd117(77)
      acd117(115)=acd117(74)*acd117(75)
      acd117(116)=acd117(72)*acd117(73)
      acd117(117)=acd117(70)*acd117(71)
      acd117(118)=acd117(40)*acd117(43)
      acd117(119)=acd117(33)*acd117(36)
      brack=2.0_ki*acd117(96)+acd117(97)+acd117(98)+acd117(99)+acd117(100)+acd1&
      &17(101)+acd117(102)+acd117(103)+acd117(104)+acd117(105)+acd117(106)+acd1&
      &17(107)+acd117(108)+acd117(109)+acd117(110)+acd117(111)+acd117(112)+acd1&
      &17(113)+acd117(114)+acd117(115)+acd117(116)+acd117(117)+acd117(118)+acd1&
      &17(119)
   end function brack_2
!---#] function brack_2:
!---#[ function brack_3:
   pure function brack_3(Q, mu2) result(brack)
      use p2_part21part21_part25part25part21_model
      use p2_part21part21_part25part25part21_kinematics
      use p2_part21part21_part25part25part21_color
      use p2_part21part21_part25part25part21_abbrevd117h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(109) :: acd117
      complex(ki) :: brack
      acd117(1)=d(iv1,iv2)
      acd117(2)=dotproduct(qshift,qshift)
      acd117(3)=abb117(35)
      acd117(4)=dotproduct(qshift,spvak1k5)
      acd117(5)=abb117(34)
      acd117(6)=dotproduct(qshift,spvak2k5)
      acd117(7)=abb117(38)
      acd117(8)=dotproduct(qshift,spvak5k1)
      acd117(9)=abb117(41)
      acd117(10)=dotproduct(qshift,spvak5k2)
      acd117(11)=abb117(28)
      acd117(12)=dotproduct(qshift,spvak2e5)
      acd117(13)=abb117(10)
      acd117(14)=dotproduct(qshift,spvae5k2)
      acd117(15)=abb117(15)
      acd117(16)=dotproduct(qshift,spval3e5)
      acd117(17)=abb117(13)
      acd117(18)=dotproduct(qshift,spvae5l3)
      acd117(19)=abb117(60)
      acd117(20)=dotproduct(qshift,spval4e5)
      acd117(21)=abb117(18)
      acd117(22)=dotproduct(qshift,spvae5l4)
      acd117(23)=abb117(44)
      acd117(24)=dotproduct(qshift,spvak5e1)
      acd117(25)=abb117(72)
      acd117(26)=dotproduct(qshift,spvae1k5)
      acd117(27)=abb117(69)
      acd117(28)=dotproduct(qshift,spvak5e5)
      acd117(29)=abb117(66)
      acd117(30)=dotproduct(qshift,spvae5k5)
      acd117(31)=abb117(56)
      acd117(32)=abb117(22)
      acd117(33)=qshift(iv1)
      acd117(34)=qshift(iv2)
      acd117(35)=spvak1k5(iv2)
      acd117(36)=spvak2k5(iv2)
      acd117(37)=spvak5k1(iv2)
      acd117(38)=spvak5k2(iv2)
      acd117(39)=spvak2e5(iv2)
      acd117(40)=spvae5k2(iv2)
      acd117(41)=spval3e5(iv2)
      acd117(42)=spvae5l3(iv2)
      acd117(43)=spval4e5(iv2)
      acd117(44)=spvae5l4(iv2)
      acd117(45)=spvak5e1(iv2)
      acd117(46)=spvae1k5(iv2)
      acd117(47)=spvak5e5(iv2)
      acd117(48)=spvae5k5(iv2)
      acd117(49)=spvak1k5(iv1)
      acd117(50)=spvak2k5(iv1)
      acd117(51)=spvak5k1(iv1)
      acd117(52)=spvak5k2(iv1)
      acd117(53)=spvak2e5(iv1)
      acd117(54)=spvae5k2(iv1)
      acd117(55)=spval3e5(iv1)
      acd117(56)=spvae5l3(iv1)
      acd117(57)=spval4e5(iv1)
      acd117(58)=spvae5l4(iv1)
      acd117(59)=spvak5e1(iv1)
      acd117(60)=spvae1k5(iv1)
      acd117(61)=spvak5e5(iv1)
      acd117(62)=spvae5k5(iv1)
      acd117(63)=abb117(30)
      acd117(64)=abb117(42)
      acd117(65)=abb117(40)
      acd117(66)=abb117(39)
      acd117(67)=abb117(36)
      acd117(68)=abb117(29)
      acd117(69)=abb117(8)
      acd117(70)=abb117(17)
      acd117(71)=abb117(12)
      acd117(72)=abb117(26)
      acd117(73)=abb117(16)
      acd117(74)=abb117(37)
      acd117(75)=abb117(31)
      acd117(76)=abb117(33)
      acd117(77)=abb117(24)
      acd117(78)=abb117(78)
      acd117(79)=abb117(80)
      acd117(80)=abb117(73)
      acd117(81)=abb117(65)
      acd117(82)=abb117(76)
      acd117(83)=abb117(71)
      acd117(84)=spvak1e5(iv1)
      acd117(85)=spvae5k1(iv2)
      acd117(86)=abb117(19)
      acd117(87)=spvak1e5(iv2)
      acd117(88)=spvae5k1(iv1)
      acd117(89)=spvae1e5(iv1)
      acd117(90)=spvae5e1(iv2)
      acd117(91)=abb117(49)
      acd117(92)=spvae1e5(iv2)
      acd117(93)=spvae5e1(iv1)
      acd117(94)=-acd117(31)*acd117(48)
      acd117(95)=-acd117(29)*acd117(47)
      acd117(96)=-acd117(27)*acd117(46)
      acd117(97)=-acd117(25)*acd117(45)
      acd117(98)=-acd117(23)*acd117(44)
      acd117(99)=-acd117(21)*acd117(43)
      acd117(100)=-acd117(11)*acd117(38)
      acd117(101)=-acd117(9)*acd117(37)
      acd117(102)=-acd117(7)*acd117(36)
      acd117(103)=-acd117(5)*acd117(35)
      acd117(104)=-acd117(42)*acd117(19)
      acd117(105)=-acd117(41)*acd117(17)
      acd117(106)=-acd117(40)*acd117(15)
      acd117(107)=-acd117(39)*acd117(13)
      acd117(108)=acd117(34)*acd117(3)
      acd117(94)=4.0_ki*acd117(108)+acd117(107)+acd117(106)+acd117(105)+acd117(&
      &104)+acd117(103)+acd117(102)+acd117(101)+acd117(100)+acd117(99)+acd117(9&
      &8)+acd117(97)+acd117(96)+acd117(94)+acd117(95)
      acd117(94)=acd117(33)*acd117(94)
      acd117(95)=acd117(3)*acd117(2)
      acd117(96)=-acd117(31)*acd117(30)
      acd117(97)=-acd117(29)*acd117(28)
      acd117(98)=-acd117(27)*acd117(26)
      acd117(99)=-acd117(25)*acd117(24)
      acd117(100)=-acd117(23)*acd117(22)
      acd117(101)=-acd117(21)*acd117(20)
      acd117(102)=-acd117(19)*acd117(18)
      acd117(103)=-acd117(17)*acd117(16)
      acd117(104)=-acd117(15)*acd117(14)
      acd117(105)=-acd117(13)*acd117(12)
      acd117(106)=-acd117(11)*acd117(10)
      acd117(107)=-acd117(9)*acd117(8)
      acd117(108)=-acd117(7)*acd117(6)
      acd117(109)=-acd117(5)*acd117(4)
      acd117(95)=acd117(109)+acd117(108)+acd117(107)+acd117(106)+acd117(105)+ac&
      &d117(104)+acd117(103)+acd117(102)+acd117(101)+acd117(100)+acd117(99)+acd&
      &117(98)+acd117(97)+acd117(96)+acd117(32)+2.0_ki*acd117(95)
      acd117(95)=acd117(1)*acd117(95)
      acd117(96)=-acd117(31)*acd117(62)
      acd117(97)=-acd117(29)*acd117(61)
      acd117(98)=-acd117(27)*acd117(60)
      acd117(99)=-acd117(25)*acd117(59)
      acd117(100)=-acd117(23)*acd117(58)
      acd117(101)=-acd117(21)*acd117(57)
      acd117(102)=-acd117(11)*acd117(52)
      acd117(103)=-acd117(9)*acd117(51)
      acd117(104)=-acd117(7)*acd117(50)
      acd117(105)=-acd117(5)*acd117(49)
      acd117(106)=-acd117(56)*acd117(19)
      acd117(107)=-acd117(55)*acd117(17)
      acd117(108)=-acd117(54)*acd117(15)
      acd117(109)=-acd117(53)*acd117(13)
      acd117(96)=acd117(109)+acd117(108)+acd117(107)+acd117(106)+acd117(105)+ac&
      &d117(104)+acd117(103)+acd117(102)+acd117(101)+acd117(100)+acd117(99)+acd&
      &117(98)+acd117(96)+acd117(97)
      acd117(96)=acd117(34)*acd117(96)
      acd117(94)=acd117(95)+acd117(96)+acd117(94)
      acd117(95)=acd117(52)*acd117(68)
      acd117(96)=acd117(61)*acd117(83)
      acd117(97)=acd117(59)*acd117(82)
      acd117(98)=acd117(57)*acd117(81)
      acd117(99)=acd117(51)*acd117(67)
      acd117(100)=acd117(53)*acd117(70)
      acd117(95)=acd117(100)+acd117(99)+acd117(98)+acd117(97)+acd117(95)+acd117&
      &(96)
      acd117(95)=acd117(42)*acd117(95)
      acd117(96)=acd117(50)*acd117(65)
      acd117(97)=acd117(62)*acd117(80)
      acd117(98)=acd117(60)*acd117(79)
      acd117(99)=acd117(58)*acd117(78)
      acd117(100)=acd117(49)*acd117(64)
      acd117(101)=acd117(54)*acd117(74)
      acd117(96)=acd117(101)+acd117(100)+acd117(99)+acd117(98)+acd117(96)+acd11&
      &7(97)
      acd117(96)=acd117(41)*acd117(96)
      acd117(97)=acd117(61)*acd117(77)
      acd117(98)=acd117(59)*acd117(76)
      acd117(99)=acd117(57)*acd117(75)
      acd117(100)=acd117(51)*acd117(66)
      acd117(101)=acd117(55)*acd117(74)
      acd117(102)=acd117(53)*acd117(69)
      acd117(97)=acd117(102)+acd117(101)+acd117(100)+acd117(99)+acd117(97)+acd1&
      &17(98)
      acd117(97)=acd117(40)*acd117(97)
      acd117(98)=acd117(62)*acd117(73)
      acd117(99)=acd117(60)*acd117(72)
      acd117(100)=acd117(58)*acd117(71)
      acd117(101)=acd117(49)*acd117(63)
      acd117(102)=acd117(56)*acd117(70)
      acd117(103)=acd117(54)*acd117(69)
      acd117(98)=acd117(103)+acd117(102)+acd117(101)+acd117(100)+acd117(98)+acd&
      &117(99)
      acd117(98)=acd117(39)*acd117(98)
      acd117(99)=acd117(38)*acd117(68)
      acd117(100)=acd117(47)*acd117(83)
      acd117(101)=acd117(45)*acd117(82)
      acd117(102)=acd117(43)*acd117(81)
      acd117(103)=acd117(37)*acd117(67)
      acd117(99)=acd117(103)+acd117(102)+acd117(101)+acd117(99)+acd117(100)
      acd117(99)=acd117(56)*acd117(99)
      acd117(100)=acd117(36)*acd117(65)
      acd117(101)=acd117(48)*acd117(80)
      acd117(102)=acd117(46)*acd117(79)
      acd117(103)=acd117(44)*acd117(78)
      acd117(104)=acd117(35)*acd117(64)
      acd117(100)=acd117(104)+acd117(103)+acd117(102)+acd117(100)+acd117(101)
      acd117(100)=acd117(55)*acd117(100)
      acd117(101)=acd117(47)*acd117(77)
      acd117(102)=acd117(45)*acd117(76)
      acd117(103)=acd117(43)*acd117(75)
      acd117(104)=acd117(37)*acd117(66)
      acd117(101)=acd117(104)+acd117(103)+acd117(101)+acd117(102)
      acd117(101)=acd117(54)*acd117(101)
      acd117(102)=acd117(48)*acd117(73)
      acd117(103)=acd117(46)*acd117(72)
      acd117(104)=acd117(44)*acd117(71)
      acd117(105)=acd117(35)*acd117(63)
      acd117(102)=acd117(105)+acd117(104)+acd117(102)+acd117(103)
      acd117(102)=acd117(53)*acd117(102)
      acd117(103)=acd117(92)*acd117(93)
      acd117(104)=acd117(89)*acd117(90)
      acd117(103)=acd117(103)+acd117(104)
      acd117(103)=acd117(91)*acd117(103)
      acd117(104)=-acd117(87)*acd117(88)
      acd117(105)=-acd117(84)*acd117(85)
      acd117(104)=acd117(104)+acd117(105)
      acd117(104)=acd117(86)*acd117(104)
      brack=2.0_ki*acd117(94)+acd117(95)+acd117(96)+acd117(97)+acd117(98)+acd11&
      &7(99)+acd117(100)+acd117(101)+acd117(102)+acd117(103)+acd117(104)
   end function brack_3
!---#] function brack_3:
!---#[ function brack_4:
   pure function brack_4(Q, mu2) result(brack)
      use p2_part21part21_part25part25part21_model
      use p2_part21part21_part25part25part21_kinematics
      use p2_part21part21_part25part25part21_color
      use p2_part21part21_part25part25part21_abbrevd117h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(80) :: acd117
      complex(ki) :: brack
      acd117(1)=d(iv1,iv2)
      acd117(2)=qshift(iv3)
      acd117(3)=abb117(35)
      acd117(4)=spvak1k5(iv3)
      acd117(5)=abb117(34)
      acd117(6)=spvak2k5(iv3)
      acd117(7)=abb117(38)
      acd117(8)=spvak5k1(iv3)
      acd117(9)=abb117(41)
      acd117(10)=spvak5k2(iv3)
      acd117(11)=abb117(28)
      acd117(12)=spvak2e5(iv3)
      acd117(13)=abb117(10)
      acd117(14)=spvae5k2(iv3)
      acd117(15)=abb117(15)
      acd117(16)=spval3e5(iv3)
      acd117(17)=abb117(13)
      acd117(18)=spvae5l3(iv3)
      acd117(19)=abb117(60)
      acd117(20)=spval4e5(iv3)
      acd117(21)=abb117(18)
      acd117(22)=spvae5l4(iv3)
      acd117(23)=abb117(44)
      acd117(24)=spvak5e1(iv3)
      acd117(25)=abb117(72)
      acd117(26)=spvae1k5(iv3)
      acd117(27)=abb117(69)
      acd117(28)=spvak5e5(iv3)
      acd117(29)=abb117(66)
      acd117(30)=spvae5k5(iv3)
      acd117(31)=abb117(56)
      acd117(32)=d(iv1,iv3)
      acd117(33)=qshift(iv2)
      acd117(34)=spvak1k5(iv2)
      acd117(35)=spvak2k5(iv2)
      acd117(36)=spvak5k1(iv2)
      acd117(37)=spvak5k2(iv2)
      acd117(38)=spvak2e5(iv2)
      acd117(39)=spvae5k2(iv2)
      acd117(40)=spval3e5(iv2)
      acd117(41)=spvae5l3(iv2)
      acd117(42)=spval4e5(iv2)
      acd117(43)=spvae5l4(iv2)
      acd117(44)=spvak5e1(iv2)
      acd117(45)=spvae1k5(iv2)
      acd117(46)=spvak5e5(iv2)
      acd117(47)=spvae5k5(iv2)
      acd117(48)=d(iv2,iv3)
      acd117(49)=qshift(iv1)
      acd117(50)=spvak1k5(iv1)
      acd117(51)=spvak2k5(iv1)
      acd117(52)=spvak5k1(iv1)
      acd117(53)=spvak5k2(iv1)
      acd117(54)=spvak2e5(iv1)
      acd117(55)=spvae5k2(iv1)
      acd117(56)=spval3e5(iv1)
      acd117(57)=spvae5l3(iv1)
      acd117(58)=spval4e5(iv1)
      acd117(59)=spvae5l4(iv1)
      acd117(60)=spvak5e1(iv1)
      acd117(61)=spvae1k5(iv1)
      acd117(62)=spvak5e5(iv1)
      acd117(63)=spvae5k5(iv1)
      acd117(64)=4.0_ki*acd117(3)
      acd117(65)=-acd117(2)*acd117(64)
      acd117(66)=acd117(4)*acd117(5)
      acd117(67)=acd117(6)*acd117(7)
      acd117(68)=acd117(8)*acd117(9)
      acd117(69)=acd117(10)*acd117(11)
      acd117(70)=acd117(12)*acd117(13)
      acd117(71)=acd117(14)*acd117(15)
      acd117(72)=acd117(16)*acd117(17)
      acd117(73)=acd117(18)*acd117(19)
      acd117(74)=acd117(20)*acd117(21)
      acd117(75)=acd117(22)*acd117(23)
      acd117(76)=acd117(24)*acd117(25)
      acd117(77)=acd117(26)*acd117(27)
      acd117(78)=acd117(28)*acd117(29)
      acd117(79)=acd117(30)*acd117(31)
      acd117(65)=acd117(79)+acd117(78)+acd117(77)+acd117(76)+acd117(75)+acd117(&
      &74)+acd117(73)+acd117(72)+acd117(71)+acd117(70)+acd117(69)+acd117(68)+ac&
      &d117(67)+acd117(65)+acd117(66)
      acd117(65)=acd117(1)*acd117(65)
      acd117(66)=-acd117(33)*acd117(64)
      acd117(67)=acd117(34)*acd117(5)
      acd117(68)=acd117(35)*acd117(7)
      acd117(69)=acd117(36)*acd117(9)
      acd117(70)=acd117(37)*acd117(11)
      acd117(71)=acd117(38)*acd117(13)
      acd117(72)=acd117(39)*acd117(15)
      acd117(73)=acd117(40)*acd117(17)
      acd117(74)=acd117(41)*acd117(19)
      acd117(75)=acd117(42)*acd117(21)
      acd117(76)=acd117(43)*acd117(23)
      acd117(77)=acd117(44)*acd117(25)
      acd117(78)=acd117(45)*acd117(27)
      acd117(79)=acd117(46)*acd117(29)
      acd117(80)=acd117(47)*acd117(31)
      acd117(66)=acd117(80)+acd117(79)+acd117(78)+acd117(77)+acd117(76)+acd117(&
      &75)+acd117(74)+acd117(73)+acd117(72)+acd117(71)+acd117(70)+acd117(69)+ac&
      &d117(68)+acd117(67)+acd117(66)
      acd117(66)=acd117(32)*acd117(66)
      acd117(64)=-acd117(49)*acd117(64)
      acd117(67)=acd117(50)*acd117(5)
      acd117(68)=acd117(51)*acd117(7)
      acd117(69)=acd117(52)*acd117(9)
      acd117(70)=acd117(53)*acd117(11)
      acd117(71)=acd117(54)*acd117(13)
      acd117(72)=acd117(55)*acd117(15)
      acd117(73)=acd117(56)*acd117(17)
      acd117(74)=acd117(57)*acd117(19)
      acd117(75)=acd117(58)*acd117(21)
      acd117(76)=acd117(59)*acd117(23)
      acd117(77)=acd117(60)*acd117(25)
      acd117(78)=acd117(61)*acd117(27)
      acd117(79)=acd117(62)*acd117(29)
      acd117(80)=acd117(63)*acd117(31)
      acd117(64)=acd117(80)+acd117(79)+acd117(78)+acd117(77)+acd117(76)+acd117(&
      &75)+acd117(74)+acd117(73)+acd117(72)+acd117(71)+acd117(70)+acd117(69)+ac&
      &d117(68)+acd117(67)+acd117(64)
      acd117(64)=acd117(48)*acd117(64)
      acd117(64)=acd117(64)+acd117(66)+acd117(65)
      brack=2.0_ki*acd117(64)
   end function brack_4
!---#] function brack_4:
!---#[ function brack_5:
   pure function brack_5(Q, mu2) result(brack)
      use p2_part21part21_part25part25part21_model
      use p2_part21part21_part25part25part21_kinematics
      use p2_part21part21_part25part25part21_color
      use p2_part21part21_part25part25part21_abbrevd117h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(10) :: acd117
      complex(ki) :: brack
      acd117(1)=d(iv1,iv2)
      acd117(2)=d(iv3,iv4)
      acd117(3)=abb117(35)
      acd117(4)=d(iv1,iv3)
      acd117(5)=d(iv2,iv4)
      acd117(6)=d(iv1,iv4)
      acd117(7)=d(iv2,iv3)
      acd117(8)=acd117(2)*acd117(1)
      acd117(9)=acd117(5)*acd117(4)
      acd117(10)=acd117(7)*acd117(6)
      acd117(8)=acd117(10)+acd117(8)+acd117(9)
      brack=8.0_ki*acd117(8)*acd117(3)
   end function brack_5
!---#] function brack_5:
!---#[ function derivative:
   function derivative(mu2,i1,i2,i3,i4) result(numerator)
      use p2_part21part21_part25part25part21_globalsl1, only: epspow
      use p2_part21part21_part25part25part21_kinematics
      use p2_part21part21_part25part25part21_abbrevd117h0
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
      qshift = k5+k4
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
end module     p2_part21part21_part25part25part21_d117h0l1d
