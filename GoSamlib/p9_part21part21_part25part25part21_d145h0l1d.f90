module     p9_part21part21_part25part25part21_d145h0l1d
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p9_part21part &
   ! &21_part25part25part21/helicity0d145h0l1d.f90
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
   integer, private :: iv4
   real(ki), dimension(4), private :: qshift
   public :: derivative
contains
!---#[ function brack_1:
   pure function brack_1(Q, mu2) result(brack)
      use p9_part21part21_part25part25part21_model
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_color
      use p9_part21part21_part25part25part21_abbrevd145h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(61) :: acd145
      complex(ki) :: brack
      acd145(1)=dotproduct(qshift,qshift)
      acd145(2)=dotproduct(qshift,spvae1e2)
      acd145(3)=abb145(25)
      acd145(4)=dotproduct(qshift,spvae2e1)
      acd145(5)=abb145(19)
      acd145(6)=dotproduct(qshift,spvae1e5)
      acd145(7)=abb145(27)
      acd145(8)=dotproduct(qshift,spvae5e1)
      acd145(9)=abb145(22)
      acd145(10)=abb145(13)
      acd145(11)=dotproduct(qshift,spvae2e5)
      acd145(12)=abb145(7)
      acd145(13)=abb145(6)
      acd145(14)=dotproduct(qshift,spvae2k1)
      acd145(15)=abb145(24)
      acd145(16)=dotproduct(qshift,spvae5k2)
      acd145(17)=abb145(20)
      acd145(18)=dotproduct(qshift,spvae2k5)
      acd145(19)=abb145(31)
      acd145(20)=abb145(23)
      acd145(21)=abb145(11)
      acd145(22)=dotproduct(qshift,spvae5e2)
      acd145(23)=abb145(16)
      acd145(24)=dotproduct(qshift,spvak1e2)
      acd145(25)=abb145(4)
      acd145(26)=dotproduct(qshift,spvak2e5)
      acd145(27)=abb145(37)
      acd145(28)=dotproduct(qshift,spvak5e2)
      acd145(29)=abb145(34)
      acd145(30)=abb145(15)
      acd145(31)=abb145(3)
      acd145(32)=abb145(26)
      acd145(33)=abb145(8)
      acd145(34)=abb145(30)
      acd145(35)=abb145(9)
      acd145(36)=abb145(21)
      acd145(37)=abb145(36)
      acd145(38)=abb145(33)
      acd145(39)=abb145(18)
      acd145(40)=abb145(12)
      acd145(41)=abb145(28)
      acd145(42)=abb145(35)
      acd145(43)=abb145(10)
      acd145(44)=abb145(32)
      acd145(45)=abb145(29)
      acd145(46)=abb145(17)
      acd145(47)=abb145(14)
      acd145(48)=abb145(5)
      acd145(49)=acd145(28)*acd145(29)
      acd145(50)=acd145(26)*acd145(27)
      acd145(51)=acd145(24)*acd145(25)
      acd145(52)=-acd145(22)*acd145(30)
      acd145(53)=-acd145(1)*acd145(5)
      acd145(54)=acd145(22)*acd145(12)
      acd145(54)=acd145(23)+acd145(54)
      acd145(54)=acd145(6)*acd145(54)
      acd145(49)=acd145(54)+acd145(53)+acd145(52)+acd145(51)+acd145(50)-acd145(&
      &31)+acd145(49)
      acd145(49)=acd145(4)*acd145(49)
      acd145(50)=acd145(18)*acd145(19)
      acd145(51)=acd145(16)*acd145(17)
      acd145(52)=acd145(14)*acd145(15)
      acd145(53)=-acd145(11)*acd145(20)
      acd145(54)=-acd145(1)*acd145(3)
      acd145(55)=acd145(11)*acd145(12)
      acd145(55)=acd145(13)+acd145(55)
      acd145(55)=acd145(8)*acd145(55)
      acd145(50)=acd145(55)+acd145(54)+acd145(53)+acd145(52)+acd145(51)-acd145(&
      &21)+acd145(50)
      acd145(50)=acd145(2)*acd145(50)
      acd145(51)=acd145(28)*acd145(38)
      acd145(52)=acd145(26)*acd145(37)
      acd145(53)=acd145(24)*acd145(36)
      acd145(54)=-acd145(1)*acd145(9)
      acd145(51)=acd145(54)+acd145(53)+acd145(52)-acd145(39)+acd145(51)
      acd145(51)=acd145(8)*acd145(51)
      acd145(52)=-acd145(18)*acd145(34)
      acd145(53)=acd145(16)*acd145(33)
      acd145(54)=acd145(14)*acd145(32)
      acd145(55)=-acd145(1)*acd145(7)
      acd145(52)=acd145(55)+acd145(54)+acd145(53)-acd145(35)+acd145(52)
      acd145(52)=acd145(6)*acd145(52)
      acd145(53)=-acd145(28)*acd145(44)
      acd145(54)=-acd145(26)*acd145(42)
      acd145(55)=-acd145(24)*acd145(40)
      acd145(56)=-acd145(22)*acd145(47)
      acd145(57)=-acd145(18)*acd145(45)
      acd145(58)=-acd145(16)*acd145(43)
      acd145(59)=-acd145(14)*acd145(41)
      acd145(60)=-acd145(11)*acd145(46)
      acd145(61)=acd145(1)*acd145(10)
      brack=acd145(48)+acd145(49)+acd145(50)+acd145(51)+acd145(52)+acd145(53)+a&
      &cd145(54)+acd145(55)+acd145(56)+acd145(57)+acd145(58)+acd145(59)+acd145(&
      &60)+acd145(61)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p9_part21part21_part25part25part21_model
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_color
      use p9_part21part21_part25part25part21_abbrevd145h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(77) :: acd145
      complex(ki) :: brack
      acd145(1)=qshift(iv1)
      acd145(2)=dotproduct(qshift,spvae1e2)
      acd145(3)=abb145(25)
      acd145(4)=dotproduct(qshift,spvae2e1)
      acd145(5)=abb145(19)
      acd145(6)=dotproduct(qshift,spvae1e5)
      acd145(7)=abb145(27)
      acd145(8)=dotproduct(qshift,spvae5e1)
      acd145(9)=abb145(22)
      acd145(10)=abb145(13)
      acd145(11)=spvae1e2(iv1)
      acd145(12)=dotproduct(qshift,qshift)
      acd145(13)=dotproduct(qshift,spvae2e5)
      acd145(14)=abb145(7)
      acd145(15)=abb145(6)
      acd145(16)=dotproduct(qshift,spvae2k1)
      acd145(17)=abb145(24)
      acd145(18)=dotproduct(qshift,spvae5k2)
      acd145(19)=abb145(20)
      acd145(20)=dotproduct(qshift,spvae2k5)
      acd145(21)=abb145(31)
      acd145(22)=abb145(23)
      acd145(23)=abb145(11)
      acd145(24)=spvae2e1(iv1)
      acd145(25)=dotproduct(qshift,spvae5e2)
      acd145(26)=abb145(16)
      acd145(27)=dotproduct(qshift,spvak1e2)
      acd145(28)=abb145(4)
      acd145(29)=dotproduct(qshift,spvak2e5)
      acd145(30)=abb145(37)
      acd145(31)=dotproduct(qshift,spvak5e2)
      acd145(32)=abb145(34)
      acd145(33)=abb145(15)
      acd145(34)=abb145(3)
      acd145(35)=spvae1e5(iv1)
      acd145(36)=abb145(26)
      acd145(37)=abb145(8)
      acd145(38)=abb145(30)
      acd145(39)=abb145(9)
      acd145(40)=spvae5e1(iv1)
      acd145(41)=abb145(21)
      acd145(42)=abb145(36)
      acd145(43)=abb145(33)
      acd145(44)=abb145(18)
      acd145(45)=spvak1e2(iv1)
      acd145(46)=abb145(12)
      acd145(47)=spvae2k1(iv1)
      acd145(48)=abb145(28)
      acd145(49)=spvak2e5(iv1)
      acd145(50)=abb145(35)
      acd145(51)=spvae5k2(iv1)
      acd145(52)=abb145(10)
      acd145(53)=spvak5e2(iv1)
      acd145(54)=abb145(32)
      acd145(55)=spvae2k5(iv1)
      acd145(56)=abb145(29)
      acd145(57)=spvae2e5(iv1)
      acd145(58)=abb145(17)
      acd145(59)=spvae5e2(iv1)
      acd145(60)=abb145(14)
      acd145(61)=-acd145(31)*acd145(32)
      acd145(62)=-acd145(29)*acd145(30)
      acd145(63)=-acd145(27)*acd145(28)
      acd145(64)=acd145(25)*acd145(33)
      acd145(65)=acd145(12)*acd145(5)
      acd145(66)=acd145(14)*acd145(25)
      acd145(66)=acd145(66)+acd145(26)
      acd145(67)=-acd145(6)*acd145(66)
      acd145(61)=acd145(67)+acd145(65)+acd145(64)+acd145(63)+acd145(62)+acd145(&
      &34)+acd145(61)
      acd145(61)=acd145(24)*acd145(61)
      acd145(62)=-acd145(20)*acd145(21)
      acd145(63)=-acd145(18)*acd145(19)
      acd145(64)=-acd145(16)*acd145(17)
      acd145(65)=acd145(13)*acd145(22)
      acd145(67)=acd145(12)*acd145(3)
      acd145(68)=acd145(14)*acd145(13)
      acd145(68)=acd145(68)+acd145(15)
      acd145(69)=-acd145(8)*acd145(68)
      acd145(62)=acd145(69)+acd145(67)+acd145(65)+acd145(64)+acd145(63)+acd145(&
      &23)+acd145(62)
      acd145(62)=acd145(11)*acd145(62)
      acd145(63)=-acd145(6)*acd145(14)
      acd145(63)=acd145(63)+acd145(33)
      acd145(63)=acd145(59)*acd145(63)
      acd145(64)=-acd145(53)*acd145(32)
      acd145(65)=-acd145(49)*acd145(30)
      acd145(67)=-acd145(45)*acd145(28)
      acd145(69)=2.0_ki*acd145(1)
      acd145(70)=acd145(5)*acd145(69)
      acd145(66)=-acd145(35)*acd145(66)
      acd145(63)=acd145(66)+acd145(70)+acd145(67)+acd145(65)+acd145(64)+acd145(&
      &63)
      acd145(63)=acd145(4)*acd145(63)
      acd145(64)=-acd145(8)*acd145(14)
      acd145(64)=acd145(64)+acd145(22)
      acd145(64)=acd145(57)*acd145(64)
      acd145(65)=-acd145(55)*acd145(21)
      acd145(66)=-acd145(51)*acd145(19)
      acd145(67)=-acd145(47)*acd145(17)
      acd145(70)=acd145(3)*acd145(69)
      acd145(68)=-acd145(40)*acd145(68)
      acd145(64)=acd145(68)+acd145(70)+acd145(67)+acd145(66)+acd145(65)+acd145(&
      &64)
      acd145(64)=acd145(2)*acd145(64)
      acd145(65)=-acd145(31)*acd145(43)
      acd145(66)=-acd145(29)*acd145(42)
      acd145(67)=-acd145(27)*acd145(41)
      acd145(68)=acd145(12)*acd145(9)
      acd145(65)=acd145(68)+acd145(67)+acd145(66)+acd145(44)+acd145(65)
      acd145(65)=acd145(40)*acd145(65)
      acd145(66)=acd145(20)*acd145(38)
      acd145(67)=-acd145(18)*acd145(37)
      acd145(68)=-acd145(16)*acd145(36)
      acd145(70)=acd145(12)*acd145(7)
      acd145(66)=acd145(70)+acd145(68)+acd145(67)+acd145(39)+acd145(66)
      acd145(66)=acd145(35)*acd145(66)
      acd145(67)=-acd145(53)*acd145(43)
      acd145(68)=-acd145(49)*acd145(42)
      acd145(70)=-acd145(45)*acd145(41)
      acd145(71)=acd145(9)*acd145(69)
      acd145(67)=acd145(71)+acd145(70)+acd145(67)+acd145(68)
      acd145(67)=acd145(8)*acd145(67)
      acd145(68)=acd145(55)*acd145(38)
      acd145(70)=-acd145(51)*acd145(37)
      acd145(71)=-acd145(47)*acd145(36)
      acd145(72)=acd145(7)*acd145(69)
      acd145(68)=acd145(72)+acd145(71)+acd145(68)+acd145(70)
      acd145(68)=acd145(6)*acd145(68)
      acd145(70)=acd145(59)*acd145(60)
      acd145(71)=acd145(57)*acd145(58)
      acd145(72)=acd145(55)*acd145(56)
      acd145(73)=acd145(53)*acd145(54)
      acd145(74)=acd145(51)*acd145(52)
      acd145(75)=acd145(49)*acd145(50)
      acd145(76)=acd145(47)*acd145(48)
      acd145(77)=acd145(45)*acd145(46)
      acd145(69)=-acd145(10)*acd145(69)
      brack=acd145(61)+acd145(62)+acd145(63)+acd145(64)+acd145(65)+acd145(66)+a&
      &cd145(67)+acd145(68)+acd145(69)+acd145(70)+acd145(71)+acd145(72)+acd145(&
      &73)+acd145(74)+acd145(75)+acd145(76)+acd145(77)
   end function brack_2
!---#] function brack_2:
!---#[ function brack_3:
   pure function brack_3(Q, mu2) result(brack)
      use p9_part21part21_part25part25part21_model
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_color
      use p9_part21part21_part25part25part21_abbrevd145h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(68) :: acd145
      complex(ki) :: brack
      acd145(1)=d(iv1,iv2)
      acd145(2)=dotproduct(qshift,spvae1e2)
      acd145(3)=abb145(25)
      acd145(4)=dotproduct(qshift,spvae2e1)
      acd145(5)=abb145(19)
      acd145(6)=dotproduct(qshift,spvae1e5)
      acd145(7)=abb145(27)
      acd145(8)=dotproduct(qshift,spvae5e1)
      acd145(9)=abb145(22)
      acd145(10)=abb145(13)
      acd145(11)=qshift(iv1)
      acd145(12)=spvae1e2(iv2)
      acd145(13)=spvae2e1(iv2)
      acd145(14)=spvae1e5(iv2)
      acd145(15)=spvae5e1(iv2)
      acd145(16)=qshift(iv2)
      acd145(17)=spvae1e2(iv1)
      acd145(18)=spvae2e1(iv1)
      acd145(19)=spvae1e5(iv1)
      acd145(20)=spvae5e1(iv1)
      acd145(21)=dotproduct(qshift,spvae2e5)
      acd145(22)=abb145(7)
      acd145(23)=abb145(6)
      acd145(24)=spvae2k1(iv2)
      acd145(25)=abb145(24)
      acd145(26)=spvae5k2(iv2)
      acd145(27)=abb145(20)
      acd145(28)=spvae2k5(iv2)
      acd145(29)=abb145(31)
      acd145(30)=spvae2e5(iv2)
      acd145(31)=abb145(23)
      acd145(32)=spvae2k1(iv1)
      acd145(33)=spvae5k2(iv1)
      acd145(34)=spvae2k5(iv1)
      acd145(35)=spvae2e5(iv1)
      acd145(36)=dotproduct(qshift,spvae5e2)
      acd145(37)=abb145(16)
      acd145(38)=spvak1e2(iv2)
      acd145(39)=abb145(4)
      acd145(40)=spvak2e5(iv2)
      acd145(41)=abb145(37)
      acd145(42)=spvak5e2(iv2)
      acd145(43)=abb145(34)
      acd145(44)=spvae5e2(iv2)
      acd145(45)=abb145(15)
      acd145(46)=spvak1e2(iv1)
      acd145(47)=spvak2e5(iv1)
      acd145(48)=spvak5e2(iv1)
      acd145(49)=spvae5e2(iv1)
      acd145(50)=abb145(26)
      acd145(51)=abb145(8)
      acd145(52)=abb145(30)
      acd145(53)=abb145(21)
      acd145(54)=abb145(36)
      acd145(55)=abb145(33)
      acd145(56)=acd145(19)*acd145(44)
      acd145(57)=acd145(14)*acd145(49)
      acd145(56)=acd145(57)+acd145(56)
      acd145(56)=acd145(4)*acd145(56)
      acd145(57)=acd145(20)*acd145(30)
      acd145(58)=acd145(15)*acd145(35)
      acd145(57)=acd145(58)+acd145(57)
      acd145(57)=acd145(2)*acd145(57)
      acd145(58)=acd145(6)*acd145(44)
      acd145(59)=acd145(14)*acd145(36)
      acd145(58)=acd145(58)+acd145(59)
      acd145(58)=acd145(18)*acd145(58)
      acd145(59)=acd145(8)*acd145(30)
      acd145(60)=acd145(15)*acd145(21)
      acd145(59)=acd145(59)+acd145(60)
      acd145(59)=acd145(17)*acd145(59)
      acd145(60)=acd145(6)*acd145(49)
      acd145(61)=acd145(19)*acd145(36)
      acd145(60)=acd145(60)+acd145(61)
      acd145(60)=acd145(13)*acd145(60)
      acd145(61)=acd145(8)*acd145(35)
      acd145(62)=acd145(20)*acd145(21)
      acd145(61)=acd145(61)+acd145(62)
      acd145(61)=acd145(12)*acd145(61)
      acd145(56)=acd145(61)+acd145(60)+acd145(59)+acd145(58)+acd145(57)+acd145(&
      &56)
      acd145(56)=acd145(22)*acd145(56)
      acd145(57)=acd145(42)*acd145(43)
      acd145(58)=acd145(40)*acd145(41)
      acd145(59)=acd145(38)*acd145(39)
      acd145(60)=-acd145(44)*acd145(45)
      acd145(61)=2.0_ki*acd145(16)
      acd145(62)=-acd145(5)*acd145(61)
      acd145(63)=acd145(14)*acd145(37)
      acd145(57)=acd145(63)+acd145(62)+acd145(60)+acd145(59)+acd145(57)+acd145(&
      &58)
      acd145(57)=acd145(18)*acd145(57)
      acd145(58)=acd145(28)*acd145(29)
      acd145(59)=acd145(26)*acd145(27)
      acd145(60)=acd145(24)*acd145(25)
      acd145(62)=-acd145(30)*acd145(31)
      acd145(63)=-acd145(3)*acd145(61)
      acd145(64)=acd145(15)*acd145(23)
      acd145(58)=acd145(64)+acd145(63)+acd145(62)+acd145(60)+acd145(58)+acd145(&
      &59)
      acd145(58)=acd145(17)*acd145(58)
      acd145(59)=acd145(43)*acd145(48)
      acd145(60)=acd145(41)*acd145(47)
      acd145(62)=acd145(39)*acd145(46)
      acd145(63)=-acd145(49)*acd145(45)
      acd145(64)=2.0_ki*acd145(11)
      acd145(65)=-acd145(5)*acd145(64)
      acd145(66)=acd145(19)*acd145(37)
      acd145(59)=acd145(66)+acd145(65)+acd145(63)+acd145(62)+acd145(59)+acd145(&
      &60)
      acd145(59)=acd145(13)*acd145(59)
      acd145(60)=acd145(29)*acd145(34)
      acd145(62)=acd145(27)*acd145(33)
      acd145(63)=acd145(25)*acd145(32)
      acd145(65)=-acd145(35)*acd145(31)
      acd145(66)=-acd145(3)*acd145(64)
      acd145(67)=acd145(20)*acd145(23)
      acd145(60)=acd145(67)+acd145(66)+acd145(65)+acd145(63)+acd145(60)+acd145(&
      &62)
      acd145(60)=acd145(12)*acd145(60)
      acd145(62)=-acd145(8)*acd145(9)
      acd145(63)=-acd145(6)*acd145(7)
      acd145(65)=-acd145(4)*acd145(5)
      acd145(66)=-acd145(2)*acd145(3)
      acd145(62)=acd145(66)+acd145(65)+acd145(63)+acd145(10)+acd145(62)
      acd145(62)=acd145(1)*acd145(62)
      acd145(63)=acd145(42)*acd145(55)
      acd145(65)=acd145(40)*acd145(54)
      acd145(66)=acd145(38)*acd145(53)
      acd145(67)=-acd145(9)*acd145(61)
      acd145(63)=acd145(67)+acd145(66)+acd145(63)+acd145(65)
      acd145(63)=acd145(20)*acd145(63)
      acd145(65)=-acd145(28)*acd145(52)
      acd145(66)=acd145(26)*acd145(51)
      acd145(67)=acd145(24)*acd145(50)
      acd145(61)=-acd145(7)*acd145(61)
      acd145(61)=acd145(61)+acd145(67)+acd145(65)+acd145(66)
      acd145(61)=acd145(19)*acd145(61)
      acd145(65)=acd145(48)*acd145(55)
      acd145(66)=acd145(47)*acd145(54)
      acd145(67)=acd145(46)*acd145(53)
      acd145(68)=-acd145(9)*acd145(64)
      acd145(65)=acd145(68)+acd145(67)+acd145(65)+acd145(66)
      acd145(65)=acd145(15)*acd145(65)
      acd145(66)=-acd145(34)*acd145(52)
      acd145(67)=acd145(33)*acd145(51)
      acd145(68)=acd145(32)*acd145(50)
      acd145(64)=-acd145(7)*acd145(64)
      acd145(64)=acd145(64)+acd145(68)+acd145(66)+acd145(67)
      acd145(64)=acd145(14)*acd145(64)
      brack=acd145(56)+acd145(57)+acd145(58)+acd145(59)+acd145(60)+acd145(61)+2&
      &.0_ki*acd145(62)+acd145(63)+acd145(64)+acd145(65)
   end function brack_3
!---#] function brack_3:
!---#[ function brack_4:
   pure function brack_4(Q, mu2) result(brack)
      use p9_part21part21_part25part25part21_model
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_color
      use p9_part21part21_part25part25part21_abbrevd145h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(33) :: acd145
      complex(ki) :: brack
      acd145(1)=d(iv1,iv2)
      acd145(2)=spvae1e2(iv3)
      acd145(3)=abb145(25)
      acd145(4)=spvae2e1(iv3)
      acd145(5)=abb145(19)
      acd145(6)=spvae1e5(iv3)
      acd145(7)=abb145(27)
      acd145(8)=spvae5e1(iv3)
      acd145(9)=abb145(22)
      acd145(10)=d(iv1,iv3)
      acd145(11)=spvae1e2(iv2)
      acd145(12)=spvae2e1(iv2)
      acd145(13)=spvae1e5(iv2)
      acd145(14)=spvae5e1(iv2)
      acd145(15)=d(iv2,iv3)
      acd145(16)=spvae1e2(iv1)
      acd145(17)=spvae2e1(iv1)
      acd145(18)=spvae1e5(iv1)
      acd145(19)=spvae5e1(iv1)
      acd145(20)=spvae2e5(iv3)
      acd145(21)=abb145(7)
      acd145(22)=spvae2e5(iv2)
      acd145(23)=spvae2e5(iv1)
      acd145(24)=spvae5e2(iv3)
      acd145(25)=spvae5e2(iv2)
      acd145(26)=spvae5e2(iv1)
      acd145(27)=-acd145(13)*acd145(17)
      acd145(28)=-acd145(12)*acd145(18)
      acd145(27)=acd145(28)+acd145(27)
      acd145(27)=acd145(24)*acd145(27)
      acd145(28)=-acd145(14)*acd145(16)
      acd145(29)=-acd145(11)*acd145(19)
      acd145(28)=acd145(29)+acd145(28)
      acd145(28)=acd145(20)*acd145(28)
      acd145(29)=-acd145(16)*acd145(22)
      acd145(30)=-acd145(11)*acd145(23)
      acd145(29)=acd145(29)+acd145(30)
      acd145(29)=acd145(8)*acd145(29)
      acd145(30)=-acd145(17)*acd145(25)
      acd145(31)=-acd145(12)*acd145(26)
      acd145(30)=acd145(30)+acd145(31)
      acd145(30)=acd145(6)*acd145(30)
      acd145(31)=-acd145(18)*acd145(25)
      acd145(32)=-acd145(13)*acd145(26)
      acd145(31)=acd145(31)+acd145(32)
      acd145(31)=acd145(4)*acd145(31)
      acd145(32)=-acd145(19)*acd145(22)
      acd145(33)=-acd145(14)*acd145(23)
      acd145(32)=acd145(32)+acd145(33)
      acd145(32)=acd145(2)*acd145(32)
      acd145(27)=acd145(32)+acd145(31)+acd145(30)+acd145(29)+acd145(28)+acd145(&
      &27)
      acd145(27)=acd145(21)*acd145(27)
      acd145(28)=acd145(9)*acd145(19)
      acd145(29)=acd145(7)*acd145(18)
      acd145(30)=acd145(5)*acd145(17)
      acd145(31)=acd145(3)*acd145(16)
      acd145(28)=acd145(31)+acd145(30)+acd145(28)+acd145(29)
      acd145(28)=acd145(15)*acd145(28)
      acd145(29)=acd145(9)*acd145(14)
      acd145(30)=acd145(7)*acd145(13)
      acd145(31)=acd145(5)*acd145(12)
      acd145(32)=acd145(3)*acd145(11)
      acd145(29)=acd145(32)+acd145(31)+acd145(29)+acd145(30)
      acd145(29)=acd145(10)*acd145(29)
      acd145(30)=acd145(8)*acd145(9)
      acd145(31)=acd145(6)*acd145(7)
      acd145(32)=acd145(4)*acd145(5)
      acd145(33)=acd145(2)*acd145(3)
      acd145(30)=acd145(33)+acd145(32)+acd145(30)+acd145(31)
      acd145(30)=acd145(1)*acd145(30)
      acd145(28)=acd145(30)+acd145(28)+acd145(29)
      brack=acd145(27)+2.0_ki*acd145(28)
   end function brack_4
!---#] function brack_4:
!---#[ function brack_5:
   pure function brack_5(Q, mu2) result(brack)
      use p9_part21part21_part25part25part21_model
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_color
      use p9_part21part21_part25part25part21_abbrevd145h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(1) :: acd145
      complex(ki) :: brack
      brack=0.0_ki
   end function brack_5
!---#] function brack_5:
!---#[ function derivative:
   function derivative(mu2,i1,i2,i3,i4) result(numerator)
      use p9_part21part21_part25part25part21_globalsl1, only: epspow
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_abbrevd145h0
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
end module     p9_part21part21_part25part25part21_d145h0l1d
