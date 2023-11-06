module     p2_part21part21_part25part25part21_d52h0l1d_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d52h0l1d_qp.f90
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
      use p2_part21part21_part25part25part21_abbrevd52h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(61) :: acd52
      complex(ki) :: brack
      acd52(1)=dotproduct(qshift,qshift)
      acd52(2)=dotproduct(qshift,spvae1e2)
      acd52(3)=abb52(16)
      acd52(4)=dotproduct(qshift,spvae2e1)
      acd52(5)=abb52(24)
      acd52(6)=dotproduct(qshift,spvae1e5)
      acd52(7)=abb52(30)
      acd52(8)=dotproduct(qshift,spvae5e1)
      acd52(9)=abb52(26)
      acd52(10)=abb52(4)
      acd52(11)=dotproduct(qshift,spvae2e5)
      acd52(12)=abb52(12)
      acd52(13)=abb52(18)
      acd52(14)=dotproduct(qshift,spvae2k1)
      acd52(15)=abb52(13)
      acd52(16)=dotproduct(qshift,spvae5k2)
      acd52(17)=abb52(37)
      acd52(18)=dotproduct(qshift,spvae2k5)
      acd52(19)=abb52(33)
      acd52(20)=abb52(2)
      acd52(21)=abb52(10)
      acd52(22)=dotproduct(qshift,spvae5e2)
      acd52(23)=abb52(31)
      acd52(24)=dotproduct(qshift,spvak1e2)
      acd52(25)=abb52(19)
      acd52(26)=dotproduct(qshift,spvak2e5)
      acd52(27)=abb52(15)
      acd52(28)=dotproduct(qshift,spvak5e2)
      acd52(29)=abb52(35)
      acd52(30)=abb52(5)
      acd52(31)=abb52(9)
      acd52(32)=abb52(28)
      acd52(33)=abb52(17)
      acd52(34)=abb52(27)
      acd52(35)=abb52(14)
      acd52(36)=abb52(23)
      acd52(37)=abb52(7)
      acd52(38)=abb52(21)
      acd52(39)=abb52(22)
      acd52(40)=abb52(3)
      acd52(41)=abb52(6)
      acd52(42)=abb52(38)
      acd52(43)=abb52(36)
      acd52(44)=abb52(34)
      acd52(45)=abb52(32)
      acd52(46)=abb52(11)
      acd52(47)=abb52(20)
      acd52(48)=abb52(8)
      acd52(49)=acd52(28)*acd52(29)
      acd52(50)=acd52(26)*acd52(27)
      acd52(51)=acd52(24)*acd52(25)
      acd52(52)=acd52(22)*acd52(30)
      acd52(53)=-acd52(1)*acd52(5)
      acd52(54)=-acd52(22)*acd52(12)
      acd52(54)=acd52(23)+acd52(54)
      acd52(54)=acd52(6)*acd52(54)
      acd52(49)=acd52(54)+acd52(53)+acd52(52)+acd52(51)+acd52(50)-acd52(31)+acd&
      &52(49)
      acd52(49)=acd52(4)*acd52(49)
      acd52(50)=acd52(18)*acd52(19)
      acd52(51)=acd52(16)*acd52(17)
      acd52(52)=acd52(14)*acd52(15)
      acd52(53)=acd52(11)*acd52(20)
      acd52(54)=-acd52(1)*acd52(3)
      acd52(55)=-acd52(11)*acd52(12)
      acd52(55)=acd52(13)+acd52(55)
      acd52(55)=acd52(8)*acd52(55)
      acd52(50)=acd52(55)+acd52(54)+acd52(53)+acd52(52)+acd52(51)-acd52(21)+acd&
      &52(50)
      acd52(50)=acd52(2)*acd52(50)
      acd52(51)=acd52(28)*acd52(38)
      acd52(52)=acd52(26)*acd52(37)
      acd52(53)=acd52(24)*acd52(36)
      acd52(54)=-acd52(1)*acd52(9)
      acd52(51)=acd52(54)+acd52(53)+acd52(52)-acd52(39)+acd52(51)
      acd52(51)=acd52(8)*acd52(51)
      acd52(52)=acd52(18)*acd52(34)
      acd52(53)=acd52(16)*acd52(33)
      acd52(54)=acd52(14)*acd52(32)
      acd52(55)=-acd52(1)*acd52(7)
      acd52(52)=acd52(55)+acd52(54)+acd52(53)-acd52(35)+acd52(52)
      acd52(52)=acd52(6)*acd52(52)
      acd52(53)=-acd52(28)*acd52(44)
      acd52(54)=-acd52(26)*acd52(42)
      acd52(55)=-acd52(24)*acd52(40)
      acd52(56)=-acd52(22)*acd52(47)
      acd52(57)=-acd52(18)*acd52(45)
      acd52(58)=-acd52(16)*acd52(43)
      acd52(59)=-acd52(14)*acd52(41)
      acd52(60)=-acd52(11)*acd52(46)
      acd52(61)=acd52(1)*acd52(10)
      brack=acd52(48)+acd52(49)+acd52(50)+acd52(51)+acd52(52)+acd52(53)+acd52(5&
      &4)+acd52(55)+acd52(56)+acd52(57)+acd52(58)+acd52(59)+acd52(60)+acd52(61)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p2_part21part21_part25part25part21_model_qp
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_color_qp
      use p2_part21part21_part25part25part21_abbrevd52h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(77) :: acd52
      complex(ki) :: brack
      acd52(1)=qshift(iv1)
      acd52(2)=dotproduct(qshift,spvae1e2)
      acd52(3)=abb52(16)
      acd52(4)=dotproduct(qshift,spvae2e1)
      acd52(5)=abb52(24)
      acd52(6)=dotproduct(qshift,spvae1e5)
      acd52(7)=abb52(30)
      acd52(8)=dotproduct(qshift,spvae5e1)
      acd52(9)=abb52(26)
      acd52(10)=abb52(4)
      acd52(11)=spvae1e2(iv1)
      acd52(12)=dotproduct(qshift,qshift)
      acd52(13)=dotproduct(qshift,spvae2e5)
      acd52(14)=abb52(12)
      acd52(15)=abb52(18)
      acd52(16)=dotproduct(qshift,spvae2k1)
      acd52(17)=abb52(13)
      acd52(18)=dotproduct(qshift,spvae5k2)
      acd52(19)=abb52(37)
      acd52(20)=dotproduct(qshift,spvae2k5)
      acd52(21)=abb52(33)
      acd52(22)=abb52(2)
      acd52(23)=abb52(10)
      acd52(24)=spvae2e1(iv1)
      acd52(25)=dotproduct(qshift,spvae5e2)
      acd52(26)=abb52(31)
      acd52(27)=dotproduct(qshift,spvak1e2)
      acd52(28)=abb52(19)
      acd52(29)=dotproduct(qshift,spvak2e5)
      acd52(30)=abb52(15)
      acd52(31)=dotproduct(qshift,spvak5e2)
      acd52(32)=abb52(35)
      acd52(33)=abb52(5)
      acd52(34)=abb52(9)
      acd52(35)=spvae1e5(iv1)
      acd52(36)=abb52(28)
      acd52(37)=abb52(17)
      acd52(38)=abb52(27)
      acd52(39)=abb52(14)
      acd52(40)=spvae5e1(iv1)
      acd52(41)=abb52(23)
      acd52(42)=abb52(7)
      acd52(43)=abb52(21)
      acd52(44)=abb52(22)
      acd52(45)=spvak1e2(iv1)
      acd52(46)=abb52(3)
      acd52(47)=spvae2k1(iv1)
      acd52(48)=abb52(6)
      acd52(49)=spvak2e5(iv1)
      acd52(50)=abb52(38)
      acd52(51)=spvae5k2(iv1)
      acd52(52)=abb52(36)
      acd52(53)=spvak5e2(iv1)
      acd52(54)=abb52(34)
      acd52(55)=spvae2k5(iv1)
      acd52(56)=abb52(32)
      acd52(57)=spvae2e5(iv1)
      acd52(58)=abb52(11)
      acd52(59)=spvae5e2(iv1)
      acd52(60)=abb52(20)
      acd52(61)=-acd52(31)*acd52(32)
      acd52(62)=-acd52(29)*acd52(30)
      acd52(63)=-acd52(27)*acd52(28)
      acd52(64)=-acd52(25)*acd52(33)
      acd52(65)=acd52(12)*acd52(5)
      acd52(66)=acd52(14)*acd52(25)
      acd52(66)=acd52(66)-acd52(26)
      acd52(67)=acd52(6)*acd52(66)
      acd52(61)=acd52(67)+acd52(65)+acd52(64)+acd52(63)+acd52(62)+acd52(34)+acd&
      &52(61)
      acd52(61)=acd52(24)*acd52(61)
      acd52(62)=-acd52(20)*acd52(21)
      acd52(63)=-acd52(18)*acd52(19)
      acd52(64)=-acd52(16)*acd52(17)
      acd52(65)=-acd52(13)*acd52(22)
      acd52(67)=acd52(12)*acd52(3)
      acd52(68)=acd52(14)*acd52(13)
      acd52(68)=acd52(68)-acd52(15)
      acd52(69)=acd52(8)*acd52(68)
      acd52(62)=acd52(69)+acd52(67)+acd52(65)+acd52(64)+acd52(63)+acd52(23)+acd&
      &52(62)
      acd52(62)=acd52(11)*acd52(62)
      acd52(63)=acd52(6)*acd52(14)
      acd52(63)=acd52(63)-acd52(33)
      acd52(63)=acd52(59)*acd52(63)
      acd52(64)=-acd52(53)*acd52(32)
      acd52(65)=-acd52(49)*acd52(30)
      acd52(67)=-acd52(45)*acd52(28)
      acd52(69)=2.0_ki*acd52(1)
      acd52(70)=acd52(5)*acd52(69)
      acd52(66)=acd52(35)*acd52(66)
      acd52(63)=acd52(66)+acd52(70)+acd52(67)+acd52(65)+acd52(64)+acd52(63)
      acd52(63)=acd52(4)*acd52(63)
      acd52(64)=acd52(8)*acd52(14)
      acd52(64)=acd52(64)-acd52(22)
      acd52(64)=acd52(57)*acd52(64)
      acd52(65)=-acd52(55)*acd52(21)
      acd52(66)=-acd52(51)*acd52(19)
      acd52(67)=-acd52(47)*acd52(17)
      acd52(70)=acd52(3)*acd52(69)
      acd52(68)=acd52(40)*acd52(68)
      acd52(64)=acd52(68)+acd52(70)+acd52(67)+acd52(66)+acd52(65)+acd52(64)
      acd52(64)=acd52(2)*acd52(64)
      acd52(65)=-acd52(31)*acd52(43)
      acd52(66)=-acd52(29)*acd52(42)
      acd52(67)=-acd52(27)*acd52(41)
      acd52(68)=acd52(12)*acd52(9)
      acd52(65)=acd52(68)+acd52(67)+acd52(66)+acd52(44)+acd52(65)
      acd52(65)=acd52(40)*acd52(65)
      acd52(66)=-acd52(20)*acd52(38)
      acd52(67)=-acd52(18)*acd52(37)
      acd52(68)=-acd52(16)*acd52(36)
      acd52(70)=acd52(12)*acd52(7)
      acd52(66)=acd52(70)+acd52(68)+acd52(67)+acd52(39)+acd52(66)
      acd52(66)=acd52(35)*acd52(66)
      acd52(67)=-acd52(53)*acd52(43)
      acd52(68)=-acd52(49)*acd52(42)
      acd52(70)=-acd52(45)*acd52(41)
      acd52(71)=acd52(9)*acd52(69)
      acd52(67)=acd52(71)+acd52(70)+acd52(67)+acd52(68)
      acd52(67)=acd52(8)*acd52(67)
      acd52(68)=-acd52(55)*acd52(38)
      acd52(70)=-acd52(51)*acd52(37)
      acd52(71)=-acd52(47)*acd52(36)
      acd52(72)=acd52(7)*acd52(69)
      acd52(68)=acd52(72)+acd52(71)+acd52(68)+acd52(70)
      acd52(68)=acd52(6)*acd52(68)
      acd52(70)=acd52(59)*acd52(60)
      acd52(71)=acd52(57)*acd52(58)
      acd52(72)=acd52(55)*acd52(56)
      acd52(73)=acd52(53)*acd52(54)
      acd52(74)=acd52(51)*acd52(52)
      acd52(75)=acd52(49)*acd52(50)
      acd52(76)=acd52(47)*acd52(48)
      acd52(77)=acd52(45)*acd52(46)
      acd52(69)=-acd52(10)*acd52(69)
      brack=acd52(61)+acd52(62)+acd52(63)+acd52(64)+acd52(65)+acd52(66)+acd52(6&
      &7)+acd52(68)+acd52(69)+acd52(70)+acd52(71)+acd52(72)+acd52(73)+acd52(74)&
      &+acd52(75)+acd52(76)+acd52(77)
   end function brack_2
!---#] function brack_2:
!---#[ function brack_3:
   pure function brack_3(Q, mu2) result(brack)
      use p2_part21part21_part25part25part21_model_qp
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_color_qp
      use p2_part21part21_part25part25part21_abbrevd52h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(68) :: acd52
      complex(ki) :: brack
      acd52(1)=d(iv1,iv2)
      acd52(2)=dotproduct(qshift,spvae1e2)
      acd52(3)=abb52(16)
      acd52(4)=dotproduct(qshift,spvae2e1)
      acd52(5)=abb52(24)
      acd52(6)=dotproduct(qshift,spvae1e5)
      acd52(7)=abb52(30)
      acd52(8)=dotproduct(qshift,spvae5e1)
      acd52(9)=abb52(26)
      acd52(10)=abb52(4)
      acd52(11)=qshift(iv1)
      acd52(12)=spvae1e2(iv2)
      acd52(13)=spvae2e1(iv2)
      acd52(14)=spvae1e5(iv2)
      acd52(15)=spvae5e1(iv2)
      acd52(16)=qshift(iv2)
      acd52(17)=spvae1e2(iv1)
      acd52(18)=spvae2e1(iv1)
      acd52(19)=spvae1e5(iv1)
      acd52(20)=spvae5e1(iv1)
      acd52(21)=dotproduct(qshift,spvae2e5)
      acd52(22)=abb52(12)
      acd52(23)=abb52(18)
      acd52(24)=spvae2k1(iv2)
      acd52(25)=abb52(13)
      acd52(26)=spvae5k2(iv2)
      acd52(27)=abb52(37)
      acd52(28)=spvae2k5(iv2)
      acd52(29)=abb52(33)
      acd52(30)=spvae2e5(iv2)
      acd52(31)=abb52(2)
      acd52(32)=spvae2k1(iv1)
      acd52(33)=spvae5k2(iv1)
      acd52(34)=spvae2k5(iv1)
      acd52(35)=spvae2e5(iv1)
      acd52(36)=dotproduct(qshift,spvae5e2)
      acd52(37)=abb52(31)
      acd52(38)=spvak1e2(iv2)
      acd52(39)=abb52(19)
      acd52(40)=spvak2e5(iv2)
      acd52(41)=abb52(15)
      acd52(42)=spvak5e2(iv2)
      acd52(43)=abb52(35)
      acd52(44)=spvae5e2(iv2)
      acd52(45)=abb52(5)
      acd52(46)=spvak1e2(iv1)
      acd52(47)=spvak2e5(iv1)
      acd52(48)=spvak5e2(iv1)
      acd52(49)=spvae5e2(iv1)
      acd52(50)=abb52(28)
      acd52(51)=abb52(17)
      acd52(52)=abb52(27)
      acd52(53)=abb52(23)
      acd52(54)=abb52(7)
      acd52(55)=abb52(21)
      acd52(56)=-acd52(19)*acd52(44)
      acd52(57)=-acd52(14)*acd52(49)
      acd52(56)=acd52(57)+acd52(56)
      acd52(56)=acd52(4)*acd52(56)
      acd52(57)=-acd52(20)*acd52(30)
      acd52(58)=-acd52(15)*acd52(35)
      acd52(57)=acd52(58)+acd52(57)
      acd52(57)=acd52(2)*acd52(57)
      acd52(58)=-acd52(6)*acd52(44)
      acd52(59)=-acd52(14)*acd52(36)
      acd52(58)=acd52(58)+acd52(59)
      acd52(58)=acd52(18)*acd52(58)
      acd52(59)=-acd52(8)*acd52(30)
      acd52(60)=-acd52(15)*acd52(21)
      acd52(59)=acd52(59)+acd52(60)
      acd52(59)=acd52(17)*acd52(59)
      acd52(60)=-acd52(6)*acd52(49)
      acd52(61)=-acd52(19)*acd52(36)
      acd52(60)=acd52(60)+acd52(61)
      acd52(60)=acd52(13)*acd52(60)
      acd52(61)=-acd52(8)*acd52(35)
      acd52(62)=-acd52(20)*acd52(21)
      acd52(61)=acd52(61)+acd52(62)
      acd52(61)=acd52(12)*acd52(61)
      acd52(56)=acd52(61)+acd52(60)+acd52(59)+acd52(58)+acd52(57)+acd52(56)
      acd52(56)=acd52(22)*acd52(56)
      acd52(57)=acd52(42)*acd52(43)
      acd52(58)=acd52(40)*acd52(41)
      acd52(59)=acd52(38)*acd52(39)
      acd52(60)=acd52(44)*acd52(45)
      acd52(61)=2.0_ki*acd52(16)
      acd52(62)=-acd52(5)*acd52(61)
      acd52(63)=acd52(14)*acd52(37)
      acd52(57)=acd52(63)+acd52(62)+acd52(60)+acd52(59)+acd52(57)+acd52(58)
      acd52(57)=acd52(18)*acd52(57)
      acd52(58)=acd52(28)*acd52(29)
      acd52(59)=acd52(26)*acd52(27)
      acd52(60)=acd52(24)*acd52(25)
      acd52(62)=acd52(30)*acd52(31)
      acd52(63)=-acd52(3)*acd52(61)
      acd52(64)=acd52(15)*acd52(23)
      acd52(58)=acd52(64)+acd52(63)+acd52(62)+acd52(60)+acd52(58)+acd52(59)
      acd52(58)=acd52(17)*acd52(58)
      acd52(59)=acd52(43)*acd52(48)
      acd52(60)=acd52(41)*acd52(47)
      acd52(62)=acd52(39)*acd52(46)
      acd52(63)=acd52(49)*acd52(45)
      acd52(64)=2.0_ki*acd52(11)
      acd52(65)=-acd52(5)*acd52(64)
      acd52(66)=acd52(19)*acd52(37)
      acd52(59)=acd52(66)+acd52(65)+acd52(63)+acd52(62)+acd52(59)+acd52(60)
      acd52(59)=acd52(13)*acd52(59)
      acd52(60)=acd52(29)*acd52(34)
      acd52(62)=acd52(27)*acd52(33)
      acd52(63)=acd52(25)*acd52(32)
      acd52(65)=acd52(35)*acd52(31)
      acd52(66)=-acd52(3)*acd52(64)
      acd52(67)=acd52(20)*acd52(23)
      acd52(60)=acd52(67)+acd52(66)+acd52(65)+acd52(63)+acd52(60)+acd52(62)
      acd52(60)=acd52(12)*acd52(60)
      acd52(62)=-acd52(8)*acd52(9)
      acd52(63)=-acd52(6)*acd52(7)
      acd52(65)=-acd52(4)*acd52(5)
      acd52(66)=-acd52(2)*acd52(3)
      acd52(62)=acd52(66)+acd52(65)+acd52(63)+acd52(10)+acd52(62)
      acd52(62)=acd52(1)*acd52(62)
      acd52(63)=acd52(42)*acd52(55)
      acd52(65)=acd52(40)*acd52(54)
      acd52(66)=acd52(38)*acd52(53)
      acd52(67)=-acd52(9)*acd52(61)
      acd52(63)=acd52(67)+acd52(66)+acd52(63)+acd52(65)
      acd52(63)=acd52(20)*acd52(63)
      acd52(65)=acd52(28)*acd52(52)
      acd52(66)=acd52(26)*acd52(51)
      acd52(67)=acd52(24)*acd52(50)
      acd52(61)=-acd52(7)*acd52(61)
      acd52(61)=acd52(61)+acd52(67)+acd52(65)+acd52(66)
      acd52(61)=acd52(19)*acd52(61)
      acd52(65)=acd52(48)*acd52(55)
      acd52(66)=acd52(47)*acd52(54)
      acd52(67)=acd52(46)*acd52(53)
      acd52(68)=-acd52(9)*acd52(64)
      acd52(65)=acd52(68)+acd52(67)+acd52(65)+acd52(66)
      acd52(65)=acd52(15)*acd52(65)
      acd52(66)=acd52(34)*acd52(52)
      acd52(67)=acd52(33)*acd52(51)
      acd52(68)=acd52(32)*acd52(50)
      acd52(64)=-acd52(7)*acd52(64)
      acd52(64)=acd52(64)+acd52(68)+acd52(66)+acd52(67)
      acd52(64)=acd52(14)*acd52(64)
      brack=acd52(56)+acd52(57)+acd52(58)+acd52(59)+acd52(60)+acd52(61)+2.0_ki*&
      &acd52(62)+acd52(63)+acd52(64)+acd52(65)
   end function brack_3
!---#] function brack_3:
!---#[ function brack_4:
   pure function brack_4(Q, mu2) result(brack)
      use p2_part21part21_part25part25part21_model_qp
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_color_qp
      use p2_part21part21_part25part25part21_abbrevd52h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(33) :: acd52
      complex(ki) :: brack
      acd52(1)=d(iv1,iv2)
      acd52(2)=spvae1e2(iv3)
      acd52(3)=abb52(16)
      acd52(4)=spvae2e1(iv3)
      acd52(5)=abb52(24)
      acd52(6)=spvae1e5(iv3)
      acd52(7)=abb52(30)
      acd52(8)=spvae5e1(iv3)
      acd52(9)=abb52(26)
      acd52(10)=d(iv1,iv3)
      acd52(11)=spvae1e2(iv2)
      acd52(12)=spvae2e1(iv2)
      acd52(13)=spvae1e5(iv2)
      acd52(14)=spvae5e1(iv2)
      acd52(15)=d(iv2,iv3)
      acd52(16)=spvae1e2(iv1)
      acd52(17)=spvae2e1(iv1)
      acd52(18)=spvae1e5(iv1)
      acd52(19)=spvae5e1(iv1)
      acd52(20)=spvae2e5(iv3)
      acd52(21)=abb52(12)
      acd52(22)=spvae2e5(iv2)
      acd52(23)=spvae2e5(iv1)
      acd52(24)=spvae5e2(iv3)
      acd52(25)=spvae5e2(iv2)
      acd52(26)=spvae5e2(iv1)
      acd52(27)=acd52(13)*acd52(17)
      acd52(28)=acd52(12)*acd52(18)
      acd52(27)=acd52(28)+acd52(27)
      acd52(27)=acd52(24)*acd52(27)
      acd52(28)=acd52(14)*acd52(16)
      acd52(29)=acd52(11)*acd52(19)
      acd52(28)=acd52(29)+acd52(28)
      acd52(28)=acd52(20)*acd52(28)
      acd52(29)=acd52(16)*acd52(22)
      acd52(30)=acd52(11)*acd52(23)
      acd52(29)=acd52(29)+acd52(30)
      acd52(29)=acd52(8)*acd52(29)
      acd52(30)=acd52(17)*acd52(25)
      acd52(31)=acd52(12)*acd52(26)
      acd52(30)=acd52(30)+acd52(31)
      acd52(30)=acd52(6)*acd52(30)
      acd52(31)=acd52(18)*acd52(25)
      acd52(32)=acd52(13)*acd52(26)
      acd52(31)=acd52(31)+acd52(32)
      acd52(31)=acd52(4)*acd52(31)
      acd52(32)=acd52(19)*acd52(22)
      acd52(33)=acd52(14)*acd52(23)
      acd52(32)=acd52(32)+acd52(33)
      acd52(32)=acd52(2)*acd52(32)
      acd52(27)=acd52(32)+acd52(31)+acd52(30)+acd52(29)+acd52(28)+acd52(27)
      acd52(27)=acd52(21)*acd52(27)
      acd52(28)=acd52(9)*acd52(19)
      acd52(29)=acd52(7)*acd52(18)
      acd52(30)=acd52(5)*acd52(17)
      acd52(31)=acd52(3)*acd52(16)
      acd52(28)=acd52(31)+acd52(30)+acd52(28)+acd52(29)
      acd52(28)=acd52(15)*acd52(28)
      acd52(29)=acd52(9)*acd52(14)
      acd52(30)=acd52(7)*acd52(13)
      acd52(31)=acd52(5)*acd52(12)
      acd52(32)=acd52(3)*acd52(11)
      acd52(29)=acd52(32)+acd52(31)+acd52(29)+acd52(30)
      acd52(29)=acd52(10)*acd52(29)
      acd52(30)=acd52(8)*acd52(9)
      acd52(31)=acd52(6)*acd52(7)
      acd52(32)=acd52(4)*acd52(5)
      acd52(33)=acd52(2)*acd52(3)
      acd52(30)=acd52(33)+acd52(32)+acd52(30)+acd52(31)
      acd52(30)=acd52(1)*acd52(30)
      acd52(28)=acd52(30)+acd52(28)+acd52(29)
      brack=acd52(27)+2.0_ki*acd52(28)
   end function brack_4
!---#] function brack_4:
!---#[ function brack_5:
   pure function brack_5(Q, mu2) result(brack)
      use p2_part21part21_part25part25part21_model_qp
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_color_qp
      use p2_part21part21_part25part25part21_abbrevd52h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(1) :: acd52
      complex(ki) :: brack
      brack=0.0_ki
   end function brack_5
!---#] function brack_5:
!---#[ function derivative:
   function derivative(mu2,i1,i2,i3,i4) result(numerator)
      use p2_part21part21_part25part25part21_globalsl1_qp, only: epspow
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_abbrevd52h0_qp
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
end module     p2_part21part21_part25part25part21_d52h0l1d_qp
