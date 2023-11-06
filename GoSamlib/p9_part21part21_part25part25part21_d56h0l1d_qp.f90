module     p9_part21part21_part25part25part21_d56h0l1d_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p9_part21part &
   ! &21_part25part25part21/helicity0d56h0l1d_qp.f90
   ! generator: buildfortran_d.py
   use p9_part21part21_part25part25part21_config, only: ki => ki_qp
   use p9_part21part21_part25part25part21_util_qp, only: cond, d => metric_tens&
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
      use p9_part21part21_part25part25part21_model_qp
      use p9_part21part21_part25part25part21_kinematics_qp
      use p9_part21part21_part25part25part21_color_qp
      use p9_part21part21_part25part25part21_abbrevd56h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(57) :: acd56
      complex(ki) :: brack
      acd56(1)=dotproduct(qshift,qshift)
      acd56(2)=dotproduct(qshift,spvae1e5)
      acd56(3)=abb56(25)
      acd56(4)=dotproduct(qshift,spvae5e1)
      acd56(5)=abb56(30)
      acd56(6)=dotproduct(qshift,spvae2e5)
      acd56(7)=abb56(27)
      acd56(8)=dotproduct(qshift,spvae5e2)
      acd56(9)=abb56(9)
      acd56(10)=abb56(13)
      acd56(11)=dotproduct(qshift,spvae2e1)
      acd56(12)=abb56(7)
      acd56(13)=abb56(6)
      acd56(14)=dotproduct(qshift,spvak1e2)
      acd56(15)=abb56(5)
      acd56(16)=dotproduct(qshift,spvak2e1)
      acd56(17)=abb56(35)
      acd56(18)=dotproduct(qshift,spvak5e1)
      acd56(19)=abb56(32)
      acd56(20)=abb56(18)
      acd56(21)=abb56(2)
      acd56(22)=dotproduct(qshift,spvae1e2)
      acd56(23)=abb56(4)
      acd56(24)=dotproduct(qshift,spvae2k1)
      acd56(25)=abb56(29)
      acd56(26)=dotproduct(qshift,spvae1k2)
      acd56(27)=abb56(34)
      acd56(28)=dotproduct(qshift,spvae1k5)
      acd56(29)=abb56(31)
      acd56(30)=abb56(15)
      acd56(31)=abb56(21)
      acd56(32)=abb56(16)
      acd56(33)=abb56(19)
      acd56(34)=abb56(26)
      acd56(35)=abb56(20)
      acd56(36)=abb56(14)
      acd56(37)=abb56(23)
      acd56(38)=abb56(28)
      acd56(39)=abb56(17)
      acd56(40)=abb56(12)
      acd56(41)=abb56(10)
      acd56(42)=abb56(24)
      acd56(43)=abb56(33)
      acd56(44)=abb56(11)
      acd56(45)=abb56(8)
      acd56(46)=abb56(3)
      acd56(47)=acd56(28)*acd56(29)
      acd56(48)=acd56(26)*acd56(27)
      acd56(49)=acd56(24)*acd56(25)
      acd56(50)=acd56(22)*acd56(30)
      acd56(51)=-acd56(1)*acd56(5)
      acd56(52)=-acd56(22)*acd56(12)
      acd56(52)=acd56(23)+acd56(52)
      acd56(52)=acd56(6)*acd56(52)
      acd56(47)=acd56(52)+acd56(51)+acd56(50)+acd56(49)+acd56(48)-acd56(31)+acd&
      &56(47)
      acd56(47)=acd56(4)*acd56(47)
      acd56(48)=acd56(18)*acd56(19)
      acd56(49)=acd56(16)*acd56(17)
      acd56(50)=acd56(14)*acd56(15)
      acd56(51)=acd56(11)*acd56(20)
      acd56(52)=-acd56(1)*acd56(3)
      acd56(53)=-acd56(11)*acd56(12)
      acd56(53)=acd56(13)+acd56(53)
      acd56(53)=acd56(8)*acd56(53)
      acd56(48)=acd56(53)+acd56(52)+acd56(51)+acd56(50)+acd56(49)-acd56(21)+acd&
      &56(48)
      acd56(48)=acd56(2)*acd56(48)
      acd56(49)=acd56(28)*acd56(38)
      acd56(50)=acd56(26)*acd56(37)
      acd56(51)=acd56(24)*acd56(36)
      acd56(52)=-acd56(1)*acd56(9)
      acd56(49)=acd56(52)+acd56(51)+acd56(50)-acd56(39)+acd56(49)
      acd56(49)=acd56(8)*acd56(49)
      acd56(50)=acd56(18)*acd56(34)
      acd56(51)=acd56(16)*acd56(33)
      acd56(52)=acd56(14)*acd56(32)
      acd56(53)=-acd56(1)*acd56(7)
      acd56(50)=acd56(53)+acd56(52)+acd56(51)-acd56(35)+acd56(50)
      acd56(50)=acd56(6)*acd56(50)
      acd56(51)=-acd56(26)*acd56(43)
      acd56(52)=-acd56(24)*acd56(41)
      acd56(53)=-acd56(22)*acd56(44)
      acd56(54)=-acd56(16)*acd56(42)
      acd56(55)=-acd56(14)*acd56(40)
      acd56(56)=-acd56(11)*acd56(45)
      acd56(57)=acd56(1)*acd56(10)
      brack=acd56(46)+acd56(47)+acd56(48)+acd56(49)+acd56(50)+acd56(51)+acd56(5&
      &2)+acd56(53)+acd56(54)+acd56(55)+acd56(56)+acd56(57)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p9_part21part21_part25part25part21_model_qp
      use p9_part21part21_part25part25part21_kinematics_qp
      use p9_part21part21_part25part25part21_color_qp
      use p9_part21part21_part25part25part21_abbrevd56h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(73) :: acd56
      complex(ki) :: brack
      acd56(1)=qshift(iv1)
      acd56(2)=dotproduct(qshift,spvae1e5)
      acd56(3)=abb56(25)
      acd56(4)=dotproduct(qshift,spvae5e1)
      acd56(5)=abb56(30)
      acd56(6)=dotproduct(qshift,spvae2e5)
      acd56(7)=abb56(27)
      acd56(8)=dotproduct(qshift,spvae5e2)
      acd56(9)=abb56(9)
      acd56(10)=abb56(13)
      acd56(11)=spvae1e5(iv1)
      acd56(12)=dotproduct(qshift,qshift)
      acd56(13)=dotproduct(qshift,spvae2e1)
      acd56(14)=abb56(7)
      acd56(15)=abb56(6)
      acd56(16)=dotproduct(qshift,spvak1e2)
      acd56(17)=abb56(5)
      acd56(18)=dotproduct(qshift,spvak2e1)
      acd56(19)=abb56(35)
      acd56(20)=dotproduct(qshift,spvak5e1)
      acd56(21)=abb56(32)
      acd56(22)=abb56(18)
      acd56(23)=abb56(2)
      acd56(24)=spvae5e1(iv1)
      acd56(25)=dotproduct(qshift,spvae1e2)
      acd56(26)=abb56(4)
      acd56(27)=dotproduct(qshift,spvae2k1)
      acd56(28)=abb56(29)
      acd56(29)=dotproduct(qshift,spvae1k2)
      acd56(30)=abb56(34)
      acd56(31)=dotproduct(qshift,spvae1k5)
      acd56(32)=abb56(31)
      acd56(33)=abb56(15)
      acd56(34)=abb56(21)
      acd56(35)=spvae2e5(iv1)
      acd56(36)=abb56(16)
      acd56(37)=abb56(19)
      acd56(38)=abb56(26)
      acd56(39)=abb56(20)
      acd56(40)=spvae5e2(iv1)
      acd56(41)=abb56(14)
      acd56(42)=abb56(23)
      acd56(43)=abb56(28)
      acd56(44)=abb56(17)
      acd56(45)=spvak1e2(iv1)
      acd56(46)=abb56(12)
      acd56(47)=spvae2k1(iv1)
      acd56(48)=abb56(10)
      acd56(49)=spvak2e1(iv1)
      acd56(50)=abb56(24)
      acd56(51)=spvae1k2(iv1)
      acd56(52)=abb56(33)
      acd56(53)=spvak5e1(iv1)
      acd56(54)=spvae1k5(iv1)
      acd56(55)=spvae1e2(iv1)
      acd56(56)=abb56(11)
      acd56(57)=spvae2e1(iv1)
      acd56(58)=abb56(8)
      acd56(59)=acd56(31)*acd56(32)
      acd56(60)=acd56(29)*acd56(30)
      acd56(61)=acd56(27)*acd56(28)
      acd56(62)=acd56(25)*acd56(33)
      acd56(63)=-acd56(12)*acd56(5)
      acd56(64)=acd56(14)*acd56(25)
      acd56(64)=acd56(64)-acd56(26)
      acd56(65)=-acd56(6)*acd56(64)
      acd56(59)=acd56(65)+acd56(63)+acd56(62)+acd56(61)+acd56(60)-acd56(34)+acd&
      &56(59)
      acd56(59)=acd56(24)*acd56(59)
      acd56(60)=acd56(20)*acd56(21)
      acd56(61)=acd56(18)*acd56(19)
      acd56(62)=acd56(16)*acd56(17)
      acd56(63)=acd56(13)*acd56(22)
      acd56(65)=-acd56(12)*acd56(3)
      acd56(66)=acd56(14)*acd56(13)
      acd56(66)=acd56(66)-acd56(15)
      acd56(67)=-acd56(8)*acd56(66)
      acd56(60)=acd56(67)+acd56(65)+acd56(63)+acd56(62)+acd56(61)-acd56(23)+acd&
      &56(60)
      acd56(60)=acd56(11)*acd56(60)
      acd56(61)=-acd56(6)*acd56(14)
      acd56(61)=acd56(61)+acd56(33)
      acd56(61)=acd56(55)*acd56(61)
      acd56(62)=acd56(32)*acd56(54)
      acd56(63)=acd56(51)*acd56(30)
      acd56(65)=acd56(47)*acd56(28)
      acd56(67)=2.0_ki*acd56(1)
      acd56(68)=-acd56(5)*acd56(67)
      acd56(64)=-acd56(35)*acd56(64)
      acd56(61)=acd56(64)+acd56(68)+acd56(65)+acd56(63)+acd56(62)+acd56(61)
      acd56(61)=acd56(4)*acd56(61)
      acd56(62)=-acd56(8)*acd56(14)
      acd56(62)=acd56(62)+acd56(22)
      acd56(62)=acd56(57)*acd56(62)
      acd56(63)=acd56(21)*acd56(53)
      acd56(64)=acd56(49)*acd56(19)
      acd56(65)=acd56(45)*acd56(17)
      acd56(68)=-acd56(3)*acd56(67)
      acd56(66)=-acd56(40)*acd56(66)
      acd56(62)=acd56(66)+acd56(68)+acd56(65)+acd56(64)+acd56(63)+acd56(62)
      acd56(62)=acd56(2)*acd56(62)
      acd56(63)=acd56(31)*acd56(43)
      acd56(64)=acd56(29)*acd56(42)
      acd56(65)=acd56(27)*acd56(41)
      acd56(66)=-acd56(12)*acd56(9)
      acd56(63)=acd56(66)+acd56(65)+acd56(64)-acd56(44)+acd56(63)
      acd56(63)=acd56(40)*acd56(63)
      acd56(64)=acd56(20)*acd56(38)
      acd56(65)=acd56(18)*acd56(37)
      acd56(66)=acd56(16)*acd56(36)
      acd56(68)=-acd56(12)*acd56(7)
      acd56(64)=acd56(68)+acd56(66)+acd56(65)-acd56(39)+acd56(64)
      acd56(64)=acd56(35)*acd56(64)
      acd56(65)=acd56(43)*acd56(54)
      acd56(66)=acd56(51)*acd56(42)
      acd56(68)=acd56(47)*acd56(41)
      acd56(69)=-acd56(9)*acd56(67)
      acd56(65)=acd56(69)+acd56(68)+acd56(65)+acd56(66)
      acd56(65)=acd56(8)*acd56(65)
      acd56(66)=acd56(38)*acd56(53)
      acd56(68)=acd56(49)*acd56(37)
      acd56(69)=acd56(45)*acd56(36)
      acd56(70)=-acd56(7)*acd56(67)
      acd56(66)=acd56(70)+acd56(69)+acd56(66)+acd56(68)
      acd56(66)=acd56(6)*acd56(66)
      acd56(68)=-acd56(57)*acd56(58)
      acd56(69)=-acd56(55)*acd56(56)
      acd56(70)=-acd56(51)*acd56(52)
      acd56(71)=-acd56(49)*acd56(50)
      acd56(72)=-acd56(47)*acd56(48)
      acd56(73)=-acd56(45)*acd56(46)
      acd56(67)=acd56(10)*acd56(67)
      brack=acd56(59)+acd56(60)+acd56(61)+acd56(62)+acd56(63)+acd56(64)+acd56(6&
      &5)+acd56(66)+acd56(67)+acd56(68)+acd56(69)+acd56(70)+acd56(71)+acd56(72)&
      &+acd56(73)
   end function brack_2
!---#] function brack_2:
!---#[ function brack_3:
   pure function brack_3(Q, mu2) result(brack)
      use p9_part21part21_part25part25part21_model_qp
      use p9_part21part21_part25part25part21_kinematics_qp
      use p9_part21part21_part25part25part21_color_qp
      use p9_part21part21_part25part25part21_abbrevd56h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(68) :: acd56
      complex(ki) :: brack
      acd56(1)=d(iv1,iv2)
      acd56(2)=dotproduct(qshift,spvae1e5)
      acd56(3)=abb56(25)
      acd56(4)=dotproduct(qshift,spvae5e1)
      acd56(5)=abb56(30)
      acd56(6)=dotproduct(qshift,spvae2e5)
      acd56(7)=abb56(27)
      acd56(8)=dotproduct(qshift,spvae5e2)
      acd56(9)=abb56(9)
      acd56(10)=abb56(13)
      acd56(11)=qshift(iv1)
      acd56(12)=spvae1e5(iv2)
      acd56(13)=spvae5e1(iv2)
      acd56(14)=spvae2e5(iv2)
      acd56(15)=spvae5e2(iv2)
      acd56(16)=qshift(iv2)
      acd56(17)=spvae1e5(iv1)
      acd56(18)=spvae5e1(iv1)
      acd56(19)=spvae2e5(iv1)
      acd56(20)=spvae5e2(iv1)
      acd56(21)=dotproduct(qshift,spvae2e1)
      acd56(22)=abb56(7)
      acd56(23)=abb56(6)
      acd56(24)=spvak1e2(iv2)
      acd56(25)=abb56(5)
      acd56(26)=spvak2e1(iv2)
      acd56(27)=abb56(35)
      acd56(28)=spvak5e1(iv2)
      acd56(29)=abb56(32)
      acd56(30)=spvae2e1(iv2)
      acd56(31)=abb56(18)
      acd56(32)=spvak1e2(iv1)
      acd56(33)=spvak2e1(iv1)
      acd56(34)=spvak5e1(iv1)
      acd56(35)=spvae2e1(iv1)
      acd56(36)=dotproduct(qshift,spvae1e2)
      acd56(37)=abb56(4)
      acd56(38)=spvae2k1(iv2)
      acd56(39)=abb56(29)
      acd56(40)=spvae1k2(iv2)
      acd56(41)=abb56(34)
      acd56(42)=spvae1k5(iv2)
      acd56(43)=abb56(31)
      acd56(44)=spvae1e2(iv2)
      acd56(45)=abb56(15)
      acd56(46)=spvae2k1(iv1)
      acd56(47)=spvae1k2(iv1)
      acd56(48)=spvae1k5(iv1)
      acd56(49)=spvae1e2(iv1)
      acd56(50)=abb56(16)
      acd56(51)=abb56(19)
      acd56(52)=abb56(26)
      acd56(53)=abb56(14)
      acd56(54)=abb56(23)
      acd56(55)=abb56(28)
      acd56(56)=-acd56(19)*acd56(44)
      acd56(57)=-acd56(14)*acd56(49)
      acd56(56)=acd56(57)+acd56(56)
      acd56(56)=acd56(4)*acd56(56)
      acd56(57)=-acd56(20)*acd56(30)
      acd56(58)=-acd56(15)*acd56(35)
      acd56(57)=acd56(58)+acd56(57)
      acd56(57)=acd56(2)*acd56(57)
      acd56(58)=-acd56(6)*acd56(44)
      acd56(59)=-acd56(14)*acd56(36)
      acd56(58)=acd56(58)+acd56(59)
      acd56(58)=acd56(18)*acd56(58)
      acd56(59)=-acd56(8)*acd56(30)
      acd56(60)=-acd56(15)*acd56(21)
      acd56(59)=acd56(59)+acd56(60)
      acd56(59)=acd56(17)*acd56(59)
      acd56(60)=-acd56(6)*acd56(49)
      acd56(61)=-acd56(19)*acd56(36)
      acd56(60)=acd56(60)+acd56(61)
      acd56(60)=acd56(13)*acd56(60)
      acd56(61)=-acd56(8)*acd56(35)
      acd56(62)=-acd56(20)*acd56(21)
      acd56(61)=acd56(61)+acd56(62)
      acd56(61)=acd56(12)*acd56(61)
      acd56(56)=acd56(61)+acd56(60)+acd56(59)+acd56(58)+acd56(57)+acd56(56)
      acd56(56)=acd56(22)*acd56(56)
      acd56(57)=acd56(42)*acd56(43)
      acd56(58)=acd56(40)*acd56(41)
      acd56(59)=acd56(38)*acd56(39)
      acd56(60)=acd56(44)*acd56(45)
      acd56(61)=2.0_ki*acd56(16)
      acd56(62)=-acd56(5)*acd56(61)
      acd56(63)=acd56(14)*acd56(37)
      acd56(57)=acd56(63)+acd56(62)+acd56(60)+acd56(59)+acd56(57)+acd56(58)
      acd56(57)=acd56(18)*acd56(57)
      acd56(58)=acd56(28)*acd56(29)
      acd56(59)=acd56(26)*acd56(27)
      acd56(60)=acd56(24)*acd56(25)
      acd56(62)=acd56(30)*acd56(31)
      acd56(63)=-acd56(3)*acd56(61)
      acd56(64)=acd56(15)*acd56(23)
      acd56(58)=acd56(64)+acd56(63)+acd56(62)+acd56(60)+acd56(58)+acd56(59)
      acd56(58)=acd56(17)*acd56(58)
      acd56(59)=acd56(43)*acd56(48)
      acd56(60)=acd56(41)*acd56(47)
      acd56(62)=acd56(39)*acd56(46)
      acd56(63)=acd56(49)*acd56(45)
      acd56(64)=2.0_ki*acd56(11)
      acd56(65)=-acd56(5)*acd56(64)
      acd56(66)=acd56(19)*acd56(37)
      acd56(59)=acd56(66)+acd56(65)+acd56(63)+acd56(62)+acd56(59)+acd56(60)
      acd56(59)=acd56(13)*acd56(59)
      acd56(60)=acd56(29)*acd56(34)
      acd56(62)=acd56(27)*acd56(33)
      acd56(63)=acd56(25)*acd56(32)
      acd56(65)=acd56(35)*acd56(31)
      acd56(66)=-acd56(3)*acd56(64)
      acd56(67)=acd56(20)*acd56(23)
      acd56(60)=acd56(67)+acd56(66)+acd56(65)+acd56(63)+acd56(60)+acd56(62)
      acd56(60)=acd56(12)*acd56(60)
      acd56(62)=-acd56(8)*acd56(9)
      acd56(63)=-acd56(6)*acd56(7)
      acd56(65)=-acd56(4)*acd56(5)
      acd56(66)=-acd56(2)*acd56(3)
      acd56(62)=acd56(66)+acd56(65)+acd56(63)+acd56(10)+acd56(62)
      acd56(62)=acd56(1)*acd56(62)
      acd56(63)=acd56(42)*acd56(55)
      acd56(65)=acd56(40)*acd56(54)
      acd56(66)=acd56(38)*acd56(53)
      acd56(67)=-acd56(9)*acd56(61)
      acd56(63)=acd56(67)+acd56(66)+acd56(63)+acd56(65)
      acd56(63)=acd56(20)*acd56(63)
      acd56(65)=acd56(28)*acd56(52)
      acd56(66)=acd56(26)*acd56(51)
      acd56(67)=acd56(24)*acd56(50)
      acd56(61)=-acd56(7)*acd56(61)
      acd56(61)=acd56(61)+acd56(67)+acd56(65)+acd56(66)
      acd56(61)=acd56(19)*acd56(61)
      acd56(65)=acd56(48)*acd56(55)
      acd56(66)=acd56(47)*acd56(54)
      acd56(67)=acd56(46)*acd56(53)
      acd56(68)=-acd56(9)*acd56(64)
      acd56(65)=acd56(68)+acd56(67)+acd56(65)+acd56(66)
      acd56(65)=acd56(15)*acd56(65)
      acd56(66)=acd56(34)*acd56(52)
      acd56(67)=acd56(33)*acd56(51)
      acd56(68)=acd56(32)*acd56(50)
      acd56(64)=-acd56(7)*acd56(64)
      acd56(64)=acd56(64)+acd56(68)+acd56(66)+acd56(67)
      acd56(64)=acd56(14)*acd56(64)
      brack=acd56(56)+acd56(57)+acd56(58)+acd56(59)+acd56(60)+acd56(61)+2.0_ki*&
      &acd56(62)+acd56(63)+acd56(64)+acd56(65)
   end function brack_3
!---#] function brack_3:
!---#[ function brack_4:
   pure function brack_4(Q, mu2) result(brack)
      use p9_part21part21_part25part25part21_model_qp
      use p9_part21part21_part25part25part21_kinematics_qp
      use p9_part21part21_part25part25part21_color_qp
      use p9_part21part21_part25part25part21_abbrevd56h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(33) :: acd56
      complex(ki) :: brack
      acd56(1)=d(iv1,iv2)
      acd56(2)=spvae1e5(iv3)
      acd56(3)=abb56(25)
      acd56(4)=spvae5e1(iv3)
      acd56(5)=abb56(30)
      acd56(6)=spvae2e5(iv3)
      acd56(7)=abb56(27)
      acd56(8)=spvae5e2(iv3)
      acd56(9)=abb56(9)
      acd56(10)=d(iv1,iv3)
      acd56(11)=spvae1e5(iv2)
      acd56(12)=spvae5e1(iv2)
      acd56(13)=spvae2e5(iv2)
      acd56(14)=spvae5e2(iv2)
      acd56(15)=d(iv2,iv3)
      acd56(16)=spvae1e5(iv1)
      acd56(17)=spvae5e1(iv1)
      acd56(18)=spvae2e5(iv1)
      acd56(19)=spvae5e2(iv1)
      acd56(20)=spvae2e1(iv3)
      acd56(21)=abb56(7)
      acd56(22)=spvae2e1(iv2)
      acd56(23)=spvae2e1(iv1)
      acd56(24)=spvae1e2(iv3)
      acd56(25)=spvae1e2(iv2)
      acd56(26)=spvae1e2(iv1)
      acd56(27)=-acd56(13)*acd56(17)
      acd56(28)=-acd56(12)*acd56(18)
      acd56(27)=acd56(28)+acd56(27)
      acd56(27)=acd56(24)*acd56(27)
      acd56(28)=-acd56(14)*acd56(16)
      acd56(29)=-acd56(11)*acd56(19)
      acd56(28)=acd56(29)+acd56(28)
      acd56(28)=acd56(20)*acd56(28)
      acd56(29)=-acd56(16)*acd56(22)
      acd56(30)=-acd56(11)*acd56(23)
      acd56(29)=acd56(29)+acd56(30)
      acd56(29)=acd56(8)*acd56(29)
      acd56(30)=-acd56(17)*acd56(25)
      acd56(31)=-acd56(12)*acd56(26)
      acd56(30)=acd56(30)+acd56(31)
      acd56(30)=acd56(6)*acd56(30)
      acd56(31)=-acd56(18)*acd56(25)
      acd56(32)=-acd56(13)*acd56(26)
      acd56(31)=acd56(31)+acd56(32)
      acd56(31)=acd56(4)*acd56(31)
      acd56(32)=-acd56(19)*acd56(22)
      acd56(33)=-acd56(14)*acd56(23)
      acd56(32)=acd56(32)+acd56(33)
      acd56(32)=acd56(2)*acd56(32)
      acd56(27)=acd56(32)+acd56(31)+acd56(30)+acd56(29)+acd56(28)+acd56(27)
      acd56(27)=acd56(21)*acd56(27)
      acd56(28)=-acd56(9)*acd56(19)
      acd56(29)=-acd56(7)*acd56(18)
      acd56(30)=-acd56(5)*acd56(17)
      acd56(31)=-acd56(3)*acd56(16)
      acd56(28)=acd56(31)+acd56(30)+acd56(28)+acd56(29)
      acd56(28)=acd56(15)*acd56(28)
      acd56(29)=-acd56(9)*acd56(14)
      acd56(30)=-acd56(7)*acd56(13)
      acd56(31)=-acd56(5)*acd56(12)
      acd56(32)=-acd56(3)*acd56(11)
      acd56(29)=acd56(32)+acd56(31)+acd56(29)+acd56(30)
      acd56(29)=acd56(10)*acd56(29)
      acd56(30)=-acd56(8)*acd56(9)
      acd56(31)=-acd56(6)*acd56(7)
      acd56(32)=-acd56(4)*acd56(5)
      acd56(33)=-acd56(2)*acd56(3)
      acd56(30)=acd56(33)+acd56(32)+acd56(30)+acd56(31)
      acd56(30)=acd56(1)*acd56(30)
      acd56(28)=acd56(30)+acd56(28)+acd56(29)
      brack=acd56(27)+2.0_ki*acd56(28)
   end function brack_4
!---#] function brack_4:
!---#[ function brack_5:
   pure function brack_5(Q, mu2) result(brack)
      use p9_part21part21_part25part25part21_model_qp
      use p9_part21part21_part25part25part21_kinematics_qp
      use p9_part21part21_part25part25part21_color_qp
      use p9_part21part21_part25part25part21_abbrevd56h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(1) :: acd56
      complex(ki) :: brack
      brack=0.0_ki
   end function brack_5
!---#] function brack_5:
!---#[ function derivative:
   function derivative(mu2,i1,i2,i3,i4) result(numerator)
      use p9_part21part21_part25part25part21_globalsl1_qp, only: epspow
      use p9_part21part21_part25part25part21_kinematics_qp
      use p9_part21part21_part25part25part21_abbrevd56h0_qp
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
      qshift = -k2+k5
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
end module     p9_part21part21_part25part25part21_d56h0l1d_qp
