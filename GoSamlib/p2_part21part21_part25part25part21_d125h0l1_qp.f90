module     p2_part21part21_part25part25part21_d125h0l1_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d125h0l1_qp.f90
   ! generator: buildfortran.py
   use p2_part21part21_part25part25part21_config, only: ki => ki_qp
   use p2_part21part21_part25part25part21_util_qp, only: cond
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   public :: numerator_ninja
contains
!---#[ function brack_1:
   pure function brack_1(Q,mu2) result(brack)
      use p2_part21part21_part25part25part21_model_qp
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_color_qp
      use p2_part21part21_part25part25part21_abbrevd125h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc125(115)
      complex(ki) :: Qspvak1e2
      complex(ki) :: Qspvae2k1
      complex(ki) :: Qspvak5e2
      complex(ki) :: Qspvae2k5
      complex(ki) :: Qspvae1e2
      complex(ki) :: Qspvae2e1
      complex(ki) :: Qspvae2e5
      complex(ki) :: Qspvae5e2
      complex(ki) :: QspQ
      complex(ki) :: Qspvak1k2
      complex(ki) :: Qspvak5k2
      complex(ki) :: Qspvae1k2
      complex(ki) :: Qspvae5k2
      complex(ki) :: Qspvak2e2
      complex(ki) :: Qspvak2k1
      complex(ki) :: Qspvak2k5
      complex(ki) :: Qspvak2e1
      complex(ki) :: Qspvak2e5
      complex(ki) :: Qspvae2k2
      complex(ki) :: Qspvak1l3
      complex(ki) :: Qspvak5l3
      complex(ki) :: Qspvae1l3
      complex(ki) :: Qspvae5l3
      complex(ki) :: Qspval3e2
      complex(ki) :: Qspval3k1
      complex(ki) :: Qspval3k5
      complex(ki) :: Qspval3e1
      complex(ki) :: Qspval3e5
      complex(ki) :: Qspvae2l3
      complex(ki) :: Qspval4e2
      complex(ki) :: Qspvae2l4
      Qspvak1e2 = dotproduct(Q,spvak1e2)
      Qspvae2k1 = dotproduct(Q,spvae2k1)
      Qspvak5e2 = dotproduct(Q,spvak5e2)
      Qspvae2k5 = dotproduct(Q,spvae2k5)
      Qspvae1e2 = dotproduct(Q,spvae1e2)
      Qspvae2e1 = dotproduct(Q,spvae2e1)
      Qspvae2e5 = dotproduct(Q,spvae2e5)
      Qspvae5e2 = dotproduct(Q,spvae5e2)
      QspQ = dotproduct(Q,Q)
      Qspvak1k2 = dotproduct(Q,spvak1k2)
      Qspvak5k2 = dotproduct(Q,spvak5k2)
      Qspvae1k2 = dotproduct(Q,spvae1k2)
      Qspvae5k2 = dotproduct(Q,spvae5k2)
      Qspvak2e2 = dotproduct(Q,spvak2e2)
      Qspvak2k1 = dotproduct(Q,spvak2k1)
      Qspvak2k5 = dotproduct(Q,spvak2k5)
      Qspvak2e1 = dotproduct(Q,spvak2e1)
      Qspvak2e5 = dotproduct(Q,spvak2e5)
      Qspvae2k2 = dotproduct(Q,spvae2k2)
      Qspvak1l3 = dotproduct(Q,spvak1l3)
      Qspvak5l3 = dotproduct(Q,spvak5l3)
      Qspvae1l3 = dotproduct(Q,spvae1l3)
      Qspvae5l3 = dotproduct(Q,spvae5l3)
      Qspval3e2 = dotproduct(Q,spval3e2)
      Qspval3k1 = dotproduct(Q,spval3k1)
      Qspval3k5 = dotproduct(Q,spval3k5)
      Qspval3e1 = dotproduct(Q,spval3e1)
      Qspval3e5 = dotproduct(Q,spval3e5)
      Qspvae2l3 = dotproduct(Q,spvae2l3)
      Qspval4e2 = dotproduct(Q,spval4e2)
      Qspvae2l4 = dotproduct(Q,spvae2l4)
      acc125(1)=abb125(7)
      acc125(2)=abb125(8)
      acc125(3)=abb125(9)
      acc125(4)=abb125(10)
      acc125(5)=abb125(11)
      acc125(6)=abb125(12)
      acc125(7)=abb125(13)
      acc125(8)=abb125(15)
      acc125(9)=abb125(16)
      acc125(10)=abb125(17)
      acc125(11)=abb125(18)
      acc125(12)=abb125(19)
      acc125(13)=abb125(20)
      acc125(14)=abb125(21)
      acc125(15)=abb125(22)
      acc125(16)=abb125(23)
      acc125(17)=abb125(24)
      acc125(18)=abb125(25)
      acc125(19)=abb125(26)
      acc125(20)=abb125(27)
      acc125(21)=abb125(28)
      acc125(22)=abb125(29)
      acc125(23)=abb125(30)
      acc125(24)=abb125(31)
      acc125(25)=abb125(32)
      acc125(26)=abb125(33)
      acc125(27)=abb125(34)
      acc125(28)=abb125(35)
      acc125(29)=abb125(36)
      acc125(30)=abb125(37)
      acc125(31)=abb125(38)
      acc125(32)=abb125(39)
      acc125(33)=abb125(40)
      acc125(34)=abb125(41)
      acc125(35)=abb125(42)
      acc125(36)=abb125(43)
      acc125(37)=abb125(44)
      acc125(38)=abb125(45)
      acc125(39)=abb125(46)
      acc125(40)=abb125(47)
      acc125(41)=abb125(48)
      acc125(42)=abb125(49)
      acc125(43)=abb125(50)
      acc125(44)=abb125(51)
      acc125(45)=abb125(52)
      acc125(46)=abb125(53)
      acc125(47)=abb125(54)
      acc125(48)=abb125(55)
      acc125(49)=abb125(56)
      acc125(50)=abb125(57)
      acc125(51)=abb125(58)
      acc125(52)=abb125(59)
      acc125(53)=abb125(60)
      acc125(54)=abb125(61)
      acc125(55)=abb125(62)
      acc125(56)=abb125(63)
      acc125(57)=abb125(64)
      acc125(58)=abb125(65)
      acc125(59)=abb125(66)
      acc125(60)=abb125(67)
      acc125(61)=abb125(68)
      acc125(62)=abb125(69)
      acc125(63)=abb125(70)
      acc125(64)=abb125(71)
      acc125(65)=abb125(72)
      acc125(66)=abb125(73)
      acc125(67)=abb125(74)
      acc125(68)=abb125(75)
      acc125(69)=abb125(76)
      acc125(70)=abb125(77)
      acc125(71)=abb125(78)
      acc125(72)=abb125(79)
      acc125(73)=abb125(80)
      acc125(74)=abb125(81)
      acc125(75)=abb125(82)
      acc125(76)=abb125(83)
      acc125(77)=abb125(84)
      acc125(78)=abb125(85)
      acc125(79)=abb125(87)
      acc125(80)=abb125(90)
      acc125(81)=abb125(94)
      acc125(82)=abb125(97)
      acc125(83)=abb125(98)
      acc125(84)=abb125(101)
      acc125(85)=Qspvak1e2*acc125(24)
      acc125(86)=Qspvak1e2*acc125(12)
      acc125(86)=acc125(29)+acc125(86)
      acc125(86)=Qspvae2k1*acc125(86)
      acc125(87)=Qspvak5e2*acc125(82)
      acc125(88)=Qspvak5e2*acc125(12)
      acc125(88)=acc125(81)+acc125(88)
      acc125(88)=Qspvae2k5*acc125(88)
      acc125(89)=Qspvae1e2*acc125(62)
      acc125(90)=Qspvae1e2*acc125(75)
      acc125(90)=acc125(8)+acc125(90)
      acc125(90)=Qspvae2e1*acc125(90)
      acc125(91)=Qspvae2e5*acc125(44)
      acc125(92)=Qspvae2e5*acc125(72)
      acc125(92)=acc125(70)+acc125(92)
      acc125(92)=Qspvae5e2*acc125(92)
      acc125(85)=acc125(92)+acc125(91)+acc125(90)+acc125(89)+acc125(88)+acc125(&
      &87)+acc125(86)+acc125(10)+acc125(85)
      acc125(85)=QspQ*acc125(85)
      acc125(86)=Qspvak1k2*acc125(77)
      acc125(87)=Qspvak5k2*acc125(74)
      acc125(88)=Qspvae1k2*acc125(54)
      acc125(89)=Qspvae5k2*acc125(34)
      acc125(90)=acc125(15)*Qspvak1k2
      acc125(90)=acc125(63)+acc125(90)
      acc125(90)=Qspvae2k1*acc125(90)
      acc125(91)=acc125(15)*Qspvak5k2
      acc125(91)=acc125(33)+acc125(91)
      acc125(91)=Qspvae2k5*acc125(91)
      acc125(92)=-Qspvae1k2*acc125(3)
      acc125(92)=acc125(30)+acc125(92)
      acc125(92)=Qspvae2e1*acc125(92)
      acc125(93)=-Qspvae5k2*acc125(65)
      acc125(93)=acc125(5)+acc125(93)
      acc125(93)=Qspvae2e5*acc125(93)
      acc125(86)=acc125(93)+acc125(92)+acc125(91)+acc125(90)+acc125(89)+acc125(&
      &88)+acc125(87)+acc125(17)+acc125(86)
      acc125(86)=Qspvak2e2*acc125(86)
      acc125(87)=Qspvak2k1*acc125(76)
      acc125(88)=Qspvak2k5*acc125(71)
      acc125(89)=Qspvak2e1*acc125(47)
      acc125(90)=Qspvak2e5*acc125(42)
      acc125(91)=acc125(15)*Qspvak2k1
      acc125(91)=acc125(27)+acc125(91)
      acc125(91)=Qspvak1e2*acc125(91)
      acc125(92)=acc125(15)*Qspvak2k5
      acc125(92)=acc125(48)+acc125(92)
      acc125(92)=Qspvak5e2*acc125(92)
      acc125(93)=-Qspvak2e1*acc125(3)
      acc125(93)=acc125(37)+acc125(93)
      acc125(93)=Qspvae1e2*acc125(93)
      acc125(94)=-Qspvak2e5*acc125(65)
      acc125(94)=acc125(36)+acc125(94)
      acc125(94)=Qspvae5e2*acc125(94)
      acc125(87)=acc125(94)+acc125(93)+acc125(92)+acc125(91)+acc125(90)+acc125(&
      &89)+acc125(88)+acc125(23)+acc125(87)
      acc125(87)=Qspvae2k2*acc125(87)
      acc125(88)=Qspvak1l3*acc125(28)
      acc125(89)=-Qspvak5l3*acc125(20)
      acc125(90)=Qspvae1l3*acc125(64)
      acc125(91)=Qspvae5l3*acc125(49)
      acc125(92)=acc125(68)*Qspvak1l3
      acc125(92)=acc125(59)+acc125(92)
      acc125(92)=Qspvae2k1*acc125(92)
      acc125(93)=acc125(68)*Qspvak5l3
      acc125(93)=acc125(84)+acc125(93)
      acc125(93)=Qspvae2k5*acc125(93)
      acc125(94)=Qspvae1l3*acc125(9)
      acc125(94)=acc125(73)+acc125(94)
      acc125(94)=Qspvae2e1*acc125(94)
      acc125(95)=Qspvae5l3*acc125(6)
      acc125(95)=acc125(45)+acc125(95)
      acc125(95)=Qspvae2e5*acc125(95)
      acc125(88)=acc125(95)+acc125(94)+acc125(93)+acc125(92)+acc125(91)+acc125(&
      &90)+acc125(89)+acc125(1)+acc125(88)
      acc125(88)=Qspval3e2*acc125(88)
      acc125(89)=Qspval3k1*acc125(57)
      acc125(90)=-Qspval3k5*acc125(50)
      acc125(91)=Qspval3e1*acc125(56)
      acc125(92)=Qspval3e5*acc125(39)
      acc125(93)=acc125(68)*Qspval3k1
      acc125(93)=acc125(2)+acc125(93)
      acc125(93)=Qspvak1e2*acc125(93)
      acc125(94)=acc125(68)*Qspval3k5
      acc125(94)=acc125(80)+acc125(94)
      acc125(94)=Qspvak5e2*acc125(94)
      acc125(95)=Qspval3e1*acc125(9)
      acc125(95)=acc125(11)+acc125(95)
      acc125(95)=Qspvae1e2*acc125(95)
      acc125(96)=Qspval3e5*acc125(6)
      acc125(96)=acc125(38)+acc125(96)
      acc125(96)=Qspvae5e2*acc125(96)
      acc125(89)=acc125(96)+acc125(95)+acc125(94)+acc125(93)+acc125(92)+acc125(&
      &91)+acc125(90)+acc125(31)+acc125(89)
      acc125(89)=Qspvae2l3*acc125(89)
      acc125(90)=acc125(21)*Qspval4e2
      acc125(91)=acc125(4)*Qspvae2l4
      acc125(92)=Qspvak1k2*acc125(13)
      acc125(93)=Qspvak1l3*acc125(66)
      acc125(94)=Qspvak2k1*acc125(19)
      acc125(95)=Qspvak2k5*acc125(41)
      acc125(96)=Qspval3k1*acc125(16)
      acc125(97)=Qspval3k5*acc125(25)
      acc125(98)=Qspvak5k2*acc125(69)
      acc125(99)=Qspvak5l3*acc125(67)
      acc125(100)=Qspvak2e1*acc125(55)
      acc125(101)=Qspvae1k2*acc125(35)
      acc125(102)=Qspvak2e5*acc125(46)
      acc125(103)=Qspvae5k2*acc125(14)
      acc125(104)=-Qspval3e1*acc125(60)
      acc125(105)=-Qspvae1l3*acc125(51)
      acc125(106)=Qspval3e5*acc125(43)
      acc125(107)=Qspvae5l3*acc125(52)
      acc125(108)=Qspvak1e2*acc125(26)
      acc125(109)=Qspvak1e2*acc125(83)
      acc125(109)=acc125(22)+acc125(109)
      acc125(109)=Qspvae2k1*acc125(109)
      acc125(110)=Qspvak5e2*acc125(78)
      acc125(111)=Qspvak5e2*acc125(83)
      acc125(111)=acc125(79)+acc125(111)
      acc125(111)=Qspvae2k5*acc125(111)
      acc125(112)=Qspvae1e2*acc125(18)
      acc125(113)=Qspvae1e2*acc125(53)
      acc125(113)=acc125(61)+acc125(113)
      acc125(113)=Qspvae2e1*acc125(113)
      acc125(114)=Qspvae2e5*acc125(40)
      acc125(115)=Qspvae2e5*acc125(32)
      acc125(115)=acc125(58)+acc125(115)
      acc125(115)=Qspvae5e2*acc125(115)
      brack=acc125(7)+acc125(85)+acc125(86)+acc125(87)+acc125(88)+acc125(89)+ac&
      &c125(90)+acc125(91)+acc125(92)+acc125(93)+acc125(94)+acc125(95)+acc125(9&
      &6)+acc125(97)+acc125(98)+acc125(99)+acc125(100)+acc125(101)+acc125(102)+&
      &acc125(103)+acc125(104)+acc125(105)+acc125(106)+acc125(107)+acc125(108)+&
      &acc125(109)+acc125(110)+acc125(111)+acc125(112)+acc125(113)+acc125(114)+&
      &acc125(115)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ subroutine numerator_ninja:
   subroutine numerator_ninja(ncut, Q_ext, mu2_ext, numerator) &
   & bind(c, name="p2_part21part21_part25part25part21_d125h0l1_qp_ninja")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki_nin
      use p2_part21part21_part25part25part21_globalsl1_qp, only: epspow
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_abbrevd125h0_qp
      implicit none
      integer(c_int), intent(in) :: ncut
      complex(ki_nin), dimension(0:3), intent(in) :: Q_ext
      complex(ki_nin), intent(in) :: mu2_ext
      complex(ki_nin), intent(out) :: numerator
      complex(ki) :: d125
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      Q(1:4)  =cmplx(real(-Q_ext(0:3),  ki_nin), aimag(-Q_ext(0:3)), ki)
      d125 = 0.0_ki
      d125 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d125, ki), aimag(d125), ki_nin)
   end subroutine numerator_ninja
   !------#] subroutine numerator_ninja:
!---#] numerator interfaces:
end module p2_part21part21_part25part25part21_d125h0l1_qp
