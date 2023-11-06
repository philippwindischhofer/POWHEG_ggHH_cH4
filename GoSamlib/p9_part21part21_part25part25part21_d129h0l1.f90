module     p9_part21part21_part25part25part21_d129h0l1
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p9_part21part &
   ! &21_part25part25part21/helicity0d129h0l1.f90
   ! generator: buildfortran.py
   use p9_part21part21_part25part25part21_config, only: ki
   use p9_part21part21_part25part25part21_util, only: cond
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   public :: numerator_ninja
contains
!---#[ function brack_1:
   pure function brack_1(Q,mu2) result(brack)
      use p9_part21part21_part25part25part21_model
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_color
      use p9_part21part21_part25part25part21_abbrevd129h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc129(115)
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
      complex(ki) :: Qspvak1l4
      complex(ki) :: Qspvak5l4
      complex(ki) :: Qspvae1l4
      complex(ki) :: Qspvae5l4
      complex(ki) :: Qspval4e2
      complex(ki) :: Qspval4k1
      complex(ki) :: Qspval4k5
      complex(ki) :: Qspval4e1
      complex(ki) :: Qspval4e5
      complex(ki) :: Qspvae2l4
      complex(ki) :: Qspval3e2
      complex(ki) :: Qspvae2l3
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
      Qspvak1l4 = dotproduct(Q,spvak1l4)
      Qspvak5l4 = dotproduct(Q,spvak5l4)
      Qspvae1l4 = dotproduct(Q,spvae1l4)
      Qspvae5l4 = dotproduct(Q,spvae5l4)
      Qspval4e2 = dotproduct(Q,spval4e2)
      Qspval4k1 = dotproduct(Q,spval4k1)
      Qspval4k5 = dotproduct(Q,spval4k5)
      Qspval4e1 = dotproduct(Q,spval4e1)
      Qspval4e5 = dotproduct(Q,spval4e5)
      Qspvae2l4 = dotproduct(Q,spvae2l4)
      Qspval3e2 = dotproduct(Q,spval3e2)
      Qspvae2l3 = dotproduct(Q,spvae2l3)
      acc129(1)=abb129(7)
      acc129(2)=abb129(8)
      acc129(3)=abb129(9)
      acc129(4)=abb129(10)
      acc129(5)=abb129(11)
      acc129(6)=abb129(12)
      acc129(7)=abb129(13)
      acc129(8)=abb129(15)
      acc129(9)=abb129(16)
      acc129(10)=abb129(17)
      acc129(11)=abb129(18)
      acc129(12)=abb129(19)
      acc129(13)=abb129(20)
      acc129(14)=abb129(21)
      acc129(15)=abb129(22)
      acc129(16)=abb129(23)
      acc129(17)=abb129(24)
      acc129(18)=abb129(25)
      acc129(19)=abb129(26)
      acc129(20)=abb129(27)
      acc129(21)=abb129(28)
      acc129(22)=abb129(29)
      acc129(23)=abb129(30)
      acc129(24)=abb129(31)
      acc129(25)=abb129(32)
      acc129(26)=abb129(33)
      acc129(27)=abb129(34)
      acc129(28)=abb129(35)
      acc129(29)=abb129(36)
      acc129(30)=abb129(37)
      acc129(31)=abb129(38)
      acc129(32)=abb129(39)
      acc129(33)=abb129(40)
      acc129(34)=abb129(41)
      acc129(35)=abb129(42)
      acc129(36)=abb129(43)
      acc129(37)=abb129(44)
      acc129(38)=abb129(45)
      acc129(39)=abb129(46)
      acc129(40)=abb129(47)
      acc129(41)=abb129(48)
      acc129(42)=abb129(49)
      acc129(43)=abb129(50)
      acc129(44)=abb129(51)
      acc129(45)=abb129(52)
      acc129(46)=abb129(53)
      acc129(47)=abb129(54)
      acc129(48)=abb129(55)
      acc129(49)=abb129(56)
      acc129(50)=abb129(57)
      acc129(51)=abb129(58)
      acc129(52)=abb129(59)
      acc129(53)=abb129(60)
      acc129(54)=abb129(61)
      acc129(55)=abb129(62)
      acc129(56)=abb129(63)
      acc129(57)=abb129(64)
      acc129(58)=abb129(65)
      acc129(59)=abb129(66)
      acc129(60)=abb129(67)
      acc129(61)=abb129(68)
      acc129(62)=abb129(69)
      acc129(63)=abb129(70)
      acc129(64)=abb129(71)
      acc129(65)=abb129(72)
      acc129(66)=abb129(73)
      acc129(67)=abb129(74)
      acc129(68)=abb129(75)
      acc129(69)=abb129(76)
      acc129(70)=abb129(77)
      acc129(71)=abb129(78)
      acc129(72)=abb129(79)
      acc129(73)=abb129(80)
      acc129(74)=abb129(81)
      acc129(75)=abb129(82)
      acc129(76)=abb129(84)
      acc129(77)=abb129(85)
      acc129(78)=abb129(87)
      acc129(79)=abb129(90)
      acc129(80)=abb129(94)
      acc129(81)=abb129(97)
      acc129(82)=abb129(98)
      acc129(83)=abb129(102)
      acc129(84)=abb129(105)
      acc129(85)=Qspvak1e2*acc129(24)
      acc129(86)=Qspvak1e2*acc129(12)
      acc129(86)=acc129(29)+acc129(86)
      acc129(86)=Qspvae2k1*acc129(86)
      acc129(87)=Qspvak5e2*acc129(81)
      acc129(88)=Qspvak5e2*acc129(12)
      acc129(88)=acc129(80)+acc129(88)
      acc129(88)=Qspvae2k5*acc129(88)
      acc129(89)=Qspvae1e2*acc129(62)
      acc129(90)=Qspvae1e2*acc129(75)
      acc129(90)=acc129(8)+acc129(90)
      acc129(90)=Qspvae2e1*acc129(90)
      acc129(91)=Qspvae2e5*acc129(44)
      acc129(92)=Qspvae2e5*acc129(72)
      acc129(92)=acc129(70)+acc129(92)
      acc129(92)=Qspvae5e2*acc129(92)
      acc129(85)=acc129(92)+acc129(91)+acc129(90)+acc129(89)+acc129(88)+acc129(&
      &87)+acc129(86)+acc129(10)+acc129(85)
      acc129(85)=QspQ*acc129(85)
      acc129(86)=Qspvak1k2*acc129(73)
      acc129(87)=Qspvak5k2*acc129(69)
      acc129(88)=Qspvae1k2*acc129(54)
      acc129(89)=Qspvae5k2*acc129(55)
      acc129(90)=acc129(15)*Qspvak1k2
      acc129(90)=acc129(66)+acc129(90)
      acc129(90)=Qspvae2k1*acc129(90)
      acc129(91)=acc129(15)*Qspvak5k2
      acc129(91)=acc129(33)+acc129(91)
      acc129(91)=Qspvae2k5*acc129(91)
      acc129(92)=-Qspvae1k2*acc129(3)
      acc129(92)=acc129(30)+acc129(92)
      acc129(92)=Qspvae2e1*acc129(92)
      acc129(93)=-Qspvae5k2*acc129(59)
      acc129(93)=acc129(5)+acc129(93)
      acc129(93)=Qspvae2e5*acc129(93)
      acc129(86)=acc129(93)+acc129(92)+acc129(91)+acc129(90)+acc129(89)+acc129(&
      &88)+acc129(87)+acc129(17)+acc129(86)
      acc129(86)=Qspvak2e2*acc129(86)
      acc129(87)=Qspvak2k1*acc129(46)
      acc129(88)=Qspvak2k5*acc129(36)
      acc129(89)=Qspvak2e1*acc129(47)
      acc129(90)=Qspvak2e5*acc129(48)
      acc129(91)=acc129(15)*Qspvak2k1
      acc129(91)=acc129(27)+acc129(91)
      acc129(91)=Qspvak1e2*acc129(91)
      acc129(92)=acc129(15)*Qspvak2k5
      acc129(92)=acc129(42)+acc129(92)
      acc129(92)=Qspvak5e2*acc129(92)
      acc129(93)=-Qspvak2e1*acc129(3)
      acc129(93)=acc129(37)+acc129(93)
      acc129(93)=Qspvae1e2*acc129(93)
      acc129(94)=-Qspvak2e5*acc129(59)
      acc129(94)=acc129(34)+acc129(94)
      acc129(94)=Qspvae5e2*acc129(94)
      acc129(87)=acc129(94)+acc129(93)+acc129(92)+acc129(91)+acc129(90)+acc129(&
      &89)+acc129(88)+acc129(23)+acc129(87)
      acc129(87)=Qspvae2k2*acc129(87)
      acc129(88)=Qspvak1l4*acc129(28)
      acc129(89)=-Qspvak5l4*acc129(20)
      acc129(90)=Qspvae1l4*acc129(64)
      acc129(91)=Qspvae5l4*acc129(49)
      acc129(92)=acc129(71)*Qspvak1l4
      acc129(92)=acc129(65)+acc129(92)
      acc129(92)=Qspvae2k1*acc129(92)
      acc129(93)=acc129(71)*Qspvak5l4
      acc129(93)=acc129(84)+acc129(93)
      acc129(93)=Qspvae2k5*acc129(93)
      acc129(94)=Qspvae1l4*acc129(9)
      acc129(94)=acc129(76)+acc129(94)
      acc129(94)=Qspvae2e1*acc129(94)
      acc129(95)=Qspvae5l4*acc129(6)
      acc129(95)=acc129(45)+acc129(95)
      acc129(95)=Qspvae2e5*acc129(95)
      acc129(88)=acc129(95)+acc129(94)+acc129(93)+acc129(92)+acc129(91)+acc129(&
      &90)+acc129(89)+acc129(1)+acc129(88)
      acc129(88)=Qspval4e2*acc129(88)
      acc129(89)=Qspval4k1*acc129(57)
      acc129(90)=-Qspval4k5*acc129(50)
      acc129(91)=Qspval4e1*acc129(56)
      acc129(92)=Qspval4e5*acc129(39)
      acc129(93)=acc129(71)*Qspval4k1
      acc129(93)=acc129(2)+acc129(93)
      acc129(93)=Qspvak1e2*acc129(93)
      acc129(94)=acc129(71)*Qspval4k5
      acc129(94)=acc129(83)+acc129(94)
      acc129(94)=Qspvak5e2*acc129(94)
      acc129(95)=Qspval4e1*acc129(9)
      acc129(95)=acc129(11)+acc129(95)
      acc129(95)=Qspvae1e2*acc129(95)
      acc129(96)=Qspval4e5*acc129(6)
      acc129(96)=acc129(79)+acc129(96)
      acc129(96)=Qspvae5e2*acc129(96)
      acc129(89)=acc129(96)+acc129(95)+acc129(94)+acc129(93)+acc129(92)+acc129(&
      &91)+acc129(90)+acc129(38)+acc129(89)
      acc129(89)=Qspvae2l4*acc129(89)
      acc129(90)=acc129(21)*Qspval3e2
      acc129(91)=acc129(4)*Qspvae2l3
      acc129(92)=Qspvak1k2*acc129(13)
      acc129(93)=Qspvak1l4*acc129(67)
      acc129(94)=Qspvak2k1*acc129(19)
      acc129(95)=Qspvak2k5*acc129(41)
      acc129(96)=Qspval4k1*acc129(16)
      acc129(97)=Qspval4k5*acc129(25)
      acc129(98)=Qspvak5k2*acc129(74)
      acc129(99)=Qspvak5l4*acc129(68)
      acc129(100)=Qspvak2e1*acc129(63)
      acc129(101)=Qspvae1k2*acc129(35)
      acc129(102)=Qspvak2e5*acc129(31)
      acc129(103)=Qspvae5k2*acc129(14)
      acc129(104)=-Qspval4e1*acc129(60)
      acc129(105)=-Qspvae1l4*acc129(51)
      acc129(106)=Qspval4e5*acc129(43)
      acc129(107)=Qspvae5l4*acc129(52)
      acc129(108)=Qspvak1e2*acc129(26)
      acc129(109)=Qspvak1e2*acc129(82)
      acc129(109)=acc129(22)+acc129(109)
      acc129(109)=Qspvae2k1*acc129(109)
      acc129(110)=Qspvak5e2*acc129(77)
      acc129(111)=Qspvak5e2*acc129(82)
      acc129(111)=acc129(78)+acc129(111)
      acc129(111)=Qspvae2k5*acc129(111)
      acc129(112)=Qspvae1e2*acc129(18)
      acc129(113)=Qspvae1e2*acc129(53)
      acc129(113)=acc129(61)+acc129(113)
      acc129(113)=Qspvae2e1*acc129(113)
      acc129(114)=Qspvae2e5*acc129(40)
      acc129(115)=Qspvae2e5*acc129(32)
      acc129(115)=acc129(58)+acc129(115)
      acc129(115)=Qspvae5e2*acc129(115)
      brack=acc129(7)+acc129(85)+acc129(86)+acc129(87)+acc129(88)+acc129(89)+ac&
      &c129(90)+acc129(91)+acc129(92)+acc129(93)+acc129(94)+acc129(95)+acc129(9&
      &6)+acc129(97)+acc129(98)+acc129(99)+acc129(100)+acc129(101)+acc129(102)+&
      &acc129(103)+acc129(104)+acc129(105)+acc129(106)+acc129(107)+acc129(108)+&
      &acc129(109)+acc129(110)+acc129(111)+acc129(112)+acc129(113)+acc129(114)+&
      &acc129(115)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ subroutine numerator_ninja:
   subroutine numerator_ninja(ncut, Q_ext, mu2_ext, numerator) &
   & bind(c, name="p9_part21part21_part25part25part21_d129h0l1_ninja")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki_nin
      use p9_part21part21_part25part25part21_globalsl1, only: epspow
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_abbrevd129h0
      implicit none
      integer(c_int), intent(in) :: ncut
      complex(ki_nin), dimension(0:3), intent(in) :: Q_ext
      complex(ki_nin), intent(in) :: mu2_ext
      complex(ki_nin), intent(out) :: numerator
      complex(ki) :: d129
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      Q(1:4)  =cmplx(real(-Q_ext(0:3),  ki_nin), aimag(-Q_ext(0:3)), ki)
      d129 = 0.0_ki
      d129 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d129, ki), aimag(d129), ki_nin)
   end subroutine numerator_ninja
   !------#] subroutine numerator_ninja:
!---#] numerator interfaces:
end module p9_part21part21_part25part25part21_d129h0l1
