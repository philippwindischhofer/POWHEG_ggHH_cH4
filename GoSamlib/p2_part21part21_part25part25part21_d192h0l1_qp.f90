module     p2_part21part21_part25part25part21_d192h0l1_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d192h0l1_qp.f90
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
      use p2_part21part21_part25part25part21_abbrevd192h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc192(108)
      complex(ki) :: Qspvak1k2
      complex(ki) :: Qspvak5k2
      complex(ki) :: Qspvae2e5
      complex(ki) :: Qspvak2e5
      complex(ki) :: Qspvae2k2
      complex(ki) :: Qspvae2l3
      complex(ki) :: Qspval3e5
      complex(ki) :: QspQ
      complex(ki) :: Qspvae5e1
      complex(ki) :: Qspvae2k1
      complex(ki) :: Qspvae2k5
      complex(ki) :: Qspvae5k2
      complex(ki) :: Qspvae5l3
      complex(ki) :: Qspvae1e2
      complex(ki) :: Qspvak2k1
      complex(ki) :: Qspvak2k5
      complex(ki) :: Qspvae5e2
      complex(ki) :: Qspvak2e2
      complex(ki) :: Qspval3e2
      complex(ki) :: Qspvae1e5
      complex(ki) :: Qspvak1e2
      complex(ki) :: Qspvak5e2
      complex(ki) :: Qspvae2e1
      Qspvak1k2 = dotproduct(Q,spvak1k2)
      Qspvak5k2 = dotproduct(Q,spvak5k2)
      Qspvae2e5 = dotproduct(Q,spvae2e5)
      Qspvak2e5 = dotproduct(Q,spvak2e5)
      Qspvae2k2 = dotproduct(Q,spvae2k2)
      Qspvae2l3 = dotproduct(Q,spvae2l3)
      Qspval3e5 = dotproduct(Q,spval3e5)
      QspQ = dotproduct(Q,Q)
      Qspvae5e1 = dotproduct(Q,spvae5e1)
      Qspvae2k1 = dotproduct(Q,spvae2k1)
      Qspvae2k5 = dotproduct(Q,spvae2k5)
      Qspvae5k2 = dotproduct(Q,spvae5k2)
      Qspvae5l3 = dotproduct(Q,spvae5l3)
      Qspvae1e2 = dotproduct(Q,spvae1e2)
      Qspvak2k1 = dotproduct(Q,spvak2k1)
      Qspvak2k5 = dotproduct(Q,spvak2k5)
      Qspvae5e2 = dotproduct(Q,spvae5e2)
      Qspvak2e2 = dotproduct(Q,spvak2e2)
      Qspval3e2 = dotproduct(Q,spval3e2)
      Qspvae1e5 = dotproduct(Q,spvae1e5)
      Qspvak1e2 = dotproduct(Q,spvak1e2)
      Qspvak5e2 = dotproduct(Q,spvak5e2)
      Qspvae2e1 = dotproduct(Q,spvae2e1)
      acc192(1)=abb192(4)
      acc192(2)=abb192(5)
      acc192(3)=abb192(6)
      acc192(4)=abb192(7)
      acc192(5)=abb192(8)
      acc192(6)=abb192(9)
      acc192(7)=abb192(10)
      acc192(8)=abb192(11)
      acc192(9)=abb192(12)
      acc192(10)=abb192(13)
      acc192(11)=abb192(14)
      acc192(12)=abb192(15)
      acc192(13)=abb192(16)
      acc192(14)=abb192(17)
      acc192(15)=abb192(18)
      acc192(16)=abb192(19)
      acc192(17)=abb192(20)
      acc192(18)=abb192(21)
      acc192(19)=abb192(22)
      acc192(20)=abb192(23)
      acc192(21)=abb192(24)
      acc192(22)=abb192(25)
      acc192(23)=abb192(26)
      acc192(24)=abb192(27)
      acc192(25)=abb192(28)
      acc192(26)=abb192(29)
      acc192(27)=abb192(30)
      acc192(28)=abb192(31)
      acc192(29)=abb192(32)
      acc192(30)=abb192(33)
      acc192(31)=abb192(34)
      acc192(32)=abb192(35)
      acc192(33)=abb192(36)
      acc192(34)=abb192(37)
      acc192(35)=abb192(38)
      acc192(36)=abb192(39)
      acc192(37)=abb192(40)
      acc192(38)=abb192(42)
      acc192(39)=abb192(43)
      acc192(40)=abb192(44)
      acc192(41)=abb192(45)
      acc192(42)=abb192(46)
      acc192(43)=abb192(47)
      acc192(44)=abb192(48)
      acc192(45)=abb192(49)
      acc192(46)=abb192(50)
      acc192(47)=abb192(51)
      acc192(48)=abb192(52)
      acc192(49)=abb192(53)
      acc192(50)=abb192(54)
      acc192(51)=abb192(55)
      acc192(52)=abb192(56)
      acc192(53)=abb192(57)
      acc192(54)=abb192(58)
      acc192(55)=abb192(59)
      acc192(56)=abb192(60)
      acc192(57)=abb192(62)
      acc192(58)=abb192(63)
      acc192(59)=abb192(64)
      acc192(60)=abb192(65)
      acc192(61)=abb192(66)
      acc192(62)=abb192(67)
      acc192(63)=abb192(68)
      acc192(64)=abb192(69)
      acc192(65)=abb192(70)
      acc192(66)=abb192(71)
      acc192(67)=abb192(72)
      acc192(68)=abb192(73)
      acc192(69)=abb192(74)
      acc192(70)=abb192(75)
      acc192(71)=abb192(76)
      acc192(72)=abb192(77)
      acc192(73)=abb192(78)
      acc192(74)=abb192(79)
      acc192(75)=abb192(80)
      acc192(76)=abb192(81)
      acc192(77)=abb192(82)
      acc192(78)=abb192(83)
      acc192(79)=abb192(84)
      acc192(80)=abb192(85)
      acc192(81)=abb192(86)
      acc192(82)=abb192(87)
      acc192(83)=abb192(88)
      acc192(84)=abb192(89)
      acc192(85)=abb192(91)
      acc192(86)=Qspvak1k2*acc192(31)
      acc192(87)=Qspvak5k2*acc192(73)
      acc192(88)=Qspvae2e5*acc192(55)
      acc192(89)=-Qspvak2e5*acc192(81)
      acc192(89)=acc192(20)+acc192(89)
      acc192(89)=Qspvae2k2*acc192(89)
      acc192(90)=Qspvae2l3*acc192(22)
      acc192(91)=-Qspvae2l3*acc192(28)
      acc192(91)=acc192(23)+acc192(91)
      acc192(91)=Qspval3e5*acc192(91)
      acc192(92)=Qspvae2e5*acc192(28)
      acc192(92)=acc192(24)+acc192(92)
      acc192(92)=QspQ*acc192(92)
      acc192(86)=acc192(92)+acc192(91)+acc192(90)+acc192(89)+acc192(88)+acc192(&
      &87)+acc192(7)+acc192(86)
      acc192(86)=Qspvae5e1*acc192(86)
      acc192(87)=Qspvak1k2*acc192(26)
      acc192(88)=Qspvak5k2*acc192(83)
      acc192(89)=Qspvae2k1*acc192(11)
      acc192(90)=Qspvae2k5*acc192(30)
      acc192(91)=Qspvae5k2*acc192(56)
      acc192(92)=Qspvae2e5*acc192(51)
      acc192(93)=Qspvak2e5*acc192(82)
      acc192(93)=acc192(15)+acc192(93)
      acc192(93)=Qspvae2k2*acc192(93)
      acc192(94)=Qspvae2l3*acc192(61)
      acc192(95)=Qspvae2l3*acc192(57)
      acc192(95)=acc192(76)+acc192(95)
      acc192(95)=Qspval3e5*acc192(95)
      acc192(96)=Qspvae5l3*acc192(71)
      acc192(97)=-Qspvae2e5*acc192(57)
      acc192(97)=acc192(50)+acc192(97)
      acc192(97)=QspQ*acc192(97)
      acc192(86)=acc192(86)+acc192(97)+acc192(96)+acc192(95)+acc192(94)+acc192(&
      &93)+acc192(92)+acc192(91)+acc192(90)+acc192(89)+acc192(88)+acc192(2)+acc&
      &192(87)
      acc192(86)=Qspvae1e2*acc192(86)
      acc192(87)=Qspvak2k1*acc192(84)
      acc192(88)=Qspvak2k5*acc192(85)
      acc192(89)=Qspvae5e2*acc192(55)
      acc192(90)=-Qspvae5k2*acc192(81)
      acc192(90)=acc192(4)+acc192(90)
      acc192(90)=Qspvak2e2*acc192(90)
      acc192(91)=Qspval3e2*acc192(10)
      acc192(92)=-Qspval3e2*acc192(28)
      acc192(92)=acc192(27)+acc192(92)
      acc192(92)=Qspvae5l3*acc192(92)
      acc192(93)=Qspvae5e2*acc192(28)
      acc192(93)=acc192(35)+acc192(93)
      acc192(93)=QspQ*acc192(93)
      acc192(87)=acc192(93)+acc192(92)+acc192(91)+acc192(90)+acc192(89)+acc192(&
      &88)+acc192(6)+acc192(87)
      acc192(87)=Qspvae1e5*acc192(87)
      acc192(88)=Qspvak2k1*acc192(53)
      acc192(89)=Qspvak2k5*acc192(39)
      acc192(90)=Qspvak1e2*acc192(38)
      acc192(91)=Qspvak5e2*acc192(66)
      acc192(92)=Qspvak2e5*acc192(36)
      acc192(93)=Qspvae5e2*acc192(54)
      acc192(94)=Qspvae5k2*acc192(5)
      acc192(94)=acc192(12)+acc192(94)
      acc192(94)=Qspvak2e2*acc192(94)
      acc192(95)=Qspval3e2*acc192(16)
      acc192(96)=Qspval3e5*acc192(75)
      acc192(97)=Qspval3e2*acc192(58)
      acc192(97)=acc192(70)+acc192(97)
      acc192(97)=Qspvae5l3*acc192(97)
      acc192(98)=-Qspvae5e2*acc192(58)
      acc192(98)=acc192(34)+acc192(98)
      acc192(98)=QspQ*acc192(98)
      acc192(87)=acc192(87)+acc192(98)+acc192(97)+acc192(96)+acc192(95)+acc192(&
      &94)+acc192(93)+acc192(92)+acc192(91)+acc192(90)+acc192(89)+acc192(14)+ac&
      &c192(88)
      acc192(87)=Qspvae2e1*acc192(87)
      acc192(88)=Qspvae2k1*acc192(9)
      acc192(89)=Qspvae2k5*acc192(63)
      acc192(90)=Qspvae5k2*acc192(59)
      acc192(91)=Qspvae2k2*acc192(52)
      acc192(92)=Qspvae2l3*acc192(77)
      acc192(93)=Qspvae5l3*acc192(69)
      acc192(94)=QspQ*acc192(48)
      acc192(88)=acc192(94)+acc192(93)+acc192(92)+acc192(91)+acc192(90)+acc192(&
      &89)+acc192(8)+acc192(88)
      acc192(88)=Qspvae1e5*acc192(88)
      acc192(89)=Qspvak1e2*acc192(60)
      acc192(90)=Qspvak5e2*acc192(64)
      acc192(91)=-Qspvak2e5*acc192(78)
      acc192(92)=Qspvak2e2*acc192(29)
      acc192(93)=Qspval3e2*acc192(62)
      acc192(94)=Qspval3e5*acc192(74)
      acc192(95)=QspQ*acc192(18)
      acc192(89)=acc192(95)+acc192(94)+acc192(93)+acc192(92)+acc192(91)+acc192(&
      &90)+acc192(17)+acc192(89)
      acc192(89)=Qspvae5e1*acc192(89)
      acc192(90)=-Qspvae2e5*acc192(49)
      acc192(91)=-Qspvae5e2*acc192(46)
      acc192(90)=acc192(91)+acc192(33)+acc192(90)
      acc192(90)=QspQ*acc192(90)
      acc192(91)=Qspvak1k2*acc192(25)
      acc192(92)=Qspvak2k1*acc192(37)
      acc192(93)=Qspvak2k5*acc192(44)
      acc192(94)=Qspvak5k2*acc192(80)
      acc192(95)=Qspvak1e2*acc192(45)
      acc192(96)=Qspvae2k1*acc192(40)
      acc192(97)=Qspvak5e2*acc192(65)
      acc192(98)=Qspvae2k5*acc192(42)
      acc192(99)=Qspvak2e5*acc192(79)
      acc192(100)=Qspvae5k2*acc192(32)
      acc192(101)=Qspvae2e5*acc192(43)
      acc192(102)=Qspvae5e2*acc192(47)
      acc192(103)=Qspvae5k2*acc192(19)
      acc192(103)=acc192(21)+acc192(103)
      acc192(103)=Qspvak2e2*acc192(103)
      acc192(104)=Qspvak2e5*acc192(68)
      acc192(104)=acc192(41)+acc192(104)
      acc192(104)=Qspvae2k2*acc192(104)
      acc192(105)=Qspval3e2*acc192(13)
      acc192(106)=Qspvae2l3*acc192(1)
      acc192(107)=Qspvae2l3*acc192(49)
      acc192(107)=acc192(72)+acc192(107)
      acc192(107)=Qspval3e5*acc192(107)
      acc192(108)=Qspval3e2*acc192(46)
      acc192(108)=acc192(67)+acc192(108)
      acc192(108)=Qspvae5l3*acc192(108)
      brack=acc192(3)+acc192(86)+acc192(87)+acc192(88)+acc192(89)+acc192(90)+ac&
      &c192(91)+acc192(92)+acc192(93)+acc192(94)+acc192(95)+acc192(96)+acc192(9&
      &7)+acc192(98)+acc192(99)+acc192(100)+acc192(101)+acc192(102)+acc192(103)&
      &+acc192(104)+acc192(105)+acc192(106)+acc192(107)+acc192(108)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ subroutine numerator_ninja:
   subroutine numerator_ninja(ncut, Q_ext, mu2_ext, numerator) &
   & bind(c, name="p2_part21part21_part25part25part21_d192h0l1_qp_ninja")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki_nin
      use p2_part21part21_part25part25part21_globalsl1_qp, only: epspow
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_abbrevd192h0_qp
      implicit none
      integer(c_int), intent(in) :: ncut
      complex(ki_nin), dimension(0:3), intent(in) :: Q_ext
      complex(ki_nin), intent(in) :: mu2_ext
      complex(ki_nin), intent(out) :: numerator
      complex(ki) :: d192
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = k2
      Q(1:4)  =cmplx(real(+Q_ext(0:3)  -qshift(:),  ki_nin), aimag(+Q_ext(0:3))&
      &, ki)
      d192 = 0.0_ki
      d192 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d192, ki), aimag(d192), ki_nin)
   end subroutine numerator_ninja
   !------#] subroutine numerator_ninja:
!---#] numerator interfaces:
end module p2_part21part21_part25part25part21_d192h0l1_qp
