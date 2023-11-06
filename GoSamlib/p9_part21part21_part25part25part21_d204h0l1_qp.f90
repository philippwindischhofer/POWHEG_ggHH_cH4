module     p9_part21part21_part25part25part21_d204h0l1_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p9_part21part &
   ! &21_part25part25part21/helicity0d204h0l1_qp.f90
   ! generator: buildfortran.py
   use p9_part21part21_part25part25part21_config, only: ki => ki_qp
   use p9_part21part21_part25part25part21_util_qp, only: cond
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   public :: numerator_ninja
contains
!---#[ function brack_1:
   pure function brack_1(Q,mu2) result(brack)
      use p9_part21part21_part25part25part21_model_qp
      use p9_part21part21_part25part25part21_kinematics_qp
      use p9_part21part21_part25part25part21_color_qp
      use p9_part21part21_part25part25part21_abbrevd204h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc204(128)
      complex(ki) :: Qspvae5k2
      complex(ki) :: Qspvae2k1
      complex(ki) :: Qspvae5l4
      complex(ki) :: Qspvae2k5
      complex(ki) :: Qspvak2e1
      complex(ki) :: Qspval4e1
      complex(ki) :: Qspvae1e5
      complex(ki) :: Qspvak2e5
      complex(ki) :: Qspvak1e2
      complex(ki) :: Qspval4e5
      complex(ki) :: Qspvak5e2
      complex(ki) :: Qspvae1k2
      complex(ki) :: Qspvae1l4
      complex(ki) :: Qspvae5e1
      complex(ki) :: QspQ
      complex(ki) :: Qspvae2e5
      complex(ki) :: Qspvae5e2
      complex(ki) :: Qspvae1e2
      complex(ki) :: Qspvae2e1
      Qspvae5k2 = dotproduct(Q,spvae5k2)
      Qspvae2k1 = dotproduct(Q,spvae2k1)
      Qspvae5l4 = dotproduct(Q,spvae5l4)
      Qspvae2k5 = dotproduct(Q,spvae2k5)
      Qspvak2e1 = dotproduct(Q,spvak2e1)
      Qspval4e1 = dotproduct(Q,spval4e1)
      Qspvae1e5 = dotproduct(Q,spvae1e5)
      Qspvak2e5 = dotproduct(Q,spvak2e5)
      Qspvak1e2 = dotproduct(Q,spvak1e2)
      Qspval4e5 = dotproduct(Q,spval4e5)
      Qspvak5e2 = dotproduct(Q,spvak5e2)
      Qspvae1k2 = dotproduct(Q,spvae1k2)
      Qspvae1l4 = dotproduct(Q,spvae1l4)
      Qspvae5e1 = dotproduct(Q,spvae5e1)
      QspQ = dotproduct(Q,Q)
      Qspvae2e5 = dotproduct(Q,spvae2e5)
      Qspvae5e2 = dotproduct(Q,spvae5e2)
      Qspvae1e2 = dotproduct(Q,spvae1e2)
      Qspvae2e1 = dotproduct(Q,spvae2e1)
      acc204(1)=abb204(4)
      acc204(2)=abb204(5)
      acc204(3)=abb204(6)
      acc204(4)=abb204(7)
      acc204(5)=abb204(8)
      acc204(6)=abb204(9)
      acc204(7)=abb204(10)
      acc204(8)=abb204(11)
      acc204(9)=abb204(12)
      acc204(10)=abb204(13)
      acc204(11)=abb204(14)
      acc204(12)=abb204(15)
      acc204(13)=abb204(16)
      acc204(14)=abb204(17)
      acc204(15)=abb204(18)
      acc204(16)=abb204(19)
      acc204(17)=abb204(20)
      acc204(18)=abb204(21)
      acc204(19)=abb204(22)
      acc204(20)=abb204(23)
      acc204(21)=abb204(24)
      acc204(22)=abb204(25)
      acc204(23)=abb204(26)
      acc204(24)=abb204(27)
      acc204(25)=abb204(28)
      acc204(26)=abb204(29)
      acc204(27)=abb204(30)
      acc204(28)=abb204(31)
      acc204(29)=abb204(32)
      acc204(30)=abb204(33)
      acc204(31)=abb204(34)
      acc204(32)=abb204(35)
      acc204(33)=abb204(36)
      acc204(34)=abb204(37)
      acc204(35)=abb204(38)
      acc204(36)=abb204(39)
      acc204(37)=abb204(40)
      acc204(38)=abb204(41)
      acc204(39)=abb204(42)
      acc204(40)=abb204(43)
      acc204(41)=abb204(44)
      acc204(42)=abb204(45)
      acc204(43)=abb204(46)
      acc204(44)=abb204(47)
      acc204(45)=abb204(48)
      acc204(46)=abb204(49)
      acc204(47)=abb204(50)
      acc204(48)=abb204(51)
      acc204(49)=abb204(52)
      acc204(50)=abb204(53)
      acc204(51)=abb204(54)
      acc204(52)=abb204(55)
      acc204(53)=abb204(56)
      acc204(54)=abb204(57)
      acc204(55)=abb204(58)
      acc204(56)=abb204(59)
      acc204(57)=abb204(60)
      acc204(58)=abb204(61)
      acc204(59)=abb204(62)
      acc204(60)=abb204(63)
      acc204(61)=abb204(64)
      acc204(62)=abb204(65)
      acc204(63)=abb204(66)
      acc204(64)=abb204(67)
      acc204(65)=abb204(68)
      acc204(66)=abb204(69)
      acc204(67)=abb204(70)
      acc204(68)=abb204(71)
      acc204(69)=abb204(72)
      acc204(70)=abb204(73)
      acc204(71)=abb204(74)
      acc204(72)=abb204(75)
      acc204(73)=abb204(76)
      acc204(74)=abb204(77)
      acc204(75)=abb204(78)
      acc204(76)=abb204(79)
      acc204(77)=abb204(80)
      acc204(78)=abb204(81)
      acc204(79)=abb204(82)
      acc204(80)=abb204(83)
      acc204(81)=abb204(84)
      acc204(82)=abb204(85)
      acc204(83)=abb204(86)
      acc204(84)=abb204(87)
      acc204(85)=abb204(88)
      acc204(86)=abb204(89)
      acc204(87)=abb204(90)
      acc204(88)=abb204(91)
      acc204(89)=abb204(92)
      acc204(90)=abb204(93)
      acc204(91)=abb204(94)
      acc204(92)=abb204(95)
      acc204(93)=abb204(96)
      acc204(94)=abb204(97)
      acc204(95)=abb204(98)
      acc204(96)=abb204(99)
      acc204(97)=abb204(100)
      acc204(98)=abb204(101)
      acc204(99)=abb204(102)
      acc204(100)=abb204(103)
      acc204(101)=abb204(104)
      acc204(102)=abb204(105)
      acc204(103)=abb204(106)
      acc204(104)=abb204(107)
      acc204(105)=abb204(108)
      acc204(106)=abb204(109)
      acc204(107)=abb204(110)
      acc204(108)=abb204(111)
      acc204(109)=abb204(112)
      acc204(110)=-Qspvae5k2*acc204(97)
      acc204(111)=Qspvae2k1*acc204(71)
      acc204(112)=-Qspvae5l4*acc204(56)
      acc204(113)=Qspvae2k5*acc204(106)
      acc204(114)=Qspvak2e1*acc204(21)
      acc204(115)=Qspval4e1*acc204(52)
      acc204(110)=acc204(115)+acc204(114)+acc204(113)+acc204(112)+acc204(111)+a&
      &cc204(72)+acc204(110)
      acc204(110)=Qspvae1e5*acc204(110)
      acc204(111)=-Qspvak2e5*acc204(81)
      acc204(112)=Qspvak1e2*acc204(55)
      acc204(113)=-Qspval4e5*acc204(88)
      acc204(114)=Qspvak5e2*acc204(108)
      acc204(115)=Qspvae1k2*acc204(67)
      acc204(116)=Qspvae1l4*acc204(44)
      acc204(111)=acc204(116)+acc204(115)+acc204(114)+acc204(113)+acc204(112)+a&
      &cc204(51)+acc204(111)
      acc204(111)=Qspvae5e1*acc204(111)
      acc204(112)=Qspvae1e5*acc204(101)
      acc204(113)=Qspvae5e1*acc204(95)
      acc204(112)=acc204(113)+acc204(68)+acc204(112)
      acc204(112)=QspQ*acc204(112)
      acc204(113)=Qspvak2e5*acc204(90)
      acc204(114)=-Qspvae5k2*acc204(70)
      acc204(115)=Qspvae2e5*acc204(69)
      acc204(116)=Qspvae5e2*acc204(74)
      acc204(117)=Qspvak1e2*acc204(39)
      acc204(118)=Qspvae2k1*acc204(62)
      acc204(119)=-Qspval4e5*acc204(87)
      acc204(120)=-Qspvae5l4*acc204(93)
      acc204(121)=Qspvak5e2*acc204(107)
      acc204(122)=Qspvae2k5*acc204(94)
      acc204(123)=-Qspvae1e2*acc204(77)
      acc204(124)=-Qspvae2e1*acc204(75)
      acc204(125)=Qspvak2e1*acc204(91)
      acc204(126)=Qspvae1k2*acc204(79)
      acc204(127)=Qspval4e1*acc204(50)
      acc204(128)=Qspvae1l4*acc204(42)
      acc204(110)=acc204(112)+acc204(111)+acc204(110)+acc204(128)+acc204(127)+a&
      &cc204(126)+acc204(125)+acc204(124)+acc204(123)+acc204(122)+acc204(121)+a&
      &cc204(120)+acc204(119)+acc204(118)+acc204(117)+acc204(116)+acc204(115)+a&
      &cc204(114)+acc204(45)+acc204(113)
      acc204(110)=QspQ*acc204(110)
      acc204(111)=-Qspvae5k2*acc204(14)
      acc204(112)=Qspvae2k1*acc204(41)
      acc204(113)=-Qspvae5l4*acc204(84)
      acc204(114)=Qspvae2k5*acc204(60)
      acc204(111)=acc204(114)+acc204(113)+acc204(112)+acc204(49)+acc204(111)
      acc204(111)=Qspval4e1*acc204(111)
      acc204(112)=Qspvae2k1*acc204(76)
      acc204(113)=Qspvae5l4*acc204(103)
      acc204(114)=Qspvae2k5*acc204(98)
      acc204(112)=acc204(114)+acc204(113)+acc204(1)+acc204(112)
      acc204(112)=Qspvak2e1*acc204(112)
      acc204(113)=Qspvae5k2*acc204(46)
      acc204(114)=Qspvae5e2*acc204(82)
      acc204(115)=Qspvae2k1*acc204(24)
      acc204(116)=Qspvae5l4*acc204(11)
      acc204(117)=Qspvae2k5*acc204(3)
      acc204(118)=-Qspvae5e2*acc204(59)
      acc204(118)=acc204(33)+acc204(118)
      acc204(118)=Qspvae2e1*acc204(118)
      acc204(111)=acc204(111)+acc204(112)+acc204(118)+acc204(117)+acc204(116)+a&
      &cc204(115)+acc204(114)+acc204(35)+acc204(113)
      acc204(111)=Qspvae1e5*acc204(111)
      acc204(112)=-Qspvak2e5*acc204(32)
      acc204(113)=Qspvak1e2*acc204(37)
      acc204(114)=-Qspval4e5*acc204(84)
      acc204(115)=Qspvak5e2*acc204(83)
      acc204(112)=acc204(115)+acc204(114)+acc204(113)+acc204(43)+acc204(112)
      acc204(112)=Qspvae1l4*acc204(112)
      acc204(113)=Qspvak1e2*acc204(99)
      acc204(114)=Qspval4e5*acc204(96)
      acc204(115)=Qspvak5e2*acc204(85)
      acc204(113)=acc204(115)+acc204(114)+acc204(26)+acc204(113)
      acc204(113)=Qspvae1k2*acc204(113)
      acc204(114)=Qspvak2e5*acc204(38)
      acc204(115)=Qspvae2e5*acc204(100)
      acc204(116)=Qspvak1e2*acc204(30)
      acc204(117)=Qspval4e5*acc204(23)
      acc204(118)=Qspvak5e2*acc204(12)
      acc204(119)=-Qspvae2e5*acc204(59)
      acc204(119)=-acc204(64)+acc204(119)
      acc204(119)=Qspvae1e2*acc204(119)
      acc204(112)=acc204(112)+acc204(113)+acc204(119)+acc204(118)+acc204(117)+a&
      &cc204(116)+acc204(115)+acc204(16)+acc204(114)
      acc204(112)=Qspvae5e1*acc204(112)
      acc204(113)=Qspvae5k2*acc204(105)
      acc204(114)=Qspvae2e5*acc204(54)
      acc204(115)=Qspvae2k1*acc204(29)
      acc204(116)=Qspvae5l4*acc204(31)
      acc204(117)=Qspvae2k5*acc204(9)
      acc204(113)=acc204(117)+acc204(116)+acc204(115)+acc204(114)+acc204(13)+ac&
      &c204(113)
      acc204(113)=Qspvae1e2*acc204(113)
      acc204(114)=Qspvak2e5*acc204(104)
      acc204(115)=Qspvae5e2*acc204(102)
      acc204(116)=Qspvak1e2*acc204(34)
      acc204(117)=Qspval4e5*acc204(36)
      acc204(118)=Qspvak5e2*acc204(22)
      acc204(114)=acc204(118)+acc204(117)+acc204(116)+acc204(115)+acc204(25)+ac&
      &c204(114)
      acc204(114)=Qspvae2e1*acc204(114)
      acc204(115)=Qspvae5k2*acc204(89)
      acc204(116)=-Qspvae2e5*acc204(36)
      acc204(117)=Qspvae2k1*acc204(92)
      acc204(118)=Qspvae5l4*acc204(66)
      acc204(119)=Qspvae2k5*acc204(53)
      acc204(115)=acc204(119)+acc204(118)+acc204(117)+acc204(116)+acc204(10)+ac&
      &c204(115)
      acc204(115)=Qspval4e1*acc204(115)
      acc204(116)=Qspvak2e5*acc204(17)
      acc204(117)=-Qspvae5e2*acc204(31)
      acc204(118)=Qspvak1e2*acc204(86)
      acc204(119)=Qspval4e5*acc204(109)
      acc204(120)=Qspvak5e2*acc204(65)
      acc204(116)=acc204(120)+acc204(119)+acc204(118)+acc204(117)+acc204(28)+ac&
      &c204(116)
      acc204(116)=Qspvae1l4*acc204(116)
      acc204(117)=Qspvae2e5*acc204(15)
      acc204(118)=Qspvae2k1*acc204(4)
      acc204(119)=Qspvae5l4*acc204(5)
      acc204(120)=Qspvae2k5*acc204(57)
      acc204(117)=acc204(120)+acc204(119)+acc204(118)+acc204(2)+acc204(117)
      acc204(117)=Qspvak2e1*acc204(117)
      acc204(118)=Qspvae5e2*acc204(8)
      acc204(119)=Qspvak1e2*acc204(80)
      acc204(120)=Qspval4e5*acc204(18)
      acc204(121)=Qspvak5e2*acc204(78)
      acc204(118)=acc204(121)+acc204(120)+acc204(119)+acc204(20)+acc204(118)
      acc204(118)=Qspvae1k2*acc204(118)
      acc204(119)=Qspvak2e5*acc204(40)
      acc204(120)=Qspvae5k2*acc204(48)
      acc204(121)=Qspvae2e5*acc204(63)
      acc204(122)=Qspvae5e2*acc204(73)
      acc204(123)=Qspvak1e2*acc204(27)
      acc204(124)=Qspvae2k1*acc204(47)
      acc204(125)=Qspval4e5*acc204(19)
      acc204(126)=Qspvae5l4*acc204(6)
      acc204(127)=Qspvak5e2*acc204(58)
      acc204(128)=Qspvae2k5*acc204(61)
      brack=acc204(7)+acc204(110)+acc204(111)+acc204(112)+acc204(113)+acc204(11&
      &4)+acc204(115)+acc204(116)+acc204(117)+acc204(118)+acc204(119)+acc204(12&
      &0)+acc204(121)+acc204(122)+acc204(123)+acc204(124)+acc204(125)+acc204(12&
      &6)+acc204(127)+acc204(128)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ subroutine numerator_ninja:
   subroutine numerator_ninja(ncut, Q_ext, mu2_ext, numerator) &
   & bind(c, name="p9_part21part21_part25part25part21_d204h0l1_qp_ninja")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki_nin
      use p9_part21part21_part25part25part21_globalsl1_qp, only: epspow
      use p9_part21part21_part25part25part21_kinematics_qp
      use p9_part21part21_part25part25part21_abbrevd204h0_qp
      implicit none
      integer(c_int), intent(in) :: ncut
      complex(ki_nin), dimension(0:3), intent(in) :: Q_ext
      complex(ki_nin), intent(in) :: mu2_ext
      complex(ki_nin), intent(out) :: numerator
      complex(ki) :: d204
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = k2-k4
      Q(1:4)  =cmplx(real(+Q_ext(0:3)  -qshift(:),  ki_nin), aimag(+Q_ext(0:3))&
      &, ki)
      d204 = 0.0_ki
      d204 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d204, ki), aimag(d204), ki_nin)
   end subroutine numerator_ninja
   !------#] subroutine numerator_ninja:
!---#] numerator interfaces:
end module p9_part21part21_part25part25part21_d204h0l1_qp
