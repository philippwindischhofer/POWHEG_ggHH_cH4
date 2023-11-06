module     p2_part21part21_part25part25part21_d196h0l1_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d196h0l1_qp.f90
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
      use p2_part21part21_part25part25part21_abbrevd196h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc196(128)
      complex(ki) :: Qspvae1k1
      complex(ki) :: Qspvae5k1
      complex(ki) :: Qspvae1l4
      complex(ki) :: Qspval3e1
      complex(ki) :: Qspvak2e1
      complex(ki) :: Qspvae1k2
      complex(ki) :: Qspvae2e5
      complex(ki) :: Qspvak1e1
      complex(ki) :: Qspvak1e5
      complex(ki) :: Qspval4e1
      complex(ki) :: Qspvae1l3
      complex(ki) :: Qspvae5e2
      complex(ki) :: QspQ
      complex(ki) :: Qspvae1e5
      complex(ki) :: Qspvae5e1
      complex(ki) :: Qspvae1e2
      complex(ki) :: Qspvae2e1
      Qspvae1k1 = dotproduct(Q,spvae1k1)
      Qspvae5k1 = dotproduct(Q,spvae5k1)
      Qspvae1l4 = dotproduct(Q,spvae1l4)
      Qspval3e1 = dotproduct(Q,spval3e1)
      Qspvak2e1 = dotproduct(Q,spvak2e1)
      Qspvae1k2 = dotproduct(Q,spvae1k2)
      Qspvae2e5 = dotproduct(Q,spvae2e5)
      Qspvak1e1 = dotproduct(Q,spvak1e1)
      Qspvak1e5 = dotproduct(Q,spvak1e5)
      Qspval4e1 = dotproduct(Q,spval4e1)
      Qspvae1l3 = dotproduct(Q,spvae1l3)
      Qspvae5e2 = dotproduct(Q,spvae5e2)
      QspQ = dotproduct(Q,Q)
      Qspvae1e5 = dotproduct(Q,spvae1e5)
      Qspvae5e1 = dotproduct(Q,spvae5e1)
      Qspvae1e2 = dotproduct(Q,spvae1e2)
      Qspvae2e1 = dotproduct(Q,spvae2e1)
      acc196(1)=abb196(6)
      acc196(2)=abb196(7)
      acc196(3)=abb196(8)
      acc196(4)=abb196(9)
      acc196(5)=abb196(10)
      acc196(6)=abb196(11)
      acc196(7)=abb196(12)
      acc196(8)=abb196(13)
      acc196(9)=abb196(14)
      acc196(10)=abb196(15)
      acc196(11)=abb196(16)
      acc196(12)=abb196(17)
      acc196(13)=abb196(18)
      acc196(14)=abb196(19)
      acc196(15)=abb196(20)
      acc196(16)=abb196(21)
      acc196(17)=abb196(22)
      acc196(18)=abb196(23)
      acc196(19)=abb196(24)
      acc196(20)=abb196(25)
      acc196(21)=abb196(26)
      acc196(22)=abb196(27)
      acc196(23)=abb196(28)
      acc196(24)=abb196(29)
      acc196(25)=abb196(30)
      acc196(26)=abb196(31)
      acc196(27)=abb196(32)
      acc196(28)=abb196(33)
      acc196(29)=abb196(34)
      acc196(30)=abb196(35)
      acc196(31)=abb196(36)
      acc196(32)=abb196(37)
      acc196(33)=abb196(38)
      acc196(34)=abb196(39)
      acc196(35)=abb196(40)
      acc196(36)=abb196(41)
      acc196(37)=abb196(42)
      acc196(38)=abb196(43)
      acc196(39)=abb196(44)
      acc196(40)=abb196(45)
      acc196(41)=abb196(46)
      acc196(42)=abb196(47)
      acc196(43)=abb196(48)
      acc196(44)=abb196(49)
      acc196(45)=abb196(50)
      acc196(46)=abb196(51)
      acc196(47)=abb196(52)
      acc196(48)=abb196(53)
      acc196(49)=abb196(54)
      acc196(50)=abb196(55)
      acc196(51)=abb196(56)
      acc196(52)=abb196(57)
      acc196(53)=abb196(58)
      acc196(54)=abb196(59)
      acc196(55)=abb196(60)
      acc196(56)=abb196(61)
      acc196(57)=abb196(62)
      acc196(58)=abb196(63)
      acc196(59)=abb196(64)
      acc196(60)=abb196(65)
      acc196(61)=abb196(66)
      acc196(62)=abb196(67)
      acc196(63)=abb196(68)
      acc196(64)=abb196(69)
      acc196(65)=abb196(70)
      acc196(66)=abb196(71)
      acc196(67)=abb196(72)
      acc196(68)=abb196(73)
      acc196(69)=abb196(74)
      acc196(70)=abb196(75)
      acc196(71)=abb196(76)
      acc196(72)=abb196(77)
      acc196(73)=abb196(78)
      acc196(74)=abb196(79)
      acc196(75)=abb196(80)
      acc196(76)=abb196(81)
      acc196(77)=abb196(82)
      acc196(78)=abb196(83)
      acc196(79)=abb196(84)
      acc196(80)=abb196(85)
      acc196(81)=abb196(86)
      acc196(82)=abb196(87)
      acc196(83)=abb196(88)
      acc196(84)=abb196(89)
      acc196(85)=abb196(90)
      acc196(86)=abb196(91)
      acc196(87)=abb196(92)
      acc196(88)=abb196(93)
      acc196(89)=abb196(94)
      acc196(90)=abb196(95)
      acc196(91)=abb196(96)
      acc196(92)=abb196(97)
      acc196(93)=abb196(98)
      acc196(94)=abb196(99)
      acc196(95)=abb196(100)
      acc196(96)=abb196(101)
      acc196(97)=abb196(102)
      acc196(98)=abb196(103)
      acc196(99)=abb196(104)
      acc196(100)=abb196(105)
      acc196(101)=abb196(106)
      acc196(102)=abb196(109)
      acc196(103)=abb196(110)
      acc196(104)=abb196(111)
      acc196(105)=abb196(112)
      acc196(106)=abb196(113)
      acc196(107)=abb196(115)
      acc196(108)=abb196(116)
      acc196(109)=abb196(117)
      acc196(110)=abb196(118)
      acc196(111)=abb196(119)
      acc196(112)=Qspvae1k1*acc196(81)
      acc196(113)=Qspvae5k1*acc196(5)
      acc196(114)=Qspvae1l4*acc196(103)
      acc196(115)=Qspval3e1*acc196(67)
      acc196(116)=Qspvak2e1*acc196(27)
      acc196(117)=Qspvae1k2*acc196(34)
      acc196(112)=acc196(117)+acc196(116)+acc196(115)+acc196(114)+acc196(113)+a&
      &cc196(56)+acc196(112)
      acc196(112)=Qspvae2e5*acc196(112)
      acc196(113)=Qspvak1e1*acc196(83)
      acc196(114)=Qspvak1e5*acc196(21)
      acc196(115)=Qspval4e1*acc196(108)
      acc196(116)=Qspvae1l3*acc196(88)
      acc196(117)=Qspvak2e1*acc196(70)
      acc196(118)=Qspvae1k2*acc196(18)
      acc196(113)=acc196(118)+acc196(117)+acc196(116)+acc196(115)+acc196(114)+a&
      &cc196(65)+acc196(113)
      acc196(113)=Qspvae5e2*acc196(113)
      acc196(114)=Qspvae2e5*acc196(90)
      acc196(115)=Qspvae5e2*acc196(84)
      acc196(114)=acc196(115)+acc196(74)+acc196(114)
      acc196(114)=QspQ*acc196(114)
      acc196(115)=Qspvae1e5*acc196(93)
      acc196(116)=Qspvae5e1*acc196(91)
      acc196(117)=Qspvak1e1*acc196(61)
      acc196(118)=Qspvae1k1*acc196(75)
      acc196(119)=Qspvak1e5*acc196(62)
      acc196(120)=Qspvae5k1*acc196(53)
      acc196(121)=Qspval4e1*acc196(106)
      acc196(122)=Qspvae1l4*acc196(101)
      acc196(123)=Qspvae1e2*acc196(41)
      acc196(124)=Qspvae2e1*acc196(31)
      acc196(125)=Qspval3e1*acc196(37)
      acc196(126)=Qspvae1l3*acc196(47)
      acc196(127)=Qspvak2e1*acc196(12)
      acc196(128)=Qspvae1k2*acc196(13)
      acc196(112)=acc196(114)+acc196(113)+acc196(112)+acc196(128)+acc196(127)+a&
      &cc196(126)+acc196(125)+acc196(124)+acc196(123)+acc196(122)+acc196(121)+a&
      &cc196(120)+acc196(119)+acc196(118)+acc196(117)+acc196(116)+acc196(38)+ac&
      &c196(115)
      acc196(112)=QspQ*acc196(112)
      acc196(113)=Qspvae1k1*acc196(50)
      acc196(114)=Qspvae5k1*acc196(77)
      acc196(115)=Qspvae1l4*acc196(111)
      acc196(113)=acc196(115)+acc196(114)+acc196(22)+acc196(113)
      acc196(113)=Qspval3e1*acc196(113)
      acc196(114)=Qspvae1k1*acc196(51)
      acc196(115)=Qspvae5k1*acc196(79)
      acc196(116)=Qspvae1l4*acc196(44)
      acc196(114)=acc196(116)+acc196(115)+acc196(1)+acc196(114)
      acc196(114)=Qspvak2e1*acc196(114)
      acc196(115)=Qspval3e1*acc196(33)
      acc196(116)=-Qspvak2e1*acc196(16)
      acc196(115)=acc196(116)+acc196(2)+acc196(115)
      acc196(115)=Qspvae1k2*acc196(115)
      acc196(116)=Qspvae5e1*acc196(69)
      acc196(117)=Qspvae1k1*acc196(98)
      acc196(118)=Qspvae5k1*acc196(73)
      acc196(119)=Qspvae1l4*acc196(102)
      acc196(120)=-Qspvae5e1*acc196(94)
      acc196(120)=acc196(97)+acc196(120)
      acc196(120)=Qspvae1e2*acc196(120)
      acc196(113)=acc196(115)+acc196(114)+acc196(113)+acc196(120)+acc196(119)+a&
      &cc196(118)+acc196(117)+acc196(30)+acc196(116)
      acc196(113)=Qspvae2e5*acc196(113)
      acc196(114)=Qspvak1e1*acc196(71)
      acc196(115)=Qspvak1e5*acc196(96)
      acc196(116)=Qspval4e1*acc196(68)
      acc196(117)=Qspvak2e1*acc196(15)
      acc196(114)=acc196(117)+acc196(116)+acc196(115)+acc196(17)+acc196(114)
      acc196(114)=Qspvae1k2*acc196(114)
      acc196(115)=Qspvak1e1*acc196(60)
      acc196(116)=Qspvak1e5*acc196(86)
      acc196(117)=Qspval4e1*acc196(109)
      acc196(115)=acc196(117)+acc196(116)+acc196(36)+acc196(115)
      acc196(115)=Qspvae1l3*acc196(115)
      acc196(116)=Qspvae1e5*acc196(72)
      acc196(117)=Qspvak1e1*acc196(100)
      acc196(118)=Qspvak1e5*acc196(64)
      acc196(119)=Qspval4e1*acc196(107)
      acc196(120)=-Qspvae1e5*acc196(94)
      acc196(120)=acc196(92)+acc196(120)
      acc196(120)=Qspvae2e1*acc196(120)
      acc196(121)=Qspvae1l3*acc196(10)
      acc196(121)=acc196(35)+acc196(121)
      acc196(121)=Qspvak2e1*acc196(121)
      acc196(114)=acc196(114)+acc196(121)+acc196(115)+acc196(120)+acc196(119)+a&
      &cc196(118)+acc196(117)+acc196(80)+acc196(116)
      acc196(114)=Qspvae5e2*acc196(114)
      acc196(115)=Qspvae5e1*acc196(46)
      acc196(116)=Qspvak1e1*acc196(19)
      acc196(117)=Qspvak1e5*acc196(89)
      acc196(118)=Qspval4e1*acc196(45)
      acc196(119)=-Qspvae2e1*acc196(29)
      acc196(120)=-Qspval3e1*acc196(57)
      acc196(121)=Qspvak2e1*acc196(52)
      acc196(115)=acc196(121)+acc196(120)+acc196(119)+acc196(118)+acc196(117)+a&
      &cc196(116)+acc196(9)+acc196(115)
      acc196(115)=Qspvae1k2*acc196(115)
      acc196(116)=Qspvae1e5*acc196(58)
      acc196(117)=Qspvae1k1*acc196(3)
      acc196(118)=Qspvae5k1*acc196(78)
      acc196(119)=Qspvae1l4*acc196(43)
      acc196(120)=-Qspvae1e2*acc196(40)
      acc196(121)=-Qspvae1l3*acc196(66)
      acc196(116)=acc196(121)+acc196(120)+acc196(119)+acc196(118)+acc196(117)+a&
      &cc196(25)+acc196(116)
      acc196(116)=Qspvak2e1*acc196(116)
      acc196(117)=Qspvae5e1*acc196(99)
      acc196(118)=Qspvak1e1*acc196(85)
      acc196(119)=Qspvak1e5*acc196(26)
      acc196(120)=Qspval4e1*acc196(23)
      acc196(117)=acc196(120)+acc196(119)+acc196(118)+acc196(4)+acc196(117)
      acc196(117)=Qspvae1e2*acc196(117)
      acc196(118)=Qspvae1e5*acc196(95)
      acc196(119)=Qspvae1k1*acc196(87)
      acc196(120)=Qspvae5k1*acc196(14)
      acc196(121)=Qspvae1l4*acc196(7)
      acc196(118)=acc196(121)+acc196(120)+acc196(119)+acc196(20)+acc196(118)
      acc196(118)=Qspvae2e1*acc196(118)
      acc196(119)=Qspvae1e5*acc196(24)
      acc196(120)=Qspvae1k1*acc196(11)
      acc196(121)=Qspvae5k1*acc196(76)
      acc196(122)=Qspvae1l4*acc196(110)
      acc196(119)=acc196(122)+acc196(121)+acc196(120)+acc196(28)+acc196(119)
      acc196(119)=Qspval3e1*acc196(119)
      acc196(120)=Qspvae5e1*acc196(8)
      acc196(121)=Qspvak1e1*acc196(63)
      acc196(122)=Qspvak1e5*acc196(82)
      acc196(123)=Qspval4e1*acc196(104)
      acc196(120)=acc196(123)+acc196(122)+acc196(121)+acc196(39)+acc196(120)
      acc196(120)=Qspvae1l3*acc196(120)
      acc196(121)=Qspvae1e5*acc196(42)
      acc196(122)=Qspvae5e1*acc196(32)
      acc196(123)=Qspvak1e1*acc196(48)
      acc196(124)=Qspvae1k1*acc196(49)
      acc196(125)=Qspvak1e5*acc196(59)
      acc196(126)=Qspvae5k1*acc196(55)
      acc196(127)=Qspval4e1*acc196(105)
      acc196(128)=Qspvae1l4*acc196(54)
      brack=acc196(6)+acc196(112)+acc196(113)+acc196(114)+acc196(115)+acc196(11&
      &6)+acc196(117)+acc196(118)+acc196(119)+acc196(120)+acc196(121)+acc196(12&
      &2)+acc196(123)+acc196(124)+acc196(125)+acc196(126)+acc196(127)+acc196(12&
      &8)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ subroutine numerator_ninja:
   subroutine numerator_ninja(ncut, Q_ext, mu2_ext, numerator) &
   & bind(c, name="p2_part21part21_part25part25part21_d196h0l1_qp_ninja")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki_nin
      use p2_part21part21_part25part25part21_globalsl1_qp, only: epspow
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_abbrevd196h0_qp
      implicit none
      integer(c_int), intent(in) :: ncut
      complex(ki_nin), dimension(0:3), intent(in) :: Q_ext
      complex(ki_nin), intent(in) :: mu2_ext
      complex(ki_nin), intent(out) :: numerator
      complex(ki) :: d196
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = -k3
      Q(1:4)  =cmplx(real(-Q_ext(0:3)  -qshift(:),  ki_nin), aimag(-Q_ext(0:3))&
      &, ki)
      d196 = 0.0_ki
      d196 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d196, ki), aimag(d196), ki_nin)
   end subroutine numerator_ninja
   !------#] subroutine numerator_ninja:
!---#] numerator interfaces:
end module p2_part21part21_part25part25part21_d196h0l1_qp
