module     p9_part21part21_part25part25part21_d198h0l1_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p9_part21part &
   ! &21_part25part25part21/helicity0d198h0l1_qp.f90
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
      use p9_part21part21_part25part25part21_abbrevd198h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc198(128)
      complex(ki) :: Qspvak1e1
      complex(ki) :: Qspvak1e2
      complex(ki) :: Qspval4e1
      complex(ki) :: Qspvae1l3
      complex(ki) :: Qspvak2e1
      complex(ki) :: Qspvae1k2
      complex(ki) :: Qspvae2e5
      complex(ki) :: Qspvae1k1
      complex(ki) :: Qspvae2k1
      complex(ki) :: Qspvae1l4
      complex(ki) :: Qspval3e1
      complex(ki) :: Qspvae5e2
      complex(ki) :: QspQ
      complex(ki) :: Qspvae1e2
      complex(ki) :: Qspvae2e1
      complex(ki) :: Qspvae1e5
      complex(ki) :: Qspvae5e1
      Qspvak1e1 = dotproduct(Q,spvak1e1)
      Qspvak1e2 = dotproduct(Q,spvak1e2)
      Qspval4e1 = dotproduct(Q,spval4e1)
      Qspvae1l3 = dotproduct(Q,spvae1l3)
      Qspvak2e1 = dotproduct(Q,spvak2e1)
      Qspvae1k2 = dotproduct(Q,spvae1k2)
      Qspvae2e5 = dotproduct(Q,spvae2e5)
      Qspvae1k1 = dotproduct(Q,spvae1k1)
      Qspvae2k1 = dotproduct(Q,spvae2k1)
      Qspvae1l4 = dotproduct(Q,spvae1l4)
      Qspval3e1 = dotproduct(Q,spval3e1)
      Qspvae5e2 = dotproduct(Q,spvae5e2)
      QspQ = dotproduct(Q,Q)
      Qspvae1e2 = dotproduct(Q,spvae1e2)
      Qspvae2e1 = dotproduct(Q,spvae2e1)
      Qspvae1e5 = dotproduct(Q,spvae1e5)
      Qspvae5e1 = dotproduct(Q,spvae5e1)
      acc198(1)=abb198(6)
      acc198(2)=abb198(7)
      acc198(3)=abb198(8)
      acc198(4)=abb198(9)
      acc198(5)=abb198(10)
      acc198(6)=abb198(11)
      acc198(7)=abb198(12)
      acc198(8)=abb198(13)
      acc198(9)=abb198(14)
      acc198(10)=abb198(15)
      acc198(11)=abb198(16)
      acc198(12)=abb198(17)
      acc198(13)=abb198(18)
      acc198(14)=abb198(19)
      acc198(15)=abb198(20)
      acc198(16)=abb198(21)
      acc198(17)=abb198(22)
      acc198(18)=abb198(23)
      acc198(19)=abb198(24)
      acc198(20)=abb198(25)
      acc198(21)=abb198(26)
      acc198(22)=abb198(27)
      acc198(23)=abb198(28)
      acc198(24)=abb198(29)
      acc198(25)=abb198(30)
      acc198(26)=abb198(31)
      acc198(27)=abb198(32)
      acc198(28)=abb198(33)
      acc198(29)=abb198(34)
      acc198(30)=abb198(35)
      acc198(31)=abb198(36)
      acc198(32)=abb198(37)
      acc198(33)=abb198(38)
      acc198(34)=abb198(39)
      acc198(35)=abb198(40)
      acc198(36)=abb198(41)
      acc198(37)=abb198(42)
      acc198(38)=abb198(43)
      acc198(39)=abb198(44)
      acc198(40)=abb198(45)
      acc198(41)=abb198(46)
      acc198(42)=abb198(47)
      acc198(43)=abb198(48)
      acc198(44)=abb198(49)
      acc198(45)=abb198(50)
      acc198(46)=abb198(51)
      acc198(47)=abb198(52)
      acc198(48)=abb198(53)
      acc198(49)=abb198(54)
      acc198(50)=abb198(55)
      acc198(51)=abb198(56)
      acc198(52)=abb198(57)
      acc198(53)=abb198(58)
      acc198(54)=abb198(59)
      acc198(55)=abb198(60)
      acc198(56)=abb198(61)
      acc198(57)=abb198(62)
      acc198(58)=abb198(63)
      acc198(59)=abb198(64)
      acc198(60)=abb198(65)
      acc198(61)=abb198(66)
      acc198(62)=abb198(67)
      acc198(63)=abb198(68)
      acc198(64)=abb198(69)
      acc198(65)=abb198(70)
      acc198(66)=abb198(71)
      acc198(67)=abb198(72)
      acc198(68)=abb198(73)
      acc198(69)=abb198(74)
      acc198(70)=abb198(75)
      acc198(71)=abb198(76)
      acc198(72)=abb198(77)
      acc198(73)=abb198(78)
      acc198(74)=abb198(79)
      acc198(75)=abb198(80)
      acc198(76)=abb198(81)
      acc198(77)=abb198(82)
      acc198(78)=abb198(83)
      acc198(79)=abb198(84)
      acc198(80)=abb198(85)
      acc198(81)=abb198(86)
      acc198(82)=abb198(87)
      acc198(83)=abb198(88)
      acc198(84)=abb198(89)
      acc198(85)=abb198(90)
      acc198(86)=abb198(91)
      acc198(87)=abb198(92)
      acc198(88)=abb198(93)
      acc198(89)=abb198(94)
      acc198(90)=abb198(96)
      acc198(91)=abb198(97)
      acc198(92)=abb198(100)
      acc198(93)=abb198(101)
      acc198(94)=abb198(103)
      acc198(95)=abb198(106)
      acc198(96)=abb198(107)
      acc198(97)=abb198(108)
      acc198(98)=abb198(109)
      acc198(99)=abb198(110)
      acc198(100)=abb198(111)
      acc198(101)=abb198(112)
      acc198(102)=abb198(113)
      acc198(103)=abb198(114)
      acc198(104)=abb198(116)
      acc198(105)=abb198(119)
      acc198(106)=abb198(120)
      acc198(107)=abb198(121)
      acc198(108)=abb198(123)
      acc198(109)=abb198(126)
      acc198(110)=abb198(127)
      acc198(111)=abb198(128)
      acc198(112)=Qspvak1e1*acc198(63)
      acc198(113)=Qspvak1e2*acc198(36)
      acc198(114)=Qspval4e1*acc198(103)
      acc198(115)=Qspvae1l3*acc198(104)
      acc198(116)=Qspvak2e1*acc198(28)
      acc198(117)=Qspvae1k2*acc198(23)
      acc198(112)=acc198(117)+acc198(116)+acc198(115)+acc198(114)+acc198(113)+a&
      &cc198(65)+acc198(112)
      acc198(112)=Qspvae2e5*acc198(112)
      acc198(113)=Qspvae1k1*acc198(86)
      acc198(114)=Qspvae2k1*acc198(47)
      acc198(115)=Qspvae1l4*acc198(102)
      acc198(116)=Qspval3e1*acc198(108)
      acc198(117)=Qspvak2e1*acc198(20)
      acc198(118)=Qspvae1k2*acc198(38)
      acc198(113)=acc198(118)+acc198(117)+acc198(116)+acc198(115)+acc198(114)+a&
      &cc198(85)+acc198(113)
      acc198(113)=Qspvae5e2*acc198(113)
      acc198(114)=Qspvae2e5*acc198(94)
      acc198(115)=Qspvae5e2*acc198(93)
      acc198(114)=acc198(115)+acc198(34)+acc198(114)
      acc198(114)=QspQ*acc198(114)
      acc198(115)=Qspvae1e2*acc198(98)
      acc198(116)=Qspvae2e1*acc198(90)
      acc198(117)=Qspvak1e1*acc198(30)
      acc198(118)=Qspvae1k1*acc198(50)
      acc198(119)=Qspvak1e2*acc198(2)
      acc198(120)=Qspvae2k1*acc198(13)
      acc198(121)=Qspval4e1*acc198(31)
      acc198(122)=Qspvae1l4*acc198(40)
      acc198(123)=Qspvae1e5*acc198(42)
      acc198(124)=Qspvae5e1*acc198(32)
      acc198(125)=Qspval3e1*acc198(95)
      acc198(126)=Qspvae1l3*acc198(88)
      acc198(127)=Qspvak2e1*acc198(18)
      acc198(128)=Qspvae1k2*acc198(16)
      acc198(112)=acc198(114)+acc198(113)+acc198(112)+acc198(128)+acc198(127)+a&
      &cc198(126)+acc198(125)+acc198(124)+acc198(123)+acc198(122)+acc198(121)+a&
      &cc198(120)+acc198(119)+acc198(118)+acc198(117)+acc198(116)+acc198(19)+ac&
      &c198(115)
      acc198(112)=QspQ*acc198(112)
      acc198(113)=Qspvak1e1*acc198(52)
      acc198(114)=Qspvak1e2*acc198(67)
      acc198(115)=Qspval4e1*acc198(79)
      acc198(116)=Qspvak2e1*acc198(9)
      acc198(113)=acc198(116)+acc198(115)+acc198(114)+acc198(15)+acc198(113)
      acc198(113)=Qspvae1k2*acc198(113)
      acc198(114)=Qspvak1e1*acc198(72)
      acc198(115)=Qspvak1e2*acc198(84)
      acc198(116)=Qspval4e1*acc198(107)
      acc198(114)=acc198(116)+acc198(115)+acc198(33)+acc198(114)
      acc198(114)=Qspvae1l3*acc198(114)
      acc198(115)=Qspvae1e2*acc198(53)
      acc198(116)=Qspvak1e1*acc198(7)
      acc198(117)=Qspvak1e2*acc198(14)
      acc198(118)=Qspval4e1*acc198(60)
      acc198(119)=Qspvae1e2*acc198(97)
      acc198(119)=acc198(78)+acc198(119)
      acc198(119)=Qspvae5e1*acc198(119)
      acc198(120)=Qspvae1l3*acc198(51)
      acc198(120)=acc198(26)+acc198(120)
      acc198(120)=Qspvak2e1*acc198(120)
      acc198(113)=acc198(113)+acc198(120)+acc198(114)+acc198(119)+acc198(118)+a&
      &cc198(117)+acc198(116)+acc198(5)+acc198(115)
      acc198(113)=Qspvae2e5*acc198(113)
      acc198(114)=Qspvae1k1*acc198(87)
      acc198(115)=Qspvae2k1*acc198(82)
      acc198(116)=Qspvae1l4*acc198(111)
      acc198(114)=acc198(116)+acc198(115)+acc198(43)+acc198(114)
      acc198(114)=Qspval3e1*acc198(114)
      acc198(115)=Qspvae1k1*acc198(4)
      acc198(116)=Qspvae2k1*acc198(73)
      acc198(117)=Qspvae1l4*acc198(74)
      acc198(115)=acc198(117)+acc198(116)+acc198(25)+acc198(115)
      acc198(115)=Qspvak2e1*acc198(115)
      acc198(116)=Qspval3e1*acc198(70)
      acc198(117)=Qspvak2e1*acc198(10)
      acc198(116)=acc198(117)+acc198(27)+acc198(116)
      acc198(116)=Qspvae1k2*acc198(116)
      acc198(117)=Qspvae2e1*acc198(64)
      acc198(118)=Qspvae1k1*acc198(21)
      acc198(119)=Qspvae2k1*acc198(24)
      acc198(120)=Qspvae1l4*acc198(49)
      acc198(121)=Qspvae2e1*acc198(97)
      acc198(121)=acc198(81)+acc198(121)
      acc198(121)=Qspvae1e5*acc198(121)
      acc198(114)=acc198(116)+acc198(115)+acc198(114)+acc198(121)+acc198(120)+a&
      &cc198(119)+acc198(118)+acc198(62)+acc198(117)
      acc198(114)=Qspvae5e2*acc198(114)
      acc198(115)=Qspvae2e1*acc198(41)
      acc198(116)=Qspvak1e1*acc198(89)
      acc198(117)=Qspvak1e2*acc198(54)
      acc198(118)=Qspval4e1*acc198(66)
      acc198(119)=-Qspvae5e1*acc198(77)
      acc198(120)=Qspval3e1*acc198(80)
      acc198(121)=Qspvak2e1*acc198(12)
      acc198(115)=acc198(121)+acc198(120)+acc198(119)+acc198(118)+acc198(117)+a&
      &cc198(116)+acc198(8)+acc198(115)
      acc198(115)=Qspvae1k2*acc198(115)
      acc198(116)=Qspvae1e2*acc198(44)
      acc198(117)=Qspvae1k1*acc198(83)
      acc198(118)=Qspvae2k1*acc198(55)
      acc198(119)=Qspvae1l4*acc198(56)
      acc198(120)=-Qspvae1e5*acc198(75)
      acc198(121)=Qspvae1l3*acc198(76)
      acc198(116)=acc198(121)+acc198(120)+acc198(119)+acc198(118)+acc198(117)+a&
      &cc198(17)+acc198(116)
      acc198(116)=Qspvak2e1*acc198(116)
      acc198(117)=Qspvae2e1*acc198(96)
      acc198(118)=Qspvak1e1*acc198(91)
      acc198(119)=Qspvak1e2*acc198(58)
      acc198(120)=Qspval4e1*acc198(100)
      acc198(117)=acc198(120)+acc198(119)+acc198(118)+acc198(39)+acc198(117)
      acc198(117)=Qspvae1e5*acc198(117)
      acc198(118)=Qspvae1e2*acc198(101)
      acc198(119)=Qspvae1k1*acc198(92)
      acc198(120)=Qspvae2k1*acc198(48)
      acc198(121)=Qspvae1l4*acc198(99)
      acc198(118)=acc198(121)+acc198(120)+acc198(119)+acc198(29)+acc198(118)
      acc198(118)=Qspvae5e1*acc198(118)
      acc198(119)=-Qspvae1e2*acc198(110)
      acc198(120)=Qspvae1k1*acc198(45)
      acc198(121)=Qspvae2k1*acc198(61)
      acc198(122)=Qspvae1l4*acc198(109)
      acc198(119)=acc198(122)+acc198(121)+acc198(120)+acc198(46)+acc198(119)
      acc198(119)=Qspval3e1*acc198(119)
      acc198(120)=-Qspvae2e1*acc198(106)
      acc198(121)=Qspvak1e1*acc198(35)
      acc198(122)=Qspvak1e2*acc198(59)
      acc198(123)=Qspval4e1*acc198(105)
      acc198(120)=acc198(123)+acc198(122)+acc198(121)+acc198(37)+acc198(120)
      acc198(120)=Qspvae1l3*acc198(120)
      acc198(121)=Qspvae1e2*acc198(71)
      acc198(122)=Qspvae2e1*acc198(69)
      acc198(123)=Qspvak1e1*acc198(11)
      acc198(124)=Qspvae1k1*acc198(22)
      acc198(125)=Qspvak1e2*acc198(1)
      acc198(126)=Qspvae2k1*acc198(3)
      acc198(127)=Qspval4e1*acc198(68)
      acc198(128)=Qspvae1l4*acc198(57)
      brack=acc198(6)+acc198(112)+acc198(113)+acc198(114)+acc198(115)+acc198(11&
      &6)+acc198(117)+acc198(118)+acc198(119)+acc198(120)+acc198(121)+acc198(12&
      &2)+acc198(123)+acc198(124)+acc198(125)+acc198(126)+acc198(127)+acc198(12&
      &8)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ subroutine numerator_ninja:
   subroutine numerator_ninja(ncut, Q_ext, mu2_ext, numerator) &
   & bind(c, name="p9_part21part21_part25part25part21_d198h0l1_qp_ninja")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki_nin
      use p9_part21part21_part25part25part21_globalsl1_qp, only: epspow
      use p9_part21part21_part25part25part21_kinematics_qp
      use p9_part21part21_part25part25part21_abbrevd198h0_qp
      implicit none
      integer(c_int), intent(in) :: ncut
      complex(ki_nin), dimension(0:3), intent(in) :: Q_ext
      complex(ki_nin), intent(in) :: mu2_ext
      complex(ki_nin), intent(out) :: numerator
      complex(ki) :: d198
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = -k3+k2-k5
      Q(1:4)  =cmplx(real(+Q_ext(0:3)  -qshift(:),  ki_nin), aimag(+Q_ext(0:3))&
      &, ki)
      d198 = 0.0_ki
      d198 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d198, ki), aimag(d198), ki_nin)
   end subroutine numerator_ninja
   !------#] subroutine numerator_ninja:
!---#] numerator interfaces:
end module p9_part21part21_part25part25part21_d198h0l1_qp
