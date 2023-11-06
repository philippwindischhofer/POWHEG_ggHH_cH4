module     p9_part21part21_part25part25part21_d200h0l1
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p9_part21part &
   ! &21_part25part25part21/helicity0d200h0l1.f90
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
      use p9_part21part21_part25part25part21_abbrevd200h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc200(128)
      complex(ki) :: Qspvae5k2
      complex(ki) :: Qspvae2k1
      complex(ki) :: Qspvae5l3
      complex(ki) :: Qspvae2k5
      complex(ki) :: Qspvak2e1
      complex(ki) :: Qspval3e1
      complex(ki) :: Qspvae1e5
      complex(ki) :: Qspvak2e5
      complex(ki) :: Qspvak1e2
      complex(ki) :: Qspval3e5
      complex(ki) :: Qspvak5e2
      complex(ki) :: Qspvae1k2
      complex(ki) :: Qspvae1l3
      complex(ki) :: Qspvae5e1
      complex(ki) :: QspQ
      complex(ki) :: Qspvae2e5
      complex(ki) :: Qspvae5e2
      complex(ki) :: Qspvae1e2
      complex(ki) :: Qspvae2e1
      Qspvae5k2 = dotproduct(Q,spvae5k2)
      Qspvae2k1 = dotproduct(Q,spvae2k1)
      Qspvae5l3 = dotproduct(Q,spvae5l3)
      Qspvae2k5 = dotproduct(Q,spvae2k5)
      Qspvak2e1 = dotproduct(Q,spvak2e1)
      Qspval3e1 = dotproduct(Q,spval3e1)
      Qspvae1e5 = dotproduct(Q,spvae1e5)
      Qspvak2e5 = dotproduct(Q,spvak2e5)
      Qspvak1e2 = dotproduct(Q,spvak1e2)
      Qspval3e5 = dotproduct(Q,spval3e5)
      Qspvak5e2 = dotproduct(Q,spvak5e2)
      Qspvae1k2 = dotproduct(Q,spvae1k2)
      Qspvae1l3 = dotproduct(Q,spvae1l3)
      Qspvae5e1 = dotproduct(Q,spvae5e1)
      QspQ = dotproduct(Q,Q)
      Qspvae2e5 = dotproduct(Q,spvae2e5)
      Qspvae5e2 = dotproduct(Q,spvae5e2)
      Qspvae1e2 = dotproduct(Q,spvae1e2)
      Qspvae2e1 = dotproduct(Q,spvae2e1)
      acc200(1)=abb200(4)
      acc200(2)=abb200(5)
      acc200(3)=abb200(6)
      acc200(4)=abb200(7)
      acc200(5)=abb200(8)
      acc200(6)=abb200(9)
      acc200(7)=abb200(10)
      acc200(8)=abb200(11)
      acc200(9)=abb200(12)
      acc200(10)=abb200(13)
      acc200(11)=abb200(14)
      acc200(12)=abb200(15)
      acc200(13)=abb200(16)
      acc200(14)=abb200(17)
      acc200(15)=abb200(18)
      acc200(16)=abb200(19)
      acc200(17)=abb200(20)
      acc200(18)=abb200(21)
      acc200(19)=abb200(22)
      acc200(20)=abb200(23)
      acc200(21)=abb200(24)
      acc200(22)=abb200(25)
      acc200(23)=abb200(26)
      acc200(24)=abb200(27)
      acc200(25)=abb200(28)
      acc200(26)=abb200(29)
      acc200(27)=abb200(30)
      acc200(28)=abb200(31)
      acc200(29)=abb200(32)
      acc200(30)=abb200(33)
      acc200(31)=abb200(34)
      acc200(32)=abb200(35)
      acc200(33)=abb200(36)
      acc200(34)=abb200(37)
      acc200(35)=abb200(38)
      acc200(36)=abb200(39)
      acc200(37)=abb200(40)
      acc200(38)=abb200(41)
      acc200(39)=abb200(42)
      acc200(40)=abb200(43)
      acc200(41)=abb200(44)
      acc200(42)=abb200(45)
      acc200(43)=abb200(46)
      acc200(44)=abb200(47)
      acc200(45)=abb200(48)
      acc200(46)=abb200(49)
      acc200(47)=abb200(50)
      acc200(48)=abb200(51)
      acc200(49)=abb200(52)
      acc200(50)=abb200(53)
      acc200(51)=abb200(54)
      acc200(52)=abb200(55)
      acc200(53)=abb200(56)
      acc200(54)=abb200(57)
      acc200(55)=abb200(58)
      acc200(56)=abb200(59)
      acc200(57)=abb200(60)
      acc200(58)=abb200(61)
      acc200(59)=abb200(62)
      acc200(60)=abb200(63)
      acc200(61)=abb200(64)
      acc200(62)=abb200(65)
      acc200(63)=abb200(66)
      acc200(64)=abb200(67)
      acc200(65)=abb200(68)
      acc200(66)=abb200(69)
      acc200(67)=abb200(70)
      acc200(68)=abb200(71)
      acc200(69)=abb200(72)
      acc200(70)=abb200(73)
      acc200(71)=abb200(74)
      acc200(72)=abb200(75)
      acc200(73)=abb200(76)
      acc200(74)=abb200(77)
      acc200(75)=abb200(78)
      acc200(76)=abb200(79)
      acc200(77)=abb200(80)
      acc200(78)=abb200(81)
      acc200(79)=abb200(82)
      acc200(80)=abb200(83)
      acc200(81)=abb200(84)
      acc200(82)=abb200(85)
      acc200(83)=abb200(86)
      acc200(84)=abb200(87)
      acc200(85)=abb200(88)
      acc200(86)=abb200(89)
      acc200(87)=abb200(90)
      acc200(88)=abb200(91)
      acc200(89)=abb200(92)
      acc200(90)=abb200(93)
      acc200(91)=abb200(94)
      acc200(92)=abb200(95)
      acc200(93)=abb200(96)
      acc200(94)=abb200(97)
      acc200(95)=abb200(98)
      acc200(96)=abb200(99)
      acc200(97)=abb200(100)
      acc200(98)=abb200(101)
      acc200(99)=abb200(102)
      acc200(100)=abb200(103)
      acc200(101)=abb200(104)
      acc200(102)=abb200(105)
      acc200(103)=abb200(106)
      acc200(104)=abb200(107)
      acc200(105)=abb200(108)
      acc200(106)=abb200(109)
      acc200(107)=abb200(110)
      acc200(108)=abb200(111)
      acc200(109)=abb200(112)
      acc200(110)=-Qspvae5k2*acc200(97)
      acc200(111)=Qspvae2k1*acc200(71)
      acc200(112)=-Qspvae5l3*acc200(56)
      acc200(113)=Qspvae2k5*acc200(106)
      acc200(114)=Qspvak2e1*acc200(21)
      acc200(115)=Qspval3e1*acc200(52)
      acc200(110)=acc200(115)+acc200(114)+acc200(113)+acc200(112)+acc200(111)+a&
      &cc200(72)+acc200(110)
      acc200(110)=Qspvae1e5*acc200(110)
      acc200(111)=-Qspvak2e5*acc200(81)
      acc200(112)=Qspvak1e2*acc200(55)
      acc200(113)=-Qspval3e5*acc200(88)
      acc200(114)=Qspvak5e2*acc200(108)
      acc200(115)=Qspvae1k2*acc200(67)
      acc200(116)=Qspvae1l3*acc200(44)
      acc200(111)=acc200(116)+acc200(115)+acc200(114)+acc200(113)+acc200(112)+a&
      &cc200(51)+acc200(111)
      acc200(111)=Qspvae5e1*acc200(111)
      acc200(112)=Qspvae1e5*acc200(101)
      acc200(113)=Qspvae5e1*acc200(95)
      acc200(112)=acc200(113)+acc200(68)+acc200(112)
      acc200(112)=QspQ*acc200(112)
      acc200(113)=Qspvak2e5*acc200(90)
      acc200(114)=-Qspvae5k2*acc200(70)
      acc200(115)=Qspvae2e5*acc200(69)
      acc200(116)=Qspvae5e2*acc200(74)
      acc200(117)=Qspvak1e2*acc200(39)
      acc200(118)=Qspvae2k1*acc200(62)
      acc200(119)=-Qspval3e5*acc200(87)
      acc200(120)=-Qspvae5l3*acc200(93)
      acc200(121)=Qspvak5e2*acc200(107)
      acc200(122)=Qspvae2k5*acc200(94)
      acc200(123)=-Qspvae1e2*acc200(77)
      acc200(124)=-Qspvae2e1*acc200(75)
      acc200(125)=Qspvak2e1*acc200(91)
      acc200(126)=Qspvae1k2*acc200(79)
      acc200(127)=Qspval3e1*acc200(50)
      acc200(128)=Qspvae1l3*acc200(42)
      acc200(110)=acc200(112)+acc200(111)+acc200(110)+acc200(128)+acc200(127)+a&
      &cc200(126)+acc200(125)+acc200(124)+acc200(123)+acc200(122)+acc200(121)+a&
      &cc200(120)+acc200(119)+acc200(118)+acc200(117)+acc200(116)+acc200(115)+a&
      &cc200(114)+acc200(45)+acc200(113)
      acc200(110)=QspQ*acc200(110)
      acc200(111)=-Qspvae5k2*acc200(14)
      acc200(112)=Qspvae2k1*acc200(41)
      acc200(113)=-Qspvae5l3*acc200(84)
      acc200(114)=Qspvae2k5*acc200(60)
      acc200(111)=acc200(114)+acc200(113)+acc200(112)+acc200(49)+acc200(111)
      acc200(111)=Qspval3e1*acc200(111)
      acc200(112)=Qspvae2k1*acc200(76)
      acc200(113)=Qspvae5l3*acc200(103)
      acc200(114)=Qspvae2k5*acc200(98)
      acc200(112)=acc200(114)+acc200(113)+acc200(1)+acc200(112)
      acc200(112)=Qspvak2e1*acc200(112)
      acc200(113)=Qspvae5k2*acc200(46)
      acc200(114)=Qspvae5e2*acc200(82)
      acc200(115)=Qspvae2k1*acc200(24)
      acc200(116)=Qspvae5l3*acc200(11)
      acc200(117)=Qspvae2k5*acc200(3)
      acc200(118)=-Qspvae5e2*acc200(59)
      acc200(118)=acc200(33)+acc200(118)
      acc200(118)=Qspvae2e1*acc200(118)
      acc200(111)=acc200(111)+acc200(112)+acc200(118)+acc200(117)+acc200(116)+a&
      &cc200(115)+acc200(114)+acc200(35)+acc200(113)
      acc200(111)=Qspvae1e5*acc200(111)
      acc200(112)=-Qspvak2e5*acc200(32)
      acc200(113)=Qspvak1e2*acc200(37)
      acc200(114)=-Qspval3e5*acc200(84)
      acc200(115)=Qspvak5e2*acc200(83)
      acc200(112)=acc200(115)+acc200(114)+acc200(113)+acc200(43)+acc200(112)
      acc200(112)=Qspvae1l3*acc200(112)
      acc200(113)=Qspvak1e2*acc200(99)
      acc200(114)=Qspval3e5*acc200(96)
      acc200(115)=Qspvak5e2*acc200(85)
      acc200(113)=acc200(115)+acc200(114)+acc200(26)+acc200(113)
      acc200(113)=Qspvae1k2*acc200(113)
      acc200(114)=Qspvak2e5*acc200(38)
      acc200(115)=Qspvae2e5*acc200(100)
      acc200(116)=Qspvak1e2*acc200(30)
      acc200(117)=Qspval3e5*acc200(23)
      acc200(118)=Qspvak5e2*acc200(12)
      acc200(119)=-Qspvae2e5*acc200(59)
      acc200(119)=-acc200(64)+acc200(119)
      acc200(119)=Qspvae1e2*acc200(119)
      acc200(112)=acc200(112)+acc200(113)+acc200(119)+acc200(118)+acc200(117)+a&
      &cc200(116)+acc200(115)+acc200(16)+acc200(114)
      acc200(112)=Qspvae5e1*acc200(112)
      acc200(113)=Qspvae5k2*acc200(105)
      acc200(114)=Qspvae2e5*acc200(54)
      acc200(115)=Qspvae2k1*acc200(29)
      acc200(116)=Qspvae5l3*acc200(31)
      acc200(117)=Qspvae2k5*acc200(9)
      acc200(113)=acc200(117)+acc200(116)+acc200(115)+acc200(114)+acc200(13)+ac&
      &c200(113)
      acc200(113)=Qspvae1e2*acc200(113)
      acc200(114)=Qspvak2e5*acc200(104)
      acc200(115)=Qspvae5e2*acc200(102)
      acc200(116)=Qspvak1e2*acc200(34)
      acc200(117)=Qspval3e5*acc200(36)
      acc200(118)=Qspvak5e2*acc200(22)
      acc200(114)=acc200(118)+acc200(117)+acc200(116)+acc200(115)+acc200(25)+ac&
      &c200(114)
      acc200(114)=Qspvae2e1*acc200(114)
      acc200(115)=Qspvae5k2*acc200(89)
      acc200(116)=-Qspvae2e5*acc200(36)
      acc200(117)=Qspvae2k1*acc200(92)
      acc200(118)=Qspvae5l3*acc200(66)
      acc200(119)=Qspvae2k5*acc200(53)
      acc200(115)=acc200(119)+acc200(118)+acc200(117)+acc200(116)+acc200(10)+ac&
      &c200(115)
      acc200(115)=Qspval3e1*acc200(115)
      acc200(116)=Qspvak2e5*acc200(17)
      acc200(117)=-Qspvae5e2*acc200(31)
      acc200(118)=Qspvak1e2*acc200(86)
      acc200(119)=Qspval3e5*acc200(109)
      acc200(120)=Qspvak5e2*acc200(65)
      acc200(116)=acc200(120)+acc200(119)+acc200(118)+acc200(117)+acc200(28)+ac&
      &c200(116)
      acc200(116)=Qspvae1l3*acc200(116)
      acc200(117)=Qspvae2e5*acc200(15)
      acc200(118)=Qspvae2k1*acc200(4)
      acc200(119)=Qspvae5l3*acc200(5)
      acc200(120)=Qspvae2k5*acc200(57)
      acc200(117)=acc200(120)+acc200(119)+acc200(118)+acc200(2)+acc200(117)
      acc200(117)=Qspvak2e1*acc200(117)
      acc200(118)=Qspvae5e2*acc200(8)
      acc200(119)=Qspvak1e2*acc200(80)
      acc200(120)=Qspval3e5*acc200(18)
      acc200(121)=Qspvak5e2*acc200(78)
      acc200(118)=acc200(121)+acc200(120)+acc200(119)+acc200(20)+acc200(118)
      acc200(118)=Qspvae1k2*acc200(118)
      acc200(119)=Qspvak2e5*acc200(40)
      acc200(120)=Qspvae5k2*acc200(48)
      acc200(121)=Qspvae2e5*acc200(63)
      acc200(122)=Qspvae5e2*acc200(73)
      acc200(123)=Qspvak1e2*acc200(27)
      acc200(124)=Qspvae2k1*acc200(47)
      acc200(125)=Qspval3e5*acc200(19)
      acc200(126)=Qspvae5l3*acc200(6)
      acc200(127)=Qspvak5e2*acc200(58)
      acc200(128)=Qspvae2k5*acc200(61)
      brack=acc200(7)+acc200(110)+acc200(111)+acc200(112)+acc200(113)+acc200(11&
      &4)+acc200(115)+acc200(116)+acc200(117)+acc200(118)+acc200(119)+acc200(12&
      &0)+acc200(121)+acc200(122)+acc200(123)+acc200(124)+acc200(125)+acc200(12&
      &6)+acc200(127)+acc200(128)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ subroutine numerator_ninja:
   subroutine numerator_ninja(ncut, Q_ext, mu2_ext, numerator) &
   & bind(c, name="p9_part21part21_part25part25part21_d200h0l1_ninja")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki_nin
      use p9_part21part21_part25part25part21_globalsl1, only: epspow
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_abbrevd200h0
      implicit none
      integer(c_int), intent(in) :: ncut
      complex(ki_nin), dimension(0:3), intent(in) :: Q_ext
      complex(ki_nin), intent(in) :: mu2_ext
      complex(ki_nin), intent(out) :: numerator
      complex(ki) :: d200
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = -k3
      Q(1:4)  =cmplx(real(-Q_ext(0:3)  -qshift(:),  ki_nin), aimag(-Q_ext(0:3))&
      &, ki)
      d200 = 0.0_ki
      d200 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d200, ki), aimag(d200), ki_nin)
   end subroutine numerator_ninja
   !------#] subroutine numerator_ninja:
!---#] numerator interfaces:
end module p9_part21part21_part25part25part21_d200h0l1
