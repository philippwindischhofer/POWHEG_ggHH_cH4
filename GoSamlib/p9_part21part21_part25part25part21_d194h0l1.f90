module     p9_part21part21_part25part25part21_d194h0l1
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p9_part21part &
   ! &21_part25part25part21/helicity0d194h0l1.f90
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
      use p9_part21part21_part25part25part21_abbrevd194h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc194(108)
      complex(ki) :: Qspvak1k2
      complex(ki) :: Qspvak5k2
      complex(ki) :: Qspvae2e5
      complex(ki) :: Qspvak2e5
      complex(ki) :: Qspvae2k2
      complex(ki) :: Qspvae2l4
      complex(ki) :: Qspval4e5
      complex(ki) :: QspQ
      complex(ki) :: Qspvae5e1
      complex(ki) :: Qspvae2k1
      complex(ki) :: Qspvae2k5
      complex(ki) :: Qspvae5k2
      complex(ki) :: Qspvae5l4
      complex(ki) :: Qspvae1e2
      complex(ki) :: Qspvak2k1
      complex(ki) :: Qspvak2k5
      complex(ki) :: Qspvae5e2
      complex(ki) :: Qspvak2e2
      complex(ki) :: Qspval4e2
      complex(ki) :: Qspvae1e5
      complex(ki) :: Qspvak1e2
      complex(ki) :: Qspvak5e2
      complex(ki) :: Qspvae2e1
      Qspvak1k2 = dotproduct(Q,spvak1k2)
      Qspvak5k2 = dotproduct(Q,spvak5k2)
      Qspvae2e5 = dotproduct(Q,spvae2e5)
      Qspvak2e5 = dotproduct(Q,spvak2e5)
      Qspvae2k2 = dotproduct(Q,spvae2k2)
      Qspvae2l4 = dotproduct(Q,spvae2l4)
      Qspval4e5 = dotproduct(Q,spval4e5)
      QspQ = dotproduct(Q,Q)
      Qspvae5e1 = dotproduct(Q,spvae5e1)
      Qspvae2k1 = dotproduct(Q,spvae2k1)
      Qspvae2k5 = dotproduct(Q,spvae2k5)
      Qspvae5k2 = dotproduct(Q,spvae5k2)
      Qspvae5l4 = dotproduct(Q,spvae5l4)
      Qspvae1e2 = dotproduct(Q,spvae1e2)
      Qspvak2k1 = dotproduct(Q,spvak2k1)
      Qspvak2k5 = dotproduct(Q,spvak2k5)
      Qspvae5e2 = dotproduct(Q,spvae5e2)
      Qspvak2e2 = dotproduct(Q,spvak2e2)
      Qspval4e2 = dotproduct(Q,spval4e2)
      Qspvae1e5 = dotproduct(Q,spvae1e5)
      Qspvak1e2 = dotproduct(Q,spvak1e2)
      Qspvak5e2 = dotproduct(Q,spvak5e2)
      Qspvae2e1 = dotproduct(Q,spvae2e1)
      acc194(1)=abb194(4)
      acc194(2)=abb194(5)
      acc194(3)=abb194(6)
      acc194(4)=abb194(7)
      acc194(5)=abb194(8)
      acc194(6)=abb194(9)
      acc194(7)=abb194(10)
      acc194(8)=abb194(11)
      acc194(9)=abb194(12)
      acc194(10)=abb194(13)
      acc194(11)=abb194(14)
      acc194(12)=abb194(15)
      acc194(13)=abb194(16)
      acc194(14)=abb194(17)
      acc194(15)=abb194(18)
      acc194(16)=abb194(19)
      acc194(17)=abb194(20)
      acc194(18)=abb194(21)
      acc194(19)=abb194(22)
      acc194(20)=abb194(23)
      acc194(21)=abb194(24)
      acc194(22)=abb194(25)
      acc194(23)=abb194(26)
      acc194(24)=abb194(27)
      acc194(25)=abb194(28)
      acc194(26)=abb194(29)
      acc194(27)=abb194(30)
      acc194(28)=abb194(31)
      acc194(29)=abb194(32)
      acc194(30)=abb194(33)
      acc194(31)=abb194(34)
      acc194(32)=abb194(35)
      acc194(33)=abb194(36)
      acc194(34)=abb194(37)
      acc194(35)=abb194(38)
      acc194(36)=abb194(39)
      acc194(37)=abb194(40)
      acc194(38)=abb194(42)
      acc194(39)=abb194(43)
      acc194(40)=abb194(44)
      acc194(41)=abb194(45)
      acc194(42)=abb194(46)
      acc194(43)=abb194(47)
      acc194(44)=abb194(48)
      acc194(45)=abb194(49)
      acc194(46)=abb194(50)
      acc194(47)=abb194(51)
      acc194(48)=abb194(52)
      acc194(49)=abb194(53)
      acc194(50)=abb194(54)
      acc194(51)=abb194(55)
      acc194(52)=abb194(56)
      acc194(53)=abb194(57)
      acc194(54)=abb194(58)
      acc194(55)=abb194(59)
      acc194(56)=abb194(60)
      acc194(57)=abb194(62)
      acc194(58)=abb194(63)
      acc194(59)=abb194(64)
      acc194(60)=abb194(65)
      acc194(61)=abb194(66)
      acc194(62)=abb194(67)
      acc194(63)=abb194(68)
      acc194(64)=abb194(69)
      acc194(65)=abb194(70)
      acc194(66)=abb194(71)
      acc194(67)=abb194(72)
      acc194(68)=abb194(73)
      acc194(69)=abb194(74)
      acc194(70)=abb194(75)
      acc194(71)=abb194(76)
      acc194(72)=abb194(77)
      acc194(73)=abb194(78)
      acc194(74)=abb194(79)
      acc194(75)=abb194(80)
      acc194(76)=abb194(81)
      acc194(77)=abb194(82)
      acc194(78)=abb194(83)
      acc194(79)=abb194(84)
      acc194(80)=abb194(85)
      acc194(81)=abb194(86)
      acc194(82)=abb194(87)
      acc194(83)=abb194(88)
      acc194(84)=abb194(89)
      acc194(85)=abb194(91)
      acc194(86)=Qspvak1k2*acc194(31)
      acc194(87)=Qspvak5k2*acc194(73)
      acc194(88)=Qspvae2e5*acc194(55)
      acc194(89)=-Qspvak2e5*acc194(81)
      acc194(89)=acc194(20)+acc194(89)
      acc194(89)=Qspvae2k2*acc194(89)
      acc194(90)=Qspvae2l4*acc194(22)
      acc194(91)=-Qspvae2l4*acc194(28)
      acc194(91)=acc194(23)+acc194(91)
      acc194(91)=Qspval4e5*acc194(91)
      acc194(92)=Qspvae2e5*acc194(28)
      acc194(92)=acc194(24)+acc194(92)
      acc194(92)=QspQ*acc194(92)
      acc194(86)=acc194(92)+acc194(91)+acc194(90)+acc194(89)+acc194(88)+acc194(&
      &87)+acc194(7)+acc194(86)
      acc194(86)=Qspvae5e1*acc194(86)
      acc194(87)=Qspvak1k2*acc194(26)
      acc194(88)=Qspvak5k2*acc194(83)
      acc194(89)=Qspvae2k1*acc194(11)
      acc194(90)=Qspvae2k5*acc194(30)
      acc194(91)=Qspvae5k2*acc194(56)
      acc194(92)=Qspvae2e5*acc194(51)
      acc194(93)=Qspvak2e5*acc194(82)
      acc194(93)=acc194(15)+acc194(93)
      acc194(93)=Qspvae2k2*acc194(93)
      acc194(94)=Qspvae2l4*acc194(61)
      acc194(95)=Qspvae2l4*acc194(57)
      acc194(95)=acc194(76)+acc194(95)
      acc194(95)=Qspval4e5*acc194(95)
      acc194(96)=Qspvae5l4*acc194(71)
      acc194(97)=-Qspvae2e5*acc194(57)
      acc194(97)=acc194(50)+acc194(97)
      acc194(97)=QspQ*acc194(97)
      acc194(86)=acc194(86)+acc194(97)+acc194(96)+acc194(95)+acc194(94)+acc194(&
      &93)+acc194(92)+acc194(91)+acc194(90)+acc194(89)+acc194(88)+acc194(2)+acc&
      &194(87)
      acc194(86)=Qspvae1e2*acc194(86)
      acc194(87)=Qspvak2k1*acc194(84)
      acc194(88)=Qspvak2k5*acc194(85)
      acc194(89)=Qspvae5e2*acc194(55)
      acc194(90)=-Qspvae5k2*acc194(81)
      acc194(90)=acc194(4)+acc194(90)
      acc194(90)=Qspvak2e2*acc194(90)
      acc194(91)=Qspval4e2*acc194(10)
      acc194(92)=-Qspval4e2*acc194(28)
      acc194(92)=acc194(27)+acc194(92)
      acc194(92)=Qspvae5l4*acc194(92)
      acc194(93)=Qspvae5e2*acc194(28)
      acc194(93)=acc194(35)+acc194(93)
      acc194(93)=QspQ*acc194(93)
      acc194(87)=acc194(93)+acc194(92)+acc194(91)+acc194(90)+acc194(89)+acc194(&
      &88)+acc194(6)+acc194(87)
      acc194(87)=Qspvae1e5*acc194(87)
      acc194(88)=Qspvak2k1*acc194(53)
      acc194(89)=Qspvak2k5*acc194(39)
      acc194(90)=Qspvak1e2*acc194(38)
      acc194(91)=Qspvak5e2*acc194(66)
      acc194(92)=Qspvak2e5*acc194(36)
      acc194(93)=Qspvae5e2*acc194(54)
      acc194(94)=Qspvae5k2*acc194(5)
      acc194(94)=acc194(12)+acc194(94)
      acc194(94)=Qspvak2e2*acc194(94)
      acc194(95)=Qspval4e2*acc194(16)
      acc194(96)=Qspval4e5*acc194(75)
      acc194(97)=Qspval4e2*acc194(58)
      acc194(97)=acc194(70)+acc194(97)
      acc194(97)=Qspvae5l4*acc194(97)
      acc194(98)=-Qspvae5e2*acc194(58)
      acc194(98)=acc194(34)+acc194(98)
      acc194(98)=QspQ*acc194(98)
      acc194(87)=acc194(87)+acc194(98)+acc194(97)+acc194(96)+acc194(95)+acc194(&
      &94)+acc194(93)+acc194(92)+acc194(91)+acc194(90)+acc194(89)+acc194(14)+ac&
      &c194(88)
      acc194(87)=Qspvae2e1*acc194(87)
      acc194(88)=Qspvae2k1*acc194(9)
      acc194(89)=Qspvae2k5*acc194(63)
      acc194(90)=Qspvae5k2*acc194(59)
      acc194(91)=Qspvae2k2*acc194(52)
      acc194(92)=Qspvae2l4*acc194(77)
      acc194(93)=Qspvae5l4*acc194(69)
      acc194(94)=QspQ*acc194(48)
      acc194(88)=acc194(94)+acc194(93)+acc194(92)+acc194(91)+acc194(90)+acc194(&
      &89)+acc194(8)+acc194(88)
      acc194(88)=Qspvae1e5*acc194(88)
      acc194(89)=Qspvak1e2*acc194(60)
      acc194(90)=Qspvak5e2*acc194(64)
      acc194(91)=-Qspvak2e5*acc194(78)
      acc194(92)=Qspvak2e2*acc194(29)
      acc194(93)=Qspval4e2*acc194(62)
      acc194(94)=Qspval4e5*acc194(74)
      acc194(95)=QspQ*acc194(18)
      acc194(89)=acc194(95)+acc194(94)+acc194(93)+acc194(92)+acc194(91)+acc194(&
      &90)+acc194(17)+acc194(89)
      acc194(89)=Qspvae5e1*acc194(89)
      acc194(90)=-Qspvae2e5*acc194(49)
      acc194(91)=-Qspvae5e2*acc194(46)
      acc194(90)=acc194(91)+acc194(33)+acc194(90)
      acc194(90)=QspQ*acc194(90)
      acc194(91)=Qspvak1k2*acc194(25)
      acc194(92)=Qspvak2k1*acc194(37)
      acc194(93)=Qspvak2k5*acc194(44)
      acc194(94)=Qspvak5k2*acc194(80)
      acc194(95)=Qspvak1e2*acc194(45)
      acc194(96)=Qspvae2k1*acc194(40)
      acc194(97)=Qspvak5e2*acc194(65)
      acc194(98)=Qspvae2k5*acc194(42)
      acc194(99)=Qspvak2e5*acc194(79)
      acc194(100)=Qspvae5k2*acc194(32)
      acc194(101)=Qspvae2e5*acc194(43)
      acc194(102)=Qspvae5e2*acc194(47)
      acc194(103)=Qspvae5k2*acc194(19)
      acc194(103)=acc194(21)+acc194(103)
      acc194(103)=Qspvak2e2*acc194(103)
      acc194(104)=Qspvak2e5*acc194(68)
      acc194(104)=acc194(41)+acc194(104)
      acc194(104)=Qspvae2k2*acc194(104)
      acc194(105)=Qspval4e2*acc194(13)
      acc194(106)=Qspvae2l4*acc194(1)
      acc194(107)=Qspvae2l4*acc194(49)
      acc194(107)=acc194(72)+acc194(107)
      acc194(107)=Qspval4e5*acc194(107)
      acc194(108)=Qspval4e2*acc194(46)
      acc194(108)=acc194(67)+acc194(108)
      acc194(108)=Qspvae5l4*acc194(108)
      brack=acc194(3)+acc194(86)+acc194(87)+acc194(88)+acc194(89)+acc194(90)+ac&
      &c194(91)+acc194(92)+acc194(93)+acc194(94)+acc194(95)+acc194(96)+acc194(9&
      &7)+acc194(98)+acc194(99)+acc194(100)+acc194(101)+acc194(102)+acc194(103)&
      &+acc194(104)+acc194(105)+acc194(106)+acc194(107)+acc194(108)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ subroutine numerator_ninja:
   subroutine numerator_ninja(ncut, Q_ext, mu2_ext, numerator) &
   & bind(c, name="p9_part21part21_part25part25part21_d194h0l1_ninja")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki_nin
      use p9_part21part21_part25part25part21_globalsl1, only: epspow
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_abbrevd194h0
      implicit none
      integer(c_int), intent(in) :: ncut
      complex(ki_nin), dimension(0:3), intent(in) :: Q_ext
      complex(ki_nin), intent(in) :: mu2_ext
      complex(ki_nin), intent(out) :: numerator
      complex(ki) :: d194
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = k2
      Q(1:4)  =cmplx(real(+Q_ext(0:3)  -qshift(:),  ki_nin), aimag(+Q_ext(0:3))&
      &, ki)
      d194 = 0.0_ki
      d194 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d194, ki), aimag(d194), ki_nin)
   end subroutine numerator_ninja
   !------#] subroutine numerator_ninja:
!---#] numerator interfaces:
end module p9_part21part21_part25part25part21_d194h0l1
