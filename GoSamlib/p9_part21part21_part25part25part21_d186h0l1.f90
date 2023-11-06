module     p9_part21part21_part25part25part21_d186h0l1
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p9_part21part &
   ! &21_part25part25part21/helicity0d186h0l1.f90
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
      use p9_part21part21_part25part25part21_abbrevd186h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc186(120)
      complex(ki) :: Qspvak1k2
      complex(ki) :: Qspvak1k5
      complex(ki) :: Qspvae5e1
      complex(ki) :: Qspvae5k1
      complex(ki) :: Qspvak1e1
      complex(ki) :: Qspvak2e1
      complex(ki) :: Qspval3e1
      complex(ki) :: Qspvae5k2
      complex(ki) :: Qspvae5l3
      complex(ki) :: QspQ
      complex(ki) :: Qspvae2e5
      complex(ki) :: Qspvak5e1
      complex(ki) :: Qspvak1e5
      complex(ki) :: Qspvak2e5
      complex(ki) :: Qspval3e5
      complex(ki) :: Qspvae1e2
      complex(ki) :: Qspvak2k1
      complex(ki) :: Qspvak5k1
      complex(ki) :: Qspvae1e5
      complex(ki) :: Qspvae1k1
      complex(ki) :: Qspvae1k2
      complex(ki) :: Qspvae1l3
      complex(ki) :: Qspvae5e2
      complex(ki) :: Qspvae1k5
      complex(ki) :: Qspvae2e1
      Qspvak1k2 = dotproduct(Q,spvak1k2)
      Qspvak1k5 = dotproduct(Q,spvak1k5)
      Qspvae5e1 = dotproduct(Q,spvae5e1)
      Qspvae5k1 = dotproduct(Q,spvae5k1)
      Qspvak1e1 = dotproduct(Q,spvak1e1)
      Qspvak2e1 = dotproduct(Q,spvak2e1)
      Qspval3e1 = dotproduct(Q,spval3e1)
      Qspvae5k2 = dotproduct(Q,spvae5k2)
      Qspvae5l3 = dotproduct(Q,spvae5l3)
      QspQ = dotproduct(Q,Q)
      Qspvae2e5 = dotproduct(Q,spvae2e5)
      Qspvak5e1 = dotproduct(Q,spvak5e1)
      Qspvak1e5 = dotproduct(Q,spvak1e5)
      Qspvak2e5 = dotproduct(Q,spvak2e5)
      Qspval3e5 = dotproduct(Q,spval3e5)
      Qspvae1e2 = dotproduct(Q,spvae1e2)
      Qspvak2k1 = dotproduct(Q,spvak2k1)
      Qspvak5k1 = dotproduct(Q,spvak5k1)
      Qspvae1e5 = dotproduct(Q,spvae1e5)
      Qspvae1k1 = dotproduct(Q,spvae1k1)
      Qspvae1k2 = dotproduct(Q,spvae1k2)
      Qspvae1l3 = dotproduct(Q,spvae1l3)
      Qspvae5e2 = dotproduct(Q,spvae5e2)
      Qspvae1k5 = dotproduct(Q,spvae1k5)
      Qspvae2e1 = dotproduct(Q,spvae2e1)
      acc186(1)=abb186(4)
      acc186(2)=abb186(5)
      acc186(3)=abb186(6)
      acc186(4)=abb186(7)
      acc186(5)=abb186(8)
      acc186(6)=abb186(9)
      acc186(7)=abb186(10)
      acc186(8)=abb186(11)
      acc186(9)=abb186(12)
      acc186(10)=abb186(13)
      acc186(11)=abb186(14)
      acc186(12)=abb186(15)
      acc186(13)=abb186(16)
      acc186(14)=abb186(17)
      acc186(15)=abb186(18)
      acc186(16)=abb186(19)
      acc186(17)=abb186(20)
      acc186(18)=abb186(21)
      acc186(19)=abb186(22)
      acc186(20)=abb186(23)
      acc186(21)=abb186(24)
      acc186(22)=abb186(25)
      acc186(23)=abb186(26)
      acc186(24)=abb186(27)
      acc186(25)=abb186(28)
      acc186(26)=abb186(29)
      acc186(27)=abb186(30)
      acc186(28)=abb186(31)
      acc186(29)=abb186(32)
      acc186(30)=abb186(33)
      acc186(31)=abb186(34)
      acc186(32)=abb186(35)
      acc186(33)=abb186(36)
      acc186(34)=abb186(37)
      acc186(35)=abb186(38)
      acc186(36)=abb186(39)
      acc186(37)=abb186(40)
      acc186(38)=abb186(41)
      acc186(39)=abb186(42)
      acc186(40)=abb186(43)
      acc186(41)=abb186(44)
      acc186(42)=abb186(45)
      acc186(43)=abb186(46)
      acc186(44)=abb186(47)
      acc186(45)=abb186(48)
      acc186(46)=abb186(49)
      acc186(47)=abb186(50)
      acc186(48)=abb186(51)
      acc186(49)=abb186(52)
      acc186(50)=abb186(53)
      acc186(51)=abb186(54)
      acc186(52)=abb186(55)
      acc186(53)=abb186(56)
      acc186(54)=abb186(57)
      acc186(55)=abb186(58)
      acc186(56)=abb186(59)
      acc186(57)=abb186(60)
      acc186(58)=abb186(61)
      acc186(59)=abb186(62)
      acc186(60)=abb186(63)
      acc186(61)=abb186(64)
      acc186(62)=abb186(65)
      acc186(63)=abb186(66)
      acc186(64)=abb186(67)
      acc186(65)=abb186(68)
      acc186(66)=abb186(69)
      acc186(67)=abb186(70)
      acc186(68)=abb186(71)
      acc186(69)=abb186(72)
      acc186(70)=abb186(73)
      acc186(71)=abb186(74)
      acc186(72)=abb186(75)
      acc186(73)=abb186(76)
      acc186(74)=abb186(80)
      acc186(75)=abb186(81)
      acc186(76)=abb186(82)
      acc186(77)=abb186(83)
      acc186(78)=abb186(84)
      acc186(79)=abb186(85)
      acc186(80)=abb186(87)
      acc186(81)=abb186(89)
      acc186(82)=abb186(90)
      acc186(83)=abb186(91)
      acc186(84)=abb186(92)
      acc186(85)=abb186(94)
      acc186(86)=abb186(95)
      acc186(87)=abb186(97)
      acc186(88)=abb186(98)
      acc186(89)=abb186(99)
      acc186(90)=abb186(100)
      acc186(91)=abb186(101)
      acc186(92)=abb186(102)
      acc186(93)=abb186(103)
      acc186(94)=abb186(104)
      acc186(95)=abb186(105)
      acc186(96)=Qspvak1k2*acc186(40)
      acc186(97)=-Qspvak1k5*acc186(35)
      acc186(98)=-Qspvae5e1*acc186(76)
      acc186(99)=-acc186(82)*Qspvae5k1
      acc186(99)=acc186(61)+acc186(99)
      acc186(99)=Qspvak1e1*acc186(99)
      acc186(100)=Qspvak2e1*acc186(69)
      acc186(101)=Qspval3e1*acc186(9)
      acc186(102)=-Qspvak2e1*acc186(78)
      acc186(102)=-acc186(55)+acc186(102)
      acc186(102)=Qspvae5k2*acc186(102)
      acc186(103)=Qspval3e1*acc186(82)
      acc186(103)=-acc186(16)+acc186(103)
      acc186(103)=Qspvae5l3*acc186(103)
      acc186(104)=acc186(82)*Qspvae5e1
      acc186(104)=acc186(2)+acc186(104)
      acc186(104)=QspQ*acc186(104)
      acc186(96)=acc186(104)+acc186(103)+acc186(102)+acc186(101)+acc186(100)+ac&
      &c186(99)+acc186(98)+acc186(97)+acc186(7)+acc186(96)
      acc186(96)=Qspvae2e5*acc186(96)
      acc186(97)=Qspvak1k5*acc186(19)
      acc186(98)=Qspvak5e1*acc186(27)
      acc186(99)=Qspvak1e5*acc186(95)
      acc186(100)=Qspvae5e1*acc186(34)
      acc186(101)=-Qspvae5k1*acc186(36)
      acc186(101)=acc186(51)+acc186(101)
      acc186(101)=Qspvak1e1*acc186(101)
      acc186(102)=Qspvak2e1*acc186(54)
      acc186(103)=Qspval3e1*acc186(85)
      acc186(104)=Qspvak2e5*acc186(57)
      acc186(105)=Qspvak2e1*acc186(80)
      acc186(105)=acc186(93)+acc186(105)
      acc186(105)=Qspvae5k2*acc186(105)
      acc186(106)=Qspval3e5*acc186(88)
      acc186(107)=Qspval3e1*acc186(36)
      acc186(107)=acc186(62)+acc186(107)
      acc186(107)=Qspvae5l3*acc186(107)
      acc186(108)=Qspvae5e1*acc186(36)
      acc186(108)=acc186(33)+acc186(108)
      acc186(108)=QspQ*acc186(108)
      acc186(96)=acc186(96)+acc186(108)+acc186(107)+acc186(106)+acc186(105)+acc&
      &186(104)+acc186(103)+acc186(102)+acc186(101)+acc186(100)+acc186(99)+acc1&
      &86(98)+acc186(12)+acc186(97)
      acc186(96)=Qspvae1e2*acc186(96)
      acc186(97)=Qspvak2k1*acc186(48)
      acc186(98)=-Qspvak5k1*acc186(74)
      acc186(99)=-Qspvae1e5*acc186(76)
      acc186(100)=-acc186(82)*Qspvak1e5
      acc186(100)=acc186(77)+acc186(100)
      acc186(100)=Qspvae1k1*acc186(100)
      acc186(101)=Qspvae1k2*acc186(59)
      acc186(102)=Qspvae1l3*acc186(46)
      acc186(103)=-Qspvae1k2*acc186(78)
      acc186(103)=-acc186(79)+acc186(103)
      acc186(103)=Qspvak2e5*acc186(103)
      acc186(104)=Qspvae1l3*acc186(82)
      acc186(104)=-acc186(49)+acc186(104)
      acc186(104)=Qspval3e5*acc186(104)
      acc186(105)=acc186(82)*Qspvae1e5
      acc186(105)=acc186(60)+acc186(105)
      acc186(105)=QspQ*acc186(105)
      acc186(97)=acc186(105)+acc186(104)+acc186(103)+acc186(102)+acc186(101)+ac&
      &c186(100)+acc186(99)+acc186(98)+acc186(32)+acc186(97)
      acc186(97)=Qspvae5e2*acc186(97)
      acc186(98)=Qspvak5k1*acc186(75)
      acc186(99)=Qspvae1k5*acc186(23)
      acc186(100)=Qspvae5k1*acc186(86)
      acc186(101)=Qspvae1e5*acc186(81)
      acc186(102)=-Qspvak1e5*acc186(83)
      acc186(102)=acc186(71)+acc186(102)
      acc186(102)=Qspvae1k1*acc186(102)
      acc186(103)=Qspvae1k2*acc186(65)
      acc186(104)=Qspvae1l3*acc186(90)
      acc186(105)=Qspvae1k2*acc186(13)
      acc186(105)=acc186(29)+acc186(105)
      acc186(105)=Qspvak2e5*acc186(105)
      acc186(106)=Qspvae5k2*acc186(67)
      acc186(107)=Qspvae1l3*acc186(83)
      acc186(107)=acc186(87)+acc186(107)
      acc186(107)=Qspval3e5*acc186(107)
      acc186(108)=Qspvae5l3*acc186(45)
      acc186(109)=Qspvae1e5*acc186(83)
      acc186(109)=acc186(42)+acc186(109)
      acc186(109)=QspQ*acc186(109)
      acc186(97)=acc186(97)+acc186(109)+acc186(108)+acc186(107)+acc186(106)+acc&
      &186(105)+acc186(104)+acc186(103)+acc186(102)+acc186(101)+acc186(100)+acc&
      &186(99)+acc186(10)+acc186(98)
      acc186(97)=Qspvae2e1*acc186(97)
      acc186(98)=Qspvae1k5*acc186(39)
      acc186(99)=Qspvae5k1*acc186(66)
      acc186(100)=Qspvae1k1*acc186(31)
      acc186(101)=-Qspvae1k2*acc186(64)
      acc186(102)=Qspvae1l3*acc186(89)
      acc186(103)=-Qspvae5k2*acc186(68)
      acc186(104)=-Qspvae5l3*acc186(72)
      acc186(105)=QspQ*acc186(70)
      acc186(98)=acc186(105)+acc186(104)+acc186(103)+acc186(102)+acc186(101)+ac&
      &c186(100)+acc186(99)+acc186(63)+acc186(98)
      acc186(98)=Qspvae2e5*acc186(98)
      acc186(99)=Qspvak5e1*acc186(43)
      acc186(100)=Qspvak1e5*acc186(50)
      acc186(101)=Qspvak1e1*acc186(26)
      acc186(102)=-Qspvak2e1*acc186(52)
      acc186(103)=Qspval3e1*acc186(92)
      acc186(104)=-Qspvak2e5*acc186(58)
      acc186(105)=Qspval3e5*acc186(84)
      acc186(106)=QspQ*acc186(56)
      acc186(99)=acc186(106)+acc186(105)+acc186(104)+acc186(103)+acc186(102)+ac&
      &c186(101)+acc186(100)+acc186(38)+acc186(99)
      acc186(99)=Qspvae5e2*acc186(99)
      acc186(100)=Qspvae1e5*acc186(73)
      acc186(101)=Qspvae5e1*acc186(41)
      acc186(100)=acc186(101)+acc186(5)+acc186(100)
      acc186(100)=QspQ*acc186(100)
      acc186(101)=Qspvak1k2*acc186(18)
      acc186(102)=Qspvak2k1*acc186(30)
      acc186(103)=Qspvak1k5*acc186(8)
      acc186(104)=Qspvak5k1*acc186(25)
      acc186(105)=Qspvak5e1*acc186(28)
      acc186(106)=Qspvae1k5*acc186(24)
      acc186(107)=Qspvak1e5*acc186(14)
      acc186(108)=Qspvae5k1*acc186(20)
      acc186(109)=Qspvae1e5*acc186(11)
      acc186(110)=Qspvae5e1*acc186(17)
      acc186(111)=-Qspvae5k1*acc186(41)
      acc186(111)=acc186(15)+acc186(111)
      acc186(111)=Qspvak1e1*acc186(111)
      acc186(112)=-Qspvak1e5*acc186(73)
      acc186(112)=acc186(21)+acc186(112)
      acc186(112)=Qspvae1k1*acc186(112)
      acc186(113)=Qspvak2e1*acc186(53)
      acc186(114)=Qspvae1k2*acc186(94)
      acc186(115)=Qspval3e1*acc186(91)
      acc186(116)=Qspvae1l3*acc186(47)
      acc186(117)=Qspvae1k2*acc186(4)
      acc186(117)=acc186(6)+acc186(117)
      acc186(117)=Qspvak2e5*acc186(117)
      acc186(118)=Qspvak2e1*acc186(1)
      acc186(118)=acc186(37)+acc186(118)
      acc186(118)=Qspvae5k2*acc186(118)
      acc186(119)=Qspvae1l3*acc186(73)
      acc186(119)=acc186(22)+acc186(119)
      acc186(119)=Qspval3e5*acc186(119)
      acc186(120)=Qspval3e1*acc186(41)
      acc186(120)=acc186(44)+acc186(120)
      acc186(120)=Qspvae5l3*acc186(120)
      brack=acc186(3)+acc186(96)+acc186(97)+acc186(98)+acc186(99)+acc186(100)+a&
      &cc186(101)+acc186(102)+acc186(103)+acc186(104)+acc186(105)+acc186(106)+a&
      &cc186(107)+acc186(108)+acc186(109)+acc186(110)+acc186(111)+acc186(112)+a&
      &cc186(113)+acc186(114)+acc186(115)+acc186(116)+acc186(117)+acc186(118)+a&
      &cc186(119)+acc186(120)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ subroutine numerator_ninja:
   subroutine numerator_ninja(ncut, Q_ext, mu2_ext, numerator) &
   & bind(c, name="p9_part21part21_part25part25part21_d186h0l1_ninja")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki_nin
      use p9_part21part21_part25part25part21_globalsl1, only: epspow
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_abbrevd186h0
      implicit none
      integer(c_int), intent(in) :: ncut
      complex(ki_nin), dimension(0:3), intent(in) :: Q_ext
      complex(ki_nin), intent(in) :: mu2_ext
      complex(ki_nin), intent(out) :: numerator
      complex(ki) :: d186
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = k2
      Q(1:4)  =cmplx(real(+Q_ext(0:3)  -qshift(:),  ki_nin), aimag(+Q_ext(0:3))&
      &, ki)
      d186 = 0.0_ki
      d186 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d186, ki), aimag(d186), ki_nin)
   end subroutine numerator_ninja
   !------#] subroutine numerator_ninja:
!---#] numerator interfaces:
end module p9_part21part21_part25part25part21_d186h0l1
