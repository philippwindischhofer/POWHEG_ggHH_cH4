module     p2_part21part21_part25part25part21_d184h0l1_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d184h0l1_qp.f90
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
      use p2_part21part21_part25part25part21_abbrevd184h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc184(120)
      complex(ki) :: Qspvak2e1
      complex(ki) :: Qspvak1e1
      complex(ki) :: Qspvak1e5
      complex(ki) :: Qspvae2k2
      complex(ki) :: Qspval3e1
      complex(ki) :: Qspvae5k2
      complex(ki) :: Qspvae2l4
      complex(ki) :: QspQ
      complex(ki) :: Qspvae5e1
      complex(ki) :: Qspvae2e5
      complex(ki) :: Qspvae1e2
      complex(ki) :: Qspvae1k2
      complex(ki) :: Qspvae1k1
      complex(ki) :: Qspvae5k1
      complex(ki) :: Qspvak2e2
      complex(ki) :: Qspvae1l3
      complex(ki) :: Qspvak2e5
      complex(ki) :: Qspval4e2
      complex(ki) :: Qspvae1e5
      complex(ki) :: Qspvae5e2
      complex(ki) :: Qspvae2e1
      Qspvak2e1 = dotproduct(Q,spvak2e1)
      Qspvak1e1 = dotproduct(Q,spvak1e1)
      Qspvak1e5 = dotproduct(Q,spvak1e5)
      Qspvae2k2 = dotproduct(Q,spvae2k2)
      Qspval3e1 = dotproduct(Q,spval3e1)
      Qspvae5k2 = dotproduct(Q,spvae5k2)
      Qspvae2l4 = dotproduct(Q,spvae2l4)
      QspQ = dotproduct(Q,Q)
      Qspvae5e1 = dotproduct(Q,spvae5e1)
      Qspvae2e5 = dotproduct(Q,spvae2e5)
      Qspvae1e2 = dotproduct(Q,spvae1e2)
      Qspvae1k2 = dotproduct(Q,spvae1k2)
      Qspvae1k1 = dotproduct(Q,spvae1k1)
      Qspvae5k1 = dotproduct(Q,spvae5k1)
      Qspvak2e2 = dotproduct(Q,spvak2e2)
      Qspvae1l3 = dotproduct(Q,spvae1l3)
      Qspvak2e5 = dotproduct(Q,spvak2e5)
      Qspval4e2 = dotproduct(Q,spval4e2)
      Qspvae1e5 = dotproduct(Q,spvae1e5)
      Qspvae5e2 = dotproduct(Q,spvae5e2)
      Qspvae2e1 = dotproduct(Q,spvae2e1)
      acc184(1)=abb184(6)
      acc184(2)=abb184(7)
      acc184(3)=abb184(8)
      acc184(4)=abb184(9)
      acc184(5)=abb184(10)
      acc184(6)=abb184(11)
      acc184(7)=abb184(12)
      acc184(8)=abb184(13)
      acc184(9)=abb184(14)
      acc184(10)=abb184(15)
      acc184(11)=abb184(16)
      acc184(12)=abb184(17)
      acc184(13)=abb184(18)
      acc184(14)=abb184(19)
      acc184(15)=abb184(20)
      acc184(16)=abb184(21)
      acc184(17)=abb184(22)
      acc184(18)=abb184(23)
      acc184(19)=abb184(24)
      acc184(20)=abb184(25)
      acc184(21)=abb184(26)
      acc184(22)=abb184(27)
      acc184(23)=abb184(28)
      acc184(24)=abb184(29)
      acc184(25)=abb184(30)
      acc184(26)=abb184(31)
      acc184(27)=abb184(32)
      acc184(28)=abb184(33)
      acc184(29)=abb184(34)
      acc184(30)=abb184(35)
      acc184(31)=abb184(36)
      acc184(32)=abb184(37)
      acc184(33)=abb184(38)
      acc184(34)=abb184(39)
      acc184(35)=abb184(40)
      acc184(36)=abb184(41)
      acc184(37)=abb184(42)
      acc184(38)=abb184(43)
      acc184(39)=abb184(44)
      acc184(40)=abb184(45)
      acc184(41)=abb184(46)
      acc184(42)=abb184(47)
      acc184(43)=abb184(48)
      acc184(44)=abb184(49)
      acc184(45)=abb184(50)
      acc184(46)=abb184(51)
      acc184(47)=abb184(52)
      acc184(48)=abb184(53)
      acc184(49)=abb184(54)
      acc184(50)=abb184(55)
      acc184(51)=abb184(56)
      acc184(52)=abb184(57)
      acc184(53)=abb184(58)
      acc184(54)=abb184(59)
      acc184(55)=abb184(60)
      acc184(56)=abb184(61)
      acc184(57)=abb184(62)
      acc184(58)=abb184(63)
      acc184(59)=abb184(64)
      acc184(60)=abb184(65)
      acc184(61)=abb184(66)
      acc184(62)=abb184(67)
      acc184(63)=abb184(68)
      acc184(64)=abb184(69)
      acc184(65)=abb184(70)
      acc184(66)=abb184(71)
      acc184(67)=abb184(72)
      acc184(68)=abb184(73)
      acc184(69)=abb184(74)
      acc184(70)=abb184(75)
      acc184(71)=abb184(76)
      acc184(72)=abb184(77)
      acc184(73)=abb184(78)
      acc184(74)=abb184(79)
      acc184(75)=abb184(80)
      acc184(76)=abb184(81)
      acc184(77)=abb184(82)
      acc184(78)=abb184(83)
      acc184(79)=abb184(84)
      acc184(80)=abb184(85)
      acc184(81)=abb184(86)
      acc184(82)=abb184(87)
      acc184(83)=abb184(88)
      acc184(84)=abb184(89)
      acc184(85)=abb184(93)
      acc184(86)=abb184(95)
      acc184(87)=abb184(96)
      acc184(88)=abb184(98)
      acc184(89)=abb184(99)
      acc184(90)=abb184(100)
      acc184(91)=abb184(101)
      acc184(92)=abb184(102)
      acc184(93)=abb184(103)
      acc184(94)=abb184(104)
      acc184(95)=abb184(105)
      acc184(96)=abb184(107)
      acc184(97)=abb184(108)
      acc184(98)=abb184(110)
      acc184(99)=abb184(112)
      acc184(100)=Qspvak2e1*acc184(82)
      acc184(101)=Qspvak1e1*acc184(14)
      acc184(102)=Qspvak1e5*acc184(88)
      acc184(103)=Qspvae2k2*acc184(96)
      acc184(104)=Qspval3e1*acc184(60)
      acc184(105)=Qspvae5k2*acc184(90)
      acc184(106)=Qspvae2l4*acc184(40)
      acc184(107)=QspQ*acc184(47)
      acc184(100)=acc184(107)+acc184(106)+acc184(105)+acc184(104)+acc184(103)+a&
      &cc184(102)+acc184(101)+acc184(20)+acc184(100)
      acc184(100)=QspQ*acc184(100)
      acc184(101)=Qspvak2e1*acc184(17)
      acc184(102)=Qspvak1e1*acc184(66)
      acc184(103)=Qspvak1e5*acc184(72)
      acc184(104)=Qspval3e1*acc184(92)
      acc184(101)=acc184(104)+acc184(103)+acc184(102)+acc184(48)+acc184(101)
      acc184(101)=Qspvae5k2*acc184(101)
      acc184(102)=Qspvak2e1*acc184(39)
      acc184(103)=Qspvak1e1*acc184(26)
      acc184(104)=Qspvak1e5*acc184(41)
      acc184(105)=Qspval3e1*acc184(91)
      acc184(102)=acc184(105)+acc184(104)+acc184(103)+acc184(22)+acc184(102)
      acc184(102)=Qspvae2l4*acc184(102)
      acc184(103)=Qspvak1e1*acc184(80)
      acc184(104)=Qspvak1e5*acc184(67)
      acc184(103)=acc184(104)+acc184(38)+acc184(103)
      acc184(103)=Qspvae2k2*acc184(103)
      acc184(104)=Qspvak2e1*acc184(15)
      acc184(105)=Qspvak1e1*acc184(5)
      acc184(106)=Qspvak1e5*acc184(70)
      acc184(107)=Qspvae2k2*acc184(97)
      acc184(107)=acc184(37)+acc184(107)
      acc184(107)=Qspval3e1*acc184(107)
      acc184(108)=Qspvae5e1*acc184(58)
      acc184(109)=-Qspvae5e1*acc184(83)
      acc184(109)=acc184(64)+acc184(109)
      acc184(109)=Qspvae2e5*acc184(109)
      acc184(100)=acc184(100)+acc184(109)+acc184(102)+acc184(101)+acc184(108)+a&
      &cc184(107)+acc184(103)+acc184(106)+acc184(105)+acc184(19)+acc184(104)
      acc184(100)=Qspvae1e2*acc184(100)
      acc184(101)=Qspvae1k2*acc184(65)
      acc184(102)=Qspvae1k1*acc184(76)
      acc184(103)=Qspvae5k1*acc184(81)
      acc184(104)=Qspvak2e2*acc184(84)
      acc184(105)=Qspvae1l3*acc184(57)
      acc184(106)=Qspvak2e5*acc184(79)
      acc184(107)=Qspval4e2*acc184(78)
      acc184(108)=QspQ*acc184(77)
      acc184(101)=acc184(108)+acc184(107)+acc184(106)+acc184(105)+acc184(104)+a&
      &cc184(103)+acc184(102)+acc184(12)+acc184(101)
      acc184(101)=QspQ*acc184(101)
      acc184(102)=Qspvae1k2*acc184(7)
      acc184(103)=Qspvae1k1*acc184(75)
      acc184(104)=Qspvae5k1*acc184(28)
      acc184(105)=Qspvae1l3*acc184(93)
      acc184(102)=acc184(105)+acc184(104)+acc184(103)+acc184(73)+acc184(102)
      acc184(102)=Qspvak2e5*acc184(102)
      acc184(103)=Qspvae1k2*acc184(45)
      acc184(104)=Qspvae1k1*acc184(29)
      acc184(105)=Qspvae5k1*acc184(27)
      acc184(106)=Qspvae1l3*acc184(74)
      acc184(103)=acc184(106)+acc184(105)+acc184(104)+acc184(8)+acc184(103)
      acc184(103)=Qspval4e2*acc184(103)
      acc184(104)=Qspvae1k1*acc184(2)
      acc184(105)=Qspvae5k1*acc184(59)
      acc184(104)=acc184(105)+acc184(25)+acc184(104)
      acc184(104)=Qspvak2e2*acc184(104)
      acc184(105)=Qspvae1k2*acc184(4)
      acc184(106)=Qspvae1k1*acc184(54)
      acc184(107)=Qspvae5k1*acc184(6)
      acc184(108)=Qspvak2e2*acc184(99)
      acc184(108)=acc184(32)+acc184(108)
      acc184(108)=Qspvae1l3*acc184(108)
      acc184(109)=Qspvae1e5*acc184(71)
      acc184(110)=-Qspvae1e5*acc184(83)
      acc184(110)=acc184(46)+acc184(110)
      acc184(110)=Qspvae5e2*acc184(110)
      acc184(101)=acc184(101)+acc184(110)+acc184(103)+acc184(102)+acc184(109)+a&
      &cc184(108)+acc184(104)+acc184(107)+acc184(106)+acc184(1)+acc184(105)
      acc184(101)=Qspvae2e1*acc184(101)
      acc184(102)=Qspvae1k2*acc184(95)
      acc184(103)=Qspvae1k1*acc184(16)
      acc184(104)=Qspvae5k1*acc184(23)
      acc184(105)=Qspvae1l3*acc184(55)
      acc184(102)=acc184(105)+acc184(104)+acc184(103)+acc184(11)+acc184(102)
      acc184(102)=Qspvae2e5*acc184(102)
      acc184(103)=Qspvak2e1*acc184(3)
      acc184(104)=Qspvak1e1*acc184(56)
      acc184(105)=Qspvak1e5*acc184(52)
      acc184(106)=Qspval3e1*acc184(89)
      acc184(103)=acc184(106)+acc184(105)+acc184(104)+acc184(50)+acc184(103)
      acc184(103)=Qspvae5e2*acc184(103)
      acc184(104)=Qspvae1e5*acc184(69)
      acc184(105)=Qspvae5e1*acc184(68)
      acc184(106)=Qspvae2e5*acc184(63)
      acc184(107)=Qspvae5e2*acc184(51)
      acc184(104)=acc184(107)+acc184(106)+acc184(105)+acc184(49)+acc184(104)
      acc184(104)=QspQ*acc184(104)
      acc184(105)=Qspvak2e1*acc184(34)
      acc184(106)=Qspvae1k2*acc184(30)
      acc184(107)=Qspvak1e1*acc184(42)
      acc184(108)=Qspvae1k1*acc184(61)
      acc184(109)=Qspvak1e5*acc184(43)
      acc184(110)=Qspvae5k1*acc184(24)
      acc184(111)=Qspvak2e2*acc184(44)
      acc184(112)=Qspvae2k2*acc184(62)
      acc184(113)=Qspval3e1*acc184(87)
      acc184(114)=Qspvae1l3*acc184(33)
      acc184(115)=-Qspvae2k2*acc184(94)
      acc184(115)=acc184(21)+acc184(115)
      acc184(115)=Qspvae1e5*acc184(115)
      acc184(116)=-Qspvak2e2*acc184(98)
      acc184(116)=acc184(13)+acc184(116)
      acc184(116)=Qspvae5e1*acc184(116)
      acc184(117)=Qspvae5e1*acc184(35)
      acc184(117)=acc184(9)+acc184(117)
      acc184(117)=Qspvak2e5*acc184(117)
      acc184(118)=Qspvae1e5*acc184(36)
      acc184(118)=acc184(18)+acc184(118)
      acc184(118)=Qspvae5k2*acc184(118)
      acc184(119)=Qspvae5e1*acc184(86)
      acc184(119)=acc184(53)+acc184(119)
      acc184(119)=Qspval4e2*acc184(119)
      acc184(120)=Qspvae1e5*acc184(85)
      acc184(120)=acc184(31)+acc184(120)
      acc184(120)=Qspvae2l4*acc184(120)
      brack=acc184(10)+acc184(100)+acc184(101)+acc184(102)+acc184(103)+acc184(1&
      &04)+acc184(105)+acc184(106)+acc184(107)+acc184(108)+acc184(109)+acc184(1&
      &10)+acc184(111)+acc184(112)+acc184(113)+acc184(114)+acc184(115)+acc184(1&
      &16)+acc184(117)+acc184(118)+acc184(119)+acc184(120)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ subroutine numerator_ninja:
   subroutine numerator_ninja(ncut, Q_ext, mu2_ext, numerator) &
   & bind(c, name="p2_part21part21_part25part25part21_d184h0l1_qp_ninja")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki_nin
      use p2_part21part21_part25part25part21_globalsl1_qp, only: epspow
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_abbrevd184h0_qp
      implicit none
      integer(c_int), intent(in) :: ncut
      complex(ki_nin), dimension(0:3), intent(in) :: Q_ext
      complex(ki_nin), intent(in) :: mu2_ext
      complex(ki_nin), intent(out) :: numerator
      complex(ki) :: d184
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = k2
      Q(1:4)  =cmplx(real(+Q_ext(0:3)  -qshift(:),  ki_nin), aimag(+Q_ext(0:3))&
      &, ki)
      d184 = 0.0_ki
      d184 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d184, ki), aimag(d184), ki_nin)
   end subroutine numerator_ninja
   !------#] subroutine numerator_ninja:
!---#] numerator interfaces:
end module p2_part21part21_part25part25part21_d184h0l1_qp
