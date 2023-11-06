module     p9_part21part21_part25part25part21_d202h0l1
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p9_part21part &
   ! &21_part25part25part21/helicity0d202h0l1.f90
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
      use p9_part21part21_part25part25part21_abbrevd202h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc202(99)
      complex(ki) :: Qspvae2k2
      complex(ki) :: Qspvak2e1
      complex(ki) :: Qspvae1k2
      complex(ki) :: Qspval3e1
      complex(ki) :: QspQ
      complex(ki) :: Qspvae2e1
      complex(ki) :: Qspvae1l3
      complex(ki) :: Qspvae1l4
      complex(ki) :: Qspvae2l3
      complex(ki) :: Qspvae5e2
      complex(ki) :: Qspvak2e2
      complex(ki) :: Qspval4e1
      complex(ki) :: Qspval3e2
      complex(ki) :: Qspvae2e5
      complex(ki) :: Qspvae1e5
      complex(ki) :: Qspvae5e1
      complex(ki) :: Qspvae1e2
      Qspvae2k2 = dotproduct(Q,spvae2k2)
      Qspvak2e1 = dotproduct(Q,spvak2e1)
      Qspvae1k2 = dotproduct(Q,spvae1k2)
      Qspval3e1 = dotproduct(Q,spval3e1)
      QspQ = dotproduct(Q,Q)
      Qspvae2e1 = dotproduct(Q,spvae2e1)
      Qspvae1l3 = dotproduct(Q,spvae1l3)
      Qspvae1l4 = dotproduct(Q,spvae1l4)
      Qspvae2l3 = dotproduct(Q,spvae2l3)
      Qspvae5e2 = dotproduct(Q,spvae5e2)
      Qspvak2e2 = dotproduct(Q,spvak2e2)
      Qspval4e1 = dotproduct(Q,spval4e1)
      Qspval3e2 = dotproduct(Q,spval3e2)
      Qspvae2e5 = dotproduct(Q,spvae2e5)
      Qspvae1e5 = dotproduct(Q,spvae1e5)
      Qspvae5e1 = dotproduct(Q,spvae5e1)
      Qspvae1e2 = dotproduct(Q,spvae1e2)
      acc202(1)=abb202(6)
      acc202(2)=abb202(7)
      acc202(3)=abb202(8)
      acc202(4)=abb202(9)
      acc202(5)=abb202(10)
      acc202(6)=abb202(11)
      acc202(7)=abb202(12)
      acc202(8)=abb202(13)
      acc202(9)=abb202(14)
      acc202(10)=abb202(15)
      acc202(11)=abb202(16)
      acc202(12)=abb202(17)
      acc202(13)=abb202(18)
      acc202(14)=abb202(19)
      acc202(15)=abb202(20)
      acc202(16)=abb202(21)
      acc202(17)=abb202(22)
      acc202(18)=abb202(23)
      acc202(19)=abb202(24)
      acc202(20)=abb202(25)
      acc202(21)=abb202(26)
      acc202(22)=abb202(27)
      acc202(23)=abb202(28)
      acc202(24)=abb202(29)
      acc202(25)=abb202(30)
      acc202(26)=abb202(31)
      acc202(27)=abb202(32)
      acc202(28)=abb202(33)
      acc202(29)=abb202(34)
      acc202(30)=abb202(35)
      acc202(31)=abb202(36)
      acc202(32)=abb202(37)
      acc202(33)=abb202(38)
      acc202(34)=abb202(39)
      acc202(35)=abb202(41)
      acc202(36)=abb202(42)
      acc202(37)=abb202(43)
      acc202(38)=abb202(44)
      acc202(39)=abb202(45)
      acc202(40)=abb202(46)
      acc202(41)=abb202(47)
      acc202(42)=abb202(48)
      acc202(43)=abb202(49)
      acc202(44)=abb202(50)
      acc202(45)=abb202(51)
      acc202(46)=abb202(52)
      acc202(47)=abb202(53)
      acc202(48)=abb202(54)
      acc202(49)=abb202(55)
      acc202(50)=abb202(57)
      acc202(51)=abb202(58)
      acc202(52)=abb202(59)
      acc202(53)=abb202(60)
      acc202(54)=abb202(62)
      acc202(55)=abb202(63)
      acc202(56)=abb202(65)
      acc202(57)=abb202(70)
      acc202(58)=abb202(71)
      acc202(59)=abb202(76)
      acc202(60)=abb202(81)
      acc202(61)=abb202(82)
      acc202(62)=abb202(83)
      acc202(63)=abb202(84)
      acc202(64)=abb202(85)
      acc202(65)=abb202(87)
      acc202(66)=abb202(90)
      acc202(67)=abb202(91)
      acc202(68)=abb202(92)
      acc202(69)=abb202(93)
      acc202(70)=abb202(94)
      acc202(71)=abb202(95)
      acc202(72)=abb202(96)
      acc202(73)=abb202(97)
      acc202(74)=acc202(17)*Qspvae2k2
      acc202(75)=acc202(21)*Qspvak2e1
      acc202(76)=acc202(23)*Qspvae1k2
      acc202(77)=Qspvae2k2*Qspvak2e1
      acc202(78)=acc202(33)*acc202(77)
      acc202(79)=acc202(34)*Qspval3e1
      acc202(80)=acc202(44)*QspQ
      acc202(81)=acc202(47)*Qspvae2e1
      acc202(82)=acc202(54)*Qspvae1l3
      acc202(83)=acc202(65)*Qspvae1l4
      acc202(84)=acc202(71)*Qspvae2l3
      acc202(74)=acc202(84)+acc202(83)+acc202(82)+acc202(81)+acc202(80)+acc202(&
      &38)+acc202(79)+acc202(78)+acc202(76)+acc202(75)+acc202(74)
      acc202(74)=Qspvae5e2*acc202(74)
      acc202(75)=acc202(11)*QspQ
      acc202(76)=acc202(18)*Qspvae1l3
      acc202(78)=acc202(25)*Qspvae1k2
      acc202(79)=Qspvak2e2*Qspvae1k2
      acc202(80)=acc202(26)*acc202(79)
      acc202(81)=acc202(35)*Qspvak2e1
      acc202(82)=acc202(49)*Qspval3e1
      acc202(83)=acc202(51)*Qspvak2e2
      acc202(84)=acc202(68)*Qspval4e1
      acc202(85)=acc202(73)*Qspval3e2
      acc202(75)=acc202(85)+acc202(84)+acc202(83)+acc202(82)+acc202(81)+acc202(&
      &80)+acc202(78)+acc202(76)+acc202(75)+acc202(5)
      acc202(75)=Qspvae2e5*acc202(75)
      acc202(76)=acc202(6)*Qspval3e1
      acc202(78)=acc202(8)*Qspvak2e1
      acc202(80)=acc202(24)*acc202(77)
      acc202(81)=acc202(37)*Qspval3e2
      acc202(82)=acc202(39)*Qspvak2e2
      acc202(83)=acc202(56)*Qspvae2e1
      acc202(84)=acc202(57)*QspQ
      acc202(85)=acc202(69)*Qspval4e1
      acc202(76)=acc202(85)+acc202(84)+acc202(83)+acc202(82)+acc202(81)+acc202(&
      &80)+acc202(78)+acc202(76)+acc202(3)
      acc202(76)=Qspvae1e5*acc202(76)
      acc202(78)=acc202(12)*acc202(79)
      acc202(80)=acc202(20)*Qspvae1l3
      acc202(81)=acc202(27)*Qspvae2k2
      acc202(82)=acc202(29)*Qspvae2l3
      acc202(83)=acc202(30)*Qspvae1k2
      acc202(84)=acc202(41)*QspQ
      acc202(85)=acc202(66)*Qspvae1l4
      acc202(78)=acc202(85)+acc202(84)+acc202(83)+acc202(82)+acc202(81)+acc202(&
      &80)+acc202(15)+acc202(78)
      acc202(78)=Qspvae5e1*acc202(78)
      acc202(80)=-acc202(58)*Qspvae5e2
      acc202(81)=-acc202(59)*Qspvae1e5
      acc202(80)=acc202(81)+acc202(80)-acc202(50)
      acc202(81)=Qspvae2e1*QspQ
      acc202(82)=Qspvae2l3*Qspval3e1
      acc202(81)=acc202(81)-acc202(82)
      acc202(80)=acc202(81)*acc202(80)
      acc202(82)=-acc202(62)*Qspvae2e5
      acc202(83)=-acc202(63)*Qspvae5e1
      acc202(82)=acc202(83)+acc202(82)-acc202(61)
      acc202(83)=Qspvae1e2*QspQ
      acc202(84)=Qspval3e2*Qspvae1l3
      acc202(83)=acc202(83)-acc202(84)
      acc202(82)=acc202(83)*acc202(82)
      acc202(84)=acc202(1)*Qspval3e1
      acc202(85)=acc202(7)*Qspvak2e1
      acc202(86)=acc202(36)*QspQ
      acc202(84)=acc202(86)+acc202(9)+acc202(85)+acc202(84)
      acc202(85)=Qspvae1e5*Qspvae5e2
      acc202(84)=acc202(85)*acc202(84)
      acc202(86)=acc202(31)*Qspvae1k2
      acc202(87)=acc202(32)*Qspvae1l3
      acc202(88)=acc202(46)*QspQ
      acc202(86)=acc202(53)+acc202(88)+acc202(87)+acc202(86)
      acc202(87)=Qspvae5e1*Qspvae2e5
      acc202(86)=acc202(87)*acc202(86)
      acc202(81)=-acc202(85)*acc202(81)
      acc202(83)=-acc202(87)*acc202(83)
      acc202(81)=acc202(83)+acc202(81)
      acc202(81)=acc202(60)*acc202(81)
      acc202(83)=acc202(85)*acc202(77)
      acc202(88)=acc202(87)*acc202(79)
      acc202(83)=acc202(83)+acc202(88)
      acc202(83)=acc202(55)*acc202(83)
      acc202(88)=acc202(28)*Qspvae5e1
      acc202(89)=acc202(42)*Qspvae2e5
      acc202(88)=acc202(43)+acc202(89)+acc202(88)
      acc202(88)=Qspvae1e2*acc202(88)
      acc202(85)=-Qspvae2e1*acc202(85)
      acc202(87)=-Qspvae1e2*acc202(87)
      acc202(85)=acc202(85)+acc202(87)
      acc202(85)=acc202(2)*acc202(85)
      acc202(87)=acc202(4)*Qspvae2k2
      acc202(89)=acc202(10)*Qspvae1k2
      acc202(90)=acc202(13)*Qspvak2e1
      acc202(91)=acc202(16)*QspQ
      acc202(77)=acc202(19)*acc202(77)
      acc202(92)=acc202(22)*Qspvae1l3
      acc202(93)=acc202(40)*Qspval3e1
      acc202(94)=acc202(45)*Qspvak2e2
      acc202(95)=acc202(48)*Qspvae2e1
      acc202(79)=acc202(52)*acc202(79)
      acc202(96)=acc202(64)*Qspvae1l4
      acc202(97)=acc202(67)*Qspval4e1
      acc202(98)=acc202(70)*Qspvae2l3
      acc202(99)=acc202(72)*Qspval3e2
      brack=acc202(14)+acc202(74)+acc202(75)+acc202(76)+acc202(77)+acc202(78)+a&
      &cc202(79)+acc202(80)+acc202(81)+acc202(82)+acc202(83)+acc202(84)+acc202(&
      &85)+acc202(86)+acc202(87)+acc202(88)+acc202(89)+acc202(90)+acc202(91)+ac&
      &c202(92)+acc202(93)+acc202(94)+acc202(95)+acc202(96)+acc202(97)+acc202(9&
      &8)+acc202(99)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ subroutine numerator_ninja:
   subroutine numerator_ninja(ncut, Q_ext, mu2_ext, numerator) &
   & bind(c, name="p9_part21part21_part25part25part21_d202h0l1_ninja")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki_nin
      use p9_part21part21_part25part25part21_globalsl1, only: epspow
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_abbrevd202h0
      implicit none
      integer(c_int), intent(in) :: ncut
      complex(ki_nin), dimension(0:3), intent(in) :: Q_ext
      complex(ki_nin), intent(in) :: mu2_ext
      complex(ki_nin), intent(out) :: numerator
      complex(ki) :: d202
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = -k3-k4
      Q(1:4)  =cmplx(real(-Q_ext(0:3)  -qshift(:),  ki_nin), aimag(-Q_ext(0:3))&
      &, ki)
      d202 = 0.0_ki
      d202 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d202, ki), aimag(d202), ki_nin)
   end subroutine numerator_ninja
   !------#] subroutine numerator_ninja:
!---#] numerator interfaces:
end module p9_part21part21_part25part25part21_d202h0l1
