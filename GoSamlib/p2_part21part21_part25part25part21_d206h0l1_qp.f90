module     p2_part21part21_part25part25part21_d206h0l1_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d206h0l1_qp.f90
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
      use p2_part21part21_part25part25part21_abbrevd206h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc206(99)
      complex(ki) :: Qspvae1l3
      complex(ki) :: Qspvak2e1
      complex(ki) :: Qspvae1k2
      complex(ki) :: Qspvae2k2
      complex(ki) :: Qspval4e1
      complex(ki) :: QspQ
      complex(ki) :: Qspvae2e1
      complex(ki) :: Qspvae2l4
      complex(ki) :: Qspvae1l4
      complex(ki) :: Qspvae5e2
      complex(ki) :: Qspvak2e2
      complex(ki) :: Qspval3e1
      complex(ki) :: Qspval4e2
      complex(ki) :: Qspvae2e5
      complex(ki) :: Qspvae1e5
      complex(ki) :: Qspvae5e1
      complex(ki) :: Qspvae1e2
      Qspvae1l3 = dotproduct(Q,spvae1l3)
      Qspvak2e1 = dotproduct(Q,spvak2e1)
      Qspvae1k2 = dotproduct(Q,spvae1k2)
      Qspvae2k2 = dotproduct(Q,spvae2k2)
      Qspval4e1 = dotproduct(Q,spval4e1)
      QspQ = dotproduct(Q,Q)
      Qspvae2e1 = dotproduct(Q,spvae2e1)
      Qspvae2l4 = dotproduct(Q,spvae2l4)
      Qspvae1l4 = dotproduct(Q,spvae1l4)
      Qspvae5e2 = dotproduct(Q,spvae5e2)
      Qspvak2e2 = dotproduct(Q,spvak2e2)
      Qspval3e1 = dotproduct(Q,spval3e1)
      Qspval4e2 = dotproduct(Q,spval4e2)
      Qspvae2e5 = dotproduct(Q,spvae2e5)
      Qspvae1e5 = dotproduct(Q,spvae1e5)
      Qspvae5e1 = dotproduct(Q,spvae5e1)
      Qspvae1e2 = dotproduct(Q,spvae1e2)
      acc206(1)=abb206(6)
      acc206(2)=abb206(7)
      acc206(3)=abb206(8)
      acc206(4)=abb206(9)
      acc206(5)=abb206(10)
      acc206(6)=abb206(11)
      acc206(7)=abb206(12)
      acc206(8)=abb206(13)
      acc206(9)=abb206(14)
      acc206(10)=abb206(15)
      acc206(11)=abb206(16)
      acc206(12)=abb206(17)
      acc206(13)=abb206(18)
      acc206(14)=abb206(19)
      acc206(15)=abb206(20)
      acc206(16)=abb206(21)
      acc206(17)=abb206(22)
      acc206(18)=abb206(23)
      acc206(19)=abb206(24)
      acc206(20)=abb206(25)
      acc206(21)=abb206(26)
      acc206(22)=abb206(27)
      acc206(23)=abb206(28)
      acc206(24)=abb206(29)
      acc206(25)=abb206(30)
      acc206(26)=abb206(31)
      acc206(27)=abb206(32)
      acc206(28)=abb206(33)
      acc206(29)=abb206(34)
      acc206(30)=abb206(35)
      acc206(31)=abb206(36)
      acc206(32)=abb206(37)
      acc206(33)=abb206(38)
      acc206(34)=abb206(39)
      acc206(35)=abb206(41)
      acc206(36)=abb206(42)
      acc206(37)=abb206(43)
      acc206(38)=abb206(44)
      acc206(39)=abb206(45)
      acc206(40)=abb206(46)
      acc206(41)=abb206(47)
      acc206(42)=abb206(48)
      acc206(43)=abb206(49)
      acc206(44)=abb206(50)
      acc206(45)=abb206(51)
      acc206(46)=abb206(52)
      acc206(47)=abb206(53)
      acc206(48)=abb206(54)
      acc206(49)=abb206(55)
      acc206(50)=abb206(57)
      acc206(51)=abb206(58)
      acc206(52)=abb206(59)
      acc206(53)=abb206(60)
      acc206(54)=abb206(62)
      acc206(55)=abb206(63)
      acc206(56)=abb206(65)
      acc206(57)=abb206(66)
      acc206(58)=abb206(67)
      acc206(59)=abb206(68)
      acc206(60)=abb206(69)
      acc206(61)=abb206(70)
      acc206(62)=abb206(71)
      acc206(63)=abb206(72)
      acc206(64)=abb206(76)
      acc206(65)=abb206(81)
      acc206(66)=abb206(82)
      acc206(67)=abb206(83)
      acc206(68)=abb206(84)
      acc206(69)=abb206(85)
      acc206(70)=abb206(87)
      acc206(71)=abb206(89)
      acc206(72)=abb206(90)
      acc206(73)=abb206(91)
      acc206(74)=acc206(17)*Qspvae1l3
      acc206(75)=acc206(21)*Qspvak2e1
      acc206(76)=acc206(23)*Qspvae1k2
      acc206(77)=Qspvae2k2*Qspvak2e1
      acc206(78)=acc206(33)*acc206(77)
      acc206(79)=acc206(34)*Qspval4e1
      acc206(80)=acc206(44)*QspQ
      acc206(81)=acc206(47)*Qspvae2e1
      acc206(82)=acc206(59)*Qspvae2k2
      acc206(83)=acc206(70)*Qspvae2l4
      acc206(84)=acc206(71)*Qspvae1l4
      acc206(74)=acc206(84)+acc206(83)+acc206(82)+acc206(81)+acc206(80)+acc206(&
      &38)+acc206(79)+acc206(78)+acc206(76)+acc206(75)+acc206(74)
      acc206(74)=Qspvae5e2*acc206(74)
      acc206(75)=acc206(11)*QspQ
      acc206(76)=acc206(18)*Qspvae1l4
      acc206(78)=acc206(25)*Qspvae1k2
      acc206(79)=Qspvak2e2*Qspvae1k2
      acc206(80)=acc206(26)*acc206(79)
      acc206(81)=acc206(35)*Qspvak2e1
      acc206(82)=acc206(51)*Qspval4e1
      acc206(83)=acc206(54)*Qspval3e1
      acc206(84)=acc206(63)*Qspvak2e2
      acc206(85)=acc206(73)*Qspval4e2
      acc206(75)=acc206(85)+acc206(84)+acc206(83)+acc206(82)+acc206(81)+acc206(&
      &80)+acc206(78)+acc206(76)+acc206(75)+acc206(5)
      acc206(75)=Qspvae2e5*acc206(75)
      acc206(76)=acc206(6)*Qspval4e1
      acc206(78)=acc206(8)*Qspvak2e1
      acc206(80)=acc206(24)*acc206(77)
      acc206(81)=acc206(37)*Qspval4e2
      acc206(82)=acc206(39)*Qspvak2e2
      acc206(83)=acc206(55)*Qspval3e1
      acc206(84)=acc206(56)*Qspvae2e1
      acc206(85)=acc206(61)*QspQ
      acc206(76)=acc206(85)+acc206(84)+acc206(83)+acc206(82)+acc206(81)+acc206(&
      &80)+acc206(78)+acc206(76)+acc206(3)
      acc206(76)=Qspvae1e5*acc206(76)
      acc206(78)=acc206(12)*acc206(79)
      acc206(80)=acc206(20)*Qspvae1l4
      acc206(81)=acc206(27)*Qspvae1l3
      acc206(82)=acc206(29)*Qspvae2l4
      acc206(83)=acc206(30)*Qspvae1k2
      acc206(84)=acc206(41)*QspQ
      acc206(85)=acc206(45)*Qspvae2k2
      acc206(78)=acc206(85)+acc206(84)+acc206(83)+acc206(82)+acc206(81)+acc206(&
      &80)+acc206(15)+acc206(78)
      acc206(78)=Qspvae5e1*acc206(78)
      acc206(80)=-acc206(62)*Qspvae5e2
      acc206(81)=-acc206(64)*Qspvae1e5
      acc206(80)=acc206(81)+acc206(80)-acc206(50)
      acc206(81)=Qspvae2e1*QspQ
      acc206(82)=Qspvae2l4*Qspval4e1
      acc206(81)=acc206(81)-acc206(82)
      acc206(80)=acc206(81)*acc206(80)
      acc206(82)=-acc206(67)*Qspvae2e5
      acc206(83)=-acc206(68)*Qspvae5e1
      acc206(82)=acc206(83)+acc206(82)-acc206(66)
      acc206(83)=Qspvae1e2*QspQ
      acc206(84)=Qspval4e2*Qspvae1l4
      acc206(83)=acc206(83)-acc206(84)
      acc206(82)=acc206(83)*acc206(82)
      acc206(84)=acc206(1)*Qspval4e1
      acc206(85)=acc206(7)*Qspvak2e1
      acc206(86)=acc206(36)*QspQ
      acc206(84)=acc206(86)+acc206(9)+acc206(85)+acc206(84)
      acc206(85)=Qspvae1e5*Qspvae5e2
      acc206(84)=acc206(85)*acc206(84)
      acc206(86)=acc206(31)*Qspvae1k2
      acc206(87)=acc206(32)*Qspvae1l4
      acc206(88)=acc206(46)*QspQ
      acc206(86)=acc206(53)+acc206(88)+acc206(87)+acc206(86)
      acc206(87)=Qspvae5e1*Qspvae2e5
      acc206(86)=acc206(87)*acc206(86)
      acc206(81)=-acc206(85)*acc206(81)
      acc206(83)=-acc206(87)*acc206(83)
      acc206(81)=acc206(83)+acc206(81)
      acc206(81)=acc206(65)*acc206(81)
      acc206(83)=acc206(85)*acc206(77)
      acc206(88)=acc206(87)*acc206(79)
      acc206(83)=acc206(83)+acc206(88)
      acc206(83)=acc206(57)*acc206(83)
      acc206(88)=acc206(28)*Qspvae5e1
      acc206(89)=acc206(42)*Qspvae2e5
      acc206(88)=acc206(43)+acc206(89)+acc206(88)
      acc206(88)=Qspvae1e2*acc206(88)
      acc206(85)=-Qspvae2e1*acc206(85)
      acc206(87)=-Qspvae1e2*acc206(87)
      acc206(85)=acc206(85)+acc206(87)
      acc206(85)=acc206(2)*acc206(85)
      acc206(87)=acc206(4)*Qspvae1l3
      acc206(89)=acc206(10)*Qspvae1k2
      acc206(90)=acc206(13)*Qspvak2e1
      acc206(91)=acc206(16)*QspQ
      acc206(77)=acc206(19)*acc206(77)
      acc206(92)=acc206(22)*Qspvae1l4
      acc206(93)=acc206(40)*Qspval4e1
      acc206(94)=acc206(48)*Qspvae2e1
      acc206(95)=acc206(49)*Qspval3e1
      acc206(79)=acc206(52)*acc206(79)
      acc206(96)=acc206(58)*Qspvae2k2
      acc206(97)=acc206(60)*Qspvak2e2
      acc206(98)=acc206(69)*Qspvae2l4
      acc206(99)=acc206(72)*Qspval4e2
      brack=acc206(14)+acc206(74)+acc206(75)+acc206(76)+acc206(77)+acc206(78)+a&
      &cc206(79)+acc206(80)+acc206(81)+acc206(82)+acc206(83)+acc206(84)+acc206(&
      &85)+acc206(86)+acc206(87)+acc206(88)+acc206(89)+acc206(90)+acc206(91)+ac&
      &c206(92)+acc206(93)+acc206(94)+acc206(95)+acc206(96)+acc206(97)+acc206(9&
      &8)+acc206(99)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ subroutine numerator_ninja:
   subroutine numerator_ninja(ncut, Q_ext, mu2_ext, numerator) &
   & bind(c, name="p2_part21part21_part25part25part21_d206h0l1_qp_ninja")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki_nin
      use p2_part21part21_part25part25part21_globalsl1_qp, only: epspow
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_abbrevd206h0_qp
      implicit none
      integer(c_int), intent(in) :: ncut
      complex(ki_nin), dimension(0:3), intent(in) :: Q_ext
      complex(ki_nin), intent(in) :: mu2_ext
      complex(ki_nin), intent(out) :: numerator
      complex(ki) :: d206
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = -k3-k4
      Q(1:4)  =cmplx(real(-Q_ext(0:3)  -qshift(:),  ki_nin), aimag(-Q_ext(0:3))&
      &, ki)
      d206 = 0.0_ki
      d206 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d206, ki), aimag(d206), ki_nin)
   end subroutine numerator_ninja
   !------#] subroutine numerator_ninja:
!---#] numerator interfaces:
end module p2_part21part21_part25part25part21_d206h0l1_qp
