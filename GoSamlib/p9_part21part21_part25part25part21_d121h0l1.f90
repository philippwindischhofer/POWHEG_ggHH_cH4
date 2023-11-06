module     p9_part21part21_part25part25part21_d121h0l1
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p9_part21part &
   ! &21_part25part25part21/helicity0d121h0l1.f90
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
      use p9_part21part21_part25part25part21_abbrevd121h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc121(89)
      complex(ki) :: Qspvae1k2
      complex(ki) :: Qspvae2k2
      complex(ki) :: Qspvak1k2
      complex(ki) :: Qspvak5k2
      complex(ki) :: Qspvak2e1
      complex(ki) :: Qspvak2e2
      complex(ki) :: Qspvae2l4
      complex(ki) :: Qspval4e2
      complex(ki) :: Qspvak2k1
      complex(ki) :: Qspvae2l3
      complex(ki) :: Qspval3e2
      complex(ki) :: Qspvae5k2
      complex(ki) :: Qspvak2k5
      complex(ki) :: Qspvak2e5
      complex(ki) :: QspQ
      complex(ki) :: Qspvae5e2
      complex(ki) :: Qspvae2e5
      complex(ki) :: Qspvae2e1
      complex(ki) :: Qspvae1e2
      complex(ki) :: Qspvae2k5
      complex(ki) :: Qspvak5e2
      complex(ki) :: Qspvae2k1
      complex(ki) :: Qspvak1e2
      Qspvae1k2 = dotproduct(Q,spvae1k2)
      Qspvae2k2 = dotproduct(Q,spvae2k2)
      Qspvak1k2 = dotproduct(Q,spvak1k2)
      Qspvak5k2 = dotproduct(Q,spvak5k2)
      Qspvak2e1 = dotproduct(Q,spvak2e1)
      Qspvak2e2 = dotproduct(Q,spvak2e2)
      Qspvae2l4 = dotproduct(Q,spvae2l4)
      Qspval4e2 = dotproduct(Q,spval4e2)
      Qspvak2k1 = dotproduct(Q,spvak2k1)
      Qspvae2l3 = dotproduct(Q,spvae2l3)
      Qspval3e2 = dotproduct(Q,spval3e2)
      Qspvae5k2 = dotproduct(Q,spvae5k2)
      Qspvak2k5 = dotproduct(Q,spvak2k5)
      Qspvak2e5 = dotproduct(Q,spvak2e5)
      QspQ = dotproduct(Q,Q)
      Qspvae5e2 = dotproduct(Q,spvae5e2)
      Qspvae2e5 = dotproduct(Q,spvae2e5)
      Qspvae2e1 = dotproduct(Q,spvae2e1)
      Qspvae1e2 = dotproduct(Q,spvae1e2)
      Qspvae2k5 = dotproduct(Q,spvae2k5)
      Qspvak5e2 = dotproduct(Q,spvak5e2)
      Qspvae2k1 = dotproduct(Q,spvae2k1)
      Qspvak1e2 = dotproduct(Q,spvak1e2)
      acc121(1)=abb121(7)
      acc121(2)=abb121(8)
      acc121(3)=abb121(9)
      acc121(4)=abb121(10)
      acc121(5)=abb121(11)
      acc121(6)=abb121(12)
      acc121(7)=abb121(13)
      acc121(8)=abb121(14)
      acc121(9)=abb121(15)
      acc121(10)=abb121(17)
      acc121(11)=abb121(18)
      acc121(12)=abb121(19)
      acc121(13)=abb121(20)
      acc121(14)=abb121(21)
      acc121(15)=abb121(22)
      acc121(16)=abb121(23)
      acc121(17)=abb121(24)
      acc121(18)=abb121(25)
      acc121(19)=abb121(26)
      acc121(20)=abb121(27)
      acc121(21)=abb121(28)
      acc121(22)=abb121(29)
      acc121(23)=abb121(30)
      acc121(24)=abb121(31)
      acc121(25)=abb121(32)
      acc121(26)=abb121(33)
      acc121(27)=abb121(34)
      acc121(28)=abb121(35)
      acc121(29)=abb121(36)
      acc121(30)=abb121(37)
      acc121(31)=abb121(38)
      acc121(32)=abb121(39)
      acc121(33)=abb121(40)
      acc121(34)=abb121(41)
      acc121(35)=abb121(42)
      acc121(36)=abb121(43)
      acc121(37)=abb121(44)
      acc121(38)=abb121(45)
      acc121(39)=abb121(46)
      acc121(40)=abb121(47)
      acc121(41)=abb121(48)
      acc121(42)=abb121(49)
      acc121(43)=abb121(50)
      acc121(44)=abb121(51)
      acc121(45)=abb121(52)
      acc121(46)=abb121(53)
      acc121(47)=abb121(54)
      acc121(48)=abb121(55)
      acc121(49)=abb121(56)
      acc121(50)=abb121(58)
      acc121(51)=abb121(59)
      acc121(52)=abb121(60)
      acc121(53)=abb121(61)
      acc121(54)=abb121(62)
      acc121(55)=abb121(63)
      acc121(56)=abb121(64)
      acc121(57)=abb121(65)
      acc121(58)=abb121(66)
      acc121(59)=abb121(67)
      acc121(60)=abb121(68)
      acc121(61)=abb121(69)
      acc121(62)=abb121(70)
      acc121(63)=abb121(71)
      acc121(64)=abb121(72)
      acc121(65)=abb121(73)
      acc121(66)=acc121(2)*Qspvae1k2
      acc121(67)=acc121(5)*Qspvae2k2
      acc121(68)=acc121(9)*Qspvak1k2
      acc121(69)=acc121(11)*Qspvak5k2
      acc121(70)=acc121(16)*Qspvak2e1
      acc121(71)=acc121(20)*Qspvak2e2
      acc121(72)=acc121(24)*Qspvae2l4
      acc121(73)=acc121(26)*Qspval4e2
      acc121(74)=acc121(42)*Qspvak2k1
      acc121(75)=acc121(50)*Qspvae2l3
      acc121(76)=acc121(53)*Qspval3e2
      acc121(77)=acc121(54)*Qspvae5k2
      acc121(78)=acc121(56)*Qspvak2k5
      acc121(79)=acc121(57)*Qspvak2e5
      acc121(66)=acc121(79)+acc121(78)+acc121(77)+acc121(76)+acc121(75)+acc121(&
      &74)+acc121(73)+acc121(72)+acc121(22)+acc121(71)+acc121(70)+acc121(69)+ac&
      &c121(68)+acc121(67)+acc121(66)
      acc121(66)=QspQ*acc121(66)
      acc121(67)=acc121(1)*Qspvae1k2
      acc121(68)=acc121(4)*Qspvae5k2
      acc121(69)=acc121(6)*Qspvae2l3
      acc121(70)=acc121(8)*Qspvae2k2
      acc121(71)=acc121(14)*Qspvae2l4
      acc121(72)=acc121(36)*Qspvak1k2
      acc121(73)=acc121(52)*Qspvak5k2
      acc121(67)=acc121(73)+acc121(72)+acc121(71)+acc121(70)+acc121(7)+acc121(6&
      &9)+acc121(68)+acc121(67)
      acc121(67)=Qspvak2e2*acc121(67)
      acc121(68)=acc121(15)*Qspvae1k2
      acc121(69)=acc121(21)*Qspvae2k2
      acc121(70)=acc121(30)*Qspvae2l4
      acc121(71)=acc121(34)*Qspvak1k2
      acc121(72)=acc121(48)*Qspvak5k2
      acc121(73)=acc121(51)*Qspvae5k2
      acc121(68)=acc121(73)+acc121(72)+acc121(35)+acc121(71)+acc121(70)+acc121(&
      &69)+acc121(68)
      acc121(68)=Qspval3e2*acc121(68)
      acc121(69)=acc121(12)*Qspval4e2
      acc121(70)=acc121(38)*Qspvak2e1
      acc121(71)=acc121(62)*Qspvak2k1
      acc121(72)=acc121(63)*Qspvak2e5
      acc121(73)=acc121(64)*Qspvak2k5
      acc121(69)=acc121(73)+acc121(72)+acc121(71)+acc121(41)+acc121(70)+acc121(&
      &69)
      acc121(69)=Qspvae2k2*acc121(69)
      acc121(70)=acc121(27)*Qspvak2e1
      acc121(71)=acc121(28)*Qspval4e2
      acc121(72)=acc121(49)*Qspvak2e5
      acc121(73)=acc121(60)*Qspvak2k1
      acc121(74)=acc121(61)*Qspvak2k5
      acc121(70)=acc121(74)+acc121(73)+acc121(72)+acc121(40)+acc121(71)+acc121(&
      &70)
      acc121(70)=Qspvae2l3*acc121(70)
      acc121(71)=acc121(45)*Qspvae5e2
      acc121(71)=acc121(71)+acc121(37)
      acc121(71)=Qspvae2e5*acc121(71)
      acc121(72)=acc121(31)*Qspvae2e1
      acc121(72)=acc121(72)+acc121(17)
      acc121(72)=Qspvae1e2*acc121(72)
      acc121(73)=-Qspvae2k5*acc121(43)
      acc121(73)=acc121(25)+acc121(73)
      acc121(73)=Qspvak5e2*acc121(73)
      acc121(74)=-Qspvae2k1*acc121(43)
      acc121(74)=acc121(19)+acc121(74)
      acc121(74)=Qspvak1e2*acc121(74)
      acc121(75)=acc121(3)*QspQ**2
      acc121(76)=acc121(10)*Qspvae1k2
      acc121(77)=acc121(18)*Qspvak2e1
      acc121(78)=acc121(23)*Qspvae2k1
      acc121(79)=acc121(29)*Qspvae2e1
      acc121(80)=acc121(32)*Qspvak1k2
      acc121(81)=acc121(33)*Qspvae5e2
      acc121(82)=acc121(39)*Qspvae2k5
      acc121(83)=acc121(44)*Qspvae2l4
      acc121(84)=acc121(46)*Qspvak5k2
      acc121(85)=acc121(47)*Qspval4e2
      acc121(86)=acc121(55)*Qspvae5k2
      acc121(87)=acc121(58)*Qspvak2e5
      acc121(88)=acc121(59)*Qspvak2k5
      acc121(89)=acc121(65)*Qspvak2k1
      brack=acc121(13)+acc121(66)+acc121(67)+acc121(68)+acc121(69)+acc121(70)+a&
      &cc121(71)+acc121(72)+acc121(73)+acc121(74)+acc121(75)+acc121(76)+acc121(&
      &77)+acc121(78)+acc121(79)+acc121(80)+acc121(81)+acc121(82)+acc121(83)+ac&
      &c121(84)+acc121(85)+acc121(86)+acc121(87)+acc121(88)+acc121(89)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ subroutine numerator_ninja:
   subroutine numerator_ninja(ncut, Q_ext, mu2_ext, numerator) &
   & bind(c, name="p9_part21part21_part25part25part21_d121h0l1_ninja")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki_nin
      use p9_part21part21_part25part25part21_globalsl1, only: epspow
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_abbrevd121h0
      implicit none
      integer(c_int), intent(in) :: ncut
      complex(ki_nin), dimension(0:3), intent(in) :: Q_ext
      complex(ki_nin), intent(in) :: mu2_ext
      complex(ki_nin), intent(out) :: numerator
      complex(ki) :: d121
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      Q(1:4)  =cmplx(real(+Q_ext(0:3),  ki_nin), aimag(+Q_ext(0:3)), ki)
      d121 = 0.0_ki
      d121 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d121, ki), aimag(d121), ki_nin)
   end subroutine numerator_ninja
   !------#] subroutine numerator_ninja:
!---#] numerator interfaces:
end module p9_part21part21_part25part25part21_d121h0l1
