module     p2_part21part21_part25part25part21_d117h0l1_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d117h0l1_qp.f90
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
      use p2_part21part21_part25part25part21_abbrevd117h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc117(77)
      complex(ki) :: Qspvak2k5
      complex(ki) :: Qspvak5k2
      complex(ki) :: Qspvak1k5
      complex(ki) :: Qspvak5k1
      complex(ki) :: Qspval4e5
      complex(ki) :: Qspvae5l4
      complex(ki) :: Qspvak5e1
      complex(ki) :: Qspvae1k5
      complex(ki) :: Qspvak5e5
      complex(ki) :: Qspvae5k5
      complex(ki) :: Qspvak2e5
      complex(ki) :: Qspvae5k2
      complex(ki) :: Qspval3e5
      complex(ki) :: Qspvae5l3
      complex(ki) :: QspQ
      complex(ki) :: Qspvak1e5
      complex(ki) :: Qspvae5k1
      complex(ki) :: Qspvae1e5
      complex(ki) :: Qspvae5e1
      Qspvak2k5 = dotproduct(Q,spvak2k5)
      Qspvak5k2 = dotproduct(Q,spvak5k2)
      Qspvak1k5 = dotproduct(Q,spvak1k5)
      Qspvak5k1 = dotproduct(Q,spvak5k1)
      Qspval4e5 = dotproduct(Q,spval4e5)
      Qspvae5l4 = dotproduct(Q,spvae5l4)
      Qspvak5e1 = dotproduct(Q,spvak5e1)
      Qspvae1k5 = dotproduct(Q,spvae1k5)
      Qspvak5e5 = dotproduct(Q,spvak5e5)
      Qspvae5k5 = dotproduct(Q,spvae5k5)
      Qspvak2e5 = dotproduct(Q,spvak2e5)
      Qspvae5k2 = dotproduct(Q,spvae5k2)
      Qspval3e5 = dotproduct(Q,spval3e5)
      Qspvae5l3 = dotproduct(Q,spvae5l3)
      QspQ = dotproduct(Q,Q)
      Qspvak1e5 = dotproduct(Q,spvak1e5)
      Qspvae5k1 = dotproduct(Q,spvae5k1)
      Qspvae1e5 = dotproduct(Q,spvae1e5)
      Qspvae5e1 = dotproduct(Q,spvae5e1)
      acc117(1)=abb117(7)
      acc117(2)=abb117(8)
      acc117(3)=abb117(9)
      acc117(4)=abb117(10)
      acc117(5)=abb117(11)
      acc117(6)=abb117(12)
      acc117(7)=abb117(13)
      acc117(8)=abb117(15)
      acc117(9)=abb117(16)
      acc117(10)=abb117(17)
      acc117(11)=abb117(18)
      acc117(12)=abb117(19)
      acc117(13)=abb117(20)
      acc117(14)=abb117(21)
      acc117(15)=abb117(22)
      acc117(16)=abb117(23)
      acc117(17)=abb117(24)
      acc117(18)=abb117(25)
      acc117(19)=abb117(26)
      acc117(20)=abb117(27)
      acc117(21)=abb117(28)
      acc117(22)=abb117(29)
      acc117(23)=abb117(30)
      acc117(24)=abb117(31)
      acc117(25)=abb117(33)
      acc117(26)=abb117(34)
      acc117(27)=abb117(35)
      acc117(28)=abb117(36)
      acc117(29)=abb117(37)
      acc117(30)=abb117(38)
      acc117(31)=abb117(39)
      acc117(32)=abb117(40)
      acc117(33)=abb117(41)
      acc117(34)=abb117(42)
      acc117(35)=abb117(43)
      acc117(36)=abb117(44)
      acc117(37)=abb117(46)
      acc117(38)=abb117(48)
      acc117(39)=abb117(49)
      acc117(40)=abb117(51)
      acc117(41)=abb117(52)
      acc117(42)=abb117(54)
      acc117(43)=abb117(56)
      acc117(44)=abb117(57)
      acc117(45)=abb117(60)
      acc117(46)=abb117(65)
      acc117(47)=abb117(66)
      acc117(48)=abb117(68)
      acc117(49)=abb117(69)
      acc117(50)=abb117(70)
      acc117(51)=abb117(71)
      acc117(52)=abb117(72)
      acc117(53)=abb117(73)
      acc117(54)=abb117(74)
      acc117(55)=abb117(76)
      acc117(56)=abb117(77)
      acc117(57)=abb117(78)
      acc117(58)=abb117(80)
      acc117(59)=Qspvak2k5*acc117(30)
      acc117(60)=Qspvak5k2*acc117(21)
      acc117(61)=Qspvak1k5*acc117(26)
      acc117(62)=Qspvak5k1*acc117(33)
      acc117(63)=Qspval4e5*acc117(11)
      acc117(64)=Qspvae5l4*acc117(36)
      acc117(65)=Qspvak5e1*acc117(52)
      acc117(66)=Qspvae1k5*acc117(49)
      acc117(67)=Qspvak5e5*acc117(47)
      acc117(68)=Qspvae5k5*acc117(43)
      acc117(69)=Qspvak2e5*acc117(4)
      acc117(70)=Qspvae5k2*acc117(8)
      acc117(71)=Qspval3e5*acc117(7)
      acc117(72)=Qspvae5l3*acc117(45)
      acc117(73)=QspQ*acc117(27)
      acc117(59)=acc117(73)+acc117(72)+acc117(71)+acc117(70)+acc117(69)+acc117(&
      &68)+acc117(67)+acc117(66)+acc117(65)+acc117(64)+acc117(63)+acc117(62)+ac&
      &c117(61)+acc117(60)+acc117(15)+acc117(59)
      acc117(59)=QspQ*acc117(59)
      acc117(60)=Qspvak2k5*acc117(32)
      acc117(61)=Qspvak1k5*acc117(34)
      acc117(62)=Qspvae5l4*acc117(57)
      acc117(63)=Qspvae1k5*acc117(58)
      acc117(64)=Qspvae5k5*acc117(53)
      acc117(65)=Qspvae5k2*acc117(29)
      acc117(60)=acc117(65)+acc117(64)+acc117(63)+acc117(62)+acc117(61)+acc117(&
      &1)+acc117(60)
      acc117(60)=Qspval3e5*acc117(60)
      acc117(61)=Qspvak5k2*acc117(22)
      acc117(62)=Qspvak5k1*acc117(28)
      acc117(63)=Qspval4e5*acc117(46)
      acc117(64)=Qspvak5e1*acc117(55)
      acc117(65)=Qspvak5e5*acc117(51)
      acc117(66)=Qspvak2e5*acc117(10)
      acc117(61)=acc117(66)+acc117(65)+acc117(64)+acc117(63)+acc117(62)+acc117(&
      &40)+acc117(61)
      acc117(61)=Qspvae5l3*acc117(61)
      acc117(62)=Qspvak5k1*acc117(31)
      acc117(63)=Qspval4e5*acc117(24)
      acc117(64)=Qspvak5e1*acc117(25)
      acc117(65)=Qspvak5e5*acc117(17)
      acc117(66)=Qspvak2e5*acc117(2)
      acc117(62)=acc117(66)+acc117(65)+acc117(64)+acc117(63)+acc117(5)+acc117(6&
      &2)
      acc117(62)=Qspvae5k2*acc117(62)
      acc117(63)=Qspvak1k5*acc117(23)
      acc117(64)=Qspvae5l4*acc117(6)
      acc117(65)=Qspvae1k5*acc117(19)
      acc117(66)=Qspvae5k5*acc117(9)
      acc117(63)=acc117(66)+acc117(65)+acc117(64)+acc117(3)+acc117(63)
      acc117(63)=Qspvak2e5*acc117(63)
      acc117(64)=Qspvak1e5*acc117(20)
      acc117(65)=-Qspvak1e5*acc117(12)
      acc117(65)=acc117(16)+acc117(65)
      acc117(65)=Qspvae5k1*acc117(65)
      acc117(66)=Qspvae1e5*acc117(35)
      acc117(67)=Qspvae1e5*acc117(39)
      acc117(67)=acc117(38)+acc117(67)
      acc117(67)=Qspvae5e1*acc117(67)
      acc117(68)=Qspvak2k5*acc117(18)
      acc117(69)=Qspvak5k2*acc117(13)
      acc117(70)=Qspvak1k5*acc117(56)
      acc117(71)=Qspvak5k1*acc117(54)
      acc117(72)=Qspval4e5*acc117(42)
      acc117(73)=Qspvae5l4*acc117(37)
      acc117(74)=Qspvak5e1*acc117(50)
      acc117(75)=Qspvae1k5*acc117(48)
      acc117(76)=Qspvak5e5*acc117(44)
      acc117(77)=Qspvae5k5*acc117(41)
      brack=acc117(14)+acc117(59)+acc117(60)+acc117(61)+acc117(62)+acc117(63)+a&
      &cc117(64)+acc117(65)+acc117(66)+acc117(67)+acc117(68)+acc117(69)+acc117(&
      &70)+acc117(71)+acc117(72)+acc117(73)+acc117(74)+acc117(75)+acc117(76)+ac&
      &c117(77)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ subroutine numerator_ninja:
   subroutine numerator_ninja(ncut, Q_ext, mu2_ext, numerator) &
   & bind(c, name="p2_part21part21_part25part25part21_d117h0l1_qp_ninja")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki_nin
      use p2_part21part21_part25part25part21_globalsl1_qp, only: epspow
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_abbrevd117h0_qp
      implicit none
      integer(c_int), intent(in) :: ncut
      complex(ki_nin), dimension(0:3), intent(in) :: Q_ext
      complex(ki_nin), intent(in) :: mu2_ext
      complex(ki_nin), intent(out) :: numerator
      complex(ki) :: d117
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = k5+k4
      Q(1:4)  =cmplx(real(+Q_ext(0:3)  -qshift(:),  ki_nin), aimag(+Q_ext(0:3))&
      &, ki)
      d117 = 0.0_ki
      d117 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d117, ki), aimag(d117), ki_nin)
   end subroutine numerator_ninja
   !------#] subroutine numerator_ninja:
!---#] numerator interfaces:
end module p2_part21part21_part25part25part21_d117h0l1_qp
