module     p9_part21part21_part25part25part21_d109h0l1_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p9_part21part &
   ! &21_part25part25part21/helicity0d109h0l1_qp.f90
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
      use p9_part21part21_part25part25part21_abbrevd109h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc109(80)
      complex(ki) :: Qspvak1e5
      complex(ki) :: Qspvae5k1
      complex(ki) :: Qspvak2e5
      complex(ki) :: Qspvae5k2
      complex(ki) :: Qspvae1e5
      complex(ki) :: Qspvae5e1
      complex(ki) :: QspQ
      complex(ki) :: Qspk2
      complex(ki) :: Qspval4l3
      complex(ki) :: Qspval4k2
      complex(ki) :: Qspval3l4
      complex(ki) :: Qspvak2l4
      complex(ki) :: Qspval3k2
      complex(ki) :: Qspvak2l3
      complex(ki) :: Qspvae5l4
      complex(ki) :: Qspval4e5
      complex(ki) :: Qspvae1l3
      complex(ki) :: Qspval3e1
      complex(ki) :: Qspvae1k2
      complex(ki) :: Qspvak2e1
      complex(ki) :: Qspval3k1
      complex(ki) :: Qspvak2k1
      complex(ki) :: Qspvak1l3
      complex(ki) :: Qspvak1k2
      Qspvak1e5 = dotproduct(Q,spvak1e5)
      Qspvae5k1 = dotproduct(Q,spvae5k1)
      Qspvak2e5 = dotproduct(Q,spvak2e5)
      Qspvae5k2 = dotproduct(Q,spvae5k2)
      Qspvae1e5 = dotproduct(Q,spvae1e5)
      Qspvae5e1 = dotproduct(Q,spvae5e1)
      QspQ = dotproduct(Q,Q)
      Qspk2 = dotproduct(Q,k2)
      Qspval4l3 = dotproduct(Q,spval4l3)
      Qspval4k2 = dotproduct(Q,spval4k2)
      Qspval3l4 = dotproduct(Q,spval3l4)
      Qspvak2l4 = dotproduct(Q,spvak2l4)
      Qspval3k2 = dotproduct(Q,spval3k2)
      Qspvak2l3 = dotproduct(Q,spvak2l3)
      Qspvae5l4 = dotproduct(Q,spvae5l4)
      Qspval4e5 = dotproduct(Q,spval4e5)
      Qspvae1l3 = dotproduct(Q,spvae1l3)
      Qspval3e1 = dotproduct(Q,spval3e1)
      Qspvae1k2 = dotproduct(Q,spvae1k2)
      Qspvak2e1 = dotproduct(Q,spvak2e1)
      Qspval3k1 = dotproduct(Q,spval3k1)
      Qspvak2k1 = dotproduct(Q,spvak2k1)
      Qspvak1l3 = dotproduct(Q,spvak1l3)
      Qspvak1k2 = dotproduct(Q,spvak1k2)
      acc109(1)=abb109(7)
      acc109(2)=abb109(8)
      acc109(3)=abb109(9)
      acc109(4)=abb109(11)
      acc109(5)=abb109(12)
      acc109(6)=abb109(13)
      acc109(7)=abb109(14)
      acc109(8)=abb109(15)
      acc109(9)=abb109(16)
      acc109(10)=abb109(17)
      acc109(11)=abb109(18)
      acc109(12)=abb109(19)
      acc109(13)=abb109(20)
      acc109(14)=abb109(21)
      acc109(15)=abb109(22)
      acc109(16)=abb109(23)
      acc109(17)=abb109(24)
      acc109(18)=abb109(25)
      acc109(19)=abb109(26)
      acc109(20)=abb109(27)
      acc109(21)=abb109(29)
      acc109(22)=abb109(30)
      acc109(23)=abb109(31)
      acc109(24)=abb109(33)
      acc109(25)=abb109(37)
      acc109(26)=abb109(39)
      acc109(27)=abb109(42)
      acc109(28)=abb109(44)
      acc109(29)=abb109(45)
      acc109(30)=abb109(46)
      acc109(31)=abb109(47)
      acc109(32)=abb109(48)
      acc109(33)=abb109(49)
      acc109(34)=abb109(50)
      acc109(35)=abb109(51)
      acc109(36)=abb109(52)
      acc109(37)=abb109(54)
      acc109(38)=abb109(55)
      acc109(39)=abb109(56)
      acc109(40)=abb109(57)
      acc109(41)=abb109(58)
      acc109(42)=abb109(59)
      acc109(43)=abb109(60)
      acc109(44)=abb109(61)
      acc109(45)=abb109(63)
      acc109(46)=abb109(65)
      acc109(47)=abb109(67)
      acc109(48)=abb109(69)
      acc109(49)=abb109(71)
      acc109(50)=abb109(73)
      acc109(51)=abb109(79)
      acc109(52)=abb109(82)
      acc109(53)=abb109(83)
      acc109(54)=abb109(87)
      acc109(55)=abb109(88)
      acc109(56)=Qspvak1e5*Qspvae5k1
      acc109(57)=Qspvak2e5*Qspvae5k2
      acc109(56)=acc109(56)-acc109(57)
      acc109(57)=acc109(51)*acc109(56)
      acc109(58)=acc109(9)*Qspvae1e5
      acc109(59)=acc109(13)*Qspvae5e1
      acc109(60)=acc109(25)*Qspvak1e5
      acc109(61)=acc109(28)*Qspvae5k1
      acc109(62)=acc109(44)*Qspvae5k2
      acc109(63)=acc109(45)*Qspvak2e5
      acc109(57)=acc109(63)+acc109(62)+acc109(61)+acc109(60)+acc109(59)+acc109(&
      &58)+acc109(6)+acc109(57)
      acc109(57)=QspQ*acc109(57)
      acc109(58)=acc109(29)*Qspvak1e5
      acc109(59)=acc109(33)*Qspvae5k1
      acc109(60)=acc109(35)*Qspvae5e1
      acc109(61)=acc109(37)*Qspvae1e5
      acc109(58)=acc109(47)+acc109(61)+acc109(60)+acc109(59)+acc109(58)
      acc109(58)=Qspk2*acc109(58)
      acc109(59)=acc109(40)*Qspvak1e5
      acc109(60)=acc109(42)*Qspvae1e5
      acc109(61)=acc109(53)*Qspvak2e5
      acc109(59)=acc109(61)+acc109(60)+acc109(59)+acc109(36)
      acc109(59)=Qspval4l3*acc109(59)
      acc109(60)=acc109(22)*Qspvae1e5
      acc109(61)=acc109(54)*Qspvak1e5
      acc109(62)=acc109(55)*Qspvak2e5
      acc109(60)=acc109(62)+acc109(61)+acc109(60)+acc109(11)
      acc109(60)=Qspval4k2*acc109(60)
      acc109(61)=acc109(2)*Qspvae5e1
      acc109(62)=acc109(3)*Qspvae5k1
      acc109(63)=acc109(30)*Qspvae5k2
      acc109(61)=acc109(63)+acc109(62)+acc109(61)+acc109(1)
      acc109(61)=Qspval3l4*acc109(61)
      acc109(62)=acc109(15)*Qspvae5e1
      acc109(63)=acc109(38)*Qspvae5k1
      acc109(64)=acc109(43)*Qspvae5k2
      acc109(62)=acc109(64)+acc109(63)+acc109(62)+acc109(8)
      acc109(62)=Qspvak2l4*acc109(62)
      acc109(63)=acc109(18)*Qspvae5e1
      acc109(64)=acc109(23)*Qspvae5k1
      acc109(63)=acc109(64)+acc109(63)+acc109(10)
      acc109(63)=Qspval3k2*acc109(63)
      acc109(64)=acc109(5)*Qspvae1e5
      acc109(65)=acc109(7)*Qspvak1e5
      acc109(64)=acc109(65)+acc109(64)+acc109(4)
      acc109(64)=Qspvak2l3*acc109(64)
      acc109(56)=-acc109(48)*acc109(56)
      acc109(65)=acc109(31)*QspQ
      acc109(65)=acc109(65)+acc109(21)
      acc109(65)=acc109(65)*Qspvae5e1
      acc109(65)=acc109(16)+acc109(65)
      acc109(65)=Qspvae1e5*acc109(65)
      acc109(66)=acc109(17)*Qspvae5e1
      acc109(67)=acc109(41)*Qspvak1e5
      acc109(68)=acc109(46)*Qspvae5k2
      acc109(69)=acc109(50)*Qspvak2e5
      acc109(70)=acc109(52)*Qspvae5k1
      acc109(71)=Qspvae5l4*acc109(39)
      acc109(72)=Qspval4e5*acc109(32)
      acc109(73)=Qspvae1l3*acc109(34)
      acc109(74)=Qspval3e1*acc109(24)
      acc109(75)=Qspvae1k2*acc109(49)
      acc109(76)=Qspvak2e1*acc109(27)
      acc109(77)=Qspval3k1*acc109(26)
      acc109(78)=Qspvak2k1*acc109(12)
      acc109(79)=Qspvak1l3*acc109(19)
      acc109(80)=Qspvak1k2*acc109(14)
      brack=acc109(20)+acc109(56)+acc109(57)+acc109(58)+acc109(59)+acc109(60)+a&
      &cc109(61)+acc109(62)+acc109(63)+acc109(64)+acc109(65)+acc109(66)+acc109(&
      &67)+acc109(68)+acc109(69)+acc109(70)+acc109(71)+acc109(72)+acc109(73)+ac&
      &c109(74)+acc109(75)+acc109(76)+acc109(77)+acc109(78)+acc109(79)+acc109(8&
      &0)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ subroutine numerator_ninja:
   subroutine numerator_ninja(ncut, Q_ext, mu2_ext, numerator) &
   & bind(c, name="p9_part21part21_part25part25part21_d109h0l1_qp_ninja")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki_nin
      use p9_part21part21_part25part25part21_globalsl1_qp, only: epspow
      use p9_part21part21_part25part25part21_kinematics_qp
      use p9_part21part21_part25part25part21_abbrevd109h0_qp
      implicit none
      integer(c_int), intent(in) :: ncut
      complex(ki_nin), dimension(0:3), intent(in) :: Q_ext
      complex(ki_nin), intent(in) :: mu2_ext
      complex(ki_nin), intent(out) :: numerator
      complex(ki) :: d109
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = k3+k5
      Q(1:4)  =cmplx(real(+Q_ext(0:3)  -qshift(:),  ki_nin), aimag(+Q_ext(0:3))&
      &, ki)
      d109 = 0.0_ki
      d109 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d109, ki), aimag(d109), ki_nin)
   end subroutine numerator_ninja
   !------#] subroutine numerator_ninja:
!---#] numerator interfaces:
end module p9_part21part21_part25part25part21_d109h0l1_qp
