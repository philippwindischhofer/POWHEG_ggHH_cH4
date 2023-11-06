module     p9_part21part21_part25part25part21_d145h0l131_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p9_part21part &
   ! &21_part25part25part21/helicity0d145h0l131_qp.f90
   ! generator: buildfortran_tn3.py
   use p9_part21part21_part25part25part21_config, only: ki => ki_qp
   use p9_part21part21_part25part25part21_util_qp, only: cond_t, d => metric_te&
   &nsor
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   integer, parameter :: ninjaidxt3mu0 = 0
   integer, parameter :: ninjaidxt2mu0 = 1
   integer, parameter :: ninjaidxt1mu0 = 2
   integer, parameter :: ninjaidxt1mu2 = 3
   integer, parameter :: ninjaidxt0mu0 = 4
   integer, parameter :: ninjaidxt0mu2 = 5
   public :: numerator_t3
contains
!---#[ subroutine brack_31:
   pure subroutine brack_31(ninjaA, ninjaE3, ninjaE4, ninjaP, brack)
      use p9_part21part21_part25part25part21_model_qp
      use p9_part21part21_part25part25part21_kinematics_qp
      use p9_part21part21_part25part25part21_color_qp
      use p9_part21part21_part25part25part21_abbrevd145h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA, ninjaE3, ninjaE4
      complex(ki), intent(in) :: ninjaP
      complex(ki), dimension(9) :: acd145
      complex(ki), dimension (0:*), intent(inout) :: brack
      acd145(1)=dotproduct(ninjaE3,spvae1e2)
      acd145(2)=dotproduct(ninjaE3,spvae5e1)
      acd145(3)=dotproduct(ninjaE3,spvae2e5)
      acd145(4)=abb145(7)
      acd145(5)=dotproduct(ninjaE3,spvae2e1)
      acd145(6)=dotproduct(ninjaE3,spvae1e5)
      acd145(7)=dotproduct(ninjaE3,spvae5e2)
      acd145(8)=-acd145(3)*acd145(2)*acd145(1)
      acd145(9)=-acd145(7)*acd145(6)*acd145(5)
      acd145(8)=acd145(8)+acd145(9)
      acd145(8)=acd145(4)*acd145(8)
      brack(ninjaidxt3mu0)=0.0_ki
      brack(ninjaidxt2mu0)=acd145(8)
   end subroutine brack_31
!---#] subroutine brack_31:
!---#[ subroutine brack_32:
   pure subroutine brack_32(ninjaA, ninjaE3, ninjaE4, ninjaP, brack)
      use p9_part21part21_part25part25part21_model_qp
      use p9_part21part21_part25part25part21_kinematics_qp
      use p9_part21part21_part25part25part21_color_qp
      use p9_part21part21_part25part25part21_abbrevd145h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA, ninjaE3, ninjaE4
      complex(ki), intent(in) :: ninjaP
      complex(ki), dimension(103) :: acd145
      complex(ki), dimension (0:*), intent(inout) :: brack
      acd145(1)=dotproduct(ninjaE3,spvae2e1)
      acd145(2)=dotproduct(ninjaE3,spvae1e5)
      acd145(3)=dotproduct(ninjaE4,spvae5e2)
      acd145(4)=abb145(7)
      acd145(5)=dotproduct(ninjaE3,spvae5e2)
      acd145(6)=dotproduct(ninjaE4,spvae1e5)
      acd145(7)=abb145(19)
      acd145(8)=dotproduct(ninjaE3,spvae1e2)
      acd145(9)=dotproduct(ninjaE3,spvae5e1)
      acd145(10)=dotproduct(ninjaE4,spvae2e5)
      acd145(11)=dotproduct(ninjaE3,spvae2e5)
      acd145(12)=dotproduct(ninjaE4,spvae5e1)
      acd145(13)=abb145(25)
      acd145(14)=dotproduct(ninjaE4,spvae1e2)
      acd145(15)=abb145(22)
      acd145(16)=dotproduct(ninjaE4,spvae2e1)
      acd145(17)=abb145(27)
      acd145(18)=dotproduct(ninjaA,ninjaE3)
      acd145(19)=dotproduct(ninjaA,spvae2e1)
      acd145(20)=dotproduct(ninjaA,spvae1e2)
      acd145(21)=dotproduct(ninjaA,spvae5e1)
      acd145(22)=dotproduct(ninjaA,spvae2e5)
      acd145(23)=dotproduct(ninjaA,spvae1e5)
      acd145(24)=dotproduct(ninjaA,spvae5e2)
      acd145(25)=dotproduct(ninjaE3,spvak1e2)
      acd145(26)=abb145(4)
      acd145(27)=abb145(16)
      acd145(28)=abb145(15)
      acd145(29)=dotproduct(ninjaE3,spvak5e2)
      acd145(30)=abb145(34)
      acd145(31)=dotproduct(ninjaE3,spvak2e5)
      acd145(32)=abb145(37)
      acd145(33)=abb145(21)
      acd145(34)=abb145(6)
      acd145(35)=abb145(23)
      acd145(36)=dotproduct(ninjaE3,spvae5k2)
      acd145(37)=abb145(20)
      acd145(38)=dotproduct(ninjaE3,spvae2k1)
      acd145(39)=abb145(24)
      acd145(40)=dotproduct(ninjaE3,spvae2k5)
      acd145(41)=abb145(31)
      acd145(42)=abb145(33)
      acd145(43)=abb145(36)
      acd145(44)=abb145(8)
      acd145(45)=abb145(26)
      acd145(46)=abb145(30)
      acd145(47)=dotproduct(ninjaA,ninjaA)
      acd145(48)=abb145(13)
      acd145(49)=dotproduct(ninjaA,spvak1e2)
      acd145(50)=dotproduct(ninjaA,spvae5k2)
      acd145(51)=dotproduct(ninjaA,spvae2k1)
      acd145(52)=dotproduct(ninjaA,spvae2k5)
      acd145(53)=dotproduct(ninjaA,spvak5e2)
      acd145(54)=dotproduct(ninjaA,spvak2e5)
      acd145(55)=abb145(3)
      acd145(56)=abb145(12)
      acd145(57)=abb145(11)
      acd145(58)=abb145(18)
      acd145(59)=abb145(17)
      acd145(60)=abb145(9)
      acd145(61)=abb145(14)
      acd145(62)=abb145(10)
      acd145(63)=abb145(28)
      acd145(64)=abb145(29)
      acd145(65)=abb145(32)
      acd145(66)=abb145(35)
      acd145(67)=acd145(11)*acd145(4)
      acd145(68)=acd145(67)*acd145(12)
      acd145(69)=acd145(9)*acd145(4)
      acd145(70)=acd145(69)*acd145(10)
      acd145(68)=acd145(68)+acd145(70)
      acd145(68)=acd145(68)*acd145(8)
      acd145(70)=acd145(5)*acd145(4)
      acd145(71)=acd145(70)*acd145(6)
      acd145(72)=acd145(2)*acd145(4)
      acd145(73)=acd145(72)*acd145(3)
      acd145(71)=acd145(71)+acd145(73)
      acd145(71)=acd145(71)*acd145(1)
      acd145(73)=acd145(2)*acd145(70)*acd145(16)
      acd145(74)=acd145(9)*acd145(67)*acd145(14)
      acd145(68)=acd145(71)+acd145(68)+acd145(73)+acd145(74)
      acd145(71)=acd145(7)*acd145(1)
      acd145(73)=acd145(13)*acd145(8)
      acd145(74)=acd145(15)*acd145(9)
      acd145(75)=acd145(17)*acd145(2)
      acd145(71)=acd145(71)+acd145(73)+acd145(74)+acd145(75)
      acd145(73)=acd145(71)-acd145(68)
      acd145(74)=acd145(26)*acd145(25)
      acd145(75)=acd145(30)*acd145(29)
      acd145(76)=acd145(32)*acd145(31)
      acd145(74)=acd145(76)+acd145(74)+acd145(75)
      acd145(75)=-acd145(23)*acd145(70)
      acd145(76)=2.0_ki*acd145(18)
      acd145(77)=acd145(7)*acd145(76)
      acd145(78)=-acd145(24)*acd145(72)
      acd145(79)=acd145(27)*acd145(2)
      acd145(80)=-acd145(28)*acd145(5)
      acd145(75)=acd145(80)+acd145(79)+acd145(78)+acd145(77)+acd145(75)+acd145(&
      &74)
      acd145(75)=acd145(1)*acd145(75)
      acd145(77)=acd145(37)*acd145(36)
      acd145(78)=acd145(39)*acd145(38)
      acd145(79)=acd145(41)*acd145(40)
      acd145(77)=acd145(79)+acd145(77)+acd145(78)
      acd145(78)=-acd145(21)*acd145(67)
      acd145(79)=acd145(13)*acd145(76)
      acd145(80)=-acd145(22)*acd145(69)
      acd145(81)=acd145(34)*acd145(9)
      acd145(82)=-acd145(35)*acd145(11)
      acd145(78)=acd145(82)+acd145(81)+acd145(80)+acd145(79)+acd145(78)+acd145(&
      &77)
      acd145(78)=acd145(8)*acd145(78)
      acd145(79)=acd145(44)*acd145(36)
      acd145(80)=acd145(45)*acd145(38)
      acd145(81)=acd145(46)*acd145(40)
      acd145(70)=acd145(70)*acd145(19)
      acd145(70)=-acd145(79)-acd145(80)+acd145(81)+acd145(70)
      acd145(79)=acd145(17)*acd145(76)
      acd145(79)=acd145(79)-acd145(70)
      acd145(79)=acd145(2)*acd145(79)
      acd145(80)=acd145(33)*acd145(25)
      acd145(81)=acd145(42)*acd145(29)
      acd145(82)=acd145(43)*acd145(31)
      acd145(67)=acd145(67)*acd145(20)
      acd145(67)=acd145(80)+acd145(81)+acd145(82)-acd145(67)
      acd145(80)=acd145(15)*acd145(76)
      acd145(80)=acd145(80)+acd145(67)
      acd145(80)=acd145(9)*acd145(80)
      acd145(75)=acd145(80)+acd145(79)+acd145(78)+acd145(75)
      acd145(67)=acd145(21)*acd145(67)
      acd145(70)=-acd145(23)*acd145(70)
      acd145(68)=-ninjaP*acd145(68)
      acd145(71)=acd145(47)*acd145(71)
      acd145(74)=acd145(19)*acd145(74)
      acd145(77)=acd145(20)*acd145(77)
      acd145(78)=ninjaP*acd145(1)
      acd145(79)=acd145(19)*acd145(76)
      acd145(78)=acd145(78)+acd145(79)
      acd145(78)=acd145(7)*acd145(78)
      acd145(79)=ninjaP*acd145(8)
      acd145(80)=acd145(20)*acd145(76)
      acd145(79)=acd145(79)+acd145(80)
      acd145(79)=acd145(13)*acd145(79)
      acd145(80)=ninjaP*acd145(9)
      acd145(81)=acd145(21)*acd145(76)
      acd145(80)=acd145(80)+acd145(81)
      acd145(80)=acd145(15)*acd145(80)
      acd145(81)=ninjaP*acd145(2)
      acd145(82)=acd145(23)*acd145(76)
      acd145(81)=acd145(81)+acd145(82)
      acd145(81)=acd145(17)*acd145(81)
      acd145(69)=-acd145(20)*acd145(69)
      acd145(82)=acd145(21)*acd145(8)
      acd145(83)=-acd145(4)*acd145(82)
      acd145(69)=acd145(69)+acd145(83)
      acd145(69)=acd145(22)*acd145(69)
      acd145(72)=-acd145(19)*acd145(72)
      acd145(83)=acd145(23)*acd145(1)
      acd145(84)=-acd145(4)*acd145(83)
      acd145(72)=acd145(72)+acd145(84)
      acd145(72)=acd145(24)*acd145(72)
      acd145(84)=acd145(19)*acd145(2)
      acd145(83)=acd145(84)+acd145(83)
      acd145(83)=acd145(27)*acd145(83)
      acd145(84)=-acd145(19)*acd145(5)
      acd145(85)=-acd145(24)*acd145(1)
      acd145(84)=acd145(84)+acd145(85)
      acd145(84)=acd145(28)*acd145(84)
      acd145(85)=acd145(20)*acd145(9)
      acd145(82)=acd145(85)+acd145(82)
      acd145(82)=acd145(34)*acd145(82)
      acd145(85)=-acd145(20)*acd145(11)
      acd145(86)=-acd145(22)*acd145(8)
      acd145(85)=acd145(85)+acd145(86)
      acd145(85)=acd145(35)*acd145(85)
      acd145(86)=acd145(26)*acd145(1)
      acd145(87)=acd145(33)*acd145(9)
      acd145(86)=acd145(86)+acd145(87)
      acd145(86)=acd145(49)*acd145(86)
      acd145(87)=acd145(37)*acd145(8)
      acd145(88)=acd145(44)*acd145(2)
      acd145(87)=acd145(87)+acd145(88)
      acd145(87)=acd145(50)*acd145(87)
      acd145(88)=acd145(39)*acd145(8)
      acd145(89)=acd145(45)*acd145(2)
      acd145(88)=acd145(88)+acd145(89)
      acd145(88)=acd145(51)*acd145(88)
      acd145(89)=acd145(41)*acd145(8)
      acd145(90)=-acd145(46)*acd145(2)
      acd145(89)=acd145(89)+acd145(90)
      acd145(89)=acd145(52)*acd145(89)
      acd145(90)=acd145(30)*acd145(1)
      acd145(91)=acd145(42)*acd145(9)
      acd145(90)=acd145(90)+acd145(91)
      acd145(90)=acd145(53)*acd145(90)
      acd145(91)=acd145(32)*acd145(1)
      acd145(92)=acd145(43)*acd145(9)
      acd145(91)=acd145(91)+acd145(92)
      acd145(91)=acd145(54)*acd145(91)
      acd145(76)=acd145(48)*acd145(76)
      acd145(92)=acd145(55)*acd145(1)
      acd145(93)=acd145(56)*acd145(25)
      acd145(94)=acd145(57)*acd145(8)
      acd145(95)=acd145(58)*acd145(9)
      acd145(96)=acd145(59)*acd145(11)
      acd145(97)=acd145(60)*acd145(2)
      acd145(98)=acd145(61)*acd145(5)
      acd145(99)=acd145(62)*acd145(36)
      acd145(100)=acd145(63)*acd145(38)
      acd145(101)=acd145(64)*acd145(40)
      acd145(102)=acd145(65)*acd145(29)
      acd145(103)=acd145(66)*acd145(31)
      acd145(67)=acd145(103)+acd145(102)+acd145(101)+acd145(100)+acd145(99)+acd&
      &145(98)+acd145(97)+acd145(96)+acd145(95)+acd145(94)+acd145(93)+acd145(92&
      &)+acd145(76)+acd145(91)+acd145(90)+acd145(89)+acd145(88)+acd145(87)+acd1&
      &45(86)+acd145(85)+acd145(82)+acd145(84)+acd145(83)+acd145(71)+acd145(72)&
      &+acd145(69)+acd145(81)+acd145(80)+acd145(79)+acd145(78)+acd145(68)+acd14&
      &5(70)+acd145(67)+acd145(77)+acd145(74)
      brack(ninjaidxt1mu0)=acd145(75)
      brack(ninjaidxt1mu2)=0.0_ki
      brack(ninjaidxt0mu0)=acd145(67)
      brack(ninjaidxt0mu2)=acd145(73)
   end subroutine brack_32
!---#] subroutine brack_32:
!---#[ subroutine numerator_t3:
   subroutine numerator_t3(ncut, a, b, c, param, deg, coeffs) &
   & bind(c, name="p9_part21part21_part25part25part21_d145h0_qp_ninja_t3")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki => ki_nin
      use p9_part21part21_part25part25part21_globalsl1_qp, only: epspow
      use p9_part21part21_part25part25part21_kinematics_qp
      use p9_part21part21_part25part25part21_abbrevd145h0_qp
      implicit none
      integer(c_int), intent(in) :: ncut, deg
      complex(ki), dimension(0:3), intent(in) :: a, b, c
      complex(ki), intent(in) :: param
      complex(ki), dimension(4) :: vecA, vecB, vecC
      complex(ki), dimension(0:*), intent(out) :: coeffs
      integer :: t1
      complex(ki), dimension(4) :: qshift
      qshift = k2
	     vecA(1:4) = + a(0:3) - qshift(1:4)
	     vecB(1:4) = + b(0:3)
	     vecC(1:4) = + c(0:3)
      if (deg.lt.0) return
      t1 = 0
      call cond_t(epspow.eq.t1,brack_31,vecA,vecB,vecC,param,coeffs)
      if (deg.le.(1+(0))) return
      call cond_t(epspow.eq.t1,brack_32,vecA,vecB,vecC,param,coeffs)
   end subroutine numerator_t3
!---#] subroutine numerator_t3:
end module     p9_part21part21_part25part25part21_d145h0l131_qp
