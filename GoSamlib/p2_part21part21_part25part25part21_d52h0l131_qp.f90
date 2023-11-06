module     p2_part21part21_part25part25part21_d52h0l131_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d52h0l131_qp.f90
   ! generator: buildfortran_tn3.py
   use p2_part21part21_part25part25part21_config, only: ki => ki_qp
   use p2_part21part21_part25part25part21_util_qp, only: cond_t, d => metric_te&
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
      use p2_part21part21_part25part25part21_model_qp
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_color_qp
      use p2_part21part21_part25part25part21_abbrevd52h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA, ninjaE3, ninjaE4
      complex(ki), intent(in) :: ninjaP
      complex(ki), dimension(9) :: acd52
      complex(ki), dimension (0:*), intent(inout) :: brack
      acd52(1)=dotproduct(ninjaE3,spvae1e2)
      acd52(2)=dotproduct(ninjaE3,spvae5e1)
      acd52(3)=dotproduct(ninjaE3,spvae2e5)
      acd52(4)=abb52(12)
      acd52(5)=dotproduct(ninjaE3,spvae2e1)
      acd52(6)=dotproduct(ninjaE3,spvae1e5)
      acd52(7)=dotproduct(ninjaE3,spvae5e2)
      acd52(8)=acd52(3)*acd52(2)*acd52(1)
      acd52(9)=acd52(7)*acd52(6)*acd52(5)
      acd52(8)=acd52(8)+acd52(9)
      acd52(8)=acd52(4)*acd52(8)
      brack(ninjaidxt3mu0)=0.0_ki
      brack(ninjaidxt2mu0)=acd52(8)
   end subroutine brack_31
!---#] subroutine brack_31:
!---#[ subroutine brack_32:
   pure subroutine brack_32(ninjaA, ninjaE3, ninjaE4, ninjaP, brack)
      use p2_part21part21_part25part25part21_model_qp
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_color_qp
      use p2_part21part21_part25part25part21_abbrevd52h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA, ninjaE3, ninjaE4
      complex(ki), intent(in) :: ninjaP
      complex(ki), dimension(103) :: acd52
      complex(ki), dimension (0:*), intent(inout) :: brack
      acd52(1)=dotproduct(ninjaE3,spvae1e2)
      acd52(2)=dotproduct(ninjaE3,spvae2e5)
      acd52(3)=dotproduct(ninjaE4,spvae5e1)
      acd52(4)=abb52(12)
      acd52(5)=dotproduct(ninjaE3,spvae5e1)
      acd52(6)=dotproduct(ninjaE4,spvae2e5)
      acd52(7)=abb52(16)
      acd52(8)=dotproduct(ninjaE4,spvae1e2)
      acd52(9)=dotproduct(ninjaE3,spvae2e1)
      acd52(10)=dotproduct(ninjaE3,spvae5e2)
      acd52(11)=dotproduct(ninjaE4,spvae1e5)
      acd52(12)=dotproduct(ninjaE3,spvae1e5)
      acd52(13)=dotproduct(ninjaE4,spvae5e2)
      acd52(14)=abb52(24)
      acd52(15)=dotproduct(ninjaE4,spvae2e1)
      acd52(16)=abb52(26)
      acd52(17)=abb52(30)
      acd52(18)=dotproduct(ninjaA,ninjaE3)
      acd52(19)=dotproduct(ninjaA,spvae1e2)
      acd52(20)=dotproduct(ninjaA,spvae2e5)
      acd52(21)=dotproduct(ninjaA,spvae2e1)
      acd52(22)=dotproduct(ninjaA,spvae5e2)
      acd52(23)=dotproduct(ninjaA,spvae5e1)
      acd52(24)=dotproduct(ninjaA,spvae1e5)
      acd52(25)=abb52(2)
      acd52(26)=dotproduct(ninjaE3,spvae2k1)
      acd52(27)=abb52(13)
      acd52(28)=abb52(18)
      acd52(29)=dotproduct(ninjaE3,spvae5k2)
      acd52(30)=abb52(37)
      acd52(31)=dotproduct(ninjaE3,spvae2k5)
      acd52(32)=abb52(33)
      acd52(33)=dotproduct(ninjaE3,spvak1e2)
      acd52(34)=abb52(19)
      acd52(35)=abb52(23)
      acd52(36)=abb52(5)
      acd52(37)=dotproduct(ninjaE3,spvak2e5)
      acd52(38)=abb52(15)
      acd52(39)=abb52(31)
      acd52(40)=dotproduct(ninjaE3,spvak5e2)
      acd52(41)=abb52(35)
      acd52(42)=abb52(28)
      acd52(43)=abb52(7)
      acd52(44)=abb52(21)
      acd52(45)=abb52(17)
      acd52(46)=abb52(27)
      acd52(47)=dotproduct(ninjaA,ninjaA)
      acd52(48)=abb52(4)
      acd52(49)=dotproduct(ninjaA,spvak1e2)
      acd52(50)=dotproduct(ninjaA,spvae2k1)
      acd52(51)=dotproduct(ninjaA,spvak2e5)
      acd52(52)=dotproduct(ninjaA,spvae5k2)
      acd52(53)=dotproduct(ninjaA,spvak5e2)
      acd52(54)=dotproduct(ninjaA,spvae2k5)
      acd52(55)=abb52(10)
      acd52(56)=abb52(11)
      acd52(57)=abb52(3)
      acd52(58)=abb52(9)
      acd52(59)=abb52(20)
      acd52(60)=abb52(6)
      acd52(61)=abb52(38)
      acd52(62)=abb52(22)
      acd52(63)=abb52(14)
      acd52(64)=abb52(36)
      acd52(65)=abb52(34)
      acd52(66)=abb52(32)
      acd52(67)=acd52(10)*acd52(4)
      acd52(68)=acd52(67)*acd52(11)
      acd52(69)=acd52(12)*acd52(4)
      acd52(70)=acd52(69)*acd52(13)
      acd52(68)=acd52(68)+acd52(70)
      acd52(68)=acd52(68)*acd52(9)
      acd52(70)=acd52(2)*acd52(4)
      acd52(71)=acd52(70)*acd52(3)
      acd52(72)=acd52(5)*acd52(4)
      acd52(73)=acd52(72)*acd52(6)
      acd52(71)=acd52(71)+acd52(73)
      acd52(71)=acd52(71)*acd52(1)
      acd52(73)=acd52(5)*acd52(70)*acd52(8)
      acd52(74)=acd52(12)*acd52(67)*acd52(15)
      acd52(68)=acd52(71)+acd52(68)+acd52(73)+acd52(74)
      acd52(71)=acd52(7)*acd52(1)
      acd52(73)=acd52(14)*acd52(9)
      acd52(74)=acd52(16)*acd52(5)
      acd52(75)=acd52(17)*acd52(12)
      acd52(71)=acd52(71)+acd52(73)+acd52(74)+acd52(75)
      acd52(73)=acd52(71)+acd52(68)
      acd52(74)=acd52(27)*acd52(26)
      acd52(75)=acd52(30)*acd52(29)
      acd52(76)=acd52(32)*acd52(31)
      acd52(74)=acd52(76)+acd52(74)+acd52(75)
      acd52(75)=acd52(23)*acd52(70)
      acd52(76)=2.0_ki*acd52(18)
      acd52(77)=acd52(7)*acd52(76)
      acd52(78)=acd52(20)*acd52(72)
      acd52(79)=acd52(25)*acd52(2)
      acd52(80)=acd52(28)*acd52(5)
      acd52(75)=acd52(79)+acd52(80)+acd52(78)+acd52(77)+acd52(75)+acd52(74)
      acd52(75)=acd52(1)*acd52(75)
      acd52(77)=acd52(34)*acd52(33)
      acd52(78)=acd52(38)*acd52(37)
      acd52(79)=acd52(41)*acd52(40)
      acd52(77)=acd52(79)+acd52(77)+acd52(78)
      acd52(78)=acd52(24)*acd52(67)
      acd52(79)=acd52(14)*acd52(76)
      acd52(80)=acd52(22)*acd52(69)
      acd52(81)=acd52(36)*acd52(10)
      acd52(82)=acd52(39)*acd52(12)
      acd52(78)=acd52(81)+acd52(82)+acd52(80)+acd52(79)+acd52(78)+acd52(77)
      acd52(78)=acd52(9)*acd52(78)
      acd52(79)=acd52(35)*acd52(33)
      acd52(80)=acd52(43)*acd52(37)
      acd52(81)=acd52(44)*acd52(40)
      acd52(70)=acd52(70)*acd52(19)
      acd52(70)=acd52(79)+acd52(80)+acd52(81)+acd52(70)
      acd52(79)=acd52(16)*acd52(76)
      acd52(79)=acd52(79)+acd52(70)
      acd52(79)=acd52(5)*acd52(79)
      acd52(80)=acd52(42)*acd52(26)
      acd52(81)=acd52(45)*acd52(29)
      acd52(82)=acd52(46)*acd52(31)
      acd52(67)=acd52(67)*acd52(21)
      acd52(67)=acd52(80)+acd52(81)+acd52(82)+acd52(67)
      acd52(80)=acd52(17)*acd52(76)
      acd52(80)=acd52(80)+acd52(67)
      acd52(80)=acd52(12)*acd52(80)
      acd52(75)=acd52(80)+acd52(79)+acd52(78)+acd52(75)
      acd52(70)=acd52(23)*acd52(70)
      acd52(67)=acd52(24)*acd52(67)
      acd52(68)=ninjaP*acd52(68)
      acd52(71)=acd52(47)*acd52(71)
      acd52(74)=acd52(19)*acd52(74)
      acd52(77)=acd52(21)*acd52(77)
      acd52(78)=ninjaP*acd52(1)
      acd52(79)=acd52(19)*acd52(76)
      acd52(78)=acd52(78)+acd52(79)
      acd52(78)=acd52(7)*acd52(78)
      acd52(79)=ninjaP*acd52(9)
      acd52(80)=acd52(21)*acd52(76)
      acd52(79)=acd52(79)+acd52(80)
      acd52(79)=acd52(14)*acd52(79)
      acd52(80)=ninjaP*acd52(5)
      acd52(81)=acd52(23)*acd52(76)
      acd52(80)=acd52(80)+acd52(81)
      acd52(80)=acd52(16)*acd52(80)
      acd52(81)=ninjaP*acd52(12)
      acd52(82)=acd52(24)*acd52(76)
      acd52(81)=acd52(81)+acd52(82)
      acd52(81)=acd52(17)*acd52(81)
      acd52(72)=acd52(19)*acd52(72)
      acd52(82)=acd52(23)*acd52(1)
      acd52(83)=acd52(4)*acd52(82)
      acd52(72)=acd52(72)+acd52(83)
      acd52(72)=acd52(20)*acd52(72)
      acd52(69)=acd52(21)*acd52(69)
      acd52(83)=acd52(24)*acd52(9)
      acd52(84)=acd52(4)*acd52(83)
      acd52(69)=acd52(69)+acd52(84)
      acd52(69)=acd52(22)*acd52(69)
      acd52(84)=acd52(19)*acd52(2)
      acd52(85)=acd52(20)*acd52(1)
      acd52(84)=acd52(84)+acd52(85)
      acd52(84)=acd52(25)*acd52(84)
      acd52(85)=acd52(19)*acd52(5)
      acd52(82)=acd52(85)+acd52(82)
      acd52(82)=acd52(28)*acd52(82)
      acd52(85)=acd52(21)*acd52(10)
      acd52(86)=acd52(22)*acd52(9)
      acd52(85)=acd52(85)+acd52(86)
      acd52(85)=acd52(36)*acd52(85)
      acd52(86)=acd52(21)*acd52(12)
      acd52(83)=acd52(86)+acd52(83)
      acd52(83)=acd52(39)*acd52(83)
      acd52(86)=acd52(34)*acd52(9)
      acd52(87)=acd52(35)*acd52(5)
      acd52(86)=acd52(86)+acd52(87)
      acd52(86)=acd52(49)*acd52(86)
      acd52(87)=acd52(27)*acd52(1)
      acd52(88)=acd52(42)*acd52(12)
      acd52(87)=acd52(87)+acd52(88)
      acd52(87)=acd52(50)*acd52(87)
      acd52(88)=acd52(38)*acd52(9)
      acd52(89)=acd52(43)*acd52(5)
      acd52(88)=acd52(88)+acd52(89)
      acd52(88)=acd52(51)*acd52(88)
      acd52(89)=acd52(30)*acd52(1)
      acd52(90)=acd52(45)*acd52(12)
      acd52(89)=acd52(89)+acd52(90)
      acd52(89)=acd52(52)*acd52(89)
      acd52(90)=acd52(41)*acd52(9)
      acd52(91)=acd52(44)*acd52(5)
      acd52(90)=acd52(90)+acd52(91)
      acd52(90)=acd52(53)*acd52(90)
      acd52(91)=acd52(32)*acd52(1)
      acd52(92)=acd52(46)*acd52(12)
      acd52(91)=acd52(91)+acd52(92)
      acd52(91)=acd52(54)*acd52(91)
      acd52(76)=acd52(48)*acd52(76)
      acd52(92)=acd52(55)*acd52(1)
      acd52(93)=acd52(56)*acd52(2)
      acd52(94)=acd52(57)*acd52(33)
      acd52(95)=acd52(58)*acd52(9)
      acd52(96)=acd52(59)*acd52(10)
      acd52(97)=acd52(60)*acd52(26)
      acd52(98)=acd52(61)*acd52(37)
      acd52(99)=acd52(62)*acd52(5)
      acd52(100)=acd52(63)*acd52(12)
      acd52(101)=acd52(64)*acd52(29)
      acd52(102)=acd52(65)*acd52(40)
      acd52(103)=acd52(66)*acd52(31)
      acd52(67)=acd52(103)+acd52(102)+acd52(101)+acd52(100)+acd52(99)+acd52(98)&
      &+acd52(97)+acd52(96)+acd52(95)+acd52(94)+acd52(93)+acd52(92)+acd52(76)+a&
      &cd52(91)+acd52(90)+acd52(89)+acd52(88)+acd52(87)+acd52(86)+acd52(83)+acd&
      &52(85)+acd52(82)+acd52(84)+acd52(71)+acd52(69)+acd52(72)+acd52(81)+acd52&
      &(80)+acd52(79)+acd52(78)+acd52(68)+acd52(67)+acd52(70)+acd52(77)+acd52(7&
      &4)
      brack(ninjaidxt1mu0)=acd52(75)
      brack(ninjaidxt1mu2)=0.0_ki
      brack(ninjaidxt0mu0)=acd52(67)
      brack(ninjaidxt0mu2)=acd52(73)
   end subroutine brack_32
!---#] subroutine brack_32:
!---#[ subroutine numerator_t3:
   subroutine numerator_t3(ncut, a, b, c, param, deg, coeffs) &
   & bind(c, name="p2_part21part21_part25part25part21_d52h0_qp_ninja_t3")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki => ki_nin
      use p2_part21part21_part25part25part21_globalsl1_qp, only: epspow
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_abbrevd52h0_qp
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
end module     p2_part21part21_part25part25part21_d52h0l131_qp
