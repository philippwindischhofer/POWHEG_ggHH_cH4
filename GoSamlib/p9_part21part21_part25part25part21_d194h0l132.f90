module     p9_part21part21_part25part25part21_d194h0l132
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p9_part21part &
   ! &21_part25part25part21/helicity0d194h0l132.f90
   ! generator: buildfortran_tn3.py
   use p9_part21part21_part25part25part21_config, only: ki
   use p9_part21part21_part25part25part21_util, only: cond_t, d => metric_tensor
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   integer, parameter :: ninjaidxt2x0mu0 = 0
   integer, parameter :: ninjaidxt1x0mu0 = 1
   integer, parameter :: ninjaidxt1x1mu0 = 2
   integer, parameter :: ninjaidxt0x0mu0 = 3
   integer, parameter :: ninjaidxt0x0mu2 = 4
   integer, parameter :: ninjaidxt0x1mu0 = 5
   integer, parameter :: ninjaidxt0x2mu0 = 6
   public :: numerator_t2
contains
!---#[ subroutine brack_21:
   pure subroutine brack_21(ninjaA0, ninjaA1, ninjaE3, ninjaE4, ninjaP, brack)
      use p9_part21part21_part25part25part21_model
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_color
      use p9_part21part21_part25part25part21_abbrevd194h0
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(21) :: acd194
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd194(1)=dotproduct(ninjaA0,ninjaE3)
      acd194(2)=dotproduct(ninjaE3,spvae1e2)
      acd194(3)=dotproduct(ninjaE3,spvae5e1)
      acd194(4)=dotproduct(ninjaE3,spvae2e5)
      acd194(5)=abb194(31)
      acd194(6)=dotproduct(ninjaE3,spvae2e1)
      acd194(7)=dotproduct(ninjaE3,spvae1e5)
      acd194(8)=dotproduct(ninjaE3,spvae5e2)
      acd194(9)=dotproduct(ninjaE3,spvae2l4)
      acd194(10)=dotproduct(ninjaE3,spval4e5)
      acd194(11)=dotproduct(ninjaE3,spvae2k2)
      acd194(12)=dotproduct(ninjaE3,spvak2e5)
      acd194(13)=abb194(86)
      acd194(14)=dotproduct(ninjaE3,spval4e2)
      acd194(15)=dotproduct(ninjaE3,spvae5l4)
      acd194(16)=dotproduct(ninjaE3,spvak2e2)
      acd194(17)=dotproduct(ninjaE3,spvae5k2)
      acd194(18)=-acd194(14)*acd194(15)
      acd194(19)=2.0_ki*acd194(1)
      acd194(20)=acd194(8)*acd194(19)
      acd194(18)=acd194(18)+acd194(20)
      acd194(20)=acd194(6)*acd194(7)
      acd194(18)=acd194(18)*acd194(20)
      acd194(21)=-acd194(9)*acd194(10)
      acd194(19)=acd194(4)*acd194(19)
      acd194(19)=acd194(21)+acd194(19)
      acd194(21)=acd194(2)*acd194(3)
      acd194(19)=acd194(19)*acd194(21)
      acd194(18)=acd194(18)+acd194(19)
      acd194(18)=acd194(5)*acd194(18)
      acd194(19)=-acd194(16)*acd194(17)*acd194(20)
      acd194(20)=-acd194(11)*acd194(12)*acd194(21)
      acd194(19)=acd194(19)+acd194(20)
      acd194(19)=acd194(13)*acd194(19)
      acd194(18)=acd194(18)+acd194(19)
      brack(ninjaidxt2x0mu0)=0.0_ki
      brack(ninjaidxt1x0mu0)=acd194(18)
      brack(ninjaidxt1x1mu0)=0.0_ki
   end subroutine brack_21
!---#] subroutine brack_21:
!---#[ subroutine brack_22:
   pure subroutine brack_22(ninjaA0, ninjaA1, ninjaE3, ninjaE4, ninjaP, brack)
      use p9_part21part21_part25part25part21_model
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_color
      use p9_part21part21_part25part25part21_abbrevd194h0
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(92) :: acd194
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd194(1)=dotproduct(ninjaE3,spvae2e1)
      acd194(2)=dotproduct(ninjaE3,spvae1e5)
      acd194(3)=dotproduct(ninjaE3,spvae5e2)
      acd194(4)=abb194(31)
      acd194(5)=dotproduct(ninjaE3,spvae1e2)
      acd194(6)=dotproduct(ninjaE3,spvae5e1)
      acd194(7)=dotproduct(ninjaE3,spvae2e5)
      acd194(8)=dotproduct(ninjaA0,ninjaA1)
      acd194(9)=dotproduct(ninjaA0,ninjaE3)
      acd194(10)=dotproduct(ninjaA1,spvae2e1)
      acd194(11)=dotproduct(ninjaA1,spvae1e5)
      acd194(12)=dotproduct(ninjaA1,spvae1e2)
      acd194(13)=dotproduct(ninjaA1,spvae5e1)
      acd194(14)=dotproduct(ninjaA1,spvae2e5)
      acd194(15)=dotproduct(ninjaA1,spvae5e2)
      acd194(16)=dotproduct(ninjaA1,spvak2e2)
      acd194(17)=dotproduct(ninjaE3,spvae5k2)
      acd194(18)=abb194(86)
      acd194(19)=dotproduct(ninjaE3,spvak2e2)
      acd194(20)=dotproduct(ninjaE3,spval4e2)
      acd194(21)=dotproduct(ninjaE3,spvae5l4)
      acd194(22)=dotproduct(ninjaA1,spvae5k2)
      acd194(23)=dotproduct(ninjaA1,spval4e2)
      acd194(24)=dotproduct(ninjaA1,spvae2k2)
      acd194(25)=dotproduct(ninjaE3,spvak2e5)
      acd194(26)=dotproduct(ninjaE3,spvae2k2)
      acd194(27)=dotproduct(ninjaE3,spvae2l4)
      acd194(28)=dotproduct(ninjaE3,spval4e5)
      acd194(29)=dotproduct(ninjaA1,spvae2l4)
      acd194(30)=dotproduct(ninjaA1,spval4e5)
      acd194(31)=dotproduct(ninjaA1,spvae5l4)
      acd194(32)=dotproduct(ninjaA1,spvak2e5)
      acd194(33)=dotproduct(ninjaA1,ninjaA1)
      acd194(34)=dotproduct(ninjaA0,ninjaA0)
      acd194(35)=dotproduct(ninjaA0,spvae2e1)
      acd194(36)=dotproduct(ninjaA0,spvae1e5)
      acd194(37)=dotproduct(ninjaA0,spvae1e2)
      acd194(38)=dotproduct(ninjaA0,spvae5e1)
      acd194(39)=dotproduct(ninjaA0,spvae2e5)
      acd194(40)=dotproduct(ninjaA0,spvae5e2)
      acd194(41)=abb194(38)
      acd194(42)=abb194(63)
      acd194(43)=abb194(27)
      acd194(44)=abb194(62)
      acd194(45)=dotproduct(ninjaA0,spvak2e2)
      acd194(46)=dotproduct(ninjaA0,spvae5k2)
      acd194(47)=dotproduct(ninjaA0,spval4e2)
      acd194(48)=dotproduct(ninjaA0,spvae2k2)
      acd194(49)=dotproduct(ninjaA0,spvae2l4)
      acd194(50)=dotproduct(ninjaA0,spval4e5)
      acd194(51)=dotproduct(ninjaA0,spvae5l4)
      acd194(52)=dotproduct(ninjaA0,spvak2e5)
      acd194(53)=abb194(7)
      acd194(54)=abb194(8)
      acd194(55)=abb194(13)
      acd194(56)=abb194(30)
      acd194(57)=abb194(59)
      acd194(58)=dotproduct(ninjaE3,spvak2k1)
      acd194(59)=abb194(89)
      acd194(60)=dotproduct(ninjaE3,spvak2k5)
      acd194(61)=abb194(91)
      acd194(62)=abb194(23)
      acd194(63)=abb194(87)
      acd194(64)=abb194(25)
      acd194(65)=abb194(26)
      acd194(66)=dotproduct(ninjaE3,spvak1k2)
      acd194(67)=abb194(34)
      acd194(68)=dotproduct(ninjaE3,spvak5k2)
      acd194(69)=abb194(78)
      acd194(70)=acd194(2)*acd194(1)
      acd194(71)=acd194(70)*acd194(4)
      acd194(72)=acd194(71)*acd194(3)
      acd194(73)=acd194(6)*acd194(5)
      acd194(74)=acd194(73)*acd194(4)
      acd194(75)=acd194(74)*acd194(7)
      acd194(72)=acd194(72)+acd194(75)
      acd194(75)=2.0_ki*acd194(9)
      acd194(76)=acd194(15)*acd194(75)
      acd194(77)=-acd194(23)*acd194(21)
      acd194(78)=-acd194(31)*acd194(20)
      acd194(76)=acd194(78)+acd194(77)+acd194(76)
      acd194(71)=acd194(71)*acd194(76)
      acd194(76)=acd194(14)*acd194(75)
      acd194(77)=-acd194(29)*acd194(28)
      acd194(78)=-acd194(30)*acd194(27)
      acd194(76)=acd194(78)+acd194(77)+acd194(76)
      acd194(74)=acd194(74)*acd194(76)
      acd194(76)=acd194(25)*acd194(18)
      acd194(77)=-acd194(24)*acd194(76)
      acd194(78)=acd194(26)*acd194(18)
      acd194(79)=-acd194(32)*acd194(78)
      acd194(77)=acd194(79)+acd194(77)
      acd194(77)=acd194(73)*acd194(77)
      acd194(79)=acd194(17)*acd194(18)
      acd194(80)=-acd194(16)*acd194(79)
      acd194(81)=acd194(19)*acd194(18)
      acd194(82)=-acd194(22)*acd194(81)
      acd194(80)=acd194(82)+acd194(80)
      acd194(80)=acd194(70)*acd194(80)
      acd194(82)=ninjaP1+2.0_ki*acd194(8)
      acd194(82)=acd194(72)*acd194(82)
      acd194(83)=acd194(21)*acd194(20)
      acd194(84)=acd194(75)*acd194(3)
      acd194(83)=acd194(83)-acd194(84)
      acd194(84)=acd194(83)*acd194(2)*acd194(4)
      acd194(85)=acd194(79)*acd194(19)
      acd194(86)=acd194(85)*acd194(2)
      acd194(84)=acd194(86)+acd194(84)
      acd194(86)=-acd194(10)*acd194(84)
      acd194(87)=acd194(83)*acd194(1)*acd194(4)
      acd194(85)=acd194(85)*acd194(1)
      acd194(85)=acd194(85)+acd194(87)
      acd194(87)=-acd194(11)*acd194(85)
      acd194(88)=acd194(75)*acd194(7)
      acd194(89)=acd194(28)*acd194(27)
      acd194(88)=acd194(88)-acd194(89)
      acd194(89)=-acd194(88)*acd194(6)*acd194(4)
      acd194(90)=acd194(76)*acd194(26)
      acd194(91)=acd194(90)*acd194(6)
      acd194(89)=acd194(91)+acd194(89)
      acd194(91)=-acd194(12)*acd194(89)
      acd194(92)=-acd194(88)*acd194(5)*acd194(4)
      acd194(90)=acd194(90)*acd194(5)
      acd194(90)=acd194(90)+acd194(92)
      acd194(92)=-acd194(13)*acd194(90)
      acd194(71)=acd194(92)+acd194(91)+acd194(86)+acd194(87)+acd194(74)+acd194(&
      &71)+acd194(82)+acd194(80)+acd194(77)
      acd194(74)=acd194(33)+ninjaP2
      acd194(74)=acd194(72)*acd194(74)
      acd194(77)=acd194(39)*acd194(75)
      acd194(80)=-acd194(49)*acd194(28)
      acd194(82)=-acd194(50)*acd194(27)
      acd194(77)=acd194(82)+acd194(80)+acd194(77)
      acd194(77)=acd194(4)*acd194(77)
      acd194(80)=acd194(43)*acd194(75)
      acd194(76)=-acd194(48)*acd194(76)
      acd194(78)=-acd194(52)*acd194(78)
      acd194(82)=acd194(62)*acd194(26)
      acd194(86)=acd194(64)*acd194(27)
      acd194(87)=acd194(65)*acd194(28)
      acd194(91)=acd194(67)*acd194(66)
      acd194(92)=acd194(69)*acd194(68)
      acd194(76)=acd194(92)+acd194(91)+acd194(87)+acd194(86)+acd194(82)+acd194(&
      &78)+acd194(76)+acd194(80)+acd194(77)
      acd194(76)=acd194(73)*acd194(76)
      acd194(77)=acd194(40)*acd194(75)
      acd194(78)=-acd194(47)*acd194(21)
      acd194(80)=-acd194(51)*acd194(20)
      acd194(77)=acd194(80)+acd194(78)+acd194(77)
      acd194(77)=acd194(4)*acd194(77)
      acd194(75)=acd194(41)*acd194(75)
      acd194(78)=-acd194(45)*acd194(79)
      acd194(79)=-acd194(46)*acd194(81)
      acd194(80)=acd194(53)*acd194(19)
      acd194(81)=acd194(55)*acd194(20)
      acd194(82)=acd194(56)*acd194(21)
      acd194(86)=acd194(59)*acd194(58)
      acd194(87)=acd194(61)*acd194(60)
      acd194(75)=acd194(87)+acd194(86)+acd194(82)+acd194(81)+acd194(80)+acd194(&
      &79)+acd194(78)+acd194(75)+acd194(77)
      acd194(75)=acd194(70)*acd194(75)
      acd194(77)=acd194(42)*acd194(83)
      acd194(78)=acd194(54)*acd194(17)*acd194(19)
      acd194(77)=acd194(78)+acd194(77)
      acd194(77)=acd194(1)*acd194(77)
      acd194(78)=ninjaP0+acd194(34)
      acd194(78)=acd194(72)*acd194(78)
      acd194(79)=-acd194(35)*acd194(84)
      acd194(80)=-acd194(36)*acd194(85)
      acd194(81)=-acd194(37)*acd194(89)
      acd194(82)=-acd194(38)*acd194(90)
      acd194(83)=-acd194(44)*acd194(5)*acd194(88)
      acd194(73)=acd194(73)*acd194(7)
      acd194(70)=acd194(70)*acd194(3)
      acd194(70)=acd194(70)+acd194(73)
      acd194(70)=acd194(57)*acd194(70)
      acd194(73)=acd194(63)*acd194(25)*acd194(26)*acd194(5)
      acd194(70)=acd194(73)+acd194(70)+acd194(83)+acd194(82)+acd194(81)+acd194(&
      &79)+acd194(80)+acd194(75)+acd194(76)+acd194(78)+acd194(77)
      brack(ninjaidxt0x0mu0)=acd194(70)
      brack(ninjaidxt0x0mu2)=acd194(72)
      brack(ninjaidxt0x1mu0)=acd194(71)
      brack(ninjaidxt0x2mu0)=acd194(74)
   end subroutine brack_22
!---#] subroutine brack_22:
!---#[ subroutine numerator_t2:
   subroutine numerator_t2(ncut, a0, a1, b, c, param, deg, coeffs) &
   & bind(c, name="p9_part21part21_part25part25part21_d194h0_ninja_t2")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki => ki_nin
      use p9_part21part21_part25part25part21_globalsl1, only: epspow
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_abbrevd194h0
      implicit none
      integer(c_int), intent(in) :: ncut, deg
      complex(ki), dimension(0:3), intent(in) :: a0, a1, b, c
      complex(ki), dimension(0:2), intent(in) :: param
      complex(ki), dimension(4) :: vecA0, vecA1, vecB, vecC
      complex(ki), dimension(0:*), intent(out) :: coeffs
      integer :: t1
      complex(ki), dimension(4) :: qshift
      qshift = k2
	     vecA0(1:4) = + a0(0:3) - qshift(1:4)
	     vecA1(1:4) = + a1(0:3)
	     vecB(1:4) = + b(0:3)
	     vecC(1:4) = + c(0:3)
      if (deg.lt.0) return
      t1 = 0
      call cond_t(epspow.eq.t1,brack_21,vecA0,vecA1,vecB,vecC,param,coeffs)
      if (deg.le.(1+(0))) return
      call cond_t(epspow.eq.t1,brack_22,vecA0,vecA1,vecB,vecC,param,coeffs)
   end subroutine numerator_t2
!---#] subroutine numerator_t2:
end module     p9_part21part21_part25part25part21_d194h0l132
