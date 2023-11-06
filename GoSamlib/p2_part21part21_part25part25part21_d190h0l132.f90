module     p2_part21part21_part25part25part21_d190h0l132
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d190h0l132.f90
   ! generator: buildfortran_tn3.py
   use p2_part21part21_part25part25part21_config, only: ki
   use p2_part21part21_part25part25part21_util, only: cond_t, d => metric_tensor
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
      use p2_part21part21_part25part25part21_model
      use p2_part21part21_part25part25part21_kinematics
      use p2_part21part21_part25part25part21_color
      use p2_part21part21_part25part25part21_abbrevd190h0
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(26) :: acd190
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd190(1)=dotproduct(ninjaA0,ninjaE3)
      acd190(2)=dotproduct(ninjaE3,spvae1e2)
      acd190(3)=dotproduct(ninjaE3,spvae2e5)
      acd190(4)=dotproduct(ninjaE3,spvae5e1)
      acd190(5)=abb190(90)
      acd190(6)=dotproduct(ninjaE3,spvae2e1)
      acd190(7)=dotproduct(ninjaE3,spvae5e2)
      acd190(8)=dotproduct(ninjaE3,spvae1e5)
      acd190(9)=dotproduct(ninjaE3,spvak2e1)
      acd190(10)=dotproduct(ninjaE3,spvae5k2)
      acd190(11)=abb190(84)
      acd190(12)=dotproduct(ninjaE3,spvak1e1)
      acd190(13)=dotproduct(ninjaE3,spvae5k1)
      acd190(14)=dotproduct(ninjaE3,spval4e1)
      acd190(15)=dotproduct(ninjaE3,spvae5l4)
      acd190(16)=dotproduct(ninjaE3,spvae1k2)
      acd190(17)=dotproduct(ninjaE3,spvak2e5)
      acd190(18)=dotproduct(ninjaE3,spvae1k1)
      acd190(19)=dotproduct(ninjaE3,spvak1e5)
      acd190(20)=dotproduct(ninjaE3,spvae1l4)
      acd190(21)=dotproduct(ninjaE3,spval4e5)
      acd190(22)=acd190(20)*acd190(21)
      acd190(23)=-acd190(18)*acd190(19)
      acd190(24)=2.0_ki*acd190(1)
      acd190(25)=acd190(8)*acd190(24)
      acd190(22)=acd190(25)+acd190(22)+acd190(23)
      acd190(23)=acd190(6)*acd190(7)
      acd190(22)=acd190(22)*acd190(23)
      acd190(25)=acd190(14)*acd190(15)
      acd190(26)=-acd190(12)*acd190(13)
      acd190(24)=acd190(4)*acd190(24)
      acd190(24)=acd190(24)+acd190(25)+acd190(26)
      acd190(25)=acd190(2)*acd190(3)
      acd190(24)=acd190(24)*acd190(25)
      acd190(22)=acd190(22)+acd190(24)
      acd190(22)=acd190(5)*acd190(22)
      acd190(23)=-acd190(16)*acd190(17)*acd190(23)
      acd190(24)=-acd190(9)*acd190(10)*acd190(25)
      acd190(23)=acd190(23)+acd190(24)
      acd190(23)=acd190(11)*acd190(23)
      acd190(22)=acd190(22)+acd190(23)
      brack(ninjaidxt2x0mu0)=0.0_ki
      brack(ninjaidxt1x0mu0)=acd190(22)
      brack(ninjaidxt1x1mu0)=0.0_ki
   end subroutine brack_21
!---#] subroutine brack_21:
!---#[ subroutine brack_22:
   pure subroutine brack_22(ninjaA0, ninjaA1, ninjaE3, ninjaE4, ninjaP, brack)
      use p2_part21part21_part25part25part21_model
      use p2_part21part21_part25part25part21_kinematics
      use p2_part21part21_part25part25part21_color
      use p2_part21part21_part25part25part21_abbrevd190h0
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(110) :: acd190
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd190(1)=dotproduct(ninjaE3,spvae1e2)
      acd190(2)=dotproduct(ninjaE3,spvae2e5)
      acd190(3)=dotproduct(ninjaE3,spvae5e1)
      acd190(4)=abb190(90)
      acd190(5)=dotproduct(ninjaE3,spvae2e1)
      acd190(6)=dotproduct(ninjaE3,spvae5e2)
      acd190(7)=dotproduct(ninjaE3,spvae1e5)
      acd190(8)=dotproduct(ninjaA0,ninjaA1)
      acd190(9)=dotproduct(ninjaA0,ninjaE3)
      acd190(10)=dotproduct(ninjaA1,spvae1e2)
      acd190(11)=dotproduct(ninjaA1,spvae2e5)
      acd190(12)=dotproduct(ninjaA1,spvae2e1)
      acd190(13)=dotproduct(ninjaA1,spvae5e1)
      acd190(14)=dotproduct(ninjaA1,spvae5e2)
      acd190(15)=dotproduct(ninjaA1,spvae1e5)
      acd190(16)=dotproduct(ninjaE3,spval4e1)
      acd190(17)=dotproduct(ninjaE3,spvae5l4)
      acd190(18)=dotproduct(ninjaE3,spvak1e1)
      acd190(19)=dotproduct(ninjaE3,spvae5k1)
      acd190(20)=dotproduct(ninjaE3,spvae5k2)
      acd190(21)=dotproduct(ninjaE3,spvak2e1)
      acd190(22)=abb190(84)
      acd190(23)=dotproduct(ninjaA1,spval4e1)
      acd190(24)=dotproduct(ninjaA1,spvae1k2)
      acd190(25)=dotproduct(ninjaE3,spvak2e5)
      acd190(26)=dotproduct(ninjaA1,spvak2e5)
      acd190(27)=dotproduct(ninjaE3,spvae1k2)
      acd190(28)=dotproduct(ninjaE3,spvae1l4)
      acd190(29)=dotproduct(ninjaE3,spval4e5)
      acd190(30)=dotproduct(ninjaE3,spvae1k1)
      acd190(31)=dotproduct(ninjaE3,spvak1e5)
      acd190(32)=dotproduct(ninjaA1,spvae5l4)
      acd190(33)=dotproduct(ninjaA1,spvak1e1)
      acd190(34)=dotproduct(ninjaA1,spvae5k1)
      acd190(35)=dotproduct(ninjaA1,spvae1l4)
      acd190(36)=dotproduct(ninjaA1,spval4e5)
      acd190(37)=dotproduct(ninjaA1,spvae5k2)
      acd190(38)=dotproduct(ninjaA1,spvak2e1)
      acd190(39)=dotproduct(ninjaA1,spvae1k1)
      acd190(40)=dotproduct(ninjaA1,spvak1e5)
      acd190(41)=dotproduct(ninjaA1,ninjaA1)
      acd190(42)=dotproduct(ninjaA0,ninjaA0)
      acd190(43)=dotproduct(ninjaA0,spvae1e2)
      acd190(44)=dotproduct(ninjaA0,spvae2e5)
      acd190(45)=dotproduct(ninjaA0,spvae2e1)
      acd190(46)=dotproduct(ninjaA0,spvae5e1)
      acd190(47)=dotproduct(ninjaA0,spvae5e2)
      acd190(48)=dotproduct(ninjaA0,spvae1e5)
      acd190(49)=abb190(5)
      acd190(50)=abb190(39)
      acd190(51)=abb190(63)
      acd190(52)=abb190(91)
      acd190(53)=dotproduct(ninjaA0,spval4e1)
      acd190(54)=dotproduct(ninjaA0,spvae1k2)
      acd190(55)=dotproduct(ninjaA0,spvak2e5)
      acd190(56)=dotproduct(ninjaA0,spvae5l4)
      acd190(57)=dotproduct(ninjaA0,spvak1e1)
      acd190(58)=dotproduct(ninjaA0,spvae5k1)
      acd190(59)=dotproduct(ninjaA0,spvae1l4)
      acd190(60)=dotproduct(ninjaA0,spval4e5)
      acd190(61)=dotproduct(ninjaA0,spvae5k2)
      acd190(62)=dotproduct(ninjaA0,spvak2e1)
      acd190(63)=dotproduct(ninjaA0,spvae1k1)
      acd190(64)=dotproduct(ninjaA0,spvak1e5)
      acd190(65)=abb190(12)
      acd190(66)=abb190(19)
      acd190(67)=dotproduct(ninjaE3,spvak1k5)
      acd190(68)=abb190(38)
      acd190(69)=abb190(82)
      acd190(70)=abb190(64)
      acd190(71)=dotproduct(ninjaE3,spvak1k2)
      acd190(72)=abb190(43)
      acd190(73)=abb190(58)
      acd190(74)=abb190(72)
      acd190(75)=abb190(87)
      acd190(76)=abb190(16)
      acd190(77)=abb190(62)
      acd190(78)=abb190(85)
      acd190(79)=abb190(49)
      acd190(80)=dotproduct(ninjaE3,spvak2k1)
      acd190(81)=abb190(51)
      acd190(82)=abb190(52)
      acd190(83)=dotproduct(ninjaE3,spvak5k1)
      acd190(84)=abb190(80)
      acd190(85)=abb190(83)
      acd190(86)=acd190(2)*acd190(1)
      acd190(87)=acd190(86)*acd190(4)
      acd190(88)=acd190(87)*acd190(3)
      acd190(89)=acd190(6)*acd190(5)
      acd190(90)=acd190(89)*acd190(4)
      acd190(91)=acd190(90)*acd190(7)
      acd190(88)=acd190(88)+acd190(91)
      acd190(91)=2.0_ki*acd190(9)
      acd190(92)=acd190(13)*acd190(91)
      acd190(93)=acd190(23)*acd190(17)
      acd190(94)=acd190(32)*acd190(16)
      acd190(95)=-acd190(33)*acd190(19)
      acd190(96)=-acd190(34)*acd190(18)
      acd190(92)=acd190(96)+acd190(95)+acd190(94)+acd190(93)+acd190(92)
      acd190(87)=acd190(87)*acd190(92)
      acd190(92)=acd190(15)*acd190(91)
      acd190(93)=acd190(35)*acd190(29)
      acd190(94)=acd190(36)*acd190(28)
      acd190(95)=-acd190(39)*acd190(31)
      acd190(96)=-acd190(40)*acd190(30)
      acd190(92)=acd190(96)+acd190(95)+acd190(94)+acd190(93)+acd190(92)
      acd190(90)=acd190(90)*acd190(92)
      acd190(92)=acd190(17)*acd190(16)
      acd190(93)=acd190(19)*acd190(18)
      acd190(94)=acd190(91)*acd190(3)
      acd190(92)=-acd190(92)+acd190(93)-acd190(94)
      acd190(93)=acd190(92)*acd190(2)*acd190(4)
      acd190(94)=acd190(21)*acd190(22)
      acd190(95)=acd190(94)*acd190(20)
      acd190(96)=acd190(95)*acd190(2)
      acd190(93)=acd190(96)+acd190(93)
      acd190(96)=-acd190(10)*acd190(93)
      acd190(97)=acd190(92)*acd190(1)*acd190(4)
      acd190(95)=acd190(95)*acd190(1)
      acd190(95)=acd190(95)+acd190(97)
      acd190(97)=-acd190(11)*acd190(95)
      acd190(98)=acd190(91)*acd190(7)
      acd190(99)=acd190(29)*acd190(28)
      acd190(100)=acd190(31)*acd190(30)
      acd190(98)=acd190(98)+acd190(99)-acd190(100)
      acd190(99)=acd190(98)*acd190(6)*acd190(4)
      acd190(100)=acd190(27)*acd190(22)
      acd190(101)=acd190(100)*acd190(25)
      acd190(102)=acd190(101)*acd190(6)
      acd190(99)=-acd190(102)+acd190(99)
      acd190(102)=acd190(12)*acd190(99)
      acd190(103)=acd190(98)*acd190(5)*acd190(4)
      acd190(101)=acd190(101)*acd190(5)
      acd190(101)=-acd190(101)+acd190(103)
      acd190(103)=acd190(14)*acd190(101)
      acd190(104)=acd190(25)*acd190(22)
      acd190(105)=-acd190(24)*acd190(104)
      acd190(106)=-acd190(26)*acd190(100)
      acd190(105)=acd190(106)+acd190(105)
      acd190(105)=acd190(89)*acd190(105)
      acd190(106)=-acd190(37)*acd190(94)
      acd190(107)=acd190(20)*acd190(22)
      acd190(108)=-acd190(38)*acd190(107)
      acd190(106)=acd190(108)+acd190(106)
      acd190(106)=acd190(86)*acd190(106)
      acd190(108)=ninjaP1+2.0_ki*acd190(8)
      acd190(108)=acd190(88)*acd190(108)
      acd190(87)=acd190(103)+acd190(102)+acd190(96)+acd190(97)+acd190(90)+acd19&
      &0(87)+acd190(108)+acd190(106)+acd190(105)
      acd190(90)=acd190(41)+ninjaP2
      acd190(90)=acd190(88)*acd190(90)
      acd190(96)=acd190(48)*acd190(91)
      acd190(97)=acd190(59)*acd190(29)
      acd190(102)=acd190(60)*acd190(28)
      acd190(103)=-acd190(63)*acd190(31)
      acd190(105)=-acd190(64)*acd190(30)
      acd190(96)=acd190(105)+acd190(103)+acd190(102)+acd190(97)+acd190(96)
      acd190(96)=acd190(4)*acd190(96)
      acd190(97)=acd190(51)*acd190(91)
      acd190(102)=-acd190(54)*acd190(104)
      acd190(100)=-acd190(55)*acd190(100)
      acd190(103)=acd190(77)*acd190(27)
      acd190(104)=-acd190(78)*acd190(25)
      acd190(105)=acd190(79)*acd190(28)
      acd190(106)=acd190(81)*acd190(80)
      acd190(108)=-acd190(82)*acd190(29)
      acd190(109)=-acd190(84)*acd190(83)
      acd190(110)=acd190(85)*acd190(30)
      acd190(96)=acd190(110)+acd190(109)+acd190(108)+acd190(106)+acd190(105)+ac&
      &d190(104)+acd190(103)+acd190(100)+acd190(102)+acd190(97)+acd190(96)
      acd190(96)=acd190(89)*acd190(96)
      acd190(97)=acd190(46)*acd190(91)
      acd190(100)=acd190(53)*acd190(17)
      acd190(102)=acd190(56)*acd190(16)
      acd190(103)=-acd190(57)*acd190(19)
      acd190(104)=-acd190(58)*acd190(18)
      acd190(97)=acd190(104)+acd190(103)+acd190(102)+acd190(100)+acd190(97)
      acd190(97)=acd190(4)*acd190(97)
      acd190(91)=acd190(49)*acd190(91)
      acd190(94)=-acd190(61)*acd190(94)
      acd190(100)=-acd190(62)*acd190(107)
      acd190(102)=acd190(65)*acd190(16)
      acd190(103)=-acd190(66)*acd190(17)
      acd190(104)=-acd190(68)*acd190(67)
      acd190(105)=acd190(70)*acd190(18)
      acd190(106)=acd190(72)*acd190(71)
      acd190(107)=-acd190(73)*acd190(20)
      acd190(108)=acd190(74)*acd190(21)
      acd190(91)=acd190(108)+acd190(107)+acd190(106)+acd190(105)+acd190(104)+ac&
      &d190(103)+acd190(102)+acd190(100)+acd190(94)+acd190(91)+acd190(97)
      acd190(91)=acd190(86)*acd190(91)
      acd190(92)=-acd190(50)*acd190(92)
      acd190(94)=acd190(75)*acd190(21)*acd190(20)
      acd190(92)=acd190(94)+acd190(92)
      acd190(92)=acd190(1)*acd190(92)
      acd190(93)=-acd190(43)*acd190(93)
      acd190(94)=-acd190(44)*acd190(95)
      acd190(95)=acd190(45)*acd190(99)
      acd190(97)=acd190(47)*acd190(101)
      acd190(98)=acd190(52)*acd190(5)*acd190(98)
      acd190(99)=ninjaP0+acd190(42)
      acd190(99)=acd190(88)*acd190(99)
      acd190(89)=acd190(89)*acd190(7)
      acd190(86)=acd190(86)*acd190(3)
      acd190(86)=-acd190(86)-acd190(89)
      acd190(86)=acd190(69)*acd190(86)
      acd190(89)=acd190(76)*acd190(27)*acd190(25)*acd190(5)
      acd190(86)=acd190(89)+acd190(86)+acd190(98)+acd190(97)+acd190(95)+acd190(&
      &93)+acd190(94)+acd190(91)+acd190(96)+acd190(99)+acd190(92)
      brack(ninjaidxt0x0mu0)=acd190(86)
      brack(ninjaidxt0x0mu2)=acd190(88)
      brack(ninjaidxt0x1mu0)=acd190(87)
      brack(ninjaidxt0x2mu0)=acd190(90)
   end subroutine brack_22
!---#] subroutine brack_22:
!---#[ subroutine numerator_t2:
   subroutine numerator_t2(ncut, a0, a1, b, c, param, deg, coeffs) &
   & bind(c, name="p2_part21part21_part25part25part21_d190h0_ninja_t2")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki => ki_nin
      use p2_part21part21_part25part25part21_globalsl1, only: epspow
      use p2_part21part21_part25part25part21_kinematics
      use p2_part21part21_part25part25part21_abbrevd190h0
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
end module     p2_part21part21_part25part25part21_d190h0l132
