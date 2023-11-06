module     p9_part21part21_part25part25part21_d186h0l132
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p9_part21part &
   ! &21_part25part25part21/helicity0d186h0l132.f90
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
      use p9_part21part21_part25part25part21_abbrevd186h0
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(26) :: acd186
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd186(1)=dotproduct(ninjaA0,ninjaE3)
      acd186(2)=dotproduct(ninjaE3,spvae1e2)
      acd186(3)=dotproduct(ninjaE3,spvae2e5)
      acd186(4)=dotproduct(ninjaE3,spvae5e1)
      acd186(5)=abb186(90)
      acd186(6)=dotproduct(ninjaE3,spvae2e1)
      acd186(7)=dotproduct(ninjaE3,spvae5e2)
      acd186(8)=dotproduct(ninjaE3,spvae1e5)
      acd186(9)=dotproduct(ninjaE3,spvak2e1)
      acd186(10)=dotproduct(ninjaE3,spvae5k2)
      acd186(11)=abb186(84)
      acd186(12)=dotproduct(ninjaE3,spvak1e1)
      acd186(13)=dotproduct(ninjaE3,spvae5k1)
      acd186(14)=dotproduct(ninjaE3,spval3e1)
      acd186(15)=dotproduct(ninjaE3,spvae5l3)
      acd186(16)=dotproduct(ninjaE3,spvae1k2)
      acd186(17)=dotproduct(ninjaE3,spvak2e5)
      acd186(18)=dotproduct(ninjaE3,spvae1k1)
      acd186(19)=dotproduct(ninjaE3,spvak1e5)
      acd186(20)=dotproduct(ninjaE3,spvae1l3)
      acd186(21)=dotproduct(ninjaE3,spval3e5)
      acd186(22)=acd186(20)*acd186(21)
      acd186(23)=-acd186(18)*acd186(19)
      acd186(24)=2.0_ki*acd186(1)
      acd186(25)=acd186(8)*acd186(24)
      acd186(22)=acd186(25)+acd186(22)+acd186(23)
      acd186(23)=acd186(6)*acd186(7)
      acd186(22)=acd186(22)*acd186(23)
      acd186(25)=acd186(14)*acd186(15)
      acd186(26)=-acd186(12)*acd186(13)
      acd186(24)=acd186(4)*acd186(24)
      acd186(24)=acd186(24)+acd186(25)+acd186(26)
      acd186(25)=acd186(2)*acd186(3)
      acd186(24)=acd186(24)*acd186(25)
      acd186(22)=acd186(22)+acd186(24)
      acd186(22)=acd186(5)*acd186(22)
      acd186(23)=-acd186(16)*acd186(17)*acd186(23)
      acd186(24)=-acd186(9)*acd186(10)*acd186(25)
      acd186(23)=acd186(23)+acd186(24)
      acd186(23)=acd186(11)*acd186(23)
      acd186(22)=acd186(22)+acd186(23)
      brack(ninjaidxt2x0mu0)=0.0_ki
      brack(ninjaidxt1x0mu0)=acd186(22)
      brack(ninjaidxt1x1mu0)=0.0_ki
   end subroutine brack_21
!---#] subroutine brack_21:
!---#[ subroutine brack_22:
   pure subroutine brack_22(ninjaA0, ninjaA1, ninjaE3, ninjaE4, ninjaP, brack)
      use p9_part21part21_part25part25part21_model
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_color
      use p9_part21part21_part25part25part21_abbrevd186h0
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(110) :: acd186
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd186(1)=dotproduct(ninjaE3,spvae1e2)
      acd186(2)=dotproduct(ninjaE3,spvae2e5)
      acd186(3)=dotproduct(ninjaE3,spvae5e1)
      acd186(4)=abb186(90)
      acd186(5)=dotproduct(ninjaE3,spvae2e1)
      acd186(6)=dotproduct(ninjaE3,spvae5e2)
      acd186(7)=dotproduct(ninjaE3,spvae1e5)
      acd186(8)=dotproduct(ninjaA0,ninjaA1)
      acd186(9)=dotproduct(ninjaA0,ninjaE3)
      acd186(10)=dotproduct(ninjaA1,spvae1e2)
      acd186(11)=dotproduct(ninjaA1,spvae2e5)
      acd186(12)=dotproduct(ninjaA1,spvae2e1)
      acd186(13)=dotproduct(ninjaA1,spvae5e1)
      acd186(14)=dotproduct(ninjaA1,spvae5e2)
      acd186(15)=dotproduct(ninjaA1,spvae1e5)
      acd186(16)=dotproduct(ninjaE3,spval3e1)
      acd186(17)=dotproduct(ninjaE3,spvae5l3)
      acd186(18)=dotproduct(ninjaE3,spvak1e1)
      acd186(19)=dotproduct(ninjaE3,spvae5k1)
      acd186(20)=dotproduct(ninjaE3,spvae5k2)
      acd186(21)=dotproduct(ninjaE3,spvak2e1)
      acd186(22)=abb186(84)
      acd186(23)=dotproduct(ninjaA1,spval3e1)
      acd186(24)=dotproduct(ninjaA1,spvae1k2)
      acd186(25)=dotproduct(ninjaE3,spvak2e5)
      acd186(26)=dotproduct(ninjaA1,spvak2e5)
      acd186(27)=dotproduct(ninjaE3,spvae1k2)
      acd186(28)=dotproduct(ninjaE3,spvae1l3)
      acd186(29)=dotproduct(ninjaE3,spval3e5)
      acd186(30)=dotproduct(ninjaE3,spvae1k1)
      acd186(31)=dotproduct(ninjaE3,spvak1e5)
      acd186(32)=dotproduct(ninjaA1,spvae5l3)
      acd186(33)=dotproduct(ninjaA1,spvak1e1)
      acd186(34)=dotproduct(ninjaA1,spvae5k1)
      acd186(35)=dotproduct(ninjaA1,spvae1l3)
      acd186(36)=dotproduct(ninjaA1,spval3e5)
      acd186(37)=dotproduct(ninjaA1,spvae5k2)
      acd186(38)=dotproduct(ninjaA1,spvak2e1)
      acd186(39)=dotproduct(ninjaA1,spvae1k1)
      acd186(40)=dotproduct(ninjaA1,spvak1e5)
      acd186(41)=dotproduct(ninjaA1,ninjaA1)
      acd186(42)=dotproduct(ninjaA0,ninjaA0)
      acd186(43)=dotproduct(ninjaA0,spvae1e2)
      acd186(44)=dotproduct(ninjaA0,spvae2e5)
      acd186(45)=dotproduct(ninjaA0,spvae2e1)
      acd186(46)=dotproduct(ninjaA0,spvae5e1)
      acd186(47)=dotproduct(ninjaA0,spvae5e2)
      acd186(48)=dotproduct(ninjaA0,spvae1e5)
      acd186(49)=abb186(5)
      acd186(50)=abb186(39)
      acd186(51)=abb186(63)
      acd186(52)=abb186(91)
      acd186(53)=dotproduct(ninjaA0,spval3e1)
      acd186(54)=dotproduct(ninjaA0,spvae1k2)
      acd186(55)=dotproduct(ninjaA0,spvak2e5)
      acd186(56)=dotproduct(ninjaA0,spvae5l3)
      acd186(57)=dotproduct(ninjaA0,spvak1e1)
      acd186(58)=dotproduct(ninjaA0,spvae5k1)
      acd186(59)=dotproduct(ninjaA0,spvae1l3)
      acd186(60)=dotproduct(ninjaA0,spval3e5)
      acd186(61)=dotproduct(ninjaA0,spvae5k2)
      acd186(62)=dotproduct(ninjaA0,spvak2e1)
      acd186(63)=dotproduct(ninjaA0,spvae1k1)
      acd186(64)=dotproduct(ninjaA0,spvak1e5)
      acd186(65)=abb186(12)
      acd186(66)=abb186(19)
      acd186(67)=dotproduct(ninjaE3,spvak1k5)
      acd186(68)=abb186(38)
      acd186(69)=abb186(82)
      acd186(70)=abb186(64)
      acd186(71)=dotproduct(ninjaE3,spvak1k2)
      acd186(72)=abb186(43)
      acd186(73)=abb186(58)
      acd186(74)=abb186(72)
      acd186(75)=abb186(87)
      acd186(76)=abb186(16)
      acd186(77)=abb186(62)
      acd186(78)=abb186(85)
      acd186(79)=abb186(49)
      acd186(80)=dotproduct(ninjaE3,spvak2k1)
      acd186(81)=abb186(51)
      acd186(82)=abb186(52)
      acd186(83)=dotproduct(ninjaE3,spvak5k1)
      acd186(84)=abb186(80)
      acd186(85)=abb186(83)
      acd186(86)=acd186(2)*acd186(1)
      acd186(87)=acd186(86)*acd186(4)
      acd186(88)=acd186(87)*acd186(3)
      acd186(89)=acd186(6)*acd186(5)
      acd186(90)=acd186(89)*acd186(4)
      acd186(91)=acd186(90)*acd186(7)
      acd186(88)=acd186(88)+acd186(91)
      acd186(91)=2.0_ki*acd186(9)
      acd186(92)=acd186(13)*acd186(91)
      acd186(93)=acd186(23)*acd186(17)
      acd186(94)=acd186(32)*acd186(16)
      acd186(95)=-acd186(33)*acd186(19)
      acd186(96)=-acd186(34)*acd186(18)
      acd186(92)=acd186(96)+acd186(95)+acd186(94)+acd186(93)+acd186(92)
      acd186(87)=acd186(87)*acd186(92)
      acd186(92)=acd186(15)*acd186(91)
      acd186(93)=acd186(35)*acd186(29)
      acd186(94)=acd186(36)*acd186(28)
      acd186(95)=-acd186(39)*acd186(31)
      acd186(96)=-acd186(40)*acd186(30)
      acd186(92)=acd186(96)+acd186(95)+acd186(94)+acd186(93)+acd186(92)
      acd186(90)=acd186(90)*acd186(92)
      acd186(92)=acd186(17)*acd186(16)
      acd186(93)=acd186(19)*acd186(18)
      acd186(94)=acd186(91)*acd186(3)
      acd186(92)=-acd186(92)+acd186(93)-acd186(94)
      acd186(93)=acd186(92)*acd186(2)*acd186(4)
      acd186(94)=acd186(21)*acd186(22)
      acd186(95)=acd186(94)*acd186(20)
      acd186(96)=acd186(95)*acd186(2)
      acd186(93)=acd186(96)+acd186(93)
      acd186(96)=-acd186(10)*acd186(93)
      acd186(97)=acd186(92)*acd186(1)*acd186(4)
      acd186(95)=acd186(95)*acd186(1)
      acd186(95)=acd186(95)+acd186(97)
      acd186(97)=-acd186(11)*acd186(95)
      acd186(98)=acd186(91)*acd186(7)
      acd186(99)=acd186(29)*acd186(28)
      acd186(100)=acd186(31)*acd186(30)
      acd186(98)=acd186(98)+acd186(99)-acd186(100)
      acd186(99)=acd186(98)*acd186(6)*acd186(4)
      acd186(100)=acd186(27)*acd186(22)
      acd186(101)=acd186(100)*acd186(25)
      acd186(102)=acd186(101)*acd186(6)
      acd186(99)=-acd186(102)+acd186(99)
      acd186(102)=acd186(12)*acd186(99)
      acd186(103)=acd186(98)*acd186(5)*acd186(4)
      acd186(101)=acd186(101)*acd186(5)
      acd186(101)=-acd186(101)+acd186(103)
      acd186(103)=acd186(14)*acd186(101)
      acd186(104)=acd186(25)*acd186(22)
      acd186(105)=-acd186(24)*acd186(104)
      acd186(106)=-acd186(26)*acd186(100)
      acd186(105)=acd186(106)+acd186(105)
      acd186(105)=acd186(89)*acd186(105)
      acd186(106)=-acd186(37)*acd186(94)
      acd186(107)=acd186(20)*acd186(22)
      acd186(108)=-acd186(38)*acd186(107)
      acd186(106)=acd186(108)+acd186(106)
      acd186(106)=acd186(86)*acd186(106)
      acd186(108)=ninjaP1+2.0_ki*acd186(8)
      acd186(108)=acd186(88)*acd186(108)
      acd186(87)=acd186(103)+acd186(102)+acd186(96)+acd186(97)+acd186(90)+acd18&
      &6(87)+acd186(108)+acd186(106)+acd186(105)
      acd186(90)=acd186(41)+ninjaP2
      acd186(90)=acd186(88)*acd186(90)
      acd186(96)=acd186(48)*acd186(91)
      acd186(97)=acd186(59)*acd186(29)
      acd186(102)=acd186(60)*acd186(28)
      acd186(103)=-acd186(63)*acd186(31)
      acd186(105)=-acd186(64)*acd186(30)
      acd186(96)=acd186(105)+acd186(103)+acd186(102)+acd186(97)+acd186(96)
      acd186(96)=acd186(4)*acd186(96)
      acd186(97)=acd186(51)*acd186(91)
      acd186(102)=-acd186(54)*acd186(104)
      acd186(100)=-acd186(55)*acd186(100)
      acd186(103)=acd186(77)*acd186(27)
      acd186(104)=-acd186(78)*acd186(25)
      acd186(105)=acd186(79)*acd186(28)
      acd186(106)=acd186(81)*acd186(80)
      acd186(108)=-acd186(82)*acd186(29)
      acd186(109)=-acd186(84)*acd186(83)
      acd186(110)=acd186(85)*acd186(30)
      acd186(96)=acd186(110)+acd186(109)+acd186(108)+acd186(106)+acd186(105)+ac&
      &d186(104)+acd186(103)+acd186(100)+acd186(102)+acd186(97)+acd186(96)
      acd186(96)=acd186(89)*acd186(96)
      acd186(97)=acd186(46)*acd186(91)
      acd186(100)=acd186(53)*acd186(17)
      acd186(102)=acd186(56)*acd186(16)
      acd186(103)=-acd186(57)*acd186(19)
      acd186(104)=-acd186(58)*acd186(18)
      acd186(97)=acd186(104)+acd186(103)+acd186(102)+acd186(100)+acd186(97)
      acd186(97)=acd186(4)*acd186(97)
      acd186(91)=acd186(49)*acd186(91)
      acd186(94)=-acd186(61)*acd186(94)
      acd186(100)=-acd186(62)*acd186(107)
      acd186(102)=acd186(65)*acd186(16)
      acd186(103)=-acd186(66)*acd186(17)
      acd186(104)=-acd186(68)*acd186(67)
      acd186(105)=acd186(70)*acd186(18)
      acd186(106)=acd186(72)*acd186(71)
      acd186(107)=-acd186(73)*acd186(20)
      acd186(108)=acd186(74)*acd186(21)
      acd186(91)=acd186(108)+acd186(107)+acd186(106)+acd186(105)+acd186(104)+ac&
      &d186(103)+acd186(102)+acd186(100)+acd186(94)+acd186(91)+acd186(97)
      acd186(91)=acd186(86)*acd186(91)
      acd186(92)=-acd186(50)*acd186(92)
      acd186(94)=acd186(75)*acd186(21)*acd186(20)
      acd186(92)=acd186(94)+acd186(92)
      acd186(92)=acd186(1)*acd186(92)
      acd186(93)=-acd186(43)*acd186(93)
      acd186(94)=-acd186(44)*acd186(95)
      acd186(95)=acd186(45)*acd186(99)
      acd186(97)=acd186(47)*acd186(101)
      acd186(98)=acd186(52)*acd186(5)*acd186(98)
      acd186(99)=ninjaP0+acd186(42)
      acd186(99)=acd186(88)*acd186(99)
      acd186(89)=acd186(89)*acd186(7)
      acd186(86)=acd186(86)*acd186(3)
      acd186(86)=-acd186(86)-acd186(89)
      acd186(86)=acd186(69)*acd186(86)
      acd186(89)=acd186(76)*acd186(27)*acd186(25)*acd186(5)
      acd186(86)=acd186(89)+acd186(86)+acd186(98)+acd186(97)+acd186(95)+acd186(&
      &93)+acd186(94)+acd186(91)+acd186(96)+acd186(99)+acd186(92)
      brack(ninjaidxt0x0mu0)=acd186(86)
      brack(ninjaidxt0x0mu2)=acd186(88)
      brack(ninjaidxt0x1mu0)=acd186(87)
      brack(ninjaidxt0x2mu0)=acd186(90)
   end subroutine brack_22
!---#] subroutine brack_22:
!---#[ subroutine numerator_t2:
   subroutine numerator_t2(ncut, a0, a1, b, c, param, deg, coeffs) &
   & bind(c, name="p9_part21part21_part25part25part21_d186h0_ninja_t2")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki => ki_nin
      use p9_part21part21_part25part25part21_globalsl1, only: epspow
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_abbrevd186h0
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
end module     p9_part21part21_part25part25part21_d186h0l132
