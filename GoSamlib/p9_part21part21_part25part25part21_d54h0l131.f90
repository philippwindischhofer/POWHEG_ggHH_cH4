module     p9_part21part21_part25part25part21_d54h0l131
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p9_part21part &
   ! &21_part25part25part21/helicity0d54h0l131.f90
   ! generator: buildfortran_tn3.py
   use p9_part21part21_part25part25part21_config, only: ki
   use p9_part21part21_part25part25part21_util, only: cond_t, d => metric_tensor
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
      use p9_part21part21_part25part25part21_model
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_color
      use p9_part21part21_part25part25part21_abbrevd54h0
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA, ninjaE3, ninjaE4
      complex(ki), intent(in) :: ninjaP
      complex(ki), dimension(9) :: acd54
      complex(ki), dimension (0:*), intent(inout) :: brack
      acd54(1)=dotproduct(ninjaE3,spvae1e2)
      acd54(2)=dotproduct(ninjaE3,spvae5e1)
      acd54(3)=dotproduct(ninjaE3,spvae2e5)
      acd54(4)=abb54(7)
      acd54(5)=dotproduct(ninjaE3,spvae2e1)
      acd54(6)=dotproduct(ninjaE3,spvae1e5)
      acd54(7)=dotproduct(ninjaE3,spvae5e2)
      acd54(8)=-acd54(3)*acd54(2)*acd54(1)
      acd54(9)=-acd54(7)*acd54(6)*acd54(5)
      acd54(8)=acd54(8)+acd54(9)
      acd54(8)=acd54(4)*acd54(8)
      brack(ninjaidxt3mu0)=0.0_ki
      brack(ninjaidxt2mu0)=acd54(8)
   end subroutine brack_31
!---#] subroutine brack_31:
!---#[ subroutine brack_32:
   pure subroutine brack_32(ninjaA, ninjaE3, ninjaE4, ninjaP, brack)
      use p9_part21part21_part25part25part21_model
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_color
      use p9_part21part21_part25part25part21_abbrevd54h0
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA, ninjaE3, ninjaE4
      complex(ki), intent(in) :: ninjaP
      complex(ki), dimension(99) :: acd54
      complex(ki), dimension (0:*), intent(inout) :: brack
      acd54(1)=dotproduct(ninjaE3,spvae1e2)
      acd54(2)=dotproduct(ninjaE3,spvae2e5)
      acd54(3)=dotproduct(ninjaE4,spvae5e1)
      acd54(4)=abb54(7)
      acd54(5)=dotproduct(ninjaE3,spvae5e1)
      acd54(6)=dotproduct(ninjaE4,spvae2e5)
      acd54(7)=abb54(23)
      acd54(8)=dotproduct(ninjaE3,spvae2e1)
      acd54(9)=dotproduct(ninjaE3,spvae5e2)
      acd54(10)=dotproduct(ninjaE4,spvae1e5)
      acd54(11)=dotproduct(ninjaE3,spvae1e5)
      acd54(12)=dotproduct(ninjaE4,spvae5e2)
      acd54(13)=abb54(30)
      acd54(14)=dotproduct(ninjaE4,spvae2e1)
      acd54(15)=abb54(25)
      acd54(16)=dotproduct(ninjaE4,spvae1e2)
      acd54(17)=abb54(9)
      acd54(18)=dotproduct(ninjaA,ninjaE3)
      acd54(19)=dotproduct(ninjaA,spvae1e2)
      acd54(20)=dotproduct(ninjaA,spvae2e1)
      acd54(21)=dotproduct(ninjaA,spvae5e2)
      acd54(22)=dotproduct(ninjaA,spvae2e5)
      acd54(23)=dotproduct(ninjaA,spvae5e1)
      acd54(24)=dotproduct(ninjaA,spvae1e5)
      acd54(25)=dotproduct(ninjaE3,spvak1e5)
      acd54(26)=abb54(5)
      acd54(27)=abb54(6)
      acd54(28)=abb54(18)
      acd54(29)=dotproduct(ninjaE3,spvak5e1)
      acd54(30)=abb54(33)
      acd54(31)=dotproduct(ninjaE3,spvak2e1)
      acd54(32)=abb54(35)
      acd54(33)=abb54(4)
      acd54(34)=abb54(29)
      acd54(35)=dotproduct(ninjaE3,spvae5k1)
      acd54(36)=abb54(26)
      acd54(37)=dotproduct(ninjaE3,spvae1k2)
      acd54(38)=abb54(34)
      acd54(39)=dotproduct(ninjaE3,spvae1k5)
      acd54(40)=abb54(32)
      acd54(41)=abb54(16)
      acd54(42)=abb54(19)
      acd54(43)=abb54(24)
      acd54(44)=abb54(14)
      acd54(45)=abb54(15)
      acd54(46)=abb54(21)
      acd54(47)=dotproduct(ninjaA,ninjaA)
      acd54(48)=abb54(13)
      acd54(49)=dotproduct(ninjaA,spvak1e5)
      acd54(50)=dotproduct(ninjaA,spvae5k1)
      acd54(51)=dotproduct(ninjaA,spvae1k2)
      acd54(52)=dotproduct(ninjaA,spvak5e1)
      acd54(53)=dotproduct(ninjaA,spvae1k5)
      acd54(54)=dotproduct(ninjaA,spvak2e1)
      acd54(55)=abb54(2)
      acd54(56)=abb54(11)
      acd54(57)=abb54(20)
      acd54(58)=abb54(12)
      acd54(59)=abb54(17)
      acd54(60)=abb54(8)
      acd54(61)=abb54(28)
      acd54(62)=abb54(10)
      acd54(63)=abb54(22)
      acd54(64)=abb54(31)
      acd54(65)=acd54(11)*acd54(4)
      acd54(66)=acd54(65)*acd54(12)
      acd54(67)=acd54(9)*acd54(4)
      acd54(68)=acd54(67)*acd54(10)
      acd54(66)=acd54(66)+acd54(68)
      acd54(66)=acd54(66)*acd54(8)
      acd54(68)=acd54(5)*acd54(4)
      acd54(69)=acd54(68)*acd54(6)
      acd54(70)=acd54(2)*acd54(4)
      acd54(71)=acd54(70)*acd54(3)
      acd54(69)=acd54(69)+acd54(71)
      acd54(69)=acd54(69)*acd54(1)
      acd54(71)=acd54(2)*acd54(68)*acd54(16)
      acd54(72)=acd54(9)*acd54(65)*acd54(14)
      acd54(66)=acd54(69)+acd54(66)+acd54(71)+acd54(72)
      acd54(69)=acd54(7)*acd54(1)
      acd54(71)=acd54(13)*acd54(8)
      acd54(72)=acd54(15)*acd54(9)
      acd54(73)=acd54(17)*acd54(2)
      acd54(69)=acd54(69)+acd54(71)+acd54(72)+acd54(73)
      acd54(71)=acd54(69)-acd54(66)
      acd54(72)=acd54(26)*acd54(25)
      acd54(73)=acd54(30)*acd54(29)
      acd54(74)=acd54(32)*acd54(31)
      acd54(72)=acd54(74)+acd54(72)+acd54(73)
      acd54(73)=-acd54(22)*acd54(68)
      acd54(74)=2.0_ki*acd54(18)
      acd54(75)=acd54(7)*acd54(74)
      acd54(76)=-acd54(23)*acd54(70)
      acd54(77)=acd54(27)*acd54(2)
      acd54(78)=-acd54(28)*acd54(5)
      acd54(73)=acd54(78)+acd54(77)+acd54(76)+acd54(75)+acd54(73)+acd54(72)
      acd54(73)=acd54(1)*acd54(73)
      acd54(75)=acd54(36)*acd54(35)
      acd54(76)=acd54(38)*acd54(37)
      acd54(77)=acd54(40)*acd54(39)
      acd54(75)=acd54(77)+acd54(75)+acd54(76)
      acd54(76)=-acd54(21)*acd54(65)
      acd54(77)=acd54(13)*acd54(74)
      acd54(78)=-acd54(24)*acd54(67)
      acd54(79)=acd54(33)*acd54(9)
      acd54(80)=-acd54(34)*acd54(11)
      acd54(76)=acd54(80)+acd54(79)+acd54(78)+acd54(77)+acd54(76)+acd54(75)
      acd54(76)=acd54(8)*acd54(76)
      acd54(77)=acd54(44)*acd54(35)
      acd54(78)=acd54(45)*acd54(37)
      acd54(79)=acd54(46)*acd54(39)
      acd54(68)=acd54(68)*acd54(19)
      acd54(68)=-acd54(77)+acd54(78)+acd54(79)-acd54(68)
      acd54(77)=acd54(17)*acd54(74)
      acd54(77)=acd54(77)+acd54(68)
      acd54(77)=acd54(2)*acd54(77)
      acd54(78)=acd54(41)*acd54(25)
      acd54(79)=acd54(42)*acd54(29)
      acd54(80)=acd54(43)*acd54(31)
      acd54(65)=acd54(65)*acd54(20)
      acd54(65)=-acd54(78)+acd54(79)+acd54(80)-acd54(65)
      acd54(78)=acd54(15)*acd54(74)
      acd54(78)=acd54(78)+acd54(65)
      acd54(78)=acd54(9)*acd54(78)
      acd54(73)=acd54(78)+acd54(77)+acd54(76)+acd54(73)
      acd54(65)=acd54(21)*acd54(65)
      acd54(68)=acd54(22)*acd54(68)
      acd54(66)=-ninjaP*acd54(66)
      acd54(69)=acd54(47)*acd54(69)
      acd54(72)=acd54(19)*acd54(72)
      acd54(75)=acd54(20)*acd54(75)
      acd54(76)=ninjaP*acd54(1)
      acd54(77)=acd54(19)*acd54(74)
      acd54(76)=acd54(76)+acd54(77)
      acd54(76)=acd54(7)*acd54(76)
      acd54(77)=ninjaP*acd54(8)
      acd54(78)=acd54(20)*acd54(74)
      acd54(77)=acd54(77)+acd54(78)
      acd54(77)=acd54(13)*acd54(77)
      acd54(78)=ninjaP*acd54(9)
      acd54(79)=acd54(21)*acd54(74)
      acd54(78)=acd54(78)+acd54(79)
      acd54(78)=acd54(15)*acd54(78)
      acd54(79)=ninjaP*acd54(2)
      acd54(80)=acd54(22)*acd54(74)
      acd54(79)=acd54(79)+acd54(80)
      acd54(79)=acd54(17)*acd54(79)
      acd54(70)=-acd54(19)*acd54(70)
      acd54(80)=acd54(22)*acd54(1)
      acd54(81)=-acd54(4)*acd54(80)
      acd54(70)=acd54(70)+acd54(81)
      acd54(70)=acd54(23)*acd54(70)
      acd54(67)=-acd54(20)*acd54(67)
      acd54(81)=acd54(21)*acd54(8)
      acd54(82)=-acd54(4)*acd54(81)
      acd54(67)=acd54(67)+acd54(82)
      acd54(67)=acd54(24)*acd54(67)
      acd54(82)=acd54(19)*acd54(2)
      acd54(80)=acd54(82)+acd54(80)
      acd54(80)=acd54(27)*acd54(80)
      acd54(82)=-acd54(19)*acd54(5)
      acd54(83)=-acd54(23)*acd54(1)
      acd54(82)=acd54(82)+acd54(83)
      acd54(82)=acd54(28)*acd54(82)
      acd54(83)=acd54(20)*acd54(9)
      acd54(81)=acd54(83)+acd54(81)
      acd54(81)=acd54(33)*acd54(81)
      acd54(83)=-acd54(20)*acd54(11)
      acd54(84)=-acd54(24)*acd54(8)
      acd54(83)=acd54(83)+acd54(84)
      acd54(83)=acd54(34)*acd54(83)
      acd54(84)=acd54(26)*acd54(1)
      acd54(85)=-acd54(41)*acd54(9)
      acd54(84)=acd54(84)+acd54(85)
      acd54(84)=acd54(49)*acd54(84)
      acd54(85)=acd54(36)*acd54(8)
      acd54(86)=-acd54(44)*acd54(2)
      acd54(85)=acd54(85)+acd54(86)
      acd54(85)=acd54(50)*acd54(85)
      acd54(86)=acd54(38)*acd54(8)
      acd54(87)=acd54(45)*acd54(2)
      acd54(86)=acd54(86)+acd54(87)
      acd54(86)=acd54(51)*acd54(86)
      acd54(87)=acd54(30)*acd54(1)
      acd54(88)=acd54(42)*acd54(9)
      acd54(87)=acd54(87)+acd54(88)
      acd54(87)=acd54(52)*acd54(87)
      acd54(88)=acd54(40)*acd54(8)
      acd54(89)=acd54(46)*acd54(2)
      acd54(88)=acd54(88)+acd54(89)
      acd54(88)=acd54(53)*acd54(88)
      acd54(89)=acd54(32)*acd54(1)
      acd54(90)=acd54(43)*acd54(9)
      acd54(89)=acd54(89)+acd54(90)
      acd54(89)=acd54(54)*acd54(89)
      acd54(74)=acd54(48)*acd54(74)
      acd54(90)=acd54(55)*acd54(1)
      acd54(91)=acd54(56)*acd54(8)
      acd54(92)=acd54(57)*acd54(9)
      acd54(93)=acd54(58)*acd54(25)
      acd54(94)=acd54(59)*acd54(2)
      acd54(95)=acd54(60)*acd54(5)
      acd54(96)=acd54(61)*acd54(11)
      acd54(97)=acd54(62)*acd54(35)
      acd54(98)=acd54(63)*acd54(29)
      acd54(99)=acd54(64)*acd54(39)
      acd54(65)=acd54(99)+acd54(98)+acd54(97)+acd54(96)+acd54(95)+acd54(94)+acd&
      &54(93)+acd54(92)+acd54(91)+acd54(90)+acd54(74)+acd54(89)+acd54(88)+acd54&
      &(87)+acd54(86)+acd54(85)+acd54(84)+acd54(83)+acd54(81)+acd54(82)+acd54(8&
      &0)+acd54(69)+acd54(67)+acd54(70)+acd54(79)+acd54(78)+acd54(77)+acd54(76)&
      &+acd54(66)+acd54(68)+acd54(65)+acd54(75)+acd54(72)
      brack(ninjaidxt1mu0)=acd54(73)
      brack(ninjaidxt1mu2)=0.0_ki
      brack(ninjaidxt0mu0)=acd54(65)
      brack(ninjaidxt0mu2)=acd54(71)
   end subroutine brack_32
!---#] subroutine brack_32:
!---#[ subroutine numerator_t3:
   subroutine numerator_t3(ncut, a, b, c, param, deg, coeffs) &
   & bind(c, name="p9_part21part21_part25part25part21_d54h0_ninja_t3")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki => ki_nin
      use p9_part21part21_part25part25part21_globalsl1, only: epspow
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_abbrevd54h0
      implicit none
      integer(c_int), intent(in) :: ncut, deg
      complex(ki), dimension(0:3), intent(in) :: a, b, c
      complex(ki), intent(in) :: param
      complex(ki), dimension(4) :: vecA, vecB, vecC
      complex(ki), dimension(0:*), intent(out) :: coeffs
      integer :: t1
      complex(ki), dimension(4) :: qshift
      qshift = -k2
	     vecA(1:4) = - a(0:3) - qshift(1:4)
	     vecB(1:4) = - b(0:3)
	     vecC(1:4) = - c(0:3)
      if (deg.lt.0) return
      t1 = 0
      call cond_t(epspow.eq.t1,brack_31,vecA,vecB,vecC,param,coeffs)
      if (deg.le.(1+(0))) return
      call cond_t(epspow.eq.t1,brack_32,vecA,vecB,vecC,param,coeffs)
   end subroutine numerator_t3
!---#] subroutine numerator_t3:
end module     p9_part21part21_part25part25part21_d54h0l131
