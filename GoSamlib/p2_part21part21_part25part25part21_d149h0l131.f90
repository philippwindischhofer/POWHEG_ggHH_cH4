module     p2_part21part21_part25part25part21_d149h0l131
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d149h0l131.f90
   ! generator: buildfortran_tn3.py
   use p2_part21part21_part25part25part21_config, only: ki
   use p2_part21part21_part25part25part21_util, only: cond_t, d => metric_tensor
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
      use p2_part21part21_part25part25part21_model
      use p2_part21part21_part25part25part21_kinematics
      use p2_part21part21_part25part25part21_color
      use p2_part21part21_part25part25part21_abbrevd149h0
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA, ninjaE3, ninjaE4
      complex(ki), intent(in) :: ninjaP
      complex(ki), dimension(9) :: acd149
      complex(ki), dimension (0:*), intent(inout) :: brack
      acd149(1)=dotproduct(ninjaE3,spvae1e2)
      acd149(2)=dotproduct(ninjaE3,spvae5e1)
      acd149(3)=dotproduct(ninjaE3,spvae2e5)
      acd149(4)=abb149(7)
      acd149(5)=dotproduct(ninjaE3,spvae2e1)
      acd149(6)=dotproduct(ninjaE3,spvae1e5)
      acd149(7)=dotproduct(ninjaE3,spvae5e2)
      acd149(8)=-acd149(3)*acd149(2)*acd149(1)
      acd149(9)=-acd149(7)*acd149(6)*acd149(5)
      acd149(8)=acd149(8)+acd149(9)
      acd149(8)=acd149(4)*acd149(8)
      brack(ninjaidxt3mu0)=0.0_ki
      brack(ninjaidxt2mu0)=acd149(8)
   end subroutine brack_31
!---#] subroutine brack_31:
!---#[ subroutine brack_32:
   pure subroutine brack_32(ninjaA, ninjaE3, ninjaE4, ninjaP, brack)
      use p2_part21part21_part25part25part21_model
      use p2_part21part21_part25part25part21_kinematics
      use p2_part21part21_part25part25part21_color
      use p2_part21part21_part25part25part21_abbrevd149h0
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA, ninjaE3, ninjaE4
      complex(ki), intent(in) :: ninjaP
      complex(ki), dimension(99) :: acd149
      complex(ki), dimension (0:*), intent(inout) :: brack
      acd149(1)=dotproduct(ninjaE3,spvae1e5)
      acd149(2)=dotproduct(ninjaE3,spvae5e2)
      acd149(3)=dotproduct(ninjaE4,spvae2e1)
      acd149(4)=abb149(7)
      acd149(5)=dotproduct(ninjaE3,spvae2e1)
      acd149(6)=dotproduct(ninjaE4,spvae5e2)
      acd149(7)=abb149(32)
      acd149(8)=dotproduct(ninjaE4,spvae1e5)
      acd149(9)=abb149(24)
      acd149(10)=dotproduct(ninjaE3,spvae2e5)
      acd149(11)=dotproduct(ninjaE3,spvae5e1)
      acd149(12)=dotproduct(ninjaE4,spvae1e2)
      acd149(13)=dotproduct(ninjaE3,spvae1e2)
      acd149(14)=dotproduct(ninjaE4,spvae5e1)
      acd149(15)=abb149(28)
      acd149(16)=dotproduct(ninjaE4,spvae2e5)
      acd149(17)=abb149(30)
      acd149(18)=dotproduct(ninjaA,ninjaE3)
      acd149(19)=dotproduct(ninjaA,spvae1e5)
      acd149(20)=dotproduct(ninjaA,spvae5e2)
      acd149(21)=dotproduct(ninjaA,spvae2e5)
      acd149(22)=dotproduct(ninjaA,spvae5e1)
      acd149(23)=dotproduct(ninjaA,spvae1e2)
      acd149(24)=dotproduct(ninjaA,spvae2e1)
      acd149(25)=abb149(3)
      acd149(26)=dotproduct(ninjaE3,spvak1e2)
      acd149(27)=abb149(16)
      acd149(28)=abb149(21)
      acd149(29)=dotproduct(ninjaE3,spvak2e1)
      acd149(30)=abb149(26)
      acd149(31)=dotproduct(ninjaE3,spvak5e1)
      acd149(32)=abb149(34)
      acd149(33)=dotproduct(ninjaE3,spvae2k1)
      acd149(34)=abb149(8)
      acd149(35)=dotproduct(ninjaE3,spvae1k2)
      acd149(36)=abb149(23)
      acd149(37)=dotproduct(ninjaE3,spvae1k5)
      acd149(38)=abb149(29)
      acd149(39)=abb149(5)
      acd149(40)=abb149(6)
      acd149(41)=abb149(18)
      acd149(42)=abb149(27)
      acd149(43)=abb149(17)
      acd149(44)=abb149(22)
      acd149(45)=abb149(36)
      acd149(46)=abb149(33)
      acd149(47)=dotproduct(ninjaA,ninjaA)
      acd149(48)=abb149(20)
      acd149(49)=dotproduct(ninjaA,spvak1e2)
      acd149(50)=dotproduct(ninjaA,spvae2k1)
      acd149(51)=dotproduct(ninjaA,spvak2e1)
      acd149(52)=dotproduct(ninjaA,spvae1k2)
      acd149(53)=dotproduct(ninjaA,spvak5e1)
      acd149(54)=dotproduct(ninjaA,spvae1k5)
      acd149(55)=abb149(25)
      acd149(56)=abb149(12)
      acd149(57)=abb149(10)
      acd149(58)=abb149(9)
      acd149(59)=abb149(15)
      acd149(60)=abb149(11)
      acd149(61)=abb149(14)
      acd149(62)=abb149(13)
      acd149(63)=abb149(19)
      acd149(64)=abb149(35)
      acd149(65)=acd149(13)*acd149(4)
      acd149(66)=acd149(65)*acd149(16)
      acd149(67)=acd149(10)*acd149(4)
      acd149(68)=acd149(67)*acd149(12)
      acd149(66)=acd149(66)+acd149(68)
      acd149(66)=acd149(66)*acd149(11)
      acd149(68)=acd149(5)*acd149(4)
      acd149(69)=acd149(68)*acd149(6)
      acd149(70)=acd149(2)*acd149(4)
      acd149(71)=acd149(70)*acd149(3)
      acd149(69)=acd149(69)+acd149(71)
      acd149(69)=acd149(69)*acd149(1)
      acd149(71)=acd149(2)*acd149(68)*acd149(8)
      acd149(72)=acd149(10)*acd149(65)*acd149(14)
      acd149(66)=acd149(69)+acd149(66)+acd149(71)+acd149(72)
      acd149(69)=acd149(7)*acd149(1)
      acd149(71)=acd149(9)*acd149(2)
      acd149(72)=acd149(15)*acd149(10)
      acd149(73)=acd149(17)*acd149(11)
      acd149(69)=acd149(69)+acd149(71)+acd149(72)+acd149(73)
      acd149(71)=acd149(69)-acd149(66)
      acd149(72)=acd149(27)*acd149(26)
      acd149(73)=acd149(30)*acd149(29)
      acd149(74)=acd149(32)*acd149(31)
      acd149(72)=acd149(74)+acd149(72)+acd149(73)
      acd149(73)=-acd149(20)*acd149(68)
      acd149(74)=2.0_ki*acd149(18)
      acd149(75)=acd149(7)*acd149(74)
      acd149(76)=-acd149(24)*acd149(70)
      acd149(77)=acd149(25)*acd149(2)
      acd149(78)=-acd149(28)*acd149(5)
      acd149(73)=acd149(78)+acd149(77)+acd149(76)+acd149(75)+acd149(73)+acd149(&
      &72)
      acd149(73)=acd149(1)*acd149(73)
      acd149(75)=acd149(44)*acd149(33)
      acd149(76)=acd149(45)*acd149(35)
      acd149(77)=acd149(46)*acd149(37)
      acd149(75)=acd149(77)+acd149(75)+acd149(76)
      acd149(76)=-acd149(21)*acd149(65)
      acd149(77)=acd149(17)*acd149(74)
      acd149(78)=-acd149(23)*acd149(67)
      acd149(79)=acd149(40)*acd149(10)
      acd149(80)=-acd149(43)*acd149(13)
      acd149(76)=acd149(80)+acd149(79)+acd149(78)+acd149(77)+acd149(76)+acd149(&
      &75)
      acd149(76)=acd149(11)*acd149(76)
      acd149(77)=acd149(34)*acd149(33)
      acd149(78)=acd149(36)*acd149(35)
      acd149(79)=acd149(38)*acd149(37)
      acd149(68)=acd149(68)*acd149(19)
      acd149(68)=-acd149(77)-acd149(78)+acd149(79)-acd149(68)
      acd149(77)=acd149(9)*acd149(74)
      acd149(77)=acd149(77)+acd149(68)
      acd149(77)=acd149(2)*acd149(77)
      acd149(78)=acd149(39)*acd149(26)
      acd149(79)=acd149(41)*acd149(29)
      acd149(80)=acd149(42)*acd149(31)
      acd149(65)=acd149(65)*acd149(22)
      acd149(65)=-acd149(78)-acd149(79)+acd149(80)-acd149(65)
      acd149(78)=acd149(15)*acd149(74)
      acd149(78)=acd149(78)+acd149(65)
      acd149(78)=acd149(10)*acd149(78)
      acd149(73)=acd149(78)+acd149(77)+acd149(76)+acd149(73)
      acd149(68)=acd149(20)*acd149(68)
      acd149(65)=acd149(21)*acd149(65)
      acd149(66)=-ninjaP*acd149(66)
      acd149(69)=acd149(47)*acd149(69)
      acd149(72)=acd149(19)*acd149(72)
      acd149(75)=acd149(22)*acd149(75)
      acd149(76)=ninjaP*acd149(1)
      acd149(77)=acd149(19)*acd149(74)
      acd149(76)=acd149(76)+acd149(77)
      acd149(76)=acd149(7)*acd149(76)
      acd149(77)=ninjaP*acd149(2)
      acd149(78)=acd149(20)*acd149(74)
      acd149(77)=acd149(77)+acd149(78)
      acd149(77)=acd149(9)*acd149(77)
      acd149(78)=ninjaP*acd149(10)
      acd149(79)=acd149(21)*acd149(74)
      acd149(78)=acd149(78)+acd149(79)
      acd149(78)=acd149(15)*acd149(78)
      acd149(79)=ninjaP*acd149(11)
      acd149(80)=acd149(22)*acd149(74)
      acd149(79)=acd149(79)+acd149(80)
      acd149(79)=acd149(17)*acd149(79)
      acd149(67)=-acd149(22)*acd149(67)
      acd149(80)=acd149(21)*acd149(11)
      acd149(81)=-acd149(4)*acd149(80)
      acd149(67)=acd149(67)+acd149(81)
      acd149(67)=acd149(23)*acd149(67)
      acd149(70)=-acd149(19)*acd149(70)
      acd149(81)=acd149(20)*acd149(1)
      acd149(82)=-acd149(4)*acd149(81)
      acd149(70)=acd149(70)+acd149(82)
      acd149(70)=acd149(24)*acd149(70)
      acd149(82)=acd149(19)*acd149(2)
      acd149(81)=acd149(82)+acd149(81)
      acd149(81)=acd149(25)*acd149(81)
      acd149(82)=-acd149(19)*acd149(5)
      acd149(83)=-acd149(24)*acd149(1)
      acd149(82)=acd149(82)+acd149(83)
      acd149(82)=acd149(28)*acd149(82)
      acd149(83)=acd149(22)*acd149(10)
      acd149(80)=acd149(83)+acd149(80)
      acd149(80)=acd149(40)*acd149(80)
      acd149(83)=-acd149(22)*acd149(13)
      acd149(84)=-acd149(23)*acd149(11)
      acd149(83)=acd149(83)+acd149(84)
      acd149(83)=acd149(43)*acd149(83)
      acd149(84)=acd149(27)*acd149(1)
      acd149(85)=-acd149(39)*acd149(10)
      acd149(84)=acd149(84)+acd149(85)
      acd149(84)=acd149(49)*acd149(84)
      acd149(85)=-acd149(34)*acd149(2)
      acd149(86)=acd149(44)*acd149(11)
      acd149(85)=acd149(85)+acd149(86)
      acd149(85)=acd149(50)*acd149(85)
      acd149(86)=acd149(30)*acd149(1)
      acd149(87)=-acd149(41)*acd149(10)
      acd149(86)=acd149(86)+acd149(87)
      acd149(86)=acd149(51)*acd149(86)
      acd149(87)=-acd149(36)*acd149(2)
      acd149(88)=acd149(45)*acd149(11)
      acd149(87)=acd149(87)+acd149(88)
      acd149(87)=acd149(52)*acd149(87)
      acd149(88)=acd149(32)*acd149(1)
      acd149(89)=acd149(42)*acd149(10)
      acd149(88)=acd149(88)+acd149(89)
      acd149(88)=acd149(53)*acd149(88)
      acd149(89)=acd149(38)*acd149(2)
      acd149(90)=acd149(46)*acd149(11)
      acd149(89)=acd149(89)+acd149(90)
      acd149(89)=acd149(54)*acd149(89)
      acd149(74)=acd149(48)*acd149(74)
      acd149(90)=acd149(55)*acd149(1)
      acd149(91)=acd149(56)*acd149(2)
      acd149(92)=acd149(57)*acd149(26)
      acd149(93)=acd149(58)*acd149(10)
      acd149(94)=acd149(59)*acd149(11)
      acd149(95)=acd149(60)*acd149(13)
      acd149(96)=acd149(61)*acd149(5)
      acd149(97)=acd149(62)*acd149(33)
      acd149(98)=acd149(63)*acd149(29)
      acd149(99)=acd149(64)*acd149(35)
      acd149(65)=acd149(99)+acd149(98)+acd149(97)+acd149(96)+acd149(95)+acd149(&
      &94)+acd149(93)+acd149(92)+acd149(91)+acd149(90)+acd149(74)+acd149(89)+ac&
      &d149(88)+acd149(87)+acd149(86)+acd149(85)+acd149(84)+acd149(83)+acd149(8&
      &0)+acd149(82)+acd149(81)+acd149(69)+acd149(70)+acd149(67)+acd149(79)+acd&
      &149(78)+acd149(77)+acd149(76)+acd149(66)+acd149(65)+acd149(68)+acd149(75&
      &)+acd149(72)
      brack(ninjaidxt1mu0)=acd149(73)
      brack(ninjaidxt1mu2)=0.0_ki
      brack(ninjaidxt0mu0)=acd149(65)
      brack(ninjaidxt0mu2)=acd149(71)
   end subroutine brack_32
!---#] subroutine brack_32:
!---#[ subroutine numerator_t3:
   subroutine numerator_t3(ncut, a, b, c, param, deg, coeffs) &
   & bind(c, name="p2_part21part21_part25part25part21_d149h0_ninja_t3")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki => ki_nin
      use p2_part21part21_part25part25part21_globalsl1, only: epspow
      use p2_part21part21_part25part25part21_kinematics
      use p2_part21part21_part25part25part21_abbrevd149h0
      implicit none
      integer(c_int), intent(in) :: ncut, deg
      complex(ki), dimension(0:3), intent(in) :: a, b, c
      complex(ki), intent(in) :: param
      complex(ki), dimension(4) :: vecA, vecB, vecC
      complex(ki), dimension(0:*), intent(out) :: coeffs
      integer :: t1
      complex(ki), dimension(4) :: qshift
      qshift = -k2+k5
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
end module     p2_part21part21_part25part25part21_d149h0l131
