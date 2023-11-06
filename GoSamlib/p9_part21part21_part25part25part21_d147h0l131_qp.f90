module     p9_part21part21_part25part25part21_d147h0l131_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p9_part21part &
   ! &21_part25part25part21/helicity0d147h0l131_qp.f90
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
      use p9_part21part21_part25part25part21_abbrevd147h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA, ninjaE3, ninjaE4
      complex(ki), intent(in) :: ninjaP
      complex(ki), dimension(9) :: acd147
      complex(ki), dimension (0:*), intent(inout) :: brack
      acd147(1)=dotproduct(ninjaE3,spvae1e2)
      acd147(2)=dotproduct(ninjaE3,spvae5e1)
      acd147(3)=dotproduct(ninjaE3,spvae2e5)
      acd147(4)=abb147(7)
      acd147(5)=dotproduct(ninjaE3,spvae2e1)
      acd147(6)=dotproduct(ninjaE3,spvae1e5)
      acd147(7)=dotproduct(ninjaE3,spvae5e2)
      acd147(8)=acd147(3)*acd147(2)*acd147(1)
      acd147(9)=acd147(7)*acd147(6)*acd147(5)
      acd147(8)=acd147(8)+acd147(9)
      acd147(8)=acd147(4)*acd147(8)
      brack(ninjaidxt3mu0)=0.0_ki
      brack(ninjaidxt2mu0)=acd147(8)
   end subroutine brack_31
!---#] subroutine brack_31:
!---#[ subroutine brack_32:
   pure subroutine brack_32(ninjaA, ninjaE3, ninjaE4, ninjaP, brack)
      use p9_part21part21_part25part25part21_model_qp
      use p9_part21part21_part25part25part21_kinematics_qp
      use p9_part21part21_part25part25part21_color_qp
      use p9_part21part21_part25part25part21_abbrevd147h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA, ninjaE3, ninjaE4
      complex(ki), intent(in) :: ninjaP
      complex(ki), dimension(99) :: acd147
      complex(ki), dimension (0:*), intent(inout) :: brack
      acd147(1)=dotproduct(ninjaE3,spvae1e2)
      acd147(2)=dotproduct(ninjaE3,spvae2e5)
      acd147(3)=dotproduct(ninjaE4,spvae5e1)
      acd147(4)=abb147(7)
      acd147(5)=dotproduct(ninjaE3,spvae5e1)
      acd147(6)=dotproduct(ninjaE4,spvae2e5)
      acd147(7)=abb147(32)
      acd147(8)=dotproduct(ninjaE4,spvae1e2)
      acd147(9)=abb147(28)
      acd147(10)=dotproduct(ninjaE3,spvae5e2)
      acd147(11)=dotproduct(ninjaE3,spvae2e1)
      acd147(12)=dotproduct(ninjaE4,spvae1e5)
      acd147(13)=dotproduct(ninjaE3,spvae1e5)
      acd147(14)=dotproduct(ninjaE4,spvae2e1)
      acd147(15)=abb147(25)
      acd147(16)=dotproduct(ninjaE4,spvae5e2)
      acd147(17)=abb147(29)
      acd147(18)=dotproduct(ninjaA,ninjaE3)
      acd147(19)=dotproduct(ninjaA,spvae1e2)
      acd147(20)=dotproduct(ninjaA,spvae2e5)
      acd147(21)=dotproduct(ninjaA,spvae5e2)
      acd147(22)=dotproduct(ninjaA,spvae2e1)
      acd147(23)=dotproduct(ninjaA,spvae5e1)
      acd147(24)=dotproduct(ninjaA,spvae1e5)
      acd147(25)=abb147(3)
      acd147(26)=dotproduct(ninjaE3,spvak1e5)
      acd147(27)=abb147(16)
      acd147(28)=abb147(21)
      acd147(29)=dotproduct(ninjaE3,spvak5e1)
      acd147(30)=abb147(23)
      acd147(31)=dotproduct(ninjaE3,spvak2e1)
      acd147(32)=abb147(36)
      acd147(33)=dotproduct(ninjaE3,spvae5k1)
      acd147(34)=abb147(8)
      acd147(35)=dotproduct(ninjaE3,spvae1k5)
      acd147(36)=abb147(14)
      acd147(37)=dotproduct(ninjaE3,spvae1k2)
      acd147(38)=abb147(27)
      acd147(39)=abb147(5)
      acd147(40)=abb147(6)
      acd147(41)=abb147(18)
      acd147(42)=abb147(24)
      acd147(43)=abb147(17)
      acd147(44)=abb147(26)
      acd147(45)=abb147(34)
      acd147(46)=abb147(35)
      acd147(47)=dotproduct(ninjaA,ninjaA)
      acd147(48)=abb147(20)
      acd147(49)=dotproduct(ninjaA,spvak1e5)
      acd147(50)=dotproduct(ninjaA,spvae5k1)
      acd147(51)=dotproduct(ninjaA,spvae1k5)
      acd147(52)=dotproduct(ninjaA,spvak5e1)
      acd147(53)=dotproduct(ninjaA,spvak2e1)
      acd147(54)=dotproduct(ninjaA,spvae1k2)
      acd147(55)=abb147(11)
      acd147(56)=abb147(12)
      acd147(57)=abb147(10)
      acd147(58)=abb147(9)
      acd147(59)=abb147(15)
      acd147(60)=abb147(30)
      acd147(61)=abb147(31)
      acd147(62)=abb147(13)
      acd147(63)=abb147(33)
      acd147(64)=abb147(19)
      acd147(65)=acd147(13)*acd147(4)
      acd147(66)=acd147(65)*acd147(16)
      acd147(67)=acd147(10)*acd147(4)
      acd147(68)=acd147(67)*acd147(12)
      acd147(66)=acd147(66)+acd147(68)
      acd147(66)=acd147(66)*acd147(11)
      acd147(68)=acd147(5)*acd147(4)
      acd147(69)=acd147(68)*acd147(6)
      acd147(70)=acd147(2)*acd147(4)
      acd147(71)=acd147(70)*acd147(3)
      acd147(69)=acd147(69)+acd147(71)
      acd147(69)=acd147(69)*acd147(1)
      acd147(71)=acd147(2)*acd147(68)*acd147(8)
      acd147(72)=acd147(10)*acd147(65)*acd147(14)
      acd147(66)=acd147(69)+acd147(66)+acd147(71)+acd147(72)
      acd147(69)=acd147(7)*acd147(1)
      acd147(71)=acd147(9)*acd147(2)
      acd147(72)=acd147(15)*acd147(10)
      acd147(73)=acd147(17)*acd147(11)
      acd147(69)=acd147(69)+acd147(71)+acd147(72)+acd147(73)
      acd147(71)=acd147(69)+acd147(66)
      acd147(72)=acd147(27)*acd147(26)
      acd147(73)=acd147(30)*acd147(29)
      acd147(74)=acd147(32)*acd147(31)
      acd147(72)=acd147(74)+acd147(72)+acd147(73)
      acd147(73)=acd147(20)*acd147(68)
      acd147(74)=2.0_ki*acd147(18)
      acd147(75)=acd147(7)*acd147(74)
      acd147(76)=acd147(23)*acd147(70)
      acd147(77)=acd147(25)*acd147(2)
      acd147(78)=acd147(28)*acd147(5)
      acd147(73)=acd147(78)+acd147(77)+acd147(76)+acd147(75)+acd147(73)+acd147(&
      &72)
      acd147(73)=acd147(1)*acd147(73)
      acd147(75)=acd147(44)*acd147(33)
      acd147(76)=acd147(45)*acd147(35)
      acd147(77)=acd147(46)*acd147(37)
      acd147(75)=acd147(77)+acd147(75)+acd147(76)
      acd147(76)=acd147(21)*acd147(65)
      acd147(77)=acd147(17)*acd147(74)
      acd147(78)=acd147(24)*acd147(67)
      acd147(79)=acd147(40)*acd147(10)
      acd147(80)=acd147(43)*acd147(13)
      acd147(76)=acd147(80)+acd147(79)+acd147(78)+acd147(77)+acd147(76)+acd147(&
      &75)
      acd147(76)=acd147(11)*acd147(76)
      acd147(77)=acd147(34)*acd147(33)
      acd147(78)=acd147(36)*acd147(35)
      acd147(79)=acd147(38)*acd147(37)
      acd147(68)=acd147(68)*acd147(19)
      acd147(68)=acd147(77)-acd147(78)+acd147(79)+acd147(68)
      acd147(77)=acd147(9)*acd147(74)
      acd147(77)=acd147(77)+acd147(68)
      acd147(77)=acd147(2)*acd147(77)
      acd147(78)=acd147(39)*acd147(26)
      acd147(79)=acd147(41)*acd147(29)
      acd147(80)=acd147(42)*acd147(31)
      acd147(65)=acd147(65)*acd147(22)
      acd147(65)=acd147(78)-acd147(79)+acd147(80)+acd147(65)
      acd147(78)=acd147(15)*acd147(74)
      acd147(78)=acd147(78)+acd147(65)
      acd147(78)=acd147(10)*acd147(78)
      acd147(73)=acd147(78)+acd147(77)+acd147(76)+acd147(73)
      acd147(68)=acd147(20)*acd147(68)
      acd147(65)=acd147(21)*acd147(65)
      acd147(66)=ninjaP*acd147(66)
      acd147(69)=acd147(47)*acd147(69)
      acd147(72)=acd147(19)*acd147(72)
      acd147(75)=acd147(22)*acd147(75)
      acd147(76)=ninjaP*acd147(1)
      acd147(77)=acd147(19)*acd147(74)
      acd147(76)=acd147(76)+acd147(77)
      acd147(76)=acd147(7)*acd147(76)
      acd147(77)=ninjaP*acd147(2)
      acd147(78)=acd147(20)*acd147(74)
      acd147(77)=acd147(77)+acd147(78)
      acd147(77)=acd147(9)*acd147(77)
      acd147(78)=ninjaP*acd147(10)
      acd147(79)=acd147(21)*acd147(74)
      acd147(78)=acd147(78)+acd147(79)
      acd147(78)=acd147(15)*acd147(78)
      acd147(79)=ninjaP*acd147(11)
      acd147(80)=acd147(22)*acd147(74)
      acd147(79)=acd147(79)+acd147(80)
      acd147(79)=acd147(17)*acd147(79)
      acd147(70)=acd147(19)*acd147(70)
      acd147(80)=acd147(20)*acd147(1)
      acd147(81)=acd147(4)*acd147(80)
      acd147(70)=acd147(70)+acd147(81)
      acd147(70)=acd147(23)*acd147(70)
      acd147(67)=acd147(22)*acd147(67)
      acd147(81)=acd147(21)*acd147(11)
      acd147(82)=acd147(4)*acd147(81)
      acd147(67)=acd147(67)+acd147(82)
      acd147(67)=acd147(24)*acd147(67)
      acd147(82)=acd147(19)*acd147(2)
      acd147(80)=acd147(82)+acd147(80)
      acd147(80)=acd147(25)*acd147(80)
      acd147(82)=acd147(19)*acd147(5)
      acd147(83)=acd147(23)*acd147(1)
      acd147(82)=acd147(82)+acd147(83)
      acd147(82)=acd147(28)*acd147(82)
      acd147(83)=acd147(22)*acd147(10)
      acd147(81)=acd147(83)+acd147(81)
      acd147(81)=acd147(40)*acd147(81)
      acd147(83)=acd147(22)*acd147(13)
      acd147(84)=acd147(24)*acd147(11)
      acd147(83)=acd147(83)+acd147(84)
      acd147(83)=acd147(43)*acd147(83)
      acd147(84)=acd147(27)*acd147(1)
      acd147(85)=acd147(39)*acd147(10)
      acd147(84)=acd147(84)+acd147(85)
      acd147(84)=acd147(49)*acd147(84)
      acd147(85)=acd147(34)*acd147(2)
      acd147(86)=acd147(44)*acd147(11)
      acd147(85)=acd147(85)+acd147(86)
      acd147(85)=acd147(50)*acd147(85)
      acd147(86)=-acd147(36)*acd147(2)
      acd147(87)=acd147(45)*acd147(11)
      acd147(86)=acd147(86)+acd147(87)
      acd147(86)=acd147(51)*acd147(86)
      acd147(87)=acd147(30)*acd147(1)
      acd147(88)=-acd147(41)*acd147(10)
      acd147(87)=acd147(87)+acd147(88)
      acd147(87)=acd147(52)*acd147(87)
      acd147(88)=acd147(32)*acd147(1)
      acd147(89)=acd147(42)*acd147(10)
      acd147(88)=acd147(88)+acd147(89)
      acd147(88)=acd147(53)*acd147(88)
      acd147(89)=acd147(38)*acd147(2)
      acd147(90)=acd147(46)*acd147(11)
      acd147(89)=acd147(89)+acd147(90)
      acd147(89)=acd147(54)*acd147(89)
      acd147(74)=acd147(48)*acd147(74)
      acd147(90)=acd147(55)*acd147(1)
      acd147(91)=acd147(56)*acd147(2)
      acd147(92)=acd147(57)*acd147(26)
      acd147(93)=acd147(58)*acd147(10)
      acd147(94)=acd147(59)*acd147(11)
      acd147(95)=acd147(60)*acd147(5)
      acd147(96)=acd147(61)*acd147(13)
      acd147(97)=acd147(62)*acd147(33)
      acd147(98)=acd147(63)*acd147(35)
      acd147(99)=acd147(64)*acd147(29)
      acd147(65)=acd147(99)+acd147(98)+acd147(97)+acd147(96)+acd147(95)+acd147(&
      &94)+acd147(93)+acd147(92)+acd147(91)+acd147(90)+acd147(74)+acd147(89)+ac&
      &d147(88)+acd147(87)+acd147(86)+acd147(85)+acd147(84)+acd147(83)+acd147(8&
      &1)+acd147(82)+acd147(80)+acd147(69)+acd147(67)+acd147(70)+acd147(79)+acd&
      &147(78)+acd147(77)+acd147(76)+acd147(66)+acd147(65)+acd147(68)+acd147(75&
      &)+acd147(72)
      brack(ninjaidxt1mu0)=acd147(73)
      brack(ninjaidxt1mu2)=0.0_ki
      brack(ninjaidxt0mu0)=acd147(65)
      brack(ninjaidxt0mu2)=acd147(71)
   end subroutine brack_32
!---#] subroutine brack_32:
!---#[ subroutine numerator_t3:
   subroutine numerator_t3(ncut, a, b, c, param, deg, coeffs) &
   & bind(c, name="p9_part21part21_part25part25part21_d147h0_qp_ninja_t3")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki => ki_nin
      use p9_part21part21_part25part25part21_globalsl1_qp, only: epspow
      use p9_part21part21_part25part25part21_kinematics_qp
      use p9_part21part21_part25part25part21_abbrevd147h0_qp
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
end module     p9_part21part21_part25part25part21_d147h0l131_qp
