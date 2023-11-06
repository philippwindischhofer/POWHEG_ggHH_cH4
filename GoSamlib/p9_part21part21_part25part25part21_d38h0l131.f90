module     p9_part21part21_part25part25part21_d38h0l131
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p9_part21part &
   ! &21_part25part25part21/helicity0d38h0l131.f90
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
      use p9_part21part21_part25part25part21_abbrevd38h0
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA, ninjaE3, ninjaE4
      complex(ki), intent(in) :: ninjaP
      complex(ki), dimension(14) :: acd38
      complex(ki), dimension (0:*), intent(inout) :: brack
      acd38(1)=dotproduct(ninjaE3,spvae2e5)
      acd38(2)=dotproduct(ninjaE3,spvae5e2)
      acd38(3)=abb38(4)
      acd38(4)=dotproduct(ninjaE3,spvae1e2)
      acd38(5)=dotproduct(ninjaE3,spvae2e1)
      acd38(6)=abb38(7)
      acd38(7)=dotproduct(ninjaE3,spvak1e2)
      acd38(8)=dotproduct(ninjaE3,spvae2k1)
      acd38(9)=abb38(10)
      acd38(10)=dotproduct(ninjaE3,spvak5e2)
      acd38(11)=dotproduct(ninjaE3,spvae2k5)
      acd38(12)=-acd38(8)*acd38(7)
      acd38(13)=-acd38(11)*acd38(10)
      acd38(12)=acd38(13)+acd38(12)
      acd38(12)=acd38(9)*acd38(12)
      acd38(13)=acd38(3)*acd38(2)*acd38(1)
      acd38(14)=acd38(6)*acd38(5)*acd38(4)
      acd38(12)=acd38(13)+acd38(14)+acd38(12)
      brack(ninjaidxt3mu0)=0.0_ki
      brack(ninjaidxt2mu0)=acd38(12)
   end subroutine brack_31
!---#] subroutine brack_31:
!---#[ subroutine brack_32:
   pure subroutine brack_32(ninjaA, ninjaE3, ninjaE4, ninjaP, brack)
      use p9_part21part21_part25part25part21_model
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_color
      use p9_part21part21_part25part25part21_abbrevd38h0
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA, ninjaE3, ninjaE4
      complex(ki), intent(in) :: ninjaP
      complex(ki), dimension(86) :: acd38
      complex(ki), dimension (0:*), intent(inout) :: brack
      acd38(1)=dotproduct(ninjaE3,spvae2e5)
      acd38(2)=dotproduct(ninjaE4,spvae5e2)
      acd38(3)=abb38(4)
      acd38(4)=dotproduct(ninjaE3,spvae5e2)
      acd38(5)=dotproduct(ninjaE4,spvae2e5)
      acd38(6)=dotproduct(ninjaE3,spvak1e2)
      acd38(7)=dotproduct(ninjaE4,spvae2k1)
      acd38(8)=abb38(10)
      acd38(9)=dotproduct(ninjaE3,spvae1e2)
      acd38(10)=dotproduct(ninjaE4,spvae2e1)
      acd38(11)=abb38(7)
      acd38(12)=dotproduct(ninjaE3,spvae2e1)
      acd38(13)=dotproduct(ninjaE4,spvae1e2)
      acd38(14)=dotproduct(ninjaE3,spvae2k1)
      acd38(15)=dotproduct(ninjaE4,spvak1e2)
      acd38(16)=dotproduct(ninjaE3,spvak5e2)
      acd38(17)=dotproduct(ninjaE4,spvae2k5)
      acd38(18)=dotproduct(ninjaE3,spvae2k5)
      acd38(19)=dotproduct(ninjaE4,spvak5e2)
      acd38(20)=abb38(6)
      acd38(21)=dotproduct(ninjaA,ninjaE3)
      acd38(22)=dotproduct(ninjaA,spvae2e5)
      acd38(23)=dotproduct(ninjaA,spvae5e2)
      acd38(24)=dotproduct(ninjaA,spvak1e2)
      acd38(25)=dotproduct(ninjaA,spvae1e2)
      acd38(26)=dotproduct(ninjaA,spvae2e1)
      acd38(27)=dotproduct(ninjaA,spvae2k1)
      acd38(28)=dotproduct(ninjaA,spvak5e2)
      acd38(29)=dotproduct(ninjaA,spvae2k5)
      acd38(30)=abb38(13)
      acd38(31)=abb38(14)
      acd38(32)=abb38(5)
      acd38(33)=abb38(17)
      acd38(34)=abb38(20)
      acd38(35)=dotproduct(ninjaE3,spvak1k2)
      acd38(36)=abb38(8)
      acd38(37)=abb38(18)
      acd38(38)=abb38(22)
      acd38(39)=abb38(16)
      acd38(40)=dotproduct(ninjaE3,spvak5k2)
      acd38(41)=abb38(11)
      acd38(42)=dotproduct(ninjaE3,spvak2k1)
      acd38(43)=abb38(15)
      acd38(44)=dotproduct(ninjaE3,spvae1k2)
      acd38(45)=abb38(19)
      acd38(46)=dotproduct(ninjaE3,spvak2e1)
      acd38(47)=abb38(21)
      acd38(48)=dotproduct(ninjaE3,spvak2k5)
      acd38(49)=abb38(23)
      acd38(50)=dotproduct(ninjaE3,spvae5k2)
      acd38(51)=abb38(24)
      acd38(52)=dotproduct(ninjaE3,spvak2e5)
      acd38(53)=abb38(26)
      acd38(54)=dotproduct(ninjaA,ninjaA)
      acd38(55)=dotproduct(ninjaA,spvak1k2)
      acd38(56)=dotproduct(ninjaA,spvak5k2)
      acd38(57)=dotproduct(ninjaA,spvak2k1)
      acd38(58)=dotproduct(ninjaA,spvae1k2)
      acd38(59)=dotproduct(ninjaA,spvak2e1)
      acd38(60)=dotproduct(ninjaA,spvak2k5)
      acd38(61)=dotproduct(ninjaA,spvae5k2)
      acd38(62)=dotproduct(ninjaA,spvak2e5)
      acd38(63)=abb38(9)
      acd38(64)=acd38(18)*acd38(19)
      acd38(65)=acd38(16)*acd38(17)
      acd38(66)=acd38(14)*acd38(15)
      acd38(67)=acd38(6)*acd38(7)
      acd38(64)=acd38(64)+acd38(65)+acd38(66)+acd38(67)
      acd38(65)=-acd38(8)*acd38(64)
      acd38(66)=acd38(12)*acd38(13)
      acd38(67)=acd38(9)*acd38(10)
      acd38(66)=acd38(66)+acd38(67)
      acd38(66)=acd38(66)*acd38(11)
      acd38(67)=acd38(4)*acd38(5)
      acd38(68)=acd38(1)*acd38(2)
      acd38(67)=acd38(67)+acd38(68)
      acd38(67)=acd38(67)*acd38(3)
      acd38(66)=acd38(20)+acd38(66)+acd38(67)
      acd38(65)=acd38(65)+acd38(66)
      acd38(67)=-acd38(18)*acd38(28)
      acd38(68)=-acd38(16)*acd38(29)
      acd38(69)=-acd38(14)*acd38(24)
      acd38(70)=-acd38(6)*acd38(27)
      acd38(67)=acd38(70)+acd38(69)+acd38(67)+acd38(68)
      acd38(67)=acd38(8)*acd38(67)
      acd38(68)=acd38(12)*acd38(25)
      acd38(69)=acd38(9)*acd38(26)
      acd38(68)=acd38(68)+acd38(69)
      acd38(68)=acd38(11)*acd38(68)
      acd38(69)=acd38(4)*acd38(22)
      acd38(70)=acd38(1)*acd38(23)
      acd38(69)=acd38(69)+acd38(70)
      acd38(69)=acd38(3)*acd38(69)
      acd38(70)=acd38(53)*acd38(52)
      acd38(71)=acd38(51)*acd38(50)
      acd38(72)=acd38(49)*acd38(48)
      acd38(73)=acd38(47)*acd38(46)
      acd38(74)=acd38(45)*acd38(44)
      acd38(75)=acd38(43)*acd38(42)
      acd38(76)=acd38(41)*acd38(40)
      acd38(77)=acd38(36)*acd38(35)
      acd38(78)=acd38(20)*acd38(21)
      acd38(79)=acd38(18)*acd38(39)
      acd38(80)=acd38(16)*acd38(38)
      acd38(81)=acd38(14)*acd38(37)
      acd38(82)=acd38(12)*acd38(34)
      acd38(83)=acd38(9)*acd38(33)
      acd38(84)=acd38(6)*acd38(32)
      acd38(85)=acd38(4)*acd38(31)
      acd38(86)=acd38(1)*acd38(30)
      acd38(67)=acd38(67)+acd38(69)+acd38(68)+acd38(86)+acd38(85)+acd38(84)+acd&
      &38(83)+acd38(82)+acd38(81)+acd38(80)+acd38(79)+2.0_ki*acd38(78)+acd38(77&
      &)+acd38(76)+acd38(75)+acd38(74)+acd38(73)+acd38(72)+acd38(70)+acd38(71)
      acd38(64)=-ninjaP*acd38(64)
      acd38(68)=-acd38(28)*acd38(29)
      acd38(69)=-acd38(24)*acd38(27)
      acd38(64)=acd38(64)+acd38(68)+acd38(69)
      acd38(64)=acd38(8)*acd38(64)
      acd38(66)=ninjaP*acd38(66)
      acd38(68)=acd38(11)*acd38(26)
      acd38(68)=acd38(68)+acd38(33)
      acd38(68)=acd38(25)*acd38(68)
      acd38(69)=acd38(3)*acd38(23)
      acd38(69)=acd38(69)+acd38(30)
      acd38(69)=acd38(22)*acd38(69)
      acd38(70)=acd38(53)*acd38(62)
      acd38(71)=acd38(51)*acd38(61)
      acd38(72)=acd38(49)*acd38(60)
      acd38(73)=acd38(47)*acd38(59)
      acd38(74)=acd38(45)*acd38(58)
      acd38(75)=acd38(43)*acd38(57)
      acd38(76)=acd38(41)*acd38(56)
      acd38(77)=acd38(36)*acd38(55)
      acd38(78)=acd38(29)*acd38(39)
      acd38(79)=acd38(28)*acd38(38)
      acd38(80)=acd38(27)*acd38(37)
      acd38(81)=acd38(26)*acd38(34)
      acd38(82)=acd38(24)*acd38(32)
      acd38(83)=acd38(23)*acd38(31)
      acd38(84)=acd38(20)*acd38(54)
      acd38(64)=acd38(64)+acd38(66)+acd38(84)+acd38(83)+acd38(82)+acd38(81)+acd&
      &38(80)+acd38(79)+acd38(78)+acd38(77)+acd38(76)+acd38(75)+acd38(74)+acd38&
      &(73)+acd38(72)+acd38(71)+acd38(63)+acd38(70)+acd38(69)+acd38(68)
      brack(ninjaidxt1mu0)=acd38(67)
      brack(ninjaidxt1mu2)=0.0_ki
      brack(ninjaidxt0mu0)=acd38(64)
      brack(ninjaidxt0mu2)=acd38(65)
   end subroutine brack_32
!---#] subroutine brack_32:
!---#[ subroutine numerator_t3:
   subroutine numerator_t3(ncut, a, b, c, param, deg, coeffs) &
   & bind(c, name="p9_part21part21_part25part25part21_d38h0_ninja_t3")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki => ki_nin
      use p9_part21part21_part25part25part21_globalsl1, only: epspow
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_abbrevd38h0
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
end module     p9_part21part21_part25part25part21_d38h0l131
