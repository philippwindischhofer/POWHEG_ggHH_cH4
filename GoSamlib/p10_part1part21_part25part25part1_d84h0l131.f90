module     p10_part1part21_part25part25part1_d84h0l131
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p10_part1part &
   ! &21_part25part25part1/helicity0d84h0l131.f90
   ! generator: buildfortran_tn3.py
   use p10_part1part21_part25part25part1_config, only: ki
   use p10_part1part21_part25part25part1_util, only: cond_t, d => metric_tensor
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
      use p10_part1part21_part25part25part1_model
      use p10_part1part21_part25part25part1_kinematics
      use p10_part1part21_part25part25part1_color
      use p10_part1part21_part25part25part1_abbrevd84h0
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA, ninjaE3, ninjaE4
      complex(ki), intent(in) :: ninjaP
      complex(ki), dimension(1) :: acd84
      complex(ki), dimension (0:*), intent(inout) :: brack
      brack(ninjaidxt3mu0)=0.0_ki
      brack(ninjaidxt2mu0)=0.0_ki
   end subroutine brack_31
!---#] subroutine brack_31:
!---#[ subroutine brack_32:
   pure subroutine brack_32(ninjaA, ninjaE3, ninjaE4, ninjaP, brack)
      use p10_part1part21_part25part25part1_model
      use p10_part1part21_part25part25part1_kinematics
      use p10_part1part21_part25part25part1_color
      use p10_part1part21_part25part25part1_abbrevd84h0
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA, ninjaE3, ninjaE4
      complex(ki), intent(in) :: ninjaP
      complex(ki), dimension(80) :: acd84
      complex(ki), dimension (0:*), intent(inout) :: brack
      acd84(1)=dotproduct(ninjaA,ninjaE3)
      acd84(2)=abb84(10)
      acd84(3)=dotproduct(ninjaE3,spvak2e2)
      acd84(4)=abb84(20)
      acd84(5)=dotproduct(ninjaE3,spvae2k2)
      acd84(6)=abb84(43)
      acd84(7)=dotproduct(ninjaE3,spvak5k2)
      acd84(8)=abb84(13)
      acd84(9)=dotproduct(ninjaE3,spval3e2)
      acd84(10)=abb84(40)
      acd84(11)=dotproduct(ninjaE3,spvak2k1)
      acd84(12)=abb84(25)
      acd84(13)=dotproduct(ninjaE3,spvae2l3)
      acd84(14)=abb84(37)
      acd84(15)=dotproduct(ninjaE3,spvae2l4)
      acd84(16)=abb84(31)
      acd84(17)=dotproduct(ninjaE3,spval4e2)
      acd84(18)=abb84(34)
      acd84(19)=abb84(7)
      acd84(20)=abb84(23)
      acd84(21)=abb84(35)
      acd84(22)=abb84(44)
      acd84(23)=abb84(18)
      acd84(24)=abb84(26)
      acd84(25)=abb84(46)
      acd84(26)=dotproduct(ninjaE3,spvae2k1)
      acd84(27)=dotproduct(ninjaE3,spvak5e2)
      acd84(28)=abb84(8)
      acd84(29)=abb84(16)
      acd84(30)=abb84(41)
      acd84(31)=abb84(22)
      acd84(32)=abb84(32)
      acd84(33)=dotproduct(ninjaA,ninjaA)
      acd84(34)=dotproduct(ninjaA,spvak2e2)
      acd84(35)=dotproduct(ninjaA,spvae2k2)
      acd84(36)=dotproduct(ninjaA,spvak5k2)
      acd84(37)=dotproduct(ninjaA,spval3e2)
      acd84(38)=dotproduct(ninjaA,spvak2k1)
      acd84(39)=dotproduct(ninjaA,spvae2l3)
      acd84(40)=dotproduct(ninjaA,spvae2l4)
      acd84(41)=dotproduct(ninjaA,spval4e2)
      acd84(42)=abb84(24)
      acd84(43)=dotproduct(ninjaA,spvae2k1)
      acd84(44)=dotproduct(ninjaA,spvak5e2)
      acd84(45)=abb84(15)
      acd84(46)=abb84(19)
      acd84(47)=abb84(9)
      acd84(48)=abb84(17)
      acd84(49)=abb84(12)
      acd84(50)=abb84(39)
      acd84(51)=abb84(21)
      acd84(52)=abb84(36)
      acd84(53)=abb84(30)
      acd84(54)=abb84(33)
      acd84(55)=acd84(4)*acd84(3)
      acd84(56)=acd84(6)*acd84(5)
      acd84(57)=acd84(8)*acd84(7)
      acd84(58)=acd84(10)*acd84(9)
      acd84(59)=acd84(12)*acd84(11)
      acd84(60)=acd84(14)*acd84(13)
      acd84(61)=acd84(16)*acd84(15)
      acd84(62)=acd84(18)*acd84(17)
      acd84(55)=acd84(59)+acd84(60)-acd84(61)+acd84(62)+acd84(55)+acd84(56)+acd&
      &84(57)+acd84(58)
      acd84(56)=4.0_ki*acd84(2)
      acd84(56)=acd84(56)*acd84(1)
      acd84(56)=acd84(55)+acd84(56)
      acd84(55)=acd84(1)*acd84(55)
      acd84(57)=acd84(2)*acd84(1)**2
      acd84(55)=2.0_ki*acd84(57)+acd84(55)
      acd84(57)=acd84(24)*acd84(11)
      acd84(58)=acd84(25)*acd84(17)
      acd84(57)=acd84(57)+acd84(58)
      acd84(58)=acd84(19)*acd84(3)
      acd84(59)=acd84(23)*acd84(9)
      acd84(58)=acd84(59)+acd84(58)+acd84(57)
      acd84(58)=acd84(5)*acd84(58)
      acd84(59)=acd84(20)*acd84(7)
      acd84(60)=acd84(22)*acd84(15)
      acd84(59)=acd84(59)+acd84(60)
      acd84(60)=acd84(3)*acd84(59)
      acd84(61)=acd84(29)*acd84(7)
      acd84(62)=acd84(30)*acd84(15)
      acd84(61)=acd84(61)+acd84(62)
      acd84(62)=acd84(9)*acd84(61)
      acd84(63)=acd84(31)*acd84(11)
      acd84(64)=acd84(32)*acd84(17)
      acd84(63)=acd84(63)+acd84(64)
      acd84(64)=acd84(21)*acd84(3)
      acd84(64)=acd84(64)+acd84(63)
      acd84(64)=acd84(13)*acd84(64)
      acd84(65)=acd84(28)*acd84(27)
      acd84(66)=acd84(26)*acd84(65)
      acd84(55)=acd84(66)+2.0_ki*acd84(55)+acd84(58)+acd84(64)+acd84(62)+acd84(&
      &60)
      acd84(58)=ninjaP+acd84(33)
      acd84(58)=acd84(58)*acd84(56)
      acd84(60)=acd84(34)*acd84(4)
      acd84(62)=acd84(35)*acd84(6)
      acd84(64)=acd84(37)*acd84(10)
      acd84(66)=acd84(39)*acd84(14)
      acd84(60)=acd84(42)+acd84(66)+acd84(64)+acd84(62)+acd84(60)
      acd84(62)=2.0_ki*acd84(1)
      acd84(60)=acd84(62)*acd84(60)
      acd84(64)=acd84(29)*acd84(9)
      acd84(66)=acd84(8)*acd84(62)
      acd84(67)=acd84(20)*acd84(3)
      acd84(64)=acd84(64)+acd84(66)+acd84(67)
      acd84(64)=acd84(36)*acd84(64)
      acd84(66)=acd84(31)*acd84(13)
      acd84(67)=acd84(12)*acd84(62)
      acd84(68)=acd84(24)*acd84(5)
      acd84(66)=acd84(66)+acd84(67)+acd84(68)
      acd84(66)=acd84(38)*acd84(66)
      acd84(67)=acd84(30)*acd84(9)
      acd84(68)=-acd84(16)*acd84(62)
      acd84(69)=acd84(22)*acd84(3)
      acd84(67)=acd84(67)+acd84(68)+acd84(69)
      acd84(67)=acd84(40)*acd84(67)
      acd84(68)=acd84(32)*acd84(13)
      acd84(62)=acd84(18)*acd84(62)
      acd84(69)=acd84(25)*acd84(5)
      acd84(62)=acd84(68)+acd84(62)+acd84(69)
      acd84(62)=acd84(41)*acd84(62)
      acd84(59)=acd84(34)*acd84(59)
      acd84(57)=acd84(35)*acd84(57)
      acd84(61)=acd84(37)*acd84(61)
      acd84(63)=acd84(39)*acd84(63)
      acd84(68)=acd84(34)*acd84(5)
      acd84(69)=acd84(35)*acd84(3)
      acd84(68)=acd84(68)+acd84(69)
      acd84(68)=acd84(19)*acd84(68)
      acd84(69)=acd84(34)*acd84(13)
      acd84(70)=acd84(39)*acd84(3)
      acd84(69)=acd84(69)+acd84(70)
      acd84(69)=acd84(21)*acd84(69)
      acd84(70)=acd84(35)*acd84(9)
      acd84(71)=acd84(37)*acd84(5)
      acd84(70)=acd84(70)+acd84(71)
      acd84(70)=acd84(23)*acd84(70)
      acd84(71)=acd84(44)*acd84(28)
      acd84(71)=acd84(47)+acd84(71)
      acd84(71)=acd84(26)*acd84(71)
      acd84(65)=acd84(43)*acd84(65)
      acd84(72)=acd84(45)*acd84(3)
      acd84(73)=acd84(46)*acd84(5)
      acd84(74)=acd84(48)*acd84(27)
      acd84(75)=acd84(49)*acd84(7)
      acd84(76)=acd84(50)*acd84(9)
      acd84(77)=acd84(51)*acd84(11)
      acd84(78)=acd84(52)*acd84(13)
      acd84(79)=acd84(53)*acd84(15)
      acd84(80)=acd84(54)*acd84(17)
      acd84(57)=acd84(80)+acd84(79)+acd84(78)+acd84(77)+acd84(76)+acd84(75)+acd&
      &84(74)+acd84(73)+acd84(72)+acd84(65)+acd84(62)+acd84(67)+acd84(66)+acd84&
      &(64)+acd84(70)+acd84(69)+acd84(68)+acd84(60)+acd84(71)+acd84(63)+acd84(6&
      &1)+acd84(57)+acd84(59)+acd84(58)
      brack(ninjaidxt1mu0)=acd84(55)
      brack(ninjaidxt1mu2)=0.0_ki
      brack(ninjaidxt0mu0)=acd84(57)
      brack(ninjaidxt0mu2)=acd84(56)
   end subroutine brack_32
!---#] subroutine brack_32:
!---#[ subroutine numerator_t3:
   subroutine numerator_t3(ncut, a, b, c, param, deg, coeffs) &
   & bind(c, name="p10_part1part21_part25part25part1_d84h0_ninja_t3")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki => ki_nin
      use p10_part1part21_part25part25part1_globalsl1, only: epspow
      use p10_part1part21_part25part25part1_kinematics
      use p10_part1part21_part25part25part1_abbrevd84h0
      implicit none
      integer(c_int), intent(in) :: ncut, deg
      complex(ki), dimension(0:3), intent(in) :: a, b, c
      complex(ki), intent(in) :: param
      complex(ki), dimension(4) :: vecA, vecB, vecC
      complex(ki), dimension(0:*), intent(out) :: coeffs
      integer :: t1
	     vecA(1:4) = + a(0:3)
	     vecB(1:4) = + b(0:3)
	     vecC(1:4) = + c(0:3)
      if (deg.lt.0) return
      t1 = 0
      call cond_t(epspow.eq.t1,brack_31,vecA,vecB,vecC,param,coeffs)
      if (deg.le.(1+(0))) return
      call cond_t(epspow.eq.t1,brack_32,vecA,vecB,vecC,param,coeffs)
   end subroutine numerator_t3
!---#] subroutine numerator_t3:
end module     p10_part1part21_part25part25part1_d84h0l131
