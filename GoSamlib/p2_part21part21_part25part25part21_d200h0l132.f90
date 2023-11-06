module     p2_part21part21_part25part25part21_d200h0l132
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d200h0l132.f90
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
      use p2_part21part21_part25part25part21_abbrevd200h0
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(1) :: acd200
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      brack(ninjaidxt2x0mu0)=0.0_ki
      brack(ninjaidxt1x0mu0)=0.0_ki
      brack(ninjaidxt1x1mu0)=0.0_ki
   end subroutine brack_21
!---#] subroutine brack_21:
!---#[ subroutine brack_22:
   pure subroutine brack_22(ninjaA0, ninjaA1, ninjaE3, ninjaE4, ninjaP, brack)
      use p2_part21part21_part25part25part21_model
      use p2_part21part21_part25part25part21_kinematics
      use p2_part21part21_part25part25part21_color
      use p2_part21part21_part25part25part21_abbrevd200h0
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(56) :: acd200
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd200(1)=dotproduct(ninjaA0,ninjaE3)
      acd200(2)=dotproduct(ninjaE3,spvae1e5)
      acd200(3)=abb200(104)
      acd200(4)=dotproduct(ninjaE3,spvae5e1)
      acd200(5)=abb200(98)
      acd200(6)=dotproduct(ninjaE3,spvae5k2)
      acd200(7)=abb200(100)
      acd200(8)=dotproduct(ninjaE3,spval3e1)
      acd200(9)=abb200(55)
      acd200(10)=dotproduct(ninjaE3,spvak2e1)
      acd200(11)=abb200(24)
      acd200(12)=dotproduct(ninjaE3,spvae2k1)
      acd200(13)=abb200(74)
      acd200(14)=dotproduct(ninjaE3,spvae5l3)
      acd200(15)=abb200(59)
      acd200(16)=dotproduct(ninjaE3,spvae2k5)
      acd200(17)=abb200(109)
      acd200(18)=dotproduct(ninjaE3,spvak2e5)
      acd200(19)=abb200(84)
      acd200(20)=dotproduct(ninjaE3,spvae1l3)
      acd200(21)=abb200(47)
      acd200(22)=dotproduct(ninjaE3,spvak1e2)
      acd200(23)=abb200(58)
      acd200(24)=dotproduct(ninjaE3,spvae1k2)
      acd200(25)=abb200(70)
      acd200(26)=dotproduct(ninjaE3,spvak5e2)
      acd200(27)=abb200(111)
      acd200(28)=dotproduct(ninjaE3,spval3e5)
      acd200(29)=abb200(91)
      acd200(30)=abb200(17)
      acd200(31)=abb200(44)
      acd200(32)=abb200(87)
      acd200(33)=abb200(63)
      acd200(34)=abb200(79)
      acd200(35)=abb200(106)
      acd200(36)=abb200(101)
      acd200(37)=dotproduct(ninjaE3,spvae2e1)
      acd200(38)=dotproduct(ninjaE3,spvae5e2)
      acd200(39)=abb200(62)
      acd200(40)=abb200(35)
      acd200(41)=abb200(40)
      acd200(42)=abb200(86)
      acd200(43)=abb200(102)
      acd200(44)=dotproduct(ninjaE3,spvae1e2)
      acd200(45)=dotproduct(ninjaE3,spvae2e5)
      acd200(46)=abb200(88)
      acd200(47)=abb200(99)
      acd200(48)=-acd200(18)*acd200(19)
      acd200(49)=-acd200(28)*acd200(29)
      acd200(50)=acd200(26)*acd200(27)
      acd200(51)=acd200(22)*acd200(23)
      acd200(52)=acd200(24)*acd200(25)
      acd200(53)=acd200(20)*acd200(21)
      acd200(54)=2.0_ki*acd200(1)
      acd200(55)=acd200(5)*acd200(54)
      acd200(48)=acd200(55)+acd200(53)+acd200(52)+acd200(51)+acd200(50)+acd200(&
      &48)+acd200(49)
      acd200(48)=acd200(48)*acd200(54)
      acd200(49)=-acd200(18)*acd200(40)
      acd200(50)=-acd200(28)*acd200(32)
      acd200(51)=acd200(26)*acd200(42)
      acd200(52)=acd200(22)*acd200(41)
      acd200(49)=acd200(52)+acd200(51)+acd200(49)+acd200(50)
      acd200(49)=acd200(20)*acd200(49)
      acd200(50)=acd200(28)*acd200(47)
      acd200(51)=acd200(26)*acd200(46)
      acd200(52)=acd200(22)*acd200(43)
      acd200(50)=acd200(52)+acd200(50)+acd200(51)
      acd200(50)=acd200(24)*acd200(50)
      acd200(51)=-acd200(39)*acd200(44)*acd200(45)
      acd200(48)=acd200(48)+acd200(49)+acd200(51)+acd200(50)
      acd200(48)=acd200(4)*acd200(48)
      acd200(49)=-acd200(6)*acd200(7)
      acd200(50)=acd200(16)*acd200(17)
      acd200(51)=-acd200(14)*acd200(15)
      acd200(52)=acd200(12)*acd200(13)
      acd200(53)=acd200(10)*acd200(11)
      acd200(55)=acd200(8)*acd200(9)
      acd200(56)=acd200(3)*acd200(54)
      acd200(49)=acd200(56)+acd200(55)+acd200(53)+acd200(52)+acd200(51)+acd200(&
      &49)+acd200(50)
      acd200(49)=acd200(49)*acd200(54)
      acd200(50)=-acd200(6)*acd200(30)
      acd200(51)=acd200(16)*acd200(33)
      acd200(52)=-acd200(14)*acd200(32)
      acd200(53)=acd200(12)*acd200(31)
      acd200(50)=acd200(53)+acd200(52)+acd200(50)+acd200(51)
      acd200(50)=acd200(8)*acd200(50)
      acd200(51)=acd200(16)*acd200(36)
      acd200(52)=acd200(14)*acd200(35)
      acd200(53)=acd200(12)*acd200(34)
      acd200(51)=acd200(53)+acd200(51)+acd200(52)
      acd200(51)=acd200(10)*acd200(51)
      acd200(52)=-acd200(39)*acd200(37)*acd200(38)
      acd200(49)=acd200(49)+acd200(50)+acd200(52)+acd200(51)
      acd200(49)=acd200(2)*acd200(49)
      acd200(48)=acd200(48)+acd200(49)
      brack(ninjaidxt0x0mu0)=acd200(48)
      brack(ninjaidxt0x0mu2)=0.0_ki
      brack(ninjaidxt0x1mu0)=0.0_ki
      brack(ninjaidxt0x2mu0)=0.0_ki
   end subroutine brack_22
!---#] subroutine brack_22:
!---#[ subroutine numerator_t2:
   subroutine numerator_t2(ncut, a0, a1, b, c, param, deg, coeffs) &
   & bind(c, name="p2_part21part21_part25part25part21_d200h0_ninja_t2")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki => ki_nin
      use p2_part21part21_part25part25part21_globalsl1, only: epspow
      use p2_part21part21_part25part25part21_kinematics
      use p2_part21part21_part25part25part21_abbrevd200h0
      implicit none
      integer(c_int), intent(in) :: ncut, deg
      complex(ki), dimension(0:3), intent(in) :: a0, a1, b, c
      complex(ki), dimension(0:2), intent(in) :: param
      complex(ki), dimension(4) :: vecA0, vecA1, vecB, vecC
      complex(ki), dimension(0:*), intent(out) :: coeffs
      integer :: t1
      complex(ki), dimension(4) :: qshift
      qshift = -k3
	     vecA0(1:4) = - a0(0:3) - qshift(1:4)
	     vecA1(1:4) = - a1(0:3)
	     vecB(1:4) = - b(0:3)
	     vecC(1:4) = - c(0:3)
      if (deg.lt.0) return
      t1 = 0
      call cond_t(epspow.eq.t1,brack_21,vecA0,vecA1,vecB,vecC,param,coeffs)
      if (deg.le.(1+(0))) return
      call cond_t(epspow.eq.t1,brack_22,vecA0,vecA1,vecB,vecC,param,coeffs)
   end subroutine numerator_t2
!---#] subroutine numerator_t2:
end module     p2_part21part21_part25part25part21_d200h0l132
