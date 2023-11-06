module     p2_part21part21_part25part25part21_d42h0l132
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d42h0l132.f90
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
      use p2_part21part21_part25part25part21_abbrevd42h0
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(7) :: acd42
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd42(1)=dotproduct(ninjaE3,spvak2e1)
      acd42(2)=dotproduct(ninjaE3,spvae1k2)
      acd42(3)=abb42(13)
      acd42(4)=dotproduct(ninjaE3,spvak5e1)
      acd42(5)=dotproduct(ninjaE3,spvae1k5)
      acd42(6)=-acd42(2)*acd42(1)
      acd42(7)=-acd42(5)*acd42(4)
      acd42(6)=acd42(6)+acd42(7)
      acd42(6)=acd42(3)*acd42(6)
      brack(ninjaidxt2x0mu0)=0.0_ki
      brack(ninjaidxt1x0mu0)=acd42(6)
      brack(ninjaidxt1x1mu0)=0.0_ki
   end subroutine brack_21
!---#] subroutine brack_21:
!---#[ subroutine brack_22:
   pure subroutine brack_22(ninjaA0, ninjaA1, ninjaE3, ninjaE4, ninjaP, brack)
      use p2_part21part21_part25part25part21_model
      use p2_part21part21_part25part25part21_kinematics
      use p2_part21part21_part25part25part21_color
      use p2_part21part21_part25part25part21_abbrevd42h0
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(38) :: acd42
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd42(1)=dotproduct(ninjaA1,spvae1k5)
      acd42(2)=dotproduct(ninjaE3,spvak5e1)
      acd42(3)=abb42(13)
      acd42(4)=dotproduct(ninjaA1,spvak5e1)
      acd42(5)=dotproduct(ninjaE3,spvae1k5)
      acd42(6)=dotproduct(ninjaA1,spvae1k2)
      acd42(7)=dotproduct(ninjaE3,spvak2e1)
      acd42(8)=dotproduct(ninjaA1,spvak2e1)
      acd42(9)=dotproduct(ninjaE3,spvae1k2)
      acd42(10)=dotproduct(ninjaA0,ninjaE3)
      acd42(11)=abb42(9)
      acd42(12)=dotproduct(ninjaA0,spvae1k5)
      acd42(13)=dotproduct(ninjaA0,spvak5e1)
      acd42(14)=dotproduct(ninjaA0,spvae1k2)
      acd42(15)=dotproduct(ninjaA0,spvak2e1)
      acd42(16)=dotproduct(ninjaE3,spvak1k2)
      acd42(17)=abb42(5)
      acd42(18)=abb42(6)
      acd42(19)=dotproduct(ninjaE3,spvak1k5)
      acd42(20)=abb42(7)
      acd42(21)=abb42(8)
      acd42(22)=abb42(11)
      acd42(23)=dotproduct(ninjaE3,spvak2k1)
      acd42(24)=abb42(12)
      acd42(25)=abb42(14)
      acd42(26)=dotproduct(ninjaE3,spvak5k1)
      acd42(27)=abb42(15)
      acd42(28)=-acd42(9)*acd42(8)
      acd42(29)=-acd42(7)*acd42(6)
      acd42(30)=-acd42(5)*acd42(4)
      acd42(31)=-acd42(2)*acd42(1)
      acd42(28)=acd42(31)+acd42(30)+acd42(28)+acd42(29)
      acd42(28)=acd42(3)*acd42(28)
      acd42(29)=-acd42(9)*acd42(15)
      acd42(30)=-acd42(7)*acd42(14)
      acd42(31)=-acd42(5)*acd42(13)
      acd42(32)=-acd42(2)*acd42(12)
      acd42(29)=acd42(32)+acd42(31)+acd42(29)+acd42(30)
      acd42(29)=acd42(3)*acd42(29)
      acd42(30)=acd42(26)*acd42(27)
      acd42(31)=acd42(23)*acd42(24)
      acd42(32)=acd42(19)*acd42(20)
      acd42(33)=acd42(16)*acd42(17)
      acd42(34)=acd42(10)*acd42(11)
      acd42(35)=acd42(9)*acd42(22)
      acd42(36)=acd42(7)*acd42(25)
      acd42(37)=acd42(5)*acd42(18)
      acd42(38)=acd42(2)*acd42(21)
      acd42(29)=acd42(29)+acd42(38)+acd42(37)+acd42(36)+acd42(35)-2.0_ki*acd42(&
      &34)+acd42(33)+acd42(32)+acd42(30)+acd42(31)
      brack(ninjaidxt0x0mu0)=acd42(29)
      brack(ninjaidxt0x0mu2)=0.0_ki
      brack(ninjaidxt0x1mu0)=acd42(28)
      brack(ninjaidxt0x2mu0)=0.0_ki
   end subroutine brack_22
!---#] subroutine brack_22:
!---#[ subroutine numerator_t2:
   subroutine numerator_t2(ncut, a0, a1, b, c, param, deg, coeffs) &
   & bind(c, name="p2_part21part21_part25part25part21_d42h0_ninja_t2")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki => ki_nin
      use p2_part21part21_part25part25part21_globalsl1, only: epspow
      use p2_part21part21_part25part25part21_kinematics
      use p2_part21part21_part25part25part21_abbrevd42h0
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
end module     p2_part21part21_part25part25part21_d42h0l132
