module     p2_part21part21_part25part25part21_d117h0l131_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d117h0l131_qp.f90
   ! generator: buildfortran_tn3.py
   use p2_part21part21_part25part25part21_config, only: ki => ki_qp
   use p2_part21part21_part25part25part21_util_qp, only: cond_t, d => metric_te&
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
      use p2_part21part21_part25part25part21_model_qp
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_color_qp
      use p2_part21part21_part25part25part21_abbrevd117h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA, ninjaE3, ninjaE4
      complex(ki), intent(in) :: ninjaP
      complex(ki), dimension(1) :: acd117
      complex(ki), dimension (0:*), intent(inout) :: brack
      brack(ninjaidxt3mu0)=0.0_ki
      brack(ninjaidxt2mu0)=0.0_ki
   end subroutine brack_31
!---#] subroutine brack_31:
!---#[ subroutine brack_32:
   pure subroutine brack_32(ninjaA, ninjaE3, ninjaE4, ninjaP, brack)
      use p2_part21part21_part25part25part21_model_qp
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_color_qp
      use p2_part21part21_part25part25part21_abbrevd117h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA, ninjaE3, ninjaE4
      complex(ki), intent(in) :: ninjaP
      complex(ki), dimension(140) :: acd117
      complex(ki), dimension (0:*), intent(inout) :: brack
      acd117(1)=dotproduct(ninjaA,ninjaE3)
      acd117(2)=abb117(35)
      acd117(3)=dotproduct(ninjaE3,spval3e5)
      acd117(4)=abb117(13)
      acd117(5)=dotproduct(ninjaE3,spvak2e5)
      acd117(6)=abb117(10)
      acd117(7)=dotproduct(ninjaE3,spvae5k2)
      acd117(8)=abb117(15)
      acd117(9)=dotproduct(ninjaE3,spvae5l4)
      acd117(10)=abb117(44)
      acd117(11)=dotproduct(ninjaE3,spvae5k5)
      acd117(12)=abb117(56)
      acd117(13)=dotproduct(ninjaE3,spvae5l3)
      acd117(14)=abb117(60)
      acd117(15)=dotproduct(ninjaE3,spval4e5)
      acd117(16)=abb117(18)
      acd117(17)=dotproduct(ninjaE3,spvak5k2)
      acd117(18)=abb117(28)
      acd117(19)=dotproduct(ninjaE3,spvak5e5)
      acd117(20)=abb117(66)
      acd117(21)=dotproduct(ninjaE3,spvak2k5)
      acd117(22)=abb117(38)
      acd117(23)=dotproduct(ninjaE3,spvae1k5)
      acd117(24)=abb117(69)
      acd117(25)=dotproduct(ninjaE3,spvak1k5)
      acd117(26)=abb117(34)
      acd117(27)=dotproduct(ninjaE3,spvak5e1)
      acd117(28)=abb117(72)
      acd117(29)=dotproduct(ninjaE3,spvak5k1)
      acd117(30)=abb117(41)
      acd117(31)=abb117(37)
      acd117(32)=abb117(78)
      acd117(33)=abb117(73)
      acd117(34)=abb117(40)
      acd117(35)=abb117(80)
      acd117(36)=abb117(42)
      acd117(37)=abb117(8)
      acd117(38)=abb117(12)
      acd117(39)=abb117(16)
      acd117(40)=abb117(17)
      acd117(41)=abb117(26)
      acd117(42)=abb117(30)
      acd117(43)=abb117(31)
      acd117(44)=abb117(24)
      acd117(45)=abb117(33)
      acd117(46)=abb117(39)
      acd117(47)=abb117(65)
      acd117(48)=abb117(29)
      acd117(49)=abb117(71)
      acd117(50)=abb117(76)
      acd117(51)=abb117(36)
      acd117(52)=dotproduct(ninjaE3,spvak1e5)
      acd117(53)=dotproduct(ninjaE3,spvae5k1)
      acd117(54)=abb117(19)
      acd117(55)=dotproduct(ninjaE3,spvae1e5)
      acd117(56)=dotproduct(ninjaE3,spvae5e1)
      acd117(57)=abb117(49)
      acd117(58)=dotproduct(ninjaA,ninjaA)
      acd117(59)=dotproduct(ninjaA,spval3e5)
      acd117(60)=dotproduct(ninjaA,spvak2e5)
      acd117(61)=dotproduct(ninjaA,spvae5k2)
      acd117(62)=dotproduct(ninjaA,spvae5l4)
      acd117(63)=dotproduct(ninjaA,spvae5k5)
      acd117(64)=dotproduct(ninjaA,spvae5l3)
      acd117(65)=dotproduct(ninjaA,spval4e5)
      acd117(66)=dotproduct(ninjaA,spvak5k2)
      acd117(67)=dotproduct(ninjaA,spvak5e5)
      acd117(68)=dotproduct(ninjaA,spvak2k5)
      acd117(69)=dotproduct(ninjaA,spvae1k5)
      acd117(70)=dotproduct(ninjaA,spvak1k5)
      acd117(71)=dotproduct(ninjaA,spvak5e1)
      acd117(72)=dotproduct(ninjaA,spvak5k1)
      acd117(73)=abb117(22)
      acd117(74)=dotproduct(ninjaA,spvak1e5)
      acd117(75)=dotproduct(ninjaA,spvae5k1)
      acd117(76)=dotproduct(ninjaA,spvae1e5)
      acd117(77)=dotproduct(ninjaA,spvae5e1)
      acd117(78)=abb117(7)
      acd117(79)=abb117(9)
      acd117(80)=abb117(11)
      acd117(81)=abb117(46)
      acd117(82)=abb117(52)
      acd117(83)=abb117(51)
      acd117(84)=abb117(54)
      acd117(85)=abb117(27)
      acd117(86)=abb117(23)
      acd117(87)=abb117(20)
      acd117(88)=abb117(57)
      acd117(89)=abb117(25)
      acd117(90)=abb117(68)
      acd117(91)=abb117(77)
      acd117(92)=abb117(70)
      acd117(93)=abb117(74)
      acd117(94)=abb117(43)
      acd117(95)=abb117(48)
      acd117(96)=acd117(10)*acd117(9)
      acd117(97)=acd117(12)*acd117(11)
      acd117(98)=acd117(16)*acd117(15)
      acd117(99)=acd117(18)*acd117(17)
      acd117(100)=acd117(20)*acd117(19)
      acd117(101)=acd117(22)*acd117(21)
      acd117(102)=acd117(24)*acd117(23)
      acd117(103)=acd117(26)*acd117(25)
      acd117(104)=acd117(28)*acd117(27)
      acd117(105)=acd117(30)*acd117(29)
      acd117(96)=acd117(101)+acd117(100)+acd117(97)+acd117(96)+acd117(98)+acd11&
      &7(99)+acd117(102)+acd117(103)+acd117(104)+acd117(105)
      acd117(97)=acd117(4)*acd117(3)
      acd117(98)=acd117(6)*acd117(5)
      acd117(99)=acd117(8)*acd117(7)
      acd117(100)=acd117(14)*acd117(13)
      acd117(97)=acd117(96)+acd117(97)+acd117(98)+acd117(99)+acd117(100)
      acd117(98)=4.0_ki*acd117(2)
      acd117(98)=acd117(98)*acd117(1)
      acd117(99)=acd117(98)+acd117(97)
      acd117(97)=acd117(1)*acd117(97)
      acd117(100)=acd117(2)*acd117(1)**2
      acd117(97)=2.0_ki*acd117(100)+acd117(97)
      acd117(100)=acd117(32)*acd117(9)
      acd117(101)=acd117(33)*acd117(11)
      acd117(102)=acd117(34)*acd117(21)
      acd117(103)=acd117(35)*acd117(23)
      acd117(104)=acd117(36)*acd117(25)
      acd117(100)=acd117(104)+acd117(100)+acd117(101)+acd117(102)+acd117(103)
      acd117(101)=acd117(3)*acd117(100)
      acd117(102)=acd117(47)*acd117(15)
      acd117(103)=acd117(48)*acd117(17)
      acd117(104)=acd117(49)*acd117(19)
      acd117(105)=acd117(50)*acd117(27)
      acd117(106)=acd117(51)*acd117(29)
      acd117(102)=acd117(106)+acd117(102)+acd117(103)+acd117(104)+acd117(105)
      acd117(103)=acd117(40)*acd117(5)
      acd117(103)=acd117(103)+acd117(102)
      acd117(103)=acd117(13)*acd117(103)
      acd117(104)=acd117(38)*acd117(9)
      acd117(105)=acd117(39)*acd117(11)
      acd117(106)=acd117(41)*acd117(23)
      acd117(107)=acd117(42)*acd117(25)
      acd117(104)=acd117(104)+acd117(105)+acd117(106)+acd117(107)
      acd117(105)=acd117(5)*acd117(104)
      acd117(106)=acd117(43)*acd117(15)
      acd117(107)=acd117(44)*acd117(19)
      acd117(108)=acd117(45)*acd117(27)
      acd117(109)=acd117(46)*acd117(29)
      acd117(106)=acd117(109)+acd117(106)+acd117(107)+acd117(108)
      acd117(107)=acd117(31)*acd117(3)
      acd117(108)=acd117(37)*acd117(5)
      acd117(107)=acd117(108)+acd117(107)+acd117(106)
      acd117(107)=acd117(7)*acd117(107)
      acd117(108)=acd117(54)*acd117(53)
      acd117(109)=-acd117(52)*acd117(108)
      acd117(110)=acd117(57)*acd117(56)
      acd117(111)=acd117(55)*acd117(110)
      acd117(97)=acd117(111)+acd117(109)+2.0_ki*acd117(97)+acd117(103)+acd117(1&
      &07)+acd117(101)+acd117(105)
      acd117(96)=acd117(98)+acd117(96)
      acd117(98)=ninjaP+acd117(58)
      acd117(96)=acd117(98)*acd117(96)
      acd117(100)=acd117(59)*acd117(100)
      acd117(101)=acd117(64)*acd117(102)
      acd117(102)=acd117(60)*acd117(104)
      acd117(103)=acd117(61)*acd117(106)
      acd117(104)=2.0_ki*acd117(1)
      acd117(105)=acd117(10)*acd117(104)
      acd117(106)=acd117(32)*acd117(3)
      acd117(107)=acd117(38)*acd117(5)
      acd117(105)=acd117(107)+acd117(105)+acd117(106)
      acd117(105)=acd117(62)*acd117(105)
      acd117(106)=acd117(12)*acd117(104)
      acd117(107)=acd117(33)*acd117(3)
      acd117(109)=acd117(39)*acd117(5)
      acd117(106)=acd117(109)+acd117(106)+acd117(107)
      acd117(106)=acd117(63)*acd117(106)
      acd117(107)=acd117(16)*acd117(104)
      acd117(109)=acd117(43)*acd117(7)
      acd117(111)=acd117(47)*acd117(13)
      acd117(107)=acd117(111)+acd117(107)+acd117(109)
      acd117(107)=acd117(65)*acd117(107)
      acd117(109)=acd117(20)*acd117(104)
      acd117(111)=acd117(44)*acd117(7)
      acd117(112)=acd117(49)*acd117(13)
      acd117(109)=acd117(112)+acd117(109)+acd117(111)
      acd117(109)=acd117(67)*acd117(109)
      acd117(111)=acd117(24)*acd117(104)
      acd117(112)=acd117(35)*acd117(3)
      acd117(113)=acd117(41)*acd117(5)
      acd117(111)=acd117(113)+acd117(111)+acd117(112)
      acd117(111)=acd117(69)*acd117(111)
      acd117(112)=acd117(26)*acd117(104)
      acd117(113)=acd117(36)*acd117(3)
      acd117(114)=acd117(42)*acd117(5)
      acd117(112)=acd117(114)+acd117(112)+acd117(113)
      acd117(112)=acd117(70)*acd117(112)
      acd117(113)=acd117(28)*acd117(104)
      acd117(114)=acd117(45)*acd117(7)
      acd117(115)=acd117(50)*acd117(13)
      acd117(113)=acd117(115)+acd117(113)+acd117(114)
      acd117(113)=acd117(71)*acd117(113)
      acd117(114)=acd117(30)*acd117(104)
      acd117(115)=acd117(46)*acd117(7)
      acd117(116)=acd117(51)*acd117(13)
      acd117(114)=acd117(116)+acd117(114)+acd117(115)
      acd117(114)=acd117(72)*acd117(114)
      acd117(115)=acd117(3)*acd117(98)
      acd117(116)=acd117(59)*acd117(104)
      acd117(115)=acd117(116)+acd117(115)
      acd117(115)=acd117(4)*acd117(115)
      acd117(116)=acd117(5)*acd117(98)
      acd117(117)=acd117(60)*acd117(104)
      acd117(116)=acd117(117)+acd117(116)
      acd117(116)=acd117(6)*acd117(116)
      acd117(117)=acd117(7)*acd117(98)
      acd117(118)=acd117(61)*acd117(104)
      acd117(117)=acd117(118)+acd117(117)
      acd117(117)=acd117(8)*acd117(117)
      acd117(98)=acd117(13)*acd117(98)
      acd117(118)=acd117(64)*acd117(104)
      acd117(98)=acd117(118)+acd117(98)
      acd117(98)=acd117(14)*acd117(98)
      acd117(118)=acd117(59)*acd117(7)
      acd117(119)=acd117(61)*acd117(3)
      acd117(118)=acd117(118)+acd117(119)
      acd117(118)=acd117(31)*acd117(118)
      acd117(119)=acd117(60)*acd117(7)
      acd117(120)=acd117(61)*acd117(5)
      acd117(119)=acd117(119)+acd117(120)
      acd117(119)=acd117(37)*acd117(119)
      acd117(120)=acd117(60)*acd117(13)
      acd117(121)=acd117(64)*acd117(5)
      acd117(120)=acd117(120)+acd117(121)
      acd117(120)=acd117(40)*acd117(120)
      acd117(121)=acd117(18)*acd117(104)
      acd117(122)=acd117(48)*acd117(13)
      acd117(121)=acd117(121)+acd117(122)
      acd117(121)=acd117(66)*acd117(121)
      acd117(122)=acd117(22)*acd117(104)
      acd117(123)=acd117(34)*acd117(3)
      acd117(122)=acd117(122)+acd117(123)
      acd117(122)=acd117(68)*acd117(122)
      acd117(123)=-acd117(75)*acd117(54)
      acd117(123)=acd117(85)+acd117(123)
      acd117(123)=acd117(52)*acd117(123)
      acd117(124)=acd117(77)*acd117(57)
      acd117(124)=acd117(94)+acd117(124)
      acd117(124)=acd117(55)*acd117(124)
      acd117(104)=acd117(73)*acd117(104)
      acd117(108)=-acd117(74)*acd117(108)
      acd117(110)=acd117(76)*acd117(110)
      acd117(125)=acd117(78)*acd117(3)
      acd117(126)=acd117(79)*acd117(5)
      acd117(127)=acd117(80)*acd117(7)
      acd117(128)=acd117(81)*acd117(9)
      acd117(129)=acd117(82)*acd117(11)
      acd117(130)=acd117(83)*acd117(13)
      acd117(131)=acd117(84)*acd117(15)
      acd117(132)=acd117(86)*acd117(53)
      acd117(133)=acd117(87)*acd117(17)
      acd117(134)=acd117(88)*acd117(19)
      acd117(135)=acd117(89)*acd117(21)
      acd117(136)=acd117(90)*acd117(23)
      acd117(137)=acd117(91)*acd117(25)
      acd117(138)=acd117(92)*acd117(27)
      acd117(139)=acd117(93)*acd117(29)
      acd117(140)=acd117(95)*acd117(56)
      acd117(96)=acd117(140)+acd117(139)+acd117(138)+acd117(137)+acd117(136)+ac&
      &d117(135)+acd117(134)+acd117(133)+acd117(132)+acd117(131)+acd117(130)+ac&
      &d117(129)+acd117(128)+acd117(127)+acd117(126)+acd117(125)+acd117(110)+ac&
      &d117(108)+acd117(104)+acd117(122)+acd117(121)+acd117(114)+acd117(113)+ac&
      &d117(112)+acd117(111)+acd117(109)+acd117(107)+acd117(106)+acd117(105)+ac&
      &d117(120)+acd117(119)+acd117(118)+acd117(98)+acd117(117)+acd117(115)+acd&
      &117(116)+acd117(101)+acd117(100)+acd117(103)+acd117(102)+acd117(96)+acd1&
      &17(124)+acd117(123)
      brack(ninjaidxt1mu0)=acd117(97)
      brack(ninjaidxt1mu2)=0.0_ki
      brack(ninjaidxt0mu0)=acd117(96)
      brack(ninjaidxt0mu2)=acd117(99)
   end subroutine brack_32
!---#] subroutine brack_32:
!---#[ subroutine numerator_t3:
   subroutine numerator_t3(ncut, a, b, c, param, deg, coeffs) &
   & bind(c, name="p2_part21part21_part25part25part21_d117h0_qp_ninja_t3")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki => ki_nin
      use p2_part21part21_part25part25part21_globalsl1_qp, only: epspow
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_abbrevd117h0_qp
      implicit none
      integer(c_int), intent(in) :: ncut, deg
      complex(ki), dimension(0:3), intent(in) :: a, b, c
      complex(ki), intent(in) :: param
      complex(ki), dimension(4) :: vecA, vecB, vecC
      complex(ki), dimension(0:*), intent(out) :: coeffs
      integer :: t1
      complex(ki), dimension(4) :: qshift
      qshift = k5+k4
	     vecA(1:4) = + a(0:3) - qshift(1:4)
	     vecB(1:4) = + b(0:3)
	     vecC(1:4) = + c(0:3)
      if (deg.lt.0) return
      t1 = 0
      call cond_t(epspow.eq.t1,brack_31,vecA,vecB,vecC,param,coeffs)
      if (deg.le.(1+(0))) return
      call cond_t(epspow.eq.t1,brack_32,vecA,vecB,vecC,param,coeffs)
   end subroutine numerator_t3
!---#] subroutine numerator_t3:
end module     p2_part21part21_part25part25part21_d117h0l131_qp
