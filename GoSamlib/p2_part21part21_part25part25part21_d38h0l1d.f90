module     p2_part21part21_part25part25part21_d38h0l1d
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d38h0l1d.f90
   ! generator: buildfortran_d.py
   use p2_part21part21_part25part25part21_config, only: ki
   use p2_part21part21_part25part25part21_util, only: cond, d => metric_tensor
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   integer, private :: iv0
   integer, private :: iv1
   integer, private :: iv2
   integer, private :: iv3
   real(ki), dimension(4), private :: qshift
   public :: derivative
contains
!---#[ function brack_1:
   pure function brack_1(Q, mu2) result(brack)
      use p2_part21part21_part25part25part21_model
      use p2_part21part21_part25part25part21_kinematics
      use p2_part21part21_part25part25part21_color
      use p2_part21part21_part25part25part21_abbrevd38h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(55) :: acd38
      complex(ki) :: brack
      acd38(1)=dotproduct(qshift,qshift)
      acd38(2)=abb38(6)
      acd38(3)=dotproduct(qshift,spvak1k2)
      acd38(4)=abb38(8)
      acd38(5)=dotproduct(qshift,spvak2k1)
      acd38(6)=abb38(15)
      acd38(7)=dotproduct(qshift,spvak2k5)
      acd38(8)=abb38(23)
      acd38(9)=dotproduct(qshift,spvak5k2)
      acd38(10)=abb38(11)
      acd38(11)=dotproduct(qshift,spvak1e2)
      acd38(12)=dotproduct(qshift,spvae2k1)
      acd38(13)=abb38(10)
      acd38(14)=abb38(5)
      acd38(15)=abb38(18)
      acd38(16)=dotproduct(qshift,spvak2e1)
      acd38(17)=abb38(21)
      acd38(18)=dotproduct(qshift,spvae1k2)
      acd38(19)=abb38(19)
      acd38(20)=dotproduct(qshift,spvak2e5)
      acd38(21)=abb38(26)
      acd38(22)=dotproduct(qshift,spvae5k2)
      acd38(23)=abb38(24)
      acd38(24)=dotproduct(qshift,spvak5e2)
      acd38(25)=dotproduct(qshift,spvae2k5)
      acd38(26)=abb38(22)
      acd38(27)=abb38(16)
      acd38(28)=dotproduct(qshift,spvae1e2)
      acd38(29)=dotproduct(qshift,spvae2e1)
      acd38(30)=abb38(7)
      acd38(31)=abb38(17)
      acd38(32)=abb38(20)
      acd38(33)=dotproduct(qshift,spvae2e5)
      acd38(34)=dotproduct(qshift,spvae5e2)
      acd38(35)=abb38(4)
      acd38(36)=abb38(13)
      acd38(37)=abb38(14)
      acd38(38)=abb38(9)
      acd38(39)=-acd38(13)*acd38(11)
      acd38(39)=-acd38(15)+acd38(39)
      acd38(39)=acd38(12)*acd38(39)
      acd38(40)=-acd38(24)*acd38(13)
      acd38(40)=-acd38(27)+acd38(40)
      acd38(40)=acd38(25)*acd38(40)
      acd38(41)=acd38(30)*acd38(28)
      acd38(41)=-acd38(32)+acd38(41)
      acd38(41)=acd38(29)*acd38(41)
      acd38(42)=acd38(35)*acd38(33)
      acd38(42)=-acd38(37)+acd38(42)
      acd38(42)=acd38(34)*acd38(42)
      acd38(43)=acd38(2)*acd38(1)
      acd38(44)=-acd38(4)*acd38(3)
      acd38(45)=-acd38(6)*acd38(5)
      acd38(46)=-acd38(8)*acd38(7)
      acd38(47)=-acd38(10)*acd38(9)
      acd38(48)=-acd38(14)*acd38(11)
      acd38(49)=-acd38(17)*acd38(16)
      acd38(50)=-acd38(19)*acd38(18)
      acd38(51)=-acd38(21)*acd38(20)
      acd38(52)=-acd38(23)*acd38(22)
      acd38(53)=-acd38(26)*acd38(24)
      acd38(54)=-acd38(31)*acd38(28)
      acd38(55)=-acd38(36)*acd38(33)
      brack=acd38(38)+acd38(39)+acd38(40)+acd38(41)+acd38(42)+acd38(43)+acd38(4&
      &4)+acd38(45)+acd38(46)+acd38(47)+acd38(48)+acd38(49)+acd38(50)+acd38(51)&
      &+acd38(52)+acd38(53)+acd38(54)+acd38(55)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p2_part21part21_part25part25part21_model
      use p2_part21part21_part25part25part21_kinematics
      use p2_part21part21_part25part25part21_color
      use p2_part21part21_part25part25part21_abbrevd38h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(63) :: acd38
      complex(ki) :: brack
      acd38(1)=qshift(iv1)
      acd38(2)=abb38(6)
      acd38(3)=spvak1k2(iv1)
      acd38(4)=abb38(8)
      acd38(5)=spvak2k1(iv1)
      acd38(6)=abb38(15)
      acd38(7)=spvak2k5(iv1)
      acd38(8)=abb38(23)
      acd38(9)=spvak5k2(iv1)
      acd38(10)=abb38(11)
      acd38(11)=spvak1e2(iv1)
      acd38(12)=dotproduct(qshift,spvae2k1)
      acd38(13)=abb38(10)
      acd38(14)=abb38(5)
      acd38(15)=spvae2k1(iv1)
      acd38(16)=dotproduct(qshift,spvak1e2)
      acd38(17)=abb38(18)
      acd38(18)=spvak2e1(iv1)
      acd38(19)=abb38(21)
      acd38(20)=spvae1k2(iv1)
      acd38(21)=abb38(19)
      acd38(22)=spvak2e5(iv1)
      acd38(23)=abb38(26)
      acd38(24)=spvae5k2(iv1)
      acd38(25)=abb38(24)
      acd38(26)=spvak5e2(iv1)
      acd38(27)=dotproduct(qshift,spvae2k5)
      acd38(28)=abb38(22)
      acd38(29)=spvae2k5(iv1)
      acd38(30)=dotproduct(qshift,spvak5e2)
      acd38(31)=abb38(16)
      acd38(32)=spvae1e2(iv1)
      acd38(33)=dotproduct(qshift,spvae2e1)
      acd38(34)=abb38(7)
      acd38(35)=abb38(17)
      acd38(36)=spvae2e1(iv1)
      acd38(37)=dotproduct(qshift,spvae1e2)
      acd38(38)=abb38(20)
      acd38(39)=spvae2e5(iv1)
      acd38(40)=dotproduct(qshift,spvae5e2)
      acd38(41)=abb38(4)
      acd38(42)=abb38(13)
      acd38(43)=spvae5e2(iv1)
      acd38(44)=dotproduct(qshift,spvae2e5)
      acd38(45)=abb38(14)
      acd38(46)=-acd38(12)*acd38(11)
      acd38(47)=-acd38(16)*acd38(15)
      acd38(48)=-acd38(27)*acd38(26)
      acd38(49)=-acd38(30)*acd38(29)
      acd38(46)=acd38(49)+acd38(48)+acd38(47)+acd38(46)
      acd38(46)=acd38(13)*acd38(46)
      acd38(47)=acd38(33)*acd38(34)
      acd38(47)=-acd38(35)+acd38(47)
      acd38(47)=acd38(32)*acd38(47)
      acd38(48)=acd38(37)*acd38(34)
      acd38(48)=-acd38(38)+acd38(48)
      acd38(48)=acd38(36)*acd38(48)
      acd38(49)=acd38(40)*acd38(41)
      acd38(49)=-acd38(42)+acd38(49)
      acd38(49)=acd38(39)*acd38(49)
      acd38(50)=acd38(44)*acd38(41)
      acd38(50)=-acd38(45)+acd38(50)
      acd38(50)=acd38(43)*acd38(50)
      acd38(51)=acd38(2)*acd38(1)
      acd38(52)=-acd38(4)*acd38(3)
      acd38(53)=-acd38(6)*acd38(5)
      acd38(54)=-acd38(8)*acd38(7)
      acd38(55)=-acd38(10)*acd38(9)
      acd38(56)=-acd38(14)*acd38(11)
      acd38(57)=-acd38(17)*acd38(15)
      acd38(58)=-acd38(19)*acd38(18)
      acd38(59)=-acd38(21)*acd38(20)
      acd38(60)=-acd38(23)*acd38(22)
      acd38(61)=-acd38(25)*acd38(24)
      acd38(62)=-acd38(28)*acd38(26)
      acd38(63)=-acd38(31)*acd38(29)
      brack=acd38(46)+acd38(47)+acd38(48)+acd38(49)+acd38(50)+2.0_ki*acd38(51)+&
      &acd38(52)+acd38(53)+acd38(54)+acd38(55)+acd38(56)+acd38(57)+acd38(58)+ac&
      &d38(59)+acd38(60)+acd38(61)+acd38(62)+acd38(63)
   end function brack_2
!---#] function brack_2:
!---#[ function brack_3:
   pure function brack_3(Q, mu2) result(brack)
      use p2_part21part21_part25part25part21_model
      use p2_part21part21_part25part25part21_kinematics
      use p2_part21part21_part25part25part21_color
      use p2_part21part21_part25part25part21_abbrevd38h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(25) :: acd38
      complex(ki) :: brack
      acd38(1)=d(iv1,iv2)
      acd38(2)=abb38(6)
      acd38(3)=spvak1e2(iv1)
      acd38(4)=spvae2k1(iv2)
      acd38(5)=abb38(10)
      acd38(6)=spvak1e2(iv2)
      acd38(7)=spvae2k1(iv1)
      acd38(8)=spvak5e2(iv1)
      acd38(9)=spvae2k5(iv2)
      acd38(10)=spvak5e2(iv2)
      acd38(11)=spvae2k5(iv1)
      acd38(12)=spvae1e2(iv1)
      acd38(13)=spvae2e1(iv2)
      acd38(14)=abb38(7)
      acd38(15)=spvae1e2(iv2)
      acd38(16)=spvae2e1(iv1)
      acd38(17)=spvae2e5(iv1)
      acd38(18)=spvae5e2(iv2)
      acd38(19)=abb38(4)
      acd38(20)=spvae2e5(iv2)
      acd38(21)=spvae5e2(iv1)
      acd38(22)=-acd38(4)*acd38(3)
      acd38(23)=-acd38(7)*acd38(6)
      acd38(24)=-acd38(9)*acd38(8)
      acd38(25)=-acd38(11)*acd38(10)
      acd38(22)=acd38(25)+acd38(24)+acd38(23)+acd38(22)
      acd38(22)=acd38(5)*acd38(22)
      acd38(23)=acd38(13)*acd38(12)
      acd38(24)=acd38(16)*acd38(15)
      acd38(23)=acd38(24)+acd38(23)
      acd38(23)=acd38(14)*acd38(23)
      acd38(24)=acd38(18)*acd38(17)
      acd38(25)=acd38(21)*acd38(20)
      acd38(24)=acd38(25)+acd38(24)
      acd38(24)=acd38(19)*acd38(24)
      acd38(25)=acd38(2)*acd38(1)
      brack=acd38(22)+acd38(23)+acd38(24)+2.0_ki*acd38(25)
   end function brack_3
!---#] function brack_3:
!---#[ function brack_4:
   pure function brack_4(Q, mu2) result(brack)
      use p2_part21part21_part25part25part21_model
      use p2_part21part21_part25part25part21_kinematics
      use p2_part21part21_part25part25part21_color
      use p2_part21part21_part25part25part21_abbrevd38h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(1) :: acd38
      complex(ki) :: brack
      brack=0.0_ki
   end function brack_4
!---#] function brack_4:
!---#[ function derivative:
   function derivative(mu2,i1,i2,i3) result(numerator)
      use p2_part21part21_part25part25part21_globalsl1, only: epspow
      use p2_part21part21_part25part25part21_kinematics
      use p2_part21part21_part25part25part21_abbrevd38h0
      implicit none
      complex(ki), intent(in) :: mu2
      integer, intent(in), optional :: i1
      integer, intent(in), optional :: i2
      integer, intent(in), optional :: i3
      complex(ki) :: numerator
      complex(ki) :: loc
      integer :: t1
      integer :: deg
      complex(ki), dimension(4), parameter :: Q = (/ (0.0_ki,0.0_ki),(0.0_ki,0.&
      &0_ki),(0.0_ki,0.0_ki),(0.0_ki,0.0_ki)/)
      qshift = -k2
      numerator = 0.0_ki
      deg = 0
      if(present(i1)) then
          iv1=i1
          deg=1
      else
          iv1=1
      end if
      if(present(i2)) then
          iv2=i2
          deg=2
      else
          iv2=1
      end if
      if(present(i3)) then
          iv3=i3
          deg=3
      else
          iv3=1
      end if
      t1 = 0
      if(deg.eq.0) then
         numerator = cond(epspow.eq.t1,brack_1,Q,mu2)
         return
      end if
      if(deg.eq.1) then
         numerator = cond(epspow.eq.t1,brack_2,Q,mu2)
         return
      end if
      if(deg.eq.2) then
         numerator = cond(epspow.eq.t1,brack_3,Q,mu2)
         return
      end if
      if(deg.eq.3) then
         numerator = cond(epspow.eq.t1,brack_4,Q,mu2)
         return
      end if
   end function derivative
!---#] function derivative:
end module     p2_part21part21_part25part25part21_d38h0l1d
