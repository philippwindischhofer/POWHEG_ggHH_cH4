module     p3_part1part21_part25part25part1_d86h0l1d
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p3_part1part2 &
   ! &1_part25part25part1/helicity0d86h0l1d.f90
   ! generator: buildfortran_d.py
   use p3_part1part21_part25part25part1_config, only: ki
   use p3_part1part21_part25part25part1_util, only: cond, d => metric_tensor
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   integer, private :: iv0
   integer, private :: iv1
   integer, private :: iv2
   integer, private :: iv3
   integer, private :: iv4
   real(ki), dimension(4), private :: qshift
   public :: derivative
contains
!---#[ function brack_1:
   pure function brack_1(Q, mu2) result(brack)
      use p3_part1part21_part25part25part1_model
      use p3_part1part21_part25part25part1_kinematics
      use p3_part1part21_part25part25part1_color
      use p3_part1part21_part25part25part1_abbrevd86h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(1) :: acd86
      complex(ki) :: brack
      acd86(1)=abb86(15)
      brack=acd86(1)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p3_part1part21_part25part25part1_model
      use p3_part1part21_part25part25part1_kinematics
      use p3_part1part21_part25part25part1_color
      use p3_part1part21_part25part25part1_abbrevd86h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(36) :: acd86
      complex(ki) :: brack
      acd86(1)=spvak2k1(iv1)
      acd86(2)=abb86(34)
      acd86(3)=spval3k1(iv1)
      acd86(4)=abb86(27)
      acd86(5)=spvak5k2(iv1)
      acd86(6)=abb86(17)
      acd86(7)=spvak5l3(iv1)
      acd86(8)=abb86(8)
      acd86(9)=spvae2k1(iv1)
      acd86(10)=abb86(26)
      acd86(11)=spvak2e2(iv1)
      acd86(12)=abb86(13)
      acd86(13)=spvae2k2(iv1)
      acd86(14)=abb86(12)
      acd86(15)=spval3e2(iv1)
      acd86(16)=abb86(30)
      acd86(17)=spvae2l3(iv1)
      acd86(18)=abb86(11)
      acd86(19)=spval4e2(iv1)
      acd86(20)=abb86(19)
      acd86(21)=spvae2l4(iv1)
      acd86(22)=abb86(37)
      acd86(23)=spvak5e2(iv1)
      acd86(24)=abb86(31)
      acd86(25)=-acd86(2)*acd86(1)
      acd86(26)=-acd86(4)*acd86(3)
      acd86(27)=-acd86(6)*acd86(5)
      acd86(28)=-acd86(8)*acd86(7)
      acd86(29)=-acd86(10)*acd86(9)
      acd86(30)=-acd86(12)*acd86(11)
      acd86(31)=-acd86(14)*acd86(13)
      acd86(32)=-acd86(16)*acd86(15)
      acd86(33)=-acd86(18)*acd86(17)
      acd86(34)=-acd86(20)*acd86(19)
      acd86(35)=-acd86(22)*acd86(21)
      acd86(36)=-acd86(24)*acd86(23)
      brack=acd86(25)+acd86(26)+acd86(27)+acd86(28)+acd86(29)+acd86(30)+acd86(3&
      &1)+acd86(32)+acd86(33)+acd86(34)+acd86(35)+acd86(36)
   end function brack_2
!---#] function brack_2:
!---#[ function brack_3:
   pure function brack_3(Q, mu2) result(brack)
      use p3_part1part21_part25part25part1_model
      use p3_part1part21_part25part25part1_kinematics
      use p3_part1part21_part25part25part1_color
      use p3_part1part21_part25part25part1_abbrevd86h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(40) :: acd86
      complex(ki) :: brack
      acd86(1)=d(iv1,iv2)
      acd86(2)=abb86(24)
      acd86(3)=spvak2k1(iv1)
      acd86(4)=spvae2k2(iv2)
      acd86(5)=abb86(32)
      acd86(6)=spvak2k1(iv2)
      acd86(7)=spvae2k2(iv1)
      acd86(8)=spvak5e2(iv2)
      acd86(9)=abb86(33)
      acd86(10)=spvak5e2(iv1)
      acd86(11)=spval3k1(iv1)
      acd86(12)=spvae2l3(iv2)
      acd86(13)=abb86(28)
      acd86(14)=spval3k1(iv2)
      acd86(15)=spvae2l3(iv1)
      acd86(16)=abb86(14)
      acd86(17)=spvak5k2(iv1)
      acd86(18)=spvak2e2(iv2)
      acd86(19)=abb86(22)
      acd86(20)=spvak5k2(iv2)
      acd86(21)=spvak2e2(iv1)
      acd86(22)=spvae2k1(iv2)
      acd86(23)=abb86(16)
      acd86(24)=spvae2k1(iv1)
      acd86(25)=spvak5l3(iv1)
      acd86(26)=spval3e2(iv2)
      acd86(27)=abb86(21)
      acd86(28)=spvak5l3(iv2)
      acd86(29)=spval3e2(iv1)
      acd86(30)=abb86(20)
      acd86(31)=abb86(29)
      acd86(32)=acd86(12)*acd86(16)
      acd86(33)=acd86(4)*acd86(9)
      acd86(34)=acd86(22)*acd86(31)
      acd86(32)=acd86(34)+acd86(32)+acd86(33)
      acd86(32)=acd86(10)*acd86(32)
      acd86(33)=acd86(15)*acd86(16)
      acd86(34)=acd86(7)*acd86(9)
      acd86(35)=acd86(24)*acd86(31)
      acd86(33)=acd86(35)+acd86(33)+acd86(34)
      acd86(33)=acd86(8)*acd86(33)
      acd86(34)=acd86(29)*acd86(28)
      acd86(35)=acd86(26)*acd86(25)
      acd86(34)=acd86(35)+acd86(34)
      acd86(34)=acd86(27)*acd86(34)
      acd86(35)=acd86(21)*acd86(20)
      acd86(36)=acd86(18)*acd86(17)
      acd86(35)=acd86(36)+acd86(35)
      acd86(35)=acd86(19)*acd86(35)
      acd86(36)=acd86(15)*acd86(14)
      acd86(37)=acd86(12)*acd86(11)
      acd86(36)=acd86(37)+acd86(36)
      acd86(36)=acd86(13)*acd86(36)
      acd86(37)=acd86(7)*acd86(6)
      acd86(38)=acd86(4)*acd86(3)
      acd86(37)=acd86(38)+acd86(37)
      acd86(37)=acd86(5)*acd86(37)
      acd86(38)=acd86(26)*acd86(30)
      acd86(39)=acd86(18)*acd86(23)
      acd86(38)=acd86(38)+acd86(39)
      acd86(38)=acd86(24)*acd86(38)
      acd86(39)=acd86(29)*acd86(30)
      acd86(40)=acd86(21)*acd86(23)
      acd86(39)=acd86(39)+acd86(40)
      acd86(39)=acd86(22)*acd86(39)
      acd86(40)=acd86(1)*acd86(2)
      brack=acd86(32)+acd86(33)+acd86(34)+acd86(35)+acd86(36)+acd86(37)+acd86(3&
      &8)+acd86(39)+2.0_ki*acd86(40)
   end function brack_3
!---#] function brack_3:
!---#[ function brack_4:
   pure function brack_4(Q, mu2) result(brack)
      use p3_part1part21_part25part25part1_model
      use p3_part1part21_part25part25part1_kinematics
      use p3_part1part21_part25part25part1_color
      use p3_part1part21_part25part25part1_abbrevd86h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(45) :: acd86
      complex(ki) :: brack
      acd86(1)=d(iv1,iv2)
      acd86(2)=spvae2k1(iv3)
      acd86(3)=abb86(23)
      acd86(4)=spvak5e2(iv3)
      acd86(5)=abb86(25)
      acd86(6)=d(iv1,iv3)
      acd86(7)=spvae2k1(iv2)
      acd86(8)=spvak5e2(iv2)
      acd86(9)=d(iv2,iv3)
      acd86(10)=spvae2k1(iv1)
      acd86(11)=spvak5e2(iv1)
      acd86(12)=spvak5k2(iv2)
      acd86(13)=spvak2e2(iv3)
      acd86(14)=abb86(18)
      acd86(15)=spvak5k2(iv3)
      acd86(16)=spvak2e2(iv2)
      acd86(17)=spvak5l3(iv2)
      acd86(18)=spval3e2(iv3)
      acd86(19)=abb86(7)
      acd86(20)=spvak5l3(iv3)
      acd86(21)=spval3e2(iv2)
      acd86(22)=spvak5k2(iv1)
      acd86(23)=spvak2e2(iv1)
      acd86(24)=spvak5l3(iv1)
      acd86(25)=spval3e2(iv1)
      acd86(26)=spvak2k1(iv2)
      acd86(27)=spvae2k2(iv3)
      acd86(28)=spvak2k1(iv3)
      acd86(29)=spvae2k2(iv2)
      acd86(30)=spval3k1(iv2)
      acd86(31)=spvae2l3(iv3)
      acd86(32)=spval3k1(iv3)
      acd86(33)=spvae2l3(iv2)
      acd86(34)=spvak2k1(iv1)
      acd86(35)=spvae2k2(iv1)
      acd86(36)=spval3k1(iv1)
      acd86(37)=spvae2l3(iv1)
      acd86(38)=-acd86(32)*acd86(33)
      acd86(39)=-acd86(30)*acd86(31)
      acd86(38)=acd86(38)+acd86(39)
      acd86(38)=acd86(11)*acd86(38)
      acd86(39)=-acd86(20)*acd86(21)
      acd86(40)=-acd86(17)*acd86(18)
      acd86(39)=acd86(39)+acd86(40)
      acd86(39)=acd86(10)*acd86(39)
      acd86(40)=-acd86(32)*acd86(37)
      acd86(41)=-acd86(31)*acd86(36)
      acd86(40)=acd86(40)+acd86(41)
      acd86(40)=acd86(8)*acd86(40)
      acd86(41)=-acd86(20)*acd86(25)
      acd86(42)=-acd86(18)*acd86(24)
      acd86(41)=acd86(41)+acd86(42)
      acd86(41)=acd86(7)*acd86(41)
      acd86(42)=-acd86(33)*acd86(36)
      acd86(43)=-acd86(30)*acd86(37)
      acd86(42)=acd86(42)+acd86(43)
      acd86(42)=acd86(4)*acd86(42)
      acd86(43)=-acd86(21)*acd86(24)
      acd86(44)=-acd86(17)*acd86(25)
      acd86(43)=acd86(43)+acd86(44)
      acd86(43)=acd86(2)*acd86(43)
      acd86(38)=acd86(43)+acd86(42)+acd86(41)+acd86(40)+acd86(38)+acd86(39)
      acd86(38)=acd86(19)*acd86(38)
      acd86(39)=-acd86(28)*acd86(29)
      acd86(40)=-acd86(26)*acd86(27)
      acd86(39)=acd86(39)+acd86(40)
      acd86(39)=acd86(11)*acd86(39)
      acd86(40)=-acd86(15)*acd86(16)
      acd86(41)=-acd86(12)*acd86(13)
      acd86(40)=acd86(40)+acd86(41)
      acd86(40)=acd86(10)*acd86(40)
      acd86(41)=-acd86(28)*acd86(35)
      acd86(42)=-acd86(27)*acd86(34)
      acd86(41)=acd86(41)+acd86(42)
      acd86(41)=acd86(8)*acd86(41)
      acd86(42)=-acd86(15)*acd86(23)
      acd86(43)=-acd86(13)*acd86(22)
      acd86(42)=acd86(42)+acd86(43)
      acd86(42)=acd86(7)*acd86(42)
      acd86(43)=-acd86(29)*acd86(34)
      acd86(44)=-acd86(26)*acd86(35)
      acd86(43)=acd86(43)+acd86(44)
      acd86(43)=acd86(4)*acd86(43)
      acd86(44)=-acd86(16)*acd86(22)
      acd86(45)=-acd86(12)*acd86(23)
      acd86(44)=acd86(44)+acd86(45)
      acd86(44)=acd86(2)*acd86(44)
      acd86(39)=acd86(44)+acd86(43)+acd86(42)+acd86(41)+acd86(39)+acd86(40)
      acd86(39)=acd86(14)*acd86(39)
      acd86(40)=-acd86(11)*acd86(9)
      acd86(41)=-acd86(8)*acd86(6)
      acd86(42)=-acd86(4)*acd86(1)
      acd86(40)=acd86(42)+acd86(41)+acd86(40)
      acd86(40)=acd86(5)*acd86(40)
      acd86(41)=-acd86(10)*acd86(9)
      acd86(42)=-acd86(7)*acd86(6)
      acd86(43)=-acd86(2)*acd86(1)
      acd86(41)=acd86(43)+acd86(42)+acd86(41)
      acd86(41)=acd86(3)*acd86(41)
      acd86(40)=acd86(41)+acd86(40)
      brack=acd86(38)+acd86(39)+2.0_ki*acd86(40)
   end function brack_4
!---#] function brack_4:
!---#[ function brack_5:
   pure function brack_5(Q, mu2) result(brack)
      use p3_part1part21_part25part25part1_model
      use p3_part1part21_part25part25part1_kinematics
      use p3_part1part21_part25part25part1_color
      use p3_part1part21_part25part25part1_abbrevd86h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(21) :: acd86
      complex(ki) :: brack
      acd86(1)=d(iv1,iv2)
      acd86(2)=spvae2k1(iv3)
      acd86(3)=spvak5e2(iv4)
      acd86(4)=abb86(9)
      acd86(5)=spvae2k1(iv4)
      acd86(6)=spvak5e2(iv3)
      acd86(7)=d(iv1,iv3)
      acd86(8)=spvae2k1(iv2)
      acd86(9)=spvak5e2(iv2)
      acd86(10)=d(iv1,iv4)
      acd86(11)=d(iv2,iv3)
      acd86(12)=spvae2k1(iv1)
      acd86(13)=spvak5e2(iv1)
      acd86(14)=d(iv2,iv4)
      acd86(15)=d(iv3,iv4)
      acd86(16)=acd86(13)*acd86(11)
      acd86(17)=acd86(9)*acd86(7)
      acd86(18)=acd86(6)*acd86(1)
      acd86(16)=acd86(18)+acd86(16)+acd86(17)
      acd86(16)=acd86(5)*acd86(16)
      acd86(17)=acd86(13)*acd86(14)
      acd86(18)=acd86(9)*acd86(10)
      acd86(19)=acd86(3)*acd86(1)
      acd86(17)=acd86(19)+acd86(17)+acd86(18)
      acd86(17)=acd86(2)*acd86(17)
      acd86(18)=acd86(9)*acd86(12)
      acd86(19)=acd86(8)*acd86(13)
      acd86(18)=acd86(18)+acd86(19)
      acd86(18)=acd86(15)*acd86(18)
      acd86(19)=acd86(12)*acd86(14)
      acd86(20)=acd86(8)*acd86(10)
      acd86(19)=acd86(19)+acd86(20)
      acd86(19)=acd86(6)*acd86(19)
      acd86(20)=acd86(12)*acd86(11)
      acd86(21)=acd86(8)*acd86(7)
      acd86(20)=acd86(20)+acd86(21)
      acd86(20)=acd86(3)*acd86(20)
      acd86(16)=acd86(17)+acd86(20)+acd86(16)+acd86(19)+acd86(18)
      brack=2.0_ki*acd86(16)*acd86(4)
   end function brack_5
!---#] function brack_5:
!---#[ function derivative:
   function derivative(mu2,i1,i2,i3,i4) result(numerator)
      use p3_part1part21_part25part25part1_globalsl1, only: epspow
      use p3_part1part21_part25part25part1_kinematics
      use p3_part1part21_part25part25part1_abbrevd86h0
      implicit none
      complex(ki), intent(in) :: mu2
      integer, intent(in), optional :: i1
      integer, intent(in), optional :: i2
      integer, intent(in), optional :: i3
      integer, intent(in), optional :: i4
      complex(ki) :: numerator
      complex(ki) :: loc
      integer :: t1
      integer :: deg
      complex(ki), dimension(4), parameter :: Q = (/ (0.0_ki,0.0_ki),(0.0_ki,0.&
      &0_ki),(0.0_ki,0.0_ki),(0.0_ki,0.0_ki)/)
      qshift = 0
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
      if(present(i4)) then
          iv4=i4
          deg=4
      else
          iv4=1
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
      if(deg.eq.4) then
         numerator = cond(epspow.eq.t1,brack_5,Q,mu2)
         return
      end if
   end function derivative
!---#] function derivative:
end module     p3_part1part21_part25part25part1_d86h0l1d
