module     p10_part1part21_part25part25part1_d84h2l1d
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p10_part1part &
   ! &21_part25part25part1/helicity2d84h2l1d.f90
   ! generator: buildfortran_d.py
   use p10_part1part21_part25part25part1_config, only: ki
   use p10_part1part21_part25part25part1_util, only: cond, d => metric_tensor
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
      use p10_part1part21_part25part25part1_model
      use p10_part1part21_part25part25part1_kinematics
      use p10_part1part21_part25part25part1_color
      use p10_part1part21_part25part25part1_abbrevd84h2
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(1) :: acd84
      complex(ki) :: brack
      acd84(1)=abb84(14)
      brack=acd84(1)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p10_part1part21_part25part25part1_model
      use p10_part1part21_part25part25part1_kinematics
      use p10_part1part21_part25part25part1_color
      use p10_part1part21_part25part25part1_abbrevd84h2
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(30) :: acd84
      complex(ki) :: brack
      acd84(1)=spvak1k2(iv1)
      acd84(2)=abb84(21)
      acd84(3)=spvak2k5(iv1)
      acd84(4)=abb84(12)
      acd84(5)=spvak1e2(iv1)
      acd84(6)=abb84(9)
      acd84(7)=spvak2e2(iv1)
      acd84(8)=abb84(19)
      acd84(9)=spvae2k2(iv1)
      acd84(10)=abb84(15)
      acd84(11)=spval3e2(iv1)
      acd84(12)=abb84(40)
      acd84(13)=spvae2l3(iv1)
      acd84(14)=abb84(36)
      acd84(15)=spval4e2(iv1)
      acd84(16)=abb84(33)
      acd84(17)=spvae2l4(iv1)
      acd84(18)=abb84(29)
      acd84(19)=spvae2k5(iv1)
      acd84(20)=abb84(17)
      acd84(21)=acd84(2)*acd84(1)
      acd84(22)=acd84(4)*acd84(3)
      acd84(23)=acd84(6)*acd84(5)
      acd84(24)=acd84(8)*acd84(7)
      acd84(25)=acd84(10)*acd84(9)
      acd84(26)=acd84(12)*acd84(11)
      acd84(27)=acd84(14)*acd84(13)
      acd84(28)=acd84(16)*acd84(15)
      acd84(29)=acd84(18)*acd84(17)
      acd84(30)=acd84(20)*acd84(19)
      brack=acd84(21)+acd84(22)+acd84(23)+acd84(24)+acd84(25)+acd84(26)+acd84(2&
      &7)+acd84(28)+acd84(29)+acd84(30)
   end function brack_2
!---#] function brack_2:
!---#[ function brack_3:
   pure function brack_3(Q, mu2) result(brack)
      use p10_part1part21_part25part25part1_model
      use p10_part1part21_part25part25part1_kinematics
      use p10_part1part21_part25part25part1_color
      use p10_part1part21_part25part25part1_abbrevd84h2
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(44) :: acd84
      complex(ki) :: brack
      acd84(1)=d(iv1,iv2)
      acd84(2)=abb84(24)
      acd84(3)=spvak1k2(iv1)
      acd84(4)=spvak2e2(iv2)
      acd84(5)=abb84(27)
      acd84(6)=spval3e2(iv2)
      acd84(7)=abb84(22)
      acd84(8)=spvak1k2(iv2)
      acd84(9)=spvak2e2(iv1)
      acd84(10)=spval3e2(iv1)
      acd84(11)=spvae2k2(iv2)
      acd84(12)=abb84(7)
      acd84(13)=spvae2l3(iv2)
      acd84(14)=abb84(18)
      acd84(15)=spvae2l4(iv2)
      acd84(16)=abb84(46)
      acd84(17)=spvae2k2(iv1)
      acd84(18)=spvae2l3(iv1)
      acd84(19)=spvae2l4(iv1)
      acd84(20)=abb84(39)
      acd84(21)=abb84(23)
      acd84(22)=spvak2k5(iv1)
      acd84(23)=abb84(25)
      acd84(24)=abb84(16)
      acd84(25)=spvak2k5(iv2)
      acd84(26)=spval4e2(iv2)
      acd84(27)=abb84(44)
      acd84(28)=spval4e2(iv1)
      acd84(29)=abb84(38)
      acd84(30)=spvak1e2(iv1)
      acd84(31)=spvae2k5(iv2)
      acd84(32)=abb84(8)
      acd84(33)=spvak1e2(iv2)
      acd84(34)=spvae2k5(iv1)
      acd84(35)=acd84(15)*acd84(16)
      acd84(36)=acd84(5)*acd84(8)
      acd84(37)=acd84(13)*acd84(14)
      acd84(38)=acd84(11)*acd84(12)
      acd84(35)=acd84(38)+acd84(37)+acd84(35)+acd84(36)
      acd84(35)=acd84(9)*acd84(35)
      acd84(36)=acd84(16)*acd84(19)
      acd84(37)=acd84(3)*acd84(5)
      acd84(38)=acd84(18)*acd84(14)
      acd84(39)=acd84(17)*acd84(12)
      acd84(36)=acd84(39)+acd84(38)+acd84(36)+acd84(37)
      acd84(36)=acd84(4)*acd84(36)
      acd84(37)=acd84(26)*acd84(27)
      acd84(38)=acd84(23)*acd84(25)
      acd84(39)=acd84(6)*acd84(20)
      acd84(37)=acd84(39)+acd84(37)+acd84(38)
      acd84(37)=acd84(17)*acd84(37)
      acd84(38)=acd84(27)*acd84(28)
      acd84(39)=acd84(22)*acd84(23)
      acd84(40)=acd84(10)*acd84(20)
      acd84(38)=acd84(40)+acd84(38)+acd84(39)
      acd84(38)=acd84(11)*acd84(38)
      acd84(39)=acd84(33)*acd84(34)
      acd84(40)=acd84(30)*acd84(31)
      acd84(39)=acd84(39)+acd84(40)
      acd84(39)=acd84(32)*acd84(39)
      acd84(40)=acd84(26)*acd84(29)
      acd84(41)=acd84(24)*acd84(25)
      acd84(40)=acd84(40)+acd84(41)
      acd84(40)=acd84(18)*acd84(40)
      acd84(41)=acd84(28)*acd84(29)
      acd84(42)=acd84(22)*acd84(24)
      acd84(41)=acd84(41)+acd84(42)
      acd84(41)=acd84(13)*acd84(41)
      acd84(42)=acd84(15)*acd84(21)
      acd84(43)=acd84(7)*acd84(8)
      acd84(42)=acd84(42)+acd84(43)
      acd84(42)=acd84(10)*acd84(42)
      acd84(43)=acd84(19)*acd84(21)
      acd84(44)=acd84(3)*acd84(7)
      acd84(43)=acd84(43)+acd84(44)
      acd84(43)=acd84(6)*acd84(43)
      acd84(44)=acd84(1)*acd84(2)
      brack=acd84(35)+acd84(36)+acd84(37)+acd84(38)+acd84(39)+acd84(40)+acd84(4&
      &1)+acd84(42)+acd84(43)+2.0_ki*acd84(44)
   end function brack_3
!---#] function brack_3:
!---#[ function brack_4:
   pure function brack_4(Q, mu2) result(brack)
      use p10_part1part21_part25part25part1_model
      use p10_part1part21_part25part25part1_kinematics
      use p10_part1part21_part25part25part1_color
      use p10_part1part21_part25part25part1_abbrevd84h2
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(45) :: acd84
      complex(ki) :: brack
      acd84(1)=d(iv1,iv2)
      acd84(2)=spvak1k2(iv3)
      acd84(3)=abb84(26)
      acd84(4)=spvak2k5(iv3)
      acd84(5)=abb84(13)
      acd84(6)=spvak2e2(iv3)
      acd84(7)=abb84(20)
      acd84(8)=spvae2k2(iv3)
      acd84(9)=abb84(28)
      acd84(10)=spval3e2(iv3)
      acd84(11)=abb84(41)
      acd84(12)=spvae2l3(iv3)
      acd84(13)=abb84(37)
      acd84(14)=spval4e2(iv3)
      acd84(15)=abb84(34)
      acd84(16)=spvae2l4(iv3)
      acd84(17)=abb84(30)
      acd84(18)=d(iv1,iv3)
      acd84(19)=spvak1k2(iv2)
      acd84(20)=spvak2k5(iv2)
      acd84(21)=spvak2e2(iv2)
      acd84(22)=spvae2k2(iv2)
      acd84(23)=spval3e2(iv2)
      acd84(24)=spvae2l3(iv2)
      acd84(25)=spval4e2(iv2)
      acd84(26)=spvae2l4(iv2)
      acd84(27)=d(iv2,iv3)
      acd84(28)=spvak1k2(iv1)
      acd84(29)=spvak2k5(iv1)
      acd84(30)=spvak2e2(iv1)
      acd84(31)=spvae2k2(iv1)
      acd84(32)=spval3e2(iv1)
      acd84(33)=spvae2l3(iv1)
      acd84(34)=spval4e2(iv1)
      acd84(35)=spvae2l4(iv1)
      acd84(36)=acd84(2)*acd84(3)
      acd84(37)=acd84(4)*acd84(5)
      acd84(38)=acd84(6)*acd84(7)
      acd84(39)=acd84(8)*acd84(9)
      acd84(40)=acd84(10)*acd84(11)
      acd84(41)=acd84(12)*acd84(13)
      acd84(42)=-acd84(14)*acd84(15)
      acd84(43)=acd84(16)*acd84(17)
      acd84(36)=acd84(43)+acd84(42)+acd84(41)+acd84(40)+acd84(39)+acd84(38)+acd&
      &84(36)+acd84(37)
      acd84(36)=acd84(1)*acd84(36)
      acd84(37)=acd84(19)*acd84(3)
      acd84(38)=acd84(20)*acd84(5)
      acd84(39)=acd84(21)*acd84(7)
      acd84(40)=acd84(22)*acd84(9)
      acd84(41)=acd84(23)*acd84(11)
      acd84(42)=acd84(24)*acd84(13)
      acd84(43)=-acd84(25)*acd84(15)
      acd84(44)=acd84(26)*acd84(17)
      acd84(37)=acd84(44)+acd84(43)+acd84(42)+acd84(41)+acd84(40)+acd84(39)+acd&
      &84(38)+acd84(37)
      acd84(37)=acd84(18)*acd84(37)
      acd84(38)=acd84(28)*acd84(3)
      acd84(39)=acd84(29)*acd84(5)
      acd84(40)=acd84(30)*acd84(7)
      acd84(41)=acd84(31)*acd84(9)
      acd84(42)=acd84(32)*acd84(11)
      acd84(43)=acd84(33)*acd84(13)
      acd84(44)=-acd84(34)*acd84(15)
      acd84(45)=acd84(35)*acd84(17)
      acd84(38)=acd84(45)+acd84(44)+acd84(43)+acd84(42)+acd84(41)+acd84(40)+acd&
      &84(39)+acd84(38)
      acd84(38)=acd84(27)*acd84(38)
      acd84(36)=acd84(38)+acd84(37)+acd84(36)
      brack=2.0_ki*acd84(36)
   end function brack_4
!---#] function brack_4:
!---#[ function brack_5:
   pure function brack_5(Q, mu2) result(brack)
      use p10_part1part21_part25part25part1_model
      use p10_part1part21_part25part25part1_kinematics
      use p10_part1part21_part25part25part1_color
      use p10_part1part21_part25part25part1_abbrevd84h2
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(10) :: acd84
      complex(ki) :: brack
      acd84(1)=d(iv1,iv2)
      acd84(2)=d(iv3,iv4)
      acd84(3)=abb84(10)
      acd84(4)=d(iv1,iv3)
      acd84(5)=d(iv2,iv4)
      acd84(6)=d(iv1,iv4)
      acd84(7)=d(iv2,iv3)
      acd84(8)=acd84(2)*acd84(1)
      acd84(9)=acd84(5)*acd84(4)
      acd84(10)=acd84(7)*acd84(6)
      acd84(8)=acd84(10)+acd84(8)+acd84(9)
      brack=8.0_ki*acd84(8)*acd84(3)
   end function brack_5
!---#] function brack_5:
!---#[ function derivative:
   function derivative(mu2,i1,i2,i3,i4) result(numerator)
      use p10_part1part21_part25part25part1_globalsl1, only: epspow
      use p10_part1part21_part25part25part1_kinematics
      use p10_part1part21_part25part25part1_abbrevd84h2
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
end module     p10_part1part21_part25part25part1_d84h2l1d
