module     pb_part21part21_part25part25_abbrevd13h0
   use pb_part21part21_part25part25_config, only: ki
   use pb_part21part21_part25part25_kinematics, only: epstensor
   use pb_part21part21_part25part25_globalsh0
   implicit none
   private
   complex(ki), dimension(46), public :: abb13
   complex(ki), public :: R2d13
   public :: init_abbrev
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
contains
   subroutine     init_abbrev()
      use pb_part21part21_part25part25_config, only: deltaOS, &
     &    logfile, debug_nlo_diagrams
      use pb_part21part21_part25part25_kinematics
      use pb_part21part21_part25part25_model
      use pb_part21part21_part25part25_color, only: TR
      use pb_part21part21_part25part25_globalsl1, only: epspow
      implicit none
      abb13(1)=sqrt(mdlMT**2)
      abb13(2)=spak2l3**(-1)
      abb13(3)=spbl3k2**(-1)
      abb13(4)=spak2l4**(-1)
      abb13(5)=spbl4k2**(-1)
      abb13(6)=mdlGC45*mdlGC7
      abb13(6)=c1*i_*abb13(6)**2
      abb13(7)=spbe2e1*abb13(6)
      abb13(8)=abb13(7)*spae1e2
      abb13(9)=5.0_ki/3.0_ki*abb13(8)
      abb13(10)=abb13(4)*abb13(5)
      abb13(11)=abb13(10)*spae2k2
      abb13(12)=abb13(11)*spbk2k1
      abb13(13)=spbe2e1*spae1k1
      abb13(14)=abb13(12)*abb13(13)
      abb13(15)=spae1e2*spbk1e1
      abb13(10)=abb13(10)*spbk2e2
      abb13(16)=abb13(15)*abb13(10)*spak1k2
      abb13(14)=abb13(14)+abb13(16)
      abb13(16)=spbl4e2*spae1l4
      abb13(17)=spae1k1*spbe2k1
      abb13(16)=abb13(16)+2.0_ki*abb13(17)
      abb13(17)=spae2k2*spbk2e1*abb13(16)
      abb13(18)=spae2l4*spbl4e1
      abb13(19)=spbk1e1*spak1e2
      abb13(18)=abb13(18)+2.0_ki*abb13(19)
      abb13(19)=spbk2e2*spae1k2*abb13(18)
      abb13(17)=abb13(17)+abb13(19)
      abb13(19)=abb13(2)*abb13(3)
      abb13(17)=abb13(17)*abb13(19)
      abb13(20)=spbl3e2*spae1l3
      abb13(21)=spbk2e1*abb13(11)*abb13(20)
      abb13(22)=spae2l3*spbl3e1
      abb13(23)=spae1k2*abb13(10)*abb13(22)
      abb13(17)=abb13(17)+abb13(21)+abb13(23)+abb13(14)
      abb13(21)=mdlMh**2
      abb13(17)=abb13(17)*abb13(21)
      abb13(23)=spae2l4*spbl4k1
      abb13(24)=abb13(13)*abb13(23)
      abb13(25)=spbl4e2*spak1l4
      abb13(26)=abb13(15)*abb13(25)
      abb13(24)=abb13(24)+abb13(26)
      abb13(26)=abb13(20)*abb13(18)
      abb13(27)=abb13(22)*abb13(16)
      abb13(28)=abb13(1)**2
      abb13(29)=abb13(28)*spbe2e1
      abb13(30)=abb13(29)*spae1e2
      abb13(17)=-2.0_ki*abb13(30)+abb13(17)+abb13(27)+abb13(26)+abb13(24)
      abb13(26)=abb13(28)*abb13(6)
      abb13(17)=abb13(17)*abb13(26)
      abb13(14)=abb13(14)*abb13(21)
      abb13(14)=-4.0_ki*abb13(30)+abb13(14)+abb13(24)
      abb13(14)=abb13(14)*abb13(6)
      abb13(8)=-2.0_ki*abb13(8)
      abb13(24)=abb13(19)*spbk2e2
      abb13(27)=abb13(24)-abb13(10)
      abb13(30)=abb13(21)*spae1k2*abb13(27)
      abb13(16)=abb13(30)+abb13(20)-abb13(16)
      abb13(20)=2.0_ki*abb13(26)
      abb13(16)=abb13(16)*abb13(20)
      abb13(19)=abb13(19)*spae2k2
      abb13(30)=abb13(19)-abb13(11)
      abb13(31)=abb13(21)*spbk2e1*abb13(30)
      abb13(18)=abb13(31)+abb13(22)-abb13(18)
      abb13(18)=abb13(18)*abb13(20)
      abb13(20)=-8.0_ki*abb13(26)
      abb13(22)=abb13(26)*spbe2e1
      abb13(31)=spae2l4*abb13(22)
      abb13(32)=spae2l4*abb13(7)
      abb13(33)=abb13(26)*spae1e2
      abb13(34)=spbl4e2*abb13(33)
      abb13(35)=spae1e2*abb13(6)
      abb13(36)=spbl4e2*abb13(35)
      abb13(25)=spbk1e1*abb13(25)
      abb13(37)=-abb13(29)+abb13(25)
      abb13(37)=spae2l3*abb13(37)
      abb13(10)=abb13(10)*spae2l3
      abb13(38)=spbk1e1*spak1k2*abb13(21)*abb13(10)
      abb13(37)=abb13(38)+abb13(37)
      abb13(37)=abb13(37)*abb13(6)
      abb13(38)=-spae2l3*abb13(7)
      abb13(39)=abb13(6)*spae2l3
      abb13(40)=spbl4e2*abb13(39)
      abb13(12)=spae1k1*abb13(21)*abb13(12)
      abb13(23)=spae1k1*abb13(23)
      abb13(28)=abb13(28)*spae1e2
      abb13(12)=-abb13(28)+abb13(23)+abb13(12)
      abb13(41)=abb13(6)*spbl3e2
      abb13(12)=abb13(12)*abb13(41)
      abb13(42)=-spbl3e2*abb13(35)
      abb13(43)=spae2l4*abb13(41)
      abb13(25)=abb13(19)*abb13(25)
      abb13(29)=-abb13(30)*abb13(29)
      abb13(25)=abb13(25)+abb13(29)
      abb13(6)=abb13(6)*abb13(21)
      abb13(21)=abb13(25)*abb13(6)
      abb13(25)=-spbe2e1*abb13(30)*abb13(6)
      abb13(19)=abb13(19)*abb13(6)
      abb13(29)=spbl4e2*abb13(19)
      abb13(11)=spbl3e2*abb13(6)*abb13(11)
      abb13(23)=abb13(24)*abb13(23)
      abb13(28)=-abb13(27)*abb13(28)
      abb13(23)=abb13(23)+abb13(28)
      abb13(23)=abb13(23)*abb13(6)
      abb13(27)=-spae1e2*abb13(27)*abb13(6)
      abb13(24)=abb13(24)*abb13(6)
      abb13(28)=spae2l4*abb13(24)
      abb13(6)=abb13(6)*abb13(10)
      abb13(10)=abb13(26)*abb13(13)
      abb13(13)=spae1k1*abb13(7)
      abb13(30)=spae1k1*abb13(41)
      abb13(44)=spae1k1*abb13(24)
      abb13(15)=abb13(26)*abb13(15)
      abb13(26)=spbk1e1*abb13(35)
      abb13(45)=spbk1e1*abb13(39)
      abb13(46)=spbk1e1*abb13(19)
      abb13(22)=spak1e2*abb13(22)
      abb13(7)=spak1e2*abb13(7)
      abb13(41)=spak1e2*abb13(41)
      abb13(24)=spak1e2*abb13(24)
      abb13(33)=spbe2k1*abb13(33)
      abb13(35)=spbe2k1*abb13(35)
      abb13(39)=spbe2k1*abb13(39)
      abb13(19)=spbe2k1*abb13(19)
      R2d13=abb13(9)
      rat2 = rat2 + R2d13
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='13' value='", &
          & R2d13, "'/>"
      end if
   end subroutine
end module pb_part21part21_part25part25_abbrevd13h0
