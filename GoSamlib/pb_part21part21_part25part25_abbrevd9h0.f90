module     pb_part21part21_part25part25_abbrevd9h0
   use pb_part21part21_part25part25_config, only: ki
   use pb_part21part21_part25part25_kinematics, only: epstensor
   use pb_part21part21_part25part25_globalsh0
   implicit none
   private
   complex(ki), dimension(36), public :: abb9
   complex(ki), public :: R2d9
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
      abb9(1)=sqrt(mdlMT**2)
      abb9(2)=spak2l3**(-1)
      abb9(3)=spbl3k2**(-1)
      abb9(4)=mdlGC7*mdlGC45
      abb9(4)=i_*c1*abb9(4)**2
      abb9(5)=spae1e2*spbe2e1*abb9(4)
      abb9(6)=-4.0_ki/3.0_ki*abb9(5)
      abb9(7)=spak2l3*spbl3e1*spbk2e2
      abb9(8)=spak1l3*spbk1e1
      abb9(9)=spbl3e2*abb9(8)
      abb9(7)=abb9(9)+abb9(7)
      abb9(9)=abb9(1)**2*abb9(4)
      abb9(10)=abb9(9)*spae1e2
      abb9(7)=abb9(10)*abb9(7)
      abb9(11)=spbl3k2*spae1l3*spae2k2
      abb9(12)=spbl3k1*spae1k1
      abb9(13)=spae2l3*abb9(12)
      abb9(11)=abb9(13)+abb9(11)
      abb9(13)=abb9(9)*spbe2e1
      abb9(11)=abb9(13)*abb9(11)
      abb9(14)=abb9(3)*abb9(2)*mdlMh**2
      abb9(15)=abb9(14)*abb9(13)
      abb9(16)=abb9(15)+2.0_ki*abb9(13)
      abb9(16)=abb9(16)*spae2k2
      abb9(17)=spbk2k1*spae1k1
      abb9(18)=abb9(17)*abb9(16)
      abb9(19)=abb9(14)*abb9(10)
      abb9(20)=abb9(19)+2.0_ki*abb9(10)
      abb9(20)=abb9(20)*spbk2e2
      abb9(21)=spak1k2*spbk1e1
      abb9(22)=abb9(21)*abb9(20)
      abb9(5)=abb9(1)**4*abb9(5)
      abb9(5)=abb9(22)-2.0_ki*abb9(5)+abb9(18)+abb9(11)+abb9(7)
      abb9(7)=-6.0_ki*spbe2e1*abb9(10)
      abb9(11)=abb9(14)+1.0_ki
      abb9(11)=abb9(9)*abb9(11)
      abb9(18)=-spae1k2*spbk2e2*abb9(11)
      abb9(22)=spbe2k1*spae1k1
      abb9(23)=2.0_ki*abb9(9)
      abb9(24)=-abb9(23)*abb9(22)
      abb9(25)=abb9(4)*spbk2e2
      abb9(26)=abb9(25)*spak2l3
      abb9(27)=abb9(12)*abb9(26)
      abb9(28)=-spae1l3*spbl3e2*abb9(9)
      abb9(18)=abb9(27)+abb9(24)+abb9(18)+abb9(28)
      abb9(22)=-abb9(4)*abb9(22)
      abb9(24)=-spae1k2*abb9(25)
      abb9(22)=abb9(24)+abb9(22)
      abb9(11)=-spbk2e1*spae2k2*abb9(11)
      abb9(24)=spak1e2*spbk1e1
      abb9(23)=-abb9(23)*abb9(24)
      abb9(27)=abb9(4)*spae2k2
      abb9(28)=abb9(27)*spbl3k2
      abb9(29)=abb9(8)*abb9(28)
      abb9(30)=-spbl3e1*spae2l3*abb9(9)
      abb9(11)=abb9(29)+abb9(23)+abb9(11)+abb9(30)
      abb9(23)=-abb9(4)*abb9(24)
      abb9(24)=-spbk2e1*abb9(27)
      abb9(23)=abb9(24)+abb9(23)
      abb9(9)=-6.0_ki*abb9(9)
      abb9(24)=-2.0_ki*abb9(4)
      abb9(29)=-spae1l3*abb9(13)
      abb9(8)=-abb9(4)*abb9(8)
      abb9(30)=-spbl3e1*abb9(10)
      abb9(12)=-abb9(4)*abb9(12)
      abb9(31)=spae2l3*abb9(13)
      abb9(32)=spbl3e2*abb9(10)
      abb9(15)=-spae1k2*abb9(15)
      abb9(14)=abb9(4)*abb9(14)
      abb9(21)=-abb9(14)*abb9(21)
      abb9(19)=-spbk2e1*abb9(19)
      abb9(17)=-abb9(14)*abb9(17)
      abb9(33)=3.0_ki*spae1k1*abb9(13)
      abb9(34)=3.0_ki*spbk1e1*abb9(10)
      abb9(13)=spak1e2*abb9(13)
      abb9(35)=spak1k2*abb9(25)
      abb9(10)=spbe2k1*abb9(10)
      abb9(36)=spbk2k1*abb9(27)
      abb9(25)=abb9(25)*spae1k1
      abb9(27)=abb9(27)*spbk1e1
      R2d9=abb9(6)
      rat2 = rat2 + R2d9
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='9' value='", &
          & R2d9, "'/>"
      end if
   end subroutine
end module pb_part21part21_part25part25_abbrevd9h0
