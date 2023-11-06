module     p2_part21part21_part25part25part21_abbrevd52h0
   use p2_part21part21_part25part25part21_config, only: ki
   use p2_part21part21_part25part25part21_kinematics, only: epstensor
   use p2_part21part21_part25part25part21_globalsh0
   implicit none
   private
   complex(ki), dimension(38), public :: abb52
   complex(ki), public :: R2d52
   public :: init_abbrev
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
contains
   subroutine     init_abbrev()
      use p2_part21part21_part25part25part21_config, only: deltaOS, &
     &    logfile, debug_nlo_diagrams
      use p2_part21part21_part25part25part21_kinematics
      use p2_part21part21_part25part25part21_model
      use p2_part21part21_part25part25part21_color, only: TR
      use p2_part21part21_part25part25part21_globalsl1, only: epspow
      implicit none
      abb52(1)=sqrt(mdlMT**2)
      abb52(2)=spbk1e1*spak1e5
      abb52(3)=abb52(2)*spae1e2
      abb52(4)=abb52(3)*spbe5e2
      abb52(5)=spae1k1*spbe5k1
      abb52(6)=abb52(5)*spbe2e1
      abb52(7)=abb52(6)*spae2e5
      abb52(4)=abb52(4)+abb52(7)
      abb52(7)=spbk2e2*spae1e2
      abb52(8)=spak2e5*abb52(7)
      abb52(9)=spbe2k1*spak1e5
      abb52(10)=spae1e2*abb52(9)
      abb52(8)=abb52(10)+abb52(8)
      abb52(8)=spbe5e1*abb52(8)
      abb52(10)=spak1e2*spbe5k1
      abb52(11)=spae1e5*spbe2e1
      abb52(12)=abb52(11)*abb52(10)
      abb52(13)=spae5k5*spbe5e1
      abb52(14)=-spbk5e2*spae1e2*abb52(13)
      abb52(15)=spbk5e5*spae1e5
      abb52(16)=-spae2k5*spbe2e1*abb52(15)
      abb52(17)=spae2k2*spbe2e1
      abb52(18)=spbe5k2*spae1e5*abb52(17)
      abb52(8)=abb52(12)+abb52(18)+abb52(16)+abb52(14)+abb52(8)+abb52(4)
      abb52(12)=c1-c2
      abb52(12)=abb52(12)*mdlGC7**3*i_*mdlGC46
      abb52(14)=-abb52(1)**3*abb52(12)
      abb52(8)=abb52(14)*abb52(8)
      abb52(16)=spbk2k1*spak1e5
      abb52(18)=-spbk5k2*spae5k5
      abb52(18)=abb52(16)+abb52(18)
      abb52(18)=abb52(18)*abb52(6)*spae2k2
      abb52(19)=spak1k2*spbe5k1
      abb52(20)=-spak2k5*spbk5e5
      abb52(20)=abb52(19)+abb52(20)
      abb52(20)=abb52(20)*abb52(3)*spbk2e2
      abb52(18)=abb52(20)+abb52(18)
      abb52(12)=abb52(1)*abb52(12)
      abb52(18)=abb52(12)*abb52(18)
      abb52(8)=abb52(18)+abb52(8)
      abb52(8)=1.0_ki/2.0_ki*abb52(8)
      abb52(4)=-1.0_ki/2.0_ki*abb52(12)*abb52(4)
      abb52(18)=abb52(14)*spbe5e1
      abb52(20)=-spae1e2*abb52(18)
      abb52(11)=-abb52(14)*abb52(11)
      abb52(21)=abb52(12)*spbk5e5
      abb52(22)=-spak2k5*abb52(21)
      abb52(23)=abb52(12)*spbe5k1
      abb52(24)=spak1k2*abb52(23)
      abb52(22)=abb52(24)+abb52(22)
      abb52(22)=abb52(7)*abb52(22)
      abb52(24)=abb52(14)*spbe5e2
      abb52(25)=spae1e2*abb52(24)
      abb52(22)=abb52(25)+abb52(22)
      abb52(22)=1.0_ki/2.0_ki*abb52(22)
      abb52(25)=1.0_ki/2.0_ki*spae1e2
      abb52(26)=-abb52(25)*abb52(12)*spbe5e2
      abb52(27)=abb52(12)*spae5k5
      abb52(28)=-spbk5k2*abb52(27)
      abb52(29)=abb52(12)*spak1e5
      abb52(30)=spbk2k1*abb52(29)
      abb52(28)=abb52(30)+abb52(28)
      abb52(28)=abb52(17)*abb52(28)
      abb52(14)=abb52(14)*spae2e5*spbe2e1
      abb52(14)=abb52(14)+abb52(28)
      abb52(14)=1.0_ki/2.0_ki*abb52(14)
      abb52(28)=1.0_ki/2.0_ki*spbe2e1
      abb52(30)=-abb52(28)*abb52(12)*spae2e5
      abb52(31)=spak2e5*spbk2e2
      abb52(32)=spbk5e2*spae5k5
      abb52(9)=-abb52(9)-abb52(31)+abb52(32)
      abb52(9)=abb52(5)*abb52(9)
      abb52(15)=spak2k5*abb52(15)
      abb52(19)=-abb52(19)*spae1e5
      abb52(15)=abb52(19)+abb52(15)
      abb52(15)=spbk2e2*abb52(15)
      abb52(9)=abb52(9)+abb52(15)
      abb52(9)=abb52(12)*abb52(9)
      abb52(15)=-spae1e5*abb52(24)
      abb52(9)=abb52(15)+abb52(9)
      abb52(9)=1.0_ki/2.0_ki*abb52(9)
      abb52(15)=1.0_ki/2.0_ki*abb52(12)
      abb52(19)=abb52(15)*spae1e5
      abb52(24)=spbe5e2*abb52(19)
      abb52(5)=abb52(12)*abb52(5)
      abb52(31)=-abb52(12)*abb52(31)
      abb52(32)=spbk5e2*abb52(27)
      abb52(33)=-spbe2k1*abb52(29)
      abb52(31)=abb52(33)+abb52(32)+abb52(31)
      abb52(31)=1.0_ki/2.0_ki*abb52(31)
      abb52(32)=spbe5k2*spae2k2
      abb52(33)=spae2k5*spbk5e5
      abb52(10)=-abb52(10)-abb52(32)+abb52(33)
      abb52(10)=abb52(2)*abb52(10)
      abb52(13)=spbk5k2*abb52(13)
      abb52(16)=-spbe5e1*abb52(16)
      abb52(13)=abb52(16)+abb52(13)
      abb52(13)=spae2k2*abb52(13)
      abb52(10)=abb52(10)+abb52(13)
      abb52(10)=abb52(12)*abb52(10)
      abb52(13)=-spae2e5*abb52(18)
      abb52(10)=abb52(13)+abb52(10)
      abb52(10)=1.0_ki/2.0_ki*abb52(10)
      abb52(13)=abb52(15)*spbe5e1
      abb52(16)=spae2e5*abb52(13)
      abb52(2)=abb52(12)*abb52(2)
      abb52(18)=-abb52(12)*abb52(32)
      abb52(32)=spae2k5*abb52(21)
      abb52(33)=-spak1e2*abb52(23)
      abb52(18)=abb52(33)+abb52(32)+abb52(18)
      abb52(18)=1.0_ki/2.0_ki*abb52(18)
      abb52(6)=abb52(15)*abb52(6)
      abb52(32)=spae5k5*abb52(6)
      abb52(27)=abb52(28)*abb52(27)
      abb52(33)=-spae5k5*abb52(13)
      abb52(3)=abb52(15)*abb52(3)
      abb52(34)=spbk5e5*abb52(3)
      abb52(21)=abb52(25)*abb52(21)
      abb52(35)=-spbk5e5*abb52(19)
      abb52(36)=-spae2k2*abb52(6)
      abb52(17)=-abb52(17)*abb52(15)
      abb52(37)=spae2k2*abb52(13)
      abb52(38)=-spbk2e2*abb52(3)
      abb52(7)=-abb52(7)*abb52(15)
      abb52(15)=spbk2e2*abb52(19)
      abb52(6)=spak1e5*abb52(6)
      abb52(28)=abb52(28)*abb52(29)
      abb52(13)=-spak1e5*abb52(13)
      abb52(3)=spbe5k1*abb52(3)
      abb52(23)=abb52(25)*abb52(23)
      abb52(19)=-spbe5k1*abb52(19)
      R2d52=0.0_ki
      rat2 = rat2 + R2d52
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='52' value='", &
          & R2d52, "'/>"
      end if
   end subroutine
end module p2_part21part21_part25part25part21_abbrevd52h0
