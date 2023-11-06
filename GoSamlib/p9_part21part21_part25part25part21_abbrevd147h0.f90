module     p9_part21part21_part25part25part21_abbrevd147h0
   use p9_part21part21_part25part25part21_config, only: ki
   use p9_part21part21_part25part25part21_kinematics, only: epstensor
   use p9_part21part21_part25part25part21_globalsh0
   implicit none
   private
   complex(ki), dimension(36), public :: abb147
   complex(ki), public :: R2d147
   public :: init_abbrev
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
contains
   subroutine     init_abbrev()
      use p9_part21part21_part25part25part21_config, only: deltaOS, &
     &    logfile, debug_nlo_diagrams
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_model
      use p9_part21part21_part25part25part21_color, only: TR
      use p9_part21part21_part25part25part21_globalsl1, only: epspow
      implicit none
      abb147(1)=1.0_ki/(-mdlMh**2+es34)
      abb147(2)=sqrt(mdlMT**2)
      abb147(3)=spbk2e1*spak2e5
      abb147(4)=-spbk5e1*spae5k5
      abb147(4)=-abb147(3)+abb147(4)
      abb147(4)=spae1e2*abb147(4)
      abb147(5)=spbk1e1*spae1e2
      abb147(6)=spak1e5*abb147(5)
      abb147(4)=abb147(6)+abb147(4)
      abb147(4)=spbe5e2*abb147(4)
      abb147(6)=spae1k2*spbe5k2
      abb147(7)=-spae1k5*spbk5e5
      abb147(7)=-abb147(6)+abb147(7)
      abb147(7)=spbe2e1*abb147(7)
      abb147(8)=spae1k1*spbe2e1
      abb147(9)=spbe5k1*abb147(8)
      abb147(7)=abb147(9)+abb147(7)
      abb147(7)=spae2e5*abb147(7)
      abb147(9)=spbe5e1*spae1e2
      abb147(10)=spbk2e2*spak2e5
      abb147(11)=abb147(10)*abb147(9)
      abb147(12)=spae1e5*spbe2e1
      abb147(13)=spae2k2*spbe5k2
      abb147(14)=abb147(13)*abb147(12)
      abb147(4)=abb147(11)+abb147(14)+abb147(7)+abb147(4)
      abb147(7)=c1-c2
      abb147(7)=abb147(1)*abb147(7)*mdlGC7**3*mdlGC43*mdlGC45
      abb147(11)=abb147(2)**3*abb147(7)
      abb147(4)=abb147(11)*abb147(4)
      abb147(7)=abb147(2)*abb147(7)
      abb147(14)=abb147(7)*spae5k5
      abb147(15)=abb147(14)*spbk5k1
      abb147(16)=abb147(15)*abb147(8)
      abb147(17)=abb147(13)*abb147(16)
      abb147(18)=abb147(7)*spbk5e5
      abb147(19)=abb147(18)*abb147(10)
      abb147(20)=spak1k5*abb147(5)*abb147(19)
      abb147(4)=abb147(20)+abb147(17)+abb147(4)
      abb147(4)=1.0_ki/2.0_ki*abb147(4)
      abb147(17)=abb147(10)*abb147(7)
      abb147(20)=abb147(17)*abb147(9)
      abb147(21)=abb147(13)*abb147(7)
      abb147(22)=abb147(21)*abb147(12)
      abb147(20)=abb147(20)+abb147(22)
      abb147(20)=1.0_ki/2.0_ki*abb147(20)
      abb147(22)=abb147(7)*spak1k2*spbe5k2
      abb147(23)=abb147(18)*spak1k5
      abb147(22)=abb147(22)+abb147(23)
      abb147(23)=abb147(5)*abb147(22)
      abb147(9)=abb147(11)*abb147(9)
      abb147(9)=abb147(9)+abb147(23)
      abb147(9)=1.0_ki/2.0_ki*abb147(9)
      abb147(23)=1.0_ki/2.0_ki*spae1e2
      abb147(24)=abb147(23)*abb147(7)
      abb147(25)=spbe5e1*abb147(24)
      abb147(12)=abb147(11)*abb147(12)
      abb147(26)=abb147(7)*spbk2k1*spak2e5
      abb147(27)=abb147(8)*abb147(26)
      abb147(12)=abb147(27)+abb147(12)+abb147(16)
      abb147(12)=1.0_ki/2.0_ki*abb147(12)
      abb147(16)=1.0_ki/2.0_ki*abb147(7)
      abb147(27)=abb147(16)*spbe2e1
      abb147(28)=spae1e5*abb147(27)
      abb147(29)=abb147(11)*spbe5e2
      abb147(30)=-spae1e2*abb147(29)
      abb147(11)=abb147(11)*spae2e5
      abb147(31)=-spbe2e1*abb147(11)
      abb147(15)=abb147(26)+abb147(15)
      abb147(15)=spbe5e2*abb147(15)
      abb147(26)=-spbe5k1*abb147(17)
      abb147(15)=abb147(26)+abb147(15)
      abb147(15)=spae1k1*abb147(15)
      abb147(26)=spae1k5*abb147(19)
      abb147(29)=spae1e5*abb147(29)
      abb147(15)=abb147(29)+abb147(26)+abb147(15)
      abb147(15)=1.0_ki/2.0_ki*abb147(15)
      abb147(26)=abb147(16)*spbe5e2
      abb147(29)=spae1e5*abb147(26)
      abb147(32)=-spbe5k1*spae1k1
      abb147(6)=abb147(6)+abb147(32)
      abb147(6)=abb147(7)*abb147(6)
      abb147(32)=spae1k5*abb147(18)
      abb147(6)=abb147(32)+abb147(6)
      abb147(6)=1.0_ki/2.0_ki*abb147(6)
      abb147(22)=spae2e5*abb147(22)
      abb147(32)=-spak1e5*abb147(21)
      abb147(22)=abb147(32)+abb147(22)
      abb147(22)=spbk1e1*abb147(22)
      abb147(14)=abb147(14)*spbk5e1
      abb147(32)=abb147(13)*abb147(14)
      abb147(11)=spbe5e1*abb147(11)
      abb147(11)=abb147(11)+abb147(32)+abb147(22)
      abb147(11)=1.0_ki/2.0_ki*abb147(11)
      abb147(22)=abb147(16)*spae2e5
      abb147(32)=spbe5e1*abb147(22)
      abb147(33)=-spak1e5*spbk1e1
      abb147(3)=abb147(3)+abb147(33)
      abb147(3)=abb147(7)*abb147(3)
      abb147(3)=abb147(14)+abb147(3)
      abb147(3)=1.0_ki/2.0_ki*abb147(3)
      abb147(14)=abb147(27)*spae5k5
      abb147(33)=-abb147(13)*abb147(14)
      abb147(34)=-spae5k5*abb147(26)
      abb147(19)=-abb147(23)*abb147(19)
      abb147(18)=abb147(18)*abb147(23)
      abb147(23)=-spbk5e5*abb147(22)
      abb147(27)=-spak2e5*abb147(27)
      abb147(35)=-spak2e5*abb147(26)
      abb147(24)=-spbe5k2*abb147(24)
      abb147(36)=-spbe5k2*abb147(22)
      abb147(8)=abb147(16)*abb147(8)
      abb147(13)=abb147(13)*abb147(8)
      abb147(26)=spae1k1*abb147(26)
      abb147(5)=abb147(16)*abb147(5)
      abb147(10)=abb147(10)*abb147(5)
      abb147(16)=spbk1e1*abb147(22)
      R2d147=0.0_ki
      rat2 = rat2 + R2d147
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='147' value='", &
          & R2d147, "'/>"
      end if
   end subroutine
end module p9_part21part21_part25part25part21_abbrevd147h0
