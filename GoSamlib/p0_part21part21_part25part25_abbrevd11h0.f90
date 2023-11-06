module     p0_part21part21_part25part25_abbrevd11h0
   use p0_part21part21_part25part25_config, only: ki
   use p0_part21part21_part25part25_kinematics, only: epstensor
   use p0_part21part21_part25part25_globalsh0
   implicit none
   private
   complex(ki), dimension(36), public :: abb11
   complex(ki), public :: R2d11
   public :: init_abbrev
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
contains
   subroutine     init_abbrev()
      use p0_part21part21_part25part25_config, only: deltaOS, &
     &    logfile, debug_nlo_diagrams
      use p0_part21part21_part25part25_kinematics
      use p0_part21part21_part25part25_model
      use p0_part21part21_part25part25_color, only: TR
      use p0_part21part21_part25part25_globalsl1, only: epspow
      implicit none
      abb11(1)=sqrt(mdlMT**2)
      abb11(2)=spak2l4**(-1)
      abb11(3)=spbl4k2**(-1)
      abb11(4)=mdlGC7*mdlGC45
      abb11(4)=i_*c1*abb11(4)**2
      abb11(5)=spae1e2*spbe2e1*abb11(4)
      abb11(6)=-4.0_ki/3.0_ki*abb11(5)
      abb11(7)=spak2l4*spbl4e1*spbk2e2
      abb11(8)=spak1l4*spbk1e1
      abb11(9)=spbl4e2*abb11(8)
      abb11(7)=abb11(9)+abb11(7)
      abb11(9)=abb11(1)**2*abb11(4)
      abb11(10)=abb11(9)*spae1e2
      abb11(7)=abb11(10)*abb11(7)
      abb11(11)=spbl4k2*spae1l4*spae2k2
      abb11(12)=spbl4k1*spae1k1
      abb11(13)=spae2l4*abb11(12)
      abb11(11)=abb11(13)+abb11(11)
      abb11(13)=abb11(9)*spbe2e1
      abb11(11)=abb11(13)*abb11(11)
      abb11(14)=abb11(3)*abb11(2)*mdlMh**2
      abb11(15)=abb11(14)*abb11(13)
      abb11(16)=abb11(15)+2.0_ki*abb11(13)
      abb11(16)=abb11(16)*spae2k2
      abb11(17)=spbk2k1*spae1k1
      abb11(18)=abb11(17)*abb11(16)
      abb11(19)=abb11(14)*abb11(10)
      abb11(20)=abb11(19)+2.0_ki*abb11(10)
      abb11(20)=abb11(20)*spbk2e2
      abb11(21)=spak1k2*spbk1e1
      abb11(22)=abb11(21)*abb11(20)
      abb11(5)=abb11(1)**4*abb11(5)
      abb11(5)=abb11(22)-2.0_ki*abb11(5)+abb11(18)+abb11(11)+abb11(7)
      abb11(7)=-6.0_ki*spbe2e1*abb11(10)
      abb11(11)=abb11(14)+1.0_ki
      abb11(11)=abb11(9)*abb11(11)
      abb11(18)=-spae1k2*spbk2e2*abb11(11)
      abb11(22)=spbe2k1*spae1k1
      abb11(23)=2.0_ki*abb11(9)
      abb11(24)=-abb11(23)*abb11(22)
      abb11(25)=abb11(4)*spbk2e2
      abb11(26)=abb11(25)*spak2l4
      abb11(27)=abb11(12)*abb11(26)
      abb11(28)=-spae1l4*spbl4e2*abb11(9)
      abb11(18)=abb11(27)+abb11(24)+abb11(18)+abb11(28)
      abb11(22)=-abb11(4)*abb11(22)
      abb11(24)=-spae1k2*abb11(25)
      abb11(22)=abb11(24)+abb11(22)
      abb11(11)=-spbk2e1*spae2k2*abb11(11)
      abb11(24)=spak1e2*spbk1e1
      abb11(23)=-abb11(23)*abb11(24)
      abb11(27)=abb11(4)*spae2k2
      abb11(28)=abb11(27)*spbl4k2
      abb11(29)=abb11(8)*abb11(28)
      abb11(30)=-spbl4e1*spae2l4*abb11(9)
      abb11(11)=abb11(29)+abb11(23)+abb11(11)+abb11(30)
      abb11(23)=-abb11(4)*abb11(24)
      abb11(24)=-spbk2e1*abb11(27)
      abb11(23)=abb11(24)+abb11(23)
      abb11(9)=-6.0_ki*abb11(9)
      abb11(24)=-2.0_ki*abb11(4)
      abb11(29)=-spae1l4*abb11(13)
      abb11(8)=-abb11(4)*abb11(8)
      abb11(30)=-spbl4e1*abb11(10)
      abb11(12)=-abb11(4)*abb11(12)
      abb11(31)=spae2l4*abb11(13)
      abb11(32)=spbl4e2*abb11(10)
      abb11(15)=-spae1k2*abb11(15)
      abb11(14)=abb11(4)*abb11(14)
      abb11(21)=-abb11(14)*abb11(21)
      abb11(19)=-spbk2e1*abb11(19)
      abb11(17)=-abb11(14)*abb11(17)
      abb11(33)=3.0_ki*spae1k1*abb11(13)
      abb11(34)=3.0_ki*spbk1e1*abb11(10)
      abb11(13)=spak1e2*abb11(13)
      abb11(35)=spak1k2*abb11(25)
      abb11(10)=spbe2k1*abb11(10)
      abb11(36)=spbk2k1*abb11(27)
      abb11(25)=abb11(25)*spae1k1
      abb11(27)=abb11(27)*spbk1e1
      R2d11=abb11(6)
      rat2 = rat2 + R2d11
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='11' value='", &
          & R2d11, "'/>"
      end if
   end subroutine
end module p0_part21part21_part25part25_abbrevd11h0
