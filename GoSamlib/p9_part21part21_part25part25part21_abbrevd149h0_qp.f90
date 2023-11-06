module     p9_part21part21_part25part25part21_abbrevd149h0_qp
   use p9_part21part21_part25part25part21_config, only: ki => ki_qp
   use p9_part21part21_part25part25part21_kinematics_qp, only: epstensor
   use p9_part21part21_part25part25part21_globalsh0_qp
   implicit none
   private
   complex(ki), dimension(36), public :: abb149
   complex(ki), public :: R2d149
   public :: init_abbrev
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
contains
   subroutine     init_abbrev()
      use p9_part21part21_part25part25part21_config, only: deltaOS, &
     &    logfile, debug_nlo_diagrams
      use p9_part21part21_part25part25part21_kinematics_qp
      use p9_part21part21_part25part25part21_model_qp
      use p9_part21part21_part25part25part21_color_qp, only: TR
      use p9_part21part21_part25part25part21_globalsl1_qp, only: epspow
      implicit none
      abb149(1)=1.0_ki/(-mdlMh**2+es34)
      abb149(2)=sqrt(mdlMT**2)
      abb149(3)=spbk5e1*spae2k5
      abb149(4)=spbk2e1*spae2k2
      abb149(4)=-abb149(3)+abb149(4)
      abb149(4)=spae1e5*abb149(4)
      abb149(5)=spbk1e1*spae1e5
      abb149(6)=spak1e2*abb149(5)
      abb149(4)=abb149(6)+abb149(4)
      abb149(4)=spbe5e2*abb149(4)
      abb149(6)=spae1k5*spbk5e2
      abb149(7)=spae1k2*spbk2e2
      abb149(7)=-abb149(6)+abb149(7)
      abb149(7)=spbe5e1*abb149(7)
      abb149(8)=spae1k1*spbe5e1
      abb149(9)=spbe2k1*abb149(8)
      abb149(7)=abb149(9)+abb149(7)
      abb149(7)=spae2e5*abb149(7)
      abb149(9)=spbe2e1*spae1e5
      abb149(10)=spbk5e5*spae2k5
      abb149(11)=-abb149(10)*abb149(9)
      abb149(12)=spae1e2*spbe5e1
      abb149(13)=spae5k5*spbk5e2
      abb149(14)=-abb149(13)*abb149(12)
      abb149(4)=abb149(11)+abb149(14)+abb149(7)+abb149(4)
      abb149(7)=c1-c2
      abb149(7)=abb149(1)*abb149(7)*mdlGC7**3*mdlGC43*mdlGC45
      abb149(11)=abb149(2)**3*abb149(7)
      abb149(4)=abb149(11)*abb149(4)
      abb149(7)=abb149(2)*abb149(7)
      abb149(14)=abb149(7)*spae2k2
      abb149(15)=abb149(14)*spbk2k1
      abb149(16)=abb149(15)*abb149(8)
      abb149(17)=abb149(13)*abb149(16)
      abb149(18)=abb149(7)*spbk2e2
      abb149(19)=abb149(18)*abb149(10)
      abb149(20)=spak1k2*abb149(5)*abb149(19)
      abb149(4)=abb149(20)+abb149(17)+abb149(4)
      abb149(4)=1.0_ki/2.0_ki*abb149(4)
      abb149(17)=abb149(10)*abb149(7)
      abb149(20)=-abb149(17)*abb149(9)
      abb149(21)=abb149(13)*abb149(7)
      abb149(22)=-abb149(21)*abb149(12)
      abb149(20)=abb149(20)+abb149(22)
      abb149(20)=1.0_ki/2.0_ki*abb149(20)
      abb149(12)=-abb149(11)*abb149(12)
      abb149(22)=abb149(7)*spbk5k1*spae2k5
      abb149(23)=-abb149(8)*abb149(22)
      abb149(12)=abb149(23)+abb149(12)+abb149(16)
      abb149(12)=1.0_ki/2.0_ki*abb149(12)
      abb149(16)=1.0_ki/2.0_ki*abb149(7)
      abb149(23)=abb149(16)*spbe5e1
      abb149(24)=-spae1e2*abb149(23)
      abb149(25)=abb149(7)*spak1k5*spbk5e2
      abb149(26)=abb149(18)*spak1k2
      abb149(25)=abb149(25)-abb149(26)
      abb149(26)=-abb149(5)*abb149(25)
      abb149(9)=-abb149(11)*abb149(9)
      abb149(9)=abb149(9)+abb149(26)
      abb149(9)=1.0_ki/2.0_ki*abb149(9)
      abb149(26)=1.0_ki/2.0_ki*spae1e5
      abb149(27)=abb149(26)*abb149(7)
      abb149(28)=-spbe2e1*abb149(27)
      abb149(15)=abb149(22)-abb149(15)
      abb149(15)=spbe5e2*abb149(15)
      abb149(22)=spbe2k1*abb149(17)
      abb149(15)=abb149(22)+abb149(15)
      abb149(15)=spae1k1*abb149(15)
      abb149(22)=spae1k2*abb149(19)
      abb149(29)=abb149(11)*spbe5e2
      abb149(30)=spae1e2*abb149(29)
      abb149(15)=abb149(30)+abb149(22)+abb149(15)
      abb149(15)=1.0_ki/2.0_ki*abb149(15)
      abb149(22)=abb149(16)*spbe5e2
      abb149(30)=spae1e2*abb149(22)
      abb149(31)=spbe2k1*spae1k1
      abb149(6)=-abb149(6)+abb149(31)
      abb149(6)=abb149(7)*abb149(6)
      abb149(31)=spae1k2*abb149(18)
      abb149(6)=abb149(31)+abb149(6)
      abb149(6)=1.0_ki/2.0_ki*abb149(6)
      abb149(25)=spae2e5*abb149(25)
      abb149(31)=spak1e2*abb149(21)
      abb149(25)=abb149(31)+abb149(25)
      abb149(25)=spbk1e1*abb149(25)
      abb149(14)=abb149(14)*spbk2e1
      abb149(31)=abb149(13)*abb149(14)
      abb149(11)=abb149(11)*spae2e5
      abb149(32)=spbe2e1*abb149(11)
      abb149(25)=abb149(32)+abb149(31)+abb149(25)
      abb149(25)=1.0_ki/2.0_ki*abb149(25)
      abb149(31)=abb149(16)*spae2e5
      abb149(32)=spbe2e1*abb149(31)
      abb149(33)=spak1e2*spbk1e1
      abb149(3)=-abb149(3)+abb149(33)
      abb149(3)=abb149(7)*abb149(3)
      abb149(3)=abb149(14)+abb149(3)
      abb149(3)=1.0_ki/2.0_ki*abb149(3)
      abb149(14)=-spae1e5*abb149(29)
      abb149(11)=-spbe5e1*abb149(11)
      abb149(29)=spae2k5*abb149(23)
      abb149(33)=-spae2k5*abb149(22)
      abb149(27)=spbk5e2*abb149(27)
      abb149(34)=-spbk5e2*abb149(31)
      abb149(23)=abb149(23)*spae2k2
      abb149(35)=-abb149(13)*abb149(23)
      abb149(36)=spae2k2*abb149(22)
      abb149(19)=-abb149(26)*abb149(19)
      abb149(18)=abb149(18)*abb149(26)
      abb149(26)=spbk2e2*abb149(31)
      abb149(8)=abb149(16)*abb149(8)
      abb149(13)=-abb149(13)*abb149(8)
      abb149(22)=spae1k1*abb149(22)
      abb149(5)=abb149(16)*abb149(5)
      abb149(10)=-abb149(10)*abb149(5)
      abb149(16)=spbk1e1*abb149(31)
      R2d149=0.0_ki
      rat2 = rat2 + R2d149
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='149' value='", &
          & R2d149, "'/>"
      end if
   end subroutine
end module p9_part21part21_part25part25part21_abbrevd149h0_qp
