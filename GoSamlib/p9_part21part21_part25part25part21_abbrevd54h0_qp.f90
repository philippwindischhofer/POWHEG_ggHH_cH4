module     p9_part21part21_part25part25part21_abbrevd54h0_qp
   use p9_part21part21_part25part25part21_config, only: ki => ki_qp
   use p9_part21part21_part25part25part21_kinematics_qp, only: epstensor
   use p9_part21part21_part25part25part21_globalsh0_qp
   implicit none
   private
   complex(ki), dimension(35), public :: abb54
   complex(ki), public :: R2d54
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
      abb54(1)=sqrt(mdlMT**2)
      abb54(2)=spbk5e5*spae2e5
      abb54(3)=-spae1k5*abb54(2)
      abb54(4)=spae1k2*spbe5k2
      abb54(5)=-spae2e5*abb54(4)
      abb54(3)=abb54(5)+abb54(3)
      abb54(3)=spbe2e1*abb54(3)
      abb54(5)=spbk5e1*spae5k5
      abb54(6)=spbk2e1*spak2e5
      abb54(7)=-abb54(6)-abb54(5)
      abb54(8)=spbe5e2*spae1e2
      abb54(7)=abb54(8)*abb54(7)
      abb54(9)=spae1e5*spbe2e1
      abb54(10)=spae2k2*spbe5k2
      abb54(11)=abb54(10)*abb54(9)
      abb54(12)=spbk2e2*spak2e5
      abb54(13)=abb54(12)*spbe5e1*spae1e2
      abb54(14)=spae1k1*spbe2e1
      abb54(15)=spbe5k1*spae2e5*abb54(14)
      abb54(16)=spbk1e1*spae1e2
      abb54(17)=spak1e5*spbe5e2*abb54(16)
      abb54(3)=abb54(17)+abb54(15)+abb54(13)+abb54(11)+abb54(7)+abb54(3)
      abb54(7)=c1-c2
      abb54(7)=abb54(7)*mdlGC7**3*i_*mdlGC46
      abb54(11)=-abb54(1)**3*abb54(7)
      abb54(3)=abb54(11)*abb54(3)
      abb54(7)=abb54(1)*abb54(7)
      abb54(15)=abb54(7)*spae5k5
      abb54(17)=spbk5k1*abb54(15)*abb54(14)
      abb54(18)=-abb54(10)*abb54(17)
      abb54(19)=abb54(7)*spbk5e5
      abb54(20)=spak1k5*abb54(19)*abb54(16)
      abb54(21)=-abb54(12)*abb54(20)
      abb54(3)=abb54(21)+abb54(18)+abb54(3)
      abb54(3)=1.0_ki/2.0_ki*abb54(3)
      abb54(13)=-abb54(7)*abb54(13)
      abb54(18)=abb54(7)*abb54(10)
      abb54(21)=-abb54(9)*abb54(18)
      abb54(13)=abb54(13)+abb54(21)
      abb54(13)=1.0_ki/2.0_ki*abb54(13)
      abb54(21)=spak1k2*spbe5k2
      abb54(22)=-abb54(21)*abb54(16)*abb54(7)
      abb54(23)=abb54(11)*spbe5e1
      abb54(24)=spae1e2*abb54(23)
      abb54(20)=abb54(22)+abb54(24)-abb54(20)
      abb54(20)=1.0_ki/2.0_ki*abb54(20)
      abb54(22)=1.0_ki/2.0_ki*spae1e2
      abb54(24)=abb54(22)*abb54(7)
      abb54(25)=-spbe5e1*abb54(24)
      abb54(26)=spbk2k1*spak2e5
      abb54(27)=-abb54(26)*abb54(7)*abb54(14)
      abb54(28)=abb54(11)*abb54(9)
      abb54(17)=abb54(27)+abb54(28)-abb54(17)
      abb54(17)=1.0_ki/2.0_ki*abb54(17)
      abb54(27)=1.0_ki/2.0_ki*abb54(7)
      abb54(9)=-abb54(9)*abb54(27)
      abb54(8)=-abb54(11)*abb54(8)
      abb54(28)=-spbe2e1*abb54(11)*spae2e5
      abb54(29)=spbe5k1*spae1k1
      abb54(30)=-spae1k5*spbk5e5
      abb54(30)=abb54(29)+abb54(30)
      abb54(30)=abb54(12)*abb54(30)
      abb54(31)=-spbk5k1*spae5k5
      abb54(26)=-abb54(26)+abb54(31)
      abb54(26)=abb54(26)*spae1k1*spbe5e2
      abb54(26)=abb54(26)+abb54(30)
      abb54(26)=abb54(7)*abb54(26)
      abb54(11)=spae1e5*spbe5e2*abb54(11)
      abb54(11)=abb54(11)+abb54(26)
      abb54(11)=1.0_ki/2.0_ki*abb54(11)
      abb54(26)=abb54(27)*spbe5e2
      abb54(30)=-spae1e5*abb54(26)
      abb54(4)=-abb54(4)+abb54(29)
      abb54(4)=abb54(7)*abb54(4)
      abb54(29)=-spae1k5*abb54(19)
      abb54(4)=abb54(29)+abb54(4)
      abb54(4)=1.0_ki/2.0_ki*abb54(4)
      abb54(29)=abb54(7)*abb54(12)
      abb54(2)=-spak1k5*abb54(2)
      abb54(21)=-abb54(21)*spae2e5
      abb54(2)=abb54(21)+abb54(2)
      abb54(2)=spbk1e1*abb54(2)
      abb54(21)=spak1e5*spbk1e1
      abb54(5)=abb54(21)-abb54(5)
      abb54(5)=abb54(10)*abb54(5)
      abb54(2)=abb54(5)+abb54(2)
      abb54(2)=abb54(7)*abb54(2)
      abb54(5)=spae2e5*abb54(23)
      abb54(2)=abb54(5)+abb54(2)
      abb54(2)=1.0_ki/2.0_ki*abb54(2)
      abb54(5)=abb54(27)*spae2e5
      abb54(23)=-spbe5e1*abb54(5)
      abb54(6)=-abb54(6)+abb54(21)
      abb54(6)=abb54(7)*abb54(6)
      abb54(15)=-spbk5e1*abb54(15)
      abb54(6)=abb54(15)+abb54(6)
      abb54(6)=1.0_ki/2.0_ki*abb54(6)
      abb54(15)=abb54(27)*spbe2e1
      abb54(21)=spae5k5*abb54(15)
      abb54(31)=abb54(10)*abb54(21)
      abb54(32)=spae5k5*abb54(26)
      abb54(19)=abb54(19)*abb54(22)
      abb54(22)=abb54(12)*abb54(19)
      abb54(33)=spbk5e5*abb54(5)
      abb54(15)=spak2e5*abb54(15)
      abb54(34)=spak2e5*abb54(26)
      abb54(24)=spbe5k2*abb54(24)
      abb54(35)=spbe5k2*abb54(5)
      abb54(14)=abb54(27)*abb54(14)
      abb54(10)=-abb54(10)*abb54(14)
      abb54(26)=-spae1k1*abb54(26)
      abb54(16)=abb54(27)*abb54(16)
      abb54(12)=-abb54(12)*abb54(16)
      abb54(5)=-spbk1e1*abb54(5)
      R2d54=0.0_ki
      rat2 = rat2 + R2d54
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='54' value='", &
          & R2d54, "'/>"
      end if
   end subroutine
end module p9_part21part21_part25part25part21_abbrevd54h0_qp
