module     p10_part1part21_part25part25part1_abbrevd86h2_qp
   use p10_part1part21_part25part25part1_config, only: ki => ki_qp
   use p10_part1part21_part25part25part1_kinematics_qp, only: epstensor
   use p10_part1part21_part25part25part1_globalsh2_qp
   implicit none
   private
   complex(ki), dimension(44), public :: abb86
   complex(ki), public :: R2d86
   public :: init_abbrev
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
contains
   subroutine     init_abbrev()
      use p10_part1part21_part25part25part1_config, only: deltaOS, &
     &    logfile, debug_nlo_diagrams
      use p10_part1part21_part25part25part1_kinematics_qp
      use p10_part1part21_part25part25part1_model_qp
      use p10_part1part21_part25part25part1_color_qp, only: TR
      use p10_part1part21_part25part25part1_globalsl1_qp, only: epspow
      implicit none
      abb86(1)=es51**(-1)
      abb86(2)=sqrt(mdlMT**2)
      abb86(3)=spak2l3**(-1)
      abb86(4)=spbl3k2**(-1)
      abb86(5)=spak2l4**(-1)
      abb86(6)=spbl4k2**(-1)
      abb86(7)=mdlGC7**3*abb86(1)*c1
      abb86(8)=mdlGC45**2
      abb86(9)=abb86(7)*abb86(8)*TR
      abb86(10)=-8.0_ki/3.0_ki*spbk5e2*spak1e2*abb86(9)
      abb86(11)=spbk5l4*spak2l4
      abb86(12)=spbk2e2*spak1e2
      abb86(13)=abb86(12)*abb86(11)
      abb86(14)=spbk5e2*spak1l4
      abb86(15)=abb86(14)*spae2k2*spbl4k2
      abb86(13)=abb86(15)+abb86(13)
      abb86(13)=abb86(4)*abb86(13)
      abb86(15)=spae2k2*spak1l3
      abb86(16)=spbk5e2*abb86(15)
      abb86(13)=abb86(16)+abb86(13)
      abb86(13)=abb86(3)*abb86(13)
      abb86(16)=abb86(5)*abb86(6)
      abb86(17)=abb86(16)*spak2l3
      abb86(17)=abb86(17)+abb86(4)
      abb86(18)=spak1e2*spbl3e2
      abb86(19)=-spbk5k2*abb86(18)*abb86(17)
      abb86(20)=abb86(4)*spbk2e2
      abb86(21)=abb86(20)*spbk5l3
      abb86(22)=spak1e2*abb86(21)
      abb86(23)=abb86(16)*spbl3k2
      abb86(23)=abb86(23)+abb86(3)
      abb86(24)=-spak1k2*spbk5e2*spae2l3*abb86(23)
      abb86(13)=abb86(19)+abb86(24)+abb86(22)+abb86(13)
      abb86(19)=mdlMh**2
      abb86(13)=abb86(13)*abb86(19)
      abb86(22)=spbk5l4*spal3l4
      abb86(24)=abb86(18)*abb86(22)
      abb86(25)=spae2l3*spbk5l3
      abb86(26)=spbk5l4*spae2l4
      abb86(27)=abb86(26)+2.0_ki*abb86(25)
      abb86(28)=spak1k2*spbk2e2
      abb86(29)=-abb86(27)*abb86(28)
      abb86(30)=spbl3e2*spak1l3
      abb86(31)=spak1l4*spbl4e2
      abb86(32)=abb86(31)+2.0_ki*abb86(30)
      abb86(33)=spbk5k2*spae2k2
      abb86(34)=-abb86(32)*abb86(33)
      abb86(35)=abb86(2)**2
      abb86(36)=abb86(35)*spbk5e2
      abb86(37)=abb86(36)*spak1e2
      abb86(38)=spae2l3*spbl4l3*abb86(14)
      abb86(13)=abb86(13)-2.0_ki*abb86(37)+abb86(34)+abb86(29)+abb86(24)+abb86(&
      &38)
      abb86(24)=2.0_ki*TR
      abb86(7)=abb86(7)*abb86(24)*abb86(8)
      abb86(8)=abb86(7)*abb86(35)
      abb86(13)=abb86(13)*abb86(8)
      abb86(24)=abb86(25)+abb86(26)
      abb86(25)=-abb86(24)*abb86(28)
      abb86(26)=abb86(30)+abb86(31)
      abb86(29)=-abb86(26)*abb86(33)
      abb86(25)=-6.0_ki*abb86(37)+abb86(25)+abb86(29)
      abb86(25)=abb86(25)*abb86(7)
      abb86(29)=abb86(3)*abb86(4)
      abb86(30)=abb86(29)*spbk2e2
      abb86(31)=abb86(16)*spbk2e2
      abb86(34)=abb86(31)+2.0_ki*abb86(30)
      abb86(37)=abb86(19)*spak1k2
      abb86(34)=abb86(34)*abb86(37)
      abb86(32)=abb86(34)+3.0_ki*abb86(28)+abb86(32)
      abb86(32)=abb86(32)*abb86(8)
      abb86(34)=abb86(31)+abb86(30)
      abb86(34)=abb86(34)*abb86(37)
      abb86(26)=abb86(34)+abb86(28)+abb86(26)
      abb86(26)=abb86(26)*abb86(7)
      abb86(34)=abb86(35)*abb86(9)
      abb86(38)=4.0_ki*abb86(34)
      abb86(14)=-abb86(14)*abb86(38)
      abb86(38)=-spak1e2*spbk5l4*abb86(38)
      abb86(39)=-abb86(28)*abb86(22)
      abb86(31)=spak2l3*abb86(31)
      abb86(20)=abb86(31)+abb86(20)
      abb86(31)=abb86(19)*spbk5k2
      abb86(40)=abb86(31)*spak1k2
      abb86(20)=abb86(20)*abb86(40)
      abb86(41)=spak1l3*abb86(36)
      abb86(20)=abb86(20)+abb86(39)-3.0_ki*abb86(41)
      abb86(20)=abb86(20)*abb86(7)
      abb86(39)=spak1l4*spbl4l3
      abb86(41)=-abb86(33)*abb86(39)
      abb86(42)=abb86(16)*spae2k2
      abb86(43)=spbl3k2*abb86(42)
      abb86(44)=abb86(3)*spae2k2
      abb86(43)=abb86(43)+abb86(44)
      abb86(40)=abb86(43)*abb86(40)
      abb86(43)=spbk5l3*abb86(35)*spak1e2
      abb86(40)=abb86(40)+abb86(41)-3.0_ki*abb86(43)
      abb86(40)=abb86(40)*abb86(7)
      abb86(23)=-abb86(23)*abb86(37)
      abb86(23)=abb86(39)+abb86(23)
      abb86(23)=abb86(23)*abb86(7)
      abb86(19)=abb86(7)*abb86(19)
      abb86(30)=-abb86(30)*abb86(11)
      abb86(16)=2.0_ki*abb86(16)+3.0_ki*abb86(29)
      abb86(16)=abb86(16)*abb86(35)
      abb86(37)=-spbk5e2*abb86(16)
      abb86(21)=abb86(37)-abb86(21)+abb86(30)
      abb86(21)=abb86(19)*spak1k2*abb86(21)
      abb86(30)=spbl4k2*abb86(4)*spak1l4
      abb86(37)=-spae2k2*abb86(30)
      abb86(15)=-abb86(15)+abb86(37)
      abb86(15)=abb86(3)*abb86(15)
      abb86(16)=-spak1e2*abb86(16)
      abb86(15)=abb86(15)+abb86(16)
      abb86(15)=abb86(19)*spbk5k2*abb86(15)
      abb86(16)=spak1l3+abb86(30)
      abb86(16)=abb86(3)*abb86(16)*abb86(19)
      abb86(30)=abb86(29)*spae2k2
      abb86(37)=-abb86(42)-2.0_ki*abb86(30)
      abb86(37)=abb86(37)*abb86(31)
      abb86(27)=abb86(37)-3.0_ki*abb86(33)-abb86(27)
      abb86(27)=abb86(27)*abb86(8)
      abb86(30)=-abb86(42)-abb86(30)
      abb86(30)=abb86(30)*abb86(31)
      abb86(24)=abb86(30)-abb86(33)-abb86(24)
      abb86(24)=abb86(24)*abb86(7)
      abb86(30)=12.0_ki*abb86(34)
      abb86(9)=4.0_ki*abb86(9)
      abb86(17)=abb86(17)*abb86(31)
      abb86(17)=-abb86(22)+abb86(17)
      abb86(17)=abb86(17)*abb86(7)
      abb86(11)=-abb86(29)*abb86(11)
      abb86(22)=abb86(4)*spbk5l3
      abb86(11)=-abb86(22)+abb86(11)
      abb86(11)=abb86(11)*abb86(19)
      abb86(8)=-abb86(18)*abb86(8)
      abb86(18)=abb86(7)*abb86(28)
      abb86(19)=abb86(29)*abb86(19)
      abb86(12)=-abb86(35)*abb86(19)*abb86(12)
      abb86(22)=abb86(28)*abb86(19)
      abb86(28)=spae2l3*abb86(7)*abb86(36)
      abb86(29)=-abb86(7)*abb86(33)
      abb86(31)=spae2k2*abb86(36)*abb86(19)
      abb86(33)=-abb86(33)*abb86(19)
      R2d86=abb86(10)
      rat2 = rat2 + R2d86
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='86' value='", &
          & R2d86, "'/>"
      end if
   end subroutine
end module p10_part1part21_part25part25part1_abbrevd86h2_qp
