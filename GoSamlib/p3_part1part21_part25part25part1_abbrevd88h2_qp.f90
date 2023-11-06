module     p3_part1part21_part25part25part1_abbrevd88h2_qp
   use p3_part1part21_part25part25part1_config, only: ki => ki_qp
   use p3_part1part21_part25part25part1_kinematics_qp, only: epstensor
   use p3_part1part21_part25part25part1_globalsh2_qp
   implicit none
   private
   complex(ki), dimension(41), public :: abb88
   complex(ki), public :: R2d88
   public :: init_abbrev
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
contains
   subroutine     init_abbrev()
      use p3_part1part21_part25part25part1_config, only: deltaOS, &
     &    logfile, debug_nlo_diagrams
      use p3_part1part21_part25part25part1_kinematics_qp
      use p3_part1part21_part25part25part1_model_qp
      use p3_part1part21_part25part25part1_color_qp, only: TR
      use p3_part1part21_part25part25part1_globalsl1_qp, only: epspow
      implicit none
      abb88(1)=es51**(-1)
      abb88(2)=sqrt(mdlMT**2)
      abb88(3)=spak2l3**(-1)
      abb88(4)=spbl3k2**(-1)
      abb88(5)=spbl4k2**(-1)
      abb88(6)=spak2l4**(-1)
      abb88(7)=mdlGC7**3*abb88(1)*c1
      abb88(8)=mdlGC45**2
      abb88(9)=abb88(7)*abb88(8)*TR
      abb88(10)=-8.0_ki/3.0_ki*spbk5e2*spak1e2*abb88(9)
      abb88(11)=abb88(6)*spae2k2
      abb88(12)=abb88(11)*spak1l4
      abb88(13)=abb88(3)*abb88(4)
      abb88(14)=abb88(13)*spbl4k2
      abb88(14)=abb88(14)+abb88(6)
      abb88(15)=-spak1k2*abb88(14)*spae2l4
      abb88(15)=abb88(15)+abb88(12)
      abb88(15)=spbk5e2*abb88(15)
      abb88(16)=spbk5l3*spak2l3
      abb88(17)=spbk2e2*spak1e2
      abb88(18)=abb88(17)*abb88(16)
      abb88(19)=spbk5e2*spak1l3
      abb88(20)=spae2k2*spbl3k2*abb88(19)
      abb88(18)=abb88(20)+abb88(18)
      abb88(18)=abb88(6)*abb88(18)
      abb88(20)=spbk5l4*abb88(17)
      abb88(18)=abb88(20)+abb88(18)
      abb88(18)=abb88(5)*abb88(18)
      abb88(20)=abb88(13)*spak2l4
      abb88(20)=abb88(20)+abb88(5)
      abb88(21)=spak1e2*spbl4e2
      abb88(22)=-spbk5k2*abb88(21)*abb88(20)
      abb88(15)=abb88(22)+abb88(18)+abb88(15)
      abb88(18)=mdlMh**2
      abb88(15)=abb88(15)*abb88(18)
      abb88(22)=spbk5l3*spal3l4
      abb88(23)=-abb88(21)*abb88(22)
      abb88(24)=spae2l4*spbk5l4
      abb88(25)=spbk5l3*spae2l3
      abb88(26)=abb88(25)+2.0_ki*abb88(24)
      abb88(27)=spak1k2*spbk2e2
      abb88(28)=-abb88(26)*abb88(27)
      abb88(29)=spbl4e2*spak1l4
      abb88(30)=spak1l3*spbl3e2
      abb88(31)=abb88(30)+2.0_ki*abb88(29)
      abb88(32)=spbk5k2*spae2k2
      abb88(33)=-abb88(31)*abb88(32)
      abb88(34)=abb88(2)**2
      abb88(35)=abb88(34)*spbk5e2
      abb88(36)=abb88(35)*spak1e2
      abb88(37)=-spae2l4*spbl4l3*abb88(19)
      abb88(15)=abb88(15)-2.0_ki*abb88(36)+abb88(33)+abb88(28)+abb88(23)+abb88(&
      &37)
      abb88(23)=2.0_ki*TR
      abb88(7)=abb88(7)*abb88(23)*abb88(8)
      abb88(8)=abb88(7)*abb88(34)
      abb88(15)=abb88(15)*abb88(8)
      abb88(23)=abb88(24)+abb88(25)
      abb88(24)=-abb88(23)*abb88(27)
      abb88(25)=abb88(29)+abb88(30)
      abb88(28)=-abb88(25)*abb88(32)
      abb88(24)=-6.0_ki*abb88(36)+abb88(24)+abb88(28)
      abb88(24)=abb88(24)*abb88(7)
      abb88(28)=abb88(5)*abb88(6)
      abb88(29)=abb88(28)*spbk2e2
      abb88(30)=abb88(13)*spbk2e2
      abb88(33)=abb88(30)+2.0_ki*abb88(29)
      abb88(36)=abb88(18)*spak1k2
      abb88(33)=abb88(33)*abb88(36)
      abb88(31)=abb88(33)+3.0_ki*abb88(27)+abb88(31)
      abb88(31)=abb88(31)*abb88(8)
      abb88(29)=abb88(30)+abb88(29)
      abb88(29)=abb88(29)*abb88(36)
      abb88(25)=abb88(29)+abb88(27)+abb88(25)
      abb88(25)=abb88(25)*abb88(7)
      abb88(29)=abb88(27)*abb88(22)
      abb88(33)=abb88(5)*spbk2e2
      abb88(30)=spak2l4*abb88(30)
      abb88(30)=abb88(30)+abb88(33)
      abb88(33)=abb88(18)*spbk5k2
      abb88(37)=abb88(33)*spak1k2
      abb88(30)=abb88(30)*abb88(37)
      abb88(38)=spak1l4*abb88(35)
      abb88(29)=abb88(30)+abb88(29)-3.0_ki*abb88(38)
      abb88(29)=abb88(29)*abb88(7)
      abb88(30)=spak1l3*spbl4l3
      abb88(38)=abb88(32)*abb88(30)
      abb88(39)=abb88(13)*spae2k2
      abb88(40)=spbl4k2*abb88(39)
      abb88(11)=abb88(40)+abb88(11)
      abb88(11)=abb88(11)*abb88(37)
      abb88(37)=spbk5l4*abb88(34)*spak1e2
      abb88(11)=abb88(11)+abb88(38)-3.0_ki*abb88(37)
      abb88(11)=abb88(11)*abb88(7)
      abb88(14)=-abb88(14)*abb88(36)
      abb88(14)=-abb88(30)+abb88(14)
      abb88(14)=abb88(14)*abb88(7)
      abb88(30)=abb88(34)*abb88(9)
      abb88(36)=4.0_ki*abb88(30)
      abb88(19)=-abb88(19)*abb88(36)
      abb88(36)=-spak1e2*spbk5l3*abb88(36)
      abb88(18)=abb88(7)*abb88(18)
      abb88(13)=2.0_ki*abb88(13)+3.0_ki*abb88(28)
      abb88(13)=abb88(13)*abb88(34)
      abb88(37)=-spbk5e2*abb88(13)
      abb88(16)=abb88(16)*abb88(6)
      abb88(16)=abb88(16)+spbk5l4
      abb88(16)=abb88(16)*abb88(5)
      abb88(38)=-spbk2e2*abb88(16)
      abb88(37)=abb88(38)+abb88(37)
      abb88(37)=abb88(18)*spak1k2*abb88(37)
      abb88(38)=spak1l3*spbl3k2
      abb88(40)=abb88(28)*spae2k2
      abb88(41)=-abb88(40)*abb88(38)
      abb88(13)=-spak1e2*abb88(13)
      abb88(12)=abb88(13)-abb88(12)+abb88(41)
      abb88(12)=abb88(18)*spbk5k2*abb88(12)
      abb88(13)=abb88(28)*abb88(38)
      abb88(38)=abb88(6)*spak1l4
      abb88(13)=abb88(38)+abb88(13)
      abb88(13)=abb88(13)*abb88(18)
      abb88(38)=-abb88(39)-2.0_ki*abb88(40)
      abb88(38)=abb88(38)*abb88(33)
      abb88(26)=abb88(38)-3.0_ki*abb88(32)-abb88(26)
      abb88(26)=abb88(26)*abb88(8)
      abb88(38)=-abb88(39)-abb88(40)
      abb88(38)=abb88(38)*abb88(33)
      abb88(23)=abb88(38)-abb88(32)-abb88(23)
      abb88(23)=abb88(23)*abb88(7)
      abb88(30)=12.0_ki*abb88(30)
      abb88(9)=4.0_ki*abb88(9)
      abb88(20)=abb88(20)*abb88(33)
      abb88(20)=abb88(22)+abb88(20)
      abb88(20)=abb88(20)*abb88(7)
      abb88(16)=-abb88(18)*abb88(16)
      abb88(8)=-abb88(21)*abb88(8)
      abb88(21)=abb88(7)*abb88(27)
      abb88(18)=abb88(28)*abb88(18)
      abb88(17)=-abb88(34)*abb88(18)*abb88(17)
      abb88(22)=abb88(27)*abb88(18)
      abb88(27)=spae2l4*abb88(7)*abb88(35)
      abb88(28)=-abb88(7)*abb88(32)
      abb88(33)=spae2k2*abb88(35)*abb88(18)
      abb88(32)=-abb88(32)*abb88(18)
      R2d88=abb88(10)
      rat2 = rat2 + R2d88
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='88' value='", &
          & R2d88, "'/>"
      end if
   end subroutine
end module p3_part1part21_part25part25part1_abbrevd88h2_qp
