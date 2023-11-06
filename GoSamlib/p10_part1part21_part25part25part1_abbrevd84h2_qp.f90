module     p10_part1part21_part25part25part1_abbrevd84h2_qp
   use p10_part1part21_part25part25part1_config, only: ki => ki_qp
   use p10_part1part21_part25part25part1_kinematics_qp, only: epstensor
   use p10_part1part21_part25part25part1_globalsh2_qp
   implicit none
   private
   complex(ki), dimension(46), public :: abb84
   complex(ki), public :: R2d84
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
      abb84(1)=es51**(-1)
      abb84(2)=sqrt(mdlMT**2)
      abb84(3)=spak2l3**(-1)
      abb84(4)=spbl3k2**(-1)
      abb84(5)=spak2l4**(-1)
      abb84(6)=spbl4k2**(-1)
      abb84(7)=mdlGC45**2
      abb84(8)=abb84(7)*TR
      abb84(9)=mdlGC7**3*abb84(1)*c1
      abb84(10)=spak1e2*abb84(9)*spbk5e2
      abb84(11)=10.0_ki/3.0_ki*abb84(8)*abb84(10)
      abb84(12)=spbk5l3*spae2l3
      abb84(13)=abb84(5)*abb84(6)
      abb84(14)=spbk2e2*abb84(13)*abb84(12)
      abb84(15)=spbk2e2*spbk5l4
      abb84(16)=abb84(3)*abb84(4)
      abb84(17)=spae2l4*abb84(16)*abb84(15)
      abb84(14)=abb84(14)+abb84(17)
      abb84(14)=spak1k2*abb84(14)
      abb84(17)=spak1l3*spbl3e2
      abb84(18)=spae2k2*abb84(13)*abb84(17)
      abb84(19)=spae2k2*spak1l4
      abb84(20)=spbl4e2*abb84(16)*abb84(19)
      abb84(18)=abb84(18)+abb84(20)
      abb84(18)=spbk5k2*abb84(18)
      abb84(14)=abb84(14)+abb84(18)
      abb84(18)=mdlMh**2
      abb84(14)=abb84(14)*abb84(18)
      abb84(15)=abb84(15)*spak2l4
      abb84(20)=abb84(2)**2
      abb84(21)=abb84(20)*spbk5e2
      abb84(22)=-2.0_ki*abb84(21)-abb84(15)
      abb84(22)=spak1e2*abb84(22)
      abb84(23)=spak1l4*spbl4e2
      abb84(24)=spak1k2*spbk2e2
      abb84(23)=-abb84(23)+2.0_ki*abb84(24)
      abb84(24)=-abb84(12)*abb84(23)
      abb84(25)=spbk5l4*spae2l4
      abb84(26)=spbk5k2*spae2k2
      abb84(25)=-abb84(25)+2.0_ki*abb84(26)
      abb84(26)=-abb84(17)*abb84(25)
      abb84(19)=abb84(19)*spbl4k2
      abb84(27)=abb84(19)*spbk5e2
      abb84(14)=abb84(14)-abb84(27)+abb84(26)+abb84(24)+abb84(22)
      abb84(7)=abb84(9)*abb84(7)
      abb84(22)=2.0_ki*TR
      abb84(22)=abb84(22)*abb84(7)
      abb84(14)=abb84(14)*abb84(20)*abb84(22)
      abb84(24)=-4.0_ki*abb84(21)-abb84(15)
      abb84(24)=spak1e2*abb84(24)
      abb84(24)=-abb84(27)+abb84(24)
      abb84(24)=abb84(24)*abb84(22)
      abb84(26)=4.0_ki*abb84(8)
      abb84(10)=-abb84(26)*abb84(10)
      abb84(8)=abb84(8)*abb84(20)*abb84(9)
      abb84(9)=4.0_ki*abb84(8)
      abb84(26)=abb84(13)-abb84(16)
      abb84(27)=abb84(18)*spak1k2
      abb84(28)=abb84(26)*abb84(27)
      abb84(29)=spbk2e2*abb84(28)
      abb84(17)=abb84(29)-abb84(17)-abb84(23)
      abb84(17)=abb84(17)*abb84(9)
      abb84(23)=abb84(22)*spak1l4
      abb84(29)=-abb84(21)*abb84(23)
      abb84(30)=-spbk5e2*abb84(23)
      abb84(31)=abb84(22)*spak1e2
      abb84(32)=abb84(31)*abb84(20)
      abb84(33)=-spbk5l4*abb84(32)
      abb84(34)=abb84(31)*spbk5l4
      abb84(35)=abb84(22)*spak1l3
      abb84(36)=abb84(15)+abb84(21)
      abb84(36)=abb84(36)*abb84(35)
      abb84(37)=spbk5e2*abb84(35)
      abb84(38)=spbk5l4*abb84(35)
      abb84(39)=abb84(22)*spbk5l3
      abb84(40)=spak1e2*abb84(20)
      abb84(40)=abb84(19)+abb84(40)
      abb84(40)=abb84(40)*abb84(39)
      abb84(41)=spbk5l3*abb84(31)
      abb84(23)=spbk5l3*abb84(23)
      abb84(15)=abb84(16)*abb84(27)*abb84(15)
      abb84(28)=abb84(28)-spak1k2
      abb84(28)=abb84(28)*spbk5e2
      abb84(42)=-abb84(28)*abb84(20)
      abb84(15)=abb84(15)+abb84(42)
      abb84(15)=abb84(15)*abb84(22)
      abb84(28)=-abb84(28)*abb84(22)
      abb84(42)=abb84(22)*abb84(16)
      abb84(43)=abb84(42)*abb84(27)
      abb84(44)=spbk5l4*abb84(43)
      abb84(45)=abb84(13)*abb84(27)
      abb84(45)=-spak1k2+abb84(45)
      abb84(39)=abb84(45)*abb84(39)
      abb84(18)=abb84(18)*spbk5k2
      abb84(45)=abb84(18)*abb84(16)
      abb84(19)=abb84(45)*abb84(19)
      abb84(26)=abb84(26)*abb84(18)
      abb84(46)=abb84(26)-spbk5k2
      abb84(46)=abb84(46)*spak1e2
      abb84(20)=-abb84(46)*abb84(20)
      abb84(19)=abb84(19)+abb84(20)
      abb84(19)=abb84(19)*abb84(22)
      abb84(20)=-abb84(46)*abb84(22)
      abb84(42)=abb84(42)*abb84(18)
      abb84(46)=spak1l4*abb84(42)
      abb84(18)=abb84(13)*abb84(18)
      abb84(18)=-spbk5k2+abb84(18)
      abb84(18)=abb84(18)*abb84(35)
      abb84(13)=abb84(13)*abb84(45)
      abb84(16)=-spbk5k2*abb84(16)
      abb84(13)=abb84(16)+abb84(13)
      abb84(7)=4.0_ki*TR*abb84(13)*abb84(7)*abb84(27)
      abb84(13)=-spae2k2*abb84(26)
      abb84(12)=abb84(13)+abb84(12)+abb84(25)
      abb84(9)=abb84(12)*abb84(9)
      abb84(8)=16.0_ki*abb84(8)
      abb84(12)=spbk2e2*abb84(32)
      abb84(13)=abb84(31)*spbk2e2
      abb84(16)=-spbk2e2*abb84(35)
      abb84(25)=-spbk2e2*abb84(43)
      abb84(22)=abb84(22)*spae2k2
      abb84(21)=-abb84(21)*abb84(22)
      abb84(26)=-spbk5e2*abb84(22)
      abb84(22)=spbk5l3*abb84(22)
      abb84(27)=spae2k2*abb84(42)
      R2d84=abb84(11)
      rat2 = rat2 + R2d84
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='84' value='", &
          & R2d84, "'/>"
      end if
   end subroutine
end module p10_part1part21_part25part25part1_abbrevd84h2_qp
