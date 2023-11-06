module     p9_part21part21_part25part25part21_abbrevd113h0
   use p9_part21part21_part25part25part21_config, only: ki
   use p9_part21part21_part25part25part21_kinematics, only: epstensor
   use p9_part21part21_part25part25part21_globalsh0
   implicit none
   private
   complex(ki), dimension(89), public :: abb113
   complex(ki), public :: R2d113
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
      abb113(1)=es12**(-1)
      abb113(2)=sqrt(mdlMT**2)
      abb113(3)=spak2l3**(-1)
      abb113(4)=spbl3k2**(-1)
      abb113(5)=spak2l4**(-1)
      abb113(6)=spbl4k2**(-1)
      abb113(7)=c2-c1
      abb113(8)=spbe2k1*spak1e2
      abb113(9)=-abb113(8)*abb113(7)
      abb113(10)=spbe5e1*spae1e5*abb113(9)
      abb113(11)=spbe2e1*spae1e2
      abb113(12)=abb113(7)*abb113(11)
      abb113(13)=spbe5k1*spak1e5
      abb113(14)=abb113(13)*abb113(12)
      abb113(10)=2.0_ki*abb113(10)+abb113(14)
      abb113(14)=mdlGC45*mdlGC7
      abb113(14)=mdlGC6*i_*abb113(14)**2
      abb113(15)=abb113(1)*abb113(14)
      abb113(10)=2.0_ki/3.0_ki*abb113(10)*abb113(15)
      abb113(16)=2.0_ki*spbk2e1
      abb113(17)=spbe5l4*spak2e5
      abb113(18)=abb113(17)*spae1l4
      abb113(19)=abb113(16)*abb113(18)
      abb113(20)=2.0_ki*spae1k2
      abb113(21)=spal4e5*spbe5k2
      abb113(22)=abb113(21)*spbl4e1
      abb113(23)=abb113(20)*abb113(22)
      abb113(19)=abb113(19)+abb113(23)
      abb113(23)=spal3e5*spbl3e1
      abb113(24)=abb113(23)*spae1k2
      abb113(25)=spae5k5*spbk5e1
      abb113(26)=abb113(25)*spae1k2
      abb113(24)=abb113(24)+abb113(26)
      abb113(27)=spbl3e1*spae1e5
      abb113(28)=abb113(27)*spak2l3
      abb113(28)=abb113(28)+abb113(24)
      abb113(28)=abb113(28)*spbe5k2
      abb113(29)=spbe5l3*spae1l3
      abb113(30)=abb113(29)*spbk2e1
      abb113(31)=spbk5e5*spae1k5
      abb113(32)=abb113(31)*spbk2e1
      abb113(30)=abb113(30)+abb113(32)
      abb113(33)=spae1l3*spbe5e1
      abb113(34)=abb113(33)*spbl3k2
      abb113(34)=abb113(34)+abb113(30)
      abb113(34)=abb113(34)*spak2e5
      abb113(28)=-abb113(19)+abb113(28)+abb113(34)
      abb113(34)=abb113(5)*abb113(6)
      abb113(28)=abb113(28)*abb113(34)
      abb113(35)=spak2l4*spae1e5
      abb113(36)=spak2e5*spae1l4
      abb113(35)=abb113(35)-abb113(36)
      abb113(36)=abb113(3)*abb113(4)
      abb113(37)=abb113(36)*spbe5l4
      abb113(35)=abb113(37)*abb113(35)*spbk2e1
      abb113(38)=spbe5e1*spae1k2*spbl4k2
      abb113(39)=spae1k2*spbl4e1
      abb113(40)=abb113(39)*spbe5k2
      abb113(38)=abb113(38)-abb113(40)
      abb113(40)=abb113(36)*spal4e5
      abb113(38)=abb113(38)*abb113(40)
      abb113(26)=abb113(26)*spbe5k2
      abb113(41)=abb113(32)*spak2e5
      abb113(26)=abb113(26)+abb113(41)
      abb113(26)=abb113(26)*abb113(36)
      abb113(28)=-abb113(35)+abb113(28)-abb113(38)-abb113(26)
      abb113(35)=mdlMh**2
      abb113(38)=abb113(8)*abb113(35)
      abb113(41)=abb113(2)**2
      abb113(42)=abb113(38)*abb113(41)
      abb113(28)=abb113(28)*abb113(42)
      abb113(43)=abb113(23)+abb113(25)
      abb113(44)=spae1l4*abb113(43)
      abb113(45)=2.0_ki*spae1l4
      abb113(46)=spal4e5*spbl4e1
      abb113(47)=abb113(45)*abb113(46)
      abb113(44)=abb113(44)-abb113(47)
      abb113(27)=abb113(27)*spal3l4
      abb113(27)=abb113(27)-abb113(44)
      abb113(27)=abb113(27)*spbe5l4
      abb113(47)=abb113(29)*abb113(25)
      abb113(48)=abb113(23)*abb113(31)
      abb113(47)=abb113(47)+abb113(48)
      abb113(48)=abb113(29)+abb113(31)
      abb113(49)=spbl4e1*abb113(48)
      abb113(33)=abb113(33)*spbl4l3
      abb113(33)=abb113(33)-abb113(49)
      abb113(33)=abb113(33)*spal4e5
      abb113(27)=abb113(47)+abb113(27)+abb113(33)
      abb113(27)=abb113(27)*abb113(8)
      abb113(33)=abb113(41)*abb113(8)
      abb113(50)=abb113(33)*spbe5e1
      abb113(51)=abb113(50)*spae1e5
      abb113(27)=abb113(27)+2.0_ki*abb113(51)
      abb113(27)=abb113(27)*abb113(41)
      abb113(27)=abb113(28)-abb113(27)
      abb113(27)=abb113(27)*abb113(7)
      abb113(28)=spak1e5*spak2l4
      abb113(52)=spak1l4*spak2e5
      abb113(28)=abb113(28)-abb113(52)
      abb113(28)=abb113(37)*spbk2k1*abb113(28)
      abb113(37)=spbe5k1*spbl4k2
      abb113(52)=spbl4k1*spbe5k2
      abb113(37)=abb113(37)-abb113(52)
      abb113(37)=abb113(40)*spak1k2*abb113(37)
      abb113(52)=spbk5e5*spak1k5
      abb113(53)=abb113(52)*spbk2k1
      abb113(54)=abb113(53)*spak2e5
      abb113(55)=spae5k5*spbk5k1
      abb113(56)=abb113(55)*spak1k2
      abb113(57)=abb113(56)*spbe5k2
      abb113(54)=abb113(54)+abb113(57)
      abb113(54)=abb113(54)*abb113(36)
      abb113(28)=abb113(54)+abb113(28)+abb113(37)
      abb113(37)=spal3e5*spbl3k1
      abb113(57)=abb113(37)*spak1k2
      abb113(56)=abb113(57)+abb113(56)
      abb113(57)=spbl3k1*spak1e5
      abb113(58)=abb113(57)*spak2l3
      abb113(58)=abb113(58)+abb113(56)
      abb113(58)=abb113(58)*spbe5k2
      abb113(59)=spbe5l3*spak1l3
      abb113(60)=abb113(59)*spbk2k1
      abb113(60)=abb113(60)+abb113(53)
      abb113(61)=spak1l3*spbe5k1
      abb113(62)=abb113(61)*spbl3k2
      abb113(62)=abb113(62)+abb113(60)
      abb113(62)=abb113(62)*spak2e5
      abb113(58)=abb113(58)+abb113(62)
      abb113(62)=abb113(21)*spbl4k1
      abb113(63)=abb113(62)*spak1k2
      abb113(64)=abb113(17)*spak1l4
      abb113(65)=abb113(64)*spbk2k1
      abb113(63)=abb113(63)+abb113(65)
      abb113(58)=-abb113(63)+1.0_ki/2.0_ki*abb113(58)
      abb113(58)=abb113(58)*abb113(34)
      abb113(28)=-abb113(58)+1.0_ki/2.0_ki*abb113(28)
      abb113(58)=abb113(41)*abb113(35)
      abb113(28)=abb113(28)*abb113(58)
      abb113(65)=spbl4k2*spak2k5
      abb113(66)=spbl4k1*spak1k5
      abb113(65)=abb113(65)-abb113(66)
      abb113(65)=abb113(65)*spbk5e5
      abb113(66)=spak2l3*spbl4k2
      abb113(67)=spak1l3*spbl4k1
      abb113(66)=abb113(66)-abb113(67)
      abb113(67)=abb113(66)*spbe5l3
      abb113(65)=abb113(65)+abb113(67)
      abb113(67)=spbe5k2*spak2l3
      abb113(61)=abb113(67)-abb113(61)
      abb113(67)=abb113(61)*spbl4l3
      abb113(67)=abb113(67)-abb113(65)
      abb113(67)=abb113(67)*spal4e5
      abb113(68)=spbl3k2*spak2k5
      abb113(69)=spbl3k1*spak1k5
      abb113(68)=abb113(68)-abb113(69)
      abb113(68)=abb113(68)*spbk5e5
      abb113(69)=abb113(68)*spal3e5
      abb113(70)=spak2l3*spbk5k2
      abb113(71)=spak1l3*spbk5k1
      abb113(70)=abb113(70)-abb113(71)
      abb113(70)=abb113(70)*spae5k5
      abb113(71)=abb113(70)*spbe5l3
      abb113(69)=abb113(69)+abb113(71)
      abb113(67)=abb113(67)+abb113(69)
      abb113(71)=spak2l4*spbk5k2
      abb113(72)=spak1l4*spbk5k1
      abb113(71)=abb113(71)-abb113(72)
      abb113(71)=abb113(71)*spae5k5
      abb113(72)=spbl3k2*spak2l4
      abb113(73)=spbl3k1*spak1l4
      abb113(72)=abb113(72)-abb113(73)
      abb113(73)=abb113(72)*spal3e5
      abb113(71)=abb113(71)+abb113(73)
      abb113(73)=spak2e5*spbl3k2
      abb113(57)=abb113(73)-abb113(57)
      abb113(73)=abb113(57)*spal3l4
      abb113(73)=abb113(73)-abb113(71)
      abb113(74)=spbl4k2*spak2l4
      abb113(75)=spbl4k1*spak1l4
      abb113(74)=abb113(74)-abb113(75)
      abb113(74)=abb113(74)*spal4e5
      abb113(73)=abb113(74)+1.0_ki/2.0_ki*abb113(73)
      abb113(73)=abb113(73)*spbe5l4
      abb113(13)=abb113(13)*abb113(41)
      abb113(67)=1.0_ki/2.0_ki*abb113(67)+abb113(73)-abb113(13)
      abb113(67)=abb113(67)*abb113(41)
      abb113(28)=abb113(28)-abb113(67)
      abb113(28)=abb113(28)*abb113(12)
      abb113(27)=abb113(28)+abb113(27)
      abb113(27)=abb113(27)*abb113(15)
      abb113(28)=abb113(36)*spae1k2
      abb113(67)=abb113(22)*abb113(28)
      abb113(73)=abb113(36)*spbk2e1
      abb113(75)=abb113(18)*abb113(73)
      abb113(26)=-abb113(26)+abb113(67)+abb113(75)
      abb113(67)=-abb113(26)*abb113(34)*abb113(9)
      abb113(75)=abb113(36)*spak1k2
      abb113(76)=abb113(62)*abb113(75)
      abb113(77)=abb113(36)*spbk2k1
      abb113(78)=abb113(64)*abb113(77)
      abb113(54)=-abb113(54)+abb113(76)+abb113(78)
      abb113(76)=-abb113(34)*abb113(7)
      abb113(78)=abb113(11)*abb113(54)*abb113(76)
      abb113(67)=2.0_ki*abb113(67)+abb113(78)
      abb113(78)=mdlMh**4
      abb113(67)=abb113(15)*abb113(78)*abb113(67)
      abb113(30)=abb113(30)*spak2e5
      abb113(24)=abb113(24)*spbe5k2
      abb113(19)=-abb113(19)+abb113(30)+abb113(24)
      abb113(19)=abb113(19)*abb113(34)
      abb113(19)=abb113(19)+abb113(26)
      abb113(19)=abb113(19)*abb113(38)
      abb113(24)=abb113(44)*spbe5l4
      abb113(26)=abb113(49)*spal4e5
      abb113(24)=-abb113(47)+abb113(26)+abb113(24)
      abb113(24)=abb113(24)*abb113(8)
      abb113(19)=abb113(19)+abb113(24)-6.0_ki*abb113(51)
      abb113(19)=abb113(19)*abb113(7)
      abb113(24)=abb113(60)*spak2e5
      abb113(26)=abb113(56)*spbe5k2
      abb113(24)=abb113(24)+abb113(26)
      abb113(24)=-abb113(63)+1.0_ki/2.0_ki*abb113(24)
      abb113(24)=abb113(24)*abb113(34)
      abb113(24)=abb113(24)+1.0_ki/2.0_ki*abb113(54)
      abb113(24)=abb113(24)*abb113(35)
      abb113(26)=abb113(65)*spal4e5
      abb113(26)=abb113(26)-abb113(69)
      abb113(30)=-abb113(74)+1.0_ki/2.0_ki*abb113(71)
      abb113(30)=abb113(30)*spbe5l4
      abb113(13)=abb113(24)-abb113(30)-3.0_ki*abb113(13)-1.0_ki/2.0_ki*abb113(2&
      &6)
      abb113(13)=-abb113(13)*abb113(12)
      abb113(13)=abb113(13)+abb113(19)
      abb113(13)=abb113(13)*abb113(15)
      abb113(19)=spbe5l4*spae1l4
      abb113(24)=abb113(29)-4.0_ki*abb113(19)+3.0_ki*abb113(31)
      abb113(24)=abb113(24)*abb113(33)
      abb113(26)=abb113(34)*spbe5k2
      abb113(29)=abb113(26)*spae1k2
      abb113(30)=abb113(28)*spbe5k2
      abb113(29)=-abb113(30)+4.0_ki*abb113(29)
      abb113(29)=abb113(29)*abb113(42)
      abb113(24)=abb113(24)-abb113(29)
      abb113(29)=abb113(7)*abb113(15)
      abb113(24)=-abb113(24)*abb113(29)
      abb113(44)=2.0_ki*abb113(15)
      abb113(47)=abb113(9)*abb113(44)
      abb113(49)=abb113(78)*abb113(26)
      abb113(51)=abb113(49)*abb113(28)*abb113(47)
      abb113(20)=abb113(20)*abb113(26)
      abb113(20)=abb113(20)-abb113(30)
      abb113(20)=abb113(20)*abb113(38)
      abb113(30)=abb113(45)*spbe5l4
      abb113(30)=abb113(30)-abb113(48)
      abb113(30)=abb113(30)*abb113(8)
      abb113(20)=abb113(30)+abb113(20)
      abb113(20)=abb113(20)*abb113(29)
      abb113(23)=abb113(23)-4.0_ki*abb113(46)+3.0_ki*abb113(25)
      abb113(23)=abb113(23)*abb113(33)
      abb113(30)=abb113(34)*spak2e5
      abb113(45)=abb113(30)*spbk2e1
      abb113(48)=abb113(73)*spak2e5
      abb113(45)=-abb113(48)+4.0_ki*abb113(45)
      abb113(45)=abb113(45)*abb113(42)
      abb113(23)=abb113(23)-abb113(45)
      abb113(23)=-abb113(23)*abb113(29)
      abb113(45)=abb113(78)*abb113(30)
      abb113(54)=abb113(45)*abb113(73)*abb113(47)
      abb113(16)=abb113(16)*abb113(30)
      abb113(16)=abb113(16)-abb113(48)
      abb113(16)=abb113(16)*abb113(38)
      abb113(43)=-abb113(43)+2.0_ki*abb113(46)
      abb113(43)=abb113(43)*abb113(8)
      abb113(16)=abb113(43)+abb113(16)
      abb113(16)=abb113(16)*abb113(29)
      abb113(29)=-6.0_ki*abb113(33)*abb113(29)
      abb113(43)=-spal4e5*abb113(7)
      abb113(48)=abb113(44)*abb113(50)*abb113(43)
      abb113(33)=abb113(33)*spae1e5
      abb113(56)=-spbe5l4*abb113(7)
      abb113(60)=abb113(44)*abb113(33)*abb113(56)
      abb113(63)=-abb113(41)*abb113(7)
      abb113(65)=abb113(63)*abb113(11)
      abb113(61)=-abb113(61)*abb113(65)
      abb113(50)=abb113(50)*spae1l3*abb113(7)
      abb113(50)=2.0_ki*abb113(50)+abb113(61)
      abb113(50)=abb113(50)*abb113(15)
      abb113(57)=-abb113(57)*abb113(65)
      abb113(33)=abb113(33)*spbl3e1*abb113(7)
      abb113(33)=2.0_ki*abb113(33)+abb113(57)
      abb113(33)=abb113(33)*abb113(15)
      abb113(57)=spbe5l3*spak2l3
      abb113(61)=spbk5e5*spak2k5
      abb113(65)=abb113(57)+3.0_ki*abb113(61)
      abb113(69)=spbe5l4*spak2l4
      abb113(65)=-2.0_ki*abb113(69)+1.0_ki/2.0_ki*abb113(65)
      abb113(65)=abb113(65)*abb113(41)
      abb113(71)=spbe5k1*abb113(75)*abb113(58)
      abb113(65)=abb113(65)+abb113(71)
      abb113(65)=-abb113(65)*abb113(12)
      abb113(42)=-abb113(42)*abb113(7)
      abb113(71)=abb113(42)*spbe5e1
      abb113(74)=-abb113(28)*abb113(71)
      abb113(65)=2.0_ki*abb113(74)+abb113(65)
      abb113(65)=abb113(65)*abb113(15)
      abb113(74)=abb113(15)*abb113(11)
      abb113(79)=abb113(7)*abb113(74)
      abb113(57)=abb113(57)+abb113(61)
      abb113(57)=-abb113(69)+1.0_ki/2.0_ki*abb113(57)
      abb113(57)=-abb113(57)*abb113(79)
      abb113(61)=spal3e5*spbl3k2
      abb113(69)=spae5k5*spbk5k2
      abb113(80)=abb113(61)+3.0_ki*abb113(69)
      abb113(81)=2.0_ki*spal4e5
      abb113(82)=abb113(81)*spbl4k2
      abb113(80)=-abb113(82)+1.0_ki/2.0_ki*abb113(80)
      abb113(80)=abb113(80)*abb113(41)
      abb113(82)=spak1e5*abb113(77)*abb113(58)
      abb113(80)=abb113(80)+abb113(82)
      abb113(80)=-abb113(80)*abb113(12)
      abb113(42)=abb113(42)*spae1e5
      abb113(73)=-abb113(73)*abb113(42)
      abb113(73)=2.0_ki*abb113(73)+abb113(80)
      abb113(73)=abb113(73)*abb113(15)
      abb113(61)=abb113(61)+abb113(69)
      abb113(69)=spal4e5*spbl4k2
      abb113(61)=-abb113(69)+1.0_ki/2.0_ki*abb113(61)
      abb113(61)=-abb113(61)*abb113(79)
      abb113(69)=abb113(63)*abb113(74)
      abb113(69)=3.0_ki*abb113(69)
      abb113(71)=abb113(44)*abb113(30)*abb113(71)
      abb113(42)=abb113(44)*abb113(26)*abb113(42)
      abb113(44)=abb113(75)*spbe5k2
      abb113(44)=1.0_ki/2.0_ki*abb113(44)
      abb113(80)=abb113(26)*spak1k2
      abb113(82)=-abb113(44)+2.0_ki*abb113(80)
      abb113(82)=abb113(82)*abb113(58)
      abb113(83)=abb113(59)+3.0_ki*abb113(52)
      abb113(84)=spbe5l4*spak1l4
      abb113(83)=-2.0_ki*abb113(84)+1.0_ki/2.0_ki*abb113(83)
      abb113(83)=abb113(83)*abb113(41)
      abb113(82)=abb113(83)-abb113(82)
      abb113(82)=abb113(82)*abb113(79)
      abb113(49)=abb113(49)*abb113(75)*abb113(79)
      abb113(44)=abb113(44)-abb113(80)
      abb113(44)=abb113(44)*abb113(35)
      abb113(52)=abb113(59)+abb113(52)
      abb113(44)=-abb113(84)+abb113(44)+1.0_ki/2.0_ki*abb113(52)
      abb113(44)=abb113(44)*abb113(79)
      abb113(52)=abb113(77)*spak2e5
      abb113(52)=1.0_ki/2.0_ki*abb113(52)
      abb113(59)=abb113(30)*spbk2k1
      abb113(80)=-abb113(52)+2.0_ki*abb113(59)
      abb113(58)=abb113(80)*abb113(58)
      abb113(80)=abb113(37)+3.0_ki*abb113(55)
      abb113(81)=abb113(81)*spbl4k1
      abb113(80)=-abb113(81)+1.0_ki/2.0_ki*abb113(80)
      abb113(80)=abb113(80)*abb113(41)
      abb113(58)=abb113(80)-abb113(58)
      abb113(58)=abb113(58)*abb113(79)
      abb113(45)=abb113(45)*abb113(77)*abb113(79)
      abb113(52)=abb113(52)-abb113(59)
      abb113(52)=abb113(52)*abb113(35)
      abb113(37)=abb113(37)+abb113(55)
      abb113(59)=spal4e5*spbl4k1
      abb113(37)=-abb113(59)+abb113(52)+1.0_ki/2.0_ki*abb113(37)
      abb113(37)=abb113(37)*abb113(79)
      abb113(46)=abb113(46)*spae1l3
      abb113(52)=abb113(25)*spae1l3
      abb113(46)=abb113(46)-abb113(52)
      abb113(46)=abb113(46)*abb113(8)*spbe5l4
      abb113(59)=abb113(38)*abb113(34)
      abb113(77)=abb113(59)*abb113(17)*spbk2e1*spae1l3
      abb113(46)=abb113(46)+abb113(77)
      abb113(46)=abb113(46)*abb113(7)
      abb113(66)=abb113(66)*spal4e5
      abb113(66)=abb113(66)-abb113(70)
      abb113(66)=abb113(66)*spbe5l4
      abb113(77)=abb113(34)*abb113(35)
      abb113(80)=abb113(77)*abb113(17)*spbk2k1*spak1l3
      abb113(66)=abb113(66)-abb113(80)
      abb113(11)=1.0_ki/2.0_ki*abb113(11)
      abb113(80)=abb113(7)*abb113(11)
      abb113(66)=abb113(66)*abb113(80)
      abb113(46)=abb113(66)+abb113(46)
      abb113(46)=abb113(46)*abb113(15)
      abb113(9)=abb113(9)*abb113(15)
      abb113(66)=spbe5l4*spae1l3*abb113(9)
      abb113(81)=abb113(11)*abb113(15)
      abb113(56)=abb113(56)*abb113(81)
      abb113(83)=spak2l3*abb113(56)
      abb113(56)=-spak1l3*abb113(56)
      abb113(36)=abb113(11)*abb113(76)*abb113(78)*abb113(36)
      abb113(78)=abb113(17)*abb113(36)
      abb113(25)=abb113(28)*abb113(25)
      abb113(39)=abb113(39)*abb113(40)
      abb113(25)=abb113(25)-abb113(39)
      abb113(39)=-abb113(38)*abb113(7)
      abb113(25)=spbe5l4*abb113(25)*abb113(39)
      abb113(55)=abb113(55)*abb113(75)
      abb113(85)=spbl4k1*abb113(40)*spak1k2
      abb113(55)=abb113(85)-abb113(55)
      abb113(85)=1.0_ki/2.0_ki*abb113(35)
      abb113(55)=abb113(55)*abb113(85)*spbe5l4
      abb113(86)=abb113(17)*abb113(41)
      abb113(55)=abb113(55)+abb113(86)
      abb113(55)=-abb113(55)*abb113(12)
      abb113(25)=abb113(55)+abb113(25)
      abb113(25)=abb113(1)*abb113(25)
      abb113(25)=abb113(78)+abb113(25)
      abb113(25)=abb113(25)*abb113(14)
      abb113(55)=abb113(39)*abb113(15)
      abb113(28)=abb113(28)*spbe5l4*abb113(55)
      abb113(35)=-abb113(35)*abb113(7)
      abb113(78)=abb113(35)*abb113(81)
      abb113(75)=-abb113(75)*spbe5l4*abb113(78)
      abb113(86)=abb113(74)*abb113(63)*spak1e5
      abb113(87)=-spbe5l4*abb113(86)
      abb113(88)=abb113(19)-abb113(31)
      abb113(89)=spal4e5*spbl3e1
      abb113(8)=abb113(8)*abb113(88)*abb113(89)
      abb113(59)=abb113(59)*abb113(21)*spae1k2*spbl3e1
      abb113(8)=abb113(8)+abb113(59)
      abb113(8)=abb113(8)*abb113(7)
      abb113(59)=abb113(72)*spbe5l4
      abb113(59)=abb113(59)-abb113(68)
      abb113(59)=abb113(59)*spal4e5
      abb113(72)=abb113(77)*abb113(21)*spak1k2*spbl3k1
      abb113(59)=abb113(59)-abb113(72)
      abb113(59)=abb113(59)*abb113(80)
      abb113(8)=abb113(59)+abb113(8)
      abb113(8)=abb113(8)*abb113(15)
      abb113(9)=abb113(89)*abb113(9)
      abb113(43)=abb113(43)*abb113(81)
      abb113(59)=spbl3k2*abb113(43)
      abb113(43)=-spbl3k1*abb113(43)
      abb113(38)=abb113(76)*abb113(38)
      abb113(31)=abb113(31)*spak2e5
      abb113(18)=abb113(31)-abb113(18)
      abb113(18)=abb113(18)*spbl3e1*abb113(38)
      abb113(11)=abb113(35)*abb113(11)
      abb113(31)=abb113(68)*spak2e5
      abb113(35)=abb113(64)*spbl3k1
      abb113(31)=abb113(35)+abb113(31)
      abb113(31)=abb113(31)*abb113(34)
      abb113(17)=spbl3k2*abb113(6)*abb113(17)
      abb113(17)=abb113(31)-abb113(17)
      abb113(17)=abb113(17)*abb113(11)
      abb113(17)=abb113(17)+abb113(18)
      abb113(17)=abb113(17)*abb113(15)
      abb113(18)=abb113(30)*spbl3e1*abb113(55)
      abb113(7)=abb113(7)*abb113(81)
      abb113(30)=abb113(77)*spak2e5
      abb113(31)=abb113(7)*abb113(30)
      abb113(34)=-spbl3k2*abb113(31)
      abb113(31)=spbl3k1*abb113(31)
      abb113(35)=abb113(21)*abb113(36)
      abb113(32)=abb113(32)*abb113(40)
      abb113(36)=abb113(40)*spbk2e1
      abb113(19)=abb113(36)*abb113(19)
      abb113(19)=abb113(32)-abb113(19)
      abb113(19)=abb113(19)*abb113(39)
      abb113(32)=abb113(53)*abb113(40)
      abb113(39)=abb113(40)*spbk2k1
      abb113(40)=abb113(39)*abb113(84)
      abb113(32)=abb113(32)-abb113(40)
      abb113(32)=abb113(32)*abb113(85)
      abb113(40)=abb113(21)*abb113(41)
      abb113(32)=abb113(32)-abb113(40)
      abb113(12)=abb113(32)*abb113(12)
      abb113(12)=abb113(12)+abb113(19)
      abb113(12)=abb113(1)*abb113(12)
      abb113(12)=abb113(35)+abb113(12)
      abb113(12)=abb113(12)*abb113(14)
      abb113(14)=abb113(36)*abb113(55)
      abb113(19)=-abb113(39)*abb113(78)
      abb113(22)=abb113(22)*spae1l3
      abb113(32)=abb113(52)*spbe5k2
      abb113(22)=abb113(22)-abb113(32)
      abb113(22)=-abb113(22)*abb113(38)
      abb113(32)=abb113(70)*spbe5k2
      abb113(35)=abb113(62)*spak1l3
      abb113(32)=abb113(35)+abb113(32)
      abb113(32)=abb113(32)*abb113(6)
      abb113(21)=abb113(21)*spak2l3
      abb113(21)=abb113(32)-abb113(21)
      abb113(11)=abb113(5)*abb113(21)*abb113(11)
      abb113(11)=abb113(11)+abb113(22)
      abb113(11)=abb113(11)*abb113(15)
      abb113(15)=abb113(26)*spae1l3*abb113(55)
      abb113(21)=abb113(77)*spbe5k2
      abb113(7)=abb113(7)*abb113(21)
      abb113(22)=-spak2l3*abb113(7)
      abb113(7)=spak1l3*abb113(7)
      abb113(21)=-abb113(21)*abb113(86)
      abb113(26)=abb113(74)*abb113(63)*spbe5k1
      abb113(32)=-spal4e5*abb113(26)
      abb113(26)=-abb113(30)*abb113(26)
      R2d113=abb113(10)
      rat2 = rat2 + R2d113
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='113' value='", &
          & R2d113, "'/>"
      end if
   end subroutine
end module p9_part21part21_part25part25part21_abbrevd113h0
