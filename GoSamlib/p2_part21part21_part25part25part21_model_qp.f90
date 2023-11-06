module     p2_part21part21_part25part25part21_model_qp
   ! Model parameters for the model: FeynRules,EWChL_UFO
   use p2_part21part21_part25part25part21_config, only: ki => ki_qp, &
   & renormalisation, reduction_interoperation, &
   & reduction_interoperation_rescue, deltaOS, &
   & nlo_prefactors, convert_to_cdr
   implicit none
   private :: ki
   private :: renormalisation, reduction_interoperation
   private :: reduction_interoperation_rescue
   private :: deltaOS, nlo_prefactors
   real(ki), parameter :: sqrt2 = &
      &1.414213562373095048801688724209698078&
      &5696718753769480731766797379_ki
   real(ki), parameter :: sqrt3 = &
      &1.732050807568877293527446341505872366&
      &9428052538103806280558069795_ki
   complex(ki) :: gauge1z = (       0.000000000000000_ki,        0.000000000000&
   &000_ki)
   complex(ki) :: gauge2z = (       0.000000000000000_ki,        0.000000000000&
   &000_ki)
   complex(ki) :: gauge5z = (       0.000000000000000_ki,        0.000000000000&
   &000_ki)
   real(ki) :: mdlaEWM1 =      127.900000000000006_ki
   real(ki),parameter :: mdlaS =        1.000000000000000_ki
   real(ki) :: mdlcb =        1.000000000000000_ki
   real(ki) :: mdlcgaga =        0.000000000000000_ki
   real(ki) :: mdlcgg =        0.000000000000000_ki
   real(ki) :: mdlcgghh =        0.000000000000000_ki
   real(ki) :: mdlchhh =        1.000000000000000_ki
   real(ki) :: mdlchhhh =        1.000000000000000_ki
   real(ki) :: mdlct =        1.000000000000000_ki
   real(ki) :: mdlctau =        1.000000000000000_ki
   real(ki) :: mdlcthh =        0.000000000000000_ki
   real(ki) :: mdlcV =        1.000000000000000_ki
   real(ki) :: mdlcVhh =        1.000000000000000_ki
   real(ki) :: mdlcZga =        0.000000000000000_ki
   real(ki),parameter :: mdlfloat1 =        2.000000000000000_ki
   real(ki),parameter :: mdlfloat2 =        4.000000000000000_ki
   real(ki),parameter :: mdlfloat3 =        3.000000000000000_ki
   real(ki),parameter :: mdlfloat4 =      125.000000000000000_ki
   real(ki),parameter :: mdlfloat5 =        8.000000000000000_ki
   real(ki),parameter :: mdlfloat6 =        0.000000000000000_ki
   real(ki),parameter :: mdlfloat7 =        6.000000000000000_ki
   real(ki) :: mdlGf =        0.000011663700000_ki
   real(ki),parameter :: mdlMB =        0.000000000000000_ki
   real(ki),parameter :: mdlMC =        0.000000000000000_ki
   real(ki),parameter :: mdlMD =        0.000000000000000_ki
   real(ki),parameter :: mdlMe =        0.000000000000000_ki
   real(ki) :: mdlMMU =        0.105660000000000_ki
   real(ki),parameter :: mdlMS =        0.000000000000000_ki
   real(ki) :: mdlMT =      173.000000000000000_ki
   real(ki) :: mdlMTA =        1.777000000000000_ki
   real(ki),parameter :: mdlMU =        0.000000000000000_ki
   real(ki) :: mdlMZ =       91.187600000000003_ki
   real(ki),parameter :: mdlWh =        0.000000000000000_ki
   real(ki),parameter :: mdlWT =        0.000000000000000_ki
   real(ki),parameter :: mdlWW =        0.000000000000000_ki
   real(ki),parameter :: mdlWZ =        0.000000000000000_ki
   real(ki),parameter :: mdlymb =        0.000000000000000_ki
   real(ki) :: mdlymc =        0.000000000000000_ki
   real(ki) :: mdlymdo =        0.000000000000000_ki
   real(ki) :: mdlyme =        0.000000000000000_ki
   real(ki) :: mdlymm =        0.000000000000000_ki
   real(ki) :: mdlyms =        0.000000000000000_ki
   real(ki) :: mdlymt =      173.000000000000000_ki
   real(ki) :: mdlymtau =        1.777000000000000_ki
   real(ki) :: mdlymup =        0.000000000000000_ki
   real(ki) :: NC =        3.000000000000000_ki
   real(ki) :: Nf =        5.000000000000000_ki
   real(ki) :: Nfgen =       -1.000000000000000_ki
   real(ki) :: mdlloop
   complex(ki) :: mdlGC53
   complex(ki) :: mdlGC52
   complex(ki) :: mdlGC51
   complex(ki) :: mdlGC50
   complex(ki) :: mdlGC56
   complex(ki) :: mdlGC55
   complex(ki) :: mdlGC54
   real(ki) :: mdlyyt
   real(ki) :: mdlyyb
   complex(ki) :: mdlGC48
   complex(ki) :: mdlGC49
   complex(ki) :: mdlGC4
   real(ki) :: mdlMW
   complex(ki) :: mdlGC40
   complex(ki) :: mdlGC41
   complex(ki) :: mdlGC42
   complex(ki) :: mdlGC43
   complex(ki) :: mdlGC44
   complex(ki) :: mdlGC45
   complex(ki) :: mdlGC46
   complex(ki) :: mdlGC47
   complex(ki) :: mdlGC7
   complex(ki) :: mdlGC6
   complex(ki) :: mdlGC5
   real(ki) :: mdlydo
   complex(ki) :: mdlGC3
   complex(ki) :: mdlGC2
   complex(ki) :: mdlGC1
   real(ki) :: mdlMh
   complex(ki) :: mdlGC9
   complex(ki) :: mdlGC8
   real(ki) :: mdlG
   real(ki) :: mdlyytau
   real(ki) :: mdlytau
   complex(ki) :: mdlCKM2x3
   complex(ki) :: mdlCKM2x2
   complex(ki) :: mdlCKM2x1
   real(ki) :: mdlee
   complex(ki) :: mdlCKM1x1
   complex(ki) :: mdlCKM1x2
   complex(ki) :: mdlCKM1x3
   real(ki) :: mdlv
   real(ki) :: mdlsw2
   complex(ki) :: mdlGC17
   complex(ki) :: mdlGC16
   complex(ki) :: mdlGC15
   complex(ki) :: mdlGC14
   complex(ki) :: mdlGC13
   complex(ki) :: mdlGC12
   complex(ki) :: mdlGC11
   complex(ki) :: mdlGC10
   real(ki) :: Nfrat
   real(ki) :: mdlvev
   complex(ki) :: mdlGC19
   complex(ki) :: mdlGC18
   real(ki) :: mdlyup
   real(ki) :: mdlsw
   complex(ki) :: mdlCKM3x2
   complex(ki) :: mdlCKM3x3
   complex(ki) :: mdlCKM3x1
   complex(ki) :: mdlGC33
   real(ki) :: mdlcw
   real(ki) :: mdlaEW
   complex(ki) :: mdlGC34
   real(ki) :: mdlgw
   complex(ki) :: mdlGC31
   complex(ki) :: mdlGC30
   real(ki) :: mdlnormh4
   complex(ki) :: mdlGC32
   complex(ki) :: mdlGC35
   real(ki) :: mdlnormh3
   complex(ki) :: mdlGC37
   complex(ki) :: mdlGC36
   complex(ki) :: mdlGC39
   complex(ki) :: mdlGC38
   real(ki) :: mdlZERO
   real(ki) :: mdlye
   real(ki) :: mdlyc
   real(ki) :: mdlyb
   real(ki) :: mdlym
   real(ki) :: mdlg1
   real(ki) :: mdlyt
   real(ki) :: mdlys
   complex(ki) :: mdlGC28
   complex(ki) :: mdlGC29
   complex(ki) :: mdlGC26
   complex(ki) :: mdlGC27
   complex(ki) :: mdlGC24
   complex(ki) :: mdlGC25
   complex(ki) :: mdlGC22
   complex(ki) :: mdlGC23
   complex(ki) :: mdlGC20
   complex(ki) :: mdlGC21
   integer, parameter, private :: line_length = 80
   integer, parameter, private :: name_length = max(9,24)
   character(len=name_length), dimension(32) :: names = (/&
      & "gauge1z  ", &
      & "gauge2z  ", &
      & "gauge5z  ", &
      & "mdlaEWM1 ", &
      & "mdlcb    ", &
      & "mdlcgaga ", &
      & "mdlcgg   ", &
      & "mdlcgghh ", &
      & "mdlchhh  ", &
      & "mdlchhhh ", &
      & "mdlct    ", &
      & "mdlctau  ", &
      & "mdlcthh  ", &
      & "mdlcV    ", &
      & "mdlcVhh  ", &
      & "mdlcZga  ", &
      & "mdlGf    ", &
      & "mdlMMU   ", &
      & "mdlMT    ", &
      & "mdlMTA   ", &
      & "mdlMZ    ", &
      & "mdlymc   ", &
      & "mdlymdo  ", &
      & "mdlyme   ", &
      & "mdlymm   ", &
      & "mdlyms   ", &
      & "mdlymt   ", &
      & "mdlymtau ", &
      & "mdlymup  ", &
      & "NC       ", &
      & "Nf       ", &
      & "Nfgen    "/)
   character(len=1), dimension(3) :: cc = (/ '#', '!', ';'/)
   private :: digit, parsereal, names, cc
contains
!---#[ print_parameter:
   ! Print current parameters / setup to stdout or output_unit
   subroutine   print_parameter(verbose,output_unit)
      implicit none
      logical, intent(in), optional :: verbose
      integer, intent(in), optional :: output_unit
      logical :: is_verbose
      integer :: unit
      real(ki), parameter :: pi = 3.14159265358979323846264&
     &3383279502884197169399375105820974944592307816406286209_ki
      is_verbose = .false.
      if(present(verbose)) then
          is_verbose = verbose
      end if
      unit = 6 ! stdout
      if(present(output_unit)) then
          unit = output_unit
      end if
   write(unit,'(A1,1x,A26)') "#", "--------- SETUP ---------"
   write(unit,'(A1,1x,A18,I2)') "#", "renormalisation = ", renormalisation
   if(convert_to_cdr) then
      write(unit,'(A1,1x,A9,A3)') "#", "scheme = ", "CDR"
   else
      write(unit,'(A1,1x,A9,A4)') "#", "scheme = ", "DRED"
   end if
   if(reduction_interoperation.eq.0) then
      write(unit,'(A1,1x,A15,A7)') "#", "reduction with ", "SAMURAI"
   else if(reduction_interoperation.eq.1) then
      write(unit,'(A1,1x,A15,A7)') "#", "reduction with ", "GOLEM95"
   else if(reduction_interoperation.eq.2) then
      write(unit,'(A1,1x,A15,A15)') "#", "reduction with ", "NINJA"
   else if(reduction_interoperation.eq.3) then
      write(unit,'(A1,1x,A15,A5)') "#", "reduction with ", "PJFRY"
   end if
   if(reduction_interoperation_rescue.ne.reduction_interoperation) then
      if(reduction_interoperation_rescue.eq.0) then
         write(unit,'(A1,1x,A15,A7)') "#", "    --> rescue ", "SAMURAI"
      else if(reduction_interoperation_rescue.eq.1) then
         write(unit,'(A1,1x,A15,A7)') "#", "    --> rescue ", "GOLEM95"
      else if(reduction_interoperation_rescue.eq.2) then
         write(unit,'(A1,1x,A15,A15)') "#", "    --> rescue ", "NINJA"
      else if(reduction_interoperation_rescue.eq.3) then
         write(unit,'(A1,1x,A15,A5)') "#", "    --> rescue ", "PJFRY"
      end if
   end if
   write(unit,'(A1,1x,A21)') "#", "--- ALL PARAMETERS ---"
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")")') "#", "gauge1z= ", gauge1z
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")")') "#", "gauge2z= ", gauge2z
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")")') "#", "gauge5z= ", gauge5z
   write(unit,'(A1,1x,A11,G23.16)') "#", "mdlaEWM1= ", mdlaEWM1
   write(unit,'(A1,1x,A11,G23.16,"const.")') "#", "mdlaS= ", mdlaS
   write(unit,'(A1,1x,A11,G23.16)') "#", "mdlcb= ", mdlcb
   write(unit,'(A1,1x,A11,G23.16)') "#", "mdlcgaga= ", mdlcgaga
   write(unit,'(A1,1x,A11,G23.16)') "#", "mdlcgg= ", mdlcgg
   write(unit,'(A1,1x,A11,G23.16)') "#", "mdlcgghh= ", mdlcgghh
   write(unit,'(A1,1x,A11,G23.16)') "#", "mdlchhh= ", mdlchhh
   write(unit,'(A1,1x,A11,G23.16)') "#", "mdlchhhh= ", mdlchhhh
   write(unit,'(A1,1x,A11,G23.16)') "#", "mdlct= ", mdlct
   write(unit,'(A1,1x,A11,G23.16)') "#", "mdlctau= ", mdlctau
   write(unit,'(A1,1x,A11,G23.16)') "#", "mdlcthh= ", mdlcthh
   write(unit,'(A1,1x,A11,G23.16)') "#", "mdlcV= ", mdlcV
   write(unit,'(A1,1x,A11,G23.16)') "#", "mdlcVhh= ", mdlcVhh
   write(unit,'(A1,1x,A11,G23.16)') "#", "mdlcZga= ", mdlcZga
   write(unit,'(A1,1x,A11,G23.16,"const.")') "#", "mdlfloat1= ", mdlfloat1
   write(unit,'(A1,1x,A11,G23.16,"const.")') "#", "mdlfloat2= ", mdlfloat2
   write(unit,'(A1,1x,A11,G23.16,"const.")') "#", "mdlfloat3= ", mdlfloat3
   write(unit,'(A1,1x,A11,G23.16,"const.")') "#", "mdlfloat4= ", mdlfloat4
   write(unit,'(A1,1x,A11,G23.16,"const.")') "#", "mdlfloat5= ", mdlfloat5
   write(unit,'(A1,1x,A11,G23.16,"const.")') "#", "mdlfloat6= ", mdlfloat6
   write(unit,'(A1,1x,A11,G23.16,"const.")') "#", "mdlfloat7= ", mdlfloat7
   write(unit,'(A1,1x,A11,G23.16)') "#", "mdlGf= ", mdlGf
   write(unit,'(A1,1x,A11,G23.16,"const.")') "#", "mdlMB= ", mdlMB
   write(unit,'(A1,1x,A11,G23.16,"const.")') "#", "mdlMC= ", mdlMC
   write(unit,'(A1,1x,A11,G23.16,"const.")') "#", "mdlMD= ", mdlMD
   write(unit,'(A1,1x,A11,G23.16,"const.")') "#", "mdlMe= ", mdlMe
   write(unit,'(A1,1x,A11,G23.16)') "#", "mdlMMU= ", mdlMMU
   write(unit,'(A1,1x,A11,G23.16,"const.")') "#", "mdlMS= ", mdlMS
   write(unit,'(A1,1x,A11,G23.16)') "#", "mdlMT= ", mdlMT
   write(unit,'(A1,1x,A11,G23.16)') "#", "mdlMTA= ", mdlMTA
   write(unit,'(A1,1x,A11,G23.16,"const.")') "#", "mdlMU= ", mdlMU
   write(unit,'(A1,1x,A11,G23.16)') "#", "mdlMZ= ", mdlMZ
   write(unit,'(A1,1x,A11,G23.16,"const.")') "#", "mdlWh= ", mdlWh
   write(unit,'(A1,1x,A11,G23.16,"const.")') "#", "mdlWT= ", mdlWT
   write(unit,'(A1,1x,A11,G23.16,"const.")') "#", "mdlWW= ", mdlWW
   write(unit,'(A1,1x,A11,G23.16,"const.")') "#", "mdlWZ= ", mdlWZ
   write(unit,'(A1,1x,A11,G23.16,"const.")') "#", "mdlymb= ", mdlymb
   write(unit,'(A1,1x,A11,G23.16)') "#", "mdlymc= ", mdlymc
   write(unit,'(A1,1x,A11,G23.16)') "#", "mdlymdo= ", mdlymdo
   write(unit,'(A1,1x,A11,G23.16)') "#", "mdlyme= ", mdlyme
   write(unit,'(A1,1x,A11,G23.16)') "#", "mdlymm= ", mdlymm
   write(unit,'(A1,1x,A11,G23.16)') "#", "mdlyms= ", mdlyms
   write(unit,'(A1,1x,A11,G23.16)') "#", "mdlymt= ", mdlymt
   write(unit,'(A1,1x,A11,G23.16)') "#", "mdlymtau= ", mdlymtau
   write(unit,'(A1,1x,A11,G23.16)') "#", "mdlymup= ", mdlymup
   write(unit,'(A1,1x,A11,G23.16)') "#", "NC   = ", NC
   write(unit,'(A1,1x,A11,G23.16)') "#", "Nf   = ", Nf
   write(unit,'(A1,1x,A11,G23.16)') "#", "Nfgen= ", Nfgen
   if(is_verbose) then
   write(unit,'(A1,1x,A11,G23.16,"calc.")') "#", "mdlloop= ", mdlloop
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC53= ", m&
   &dlGC53
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC52= ", m&
   &dlGC52
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC51= ", m&
   &dlGC51
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC50= ", m&
   &dlGC50
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC56= ", m&
   &dlGC56
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC55= ", m&
   &dlGC55
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC54= ", m&
   &dlGC54
   write(unit,'(A1,1x,A11,G23.16,"calc.")') "#", "mdlyyt= ", mdlyyt
   write(unit,'(A1,1x,A11,G23.16,"calc.")') "#", "mdlyyb= ", mdlyyb
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC48= ", m&
   &dlGC48
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC49= ", m&
   &dlGC49
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC4= ", md&
   &lGC4
   write(unit,'(A1,1x,A11,G23.16,"calc.")') "#", "mdlMW= ", mdlMW
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC40= ", m&
   &dlGC40
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC41= ", m&
   &dlGC41
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC42= ", m&
   &dlGC42
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC43= ", m&
   &dlGC43
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC44= ", m&
   &dlGC44
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC45= ", m&
   &dlGC45
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC46= ", m&
   &dlGC46
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC47= ", m&
   &dlGC47
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC7= ", md&
   &lGC7
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC6= ", md&
   &lGC6
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC5= ", md&
   &lGC5
   write(unit,'(A1,1x,A11,G23.16,"calc.")') "#", "mdlydo= ", mdlydo
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC3= ", md&
   &lGC3
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC2= ", md&
   &lGC2
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC1= ", md&
   &lGC1
   write(unit,'(A1,1x,A11,G23.16,"calc.")') "#", "mdlMh= ", mdlMh
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC9= ", md&
   &lGC9
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC8= ", md&
   &lGC8
   write(unit,'(A1,1x,A11,G23.16,"calc.")') "#", "mdlG = ", mdlG
   write(unit,'(A1,1x,A11,G23.16,"calc.")') "#", "mdlyytau= ", mdlyytau
   write(unit,'(A1,1x,A11,G23.16,"calc.")') "#", "mdlytau= ", mdlytau
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlCKM2x3= ",&
   & mdlCKM2x3
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlCKM2x2= ",&
   & mdlCKM2x2
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlCKM2x1= ",&
   & mdlCKM2x1
   write(unit,'(A1,1x,A11,G23.16,"calc.")') "#", "mdlee= ", mdlee
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlCKM1x1= ",&
   & mdlCKM1x1
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlCKM1x2= ",&
   & mdlCKM1x2
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlCKM1x3= ",&
   & mdlCKM1x3
   write(unit,'(A1,1x,A11,G23.16,"calc.")') "#", "mdlv = ", mdlv
   write(unit,'(A1,1x,A11,G23.16,"calc.")') "#", "mdlsw2= ", mdlsw2
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC17= ", m&
   &dlGC17
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC16= ", m&
   &dlGC16
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC15= ", m&
   &dlGC15
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC14= ", m&
   &dlGC14
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC13= ", m&
   &dlGC13
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC12= ", m&
   &dlGC12
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC11= ", m&
   &dlGC11
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC10= ", m&
   &dlGC10
   write(unit,'(A1,1x,A11,G23.16,"calc.")') "#", "Nfrat= ", Nfrat
   write(unit,'(A1,1x,A11,G23.16,"calc.")') "#", "mdlvev= ", mdlvev
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC19= ", m&
   &dlGC19
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC18= ", m&
   &dlGC18
   write(unit,'(A1,1x,A11,G23.16,"calc.")') "#", "mdlyup= ", mdlyup
   write(unit,'(A1,1x,A11,G23.16,"calc.")') "#", "mdlsw= ", mdlsw
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlCKM3x2= ",&
   & mdlCKM3x2
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlCKM3x3= ",&
   & mdlCKM3x3
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlCKM3x1= ",&
   & mdlCKM3x1
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC33= ", m&
   &dlGC33
   write(unit,'(A1,1x,A11,G23.16,"calc.")') "#", "mdlcw= ", mdlcw
   write(unit,'(A1,1x,A11,G23.16,"calc.")') "#", "mdlaEW= ", mdlaEW
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC34= ", m&
   &dlGC34
   write(unit,'(A1,1x,A11,G23.16,"calc.")') "#", "mdlgw= ", mdlgw
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC31= ", m&
   &dlGC31
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC30= ", m&
   &dlGC30
   write(unit,'(A1,1x,A11,G23.16,"calc.")') "#", "mdlnormh4= ", mdlnormh4
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC32= ", m&
   &dlGC32
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC35= ", m&
   &dlGC35
   write(unit,'(A1,1x,A11,G23.16,"calc.")') "#", "mdlnormh3= ", mdlnormh3
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC37= ", m&
   &dlGC37
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC36= ", m&
   &dlGC36
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC39= ", m&
   &dlGC39
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC38= ", m&
   &dlGC38
   write(unit,'(A1,1x,A11,G23.16,"calc.")') "#", "mdlZERO= ", mdlZERO
   write(unit,'(A1,1x,A11,G23.16,"calc.")') "#", "mdlye= ", mdlye
   write(unit,'(A1,1x,A11,G23.16,"calc.")') "#", "mdlyc= ", mdlyc
   write(unit,'(A1,1x,A11,G23.16,"calc.")') "#", "mdlyb= ", mdlyb
   write(unit,'(A1,1x,A11,G23.16,"calc.")') "#", "mdlym= ", mdlym
   write(unit,'(A1,1x,A11,G23.16,"calc.")') "#", "mdlg1= ", mdlg1
   write(unit,'(A1,1x,A11,G23.16,"calc.")') "#", "mdlyt= ", mdlyt
   write(unit,'(A1,1x,A11,G23.16,"calc.")') "#", "mdlys= ", mdlys
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC28= ", m&
   &dlGC28
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC29= ", m&
   &dlGC29
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC26= ", m&
   &dlGC26
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC27= ", m&
   &dlGC27
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC24= ", m&
   &dlGC24
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC25= ", m&
   &dlGC25
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC22= ", m&
   &dlGC22
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC23= ", m&
   &dlGC23
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC20= ", m&
   &dlGC20
   write(unit,'(A1,1x,A11,"(",G23.16,G23.16,")"," calc.")') "#", "mdlGC21= ", m&
   &dlGC21
   end if
   write(unit,'(A1,1x,A25)') "#", "-------------------------"
   end subroutine
!---#] print_parameter:
   function     digit(ch, lnr) result(d)
      implicit none
      character(len=1), intent(in) :: ch
      integer, intent(in) :: lnr
      integer :: d
      d = -1
      select case(ch)
         case('0')
            d = 0
         case('1')
            d = 1
         case('2')
            d = 2
         case('3')
            d = 3
         case('4')
            d = 4
         case('5')
            d = 5
         case('6')
            d = 6
         case('7')
            d = 7
         case('8')
            d = 8
         case('9')
            d = 9
         case default
            write(*,'(A21,1x,I5)') 'invalid digit in line', lnr
         end select
   end function digit
   function     parsereal(str, ierr, lnr) result(num)
      implicit none
      character(len=*), intent(in) :: str
      integer, intent(out) :: ierr
      integer, intent(in) :: lnr
      integer, dimension(0:3,0:4), parameter :: DFA = &
      & reshape( (/1,  1,  2, -1,   & ! state = 0
      &            1, -1,  2,  3,   & ! state = 1
      &            2, -1, -1,  3,   & ! state = 2
      &            4,  4, -1, -1,   & ! state = 3
      &            4, -1, -1, -1/), (/4, 5/))
      real(ki) :: num
      integer :: i, expo, ofs, state, cclass, d, s1, s2
      num = 0.0_ki
      expo = 0
      state = 0
      ofs = 0
      s1 = 1
      s2 = 1
      d = -1
      cclass = -1
      do i=1,len(str)
         select case(str(i:i))
         case('_', '''', ' ')
            cycle
         case('+', '-')
            cclass = 1
         case('.')
            cclass = 2
         case('E', 'e')
            cclass = 3
         case default
            cclass = 0
            d = digit(str(i:i), lnr)
            if (d .eq. -1) then
               ierr = 1
               return
            end if
         end select
         if (cclass .eq. 0) then
            select case(state)
            case(0, 1)
               num = 10.0_ki * num + d
            case(2)
               num = 10.0_ki * num + d
               ofs = ofs - 1
            case(4)
               expo = 10 * expo + d
            end select
         elseif ((cclass .eq. 1) .and. (str(i:i) .eq. '-')) then
            if (state .eq. 0) then
               s1 = -1
            else
               s2 = -1
            endif
         end if
         ! Advance in the DFA
         state = DFA(cclass, state)
         if (state < 0) then
            write(*,'(A21,1x,A1,1x,A7,I5)') 'invalid position for', &
            & str(i:i), 'in line', lnr
            ierr = 1
            return
         end if
      end do
      num = s1 * num * 10.0_ki**(ofs + s2 * expo)
      ierr = 0
   end function parsereal
   subroutine     parseline(line,stat,line_number)
      implicit none
      character(len=*), intent(in) :: line
      integer, intent(out) :: stat
      integer, intent(in), optional :: line_number
      character(len=line_length) :: rvalue, ivalue, value
      character(len=name_length) :: name
      real(ki) :: re, im
      integer :: idx, icomma, idx1, idx2, lnr, nidx, ierr, pdg
      if(present(line_number)) then
         lnr = line_number
      else
         lnr = 0
      end if
      idx = scan(line, '=', .false.)
      if (idx .eq. 0) then
         if(present(line_number)) then
            write(*,'(A13,1x,I5)') 'error at line', line_number
         else
            write(*,'(A18)') 'error in parseline'
         end if
         stat = 1
         return
      end if
      name = adjustl(line(1:idx-1))
      value = adjustl(line(idx+1:len(line)))
      idx = scan(value, ',', .false.)
      if (name .eq. "renormalisation") then
         re = parsereal(value, ierr, lnr)
         if (ierr .ne. 0) then
            stat = 1
            return
         end if
         renormalisation = int(re)
      elseif (name .eq. "nlo_prefactors") then
         re = parsereal(value, ierr, lnr)
         if (ierr .ne. 0) then
            stat = 1
            return
         end if
         nlo_prefactors = int(re)
      elseif (name .eq. "deltaOS") then
         re = parsereal(value, ierr, lnr)
         if (ierr .ne. 0) then
            stat = 1
            return
         end if
         deltaOS = int(re)
      elseif (name .eq. "reduction_interoperation") then
         re = parsereal(value, ierr, lnr)
         if (ierr .ne. 0) then
            stat = 1
            return
         end if
         reduction_interoperation = int(re)
      elseif (any(names .eq. name)) then
         do nidx=1,size(names)
            if (names(nidx) .eq. name) exit
         end do
         if (idx .gt. 0) then
            rvalue = value(1:idx-1)
            ivalue = value(idx+1:len(value))
            re = parsereal(rvalue, ierr, lnr)
            if (ierr .ne. 0) then
               stat = 1
               return
            end if
            im = parsereal(ivalue, ierr, lnr)
            if (ierr .ne. 0) then
               stat = 1
               return
            end if
         else
            re = parsereal(value, ierr, lnr)
            if (ierr .ne. 0) then
               stat = 1
               return
            end if
            im = 0.0_ki
         end if
         select case (nidx)
         case(1)
            gauge1z = cmplx(re, im, ki)
         case(2)
            gauge2z = cmplx(re, im, ki)
         case(3)
            gauge5z = cmplx(re, im, ki)
         case(4)
            mdlaEWM1 = re
         case(5)
            mdlcb = re
         case(6)
            mdlcgaga = re
         case(7)
            mdlcgg = re
         case(8)
            mdlcgghh = re
         case(9)
            mdlchhh = re
         case(10)
            mdlchhhh = re
         case(11)
            mdlct = re
         case(12)
            mdlctau = re
         case(13)
            mdlcthh = re
         case(14)
            mdlcV = re
         case(15)
            mdlcVhh = re
         case(16)
            mdlcZga = re
         case(17)
            mdlGf = re
         case(18)
            mdlMMU = re
         case(19)
            mdlMT = re
         case(20)
            mdlMTA = re
         case(21)
            mdlMZ = re
         case(22)
            mdlymc = re
         case(23)
            mdlymdo = re
         case(24)
            mdlyme = re
         case(25)
            mdlymm = re
         case(26)
            mdlyms = re
         case(27)
            mdlymt = re
         case(28)
            mdlymtau = re
         case(29)
            mdlymup = re
         case(30)
            NC = re
         case(31)
            Nf = re
         case(32)
            Nfgen = re
         end select
      elseif (name(1:8).eq."frblock(") then
         idx = scan(name, ')', .false.)
         if (idx.eq.0) then
            if(present(line_number)) then
               write(*,'(A13,1x,I5)') 'error at line', line_number
            else
               write(*,'(A18)') 'error in parseline'
            end if
            stat = 1
            return
         endif
         read(name(9:idx-1),*, iostat=ierr) pdg
         if (ierr.ne.0) then
            write(*,*) "Not an integer:", name(9:idx-1)
            if(present(line_number)) then
               write(*,'(A13,1x,I5)') 'error at line', line_number
            else
               write(*,'(A18)') 'error in parseline'
            end if
            stat = 1
            return
         end if
         select case(pdg)
            case(1)
               mdlcV = parsereal(value, ierr, lnr)
            case(2)
               mdlct = parsereal(value, ierr, lnr)
            case(3)
               mdlcb = parsereal(value, ierr, lnr)
            case(4)
               mdlctau = parsereal(value, ierr, lnr)
            case(5)
               mdlcgaga = parsereal(value, ierr, lnr)
            case(6)
               mdlcgg = parsereal(value, ierr, lnr)
            case(7)
               mdlcZga = parsereal(value, ierr, lnr)
            case(8)
               mdlcgghh = parsereal(value, ierr, lnr)
            case(9)
               mdlcVhh = parsereal(value, ierr, lnr)
            case(10)
               mdlcthh = parsereal(value, ierr, lnr)
            case(11)
               mdlchhh = parsereal(value, ierr, lnr)
            case(12)
               mdlchhhh = parsereal(value, ierr, lnr)
            case default
               write(*,'(A20,1x,I10)') "Cannot set frblock for code:", pdg
               stat = 1
               return
         end select
      elseif (name(1:9).eq."sminputs(") then
         idx = scan(name, ')', .false.)
         if (idx.eq.0) then
            if(present(line_number)) then
               write(*,'(A13,1x,I5)') 'error at line', line_number
            else
               write(*,'(A18)') 'error in parseline'
            end if
            stat = 1
            return
         endif
         read(name(10:idx-1),*, iostat=ierr) pdg
         if (ierr.ne.0) then
            write(*,*) "Not an integer:", name(10:idx-1)
            if(present(line_number)) then
               write(*,'(A13,1x,I5)') 'error at line', line_number
            else
               write(*,'(A18)') 'error in parseline'
            end if
            stat = 1
            return
         end if
         select case(pdg)
            case(1)
               mdlaEWM1 = parsereal(value, ierr, lnr)
            case(2)
               mdlGf = parsereal(value, ierr, lnr)
            case default
               write(*,'(A20,1x,I10)') "Cannot set sminputs for code:", pdg
               stat = 1
               return
         end select
      elseif (name(1:7).eq."yukawa(") then
         idx = scan(name, ')', .false.)
         if (idx.eq.0) then
            if(present(line_number)) then
               write(*,'(A13,1x,I5)') 'error at line', line_number
            else
               write(*,'(A18)') 'error in parseline'
            end if
            stat = 1
            return
         endif
         read(name(8:idx-1),*, iostat=ierr) pdg
         if (ierr.ne.0) then
            write(*,*) "Not an integer:", name(8:idx-1)
            if(present(line_number)) then
               write(*,'(A13,1x,I5)') 'error at line', line_number
            else
               write(*,'(A18)') 'error in parseline'
            end if
            stat = 1
            return
         end if
         select case(pdg)
            case(1)
               mdlymdo = parsereal(value, ierr, lnr)
            case(2)
               mdlymup = parsereal(value, ierr, lnr)
            case(3)
               mdlyms = parsereal(value, ierr, lnr)
            case(4)
               mdlymc = parsereal(value, ierr, lnr)
            case(6)
               mdlymt = parsereal(value, ierr, lnr)
            case(11)
               mdlyme = parsereal(value, ierr, lnr)
            case(13)
               mdlymm = parsereal(value, ierr, lnr)
            case(15)
               mdlymtau = parsereal(value, ierr, lnr)
            case default
               write(*,'(A20,1x,I10)') "Cannot set yukawa for code:", pdg
               stat = 1
               return
         end select
      elseif (name(1:5).eq."mass(") then
         idx = scan(name, ')', .false.)
         if (idx.eq.0) then
            if(present(line_number)) then
               write(*,'(A13,1x,I5)') 'error at line', line_number
            else
               write(*,'(A18)') 'error in parseline'
            end if
            stat = 1
            return
         endif
         read(name(6:idx-1),*, iostat=ierr) pdg
         if (ierr.ne.0) then
            write(*,*) "Not an integer:", name(6:idx-1)
            if(present(line_number)) then
               write(*,'(A13,1x,I5)') 'error at line', line_number
            else
               write(*,'(A18)') 'error in parseline'
            end if
            stat = 1
            return
         end if
         select case(pdg)
            case(23)
               mdlMZ = parsereal(value, ierr, lnr)
            case(13)
               mdlMMU = parsereal(value, ierr, lnr)
            case(6)
               mdlMT = parsereal(value, ierr, lnr)
            case(15)
               mdlMTA = parsereal(value, ierr, lnr)
            case default
               write(*,'(A20,1x,I10)') "Cannot set mass for code:", pdg
               stat = 1
               return
         end select
      elseif (name(1:2).eq."m(" .or. name(1:2).eq."w(") then
         idx = scan(name, ')', .false.)
         if (idx.eq.0) then
            if(present(line_number)) then
               write(*,'(A13,1x,I5)') 'error at line', line_number
            else
               write(*,'(A18)') 'error in parseline'
            end if
            stat = 1
            return
         endif
         read(name(3:idx-1),*, iostat=ierr) pdg
         if (ierr.ne.0) then
            write(*,*) "pdg is not an integer:", name(3:idx-1)
            if(present(line_number)) then
               write(*,'(A13,1x,I5)') 'error at line', line_number
            else
               write(*,'(A18)') 'error in parseline'
            end if
            stat = 1
            return
         end if
         if (name(1:1).eq."m") then
            ! set mass according to PDG code
            select case(pdg)
            case(23)
               mdlMZ = parsereal(value, ierr, lnr)
            case(13)
               mdlMMU = parsereal(value, ierr, lnr)
            case(6)
               mdlMT = parsereal(value, ierr, lnr)
            case(15)
               mdlMTA = parsereal(value, ierr, lnr)
            case default
               write(*,'(A20,1x,I10)') "Cannot set mass for PDG code:", pdg
               stat = 1
               return
            end select
         else
            ! set width according to PDG code
            select case(pdg)
            case default
               write(*,'(A20,1x,I10)') "Cannot set width for PDG code:", pdg
               stat = 1
               return
            end select
         endif
      else
         write(*,'(A20,1x,A20)') 'Unrecognized option:', name
         stat = 1
         return
      end if
      stat = 0
   end subroutine parseline
   subroutine     parse(aunit)
      implicit none
      integer, intent(in) :: aunit
      character(len=line_length) :: line
      integer :: ios, lnr
      lnr = 0
      loop1: do
         read(unit=aunit,fmt='(A80)',iostat=ios) line
         if(ios .ne. 0) exit
         lnr = lnr + 1
         line = adjustl(line)
         if (line .eq. '') cycle loop1
         if (any(cc .eq. line(1:1))) cycle loop1
         call parseline(line,ios,lnr)
         if(ios .ne. 0) then
            write(*,'(A44,I2,A1)') &
            & 'Error while reading parameter file in parse(', aunit, ')'
         end if
      end do loop1
   end subroutine parse
!---#[ SLHA READER:
   subroutine     read_slha(ch, ierr)
      implicit none
      integer, intent(in) :: ch
      integer, intent(out), optional :: ierr
      integer :: lnr, i, l, ofs, ios
      character(len=255) :: line
      integer :: block
      ofs = iachar('A') - iachar('a')
      lnr = 0
      loop1: do
         read(unit=ch,fmt='(A80)',iostat=ios) line
         if(ios .ne. 0) exit
         lnr = lnr + 1
         i = scan(line, '#', .false.)
         if (i .eq. 0) then
            l = len_trim(line)
         else
            l = i - 1
         end if
         if (l .eq. 0) cycle loop1
         ucase: do i = 1, l
            if (line(i:i) >= 'a' .and. line(i:i) <= 'z') then
               line(i:i) = achar(iachar(line(i:i))+ofs)
            end if
         end do ucase
         if (line(1:1) .eq. 'B') then
            if (line(1:5) .eq. 'BLOCK') then
               line = adjustl(line(6:l))
               do i=1,l
                 if (line(i:i) <= ' ') exit
               end do
               l = i
           if ("FRBLOCK" .eq. line(1:l)) then
                  block = 0
           else           if ("SMINPUTS" .eq. line(1:l)) then
                  block = 1
           else           if ("YUKAWA" .eq. line(1:l)) then
                  block = 2
           else           if ("MASS" .eq. line(1:l)) then
                  block = 3
               else
                  block = -1
               end if
            else
               write(*,'(A37,I5)') "Illegal statement in SLHA file, line ", lnr
               if (present(ierr)) ierr = 1
               return
            end if
         else
            ! read a parameter line
            select case(block)
            case(0)
               call read_slha_block_frblock(line(1:l), i)
               if (i .ne. 0) then
                  if (present(ierr)) ierr = 1
                  write(*,'(A44,I5)') &
                  & "Unrecognized line format in SLHA file, line ", lnr
                  return
               end if
            case(1)
               call read_slha_block_sminputs(line(1:l), i)
               if (i .ne. 0) then
                  if (present(ierr)) ierr = 1
                  write(*,'(A44,I5)') &
                  & "Unrecognized line format in SLHA file, line ", lnr
                  return
               end if
            case(2)
               call read_slha_block_yukawa(line(1:l), i)
               if (i .ne. 0) then
                  if (present(ierr)) ierr = 1
                  write(*,'(A44,I5)') &
                  & "Unrecognized line format in SLHA file, line ", lnr
                  return
               end if
            case(3)
               call read_slha_block_mass(line(1:l), i)
               if (i .ne. 0) then
                  if (present(ierr)) ierr = 1
                  write(*,'(A44,I5)') &
                  & "Unrecognized line format in SLHA file, line ", lnr
                  return
               end if
            case default
               cycle loop1
            end select
         end if
      end do loop1
      if (present(ierr)) ierr = 0
   end subroutine read_slha
   subroutine read_slha_block_frblock(line, ierr)
      implicit none
      character(len=*), intent(in) :: line
      integer, intent(out), optional :: ierr
      integer :: idx1,ioerr
      real(ki) :: value
      read(line,*,iostat=ioerr) idx1, value
      if (ioerr .ne. 0) then
         if (present(ierr)) ierr = 1
         return
      end if
      select case(idx1)
      case(1)
         mdlcV = value
      case(2)
         mdlct = value
      case(3)
         mdlcb = value
      case(4)
         mdlctau = value
      case(5)
         mdlcgaga = value
      case(6)
         mdlcgg = value
      case(7)
         mdlcZga = value
      case(8)
         mdlcgghh = value
      case(9)
         mdlcVhh = value
      case(10)
         mdlcthh = value
      case(11)
         mdlchhh = value
      case(12)
         mdlchhhh = value
      end select
      if (present(ierr)) ierr = 0
   end subroutine read_slha_block_frblock
   subroutine read_slha_block_sminputs(line, ierr)
      implicit none
      character(len=*), intent(in) :: line
      integer, intent(out), optional :: ierr
      integer :: idx1,ioerr
      real(ki) :: value
      read(line,*,iostat=ioerr) idx1, value
      if (ioerr .ne. 0) then
         if (present(ierr)) ierr = 1
         return
      end if
      select case(idx1)
      case(1)
         mdlaEWM1 = value
      case(2)
         mdlGf = value
      end select
      if (present(ierr)) ierr = 0
   end subroutine read_slha_block_sminputs
   subroutine read_slha_block_yukawa(line, ierr)
      implicit none
      character(len=*), intent(in) :: line
      integer, intent(out), optional :: ierr
      integer :: idx1,ioerr
      real(ki) :: value
      read(line,*,iostat=ioerr) idx1, value
      if (ioerr .ne. 0) then
         if (present(ierr)) ierr = 1
         return
      end if
      select case(idx1)
      case(1)
         mdlymdo = value
      case(2)
         mdlymup = value
      case(3)
         mdlyms = value
      case(4)
         mdlymc = value
      case(6)
         mdlymt = value
      case(11)
         mdlyme = value
      case(13)
         mdlymm = value
      case(15)
         mdlymtau = value
      end select
      if (present(ierr)) ierr = 0
   end subroutine read_slha_block_yukawa
   subroutine read_slha_block_mass(line, ierr)
      implicit none
      character(len=*), intent(in) :: line
      integer, intent(out), optional :: ierr
      integer :: idx1,ioerr
      real(ki) :: value
      read(line,*,iostat=ioerr) idx1, value
      if (ioerr .ne. 0) then
         if (present(ierr)) ierr = 1
         return
      end if
      select case(idx1)
      case(23)
         mdlMZ = value
      case(13)
         mdlMMU = value
      case(6)
         mdlMT = value
      case(15)
         mdlMTA = value
      end select
      if (present(ierr)) ierr = 0
   end subroutine read_slha_block_mass
!---#] SLHA READER:
!---#[ subroutine set_parameter
   recursive subroutine set_parameter(name, re, im, ierr)
      implicit none
      real(ki), parameter :: pi = 3.14159265358979323846264&
     &3383279502884197169399375105820974944592307816406286209_ki
      character(len=*), intent(in) :: name
      real(ki), intent(in) :: re, im
      integer, intent(out) :: ierr
      integer :: err, pdg, nidx, idx
      complex(ki) :: tmp
      logical :: must_be_real
      must_be_real = .false.
      ierr = 1 ! OK
if (name(1:5).eq."mass(" .and. len_trim(name)>=7) then
         idx = scan(name,")",.false.)
         if (idx.eq.0) then
            idx=len_trim(name)+1
         end if
         read(name(6:idx-1),*, iostat=err) pdg
         if (err.ne.0) then
            ierr=0 !FAIL
            return
         end if
         must_be_real = .true.
         select case(pdg)
            case(23)
               mdlMZ = re
            case(13)
               mdlMMU = re
            case(6)
               mdlMT = re
            case(15)
               mdlMTA = re
            case default
               write(*,'(A20,1x,I10)') "Cannot set mass for PDG code:", pdg
               ierr = 0
               return
            end select
     elseif (len_trim(name)>=8 .and. name(1:6).eq."width(") then
         idx = scan(name,")",.false.)
         if (idx.eq.0) then
            idx=len_trim(name)+1
         end if
         read(name(7:idx-1),*, iostat=err) pdg
         if (err.ne.0) then
            ierr=0 !FAIL
            return
         end if
         must_be_real = .true.
         select case(pdg)
            case default
               write(*,'(A20,1x,I10)') "Cannot set width for PDG code:", pdg
               ierr = 0 !FAIL
               return
            end select
     elseif (name .eq. "FRBLOCK&&1") then
               must_be_real = .true.
               mdlcV = re
     elseif (name .eq. "FRBLOCK&&2") then
               must_be_real = .true.
               mdlct = re
     elseif (name .eq. "FRBLOCK&&3") then
               must_be_real = .true.
               mdlcb = re
     elseif (name .eq. "FRBLOCK&&4") then
               must_be_real = .true.
               mdlctau = re
     elseif (name .eq. "FRBLOCK&&5") then
               must_be_real = .true.
               mdlcgaga = re
     elseif (name .eq. "FRBLOCK&&6") then
               must_be_real = .true.
               mdlcgg = re
     elseif (name .eq. "FRBLOCK&&7") then
               must_be_real = .true.
               mdlcZga = re
     elseif (name .eq. "FRBLOCK&&8") then
               must_be_real = .true.
               mdlcgghh = re
     elseif (name .eq. "FRBLOCK&&9") then
               must_be_real = .true.
               mdlcVhh = re
     elseif (name .eq. "FRBLOCK&&10") then
               must_be_real = .true.
               mdlcthh = re
     elseif (name .eq. "FRBLOCK&&11") then
               must_be_real = .true.
               mdlchhh = re
     elseif (name .eq. "FRBLOCK&&12") then
               must_be_real = .true.
               mdlchhhh = re
     elseif (name .eq. "SMINPUTS&&1") then
               must_be_real = .true.
               mdlaEWM1 = re
     elseif (name .eq. "SMINPUTS&&2") then
               must_be_real = .true.
               mdlGf = re
     elseif (name .eq. "YUKAWA&&1") then
               must_be_real = .true.
               mdlymdo = re
     elseif (name .eq. "YUKAWA&&2") then
               must_be_real = .true.
               mdlymup = re
     elseif (name .eq. "YUKAWA&&3") then
               must_be_real = .true.
               mdlyms = re
     elseif (name .eq. "YUKAWA&&4") then
               must_be_real = .true.
               mdlymc = re
     elseif (name .eq. "YUKAWA&&6") then
               must_be_real = .true.
               mdlymt = re
     elseif (name .eq. "YUKAWA&&11") then
               must_be_real = .true.
               mdlyme = re
     elseif (name .eq. "YUKAWA&&13") then
               must_be_real = .true.
               mdlymm = re
     elseif (name .eq. "YUKAWA&&15") then
               must_be_real = .true.
               mdlymtau = re
     elseif (name .eq. "MASS&&23") then
               must_be_real = .true.
               mdlMZ = re
     elseif (name .eq. "MASS&&13") then
               must_be_real = .true.
               mdlMMU = re
     elseif (name .eq. "MASS&&6") then
               must_be_real = .true.
               mdlMT = re
     elseif (name .eq. "MASS&&15") then
               must_be_real = .true.
               mdlMTA = re
      elseif (name .eq. "renormalisation") then
          if ( real(int(re),ki) == re .and. im == 0.0_ki ) then
             renormalisation = int(re)
          else
             ierr=0 !FAIL
          end if
      elseif (name .eq. "nlo_prefactors") then
         if ( real(int(re),ki) == re .and. im == 0.0_ki ) then
            nlo_prefactors = int(re)
         else
            ierr=0 !FAIL
         end if
      elseif (name .eq. "deltaOS") then
         if ( real(int(re),ki) == re .and. im == 0.0_ki ) then
            deltaOS = int(re)
         else
            ierr=0 !FAIL
         end if
      elseif (name .eq. "reduction_interoperation") then
         if ( real(int(re),ki) == re .and. im == 0.0_ki ) then
            reduction_interoperation = int(re)
         else
            ierr=0 !FAIL
         end if
     elseif (any(names .eq. name)) then
         do nidx=1,size(names)
            if (names(nidx) .eq. name) exit
         end do
         select case (nidx)
         case(1)
            gauge1z = cmplx(re, im, ki)
         case(2)
            gauge2z = cmplx(re, im, ki)
         case(3)
            gauge5z = cmplx(re, im, ki)
         case(4)
            mdlaEWM1 = re
            must_be_real=.true.
         case(5)
            mdlcb = re
            must_be_real=.true.
         case(6)
            mdlcgaga = re
            must_be_real=.true.
         case(7)
            mdlcgg = re
            must_be_real=.true.
         case(8)
            mdlcgghh = re
            must_be_real=.true.
         case(9)
            mdlchhh = re
            must_be_real=.true.
         case(10)
            mdlchhhh = re
            must_be_real=.true.
         case(11)
            mdlct = re
            must_be_real=.true.
         case(12)
            mdlctau = re
            must_be_real=.true.
         case(13)
            mdlcthh = re
            must_be_real=.true.
         case(14)
            mdlcV = re
            must_be_real=.true.
         case(15)
            mdlcVhh = re
            must_be_real=.true.
         case(16)
            mdlcZga = re
            must_be_real=.true.
         case(17)
            mdlGf = re
            must_be_real=.true.
         case(18)
            mdlMMU = re
            must_be_real=.true.
         case(19)
            mdlMT = re
            must_be_real=.true.
         case(20)
            mdlMTA = re
            must_be_real=.true.
         case(21)
            mdlMZ = re
            must_be_real=.true.
         case(22)
            mdlymc = re
            must_be_real=.true.
         case(23)
            mdlymdo = re
            must_be_real=.true.
         case(24)
            mdlyme = re
            must_be_real=.true.
         case(25)
            mdlymm = re
            must_be_real=.true.
         case(26)
            mdlyms = re
            must_be_real=.true.
         case(27)
            mdlymt = re
            must_be_real=.true.
         case(28)
            mdlymtau = re
            must_be_real=.true.
         case(29)
            mdlymup = re
            must_be_real=.true.
         case(30)
            NC = re
            must_be_real=.true.
         case(31)
            Nf = re
            must_be_real=.true.
         case(32)
            Nfgen = re
            must_be_real=.true.
         end select
      elseif (name .eq. "alpha") then
            mdlaEWM1 = 1._ki/re
            must_be_real=.true.
      elseif (name .eq. "alphaEW") then
            mdlaEWM1 = 1._ki/re
            must_be_real=.true.
      elseif (name .eq. "GF") then
            mdlGf = re
            must_be_real=.true.
     else
         if (name(1:3) /= "mdl") then
            call set_parameter("mdl" // name(4:),re,im,ierr)
            return
         end if
         ierr = 0 !FAIL / UNKNOWN
     end if
     if (must_be_real .and. im /= 0.0_ki .and. ierr.eq.1) then
        ierr = 0 ! FAIL
     end if
     call init_functions()
      ! TODO init_color
   end subroutine
!---#] subroutine set_parameter
!---#[ subroutine init_functions:
   subroutine     init_functions()
      implicit none
      complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
      real(ki), parameter :: pi = 3.14159265358979323846264&
     &3383279502884197169399375105820974944592307816406286209_ki
      mdlloop=16.0_ki*pi**2
      mdlMh=mdlfloat4
      mdlG=2.0_ki*sqrt(pi)*sqrt(mdlaS)
      mdlGC7=i_*mdlG
      mdlGC6=-mdlG
      mdlGC8=i_*mdlG**2
      mdlCKM2x3=0.0_ki
      mdlCKM2x2=1.0_ki
      mdlCKM2x1=0.0_ki
      mdlCKM1x1=1.0_ki
      mdlCKM1x2=0.0_ki
      mdlCKM1x3=0.0_ki
      Nfrat=ifpos(Nfgen,Nf*Nfgen**(-1),1.0_ki)
      mdlCKM3x2=0.0_ki
      mdlCKM3x3=1.0_ki
      mdlCKM3x1=0.0_ki
      mdlaEW=mdlaEWM1**(-1)
      mdlMW=sqrt(mdlfloat1**(-1)*mdlMZ**2+sqrt(mdlfloat2**(-1)*mdlMZ**4-1.0_ki/&
      &(sqrt(2.0_ki))*pi*mdlaEW*mdlGf**(-1)*mdlMZ**2))
      mdlee=2.0_ki*sqrt(pi)*sqrt(mdlaEW)
      mdlGC4=i_*mdlee
      mdlGC5=i_*mdlee**2
      mdlGC3=-i_*mdlee
      mdlGC2=2.0_ki*i_*mdlee*mdlfloat3**(-1)
      mdlGC1=-i_*mdlee*mdlfloat3**(-1)
      mdlsw2=1.0_ki-mdlMW**2*mdlMZ**(-2)
      mdlsw=sqrt(mdlsw2)
      mdlGC53=1.0_ki/(sqrt(2.0_ki))*conjg(mdlCKM2x3)*i_*mdlee*mdlsw**(-1)
      mdlGC52=1.0_ki/(sqrt(2.0_ki))*conjg(mdlCKM2x2)*i_*mdlee*mdlsw**(-1)
      mdlGC51=1.0_ki/(sqrt(2.0_ki))*conjg(mdlCKM2x1)*i_*mdlee*mdlsw**(-1)
      mdlGC50=1.0_ki/(sqrt(2.0_ki))*conjg(mdlCKM1x3)*i_*mdlee*mdlsw**(-1)
      mdlGC56=1.0_ki/(sqrt(2.0_ki))*conjg(mdlCKM3x3)*i_*mdlee*mdlsw**(-1)
      mdlGC55=1.0_ki/(sqrt(2.0_ki))*conjg(mdlCKM3x2)*i_*mdlee*mdlsw**(-1)
      mdlGC54=1.0_ki/(sqrt(2.0_ki))*conjg(mdlCKM3x1)*i_*mdlee*mdlsw**(-1)
      mdlGC48=1.0_ki/(sqrt(2.0_ki))*conjg(mdlCKM1x1)*i_*mdlee*mdlsw**(-1)
      mdlGC49=1.0_ki/(sqrt(2.0_ki))*conjg(mdlCKM1x2)*i_*mdlee*mdlsw**(-1)
      mdlv=2.0_ki*mdlMW*mdlee**(-1)*mdlsw
      mdlyyt=mdlv**(-1)*mdlymt
      mdlyyb=mdlv**(-1)*mdlymb
      mdlGC40=-2.0_ki*i_*mdlloop**(-1)*mdlG**4*mdlv**(-1)*mdlcgg
      mdlGC41=2.0_ki*i_*mdlMW**2*mdlv**(-1)*mdlcV
      mdlGC42=2.0_ki*i_*mdlv**(-1)*mdlcV*mdlMZ**2
      mdlGC44=-i_*mdlyyb*mdlcb
      mdlGC45=-i_*mdlyyt*mdlct
      mdlGC46=-2.0_ki*i_*mdlv**(-1)*mdlyyt*mdlcthh
      mdlyytau=mdlv**(-1)*mdlymtau
      mdlGC47=-i_*mdlyytau*mdlctau
      mdlGC17=1.0_ki/(sqrt(2.0_ki))*i_*mdlCKM2x2*mdlee*mdlsw**(-1)
      mdlGC16=1.0_ki/(sqrt(2.0_ki))*i_*mdlCKM2x1*mdlee*mdlsw**(-1)
      mdlGC15=1.0_ki/(sqrt(2.0_ki))*i_*mdlCKM1x3*mdlee*mdlsw**(-1)
      mdlGC14=1.0_ki/(sqrt(2.0_ki))*i_*mdlCKM1x2*mdlee*mdlsw**(-1)
      mdlGC13=1.0_ki/(sqrt(2.0_ki))*i_*mdlCKM1x1*mdlee*mdlsw**(-1)
      mdlGC12=1.0_ki/(sqrt(2.0_ki))*i_*mdlee*mdlsw**(-1)
      mdlGC10=-i_*mdlee**2*mdlsw**(-2)
      mdlGC19=1.0_ki/(sqrt(2.0_ki))*i_*mdlCKM3x1*mdlee*mdlsw**(-1)
      mdlGC18=1.0_ki/(sqrt(2.0_ki))*i_*mdlCKM2x3*mdlee*mdlsw**(-1)
      mdlGC33=-4.0_ki*i_*mdlloop**(-1)*mdlG**4*mdlv**(-2)*mdlcgghh
      mdlcw=sqrt(1.0_ki-mdlsw2)
      mdlGC11=i_*mdlee**2*mdlsw**(-2)*mdlcw**2
      mdlGC24=i_*mdlee*mdlsw*mdlcw**(-1)*mdlfloat3**(-1)
      mdlGC34=2.0_ki*i_*mdlMW**2*mdlv**(-2)*mdlcVhh
      mdlgw=mdlee*mdlsw**(-1)
      mdlGC31=4.0_ki*i_*mdlloop**(-1)*mdlG**2*mdlv**(-2)*mdlcgghh
      mdlGC30=i_*mdlee*mdlsw**(-1)*mdlcw*mdlfloat1**(-1)+i_*mdlee*mdlsw*mdlcw**&
      &(-1)*mdlfloat1**(-1)
      mdlnormh4=mdlMh**2*mdlv**(-2)*mdlfloat5**(-1)
      mdlGC9=-24.0_ki*i_*mdlnormh4*mdlchhhh
      mdlGC32=4.0_ki*mdlloop**(-1)*mdlG**3*mdlv**(-2)*mdlcgghh
      mdlGC35=2.0_ki*i_*mdlv**(-2)*mdlcVhh*mdlMZ**2
      mdlnormh3=mdlMh**2*mdlv**(-2)*mdlfloat1**(-1)
      mdlGC43=-6.0_ki*i_*mdlv*mdlnormh3*mdlchhh
      mdlGC37=2.0_ki*i_*mdlloop**(-1)*mdlee**2*mdlv**(-1)*mdlcw**(-1)*mdlcZga
      mdlGC36=4.0_ki*i_*mdlloop**(-1)*mdlee**2*mdlv**(-1)*mdlcgaga
      mdlGC39=2.0_ki*mdlloop**(-1)*mdlG**3*mdlv**(-1)*mdlcgg
      mdlGC38=2.0_ki*i_*mdlloop**(-1)*mdlG**2*mdlv**(-1)*mdlcgg
      mdlZERO=mdlfloat6
      mdlg1=mdlee*mdlcw**(-1)
      mdlGC28=i_*mdlee*mdlsw**(-1)*mdlcw*mdlfloat1**(-1)-i_*mdlee*mdlsw*mdlcw**&
      &(-1)*mdlfloat7**(-1)
      mdlGC29=-i_*mdlee*mdlsw**(-1)*mdlcw*mdlfloat1**(-1)+i_*mdlee*mdlsw*mdlcw*&
      &*(-1)*mdlfloat1**(-1)
      mdlGC26=i_*mdlee*mdlsw*mdlcw**(-1)
      mdlGC27=-i_*mdlee*mdlsw**(-1)*mdlcw*mdlfloat1**(-1)-i_*mdlee*mdlsw*mdlcw*&
      &*(-1)*mdlfloat7**(-1)
      mdlvev=2.0_ki*mdlMW*mdlee**(-1)*mdlsw
      mdlydo=sqrt(2.0_ki)*mdlvev**(-1)*mdlymdo
      mdlytau=sqrt(2.0_ki)*mdlvev**(-1)*mdlymtau
      mdlyup=sqrt(2.0_ki)*mdlvev**(-1)*mdlymup
      mdlye=sqrt(2.0_ki)*mdlvev**(-1)*mdlyme
      mdlyc=sqrt(2.0_ki)*mdlvev**(-1)*mdlymc
      mdlyb=sqrt(2.0_ki)*mdlvev**(-1)*mdlymb
      mdlym=sqrt(2.0_ki)*mdlvev**(-1)*mdlymm
      mdlyt=sqrt(2.0_ki)*mdlvev**(-1)*mdlymt
      mdlys=sqrt(2.0_ki)*mdlvev**(-1)*mdlyms
      mdlGC25=-2.0_ki*i_*mdlee*mdlsw*mdlcw**(-1)*mdlfloat3**(-1)
      mdlGC22=i_*mdlee*mdlsw**(-1)*mdlcw
      mdlGC23=-2.0_ki*i_*mdlee**2*mdlsw**(-1)*mdlcw
      mdlGC20=1.0_ki/(sqrt(2.0_ki))*i_*mdlCKM3x2*mdlee*mdlsw**(-1)
      mdlGC21=1.0_ki/(sqrt(2.0_ki))*i_*mdlCKM3x3*mdlee*mdlsw**(-1)
end subroutine init_functions
!---#] subroutine init_functions:
!---#[ utility functions for model initialization:
   pure function ifpos(x0, x1, x2)
      implicit none
      real(ki), intent(in) :: x0, x1, x2
      real(ki) :: ifpos
      if (x0 > 0.0_ki) then
         ifpos = x1
      else
         ifpos = x2
      endif
   end  function ifpos
   pure function sort4(m1, m2, m3, m4, n)
      implicit none
      real(ki), intent(in) :: m1, m2, m3, m4
      integer, intent(in) :: n
      real(ki) :: sort4
      real(ki), dimension(4) :: m
      logical :: f
      integer :: i
      real(ki) :: tmp
      m(1) = m1
      m(2) = m2
      m(3) = m3
      m(4) = m4
      ! Bubble Sort
      do
         f = .false.
         do i=1,3
            if (abs(m(i)) .gt. abs(m(i+1))) then
               tmp = m(i)
               m(i) = m(i+1)
               m(i+1) = tmp
               f = .true.
            end if
         end do
         if (.not. f) exit
      end do
      sort4 = m(n)
   end  function sort4
!---#] utility functions for model initialization:
end module p2_part21part21_part25part25part21_model_qp
