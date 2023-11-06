module     p2_part21part21_part25part25part21_matrix
   use p2_part21part21_part25part25part21_util, only: square
   use p2_part21part21_part25part25part21_util_qp, only: square_qp => square
   use p2_part21part21_part25part25part21_config, only: ki, ki_qp, &
     & include_helicity_avg_factor, include_color_avg_factor, &
     & debug_lo_diagrams, debug_nlo_diagrams, &
     & include_symmetry_factor, &
     & PSP_check, PSP_verbosity, PSP_rescue, PSP_chk_th1, &
     & PSP_chk_th2, PSP_chk_th3, PSP_chk_kfactor, reduction_interoperation, &
     & PSP_chk_li1, PSP_chk_li2, PSP_chk_li3, PSP_chk_li4, &
     & reduction_interoperation_rescue, convert_to_cdr
   use p2_part21part21_part25part25part21_kinematics, only: &
       in_helicities, symmetry_factor, num_legs, &
       lo_qcd_couplings, corrections_are_qcd, num_light_quarks, num_gluons
   use p2_part21part21_part25part25part21_model, only: Nf, NC, sqrt2, init_functions
   use p2_part21part21_part25part25part21_model_qp, only: Nf_qp => Nf, NC_qp => NC, sqrt2_qp => sqrt2, &
     & init_functions_qp => init_functions
   use p2_part21part21_part25part25part21_color, only: TR, CA, CF, numcs, &
     & incolors, init_color
   use p2_part21part21_part25part25part21_color_qp, only: TR_qp => TR, CA_qp => CA, CF_qp => CF, &
     & incolors_qp => incolors, init_color_qp => init_color
   use p2_part21part21_part25part25part21_amplitudeh0, only: samplitudeh0l1 => samplitude, &
     &   finite_renormalisation0 => finite_renormalisation
   use p2_part21part21_part25part25part21_amplitudeh0_qp, only: samplitudeh0l1_qp => samplitude, &
     &   finite_renormalisation0_qp => finite_renormalisation
   use p2_part21part21_part25part25part21_dipoles, only: insertion_operator, insertion_operator_qed
   use p2_part21part21_part25part25part21_dipoles_qp, only: insertion_operator_qp => insertion_operator, &
     & insertion_operator_qed_qp => insertion_operator_qed

   implicit none
   save

   private

   integer :: banner_ch = 6

   public :: initgolem, exitgolem, samplitude
   public :: samplitudel0, samplitudel1
   public :: ir_subtraction, color_correlated_lo2, spin_correlated_lo2
   public :: OLP_color_correlated, OLP_spin_correlated_lo2




contains
   !---#[ subroutine banner:
   subroutine     banner()
      implicit none

      character(len=74) :: frame = "+" // repeat("-", 72) // "+"

      if (banner_ch .le. 0) return

      write(banner_ch,'(A74)') frame
      write(banner_ch,'(A74)') "|   __   __   ___   __   __  __                   GoSam                  |"
      write(banner_ch,'(A74)') "|  / _) /  \ / __) (  ) (  \/  )          An Automated One-Loop          |"
      write(banner_ch,'(A74)') "| ( (/\( () )\__ \ /__\  )    (          Matrix Element Generator        |"
      write(banner_ch,'(A74)') "|  \__/ \__/ (___/(_)(_)(_/\/\_)        Version 2.0.4 Rev: 1fe5178       |"
      write(banner_ch,'(A74)') "|                                                                        |"
      write(banner_ch,'(A74)') "|                                 (c) The GoSam Collaboration 2011-2016  |"
      write(banner_ch,'(A74)') "|                                                                        |"
      write(banner_ch,'(A74)') "|  AUTHORS:                                                              |"
      write(banner_ch,'(A74)') "|  * Nicolas Greiner                   <greiner@mpp.mpg.de>              |"
      write(banner_ch,'(A74)') "|  * Gudrun Heinrich                   <gudrun@mpp.mpg.de>               |"
      write(banner_ch,'(A74)') "|  * Stephan Jahn                      <sjahn@mpp.mpg.de>                |"
      write(banner_ch,'(A74)') "|  * Gionata Luisoni                   <luisonig@mpp.mpg.de>             |"
      write(banner_ch,'(A74)') "|  * Pierpaolo Mastrolia               <Pierpaolo.Mastrolia@cern.ch>     |"
      write(banner_ch,'(A74)') "|  * Giovanni Ossola                   <gossola@citytech.cuny.edu>       |"
      write(banner_ch,'(A74)') "|  * Tiziano Peraro                    <peraro@mpp.mpg.de>               |"
      write(banner_ch,'(A74)') "|  * Johannes Schlenk                  <jschlenk@mpp.mpg.de>             |"
      write(banner_ch,'(A74)') "|  * Johann Felix von Soden-Fraunhofen <jfsoden@mpp.mpg.de>              |"
      write(banner_ch,'(A74)') "|  * Francesco Tramontano              <Francesco.Tramontano@cern.ch>    |"
      write(banner_ch,'(A74)') "|                                                                        |"
      write(banner_ch,'(A74)') "|  FORMER AUTHORS:                                                       |"
      write(banner_ch,'(A74)') "|  * Gavin Cullen                      <gavin.cullen@desy.de>            |"
      write(banner_ch,'(A74)') "|  * Hans van Deurzen                  <hdeurzen@mpp.mpg.de>             |"
      write(banner_ch,'(A74)') "|  * Edoardo Mirabella                 <mirabell@mpp.mpg.de>             |"
      write(banner_ch,'(A74)') "|  * Joscha Reichel                    <joscha@mpp.mpg.de>               |"
      write(banner_ch,'(A74)') "|  * Thomas Reiter                     <reiterth@mpp.mpg.de>             |"
      write(banner_ch,'(A74)') "|                                                                        |"
      write(banner_ch,'(A74)') "|  This program is free software: you can redistribute it and/or modify  |"
      write(banner_ch,'(A74)') "|  it under the terms of the GNU General Public License either           |"
      write(banner_ch,'(A74)') "|  version 3, or (at your option) any later version.                     |"
      write(banner_ch,'(A74)') "|                                                                        |"
      write(banner_ch,'(A74)') "|  Scientific publications prepared using the present version of         |"
      write(banner_ch,'(A74)') "|  GoSam or any modified version of it or any code linking to GoSam      |"
      write(banner_ch,'(A74)') "|  or parts of it should make a clear reference to the publication:      |"
      write(banner_ch,'(A74)') "|                                                                        |"
      write(banner_ch,'(A74)') "|      G. Cullen et al.,                                                 |"
      write(banner_ch,'(A74)') "|      ``GoSam-2.0: a tool for automated one-loop calculations           |"
      write(banner_ch,'(A74)') "|                        within the Standard Model and Beyond'',         |"
      write(banner_ch,'(A74)') "|      Eur. Phys. J. C 74 (2014) 8,  3001                                |"
      write(banner_ch,'(A74)') "|      [arXiv:1404.7096 [hep-ph]].                                       |"
      write(banner_ch,'(A74)') frame

      banner_ch = 0
   end subroutine banner
   !---#] subroutine banner:

   !---#[ subroutine initgolem :
   subroutine     initgolem(is_first,stage,rndseed)
      implicit none
      logical, optional, intent(in) :: is_first
      integer, optional, intent(in) :: stage
      integer, optional, intent(in) :: rndseed
      logical :: init_third_party
      logical :: file_exists, dir_exists
      integer i, j
      character(len=50) :: file_name
      character(len=9)  :: dir_name = "BadPoints"
      character(len=6)  :: file_numb
      character(len=9)  :: file_pre = "gs_badpts"
      character(len=3)  :: file_ext = "log"
      character(len=1)  :: cstage
      character(len=4)  :: crndseed
      i = 1
      file_exists =.true.

      if(present(is_first)) then
         init_third_party = is_first
      else
         init_third_party = .true.
      end if

      if(.not. corrections_are_qcd) then
         PSP_check = .false.
      end if
      if (init_third_party) then
      ! call our banner
      call banner()
      if(PSP_check.and.PSP_rescue.and.PSP_verbosity) then
         inquire(file=dir_name, exist=dir_exists)
         if(.not. dir_exists) then
            call system('mkdir BadPoints')
         end if
         if(present(stage)) then
            write(cstage,'(i1)') stage
            write(crndseed,'(i4)') rndseed
            do j=1,4
               if(crndseed(j:j).eq.' ') crndseed(j:j)='0'
            enddo
            file_name = dir_name//"/"//file_pre//"-"//cstage//"-"//crndseed//"."//file_ext
            open(unit=42, file=file_name, status='replace', action='write')
            write(42,'(A22)') "<?xml version='1.0' ?>"
            write(42,'(A5)')  "<run>"
         else
            do while(file_exists)
               write(file_numb, '(I6.1)') i
               file_name = dir_name//"/"//file_pre//trim(adjustl(file_numb))//"."//file_ext
               inquire(file=file_name, exist=file_exists)
               if(file_exists) then
                  write(*,*) "File ", file_name, " already exists!"
                  i = i+1
               else
                  write(*,*) "Bad points stored in file: ", file_name
                  open(unit=42, file=file_name, status='unknown', action='write')
                  write(42,'(A22)') "<?xml version='1.0' ?>"
                  write(42,'(A5)')  "<run>"
               end if
            enddo
         end if
      end if
      end if

      call init_functions()
      call init_color()

      call init_functions_qp()
      call init_color_qp()

   end subroutine initgolem
   !---#] subroutine initgolem :
   !---#[ subroutine exitgolem :
   subroutine     exitgolem(is_last)
      implicit none
      logical, optional, intent(in) :: is_last

      logical :: exit_third_party

      if(present(is_last)) then
         exit_third_party = is_last
      else
         exit_third_party = .true.
      end if
      if (exit_third_party) then
         if(PSP_check.and.PSP_rescue.and.PSP_verbosity) then
            write(42,'(A6)')  "</run>"
            close(unit=42)
         endif
      end if
   end subroutine exitgolem
   !---#] subroutine exitgolem :

   !---#[ subroutine samplitude :
   subroutine     samplitude(vecs, scale2, amp, prec, ichecked, ok, h)
      use p2_part21part21_part25part25part21_kinematics_qp, only: adjust_kinematics_qp => adjust_kinematics
      use p2_part21part21_part25part25part21_model
      implicit none
      real(ki), dimension(5, 4), intent(in) :: vecs
      real(ki_qp), dimension(5, 4) :: vecs_qp
      real(ki), dimension(5, 4) :: vecsrot, vecsrot2
      real(ki), intent(in) :: scale2
      real(ki), dimension(1:4), intent(out) :: amp
      real(ki_qp), dimension(1:4) :: amp_qp
      real(ki_qp), dimension(2:3) :: irp_qp
      real(ki_qp) :: scale2_qp, rat2_qp
      real(ki), dimension(0:2) :: scales2
      real(ki), dimension(1:4) :: ampdef, amprot, ampres, ampresrot
      real(ki) :: rat2, kfac, zero, angle
      real(ki), dimension(2:3) :: irp
      integer, intent(out) :: prec
      integer, intent(out) :: ichecked
      logical, intent(out), optional :: ok
      integer, intent(in), optional :: h
      integer spprec1, fpprec1, spprec2, fpprec2
      integer tmp_red_int, icheck, i, irot
      ampdef=0.0_ki
      amprot=0.0_ki
      ampres=0.0_ki
      ampresrot=0.0_ki
      icheck = 1
      angle = 1.234_ki
      spprec1 = 18
      spprec2 = 18
      fpprec1 = 18
      fpprec2 = 18
      prec = 20 ! Start with unrealistically high value
      scales2(:) = (/0.0_ki, &
     &              mdlMh, &
     &              mdlMh/)
      tmp_red_int = reduction_interoperation
      if(reduction_interoperation.eq.4) then
         PSP_check=.true.
         PSP_rescue=.true.
         icheck = 3
      else
         call samplitudel01(vecs, scale2, ampdef, rat2, ok, h)
         amp = ampdef
      endif
      ! RESCUE SYSTEM
      if(PSP_check) then
         ! CHECK ON ROTATED AMPLITUDE
         do irot = 1,5
            vecsrot(irot,1) = vecs(irot,1)
            vecsrot(irot,2) = vecs(irot,2)*Cos(angle)-vecs(irot,3)*Sin(angle)
            vecsrot(irot,3) = vecs(irot,2)*Sin(angle)+vecs(irot,3)*Cos(angle)
            vecsrot(irot,4) = vecs(irot,4)
         enddo
         do irot = 1,5
            vecsrot2(irot,1) = vecsrot(irot,1)
            vecsrot2(irot,2) = vecsrot(irot,2)*Cos(angle)-vecsrot(irot,4)*Sin(angle)
            vecsrot2(irot,3) = vecsrot(irot,3)
            vecsrot2(irot,4) = vecsrot(irot,2)*Sin(angle)+vecsrot(irot,4)*Cos(angle)
         enddo
         !call adjust_kinematics(vecsrot)
         call samplitudel01(vecsrot2, scale2, ampresrot, rat2, ok, h)
         if((ampresrot(2)-ampdef(2)) .ne. 0.0_ki) then
            fpprec1 = -int(log10(abs((ampresrot(2)-ampdef(2))/((ampresrot(2)+ampdef(2))/2.0_ki))))
         else
            fpprec1 = 16
         endif
         kfac = 0.0_ki
         if(fpprec1.lt.PSP_chk_li2) then                                       ! RESCUE
            icheck=3
            fpprec1=-10        ! Set -10 as finite part precision
         endif
         prec = min(spprec1,fpprec1)
         if(icheck.eq.3.and.PSP_rescue) then
            !icheck = 1
            reduction_interoperation = reduction_interoperation_rescue
            scale2_qp = real(scale2,ki_qp)
            !call refine_momenta_to_qp(5,vecs,vecs_qp,2+1,scales2)
            vecs_qp = vecs
            call adjust_kinematics_qp(vecs_qp)
            call samplitudel01_qp(vecs_qp, scale2_qp, amp_qp, rat2_qp, ok, h)
            ampres = real(amp_qp,ki)
            amp=ampres
            ! poles should be zero for loop-induced processes
            if(ampres(2) .ne. 0.0_ki .and. ampres(3) .ne. 0.0_ki) then
               spprec2 = -int(log10(abs(ampres(3)/ampres(2))))
            else
               spprec2 = 16
            endif
            kfac = 0.0_ki
            if(spprec2.lt.PSP_chk_li4) then ! DISCARD
               icheck=3
               fpprec2=-10        ! Set -10 as finite part precision
            endif
            !if(icheck.eq.2) then
            !   do irot = 1,5
            !      vecsrot_qp(irot,1) = vecs_qp(irot,1)
            !      vecsrot_qp(irot,2) = vecs_qp(irot,2)*Cos(angle_qp)-vecs_qp(irot,3)*Sin(angle_qp)
            !      vecsrot_qp(irot,3) = vecs_qp(irot,2)*Sin(angle_qp)+vecs_qp(irot,3)*Cos(angle_qp)
            !      vecsrot_qp(irot,4) = vecs_qp(irot,4)
            !   enddo
            !   call adjust_kinematics_qp(vecsrot_qp)
            !   call samplitudel01_qp(vecsrot_qp, scale2_qp, ampresrot_qp, rat2_qp, ok, h)
            !   if((ampresrot_qp(2)-ampres(2)) .ne. 0.0_ki) then
            !      fpprec2 = -int(log10(abs((ampresrot_qp(2)-ampres(2))/((ampresrot_qp(2)+ampres(2))/2.0_ki))))
            !   else
            !      fpprec2 = 16
            !   endif
            !   if(fpprec2.ge.PSP_chk_li3) icheck=1                         ! ACCEPTED
            !   if(fpprec2.lt.PSP_chk_li3) icheck=3                         ! DISCARD
            !endif
            reduction_interoperation = tmp_red_int
            prec = min(spprec2,fpprec2)
         endif

         if(icheck.eq.3.and.PSP_verbosity) then
            write(42,'(2x,A7)')"<event>"
            write(42,'(4x,A15,A34,A3)') &
                 &  "<process name='","p2_part21part21_part25part25part21","'/>"
            write(42,'(4x,A21,I2.1,A7,I2.1,A7,I2.1,A7,I2.1,A3)') &
                 &  "<PSP_thresholds li1='", PSP_chk_li1, &
                 &                "' li2='", PSP_chk_li2, &
                 &                "' li3='", PSP_chk_li3, &
                 &                "' li4='", PSP_chk_li4,"'/>"
            write(42,'(4x,A16,D23.16,A3)') &
                 &  "<PSP_kfaktor k='", PSP_chk_kfactor,"'/>"
            write(42,'(4x,A15,I3.1,A6,I3.1,A3)') &
                 &  "<PSP_prec1 sp='", spprec1, "' fp='", fpprec1, "'/>"
            write(42,'(4x,A15,I3.1,A6,I3.1,A3)') &
                 &  "<PSP_prec2 sp='", spprec2, "' fp='", fpprec2, "'/>"
            write(42,'(4x,A10,D23.16,A3)') &
                 &  "<born LO='", ampdef(1), "'/>"
            write(42,'(4x,A10,D23.16,A3)') &
                 &  "<rat2 r2='", rat2, "'/>"
            write(42,'(4x,A15,D23.16,A6,D23.16,A6,D23.16,A3)') &
                 &  "<amp       sp='", ampdef(3)   ,"' ir='", irp(2),"' fp='", ampdef(2)   ,"'/>"
            write(42,'(4x,A15,D23.16,A6,D23.16,A6,D23.16,A3)') &
                 &  "<amprot    sp='", amprot(3)   ,"' ir='", irp(2),"' fp='", amprot(2)   ,"'/>"
            write(42,'(4x,A15,D23.16,A6,D23.16,A6,D23.16,A3)') &
                 &  "<ampres    sp='", ampres(3)   ,"' ir='", irp(2),"' fp='", ampres(2)   ,"'/>"
            write(42,'(4x,A15,D23.16,A6,D23.16,A6,D23.16,A3)') &
                 &  "<ampresrot sp='", ampresrot(3),"' ir='", irp(2),"' fp='", ampresrot(2),"'/>"
            write(42,'(4x,A9)') "<momenta>"
            do i=1,5
               write(42,'(8x,A8,3(D23.16,A6),D23.16,A3)') "<mom e='", vecs(i,1), &
                    &  "' px='", vecs(i,2), &
                    &  "' py='", vecs(i,3), &
                    &  "' pz='", vecs(i,4), "'/>"
            enddo
            write(42,'(4x,A10)')"</momenta>"
            write(42,'(2x,A8)')"</event>"
         endif
      else
         prec = 20 ! If PSP_check is off, precision is set to unrealistic value = 20.
      end if
      ichecked = icheck
 end subroutine samplitude
   !---#] subroutine samplitude :

   !---#[ subroutine samplitudel01 :
   subroutine     samplitudel01(vecs, scale2, amp, rat2, ok, h)
      use p2_part21part21_part25part25part21_config, only: &
         & debug_lo_diagrams, debug_nlo_diagrams, logfile, deltaOS, &
         & renormalisation, renorm_beta, renorm_mqwf, renorm_decoupling, &
         & renorm_logs, renorm_mqse, renorm_yukawa, nlo_prefactors
      use p2_part21part21_part25part25part21_kinematics, only: &
         & inspect_kinematics, init_event
      use p2_part21part21_part25part25part21_model
      use p2_part21part21_part25part25part21_dipoles, only: pi
      implicit none
      real(ki), dimension(5, 4), intent(in) :: vecs
      real(ki), intent(in) :: scale2
      real(ki), dimension(4), intent(out) :: amp
      real(ki), intent(out) :: rat2
      logical, intent(out), optional :: ok
      integer, intent(in), optional :: h
      real(ki) :: nlo_coupling

      complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)

      ! Number of heavy quark flavours in loops.
      real(ki), parameter :: NFh = 0.0_ki

      logical :: my_ok

      ! used for m=0 QCD renormalisation
      real(ki) :: beta0

      if(corrections_are_qcd) then
         nlo_coupling = mdlG*mdlG
      else
         nlo_coupling = mdlee*mdlee
      end if

      if(debug_lo_diagrams .or. debug_nlo_diagrams) then
         call init_event(vecs)
         write(logfile,'(A7)') "<event>"
         call inspect_kinematics(logfile)
      end if


      amp(1)   = 0.0_ki
      select case (renormalisation)
      case (0)
         ! no renormalisation
         deltaOS = 0.0_ki
      case (1)
         ! fully renormalized
         if(renorm_mqse) then
            deltaOS = 1.0_ki
         else
            deltaOS = 0.0_ki
         end if
      case (2)
         ! massive quark counterterms only
         deltaOS = 1.0_ki
      case default
         ! not implemented
         print*, "In p2_part21part21_part25part25part21_matrix:"
         print*, "  invalid value for renormalisation=", renormalisation
         stop
      end select

      if (present(h)) then
         amp((/4,3,2/)) = samplitudel1(vecs, scale2, my_ok, rat2, h)/nlo_coupling/nlo_coupling
      else
         amp((/4,3,2/)) = samplitudel1(vecs, scale2, my_ok, rat2)/nlo_coupling/nlo_coupling
      end if
      select case (renormalisation)
      case (0)
         ! no renormalisation
      case (1)
         ! fully renormalized
         ! No tree level present
      case (2)
         ! massive quark counterterms only
      case default
         ! not implemented
         print*, "In p2_part21part21_part25part25part21_matrix:"
         print*, "  invalid value for renormalisation=", renormalisation
         stop
      end select
      if (convert_to_cdr) then
         ! Scheme conversion for infrared structure
         ! Reference:
         ! S. Catani, M. H. Seymour, Z. Trocsanyi,
         ! ``Regularisation scheme independence and unitarity
         !   in QCD cross-sections,''
         ! Phys.Rev. D 55 (1997) 6819
         ! arXiv:hep-ph/9610553
         amp(2) = amp(2) - amp(1) * (&
           &          num_light_quarks * 0.5_ki * CF &
           &        + num_gluons * 1.0_ki/6.0_ki * CA)
      end if
      if (present(ok)) ok = my_ok

      if(debug_lo_diagrams .or. debug_nlo_diagrams) then
         write(logfile,'(A25,E24.16,A3)') &
            & "<result kind='lo' value='", amp(1), "'/>"
         write(logfile,'(A33,E24.16,A3)') &
            & "<result kind='nlo-finite' value='", amp(2), "'/>"
         write(logfile,'(A33,E24.16,A3)') &
            & "<result kind='nlo-single' value='", amp(3), "'/>"
         write(logfile,'(A33,E24.16,A3)') &
            & "<result kind='nlo-double' value='", amp(4), "'/>"
         if(my_ok) then
            write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
         else
            write(logfile,'(A29)') "<flag name='ok' status='no'/>"
         end if
         write(logfile,'(A8)') "</event>"
      end if
      select case(nlo_prefactors)
      case(0)
         ! The result is already in its desired form
      case(1)
         ! loop-induced
         amp(2:4) = amp(2:4) * nlo_coupling * nlo_coupling
      case(2)
         ! loop-induced
         amp(2:4) = amp(2:4) * (nlo_coupling / 8.0_ki / pi / pi)**2
      end select
   end subroutine samplitudel01
   !---#] subroutine samplitudel01 :
   !---#[ function samplitudel0 :
   function     samplitudel0(vecs, h) result(amp)
      use p2_part21part21_part25part25part21_config, only: logfile
      use p2_part21part21_part25part25part21_kinematics, only: init_event
      implicit none
      real(ki), dimension(5, 4), intent(in) :: vecs
      integer, optional, intent(in) :: h
      real(ki) :: amp, heli_amp
      complex(ki), dimension(numcs) :: color_vector
      logical, dimension(0:7) :: eval_heli
      real(ki), dimension(5, 4) :: pvecs

      if (present(h)) then
         eval_heli(:) = .false.
         eval_heli(h) = .true.
      else
         eval_heli(:) = .true.
      end if

      amp = 0.0_ki
   end function samplitudel0
   !---#] function samplitudel0 :
   !---#[ function samplitudel1 :
   function     samplitudel1(vecs,scale2,ok,rat2,h) result(amp)
      use p2_part21part21_part25part25part21_config, only: &
         & debug_nlo_diagrams, logfile, renorm_gamma5
      use p2_part21part21_part25part25part21_kinematics, only: init_event
      implicit none
      real(ki), dimension(5, 4), intent(in) :: vecs
      logical, intent(out) :: ok
      real(ki), intent(in) :: scale2
      real(ki), intent(out) :: rat2
      integer, optional, intent(in) :: h
      real(ki), dimension(5, 4) :: pvecs
      real(ki), dimension(-2:0) :: amp, heli_amp
      complex(ki), dimension(numcs,-2:0) :: colorvec
      integer :: c
      logical :: my_ok
      logical, dimension(0:7) :: eval_heli
      real(ki) :: fr, rational2

      amp(:) = 0.0_ki
      rat2 = 0.0_ki
      ok = .true.

      if (present(h)) then
         eval_heli(:) = .false.
         eval_heli(h) = .true.
      else
         eval_heli(:) = .true.
      end if


      if (eval_heli(0)) then
         if(debug_nlo_diagrams) then
            write(logfile,*) "<helicity index='0'>"
         end if
         !---#[ reinitialize kinematics:
         pvecs(1,:) = vecs(1,:)
         pvecs(2,:) = vecs(2,:)
         pvecs(3,:) = vecs(3,:)
         pvecs(4,:) = vecs(4,:)
         pvecs(5,:) = vecs(5,:)
         call init_event(pvecs, -1, -1, -1)
            !---#] reinitialize kinematics:
         do c=1,numcs
            colorvec(c,:) = samplitudeh0l1(real(scale2,ki),my_ok,rational2,c)
         end do
         heli_amp( 0) = square(colorvec(:, 0))
         heli_amp(-1) = square(colorvec(:,-1))
         heli_amp(-2) = square(colorvec(:,-2))
      
         if (corrections_are_qcd .and. renorm_gamma5) then
            !---#[ reinitialize kinematics:
            pvecs(1,:) = vecs(1,:)
            pvecs(2,:) = vecs(2,:)
            pvecs(3,:) = vecs(3,:)
            pvecs(4,:) = vecs(4,:)
            pvecs(5,:) = vecs(5,:)
            call init_event(pvecs, -1, -1, -1)
            !---#] reinitialize kinematics:
            fr = finite_renormalisation0(real(scale2,ki))
            heli_amp(0) = heli_amp(0) + fr
         end if
         ok = ok .and. my_ok
         amp = amp + heli_amp
         rat2 = rat2 + rational2

         if(debug_nlo_diagrams) then
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-finite' value='", heli_amp(0), "'/>"
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-single' value='", heli_amp(-1), "'/>"
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-double' value='", heli_amp(-2), "'/>"
            if (corrections_are_qcd .and. renorm_gamma5) then
               write(logfile,'(A30,E24.16,A3)') &
                   & "<result kind='fin-ren' value='", fr, "'/>"
            end if
            if(my_ok) then
               write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
            else
               write(logfile,'(A29)') "<flag name='ok' status='no'/>"
            end if
            write(logfile,*) "</helicity>"
         end if
      end if
      if (eval_heli(1)) then
         if(debug_nlo_diagrams) then
            write(logfile,*) "<helicity index='1'>"
         end if
         !---#[ reinitialize kinematics:
         pvecs(1,:) = vecs(1,:)
         pvecs(2,:) = vecs(2,:)
         pvecs(3,:) = vecs(4,:)
         pvecs(4,:) = vecs(3,:)
         pvecs(5,:) = vecs(5,:)
         call init_event(pvecs, +1, -1, -1)
            !---#] reinitialize kinematics:
         do c=1,numcs
            colorvec(c,:) = samplitudeh0l1(real(scale2,ki),my_ok,rational2,c)
         end do
         heli_amp( 0) = square(colorvec(:, 0))
         heli_amp(-1) = square(colorvec(:,-1))
         heli_amp(-2) = square(colorvec(:,-2))
      
         if (corrections_are_qcd .and. renorm_gamma5) then
            !---#[ reinitialize kinematics:
            pvecs(1,:) = vecs(1,:)
            pvecs(2,:) = vecs(2,:)
            pvecs(3,:) = vecs(4,:)
            pvecs(4,:) = vecs(3,:)
            pvecs(5,:) = vecs(5,:)
            call init_event(pvecs, +1, -1, -1)
            !---#] reinitialize kinematics:
            fr = finite_renormalisation0(real(scale2,ki))
            heli_amp(0) = heli_amp(0) + fr
         end if
         ok = ok .and. my_ok
         amp = amp + heli_amp
         rat2 = rat2 + rational2

         if(debug_nlo_diagrams) then
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-finite' value='", heli_amp(0), "'/>"
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-single' value='", heli_amp(-1), "'/>"
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-double' value='", heli_amp(-2), "'/>"
            if (corrections_are_qcd .and. renorm_gamma5) then
               write(logfile,'(A30,E24.16,A3)') &
                   & "<result kind='fin-ren' value='", fr, "'/>"
            end if
            if(my_ok) then
               write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
            else
               write(logfile,'(A29)') "<flag name='ok' status='no'/>"
            end if
            write(logfile,*) "</helicity>"
         end if
      end if
      if (eval_heli(2)) then
         if(debug_nlo_diagrams) then
            write(logfile,*) "<helicity index='2'>"
         end if
         !---#[ reinitialize kinematics:
         pvecs(1,:) = vecs(1,:)
         pvecs(2,:) = vecs(2,:)
         pvecs(3,:) = vecs(4,:)
         pvecs(4,:) = vecs(3,:)
         pvecs(5,:) = vecs(5,:)
         call init_event(pvecs, -1, +1, -1)
            !---#] reinitialize kinematics:
         do c=1,numcs
            colorvec(c,:) = samplitudeh0l1(real(scale2,ki),my_ok,rational2,c)
         end do
         heli_amp( 0) = square(colorvec(:, 0))
         heli_amp(-1) = square(colorvec(:,-1))
         heli_amp(-2) = square(colorvec(:,-2))
      
         if (corrections_are_qcd .and. renorm_gamma5) then
            !---#[ reinitialize kinematics:
            pvecs(1,:) = vecs(1,:)
            pvecs(2,:) = vecs(2,:)
            pvecs(3,:) = vecs(4,:)
            pvecs(4,:) = vecs(3,:)
            pvecs(5,:) = vecs(5,:)
            call init_event(pvecs, -1, +1, -1)
            !---#] reinitialize kinematics:
            fr = finite_renormalisation0(real(scale2,ki))
            heli_amp(0) = heli_amp(0) + fr
         end if
         ok = ok .and. my_ok
         amp = amp + heli_amp
         rat2 = rat2 + rational2

         if(debug_nlo_diagrams) then
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-finite' value='", heli_amp(0), "'/>"
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-single' value='", heli_amp(-1), "'/>"
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-double' value='", heli_amp(-2), "'/>"
            if (corrections_are_qcd .and. renorm_gamma5) then
               write(logfile,'(A30,E24.16,A3)') &
                   & "<result kind='fin-ren' value='", fr, "'/>"
            end if
            if(my_ok) then
               write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
            else
               write(logfile,'(A29)') "<flag name='ok' status='no'/>"
            end if
            write(logfile,*) "</helicity>"
         end if
      end if
      if (eval_heli(3)) then
         if(debug_nlo_diagrams) then
            write(logfile,*) "<helicity index='3'>"
         end if
         !---#[ reinitialize kinematics:
         pvecs(1,:) = vecs(1,:)
         pvecs(2,:) = vecs(2,:)
         pvecs(3,:) = vecs(4,:)
         pvecs(4,:) = vecs(3,:)
         pvecs(5,:) = vecs(5,:)
         call init_event(pvecs, +1, +1, -1)
            !---#] reinitialize kinematics:
         do c=1,numcs
            colorvec(c,:) = samplitudeh0l1(real(scale2,ki),my_ok,rational2,c)
         end do
         heli_amp( 0) = square(colorvec(:, 0))
         heli_amp(-1) = square(colorvec(:,-1))
         heli_amp(-2) = square(colorvec(:,-2))
      
         if (corrections_are_qcd .and. renorm_gamma5) then
            !---#[ reinitialize kinematics:
            pvecs(1,:) = vecs(1,:)
            pvecs(2,:) = vecs(2,:)
            pvecs(3,:) = vecs(4,:)
            pvecs(4,:) = vecs(3,:)
            pvecs(5,:) = vecs(5,:)
            call init_event(pvecs, +1, +1, -1)
            !---#] reinitialize kinematics:
            fr = finite_renormalisation0(real(scale2,ki))
            heli_amp(0) = heli_amp(0) + fr
         end if
         ok = ok .and. my_ok
         amp = amp + heli_amp
         rat2 = rat2 + rational2

         if(debug_nlo_diagrams) then
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-finite' value='", heli_amp(0), "'/>"
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-single' value='", heli_amp(-1), "'/>"
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-double' value='", heli_amp(-2), "'/>"
            if (corrections_are_qcd .and. renorm_gamma5) then
               write(logfile,'(A30,E24.16,A3)') &
                   & "<result kind='fin-ren' value='", fr, "'/>"
            end if
            if(my_ok) then
               write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
            else
               write(logfile,'(A29)') "<flag name='ok' status='no'/>"
            end if
            write(logfile,*) "</helicity>"
         end if
      end if
      if (eval_heli(4)) then
         if(debug_nlo_diagrams) then
            write(logfile,*) "<helicity index='4'>"
         end if
         !---#[ reinitialize kinematics:
         pvecs(1,:) = vecs(1,:)
         pvecs(2,:) = vecs(2,:)
         pvecs(3,:) = vecs(4,:)
         pvecs(4,:) = vecs(3,:)
         pvecs(5,:) = vecs(5,:)
         call init_event(pvecs, -1, -1, +1)
            !---#] reinitialize kinematics:
         do c=1,numcs
            colorvec(c,:) = samplitudeh0l1(real(scale2,ki),my_ok,rational2,c)
         end do
         heli_amp( 0) = square(colorvec(:, 0))
         heli_amp(-1) = square(colorvec(:,-1))
         heli_amp(-2) = square(colorvec(:,-2))
      
         if (corrections_are_qcd .and. renorm_gamma5) then
            !---#[ reinitialize kinematics:
            pvecs(1,:) = vecs(1,:)
            pvecs(2,:) = vecs(2,:)
            pvecs(3,:) = vecs(4,:)
            pvecs(4,:) = vecs(3,:)
            pvecs(5,:) = vecs(5,:)
            call init_event(pvecs, -1, -1, +1)
            !---#] reinitialize kinematics:
            fr = finite_renormalisation0(real(scale2,ki))
            heli_amp(0) = heli_amp(0) + fr
         end if
         ok = ok .and. my_ok
         amp = amp + heli_amp
         rat2 = rat2 + rational2

         if(debug_nlo_diagrams) then
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-finite' value='", heli_amp(0), "'/>"
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-single' value='", heli_amp(-1), "'/>"
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-double' value='", heli_amp(-2), "'/>"
            if (corrections_are_qcd .and. renorm_gamma5) then
               write(logfile,'(A30,E24.16,A3)') &
                   & "<result kind='fin-ren' value='", fr, "'/>"
            end if
            if(my_ok) then
               write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
            else
               write(logfile,'(A29)') "<flag name='ok' status='no'/>"
            end if
            write(logfile,*) "</helicity>"
         end if
      end if
      if (eval_heli(5)) then
         if(debug_nlo_diagrams) then
            write(logfile,*) "<helicity index='5'>"
         end if
         !---#[ reinitialize kinematics:
         pvecs(1,:) = vecs(1,:)
         pvecs(2,:) = vecs(2,:)
         pvecs(3,:) = vecs(4,:)
         pvecs(4,:) = vecs(3,:)
         pvecs(5,:) = vecs(5,:)
         call init_event(pvecs, +1, -1, +1)
            !---#] reinitialize kinematics:
         do c=1,numcs
            colorvec(c,:) = samplitudeh0l1(real(scale2,ki),my_ok,rational2,c)
         end do
         heli_amp( 0) = square(colorvec(:, 0))
         heli_amp(-1) = square(colorvec(:,-1))
         heli_amp(-2) = square(colorvec(:,-2))
      
         if (corrections_are_qcd .and. renorm_gamma5) then
            !---#[ reinitialize kinematics:
            pvecs(1,:) = vecs(1,:)
            pvecs(2,:) = vecs(2,:)
            pvecs(3,:) = vecs(4,:)
            pvecs(4,:) = vecs(3,:)
            pvecs(5,:) = vecs(5,:)
            call init_event(pvecs, +1, -1, +1)
            !---#] reinitialize kinematics:
            fr = finite_renormalisation0(real(scale2,ki))
            heli_amp(0) = heli_amp(0) + fr
         end if
         ok = ok .and. my_ok
         amp = amp + heli_amp
         rat2 = rat2 + rational2

         if(debug_nlo_diagrams) then
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-finite' value='", heli_amp(0), "'/>"
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-single' value='", heli_amp(-1), "'/>"
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-double' value='", heli_amp(-2), "'/>"
            if (corrections_are_qcd .and. renorm_gamma5) then
               write(logfile,'(A30,E24.16,A3)') &
                   & "<result kind='fin-ren' value='", fr, "'/>"
            end if
            if(my_ok) then
               write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
            else
               write(logfile,'(A29)') "<flag name='ok' status='no'/>"
            end if
            write(logfile,*) "</helicity>"
         end if
      end if
      if (eval_heli(6)) then
         if(debug_nlo_diagrams) then
            write(logfile,*) "<helicity index='6'>"
         end if
         !---#[ reinitialize kinematics:
         pvecs(1,:) = vecs(1,:)
         pvecs(2,:) = vecs(2,:)
         pvecs(3,:) = vecs(4,:)
         pvecs(4,:) = vecs(3,:)
         pvecs(5,:) = vecs(5,:)
         call init_event(pvecs, -1, +1, +1)
            !---#] reinitialize kinematics:
         do c=1,numcs
            colorvec(c,:) = samplitudeh0l1(real(scale2,ki),my_ok,rational2,c)
         end do
         heli_amp( 0) = square(colorvec(:, 0))
         heli_amp(-1) = square(colorvec(:,-1))
         heli_amp(-2) = square(colorvec(:,-2))
      
         if (corrections_are_qcd .and. renorm_gamma5) then
            !---#[ reinitialize kinematics:
            pvecs(1,:) = vecs(1,:)
            pvecs(2,:) = vecs(2,:)
            pvecs(3,:) = vecs(4,:)
            pvecs(4,:) = vecs(3,:)
            pvecs(5,:) = vecs(5,:)
            call init_event(pvecs, -1, +1, +1)
            !---#] reinitialize kinematics:
            fr = finite_renormalisation0(real(scale2,ki))
            heli_amp(0) = heli_amp(0) + fr
         end if
         ok = ok .and. my_ok
         amp = amp + heli_amp
         rat2 = rat2 + rational2

         if(debug_nlo_diagrams) then
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-finite' value='", heli_amp(0), "'/>"
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-single' value='", heli_amp(-1), "'/>"
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-double' value='", heli_amp(-2), "'/>"
            if (corrections_are_qcd .and. renorm_gamma5) then
               write(logfile,'(A30,E24.16,A3)') &
                   & "<result kind='fin-ren' value='", fr, "'/>"
            end if
            if(my_ok) then
               write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
            else
               write(logfile,'(A29)') "<flag name='ok' status='no'/>"
            end if
            write(logfile,*) "</helicity>"
         end if
      end if
      if (eval_heli(7)) then
         if(debug_nlo_diagrams) then
            write(logfile,*) "<helicity index='7'>"
         end if
         !---#[ reinitialize kinematics:
         pvecs(1,:) = vecs(1,:)
         pvecs(2,:) = vecs(2,:)
         pvecs(3,:) = vecs(4,:)
         pvecs(4,:) = vecs(3,:)
         pvecs(5,:) = vecs(5,:)
         call init_event(pvecs, +1, +1, +1)
            !---#] reinitialize kinematics:
         do c=1,numcs
            colorvec(c,:) = samplitudeh0l1(real(scale2,ki),my_ok,rational2,c)
         end do
         heli_amp( 0) = square(colorvec(:, 0))
         heli_amp(-1) = square(colorvec(:,-1))
         heli_amp(-2) = square(colorvec(:,-2))
      
         if (corrections_are_qcd .and. renorm_gamma5) then
            !---#[ reinitialize kinematics:
            pvecs(1,:) = vecs(1,:)
            pvecs(2,:) = vecs(2,:)
            pvecs(3,:) = vecs(4,:)
            pvecs(4,:) = vecs(3,:)
            pvecs(5,:) = vecs(5,:)
            call init_event(pvecs, +1, +1, +1)
            !---#] reinitialize kinematics:
            fr = finite_renormalisation0(real(scale2,ki))
            heli_amp(0) = heli_amp(0) + fr
         end if
         ok = ok .and. my_ok
         amp = amp + heli_amp
         rat2 = rat2 + rational2

         if(debug_nlo_diagrams) then
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-finite' value='", heli_amp(0), "'/>"
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-single' value='", heli_amp(-1), "'/>"
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-double' value='", heli_amp(-2), "'/>"
            if (corrections_are_qcd .and. renorm_gamma5) then
               write(logfile,'(A30,E24.16,A3)') &
                   & "<result kind='fin-ren' value='", fr, "'/>"
            end if
            if(my_ok) then
               write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
            else
               write(logfile,'(A29)') "<flag name='ok' status='no'/>"
            end if
            write(logfile,*) "</helicity>"
         end if
      end if
      if (include_helicity_avg_factor) then
         amp = amp / real(in_helicities, ki)
      end if
      if (include_color_avg_factor) then
         amp = amp / incolors
      end if
      if (include_symmetry_factor) then
         amp = amp / real(symmetry_factor, ki)
      end if
   end function samplitudel1
   !---#] function samplitudel1 :
   !---#[ subroutine ir_subtraction :
   subroutine     ir_subtraction(vecs,scale2,amp,h)
      use p2_part21part21_part25part25part21_config, only: &
         & nlo_prefactors
      use p2_part21part21_part25part25part21_dipoles, only: pi
      use p2_part21part21_part25part25part21_kinematics, only: &
         & init_event, corrections_are_qcd
      use p2_part21part21_part25part25part21_model
      implicit none
      real(ki), dimension(5, 4), intent(in) :: vecs
      real(ki), intent(in) :: scale2
      integer, optional, intent(in) :: h
      real(ki), dimension(2), intent(out) :: amp
      real(ki), dimension(2) :: heli_amp
      real(ki), dimension(5, 4) :: pvecs
      complex(ki), dimension(numcs,numcs,2) :: oper
      complex(ki), dimension(numcs) :: color_vectorl0, pcolor
      logical, dimension(0:7) :: eval_heli
      real(ki) :: nlo_coupling

      if (present(h)) then
         eval_heli(:) = .false.
         eval_heli(h) = .true.
      else
         eval_heli(:) = .true.
      end if

      if(corrections_are_qcd) then
         nlo_coupling = mdlG*mdlG
      else
         nlo_coupling = mdlee*mdlee
      end if

      if (corrections_are_qcd) then
        oper = insertion_operator(real(scale2,ki), vecs)
      else
        oper = insertion_operator_qed(real(scale2,ki), vecs)
      endif
      amp(:) = 0.0_ki
      select case(nlo_prefactors)
      case(0)
         ! The result is already in its desired form
      case(1)
         amp(:) = amp(:) * nlo_coupling
      case(2)
         amp(:) = amp(:) * nlo_coupling / 8.0_ki / pi / pi
      end select
   end subroutine ir_subtraction
   !---#] subroutine ir_subtraction :
   !---#[ subroutine samplitudel01_qp :
   subroutine     samplitudel01_qp(vecs, scale2, amp, rat2, ok, h)
      use p2_part21part21_part25part25part21_config, only: &
         & debug_lo_diagrams, debug_nlo_diagrams, logfile, deltaOS, &
         & renormalisation, renorm_beta, renorm_mqwf, renorm_decoupling, &
         & renorm_logs, renorm_mqse, renorm_yukawa, nlo_prefactors
      use p2_part21part21_part25part25part21_kinematics_qp, only: &
         & inspect_kinematics, init_event
      use p2_part21part21_part25part25part21_model_qp
      use p2_part21part21_part25part25part21_dipoles_qp, only: pi
      implicit none
      real(ki_qp), dimension(5, 4), intent(in) :: vecs
      real(ki_qp), intent(in) :: scale2
      real(ki_qp), dimension(4), intent(out) :: amp
      real(ki_qp), intent(out) :: rat2
      logical, intent(out), optional :: ok
      integer, intent(in), optional :: h
      real(ki_qp) :: nlo_coupling

      complex(ki_qp), parameter :: i_ = (0.0_ki_qp, 1.0_ki_qp)

      ! Number of heavy quark flavours in loops.
      real(ki_qp), parameter :: NFh_qp = 0.0_ki_qp

      logical :: my_ok

      ! used for m=0 QCD renormalisation
      real(ki_qp) :: beta0

      if(corrections_are_qcd) then
         nlo_coupling = mdlG*mdlG
      else
         nlo_coupling = mdlee*mdlee
      end if

      if(debug_lo_diagrams .or. debug_nlo_diagrams) then
         call init_event(vecs)
         write(logfile,'(A7)') "<event>"
         call inspect_kinematics(logfile)
      end if

      
      amp(1)   = 0.0_ki_qp
      select case (renormalisation)
      case (0)
         ! no renormalisation
         deltaOS = 0.0_ki_qp
      case (1)
         ! fully renormalized
         if(renorm_mqse) then
            deltaOS = 1.0_ki_qp
         else
            deltaOS = 0.0_ki_qp
         end if
      case (2)
         ! massive quark counterterms only
         deltaOS = 1.0_ki_qp
      case default
         ! not implemented
         print*, "In p2_part21part21_part25part25part21_matrix:"
         print*, "  invalid value for renormalisation=", renormalisation
         stop
      end select

      if (present(h)) then
         amp((/4,3,2/)) = samplitudel1_qp(vecs, scale2, my_ok, rat2, h)/nlo_coupling/nlo_coupling
      else
         amp((/4,3,2/)) = samplitudel1_qp(vecs, scale2, my_ok, rat2)/nlo_coupling/nlo_coupling
      end if
      select case (renormalisation)
      case (0)
         ! no renormalisation
      case (1)
         ! fully renormalized
         ! No tree level present
      case (2)
         ! massive quark counterterms only
      case default
         ! not implemented
         print*, "In p2_part21part21_part25part25part21_matrix:"
         print*, "  invalid value for renormalisation=", renormalisation
         stop
      end select
      if (convert_to_cdr) then
         ! Scheme conversion for infrared structure
         ! Reference:
         ! S. Catani, M. H. Seymour, Z. Trocsanyi,
         ! ``Regularisation scheme independence and unitarity
         !   in QCD cross-sections,''
         ! Phys.Rev. D 55 (1997) 6819
         ! arXiv:hep-ph/9610553
         amp(2) = amp(2) - amp(1) * (&
           &          num_light_quarks * 0.5_ki_qp * CF_qp &
           &        + num_gluons * 1.0_ki_qp/6.0_ki_qp * CA_qp)
      end if
      if (present(ok)) ok = my_ok

      if(debug_lo_diagrams .or. debug_nlo_diagrams) then
         write(logfile,'(A25,E24.16,A3)') &
            & "<result kind='lo' value='", amp(1), "'/>"
         write(logfile,'(A33,E24.16,A3)') &
            & "<result kind='nlo-finite' value='", amp(2), "'/>"
         write(logfile,'(A33,E24.16,A3)') &
            & "<result kind='nlo-single' value='", amp(3), "'/>"
         write(logfile,'(A33,E24.16,A3)') &
            & "<result kind='nlo-double' value='", amp(4), "'/>"
         if(my_ok) then
            write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
         else
            write(logfile,'(A29)') "<flag name='ok' status='no'/>"
         end if
         write(logfile,'(A8)') "</event>"
      end if
      select case(nlo_prefactors)
      case(0)
         ! The result is already in its desired form
      case(1)
         ! loop-induced
         amp(2:4) = amp(2:4) * nlo_coupling * nlo_coupling
      case(2)
         ! loop-induced
         amp(2:4) = amp(2:4) * (nlo_coupling / 8.0_ki_qp / pi / pi)**2
      end select
   end subroutine samplitudel01_qp
   !---#] subroutine samplitudel01_qp :
   !---#[ function samplitudel0_qp :
   function     samplitudel0_qp(vecs, h) result(amp)
      use p2_part21part21_part25part25part21_config, only: logfile
      use p2_part21part21_part25part25part21_kinematics_qp, only: init_event
      implicit none
      real(ki_qp), dimension(5, 4), intent(in) :: vecs
      integer, optional, intent(in) :: h
      real(ki_qp) :: amp, heli_amp
      complex(ki_qp), dimension(numcs) :: color_vector
      logical, dimension(0:7) :: eval_heli
      real(ki_qp), dimension(5, 4) :: pvecs

      if (present(h)) then
         eval_heli(:) = .false.
         eval_heli(h) = .true.
      else
         eval_heli(:) = .true.
      end if

      amp = 0.0_ki_qp
   end function samplitudel0_qp
   !---#] function samplitudel0_qp :
   !---#[ function samplitudel1_qp :
   function     samplitudel1_qp(vecs,scale2,ok,rat2,h) result(amp)
      use p2_part21part21_part25part25part21_config, only: &
         & debug_nlo_diagrams, logfile, renorm_gamma5
      use p2_part21part21_part25part25part21_kinematics_qp, only: init_event
      implicit none
      real(ki_qp), dimension(5, 4), intent(in) :: vecs
      logical, intent(out) :: ok
      real(ki_qp), intent(in) :: scale2
      real(ki_qp), intent(out) :: rat2
      integer, optional, intent(in) :: h
      real(ki_qp), dimension(5, 4) :: pvecs
      real(ki_qp), dimension(-2:0) :: amp, heli_amp
      complex(ki_qp), dimension(numcs,-2:0) :: colorvec
      integer :: c
      logical :: my_ok
      logical, dimension(0:7) :: eval_heli
      real(ki_qp) :: fr, rational2

      amp(:) = 0.0_ki_qp
      rat2 = 0.0_ki_qp
      ok = .true.

      if (present(h)) then
         eval_heli(:) = .false.
         eval_heli(h) = .true.
      else
         eval_heli(:) = .true.
      end if


      if (eval_heli(0)) then
         if(debug_nlo_diagrams) then
            write(logfile,*) "<helicity index='0'>"
         end if
         !---#[ reinitialize kinematics:
         pvecs(1,:) = vecs(1,:)
         pvecs(2,:) = vecs(2,:)
         pvecs(3,:) = vecs(3,:)
         pvecs(4,:) = vecs(4,:)
         pvecs(5,:) = vecs(5,:)
         call init_event(pvecs, -1, -1, -1)
            !---#] reinitialize kinematics:
         do c=1,numcs
            colorvec(c,:) = samplitudeh0l1_qp(real(scale2,ki_qp),my_ok,rational2,c)
         end do
         heli_amp( 0) = square_qp(colorvec(:, 0))
         heli_amp(-1) = square_qp(colorvec(:,-1))
         heli_amp(-2) = square_qp(colorvec(:,-2))
      
         if (corrections_are_qcd .and. renorm_gamma5) then
            !---#[ reinitialize kinematics:
            pvecs(1,:) = vecs(1,:)
            pvecs(2,:) = vecs(2,:)
            pvecs(3,:) = vecs(3,:)
            pvecs(4,:) = vecs(4,:)
            pvecs(5,:) = vecs(5,:)
            call init_event(pvecs, -1, -1, -1)
            !---#] reinitialize kinematics:
            fr = finite_renormalisation0_qp(real(scale2,ki_qp))
            heli_amp(0) = heli_amp(0) + fr
         end if
         ok = ok .and. my_ok
         amp = amp + heli_amp
         rat2 = rat2 + rational2

         if(debug_nlo_diagrams) then
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-finite' value='", heli_amp(0), "'/>"
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-single' value='", heli_amp(-1), "'/>"
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-double' value='", heli_amp(-2), "'/>"
            if (corrections_are_qcd .and. renorm_gamma5) then
               write(logfile,'(A30,E24.16,A3)') &
                   & "<result kind='fin-ren' value='", fr, "'/>"
            end if
            if(my_ok) then
               write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
            else
               write(logfile,'(A29)') "<flag name='ok' status='no'/>"
            end if
            write(logfile,*) "</helicity>"
         end if
      end if
      if (eval_heli(1)) then
         if(debug_nlo_diagrams) then
            write(logfile,*) "<helicity index='1'>"
         end if
         !---#[ reinitialize kinematics:
         pvecs(1,:) = vecs(1,:)
         pvecs(2,:) = vecs(2,:)
         pvecs(3,:) = vecs(4,:)
         pvecs(4,:) = vecs(3,:)
         pvecs(5,:) = vecs(5,:)
         call init_event(pvecs, +1, -1, -1)
            !---#] reinitialize kinematics:
         do c=1,numcs
            colorvec(c,:) = samplitudeh0l1_qp(real(scale2,ki_qp),my_ok,rational2,c)
         end do
         heli_amp( 0) = square_qp(colorvec(:, 0))
         heli_amp(-1) = square_qp(colorvec(:,-1))
         heli_amp(-2) = square_qp(colorvec(:,-2))
      
         if (corrections_are_qcd .and. renorm_gamma5) then
            !---#[ reinitialize kinematics:
            pvecs(1,:) = vecs(1,:)
            pvecs(2,:) = vecs(2,:)
            pvecs(3,:) = vecs(4,:)
            pvecs(4,:) = vecs(3,:)
            pvecs(5,:) = vecs(5,:)
            call init_event(pvecs, +1, -1, -1)
            !---#] reinitialize kinematics:
            fr = finite_renormalisation0_qp(real(scale2,ki_qp))
            heli_amp(0) = heli_amp(0) + fr
         end if
         ok = ok .and. my_ok
         amp = amp + heli_amp
         rat2 = rat2 + rational2

         if(debug_nlo_diagrams) then
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-finite' value='", heli_amp(0), "'/>"
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-single' value='", heli_amp(-1), "'/>"
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-double' value='", heli_amp(-2), "'/>"
            if (corrections_are_qcd .and. renorm_gamma5) then
               write(logfile,'(A30,E24.16,A3)') &
                   & "<result kind='fin-ren' value='", fr, "'/>"
            end if
            if(my_ok) then
               write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
            else
               write(logfile,'(A29)') "<flag name='ok' status='no'/>"
            end if
            write(logfile,*) "</helicity>"
         end if
      end if
      if (eval_heli(2)) then
         if(debug_nlo_diagrams) then
            write(logfile,*) "<helicity index='2'>"
         end if
         !---#[ reinitialize kinematics:
         pvecs(1,:) = vecs(1,:)
         pvecs(2,:) = vecs(2,:)
         pvecs(3,:) = vecs(4,:)
         pvecs(4,:) = vecs(3,:)
         pvecs(5,:) = vecs(5,:)
         call init_event(pvecs, -1, +1, -1)
            !---#] reinitialize kinematics:
         do c=1,numcs
            colorvec(c,:) = samplitudeh0l1_qp(real(scale2,ki_qp),my_ok,rational2,c)
         end do
         heli_amp( 0) = square_qp(colorvec(:, 0))
         heli_amp(-1) = square_qp(colorvec(:,-1))
         heli_amp(-2) = square_qp(colorvec(:,-2))
      
         if (corrections_are_qcd .and. renorm_gamma5) then
            !---#[ reinitialize kinematics:
            pvecs(1,:) = vecs(1,:)
            pvecs(2,:) = vecs(2,:)
            pvecs(3,:) = vecs(4,:)
            pvecs(4,:) = vecs(3,:)
            pvecs(5,:) = vecs(5,:)
            call init_event(pvecs, -1, +1, -1)
            !---#] reinitialize kinematics:
            fr = finite_renormalisation0_qp(real(scale2,ki_qp))
            heli_amp(0) = heli_amp(0) + fr
         end if
         ok = ok .and. my_ok
         amp = amp + heli_amp
         rat2 = rat2 + rational2

         if(debug_nlo_diagrams) then
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-finite' value='", heli_amp(0), "'/>"
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-single' value='", heli_amp(-1), "'/>"
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-double' value='", heli_amp(-2), "'/>"
            if (corrections_are_qcd .and. renorm_gamma5) then
               write(logfile,'(A30,E24.16,A3)') &
                   & "<result kind='fin-ren' value='", fr, "'/>"
            end if
            if(my_ok) then
               write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
            else
               write(logfile,'(A29)') "<flag name='ok' status='no'/>"
            end if
            write(logfile,*) "</helicity>"
         end if
      end if
      if (eval_heli(3)) then
         if(debug_nlo_diagrams) then
            write(logfile,*) "<helicity index='3'>"
         end if
         !---#[ reinitialize kinematics:
         pvecs(1,:) = vecs(1,:)
         pvecs(2,:) = vecs(2,:)
         pvecs(3,:) = vecs(4,:)
         pvecs(4,:) = vecs(3,:)
         pvecs(5,:) = vecs(5,:)
         call init_event(pvecs, +1, +1, -1)
            !---#] reinitialize kinematics:
         do c=1,numcs
            colorvec(c,:) = samplitudeh0l1_qp(real(scale2,ki_qp),my_ok,rational2,c)
         end do
         heli_amp( 0) = square_qp(colorvec(:, 0))
         heli_amp(-1) = square_qp(colorvec(:,-1))
         heli_amp(-2) = square_qp(colorvec(:,-2))
      
         if (corrections_are_qcd .and. renorm_gamma5) then
            !---#[ reinitialize kinematics:
            pvecs(1,:) = vecs(1,:)
            pvecs(2,:) = vecs(2,:)
            pvecs(3,:) = vecs(4,:)
            pvecs(4,:) = vecs(3,:)
            pvecs(5,:) = vecs(5,:)
            call init_event(pvecs, +1, +1, -1)
            !---#] reinitialize kinematics:
            fr = finite_renormalisation0_qp(real(scale2,ki_qp))
            heli_amp(0) = heli_amp(0) + fr
         end if
         ok = ok .and. my_ok
         amp = amp + heli_amp
         rat2 = rat2 + rational2

         if(debug_nlo_diagrams) then
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-finite' value='", heli_amp(0), "'/>"
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-single' value='", heli_amp(-1), "'/>"
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-double' value='", heli_amp(-2), "'/>"
            if (corrections_are_qcd .and. renorm_gamma5) then
               write(logfile,'(A30,E24.16,A3)') &
                   & "<result kind='fin-ren' value='", fr, "'/>"
            end if
            if(my_ok) then
               write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
            else
               write(logfile,'(A29)') "<flag name='ok' status='no'/>"
            end if
            write(logfile,*) "</helicity>"
         end if
      end if
      if (eval_heli(4)) then
         if(debug_nlo_diagrams) then
            write(logfile,*) "<helicity index='4'>"
         end if
         !---#[ reinitialize kinematics:
         pvecs(1,:) = vecs(1,:)
         pvecs(2,:) = vecs(2,:)
         pvecs(3,:) = vecs(4,:)
         pvecs(4,:) = vecs(3,:)
         pvecs(5,:) = vecs(5,:)
         call init_event(pvecs, -1, -1, +1)
            !---#] reinitialize kinematics:
         do c=1,numcs
            colorvec(c,:) = samplitudeh0l1_qp(real(scale2,ki_qp),my_ok,rational2,c)
         end do
         heli_amp( 0) = square_qp(colorvec(:, 0))
         heli_amp(-1) = square_qp(colorvec(:,-1))
         heli_amp(-2) = square_qp(colorvec(:,-2))
      
         if (corrections_are_qcd .and. renorm_gamma5) then
            !---#[ reinitialize kinematics:
            pvecs(1,:) = vecs(1,:)
            pvecs(2,:) = vecs(2,:)
            pvecs(3,:) = vecs(4,:)
            pvecs(4,:) = vecs(3,:)
            pvecs(5,:) = vecs(5,:)
            call init_event(pvecs, -1, -1, +1)
            !---#] reinitialize kinematics:
            fr = finite_renormalisation0_qp(real(scale2,ki_qp))
            heli_amp(0) = heli_amp(0) + fr
         end if
         ok = ok .and. my_ok
         amp = amp + heli_amp
         rat2 = rat2 + rational2

         if(debug_nlo_diagrams) then
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-finite' value='", heli_amp(0), "'/>"
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-single' value='", heli_amp(-1), "'/>"
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-double' value='", heli_amp(-2), "'/>"
            if (corrections_are_qcd .and. renorm_gamma5) then
               write(logfile,'(A30,E24.16,A3)') &
                   & "<result kind='fin-ren' value='", fr, "'/>"
            end if
            if(my_ok) then
               write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
            else
               write(logfile,'(A29)') "<flag name='ok' status='no'/>"
            end if
            write(logfile,*) "</helicity>"
         end if
      end if
      if (eval_heli(5)) then
         if(debug_nlo_diagrams) then
            write(logfile,*) "<helicity index='5'>"
         end if
         !---#[ reinitialize kinematics:
         pvecs(1,:) = vecs(1,:)
         pvecs(2,:) = vecs(2,:)
         pvecs(3,:) = vecs(4,:)
         pvecs(4,:) = vecs(3,:)
         pvecs(5,:) = vecs(5,:)
         call init_event(pvecs, +1, -1, +1)
            !---#] reinitialize kinematics:
         do c=1,numcs
            colorvec(c,:) = samplitudeh0l1_qp(real(scale2,ki_qp),my_ok,rational2,c)
         end do
         heli_amp( 0) = square_qp(colorvec(:, 0))
         heli_amp(-1) = square_qp(colorvec(:,-1))
         heli_amp(-2) = square_qp(colorvec(:,-2))
      
         if (corrections_are_qcd .and. renorm_gamma5) then
            !---#[ reinitialize kinematics:
            pvecs(1,:) = vecs(1,:)
            pvecs(2,:) = vecs(2,:)
            pvecs(3,:) = vecs(4,:)
            pvecs(4,:) = vecs(3,:)
            pvecs(5,:) = vecs(5,:)
            call init_event(pvecs, +1, -1, +1)
            !---#] reinitialize kinematics:
            fr = finite_renormalisation0_qp(real(scale2,ki_qp))
            heli_amp(0) = heli_amp(0) + fr
         end if
         ok = ok .and. my_ok
         amp = amp + heli_amp
         rat2 = rat2 + rational2

         if(debug_nlo_diagrams) then
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-finite' value='", heli_amp(0), "'/>"
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-single' value='", heli_amp(-1), "'/>"
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-double' value='", heli_amp(-2), "'/>"
            if (corrections_are_qcd .and. renorm_gamma5) then
               write(logfile,'(A30,E24.16,A3)') &
                   & "<result kind='fin-ren' value='", fr, "'/>"
            end if
            if(my_ok) then
               write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
            else
               write(logfile,'(A29)') "<flag name='ok' status='no'/>"
            end if
            write(logfile,*) "</helicity>"
         end if
      end if
      if (eval_heli(6)) then
         if(debug_nlo_diagrams) then
            write(logfile,*) "<helicity index='6'>"
         end if
         !---#[ reinitialize kinematics:
         pvecs(1,:) = vecs(1,:)
         pvecs(2,:) = vecs(2,:)
         pvecs(3,:) = vecs(4,:)
         pvecs(4,:) = vecs(3,:)
         pvecs(5,:) = vecs(5,:)
         call init_event(pvecs, -1, +1, +1)
            !---#] reinitialize kinematics:
         do c=1,numcs
            colorvec(c,:) = samplitudeh0l1_qp(real(scale2,ki_qp),my_ok,rational2,c)
         end do
         heli_amp( 0) = square_qp(colorvec(:, 0))
         heli_amp(-1) = square_qp(colorvec(:,-1))
         heli_amp(-2) = square_qp(colorvec(:,-2))
      
         if (corrections_are_qcd .and. renorm_gamma5) then
            !---#[ reinitialize kinematics:
            pvecs(1,:) = vecs(1,:)
            pvecs(2,:) = vecs(2,:)
            pvecs(3,:) = vecs(4,:)
            pvecs(4,:) = vecs(3,:)
            pvecs(5,:) = vecs(5,:)
            call init_event(pvecs, -1, +1, +1)
            !---#] reinitialize kinematics:
            fr = finite_renormalisation0_qp(real(scale2,ki_qp))
            heli_amp(0) = heli_amp(0) + fr
         end if
         ok = ok .and. my_ok
         amp = amp + heli_amp
         rat2 = rat2 + rational2

         if(debug_nlo_diagrams) then
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-finite' value='", heli_amp(0), "'/>"
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-single' value='", heli_amp(-1), "'/>"
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-double' value='", heli_amp(-2), "'/>"
            if (corrections_are_qcd .and. renorm_gamma5) then
               write(logfile,'(A30,E24.16,A3)') &
                   & "<result kind='fin-ren' value='", fr, "'/>"
            end if
            if(my_ok) then
               write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
            else
               write(logfile,'(A29)') "<flag name='ok' status='no'/>"
            end if
            write(logfile,*) "</helicity>"
         end if
      end if
      if (eval_heli(7)) then
         if(debug_nlo_diagrams) then
            write(logfile,*) "<helicity index='7'>"
         end if
         !---#[ reinitialize kinematics:
         pvecs(1,:) = vecs(1,:)
         pvecs(2,:) = vecs(2,:)
         pvecs(3,:) = vecs(4,:)
         pvecs(4,:) = vecs(3,:)
         pvecs(5,:) = vecs(5,:)
         call init_event(pvecs, +1, +1, +1)
            !---#] reinitialize kinematics:
         do c=1,numcs
            colorvec(c,:) = samplitudeh0l1_qp(real(scale2,ki_qp),my_ok,rational2,c)
         end do
         heli_amp( 0) = square_qp(colorvec(:, 0))
         heli_amp(-1) = square_qp(colorvec(:,-1))
         heli_amp(-2) = square_qp(colorvec(:,-2))
      
         if (corrections_are_qcd .and. renorm_gamma5) then
            !---#[ reinitialize kinematics:
            pvecs(1,:) = vecs(1,:)
            pvecs(2,:) = vecs(2,:)
            pvecs(3,:) = vecs(4,:)
            pvecs(4,:) = vecs(3,:)
            pvecs(5,:) = vecs(5,:)
            call init_event(pvecs, +1, +1, +1)
            !---#] reinitialize kinematics:
            fr = finite_renormalisation0_qp(real(scale2,ki_qp))
            heli_amp(0) = heli_amp(0) + fr
         end if
         ok = ok .and. my_ok
         amp = amp + heli_amp
         rat2 = rat2 + rational2

         if(debug_nlo_diagrams) then
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-finite' value='", heli_amp(0), "'/>"
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-single' value='", heli_amp(-1), "'/>"
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-double' value='", heli_amp(-2), "'/>"
            if (corrections_are_qcd .and. renorm_gamma5) then
               write(logfile,'(A30,E24.16,A3)') &
                   & "<result kind='fin-ren' value='", fr, "'/>"
            end if
            if(my_ok) then
               write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
            else
               write(logfile,'(A29)') "<flag name='ok' status='no'/>"
            end if
            write(logfile,*) "</helicity>"
         end if
      end if
      if (include_helicity_avg_factor) then
         amp = amp / real(in_helicities, ki_qp)
      end if
      if (include_color_avg_factor) then
         amp = amp / incolors
      end if
      if (include_symmetry_factor) then
         amp = amp / real(symmetry_factor, ki_qp)
      end if
   end function samplitudel1_qp
   !---#] function samplitudel1_qp :
   !---#[ subroutine ir_subtraction_qp :
   subroutine     ir_subtraction_qp(vecs,scale2,amp,h)
      use p2_part21part21_part25part25part21_config, only: &
         & nlo_prefactors
      use p2_part21part21_part25part25part21_dipoles_qp, only: pi
      use p2_part21part21_part25part25part21_kinematics_qp, only: &
         & init_event, corrections_are_qcd
      use p2_part21part21_part25part25part21_model_qp
      implicit none
      real(ki_qp), dimension(5, 4), intent(in) :: vecs
      real(ki_qp), intent(in) :: scale2
      integer, optional, intent(in) :: h
      real(ki_qp), dimension(2), intent(out) :: amp
      real(ki_qp), dimension(2) :: heli_amp
      real(ki_qp), dimension(5, 4) :: pvecs
      complex(ki_qp), dimension(numcs,numcs,2) :: oper
      complex(ki_qp), dimension(numcs) :: color_vectorl0, pcolor
      logical, dimension(0:7) :: eval_heli
      real(ki_qp) :: nlo_coupling

      if (present(h)) then
         eval_heli(:) = .false.
         eval_heli(h) = .true.
      else
         eval_heli(:) = .true.
      end if

      if(corrections_are_qcd) then
         nlo_coupling = mdlG*mdlG
      else
         nlo_coupling = mdlee*mdlee
      end if

      if (corrections_are_qcd) then
        oper = insertion_operator_qp(real(scale2,ki_qp), vecs)
      else
        oper = insertion_operator_qed_qp(real(scale2,ki_qp), vecs)
      endif
      amp(:) = 0.0_ki_qp
      select case(nlo_prefactors)
      case(0)
         ! The result is already in its desired form
      case(1)
         amp(:) = amp(:) * nlo_coupling
      case(2)
         amp(:) = amp(:) * nlo_coupling / 8.0_ki_qp / pi / pi
      end select
   end subroutine ir_subtraction_qp
   !---#] subroutine ir_subtraction_qp :
   !---#[ color correlated ME :
   pure subroutine color_correlated_lo(color_vector,res)
      use p2_part21part21_part25part25part21_color, only: T1T1, &
      & T1T2, &
      & T1T5, &
      & T2T2, &
      & T2T5, &
      & T5T5
      implicit none
      complex(ki), dimension(numcs), intent(in) :: color_vector
      real(ki), dimension(num_legs,num_legs), intent(out) :: res
      res(:,:)=0.0_ki
      res(1,1) = square(color_vector,T1T1)
      res(1,1) = res(1,1)
      res(1,2) = square(color_vector,T1T2)
      res(2,1) = res(1,2)
      res(1,5) = square(color_vector,T1T5)
      res(5,1) = res(1,5)
      res(2,2) = square(color_vector,T2T2)
      res(2,2) = res(2,2)
      res(2,5) = square(color_vector,T2T5)
      res(5,2) = res(2,5)
      res(5,5) = square(color_vector,T5T5)
      res(5,5) = res(5,5)
   end subroutine color_correlated_lo

   subroutine     color_correlated_lo2(vecs,borncc)
      use p2_part21part21_part25part25part21_kinematics, only: init_event
      implicit none
      real(ki), dimension(num_legs, 4), intent(in) :: vecs
      real(ki), dimension(num_legs,num_legs), intent(out) :: borncc
      real(ki), dimension(num_legs,num_legs) :: borncc_heli
      real(ki), dimension(num_legs, 4) :: pvecs
      complex(ki), dimension(numcs) :: color_vector

      borncc(:,:) = 0.0_ki
   end subroutine color_correlated_lo2


   pure subroutine OLP_color_correlated_lo(color_vector,res)
      use p2_part21part21_part25part25part21_color, only: T1T1, &
      & T1T2, &
      & T1T5, &
      & T2T2, &
      & T2T5, &
      & T5T5
      implicit none
      complex(ki), dimension(numcs), intent(in) :: color_vector
      real(ki), dimension(num_legs*(num_legs-1)/2), intent(out) :: res
      res(:)=0.0_ki   
      res(1) = square(color_vector,T1T2)  
      res(7) = square(color_vector,T1T5)    
      res(8) = square(color_vector,T2T5)   
   end subroutine OLP_color_correlated_lo


   subroutine OLP_color_correlated(vecs,ampcc)
      use p2_part21part21_part25part25part21_kinematics, only: init_event
      implicit none
      real(ki), dimension(num_legs, 4), intent(in) :: vecs
      real(ki), dimension(num_legs*(num_legs-1)/2), intent(out) :: ampcc
      real(ki), dimension(num_legs,num_legs) :: borncc
      real(ki), dimension(num_legs*(num_legs-1)/2) :: ampcc_heli
      real(ki), dimension(num_legs, 4) :: pvecs
      complex(ki), dimension(numcs) :: color_vector
      ampcc(:) = 0.0_ki
      !---#[ reinitialize kinematics:
      pvecs(1,:) = vecs(1,:)
      pvecs(2,:) = vecs(2,:)
      pvecs(3,:) = vecs(3,:)
      pvecs(4,:) = vecs(4,:)
      pvecs(5,:) = vecs(5,:)
      call init_event(pvecs, -1, -1, -1)
      !---#] reinitialize kinematics:
!      color_vector = amplitude0l0()
      color_vector = 0
      call OLP_color_correlated_lo(color_vector,ampcc_heli)
      ampcc(:) = ampcc(:) + ampcc_heli(:)
      !---#[ reinitialize kinematics:
      pvecs(1,:) = vecs(1,:)
      pvecs(2,:) = vecs(2,:)
      pvecs(3,:) = vecs(4,:)
      pvecs(4,:) = vecs(3,:)
      pvecs(5,:) = vecs(5,:)
      call init_event(pvecs, +1, -1, -1)
      !---#] reinitialize kinematics:
!      color_vector = amplitude0l0()
      color_vector = 0
      call OLP_color_correlated_lo(color_vector,ampcc_heli)
      ampcc(:) = ampcc(:) + ampcc_heli(:)
      !---#[ reinitialize kinematics:
      pvecs(1,:) = vecs(1,:)
      pvecs(2,:) = vecs(2,:)
      pvecs(3,:) = vecs(4,:)
      pvecs(4,:) = vecs(3,:)
      pvecs(5,:) = vecs(5,:)
      call init_event(pvecs, -1, +1, -1)
      !---#] reinitialize kinematics:
!      color_vector = amplitude0l0()
      color_vector = 0
      call OLP_color_correlated_lo(color_vector,ampcc_heli)
      ampcc(:) = ampcc(:) + ampcc_heli(:)
      !---#[ reinitialize kinematics:
      pvecs(1,:) = vecs(1,:)
      pvecs(2,:) = vecs(2,:)
      pvecs(3,:) = vecs(4,:)
      pvecs(4,:) = vecs(3,:)
      pvecs(5,:) = vecs(5,:)
      call init_event(pvecs, +1, +1, -1)
      !---#] reinitialize kinematics:
!      color_vector = amplitude0l0()
      color_vector = 0
      call OLP_color_correlated_lo(color_vector,ampcc_heli)
      ampcc(:) = ampcc(:) + ampcc_heli(:)
      !---#[ reinitialize kinematics:
      pvecs(1,:) = vecs(1,:)
      pvecs(2,:) = vecs(2,:)
      pvecs(3,:) = vecs(4,:)
      pvecs(4,:) = vecs(3,:)
      pvecs(5,:) = vecs(5,:)
      call init_event(pvecs, -1, -1, +1)
      !---#] reinitialize kinematics:
!      color_vector = amplitude0l0()
      color_vector = 0
      call OLP_color_correlated_lo(color_vector,ampcc_heli)
      ampcc(:) = ampcc(:) + ampcc_heli(:)
      !---#[ reinitialize kinematics:
      pvecs(1,:) = vecs(1,:)
      pvecs(2,:) = vecs(2,:)
      pvecs(3,:) = vecs(4,:)
      pvecs(4,:) = vecs(3,:)
      pvecs(5,:) = vecs(5,:)
      call init_event(pvecs, +1, -1, +1)
      !---#] reinitialize kinematics:
!      color_vector = amplitude0l0()
      color_vector = 0
      call OLP_color_correlated_lo(color_vector,ampcc_heli)
      ampcc(:) = ampcc(:) + ampcc_heli(:)
      !---#[ reinitialize kinematics:
      pvecs(1,:) = vecs(1,:)
      pvecs(2,:) = vecs(2,:)
      pvecs(3,:) = vecs(4,:)
      pvecs(4,:) = vecs(3,:)
      pvecs(5,:) = vecs(5,:)
      call init_event(pvecs, -1, +1, +1)
      !---#] reinitialize kinematics:
!      color_vector = amplitude0l0()
      color_vector = 0
      call OLP_color_correlated_lo(color_vector,ampcc_heli)
      ampcc(:) = ampcc(:) + ampcc_heli(:)
      !---#[ reinitialize kinematics:
      pvecs(1,:) = vecs(1,:)
      pvecs(2,:) = vecs(2,:)
      pvecs(3,:) = vecs(4,:)
      pvecs(4,:) = vecs(3,:)
      pvecs(5,:) = vecs(5,:)
      call init_event(pvecs, +1, +1, +1)
      !---#] reinitialize kinematics:
!      color_vector = amplitude0l0()
      color_vector = 0
      call OLP_color_correlated_lo(color_vector,ampcc_heli)
      ampcc(:) = ampcc(:) + ampcc_heli(:)

      
      if (include_helicity_avg_factor) then
         ampcc = ampcc / real(in_helicities, ki)
      end if
      if (include_color_avg_factor) then
         ampcc = ampcc / incolors
      end if
      if (include_symmetry_factor) then
         ampcc = ampcc / real(symmetry_factor, ki)
      end if


   end subroutine OLP_color_correlated


   !---#] color correlated ME :
   !---#[ spin correlated ME :
   subroutine spin_correlated_lo2(vecs, bornsc)
      use p2_part21part21_part25part25part21_kinematics
      implicit none
      real(ki), dimension(num_legs, 4), intent(in) :: vecs
      real(ki), dimension(num_legs,4,4) :: bornsc
      real(ki), dimension(num_legs, 4) :: pvecs
      complex(ki), dimension(4,4) :: tens
      complex(ki) :: pp, pm, mp, mm

      bornsc(:,:,:) = 0.0_ki
      !---#[ Initialize helicity amplitudes :
   end subroutine spin_correlated_lo2




   subroutine OLP_spin_correlated_lo2(vecs, ampsc)
      use p2_part21part21_part25part25part21_kinematics
      implicit none
      real(ki), dimension(num_legs, 4), intent(in) :: vecs
      real(ki), dimension(2*num_legs*num_legs) :: ampsc
      real(ki), dimension(num_legs, 4) :: pvecs
      integer :: i
      complex(ki) :: pm, mp
      complex(ki), dimension(numcs) :: heli_amp0
      complex(ki), dimension(numcs) :: heli_amp1
      complex(ki), dimension(numcs) :: heli_amp2
      complex(ki), dimension(numcs) :: heli_amp3
      complex(ki), dimension(numcs) :: heli_amp4
      complex(ki), dimension(numcs) :: heli_amp5
      complex(ki), dimension(numcs) :: heli_amp6
      complex(ki), dimension(numcs) :: heli_amp7
      complex(ki), dimension(4) :: eps1
      complex(ki), dimension(4) :: eps2
      complex(ki), dimension(4) :: eps5
      complex(ki), dimension(numcs,-2:0) :: colorvec
      integer :: c
      logical :: my_ok
      real(ki) :: rational2, scale2

      ampsc(:) = 0.0_ki
      !---#[ Initialize helicity amplitudes :
      !---#[ reinitialize kinematics:
      pvecs(1,:) = vecs(1,:)
      pvecs(2,:) = vecs(2,:)
      pvecs(3,:) = vecs(3,:)
      pvecs(4,:) = vecs(4,:)
      pvecs(5,:) = vecs(5,:)
      call init_event(pvecs, -1, -1, -1)
      !---#] reinitialize kinematics:
      ! For loop induced diagrams the scale should not matter
      scale2 = 100.0_ki
      do c=1,numcs
         colorvec(c,:) = samplitudeh0l1(real(scale2,ki),my_ok,rational2,c)
      end do
      heli_amp0 = colorvec(:, 0)
      !---#[ reinitialize kinematics:
      pvecs(1,:) = vecs(1,:)
      pvecs(2,:) = vecs(2,:)
      pvecs(3,:) = vecs(4,:)
      pvecs(4,:) = vecs(3,:)
      pvecs(5,:) = vecs(5,:)
      call init_event(pvecs, +1, -1, -1)
      !---#] reinitialize kinematics:
      ! For loop induced diagrams the scale should not matter
      scale2 = 100.0_ki
      do c=1,numcs
         colorvec(c,:) = samplitudeh0l1(real(scale2,ki),my_ok,rational2,c)
      end do
      heli_amp1 = colorvec(:, 0)
      !---#[ reinitialize kinematics:
      pvecs(1,:) = vecs(1,:)
      pvecs(2,:) = vecs(2,:)
      pvecs(3,:) = vecs(4,:)
      pvecs(4,:) = vecs(3,:)
      pvecs(5,:) = vecs(5,:)
      call init_event(pvecs, -1, +1, -1)
      !---#] reinitialize kinematics:
      ! For loop induced diagrams the scale should not matter
      scale2 = 100.0_ki
      do c=1,numcs
         colorvec(c,:) = samplitudeh0l1(real(scale2,ki),my_ok,rational2,c)
      end do
      heli_amp2 = colorvec(:, 0)
      !---#[ reinitialize kinematics:
      pvecs(1,:) = vecs(1,:)
      pvecs(2,:) = vecs(2,:)
      pvecs(3,:) = vecs(4,:)
      pvecs(4,:) = vecs(3,:)
      pvecs(5,:) = vecs(5,:)
      call init_event(pvecs, +1, +1, -1)
      !---#] reinitialize kinematics:
      ! For loop induced diagrams the scale should not matter
      scale2 = 100.0_ki
      do c=1,numcs
         colorvec(c,:) = samplitudeh0l1(real(scale2,ki),my_ok,rational2,c)
      end do
      heli_amp3 = colorvec(:, 0)
      !---#[ reinitialize kinematics:
      pvecs(1,:) = vecs(1,:)
      pvecs(2,:) = vecs(2,:)
      pvecs(3,:) = vecs(4,:)
      pvecs(4,:) = vecs(3,:)
      pvecs(5,:) = vecs(5,:)
      call init_event(pvecs, -1, -1, +1)
      !---#] reinitialize kinematics:
      ! For loop induced diagrams the scale should not matter
      scale2 = 100.0_ki
      do c=1,numcs
         colorvec(c,:) = samplitudeh0l1(real(scale2,ki),my_ok,rational2,c)
      end do
      heli_amp4 = colorvec(:, 0)
      !---#[ reinitialize kinematics:
      pvecs(1,:) = vecs(1,:)
      pvecs(2,:) = vecs(2,:)
      pvecs(3,:) = vecs(4,:)
      pvecs(4,:) = vecs(3,:)
      pvecs(5,:) = vecs(5,:)
      call init_event(pvecs, +1, -1, +1)
      !---#] reinitialize kinematics:
      ! For loop induced diagrams the scale should not matter
      scale2 = 100.0_ki
      do c=1,numcs
         colorvec(c,:) = samplitudeh0l1(real(scale2,ki),my_ok,rational2,c)
      end do
      heli_amp5 = colorvec(:, 0)
      !---#[ reinitialize kinematics:
      pvecs(1,:) = vecs(1,:)
      pvecs(2,:) = vecs(2,:)
      pvecs(3,:) = vecs(4,:)
      pvecs(4,:) = vecs(3,:)
      pvecs(5,:) = vecs(5,:)
      call init_event(pvecs, -1, +1, +1)
      !---#] reinitialize kinematics:
      ! For loop induced diagrams the scale should not matter
      scale2 = 100.0_ki
      do c=1,numcs
         colorvec(c,:) = samplitudeh0l1(real(scale2,ki),my_ok,rational2,c)
      end do
      heli_amp6 = colorvec(:, 0)
      !---#[ reinitialize kinematics:
      pvecs(1,:) = vecs(1,:)
      pvecs(2,:) = vecs(2,:)
      pvecs(3,:) = vecs(4,:)
      pvecs(4,:) = vecs(3,:)
      pvecs(5,:) = vecs(5,:)
      call init_event(pvecs, +1, +1, +1)
      !---#] reinitialize kinematics:
      ! For loop induced diagrams the scale should not matter
      scale2 = 100.0_ki
      do c=1,numcs
         colorvec(c,:) = samplitudeh0l1(real(scale2,ki),my_ok,rational2,c)
      end do
      heli_amp7 = colorvec(:, 0)
      !---#] Initialize helicity amplitudes :

       
      !---#[ pair 12 :

      mp  = 0.0_ki &
      &          + square_1_2_sc(heli_amp0,heli_amp1) &
      &          + square_1_2_sc(heli_amp2,heli_amp3) &
      &          + square_1_2_sc(heli_amp4,heli_amp5) &
      &          + square_1_2_sc(heli_amp6,heli_amp7)


      ampsc(2*(1-1)+2*(2-1)*num_legs+1)   = ampsc(2*(1-1)+2*(2-1)*num_legs +1) + real(mp, ki)
      ampsc(2*(1-1)+2*(2-1)*num_legs+2) = ampsc(2*(1-1)+2*(2-1)*num_legs + 2)  + real(aimag(mp),ki)

      !---#] pair 12 :
      
      !---#[ pair 15 :

      mp  = 0.0_ki &
      &          + square_1_5_sc(heli_amp0,heli_amp1) &
      &          + square_1_5_sc(heli_amp2,heli_amp3) &
      &          + square_1_5_sc(heli_amp4,heli_amp5) &
      &          + square_1_5_sc(heli_amp6,heli_amp7)


      ampsc(2*(1-1)+2*(5-1)*num_legs+1)   = ampsc(2*(1-1)+2*(5-1)*num_legs +1) + real(mp, ki)
      ampsc(2*(1-1)+2*(5-1)*num_legs+2) = ampsc(2*(1-1)+2*(5-1)*num_legs + 2)  + real(aimag(mp),ki)

      !---#] pair 15 :
      
      !---#[ pair 21 :

      mp  = 0.0_ki &
      &          + square_2_1_sc(heli_amp0,heli_amp2) &
      &          + square_2_1_sc(heli_amp1,heli_amp3) &
      &          + square_2_1_sc(heli_amp4,heli_amp6) &
      &          + square_2_1_sc(heli_amp5,heli_amp7)


      ampsc(2*(2-1)+2*(1-1)*num_legs+1)   = ampsc(2*(2-1)+2*(1-1)*num_legs +1) + real(mp, ki)
      ampsc(2*(2-1)+2*(1-1)*num_legs+2) = ampsc(2*(2-1)+2*(1-1)*num_legs + 2)  + real(aimag(mp),ki)

      !---#] pair 21 :
       
      !---#[ pair 25 :

      mp  = 0.0_ki &
      &          + square_2_5_sc(heli_amp0,heli_amp2) &
      &          + square_2_5_sc(heli_amp1,heli_amp3) &
      &          + square_2_5_sc(heli_amp4,heli_amp6) &
      &          + square_2_5_sc(heli_amp5,heli_amp7)


      ampsc(2*(2-1)+2*(5-1)*num_legs+1)   = ampsc(2*(2-1)+2*(5-1)*num_legs +1) + real(mp, ki)
      ampsc(2*(2-1)+2*(5-1)*num_legs+2) = ampsc(2*(2-1)+2*(5-1)*num_legs + 2)  + real(aimag(mp),ki)

      !---#] pair 25 :
      
      !---#[ pair 51 :

      mp  = 0.0_ki &
      &          + square_5_1_sc(heli_amp0,heli_amp4) &
      &          + square_5_1_sc(heli_amp1,heli_amp5) &
      &          + square_5_1_sc(heli_amp2,heli_amp6) &
      &          + square_5_1_sc(heli_amp3,heli_amp7)


      ampsc(2*(5-1)+2*(1-1)*num_legs+1)   = ampsc(2*(5-1)+2*(1-1)*num_legs +1) + real(mp, ki)
      ampsc(2*(5-1)+2*(1-1)*num_legs+2) = ampsc(2*(5-1)+2*(1-1)*num_legs + 2)  + real(aimag(mp),ki)

      !---#] pair 51 :
      
      !---#[ pair 52 :

      mp  = 0.0_ki &
      &          + square_5_2_sc(heli_amp0,heli_amp4) &
      &          + square_5_2_sc(heli_amp1,heli_amp5) &
      &          + square_5_2_sc(heli_amp2,heli_amp6) &
      &          + square_5_2_sc(heli_amp3,heli_amp7)


      ampsc(2*(5-1)+2*(2-1)*num_legs+1)   = ampsc(2*(5-1)+2*(2-1)*num_legs +1) + real(mp, ki)
      ampsc(2*(5-1)+2*(2-1)*num_legs+2) = ampsc(2*(5-1)+2*(2-1)*num_legs + 2)  + real(aimag(mp),ki)

      !---#] pair 52 :
       

      

      if (include_helicity_avg_factor) then
         ampsc = ampsc / real(in_helicities, ki)
      end if
      if (include_color_avg_factor) then
         ampsc = ampsc / incolors
      end if
      if (include_symmetry_factor) then
         ampsc = ampsc / real(symmetry_factor, ki)
      end if
   end subroutine OLP_spin_correlated_lo2
   !---#] spin correlated ME :


   !---#[ construct polarisation tensor :
   pure subroutine construct_polarization_tensor(eps1, eps2, tens)
      implicit none
      complex(ki), dimension(0:3), intent(in) :: eps1, eps2
      complex(ki), dimension(0:3,0:3), intent(out) :: tens

      integer :: mu, nu

      do mu = 0,3
         do nu = 0, 3
            tens(mu,nu) = eps1(mu) * eps2(nu)
         end do
      end do
   end  subroutine construct_polarization_tensor
   !---#] construct polarisation tensor :

   pure function square_0l_0l_sc(color_vector1, color_vector2) result(amp)
      use p2_part21part21_part25part25part21_color, only: cmat => CC
      implicit none
      complex(ki), dimension(numcs), intent(in) :: color_vector1, color_vector2
      complex(ki) :: amp
      complex(ki), dimension(numcs) :: v1, v2

      v1 = matmul(cmat, color_vector2)
      v2 = conjg(color_vector1)
      amp = sum(v1(:) * v2(:))
   end function  square_0l_0l_sc


    
   pure function square_1_2_sc(color_vector1, color_vector2) result(amp)
      use p2_part21part21_part25part25part21_color, only: cmat => T1T2
      
      implicit none
      complex(ki), dimension(numcs), intent(in) :: color_vector1, color_vector2
      complex(ki) :: amp
      complex(ki), dimension(numcs) :: v1, v2

      v1 = matmul(cmat, color_vector2)
      v2 = conjg(color_vector1)
      amp = sum(v1(:) * v2(:))
   end function  square_1_2_sc
      
   pure function square_1_5_sc(color_vector1, color_vector2) result(amp)
      use p2_part21part21_part25part25part21_color, only: cmat => T1T5
      
      implicit none
      complex(ki), dimension(numcs), intent(in) :: color_vector1, color_vector2
      complex(ki) :: amp
      complex(ki), dimension(numcs) :: v1, v2

      v1 = matmul(cmat, color_vector2)
      v2 = conjg(color_vector1)
      amp = sum(v1(:) * v2(:))
   end function  square_1_5_sc
      
   pure function square_2_1_sc(color_vector1, color_vector2) result(amp)
      use p2_part21part21_part25part25part21_color, only: cmat => T1T2
      
      implicit none
      complex(ki), dimension(numcs), intent(in) :: color_vector1, color_vector2
      complex(ki) :: amp
      complex(ki), dimension(numcs) :: v1, v2

      v1 = matmul(cmat, color_vector2)
      v2 = conjg(color_vector1)
      amp = sum(v1(:) * v2(:))
   end function  square_2_1_sc
       
   pure function square_2_5_sc(color_vector1, color_vector2) result(amp)
      use p2_part21part21_part25part25part21_color, only: cmat => T2T5
      
      implicit none
      complex(ki), dimension(numcs), intent(in) :: color_vector1, color_vector2
      complex(ki) :: amp
      complex(ki), dimension(numcs) :: v1, v2

      v1 = matmul(cmat, color_vector2)
      v2 = conjg(color_vector1)
      amp = sum(v1(:) * v2(:))
   end function  square_2_5_sc
      
   pure function square_5_1_sc(color_vector1, color_vector2) result(amp)
      use p2_part21part21_part25part25part21_color, only: cmat => T1T5
      
      implicit none
      complex(ki), dimension(numcs), intent(in) :: color_vector1, color_vector2
      complex(ki) :: amp
      complex(ki), dimension(numcs) :: v1, v2

      v1 = matmul(cmat, color_vector2)
      v2 = conjg(color_vector1)
      amp = sum(v1(:) * v2(:))
   end function  square_5_1_sc
      
   pure function square_5_2_sc(color_vector1, color_vector2) result(amp)
      use p2_part21part21_part25part25part21_color, only: cmat => T2T5
      
      implicit none
      complex(ki), dimension(numcs), intent(in) :: color_vector1, color_vector2
      complex(ki) :: amp
      complex(ki), dimension(numcs) :: v1, v2

      v1 = matmul(cmat, color_vector2)
      v2 = conjg(color_vector1)
      amp = sum(v1(:) * v2(:))
   end function  square_5_2_sc
       

   

   

end module p2_part21part21_part25part25part21_matrix
