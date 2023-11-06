module     pb_part21part21_part25part25_groups
   use precision_golem, only: ki_gol => ki
   use tens_rec
   use pb_part21part21_part25part25_config, only: ki 
   implicit none
   save

   private
!---#[ tensor coefficients for golem95:
   !-----#[ tensor coefficients group 0:
   type tensrec_info_group0
      type(coeff_type_4) :: coeffs_13
   end type

   public :: tensrec_info_group0
   !-----#] tensor coefficients group 0:
   !-----#[ tensor coefficients group 1:
   type tensrec_info_group1
      type(coeff_type_3) :: coeffs_2
      type(coeff_type_3) :: coeffs_7
      type(coeff_type_4) :: coeffs_11
   end type

   public :: tensrec_info_group1
   !-----#] tensor coefficients group 1:
   !-----#[ tensor coefficients group 2:
   type tensrec_info_group2
      type(coeff_type_4) :: coeffs_9
   end type

   public :: tensrec_info_group2
   !-----#] tensor coefficients group 2:
!---#] tensor coefficients for golem95:
   integer :: prev_ls = -1

   interface contract_golem95
      module procedure contract_tensor_coefficients_group_0
      module procedure contract_tensor_coefficients_group_1
      module procedure contract_tensor_coefficients_group_2
   end interface

   public :: contract_golem95
   public :: tear_down_golem95
   public :: ninja_exit
contains
!---#[ contract tensor coefficients golem95:
!-----#[ function contract_tensor_coefficients_group_0:
function     contract_tensor_coefficients_group_0(coeffs) result(amp)
   use matrice_s, only: allocation_s, deallocation_s, s_mat, set_ref, &
                      & s_mat_c, b_ref, preparesmatrix
   use parametre, only: rmass_or_cmass_par, cmass
   use cache, only: allocate_cache, clear_cache, reset_cache
   use array, only: packb
   use tens_comb
   use form_factor_1p, only: a10
   use form_factor_2p, only: a20
   use form_factor_3p, only: a30
   use form_factor_4p, only: a40
   use form_factor_type, only: form_factor, operator(+), operator(-)
   use pb_part21part21_part25part25_config, only: debug_nlo_diagrams, logfile
   use pb_part21part21_part25part25_kinematics, only: k1, k2, k3, k4, es12, es3, es23, es4
   use pb_part21part21_part25part25_model
   implicit none
   type(tensrec_info_group0), intent(in) :: coeffs
   type(form_factor) :: amp, dbg_amp

   integer :: b_set
   real(ki_gol), dimension(4,0:3) :: rmomenta
   logical :: ev_diagram

   if(prev_ls.ne.4) then
      if(prev_ls > 0) then
         !------#[ call sequence of exitgolem95():
         rmass_or_cmass_par = cmass
         nullify(s_mat)
         call deallocation_s()
         call clear_cache()
         !------#] call sequence of exitgolem95():
      end if

      !------#[ call sequence of initgolem95():
      rmass_or_cmass_par = cmass
      call allocation_s(4)
      set_ref = (/1, 2, 3, 4/)
      b_ref = packb(set_ref)
      call allocate_cache(4)
      s_mat => s_mat_c
      !------#] call sequence of initgolem95():
      prev_ls = 4
   !else
   !   reset_cache() is called by preparesmatrix()
   !   call reset_cache()
   end if
   s_mat(1,1)=real(-2.0_ki*mdlMT**2, ki_gol)
   s_mat(1,2)=real(-2.0_ki*mdlMT**2+mdlMh**2, ki_gol)
   s_mat(2,1)=s_mat(1,2)
   s_mat(1,3)=real(-2.0_ki*mdlMT**2-es12+2.0_ki*mdlMh**2-es23, ki_gol)
   s_mat(3,1)=s_mat(1,3)
   s_mat(1,4)=real(-2.0_ki*mdlMT**2, ki_gol)
   s_mat(4,1)=s_mat(1,4)
   s_mat(2,2)=real(-2.0_ki*mdlMT**2, ki_gol)
   s_mat(2,3)=real(-2.0_ki*mdlMT**2, ki_gol)
   s_mat(3,2)=s_mat(2,3)
   s_mat(2,4)=real(-2.0_ki*mdlMT**2+es23, ki_gol)
   s_mat(4,2)=s_mat(2,4)
   s_mat(3,3)=real(-2.0_ki*mdlMT**2, ki_gol)
   s_mat(3,4)=real(-2.0_ki*mdlMT**2+mdlMh**2, ki_gol)
   s_mat(4,3)=s_mat(3,4)
   s_mat(4,4)=real(-2.0_ki*mdlMT**2, ki_gol)
   call preparesmatrix()
   rmomenta(1,:) = real(-k2+k4, ki_gol)
   rmomenta(2,:) = real(-k2, ki_gol)
   rmomenta(3,:) = 0.0_ki_gol
   rmomenta(4,:) = real(-k3, ki_gol)
   !-------#[ Diagram 13:
      if(debug_nlo_diagrams) then
         write(logfile,*) "<diagram index='13'>"
         dbg_amp = 0.0_ki_gol
      end if
      b_set = 0

      amp = + (contract4_4(coeffs%coeffs_13, rmomenta, b_set))
      if(debug_nlo_diagrams) then
         dbg_amp = amp - dbg_amp
         write(logfile,100) &
            & "<result kind='nlo-finite' re='", real(dbg_amp%C, ki), &
            & "' im='", aimag(dbg_amp%C), "'/>"
         write(logfile,100) &
            & "<result kind='nlo-single' re='", real(dbg_amp%B, ki), &
            & "' im='", aimag(dbg_amp%B), "'/>"
         write(logfile,100) &
            & "<result kind='nlo-double' re='", real(dbg_amp%A, ki), &
            & "' im='", aimag(dbg_amp%A), "'/>"
         write(logfile,*) "</diagram>"
      end if
   !-------#] Diagram 13:

   100 format (A30,E24.16,A6,E24.16,A3)
end function contract_tensor_coefficients_group_0
!-----#] function contract_tensor_coefficients_group_0:
!-----#[ function contract_tensor_coefficients_group_1:
function     contract_tensor_coefficients_group_1(coeffs) result(amp)
   use matrice_s, only: allocation_s, deallocation_s, s_mat, set_ref, &
                      & s_mat_c, b_ref, preparesmatrix
   use parametre, only: rmass_or_cmass_par, cmass
   use cache, only: allocate_cache, clear_cache, reset_cache
   use array, only: packb
   use tens_comb
   use form_factor_1p, only: a10
   use form_factor_2p, only: a20
   use form_factor_3p, only: a30
   use form_factor_4p, only: a40
   use form_factor_type, only: form_factor, operator(+), operator(-)
   use pb_part21part21_part25part25_config, only: debug_nlo_diagrams, logfile
   use pb_part21part21_part25part25_kinematics, only: k1, k2, k3, k4, es12, es3, es23, es4
   use pb_part21part21_part25part25_model
   implicit none
   type(tensrec_info_group1), intent(in) :: coeffs
   type(form_factor) :: amp, dbg_amp

   integer :: b_set
   real(ki_gol), dimension(4,0:3) :: rmomenta
   logical :: ev_diagram

   if(prev_ls.ne.4) then
      if(prev_ls > 0) then
         !------#[ call sequence of exitgolem95():
         rmass_or_cmass_par = cmass
         nullify(s_mat)
         call deallocation_s()
         call clear_cache()
         !------#] call sequence of exitgolem95():
      end if

      !------#[ call sequence of initgolem95():
      rmass_or_cmass_par = cmass
      call allocation_s(4)
      set_ref = (/1, 2, 3, 4/)
      b_ref = packb(set_ref)
      call allocate_cache(4)
      s_mat => s_mat_c
      !------#] call sequence of initgolem95():
      prev_ls = 4
   !else
   !   reset_cache() is called by preparesmatrix()
   !   call reset_cache()
   end if
   s_mat(1,1)=real(-2.0_ki*mdlMT**2, ki_gol)
   s_mat(1,2)=real(-2.0_ki*mdlMT**2+mdlMh**2, ki_gol)
   s_mat(2,1)=s_mat(1,2)
   s_mat(1,3)=real(-2.0_ki*mdlMT**2+es12, ki_gol)
   s_mat(3,1)=s_mat(1,3)
   s_mat(1,4)=real(-2.0_ki*mdlMT**2, ki_gol)
   s_mat(4,1)=s_mat(1,4)
   s_mat(2,2)=real(-2.0_ki*mdlMT**2, ki_gol)
   s_mat(2,3)=real(-2.0_ki*mdlMT**2+mdlMh**2, ki_gol)
   s_mat(3,2)=s_mat(2,3)
   s_mat(2,4)=real(-2.0_ki*mdlMT**2+es23, ki_gol)
   s_mat(4,2)=s_mat(2,4)
   s_mat(3,3)=real(-2.0_ki*mdlMT**2, ki_gol)
   s_mat(3,4)=real(-2.0_ki*mdlMT**2, ki_gol)
   s_mat(4,3)=s_mat(3,4)
   s_mat(4,4)=real(-2.0_ki*mdlMT**2, ki_gol)
   call preparesmatrix()
   rmomenta(1,:) = real(-k3-k4, ki_gol)
   rmomenta(2,:) = real(-k3, ki_gol)
   rmomenta(3,:) = 0.0_ki_gol
   rmomenta(4,:) = real(-k2, ki_gol)
   !-------#[ Diagram 2:
      if(debug_nlo_diagrams) then
         write(logfile,*) "<diagram index='2'>"
         dbg_amp = 0.0_ki_gol
      end if
      b_set = packb((/2/))

      amp = + (contract3_3(coeffs%coeffs_2, rmomenta, b_set))
      if(debug_nlo_diagrams) then
         dbg_amp = amp - dbg_amp
         write(logfile,100) &
            & "<result kind='nlo-finite' re='", real(dbg_amp%C, ki), &
            & "' im='", aimag(dbg_amp%C), "'/>"
         write(logfile,100) &
            & "<result kind='nlo-single' re='", real(dbg_amp%B, ki), &
            & "' im='", aimag(dbg_amp%B), "'/>"
         write(logfile,100) &
            & "<result kind='nlo-double' re='", real(dbg_amp%A, ki), &
            & "' im='", aimag(dbg_amp%A), "'/>"
         write(logfile,*) "</diagram>"
      end if
   !-------#] Diagram 2:
   !-------#[ Diagram 7:
      if(debug_nlo_diagrams) then
         write(logfile,*) "<diagram index='7'>"
         dbg_amp = amp
      end if
      b_set = packb((/2/))

      amp = amp + (contract3_3(coeffs%coeffs_7, rmomenta, b_set))
      if(debug_nlo_diagrams) then
         dbg_amp = amp - dbg_amp
         write(logfile,100) &
            & "<result kind='nlo-finite' re='", real(dbg_amp%C, ki), &
            & "' im='", aimag(dbg_amp%C), "'/>"
         write(logfile,100) &
            & "<result kind='nlo-single' re='", real(dbg_amp%B, ki), &
            & "' im='", aimag(dbg_amp%B), "'/>"
         write(logfile,100) &
            & "<result kind='nlo-double' re='", real(dbg_amp%A, ki), &
            & "' im='", aimag(dbg_amp%A), "'/>"
         write(logfile,*) "</diagram>"
      end if
   !-------#] Diagram 7:
   !-------#[ Diagram 11:
      if(debug_nlo_diagrams) then
         write(logfile,*) "<diagram index='11'>"
         dbg_amp = amp
      end if
      b_set = 0

      amp = amp + (contract4_4(coeffs%coeffs_11, rmomenta, b_set))
      if(debug_nlo_diagrams) then
         dbg_amp = amp - dbg_amp
         write(logfile,100) &
            & "<result kind='nlo-finite' re='", real(dbg_amp%C, ki), &
            & "' im='", aimag(dbg_amp%C), "'/>"
         write(logfile,100) &
            & "<result kind='nlo-single' re='", real(dbg_amp%B, ki), &
            & "' im='", aimag(dbg_amp%B), "'/>"
         write(logfile,100) &
            & "<result kind='nlo-double' re='", real(dbg_amp%A, ki), &
            & "' im='", aimag(dbg_amp%A), "'/>"
         write(logfile,*) "</diagram>"
      end if
   !-------#] Diagram 11:

   100 format (A30,E24.16,A6,E24.16,A3)
end function contract_tensor_coefficients_group_1
!-----#] function contract_tensor_coefficients_group_1:
!-----#[ function contract_tensor_coefficients_group_2:
function     contract_tensor_coefficients_group_2(coeffs) result(amp)
   use matrice_s, only: allocation_s, deallocation_s, s_mat, set_ref, &
                      & s_mat_c, b_ref, preparesmatrix
   use parametre, only: rmass_or_cmass_par, cmass
   use cache, only: allocate_cache, clear_cache, reset_cache
   use array, only: packb
   use tens_comb
   use form_factor_1p, only: a10
   use form_factor_2p, only: a20
   use form_factor_3p, only: a30
   use form_factor_4p, only: a40
   use form_factor_type, only: form_factor, operator(+), operator(-)
   use pb_part21part21_part25part25_config, only: debug_nlo_diagrams, logfile
   use pb_part21part21_part25part25_kinematics, only: k1, k2, k3, k4, es12, es3, es23, es4
   use pb_part21part21_part25part25_model
   implicit none
   type(tensrec_info_group2), intent(in) :: coeffs
   type(form_factor) :: amp, dbg_amp

   integer :: b_set
   real(ki_gol), dimension(4,0:3) :: rmomenta
   logical :: ev_diagram

   if(prev_ls.ne.4) then
      if(prev_ls > 0) then
         !------#[ call sequence of exitgolem95():
         rmass_or_cmass_par = cmass
         nullify(s_mat)
         call deallocation_s()
         call clear_cache()
         !------#] call sequence of exitgolem95():
      end if

      !------#[ call sequence of initgolem95():
      rmass_or_cmass_par = cmass
      call allocation_s(4)
      set_ref = (/1, 2, 3, 4/)
      b_ref = packb(set_ref)
      call allocate_cache(4)
      s_mat => s_mat_c
      !------#] call sequence of initgolem95():
      prev_ls = 4
   !else
   !   reset_cache() is called by preparesmatrix()
   !   call reset_cache()
   end if
   s_mat(1,1)=real(-2.0_ki*mdlMT**2, ki_gol)
   s_mat(1,2)=real(-2.0_ki*mdlMT**2+mdlMh**2, ki_gol)
   s_mat(2,1)=s_mat(1,2)
   s_mat(1,3)=real(-2.0_ki*mdlMT**2+es12, ki_gol)
   s_mat(3,1)=s_mat(1,3)
   s_mat(1,4)=real(-2.0_ki*mdlMT**2, ki_gol)
   s_mat(4,1)=s_mat(1,4)
   s_mat(2,2)=real(-2.0_ki*mdlMT**2, ki_gol)
   s_mat(2,3)=real(-2.0_ki*mdlMT**2+mdlMh**2, ki_gol)
   s_mat(3,2)=s_mat(2,3)
   s_mat(2,4)=real(-2.0_ki*mdlMT**2-es12+2.0_ki*mdlMh**2-es23, ki_gol)
   s_mat(4,2)=s_mat(2,4)
   s_mat(3,3)=real(-2.0_ki*mdlMT**2, ki_gol)
   s_mat(3,4)=real(-2.0_ki*mdlMT**2, ki_gol)
   s_mat(4,3)=s_mat(3,4)
   s_mat(4,4)=real(-2.0_ki*mdlMT**2, ki_gol)
   call preparesmatrix()
   rmomenta(1,:) = real(-k3-k4, ki_gol)
   rmomenta(2,:) = real(-k4, ki_gol)
   rmomenta(3,:) = 0.0_ki_gol
   rmomenta(4,:) = real(-k2, ki_gol)
   !-------#[ Diagram 9:
      if(debug_nlo_diagrams) then
         write(logfile,*) "<diagram index='9'>"
         dbg_amp = 0.0_ki_gol
      end if
      b_set = 0

      amp = + (contract4_4(coeffs%coeffs_9, rmomenta, b_set))
      if(debug_nlo_diagrams) then
         dbg_amp = amp - dbg_amp
         write(logfile,100) &
            & "<result kind='nlo-finite' re='", real(dbg_amp%C, ki), &
            & "' im='", aimag(dbg_amp%C), "'/>"
         write(logfile,100) &
            & "<result kind='nlo-single' re='", real(dbg_amp%B, ki), &
            & "' im='", aimag(dbg_amp%B), "'/>"
         write(logfile,100) &
            & "<result kind='nlo-double' re='", real(dbg_amp%A, ki), &
            & "' im='", aimag(dbg_amp%A), "'/>"
         write(logfile,*) "</diagram>"
      end if
   !-------#] Diagram 9:

   100 format (A30,E24.16,A6,E24.16,A3)
end function contract_tensor_coefficients_group_2
!-----#] function contract_tensor_coefficients_group_2:
!---#] contract tensor coefficients golem95:
!---#[ subroutine tear_down_golem95:
subroutine     tear_down_golem95()
   use matrice_s, only: deallocation_s, s_mat
   use parametre, only: rmass_or_cmass_par, cmass
   use cache, only: clear_cache
   implicit none
   if(prev_ls.gt.0) then
      !------#[ call sequence of exitgolem95():
      rmass_or_cmass_par = cmass
      nullify(s_mat)
      call deallocation_s()
      call clear_cache()
      !------#] call sequence of exitgolem95():
      prev_ls = -1
   end if
end subroutine tear_down_golem95
!---#] subroutine tear_down_golem95:
!---#[ subroutine ninja_exit:
subroutine ninja_exit()
  use ninjago_module, only: ninja_clear_integral_cache
  implicit none
  !------#[ call ninja_clear_integral_cache():
  call ninja_clear_integral_cache()
  !------#] call ninja_clear_integral_cache():
end subroutine ninja_exit
!---#] subroutine ninja_exit:
end module pb_part21part21_part25part25_groups
