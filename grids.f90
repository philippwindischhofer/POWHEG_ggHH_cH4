subroutine initgrids(cHHH,ct,ctt,cg,cgg)
  use, intrinsic :: iso_c_binding
  implicit none
  real(c_double) :: s, t, result, expected
  real * 8 cHHH,ct,ctt,cg,cgg
  character(len=80) :: gridname
  character(len=100,kind=C_CHAR) :: c_gridname
  type(c_ptr) :: grid
  common/cbgrid/grid
  integer mtdep
  common/cmtdep/mtdep

  
  interface
     subroutine python_initialize() bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
     end subroutine python_initialize

     subroutine python_decref(grid) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), intent(in) :: grid
     end subroutine python_decref

     subroutine python_finalize() bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
     end subroutine python_finalize

     subroutine python_printinfo() bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
     end subroutine python_printinfo

     subroutine combine_grids(grid_temp, cHHH, ct, ctt, cg, cgg) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       real(c_double), intent(in), value :: cHHH, ct, ctt, cg, cgg
       character(kind=c_char) :: grid_temp
     end subroutine combine_grids

     type(c_ptr) function grid_initialize(grid_name) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       character(kind=c_char) :: grid_name(*)
     end function grid_initialize

     real(c_double) function grid_virt(grid, s, t) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: grid
       real(c_double), intent(in), value :: s,t
     end function grid_virt
  end interface

  call python_initialize
  call python_printinfo

  if (mtdep.ne.5) then
    write(gridname, "(A10,SP,ES11.4,A,ES11.4,A,ES11.4,A,ES11.4,A,ES11.4,A5)") &
           "Virt_full_", cHHH, "_", ct,"_", ctt, "_", cg, "_", cgg, ".grid"
    c_gridname = TRIM(gridname)//C_NULL_CHAR

    call combine_grids(c_gridname,cHHH,ct,ctt,cg,cgg)
  else
    c_gridname = TRIM("grid.in")//C_NULL_CHAR
  endif


  grid = grid_initialize(c_gridname)

  ! s = 2.56513D6
  ! t = -482321.D0
  ! result = grid_virt(grid,s,t)

  ! write(*,*) 'Sent: ', s,' ', t
  ! write(*,*) 'Received: ', result

end subroutine initgrids

function gridvirt(s,t)
  use, intrinsic :: iso_c_binding
  implicit none
  real(c_double) :: s, t, result, expected
  real(c_double) :: gridvirt
  type(c_ptr) :: grid
  common/cbgrid/grid

  interface
     subroutine python_initialize() bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
     end subroutine python_initialize

     subroutine python_decref(grid) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), intent(in) :: grid
     end subroutine python_decref

     subroutine python_finalize() bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
     end subroutine python_finalize

     subroutine python_printinfo() bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
     end subroutine python_printinfo

     type(c_ptr) function grid_initialize(grid_name) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       character(kind=c_char) :: grid_name(*)
     end function grid_initialize

     real(c_double) function grid_virt(grid, s, t) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: grid
       real(c_double), intent(in), value :: s,t
     end function grid_virt
  end interface
  
  gridvirt = grid_virt(grid, s, t)

end function gridvirt
  
subroutine teardowngrids()
    use, intrinsic :: iso_c_binding
  implicit none
  real(c_double) :: s, t, result, expected
  type(c_ptr) :: grid

  interface
     subroutine python_initialize() bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
     end subroutine python_initialize

     subroutine python_decref(grid) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), intent(in) :: grid
     end subroutine python_decref

     subroutine python_finalize() bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
     end subroutine python_finalize

     subroutine python_printinfo() bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
     end subroutine python_printinfo

     type(c_ptr) function grid_initialize(grid_name) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       character(kind=c_char) :: grid_name(*)
     end function grid_initialize

     real(c_double) function grid_virt(grid, s, t) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: grid
       real(c_double), intent(in), value :: s,t
     end function grid_virt
  end interface

  call python_decref(grid)
  call python_finalize

end subroutine teardowngrids
