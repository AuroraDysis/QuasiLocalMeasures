#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



! TODO:
! instead of interpolating to the symmetry points, copy them



! A convenient shortcut
#define P(x) CCTK_PointerTo(x)



subroutine qlm_interpolate (CCTK_ARGUMENTS, hn)
  use cctk
  use qlm_variables
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: hn
  
  CCTK_INT,  parameter :: izero = 0
  integer,   parameter :: ik = kind(izero)
  integer,   parameter :: sk = kind(interpolator)
  CCTK_REAL, parameter :: one = 1
  
  integer      :: len_coordsystem
  integer      :: len_interpolator
  integer      :: len_interpolator_options
  
  character    :: fort_coordsystem*100
  character    :: fort_interpolator*100
  character    :: fort_interpolator_options*1000
  
  integer      :: nvars
  
  integer      :: coord_handle
  integer      :: interp_handle
  integer      :: options_table
  
  integer      :: ninputs
  integer      :: noutputs
  
  CCTK_REAL, allocatable :: xcoord(:,:)
  CCTK_REAL, allocatable :: ycoord(:,:)
  CCTK_REAL, allocatable :: zcoord(:,:)
  
  integer      :: ind_gxx, ind_gxy, ind_gxz, ind_gyy, ind_gyz, ind_gzz
  integer      :: ind_kxx, ind_kxy, ind_kxz, ind_kyy, ind_kyz, ind_kzz
  integer      :: ind_alpha
  integer      :: ind_betax, ind_betay, ind_betaz
  
  integer      :: coord_type
  CCTK_POINTER :: coords(3)
  CCTK_INT     :: inputs(16)
  CCTK_INT     :: output_types(88)
  CCTK_POINTER :: outputs(88)
  CCTK_INT     :: operand_indices(88)
  CCTK_INT     :: operation_codes(88)
  integer      :: npoints
  
  character    :: msg*1000
  
  integer      :: ni, nj
  
  integer      :: ierr
  
  
  
  if (veryverbose/=0) then
     call CCTK_INFO ("Interpolating 3d grid functions")
  end if
  
  
  
  if (shift_state==0) then
     call CCTK_WARN (0, "The shift must have storage")
  end if
  
  
  
  ! Get coordinate system
  call CCTK_FortranString &
       (len_coordsystem, int(coordsystem,sk), fort_coordsystem)
  call CCTK_CoordSystemHandle (coord_handle, fort_coordsystem)
  if (coord_handle<0) then
     write (msg, '("The coordinate system """, a, """ does not exist")') &
          trim(fort_coordsystem)
     call CCTK_WARN (0, msg)
  end if
  
  ! Get interpolator
  call CCTK_FortranString &
       (len_interpolator, int(interpolator,sk), fort_interpolator)
  call CCTK_InterpHandle (interp_handle, fort_interpolator)
  if (interp_handle<0) then
     write (msg, '("The interpolator """,a,""" does not exist")') &
          trim(fort_interpolator)
     call CCTK_WARN (0, msg)
  end if
  
  ! Get interpolator options
  call CCTK_FortranString &
       (len_interpolator_options, int(interpolator_options,sk), &
       fort_interpolator_options)
  call Util_TableCreateFromString (options_table, fort_interpolator_options)
  if (options_table<0) then
     write (msg, '("The interpolator_options """,a,""" have a wrong syntax")') &
          trim(fort_interpolator_options)
     call CCTK_WARN (0, msg)
  end if
  
  
  
  if (hn > 0) then
     
     ni = qlm_ntheta(hn)
     nj = qlm_nphi(hn)
     
     allocate (xcoord(ni,nj))
     allocate (ycoord(ni,nj))
     allocate (zcoord(ni,nj))
     
     xcoord(:,:) = qlm_x(:ni,:nj,hn)
     ycoord(:,:) = qlm_y(:ni,:nj,hn)
     zcoord(:,:) = qlm_z(:ni,:nj,hn)
     
  end if
  
  
  
  ! TODO: check the excision mask
  
  ! Get variable indices
  call CCTK_VarIndex (ind_gxx  , "ADMBase::gxx"  )
  call CCTK_VarIndex (ind_gxy  , "ADMBase::gxy"  )
  call CCTK_VarIndex (ind_gxz  , "ADMBase::gxz"  )
  call CCTK_VarIndex (ind_gyy  , "ADMBase::gyy"  )
  call CCTK_VarIndex (ind_gyz  , "ADMBase::gyz"  )
  call CCTK_VarIndex (ind_gzz  , "ADMBase::gzz"  )
  call CCTK_VarIndex (ind_kxx  , "ADMBase::kxx"  )
  call CCTK_VarIndex (ind_kxy  , "ADMBase::kxy"  )
  call CCTK_VarIndex (ind_kxz  , "ADMBase::kxz"  )
  call CCTK_VarIndex (ind_kyy  , "ADMBase::kyy"  )
  call CCTK_VarIndex (ind_kyz  , "ADMBase::kyz"  )
  call CCTK_VarIndex (ind_kzz  , "ADMBase::kzz"  )
  call CCTK_VarIndex (ind_alpha, "ADMBase::alp"  )
  call CCTK_VarIndex (ind_betax, "ADMBase::betax")
  call CCTK_VarIndex (ind_betay, "ADMBase::betay")
  call CCTK_VarIndex (ind_betaz, "ADMBase::betaz")
  
  
  
  ! Set up the interpolator arguments
  coord_type = CCTK_VARIABLE_REAL
  if (hn > 0) then
     npoints = ni * nj
     coords(:) = (/ P(xcoord), P(ycoord), P(zcoord) /)
  else
     npoints = 0
     coords(:) = CCTK_NullPointer()
  end if
  
  inputs = (/ &
       ind_gxx, ind_gxy, ind_gxz, ind_gyy, ind_gyz, ind_gzz, &
       ind_kxx, ind_kxy, ind_kxz, ind_kyy, ind_kyz, ind_kzz, &
       ind_alpha, &
       ind_betax, ind_betay, ind_betaz /)
  call CCTK_NumVars (nvars)
  if (nvars < 0) call CCTK_WARN (0, "internal error")
  if (any(inputs /= -1 .and. (inputs < 0 .or. inputs >= nvars))) then
     call CCTK_WARN (0, "internal error")
  end if
  
  operand_indices = (/ &
       00, 01, 02, 03, 04, 05, & ! g_ij
       00, 01, 02, 03, 04, 05, & ! g_ij,k
       00, 01, 02, 03, 04, 05, &
       00, 01, 02, 03, 04, 05, &
       00, 01, 02, 03, 04, 05, & ! g_ij,kl
       00, 01, 02, 03, 04, 05, &
       00, 01, 02, 03, 04, 05, &
       00, 01, 02, 03, 04, 05, &
       00, 01, 02, 03, 04, 05, &
       00, 01, 02, 03, 04, 05, &
       06, 07, 08, 09, 10, 11, & ! K_ij
       06, 07, 08, 09, 10, 11, & ! K_ij,k
       06, 07, 08, 09, 10, 11, &
       06, 07, 08, 09, 10, 11, &
       12, &                    ! alp
       13, 14, 15 /)            ! beta^i
  
  operation_codes = (/ &
       0, 0, 0, 0, 0, 0, &      ! g_ij
       1, 1, 1, 1, 1, 1, &      ! g_ij,k
       2, 2, 2, 2, 2, 2, &
       3, 3, 3, 3, 3, 3, &
       11, 11, 11, 11, 11, 11, & ! g_ij,kl
       12, 12, 12, 12, 12, 12, &
       13, 13, 13, 13, 13, 13, &
       22, 22, 22, 22, 22, 22, &
       23, 23, 23, 23, 23, 23, &
       33, 33, 33, 33, 33, 33, &
       0, 0, 0, 0, 0, 0, &      ! K_ij
       1, 1, 1, 1, 1, 1, &      ! K_ij,k
       2, 2, 2, 2, 2, 2, &
       3, 3, 3, 3, 3, 3, &
       0, &                     ! alp
       0, 0, 0 /)               ! beta^i
  
  output_types(:) = CCTK_VARIABLE_REAL
  if (hn > 0) then
     outputs = (/ &
          P(qlm_gxx), P(qlm_gxy), P(qlm_gxz), P(qlm_gyy), P(qlm_gyz), P(qlm_gzz), &
          P(qlm_dgxxx), P(qlm_dgxyx), P(qlm_dgxzx), P(qlm_dgyyx), P(qlm_dgyzx), P(qlm_dgzzx), &
          P(qlm_dgxxy), P(qlm_dgxyy), P(qlm_dgxzy), P(qlm_dgyyy), P(qlm_dgyzy), P(qlm_dgzzy), &
          P(qlm_dgxxz), P(qlm_dgxyz), P(qlm_dgxzz), P(qlm_dgyyz), P(qlm_dgyzz), P(qlm_dgzzz), &
          P(qlm_ddgxxxx), P(qlm_ddgxyxx), P(qlm_ddgxzxx), P(qlm_ddgyyxx), P(qlm_ddgyzxx), P(qlm_ddgzzxx), &
          P(qlm_ddgxxxy), P(qlm_ddgxyxy), P(qlm_ddgxzxy), P(qlm_ddgyyxy), P(qlm_ddgyzxy), P(qlm_ddgzzxy), &
          P(qlm_ddgxxxz), P(qlm_ddgxyxz), P(qlm_ddgxzxz), P(qlm_ddgyyxz), P(qlm_ddgyzxz), P(qlm_ddgzzxz), &
          P(qlm_ddgxxyy), P(qlm_ddgxyyy), P(qlm_ddgxzyy), P(qlm_ddgyyyy), P(qlm_ddgyzyy), P(qlm_ddgzzyy), &
          P(qlm_ddgxxyz), P(qlm_ddgxyyz), P(qlm_ddgxzyz), P(qlm_ddgyyyz), P(qlm_ddgyzyz), P(qlm_ddgzzyz), &
          P(qlm_ddgxxzz), P(qlm_ddgxyzz), P(qlm_ddgxzzz), P(qlm_ddgyyzz), P(qlm_ddgyzzz), P(qlm_ddgzzzz), &
          P(qlm_kxx), P(qlm_kxy), P(qlm_kxz), P(qlm_kyy), P(qlm_kyz), P(qlm_kzz), &
          P(qlm_dkxxx), P(qlm_dkxyx), P(qlm_dkxzx), P(qlm_dkyyx), P(qlm_dkyzx), P(qlm_dkzzx), &
          P(qlm_dkxxy), P(qlm_dkxyy), P(qlm_dkxzy), P(qlm_dkyyy), P(qlm_dkyzy), P(qlm_dkzzy), &
          P(qlm_dkxxz), P(qlm_dkxyz), P(qlm_dkxzz), P(qlm_dkyyz), P(qlm_dkyzz), P(qlm_dkzzz), &
          P(qlm_alpha), &
          P(qlm_betax), P(qlm_betay), P(qlm_betaz) /)
  else
     outputs(:) = CCTK_NullPointer()
  end if
  
  
  
  ninputs = size(inputs)
  noutputs = size(outputs)
  
  
  
#if 0
  ! Poison the output variables
#define poison -42
  qlm_gxx = poison
  qlm_gxy = poison
  qlm_gxz = poison
  qlm_gyy = poison
  qlm_gyz = poison
  qlm_gzz = poison
  qlm_dgxxx = poison
  qlm_dgxyx = poison
  qlm_dgxzx = poison
  qlm_dgyyx = poison
  qlm_dgyzx = poison
  qlm_dgzzx = poison
  qlm_dgxxy = poison
  qlm_dgxyy = poison
  qlm_dgxzy = poison
  qlm_dgyyy = poison
  qlm_dgyzy = poison
  qlm_dgzzy = poison
  qlm_dgxxz = poison
  qlm_dgxyz = poison
  qlm_dgxzz = poison
  qlm_dgyyz = poison
  qlm_dgyzz = poison
  qlm_dgzzz = poison
  qlm_ddgxxxx = poison
  qlm_ddgxyxx = poison
  qlm_ddgxzxx = poison
  qlm_ddgyyxx = poison
  qlm_ddgyzxx = poison
  qlm_ddgzzxx = poison
  qlm_ddgxxxy = poison
  qlm_ddgxyxy = poison
  qlm_ddgxzxy = poison
  qlm_ddgyyxy = poison
  qlm_ddgyzxy = poison
  qlm_ddgzzxy = poison
  qlm_ddgxxxz = poison
  qlm_ddgxyxz = poison
  qlm_ddgxzxz = poison
  qlm_ddgyyxz = poison
  qlm_ddgyzxz = poison
  qlm_ddgzzxz = poison
  qlm_ddgxxyy = poison
  qlm_ddgxyyy = poison
  qlm_ddgxzyy = poison
  qlm_ddgyyyy = poison
  qlm_ddgyzyy = poison
  qlm_ddgzzyy = poison
  qlm_ddgxxyz = poison
  qlm_ddgxyyz = poison
  qlm_ddgxzyz = poison
  qlm_ddgyyyz = poison
  qlm_ddgyzyz = poison
  qlm_ddgzzyz = poison
  qlm_ddgxxzz = poison
  qlm_ddgxyzz = poison
  qlm_ddgxzzz = poison
  qlm_ddgyyzz = poison
  qlm_ddgyzzz = poison
  qlm_ddgzzzz = poison
  qlm_kxx = poison
  qlm_kxy = poison
  qlm_kxz = poison
  qlm_kyy = poison
  qlm_kyz = poison
  qlm_kzz = poison
  qlm_dkxxx = poison
  qlm_dkxyx = poison
  qlm_dkxzx = poison
  qlm_dkyyx = poison
  qlm_dkyzx = poison
  qlm_dkzzx = poison
  qlm_dkxxy = poison
  qlm_dkxyy = poison
  qlm_dkxzy = poison
  qlm_dkyyy = poison
  qlm_dkyzy = poison
  qlm_dkzzy = poison
  qlm_dkxxz = poison
  qlm_dkxyz = poison
  qlm_dkxzz = poison
  qlm_dkyyz = poison
  qlm_dkyzz = poison
  qlm_dkzzz = poison
  qlm_alpha = poison
  qlm_betax = poison
  qlm_betay = poison
  qlm_betaz = poison
#endif
  
  
  
  ! Call the interpolator
  call Util_TableSetIntArray &
       (ierr, options_table, noutputs, &
       operand_indices, "operand_indices")
  if (ierr /= 0) call CCTK_WARN (0, "internal error")
  call Util_TableSetIntArray &
       (ierr, options_table, noutputs, &
       operation_codes, "operation_codes")
  if (ierr /= 0) call CCTK_WARN (0, "internal error")
  
  call CCTK_InterpGridArrays &
       (ierr, cctkGH, 3, &
       interp_handle, options_table, coord_handle, &
       npoints, coord_type, coords, &
       ninputs, inputs, &
       noutputs, output_types, outputs)
  
  if (ierr /= 0) then
     if (hn > 0) then
        qlm_calc_error(hn) = 1
     end if
     call CCTK_WARN (1, "Interpolator failed")
     return
  end if
  
  
  
  ! Unpack the variables
  if (hn > 0) then
     
     call unpack (qlm_gxx         , ni, nj)
     call unpack (qlm_gxy         , ni, nj)
     call unpack (qlm_gxz         , ni, nj)
     call unpack (qlm_gyy         , ni, nj)
     call unpack (qlm_gyz         , ni, nj)
     call unpack (qlm_gzz         , ni, nj)
     call unpack (qlm_dgxxx       , ni, nj)
     call unpack (qlm_dgxyx       , ni, nj)
     call unpack (qlm_dgxzx       , ni, nj)
     call unpack (qlm_dgyyx       , ni, nj)
     call unpack (qlm_dgyzx       , ni, nj)
     call unpack (qlm_dgzzx       , ni, nj)
     call unpack (qlm_dgxxy       , ni, nj)
     call unpack (qlm_dgxyy       , ni, nj)
     call unpack (qlm_dgxzy       , ni, nj)
     call unpack (qlm_dgyyy       , ni, nj)
     call unpack (qlm_dgyzy       , ni, nj)
     call unpack (qlm_dgzzy       , ni, nj)
     call unpack (qlm_dgxxz       , ni, nj)
     call unpack (qlm_dgxyz       , ni, nj)
     call unpack (qlm_dgxzz       , ni, nj)
     call unpack (qlm_dgyyz       , ni, nj)
     call unpack (qlm_dgyzz       , ni, nj)
     call unpack (qlm_dgzzz       , ni, nj)
     call unpack (qlm_ddgxxxx     , ni, nj)
     call unpack (qlm_ddgxyxx     , ni, nj)
     call unpack (qlm_ddgxzxx     , ni, nj)
     call unpack (qlm_ddgyyxx     , ni, nj)
     call unpack (qlm_ddgyzxx     , ni, nj)
     call unpack (qlm_ddgzzxx     , ni, nj)
     call unpack (qlm_ddgxxxy     , ni, nj)
     call unpack (qlm_ddgxyxy     , ni, nj)
     call unpack (qlm_ddgxzxy     , ni, nj)
     call unpack (qlm_ddgyyxy     , ni, nj)
     call unpack (qlm_ddgyzxy     , ni, nj)
     call unpack (qlm_ddgzzxy     , ni, nj)
     call unpack (qlm_ddgxxxz     , ni, nj)
     call unpack (qlm_ddgxyxz     , ni, nj)
     call unpack (qlm_ddgxzxz     , ni, nj)
     call unpack (qlm_ddgyyxz     , ni, nj)
     call unpack (qlm_ddgyzxz     , ni, nj)
     call unpack (qlm_ddgzzxz     , ni, nj)
     call unpack (qlm_ddgxxyy     , ni, nj)
     call unpack (qlm_ddgxyyy     , ni, nj)
     call unpack (qlm_ddgxzyy     , ni, nj)
     call unpack (qlm_ddgyyyy     , ni, nj)
     call unpack (qlm_ddgyzyy     , ni, nj)
     call unpack (qlm_ddgzzyy     , ni, nj)
     call unpack (qlm_ddgxxyz     , ni, nj)
     call unpack (qlm_ddgxyyz     , ni, nj)
     call unpack (qlm_ddgxzyz     , ni, nj)
     call unpack (qlm_ddgyyyz     , ni, nj)
     call unpack (qlm_ddgyzyz     , ni, nj)
     call unpack (qlm_ddgzzyz     , ni, nj)
     call unpack (qlm_ddgxxzz     , ni, nj)
     call unpack (qlm_ddgxyzz     , ni, nj)
     call unpack (qlm_ddgxzzz     , ni, nj)
     call unpack (qlm_ddgyyzz     , ni, nj)
     call unpack (qlm_ddgyzzz     , ni, nj)
     call unpack (qlm_ddgzzzz     , ni, nj)
     call unpack (qlm_kxx         , ni, nj)
     call unpack (qlm_kxy         , ni, nj)
     call unpack (qlm_kxz         , ni, nj)
     call unpack (qlm_kyy         , ni, nj)
     call unpack (qlm_kyz         , ni, nj)
     call unpack (qlm_kzz         , ni, nj)
     call unpack (qlm_dkxxx       , ni, nj)
     call unpack (qlm_dkxyx       , ni, nj)
     call unpack (qlm_dkxzx       , ni, nj)
     call unpack (qlm_dkyyx       , ni, nj)
     call unpack (qlm_dkyzx       , ni, nj)
     call unpack (qlm_dkzzx       , ni, nj)
     call unpack (qlm_dkxxy       , ni, nj)
     call unpack (qlm_dkxyy       , ni, nj)
     call unpack (qlm_dkxzy       , ni, nj)
     call unpack (qlm_dkyyy       , ni, nj)
     call unpack (qlm_dkyzy       , ni, nj)
     call unpack (qlm_dkzzy       , ni, nj)
     call unpack (qlm_dkxxz       , ni, nj)
     call unpack (qlm_dkxyz       , ni, nj)
     call unpack (qlm_dkxzz       , ni, nj)
     call unpack (qlm_dkyyz       , ni, nj)
     call unpack (qlm_dkyzz       , ni, nj)
     call unpack (qlm_dkzzz       , ni, nj)
     call unpack (qlm_alpha       , ni, nj)
     call unpack (qlm_betax       , ni, nj)
     call unpack (qlm_betay       , ni, nj)
     call unpack (qlm_betaz       , ni, nj)
     
     
     
#if 0
     ! Check for poison
     if (any(qlm_gxx == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_gxy == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_gxz == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_gyy == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_gyz == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_gzz == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_dgxxx == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_dgxyx == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_dgxzx == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_dgyyx == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_dgyzx == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_dgzzx == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_dgxxy == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_dgxyy == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_dgxzy == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_dgyyy == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_dgyzy == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_dgzzy == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_dgxxz == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_dgxyz == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_dgxzz == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_dgyyz == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_dgyzz == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_dgzzz == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_ddgxxxx == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_ddgxyxx == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_ddgxzxx == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_ddgyyxx == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_ddgyzxx == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_ddgzzxx == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_ddgxxxy == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_ddgxyxy == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_ddgxzxy == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_ddgyyxy == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_ddgyzxy == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_ddgzzxy == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_ddgxxxz == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_ddgxyxz == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_ddgxzxz == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_ddgyyxz == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_ddgyzxz == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_ddgzzxz == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_ddgxxyy == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_ddgxyyy == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_ddgxzyy == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_ddgyyyy == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_ddgyzyy == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_ddgzzyy == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_ddgxxyz == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_ddgxyyz == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_ddgxzyz == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_ddgyyyz == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_ddgyzyz == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_ddgzzyz == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_ddgxxzz == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_ddgxyzz == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_ddgxzzz == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_ddgyyzz == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_ddgyzzz == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_ddgzzzz == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_kxx == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_kxy == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_kxz == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_kyy == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_kyz == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_kzz == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_dkxxx == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_dkxyx == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_dkxzx == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_dkyyx == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_dkyzx == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_dkzzx == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_dkxxy == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_dkxyy == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_dkxzy == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_dkyyy == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_dkyzy == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_dkzzy == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_dkxxz == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_dkxyz == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_dkxzz == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_dkyyz == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_dkyzz == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_dkzzz == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_alpha == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_betax == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_betay == poison)) call CCTK_WARN (0, "poison found")
     if (any(qlm_betaz == poison)) call CCTK_WARN (0, "poison found")
#endif
     
  end if
  
  
  
  ! Free interpolator options
  call Util_TableDestroy (ierr, options_table)
  
  
  
  if (hn > 0) then
     
     qlm_have_valid_data(hn) = 1
     
     deallocate (xcoord)
     deallocate (ycoord)
     deallocate (zcoord)
     
  end if
  
  
  
contains
  
  subroutine pack (arr, ni, nj)
    integer,   intent(in)    :: ni, nj
    CCTK_REAL, intent(inout) :: arr(:,:)
    CCTK_REAL :: tmp(ni,nj)
    tmp(:,:) = arr(:ni, :nj)
    call copy (arr, tmp, size(tmp))
  end subroutine pack
  
  subroutine unpack (arr, ni, nj)
    integer,   intent(in)    :: ni, nj
    CCTK_REAL, intent(inout) :: arr(:,:)
    CCTK_REAL :: tmp(ni,nj)
    call copy (tmp, arr, size(tmp))
    arr(:ni, :nj) = tmp(:,:)
    arr(ni+1:, :nj) = 0
    arr(:, nj+1:) = 0
  end subroutine unpack
  
  subroutine copy (a, b, n)
    integer,   intent(in)  :: n
    CCTK_REAL, intent(out) :: a(n)
    CCTK_REAL, intent(in)  :: b(n)
    a = b
  end subroutine copy
  
end subroutine qlm_interpolate
