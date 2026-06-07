module hompack_global_legacy
!! Global allocatable arrays used for the sparse matrix data structures
!! and by the polynomial system solver.  Used by module homotopy.

   use hompack_kinds, only: dp
   implicit none

   integer, dimension(:), allocatable :: colpos, ipar, rowpos
   real(dp), dimension(:), allocatable :: par, pp, qrsparse

end module hompack_global_legacy
