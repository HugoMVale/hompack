module hompack90_global
!! Global allocatable arrays used for the sparse matrix data structures
!! and by the polynomial system solver.  Used by module homotopy.

   use hompack_kinds, only: dp

   integer, dimension(:), allocatable :: colpos, ipar, rowpos
   real(dp), dimension(:), allocatable :: par, pp, qrsparse

end module hompack90_global
