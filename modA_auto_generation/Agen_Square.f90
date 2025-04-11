module Agen
  implicit none

  integer :: modA_npols
  parameter ( modA_npols = 1 ) !number of polygons

  !number of vertices (or borders) of each polygon
  integer :: modA_nvpols(modA_npols)
  data modA_nvpols(1:modA_npols) / 5 /

  !Acquired coordinates of each boundary points
  real(kind=8) :: modA_putx(1, 5),  modA_puty(1, 5)
  data modA_putx(1,1:5) / &
       -0.5000d0, 0.5000d0, 0.5000d0, -0.5000d0, -0.5000d0 /

  data modA_puty(1,1:5) / &
       -0.5000d0, -0.5000d0, 0.5000d0, 0.5000d0, -0.5000d0 /

  !Edges ID
  integer :: modA_ploedges(1, 4)
  data modA_ploedges(1,1:4) / &
       0, 0, 0, 0 /

  ! Rectangle containing A
  real(kind=8) :: modA_dbl(2), modA_dtr(2)
  data modA_dbl(1:2) / -0.6000d0, -0.6000d0 /, &
       modA_dtr(1:2) / 0.6000d0, 0.6000d0 /

end module Agen
