!program prog
!    implicit none
!    real(8) :: val
!    real(8) :: x
!    x = 0.0_8
!    val = ki3(x)
!    write(*,*) val
!
!contains
    real(8) function ki3(x)
        use iso_c_binding
        implicit none
        
        real(8) :: x
        real(8) :: y(1)
        integer :: n
        integer :: kode
        integer :: m
        integer :: nz
        integer :: ierr

        n = 3
        m = 1
        kode = 1

        call dbskin(x, n, kode, m, y, nz, ierr)

        ki3 = y(1)
        return
    end function ki3
!end program
