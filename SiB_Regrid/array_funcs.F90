MODULE ARRAY_FUNCS

  IMPLICIT NONE

CONTAINS

  SUBROUTINE AF_transpose(arr)
    !--------------------------------------------------
    !DESC: transpose a 2-D array
    !--------------------------------------------------
    !INPUT:
    !arr: a two dimensional Fortran array
    !--------------------------------------------------
    !OUTPUT
    !the input array transposed
    !--------------------------------------------------
    !AUTHOR: Timothy W. Hilton, UC Merced, 22 Dec 2014
    !--------------------------------------------------

    real, dimension(:,:), allocatable, intent(inout) :: arr
    real, dimension(:,:), allocatable :: tmp
    integer i, j, nrow, ncol, ierr

    nrow = size(arr, 1)
    ncol = size(arr, 2)

    allocate(tmp(ncol, nrow), STAT=ierr)

    DO i=1, nrow
       DO j=1, ncol
          tmp(j,i) = arr(i,j)
       ENDDO
    ENDDO
    deallocate(arr)
    allocate(arr(ncol, nrow), STAT=ierr)
    arr(:,:) = tmp(:,:)

  END SUBROUTINE AF_transpose

!--------------------------------------------------

  SUBROUTINE fliplr(arr)
    !--------------------------------------------------
    !DESC: flip a 2-D array in the left-right dimension
    !--------------------------------------------------
    !INPUT:
    !arr: a two dimensional Fortran array
    !--------------------------------------------------
    !OUTPUT
    !the input array, flipped horizontally
    !--------------------------------------------------
    !AUTHOR: Timothy W. Hilton, UC Merced, 22 Dec 2014
    !--------------------------------------------------

    real, dimension(:,:), intent(inout) :: arr
    real, dimension(:,:), allocatable :: tmp
    integer i, j, nrow, ncol, ierr

    nrow = size(arr, 1)
    ncol = size(arr, 2)

    allocate(tmp(nrow, ncol), STAT=ierr)

    tmp(:,:) = arr(:, ncol:1:-1)
    arr = tmp
    deallocate(tmp)

  END SUBROUTINE fliplr

!--------------------------------------------------

  SUBROUTINE flipud(arr)
    !--------------------------------------------------
    !DESC: flip a 2-D array in the up-down dimension
    !--------------------------------------------------
    !INPUT:
    !arr: a two dimensional Fortran array
    !--------------------------------------------------
    !OUTPUT
    !the input array, flipped vertically
    !--------------------------------------------------
    !AUTHOR: Timothy W. Hilton, UC Merced, 22 Dec 2014
    !--------------------------------------------------

    real, dimension(:,:), intent(inout) :: arr
    real, dimension(:,:), allocatable :: tmp
    integer i, j, nrow, ncol, ierr

    nrow = size(arr, 1)
    ncol = size(arr, 2)

    allocate(tmp(nrow, ncol), STAT=ierr)

    tmp(:,:) = arr(nrow:1:-1, :)
    arr = tmp
    deallocate(tmp)

  END SUBROUTINE flipud

!--------------------------------------------------
  SUBROUTINE write_array(arr)

    real, dimension(:,:), intent(in) :: arr
    integer i, j, nrow, ncol

    nrow = size(arr, 1)
    ncol = size(arr, 2)

    DO i = 1, nrow
       DO j = 1, ncol
          write(*, '(e15.4)', advance='no') arr(i,j)
       ENDDO
       write(*,*)
    ENDDO

  END SUBROUTINE write_array

END MODULE ARRAY_FUNCS

!============================================================
