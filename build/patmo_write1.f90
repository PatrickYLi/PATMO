program create_netcdf
    use netcdf
    implicit none
    integer :: status, ncid

    status = nf90_create('data.nc', NF90_NETCDF4, ncid)
    call check(status, 'open')

    status = nf90_close(ncid)
    call check(status, 'close')

contains

    subroutine check(status, operation)
        ....
    end subroutine check
end program create_netcdf  
