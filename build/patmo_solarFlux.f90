!This program is aim to give the output of solar flux at different altitudes
!*****Hope it works :)*****
subroutine patmo_dumpPhotoFlux(fname)
    use patmo_commons
    use patmo_parameters
    use patmo_utils
    implicit none
    character(len=*),intent(in)::fname
    real*8::zkm(cellsNumber)
    integer::i

    !upper layer opacity is zero
    tauAll(:,cellsNumber) = 0d0
    !loop on cells
    do j=cellsNumber-1,1,-1
       sumxn(:) = 0d0
       !loop on reactions
       do i=1,photoReactionsNumber
          sumxn(:) = sumxn(:) + xsecAll(:,i) * nall(j,photoPartnerIndex(i))
       end do
       tauAll(:,j) = tauAll(:,j+1) + gridSpace(j) * sumxn(:)

    open(70,file=trim(fname),status="replace")
    write(70,*)"altitude/km, 1, 10, 20, 30, 40, 50, 60"

    !output of solar flux at diff. altitudes
       
       write(70,*) &
       photoFlux 

    end do
    
    write(70,*)
    close(70)

 end subroutine patmo_dumpPhotoFlux