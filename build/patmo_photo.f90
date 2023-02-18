module patmo_photo
contains

  !**************
  subroutine loadPhotoMetric(fname)
    use patmo_parameters
    implicit none
    character(len=*),intent(in)::fname
    character::aout
    real*8::rout(4)
    integer::ios,i

    !open metric file (energy: mid/eV, span/eV, left/eV, right/eV)
    open(22,file=trim(fname),status="old",iostat=ios)
    if(ios/=0) then
      print *,"ERROR: problem loading "//trim(fname)
      stop
    end if

    read(22,'(a)') aout

    !loop on bins to read (already interpolated by python)
    do i=1,photoBinsNumber
      read(22,*,iostat=ios) rout(:)
      if(ios/=0) then
        print *,"ERROR: problem while reading "//trim(fname)
        stop
      end if
      !load energy mid, span, left, right (eV)
      energyMid(i) = rout(1)
      energySpan(i) = rout(2)
      energyLeft(i) = rout(3)
      energyRight(i) = rout(4)
    end do

    close(22)

  end subroutine loadPhotoMetric

  !***************
  subroutine loadAllPhotoXsecs()
    use patmo_parameters
    implicit none

    !36SO2 -> 36SO + O
    call loadPhotoXsec("xsecs/36SO2__36SO_O.dat",1)
    !36SO2 -> 36SO2_1 + O
    call loadPhotoXsec("xsecs/36SO2__36SO2_1_O.dat",2)
    !36SO2 -> 36SO2_3 + O
    call loadPhotoXsec("xsecs/36SO2__36SO2_3_O.dat",3)
    !36SO3 -> 36SO2 + O
    call loadPhotoXsec("xsecs/36SO3__36SO2_O.dat",4)
    !32COS -> CO + 32S
    call loadPhotoXsec("xsecs/32COS__CO_32S.dat",5)
    !34SO2 -> 34SO2_3 + O
    call loadPhotoXsec("xsecs/34SO2__34SO2_3_O.dat",6)
    !34SO2 -> 34SO + O
    call loadPhotoXsec("xsecs/34SO2__34SO_O.dat",7)
    !34SO2 -> 34SO2_1 + O
    call loadPhotoXsec("xsecs/34SO2__34SO2_1_O.dat",8)
    !33CS2 -> 33CS + 33S
    call loadPhotoXsec("xsecs/33CS2__33CS_33S.dat",9)
    !O3 -> O2 + O
    call loadPhotoXsec("xsecs/O3__O2_O.dat",10)
    !O2 -> O + O
    call loadPhotoXsec("xsecs/O2__O_O.dat",11)
    !34COS -> CO + 34S
    call loadPhotoXsec("xsecs/34COS__CO_34S.dat",12)
    !36H2S -> 36SH + H
    call loadPhotoXsec("xsecs/36H2S__36SH_H.dat",13)
    !33SO3 -> 33SO2 + O
    call loadPhotoXsec("xsecs/33SO3__33SO2_O.dat",14)
    !33H2S -> 33SH + H
    call loadPhotoXsec("xsecs/33H2S__33SH_H.dat",15)
    !32SO2 -> 32SO2_3 + O
    call loadPhotoXsec("xsecs/32SO2__32SO2_3_O.dat",16)
    !32SO2 -> 32SO + O
    call loadPhotoXsec("xsecs/32SO2__32SO_O.dat",17)
    !32SO2 -> 32SO2_1 + O
    call loadPhotoXsec("xsecs/32SO2__32SO2_1_O.dat",18)
    !34CS2 -> 34CS + 34S
    call loadPhotoXsec("xsecs/34CS2__34CS_34S.dat",19)
    !32SO3 -> 32SO2 + O
    call loadPhotoXsec("xsecs/32SO3__32SO2_O.dat",20)
    !34SO -> 34S + O
    call loadPhotoXsec("xsecs/34SO__34S_O.dat",21)
    !32SO -> 32S + O
    call loadPhotoXsec("xsecs/32SO__32S_O.dat",22)
    !33COS -> CO + 33S
    call loadPhotoXsec("xsecs/33COS__CO_33S.dat",23)
    !33SO2 -> 33SO2_3 + O
    call loadPhotoXsec("xsecs/33SO2__33SO2_3_O.dat",24)
    !33SO2 -> 33SO + O
    call loadPhotoXsec("xsecs/33SO2__33SO_O.dat",25)
    !33SO2 -> 33SO2_1 + O
    call loadPhotoXsec("xsecs/33SO2__33SO2_1_O.dat",26)
    !34H2S -> 34SH + H
    call loadPhotoXsec("xsecs/34H2S__34SH_H.dat",27)
    !32CS2 -> 32CS + 32S
    call loadPhotoXsec("xsecs/32CS2__32CS_32S.dat",28)
    !33SO -> 33S + O
    call loadPhotoXsec("xsecs/33SO__33S_O.dat",29)
    !36SO -> 36S + O
    call loadPhotoXsec("xsecs/36SO__36S_O.dat",30)
    !36CS2 -> 36CS + 36S
    call loadPhotoXsec("xsecs/36CS2__36CS_36S.dat",31)
    !32H2S -> 32SH + H
    call loadPhotoXsec("xsecs/32H2S__32SH_H.dat",32)
    !36COS -> CO + 36S
    call loadPhotoXsec("xsecs/36COS__CO_36S.dat",33)
    !34SO3 -> 34SO2 + O
    call loadPhotoXsec("xsecs/34SO3__34SO2_O.dat",34)

  end subroutine loadAllPhotoXsecs

  !***************
  subroutine loadPhotoXsec(fname,index)
    use patmo_commons
    use patmo_parameters
    implicit none
    character(len=*),intent(in)::fname
    character(len=100)::aout
    integer,intent(in)::index
    integer::ios,i
    real*8::rout(2)

    !open xsec file
    open(22,file=trim(fname),status="old",iostat=ios)
    if(ios/=0) then
      print *,"ERROR: problem loading "//trim(fname)
      stop
    end if

    !skip header (one line)
    read(22,'(a)') aout

    !loop on bins to read (already interpolated by python)
    do i=1,photoBinsNumber
      read(22,*,iostat=ios) rout(:)
      if(ios/=0) then
        print *,"ERROR: problem while reading "//trim(fname)
        stop
      end if
      !load xsecs for all cells (cm2)
      xsecAll(i,index) = rout(2)
    end do

    close(22)
  end subroutine loadPhotoXsec

  !**************************
  function fluxBB(x,Tbb)
    use patmo_constants
    implicit none
    real*8,intent(in)::x,Tbb
    real*8::fluxBB,xexp

    !exponent
    xexp = x/kboltzmann_eV/Tbb

    !default value
    fluxBB = 0d0

    !limit exp overflow
    if(xexp<3d2.and.x>1d-10) then
      fluxBB = 2d0*x**3/planck_eV**2/clight**2 &
          / (exp(xexp)-1d0)
    end if

  end function fluxBB

end module patmo_photo
