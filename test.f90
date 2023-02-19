program test
  use patmo
  use patmo_commons
  use patmo_constants
  use patmo_parameters
  implicit none
  real*8::dt,x(speciesNumber),t,tend,imass
    real*8::heff(chemSpeciesNumber)
  integer::icell,i,j

  !init photochemistry
  call patmo_init()

  !load temperature and density profile
  call patmo_loadInitialProfile("profile.dat",unitH="km",unitX="1/cm3")

  !set BB flux (default Sun@1AU)
  call patmo_setFluxBB()
  call patmo_setGravity(9.8d2)
  call patmo_setEddyKzzAll(1.0d5)
  
  !read initial value 
  call patmo_dumpHydrostaticProfile("hydrostat.out")
  k_rain(:,:) = 0d0
    !calculate rainout
      !do j=1,chemSpeciesNumber
      call computeRainout(1,2.0d-3)   !OCS
      call computeRainout(12,5.0d-2) !CS2
   ! end do
  
  !get initial mass, g/cm3
  imass = patmo_getTotalMass()
 !! print *,"mass:",imass
  !dt = 1.0
  dt = secondsPerYear * 1d-5
  tend = secondsPerYear * 10d0
  t = 0d0

  !loop on time
  do
     dt = dt * 1.1
 	 !nall(1,patmo_idx_COS) = nall(1,patmo_idx_COS) + 600d0 * dt
     !nall(1,patmo_idx_CS2) = nall(1,patmo_idx_CS2) + 5d0 * dt
   	 
	call patmo_run(dt)
	
	!nall(1,patmo_idx_COS) = nall(1,patmo_idx_COS) - 0.01/1.0d5*nall(1,patmo_idx_COS)
    
	t = t + dt

     call patmo_dumpDensityToFile(35,t,patmo_idx_COS)
	 call patmo_dumpDensityToFile(36,t,patmo_idx_CS2)
	 call patmo_dumpDensityToFile(38,t,patmo_idx_CS)


	 write(37,*) t, krate(:,18)   !OCS!
    
	print '(F11.2,a2)',t/tend*1d2," %"
     if(t>=tend) exit
  end do

  !dump final hydrostatic equilibrium
  call patmo_dumpHydrostaticProfile("hydrostatEnd.out")

end program test

!**************
  !  subroutine read_Henry(fname)
  !    use patmo_commons
  !    use patmo_parameters
  !    implicit none

	  ! read file 
  !    open(22,file=trim(fname),status="old",iostat=ios)
          !check for file opening
  !        if(ios/=0) then
  !           print *,"ERROR: problem while opening ",trim(fname)
  !           stop
  !        end if
	!do j=1,chemSpeciesNumber
    ! 	read(22,*,iostat=ios) datar(j,2)
    !	heff(j)=datar(j,2)
	!end do
    
!	close(22)
!	end subroutine read_Henry
!**************
    subroutine computeRainout(i,heff)
      use patmo_commons
      use patmo_constants
      use patmo_parameters
      implicit none
      real*8::gamma(cellsNumber)  ! Precipation and Nonprecipitation time 
    !  real*8::gam15(cellsNumber)  ! Precipation and Nonprecipitation time
     ! real*8::gam8(cellsNumber)   ! Precipation and Nonprecipitation time
      real*8::wh2o(cellsNumber)   ! Rate of wet removal
      real*8::rkj(cellsNumber,chemSpeciesNumber)    ! Average Remeval Frequency	  
      real*8::y(cellsNumber),fz(cellsNumber),wl,qj(cellsNumber,chemSpeciesNumber)
	  real*8::heff !Henry's Constant
    !  real*8::k_rain(cellsNumber,chemSpeciesNumber)
      real*8::zkm(cellsNumber),temp(cellsNumber),gam15,gam8
    !  real*8::height(cellsNumber),TgasAll(cellsNumber)
	  integer::i,j
		   !Gas constant
   real*8,parameter::Rgas_atm = 1.36d-22 !Latm/K/mol
   !  real*8,parameter::Rgas_atm = 8.20574d-2 !Latm/K/mol	   !Gas constant

      gam15 = 8.64d5/2d0
      gam8 = 7.0d6/2d0
     !ZKM = 5.0
      wl = 1d0
     !HEFF = 1000 H2SO4
	 !TEMP = 200
	 !R = 1.36E-22
	 !AV = 6.02E+23
	 zkm(:) = height(:)/1d5
	 temp(:) = TgasAll(:)
	 k_rain(:,:) = 0d0

    do j=1,17

	   !FIND APPROPRIATE GAMMA
        if (zkm(j).LE.1.51d0) then
           gamma(j) = gam15
        else if (zkm(j).LT.8d0) then
           gamma(j) = gam15 + (gam8-gam15)*((zkm(j)-1.5d0)/6.5d0)
        else
           gamma(j) = gam8
        end if 

       !FIND WH2O
        if (zkm(j).LE.1d0) then
           y(j) = 11.35d0 + 1d-1*zkm(j)
        else
           y(j) = 11.5444d0 - 0.085333d0*zkm(j) - 9.1111d-03*zkm(j)*zkm(j)
        end if
        wh2o(j) = 10d0**y(j)

	   !FIND F(Z)
        if (zkm(j).LE.1.51d0) then
           fz(j) = 1d-1
        else
           fz(j) = 0.16615d0 - 0.04916d0*zkm(j) + 3.37451d-3*zkm(j)*zkm(j)
        end if 
		
	   !Raintout rates
        rkj(j,i) = wh2o(j)/55d0 /(av*wl*1.0d-9 + 1d0/(heff*Rgas_atm*temp(j)))
	    qj(j,i) = 1d0 - fz(j) + fz(j)/(gamma(j)*rkj(j,i)) * (1d0 - EXP(-rkj(j,i)*gamma(j)))
        k_rain(j,i) = (1d0 - EXP(-rkj(j,i)*gamma(j)))/(gamma(j)*qj(j,i))
  		!  write(*,*) k_rain(j,1)
		!  write(*,*) k_rain(j,3)
	end do
	
	
	
	 if (i==1) then
    !Output 
	do j=1,17
      open  (20,file="Rainout-OCS.txt")
      write (20,*) 'GAMMA', gamma(j)
      write (20,*) 'ZKM', zkm(j)
      write (20,*) 'WH2O', wh2o(j)
	  write (20,*) 'RKJ', rkj(j,i)
 	  write (20,*) 'QJ', qj(j,i)
	  write (20,*) 'K_RAIN',  k_rain(j,i)
     end do
	 end if

	 end subroutine computeRainout
	!**************