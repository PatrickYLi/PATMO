module patmo
contains

  !*************
  !initialize all
  subroutine patmo_init()
    use patmo_photo
    use patmo_parameters
    use patmo_utils

    !load photo cross-sections
    call loadAllPhotoXsecs()
    !load photo metrics (i.e. binning)
    call loadPhotoMetric("xsecs/photoMetric.dat")

    !init default photoflux
    photoFlux(:) = 0d0

    !init cumulative flux for entropy production
    ! (i.e. integrated with time)
    cumulativeFlux(:,:) = 0d0

    !set reactions rate to zero by default
    krate(:,:) = 0d0

    !load verbatim reactions
    call loadReactionsVerbatim()

  end subroutine patmo_init

  !**************
  !run model for a time-step
  subroutine patmo_run(dt)
    !use dvode_f90_m
    use patmo_parameters
    use patmo_commons
    use patmo_ode
    use patmo_jacobian
    use patmo_sparsity
    use patmo_rates
    use patmo_photoRates
    use patmo_reverseRates
    use patmo_utils
    implicit none
    real*8,intent(in)::dt
    real*8::atol(neqAll),rtol(neqAll)
    real*8::tstart,tend,n(neqAll)
    real*8::sumxn(photoBinsNumber),m(speciesNumber)
    !type(VODE_OPTS)::OPTIONS
    integer::istate,itask,i,j

    integer,parameter::meth=2
    integer,parameter::lwm=2*neqAll**2 + 2*neqAll &
         + (neqAll**2+10*neqAll)/2
    integer::itol,iopt,mf,lrw,liw
    integer::iwork(20+9*neqAll+LWM),neq(1)
    real*8::rwork(20+neqAll*6+3*neqAll+lwm)

    lrw = size(rwork)
    liw = size(iwork)

    iwork(:) = 0
    rwork(:) = 0d0

    atol(:) = 1d-10 !absolute tolerances (array)
    rtol(:) = 1d-4 !relative tolerances (array)

    !computes sparsity if not already done
    if(nonZeroElements==0) then
       call computeSparsity()
    end if

    itol = 4
    istate = 1
    itask = 1
    iopt = 0
    MF = 222
    call xsetf(0)

    !set solver options (DVODE_f90)
    !OPTIONS = SET_OPTS(SPARSE_J=.true., ABSERR_VECTOR=ATOL(:), &
    !     RELERR_VECTOR=RTOL(:), MXSTEP=100000, &
    !     USER_SUPPLIED_SPARSITY = .true., &
    !     MA28_RPS = .true., &
    !     USER_SUPPLIED_JACOBIAN = .false.)

    !set the sparsity structure (DVODE_f90)
    !CALL USERSETS_IAJA(iaSparsity, size(iaSparsity), &
    !     jaSparsity, size(jaSparsity))

    tstart = 0d0
    tend = dt

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
    end do

    !unroll chemistry
    do i=1,speciesNumber
       n((i-1)*cellsNumber+1:(i*cellsNumber)) = nall(:,i)
    end do
    !unroll Tgas
    n((positionTgas-1)*cellsNumber+1:(positionTgas*cellsNumber)) &
         = TgasAll(:)

    !compute rates
    call computeRates(TgasAll(:))
    call computePhotoRates(tauAll(:,:))
    call computeReverseRates(TgasAll(:))

    !compute tot density
    ntotAll(:) = sum(nall(:,1:chemSpeciesNumber),2)

    !compute mean molecular mass of the whole atmosphere
    ! (averaged between layers)
    m(:) = getSpeciesMass()
    meanMolecularMass = 0d0
    do i=1,cellsNumber
       meanMolecularMass = meanMolecularMass &
            + sum(m(1:chemSpeciesNumber) &
            * nAll(i,1:chemSpeciesNumber)) &
            / ntotAll(i) / cellsNumber
    end do

    !call the solver (DVODE_f90)
    !CALL DVODE_F90(fex, &
    !     neqAll, n(:), &
    !     tstart, tend, ITASK, ISTATE, OPTIONS, &
    !     jex)

    neq(:) = neqAll

    !loop until istate=2 or istate=error
    do
       CALL DLSODES(fex, NEQ(:), n(:), tstart, dt, ITOL, RTOL, ATOL,&
            ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JES, MF)
       !recompute sparsity if required
       if(istate==-5.or.istate==-4) then
          istate = 3
          cycle
       end if
       !loop when max iteration reached
       if(istate/=-1) exit
       istate = 1
    end do

    !check output state
    if(istate/=2) then
       print *,"ERROR: istate=",istate
       stop
    end if

    !avoid negative species
    do i=1,neqAll
       n(i) = max(n(i),0d0)
    end do

    !roll chemistry
    do i=1,speciesNumber
       nall(:,i) = n((i-1)*cellsNumber+1:(i*cellsNumber))
    end do
    !roll Tgas
    TgasAll(:) = n((positionTgas-1)*cellsNumber+1:(positionTgas*cellsNumber))

  end subroutine patmo_run

  !****************
  !dump the histogram of the connection degree to ifile, for the
  ! cell icell, at a given time (or any other independent varibles)
  subroutine patmo_dumpWeightedDegreeHistogram(ifile,icell,time)
    use patmo_utils
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::time
    integer,intent(in)::ifile,icell
    integer::hist(speciesNumber),i

    !get histogram
    hist(:) = getDegreeHistogram(nAll(:,:),icell)
    !write histogram to file
    do i=1,speciesNumber
       write(ifile,*) time, i, hist(i)
    end do
    write(ifile,*)

  end subroutine patmo_dumpWeightedDegreeHistogram

  !**************
  subroutine patmo_printBestFluxes(icell,bestFluxesNumber)
    use patmo_commons
    use patmo_utils
    use patmo_parameters
    implicit none
    integer,intent(in)::bestFluxesNumber,icell
    integer::idx(bestFluxesNumber),i
    real*8::flux(reactionsNumber)

    !get fluxes
    flux(:) = getFlux(nAll(:,:),icell)

    idx(:) = getBestFluxIdx(icell,bestFluxesNumber)
    print *,"*************"
    do i=1,bestFluxesNumber
       print *,idx(i),trim(reactionsVerbatim(idx(i))),flux(idx(i))
    end do

  end subroutine patmo_printBestFluxes

  !**************
  !compute cumulative flux for entropy production
  subroutine patmo_computeEntropyProductionFlux(dt)
    use patmo_utils
    implicit none
    real*8,intent(in)::dt

    call computeEntropyProductionFlux(dt)

  end subroutine patmo_computeEntropyProductionFlux

  !***************
  function patmo_getEntropyProduction(timeInterval)
    use patmo_utils
    implicit none
    real*8,intent(in)::timeInterval
    real*8::patmo_getEntropyProduction

    patmo_getEntropyProduction = getEntropyProduction(timeInterval)

  end function patmo_getEntropyProduction

  !***************
  !Assume a black-body flux, with starTbb (K), starRadius (Rsun)
  ! starDistance (AU). Default is Sun at 1AU
  subroutine patmo_setFluxBB(starTbb,starRadius,starDistance)
    use patmo_commons
    use patmo_parameters
    use patmo_constants
    use patmo_photo
    implicit none
    real*8,optional,intent(in)::starTbb,starRadius,starDistance
    real*8,parameter::AU2cm=1.496d13 !AU->cm
    real*8,parameter::Rsun2cm=6.963d10 !Rsun->cm
    real*8::Tbb,rstar,dstar
    integer::i

    !default is Sun
    Tbb = 5.777d3 !K
    rstar = Rsun2cm !cm
    dstar = AU2cm !cm

    !check optional parameters
    if(present(starTbb)) Tbb = starTbb
    if(present(starRadius)) rstar = starRadius*Rsun2cm
    if(present(starTbb)) dstar = starDistance*AU2cm

    !integrate flux
    do i=1,photoBinsNumber
       photoFlux(i) = (fluxBB(energyLeft(i),Tbb)+fluxBB(energyRight(i),Tbb)) &
            * energySpan(i)/2d0
    end do

    !scale geometric flux
    photoFlux(:) =  pi*rstar**2/dstar**2 * photoFlux(:)
		open(58,file="solar_flux.txt",status="old")
		do i=1,photoBinsNumber
			read(58,*) photoFlux(i)
		end do
	close(58)

  end subroutine patmo_setFluxBB

  !***************
  !dump opacity to file fname using unitEenergy as unit for
  ! energy (eV or mircron, eV default).
  ! File format is energy,layer,opacity
  subroutine patmo_dumpOpacity(fname,unitEnergy)
    use patmo_commons
    use patmo_constants
    use patmo_parameters
    implicit none
    character(len=*),optional,intent(in)::unitEnergy
    character(len=*),intent(in)::fname
    character(len=100)::unitE
    integer::i,j

    unitE = "eV"
    if(present(unitEnergy)) unitE = trim(unitEnergy)

    open(22,file=trim(fname),status="replace")
    if(trim(unitE)=="eV") then
       !loop on energy
       do i=1,photoBinsNumber
          !loop on cells
          do j=1,cellsNumber
             write(22,*) energyMid(i),height(j)/1d5,tauAll(i,j)
          end do
          write(22,*)
       end do
    else if(trim(unitE)=="micron") then
       !loop on energy
       do i=photoBinsNumber,1,-1
          !loop on cells
          do j=1,cellsNumber
             !h*c/E -> cm -> micron
             write(22,*) 1d4*planck_eV*clight/energyMid(i),height(j)/1d5,&
                  tauAll(i,j)
          end do
          write(22,*)
       end do
    else
       print *,"ERROR: unknown unit "//trim(unitEnergy)
       stop
    end if
    close(22)

  end subroutine patmo_dumpOpacity

  !***************
  !find hydrostatic equilbrium knowing the pressure at ground
  ! (pground), using dp/dz = -mu*p*g/k/T
  ! Pressure unit is defined in unitP, default dyne
  subroutine patmo_hydrostaticEquilibrium(pground,unitP)
    use patmo_commons
    use patmo_constants
    use patmo_parameters
    use patmo_utils
    implicit none
    real*8,intent(in)::pground
    real*8::p,zold,dz,ntot,n(speciesNumber)
    integer::i
    character(len=*),optional,intent(in)::unitP
    character(len=50)::units

    !optional argument units
    if(present(unitP)) then
       units = trim(unitP)
    else
       !default
       units = "dyne"
    end if

    !convert initial pressure to dyne/cm2 if necessary
    if(trim(units)=="dyne") then
       p = pground
    elseif(trim(units)=="dyne/cm2") then
       p = pground
    elseif(trim(units)=="atm") then
       p = pground*1.013250d6
    elseif(trim(units)=="mbar") then
       p = pground*1d3
    elseif(trim(units)=="bar") then
       p = pground*1d6
    else
       !error if units unknonw
       print *,"ERROR: unknown pressure unit for hydrostatic eq ",trim(units)
       stop
    end if

    !initial conditions
    zold = 0d0
    !loop on cells
    do i=1,cellsNumber
       !difference in height
       dz = height(i)-zold
       !temp array
       n(:) = nall(i,:)
       !compute difference in pressure
       p = p - getMeanMass(n(:)) * p / kboltzmann &
            / TgasAll(i) * gravity * dz
       !total number density p=n*k*T
       ntot = p/kboltzmann/TgasAll(i)
       !resacale abundances depending on total pressure
       nall(i,1:chemSpeciesNumber) = &
            n(1:chemSpeciesNumber) &
            / sum(n(1:chemSpeciesNumber))*ntot
       !store old height
       zold = height(i)
    end do

  end subroutine patmo_hydrostaticEquilibrium

  !***************
  !dump hydrostatic pressure profile to fname.
  ! Format: h/km, p/mbar, Tgas/K
  subroutine patmo_dumpHydrostaticProfile(fname)
    use patmo_commons
    use patmo_parameters
    use patmo_constants
    implicit none
    integer::i
    character(len=*),intent(in)::fname
    real*8::ntot

    open(22,file=trim(fname),status="replace")
    write(22,*) "#hydrostatic equilibrium dump"
    write(22,*) "#alt/km p/mbar Tgas/K"
    !loop on cells
    do i=1,cellsNumber
       ntot = sum(nall(i,1:chemSpeciesNumber))
       write(22,*) height(i)/1d5,ntot*kboltzmann*TgasAll(i)/1d3,TgasAll(i)
    end do
    close(33)
    print *,"Hydrostatic equilibrium dumped in ",trim(fname)

  end subroutine patmo_dumpHydrostaticProfile

  !****************
  !load initial profile (density, etc...) from fname.
  ! Height in unitH, species in unitX
  subroutine patmo_loadInitialProfile(fname,unitH,unitX,defaultDensity)
    use patmo_parameters
    use patmo_commons
    use patmo_utils
    implicit none
    character(len=*),intent(in)::fname
    character(len=*),optional,intent(in)::unitH,unitX
    character(len=50)::units,unitsX
    character(len=50),allocatable::arow(:)
    real*8,optional,intent(in)::defaultDensity
    real*8,allocatable::x(:),rout(:)
    real*8::zold,defaultN
    integer::ios,i,idx,j,nonZero,offset
    logical::firstRow

    units = "cm"
    !optional argument units height
    if(present(unitH)) units = trim(unitH)

    unitsX = "1/cm3"
    !optional argument units chemical species
    if(present(unitX)) unitsX = trim(unitX)

    defaultN = 0d0
    !optional argument for default density
    if(present(defaultDensity)) defaultN = defaultDensity

    !read file
    print *,"reading ",trim(fname)
    open(22,file=trim(fname),status="old",iostat=ios)
    !check for file opening
    if(ios/=0) then
       print *,"ERROR: problem while opening ",trim(fname)
       stop
    end if

    !read until comment is found
    do
       read(22,*,iostat=ios) offset,nonZero
       if(ios==0) exit
    end do
    allocate(arow(offset+nonZero))
    allocate(x(nonZero))
    allocate(rout(offset))
    read(22,*) arow(:)

    !set default abundance
    nall(:,1:chemSpeciesNumber) = defaultN
    !set not chemial species to zero
    nall(:,chemSpeciesNumber+1:speciesNumber) = 0d0
    !loop on cells (file lines have to be the same number)
    do j=1,cellsNumber
       !read data+chemistry
       read(22,*,iostat=ios) rout(:),x(:)
       if(ios/=0) then
          print *,"ERROR: problem while reading ",trim(fname)
          if(j>1) print *,&
               "(could be less file lines than declared lines number)"
          stop
       end if

       !loop on data accoding to header
       do i=1,offset
          if(trim(arow(i))=="alt") then
             height(j) = rout(i)
          elseif(trim(arow(i))=="Tgas") then
             TgasAll(j) = rout(i)
          elseif(trim(arow(i))=="Dzz") then
             diffusionDzz(j) = rout(i)
          elseif(trim(arow(i))=="Kzz") then
             eddyKzz(j) = rout(i)
          elseif(trim(arow(i))=="index") then
             continue
          elseif(trim(arow(i))=="idx") then
             continue
          elseif(trim(arow(i))=="dummy") then
             continue
          else
             print *,"ERROR: unknown header element: ",trim(arow(i))
             stop
          end if
       end do

       !load species into common array
       do i=1,nonZero
          idx = getSpeciesIndex(arow(offset+i),error=.false.)
          if(idx/=-1) nall(j,idx) = x(i)
       end do

       !convert units if necessary
       if(trim(unitsX)=="ppbv") then
          nall(j,:) = nall(j,:)/sum(nall(j,1:chemSpeciesNumber))
       elseif(trim(unitsX)=="1/cm3") then
          continue
       else
          print *,&
               "ERROR: unknown chemical abundance units while reading profile",&
               trim(fname),trim(units)
          stop
       end if

    end do
    close(22)
    deallocate(arow)
    deallocate(x)
    deallocate(rout)

    !convert units if necessary
    if(trim(units)=="km") then
       height(:) = height(:)*1d5
    elseif(trim(units)=="cm") then
       continue
    else
       print *,"ERROR: unknown units while reading profile", &
            trim(fname),trim(units)
       stop
    end if

    !store inverse grid space squared, 1/dz**2, and dz
    zold = 0d0
    do j=1,cellsNumber
       idh2(j) = 1d0/(height(j)-zold)**2
       gridSpace(j) = (height(j)-zold)
       zold = height(j)
    end do

  end subroutine patmo_loadInitialProfile

  !****************
  !return total mass in g/cm3
  function patmo_getTotalMass()
    use patmo_commons
    use patmo_parameters
    use patmo_utils
    implicit none
    integer::icell
    real*8::patmo_getTotalMass
    real*8::m(speciesNumber)

    m(:) = getSpeciesMass()

    patmo_getTotalMass = 0d0
    do icell=1,cellsNumber
       patmo_getTotalMass = patmo_getTotalMass &
            + sum(m(1:chemSpeciesNumber) &
            * nall(icell,1:chemSpeciesNumber))
    end do

  end function patmo_getTotalMass

!***************************
function patmo_getTotalMassNuclei_C()
 use patmo_utils
 implicit none
 real*8::patmo_getTotalMassNuclei_C

patmo_getTotalMassNuclei_C = getTotalMassNuclei_C() 

end function

!***************************
function patmo_getTotalMassNuclei_H()
 use patmo_utils
 implicit none
 real*8::patmo_getTotalMassNuclei_H

patmo_getTotalMassNuclei_H = getTotalMassNuclei_H() 

end function

!***************************
function patmo_getTotalMassNuclei_O()
 use patmo_utils
 implicit none
 real*8::patmo_getTotalMassNuclei_O

patmo_getTotalMassNuclei_O = getTotalMassNuclei_O() 

end function

!***************************
function patmo_getTotalMassNuclei_N()
 use patmo_utils
 implicit none
 real*8::patmo_getTotalMassNuclei_N

patmo_getTotalMassNuclei_N = getTotalMassNuclei_N() 

end function

!***************************
function patmo_getTotalMassNuclei_S()
 use patmo_utils
 implicit none
 real*8::patmo_getTotalMassNuclei_S

patmo_getTotalMassNuclei_S = getTotalMassNuclei_S() 

end function

!***************************
function patmo_getTotalMassNuclei_CS2E()
 use patmo_utils
 implicit none
 real*8::patmo_getTotalMassNuclei_CS2E

patmo_getTotalMassNuclei_CS2E = getTotalMassNuclei_CS2E() 

end function



  !***************
  !set uniform grid spacing, cm
  subroutine patmo_setGridSpacing(dz)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::dz
    real*8::zold
    integer::j

    gridSpace(:) = dz
    !store inverse grid space squared, 1/dz**2, and height
    zold = 0d0
    do j=1,cellsNumber
       idh2(j) = 1d0/gridSpace(j)**2
       height(j) = zold
       zold = zold + gridSpace(j)
    end do

  end subroutine patmo_setGridSpacing

  !***************
  !set thermal diffusion
  subroutine patmo_setThermalDiffusion(alpha)
    use patmo_parameters
    implicit none
    real*8,intent(in)::alpha

    thermalDiffusionFactor = alpha

  end subroutine patmo_setThermalDiffusion

  !***************
  !set eddy Kzz coefficient of icell layer
  subroutine patmo_setEddyKzz(icell,kzz)
    use patmo_parameters
    implicit none
    real*8,intent(in)::kzz
    integer,intent(in)::icell

    eddyKzz(icell) = kzz

  end subroutine patmo_setEddyKzz

  !***************
  !set eddy Kzz, same for all layers
  subroutine patmo_setEddyKzzAll(kzz)
    use patmo_parameters
    implicit none
    real*8,intent(in)::kzz

    eddyKzz(:) = kzz

  end subroutine patmo_setEddyKzzAll

  !***************
  !set diffusion Dzz for layer icell
  subroutine patmo_setDiffusionDzz(icell,dzz)
    use patmo_parameters
    implicit none
    real*8,intent(in)::dzz
    integer,intent(in)::icell

    diffusionDzz(icell) = dzz

  end subroutine patmo_setDiffusionDzz

  !***************
  !set diffusion Dzz, same for all layers
  subroutine patmo_setDiffusionDzzAll(dzz)
    use patmo_parameters
    implicit none
    real*8,intent(in)::dzz

    diffusionDzz(:) = dzz

  end subroutine patmo_setDiffusionDzzAll

  !***************
  !append density of species idx to file number ifile
  subroutine patmo_dumpDensityToFile(ifile,time,idx)
    use patmo_commons
    use patmo_parameters
    implicit none
    integer,intent(in)::ifile,idx
    real*8,intent(in)::time
    integer::i,j

    do i=1,cellsNumber
       write(ifile,'(E17.8,I8,E17.8)') time, i, nall(i,idx)
    end do
    write(ifile,*)

  end subroutine patmo_dumpDensityToFile

  !****************
  !dump all mixing rations to file (one column one species)
  ! first column is layer number
  subroutine patmo_dumpAllMixingRatioToFile(fname)
    use patmo_commons
    use patmo_parameters
    use patmo_utils
    implicit none
    character(len=*),intent(in)::fname
    character(len=500)::header
    character(len=maxNameLength)::names(speciesNumber)
    integer::i

    names(:) = getSpeciesNames()
    !prepare header (species names)
    header = "#layer"
    do i=1,chemSpeciesNumber
       header = trim(header)//" "//names(i)
    end do

    !open file to dump mixing ratios
    open(67,file=trim(fname),status="replace")
    !write file header (species names)
    write(67,*) trim(header)
    !write mixing ratios
    do i=1,cellsNumber
       write(67,'(I5,999E17.8e3)') i,nall(i,1:chemSpeciesNumber)
    end do
    close(67)

  end subroutine patmo_dumpAllMixingRatioToFile

  !***************
  !append mixing ration of species idx to file number ifile
  subroutine patmo_dumpMixingRatioToFile(ifile,time,idx)
    use patmo_commons
    use patmo_parameters
    implicit none
    integer,intent(in)::ifile,idx
    real*8,intent(in)::time
    integer::i,j

    do i=1,cellsNumber
       write(ifile,'(E17.8,I8,E17.8)') time, i, nall(i,idx) &
            / sum(nall(i,1:chemSpeciesNumber))
    end do
    write(ifile,*)

  end subroutine patmo_dumpMixingRatioToFile

  !****************
  !set gravity in cm/s2
  subroutine patmo_setGravity(g)
    use patmo_parameters
    implicit none
    real*8,intent(in)::g

    gravity = g

  end subroutine patmo_setGravity

  !****************
  !set chemistry of layer icell
  subroutine patmo_setChemistry(icell,n)
    use patmo_commons
    use patmo_parameters
    implicit none
    integer,intent(in)::icell
    real*8,intent(in)::n(speciesNumber)

    nall(icell,:) = n(:)

  end subroutine patmo_setChemistry

  !****************
  !set the same chemistry for all the layers
  subroutine patmo_setChemistryAll(n)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::n(speciesNumber)
    integer::icell

    do icell=1,cellsNumber
       nall(icell,:) = n(:)
    end do

  end subroutine patmo_setChemistryAll

  !**************
  !set Tgas for layer icell
  subroutine patmo_setTgas(icell,Tgas)
    use patmo_commons
    use patmo_parameters
    implicit none
    integer,intent(in)::icell
    real*8,intent(in)::Tgas

    TgasAll(icell) = Tgas

  end subroutine patmo_setTgas

  !**************
  !set the same Tgas for all layers
  subroutine patmo_setTgasAll(Tgas)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::Tgas

    TgasAll(:) = Tgas

  end subroutine patmo_setTgasAll

  !**************
  !get density of species idx_species at layer icell
  function patmo_getDensity(icell,idx_species)
    use patmo_commons
    use patmo_parameters
    implicit none
    integer,intent(in)::icell,idx_species
    real*8::patmo_getDensity

    patmo_getDensity = nall(icell,idx_species)

  end function patmo_getDensity

  !**************
  !get Tgas of layer icell
  function patmo_getTgas(icell)
    use patmo_commons
    use patmo_parameters
    implicit none
    integer,intent(in)::icell
    real*8::patmo_getTgas

    patmo_getTgas = TgasAll(icell)

  end function patmo_getTgas

  !**************
!dump J-Values
  subroutine patmo_dumpJValue(fname)
   use patmo_commons
   use patmo_constants
   use patmo_parameters
   implicit none
   character(len=*),intent(in)::fname
   integer::i

   open(22,file=trim(fname),status="replace")
   write(22,*)"altitude/km, &
   5.CS2+OH->SCSOH,19.SO+O2->SO2+O,CS2->CS2E"
   !loop on cells
   do i=1,cellsNumber
       write(22,*) i, &
     krate(i,5), krate(i,19), krate(i,34)
      
   end do
   write(22,*)
   close(22)

 end subroutine patmo_dumpJValue

  !**************
 !dump all reaction rates
 subroutine patmo_dumpAllRates(fname)
   use patmo_commons
   use patmo_constants
   use patmo_parameters
   implicit none
   character(len=*),intent(in)::fname
   integer::i

   open(22,file=trim(fname),status="replace")
   write(22,*) "altitude/km, &
   1.CS2E+O2->CS2, 2.CS2E+N2->CS2, 3.CS2E+O2->CS+SO2, 4.CS2+OH->COS+SH, 5.CS2+OH->SCSOH, &
   6.SCSOH+O2->COS+HSO2, 7.CS2+O->CS+SO, 8.CS2+O->COS+S, 9.CS2+O->S2+CO, 10.CS+O->S+CO, &
   11.CS+O2->COS+O, 12.CS+O2->SO+CO, 13.CS+O3->COS+O2, 14.S+O2->SO+O, 15.S+O3->SO+O2, &
   16.S+OH->SO2+H, 17.S2+O->S+SO, 18.SO+O2->SO2+O, 19.SO+O3->SO2+O2, 20.SO+OH->SO2+H, &
   21.SH+O->SO+H, 22.SH+O2->SO+OH, 23.SH+O3->HSO+O2, 24.HSO+O2->SO2+OH, 25.HSO+O3->SO2+SH, &
   26.HSO2+O2->SO2+HO2, 27.COS+hv->CO+S, 28.SO2+hv->SO+O, 29.SO+hv->S+O, 30.CS2+hv->CS+S, &
   31.CS2+hv->CS2E "
   !loop on cells
   do i=1,cellsNumber
      write(22,*) i, &
      krate(i,1)*nall(i,patmo_idx_CS2E)*nall(i,patmo_idx_O2), &  !CS2E+O2->CS2 pesudo-1st order reaction
      krate(i,2)*nall(i,patmo_idx_CS2E)*nall(i,patmo_idx_N2), &  !CS2E+N2->CS2 pesudo-1st order reaction
      krate(i,3)*nall(i,patmo_idx_CS2E)*nall(i,patmo_idx_O2), &  !CS2E+O2->CS+SO2 
      krate(i,4)*nall(i,patmo_idx_CS2)*nall(i,patmo_idx_OH), &   !CS2+OH->COS+SH
      krate(i,5)*nall(i,patmo_idx_CS2)*nall(i,patmo_idx_OH), &   !CS2+OH->SCSOH

   !   krate(i,6)*nall(i,patmo_idx_SCSOH), &                      !SCSOH->CS2+OH
      krate(i,7)*nall(i,patmo_idx_SCSOH)*nall(i,patmo_idx_O2), & !SCSOH+O2->COS+HSO2 
      krate(i,8)*nall(i,patmo_idx_CS2)*nall(i,patmo_idx_O), &    !CS2+O->CS+SO
      krate(i,9)*nall(i,patmo_idx_CS2)*nall(i,patmo_idx_O), &    !CS2+O->COS+S
      krate(i,10)*nall(i,patmo_idx_CS2)*nall(i,patmo_idx_O), &   !CS2+O->S2+CO

      krate(i,11)*nall(i,patmo_idx_CS)*nall(i,patmo_idx_O), &    !CS+O->S+CO
      krate(i,12)*nall(i,patmo_idx_CS)*nall(i,patmo_idx_O2), &   !CS+O2->COS+O 
      krate(i,13)*nall(i,patmo_idx_CS)*nall(i,patmo_idx_O2), &   !CS+O2->SO+CO 
      krate(i,14)*nall(i,patmo_idx_CS)*nall(i,patmo_idx_O3), &   !CS+O3->COS+O2 
      krate(i,15)*nall(i,patmo_idx_S)*nall(i,patmo_idx_O2), &    !S+O2->SO+O 
 
      krate(i,16)*nall(i,patmo_idx_S)*nall(i,patmo_idx_O3), &    !S+O3->SO+O2 
      krate(i,17)*nall(i,patmo_idx_S)*nall(i,patmo_idx_OH), &    !S+OH->SO2+H
      krate(i,18)*nall(i,patmo_idx_S2)*nall(i,patmo_idx_O), &    !S2+O->S+SO
      krate(i,19)*nall(i,patmo_idx_SO)*nall(i,patmo_idx_O2), &   !SO+O2->SO2+O 
      krate(i,20)*nall(i,patmo_idx_SO)*nall(i,patmo_idx_O3), &   !SO+O3->SO2+O2 

      krate(i,21)*nall(i,patmo_idx_SO)*nall(i,patmo_idx_OH), &   !SO+OH->SO2+H
      krate(i,22)*nall(i,patmo_idx_SH)*nall(i,patmo_idx_O), &    !SH+O->SO+H
      krate(i,23)*nall(i,patmo_idx_SH)*nall(i,patmo_idx_O2), &   !SH+O2->SO+OH 
      krate(i,24)*nall(i,patmo_idx_SH)*nall(i,patmo_idx_O3), &   !SH+O3->HSO+O2 
      krate(i,25)*nall(i,patmo_idx_HSO)*nall(i,patmo_idx_O2), &  !HSO+O2->SO2+OH 

      krate(i,26)*nall(i,patmo_idx_HSO)*nall(i,patmo_idx_O3), &  !HSO+O3->SO2+SH 
      krate(i,27)*nall(i,patmo_idx_HSO2)*nall(i,patmo_idx_O2), & !HSO2+O2->SO2+HO2 
      krate(i,28)*nall(i,patmo_idx_COS), &                       !COS+hv->CO+S
      krate(i,29)*nall(i,patmo_idx_SO2), &                       !SO2+hv->SO+O
      krate(i,30)*nall(i,patmo_idx_SO), &                        !SO+hv->S+O

      krate(i,33)*nall(i,patmo_idx_CS2), &                       !CS2+hv->CS+S
      krate(i,34)*nall(i,patmo_idx_CS2)                          !CS2+hv->CS2E

   !   krate(i,1)*nall(i,patmo_idx_CS2E), &                       !CS2E+O2->CS2 pesudo-1st order reaction
   !   krate(i,2)*nall(i,patmo_idx_CS2E), &                       !CS2E+N2->CS2 pesudo-1st order reaction
   !   krate(i,3)*nall(i,patmo_idx_CS2E), &                       !CS2E+O2->CS+SO2 pesudo-1st order reaction
   !   krate(i,4)*nall(i,patmo_idx_CS2)*nall(i,patmo_idx_OH), &   !CS2+OH->COS+SH
   !   krate(i,5)*nall(i,patmo_idx_CS2)*nall(i,patmo_idx_OH), &   !CS2+OH->SCSOH
!
   !   krate(i,6)*nall(i,patmo_idx_SCSOH), &                      !SCSOH->CS2+OH
   !   krate(i,7)*nall(i,patmo_idx_SCSOH), &                      !SCSOH+O2->COS+HSO2 pesudo-1st order reaction
   !   krate(i,8)*nall(i,patmo_idx_CS2)*nall(i,patmo_idx_O), &    !CS2+O->CS+SO
   !   krate(i,9)*nall(i,patmo_idx_CS2)*nall(i,patmo_idx_O), &    !CS2+O->COS+S
   !   krate(i,10)*nall(i,patmo_idx_CS2)*nall(i,patmo_idx_O), &   !CS2+O->S2+CO
!
   !   krate(i,11)*nall(i,patmo_idx_CS)*nall(i,patmo_idx_O), &    !CS+O->S+CO
   !   krate(i,12)*nall(i,patmo_idx_CS), &                        !CS+O2->COS+O pesudo-1st order reaction
   !   krate(i,13)*nall(i,patmo_idx_CS), &                        !CS+O2->SO+CO pesudo-1st order reaction
   !   krate(i,14)*nall(i,patmo_idx_CS), &                        !CS+O3->COS+O2 pesudo-1st order reaction
   !   krate(i,15)*nall(i,patmo_idx_S), &                         !S+O2->SO+O pesudo-1st order reaction
 !
   !   krate(i,16)*nall(i,patmo_idx_S), &                         !S+O3->SO+O2 pesudo-1st order reaction
   !   krate(i,17)*nall(i,patmo_idx_S)*nall(i,patmo_idx_OH), &    !S+OH->SO2+H
   !   krate(i,18)*nall(i,patmo_idx_S2)*nall(i,patmo_idx_O), &    !S2+O->S+SO
   !   krate(i,19)*nall(i,patmo_idx_SO), &                        !SO+O2->SO2+O pesudo-1st order reaction
   !   krate(i,20)*nall(i,patmo_idx_SO), &                        !SO+O3->SO2+O2 pesudo-1st order reaction
!
   !   krate(i,21)*nall(i,patmo_idx_SO)*nall(i,patmo_idx_OH), &   !SO+OH->SO2+H
   !   krate(i,22)*nall(i,patmo_idx_SH)*nall(i,patmo_idx_O), &    !SH+O->SO+H
   !   krate(i,23)*nall(i,patmo_idx_SH), &                        !SH+O2->SO+OH pesudo-1st order reaction
   !   krate(i,24)*nall(i,patmo_idx_SH), &                        !SH+O3->HSO+O2 pesudo-1st order reaction
   !   krate(i,25)*nall(i,patmo_idx_HSO), &                       !HSO+O2->SO2+OH pesudo-1st order reaction
!
   !   krate(i,26)*nall(i,patmo_idx_HSO), &                       !HSO+O3->SO2+SH pesudo-1st order reaction
   !   krate(i,27)*nall(i,patmo_idx_HSO2), &                      !HSO2+O2->SO2+HO2 pesudo-1st order reaction
   !   krate(i,28)*nall(i,patmo_idx_COS), &                       !COS+hv->CO+S
   !   krate(i,29)*nall(i,patmo_idx_SO2), &                       !SO2+hv->SO+O
   !   krate(i,30)*nall(i,patmo_idx_SO), &                        !SO+hv->S+O
!
   !   krate(i,31)*nall(i,patmo_idx_CS2), &                       !CS2+hv->CS+S
   !   krate(i,32)*nall(i,patmo_idx_CS2)                          !CS2+hv->CS2E

   end do
   write(22,*)
   close(22)

end subroutine patmo_dumpAllRates

!**************
 !dump all reverse reaction rates
subroutine patmo_dumpAllReverseRates(fname)
   use patmo_commons
   use patmo_constants
   use patmo_parameters
   implicit none
   character(len=*),intent(in)::fname
   integer::i

   open(22,file=trim(fname),status="replace")
   write(22,*) "altitude/km, &
   1.CS2+OH->COS+SH,    2.CS2+OH->SCSOH,     3.SCSOH+O2->COS+HSO2,   4.CS2+O->CS+SO,   5.CS2+O->COS+S, &
   6.CS2+O->S2+CO,      7.CS+O->S+CO,        8.CS+O2->COS+O,         9.CS+O2->SO+CO,   10.CS+O3->COS+O2,  &
   11.S+O2->SO+O,       12.S+O3->SO+O2,      13.S+OH->SO2+H,         14.S2+O->S+SO,    15.SO+O2->SO2+O, &
   16.SO+O3->SO2+O2,    17.SO+OH->SO2+H,     18.SH+O->SO+H,          19.SH+O2->SO+OH,20.SH+O3->HSO+O2,&
   21.HSO+O2->SO2+OH,   22.HSO+O3->SO2+SH,   23.HSO2+O2->SO2+HO2"
   
   !loop on cells
   do i=1,cellsNumber
      write(22,*) i, &
      !krate(i,)*nall(i,patmo_idx_)*nall(i,patmo_idx_), &  !CS2E+O2->CS2 pesudo-1st order reaction
      !krate(i,)*nall(i,patmo_idx_)*nall(i,patmo_idx_), &  !CS2E+N2->CS2 pesudo-1st order reaction
      !krate(i,)*nall(i,patmo_idx_)*nall(i,patmo_idx_), &  !CS2E+O2->CS+SO2 pesudo-1st order reaction
      krate(i,38)*nall(i,patmo_idx_COS)*nall(i,patmo_idx_SH), &   !CS2+OH->COS+SH
      krate(i,39)*nall(i,patmo_idx_SCSOH), &   !CS2+OH->SCSOH

   !   krate(i,6)*nall(i,patmo_idx_SCSOH), &                      !SCSOH->CS2+OH
      krate(i,41)*nall(i,patmo_idx_COS)*nall(i,patmo_idx_SCSOH), &    !SCSOH+O2->COS+HSO2 pesudo-1st order reaction
      krate(i,42)*nall(i,patmo_idx_CS)*nall(i,patmo_idx_SO), &!CS2+O->CS+SO
      krate(i,43)*nall(i,patmo_idx_COS)*nall(i,patmo_idx_S), &!CS2+O->COS+S
      krate(i,44)*nall(i,patmo_idx_S2)*nall(i,patmo_idx_O), &!CS2+O->S2+CO

      krate(i,45)*nall(i,patmo_idx_S)*nall(i,patmo_idx_CO), &!CS+O->S+CO
      krate(i,46)*nall(i,patmo_idx_COS)*nall(i,patmo_idx_O), &!CS+O2->COS+O pesudo-1st order reaction
      krate(i,47)*nall(i,patmo_idx_SO)*nall(i,patmo_idx_CO), &!CS+O2->SO+CO pesudo-1st order reaction
      krate(i,48)*nall(i,patmo_idx_COS)*nall(i,patmo_idx_O2), &!CS+O3->COS+O2 pesudo-1st order reaction
      krate(i,49)*nall(i,patmo_idx_SO)*nall(i,patmo_idx_O), &!S+O2->SO+O pesudo-1st order reaction
 
      krate(i,50)*nall(i,patmo_idx_SO)*nall(i,patmo_idx_O2), &!S+O3->SO+O2 pesudo-1st order reaction
      krate(i,51)*nall(i,patmo_idx_SO2)*nall(i,patmo_idx_H), &!S+OH->SO2+H
      krate(i,52)*nall(i,patmo_idx_S)*nall(i,patmo_idx_SO), &!S2+O->S+SO
      krate(i,53)*nall(i,patmo_idx_SO2)*nall(i,patmo_idx_O), &!SO+O2->SO2+O pesudo-1st order reaction
      krate(i,54)*nall(i,patmo_idx_SO2)*nall(i,patmo_idx_O2), &!SO+O3->SO2+O2 pesudo-1st order reaction

      krate(i,55)*nall(i,patmo_idx_SO2)*nall(i,patmo_idx_H), &!SO+OH->SO2+H
      krate(i,56)*nall(i,patmo_idx_SO)*nall(i,patmo_idx_H), &!SH+O->SO+H
      krate(i,57)*nall(i,patmo_idx_SO)*nall(i,patmo_idx_OH), &!SH+O2->SO+OH pesudo-1st order reaction
      krate(i,58)*nall(i,patmo_idx_HSO)*nall(i,patmo_idx_O2), &!SH+O3->HSO+O2 pesudo-1st order reaction
      krate(i,59)*nall(i,patmo_idx_SO2)*nall(i,patmo_idx_OH), &!HSO+O2->SO2+OH pesudo-1st order reaction

      krate(i,60)*nall(i,patmo_idx_SO2)*nall(i,patmo_idx_SH), &!HSO+O3->SO2+SH pesudo-1st order reaction
      krate(i,61)*nall(i,patmo_idx_SO2)*nall(i,patmo_idx_HO2)!HSO2+O2->SO2+HO2 pesudo-1st order reaction
      
      !photochemistry no reverse reaction
      !krate(i,)*nall(i,patmo_idx_)*nall(i,patmo_idx_), &!COS+hv->CO+S
      !krate(i,)*nall(i,patmo_idx_)*nall(i,patmo_idx_), &!SO2+hv->SO+O
      !krate(i,)*nall(i,patmo_idx_)*nall(i,patmo_idx_), &!SO+hv->S+O

      !krate(i,)*nall(i,patmo_idx_)*nall(i,patmo_idx_), &!CS2+hv->CS+S
      !krate(i,)*nall(i,patmo_idx_)*nall(i,patmo_idx_), &!CS2+hv->CS2E
     


   end do
   write(22,*)
   close(22)

end subroutine patmo_dumpAllReverseRates

end module patmo
