module patmo_volc
  implicit none
  private

  type volcano_event
     real*8::startTime
     real*8::duration
     real*8::center
     real*8::sigma
     real*8::so2ColumnFlux
     real*8::ashTau550
     real*8::ashLifetime
     real*8::ashSettling
     real*8::ashWavelengthExp
  end type volcano_event

  logical::volcanoEnabled = .false.
  logical::warnedMissingSO2 = .false.
  integer::volcanoEventsNumber = 0
  integer::so2IndexCache = -2
  real*8::volcanoElapsedTime = 0d0
  type(volcano_event),allocatable::events(:)

  public::patmo_volc_loadEvents
  public::patmo_volc_reset
  public::patmo_volc_setEnabled
  public::patmo_volc_setTime
  public::patmo_volc_getTime
  public::patmo_volc_advanceTime
  public::patmo_volc_injectSO2Column
  public::patmo_volc_addSources
  public::patmo_volc_applyAshOpacity
  public::patmo_volc_dumpState

contains

  !***************
  subroutine patmo_volc_reset()
    implicit none

    if(allocated(events)) deallocate(events)
    volcanoEnabled = .false.
    warnedMissingSO2 = .false.
    volcanoEventsNumber = 0
    so2IndexCache = -2
    volcanoElapsedTime = 0d0

  end subroutine patmo_volc_reset

  !***************
  subroutine patmo_volc_setEnabled(enabled)
    implicit none
    logical,intent(in)::enabled

    volcanoEnabled = enabled

  end subroutine patmo_volc_setEnabled

  !***************
  subroutine patmo_volc_setTime(timeSeconds)
    implicit none
    real*8,intent(in)::timeSeconds

    volcanoElapsedTime = max(timeSeconds,0d0)

  end subroutine patmo_volc_setTime

  !***************
  function patmo_volc_getTime()
    implicit none
    real*8::patmo_volc_getTime

    patmo_volc_getTime = volcanoElapsedTime

  end function patmo_volc_getTime

  !***************
  subroutine patmo_volc_advanceTime(dt)
    implicit none
    real*8,intent(in)::dt

    volcanoElapsedTime = max(volcanoElapsedTime + dt,0d0)

  end subroutine patmo_volc_advanceTime

  !***************
  ! Preferred event file format, one event per line:
  ! event_id=name start_day=... duration_day=... plume_center_km=...
  ! plume_sigma_km=... so2_flux_cm2_s=... ash_tau_550=...
  ! ash_lifetime_day=... ash_settling_cm_s=... ash_lambda_exponent=...
  ! Legacy numeric rows are still accepted in the same column order.
  subroutine patmo_volc_loadEvents(fname)
    use patmo_constants
    implicit none
    character(len=*),intent(in)::fname
    character(len=512)::line
    integer::ios,unitEvent,i,commentPos
    real*8::startDay,durationDay,centerKm,sigmaKm
    real*8::so2ColumnFlux,ashTau550,ashLifetimeDay
    real*8::ashSettling,ashWavelengthExp
    logical::needsSO2

    call patmo_volc_reset()

    unitEvent = 91
    open(unitEvent,file=trim(fname),status="old",iostat=ios)
    if(ios/=0) then
       print *,"ERROR: problem while opening volcano event file ",trim(fname)
       stop
    end if

    volcanoEventsNumber = 0
    do
       read(unitEvent,'(A)',iostat=ios) line
       if(ios/=0) exit
       commentPos = patmo_volc_commentStart(line)
       if(commentPos>0) line = line(:commentPos-1)
       if(len_trim(line)==0) cycle
       volcanoEventsNumber = volcanoEventsNumber + 1
    end do
    close(unitEvent)

    if(volcanoEventsNumber<=0) then
       print *,"WARNING: volcano event file has no active event: ",trim(fname)
       return
    end if

    allocate(events(volcanoEventsNumber))
    open(unitEvent,file=trim(fname),status="old",iostat=ios)
    if(ios/=0) then
       print *,"ERROR: problem while reopening volcano event file ",trim(fname)
       stop
    end if

    i = 0
    do
       read(unitEvent,'(A)',iostat=ios) line
       if(ios/=0) exit
       commentPos = patmo_volc_commentStart(line)
       if(commentPos>0) line = line(:commentPos-1)
       if(len_trim(line)==0) cycle

       startDay = 0d0
       durationDay = 0d0
       centerKm = 0d0
       sigmaKm = 1d0
       so2ColumnFlux = 0d0
       ashTau550 = 0d0
       ashLifetimeDay = 0d0
       ashSettling = 0d0
       ashWavelengthExp = 0d0

       if(index(line,"=")>0) then
          call patmo_volc_parseKeywordEvent(line,startDay,durationDay, &
               centerKm,sigmaKm,so2ColumnFlux,ashTau550,ashLifetimeDay, &
               ashSettling,ashWavelengthExp,ios)
       else
          read(line,*,iostat=ios) startDay,durationDay,centerKm,sigmaKm, &
               so2ColumnFlux,ashTau550,ashLifetimeDay,ashSettling, &
               ashWavelengthExp
       end if
       if(ios/=0) then
          print *,"ERROR: malformed volcano event row:"
          print *,trim(line)
          stop
       end if

       i = i + 1
       events(i)%startTime = startDay * secondsPerDay
       events(i)%duration = max(durationDay,0d0) * secondsPerDay
       events(i)%center = centerKm * 1d5
       events(i)%sigma = max(abs(sigmaKm) * 1d5,1d0)
       events(i)%so2ColumnFlux = max(so2ColumnFlux,0d0)
       events(i)%ashTau550 = max(ashTau550,0d0)
       events(i)%ashLifetime = ashLifetimeDay * secondsPerDay
       events(i)%ashSettling = max(ashSettling,0d0)
       events(i)%ashWavelengthExp = ashWavelengthExp
    end do
    close(unitEvent)

    volcanoEnabled = .true.
    needsSO2 = .false.
    do i=1,volcanoEventsNumber
       if(events(i)%so2ColumnFlux>0d0) needsSO2 = .true.
    end do
    if(needsSO2) i = patmo_volc_getSO2Index(required=.true.)

    print *,"Loaded volcano events from ",trim(fname),": ",volcanoEventsNumber

  end subroutine patmo_volc_loadEvents

  !***************
  subroutine patmo_volc_injectSO2Column(centerKm,sigmaKm,columnSO2)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::centerKm,sigmaKm,columnSO2
    real*8::w(cellsNumber),norm,center,sigma
    integer::j,idxSO2

    if(columnSO2<=0d0) return

    idxSO2 = patmo_volc_getSO2Index(required=.true.)
    center = centerKm * 1d5
    sigma = max(abs(sigmaKm)*1d5,1d0)
    call patmo_volc_layerWeights(center,sigma,w,norm)
    if(norm<=0d0) return

    do j=1,cellsNumber
       nall(j,idxSO2) = nall(j,idxSO2) + columnSO2 * w(j) / norm
    end do

  end subroutine patmo_volc_injectSO2Column

  !***************
  subroutine patmo_volc_addSources(tlocal,n,dn)
    use patmo_commons
    implicit none
    real*8,intent(in)::tlocal
    real*8,intent(in)::n(cellsNumber,speciesNumber)
    real*8,intent(inout)::dn(cellsNumber,speciesNumber)
    real*8::absoluteTime,w(cellsNumber),norm
    integer::i,j,idxSO2

    if(.not.volcanoEnabled) return
    if(volcanoEventsNumber<=0) return

    idxSO2 = patmo_volc_getSO2Index(required=.false.)
    if(idxSO2<0) then
       if(.not.warnedMissingSO2) then
          print *,"WARNING: volcano SO2 source requested, but SO2 is absent."
          warnedMissingSO2 = .true.
       end if
       return
    end if

    absoluteTime = volcanoElapsedTime + max(tlocal,0d0)
    do i=1,volcanoEventsNumber
       if(.not.patmo_volc_so2Active(i,absoluteTime)) cycle
       if(events(i)%so2ColumnFlux<=0d0) cycle
       call patmo_volc_layerWeights(events(i)%center,events(i)%sigma,w,norm)
       if(norm<=0d0) cycle
       do j=1,cellsNumber
          dn(j,idxSO2) = dn(j,idxSO2) + events(i)%so2ColumnFlux * w(j) / norm
       end do
    end do

  end subroutine patmo_volc_addSources

  !***************
  subroutine patmo_volc_applyAshOpacity(tau)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(inout)::tau(photoBinsNumber,cellsNumber)
    real*8::absoluteTime,w(cellsNumber),norm
    real*8::center,age,ashFactor,columnTau(photoBinsNumber)
    real*8::layerFraction,scale
    integer::i,j,k

    if(.not.volcanoEnabled) return
    if(volcanoEventsNumber<=0) return

    absoluteTime = volcanoElapsedTime
    do k=1,volcanoEventsNumber
       ashFactor = patmo_volc_ashFactor(k,absoluteTime)
       if(ashFactor<=0d0) cycle
       if(events(k)%ashTau550<=0d0) cycle

       age = max(absoluteTime - events(k)%startTime,0d0)
       center = events(k)%center - events(k)%ashSettling * age
       call patmo_volc_layerWeights(center,events(k)%sigma,w,norm)
       if(norm<=0d0) cycle

       columnTau(:) = 0d0
       do j=cellsNumber-1,1,-1
          layerFraction = w(j) * max(gridSpace(j),1d0) / norm
          do i=1,photoBinsNumber
             scale = patmo_volc_ashSpectralScale(i,events(k)%ashWavelengthExp)
             columnTau(i) = columnTau(i) + events(k)%ashTau550 &
                  * ashFactor * layerFraction * scale
             tau(i,j) = tau(i,j) + columnTau(i)
          end do
       end do
    end do

  end subroutine patmo_volc_applyAshOpacity

  !***************
  subroutine patmo_volc_dumpState(fname)
    use patmo_commons
    use patmo_parameters
    use patmo_constants
    implicit none
    character(len=*),intent(in)::fname
    real*8::source(cellsNumber),ashLocal(cellsNumber),ashColumn(cellsNumber)
    real*8::w(cellsNumber),norm,absoluteTime,ashFactor,age,center
    real*8::layerFraction
    integer::i,j,idxSO2,unitOut

    source(:) = 0d0
    ashLocal(:) = 0d0
    ashColumn(:) = 0d0
    absoluteTime = volcanoElapsedTime

    if(volcanoEnabled.and.volcanoEventsNumber>0) then
       idxSO2 = patmo_volc_getSO2Index(required=.false.)
       do i=1,volcanoEventsNumber
          if(idxSO2>0.and.patmo_volc_so2Active(i,absoluteTime)) then
             call patmo_volc_layerWeights(events(i)%center,events(i)%sigma,w,norm)
             if(norm>0d0) then
                do j=1,cellsNumber
                   source(j) = source(j) + events(i)%so2ColumnFlux * w(j) / norm
                end do
             end if
          end if

          ashFactor = patmo_volc_ashFactor(i,absoluteTime)
          if(ashFactor>0d0.and.events(i)%ashTau550>0d0) then
             age = max(absoluteTime - events(i)%startTime,0d0)
             center = events(i)%center - events(i)%ashSettling * age
             call patmo_volc_layerWeights(center,events(i)%sigma,w,norm)
             if(norm>0d0) then
                do j=1,cellsNumber
                   layerFraction = w(j) * max(gridSpace(j),1d0) / norm
                   ashLocal(j) = ashLocal(j) + events(i)%ashTau550 &
                        * ashFactor * layerFraction
                end do
             end if
          end if
       end do
    end if

    do j=cellsNumber-1,1,-1
       ashColumn(j) = ashColumn(j+1) + ashLocal(j)
    end do

    unitOut = 92
    open(unitOut,file=trim(fname),status="replace")
    write(unitOut,'(A)') "# altitude_km so2_source_cm-3_s-1 ash_local_tau550 ash_column_tau550"
    do j=1,cellsNumber
       write(unitOut,'(4E17.8E3)') height(j)/1d5,source(j),ashLocal(j),ashColumn(j)
    end do
    close(unitOut)

  end subroutine patmo_volc_dumpState

  !***************
  function patmo_volc_getSO2Index(required)
    use patmo_utils
    implicit none
    logical,intent(in)::required
    integer::patmo_volc_getSO2Index

    if(so2IndexCache==-2) so2IndexCache = getSpeciesIndex("SO2",error=.false.)
    if(required.and.so2IndexCache<0) then
       print *,"ERROR: volcano module requires SO2 in the reaction network."
       stop
    end if
    patmo_volc_getSO2Index = so2IndexCache

  end function patmo_volc_getSO2Index

  !***************
  function patmo_volc_so2Active(eventIndex,timeSeconds)
    implicit none
    integer,intent(in)::eventIndex
    real*8,intent(in)::timeSeconds
    logical::patmo_volc_so2Active

    patmo_volc_so2Active = .false.
    if(eventIndex<1.or.eventIndex>volcanoEventsNumber) return
    if(events(eventIndex)%duration<=0d0) return
    if(timeSeconds<events(eventIndex)%startTime) return
    if(timeSeconds>events(eventIndex)%startTime+events(eventIndex)%duration) return
    patmo_volc_so2Active = .true.

  end function patmo_volc_so2Active

  !***************
  function patmo_volc_ashFactor(eventIndex,timeSeconds)
    implicit none
    integer,intent(in)::eventIndex
    real*8,intent(in)::timeSeconds
    real*8::patmo_volc_ashFactor,age,decayAge

    patmo_volc_ashFactor = 0d0
    if(eventIndex<1.or.eventIndex>volcanoEventsNumber) return
    if(timeSeconds<events(eventIndex)%startTime) return

    age = timeSeconds - events(eventIndex)%startTime
    if(age<=events(eventIndex)%duration) then
       patmo_volc_ashFactor = 1d0
       return
    end if

    if(events(eventIndex)%ashLifetime<0d0) then
       patmo_volc_ashFactor = 1d0
    elseif(events(eventIndex)%ashLifetime>0d0) then
       decayAge = age - events(eventIndex)%duration
       patmo_volc_ashFactor = exp(-decayAge/events(eventIndex)%ashLifetime)
    else
       patmo_volc_ashFactor = 0d0
    end if

  end function patmo_volc_ashFactor

  !***************
  subroutine patmo_volc_layerWeights(center,sigma,w,norm)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::center,sigma
    real*8,intent(out)::w(cellsNumber),norm
    real*8::arg,dz,bestDistance,distance,sig
    integer::j,bestCell

    sig = max(abs(sigma),1d0)
    norm = 0d0
    do j=1,cellsNumber
       arg = (height(j)-center)/sig
       if(abs(arg)>40d0) then
          w(j) = 0d0
       else
          w(j) = exp(-0.5d0*arg*arg)
       end if
       dz = max(gridSpace(j),1d0)
       norm = norm + w(j) * dz
    end do

    if(norm>0d0) return

    bestCell = 1
    bestDistance = abs(height(1)-center)
    do j=2,cellsNumber
       distance = abs(height(j)-center)
       if(distance<bestDistance) then
          bestDistance = distance
          bestCell = j
       end if
    end do
    w(:) = 0d0
    w(bestCell) = 1d0
    norm = max(gridSpace(bestCell),1d0)

  end subroutine patmo_volc_layerWeights

  !***************
  function patmo_volc_ashSpectralScale(ibin,wavelengthExp)
    use patmo_commons
    use patmo_parameters
    use patmo_constants
    implicit none
    integer,intent(in)::ibin
    real*8,intent(in)::wavelengthExp
    real*8::patmo_volc_ashSpectralScale,lambdaNm

    patmo_volc_ashSpectralScale = 1d0
    if(abs(wavelengthExp)<=1d-99) return
    if(energyMid(ibin)<=0d0) return

    lambdaNm = 1d7 * planck_eV * clight / energyMid(ibin)
    if(lambdaNm<=0d0) return

    patmo_volc_ashSpectralScale = (lambdaNm/550d0)**(-wavelengthExp)
    patmo_volc_ashSpectralScale = min(max(patmo_volc_ashSpectralScale,1d-6),1d6)

  end function patmo_volc_ashSpectralScale

  !***************
  function patmo_volc_commentStart(line)
    implicit none
    character(len=*),intent(in)::line
    integer::patmo_volc_commentStart,hashPos,bangPos

    hashPos = index(line,"#")
    bangPos = index(line,"!")
    if(hashPos>0.and.bangPos>0) then
       patmo_volc_commentStart = min(hashPos,bangPos)
    elseif(hashPos>0) then
       patmo_volc_commentStart = hashPos
    elseif(bangPos>0) then
       patmo_volc_commentStart = bangPos
    else
       patmo_volc_commentStart = 0
    end if

  end function patmo_volc_commentStart

  !***************
  subroutine patmo_volc_parseKeywordEvent(line,startDay,durationDay, &
       centerKm,sigmaKm,so2ColumnFlux,ashTau550,ashLifetimeDay, &
       ashSettling,ashWavelengthExp,ios)
    implicit none
    character(len=*),intent(in)::line
    real*8,intent(inout)::startDay,durationDay,centerKm,sigmaKm
    real*8,intent(inout)::so2ColumnFlux,ashTau550,ashLifetimeDay
    real*8,intent(inout)::ashSettling,ashWavelengthExp
    integer,intent(out)::ios
    character(len=512)::work
    character(len=80)::key,value
    integer::eqPos,keyStart,keyEnd,valueStart,valueEnd,nline,readIos
    logical::hasStart,hasDuration,hasCenter,hasSigma,hasSource
    real*8::tmp

    ios = 0
    hasStart = .false.
    hasDuration = .false.
    hasCenter = .false.
    hasSigma = .false.
    hasSource = .false.

    work = adjustl(line)
    do eqPos=1,len(work)
       if(work(eqPos:eqPos)==",".or.work(eqPos:eqPos)==char(9)) work(eqPos:eqPos) = " "
    end do
    nline = len_trim(work)

    do
       eqPos = index(work,"=")
       if(eqPos<=0) exit

       keyEnd = eqPos - 1
       do while(keyEnd>=1.and.work(keyEnd:keyEnd)==" ")
          keyEnd = keyEnd - 1
       end do

       keyStart = keyEnd
       do while(keyStart>1.and.work(keyStart-1:keyStart-1)/=" ")
          keyStart = keyStart - 1
       end do

       valueStart = eqPos + 1
       do while(valueStart<=nline.and.work(valueStart:valueStart)==" ")
          valueStart = valueStart + 1
       end do

       valueEnd = valueStart
       do while(valueEnd<=nline)
          if(work(valueEnd:valueEnd)==" ".or.work(valueEnd:valueEnd)==",") exit
          valueEnd = valueEnd + 1
       end do
       valueEnd = valueEnd - 1

       key = " "
       value = " "
       if(keyStart<=keyEnd) key = patmo_volc_lower(work(keyStart:keyEnd))
       if(valueStart<=valueEnd) value = adjustl(work(valueStart:valueEnd))

       select case(trim(key))
       case("event_id","event","name","scenario")
          continue
       case("start_day","start")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          startDay = tmp
          hasStart = .true.
       case("duration_day","duration")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          durationDay = tmp
          hasDuration = .true.
       case("plume_center_km","center_km","injection_center_km")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          centerKm = tmp
          hasCenter = .true.
       case("plume_sigma_km","sigma_km","injection_sigma_km")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          sigmaKm = tmp
          hasSigma = .true.
       case("so2_flux_cm2_s","so2_column_flux_cm2_s","so2_column_flux")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          so2ColumnFlux = tmp
          hasSource = .true.
       case("ash_tau_550","ash_optical_depth_550","ash_aod_550")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          ashTau550 = tmp
          hasSource = .true.
       case("ash_lifetime_day","ash_decay_day")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          ashLifetimeDay = tmp
       case("ash_settling_cm_s","ash_fall_speed_cm_s")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          ashSettling = tmp
       case("ash_lambda_exponent","ash_wavelength_exp","ash_angstrom_exp")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          ashWavelengthExp = tmp
       case default
          print *,"WARNING: unknown volcano event key ignored: ",trim(key)
       end select

       if(valueEnd>=nline) exit
       work = adjustl(work(valueEnd+1:))
       nline = len_trim(work)
    end do

    if(.not.(hasStart.and.hasDuration.and.hasCenter.and.hasSigma.and.hasSource)) then
       ios = 1
    end if

  end subroutine patmo_volc_parseKeywordEvent

  !***************
  function patmo_volc_lower(text)
    implicit none
    character(len=*),intent(in)::text
    character(len=len(text))::patmo_volc_lower
    integer::i,ich

    do i=1,len(text)
       ich = iachar(text(i:i))
       if(ich>=iachar("A").and.ich<=iachar("Z")) then
          patmo_volc_lower(i:i) = achar(ich + iachar("a") - iachar("A"))
       else
          patmo_volc_lower(i:i) = text(i:i)
       end if
    end do

  end function patmo_volc_lower

end module patmo_volc
