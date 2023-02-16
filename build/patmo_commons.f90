module patmo_commons
  implicit none

  integer,parameter::reactionsNumber = 61
  integer,parameter::chemReactionsNumber = 27
  integer,parameter::photoReactionsNumber = 7
  integer,parameter::reverseReactionsNumber = 27
  integer,parameter::chemSpeciesNumber = 20
  integer,parameter::speciesNumber = 22
  integer,parameter::positionTgas = 21
  integer,parameter::positionDummy = 22
  integer,parameter::cellsNumber = 60
  integer,parameter::photoBinsNumber = 4440
  integer,parameter::patmo_idx_COS = 1
  integer,parameter::patmo_idx_CO = 2
  integer,parameter::patmo_idx_SCSOH = 3
  integer,parameter::patmo_idx_OH = 4
  integer,parameter::patmo_idx_S2 = 5
  integer,parameter::patmo_idx_SH = 6
  integer,parameter::patmo_idx_HSO = 7
  integer,parameter::patmo_idx_O = 8
  integer,parameter::patmo_idx_S = 9
  integer,parameter::patmo_idx_SO2 = 10
  integer,parameter::patmo_idx_SO = 11
  integer,parameter::patmo_idx_CS2E = 12
  integer,parameter::patmo_idx_H = 13
  integer,parameter::patmo_idx_CS = 14
  integer,parameter::patmo_idx_HO2 = 15
  integer,parameter::patmo_idx_N2 = 16
  integer,parameter::patmo_idx_HSO2 = 17
  integer,parameter::patmo_idx_O2 = 18
  integer,parameter::patmo_idx_O3 = 19
  integer,parameter::patmo_idx_CS2 = 20

  integer,parameter::chemReactionsOffset = 0
  integer,parameter::photoReactionsOffset = chemReactionsNumber
  integer,parameter::reverseReactionsOffset = &
      photoReactionsOffset + photoReactionsNumber

  integer,parameter::neqAll = speciesNumber*cellsNumber
  integer,parameter::maxNameLength = 50

  integer,dimension(photoReactionsNumber)::photoPartnerIndex = (/patmo_idx_COS,patmo_idx_SO2,patmo_idx_SO,patmo_idx_O2,patmo_idx_O3,patmo_idx_CS2,patmo_idx_CS2/)

  integer,parameter,dimension(reactionsNumber)::indexReactants2 = (/patmo_idx_O2,&
      patmo_idx_N2,&
      patmo_idx_O2,&
      patmo_idx_OH,&
      patmo_idx_OH,&
      positionDummy,&
      patmo_idx_O2,&
      patmo_idx_O,&
      patmo_idx_O,&
      patmo_idx_O,&
      patmo_idx_O,&
      patmo_idx_O2,&
      patmo_idx_O2,&
      patmo_idx_O3,&
      patmo_idx_O2,&
      patmo_idx_O3,&
      patmo_idx_OH,&
      patmo_idx_O,&
      patmo_idx_O2,&
      patmo_idx_O3,&
      patmo_idx_OH,&
      patmo_idx_O,&
      patmo_idx_O2,&
      patmo_idx_O3,&
      patmo_idx_O2,&
      patmo_idx_O3,&
      patmo_idx_O2,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      patmo_idx_SO2,&
      patmo_idx_SH,&
      positionDummy,&
      patmo_idx_OH,&
      patmo_idx_HSO2,&
      patmo_idx_SO,&
      patmo_idx_S,&
      patmo_idx_CO,&
      patmo_idx_CO,&
      patmo_idx_O,&
      patmo_idx_CO,&
      patmo_idx_O2,&
      patmo_idx_O,&
      patmo_idx_O2,&
      patmo_idx_H,&
      patmo_idx_SO,&
      patmo_idx_O,&
      patmo_idx_O2,&
      patmo_idx_H,&
      patmo_idx_H,&
      patmo_idx_OH,&
      patmo_idx_O2,&
      patmo_idx_OH,&
      patmo_idx_SH,&
      patmo_idx_HO2/)
  integer,parameter,dimension(reactionsNumber)::indexReactants1 = (/patmo_idx_CS2E,&
      patmo_idx_CS2E,&
      patmo_idx_CS2E,&
      patmo_idx_CS2,&
      patmo_idx_CS2,&
      patmo_idx_SCSOH,&
      patmo_idx_SCSOH,&
      patmo_idx_CS2,&
      patmo_idx_CS2,&
      patmo_idx_CS2,&
      patmo_idx_CS,&
      patmo_idx_CS,&
      patmo_idx_CS,&
      patmo_idx_CS,&
      patmo_idx_S,&
      patmo_idx_S,&
      patmo_idx_S,&
      patmo_idx_S2,&
      patmo_idx_SO,&
      patmo_idx_SO,&
      patmo_idx_SO,&
      patmo_idx_SH,&
      patmo_idx_SH,&
      patmo_idx_SH,&
      patmo_idx_HSO,&
      patmo_idx_HSO,&
      patmo_idx_HSO2,&
      patmo_idx_COS,&
      patmo_idx_SO2,&
      patmo_idx_SO,&
      patmo_idx_O2,&
      patmo_idx_O3,&
      patmo_idx_CS2,&
      patmo_idx_CS2,&
      patmo_idx_CS2,&
      patmo_idx_CS2,&
      patmo_idx_CS,&
      patmo_idx_COS,&
      patmo_idx_SCSOH,&
      patmo_idx_CS2,&
      patmo_idx_COS,&
      patmo_idx_CS,&
      patmo_idx_COS,&
      patmo_idx_S2,&
      patmo_idx_S,&
      patmo_idx_COS,&
      patmo_idx_SO,&
      patmo_idx_COS,&
      patmo_idx_SO,&
      patmo_idx_SO,&
      patmo_idx_SO2,&
      patmo_idx_S,&
      patmo_idx_SO2,&
      patmo_idx_SO2,&
      patmo_idx_SO2,&
      patmo_idx_SO,&
      patmo_idx_SO,&
      patmo_idx_HSO,&
      patmo_idx_SO2,&
      patmo_idx_SO2,&
      patmo_idx_SO2/)
  integer,parameter,dimension(reactionsNumber)::indexProducts2 = (/positionDummy,&
      positionDummy,&
      patmo_idx_SO2,&
      patmo_idx_SH,&
      positionDummy,&
      patmo_idx_OH,&
      patmo_idx_HSO2,&
      patmo_idx_SO,&
      patmo_idx_S,&
      patmo_idx_CO,&
      patmo_idx_CO,&
      patmo_idx_O,&
      patmo_idx_CO,&
      patmo_idx_O2,&
      patmo_idx_O,&
      patmo_idx_O2,&
      patmo_idx_H,&
      patmo_idx_SO,&
      patmo_idx_O,&
      patmo_idx_O2,&
      patmo_idx_H,&
      patmo_idx_H,&
      patmo_idx_OH,&
      patmo_idx_O2,&
      patmo_idx_OH,&
      patmo_idx_SH,&
      patmo_idx_HO2,&
      patmo_idx_S,&
      patmo_idx_O,&
      patmo_idx_O,&
      patmo_idx_O,&
      patmo_idx_O,&
      patmo_idx_S,&
      positionDummy,&
      patmo_idx_O2,&
      patmo_idx_N2,&
      patmo_idx_O2,&
      patmo_idx_OH,&
      patmo_idx_OH,&
      positionDummy,&
      patmo_idx_O2,&
      patmo_idx_O,&
      patmo_idx_O,&
      patmo_idx_O,&
      patmo_idx_O,&
      patmo_idx_O2,&
      patmo_idx_O2,&
      patmo_idx_O3,&
      patmo_idx_O2,&
      patmo_idx_O3,&
      patmo_idx_OH,&
      patmo_idx_O,&
      patmo_idx_O2,&
      patmo_idx_O3,&
      patmo_idx_OH,&
      patmo_idx_O,&
      patmo_idx_O2,&
      patmo_idx_O3,&
      patmo_idx_O2,&
      patmo_idx_O3,&
      patmo_idx_O2/)
  integer,parameter,dimension(reactionsNumber)::indexProducts1 = (/patmo_idx_CS2,&
      patmo_idx_CS2,&
      patmo_idx_CS,&
      patmo_idx_COS,&
      patmo_idx_SCSOH,&
      patmo_idx_CS2,&
      patmo_idx_COS,&
      patmo_idx_CS,&
      patmo_idx_COS,&
      patmo_idx_S2,&
      patmo_idx_S,&
      patmo_idx_COS,&
      patmo_idx_SO,&
      patmo_idx_COS,&
      patmo_idx_SO,&
      patmo_idx_SO,&
      patmo_idx_SO2,&
      patmo_idx_S,&
      patmo_idx_SO2,&
      patmo_idx_SO2,&
      patmo_idx_SO2,&
      patmo_idx_SO,&
      patmo_idx_SO,&
      patmo_idx_HSO,&
      patmo_idx_SO2,&
      patmo_idx_SO2,&
      patmo_idx_SO2,&
      patmo_idx_CO,&
      patmo_idx_SO,&
      patmo_idx_S,&
      patmo_idx_O,&
      patmo_idx_O2,&
      patmo_idx_CS,&
      patmo_idx_CS2E,&
      patmo_idx_CS2E,&
      patmo_idx_CS2E,&
      patmo_idx_CS2E,&
      patmo_idx_CS2,&
      patmo_idx_CS2,&
      patmo_idx_SCSOH,&
      patmo_idx_SCSOH,&
      patmo_idx_CS2,&
      patmo_idx_CS2,&
      patmo_idx_CS2,&
      patmo_idx_CS,&
      patmo_idx_CS,&
      patmo_idx_CS,&
      patmo_idx_CS,&
      patmo_idx_S,&
      patmo_idx_S,&
      patmo_idx_S,&
      patmo_idx_S2,&
      patmo_idx_SO,&
      patmo_idx_SO,&
      patmo_idx_SO,&
      patmo_idx_SH,&
      patmo_idx_SH,&
      patmo_idx_SH,&
      patmo_idx_HSO,&
      patmo_idx_HSO,&
      patmo_idx_HSO2/)

end module patmo_commons
