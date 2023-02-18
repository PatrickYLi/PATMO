module patmo_budget
        use patmo_commons
         implicit none
         integer,parameter::nprocess = 56
         integer,parameter::patmo_npr_chemprod = 1  
         integer,parameter::patmo_npr_chemloss = 2 
         integer,parameter::patmo_npr_diffup = 3
         integer,parameter::patmo_npr_diffdown = 4
         integer,parameter::patmo_npr_emission = 5
         integer,parameter::patmo_npr_drydep = 6
         integer,parameter::patmo_npr_wetdep = 7
         integer,parameter::patmo_npr_SO2_OH = 8
         integer,parameter::patmo_npr_SO3_H2O = 9 
         integer,parameter::patmo_npr_photo = 10
         integer,parameter::patmo_npr_CO2_SH = 11
         integer,parameter::patmo_npr_CO_SO = 12
         integer,parameter::patmo_npr_COS_OH = 13
         integer,parameter::patmo_npr_COS_O = 14
         integer,parameter::patmo_npr_COS_hv = 15
        
         integer,parameter::patmo_npr_16 = 16
         integer,parameter::patmo_npr_17 = 17  
         integer,parameter::patmo_npr_18 = 19 
         integer,parameter::patmo_npr_20 = 20
         integer,parameter::patmo_npr_21 = 21
         integer,parameter::patmo_npr_22 = 22
         integer,parameter::patmo_npr_23 = 23
         integer,parameter::patmo_npr_24 = 24
         integer,parameter::patmo_npr_25 = 26
         integer,parameter::patmo_npr_27 = 28 
         integer,parameter::patmo_npr_29 = 29
         integer,parameter::patmo_npr_30 = 30
         integer,parameter::patmo_npr_31 = 31
         integer,parameter::patmo_npr_32 = 32
         integer,parameter::patmo_npr_33 = 33
         integer,parameter::patmo_npr_34 = 34


         integer,parameter::patmo_npr_35 = 35
         integer,parameter::patmo_npr_36 = 36  
         integer,parameter::patmo_npr_37 = 37 
         integer,parameter::patmo_npr_38 = 38
         integer,parameter::patmo_npr_39 = 39
         integer,parameter::patmo_npr_40 = 40
         integer,parameter::patmo_npr_41 = 41
         integer,parameter::patmo_npr_42 = 42
         integer,parameter::patmo_npr_43 = 43
         integer,parameter::patmo_npr_44 = 44 
         integer,parameter::patmo_npr_45 = 45
         integer,parameter::patmo_npr_46 = 46
         integer,parameter::patmo_npr_47 = 47
         integer,parameter::patmo_npr_48 = 48
         integer,parameter::patmo_npr_49 = 49
         integer,parameter::patmo_npr_50 = 50
         integer,parameter::patmo_npr_51 = 51 
         integer,parameter::patmo_npr_52 = 52
         integer,parameter::patmo_npr_53 = 53
         integer,parameter::patmo_npr_54 = 54
         integer,parameter::patmo_npr_55 = 55
         integer,parameter::patmo_npr_56 = 56
         integer,parameter::patmo_npr_57 = 57
         
         integer,parameter::patmo_npr_58 = 58
         integer,parameter::patmo_npr_59 = 59
         integer,parameter,dimension(nprocess)::indexProcess3 = (/&
         patmo_npr_chemprod,&
         patmo_npr_chemloss,&
         patmo_npr_diffup,&
         patmo_npr_diffdown,&
         patmo_npr_emission,&
         patmo_npr_drydep,&
         patmo_npr_wetdep,&
         patmo_npr_SO2_OH,&
         patmo_npr_SO3_H2O,&
         patmo_npr_photo,&
         patmo_npr_CO2_SH,&
         patmo_npr_CO_SO,&
         patmo_npr_COS_OH,&
         patmo_npr_COS_O,&
         patmo_npr_COS_hv,&
        patmo_npr_16,& 
        patmo_npr_17,& 
        patmo_npr_18,& 
        patmo_npr_20,& 
        patmo_npr_21,& 
        patmo_npr_22,& 
        patmo_npr_23,& 
        patmo_npr_24,& 
        patmo_npr_25,& 
        patmo_npr_27,& 
        patmo_npr_29,& 
        patmo_npr_30,& 
        patmo_npr_31,& 
        patmo_npr_32,& 
        patmo_npr_33,& 
        patmo_npr_34,&             
        patmo_npr_35,& 
        patmo_npr_36,& 
        patmo_npr_37,& 
        patmo_npr_38,& 
        patmo_npr_39,& 
        patmo_npr_40,& 
        patmo_npr_41,& 
        patmo_npr_42,& 
        patmo_npr_43,& 
        patmo_npr_44,& 
        patmo_npr_45,& 
        patmo_npr_46,& 
        patmo_npr_47,& 
        patmo_npr_48,& 
        patmo_npr_49,& 
        patmo_npr_50,& 
        patmo_npr_51,& 
        patmo_npr_52,& 
        patmo_npr_53,& 
        patmo_npr_54,& 
        patmo_npr_55,& 
        patmo_npr_56,& 
        patmo_npr_57,&
        patmo_npr_58,&
        patmo_npr_59/)
          real*8::budget(cellsNumber,speciesNumber,nprocess)
 end module patmo_budget
