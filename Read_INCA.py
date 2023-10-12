import numpy as np

MeasureName = {"prs_tur_down":u"prs_tur_down", 'ang_bpa[0]':u"ang_bpa[0]", "ti_main1_tu_tmp":u"ti_main1_tu_tmp", 'prs_bpa':u"prs_bpa", 'fup_mes':u"fup_mes", "prs_thr_up":u"prs_thr_up",  "prs_thr_up_sp":u"prs_thr_up_sp", 'prs_sp_bpa':u"prs_sp_bpa", "poil":u"poil", 'soi_main1_chk_tu[2]':u"soi_main1_chk_tu[2]", 
               "soi_main1_chk_tu[4]":u"soi_main1_chk_tu[4]", "soi_main1_chk_tu[3]":u"soi_main1_chk_tu[3]", 'prs_eg[0]':u"prs_eg[0]", 'soi_main1_chk_tu[5]':u"soi_main1_chk_tu[5]", "prs_eg_mes[0]":u"prs_eg_mes[0]", "prs_cmpr_up":u"prs_cmpr_up", 'prs_ex':u"prs_ex",  "soi_main1_chk_tu[1]":u"soi_main1_chk_tu[1]", "ti_main1_bas":u"ti_main1_bas", "soi_main1_st":u"soi_main1_st",
               'map':u"map", 'prs_cmpr_down':u"prs_cmpr_down", "prs_egrv_up_mdl":u"prs_egrv_up_mdl", "prs_up_cmp_prot[0]":u"prs_up_cmp_prot[0]",  "amp_mes":u"amp_mes", 'prs_bpa_dif':u"prs_bpa_dif", "egbp_mes[0]":u"egbp_mes[0]", "pfu_l_mes":u"pfu_l_mes", "map_mes":u"map_mes", 'pfu_mes':u"pfu_mes",
               'ti_main1_tu':u"ti_main1_tu", "ti_main1_tu_[0]":u"ti_main1_tu_[0]", 'prs_thr_down':u"prs_thr_down", "soi_main1":u"soi_main1", "ti_main1_tu_[3]":u"ti_main1_tu_[3]", 'n':u"n", 'ti_main1_tu_[1]':u"ti_main1_tu_[1]", "ti_main1_tu_[2]":u"ti_main1_tu_[2]", 'soi_main1_chk_tu[0]':u"soi_main1_chk_tu[0]", "ti_main1_tu_[4]":u"ti_main1_tu_[4]",
               "fup":u"fup",  "ti_main1_tu_[5]":u"ti_main1_tu_[5]", 'soi_main1_tu':u"soi_main1_tu", 'ang_thr[0]':u"ang_thr[0]", "soi_main1_bas":u"soi_main1_bas", 'fup_sp':u"fup_sp" }

def Read_INCA(ExperimentObj, MeasureName): 

    INCA_Data = np.zeros((1, len(MeasureName)))

    INCA_Data[0,0] = (MeasureValues(ExperimentObj, list(MeasureName.keys())[0]))/12.0601766654398 # prs_tur_down
    INCA_Data[0,1] = (MeasureValues(ExperimentObj, list(MeasureName.keys())[1]))/327.68 # ang_bpa[0]
    INCA_Data[0,2] = (MeasureValues(ExperimentObj, list(MeasureName.keys())[2]))/1250.0 # ti_main1_tu_tmp
    INCA_Data[0,3] = (MeasureValues(ExperimentObj, list(MeasureName.keys())[3]))/12.0601766654398 # prs_bpa
    INCA_Data[0,4] = (MeasureValues(ExperimentObj, list(MeasureName.keys())[4]))/327.680000069121 # fup_mes
    INCA_Data[0,5] = (MeasureValues(ExperimentObj, list(MeasureName.keys())[5]))/12.0601766654398 # prs_thr_up
    INCA_Data[0,6] = (MeasureValues(ExperimentObj, list(MeasureName.keys())[6]))/12.0601766654398 # prs_thr_up_sp
    INCA_Data[0,7] = (MeasureValues(ExperimentObj, list(MeasureName.keys())[7]))/12.0601766654398 # prs_sp_bpa
    INCA_Data[0,8] = (MeasureValues(ExperimentObj, list(MeasureName.keys())[8]))/6.03008833271991 # poil
    INCA_Data[0,9] = (-5120 + MeasureValues(ExperimentObj, list(MeasureName.keys())[9]))/42.6666666666667 # soi_main1_chk_tu[2]
    INCA_Data[0,10] = (-5120 + MeasureValues(ExperimentObj, list(MeasureName.keys())[10]))/42.6666666666667 # soi_main1_chk_tu[4]
    INCA_Data[0,11] = (-5120 + MeasureValues(ExperimentObj, list(MeasureName.keys())[11]))/42.6666666666667 # soi_main1_chk_tu[3]
    INCA_Data[0,12] = (MeasureValues(ExperimentObj, list(MeasureName.keys())[12]))/6.03008833271991 # prs_eg[0]
    INCA_Data[0,13] = (-5120 + MeasureValues(ExperimentObj, list(MeasureName.keys())[13]))/42.6666666666667 # soi_main1_chk_tu[5]
    INCA_Data[0,14] = (MeasureValues(ExperimentObj, list(MeasureName.keys())[14]))/6.03008833271991 # prs_eg_mes[0]
    INCA_Data[0,15] = (MeasureValues(ExperimentObj, list(MeasureName.keys())[15]))/12.0601766654398 # prs_cmpr_up
    INCA_Data[0,16] = (MeasureValues(ExperimentObj, list(MeasureName.keys())[16]))/12.0601766654398 # prs_ex
    INCA_Data[0,17] = (-5120 + MeasureValues(ExperimentObj, list(MeasureName.keys())[17]))/42.6666666666667 # soi_main1_chk_tu[1]  
    INCA_Data[0,18] = (MeasureValues(ExperimentObj, list(MeasureName.keys())[18]))/1250.0 # ti_main1_bas
    INCA_Data[0,19] = (-5120 + MeasureValues(ExperimentObj, list(MeasureName.keys())[19]))/42.6666666666667 # soi_main1_st
    INCA_Data[0,20] = (MeasureValues(ExperimentObj, list(MeasureName.keys())[20]))/12.0601766654398 # map
    INCA_Data[0,21] = (MeasureValues(ExperimentObj, list(MeasureName.keys())[21]))/12.0601766654398 # prs_cmpr_down
    INCA_Data[0,22] = (MeasureValues(ExperimentObj, list(MeasureName.keys())[22]))/12.0601766654398 # prs_egrv_up_mdl
    INCA_Data[0,23] = (MeasureValues(ExperimentObj, list(MeasureName.keys())[23]))/12.0601766654398 # prs_up_cmp_prot[0]
    INCA_Data[0,24] = (MeasureValues(ExperimentObj, list(MeasureName.keys())[24]))/12.0601766654398 # amp_mes
    INCA_Data[0,25] = (4.00177668780088e-11 +  MeasureValues(ExperimentObj, list(MeasureName.keys())[25]))/12.0601766654398 # prs_bpa_dif
    INCA_Data[0,26] = (MeasureValues(ExperimentObj, list(MeasureName.keys())[26]))/6.03008833271991 # egbp_mes[0]
    INCA_Data[0,27] = (MeasureValues(ExperimentObj, list(MeasureName.keys())[27]))/6.03008833271991 # pfu_l_mes
    INCA_Data[0,28] = (MeasureValues(ExperimentObj, list(MeasureName.keys())[28]))/12.0601766654398 # map_mes
    INCA_Data[0,29] = (MeasureValues(ExperimentObj, list(MeasureName.keys())[29]))/327.680000069121 # pfu_mes
    INCA_Data[0,30] = (MeasureValues(ExperimentObj, list(MeasureName.keys())[30]))/1250.0 # ti_main1_tu
    INCA_Data[0,31] = (MeasureValues(ExperimentObj, list(MeasureName.keys())[31]))/1250.0 # ti_main1_tu_[0]
    INCA_Data[0,32] = (MeasureValues(ExperimentObj, list(MeasureName.keys())[32]))/12.0601766654398 # prs_thr_down
    INCA_Data[0,33] = (-5120 + MeasureValues(ExperimentObj, list(MeasureName.keys())[33]))/42.6666666666667 # soi_main1
    INCA_Data[0,34] = (MeasureValues(ExperimentObj, list(MeasureName.keys())[34]))/1250.0 # ti_main1_tu_[3]
    INCA_Data[0,35] = (MeasureValues(ExperimentObj, list(MeasureName.keys())[35])) # n
    INCA_Data[0,36] = (MeasureValues(ExperimentObj, list(MeasureName.keys())[36]))/1250.0 # ti_main1_tu_[1]
    INCA_Data[0,37] = (MeasureValues(ExperimentObj, list(MeasureName.keys())[37]))/1250.0 # ti_main1_tu_[2]
    INCA_Data[0,38] = (-5120 + MeasureValues(ExperimentObj, list(MeasureName.keys())[38]))/42.6666666666667 # soi_main1_chk_tu[0]
    INCA_Data[0,39] = (MeasureValues(ExperimentObj, list(MeasureName.keys())[39]))/1250.0 # ti_main1_tu_[4]
    INCA_Data[0,40] = (MeasureValues(ExperimentObj, list(MeasureName.keys())[40]))/327.680000069121 # fup
    INCA_Data[0,41] = (MeasureValues(ExperimentObj, list(MeasureName.keys())[41]))/1250.0 # ti_main1_tu_[5]
    INCA_Data[0,42] = (-5120 + MeasureValues(ExperimentObj, list(MeasureName.keys())[42]))/42.6666666666667 # soi_main1_tu
    INCA_Data[0,43] = (MeasureValues(ExperimentObj, list(MeasureName.keys())[43]))/327.68 # ang_thr[0]
    INCA_Data[0,44] = (-5120 + MeasureValues(ExperimentObj, list(MeasureName.keys())[44]))/42.6666666666667 # soi_main1_bas
    INCA_Data[0,45] = (MeasureValues(ExperimentObj, list(MeasureName.keys())[45]))/327.680000069121 # fup_sp
    
    return INCA_Data

# Read the variable (MEASDOUBLENAME1)
def MeasureValues(ExperimentObj, MEASDOUBLENAME1):
    MeasureElementObj_Double_M1  = ExperimentObj.GetMeasureElement(MEASDOUBLENAME1)
    MeasureValueObj_Double_M1  = MeasureElementObj_Double_M1.GetValue() 
    Value_M1 = MeasureValueObj_Double_M1.GetImplValue()
    return Value_M1