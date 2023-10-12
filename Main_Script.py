from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler
from autograd.scipy.special import logsumexp
import sys
import os
from win32com import client
import time as time_measure
import numpy as np
import pandas as pd
from canlib import canlib
import cantools
from pprint import pprint
import can
from Controller_PI import control_loop
from Read_INCA import *
import pickle
import csv
from collections import namedtuple
gaussian = namedtuple('Gaussian', ['mean', 'var'])
gaussian.__repr__ = lambda s: f'𝒩(μ={s[0]:.3f}, 𝜎²={s[1]:.3f})'

#%%
#------------------------------------------------------------------------------
TOOLADD = u"Inca.Inca.7.2"# INCA V7.x COM address
DBPATH = u"D:\ETASData\INCA7.2\Database\ONR"		# path and name of used database
TFNAME  = u"DEFAULT"                 # name of top folder
SFNAME_EXP  = u"Experiment"         # name of sub folder for Experiments
EXPNAME = u"I334 - 9-11-2018"					 # name of experiment element
SFNAME_PRJ = u"Project"             # name of sub folder for ECU-Projects
PRJNAME = u"NS660000"     # name of ECU-Project Element
SFNAME_WRK = u"Workspace"           # name of sub folder for Workspaces
WSNAME  = u"ONR"       # name of workspace element
WORKINGDIRECTORY = u"D:\ETASData\INCA7.2\Measure"
SYS_NAME = u"TS-Testsystem"          # Global System Name
DEV_NAME = u"ETK-Testdevice"         # Global Device Name
DatasetStr = u"NS660000"
RASTER = u"10ms"
Calibration_Map = "Calibration_Map.csv"

# -----------------------------------------------------------------------------
# Data_Raster = 0.036 # Seconds
Data_Raster = 0.5 # Seconds --> This decides at what interval we are collecting data and changing ECU inputs.
TIME = 60*13# sec
Data_Pts = int(TIME/Data_Raster)

# These are the hard constraints we applied at the end while writing inputs to ECU. 
# This is to make sure that the Python API shold not give any infesible inputs to engine.
InpConstr = {'maxDur':1.0,'minDur':-0.7,'maxPres':12,'minPres':3, 'maxSOI':7,'minSOI':-10}

# -----------------------------------------------------------------------------
# Creating the dict for the name or values we want to set or measure
SetBool = {'lc_fup_sp_man':u"lc_fup_sp_man", 'lc_ti_main1_tu_on':u'lc_ti_main1_tu_on', 'lc_soi_main1_tu_on':u'lc_soi_main1_tu_on'}
SetName = {'c_fup_sp_man':u"c_fup_sp_man", 'c_ti_main1_add_tu[0]':u'c_ti_main1_add_tu[0]', 'c_ti_main1_add_tu[1]':u'c_ti_main1_add_tu[1]', 'c_ti_main1_add_tu[2]':u'c_ti_main1_add_tu[2]', 'c_ti_main1_add_tu[3]':u'c_ti_main1_add_tu[3]', 'c_ti_main1_add_tu[4]':u'c_ti_main1_add_tu[4]', 'c_ti_main1_add_tu[5]':u'c_ti_main1_add_tu[5]', 'c_soi_main1_add_tu':u"c_soi_main1_add_tu"}
SetValue = {'c_ang_sp_man_thr[0]':2}

MeasureName = {"prs_tur_down":u"prs_tur_down", 'ang_bpa[0]':u"ang_bpa[0]", "ti_main1_tu_tmp":u"ti_main1_tu_tmp", 'prs_bpa':u"prs_bpa", 'fup_mes':u"fup_mes", "prs_thr_up":u"prs_thr_up",  "prs_thr_up_sp":u"prs_thr_up_sp", 'prs_sp_bpa':u"prs_sp_bpa", "poil":u"poil", 'soi_main1_chk_tu[2]':u"soi_main1_chk_tu[2]", 
               "soi_main1_chk_tu[4]":u"soi_main1_chk_tu[4]", "soi_main1_chk_tu[3]":u"soi_main1_chk_tu[3]", 'prs_eg[0]':u"prs_eg[0]", 'soi_main1_chk_tu[5]':u"soi_main1_chk_tu[5]", "prs_eg_mes[0]":u"prs_eg_mes[0]", "prs_cmpr_up":u"prs_cmpr_up", 'prs_ex':u"prs_ex",  "soi_main1_chk_tu[1]":u"soi_main1_chk_tu[1]", "ti_main1_bas":u"ti_main1_bas", "soi_main1_st":u"soi_main1_st",
               'map':u"map", 'prs_cmpr_down':u"prs_cmpr_down", "prs_egrv_up_mdl":u"prs_egrv_up_mdl", "prs_up_cmp_prot[0]":u"prs_up_cmp_prot[0]",  "amp_mes":u"amp_mes", 'prs_bpa_dif':u"prs_bpa_dif", "egbp_mes[0]":u"egbp_mes[0]", "pfu_l_mes":u"pfu_l_mes", "map_mes":u"map_mes", 'pfu_mes':u"pfu_mes",
               'ti_main1_tu':u"ti_main1_tu", "ti_main1_tu_[0]":u"ti_main1_tu_[0]", 'prs_thr_down':u"prs_thr_down", "soi_main1":u"soi_main1", "ti_main1_tu_[3]":u"ti_main1_tu_[3]", 'n':u"n", 'ti_main1_tu_[1]':u"ti_main1_tu_[1]", "ti_main1_tu_[2]":u"ti_main1_tu_[2]", 'soi_main1_chk_tu[0]':u"soi_main1_chk_tu[0]", "ti_main1_tu_[4]":u"ti_main1_tu_[4]",
               "fup":u"fup",  "ti_main1_tu_[5]":u"ti_main1_tu_[5]", 'soi_main1_tu':u"soi_main1_tu", 'ang_thr[0]':u"ang_thr[0]", "soi_main1_bas":u"soi_main1_bas", 'fup_sp':u"fup_sp" }

CaliBoolON = 1
CaliBoolOFF = 0
# Create empty list and array
DoubleValue_Bool = []
DoubleValue_S1 = []
INCA_PUMA_Data = np.zeros((int(Data_Pts), 31+len(MeasureName)))
Torque_filter_Array = np.zeros((int(Data_Pts),1))

U_Prs_App_Array = np.zeros((int(Data_Pts), 1))
U_Dur_App_Array = np.zeros((int(Data_Pts), 1))
U_SOI_App_Array = np.zeros((int(Data_Pts), 1))

#------------------------------------------------------------------------------
# PUMA initial things
db = cantools.database.load_file('PUMA.dbc') # A database file: It primarily includes the required conversions from data acquisition system to python. 
# Initialization
channel_number = 0

# Specific CANlib channel number may be specified as first argument
if len(sys.argv) == 2:
    channel_number = int(sys.argv[1])

print("Opening channel %d" % (channel_number))
# Open CAN channel, virtual channels are considered ok to use
ch = canlib.openChannel(channel_number, canlib.canOPEN_ACCEPT_VIRTUAL, canlib.Open.EXCLUSIVE)

print("Setting bitrate to 500 kb/s")
ch.setBusParams(canlib.canBITRATE_500K)
ch.setBusOutputControl(canlib.Driver.SILENT)
ch.busOn()

# Start listening for messages
finished = False
# -----------------------------------------------------------------------------

def main():
    out_model_path = 'All/'
    
    ## Input the plan and create torque setpoints
    plan = pd.read_csv('ToruqeProfile_plan.csv')
    
    change = plan['torque'].ne(plan['torque'].shift())
    plan_change_var = plan.loc[change,['time','torque']]
       
    ## Define speed Setpoint
    Speed_Setpoint = 1400
    #initial torque setpoint
    Torque_Setpoint = plan_change_var['torque'].iloc[0]  
    
    step_size = Data_Raster
    t_offset=int(2/step_size)
    time=0
    
    #Control on/off
    pi_c = 1  ##PI Controller 0-off, 1-on
    pi_esc = 0  #PI+ES Controller 0-off, 1-0n
    esc_type = 0 ## esc_type - 0 (normal), 1 (butterworth)
    
    # Read the CSV file into a dataframe
    df = pd.read_csv(Calibration_Map)
    
    # Get Nominal control input maps 
    nominal_Pic = []
    nominal_InjDur = []
    nominal_SOI = []
    
    
    rows = df[(df["Speed"] == Speed_Setpoint) & (df["Torque"] == Torque_Setpoint)]
    output_vars=["Inj Pressure","Inj Duration","SOI"]
    nominal_control_inputs = rows[output_vars].values.tolist()
    nominal_Pic = np.array(10*1e6)
    nominal_InjDur  = np.array(1)
    nominal_SOI  = np.array(360)
    
   ## PI Controller parameters
    # Injection pressure
    KP_PRS=0*1e1
    KI_PRS=0*1e3
    # Injection Duration
    KP_DUR=0*1e-3
    KI_DUR=0*1e-3
    ## Start of Injection
    KP_SOI=0*1e-3
    KI_SOI=0*1e-2
    
    ## Store as structure
    # PI Gains
    PI_Gains = np.array([('KP_PRS', KP_PRS), ('KI_PRS',  KI_PRS),('KP_DUR',  KP_DUR),('KI_DUR',  KI_DUR),('KP_SOI',  KP_SOI),('KI_SOI',  KI_SOI)],
                 dtype=[('name', 'U10'), ('value',np.float32)])
    
    # Baseline_Injection_Parameters
    Baseline_Inj_Parameters = np.array([('Inj_Prs', nominal_Pic), ('Inj_Dur', nominal_InjDur),('SOI',  nominal_SOI )],
                 dtype=[('name', 'U10'),('value',np.float32)])
    
    # Control_Configuration
    Control_Config = np.array([('Torque_Setpoint', Torque_Setpoint),('Speed_Setpoint', Speed_Setpoint), ('ESC_ON_OFF', pi_esc ),('PI_ON_OFF',  pi_c ), ('ESC_TYPE',esc_type)],
                 dtype=[('name', 'U10'), ('value',np.float32)])
    
    ##Variables to store results
    Torque=[]
    Fuel=[]
    sim_time=[]
    
    SOI_var=[]
    Jobj_HP=[]
    Jobj_LP=[]
    SOI_mean=[]
    Prs_int=[]
    Dur_int=[]
    SOI_int=[]
    Prs_PI=[]
    Dur_PI=[]
    SOI_PI=[]
    
    U_Prs_App=[]
    U_Dur_App=[]
    U_SOI_App=[]
    
    Jobj_HP.append(0)
    Jobj_LP.append(0)
    SOI_var.append(0)
    SOI_mean.append(0)
    Prs_int.append(0)
    Dur_int.append(0)
    SOI_int.append(0)
    Prs_PI.append(0)
    Dur_PI.append(0)
    SOI_PI.append(0)
        
    U_Prs_App.append(Baseline_Inj_Parameters[0][1])
    U_Dur_App.append(Baseline_Inj_Parameters[1][1])
    U_SOI_App.append(Baseline_Inj_Parameters[2][1])
    
    Torque_window=[]
    Torque_filter=[]
    
    inj_prs_init=nominal_Pic
    inj_dur_init=nominal_InjDur
    soi_init=nominal_SOI
    
    iter_loop=1
    
    def predict(pos, movement):
        return gaussian(pos.mean + movement.mean, pos.var + movement.var)

    def gaussian_multiply(g1, g2):
        mean = (g1.var * g2.mean + g2.var * g1.mean) / (g1.var + g2.var)
        variance = (g1.var * g2.var) / (g1.var + g2.var)
        return gaussian(mean, variance)

    def update(prior, likelihood):
        posterior = gaussian_multiply(likelihood, prior)
        return posterior
##--------------------------------------------------------------------------------------------------
    print (u"----- Start Experiment -----")

    # Get connection to INCA and open it, if it was closed, exit script if fails object dispatch
    IncaObj = client.Dispatch(TOOLADD)
    if IncaObj is None:
        print (u"Unable to connect to ", TOOLADD)
        IncaObj.DisconnectFromTool()
        sys.exit(1)
    else:
        print (u"Connected to %s\n" % TOOLADD)

    DBObj = IncaObj.GetCurrentDataBase()
    if DBObj is None:
        print (u"Database not opened!")
        exit()

    ExperimentObj = MyExperiment(IncaObj, DBObj)

    # Is current experiment online?
    if ExperimentObj.IsIncaOnlineExperiment():
        print (u"Experiment %s is active!\n" % EXPNAME)
    else:
        print(u"Experiment %s is not active!\n" % EXPNAME)
        exit()

    # Start Measuring
    if ExperimentObj.StartMeasurement():
        print ("Measurement started!")
        
        PUMA_Data = []
        xloop = 0
        ExperimentObj = IncaObj.GetOpenedExperiment()
        
        while len(PUMA_Data) < Data_Pts:
            Start_Time = time_measure.process_time()
            
            try:
                # Store/get PUMA data
                frame = ch.read(timeout=1)
                CAN_Data1 = db.decode_message(frame.id, frame.data)
                frame = ch.read(timeout=1)
                CAN_Data2 = db.decode_message(frame.id, frame.data)
                frame = ch.read(timeout=1)
                CAN_Data3 = db.decode_message(frame.id, frame.data)
                frame = ch.read(timeout=1)
                CAN_Data4 = db.decode_message(frame.id, frame.data)
                frame = ch.read(timeout=1)
                CAN_Data5 = db.decode_message(frame.id, frame.data)
                frame = ch.read(timeout=1)
                CAN_Data6 = db.decode_message(frame.id, frame.data)
                frame = ch.read(timeout=1)
                CAN_Data7 = db.decode_message(frame.id, frame.data)
                frame = ch.read(timeout=1)
                CAN_Data8 = db.decode_message(frame.id, frame.data)
                Current_Data = list(CAN_Data1.values()) + list(CAN_Data2.values()) +list(CAN_Data3.values()) + list(CAN_Data4.values()) + list(CAN_Data5.values()) + list(CAN_Data6.values()) +list(CAN_Data7.values()) + list(CAN_Data8.values())
                PUMA_Data.append(Current_Data)
                INCA_PUMA_Data[xloop,0:31] = np.array(PUMA_Data[-1])

                # Store/get INCA data
                INCA_PUMA_Data[xloop,31:len(INCA_PUMA_Data[0,:])] = Read_INCA(ExperimentObj, MeasureName)
                
    ##-----------------------------------------------------------------------------------------------------
    
                # COLTROLLER PIECE (Torque Error) if we detect a fault and we need to change the input parameters
                
                #Get the Torque
                Engine_Torque = INCA_PUMA_Data[xloop,30]
                Torque_window.append(Engine_Torque)
                window=int(Speed_Setpoint/30)
                process_var = .001**2
                
                if len(Torque_window)<window:
                    torque_std = np.std(Torque_window)
                    torque_mean= np.mean(Torque_window)
                    # x = gaussian(Torque_Setpoint, 150) # initial state
                    x = gaussian(torque_mean, 90) # initial state
                    process_model = gaussian(0., process_var)
                    ps = []
                    estimates = []
                    zs=Torque_window
    
                    for z in zs:
                        prior = predict(x, process_model)
                        x = update(prior, gaussian(z, torque_std**2))
                        estimates.append(x.mean)
                        ps.append(x.var)
                else:
                    Torque_window1 = Torque_window[len(Torque_window)-window:len(Torque_window)] 
                    torque_mean=np.mean(Torque_window1)
                    torque_std = np.std(Torque_window1)
                    # x = gaussian(Torque_Setpoint, 150) # initial state
                    x = gaussian(torque_mean, 90) # initial state
                    process_model = gaussian(0., process_var)
                    ps = []
                    estimates = []
                    zs=Torque_window1
    
                    for z in zs:
                        prior = predict(x, process_model)
                        x = update(prior, gaussian(z, torque_std**2))
                        estimates.append(x.mean)
                        ps.append(x.var)
                
                Engine_Torque=estimates[-1]
                Torque_filter_Array[xloop,0] = Engine_Torque
                # Control piece if we detect a fault and we need to change the input parameters
                #Time Information
                Time=np.array([time,t_offset,step_size,iter_loop])
    
                # Controller past state 
                Controller_Data=np.zeros(10)
                Controller_Data[0]=Jobj_HP[iter_loop-1]
                Controller_Data[1]=Jobj_LP[iter_loop-1]
                Controller_Data[2]=SOI_var[iter_loop-1]
                Controller_Data[3]=SOI_mean[iter_loop-1]
                Controller_Data[4]=Prs_int[iter_loop-1]
                Controller_Data[5]=Dur_int[iter_loop-1]
                Controller_Data[6]=SOI_int[iter_loop-1]
                Controller_Data[7]=Prs_PI[iter_loop-1]
                Controller_Data[8]=Dur_PI[iter_loop-1]
                Controller_Data[9]=SOI_PI[iter_loop-1]
                
                if np.around(time, decimals=4)>t_offset:
                    u_prs, u_dur, u_soi, Controller_Data_next = control_loop(PI_Gains,ESC_Parameters,Baseline_Inj_Parameters,Control_Config,Engine_Torque,Time,Controller_Data)
                    Jobj_HP.append(Controller_Data_next[0])
                    Jobj_LP.append(Controller_Data_next[1])
                    SOI_var.append(Controller_Data_next[2])
                    SOI_mean.append(Controller_Data_next[3])
                    Prs_int.append(Controller_Data_next[4])
                    Dur_int.append(Controller_Data_next[5])
                    SOI_int.append(Controller_Data_next[6])
                    Prs_PI.append(Controller_Data_next[7])
                    Dur_PI.append(Controller_Data_next[8])
                    SOI_PI.append(Controller_Data_next[9])
    
                    U_Prs_App.append(u_prs)
                    U_Dur_App.append(u_dur)
                    U_SOI_App.append(u_soi)
                else:
                    u_prs=(Baseline_Inj_Parameters[0][1])
                    u_dur=(Baseline_Inj_Parameters[1][1])
                    u_soi=(Baseline_Inj_Parameters[2][1])
    
                    Jobj_HP.append(0)
                    Jobj_LP.append(0)
                    SOI_var.append(0)
                    SOI_mean.append(0)
                    Prs_int.append(0)
                    Dur_int.append(0)
                    SOI_int.append(0)
                    Prs_PI.append(0)
                    Dur_PI.append(0)
                    SOI_PI.append(0)
    
                    U_Prs_App.append(u_prs)
                    U_Dur_App.append(u_dur)
                    U_SOI_App.append(u_soi)
    
                    pass
                U_Prs_App_Array[xloop,0] = u_prs
                U_Dur_App_Array[xloop,0] = u_dur
                U_SOI_App_Array[xloop,0] = u_soi
                
                Torque.append(Engine_Torque)
                sim_time.append(time)
                Fuel.append(INCA_PUMA_Data[xloop,25])
                
                
                
                iter_loop=iter_loop+1
                time=time+step_size
    ##----------------------------------------------------------------------------------------------------------------                

                # Conversions to match with required input format for ECU.
                                        # SOI_main1_st (from INCA data)
                u_soi = (u_soi - 360) - (2.8603) # Need to check with Brian
                u_dur = u_dur - 1.2
                u_prs = u_prs/1e6
                
                # Add Pressure Shift Faults [This is done for experiment purpose to demonstrate fault mitigation].
                if np.around(time, decimals=4) >= 300:
                    u_prs= u_prs + (2)
                if np.around(time, decimals=4) >= 660:
                    u_prs= u_prs + (1)

                # Constraints for control inputs
                while( (u_dur > InpConstr['maxDur']) or (u_dur < InpConstr['minDur'])):
                      if u_dur > InpConstr['maxDur']:
                          u_dur = InpConstr['maxDur']
                      if u_dur < InpConstr['minDur']:
                          u_dur = InpConstr['minDur']
                
                while( (u_prs > InpConstr['maxPres']) or (u_prs < InpConstr['minPres'])):
                      if u_prs > InpConstr['maxPres']:
                          u_prs = InpConstr['maxPres']
                      if u_prs < InpConstr['minPres']:
                          u_prs = InpConstr['minPres']
                
                while( (u_soi > InpConstr['maxSOI']) or (u_soi < InpConstr['minSOI'])):
                      if u_soi > InpConstr['maxSOI']:
                          u_soi = InpConstr['maxSOI']
                      if u_soi < InpConstr['minSOI']:
                          u_soi = InpConstr['minSOI']
                    
                print('Pressure: {} Duration: {} SOI {} Time{}'.format(round(u_prs,2),round(u_dur,2),round(u_soi,2), time))    
                print('Filtered Torque: {} Torque: {}'.format(round(Torque_filter_Array[xloop,0],3),round(INCA_PUMA_Data[xloop,30],2)))
                
                # Setting/writing the control parameter in ECU
                if Current_Data: # If the list is not empty
                    
                    # Injection Duration setting
                    DoubleValue_Bool.append(SetValues(ExperimentObj, SetBool['lc_ti_main1_tu_on'], CaliBoolON))
                    DoubleValue_S1.append(SetValues(ExperimentObj, SetName['c_ti_main1_add_tu[0]'], u_dur))
                    DoubleValue_S1.append(SetValues(ExperimentObj, SetName['c_ti_main1_add_tu[1]'], u_dur))
                    DoubleValue_S1.append(SetValues(ExperimentObj, SetName['c_ti_main1_add_tu[2]'], u_dur))
                    DoubleValue_S1.append(SetValues(ExperimentObj, SetName['c_ti_main1_add_tu[3]'], u_dur))
                    DoubleValue_S1.append(SetValues(ExperimentObj, SetName['c_ti_main1_add_tu[4]'], u_dur))
                    DoubleValue_S1.append(SetValues(ExperimentObj, SetName['c_ti_main1_add_tu[5]'], u_dur))
        
                    # Injection Pressure setting
                    DoubleValue_S1.append(SetValues(ExperimentObj, SetName['c_fup_sp_man'], u_prs))
                    DoubleValue_Bool.append(SetValues(ExperimentObj, SetBool['lc_fup_sp_man'], CaliBoolON))
                
                    # # Start Of Injection setting
                    DoubleValue_S1.append(SetValues(ExperimentObj, SetName['c_soi_main1_add_tu'], u_soi))
                    DoubleValue_Bool.append(SetValues(ExperimentObj, SetBool['lc_soi_main1_tu_on'], CaliBoolON))
                    
                Current_Data = []
                xloop = xloop +1
                ctrl = canlib.IOControl(0)
                ctrl.flush_rx_buffer
                DataPt_Time = time_measure.process_time() - Start_Time
                
            except(canlib.canNoMsg) as ex:
                ctrl = canlib.IOControl(0)
                ctrl.flush_rx_buffer
                None

# -----------------------------------------------------------------------------------------------------------------
    # Stop Measuring
    if ExperimentObj.StopMeasurement():
        print ("Measurement stopped!")

    # Clean up
    if IncaObj.DisconnectFromTool():
        del IncaObj
        print ("\nDisconnected from %s!" % TOOLADD)
        print ("\n-------- End ---------\n")
    else:
        print ("Could not disconnect from %s!!!...\n\n" % TOOLADD)

    # Channel teardown
    ch.busOff()
    ch.close()

#------------------------------------------------------------------------------
def MyExperiment(myIncaObj, myDBObj):
    myExperimentObj = myIncaObj.GetOpenedExperiment()
	# If closed, open the current experiment assigned to the defined workspace of the opened database
    if myExperimentObj is None:
        print ("Experiment is not opened.")
        TFolderObj = myDBObj.GetFolder(TFNAME)
        SFolderExpObj = TFolderObj.GetSubFolder(SFNAME_EXP)
        ExperimentItemObj = SFolderExpObj.GetComponent(EXPNAME)
        SFolderWrkObj = TFolderObj.GetSubFolder(SFNAME_WRK)
        WorkspaceItemObj = SFolderWrkObj.GetComponent(WSNAME)
        ExperimentItemObj.SetHardwareConfiguration(WorkspaceItemObj)
        if ExperimentItemObj.SetHardwareConfiguration(WorkspaceItemObj):
            ExperimentViewObj = ExperimentItemObj.OpenExperiment()
            myExperimentObj = ExperimentViewObj.GetExperiment()
            print ("Experiment opened!")
        else:
            print ("Experiment not connected to workspace!")
            exit()
    return myExperimentObj

# Set the variable or boolean (CALIDOUBLENAME) to value (SetValue)
def SetValues(ExperimentObj, CALIDOUBLENAME, SetValue):
    CalibrationElementObj_Double = ExperimentObj.GetCalibrationElement(CALIDOUBLENAME)
    CalibrationValueObj_Double = CalibrationElementObj_Double.GetValue()
    SetValue_C1 = CalibrationValueObj_Double.SetDoublePhysValue(SetValue)
    DoubleValue_C1 = CalibrationValueObj_Double.GetDoublePhysValue()
    return DoubleValue_C1

# Read the variable (MEASDOUBLENAME1)
def MeasureValues(ExperimentObj, MEASDOUBLENAME1):
    MeasureElementObj_Double_M1  = ExperimentObj.GetMeasureElement(MEASDOUBLENAME1)
    MeasureValueObj_Double_M1  = MeasureElementObj_Double_M1.GetValue() 
    Value_M1 = MeasureValueObj_Double_M1.GetImplValue()
    return Value_M1

main()

#%% Storing the INCA PUMA data in csv files with variable heading 
PUMA_Name_List = ['T_EXH_CYL_1', 'T_EXH_CYL_2', 'T_EXH_CYL_3', 'T_EXH_CYL_4', 'T_EXH_MNF', 'T_TURB_OUT', 'T_EXH_CYL_5', 'T_EXH_CYL_6', 'P_COMP_IN', 'P_COMP_OUT', 'P_EXH_MAN', 'P_TURB_OUT',
                   'T_ENG_COOL_IN', 'T_ENG_COOL_OUT', 'T_ENG_OIL', 'P_ENG_COOL_IN', 'P_ENG_COOL_OUT', 'P_ENG_OIL', 'P_Fuel_Rail', 'T_ENG_AIR_down_CAC', 'T_COMP_IN', 'T_COMP_OUT', 
                   'T_INT_MAN', 'ENG_MF_LFE', 'ENGINE_COOL_VF', 'FB_VAL', 'P_BACKPRESSURE', 'P_BARO', 'SPEED', 'T_CELL', 'TORQUE']

INCA_Name_List = list(MeasureName.keys())
MeasureName_List = PUMA_Name_List + INCA_Name_List

df_INCA_PUMA_Data = pd.DataFrame(INCA_PUMA_Data)
df_INCA_PUMA_Data.columns= MeasureName_List
df_INCA_PUMA_Data.to_csv("ExperimentData_INCA_PUMA.csv")
