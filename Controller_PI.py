import matplotlib.pyplot as plt
import numpy as np
import csv
import math
import time
import pickle
# from scipy.optimize import least_squares
# from scipy.signal import medfilt
from matplotlib import rcParams
# from scipy import signal

def control_loop(PI_Gains,ESC_Parameters,Baseline_Inj_Parameters,Control_Config,Engine_Torque,Time,Controller_Data):

    ## Injection Pressure
    KP_PRS=PI_Gains[0][1]
    KI_PRS=PI_Gains[1][1]
    
    ## Injection Duration
    KP_DUR=PI_Gains[2][1]
    KI_DUR=PI_Gains[3][1]

    ## Start of Injection
    KP_SOI=PI_Gains[4][1]
    KI_SOI=PI_Gains[5][1]

    ## Extremum Seeking Controller Parameters Assignment
    K_ESC    = ESC_Parameters[0][1]

    ## Perturbation Frequency
    FP       = ESC_Parameters[1][1]
    OMEGA    = ESC_Parameters[2][1]

    ##Filter Frequency
    OMEGA_HP = ESC_Parameters[3][1]
    OMEGA_LP = ESC_Parameters[4][1]

    ##Perturbation Amplitude
    AM       = ESC_Parameters[5][1]
    BM       = ESC_Parameters[6][1]
    PHASE    = ESC_Parameters[7][1]

    ## Reference Setpoints
    M_e_setp=Control_Config[0][1]

    ## Baseline_Inj_Parameters
    inj_prs_init=np.array(Baseline_Inj_Parameters[0][1])
    inj_dur_init=np.array(Baseline_Inj_Parameters[1][1])
    soi_init=np.array(Baseline_Inj_Parameters[2][1])
    
    ##Control_Config
    pi_esc = Control_Config[2][1]
    pi_c   = Control_Config[3][1]
    
    ## esc_type - 0 (normal), 1 (butterworth)
    esc_type = int(Control_Config[4][1])
    
    ## Time
    time=Time[0]
    t_offset=Time[1]
    step_size=Time[2]
    iteration=Time[3]
    freq=1/step_size

    ## High pass filter (Butterworth filter)
    butterorder = 1
    ys = np.zeros(butterorder+1)
    HPF = np.zeros(butterorder+1)

    ## Engine Torque (Engine_Torque)
    M_e_all=Engine_Torque

    Controller_Data_prev=Controller_Data
    Controller_Data_next=np.zeros(10)

    uapp_prs=0
    uapp_dur=0
    uapp_soi=0

    M_e=Engine_Torque
    
    ## Torque Error
    M_err=(M_e_setp-M_e)/6
    
    ##Integral Control (Memory Assignment)
    prs_int = Controller_Data_prev[4]
    dur_int = Controller_Data_prev[5]
    
    
    prs_int=prs_int+(step_size*M_err)
    dur_int=prs_int+(step_size*M_err)

    prs_pro = KP_PRS*M_err + KI_PRS*prs_int
    dur_pro = KP_DUR*M_err + KI_DUR*dur_int
    
    
    Controller_Data_next[4]=prs_int
    Controller_Data_next[5]=dur_int

    Controller_Data_next[7]=prs_pro
    Controller_Data_next[8]=dur_pro
    

    if pi_c==1:
        soi_int=Controller_Data_prev[6]
        soi_int=soi_int+(step_size*M_err)
        soi_pro=KP_SOI*M_err + KI_SOI*soi_int
        Controller_Data_next[6]=soi_int
        Controller_Data_next[9]=soi_pro

    elif pi_esc==1:

        ## Extremum Seeking Control
        osc_am=AM*math.sin(OMEGA*(time-t_offset+PHASE))
        osc_bm=BM*math.sin(OMEGA*(time-t_offset+PHASE))

        ## Prev SOI Value
        Qopt = Controller_Data_prev[3]

        if esc_type ==1:

            for j in range(len(ys)):
                ys[j] = M_e/6
                HPF[j]= Controller_Data_prev[0]
            pass              

             ## High pass filter (Butterworth filter)
            butterfreq =OMEGA/(2*np.pi) # in Hz
            butterfreq = butterfreq/(freq/2) # normalize to Nyquist frequency
            b,a = signal.butter(butterorder,butterfreq,'highpass')

            for k in range(butterorder):
                ys[k] = ys[k+1]
                HPF[k] = HPF[k+1]
                
            ys[butterorder] = M_e/6

            HPFnew = 0
            
            for k in range(butterorder+1):
                HPFnew = HPFnew + b[k]*ys[butterorder-k]
            for k in range(1,butterorder+1):
                HPFnew = HPFnew - a[k]*HPF[butterorder-k]
                
            HPF[butterorder] = HPFnew
            
            xi = AM*HPFnew*math.sin(OMEGA*(time-t_offset+PHASE))
            
            Qopt = Qopt + xi*K_ESC*step_size
            
            soi_sys= Qopt + BM*math.sin(OMEGA*(time-t_offset+PHASE))
            
            J_hp=HPFnew
            
            J_lp=0
            
            soi_pro=soi_sys
            
            pass

        elif esc_type==0:
            
            J_hp_old=Controller_Data_prev[0]
            J_lp_old=Controller_Data_prev[1]
            J_hp=(1-(step_size*OMEGA_HP))*J_hp_old+(step_size*OMEGA_HP*(M_e/6))
            J_lp=(1-(step_size*OMEGA_LP))*J_lp_old+(((M_e/6)-J_hp_old)*(step_size*OMEGA_LP*osc_am))
            Qopt=Qopt+(step_size*K_ESC*J_lp_old) 
            soi_sys=Qopt+osc_bm
            soi_pro=soi_sys
            pass

        Controller_Data_next[0]=J_hp
        Controller_Data_next[1]=J_lp
        Controller_Data_next[2]=soi_sys
        Controller_Data_next[3]=Qopt
        pass

    uapp_prs=prs_pro+inj_prs_init
    uapp_dur=dur_pro+inj_dur_init
    uapp_soi=soi_pro+soi_init

    return uapp_prs, uapp_dur, uapp_soi, Controller_Data_next
