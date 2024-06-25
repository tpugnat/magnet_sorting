from datetime import datetime
import os
from os import system
import sys
from subprocess import Popen
sys.path.append("omc3")
from omc3.scripts.fake_measurement_from_model import generate as fake_measurement
from omc3 import global_correction
from omc3.response_creator import create_response_entrypoint
from pairs import CorrectabilityError, Pairs, mask_monitor
from q1_errors import Q1Pairs
from q2_errors import Q2Pairs
from magnet_errors import *
from pathlib import Path
import tfs
import numpy as np
import pandas as pd
from typing import List, Union, Tuple


# Remove FutureWarning
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)



# path to madx executable
#MADX = "/home/awegsche/programs/madx/madx-gnu64"
#MADX = "/home/thpugnat/Documents/CERN/madx"
MADX = "madx"

# path to twiss output
MODEL_NUMBER= "_lxplus1"
MODEL_TWISS = f"model{MODEL_NUMBER}/twiss_err_b1.tfs"

# variable categories for global corrections
# these are the magnet circuits used to correct, check that they match the triplets you want to correct
# e.g. correcting IR1 with IR2 circuits will not yield good results
VAR_CATS = [
        'kqx2a.l1', 'kqx2a.r1',
        'kqx2a.l5', 'kqx2a.r5',
        'kqx2b.l1', 'kqx2b.r1',
        'kqx2b.l5', 'kqx2b.r5',
        'kqx1.l1', 'kqx1.r1',
        'kqx1.l5', 'kqx1.r5',
        #'ktqx1.l1', 'ktqx1.r1',
        #'MQY',
        ]


# TYPE_SCORE can be either "RMS_XY", "RMS_X", "MAXABS_XY", "MAXABS_X"
TYPE_SCORE = 'RMS_XY'; 
# TYPE_ERROR = 'REAL'; 
TYPE_RAND  = 'Gaussian'; 
FLAG_WITH_calibration = True; #False


# turn off certain parts of the sim, leave all at `TRUE` for more than one simulation or you will
# get the same result over and over
DO_MADX = True
DO_FR   = True
DO_CORR = True

# error specification, this has to be kept up to date with estimations from Massimo and Ezio 
# (initially AMP_REAL_ERROR = 50 , AMP_MEAS_ERROR = 2)
AMP_MEAS_ERROR = 10; # 25; # 15;
AMP_CALI_ERROR = 0; #5; #
AMP_PRES_ERROR = 0; #1; #
STAGE = 2;


# maximum number of simulations, `Ctrl-C` stops the program early
MAX_SIMS = 1000 #100 #1 #1000

NUMB_PERMUT_CROSSSCORING = 100

# some more constants
WRITE_SUMMARY_REFRESH_FREQUENCY = 500
#SUMM_PAIRING = "summ_pairing.tfs"
#SUMM_SUM = "summ_sum.tfs"
#SUMM_SUMBETA = "summ_sumbeta.tfs"
#SUMM_SUMBBEATING = "summ_sumbbeating.tfs"
#SUMM_SUMBBEATING_V2 = "summ_sumbbeating_v2.tfs"

if FLAG_WITH_calibration:
    NAME_SAVEFILE_FORMAT = f"{TYPE_SCORE}_dist_m{AMP_MEAS_ERROR}u_{TYPE_RAND}_c{AMP_CALI_ERROR}u_p{AMP_PRES_ERROR}u_nbsimu{MAX_SIMS}_wcor_Phase{STAGE}"
else:
    NAME_SAVEFILE_FORMAT = f"{TYPE_SCORE}_dist_m{AMP_MEAS_ERROR}u_{TYPE_RAND}_c{AMP_CALI_ERROR}u_p{AMP_PRES_ERROR}u_nbsimu{MAX_SIMS}_Phase{STAGE}"
SUMM_SUMBBEATING_V2 = f"summ_score_{NAME_SAVEFILE_FORMAT}.tfs"
#SUMM_DIFF = "summ_diff.tfs"
#SUMM_BBEAT = "summ_bbeat.tfs"

# some flags
FLAG_DEBUG=False





def main():
    """
    """
    if not os.path.exists(f"model/"):
        raise FileNotFoundError("could not find model/")
    
    if not os.path.exists(f"model{MODEL_NUMBER}/"):
        os.makedirs(f"model{MODEL_NUMBER}/")
        cwd = os.getcwd()
        os.symlink(Path(cwd,'model/acc-models-lhc'), f"model{MODEL_NUMBER}/acc-models-lhc")
        os.symlink(Path(cwd,'macros'), f"model{MODEL_NUMBER}/macros")
        system(f"touch model{MODEL_NUMBER}/knobs.madx")
        
    generate_optic()

    summ = Summary()

    try:
        optic_without_error = tfs.read('model1_ft/twiss_elements.dat')
    except:
        print('WARNING: No twiss file found for the sorting of the magnet. The betas will be set to 1.')
        optic_without_error = None

    for i in range(MAX_SIMS):
        summ.final_seed = (MAX_SIMS==i+1)
        do_analysis(summ, optic_without_error=optic_without_error)
        print(summ)

        if FLAG_DEBUG:
            with open(f"sim_summary{datetime.now()}.txt", "w") as summfile:
                summfile.write(str(summ))
                summfile.write("\n")
                
    with open(f"sim_summary_{NAME_SAVEFILE_FORMAT}.txt", "w") as summfile:
                summfile.write(str(summ))
                summfile.write("\n")



#WRITE_SUMMARY_COUNT = 0
# ---- summary (helpfull to check the global state of the simulations) -----------------------------
def write_summary(parameters, filename, summ):
    print("concatenate summary table")
    index = len(summ.WRITE_SUMMARY_TABLE.index)
    summ.WRITE_SUMMARY_TABLE = pd.concat([summ.WRITE_SUMMARY_TABLE, pd.DataFrame(parameters, index=[index])], ignore_index=True)
    
    if summ.final_seed or len(summ.WRITE_SUMMARY_TABLE) >= WRITE_SUMMARY_REFRESH_FREQUENCY:
        print("writing summary tfs")
        tfs_summary = tfs.read(filename) if os.path.exists(filename) else pd.DataFrame()
        
        if "SEED" in tfs_summary.columns.tolist() and "SEED" in summ.WRITE_SUMMARY_TABLE.columns.tolist():
            if max(tfs_summary.SEED) >= min(summ.WRITE_SUMMARY_TABLE.SEED):
                summ.WRITE_SUMMARY_TABLE.SEED += max(tfs_summary.SEED) - min(summ.WRITE_SUMMARY_TABLE.SEED) + 1

        tfs.write(filename,
                  pd.concat([tfs_summary, summ.WRITE_SUMMARY_TABLE], ignore_index=True)
                  )
        summ.WRITE_SUMMARY_TABLE = pd.DataFrame()
        
        
    
    #print("writing summary tfs")
    #tfs_summary = tfs.read(filename) if os.path.exists(filename) else pd.DataFrame()
    #index = len(tfs_summary.index)

    #tfs.write(filename,
    #          pd.concat([tfs_summary, pd.DataFrame(parameters, index=[index])])
    #          )



# ---- summary (helpfull to check the global state of the simulations) -----------------------------
def write_permutation_table(q1_errors, q2_errors, summ):
    if True: #not os.path.exists(f"permutation_s{summ.total():d}.tfs"):
        q1_best_comb = q1_errors.selected_permutation
        q2_best_comb = q2_errors.selected_permutation

        params = {}
        q2_errors.log_strengths(params)
        q1_errors.log_strengths(params)

        keys = [*params.keys(), "DBETX_MIN", "DBETY_MIN", "DBETX_RMS", "DBETY_RMS", "DBETX_MAX", "DBETY_MAX", "PERM_Q1", "PERM_Q2"]

        print("writing permutation tfs")
        tfs_permutation = pd.DataFrame(0.0,index=range(len(q1_errors.permutations)+len(q2_errors.permutations)),columns=keys)

        q2_BBETX, q2_BBETY = q2_errors.get_generated_betabeating_real(with_calibration=False)
        q2_BBETX_wcor, q2_BBETY_wcor = q2_errors.get_generated_betabeating_real(with_calibration=True)
        for i in range(len(q1_errors.permutations)):
            q1_errors.selected_permutation = i

            q1_BBETX, q1_BBETY = q1_errors.get_generated_betabeating_real(with_calibration=False)
            q1_BBETX_wcor, q1_BBETY_wcor = q1_errors.get_generated_betabeating_real(with_calibration=True)

            params = {}
            q1_errors.log_strengths(params)
            q2_errors.log_strengths(params)
            params["DBETX_MIN"] = np.min(q1_BBETX+q2_BBETX)
            params["DBETY_MIN"] = np.min(q1_BBETY+q2_BBETY)
            params["DBETX_RMS"] =    rms(q1_BBETX+q2_BBETX)
            params["DBETY_RMS"] =    rms(q1_BBETY+q2_BBETY)
            params["DBETX_MAX"] = np.max(q1_BBETX+q2_BBETX)
            params["DBETY_MAX"] = np.max(q1_BBETY+q2_BBETY)
            
            params["DBETX_MIN_wcor"] = np.min(q1_BBETX_wcor+q2_BBETX_wcor)
            params["DBETY_MIN_wcor"] = np.min(q1_BBETY_wcor+q2_BBETY_wcor)
            params["DBETX_RMS_wcor"] =    rms(q1_BBETX_wcor+q2_BBETX_wcor)
            params["DBETY_RMS_wcor"] =    rms(q1_BBETY_wcor+q2_BBETY_wcor)
            params["DBETX_MAX_wcor"] = np.max(q1_BBETX_wcor+q2_BBETX_wcor)
            params["DBETY_MAX_wcor"] = np.max(q1_BBETY_wcor+q2_BBETY_wcor)
            
            
            params["PERM_Q1"]   = q1_errors.selected_permutation
            params["PERM_Q2"]   = q2_errors.selected_permutation

            for kk in params.keys():
                tfs_permutation.loc[i,kk]=params[kk]

        q1_errors.selected_permutation = q1_best_comb
        q2_errors.selected_permutation = q2_best_comb

        q1_BBETX, q1_BBETY = q1_errors.get_generated_betabeating_real(with_calibration=False)
        q1_BBETX_wcor, q1_BBETY_wcor = q1_errors.get_generated_betabeating_real(with_calibration=True)
        for i in range(len(q2_errors.permutations)):
            q2_errors.selected_permutation = i

            q2_BBETX, q2_BBETY = q2_errors.get_generated_betabeating_real(with_calibration=False)
            q2_BBETX_wcor, q2_BBETY_wcor = q2_errors.get_generated_betabeating_real(with_calibration=True)

            params = {}
            q1_errors.log_strengths(params)
            q2_errors.log_strengths(params)
            params["DBETX_MIN"] = np.min(q1_BBETX+q2_BBETX)
            params["DBETY_MIN"] = np.min(q1_BBETY+q2_BBETY)
            params["DBETX_RMS"] =    rms(q1_BBETX+q2_BBETX)
            params["DBETY_RMS"] =    rms(q1_BBETY+q2_BBETY)
            params["DBETX_MAX"] = np.max(q1_BBETX+q2_BBETX)
            params["DBETY_MAX"] = np.max(q1_BBETY+q2_BBETY)
            
            params["DBETX_MIN_wcor"] = np.min(q1_BBETX_wcor+q2_BBETX_wcor)
            params["DBETY_MIN_wcor"] = np.min(q1_BBETY_wcor+q2_BBETY_wcor)
            params["DBETX_RMS_wcor"] =    rms(q1_BBETX_wcor+q2_BBETX_wcor)
            params["DBETY_RMS_wcor"] =    rms(q1_BBETY_wcor+q2_BBETY_wcor)
            params["DBETX_MAX_wcor"] = np.max(q1_BBETX_wcor+q2_BBETX_wcor)
            params["DBETY_MAX_wcor"] = np.max(q1_BBETY_wcor+q2_BBETY_wcor)
            
            #params["PERM_Q1"]   = ",".join([str(i) for i in q1_errors.permutations[q1_errors.selected_permutation]])
            #params["PERM_Q2"]   = ",".join([str(i) for i in q2_errors.permutations[q2_errors.selected_permutation]])
            params["PERM_Q1"]   = q1_errors.selected_permutation
            params["PERM_Q2"]   = q2_errors.selected_permutation

            for kk in params.keys():
                tfs_permutation.loc[i+len(q1_errors.permutations),kk]=params[kk]

        q1_errors.selected_permutation = q1_best_comb
        q2_errors.selected_permutation = q2_best_comb
        tfs_permutation['SEED'] = summ.total()
        
        tfs.write(f"permutation_s{summ.total():d}_{NAME_SAVEFILE_FORMAT}.tfs", tfs_permutation)



# ---- summary (helpfull to check the global state of the simulations) -----------------------------
def write_betabeating_table(q1_errors, q2_errors,optic):
    if not os.path.exists("comparison_betabeating.tfs"):
        # compare
        err_twiss = tfs.read(MODEL_TWISS)
        #model_twiss = tfs.read("model1/twiss.dat")
        model_twiss = tfs.read("model1_ft/twiss.dat")
    
        mask_err_monitor   = mask_monitor(  err_twiss)
        mask_model_monitor = mask_monitor(model_twiss)

        q1_BBETX, q1_BBETY = q1_errors.get_generated_betabeating_real(with_calibration=False)
        q2_BBETX, q2_BBETY = q2_errors.get_generated_betabeating_real(with_calibration=False)
    
        params = {}
        params["NAME_MADX"] = err_twiss.NAME[mask_err_monitor]
        params["S"]    = err_twiss.S[mask_err_monitor]
        params["BBETX_MADX"] = err_twiss["BETX"][mask_err_monitor]/model_twiss["BETX"][mask_model_monitor]-1
        params["BBETY_MADX"] = err_twiss["BETY"][mask_err_monitor]/model_twiss["BETY"][mask_model_monitor]-1
    
        params["NAME_TPUG"]  = optic.NAME[mask_monitor(optic)].values
        params["BBETX_TPUG"] = (q1_BBETX+q2_BBETX).values
        params["BBETY_TPUG"] = (q1_BBETY+q2_BBETY).values
        
        print(f'\n\nparams["BBETX_MADX"]=')
        print(params["BBETX_MADX"])
        print(type(params["BBETX_MADX"]))
        print(len(params["BBETX_MADX"]))
        
        print(f'\n\nparams["BBETX_TPUG"]=')
        print(params["BBETX_TPUG"])
        print(type(params["BBETX_TPUG"]))
        print(len(params["BBETX_TPUG"]))
        
        print('\n\nMissing BPM:')
        print([ bpm for bpm in optic.NAME[mask_monitor(optic)] if bpm not in err_twiss.NAME ])
        print(len([ bpm for bpm in optic.NAME[mask_monitor(optic)] if bpm not in err_twiss.NAME ]))
        
        
        print(f'\n\nmodel1_ft/twiss_elements.dat   Qx, Qy = {optic.headers["Q1"]}, {optic.headers["Q2"]}')
        print(f'model1_ft/twiss.dat            Qx, Qy = {optic.headers["Q1"]}, {optic.headers["Q2"]}')
        print(f'{MODEL_TWISS:<25}   Qx, Qy = {optic.headers["Q1"]}, {optic.headers["Q2"]}\n\n')
    
        tfs.write("comparison_betabeating.tfs", pd.DataFrame(params))
        


class Summary():
    def __init__(self) -> None:
        self.not_correctable = 0
        self.low_bbeat = 0
        self.general_failed = 0
        self.passed = 0
        self.previous_study = 0
        self.WRITE_SUMMARY_TABLE = pd.DataFrame()
        self.final_seed = False
        

    def __repr__(self) -> str:
        total = self.total()
        return f"""
Simulation Summary
==================
not correctable : {self.not_correctable:5} ({self.not_correctable/total * 100.0:4.0f}%)
low beta beat   : {self.low_bbeat:5} ({self.low_bbeat/total * 100.0:4.0f}%)
general failed  : {self.general_failed:5} ({self.general_failed/total * 100.0:4.0f}%)
passed          : {self.passed:5} ({self.passed/total * 100.0:4.0f}%)
previous study  : {self.previous_study:5} ({self.previous_study/total * 100.0:4.0f}%)
--------------------------------
total           : {total:5}
"""
    def __str__(self) -> str:
        return self.__repr__()
        
    def total(self) -> int:
        return self.not_correctable + self.general_failed + self.passed + self.low_bbeat + self.previous_study



# ---- ANALYSIS ------------------------------------------------------------------------------------
def do_analysis(summ: Summary, optic_without_error = None):
    """ runs a full analysis, e.g. 
    - runs madx and measures bbeat before and after corrs with initial distribution
    - does the sorting, according to certain criteria
    - runs madx and corrs again after sorting
    """
    
    # Fix the seed random generator per seed
    #np.random.seed(summ.total())
    
#     try:
#         #optic_without_error = tfs.read('model1/twiss_elements.dat')
#         optic_without_error = tfs.read('model1_ft/twiss_elements.dat')
#         #print(optic_without_error)
#         #print(optic_without_error.columns)
#     except:
#         print('WARNING: No twiss file found for the sorting of the magnet. The betas will be set to 1.')
#         optic_without_error = None

    #q2_errors = Q2Pairs(10,2, optic = optic_without_error)
    #q1_errors = Q1Pairs(10,2, optic = optic_without_error)
    
    cerror = uniform(AMP_CALI_ERROR)
    cerror = AMP_CALI_ERROR*cerror/abs(cerror) if (cerror != 0) else AMP_CALI_ERROR
    q2_errors = Q2Pairs(meas_error = AMP_MEAS_ERROR, cali_error = cerror, pres_error = 3*AMP_PRES_ERROR, 
                        stage=STAGE, optic = optic_without_error)
    
    cerror = uniform(AMP_CALI_ERROR)
    cerror = AMP_CALI_ERROR*cerror/abs(cerror) if (cerror != 0) else AMP_CALI_ERROR
    q1_errors = Q1Pairs(meas_error = AMP_MEAS_ERROR, cali_error = cerror, pres_error = 3*AMP_PRES_ERROR, 
                        stage=STAGE, optic = optic_without_error)
    
    cbbetx = pd.concat([q1_errors.monitor_responce_bbetx,q2_errors.monitor_responce_bbetx], axis=1)
    cbbety = pd.concat([q1_errors.monitor_responce_bbety,q2_errors.monitor_responce_bbety], axis=1)

    q1_errors.monitor_responce_bbetx_rms2 = q2_errors.monitor_responce_bbetx_rms2 = cbbetx.T.dot(cbbetx)/len(cbbetx)
    q1_errors.monitor_responce_bbety_rms2 = q2_errors.monitor_responce_bbety_rms2 = cbbety.T.dot(cbbety)/len(cbbety)
    
    if TYPE_SCORE == "RMS_XY":
        #q1_errors.score_def_fn = q1_errors.score_rms_betabeatingXY
        #q2_errors.score_def_fn = q2_errors.score_rms_betabeatingXY
        if FLAG_WITH_calibration:
#             if TYPE_ERROR == "REAL":
#                 q1_errors.score_def_fn = q1_errors.score_rms2_betabeatingXY_wicorrb2_real
#                 q2_errors.score_def_fn = q2_errors.score_rms2_betabeatingXY_wicorrb2_real
#             elif TYPE_ERROR == "MEAS":
            q1_errors.score_def_fn = q1_errors.score_rms2_betabeatingXY_wicorrb2_meas
            q2_errors.score_def_fn = q2_errors.score_rms2_betabeatingXY_wicorrb2_meas
#             else:
#                 raise TypeError(f"TYPE_ERROR for the choise of sorting must be either REAL or MEAS")
        else:
#             if TYPE_ERROR == "REAL":
#                 q1_errors.score_def_fn = q1_errors.score_rms2_betabeatingXY_nocorrb2_real
#                 q2_errors.score_def_fn = q2_errors.score_rms2_betabeatingXY_nocorrb2_real
#             elif TYPE_ERROR == "MEAS":
            q1_errors.score_def_fn = q1_errors.score_rms2_betabeatingXY_nocorrb2_meas
            q2_errors.score_def_fn = q2_errors.score_rms2_betabeatingXY_nocorrb2_meas
#             else:
#                 raise TypeError(f"TYPE_ERROR for the choise of sorting must be either REAL or MEAS")
    elif TYPE_SCORE == "RMS_X":
        q1_errors.score_def_fn = q1_errors.score_rms_betabeatingX
        q2_errors.score_def_fn = q2_errors.score_rms_betabeatingX
    elif TYPE_SCORE == "MAXABS_XY":
        q1_errors.score_def_fn = q1_errors.score_maxabs_betabeatingXY
        q2_errors.score_def_fn = q2_errors.score_maxabs_betabeatingXY
    elif TYPE_SCORE == "MAXABS_X":
        q1_errors.score_def_fn = q1_errors.score_maxabs_betabeatingX
        q2_errors.score_def_fn = q2_errors.score_maxabs_betabeatingX
    else:
        raise TypeError(f'TYPE_SCORE must be either "RMS_XY", "RMS_X", "MAXABS_XY", "MAXABS_X": {TYPE_SCORE=}!')

    # initial distribution, the try ... except block is to remove cases that would either not result
    # in a stable machine anyways or would not be corrected by us because the bbeat is already fine
    fail = False
    #if True:
    try:
        check1x, check1y, err1x, err1y, diff1x, diff1y = do_sim(q1_errors, q2_errors)

        #if err1x < 0.05:
        #    summ.low_bbeat = summ.low_bbeat + 1
        #    print("Low initial beta beating, don't try to improve")
        #    return

        summ.passed = summ.passed + 1

    except CorrectabilityError as e:
        summ.not_correctable = summ.not_correctable + 1 
        print("test case not correctable")
        print(e)
        fail = True
        #return
        
    except Exception as e:
        summ.general_failed = summ.general_failed + 1
        print(e)
        fail = True
        #return

    def sort_and_sim(summ,summ_filename):
        """ applies the sorting, then simulates the sorted lattice, measures beta beating before
        and after corrections
        """

        # ---- do simulations, write results to df ------------------------------------------------
        try:
            check2x, check2y, err2x, err2y, diff2x, diff2y = do_sim(q1_errors, q2_errors)
        except:
            check2x = check2y = err2x = err2y = diff2x = diff2y = 0

        params = {}
        q2_errors.log_strengths(params)
        q1_errors.log_strengths(params)
        params["BBEATX"] = 0 if fail else err1x
        params["BBEATY"] = 0 if fail else err1y
        params["BBEATX_AFTER"] = err2x
        params["BBEATY_AFTER"] = err2y
        params["CORRX"] = 0 if fail else diff1x
        params["CORRY"] = 0 if fail else diff1y
        params["CORRX_AFTER"] = diff2x
        params["CORRY_AFTER"] = diff2y
        params["CHECKX"] = 0 if fail else check1x
        params["CHECKY"] = 0 if fail else check1y
        params["CHECKX_AFTER"] = check2x
        params["CHECKY_AFTER"] = check2y
        # add permutation
        #params["PERM_Q1"] = ",".join([str(i) for i in q1_errors.permutations[q1_errors.selected_permutation]])
        #params["PERM_Q2"] = ",".join([str(i) for i in q2_errors.permutations[q2_errors.selected_permutation]])
        params["PERM_Q1"] = q1_errors.selected_permutation
        params["PERM_Q2"] = q2_errors.selected_permutation
        params["SEED"]    = summ.total()
        
        params = {**params,'init_id_q1':q1_errors.initial_permutation['id'] , 'init_id_q2':q2_errors.initial_permutation['id'] , \
                  **{'init_'+kk:vv for kk,vv in q1_errors.initial_permutation.items() if kk != 'id'} }
        
        params = {**params,'baln_id_q1':q1_errors.best_permutation_alone['id'] , 'baln_id_q2':q2_errors.best_permutation_alone['id'] , \
                  **{'baln_'+kk:vv for kk,vv in q1_errors.best_permutation_alone.items() if kk != 'id'} }
        
        params = {**params,'bbth_id_q1':q1_errors.best_permutation_both['id'] , 'bbth_id_q2':q2_errors.best_permutation_both['id'] , \
                  **{'bbth_'+kk:vv for kk,vv in q1_errors.best_permutation_both.items() if kk != 'id'} }
        
        params = {**params,'wrst_id_q1':q1_errors.worst_permutation_both['id'] , 'wrst_id_q2':q2_errors.worst_permutation_both['id'] , \
                  **{'wrst_'+kk:vv for kk,vv in q1_errors.worst_permutation_both.items() if kk != 'id'} }

        write_summary(params, summ_filename, summ)

        # ---- end of sort_and_sim ----------------------------------------------------------------

    #q1_errors.sort_betabeating()
    #q2_errors.sort_betabeating()
    #sort_and_sim(summ,SUMM_SUMBBEATING)
    try:
        # ---- sorting based on initial beta beating ---------------------------------------------

        def sort_on_bbeat_q1(pairs: Q1Pairs):
            return run_madx_for_sorting(pairs, q2_errors)

        def sort_on_bbeat_q2(pairs: Q2Pairs):
            return run_madx_for_sorting(q1_errors, pairs)
        
        q1_errors.initial_permutation, q2_errors.initial_permutationq = prepare_data_for_analysis(q1_errors, q2_errors)
        
#         (q1_BETX, q1_BETY) = q1_errors.get_generated_betabeating_real(with_calibration=False)
#         (q2_BETX, q2_BETY) = q2_errors.get_generated_betabeating_real(with_calibration=False)
#         BETX = q1_BETX + q2_BETX
#         BETY = q1_BETY + q2_BETY

#         q1_errors.initial_permutation['id'] = q1_errors.selected_permutation
#         q2_errors.initial_permutation['id'] = q2_errors.selected_permutation
        
#         q2_errors.initial_permutation['meanBETX']  = q1_errors.initial_permutation['meanBETX']  = np.average(BETX)
#         q2_errors.initial_permutation['meanBETY']  = q1_errors.initial_permutation['meanBETY']  = np.average(BETY)
        
#         q2_errors.initial_permutation['rmsBETX']  = q1_errors.initial_permutation['rmsBETX']  = rms(BETX)
#         q2_errors.initial_permutation['rmsBETY']  = q1_errors.initial_permutation['rmsBETY']  = rms(BETY)
        
#         q2_errors.initial_permutation['rmsBETXY'] = q1_errors.initial_permutation['rmsBETXY'] = np.sqrt(rms(BETX)**2 + rms(BETY)**2)
        
#         q2_errors.initial_permutation['maxBETX']  = q1_errors.initial_permutation['maxBETX']  = max(abs(BETX))
#         q2_errors.initial_permutation['maxBETY']  = q1_errors.initial_permutation['maxBETY']  = max(abs(BETY))
        
#         q2_errors.initial_permutation['maxBETXY'] = q1_errors.initial_permutation['maxBETXY'] = np.sqrt(max(abs(BETX))**2 + max(abs(BETY))**2)

        
#         (q1_DQX, q1_DQY) = q1_errors.get_generated_tuneshift_real(with_calibration=False)
#         (q2_DQX, q2_DQY) = q2_errors.get_generated_tuneshift_real(with_calibration=False)
#         DQX = q1_DQX + q2_DQX
#         DQY = q1_DQY + q2_DQY
#         q2_errors.initial_permutation['DQX']  = q1_errors.initial_permutation['DQX']  = DQX
#         q2_errors.initial_permutation['DQY']  = q1_errors.initial_permutation['DQY']  = DQY
        
        
#         (q1_BETX, q1_BETY) = q1_errors.get_generated_betabeating_real(with_calibration=True)
#         (q2_BETX, q2_BETY) = q2_errors.get_generated_betabeating_real(with_calibration=True)
#         BETX = q1_BETX + q2_BETX
#         BETY = q1_BETY + q2_BETY
        
#         q2_errors.initial_permutation['meanBETX_wcor']  = q1_errors.initial_permutation['meanBETX_wcor']  = np.average(BETX)
#         q2_errors.initial_permutation['meanBETY_wcor']  = q1_errors.initial_permutation['meanBETY_wcor']  = np.average(BETY)
        
#         q2_errors.initial_permutation['rmsBETX_wcor']  = q1_errors.initial_permutation['rmsBETX_wcor']  = rms(BETX)
#         q2_errors.initial_permutation['rmsBETY_wcor']  = q1_errors.initial_permutation['rmsBETY_wcor']  = rms(BETY)
        
#         q2_errors.initial_permutation['rmsBETXY_wcor'] = q1_errors.initial_permutation['rmsBETXY_wcor'] = np.sqrt(rms(BETX)**2 + rms(BETY)**2)
        
#         q2_errors.initial_permutation['maxBETX_wcor']  = q1_errors.initial_permutation['maxBETX_wcor']  = max(abs(BETX))
#         q2_errors.initial_permutation['maxBETY_wcor']  = q1_errors.initial_permutation['maxBETY_wcor']  = max(abs(BETY))
        
#         q2_errors.initial_permutation['maxBETXY_wcor'] = q1_errors.initial_permutation['maxBETXY_wcor'] = np.sqrt(max(abs(BETX))**2 + max(abs(BETY))**2)

        
#         (q1_DQX, q1_DQY) = q1_errors.get_generated_tuneshift_real(with_calibration=True)
#         (q2_DQX, q2_DQY) = q2_errors.get_generated_tuneshift_real(with_calibration=True)
#         DQX = q1_DQX + q2_DQX
#         DQY = q1_DQY + q2_DQY
#         q2_errors.initial_permutation['DQX_wcor']  = q1_errors.initial_permutation['DQX_wcor']  = DQX
#         q2_errors.initial_permutation['DQY_wcor']  = q1_errors.initial_permutation['DQY_wcor']  = DQY
        
        # q1_errors.sort(sort_on_bbeat_q1)
        # q2_errors.sort(sort_on_bbeat_q2)
        # sort_and_sim(summ,SUMM_BBEAT)

        # # ---- naive sorting methods (but fast, we can crunch through 80k combinations in < 1s ----
        # # ---- this sorts on the difference of errors --------------------------------------------
        # q1_errors.sort_diff()
        # q2_errors.sort_diff()
        # sort_and_sim(summ,SUMM_DIFF)

        # ---- this sorts on the sum of errors (because nearby magnets should self cancel, right? -
        #q1_errors.sort_sum()
        #q2_errors.sort_sum()
        #sort_and_sim(summ,SUMM_SUM)

        # ---- this sorts on the sum of errors weigthed by the betas ------------------------------
        #q1_errors.sort_beta()
        #q2_errors.sort_beta()
        #sort_and_sim(summ,SUMM_SUMBETA)

        # ---- this sorts on the sum of errors (because nearby magnets should self cancel, right? -
        #list_bestscore_q1 = q1_errors.sort_betabeatingXY()
        #list_bestscore_q2 = q2_errors.sort_betabeatingXY()
        list_bestscore_q2 = q2_errors.sort_def_fn() #with_calibration=FLAG_WITH_calibration)
        list_bestscore_q1 = q1_errors.sort_def_fn(other=q2_errors) #,with_calibration=FLAG_WITH_calibration)
        
        q1_errors.best_permutation_alone, q2_errors.best_permutation_alone = prepare_data_for_analysis(q1_errors, q2_errors)
        
#         (q1_BETX, q1_BETY) = q1_errors.get_generated_betabeating_real(with_calibration=False)
#         (q2_BETX, q2_BETY) = q2_errors.get_generated_betabeating_real(with_calibration=False)
#         BETX = q1_BETX + q2_BETX
#         BETY = q1_BETY + q2_BETY

#         q1_errors.best_permutation_alone['id'] = q1_errors.selected_permutation
#         q2_errors.best_permutation_alone['id'] = q2_errors.selected_permutation
        
#         q2_errors.best_permutation_alone['meanBETX']  = q1_errors.best_permutation_alone['meanBETX']  = np.average(BETX)
#         q2_errors.best_permutation_alone['meanBETY']  = q1_errors.best_permutation_alone['meanBETY']  = np.average(BETY)
        
#         q2_errors.best_permutation_alone['rmsBETX']  = q1_errors.best_permutation_alone['rmsBETX']  = rms(BETX)
#         q2_errors.best_permutation_alone['rmsBETY']  = q1_errors.best_permutation_alone['rmsBETY']  = rms(BETY)
        
#         q2_errors.best_permutation_alone['rmsBETXY'] = q1_errors.best_permutation_alone['rmsBETXY'] = np.sqrt(rms(BETX)**2 + rms(BETY)**2)
        
#         q2_errors.best_permutation_alone['maxBETX']  = q1_errors.best_permutation_alone['maxBETX']  = max(abs(BETX))
#         q2_errors.best_permutation_alone['maxBETY']  = q1_errors.best_permutation_alone['maxBETY']  = max(abs(BETY))
        
#         q2_errors.best_permutation_alone['maxBETXY'] = q1_errors.best_permutation_alone['maxBETXY'] = np.sqrt(max(abs(BETX))**2 + max(abs(BETY))**2)
        
#         (q1_DQX, q1_DQY) = q1_errors.get_generated_tuneshift_real(with_calibration=False)
#         (q2_DQX, q2_DQY) = q2_errors.get_generated_tuneshift_real(with_calibration=False)
#         DQX = q1_DQX + q2_DQX
#         DQY = q1_DQY + q2_DQY
#         q2_errors.best_permutation_alone['DQX']  = q1_errors.best_permutation_alone['DQX']  = DQX
#         q2_errors.best_permutation_alone['DQY']  = q1_errors.best_permutation_alone['DQY']  = DQY

        
#         (q1_BETX, q1_BETY) = q1_errors.get_generated_betabeating_real(with_calibration=True)
#         (q2_BETX, q2_BETY) = q2_errors.get_generated_betabeating_real(with_calibration=True)
#         BETX = q1_BETX + q2_BETX
#         BETY = q1_BETY + q2_BETY
        
#         q2_errors.best_permutation_alone['meanBETX_wcor']  = q1_errors.best_permutation_alone['meanBETX_wcor']  = np.average(BETX)
#         q2_errors.best_permutation_alone['meanBETY_wcor']  = q1_errors.best_permutation_alone['meanBETY_wcor']  = np.average(BETY)
        
#         q2_errors.best_permutation_alone['rmsBETX_wcor']  = q1_errors.best_permutation_alone['rmsBETX_wcor']  = rms(BETX)
#         q2_errors.best_permutation_alone['rmsBETY_wcor']  = q1_errors.best_permutation_alone['rmsBETY_wcor']  = rms(BETY)
        
#         q2_errors.best_permutation_alone['rmsBETXY_wcor'] = q1_errors.best_permutation_alone['rmsBETXY_wcor'] = np.sqrt(rms(BETX)**2 + rms(BETY)**2)
        
#         q2_errors.best_permutation_alone['maxBETX_wcor']  = q1_errors.best_permutation_alone['maxBETX_wcor']  = max(abs(BETX))
#         q2_errors.best_permutation_alone['maxBETY_wcor']  = q1_errors.best_permutation_alone['maxBETY_wcor']  = max(abs(BETY))
        
#         q2_errors.best_permutation_alone['maxBETXY_wcor'] = q1_errors.best_permutation_alone['maxBETXY_wcor'] = np.sqrt(max(abs(BETX))**2 + max(abs(BETY))**2)
        
#         (q1_DQX, q1_DQY) = q1_errors.get_generated_tuneshift_real(with_calibration=True)
#         (q2_DQX, q2_DQY) = q2_errors.get_generated_tuneshift_real(with_calibration=True)
#         DQX = q1_DQX + q2_DQX
#         DQY = q1_DQY + q2_DQY
#         q2_errors.best_permutation_alone['DQX_wcor']  = q1_errors.best_permutation_alone['DQX_wcor']  = DQX
#         q2_errors.best_permutation_alone['DQY_wcor']  = q1_errors.best_permutation_alone['DQY_wcor']  = DQY

        #sort_and_sim(summ,SUMM_SUMBBEATING)

        # ---- this sorts on the sum of errors (because nearby magnets should self cancel, right? -
        q1_errors, q2_errors = sort_betabeating_v2(q1_errors, list_bestscore_q1[:NUMB_PERMUT_CROSSSCORING], 
                                                   q2_errors, list_bestscore_q2[:NUMB_PERMUT_CROSSSCORING])

        # ---- this sorts on the sum of errors (because nearby magnets should self cancel, right? -
        q1_errors, q2_errors, list_wrstscore_q1, list_wrstscore_q2 = sort_betabeating_worst_v2(q1_errors, list_bestscore_q1[-NUMB_PERMUT_CROSSSCORING:], 
                                                         q2_errors, list_bestscore_q2[-NUMB_PERMUT_CROSSSCORING:])
        sort_and_sim(summ,SUMM_SUMBBEATING_V2)

    #except:
    #    return
    #except Exception as e:
    except CorrectabilityError as e:
        print('\n# '+'='*60)
        print('#   * BUG IN SORTING')
        print(e)
        print('# '+'='*60+'\n')
        return
    
    #write_betabeating_table(q1_errors, q2_errors,optic_without_error)
    if summ.total() <=1:
        write_permutation_table(q1_errors, q2_errors, summ)
    #if summ.total() <=100:
    #    write_permutation_table(q1_errors, q2_errors, summ)
    #if summ.total() < 2:
    #    scan_limit_MADX(q1_errors, q2_errors, list_wrstscore_q1, list_wrstscore_q2)



def generate_optic():
    if not os.path.exists("model1_ft/hllhc_lhcb1.seq"): 
        # remove twiss output, its existance after running madx indicates success
        system(f"cp job.hl16_nominal.madx model{MODEL_NUMBER}/run_job.hl16_nominal.madx")

        # remove twiss output, its existance after running madx indicates success
        system(f"rm {MODEL_TWISS}")
        system(f"touch model{MODEL_NUMBER}/knobs.madx")

        with open(os.devnull, "w") as devnull:
            #Popen([MADX, "run_job.inj.madx"], stdout=devnull, cwd="model").wait()
            Popen([MADX, "run_job.hl16_nominal.madx"], stdout=devnull, cwd=f"model{MODEL_NUMBER}").wait()
        system(f"rm model{MODEL_NUMBER}*.*")
        system(f"touch model{MODEL_NUMBER}/knobs.madx")

    if not os.path.exists("model1_ft/acc-models-lhc") or not os.path.exists("model1_ft/macros"): 
        cwd = os.getcwd()
        os.symlink(Path(cwd,'model/acc-models-lhc'), "model1_ft/acc-models-lhc")
        os.symlink(Path(cwd,'macros'), "model1_ft/macros")
        system(f"touch model1_ft/knobs.madx")


            
def run_madx(q1_errors: Q1Pairs, q2_errors: Q2Pairs):
    print("running madx")

    q2_errors.write_errors_to_file(f"model{MODEL_NUMBER}/errors_Q2.madx")
    q1_errors.write_errors_to_file(f"model{MODEL_NUMBER}/errors_Q1.madx")

    if not q1_errors.check_correctability():
        raise CorrectabilityError([x.real_error for x in q1_errors.cold_masses])
    if not q2_errors.check_correctability():
        raise CorrectabilityError([x.real_error for x in q2_errors.cold_masses])

    #template_file = open("job.inj.madx", "r")
    #jobfile = open(f"model{MODEL_NUMBER}/run_job.inj.madx", "w")
    template_file = open("job.hl16_nominal_thin.madx", "r")
    jobfile = open(f"model{MODEL_NUMBER}/run_job.hl16_nominal_thin.madx", "w")
    jobfile.write(template_file.read().replace("_TRACK_", "0")) # no tracking
    template_file.close()
    jobfile.close()

    # remove twiss output, its existance after running madx indicates success
    system(f"rm {MODEL_TWISS}")

    with open(os.devnull, "w") as devnull:
        #Popen([MADX, "run_job.inj.madx"], stdout=devnull, cwd="model").wait()
        Popen([MADX, "run_job.hl16_nominal_thin.madx"], stdout=devnull, cwd=f"model{MODEL_NUMBER}").wait()

    if not os.path.exists(f"{MODEL_TWISS}"):
        raise RuntimeError("twiss_err_b1.tfs does not exist")



def run_madx_for_sorting(q1_errors: Q1Pairs, q2_errors: Q2Pairs) -> float:

    run_madx(q1_errors, q2_errors)

    # error beta beating
    err_twiss = tfs.read(MODEL_TWISS)
    #model_twiss = tfs.read("model1/twiss.dat")
    model_twiss = tfs.read("model1_ft/twiss.dat")
    
    mask_err_monitor   = mask_monitor(  err_twiss)
    mask_model_monitor = mask_monitor(model_twiss)
    return rms(err_twiss["BETX"][mask_err_monitor]/model_twiss["BETX"][mask_model_monitor]-1), rms(err_twiss["BETY"][mask_err_monitor]/model_twiss["BETY"][mask_model_monitor]-1)



def do_sim(q1_errors: Q1Pairs, q2_errors: Q2Pairs) -> Tuple[float, float, float]:
    """
    runs madx and simulates corrections
    returns:
        (checkx, checky, errx, erry, diffx, diffy)
        - err = bbeat of the virgin lattice
        - check = bbeat reconstruction with the result from global_corrections
        - diff = difference between the two

    Note:
        One could (should) apply the corrections (= negative of the reconstructed errors) to the
        virgin lattice, but this is not implemented yet.
    """
    print("running simulation")

    if DO_MADX:
        run_madx(q1_errors, q2_errors)

    if MODEL_NUMBER:
        if not os.path.exists(f"model1_ft{MODEL_NUMBER}"):
            os.makedirs(f"model1_ft{MODEL_NUMBER}/")
            cwd = os.getcwd()
            os.symlink(Path(cwd,'model/acc-models-lhc'), f"model1_ft{MODEL_NUMBER}/acc-models-lhc")
            os.symlink(Path(cwd,'macros'), f"model1_ft{MODEL_NUMBER}/macros")
            system(f"touch model1_ft{MODEL_NUMBER}/knobs.madx")
        system(f"cp model1_ft/*.* model1_ft{MODEL_NUMBER}/")
            

    accel_params = dict(
        accel="lhc",
        year="hl16", # to be checked
        beam=1,
        #model_dir=Path("model1").absolute(),
        model_dir=Path(f"model1_ft{MODEL_NUMBER}").absolute(),
            )


    if DO_FR:
        print("creating response entrypoint")
        create_response_entrypoint(
            #outfile_path=Path("model1/FullResponse.h5"),
            outfile_path=Path(f"model1_ft{MODEL_NUMBER}/FullResponse.h5"),
            creator="madx",
            optics_params=['PHASEX', 'PHASEY', 'BETX', 'BETY', 'Q'],
            variable_categories=VAR_CATS,
            delta_k=1.0e-7,
            **accel_params 
        )


    # fake measurement
    fake_measurement(twiss=MODEL_TWISS,
                     #model="model1/twiss.dat",
                     model=f"model1_ft{MODEL_NUMBER}/twiss.dat",
                     outputdir=f"fake_measurements{MODEL_NUMBER}")


    if DO_CORR:
        print("running global correction")
        global_correction.global_correction_entrypoint(
            meas_dir=f"fake_measurements{MODEL_NUMBER}",
            output_dir=f"global_corrections{MODEL_NUMBER}",
            #fullresponse_path="model1/FullResponse.h5",
            fullresponse_path=f"model1_ft{MODEL_NUMBER}/FullResponse.h5",
            iterations=1,
            optics_params=['PHASEX', 'PHASEY', 'BETX', 'BETY', 'Q'],
            variable_categories=VAR_CATS,
            **accel_params,
                )
        #system("cp scripts_global_corrs/rerun.job.madx global_corrections/rerun.job.madx")
        system(f"cp scripts_global_corrs/rerun.job.hl16_nominal_thin.madx global_corrections{MODEL_NUMBER}/rerun.job.madx")
        if MODEL_NUMBER:
            system(f"mv global_corrections{MODEL_NUMBER}/rerun.job.madx global_corrections{MODEL_NUMBER}/rerun_ref.job.madx")
            template_file = open(f"global_corrections{MODEL_NUMBER}/rerun_ref.job.madx", "r")
            jobfile = open(f"global_corrections{MODEL_NUMBER}/rerun.job.madx", "w")
            jobfile.write(template_file.read().replace("../global_corrections/", f"../global_corrections{MODEL_NUMBER}/")) # no tracking
            template_file.close()
            jobfile.close()
        with open(os.devnull, "w") as devnull:
            Popen([MADX, "rerun.job.madx"], cwd=f"global_corrections{MODEL_NUMBER}", stdout=devnull).wait()




    # compare
    check_twiss = tfs.read(f"global_corrections{MODEL_NUMBER}/twiss_global_b1.tfs")
    err_twiss   = tfs.read(MODEL_TWISS)
    #model_twiss = tfs.read("model1/twiss.dat")
    model_twiss = tfs.read(f"model1_ft{MODEL_NUMBER}/twiss.dat")
    
    mask_check_monitor = mask_monitor(check_twiss)
    mask_err_monitor   = mask_monitor(  err_twiss)
    mask_model_monitor = mask_monitor(model_twiss)

    #plt.plot(check_twiss["S"], check_twiss["BETX"]/model_twiss["BETX"]-1, label="check")
    #plt.plot(err_twiss["S"], err_twiss["BETX"]/model_twiss["BETX"]-1, label="err")
    #plt.legend()

    #plt.show()

    # reconstructed beta beating
    rms_checkx = rms(check_twiss["BETX"][mask_check_monitor]/model_twiss["BETX"][mask_model_monitor]-1)
    rms_checky = rms(check_twiss["BETY"][mask_check_monitor]/model_twiss["BETY"][mask_model_monitor]-1)

    # error beta beating
    rms_errx = rms(err_twiss["BETX"][mask_err_monitor]/model_twiss["BETX"][mask_model_monitor]-1)
    rms_erry = rms(err_twiss["BETY"][mask_err_monitor]/model_twiss["BETY"][mask_model_monitor]-1)

    # difference (how close do we get?)
    rms_diffx = rms(check_twiss["BETX"][mask_check_monitor]/err_twiss["BETX"][mask_err_monitor]-1)
    rms_diffy = rms(check_twiss["BETY"][mask_check_monitor]/err_twiss["BETY"][mask_err_monitor]-1)

    print(f"rms check: ({rms_checkx},{rms_checky}), rms err: ({rms_errx},{rms_erry})")
    print(f"rms diff: ({rms_diffx},{rms_diffy})")

    return rms_checkx, rms_checky, rms_errx, rms_erry, rms_diffx, rms_diffy
    
    
    
def sort_betabeating_v2(q1_errors: Q1Pairs, list_q1: np.array, q2_errors: Q2Pairs, list_q2: np.array):
    len_lq1 = len(list_q1)
    len_lq2 = len(list_q2)
    
    print(f"sorting both group of magnet at the same time") 
    print(f"searching for best combination in {len_lq1*len_lq2} permutations")
    print(f"using the diff method")

    #score = q1_errors.score_rms_betabeatingXY(q1_errors, q2_errors)
    score = q1_errors.score_def_fn(q1_errors, q2_errors) #, with_calibration=FLAG_WITH_calibration)
    print(f" -- initial score: {score}")

    next_progress = 0.1
    best_comb = ( q1_errors.selected_permutation , q2_errors.selected_permutation )
        
    for i1 in range(len_lq1):
        for i2 in range(len_lq2):
            if ( i1*len_lq2 + i2 ) / ( len_lq1*len_lq2 ) > next_progress:
                print(f"progress: {100* ( i1*len_lq2 + i2 ) / ( len_lq1*len_lq2 ):.0f} %")
                next_progress += 0.1

            q1_errors.selected_permutation = int(list_q1[i1])
            q2_errors.selected_permutation = int(list_q2[i2])

            #sum = q1_errors.score_rms_betabeatingXY(q1_errors, q2_errors)
            sum = q1_errors.score_def_fn(q1_errors, q2_errors) #, with_calibration=FLAG_WITH_calibration)

            if sum < score:
                score = sum
                best_comb = ( q1_errors.selected_permutation , q2_errors.selected_permutation )
                
    q1_errors.selected_permutation = best_comb[0]
    q2_errors.selected_permutation = best_comb[1]
    
    q1_errors.best_permutation_both, q2_errors.best_permutation_both = prepare_data_for_analysis(q1_errors, q2_errors)
    
#     (q1_BETX, q1_BETY) = q1_errors.get_generated_betabeating_real(with_calibration=False)
#     (q2_BETX, q2_BETY) = q2_errors.get_generated_betabeating_real(with_calibration=False)
#     BETX = q1_BETX + q2_BETX
#     BETY = q1_BETY + q2_BETY

#     q1_errors.best_permutation_both['id'] = q1_errors.selected_permutation
#     q2_errors.best_permutation_both['id'] = q1_errors.selected_permutation
    
#     q2_errors.best_permutation_both['meanBETX']  = q1_errors.best_permutation_both['meanBETX']  = np.average(BETX)
#     q2_errors.best_permutation_both['meanBETY']  = q1_errors.best_permutation_both['meanBETY']  = np.average(BETY)
    
#     q2_errors.best_permutation_both['rmsBETX']  = q1_errors.best_permutation_both['rmsBETX']  = rms(BETX)
#     q2_errors.best_permutation_both['rmsBETY']  = q1_errors.best_permutation_both['rmsBETY']  = rms(BETY)
    
#     q2_errors.best_permutation_both['rmsBETXY'] = q1_errors.best_permutation_both['rmsBETXY'] = np.sqrt(rms(BETX)**2 + rms(BETY)**2)
    
#     q2_errors.best_permutation_both['maxBETX']  = q1_errors.best_permutation_both['maxBETX']  = max(abs(BETX))
#     q2_errors.best_permutation_both['maxBETY']  = q1_errors.best_permutation_both['maxBETY']  = max(abs(BETY))
    
#     q2_errors.best_permutation_both['maxBETXY'] = q1_errors.best_permutation_both['maxBETXY'] = np.sqrt(max(abs(BETX))**2 + max(abs(BETY))**2)

#     (q1_DQX, q1_DQY) = q1_errors.get_generated_tuneshift_real(with_calibration=False)
#     (q2_DQX, q2_DQY) = q2_errors.get_generated_tuneshift_real(with_calibration=False)
#     DQX = q1_DQX + q2_DQX
#     DQY = q1_DQY + q2_DQY
#     q2_errors.best_permutation_both['DQX']  = q1_errors.best_permutation_both['DQX']  = DQX
#     q2_errors.best_permutation_both['DQY']  = q1_errors.best_permutation_both['DQY']  = DQY
    
#     (q1_BETX, q1_BETY) = q1_errors.get_generated_betabeating_real(with_calibration=True)
#     (q2_BETX, q2_BETY) = q2_errors.get_generated_betabeating_real(with_calibration=True)
#     BETX = q1_BETX + q2_BETX
#     BETY = q1_BETY + q2_BETY
    
#     q2_errors.best_permutation_both['meanBETX_wcor']  = q1_errors.best_permutation_both['meanBETX_wcor']  = np.average(BETX)
#     q2_errors.best_permutation_both['meanBETY_wcor']  = q1_errors.best_permutation_both['meanBETY_wcor']  = np.average(BETY)
    
#     q2_errors.best_permutation_both['rmsBETX_wcor']  = q1_errors.best_permutation_both['rmsBETX_wcor']  = rms(BETX)
#     q2_errors.best_permutation_both['rmsBETY_wcor']  = q1_errors.best_permutation_both['rmsBETY_wcor']  = rms(BETY)
    
#     q2_errors.best_permutation_both['rmsBETXY_wcor'] = q1_errors.best_permutation_both['rmsBETXY_wcor'] = np.sqrt(rms(BETX)**2 + rms(BETY)**2)
    
#     q2_errors.best_permutation_both['maxBETX_wcor']  = q1_errors.best_permutation_both['maxBETX_wcor']  = max(abs(BETX))
#     q2_errors.best_permutation_both['maxBETY_wcor']  = q1_errors.best_permutation_both['maxBETY_wcor']  = max(abs(BETY))
    
#     q2_errors.best_permutation_both['maxBETXY_wcor'] = q1_errors.best_permutation_both['maxBETXY_wcor'] = np.sqrt(max(abs(BETX))**2 + max(abs(BETY))**2)
    
#     (q1_DQX, q1_DQY) = q1_errors.get_generated_tuneshift_real(with_calibration=True)
#     (q2_DQX, q2_DQY) = q2_errors.get_generated_tuneshift_real(with_calibration=True)
#     DQX = q1_DQX + q2_DQX
#     DQY = q1_DQY + q2_DQY
#     q2_errors.best_permutation_both['DQX_wcor']  = q1_errors.best_permutation_both['DQX_wcor']  = DQX
#     q2_errors.best_permutation_both['DQY_wcor']  = q1_errors.best_permutation_both['DQY_wcor']  = DQY
 
    print(f"final score: {score}")
    return q1_errors, q2_errors
    
    
    
def sort_betabeating_worst_v2(q1_errors: Q1Pairs, list_q1: np.array, q2_errors: Q2Pairs, list_q2: np.array):
    len_lq1 = len(list_q1)
    len_lq2 = len(list_q2)
    
    print(f"sorting both group of magnet at the same time") 
    print(f"searching for worst combination in {len_lq1*len_lq2} permutations")
    print(f"using the diff method")

    #score = q1_errors.score_rms_betabeatingXY(q1_errors, q2_errors)
    score = q1_errors.score_def_fn(q1_errors, q2_errors) #, with_calibration=False)
    print(f" -- initial score: {score}")

    next_progress = 0.1
    best_comb  = ( q1_errors.selected_permutation , q2_errors.selected_permutation )
    worst_comb = ( q1_errors.selected_permutation , q2_errors.selected_permutation )
    
    list_score = np.empty( len_lq1 * len_lq2 )
    list_idxq1 = np.empty( len_lq1 * len_lq2 )
    list_idxq2 = np.empty( len_lq1 * len_lq2 )
        
    for i1 in range(len_lq1):
        for i2 in range(len_lq2):
            ii = i1*len_lq2 + i2
            if ( ii ) / ( len_lq1*len_lq2 ) > next_progress:
                print(f"progress: {100* ( ii ) / ( len_lq1*len_lq2 ):.0f} %")
                next_progress += 0.1

            q1_errors.selected_permutation = int(list_q1[i1])
            q2_errors.selected_permutation = int(list_q2[i2])

            #sum = q1_errors.score_rms_betabeatingXY(q1_errors, q2_errors)
            sum = q1_errors.score_def_fn(q1_errors, q2_errors) #, with_calibration=False)
            
            list_score[ii] = sum
            list_idxq1[ii] = q1_errors.selected_permutation
            list_idxq2[ii] = q2_errors.selected_permutation

            if sum > score:
                score = sum
                worst_comb = ( q1_errors.selected_permutation , q2_errors.selected_permutation )
                
    q1_errors.selected_permutation = worst_comb[0]
    q2_errors.selected_permutation = worst_comb[1]
    
    q1_errors.worst_permutation_both, q2_errors.worst_permutation_both = prepare_data_for_analysis(q1_errors, q2_errors)
    
#     (q1_BETX, q1_BETY) = q1_errors.get_generated_betabeating_real(with_calibration=False)
#     (q2_BETX, q2_BETY) = q2_errors.get_generated_betabeating_real(with_calibration=False)
#     BETX = q1_BETX + q2_BETX
#     BETY = q1_BETY + q2_BETY

#     q1_errors.worst_permutation_both['id'] = q1_errors.selected_permutation
#     q2_errors.worst_permutation_both['id'] = q2_errors.selected_permutation
    
#     q2_errors.worst_permutation_both['meanBETX'] = q1_errors.worst_permutation_both['meanBETX'] = np.average(BETX)
#     q2_errors.worst_permutation_both['meanBETY'] = q1_errors.worst_permutation_both['meanBETY'] = np.average(BETY)

#     q2_errors.worst_permutation_both['rmsBETX']  = q1_errors.worst_permutation_both['rmsBETX']  = rms(BETX)
#     q2_errors.worst_permutation_both['rmsBETY']  = q1_errors.worst_permutation_both['rmsBETY']  = rms(BETY)
    
#     q2_errors.worst_permutation_both['rmsBETXY'] = q1_errors.worst_permutation_both['rmsBETXY'] = np.sqrt(rms(BETX)**2 + rms(BETY)**2)
    
#     q2_errors.worst_permutation_both['maxBETX']  = q1_errors.worst_permutation_both['maxBETX']  = max(abs(BETX))
#     q2_errors.worst_permutation_both['maxBETY']  = q1_errors.worst_permutation_both['maxBETY']  = max(abs(BETY))
    
#     q2_errors.worst_permutation_both['maxBETXY'] = q1_errors.worst_permutation_both['maxBETXY'] = np.sqrt(max(abs(BETX))**2 + max(abs(BETY))**2)
    
#     (q1_DQX, q1_DQY) = q1_errors.get_generated_tuneshift_real(with_calibration=False)
#     (q2_DQX, q2_DQY) = q2_errors.get_generated_tuneshift_real(with_calibration=False)
#     DQX = q1_DQX + q2_DQX
#     DQY = q1_DQY + q2_DQY
#     q2_errors.worst_permutation_both['DQX']  = q1_errors.worst_permutation_both['DQX']  = DQX
#     q2_errors.worst_permutation_both['DQY']  = q1_errors.worst_permutation_both['DQY']  = DQY
 
    
#     (q1_BETX, q1_BETY) = q1_errors.get_generated_betabeating_real(with_calibration=True)
#     (q2_BETX, q2_BETY) = q2_errors.get_generated_betabeating_real(with_calibration=True)
#     BETX = q1_BETX + q2_BETX
#     BETY = q1_BETY + q2_BETY
    
#     q2_errors.worst_permutation_both['meanBETX_wcor'] = q1_errors.worst_permutation_both['meanBETX_wcor'] = np.average(BETX)
#     q2_errors.worst_permutation_both['meanBETY_wcor'] = q1_errors.worst_permutation_both['meanBETY_wcor'] = np.average(BETY)

#     q2_errors.worst_permutation_both['rmsBETX_wcor']  = q1_errors.worst_permutation_both['rmsBETX_wcor']  = rms(BETX)
#     q2_errors.worst_permutation_both['rmsBETY_wcor']  = q1_errors.worst_permutation_both['rmsBETY_wcor']  = rms(BETY)
    
#     q2_errors.worst_permutation_both['rmsBETXY_wcor'] = q1_errors.worst_permutation_both['rmsBETXY_wcor'] = np.sqrt(rms(BETX)**2 + rms(BETY)**2)
    
#     q2_errors.worst_permutation_both['maxBETX_wcor']  = q1_errors.worst_permutation_both['maxBETX_wcor']  = max(abs(BETX))
#     q2_errors.worst_permutation_both['maxBETY_wcor']  = q1_errors.worst_permutation_both['maxBETY_wcor']  = max(abs(BETY))
    
#     q2_errors.worst_permutation_both['maxBETXY_wcor'] = q1_errors.worst_permutation_both['maxBETXY_wcor'] = np.sqrt(max(abs(BETX))**2 + max(abs(BETY))**2)
            
#     (q1_DQX, q1_DQY) = q1_errors.get_generated_tuneshift_real(with_calibration=True)
#     (q2_DQX, q2_DQY) = q2_errors.get_generated_tuneshift_real(with_calibration=True)
#     DQX = q1_DQX + q2_DQX
#     DQY = q1_DQY + q2_DQY
#     q2_errors.worst_permutation_both['DQX_wcor']  = q1_errors.worst_permutation_both['DQX_wcor']  = DQX
#     q2_errors.worst_permutation_both['DQY_wcor']  = q1_errors.worst_permutation_both['DQY_wcor']  = DQY
 
    q1_errors.selected_permutation = best_comb[0]
    q2_errors.selected_permutation = best_comb[1]
                
    print(f"final score: {score}")
    return q1_errors, q2_errors, list_idxq1[np.argsort(list_score)], list_idxq2[np.argsort(list_score)]



# def scan_limit_MADX(q1_errors, q2_errors, list_wrstscore_q1, list_wrstscore_q2):
#     list_idxq1 = list_wrstscore_q1[-1::-1]
#     list_idxq2 = list_wrstscore_q2[-1::-1]
#     next_progress = 0.1
#     best_comb  = ( q1_errors.selected_permutation , q2_errors.selected_permutation )
    
#     nb_good_run = 0
#     border_ii = ii = 0
#     #while nb_good_run < 5 and ii < len(list_idxq1):
#     for ii in range(len(list_idxq1)):
#         print(f"{ii} = ({list_idxq1[ii]} , {list_idxq2[ii]})")
#         q1_errors.selected_permutation = int(list_idxq1[ii])
#         q2_errors.selected_permutation = int(list_idxq2[ii])
        
#         status = 0.0
#         try:
#             do_sim(q1_errors, q2_errors)
#             nb_good_run += 1
#             status = +1.0
#             print(f"{ii} = ({q1_errors.selected_permutation} , {q2_errors.selected_permutation})    ->   GOOD!")
        
#         #except CorrectabilityError as e:
#         #    nb_good_run = 0
#         except Exception as e:
#             nb_good_run = 0
#             border_ii = ii
#             status = -1.0
#             print(e)
#             print(f"{ii} = ({q1_errors.selected_permutation} , {q2_errors.selected_permutation})  ->   BAD!")
            
#         (q1_BETX, q1_BETY) = q1_errors.get_generated_betabeating_real()
#         (q2_BETX, q2_BETY) = q2_errors.get_generated_betabeating_real()
#         BETX = q1_BETX + q2_BETX
#         BETY = q1_BETY + q2_BETY
    
#         params = {}
#         q2_errors.log_strengths(params)
#         q1_errors.log_strengths(params)
    
#         params['meanBETX'] = np.average(BETX)
#         params['meanBETY'] = np.average(BETY)

    
#         params['rmsBETX'] = rms(BETX)
#         params['rmsBETY'] = rms(BETY)
    
#         params['maxBETX'] = max(abs(BETX))
#         params['maxBETY'] = max(abs(BETY))
#         params['status']  = status
    
#         write_summary(params, 'summ_beta_wall_madx.tfs')

#     q1_errors.selected_permutation , q2_errors.selected_permutation = best_comb
#     return q1_errors, q2_errors



def rms(array):
    """ root mean square """
    return np.sqrt(np.mean(array**2))



def prepare_data_for_analysis(q1_errors, q2_errors):
    
    analysis_data_q1 = {'id': q1_errors.selected_permutation}
    analysis_data_q2 = {'id': q2_errors.selected_permutation}
    
    #(q1_BETX, q1_BETY) = q1_errors.get_generated_betabeating_meas(with_calibration=False)
    #(q2_BETX, q2_BETY) = q2_errors.get_generated_betabeating_meas(with_calibration=False)
    #BETX = q1_BETX + q2_BETX
    #BETY = q1_BETY + q2_BETY
    
    #analysis_data_q2['meanBETX'] = analysis_data_q1['meanBETX'] = np.average(BETX)
    #analysis_data_q2['meanBETY'] = analysis_data_q1['meanBETY'] = np.average(BETY)

    #analysis_data_q2['rmsBETX']  = analysis_data_q1['rmsBETX']  = rms(BETX)
    #analysis_data_q2['rmsBETY']  = analysis_data_q1['rmsBETY']  = rms(BETY)
    
    #analysis_data_q2['rmsBETXY'] = analysis_data_q1['rmsBETXY'] = np.sqrt(rms(BETX)**2 + rms(BETY)**2)
    
    #analysis_data_q2['maxBETX']  = analysis_data_q1['maxBETX']  = max(abs(BETX))
    #analysis_data_q2['maxBETY']  = analysis_data_q1['maxBETY']  = max(abs(BETY))
    
    ##analysis_data_q2['maxBETXY'] = analysis_data_q1['maxBETXY'] = np.sqrt(max(abs(BETX))**2 + max(abs(BETY))**2)
    
    #(q1_DQX, q1_DQY) = q1_errors.get_generated_tuneshift_meas(with_calibration=False)
    #(q2_DQX, q2_DQY) = q2_errors.get_generated_tuneshift_meas(with_calibration=False)
    #DQX = q1_DQX + q2_DQX
    #DQY = q1_DQY + q2_DQY
    #analysis_data_q2['DQX']  = analysis_data_q1['DQX']  = DQX
    #analysis_data_q2['DQY']  = analysis_data_q1['DQY']  = DQY
 
    
    #(q1_BETX, q1_BETY) = q1_errors.get_generated_betabeating_meas(with_calibration=True)
    #(q2_BETX, q2_BETY) = q2_errors.get_generated_betabeating_meas(with_calibration=True)
    #BETX = q1_BETX + q2_BETX
    #BETY = q1_BETY + q2_BETY
    
    #analysis_data_q2['meanBETX_wcor'] = analysis_data_q1['meanBETX_wcor'] = np.average(BETX)
    #analysis_data_q2['meanBETY_wcor'] = analysis_data_q1['meanBETY_wcor'] = np.average(BETY)

    #analysis_data_q2['rmsBETX_wcor']  = analysis_data_q1['rmsBETX_wcor']  = rms(BETX)
    #analysis_data_q2['rmsBETY_wcor']  = analysis_data_q1['rmsBETY_wcor']  = rms(BETY)
    
    #analysis_data_q2['rmsBETXY_wcor'] = analysis_data_q1['rmsBETXY_wcor'] = np.sqrt(rms(BETX)**2 + rms(BETY)**2)
    
    #analysis_data_q2['maxBETX_wcor']  = analysis_data_q1['maxBETX_wcor']  = max(abs(BETX))
    #analysis_data_q2['maxBETY_wcor']  = analysis_data_q1['maxBETY_wcor']  = max(abs(BETY))
    
    ##analysis_data_q2['maxBETXY_wcor'] = analysis_data_q1['maxBETXY_wcor'] = np.sqrt(max(abs(BETX))**2 + max(abs(BETY))**2)
            
    #(q1_DQX, q1_DQY) = q1_errors.get_generated_tuneshift_meas(with_calibration=True)
    #(q2_DQX, q2_DQY) = q2_errors.get_generated_tuneshift_meas(with_calibration=True)
    #DQX = q1_DQX + q2_DQX
    #DQY = q1_DQY + q2_DQY
    #analysis_data_q2['DQX_wcor']  = analysis_data_q1['DQX_wcor']  = DQX
    #analysis_data_q2['DQY_wcor']  = analysis_data_q1['DQY_wcor']  = DQY
    
    
    for lab,withcal in zip([0,1],[False,True]):
        (q1_BETX, q1_BETY) = q1_errors.get_generated_betabeating_meas(with_calibration=withcal)
        (q2_BETX, q2_BETY) = q2_errors.get_generated_betabeating_meas(with_calibration=withcal)
        BETX = q1_BETX + q2_BETX
        BETY = q1_BETY + q2_BETY

        (q1_DQX, q1_DQY) = q1_errors.get_generated_tuneshift_meas(with_calibration=withcal)
        (q2_DQX, q2_DQY) = q2_errors.get_generated_tuneshift_meas(with_calibration=withcal)
        DQX = q1_DQX + q2_DQX
        DQY = q1_DQY + q2_DQY

        analysis_data_q2[f'meanBETX_TEmeas_TC{lab}'] = analysis_data_q1[f'meanBETX_TEmeas_TC{lab}'] = np.average(BETX)
        analysis_data_q2[f'meanBETY_TEmeas_TC{lab}'] = analysis_data_q1[f'meanBETY_TEmeas_TC{lab}'] = np.average(BETY)

        analysis_data_q2[f'rmsBETX_TEmeas_TC{lab}']  = analysis_data_q1[f'rmsBETX_TEmeas_TC{lab}']  = rms(BETX)
        analysis_data_q2[f'rmsBETY_TEmeas_TC{lab}']  = analysis_data_q1[f'rmsBETY_TEmeas_TC{lab}']  = rms(BETY)

        analysis_data_q2[f'rmsBETXY_TEmeas_TC{lab}'] = analysis_data_q1[f'rmsBETXY_TEmeas_TC{lab}'] = np.sqrt(rms(BETX)**2 + rms(BETY)**2)

        analysis_data_q2[f'maxBETX_TEmeas_TC{lab}']  = analysis_data_q1[f'maxBETX_TEmeas_TC{lab}']  = max(abs(BETX))
        analysis_data_q2[f'maxBETY_TEmeas_TC{lab}']  = analysis_data_q1[f'maxBETY_TEmeas_TC{lab}']  = max(abs(BETY))

        #analysis_data_q2['maxBETXY_TEmeas_TC{lab}'] = analysis_data_q1['maxBETXY_TEmeas_TC{lab}'] = np.sqrt(max(abs(BETX))**2 + max(abs(BETY))**2)
        analysis_data_q2[f'DQX_TEmeas_TC{lab}']  = analysis_data_q1[f'DQX_TEmeas_TC{lab}']  = DQX
        analysis_data_q2[f'DQY_TEmeas_TC{lab}']  = analysis_data_q1[f'DQY_TEmeas_TC{lab}']  = DQY
    
    
    for lab,withcal in zip([0,1,2],[False,True,True]):
        (q1_BETX, q1_BETY) = q1_errors.get_generated_betabeating_real(with_calibration=withcal, type_pair_calibration=lab)
        (q2_BETX, q2_BETY) = q2_errors.get_generated_betabeating_real(with_calibration=withcal, type_pair_calibration=lab)
        BETX = q1_BETX + q2_BETX
        BETY = q1_BETY + q2_BETY

        (q1_DQX, q1_DQY) = q1_errors.get_generated_tuneshift_real(with_calibration=withcal, type_pair_calibration=lab)
        (q2_DQX, q2_DQY) = q2_errors.get_generated_tuneshift_real(with_calibration=withcal, type_pair_calibration=lab)
        DQX = q1_DQX + q2_DQX
        DQY = q1_DQY + q2_DQY

        analysis_data_q2[f'meanBETX_TEreal_TC{lab}'] = analysis_data_q1[f'meanBETX_TEreal_TC{lab}'] = np.average(BETX)
        analysis_data_q2[f'meanBETY_TEreal_TC{lab}'] = analysis_data_q1[f'meanBETY_TEreal_TC{lab}'] = np.average(BETY)

        analysis_data_q2[f'rmsBETX_TEreal_TC{lab}']  = analysis_data_q1[f'rmsBETX_TEreal_TC{lab}']  = rms(BETX)
        analysis_data_q2[f'rmsBETY_TEreal_TC{lab}']  = analysis_data_q1[f'rmsBETY_TEreal_TC{lab}']  = rms(BETY)

        analysis_data_q2[f'rmsBETXY_TEreal_TC{lab}'] = analysis_data_q1[f'rmsBETXY_TEreal_TC{lab}'] = np.sqrt(rms(BETX)**2 + rms(BETY)**2)

        analysis_data_q2[f'maxBETX_TEreal_TC{lab}']  = analysis_data_q1[f'maxBETX_TEreal_TC{lab}']  = max(abs(BETX))
        analysis_data_q2[f'maxBETY_TEreal_TC{lab}']  = analysis_data_q1[f'maxBETY_TEreal_TC{lab}']  = max(abs(BETY))

        #analysis_data_q2['maxBETXY_TEreal_TC{lab}'] = analysis_data_q1['maxBETXY_TEreal_TC{lab}'] = np.sqrt(max(abs(BETX))**2 + max(abs(BETY))**2)
        analysis_data_q2[f'DQX_TEreal_TC{lab}']  = analysis_data_q1[f'DQX_TEreal_TC{lab}']  = DQX
        analysis_data_q2[f'DQY_TEreal_TC{lab}']  = analysis_data_q1[f'DQY_TEreal_TC{lab}']  = DQY
    
    
    return analysis_data_q1,analysis_data_q2






main()
