from datetime import datetime
import os
from os import system
from subprocess import Popen
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

# path to madx executable
#MADX = "/home/awegsche/programs/madx/madx-gnu64"
MADX = "/home/thpugnat/Documents/CERN/madx"

# path to twiss output
MODEL_TWISS = "model/twiss_err_b1.tfs"

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


# turn off certain parts of the sim, leave all at `TRUE` for more than one simulation or you will
# get the same result over and over
DO_MADX = True
DO_FR = True
DO_CORR = True

# error specification, this has to be kept up to date with estimations from Massimo and Ezio
AMP_REAL_ERROR = 50
AMP_MEAS_ERROR = 2

# maximum number of simulations, `Ctrl-C` stops the program early
MAX_SIMS = 100 #1000

NUMB_PERMUT_CROSSSCORING = 100

# some more constants
SUMM_PAIRING = "summ_pairing.tfs"
SUMM_SUM = "summ_sum.tfs"
SUMM_SUMBETA = "summ_sumbeta.tfs"
SUMM_SUMBBEATING = "summ_sumbbeating.tfs"
SUMM_SUMBBEATING_V2 = "summ_sumbbeating_v2.tfs"

SUMM_DIFF = "summ_diff.tfs"
SUMM_BBEAT = "summ_bbeat.tfs"

# some flags
FLAG_DEBUG=False

def main():
    """
    """

    summ = Summary()

    for i in range(MAX_SIMS):
        do_analysis(summ)
        print(summ)

        if FLAG_DEBUG:
            with open(f"sim_summary{datetime.now()}.txt", "w") as summfile:
                summfile.write(str(summ))
                summfile.write("\n")



# ---- summary (helpfull to check the global state of the simulations) -----------------------------
def write_summary(parameters, filename):
    print("writing summary tfs")
    tfs_summary = tfs.read(filename) if os.path.exists(filename) else pd.DataFrame()
    index = len(tfs_summary.index)

    tfs.write(filename,
              pd.concat([tfs_summary, pd.DataFrame(parameters, index=[index])])
              )



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
        tfs_permutation = pd.DataFrame(0,index=range(len(q1_errors.permutations)+len(q2_errors.permutations)),columns=keys)

        q2_BBETX, q2_BBETY = q2_errors.get_generated_betabeating()
        for i in range(len(q1_errors.permutations)):
            q1_errors.selected_permutation = i

            q1_BBETX, q1_BBETY = q1_errors.get_generated_betabeating()

            params = {}
            q1_errors.log_strengths(params)
            q2_errors.log_strengths(params)
            params["DBETX_MIN"] = np.min(q1_BBETX+q2_BBETX)
            params["DBETY_MIN"] = np.min(q1_BBETY+q2_BBETY)
            params["DBETX_RMS"] = rms(q1_BBETX+q2_BBETX)
            params["DBETY_RMS"] = rms(q1_BBETY+q2_BBETY)
            params["DBETX_MAX"] = np.max(q1_BBETX+q2_BBETX)
            params["DBETY_MAX"] = np.max(q1_BBETY+q2_BBETY)
            params["PERM_Q1"]   = ",".join([str(i) for i in q1_errors.permutations[q1_errors.selected_permutation]])
            params["PERM_Q2"]   = ",".join([str(i) for i in q2_errors.permutations[q2_errors.selected_permutation]])

            for kk in params.keys():
                tfs_permutation.loc[i,kk]=params[kk]

        q1_errors.selected_permutation = q1_best_comb
        q2_errors.selected_permutation = q2_best_comb

        q1_BBETX, q1_BBETY = q1_errors.get_generated_betabeating()
        for i in range(len(q2_errors.permutations)):
            q2_errors.selected_permutation = i

            q2_BBETX, q2_BBETY = q2_errors.get_generated_betabeating()

            params = {}
            q1_errors.log_strengths(params)
            q2_errors.log_strengths(params)
            params["DBETX_MIN"] = np.min(q1_BBETX+q2_BBETX)
            params["DBETY_MIN"] = np.min(q1_BBETY+q2_BBETY)
            params["DBETX_RMS"] = rms(q1_BBETX+q2_BBETX)
            params["DBETY_RMS"] = rms(q1_BBETY+q2_BBETY)
            params["DBETX_MAX"] = np.max(q1_BBETX+q2_BBETX)
            params["DBETY_MAX"] = np.max(q1_BBETY+q2_BBETY)
            params["PERM_Q1"]   = ",".join([str(i) for i in q1_errors.permutations[q1_errors.selected_permutation]])
            params["PERM_Q2"]   = ",".join([str(i) for i in q2_errors.permutations[q2_errors.selected_permutation]])

            for kk in params.keys():
                tfs_permutation.loc[i+len(q1_errors.permutations),kk]=params[kk]


        q1_errors.selected_permutation = q1_best_comb
        q2_errors.selected_permutation = q2_best_comb

        #if os.path.exists("permutation.tfs"):
        #    tfs_file = tfs.read("permutation.tfs")
        #    seed = max(tfs_file.SEED) + 1
        #else:
        #    tfs_file = pd.DataFrame()
        #    seed = 0
                
        tfs_permutation['SEED'] = summ.total()
        tfs.write(f"permutation_s{summ.total():d}.tfs", tfs_permutation)



# ---- summary (helpfull to check the global state of the simulations) -----------------------------
def write_betabeating_table(q1_errors, q2_errors,optic):
    if not os.path.exists("comparison_betabeating.tfs"):
        # compare
        err_twiss = tfs.read(MODEL_TWISS)
        model_twiss = tfs.read("model1/twiss.dat")
    
        mask_err_monitor   = mask_monitor(  err_twiss)
        mask_model_monitor = mask_monitor(model_twiss)

        q1_BBETX, q1_BBETY = q1_errors.get_generated_betabeating()
        q2_BBETX, q2_BBETY = q2_errors.get_generated_betabeating()
    
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
        
        
        print(f'\n\nmodel1/twiss_elements.dat   Qx, Qy = {optic.headers["Q1"]}, {optic.headers["Q2"]}')
        print(f'model1/twiss.dat            Qx, Qy = {optic.headers["Q1"]}, {optic.headers["Q2"]}')
        print(f'{MODEL_TWISS:<25}   Qx, Qy = {optic.headers["Q1"]}, {optic.headers["Q2"]}\n\n')
    
        tfs.write("comparison_betabeating.tfs", pd.DataFrame(params))
        


class Summary():
    def __init__(self) -> None:
        self.not_correctable = 0
        self.low_bbeat = 0
        self.general_failed = 0
        self.passed = 0

    def __repr__(self) -> str:
        total = self.total()
        return f"""
Simulation Summary
==================
not correctable : {self.not_correctable:5} ({self.not_correctable/total * 100.0:4.0f}%)
low beta beat   : {self.low_bbeat:5} ({self.low_bbeat/total * 100.0:4.0f}%)
general failed  : {self.general_failed:5} ({self.general_failed/total * 100.0:4.0f}%)
passed          : {self.passed:5} ({self.passed/total * 100.0:4.0f}%)
--------------------------------
total           : {total:5}
"""
    def __str__(self) -> str:
        return self.__repr__()
        
    def total(self) -> int:
    	return self.not_correctable + self.general_failed + self.passed + self.low_bbeat



# ---- ANALYSIS ------------------------------------------------------------------------------------
def do_analysis(summ: Summary):
    """ runs a full analysis, e.g. 
    - runs madx and measures bbeat before and after corrs with initial distribution
    - does the sorting, according to certain criteria
    - runs madx and corrs again after sorting
    """
    
    try:
        optic_without_error = tfs.read('model1/twiss_elements.dat')
        #print(optic_without_error)
        #print(optic_without_error.columns)
    except:
        print('WARNING: No twiss file found for the sorting of the magnet. The betas will be set to 1.')
        optic_without_error = None

    q2_errors = Q2Pairs(10,2, optic = optic_without_error)
    q1_errors = Q1Pairs(10,2, optic = optic_without_error)

    # initial distribution, the try ... except block is to remove cases that would either not result
    # in a stable machinge anyways or would not be corrected by us because the bbeat is already fine
    try:
        check1x, check1y, err1x, err1y, diff1x, diff1y = do_sim(q1_errors, q2_errors)

    except CorrectabilityError as e:
        summ.not_correctable = summ.not_correctable + 1 
        print("test case not correctable")
        print(e)
        return
    except Exception as e:
        summ.general_failed = summ.general_failed + 1
        print(e)
        return

    #if err1x < 0.05:
    #    summ.low_bbeat = summ.low_bbeat + 1
    #    print("Low initial beta beating, don't try to improve")
    #    return

    summ.passed = summ.passed + 1

    def sort_and_sim(summ_filename):
        """ applies the sorting, then simulates the sorted lattice, measures beta beating before
        and after corrections
        """

        # ---- do simulations, write results to df ------------------------------------------------
        check2x, check2y, err2x, err2y, diff2x, diff2y = do_sim(q1_errors, q2_errors)

        params = {}
        q2_errors.log_strengths(params)
        q1_errors.log_strengths(params)
        params["BBEATX"] = err1x
        params["BBEATY"] = err1y
        params["BBEATX_AFTER"] = err2x
        params["BBEATY_AFTER"] = err2y
        params["CORRX"] = diff1x
        params["CORRY"] = diff1y
        params["CORRX_AFTER"] = diff2x
        params["CORRY_AFTER"] = diff2y
        params["CHECKX"] = check1x
        params["CHECKY"] = check1y
        params["CHECKX_AFTER"] = check2x
        params["CHECKY_AFTER"] = check2y
        # add permutation
        params["PERM_Q1"] = ",".join([str(i) for i in q1_errors.permutations[q1_errors.selected_permutation]])
        params["PERM_Q2"] = ",".join([str(i) for i in q2_errors.permutations[q2_errors.selected_permutation]])

        write_summary(params, summ_filename)

        # ---- end of sort_and_sim ----------------------------------------------------------------

    #q1_errors.sort_betabeating()
    #q2_errors.sort_betabeating()
    #sort_and_sim(SUMM_SUMBBEATING)
    try:
        # ---- sorting based on initial beta beating ---------------------------------------------

        def sort_on_bbeat_q1(pairs: Q1Pairs):
            return run_madx_for_sorting(pairs, q2_errors)

        def sort_on_bbeat_q2(pairs: Q2Pairs):
            return run_madx_for_sorting(q1_errors, pairs)

        # q1_errors.sort(sort_on_bbeat_q1)
        # q2_errors.sort(sort_on_bbeat_q2)
        # sort_and_sim(SUMM_BBEAT)

        # # ---- naive sorting methods (but fast, we can crunch through 80k combinations in < 1s ----
        # # ---- this sorts on the difference of errors --------------------------------------------
        # q1_errors.sort_diff()
        # q2_errors.sort_diff()
        # sort_and_sim(SUMM_DIFF)

        # ---- this sorts on the sum of errors (because nearby magnets should self cancel, right? -
        #q1_errors.sort_sum()
        #q2_errors.sort_sum()
        #sort_and_sim(SUMM_SUM)

        # ---- this sorts on the sum of errors weigthed by the betas ------------------------------
        #q1_errors.sort_beta()
        #q2_errors.sort_beta()
        #sort_and_sim(SUMM_SUMBETA)

        # ---- this sorts on the sum of errors (because nearby magnets should self cancel, right? -
        list_bestscore_q1 = q1_errors.sort_betabeating()
        list_bestscore_q2 = q2_errors.sort_betabeating()
        sort_and_sim(SUMM_SUMBBEATING)

        # ---- this sorts on the sum of errors (because nearby magnets should self cancel, right? -
        q1_errors, q2_errors = sort_betabeating_v2(q1_errors, list_bestscore_q1[:NUMB_PERMUT_CROSSSCORING], 
                                                   q2_errors, list_bestscore_q2[:NUMB_PERMUT_CROSSSCORING])
        sort_and_sim(SUMM_SUMBBEATING_V2)

    #except:
    #    return
    except Exception as e:
        print('\n# '+'='*60)
        print('#   * BUG IN SORTING')
        print(e)
        print('# '+'='*60+'\n')
        return
    
    #write_betabeating_table(q1_errors, q2_errors,optic_without_error)
    write_permutation_table(q1_errors, q2_errors, summ)



def run_madx(q1_errors: Q1Pairs, q2_errors: Q2Pairs):
    print("running madx")

    q2_errors.write_errors_to_file("model/errors_Q2.madx")
    q1_errors.write_errors_to_file("model/errors_Q1.madx")

    if not q1_errors.check_correctability():
        raise CorrectabilityError([x.real_error for x in q1_errors.cold_masses])
    if not q2_errors.check_correctability():
        raise CorrectabilityError([x.real_error for x in q2_errors.cold_masses])

    template_file = open("job.inj.madx", "r")
    jobfile = open("model/run_job.inj.madx", "w")
    jobfile.write(template_file.read().replace("_TRACK_", "0")) # no tracking
    template_file.close()
    jobfile.close()

    # remove twiss output, its existance after running madx indicates success
    system(f"rm {MODEL_TWISS}")

    with open(os.devnull, "w") as devnull:
        Popen([MADX, "run_job.inj.madx"], stdout=devnull, cwd="model").wait()

    if not os.path.exists(f"{MODEL_TWISS}"):
        raise RuntimeError("twiss_err_b1.tfs does not exist")


def run_madx_for_sorting(q1_errors: Q1Pairs, q2_errors: Q2Pairs) -> float:

    run_madx(q1_errors, q2_errors)

    # error beta beating
    err_twiss = tfs.read(MODEL_TWISS)
    model_twiss = tfs.read("model1/twiss.dat")
    
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

    accel_params = dict(
        accel="lhc",
        year="hl16", # to be checked
        beam=1,
        model_dir=Path("model1").absolute(),
            )


    if DO_FR:
        print("creating response entrypoint")
        create_response_entrypoint(
            outfile_path=Path("model1/FullResponse.h5"),
            creator="madx",
            optics_params=['PHASEX', 'PHASEY', 'BETX', 'BETY', 'Q'],
            variable_categories=VAR_CATS,
            delta_k=1.0e-7,
            **accel_params 
        )


    # fake measurement
    fake_measurement(twiss=MODEL_TWISS,
                     model="model1/twiss.dat",
                     outputdir="fake_measurements")


    if DO_CORR:
        print("running global correction")
        global_correction.global_correction_entrypoint(
            meas_dir="fake_measurements",
            output_dir="global_corrections",
            fullresponse_path="model1/FullResponse.h5",
            iterations=1,
            optics_params=['PHASEX', 'PHASEY', 'BETX', 'BETY', 'Q'],
            variable_categories=VAR_CATS,
            **accel_params,
                )
        system("cp scripts_global_corrs/rerun.job.madx global_corrections/rerun.job.madx")
        with open(os.devnull, "w") as devnull:
            Popen([MADX, "rerun.job.madx"], cwd="global_corrections", stdout=devnull).wait()




    # compare
    check_twiss = tfs.read("global_corrections/twiss_global_b1.tfs")
    err_twiss   = tfs.read(MODEL_TWISS)
    model_twiss = tfs.read("model1/twiss.dat")
    
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

    score = q1_errors.score_rms_betabeating(q1_errors, q2_errors)
    print(f" -- initial score: {score}")

    next_progress = 0.1
    best_comb = ( q1_errors.selected_permutation , q2_errors.selected_permutation )
        
    for i1 in range(len_lq1):
        for i2 in range(len_lq2):
            if ( i1*len_lq2 + i2 ) / ( len_lq1*len_lq2 ) > next_progress:
                print(f"progress: {100* ( i1*len_lq2 + i2 ) / ( len_lq1*len_lq2 ):.0f} %")
                next_progress += 0.1

            q1_errors.selected_permutation = list_q1[i1]
            q2_errors.selected_permutation = list_q2[i2]

            sum = q1_errors.score_rms_betabeating(q1_errors, q2_errors)

            if sum < score:
                score = sum
                best_comb = ( q1_errors.selected_permutation , q2_errors.selected_permutation )
                
    q1_errors.selected_permutation = best_comb[0]
    q2_errors.selected_permutation = best_comb[1]
                
    print(f"final score: {score}")
    return q1_errors, q2_errors


def rms(array):
    """ root mean square """
    return np.sqrt(np.mean(array**2))

main()
