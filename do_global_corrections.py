from datetime import datetime
import os
from os import system
from subprocess import Popen
from omc3.scripts.fake_measurement_from_model import generate as fake_measurement
from omc3 import global_correction
from omc3.response_creator import create_response_entrypoint
from pairs import CorrectabilityError, Pairs
from q1_errors import Q1Pairs
from q2_errors import Q2Pairs
from magnet_errors import *
from pathlib import Path
import tfs
import numpy as np
import pandas as pd
from typing import List, Union, Tuple

# path to madx executable
MADX = "/home/awegsche/programs/madx/madx-gnu64"

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
MAX_SIMS = 1000

# some more constants
SUMM_PAIRING = "summ_pairing.tfs"
SUMM_SUM = "summ_sum.tfs"
SUMM_DIFF = "summ_diff.tfs"
SUMM_BBEAT = "summ_bbeat.tfs"


def main():
    """
    """

    summ = Summary()

    for i in range(MAX_SIMS):
        do_analysis(summ)
        print(summ)

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

class Summary():
    def __init__(self) -> None:
        self.not_correctable = 0
        self.low_bbeat = 0
        self.general_failed = 0
        self.passed = 0

    def __repr__(self) -> str:
        total = self.not_correctable + self.general_failed + self.passed + self.low_bbeat
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



# ---- ANALYSIS ------------------------------------------------------------------------------------
def do_analysis(summ: Summary):
    """ runs a full analysis, e.g. 
    - runs madx and measures bbeat before and after corrs with initial distribution
    - does the sorting, according to certain criteria
    - runs madx and corrs again after sorting
    """

    q2_errors = Q2Pairs(10,2)
    q1_errors = Q1Pairs(10,2)

    # initial distribution, the try ... except block is to remove cases that would either not result
    # in a stable machinge anyways or would not be corrected by us because the bbeat is already fine
    try:
        check1, err1, diff1 = do_sim(q1_errors, q2_errors)

    except CorrectabilityError as e:
        summ.not_correctable = summ.not_correctable + 1 
        print("test case not correctable")
        print(e)
        return
    except Exception as e:
        summ.general_failed = summ.general_failed + 1
        print(e)
        return

    if err1 < 0.05:
        summ.low_bbeat = summ.low_bbeat + 1
        print("Low initial beta beating, don't try to improve")
        return

    summ.passed = summ.passed + 1

    def sort_and_sim(summ_filename):
        """ applies the sorting, then simulates the sorted lattice, measures beta beating before
        and after corrections
        """

        # ---- do simulations, write results to df ------------------------------------------------
        check2, err2, diff2 = do_sim(q1_errors, q2_errors)

        params = {}
        q2_errors.log_strengths(params)
        q1_errors.log_strengths(params)
        params["BBEAT"] = err1
        params["BBEAT_AFTER"] = err2
        params["CORR"] = diff1
        params["CORR_AFTER"] = diff2
        params["CHECK"] = check1
        params["CHECK_AFTER"] = check2

        write_summary(params, summ_filename)

        # ---- end of sort_and_sim ----------------------------------------------------------------

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
        q1_errors.sort_sum()
        q2_errors.sort_sum()
        sort_and_sim(SUMM_SUM)

    except:
        return


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
    return rms(err_twiss["BETX"]/model_twiss["BETX"]-1)


def do_sim(q1_errors: Q1Pairs, q2_errors: Q2Pairs) -> Tuple[float, float, float]:
    """
    runs madx and simulates corrections
    returns:
        (err, check, diff)
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
    err_twiss = tfs.read(MODEL_TWISS)
    model_twiss = tfs.read("model1/twiss.dat")

    #plt.plot(check_twiss["S"], check_twiss["BETX"]/model_twiss["BETX"]-1, label="check")
    #plt.plot(err_twiss["S"], err_twiss["BETX"]/model_twiss["BETX"]-1, label="err")
    #plt.legend()

    #plt.show()

    # reconstructed beta beating
    rms_check = rms(check_twiss["BETX"]/model_twiss["BETX"]-1)

    # error beta beating
    rms_err = rms(err_twiss["BETX"]/model_twiss["BETX"]-1)

    # difference (how close do we get?)
    rms_diff = rms(check_twiss["BETX"]/err_twiss["BETX"]-1)

    print(f"rms check: {rms_check}, rms err: {rms_err}")
    print(f"rms diff: {rms_diff}")

    return rms_check, rms_err, rms_diff


def rms(array):
    """ root mean square """
    return np.sqrt(np.mean(array**2))

main()
