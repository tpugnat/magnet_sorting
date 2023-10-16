from os import system
import os
import pandas as pd
from omc3.hole_in_one import hole_in_one_entrypoint
from subprocess import Popen
import sortq2withuncertainty
import shutil
import pathlib
import traceback
import re
import argparse

import tfs
from matplotlib import pyplot as plt
import numpy as np
import plot_bb

MADX = "/home/awegsche/programs/madx/madx-gnu64"

RERUN_MADX = True
DO_ANALYSIS = False
DO_TRACK = "1" if DO_ANALYSIS else "0"

SUMMARY_FILENAME = "meas_bbs_ramp.txt"
SUMMARY_TFS = "summary.tfs"
SUMMARY_FILENAME_MOD = "modl_bbs_ramp.txt"


TRACKONE = "witherrorsone"

BPMR_RE = re.compile("(bpm\\S*)")

def main():
    args = parse()

    if args.reset:
        print("DELETING summary files")
        system(f"rm {SUMMARY_FILENAME}")
        system(f"rm {SUMMARY_FILENAME_MOD}")
    for i in range(args.num):
        run_sim()


def run_sim():
    model = "model1"
    files = [TRACKONE]
    outputdir = "witherrors"

    if RERUN_MADX:
        try:
            system(f"rm {TRACKONE}")
        except:
            pass

        parameters = sortq2withuncertainty.main()

        template_file = open("job.inj.madx", "r")
        jobfile = open("run_job.inj.madx", "w")
        jobfile.write(template_file.read().replace("_TRACK_", DO_TRACK))
        template_file.close()
        jobfile.close()

        system(f"rm twiss_err_b1.tfs")

        Popen([MADX, "run_job.inj.madx"]).wait()

        # get beta_beating from twiss
        if os.path.exists("twiss_err_b1.tfs"):
            twiss = tfs.read("twiss_err_b1.tfs")
            twiss_model = tfs.read("model1/twiss.dat")

            bbx = _rms(twiss_model["BETX"]/twiss["BETX"]-1)
            bby = _rms(twiss_model["BETY"]/twiss["BETY"]-1)

            parameters["BETABEAT_X"] = bbx
            parameters["BETABEAT_Y"] = bby
            parameters["NOTES"] = "b2 errors only"
        else:

            parameters["BETABEAT_X"] = float("inf")
            parameters["BETABEAT_Y"] = float("inf")
            parameters["NOTES"] = "twiss failed, b2 errors only"

        parameters["OPTICS"] = "inj"

        if DO_ANALYSIS:
            print("done running madx, converting bpm names to uppercase")
            trackone = open(TRACKONE, "r")
            lines = [BPMR_RE.sub(lambda match: match.group(1).upper(), l) for l in trackone]
            trackone.close()

            for l in lines[55:65]:
                print(l)

            system(f"rm {TRACKONE}")
            trackone = open(TRACKONE, "w")
            trackone.writelines(lines)
            trackone.close()
            print("done converting names")

m       print("writing summary tfs")
        tfs_summary = tfs.read(SUMMARY_TFS) if os.path.exists(SUMMARY_TFS) else pd.DataFrame()
        index = len(tfs_summary.index)

        tfs.write(SUMMARY_TFS,
                  pd.concat([tfs_summary, pd.DataFrame(parameters, index=[index])])
                  )


    if not pathlib.Path("witherrorsone").exists():
        with open(SUMMARY_FILENAME, "a") as bbs:
            bbs.write(f"failed madx\n")
        return


    try:
        if DO_ANALYSIS:
            print("run hole in one")
            hole_in_one_entrypoint(
                    # harpy
                    harpy=True,
                    files=files,
                    tbt_datatype="ptc",
                    tunes=[0.31,0.32,0.0],
                    turn_bits=16,
                    to_write=['lin', 'bpm_summary', 'spectra'],
                    clean=True,
                    # optics
                    optics=True,
                    calibrationdir="/afs/cern.ch/eng/sl/lintrack/LHC_commissioning2017/Calibration_factors_2017/Calibration_factors_2017_beam1",
                    model_dir=model,
                    accel="lhc",
                    year="hllhc1.3",
                    beam=1,
                    three_bpm_method=True,
                    compensation="none",
                    outputdir=outputdir,
                    #lobster_rescaling=True,
                    )
            calc_bb()
            plot_bb.main(SUMMARY_FILENAME_MOD)
    except Exception as e:
        print(traceback.format_exc())
        with open(SUMMARY_FILENAME, "a") as bbs:
            bbs.write(f"failed omc3\n")
        return


def _rms(arr):
    return np.sqrt(np.mean(np.square(arr)))

def calc_bb():
    try:
        betx = tfs.read_tfs("witherrors/beta_phase_x.tfs")
    except:
        return

    errbx = betx["BETX"]/betx["BETXMDL"]-1

    rms_err = _rms(errbx**2)

    print("Beta Beating")
    print(f"from omc3 = {rms_err*100:.3f}%")

    with open(SUMMARY_FILENAME, "a") as bbs:
        bbs.write(f"{rms_err}\n")


def parse() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reset", action="store_true",
                        help="deletes summary files")
    parser.add_argument("--num", "-n", dest="num", type=int, default=1,
                        help="number of simulation runs. Defaul=1")
    return parser.parse_args()

main()
