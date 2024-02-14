from pandas.core.frame import itertools
from magnet_errors import *
from pairs import Pairs, mask_monitor
from typing import Tuple
import numpy as np
import pandas as pd

#MAGNETS = ["MQXB.A2L1", "MQXB.B2L1", "MQXB.A2R1", "MQXB.B2R1",
#"MQXB.A2L5", "MQXB.B2L5", "MQXB.A2R5", "MQXB.B2R5"]
#HL-LHC

STR_FLE = "strengths.tfs"
MAGNETS = ["MQXFB.A2L1", "MQXFB.B2L1", "MQXFB.A2R1", "MQXFB.B2R1",
           "MQXFB.A2L5", "MQXFB.B2L5", "MQXFB.A2R5", "MQXFB.B2R5"]
#MAGNETS = []

AMP_REAL_ERROR = 10  # Full range of integrated gradient error in units
AMP_MEAS_ERROR = 2  # Random error of measurement


class Q2Pairs(Pairs):
    NAME = "Q2"
    def __init__(self, real_error, meas_error, optic = None):
        self.cold_masses = [MagnetError(real_error , meas_error) for _ in range(8)]
        self.selected_permutation = 0
        self.permutations = list(itertools.permutations(range(8)))
        if optic is None:
            #self.BETX = np.ones(len(MAGNETS))
            #self.BETY = np.ones(len(MAGNETS))
            
            self.monitor_responce_betx = None
            self.monitor_responce_bety = None
        else:
            #idx = [optic.index[optic["NAME"] == nn][0] for nn in MAGNETS]
            #print(" - Q1-Q3:")
            #print(idx)
            #print(optic["NAME"][idx])
            #self.BETX = optic["BETX"][idx].values
            #self.BETY = optic["BETY"][idx].values
            
            Qx=optic.headers['Q1']
            Qy=optic.headers['Q2']
            msk_monitor = mask_monitor(optic)
            msk_magnets = {nn:optic["NAME"] == nn for nn in MAGNETS}
            
            self.monitor_responce_bbetx = pd.DataFrame(0,index=optic[msk_monitor].index,columns=MAGNETS)
            self.monitor_responce_bbety = pd.DataFrame(0,index=optic[msk_monitor].index,columns=MAGNETS)
            
            for nn in MAGNETS:
                for nnn in [nop for nop in optic.NAME if f"{nn}.." in nop]:
                    idx = optic.index[optic["NAME"] == nnn][0]
                    self.monitor_responce_bbetx.loc[:,nn] += -optic.K1L[idx] * optic.BETX[idx] * np.cos(2*np.pi*( 2*abs(optic.MUX[idx]-optic.MUX[msk_monitor]) - Qx ))
                    self.monitor_responce_bbety.loc[:,nn] +=  optic.K1L[idx] * optic.BETY[idx] * np.cos(2*np.pi*( 2*abs(optic.MUY[idx]-optic.MUY[msk_monitor]) - Qy ))
            self.monitor_responce_bbetx = self.monitor_responce_bbetx * 1.0e-4 / (2*np.sin(2*np.pi*Qx))
            self.monitor_responce_bbety = self.monitor_responce_bbety * 1.0e-4 / (2*np.sin(2*np.pi*Qy))
        #self.sign_phase = [-1 if mm[-2:] in ['L1','R5'] else +1 for mm in MAGNETS]

    def positions(self):
        return self.permutations[self.selected_permutation]

    def get_magnet_count(self):
        return len(MAGNETS)

    def get_magnet(self, index):
        """ returns a tuple (magnet_name, magnet_error) """
        return (MAGNETS[index], self[index])

    #def get_magnet_twiss(self, index: int):
    #    """ returns a tuple (magnet_name, betax, betay, sign) """
    #    return (MAGNETS[index], self.BETX[index], self.BETY[index], self.sign_phase[index])

    def get_generated_betabeating(self):
        """ returns a tuple (Dbetax, Dbetay) """
        bbetx=0 ; bbety=0;
        if self.monitor_responce_bbetx is not None:
            for ii in range(self.get_magnet_count()):
                bbetx+= self.get_magnet(ii)[1].real_error*self.monitor_responce_bbetx[self.get_magnet(ii)[0]]
                bbety+= self.get_magnet(ii)[1].real_error*self.monitor_responce_bbety[self.get_magnet(ii)[0]]
        return (bbetx, bbety)

    def get_pair_count(self):
        return len(MAGNETS) // 2

    def get_pair(self, index: int) -> Tuple[MagnetError, MagnetError]:
        return (self.get_magnet(2*index)[1], self.get_magnet(2*index+1)[1])

    #def get_pair_betx(self, index: int) -> Tuple[float, float]:
    #    return (self.get_magnet_twiss(2*index)[1], self.get_magnet_twiss(2*index+1)[1])

    #def get_pair_bety(self, index: int) -> Tuple[float, float]:
    #    return (self.get_magnet_twiss(2*index)[2], self.get_magnet_twiss(2*index+1)[2])

    #def get_pair_sign(self, index: int) -> Tuple[float, float]:
    #    return (self.get_magnet_twiss(2*index)[3], self.get_magnet_twiss(2*index+1)[3])

    def get_pair_names(self, index: int) -> Tuple[str, str]:
        return (self.get_magnet(2*index)[0], self.get_magnet(2*index+1)[0])

    def __getitem__(self, index):
        return self.cold_masses[self.positions()[index]]

