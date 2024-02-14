from magnet_errors import *
import itertools
from pairs import Pairs, mask_monitor
from typing import Tuple
import numpy as np
import pandas as pd

MAGNETS = ["MQXA.1L1", "MQXA.1R1", "MQXA.1L5", "MQXA.1R5"]
# HL-LHC
MAGNETS = [
        "MQXFA.A1L1", "MQXFA.B1L1", "MQXFA.A1R1", "MQXFA.B1R1",
        "MQXFA.A1L5", "MQXFA.B1L5", "MQXFA.A1R5", "MQXFA.B1R5",
        "MQXFA.A3L1", "MQXFA.B3L1", "MQXFA.A3R1", "MQXFA.B3R1",
        "MQXFA.A3L5", "MQXFA.B3L5", "MQXFA.A3R5", "MQXFA.B3R5",
        ]
#MAGNETS = ["MQXFA.A1L5", "MQXFA.B1L5",]

AMP_REAL_ERROR = 50  # Full range of integrated gradient error in units
AMP_MEAS_ERROR = 2  # Random error of measurement

class Q1Pairs(Pairs):
    """ Actually Q1 and Q3, because those are "pairable" """

    NAME = "Q1"

    def __init__(self,
                 real_error: float = AMP_REAL_ERROR,
                 meas_error: float = AMP_MEAS_ERROR,
                 stage: int = 1,
                 optic = None):
        self.cold_masses = [MagnetError(real_error , meas_error) for _ in range(self.get_magnet_count())]
        self.stage = stage
        self.selected_permutation = 0
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
            
            Qx = optic.headers['Q1']
            Qy = optic.headers['Q2']
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
        if stage == 1 or stage == 2:
            self.permutations = list(itertools.permutations(range(self.get_magnet_count()//2)))
        elif stage == 3:
            raise NotImplementedError("slhfe")

    def positions(self):
        return self.permutations[self.selected_permutation]

    def get_magnet_count(self):
        return len(MAGNETS)

    def get_magnet(self, index: int):
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
                bbetx += self.get_magnet(ii)[1].real_error*self.monitor_responce_bbetx[self.get_magnet(ii)[0]]
                bbety += self.get_magnet(ii)[1].real_error*self.monitor_responce_bbety[self.get_magnet(ii)[0]]
        return (bbetx, bbety)

    def get_pair_count(self):
        return len(MAGNETS) // 2

    def get_pair(self, index: int) -> Tuple[MagnetError, MagnetError]:
        assert index < self.get_pair_count()
        return (self.get_magnet(index)[1], self.get_magnet(index+self.get_pair_count())[1])

    #def get_pair_betx(self, index: int) -> Tuple[float, float]:
    #    assert index < self.get_pair_count()
    #    return (self.get_magnet_twiss(index)[1], self.get_magnet_twiss(index+self.get_pair_count())[1])

    #def get_pair_bety(self, index: int) -> Tuple[float, float]:
    #    assert index < self.get_pair_count()
    #    return (self.get_magnet_twiss(index)[2], self.get_magnet_twiss(index+self.get_pair_count())[2])

    #def get_pair_sign(self, index: int) -> Tuple[float, float]:
    #    assert index < self.get_pair_count()
    #    return (self.get_magnet_twiss(index)[3], self.get_magnet_twiss(index+self.get_pair_count())[3])

    def get_pair_names(self, index: int) -> Tuple[str, str]:
        assert index < self.get_pair_count()
        return (self.get_magnet(index)[0], self.get_magnet(index+self.get_pair_count())[0])

    def __getitem__(self, index: int):
        # permuted bucket
        bucket = index//2
        # even / odd
        position = index % 2

        return self.cold_masses[self.positions()[bucket]*2 + position]
