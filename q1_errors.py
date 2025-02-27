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
                 #real_error: float = AMP_REAL_ERROR,
                 meas_error: float = 0,
                 cali_error: float = 0,
                 pres_error: float = 0,
                 stage: int = 1,
                 optic = None):
        # Set magnets to sort
        self.cold_masses = [MagnetError(ampl_meas=meas_error, ampl_cali=cali_error, ampl_prec=pres_error) for _ in range(self.get_magnet_count())]
        
        # Set pair callibrations
        self.pairs_average_meas     = pd.DataFrame([], index=range(self.get_magnet_count()), 
                                                       columns=range(self.get_magnet_count()) )
        self.pairs_calibration_meas = pd.DataFrame([], index=range(self.get_magnet_count()), 
                                                       columns=range(self.get_magnet_count()) )
        
        self.pairs_average_real     = pd.DataFrame([], index=range(self.get_magnet_count()), 
                                                       columns=range(self.get_magnet_count()) )
        self.pairs_calibration_real = pd.DataFrame([], index=range(self.get_magnet_count()), 
                                                       columns=range(self.get_magnet_count()) )
        
        self.update_pair_calibration()
#         for mmA in range(self.get_magnet_count()):
#             for mmB in range(mmA,self.get_magnet_count()):
#                 self.pairs_average.loc[mmA,mmB]=5e-1*(self.cold_masses[mmA].meas_error +
#                                                       self.cold_masses[mmB].meas_error)
#                 self.pairs_average.loc[mmB,mmA]=self.pairs_average.loc[mmA,mmB] 
#                 self.pairs_calibration.loc[mmA,mmB]=1/(1+1e-4*self.pairs_average.loc[mmA,mmB])
#                 self.pairs_calibration.loc[mmB,mmA]=self.pairs_calibration.loc[mmA,mmB]
        
        # Dict for step by step saves
        self.initial_permutation    = {'id':0}
        self.best_permutation_alone = {'id':0}
        self.best_permutation_both  = {'id':0}
        self.worst_permutation_both = {'id':0}
#         self.beta_wall_MADX         = {'id':0}
        
#         self.initial_permutation    = {'id':0, 'rmsBETX':0, 'rmsBETY':0, 'meanBETX':0, 'meanBETY':0, 
#                                        'rmsBETXY':0, 'maxBETX':0, 'maxBETY':0, 'maxBETXY':0}
#         self.best_permutation_alone = {'id':0, 'rmsBETX':0, 'rmsBETY':0, 'meanBETX':0, 'meanBETY':0, 
#                                        'rmsBETXY':0, 'maxBETX':0, 'maxBETY':0, 'maxBETXY':0}
#         self.best_permutation_both  = {'id':0, 'rmsBETX':0, 'rmsBETY':0, 'meanBETX':0, 'meanBETY':0, 
#                                        'rmsBETXY':0, 'maxBETX':0, 'maxBETY':0, 'maxBETXY':0}
#         self.worst_permutation_both = {'id':0, 'rmsBETX':0, 'rmsBETY':0, 'meanBETX':0, 'meanBETY':0, 
#                                        'rmsBETXY':0, 'maxBETX':0, 'maxBETY':0, 'maxBETXY':0}
#         self.beta_wall_MADX         = {'id':0, 'rmsBETX':0, 'rmsBETY':0, 'meanBETX':0, 'meanBETY':0, 
#                                        'rmsBETXY':0, 'maxBETX':0, 'maxBETY':0, 'maxBETXY':0}
        
        # Default score function
        self.score_def_fn = None
        
        # Responce matrices from magnet error to beta-beating and tune shift
        if optic is None:
            self.monitor_responce_betx = None
            self.monitor_responce_bety = None
            self.responce_Qx = None
            self.responce_Qy = None
            
        else:
            Qx = optic.headers['Q1']
            Qy = optic.headers['Q2']
            msk_monitor = mask_monitor(optic)
            msk_magnets = {nn:optic["NAME"] == nn for nn in MAGNETS}
            
            self.monitor_responce_bbetx = pd.DataFrame(0,index=optic[msk_monitor].index,columns=MAGNETS)
            self.monitor_responce_bbety = pd.DataFrame(0,index=optic[msk_monitor].index,columns=MAGNETS)
            
            self.responce_Qx = pd.DataFrame(0,index=[0],columns=MAGNETS)
            self.responce_Qy = pd.DataFrame(0,index=[0],columns=MAGNETS)
            
            for nn in MAGNETS:
                for nnn in [nop for nop in optic.NAME if f"{nn}.." in nop]:
                    idx = optic.index[optic["NAME"] == nnn][0]
                    self.monitor_responce_bbetx.loc[:,nn] += -optic.K1L[idx] * optic.BETX[idx] * np.cos(2*np.pi*( 2*abs(optic.MUX[idx]-optic.MUX[msk_monitor]) - Qx ))
                    self.monitor_responce_bbety.loc[:,nn] +=  optic.K1L[idx] * optic.BETY[idx] * np.cos(2*np.pi*( 2*abs(optic.MUY[idx]-optic.MUY[msk_monitor]) - Qy ))
                    
                    self.responce_Qx.loc[:,nn] += -optic.K1L[idx] * optic.BETX[idx]
                    self.responce_Qy.loc[:,nn] +=  optic.K1L[idx] * optic.BETY[idx]
                    
            self.monitor_responce_bbetx = self.monitor_responce_bbetx * 1.0e-4 / (2*np.sin(2*np.pi*Qx))
            self.monitor_responce_bbety = self.monitor_responce_bbety * 1.0e-4 / (2*np.sin(2*np.pi*Qy))
            
            self.responce_Qx = self.responce_Qx * 1.0e-4 / (4*np.pi)
            self.responce_Qy = self.responce_Qy * 1.0e-4 / (4*np.pi)
            
        self.monitor_responce_bbetx_rms2 = None
        self.monitor_responce_bbety_rms2 = None
            
        # List allowed permutations
        self.stage = stage
        self.selected_permutation = 0
        if stage == 1 or stage == 2 or stage == 4:
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
        
    def set_magnets_meas_strength(self,strength_and_errors):
        for idx,var in enumerate(strength_and_errors):
            self.cold_masses[idx].set_meas_strength(*var)
        self.update_pair_calibration()
        
    def set_meas_strength_no_cold_incertitudes_for_sorting(self,strength_and_errors):
        for idx,var in enumerate(strength_and_errors):
            self.cold_masses[idx].set_meas_strength_no_cold_incertitudes_for_sorting(*var)
        self.update_pair_calibration()
        
    def set_meas_strength_no_cold_nor_unknown_incertitudes_for_sorting(self,strength_and_errors):
        for idx,var in enumerate(strength_and_errors):
            self.cold_masses[idx].set_meas_strength_no_cold_nor_unknown_incertitudes_for_sorting(*var)
        self.update_pair_calibration()
        
    def set_magnets_real_strength(self,strength_and_errors):
        for idx,var in enumerate(strength_and_errors):
            self.cold_masses[idx].set_real_strength(*var)
        self.update_pair_calibration()

    #def get_magnet_twiss(self, index: int):
    #    """ returns a tuple (magnet_name, betax, betay, sign) """
    #    return (MAGNETS[index], self.BETX[index], self.BETY[index], self.sign_phase[index])

    def get_pair_count(self):
        return len(MAGNETS) // 2

    def get_pair(self, index: int) -> Tuple[MagnetError, MagnetError]:
        assert index < self.get_pair_count()
        return (self.get_magnet(index)[1], self.get_magnet(index+self.get_pair_count())[1])

    def get_pair_names_magnet(self, index: int) -> Tuple[Tuple, Tuple]:
        assert index < self.get_pair_count()
        #return (self.get_magnet(index), self.get_magnet(index+self.get_pair_count()))
        return (self.get_magnet(2*index), self.get_magnet(2*index+1))

    def get_pair_names(self, index: int) -> Tuple[str, str]:
        assert index < self.get_pair_count()
        #return (self.get_magnet(index)[0], self.get_magnet(index+self.get_pair_count())[0])
        return (self.get_magnet(2*index)[0], self.get_magnet(2*index+1)[0])
    
    def get_pair_index(self, index: int) -> Tuple[int, int]:
        return (self.positions()[index]*2, self.positions()[index]*2+1)

    def __getitem__(self, index: int):
        # permuted bucket
        bucket = index//2
        # even / odd
        position = index % 2

        return self.cold_masses[self.positions()[bucket]*2 + position]
