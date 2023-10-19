from pandas.core.frame import itertools
from magnet_errors import *
from pairs import Pairs
from typing import Tuple

#MAGNETS = ["MQXB.A2L1", "MQXB.B2L1", "MQXB.A2R1", "MQXB.B2R1",
#"MQXB.A2L5", "MQXB.B2L5", "MQXB.A2R5", "MQXB.B2R5"]
#HL-LHC

STR_FLE = "strengths.tfs"
MAGNETS = ["MQXFB.A2L1", "MQXFB.B2L1", "MQXFB.A2R1", "MQXFB.B2R1",
           "MQXFB.A2L5", "MQXFB.B2L5", "MQXFB.A2R5", "MQXFB.B2R5"]

AMP_REAL_ERROR = 10  # Full range of integrated gradient error in units
AMP_MEAS_ERROR = 2  # Random error of measurement


class Q2Pairs(Pairs):
    NAME = "Q2"
    def __init__(self, real_error, meas_error):
        self.cold_masses = [MagnetError(real_error , meas_error) for _ in range(8)]
        self.selected_permutation = 0
        self.permutations = list(itertools.permutations(range(8)))

    def positions(self):
        return self.permutations[self.selected_permutation]

    def get_magnet_count(self):
        return len(MAGNETS)

    def get_pair_count(self):
        return len(MAGNETS) // 2

    def get_pair(self, index: int) -> Tuple[MagnetError, MagnetError]:
        return (self.get_magnet(2*index)[1], self.get_magnet(2*index+1)[1])

    def get_magnet(self, index):
        """ returns a tuple (magnet_name, magnet_error) """
        return (MAGNETS[index], self[index])

    def __getitem__(self, index):
        return self.cold_masses[self.positions()[index]]

