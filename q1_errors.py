from magnet_errors import *
import itertools
from pairs import Pairs
from typing import Tuple
import numpy as np

MAGNETS = ["MQXA.1L1", "MQXA.1R1", "MQXA.1L5", "MQXA.1R5"]
# HL-LHC
MAGNETS = [
        "MQXFA.A1L1", "MQXFA.B1L1", "MQXFA.A1R1", "MQXFA.B1R1",
        "MQXFA.A1L5", "MQXFA.B1L5", "MQXFA.A1R5", "MQXFA.B1R5",
        "MQXFA.A3L1", "MQXFA.B3L1", "MQXFA.A3R1", "MQXFA.B3R1",
        "MQXFA.A3L5", "MQXFA.B3L5", "MQXFA.A3R5", "MQXFA.B3R5",
        ]

AMP_REAL_ERROR = 50  # Full range of integrated gradient error in units
AMP_MEAS_ERROR = 2  # Random error of measurement

class Q1Pairs(Pairs):
    """ Actually Q1 and Q3, because those are "pairable" """

    NAME = "Q1"

    def __init__(self,
                 real_error: float = AMP_REAL_ERROR,
                 meas_error: float = AMP_MEAS_ERROR,
                 stage: int = 1):
        self.cold_masses = [MagnetError(real_error , meas_error) for _ in range(self.get_magnet_count())]
        self.stage = stage
        self.selected_permutation = 0
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

    def get_pair_count(self):
        return len(MAGNETS) // 2

    def get_pair(self, index: int) -> Tuple[MagnetError, MagnetError]:
        assert index < self.get_pair_count()
        return (self.get_magnet(index)[1], self.get_magnet(index+self.get_pair_count())[1])

    def get_pair_names(self, index: int) -> Tuple[str, str]:
        assert index < self.get_pair_count()
        return (self.get_magnet(index)[0], self.get_magnet(index+self.get_pair_count())[0])

    def __getitem__(self, index: int):
        # permuted bucket
        bucket = index//2
        # even / odd
        position = index % 2

        return self.cold_masses[self.positions()[bucket]*2 + position]

