from magnet_errors import *
from typing import Tuple

class Pairs:
    def __init__(self):
        self.cold_masses = []
        self.selected_permutation = 0
        self.permutations = []

    
    def positions(self):
        """ returns a list of permutated indices """
        raise NotImplementedError()

    def get_magnet_count(self):
        """ returns the number of magnets """
        raise NotImplementedError()

    def get_magnet(self, index):
        """ returns a tuple (magnet_name, magnet_error) """
        raise NotImplementedError()

    def get_pair_count(self):
        """ returns the number of pairs """
        raise NotImplementedError()

    def get_pair(self, index: int) -> Tuple["MagnetError", "MagnetError"]:
        """ returns a natural pair (magnet1, magnet2)
        where 'natural' means, A/B parts of Q2, whereas the whole Q1 is paired with Q3:
            Q2a and Q2b
            Q1ab Q3ab
        """
        raise NotImplementedError()

    def write_errors_to_file(self, filename):
        """ writes the errors to a file """
        with open(filename, 'w') as error_file:
            for i in range(self.get_magnet_count()):
                (name, error) = self.get_magnet(i)
                MagnetPair.write_error_to_file(error_file, name, error.real_error)


    def log_strengths(self, parameters):
        for i in range(self.get_magnet_count()):
            (name, error) = self.get_magnet(i)
            parameters[name] = error.real_error
