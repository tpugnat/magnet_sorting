from magnet_errors import *
from typing import Tuple
import numpy as np

class Pairs:
    NAME = None

    # ---- Init and accessors ----------------------------------------------------------------------
    def __init__(self):
        self.cold_masses = []
        self.selected_permutation = 0
        self.permutations = []
        self.stage = 1

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

    # ---- Output to files -------------------------------------------------------------------------
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

    # ---- Sorting ---------------------------------------------------------------------------------
    @staticmethod
    def score_diff(errors: "Pairs"):
        """
        Sorts according to the difference between pair elements
        (this favors errors that are close to each other)
        """

        return np.sum([
            (errors.get_pair(i)[0].real_error - errors.get_pair(i)[1].real_error)**2
            for i in range(errors.get_pair_count())])

    @staticmethod
    def score_sum(errors):
        """
        Sorts according to sum of pair elements
        (this favors errors that naturally cancel each other)
        """

        return np.sum([
            (errors.get_pair(i)[0].real_error + errors.get_pair(i)[1].real_error)**2
            for i in range(errors.get_pair_count())])

    def sort(self, sort_fn):
        """
        Performs the sorting.

        sort_fn should take a Pairs object and return a score
        """

        print(f"sorting {self.NAME}") 
        print(f"searching for best combination in {len(self.permutations)} permutations")
        print(f"using the diff method")

        score = sort_fn(self)
        print(f" -- initial score: {score}")

        next_progress = 0.1
        best_comb = 0
        
        for i in range(len(self.permutations)):
            if i / len(self.permutations) > next_progress:
                print(f"progress: {100* i / len(self.permutations):.0f} %")
                next_progress += 0.1

            self.selected_permutation = i

            sum = sort_fn(self)

            if sum < score:
                score = sum
                best_comb = i
                
        self.selected_permutation = best_comb
        print(f"final score: {score}")

    def sort_diff(self):
        """ Convenience method for performing the sorting according to difference"""
        self.sort(self.score_diff)

    def sort_sum(self):
        """ Convenience method for performing the sorting according to sum"""
        self.sort(self.score_sum)

