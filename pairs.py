from magnet_errors import *
from typing import Tuple
import numpy as np

# threshold on correction strengths, if the magnet error is below this value (relative for now)
# we won't (be able) to do corrections on it
CORRECTABILITY_THRESHOLD = 1.0

def mask_monitor(optic):
    return (optic.KEYWORD == 'MONITOR') & (~optic.NAME.str.contains('BPT')) & (~optic.NAME.str.contains('_ITLK')) & (~optic.NAME.str.contains('_DOROS'))


def rms(array):
    """ root mean square """
    return np.sqrt(np.mean(array**2))
    

class Pairs:
    NAME = None

    # ---- Init and accessors ----------------------------------------------------------------------
    def __init__(self):
        self.cold_masses = []
        self.selected_permutation = 0
        self.initial_permutation    = {'id':0, 'rmsBETX':0, 'rmsBETY':0, 'rmsBETXY':0, 'maxBETX':0, 'maxBETY':0, 'maxBETXY':0}
        self.best_permutation_alone = {'id':0, 'rmsBETX':0, 'rmsBETY':0, 'rmsBETXY':0, 'maxBETX':0, 'maxBETY':0, 'maxBETXY':0}
        self.best_permutation_both  = {'id':0, 'rmsBETX':0, 'rmsBETY':0, 'rmsBETXY':0, 'maxBETX':0, 'maxBETY':0, 'maxBETXY':0}
        self.worst_permutation_both = {'id':0, 'rmsBETX':0, 'rmsBETY':0, 'rmsBETXY':0, 'maxBETX':0, 'maxBETY':0, 'maxBETXY':0}
        self.beta_wall_MADX         = {'id':0, 'rmsBETX':0, 'rmsBETY':0, 'meanBETX':0, 'meanBETY':0, 'rmsBETXY':0, 'maxBETX':0, 'maxBETY':0, 'maxBETXY':0}
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

    def get_pair_betx(self, index: int) -> Tuple[float, float]:
        raise NotImplementedError()

    def get_pair_bety(self, index: int) -> Tuple[float, float]:
        raise NotImplementedError()

    def get_pair_names(self, index: int) -> Tuple[str, str]:
        """ returns the names of the magnets of a natural pair (magnet1, magnet2)
        where 'natural' means, A/B parts of Q2, whereas the whole Q1 is paired with Q3:
            Q2a and Q2b
            Q1ab Q3ab
        """
        raise NotImplementedError()

    def check_correctability(self):
        """ checks if all values are above the `CORRECTABILITY_THRESHOLD`

            returns `True` if at least one cold mass error is above threshold
        """
        errs = [e.real_error for e in self.cold_masses]
        return np.any(np.abs(errs) > CORRECTABILITY_THRESHOLD)

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


    #@staticmethod
    #def score_sum_with_beta_as_weight(errors):
    #    """
    #    Sorts according to sum of pair elements
    #    (this favors errors that naturally cancel each other)
    #    """

    #    return np.sqrt(np.sum([
    #        (errors.get_pair(i)[0].real_error * errors.get_pair_betx(i)[0]  + errors.get_pair(i)[1].real_error * errors.get_pair_betx(i)[1])**2 
    #        for i in range(errors.get_pair_count())])**2 + np.sum([
    #        (errors.get_pair(i)[0].real_error * errors.get_pair_bety(i)[0]  + errors.get_pair(i)[1].real_error * errors.get_pair_bety(i)[1])**2 
    #        for i in range(errors.get_pair_count())])**2)


    #@staticmethod
    #def score_sum_betabeating(errors):
    #    """
    #    Sorts according to sum of pair elements
    #    (this favors errors that naturally cancel each other)
    #    """

    #    score_bbx = np.sum([
    #        ( errors.get_pair(i)[0].real_error*errors.get_pair_betx(i)[0]*errors.get_pair_sign(i)[0] 
    #        + errors.get_pair(i)[1].real_error*errors.get_pair_betx(i)[1]*errors.get_pair_sign(i)[1])
    #        for i in range(errors.get_pair_count())
    #        ])

    #    score_bby = np.sum([
    #        ( errors.get_pair(i)[0].real_error*errors.get_pair_bety(i)[0]*errors.get_pair_sign(i)[0] 
    #        + errors.get_pair(i)[1].real_error*errors.get_pair_bety(i)[1]*errors.get_pair_sign(i)[1])
    #        for i in range(errors.get_pair_count())
    #        ])

    #    return np.sqrt( score_bbx**2 + score_bby**2 )


    @staticmethod
    def score_maxabs_betabeatingX(errors,errors_bis=None,with_correction=True):
        """
        Sorts according to sum of pair elements
        (this favors errors that naturally cancel each other)
        """
        
        bbeating = errors.get_generated_betabeating(with_correction=with_correction)
        bbeating_bis = [0,0]
        if errors_bis is not None:
            bbeating_bis = errors_bis.get_generated_betabeating(with_correction=with_correction)
        
        score_bbx = np.max(abs(bbeating[0]+bbeating_bis[0]))
        #score_bby = np.max(abs(bbeating[1]+bbeating_bis[1]))

        return score_bbx


    @staticmethod
    def score_rms_betabeatingX(errors,errors_bis=None,with_correction=True):
        """
        Sorts according to sum of pair elements
        (this favors errors that naturally cancel each other)
        """
        
        bbeating = errors.get_generated_betabeating(with_correction=with_correction)
        bbeating_bis = [0,0]
        if errors_bis is not None:
            bbeating_bis = errors_bis.get_generated_betabeating(with_correction=with_correction)
        
        score_bbx = rms(bbeating[0]+bbeating_bis[0])
        #score_bby = rms(bbeating[1]+bbeating_bis[1])

        return score_bbx
        
        


    @staticmethod
    def score_maxabs_betabeatingXY(errors,errors_bis=None,with_correction=True):
        """
        Sorts according to sum of pair elements
        (this favors errors that naturally cancel each other)
        """
        
        bbeating = errors.get_generated_betabeating(with_correction=with_correction)
        bbeating_bis = [0,0]
        if errors_bis is not None:
            bbeating_bis = errors_bis.get_generated_betabeating(with_correction=with_correction)
        
        score_bbx = np.max(abs(bbeating[0]+bbeating_bis[0]))
        score_bby = np.max(abs(bbeating[1]+bbeating_bis[1]))

        return np.sqrt( score_bbx**2 + score_bby**2 )


    @staticmethod
    def score_max_betabeatingXY(errors,errors_bis=None,with_correction=True):
        """
        Sorts according to sum of pair elements
        (this favors errors that naturally cancel each other)
        """
        
        bbeating = errors.get_generated_betabeating(with_correction=with_correction)
        bbeating_bis = [0,0]
        if errors_bis is not None:
            bbeating_bis = errors_bis.get_generated_betabeating(with_correction=with_correction)
        
        score_bbx = np.max(bbeating[0]+bbeating_bis[0])
        score_bby = np.max(bbeating[1]+bbeating_bis[1])

        return np.sqrt( score_bbx**2 + score_bby**2 )


    @staticmethod
    def score_min_betabeatingXY(errors,errors_bis=None,with_correction=True):
        """
        Sorts according to sum of pair elements
        (this favors errors that naturally cancel each other)
        """
        
        bbeating = errors.get_generated_betabeating(with_correction=with_correction)
        bbeating_bis = [0,0]
        if errors_bis is not None:
            bbeating_bis = errors_bis.get_generated_betabeating(with_correction=with_correction)
        
        score_bbx = np.min(bbeating[0]+bbeating_bis[0])
        score_bby = np.min(bbeating[1]+bbeating_bis[1])

        return np.sqrt( score_bbx**2 + score_bby**2 )


    @staticmethod
    def score_rms_betabeatingXY(errors,errors_bis=None,with_correction=True):
        """
        Sorts according to sum of pair elements
        (this favors errors that naturally cancel each other)
        """
        
        bbeating = errors.get_generated_betabeating(with_correction=with_correction)
        bbeating_bis = [0,0]
        if errors_bis is not None:
            bbeating_bis = errors_bis.get_generated_betabeating(with_correction=with_correction)
        
        score_bbx = rms(bbeating[0]+bbeating_bis[0])
        score_bby = rms(bbeating[1]+bbeating_bis[1])

        return np.sqrt( score_bbx**2 + score_bby**2 )


    def sort(self, sort_fn,with_correction=True):
        """
        Performs the sorting, i.e. brute force calculates the score for each permutation and selects
        the best one.

        sort_fn should take a Pairs object and return a score

        Example:
            ``` python

            pairs = Pairs()

            def sort_fn(pairs):
                return np.sum([
                    (errors.get_pair(i)[0].real_error - errors.get_pair(i)[1].real_error)**2
                    for i in range(errors.get_pair_count())])

            pairs.sort(sort_fn)
            ```

        """
        print(f"sorting {self.NAME}") 
        print(f"searching for best combination in {len(self.permutations)} permutations")
        print(f"using the diff method")

        score = sort_fn(self,with_correction=with_correction)
        print(f" -- initial score: {score}")

        next_progress = 0.1
        best_comb = 0
        
        list_score = np.empty(len(self.permutations))
        for i in range(len(self.permutations)):
            if i / len(self.permutations) > next_progress:
                print(f"progress: {100* i / len(self.permutations):.0f} %")
                next_progress += 0.1

            self.selected_permutation = i

            sum = sort_fn(self,with_correction=with_correction)
            
            list_score[i] = sum

            if sum < score:
                score = sum
                best_comb = i
                
        self.selected_permutation = best_comb
        
        print(f"final score: {score}")
        return np.argsort(list_score)

    def sort_diff(self):
        """ Convenience method for performing the sorting according to difference"""
        return self.sort(self.score_diff)

    def sort_sum(self):
        """ Convenience method for performing the sorting according to sum"""
        return self.sort(self.score_sum)

    #def sort_beta(self):
    #    """ Convenience method for performing the sorting according to sum with respect to the betas"""
    #    return self.sort(self.score_sum_with_beta_as_weight)

    def sort_def_fn(self,with_correction=True):
        """ Convenience method for performing the sorting according to beta-beating expression"""
        return self.sort(self.score_def_fn,with_correction=with_correction)

    def sort_rmsbetabeatingXY(self,with_correction=True):
        """ Convenience method for performing the sorting according to beta-beating expression"""
        return self.sort(self.score_rms_betabeatingXY,with_correction=with_correction)

    def sort_maxbetabeatingXY(self,with_correction=True):
        """ Convenience method for performing the sorting according to beta-beating expression"""
        return self.sort(self.score_maxabs_betabeatingXY,with_correction=with_correction)

    def sort_rmsbetabeatingX(self,with_correction=True):
        """ Convenience method for performing the sorting according to beta-beating expression"""
        return self.sort(self.score_rms_betabeatingX,with_correction=with_correction)

    def sort_maxbetabeatingX(self,with_correction=True):
        """ Convenience method for performing the sorting according to beta-beating expression"""
        return self.sort(self.score_maxabs_betabeatingX,with_correction=with_correction)


class CorrectabilityError(Exception):
    def __init__(self, errs) -> None:
        self.errs = errs

    def __repr__(self) -> str:
        return (f"the following errors are all below the threshold {CORRECTABILITY_THRESHOLD}:\n"
                + ", ".join([str(x) for x in self.errs]))

    def __str__(self) -> str:
        return self.__repr__()

