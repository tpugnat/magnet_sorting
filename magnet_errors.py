from random import random

# ---- helper functions ----------------------------------------------------------------------------

# what does this do?
def sortOnSum(val):
    return sum(val)

def rand(amplitude: float):
    return (random() - 0.5) * 2 * amplitude
    #return normal(0, amplitude)


SYST_ERROR = rand(10.0) # std deviation of systematic errors


class MagnetPair:
    """pair of errors for a pair of magnets

    Magnets inside pair can be swapped by calling pair.swap_magnets()
    And entire pairs can be swapped by calling MagnetPair.swap_pairs(pair1, pair2)
    """
    def __init__(self,
                 amp_real: float,
                 amp_meas: float,
                 name_a: str,
                 name_b: str
                 ) -> None:
        self.A = MagnetError(amp_real, amp_meas)
        self.B = MagnetError(amp_real, amp_meas)
        self.name_A = name_a
        self.name_B = name_b

    def powered(self) -> float:
        return (self.A.believed_error() + self.B.believed_error()) / 2

    def powered_a(self) -> float:
        """returns the pairpowered value for magnet A as in Hectors script
        to get the error remaining after corrections, use `-0.5 * pair.powered_a()`
        """
        return -self.A.real_error + self.powered()

    def powered_b(self) -> float:
        """returns the pairpowered value for magnet B as in Hectors script
        to get the error remaining after corrections, use `-0.5 * pair.powered_b()`
        """
        return -self.B.real_error + self.powered()

    def swap_magnest(self):
        temp = self.A
        self.A = self.B
        self.B = temp

    def log_strengths(self, str_map):
        str_map[self.name_A] = self.A.real_error
        str_map[self.name_B] = self.B.real_error

    def write_real_error(self, error_file):
        # error_file.write(qgrouped[-1] '\n')
        MagnetPair.write_error_to_file(error_file, self.name_A, self.A.real_error)
        MagnetPair.write_error_to_file(error_file, self.name_B, self.B.real_error)
        error_file.write('\n')

    def write_corr_error(self, error_file):
        """writes the corrected error,
        i.e. the error remaining in the machine after corrections
        """
        # error_file.write(qgrouped[-1] '\n')
        MagnetPair.write_error_to_file(error_file, self.name_A, - 0.5 * self.powered_a())
        MagnetPair.write_error_to_file(error_file, self.name_B, - 0.5 * self.powered_b())
        error_file.write('\n')

    def __getitem__(self, index):
        if index == 0:
            return self.A
        elif index == 1:
            return self.B
        else:
            raise IndexError

    def __setitem__(self, index, item):
        if index == 0:
            self.A = item
        elif index == 1:
            self.B = item
        else:
            raise IndexError

    @staticmethod
    def write_all_errors_to_file(filename: str, errors):
        with open(filename, 'w') as error_file:
            for pair in errors:
                pair.write_real_error(error_file)

    @staticmethod
    def write_error_to_file(error_file, name, value):
        error_file.write('select, flag=error, clear;\n')
        error_file.write(f'select, flag=error, PATTERN={name};\n' )
        error_file.write('Efcomp, radius = 1, order= 1, dknr={0, ' + str(unit(value)) + '};\n')
        error_file.write('\n')


    @staticmethod
    def swap_pairs(pair1: "MagnetPair", pair2: "MagnetPair"):
        temp = pair1
        pair1 = pair2
        pair2 = temp


def unit(x: float) -> float:
    return x * 1.0e-4

class MagnetError:
    def __init__(self, ampl_real: float, ampl_meas: float) -> None:
        self.real_error: float = rand(ampl_real)
        self.meas_error: float = rand(ampl_meas)

    def believed_error(self) -> float:
        return self.real_error + self.meas_error + SYST_ERROR
