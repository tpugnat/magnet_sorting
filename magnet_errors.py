from random import random
import numpy as np

# ---- helper functions ----------------------------------------------------------------------------

# Set the seed
#print(f"Numpy random generator seed: {np.random.seed()}\n")
TYPE_RAND = 'Gaussian';

# what does this do?
def sortOnSum(val):
    return sum(val)

def rand(amplitude: float, cut_3sigma= False):
    if TYPE_RAND == "Uniform":
        #print(f"Uniform distribution pm{amplitude}")
        ran = (random() - 0.5) * 2 * amplitude
        #ran = (np.random.random() - 0.5) * 2 * amplitude
    if TYPE_RAND == "Gaussian":
        #print(f"Gaussian distribution sigma={amplitude/3}  (3sigma cut: {cut_3sigma})")
        ran = (amplitude/3) * np.random.normal()
        if cut_3sigma:
            while abs(ran)> 3*(amplitude/3):
                ran = (amplitude/3) * np.random.normal()
    return ran
    #return normal(0, amplitude)
    

def uniform(amplitude: float):
    ran = (random() - 0.5) * 2 * amplitude
    return ran


def gauss3sc(amplitude: float, cut_3sigma= False):
    ran = (amplitude/3) * np.random.normal()
    if cut_3sigma:
        while abs(ran)> 3*(amplitude/3):
            ran = (amplitude/3) * np.random.normal()
    return ran
    


#SYST_ERROR = rand(10.0) # std deviation of systematic errors
SYST_ERROR = 10.0


# class MagnetPair:
#     """pair of errors for a pair of magnets

#     Magnets inside pair can be swapped by calling pair.swap_magnets()
#     And entire pairs can be swapped by calling MagnetPair.swap_pairs(pair1, pair2)
#     """
#     def __init__(self,
#                  amp_real: float,
#                  amp_meas: float,
#                  name_a: str,
#                  name_b: str
#                  ) -> None:
#         self.A = MagnetError(amp_real, amp_meas)
#         self.B = MagnetError(amp_real, amp_meas)
#         self.name_A = name_a
#         self.name_B = name_b

#     def powered(self) -> float:
#         return (self.A.believed_error() + self.B.believed_error()) / 2

#     def powered_a(self) -> float:
#         """returns the pairpowered value for magnet A as in Hectors script
#         to get the error remaining after corrections, use `-0.5 * pair.powered_a()`
#         """
#         return -self.A.real_error + self.powered()

#     def powered_b(self) -> float:
#         """returns the pairpowered value for magnet B as in Hectors script
#         to get the error remaining after corrections, use `-0.5 * pair.powered_b()`
#         """
#         return -self.B.real_error + self.powered()

#     def swap_magnest(self):
#         temp = self.A
#         self.A = self.B
#         self.B = temp

#     def log_strengths(self, str_map):
#         str_map[self.name_A] = self.A.real_error
#         str_map[self.name_B] = self.B.real_error

#     def write_real_error(self, error_file):
#         # error_file.write(qgrouped[-1] '\n')
#         MagnetPair.write_error_to_file(error_file, self.name_A, self.A.real_error)
#         MagnetPair.write_error_to_file(error_file, self.name_B, self.B.real_error)
#         error_file.write('\n')

#     def write_corr_error(self, error_file):
#         """writes the corrected error,
#         i.e. the error remaining in the machine after corrections
#         """
#         # error_file.write(qgrouped[-1] '\n')
#         MagnetPair.write_error_to_file(error_file, self.name_A, - 0.5 * self.powered_a())
#         MagnetPair.write_error_to_file(error_file, self.name_B, - 0.5 * self.powered_b())
#         error_file.write('\n')

#     def __getitem__(self, index):
#         if index == 0:
#             return self.A
#         elif index == 1:
#             return self.B
#         else:
#             raise IndexError

#     def __setitem__(self, index, item):
#         if index == 0:
#             self.A = item
#         elif index == 1:
#             self.B = item
#         else:
#             raise IndexError

#     @staticmethod
#     def write_all_errors_to_file(filename: str, errors):
#         with open(filename, 'w') as error_file:
#             for pair in errors:
#                 pair.write_real_error(error_file)

#     @staticmethod
#     def write_error_to_file(error_file, name, value):
#         error_file.write('select, flag=error, clear;\n')
#         error_file.write(f'select, flag=error, PATTERN={name};\n' )
#         #error_file.write('Efcomp, radius = 1, order= 1, dknr={0, ' + str(unit(value)) + '};\n')
#         error_file.write('Efcomp, radius = 0.050, order= 1, dknr={0, ' + str(unit(value)) + '};\n')
#         error_file.write('\n')


#     @staticmethod
#     def swap_pairs(pair1: "MagnetPair", pair2: "MagnetPair"):
#         temp = pair1
#         pair1 = pair2
#         pair2 = temp


# def unit(x: float) -> float:
#     return x * 1.0e-4

class MagnetError:
    def __init__(self, ampl_meas=0, ampl_cali=0, ampl_prec=0) -> None:
        #self.real_error: float = rand(ampl_real, cut_3sigma= True)
        self.meas_error: float = rand(ampl_meas, cut_3sigma= True)
        self.cali_error: float = ampl_cali
        self.pres_error: float = gauss3sc(ampl_prec, cut_3sigma= True)

    #def believed_error(self) -> float:
    #    return self.real_error + self.meas_error + SYST_ERROR

    @property
    def real_error(self) -> float:
        return self.meas_error - self.cali_error - self.pres_error
