# Intro

## Setup

The script `bootstrap.sh` should get you up and running

``` shell
chmod u+x bootstrap.sh
/bootstrap.sh
```

## Structure

- `do_global_corrections.py` will run a full set of simulations
  go into the file and check it out
- `pairs.py`, `q1_errors.py` `q2_errors.py` and `magnet_errors.py` define helpful classes for
  the setup of magnet errors
- `analysis_pairing.py` does the statistics on many simulations
- `job*.madx` are various job files to create the madx lattice and produce model and meas output
- the folder `macros` contains our omc3 macros for lattice and model creation
- some folders are created by `bootstrap.sh`

## Quick links

- [History](#history) gives an overview over studies already done
- [FiDeL homepage](https://lhc-div-mms.web.cern.ch/tests/MAG/Fidel/)
- [Magnet sorting Decision](https://edms.cern.ch/document/LHC-LM-ED-0001/1.0/approvalAndComments)
- [MADX Lattice](#madx-lattice) info

## Current progress

- defined steps to take with Massimo
- [Error distribution](#error-distribution) defines how the errors are created
- run many [Simulations](#simulations) with different magnet errors
- get global corrections from [Full Response](#full-response)

### Many simulations result

- some correlation visible for high bbeat
- low bbeat at Q2a + Q2b less significant

### Global corrections sims

- shows good improvement for bbeat by sorting and pairing
- shows also good improvement for corr
- both seem correlated (luckily, now we can fix both in one go)

###  run simulations with 'idealised' corrections: 

- i.e. calc corrections by hand / using mad-x
- plug I into transfer functions
- repeat tracking and measure betabeat
  
### History

- Jaime: study on correctability, long time ago. showed that a tiny deviation from optimal could
  already be a problem.
- Global correction: not good for triplet (global factor?)
- Hectors work on magnet sorting:
  [first talk](https://indico.cern.ch/event/966436/),
  [second talk](https://indico.cern.ch/event/981399/contributions/4133558/attachments/2155655/3635867/HL_inner_triplet_prediction.pdf) (with ML),
  [third talk](https://indico.cern.ch/event/1019508/contributions/4278837/attachments/2211355/3742617/Reinforcement_learning_for_IR_local_optics_correction%20(1).pdf) (more ML),

### Possible Tasks

- create a global corrections script, taking a global factor into account
- corrections from twiss
- add parameter to FR? YES
- check which correction method is the best (SbS, global+factor, APJ?, ML?)    
- decide if sorting is needed after all

### Parameters

- measured magnet strengths
- aperture?
- beta beating
- required correction strengths (pc limits?, K limits?)

### Immediate tasks:

- check corrections.py - can we get corrections from twiss?
  A: YES   
- start some simulations, get overview of beta beating
  A: 6-10%

### result from meeting with massimo

#### Constraints

- probably only Q2 pairing possible,
- possibly Q1, Q2, Q3 sorting between different IPs same side
  
#### Notation

Pairing: A - B
Sorting: Q1IP1 - Q1IP5

#### Merit

- sorting: could improve beta beating
- pairing: could improve correctability
  
  Aperture / tolerances: try simulating sorting for aperture

# Work Details

## Error distribution

The error distributions are created and managed by subclasses of [`Pairs`](pairs.py).

Example:

``` python
q1_errors = Q1Pairs(real_error=10, meas_error=2, stage=1)
q1_errors.write_errors_to_file("model/errors_Q1.madx")
```

This creates a set of error distributions `q1_errors` and writes the errors of the initial
distribution to a file `model/errors_Q1.madx`.

Now we can perform sorting:
```python
q1_errors.sort_sum()
q1_errors.write_errors_to_file("model/errors_Q1_sorted.madx")
```

Or define our own sorting function:

```python
def sort_fn(pairs):
    return np.sum([
        (errors.get_pair(i)[0].real_error - errors.get_pair(i)[1].real_error)**2
        for i in range(errors.get_pair_count())])

q1_errors.sort(sort_fn)
q1_errors.write_errors_to_file("model/errors_Q1_sorted.madx")
```




# MADX lattice


root: [acc-models-lhc](#acc-models-lhc)



## Triplet


Q1ab, Q2a, Q2b, Q3ab

defined in [acc-models-lhc/hllhc_sequence.madx:420](#acc-models-lhchllhcsequencemadx420)


###  Powering

task: please check if Q2a/b might be powered individually

