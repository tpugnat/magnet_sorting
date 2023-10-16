#!/bin/bash

python3 -m venv venv
. venv/bin/activate

pip install omc3 matplotlib numpy

cd model
ln -s ../macros macros
git clone ssh://git@gitlab.cern.ch:7999/acc-models/acc-models-lhc.git
cd acc-models-lhc
git checkout hl16
cd ../..


