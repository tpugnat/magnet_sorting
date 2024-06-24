#!/bin/bash


MAGNET_SORTING="$(pwd -P)"


# Manage lxplus linux version
locfunc=`cat /etc/*-release | grep -i "VERSION_ID="`
if [[ $locfunc == *"7"* ]]; then
    echo " * LOAD CENTOS7 ENV"
    source /cvmfs/sft.cern.ch/lcg/views/LCG_104a/x86_64-centos7-gcc11-opt/setup.sh
    LINUX_VERSION=7
elif [[ $locfunc == *"8"* ]]; then
    echo " * LOAD CENTOS8 ENV"
    source /cvmfs/sft.cern.ch/lcg/views/LCG_104a/x86_64-centos8-gcc11-opt/setup.sh
    LINUX_VERSION=8
else
    echo " * LOAD EL9 ENV"
    source /cvmfs/sft.cern.ch/lcg/views/LCG_104a/x86_64-el9-gcc11-opt/setup.sh
    LINUX_VERSION=9
fi

# --- setup python ---------------------------------------------------------------------------------
if [[ ! -d venv${LINUX_VERSION} ]]; then
    python3 -m venv venv${LINUX_VERSION} #--no-site-packages func-adl-servicex servicex
fi
. venv${LINUX_VERSION}/bin/activate

pip install matplotlib numpy
#pip install 'pandas~=1.0'
#pip install 'idna==2.10'

if [[ ! -d omc3 ]]; then
    git clone git@github.com:pylhc/omc3.git

    cd omc3
    git checkout enhancement/292/modelcreation

    cd ..
fi
#pip install -e omc3

# --- setup model output ---------------------------------------------------------------------------
mkdir model

cd model
ln -s ${MAGNET_SORTING}/macros macros
git clone ssh://git@gitlab.cern.ch:7999/acc-models/acc-models-lhc.git
cd acc-models-lhc
git checkout hl16
touch knobs.madx
cd ../..

mkdir model1
cd model1
ln -s ${MAGNET_SORTING}/macros macros
ln -s ${MAGNET_SORTING}/model/acc-models-lhc acc-models-lhc
touch knobs.madx
cd ..

echo "pwd = ${MAGNET_SORTING}"
