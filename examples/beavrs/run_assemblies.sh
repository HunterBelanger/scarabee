#!/usr/bin/bash

export SCARABEE_ND_LIBRARY="/home/hunter/Documents/nuclear_data/scarabee/endf8_shem281.h5"

python assembly_16_0.py
wait
echo

python assembly_24_0.py
wait
echo

python assembly_24_12.py
wait
echo

python assembly_24_16.py
wait
echo

python assembly_31_0.py
wait
echo

python assembly_31_6.py
wait
echo

python assembly_31_15.py
wait
echo

python assembly_31_16.py
wait
echo

python assembly_31_20.py
wait
echo

python core.py
