#!/bin/bash

PYTHIA8_BASE=$(scram tool info pythia8 | grep "PYTHIA8_BASE=" | cut -d'=' -f2)
FASTJET_BASE=$(scram tool info fastjet | grep "FASTJET_BASE=" | cut -d'=' -f2)

g++ -o sps_matching sps_matching.cpp \
    -I$PYTHIA8_BASE/include \
    -I$PYTHIA8_BASE/include/Pythia8Plugins \
    -I$FASTJET_BASE/include \
    -L$PYTHIA8_BASE/lib -lpythia8 \
    -L$FASTJET_BASE/lib -lfastjet \
    `root-config --cflags --libs` \
    -std=c++17 \
    -O2
