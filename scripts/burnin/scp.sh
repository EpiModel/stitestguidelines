#!/bin/bash

## MOX

# Send
scp est/*.rda mox:/gscratch/csde/sjenness/stitnt/est
scp scripts/burnin/*.* mox:/gscratch/csde/sjenness/stitnt

# Receive
scp mox:/gscratch/csde/sjenness/stitnt/data/*.rda data/


## IKT

# Send
scp est/*.rda hyak:/gscratch/csde/sjenness/stitnt/est
scp scripts/burnin/*.* hyak:/gscratch/csde/sjenness/stitnt

# Receive
scp hyak:/gscratch/csde/sjenness/stitnt/data/*.rda data/
