#!/bin/bash

# Send
scp est/*.rda mox:/gscratch/csde/sjenness/stitnt/est
scp scripts/burnin/abc/data/abc.prep.rda mox:/gscratch/csde/sjenness/stitnt/data
scp scripts/burnin/abc/*.* mox:/gscratch/csde/sjenness/stitnt

# Receive
scp mox:/gscratch/csde/sjenness/stitnt/data/*.rda scripts/burnin/abc/data/
scp mox:/gscratch/csde/sjenness/stitnt/data/abc.wave27.rda scripts/burnin/abc/data/


cp scripts/burnin/abc/data/abc.wave27.rda est/
