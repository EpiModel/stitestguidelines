#!/bin/bash

# Send
# scp est/*.rda hyak:/gscratch/csde/sjenness/stitnt/est
scp scripts/burnin/*.* hyak:/gscratch/csde/sjenness/stitnt

# Receive
scp hyak:/gscratch/csde/sjenness/stitnt/data/*.rda data/
