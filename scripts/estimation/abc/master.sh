#!/bin/bash

cd scripts/estimation/abc
scp -r * hyak:/gscratch/csde/sjenness/stitnt
scp hyak:/gscratch/csde/sjenness/stitnt/data/*.rda data/

# Send
scp est/fit.20k.rda mox:/gscratch/csde/sjenness/stitnt/est
scp scripts/estimation/abc/data/abc.prep.rda mox:/gscratch/csde/sjenness/stitnt/data
scp scripts/estimation/abc/*.* mox:/gscratch/csde/sjenness/stitnt

# Receive
scp mox:/gscratch/csde/sjenness/stitnt/data/abc.wave*.rda scripts/estimation/abc/data/

sbatch -p ckpt -A csde-ckpt --array=1-72 --job-name=wave0 --ntasks-per-node=28 --mem=100G --time=1:00:00 --export=ALL,wave=0 runsim.sh
sbatch -p ckpt -A csde-ckpt --array=1-58 --depend=afterany:582211 --job-name=wave1 --ntasks-per-node=28 --mem=100G --time=1:00:00 --export=ALL,wave=1 runsim.sh
sbatch -p ckpt -A csde-ckpt --array=1-58 --depend=afterany:582283 --job-name=wave2 --ntasks-per-node=28 --mem=100G --time=1:00:00 --export=ALL,wave=2 runsim.sh
sbatch -p ckpt -A csde-ckpt --array=1-58 --depend=afterany:582295 --job-name=wave3 --ntasks-per-node=28 --mem=100G --time=1:00:00 --export=ALL,wave=3 runsim.sh
sbatch -p ckpt -A csde-ckpt --array=1-58 --depend=afterany:582306 --job-name=wave4 --ntasks-per-node=28 --mem=100G --time=1:00:00 --export=ALL,wave=4 runsim.sh
sbatch -p ckpt -A csde-ckpt --array=1-58 --depend=afterany:582345 --job-name=wave5 --ntasks-per-node=28 --mem=100G --time=1:00:00 --export=ALL,wave=5 runsim.sh

sbatch -p ckpt -A csde-ckpt --array=1-58 --job-name=wave6 --ntasks-per-node=28 --mem=100G --time=1:00:00 --export=ALL,wave=6 runsim.sh
sbatch -p ckpt -A csde-ckpt --array=1-58 --depend=afterany:582579 --job-name=wave7 --ntasks-per-node=28 --mem=100G --time=1:00:00 --export=ALL,wave=7 runsim.sh
sbatch -p ckpt -A csde-ckpt --array=1-58 --depend=afterany:582637 --job-name=wave8 --ntasks-per-node=28 --mem=100G --time=1:00:00 --export=ALL,wave=8 runsim.sh
sbatch -p ckpt -A csde-ckpt --array=1-58 --depend=afterany:582638 --job-name=wave9 --ntasks-per-node=28 --mem=100G --time=1:00:00 --export=ALL,wave=9 runsim.sh
sbatch -p ckpt -A csde-ckpt --array=1-58 --depend=afterany:582639 --job-name=wave10 --ntasks-per-node=28 --mem=100G --time=1:00:00 --export=ALL,wave=10 runsim.sh

sbatch -p ckpt -A csde-ckpt --array=1-58 --job-name=wave11 --ntasks-per-node=28 --mem=100G --time=1:00:00 --export=ALL,wave=11 runsim.sh
sbatch -p ckpt -A csde-ckpt --array=1-58 --depend=afterany:582879 --job-name=wave12 --ntasks-per-node=28 --mem=100G --time=1:00:00 --export=ALL,wave=12 runsim.sh
sbatch -p ckpt -A csde-ckpt --array=1-58 --depend=afterany:582890 --job-name=wave13 --ntasks-per-node=28 --mem=100G --time=1:00:00 --export=ALL,wave=13 runsim.sh
sbatch -p ckpt -A csde-ckpt --array=1-58 --depend=afterany:582901 --job-name=wave14 --ntasks-per-node=28 --mem=100G --time=1:00:00 --export=ALL,wave=14 runsim.sh
sbatch -p ckpt -A csde-ckpt --array=1-58 --depend=afterany:582912 --job-name=wave15 --ntasks-per-node=28 --mem=100G --time=1:00:00 --export=ALL,wave=15 runsim.sh


sbatch --dependency=$(squeue --noheader --format %i --name <JOB_NAME>)
