#!/bin/bash

#SBATCH -p singlepe
#SBATCH -t 2:00:00
#SBATCH -n 1
#SBATCH --gid=mod3eval
#SBATCH --account=mod3eval
#SBATCH -J omitocmaq
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.out

#source ~/.bashrc


satdir=/work/ROMO/satellite/OMNO2d/aura.gesdisc.eosdis.nasa.gov/data/Aura_OMI_Level2/OMNO2.003/2016/
pixdir=/work/MOD3EVAL/jeast/satellite/OMPIXCOR/aura.gesdisc.eosdis.nasa.gov/data/Aura_OMI_Level2/OMPIXCOR.003/2016/

satf=${satdir}/183/OMI-Aura_L2-OMNO2_2016m0701t0532-o63623_v003-2018m0617t030115.he5
pixelcornerf=${pixdir}/183/OMI-Aura_L2-OMPIXCOR_2016m0701t0532-o63623_v003-2018m0301t214028.he5

concf=/work/MOD3EVAL/jeast/scratch/regrid_conservative/CCTM_CONC_v532_WWLLN_20160701.nc
met2d=/work/MOD3EVAL/jeast/scratch/regrid_conservative/METCRO2D.108NHEMI2.44L.20160701
met3d=/work/MOD3EVAL/jeast/scratch/regrid_conservative/METCRO3D.108NHEMI2.44L.20160701
outpath=./

./regrid_cmaq_filter.sh $pixelcornerf $satf

python processOMI.py $satf dat_rgrcmaq_nco_con.nc $concf $met2d $met3d $outpath
