#!/bin/bash
#inf=OMI-Aura_L2-OMPIXCOR_2016m0701t1843-o63631_v003-2018m0301t214026.he5
inf=${1}
tmp=${inf}.tmp
outf=omino2.nc

# get lat/lon vars NO2 uses VIS channel
ncks -O -G : -g OMI\ Ground\ Pixel\ Corners\ VIS $inf $outf
ncks -O -v TiledCornerLatitude,TiledCornerLongitude $outf ${outf}.1
ncks -O -v Latitude,Longitude $outf ${outf}.2
ncrename -O -d phony_dim_12,nXtrack -d phony_dim_13,Ncorners -d phony_dim_14,nTimes  ${outf}.1
ncrename -O -d phony_dim_15,nTimes -d phony_dim_16,nXtrack  ${outf}.2
ncks -A ${outf}.2 ${outf}.1
mv ${outf}.1 ${outf}
rm ${outf}.2

# rename lat
ncrename -O -v Latitude,lat ${outf}
ncrename -O -a lat@Units,units  ${outf}
ncatted -h -a long_name,lat,c,c,latitude -a bounds,lat,c,c,lat_bnds -a units,lat,o,c,degrees_north  ${outf}

# shift and rename lon
ncrename -O -v Longitude,lon ${outf}
ncrename -O -a lon@Units,units  ${outf}
ncatted -h -a long_name,lon,c,c,longitude -a bounds,lon,c,c,lon_bnds -a units,lon,o,c,degrees_east  ${outf}

# rename bnds
ncrename -O -v TiledCornerLatitude,lat_bnds -v TiledCornerLongitude,lon_bnds ${outf}

# shift dim order
ncpdq -O -a nTimes,nXtrack,Ncorners  ${outf} ${outf}

# add data
#inf2=OMI-Aura_L2-OMNO2_2016m0701t1843-o63631_v003-2019m0819t165356.he5
inf2=${2}
tmpf=omino2.tmp

ncks -O -G : -g HDFEOS/SWATHS/ColumnAmountNO2 $inf2 $tmpf
ncks -O -v AmfTrop,ColumnAmountNO2,ColumnAmountNO2Trop,XTrackQualityFlags,CloudFraction,VcdQualityFlags,ScatteringWeight,ScatteringWtPressure,TropopausePressure $tmpf ${tmpf}.1
ncrename -O -d phony_dim_0,nTimes -d phony_dim_1,nXtrack -d phony_dim_2,NumSWLevels ${tmpf}.1 ${tmpf}.1
ncpdq -O -a NumSWLevels,nTimes,nXtrack ${tmpf}.1 ${tmpf}.1

ncks -O -v SolarZenithAngle,ViewingZenithAngle $tmpf ${tmpf}.2
ncrename -O -d phony_dim_7,nTimes -d phony_dim_5,nXtrack ${tmpf}.2 ${tmpf}.2

#ncks -A -v ColumnAmountNO2,XTrackQualityFlags,CloudFraction,VcdQualityFlags,ScatteringWeight,ScatteringWtPressure,SolarZenithAngle $tmpf $outf
ncks -A ${tmpf}.1 $outf
ncks -A ${tmpf}.2 $outf
rm $tmpf
rm ${tmpf}.1
rm ${tmpf}.2

# edit XTQF for masking
slicef=omino2.nc
echo vcdqualityflag
ncap2 -O -s 'XTrackQualityFlags=float(XTrackQualityFlags)' $slicef $slicef
cat > vcdqf.nco << EOF
*validcol[\$nTimes,\$nXtrack]=0;
where(XTrackQualityFlags==0 && VcdQualityFlags==0 && SolarZenithAngle<85 && CloudFraction<300) validcol=1;
validcol.ram_write();
EOF
ncap2 -O -S vcdqf.nco  $slicef $slicef
ncap2 -O -s 'where(ColumnAmountNO2<0) ColumnAmountNO2=ColumnAmountNO2.get_miss()' $slicef $slicef
ncap2 -O -s 'where(validcol==0) ColumnAmountNO2=ColumnAmountNO2.get_miss()' $slicef $slicef
ncap2 -O -s 'where(ColumnAmountNO2Trop<0) ColumnAmountNO2Trop=ColumnAmountNO2Trop.get_miss()' $slicef $slicef
ncap2 -O -s 'where(validcol==0) ColumnAmountNO2Trop=ColumnAmountNO2Trop.get_miss()' $slicef $slicef
echo done quality filter

# make cmaq grid file
qgridin=108NHEMI2_grid.nc
qgridout=cmaqgrid.nc
ncks -O -v latitude,longitude,latitude_bounds,longitude_bounds $qgridin $qgridout

satf=omino2.nc

# try to infer sat grid
#ncremap --mask_src=validcol --msk_src=XTrackQualityFlags -d $satf -g grd.nc
ncremap --mask_src=validcol -d $satf -g grd.nc

# try to infer cmaq grid
ncremap -d $qgridout -g qgrd.nc

# make map file from sat to global grid
ncremap --devnull=No -a $method -s grd.nc -g qgrd.nc -m qmap.nc

# regrid sat to global
ncremap --preserve=mean -m qmap.nc $satf dat_rgrcmaq_${method}.nc

# clean
rm vcdqf.nco omino2.nc cmaqgrid.nc grd.nc qgrd.nc qmap.nc
