#!/usr/bin/env python
# coding: utf-8

# ### Process OMI VCD for model comparison after regridding
#     James East
#     2021-09-27
#     
# Needs to do:
# * Get model data at right time step (nearest neighbor)
# * Get model pressure levels
# * interpolate SW to model, accounting for tropopause
# * Compute model VCD (simple) and partial columns
#     * $\Omega_V^M$[ppm] = $\frac{MW_{air}*Av*1e-4}{1e6} \sum_{surface}^{tp} \frac{\rho (z)}{\Delta H (z)}*C_{NO_2}(z)$
# * Apply SW to model partial columns and sum for model SCD
#     * $\Omega_S^M$ = $\sum_{surface}^{tp} \Omega_V^M(z)*w(z)$ 
# * Compute model AMF (see Cooper et al. 2020 ACP)
#     * $\text{AMF}_M=M_G * \sum_{surface}^{tp}w(z)\frac{n(z)}{\Omega_V^M}$
#     * $M_G=\sec \theta_o + \sec \theta$
#     * $\theta_o$ is solar zenith angle
#     * $\theta$ is viewing zenith angle
#     * see OMNO2G README p. 28 (https://aura.gesdisc.eosdis.nasa.gov/data/Aura_OMI_Level2/OMNO2.003/doc/README.OMNO2.pdf)
# * Compute satellite VCD with CMAQ profile (see Duncan et al. 2014 Section S.5 and Palmer et al. 2001 eqs 13-15 and discussion)
#     * $\Omega_V^{O'} = \frac{\Omega_V^O * \text{AMF}_{\text{OMI}}}{\text{AMF}_M}$
# * Add OMI VCD with CMAQ profile to file
# * Add CMAQ VCD & SCD to file
# * Change dims to match IOAPI expectations for easier future comparison
# * 1 retrieval at a time for now, loop later

# Assumptions:
# * Tropopause cutoff based on OMI tropopause definition
# * All sat fields regridded to CMAQ grid first
# * Surface pressure set to CMAQ for SW application
# * Cloud filtering based on OMI cloud fraction



import xarray as xr
import numpy as np
import PseudoNetCDF as pnc
import cftime
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import sys




# original sat file: OMI-Aura_L2-OMNO2_2016m0701t1843-o63631_v003-2019m0819t165356.he5
# horizontal regridded sat file: dat_rgrcmaq_nco_con.nc
# CMAQ file: CCTM_CONC_v532_WWLLN_20160701.nc

# input files:
satfile = sys.argv[1] 
regrid_file = sys.argv[2] #'./dat_rgrcmaq_nco_con.nc'
concfile = sys.argv[3] #'./CCTM_CONC_v532_WWLLN_20160701.nc'
met2d = sys.argv[4]
met3d = sys.argv[5]
outpath = sys.argv[6]

def get_model_tidx(satf):
    '''get sat time brute force way
    returns indices of TSTEPs to use
    use higher index, we only look at northern hemisphere'''
    
    omiraw=xr.open_dataset(satf,group='/HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields')
    otime=omiraw['Time'].rename({'phony_dim_7':'nTimes'})
    t = cftime.num2date(otime[:]-32-4, units='seconds since 1993-01-01T00:00:00Z')

    timevf = pnc.PseudoNetCDFFile()
    timevf.createDimension('time', 25)
    timev = timevf.createVariable('time', 'd', ('time',))
    timev.units = f'hours since {t[0].strftime("%Y-%m-%d 00+0000Z")}'
    timev[:] = np.arange(25)
    # return tidx for sat swath shape
    tidx = cftime.date2index(t, nctime=timev, select='nearest')[:, None].repeat(omiraw.phony_dim_5.shape, 1)
    return t,np.unique(tidx)[-1]

def makevar(omirgr,var,name,units):
    gridvar = np.full(satvalid.shape, np.nan)
    gridvar[satvalid] = var
    outvar = xr.DataArray(
        data = gridvar,
        coords = omirgr.ColumnAmountNO2Trop.coords,
        dims = omirgr.ColumnAmountNO2Trop.dims,
        name = name,
        attrs = {
            'units':units
        }
    )
    return outvar

def droplays2(arr1, arr2):
    '''If repeated values at bottom of arr, get rid of repeats
    and append 0 to the end of the arr
    Also changes arr2 in same way, based only on arr1 condition'''
    if arr1[0]==arr1[1]:
        arr2 = np.append(np.delete(arr2,0),0)
        arr1, arr2 = droplays2(np.append(np.delete(arr1, 0), 0), arr2)
    else:
        pass
    return arr1, arr2

def interpno2levs2omi(no2in): 
    no2in_cum = np.insert(np.cumsum(no2in, axis=-1),0,0,axis=-1) # cumulative no2 at edges
    no2out = np.zeros(opresi.shape)
    for iloc in range(no2in.shape[0]):
        ifunc1 = interp1d(LAYi,no2in_cum[iloc,:],fill_value='extrapolate')
        no2out[iloc,:] = ifunc1(olidx[iloc])
    return no2out


t,tidx=get_model_tidx(satfile)

gridname = '108NHEMI2'
satdate=satfile.split('-')[2].split('_')[-1]
outfname = f'{outpath}/./OMNO2-L3_{gridname}_{satdate}.nc'

# open regrid file and get valid cells
omirgr = xr.open_dataset(regrid_file)
satvalid = omirgr['validcol'].values.astype(bool)
jvalid,ivalid=np.where(satvalid) # i=(COL,nXtrack), j=(ROW,nTimes)
satvcdtrop = omirgr['ColumnAmountNO2Trop'].values[jvalid,ivalid]
sw = omirgr['ScatteringWeight'].values[:,jvalid,ivalid]
sw = np.transpose(sw,[1,0])
tppres = omirgr['TropopausePressure'].values[jvalid,ivalid]*100 #convert hPa --> Pa
swpres = omirgr['ScatteringWtPressure'].values*100 #convert hPa --> Pa
satamftrop = omirgr['AmfTrop'].values[jvalid,ivalid]
sza = omirgr['SolarZenithAngle'].values[jvalid,ivalid]
vza = omirgr['ViewingZenithAngle'].values[jvalid,ivalid]

# apply filter to ColumnAmountNO2
omirgr['ColumnAmountNO2Trop'].values = np.where(
    satvalid,
    omirgr['ColumnAmountNO2Trop'].values,
    omirgr['ColumnAmountNO2Trop'].MissingValue
)

omirgr['ColumnAmountNO2'].values = np.where(
    satvalid,
    omirgr['ColumnAmountNO2'].values,
    omirgr['ColumnAmountNO2'].MissingValue
)

# get model files
qf=xr.open_dataset(concfile)

m2f = xr.open_dataset(met2d)

m3f = xr.open_dataset(met3d)

gf = pnc.pncopen('/home/jeast/GRIDDESC',format='griddesc')

# Get model vars
qno2 = qf.variables['NO2'].values[tidx, :, jvalid, ivalid]
qdens = m3f.variables['DENS'].values[tidx, :, jvalid, ivalid]
qzf = m3f.variables['ZF'].values[tidx, :, jvalid, ivalid]
qprsfc=m2f.variables['PRSFC'].values[tidx, :, jvalid, ivalid]

# Calc model pressure
# cmaq pres edges, valid only
qpedges = (
    (qprsfc - gf.VGTOP) *
    qf.VGLVLS[None,:] + qf.VGTOP
)
qdp = -np.diff(qpedges, axis=1) / 100 #hPa

# model dz
qdz = qzf.copy()
qdz[..., 1:] = qzf[..., 1:] - qzf[..., :-1]

# Calc model VCD and partial columns
MWair = 0.0289628
nair = qdens / MWair * qdz # mole / m2
qno2prof = qno2 / 1e6 * nair * 6.022e23 * 1e-4 # ppm_i * mole_air/m2 * Av * m2/cm2 = molec_i / cm2
qvcd_simple = qno2prof.sum(-1) #not accounting for tropopause

# Interpolate SW to model and apply in tropopause

# make empty array for pressure interfaces
olevsn = swpres.shape[-1]
opresi = np.zeros(satvcdtrop.shape)
opresi = opresi[:,None].repeat(olevsn+1,axis=-1) # sat interface pressures

# Set OMI scattering weight interface P like GSI (omimod.f90)
opresi[:,0] = qprsfc[:,0] # set to cmaq lowest
opresi[:,-1] = swpres.min() # top P 0
for lev in range(swpres.size-1):
    opresi[:,lev+1] = 0.5*(swpres[lev]+swpres[lev+1]) # interface Ps avg of midpt Pressures

# Set pressures like GSI 
for k in range(1, swpres.size):
    opresi[:,k] = np.where(
        opresi[:,k] < qpedges[:,0], # where omi P less than lowest lev cmaq P
        opresi[:,k], # set to TM4 P
        qpedges[:,0] # other set to lowest lev cmaq P
    )


# Want bottom pressure to be CMAQ pressure,
# Then next pressure to be lowest OMI pressure
# which is less than CMAQ bottom pressure

npts = opresi.shape[0]
for iloc in range(npts):
    opresi[iloc,:], sw[iloc,:] = droplays2(opresi[iloc,:], sw[iloc,:])

# array to fill with OMI layer indices relative to CMAQ lays
olidx = np.zeros(opresi.shape)

# interpolate omi pressures to get exact lev indices relative to cmaq lev indices
krange = np.arange(qpedges.shape[1]) #CMAQ LAY+1 indices

npts = opresi.shape[0]
for iloc in range(npts):
    ifunc = interp1d(qpedges[iloc], krange, fill_value='extrapolate')
    olidx[iloc,:] = ifunc(opresi[iloc])


# Interpolate to OMI levels
LAYi=np.arange(qpedges.shape[1]) #CMAQ LAY+1 indices

# interpolate integrated profile
qno2_olevsi = interpno2levs2omi(qno2prof)

# get difference bt integrated levels for profile
qno2_olevs_prof = np.diff(qno2_olevsi, axis=-1)


# Apply scattering weights
qno2profsw = np.where(
    (opresi[:,:-1] < (tppres[:,None])),
    0, # zero above TM4 tropopause
    qno2_olevs_prof * sw # model profile * SW in troposphere
)

qno2prof = np.where(
    (opresi[:,:-1] < (tppres[:,None])),
    0, # zero above TM4 tropopause
    qno2_olevs_prof # model profile
)

qscd_sw = qno2profsw.sum(-1)
qvcd = qno2prof.sum(-1)


# compute model AMF
amfgeo = 1/np.cos(np.deg2rad(vza)) + 1/np.cos(np.deg2rad(sza))
qamf = amfgeo/qvcd*qscd_sw
qscd = qvcd*qamf
#qamf = qscd / qvcd

# compute OMI VCD with CMAQ profile
satvcdqprof = satvcdtrop * satamftrop / qamf

# Set up file for saving to disk
keepvars = [
    satvcdqprof,
    qamf,
    qscd,
    qvcd
]

keepname = [
    'OMNO2VCDTrop_modprof',
    'Mod_amf',
    'ModNO2SCDTrop',
    'ModNO2VCDTrop'
]

keepunits = [
    'molecules cm-2',
    'molecules cm-2',
    'molecules cm-2',
    'unitless'
]


omiout = xr.merge([makevar(omirgr,var,name,units) for (var,name,units) in zip(keepvars,keepname,keepunits)]+[omirgr])

omiout=omiout.rename(
    {
        'nXtrack':'COL',
        'nTimes':'ROW',
        'NumSWLevels':'LAY'
    }
).expand_dims({'TSTEP':1})
omiout = omiout.transpose('TSTEP','LAY','ROW','COL','nvertices')


# Assign attrs to file
keepattrs = [
    'IOAPI_VERSION','EXEC_ID','FTYPE','CDATE','CTIME',
    'WDATE','WTIME','SDATE','STIME','TSTEP','NTHIK',
    'NCOLS','NROWS','GDTYP','P_ALP','P_BET','P_GAM',
    'XCENT','YCENT','XORIG','YORIG','XCELL','YCELL'
]

qattrs={attkey:qf.attrs[attkey] for attkey in keepattrs}

satattrs = {
    'Vertical_grid_note':'Vertical grid inherited from OMI NO2 L2 file',
    'Satellite_file':satfile.split('/')[-1],
    'STIME':np.int32(tidx)
}

omiout = omiout.assign_attrs(**qattrs)
omiout = omiout.assign_attrs(**satattrs)
omiout = omiout.assign_attrs(omirgr.attrs)

omiout.to_netcdf(outfname)
