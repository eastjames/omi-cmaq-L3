### Process OMI and CMAQ NO2 columns
Produces apples-to-apples comparison between model and satellite columns

### Process OMI VCD for model comparison after regridding
    James East
    2021-09-27
    
Needs to do:
* Get model data at right time step (nearest neighbor)
* Get model pressure levels
* interpolate SW to model, accounting for tropopause
* Compute model VCD (simple) and partial columns
    * $\Omega_V^M$[ppm] = $\frac{MW_{air}*Av*1e-4}{1e6} \sum_{surface}^{tp} \frac{\rho (z)}{\Delta H (z)}*C_{NO_2}(z)$
* Apply SW to model partial columns and sum for model SCD
    * $\Omega_S^M$ = $\sum_{surface}^{tp} \Omega_V^M(z)*w(z)$ 
* Compute model AMF (see Cooper et al. 2020 ACP)
    * $\text{AMF}_M=M_G * \sum_{surface}^{tp}w(z)\frac{n(z)}{\Omega_V^M}$
    * $M_G=\sec \theta_o + \sec \theta$
    * $\theta_o$ is solar zenith angle
    * $\theta$ is viewing zenith angle
    * see OMNO2G README p. 28 (https://aura.gesdisc.eosdis.nasa.gov/data/Aura_OMI_Level2/OMNO2.003/doc/README.OMNO2.pdf)
* Compute satellite VCD with CMAQ profile (see Duncan et al. 2014 Section S.5 and Palmer et al. 2001 eqs 13-15 and discussion)
    * $\Omega_V^{O'} = \frac{\Omega_V^O * \text{AMF}_{\text{OMI}}}{\text{AMF}_M}$
* Add OMI VCD with CMAQ profile to file
* Add CMAQ VCD & SCD to file
* Change dims to match IOAPI expectations for easier future comparison
* 1 retrieval at a time for now, loop later

Assumptions:
* Tropopause cutoff based on OMI tropopause definition
* All sat fields regridded to CMAQ grid first
* Surface pressure set to CMAQ for SW application
* Cloud filtering based on OMI cloud fraction
* Conservative horizontal regridding (requires nco >= 5.0.1)
