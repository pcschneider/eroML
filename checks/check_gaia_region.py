from astropy.io import fits as pyfits
import matplotlib.pyplot as plt
import healpy as hp
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord

nside=32
hpix = 1210
step = 1


resol = hp.nside2resol(nside, arcmin=True)
pix_center = hp.pix2ang(nside, hpix, nest=True)
ra, dec = pix_center[1]/np.pi*180, 90. - pix_center[0]/np.pi*180

ww= u.Quantity(resol, u.arcmin) # + 2*edge*u.arcmin
hh = u.Quantity(resol, u.arcmin)# + 2*edge*u.arcmin

width = ww.to(u.degree).value
height = hh.to(u.degree).value

query_str = str("SELECT gaia_source.source_id,gaia_source.ra,gaia_source.ra_error,gaia_source.dec,gaia_source.dec_error,gaia_source.parallax,gaia_source.parallax_error,gaia_source.phot_g_mean_mag,gaia_source.bp_rp,gaia_source.radial_velocity,gaia_source.radial_velocity_error,gaia_source.phot_bp_mean_mag,gaia_source.phot_rp_mean_mag,gaia_source.phot_g_mean_flux_over_error,gaia_source.phot_rp_mean_flux_over_error,gaia_source.phot_bp_mean_flux_over_error,gaia_source.phot_variable_flag,gaia_source.teff_val,gaia_source.a_g_val,gaia_source.visibility_periods_used,gaia_source.astrometric_chi2_al,gaia_source.astrometric_n_good_obs_al, gaia_source.astrometric_excess_noise,gaia_source.phot_bp_rp_excess_factor  FROM gaiadr2.gaia_source WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),BOX('ICRS', %f, %f, %f, %f))=1;" % (ra, dec, width,height))   

print("Query1: \n",query_str)
print()

b = hp.boundaries(nside, hpix, nest=True, step=step).transpose()
v = hp.pixelfunc.vec2ang(b)
x, y = v[1]/np.pi*180, 90.-v[0]/np.pi*180  




tpl = (x[0], y[0], x[1], y[1], x[2], y[2], x[3], y[3])

query_str = str("SELECT gaia_source.source_id,gaia_source.ra,gaia_source.ra_error,gaia_source.dec,gaia_source.dec_error,gaia_source.parallax,gaia_source.parallax_error,gaia_source.phot_g_mean_mag,gaia_source.bp_rp,gaia_source.radial_velocity,gaia_source.radial_velocity_error,gaia_source.phot_bp_mean_mag,gaia_source.phot_rp_mean_mag,gaia_source.phot_g_mean_flux_over_error,gaia_source.phot_rp_mean_flux_over_error,gaia_source.phot_bp_mean_flux_over_error,gaia_source.phot_variable_flag,gaia_source.teff_val,gaia_source.a_g_val,gaia_source.visibility_periods_used,gaia_source.astrometric_chi2_al,gaia_source.astrometric_n_good_obs_al, gaia_source.astrometric_excess_noise,gaia_source.phot_bp_rp_excess_factor  FROM gaiadr2.gaia_source WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),POLYGON('ICRS', %f, %f, %f, %f, %f, %f, %f, %f))=1;" % tpl)   
print("Query2: \n",query_str)
print()



query_str = str("SELECT gaia_source.source_id,gaia_source.ra,gaia_source.ra_error,gaia_source.dec,gaia_source.dec_error,gaia_source.parallax,gaia_source.parallax_error,gaia_source.phot_g_mean_mag,gaia_source.bp_rp,gaia_source.radial_velocity,gaia_source.radial_velocity_error,gaia_source.phot_bp_mean_mag,gaia_source.phot_rp_mean_mag,gaia_source.phot_g_mean_flux_over_error,gaia_source.phot_rp_mean_flux_over_error,gaia_source.phot_bp_mean_flux_over_error,gaia_source.phot_variable_flag,gaia_source.teff_val,gaia_source.a_g_val,gaia_source.visibility_periods_used,gaia_source.astrometric_chi2_al,gaia_source.astrometric_n_good_obs_al, gaia_source.astrometric_excess_noise,gaia_source.phot_bp_rp_excess_factor  FROM gaiadr2.gaia_source WHERE INTERSECTS(CIRCLE('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec, 0.05),POLYGON('ICRS', %f, %f, %f, %f, %f, %f, %f, %f))=1;" % tpl)   
print("Query3: \n",query_str)
print()

