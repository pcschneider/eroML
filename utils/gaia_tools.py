import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
from astroquery.gaia import Gaia
from astropy.io import fits as pyfits
from astropy.io.votable import parse
from astropy.io.votable import parse_single_table

def vo2fits(ifn, ofn, overwrite=False):
    """
    Convert VO table (file) into fits-file
    """
    hdu = pyfits.PrimaryHDU()       
    print("Reading: ",ifn)
    votable = parse_single_table(ifn).to_table()
    print("Rows:",len(votable), " Columns: ",len(votable.columns))
    cols = []
    for f in votable.columns:
        cc = votable[f]
        
        dtype = cc.dtype
        #print(f, dtype)#, cc._format, , type(f), type(cc), dir(cc))
        if dtype == "int16" or dtype=="int32" or dtype=="int64":
            fmt = "K"
        elif dtype == "float32":
            fmt = "E"
        elif dtype == "float64":   
            fmt = "D"
        elif f=="phot_variable_flag":     
            fmt = "10A"            
        elif dtype == "object":
            fmt = "L"
        else: fmt = None
        mm = cc.torecords()
        print("Adding column \'%s\' with dtype=%s using fits-format=%s." % (f,dtype,fmt))
        c = pyfits.Column(name=f, array=cc.data, format=fmt)
        cols.append(c)

    cc = pyfits.ColDefs(cols)
    xx = pyfits.BinTableHDU.from_columns(cc)
    hdul = pyfits.HDUList([hdu, xx])
    
    qc = quality_filter(hdul)
    c = pyfits.Column(name="Quality", array=qc, format="L")
    cols.append(c)
    cc = pyfits.ColDefs(cols)
    xx = pyfits.BinTableHDU.from_columns(cc)
    hdul = pyfits.HDUList([hdu, xx])
    
    hdul.writeto(ofn, overwrite=overwrite)
    print("Written: \'%s\'" % ofn)        


def quality_filter(ff, filter_Nr=0):
    """
    Determine Gaia-sources for which certain filter criteria are fullfilled
    
    Parameters
    ----------
    ff : pyfits.HDUList
    filter_Nr : int
        There are three filters defined (0, 1, 2). Filter-Nr 0 corresponds to the tightest contrains and equals the criteria 
        defined in <> to obtain a well-defined HR diagram.
    """
    d = ff[1].data
    N = len(d["source_id"])
    tmp = np.array([N*[1],np.exp(-0.4*(d["phot_g_mean_mag"]-19.5))])
    
    tmp = np.max(tmp.T, axis=1)
    #print(tmp)
    if filter_Nr==0:
        gi = np.where((d["parallax"]/d["parallax_error"] > 10) & (d["phot_g_mean_flux_over_error"]>50) &\
            (d["phot_rp_mean_flux_over_error"]>20) & (d["phot_bp_mean_flux_over_error"]>20) &\
            (d["phot_bp_rp_excess_factor"] < 1.3+0.06 *(d["phot_bp_mean_mag"]-d["phot_rp_mean_mag"])**2) &\
            (d["phot_bp_rp_excess_factor"] > 1.0+0.015*(d["phot_bp_mean_mag"]-d["phot_rp_mean_mag"])**2) &\
            (d["visibility_periods_used"]>8) &\
            (d["astrometric_chi2_al"]/(d["astrometric_n_good_obs_al"]-5)<1.44*tmp))[0]
    elif filter_Nr==1:
        gi = np.where((d["parallax"]/d["parallax_error"] > 3) & (d["phot_g_mean_flux_over_error"]>50) &\
            (d["phot_rp_mean_flux_over_error"]>20) & (d["phot_bp_mean_flux_over_error"]>20) &\
            (d["phot_bp_rp_excess_factor"] < 1.3+0.06 *(d["phot_bp_mean_mag"]-d["phot_rp_mean_mag"])**2) &\
            (d["phot_bp_rp_excess_factor"] > 1.0+0.015*(d["phot_bp_mean_mag"]-d["phot_rp_mean_mag"])**2) &\
            (d["visibility_periods_used"]>8) &\
            (d["astrometric_chi2_al"]/(d["astrometric_n_good_obs_al"]-5)<1.44*tmp))[0]
    elif filter_Nr==2:
        gi = np.where((d["parallax"]/d["parallax_error"] > 3) & (d["phot_g_mean_flux_over_error"]>30) &\
            (d["phot_rp_mean_flux_over_error"]>10) & (d["phot_bp_mean_flux_over_error"]>10) &\
            (d["phot_bp_rp_excess_factor"] < 1.3+0.06 *(d["phot_bp_mean_mag"]-d["phot_rp_mean_mag"])**2) &\
            (d["phot_bp_rp_excess_factor"] > 1.0+0.015*(d["phot_bp_mean_mag"]-d["phot_rp_mean_mag"])**2) &\
            (d["visibility_periods_used"]>8) &\
            (d["astrometric_chi2_al"]/(d["astrometric_n_good_obs_al"]-5)<1.44*tmp))[0]
    q = np.zeros(N)
    q[gi] = 1
    return q

def add_quality_column(fn, ofile, colname="quality", overwrite=False, filterNr=0):
    ff = pyfits.open(fn)
    cols = ff[1].columns
    print(cols) 
    q = quality_filter(ff, filterNr=filterNr)
    c = pyfits.Column(name=colname, array=q, format="L")
    print("#rows: %i, good quality: %i, fraction: %f" % (len(ff[1].data[cols[0].name]), np.sum(q),np.sum(q)/len(ff[1].data[cols[0].name]) ))
    cols.add_col(c)
    hdu = copy.copy(ff[0]) 
    hdx = pyfits.BinTableHDU.from_columns(cols)
    hdul = pyfits.HDUList([hdu, hdx])
    hdul.writeto(ofile, overwrite=overwrite)    #hdul = pyfits.HDUList([hdu, ff[1])
    


def gaia4ero(ifn, ofn=None, ext=1, radec_cols=("RA", "DEC"), verbose=10, keep_VO=False, overwrite=False):
    """
    Download Gaia sources for the sky region covered by an eROSITA source catalog
    
    The sky region downloaded is a "square" in coordinates-VALUES, not on the sky.
    
    Parameters
    ----------
    ifn : str
        Gaia source catalog filename
    ext : int
        Extension in fits-file
    radec_cols : tuple
        Column names for RA and Dec entries
    keep_VO : boolean
        The temporary VO-file is kept if True
    """
    ff = pyfits.open(ifn)
    ra, dec = ff[ext].data[radec_cols[0]], ff[ext].data[radec_cols[1]]
   
    width = max(ra) - min(ra)
    height = max(dec) - min(dec)
    RA_center = min(ra)+width/2
    Dec_center = min(dec)+height/2
    
    width/=100
    height/=100
    
    if verbose>0:
        print("Downloading Gaia sources around RA, Dec = (",RA_center, Dec_center,") with width, height = (", width, height,")")

    #exit()
    ww= u.Quantity(width, u.deg)
    hh = u.Quantity(height, u.deg)

    coord = SkyCoord(ra=RA_center, dec=Dec_center, unit=(u.degree, u.degree), frame='icrs')
    #r = Gaia.query_object_async(coordinate=coord, width=width, height=height, dump_to_file=True)

    qc=False
        
    if qc:
        query_str = str("SELECT gaia_source.source_id,gaia_source.ra,gaia_source.ra_error,gaia_source.dec,gaia_source.dec_error,gaia_source.parallax,gaia_source.parallax_error,gaia_source.phot_g_mean_mag,gaia_source.bp_rp,gaia_source.radial_velocity,gaia_source.radial_velocity_error,gaia_source.phot_bp_mean_mag,gaia_source.phot_rp_mean_mag,gaia_source.phot_g_mean_flux_over_error,gaia_source.phot_rp_mean_flux_over_error,gaia_source.phot_bp_mean_flux_over_error,gaia_source.phot_variable_flag,gaia_source.teff_val,gaia_source.a_g_val,gaia_source.visibility_periods_used,gaia_source.astrometric_chi2_al,gaia_source.astrometric_n_good_obs_al, gaia_source.astrometric_excess_noise,gaia_source.phot_bp_rp_excess_factor  FROM gaiadr2.gaia_source WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),BOX('ICRS', %f, %f, %f, %f))=1 AND parallax_over_error > 10 AND phot_g_mean_flux_over_error>50 AND phot_rp_mean_flux_over_error>20 AND phot_bp_mean_flux_over_error>20 AND phot_bp_rp_excess_factor < 1.3+0.06*power(phot_bp_mean_mag-phot_rp_mean_mag,2) AND phot_bp_rp_excess_factor > 1.0+0.015*power(phot_bp_mean_mag-phot_rp_mean_mag,2) AND visibility_periods_used>8 AND astrometric_chi2_al/(astrometric_n_good_obs_al-5)<1.44*greatest(1,exp(-0.4*(phot_g_mean_mag-19.5)));" % (RA_center, Dec_center, width*1.01,height*1.01))
    else:
        query_str = str("SELECT gaia_source.source_id,gaia_source.ra,gaia_source.ra_error,gaia_source.dec,gaia_source.dec_error,gaia_source.parallax,gaia_source.parallax_error,gaia_source.phot_g_mean_mag,gaia_source.bp_rp,gaia_source.radial_velocity,gaia_source.radial_velocity_error,gaia_source.phot_bp_mean_mag,gaia_source.phot_rp_mean_mag,gaia_source.phot_g_mean_flux_over_error,gaia_source.phot_rp_mean_flux_over_error,gaia_source.phot_bp_mean_flux_over_error,gaia_source.phot_variable_flag,gaia_source.teff_val,gaia_source.a_g_val,gaia_source.visibility_periods_used,gaia_source.astrometric_chi2_al,gaia_source.astrometric_n_good_obs_al, gaia_source.astrometric_excess_noise,gaia_source.phot_bp_rp_excess_factor  FROM gaiadr2.gaia_source WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),BOX('ICRS', %f, %f, %f, %f))=1;" % (RA_center, Dec_center, width*1.01,height*1.01))   
   
    if verbose>3:
        print(query_str)
    job = Gaia.launch_job_async(query_str, dump_to_file=True)
    #job = Gaia.launch_job_async("SELECT  gaia_source.source_id,gaia_source.ra,gaia_source.ra_error,gaia_source.dec,gaia_source.dec_error,gaia_source.parallax,gaia_source.parallax_error,gaia_source.phot_g_mean_mag,gaia_source.bp_rp,gaia_source.radial_velocity,gaia_source.radial_velocity_error,gaia_source.phot_variable_flag,gaia_source.teff_val,gaia_source.a_g_val FROM gaiadr2.gaia_source WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),BOX('ICRS',136.03142581392405,1.5018720027676264,18.822564537570514,8.516702022938103))=1;", dump_to_file=True)

    if verbose>3: print(job)
    r = job.get_results()
    
    tmp_fn = job.outputFile
    if ofn is not None:
        vo2fits(tmp_fn, ofn, overwrite=overwrite)
        
    if keep_VO == False:
        import os
        os.remove(tmp_fn)
