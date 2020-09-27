import astropy.units as u
from astropy.coordinates import SkyCoord
from eroML.ensemble import fits_support
import numpy as np
from astroquery.gaia import Gaia
from astropy.io import fits as pyfits
from astropy.io.votable import parse
from astropy.io.votable import parse_single_table
import tempfile
import copy
import os, re
import healpy as hp
import logging
import glob
#from .enrich import enrich_Gaia

logger = logging.getLogger('eroML')

def get_alternate_gaia_file(glob_str):
    """
    Return 
    """
    print("CHECK checking ",glob_str)
    #glob_str = pre+"*"+post
    fnames = glob.glob(glob_str)
    print("fnames: ",fnames)
    if len(fnames)==0:
        return None
    else:
        good = None
        for fn in fnames:
            try: 
                ff = pyfits.open(fn)
                good = fn
            except:
                continue
        print(" good: ",good)    
        return good    

def download_Gaia_tiles(outdir=".", prefix="Gaia", idx=None, nside=None, overwrite=False, verbose=1, edge=3., keep_VO=False, check_alternate=True, Glim=23):
    """
    Download all Gaia sources
    """

    if idx is None:
        idx = range(hp.nside2npix(nside))
    logger.info("Number of hpix: %i" % len(idx))   
    for j, i in enumerate(idx):
        
        pre  = outdir+"/"+prefix+"_nside"+str(nside)+"_"+str(i)
        post = ".fits"
        ofn = pre+post
                
        logger.debug("Using \'%s\' for hpix=%i  (pix %i of %i)" % (ofn, i, j+1, len(idx)))
        if os.path.exists(ofn):
            logger.debug("  Skipping...")
            continue
        
        if check_alternate:
            glob_str = outdir+"/Gaia_*"+"_nside"+str(nside)+"_"+str(i)+".fits"
            afn = get_alternate_gaia_file(glob_str)
            if afn: 
                os.symlink(afn, ofn)
                continue
        
        download_one_Gaia_polytile(ofn, i, nside, overwrite=overwrite, verbose=verbose, edge=edge, keep_VO=keep_VO,)
        add_standard_cols(ofn, overwrite=True)
        #add_quality_column(ofn, ofn, overwrite=True)



def Gaia_tile_loop(idx, prefix=None, postfix=None, filterNr=3):
    """
    Loop through ero tiles and enrich them
    """
    from .enrich import enrich_Gaia
    logger.info("Enrichting %i Gaia source tiles." % len(idx))
    
    p = re.compile("(\S*)/Gaia_(\S*)_nside(\d+)_(\S*).fits")

    for j, i in enumerate(idx):
        fn = prefix+str(i)+postfix+'.fits'
        logger.debug("Enriching Gaia tile: %s. (file %i/%i; filterNr=%i) " % (fn, j+1, len(idx), filterNr)) 

        if not os.path.exists(fn):
            
            mm = p.match(fn)
            #print(mm, mm.group(0), mm.group(2))
            if mm:
                found = mm.group(3)
                nside = int(found)
            else:
                logger.warning("Cannot find %s or alternative." % fn)
                raise FileNotFoundError("Cannot find '%s' or alternative." % fn)

            glob_str =  mm.group(1)+"/Gaia_*"+"_nside"+str(nside)+"_"+str(i)+".fits"
            #print("GLOB:", glob_str)
            afn = get_alternate_gaia_file(glob_str)    
            #print("AFN", afn)
            if afn:
                logger.info("Using %s instead of %s." % (afn, fn))
                os.symlink(afn, fn)
                
            
        enrich_Gaia(fn, filterNr=filterNr)    


def get_larger_poly(x, y, around=3):
    """
    Parameters
    ----------
    x, y : np.array of size 4 (float), 
        RA, Dec of the polynom
    around : float
        increase size by this amount (in arcmin)
    """
    
    minxi, maxxi = np.argmin(x), np.argmax(x)
    minyi, maxyi = np.argmin(y), np.argmax(y)
    
    co = SkyCoord(ra=x[minxi], dec=y[minxi], unit=(u.degree, u.degree))
    co1 = co.directional_offset_by(270*u.degree, around*u.arcmin)
    x[minxi] = co1.ra.degree
    y[minxi] = co1.dec.degree
    
    co = SkyCoord(ra=x[maxxi], dec=y[maxxi], unit=(u.degree, u.degree))
    co1 = co.directional_offset_by(90*u.degree, around*u.arcmin)
    x[maxxi] = co1.ra.degree
    y[maxxi] = co1.dec.degree
    
    co = SkyCoord(ra=x[minyi], dec=y[minyi], unit=(u.degree, u.degree))
    co1 = co.directional_offset_by(180*u.degree, around*u.arcmin)
    x[minyi] = co1.ra.degree
    y[minyi] = co1.dec.degree
    
    co = SkyCoord(ra=x[maxyi], dec=y[maxyi], unit=(u.degree, u.degree))
    co1 = co.directional_offset_by(0*u.degree, around*u.arcmin)
    x[maxyi] = co1.ra.degree
    y[maxyi] = co1.dec.degree
    
    return x, y


def download_one_Gaia_polytile(ofn, hpix, nside, overwrite=False, verbose=1, edge=3., keep_VO=False, Glim=23):
    """
    Download Gaia sources for specific healpix index
    
    Parameters
    -----------
    """
    step=1
    resol = hp.nside2resol(nside, arcmin=True)
    pix_center = hp.pix2ang(nside, hpix, nest=True)
    ra, dec = pix_center[1]/np.pi*180, 90. - pix_center[0]/np.pi*180
    logger.log(10, "pos: "+str(pix_center)+ " resol: "+str(resol)+" RA, Dec: "+str(ra)+" "+str(dec))
    logger.log(10, "pix_center: %f, %f " % (90-pix_center[0]*180/np.pi, pix_center[1]*180/np.pi))
    
    #ww= u.Quantity(resol, u.arcmin) + 2*edge*u.arcmin
    #hh = u.Quantity(resol, u.arcmin) + 2*edge*u.arcmin

    #width = ww.to(u.degree).value
    #height = hh.to(u.degree).value
    
    
    b = hp.boundaries(nside, hpix, nest=True, step=step).transpose()
    v = hp.pixelfunc.vec2ang(b)
    x0, y0 = v[1]/np.pi*180, 90.-v[0]/np.pi*180  

    x, y = get_larger_poly(x0, y0, around=3*edge)
     

    tpl = (x[0], y[0], x[1], y[1], x[2], y[2], x[3], y[3])

    query_str = str("SELECT gaia_source.source_id,gaia_source.ra,gaia_source.pmra, gaia_source.pmdec, gaia_source.ref_epoch, gaia_source.ra_error,gaia_source.dec,gaia_source.dec_error,gaia_source.parallax,gaia_source.parallax_error,gaia_source.phot_g_mean_mag,gaia_source.bp_rp,gaia_source.radial_velocity,gaia_source.radial_velocity_error,gaia_source.phot_bp_mean_mag,gaia_source.phot_rp_mean_mag,gaia_source.phot_g_mean_flux_over_error,gaia_source.phot_rp_mean_flux_over_error,gaia_source.phot_bp_mean_flux_over_error,gaia_source.phot_variable_flag,gaia_source.teff_val,gaia_source.a_g_val,gaia_source.visibility_periods_used,gaia_source.astrometric_chi2_al,gaia_source.astrometric_n_good_obs_al, gaia_source.astrometric_excess_noise,gaia_source.phot_bp_rp_excess_factor  FROM gaiadr2.gaia_source WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),POLYGON('ICRS', %f, %f, %f, %f, %f, %f, %f, %f))=1;" % tpl) 
       
    logger.log(5, "query_str: "+query_str)
    #return    
    job = Gaia.launch_job_async(query_str, dump_to_file=True)

    if verbose>3: logger.log(0, " (job)"+str(job))
    r = job.get_results()
    
    tmp_fn = job.outputFile
    if ofn is not None:
        vo2fits(tmp_fn, ofn, overwrite=overwrite)
        
    if keep_VO == False:
        os.remove(tmp_fn)
             
         
def download_one_Gaia_square_tile(ofn, hpix, nside, overwrite=False, verbose=1, edge=3., keep_VO=False, Glim=23):
    """
    Download Gaia sources within square centered on specific healpix index
    
    Parameters
    -----------
    """
    resol = hp.nside2resol(nside, arcmin=True)
    pix_center = hp.pix2ang(nside, hpix, nest=True)
    ra, dec = pix_center[1]/np.pi*180, 90. - pix_center[0]/np.pi*180
    logger.log(10, "pos: "+str(pix_center)+ " resol: "+str(resol)+" RA, Dec: "+str(ra)+" "+str(dec))
    logger.log(10, "pix_center: %f, %f " % (90-pix_center[0]*180/np.pi, pix_center[1]*180/np.pi))
    
    ww= u.Quantity(resol, u.arcmin) + 2*edge*u.arcmin
    hh = u.Quantity(resol, u.arcmin) + 2*edge*u.arcmin

    width = ww.to(u.degree).value
    height = hh.to(u.degree).value
    
    query_str = str("SELECT gaia_source.source_id,gaia_source.ra,gaia_source.ra_error,gaia_source.dec,gaia_source.pmra, gaia_source.pmdec, gaia_source.ref_epoch,gaia_source.dec_error,gaia_source.parallax,gaia_source.parallax_error,gaia_source.phot_g_mean_mag,gaia_source.bp_rp,gaia_source.radial_velocity,gaia_source.radial_velocity_error,gaia_source.phot_bp_mean_mag,gaia_source.phot_rp_mean_mag,gaia_source.phot_g_mean_flux_over_error,gaia_source.phot_rp_mean_flux_over_error,gaia_source.phot_bp_mean_flux_over_error,gaia_source.phot_variable_flag,gaia_source.teff_val,gaia_source.a_g_val,gaia_source.visibility_periods_used,gaia_source.astrometric_chi2_al,gaia_source.astrometric_n_good_obs_al, gaia_source.astrometric_excess_noise,gaia_source.phot_bp_rp_excess_factor  FROM gaiadr2.gaia_source WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),BOX('ICRS', %f, %f, %f, %f))=1;" % (ra, dec, width,height))   
       
    logger.log(5, "query_str: "+query_str)
    #return    
    job = Gaia.launch_job_async(query_str, dump_to_file=True)

    if verbose>3: logger.log(0, " (job)"+str(job))
    r = job.get_results()
    
    tmp_fn = job.outputFile
    if ofn is not None:
        vo2fits(tmp_fn, ofn, overwrite=overwrite)
        
    if keep_VO == False:
        os.remove(tmp_fn)
        
    
    
    
def vo2fits(ifn, ofn, verbose=1, overwrite=False):
    """
    Convert VO table (file) into fits-file
    """
    hdu = pyfits.PrimaryHDU()       
    logger.debug("Reading: %s" % ifn)
    votable = parse_single_table(ifn).to_table()
    logger.debug("   Rows: %i Columns: %i " % (len(votable), len(votable.columns)))
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
            fmt = "25A"            
        elif dtype == "object":
            fmt = "L"
        else: fmt = None
        mm = cc.torecords()
        if verbose>1: print("gaia_tools::vo2fits - Adding column \'%s\' with dtype=%s using fits-format=%s." % (f,dtype,fmt))
        if f=="ra":
            c = pyfits.Column(name="RA", array=cc.data, format=fmt)
        elif f=="dec":
            c = pyfits.Column(name="Dec", array=cc.data, format=fmt)
        elif f=="pmra":
            c = pyfits.Column(name="pm_RA", array=cc.data, format=fmt)
        elif f=="pmdec":
            c = pyfits.Column(name="pm_Dec", array=cc.data, format=fmt)            
        else:
            c = pyfits.Column(name=f, array=cc.data, format=fmt)
        cols.append(c)

    cc = pyfits.ColDefs(cols)
    xx = pyfits.BinTableHDU.from_columns(cc)
    hdul = pyfits.HDUList([hdu, xx])
    
    #qc = quality_filter(hdul)
    #c = pyfits.Column(name="Gaia_Quality", array=qc, format="I")
    #cols.append(c)
    cc = pyfits.ColDefs(cols)
    xx = pyfits.BinTableHDU.from_columns(cc)
    hdul = pyfits.HDUList([hdu, xx])
    
    hdul.writeto(ofn, overwrite=overwrite)
    logger.debug("   Written: \'%s\'" % ofn)        


    

def quality_filter(ff, filter_Nr=3):
    """
    Determine Gaia-sources for which certain filter criteria are fullfilled
    
    Parameters
    ----------
    ff : pyfits.HDUList
    filter_Nr : int
        There are three filters defined (0, 1, 2, 3). Filter-Nr 0 corresponds to the tightest contrains and equals the criteria 
        defined in <> to obtain a well-defined HR diagram.
    """
    #filter_Nr = 3
    print("using Filter: ",filter_Nr)
    
    try:
        ff["parallax"]
        is_array=True
    except:
        is_array=False
        
        
    if not is_array:
        d = ff[1].data
    else:
        d = ff
        
    N = len(d["srcID"])
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
    elif filter_Nr==3:
       gi = np.where((d["parallax"]/d["parallax_error"] > 3) & (d["phot_g_mean_flux_over_error"]>30) &\
            (d["phot_rp_mean_flux_over_error"]>10) & (d["phot_bp_mean_flux_over_error"]>10) &\
            (d["phot_bp_rp_excess_factor"] < 1.3+0.06 *(d["phot_bp_mean_mag"]-d["phot_rp_mean_mag"])**2) &\
            (d["phot_bp_rp_excess_factor"] > 1.0+0.015*(d["phot_bp_mean_mag"]-d["phot_rp_mean_mag"])**2) &\
            (d["visibility_periods_used"]>4) &\
            (d["astrometric_chi2_al"]/(d["astrometric_n_good_obs_al"]-5)<1.44*tmp))[0]

        
    q = np.zeros(N)
    q[gi] = 1
    logger.debug("Number of sources fullfilling Gaia quality criterium: %i (%f)" % (np.sum(q), np.sum(q)/N))
    return q

@fits_support
def add_quality_column(e, colname="Gaia_quality", filter_Nr=2):
    q = quality_filter(e.array, filter_Nr=filter_Nr).astype(int)
    e.add_col(colname, q)
    return e
    #c = pyfits.Column(name=colname, array=q, format="L")
    #if verbose>1: print("gaia_tools::add_quality_column - #rows: %i, good quality: %i, fraction: %f" % (len(ff[1].data[cols[0].name]), np.sum(q),np.sum(q)/len(ff[1].data[cols[0].name]) ))
    #cols.add_col(c)
    #hdu = copy.copy(ff[0]) 
    ##hdu = ff[0]
    ##hdu = pyfits.PrimaryHDU()
    #hdx = pyfits.BinTableHDU.from_columns(cols)
    #hdul = pyfits.HDUList([hdu, hdx])
    #hdul.writeto(ofile, overwrite=overwrite)
    #ff.close()


    
def add_standard_cols(ifn, overwrite=True):
    """
    """
    ff = pyfits.open(ifn)
    cols = ff[1].columns
    srcIDs = ff[1].data["source_id"].astype(str)
    c = pyfits.Column(name="srcID", array=srcIDs, format="38A")
    
    
    #print("#rows: %i, good quality: %i, fraction: %f" % (len(ff[1].data[cols[0].name]), np.sum(q),np.sum(q)/len(ff[1].data[cols[0].name]) ))
    cols.add_col(c)
    hdu = copy.copy(ff[0]) 
    #hdu = ff[0]
    hdx = pyfits.BinTableHDU.from_columns(cols)
    hdul = pyfits.HDUList([hdu, hdx])
    hdul.writeto(ifn, overwrite=overwrite)    #hdul = pyfits.HDUList([hdu, ff[1])
    ff.close()
    
 
def get_gaia(ifn, ofn, overwrite=False, verbose=1):
    """
    Download Gaia sources and enrich file, keep only relevant columns.
    """    
    if os.path.exists(ofn):
        if verbose>0: print("gaia_tools::get_gaia - ",ofn," already exists. Assuming it is the correct file.")
        return
    fh, tmp_fn = tempfile.mkstemp(dir='.', suffix='.fits')
    gaia4ero(ifn, ofn=tmp_fn, verbose=verbose)
    add_standard_cols(tmp_fn, overwrite=True)
    add_quality_column(tmp_fn, ofn, overwrite=True)
    os.close(fh)
    os.remove(tmp_fn)


def gaia4ero(ifn, ofn=None, ext=1, radec_cols=("RA", "DEC"), verbose=1, keep_VO=False, overwrite=False):
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
    if verbose>3: print("gaia_tools::gaia4ero - Reading: ",ifn)
    ff = pyfits.open(ifn)
    ra, dec = ff[ext].data[radec_cols[0]], ff[ext].data[radec_cols[1]]
   
    width = max(ra) - min(ra)
    height = max(dec) - min(dec)
    RA_center = min(ra)+width/2
    Dec_center = min(dec)+height/2
    
    #width/=100
    #height/=100
    
    if verbose>0:
        print("gaia_tools::gaia4ero - Downloading Gaia sources around RA, Dec = (",RA_center, Dec_center,") with width, height = (", width, height,")")

    #exit()
    ww= u.Quantity(width, u.deg) + 10.*u.arcmin
    hh = u.Quantity(height, u.deg) + 10.*u.arcmin

    width = ww.to(u.degree).value
    height = hh.to(u.degree).value

    #width/=30
    #height/=30

    coord = SkyCoord(ra=RA_center, dec=Dec_center, unit=(u.degree, u.degree), frame='icrs')
    #r = Gaia.query_object_async(coordinate=coord, width=width, height=height, dump_to_file=True)

    qc=False
        
    if qc:
        query_str = str("SELECT gaia_source.source_id,gaia_source.ra,gaia_source.ra_error,gaia_source.dec,gaia_source.dec_error,gaia_source.parallax,gaia_source.parallax_error,gaia_source.phot_g_mean_mag,gaia_source.bp_rp,gaia_source.radial_velocity,gaia_source.radial_velocity_error,gaia_source.phot_bp_mean_mag,gaia_source.phot_rp_mean_mag,gaia_source.phot_g_mean_flux_over_error,gaia_source.phot_rp_mean_flux_over_error,gaia_source.phot_bp_mean_flux_over_error,gaia_source.phot_variable_flag,gaia_source.teff_val,gaia_source.a_g_val,gaia_source.visibility_periods_used,gaia_source.astrometric_chi2_al,gaia_source.astrometric_n_good_obs_al, gaia_source.astrometric_excess_noise,gaia_source.phot_bp_rp_excess_factor  FROM gaiadr2.gaia_source WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),BOX('ICRS', %f, %f, %f, %f))=1 AND parallax_over_error > 10 AND phot_g_mean_flux_over_error>50 AND phot_rp_mean_flux_over_error>20 AND phot_bp_mean_flux_over_error>20 AND phot_bp_rp_excess_factor < 1.3+0.06*power(phot_bp_mean_mag-phot_rp_mean_mag,2) AND phot_bp_rp_excess_factor > 1.0+0.015*power(phot_bp_mean_mag-phot_rp_mean_mag,2) AND visibility_periods_used>8 AND astrometric_chi2_al/(astrometric_n_good_obs_al-5)<1.44*greatest(1,exp(-0.4*(phot_g_mean_mag-19.5)));" % (RA_center, Dec_center, width,height))
    else:
        query_str = str("SELECT gaia_source.source_id,gaia_source.ra,gaia_source.ra_error,gaia_source.dec,gaia_source.dec_error,gaia_source.parallax,gaia_source.parallax_error,gaia_source.phot_g_mean_mag,gaia_source.bp_rp,gaia_source.radial_velocity,gaia_source.radial_velocity_error,gaia_source.phot_bp_mean_mag,gaia_source.phot_rp_mean_mag,gaia_source.phot_g_mean_flux_over_error,gaia_source.phot_rp_mean_flux_over_error,gaia_source.phot_bp_mean_flux_over_error,gaia_source.phot_variable_flag,gaia_source.teff_val,gaia_source.a_g_val,gaia_source.visibility_periods_used,gaia_source.astrometric_chi2_al,gaia_source.astrometric_n_good_obs_al, gaia_source.astrometric_excess_noise,gaia_source.phot_bp_rp_excess_factor  FROM gaiadr2.gaia_source WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),BOX('ICRS', %f, %f, %f, %f))=1;" % (RA_center, Dec_center, width,height))   
   
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
        os.remove(tmp_fn)
        
    ff.close()
