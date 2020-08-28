from astropy.io import fits as pyfits
import numpy as np
try:
    import matplotlib.pyplot as plt
except:
    pass    

def generate():
    pass

def prepare_classify(ifn, extension=1, ofn=None, overwrite=False, verbose=1):
    """
    Keep only relevant columns
    
    Parameters
    ----------
    
    """
    relevant_cols = ["srcID", "Fx","Fg","bp_rp","offset_sig","parallax","eligible_sky_density"]
    
    ff = pyfits.open(ifn)
    
    cols = []
    columns = {col.name:col.format for col in ff[extension].columns}
    
    col = pyfits.Column(name="srcID", array=ff[extension].data["srcID"], format=columns["srcID"])    
    cols.append(col)
    
    arr = np.log10(ff[extension].data["Fx"])+13
    print("log Fx",np.mean(arr), np.median(arr), np.std(arr))
    col = pyfits.Column(name="logFx", array=arr, format=columns["Fx"])    
    cols.append(col)
    
    arr = np.log10(ff[extension].data["Fg"])+12
    print("log Fg",np.mean(arr), np.median(arr), np.std(arr))
    col = pyfits.Column(name="logFg", array=arr, format=columns["Fg"])    
    cols.append(col)
    
    arr = np.log10(ff[extension].data["Fx"]) - np.log10(ff[extension].data["Fg"])+4
    print("log FxFg",np.mean(arr), np.median(arr), np.std(arr))   
    col = pyfits.Column(name="logFxFg", array=arr, format=columns["Fg"])    
    cols.append(col)
    
    
    arr = ff[extension].data["bp_rp"]
    print("bp_rp",np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
    col = pyfits.Column(name="bp_rp", array=arr, format=columns["bp_rp"])    
    cols.append(col)
    
    dst = ff[extension].data["match_dist"]
    dst = ff[extension].data["fake_match_dist"]
    err = ff[extension].data["RADEC_ERR"]
    val = np.zeros(len(dst))
    for i, (x, s) in enumerate(zip(dst, err)):
        xxx = np.arange(0,x,0.001)
        y = len(xxx)*[0.1]
        y = xxx/s**2 * np.exp(-xxx**2/(2*s**2))
        val[i] = 1- np.trapz(y, xxx)
    #print(val)
    #import matplotlib.pyplot as plt
    plt.hist(ff[extension].data["offset_sig"])
    plt.show()
    plt.scatter(ff[extension].data["offset_sig"], val)
    plt.show()
    print("pos",np.nanmean(val), np.nanmedian(val), np.nanstd(val))
    col = pyfits.Column(name="pos", array=val, format=columns["match_dist"])    
    cols.append(col)
    
    arr = np.log10(ff[extension].data["parallax"])
    gi = np.where(np.isnan(arr))[0]
    arr[gi] = -1
    print("plx",np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
    col = pyfits.Column(name="log_plx", array=arr, format=columns["parallax"])    
    cols.append(col)
    
    arr = np.log10(ff[extension].data["eligible_sky_density"])
    print("sky_density",np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
    col = pyfits.Column(name="log_sk", array=arr, format=columns["eligible_sky_density"])    
    cols.append(col)
        
    arr = ff[extension].data["category"]
    print("category",np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
    col = pyfits.Column(name="category", array=arr, format=columns["category"])    
    cols.append(col)
         
    hdu = pyfits.PrimaryHDU()    
    cc = pyfits.ColDefs(cols)
    xx = pyfits.BinTableHDU.from_columns(cc)
    hdul = pyfits.HDUList([hdu, xx])
    print("XXXXXX")
    if ofn is not None:
        hdul.writeto(ofn, overwrite=overwrite)        
        if verbose>0:
            print("datasets::prep_classify - Written ",len(dst)," objects with ",len(cols)," properties to ",ofn)
    ff.close()            
    return hdul
