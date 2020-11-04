from astropy.io import fits as pyfits
import numpy as np
try:
    import matplotlib.pyplot as plt
except:
    pass    

def generate():
    pass

def gen_real_pos_offset(N, sigma=1.):
        
    rnd = np.random.rand(N)
    rndx = np.sqrt(2*sigma**2 * (-np.log(1-rnd)))

    return rndx

def gen_random_pos_offset(N, dens=1.):

    rnd = np.random.rand(N)
    rndx = np.sqrt(-np.log(1-rnd)/(np.pi*dens))
    return rndx


def repeat_fits(hdu, multi=10):
    """
    """
    NN0 = len(hdu.data["srcID"])
    cols = []
    for c in hdu.columns:
        #print(c)
        col = pyfits.Column(name=c.name, array=np.repeat(hdu.data[c.name], multi), format=c.format)    
        cols.append(col)
    cc = pyfits.ColDefs(cols)    
    xx = pyfits.BinTableHDU.from_columns(cc)
    xx.header = hdu.header
    #print("repeat: ",NN0, len(xx.data["srcID"]))
    return xx
    

def prepare_classify(ifn, extension=1, ofn=None, overwrite=False, verbose=1, display=False):
    """
    Keep only relevant columns
    
    Parameters
    ----------
    
    """
    relevant_cols = ["srcID", "Fx","Fg","bp_rp","offset_sig","parallax","eligible_sky_density"]
    
    ff = pyfits.open(ifn)
    
    if "category" in ff[extension].data.columns.names:
        ff[extension] = repeat_fits(ff[extension], multi=10)
    
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
    err = ff[extension].data["RADEC_ERR"]
    
    if "category" in ff[extension].data.columns.names:
        dist = 1000/ff[extension].data["parallax"] 
        Lx = 4*np.pi*(3.1e18*dist)**2  * ff[extension].data["Fx"]      
        ci = np.where(ff[extension].data["category"] == 0)[0]
        gi = np.where(Lx[ci] > 2e31)[0]
        print("all: ",len(ci)," Lx filter: ",len(gi), "(max Lx: ",max(Lx[ci]),")")
        ff[extension].data["category"][ci[gi]] = 1
    
    
    dens = ff[extension].data['eligible_sky_density'] / 3600
    
    
    
    if "category" in ff[extension].data.columns.names:
        # random matches
        gi = np.where(ff[extension].data["category"] == 2)[0] 
        dst[gi] = gen_random_pos_offset(len(gi), dens = dens[gi])
        
        # real matches
        gi = np.where(ff[extension].data["category"] < 2)[0] 
        dst[gi] = gen_real_pos_offset(len(gi), sigma = err[gi])
    else:
        print("No category column.")
    
    
    #val = np.zeros(len(dst))
    #for i, (x, s) in enumerate(zip(dst, err)):
        #xxx = np.arange(0,x,0.001)
        #y = len(xxx)*[0.1]
        #y = xxx/s**2 * np.exp(-xxx**2/(2*s**2))
        #val[i] = 1- np.trapz(y, xxx)
    val = 1 - np.exp(-dst**2/(2*err**2))        
    print(val)
    #import matplotlib.pyplot as plt
    if "category" in ff[extension].data.columns.names:
        display=True
    if display:
        plt.hist(ff[extension].data["offset_sig"])
        plt.title(ifn)
        plt.xlabel("Offset sig")
        plt.show()
        plt.scatter(ff[extension].data["offset_sig"], (val*6)**2)
        plt.xlabel("offset sig")
        plt.ylabel("pos")
        plt.title(ifn)
        plt.show()
    print("pos",np.nanmean(val), np.nanmedian(val), np.nanstd(val))
    col = pyfits.Column(name="pos", array=(val*6)**2, format=columns["match_dist"])
    col = pyfits.Column(name="pos", array=(val*6)**2, format=columns["match_dist"])    
    cols.append(col)
    
    
    
    
    arr = np.log10(ff[extension].data["parallax"])
    gi = np.where(np.isnan(arr))[0]
    arr[gi] = -1
    print("plx",np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
    col = pyfits.Column(name="log_plx", array=arr*5, format=columns["parallax"])    
    cols.append(col)
    
    arr = np.log10(ff[extension].data["eligible_sky_density"])
    print("sky_density",np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
    col = pyfits.Column(name="log_sk", array=arr, format=columns["eligible_sky_density"])    
    cols.append(col)

    arr = ff[extension].data["eligible_sky_density"]
    print("sky_density",np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
    col = pyfits.Column(name="skd", array=arr, format=columns["eligible_sky_density"])    
    cols.append(col)
        
    if "category" in ff[extension].data.columns.names:
        sigma = np.sqrt(-2*np.log(1-val))
        print(sigma)
        col = pyfits.Column(name="offset_sig", array=sigma, format=columns["offset_sig"])    
        cols.append(col)
        plt.scatter(np.exp(val), sigma)
        plt.xlabel("pos")
        plt.ylabel("offset_sig")
        plt.title(ifn)
        plt.show()
    else:
        arr = ff[extension].data["offset_sig"]
        print("offset_sig",np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
        col = pyfits.Column(name="offset_sig", array=arr, format=columns["offset_sig"])    
        cols.append(col)
        
    
    
    arr = ff[extension].data["NN"]
    print("NN",np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
    col = pyfits.Column(name="NN", array=arr, format=columns["NN"])    
    cols.append(col)        
        
    arr = ff[extension].data["offset_sig"]
    print("offset_sig",np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
    col = pyfits.Column(name="offset_sig", array=arr, format=columns["offset_sig"])    
    cols.append(col)
    
    arr = ff[extension].data["NN"]
    print("NN",np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
    col = pyfits.Column(name="NN", array=arr, format=columns["NN"])    
    cols.append(col)        
        
    if "category" in ff[extension].data.columns.names:
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
