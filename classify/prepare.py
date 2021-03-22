from astropy.io import fits as pyfits
import numpy as np
from eroML.positions import gen_random_pos_offset, gen_real_pos_offset
from eroML.positions import calc_sigma_from_RADEC_ERR
from eroML.classify import training_filter

try:
    import matplotlib.pyplot as plt
except:
    pass    

def generate():
    pass

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
    return xx
    

def preprocess(ifn, extension=1, ofn=None, overwrite=False, verbose=1, display=False):
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
    col = pyfits.Column(name="logFx", array=5*arr, format=columns["Fx"])    
    cols.append(col)
    
    arr = np.log10(ff[extension].data["Fg"])+12
    print("log Fg",np.mean(arr), np.median(arr), np.std(arr))
    col = pyfits.Column(name="logFg", array=5*arr, format=columns["Fg"])    
    cols.append(col)
    
    arr = np.log10(ff[extension].data["Fx"]) - np.log10(ff[extension].data["Fg"])+4
    print("log FxFg",np.mean(arr), np.median(arr), np.std(arr))   
    col = pyfits.Column(name="logFxFg", array=5*arr, format=columns["Fg"])    
    cols.append(col)
    
    
    arr = ff[extension].data["bp_rp"]
    print("bp_rp",np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
    col = pyfits.Column(name="bp_rp", array=arr, format=columns["bp_rp"])    
    cols.append(col)
    
    dst = ff[extension].data["match_dist"]
    #dst = ff[extension].data["fake_match_dist"]
    err = ff[extension].data["RADEC_ERR"]
    
    
    if "category" in ff[extension].data.columns.names:
        #dist = 1000/ff[extension].data["parallax"] 
        #Lx = 4*np.pi*(3.1e18*dist)**2  * ff[extension].data["Fx"]      
        #ci = np.where(ff[extension].data["category"] == 0)[0]
        #gi = np.where(Lx[ci] > 2e31)[0]
        #print("all: ",len(ci)," Lx filter: ",len(gi), "(max Lx: ",max(Lx[ci]),")")
        #good = training_filter(
        dct = {}
        for kw in ["parallax","Fx", "FxFg", "bp_rp"]:
            dct[kw] = ff[extension].data[kw]
        ai = training_filter(dct)    
        gi = np.where(ai==1)[0]
        ff[extension].data["category"][gi] = 1
    
    
    dens = ff[extension].data['eligible_sky_density'] / 3600
    
    if "category" in ff[extension].data.columns.names:
        # random matches
        gi = np.where(ff[extension].data["category"] == 2)[0] 
        NN = ff[extension].data["NN"]
        sgm = 1/np.sqrt(dens[gi]*2*np.pi)
        dst[gi] = gen_real_pos_offset(len(gi), sigma = sgm)
        #dst[gi] = gen_random_pos_offset(dens = dens[gi], NN=NN[gi])
        
        
        bns = np.histogram(dst[gi], range=(0, 100), bins=30, density=True)

        bnx = (bns[1][1:] + bns[1][0:-1])/2

        bny = bns[0]

        plt.hist(dst[gi], bins=30, range=(0, 100), density=True)
        plt.bar(bnx, height=bny, color='r')
        
        
        M = len(bnx)
        N = len(gi)
        xx = bnx.repeat(N).reshape((M,N))
        sigma = sgm
        zz = xx/sigma**2*np.exp(-xx**2/(2*sigma**2))        
        mm = np.mean(zz, axis=1)
        plt.plot(bnx, mm)
        plt.title("Random")
        plt.xlabel("Radius (arcsec)")
        plt.show()
        
        # real matches
        gi = np.where(ff[extension].data["category"] < 2)[0] 
        dst[gi] = gen_real_pos_offset(len(gi), sigma = err[gi])
        
        plt.hist(dst[gi], bins=30, range=(0, 30), density=True)
        
        bnx = np.linspace(0,30, 1000)
        M = len(bnx)
        N = len(gi)
        xx = bnx.repeat(N).reshape((M,N))
        sigma = err[gi]
        zz = xx/sigma**2*np.exp(-xx**2/(2*sigma**2))        
        mm = np.mean(zz, axis=1)
        plt.plot(bnx, mm)

        plt.title("Real")
        plt.xlabel("Radius (arcsec)")
        plt.show()
        
        print("dst", dst)
        
        col = pyfits.Column(name="match_dist", array=dst, format=columns["match_dist"])    
        cols.append(col)
        col = pyfits.Column(name="sigma_r", array=ff[extension].data["RADEC_ERR"], format=columns["RADEC_ERR"])    
        cols.append(col)
        col = pyfits.Column(name="offset_sig", array=dst/err, format=columns["offset_sig"])    
        cols.append(col)
            
        arr = ff[extension].data["eligible_sky_density"] * dst**2 * 4*np.pi/3600
        print("expected_rnd",np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
        col = pyfits.Column(name="expected_rnd", array=arr*10, format=columns["match_dist"])    
        cols.append(col)        
            

        
    else:
        print("No category column.")
        col = pyfits.Column(name="match_dist", array=ff[extension].data["match_dist"], format=columns["match_dist"])    
        cols.append(col)
        col = pyfits.Column(name="sigma_r", array=0.61*ff[extension].data["RADEC_ERR"], format=columns["RADEC_ERR"])    
        cols.append(col)
        arr = ff[extension].data["offset_sig"]
        print("offset_sig",np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
        col = pyfits.Column(name="offset_sig", array=arr, format=columns["offset_sig"])    
        cols.append(col)
    
        arr = ff[extension].data["eligible_sky_density"]
        print("eligible_sky_density",np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr), min(arr), max(arr))
        
        arr = ff[extension].data["eligible_sky_density"] * ff[extension].data["match_dist"]**2 * 4*np.pi/3600
        print("expected_rnd",np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
        col = pyfits.Column(name="expected_rnd", array=arr*10, format=columns["match_dist"])    
        cols.append(col)        
        

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
            plt.scatter(ff[extension].data["offset_sig"], val)
            plt.xlabel("offset sig")
            plt.ylabel("pos")
            plt.title(ifn)
            plt.show()
    print("pos",np.nanmean(val), np.nanmedian(val), np.nanstd(val))
    col = pyfits.Column(name="pos", array=val, format=columns["match_dist"])
    #col = pyfits.Column(name="pos", array=(val*6)**2, format=columns["match_dist"])    
    cols.append(col)
    
    
    
    
    arr = np.log10(ff[extension].data["parallax"])
    gi = np.where(np.isnan(arr))[0]
    arr[gi] = -1
    print("plx",np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
    col = pyfits.Column(name="log_plx", array=arr*15, format=columns["parallax"])    
    cols.append(col)
    
    arr = np.log10(ff[extension].data["eligible_sky_density"])
    print("log sky_density",np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
    col = pyfits.Column(name="log_skd", array=arr, format=columns["eligible_sky_density"])    
    cols.append(col)

    arr = ff[extension].data["eligible_sky_density"]
    print("sky_density",np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
    col = pyfits.Column(name="skd", array=arr*10, format=columns["eligible_sky_density"])    
    cols.append(col)
        
    #if "category" in ff[extension].data.columns.names:
        #sigma = np.sqrt(-2*np.log(1-val))
        #print(sigma)
        #col = pyfits.Column(name="offset_sig", array=sigma, format=columns["offset_sig"])    
        #cols.append(col)
        #plt.scatter(np.exp(val), sigma)
        #plt.xlabel("pos")
        #plt.ylabel("offset_sig")
        #plt.title(ifn)
        #plt.show()
    #else:
        #arr = ff[extension].data["offset_sig"]
        #print("offset_sig",np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
        #col = pyfits.Column(name="offset_sig", array=arr, format=columns["offset_sig"])    
        #cols.append(col)
        
    
    
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

def prepare_training(ifn, ofn, overwrite=True, verbose=1):
    """
    """
    ext = 1
    ff = pyfits.open(ifn)
    cols = []
    column_formats = {col.name:col.format for col in ff[ext].columns}

    Lx = 4*np.pi*ff[ext].data["Fx"]* (1000/ff[ext].data["parallax"] * 3.1e18)**2
    gi = np.where((ff[ext].data["RADEC_ERR"] < 10) & (ff[ext].data["match_dist"]<25) &  (ff[ext].data["match_dist"]<4*ff[ext].data["RADEC_ERR"]) & (ff[ext].data['phot_g_mean_mag']>5.5))[0]

    if "category" in ff[ext].data.columns.names:
        print("Before'training_filter': ",len(gi))

        dct = {}
        for kw in ["parallax","Fx", "FxFg", "bp_rp"]:
            dct[kw] = ff[ext].data[kw][gi]
        ai = training_filter(dct)    
        pi = np.where(ai==0)[0]
        
        #gi = gi[pi]
        print("After 'training_filter': ",len(gi))
        #ff[extension].data["category"][gi] = 1
        
    for colname in ["srcID", "srcID_NN", "original_srcID"]:        
        arr = ff[ext].data[colname][gi]
        col = pyfits.Column(name=colname , array=arr, format=column_formats[colname])    
        cols.append(col)

    colname = "match_dist"
    arr = ff[ext].data[colname][gi]
    print(colname ,np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
    col = pyfits.Column(name=colname , array=arr, format=column_formats[colname])    
    cols.append(col)
    
    colname = "RADEC_ERR"
    arr = ff[ext].data[colname][gi]
    print(colname ,np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
    col = pyfits.Column(name=colname , array=arr, format=column_formats[colname])    
    cols.append(col)
    col = pyfits.Column(name="RADEC_sigma" , array=calc_sigma_from_RADEC_ERR(arr), format=column_formats[colname])    
    cols.append(col)
    
    arr0 = ff[ext].data["match_dist"][gi]
    arr1 = calc_sigma_from_RADEC_ERR(ff[ext].data["RADEC_ERR"][gi])
    col = pyfits.Column(name="offset_in_sigma" , array=arr0/arr1, format=column_formats["match_dist"])    
    cols.append(col)
    
    arr0 = ff[ext].data["match_dist"][gi]
    arr1 = ff[ext].data["eligible_sky_density"][gi]
    col = pyfits.Column(name="expected_rnd" , array=4*np.pi*arr0*arr1, format=column_formats["match_dist"])    
    cols.append(col)
    
    
    colname = "eligible_sky_density"
    arr = ff[ext].data[colname][gi]
    print(colname ,np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
    col = pyfits.Column(name="skd" , array=arr, format=column_formats[colname])    
    cols.append(col)
    
    arr = np.log10(ff[ext].data["eligible_sky_density"][gi])
    print("log sky_density",np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
    col = pyfits.Column(name="log_skd", array=arr, format=column_formats["eligible_sky_density"])    
    cols.append(col)

    colname = "category"
    if colname in ff[ext].data.columns.names:
        arr = ff[ext].data[colname][gi]
        print("Category: ")
        for cat in np.unique(arr):
            ci = np.where(arr==cat)[0]
            print("   #cat=",cat,":", len(ci))
        col = pyfits.Column(name=colname , array=arr, format=column_formats[colname])    
        cols.append(col)


    colname = "bp_rp"
    arr = ff[ext].data[colname][gi]
    print(colname ,np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
    col = pyfits.Column(name=colname , array=arr, format=column_formats[colname])    
    cols.append(col)
    
    colname = "Fx"
    arr = ff[ext].data[colname][gi]
    print(colname ,np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
    col = pyfits.Column(name=colname , array=np.log10(arr)+13, format=column_formats[colname])    
    cols.append(col)
    
    colname = "Fg"
    arr = ff[ext].data[colname][gi]
    print(colname ,np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
    col = pyfits.Column(name=colname , array=np.log10(arr)+11, format=column_formats[colname])    
    cols.append(col)
    
    colname = "FxFg"
    arr = ff[ext].data[colname][gi]
    print(colname ,np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
    if "category" in ff[ext].data.columns.names:
        cat = ff[ext].data["category"][gi]
        ri = np.where((cat >0) & (arr < 1e-5))[0]
        #arr[ri] = 0.01
    col = pyfits.Column(name=colname , array=np.log10(arr), format=column_formats[colname])    
    cols.append(col)
    
    colname = "parallax"
    arr = ff[ext].data[colname][gi]
    print(colname ,np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
    col = pyfits.Column(name="log_plx" , array=np.log10(arr), format=column_formats[colname])    
    cols.append(col)

    colname = "parallax"
    arr = ff[ext].data[colname][gi]
    print(colname ,np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
    col = pyfits.Column(name="plx", array=arr, format=column_formats[colname])    
    cols.append(col)

    colname = "parallax_error"
    arr = ff[ext].data[colname][gi]
    print(colname ,np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
    col = pyfits.Column(name="plx_error", array=arr, format=column_formats[colname])    
    cols.append(col)

    colname = "phot_g_mean_mag"
    arr = ff[ext].data[colname][gi]
    print(colname ,np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
    col = pyfits.Column(name="Gmag", array=arr, format=column_formats[colname])    
    cols.append(col)

    #colname = "phot_g_mean_mag"
    #arr = ff[ext].data[colname][gi]
    #print(colname ,np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
    #col = pyfits.Column(name="Gmag", array=arr, format=column_formats[colname])    
    #cols.append(col)

    colname = "NN"
    arr = ff[ext].data[colname][gi]
    print(colname ,np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
    col = pyfits.Column(name=colname , array=arr, format=column_formats[colname])    
    cols.append(col)
    
    colname = "RA"
    arr = ff[ext].data[colname][gi]
    print(colname ,np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
    col = pyfits.Column(name=colname , array=arr, format=column_formats[colname])    
    cols.append(col)
    
    colname = "Dec"
    arr = ff[ext].data[colname][gi]
    print(colname ,np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
    col = pyfits.Column(name=colname , array=arr, format=column_formats[colname])    
    cols.append(col)
    
    colname = "RA_NN"
    arr = ff[ext].data[colname][gi]
    print(colname ,np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
    col = pyfits.Column(name="gaia_RA" , array=arr, format=column_formats[colname])    
    cols.append(col)
    
    colname = "Dec_NN"
    arr = ff[ext].data[colname][gi]
    print(colname ,np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
    col = pyfits.Column(name="gaia_Dec", array=arr, format=column_formats[colname])    
    cols.append(col)
    print("len ",len(arr))
    #colname = "Dec_NN"
    #arr = ff[ext].data[colname][gi]
    #print(colname ,np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
    #col = pyfits.Column(name=colname , array=arr, format=column_formats[colname])    
    #cols.append(col)
    
        
    hdu = pyfits.PrimaryHDU()    
    cc = pyfits.ColDefs(cols)    
    xx = pyfits.BinTableHDU.from_columns(cc)
    dst = len(xx.data[cc[0].name])
    hdul = pyfits.HDUList([hdu, xx])
    print("XXXXXX")
    if ofn is not None:
        hdul.writeto(ofn, overwrite=overwrite)        
        if verbose>0:
            print("datasets::prep_classify - Written ",dst," objects with ",len(cols)," properties to ",ofn)
    ff.close()            
    return hdul

if __name__ == "__main__":
    #prepare_training("train.fits", "train_preprocessed.fits")
    #prepare_training("train2.fits", "train_preprocessed2.fits")
    #prepare_training("../ero_data/merged_random_eFEDS_EDR3.fits", "random4classify_eFEDS.fits")
    #prepare_training("../ero_data/merged_major_eFEDS_EDR3.fits", "major4classify_eFEDS.fits")
    prepare_training("../ero_data/merged_random_eFEDS_EDR3_HamStar.fits", "random4classify_eFEDS_HamStar.fits")
    prepare_training("../ero_data/merged_major_eFEDS_EDR3_HamStar.fits", "major4classify_eFEDS_HamStar.fits")
    prepare_training("train_HamStar.fits", "train_HamStar_preprocessed.fits")
    #prepare_training("train.fits", "train_preprocessed.fits")
    #prepare_training("../ero_data/merged_random_eFEDS_EDR3.fits", "random4classify_eFEDS.fits")
    #prepare_training("../ero_data/merged_major_eFEDS_EDR3.fits", "major4classify_eFEDS.fits")
    

