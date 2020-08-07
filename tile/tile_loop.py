from eroML.tile import Tile, add_healpix_col, extract_healpix_range
import eroML.config as conf
import numpy as np
import glob

def loop():
    fn = conf.ero_fn
    ofn = conf.ero_fn_hp
    NSIDE = conf.hp_nside
    print("using fn=",fn)
    rID = 1
    hps = add_healpix_col(fn, ofn=ofn, nside=NSIDE, overwrite=True)
    #hps = [500,501]
    #gi = np.where(hps==2814)[0][0]
    #print(gi)
    for i in hps:

        glob_str = "../ero_data/tile_hp"+str(i)+"_nside"+str(NSIDE)+"_rID"+str(rID)+"_training.fits"
        fnames = glob.glob(glob_str)
        print(glob_str)
        print(fnames)
        if len(fnames)>0: continue
        
        print()
        print()
        print(i)
        extract_healpix_range(ofn, "tmp_hp.fits", overwrite=True, min_index=i, max_index=i)
        t = Tile()
        t.populate_filenames("tmp_hp.fits", fn_prefix="../ero_data/tile_hp"+str(i)+"_nside"+str(NSIDE)+"_rID"+str(rID))
        #t.populate_filenames("tmp_hp.fits")
        t.prepare_data()
        t.generate_sets()
        print()
        print()        
#exit()
