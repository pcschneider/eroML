from eroML.tile import Tile, add_healpix_col, extract_healpix_range
import eroML.config as conf

def loop():
    fn = conf.ero_fn
    ofn = conf.ero_fn_hp
    NSIDE = conf.hp_nside
    print("using fn=",fn)
    rID = 1
    hps = add_healpix_col(fn, ofn=ofn, nside=NSIDE, overwrite=True)
    #hps = [500,501]
    for i in hps[100:102]:
        
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
