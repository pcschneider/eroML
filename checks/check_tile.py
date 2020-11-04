from eroML.tile import Tile, add_healpix_col, extract_healpix_range

add_healpix_col("tmp_efeds.fits", ofn="tmp_efeds_hp.fits", nside=16, overwrite=True)
extract_healpix_range("tmp_efeds_hp.fits", "tst_hp.fits", overwrite=True, min_index=1400, max_index=1400)
#exit()

t = Tile()
t.populate_filenames("tst_hp.fits")
#t.populate_filenames("tmp_efeds.fits")
t.prepare_data()
t.generate_sets()
