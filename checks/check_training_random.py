from utils.datasets import training_random_set

rel_path = "../ero_data/"
ifn1 = rel_path+"tile_hp8711_nside32_rID1_training.fits"
ifn2 = rel_path+"tile_hp8711_nside32_rID1_random.fits"
ofn = "test_tr.fits"
training_random_set(ifn1, ifn2, ofn)
