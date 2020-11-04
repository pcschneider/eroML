from eroML.ensemble import Ensemble
from eroML.ensemble.tools import from_fits, to_fits
import configparser

fn = "config.ini"
config = configparser.ConfigParser()
config.read("config.ini")
ero_fn = config["Sources"]["ero_filename"]

e0 = from_fits(ero_fn, mapper={"detUID":"srcID", "DEC":"Dec"}, maxN=100)
e1 = from_fits("../eFEDS/SrcCat_V2T_shifted.fits", mapper={"detUID":"srcID", "DEC":"Dec"}, maxN=10)

e0.merge_add(e1, NN=3)

to_fits(e0, "test.fits", overwrite=True)
