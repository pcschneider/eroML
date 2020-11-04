from eroML import *
import configparser

fn = "config.ini"
config = configparser.ConfigParser()
config.read("config.ini")
ero_fn = config["Sources"]["ero_filename"]

print(ero_fn)
