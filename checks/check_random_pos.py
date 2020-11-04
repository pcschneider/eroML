from eroML.tools import fake_positions
#from eroML.config import read_config
from eroML.utils import setup_logger
from eroML.config import *

read_config("tt2.ini")
#print(config)
logger = setup_logger(config["General"]["main_log_file"], config["General"]["debug_log_file"])

fake_positions(cconfig=config)
