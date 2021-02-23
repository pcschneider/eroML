import numpy as np
from math import floor, ceil
from eroML.utils import Dist_model
import matplotlib.pyplot as plt
from scipy.stats import norm
from eroML.config import *
import argparse
import os
from astropy.io import fits as pyfits
from eroML.positions import gen_random_pos_offset, gen_real_pos_offset, random4dens, Nexp4dens


def generate_simu_data(mfn, ofn="test.dat", N=1000, rnd_factor=1,\
    dens_scaling=1.03, overwrite=False, key=1, max_dist=60):
    """
    
    Parameters
    ----------
    mfn : str
        Filename for major-file
    ofn : str
        Filename for file with simulated offsets (will be generated)
    N : int 
        Number of real sources
    rnd_factor : float
        Scaling factor for the generation of random sources
    dens_scaling : float
        Scaling of the calculated sky density for simulating random sources
    max_dist : float
        Return sources expected within ``max_dist'' (unit arcsec)
    overwrite : bool
        Overwrite `ofn` if it exists
    key : int or float
        Appended to the data file to identify this particular simulation
    """

    print("Using major file \'%s\' to simulate %i real sources (rnd_factor=%6.2f); ofn=\'%s\' (overwrite=%i)." % (mfn, N, rnd_factor, ofn, overwrite))

    if os.path.isfile(ofn) and overwrite==False:
        print("ofn=\'%s\' exists and \'overwrite\'=%i. Aborting..." % (ofn, overwrite))

    ff = pyfits.open(mfn)
    ffd = ff[1].data
    ffN = len(ffd["srcID"])
    gi = np.where(ffd["NN"]==1)[0]
    
    ffNu = len(gi)
    print("Number of sources in \'%s\': %i (NN=1: %i)" % (mfn, ffN, ffNu))
    if N>ffNu:
        print("WARNING - Number of simulated sources larger than number of sources in major-file.")
        
    i = np.random.choice(gi, size=N)
    try:
        SIG = ffd["sigma_r"][i]
        print("XXX")
    except KeyError:
        SIG = ffd["RADEC_ERR"][i]
    except:
        print("Cannot read positional error from \'%s\', aborting..." % mfn)
    #SIG[SIG>11] = 11
    print("Mean, median \'sigma_r\': %6.3f, %6.3f [arcsec]" % (np.mean(SIG), np.median(SIG)))

    try:
        sk = ffd["skd"] / 10 # Because skd is scaled by a factor of ten in "prepare.py"
    except KeyError:
        sk = ffd["eligible_sky_density"]
    except:
        print("Cannot read sky density from \'%s\', aborting..." % mfn)        
    print("Mean, median \'skd\': %6.3f, %6.3f [eligible sources/arcmin^2]" % (np.mean(sk), np.median(sk)))
    
    md = ffd["match_dist"][i]
    print("Mean, median \'match_dist\': %6.3f, %6.3f [eligible sources/arcmin^2]" % (np.mean(md), np.median(md)))
    
    Nrnd = np.sum((sk[i]/3600)*np.pi*md**2)
    print("Resulting in %i expected random sources (Nrandom/Nreal: %6.3f), simulating %i random sources." % (Nrnd, Nrnd/N, rnd_factor*Nrnd))
    Nrnd = round(N)


    #rnd_offs = gen_random_pos_offset(dens=sk)
    real_offs = gen_real_pos_offset(sigma=SIG*0.6)
    
    dens = sk[i]
    Nexp = Nexp4dens(dens*dens_scaling, max_dist=max_dist)
   
    Ndens = floor(Nexp * rnd_factor)
    #ggg = np.where()
    ii = np.random.choice(len(sk), Ndens)
    dens = sk[ii]
    print("len dens ", len(dens))
    
    rand_offs = random4dens(dens=dens* dens_scaling)
    
    
    print(np.shape(rand_offs))
    #for j in range(4):
        #print(j, len(rand_offs[j]))
    #sk_simu = np.repeat(sk,3)
    idx = rand_offs[2].astype(int)
    print(idx)
    sigi = np.random.choice(ffN, len(rand_offs[0]))
    #sig = ffd[
    #i = np.random.choice(gi, size=N)
    try:
        sig = ffd["sigma_r"][sigi]
        print("XXX")
    except KeyError:
        sig = ffd["RADEC_ERR"][sigi]
        
    sig_simu = np.concatenate((SIG,sig))


    #pos_off = np.zeros((Nreal+Nrnd)
    #skdens = np.zeros((Nreal+Nrnd)*len(sk_array))
    cls = np.zeros(len(sig_simu))
    cls[0:N] = 0
    cls[N::] = 1

    #tmp0 = gen_real_pos_offset(N, sig=SIG)
    #tmp1 = gen_random_pos_offset(Nrnd, dens=sk_simu)
    pos_off = np.concatenate((real_offs, rand_offs[0]))
    skdens = np.concatenate((sk[i], rand_offs[1]))
    sigout = sig_simu
    nth = np.ones(len(sigout))
    nth[N:] = rand_offs[3]
    k = np.repeat([key], len(nth))
    
    oo = np.transpose([sigout, pos_off, skdens, nth, cls, k])
    print(oo[0].shape, oo[1].shape, oo[2].shape, oo[3].shape)
    
    np.savetxt(ofn, oo)
    return
  

if __name__ == "__main__":
    unlogged = []
    
    parser = argparse.ArgumentParser(usage='%(prog)s [options]\n\n major-filename overrides the options in config-filename')
    parser.add_argument("--conf", nargs='?', default=None, help="Config-filename")
    parser.add_argument("--major", nargs='?', default=None, help="major-filename")
    parser.add_argument("--ofn", nargs='?', default="test.dat", help="file with simulated offsets.")
    parser.add_argument("--rnd_factor", nargs='?', type=float, default=1.0, help="Factor for increasing the Number of simulated random sources.")
    parser.add_argument("-o", dest='overwrite', action='store_true', default=False, help="Overwrite file with simulated data")
    
    parser.add_argument("N", nargs='?', type=int, default=1000, help="Number of real random sources")
    
    args = parser.parse_args()
    
    if args.major is not None:
        mfn = args.major
        if os.path.isfile(mfn) is False:
            print("Major file \'%s\' does not exist. Aborting..." % mfn)
            exit()
    
        generate_simu_data(mfn, ofn=args.ofn, N=args.N, rnd_factor=args.rnd_factor, overwrite=args.overwrite)
        
    elif args.conf is not None:
        unlogged.append(("info", "Using custom config-file: `%s` " % args.conf))        
        tmp = read_config(args.conf)
        unlogged.append(tmp)
        for ll in unlogged:
            print(ll[0].upper(), ": ", ll[1])
       
        mfn = config["Classification"]["major_filename"]
        if os.path.isfile(mfn) is False:
            print("Major file \'%s\' does not exist. Aborting..." % mfn)
            exit()
            
        generate_simu_data(mfn, ofn=args.ofn, N=args.N, rnd_factor=args.rnd_factor, overwrite=args.overwrite)
    else:
        print("Nothing to do...")
