import numpy as np
from eroML.utils import Dist_model
import matplotlib.pyplot as plt
from scipy.stats import norm
from eroML.config import *
import argparse
import os
from astropy.io import fits as pyfits
from eroML.positions import gen_random_pos_offset, gen_real_pos_offset


def generate_simu_data(mfn, ofn="test.dat", N=1000, rnd_factor=1,\
    dens_scaling=1.06, overwrite=False):

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
    except KeyError:
        SIG = ffd["RADEC_ERR"][i]
    else:
        print("Cannot read positional error from \'%s\', aborting..." % mfn)
    #SIG[SIG>11] = 11
    print("Mean, median \'sigma_r\': %6.3f, %6.3f [arcsec]" % (np.mean(SIG), np.median(SIG)))

    try:
        sk = ffd["skd"][i] / 10 # Because skd is scaled by a factor of ten in "prepare.py"
    except:
        sk = ffd["eligible_sky_density"][i]
    else:
        print("Cannot read sky density from \'%s\', aborting..." % mfn)        
    print("Mean, median \'skd\': %6.3f, %6.3f [eligible sources/arcmin^2]" % (np.mean(sk), np.median(sk)))
    
    md = ffd["match_dist"][i]
    print("Mean, median \'match_dist\': %6.3f, %6.3f [eligible sources/arcmin^2]" % (np.mean(md), np.median(md)))
    
    Nrnd = np.sum((sk/3600)*np.pi*md**2)
    print("Resulting in %i expected random sources (Nrandom/Nreal: %6.3f), simulating %i random sources." % (Nrnd, Nrnd/N, rnd_factor*Nrnd))
    Nrnd = round(N)

    ##sk_array = [0.0003, 0.0001, 0.001, 0.003, 0.01, 0.03, 0.1]
    #Nreal=N
    #Nrnd = 5000
    #SIG = 4
    #SIG = np.random.rand(Nreal)*20+1
    ##sk = np.random.rand(Nreal)**2*0.001+0.0001 # (mean eFEDS: 0.00014)
    #sk = np.random.rand(Nreal)*0.001+0.0001 # (mean eFEDS: 0.00014)
    ##sk = Nreal*[0.004]

    #rnd_offs = gen_random_pos_offset(dens=sk)
    real_offs = gen_real_pos_offset(sigma=SIG*0.6)
    rand_offs = gen_random_pos_offset(dens=np.repeat(sk/3600, rnd_factor)* dens_scaling)
    print(np.shape(rand_offs))
    for i in range(4):
        print(i, len(rand_offs[i]))
    #sk_simu = np.repeat(sk,3)
    idx = rand_offs[2].astype(int)
    print(idx)
    sig = np.repeat(SIG, rnd_factor)[idx]
    sig_simu = np.concatenate((SIG,sig))



    #Nrnd = len(rand_offs)
    #print("Simulating %i and %i real and random sources, respectively." % (Nreal, Nrnd))
    #print(len(sk_simu))


    #pos_off = np.zeros((Nreal+Nrnd)
    #skdens = np.zeros((Nreal+Nrnd)*len(sk_array))
    cls = np.zeros(len(sig_simu))
    cls[0:N] = 0
    cls[N::] = 1

    #tmp0 = gen_real_pos_offset(N, sig=SIG)
    #tmp1 = gen_random_pos_offset(Nrnd, dens=sk_simu)
    pos_off = np.concatenate((real_offs, rand_offs[0]))
    skdens = np.concatenate((sk, rand_offs[1]*3600))
    sigout = sig_simu
    nth = np.ones(len(sigout))
    nth[N:] = rand_offs[3]
    
    oo = np.transpose([sigout, pos_off, skdens, nth, cls])
    print(oo[0].shape, oo[1].shape, oo[2].shape, oo[3].shape)
    
    np.savetxt("offs.dat", oo)
    return
    #exit()
    #print(offs)
    #plt.hist(offs)

    #dm = Dist_model()

    ##"fraction","sig","dens","N", "err_scaling"]


    #dm["fraction"] = 1.0
    #dm["N"] = 1
    #dm["sig"] = SIG
    #dm["err_scaling"] = 1.

    #x = np.linspace(0, 5*SIG, 1000)
    #y = dm.evaluate(x)
    #plt.plot(x,y*N/2.4)
    #plt.show()
    ##plt.plot(x,np.cumsum(y)/np.sum(y))
    ##plt.show()


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
