import numpy as np
from eroML.utils import Dist_model
import matplotlib.pyplot as plt
from scipy.stats import norm
from eroML.config import *
import argparse
import os
from astropy.io import fits as pyfits

def gen_real_pos_offset(N, sig=1.):
        
    rnd = np.random.rand(N)
    rndx = np.sqrt(2*sig**2 * (-np.log(1-rnd)))

    return rndx


#def gen_random_pos_offset(Nrnd, dens=1.):
    #"""
    #Parameters
    #-----------
    #dens : array of float
        #Density (arcsec-2) of stars

    #Returns
    #-------
    #dens_simu, offs
    #"""
    #if type(dens) != float:
        #dens = np.array(dens)
    #max_dist = np.repeat(np.sqrt(NN/np.pi/dens), NN)
    #dd = np.repeat(dens,NN)    
    #N = len(dens)
    #print(N, np.shape(max_dist))
    #offs = max_dist*np.random.rand(N*NN)**0.5
    #print(np.shape(offs))
    #return offs



def gen_random_pos_offset(dens=1., nthN=3):
    """
    Parameters
    -----------
    dens : array of float
        Density (arcsec-2) of stars

    Returns
    -------
    dens_simu, offs
    """
    if type(dens) != float:
        dens = np.array(dens)
    N = len(dens)
    NN = 20
    max_dist = np.repeat(np.sqrt(NN/np.pi/dens), NN).reshape((N,NN))
    #print(np.shape(max_dist), max_dist)
    dd = np.repeat(dens,NN).reshape((N,NN))
    gg = np.repeat(np.arange(len(dens)),NN).reshape((N,NN))
    print(np.shape(gg), np.shape(dd))
    #print(N, np.shape(max_dist), max(max_dist))
    offs = max_dist*np.random.rand(N,NN)**0.5
    plt.hist(offs.flatten(), range=(0, 2.3), bins=30)
    #print(offs[0,:])
    plt.show()
    #print(np.shape(offs))
    nth = np.argsort(offs, axis=1)
    #offs = offs[nth]
    offs = np.sort(offs, axis=1)
    #nth = np.tile(np.arange(NN),N).flatten()[np.argsort(offs).flatten()]+1
    #print(nth[0,:])
    #print(offs[0,nth[0,:]-1])
    #offs = offs.flatten()
    #gi = np.where(nth.flatten() == 1)
    #print(offs)
    #print(nth)
    #print(offs[nth])
    #print(np.shape(offs))
    return np.array([offs.flatten(), dd.flatten(), gg.flatten()])




def gen_random_pos_offset2(dens=1., max_dist=70, sigma_in=None):
    """
    Parameters
    -----------
    dens : array of float
        Density (arcsec-2) of stars

    Returns
    -------
    dens_simu, offs
    """
    if type(dens) != float:
        dens = np.array(dens)
    NN = 3*len(dens)
    N_simu = np.random.poisson(lam=NN)    
    dens_simu = np.repeat(dens, N_simu)
    sigma_simu = np.repeat(sigma_in, N_simu)
    offs = max_dist*np.random.rand(np.sum(N_simu))**0.5
    return sigma_simu, dens_simu, offs



def generate_simu_data(mfn, ofn="test.dat", N=1000, rnd_factor=1,  overwrite=False):

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
    real_offs = gen_real_pos_offset(N, sig=SIG)
    rand_offs = gen_random_pos_offset(Nrnd, dens=sk)
    #sk_simu = np.repeat(sk,3)
    #sig_simu = np.repeat(SIG,3)



    #Nrnd = len(rand_offs)
    #print("Simulating %i and %i real and random sources, respectively." % (Nreal, Nrnd))
    #print(len(sk_simu))


    #pos_off = np.zeros((Nreal+Nrnd)
    #skdens = np.zeros((Nreal+Nrnd)*len(sk_array))
    cls = np.zeros(N+Nrnd)
    cls[0:N] = 0
    cls[N::] = 1

    #tmp0 = gen_real_pos_offset(N, sig=SIG)
    #tmp1 = gen_random_pos_offset(Nrnd, dens=sk_simu)
    pos_off = np.concatenate((real_offs, rand_offs))
    skdens = np.concatenate((sk, sk_simu))
    sigout = np.concatenate((SIG, sig_simu))   
        
    np.savetxt("offs.dat", np.transpose([sigout, pos_off, skdens*10000, cls]))
    exit()
    print(offs)
    plt.hist(offs)

    dm = Dist_model()

    #"fraction","sig","dens","N", "err_scaling"]


    dm["fraction"] = 1.0
    dm["N"] = 1
    dm["sig"] = SIG
    dm["err_scaling"] = 1.

    x = np.linspace(0, 5*SIG, 1000)
    y = dm.evaluate(x)
    plt.plot(x,y*N/2.4)
    plt.show()
    #plt.plot(x,np.cumsum(y)/np.sum(y))
    #plt.show()


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