# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 15:12:08 2024

@author: cammy
"""

from os import path
import matplotlib.pyplot as plt
import numpy as np

# All packages below are available at https://pypi.org/user/micappe/

from plotbin.plot_velfield import plot_velfield
from plotbin.symmetrize_velfield import symmetrize_velfield

import jampy as jam_package
from jampy.jam_axi_proj import jam_axi_proj
from jampy.mge_half_light_isophote import mge_half_light_radius
from jampy.mge_radial_mass import mge_radial_mass

from adamet.adamet import adamet
from adamet.corner_plot import corner_plot
from astropy.io import fits
from mgefit.mge_fit_1d import mge_fit_1d
from pafit.fit_kinematic_pa import fit_kinematic_pa

import emcee

from multiprocessing import Pool
import time
import corner

def dark_halo_gNFW_MGE(gamma,lg_rbreak,lg_rho_s):
    start=time.time()
    rbreak,rho_s=10**lg_rbreak,10**lg_rho_s
    rbreak /= pc
    
    n = 300     # Number of values to sample the gNFW profile for the MGE fit
    r = np.geomspace(1, rbreak*10, n)
    
    
    
    rho=(rho_s)/(((r/rbreak)**gamma)*((1+(r/rbreak)**2)**((3-gamma)/2.)))
    
    #rho = (r/rbreak)**gamma * (0.5 + 0.5*r/rbreak)**(-gamma - 3)  # rho=1 at r=rbreak
    
    
    m = mge_fit_1d(r, rho, ngauss=20, quiet=1, plot=0)
    #plt.pause(1)
    surf_dm, sigma_dm = m.sol
    
    
    
    qobs_dm= np.ones_like(surf_dm)
    
    #surf_pot_dm=surf_dm/(sigma_dm*np.sqrt(2*np.pi))
    
    end = time.time()
    serial_time = end - start
    #print("dark mge took {0:.1f} seconds".format(serial_time))
    return surf_dm,sigma_dm, qobs_dm

def total_mass_mge(surf_mass,surf_lum, sigma_lum, qobs_lum, gamma, lg_rbreak,lg_rho_s, inc, alpha):
    
    start=time.time()
    
    surf_dm,sigma_dm,qobs_dm= dark_halo_gNFW_MGE(gamma,lg_rbreak,lg_rho_s)
    #print(surf_lum,'wjuehfuheu')
    surf_mass *= alpha
    
    if any(np.isnan(surf_lum))==1 or any(np.isnan(surf_dm))==1:
        print('shit',surf_lum,surf_dm)
    
    surf_pot = np.append(surf_mass,surf_dm)
    sigma_pot= np.append(sigma_lum,sigma_dm)
    qobs_pot= np.append(qobs_lum,qobs_dm)
    
    
    #resets Surf lum as it is mulitplied by alpha and then does not get
    #Set again before running the next MC run
    #surf_lum=get_lum_mge()[0]
    end = time.time()
    serial_time = end - start
    #print("total mass mge took {0:.1f} seconds".format(serial_time))
    
    return surf_pot,sigma_pot,qobs_pot






def JAM_model(pars,surf_mass=None, surf_lum=None, sigma_lum=None, qobs_lum=None, distance=None,
              xbin=None, ybin=None, sigmapsf=None, normpsf=None,
              rms=None, erms=None, pixsize=None,goodbins=None):
    
    start=time.time()
    #defines the parameters used
    alpha,lg_rbreak,lg_rho_s,gamma,lg_mbh = pars
    
    #print(pars,'check')
    #print(rms,'check')
    inc=60
    mbh=10**lg_mbh
    #Need to fit the dark matter density fit mge 1d fit to add the surface 
    #density to the mass mge of the luminous part. The dark matter density is 
    #paramerterised by lg_rbreak, lg_rho_s and gamma.
    
    #total mass mge appends the dark mge's to the mass mge's to get a 1D array
    #of all the mge's
    
    surf_mass=get_lum_mge()[3]
    
    #surf_lum *= alpha
    
    surf_pot,sigma_pot,qobs_pot = total_mass_mge(surf_mass,surf_lum, sigma_lum, qobs_lum, gamma, lg_rbreak,lg_rho_s, inc, alpha)
    #Note: In the example this is where surf pot is scaled by M/L including the dark matter mge's,
    #This is not fully explained but is said to keep the fraction of the dark matter the same, this
    #also scales Mbh
    
    
    #print(surf_lum)
    #Runs the actual JAM model to get the model RMS out (other moments can be selected)
    jam = jam_axi_proj(surf_lum, sigma_lum, qobs_lum, surf_pot, sigma_pot, qobs_pot,
                       inc, mbh, distance, xbin, ybin, plot=0, pixsize=0, quiet=1,
                       align='cyl', data=rms, errors=erms, ml=None,goodbins=goodbins,interp=False)
    
    
    
    Vrms_model=jam.model
    
    if any(np.isnan(Vrms_model))==1:
        print('fuck',surf_lum,Vrms_model,erms,pars)
    
    end = time.time()
    serial_time = end - start
    #print("JAM_model took {0:.1f} seconds".format(serial_time))
    return Vrms_model







    
def lnlike(pars,surf_mass=None, surf_lum=None, sigma_lum=None, qobs_lum=None, distance=None,
              xbin=None, ybin=None, sigmapsf=None, normpsf=None,
              rms=None, erms=None, pixsize=None,goodbins=None):
    start=time.time()
    
    Vrms_model=JAM_model(pars, surf_mass,surf_lum, sigma_lum, qobs_lum, distance,
              xbin, ybin, sigmapsf, normpsf, 
              rms, erms, pixsize,goodbins)
    
    lnlike= -0.5*sum(((rms-Vrms_model)/erms)**2)
    
    if np.isnan(lnlike)==1:
        print('fuck',surf_lum,Vrms_model,pars)
    end = time.time()
    serial_time = end - start
    #print("lnlike took {0:.1f} seconds".format(serial_time))
    
    return lnlike

def lnprior(pars):
    start=time.time()
    check=np.zeros(len(pars))
    #alpha,lg_rbreak,lg_rho_s,gamma,lg_mbh = pars
    
    for i in range(len(pars)):
        if pars[i]>=bounds[0][i] and pars[i]<=bounds[1][i]:
            check[i]=0
        else:
            check[i]=1
    
    if check.any()==1:
        lp=-np.inf
    else:
        lp=0.0
    
    end = time.time()
    serial_time = end - start
    #print("lnprior took {0:.1f} seconds".format(serial_time))
    return lp
    
def lnprob(pars, surf_mass=None, surf_lum=None, sigma_lum=None, qobs_lum=None, distance=None,
              xbin=None, ybin=None, sigmapsf=None, normpsf=None, 
              rms=None, erms=None, pixsize=None,goodbins=None):
    #print(surf_lum,'check lnlike')
    start = time.time()
    lp= lnprior(pars)
    if lp != 0.0:
        #print('fuuuuckkk')
        return -np.inf
    else:
        return lp + lnlike(pars,surf_mass,surf_lum, sigma_lum, qobs_lum, distance,
              xbin, ybin, sigmapsf, normpsf, 
              rms, erms, pixsize,goodbins)
    end = time.time()
    serial_time = end - start
    #print("lnprob took {0:.1f} seconds".format(serial_time))

def get_fitz():
    hdul = fits.open("PXF_bin_MS_NGC2974_r5_idl.fits")
    hdul.info()
    #print(hdul[1].header)
    data=hdul[1].data
    xbin,ybin, V, sigma,flux,Verr, serr = hdul[1].data['XS'],hdul[1].data['YS'],hdul[1].data['VPXF'],hdul[1].data['SPXF'],hdul[1].data['FLUX'],hdul[1].data['EVPXF'],hdul[1].data['ESPXF']
    Vhel,distance=1887.,20.89
    rms=np.sqrt((V-Vhel)**2+sigma**2)
    erms=(1./rms)*np.sqrt((V*Verr)**2+(sigma*serr)**2)
    return xbin,ybin,rms,erms,distance,V

def get_lum_mge():
    surf_lum=np.asarray([4276.01,7782.37,2853.55,3171.34,220.000,970.160,252.150])
    sigma_lum=np.asarray([0.54153,0.88097,1.44526,3.81993,6.64704,10.7437,28.4453])
    qobs_lum=np.asarray([0.83144,0.82501 ,0.94271 ,0.67267 ,0.99990 ,0.55375,0.61238])
    
    surf_mass=np.asarray([16208.47,26366.23,13148.71,11329.50,1966.17, 2890.09,778.71])
    return surf_lum,sigma_lum,qobs_lum,surf_mass

def transform(x,y,angle):
    theta=np.radians(angle)
    xm=x*np.cos(theta)-y*np.sin(theta)
    ym=x*np.sin(theta)+y*np.cos(theta)
    return xm,ym


def dark_halo_gNFW_MGE(gamma,lg_rbreak,lg_rho_s):
    start=time.time()
    rbreak,rho_s=10**lg_rbreak,10**lg_rho_s
    rbreak /= pc
    
    n = 300     # Number of values to sample the gNFW profile for the MGE fit
    r = np.geomspace(1, rbreak*10, n)
    
    
    
    rho=(rho_s)/(((r/rbreak)**gamma)*((1+(r/rbreak)**2)**((3-gamma)/2.)))
    
    #rho = (r/rbreak)**gamma * (0.5 + 0.5*r/rbreak)**(-gamma - 3)  # rho=1 at r=rbreak
    
    
    m = mge_fit_1d(r, rho, ngauss=20, quiet=1, plot=0)
    #plt.pause(1)
    surf_dm, sigma_dm = m.sol
    
    
    
    qobs_dm= np.ones_like(surf_dm)
    
    #surf_pot_dm=surf_dm/(sigma_dm*np.sqrt(2*np.pi))
    
    end = time.time()
    serial_time = end - start
    #print("dark mge took {0:.1f} seconds".format(serial_time))
    return surf_dm,sigma_dm, qobs_dm

def total_mass_mge(surf_mass,surf_lum, sigma_lum, qobs_lum, gamma, lg_rbreak,lg_rho_s, inc, alpha):
    
    start=time.time()
    
    surf_dm,sigma_dm,qobs_dm= dark_halo_gNFW_MGE(gamma,lg_rbreak,lg_rho_s)
    #print(surf_lum,'wjuehfuheu')
    surf_mass *= alpha
    
    if any(np.isnan(surf_lum))==1 or any(np.isnan(surf_dm))==1:
        print('shit',surf_lum,surf_dm)
    
    surf_pot = np.append(surf_mass,surf_dm)
    sigma_pot= np.append(sigma_lum,sigma_dm)
    qobs_pot= np.append(qobs_lum,qobs_dm)
    
    
    #resets Surf lum as it is mulitplied by alpha and then does not get
    #Set again before running the next MC run
    #surf_lum=get_lum_mge()[0]
    end = time.time()
    serial_time = end - start
    #print("total mass mge took {0:.1f} seconds".format(serial_time))
    
    return surf_pot,sigma_pot,qobs_pot






def JAM_model(pars,surf_mass=None, surf_lum=None, sigma_lum=None, qobs_lum=None, distance=None,
              xbin=None, ybin=None, sigmapsf=None, normpsf=None,
              rms=None, erms=None, pixsize=None,goodbins=None):
    
    start=time.time()
    #defines the parameters used
    alpha,lg_rbreak,lg_rho_s,gamma,lg_mbh = pars
    
    #print(pars,'check')
    #print(rms,'check')
    inc=60
    mbh=10**lg_mbh
    #Need to fit the dark matter density fit mge 1d fit to add the surface 
    #density to the mass mge of the luminous part. The dark matter density is 
    #paramerterised by lg_rbreak, lg_rho_s and gamma.
    
    #total mass mge appends the dark mge's to the mass mge's to get a 1D array
    #of all the mge's
    
    surf_mass=get_lum_mge()[3]
    
    #surf_lum *= alpha
    
    surf_pot,sigma_pot,qobs_pot = total_mass_mge(surf_mass,surf_lum, sigma_lum, qobs_lum, gamma, lg_rbreak,lg_rho_s, inc, alpha)
    #Note: In the example this is where surf pot is scaled by M/L including the dark matter mge's,
    #This is not fully explained but is said to keep the fraction of the dark matter the same, this
    #also scales Mbh
    
    
    #print(surf_lum)
    #Runs the actual JAM model to get the model RMS out (other moments can be selected)
    jam = jam_axi_proj(surf_lum, sigma_lum, qobs_lum, surf_pot, sigma_pot, qobs_pot,
                       inc, mbh, distance, xbin, ybin, plot=0, pixsize=0, quiet=1,
                       align='cyl', data=rms, errors=erms, ml=None,goodbins=goodbins,interp=False)
    
    
    
    Vrms_model=jam.model
    
    if any(np.isnan(Vrms_model))==1:
        print('fuck',surf_lum,Vrms_model,erms,pars)
    
    end = time.time()
    serial_time = end - start
    #print("JAM_model took {0:.1f} seconds".format(serial_time))
    return Vrms_model







    
def lnlike(pars,surf_mass=None, surf_lum=None, sigma_lum=None, qobs_lum=None, distance=None,
              xbin=None, ybin=None, sigmapsf=None, normpsf=None,
              rms=None, erms=None, pixsize=None,goodbins=None):
    start=time.time()
    
    Vrms_model=JAM_model(pars, surf_mass,surf_lum, sigma_lum, qobs_lum, distance,
              xbin, ybin, sigmapsf, normpsf, 
              rms, erms, pixsize,goodbins)
    
    lnlike= -0.5*sum(((rms-Vrms_model)/erms)**2)
    
    if np.isnan(lnlike)==1:
        print('fuck',surf_lum,Vrms_model,pars)
    end = time.time()
    serial_time = end - start
    #print("lnlike took {0:.1f} seconds".format(serial_time))
    
    return lnlike

def lnprior(pars):
    start=time.time()
    check=np.zeros(len(pars))
    #alpha,lg_rbreak,lg_rho_s,gamma,lg_mbh = pars
    
    for i in range(len(pars)):
        if pars[i]>=bounds[0][i] and pars[i]<=bounds[1][i]:
            check[i]=0
        else:
            check[i]=1
    
    if check.any()==1:
        lp=-np.inf
    else:
        lp=0.0
    
    end = time.time()
    serial_time = end - start
    #print("lnprior took {0:.1f} seconds".format(serial_time))
    return lp
    
def lnprob(pars, surf_mass=None, surf_lum=None, sigma_lum=None, qobs_lum=None, distance=None,
              xbin=None, ybin=None, sigmapsf=None, normpsf=None, 
              rms=None, erms=None, pixsize=None,goodbins=None):
    #print(surf_lum,'check lnlike')
    start = time.time()
    lp= lnprior(pars)
    if lp != 0.0:
        #print('fuuuuckkk')
        return -np.inf
    else:
        return lp + lnlike(pars,surf_mass,surf_lum, sigma_lum, qobs_lum, distance,
              xbin, ybin, sigmapsf, normpsf, 
              rms, erms, pixsize,goodbins)
    end = time.time()
    serial_time = end - start
    #print("lnprob took {0:.1f} seconds".format(serial_time))

def get_fitz():
    hdul = fits.open("PXF_bin_MS_NGC2974_r5_idl.fits")
    hdul.info()
    #print(hdul[1].header)
    data=hdul[1].data
    xbin,ybin, V, sigma,flux,Verr, serr = hdul[1].data['XS'],hdul[1].data['YS'],hdul[1].data['VPXF'],hdul[1].data['SPXF'],hdul[1].data['FLUX'],hdul[1].data['EVPXF'],hdul[1].data['ESPXF']
    Vhel,distance=1887.,20.89
    rms=np.sqrt((V-Vhel)**2+sigma**2)
    erms=(1./rms)*np.sqrt((V*Verr)**2+(sigma*serr)**2)
    return xbin,ybin,rms,erms,distance,V

def get_lum_mge():
    surf_lum=np.asarray([4276.01,7782.37,2853.55,3171.34,220.000,970.160,252.150])
    sigma_lum=np.asarray([0.54153,0.88097,1.44526,3.81993,6.64704,10.7437,28.4453])
    qobs_lum=np.asarray([0.83144,0.82501 ,0.94271 ,0.67267 ,0.99990 ,0.55375,0.61238])
    
    surf_mass=np.asarray([16208.47,26366.23,13148.71,11329.50,1966.17, 2890.09,778.71])
    return surf_lum,sigma_lum,qobs_lum,surf_mass

def transform(x,y,angle):
    theta=np.radians(angle)
    xm=x*np.cos(theta)-y*np.sin(theta)
    ym=x*np.sin(theta)+y*np.cos(theta)
    return xm,ym



xbin,ybin,rms,erms,distance,V=get_fitz()

surf_lum,sigma_lum,qobs_lum,surf_mass=get_lum_mge()


pixsize = 0.8# spaxel size in arcsec (before Voronoi binning)
sigmapsf = 2.6/2.355      # sigma PSF in arcsec (=FWHM/2.355)
normpsf = 1
pc = distance*np.pi/0.648


alpha0 = 1.8     # 
lg_rbreak0 = 4.7   # 
lg_rho_s = -2.25  # 
gamma0 = 0.      #
lg_mbh0= 9.5     #

#labels={'alpha','rbreak','density','gamma','Mbh'}


labels=[r"$alpha$",r"$rbreak$",r"$density$",r"$\gamma$",r"$Mbh$"]

vel_corr=V-np.median(V)

angBest, angErr, vSyst = fit_kinematic_pa(xbin, ybin, vel_corr, debug=False, plot=False,quiet=1)

xbin,ybin = transform(xbin,ybin,angBest)
goodbins = np.isfinite(xbin)
loc=0
for i in range(0,len(xbin)):
    if xbin[i]==0 and ybin[i]==0:
        print('0 at',i)
    else:
        #print(i)
        loc = np.append(loc,i)

print(loc)
xbin=xbin[loc]
ybin=ybin[loc]
rms=rms[loc]
erms=erms[loc]



#means = np.random.rand(ndim)
#print(means)
kwargs = {'surf_lum': surf_lum, 'sigma_lum': sigma_lum, 'qobs_lum': qobs_lum,
              'distance': distance, 'xbin': xbin, 'ybin': ybin, 'sigmapsf': sigmapsf,
              'normpsf': normpsf, 'rms': rms, 'erms': erms, 'pixsize': pixsize}

args=[surf_mass,surf_lum, sigma_lum, qobs_lum, distance,
              xbin, ybin, sigmapsf, normpsf, 
              rms, erms, pixsize,goodbins]

#print('huuhu',[surf_lum, sigma_lum, qobs_lum, distance,
#              xbin, ybin, sigmapsf, normpsf, 
#              rms, erms, pixsize])


nwalkers=50

niter=1000

initial = np.asarray([alpha0, lg_rbreak0, lg_rho_s, gamma0, lg_mbh0])

bounds = np.asarray([[1., 3., -5., 0.,6.], [2., 5., -1,1.1, 11.]])

ndim = len(initial)
p0 = [np.array(initial) + 1e-7 * np.random.randn(ndim) for i in range(nwalkers)]


def main(p0,nwalkers,niter,ndim,lnprob,args):
    
    #print('dieejijei',surf_lum)
    
    # Set up the backend
# Don't forget to clear it in case the file already exists
    filename = "tutorial.h5"
    backend = emcee.backends.HDFBackend(filename)
    backend.reset(nwalkers, ndim)

    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=args,backend=backend)
    start = time.time()
    print("Running burn-in...")
    p0, _, _ = sampler.run_mcmc(p0, 10,progress=False)
    sampler.reset()
    end = time.time()
    serial_time = end - start
    print("Serial took {0:.1f} seconds".format(serial_time))
        
        
        
        
    start = time.time()
    print("Running production...")
    pos, prob, state = sampler.run_mcmc(p0, niter,progress=False)
    end = time.time()
    serial_time = end - start
    print("Serial took {0:.1f} seconds".format(serial_time))
    return sampler, pos, prob, state

sampler, pos, prob, state = main(p0,nwalkers,niter,ndim,lnprob,args)



samples = sampler.flatchain

theta_max  = samples[np.argmax(sampler.flatlnprobability)]
print(theta_max)
best_fit_model = JAM_model(theta_max,*args)



fig,axs =plt.subplots(1,3,figsize=(19.5,8))
norm = plt.cm.colors.Normalize(vmin = min(rms), vmax = max(rms))
cmap = 'hot'
sc = plt.cm.ScalarMappable(norm = norm, cmap =  cmap)
norm2 = plt.cm.colors.Normalize(vmin = min(rms-best_fit_model), vmax = max(rms-best_fit_model))
sc2 = plt.cm.ScalarMappable(norm = norm2, cmap =  cmap)

a=axs[0].scatter(xbin,ybin,c=rms,norm=norm,cmap=cmap)
axs[0].set_title('rms data')
#plt.colorbar()
#plt.show()
b=axs[1].scatter(xbin,ybin,c=best_fit_model,norm=norm,cmap=cmap)
fig.colorbar(sc, ax = axs[1])
axs[1].set_title('rms Model')
#plt.show()
axs[2].scatter(xbin,ybin,c=abs(rms-best_fit_model),cmap=cmap)
axs[2].set_title('residuals')
fig.colorbar(sc2, ax = axs[2])
plt.savefig('NGC2974_model_data_rms.png')
plt.show()
fig = corner.corner(samples,show_titles=True,labels=labels,plot_datapoints=True,quantiles=[0.16, 0.5, 0.84])
plt.savefig('corner_plot_NGC2974.png')



