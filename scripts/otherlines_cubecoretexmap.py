import numpy as np
import astropy.units as u
from spectral_cube import SpectralCube as sc
import matplotlib.pyplot as plt
from astroquery.splatalogue import utils, Splatalogue
from astroquery.linelists.cdms import CDMS
import scipy.constants as cnst
from astropy.io import fits
import glob
import radio_beam
import os
from astropy.modeling import models, fitting
import time
import pdb
import pickle
from astropy.wcs import WCS
import matplotlib as mpl
import copy
from astropy import coordinates
from spectral_cube import BooleanArrayMask
from astropy.nddata import Cutout2D
from spectral_cube.io.casa_masks import make_casa_mask
from utilities import *
from math import log10, floor
from manualspeciesabundances import *
from pyspeckit.spectrum.models import lte_molecule

Splatalogue.QUERY_URL= 'https://splatalogue.online/c_export.php'

def round_to_1(x):
    return round(x, -int(floor(log10(abs(x)))))

'''This wing of the script takes in continuum-subtracted cubes, cuts out a subcube around a region of interest based on a DS9 region, and converts the subcubes into brightness temperature (K) units'''
print('Cube-->Core-->Tex start\n')
print('Begin Jy/beam-to-K and region subcube conversion\n')

source='DSi'
print(f'Source: {source}\n')
fields={'SgrB2S':1,'DSi':10,'DSii':10,'DSiii':10,'DSiv':10,'DSv':10,'DSVI':2,'DSVII':3,'DSVIII':3,'DSIX':7,'DSX':7,'DSXI':8}
fnum=fields[source]
molecule=' C2H5OH '
nospace_molecule=molecule.replace(' ','')

#inpath="/orange/adamginsburg/sgrb2/d.jeff/data/field10originalimages/"
inpaths={1:'/orange/adamginsburg/sgrb2/d.jeff/data/OctReimage_K/',10:"/orange/adamginsburg/sgrb2/d.jeff/data/field10originals_K/",2:"/orange/adamginsburg/sgrb2/d.jeff/data/field2originals_K/",3:"/orange/adamginsburg/sgrb2/d.jeff/data/field3originals_K/",7:"/orange/adamginsburg/sgrb2/d.jeff/data/field7originals_K/",8:"/orange/adamginsburg/sgrb2/d.jeff/data/field8originals_K/"}
inpath=inpaths[fnum]#'/blue/adamginsburg/d.jeff/imaging_results/data/OctReimage/'
beamcubes=glob.glob(inpath+'*.fits')
homes={1:'/orange/adamginsburg/sgrb2/d.jeff/products/OctReimage_K/',10:"/orange/adamginsburg/sgrb2/d.jeff/products/field10originals_K/",2:"/orange/adamginsburg/sgrb2/d.jeff/products/field2originals_K/",3:"/orange/adamginsburg/sgrb2/d.jeff/products/field3originals_K/",7:"/orange/adamginsburg/sgrb2/d.jeff/products/field7originals_K/",8:"/orange/adamginsburg/sgrb2/d.jeff/products/field8originals_K/"}
home=homes[fnum]#'/blue/adamginsburg/d.jeff/imaging_results/products/OctReimage/'
cubes=glob.glob(home+'*pbcor_line.fits')
sourceregs={'SgrB2S':'fk5; box(266.8353410, -28.3962005, 0.0016806, 0.0016806)','DSi':'fk5; box(266.8316387, -28.3971867, 0.0010556, 0.0010556)','DSii':'fk5; box(266.8335363, -28.3963159, 0.0006389, 0.0006389)','DSiii':'fk5; box(266.8332758, -28.3969270, 0.0006944, 0.0006944)','DSiv':'fk5; box(266.8323834, -28.3954424, 0.0009000, 0.0009000)','DSv':'fk5; box(266.8321331, -28.3976585, 0.0005556, 0.0005556)','DSVI':'fk5; box(266.8380037, -28.4050741, 0.0017361, 0.0017361)','DSVII':'fk5; box(266.8426074, -28.4094401, 0.0020833, 0.0020833)', 'DSVIII':'fk5; box(266.8418408, -28.4118242, 0.0014028, 0.0014028)','DSIX':'fk5; box(266.8477371, -28.4311386, 0.0009583, 0.0009583)','DSX':'fk5; box(266.8452950, -28.4282608, 0.0017083, 0.0017083)','DSXI':'fk5; box(266.8404733, -28.4286378, 0.0013194, 0.0013194)'}
region=sourceregs[source]
outpath_base=f'/blue/adamginsburg/d.jeff/SgrB2DSminicubes/{source}/'
outstatpath_end={1:'OctReimage_K/',10:'field10originals_K/',2:'field2originals_K/',3:'field3originals_K/',7:'field7originals_K/',8:'field8originals_K/'}
outpath=outpath_base+outstatpath_end[fnum]
statfixpath_base='/blue/adamginsburg/d.jeff/SgrB2DSstatcontfix/'
statfixpath=statfixpath_base+outstatpath_end[fnum]

regionparams=[float(val) for val in region[9:(len(region)-1)].split(', ')]

'''This wing of the code runs the linelooper LTE modeling and kinetic temperature determination on the newly created region-specific subcubes'''

print('Begin core cube to Tex map process\n') 

'''Collect constants and CH3OH-specific quantum parameters'''
print('Setting constants')
c=cnst.c*u.m/u.s
k=cnst.k*u.J/u.K
h=cnst.h*u.J*u.s
sigma_sb=cnst.sigma*u.W/((u.m)**(2)*(u.K)**(4))
b_0=rotationalconstants[molecule][0]
a_0=rotationalconstants[molecule][1]
c_0=rotationalconstants[molecule][2]
m=b_0**2/(a_0*c_0)
R_i=1
f=1
Tbg=2.7355*u.K

reorgpath=f'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field{fnum}/CH3OH/{source}/'+ch3oh_sourcedict[source]
mastertxttablepath=reorgpath+'mastereuksqnsfreqsdegens.fits'
fwhmpath=glob.glob(reorgpath+'*fwhm*')[0]
nch3ohpath=reorgpath+'bootstrap_ntot_intstd_boostrap1000_nonegativeslope.fits'
picklepath=reorgpath+f'{molecule}linesdict.obj'
trotmappath=reorgpath+'bootstrap_texmap_3sigma_allspw_withnans_weighted.fits'
fluxerrorpath=reorgpath+'errorimgs/std/*.fits'
ch3oh_trotpath=reorgpath+'bootstrap_texmap_3sigma_allspw_withnans_weighted.fits'
pixcoords=pixdict[source]

c2h5oh_vlsrs={'DSi':(55.3*u.km/u.s),'DSii':(49.3*u.km/u.s)}
#dopplershifts={'SgrB2S':0.00022829138061883716,'DSi':0.0001842772437139578,'DSii':0.00016236367659115043,'DSiii':0.00017500261911843952,'DSiv':0.00018225233186845314,'DSv':0.0001838576164010067,'DSVI':0.0001661613132158407,'DSVII':0.00016320118280935546,'DSVIII':0.0001661546432045067,'DSIX':0.00015453732389175085,'DSX':0.00016375916278648755}#New S z from 5_-0 centroid taken from line spatially integrated over S

vlsr=c2h5oh_vlsrs[source]
z=vlsr/(c.to('km s-1'))
print(f'Doppler shift: {z} / {vlsr}\n')

print('Setting input LTE parameters')
sourcefwhm=fits.getdata(fwhmpath)*u.km/u.s
trotmap=fits.getdata(trotmappath)*u.K
measTrot=trotmap[pixcoords[0],pixcoords[1]]
measlinewidth=sourcefwhm[pixcoords[0],pixcoords[1]]
measTrot=trotmap[pixcoords[0],pixcoords[1]]#{'SgrB2S':215*u.K,'DSi':300*u.K,'DSii':150*u.K,'DSiii':150*u.K,'DSiv':150*u.K,'DSv':150*u.K,'DSVI':300*u.K,'DSVII':200*u.K,'DSVIII':175*u.K,'DSIX':150*u.K}
testT=measTrot#trotdict[source]#500*u.K
#ntotdict=sourcecolumns[source][molecule]#{'SgrB2S':2e14*u.cm**-2,'DSi':1e14*u.cm**-2,'DSii':1e14*u.cm**-2,'DSiii':1e14*u.cm**-2,'DSiv':1e14*u.cm**-2,'DSv':1e14*u.cm**-2,'DSVI':1e14*u.cm**-2,'DSVII':1e14*u.cm**-2,'DSVIII':1e14*u.cm**-2,'DSIX':1e14*u.cm**-2}
testntot=sourcecolumns[source][molecule]
print(f'Input Tex: {testT}\nInput Ntot: {testntot}')
    
'''Gathers beam data from moment map headers'''    
def beamer(momentmap):
    print(f'in beamer function: {momentmap}')
    hdu=fits.getheader(momentmap)
    momentbeam=radio_beam.Beam.from_fits_header(hdu).value
    return momentbeam*u.sr

'''Reorders Splatalogue table parameters to match the glob.glob filename order'''
def unscrambler(filenames,sliced_qns,linelist):
    #print('Start unscrambler')
    unscrambled_qns=[]
    unscrambled_freqs=[]
    unscrambled_eus=[]
    unscrambled_degs=[]
    unscrambled_aijs=[]
    tempfiles=np.copy(filenames)
    for i in range(len(filenames)):
        print(f'filename: {filenames[i]}')
        tempfiles[i]=tempfiles[i].replace('.fits','')
        for j in range(len(sliced_qns)):
            #print(f'sliced_qns: {sliced_qns[j]}')
            #print(f'comparison qns: {tempfiles[i][55:]}')
            comp=(sliced_qns[j]==tempfiles[i][79:])
            if comp==True:
                print(f'{sliced_qns[j]} == {tempfiles[i][79:]}\comp==True')
                unscrambled_qns.append(sliced_qns[j])
                unscrambled_freqs.append(linelist[j])
                unscrambled_eus.append(mastereuks[j]/u.K)
                unscrambled_degs.append(masterdegens[j])
                unscrambled_aijs.append(masteraijs[j])
                break
            else: 
                print(f'{sliced_qns[j]} != {tempfiles[i][79:]}')
    return unscrambled_qns,unscrambled_freqs,unscrambled_eus,unscrambled_degs,unscrambled_aijs

   
def brightnessTandintensities(fluxdict):
    intensitydict={}
    t_bright={}
    dictkeys=fluxdict.keys()
    for key in dictkeys:
        temptransdict=fluxdict[key]
        temptransdictkeys=list(temptransdict.keys())
        print(f'Transition keys in brightnessTandintensities: {temptransdictkeys}')
        for i in range(len(temptransdictkeys)):
            if 'restfreq' in temptransdictkeys[i]:
                continue
            else:
                velflux_T=temptransdict[temptransdictkeys[i]]['flux']
                intensitydict.update({temptransdictkeys[i]:velflux_T})
                temp=velflux_T/measlinewidth#***What's going on here? Should this be fixed because it uses one uniform linewidth
                t_bright.update({temptransdictkeys[i]:temp})
                d_velfluxT=(temptransdict[temptransdictkeys[i]]['stddev'])#/temp)*velflux_T
                intensityerror.append(d_velfluxT)
                #pdb.set_trace()
                
    return intensitydict,t_bright
                
def jupperfinder(quan_nums):
    j_upper=[]
    k_upper=[]
    for i in range(len(quan_nums)):
        for j in range(len(quan_nums[i])):
            comp=quan_nums[i][j].isdigit()
            if comp==False:
                appendage=quan_nums[i][:(j)]
                j_upper.append(int(appendage))
                for k in range(1,len(quan_nums[i][j:])):
                    secondary=quan_nums[i][j+k]
                    if k == 1:
                        if secondary=='-':
                            continue
                    elif secondary.isdigit()==False:
                        appendage=quan_nums[i][(j+1):(j+k)]
                        k_upper.append(int(appendage))
                        break
                break
    return j_upper,k_upper
    
def N_u(nu,Aij,velocityintegrated_intensity_K,velint_intK_err):#(ntot,qrot,gu,eu_J,T_ex):#taken from pyspeckit documentation https://pyspeckit.readthedocs.io/en/latest/lte_molecule_model.html?highlight=Aij#lte-molecule-model
    nuppercalc=((8*np.pi*k*nu**2)/(h*c**3*Aij))*velocityintegrated_intensity_K
    nuppererr=((8*np.pi*k*nu**2)/(h*c**3*Aij))*velint_intK_err#(velint_intK_err.to('K km s-1')/velocityintegrated_intensity_K.to('K km s-1'))*nuppercalc
    return nuppercalc,nuppererr#ntot/(qrot*np.exp(eu_J/(k*T_ex)))
    
def S_j(j_upper,k_upper):#Works for symmetric tops
    return (j_upper**2-k_upper**2)/(j_upper*(2*j_upper+1))
    
'''Loop through a given list of lines (in Hz), computing and saving moment0 maps of the entered data cube'''
def linelooplte(line_list,line_width,iterations,quantum_numbers):
    print('\ncubelooperLTE...')
    print('Grab cube and reference pixel')
    targetspec_K=cube[:,pixycrd,pixxcrd]
    cubebeams=(cube.beams.value)*u.sr/u.beam
    print('Compute cube brightness temperature stddev')
    targetspecK_stddev=stddata[stdpixycrd,stdpixxcrd]
    transitionbeamlist=[]
    transitionfluxlist=[]
    for i in range(iterations):
        print(f'\nStart {quantum_numbers[i]} moment0 procedure')
        temptransdict={}
        line=line_list[i]
        restline=line*(1+z)
        line_sigma=line_width#/(2*np.sqrt(2*np.log(2))),Widening to try to accomodate more of the velocity field/increase the S/N on the eastern side of DS1
        line_sigma_freq=velocitytofreq(line_sigma,line)
        nu_upper=line+line_sigma_freq
        nu_lower=line-line_sigma_freq
        print(f'Make spectral slab between {nu_lower} and {nu_upper}')
        slab=cube.spectral_slab(nu_upper,nu_lower)
        oldstyleslab=cube.spectral_slab((nu_upper-nu_offset),(nu_lower+nu_offset))
        peakchannel=slab.closest_spectral_channel(line)
        print(f'Peak channel: {peakchannel}')
        slabbeams=(slab.beams.value)*u.sr/u.beam
        slab_K=slab[:,pixycrd,pixxcrd]
        mulu2=(mulu(aijs[i],line)).to('cm5 g s-2')
        peak_amplitude=slab_K[peakchannel]#slab_K.max(axis=0)
        if peak_amplitude == slab_K.max(axis=0):
            print('Peak amplitude == Max amplitude in slab')
        else:
            print('Other bright line in slab')
            print(f'Max brightness in slab: {slab_K.max(axis=0)}\n')
        
        phi_nu=lineprofile(sigma=line_sigma_freq,nu_0=restline,nu=restline)
        est_nupper=nupper_estimated(testntot,degeneracies[i],qrot_partfunc,eujs[i],testT).to('cm-2')
        intertau=lte_molecule.line_tau(testT, testntot, qrot_partfunc, degeneracies[i], restline, eujs[i], aijs[i])
        est_tau=(intertau*phi_nu).to('')#opticaldepth(aijs[i],restline,testT,est_nupper,originallinewidth).to('')
        trad=t_rad(tau_nu=est_tau,ff=f,nu=restline,T_ex=testT).to('K')
        
        print('LTE params calculated')
        #print(f'tbthick: {tbthick}\n targetspecK_stddev: {targetspecK_stddev}\n peak_amplitude: {peak_amplitude}')
        print(f'est_nupper: {est_nupper}\n est_tau: {est_tau}\n trad: {trad}')
        #if tbthick >= targetspecK_stddev:
        #    print(f'\n est tau from data at {testT}: {(peak_amplitude/trad).to("")}')
        print('Slicing quantum numbers')
        transition=qn_replace(quantum_numbers[i])
        moment0filename=home+f'{nospace_molecule}~'+transition+'_raw.fits'
        maskedmom0fn=home+f'{nospace_molecule}~'+transition+'_masked.fits'
        maskresidualfn=home+f'{nospace_molecule}~'+transition+'_residual.fits'
        maskedmom0errfn=home+f'{nospace_molecule}~'+transition+'_error.fits'
        slabfilename=slabpath+f'{nospace_molecule}~'+transition+'_slab.fits'
        maskedslabfn=slabpath+f'{nospace_molecule}~'+transition+'_maskedslab.fits'
        maskfn=slabpath+f'{nospace_molecule}~'+transition+'_mask.fits'
        peakintfn=home+f'{nospace_molecule}~'+transition+'_peakint.fits'
        #print('Done')
        #print('Moment 0')
        if os.path.isfile(maskedmom0fn):
            print(f'{moment0filename} already exists.')
            isfilemom0=fits.getdata(maskedmom0fn)*u.K*u.km/u.s
            #isfilepixflux=isfilemom0[pixycrd,pixxcrd]
            isfilebeam=beamer(maskedmom0fn)
            isfilestdflux=stddata#fits.getdata(f'{stdpath}{images[imgnum]}fluxstd.fits')*u.K#This is confusing, notation-wise, but I'm leaving it this way for now since it's consistent between the two forks in the loop. For future reference: isfilestdflux is the error on the measured brightnesstemp in K, whereas isfilemom0 pulls from the moment0 maps and is in K km/s
            isfilemom0err=fits.getdata(maskedmom0errfn)*u.K*u.km/u.s
            temptransdict.update([('freq',restline),('flux',isfilemom0),('stddev',isfilestdflux),('beam',isfilebeam),('euk',euks[i]),('eujs',eujs[i]),('degen',degeneracies[i]),('aij',aijs[i]),('filename',moment0filename),('shift_freq',line),('mom0err',isfilemom0err)])
            transitiondict.update({transition:temptransdict})
            masterslicedqns.append(transition)
            mastereuks.append(euks[i].value)
            masterstddevs.append(targetspecK_stddev)
            masterqns.append(quantum_numbers[i])
            masterlines.append(line_list[i].value)
            print('\nDictionaries populated for this transition.')
            if os.path.isfile(maskedslabfn):
                print('Masked slab already exists...\n')
                pass
            else:
                slab.write(maskedslabfn)
                print(f'Slab written to {slabfilename}. Proceeding...\n')
            if os.path.isfile(peakintfn):
                print('Peak intensity file already exists.')
                pass
            else:
                print('Peak intensity procedure starting')
                print('Creating masked slab')
                slab=slab.with_spectral_unit((u.km/u.s),velocity_convention='radio',rest_value=lines[i])
                slab3sigmamask=slab > (3*stdcutout.data)
                slab=slab.with_mask(slab3sigmamask)
                slabspecax=slab.spectral_axis
                slabmom1=slab.moment1()
                slabfwhm=slab.linewidth_fwhm()#(7*u.MHz/line)*c.to('km s-1')#
                cubemask=(slabspecax[:,None,None] < (velocityfield_representative + fwhm_representative)[None,:,:]) & (slabspecax[:,None,None] > (velocityfield_representative - fwhm_representative)[None,:,:])
                maskedslab=slab.with_mask(cubemask)
                print('Computing peak intensity image')
                maskedpeakint=maskedslab.max(axis=0)
                print(f'Saving to {peakintfn}')
                maskedpeakint.write(peakintfn)
                print('Peak intensity image saved.\n')
            for moment in [1,2]:
                slab=slab.with_spectral_unit((u.km/u.s),velocity_convention='radio',rest_value=lines[i])
                momentnfilenames=[(sourcepath+f'mom1/{nospace_molecule}~'+transition+'.fits'),(sourcepath+f'mom2/{nospace_molecule}~'+transition+'_var.fits')]
                fwhmfilename=sourcepath+f'mom2/{nospace_molecule}~'+transition+'_fwhm.fits'
                if os.path.isfile(momentnfilenames[moment-1]):
                    print(f'{transition} moment{moment} file already exists.\n')
                    continue
                elif moment == 1:
                    print(f'Computing moment 1 and saving to {momentnfilenames[moment-1]}\n')
                    slabmom1=slab.moment1()
                    slabmom1.write(momentnfilenames[moment-1])
                elif moment == 2:
                    print(f'Computing moment 2 and saving to {momentnfilenames[moment-1]}\n')
                    slabmom2=slab.moment2()
                    slabfwhm=slab.linewidth_fwhm()
                    slabmom2.write(momentnfilenames[moment-1])
                    slabfwhm.write(fwhmfilename)
            pass
        elif trad >= 3*targetspecK_stddev and peak_amplitude >= 3* targetspecK_stddev:#*u.K:
            slab=slab.with_spectral_unit((u.km/u.s),velocity_convention='radio',rest_value=lines[i])
            
            vel_lineprofilesigma=measlinewidth/(2*np.sqrt(2*np.log(2)))
            modelline=models.Gaussian1D(mean=(0*u.km/u.s), stddev=vel_lineprofilesigma, amplitude=trad)
            search_for_fwhm=np.abs(np.ones(len(slab.spectral_axis.value))-np.abs(slab.spectral_axis.value/vel_lineprofilesigma).value)
            
            fwhm_channel=np.where(search_for_fwhm == min(search_for_fwhm))[0][0]#Find where difference is minimized and set as fwhm channel
            central_channel=np.where(slab.spectral_axis.value == min(slab.spectral_axis.value))[0][0]#Find central channel of slab
            channelrange_fwhm=np.abs(central_channel-fwhm_channel)#Compute channel range for chi squared
            
            slabcutout_forchisquare=slab_K[(central_channel-channelrange_fwhm):(central_channel+channelrange_fwhm)]#Grabs the center channel and the fwhm range
            velocityrange_chisquared=slab.spectral_axis[(central_channel-channelrange_fwhm):(central_channel+channelrange_fwhm)]
            
            length_chisquareslab=len(slabcutout_forchisquare)
            #velocityrange_chisquared=np.linspace(-(measlinewidth/2),(measlinewidth/2),length_slab)#Set FWHM velocity range over which chi-squared will be computed
            chisquared=np.sum((slabcutout_forchisquare-modelline(velocityrange_chisquared))**2/(targetspecK_stddev*np.ones(length_chisquareslab))).value/length_chisquareslab
            degrees_of_freedom=length_chisquareslab-2#we would be fitting Trot and Ntot
            goodness_of_fit=np.abs(np.sqrt(2*chisquared)-np.sqrt(2*degrees_of_freedom-1))
            #pdb.set_trace()
            if transition in excludedlines[source] or goodness_of_fit >= 3:
                print(f'\nExcluded line detected: {quantum_numbers[i]}, E_U: {euks[i]}, Freq: {line.to("GHz")}')
                print(f'Goodness of fit: {goodness_of_fit}')
                sigma1mom0=stdcutout.data*measlinewidth
                sigma1hdu=fits.PrimaryHDU(sigma1mom0.value)
                sigma1hdu.header=fits.open(rep_mom1)[0].header
                sigma1hdu.header['RESTFRQ']=line.value
                sigma1hdul=fits.HDUList([sigma1hdu])
                print(f'1sigma moment0 of shape ({cube.shape[1]},{cube.shape[2]}) populated\nTarget pix flux value: {stddata[pixycrd,pixxcrd]*measlinewidth}')
                moment0beam=slab.beams
                maskslabmom0=sigma1mom0
                #maskslabmom0.write((maskedmom0fn))
                #pdb.set_trace()
                sigma1hdul.writeto(moment0filename,overwrite=True)
                print(f'Saved to {moment0filename}')
                kkmsstdarray=maskslabmom0
                pass
            else:
                print('Commence moment0 procedure\n')
                #cubemask=BooleanArrayMask(mask=cubemaskarray,wcs=slab.wcs)
                print(f'Create {quantum_numbers[i]} spatial-velocity mask')
                slab3sigmamask=slab > (3*stdcutout.data)
                slab=slab.with_mask(slab3sigmamask)
                slabspecax=slab.spectral_axis
                slabmom1=slab.moment1()
                slabfwhm=slab.linewidth_fwhm()#(7*u.MHz/line)*c.to('km s-1')#
                cubemask=(slabspecax[:,None,None] < (velocityfield_representative + fwhm_representative)[None,:,:]) & (slabspecax[:,None,None] > (velocityfield_representative - fwhm_representative)[None,:,:])
                oldstyleslab=oldstyleslab.with_spectral_unit((u.km/u.s),velocity_convention='radio',rest_value=lines[i])
                #if imgnum > 0:
                #    pdb.set_trace()
                print('Masking spectral slab')
                maskedslab=slab.with_mask(cubemask)
                #momstart=time.time()
                print('Unmasked moment0 computing...\n')
                slabmom0=oldstyleslab.moment0()
                print('Masked moment0 computing...\n')
                maskslabmom0=maskedslab.moment0()
                #momend=time.time()-momstart
                #print(f'{quantum_numbers[i]} elapsed time: {time.strftime("%H:%M:%S", time.gmtime(momend))}')
                print('\nComputing masking residuals')
                mom0maskresiduals=maskslabmom0-slabmom0
                print('\nSaving...')
                #name='test'+str(i)
                slabmom0.write((moment0filename),overwrite=True)
                maskslabmom0.write((maskedmom0fn))
                mom0maskresiduals.write((maskresidualfn))
                #maskboolarr=BooleanArrayMask(mask=cubemask,wcs=slab.wcs)
                #make_casa_mask(maskedslab,maskfn,append_to_image=False,add_stokes=False)
                moment0beam=slabmom0.beam.value*u.sr
                kkmsstdarray=stdcutout.data*fwhm_representative
                if os.path.isfile(slabfilename):
                    print(f'Spectral slab {slabfilename} already exists.\nProceeding...\n')
                    pass
                else:
                    slab.write(slabfilename)
                    print(f'Slab written to {slabfilename}.')
                    maskedslab.write(maskedslabfn)
                    print(f'Masked slab written to {maskedslabfn}.\n')
                    print('Creating moment0 error Primary HDU')
                    kkmshdu=fits.PrimaryHDU(kkmsstdarray.value)
                    kkmshdu.header=maskedslab.header
                    kkmshdu.header['BUNIT']='K km s-1'
                    print('Wrapping moment0 error Primary HDU in HDUList')
                    kkmshdul=fits.HDUList([kkmshdu])
                    print(f'Writing moment0 error to {maskedmom0errfn}')
                    kkmshdul.writeto(maskedmom0errfn,overwrite=True)
                if os.path.isfile(peakintfn):
                    print('Peak intensity file already exists.')
                    pass
                else:
                    print('Peak intensity procedure starting')
                    print('Computing peak intensity image')
                    maskedpeakint=maskedslab.max(axis=0)
                    print(f'Saving to {peakintfn}')
                    maskedpeakint.write(peakintfn)
                    print('Peak intensity image saved.\n')
                for moment in [1,2]:
                    moment1filename=sourcepath+f'mom1/{nospace_molecule}~'+transition+'.fits'
                    moment2filename=sourcepath+f'mom2/{nospace_molecule}~'+transition+'_var.fits'
                    fwhmfilename=sourcepath+f'mom2/{nospace_molecule}~'+transition+'_fwhm.fits'
                    if moment == 1:
                        print(f'Computing moment 1 and saving to {moment1filename}\n')
                        #slabmom1=slab.moment1()
                        slabmom1.write(moment1filename)
                    elif moment == 2:
                        print(f'Computing moment 2s and saving to {moment2filename} and {fwhmfilename}\n')
                        slabmom2=slab.moment2()
                        slabfwhm=slab.linewidth_fwhm()
                        slabmom2.write(moment2filename)
                        slabfwhm.write(fwhmfilename)
                pass
            temptransdict.update([('freq',restline),('flux',maskslabmom0),('stddev',targetspecK_stddev),('beam',moment0beam),('euk',euks[i]),('eujs',eujs[i]),('degen',degeneracies[i]),('aij',aijs[i]),('filename',moment0filename),('shift_freq',line),('mom0err',kkmsstdarray)])
            transitiondict.update({transition:temptransdict})
            mastereuks.append(euks[i].value)
            masterstddevs.append(targetspecK_stddev)
            masterqns.append(quantum_numbers[i])
            masterlines.append(line_list[i].value)
            print(f'{quantum_numbers[i]} calculations complete.\n')
            pass
        else:
            if not trad >= 3*targetspecK_stddev and peak_amplitude >= 3* targetspecK_stddev:
                print('LTE Model max brightnessTemp below 3sigma threshold')
                print(f'{quantum_numbers[i]} skipped, possible contamination\n')
                pass
            elif trad >= 3*targetspecK_stddev and not peak_amplitude >= 3* targetspecK_stddev:
                print(f'Line amplitude ({peak_amplitude}) less than 3 sigma criterion ({3*targetspecK_stddev})')
                print(f'{quantum_numbers[i]} skipped\n')
            elif not trad >= 3*targetspecK_stddev and not peak_amplitude >= 3* targetspecK_stddev:
                print('3 sigma LTE model and 3 sigma amplitude criteria not met')
                print(f'{quantum_numbers[i]} skipped\n')
                pass
    spectraKdict.update({images[imgnum]:targetspec_K})
    print('lines looped.\n')
    
#Check what linelist molecule uses and pull molecular parameters from linelist
linelist=linelistdict[molecule]
if linelist=='CDMS':
        if molecule in cdmsnamelist.keys():
            queryname=cdmsnamelist[molecule]
        elif molecule == 'CH3NCO, vb = 0':
            queryname=cdms_get_molecule_name('CH3NCO, vb=0')
        else:
            param_name=molecule.replace(' ','')
            queryname=cdms_get_molecule_name(param_name)#cdmsnamelist[molecule]

        if queryname in cdmsproblemchildren:
            freqs, aij, deg, EU, qrot = get_molecular_parameters(queryname,catalog='CDMS',
                                                                           fmin=216*u.GHz, fmax=(218*u.GHz+100*u.GHz),)
        else:
            freqs, aij, deg, EU, qrot = get_molecular_parameters(queryname,catalog='CDMS',fmin=216*u.GHz,fmax=(218*u.GHz+50*u.GHz))
if linelist == 'JPL':
    if molecule in jplnamelist.keys():
        queryname=jplnamelist[molecule]
    else:
        queryname=molecule.replace(' ','')
    freqs, aij, deg, EU, qrot = get_molecular_parameters(queryname,catalog='JPL',
                                                 fmin=216*u.GHz,
                                                 fmax=234*u.GHz,)

#Check if molecule has partition function out to hot core temperatures and, if it doesn't, extrapolate using line
if molecule in incompleteqrot:
    print(f'{molecule} has an incomplete partition function')
    print('Estimating by linear fit to log-log Qrot/T relation')
    poly=Linear1D(slope=150, intercept=10)
    fitter=fitting.LinearLSQFitter()
    fitinput_xvalues=np.linspace(3,300,1000)*u.K
    power_law_fit=fitter(poly,np.log10(fitinput_xvalues.value),np.log10(qrot(fitinput_xvalues)))
    logintercept=10**power_law_fit.intercept
    logTs=logintercept*fitinput_xvalues.value**power_law_fit.slope

qrot_partfunc=fit_qrot(logintercept,measTrot,power_law_fit)

incubes=glob.glob(outpath+"*pbcor_line.fits")#'/blue/adamginsburg/d.jeff/imaging_results/field1core1box2/*.fits')

images=['spw0','spw1','spw2','spw3']

datacubes=[]

for spew in images:
    for f1 in incubes:
        if spew in f1:
            datacubes.append(f1)
            continue
    
assert 'spw0' in datacubes[0], 'Cube list out of order'

stdhomedict={1:'/orange/adamginsburg/sgrb2/d.jeff/products/OctReimage_K/',10:'/orange/adamginsburg/sgrb2/d.jeff/products/field10originals_K/',2:'/orange/adamginsburg/sgrb2/d.jeff/products/field2originals_K/',3:'/orange/adamginsburg/sgrb2/d.jeff/products/field3originals_K/',7:'/orange/adamginsburg/sgrb2/d.jeff/products/field7originals_K/'}
stdhome=stdhomedict[fnum]

#cubemaskarray=maskeddatacube.get_mask_array()
c2h5oh_sourcelocs={'DSi':'/chisquare_goodnessoffit_3_2of3contamsremoved/'}

sourcelocs=c2h5oh_sourcelocs

representativelines={'DSi':'14_0&14-13_1&13'}
representativelws=measlinewidth#{'SgrB2S':(5*u.km/u.s),'DSi':(3*u.km/u.s),'DSii':(3*u.km/u.s),'DSiii':(3*u.km/u.s),'DSiv':(4*u.km/u.s),'DSv':(4*u.km/u.s),'DSVI':(3*u.km/u.s),'DSVII':(2.5*u.km/u.s),'DSVIII':(2.5*u.km/u.s),'DSIX':(5*u.km/u.s),'DSX':(4*u.km/u.s)}#{'SgrB2S':8*u.MHz,'DSi':3.6*u.MHz}#11MHz for ~10 km/s
representativecubes={'SgrB2S':2,'DSi':2,'DSii':1,'DSiii':2,'DSiv':0,'DSv':1,'DSVI':1,'DSVII':1,'DSVIII':1,'DSIX':1,'DSX':1}#spwnumber

sourcepath=f'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field{fnum}/{nospace_molecule}/{source}/'+sourcelocs[source]
nupperpath=sourcepath+'nuppers/'
stdpath=sourcepath+'errorimgs/std/'
slabpath=sourcepath+'spectralslabs/km_s/'
mom0path=sourcepath+'mom0/'
rotdiagpath=sourcepath+'pixelwiserotationaldiagrams/'
figpath=sourcepath+'figures/'

overleafpath="/blue/adamginsburg/d.jeff/repos/C2H5OH_SgrB2DS/figures/"

picklepath=sourcepath+f'{nospace_molecule}linesdict.obj'

if os.path.isdir(slabpath):
    print(f'Source path directory tree {sourcepath} already exists.\n')
    if os.path.isdir(sourcepath+'mom1/'):
        print('Moment 1/2 directories already exist.')
    else:
        for moment in [1,2]:
            momnpath=sourcepath+f'mom{moment}/'
            print(f'Creating moment {moment} directory at {momnpath}')
            os.mkdir(momnpath)
else:
    print(f'Making source path {sourcepath}')
    os.makedirs(sourcepath)
    print(f'Making nupper folder {nupperpath}')
    os.mkdir(nupperpath)
    print(f'Making error folder {stdpath}')
    os.makedirs(stdpath)
    print(f'Making spectral slab folder {slabpath}\n')
    os.makedirs(slabpath)
    for moment in [0,1,2]:
            momnpath=sourcepath+f'mom{moment}/'
            print(f'Creating moment {moment} directory at {momnpath}')
            os.mkdir(momnpath)
    #print(f'Making mom0 folder {mom0path}')
    #os.mkdir(mom0path)
    print(f'Making rotational diagram folder')
    os.mkdir(rotdiagpath)
    print(f'Making figures folder')
    os.mkdir(figpath)

spwdict={}
kstddict={}
kkmsstddict={}
spectraKdict={}

masterlines=[]
masterqns=[]
mastereuks=[]
mastereujs=[]
masterdegens=[]
masterlog10aijs=[]
masteraijs=[]
masterslicedqns=[]
masterrestfreqs=[]

masterfluxes=[]
masterbeams=[]
masterstddevs=[]

catdir=CDMS.get_species_table()
catdir_c2h5oh=catdir[catdir['TAG'] == 46524]
catdir_qrot300=10**catdir_c2h5oh['lg(Q(300))']

excludedlines={'DSi':['30_3&27-30_2&28&anti','19_5&15-19_4&16&anti']}
restfreq_representativeline={'SgrB2S':'','DSi':230.9913834*u.GHz,'DSii':'','DSiii':'','DSiv':'','DSv':'','DSVI':'','DSVII':'','DSVIII':'','DSIX':'','DSX':''}#All taken from Splatalogue
representative_filename_base=sourcepath+representativelines[source]+'repline_'
rep_mom1=representative_filename_base+'mom1.fits'
rep_fwhm=representative_filename_base+'fwhm.fits'
rep_slab=representative_filename_base+'slab.fits'
rep_maskedslab=representative_filename_base+'maskedslab.fits'

regiondims=regionparams[3]

if os.path.isfile(rep_mom1):
    print(f'{source} representative line objects already exist.')
    print(f'Grabbing spectralslab from {rep_slab}')
    spectralslab_rep_3sigma=fits.getdata(rep_maskedslab)*u.K
    print(f'Grabbing mom1 from {rep_mom1}')
    velocityfield_representative=fits.getdata(rep_mom1)*u.km/u.s
    print(f'Grabbing fwhm from {rep_fwhm}\n')
    fwhm_representative=fits.getdata(rep_fwhm)*u.km/u.s
else:
    print(f'{representativelines[source]} representative line objects will be computed for {source}.\n')
    pass

for imgnum in range(len(datacubes)):
    print(f'Accessing data cube {datacubes[imgnum]}')
    assert images[imgnum] in datacubes[imgnum], f'{images[imgnum]} not in filename {datacubes[imgnum]}'
    home=sourcepath+'mom0/'#f'{images[imgnum]}/'#Make sure to include slash after path
    readstart=time.time()
    cube=sc.read(datacubes[imgnum])
    
    header=fits.getheader(datacubes[imgnum])
    
    stdimage=fits.open(stdhome+images[imgnum]+'minimize.image.pbcor_noise.fits')
    stdcellsize=(np.abs(stdimage[0].header['CDELT1']*u.deg)).to('arcsec')
    stdcutoutsize=round(((float(regiondims)*u.deg)/stdcellsize).to('').value)
    stddata=stdimage[0].data*u.K
    
    print('Acquiring cube rest frequency and computing target pixel coordinates')
    spwrestfreq=header['RESTFRQ']*u.Hz
    masterrestfreqs.append(spwrestfreq)
    
    freqs=cube.spectral_axis#Hz
    freqflip=False
    if freqs[1] < freqs[0]:
        freqs=freqs[::-1]
        freqflip=True
        print('Corrected decreasing frequency axis')
    else:
        pass
        
    velcube=cube.with_spectral_unit((u.km/u.s),velocity_convention='radio',rest_value=spwrestfreq)
    cube_unmasked=velcube.unmasked_data
    
    targetworldcrds={'SgrB2S':[[0,0,0],[2.66835339e+02, -2.83961660e+01, 0]], 'DSi':[[0,0,0],[266.8316149,-28.3972040,0]], 'DSii':[[0,0,0],[266.8335363,-28.3963158,0]],'DSiii':[[0,0,0],[266.8332758,-28.3969269,0]],'DSiv':[[0,0,0],[266.8323834, -28.3954424,0]],'DSv':[[0,0,0],[266.8321331, -28.3976585, 0]],'DSVI':[[0,0,0],[266.8380037, -28.4050741,0]],'DSVII':[[0,0,0],[266.8426074, -28.4094401,0]],'DSVIII':[[0,0,0],[266.8418408, -28.4118242, 0]],'DSIX':[[0,0,0],[266.8477371, -28.4311386,0]],'DSX':[[0,0,0],[266.8452950, -28.4282608,0]]}
    cube_w=cube.wcs
    stdwcs=WCS(stdimage[0].header)
    
    targetworldcrd=targetworldcrds[source]
    targetpixcrd=cube_w.all_world2pix(targetworldcrd,1,ra_dec_order=True)
    fullsize_targetpixcrd=stdwcs.wcs_world2pix(targetworldcrd,1,ra_dec_order=True)
    stdpixxcrd,stdpixycrd=int(round(fullsize_targetpixcrd[1][0])),int(round(fullsize_targetpixcrd[1][1]))
    print(f'Stddev position - x: {stdpixxcrd}/y: {stdpixycrd}')
    
    assert stdpixxcrd >= 0 and stdpixycrd >= 0, 'Negative std pixel coords'
    
    pixxcrd,pixycrd=int(round(targetpixcrd[1][0])),int(round(targetpixcrd[1][1]))
    print(f'Flux position - x: {pixxcrd}/y: {pixycrd}')
    
    assert pixxcrd >= 0 and pixycrd >= 0, 'Negative pixel coords'
    
    stdcutout=Cutout2D(stddata,(stdpixxcrd,stdpixycrd),(stdcutoutsize))
    assert np.shape(stdcutout)[0]==cube.shape[1], 'Standard deviation cutout size mismatch'
    
    '''Create representative line to use as masking template'''
    if not os.path.isfile(rep_mom1):
        print(f'Creating {representativelines[source]} representative line data objects for {source}')
        print(f'Opening file {datacubes[representativecubes[source]]}')
        reffreq_repline=restfreq_representativeline[source]/(1+z)
        repcube=sc.read(datacubes[representativecubes[source]])
        
        assert reffreq_repline > min(repcube.spectral_axis) and reffreq_repline < max(repcube.spectral_axis), f'Incorrect representative cube chosen.\nCube limits: {min(repcube.spectral_axis)},{max(repcube.spectral_axis)}\nRepresentative line reference frequency: {reffreq_repline}'
        
        repstdmain=fits.open(stdhome+f'spw{representativecubes[source]}minimize.image.pbcor_noise.fits')
        repstdmain_data=repstdmain[0].data*u.K
        repstdmain_wcs=WCS(repstdmain[0].header)
        repstdmain_cellsize=(np.abs(repstdmain[0].header['CDELT1']*u.deg)).to('arcsec')
        repstdmain_pixcrd=repstdmain_wcs.wcs_world2pix(targetworldcrd,1,ra_dec_order=True)
        repstdmain_xcrd,repstdmain_ycrd=int(round(repstdmain_pixcrd[1][0])),int(round(repstdmain_pixcrd[1][1]))
        
        assert repstdmain_xcrd and repstdmain_ycrd > 0, 'Negative representative std pixel coords'
        
        repstdcutoutsize=round(((float(regiondims)*u.deg)/stdcellsize).to('').value)
        repstdcutout=Cutout2D(repstdmain_data,(repstdmain_xcrd,repstdmain_ycrd), repstdcutoutsize)
        upperfreq=reffreq_repline+velocitytofreq(representativelws,reffreq_repline)
        lowerfreq=reffreq_repline-velocitytofreq(representativelws,reffreq_repline)
        print(f'Creating spectral slab between {lowerfreq} and {upperfreq}')
        spectralslab_representative=repcube.spectral_slab(lowerfreq,upperfreq)
        spectralslab_representative=spectralslab_representative.with_spectral_unit((u.km/u.s),velocity_convention='radio',rest_value=reffreq_repline)
        rep_3sigmamask=spectralslab_representative >= (3*repstdcutout.data)
        spectralslab_rep_3sigma=spectralslab_representative.with_mask(rep_3sigmamask)
        
        print('\nComputing moment1')
        velocityfield_representative=spectralslab_rep_3sigma.moment1()
        print('\nComputing fwhm')
        fwhm_representative=spectralslab_rep_3sigma.linewidth_fwhm()
        #pdb.set_trace()
        print('Writing objects to file')
        spectralslab_representative.write(rep_slab)
        spectralslab_rep_3sigma.write(rep_maskedslab)
        velocityfield_representative.write(rep_mom1)
        fwhm_representative.write(rep_fwhm)
        print('Continuing to line loops\n')
    else:
        pass

    #pdb.set_trace()
    
    freq_min=freqs[0]*(1+z)#215*u.GHz
    #print(freq_max)
    freq_max=freqs[(len(freqs)-1)]*(1+z)#235*u.GHz
    
    assert freq_max > freq_min, 'Inverted spectral axis'
    print('Passed increasing spectral axis check')
    #print(freq_min)
    linelist=['CDMS']
    
    print('Peforming Splatalogue queries')
    maintable = CDMS.query_lines(min_frequency=freq_min,max_frequency=freq_max,min_strength=-500,molecule='046524 C2H5OH,v=0',get_query_payload=False)#utils.minimize_table(Splatalogue.query_lines(freq_min, freq_max, chemical_name=molecule,
                                    #energy_max=1840, energy_type='eu_k',
                                    #line_lists=linelist,
                                    #show_upper_degeneracy=True))
    '''Needed for upper state degeneracies'''                                
    sparetable=CDMS.query_lines(min_frequency=freq_min,max_frequency=freq_max,min_strength=-500,molecule='046524 C2H5OH,v=0',get_query_payload=False)#Splatalogue.query_lines(freq_min, freq_max, chemical_name=molecule,
                                    #energy_max=1840, energy_type='eu_k',
                                    #line_lists=linelist,
                                    #show_upper_degeneracy=True)
                                    
    
    print('Gathering CDMS table parameters')
    nus=maintable['FREQ']
    lines=nus/(1+z)#Redshifted to source
    
    elo_lambda=(1/maintable['ELO'].data)*u.cm
    elo_K=(((h*c)/elo_lambda)/k).to('K')
    elo_J=(elo_K*k).to('J')
    deltae=((h*maintable['FREQ'])/k).to('K')
    euks=elo_K+deltae#maintable['EU_K']*u.K
    eujs=(euks*k).to('J')
    degeneracies=maintable['GUP']
    
    ju=maintable['Ju']
    jl=maintable['Jl']
    ku1=maintable['Ku']
    ku2=maintable['vu']
    kl1=maintable['Kl']
    kl2=maintable['vl']
    qns=[]
    assert len(ju)==len(maintable) and len(jl)==len(maintable)
    for jupper,jlower,kupper1,kupper2,klower1,klower2 in zip(ju,jl,ku1,ku2,kl1,kl2):
        tempqn=f'{jupper}({kupper1},{kupper2})-{jlower}({klower1},{klower2})'
        qns.append(tempqn)
    
    log10cdmsfluxes=maintable['LGINT']
    cdmsfluxes=10**log10cdmsfluxes
    aijs=pickett_aul(cdmsfluxes,nus,degeneracies,elo_J,eujs,catdir_qrot300,T=300*u.K)
    
    #pdb.set_trace()
    
    singlecmpntwidth=velocitytofreq(measlinewidth,spwrestfreq).to('GHz')#(0.00485/8)*u.GHz
    linewidth=representativelws#10*u.km/u.s#8*u.MHz
    linewidth_freq=velocitytofreq(linewidth,restfreq_representativeline[source])
    oldwideslabwidth=(15.15*u.MHz)
    originallinewidth=(11231152.36688232*u.Hz/2)#0.005*u.GHz#e###0.5*0.0097*u.GHz#from small line @ 219.9808GHz# 0.0155>>20.08km/s 
    nu_offset=oldwideslabwidth-originallinewidth
    #linewidth_vel=measlinewidth
    
    pixeldict={}
    transitiondict={}
    linelooplte(lines,linewidth,len(lines),qns)
    spwdict.update([(images[imgnum],transitiondict)])
    tempkeys=list(spwdict[images[imgnum]].keys())
    
    kstdimgpath=stdpath+f'{images[imgnum]}fluxstd.fits'
    if os.path.isfile(kstdimgpath):
        print(f'{images[imgnum]} brightness std image already exists')
        spwstdarray=fits.getdata(kstdimgpath)*u.K#stdcutout.data
        #kkmsstdarray=fits.getdata(kkmsstdimgpath)*u.K*u.km/u.s
        #print(f'Retrieved integrated intensity std data from {kkmsstdimgpath}\n')
    else:
        print(f'Start {images[imgnum]} std calculations')
        spwstdarray=stdcutout.data
        print('Set Primary HDU')
        khdu=fits.PrimaryHDU(spwstdarray.value)
        '''This transmoment0 file has intensity (K km/s) units'''
        if len(tempkeys) == 0:
            print(f'No transitions detected in this spw ({images[imgnum]})')
            transmomslab=cube.spectral_slab((lines[0]-linewidth_freq),(lines[0]+linewidth_freq))
            transmoment0=transmomslab.moment0()
            transmom0header=transmoment0.header
            print(f'Set transmoment0 to moment0 from {(lines[0]+linewidth_freq).to("GHz")} to {(lines[0]-linewidth_freq).to("GHz")}')
        else:
            transmoment0=fits.open(spwdict[images[imgnum]][tempkeys[0]]['filename'])
            transmom0header=transmoment0[0].header
            print(f'Set header from {spwdict[images[imgnum]][tempkeys[0]]["filename"]}')
        khdu.header=transmom0header
        khdu.header['BUNIT']='K'
        print('Wrapping Primary HDU in HDUList')
        khdul=fits.HDUList([khdu])
        print(f'Writing to {kstdimgpath}')
        khdul.writeto(kstdimgpath,overwrite=True)
        print(f'{images[imgnum]} std calculations complete.\n')
    #transitiondict.update({'restfreq':spwrestfreq})
    #,('pixel_0',(pixycrd,pixxcrd))])
    kstddict.update([(images[imgnum],spwstdarray)])
    #kkmsstddict.update([(images[imgnum],kkmsstdarray)])
    print(f'Finished loop for {images[imgnum]}\n')
    #masterqns.append(slicedqns)
    #pdb.set_trace()

if os.path.isfile(picklepath):
    print(f'pickle {picklepath} already exists.')
else:
    print('Saving dictionary pickle...')
    f=open(picklepath,'wb')
    pickle.dump(spwdict,f)
    f.close()
    print(f'Dictionary pickle saved at {picklepath}')

print('Computing K km/s intensities and K brightness temperatures')
intensityerror=[]
intensities,t_brights=brightnessTandintensities(spwdict)

print(intensityerror)

print('Begin fitting procedure\nCompute N_uppers')
spwdictkeys=spwdict.keys()
print(f'spwdictkeys: {spwdictkeys}')
testyshape=np.shape(cube)[1]#60
testxshape=np.shape(cube)[2]#60
testzshape=len(mastereuks)
nugsmap=np.empty(shape=(testyshape,testxshape,testzshape))
nugserrormap=np.empty(shape=(testyshape,testxshape,testzshape))
orderedeuks=[]
ordereddegens=[]
print(f'Begin pixel loops of shape ({testyshape},{testxshape})')
pixelzcoord_nupper=0
pixelzcoord_nuperr=0
master_transkeys=[]
spws_with_detections=[]
for key in spwdictkeys:
    transdict=spwdict[key]
    #print(f'transdict: {transdict}')
    transitionkeys=list(spwdict[key])
    master_transkeys.append(transitionkeys)
    if len(transitionkeys) > 0:
        spws_with_detections.append(transdict)
    else:
        pass
    #print(f'transitionkeys: {transitionkeys}')
    for transkey in range(len(transitionkeys)):#Need to figure out way to store the n_us per pixel, per moment map. possibly append in 3D array
        print(f'Transition: {transitionkeys[transkey]}/Nupper array z-coord: {pixelzcoord_nupper}')
        nupperimage_filepath=nupperpath+f'{nospace_molecule}~'+transitionkeys[transkey]+'data.fits'
        nuperrorimage_filepath=nupperpath+f'{nospace_molecule}~'+transitionkeys[transkey]+'error.fits'
        orderedeuks.append(transdict[transitionkeys[transkey]]['euk'])
        ordereddegens.append(transdict[transitionkeys[transkey]]['degen'])
        
        nupperimgexists=False
        nuperrorimgexists=False
        if os.path.isfile(nupperimage_filepath):
            print(f'{nupperimage_filepath} already exists.\nPopulating nuppers array...\n')
            tempnupper=fits.getdata(nupperimage_filepath)
            nupperimgexists=True
            nugsmap[:,:,pixelzcoord_nupper]=tempnupper
            pixelzcoord_nupper+=1
        if os.path.isfile(nuperrorimage_filepath):
            print(f'{nuperrorimage_filepath} already exists\nPopulating nupper error array...\n')
            tempnuerr=fits.getdata(nuperrorimage_filepath)
            nuperrorimgexists=True
            nugserrormap[:,:,pixelzcoord_nuperr]=tempnuerr
            pixelzcoord_nuperr+=1
        elif not nupperimgexists or not nuperrorimgexists:
            for y in range(testyshape):
                print(f'Row {y} Looping')
                n_us=[]#np.empty(np.shape(intensityerror))
                n_uerr=[]#np.empty(np.shape(intensityerror))
                for x in range(testxshape):
                    if nupperimgexists:
                        n_us.append((tempnupper[y,x])/transdict[transitionkeys[transkey]]['degen'])
                    if nuperrorimgexists:
                        n_uerr.append((tempnuerr[y,x])/transdict[transitionkeys[transkey]]['degen'])
                    else:
                        tempnupper,tempnuerr=N_u(transdict[transitionkeys[transkey]]['freq'],transdict[transitionkeys[transkey]]['aij'],intensities[transitionkeys[transkey]][y,x],transdict[transitionkeys[transkey]]['mom0err'][y,x])#kkmsstddict[key][y,x])
                        n_us.append((tempnupper.to('cm-2')*u.cm**2)/transdict[transitionkeys[transkey]]['degen'])
                        n_uerr.append((tempnuerr.to('cm-2')*u.cm**2)/transdict[transitionkeys[transkey]]['degen'])  
                nugsmap[y,:,pixelzcoord_nupper]=n_us
                nugserrormap[y,:,pixelzcoord_nuperr]=n_uerr
            if not nupperimgexists:
                nupperimgdata=nugsmap[:,:,pixelzcoord_nupper]
                primaryhdu=fits.PrimaryHDU(nupperimgdata)
                transmoment0=fits.open(transdict[transitionkeys[transkey]]['filename'])
                transmom0header=transmoment0[0].header
                primaryhdu.header=transmom0header
                primaryhdu.header['BTYPE']='Upper-state column density'
                primaryhdu.header['BUNIT']='cm-2'
                hdul=fits.HDUList([primaryhdu])
                hdul.writeto(nupperimage_filepath,overwrite=True)
            if not nuperrorimgexists:
                nuperrorimgdata=nugserrormap[:,:,pixelzcoord_nuperr]
                primaryhduerr=fits.PrimaryHDU(nuperrorimgdata)
                transmoment0=fits.open(transdict[transitionkeys[transkey]]['filename'])
                transmom0header=transmoment0[0].header
                primaryhduerr.header=transmom0header
                primaryhduerr.header['BTYPE']='Upper-state column density'
                primaryhduerr.header['BUNIT']='cm-2'
                hdulerr=fits.HDUList([primaryhduerr])
                hdulerr.writeto(nuperrorimage_filepath)
            pixelzcoord_nupper+=1
            pixelzcoord_nuperr+=1
            #pdb.set_trace()
        #print(n_us)
        #print(n_uerr)
    print('pixels looped, Nupper calcs complete\n')  

'''
for key in spwdictkeys:
    transitionkeys=list(spwdict[key])
    for transkey in range(len(transitionkeys)):
        nupperimage_filepath=filepath2+'CH3OH~'+transitionkeys[transkey]+'.fits'
        nuperrorimage_filepath=filepath2+'CH3OH~'+transitionkeys[transkey]+'error.fits'
        if os.path.isfile(nupperimage_filepath):
            print(f'{nupperimage_filepath} already exists')
        elif os.path.isfile(nupperimage_filepath)==False:
            nupperimgdata=nugsmap[:,:,transkey]
            primaryhdu=fits.PrimaryHDU(nupperimgdata)
            primaryhdu.header['UNIT']='cm-2'
            hdul.writeto(nupperimage_filepath,overwrite=True)
            hdul=fits.HDUList([primaryhdu])
        elif os.path.isfile(nuperrorimage_filepath):
            print(f'{nuperrorimage_filepath} already exists')
        elif os.path.isfile(nuperrorimage_filepath)==False:
            primaryhduerr=fits.PrimaryHDU(nuperrorimgdata)
            primaryhduerr.header['UNIT']='cm-2'
            hdulerr=fits.HDUList([primaryhduerr])
            hdulerr.writeto(nuperrorimage_filepath)
'''

print('Setting up and executing model fit')
texmap=np.empty((testyshape,testxshape))
ntotmap=np.empty((testyshape,testxshape))

ntoterrmap=np.empty((testyshape,testxshape))
texerrormap=np.empty((testyshape,testxshape))

texsigclipmap=np.empty((testyshape,testxshape))
ntotsigclipmap=np.zeros((testyshape,testxshape))
texsnrmap=np.empty((testyshape,testxshape))
ntotsnrmap=np.zeros((testyshape,testxshape))
numtransmap=np.empty((testyshape,testxshape))
degensforfit=[]
snr=3

fitdict={}

for y in range(testyshape):
    print(f'Start Row {y} fitting')
    for x in range(testxshape):
        rotdiagfilename=rotdiagpath+f'rotdiag_pixel_{y}-{x}_weighted2_sigfigs.png'
        linemod=models.Linear1D(slope=1.0,intercept=14)
        fit=fitting.LinearLSQFitter()
        nupperstofit=[]
        eukstofit=[]
        nuperrors=[]
        qnsfitted=[]
        excludedqns=[]
        excludednuppers=[]
        excludedeuks=[]
        for zed in range(testzshape):
            if nugsmap[y,x,zed] <= 0 or np.isnan(nugsmap[y,x,zed]):
                continue
            elif nugsmap[y,x,zed]/nugserrormap[y,x,zed] == 1:
                #print(f'Excluded line detected: {masterqns[z]}')
                #print('Appending to exclusion lists')
                #excludedqns.append(masterqns[zed])
                #excludednuppers.append(nugsmap[y,x,zed])
                #excludedeuks.append(mastereuks[zed])
                pass
            else:
                nupperstofit.append(nugsmap[y,x,zed])
                eukstofit.append(mastereuks[zed])
                nuperrors.append(nugserrormap[y,x,zed])
                qnsfitted.append(masterqns[zed])
                degensforfit.append(ordereddegens[zed])
        numtransmap[y,x]=len(nupperstofit)        
        if len(nupperstofit)==0:
            obsTex=np.nan
            obsNtot=np.nan
            texmap[y,x]=obsTex
            ntotmap[y,x]=obsNtot
            texsnrmap[y,x]=np.nan
            texsigclipmap[y,x]=obsTex
            texerrormap[y,x]=np.nan
            ntoterrmap[y,x]=np.nan
            ntotsigclipmap[y,x]=np.nan
            ntotsnrmap[y,x]=np.nan
        else:
            #log10nuerr=[]
            errstofit=[]
            log10variances=[]
            log10nuerr=[]
            
            for num in range(len(nupperstofit)):
                templog10=(1/nupperstofit[num])*nuperrors[num]
                temperrfit=1/templog10
                #log10nuerr.append(templog10)
                errstofit.append(temperrfit)
                log10variances.append(templog10**2)
                log10nuerr.append(templog10)
                
            
            #pdb.set_trace()
            fit_lin=fit(linemod,eukstofit,np.log10(nupperstofit),weights=errstofit)
            linemod_euks=np.linspace(min(eukstofit),max(mastereuks),100)
            #print('Model fit complete')
            #print('Compute obsTex and obsNtot')
            obsTrot=-np.log10(np.e)/(fit_lin.slope)
            obsNtot=qrot_partfunc*10**(fit_lin.intercept)#qrot_partfunc*10**(np.log10(nupperstofit[0])+fit_lin.slope*eukstofit[0])
            
            A=np.stack((eukstofit,np.ones_like(eukstofit)),axis=1)
            C=np.diagflat(log10variances)
            atc_1a=np.dot(np.dot(A.T, np.linalg.inv(C)), A)
            if np.linalg.det(atc_1a) == 0:
                print(f'Singular C matrix detected in pixel {y,x}')
                m_unc=np.nan
            else:
                covmat = np.linalg.inv(atc_1a)
                m_unc = covmat[0,0]**0.5
                b_unc = covmat[1,1]**0.5
            
            dobsTrot=np.abs(np.abs(m_unc/fit_lin.slope)*obsTrot*u.K)
            dobsNtot=np.abs(qrot_partfunc*10**(fit_lin.intercept)*(np.log(10)*b_unc))*u.cm**-2#np.sqrt((qrot_partfunc*10**(np.log10(nupperstofit[0])+fit_lin.slope*eukstofit[0])*np.log(10)*eukstofit[0]*m_unc)**2+(qrot_partfunc*10**(np.log10(nupperstofit[0])+fit_lin.slope*eukstofit[0])*(1/(nupperstofit[0]*np.log(10)))*nuperrors[0])**2)*u.cm**-2
            
            sigTrot=(obsTrot*u.K/dobsTrot).to('')
            sigNtot=(obsNtot*u.cm**-2/dobsNtot).to('')
            
            texmap[y,x]=obsTrot
            ntotmap[y,x]=obsNtot
            texerrormap[y,x]=dobsTrot.to('K').value
            ntoterrmap[y,x]=dobsNtot.value
            texsnrmap[y,x]=sigTrot
            ntotsnrmap[y,x]=sigNtot
            
            if sigTrot >= snr:
                texsigclipmap[y,x]=obsTrot
            else:
                texsigclipmap[y,x]=np.nan
                
            if sigNtot >= snr:
                '''
                if obsNtot >= 1e29 or dobsNtot.value <= 1:
                    ntotsigclipmap[y,x]=np.nan
                else:
                '''
                ntotsigclipmap[y,x]=obsNtot
            #elif np.isnan(dobsNtot) or dobsNtot == 0:
            #     ntotsigclipmap[y,x]=np.nan
            else:
                ntotsigclipmap[y,x]=np.nan

            #print('Begin plotting')
            tk='$T_{rot}$'
            ntot='log$_{10}(N_{tot})$'
            cm2='cm$^{-2}$'
            #strdobsntot=str(dobsNtot.value)[0]
            if not np.isfinite(dobsNtot.value) or dobsNtot <= 0 or obsNtot == 0:
                val_dntot='nan'#np.nan
                if not np.isfinite(obsNtot) or obsNtot <= 0:
                    val_ntot=np.nan
                else:
                    val_ntot=round(np.log10(obsNtot))
            else:
                val_dntot=round_to_1(dobsNtot.value/(1*10**int(np.log10(obsNtot))))
                val_ntot=round(np.log10(obsNtot),(len(str(round(sigNtot.value)))-1))

            if not np.isfinite(dobsTrot.value):
                dobsTrot=1e5*u.K

            if y >= 30 and y <=50:
                if x >= 30 and x <= 50:
                    plt.ioff()
                    plt.figure()
                    plt.errorbar(eukstofit,np.log10(nupperstofit),yerr=log10nuerr,fmt='o')
                    plt.plot(linemod_euks,fit_lin(linemod_euks),label=(f'{tk}: {int(round(obsTrot, 0))} $\pm$ {int(round(dobsTrot.value,0))} K\n{ntot}: {val_ntot} $\pm$ {val_dntot} '))
                    plt.xlabel(r'$E_u$ (K)',fontsize=14)
                    plt.ylabel(r'log$_{10}$($N_u$/$g_u$)',fontsize=14)
                    plt.legend()
                    plt.savefig((rotdiagfilename),dpi=100,overwrite=True)
                    plt.close()
                    print('Done.')
                    print(f'Diagram saved to {rotdiagfilename}')
            
plt.ion()
detectnum=0
transmaskarr=np.ma.masked_where(numtransmap<=detectnum,texsigclipmap)
transmasktexmap=transmaskarr.filled(fill_value=np.nan)#np.array(np.ma.masked_where(numtransmap<detectnum,texsigclipmap))
            
transmoment0=fits.open(spws_with_detections[0][master_transkeys[0][0]]['filename'])
transmom0header=transmoment0[0].header

primaryhdutex=fits.PrimaryHDU(texmap)
primaryhdutex.header=transmom0header
primaryhdutex.header['BTYPE']='Excitation temperature'
primaryhdutex.header['BUNIT']='K'
hdultex=fits.HDUList([primaryhdutex])
print(f'Saving raw temperature map at {sourcepath+"texmap_allspw_withnans_weighted.fits"}\n')
hdultex.writeto(sourcepath+'texmap_allspw_withnans_weighted.fits',overwrite=True)

primaryhduntot=fits.PrimaryHDU(ntotmap)
primaryhduntot.header=transmom0header
primaryhduntot.header['BTYPE']='Total column density'
primaryhduntot.header['BUNIT']='cm-2'
hdulntot=fits.HDUList([primaryhduntot])
print(f'Saving raw ntot map at {sourcepath+"ntotmap_allspw_withnans_weighted_useintercept.fits"}\n')
hdulntot.writeto(sourcepath+'ntotmap_allspw_withnans_weighted_useintercept.fits',overwrite=True)

primaryhdutexerr=fits.PrimaryHDU(texerrormap)
primaryhdutexerr.header=transmom0header
primaryhdutexerr.header['BTYPE']='Excitation temperature'
primaryhdutexerr.header['BUNIT']='K'
hdultexerror=fits.HDUList([primaryhdutexerr])
print(f'Saving temperature error map at {sourcepath+"texmap_error_allspw_withnans_weighted.fits"}\n')
hdultexerror.writeto(sourcepath+'texmap_error_allspw_withnans_weighted.fits',overwrite=True)

primaryhdutexclip=fits.PrimaryHDU(texsigclipmap)
primaryhdutexclip.header=transmom0header
primaryhdutexclip.header['BTYPE']='Excitation temperature'
primaryhdutexclip.header['BUNIT']='K'
hdultexclip=fits.HDUList([primaryhdutexclip])
nsigmatexpath=sourcepath+f'texmap_{snr}sigma_allspw_withnans_weighted.fits'
print(f'Saving {snr}sigma temperature map at {nsigmatexpath}')
hdultexclip.writeto(nsigmatexpath,overwrite=True)

primaryhdutexsnr=fits.PrimaryHDU(texsnrmap)
primaryhdutexsnr.header=transmom0header
primaryhdutexsnr.header['BTYPE']='Excitation temperature SNR'
primaryhdutexsnr.header['BUNIT']=''
hdultexsnr=fits.HDUList([primaryhdutexsnr])
print(f'Saving {snr}sigma temperature map at {sourcepath+"texmap_{snr}sigma_allspw_withnans_weighted.fits"}\n')
hdultexsnr.writeto(sourcepath+'texmap_snr_allspw_weighted.fits',overwrite=True)

primaryhdunumtrans=fits.PrimaryHDU(numtransmap)
primaryhdunumtrans.header=transmom0header
primaryhdunumtrans.header['BTYPE']='Number CH3OH 3sigma Detected Transitions'
primaryhdunumtrans.header['BUNIT']=''
hdulnumtrans=fits.HDUList([primaryhdunumtrans])
print(f'Saving number of {snr}sigma detected CH3OH lines map at {sourcepath+"ch3ohdetections_{snr}sigma_allspw_withnans_weighted.fits"}\n')
hdulnumtrans.writeto(sourcepath+f"ch318ohdetections{detectnum}_{snr}sigma_allspw_withnans_weighted.fits",overwrite=True)

primaryhdutransmasktex=fits.PrimaryHDU(transmasktexmap)
primaryhdutransmasktex.header=transmom0header
primaryhdutransmasktex.header['BTYPE']='Excitation temperature'
primaryhdutransmasktex.header['BUNIT']='K'
hdultransmasktex=fits.HDUList([primaryhdutransmasktex])
nsigmatransmaskedpath=sourcepath+f"texmap_{detectnum}transmask_{snr}sigma_allspw_withnans_weighted.fits"
print(f'Saving {detectnum} transition masked temperature map at {nsigmatransmaskedpath}\n')
hdultransmasktex.writeto(nsigmatransmaskedpath,overwrite=True)

primaryhduntoterr=fits.PrimaryHDU(ntoterrmap)
primaryhduntoterr.header=transmom0header
primaryhduntoterr.header['BTYPE']='Total column density error'
primaryhduntoterr.header['BUNIT']='cm-2'
hdulntoterr=fits.HDUList([primaryhduntoterr])
ntoterrpath=sourcepath+f"ntoterrmap_allspw_withnans_weighted_useintercept.fits"
print(f'Saving ntoterr map at {ntoterrpath}\n')
hdulntoterr.writeto(ntoterrpath,overwrite=True)

primaryhduntotsig=fits.PrimaryHDU(ntotsigclipmap)
primaryhduntotsig.header=transmom0header
primaryhduntotsig.header['BTYPE']='Total column density'
primaryhduntotsig.header['BUNIT']='cm-2'
hdulntotsig=fits.HDUList([primaryhduntotsig])
ntotsigpath=sourcepath+f"ntotmap_allspw_withnans_weighted_useintercept_{snr}sigma.fits"
print(f'Saving sigmaclip ntot map at {ntotsigpath}\n')
hdulntotsig.writeto(ntotsigpath,overwrite=True)

primaryhduntotsnr=fits.PrimaryHDU(ntotsnrmap)
primaryhduntotsnr.header=transmom0header
primaryhduntotsnr.header['BTYPE']='Total column density SNR'
primaryhduntotsnr.header['BUNIT']='cm-2/cm-2'
hdulntotsnr=fits.HDUList([primaryhduntotsnr])
ntotsnrpath=sourcepath+f"ntotmap_snr_allspw_withnans_weighted_useintercept.fits"
print(f'Saving ntot snr map at {ntotsnrpath}\n')
hdulntotsnr.writeto(ntotsnrpath,overwrite=True)

nugs_swapaxis2toaxis0=np.swapaxes(nugsmap,0,2)
nugserr_swapaxis2toaxis0=np.swapaxes(nugserrormap,0,2)

nugs_swapxwithy=np.swapaxes(nugs_swapaxis2toaxis0,1,2)
nugserr_swapxwithy=np.swapaxes(nugserr_swapaxis2toaxis0,1,2)

nugscube=fits.PrimaryHDU(nugs_swapxwithy)
nugserrcube=fits.PrimaryHDU(nugserr_swapxwithy)

nugscube.header=transmom0header
nugserrcube.header=transmom0header

nugscube.header['BUNIT']='cm-2'
nugserrcube.header['BUNIT']='cm-2'

nugserrcube.header['BTYPE']='Upper-state column density error'
nugscube.header['BTYPE']='Upper-state column density'

nugshdul=fits.HDUList([nugscube])
nugserrhdul=fits.HDUList([nugserrcube])

print('Saving alltransition and E_U(K) lists\n')
nugshdul.writeto(sourcepath+'alltransitions_nuppers.fits',overwrite=True)
nugserrhdul.writeto(sourcepath+'alltransitions_nupper_error.fits',overwrite=True)

eukqns=QTable(columns=[mastereuks,masterqns,masterlines,ordereddegens], names=['Eupper','QNs','Reference Frequency','Degeneracy'], descriptions=['','',f'z={vlsr}',''])#np.column_stack((mastereuks,masterqns,masterlines,ordereddegens))
eukqns.writeto(sourcepath+'mastereuksqnsfreqsdegens.fits')
#np.savetxt(sourcepath+'mastereuksqnsfreqsdegens.txt',eukqns,fmt='%s',header=f'Methanol transitions, excitation temperatures, and degeneracies used in this folder. Temperatures in units of K, frequencies are redshifted ({z}/{(z*c).to("km s-1")}) and in Hz.\nExcluded lines: {excludedlines[source]}')

'''This wing of the code plots up the temperature and ntot maps'''

print('Begin plotting procedure.\n')

def make_scalebar(ax, left_side, length, color='black', linestyle='-', label='',
                  fontsize=12, text_offset=0.1*u.arcsec, coordsys='icrs'):
    axlims = ax.axis()
    lines = ax.plot(u.Quantity([left_side.ra, left_side.ra-length]),
                    u.Quantity([left_side.dec]*2),
                    color=color, linestyle=linestyle, marker=None,
                    transform=ax.get_transform(coordsys),
                   zorder=3)
    txt = ax.text((left_side.ra-length/2).to(u.deg).value,
                  (left_side.dec+text_offset).to(u.deg).value,
                  label,
                  verticalalignment='bottom',
                  horizontalalignment='center',
                  transform=ax.get_transform(coordsys),
                  color=color,
                  fontsize=fontsize,
                 zorder=2,bbox=dict(facecolor='white', alpha=0.6))
    ax.axis(axlims)
    return lines,txt

colormap= copy.copy(mpl.cm.get_cmap("inferno"))
colormap.set_bad('black')
dGC=8.34*u.kpc#per Meng et al. 2019 https://www.aanda.org/articles/aa/pdf/2019/10/aa35920-19.pdf

plottexhdu=fits.open(nsigmatransmaskedpath)[0]

saveimghome=overleafpath+f'{source}/{sourcelocs[source]}/'
saveimgpath=saveimghome+f'texmap_{snr}sigma_allspw_withnans_weighted.png'

if not os.path.exists(saveimghome):
    os.makedirs(saveimghome)
else:
    print(f"Figure directory {saveimghome} already exists.")
    pass

plottexwcs=WCS(plottexhdu)

sliced=['x','y']
plt.figure()
ax=plt.subplot(projection=plottexwcs,slices=sliced)
plt.rcParams['figure.dpi'] = 300

ra=ax.coords[0]
dec=ax.coords[1]

vmaxdict={}#{'SgrB2S':525,'DSi':320,'DSii':224,'DSiii':300,'DSiv':312,'DSv':280,'DSVI':377,'DSVII':248,'DSVIII':225,'DSIX':215}
if source in vmaxdict.keys():
    sourcevmax=vmaxdict[source]
    plottedtex=ax.imshow(plottexhdu.data,vmax=sourcevmax,cmap=colormap)#,vmax=1000,vmin=10)
else:
    plottedtex=ax.imshow(plottexhdu.data,cmap=colormap)


scaledict={'SgrB2S':5000*u.AU,'DSi':5000*u.AU,'DSii':2000*u.AU,'DSiii':2000*u.AU,'DSiv':2000*u.AU,'DSv':2000*u.AU,'DSVI':5000*u.AU,'DSVII':5000*u.AU,'DSVIII':5000*u.AU,'DSIX':5000*u.AU,'DSX':5000*u.AU}
scale=scaledict[source]
lenn=np.arctan(scale/dGC)

if source == 'DSi':
    print(f'Scalebar source : DSi')
    make_scalebar(ax, coordinates.SkyCoord('17:47:19.5180 -28:23:51.359', unit=(u.hour, u.deg), 
                                       frame='icrs'),
              length=lenn, label=f'{int(scale.value)} AU', fontsize=12, 
              text_offset=0.02*u.arcsec)
    crd=coordinates.SkyCoord('17:47:19.4976 -28:23:51.384', unit=(u.hour,u.deg), frame='icrs')

elif source == 'SgrB2S':
    print(f'Scalebar source: SgrB2S')
    make_scalebar(ax, coordinates.SkyCoord('17:47:20.6618 -28:23:48.734', unit=(u.hour, u.deg), 
                                       frame='icrs'),
              length=lenn, label=f'{int(scale.value)} AU', fontsize=12, 
              text_offset=0.02*u.arcsec)
    crd=coordinates.SkyCoord('17:47:20.6590 -28:23:48.772', unit=(u.hour,u.deg), frame='icrs')
elif source == 'DSii':
    print(f'Scalebar source: DSii')
    make_scalebar(ax, coordinates.SkyCoord('17:47:20.1115 -28:23:47.596', unit=(u.hour, u.deg), 
                                       frame='icrs'),
              length=lenn, label=f'{int(scale.value)} AU', fontsize=12, 
              text_offset=0.02*u.arcsec)
    crd=coordinates.SkyCoord('17:47:20.1087 -28:23:48.634', unit=(u.hour,u.deg), frame='icrs')
    pass
elif source == 'DSiii':
    print(f'Scalebar source: DSiii')
    make_scalebar(ax, coordinates.SkyCoord('17:47:20.0577 -28:23:49.912', unit=(u.hour, u.deg), 
                                       frame='icrs'),
              length=lenn, label=f'{int(scale.value)} AU', fontsize=12, 
              text_offset=0.02*u.arcsec)
    crd=coordinates.SkyCoord('17:47:20.0549 -28:23:49.950', unit=(u.hour,u.deg), frame='icrs')
    
else:
    print(f'Need to make {source} scalebar!!!')
    pass

ra.set_axislabel('RA (J2000)',fontsize=14,minpad=0.9)
ra.set_ticklabel(exclude_overlapping=True)
dec.set_axislabel('Dec (J2000)',fontsize=14,minpad=-0.7)
ax.tick_params(fontsize=14)
cbar1=plt.colorbar(plottedtex,pad=0,label='K')#plt.colorbar()
plt.savefig(saveimgpath,overwrite=True)
plt.show()

print('Cube-->core-->Texmap complete.')
