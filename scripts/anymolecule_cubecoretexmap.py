from spectral_cube import SpectralCube as sc
import matplotlib.pyplot as plt
from astropy.io import fits
import glob
import radio_beam
import os
from astropy.modeling import models, fitting
import pdb
from astropy.wcs import WCS
import matplotlib as mpl
import copy
from astropy import coordinates
from spectral_cube import BooleanArrayMask
from astropy.nddata import Cutout2D
from spectral_cube.io.casa_masks import make_casa_mask
from utilities import *
from manualspeciesabundances import *
from pyspeckit.spectrum.models import lte_molecule
from astropy.table import QTable
from astropy.stats import bootstrap
import astropy.stats
import importlib
from pathlib import Path
import sys
sys.path.append('RotationalTemperatureParameters/')

reference_molecule='CH3OH'
molecule=' C2H5OH '
nospace_molecule=molecule.replace(' ','')
target_molecule_module = importlib.import_module(nospace_molecule)
reference_molecule_module = importlib.import_module(reference_molecule)
# Access attributes from target_molecule_module as needed, e.g.:
# target_molecule_module.some_function()

Splatalogue.QUERY_URL= 'https://splatalogue.online/c_export.php'

'''This wing of the script takes in continuum-subtracted cubes, cuts out a subcube around a region of interest based on a DS9 region, and converts the subcubes into brightness temperature (K) units'''
print('Cube-->Core-->Tex start\n')
print('Begin Jy/beam-to-K and region subcube conversion\n')

source='DSi' #Change this to the source you want to analyze
print(f'Source: {source}\n')
fnum=fields[source]

inpaths=orange_stacontsubpaths
inpath=inpaths[fnum]
beamcubes=glob.glob(inpath+'*.fits')
homes=orange_stdhomedict
home=homes[fnum]
cubes=glob.glob(home+'*pbcor_line.fits')
region=sourceregs[source]
outstatpath_end={1:'OctReimage_K/',10:'field10originals_K/',2:'field2originals_K/',3:'field3originals_K/',7:'field7originals_K/',8:'field8originals_K/'}
outpath=Path("/orange/adamginsburg/sgrb2/2017.1.00114.S/desmond/SgrB2DSminicubes/")/source/outstatpath_end[fnum]
statfixpath=str(Path('/orange/adamginsburg/sgrb2/2017.1.00114.S/desmond/SgrB2DSstatcontfix/')/outstatpath_end[fnum])


regionparams=[float(val) for val in region[9:(len(region)-1)].split(', ')]

'''This wing of the code runs the linelooper LTE modeling and kinetic temperature determination on the newly created region-specific subcubes'''

print('Begin core cube to Tex map process\n') 

'''Collect constants and quantum parameters'''
print('Setting constants')
b_0=target_molecule_module.rotationalconstants[0]
a_0=target_molecule_module.rotationalconstants[1]
c_0=target_molecule_module.rotationalconstants[2]
m=b_0**2/(a_0*c_0)

#Define path names for methanol/parameters used for initial guesses
reorgpath=f'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field{fnum}/CH3OH/{source}/{reference_molecule_module.sourcelocs[source]}'
mastertxttablepath=reorgpath+'mastereuksqnsfreqsdegens.fits'
fwhmpath=glob.glob(reorgpath+'*fwhm*')[0]
trotmappath=reorgpath+'bootstrap_texmap_3sigma_allspw_withnans_weighted.fits'
pixcoords=pixdict[source]

vlsr=target_molecule_module.dopplershifts[source]
z=vlsr/(c.to('km s-1'))
print(f'Doppler shift: {z} / {vlsr}\n')

print('Setting input LTE parameters')
sourcefwhm=fits.getdata(fwhmpath)*u.km/u.s
trotmap=fits.getdata(trotmappath)*u.K
measTrot=trotmap[pixcoords[0],pixcoords[1]]
measlinewidth=sourcefwhm[pixcoords[0],pixcoords[1]]
    
cube_edge_boundary=10*u.MHz
'''Loop through a given list of lines (in Hz), computing and saving moment0 maps of the entered data cube'''
def linelooplte(line_list,line_fwhm):
    print('\nlinelooperLTE...')
    for line in line_list:
        unsliced_transition=line['QNs']
        if unsliced_transition in target_molecule_module.excludedlines[source]:
            print(f'Skipping excluded line: {unsliced_transition}')
            excludedline_index=np.where(safelines['QNs']==unsliced_transition)[0][0]
            safelines.remove_row(excludedline_index)
            #pdb.set_trace()
            continue
        else:
            print(f'\nStart {unsliced_transition} moment0 procedure')
            reffreq=line['ReferenceFrequency']
            line_fwhm_freq=velocitytofreq(line_fwhm,reffreq)
            nu_upper=reffreq+(line_fwhm_freq/2)#Am testing increasing the widths of the slabs, just to make sure we're not unneccessarily omitting signal. Reset back to dividing by 2 if this doesn't work
            nu_lower=reffreq-(line_fwhm_freq/2)#Am testing increasing the widths of the slabs, just to make sure we're not unneccessarily omitting signal. Reset back to dividing by 2 if this doesn't work
            if (max(cube.spectral_axis) - nu_upper) <= cube_edge_boundary or (nu_lower - min(cube.spectral_axis)) <= cube_edge_boundary:
                print(f'{unsliced_transition} is too close to the cube edges. Skipping\n')
                continue
            print(f'Make spectral slab between {nu_lower} and {nu_upper}')
            slab=cube.spectral_slab(nu_upper,nu_lower)
            #pdb.set_trace()
            print('Slicing quantum numbers')
            transition=stringmanipulationdict[molecule](line['QNs'])#qn_replace(quantum_numbers[i])
            moment0filename = home/f"{nospace_molecule}~{transition}_raw.fits"
            maskedmom0fn = home/f"{nospace_molecule}~{transition}_masked.fits"
            maskedmom0errfn = home/f"{nospace_molecule}~{transition}_error.fits"
            slabfilename = slabpath/f"{nospace_molecule}~{transition}_slab.fits"
            maskedslabfn = slabpath/f"{nospace_molecule}~{transition}_maskedslab.fits"
            peakintfn = home/f"{nospace_molecule}~{transition}_peakint.fits"
            moment1filename = sourcepath/'mom1'/f"{nospace_molecule}~{transition}.fits"
            moment2filename = sourcepath/'mom2'/f"{nospace_molecule}~{transition}_var.fits"
            fwhmfilename = sourcepath/'mom2'/f"{nospace_molecule}~{transition}_fwhm.fits"
            masterslicedqns.append(transition)
            #Check if the moment0 file already exists, if it does, skip the rest of the loop
            if maskedmom0fn.exists():
                print(f'{moment0filename} already exists.')
                isfile_hdu=fits.open(maskedmom0fn)
                isfileheader=isfile_hdu[0].header
                isfilemom0=fits.getdata(maskedmom0fn)*u.K*u.km/u.s
                isfilemom0err=fits.getdata(maskedmom0errfn)*u.K*u.km/u.s
                mastermom0s.update({unsliced_transition:isfilemom0})
                #masterslicedqns.append(transition)
                mastermom0errors.update({unsliced_transition:isfilemom0err})
                headers_for_output_maps.update({unsliced_transition:isfileheader})
                print('\nDictionaries populated for this transition.')
            else: 
                print('Commence moment0 procedure\n')
                print(f'Create {line['QNs']} spatial-velocity mask')
    
                slab=slab.with_spectral_unit((u.km/u.s),velocity_convention='radio',rest_value=reffreq)
                slab3sigmamask=slab > (3*stdcutout.data)
                slab=slab.with_mask(slab3sigmamask)
                slabspecax=slab.spectral_axis
                slabmom1=slab.moment1()
                slabmom2=slab.moment2()
                slabfwhm=slab.linewidth_fwhm()
                cubemask=(slabspecax[:,None,None] < (velocityfield_representative + fwhm_representative)[None,:,:]) & (slabspecax[:,None,None] > (velocityfield_representative - fwhm_representative)[None,:,:])
                #oldstyleslab=oldstyleslab.with_spectral_unit((u.km/u.s),velocity_convention='radio',rest_value=lines[i])
    
                print('Masking spectral slab')
                maskedslab=slab.with_mask(cubemask)
                print('Computing unmasked moment0...\n')
                slabmom0=slab.moment0()
                print('Computing masked moment0...\n')
                contmom0=reprojcont_K*slabfwhm#continuum sanity check
                if len(doublet) > 0:
                    pdb.set_trace()
                    for path_transitiontable in doublet: #Find table containing the transition's self-contaminants
                        transitiontable=QTable.read(path_transitiontable)
                        if transition in transitiontable['QNs']:
                            targetdoublettransition=np.where(transition in transitiontable['QNs'])[0] #Find the target transition
                            all_companions_combined_flux=np.sum(doublettable['ModelBrightness']) #Compute the combined flux of the target and its companions
                            targettransition_modelflux=transitiontable['ModelBrightness'][targetdoublettransition] #Select the target transition's flux
                            target_flux_to_total_flux_ratio=targettransition_modelflux/all_companions_combined_flux #Compute the ratio of the target transition's flux to the total flux of the target and its companions
                            maskslabmom0=(maskedslab.moment0()+contmom0)/target_flux_to_total_flux_ratio #Scale the measured flux by the ratio of the target transition's flux to the total flux of the target and its companions
                            print(f'\nDoublet line identified: {unsliced_transition}')
                            print(f'Value scaled by factor of {target_flux_to_total_flux_ratio} to compensate for line blending.\n')
                            sys.exit()
                        else:
                            maskslabmom0=maskedslab.moment0()+contmom0
                else:
                    maskslabmom0=maskedslab.moment0()+contmom0
                print('Computing peak intensity')
                maskedpeakint=maskedslab.max(axis=0)
    
                print('\nSaving...')
                slabmom0.write(moment0filename,overwrite=True)
                maskslabmom0.write(maskedmom0fn)
                kkmsstdarray=stdcutout.data*fwhm_representative
    
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
                kkmshdul.writeto(maskedmom0errfn)
                print(f'Saving moment1 to {moment1filename}\n')
                slabmom1.write(moment1filename)
                print(f'Saving moment 2s to {moment2filename} and {fwhmfilename}\n')
                slabmom2.write(moment2filename)
                slabfwhm.write(fwhmfilename)
                print(f'Saving peak intensity to {peakintfn}')
                maskedpeakint.write(peakintfn)
                
                mastermom0s.update({unsliced_transition:maskslabmom0})
                mastermom0errors.update({unsliced_transition:kkmsstdarray})
                headers_for_output_maps.update({unsliced_transition:maskslabmom0.header})
                print(f'{line['QNs']} calculations complete.\n')
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

incubes=glob.glob(str(outpath/"*pbcor_line.fits"))
incubes.sort()
images=['spw0','spw1','spw2','spw3']
datacubes=incubes
assert 'spw0' in datacubes[0], 'Cube list out of order'

stdhome=orange_stdhomedict[fnum]

targetworldcrd=targetworldcrds[source]

contpath=reorgpath+'reprojectedcontinuum.fits'
reprojcontfits=fits.open(contpath)
reprojcont=reprojcontfits[0].data*u.Jy
reprojcontrestfreq=225*u.GHz#manual addition 11/9/2022, wiggle room w/i GHz
cntmbeam=radio_beam.Beam.from_fits_header(reprojcontfits[0].header)
reprojcont_K=reprojcont.to('K',cntmbeam.jtok_equiv(reprojcontrestfreq))

representativelws=measlinewidth

sourcepath=Path(f'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field{fnum}/{nospace_molecule}/{source}/{target_molecule_module.sourcelocs[source]}')
nupperpath=sourcepath/'nuppers/'
stdpath=sourcepath/'errorimgs/std/'
slabpath=sourcepath/'spectralslabs/km_s/'
mom0path=sourcepath/'mom0/'
rotdiagpath=sourcepath/'pixelwiserotationaldiagrams/'
figpath=sourcepath/'figures/'

overleafpath=Path("/blue/adamginsburg/d.jeff/repos/C2H5OH_SgrB2DS/figures/")

#Create the file tree for all the output files
if slabpath.is_dir():
    print(f'Source path directory tree {sourcepath} already exists.\n')
    mom1path=sourcepath/'mom1/'
    mom2path=sourcepath/'mom2/'
    if mom1path.is_dir() and mom2path.is_dir():
        print('Moment 1 and 2 directories already exist.')
    else:
        print(f'Creating moment 1 directory at {mom1path}')
        mom1path.mkdir()
        print(f'Creating moment 2 directory at {mom2path}')
        mom2path.mkdir()
else:
    print(f'Making source path {sourcepath}')
    os.makedirs(sourcepath)
    print(f'Making nupper folder {nupperpath}')
    nupperpath.mkdir()
    print(f'Making error folder {stdpath}')
    os.makedirs(stdpath)
    print(f'Making spectral slab folder {slabpath}\n')
    os.makedirs(slabpath)
    
    mom0path=sourcepath/'mom0/'
    mom1path=sourcepath/'mom1/'
    mom2path=sourcepath/'mom2/'
    print(f'Making moment 0 folder {mom0path}')
    mom0path.mkdir()
    print(f'Creating moment 1 directory at {mom1path}')
    mom1path.mkdir()
    print(f'Creating moment 2 directory at {mom2path}')
    mom2path.mkdir()

    print(f'Making rotational diagram folder')
    rotdiagpath.mkdir()
    print(f'Making figures folder')
    figpath.mkdir()

masterslicedqns=[]

mastermom0s={}
mastermom0errors={}
headers_for_output_maps={}

#Check which linelist the molecule is in and set the master catdir and catdirtag
if linelistdict[molecule] == 'CDMS':
    print(f'Using CDMS linelist for {molecule}')
    master_catdir=CDMS.get_species_table()
    if target_molecule_module.catdirtag:
        catdirtag=target_molecule_module.catdirtag
        print(f'Catdir tag: {catdirtag}')
    else:
        print(f'No catdirtag found for {molecule} ')
        sys.exit()
elif linelistdict[molecule] == 'JPL':
    print(f'Using JPL linelist for {molecule}')
    master_catdir=JPLSpec.get_species_table()
    if target_molecule_module.catdirtag:
        catdirtag=target_molecule_module.catdirtag
        print(f'Catdir tag: {catdirtag}')
    else:
        print(f'No catdirtag found for {molecule}')
        sys.exit()

#Select the target molecule's catdir table and extract the partition function at 300 K
target_molecule_catdir=master_catdir[master_catdir['tag'] == catdirtag]
catdir_qrot300=10**target_molecule_catdir['lg(Q(300))']

#Set up paths for representative line objects
representative_filename_base=sourcepath/f'{target_molecule_module.representativelines[source]}repline_'
rep_mom1=f'{representative_filename_base}mom1.fits'
rep_fwhm=f'{representative_filename_base}fwhm.fits'
rep_slab=f'{representative_filename_base}slab.fits'
rep_maskedslab=f'{representative_filename_base}maskedslab.fits'

#Select the region that defines the dimensions of the output maps
regiondims=regionparams[3]

#Set up paths for safelines and doublets
safelinemol=molecule.replace(' ','')
linemodelpath=linemodelhome/linemodelversion/source
safelinepathbase=linemodelpath/'safelines/'
safelinepath=safelinepathbase/f'{safelinemol}.fits'
doubletpath=safelinepathbase/f'{nospace_molecule}/doublets/*.fits'
safelines=QTable.read(safelinepath)
doublet=glob.glob(str(doubletpath))

# Check if representative line objects already exist and populate variables if they do
if os.path.isfile(rep_mom1):
    print(f'{source} representative line objects already exist.')
    print(f'Grabbing spectralslab from {rep_slab}')
    spectralslab_rep_3sigma=fits.getdata(rep_maskedslab)*u.K
    print(f'Grabbing mom1 from {rep_mom1}')
    velocityfield_representative=fits.getdata(rep_mom1)*u.km/u.s
    print(f'Grabbing fwhm from {rep_fwhm}\n')
    fwhm_representative=fits.getdata(rep_fwhm)*u.km/u.s
else:
    print(f'{target_molecule_module.representativelines[source]} representative line objects will be computed for {source}.\n')
    pass

#Loop over datacubes and create maps of observational parameters
for imgnum in range(len(datacubes)):
    print(f'Accessing data cube {datacubes[imgnum]}')
    assert images[imgnum] in datacubes[imgnum], f'{images[imgnum]} not in filename {datacubes[imgnum]}'
    home=sourcepath/'mom0/'
    cube=sc.read(datacubes[imgnum])
    header=fits.getheader(datacubes[imgnum])
    
    #The stdimage comes from the full-sized image cube, we need to make a cutout around the hot core
    stdimage=fits.open(stdhome+images[imgnum]+'minimize.image.pbcor_noise.fits')
    stdcellsize=(np.abs(stdimage[0].header['CDELT1']*u.deg)).to('arcsec')
    stdcutoutsize=round(((float(regiondims)*u.deg)/stdcellsize).to('').value)
    stddata=stdimage[0].data*u.K
    
    print('Acquiring cube rest frequency and computing target pixel coordinates')
    spwrestfreq=header['RESTFRQ']*u.Hz
    freqs=cube.spectral_axis#Hz
    freqflip=False
    if freqs[1] < freqs[0]:
        freqs=freqs[::-1]
        freqflip=True
        print('Corrected decreasing frequency axis')
    
    assert max(freqs) > min(freqs), 'Inverted spectral axis'
    print('Passed increasing spectral axis check')
    
    cube_w=cube.wcs
    stdwcs=WCS(stdimage[0].header)
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
        print(f'Creating {target_molecule_module.representativelines[source]} representative line data objects for {source}')
        print(f'Opening file {datacubes[target_molecule_module.representativecubes[source]]}')
        reffreq_repline=target_molecule_module.restfreq_representativeline[source]/(1+z)
        repcube=sc.read(datacubes[target_molecule_module.representativecubes[source]])
        assert reffreq_repline > min(repcube.spectral_axis) and reffreq_repline < max(repcube.spectral_axis), f'Incorrect representative cube chosen.\nCube limits: {min(repcube.spectral_axis)},{max(repcube.spectral_axis)}\nRepresentative line reference frequency: {reffreq_repline}'
        
        repstdmain=fits.open(stdhome+f'spw{target_molecule_module.representativecubes[source]}minimize.image.pbcor_noise.fits')
        repstdmain_data=repstdmain[0].data*u.K
        repstdmain_wcs=WCS(repstdmain[0].header)
        repstdmain_cellsize=(np.abs(repstdmain[0].header['CDELT1']*u.deg)).to('arcsec')
        repstdmain_pixcrd=repstdmain_wcs.wcs_world2pix(targetworldcrd,1,ra_dec_order=True)
        repstdmain_xcrd,repstdmain_ycrd=int(round(repstdmain_pixcrd[1][0])),int(round(repstdmain_pixcrd[1][1]))
        
        assert repstdmain_xcrd and repstdmain_ycrd > 0, 'Negative representative std pixel coords'
        
        repstdcutoutsize=round(((float(regiondims)*u.deg)/stdcellsize).to('').value)
        repstdcutout=Cutout2D(repstdmain_data,(repstdmain_xcrd,repstdmain_ycrd), repstdcutoutsize)
        upperfreq=reffreq_repline+velocitytofreq(representativelws,reffreq_repline)#Pretty sure I neglected to divide the width by 2 here, may need to do so later on
        lowerfreq=reffreq_repline-velocitytofreq(representativelws,reffreq_repline)#Pretty sure I neglected to divide the width by 2 here, may need to do so later on
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

    continuuminpix=reprojcont_K[pixcoords[0],pixcoords[1]]

    linesbelowmax=safelines['ReferenceFrequency']<=max(cube.spectral_axis)
    linesabovemin=safelines['ReferenceFrequency']>=min(cube.spectral_axis)
    lines_in_spw=safelines[linesbelowmax*linesabovemin]
    
    linewidth=fwhm_representative[pixcoords[0],pixcoords[1]]#representativelws
    linewidth_freq=velocitytofreq(linewidth,target_molecule_module.restfreq_representativeline[source])

    linelooplte(lines_in_spw,linewidth)
    print(f'Finished loop for {images[imgnum]}\n')

'''
if os.path.isfile(picklepath):
    print(f'pickle {picklepath} already exists.')
else:
    print('Saving dictionary pickle...')
    f=open(picklepath,'wb')
    pickle.dump(spwdict,f)
    f.close()
    print(f'Dictionary pickle saved at {picklepath}')
'''

print('Begin fitting routines\nCompute N_uppers')
nupper_gmaps={}
error_nupper_gmaps={}
master_nupperfilepaths=[]
master_errornupperfilepaths=[]
for count, transition in enumerate(safelines):
    current_qns=transition['QNs']
    current_deg=transition['Degeneracy']
    print(f'Transition: {current_qns}')
    nupperimage_filepath=nupperpath/f'{nospace_molecule}~{masterslicedqns[count]}data.fits'
    nuperrorimage_filepath=nupperpath/f'{nospace_molecule}~{masterslicedqns[count]}error.fits'
    
    nupperimgexists=False
    nuperrorimgexists=False
    if nupperimage_filepath.exists():
        print(f'{nupperimage_filepath} already exists. Populating Nupper/g dictionary.\n')
        nupper_g=fits.getdata(nupperimage_filepath)*u.cm**-2
        error_nupperg=fits.getdata(nuperrorimage_filepath)*u.cm**-2
        nupper_gmaps.update({current_qns:nupper_g})
        error_nupper_gmaps.update({current_qns:error_nupperg})
        nupperimgexists=True
        nuperrorimgexists=True
        master_nupperfilepaths.append(str(nupperimage_filepath))
        master_errornupperfilepaths.append(str(nuperrorimage_filepath))
    elif not nupperimgexists:
        nupper,error_nupper=N_u(transition['RestFrequency'],transition['Aij'],mastermom0s[current_qns],
                                        mastermom0errors[current_qns]) 
        nupperg=nupper/current_deg
        error_nupperg=error_nupper/current_deg
        nupper_gmaps.update({current_qns:nupperg})
        error_nupper_gmaps.update({current_qns:error_nupperg})
        if not nupperimgexists:
            nupperimgdata=nupperg.value#nugsmap[:,:,pixelzcoord_nupper]
            primaryhdu=fits.PrimaryHDU(nupperimgdata)
            templatemom0header=headers_for_output_maps[current_qns]
            primaryhdu.header=templatemom0header
            primaryhdu.header['BTYPE']='Upper-state column density'
            primaryhdu.header['BUNIT']='cm-2'
            hdul=fits.HDUList([primaryhdu])
            hdul.writeto(nupperimage_filepath,overwrite=True)
            master_nupperfilepaths.append(str(nupperimage_filepath))
        if not nuperrorimgexists:
            nuperrorimgdata=error_nupperg.value#nugserrormap[:,:,pixelzcoord_nuperr]
            primaryhduerr=fits.PrimaryHDU(nuperrorimgdata)
            templatemom0header=headers_for_output_maps[current_qns]
            primaryhduerr.header=templatemom0header
            primaryhduerr.header['BTYPE']='Upper-state column density error'
            primaryhduerr.header['BUNIT']='cm-2'
            hdulerr=fits.HDUList([primaryhduerr])
            hdulerr.writeto(nuperrorimage_filepath)
            master_errornupperfilepaths.append(str(nuperrorimage_filepath))

print('Nupper calcs complete\n')  
print('Setting up and executing model fit')
yshape=stdcutout.shape[0]
xshape=stdcutout.shape[1]
texmap=np.empty((yshape,xshape))
texerrormap=np.empty((yshape,xshape))

ntotmap=np.empty((yshape,xshape))
ntoterrmap=np.empty((yshape,xshape))

texsigclipmap=np.empty((yshape,xshape))
ntotsigclipmap=np.zeros((yshape,xshape))
texsnrmap=np.empty((yshape,xshape))
ntotsnrmap=np.zeros((yshape,xshape))
numtransmap=np.empty((yshape,xshape))
snr=3

print(f'Starting rotational {molecule} diagram loops')
centers={'SgrB2S':[80,50],'DSi':[35,50],'DSii':[15,30],'DSiii':[20,32],'DSiv':[20,40],'DSv':[15,25],'DSVI':[55,66],'DSVII':[70,80],'DSVIII':[40,60]}
bootstraps=True
midloop_rotdiagrams=True

for y in range(yshape):
    print(f'Starting row {y} of {yshape}') 
    for x in range(xshape):
        linemod=models.Linear1D(slope=1.0,intercept=14)
        fit=fitting.LinearLSQFitter()
        nuppergstofit=[]
        eukstofit=[]
        nuppergerrorstofit=[]
        qnstofit=[]
        for count, transition in enumerate(nupper_gmaps.keys()):
            if nupper_gmaps[transition][y,x] >= 0 and np.isfinite(nupper_gmaps[transition][y,x]):#or (nupper_gmaps[transition][y,x]/error_nupper_gmaps[transition][y,x]) < snr:
                nuppergstofit.append(nupper_gmaps[transition][y,x].value)
                eukstofit.append(safelines['Eupper'][count].value)
                nuppergerrorstofit.append(error_nupper_gmaps[transition][y,x].value)
                qnstofit.append(safelines['QNs'][count])
        numtransmap[y,x]=len(nuppergstofit)
        if len(nuppergstofit)==0:
            texmap[y,x]=np.nan
            ntotmap[y,x]=np.nan
            texsnrmap[y,x]=np.nan
            texsigclipmap[y,x]=np.nan
            texerrormap[y,x]=np.nan
            ntoterrmap[y,x]=np.nan
            ntotsigclipmap[y,x]=np.nan
            ntotsnrmap[y,x]=np.nan
        else:
            log10nuerr=[]
            weightstofit=[]
            log10variances=[]
            for nuppergtofit,nuppergerrortofit in zip(nuppergstofit,nuppergerrorstofit):
                convert_linear_error_to_log10_error=nuppergerrortofit/nuppergtofit
                tempweight=1/convert_linear_error_to_log10_error
                log10nuerr.append(convert_linear_error_to_log10_error)
                weightstofit.append(tempweight)
                log10variances.append(convert_linear_error_to_log10_error**2)
            
            fit_lin=fit(linemod,eukstofit,np.log10(nuppergstofit),weights=weightstofit)
            linemod_euks=np.linspace(min(eukstofit),max(safelines['Eupper'].value),100)#Goes out to the max Eupper because these values are used for plotting
            obsTrot=-np.log10(np.e)/(fit_lin.slope)
            obsNtot=qrot_partfunc*10**(fit_lin.intercept)
            if bootstraps:
                bslist=[]
                for nug,euk,weight,var in zip(np.log10(nuppergstofit),eukstofit,weightstofit,log10variances):#10
                    bslist.append((nug,euk,weight,var))
                
                numboots=1000
                bootresult=bootstrap(np.array(bslist),numboots)
                
                bootlines=[]
                bootTrots=[]
                bootNtots=[]#The fitted ntots
                bootIntercepts=[]#The fitted interecepts
                bootSlopes=[]
                for boot in bootresult:
                    tempfit=fit(linemod,boot[:,1],boot[:,0],weights=boot[:,2])
                    bootlines.append(tempfit)
                for line in bootlines:
                    tempbootTrot=-np.log10(np.e)/(line.slope)
                    tempbootslope=line.slope.value
                    tempbootNtot=qrot_partfunc*10**(line.intercept)
                    tempbootintNtot=line.intercept.value
                    if line.slope >= 0:#tempbootTrot >= 1000 or tempbootTrot <= 0:
                        continue
                    else:
                        bootTrots.append(tempbootTrot)
                        bootNtots.append(tempbootNtot)
                        bootIntercepts.append(tempbootintNtot)
                        bootSlopes.append(tempbootslope)
                
                bootstrap_error_Trot=astropy.stats.mad_std(bootTrots)
                slope_bootstd=astropy.stats.mad_std(bootSlopes)
                madstd_interceptfits=astropy.stats.mad_std(bootIntercepts)
                bootstrap_error_Ntot=astropy.stats.mad_std(bootNtots)

                dobsTrot=bootstrap_error_Trot
                dobsNtot=bootstrap_error_Ntot
                sigTrot=(obsTrot/dobsTrot)
                sigNtot=(obsNtot/dobsNtot)
            '''
            else:
                A=np.stack((eukstofit,np.ones_like(eukstofit)),axis=1)
                C=np.diagflat(log10variances)
                atc_1a=np.dot(np.dot(A.T, np.linalg.inv(C)), A)
                if np.linalg.det(atc_1a) == 0:
                    #print(f'Singular C matrix detected in pixel {y,x}')
                    m_unc=np.nan
                else:
                    covmat = np.linalg.inv(atc_1a)
                    m_unc = covmat[0,0]**0.5
                    b_unc = covmat[1,1]**0.5
                
                dobsTrot=np.abs(np.abs(m_unc/fit_lin.slope)*obsTrot*u.K)
                dobsNtot=np.abs(qrot_partfunc*10**(fit_lin.intercept)*(np.log(10)*b_unc))*u.cm**-2
                sigTrot=(obsTrot*u.K/dobsTrot).to('')
                sigNtot=(obsNtot*u.cm**-2/dobsNtot).to('')
            '''
            texmap[y,x]=obsTrot
            ntotmap[y,x]=obsNtot
            texerrormap[y,x]=dobsTrot
            ntoterrmap[y,x]=dobsNtot
            texsnrmap[y,x]=sigTrot
            ntotsnrmap[y,x]=sigNtot
                
            if sigTrot >= snr:
                texsigclipmap[y,x]=obsTrot
            else:
                texsigclipmap[y,x]=np.nan
                
            if sigNtot >= snr:
                ntotsigclipmap[y,x]=obsNtot
            else:
                ntotsigclipmap[y,x]=np.nan

            if midloop_rotdiagrams:
                tk='$T_{rot}$'
                ntot='log$_{10}(N_{tot})$'
                cm2='cm$^{-2}$'

                if np.isfinite(obsTrot) and np.isfinite(obsNtot) and np.isfinite(dobsTrot) and np.isfinite(dobsNtot):
                    if y >= centers[source][0] and y <=centers[source][1]:
                        if x >= centers[source][0] and x <= centers[source][1]:
                            plot_error_Ntot=round((dobsNtot/obsNtot),1)
                            val_ntot=round(np.log10(obsNtot),1)
                            rotdiagfilename=rotdiagpath/f'rotdiag_pixel_{y}-{x}.png'
                            plt.ioff()
                            plt.figure()
                            plt.errorbar(eukstofit,np.log10(nuppergstofit),yerr=convert_linear_error_to_log10_error,
                                     fmt='o')
                            plt.plot(linemod_euks,fit_lin(linemod_euks),label=(f'{tk}: {int(obsTrot)} $\pm$ {int(dobsTrot)} K\n{ntot}: {val_ntot} $\pm$ {plot_error_Ntot}'))
                            
                            for linmod in bootlines:
                                plt.plot(linemod_euks,linmod(linemod_euks),color='black',alpha=0.03,zorder=0)
                                
                            plt.xlabel(r'$E_u$ (K)',fontsize=14)
                            plt.ylabel(r'log$_{10}$($N_u$/$g_u$)',fontsize=14)
                            plt.legend()
                            plt.savefig((rotdiagfilename),dpi=100)
                            plt.close()
                            print(f'Diagram saved to {rotdiagfilename}')

plt.ion()
#detectnum=0
#transmaskarr=np.ma.masked_where(numtransmap<=detectnum,texsigclipmap)
#transmasktexmap=transmaskarr.filled(fill_value=np.nan)#np.array(np.ma.masked_where(numtransmap<detectnum,texsigclipmap))
            
templatemom0header=headers_for_output_maps['14(1,14)(0)-13(1,13)(0)']

primaryhdutex=fits.PrimaryHDU(texmap)
primaryhdutex.header=templatemom0header
primaryhdutex.header['BTYPE']='Excitation temperature'
primaryhdutex.header['BUNIT']='K'
hdultex=fits.HDUList([primaryhdutex])
print(f'Saving raw temperature map at {sourcepath/"texmap_allspw_withnans_weighted.fits"}\n')
hdultex.writeto(sourcepath/'texmap_allspw_withnans_weighted.fits',overwrite=True)

primaryhduntot=fits.PrimaryHDU(ntotmap)
primaryhduntot.header=templatemom0header
primaryhduntot.header['BTYPE']='Total column density'
primaryhduntot.header['BUNIT']='cm-2'
hdulntot=fits.HDUList([primaryhduntot])
print(f'Saving raw ntot map at {sourcepath/"ntotmap_allspw_withnans_weighted_useintercept.fits"}\n')
hdulntot.writeto(sourcepath/'ntotmap_allspw_withnans_weighted_useintercept.fits',overwrite=True)

primaryhdutexerr=fits.PrimaryHDU(texerrormap)
primaryhdutexerr.header=templatemom0header
primaryhdutexerr.header['BTYPE']='Excitation temperature'
primaryhdutexerr.header['BUNIT']='K'
hdultexerror=fits.HDUList([primaryhdutexerr])
print(f'Saving temperature error map at {sourcepath/"texmap_error_allspw_withnans_weighted.fits"}\n')
hdultexerror.writeto(sourcepath/'texmap_error_allspw_withnans_weighted.fits',overwrite=True)

primaryhdutexclip=fits.PrimaryHDU(texsigclipmap)
primaryhdutexclip.header=templatemom0header
primaryhdutexclip.header['BTYPE']='Excitation temperature'
primaryhdutexclip.header['BUNIT']='K'
hdultexclip=fits.HDUList([primaryhdutexclip])
nsigmatexpath=sourcepath/f'texmap_{snr}sigma_allspw_withnans_weighted.fits'
print(f'Saving {snr}sigma temperature map at {nsigmatexpath}')
hdultexclip.writeto(nsigmatexpath,overwrite=True)

primaryhdutexsnr=fits.PrimaryHDU(texsnrmap)
primaryhdutexsnr.header=templatemom0header
primaryhdutexsnr.header['BTYPE']='Excitation temperature SNR'
primaryhdutexsnr.header['BUNIT']=''
hdultexsnr=fits.HDUList([primaryhdutexsnr])
print(f'Saving {snr}sigma temperature map at {sourcepath}texmap_snr_allspw_weighted.fits\n')
hdultexsnr.writeto(sourcepath/'texmap_snr_allspw_weighted.fits',overwrite=True)

primaryhdunumtrans=fits.PrimaryHDU(numtransmap)
primaryhdunumtrans.header=templatemom0header
primaryhdunumtrans.header['BTYPE']='Number CH3OH 3sigma Detected Transitions'
primaryhdunumtrans.header['BUNIT']=''
hdulnumtrans=fits.HDUList([primaryhdunumtrans])
print(f'Saving number of {snr}sigma detected CH3OH lines map at {sourcepath/"ch3ohdetections_{snr}sigma_allspw_withnans_weighted.fits"}\n')
hdulnumtrans.writeto(sourcepath/f"{nospace_molecule}detections_{snr}sigma_allspw_withnans_weighted.fits",overwrite=True)

'''
primaryhdutransmasktex=fits.PrimaryHDU(transmasktexmap)
primaryhdutransmasktex.header=templatemom0header
primaryhdutransmasktex.header['BTYPE']='Excitation temperature'
primaryhdutransmasktex.header['BUNIT']='K'
hdultransmasktex=fits.HDUList([primaryhdutransmasktex])
nsigmatransmaskedpath=sourcepath+f"texmap_{detectnum}transmask_{snr}sigma_allspw_withnans_weighted.fits"
print(f'Saving {detectnum} transition masked temperature map at {nsigmatransmaskedpath}\n')
hdultransmasktex.writeto(nsigmatransmaskedpath,overwrite=True)
'''
primaryhduntoterr=fits.PrimaryHDU(ntoterrmap)
primaryhduntoterr.header=templatemom0header
primaryhduntoterr.header['BTYPE']='Total column density error'
primaryhduntoterr.header['BUNIT']='cm-2'
hdulntoterr=fits.HDUList([primaryhduntoterr])
ntoterrpath=sourcepath/f"ntoterrmap_allspw_withnans_weighted_useintercept.fits"
print(f'Saving ntoterr map at {ntoterrpath}\n')
hdulntoterr.writeto(ntoterrpath,overwrite=True)

primaryhduntotsig=fits.PrimaryHDU(ntotsigclipmap)
primaryhduntotsig.header=templatemom0header
primaryhduntotsig.header['BTYPE']='Total column density'
primaryhduntotsig.header['BUNIT']='cm-2'
hdulntotsig=fits.HDUList([primaryhduntotsig])
ntotsigpath=sourcepath/f"ntotmap_allspw_withnans_weighted_useintercept_{snr}sigma.fits"
print(f'Saving sigmaclip ntot map at {ntotsigpath}\n')
hdulntotsig.writeto(ntotsigpath,overwrite=True)

primaryhduntotsnr=fits.PrimaryHDU(ntotsnrmap)
primaryhduntotsnr.header=templatemom0header
primaryhduntotsnr.header['BTYPE']='Total column density SNR'
primaryhduntotsnr.header['BUNIT']='cm-2/cm-2'
hdulntotsnr=fits.HDUList([primaryhduntotsnr])
ntotsnrpath=sourcepath/f"ntotmap_snr_allspw_withnans_weighted_useintercept.fits"
print(f'Saving ntot snr map at {ntotsnrpath}\n')
hdulntotsnr.writeto(ntotsnrpath,overwrite=True)

print('Saving table of used lines\n')
safelines.add_columns([master_nupperfilepaths,master_errornupperfilepaths],names=['MeasuredNupperPath','MeasuredNupperErrorPath'])
eukqns=safelines#QTable([mastereuks,masterqns,masterlines,masterdegens], names=['Eupper','QNs','Reference Frequency','Degeneracy'], descriptions=['','',f'z={vlsr}',''])#np.column_stack((mastereuks,masterqns,masterlines,ordereddegens))
eukqns.write(sourcepath/'table_of_used_lines_and_parameters.fits',overwrite=True)
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

plottexhdu=fits.open(nsigmatexpath)[0]

saveimghome=overleafpath/f'{source}/{target_molecule_module.sourcelocs[source]}'
saveimgpath=saveimghome/f'texmap_{snr}sigma_allspw_withnans_weighted.png'

if not saveimghome.exists():
    print(f'Creating figure directory {saveimghome}')
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
plt.savefig(saveimgpath)
plt.show()

print('Cube-->core-->Texmap complete.')
