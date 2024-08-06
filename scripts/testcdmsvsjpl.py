import numpy as np
from pyspeckit.spectrum.models import lte_molecule
from pyspeckit.spectrum.models.lte_molecule import get_molecular_parameters
from astroquery.linelists.cdms import CDMS
import astropy.units as u
import matplotlib.pyplot as plt
from astroquery.splatalogue import utils, Splatalogue
import scipy.constants as cnst
from astropy.io import fits
import glob
from astropy.modeling import models#Fittable1DModel, Parameter, fitting
from utilities import *#Q_rot_asym,mulu,vradio,t_rad,nupper_estimated,opticaldepth,qngrabber
import matplotlib as mpl
import pdb
import sys
from manualspeciesabundances import *
from spectral_cube import SpectralCube as sc
import radio_beam

Splatalogue.QUERY_URL= 'https://splatalogue.online/c_export.php'

mpl.interactive(True)

plt.close('all')

def lineprofile(sigma,nu_0,nu):
    return (1/(np.sqrt(2*np.pi)*sigma))*np.exp(-(nu-nu_0)**2/(2*sigma**2))

def cdms_get_molecule_name(my_molecule_name, **kwargs):
    basename = dict(CDMS.query_lines(min_frequency=1*u.GHz, max_frequency=500*u.GHz, molecule=my_molecule_name, parse_name_locally=True, get_query_payload=True, **kwargs))['Molecules']
    return " ".join(basename.split(" ")[1:])

'''Collect constants for N_tot and N_upper calculations'''

source='DSi'
fnum=fields[source]
dpi={0:150,1:300}
mode=dpi[0]

c=cnst.c*u.m/u.s
k=cnst.k*u.J/u.K
h=cnst.h*u.J*u.s
sigma_sb=cnst.sigma*u.W/((u.m)**(2)*(u.K)**(4))
b_0=24679.98*u.MHz
a_0=127484*u.MHz
c_0=23769.70*u.MHz
m=b_0**2/(a_0*c_0)
Tbg=2.7355*u.K
cdms_catdir=CDMS.get_species_table()
jpl_catdir=JPLSpec.get_species_table()

trotdict={'SgrB2S':300*u.K,'DSi':300*u.K,'DSii':150*u.K,'DSiii':150*u.K,'DSiv':150*u.K,'DSv':100*u.K,'DSVI':300*u.K,'DSVII':200*u.K,'DSVIII':215*u.K,'DSIX':150*u.K,'DS10':150*u.K}

testT=trotdict[source]
#qrot_partfunc=partfunc(testT)#Q_rot_asym(testT).to('')

R_i=1
kappa=((2*b_0)-a_0-c_0)/(a_0-c_0)
f=1


dopplershifts={'SgrB2S':0.000228,'DSi':0.0001865,'DSii':0.000163,'DSiii':0.00017500261911843952,'DSiv':0.00018225233186845314,'DSv':0.0001838576164010067,'DSVI':0.0001661613132158407,'DSVII':0.00016320118280935546,'DSVIII':0.0001662062062062062,'DSIX':0.00015453732389175085,'DS10':0.00015794099431186572}#:0.000190713}/old doppler S: 0.0002306756533745274/0.00015954965399894244/0.00016236367659115043
pixdict={'SgrB2S':(70,59),'DSi':(36,42),'DSii':(22,24),'DSiii':(24,24),'DSiv':(32,31),'DSv':(19,19),'DSVI':(62,62),'DSVII':(75,75),'DSVIII':(50,50),'DSIX':(34,32),'DS10':(35,35)}#SgrB2S:61,64, DSIX:(34,35)
targetpix=pixdict[source]

s_othermol_dshift_v={' CH3CHO ':67.45330305*u.km/u.s,' C2H5OH ':67.45330305*u.km/u.s,' CH3OCHO ':67.45330305*u.km/u.s,' C(18)O ':69.551850256*u.km/u.s,' 13CH3OH ':67.5*u.km/u.s,' SO ':70.5*u.km/u.s}#' CH3OH ':68352.680424
ds2_othermol_dshift_v={' CH3OCHO ':49*u.km/u.s,' CH3CHO ':49*u.km/u.s,' C2H5OH ':49.3*u.km/u.s}#47831.782945392486 m / s
ds5_othermol_dshift_v={}
othermol_dopplershift={' CH3CHO ':0.000225,' C2H5OH ':0.000225,' CH3OCHO ':0.000225,' C(18)O ':0.000232}
ds9_othermol_dshift_v={}

sourceothers={'SgrB2S':s_othermol_dshift_v,'DSi':{},'DSii':ds2_othermol_dshift_v,'DSiii':{},'DSiv':{},'DSv':ds5_othermol_dshift_v,'DSVI':{},'DSVII':{},'DSVIII':{},'DSIX':ds9_othermol_dshift_v,'DS10':{}}
othermol_dshift_v=sourceothers[source]

z=dopplershifts[source]
z_vel=z*c

assemblecubepath=minicube_base+f'{source}/'+minicube_end[fnum]+'*spw*.fits'
print(f'Collecting spectra from {assemblecubepath}')
spectra=glob.glob(assemblecubepath)
spectra.sort()

homedict={'SgrB2S':'/blue/adamginsburg/d.jeff/XCLASS2021/files/SgrB2S/OctReimage_K/','DSi':'/blue/adamginsburg/d.jeff/XCLASS2021/files/DSi/field10originals_K/','DSii':'/blue/adamginsburg/d.jeff/XCLASS2021/files/DSii/field10originals_K/','DSiii':'/aug2023qrotfix/','DSiv':'/aug2023qrotfix/','DSv':f'/aug2023qrotfix/','DSVI':'/aug2023qrotfix/','DSVII':f'/aug2023qrotfix/','DSVIII':f'/aug2023qrotfix/','DSIX':f'/aug2023qrotfix/'}
sourcepath=ch3oh_sourcedict[source]
origsourcepath=f'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field{fnum}/CH3OH/{source}{sourcepath}'

texmappath=origsourcepath+'bootstrap_texmap_3sigma_allspw_withnans_weighted.fits'
    
fwhmpath=glob.glob(origsourcepath+'*fwhm*')[0]
nch3ohpath=glob.glob(origsourcepath+'ntotmap_allspw_withnans_weighted_useintercept_3sigma.fits')[0]
contpath=origsourcepath+'reprojectedcontinuum.fits'

texmapdata=fits.getdata(texmappath)*u.K
fwhmmap=fits.getdata(fwhmpath)*u.km/u.s
nch3ohmap=fits.getdata(nch3ohpath)*u.cm**-2

testT=texmapdata[targetpix[0],targetpix[1]]#350*u.K
fwhm_at_pix=fwhmmap[targetpix[0],targetpix[1]]
nch3oh_at_pix=nch3ohmap[targetpix[0],targetpix[1]]
#print(fwhm_at_pix)

reprojcontfits=fits.open(contpath)
reprojcont=reprojcontfits[0].data*u.Jy
reprojcontrestfreq=225*u.GHz#manual addition 11/9/2022, wiggle room w/i GHz
cntmbeam=radio_beam.Beam.from_fits_header(reprojcontfits[0].header)
reprojcont_K=reprojcont.to('K',cntmbeam.jtok_equiv(reprojcontrestfreq))
continuumlevel=reprojcont_K[targetpix[0],targetpix[1]]

stds=glob.glob(origsourcepath+'errorimgs/std/*.fits')
stds.sort()
images=['spw0','spw1','spw2','spw3']
    
assert 'spw0' in spectra[0] and 'spw0' in stds[0], 'List out of order'

print('Spectra and stds are sequential order')
linewidth=fwhm_at_pix#2.5*u.km/u.s#2.5 km/s is ideal for DSVI
print(f'Absolute model line width: {linewidth}\n')

sourcecolumns={'SgrB2S':sgrb2scolumns,'DSi':dsicolumns, 'DSii':ds2columns,'DSiii':ds3columns,'DSiv':ds4columns,
               'DSv':ds5columns,'DSVI':ds6columns,'DSVII':ds7columns,'DSVIII':ds8columns,'DSIX':ds9columns,'DS10':ds10columns}

columndict=sourcecolumns[source]
#plt.rcParams['figure.dpi'] = 300
#plt.figure(1, figsize=(30,10))
molcolors=['red','cyan','orange','brown','deepskyblue','darkviolet','yellow','pink','gold','darkkhaki','silver','blue','lime','blue','grey','plum','fuchsia','darkcyan','magenta','deeppink','gold','palegreen','goldenrod','indigo','dodgerblue','mediumpurple','yellow','red','grey']
spwmoldict={}
dummylist=[]
p1firstmolline={}#list(np.ones(len(columndict.keys())))
p2firstmolline={}
p3firstmolline={}
p4firstmolline={}
firstmolline=p1firstmolline
plotspecpad=0.005*u.GHz
n=1

cdms_catdir_qrot_temps=np.array([1000,500,300,225,150,75,37.5,18.75,9.375,5,2.725])
jpl_catdir_qrot_temps=[300, 225, 150, 75, 37.5, 18.75, 9.375]

aceswidth=[[85.965625,86.434375],[86.665625,87.134375],[89.159231,89.217821],[87.895942,87.954532],[97.6625,99.5375],[99.5625,101.4375]]

linedetections={}
for m in columndict.keys():
    p1firstmolline.update({m:1})
    p2firstmolline.update({m:1})
    p3firstmolline.update({m:1})
    p4firstmolline.update({m:1})

for spectrum, img, stdimage in zip(spectra,images,stds):
    
    print('Getting ready - '+img)
    plt.rcParams['figure.dpi'] = mode
    plt.figure(figsize=(20,10))
    n+=1
    if img == 'spw1':
        firstmolline=p2firstmolline
    if img == 'spw2':
        firstmolline=p3firstmolline
    if img == 'spw3':
        firstmolline=p4firstmolline
        '''
        plt.xlim(xmin=(p1minfreq-plotspecpad).value,xmax=(p1maxfreq+plotspecpad).value)
        plt.xlabel(r'$\nu$ (Hz)',fontsize=16)
        plt.ylabel('T$_b$ (K)',fontsize=16)
        plt.ylim(ymax=100)
        plt.tick_params(labelsize=13)
        plt.legend()
        plt.tight_layout()
        plt.show()
        '''
    
    cube=sc.read(spectrum)
    freqs=cube.spectral_axis
    data=cube[:,targetpix[0],targetpix[1]]+continuumlevel
    error=fits.getdata(stdimage)[targetpix[0],targetpix[1]]*u.K

    #freqs=(spec[:,0]*u.MHz).to('GHz')#cube.spectral_axis
    #data=spec[:,1]*u.K
    freqflip=False
    if freqs[0] > freqs[1]:
        freqs=freqs[::-1]
        data=data[::-1]
        freqflip=True
        print('Corrected decreasing frequency axis')
    else:
        pass
    
    freq_min=freqs[0]#*(1+z)#215*u.GHz
    freq_max=freqs[(len(freqs)-1)]#*(1+z)#235*u.GHz
    
    assert freq_max > freq_min, 'Decreasing frequency axis'
    
    print('Plotting model spectra')
    plt.plot(freqs.to('GHz').value,data.value,drawstyle='steps-mid',color='black')
    
    '''Generate methanol table for use during contaminant search'''
    #pdb.set_trace()
    Jfreqs, Jaij, Jdeg, JEU, qrot = get_molecular_parameters('CH3OH, vt=0-2', catalog='CDMS', fmin=freq_min, fmax=freq_max)
    
    qrot_partfunc=qrot(testT)
    print('Gathering CDMS table parameters')
    catdir_ch3oh=cdms_catdir[cdms_catdir['NAME'] == 'CH3OH, vt=0-2']
    mcatdir_qrot300=10**catdir_ch3oh['lg(Q(300))']
    methanol_table=CDMS.query_lines(min_frequency=freq_min,max_frequency=freq_max,min_strength=-500,molecule='032504 CH3OH, vt=0-2',get_query_payload=False)
    #Splatalogue.query_lines(freq_min, freq_max, chemical_name=' CH3OH ', energy_max=1840, energy_type='eu_k', line_lists=[linelistdict[' CH3OH ']], show_upper_degeneracy=True)
    mnus=methanol_table['FREQ']
    mlines=(mnus/(1+z)).to('GHz')#Redshifted to source
    
    melo_lambda=(1/methanol_table['ELO'].data)*u.cm
    melo_K=(((h*c)/melo_lambda)/k).to('K')
    melo_J=(melo_K*k).to('J')
    mdeltae=((h*methanol_table['FREQ'])/k).to('K')
    meuks=melo_K+mdeltae#maintable['EU_K']*u.K
    meujs=(meuks*k).to('J')
    mdegs=methanol_table['GUP']
    #Assembles the QNs for the molecule, not necessary in current implementation
    ju=methanol_table['Ju']
    jl=methanol_table['Jl']
    ku1=methanol_table['Ku']
    ku2=methanol_table['vu']
    kl1=methanol_table['Kl']
    kl2=methanol_table['vl']
    mqns=[]
    assert len(ju)==len(methanol_table) and len(jl)==len(methanol_table)
    for jupper,jlower,kupper1,kupper2,klower1,klower2 in zip(ju,jl,ku1,ku2,kl1,kl2):
        tempqn=f'{jupper}({kupper1},{kupper2})-{jlower}({klower1},{klower2})'
        mqns.append(tempqn)
    mlog10cdmsfluxes=methanol_table['LGINT']
    mcdmsfluxes=10**mlog10cdmsfluxes
    maijs=pickett_aul(mcdmsfluxes,mnus,mdegs,melo_J,meujs,mcatdir_qrot300,T=300*u.K)
    
    '''Create background model for the methanol lines and other species'''
    baseline=models.Linear1D(slope=(0*(u.K/u.Hz)),intercept=0*u.K)
    baseline.bounding_box=(freqs[0],freqs[(len(freqs)-1)])
    methmodelspec=baseline
    compositebaseline=baseline
    plot=np.linspace(freqs[0],freqs[(len(freqs)-1)],np.shape(data)[0]).to('GHz')
    modeldict={}
    #sys.exit()
    '''Generate species table for contaminant search'''
    molecule=' CH3OH '
    molname='CH3OH'
    cCfreqs, cCaij, cCdeg, cCEU, c_qrot = get_molecular_parameters(molname,catalog='JPL', fmin=freq_min, fmax=(freq_max+100*u.GHz),)
    species_catdir=jpl_catdir[jpl_catdir['NAME'] == molname]
    species_catdirtag=str(species_catdir['TAG'][0])
    scatdir_qrot300=10**species_catdir['QLOG1']
    jplname=f'{species_catdirtag} {molname}'
    species_table= JPLSpec.query_lines(min_frequency=freq_min,max_frequency=freq_max,min_strength=-500,molecule=jplname,get_query_payload=False)
    
    modelspec=baseline
    cnus=species_table['FREQ']
    celo_lambda=(1/species_table['ELO'].data)*u.cm
    celo_K=(((h*c)/celo_lambda)/k).to('K')
    celo_J=(celo_K*k).to('J')
    cdeltae=((h*species_table['FREQ'])/k).to('K')
    ceuks=celo_K+cdeltae#maintable['EU_K']*u.K
    ceujs=(ceuks*k).to('J')
    cdegs=species_table['GUP']
    clog10jplfluxes=species_table['LGINT']
    cjplfluxes=10**clog10jplfluxes
    caijs=pickett_aul(cjplfluxes,cnus,cdegs,celo_J,ceujs,scatdir_qrot300,T=300*u.K)
    cupqn=species_table['QN\'']
    clowqn=species_table['QN\"']
    cqns=[]
    for i, j in zip(cupqn,clowqn):
        i=i.replace(i[(len(i)-2):],'')
        i=i.replace(' ','.')
        j=j.replace(j[(len(j)-2):],'')
        j=j.replace(' ','.')
        a=f'{i}-{j}'
        cqns.append(a)
    c_qrot_partfunc=c_qrot(testT)
    clines=cnus/(1+z)
    
    cntot=nch3oh_at_pix#columndict[molecule]
        
    print(f'Begin model loops for JPL CH3OH')

    tempdetections={}
    jplbrights=[]
    for line,deg,euj,aij,qn in zip(clines,cdegs,ceujs,caijs,cqns):
        if isinstance(line,float):
            line=line*u.MHz
        restline=line*(1+z)#*u.MHz
        est_nupper=nupper_estimated(cntot,deg,c_qrot_partfunc,euj,testT).to('cm-2')
        modlinewidth=velocitytofreq(linewidth,line)#*u.MHz)
        lineprofilesigma=modlinewidth/2*np.sqrt(2*np.log(2))
        phi_nu=lineprofile(sigma=lineprofilesigma,nu_0=restline,nu=restline)
        intertau=lte_molecule.line_tau(testT, cntot, c_qrot_partfunc, deg, restline, euj, aij) 
        est_tau=(intertau*phi_nu).to('')#/modlinewidth
        trad=t_rad(tau_nu=est_tau,ff=f,nu=restline,T_ex=testT).to('K')
        #print(f'{qn} - {trad} - {np.log10(est_nupper.value)} - {deg} - {aij} - {euj} - {line} - {modlinewidth}')
        if trad >= 3*error:
            modelline=models.Gaussian1D(mean=line, stddev=modlinewidth, amplitude=trad)
            modelspec+=modelline
            compositebaseline+=modelline
            tempdetections.update({qn:True})
            jplbrights.append(trad.value)
        else:
            tempdetections.update({qn:False})
            continue
    linedetections.update({molecule:tempdetections})
    first=list(firstmolline.keys())[1]
    if firstmolline[first] and True in linedetections[molecule].values():
        plt.plot(freqs.to('GHz').value,modelspec(freqs.to('GHz')),drawstyle='steps-mid',color='purple',label='JPL')
        firstmolline[first]=0
        #print('hello')
        #sys.exit()
    elif True in linedetections[molecule].values():
        plt.plot(freqs.to('GHz').value,modelspec(freqs.to('GHz')),drawstyle='steps-mid',color=hue)
        #print('helloagain')
    cdmsbrights=[]
    if ' CH3OH ' in columndict:
        print('Begin CDMS CH3OH modeling\n')
        tempmdetections={}
        for line,deg,euj,aij,qn in zip(mlines,mdegs,meujs,maijs,mqns):
            #print(f'Transition: {qn} @ {line.to("GHz")}')
            restline=line*(1+z)
            modlinewidth=velocitytofreq(linewidth,line)
            lineprofilesigma=modlinewidth/2*np.sqrt(2*np.log(2))
            phi_nu=lineprofile(sigma=lineprofilesigma,nu_0=restline,nu=restline)
            
            methntot=nch3oh_at_pix#columndict[' CH3OH ']
            est_nupper=nupper_estimated(methntot,deg,qrot_partfunc,euj,testT).to('cm-2')
            intertau=lte_molecule.line_tau(testT, methntot, qrot_partfunc, deg, restline, euj, aij)
            est_tau=(intertau*phi_nu).to('')
            #print(f'Estimated tau: {"{:.3f}".format(est_tau)}')
            trad=t_rad(tau_nu=est_tau,ff=f,nu=restline,T_ex=testT).to('K')
            if trad >= 3*error:
                #print(f'Estimated brightness: {"{:.3f}".format(trad)}')
                #print(f'Model linewidth (Hz): {modlinewidth}')
                modelline=models.Gaussian1D(mean=line, stddev=modlinewidth, amplitude=trad)
                #modelgaus+=modelline
                methmodelspec+=modelline
                compositebaseline+=modelline
                tempmdetections.update({qn:True})
                cdmsbrights.append(trad.value)
            #elif qn == '16(6)-16(7)E1vt=1':
                #pdb.set_trace()
            else:
                #print(f'{qn} line brightness ({"{:.3f}".format(trad)}) below 3sigma threshold ({3*error})')
                tempmdetections.update({qn:False})
                continue
        linedetections.update({' CH3OH ':tempmdetections})
        if firstmolline[' CH3OH ']:
            plt.plot(freqs.to('GHz').value,methmodelspec(freqs.to('GHz')),drawstyle='steps-mid',linestyle='--',color='green',label='CDMS')
            firstmolline[' CH3OH ']=0
            print('yay')
        else:
            plt.plot(freqs.to('GHz').value,methmodelspec(freqs.to('GHz')),drawstyle='steps-mid',linestyle='--',color='green')
            print('yayy')
    
    '''
    print('Overplotting axvlines and transition annotations')
    for line,qn,detected in zip(clines,cqns,linedetections):
        if detected:
            plt.axvline(x=line.value,linestyle='--',color='yellow',ymin=0.25)
            plt.annotate(qn, (line.value, 0), (line.value-0.002,40),rotation=90)
        else:
            continue
    
    for mline,mqn,detected in zip(mlines,mqns,mdetections):
        if detected:
            plt.axvline(x=mline.value,linestyle='--',color='pink',ymin=0.25)
            plt.annotate(mqn, (mline.value, 0), (line.value-0.002,40),rotation=90)
        else:
            continue
    '''
    '''
    if img == 'spw3':
        p2maxfreq=max(freqs)
        plt.xlim(xmin=(p2minfreq-plotspecpad).value,xmax=(p2maxfreq+plotspecpad).value)
    '''
    if img == '':
        plt.xlabel(r'$\nu$ (GHz)',fontsize=16)
        plt.ylabel('T$_b$ (K)',fontsize=16)
        plt.xlim(xmin=231.7,xmax=232)
        plt.ylim(ymax=100,ymin=-20)
        plt.tick_params(labelsize=13)
        plt.tight_layout()
        plt.legend()
        #plt.savefig(fr'C:/Users/desmond/Desktop/submission_bootstrap_{source}_{img}_qrotfix_compositespectra.pdf')
        plt.show()
        sys.exit()
    else:
        plt.xlabel(r'$\nu$ (GHz)',fontsize=16)
        plt.ylabel('T$_b$ (K)',fontsize=16)
        plt.xlim(xmin=(min(freqs.to('GHz'))-plotspecpad).value,xmax=(max(freqs.to('GHz'))+plotspecpad).value)
        plt.ylim(ymax=100)
        plt.tick_params(labelsize=13)
        plt.tight_layout()
        plt.legend()
        #plt.savefig(f'../plots/HotCoreSpectra/{source}_{img}_individualizedspectra.pdf')
        plt.show()
        
        #plt.rcParams['figure.dpi'] = mode
        plt.figure(figsize=(20,10))
        plt.plot(freqs.to('GHz').value,data.value,drawstyle='steps-mid',color='black',label='Data')
        plt.plot(freqs.to('GHz').value,compositebaseline(freqs.to('GHz')),drawstyle='steps-mid',color='red',label='Composite Model')
        plt.xlabel(r'$\nu$ (GHz)',fontsize=16)
        plt.ylabel('T$_b$ (K)',fontsize=16)
        plt.xlim(xmin=(min(freqs.to('GHz'))-plotspecpad).value,xmax=(max(freqs.to('GHz'))+plotspecpad).value)
        plt.ylim(ymax=100)
        plt.tick_params(labelsize=13)
        plt.tight_layout()
        plt.legend()
        #plt.savefig(f'../plots/HotCoreSpectra/{source}_{img}_compositespectra.pdf')
        plt.show()
        if len(jplbrights)==len(cdmsbrights):
            print(f'JPL/CDMS brightness ratios{np.array(jplbrights)/np.array(cdmsbrights)}')