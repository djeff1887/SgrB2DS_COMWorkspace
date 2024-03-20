import astropy.units as u
import scipy.constants as cnst
import numpy as np
from pyspeckit.spectrum.models.lte_molecule import get_molecular_parameters
import os
from astroquery.splatalogue import Splatalogue
from astroquery.linelists.cdms import CDMS
from astroquery.jplspec import JPLSpec
from astropy.modeling.functional_models import Linear1D
from astropy.modeling import fitting

c=cnst.c*u.m/u.s
k=cnst.k*u.J/u.K
h=cnst.h*u.J*u.s
sigma_sb=cnst.sigma*u.W/((u.m)**(2)*(u.K)**(4))
R_i=1
f=1
Tbg=2.7355*u.K
mu_a=(0.896e-18*u.statC*u.cm).to('cm(3/2) g(1/2) s-1 cm')

incompleteqrot=[' C2H5OH ']

linelistdict={' CH3OH ':'JPL',' CH3OCHO ':'JPL',' CH3CHO ':'JPL',' C2H5OH ':'CDMS',' CH3OCH3 ':'JPL',' DCN ':'JPL',' OCS ':'CDMS',' 13CH3OH ':'CDMS',' H2CO ':'CDMS',' HC3N ':'CDMS',' C(18)O ':'CDMS',' 13CS ':'CDMS',' SO2 ':'CDMS',' NH2CHO ':'JPL',' HNCO ':'CDMS',' SO ':'CDMS', ' SiO ':'CDMS',' H2S ':'CDMS',' c-HCCCH ':'CDMS', ' HC3N v7=1 ':'CDMS',' H213CO ':'CDMS',' 13CH3CN ':'CDMS',' CH3COOH ':'CDMS',' t-HCOOH ':'CDMS',' CH3O13CHO ':'TopModel',' HNO3 ':'JPL',' CH3O13CHO, vt = 0, 1':'CDMS','NH2CN':'JPL','CH2CHCN':'CDMS','CH3OCHO v=1':'JPL','18OCS':'CDMS','CH3NCO':'CDMS','CH3CH2CN':'CDMS','NH2CO2CH3 v=1 ':'JPL',' HOCN ':'CDMS',' 13CH3CCH ':'JPL'}

jplnamelist={' CH3OCHO ':'CH3OCHO',' CH3CHO ':'CH3CHO',' CH3OH ':'CH3OH'}
cdmsnamelist={' 13CH3OH ':'C-13-H3OH, vt=0,1',' C(18)O ':'CO-18',' 13CS ':'C-13-S, v=0,1',' NH2CHO ':'HC(O)NH2, v=0',' c-HCCCH ':'c-C3H2',' HC3N v7=1 ':'HC3N, v7=1',' H213CO ':'H2C-13-O',' 13CH3CN ':'C-13-H3CN, v=0',' 18OCS ':'O-18-CS',' N17O ':'N-15-O-17',' CH3CH2CN ':'C2H5CN, v=0', ' 13CH3OCH3 ':'C-13-H3OCH3',}

cdmsproblemchildren=[' OCS ',' 13CS ',' C(18)O ',' HNCO ',' SO ',' HC3N ',' CH3NCO, vb=0 ',' CH3CH2CN ',' 13CH3CN ']
problemchildren2=[' CH3NCO, vb=0 ']

'''Hot core locations'''

fields={'SgrB2S':1,'DSi':10,'DSii':10,'DSiii':10,'DSiv':10,'DSv':10,'DSVI':2,'DSVII':3,'DSVIII':3,'DSIX':7}
orange_stacontsubpaths={1:'/orange/adamginsburg/sgrb2/d.jeff/data/OctReimage_K/',10:"/orange/adamginsburg/sgrb2/d.jeff/data/field10originals_K/",2:"/orange/adamginsburg/sgrb2/d.jeff/data/field2originals_K/",3:"/orange/adamginsburg/sgrb2/d.jeff/data/field3originals_K/",7:"/orange/adamginsburg/sgrb2/d.jeff/data/field7originals_K/",8:"/orange/adamginsburg/sgrb2/d.jeff/data/field8originals_K/"}
orange_stdhomedict={1:'/orange/adamginsburg/sgrb2/d.jeff/products/OctReimage_K/',10:'/orange/adamginsburg/sgrb2/d.jeff/products/field10originals_K/',2:'/orange/adamginsburg/sgrb2/d.jeff/products/field2originals_K/',3:'/orange/adamginsburg/sgrb2/d.jeff/products/field3originals_K/',7:'/orange/adamginsburg/sgrb2/d.jeff/products/field7originals_K/'}
minicube_base=f'/orange/adamginsburg/sgrb2/2017.1.00114.S/desmond/SgrB2DSminicubes/'
minicube_end={1:'OctReimage_K/',10:'field10originals_K/',2:'field2originals_K/',3:'field3originals_K/',7:'field7originals_K/',8:'field8originals_K/'}

dataversion={'ch3oh':'pacman_sep2023revolution'}
datadirs={'ch3oh':f'/blue/adamginsburg/d.jeff/imaging_results/SgrB2DS-CH3OH/{dataversion}/'}

rotationalconstants={' CH3OH ':[24679.98*u.MHz,127484*u.MHz,23769.70*u.MHz],' C2H5OH ':[(34.89170*u.GHz).to('MHz'),(9.35065*u.GHz).to('MHz'),(8.13520*u.GHz).to('MHz')]}#B,A,C


ch3oh_sourcedict={'SgrB2S':'/sep2023-5removelasttorsionalline/','DSi':'/sep2023-5addvt2linesbackin/','DSii':'/sep2023-2widerrefslab/','DSiii':'/sep2023-3vt2doublet/','DSiv':'/sep2023-4nextinline/','DSv':f'/sep2023phi_nu&doublet/','DSVI':'/sep2023-2removenewvt1line/','DSVII':f'/sep2023phi_nu&doublet/','DSVIII':f'/sep2023phi_nu&doublet/','DSIX':f'/sep2023phi_nu&doublet/'}
ch3oh_excludedlines={'SgrB2S':['7_6-7_7E1vt1','14_6-14_7E1vt1','11_6-11_7E1vt1','15_6-15_7E1vt1','9_6-9_7E1vt1','13_6-13_7E1vt1','12_6-12_7E1vt1','8_6-8_7E1vt1','16_6-16_7E1vt1','10_6-10_7E1vt1'],'DSi':['11_6-11_7E1vt1','25_3-24_4E1vt0','23_5-22_6E1vt0','14_6-14_7E1vt1','7_6-7_7E1vt1','15_6-15_7E1vt1','16_6-16_7E1vt1','9_6-9_7E1vt1','10_6-10_7E1vt1','11_6-11_7E1vt1','12_6-12_7E1vt1','13_6-13_7E1vt1'],'DSii':['7_6-7_7E1vt1','9_6-9_7E1vt1','14_6-14_7E1vt1','10_6-10_7E1vt1','13_6-13_7E1vt1','11_6-11_7E1vt1','23_5-22_6E1vt0'],'DSiii':'','DSiv':['8_6-8_7E1vt1','7_6-7_7E1vt1','9_6-9_7E1vt1','10_6-10_7E1vt1','11_6-11_7E1vt1','12_6-12_7E1vt1','13_6-13_7E1vt1','14_6-14_7E1vt1','6_1--7_2-vt1'],'DSv':'','DSVI':["6_1--7_2-vt1",'14_6-14_7E1vt1','10_6-10_7E1vt1','9_6-9_7E1vt1','11_6-11_7E1vt1','13_6-13_7E1vt1','12_6-12_7E1vt1','13_3--14_4-vt2','13_3+-14_4+vt2','7_6-7_7E1vt1','16_6-16_7E1vt1','8_6-8_7E1vt1','17_6-17_7E1vt1'],'DSVII':["6_1--7_2-vt1"],'DSVIII':'','DSIX':''}


c2h5oh_dopplershifts={'DSi':(55.3*u.km/u.s)}
ch3oh_dopplershifts={'SgrB2S':0.00023099669803283718,'DSi':0.00018761288466593936,'DSii':0.00016236367659115043,'DSiii':0.000176,'DSiv':0.00018225233186845314,'DSv':0.0001838576164010067,'DSVI':0.0001661613132158407,'DSVII':0.00016236727257136008,'DSVIII':0.0001661546432045067,'DSIX':0.00015787296484373237}

masterdopplershifts={' CH3OH ':ch3oh_dopplershifts,' C2H5OH ':c2h5oh_dopplershifts}

excludedlines={'ch3oh':ch3oh_excludedlines}
dopplershifts={'ch3oh':ch3oh_dopplershifts}
sourcedict={'ch3oh':ch3oh_sourcedict}

targetworldcrds={'SgrB2S':[[0,0,0],[266.8351718,-28.3961210, 0]], 'DSi':[[0,0,0],[266.8316149,-28.3972040,0]], 'DSii':[[0,0,0],[266.8335363,-28.3963158,0]],'DSiii':[[0,0,0],[266.8332758,-28.3969269,0]],'DSiv':[[0,0,0],[266.8323834, -28.3954424,0]],'DSv':[[0,0,0],[266.8321331, -28.3976585, 0]],'DSVI':[[0,0,0],[266.8380037, -28.4050741,0]],'DSVII':[[0,0,0],[266.8426074, -28.4094401,0]],'DSVIII':[[0,0,0],[266.8418408, -28.4118242, 0]],'DSIX':[[0,0,0],[266.8477371, -28.4311386,0]],'DS10':[[0,0,0],[266.8373798, -28.4009340,0]],'DS11':[[0,0,0],[266.8374572, -28.3996894, 0]],'DSX':[[0,0,0],[266.8452950, -28.4282608,0]]}
pixdict={'SgrB2S':(73,54),'DSi':(36,40),'DSii':(22,24),'DSiii':(24,24),'DSiv':(32,31),'DSv':(19,19),'DSVI':(62,62),'DSVII':(75,75),'DSVIII':(50,50),'DSIX':(34,35)}
tblconversion={'SgrB2S':'SgrB2S','DSi':'DS1','DSii':'DS2','DSiii':'DS3','DSiv':'DS4','DSv':'DS5','DSVI':'DS6','DSVII':'DS7','DSVIII':'DS8','DSIX':'DS9'}

'''Line Queries'''
def cdms_get_molecule_name(my_molecule_name, **kwargs):
    basename = dict(CDMS.query_lines(min_frequency=1*u.GHz, max_frequency=500*u.GHz, molecule=my_molecule_name, parse_name_locally=True, get_query_payload=True, **kwargs))['Molecules']
    return " ".join(basename.split(" ")[1:])

def jpl_get_molecule_name(my_molecule_name, **kwargs):
    basename = dict(JPL.query_lines(min_frequency=1*u.GHz, max_frequency=500*u.GHz, molecule=my_molecule_name, parse_name_locally=True, get_query_payload=True, **kwargs))['Molecules']
    return " ".join(basename.split(" ")[1:])

def pickett_aul(intensity,nu,g,elower,eupper,q,T=300*u.K):
    return (intensity*(nu**2)*(q/g)*((np.exp(-elower/(k*T))-np.exp(-eupper/(k*T)))**-1)*2.7964e-16).value*u.Hz

'''LTE Analysis'''
#def Tbthick(ntot,nu,line_width,mulu_2,g,q,eu_J,T_ex):
#    print(f'ntot: {ntot}, nu: {nu},line_width: {line_width},mulu_2: {mulu_2},g: {g},q: {q},eu_J: {eu_J},T_ex: {T_ex}')
#    return (1-np.exp(((-8*np.pi**3*mulu_2*R_i*g)/(3*h*q*line_width))*((np.exp((h*nu)/(k*T_ex))-1)/np.exp((eu_J)/(k*T_ex)))*ntot))*(f*(rjequivtemp(nu,T_ex)-rjequivtemp(nu,Tbg)))
    
def mulu(aij,nu):#Rearranged from Eq 11 (Magnum & Shirley 2015), returns product in units of cm5 g s-2
    return (3*h*c**3*aij)/(64*np.pi**4*nu**3)

'''Converts given line list in frequency to radio velocity'''
def vradio(frequency,rest_freq):
    velocity=c.to(u.km/u.s)*(1-((rest_freq-frequency)/rest_freq))
    return velocity.to('km s-1')
    
'''Converts given velocity width to frequency width'''
def velocitytofreq(velocity,ref_freq):
    frequency=((velocity)/c.to(u.km/u.s))*ref_freq
    return frequency.to('MHz')
    
def freqtovelocity(freq,ref_freq):
    velocity=(freq/ref_freq)*c.to('km s-1')
    return velocity.to('km s-1')

'''Compute Rayleigh-Jean equivalent brightness temperature per M&S 2015 Eq. 24'''
def rjequivtemp(nu,T_ex):
    return ((h*nu)/k)/(np.exp((h*nu)/(k*T_ex))-1)
    
'''Compute upper-state column density and associated error of transition of interest'''      
def N_u(nu,Aij,velocityintegrated_intensity_K,velint_intK_err):#(ntot,qrot,gu,eu_J,T_ex):#taken from pyspeckit documentation https://pyspeckit.readthedocs.io/en/latest/lte_molecule_model.html?highlight=Aij#lte-molecule-model
    nuppercalc=((8*np.pi*k*nu**2)/(h*c**3*Aij))*velocityintegrated_intensity_K
    nuppererr=((8*np.pi*k*nu**2)/(h*c**3*Aij))*velint_intK_err#(velint_intK_err.to('K km s-1')/velocityintegrated_intensity_K.to('K km s-1'))*nuppercalc
    return nuppercalc.to('cm-2'),nuppererr.to('cm-2')#ntot/(qrot*np.exp(eu_J/(k*T_ex)))
    
def KtoJ(T):#Convert from excitation temperature (Kelvin) to energy (Joules)
    return k*T

def JtoK(E):
    return (E/k).to('K')
    
def qngrabber(nums):
    temp=nums.split('(')
    temp2=temp[1].split(',')
    jupper=int(temp[0])
    if linelist == 'JPL':
        temp3=temp2[0].split(')')
        kupper=temp3[0]
        if 'a' in kupper:#What are these things?
            kupper=0
        else:
            kupper=int(temp3[0])
    else:
        kupper=int(temp2[0])
        
    return jupper, kupper
        
def t_rad(tau_nu, ff, nu, T_ex):#Radiation temperature per M&S (2015) Eq 27
    return ff*(1-np.exp(-tau_nu))*(rjequivtemp(nu, T_ex)-rjequivtemp(nu,Tbg))
    
def nupper_estimated(n_tot,g,q,euj,tex):#LTE model upper-state column density (nupper) using a fiducial/test total column density (ntot) and excitation temperature (Tex)
    return n_tot*(g/q)*np.exp(-euj/(k*tex))
    
def opticaldepth(aij,nu,T_ex,nupper,lw):#LTE model optical depth using nupper from the above equation
    return (c**2/(8*np.pi*nu**2*lw))*aij*nupper*np.exp((h*nu)/(k*T_ex))
    
'''Replaces unwanted characters from the QN table for use in file names'''
def qn_replace(string):
    string=string.replace('=','')
    string=string.replace('(','_')
    string=string.replace(')','')
    string=string.replace(',','&')
    return string
    
def lineprofile(sigma,nu_0,nu):
    return (1/(np.sqrt(2*np.pi)*sigma))*np.exp(-(nu-nu_0)**2/(2*sigma**2))

def Ctau(tau):
    return tau/(1-np.exp(-tau))

def fit_qrot(input_logintercept=0,input_temperature=301*u.K,input_powerlawfit=0):
    return input_logintercept*input_temperature.value**input_powerlawfit.slope
