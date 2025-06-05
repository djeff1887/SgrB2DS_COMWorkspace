import astropy.units as u
# This file contains the parameters for the C2H5OH analysis in the SgrB2DS project.
# The parameters include source names, source locations, Doppler shifts, representative lines, representative cubes,
# excluded lines, and rest frequencies of the representative line.
catdirtag=32504

rotationalconstants=[24679.98*u.MHz,127484*u.MHz,23769.70*u.MHz]#B,A,C

sourcelocs = {'SgrB2S':'/sep2023-5removelasttorsionalline/','DSi':'/sep2023-5addvt2linesbackin/',
              'DSii':'/sep2023-2widerrefslab/','DSiii':'/sep2023-3vt2doublet/',
              'DSiv':'/sep2023-4nextinline/','DSv':f'/sep2023phi_nu&doublet/',
              'DSVI':'/sep2023-2removenewvt1line/','DSVII':f'/sep2023phi_nu&doublet/',
              'DSVIII':f'/sep2023phi_nu&doublet/','DSIX':f'/sep2023phi_nu&doublet/'}

dopplershifts = {
    'SgrB2S': 0.00023099669803283718,
    'DSi': 0.00018761288466593936,
    'DSii': 0.00016236367659115043,
    'DSiii': 0.000176,
    'DSiv': 0.00018225233186845314,
    'DSv': 0.0001838576164010067,
    'DSVI': 0.0001661613132158407,
    'DSVII': 0.00016236727257136008,
    'DSVIII': 0.0001661546432045067,
    'DSIX': 0.00015787296484373237
}

representativelines = {
    'SgrB2S': '20_1-20_0vt=0',
    'DSi': '8_1-7_0vt=0',
    'DSii': '8_1-7_0vt=0',
    'DSiii': '10_2--9_3-vt=0',
    'DSiv': '20_1-20_0vt=0',
    'DSv': '8_1-7_0vt=0',
    'DSVI': '8_1-7_0vt=0',
    'DSVII': '8_1-7_0vt=0',
    'DSVIII': '8_1-7_0vt=0',
    'DSIX': '8_1-7_0vt=0',
    'DS10': '10_2--9_3-vt=0',
    'DS11': '8_1-7_0vt=0',
    'DSX': '8_1-7_0vt=0'
}

representativecubes={'SgrB2S':0,'DSi':1,'DSii':1,'DSiii':2,'DSiv':0,'DSv':1,'DSVI':1,'DSVII':1,
                     'DSVIII':1,'DSIX':1,'DS10':2,'DS11':1,'DSX':1}#spwnumber


excludedlines = {
    'SgrB2S': [
        '7_6-7_7E1vt1', '14_6-14_7E1vt1', '11_6-11_7E1vt1', '15_6-15_7E1vt1',
        '9_6-9_7E1vt1', '13_6-13_7E1vt1', '12_6-12_7E1vt1', '8_6-8_7E1vt1',
        '16_6-16_7E1vt1', '10_6-10_7E1vt1'
    ],
    'DSi': [
        '11_6-11_7E1vt1', '25_3-24_4E1vt0', '23_5-22_6E1vt0', '14_6-14_7E1vt1',
        '7_6-7_7E1vt1', '15_6-15_7E1vt1', '16_6-16_7E1vt1', '9_6-9_7E1vt1',
        '10_6-10_7E1vt1', '11_6-11_7E1vt1', '12_6-12_7E1vt1', '13_6-13_7E1vt1'
    ],
    'DSii': [
        '7_6-7_7E1vt1', '9_6-9_7E1vt1', '14_6-14_7E1vt1', '10_6-10_7E1vt1',
        '13_6-13_7E1vt1', '11_6-11_7E1vt1', '23_5-22_6E1vt0'
    ],
    'DSiii': '',
    'DSiv': [
        '8_6-8_7E1vt1', '7_6-7_7E1vt1', '9_6-9_7E1vt1', '10_6-10_7E1vt1',
        '11_6-11_7E1vt1', '12_6-12_7E1vt1', '13_6-13_7E1vt1', '14_6-14_7E1vt1',
        '6_1--7_2-vt1'
    ],
    'DSv': '',
    'DSVI': [
        "6_1--7_2-vt1", '14_6-14_7E1vt1', '10_6-10_7E1vt1', '9_6-9_7E1vt1',
        '11_6-11_7E1vt1', '13_6-13_7E1vt1', '12_6-12_7E1vt1', '13_3--14_4-vt2',
        '13_3+-14_4+vt2', '7_6-7_7E1vt1', '16_6-16_7E1vt1', '8_6-8_7E1vt1',
        '17_6-17_7E1vt1'
    ],
    'DSVII': ["6_1--7_2-vt1"],
    'DSVIII': '',
    'DSIX': ''
}

restfreq_representativeline={'SgrB2S':217.88650400*u.GHz,'DSi':220.07856100*u.GHz,
                             'DSii':220.07856100*u.GHz,'DSiii':231.28111000*u.GHz,
                             'DSiv':217.88650400*u.GHz,'DSv':220.07856100*u.GHz,
                             'DSVI':220.07856100*u.GHz,'DSVII':220.07856100*u.GHz,
                             'DSVIII':220.07856100*u.GHz,'DSIX':220.07856100*u.GHz,
                             'DS10':231.28111000*u.GHz,'DS11':220.07856100*u.GHz,
                             'DSX':220.07856100*u.GHz}#All taken from Splatalogue;  oldS 218.44006300#All taken from Splatalogue
