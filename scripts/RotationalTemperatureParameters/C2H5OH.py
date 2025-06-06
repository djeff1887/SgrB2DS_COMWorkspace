import astropy.units as u
# This file contains the parameters for the C2H5OH analysis in the SgrB2DS project.
# The parameters include source names, source locations, Doppler shifts, representative lines, representative cubes,
# excluded lines, and rest frequencies of the representative line.
catdirtag=46524

rotationalconstants=[(34.89170*u.GHz).to('MHz'),(9.35065*u.GHz).to('MHz'),(8.13520*u.GHz).to('MHz')]

sourcelocs = {
    'SgrB2S': '/mar2025_2_removesDS2exclusions/',
    'DSi': '/COPILOTREFACTORING_first_attempt_at_generalizing_script/',#'/oct2024_1_removesDS2exclusions/',
    'DSii': '/oct2024_1_removeproblemlines/',
    'DSiii': '/dec2024_3_try-close-to-FWZI/',
    'DSiv': '/nov2024_1_firstrun_removesDS2exclusions/',
    'DSv': '/mar2025_1_removesDS2exclusions/',
    'DSVI': '/nov2024_1_removesDS2exclusions/',
    'DSVII': '/apr2025_2_removesDS2exclusions/',
    'DSVIII': '/may2025_1_removesDS2exclusions/'
}

dopplershifts={'SgrB2S':(66.2*u.km/u.s),'DSi':(55.3*u.km/u.s),'DSii':(49.3*u.km/u.s),'DSiii':(52.267*u.km/u.s),
                      'DSiv':(54.648*u.km/u.s),'DSv':(56.7*u.km/u.s),'DSVI':(51*u.km/u.s),'DSVII':(47.22*u.km/u.s),
                      'DSVIII':(51.6379*u.km/u.s)}#all taken from peaks of representative lines

representativelines={'SgrB2S':'14.0.14_2-13.1.13_2','DSi':'14.0.14_2-13.1.13_2','DSii':'14.0.14_2-13.1.13_2','DSiii':'14.0.14_2-13.1.13_2',
                     'DSiii':'14.0.14_2-13.1.13_2','DSiv':'14.0.14_2-13.1.13_2','DSv':'16.5.11_2-16.4.12_2','DSVI':'14.0.14_2-13.1.13_2',
                     'DSVII':'14.0.14_2-13.1.13_2','DSVIII':'14.0.14_2-13.1.13_2'}

representativecubes={'SgrB2S':2,'DSi':2,'DSii':2,'DSiii':2,'DSiv':2,'DSv':2,'DSVI':2,'DSVII':2,'DSVIII':2,'DSIX':'','DSX':''}

excludedlines={'SgrB2S':['35(4,31)(2)-35(3,32)(2)','45(6,39)(2)-45(5,40)(2)','21(5,17)(2)-21(4,18)(2)','20(5,16)(2)-20(4,17)(2)','13(5,8)(2)-13(4,9)(2)','22(5,18)(2)-22(4,19)(2)'],
               'DSi':['13(5,8)(2)-13(4,9)(2)','21(5,17)(2)-21(4,18)(2)','20(5,16)(2)-20(4,17)(2)','22(5,18)(2)-22(4,19)(2)','13(2,11)(0)-12(2,10)(0)'],
               'DSii':['13(5,8)(2)-13(4,9)(2)','21(5,17)(2)-21(4,18)(2)','20(5,16)(2)-20(4,17)(2)','22(5,18)(2)-22(4,19)(2)'],
               'DSiii':['13(5,8)(2)-13(4,9)(2)','21(5,17)(2)-21(4,18)(2)','20(5,16)(2)-20(4,17)(2)','22(5,18)(2)-22(4,19)(2)','13(2,11)(0)-12(2,10)(0)'],
               'DSiv':['13(5,8)(2)-13(4,9)(2)','21(5,17)(2)-21(4,18)(2)','20(5,16)(2)-20(4,17)(2)','22(5,18)(2)-22(4,19)(2)','13(2,11)(0)-12(2,10)(0)'],
               'DSv':['13(5,8)(2)-13(4,9)(2)','21(5,17)(2)-21(4,18)(2)','20(5,16)(2)-20(4,17)(2)','22(5,18)(2)-22(4,19)(2)'],
               'DSVI':['13(5,8)(2)-13(4,9)(2)','21(5,17)(2)-21(4,18)(2)','20(5,16)(2)-20(4,17)(2)','22(5,18)(2)-22(4,19)(2)','13(2,11)(0)-12(2,10)(0)',],
              'DSVII':['13(5,8)(2)-13(4,9)(2)','21(5,17)(2)-21(4,18)(2)','20(5,16)(2)-20(4,17)(2)','22(5,18)(2)-22(4,19)(2)'],
              'DSVIII':['13(5,8)(2)-13(4,9)(2)','21(5,17)(2)-21(4,18)(2)','20(5,16)(2)-20(4,17)(2)','22(5,18)(2)-22(4,19)(2)']}# '35(4,31)(2)-35(3,32)(2)']}

restfreq_representativeline={'SgrB2S':230.9913834*u.GHz,'DSi':230.9913834*u.GHz,'DSii':230.9913834*u.GHz,
                             'DSiii':230.9913834*u.GHz,'DSiv':230.9913834*u.GHz,'DSv':230.9537832*u.GHz,
                             'DSVI':230.9913834*u.GHz,'DSVII':230.9913834*u.GHz,'DSVIII':230.9913834*u.GHz,
                             'DSIX':'','DSX':''}#All taken from Splatalogue
