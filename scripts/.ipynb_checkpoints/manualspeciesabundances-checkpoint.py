import astropy.units as u

sgrb2scolumns={
    ' CH3OH ':2.1e18*u.cm**-2,' CH3OCHO ':6e17*u.cm**-2, ' CH3CHO ':2e17*u.cm**-2,' C2H5OH ':2e18*u.cm**-2,' CH3OCH3 ':8e17*u.cm**-2,
    ' DCN ':5e15*u.cm**-2, ' OCS ':9e16*u.cm**-2,' 13CH3OH ':9e17*u.cm**-2,' H2CO ':7e17*u.cm**-2,' HC3N ':4e15*u.cm**-2, ' C(18)O ':1.3e18*u.cm**-2,
    ' 13CS ':6e15*u.cm**-2,' SO2 ':2e17*u.cm**-2,' NH2CHO ':5e16*u.cm**-2,' HNCO ':3e17*u.cm**-2,' SO ':9e16*u.cm**-2,' SiO ':2e15*u.cm**-2,
    ' H2S ':2e17*u.cm**-2,' c-HCCCH ':3e16*u.cm**-2, 'HC3N v7=1':5e15*u.cm**-2,' H213CO ':5e16*u.cm**-2,' 13CH3CN ':1e16*u.cm**-2,
    ' CH2CHCN ':3e16*u.cm**-2,' 18OCS ':1e17*u.cm**-2,' CH3NCO, vb = 0 ':5e17*u.cm**-2,' CH3CH2CN ':5e16*u.cm**-2, ' CH3COCH3 ':6e17*u.cm**-2,}#' 13CH3CCH ':1e17*u.cm**-2}#' NH2CN ':1e15*u.cm**-2,}#' 13CH3OCH3 ':1e14*u.cm**-2}#'CH3OCHO v=1':1e17*u.cm**-2}#' CH3O13CHO ':1e14*u.cm**-2,' H2CCO ':1e16*u.cm**-2,}#' H2CS ':1e18*u.cm**-2,' CH3(18)OH ':2.5e16*u.cm**-2,' NH2CN ':1e14*u.cm**-2, ' NH2D ':1e15*u.cm**-2, ' t-HCOOH ':5e14*u.cm**-2,' SiS ':1e14*u.cm**-2, ' c-HCCCH ':2.5e15*u.cm**-2, 'Acetone':6e13*u.cm**-2,' CH3C(15)N ':3e13*u.cm**-2,' SiN ':2e15*u.cm**-2, ' CH3NH2 ':9e15*u.cm**-2,}#' HOONO ':5e15*u.cm**-2,' CH3COOH ':2e15*u.cm**-2,
#CDMS - ' CH3OH ':1.2e17*u.cm**-2,

dsicolumns = {' CH3OH ': 1.48e18 * u.cm**-2, ' CH3OCHO ': 2e17 * u.cm**-2, ' CH3CHO ': 8e16 * u.cm**-2,
    ' C2H5OH ': 4e17 * u.cm**-2, ' CH3OCH3 ': 9e16 * u.cm**-2, ' DCN ': 6e14 * u.cm**-2,
    ' OCS ': 6e16 * u.cm**-2, ' 13CH3OH ': 4e17 * u.cm**-2, ' H2CO ': 9e16 * u.cm**-2,
    ' HC3N ': 2e15 * u.cm**-2, ' C(18)O ': 1e18 * u.cm**-2, ' 13CS ': 4e15 * u.cm**-2,
    ' SO2 ': 3e16 * u.cm**-2, ' NH2CHO ': 1e16 * u.cm**-2, ' HNCO ': 3e16 * u.cm**-2,
    ' SO ': 1e16 * u.cm**-2, ' SiO ': 1e15 * u.cm**-2, ' H2S ': 9e16 * u.cm**-2,
    ' c-HCCCH ': 5e15 * u.cm**-2, 'HC3N v7=1': 2e15 * u.cm**-2, ' H213CO ': 7e15 * u.cm**-2,
    ' 13CH3CN ': 4e15 * u.cm**-2, ' CH2CHCN ': 2e16 * u.cm**-2, ' 18OCS ': 5e15 * u.cm**-2,
    ' CH3NCO, vb = 0 ': 6e16 * u.cm**-2, ' CH3CH2CN ': 3e16 * u.cm**-2, ' CH3COCH3 ': 3e17 * u.cm**-2,
}

ds2columns = {' CH3OH ': 6e17 * u.cm**-2, ' CH3OCHO ': 6e16 * u.cm**-2, ' CH3CHO ': 3e16 * u.cm**-2,
    ' C2H5OH ': 9e16 * u.cm**-2, ' CH3OCH3 ': 1e17 * u.cm**-2, ' DCN ': 5e14 * u.cm**-2,
    ' OCS ': 3e16 * u.cm**-2, ' 13CH3OH ': 8e16 * u.cm**-2, ' H2CO ': 1e17 * u.cm**-2,
    ' HC3N ': 2e15 * u.cm**-2, ' C(18)O ': 6e17 * u.cm**-2, ' 13CS ': 1.5e15 * u.cm**-2,
    ' SO2 ': 2.5e16 * u.cm**-2, ' NH2CHO ': 5e15 * u.cm**-2, ' HNCO ': 2e16 * u.cm**-2,
    ' SO ': 1.2e16 * u.cm**-2, ' SiO ': 5e14 * u.cm**-2, ' H2S ': 6e16 * u.cm**-2,
    ' c-HCCCH ': 4e15 * u.cm**-2, 'HC3N v7=1': 3e15 * u.cm**-2, ' H213CO ': 7e15 * u.cm**-2,
    ' 13CH3CN ': 2e15 * u.cm**-2, ' CH2CHCN ': 8e15 * u.cm**-2, ' 18OCS ': 6e15 * u.cm**-2,
    ' CH3NCO, vb = 0 ': 1.5e16 * u.cm**-2, ' CH3CH2CN ': 7e15 * u.cm**-2, ' CH3COCH3 ': 5e16 * u.cm**-2,
}

ds3columns = {' CH3OH ': 7e17 * u.cm**-2, ' CH3OCHO ': 3e17 * u.cm**-2, ' CH3CHO ': 7e16 * u.cm**-2,
    ' C2H5OH ': 5e17 * u.cm**-2, ' CH3OCH3 ': 3e17 * u.cm**-2, ' DCN ': 9e14 * u.cm**-2,
    ' OCS ': 3.5e16 * u.cm**-2, ' 13CH3OH ': 3e17 * u.cm**-2, ' H2CO ': 8e16 * u.cm**-2,
    ' HC3N ': 3e15 * u.cm**-2, ' C(18)O ': 1e18 * u.cm**-2, ' 13CS ': 1e15 * u.cm**-2,
    ' SO2 ': 2.5e16 * u.cm**-2, ' NH2CHO ': 9e15 * u.cm**-2, ' HNCO ': 2e16 * u.cm**-2,
    ' SO ': 1e16 * u.cm**-2, ' SiO ': 7e14 * u.cm**-2, ' H2S ': 8e16 * u.cm**-2,
    ' c-HCCCH ': 2e16 * u.cm**-2, 'HC3N v7=1': 2e15 * u.cm**-2, ' H213CO ': 8e15 * u.cm**-2,
    ' 13CH3CN ': 3e15 * u.cm**-2, ' CH2CHCN ': 7e15 * u.cm**-2, ' 18OCS ': 2e16 * u.cm**-2,
    ' CH3CH2CN ': 9e15 * u.cm**-2, ' CH3COCH3 ': 4e17 * u.cm**-2,
}

ds4columns = {' CH3OH ': 9e17 * u.cm**-2, ' CH3OCHO ': 5e17 * u.cm**-2, ' CH3CHO ': 7e16 * u.cm**-2,
    ' C2H5OH ': 9e17 * u.cm**-2, ' CH3OCH3 ': 5e17 * u.cm**-2, ' DCN ': 2e15 * u.cm**-2,
    ' OCS ': 1e17 * u.cm**-2, ' 13CH3OH ': 5e17 * u.cm**-2, ' H2CO ': 4e17 * u.cm**-2,
    ' HC3N ': 3e15 * u.cm**-2, ' C(18)O ': 2e18 * u.cm**-2, ' 13CS ': 2e15 * u.cm**-2,
    ' SO2 ': 4e16 * u.cm**-2, ' NH2CHO ': 2e16 * u.cm**-2, ' HNCO ': 4e16 * u.cm**-2,
    ' SO ': 3e16 * u.cm**-2, ' SiO ': 5e14 * u.cm**-2, ' H2S ': 1.8e17 * u.cm**-2,
    ' c-HCCCH ': 2e16 * u.cm**-2, 'HC3N v7=1': 1.5e15 * u.cm**-2, ' H213CO ': 2e16 * u.cm**-2,
    ' 13CH3CN ': 8e15 * u.cm**-2, ' CH2CHCN ': 1.5e16 * u.cm**-2, ' 18OCS ': 3e16 * u.cm**-2,
    ' CH3NCO, vb = 0 ': 9e16 * u.cm**-2, ' CH3CH2CN ': 2e16 * u.cm**-2, ' CH3COCH3 ': 6e17 * u.cm**-2,
}

ds5columns = {
    ' CH3OH ': 1.8e17 * u.cm**-2, ' CH3OCHO ': 4e16 * u.cm**-2, ' CH3CHO ': 1e16 * u.cm**-2,
    ' C2H5OH ': 2e17 * u.cm**-2, ' CH3OCH3 ': 8e16 * u.cm**-2, ' DCN ': 5e14 * u.cm**-2,
    ' OCS ': 1e16 * u.cm**-2, ' 13CH3OH ': 1e17 * u.cm**-2, ' H2CO ': 1e17 * u.cm**-2,
    ' HC3N ': 8e14 * u.cm**-2, ' C(18)O ': 1e18 * u.cm**-2, ' 13CS ': 1e15 * u.cm**-2,
    ' SO2 ': 9e15 * u.cm**-2, ' NH2CHO ': 3e15 * u.cm**-2, ' HNCO ': 1e16 * u.cm**-2,
    ' SO ': 2e16 * u.cm**-2, ' SiO ': 5e14 * u.cm**-2, ' H2S ': 4e16 * u.cm**-2,
    ' c-HCCCH ': 5e14 * u.cm**-2, 'HC3N v7=1': 5e14 * u.cm**-2, ' H213CO ': 5e15 * u.cm**-2,
    ' 13CH3CN ': 8e14 * u.cm**-2, ' CH2CHCN ': 1e15 * u.cm**-2, ' 18OCS ': 6e15 * u.cm**-2,
    ' CH3CH2CN ': 5e15 * u.cm**-2,
}#' CH3NCO, vb = 0 ':1e16*u.cm**-2

ds6columns = {
    ' CH3OH ': 3e18 * u.cm**-2, ' CH3OCHO ': 4e17 * u.cm**-2, ' CH3CHO ': 1.3e17 * u.cm**-2,
    ' C2H5OH ': 8e17 * u.cm**-2, ' CH3OCH3 ': 4e17 * u.cm**-2, ' DCN ': 3.5e15 * u.cm**-2,
    ' OCS ': 7e16 * u.cm**-2, ' 13CH3OH ': 7e17 * u.cm**-2, ' H2CO ': 7e17 * u.cm**-2,
    ' HC3N ': 1e16 * u.cm**-2, ' C(18)O ': 7e17 * u.cm**-2, ' 13CS ': 4e15 * u.cm**-2,
    ' SO2 ': 9e16 * u.cm**-2, ' NH2CHO ': 3e16 * u.cm**-2, ' HNCO ': 6e16 * u.cm**-2,
    ' SO ': 2e16 * u.cm**-2, ' SiO ': 1e15 * u.cm**-2, ' H2S ': 7e16 * u.cm**-2,
    ' c-HCCCH ': 3e16 * u.cm**-2, 'HC3N v7=1': 3e15 * u.cm**-2, ' H213CO ': 1e16 * u.cm**-2,
    ' 13CH3CN ': 9e15 * u.cm**-2, ' CH2CHCN ': 2.5e16 * u.cm**-2, ' 18OCS ': 4e16 * u.cm**-2,
    ' CH3NCO, vb = 0 ': 2.1e17 * u.cm**-2, ' CH3CH2CN ': 7e16 * u.cm**-2, ' CH3COCH3 ': 6e17 * u.cm**-2,
}

ds7columns = {
    ' CH3OH ': 5.3e17 * u.cm**-2, ' CH3OCHO ': 2e17 * u.cm**-2, ' CH3CHO ': 4e16 * u.cm**-2,
    ' C2H5OH ': 2e17 * u.cm**-2, ' CH3OCH3 ': 2e17 * u.cm**-2, ' DCN ': 6e14 * u.cm**-2,
    ' OCS ': 2e16 * u.cm**-2, ' 13CH3OH ': 1e17 * u.cm**-2, ' H2CO ': 5e16 * u.cm**-2,
    ' HC3N ': 1e15 * u.cm**-2, ' C(18)O ': 8e17 * u.cm**-2, ' 13CS ': 2e15 * u.cm**-2,
    ' SO2 ': 4e15 * u.cm**-2, ' NH2CHO ': 5e15 * u.cm**-2, ' HNCO ': 1.5e16 * u.cm**-2,
    ' SO ': 5e15 * u.cm**-2, ' SiO ': 5e14 * u.cm**-2, ' H2S ': 4e16 * u.cm**-2,
    ' c-HCCCH ': 5e15 * u.cm**-2, 'HC3N v7=1': 2e15 * u.cm**-2, ' H213CO ': 6e15 * u.cm**-2,
    ' 13CH3CN ': 2e15 * u.cm**-2, ' CH2CHCN ': 7e15 * u.cm**-2, ' 18OCS ': 8e15 * u.cm**-2,
    ' CH3NCO, vb = 0 ': 4e16 * u.cm**-2, ' CH3CH2CN ': 8e15 * u.cm**-2, ' CH3COCH3 ': 1.5e17 * u.cm**-2,
}

ds8columns = {
    ' CH3OH ': 5e17 * u.cm**-2, ' CH3OCHO ': 2e17 * u.cm**-2, ' CH3CHO ': 3e16 * u.cm**-2,
    ' C2H5OH ': 2e17 * u.cm**-2, ' CH3OCH3 ': 2e17 * u.cm**-2, ' DCN ': 7e14 * u.cm**-2,
    ' OCS ': 3.5e16 * u.cm**-2, ' 13CH3OH ': 8e16 * u.cm**-2, ' H2CO ': 1e17 * u.cm**-2,
    ' HC3N ': 9e14 * u.cm**-2, ' C(18)O ': 6e17 * u.cm**-2, ' 13CS ': 1e15 * u.cm**-2,
    ' SO2 ': 3e15 * u.cm**-2, ' NH2CHO ': 7e15 * u.cm**-2, ' HNCO ': 1e16 * u.cm**-2,
    ' SO ': 5e15 * u.cm**-2, ' SiO ': 2e14 * u.cm**-2, ' H2S ': 5e16 * u.cm**-2,
    ' c-HCCCH ': 4e15 * u.cm**-2, 'HC3N v7=1': 1e15 * u.cm**-2, ' H213CO ': 9e15 * u.cm**-2,
    ' 13CH3CN ': 2e15 * u.cm**-2, ' CH2CHCN ': 6e15 * u.cm**-2, ' 18OCS ': 7e15 * u.cm**-2,
    ' CH3NCO, vb = 0 ': 4e16 * u.cm**-2, ' CH3CH2CN ': 7e15 * u.cm**-2, ' CH3COCH3 ': 1.5e17 * u.cm**-2,
}

ds9columns = {
    ' CH3OH ': 2.5e17 * u.cm**-2, ' CH3OCHO ': 7e16 * u.cm**-2, ' CH3CHO ': 1e16 * u.cm**-2,
    ' C2H5OH ': 5e16 * u.cm**-2, ' CH3OCH3 ': 7e16 * u.cm**-2, ' DCN ': 3e14 * u.cm**-2,
    ' OCS ': 2e16 * u.cm**-2, ' 13CH3OH ': 6e16 * u.cm**-2, ' H2CO ': 9e15 * u.cm**-2,
    ' HC3N ': 9e14 * u.cm**-2, ' C(18)O ': 2e17 * u.cm**-2, ' 13CS ': 6e14 * u.cm**-2,
    ' SO2 ': 2.5e15 * u.cm**-2, ' NH2CHO ': 2e15 * u.cm**-2, ' HNCO ': 5e15 * u.cm**-2,
    ' SO ': 7e15 * u.cm**-2, ' SiO ': 9e13 * u.cm**-2, ' H2S ': 4e16 * u.cm**-2,
    ' c-HCCCH ': 1e15 * u.cm**-2, 'HC3N v7=1': 7e14 * u.cm**-2, ' H213CO ': 2.5e15 * u.cm**-2,
    ' 13CH3CN ': 7e14 * u.cm**-2, ' CH3CH2CN ': 4e15 * u.cm**-2, ' CH2CHCN ': 2.3e15 * u.cm**-2,
}

ds10columns = {
    ' CH3OH ': 6e17 * u.cm**-2, ' CH3OCHO ': 3e15 * u.cm**-2, ' CH3CHO ': 9e14 * u.cm**-2,
    ' C2H5OH ': 9e15 * u.cm**-2, ' CH3OCH3 ': 2e15 * u.cm**-2, ' DCN ': 3e15 * u.cm**-2,
    ' OCS ': 1.5e17 * u.cm**-2, ' 13CH3OH ': 1e17 * u.cm**-2, ' H2CO ': 6e16 * u.cm**-2,
    ' HC3N ': 5e15 * u.cm**-2, ' C(18)O ': 3e19 * u.cm**-2, ' 13CS ': 1e16 * u.cm**-2,
    ' SO2 ': 6e15 * u.cm**-2, ' NH2CHO ': 5e14 * u.cm**-2, ' HNCO ': 5e15 * u.cm**-2,
    ' SO ': 5e16 * u.cm**-2, ' SiO ': 4e15 * u.cm**-2, ' H2S ': 4e17 * u.cm**-2,
    ' c-HCCCH ': 9e14 * u.cm**-2, 'HC3N v7=1': 5e15 * u.cm**-2, ' H213CO ': 1e16 * u.cm**-2,
    ' 13CH3CN ': 1.5e15 * u.cm**-2,
}

sourcecolumns={'SgrB2S':sgrb2scolumns,'DSi':dsicolumns, 'DSii':ds2columns,'DSiii':ds3columns,'DSiv':ds4columns,
               'DSv':ds5columns,'DSVI':ds6columns,'DSVII':ds7columns,'DSVIII':ds8columns,'DSIX':ds9columns,'DS10':ds10columns}