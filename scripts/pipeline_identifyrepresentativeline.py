from utilities import *
import matplotlib.pyplot as plt
from astropy.table import QTable

# This script uses 1D gaussian fitting to identify the representative line for a given molecule.
# It is assumed that the molecule has been previously identified in the pipeline, but this may be changed to import from another file in the future
molecule='CH3CH2CN'
safelinetable=QTable.read(f'../linemodels/firstrelease/DSi/{molecule}.fits')

print(safelinetable)
