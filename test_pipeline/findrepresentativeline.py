# %%
from utilities import *
import matplotlib.pyplot as plt
from astropy.table import QTable
import glob
import sys
from astropy.modeling.fitting import LMLSQFitter
from astropy.modeling.functional_models import Gaussian1D
from collections import Counter
plt.rcParams['figure.dpi'] = 300

# %%
# This script uses 1D gaussian fitting to identify the representative line for a given molecule.
# It is assumed that the molecule has been previously identified in the pipeline, but this may be changed to import from another file in the future

def run(logger, results_dir, source='DSi', molecule='C2H5OH', molecule_with_spaces=' C2H5OH ', num_best_lines=5):
    # Set logger, input model line directory, and output directory
    model_line_dir = Path(f'../../linemodels/firstrelease/{source}/')
    representative_line_output_dir=results_dir/'representativelinetests/'
    logger.info(f'Creating directory for results of representative line measurement: {representative_line_output_dir}')

    representative_line_plots_home=representative_line_output_dir/'plots/'
    representative_line_plots_home.mkdir(parents=True, exist_ok=True)  # Create the directory if it doesn't exist
    
    #Open the table of minimally contaminated, minimally self-contaminated lines
    safelinetable=QTable.read(model_line_dir / f'{molecule}.fits')

    #Sort the table by predicted line brightness
    safelinetable.sort('ModelBrightness', reverse=True)
    targetlines= safelinetable[:(num_best_lines)]
    print(targetlines)

    # %%
    representative_spectra_paths=glob.glob(f'{representative_spectra_home}{source}/*.fits')
    representative_spectra_paths.sort()
    print(representative_spectra_paths)

    # %%
    # Loop through data cubes and identify lines that are present within
    # Then fit a 1D gaussian to lines that are present
    model_fit_rows=[]
    for representative_spectra_path in representative_spectra_paths:
        # Read in the respresentative spectra for the source
        representative_spectra = QTable.read(representative_spectra_path)
        targetlines_in_cube = targetlines[(targetlines['ReferenceFrequency'] > min(representative_spectra['Frequency'])) & 
                                        (targetlines['ReferenceFrequency'] < max(representative_spectra['Frequency']))]
        if len(targetlines_in_cube) == 0: #Skip if no lines are present in the cube
            continue
        else:
            # Fit the lines with a 1D gaussian
            for targetline in targetlines_in_cube:
                print(f'Found {targetline['QNs']} between {min(representative_spectra['Frequency']).to("GHz")} and {max(representative_spectra['Frequency']).to("GHz")} in {representative_spectra_path}')
                fitter = LMLSQFitter()  
                # Set the plotting range to be 10 km/s around the line
                # This is arbitrary, but should be sufficient in most cases
                plotting_width=velocitytofreq(10*u.km/u.s, targetline['ReferenceFrequency'])
                min_plot_frequency = targetline['ReferenceFrequency'] - plotting_width
                max_plot_frequency = targetline['ReferenceFrequency'] + plotting_width
                frequency_to_plot= representative_spectra['Frequency'][(representative_spectra['Frequency'] > min_plot_frequency) & 
                                                                        (representative_spectra['Frequency'] < max_plot_frequency)]
                plotting_indices=np.where([representative_spectra['Frequency'] == x for x in frequency_to_plot])[1]
                flux_for_plot = representative_spectra['BrightnessTemperature'][plotting_indices].value

                # Make sure to calculate the linewidth needed for fitting in frequency space
                model_line_width= velocitytofreq(4*u.km/u.s, targetline['ReferenceFrequency'])
                model_line_stddev = model_line_width.to('GHz').value / (2 * np.sqrt(2 * np.log(2)))
                min_fit_frequency= targetline['ReferenceFrequency'] - model_line_width
                max_fit_frequency=targetline['ReferenceFrequency'] + model_line_width
                # Set the flux and frequency range for the fit
                frequency_to_fit = representative_spectra['Frequency'][(representative_spectra['Frequency']>min_fit_frequency) & (representative_spectra['Frequency']<max_fit_frequency)]
                target_indices=np.where([representative_spectra['Frequency'] == x for x in frequency_to_fit])[1]
                flux_to_fit = representative_spectra['BrightnessTemperature'][target_indices].value
                # Fit the gaussian model
                gaussian_model = Gaussian1D(amplitude=targetline['ModelBrightness'].value, mean=targetline['ReferenceFrequency'].value, stddev=model_line_stddev)
                fitted_model = fitter(gaussian_model, frequency_to_fit.to('GHz').value, flux_to_fit)

                cov=fitter.fit_info['param_cov']  # This will give the covariance matrix of the fit parameters
                errors_on_fit= np.sqrt(np.diag(cov)) if cov is not None else None #Convert covariance matrix to 1 sigma errors
                if errors_on_fit is None:
                    print("Covariance matrix not available. Cannot calculate errors on fit parameters.")
                    continue
                else:
                    print(f'Fitted parameters: {fitted_model.parameters}')
                    print(f'Errors: {errors_on_fit}')

                    string_manip_function=stringmanipulationdict[molecule_with_spaces]
                    qns_for_figpath=string_manip_function(targetline['QNs'])
                    figpath=representative_line_plots_home / f'{qns_for_figpath}_fitted.png'
                    #Append the results to the error table
                    model_fit_rows.append([str(targetline['QNs']),fitted_model.amplitude.value*u.K,errors_on_fit[0]*u.K, 
                                        fitted_model.mean.value*u.GHz, errors_on_fit[1]*u.GHz, 
                                        fitted_model.stddev.value*u.GHz, errors_on_fit[2]*u.GHz])
                    #Plot output
                    testplot=plt.plot(frequency_to_plot.to('GHz'), flux_for_plot,drawstyle='steps-mid', color='black',label='Data')
                    testgaussian=plt.plot(frequency_to_plot.to('GHz'), fitted_model(frequency_to_plot.to('GHz').value), color='red',label='Fit')
                    plt.xlabel(r'$\nu$ [GHz]')
                    plt.ylabel(r'$T_b$ [K]')
                    plt.legend()
                    plt.savefig(figpath)
                    plt.show()
                    #sys.exit()
    #spdb.set_trace()
    model_fit_table=QTable([np.array(model_fit_rows).T],names=['QNs','Brightness','Error_Brightness', 'Frequency', 'Error_Frequency', 'Sigma', 'Error_Sigma'],
           dtype=['str', 'float', 'float', 'float', 'float', 'float', 'float'],
           units=['', 'K', 'K', 'GHz', 'GHz', 'GHz', 'GHz'])
    
    print(model_fit_table)

    # %%
    #Find the transition that has the lowest error in the most columns
    errors_for_best_line_selection=np.vstack([model_fit_table['Error_Brightness'].value,model_fit_table['Error_Frequency'].value,model_fit_table['Error_Sigma'].value]).T  # Stack the errors in a 2D array
    # This will give you a 2D array where each row corresponds to a line and each column corresponds to an error type (Brightness, Frequency, Sigma)
    min_indices= np.argmin(errors_for_best_line_selection, axis=0)  # Get the indices of the minimum values in each column
    min_counts=Counter(min_indices) #Note, this is a counter. Supplying it with indices will count how many times each index appears in the min_indices array

    best_index=min_counts.most_common(1)[0][0]  # Get the index of the most common minimum index. most_common(1) returns a list of tuples, where the first element is the index and the second element is the count
    #Note, most_common(0) returns an empty list because it specifies how many most common elements to return. If you want to return all elements, you can use most_common() without any arguments.
    best_line=model_fit_table[best_index]
    # Grab the best line from the input table of lines
    input_table_best_line=targetlines[best_index]

    print(input_table_best_line)
    return input_table_best_line



