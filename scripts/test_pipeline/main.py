import logging
from pathlib import Path
import sys
import importlib
sys.path.append('../')  # Adjust the path to include the parent directory

utilities = importlib.import_module('utilities')
from utilities import *

def main():
    # Set up logging
    logging.basicConfig(level=logging.INFO)
    current_dir = Path(__file__).parent
    logging.info(f"Main directory: {current_dir}")
    # Set necessary variables for analysis, will probably outsource these to another script as the structure develops
    source='DSi'
    fieldnumber=fields[source]
    molecule='C2H5OH'
    results_dir = Path(f"results/{source}/")

    # Define paths for data and results
    model_line_dir = Path(f'../../linemodels/firstrelease/{source}/')
    representative_line_output_dir=results_dir/"representativelinetests/"
    logging.info(f'Creating directory for results of representative line measurement: {representative_line_output_dir}')
    results_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Identify representative lines
    num_best_lines=5
    logging.info(f"Selecting {molecule} representative line from {num_best_lines} candidate best lines in {model_line_dir}")
    import findrepresentativeline
    findrepresentativeline.run(model_line_dir, representative_line_output_dir, source, molecule, num_best_lines)

    logging.info("Pipeline complete.")

if __name__ == "__main__":
    main()