import logging
from pathlib import Path
import sys
import importlib
import pdb
sys.path.append('../')  # Adjust the path to include the parent directory

utilities = importlib.import_module('utilities')
from utilities import *

def main():
    # Set up logging
    logging.basicConfig(level=logging.INFO)
    logger= logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    current_dir = Path(__file__).parent
    logging.info(f"Main directory: {current_dir}")
    # Set necessary variables for analysis, will probably outsource these to another script as the structure develops
    source='DSi'
    fieldnumber=fields[source]
    molecule='CH3CH2CN'
    molecule_with_spaces=f' {molecule} ' #needed for some dictionaries
    results_dir = Path(f"results/{source}/{molecule}")
    logging.info(f'Creating home directory for all results: {results_dir}')
    results_dir.mkdir(parents=True, exist_ok=True)

    ### Step 1: Identify representative lines
    
    #Set number of candidate lines for fitting
    num_best_lines=5
    
    import findrepresentativeline
    print(type(results_dir))
    source_representativeline=findrepresentativeline.run(logger, results_dir, source, molecule, molecule_with_spaces,
     num_best_lines)
    print(source_representativeline)

    logging.info("Pipeline complete.")

if __name__ == "__main__":
    main()