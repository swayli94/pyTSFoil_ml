'''

'''
import os
import sys
import json

path = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(path, '..', '..'))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

import numpy as np
import time
from typing import Dict, Any
import multiprocessing as mp

from cst_modeling.foil import cst_foil 
from pytsfoil_ibl import PyTSFoil


N_PROCESSES = 40
N_TEST_AIRFOILS = 40


airfoil_database_path = os.path.join(project_root, 'data', 'airfoil_database.json')
os.makedirs(os.path.join(path, 'results'), exist_ok=True)


def load_airfoil_database() -> Dict[str, Any]:
    '''
    Load the airfoil database from the JSON file.
    
    Covert data to int, float, and list.
    
    Returns:
    --------
    Dict[str, Any]: Airfoil database
    '''
    t0 = time.time()
    with open(airfoil_database_path, 'r') as f:
        airfoil_database = json.load(f)
        
    t1 = time.time()
    print(f"Loaded {len(airfoil_database)} airfoils in {t1-t0:.2f} seconds")
    
    return airfoil_database


def run_pytsfoil(params: Dict[str, Any]) -> Dict[str, Any]:
    """
    Worker function to run a single PyTSFoil analysis in a separate process.
    
    Parameters:
    -----------
    params: Dict[str, Any]
        Dictionary containing:
        - 'case_name': Name identifier for this case
        - 'cst_u', 'cst_l': CST parameters
        - 'config': Configuration parameters for PyTSFoil
        
    Returns:
    --------
    Dict[str, Any]: Results dictionary with case info and outcomes
    """
    try:
        start_time = time.time()
        
        xx, yu, yl, _, _ =cst_foil(nn=201, cst_u=params['cst_u'], cst_l=params['cst_l'])
        x = np.concatenate((xx[::-1], xx[1:]))
        y = np.concatenate((yu[::-1], yl[1:]))
        airfoil_coordinates = np.column_stack((x, y))
        
        # Create PyTSFoil instance (each process gets its own isolated Fortran data)
        pytsfoil = PyTSFoil(airfoil_coordinates=airfoil_coordinates)
        
        # Set configuration
        pytsfoil.set_config(**params['config'])
        
        pytsfoil.set_config(
            flag_output=False,
            flag_output_summary=False,
            flag_output_shock=False,
            flag_output_field=False,
            flag_print_info=False
            )
        
        # Run analysis
        pytsfoil.run()
        
        elapsed_time = time.time() - start_time
        
        # Extract Mach distribution data for plotting
        results = {
            'xx': pytsfoil.mesh['xx'].copy(),
            'cpu': pytsfoil.data_summary['cpu'].copy(),
            'cpl': pytsfoil.data_summary['cpl'].copy(),
            'mau': pytsfoil.data_summary['mau'].copy(),
            'mal': pytsfoil.data_summary['mal'].copy(),
            'cl': pytsfoil.data_summary['cl'],
            'cd': pytsfoil.data_summary['cd'],
            'cm': pytsfoil.data_summary['cm']
        }
        
        # Save results to JSON file
        results_save = params.copy()
        results_save.update({
            'elapsed_time': elapsed_time,
            'cl': pytsfoil.data_summary['cl'],
            'cd': pytsfoil.data_summary['cd'],
            'cm': pytsfoil.data_summary['cm'],
            'xx': pytsfoil.mesh['xx'].copy().tolist(),
            'cpu': pytsfoil.data_summary['cpu'].copy().tolist(),
            'cpl': pytsfoil.data_summary['cpl'].copy().tolist(),
            'mau': pytsfoil.data_summary['mau'].copy().tolist(),
            'mal': pytsfoil.data_summary['mal'].copy().tolist(),
        })
        with open(os.path.join(path, 'results', f'{params['i_case']}.json'), 'w') as f:
            json.dump(results_save, f)
        
        # Get results
        return {
            'i_case': params['i_case'],
            'success': True,
            'elapsed_time': elapsed_time,
            'results': results,
            'message': f"Analysis completed successfully"
        }
        
    except Exception as e:
        print(f"Error in case {params['i_case']}: {str(e)}")
        return {
            'i_case': params['i_case'],
            'success': False,
            'elapsed_time': 0,
            'error': str(e)
        }


def load_results(file_name: str) -> Dict[str, Any]:
    '''
    Load the results from the JSON files.
    
    Returns:
    --------
    Dict[str, Any]: Results dictionary
    '''
    results = {}
    
    with open(file_name, 'r') as f:
        results = json.load(f)
        
        results['Aoa'] = float(results['config']['ALPHA'])
        results['Mach'] = float(results['config']['EMACH'])
        # results['Re'] = float(results['config']['REYNLD'])

    return results


if __name__ == "__main__":
    
    print('Working directory:', path)
    
    airfoil_database = load_airfoil_database()
    
    if N_TEST_AIRFOILS > 0:
        airfoil_database = airfoil_database[:N_TEST_AIRFOILS]
        print(f"Running {N_TEST_AIRFOILS} test airfoils")
    
    # Define multiple cases to run in parallel
    # Example: varying angle of attack and Mach number
    base_config = {
        'MAXIT': 9999,
        'NWDGE': 2,
        'n_point_x': 200,
        'n_point_y': 80,
        'n_point_airfoil': 100,
        'EPS': 0.2,
        'CVERGE': 1e-6,
    }
    
    # Create multiple parameter combinations
    cases = []
    i_case = 0
    for entry in airfoil_database:
        case_config = base_config.copy()
        case_config.update({
            'ALPHA': entry['AoA'],
            'EMACH': entry['Ma'],
            'REYNLD': entry['Re'],
        })
        cases.append({
            'i_case': i_case,
            'ID': entry['ID'],
            'cst_u': entry['cst_u'],
            'cst_l': entry['cst_l'],
            'config': case_config
        })
        i_case += 1
        
    print(f"Running {len(cases)} cases in parallel...")
    
    # Run analyses in parallel using multiprocessing
    start_time = time.time()
    
    # Use multiprocessing.Pool for parallel execution
    # Each process gets its own isolated PyTSFoil/Fortran data
    with mp.Pool(processes=N_PROCESSES) as pool:
        results = pool.map(run_pytsfoil, cases)
            
    total_time = time.time() - start_time
    
    # Print summary of results
    print(f"\n{'='*60}")
    print("MULTIPROCESSING ANALYSIS SUMMARY")
    print(f"{'='*60}")
    print(f"Total cases: {len(cases)}")
    print(f"Total time: {total_time:.2f}s")
    time_cases = [r['elapsed_time'] for r in results]
    print(f"Average time per case: {sum(time_cases)/len(cases):.2f}s")
    
    successful_cases = [r for r in results if r['success']]
    failed_cases = [r for r in results if not r['success']]
    
    print(f"Successful: {len(successful_cases)}")
    print(f"Failed: {len(failed_cases)}")
    



