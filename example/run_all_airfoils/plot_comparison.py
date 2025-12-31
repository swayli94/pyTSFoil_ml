'''
Random select N results, compare the Cp and Mw distribution between:

- data/airfoil_database.json
- results/*.json


'''
import os
import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, Any

path = os.path.dirname(os.path.abspath(__file__))

from run_airfoils import load_airfoil_database, load_results

results_dir = os.path.join(path, 'results')
figure_dir = os.path.join(path, 'figures')
os.makedirs(figure_dir, exist_ok=True)


N_SELECT = 10


def plot_comparison(i_case: int, airfoil_entry: Dict[str, Any], results: Dict[str, Any]):
    '''
    Plot the comparison between the airfoil database and the results.
    '''
    
    tsf_xx = np.array(results['xx'])
    tsf_mask = (tsf_xx >= 0) & (tsf_xx <= 1)
    tsf_xx = tsf_xx[tsf_mask]
    tsf_cpu = np.array(results['cpu'])[tsf_mask]
    tsf_cpl = np.array(results['cpl'])[tsf_mask]
    tsf_mau = np.array(results['mau'])[tsf_mask]
    tsf_mal = np.array(results['mal'])[tsf_mask]
    
    name = f'{airfoil_entry["ID"]}-{airfoil_entry["AoA"]:.2f}-{airfoil_entry["Ma"]:.2f}'

    fig, ax = plt.subplots(1, 2, figsize=(12, 6))
    
    ax[0].plot(airfoil_entry['x'], airfoil_entry['cpu'], 'b-', linewidth=2, label='RANS')
    ax[0].plot(airfoil_entry['x'], airfoil_entry['cpl'], 'b-', linewidth=2, label=None)
    ax[0].plot(tsf_xx, tsf_cpu, 'r-', linewidth=2, label='pyTSFoil')
    ax[0].plot(tsf_xx, tsf_cpl, 'r-', linewidth=2, label=None)
    ax[0].legend()
    ax[0].set_xlabel('X/c')
    ax[0].set_ylabel('Cp')
    ax[0].invert_yaxis()

    ax[1].plot(tsf_xx, tsf_mau, 'r-', linewidth=2, label='pyTSFoil')
    ax[1].plot(tsf_xx, tsf_mal, 'r-', linewidth=2, label=None)
    ax[1].legend()
    ax[1].set_xlabel('X/c')
    ax[1].set_ylabel('Mach Number')

    plt.suptitle(name)
    plt.tight_layout()

    plt.savefig(os.path.join(figure_dir, f'comparison-{i_case}.png'), dpi=100)
    plt.close()


if __name__ == "__main__":
    
    airfoil_database = load_airfoil_database()    
    
    indices = np.random.choice(len(airfoil_database), size=N_SELECT, replace=False)
    indices = [0, 5, 11, 19, 24]
    
    for i_case in indices:
        
        print(f'# {i_case}')
        
        airfoil_entry = airfoil_database[i_case]
        
        results_file = os.path.join(results_dir, f'{i_case}.json')
        results = load_results(results_file)
        
        plot_comparison(i_case, airfoil_entry, results)

