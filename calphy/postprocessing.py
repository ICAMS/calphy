import os
import numpy as np
import yaml
import matplotlib.pyplot as plt
import warnings

def read_report(folder):
    """
    Read the finished calculation report

    Parameters
    ----------
    folder: string
        folder from which calculation is to be read
    
    Returns
    -------
    data: dict
        dictionary with results

    """
    repfile = os.path.join(folder, "report.yaml")
    if not os.path.exists(repfile):
        raise FileNotFoundError(f"file {repfile} not found")

    with open(repfile, 'r') as fin:
        data = yaml.safe_load(fin)
    return data

def _extract_error(errfile):
    error_code = None
    if os.path.exists(errfile):
        with open(errfile, 'r') as fin:
            for line in fin:
                if 'calphy.errors' in line:
                    break
            try:
                error_code = line.split(':')[0].split('.')[-1]
            except:
                pass
    return error_code
    
def gather_results(mainfolder):
    """
    Gather results from all subfolders in a given folder into a Pandas DataFrame

    Parameters
    ----------
    mainfolder: string
        folder where calculations are stored
    
    Returns
    -------
    df: pandas DataFrame
        DataFrame with results
    """
    try:
        import pandas as pd
    except ImportError:
        raise ImportError('Please install pandas to use this function')
     
    datadict = {}
    datadict['mode'] = []
    datadict['status'] = []
    datadict['temperature'] = []
    datadict['pressure'] = []
    datadict['free_energy'] = []
    datadict['reference_phase'] = []
    datadict['error_code'] = []
    datadict['composition'] = []
    datadict['calculation'] = []
    
    folders = next(os.walk(mainfolder))[1]
    for folder in folders:
        #adjust for pyiron folder, see
        if folder.split('_')[-1] == 'hdf5':
            #this could be a pyiron calc
            withouthdf = folder.split('_hdf5')[0]
            folder = f'{folder}/{withouthdf}'
            
        inpfile = os.path.join(mainfolder, folder, 'input_file.yaml')
        #print(inpfile)
        if not os.path.exists(inpfile):
            continue;
        
        #ok, valid calculation, try to parse input file to get info
        with open(inpfile, 'r') as fin:
            inp = yaml.safe_load(fin)
        #grab the first calculation
        inp = inp['calculations'][0]
        #mode
        mode = inp['mode']
        datadict['mode'].append(mode)
        datadict['temperature'].append(inp['temperature'])
        datadict['pressure'].append(inp['pressure'])
        datadict['reference_phase'].append(inp['reference_phase'])
        datadict['composition'].append(None)
        datadict['calculation'].append(folder)
    
        #check output file
        outfile = os.path.join(mainfolder, folder, 'report.yaml')
        datadict['error_code'].append(None)
        
        #print(inpfile)
        if not os.path.exists(outfile):
            datadict['status'].append('False')
            datadict['free_energy'].append(np.NaN)
            #check if error file is found
            errfile = os.path.join(os.getcwd(), mainfolder, folder+'.sub.err')
            datadict['error_code'][-1] = _extract_error(errfile)
            continue;
    
        if mode in ['fe', 'alchemy', 'composition_scaling']:
            datadict['status'].append('True')
        
        #ok, valid calculation, try to parse input file to get info
        with open(outfile, 'r') as fin:
            out = yaml.safe_load(fin)
    
        datadict['free_energy'].append(out['results']['free_energy'])
        
        #add normal composition
        el_arr = np.array(out['input']['element'].split(' ')).astype(str)
        comp_arr = np.array(out['input']['concentration'].split(' ')).astype(float)
        composition = {x:y for x,y in zip(el_arr, comp_arr)}
        datadict['composition'][-1] = composition

        if mode == 'composition_scaling':
            #we need to update composition
            compdict = inp['composition_scaling']['output_chemical_composition']
            maxatoms = np.sum([val for key, val in compdict.items()])
            for key, val in compdict.items():
                compdict[key] = val/maxatoms
            datadict['composition'][-1] = compdict

        #parse extra info
        if mode in ['ts', 'tscale']:
            datafile = os.path.join(os.getcwd(), mainfolder, folder, 'temperature_sweep.dat')
            if os.path.exists(datafile):
                datadict['status'].append('True')
                t, f = np.loadtxt(datafile, unpack=True, usecols=(0,1))
                datadict['temperature'][-1] = t
                datadict['free_energy'][-1] = f
            else:
                datadict['status'].append('False')
                errfile = os.path.join(os.getcwd(), mainfolder, folder+'.sub.err')
                datadict['error_code'][-1] = _extract_error(errfile)

    df = pd.DataFrame(data=datadict)
    return df

def find_transition_temperature(folder1, folder2, fit_order=4, plot=True):
    """
    Find transition temperature where free energy of two phases are equal.

    Parameters
    ----------
    folder1: string
        directory with temperature scale calculation

    folder2: string
        directory with temperature scale calculation

    fit_order: int, optional
        default 4. Order for polynomial fit of temperature vs free energy
    
    plot: bool, optional
        default True. Plot the results.
    """
    file1 = os.path.join(folder1, 'temperature_sweep.dat')
    file2 = os.path.join(folder2, 'temperature_sweep.dat')
    if not os.path.exists(file1):
        raise FileNotFoundError(f'{file1} does not exist')
    if not os.path.exists(file2):
        raise FileNotFoundError(f'{file2} does not exist')

    t1, f1 = np.loadtxt(file1, unpack=True, usecols=(0,1))
    t2, f2 = np.loadtxt(file2, unpack=True, usecols=(0,1))

    #do some fitting to determine temps
    t1min = np.min(t1)
    t2min = np.min(t2)
    t1max = np.max(t1)
    t2max = np.max(t2)

    tmin = np.min([t1min, t2min])
    tmax = np.max([t1max, t2max])

    #warn about extrapolation
    if not t1min == t2min:
        warnings.warn(f'free energy is being extrapolated!')
    if not t1max == t2max:
        warnings.warn(f'free energy is being extrapolated!')

    #now fit
    f1fit = np.polyfit(t1, f1, fit_order)
    f2fit = np.polyfit(t2, f2, fit_order)

    #reevaluate over the new range
    fit_t = np.arange(tmin, tmax+1, 1)
    fit_f1 = np.polyval(f1fit, fit_t)
    fit_f2 = np.polyval(f2fit, fit_t)

    #now evaluate the intersection temp
    arg = np.argsort(np.abs(fit_f1-fit_f2))[0]
    transition_temp = fit_t[arg]

    #warn if the temperature is shady
    if np.abs(transition_temp-tmin) < 1E-3:
        warnings.warn('It is likely there is no intersection of free energies')
    elif np.abs(transition_temp-tmax) < 1E-3:
        warnings.warn('It is likely there is no intersection of free energies')

    #plot
    if plot:
        c1lo = '#ef9a9a'
        c1hi = '#b71c1c'
        c2lo = '#90caf9'
        c2hi = '#0d47a1'

        plt.plot(fit_t, fit_f1, color=c1lo, label=f'{folder1} fit')
        plt.plot(fit_t, fit_f2, color=c2lo, label=f'{folder2} fit')
        plt.plot(t1, f1, color=c1hi, label=folder1, ls='dashed')
        plt.plot(t2, f2, color=c2hi, label=folder2, ls='dashed')
        plt.axvline(transition_temp, ls='dashed', c='#37474f')
        plt.ylabel('Free energy (eV/atom)')
        plt.xlabel('Temperature (K)')
        plt.legend(frameon=False)
    return transition_temp