import os
import numpy as np
import yaml
import matplotlib.pyplot as plt
import warnings
import pandas as pd

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
    try:
        if os.path.exists(errfile):
            with open(errfile, 'r') as fin:
                for line in fin:
                    if 'calphy.errors' in line:
                        break
                    error_code = line.split(':')[0].split('.')[-1]
    except:
        pass
    return error_code
    
def gather_results(mainfolder, reduce_composition=True, 
    extract_phase_prefix=False):
    """
    Gather results from all subfolders in a given folder into a Pandas DataFrame

    Parameters
    ----------
    mainfolder: string
        folder where calculations are stored
    
    reduce_composition: bool
        If True, per species composition arrays are added.
        Might be redundant.
    
    extract_phase_prefix: bool
        Should be used in conjuction with phase diagram mode. 
        Extracts the prefix and add it as a phase_name column.
    
    Returns
    -------
    df: pandas DataFrame
        DataFrame with results
    """
    try:
        import pandas as pd
    except ImportError:
        raise ImportError('Please install pandas to use this function')

    unique_elements = []
    datadict = {}
    datadict['calculation_mode'] = []
    datadict['status'] = []
    datadict['temperature'] = []
    datadict['pressure'] = []
    datadict['free_energy'] = []
    datadict['reference_phase'] = []
    datadict['error_code'] = []
    datadict['composition'] = []
    datadict['calculation'] = []
    datadict['ideal_entropy'] = []
    datadict['phase_name'] = []
    datadict['reference_composition'] = []
    
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
        datadict['calculation_mode'].append(mode)
        datadict['temperature'].append(inp['temperature'])
        datadict['pressure'].append(inp['pressure'])
        datadict['reference_phase'].append(inp['reference_phase'])
        datadict['phase_name'].append(inp['phase_name'])
        datadict['reference_composition'].append(inp['reference_composition'])
        datadict['composition'].append(None)
        datadict['ideal_entropy'].append(0)
        datadict['calculation'].append(folder)
    
        #check output file
        outfile = os.path.join(mainfolder, folder, 'report.yaml')
        datadict['error_code'].append(None)
        
        #print(inpfile)
        if not os.path.exists(outfile):
            datadict['status'].append('False')
            datadict['free_energy'].append(np.nan)
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
            el_arr = list(compdict.keys())

            #we also need to update entropy
            if 'entropy_contribution' in out['results'].keys():
                datadict['ideal_entropy'][-1] = -1*out['results']['entropy_contribution']

        for el in el_arr:
            if el not in unique_elements:
                unique_elements.append(el)

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

    if reduce_composition:
        unique_element_dict = {x: [] for x in unique_elements}
        for x in datadict['composition']:
            if x is None:
                for key in unique_element_dict.keys():
                    unique_element_dict[key].append(0)
            else:
                for el in unique_elements:
                    if el in x.keys():
                        unique_element_dict[el].append(x[el])
                    else:
                        unique_element_dict[el].append(0)
        #add the keys to datadict
        for key, val in unique_element_dict.items():
            datadict[key] = val
    
    if not extract_phase_prefix:
        del datadict['phase_name']

    df = pd.DataFrame(data=datadict)
    return df

def _entropy(compC, comp_init=[1, 0]):
    """
    l1: initial concentration [Au, Cu]
    l2: final concentration [Au, Cu]
    """
    compA = 1 - compC
    def _log(val):
        if val == 0:
            return 0
        else:
            return np.log(val)
    dA = compA*_log(compA) - comp_init[0]*_log(comp_init[0])
    dC = compC*_log(compC) - comp_init[1]*_log(comp_init[1])
    return kb*(dA+dC)
    
def clean_df(df, reference_element, combine_direct_calculations=False, smooth=False):
    """
    Clean a parsed dataframe and drop unnecessary columns. This gets it ready for further processing
    Note that `gather_results` should be run with `reduce_composition` and `extract_phase_name` for this to work.


    Parameters
    ----------
    df: DataFrame
        dataframe parsed by `gather_results` with `reduce_composition=True`.
    
    reference_element: str
        reference element from the compositions, which will be renamed to `composition`
    
    combine_direct_calculations: bool, optional
        If True, combine direct calculations by fitting to produce temperature and free energy arrays
        If used, an extra column `error` with RMSE of the fitting is also created
    
    smooth : bool, optional
        If True, smooth the F(T) data using the thermodynamic basis
        ``[1, T, T ln T, T²]``.  If False (default), return the raw data
        points without smoothing.

    Returns
    -------
    df: DataFrame
        combined, finished DataFrame
    """

    if "phase_name" not in df.keys():
        raise ValueError("phase_name key is not found, maybe add it?")

    df = df.loc[df.status=='True']
    df = df.drop(labels=['status', 'pressure', 'reference_phase', 
                 'error_code', 'composition', 'calculation'], axis='columns')

    phases = df.groupby(df.phase_name)
    phases = [phases.get_group(x) for x in phases.groups]

    df_dict = {}
    
    for phase in phases:
        if combine_direct_calculations:
            gb = phase.groupby(by=reference_element)
            gbs = [gb.get_group(x) for x in gb.groups]

            fes = []
            tes = []
            errors = []
            comps = []
            mode_list = []
            entropies = []
            is_refs = []

            for exdf in gbs:
                temps = np.array(exdf.temperature.values)
                fe = np.array(exdf.free_energy.values)
                modes = np.array(exdf.calculation_mode.values)
                entropy = np.array(exdf.ideal_entropy.values)
                comp_ref = np.array(exdf.reference_composition.values)

                unique_modes = np.unique(modes)
                if len(unique_modes)>1:
                    warnings.warn("mixing calculations from more than one mode!")
                unique_mode = unique_modes[0]

                #REMEMBER TO SORT EVERYTHING
                args = np.argsort(temps)
                temps = temps[args]
                fe = fe[args] 
                #entropy = entropy[args]        
                #print(fe, entropy)
                #print(len(fe), len(entropy))    
                
                if smooth:
                    # Thermodynamic basis: F(T) = a + b*T + c*T*ln(T) + d*T²
                    T_ = temps.astype(float)
                    basis = np.column_stack([np.ones_like(T_), T_,
                                            T_ * np.log(T_), T_**2])
                    coeffs_t, _, _, _ = np.linalg.lstsq(basis, fe,
                                                        rcond=None)

                    fe_eval = basis @ coeffs_t
                    error = float(np.sqrt(np.mean((fe_eval - fe)**2)))

                    temp_arr = np.arange(temps.min(), temps.max()+1,
                                         1).astype(float)
                    basis_arr = np.column_stack([
                        np.ones_like(temp_arr), temp_arr,
                        temp_arr * np.log(temp_arr), temp_arr**2])
                    fe_arr = basis_arr @ coeffs_t
                else:
                    fe_arr = fe
                    temp_arr = temps
                    error = 0

                fes.append(fe_arr)
                tes.append(temp_arr)
                errors.append(error)
                entropies.append(entropy[0])
                comps.append(float(exdf[reference_element].values[0]))
                # The fe job IS the reference by definition; composition_scaling jobs are not.
                is_refs.append(unique_mode == 'fe')

                mode_list.append(unique_mode)
            
            #replace df
            df = pd.DataFrame(data={'temperature':tes, 'free_energy': fes, 
                'error':errors, reference_element:comps, 'ideal_entropy': entropies,
                'calculation_mode': mode_list, "is_reference":is_refs})
        
        df = df.rename(columns={reference_element:'composition'})
        df_dict[phase.phase_name.values[0]] = df
    return df_dict

def fix_composition_scaling(dfdict, correct_entropy=True, add_ideal_entropy=False):
    """
    Correct composition-scaling free energies by adding the reference
    free energy and (optionally) subtracting the ideal-entropy term.

    Parameters
    ----------
    dfdict : dict of DataFrame
        Output from ``clean_df``.
    correct_entropy : bool
        If True **and** *add_ideal_entropy* is True, subtract
        ``T * S_ideal`` from the free energy.
    add_ideal_entropy : bool
        Controls whether the ideal-entropy correction is applied.
    """
    for key, val in dfdict.items():
        x = val
        ref_row = x.loc[x.is_reference == True]
        ref_fe   = np.asarray(ref_row.free_energy.values[0], dtype=float)
        ref_temp = np.asarray(ref_row.temperature.values[0], dtype=float)

        # Fit reference F(T) with thermodynamic basis
        basis_ref = np.column_stack([np.ones_like(ref_temp), ref_temp,
                                     ref_temp * np.log(ref_temp),
                                     ref_temp**2])
        ref_coeffs, _, _, _ = np.linalg.lstsq(basis_ref, ref_fe, rcond=None)

        for index, row in x.iterrows():
            if (not row.is_reference) and (row.calculation_mode == 'composition_scaling'):
                T_row = np.asarray(row.temperature, dtype=float)
                basis_row = np.column_stack([np.ones_like(T_row), T_row,
                                             T_row * np.log(T_row),
                                             T_row**2])
                ref_at_T = basis_row @ ref_coeffs

                if correct_entropy and add_ideal_entropy:
                    x.at[index, 'free_energy'] = (
                        row.free_energy
                        - row.temperature * row.ideal_entropy
                        + ref_at_T)
                else:
                    x.at[index, 'free_energy'] = row.free_energy + ref_at_T

        dfdict[key] = x
    return dfdict

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