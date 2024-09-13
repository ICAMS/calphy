import os
import numpy as np
import yaml

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