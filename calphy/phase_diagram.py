import numpy as np
from tqdm.notebook import trange
import pandas as pd
import matplotlib.pyplot as plt
import warnings
import itertools
from itertools import combinations
import math
import copy
import os
from calphy.composition_transformation import CompositionTransformation
import yaml
import matplotlib.patches as mpatches

from calphy.integrators import kb

from scipy.spatial import ConvexHull
from scipy.interpolate import splrep, splev
from scipy.optimize import curve_fit


colors = ['#a6cee3','#1f78b4','#b2df8a',
'#33a02c','#fb9a99','#e31a1c',
'#fdbf6f','#ff7f00','#cab2d6',
'#6a3d9a','#ffff99','#b15928']

matcolors = {
	"amber": {
		50 : '#fff8e1',
		100 : '#ffecb3',
		200 : '#ffe082',
		300 : '#ffd54f',
		400 : '#ffca28',
		500 : '#ffc107',
		600 : '#ffb300',
		700 : '#ffa000',
		800 : '#ff8f00',
		900 : '#ff6f00',
	},
	"blue_grey": {
		50 : '#ECEFF1',
		100 : '#CFD8DC',
		200 : '#B0BEC5',
		300 : '#90A4AE',
		400 : '#78909C',
		500 : '#607D8B',
		600 : '#546E7A',
		700 : '#455A64',
		800 : '#37474F',
		900 : '#263238',
	},
	"blue": {
		50 : '#E3F2FD',
		100 : '#BBDEFB',
		200 : '#90CAF9',
		300 : '#64B5F6',
		400 : '#42A5F5',
		500 : '#2196F3',
		600 : '#1E88E5',
		700 : '#1976D2',
		800 : '#1565C0',
		900 : '#0D47A1',
	},
	"brown": {
		50 : '#EFEBE9',
		100 : '#D7CCC8',
		200 : '#BCAAA4',
		300 : '#A1887F',
		400 : '#8D6E63',
		500 : '#795548',
		600 : '#6D4C41',
		700 : '#5D4037',
		800 : '#4E342E',
		900 : '#3E2723',
	},
	"cyan": {
		50 : '#E0F7FA',
		100 : '#B2EBF2',
		200 : '#80DEEA',
		300 : '#4DD0E1',
		400 : '#26C6DA',
		500 : '#00BCD4',
		600 : '#00ACC1',
		700 : '#0097A7',
		800 : '#00838F',
		900 : '#006064',
	},
	"deep_orange": {
		50 : '#FBE9E7',
		100 : '#FFCCBC',
		200 : '#FFAB91',
		300 : '#FF8A65',
		400 : '#FF7043',
		500 : '#FF5722',
		600 : '#F4511E',
		700 : '#E64A19',
		800 : '#D84315',
		900 : '#BF360C',
	},
	"deep_purple": {
		50 : '#EDE7F6',
		100 : '#D1C4E9',
		200 : '#B39DDB',
		300 : '#9575CD',
		400 : '#7E57C2',
		500 : '#673AB7',
		600 : '#5E35B1',
		700 : '#512DA8',
		800 : '#4527A0',
		900 : '#311B92',
	},
	"green": {
		50 : '#E8F5E9',
		100 : '#C8E6C9',
		200 : '#A5D6A7',
		300 : '#81C784',
		400 : '#66BB6A',
		500 : '#4CAF50',
		600 : '#43A047',
		700 : '#388E3C',
		800 : '#2E7D32',
		900 : '#1B5E20',
	},
	"grey": {
		50 : '#FAFAFA',
		100 : '#F5F5F5',
		200 : '#EEEEEE',
		300 : '#E0E0E0',
		400 : '#BDBDBD',
		500 : '#9E9E9E',
		600 : '#757575',
		700 : '#616161',
		800 : '#424242',
		900 : '#212121',
	},
	"indigo": {
		50 : '#E8EAF6',
		100 : '#C5CAE9',
		200 : '#9FA8DA',
		300 : '#7986CB',
		400 : '#5C6BC0',
		500 : '#3F51B5',
		600 : '#3949AB',
		700 : '#303F9F',
		800 : '#283593',
		900 : '#1A237E',
	},
	"light_blue": {
		50 : '#E1F5FE',
		100 : '#B3E5FC',
		200 : '#81D4FA',
		300 : '#4FC3F7',
		400 : '#29B6F6',
		500 : '#03A9F4',
		600 : '#039BE5',
		700 : '#0288D1',
		800 : '#0277BD',
		900 : '#01579B',
	},
	"light_green": {
		50 : '#F1F8E9',
		100 : '#DCEDC8',
		200 : '#C5E1A5',
		300 : '#AED581',
		400 : '#9CCC65',
		500 : '#8BC34A',
		600 : '#7CB342',
		700 : '#689F38',
		800 : '#558B2F',
		900 : '#33691E',
	},
	"lime": {
		50 : '#F9FBE7',
		100 : '#F0F4C3',
		200 : '#E6EE9C',
		300 : '#DCE775',
		400 : '#D4E157',
		500 : '#CDDC39',
		600 : '#C0CA33',
		700 : '#AFB42B',
		800 : '#9E9D24',
		900 : '#827717',
	},
	"orange": {
		50 : '#FFF3E0',
		100 : '#FFE0B2',
		200 : '#FFCC80',
		300 : '#FFB74D',
		400 : '#FFA726',
		500 : '#FF9800',
		600 : '#FB8C00',
		700 : '#F57C00',
		800 : '#EF6C00',
		900 : '#E65100',
	},
	"pink": {
		50 : '#FCE4EC',
		100 : '#F8BBD0',
		200 : '#F48FB1',
		300 : '#F06292',
		400 : '#EC407A',
		500 : '#E91E63',
		600 : '#D81B60',
		700 : '#C2185B',
		800 : '#AD1457',
		900 : '#880E4F',
	},
	"purple": {
		50 : '#F3E5F5',
		100 : '#E1BEE7',
		200 : '#CE93D8',
		300 : '#BA68C8',
		400 : '#AB47BC',
		500 : '#9C27B0',
		600 : '#8E24AA',
		700 : '#7B1FA2',
		800 : '#6A1B9A',
		900 : '#4A148C',
	},
	"red": {
		50 : '#FFEBEE',
		100 : '#FFCDD2',
		200 : '#EF9A9A',
		300 : '#E57373',
		500 : '#F44336',
		600 : '#E53935',
		700 : '#D32F2F',
		800 : '#C62828',
		900 : '#B71C1C',
	},
	"teal": {
		50 : '#E0F2F1',
		100 : '#B2DFDB',
		200 : '#80CBC4',
		300 : '#4DB6AC',
		400 : '#26A69A',
		500 : '#009688',
		600 : '#00897B',
		700 : '#00796B',
		800 : '#00695C',
		900 : '#004D40',
	},
	"yellow": {
		50 : '#FFFDE7',
		100 : '#FFF9C4',
		200 : '#FFF59D',
		300 : '#FFF176',
		400 : '#FFEE58',
		500 : '#FFEB3B',
		600 : '#FDD835',
		700 : '#FBC02D',
		800 : '#F9A825',
		900 : '#F57F17',
	}
}

def fix_data_file(datafile, nelements):
    """
    Change the atom types keyword in the structure file
    """
    lines = []
    with open(datafile, 'r') as fin:
        for line in fin:
            if 'atom types' in line:
                lines.append(f'{nelements} atom types\n')
            else:
                lines.append(line)
    outfile = datafile + 'mod.data'
    with open(outfile, 'w') as fout:
        for line in lines:
            fout.write(line)
    return outfile

class CScale:
    def __init__(self):
        self._input_chemical_composition = None
        self._output_chemical_composition = None
        self.restrictions = []

    @property
    def input_chemical_composition(self):
        return self._input_chemical_composition

    @property
    def output_chemical_composition(self):
        return self._output_chemical_composition


class SimpleCalculation:
    """
    Simple calc class 
    """
    def __init__(self, lattice,
                element,
                input_chemical_composition,
                output_chemical_composition):
        self.lattice = lattice
        self.element = element
        self.composition_scaling = CScale()
        self.composition_scaling._input_chemical_composition = input_chemical_composition
        self.composition_scaling._output_chemical_composition = output_chemical_composition


def prepare_inputs_for_phase_diagram(inputyamlfile, calculation_base_name=None):
    with open(inputyamlfile, 'r') as fin:
        data = yaml.safe_load(fin)

    if calculation_base_name is None:
        calculation_base_name = inputyamlfile    
    
    for phase in data['phases']:
        phase_reference_state = phase['reference_phase']
        phase_name = phase['phase_name']

        comps = phase['composition']
        reference_element = comps["reference_element"]
        use_composition_scaling = bool(comps["use_composition_scaling"])
        if str(phase_reference_state) == 'liquid':
            use_composition_scaling = False

        other_element_list = copy.deepcopy(phase['element'])
        other_element_list.remove(reference_element)
        other_element = other_element_list[0]

        #convert to list if scalar
        if not isinstance(comps['range'], list):
            comps["range"] = [comps["range"]]
        if len(comps["range"]) == 2: 
            comp_arr = np.arange(comps['range'][0], comps['range'][-1], comps['interval'])
            last_val = comps['range'][-1]
            if last_val not in comp_arr:
                comp_arr = np.append(comp_arr, last_val)
            ncomps = len(comp_arr)
            is_reference = np.abs(comp_arr-comps['reference']) < 1E-5
        elif len(comps["range"]) == 1:
            ncomps = 1
            comp_arr = [comps["range"][0]]
            is_reference = [True]
        else:
            raise ValueError("Composition range should be scalar of list of two values!")

        temps = phase["temperature"]
        if not isinstance(temps['range'], list):
            temps["range"] = [temps["range"]]
        if len(temps["range"]) == 2: 
            ntemps = int((temps['range'][-1]-temps['range'][0])/temps['interval'])+1
            temp_arr = np.linspace(temps['range'][0], temps['range'][-1], ntemps, endpoint=True)
        elif len(temps["range"]) == 1:
            ntemps = 1
            temp_arr = [temps["range"][0]]
        else:
            raise ValueError("Temperature range should be scalar of list of two values!")

        all_calculations = []

        for count, comp in enumerate(comp_arr):
            #check if ref comp equals given comp
            if is_reference[count]:
                #copy the dict
                calc = copy.deepcopy(phase)

                #pop extra keys which are not needed
                #we dont kick out phase_name
                extra_keys = ['composition', 'monte_carlo']
                for key in extra_keys:
                    _ = calc.pop(key, None)
                
                #update file if needed
                outfile = fix_data_file(calc['lattice'], len(calc['element']))
                
                #add ref phase, needed
                calc['reference_phase'] = str(phase_reference_state)
                calc['reference_composition'] = comps['reference']
                calc['mode'] = str('fe')
                calc['folder_prefix'] = f'{phase_name}-{comp:.2f}'
                calc['lattice'] = str(outfile)
                
                #now we need to run this for different temp
                for temp in temp_arr:
                    calc_for_temp = copy.deepcopy(calc)
                    calc_for_temp['temperature'] = int(temp)
                    all_calculations.append(calc_for_temp)
            else:
                #off stoichiometric
                #copy the dict
                calc = copy.deepcopy(phase)

                #first thing first, we need to calculate the number of atoms
                #we follow the convention that composition is always given with the second species
                n_atoms = np.sum(calc['composition']['number_of_atoms'])
                
                #find number of atoms of second species
                output_chemical_composition = {}
                n_species_b = int(np.round(comp*n_atoms, decimals=0))
                output_chemical_composition[reference_element] = n_species_b

                n_species_a = int(n_atoms-n_species_b)
                output_chemical_composition[other_element] = n_species_a

                if n_species_a == 0:
                    raise ValueError("Please add pure phase as a new entry!")
                #create input comp dict and output comp dict
                input_chemical_composition = {element:number for element, number in zip(calc['element'],
                                                                    calc['composition']['number_of_atoms'])}

                #good, now we need to write such a structure out; likely better to use working directory for that
                folder_prefix = f'{phase_name}-{comp:.2f}'
                calc['reference_composition'] = comps['reference']
                #if solid, its very easy; kinda
                #if calc['reference_phase'] == 'solid':
                if use_composition_scaling:
                    #this is solid , and comp scale is turned on
                    #pop extra keys which are not needed
                    #we dont kick out phase_name
                    extra_keys = ['composition', 'reference_phase']
                    for key in extra_keys:
                        _ = calc.pop(key, None)
                    
                    #just submit comp scales
                    #add ref phase, needed
                    calc['mode'] = str('composition_scaling')
                    calc['folder_prefix'] = folder_prefix
                    calc['composition_scaling'] = {}
                    calc['composition_scaling']['output_chemical_composition'] = output_chemical_composition
                    
                else:
                    #manually create a mixed structure - not that the pair style is always ok :)
                    
                    outfile = os.path.join(os.getcwd(), os.path.basename(calc['lattice'])+folder_prefix+'.comp.mod')
                    #print(f'finding comp trf from {input_chemical_composition} to {output_chemical_composition}')
                    #write_structure(calc['lattice'], input_chemical_composition, output_chemical_composition, outfile)
            
                    simplecalc = SimpleCalculation(calc['lattice'], 
                                    calc["element"],
                                    input_chemical_composition,
                                    output_chemical_composition)
                    compsc = CompositionTransformation(simplecalc)
                    compsc.write_structure(outfile)
            
                    #pop extra keys which are not needed
                    #we dont kick out phase name
                    extra_keys = ['composition']
                    for key in extra_keys:
                        _ = calc.pop(key, None)

                    #add ref phase, needed
                    calc['mode'] = str('fe')
                    calc['folder_prefix'] = folder_prefix
                    calc['lattice'] = str(outfile)
            
            #now we need to run this for different temp
            for temp in temp_arr:
                calc_for_temp = copy.deepcopy(calc)
                calc_for_temp['temperature'] = int(temp)
                all_calculations.append(calc_for_temp)       
            
        #finish and write up the file
        output_data = {"calculations": all_calculations}
        for rep in ['.yml', '.yaml']:
            calculation_base_name = calculation_base_name.replace(rep, '')

        outfile_phase = phase_name + '_' + calculation_base_name + ".yaml"
        with open(outfile_phase, 'w') as fout:
            yaml.safe_dump(output_data, fout)
        print(f'Total {len(all_calculations)} calculations found for phase {phase_name}, written to {outfile_phase}')


def _get_temp_arg(tarr, temp, threshold=1E-1):
    if tarr is None:
        return None
    arg = np.argsort(np.abs(tarr-temp))[0]

    th = np.abs(tarr-temp)[arg] 
    if th > threshold:
        arg = None
    return arg

def _is_val_ok(val):
    if val is None:
        return False
    elif math.isnan(val):
        return False
    else:
        return True

def _get_fe_at_args(arr, args):
    fes = []
    for count, x in enumerate(args):
        if _is_val_ok(x):
            fes.append(arr[count][int(x)])
        else:
            fes.append(None)
    return fes

def _calculate_configurational_entropy(x, correction=0):
    if correction == 0:
        s = np.array([(c*np.log(c) + (1-c)*np.log(1-c)) if 1 > c > 0 else 0 for c in x])
    else:
        arg = np.argsort(np.abs(x-correction))[0]
        left_side = x[:arg+1]
        right_side = x[arg:]

        if len(left_side)>0:
            left_side = left_side/left_side[-1]
            s_left = np.array([(c*np.log(c) + (1-c)*np.log(1-c)) if 1 > c > 0 else 0 for c in left_side])
        
        if len(right_side)>0:
            right_side = right_side - right_side[0]
            right_side = right_side/right_side[-1]
            s_right = np.array([(c*np.log(c) + (1-c)*np.log(1-c)) if 1 > c > 0 else 0 for c in right_side])
        
        if len(left_side) == 0:
            return s_right
        elif len(right_side) == 0:
            return s_left
        else:
            return np.concatenate((s_left, s_right[1:]))
    return -s

def _get_free_energy_fit(composition, 
                        free_energy, 
                        fit_order=5,
                        end_weight=3,
                        end_indices=4):
    """
    Create splines for free energy, and return them
    """
    weights = np.ones_like(free_energy)
    weights[0:end_indices] = end_weight
    weights[-end_indices:] = end_weight
    fit = np.polyfit(composition, free_energy, fit_order, w=weights)
    return fit

def get_phase_free_energy(df, phase, temp, 
                          composition_interval=(0, 1),
                          ideal_configurational_entropy=False,
                          entropy_correction=0.0,
                          fit_order=5,
                          composition_grid=10000,
                          composition_cutoff=None,
                          reset_value=1,
                          plot=False,
                          end_weight=3,
                          end_indices=4):
    """
    Get the free energy of a phase as a function of composition.

    Parameters
    ----------
    df: Pandas dataframe
        Dataframe consisting of values from simulation. Should contain at least columns composition, phase, `free_energy` and `temperature`.
        `energy_free` and `temperature` should be arrays of equal length, generally an output from reversible scaling calculation.

    phase: str
        phase for which calculation is to be done. Should be present in `df`.

    temp: float
        temperature at which the free energy curves are to be calculated.

    composition_interval: tuple, optional
        If provided, this composition interval is considered. Default (0, 1)

    ideal_configuration_entropy: bool, optional\
        If True, add the ideal configurational entropy. See Notes. Default False.

    entropy_correction: float, optional.
        The composition of the ordered phase. See Notes. Default None.

    fit_order: int, optional
        Order of the polynomial fit used for fitting free energy as a function of composition. Default 5.

    composition_grid: int, optional
        Number of composition points to be used for fitting. Default 10000.

    composition_cutoff: float, optional
        term for correcting incomplete data. If two consecutive composition values are separated by more than `composition_cutoff`,
        it is reset to `reset_value`. Default None.

    reset_value: float, optional
        see above. Default 1.

    plot: bool, optional
        If True, plot the calculated free energy curves.

    Returns
    -------
    result_dict: dict
        contains keys: "phase", "temperature", "composition", "free_energy", and "entropy".

    Notes
    -----
    To be added
    """
    df_phase = df.loc[df['phase']==phase]
    #drop Nones
    df_phase = df_phase.sort_values(by="composition")
    df_phase = df_phase[(df_phase['composition'] >= composition_interval[0]) & (df_phase['composition'] <= composition_interval[1])]
    
    composition = df_phase['composition'].values
    args = df_phase["temperature"].apply(_get_temp_arg, args=(temp,))
    fes = _get_fe_at_args(df_phase["free_energy"].values, args)
    
    #print(fes)
    #filter out None values
    composition = np.array([composition[count] for count, x in enumerate(fes) if x is not None])
    fes = np.array([x for x in fes if x is not None])

    if (len(fes)==0) or (fes is None):
        warnings.warn("Some temperatures could not be found!")
    else:
        if ideal_configurational_entropy:
            entropy_term = kb*temp*_calculate_configurational_entropy(composition, 
                                                                     correction=entropy_correction) 
            fes = fes - entropy_term
        else:
            entropy_term = []

        fe_fit = _get_free_energy_fit(composition, fes, fit_order=fit_order,
                                            end_weight=end_weight,
                                            end_indices=end_indices)
        compfine = np.linspace(np.min(composition), np.max(composition), composition_grid)
        
        #now fit on the comp grid again
        fe = np.polyval(fe_fit, compfine)

        if composition_cutoff is not None:
            distances = [np.min(np.abs(c-composition)) for c in compfine]
            filters = [x for x in range(len(distances)) if distances[x] > composition_cutoff]
            fe[filters] = reset_value

        if plot:
            plt.scatter(composition, fes, s=4, label=f'{phase}-calc.', color="#e57373")
            plt.plot(compfine, fe, label=f'{phase}-fit', color="#b71c1c")
            plt.xlabel("x")
            plt.ylabel("F (eV/atom)")
            plt.legend()
        
        return {"phase":phase, "temperature": temp, "composition": compfine, 
                "free_energy": fe, "entropy": entropy_term}
    return None


def get_free_energy_mixing(dict_list, threshold=1E-3):
    """
    Input is a list of dictionaries

    Get free energy of mixing by subtracting end member values.
    End members are chosen automatically.
    """
    dict_list = np.atleast_1d(dict_list)

    dict_list = np.array([dct for dct in dict_list if dct is not None])

    #we have to get min_comp from all possible values
    min_comp = np.min([np.min(d["composition"]) for d in dict_list])
    max_comp = np.max([np.max(d["composition"]) for d in dict_list])
    
    #now left ref will be min fe value from all dicts, corresponds to min_comp
    min_fe = []
    max_fe = []
    for d in dict_list:
        diff = np.abs(d["composition"]-min_comp)
        arg = np.argsort(diff)[0]
        if diff[arg] < threshold:
            min_fe.append(d["free_energy"][arg])
        diff = np.abs(d["composition"]-max_comp)
        arg = np.argsort(diff)[0]
        if diff[arg] < threshold:
            max_fe.append(d["free_energy"][arg])
    
    #lists are grabbed, now get the references
    left_ref = np.min(min_fe)
    right_ref = np.min(max_fe)
    
    #print(left_ref, right_ref)
    #now once again, loop through, and add the diff
    for d in dict_list:
        #adjust ref based on composition demands
        scaled_comp = d["composition"]/max_comp
        right_ref_scaled = right_ref*scaled_comp
        left_ref_scaled = left_ref*(1-scaled_comp)
        
        #print(d["free_energy"][-1])
        #print((right_ref_scaled + left_ref_scaled)[-1])
        ref = d["free_energy"] - (right_ref_scaled + left_ref_scaled)
        d["free_energy_mix"] = ref
    return dict_list    

def create_color_list(phases):    
    combinations_list = ['-'.join(pair) for pair in combinations(phases, 2)]
    same_element_pairs = ['-'.join([item, item]) for item in phases]
    final_combinations = same_element_pairs + combinations_list

    color_dict = {}

    color_keys = list(matcolors.keys())
    int_keys = list(matcolors['red'].keys())

    for count, combination in enumerate(final_combinations):
        index = count%len(color_keys)
        second_index = -1
        color_hex = matcolors[color_keys[int(index)]][int_keys[second_index]]
        color_dict[combination] = color_hex
        raw = combination.split('-')
        if raw[0] != raw[1]:
            reversecombo = f'{raw[1]}-{raw[0]}'
            color_dict[reversecombo] = color_hex
    return color_dict

def get_tangent_type(dict_list, tangent, energy):
    left_c = tangent[0]
    right_c = tangent[1]
    
    left_e = energy[0]
    right_e = energy[1]
    
    left_phase = None
    right_phase = None
    
    left_values = []
    left_phases = []
    right_values = []
    right_phases = []

    for d in dict_list:
        diff = np.abs(left_c - d["composition"])
        arg = np.argsort(diff)[0]
        if diff[arg] < 1E-5:
            a = np.abs(left_e - d["free_energy_mix"][arg])
            left_values.append(a)
            left_phases.append(d["phase"])
        diff = np.abs(right_c - d["composition"])
        arg = np.argsort(diff)[0]
        if diff[arg] < 1E-5:
            a = np.abs(right_e - d["free_energy_mix"][arg])
            right_values.append(a)
            right_phases.append(d["phase"])
    
    #now check min values
    left_min_arg = np.argmin(left_values)
    if left_values[left_min_arg] < 1E-5:
        #this is ok
        left_phase = left_phases[left_min_arg]
    
    right_min_arg = np.argmin(right_values)
    if right_values[right_min_arg] < 1E-5:
        #this is ok
        right_phase = right_phases[right_min_arg]
    
    phase_str = f'{left_phase}-{right_phase}'
    return phase_str
            
    
def get_common_tangents(dict_list, 
                        peak_cutoff=0.01, 
                        plot=False, 
                        remove_self_tangents_for=[]):
    """
    Get common tangent constructions using convex hull method
    """
    points = np.vstack([np.column_stack((d["composition"], 
        d["free_energy_mix"])) for d in dict_list]) 
    
    #if color_dict is None:
    #    color_dict = create_color_list(dict_list) 
    
    #make common tangent constructions
    #term checks if two different phases are stable at the end points, then common tangent is needed
    hull = ConvexHull(points)
    convex_points = []
    convex_x = []
    for simplex in hull.simplices:
        ind = points[simplex, 1]<=0.0
        if all(ind):            
            convex_points.extend(points[simplex, 1][ind])
            convex_x.extend(points[simplex, 0][ind])
    
    dist = np.diff(np.sort(convex_x))
    dist = np.where(dist>peak_cutoff)[0]
    sargs = np.argsort(convex_x)
    convex_x = np.array(convex_x)
    convex_points = np.array(convex_points)

    tangents = []
    energies = []    
    tangent_types = []
    phases = []
    
    for d in dist:
        t = [convex_x[sargs][d], convex_x[sargs][d+1]]
        e = [convex_points[sargs][d], convex_points[sargs][d+1]]
        phase_str = get_tangent_type(dict_list, t, e)
        
        remove = False
        ps = phase_str.split('-')
        if ps[0] == ps[1]:
            if ps[0] in remove_self_tangents_for:
                remove = True

        if not remove:
            tangents.append(t)
            energies.append(e)
            tangent_types.append(phase_str)
            phases.append(phase_str.split("-"))
    
    if plot:
        for d in dict_list:
            plt.plot(d["composition"], d["free_energy_mix"], color=colors[np.random.randint(len(colors))])
        for t, e in zip(tangents, energies):
            plt.plot(t, e, color="black", ls="dashed")
        plt.ylim(top=0.0)
    
    return np.array(tangents), np.array(energies), np.array(tangent_types), np.array(phases)


def plot_phase_diagram(tangents, temperature,
    tangent_types,
    phases,
    edgecolor="#37474f",
    linewidth=1,
    linestyle='-'):
    
    #get a phase list
    color_dict = create_color_list(phases) 
    minimal_color_dict = {}
    color_list = []
    for key, val in color_dict.items():
        if val not in color_list:
            color_list.append(val)
            minimal_color_dict[key] = val

    legend_patches = [mpatches.Patch(color=color, label=label) for label, color in minimal_color_dict.items()]

    fig, ax = plt.subplots(edgecolor=edgecolor)

    for count, x in enumerate(tangents):
        for c, a in enumerate(x):
            ax.plot(np.array(a), 
                     [temperature[count], temperature[count]], 
                     linestyle,
                     lw=linewidth,
                     c=color_dict[tangent_types[count][c]],
                     )
    ax.legend(handles=legend_patches, loc='center left', bbox_to_anchor=(1, 0.5))
    return fig

