import numpy as np
from tqdm.notebook import trange
import pandas as pd
import matplotlib.pyplot as plt
import warnings
import itertools
import math

from calphy.integrators import kb

from scipy.spatial import ConvexHull
from scipy.interpolate import splrep, splev

colors = ['#a6cee3','#1f78b4','#b2df8a',
'#33a02c','#fb9a99','#e31a1c',
'#fdbf6f','#ff7f00','#cab2d6',
'#6a3d9a','#ffff99','#b15928']


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
                          plot=False):
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

        fe_fit = _get_free_energy_fit(composition, fes, fit_order=fit_order)
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

def create_color_list(dict_list, color_list=None):
    if color_list is None:
        color_list = colors
    
    phase_list = [d["phase"] for d in dict_list]
    combinations = list(itertools.combinations_with_replacement(phase_list, 2))
    
    if len(combinations) > len(color_list):
        raise ValueError(f'need {len(combinations)} colors, please provide using color_list=')
    
    color_dict = {}
    for count, combo in enumerate(combinations):
        color_dict[f'{combo[0]}-{combo[1]}'] = color_list[count]
        #add reverse just in case
        color_dict[f'{combo[1]}-{combo[0]}'] = color_list[count]
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
                        remove_self_tangents_for=[],
                        color_dict=None):
    """
    Get common tangent constructions using convex hull method
    """
    points = np.vstack([np.column_stack((d["composition"], 
        d["free_energy_mix"])) for d in dict_list]) 
    
    if color_dict is None:
        color_dict = create_color_list(dict_list) 
    
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
    tangent_colors = []
    phases = []
    
    for d in dist:
        t = [convex_x[sargs][d], convex_x[sargs][d+1]]
        e = [convex_points[sargs][d], convex_points[sargs][d+1]]
        phase_str = get_tangent_type(dict_list, t, e)
        
        if phase_str not in remove_self_tangents_for:
            tangents.append(t)
            energies.append(e)
            tangent_colors.append(color_dict[phase_str])
            phases.append(phase_str.split("-"))
    
    if plot:
        for d in dict_list:
            plt.plot(d["composition"], d["free_energy_mix"], color=colors[np.random.randint(len(colors))])
        for t, e in zip(tangents, energies):
            plt.plot(t, e, color="black", ls="dashed")
        plt.ylim(top=0.0)
    
    return np.array(tangents), np.array(energies), np.array(tangent_colors), color_dict, np.array(phases)


def plot_phase_diagram(tangents, temperature,
    colors,
    edgecolor="#37474f",
    linewidth=1,
    linestyle='-'):
    
    fig, ax = plt.subplots(edgecolor=edgecolor)

    for count, x in enumerate(tangents):
        for c, a in enumerate(x):
            ax.plot(np.array(a), 
                     [temperature[count], temperature[count]], 
                     linestyle,
                     lw=linewidth,
                     c=colors[count][c],
                     )
    return fig

