import numpy as np
from tqdm.notebook import trange
import pandas as pd
import matplotlib.pyplot as plt
import warnings
import itertools

from calphy.integrators import kb

from scipy.spatial import ConvexHull
from scipy.interpolate import splrep, splev

colors = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928']

def get_free_energy_at(d, phase, comp, temp, threshold=1E-1):
    """
    Extract free energy at given temperature
    """
    tarr = np.array(d[phase]["%.2f"%comp]["temperature"])
    arg = np.argsort(np.abs(tarr-temp))[0]
    th = np.abs(tarr-temp)[arg] 
    if th > threshold:
        val = None
    else:
        val = d[phase]["%.2f"%comp]["free_energy"][arg] 
    return val
    
def calculate_configurational_entropy(x, correction=0):
    """
    Calculate configurational entropy
    """
    if correction == 0:
        s = np.array([(c*np.log(c) + (1-c)*np.log(1-c)) if 1 > c > 0 else 0 for c in x])
    else:
        arg = np.argsort(np.abs(x-correction))[0]
        left_side = x[:arg+1]
        right_side = x[arg:]
        #print(len(left_side))
        #print(left_side)
        #print(len(right_side))
        #print(right_side)

        if len(left_side)>0:
            left_side = left_side/left_side[-1]
            s_left = np.array([(c*np.log(c) + (1-c)*np.log(1-c)) if 1 > c > 0 else 0 for c in left_side])
            
        #correct to zero
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

#def get_free_energy_splines(composition, free_energy, k=3):
#    """
#    Create splines for free energy, and return them
#    """
#    return splrep(comp, fes, k=3)

def get_free_energy_fit(composition, free_energy, fit_order=5):
    """
    Create splines for free energy, and return them
    """
    weights = np.ones_like(free_energy)
    weights[0:4] = 3
    weights[-4:] = 3
    fit = np.polyfit(composition, free_energy, fit_order, w=weights)
    return fit


def get_phase_free_energy(data, phase, temp, 
                          ideal_configurational_entropy=False,
                          entropy_correction=0.0,
                          composition_grid=10000,
                          fit_order=5,
                          plot=False,
                          composition_interval=(0, 1),
                          composition_cutoff=None,
                          reset_value=1):
    """
    Extract free energy for given phase
    """
    comporg = list(data[phase]["composition"])
    fes = []
    comp = []
    
    for c in comporg:
        if (composition_interval[0] <= c <= composition_interval[1]):
            f = get_free_energy_at(data, phase, float(c), temp)
            if f is not None:
                fes.append(f)
                comp.append(c)
    
    fes = np.array(fes)
    comp = np.array(comp)
    
    if (len(fes)==0) or (fes is None):
        warnings.warn("Some temperatures could not be found!")
    else:  
        if ideal_configurational_entropy:
            entropy_term = kb*temp*calculate_configurational_entropy(comp, correction=entropy_correction) 
            fes = fes - entropy_term
        else:
            entropy_term = []

        fe_fit = get_free_energy_fit(comp, fes, fit_order=fit_order)
        compfine = np.linspace(np.min(comp), np.max(comp), composition_grid)

        fe = np.polyval(fe_fit, compfine)
        
        #fix missing values; assign +0.01 to all values which are not within vicinity
        if composition_cutoff is not None:
            #so we go along composition, see if there are points with no adjacent comp values, ignore them
            distances = [np.min(np.abs(c-comp)) for c in compfine]
            filters = [x for x in range(len(distances)) if distances[x] > composition_cutoff]
            fe[filters] = reset_value
            
        if plot:
            plt.scatter(comp, fes, s=4, label=f'{phase}-calc.', color=colors[np.random.randint(len(colors))])
            plt.plot(compfine, fe, label=f'{phase}-fit', color=colors[np.random.randint(len(colors))])
            plt.xlabel("x")
            plt.ylabel("F (eV/atom)")
            plt.legend()
            #plt.ylim(top=0.0)
        return {"phase":phase, "temperature": temp, "composition": compfine, 
                "free_energy": fe, "entropy": entropy_term}
    return None


def get_free_energy_mixing(dict_list, threshold=1E-3):
    """
    Input is a list of dictionaries
    """
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
            
    
def get_common_tangents(dict_list, peak_cutoff=0.01, plot=False, 
                        remove_self_tangents_for=[],
                        color_dict=None):
    """
    Get common tangent constructions using convex hull method
    """
    points = np.vstack([np.column_stack((d["composition"], d["free_energy_mix"])) for d in dict_list]) 
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
    
    for d in dist:
        t = [convex_x[sargs][d], convex_x[sargs][d+1]]
        e = [convex_points[sargs][d], convex_points[sargs][d+1]]
        phase_str = get_tangent_type(dict_list, t, e)
        
        if phase_str not in remove_self_tangents_for:
            tangents.append(t)
            energies.append(e)
            tangent_colors.append(color_dict[phase_str])
    
    if plot:
        for d in dict_list:
            plt.plot(d["composition"], d["free_energy_mix"], color=colors[np.random.randint(len(colors))])
        for t, e in zip(tangents, energies):
            plt.plot(t, e, color="black", ls="dashed")
        plt.ylim(top=0.0)
    return np.array(tangents), np.array(energies), np.array(tangent_colors), color_dict

