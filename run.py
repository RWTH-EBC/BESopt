# -*- coding: utf-8 -*-
"""
Created on Tue Jun 06 17:51:36 2017

@author: tsz
"""

from __future__ import division

import numpy as np
import parse_inputs
import clustering_medoid as clustering
import building_optimization_decomp as opti
import pickle
import shutil
import time

def _optimize(eco, devs, clustered, par, emissions_max, options):
    max_costs = 1e5
    emis = 0
    opt_force = 0
    mfocus = 0
    
    time_begin = time.time()
    
    device_keys = ["boiler", "chp", "hp"]
    for dev in device_keys:
        list_n = np.sort(devs[dev].keys())
        for n in list_n:
            forced_devs = {dev: n}
            
            (cur_costs, cur_emis) = opti.compute(eco, devs, clustered, par, emissions_max, options, max_costs-0.01, forced_devs, mfocus)
            
            if cur_costs <= max_costs:
                max_costs = cur_costs
                emis = cur_emis
                opt_force = forced_devs
                mfocus = 3
    
    time_end = time.time()
    time_used = time_end - time_begin
    
    return (max_costs, emis, time_used)

def run_building(house="03", a_max=40):
    #%% Read inputs
    raw_inputs = {}
    
    # Scale inputs to kW
    raw_inputs["electricity"] = np.loadtxt("raw_inputs/building_"+house+"/electricity.csv")
    raw_inputs["dhw"]         = np.loadtxt("raw_inputs/building_"+house+"/dhw.csv")
    raw_inputs["sh"]          = np.loadtxt("raw_inputs/building_"+house+"/space_heating.csv")
    raw_inputs["solar_irrad"] = np.loadtxt("raw_inputs/building_"+house+"/solar_rad_35deg.csv") / 1000
    raw_inputs["solar_irrad"] = np.maximum(raw_inputs["solar_irrad"], 0)
    raw_inputs["temperature"] = np.loadtxt("raw_inputs/building_"+house+"/temperature.csv")
    
    design_heat_load = (max(raw_inputs["dhw"]) + max(raw_inputs["sh"])) * 1.2 # in kW
    
    #%% Clustering
    inputs_clustering = np.array([raw_inputs["electricity"], 
                                  raw_inputs["dhw"],
                                  raw_inputs["sh"],
                                  raw_inputs["solar_irrad"],
                                  raw_inputs["temperature"]])
    
    number_clusters = 5
    (inputs, nc, z) = clustering.cluster(inputs_clustering, 
                                         number_clusters=number_clusters,
                                         norm=2,
                                         mip_gap=0.0,
                                         weights=[1,1,1,1,0])
    
    
    # Determine time steps per day
    len_day = int(inputs_clustering.shape[1] / 365)
    
    clustered = {}
    clustered["electricity"] = inputs[0]
    clustered["dhw"]         = inputs[1]
    clustered["sh"]          = inputs[2]
    clustered["solar_irrad"] = inputs[3]
    clustered["temperature"] = inputs[4]
    clustered["design_heat_load"] = design_heat_load
    clustered["weights"]     = nc
    clustered["z"]           = z
    
    #%% Read devices
    devs = parse_inputs.read_devices(timesteps=len_day, days=number_clusters,
                             temperature_ambient=clustered["temperature"],
                             temperature_flow=35, 
                             temperature_design=-12, 
                             solar_irradiation=clustered["solar_irrad"], 
                             days_per_cluster = nc)
    
    (eco, par, devs) = parse_inputs.read_economics(devs)
    par = parse_inputs.compute_parameters(par, number_clusters, len_day)
    
    par["A_max"] = a_max
    par["time_limit"] = 12 * 3600
    
    emissions_max = 1000 # Tons of CO2
    
    #%% Store inputs
    # These inputs (and summarizing results) remain unchanged for all scenarios
    filename = "results/building_"+house+"/inputs.pkl"
    with open(filename, "wb") as f_in:
            pickle.dump("", f_in, pickle.HIGHEST_PROTOCOL)
            pickle.dump("", f_in, pickle.HIGHEST_PROTOCOL)
            pickle.dump(eco, f_in, pickle.HIGHEST_PROTOCOL)
            pickle.dump(devs, f_in, pickle.HIGHEST_PROTOCOL)
            pickle.dump(clustered, f_in, pickle.HIGHEST_PROTOCOL)
            pickle.dump(par, f_in, pickle.HIGHEST_PROTOCOL)
    
    #%% Compute reference result (only boilers + grid)
    options={"filename_results" : "results/building_"+house+"/test_costmin.pkl",
             "EEG": True,
             "KfW": True,
             "KWKG": True,
             "HP tariff": True,
             "nZEB": False,
             "scenario": "free",
             "opt_costs": True,
             "store_start_vals": True,
             "load_start_vals": False,
             "filename_start_vals": "start_values/sv_"+house+".csv"}
    
    (cost_ref, emis_ref, time_normal) = _optimize(eco, devs, clustered, par, emissions_max, options)
    
    options={"filename_results" : "results/building_"+house+"/test_co2red.pkl",
             "EEG": True,
             "KfW": True,
             "KWKG": True,
             "HP tariff": True,
             "nZEB": False,
             "scenario": "free",
             "opt_costs": True,
             "store_start_vals": True,
             "load_start_vals": True,
             "filename_start_vals": "start_values/test_"+house+".csv"}
    
    (cost_red, emis_red, time_reduced) = _optimize(eco, devs, clustered, par, 0.8 * emis_ref, options)

    
    filename = "results/building_"+house+"/time_test.pkl"
    with open(filename, "wb") as f_in:
            pickle.dump(time_normal, f_in, pickle.HIGHEST_PROTOCOL)
            pickle.dump(time_reduced, f_in, pickle.HIGHEST_PROTOCOL)
    
if __name__ == "__main__":
    houses = ("03", "07", "10")
    a_max = {"03": 40,
             "07": 100,
             "10": 170}
    
    for house in houses:
        run_building(house, a_max[house])
