#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 01 10:35:56 2015

@author: tsz
"""

from __future__ import division

import gurobipy as gp
import numpy as np
import math
import pickle

def compute(eco, devs, demands, params, max_emi, options, max_costs, forced_devs, mfocus=0):
    """
    Compute the optimal building energy system consisting of pre-defined 
    devices (devs) for a given building (demands) under given economic and 
    other parameters (eco and params).
    """
    
    # Extract parameters
    dt = params["dt"]
    
    # Define subsets
    heater = ("boiler", "chp", "eh", "hp")
    heater_ehdhw = ("boiler", "chp", "eh", "hp", "eh_dhw")
    storage = ("bat", "tes")
    solar = ("pv", "stc")
    
    time_steps = range(params["time_steps"])
    days = range(params["days"])

#%%    Start of the Model: Declarations of the model's variables
    try:
        model = gp.Model("Design computation")
        
        # Define variables
        # Costs and Revenues
        c_inv   = {dev: model.addVar(vtype="C", name="c_inv_"+dev)
                for dev in devs.keys()}
        c_om    = {dev: model.addVar(vtype="C", name="c_om_"+dev)
                for dev in devs.keys()}
        c_dem   = {dev: model.addVar(vtype="C", name="c_dem_"+dev)
                for dev in ("boiler", "chp", "grid_hou", "grid_hp")}
        c_fix   = {dev: model.addVar(vtype="C", name="c_fix_"+dev)
                for dev in ("el", "gas")}
        c_eeg   = {dev: model.addVar(vtype="C", name="c_eeg_"+dev)
                for dev in ("pv", "chp")}
        revenue = {dev: model.addVar(vtype="C", name="revenue_"+dev)
                for dev in ("chp", "pv")}
        subsidy = {dev: model.addVar(vtype="C", name="subsidy_"+dev)
                for dev in ("chp","bat")}
        sub     = {dev: model.addVar(vtype="C", name="sub_"+dev)
                for dev in ("micro","large")}

        # chp subsidy variables
        x_chp = {}
        for dev in ("lump","var"):
            sub[dev,"micro"]    = model.addVar(vtype="C", name="sub_"+str(dev)+"_micro")
            x_chp[dev,"micro"]  = model.addVar(vtype="B", name="x_chp_"+str(dev)+"_micro")

        x_eeg_levy = model.addVar(vtype="B", name="x_eeg_levy")
        x_eeg_10 = model.addVar(vtype="B", name="x_eeg_10")
        x_eeg_40 = model.addVar(vtype="B", name="x_eeg_40")
        x_eeg_100 = model.addVar(vtype="B", name="x_eeg_100")
        c_eeg_lin = {dev: model.addVar(vtype="C", name="c_eeg_lin_"+dev)
                for dev in ("pv", "chp")}
        rev_pv_10 = model.addVar(vtype="C", name="rev_pv_10")
        rev_pv_40 = model.addVar(vtype="C", name="rev_pv_40")
        rev_pv_100 = model.addVar(vtype="C", name="rev_pv_100")
        p_sell_total_pv = model.addVar(vtype="C", name="p_sell_total_pv")

        # SOC, power, heat and energy
        soc = {}
        power = {}
        heat = {}
        energy = {}
        for d in days: # All days
            for t in time_steps: # All time steps of all days
                timetag = "_"+str(d)+"_"+str(t)
                for dev in storage: # All storage devices
                    for n in devs[dev].keys(): # All listed types of each dev
                        soc[dev,n,d,t] = model.addVar(vtype="C",
                                         name="SOC_"+dev+"_"+str(n)+timetag,
                                         lb=0, ub=1)
                
                for dev in (heater_ehdhw+solar):
                    for n in devs[dev].keys():
                        power[dev,n,d,t] = model.addVar(vtype="C",
                                           name="P_"+dev+"_"+str(n)+timetag)
                        heat[dev,n,d,t] = model.addVar(vtype="C",
                                           name="Q_"+dev+"_"+str(n)+timetag)

                for dev in heater_ehdhw:
                    for n in devs[dev].keys():
                        energy[dev,n,d,t] = model.addVar(vtype="C",
                                            name="E_"+dev+"_"+str(n)+timetag)

        # Weights for the interpolation, data type depends on mudulation mode
        lin = {}
        number_nodes = {}
        number_nodes_hp = {}
        for dev in ("boiler", "chp", "hp"):
            if dev in ("chp", "boiler"):
                for n in devs[dev].keys():
                    number_nodes[dev,n] = len(devs[dev][n]["Q_heat"])
                    for d in days:
                        for t in time_steps:
                            for i in range(number_nodes[dev,n]):
                                tag = "_"+str(dev)+"_"+str(n)+"_"+str(d)+"_"+str(t)+"_"+str(i)
                                lin[dev,n,d,t,i] = model.addVar(vtype="C", name="lin_"+tag)
                
            dev = "hp"
            for n in devs[dev].keys():
                m= devs[dev][n]["Mod"]
                if m == 'stageless':
                    for d in days: # All days
                        for t in time_steps: # All time steps of all days
                            number_nodes_hp[dev,n,d,t] = len(devs[dev][n]["Q_heat"][d,t])
                            for i in range(number_nodes_hp[dev,n,d,t]):
                                tag = "_"+str(dev)+"_"+str(n)+"_"+str(d)+"_"+str(t)+"_"+str(i)
                                lin[dev,n,d,t,i] = model.addVar(vtype="C", name="lin_"+tag)
                if m == 'stage':
                    for d in days: # All days
                        for t in time_steps: # All time steps of all days
                            number_nodes_hp[dev,n,d,t] = len(devs[dev][n]["Q_heat"][d,t])
                            for i in range(number_nodes_hp[dev,n,d,t]):
                                tag = "_"+str(dev)+"_"+str(n)+"_"+str(d)+"_"+str(t)+"_"+str(i)
                                lin[dev,n,d,t,i] = model.addVar(vtype="B", name="lin_"+tag)

        # Storage initial SOC's
        soc_init = {}
        for dev in storage:
            for n in devs[dev].keys():
                for d in days:
                    tag = dev + "_" + str(n) + "_" + str(d)
                    soc_init[dev,n,d] = model.addVar(vtype="C", lb=0, ub=1,
                                                     name="SOC_init_"+tag)

        # Storage charging and discharging
        ch = {}
        dch = {}
        for dev in storage:
            for n in devs[dev].keys():
                for d in days:
                    for t in time_steps:
                        timetag = "_"+str(d)+"_"+str(t)

                        ch[dev,n,d,t] = model.addVar(vtype="C",
                                        name="ch_"+dev+"_"+str(n)+timetag)
                        dch[dev,n,d,t] = model.addVar(vtype="C",
                                         name="dch_"+dev+"_"+str(n)+timetag)
                    
        # Electricity imports, sold, self-used and transferred (to heat pump) electricity
        p_grid  = {}
        p_use   = {}
        p_sell  = {}
        p_hp    = {}
        for d in days:
            for t in time_steps:
                timetag = "_"+str(d)+"_"+str(t)
                
                p_grid["grid_hou",d,t] = model.addVar(vtype="C", name="p_grid_hou"+timetag)
                p_grid["grid_hp",d,t]  = model.addVar(vtype="C", name="p_grid_hp"+timetag)
                
                # Note: bat is referring to the discharge power
                for dev in ("pv", "bat"):
                    p_use[dev,d,t]  = model.addVar(vtype="C", 
                                                  name="P_use_"+dev+timetag)
                    p_sell[dev,d,t] = model.addVar(vtype="C", 
                                                   name="P_sell_"+dev+timetag)
                    p_hp[dev,d,t]   = model.addVar(vtype="C", 
                                                   name="P_hp_"+dev+timetag)
                for n in devs["chp"].keys():
                    dev = "chp"
                    p_hp[dev,n,d,t] = model.addVar(vtype="C", name="P_hp_"+dev+"_"+str(n)+timetag)
                    p_use[dev,n,d,t] = model.addVar(vtype="C", name="P_use_"+dev+"_"+str(n)+timetag)
                    p_sell[dev,n,d,t] = model.addVar(vtype="C", name="P_sell_"+dev+"_"+str(n)+timetag)

        # Amount of gas consumed
        gas_tariffs = eco["gas"].keys()
        gas_tariffs.remove("energy_tax")
        G = {}
        for tar in gas_tariffs:
            G[tar] = {}
            for dev in ("boiler","chp"):
                for n in xrange(len(eco["gas"][tar]["lb"])):
                    G[tar][dev,n] = model.addVar(vtype="C", name="G_"+dev+"_"+str(n))
        G_total = {dev: model.addVar(vtype="C", name="G_total_"+dev)
                  for dev in ("boiler","chp")}
        
        # Amount of electricity consumed
        el_tariffs = eco["el"].keys()
        El = {}
        for tar in el_tariffs:
            El[tar] = {}
            for dev in ("grid_hou","grid_hp"):
                for n in xrange(len(eco["el"][tar]["lb"])):
                    El[tar][dev,n] = model.addVar(vtype="C", name="El_"+dev+"_"+str(n))
        El_total = {dev: model.addVar(vtype="C", name="El_total_"+dev)
                    for dev in ("grid_hou","grid_hp")}
                                       
        # Split EH for HP tariff
        eh_split = {}
        for d in days:
            for t in time_steps:
                timetag = "_"+str(d)+"_"+str(t)
                
                eh_split["eh_w/o_hp",d,t] = model.addVar(vtype="C", name="p_eh_w/o_hp"+timetag)
                eh_split["eh_w/_hp",d,t]  = model.addVar(vtype="C", name="p_eh_w/_hp"+timetag)
        
        lim_ch = {}
        lim_sell = {}
        for d in days:
            for t in time_steps:
                lim_ch[d,t] = model.addVar(vtype="B", name="lim_ch"+"_"+str(d)+"_"+str(t))
                lim_sell[d,t] = model.addVar(vtype="B", name="lim_sell"+"_"+str(d)+"_"+str(t))
       
        # Activation and purchase decision variables
        x = {}  # Purchase (all devices)
        y = {}  # Activation (heaters)
        z = {}  # Number of modules (PV, STC)
        x_tariff = {"gas":{}, "el":{}}   # Tariffs
        
        for dev in devs.keys():
            for n in devs[dev].keys():
                x[dev,n] = model.addVar(vtype="B", name="x_"+dev+"_"+str(n))
        
        # All tariff gradations
        for tar in gas_tariffs:
            x_tariff["gas"][tar] = {}
            for n in range(len(eco["gas"][tar]["lb"])):            
                x_tariff["gas"][tar][n] = model.addVar(vtype="B", name="x_tariff_"+tar+"_"+str(n))
        for tar in el_tariffs:
            x_tariff["el"][tar] = {}
            for n in range(len(eco["el"][tar]["lb"])):            
                x_tariff["el"][tar][n] = model.addVar(vtype="B", name="x_tariff_"+tar+"_"+str(n))
        
        # General tariff decision variables
        x_gas = {tar: model.addVar(vtype="B", name="x_"+tar)
                for tar in gas_tariffs}
        x_el = {tar: model.addVar(vtype="B", name="x_"+tar)
                for tar in el_tariffs}

        for d in days:
            for t in time_steps:
                timetag = "_"+str(d)+"_"+str(t)
                for dev in heater_ehdhw: # All heating devices
                    for n in devs[dev].keys(): # All listed types of each device
                        y[dev,n,d,t] = model.addVar(vtype="B",
                                       name="y_"+dev+"_"+str(n)+timetag)
        for dev in solar: # All solar devices
            for n in devs[dev].keys():
                max_modules = params["A_max"] / devs[dev][n]["area"]
                z[dev,n] = model.addVar(vtype="I", name="z_"+dev+"_"+str(n),
                                        lb=0, ub=math.floor(max_modules))
        
        y_stc = {}
        for d in days:
            for t in time_steps:
                timetag = "_"+str(d)+"_"+str(t)
                y_stc[d,t] = model.addVar(vtype="B", name="y_stc"+timetag)
                                
        # Variables for linearization
        lin_pv_bat = {}
        for n in devs["pv"].keys():
            lin_pv_bat[n] = model.addVar(vtype="I", name="lin_pv_bat_"+str(n))
        lin_bat_sub = model.addVar(vtype="C", name="lin_bat_sub")
        
        # Definition of emissions over all devices that are used
        emission = model.addVar(vtype="C", name= "CO2_emission", lb=-gp.GRB.INFINITY)
        
        #######################################################
        # Define total costs
        #####################################################
        
        c_total =   (sum(c_inv[key]  for key in c_inv.keys())
                   + sum(c_om[key]    for key in c_om.keys())
                   + sum(c_dem[key]   for key in c_dem.keys())
                   + sum(c_fix[key]   for key in c_fix.keys())
                   + sum(c_eeg[key]   for key in c_eeg.keys())
                   - sum(revenue[key] for key in revenue.keys())
                   - sum(subsidy[key] for key in subsidy.keys()))
#%%     Loading variables into the model, definition of objective        
        # Update
        model.update()

    ###############################################################################    
    # Objective: Minimize investments, service, demand, metering costs (less 
    #   generated revenues, or minimize emissions)

        if options["opt_costs"]:
            model.setObjective(c_total, gp.GRB.MINIMIZE)
        else:
            model.setObjective(emission, gp.GRB.MINIMIZE)
           
    ###############################################################################               
        
        model.addConstr(c_total <= max_costs)

#%%     Economic constraints: costs, revenues, subsidies
        # Economic constraints
        # Investment costs
        for dev in ("boiler", "chp", "hp", "eh", "eh_dhw", "tes", "bat", "inv"):
            model.addConstr(c_inv[dev] == eco["crf"] * eco["tax"] * 
                sum(x[dev,n] * (1-devs[dev][n]["rval"]) * devs[dev][n]["c_inv"]
                    for n in devs[dev].keys()),
                name="Investment_costs_"+dev)
        for dev in ("pv", "stc"):
            model.addConstr(c_inv[dev] == eco["crf"] * eco["tax"] * 
                sum(z[dev,n] * (1-devs[dev][n]["rval"]) * devs[dev][n]["c_inv"]
                    for n in devs[dev].keys()),
                name="Investment_costs_"+dev)

        # Operation and maintenance
        for dev in ("boiler", "hp", "eh", "eh_dhw", "tes", "bat", "inv"):
            model.addConstr(c_om[dev] == eco["b"]["infl"] * eco["crf"] 
                * eco["tax"] * 
                sum(x[dev,n] * devs[dev][n]["c_om"] for n in devs[dev].keys()),
                name="O&M_costs_"+dev)
        dev = "chp"
        model.addConstr(c_om[dev] == eco["b"]["infl"] * eco["crf"] * eco["tax"] * dt *
                sum(sum(sum(power[dev,n,d,t]
                            for t in time_steps)
                        for d in days)
                    * devs[dev][n]["c_om"] for n in devs[dev].keys()),
                name="O&M_costs_"+dev)
        for dev in ("pv", "stc"):
            model.addConstr(c_om[dev] ==  eco["b"]["infl"] * eco["crf"]
                * eco["tax"] *
                sum(z[dev,n] * devs[dev][n]["c_om"] for n in devs[dev].keys()),
                name="O&M_costs_"+dev)
        
        # Demand related costs (gas)
        # Exactly one gas tariff if at least one chp or boiler is installed
        for dev in ("chp", "boiler"):
            model.addConstr(sum(x_gas[key] for key in x_gas.keys()) >= 
                            sum(x[dev,n] for n in devs[dev].keys()),
                            name="single_gas_tariff_"+dev)
            
        model.addConstr(sum(x_gas[key] for key in x_gas.keys()) <= 1,
                        name="single_gas_tariff_overall")

        for tar in gas_tariffs:
            # Choose one tariff level for the dertermined gas tariff
            model.addConstr(sum(x_tariff["gas"][tar][n] for n in x_tariff["gas"][tar].keys()) == x_gas[tar],
                            name="gas_levels_"+tar+"_"+str(n))
            # The tariff level is restricted by the consumed gas amount
            for n in x_tariff["gas"][tar].keys():
                model.addConstr(x_tariff["gas"][tar][n] * eco["gas"][tar]["lb"][n] <= 
                                (G[tar]["boiler",n] + G[tar]["chp",n]) * 0.001,
                                name="gas_level_lb"+tar+"_"+str(n))
                model.addConstr(x_tariff["gas"][tar][n] * eco["gas"][tar]["ub"][n] >= 
                                (G[tar]["boiler",n] + G[tar]["chp",n]) * 0.001,
                                name="gas_level_ub"+tar+"_"+str(n))
                
        # Divide because of energy tax
        for dev in ("boiler","chp"):   
            # Total amount of gas used
            model.addConstr(G_total[dev] == sum(sum(G[tar][dev,n] 
                                            for n in x_tariff["gas"][tar].keys())
                                            for tar in gas_tariffs))
            model.addConstr(G_total[dev] == sum(demands["weights"][d] * dt * sum(sum(energy[dev,n,d,t] 
                                    for n in devs[dev].keys()) 
                                    for t in time_steps)
                                    for d in days))
            # Variable gas costs
            if dev == "chp":
                model.addConstr(c_dem[dev] == sum(sum(
                        G[tar][dev,n] * (eco["gas"][tar]["var"][n] - eco["gas"]["energy_tax"])
                        for n in x_tariff["gas"][tar].keys())
                        for tar in gas_tariffs) * eco["b"]["gas"] * eco["crf"],
                        name="c_dem_"+dev)
            else: 
                model.addConstr(c_dem[dev] == sum(sum(
                        G[tar][dev,n] * eco["gas"][tar]["var"][n]
                        for n in x_tariff["gas"][tar].keys())
                        for tar in gas_tariffs) * eco["b"]["gas"] * eco["crf"],
                        name="c_dem_"+dev)
                            
        # Fixed costs for gas administration
        model.addConstr(c_fix["gas"] == sum(sum(x_tariff["gas"][tar][n] * eco["gas"][tar]["fix"][n]
                        for n in x_tariff["gas"][tar].keys())
                        for tar in gas_tariffs),
                        name="c_fix_gas")
        
        # Demand related costs (electricity)
        # Choose one tariff for general electricity purchase
        non_hp_tariffs = sum(x_el[tar] for tar in x_el.keys() if eco["el"][tar]["hp"] == 0)
        model.addConstr(non_hp_tariffs == 1,
                        name="single_el_tariff")

        # If a HP is installed, the HP tariff is available
        hp_tariffs = sum(x_el[tar] for tar in x_el.keys() if eco["el"][tar]["hp"] == 1)
        if options["HP tariff"]:
            # Allow special heat pump tariffs
            model.addConstr(hp_tariffs <= sum(x["hp",n] for n in devs["hp"].keys()),
                            name="optional_hp_tariff")
#            model.addConstr(hp_tariffs == sum(x["hp",n] for n in devs["hp"].keys()),
#                            name="optional_hp_tariff")
        else:
            # Prohibit special heat pump tariffs
            model.addConstr(hp_tariffs <= 0, name="optional_hp_tariff")

        # grid_hou electricity cannot be purchased with el_hp tariff
        for tar in x_el.keys():
            if eco["el"][tar]["hp"] == 1:
                for n in x_tariff["el"][tar].keys():
                    model.addConstr(El[tar]["grid_hou",n] == 0)
        
        for tar in el_tariffs:
            # Choose one tariff level for the dertermined el tariff
            model.addConstr(sum(x_tariff["el"][tar][n] for n in x_tariff["el"][tar].keys()) == x_el[tar])
            # The tariff level is restricted by the consumed el amount
            for n in x_tariff["el"][tar].keys():
                model.addConstr(x_tariff["el"][tar][n] * eco["el"][tar]["lb"][n] <= 
                (El[tar]["grid_hou",n] + El[tar]["grid_hp",n]) * 0.001,
                name="el_level_lb_"+tar+"_"+str(n))
                model.addConstr(x_tariff["el"][tar][n] * eco["el"][tar]["ub"][n] >= 
                (El[tar]["grid_hou",n] + El[tar]["grid_hp",n]) * 0.001,
                name="el_level_ub_"+tar+"_"+str(n))
            
        # Devide because of optional HP tariff
        for dev in ("grid_hou", "grid_hp"):
            # Total amount of electricity used
            model.addConstr(El_total[dev] == sum(sum(El[tar][dev,n] 
                                            for n in x_tariff["el"][tar].keys())
                                            for tar in el_tariffs))
            model.addConstr(El_total[dev] == sum(demands["weights"][d] * sum(p_grid[dev,d,t] 
                                            for t in time_steps)
                                            for d in days) * dt)
        
            # Variable costs
            model.addConstr(c_dem[dev] == sum(sum(
                    El[tar][dev,n] * eco["el"][tar]["var"][n]
                    for n in x_tariff["el"][tar].keys())
                    for tar in el_tariffs) * eco["b"]["el"] * eco["crf"],
                    name="c_dem_"+dev)
        
        # Fixed costs for el administration
        model.addConstr(c_fix["el"] == sum(sum(x_tariff["el"][tar][n] * eco["el"][tar]["fix"][n]
                        for n in x_tariff["el"][tar].keys())
                        for tar in el_tariffs),
                        name="c_fix_el")
        
        # EEG-Levy
        M = 10000
        for dev in ("pv", "chp"):
            model.addConstr(c_eeg[dev] <= x_eeg_levy * M)
            
            model.addConstr(c_eeg_lin[dev] - c_eeg[dev] >= 0)
            model.addConstr(c_eeg_lin[dev] - c_eeg[dev] <= (1 - x_eeg_levy) * M)

        dev = "pv"
        model.addConstr(c_eeg_lin[dev] == eco["eeg_levy"] * dt * 
                                      sum(demands["weights"][d] * sum(p_use[dev,d,t] + p_hp[dev,d,t]
                                                                  for t in time_steps) 
                                      for d in days),
                        "eeg levy " + dev)

        dev = "chp"
        model.addConstr(c_eeg_lin[dev] == eco["eeg_levy"] * dt * 
                                      sum(demands["weights"][d] * sum(sum(p_use[dev,n,d,t] + p_hp[dev,n,d,t]
                                                                      for n in devs[dev].keys())
                                                                  for t in time_steps)
                                      for d in days),
                        "eeg levy " + dev)

        # EEG levy for PV
        model.addConstr(sum(devs["pv"][n]["p_nom"] * z["pv",n] for n in devs["pv"].keys()) <= 
                        10 + 10000 * x_eeg_levy)
        # EEG levy for CHP
        model.addConstr(sum(x["chp",n] * devs["chp"][n]["p_nom"] for n in devs["chp"].keys()) <= 
                        10 + 10000 * x_eeg_levy)

        # Revenues from sold electricity
        dev = "pv"
        if options["EEG"]:
            # Remuneration as with EEG
            model.addConstr(sum(devs[dev][n]["p_nom"] * z[dev,n] for n in devs[dev].keys())
                         <= 10 * x_eeg_10 + 40 * x_eeg_40 + 750 * x_eeg_100)
            model.addConstr(x_eeg_10 + x_eeg_40 + x_eeg_100 <= sum(x["pv",n] for n in devs["pv"].keys()))

            model.addConstr(p_sell_total_pv == sum(demands["weights"][d] * sum(p_sell[dev,d,t] 
                                            for t in time_steps)
                                            for d in days) * dt)

            model.addConstr(p_sell_total_pv == rev_pv_10 + rev_pv_40 + rev_pv_100)

            model.addConstr(revenue[dev] == eco["b"]["eex"] * eco["crf"] *
                    (rev_pv_10 * eco["sell",dev,"10"] + 
                     rev_pv_40 * eco["sell",dev,"40"] +
                     rev_pv_100 * eco["sell",dev,"100"]),
                    name="Feed_in_rev_"+dev)
            
            M = 1000
            model.addConstr(1.0 / M * rev_pv_10 <= x_eeg_10 * M)
            model.addConstr(1.0 / M * rev_pv_40 <= x_eeg_40 * M)
            model.addConstr(1.0 / M * rev_pv_100 <= x_eeg_100 * M)

        else:
            # Market conform remuneration (similar to CHP)
            model.addConstr(revenue[dev] == eco["b"]["eex"] * eco["crf"] * dt *
                sum(demands["weights"][d] *
                    sum(p_sell[dev,d,t] * eco["sell","chp"] for t in time_steps)
                for d in days),
                name="Feed_in_rev_"+dev)
            
        dev = "chp"
        model.addConstr(revenue[dev] == eco["b"]["eex"] * eco["crf"] * dt *
            sum(demands["weights"][d] *
                sum(sum(p_sell[dev,n,d,t] for n in devs[dev].keys()) * 
                eco["sell",dev] for t in time_steps)
            for d in days),
            name="Feed_in_rev_"+dev)

        # Subsidies
        # CHP
        dev = "chp"
        if options["KWKG"]:
            # Divide into micro and bigger sized chp units
            model.addConstr(subsidy[dev] == sub["micro"] + sub["large"])
            p_micro = 2 # kW
            
            # Subsidies for bigger sized chp units
            model.addConstr(sub["large"] == eco["crf"] * eco["b"]["eex"] * dt *
                            sum(sum(demands["weights"][d] * sum(
                            eco["sub_chp"]["self"] * (p_use[dev,n,d,t] + p_hp[dev,n,d,t]) +
                            eco["sub_chp"]["sell"] * p_sell[dev,n,d,t] for t in time_steps) for d in days) 
                            for n in devs[dev].keys() if devs[dev][n]["p_nom"] > p_micro),
                            name="sub_chp_large")
                            
            # Subsidies for micro chp units
            # For micro chp units either lump or variable subsidies (maximum) 
            # according to http://www.fico.com/en/node/8140?file=5125 p. 7
            max_lump = eco["crf"] * eco["sub_chp"]["t"] * eco["sub_chp"]["lump"] * p_micro
            max_var  = eco["crf"] * eco["b"]["eex"] * 8760 * eco["sub_chp"]["sell"] * p_micro
            
            model.addConstr(sub["lump","micro"] == eco["crf"] * eco["sub_chp"]["t"] * eco["sub_chp"]["lump"] *
                            sum(x[dev,n] * devs[dev][n]["p_nom"] 
                            for n in devs[dev].keys() if devs[dev][n]["p_nom"] <= p_micro),
                            name="sub_chp_micro_lump")
            model.addConstr(sub["var","micro"]  == eco["crf"] * eco["b"]["eex"] * dt *
                            sum(sum(demands["weights"][d] * sum(
                            eco["sub_chp"]["self"] * (p_use[dev,n,d,t] + p_hp[dev,n,d,t]) +
                            eco["sub_chp"]["sell"] * p_sell[dev,n,d,t] for t in time_steps) for d in days) 
                            for n in devs[dev].keys() if devs[dev][n]["p_nom"] <= p_micro),
                            name="sub_chp_micro_var")
            
            model.addConstr(sub["micro"] >= sub["lump","micro"])
            model.addConstr(sub["micro"] >= sub["var","micro"])
            model.addConstr(sub["micro"] <= sub["lump","micro"] + max_lump * (1 - x_chp["lump","micro"]))
            model.addConstr(sub["micro"] <= sub["var","micro"]  + max_var * (1 - x_chp["var","micro"]))
#            model.addConstr(x_chp["lump","micro"] + x_chp["var","micro"] == sum(x[dev,n] 
#                            for n in devs[dev].keys() if devs[dev][n]["p_nom"] <= p_micro))
            model.addConstr(x_chp["lump","micro"] + x_chp["var","micro"] == 1)
        else:
            model.addConstr(subsidy[dev] == 0)

        # Battery
        # abbrev.
        if options["KfW"]:
            pv_pow = sum(devs["pv"][n]["p_nom"]*z["pv",n] for n in devs["pv"].keys())
            model.addConstr(subsidy["bat"] <= eco["crf"] * eco["sub_bat_max"] * lin_bat_sub * 0.19,
                            name="Bat_Subsidies_1")
            model.addConstr(subsidy["bat"] <= 
                            (c_inv["pv"] + c_inv["bat"] - eco["crf"] * eco["sub_bat"] * lin_bat_sub) * 0.19,
                            name="Bat_Subsidies_2")
            # linearization constraints because of pv_pow * x_bat
            M = max(devs["pv"][n]["p_nom"] * params["A_max"] / devs["pv"][n]["area"] for n in devs["pv"].keys())
            model.addConstr(lin_bat_sub <= M * sum(x["bat",n] for n in devs["bat"].keys()))
            model.addConstr(pv_pow - lin_bat_sub >= 0)
            model.addConstr(pv_pow - lin_bat_sub <= (1 - sum(x["bat",n] for n in devs["bat"].keys())) * M)
        else:
            model.addConstr(subsidy["bat"] == 0)
#%%     Technical constraints: logical, operation, balances               
        # Technical constraints
        # Maximum of one device per type
        for dev in devs.keys():
            model.addConstr(sum(x[dev,n] for n in devs[dev].keys()) <= 1,
                            name="One_dev_per_type_"+dev)
        
        # Maximum area constraints
        for dev in solar:
            for n in devs[dev].keys():
                max_modules = params["A_max"] / devs[dev][n]["area"]
                model.addConstr(z[dev,n] <= x[dev,n] * math.floor(max_modules),
                                name="Max_specific_area_"+dev+"_"+str(n))
        # Shared roof area
        model.addConstr(sum(sum(z[dev,n] * devs[dev][n]["area"] 
                            for n in devs[dev].keys()) 
                        for dev in solar) <= params["A_max"],
                        name="Maximum_total_area")
        
        # Devices can be switched on only if they have been purchased
        for dev in heater_ehdhw:
            for n in devs[dev].keys():
                model.addConstr(sum( sum(y[dev,n,d,t] for t in time_steps)
                                for d in days) 
                      <= params["time_steps"] * params["days"] * x[dev,n],
                      name="Activation_"+dev+"_"+str(n))
        
        # Devices operation
        # Heat output between mod_lvl*Q_nom and Q_nom (P_nom for heat pumps)
        # Power and Energy directly result from Heat output
        for dev in heater_ehdhw:
            for n in devs[dev].keys():
                for d in days:
                    for t in time_steps:
                        # Abbreviations
                        timetag = "_"+str(d)+"_"+str(t)
                    
                        if dev == "hp":
                            tag = str(dev) + "_" + str(n) + timetag
                            model.addConstr(heat[dev,n,d,t] == sum(lin[dev,n,d,t,i] * devs[dev][n]["Q_heat"][d,t][i]
                                                               for i in range(number_nodes_hp[dev,n,d,t])),
                                                                "Q_th interpolation "+tag)
                            model.addConstr(power[dev,n,d,t] == sum(lin[dev,n,d,t,i] * devs[dev][n]["P_el"][d,t][i]
                                                                for i in range(number_nodes_hp[dev,n,d,t])),
                                                                "P_el interpolation "+tag)
                            model.addConstr(energy[dev,n,d,t] == sum(lin[dev,n,d,t,i] * devs[dev][n]["Q_gas"][d,t]
                                                                for i in range(number_nodes_hp[dev,n,d,t])),
                                                                "Q_gas interpolation "+tag)
                            model.addConstr(y[dev,n,d,t] == sum(lin[dev,n,d,t,i] for i in range(number_nodes_hp[dev,n,d,t])),
                                                                "activation "+tag)
                            model.addSOS(gp.GRB.SOS_TYPE2, [lin[dev,n,d,t,i] for i in range(number_nodes_hp[dev,n,d,t])])

                        elif dev in ("boiler", "chp"):
                            # Add constraints
                            # The constraints force an interpolation that computes the electricity 
                            # generation and gas consumption for a given heat output
                            # for key in values

                            tag = str(dev) + "_" + str(n) + timetag
                            model.addConstr(heat[dev,n,d,t] == sum(lin[dev,n,d,t,i] * devs[dev][n]["Q_heat"][i]
                                                               for i in range(number_nodes[dev,n])),
                                                                "Q_th interpolation "+tag)
                            model.addConstr(power[dev,n,d,t] == sum(lin[dev,n,d,t,i] * devs[dev][n]["P_el"][i]
                                                                for i in range(number_nodes[dev,n])),
                                                                "P_el interpolation "+tag)
                            model.addConstr(energy[dev,n,d,t] == sum(lin[dev,n,d,t,i] * devs[dev][n]["Q_gas"][i]
                                                               for i in range(number_nodes[dev,n])),
                                                               "Q_gas interpolation "+tag)
                            model.addConstr(y[dev,n,d,t] == sum(lin[dev,n,d,t,i] for i in range(number_nodes[dev,n])),
                                                                "activation "+tag)
                            model.addSOS(gp.GRB.SOS_TYPE2, [lin[dev,n,d,t,i] for i in range(number_nodes[dev,n])])
                                
                        elif dev in ("eh", "eh_dhw"):
                            mod_lvl = devs[dev][n]["mod_lvl"]
                            eta = devs[dev][n]["eta"][d,t]
                            omega = devs[dev][n]["omega"][d,t]
                            #get maximum and minimum heat option 
                            q_nom = devs[dev][n]["Q_nom"]
                            model.addConstr(heat[dev,n,d,t] 
                                <= q_nom * y[dev,n,d,t],
                                name="Max_heat_output_"+dev+"_"+str(n)+timetag)
                            model.addConstr(heat[dev,n,d,t]
                                >= q_nom * mod_lvl * y[dev,n,d,t],
                                name="Min_heat_output_"+dev+"_"+str(n)+timetag)
                            #interpolate the power with the heat
                            model.addConstr(power[dev,n,d,t]
                                == 1/eta * heat[dev,n,d,t],
                                name="Power_equation_"+dev+"_"+str(n)+timetag)
                            #interpolate the gas energy with the heat (or the power)        
                            model.addConstr(energy[dev,n,d,t] == 0,
                                name="Energy_equation_"+dev+"_"+str(n)+timetag)
        
                                
        # Solar components
        eta_inv = np.mean([devs["inv"][n]["eta"] for n in devs["inv"].keys()])
        for dev in solar:
            for n in devs[dev].keys():
                a = devs[dev][n]["area"]
                for d in days:
                    for t in time_steps:
                        timetag = "_"+str(d)+"_"+str(t)
                        eta_th = devs[dev][n]["eta_th"][d][t]
                        eta_el = devs[dev][n]["eta_el"][d][t]
                        solar_irrad = demands["solar_irrad"][d][t]

                        model.addConstr(heat[dev,n,d,t] <= 
                              eta_th * a * z[dev,n] * solar_irrad,
                              name="Solar_thermal_"+dev+"_"+str(n)+timetag)
                        model.addConstr(power[dev,n,d,t] <= 
                              eta_el * a * z[dev,n] * solar_irrad * eta_inv,
                              name="Solar_electrical_"+dev+"_"+str(n)+timetag)
                         
                        # Bestrict sold electricity from PV to 70% of the rated power without
                        # battery storage and 50% with storage system
                        if dev == "pv":
                            # linearization constraints
                            M = math.floor(params["A_max"] / a)
                            model.addConstr(lin_pv_bat[n] <= M * sum(x["bat",m] for m in devs["bat"].keys()))
                            model.addConstr(z[dev,n] - lin_pv_bat[n] >= 0)
                            model.addConstr(z[dev,n] - lin_pv_bat[n] <= (1-sum(x["bat",m] for m in devs["bat"].keys()))*M)
        
        dev = "pv"
        if options["EEG"] and options["KfW"]:
            # With EEG and KfW --> cap PV exports at 70 or 50%
            for d in days:
                for t in time_steps:
                    model.addConstr(p_sell[dev,d,t] <= sum(0.7*devs[dev][n]["p_nom"]*(z[dev,n]-lin_pv_bat[n])
                                                     + lin_pv_bat[n]*0.5*devs[dev][n]["p_nom"] for n in devs[dev].keys()),
                                                    name='restrict sold elec')
        elif options["EEG"] and not options["KfW"]:
            # No KfW support --> always cap PV exports at 70%
            for d in days:
                for t in time_steps:
                    model.addConstr(p_sell[dev,d,t] <= sum(0.7*devs[dev][n]["p_nom"]*(z[dev,n]-lin_pv_bat[n])
                                                     + lin_pv_bat[n]*0.7*devs[dev][n]["p_nom"] for n in devs[dev].keys()),
                                                    name='restrict sold elec')
        elif not options["EEG"] and options["KfW"]:
            # No EEG, but with KfW support --> cap PV exports at 50% with BAT, and no limit without BAT
            for d in days:
                for t in time_steps:
                    model.addConstr(p_sell[dev,d,t] <= sum(devs[dev][n]["p_nom"]*(z[dev,n]-lin_pv_bat[n])
                                                     + lin_pv_bat[n]*0.5*devs[dev][n]["p_nom"] for n in devs[dev].keys()),
                                                    name='restrict sold elec')
        else:
            # Neither EEG nor KfW --> do not limit PV exports
            for d in days:
                for t in time_steps:
                    model.addConstr(p_sell[dev,d,t] <= sum(devs[dev][n]["p_nom"]*z[dev,n] for n in devs[dev].keys()),
                                                    name='restrict sold elec')

        # Inverter sizing
        pv_nominal = sum(z["pv",n] * devs["pv"][n]["p_nom"] for n in devs["pv"].keys())
        inv_power = sum(x["inv",n] * devs["inv"][n]["P_nom_DC"] for n in devs["inv"].keys())
        model.addConstr(pv_nominal <= inv_power, name="Sizing_Inverter")

        # Storage equations
        # TES
        dev = "tes"
        for n in devs[dev].keys():
            k_loss = devs[dev][n]["k_loss"]
            max_cap = (params["rho_w"] * params["c_w"] * 
                       devs[dev][n]["volume"] * params["dT_max"] /
                       3600000)
            eta_ch = devs[dev][n]["eta_ch"]
            eta_dch = devs[dev][n]["eta_dch"]
            for d in days:
                for t in time_steps:
                    if t == 0:
                        if np.max(demands["weights"]) == 1:
                            if d == 0:
                               soc_prev = soc_init[dev,n,d]
                            else:
                               soc_prev = soc[dev,n,d-1,params["time_steps"]-1]
                        else:
                            soc_prev = soc_init[dev,n,d]
                    else:
                        soc_prev = soc[dev,n,d,t-1]
                    
                    charge = eta_ch / max_cap * ch[dev,n,d,t]
                    discharge = 1 / eta_dch / max_cap * dch[dev,n,d,t]
                    
                    timetag = "_" + str(d) + "_" + str(t)
                    model.addConstr(soc[dev,n,d,t] == (1-k_loss) * soc_prev + 
                                    dt * (charge - discharge),
                                    name="Storage_bal_"+dev+"_"+str(n)+timetag)
                    
                    model.addConstr(ch[dev,n,d,t] <= 
                                    x[dev,n] * 5 * demands["design_heat_load"],
                                    name="Q_ch_"+str(n)+timetag)
                    
                    model.addConstr(dch[dev,n,d,t] <= 
                                    x[dev,n] * 5 * demands["design_heat_load"],
                                    name="Q_dch_"+str(n)+timetag)

        dev = "bat"
        for n in devs[dev].keys():
            k_loss = devs[dev][n]["k_loss"]
            max_cap = devs[dev][n]["cap"]
            for d in days:
                for t in time_steps:
                    if t == 0:
                        if np.max(demands["weights"]) == 1:
                            if d == 0:
                               soc_prev = soc_init[dev,n,d]
                            else:
                               soc_prev = soc[dev,n,d-1,params["time_steps"]-1]
                        else:
                            soc_prev = soc_init[dev,n,d]
                    else:
                        soc_prev = soc[dev,n,d,t-1]

                    timetag = "_" + str(d) + "_" + str(t)
        
                    model.addConstr(soc[dev,n,d,t] == (1-k_loss) * soc_prev + 
                        dt / max_cap * 
                        (devs[dev][n]["eta"] * ch[dev,n,d,t] - dch[dev,n,d,t]),
                        name="Storage_balance_"+dev+"_"+str(n)+timetag)
        
                    model.addConstr(ch[dev,n,d,t] <= 
                                    x[dev,n] * devs[dev][n]["P_ch_max"],
                                    name="P_ch_max_"+str(n)+timetag)
        
                    model.addConstr(dch[dev,n,d,t] <= 
                                    x[dev,n] * devs[dev][n]["P_dch_max"],
                                    name="P_dch_max_"+str(n)+timetag)
        
        # SOC limits and repetitions
        for dev in storage:
            for n in devs[dev].keys():
                for d in range(params["days"]):
                    # Repetitions
                    if np.max(demands["weights"]) > 1:
                        model.addConstr(soc[dev,n,d,params["time_steps"]-1] ==
                                        soc_init[dev,n,d],
                                        name="repetitions_"+dev+"_"+str(n))

                    # SOC limits
                    model.addConstr(soc_init[dev,n,d] <= x[dev,n],
                          name="soc_init_lim_"+dev+"_"+str(n)+"_"+str(d))
                    for t in range(params["time_steps"]):
                        timetag = "_" + str(d) + "_" + str(t)
                        model.addConstr(soc[dev,n,d,t] <= x[dev,n], 
                                 name="soc_limit_"+dev+"_"+str(n)+timetag)
        
        # Electricity balance
        for d in days:
            for t in time_steps:
                # elec balance for components without hp-tariff (p_use["bat"] referring to discharge)
                model.addConstr(
                        demands["electricity"][d,t]
                        + eh_split["eh_w/o_hp",d,t]
                        + sum(ch["bat",n,d,t] for n in devs["bat"].keys())
                        == p_grid["grid_hou",d,t]
                        + sum(p_use[dev,d,t] for dev in ("pv","bat"))
                        + sum(p_use["chp",n,d,t] for n in devs["chp"].keys()),
                        name="Electricity_balance_w/o_HPtariff"+str(d)+"_"+str(t))
                # elec balance for components with hp-tariff (p_hp["bat"] referring to discharge)
                model.addConstr(
                        sum(power["hp",n,d,t] for n in devs["hp"].keys())
                        + eh_split["eh_w/_hp",d,t]
                        == p_grid["grid_hp",d,t]
                        + sum(p_hp[dev,d,t] for dev in ("pv","bat"))
                        + sum(p_hp["chp",n,d,t] for n in devs["chp"].keys()),
                        name="Electricity_balance_w/_HPtariff"+str(d)+"_"+str(t))
                        
        # Thermal balance
        dev = "tes"
        for d in days:
            for t in time_steps:
                model.addConstr(sum(ch[dev,n,d,t] for n in devs[dev].keys()) 
                    == sum(sum(heat[dv,nb,d,t] 
                       for nb in devs[dv]) 
                       for dv in (heater_ehdhw+solar)),
                       name="Thermal_max_charge_"+str(d)+"_"+str(t))
                model.addConstr(sum(dch[dev,n,d,t] for n in devs[dev].keys())
                    == demands["dhw"][d,t] + demands["sh"][d,t], 
                       name="Thermal_max_discharge_"+str(d)+"_"+str(t))
        
        # High temperature heat generators have to cover at least X% of DHW demand
        dem_dhw = np.sum(demands["dhw"], axis=1)
        for d in days:
            sum_boiler = sum(sum(heat["boiler",n,d,t] for n in devs["boiler"].keys()) for t in time_steps)
            sum_chp = sum(sum(heat["chp",n,d,t] for n in devs["chp"].keys()) for t in time_steps)
            sum_eh_dhw = sum(sum(heat["eh_dhw",n,d,t] for n in devs["eh_dhw"].keys()) for t in time_steps)
            sum_eh_sh = sum(sum(heat["eh",n,d,t] for n in devs["eh"].keys()) for t in time_steps)
            
            sum_stc = sum(x["stc",n] for n in devs["stc"].keys())
            sum_hp = sum(x["hp",n] for n in devs["hp"].keys())
            
            model.addConstr(sum_boiler + sum_chp + sum_eh_dhw + sum_eh_sh >=
                            (1-params["partial_low_temp_dhw"]) * dem_dhw[d] * sum_stc)
            model.addConstr(sum_boiler + sum_chp + sum_eh_dhw + sum_eh_sh >=
                            (1-params["partial_low_temp_dhw"]) * dem_dhw[d] * sum_hp)
        
                
        # Split CHP and PV generation and bat discharge Power into self-consumed, sold and transferred powers
        for d in days:
            for t in time_steps:
                dev = "bat"
                model.addConstr(p_sell[dev,d,t] + p_use[dev,d,t] + p_hp[dev,d,t] == 
                        sum(dch[dev,n,d,t] for n in devs[dev].keys()),
                        name="power=sell+use+hp_"+dev+"_"+str(d)+"_"+str(t))
                dev = "pv"
                model.addConstr(p_sell[dev,d,t] + p_use[dev,d,t] + p_hp[dev,d,t] == 
                        sum(power[dev,n,d,t] for n in devs[dev].keys()),
                        name="power=sell+use+hp_"+dev+"_"+str(d)+"_"+str(t))
                dev = "chp"
                for n in devs[dev].keys():
                    model.addConstr(p_sell[dev,n,d,t] + p_use[dev,n,d,t] + p_hp[dev,n,d,t] == 
                            power[dev,n,d,t],
                            name="power=sell+use+hp_"+dev+"_"+str(d)+"_"+str(t))
                                
        # Split EH power consumption into cases with and without heat pump installed
        dev = "eh"
        # efficiency eta = 1 => P_nom = Q_nom
        p_nom = (max([devs["eh"][n]["Q_nom"] for n in devs["eh"].keys()]) + 
                 max([devs["eh_dhw"][n]["Q_nom"] for n in devs["eh_dhw"].keys()]))

        for d in days:
            for t in time_steps:
                model.addConstr(eh_split["eh_w/o_hp",d,t] + eh_split["eh_w/_hp",d,t]==
                                sum(power["eh",n,d,t] for n in devs["eh"].keys()) + 
                                sum(power["eh_dhw",n,d,t] for n in devs["eh_dhw"].keys()))

                for n in devs[dev].keys():
                    model.addConstr(eh_split["eh_w/_hp",d,t] <=
                                    sum(x["hp",k] for k in devs["hp"].keys()) * p_nom)
                    model.addConstr(eh_split["eh_w/o_hp",d,t] <=
                                    (1-sum(x["hp",k] for k in devs["hp"].keys())) * p_nom)
                            
        # Heat pump's operation depends on storage temperature
        dt_max_sto = params["dT_max"]
        dev = "hp"
        for n in devs[dev].keys():
            for d in days:
                for t in time_steps:
                    timetag = "_"+str(d)+"_"+str(t)
                    sum_soc = sum(soc["tes",n_tes,d,t] 
                                  for n_tes in devs["tes"].keys())
                    model.addConstr(sum_soc <= 1 - y[dev,n,d,t] * 
                                    (1 - devs[dev][n]["dt_max"] / dt_max_sto),
                                    name="Heat_pump_act_"+str(n)+timetag)

        # STC operation is also coupled with storage temperature
        dev = "stc"
#        y_stc[d,t]
        for d in days:
            for t in time_steps:
                timetag = "_"+str(d)+"_"+str(t)
                sum_soc = sum(soc["tes",n_tes,d,t] for n_tes in devs["tes"].keys())
                sum_heat = sum(heat[dev,n,d,t] for n in devs[dev].keys())
                model.addConstr(sum_heat <= y_stc[d,t] * 1000)
                model.addConstr(sum_soc <= 1 - y_stc[d,t] * 
                                    (1 - devs["hp"][1]["dt_max"] / dt_max_sto),
                                    name="stc_act_"+timetag)
                
        # Design heat load
        cap_chp = sum(np.max(devs["chp"][n]["Q_heat"]) * x["chp",n] 
                      for n in devs["chp"].keys())
        cap_boil = sum(np.max(devs["boiler"][n]["Q_heat"]) * x["boiler",n] 
                      for n in devs["boiler"].keys())
        cap_eh = sum(devs["eh"][n]["Q_nom"] * x["eh",n]
                      for n in devs["eh"].keys())
#        cap_eh_dhw = sum(devs["eh_dhw"][n]["Q_nom"] * x["eh_dhw",n]
#                      for n in devs["eh_dhw"].keys())
        cap_hp = sum(devs["hp"][n]["Q_nom"] * x["hp",n] 
                      for n in devs["hp"].keys())
        model.addConstr(cap_chp + cap_boil + cap_eh + cap_hp #+ cap_eh_dhw
                        >= demands["design_heat_load"],
                        name="Design_heat_load")

#%% CO2 emissions
        emissions_gas = sum(sum((G[tar]["boiler",n] + G[tar]["chp",n]) * eco["gas"][tar]["emi"]
                            for n in x_tariff["gas"][tar].keys())
                            for tar in gas_tariffs)
        
        emissions_grid = sum(eco["el"][tar]["emi"]*sum((El[tar]["grid_hou",n] + El[tar]["grid_hp",n])
                            for n in x_tariff["el"][tar].keys()) 
                            for tar in el_tariffs)
               
        emissions_feedin = eco["emi_el_mix"] * sum(demands["weights"][d] * dt * sum(sum(p_sell["chp",n,d,t] 
                                   for n in devs["chp"].keys()) + p_sell["pv",d,t] + p_sell["bat",d,t]
                                   for t in time_steps)
                                   for d in days)
        
        model.addConstr(emission == emissions_gas + emissions_grid - emissions_feedin)
        model.addConstr(0.001 * emission <= max_emi)   

#%% Set start values and branching priority
        if options["load_start_vals"]:
            with open(options["filename_start_vals"], "r") as fin:
                for line in fin:
                    line_split = line.replace("\"", "").split()
                    (model.getVarByName(line_split[0])).Start = float(line_split[1])

        for key in x.keys():
            x[key].BranchPriority = 100
            
#%% Define scenarios:
        model.addConstr(sum(x["tes",n] for n in devs["tes"].keys()) == 1)
        
        for d in days:
            for t in time_steps:
                model.addConstr(p_sell["bat",d,t] == 0)
                model.addConstr(sum(ch["bat",n,d,t] for n in devs["bat"].keys()) <= lim_ch[d,t] * 1000)
                model.addConstr(p_grid["grid_hou",d,t] <= (1-lim_ch[d,t]) * 1000)
                
                model.addConstr(p_sell["pv",d,t] + sum(p_sell["chp",n,d,t] for n in devs["chp"].keys()) <= lim_sell[d,t] * 1000)
                model.addConstr(p_grid["grid_hou",d,t] <= (1-lim_sell[d,t]) * 1000)
                

        if options["scenario"] == "reference":
            model.addConstr(sum(x["chp",n] for n in devs["chp"].keys()) == 0)
            model.addConstr(sum(x["hp",n] for n in devs["hp"].keys()) == 0)
            model.addConstr(sum(x["pv",n] for n in devs["pv"].keys()) == 0)
            model.addConstr(sum(x["stc",n] for n in devs["stc"].keys()) == 0)
            model.addConstr(sum(x["eh",n] for n in devs["eh"].keys()) == 0)
            model.addConstr(sum(x["bat",n] for n in devs["bat"].keys()) == 0)
            model.addConstr(sum(x["boiler",n] for n in devs["boiler"].keys()) == 1)
        elif options["scenario"] == "pv":
            model.addConstr(sum(z["pv",n] for n in devs["pv"].keys()) >= 1)
            model.addConstr(sum(x["bat",n] for n in devs["bat"].keys()) == 0)
        elif options["scenario"] == "kfw":
            model.addConstr(sum(z["pv",n] for n in devs["pv"].keys()) >= 1)
            model.addConstr(sum(x["bat",n] for n in devs["bat"].keys()) == 1)
        elif options["scenario"] == "chp":
            model.addConstr(sum(x["chp",n] for n in devs["chp"].keys()) == 1)
        elif options["scenario"] == "hp":
            model.addConstr(sum(x["hp",n] for n in devs["hp"].keys()) == 1)
        elif options["scenario"] == "free":
            pass

            
        # forced devices
        for dev in forced_devs.keys():
            n = forced_devs[dev]
            model.addConstr(x[dev,n] == 1)
            fdev = str(dev) + str(n)
            

#%%     Set model parameters, execute calculation and get results
        # Set solver parameters
        model.Params.TimeLimit  = params["time_limit"]
        model.Params.MIPGap     = params["mip_gap"]
        model.Params.MIPFocus   = mfocus
        model.Params.NumericFocus = 3
        
        # Execute calculation
        model.optimize()
        
        
        if model.Status in (3,4) or model.SolCount == 0: # "INFEASIBLE" or "INF_OR_UNBD"
            if model.SolCount == 0:
                with open(options["filename_results"]+"unfinished", "a") as fout:
                    fout.write(fdev+"\n")
                
            return (gp.GRB.INFINITY, 0)
        else:
            
            # purchase
            res_x = {}
            for dev in devs.keys():
                res_x[dev] = {n : x[dev,n].X for n in devs[dev]}

            # Operation
            res_y = {}
            for dev in heater_ehdhw:
                res_y[dev] = {}
                for n in devs[dev].keys():
                    res_y[dev][n] = np.array([[y[dev,n,d,t].X for t in time_steps] for d in days])

            # Solar modules
            res_z = {}
            for dev in solar:
                res_z[dev] = {n : z[dev,n].X for n in devs[dev]}
                
            # Tariffs
            res_x_tariff = {"el":{}, "gas":{}}
            for key in ("el", "gas"):
                for tar in x_tariff[key].keys():
                    res_x_tariff[key][tar] = {}
                    for dev in x_tariff[key][tar].keys():
                        res_x_tariff[key][tar][dev] = x_tariff[key][tar][dev].X

            res_x_gas    = {tar: x_gas[tar].X for tar in x_gas.keys()}
            res_x_el     = {tar: x_el[tar].X  for tar in x_el.keys()}

            # heat and electricity output
            res_power = {}
            res_heat  = {}
            for dev in (heater_ehdhw+solar):
                res_power[dev] = {}
                res_heat[dev]  = {}
                for n in devs[dev].keys():
                    res_power[dev][n] = np.array([[power[dev,n,d,t].X for t in time_steps] for d in days])
                    res_heat[dev][n]  = np.array([[heat[dev,n,d,t].X  for t in time_steps] for d in days])

            res_energy = {}
            for dev in heater_ehdhw:
                res_energy[dev] = {}
                for n in devs[dev].keys():
                    res_energy[dev][n] = np.array([[energy[dev,n,d,t].X for t in time_steps] for d in days])
            
            # Gas/El consumption        
            res_G_total  = {dev: G_total[dev].X  for dev in G_total.keys()}
            res_El_total = {dev: El_total[dev].X for dev in El_total.keys()}
            res_G  = {}
            res_El = {}
            for tar in gas_tariffs:
                res_G[tar] = {}
                for dev in ("boiler","chp"):
                    for n in x_tariff["gas"][tar].keys():
                        res_G[tar][dev,n] = G[tar][dev,n].X
            for tar in el_tariffs:
                res_El[tar] = {}
                for dev in ("grid_hou","grid_hp"):
                    for n in x_tariff["el"][tar].keys():
                        res_El[tar][dev,n] = El[tar][dev,n].X

            # State of charge for storage systems
            res_soc = {}
            for dev in storage:
                res_soc[dev] = {}
                for n in devs[dev].keys():
                    res_soc[dev][n] = np.array([[soc[dev,n,d,t].X for t in time_steps] for d in days])
        
            # Purchased power from the grid for either feeding a hp tariff component or a different (standard/eco tariff)
            res_p_grid          = {}
            res_p_grid["house"] = np.array([[p_grid["grid_hou",d,t].X for t in time_steps] for d in days])
            res_p_grid["hp"]    = np.array([[p_grid["grid_hp",d,t].X  for t in time_steps] for d in days])

            # Charge and discharge power for storage
            res_ch  = {}
            res_dch = {}
            for dev in storage:
                res_ch[dev]  = {}
                res_dch[dev] = {}
                for n in devs[dev].keys():
                    res_ch[dev][n]  = np.array([[ch[dev,n,d,t].X  for t in time_steps] for d in days])
                    res_dch[dev][n] = np.array([[dch[dev,n,d,t].X for t in time_steps] for d in days])

            # Power going from an electricity offering component to the demand/the grid/a hp tariff component
            res_p_use  = {}
            res_p_sell = {}
            res_p_hp   = {}
            for dev in ("pv", "bat"):
                res_p_use[dev]  = np.array([[p_use[dev,d,t].X  for t in time_steps] for d in days])
                res_p_sell[dev] = np.array([[p_sell[dev,d,t].X for t in time_steps] for d in days])
                res_p_hp[dev]   = np.array([[p_hp[dev,d,t].X   for t in time_steps] for d in days])
            dev = "chp"
            res_p_use[dev]  = np.array([[[p_use[dev,n,d,t].X  for t in time_steps] for d in days]
                                                              for n in devs[dev].keys()])
            res_p_sell[dev] = np.array([[[p_sell[dev,n,d,t].X for t in time_steps] for d in days]
                                                              for n in devs[dev].keys()])
            res_p_hp[dev]   = np.array([[[p_hp[dev,n,d,t].X   for t in time_steps] for d in days]
                                                              for n in devs[dev].keys()])
            
            # Costs
            res_c_inv   = {dev: c_inv[dev].X    for dev in c_inv.keys()}
            res_c_om    = {dev: c_om[dev].X     for dev in c_om.keys()}
            res_c_dem   = {dev: c_dem[dev].X    for dev in c_dem.keys()}
            res_c_fix   = {dev: c_fix[dev].X    for dev in c_fix.keys()}
            res_c_eeg   = {dev: c_eeg[dev].X    for dev in c_eeg.keys()}
            res_rev     = {dev: revenue[dev].X  for dev in revenue.keys()}
            res_sub     = {dev: subsidy[dev].X  for dev in subsidy.keys()}
            
            res_c_total = (sum(res_c_inv[key]  for key in c_inv.keys())
                           + sum(res_c_om[key]    for key in c_om.keys())
                           + sum(res_c_dem[key]   for key in c_dem.keys())
                           + sum(res_c_fix[key]   for key in c_fix.keys())
                           + sum(res_c_eeg[key]   for key in c_eeg.keys())
                           - sum(res_rev[key] for key in revenue.keys())
                           - sum(res_sub[key] for key in subsidy.keys()))  
            
            res_soc_init = {}
            for dev in storage:
                res_soc_init[dev] = {}
                for n in devs[dev].keys():
                    res_soc_init[dev][n] = np.array([soc_init[dev,n,d].X for d in days])

            # Emissions 
            res_emission_max = max_emi
            res_emission = emission.X / 1000

            if options["store_start_vals"]:
                with open(options["filename_start_vals"], "w") as fout:
                    for var in model.getVars():
                        if var.VType == "B":
                            fout.write(var.VarName + "\t" + str(int(var.X)) + "\n")

            # Save results 
            with open(options["filename_results"], "wb") as fout:
                pickle.dump(res_x, fout, pickle.HIGHEST_PROTOCOL)
                pickle.dump(res_y, fout, pickle.HIGHEST_PROTOCOL)
                pickle.dump(res_z, fout, pickle.HIGHEST_PROTOCOL)
                pickle.dump(res_x_tariff, fout, pickle.HIGHEST_PROTOCOL)
                pickle.dump(res_x_gas, fout, pickle.HIGHEST_PROTOCOL)
                pickle.dump(res_x_el, fout, pickle.HIGHEST_PROTOCOL)
                pickle.dump(res_power, fout, pickle.HIGHEST_PROTOCOL)
                pickle.dump(res_heat, fout, pickle.HIGHEST_PROTOCOL)
                pickle.dump(res_energy, fout, pickle.HIGHEST_PROTOCOL)
                pickle.dump(res_p_grid, fout, pickle.HIGHEST_PROTOCOL)
                pickle.dump(res_G, fout, pickle.HIGHEST_PROTOCOL)
                pickle.dump(res_G_total, fout, pickle.HIGHEST_PROTOCOL)
                pickle.dump(res_El, fout, pickle.HIGHEST_PROTOCOL)
                pickle.dump(res_El_total, fout, pickle.HIGHEST_PROTOCOL)
                pickle.dump(res_soc, fout, pickle.HIGHEST_PROTOCOL)
                pickle.dump(res_soc_init, fout, pickle.HIGHEST_PROTOCOL)
                pickle.dump(res_ch, fout, pickle.HIGHEST_PROTOCOL)
                pickle.dump(res_dch, fout, pickle.HIGHEST_PROTOCOL)
                pickle.dump(res_p_use, fout, pickle.HIGHEST_PROTOCOL)
                pickle.dump(res_p_sell, fout, pickle.HIGHEST_PROTOCOL)
                pickle.dump(res_p_hp, fout, pickle.HIGHEST_PROTOCOL)
                pickle.dump(res_c_inv, fout, pickle.HIGHEST_PROTOCOL)
                pickle.dump(res_c_om, fout, pickle.HIGHEST_PROTOCOL)
                pickle.dump(res_c_dem, fout, pickle.HIGHEST_PROTOCOL)
                pickle.dump(res_c_fix, fout, pickle.HIGHEST_PROTOCOL)
                pickle.dump(res_c_eeg, fout, pickle.HIGHEST_PROTOCOL)
                pickle.dump(res_c_total, fout, pickle.HIGHEST_PROTOCOL)
                pickle.dump(res_rev, fout, pickle.HIGHEST_PROTOCOL)
                pickle.dump(res_sub, fout, pickle.HIGHEST_PROTOCOL)
                pickle.dump(res_emission, fout, pickle.HIGHEST_PROTOCOL)  
                pickle.dump(res_emission_max, fout, pickle.HIGHEST_PROTOCOL)
                pickle.dump(model.ObjVal, fout, pickle.HIGHEST_PROTOCOL)
                pickle.dump(model.Runtime, fout, pickle.HIGHEST_PROTOCOL)
                pickle.dump(model.MIPGap, fout, pickle.HIGHEST_PROTOCOL)
                
            # Return results            
            return(res_c_total, res_emission)

    except gp.GurobiError as e:
        print("")        
        print("Error: "+e.message)
