#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 20 18:42:33 2015

@author: Thomas
"""
from __future__ import division

import xlrd
import numpy as np


interp = lambda x1, x2, y1, y2, x_v: y1 + (x_v-x1) / (x2-x1) * (y2-y1)


def _parse_boiler(str):
    line = str.replace("[","").replace("]","").replace(" ", "")
    line_split = line.split(";")
    
    return np.array([float(val) for val in line_split])


def _parse_hp(str):
    data = {}

    str_split1 = str.replace("'","").split(":{")
    
    t_flow = [int(str_split1[1])]
    data_other = []
    for i in range(2,len(str_split1)-1):
        temp = str_split1[i].split("},")
        t_flow.append(int(temp[1]))
        data_other.append(temp[0])
    data_other.append(str_split1[-1].replace("}}",""))
    
    for i in range(len(t_flow)):
        data[t_flow[i]] = {}
        
        # Handle data_other
        temp1 = data_other[i].split("],")
        for j in range(len(temp1)):
            temp2 = temp1[j].replace("[","").replace("]","")
            
            temp_split = temp2.split(":")
            
            data[t_flow[i]][int(temp_split[0])] = [float(temp3) for temp3 in temp_split[1].split(",")]
    
    return data


def _interp_flow(data, t_flow):
    sorted_keys = np.sort(data.keys())
    
    if sorted_keys[0] >= t_flow:
        result = data[sorted_keys[0]]
    elif sorted_keys[-1] <= t_flow:
        result = data[sorted_keys[-1]]
    else:
        i = 1
        while i < len(sorted_keys) and sorted_keys[i] <= t_flow:
            i += 1
        
        t_ambs = data[sorted_keys[i]].keys()
    
        result = {}
    
        for t_amb in t_ambs:
            result[t_amb] = interp(sorted_keys[i-1], 
                                   sorted_keys[i], 
                                   np.array(data[sorted_keys[i-1]][t_amb]),
                                   np.array(data[sorted_keys[i]][t_amb]),
                                   t_flow)

    return result
    
def _interp_amb(data, t_amb):
    result = {}
    
    sorted_keys = np.sort(data.keys())
    for d in range(t_amb.shape[0]):
        for t in range(t_amb.shape[1]):
            if sorted_keys[0] >= t_amb[d,t]:
                result[d,t] = data[sorted_keys[0]]
            elif sorted_keys[-1] <= t_amb[d,t]:
                result[d,t] = data[sorted_keys[-1]]
            else: 
                i = 1
                while i < len(sorted_keys) and sorted_keys[i] <= t_amb[d,t]:
                    i += 1
                
                result[d,t] = interp(sorted_keys[i-1],
                                     sorted_keys[i], 
                                     np.array(data[sorted_keys[i-1]]),
                                     np.array(data[sorted_keys[i]]),
                                     t_amb[d,t])

    return result    

def _interp_hp(p_hp, q_hp, t_flow, temperature_ambient):
    # Interpolation flow temperature
    p_hp_int = _interp_flow(p_hp, t_flow)
    q_hp_int = _interp_flow(q_hp, t_flow)
    
    
    # Interpolation ambient temperature
    res_p_hp = _interp_amb(p_hp_int, temperature_ambient)
    res_q_hp = _interp_amb(q_hp_int, temperature_ambient)
    
    return (res_p_hp, res_q_hp)


def _parse_tariffs(fixed, variable):
    """
    """

    def _parse_tariff(tariff):
        """
        """
        lb  = []
        ub  = []
        tar = []
        
        tariff = tariff.replace(" ", "")
        tariff = tariff.replace("[", "")
        tariff = tariff.replace("]", "")
        split  = tariff.split(";")
        
        for val in split:
            (value_range, value) = val.split(":")
            (temp_lb, temp_ub)   = value_range.split("-")
            lb.append(float(temp_lb.replace(",",".")))
            ub.append(float(temp_ub.replace(",",".")))
            tar.append(float(value.replace(",",".")))
        
        return (lb, ub, tar)
    
    (lb, ub, fix) = _parse_tariff(fixed)
    (lb, ub, var) = _parse_tariff(variable)
    
    return (lb, ub, fix, var)

def read_economics(devices, filename="further_parameters.xlsx"):
    """
    Read in economic parameters and update residual values of devices.
    
    Parameters
    ----------
    devices : dictionary
        All device specific characteristics.
    filename : string, optional
        Excel-file with the economic and other parameters.
    
    Returns
    -------
    eco : dictionary
        Information on economic parameters.
    par : dictionary
        All non-economic and non-technical parameters.
    devices : dictionary
        All device specific characteristics.
    """
    book = xlrd.open_workbook(filename)
    
    sheet_eco = book.sheet_by_name("gen_economics")
    sheet_par = book.sheet_by_name("further_parameters")
    sheet_gas = book.sheet_by_name("gas_economics")
    sheet_el  = book.sheet_by_name("el_economics")
    
    eco = {}
    par = {}
    
    # Economics
    t_calc = sheet_eco.cell_value(1,1)
    eco["t_calc"] = t_calc
    eco["tax"]    = sheet_eco.cell_value(2,1)
    eco["rate"]   = sheet_eco.cell_value(3,1)
    eco["q"]      = 1 + eco["rate"]
    eco["crf"]    = ((eco["q"] ** eco["t_calc"] * eco["rate"]) / 
                     (eco["q"] ** eco["t_calc"] - 1))
    
    # Current price change rates
    eco["prChange"] = {}
    eco["prChange"]["el"]   = sheet_eco.cell_value(6,1)
    eco["prChange"]["gas"]  = sheet_eco.cell_value(7,1)
    eco["prChange"]["eex"]  = sheet_eco.cell_value(8,1)
    eco["prChange"]["infl"] = sheet_eco.cell_value(9,1)

    pC = eco["prChange"]
    eco["b"] = {key: ((1 - (pC[key] / eco["q"]) ** eco["t_calc"]) / 
                      (eco["q"] - pC[key]))
                for key in pC.keys()}
    
    # Prices and tariff amount gradations
    eco["el"]                 = {}
    eco["gas"]                = {}
    eco["sub_chp"]            = {}
    
    # Subsidies
    eco["sub_chp"]["self"]    = sheet_eco.cell_value(12,1)
    eco["sub_chp"]["sell"]    = sheet_eco.cell_value(13,1)
    eco["sub_chp"]["lump"]    = sheet_eco.cell_value(14,1)
    eco["sub_chp"]["t"]       = sheet_eco.cell_value(15,1)
    
    eco["sub_bat_max"]        = sheet_eco.cell_value(16,1)
    eco["sub_bat"]            = sheet_eco.cell_value(17,1)
    # Feed-ins
    eco["sell", "chp"]        = sheet_eco.cell_value(18,1)
    eco["sell", "pv", "10"]   = sheet_eco.cell_value(19,1)
    eco["sell", "pv", "40"]   = sheet_eco.cell_value(21,1)
    eco["sell", "pv", "100"]  = sheet_eco.cell_value(22,1)
    eco["eeg_levy"]           = sheet_eco.cell_value(20,1)
    
    eco["gas"]["energy_tax"]  = sheet_gas.cell_value(1,1)  # in €/kWh
    
    # Read gas prices
    n = 2
    while (not sheet_gas.cell_value(n,0) == "Name") and n < sheet_gas.nrows:
        n += 1
    for i in range(n+1, sheet_gas.nrows):
        (lb, ub, fix, var) = _parse_tariffs(sheet_gas.cell_value(i,1),
                                            sheet_gas.cell_value(i,2))
        eco["gas"][sheet_gas.cell_value(i,0)] = {"lb": lb, "ub": ub, 
                                                 "fix": fix, "var":var}
        eco["gas"][sheet_gas.cell_value(i,0)]["emi"] = float(sheet_gas.cell_value(i,3))
    # Read el prices
    n = 2
    while (not sheet_el.cell_value(n,0) == "Name") and n < sheet_el.nrows:
        n += 1
    for i in range(n+1, sheet_el.nrows):
        (lb, ub, fix, var) = _parse_tariffs(sheet_el.cell_value(i,1),
                                            sheet_el.cell_value(i,2))
        eco["el"][sheet_el.cell_value(i,0)] = {"lb": lb, "ub": ub, 
                                               "fix": fix, "var":var}
        eco["el"][sheet_el.cell_value(i,0)]["emi"] = float(sheet_el.cell_value(i,3))
        eco["el"][sheet_el.cell_value(i,0)]["hp"] = int(sheet_el.cell_value(i,4))
        
    # Determine residual values
    for dev in devices.keys():
        for number in devices[dev].keys():
            t_life = devices[dev][number]["T_op"]
            rval = (t_life - t_calc) / t_life / (eco["q"] ** t_calc)
            devices[dev][number]["rval"] = rval
    
    # Further parameters
    par["A_max"]      = sheet_par.cell_value(2,1)
    par["mip_gap"]    = sheet_par.cell_value(3,1)
    par["time_limit"] = sheet_par.cell_value(4,1)
    par["rho_w"]      = sheet_par.cell_value(7,1)
    par["c_w"]        = sheet_par.cell_value(8,1)
    par["dT_max"]     = sheet_par.cell_value(9,1)
    
    eco["emi_el_mix"] = sheet_par.cell_value(10,1)
    par["partial_low_temp_dhw"] = sheet_par.cell_value(11,1)

    return (eco, par, devices)
    
def compute_parameters(par, number_clusters, len_day):
    """
    Add number of days, time steps per day and temporal discretization to par.
    
    Parameters
    ----------
    par : dictionary
        Dictionary which holds non-device-characteristic and non-economic 
        parameters.
    number_clusters : integer
        Number of allowed clusters.
    len_day : integer
        Time steps per day
    """
    par["days"] = number_clusters
    par["time_steps"] = len_day
    par["dt"] = 24 / len_day
    
    return par
    
    
def read_devices(timesteps, days, 
                 temperature_ambient, temperature_flow, temperature_design,
                 solar_irradiation,
                 days_per_cluster,
                 filename="devices.xlsx"):
    """
    Read all devices from a given file.
    
    Parameters
    ----------
    timesteps : integer
        Number of time steps per typical day
    days : integer
        Number of typical days
    temperature_ambient : array_like
        2-dimensional array [days, timesteps] with the ambient temperature in 
        degree Celsius
    temperature_flow : float or array_like
        Required flow temperature in degree Celsius. Either as float value or
        as 2-dimensional array [days, timesteps]
    temperature_design : float
        Nominal design temperature in degree Celsius (-12 for Aachen, Germany)
    solar_irradiation : array_like
        Solar irradiation in Watt per square meter on the tilted areas on 
        which STC or PV will be installed.
    filename : string, optional
        Path to the *.xlsx file containing all available devices
    
    Return
    ------
    results : dictionary
        Dictionary containing the information for each device specified in 
        the given input file.
    """
    # Initialize results
    results = {}
    
    # Open work book
    book = xlrd.open_workbook(filename)
    
    # Get all available sheets
    available_sheets = book.sheet_names()
    
    # Iterate over all sheets
    for dev in available_sheets:
        # Read each sheet
        results[dev] = _read_sheet(book.sheet_by_name(dev), dev, 
                                   timesteps, days, 
                                   temperature_ambient, temperature_flow,
                                   temperature_design,
                                   solar_irradiation,
                                   days_per_cluster)
    
    return results

def _read_sheet(sheet, device, timesteps, days, 
                temperature_ambient, temperature_flow, temperature_design,
                solar_irradiation, days_per_cluster):
    """
    Parameters
    ----------
    sheet : sheet-object
        Sheet of the workbook containing all available devices
    device : string
        - `"boiler"`    : Boiler
        - `"chp"`       : CHP unit
        - `"hp"`        : Heat pump
        - `"eh"`        : Electrical heater
        - `"pv"`        : Photovoltaic modules
        - `"stc"`       : Solar thermal collectors
        - `"tes"`       : Thermal energy storage units
        - `"bat"`       : Battery units
    timesteps : integer
        Number of time steps per typical day
    days : integer
        Number of typical days
    temperature_ambient : array_like
        2-dimensional array [days, timesteps] with the ambient temperature in 
        degree Celsius
    temperature_flow : float or array_like
        Required flow temperature in degree Celsius. Either as float value or
        as 2-dimensional array [days, timesteps]
    temperature_design : float
        Nominal design temperature in degree Celsius (-12 for Aachen, Germany)
    solar_irradiation : array_like
        Solar irradiation in Watt per square meter on the tilted areas on 
        which STC or PV will be installed.
        
    Implemented characteristics
    ---------------------------
    - eta = Q/P
    - omega = (Q+P) / E
    """
    # Initialize results
    results = {}
    
    # Define infinity
    infinity = np.inf
    

    # Read all rows but the headers:
    for row in range(1, sheet.nrows):
        # Create new dictionary for current entry. Add common inputs.
        current_results = {}
        # Handle each device separately
        if device in ("eh", "eh_dhw"):
            current_results["Q_nom"]   = sheet.cell_value(row, 1)
            current_results["mod_lvl"] = sheet.cell_value(row, 2)
            current_results["c_inv"]   = sheet.cell_value(row, 3)
            current_results["c_om"]    = sheet.cell_value(row, 4)
            current_results["T_op"]    = sheet.cell_value(row, 5)
            
            # Electrical heater has a given electrical efficiency
            eta_el = sheet.cell_value(row, 6)
            
            # Electrical heaters do not require any gas energy (E), therefore
            # omega has to be large
            current_results["omega"] = np.ones((days, timesteps)) * infinity
            
            # The effiency is defined like eta (Q/P)
            current_results["eta"]   = np.ones((days, timesteps)) * eta_el
            
        elif device == "boiler":
            #each boiler has multiple performance points that are saved as an array
            current_results["Q_heat"]   = _parse_boiler(sheet.cell_value(row, 1))
            current_results["Q_gas"]    = _parse_boiler(sheet.cell_value(row, 2))
            current_results["c_inv"]    = sheet.cell_value(row, 3)
            current_results["c_om"]     = sheet.cell_value(row, 4)
            current_results["T_op"]     = sheet.cell_value(row, 5)
            #the elerctricity is set to zero for all performance points
            current_results["P_el"]     = np.zeros((len(current_results["Q_heat"])))
            current_results["Mod"]      = np.str("stageless")
            
            current_results["q_nom"]    = np.max(current_results["Q_heat"])
    

        elif device == "chp":
            #each chp has multiple performance points that are saved as an array
            current_results["Q_heat"] = _parse_boiler(sheet.cell_value(row, 1))
            current_results["P_el"]   = _parse_boiler(sheet.cell_value(row, 2))
            current_results["Q_gas"]  = _parse_boiler(sheet.cell_value(row, 3))
            current_results["c_inv"]  = sheet.cell_value(row, 5)
            current_results["c_om"]   = sheet.cell_value(row, 6)
            current_results["T_op"]   = sheet.cell_value(row, 7)
            current_results["Mod"]    = sheet.cell_value(row, 4)
          
            current_results["p_nom"]   = np.max(current_results["P_el"])
            current_results["q_nom"]   = np.max(current_results["Q_heat"])

            
        elif device == "hp":
            """
            create set of operation points given by manufacturer that enable 
            interpolation between ambient temperature and 2 flow temperatures
            
            data_nom/min_pi: data set that consits of ambient temperatures with
            corresponding heat,COP and power values
            for a certain flow temperature and a modulation level(nom or min)
         
            """
            current_results["Mod"]     = sheet.cell_value(row, 3)
            current_results["c_inv"]   = sheet.cell_value(row, 4)
            current_results["c_om"]    = sheet.cell_value(row, 5)
            current_results["T_op"]    = sheet.cell_value(row, 6)
            current_results["dt_max"]  = sheet.cell_value(row, 7)
            current_results["Q_gas"]   = np.zeros_like(temperature_ambient)
            
            Q_HP = _parse_hp(sheet.cell_value(row,1))
            P_HP = _parse_hp(sheet.cell_value(row,2))
            
            p_hp_t, q_hp_t = _interp_hp(P_HP, Q_HP, temperature_flow, temperature_ambient)

            current_results["P_el"] = p_hp_t
            current_results["Q_heat"] = q_hp_t

            # Compute nominal heat output at the design temperature and HT flow
            q_nom_int = _interp_flow(Q_HP, temperature_flow)
            current_results["Q_nom"] = max(_interp_amb(q_nom_int, np.array([[temperature_design]]))[0,0])

            
        elif device == "stc":
            current_results["c_inv"] = sheet.cell_value(row, 1)
            current_results["c_om"]  = sheet.cell_value(row, 2)
            current_results["T_op"]  = sheet.cell_value(row, 3)
            current_results["area"]  = sheet.cell_value(row, 4)

            zero_loss    = sheet.cell_value(row, 5) # Optical efficiency
            first_order  = sheet.cell_value(row, 6) # Linear temperature losses
            second_order = sheet.cell_value(row, 7) # Quadratic temp. losses
            
            current_results["eta_el"] = np.zeros((days, timesteps))
            temp_diff = temperature_flow - temperature_ambient
            # Compute the thermal efficiency based on the optical, linear and
            # quadratic temperature loss coefficients
            eta_th = np.zeros_like(solar_irradiation)

            eta_th[solar_irradiation>0] = (zero_loss
                     - first_order * temp_diff[solar_irradiation>0] / solar_irradiation[solar_irradiation>0]
                     - second_order * temp_diff[solar_irradiation>0]**2 / solar_irradiation[solar_irradiation>0])
    
            eta_th[solar_irradiation <= 0.00001] = 0
            current_results["eta_th"] = np.maximum(eta_th, 0)
            
        elif device == "pv":
            current_results["c_inv"] = sheet.cell_value(row, 2)
            current_results["c_om"]  = sheet.cell_value(row, 3)
            current_results["T_op"]  = sheet.cell_value(row, 4)
            current_results["area"]  = sheet.cell_value(row, 5)

            p_NOCT = sheet.cell_value(row, 1)
            i_NOCT = 0.8 # kW / m2
            t_NOCT = sheet.cell_value(row, 6)
            gamma  = sheet.cell_value(row, 7)
            current_results["p_nom"] = sheet.cell_value(row,8)
            
            # Interpolate cell temperature.
            # Without solar irradiation, the cell temperature has to be equal
            # to the ambient temperature. At NOCT irradiation, the cell's 
            # temperature has to be equal to t_NOCT
            t_cell = (temperature_ambient + solar_irradiation / i_NOCT * 
                                            (t_NOCT - temperature_ambient))
            eta_NOCT = p_NOCT / (current_results["area"] * i_NOCT)
            # Compute electrical efficiency of the cell
            eta_el = eta_NOCT * (1 + gamma / 100 * (t_cell - t_NOCT))
            
            current_results["eta_th"] = np.zeros((days, timesteps))
            current_results["eta_el"] = eta_el
        
        elif device == "tes":
            current_results["c_inv"]    = sheet.cell_value(row, 1)
            current_results["c_om"]     = sheet.cell_value(row, 2)
            current_results["T_op"]     = sheet.cell_value(row, 3)
            current_results["eta_ch"]   = sheet.cell_value(row, 6)
            current_results["eta_dch"]  = sheet.cell_value(row, 7)
            
            standby_losses = sheet.cell_value(row, 4) # in kWh / 24h
            volume = sheet.cell_value(row, 5)         # in m3
            
            temp_diff_norm = 45 # K
            heat_cap       = 4180 # J/(kgK)
            density        = 1000 # kg/m3
            energy_content = volume * temp_diff_norm * heat_cap * density # J
            
            k_loss_day = 1 - standby_losses / energy_content * 3600*1000
            current_results["k_loss"] = 1 - (k_loss_day ** (1 / timesteps))
            current_results["volume"] = volume
        
        elif device == "bat":
            current_results["c_inv"]     = sheet.cell_value(row, 1)
            current_results["c_om"]      = sheet.cell_value(row, 2)
            current_results["T_op"]      = sheet.cell_value(row, 3)
            current_results["cap"]       = sheet.cell_value(row, 4)
            current_results["eta"]       = sheet.cell_value(row, 5)
            current_results["P_ch_max"]  = sheet.cell_value(row, 6)
            current_results["P_dch_max"] = sheet.cell_value(row, 7)
            
            current_results["k_loss"]    = 0
        
        elif device == "inv":
            current_results["P_nom_DC"]  = sheet.cell_value(row, 1)
            current_results["c_inv"]     = sheet.cell_value(row, 2)
            current_results["c_om"]      = sheet.cell_value(row, 3)
            current_results["T_op"]      = sheet.cell_value(row, 4)
            current_results["eta"]       = sheet.cell_value(row, 5)
            
        results[row] = current_results
        
        
        
    return (results)

if __name__ == "__main__":
    # Short example
    timesteps = 24
    days = 5
    # Random temperatures between -10 and +20 degC:
    temperature_ambient = np.random.rand(days, timesteps) * 30 - 10
    
    temperature_design = -12 # Aachen
    
    solar_irradiation = np.random.rand(days, timesteps)
        
    (devs) = read_devices(timesteps, days, temperature_ambient, 
                        temperature_flow=55,
                        temperature_design=temperature_design,
                        solar_irradiation=solar_irradiation,
                        days_per_cluster=[1]*days)

    (eco, par, devs) = read_economics(devs)
    
   