# BESopt #

BESopt (Building Energy System design, sizing and operation OPTimization) is an optimization model for the design and operation of building energy systems.

Additional information on this repository can be found in the corresponding publication: <http://www.sciencedirect.com/science/article/pii/S0360544217314007>


## Required software ##

BESopt is written entirely in **Python 2.7** but should also be compatible with **Python 3.X**.

The following Python packages are currently in use:
- numpy
- xlrd
- math

Further, the optimization requires the [Gurobi optimizer](http://www.gurobi.com/), which is freely available for academic applications.

## How to run these scripts ##

In order to run this tool, parameterize your devices and parameters tables (devices.xlsx and further_parameters.xlsx), provide hourly inputs in raw_inputs\your_building\ . Finally, adjust and execute run.py. The adjustments mainly concern adding your_building to the houses over which the script is run (block if __name__ == "__main__")

The run.py script loops over all listed buildings and executes the optimization model building_optimization_decomp.py multiple times to create a multi-objective optimization.
The results are stored in the results folder. Furthermore, each of these multi-objective runs creates a new set of start values that are stored in the start_values folder.
