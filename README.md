# BESopt #

BESopt (Building Energy System design, sizing and operation OPTimization) is an optimization model for the design and operation of building energy systems.

Additional information on this repository can be found in the corresponding publication: <http://www.sciencedirect.com/science/article/pii/S0360544217314007>


## Required software ##

BESopt is written entirely in **Python 2.7** but should also be compatible with **Python 3.4**.

The following Python packages are currently in use:
- numpy
- xlrd
- math

Further, the optimization requires the [Gurobi optimizer](http://www.gurobi.com/), which is freely available for academic applications.

In order to run this tool, just parameterize your devices and parameters lists (devices.xlsx and further_parameters.xlsx), provide hourly inputs in raw_inputs. Finally, adjust and execute run.py.