# Beo: Model heat flow and (U-Th)/He thermochronology in a hydrothermal system




*Elco Luijendijk, McGill University & Göttingen University*


# Introduction

Beo is a model of heat flow in hot springs and hydrothermal systems. The model code uses the generic finite element code escript (https://launchpad.net/escript-finley) to solve the advective and conductive heat flow in a 2D cross-section of the subsurface. The resulting temperature history is used to calculate two low-temperature thermochronometers, apatite fission track data and apatite (U-Th)/He data. The modeled values of these thermochronometers can then be compared to measured values. Beo also support automated model runs to explore which parameter values like fluid fluxes, fault geometry and age of the hydrothermal system best match the thermochronometer data, as well as present-day temperature data. 

 

# Downloading Beo

* click the download link on the right for a zip file of the source code
* or clone the repository

# Installing & running 

* Install Escript

    - get the code here: https://launchpad.net/escript-finley
    - an installation guide can be found here: http://esys.geocomp.uq.edu.au/docs

* Unzip the beo source code  
* Navigate to the directory where you have installed escript. Go to the subdirectory bin (``python-escript/bin/``) and run Beo by executing the following command from the command line:
	
````bash
./run-escript beo_dir/beo.py
````	

Where ``beo_dir`` is the directory where you have saved Beo.

Alternatively use the command 

````bash
./run-escript -e
````

This will show you three lines that define environment variables that your system needs to be able to find the location of escript. Add these lines to your .bashrc (Ubuntu linux) or profile file in your home directory. After adding these lines and logging out and in again, you can start beo by going to the folder where the beo code is located (so not to the escript/bin directory) and start beo.py like any regular python code:

````bash
python beo.py
````


# Required modules

* numpy:  http://sourceforge.net/projects/numpy/files/NumPy/
* For creating figures (optional):

    - scipy: http://sourceforge.net/projects/scipy/files/scipy/
    - matplotlib: http://matplotlib.org/downloads.html

These modules are available as standalone packages. For mac and windows an easy way to get a working version of python and these modules is to install the Enthought Canopy (https://www.enthought.com/products/canopy) or pythonxy (https://code.google.com/p/pythonxy/wiki/Welcome).

Beo was tested on Ubuntu 14.04 and 16.04 


# Model input & output

## Input:

All model input is contained in the file ``model_parameters.py`` located in the directory model_parameters. See the comments in this file for a brief explanation of what each parameter does.


##Multiple model runs

Optionally you can start automated runs to test a range of parameter combinations. This is useful for automated sensitivity or uncertainty analysis. 

You can start a series of model runs using the script ``run_multiple_models.py``. Once you start this script, any parameter value that is included in the file ``model_parameters/parameter_ranges.py`` will be used as input for a single model run. All results will be stored an written to a comma separated (.csv) file names ``model_output/model_params_and_results_x_runs.csv``. 

You can include any model parameter in the automated runs. Simply copy a parameter from the ``model_parameters.py`` script into the ``parameter_ranges.py`` script and add _s to the parameter name. For instance to test multiple values of the thermal gradient, add `thermal_gradient_s = [0.03, 0.04, 0.05]` to test the effect of geothermal gradients of 0.03, 0.04 and 0.05 C/m, respectively.

There are two options for running multiple model runs. The default is a sensitivity run. In each model run a single parameter will be changed, while all other parameters are kept constant at the default value specified in the ``model_parameters.py`` file. Alternatively you can test all parameter combinations by changing `parameter_combinations = False` to `parameter_combinations = True`. Note that this will generate a lot of model runs, testing ten parameter values for two parameters each will generate 10*10 = 100 model runs, for three parameters this increase to a 1000 model runs, etc...


## Output

* After each model run, the modeled temperature field and (U-Th)/He data are stored in the directory ``model_output`` as a .pck file, which can be read using the Python's pickle module. 
* In addition, Beo saves a comma separated file containing the model parameters and a summary of the results for each model run and each timestep in the same directory.


##Making figures

Currently there are two Python scripts that will generate figures of the model output:

* The script ``make_figures.py`` will make a single figure of the final temperature field for output files (with extension .pck) found in the directory ``model_output``. After running this script you will be prompted to select the output file that you would lik a figure of.
* ``make_figure_2models.py`` will make a figure with two panels containing the modeled temperatures and (U-Th)/He data for two model runs. THe model runs are specified in the script itself, for example like this: ``files = ['model_output/T_field_model_run_14_(-3000.0, 0.03).pck', 'model_output/T_field_model_run_32_(-6000.0, 0.03).pck']``


