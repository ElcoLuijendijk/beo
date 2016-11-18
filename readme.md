# Beo: Model heat flow and (U-Th)/He thermochronology in a hydrothermal system




:Authors: Elco Luijendijk, McGill University & Goettingen University, 2013-2014
 

# Downloading Beo

* click the download link on the right for a zip file of the source code
* or clone the repository

# Installing & running 

* Install Escript

    - get the code here: https://launchpad.net/escript-finley
    - an installation guide can be found here: http://esys.geocomp.uq.edu.au/docs

* Unzip the beo source code  
* Navigate to the directory where you have installed escript. Go to the subdirectory bin (``python-escript/bin/``) and run Beo by executing the following command from the command line:
	
	>>> ./run-escript beo_dir/beo.py
	
	Where ``beo_dir`` is the directory where you have saved Beo.

# Required modules

* numpy:  http://sourceforge.net/projects/numpy/files/NumPy/
* For creating figures (optional):

    - scipy: http://sourceforge.net/projects/scipy/files/scipy/
    - matplotlib: http://matplotlib.org/downloads.html

These modules are available as standalone packages. For mac and windows an easy way to get a working version of python and these modules is to install the Enthought Canopy (https://www.enthought.com/products/canopy) or pythonxy (https://code.google.com/p/pythonxy/wiki/Welcome).

Beo was tested on Ubuntu 14.04 and 16.04 


# Model input & output

## Input:
All model input is contained in 2 files located in the directory ``model_parameters``. 

* ``model_parameters.py`` : contains all model parameters and options. see comments in this file for an explanation
* optionally you can include a file ``parameter_ranges.py`` that contains all lists of parameters values for any parameter that is included in the ``model_parameters.py`` file. Beo will run any parameter combination automatically. This is useful for automated sensitivity or uncertainty analysis. 


## Output

* After each model run, the modeled temperature field and (U-Th)/He data are stored in the directory ``model_output`` as a .pck file, which can be read using the Python's pickle module. 
* In addition, Beo saves a comma separated file containing the model parameters and a summary of the results for each model run in the same directory.


##Making figures
Currently there are two Python scripts that will generate figures of the model output:

* The script ``make_figures.py`` will make a single figure of the final temperature field for each output file (with extension .pck) found in the directory ``model_output``
* ``make_figure_2models.py`` will make a figure with two panels containing the modeled temperatures and (U-Th)/He data for two model runs. THe model runs are specified in the script itself, for example like this: ``files = ['model_output/T_field_model_run_14_(-3000.0, 0.03).pck', 'model_output/T_field_model_run_32_(-6000.0, 0.03).pck']``


