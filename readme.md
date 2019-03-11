# Beo: Model heat flow and (U-Th)/He thermochronology in a hydrothermal system

*Elco Luijendijk, University of Göttingen*


# Introduction

Beo is a model of heat flow in hot springs and hydrothermal systems. The model code uses the generic finite element code escript (https://launchpad.net/escript-finley) to solve the advective and conductive heat flow equations in a 2D cross-section of the subsurface. The resulting temperature history is used to calculate the low-temperature thermochronometer apatite (U-Th)/He. The modeled values of this thermochronometers can then be compared to measured values. Beo support automated model runs to explore which parameter values like fluid fluxes, fault geometry and age of the hydrothermal system best match the thermochronometer data, as well as present-day spring temperature data or temperature records in nearby boreholes. 

 

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

This will show you three lines that define environment variables that your system needs to be able to find the location of escript. Add these lines to your .bashrc (Ubuntu linux) or profile file in your home directory. After adding these lines and logging out and in again, you can start beo by going to the directory where the beo code is located (so not to the escript/bin directory) and start beo.py like any regular python code:

````bash
python beo.py model_parameters/model_parameters.py
````

where ``model_parameters/model_parameters.py`` is a file containing all model parameters. An example file called ``model_parameters.py`` is located in the directory ``model_parameters``.


# Required modules

* numpy:  http://www.numpy.org/
* For creating figures (optional):

    - scipy: http://scipy.org/scipylib/index.html
    - matplotlib: http://matplotlib.org/downloads.html

These modules are available as standalone packages. For mac and windows an easy way to get a working version of python and these modules is to install a full Python environemnt like Anaconda (https://www.anaconda.com/), Enthought Canopy (https://www.enthought.com/products/canopy) or pythonxy (https://code.google.com/p/pythonxy/wiki/Welcome).

Note that Beo includes an option to calculate apatite (U-Th)/He ages using the RDAAM model (Flowers et al. 2009, Geochimica et Cosmochimica Acta 73(8)). The implementation of the RDAAM model uses a piece of Fortran code to speed up the model. To enable the RDAAM model you first need to compile the Fortran code using f2py, which is normally inlcuded with Numpy (https://docs.scipy.org/doc/numpy/f2py/). To do so, navigate to the subdirectory lib and run the following command:

``f2py -c calculate_reduced_AFT_lengths.f90 -m calculate_reduced_AFT_lengths``

If all goes well there should now be a file called ``calculate_reduced_AFT_lengths.so`` in the subdirectory lib, which will be used by Beo to calculate radiation damage and the diffusivity of helium in apatites. Currently new updates are planned to convert the fortran code to python, which will remove this requirement.  

Beo was tested on Ubuntu 14.04 and 16.04 


# Manual and publication

You can find a manual for beo in the directory manual. The manual contains more background on the model, explanation of how surface heat flow is modeled and a detailed list and explanation of the model paramameters. More information on the model code can be found in a discussion paper that is currently under review: https://doi.org/10.5194/gmd-2018-341. See the bottom of this readme for the full reference. 

The example model runs for the Baden & Schinznach hot springs that are discussed in this paper can be reproduced by using one of the parameter files located in the directory ``example_input_files``.


# Model input & output

## Input:

All model input parameters are contained in a single Python file. An example file can be found in ``model_parameters.py`` located in the directory model_parameters. The class ``ModelParameters`` contains all parameters needed for a single model run. See the manual for an explanation of the model parameters.


## Multiple model runs

Optionally you can start automated runs to test a range of parameter combinations. This is useful for automated sensitivity or uncertainty analysis. 

THe model input file contains a class called ``ParameterRanges``. Any parameter value that is included in this class will be used as input for a single model run. All results will be stored and written to a comma separated (.csv) file names ``model_output/model_params_and_results_x_runs.csv``. 

You can include any model parameter in the automated runs. Simply copy a parameter from the ``ModelParameters`` class into the ``ParameterRanges`` class and add _s to the parameter name. For instance to test multiple values of the thermal gradient, add `thermal_gradient_s = [0.03, 0.04, 0.05]` to test the effect of geothermal gradients of 0.03, 0.04 and 0.05 C/m, respectively.

There are two options for running multiple model runs. The default is a sensitivity run. In each model run a single parameter will be changed, while all other parameters are kept constant at the default value specified in the ``ModelParameters`` class. Alternatively you can test all parameter combinations by changing `parameter_combinations = False` to `parameter_combinations = True`. Note that this will generate a lot of model runs, testing ten parameter values for two parameters each will generate 10*10 = 100 model runs, for three parameters this increase to a 1000 model runs, etc...


## Output

* After each model run, the modeled temperature field and (U-Th)/He data are stored in the directory ``model_output`` as a .pck file, which can be read using Python's pickle module. 
* Beo saves a comma separated file containing the input parameters and a summary of the results for each model run and each timestep in the same directory.
* Beo also contains an option to save modeled temperature and advective flux to a VTK file, which can be used for visualization using software such as Paraview and Visit.


## Making figures

The script ``make_figures.py`` will make a single figure of the final temperature field for output files (with extension .pck) found in the directory ``model_output``. After running this script you will be prompted to select the output file that you would like a figure of. The file ``model_parameters/figure_params.py`` contains a number of parameters that control the figure, such as which timeslices to show, the min. and max. coordinates of the area to show, etc.. THe resulting figure is saved as a .png file in the same directory as the model output file.


# License

Beo is distributed under the GNU General Public License, version 3:
http://www.gnu.org/copyleft/gpl.html

A copy of this license is distributed along with the source code, see LICENSE.txt


# Reference

Please cite the following paper if you publish work that uses Beo:

Luijendijk, E.: Beo v1.0: Numerical model of heat flow and low-temperature thermochronology in hydrothermal systems, Geosci. Model Dev. Discuss., https://doi.org/10.5194/gmd-2018-341, in review, 2019. 

The paper is freely available and can be found here: https://doi.org/10.5194/gmd-2018-341