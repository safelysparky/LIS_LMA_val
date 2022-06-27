# Validation of  peformance of Lightning Imaging Sensor onboard International Space Station (ISS-LIS) using lightning mapping array (LMA) datasets

## Dependences other than some standard packages (e.g., numpy, pandas, shapely)
pyltg package for some handy functions for LIS data
<https://www.nsstc.uah.edu/users/phillip.bitzer/python_doc/pyltg/readme_link.html#installing>

astral for determinting the sunrise/sunset time at a give location
<https://astral.readthedocs.io/en/latest/>

pymap3d for (lat,lon,alt) <==> (x,y,z) coordinates conversion
<https://github.com/geospace-code/pymap3d>



## Overview of this repo:
Scripts in thie repo serves to 
1. [Save_passover_LIS_filenames.py](https://github.com/safelysparky/LIS_LMA_val/blob/main/Save_passover_LIS_filenames.py) find LIS .nc files that contain LIS flashes over a LMA network. 
2. [LIS_LMA_match.py](https://github.com/safelysparky/LIS_LMA_val/blob/main/LIS_LMA_match.py) find LIS events matching the LMA flashes, and save all matches in a dictionary.
3. [LIS_DE_analysis_against_LMA.py](https://github.com/safelysparky/LIS_LMA_val/blob/main/LIS_DE_analysis_against_LMA.py) calculate the detection efficiency of LIS using LMA flashes as ground-truth.

## Further info regarding each script:
1. [Save_passover_LIS_filenames.py](https://github.com/safelysparky/LIS_LMA_val/blob/main/Save_passover_LIS_filenames.py) looks for .nc files (in a specified folder) containing LIS flashes over the LMA network given the LMA center coordinates and the searching radius in km. Paths of .nc files with LIS flashes over the LMA will be saved in a .txt file.
2. [LIS_LMA_match.py](https://github.com/safelysparky/LIS_LMA_val/blob/main/LIS_LMA_match.py) reads the .nc files whose filenames are saved in a .txt file by the [Save_passover_LIS_filenames.py](https://github.com/safelysparky/LIS_LMA_val/blob/main/Save_passover_LIS_filenames.py). During the period of LIS passover the LMA, the correspoidng LMA file(s) (in a specified folder) will be searched and loaded (if found). LMA L1 data will be sorted into flashes using a spatial and temporal threshold. There is addtional option to include ENTLN CG datasets to determine the type (IC/CG) of each LMA flash. Each of the sorted lma flashes will be further checked if they are within the FOV of LIS. If yes, the corresponding LIS events will be searched and saved (if any). A lma flash is considered detected by the LIS if any of the LIS event pixel intersect with the LMA convext hull (with some buffer to account for the LIS location offset). Each of the LMA flashes's L1 data will be saved in a dictionary, along with coreesponding LIS events (if detected by LIS), and also ENTLN CG data (if any). 
3. [LIS_DE_analysis_against_LMA.py](https://github.com/safelysparky/LIS_LMA_val/blob/main/LIS_DE_analysis_against_LMA.py) To be finished...
## Contact
Please contact Yanan Zhu <yz0022@uah.edu> if you found any issue or have any questions. 