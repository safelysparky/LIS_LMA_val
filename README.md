# Validation of  peformance of Lightning Imaging Sensor onboard International Space Station (ISS-LIS) using lightning mapping array (LMA) datasets

## Dependences
pyltg package developed by Dr. Phillip Bitzer for some handy functions for LIS data
<https://www.nsstc.uah.edu/users/phillip.bitzer/python_doc/pyltg/readme_link.html#installing>

astral for determinting the sunrise/sunset time at give location
<https://astral.readthedocs.io/en/latest/>

pymap3d for (lat,lon,alt) ==> (x,y,z) coordinates conversion
<https://github.com/geospace-code/pymap3d>



## Overview of this repo:
Scripts in thie repo serves to 
1. [Save_passover_LIS_filenames.py](https://github.com/safelysparky/LIS_LMA_val/blob/main/Save_passover_LIS_filenames.py) find LIS .nc files that contain LIS flashes over a LMA network. 
2. [LIS_LMA_match.py](https://github.com/safelysparky/LIS_LMA_val/blob/main/LIS_LMA_match.py) find LIS events matching the LMA flashes, and save all matches in a dictionary.
3. [LIS_DE_analysis_against_LMA.py](https://github.com/safelysparky/LIS_LMA_val/blob/main/LIS_DE_analysis_against_LMA.py) calculate the detection efficiency of LIS using LMA flashes as ground-truth.

## Further info regarding each script:
1. [Save_passover_LIS_filenames.py](https://github.com/safelysparky/LIS_LMA_val/blob/main/Save_passover_LIS_filenames.py) looks for .nc files (in a specified folder) containing LIS flashes over the LMA network given the LMA center coordinates and the searching radius in km. Paths of .nc files with LIS flashes over the LMA will be saved in a .txt file.
2. [LIS_LMA_match.py](https://github.com/safelysparky/LIS_LMA_val/blob/main/LIS_LMA_match.py) reads the .nc files whose filenames are saved by the [Save_passover_LIS_filenames.py](https://github.com/safelysparky/LIS_LMA_val/blob/main/Save_passover_LIS_filenames.py).
## Contact
Please contact Yanan Zhu <yz0022@uah.edu> if you found any issue or have any questions. 