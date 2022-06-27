# Validation of  peformance of Lightning Imaging Sensor onboard International Space Station (ISS-LIS) using lightning mapping array (LMA) datasets

## Dependences
You need to download the pyltg package developed by Dr. Phillip Bitzer 
<https://www.nsstc.uah.edu/users/phillip.bitzer/python_doc/pyltg/readme_link.html#installing>

## Overview of this repo:
Scripts in thie repo serves to 
1. [Save_passover_LIS_filenames.py](https://github.com/safelysparky/LIS_LMA_val/blob/main/Save_passover_LIS_filenames.py) find LIS .nc files that contain LIS flashes over a LMA network 
2. find LIS events matching the flashes detected by the LMA, and save all matches in a dictionary (LIS_LMA_match.py)
3. calculate the detection efficiency of LIS using LMA flashes as ground-truth. (LIS_DE_analysis_against_LMA.py)

## Further info regarding each script:
1. [Save_passover_LIS_filenames.py](https://github.com/safelysparky/LIS_LMA_val/blob/main/Save_passover_LIS_filenames.py) looks for .nc files (in a specified folder) containing LIS flashes over the LMA network given the LMA center coordinates and the searching radius in km. Paths of .nc files with LIS flashes over the LMA will be saved in a .txt file.
## Contact
Please contact Yanan Zhu <yz0022@uah.edu> if you found any issue or have any questions. 