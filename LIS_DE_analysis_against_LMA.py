import numpy as np
import matplotlib.pyplot as plt
import pickle

def load_obj(name):
    with open(name, 'rb') as f:
        return pickle.load(f)
    
    
M=load_obj('E:/OKLMA_LIS/OKLMA_LIS_matches.pkl')

num_LIS_detection=0
for i, m in M.items():
    if m['LIS detection'] is True:
        num_LIS_detection+=1

DE=np.around(num_LIS_detection/len(M),2)

print(f"LIS detected {num_LIS_detection}/{len(M)} LMA flashes, DE is {DE}")

