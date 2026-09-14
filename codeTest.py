import numpy as np 
import pandas as pd 
from matpotlib import pyplot as plt

for k in range(5):
    x = np.linspace(0, 10, 100)
    y = np.sin(x + k)
    plt.plot(x, y, label=f'Sin wave {k}')