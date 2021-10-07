# -*- coding: utf-8 -*-
"""
Created on Fri Jul 30 22:01:25 2021

@author: 25016
"""

import pandas as pd 
import numpy as np
PE = np.array(np.genfromtxt('/Users/25016/Desktop/aaaa/PE.csv', delimiter=','))
TE = np.array(np.genfromtxt('/Users/25016/Desktop/aaaa/TE.csv', delimiter=','))
ratio = pd.DataFrame(PE/TE).to_csv('/Users/25016/Desktop/aaaa/ratio.csv')
