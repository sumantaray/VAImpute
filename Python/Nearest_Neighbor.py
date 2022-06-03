#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy
import scipy
from sklearn.neighbors import LSHForest
#numpy.set_printoptions(threshold=numpy.nan)
from numpy import genfromtxt
from scipy import io, sparse
import scipy.stats
from scipy.io import mmread
from multiprocessing import Pool


# In[ ]:


Xnew = genfromtxt('~path/data/data_process.csv',delimiter=",")


# In[ ]:


lshf = LSHForest(n_estimators=30, random_state=42)
lshf.fit(sparse.coo_matrix(X))
distances, indices = lshf.kneighbors(X, n_neighbors=11)


# In[ ]:


numpy.savetxt('~path/data/NNcell_data.csv', indices, delimiter=',') 

