
# coding: utf-8

# ## How to use the UniRep mLSTM "babbler". This version demonstrates the 64-unit and the 1900-unit architecture. 
# 
# We recommend getting started with the 64-unit architecture as it is easier and faster to run, but has the same interface as the 1900-unit one.

# Use the 64-unit or the 1900-unit model?

# In[1]:


USE_FULL_1900_DIM_MODEL = False # if True use 1900 dimensional model, else use 64 dimensional one.


# ## Setup

# In[2]:


# to allow autoreload of utils.py
get_ipython().magic('load_ext autoreload')
get_ipython().magic('autoreload 2')

import utils # our functions
import tensorflow as tf
import numpy as np

# Set seeds
tf.set_random_seed(42)
np.random.seed(42)

if USE_FULL_1900_DIM_MODEL:
    # Sync relevant weight files
    get_ipython().system('aws s3 sync --no-sign-request --quiet s3://unirep-public/1900_weights/ 1900_weights/')
    
    # Import the mLSTM babbler model
    from unirep import babbler1900 as babbler
    
    # Where model weights are stored.
    MODEL_WEIGHT_PATH = "./1900_weights"
    
else:
    # Sync relevant weight files
    get_ipython().system('aws s3 sync --no-sign-request --quiet s3://unirep-public/64_weights/ 64_weights/')
    
    # Import the mLSTM babbler model
    from unirep import babbler64 as babbler
    
    # Where model weights are stored.
    MODEL_WEIGHT_PATH = "./64_weights"


# ## Data formatting and management

# Initialize UniRep, also referred to as the "babbler" in our code. You need to provide the batch size you will use and the path to the weight directory.

# In[3]:


batch_size = 12
b = babbler(batch_size=batch_size, model_path=MODEL_WEIGHT_PATH)


# ## Code

# Initialize database and save in binary files with numpy.save("path to files", data to save)
# Once the binary files are created, there is no need to execute this code

# In[ ]:


total_start_time = time.time()
classes_avg, classes_concat = dic_init()
total_elapsed_time = time.time() - total_start_time
print(total_elapsed_time)
np.save("dataset/avg/data_avg.npy", classes_avg)
np.save("dataset/concat/data_concat.npy", classes_concat)

dist_intra_avg = get_dist_intra(classes_avg)
dist_intra_concat = get_dist_intra(classes_concat)
np.save("dataset/avg/dist_intra_avg.npy", dist_intra_avg)
np.save("dataset/concat/dist_intra_concat.npy", dist_intra_concat)

dist_extra_avg = get_dist_extra(classes_avg)
dist_extra_concat = get_dist_extra(classes_concat)
np.save("dataset/avg/dist_extra_avg.npy", dist_extra_avg)
np.save("dataset/concat/dist_extra_concat.npy", dist_extra_concat)

_, stats_intra_avg = get_dist_extra(classes_avg)
_, stats_intra_concat = get_dist_extra(classes_concat)
np.save("dataset/avg/stats_intra.npy", stats_intra_avg)
np.save("dataset/concat/stats_intra.npy", stats_intra_concat)

_, stats_extra_avg = get_dist_extra(classes_avg)
_, stats_extra_concat = get_dist_extra(classes_concat)
np.save("dataset/avg/stats_extra.npy", stats_extra_avg)
np.save("dataset/concat/stats_extra.npy", stats_extra_concat)


# Load data with numpy.load("path to binary file containing data")

# In[4]:


classes_avg = np.load("dataset/avg/data_avg.npy")[()]
dist_intra_avg = np.load("dataset/avg/dist_intra_avg.npy")[()]
dist_extra_avg = np.load("dataset/avg/dist_extra_avg.npy")[()]

classes_concat = np.load("dataset/concat/data_concat.npy")[()]
dist_intra_concat = np.load("dataset/concat/dist_intra_concat.npy")[()]
dist_extra_concat = np.load("dataset/concat/dist_extra_concat.npy")[()]


# In[9]:


stat_intra_avg = np.load("dataset/avg/stats_intra.npy")
stat_intra_concat = np.load("dataset/concat/stats_intra.npy")

stat_extra_avg = np.load("dataset/avg/stats_extra.npy")
stat_extra_concat = np.load("dataset/concat/stats_extra.npy")


# In[11]:


print(stat_intra_avg)
print(stat_extra_avg)
print(stat_intra_concat)
print(stat_extra_concat)


# In[5]:


utils.histo(dist_intra_avg, dist_extra_avg, avg = True)


# In[6]:


utils.histo(dist_intra_concat, dist_extra_concat, avg = False)


# In[64]:


utils.test()

