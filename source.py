
# coding: utf-8

# ## How to use the UniRep mLSTM "babbler". This version demonstrates the 64-unit and the 1900-unit architecture. 
# 
# We recommend getting started with the 64-unit architecture as it is easier and faster to run, but has the same interface as the 1900-unit one.

# Use the 64-unit or the 1900-unit model?

# In[ ]:


USE_FULL_1900_DIM_MODEL = False # if True use 1900 dimensional model, else use 64 dimensional one.


# ## Setup

# In[ ]:


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

# In[ ]:


batch_size = 12
b = babbler(batch_size=batch_size, model_path=MODEL_WEIGHT_PATH)


# ## Code

# Initialize database and save in binary files with numpy.save("path to files", data to save)
# Once the binary files are created, there is no need to execute this code

# In[ ]:


total_start_time = time.time()
classes_avg, classes_concat = utils.dic_init()
total_elapsed_time = time.time() - total_start_time
print(total_elapsed_time)
np.save("dataset/avg/data_avg.npy", classes_avg)
np.save("dataset/concat/data_concat.npy", classes_concat)

dist_intra_avg, stats_intra_avg = utils.get_dist_intra(classes_avg)
dist_intra_concat, stats_intra_concat = utils.get_dist_intra(classes_concat)
np.save("dataset/avg/dist_intra_avg.npy", dist_intra_avg)
np.save("dataset/concat/dist_intra_concat.npy", dist_intra_concat)
np.save("dataset/avg/stats_intra.npy", stats_intra_avg)
np.save("dataset/concat/stats_intra.npy", stats_intra_concat)

dist_extra_avg, stats_extra_avg = utils.get_dist_extra(classes_avg)
dist_extra_concat, stats_extra_concat = utils.get_dist_extra(classes_concat)
np.save("dataset/avg/dist_extra_avg.npy", dist_extra_avg)
np.save("dataset/concat/dist_extra_concat.npy", dist_extra_concat)
np.save("dataset/avg/stats_extra.npy", stats_extra_avg)
np.save("dataset/concat/stats_extra.npy", stats_extra_concat)

utils.seuil_init()


# Load data with numpy.load("path to binary file containing data")

# In[ ]:


classes_avg = np.load("dataset/avg/data_avg.npy")[()]
dist_intra_avg = np.load("dataset/avg/dist_intra_avg.npy")[()]
dist_extra_avg = np.load("dataset/avg/dist_extra_avg.npy")[()]
stat_intra_avg = np.load("dataset/avg/stats_intra.npy")
seuil_avg = np.load("dataset/avg/seuil.npy")[()]
seuil_avg2 = np.load("dataset/avg/seuil2.npy")[()]
stat_extra_avg = np.load("dataset/avg/stats_extra.npy")

classes_concat = np.load("dataset/concat/data_concat.npy")[()]
dist_intra_concat = np.load("dataset/concat/dist_intra_concat.npy")[()]
dist_extra_concat = np.load("dataset/concat/dist_extra_concat.npy")[()]
stat_intra_concat = np.load("dataset/concat/stats_intra.npy")
stat_extra_concat = np.load("dataset/concat/stats_extra.npy")


# In[ ]:


utils.histo(dist_intra_avg, dist_extra_avg, avg = True)


# In[ ]:


utils.histo(dist_intra_concat, dist_extra_concat, avg = False)


# In[ ]:


res = np.inf
for distance in dist_intra_avg["a.1.1.1"].values():
    print(distance[1])
    if distance[1] < res:
        res = distance[1]
print(res)
print("-----------")
seuil_intra = np.inf
for protein, distance in dist_intra_avg["a.1.1.1"].values():
    print(distance)
    if distance < seuil_intra:
        seuil_intra = distance
print(seuil_intra)


# In[ ]:


utils.seuil_init()


# In[ ]:


utils.seuil_init2()


# In[ ]:


utils.psiblastCode(classes_avg, seuil_avg2)


# In[ ]:


print("INTRA")
x = 0
print(dist_intra_avg["a.104.1.1"])
for key, val in dist_intra_avg["a.104.1.1"].items():
    print(val[1])
    if val[1] > x:
        x = val[1]
print("max :", x)

print("EXTRA")
y = 0
for key, val in dist_extra_avg["a.104.1.1"].items():
    print(val[1])
    if val[1] > y:
        y = val[1]
print("max :",y)

