
# coding: utf-8

# ## How to use the UniRep mLSTM "babbler". This version demonstrates the 64-unit and the 1900-unit architecture. 
# 
# We recommend getting started with the 64-unit architecture as it is easier and faster to run, but has the same interface as the 1900-unit one.

# Use the 64-unit or the 1900-unit model?

# In[ ]:


USE_FULL_1900_DIM_MODEL = False # if True use 1900 dimensional model, else use 64 dimensional one.


# ## Setup

# In[ ]:


import tensorflow as tf
import numpy as np
from scipy.spatial import distance
import matplotlib.pyplot as plt
import time

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


# In[ ]:


def get_prot_seq(file_name):
    start_time_seq = time.time()

    f = open("dataset/fastas/" + file_name + ".fasta", "r") # Retriving the file containing the sequence
    next(f) # Skipping the first line (containing the protein's name)
    seq = ""
    for line in f: # Retriving the sequence
        tmp = line.rstrip()    # Deleting "\n"
        seq += tmp
    f.close
    elapsed_time_seq = time.time() - start_time_seq
    print("fichier lu, temps : ", elapsed_time_seq)
    return seq

def get_avg_vec(seq):
    start_time_avg = time.time()
    avg_vec = b.get_rep(seq)[0] # Vector 1 : avg
    elapsed_time_avg = time.time() - start_time_avg
    print("vecteur avg, temps : ", elapsed_time_avg)
    return avg_vec

def get_concat_vec(seq):
    start_time_concat = time.time()
    seq_vec = b.get_rep(seq)
    seq_vec = np.reshape(seq_vec, 192) # Concatenation of all three vectors
    elapsed_time_concat = time.time() - start_time_concat
    print("vecteur concat, temps : ", elapsed_time_concat)
    return seq_vec

def get_classe(searched_protein): # Returning the protein's category (the key in the level 0 dictionnary)
    for classe, protein_list in classes.items(): # Browsing the category dictionnary (level 0)
        for protein_name, seq in protein_list.items(): # Browsing the protein dictionnary (level 1)
            if protein_name == searched_protein:
                return classe


def dic_init(): # Initializing nested dictionnaries all proteins and their vector (avg or concatenated)
    classes_avg = dict()
    classes_concat = dict()
    f = open("partialProtein.list", "r")
    cpt = 0 # Number of protein already processed
    start_time = time.time()
    for line in f: # Browsing all protein
        infos = line.split()
        protein = infos[0]    # Protein name
        classe = infos[1]     # Protein category
        if classe not in classes_avg: # Adding new category key if it doesn't exist
            classes_avg[classe] = dict()
        if classe not in classes_concat: # Adding new category key if it doesn't exist
            classes_concat[classe] = dict()
        prot_seq = get_prot_seq(protein) # Retrieving protein sequence
        classes_avg[classe][protein] = get_avg_vec(prot_seq) # Retrieving and stocking avg vector
        classes_concat[classe][protein] = get_concat_vec(prot_seq) # Retrieving and stocking concat vector
        cpt += 1
        if cpt % 100 == 0: # Periodical save every 100 protein processed
            print("Nombre proteines lues : ", cpt)
            np.save("database/next_batch/data_avg" + str(cpt) + ".npy", classes_avg)
            np.save("database/next_batch/data_concat" + str(cpt) + ".npy", classes_concat)
            elapsed_time = time.time() - start_time
            print(elapsed_time)
            start_time = time.time() # Reset timer
    np.save("database/next_batch/data_avg.npy", classes_avg) # Saving whole database
    np.save("database/next_batch/data_concat.npy", classes_concat)
    return classes_avg, classes_concat

def get_dist_intra(protein_dict): # Initializing a dictionnary containning the shortest euclidian distance between proteins of the same category
    dist_intra = dict()
    for classe, protein_list in protein_dict.items():
        if classe not in dist_intra: # Adding new category key if it doesn't exist
            dist_intra[classe] = dict()
        for protein_a, vec_a in protein_list.items():
            dist_intra[classe][protein_a] = (None, np.inf)
            for protein_b, vec_b in protein_list.items():
                if protein_a == protein_b:
                    continue
                dist = distance.euclidean(vec_a, vec_b)
                if dist < dist_intra[classe][protein_a][1]:
                    dist_intra[classe][protein_a] = (protein_b, dist)
    return dist_intra

def get_dist_extra(protein_dict): # Initializing a dictionnary containning the shortest euclidian distance between proteins of different category
    dist_extra = dict()
    for classe_a, protein_list_a in protein_dict.items():
        if classe_a not in dist_extra: # Adding new category key if it doesn't exist
            dist_extra[classe_a] = dict()
        for protein_a, vec_a in protein_list_a.items():
            dist_extra[classe_a][protein_a] = (None, np.inf)
            for classe_b, protein_list_b in protein_dict.items():
                if classe_a == classe_b:
                    continue
                for protein_b, vec_b in protein_list_b.items():
                    dist = distance.euclidean(vec_a, vec_b)
                    if dist < dist_extra[classe_a][protein_a][1]:
                        dist_extra[classe_a][protein_a] = (protein_b, dist)
    return dist_extra
                    
def histo(dist_intra, dist_extra):
    x_intra = []
    y_intra = []
    x_extra = []
    y_extra = []

    for classe, protein_list in dist_intra.items(): # Retrieving smallest dist value in dist_intra for each protein
        classe_dist = np.inf
        for protein, val in protein_list.items():
            if(val[1] < classe_dist):
                classe_dist = val[1]
        if(classe_dist != np.inf):
            x_intra.append(classe_dist)
            y_intra.append(len(protein_list))

    for classe, protein_list in dist_extra.items(): # Retrieving smallest dist value in dist_extra for each protein
        classe_dist = np.inf
        for protein, val in protein_list.items():
            if(val[1] < classe_dist):
                classe_dist = val[1]
        if(classe_dist != np.inf):
            x_extra.append(classe_dist)
            y_extra.append(len(protein_list))

    plt.bar(x_intra,y_intra,align='center', alpha = 0.7, width = 0.01, label='intra')
    plt.bar(x_extra,y_extra,align='center', alpha = 0.7, width = 0.01, label='extra')
    plt.xlabel('Distance')
    plt.ylabel('Nb Sequence')
    plt.legend(loc='upper right')
    plt.show()


# In[ ]:


total_start_time = time.time()
classes_avg, classes_concat = dic_init()
total_elapsed_time = time.time() - total_start_time
print(total_elapsed_time)


# In[ ]:


classes_avg = np.load("database/data_avg1294.npy")[()]
classes_concat = np.load("database/data_concat1294.npy")[()]


# In[ ]:


dist_intra_avg = get_dist_intra(classes_avg)
dist_intra_concat = get_dist_intra(classes_concat)


# In[ ]:


dist_extra_avg = get_dist_extra(classes_avg)
dist_extra_concat = get_dist_extra(classes_concat)


# In[ ]:


histo(dist_intra_avg, dist_extra_avg)


# In[ ]:


histo(dist_intra_concat, dist_extra_concat)


# In[ ]:


cpt = 0
for classe, protein_list in classes_avg.items():
    for protein_a, vec_a in protein_list.items():
        cpt += 1
print(cpt)

