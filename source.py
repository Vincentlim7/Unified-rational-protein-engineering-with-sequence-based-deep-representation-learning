
# coding: utf-8

# ## How to use the UniRep mLSTM "babbler". This version demonstrates the 64-unit and the 1900-unit architecture. 
# 
# We recommend getting started with the 64-unit architecture as it is easier and faster to run, but has the same interface as the 1900-unit one.

# Use the 64-unit or the 1900-unit model?

# In[1]:


USE_FULL_1900_DIM_MODEL = False # if True use 1900 dimensional model, else use 64 dimensional one.


# ## Setup

# In[4]:


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

# In[5]:


batch_size = 12
b = babbler(batch_size=batch_size, model_path=MODEL_WEIGHT_PATH)


# In[18]:


from scipy.spatial import distance
# PLUS BESOIN DE get_rd_prot() JE CROIS
def get_rd_prot():    # Retourne une ligne aleatoire du fichier fullProtein.list
    lines = open("fullProtein.list").read().splitlines()
    return np.random.choice(lines).split() # return ("protein name", "protein class)

def get_prot_seq(file_name):
    f = open("dataset/fastas/" + file_name + ".fasta", "r")
    next(f)
    seq = ""
    for line in f:
        tmp = line.rstrip()    # Supprimer le "\n"
        seq += tmp
    f.close
    return seq

def get_avg_vec(seq):
    avg_vec = b.get_rep(seq)[0]
    return avg_vec

def get_concat_vec(seq):
    avg_vec = b.get_rep(seq)[0]
    fnl_hid_vec = b.get_rep(seq)[1]
    fnl_cell_vec = b.get_rep(seq)[2]
    seq_vec = np.concatenate((avg_vec, fnl_hid_vec, fnl_cell_vec))
    return seq_vec

def get_classe(searched_protein):
    for classe, protein_list in classes.items():
        for protein_name, seq in protein_list.items():
            if protein_name == searched_protein:
                return classe


def dic_init(avg = True):
    classes = dict()
    f = open("partialProtein.list", "r")
    for line in f:
        infos = line.split()
        protein = infos[0]    # Protein name
        classe = infos[1]     # Protein class
        if classe not in classes:
            classes[classe] = dict()
        if avg:
            classes[classe][protein] = get_avg_vec(get_prot_seq(protein))
        else:
            classes[classe][protein] = get_concat_vec(get_prot_seq(protein))
    return classes

def dist_intra(protein_dict): # Fonctionne 
    dist_intra = dict()
    for classe, protein_list in protein_dict.items():
        if classe not in dist_intra:
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

def dist_extra(protein_dict): # A CODER 
    dist_extra = dict()
    for classe, protein_list in protein_dict.items():
        for protein_a, seq_a in protein_list.items():
            dist_extra[protein_a] = (None, np.inf)
            for protein_b, seq_b in protein_list.items():
                if protein_a == protein_b:
                    continue
                dist = distance.euclidean(seq_a, seq_b)
                if dist < dist_extra[protein_a][1]:
                    dist_extra[protein_a] = (protein_b, dist)
    return dist_extra
                    
                


# In[69]:


# EXEMPLE VECTEURS
seq = "MRKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFARYPDHMKQHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
print(get_avg_vec(seq))
print(get_concat_vec(seq))


# In[9]:


classes = dic_init()
print(classes)


# In[15]:


# EXEMPLE DISTANCE EUCLIDIENNE
prot_a = classes["a.1.1.4"]["d1kr7a_"]
prot_b = classes["a.1.1.1"]["d1ngka_"]
print(prot_a)
print(prot_b)
print(distance.euclidean(prot_a, prot_b))


# In[19]:


dist_intra = dist_intra(classes)
print(dist_intra)


# In[59]:


from scipy.spatial import distance

f = open("fullProtein.list")

line = next(f).split()
studied_prot = get_seq(line[0])

same_class = []
rd_class = []

nb_seq_class = 5


for i in range (nb_seq_class): #Recuperer les noms de fichiers
    line = next(f).split()
    same_class.append(line[0])
    rd_protein =random_line("fullProtein.list")
    rd_class.append(line[0])
    
f.close

print(same_class)
print(rd_class)

seq_list1 = [] # Recuperer les sequences
dist1 = []
for i in range(nb_seq_class):
    prot = get_seq(same_class[0])
    seq_list1.append(prot)
    #dist1.append(distance.euclidean(studied_prot[1], prot[1]))

seq_list2 = []
dist2 = []
for i in range(nb_seq_class):
    prot = get_seq(rd_class[0])
    seq_list2.append(prot)
    #dist2.append(distance.euclidean(studied_prot[1], prot[1]))

print(seq_list1)

"""rd_protein = random_line("fullProtein.list")
print(random_protein)

file_name = "dataset/fastas/" + protein[0] + ".fasta"
rd_file_name = "dataset/fastas/" + rd_protein[0] + ".fasta"

seq, seq_vector = get_seq(file_name)
rd_seq, rd_seq_vector = get_seq(rd_file_name)

print("Studied sequence :")
print(seq + "\n")
print("Sequence size : " + str(len(seq)) + "\n")
print("Vector representation of the sequence :\n" + str(seq_vector) + "\n")
print("Vector size : " + str(len(seq_vector)) + "\n")
print("Is seq a valide sequence ?\n" + str(b.is_valid_seq(seq)) + "\n")
"""

