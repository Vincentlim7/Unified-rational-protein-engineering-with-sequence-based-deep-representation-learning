import numpy as np
from scipy.spatial import distance
import matplotlib.pyplot as plt
import time
import os

    
def get_prot_seq(file_name):
    f = open("dataset/fastas/" + file_name + ".fasta", "r") # Retriving the file containing the sequence
    next(f) # Skipping the first line (containing the protein's name)
    seq = ""
    for line in f: # Retriving the sequence
        tmp = line.rstrip()    # Deleting "\n"
        seq += tmp
    f.close
    
    return seq

def get_vecs(seq, b):
    vecs = b.get_rep(seq)
    avg_vec = vecs[0]
    concat_vec = np.reshape(vecs, 192) # Concatenation of all three vectors
    return avg_vec, concat_vec

def get_classe(searched_protein, protein_dict): # Returning the protein's category (the key in the level 0 dictionnary)
    for classe, protein_list in protein_dict.items(): # Browsing the category dictionnary (level 0)
        for protein_name, seq in protein_list.items(): # Browsing the protein dictionnary (level 1)
            if protein_name == searched_protein:
                return classe


def dic_init(b): # Initializing nested dictionnaries all proteins and their vector (avg or concatenated)
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
        classes_avg[classe][protein], classes_concat[classe][protein] = get_vecs(prot_seq, b) # Retrieving and stocking vectors
        cpt += 1
        if cpt % 100 == 0: # Periodical save every 100 protein processed
            print("Nombre proteines lues : ", cpt)
            elapsed_time = time.time() - start_time
            print(elapsed_time)
            start_time = time.time() # Reset timer
            np.save("database/avg/next_batch/data_avg" + str(cpt) + ".npy", classes_avg)
            np.save("database/concat/next_batch/data_concat" + str(cpt) + ".npy", classes_concat)
    np.save("database/avg/next_batch/data_avg.npy", classes_avg) # Saving whole database
    np.save("database/concat/next_batch/data_concat.npy", classes_concat)
    f.close
    return classes_avg, classes_concat

def get_dist_intra(protein_dict): # Initializing a dictionnary containning the shortest euclidian distance between proteins of the same category
    dist_intra = dict()
    allDist = []    # All distance value, used to compute mean and std
    for classe, protein_list in protein_dict.items():
        if classe not in dist_intra: # Adding new category key if it doesn't exist
            dist_intra[classe] = dict()
        for protein_a, vec_a in protein_list.items():
            dist_intra[classe][protein_a] = (None, np.inf) # Adding new protein key if it doesn't exist
            for protein_b, vec_b in protein_list.items():
                if protein_a == protein_b:
                    continue
                dist = distance.euclidean(vec_a, vec_b)
                allDist.append(dist)
                if dist < dist_intra[classe][protein_a][1]: # Checking if the new value is smaller then the one already saved
                    dist_intra[classe][protein_a] = (protein_b, dist)
    return dist_intra, (np.mean(allDist), np.std(allDist))

def get_dist_extra(protein_dict): # Initializing a dictionnary containning the shortest euclidian distance between proteins of different category
    dist_extra = dict()
    allDist = []
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
                    allDist.append(dist)
                    if dist < dist_extra[classe_a][protein_a][1]:
                        dist_extra[classe_a][protein_a] = (protein_b, dist)
    return dist_extra, (np.mean(allDist), np.std(allDist))
                    
def histo(dist_intra, dist_extra, avg):
    x_intra = []
    x_extra = []
    
    # Retrieve datas
    for classe, protein_list in dist_intra.items(): # Retrieving smallest dist value in dist_intra for each protein
        for protein, val in protein_list.items():
            if(val[1] != np.inf):     # val = (protein, dist)
                x_intra.append(val[1])

    for classe, protein_list in dist_extra.items(): # Retrieving smallest dist value in dist_extra for each protein
        for protein, val in protein_list.items():
            if(val[1] != np.inf):     # val = (protein, dist)
                x_extra.append(val[1])
    
    # Plot histogram
    if avg:
        plt.title("Distance Euclidienne avec vecteurs avg")
        plt.hist([x_intra, x_extra], bins=100, label=['intra', 'extra'])
        
    else:
        plt.title("Distance Euclidienne avec vecteurs concat")
        plt.hist([x_intra, x_extra], bins=100, label=['intra', 'extra'])
    
    
    print("distance intra :")
    print("\tmin :", min(x_intra))
    print("\tmax :", max(x_intra))
    print("\tmoyenne :", np.mean(x_intra))
    print("\tecart-type :", np.std(x_intra))
    print("distance extra :")
    print("\tmin :", min(x_extra))
    print("\tmax :", max(x_extra))
    print("\tmoyenne :", np.mean(x_extra))
    print("\tecart-type :", np.std(x_extra))
    plt.xlabel('Distance euclidienne')
    plt.ylabel('Nb Sequence')
    plt.legend(loc='upper right')
    plt.show()

def seuil_init(): # max(max_intra, max_extra)
    dist_intra = np.load("dataset/avg/dist_intra_avg.npy")[()]
    dist_extra = np.load("dataset/avg/dist_extra_avg.npy")[()]
    seuil = dict()
    for classe, dist_lis in dist_intra.items():
        max_intra = 0
        max_extra = 0
        for protein, distance in dist_intra[classe].values():
            if distance > max_intra:
                max_intra = distance
        for protein, distance in dist_extra[classe].values():
            if distance > max_extra:
                max_extra = distance
        seuil[classe] = max(max_intra, max_extra)
    np.save("dataset/avg/seuil.npy", seuil)

def seuil_init2(): # max(max_intra, min_extra)
    dist_intra = np.load("dataset/avg/dist_intra_avg.npy")[()]
    dist_extra = np.load("dataset/avg/dist_extra_avg.npy")[()]
    seuil = dict()
    for classe, dist_lis in dist_intra.items():
        max_intra = 0
        min_extra = np.inf
        for protein, distance in dist_intra[classe].values():
            if distance == np.inf:
                print("DISTANCE = INF")
            if distance > max_intra:
                max_intra = distance
        for protein, distance in dist_extra[classe].values():
            if distance < min_extra:
                min_extra = distance
        seuil[classe] = max(max_intra, min_extra)
    np.save("dataset/avg/seuil2.npy", seuil)

# dict 1 classe --> dict 2
# dict 2 protein --> vecteur
def psiblastCode(classe_dict, seuil_avg):
    start_time = time.time()
    cpt = 0     # number of familly ignored
    test_dataset = open("./dataset/test_dataset.list", "r")
    test_cpt = 0
    for test_line in test_dataset: # Retriving the sequence
        test_prot = test_line.split()
        test_prot_name = test_prot[0]
        test_prot_class = test_prot[1]
        # print("test : ", test_prot_name)
        test_cpt += 1

        if test_cpt % 50 == 0:
            elapsed_time = time.time() - start_time
            print("50 sequences test traitee en : ", elapsed_time, "s")
            start_time = time.time() # Reset timer
            print(test_cpt, "sequences test traitee")

        valid_class = None      # Used to skip unneeded comparison
        previous_class = None   # Used to check if a family is never under the threshold
        train_dataset = open("./dataset/train_dataset.list", "r")

        train_cpt = 0
        for train_line in train_dataset:
            train_prot = train_line.split()
            train_prot_name = train_prot[0]
            train_prot_class = train_prot[1]
            train_cpt += 1
            
            if train_prot_class != previous_class and previous_class != valid_class:  # Check if we are comparing with another family and the previous one wasn't valid
                # print("-------------\n","La famille", previous_class, "a ete ignore","\n-------------\n")
                cpt += 1

            # if train_prot_class != previous_class:
            #     print("Valeur du seuil pour la famille", train_prot_class, ":", seuil_avg[train_prot_class], "\n-------------")

            if train_prot_class == valid_class:
                os.system('bash scripts/psiblast_script.sh ' + test_prot_name + ' ' + train_prot_name)
                continue
            
            dist = distance.euclidean(classe_dict[test_prot_class][test_prot_name], classe_dict[train_prot_class][train_prot_name])
            # print("dist(", test_prot_name, ",", train_prot_name, ")",  ":", dist)
            if dist < seuil_avg[train_prot_class]:
                # print("Lancement de psiblast pour", test_prot_name, "et la famille", train_prot_class)
                valid_class = train_prot_class # Current examined class is approved for psiblast 
                os.system('bash scripts/psiblast_script.sh ' + test_prot_name + ' ' + train_prot_name)
            
            previous_class = train_prot_class
        if previous_class != valid_class:  # Check if we are comparing with another family and the previous one wasn't valid
            # print("-------------\n","La famille", previous_class, "a ete ignore","\n-------------")
            cpt += 1
        train_dataset.close
    test_dataset.close
    print("Nombre de famille ignorees : ", cpt)
    np.save("dataset/psiblast_res/comp_gain.npy", cpt)