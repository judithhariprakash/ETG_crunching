
# coding: utf-8

# In[1]:


cd /mnt/sdb/ETG/Data_Toolwise/JEME_data/


# In[2]:


import collections
import pandas as pd
import matplotlib.pyplot as plt
import itertools 


# # Each line of the file has three columns:
# -  Enhancer coordinate
# -  Target gene(s) information. Genes sharing the TSS are delimited by '%'
# -  Confidence score for interaction
# 
# A single enhancer can be on multiple lines as it can interact with multiple TSS

# In[3]:


def read_lasso_fn(fn):
    """
    A function that reads the flat format lasso file.
    It returns a dictionary whose keys are enhancer coordinates
    and value a list of lists which contain gene names.
    """
    etp = {}
    unique_enhancers = []
    unique_genes = []
    total_etp = {}
    with open(fn) as handle:
        for line in handle:
            columns = line.rstrip("\n").split(",")
            enhancer = columns[0]
            genes = [x.split('$')[1] for x in columns[1].split('%')]
            if enhancer not in etp:
                etp[enhancer] = []
            etp[enhancer].extend(genes)
    return etp

def make_dict_barplot(input_dict, color='purple', figsize=(7,5)):
    """
    Takes a dictionary and makes a barplot by using keys for
    x axis and values are for y axis
    Returns an axis object which can be used to manipulate the
    figure.
    """
    keys = sorted(input_dict.keys())
    values = [input_dict[x] for x in keys]
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    ax.bar(keys, values, color=color)
    ax.set_xticks(keys)
    return ax


# In[4]:


files = [28, 117, 86, 17, 126, 112, 94]
names = {28:"HMEC_v", 117:"HMEC", 86:"Fetal_lung", 17: "IMR90", 126:"Lung_fibroblast", 112: "A549", 94: "Adult_lung" }
enhancers_for_all = {}
number_table = {}
print ("Unique_enhancers in JEME_prediction")
for i in files:
    et_pairs = read_lasso_fn("fantom5_lasso/fantom5_lasso.%d.csv" % i )
    print (str(names[i])+": "+ str(len(et_pairs)))
    enhancers_for_all[names[i]] = [key for key in et_pairs.keys()]
    number_table[names[i]] = len(et_pairs)


# In[5]:


HMEC_v = set(enhancers_for_all['HMEC_v'])
HMEC= set(enhancers_for_all['HMEC'])
Fetal_lung= set(enhancers_for_all['Fetal_lung'])
IMR90= set(enhancers_for_all['IMR90'])
Lung_fibroblast= set(enhancers_for_all['Lung_fibroblast'])
A549= set(enhancers_for_all['A549'])
Adult_lung= set(enhancers_for_all['Adult_lung'])
union_enhancers = set.union(HMEC_v, HMEC, Fetal_lung, IMR90, Lung_fibroblast, A549, Adult_lung)
number_table.update({'union': len(union_enhancers)})
intersection_enhancers = set.intersection(HMEC_v, HMEC, Fetal_lung, IMR90, Lung_fibroblast, A549, Adult_lung)
number_table.update({'intersection': len(intersection_enhancers)})


# In[7]:


ax = make_dict_barplot(number_table, figsize=(12, 9))
ax.set_xlabel('No. of enhancers', fontsize=15)
ax.set_ylabel('Cell line', fontsize=15)
plt.title("Unique enhancers predicted by JEME", fontsize=15)
plt.show()
#plt.savefig("/media/judith/MyStorage/ETG/Result/Sample_specific_unique_enhancers.jpeg")
