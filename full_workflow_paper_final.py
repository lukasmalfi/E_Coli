#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
from ete3 import Tree, faces, TreeStyle, TextFace, CircleFace
import os
import holoviews as hv
from bokeh.plotting import figure, output_file, show
from bokeh.models import ColumnDataSource, Label, LabelSet, Range1d
from bokeh.io import export_png
hv.extension('bokeh')


# In[2]:


#load the gene presence absence matrix from Horesh et al (https://microbiology.figshare.com/articles/dataset/A_comprehensive_and_high-quality_collection_of_E_coli_genomes_and_their_genes/13270073)
f_test=pd.read_csv("F4_complete_presence_absence.csv", sep=",")


# In[3]:


#Drop the row with lsrB_1, not being an actual lsrB 
f_test.drop(labels=[12470],axis=0,inplace=True)


# In[4]:


#Filter out all genes of unknown function, idetified by containing "group"
f_tst_no_grp=f_test[~f_test["Strain"].str.contains("group")]


# In[5]:


#split up the gene names (e.g lsrB_1_1 in ) with _ and use the first part- the actual gene names
for num,gnname in enumerate(f_tst_no_grp["Strain"].str.split("_")):
    f_tst_no_grp.iloc[num,0]=f_tst_no_grp.iloc[num,0].split("_")[0]


# In[6]:


#create dictionary later used to sum up the presence/absence of the same gene
myDict = {}
for i in f_tst_no_grp.columns:
    if i != "Strain":
        myDict.update(dict.fromkeys([i], "sum"))


# In[7]:


#new dataframe, add up the presence of genes with the same short name
df_new = f_tst_no_grp.groupby(f_tst_no_grp['Strain']).aggregate(myDict)

#set maximum value of 1 (gene is present)
df_new_clipped=df_new.clip(upper=1)

#also the lineage information got clipped to 1, so restoring it to the true value
df_new_clipped.loc["Lineage"]=df_new.loc["Lineage"]


# In[8]:


#calculate frequencies based on the new file, unified and w/o lsrB_1, first transpose the matrix
df_new_T=df_new_clipped.T



#calculate the sum of the present genes and divide by the total number of genomes in each lineage
df_sum=df_new_T.groupby(["Lineage"]).sum()

df_nr=df_new_T.groupby(["Lineage"]).count()

df_new_fr=df_sum/df_nr

#transpose back
freq_lsr=df_new_fr.T


# In[9]:


#create list for the gene names and correlations relative to the lsrB gene, nan corresponds to a gene present in all genomes
pearsons=[]
for i in range(len(freq_lsr)): 
    pearsons.append(str(freq_lsr.iloc[i].name)+ " " +str(freq_lsr.loc["lsrB"].astype(float).corr(freq_lsr.iloc[i].astype(float))))


# In[10]:


#convert the list into a dataframe and sort it in descending order
df_pear=pd.DataFrame(columns=["corr"])

df_pear["corr"]=pearsons

df_corr_split=df_pear["corr"].str.split(" ",expand=True)

df_corr_split.columns=["Gene","Pearson"]


#df_corr_split[df_corr_split["Jaccard"]=="nan"]

df_fin_lsrb=df_corr_split.sort_values(["Pearson"],ignore_index=True, na_position="last", ascending=False,)


# In[11]:


#export the dataframe and manually correct the sorting and annotate it via String E.Coli K12 -> Suppl. Table 2


# In[12]:




#set the environment as offscreen, important when working on a remote cluster

os.environ['QT_QPA_PLATFORM']='offscreen'


# In[13]:


#load phylogenetic tree
t = Tree("tree_50.nwk")


# In[14]:


ts = TreeStyle()


# In[15]:


#delete the nodes for which we have no genes
A = t.search_nodes(name="21")[0]
B = t.search_nodes(name="49")[0]
C = t.search_nodes(name="43")[0]

A.delete()
B.delete()
C.delete()


# In[16]:


#select the genes of interest
freq_lsrb=freq_lsr.loc["lsrB"]
freq_lsrg=freq_lsr.loc["lsrG"]
freq_frla=freq_lsr.loc["frlA"]


# In[17]:


#Create dictionary with lineage as key and the percentage as value.
myDict_lsrb = {}
for n,ele in enumerate(freq_lsrb.keys()):
    myDict_lsrb[str(ele)]=np.round(freq_lsrb.values[n],3)
myDict_frla = {}
for n,ele in enumerate(freq_frla.keys()):
    myDict_frla[str(ele)]=np.round(freq_frla.values[n],3)
myDict_lsrg = {}
for n,ele in enumerate(freq_lsrg.keys()):
    myDict_lsrg[str(ele)]=np.round(freq_lsrg.values[n],3)


# In[18]:


#function to create tree with lsrB and frlA
def mylayout_lsrbfrla(node):
    # If node is a leaf, add the nodes name and a its scientific
    # name
         
    #annotate nodes with color and shape depending on which gene and frequency it corresponds to
    if node.is_leaf():
        
        node.img_style["size"] = 9
        node.img_style["shape"] = "square"
        
        for i, name in enumerate(set(node.get_leaf_names())):
            if myDict_lsrb[name]>0.95:
                node.img_style["fgcolor"] = "green"
                
            elif myDict_lsrb[name]<0.05:
                 node.img_style["fgcolor"] = "red"
            else: 
                node.img_style["fgcolor"] = "orange"
        for i, name in enumerate(set(node.get_leaf_names())):
            if myDict_frla[name]>0.95:
                C = CircleFace(radius=6, color="green", style="circle")
              
            elif myDict_frla[name]<0.05:
                C = CircleFace(radius=6, color="red", style="circle")
               
            else: 
                C = CircleFace(radius=6, color="orange", style="circle")
               
            
            faces.add_face_to_node(C, node, 3, position="branch-right")
    
        
        # And place as a float face over the tree
     
    #If node is an internal node
    else:
        # Sets the style of internal nodes
        node.img_style["size"] = 6
        node.img_style["shape"] = "circle"
        node.img_style["fgcolor"] = "#000000" 

ts = TreeStyle()

ts.layout_fn = mylayout_lsrbfrla

ts.show_branch_support = True

ts.title.add_face(TextFace("lsrB/frlA", fsize=15), column=0)

t.render("lsrb_frla.svg", tree_style=ts, w=200, units="mm")


# In[19]:


#function to create tree with lsrB and lsrG
def mylayout_lsrblsrg(node):
    # If node is a leaf, add the nodes name and 

    if node.is_leaf():
        
        #annotate nodes with color and shape depending on which gene and frequency it corresponds to
        node.img_style["size"] = 9
        node.img_style["shape"] = "square"
        
        for i, name in enumerate(set(node.get_leaf_names())):
            if myDict_lsrb[name]>0.95:
                node.img_style["fgcolor"] = "green"
                
            elif myDict_lsrb[name]<0.05:
                 node.img_style["fgcolor"] = "red"
            else: 
                node.img_style["fgcolor"] = "orange"
        for i, name in enumerate(set(node.get_leaf_names())):
            #>95% of genomes have the gene: green, <5%=red, in the middle=yellow
            if myDict_lsrg[name]>0.95:
                C = CircleFace(radius=6, color="green", style="circle")
               
            elif myDict_lsrg[name]<0.05:
                C = CircleFace(radius=6, color="red", style="circle")
               
            else: 
                C = CircleFace(radius=6, color="orange", style="circle")
              
            faces.add_face_to_node(C, node, 3, position="branch-right")
       
        # And place as a float face over the tree
      

        # text faces support multiline. We add a text face
        # with the whole description of each leaf.
       # descFace = faces.TextFace(code2desc[node.name], fsize=10)
       # descFace.margin_top = 10
        #descFace.margin_bottom = 10
        #descFace.border.margin = 1

        # Note that this faces is added in "aligned" mode
        #faces.add_face_to_node(descFace, node, column=0, aligned=True)

        # Sets the style of leaf nodes
        
    #If node is an internal node
    else:
        # Sets the style of internal nodes
        node.img_style["size"] = 6
        node.img_style["shape"] = "circle"
        node.img_style["fgcolor"] = "#000000" 

ts = TreeStyle()

ts.layout_fn = mylayout_lsrblsrg

ts.show_branch_support = True

ts.title.add_face(TextFace("lsrB/lsrG", fsize=15), column=0)

#ts.rotation = 180
#ts.show_scale=False 

t.render("lsrb_lsrg.svg", tree_style=ts, w=200, units="mm")


# In[20]:


#combine and refine the two images manually, add legend -> Figure 3


# In[21]:


#Next, generation of ext. data figure 5


# In[22]:


#select rows with lineages contaning >10% of both frlA and lsrB
with_lsrb_lins_1=freq_lsr.T[freq_lsr.loc["lsrB"]>0.1 ]
with_lsrb_lins=with_lsrb_lins_1[with_lsrb_lins_1["frlA"]>0.1]
#select rows with lineages contaning <10% of both frlA and lsrB
no_lsrb_lins_1=freq_lsr.T[freq_lsr.loc["lsrB"]<=0.1]
no_lsrb_lins=no_lsrb_lins_1[no_lsrb_lins_1["frlA"]<=0.1]


# In[23]:


#lets calculate branch lenghts of tree and try to get together with the other similarity measurements


# In[24]:


#get the names of all lineages with lsrb and frla
lineages_with=with_lsrb_lins.iloc[:,0].index

#and same for without
lineages_no=no_lsrb_lins.iloc[:,0].index


# In[25]:


#get lineages that contain either frlA or lsrB
partly=[]
for lin in freq_lsr.columns:
    if lin not in lineages_with and lin not in lineages_no:
        partly.append(lin)

#dataframe with lineages that contain either lsrb or frla
partly_lsrb_lins=freq_lsr.T.loc[partly]


# In[26]:


#now create groups as specified in ext data figure 5 and compare them to eachother


# In[27]:


#get values for the genomic distance

#with vs with
values_with=[]
for r1 in range(len(with_lsrb_lins)-1):

    for r2 in range(r1+1,len(with_lsrb_lins)):
          #minimum of each gene frequency
        mini=np.minimum(with_lsrb_lins.iloc[r1],with_lsrb_lins.iloc[r2])
            #maximum of each gene frequency
        maxi=np.maximum(with_lsrb_lins.iloc[r1],with_lsrb_lins.iloc[r2])
            #calculate ratio/final value for genomic similarity of genes with a known function
        values_with.append(sum(mini)/sum(maxi))

#without vs without
values_no=[]
for r1 in range(len(no_lsrb_lins)-1):

    for r2 in range(r1+1,len(no_lsrb_lins)):
   
        mini=np.minimum(no_lsrb_lins.iloc[r1],no_lsrb_lins.iloc[r2])
        maxi=np.maximum(no_lsrb_lins.iloc[r1],no_lsrb_lins.iloc[r2])
        values_no.append(sum(mini)/sum(maxi))

#with vs without
values_against=[]
for r1 in range(len(no_lsrb_lins)):

    for r2 in range(len(with_lsrb_lins)):
        if r1!=r2:
            mini=np.minimum(no_lsrb_lins.iloc[r1],with_lsrb_lins.iloc[r2])
            maxi=np.maximum(no_lsrb_lins.iloc[r1],with_lsrb_lins.iloc[r2])
            values_against.append(sum(mini)/sum(maxi))

#with either or vs with either or        
values_partly=[]
for r1 in range(len(partly_lsrb_lins)-1):

    for r2 in range(r1+1,len(partly_lsrb_lins)):
    #for c1 in range(len(with_lsrb_lins.columns)):
        mini=np.minimum(partly_lsrb_lins.iloc[r1],partly_lsrb_lins.iloc[r2])
        maxi=np.maximum(partly_lsrb_lins.iloc[r1],partly_lsrb_lins.iloc[r2])
        values_partly.append(sum(mini)/sum(maxi))
        
        
#all against all
values_all=[]
for r1 in range(len(freq_lsr.T)-1):

    for r2 in range(r1+1,len(freq_lsr.T)):
        if r1!=r2:
  
            mini=np.minimum(freq_lsr.T.iloc[r1],freq_lsr.T.iloc[r2])
            maxi=np.maximum(freq_lsr.T.iloc[r1],freq_lsr.T.iloc[r2])
            values_all.append(sum(mini)/sum(maxi))
            


# In[28]:


#get the tree distances as a measure for their phylogenetic distance

#with vs with
dists_with=[]
for r1 in range(len(lineages_with)-1):

    for r2 in range(r1+1,len(lineages_with)):
    #for c1 in range(len(with_lsrb_lins.columns)):
        
        dist=t.get_distance(str(lineages_with[r1]),str(lineages_with[r2]))
        dists_with.append(dist)

dists_no=[]
#without vs without
for r1 in range(len(lineages_no)-1):

    for r2 in range(r1+1,len(lineages_no)):
    #for c1 in range(len(with_lsrb_lins.columns)):
        
        dist=t.get_distance(str(lineages_no[r1]),str(lineages_no[r2]))
        dists_no.append(dist)

#with vs without
dists_against=[]
for r1 in range(len(lineages_no)):

    for r2 in range(len(lineages_with)):
    #for c1 in range(len(with_lsrb_lins.columns)):
        if r1!=r2:
            dist=t.get_distance(str(lineages_no[r1]),str(lineages_with[r2]))
            dists_against.append(dist)

#all against all

dists_maybe=[]
#without vs without
for r1 in range(len(partly)-1):

    for r2 in range(r1+1,len(partly)):
    #for c1 in range(len(with_lsrb_lins.columns)):
        
        dist=t.get_distance(str(partly[r1]),str(partly[r2]))
        dists_maybe.append(dist)
dists_all=[]
for r1 in range(len(freq_lsr.T)-1):

    for r2 in range(r1+1,len(freq_lsr.T)):
        if r1!=r2:
            dist=t.get_distance(str(freq_lsr.columns[r1]),str(freq_lsr.columns[r2]))
            dists_all.append(dist)
            


# In[29]:


#plot the phylogenetic distance of each lineage-pair against the genomic simlarioty, and color them according to the group they belong to
p = figure(plot_width=600, plot_height=600)



p.circle( y=dists_all, x=values_all, radius=0.0005, color='grey')


p.circle( y=dists_with, x=values_with, radius=0.0005, color='blue')


p.circle( y=dists_no, x=values_no, radius=0.0005, color='red')


p.circle( y=dists_against, x=values_against, radius=0.0005, color='black')

p.circle( y=dists_maybe, x=values_partly, radius=0.0005, color='orange')



p.yaxis.axis_label = 'Evolutionary Distance'
p.xaxis.axis_label = 'Functional similarity'

show(p)
#export_png(p,filename="Ext_Fig_5.png")

