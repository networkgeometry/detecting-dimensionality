# Python library to infer the dimension of a network given a set of feature maps of its surrogates (generated with create_feats.sh script). 
# The code requires a folder with surrogates organized by dimension, as obtained with create_feats.sh script. 
#
# This main function of this library is dimension(surrogate_set,network_features,predictors,maxk) that recieves a set of feature maps of surrogates of a network, the feature map of the network, a set of predictors (a subset of ['triangles', 'squares','pentagons']) and a maximum value of k to explore (maxk) and returns the infered dimension for the given network, the value of k and the accuracy for the kNN method. surrogate_set indicates the name of te folder with surrogates organized by dimension, as obtained with create_feats.sh script. network_features indicates the name of a file containing the featuremaps of the network obtained with cyclesmap fortran script.

import math
import os
import pandas as pd
import numpy as np
import sklearn
from sklearn.neighbors import KNeighborsClassifier
from sklearn import metrics

def evalua(X,y,k):

    reg = KNeighborsClassifier(weights="distance",n_neighbors=k)
    predicteds = []
    reals=[]
    X = X.set_index(np.array(list(range(0,len(X)))))
    y = y.set_index(np.array(list(range(0,len(y)))))
    for n in range(0,len(X)):
        example_x = X.iloc[n]
        example_y = y.iloc[n]
        train_x = X.drop(n)
        train_y = y.drop(n)
        reg.fit(train_x,train_y)
        reals.append(example_y)
        predicteds.append(reg.predict([example_x]))
    res = metrics.accuracy_score(reals,predicteds)
    return res

def itera(data,reg,predictors,outcome,nc):
    X, y = data[predictors],data[outcome]

    n_features = X.shape[1]
    m=0
    bestk=-1
    accs=[]
    for k in range(1,nc+1):
        res = evalua(X,y,k)
        accs.append(res)
        if res >= m:
            bestk = k
            bestacc = res 
    return bestk,bestacc,accs

def read_dimensions(surrogate_set,network_features):
    d = surrogate_set
    df = pd.DataFrame(columns=['D','B','iter','triangles', 'squares', 'squaresP','pentagons', 'pentagonsP',"nnodes","nedges"])
    idx=0
    for root, dirs, files in os.walk(d):
        for name in files:
            a = pd.read_csv(os.path.join(root, name),names=['triangles', 'squares', 'squaresP','pentagons', 'pentagonsP',"nnodes","nedges","index"])
            a['D'] = float(name.split("_")[-3])
            a['iter'] = float(name.split("_")[-1].split(".")[0])
            a['B'] = float(name.split("_")[-2])/100
            df = df.append(a)
            idx = idx+1
    try:    
        g = pd.read_csv(network_features,names=['triangles', 'squares', 'pentagons',"nnodes","nedges"])
    except:
        g = False
    df = df.astype(float)
    return df,g

def dimension(surrogate_set,network_features,predictors,maxk):
    dataset = surrogate_set
    try:
        dfSD,g = read_dimensions(surrogate_set,network_features)
        k,acc,accs = itera(dfSD,KNeighborsClassifier(weights="distance"),predictors,["D"],maxk)   
        model = KNeighborsClassifier(weights="distance",n_neighbors=k)
        model.fit(dfSD[predictors],dfSD["D"]) 
        d = model.predict(g[predictors])[0]
    except:
        if(os.path.exists("SDfeats/"+dataset+"/D*")):
            print(dataset+"/D* folder found, assigning dimension 1.")
            d,acc = 1,1        
        else:
            print(dataset+"/D* folder not found or real network features file not found or some value in feats is NaN.")
            d,acc = -1,-1
    print("Dimension detected:"+str(d))
    print("Accuracy:"+str(acc))
    print("Number of neighbors used in the estimation (k):"+str(k))
    return d,acc,k
