# Detecting the ultra low dimensionality of real networks

This repository contains four folders:

1. cyclesmap: calculations of chordless cycles in a given network.
2. create_SD: generation of $\mathbb{S}^D$ surrogates of a given network.
3. create_feats: calculation of chordless cycles in a set of surrogates.
4. dimension: detection of optimal dimension of a network. 

The workflow to detect the optimal dimension of a given network is:

1. Calculation of chordless cycles of the network:
    ```sh
    $ ./cyclesmap/cyclesmap network network_features 
    ```
2. Generation of surrogates of the network:
    ```sh
    $ ./create_SD/create_SD.sh network resolution n_poll wsize nrealizations maxD
    ```
3. calculation of the chordless cycles of the surrogates:
    ```sh
    $ ./crete_feats/create_feats.sh SDnets/network SDfeats/network
    ```
4. Detection of optimal dimension using Python 3.x (from dimension folder):
    ```sh
    $ python -c'import dimension; dimension.dimension(SDfeats/network,network_features,["triangles", "squares","pentagons"],maxk)'
    ```
    
The following is a brief description of the codes contained in each folder (more information can be found in the corresponding sh and py files).

# cyclesmap

Contains the program cyclesmap to calculate the number of chordless cycles of length 3,4 and 5 of a given network. 

Parameters:

- The name of a file contaning an edge list of a network.

# create_SD

Contains the script create_SD.sh to generate $\mathbb{S}^D$  surrogates of a given network. 

The script requires an edgelist of the given network and a parameter setting.

Parameters:

- network: name of the network (edgelist file must be located in RealNets folder with edge extension and features file in RealFeats with csv extension)
- resolution: number of surrogates per dimension
- n_poll: number of points used to infer the relation T vs. Beta
- wsize: size of the clustering interval in which to create the surrogates
- nrealizations: number of realizations per random Beta value (default=1)  
- maxD: the script will generate surrogates from D=1 to D=maxD

The resulting surrogates will be placed in SDnets folder. Some folders will be created during the process for calculation purposes.

# create_feats 

Contains the script create_feats.sh to calculate chordless cycles from a set of surrogates.

Parameters:

- surrogates folder: name of the folder containing surrogates (they should be organized by dimension, as obtained with create_SD.sh script)
- results folder: name of the destination folder to store the results

# dimension

This folder contains the dimension.py Python library which provides functions to infer the optimal dimension of a network given a set of features (chordless cycles counts) of its surrogates.  The code requires a folder with surrogate features organized by dimension, as obtained with create_feats.sh script. 

The main function of this library is dimension whith parameters:

- surrogate_set: the name of a folder with surrogate features organized by dimension, as obtained with create_feats.sh
- network_features: the name of a file containing network features obtained with cyclesmap script
- predictors: a set of predictors (a subset of ['triangles', 'squares','pentagons']) 
- maxk: a maximum value of k to explore (maxk) 

This function and returns the optimal dimension for the given network, the value of k and the accuracy for the kNN method. 

----------------------------------------------------

# Publication

Detecting the ultra low dimensionality of real networks
Pedro Almagro, Marian Boguna, M. Angeles Serrano
arXiv:2110.14507,  	
https://doi.org/10.48550/arXiv.2110.14507

----------------------------------------------------

Pedro Almagro, Marián Boguñá, M. Ángeles Serrano

Universitat de Barcelona | Universidad de Sevilla

September 7, 2022

Questions related to code: palmagro@us.es
Questions related to calculation of cycles: marian.serrano@ub.edu



