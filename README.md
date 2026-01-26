# INV-2022-SAEZ
Matlab code for research published in Statistical Normalization for a Guided Clustering Type-2 Fuzzy System for WSN

A. J. Yuste-Delgado, J. C. Cuevas-Martínez and A. Triviño-Cabrera, "Statistical Normalization for a Guided Clustering Type-2 Fuzzy System for WSN," in IEEE Sensors Journal, vol. 22, no. 6, pp. 6187-6195, 15 March15, 2022, doi: [10.1109/JSEN.2022.3150066](https://doi.org/10.1109/JSEN.2022.3150066).
keywords: {Wireless sensor networks;Fuzzy systems;Clustering algorithms;Sensors;Input variables;Clustering methods;Base stations;Wireless sensor networks;clustering;dynamic;normalization;Type-2 fuzzy systems},

# Abstract:
One of the main concerns in Wireless Sensor Networks is the efficient energy management of the nodes. Hierarchical techniques such as clustering have been developed in an effort to solve this problem. In this paper we present a smart evolution of a distributed clustering method that uses a turn-based scheduling cluster head selection process based on an interval Type-2 fuzzy system. The method we propose offers four main improvements. First, the setup process guided by the Base Station is adapted to tune the skip parameter during the network lifetime, which controls how many rounds the clusters are not updated. Second, the normalization of the fuzzy system input variables is carefully performed based on a statistical analysis to reduce the effect of fluctuations in edge values. Third, the value of the coefficient applied to the output of the inner Type-2 fuzzy system is updated to balance the number of cluster heads at early stages. Finally, only the strongest candidate nodes, those with the highest probability, are selected to become cluster heads. The proposed design and scheduling aim to achieve low-energy processing in the nodes. When our proposed techniques are applied, they give better results compared with other similar approaches.

# How to run
File: saez_best_v3_variables_de_entrada_simples.m
It is the main file with the proposed algorithm.

## INPUTS:

* x: Vector with the x-coordinates of the sensors
* y: Vector with the y-coordinates of the sensors
* x_max: Maximum x-coordinate
* y_max: Maximum y-coordinate
* n: Number of nodes in the experiment
* stop: Number of deaths required to reach the end, expressed as a percentage of n
* p: Variable p from the LEACH method
* EBx: x-coordinate of the base station
* EBy: y-coordinate of the base station
* Eo: initial energy of nodes
* ETX,ERX,Efs,Eamp,EDA: parameters of energy model 
* packetSize,controlPacketSize: parameters of communication packets
* Eminima: Minimum energy before a node is considered dead

## OUTPUTS:
* FND: Firts node dies
* HND: Halfnode die
* LND: Last node die 
* noCH: number of rounds in which no CHs were chosen
* noCHFND: number of rounds in which no CHs were chosen before FND
* noCHHND: number of rounds in which no CHs were chosen before HND

## FILE: matriz_saez_v3.mat
Arrays with the sampled values ​​of the fuzzy interval type 2.

* FILE: variables_menos_x_y_ebx_eby.may
A* ll the variables of the energy and communications model necessary to launch experiments are in place; only the position of the sensors and the base station are missing.

## FILE: escenario_centro.m
Coordinates with the location of the sensors.

## FILE: pruebas_escenarios_centro.m
File that calculates the algorithm's output for 15 different scenarios, with the base station in the center.
