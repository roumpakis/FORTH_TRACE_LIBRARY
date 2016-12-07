# FORTH-TRACE library

# Table of contents
1. [Introduction](#introduction)
2. [Modules Description](#modules)
3. [Datasets](#datasets)
4. [Useful methods](#useful)
5. [Installation Instructions](#execution)
6. [Licence](#licence)

## Introduction <a name="introduction"></a>
FORTH-TRACE library, as part of the FORTH-TRACE benchmark framework, is a Java implemented software tool that covers the entire feature-level fusion chain applied in Human Activity Recognition (HAR) applications.
Feature selection is considered an important step of filtering redundant information, prior to applying any classification technique, especially for data streams characterized by high correlation.  

The FORTH-TRACE library consists of the following modules.

*  **Acquisition**
*  **Segmentation**
*  **Feature Extraction**
*  **Feature Selection**

## Modules Description <a name="modules"></a>
* **Acquisition:** Collect data from wearable sensors and convert them into the appropriate data structures. This module supports multiple types of HAR datasets.
* **Segmentation:** Split data into sliding windows of observations. The key parameters of the Segmentation component is the size _w_ pf the sliding windows, as well as the overlapping factor _a_.
* **Feature Extraction:** Extract _M_ features for each data segment of each sensor modality. Currently, we support 11 types of statistical features and the Pairwise Correlation among sensor channels, listed below. 
    *  Mean
    *  Median
    *  Standard Deviation
    *  Variance
    *  Root Mean Square
    *  Skewness
    *  Kurtosis
    *  Interquartile Range
    *  Mean Crossing Rate
    *  Zero Crossing Rate
    *  Spectral Entropy
    *  Pairwise Correlation
    
* **Feature Selection:** Contains a pool of available Feature Selection Algorithms (FSA) in order to calculate the 
_R (R <<M)_ representative features for each segment. Currently 4 types of FSA are implemented. 
    * Supervised (Sequential Floating Forward Selection, SFFS)
    * Unsupervised (Feature Selection based on Feature Similarity, FSSA)
    * Ranking-based (Relief-F) 
    * Unsupervised graph-based (Graph Clustering with Node Centrality, GCNC)

## Datasets <a name="datasets"></a>
The FORTH-TRACE library is designed to support any given HAR dataset.
Our implementation is tested using two HAR datasets.
* HAPT dataset -- single location - multiple modalities (Available in http://archive.ics.uci.edu/ml/datasets/Smartphone-Based+Recognition+of+Human+Activities+and+Postural+Transitions)
* FORTH-TRACE dataset -- multiple locations - multiple modalities (Available at https://github.com/spl-icsforth/FORTH_TRACE_DATASET)

## Useful Methods <a name="useful"></a>
Besides the methods used that implement the Machine Learning pipeline, this project contains methods to write the feature matrixes and the experiment statistics into .csv for post-processing evaluation purposes.
These methods are located in the IO package of the project. It is recommended to save the .csv files in the files/ directory of the project.

Furthermore, the Main.java file contains full experiment implementation for the following scenarios.

* HAPT dataset
* FORTH-TRACE dataset (collective scenario -- all devices of 10 random participants)
* FORTH-TRACE dataset (single-device scenario -- for each separate device of 10 random participants) 

## Installation Instructions <a name="execution"></a>
1. Download the project source files.
2. Import the project to Eclipse or other IDE.
3. Add the following libraries into the build path (located in the libs/ directory of the project):
    * commons-math3-3.5.jar
    * javaml-0.1.7.jar
    * weka.jar
4. Specify the path of your dataset files in the Main.java file of the project.


## License <a name="licence"></a>
Use of this source code in publications must be acknowledged by referencing the following publications:

* Katerina Karagiannaki, Athanasia Panousopoulou, Panagiotis Tsakalides. A Benchmark Study on Feature Selection for Human Activity Recognition. ACM International Joint Conference on Pervasive and Ubiquitous Computing (UbiComp), ACM, 2016.