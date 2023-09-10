# Hyperspectral Unmixing using K-means and FCLS

## Group Number - 16

### Group Members - 
 (1) Krisha Shah - 21d070039

 (2) Joel Anto Paul - 210070037
 
 (3) Shravani Kode - 210070086

### This projects contain a MATLAB spectral unmixing based on:

- K-means algorithm in combination with Principal Component Analysis for dimensionality reduction

- Fully-constrained least squares (FCLS) using endmembers of data matrix

### Add the MATLAB data file and MATLAB ground truth file to the 'data' directory.

### The dataset can be downloaded from: https://www.ehu.eus/ccwintco/index.php?title=Hyperspectral_Remote_Sensing_Scenes

### User selected parameters:

- number of components to be returned from PCA and number of clusters to be initialised in K-means

- number of endmembers to take in FCLS

###  The number of endmembers in FCLS we took were 6 for Salinas, 8 for SalinasA, 16 for Indian Pines, and 17 for KSC.

### The number of principal components in K-means we took were 16 for Salinas, 6 for SalinasA, 16 for Indian Pines and 13 for KSC.

### The 'kmeans_canberra.mlx' is a MATLAB live script which performs K-means clustering on Salinas, SalinasA, Indian Pines, and  Kennedy Space Center (KSC) using canberra distance.

### The 'kmeans_euclidean.mlx' is a MATLAB live script which performs K-means clustering on Salinas, SalinasA, Indian Pines, and  Kennedy Space Center (KSC) using euclidean distance.

### The 'fcls.mlx' is a MATLAB live script which performs abundance mapping on Salinas, SalinasA, Indian Pines, and  Kennedy Space Center (KSC) using FCLS.

### The results we got can be found in 'results' directory.

### Functions in 'kmeans_xxx.m' :

- Contains ' function [theta,bel,J]=k_means(X,theta) ' : 

    * This function implements the k-means algorithm, which requires as input the number of clusters underlying the data set. The algorithm starts with an initial estimation of the cluster representatives and iteratively tries to move them into regions that are "dense" in data vectors, so that a suitable cost function is minimized. This is achieved by performing (usually) a few passes on the data set. The algorithm terminates when the values of the cluster representatives remain unaltered between two successive iterations.

    * INPUT ARGUMENTS:
    
        + X : l x N matrix, each column of which corresponds to an l-dimensional data vector.
        + theta :   a matrix, whose columns contain the l-dimensional (mean) representatives of the clusters.

    * OUTPUT ARGUMENTS:
        * theta : a matrix, whose columns contain the final estimations of the representatives of the clusters.
        * bel : N-dimensional vector, whose i-th element contains the cluster label for the i-th data vector.
        * J : the value of the cost function (sum of squared Euclidean distances of each data vector from its closest parameter vector) that corresponds to the estimated clustering.

### Functions in 'fcls.m' :

- Contains ' function [M] = hyperConvert2d(M) ' :
    * Converts a 3D HSI cube (m x n x p) to a 2D matrix of points (p X N) where N = mn.

    * Inputs : M - 3D HSI cube (m x n x p)
    * Outputs : M - 2D data matrix (p x N)

- Contains ' function [img] = hyperConvert3d(img, h, w, numBands) ' :
    * Converts a 2D matrix (p x N) to a 3D data cube (m x n x p) where N = m * n.

    * Inputs : M - 2D data matrix (p x N)
    * Outputs : M - 3D data cube (m x n x p)

- Contains ' function [ X ] = hyperFcls( M, U ) ' :
    * performs fully constrained least squares of each pixel in M using the endmember signatures of U.  Fully constrained least squares is least squares with the abundance sum-to-one constraint (ASC) and the abundance nonnegative constraint (ANC).

    * Inputs :
        * M - HSI data matrix (p x N)
        * U - Matrix of endmembers (p x q)
    * Outputs :
        * X - Abundance maps (q x N)

## Refrences

- https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=774644

- http://www.planetary.brown.edu/pdfs/3096.pdf

- https://www.researchgate.net/publication/224242282_Fully_Constrained_Least_Squares_Spectral_Unmixing_by_Simplex_Projection

- https://in.mathworks.com/help/images/ref/countendmembershfc.html

- https://in.mathworks.com/help/images/ref/nfindr.html#mw_c88ca6df-5bfc-4fb2-ba9b-dac2cfc5b6cb_sep_mw_15786f7b-d8ef-430a-b053-ab3a5274323c
