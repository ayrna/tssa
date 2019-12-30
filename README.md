![AYRNA TSSA logo](logo.png)
<!-- TOC depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:1 -->

- [AYRNA TSSA](#ayrna-tssa)
- [Cite AYRNA TSSA](#cite-ayrna-tssa)
- [Structure and algorithms included](#structure-and-algorithms-included)
  - [Algorithms](#algorithms)
  - [Algorithms eval fast](#algorithms-evals-fast)
  - [Reporter](#reporter)
  - [Source code](#source-code)
  - [Time series](#time-series)
- [Brief explanation of how to use](#brief-explanation-of-how-to-use)
  - [Usage](#usage)
  - [Analysis of the results](#analysis-of-the-results)
- [External software](#external-software)
- [References](#references)

<!-- /TOC -->

# AYRNA TSSA
AYRNA TSSA (Time Series Segmentation Algorithms) is a set of algorithms implemented in MATLAB which integrates a wide range of time series segmentation algorithms. It has been developed by ["AYRNA Research Group"](http://www.uco.es/grupos/ayrna/index.php/en) resulting in several publications in international journals and conferences. 

There exist two main objectives which time series segmentation tries to satisfy. On the one hand, this procedure is applied to discover similarities between segments. On the other hand, time series segmentation is also applied with the objective of simplifying the time series, that is, replacing segments by simple model descriptions. AYRNA TSSA contains algorithms which tackle these objectives separately or at the same time.

## Copyright
This software is released under the The GNU General Public License v3.0 licence available at [http://www.gnu.org/licenses/gpl-3.0.html](http://www.gnu.org/licenses/gpl-3.0.html).

Please see the license [file](LICENSE) available in the repository and the headers of the documents for the copyright. 

# Cite AYRNA TSSA
If you use AYRNA TSSA files or algorithms please cite the works which are indicated in the header of each file. For example, if you use GMOTSS algorithm (see section [Structure and algorithms included](#structure-and-algorithms-included)), the header is:

```
%GMOTSS Genetic Multiobjective Time series segmentation [1]
%
%   GMOTSS methods:
%      runAlgorithm               - runs the corresponding algorithm (GMOTSS in [1])
%      saveInformation            - specific information of the algorithm
%      saveAll                    - save all information of the algorithm
%
%   References:
%     [1] A.M. Durán-Rosal, P.A. Gutiérrez, F.J. Martínez-Estudillo and C. Hervás-Martínez.
%         "Simultaneous optimisation of clustering quality and approximation error
%         for time series segmentation", Information Sciences, Vol. 442-443, May, 2018, pp. 186-201.
%         https://doi.org/10.1016/j.ins.2018.02.041
%
%   This file is part of TSSA: https://github.com/ayrna/tssa
%   Original authors: Antonio M. Duran Rosal, Pedro A. Gutierrez Peña
%   Citation: If you use this code, please cite the associated paper [1]
%   Copyright:
%       This software is released under the The GNU General Public License v3.0 licence
%       available at http://www.gnu.org/licenses/gpl-3.0.html 
```

In the case of source code and reporter files (see section [Structure and algorithms included](#structure-and-algorithms-included)), you will find the following kind of header:

```
%   COPYRIGHT
%   This file is part of TSSA: https://github.com/ayrna/tssa
%   Original authors: Antonio M. Duran Rosal, Pedro A. Gutierrez
%   Copyright:
%       This software is released under the The GNU General Public License v3.0 licence
%       available at http://www.gnu.org/licenses/gpl-3.0.html
%   Citation: If you use this code, please cite the following paper:
%     [1] Citation
%
% Function description
```

For more information about the papers and our research group, please visit [Learning and Artificial Neural Networks (AYRNA) website](http://www.uco.es/grupos/ayrna/index.php/en) at [University of Córdoba](http://www.uco.es/) (Spain).


# Structure and algorithms included
The repository is organised as follows:

## Algorithms
The [algorithms](algorithms) folder includes the MATLAB algorithms (and their hybrid versions) for time series segmentation using the number of iterations as a stop criterion.

 - [ACROTSS](algorithms/ACROTSS): Coral reefs optimisation algorithm for the reduction of the number of points of time series [1].
 - [ASCROTSS](algorithms/ASCROTSS): Statistically-driven Coral reefs optimisation algorithm for the reduction of the number of points of time series [1].
 - [ATSS](algorithms/ATSS): Genetic algorithm for the reduction of the number of points of time series [1].
 - [BBePSOTSS](algorithms/BBePSOTSS): Exploitation barebones Particle Swarm Optimisation algorithm for the reduction of the number of points of time series [2].
 - [BHTSS](algorithms/BHTSS): Hybrid algorithm using a likelihood-based segmentation assuming a beta distribution for detecting similarities between segments [3].
 - [CROTSS](algorithms/CROTSS): Coral reefs optimisation algorithm for detecting similarities between segments [4].
 - [DBBePSOTSS](algorithms/DBBePSOTSS): Dynamic Exploitation barebones Particle Swarm Optimisation algorithm for the reduction of the number of points of time series [2].
 - [EvolTSS](algorithms/EvolTSS): Genetic algorithm for detecting similarities between segments [5].
 - [GMOTSS](algorithms/GMOTSS): Multiobjective evolutionary algorithm for the reduction of the number of points and for detecting similarities between segments in the same algorithm [6].
 - [NHTSS](algorithms/NHTSS): Hybrid algorithm using a likelihood-based segmentation assuming a normal distribution for detecting similarities between segments [7].
 - [PSOTSS](algorithms/PSOTSS): Exploitation barebones Particle Swarm Optimisation algorithm for the reduction of the number of points of time series [2].
 - [TRADTSS](algorithms/TRADTSS): Traditional algorithms for the reduction of the number of points of time series extracted from [8].

## Algorithms evals fast
The [algorithms_evals_fast](algorithms_evals_fast) folder includes a reimplemented version of some previous algorithms. The algorithms use the number of evaluations as a stop criterion and a new function to evaluate the solution, which significantly decreases the computational time. Also, this folder includes two new algorithms and their hybrid versions:

 - [AMCROTSS](algorithms_evals_fast/AMCROTSS): Memetic Coral reefs optimisation algorithm for the reduction of the number of points of time series [9].
 - [WBBePSOTSS](algorithms_evals_fast/WBBePSOTSS): Weighted Exploitation barebones Particle Swarm Optimisation algorithm for the reduction of the number of points of time series [10].
 
## Reporter
The [reporter](reporter) folder includes the functions to generate the final reports at the end of the execution of the algorithms. It also consists of two externals tools with their corresponding license of use (see section [External software](#external-software)).

## Source code
The [source_code](source_code) folder saves all the functions needed to execute the different algorithms available in the repository. It also includes external software with its corresponding license to use (see section [External software](#external-software)).

## Time series
In the [time_series](time_series) folder you must include your own time series databases. We provided two examples of time series, which are:
- [Mallat](time_series/MALLAT_.txt): is extracted from [The UCR Time Series Classification Archive](https://www.cs.ucr.edu/~eamonn/time_series_data/).
- [Donoho-Johnstone](time_series/Donoho-Johnstone.txt): is extracted from a [benchmark repository](https://sites.google.com/site/icdmmdl/).


# Brief explanation of how to use
For each algorithm, we have three files:
- "Algorithm.m": this file defines the algorithm.
- "masterSaveAll.m": this file includes a function which saves a summary of a set of runs given the stochastic nature of the algorithms implemented in AYRNA TSS.
- "masterExperimenter.m": this script runs the algorithm given the setting of parameters.

## Usage
1. Click on the selected algorithm.
* For the first time, you must create a folder called "reports"
2. Open masterExperimenter.m file.
3. Configure masterExperimenter.m by setting the experimental parameters. (Note that you should change the files array with your time series).
4. Run masterExperimenter.m file.

## Analysis of the results
1. Once you execute masterExperimenter.m script, a folder (called *time_folder* in the following steps) with the current date and hour is created in the report/ folder.
2. For each run, the algorithm creates a folder (*run_folder* in the following step) inside the *time_folder* with the name equal to the number of the execution.
3. Inside each *run_folder*, and depending on the selected algorithm, we can found different files such as the plotted time series, information of segments, fitness, etc. for the execution.
4. When all executions are finished, masterSaveAll.m is executed and it creates two files in *time_folder*, which are:
- Summary information saved in resultsMultipleRunnings.csv file.
- Summary models saved in informationMultipleRunnings.mat file.

# External software
AYRNA TSSA makes use of the following external software implementations:
- [kmeans](source_code/kmeans): K-means algorithm that we have modified to be available with our TSS algorithms. Please see the [license](source_code/kmeans/license.txt) file. 
- [export_fig](reporter/external_tools/export_fig): Implementation to generate pdf outputs of graphics. Please see the [license](reporter/external_tools/export_fig/license.txt) file.
- [plot2svg](reporter/external_tools/plot2svg): Implementation to generate svg outputs of graphics. Please see the [license](reporter/external_tools/plot2svg/license.txt) file.
- [pdist2](source_code/pdist2.m): It calculates the distance between a set of vectors. Please see the copyright available in the header of the [file](source_code/pdist2.m).

# References
- [1] A.M. Durán-Rosal, P.A. Gutiérrez, S. Salcedo-Sanz and C. Hervás-Martínez. "A statistically-driven Coral Reef Optimization algorithm for optimal size reduction of time series", Applied Soft Computing, Vol. 63. 2018, pp. 139-153. (https://doi.org/10.1016/j.asoc.2017.11.037)
- [2] A.M. Durán-Rosal, P.A. Gutiérrez, Á. Carmona-Poyato and C. Hervás-Martínez. "A hybrid dynamic exploitation barebones particle swarm optimisation algorithm for time series segmentation", Neurocomputing, Vol. 353, August, 2019, pp. 45-55. (https://doi.org/10.1016/j.neucom.2018.05.129)
- [3] A.M. Durán-Rosal, J.C. Fernández, P.A. Gutiérrez and C. Hervás-Martínez. "Detection and prediction of segments containing extreme significant wave heights", Ocean Engineering, Vol. 142, September, 2017, pp. 268-279.(https://doi.org/10.1016/j.oceaneng.2017.07.009)
- [4]  A.M. Durán-Rosal, D. Guijo-Rubio, P.A. Gutiérrez, S. Salcedo-Sanz and C. Hervás-Martínez. "A coral reef optimization algorithm for wave height time series segmentation problems". International Work-Conference on Artificial and Natural Neural Networks (IWANN2017). 14th-16th June. 2017. Cádiz (Spain). LNCS, vol. 10305. pp. 673-684. (https://doi.org/10.1007/978-3-319-59153-7_58)
- [5] M. Pérez-Ortiz, A.M. Durán-Rosal, P.A. Gutiérrez, J.Sánchez-Monedero, A.Nikolaou, F.Fernández-Navarro, C.Hervás-Martínez. "On the use of evolutionary time series analysis for segmenting paleoclimate data", Neurocomputing, Vol. 326-327, January, 2019, pp. 3-14 (https://doi.org/10.1016/j.neucom.2016.11.101)
- [6] A.M. Durán-Rosal, P.A. Gutiérrez, F.J. Martínez-Estudillo and C. Hervás-Martínez. "Simultaneous optimisation of clustering quality and approximation error for time series segmentation", Information Sciences, Vol. 442-443, May, 2018, pp. 186-201. (https://doi.org/10.1016/j.ins.2018.02.041)
- [7] A.M. Durán-Rosal, M. de la Paz Marín, P.A. Gutiérrez and C. Hervás-Martínez. "Identifying market behaviours using European Stock Index time series by a hybrid segmentation algorithm", Neural Processing Letters, Vol. 46, December, 2017, pp. 767–790. (https://doi.org/10.1007/s11063-017-9592-8)
- [8] E. Keogh, S. Chu, D. Hart and M. Pazzani. "Segmenting time series: A survey and novel approach", In Data mining in time series databases, 2004, pp.1-21. (https://doi.org/10.1142/9789812565402_0001)
- [9] A.M. Durán-Rosal, P.A. Gutiérrez, S. Salcedo-Sanz and C. Hervás-Martínez. "Dynamical Memetization in Coral Reef Optimization Algorithms for Optimal Time Series Approximation", Progress in Artificial Intelligence, Vol. 8, June, 2019, pp. 253-262. (https://doi.org/10.1007/s13748-019-00176-0)
- [10] A.M. Durán-Rosal, D. Guijo-Rubio, P.A. Gutiérrez and C. Hervás-Martínez. "Hybrid Weighted Barebones Exploiting Particle Swarm Optimization Algorithm for Time Series Representation". BIOMA2018. 16th-18th May. 2018. Paris (France). LNCS, vol. 10835. pp. 126-137. (https://doi.org/10.1007/978-3-319-91641-5_11)	
