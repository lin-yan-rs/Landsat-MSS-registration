This repository contains C code for registering Landsat Multispectral Scanner (MSS) L1TP time series to Landsat 8/9 Operational Land Imager (OLI) L1TP reference images.

Features: 
    (i) many-to-many matching strategy to better ensure mutual registration of MSS-MSS images and MSS-OLI images; 
    (ii) radial-basis-function (RBF) to accommodate geometric distortions in MSS images; 
    (iii) least-squares adjustment to provide optimal solution to the many-to-many matchings;
    (iv) expert implementation of the least-square matching algorithm providing high matching accuracy (better than 0.01 pixels)

Prerequisites: 
    (i) Windows Visual Studio 2019 (or compatible versions); 
    (ii) MSS and OLI images saved as separate binary files as those in the test data included in the repository (containing two Landdsat 1 MSS L1TP images and two Landsat 9 L1TP images).

How to run the program with the test data: 
    (i) create a new project in Windows Visual Studio; 
    (ii) add the source files (.cpp and .h) to the project; 
    (iii) open the "042036_list_all.txt" ASCII file in the "testdata" folder, update the file paths for the four binary images in the "testdata/binary" folder; 
    (iv) open the "LSRegistration_main.cpp" file, update variables to specify the paths to the test data folder and the "042036_list_all.txt" file; 
    (v) run the program in Visual Studio (Release mode: Press Alt + F5 to run; Debug mode: Press F5 to run). 

The program generates the following outputs:
    (i) a stack of registered images; 
    (ii) a stack of the original unregistered images; 
    (iii) an ASCII file containing registration parameters for each MSS image and the registration accuracy statistics; 
    (iv) intermediate results of dense least-squares matching between every pair of MSS-MSS and MSS-OLI images (in "matchings" folder).

Reference outputs: 
the results from a successful run of the software, for reference, can be accessed at https://drive.google.com/file/d/1vqpbgicTdABiiypezZ7hKVLUHzljHyJK/view?usp=sharing

Please cite the following paper for software usage: Yan, L., Roy, D.P. Using Landsat 8 and 9 Operational Land Imager (OLI) data to characterize geometric distortion and improve geometric correction of Landsat Multispectral Scanner (MSS) imagery. Remote Sensing of Environment. In Review.
