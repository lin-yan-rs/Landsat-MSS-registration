This repository contains C code for registering Landsat Multispectral Scanner (MSS) L1TP time series to Landsat 8/9 Operational Land Imager (OLI) L1TP reference images.

Prerequisites: 
    (i) Windows Visual Studio 2019 (or compatible versions); 
    (ii) test data included in the repository ("testdata" folder containing two Landdsat 1 MSS L1TP images and two Landsat 9 L1TP images).

How to run the program: 
        (i) create a new project in Windows Visual Studio; 
        (ii) add the provided source files (.c and .h) to the project; 
        (iii) open the 042036_list_all.txt ASCII file in the testdata folder, update the file paths for the four binary images in the "testdata/binary" folder; 
        (iv) open the LSRegistration_main.cpp file, update variables to specify the paths to the test data folder and the 042036_list_all.txt file; 
        (v) run the Program in Visual Studio (Release mode: Press Alt + F5 to run; Debug mode: Press F5 to run). 

The program generates the following outputs:
    (i) registered images: a stack of registered images; 
    (ii) unregistered images: a stack of the original, unregistered images; 
    (iii) registration parameters and accuracy: An ASCII file containing registration parameters for each MSS image and the registration accuracy statistics; 
    (iv) matching results: dense-matching resutls between every two images.

Reference outputs
The results from a successful run of the software, for reference, can be accessed at https://drive.google.com/file/d/1vqpbgicTdABiiypezZ7hKVLUHzljHyJK/view?usp=sharing

If you use this software, please cite the following paper: Yan, L., Roy, D.P. Using Landsat 8 and 9 Operational Land Imager (OLI) data to characterize geometric distortion and improve geometric correction of Landsat Multispectral Scanner (MSS) imagery.
