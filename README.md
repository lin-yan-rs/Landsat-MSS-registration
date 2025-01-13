This repository contains C code for registering Landsat Multispectral Scanner (MSS) time series to Landsat 8/9 Operational Land Imager (OLI) reference images.

Prerequisites
    Windows Visual Studio 2019 (or compatible versions)
    Test data included in the repository ("testdata" folder).

How to run the program
        Create a new project in Windows Visual Studio 2019.
        Add the provided source files (.c and .h) to the project.
        Open the 042036_list_all.txt ASCII file in the testdata folder, update the file paths for the four binary images in the "testdata/binary" folder.
        Open the LSRegistration_main.cpp file, update variables to specify the paths to the test data folder and the 042036_list_all.txt file.
        Run the Program in Visual Studio (Release mode: Press Alt + F5 to run; Debug mode: Press F5 to run).

The program generates the following outputs:
    Registered images: a stack of registered images.
    Unregistered images: a stack of the original, unregistered images.
    Registration parameters and accuracy: An ASCII file containing registration parameters for each MSS image and the registration accuracy statistics.
    Matching results: dense-matching resutls between every two images.

Reference outputs
The results from a successful run of the software, for reference, can be accessed at: https://drive.google.com/file/d/1vqpbgicTdABiiypezZ7hKVLUHzljHyJK/view?usp=sharing

If you use this software, please cite the following paper: Yan, L., Roy, D.P. Using Landsat 8 and 9 Operational Land Imager (OLI) data to characterize geometric distortion and improve geometric correction of Landsat Multispectral Scanner (MSS) imagery.
