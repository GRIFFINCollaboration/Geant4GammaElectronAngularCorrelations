# Geant4GammaElectronAngularCorrelations
The angular correlations extension was modified to include internal conversion.

To install the Geant4 source files and the detectorSimulations source files, ```rsyncit.sh``` is provided.

1. open ```rsyncit.sh``` with a text editor and replace the paths to the original source directories if necessary

2. run ```chmod u+x rsyncit.sh``` if the file is not executable

3. run ```./rsyncit.sh```

4. run ```make clean``` and then ```make``` in your Geant4 build directory
<br><hr><br>

FOR A PURE SIMULATION:

5. run ```make``` in ```rdecay01```

6. create a new build directory

7. run 
   ``` 
   cmake -DGeant4_DIR=/opt/geant4.10.01.p03-install/lib64/Geant4-10.1.3/ ~/Geant4GammaElectronAngularCorrelations/rdecay01/
   ```
   The above command assumes that your Geant4 install directory is named ```geant4.10.01.p03-install``` like mine was, and that you've checked out this repository into your home directory. You may also have ```cmake3``` installed, as opposed to ```cmake```.

8. run ```make``` in this build directory

*** repeat steps 4 (without ```make clean```), 5, and 8 for every change to the Geant4/pure simulation source code.
<br><hr><br>

FOR A DETECTOR SIMULATION

9. run ```make``` in your ```detectorSimulations_v10``` directory

10. create a new build directory

11. run ```cmake -DGeant4_DIR=/opt/geant4.10.01.p03-install/lib64/Geant4-10.1.3/ ~/detectorSimulations_v10/``` 
    [see notes in step 7].
    
12. run ```make``` in this build directory

*** repeat steps 4 (without ```make clean```), 9, and 12 for every change to the Geant4/detector simulation source code.
