### Calculate an UHE neutrino induced underwater acoustic pulse at a specific position                      

The current python version has been implemented by:                                                 
                                                                                                     
> Dimitris Stavropoulos, NCSR "Demokritos", dstavropoulos@inp.demokritos.gr          
> Vasilis Tsourapis, NCSR "Demokritos", tsourapis@inp.demokritos.gr                     
       
Original codes implemented by the ACoRNE collaboration in MATLAB.
Description contained partially in https://doi.org/10.1016/j.astropartphys.2007.08.001 and https://doi.org/10.1016/j.nima.2009.05.009.   
The original MATLAB codes are found in https://www.hep.shef.ac.uk/research/acorne/index.php.                                                  

#### INPUTS: 
- CORSIKA-IW .dat file                                                                      
- Log10(shower_E/GeV)                                                                       
- Observer position [x<sub>0</sub>,y<sub>0</sub>,z<sub>0</sub>] in the shower cylinder's frame of reference                  
                                                                                               
#### OUTPUT: 
- Pressure time-series at the observer's position (plotted in a pop-up window)                 
                                                                                                     
Execute as: 
```
python3  acoustic_nu_pulse.py  cg_dEdr_13868_10.95_730.dat  10.95  [400,0,6]
```
where:
- `cg_dEdr_13868_10.95_730.dat` is the CORSIKA-IW file
- `10.95` is the Log10(shower_E/GeV)
- `[400,0,6]` is the [x<sub>0</sub>,y<sub>0</sub>,z<sub>0</sub>] position to calculate the pulse    

#### Additional Notes
1. The original ACoRNE machinery also provides the option to use parameterised showers as input (several parameterization models are included); In this work only simulated showers with the CORSIKA-IW can be parsed as input.
2. The underwater sound attenuation model followed in this work is the same as in https://doi.org/10.1016/j.nima.2009.05.009.
3. A convention for a 200m long cylinder for the shower simulation is followed. The input CORSIKA-IW .dat file is parsed as generated, and it is being internally modified during the execution. As a result, an additional .dat file is created that meets this convention of the 200m long cylinder.


