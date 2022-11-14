### Code to calculate an UHE neutrino induced acoustic pulse at a specific position                      

The current python version has been implemented by:                                                 
                                                                                                     
> Dimitris Stavropoulos, NCSR "Demokritos", dstavropoulos@inp.demokritos.gr          
> Vasilis Tsourapis, NCSR "Demokritos", tsourapis@inp.demokritos.gr                     
       
Original codes implemented by the ACoRNE collaboration in MATLAB.
Description contained partially in https://doi.org/10.1016/j.astropartphys.2007.08.001 and https://doi.org/10.1016/j.nima.2009.05.009   
The original MATLAB codes are found in https://www.hep.shef.ac.uk/research/acorne/index.php                                                  

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
