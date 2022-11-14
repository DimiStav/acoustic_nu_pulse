#=======================================================================================================
# Code to calculate an UHE neutrino induced acoustic pulse at a specific position                      =
#                                                                                                      =
# The current python version has been implemented by:                                                  = 
#                                                                                                      =
#               Dimitris Stavropoulos, NCSR "Demokritos", dstavropoulos@inp.demokritos.gr              =
#               Vasilis Tsourapis, NCSR "Demokritos", tsourapis@inp.demokritos.gr                      =
#                                                                                                      =
# Original codes implemented by the ACoRNE collaboration in MATLAB; Description contained partially in =
# https://doi.org/10.1016/j.astropartphys.2007.08.001 and https://doi.org/10.1016/j.nima.2009.05.009   =
# The original MATLAB codes are found in https://www.hep.shef.ac.uk/research/acorne/index.php          =
#                                                                                                      =
#------------------------------------------------------------------------------------------------------=
#                                                                                                      =
# INPUTS: 1) CORSIKA-IW .dat file                                                                      =
#         2) Log10(shower_E/GeV)                                                                       =
#         3) Observer position [x0,y0,z0] in the shower cylinder's frame of reference                  =
#                                                                                                      =
# OUTPUT: Pressure time-series at the observer's position (plotted in a pop-up window)                 =
#                                                                                                      =
# EXAMPLE: python3 acoustic_nu_pulse.py cg_dEdr_13868_10.95_730.dat 10.95 [400,0,6]                    =
#                                                                                                      =
#          where cg_dEdr_13868_10.95_730.dat is the CORSIKA-IW file                                    =
#                                      10.95 is the Log10(shower_E/GeV)                                =
#                                  [400,0,6] is the [x0,y0,z0] position to calculate the pulse         =
#                                                                                                      =
# Please contact in case of bugs.                                                                      =
#=======================================================================================================

import sys
import numpy as np
import scipy.interpolate as sip
import matplotlib.pyplot as plt
import time

start_time = time.time() #start timer


def histc( X , bins ):

    # Define the bin edges with a vector, where the first element is the left edge of the first bin, 
    # and the last element is the right edge of the last bin
    
    clean_X = np.array(X)[ ( X >= np.array(bins)[0] ) & ( X <= np.array(bins)[-1] ) ]
    
    map_to_bins = np.digitize(clean_X,np.array(bins))
    r = np.zeros(np.array(bins).shape)
    
    for i in map_to_bins: r[i-1] += 1
        
    return [ r , map_to_bins ]


def MCGEn( jpri , lscale , rscale , N , imethod ):

    # MCGEn Table based Monte Carlo Generator
    # Inputs 
    # jpri - 2D-Histogam whose statistics we wish to mimic:  size mxn
    # lscale - bin edges of the rows size: (m+1)*1
    # rscale - bin edges of the columns size: 1*(n+1)
    # n - the number of points to be produced 
    
    # Longitudinal distribution of energy by adding all matrix elements in row (radial bins).
    long = np.sum(jpri,axis=1) + 1e-10*np.sum(jpri) # add 1e-10 to ensure monotonic increase (for the interpolation)
    
    per_sum = np.cumsum( long ) / np.sum( long ) 
    X_values = np.concatenate( ( zero, per_sum ) ) # normalized cumulative distribution
    Y_values = lscale # Longitudinal bin edges
    
    inter_querry = np.random.random( ( int(N) , 1 ) ) # x-coordinates for the interpolated points
  
    lpoints_rand_order = imethod( X_values , Y_values , inter_querry ) # interpolation
    lpoints = sorted( lpoints_rand_order )
 
    rpoints = np.zeros( ( int(N) , 1 ) ) # to store interpolated radial coordinates
    
    Throw = histc( lpoints , lscale )
    nthrowa = Throw[0]
    spos = 0
    
    for i in range( len(nthrowa) - 1 ):

        nthrow = nthrowa[i]
        
        if nthrow > 0:
            
            radial = jpri[i] + 1e-8
            
            R_per_sum = np.cumsum( radial ) / np.sum( radial ) 
            R_X_values = np.concatenate( ( zero , R_per_sum ) )
            R_Y_values = rscale 
            
            R_inter_querry = np.random.random( ( int(nthrow) , 1 ) )
    
            rpoint = imethod( R_X_values , R_Y_values , R_inter_querry )
            
            rpoints[spos:spos+int(nthrow)] = rpoint
            spos = spos + int(nthrow)
    
    lpoints_list = []
    rpoints_list = []
    
    for i in range(int(N)):
        
        lpoints_list.append( lpoints[i][0] )
        rpoints_list.append( rpoints[i][0] )
        
        
    return [ np.array(lpoints_list) , np.array(rpoints_list) ]


def pol2cart(rho, phi, z):

    # convert cylindrical coordinates to cartesian
    
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    
    return np.array( [ np.array(x) , np.array(y) , np.array(z) ] )


def atten( atten_case , dkm , f_axis , c ):

    # ATTEN_FNA calculates attenuation in Sea water.
    # inputs f  (Hz) and dkm (distance in km)
    # outputs attenuation as a fraction of the pressure remaing at the particular frequency.
    # Case 1 Learned's Technique based on Lehtinen et al. w0=2.5e10.
    # Case 2 Polynomial fit to Fig 1 Lehtinen et al
    # Case 3 Neiss and Bertin's Complex attenuation
    # Case 4 Ainslie and McColm J. Acoust. Soc. AM. 103 (3) 1671 1998
    # Case 5 No attenuation
    # Case 6 Francois and Garrison
    # Case 7 ACoRNe (Combined A&Mc and N&B);
    # Case 8 Fischer & Simmons (Note there is a bug in this).

    omega = 2.5e10
    
    if atten_case == 1:

        attenb = (1e4 /  np.log(10) / omega / c )* (2*np.pi* f_axis) ** 2 
        
        y = 10 ** ( -(attenb * dkm)/20 )
        return y , attenb
        
       #attendb= 1e4./log(10)/omega/c.*(2*pi*f).^2;
        #y=10.^-((attendb*dkm)/20);
        
    elif atten_case == 7:
        
        T=15
        S=37
        pH=7.9
        z=2
        w=2*np.pi*f_axis
        s = 1j*w
        w_B = 2*np.pi*0.78e3*np.sqrt(S/35)*np.exp(T/26)       
        w_Mg= 2*np.pi*42e3*np.exp(T/17)
        K_B=0.106e-3*np.exp((pH-8)/0.56)/2/np.pi
        K_Mg= 0.52e-3*(1+T/43)*(S/35)*np.exp(-z/6)/2/np.pi     
        K_W=0.00049e-6*np.exp(-T/27+z/17)/4/np.pi**2
        attendb=K_B*w_B*s/(s+w_B) +  K_Mg*w_Mg*s/(s+w_Mg)+ K_W*w**2
        y=10**( -( (attendb*dkm)/20 ) )
        return y, attendb
    
    else:
        pass


def kernel( points , Do , energy , atten_case , nr , fs , c , m=-999 ):

    # UltraFast Acoustic Integral function
    # Inputs Points [x,y,z] where z is oriented along the axis of the shower,
    # units (m) size  n x 3 where n is typically about 10^6. Note the value of n is very
    # distance dependent. At 100m c 10^8 points are needed to give an ultra-clean pulse.
    # At 10km c 10^4 points are sufficient.
    # Do: is the position of the observer [x0,y0,z0]
   
    # energy: is total energy of the shower 10^energy GeV. i.e. if energy = 11
    # then the total shower energy is 10^20 eV.
   
    # atten; is a flag and is passed to attenfna.
   
    # nr: rotational symmetry is exploited by rotating the shower axially. a
    # value of 100 is typical. (default 1)
    # fs: sampling frequency (default 1MHz)
    # m: mean distance from observer to shower. Calculated if not supplied
    # OUTPUTS
    # p the pulse (sampling rate 1MHz default) Note: zero time is at the mean shower
    # transit time ignoring complex attenuation
    # pw the FFT of the pulse (sampling rate 1MHz default)
    # Exyz is a scaled version of the Velocity Potential
        
    # constants
    PI = np.pi 
    beta = 2e-4 # thermal expansion coefficient
    Cp = 3.8e3 # specific heat capacity in J/(kg*K) 
    
    RR = abs(bin_a) + abs(bin_b) + 1  
    time_bins = np.linspace( bin_a , bin_b , RR )
    time_axis = time_bins / fs
    
    f1 = np.linspace( 0 , abs(bin_a) , abs(bin_a)+1 )
    f2 = np.linspace( -abs(bin_b) , -1 , abs(bin_b) )
    frequency_bins = np.concatenate( ( f1 , f2 ) )
    frequency_axis = ( frequency_bins / RR)  *  fs

    diff_filt = 1j * 2 * PI * frequency_axis * np.fft.fftshift( np.blackman( RR ) ) # A Blackman window is used to smooth the integral and is optional
    
    Exyz = np.zeros(RR) # create velocity potential array
    nmc = np.size(points[0]) # determine size of integration ( points[0] = lpoints )

    Do = [Do]*nmc # to calculate the distance from the observer for each point
    
    # supposing rotational symmetry, the total shower simulation consist of [lpoints,rpoints] rotated to nr positions along phi
    phi = np.linspace( 0 , 360 , 1+nr )
    phi = phi[0:-1] # avoid double counting at zero
    
    # loop nr times over phi
    for ang in phi:
        
        ang = (ang/180.0)*PI # convert to radians
        points_rot =  points.T @ np.array( [ [ np.cos(ang) , np.sin(ang) , 0 ] , [ -np.sin(ang) , np.cos(ang) , 0 ] , [ 0 , 0 , 1 ] ] ) # rotate [lpoints,rpoints]

        d2 = ( points_rot - Do )**2 
        distances = np.sqrt( np.sum( d2 , axis=1 ) ) # points distances from the observer
        
        # Using the last argument, the shower-observer distance can be set manually.
        # The m is constantly set as m=-999 in the arguments; Thus the shower-observer distance
        # is extracted as the average of the points distances from the observer. 
        if m == -999:
            m = sum(distances) / nmc
            
        his = histc( distances/c , time_axis + m/c ) 
        
        Exyz = Exyz + his[0] # unscaled velocity potential
         
    Exyz = ( Exyz / len(phi) ) / nmc # normalize

    # scale and divide by distance
    Exyz_scale = ( ( ( ( Exyz * beta / Cp ) / 4 ) / PI ) * ( 10**(9+energy) * 1.6e-19) ) / ( np.array( [m]*len(time_axis) ) + np.array(time_axis) * c ) 
    
    Iw = np.fft.fft( Exyz_scale ) # FFT

    At = atten( atten_case , m*1e-3 , frequency_axis , c ) # account for the attenuation
    
    Pw = Iw * diff_filt * At[0] 
       
    P = np.real( np.fft.ifft( Pw * fs ) )  # Inverse FFT

    return P, Exyz



# === MAIN CODE ===

# parse arguments
corsika_file = sys.argv[1] # .dat file as generated by CORSIKA-IW
energy = float(sys.argv[2]) # Log10(E/GeV)
obs = sys.argv[3] # [x0,y0,z0], to the frame of reference of the shower cylinder
obs_spl = obs[1:-1].split(',')
Do = [ float(obs_spl[0]) , float(obs_spl[1]) , float(obs_spl[2]) ]

# modify the input CORSIKA file to keep a 200m-long cylinder, also remove the last radial overflow column
cut1000_file = []
cntr = 1
with open( corsika_file  , 'r' ) as file:    
    lines = file.readlines()
    for line in lines:
        if cntr <= 1000:
            keep = line[6:-10]+'\n'
            cut1000_file.append(keep)
            cntr += 1

corsika_file_cut1000 = corsika_file[0:-4] + '_cut1000.dat'

with open( corsika_file_cut1000 , 'w') as file:
    file.writelines(cut1000_file)

corsika_shower = np.genfromtxt( corsika_file_cut1000 , dtype=None)
tsmc_new = corsika_shower


zsc = np.linspace( 9.756 , 19502.244 , 1000 ) # 200m-long cylinder
#zsc = np.linspace( 0.09756 , 195.02244 , 100 ) # comment out in case of a 20m-long cylinder
rsc1 = np.linspace(0.5, 9.5, 10)
rsc2 = np.linspace(15, 105, 10)
rsc = np.concatenate( ( rsc1 , rsc2 ) ) #radial bin centers provided for range [0, 10] with step 1. So, 10 partitions.

#make bin arrays with edges (not centers)
zero = np.array([0])
zsc_new = zsc + 10
lscale = np.concatenate( (zero, zsc_new)  )

radial_add_array = np.concatenate( (0.5*np.ones(10), 5*np.ones(10) ) )
rsc_new = rsc + radial_add_array
rscale = np.concatenate( (zero, rsc_new) )

jpri = tsmc_new
imethod = sip.pchip_interpolate 
N = 1e6 # number of MC points along the shower L-r profile
print("===> Generating MC points according to the deposited energy density...")
points_cylindrical = MCGEn( jpri , lscale , rscale , N , imethod ) # generate points randomly according to the deposited energy density
points = pol2cart( points_cylindrical[1] , np.random.random(int(N))*2*np.pi , points_cylindrical[0] ) # convert to cartesian coordinates
points = np.array(points) * 0.01 # cm to m

# for the binning of histograms
bin_a = -512
bin_b = 511
RR = abs(bin_a) + abs(bin_b) + 1

# constants - kernel
nr = 10  # number of phi angles to rotate points  
c = 1500 # speed of sound in m/s
fs = 1e6 # sampling frequency in Hz
atten_case = 7 # attenuation model as in doi.org/10.1016/j.nima.2009.05.009
print("===> Integration over the generated points along the shower...")
zzk = kernel( points , Do , energy , atten_case , nr , fs , c ) 
pulse = zzk[0]    

# binned time axis
time_bins = np.linspace( bin_a , bin_b , RR ) 
time_axis = time_bins / fs 

print("--- %s seconds ---" % (time.time() - start_time)) # stop timer

plt.grid(axis='both')
plt.plot(time_axis, pulse)
plt.show()



