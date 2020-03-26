import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as si
import Experiments as Ex


fc = 76e9
c = 3e8
lamb = c/fc
tm = 256e-6
bw = 1e9
k = bw/tm
fs = 2e6   #sample rate after stretch proccessing
fss = 2e9  #sample rate to generate waveforms
t = np.arange(0, tm, 1/fss)
TxPower = 1e-2 #/watts  10dBm
txGain = 17 #change to function of angle later maybe

rxGain = np.loadtxt('Gain from datasheet.csv',delimiter=',')
rxGain_H= rxGain[:,0:2]
rxGain_H[:,1] = 10**(rxGain_H[:,1]/10)      #convert from db to get non negative values
rxGain_V = rxGain[:,::2]
rxGain_V[:,1] = 10**(rxGain_V[:,1]/10)

n1 = 1 #refractive indicies for fresenel
n2 = 2.5 #concrete floor https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=504907

def Calc_Gain(theta):
    '''currently treat horozontal as fixed 15'''
    idx = np.abs(rxGain_V[:,0] - np.degrees(theta)).argmin()
    gain = 10**(15/10) * rxGain_V[idx,1]
    return gain



def add_zeros(array):
    '''Add zero values before and after data to create a 
        rectangular pulse giving sinc functions after fft'''

    a = np.zeros((1024,len(array.T)))
    b = np.append(a,array,0)
    new = np.append(b,a,0)          #for range sidelobes

    a = np.zeros((len(new),100))
    b = np.append(a,new,1)
    new = np.append(b,a,1)          #for angle sidelobes

    return new


def RAngCalc(data,sweep_slope,fs):
    '''create Range angle map for a given data set'''

    Angle1 = np.linspace(-1,1,len(data.T),endpoint=True)
    Angle1 = np.degrees(np.arcsin(Angle1))                  #arcsin as the array factor is proportional to the sin

    Freq = np.linspace(0,fs,len(data),endpoint=True)   
    Range = (Freq*3e8) / (2*sweep_slope)

    R = np.fft.fft(data,None,0)

    A = np.fft.fft(data,None,1)
    A = np.fft.fftshift(A,1)

    ReflecMap = np.fft.fft(R,None,1)
    ReflecMap = np.fft.fftshift(ReflecMap,1)        #create range and angle map for input data

    
    interp_ratio = 1#20                                #ratio of initial array elements vs interpolated elements
    
    Angle = np.linspace(-1,1,len(data.T)*interp_ratio,endpoint=True) #create new angle vector to plot against
    Angle = np.degrees(np.arcsin(Angle))                    #to go from psi to angle representations
    
    interp_mesh = np.array(np.meshgrid(Range, Angle, indexing='ij'))
    interp_points = np.rollaxis(interp_mesh, 0, 3).reshape((-1, 2)) #create list of cooridnates to interpolate at

    A = si.interpn((Range,Angle1),A,interp_points)
    A = np.reshape(A,(len(Range), len(Angle)))

    ReflecMap = si.interpn((Range,Angle1),ReflecMap,interp_points)
    ReflecMap = np.reshape(ReflecMap,(len(Range), len(Angle)))
    
    return [Range, Angle, R, A, ReflecMap]

def delayedsig(dist):
    '''delays the initial signal by time corresponding to 
       distance signal travelled'''  

    delay = dist/c
    delsig = np.exp(np.pi*1j*(fc*np.subtract(t,delay)+k*(np.subtract(t,delay)**2)))
    
    return delsig


def data_set():

    num = int(input('Pick experiment #: '))
    if num == 1:
        data,Reflector = Ex.fan_on()
        data = add_zeros(data)

    elif num == 2:
        data,Reflector = Ex.fan_off()
        data = add_zeros(data)

    elif num == 3:
        data,Reflector = Ex.corner2()
        data = add_zeros(data)

    elif num == 4:
        data,Reflector = Ex.big1()
        data = add_zeros(data)

    elif num == 5:
        data,Reflector = Ex.big2()
        data = add_zeros(data)

    elif num == 6:
        data, Reflector = Ex.test()
        data = add_zeros(data)

    else:
        data,Reflector = Ex.fan_on()
        data = add_zeros(data)

    return data, Reflector


''' 
Geometric 3-D Simulation

Signals are delayed and has contributions from
reflected rays
'''

def SimData(array_size,multi=0):

    RxArray = [[(i-(array_size-1)/2)*lamb/2,0,0.98] for i in range(array_size)]   #create Uniform linear spaced array

    Temp = np.zeros((len(t),len(RxArray)),dtype=complex)                #to store data before mixing/pulse compression?
    RxSig = np.zeros((int(len(t)/(fss/fs)),len(RxArray)),dtype=complex) #smaller array to store post compression signal


    for Obj in Reflector:  #loop through all objects and array elements

        numerator =  (TxPower * txGain * Obj.rcs * (lamb**2)) / ((4*np.pi)**3)   #common numerator used for determining power           
        for i in range(len(RxArray)):
       
            Y1 = (RxArray[i][1] * Obj.z + Obj.y * RxArray[i][2]) / (Obj.z + RxArray[i][2])              #find x,y,z position of reflection point could move outside as doesnt change much for each array element
            X1 = (RxArray[i][0] + Obj.x) /2
            Z1 = 0

            dist = np.sqrt( (RxArray[i][0]-Obj.x)**2 + (RxArray[i][1]-Obj.y)**2 + (RxArray[i][2]-Obj.z)**2 )       #Direct Ray distance
            dist1 = np.sqrt((RxArray[i][0] - X1)**2 + (RxArray[i][1] - Y1)**2 + (RxArray[i][2] - Z1)**2)           #array to ground distance
            dist2 = np.sqrt((Obj.x - X1)**2 + (Obj.y - Y1)**2 + (Obj.z - Z1)**2)                                   #ground to object distance
            
            incidence = np.abs(np.arctan((Y1-Obj.y)/Obj.z))                                                                    #angle of incidence wrt the normal
            #Ref = np.abs((n1*np.sqrt(1-(n1*np.sin(incidence)/n2)**2) - n2*np.cos(incidence))/(n1*np.sqrt(1-(n1*np.sin(incidence)/n2)**2) + n2*np.cos(incidence)))**2 #assuming TE polarised need to check actually is s i think because from radar stats horo and vertical beamwidths
            Ref = np.abs((n1*np.cos(incidence) - n2*np.sqrt(1-(n1*np.sin(incidence)/n2)**2) )/(n1*np.cos(incidence) + n2*np.sqrt(1-(n1*np.sin(incidence)/n2)**2) ))**2 #TM     

            theta = np.arctan(Obj.z/(Y1 - Obj.y))
            Gain = Calc_Gain(theta)
            LOS_Gain = Calc_Gain(0)
            if multi == 0:
                Temp[:,i] += (delayedsig(2*dist) * np.sqrt((LOS_Gain * numerator) / (dist**4))     )                                  #sum each rays signal at each receiver

            else:
                Temp[:,i] += (delayedsig(2*dist) * np.sqrt((LOS_Gain * numerator)  / (dist**4))
                +( delayedsig(2*(dist1+dist2)) * np.sqrt((Ref**2 * Gain * numerator) / (dist1+dist2)**4))                        #then multiply by power received     
                + (delayedsig((dist+dist1+dist2)) * np.sqrt((Ref * Gain * numerator) / (dist**2 * (dist1 + dist2)**2)))
                + (delayedsig((dist+dist1+dist2)) * np.sqrt((Ref * LOS_Gain * numerator) / (dist**2 * (dist1 + dist2)**2))))      #receive straight on one
       
    
    mix = np.conj(Temp) * np.sqrt(TxPower)*np.exp(np.pi*1j*(fc*t+k*t**2))[:, np.newaxis]  #stretch process with original signal          
    RxSig[:,:] = mix[0:-1:int(fss/fs)][:]#Reduce sampling rate to 2e6 after stretch
    RxSig = add_zeros(RxSig)
    
    return RxSig    





data, Reflector = data_set()
multiple_rays = int(input('Enter 1 to include reflected rays, 0 to not: '))

RxSig = SimData(29,multi=multiple_rays)


'''
Plotting Functions Below
'''

##Create Polar Subplots and Adjust to correct orientation and limits

fig, ax = plt.subplots(1,2,subplot_kw=dict(projection='polar'))
ax[0].set_theta_zero_location("N"), ax[0].set_theta_direction(-1), ax[0].set_ylim([0, 5])   #create polar plots for comparing real data to simulation
ax[0].set_xlim([-np.pi/4, np.pi/4]), ax[0].set_xlabel('Angle /$^{\circ}$'), ax[0].set_ylabel('Range /m'), ax[0].set_title('Simulated Data')#, ax[0].grid(False), 

ax[1].set_theta_zero_location("N"), ax[1].set_theta_direction(-1), ax[1].set_ylim([0, 5])
ax[1].set_xlim([-np.pi/4, np.pi/4]), ax[1].set_xlabel('Angle /$^{\circ}$'), ax[1].set_ylabel('Range /m'), ax[1].set_title('Experimental Data')#, ax[1].grid(False)

levels = np.linspace(-35,0,100) # create contour levels for decibels


## Sim Data Plotting/Analysis
[rs,asim,Ran,Ang,zs] = RAngCalc(RxSig,k,fs)           #do data processing

cm = ax[0].contourf(np.radians(-asim),rs,20*np.log10(np.abs(zs)/np.amax(np.abs(zs))), levels=levels,cmap='jet', extend='both')
cb = fig.colorbar(cm,ax=ax[0],shrink=0.5)


## Real Data Plotting/Analysis
[r,a,Ra,An,z] = RAngCalc(data,k,fs) #do same as above but with the real data

cm1 = ax[1].contourf(np.radians(a), r,20*np.log10(np.abs(z)/np.amax(np.abs(z))), levels=levels, cmap='jet', extend='both') #
cb1 = fig.colorbar(cm1,ax=ax[1],shrink=0.5)



## Plot Range and Angle Profiles for Sim and Real Data for verification and troubleshooting

fig, (ax1, ax2) = plt.subplots(1, 2)            #plot range angle and contour profiles
fig2, (ax3, ax4) = plt.subplots(1, 2)  
                                                  
ax1.set_xlim([0.5, 5])
ax1.set_xlabel('Range /m')
ax1.set_title('Range Profile')

ax2.set_xlim([-45, 45])
ax2.set_xlabel('Signal AOA /$^{\circ}$')
ax2.set_title('Angular Profile')

ax3.set_xlim([0.5, 5])
ax3.set_xlabel('Range /m')

ax4.set_xlim([-45, 45])
ax4.set_xlabel('Signal AOA /$^{\circ}$')


ax1.plot(rs,np.abs(Ran[:,int(len(Ran.T) / 2)]))  # int function just finds middle element along array 

ax2.plot(-asim,np.abs(Ang[int(len(Ang) / 2),:])) # this is to avoid the array elements which are just 0

ax3.plot(r,np.abs(Ra[:,int(len(Ra.T) / 2)]))

ax4.plot(a,np.abs(An[int(len(An) / 2),:])) 

plt.show()
