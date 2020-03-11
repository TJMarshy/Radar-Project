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
rxGain = 15
n1 = 1 #refractive indicies for fresenel
n2 = 2.5 #concrete floor https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=504907


def add_zeros(array):
    '''Add zero values before and after data to create a 
        rectangular pulse giving sinc functions after fft'''

    a = np.zeros((256,len(array.T)))
    b = np.append(a,array,0)
    new = np.append(b,a,0)          #for range sidelobes

    a = np.zeros((len(new),10))
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

    
    interp_ratio = 20                                #ratio of initial array elements vs interpolated elements
    
    Angle = np.linspace(-1,1,29*interp_ratio,endpoint=True) #create new angle vector to plot against
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

data,Reflector = Ex.fan_on()
#data,Reflector = Ex.fan_off()
#data,Reflector = Ex.corner2()
#data,Reflector = Ex.big1()
#data,Reflector = Ex.big2()

data = add_zeros(data)


# Geom Sim 3D
# assume Tx element is a origin but assuming target radiates isotropically
# atm so doesnt matter 

RxArray = [[i*lamb/2,0,0.98] for i in range(29)]   #create Uniform linear spaced array

Temp = np.zeros((len(t),len(RxArray)),dtype=complex)                #to store data before mixing/pulse compression?
RxSig = np.zeros((int(len(t)/(fss/fs)),len(RxArray)),dtype=complex) #smaller array to store post compression signal


for Obj in Reflector:  #loop through all objects and array elements

    numerator =  (TxPower * (txGain*rxGain) * Obj.rcs * (lamb**2)) / ((4*np.pi)**3)   #common numerator used for determining power           
    for i in range(len(RxArray)):
       
        X1 = (RxArray[i][0] * Obj.z + Obj.x * RxArray[i][2]) / (Obj.z + RxArray[i][2])              #find x,y,z position of reflection point
        Y1 = (RxArray[i][1] + Obj.y) /2
        Z1 = 0

        dist = np.sqrt( (RxArray[i][0]-Obj.x)**2 + (RxArray[i][1]-Obj.y)**2 + (RxArray[i][2]-Obj.z)**2 )       #Direct Ray distance
        dist1 = np.sqrt((RxArray[i][0] - X1)**2 + (RxArray[i][1] - Y1)**2 + (RxArray[i][2] - Z1)**2)           #array to ground distance
        dist2 = np.sqrt((Obj.x - X1)**2 + (Obj.y - Y1)**2 + (Obj.z - Z1)**2)                                   #ground to object distance

        theta = np.arctan((X1-Obj.x)/Obj.z)                                                                    #angle of incidence wrt the normal

        Rp = np.abs((n1*np.sqrt(1-(n1*np.sin(theta)/n2)**2) - n2*np.cos(theta))/(n1*np.sqrt(1-(n1*np.sin(theta)/n2)**2) + n2*np.cos(theta)))**2#assuming p polarised need to check
         
    
        Temp[:,i] += delayedsig(2*dist) * np.sqrt(numerator / (dist**4))                                       #sum each rays signal at each receiver
        + delayedsig(2*(dist1+dist2)) * np.sqrt((Rp**2 * numerator) / (dist1+dist2)**4)                        #then multiply by power received     
        + 2 * delayedsig((dist+dist1+dist2)) * np.sqrt((Rp * numerator) / (dist**2 * (dist1 + dist2)**2))      
       


                      

mix = np.conj(Temp) * np.sqrt(TxPower)*np.exp(np.pi*1j*(fc*t+k*t**2))[:, np.newaxis]  #stretch process with original signal          
RxSig[:,:] = mix[0:-1:int(fss/fs)][:]#Reduce sampling rate to 2e6 after stretch
RxSig = add_zeros(RxSig)

fig, ax = plt.subplots(1,2,subplot_kw=dict(projection='polar'))
ax[0].set_theta_zero_location("N"), ax[0].set_theta_direction(-1), ax[0].set_ylim([0, 5])   #create polar plots for comparing real data to simulation
ax[0].set_xlim([-np.pi/4, np.pi/4])#, ax[0].grid(False)

ax[1].set_theta_zero_location("N"), ax[1].set_theta_direction(-1), ax[1].set_ylim([0, 5])
ax[1].set_xlim([-np.pi/4, np.pi/4])#, ax[1].grid(False)



## Sim Data Plotting/Analysis

[rs,asim,Ran,Ang,zs] = RAngCalc(RxSig,k,fs)           #do data processing

levels = np.linspace(-30,0,100) #for deciblels
cm = ax[0].contourf(np.radians(-asim),rs,20*np.log10(np.abs(zs)/np.amax(np.abs(zs))), levels=levels,cmap='jet', extend='both')  #need -angle not sure why yet
cb = fig.colorbar(cm,ax=ax[0],shrink=0.5)

## Real Data Plotting/Analysis

[r,a,Ra,An,z] = RAngCalc(data,k,fs) #do same as above but with the real data

cm1 = ax[1].contourf(np.radians(a), r,20*np.log10(np.abs(z)/np.amax(np.abs(z))), levels=levels, cmap='jet', extend='both') #
cb1 = fig.colorbar(cm1,ax=ax[1],shrink=0.5)









#uncomment to plot range and angle maps and stuff / non polar contours etc....

fig, (ax1, ax2, ax3) = plt.subplots(1, 3)            #plot range angle and contour profiles
                                                  # N.B. have to reverse angle axes not sure why yet
ax1.set_xlim([0, 10])
ax2.set_xlim([-45, 45])
ax3.set_ylim([0, 10])


for i in range(1):
    ax1.plot(rs,np.abs(Ran[:,i+10]))
    ax2.plot(-asim,np.abs(Ang[i+280,:])) #280 is to get to points after all the zeros


   
cm = ax3.contourf(asim,rs,np.abs(zs),levels=levels)
cb = fig.colorbar(cm)


fig2, (ax5, ax6, ax7) = plt.subplots(1, 3)  

ax5.set_xlim([0, 5])
ax6.set_xlim([-45, 45])
ax7.set_ylim([0, 5])


for i in range(1):
    ax5.plot(r,np.abs(Ra[:,i+10]))
    ax6.plot(a,np.abs(An[i+300,:])) #300 is to get to points after all the zeros

cm2 = ax7.contourf(-a,r,np.abs(z),levels=levels)
cb2 = fig2.colorbar(cm2)



plt.show()












'''
    interp_ratio = 4
    Angle = np.linspace(-1,1,29*interp_ratio,endpoint=True)
    Angle = np.degrees(np.arcsin(Angle))  
    

    interp_mesh = np.array(np.meshgrid(Range, Angle, indexing='ij'))
    interp_points = np.rollaxis(interp_mesh, 0, 3).reshape((-1, 2)) #create list of cooridnates

    f = si.RegularGridInterpolator((Range,Angle1), A)
    A = f(interp_points)
    A = np.reshape(A,(len(Range), len(Angle)))

    f = si.RegularGridInterpolator((Range,Angle1), ReflecMap)
    ReflecMap = f(interp_points)
    ReflecMap = np.reshape(ReflecMap,(len(Range), len(Angle)))
'''