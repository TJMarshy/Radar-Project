import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import scipy.interpolate as si

data = sio.loadmat('experiment_1_fan_on.mat')

data = data['DATAr'] *1e-6   #just to get decibels of real data to 0 , no understanding yet


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
txGain = 17 #change to function of angle later
rxGain = 15

def RAngCalc(data,sweep_slope,fs):
    R = np.fft.fft(data,None,0)

    A = np.fft.fft(data,None,1)
    A = np.fft.fftshift(A,1)

    ReflecMap = np.fft.fft(R,None,1)
    ReflecMap = np.fft.fftshift(ReflecMap,1)

    Angle = np.linspace(-1,1,len(data.T),endpoint=True)
    Angle = np.degrees(np.arcsin(Angle))

    Freq = np.linspace(0,fs,len(data),endpoint=True)   
    Range = (Freq*3e8) / (2*sweep_slope)

    return [Range, Angle, R, A, ReflecMap]

def delayedsig(dist,Obj):   #NEED TO INCLUDE DIFFERENT DISTANCES LEADING TO DIFFERENT POWER AND GROUND RADAR CROSS SECTION TYPE STUFF

    delay = dist/c
    delsig = np.exp(np.pi*1j*(fc*np.subtract(t,delay)+k*(np.subtract(t,delay)**2)))
    P_received =  (TxPower * (txGain*rxGain) * Obj.rcs * (lamb**2)) / ((4*np.pi)**3 * dist**4) #https://www.radartutorial.eu/01.basics/The%20Radar%20Range%20Equation.en.html

    sig = np.sqrt(P_received) * delsig
    
    return sig


class Object():
    '''Define object to give its parameters'''

    def __init__(self, x, y, z, rcs):

        self.x = x
        self.y = y
        self.z = z
        self.rcs = rcs

Corn1, Corn2, Fan, Sphere = Object(0,3.87,0.83,5), Object(-0.58,3.87,0.83,5), Object(1.25, 1.94, .95, 1), Object(-1.10, 1.35, 0.89, 1)
Reflector = [Corn1,Corn2, Fan, Sphere]  #each element = [x,y,rcs]  maybe z later #maybe class aswell


# Geom Sim 2D
# assume Tx element is a origin but assuming target radiates isotropically
# atm so doesnt matter 
#RxArray = np.zeros((29,2))
RxArray = [[i*lamb/2,0,0.98] for i in range(29)]
Temp = np.zeros((len(t),29),dtype=complex)
RxSig = np.zeros((int(len(t)/(fss/fs)),29),dtype=complex)


for Obj in Reflector:
    for i in range(29):

        dist = np.sqrt( (RxArray[i][0]-Obj.x)**2 + (RxArray[i][1]-Obj.y)**2 + (RxArray[i][2]-Obj.z)**2 )       #Direct Ray
        X1 = (RxArray[i][0] * Obj.z + Obj.x * RxArray[i][2]) / (Obj.z + RxArray[i][2])              #find x,y,z position of reflection point
        Y1 = (RxArray[i][1] + Obj.y) /2
        Z1 = 0

        dist1 = np.sqrt((RxArray[i][0] - X1)**2 + (RxArray[i][1] - Y1)**2 + (RxArray[i][2] - Z1)**2)  #array to ground
        dist2 = np.sqrt((Obj.x - X1)**2 + (Obj.y - Y1)**2 + (Obj.z - Z1)**2)                          #ground to object

        
           
        Temp[:,i] += delayedsig(2*dist,Obj) + delayedsig(2*(dist1+dist2),Obj) + delayedsig((dist+dist1+dist2),Obj) + delayedsig((dist+dist1+dist2),Obj)
        #REMEMEBER THE ORDER PROBABLY MATTER EXPLORE MORE


                      

mix = np.conj(Temp) * np.sqrt(TxPower)*np.exp(np.pi*1j*(fc*t+k*t**2))[:, np.newaxis]  #stretch process with original signal          
RxSig[:,:] = mix[0:-1:int(fss/fs)][:]#Reduce sampling rate to 2e6 after stretch


fig, ax = plt.subplots(1,2,subplot_kw=dict(projection='polar'))
ax[0].set_theta_zero_location("N"), ax[0].set_theta_direction(-1), ax[0].set_ylim([0, 5])
ax[0].set_xlim([-np.pi/4, np.pi/4])#, ax[0].grid(False)

ax[1].set_theta_zero_location("N"), ax[1].set_theta_direction(-1), ax[1].set_ylim([0, 5])
ax[1].set_xlim([-np.pi/4, np.pi/4])#, ax[1].grid(False)



## Sim Data
[rs,asim,Ran,Ang,zs] = RAngCalc(RxSig,k,fs)           #do data processing

#levels = np.linspace(0,1e4,100) 
levels = np.linspace(-30,0,100)
cm = ax[0].contourf(np.radians(-asim),rs,20*np.log10(np.abs(zs)/np.sqrt(TxPower)), levels=levels,cmap='jet', extend='both')   #
cb = fig.colorbar(cm,ax=ax[0],shrink=0.5)

# Real Data
[r,a,Ra,An,z] = RAngCalc(data,k,fs) #do same as above but with the real data

#levels2 = np.linspace(0,2e4,100)
levels2 = np.linspace(-30,0,100)
cm1 = ax[1].contourf(np.radians(a), r,20*np.log10(np.abs(z)/np.sqrt(TxPower)), levels=levels2, cmap='jet', extend='both') #
cb1 = fig.colorbar(cm1,ax=ax[1],shrink=0.5)





plt.show()



'''
fig, (ax1, ax2, ax3) = plt.subplots(1, 3)            #plot range angle and contour profiles
                                                  # N.B. have to reverse angle axes not sure why yet
ax1.set_xlim([0, 10])
ax2.set_xlim([-45, 45])
ax3.set_ylim([0, 10])


for i in range(29):
    ax1.plot(rs,np.abs(Ran[:,i]))
    ax2.plot(-asim,np.abs(Ang[i,:]))


   
cm = ax3.contourf(-asim,rs,np.abs(zs),levels=levels)
cb = fig.colorbar(cm)
'''

'''
fig2, (ax5, ax6, ax7) = plt.subplots(1, 3)  

ax5.set_xlim([0, 5])
ax6.set_xlim([-45, 45])
ax7.set_ylim([0, 5])


for i in range(29):
    ax5.plot(r,np.abs(Ra[:,i]))
    ax6.plot(-a,np.abs(An[i,:]))

cm2 = ax7.contourf(-a,r,np.abs(z),levels=levels)
cb2 = fig2.colorbar(cm2)

'''