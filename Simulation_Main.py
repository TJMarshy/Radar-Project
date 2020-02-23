import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio

data = sio.loadmat('experiment_1_fan_on.mat')
data = data['DATAr']


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
Gain = 1 #change to function of angle later

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

# Geom Sim 2D
# assume Tx element is a origin but assuming target radiates isotropically
# atm so doesnt matter 
RxArray = np.zeros((29,2))
Temp = np.zeros((len(t),29),dtype=complex)
RxSig = np.zeros((int(len(t)/(fss/fs)),29),dtype=complex)

Reflector = [[0,3.87,1],[-0.58,3.87,1],[1.25,1.94,1],[-1.10,1.35,1]]  #each element = [x,y,rcs]  maybe z later #maybe class aswell

for Obj in Reflector:

    for i in range(29):
        RxArray[i] = [i*lamb/2,0] #seperate 29 elements by distance lamb/2 
        dist = (np.sqrt((RxArray[i][0]-Obj[0])**2 + (RxArray[i][1]-Obj[1])**2)) #find distance between reflector and array element
        tdelay = (2*dist)/c
        P_received =  (TxPower * (Gain**2) * Obj[2] * (lamb**2)) / ((4*np.pi)**3 * dist**4)   #https://www.radartutorial.eu/01.basics/The%20Radar%20Range%20Equation.en.html
        print(P_received)
        Temp[:,i] = np.add(Temp[:,i],np.exp(np.pi*1j*(fc*np.subtract(t,tdelay)+k*(np.subtract(t,tdelay)**2)))) # produce delayed chirp

                      

mix = np.conj(Temp) * np.exp(np.pi*1j*(fc*t+k*t**2))[:, np.newaxis]  #stretch process with original signal          
RxSig[:,:] = mix[0:-1:int(fss/fs)][:]#Reduce sampling rate to 2e6 after stretch


fig, ax = plt.subplots(1,2,subplot_kw=dict(projection='polar'))
ax[0].set_theta_zero_location("N"), ax[0].set_theta_direction(-1), ax[0].set_ylim([0, 5])
ax[0].set_xlim([-np.pi/4, np.pi/4]), ax[0].grid(False)

ax[1].set_theta_zero_location("N"), ax[1].set_theta_direction(-1), ax[1].set_ylim([0, 5])
ax[1].set_xlim([-np.pi/4, np.pi/4]), ax[1].grid(False)



## Sim Data
[rs,asim,Ran,Ang,zs] = RAngCalc(RxSig,k,fs)           #do data processing

levels = np.linspace(0,1e4,100) 
cm = ax[0].contourf(np.radians(-asim),rs,np.abs(zs),levels=levels,cmap='plasma')
cb = fig.colorbar(cm,ax=ax[0],shrink=0.5)

# Real Data
[r,a,Ra,An,z] = RAngCalc(data,k,fs) #do same as above but with the real data

levels2 = np.linspace(0,2e4,100)
cm1 = ax[1].contourf(np.radians(a),r,np.abs(z),levels=levels2, cmap='plasma')
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