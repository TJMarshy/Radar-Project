import Experiments as Ex
import numpy as np
import matplotlib.pyplot as plt


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

def phased_sig(dist):
    '''delays the initial signal by pahse corresponding to 
       distance signal travelled'''  

    delsig = np.exp((1j*(dist)*2*np.pi)/lamb)
    
    return delsig

def dists(ran):
    Obj = Ex.Object(0,ran,0.83,1.72467)  #real
    Rx = Ex.Object(0,0,0.98,1)

    #Obj = Ex.Object(0,ran,.2,1)       #fake
    #Rx = Ex.Object(0,0,.2,1)


    Y1 = (Rx.y * Obj.z + Obj.y * Rx.z) / (Obj.z + Rx.z)              #find x,y,z position of reflection point
    X1 = (Rx.x + Obj.x) /2
    Z1 = 0

    dist = np.sqrt( (Rx.x-Obj.x)**2 + (Rx.y-Obj.y)**2 + (Rx.z-Obj.z)**2 )       #Direct Ray distance
    dist1 = np.sqrt((Rx.x - X1)**2 + (Rx.y - Y1)**2 + (Rx.z - Z1)**2)           #array to ground distance
    dist2 = np.sqrt((Obj.x - X1)**2 + (Obj.y - Y1)**2 + (Obj.z - Z1)**2)          #ground to object distance

    return [dist,dist1,dist2]
                             



ran = np.linspace(0.2,50,10000)

amp = (phased_sig(2*dists(ran)[0]) + phased_sig(2*(dists(ran)[1]+ dists(ran)[2])) + 2*phased_sig(dists(ran)[0] + dists(ran)[1] + dists(ran)[2]))/3


plt.plot(ran,20*np.log10(np.abs(amp)/3))
plt.xlabel('Range /m')
plt.ylabel('Normalised Amplitude /dB')
plt.title('Received Signal Power for 4-Ray Model, with Real Positions')
plt.show()