
import scipy.io as sio


class Object():
    '''Define object to give its parameters'''

    def __init__(self, x, y, z, rcs):

        self.x = x
        self.y = y
        self.z = z
        self.rcs = rcs

'''
Functions below, load data from a given experiment
and the corresponding points to simulate
'''

def fan_on():
    data = sio.loadmat('experiment_1_fan_on.mat')
    data = data['DATAr']

    Corn1, Corn2, Fan, Sphere = Object(0,3.87,0.83,1.72467), Object(-0.58,3.87,0.83,1.72467), Object(1.25, 1.94, .95, 0.5), Object(-1.10, 1.35, 0.89, 0.01337)
    
    return data , [Corn1,Corn2, Fan, Sphere]  #create array of all objects with their x,y,z and rcs values need value for fan off/on
    

def fan_off():
    data = sio.loadmat('experiment_2_fan_off.mat')
    data = data['DATAr']

    Corn1, Corn2, Fan, Sphere = Object(0,3.87,0.83,1.72467), Object(-0.58,3.87,0.83,1.72467), Object(1.25, 1.94, .95, 2.587), Object(-1.10, 1.35, 0.89, 0.01337)
    
    return data , [Corn1,Corn2, Fan, Sphere]  #create array of all objects with their x,y,z and rcs values need value for fan off/on

def corner2():
    data = sio.loadmat('experiment_3_2cornerref.mat')
    data = data['DATAr']

    Corn1, Corn2 = Object(0,3.87,0.83,1.72467), Object(-0.58,3.87,0.83,1.72467)
    
    return data , [Corn1,Corn2] #create array of all objects with their x,y,z and rcs values need value for fan off/on


def big1():
    data = sio.loadmat('experiment_4_bigreflector.mat')
    data = data['DATAr']

    Corn1, Corn2 = Object(0,3.87,0.83,1.72467), Object(-0.58,3.87,0.83,1.72467)

    return data, [Corn1,Corn2]


def big2():
    data = sio.loadmat('experiment_5_bigerflector.mat')
    data = data['DATAr']
 
    Corn1, Corn2 = Object(0,2,0.83,1.72467), Object(-0.58,2,0.83,1.72467)

    return data, [Corn1,Corn2]


def test():

    data = sio.loadmat('experiment_5_bigerflector.mat')
    data = data['DATAr']

    Ob, Ob2 = Object(0,3,0.83,1) ,Object(0.217,3,0.83,1)

    return data, [Ob, Ob2]