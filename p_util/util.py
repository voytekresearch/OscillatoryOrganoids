import numpy as np
import os
import re
import matplotlib.pyplot as plt
import OpenEphys
import glob
from matplotlib.ticker import MaxNLocator

def getFiles(dataLocation):
    '''
    Returns all names of files in the specified folder that resemble OpenEphys save files.
    Gather channels 1-9 first, then 10-32 and add them together. Because sorted() will not sort 1-32 properly.
    '''
    files = sorted([f for f in os.listdir(dataLocation) if re.search(r'.*CH[0-9].continuous$', f)])
    return files+sorted([f for f in os.listdir(dataLocation) if re.search(r'.*CH[0-9][0-9].continuous$', f)])
def getChannels(dataLocation):
    '''
    Gets all the channels from designated folder and returns them to a numpy array.
    '''
    files=getFiles(dataLocation)
    channels=[]
    for i, x in enumerate(files):
        raw=OpenEphys.load(dataLocation+'/'+x)
        channels.append(np.array(raw['data']))
    return channels
def getAverages(dataLocation, method='mean'):
    '''
    Gets averages of each channel and organizes it in a array by relative location on the shank.

    Parameters
    ----------
    method : string
        Put mean to return the means of all channels, and std first standard deviations. Defaults to mean.
    '''
    channelMap=[[5,4,6,3,7,2,8,1],[13,12,14,11,15,10,16,9],[21,20,22,19,23,18,24,17],[29,28,30,27,31,26,32,25]]# This map is organized by shank, and then by height on shank (e.g 5,13,21,29 being the lowest)
    channels=getChannels(dataLocation)
    shankAvgs=[]
    for i, x in enumerate(channelMap):
        channelAverages=[]
        for v in x:
            if method=='mean':
                processedChannel=np.mean(channels[v-1])
            if method=='std':
                processedChannel=np.std(channels[v-1])
            channelAverages.append(processedChannel)
        shankAvgs.append(np.array(channelAverages))
    return shankAvgs
'''
def getChannelDict():
    channelMap= {
            shank_1: {
                left: {
                    CH_4:4,
                    CH_3:3,
                    CH_2:2,
                    CH_1:1
                    },
                right: {
                    CH_5:5,
                    CH_6:6,
                    CH_7:7,
                    CH_8:8
                    }
                },
            shank_2: {
                left: {
                    CH_12:12,
                    CH_11:11,
                    CH_10:10,
                    CH_9:9
                    },
                right: {
                    CH_13:13,
                    CH_14:14,
                    CH_15:15,
                    CH_16:16
                    }
                },
            shank_3: {
                left:{
                    CH_20:20,
                    CH_19:19,
                    CH_18:18,
                    CH_17:17
                    }
                right:{
                    CH_21:21,
                    CH_22:22,
                    CH_23:23,
                    CH_24:24
                    }
                },
            shank_4: {
                left:[28,27,26,25],
                right:[29,30,31,32]
                }
            }
'''
