# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 16:48:28 2022

@author: vankronp

Modified by ertmana May 2022 
Modified by MarIniOnz July 2022

"""
#%matplotlib qt5

from IPython import get_ipython
get_ipython().magic('%reset -f') 

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import gettingTriggers

root_NP = # directory

plotting = 1 #1

##########################

filenames = #where the data is
 
keys = ['spikes', 'clusters', 'triggers','depths']
soundnames = ['FRA3.5to60', '8F', '4F4F', '4F4R', '4F4FRC', '4FSNT', '4FS', '4F4F6R', 'AEM','RE', 'NRE', 'OE', 'NRUE', 'OUE', 'RUE','APM', 'AEMSD', 'RESD', 'NRESD', 'OESD', 'NRUESD', 'OUESD', 'RUESD', 'APMSD', 'PS6'] 

sounds = np.arange(len(soundnames)) # num of total sequences (including FRAs)
something_else = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] #ordered indices for FRAs and FSs, each starting at 0 #mostly will be 0s, basically how many times everthing got played. 

spike_dict = dict.fromkeys(keys)
for key in keys:
    spike_dict[key] = dict.fromkeys(soundnames)

##############################################################################

#mi, filename = 0, filenames[0]
for mi, filename in enumerate(filenames):
    
    spike_times = np.load(root_NP + filename + '\spike_times.npy')
    spike_clusters = np.load(root_NP + filename + r'\spike_clusters.npy')
    #spike_depths = np.load(root_NP + filename + r'\spike_depths.npy') # Comment this if this file does not exist with your version of processing
    
    # plt.figure() # commented
    # plt.scatter(spike_times, spike_clusters, s=1, marker='*', c='k', alpha=0.1) # commented
     
    binFullPath_AP = Path(root_NP + filename + '\\' + filename[-23:-6] + '_t0.imec0.ap.bin')
    binFullPath_NI = Path(root_NP + filename[:-6] + '_t0.nidq.bin')
        
    tstrigNPAP, triggers_all_NP, FRA_indices, FS_indices, sRate_AP, sRate_NI = gettingTriggers_Martin.getTrigsMartin(binFullPath_AP, binFullPath_NI, soundnames, Alexandras_lengths, plotting)
    
    for i, sound in enumerate(sounds):

        spike_times_adj = (spike_times/sRate_AP)
        
        first_trig_NP = tstrigNPAP[sound]
        
        spikes_in_trial= spike_times_adj[np.where(np.logical_and(spike_times_adj >= first_trig_NP-0.1, spike_times_adj <= first_trig_NP+Alexandras_lengths[i]))[0]]
        spike_clusters_trial = spike_clusters[np.where(np.logical_and(spike_times_adj >= first_trig_NP-0.1, spike_times_adj <= first_trig_NP+Alexandras_lengths[i]))[0]]
        
        spike_dict[keys[0]][soundnames[i]]=spikes_in_trial
        spike_dict[keys[1]][soundnames[i]]=spike_clusters_trial# We commented this
        
        triggers = (triggers_all_NP[soundnames[i]]-1000)/sRate_NI + tstrigNPAP[sounds[i]]
        
        spike_dict[keys[2]][soundnames[i]] = triggers
    # scipy.io.savemat(root_NP + filename + '\spikes_by_sound.mat', spike_dict)
        
        
