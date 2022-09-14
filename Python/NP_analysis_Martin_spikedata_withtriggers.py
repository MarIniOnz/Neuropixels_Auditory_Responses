# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 16:48:28 2022

@author: vankronp

Modified by ertmana May 2022 

"""

#cd S:\Fakultaet\MFZ\NWFZ\AGdeHoz\Martni\Codes\Python\
#%matplotlib qt5

from IPython import get_ipython
get_ipython().magic('%reset -f') 


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import gettingTriggers_Martin
# import scipy.io


root_NP = r'S:\Fakultaet\MFZ\NWFZ\AG-deHoz-Scratch\Neuropixels'

plotting = 1 #1

#1702
#----------------------------------------------------------------------------
# filenames = [r'\\26_5_2022\m1702_26_05_22_g0\m1702_26_05_22_g0_imec0'] #where the data is   

filenames = [r'\\15_10_2020\GG_M609_b__g0_processed\GG_M609_b__g0_imec0\imec0_ks2_orig'] #where the data is
 
keys = ['spikes', 'clusters', 'triggers','depths']
# soundnames = ['FRA3.5to60', '8F', '4F4F', '4F4R', '4F4FRC', '4FSNT', '4FS', '4F4F6R', 'AEM','RE', 'NRE', 'OE', 'NRUE', 'OUE', 'RUE','APM', 'AEMSD', 'RESD', 'NRESD', 'OESD', 'NRUESD', 'OUESD', 'RUESD', 'APMSD', 'PS6'] #1701 25.06.22 #should FRA2 be called FRA1 bc they're the same sound? 
soundnames = ['FRA3.5to60', '8F', '4F4F', '4F4R', '4F4FRC', '4FSNT', '4FS', '4F4F6R', 'AEM','RE', 'NRE', 'OE', 'NRUE', 'OUE', 'RUE','APM', 'AEMSD', 'RESD', 'NRESD', 'OESD', 'NRUESD', 'OUESD', 'RUESD', 'APMSD', 'PS6'] 

Alexandras_sounds = np.arange(len(soundnames)) # num of total sequences (including FRAs)
Alexandras_lengths = np.array([530, 133, 133, 133, 133, 133, 133, 117, 209, 69, 69, 69, 23, 23, 23, 279, 209, 69, 69, 69, 23, 23, 23, 279, 310])
something_else = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] #ordered indices for FRAs and FSs, each starting at 0 #mostly will be 0s, basically how many times everthing got played. 


spike_dict = dict.fromkeys(keys)
for key in keys:
    spike_dict[key] = dict.fromkeys(soundnames)

##############################################################################

#mi, filename = 0, filenames[0]
for mi, filename in enumerate(filenames):
    
    spike_times = np.load(root_NP + filename + '\spike_times.npy')
    spike_clusters = np.load(root_NP + filename + r'\spike_clusters.npy')
    #spike_depths = np.load(root_NP + filename + r'\spike_depths.npy') # We commented this, we Do not have this file
    
    # plt.figure() # commented
    # plt.scatter(spike_times, spike_clusters, s=1, marker='*', c='k', alpha=0.1) # commented
     
    binFullPath_AP = Path(root_NP + filename + '\\' + filename[-23:-6] + '_t0.imec0.ap.bin')
    binFullPath_NI = Path(root_NP + filename[:-6] + '_t0.nidq.bin')
        
    # tstrigNPAP, tstrigNP_FRA_all, FRA_indices, triggers_allNP, sRate_AP, sRate_NI = Martin_trying_Hard_with_poppy.getTrigsPoppyMartin(binFullPath_AP, binFullPath_NI, soundnames, Alexandras_lengths, plotting)
    # tstrigNPAP, tstrigNP_FRA_all, tstrigNP_FS_all, sync_offset, sRate_AP, sRate_NI = Spike_GLX_Poppy_Analog_Trigger_loading_function.getTrigs(binFullPath_AP, binFullPath_NI, plotting)
    tstrigNPAP, triggers_all_NP, FRA_indices, FS_indices, sRate_AP, sRate_NI = gettingTriggers_Martin.getTrigsMartin(binFullPath_AP, binFullPath_NI, soundnames, Alexandras_lengths, plotting)
    
    for i, sound in enumerate(Alexandras_sounds):

        # spike_times_adj = (spike_times/sRate_AP)+sync_offset[sound]
        spike_times_adj = (spike_times/sRate_AP)
        
        first_trig_NP = tstrigNPAP[sound]
        
        spikes_in_trial= spike_times_adj[np.where(np.logical_and(spike_times_adj >= first_trig_NP-0.1, spike_times_adj <= first_trig_NP+Alexandras_lengths[i]))[0]]
        spike_clusters_trial = spike_clusters[np.where(np.logical_and(spike_times_adj >= first_trig_NP-0.1, spike_times_adj <= first_trig_NP+Alexandras_lengths[i]))[0]]
        # spike_depths_trial= spike_depths[np.where(np.logical_and(spike_times_adj >= first_trig_NP-0.1, spike_times_adj <= first_trig_NP+Alexandras_lengths[i]))[0]] # We Do not have this file

        
        spike_dict[keys[0]][soundnames[i]]=spikes_in_trial
        spike_dict[keys[1]][soundnames[i]]=spike_clusters_trial# We commented this
        
        # if Alexandras_lengths[i] < 180:
        #     triggers = (tstrigNP_FS_all[something_else[i]]-1000)/sRate_NI + tstrigNPAP[Alexandras_sounds[i]]
        # else:
        #     triggers = (tstrigNP_FRA_all[something_else[i]]-1000)/sRate_NI + tstrigNPAP[Alexandras_sounds[i]]
        
        triggers = (triggers_all_NP[soundnames[i]]-1000)/sRate_NI + tstrigNPAP[Alexandras_sounds[i]]
        
        spike_dict[keys[2]][soundnames[i]] = triggers
        #spike_dict[keys[3]][soundnames[i]] = spike_depths_trial# We commented this
        
    # scipy.io.savemat(root_NP + filename + '\spikes_by_sound.mat', spike_dict)
        
        