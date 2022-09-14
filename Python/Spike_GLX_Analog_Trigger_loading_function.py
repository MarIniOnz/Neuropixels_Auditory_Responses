# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 16:44:46 2021

@author: vankronp

Modified by ertmana May 2022
Modified by MarIniOnz July 2022

reading Analog triggers from AP bin file and from Analog bin file

"""

def getTrigs(binFullPath_AP, binFullPath_NI, plotting):
    
    from pathlib import Path
    import numpy as np
    import matplotlib.pyplot as plt
    
    #necessary functions:
    def readMeta(binFullPath):
        metaName = binFullPath.stem + ".meta"
        metaPath = Path(binFullPath.parent / metaName)
        metaDict = {}
        if metaPath.exists():
            # print("meta file present")
            with metaPath.open() as f:
                mdatList = f.read().splitlines()
                # convert the list entries into key value pairs
                for m in mdatList:
                    csList = m.split(sep='=')
                    if csList[0][0] == '~':
                        currKey = csList[0][1:len(csList[0])]
                    else:
                        currKey = csList[0]
                    metaDict.update({currKey: csList[1]})
        else:
            print("no meta file")
        return(metaDict)
    
    
    # Return sample rate as python float.
    # On most systems, this will be implemented as C++ double.
    # Use python command sys.float_info to get properties of float on your system.
    #
    def SampRate(meta):
        if meta['typeThis'] == 'imec':
            srate = float(meta['imSampRate'])
        else:
            srate = float(meta['niSampRate'])
        return(srate)
    
    
    # Return a multiplicative factor for converting 16-bit file data
    # to voltage. This does not take gain into account. The full
    # conversion with gain is:
    #         dataVolts = dataInt * fI2V / gain
    # Note that each channel may have its own gain.
    #
    def Int2Volts(meta):
        if meta['typeThis'] == 'imec':
            if 'imMaxInt' in meta:
                maxInt = int(meta['imMaxInt'])
            else:
                maxInt = 512
            fI2V = float(meta['imAiRangeMax'])/maxInt
        else:
            fI2V = float(meta['niAiRangeMax'])/32768
        return(fI2V)
    
    
    # Return array of original channel IDs. As an example, suppose we want the
    # imec gain for the ith channel stored in the binary data. A gain array
    # can be obtained using ChanGainsIM(), but we need an original channel
    # index to do the lookup. Because you can selectively save channels, the
    # ith channel in the file isn't necessarily the ith acquired channel.
    # Use this function to convert from ith stored to original index.
    # Note that the SpikeGLX channels are 0 based.
    #
    def OriginalChans(meta):
        if meta['snsSaveChanSubset'] == 'all':
            # output = int32, 0 to nSavedChans - 1
            chans = np.arange(0, int(meta['nSavedChans']))
        else:
            # parse the snsSaveChanSubset string
            # split at commas
            chStrList = meta['snsSaveChanSubset'].split(sep=',')
            chans = np.arange(0, 0)  # creates an empty array of int32
            for sL in chStrList:
                currList = sL.split(sep=':')
                if len(currList) > 1:
                    # each set of contiguous channels specified by
                    # chan1:chan2 inclusive
                    newChans = np.arange(int(currList[0]), int(currList[1])+1)
                else:
                    newChans = np.arange(int(currList[0]), int(currList[0])+1)
                chans = np.append(chans, newChans)
        return(chans)
    
    
    # Return counts of each nidq channel type that composes the timepoints
    # stored in the binary file.
    #
    def ChannelCountsNI(meta):
        chanCountList = meta['snsMnMaXaDw'].split(sep=',')
        MN = int(chanCountList[0])
        MA = int(chanCountList[1])
        XA = int(chanCountList[2])
        DW = int(chanCountList[3])
        return(MN, MA, XA, DW)
    
    
    # Return counts of each imec channel type that composes the timepoints
    # stored in the binary files.
    #
    def ChannelCountsIM(meta):
        chanCountList = meta['snsApLfSy'].split(sep=',')
        AP = int(chanCountList[0])
        LF = int(chanCountList[1])
        SY = int(chanCountList[2])
        return(AP, LF, SY)
    
    
    # Return gain for ith channel stored in nidq file.
    # ichan is a saved channel index, rather than the original (acquired) index.
    #
    def ChanGainNI(ichan, savedMN, savedMA, meta):
        if ichan < savedMN:
            gain = float(meta['niMNGain'])
        elif ichan < (savedMN + savedMA):
            gain = float(meta['niMAGain'])
        else:
            gain = 1    # non multiplexed channels have no extra gain
        return(gain)
    
    
    # Return gain for imec channels.
    # Index into these with the original (acquired) channel IDs.
    #
    def ChanGainsIM(meta):
        imroList = meta['imroTbl'].split(sep=')')
        # One entry for each channel plus header entry,
        # plus a final empty entry following the last ')'
        nChan = len(imroList) - 2
        APgain = np.zeros(nChan)        # default type = float
        LFgain = np.zeros(nChan)
        if 'imDatPrb_type' in meta:
            probeType = meta['imDatPrb_type']
        else:
            probeType = 0
        if (probeType == 21) or (probeType == 24):
            # NP 2.0; APGain = 80 for all AP
            # return 0 for LFgain (no LF channels)
            APgain = APgain + 80
        else:
            # 3A, 3B1, 3B2 (NP 1.0)
            for i in range(0, nChan):
                currList = imroList[i+1].split(sep=' ')
                APgain[i] = currList[3]
                LFgain[i] = currList[4]
        return(APgain, LFgain)
    
    
    # Having accessed a block of raw nidq data using makeMemMapRaw, convert
    # values to gain-corrected voltage. The conversion is only applied to the
    # saved-channel indices in chanList. Remember, saved-channel indices are
    # in the range [0:nSavedChans-1]. The dimensions of dataArray remain
    # unchanged. ChanList examples:
    # [0:MN-1]  all MN channels (MN from ChannelCountsNI)
    # [2,6,20]  just these three channels (zero based, as they appear in SGLX).
    #
    def GainCorrectNI(dataArray, chanList, meta):
        MN, MA, XA, DW = ChannelCountsNI(meta)
        fI2V = Int2Volts(meta)
        # print statements used for testing...
        # print("NI fI2V: %.3e" % (fI2V))
        # print("NI ChanGainNI: %.3f" % (ChanGainNI(0, MN, MA, meta)))
    
        # make array of floats to return. dataArray contains only the channels
        # in chanList, so output matches that shape
        convArray = np.zeros(dataArray.shape, dtype=float)
        for i in range(0, len(chanList)):
            j = chanList[i]             # index into timepoint
            conv = fI2V/ChanGainNI(j, MN, MA, meta)
            # dataArray contains only the channels in chanList
            convArray[i, :] = dataArray[i, :] * conv
        return(convArray)
    
    
    # Having accessed a block of raw imec data using makeMemMapRaw, convert
    # values to gain corrected voltages. The conversion is only applied to
    # the saved-channel indices in chanList. Remember saved-channel indices
    # are in the range [0:nSavedChans-1]. The dimensions of the dataArray
    # remain unchanged. ChanList examples:
    # [0:AP-1]  all AP channels
    # [2,6,20]  just these three channels (zero based)
    # Remember that for an lf file, the saved channel indices (fetched by
    # OriginalChans) will be in the range 384-767 for a standard 3A or 3B probe.
    #
    def GainCorrectIM(dataArray, chanList, meta):
        # Look up gain with acquired channel ID
        chans = OriginalChans(meta)
        APgain, LFgain = ChanGainsIM(meta)
        nAP = len(APgain)
        nNu = nAP * 2
    
        # Common conversion factor
        fI2V = Int2Volts(meta)
    
        # make array of floats to return. dataArray contains only the channels
        # in chanList, so output matches that shape
        convArray = np.zeros(dataArray.shape, dtype='float')
        for i in range(0, len(chanList)):
            j = chanList[i]     # index into timepoint
            k = chans[j]        # acquisition index
            if k < nAP:
                conv = fI2V / APgain[k]
            elif k < nNu:
                conv = fI2V / LFgain[k - nAP]
            else:
                conv = 1
            # The dataArray contains only the channels in chanList
            convArray[i, :] = dataArray[i, :]*conv
        return(convArray)
    
    
    def makeMemMapRaw(binFullPath, meta):
        nChan = int(meta['nSavedChans'])
        nFileSamp = int(int(meta['fileSizeBytes'])/(2*nChan))
        print("nChan: %d, nFileSamp: %d" % (nChan, nFileSamp))
        rawData = np.memmap(binFullPath, dtype='int16', mode='r',
                            shape=(nChan, nFileSamp), offset=0, order='F')
        return(rawData)
    
    
    # Return an array [lines X timepoints] of uint8 values for a
    # specified set of digital lines.
    #
    # - dwReq is the zero-based index into the saved file of the
    #    16-bit word that contains the digital lines of interest.
    # - dLineList is a zero-based list of one or more lines/bits
    #    to scan from word dwReq.
    #
    def ExtractDigital(rawData, firstSamp, lastSamp, dwReq, dLineList, meta):
        # Get channel index of requested digital word dwReq
        if meta['typeThis'] == 'imec':
            AP, LF, SY = ChannelCountsIM(meta)
            if SY == 0:
                print("No imec sync channel saved.")
                digArray = np.zeros((0), 'uint8')
                return(digArray)
            else:
                digCh = AP + LF + dwReq
        else:
            MN, MA, XA, DW = ChannelCountsNI(meta)
            if dwReq > DW-1:
                print("Maximum digital word in file = %d" % (DW-1))
                digArray = np.zeros((0), 'uint8')
                return(digArray)
            else:
                digCh = MN + MA + XA + dwReq
    
        selectData = np.ascontiguousarray(rawData[digCh, firstSamp:lastSamp+1], 'int16')
        nSamp = lastSamp-firstSamp + 1
    
        # unpack bits of selectData; unpack bits works with uint8
        # original data is int16
        bitWiseData = np.unpackbits(selectData.view(dtype='uint8'))
        # output is 1-D array, nSamp*16. Reshape and transpose
        bitWiseData = np.transpose(np.reshape(bitWiseData, (nSamp, 16)))
    
        nLine = len(dLineList)
        digArray = np.zeros((nLine, nSamp), 'uint8')
        for i in range(0, nLine):
            byteN, bitN = np.divmod(dLineList[i], 8)
            targI = byteN*8 + (7 - bitN)
            digArray[i, :] = bitWiseData[targI, :]
        return(digArray)
    
    
    binFullPath = binFullPath_NI 
    meta = readMeta(binFullPath)
    sRate_NI = SampRate(meta)
    rawData_NI = makeMemMapRaw(binFullPath, meta)
    
    # array of times for plot
    
    chanList = [0, 2] 
    selectData_NI = rawData_NI[chanList, :]
    NI_convData = 1e3*GainCorrectNI(selectData_NI, chanList, meta) ##crashing here 1701 day2
    tstrigtrace = NI_convData[0, :]
    #plt.plot(tstrigtrace[::100]) plotting when weird triggers occur
    triggs = np.where(tstrigtrace > 2000)[0] #the 2000 is determining the trigger threshold
    
    tstrigNP = np.where(np.diff(triggs)>90000)[0] #90000 before #this is the time between blocks i think #179600?
    
    tstrigNP = tstrigNP+1 #
    tstrigNP = np.insert(tstrigNP, 0, 0) #first trigger of every sound that it detected
    
    #tstrigNP_FRA = np.array([tstrigNP[0], tstrigNP[-1]]) # brackets are indicies of the sounds you want it to extract
    tstrigNP_FRA = np.array([tstrigNP[0]]) 
    
    #tstrigNP_FRA_all = np.zeros((2, 3700)) #for 1701 25.05.22. #freqs = 74, levels = 10, reps = 5, 74*10*5 = 3700 #how many triggers it should look for in each file #"2" - how many files, other number freqXlvlXreps
    #tstrigNP_FRA_all = np.zeros((1, 2295)) #for 1701,1702,1703 26.05.22. #freqs = 32 levels = 10 resps = 5,  32*10*5 = 2295 ... this is not true, its 1600
    tstrigNP_FRA_all = np.zeros((1, 1600)) #for 1701,1702,1703 26.05.22. #freqs = 32 levels = 10 resps = 5,  32*10*5 = 2295
    ##? both 1702 and 1703 say --- ValueError: could not broadcast input array from shape (1600,) into shape (2295,) ???
    #tstrigNP_FRA_all = np.zeros((1, 620)) #490 FS? 620?
    
    time_before_trig = 1000 #for baseline activity
    
    #time_after_trig = 12635000 # 1701 day 25.05.22 #7760000 #to do with the length of the file you figured out before #starts from 0 for each sound #length of sound X sRate_NI (then round up) #146700-20350 = 126350. times by 100? #ask philipp #yes
    #time_after_trig = 5405000 #for 1701 26.05.22. #2400 56450. 56450-2400 = 54050. x100 = 5405000
    time_after_trig = 5400000 #for 1702 26.05.22. #1150 55150. 55150 - 1150 = 54000. x100 = 5400000
    #time_after_trig = 5430000 #for 1703 26.05.22. #980 55280. 55280 - 980 = 54300. x100 = 5430000

    #time_after_trig = 2065000 ##fs day 28.05.22 #1100 21750, 21750 - 1100 = 20650
#getting all triggers of each sound presentation
    
    for i, trigger in enumerate(tstrigNP_FRA):
        if i == 0:
            #trigtrace_temp = tstrigtrace[(triggs[trigger]-time_before_trig):triggs[tstrigNP[1]-1]] # added the parenthesis
            trigtrace_temp = tstrigtrace[triggs[trigger]-time_before_trig:triggs[tstrigNP[1]-1]]
        else:
            #trigtrace_temp = tstrigtrace[(triggs[trigger]-time_before_trig):triggs[trigger]+time_after_trig]
            trigtrace_temp = tstrigtrace[triggs[trigger]-time_before_trig:triggs[trigger]+time_after_trig]
        triggs_trial = np.where(trigtrace_temp > 2000)[0]   #the inter trigger interval, so find what the ITI is, x60, x sRate_NI (then round up, for example, 1998 -> 2000)
        last_trigs = np.where(np.diff(triggs_trial)>2000)[0]
        tstrigNP_FRA_all[i,:] = triggs_trial[np.concatenate(([0], last_trigs+1))]
        
        if plotting == 1:
        
            fig, ax = plt.subplots()
            ax.plot(trigtrace_temp)
            ax.plot(triggs_trial[np.concatenate(([0], last_trigs+1))], np.ones(len(triggs_trial[np.concatenate(([0], last_trigs+1))])), 'ro')
            
            
            
    tstrigNP_FS = np.array([tstrigNP[1], tstrigNP[2], tstrigNP[3]])
    tstrigNP_FS_all = np.zeros((3, 510)) #freqs = 51, levels = 10, 51*9*10 = 510
    
    time_before_trig = 1000
    time_after_trig = 1700000 #1698300
#getting all triggers of each sound presentation
    for i, trigger in enumerate(tstrigNP_FS):
        
        trigtrace_temp = tstrigtrace[triggs[trigger]-time_before_trig:triggs[trigger]+time_after_trig]
        
        triggs_trial = np.where(trigtrace_temp > 2000)[0]
        last_trigs = np.where(np.diff(triggs_trial)>2000)[0]
        tstrigNP_FS_all[i,:] = triggs_trial[np.concatenate(([0], last_trigs+1))]
        
        if plotting == 1:
        
            fig, ax = plt.subplots()
            ax.plot(trigtrace_temp)
            ax.plot(tstrigNP_FS_all[i,:], np.ones(len(tstrigNP_FS_all[i,:])), 'ro')

#calculating difference of synctrace for each trigger/window
    binFullPath = binFullPath_AP
    chanList = [384] #384

    # Read in metadata; returns a dictionary with string for values
    meta = readMeta(binFullPath)
    sRate_AP = SampRate(meta)
    rawData = makeMemMapRaw(binFullPath, meta)
    
    tstrigNPAP = triggs[tstrigNP]/sRate_NI

    
    sync_offset = np.zeros(len(tstrigNPAP))
    
    for i, trig_time in enumerate(tstrigNPAP):
        
        
        firstSamp_AP = int(np.round(sRate_AP*(trig_time)))
        lastSamp_AP = int(np.round(sRate_AP*(trig_time+100))) #first trig from each sound then 9 seconds, can increase to 12 13 14 15 minutes (will take ages). check how diff the first and last 10-20 sec or so
        
        firstSamp_NI = int(np.round(sRate_NI*(trig_time)))
        lastSamp_NI = int(np.round(sRate_NI*(trig_time+100)))
        
        # if binFullPath_AP == Path('S:/Fakultaet/MFZ/NWFZ/AG-deHoz-Scratch/Neuropixels/25_5_2022/m1701_25_05_22_g0/m1701_25_05_22_g0_imec0/m1701_25_05_22_g0_t0.imec0.ap.bin') and i == 23: #what is the i==23 doing?
        #     #\25_5_2022\m1701_25_05_22_g0\m1701_25_05_22_g0_imec0 #not sure what this should be
        #     firstSamp_AP = firstSamp_AP+1000
        #     firstSamp_NI = firstSamp_NI+1000
        # array of times for plot
        tDat = np.arange(0, np.size(rawData[0,firstSamp_AP:lastSamp_AP],0))
        AP_tDat = tDat/sRate_AP     # plot time axis in sec
        
        selectData = rawData[chanList, firstSamp_AP:lastSamp_AP] # All 64, should it be?
        
        # print("NI channel counts: %d, %d, %d, %d" % (MN, MA, XA, DW))
        # apply gain correction and convert to mV
        AP_convData = 1e3*GainCorrectIM(selectData, chanList, meta) # 
        
        #path to the nidaq.bin file        
        # array of times for plot
        tDat = np.arange(0, np.size(rawData_NI[0,firstSamp_NI:lastSamp_NI],0))
        NI_tDat = tDat/sRate_NI    # plot time axis in sec
        
        # binFullPath = binFullPath_NI
        # chanList = [0,2]
    
        # # Read in metadata; returns a dictionary with string for values
        # meta = readMeta(binFullPath)
        
        # selectData= rawData_NI[chanList, firstSamp:lastSamp]
        # # selectData = rawData[chanList, :]
        
        # print("NI channel counts: %d, %d, %d, %d" % (MN, MA, XA, DW))
        # apply gain correction and convert to mV
        # NI2_convData = 1e3*GainCorrectNI(selectData, chanList, meta)
        
        
        #finding up and down times of sync traces
        NI_synctrace = NI_convData[1, firstSamp_NI:lastSamp_NI]#*(np.max(AP_convData)/np.max(NI_convData[1,:]))
        AP_synctrace = AP_convData[0, :] 
        
        if len(NI_synctrace) == 0:
            print('No synctrace for analog triggers for {} file (sorted by time), skipping...')
            continue
        if len(AP_synctrace) == 0:
            print('No synctrace for spike trace for {} file (sorted by time), skipping...')
            continue
        
        NI_triggs = np.where(NI_synctrace > 0.04)[0]
        AP_triggs = np.where(AP_synctrace > np.max(AP_synctrace)*0.04/np.max(NI_synctrace))[0]
        
        NI_tstrig = np.where(np.diff(NI_triggs)>1000)[0]
        NI_tstrig_check = np.where(np.diff(NI_tstrig)< 4990)[0]
        
        if len(NI_tstrig_check) != 0:
            print('weird trigger detected in sync pulse (NI), needs cheching')
            NI_tstrig = np.delete(NI_tstrig, NI_tstrig_check[1])
        
        AP_tstrig = np.where(np.diff(AP_triggs)>1000)[0]
        AP_tstrig_check = np.where(np.diff(AP_tstrig)< 14990)[0]
        
        if len(AP_tstrig_check) != 0:
            print('weird trigger detected in sync pulse (AP), needs cheching')
            AP_tstrig = np.delete(AP_tstrig, AP_tstrig_check[1])
        
        NI_tstrig = NI_triggs[NI_tstrig]
        AP_tstrig = AP_triggs[AP_tstrig]
        
        # Plotting sync channels plus "trigger"-times
        if plotting == 1:
            fig, ax = plt.subplots()
            # ax.plot(NI_tDat, NI_convData[0, :]*(-1))
            ax.plot(NI_tDat, NI_synctrace*(np.max(AP_synctrace)/np.max(NI_synctrace)))
            ax.plot(AP_tDat, AP_synctrace)
            ax.plot(NI_tDat[NI_tstrig], np.ones(len(NI_tstrig)), 'ro')
            ax.plot(AP_tDat[AP_tstrig], np.ones(len(AP_tstrig)), 'go')
            plt.show()
        
        #sync_offset is the time that the Neuropixel channels are shifted to the back against the analog NIDAQ channel
        sync_offset[i] = np.mean(NI_tDat[NI_tstrig]-AP_tDat[AP_tstrig]) #in seconds
        
        sync_offset_range = np.max(NI_tDat[NI_tstrig]-AP_tDat[AP_tstrig])-np.min(NI_tDat[NI_tstrig]-AP_tDat[AP_tstrig])
        
        if sync_offset_range > 0.0002:
            print('philipps conscious speaking: "differences between syncwave triggers are high!", recording number: {}'.format(i))
        if sync_offset[i] > 0.02:
            print('philipps conscious speaking: "differences between syncwave triggers are high!", recording number: {}'.format(i))
                    
        
    return(tstrigNPAP, tstrigNP_FRA_all, tstrigNP_FS_all, sync_offset, sRate_AP, sRate_NI)
    
