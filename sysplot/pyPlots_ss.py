import os
import re
import numpy as np
import matplotlib.pyplot as plt

TiologPath = r'X:\Case and FA\Phoenix\MBUSE-1241\tiotest\AHIT12ms\256G'

os.chdir(TiologPath)
loglist = [log for log in os.listdir(TiologPath) if os.path.isfile(log) and 'swsr' in log]
# print(loglist)

def ParseTiologName(TiologName):
    reStr_swsrfilename = r'.*ss_swsr_tl(\d+)kb_(\S+)_(\d+).txt'
    match_swsrfile = re.match(reStr_swsrfilename,TiologName)
    if match_swsrfile:
        swsrTL = int(match_swsrfile.group(1))
        swsrOpt = match_swsrfile.group(2)
        swsrRound = int(match_swsrfile.group(3))
    else:
        print('Tiolog name %s is not expected!'% TiologName)
    return swsrRound, swsrTL, swsrOpt

lstSWperf = [[] for row in range(len(loglist))]
lstSRperf = [[] for row in range(len(loglist))]

def ParsePerf(OrderNum, Tiolog):
    reStr_TioSWperf = r'\|\s+Write\s+.*\|\s+(\d+.\d+)\s+MB/s'
    reStr_TioSRperf = r'\|\s+Read\s+.*\|\s+(\d+.\d+)\s+MB/s'

    with open(Tiolog,newline='') as file:
        lines = file.readlines()

        for i in range(len(lines)):
            match_SWperf  = re.match(reStr_TioSWperf, lines[i])
            match_SRperf  = re.match(reStr_TioSRperf, lines[i])
            if match_SWperf:
                global lstSWperf
                lstSWperf[OrderNum].append(float(match_SWperf.group(1)))
            if match_SRperf:
                global lstSRperf
                lstSRperf[OrderNum].append(float(match_SRperf.group(1)))
        #print('SWperf rows:\t%d' % len(lstSWperf[OrderNum]))
        #print('SRperf rows:\t%d' % len(lstSRperf[OrderNum]))
        return ParseTiologName(Tiolog)

plt.figure(figsize=(6, 7))
plt.title('CS FW 0001 256G Seq Write Read 1G Cumulative test\n AHIT12ms on YJ8250')
plt.xlim(1, 218)
plt.ylim(0, 1800)
plt.xlabel('GB')
plt.ylabel('MBps')

lstMarker = ['.', 'o', '^', 's']
lstColor = ['gold', 'purple', 'lightgreen', 'royalblue' ]


i = 0
for log in loglist:
    #print('\n')
    #print(log)
    swsrRound, swsrTL, swsrOpt = ParsePerf(i, log)
    #print(lstSWperf[i])
    Pmarker = lstMarker[swsrRound]
    Pcolor = lstColor[int(i/3)]
    plt.plot(lstSWperf[i], c=Pcolor, linewidth=1, alpha=0.2, marker=Pmarker, label='SW%d_TL%dKB_%s'%(swsrRound,swsrTL,swsrOpt))
    plt.plot(lstSRperf[i], c=Pcolor, linewidth=1, alpha=0.2, marker=Pmarker, label='SR%d_TL%dKB_%s'%(swsrRound,swsrTL,swsrOpt))
    plt.legend(bbox_to_anchor=(1.02,0), loc=3, borderaxespad=1)
    i = i+1

import re
import numpy as np
import matplotlib.pyplot as plt

def ParseFtrace(ftracelog, MB):

    reStr_Write10_ufshcdcmd = r'.*(\d+.\d+): ufshcd_command:.*scsi_send.*cmd: 0x2a lba: (\d+)\s+size:\s+(\d+)'
    reStr_Read10_ufshcdcmd = r'.*(\d+.\d+): ufshcd_command:.*scsi_send.*cmd: 0x28 lba: (\d+)\s+size:\s+(\d+)'
    reStr_Unmap_blockrq = r'.*\s+(\d+.\d+):\s+block_rq_complete:\s+8,0\s+D\s+\(\)\s+(\d+)\s+\+\s+(\d+)\s+\[0'

    with open(ftracelog,newline='') as file:
        lines = file.readlines(1024*1024*MB)
        #print('Total row: {}'.format(len(lines)))
        #print('First row: {}'.format(lines[0]))
        #print('Last  row: {}'.format(lines[-1]))
        lstWrite10StartLba = []
        lstWrite10TL = []
        lstWrite10ExtLba = []
        lstWrite10Num = []
        lstRead10StartLba = []
        lstRead10TL = []
        lstRead10ExtLba = []
        lstRead10Num = []
        lstUnmapStartLba = []
        lstUnmapTL = []
        lstUnmapExtLba = []
        lstUnmapNum = []
     
        for i in range(len(lines)):
            match_Write10_ufshcdcmd = re.match(reStr_Write10_ufshcdcmd, lines[i])
            match_Read10_ufshcdcmd = re.match(reStr_Read10_ufshcdcmd, lines[i])
            match_Unmap_blockrq = re.match(reStr_Unmap_blockrq, lines[i])
            if match_Write10_ufshcdcmd:
                lstWrite10StartLba.append(int(int(match_Write10_ufshcdcmd.group(2))/8))
                lstWrite10TL.append(int(int(match_Write10_ufshcdcmd.group(3))/4096))
                lstWrite10Num.append(i)
                if lstWrite10TL[-1] > 1:
                    lstWrite10ExtLba.append([i, lstWrite10StartLba[-1]+lstWrite10TL[-1]-1])
            if match_Read10_ufshcdcmd:
                lstRead10StartLba.append(int(int(match_Read10_ufshcdcmd.group(2))/8))
                lstRead10TL.append(int(int(match_Read10_ufshcdcmd.group(3))/4096))
                lstRead10Num.append(i)
                if lstRead10TL[-1] > 1:
                    lstRead10ExtLba.append([i, lstRead10StartLba[-1]+lstRead10TL[-1]-1])       
            if match_Unmap_blockrq:
                lstUnmapStartLba.append(int(int(match_Unmap_blockrq.group(2))/8))
                lstUnmapTL.append(int(int(match_Unmap_blockrq.group(3))/8))
                lstUnmapNum.append(i)
                if lstUnmapTL[-1] > 1:
                    lstUnmapExtLba.append([i, lstUnmapStartLba[-1]+lstUnmapTL[-1]-1])
        print('\n%s' % ftracelog)
        print('Write10 commands:\t%d' % len(lstWrite10StartLba))
        print('Read10 commands:\t%d' % len(lstRead10StartLba))
        print('Unmap commands:\t\t%d' % len(lstUnmapStartLba))
        return [lstWrite10Num, lstWrite10StartLba, lstWrite10TL, lstWrite10ExtLba, lstRead10Num, lstRead10StartLba, lstRead10TL, lstRead10ExtLba, lstUnmapNum, lstUnmapStartLba, lstUnmapTL, lstUnmapExtLba]


def lst2npXY(Values2D, outputX=True, outputY=True):
   lstX = []
   lstY = []
   for val in Values2D:
       lstX.append(val[0])
       lstY.append(val[1])
   if outputX and outputY:
       return np.array(lstX), np.array(lstY)
   if outputX:
       return np.array(lstX)
   if outputY:
       return np.array(lstY)

ftracelog = r'X:\Case and FA\Phoenix\MBUSE-1241\ftrace\ParsedByGB\TraceTrace_Run1_GB1.txt'
#ftracePath = r'X:\Case and FA\Phoenix\MBUSE-1241\ftrace\ParsedByGB'
ftracePath = r'X:\Case and FA\Phoenix\MBUSE-1241\tiotest\CSL03\_n_outputfile'

os.chdir(ftracePath)
loglist = [log for log in os.listdir(ftracePath) if os.path.isfile(log) and '_' in log][:1]
print(loglist)


glstWrite10StartLba = []
glstWrite10TL = []
glstWrite10ExtLba = []
glstWrite10Num = []
glstRead10StartLba = []
glstRead10TL = []
glstRead10ExtLba = []
glstRead10Num = []
glstUnmapStartLba = []
glstUnmapTL = []
glstUnmapExtLba = []
glstUnmapNum = []

for log in loglist:
    lstFtrace = ParseFtrace(log, 50)

    glstWrite10Num.extend(lstFtrace[0])
    glstWrite10StartLba.extend(lstFtrace[1])
    glstWrite10TL.extend(lstFtrace[2])
    glstWrite10ExtLba.extend(lstFtrace[3])
    glstRead10Num.extend(lstFtrace[4])
    glstRead10StartLba.extend(lstFtrace[5])
    glstRead10TL.extend(lstFtrace[6])
    glstRead10ExtLba.extend(lstFtrace[7])
    glstUnmapNum.extend(lstFtrace[8])
    glstUnmapStartLba.extend(lstFtrace[9])
    glstUnmapTL.extend(lstFtrace[10])
    glstUnmapExtLba.extend(lstFtrace[11])


#plt.figure(figsize=(20, 15))
plt.title('LBA Dist.')
#plt.ylim(1e7,2e7)
plt.xlabel('Commands')
plt.ylabel('LBA')


plt.scatter(glstWrite10Num, glstWrite10StartLba, c='r', s=1, alpha=0.05, label='Write10')
plt.scatter(glstRead10Num, glstRead10StartLba, c='b', s=1, alpha=0.05, label='Read10')
plt.scatter(glstUnmapNum, glstUnmapStartLba, c='y', s=1, alpha=0.05, label='Unmap')

npXval, npYval = lst2npXY(glstWrite10ExtLba)
plt.scatter(npXval, npYval, c='r', s=1, alpha=0.0, marker='+')
npXval, npYval = lst2npXY(glstRead10ExtLba)
plt.scatter(npXval, npYval, c='b', s=1, alpha=0.0, marker='+')
npXval, npYval = lst2npXY(glstUnmapExtLba)
plt.scatter(npXval, npYval, c='y', s=1, alpha=0.0, marker='+')

plt.legend(bbox_to_anchor=(1.02,0), loc=3, borderaxespad=1)

