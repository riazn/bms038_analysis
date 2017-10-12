#! /usr/bin/python
import numpy as np
import sys
import os
import random
from string import split

###################### linepick_7.7.2.py #########################
# ---------- Jennifer S. Sims, IPOP/MSKCC, Oct. 2017 ----------- #
# -------- Version accompanies Riaz et al., Cell (2017) -------- #
# --------- https://doi.org/10.1016/j.cell.2017.09.028 --------- #
################### WRITTEN FOR PYTHON 2.7.5 #####################

# NOTES: #########################################################
# - Searches all subdirectories ending in "_v7.7", the suffix given by JS_entropy_v7.7
# - Output files are named "linepick_entropy.out" and "linepick_divergence.out" -- IF YOU RE-RUN THIS PROGRAM, PREVIOUS OUTPUT FILE VERSIONS WILL BE OVERWRITTEN!


### READ *H_CL_JS.out FILES AND EXTRACT LINES ###

dirlist = []
entropies = []
divergences = []
ctrls = []

mainhead = []
divhead = []
ctrlhead = []

for folder in os.listdir('.'):
	if folder.endswith('_v7.7'):														# AS CREATED BY JS_entropy_v7.7.py
		dirlist.append(folder)
		dirname = folder
		for file in os.listdir(dirname+'/'):
			if file.endswith('H_CL_JS.out'):
				input = open(dirname+'/'+file,'r')
				for line in input.readlines():
					llist = split(line)
					item = str(llist[0])
					if item.startswith('Pt'):											# Single-repertoire data line; Modify searched prefix to match the the string the samplenames start with as needed
						entropies.append(line)
					elif item.startswith('JS_Pt'):										# Paired-reperotire comparison data line; Modify searched prefix to match the the string the samplenames start with as needed
						divergences.append(line)
					elif item.startswith('CTRL_Pt') or item.startswith('CTRL_VAR'):		# Sampling control repertoire data lines; Modify searched prefix to match the the string the samplenames start with as needed
						ctrls.append(line)
					elif item == 'FileName':
						mainhead = line
					elif item == 'JS_Divergence':
						divhead = line
					elif item == 'CTRL_Name':
						ctrlhead = line
					else:
						pass
				input.close()

print 'Divergence of',len(divergences),'sample pairs'
print 'Entropy of',len(entropies),'samples'
print len(ctrls)/2,'control repertoires'

### PRODUCE OUTPUT FILE FOR SINGLE-REPERTOIRE DATA ###

outfile1 = 'linepick_entropy.out'

output1 = open(outfile1,'w')

output1.write('%(mainhead)s' % vars()) 

for x in range(0,len(divergences)):
	y = 2*x
	Hlist1 = split(entropies[y])
	Hlist2 = split(entropies[y+1])
	Clist = split(ctrls[y])
	CVARlist = split(ctrls[y+1])
# WRITE ALL DATA ABOUT THE FIRST REAL REPERTOIRE
	pt1 = Hlist1[0]
	pt2 = Hlist1[1]
	pt3 = Hlist1[2]	
	pt4 = Hlist1[3]
	pt5 = Hlist1[4]
	pt6 = Hlist1[5]
	pt7 = Hlist1[6]
	pt8 = Hlist1[7]
	pt9 = Hlist1[8]
	pt10 = Hlist1[9]
	pt11 = Hlist1[10]
	pt12 = Hlist1[11]
	pt13 = Hlist1[12]	
	output1.write('%(pt1)s\t%(pt2)s\t%(pt3)s\t%(pt4)s\t%(pt5)s\t%(pt6)s\t%(pt7)s\t%(pt8)s\t%(pt9)s\t%(pt10)s\t%(pt11)s\t%(pt12)s\t%(pt13)s\n' % vars())
# WRITE ALL DATA ABOUT THE SECOND REAL REPERTOIRE
	pt14 = Hlist2[0]
	pt15 = Hlist2[1]
	pt16 = Hlist2[2]
	pt17 = Hlist2[3]
	pt18 = Hlist2[4]
	pt19 = Hlist2[5]
	pt20 = Hlist2[6]
	pt21 = Hlist2[7]
	pt22 = Hlist2[8]
	pt23 = Hlist2[9]	
	pt24 = Hlist2[10]
	pt25 = Hlist2[11]
	pt26 = Hlist2[12]
	output1.write('%(pt14)s\t%(pt15)s\t%(pt16)s\t%(pt17)s\t%(pt18)s\t%(pt19)s\t%(pt20)s\t%(pt21)s\t%(pt22)s\t%(pt23)s\t%(pt24)s\t%(pt25)s\t%(pt26)s\n' % vars())

# WRITE ALL AVERAGE DATA ABOUT THE SAMPLING CONTROL REPERTOIRE
	pt27 = Clist[0]
	pt28 = Clist[1]
	pt29 = Clist[2]
	pt30 = Clist[3]
	pt31 = Clist[4]
	pt32 = Clist[5]
	pt33 = Clist[6]	
	pt34 = Clist[7]
	pt35 = Clist[8]
	pt36 = Clist[9]
	pt37 = Clist[10]
	pt38 = Clist[11]
	pt39 = Clist[12]
	output1.write('%(pt27)s\t%(pt28)s\t%(pt29)s\t%(pt30)s\t%(pt31)s\t%(pt32)s\t%(pt33)s\t%(pt34)s\t%(pt35)s\t%(pt36)s\t%(pt37)s\t%(pt38)s\t%(pt39)s\n' % vars())
# WRITE ALL DATA ABOUT THE VARIANCE IN SAMPLING CONTROL REPERTOIRE
	pt40 = CVARlist[0]
	pt41 = CVARlist[1]
	pt42 = CVARlist[2]
	pt43 = CVARlist[3]
	pt44 = CVARlist[4]
	pt45 = CVARlist[5]
	pt46 = CVARlist[6]
	pt47 = CVARlist[7]
	pt48 = CVARlist[8]
	pt49 = CVARlist[9]
	pt50 = CVARlist[10]
	pt51 = CVARlist[11]
	pt52 = CVARlist[12]
	output1.write('%(pt40)s\t%(pt41)s\t%(pt42)s\t%(pt43)s\t%(pt44)s\t%(pt45)s\t%(pt46)s\t%(pt47)s\t%(pt48)s\t%(pt49)s\t%(pt50)s\t%(pt51)s\t%(pt52)s\n' % vars())
output1.close()

### PRODUCE OUTPUT FILE FOR PAIRED-REPERTOIRE COMPARISON DATA ###

outfile2 = 'linepick_divergence.out'

output2 = open(outfile2,'w')

output2.write('%(divhead)s' % vars() )

for line in divergences:
	llist = split(line)
	pt1 = llist[0]
	pt2 = llist[1]
	pt3 = llist[2]
	pt4 = llist[3]
	pt5 = llist[4]
	pt6 = llist[5]
	pt7 = llist[6]
	pt8 = llist[7]
	pt9 = llist[8]
	pt10 = llist[9]
	pt11 = llist[10]
	pt12 = llist[11]
	pt13 = llist[12]
	output2.write('%(pt1)s\t%(pt2)s\t%(pt3)s\t%(pt4)s\t%(pt5)s\t%(pt6)s\t%(pt7)s\t%(pt8)s\t%(pt9)s\t%(pt10)s\t%(pt11)s\t%(pt12)s\t%(pt13)s\n' % vars())

output2.close()

print 'Lines picked.'
