#! /usr/bin/python
import numpy as np
import sys
import os
import random
from string import split

########################## JS_aaCDR3perVJ_v2.5.py ##############################
# ----------------- Jennifer S. Sims, IPOP/MSKCC, Oct. 2017 ------------------ #
# ----- Jennifer S. Sims, Columbia University Medical Center, Feb. 2015 ------ #
# ----------- Annotated to accompany Riaz et al., Cell (2017) ---------------- #
# -------------- https://doi.org/10.1016/j.cell.2017.09.028 ------------------ #
# ----- Later version of functions described in Sims et al., PNAS (2016) ----- #
# -------------- http://www.ncbi.nlm.nih.gov/pubmed/27261081 ----------------- #
######################### WRITTEN FOR PYTHON 2.7.5 #############################

# NOTES: #######################################################################
# - This program identifies all unique VJ combinations, and calculates the number of unique aaCDR3, the entropy of aaCDR3, and the normalized entropy (evenness) of aaCDR3 for each VJ combination
# - Later versions of this program include the actual tabulation of these aaCDR3 and their respective abundances; this version does not. 
# - Iterates through all files of a specified suffix filetype in the working directory.  To make it one file at a time, extract the code between lines 46-131.
# - Remember, the entropy of 1 is 0!!!

###################### BEGIN FUNCTIONS ########################

def entropy(list1):
	pvalues = np.array(list1).astype(float)
	px = pvalues/pvalues.sum()
	px = px[px.nonzero()]
	Hp = sum(-px*np.log2(px))
	Hpmax = np.log2(float(len(list1)))

	return Hp,Hpmax,Hp/Hpmax,len(pvalues)

####################### END FUNCTIONS #########################


filetype = sys.argv[1]			# This is flexible in order to eliminate however much of the filename you wish.  For example, use "only.productive.tsv" if applicable because you used JS_splitter_Imsqv4.5.py.

filelist = []					# list of files to process
filestrip = []					# list of files, with suffix stripped -- this will be concatenated with "_aaCDR3perVJ.tsv" to create the outfile name

for file in os.listdir('.'):
	if file.endswith(filetype):
		filelist.append(file)

dFILE = {}						# dictionary for retaining the line-item data on each VJ combo from each file.  Not necessary in this version -- used in later versions to combine data for a given VJcombo as it occurs in different repertoire files.
for infile in filelist:
	instrip = infile.rstrip(filetype)	# a two-step stripping process for filenames that were already using the same delimiters...
	instrip = instrip.rstrip('._')		#
	filestrip.append(instrip)

	############## READ IN THE CLONOTYPES FROM THE FILE ################
	input = open(infile,'r')
	dTOT = {}					# dictionary for totCDR3s and their abundances
	dVJ = {}					# dictionary for VJcombos and their abundances
	k1 = 0 						# iterable for lines processed
	m1 = 0 						# iterable for unique totCDR3s included
	ct1 = 0 					# iterable for reads added
	for line in input.readlines():
		if k1 != 0:
			llist = split(line)
			aa = llist[4]
			combo = llist[0]
			reads = int(llist[1])
			ct1 = ct1 + reads
			try:
				aalist = dVJ[combo]
				if aa not in aalist:		# addition of unique aa
					aalist.append(aa)
				else:
					pass
				dVJ[combo] = aalist
			except KeyError:				# addition of new VJcombo and initialization of aa list
				aalist = []
				aalist.append(aa)
				dVJ[combo] = aalist
			totid = aa,combo
			try:
				dTOT[totid] = dTOT[totid] + reads
				m1 = m1
			except KeyError:
				dTOT[totid] = reads
				m1 = m1+1
		k1 = k1+1
	input.close()

	########## FOR EACH VJcombo, PROCESS THE aaCDR3 STATISTICS ############
	aa_data = []				# list into which the VJcombo and its associate aaCDR3 repertoire data will be indexed
	for combo1 in dVJ.keys():
		aalist1 = dVJ[combo1]
		aa_reads = []
		for aa1 in aalist1:
			totid1 = aa1,combo1
			try:
				reads = int(dTOT[totid1])
				aa_reads.append(reads)
			except KeyError:
				print totid1
				pass
			except TypeError:
				print totid1
		ans = entropy(aa_reads)				# ENTROPY FUNCTION HERE
		aa_data.append([combo1,len(aalist1),ans[0],ans[2],ans[1]])
	temp1 = []	#len
	temp2 = []	#h
	temp3 = []	#h/hmax
	temp4 = []	#hmax
	for y in range(0,len(aa_data)):
		data = aa_data[y]
		temp1.append(data[1])
		temp2.append(data[2])
#		if data[3] > 0:
		temp3.append(float(data[3]))
		temp4.append(data[4])
	dFILE[instrip] = aa_data,np.median(temp1),np.mean(temp2)
	
	################## PRINT THE OUTFILE ######################

	outfile1 = str(instrip) + '_aaCDR3perVJ.tsv'

	output1 = open(outfile1,'w')
	h1 = 'VJcombo'
	h2 = 'Number_aaCDR3'
	h3 = 'H_aaCDR3'
	h4 = 'Hnorm_aaCDR3'
	output1.write('%(h1)s\t%(h2)s\t%(h3)s\t%(h4)s\n' % vars())	

	for i in range(0,len(aa_data)):
		combodata = aa_data[i]
		pt1 = combodata[0]
		pt2 = combodata[1]
		pt3 = combodata[2]
		if str(combodata[3]).find('nan') != -1:
			pt4 = 0
		else:
			pt4 = combodata[3]
		output1.write('%(pt1)s\t%(pt2)s\t%(pt3)s\t%(pt4)s\n' % vars())
	output1.close()

	##### THIS IS A FRIVOLOUS FUNCTION BUT IT'S A LOT MORE FUN THAN PRINTING 'DONE' ALL THE TIME #####
	batmansays = ['BOOM!','KAPOW!','SHAZZAM!','WHAMMY!','VRONK!','SPLAT!','BANG!','WHAP!','ZOWIE!','SPLAT!','BAM!(just kidding)','CLANK!(just kidding)','AWK!(just kidding)']
	print instrip,'Counts:',ct1,'VJcombos:',len(aa_data),'...'+str(random.choice(batmansays))

