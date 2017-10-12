#! /usr/bin/python
import numpy as np
import sys
from string import split
import os

################## JS_splitter_Imsqv4.5.py ######################
# ---------- Jennifer S. Sims, IPOP/MSKCC, Oct. 2017 ----------- #
# -------- Version accompanies Riaz et al., Cell (2017) -------- #
# --------- https://doi.org/10.1016/j.cell.2017.09.028 --------- #
################# WRITTEN FOR PYTHON 2.7.5 #######################


splittype = sys.argv[1]		#FILTER TYPE:  "only" or "all"
runfile = sys.argv[2]		# filename -- should be an Immunoseq v3.0 *.tsv exported from the Adaptive portal
	
runstrip = runfile.rstrip('.tsv')		# removing the suffix to obtain the sample name

# Dictionaries of the repertoire at different resolutions:  keys = clonotypes, values = reads
dCDR3 = {}
dVJ = {}
dTOT = {}

# Dictionary in which keys = the whole TCR identifier, values = reads
dLINE = {}

ct = 0			# iterable for cumulative reads

k = 0			# iterable for lines read
i = 0			# iterable for aaCDR3 added
l = 0			# iterable for VJcombos added
m = 0			# iterable for totCDR3 added
n = 0			# iterable for lines added
input = open(runfile,'r')
for line in input.readlines():
	if k != 0:
		llist = split(line,'\t')
		reads = int(llist[2])
		ct = ct + reads
		ntseq = str(llist[0])

# LOGIC FOR DISCERNING V AND J CASSETTES FROM ADAPTIVE *.TSV TO THE RESOLUTION OF GENE NAME IS AS INCLUSIVE AS POSSIBLE, TRYING FIRST "GeneName" MATCHES, THEN "GeneNameTies", THEN "FamilyName", THEN "FamilyTies", SUCH THAT ANY POSITIVE MAPPING RETAINS THE READS
		v1 = str(llist[7])									# Immunoseq Col7 = "vGeneName", first choice for Vcassette ID
		if v1.startswith('TCR'):
			vname = v1
		elif v1 == 'unresolved':
			vcass1 = split(llist[10],',')					# Immunoseq Col10 = "vGeneNameTies", a list that contains potential matches
			v2 = str(vcass1[0].strip('"'))					# First hit in list
			v2 = str(v2.rstrip('"'))
			if v2.startswith('TCR'):
				vname=v2
			elif v2 == 'unresolved':
				v3 = str(llist[6].strip('"'))				# Immunoseq Col6 = "vFamilyName", in the absence of a positive gene ID or potential ties, use the family name
				v3 = str(v3.rstrip('"'))
				if v3.startswith('TCR'):
					vname = v3
				elif v3 == 'unresolved':
					vcass2 = split(llist[9],',')			# Immunoseq Col9 = "vFamilyTies", in the absence of all of the above, a list that contains potential families
					v4 = str(vcass2[0],',').strip('"')		# First hit in list
					v4 = str(v4.rstrip('"'))
					vname = v4
		j1 = str(llist[21])									# Immunoseq Col21 =  "jGeneName", first choice for Jcassette ID
		if j1.startswith('TCR'):
			jname = j1
		elif j1 == 'unresolved':
			jcass1 = split(llist[24],',')					# Immunoseq Col24 = "jGeneNameTies", a list that contains potential matches
			j2 = str(jcass1[0].strip('"'))					# First hit in list
			j2 = str(j2.rstrip('"'))
			if j2.startswith('TCR'):
				jname=j2
			elif j2 == 'unresolved':
				j3 = str(llist[20].strip('"'))				# Immunoseq Col20 = "jFamilyName", in the absence of a positive gene ID or potential ties, use the family name
				j3 = str(j3.rstrip('"'))
				if j3.startswith('TCR'):
					jname = j3
				elif j3 == 'unresolved':
					jcass2 = split(llist[23],',')			# Immunoseq Col23 = "jFamilyTies", in the absence of all of the above, a list that contains potential families
					j4 = str(jcass2[0],',').strip('"')		# First hit in list
					j4 = str(j4.rstrip('"'))
					jname = j4
		combo = vname + '.' + jname							# VJcombo name is the join of the chosen designation for the Vcassette and Jcassette
		aa = llist[1]
		if aa == '':										# If not amino acid translation, DO NOT ADD ANY OF THIS TO ANY OF THE DICTIONARIES
			pass

# ADD THESE IDENTIFIERS TO THE DICIONARIES
		else:
			try:
				dCDR3[aa] = dCDR3[aa] + reads
				i = i
			except KeyError:
				dCDR3[aa] = reads
				i = i+1
			try:
				dVJ[combo] = dVJ[combo] + reads
				l = l
			except KeyError:
				dVJ[combo] = reads
				l = l+1
			totid = aa,combo
			try:
				dTOT[totid] = dTOT[totid] + reads
				m = m
			except KeyError:
				dTOT[totid] = reads
				m = m+1
			linedata = combo,vname,jname,aa,ntseq			# concatenated whole-line identifier
			try:
				dLINE[linedata] = dLINE[linedata] + reads
				n = n
			except KeyError:
				dLINE[linedata] = reads
				n=n+1
	k = k+1
input.close()
print k
print n
print 'totCDR3',m
print 'VJ',l
print 'aaCDR3',i

# ----------- aaCDR3 SPLIT OUTPUT FILES ------------------------------

dirname1 = 'aaCDR3'
cmd1 = 'mkdir %(dirname1)s' % vars()
os.system(cmd1)

if splittype == 'only':
	outfile1 = 'aaCDR3/' + runstrip + '_only.aaCDR3.out'

	output1 = open(outfile1,'w')
	h1 = 'aaCDR3_productive'
	h2 = runstrip
	output1.write('%(h1)s\t%(h2)s\n' % vars())
	for aa in dCDR3.keys():
		pt1 = aa
		if pt1.startswith('C'):
			if pt1.find('*') == -1:
				if pt1 != '':
					try:
						pt2 = int(dCDR3[aa])
					except KeyError:
						pt2 = ''
					output1.write('%(pt1)s\t%(pt2)s\n' % vars())
				else:
					pass
			else:
				pass
		else:
			pass
	output1.close()

elif splittype == 'all':
	outfile1 = runstrip + '_all.aaCDR3.out'

	output1 = open(outfile1,'w')
	h1 = 'aaCDR3'
	h2 = runstrip
	output1.write('%(h1)s\t%(h2)s\n' % vars())
	for aa in dCDR3.keys():
		pt1 = aa
		try:
			pt2 = int(dCDR3[aa])
		except KeyError:
			pt2 = ''
		output1.write('%(pt1)s\t%(pt2)s\n' % vars())
	output1.close()

#---------- VJ SPLIT OUTPUT FILE ------------------------------------------

dirname2 = 'VJ'
cmd2 = 'mkdir %(dirname2)s' % vars()
os.system(cmd2)

outfile2 = 'VJ/' + runstrip + '_VJ.out'

output2 = open(outfile2,'w')
h1 = 'VJcombo'
h2 = runstrip
output2.write('%(h1)s\t%(h2)s\n' % vars())

for combo in dVJ.keys():
	pt1 = combo
	try:
		pt2 = int(dVJ[combo])
	except KeyError:
		pt2 = ''
	output2.write('%(pt1)s\t%(pt2)s\n' % vars())
output2.close()

#---------- totCDR3 SPLIT OUTPUT FILES -------------------------------------------

dirname3 = 'totCDR3'
cmd3 = 'mkdir %(dirname3)s' % vars()
os.system(cmd3)

if splittype == 'only':

	outfile3 = 'totCDR3/' + runstrip + '_only.totCDR3.out'

	output3 = open(outfile3,'w')
	h1 = 'aaCDR3_productive'
	h2 = 'VJcombo'
	h3 = runstrip
	output3.write('%(h1)s\t%(h2)s\t%(h3)s\n' % vars())
	for totid in dTOT.keys():
		pt1 = totid[0]
		pt2 = totid[1]
		if pt1.startswith('C'):
			if pt1.find('*') == -1:
				if pt1 != '':
					try:
						pt3 = int(dTOT[totid])
					except KeyError:
						pt3 = ''
					output3.write('%(pt1)s\t%(pt2)s\t%(pt3)s\n' % vars())
				else:
					pass
			else:
				pass
		else:
			pass
	output3.close()

elif splittype == 'all':
	outfile3 = 'totCDR3/' + runstrip + '_all.totCDR3.out'

	output3 = open(outfile3,'w')
	h1 = 'aaCDR3'
	h2 = 'VJcombo'
	h3 = runstrip
	output3.write('%(h1)s\t%(h2)s\t%(h3)s\n' % vars())
	for totid in dTOT.keys():
		pt1 = totid[0]
		pt2 = totid[1]
		try:
			pt3 = int(dTOT[totid])
		except KeyError:
			pt3 = ''
		output3.write('%(pt1)s\t%(pt2)s\t%(pt3)s\n' % vars())
	output3.close()

# --------------- RECONSTITUTED WHOLE-REPERTOIRE FILES --------------------------
# --------------- No filter

if splittype == 'all':

	outfile4 = runstrip + '_all.productive.tsv'

	output4 = open(outfile4,'w')
	h1 = 'VJcombo' 
	h2 = 'Counts'
	h3 = 'Vcassette'
	h4 = 'Jcassette'
	h5 = 'aaCDR3'
	h6 = 'ntCDR3'

	output4.write('%(h1)s\t%(h2)s\t%(h3)s\t%(h4)s\t%(h5)s\t%(h6)s\n' % vars())
	for linedata in dLINE.keys():
		pt1 = str(linedata[0])
		pt2 = int(dLINE[linedata])
		pt3 = str(linedata[1])
		pt4 = str(linedata[2])
		pt5 = str(linedata[3])
		pt6 = str(linedata[4])
		output4.write('%(pt1)s\t%(pt2)s\t%(pt3)s\t%(pt4)s\t%(pt5)s\t%(pt6)s\n' % vars())
	output4.close()

# --------------- "only" filter

elif splittype == 'only':

	outfile5 = runstrip + '_only.productive.tsv'

	output5 = open(outfile5,'w')
	h1 = 'VJcombo' 
	h2 = 'Counts'
	h3 = 'Vcassette'
	h4 = 'Jcassette'
	h5 = 'aaCDR3_filtered'
	h6 = 'ntCDR3'

	output5.write('%(h1)s\t%(h2)s\t%(h3)s\t%(h4)s\t%(h5)s\t%(h6)s\n' % vars())
	for linedata in dLINE.keys():
		pt1 = str(linedata[0])
		pt2 = int(dLINE[linedata])
		pt3 = str(linedata[1])
		pt4 = str(linedata[2])
		pt5 = str(linedata[3])
		pt6 = str(linedata[4])
		if pt5.startswith('C'):
			if pt5.find('*') == -1:
				if pt5 != '':
					output5.write('%(pt1)s\t%(pt2)s\t%(pt3)s\t%(pt4)s\t%(pt5)s\t%(pt6)s\n' % vars())
				else:
					pass
			else:
				pass
		else:
			pass
	output5.close()

print 'Done'