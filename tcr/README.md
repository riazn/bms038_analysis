The set of Python code described in <A HREF="https://doi.org/10.1016/j.cell.2017.09.028">Riaz et al. (2017)</A> performs three basic operations on the tabulated repertoires – containing unique CDR3 clonotypes and their respective abundances – provided as output by the Immunoseq v2.0 (Adaptive Biotechnologies) proprietary processing pipeline (compressed as <A HREF="https://github.com/riazn/bms038_analysis/tree/master/data/tcr/">/data/tcr/Immunoseq.tgz</A>):
1)	Filtration of unproductive CDR3 sequences
2)	Parsing of CDR3 abundance data to the resolutions of:  V-J combinations (VJ), amino acid CDR3 motifs (aaCDR3), and VJ/aaCDR3-unique clonotypes (totCDR3)
3)	Calculation of entropy and other diversity statistics on each repertoire as described at each of these resolutions

The output of these functions are provided in a Microsoft Excel sheet (<A HREF="https://github.com/riazn/bms038_analysis/tree/master/tcr/TCRmetrics_Riaz2017.xlsx">TCRmetrics_Riaz2017.xlsx</A>) and its .csv version <A HREF="https://github.com/riazn/bms038_analysis/tree/master/data/tcr/TCRmetrics.csv">TCRmetrics.csv</A>.

An additional R script (<A HREF="https://github.com/riazn/bms038_analysis/tree/master/tcr/createTCRfigures.r">"createTCRfigures.R"</A> is provided to automate the graphical output in the manuscript figures.

<h1>Python functions</h1>

The pipeline employs three primary Python programs:
<P><tt><B>JS_splitter_Imsqv4.5.py</B> [filtertype] [filename] </tt>– parses Immunoseq (*.tsv) output containing 46 columns into a repertoire tabulation and resolutions (output files <tt>[samplename]_VJ.out</tt>, <tt>[samplename]_aaCDR3.out</tt>, <tt>[samplename]_totCDR3.out</tt>, <tt>[samplename]_ntCDR3.out</tt>).

<P><tt><B>JS_entropy_v7.7.py</B> [filename1] [filename2] </tt>– performs entropy calculations and entropy comparisons on the “complete” repertoire files pairwise, to produce output file <tt>[samplename]_H_CL_JS.out</tt>.  This code is functionally identical to that by the same name in https://github.com/ShenLab/Repertoire, as published by <A HREF="https://www.ncbi.nlm.nih.gov/pubmed/27261081">Sims et al. (2016)</A>.

<P><tt><B>JS_aaCDR3perVJ_v2.5.py</B> [filename] </tt>– parses the amino acid CDR3 motifs encoded by each V-J combination, and outputs the number, entropy (H), and normalized entropy (Hnorm) of the aaCDR3 repertoire encoded by each V-J combination.

and one file handling script:
<P><tt><B>linepick_v7.7.2.py</B> </tt>– which simply extracts single-sample statistics and paired-sample statistics from the output of JS_entropy_v7.7.py and assembles them into summary files for groupwise analysis.  Current version is filename dependent and can be replaced by a more flexible input function, or shell scripting.
  
Additionally, we used:
<P><tt><B>multiple_joint_kde.py</B></tt>, available through the <A HREF="http://seaborn.pydata.org/">Seaborn</A> visualization library for python (<A HREF="https://github.com/mwaskom/seaborn">https://github.com/mwaskom/seaborn</A>).

<h1>Implementation</h1>
Place the python scripts in the working directory containing tabulated repertoire files, such as those provided for export by Adaptive for the Immunoseq sequencing analysis. These versions are formatted to accommodate the columns associated with their v3.0 export(46 columns).
<BR>
<P><h2>I. <tt>JS_splitter_Imsqv4.5.py [filtertype] [filename]</tt></h2>
<P><h3>Input:</h3>
<P><tt>filtertype </tt>– The filter “only” removes any CDR3 sequences which:
A)	contain stop codons (“*”) in the amino acid translation
B)	do not contain an amino acid translation
C)	contain ambiguous nucleotides (“N”, “R”, “W”, “Y”) in the nucleotide sequence
D)	do not contain mapped V or J cassettes
The filter “all” carries such CDR3s through into the parsed repertoires, including those reads among the listed clonotypes, and in statistical calculations on the repertoire.  The “only” option was used to produce the repertoires analyzed in <A HREF="https://doi.org/10.1016/j.cell.2017.09.028">Riaz et al. (2017)</A>.
filename – adaptive-exported repertoire file
Further input instructions are in the comments on the input argument lines.

<P><h3>Output:</h3>
<P>Program strips the input filename of its suffix to provide the "samplename" base for the output filename, producing <tt>[samplename]_only.productive.tsv</tt> and directories titled <tt>aaCDR3</tt>, <tt>totCDR3</tt>, and <tt>VJ</tt>, to which thusly parsed repertoires from all samples will be written as tab-delimited text files (*.out).

In <tt>[samplename]_only.productive.tsv</tt>, data columns are as follows:
<UL>
<LI>Col0:  VJcombo – V-J cassette combination derived from Col2 and Col3
<LI>Col1:  Counts – Number of reads
<LI>Col2:  Vcassette – Identity of mapped V cassette
<LI>Col3:  Jcassette – Identity of mapped J cassette
<LI>Col4:  aaCDR3 – amino acid translation of CDR3 motif
<LI>Col5:  ntCDR3 - nucleotides comprising CDR3 motif
</UL>

<h3>Example:</h3>
<TABLE><TR><TD>
<tt>python JS_splitter_Imsqv4.5.py only Pt3_pre.tsv</tt>
<h3>.<BR>.<BR>.</h3>
<tt>Pt3_pre_only.productive.tsv</tt><BR>
<tt>VJ/Pt3_pre_VJ.out</tt><BR>
<tt>aaCDR3/Pt3_pre_aaCDR3.out</tt><BR>
<tt>totCDR3/Pt3_pre_totCDR3.out</tt><BR>
</TD></TR></TABLE

<P><h2>II. <tt>JS_entropy_v7.7.py [filename1] [filename2]</tt></h2>
<P><h3>Input:</h3>
After running JS_splitter_Imsqv4.5.py, in the working directory, each input sample file is represented by <tt>[samplename]_only.productive.tsv</tt>.  For each pair of samples, <tt>filename1</tt> and <tt>filename2</tt> are therefore now represented by <tt>filename1 = [samplename1]_only.productive.tsv</tt> and <tt>filename2 = [samplename1]_only.productive.tsv</tt>.

<P><h3>Output:</h3>
A new directory is created and named <tt>[samplename1]_[samplename2]</tt>.  In this directory, paired repertoire files (frequencies of each TCR in each of the two samples) will be produced for each resolution, and a file of individual and comparative entropy-based statistics, including a sample size control for the Jensen-Shannon Divergence Metric (JSM) which samples repertoire 1 to the number of clonotypes as repertoire 2.  For more details on this operation, see <A HREF="https://www.ncbi.nlm.nih.gov/pubmed/27261081">Sims et al. (2016)</A>.  Note that Evenness = 1-Clonality = Hnorm = H/Hmax, where H is entropy, Hmax (maximum entropy) = log2(N), therefore Hnorm signifies normalized entropy.

<P>The output file <tt>[samplename1]_[samplename2]_H_CL_JS.out</tt> contains:
<P>Line0: headers for single-repertoire output
<P>Line1 and Line2:
<UL>
<LI>Col0 = samplename1 or samplename2
<LI>Col1 = Hcdr3 (aaCDR3 entropy)
<LI>Col2 = Hvj (VJ entropy)
<LI>Col3 = Htot (totCDR3 entropy)
<LI>Col4 = CLcdr3 (aaCDR3 clonality)
<LI>Col5 = CLvj (VJ clonality)
<LI>Col6 = CLtot (totCDR3 clonality)
<LI>Col7 = Hcdr3_max (maximum aaCDR3 entropy; log2(N) where N is aaCDR3 population size)
<LI>Col8 = Hvj_max (maximum VJ entropy; log2(N) where N is VJ population size)
<LI>Col9 = Htot_max (maximum totCDR3 entropy; log2(N) where N is totCDR3 population size)
<LI>Col10 = Num_CDR3 (number of unique aaCDR3)
<LI>Col11 = Num_VJ (number of unique VJ combinations)
<LI>Col12 = Num_totCDR3 (number of unique totCDR3 clonotypes)
</UL>
<P>Line3: headers for paired-repertoire output
<P>Line4:
<UL>
<LI>Col0 = “[samplename1]_[samplename2]”
<LI>Col1 = JScdr3 (Jensen-Shannon divergence of aaCDR3 repertoires)
<LI>Col2 = JSvj (Jensen-Shannon divergence of VJ repertoires)
<LI>Col3 = JStotCDR3 (Jensen-Shannon divergence of totCDR3 repertoires)
<LI>Col4 = cdr3_COMBINEDnum (total number of unique aaCDR3s in both repertoires)
<LI>Col5 = vj_COMBINEDnum (total number of unique VJ combinations in both repertoires)
<LI>Col6 = tot_COMBINEDnum (total number of unique totCDR3s in both repertoires)
<LI>Col7 = JScdr3_ctrl (aaCDR3 Jensen-Shannon divergence of samplename1 and samplename1’)*
<LI>Col8 = JSvj_ctrl (VJ Jensen-Shannon divergence of samplename1 and samplename1’)*
<LI>Col9 = JStot_ctrl (totCDR3 Jensen-Shannon divergence of samplename1 and samplename1’)*
<LI>Col10 = JScdr3_ctrlVAR (variance in aaCDR3 Jensen-Shannon divergence of samplename1 and samplename1’)*
<LI>Col11 = JSvj_ctrlVAR (variance in VJ Jensen-Shannon divergence of samplename1 and samplename1’)*
<LI>Col12 = JStot_ctrlVAR (variance in totCDR3 Jensen-Shannon divergence of samplename1 and samplename1’)*
</UL>
* Sampling control (samplename1’) is created by downsampling the specified repertoire {aaCDR3, VJ, totCDR3} of samplename1 to the number of unique members in the equivalent repertoire of samplename2.  Divergence calculations are then performed as JSD(samplename1,samplename1’).  This process is performed 10 times, with the results in Line4,Col7-9 representing the mean value over these iterations, and Line4,Col10-12 representing the variance in those values.

<P>Line5: headers for single-repertoire statistics of sampling control (samplename1’)<BR>
<P>Line6:
<UL>
<LI>Col0 = “CTRL_[samplename1]_[samplename2]”
<LI>Col1 = Hcdr3_ctrl (entropy of samplename1aaCDR3’)
<LI>Col2 = Hvj_ctrl (entropy of samplename1VJ’)
<LI>Col3 = Htot_ctrl (entropy of samplename1totCDR3’)
<LI>Col4 = CLcdr3_ctrl (clonality of samplename1aaCDR3’)
<LI>Col5 = CLvj_ctrl (clonality of samplename1VJ’)
<LI>Col6 = CLtot_ctrl (clonality of samplename1totCDR3’)
<LI>Col7 = Hcdr3_max_ctrl (maximum entropy of samplename1aaCDR3’)
<LI>Col8 = Hvj_max_ctrl (maximum entropy of samplename1VJ’)
<LI>Col9 = Htot_max_ctrl (maximum entropy of samplename1totCDR3’)
<LI>Col10 = Num_CDR3_ctrl (number of unique aaCDR3 in samplename1aaCDR3’)
<LI>Col11 = Num_VJ_ctrl (number of unique VJ in samplename1VJ’)
<LI>Col12 = Num_totCDR3_ctrl (number of unique totCDR3 in samplename1totCDR3’)
</UL>

<h3>Example:</h3>
<TABLE><TR><TD>
<tt>python JS_entropy_v7.7.py Pt3_pre.tsv Pt3_on.tsv</tt>
<h3>.<BR>.<BR>.</h3>
<tt>Pt3_pre_Pt3_on/Pt3_pre_Pt3_on_H_CL_JS.out</tt><BR>
<tt>Pt3_pre_Pt3_on/Pt3_pre_Pt3_on_cdr3.out</tt><BR>
<tt>Pt3_pre_Pt3_on/Pt3_pre_Pt3_on_VJ.out</tt><BR>
<tt>Pt3_pre_Pt3_on/Pt3_pre_Pt3_on_tot.out</tt>
</TD></TR></TABLE>
  
In addition to these statistical outputs, the minimum number (N) or fraction (D) of unique clonotypes constituting 10, 25, 50, or 90 percent of the repertoire (e.g. D25) were calculated post-hoc from the tabulated repertoire files, when sorted in Microsoft Excel.

<P><h2>III.<tt>linepick_7.7.2.py</tt></h2>
<h3>Input:</h3>
Run from the top-level working directory, this program extracts the single-repertoire and paired-repertoire entropy statistics lines from all files *_H_CL_JS.out found inside directories of the same prefix as those files.

<h3>Output:</h3>
<tt>linepick_entropy.out</tt><BR>
<tt>linepick_divergence.out</tt><BR>
<P>All lines are extracted from the *_H_CL_JS.out file described in Col0.
<P>These files are suitable tables for all downstream and graphical analyses of the repertoires.  Values from Col10 (Num_CDR3) and Col4 (CLcdr3) were used to generate Figure 5B. 

<P><h2>IV. <tt>JS_aaCDR3perVJ_v2.5.py [filename]</tt></h2>
<h3>Input:</h3>
After running <tt>JS_splitter_Imsqv4.5.py</tt>, in the working directory, each input sample file is represented by <tt>[samplename]_only.productive.tsv</tt>.  These are the expected [filename] inputs.

<h3>Output:</h3>
Produces tabulated output file <tt>[samplename]_aaCDR3perVJ.tsv</tt>, with columns as follows:
<UL>
<LI>Col0:  VJcombo (name of V-J cassette combination)
<LI>Col1:  Number_aaCDR3 (number of aaCDR3 encoded by VJcombo in Col0)
<LI>Col2:  H_aaCDR3 (entropy of aaCDR3 encoded by VJcombo in Col0)
<LI>Col3:  Hnorm_aaCDR3 (entropy of aaCDR3 encoded by VJcombo in Col0, divided by maximum entropy; log2(N) where N is aaCDR3 population size)
</UL>

<h3>Example:</h3>
<TABLE><TR><TD>
<P><tt>python JS_aaCDR3perVJ_v2.5.py Pt3_pre_only.productive.tsv</tt>
<h3>.<BR>.<BR>.</h3>
<tt>Pt3_pre_only_aaCDR3perVJ.tsv</tt>
</TD></TR></TABLE>
<P>The summary statistics from these output files were used to generate Figures 5C and 5D in <A HREF="https://doi.org/10.1016/j.cell.2017.09.028">Riaz et al. (2017)</A>.  Figured 5C and 5D were generated using the median Hnorm of the aaCDR3 per VJ, as calculated post-hoc in Microsoft Excel from Col1 and Col3 of the *_aaCDR3perVJ.tsv output files as follows:  For each patient, the VJcombos and associated data from "pre" and "on" were merged, replacing the VJcombo identities in Col0 with "Pre" or "On" per the originating repertoire for that line.  This file was used as input for <A HREF="https://github.com/mwaskom/seaborn/blob/master/examples/multiple_joint_kde.py">multiple_joint_kde.py</A>, with "Pre" and "On" as the subsets from which to generate the kernel density estimate plots in Figure 5D.  Note that the Seaborn data visualization library must be installed (https://seaborn.pydata.org/index.html). 

<h1>Figure Creation</h1>
<p><h2>createTCRfigures.r</h2>
This script will generate Figure 5B and 5C using data generated from scripts above. 
See file for additional details.
