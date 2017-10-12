# Assessment of on-therapy neoantigen depletion

To evaluate the possibility of selective depletion of putative neoantigens within each individual patient, the pre-therapy number of neoantigens per synonymous mutation was determined, and the expected number of neoantigens on-therapy was computed using the measured number of on-therapy synonymous mutations.

To test whether the observed number of neoantigens was different from the expected count in each on-therapy tumor sample, we empirically calculated the following baseline rates pre-therapy:
1)	The average number of missense mutations per silent mutation across tumors.
2)	The number of neoantigens per missense mutation. Of note, since computational tools that predict neoantigens are dependent on patient-specific HLA class I alleles and some HLA-I alleles can bind greater numbers of peptides than others, we calculated this rate specific to the patient, rather than an average across all patients.
Using these two baseline rates, we used the number of silent mutations in each on-therapy tumor sample to calculate the expected number of neoantigens. Hence, the expected number of neoantigens in an on-therapy tumor sample is equal to the number of silent mutations in an on-therapy tumor sample multiplied by the number of missense mutations per silent mutation pre-therapy, multiplied by the number of neoantigens per missense mutation in the patient?s tumor sample pre-therapy. The chi-squared test was used to determine whether the observed number was different from the expected value. 


Note: The table 'Mut_NeoAg_Expected_vs_Observed_BMS038..csv' contains the computed values of the different rates