# use-HHT-for-1A0535-262
This is a record of my work on QPO phase resolved spectral analysis of 1A 0535+262 by HHT method.
This work and scripts is based on HXMT data analysis. You can also fallow the same logic for other data.

For the first step, you need to use the HXMT pipeline to do the general data analysis.

You can find the files in Workflow_Tutorial. 
1. Use vmd method to do HHT (see vmd_phase_pulsar_new.ipynb), and get the GTI files for each QPO phases you interest. 

2. Use phase_resolved_HXMT_gti.py to do QPO phase resolve.

3. Copy bkg files of from phase average results. See .cp files (cp0.sh and cp.sh).

4. Combine data. (cp_spec.sh,makeL_phase.sh, make_Comb_phase.sh, grppha_combine.sh, combine1_phase.sh, combine2_phase.sh, combine3_phase.sh) 
