# C-Trap-Viscoelastic-Measurements
operating, passing commands to, and analyzing oscillation experiments with the Lumicks C-Trap

Oscillating C-trap.py: uses the bluelake python environment to pass commands to the C-trap mirror which controls bead movement

opening hfivefiles.m: after data has been exported from bluelake, takes the h5 files and parses data into .txt files for future analysis

Step1_get_TF_PSx_PSy_all_files.m: takes the .txt files containing trap position and PSD voltage to calculate TF, PSx, and PSy

OTC_nosqwv.m: takes the output TF and PS from the previous script and performs an optimization algorithm to determine the relevant viscoelastic parameters based on the model from Mizuno et al. 2008. Then uses this information to plot storage and loss modulus, among other things

TF, PS data for OTC: some sample data compiled using Step1_get_TF_PSx_PSy_all_files.m. This data can be called by OTC_nosqwv.m to familiarize oneself with doing Optical Trap Calibration. 

Additional Scripts: the main scripts earlier often call separate functions, particularly when calculating TF, PS, and performing least squares optimization. Such scripts can be found here, among others.

two trapped beads with oscillation.mp4: video showing what the experiment looks like under the brighfield objective.
