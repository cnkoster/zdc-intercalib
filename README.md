# zdc-intercalib

Code to run ZDC tower inter-calibration in A-A collisions

STEPS

[1] one needs to run the task zdc-task-intercalib over a sample of A02D data over the GRID

[2] the obtained output includes some control histograms (if the configurable variable writeHistos in the task is flagged (default is true) and a table containing 5 floats for each ZN, corresponsing to: PMC, PMQ1, PMQ2, PMQ3, PMQ4 light outputs (entries of the table are spcified in ZDCInterCalib.h class). The output table needs to be merged using aod-merger.

[3] the macro interCalibZN.C can then be used to read back the light outputs from the table and run the minimization
NB -> Better to run it compiled: root[0] .x interCalibZN.C+
NB -> One needs to properly set the root file name and the directory name inside the macro (the directory name can be seen opening the root file)
