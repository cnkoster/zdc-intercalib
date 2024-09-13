export OPTIONS="-b --configuration json://config.json --aod-file AO2D544124.root --aod-writer-keep="AOD/ZDCIC/0""

o2-analysis-timestamp ${OPTIONS} | 
o2-analysis-event-selection ${OPTIONS} |  
o2-analysis-track-propagation ${OPTIONS} |  
o2-analysis-mm-zdc-task-intercalib ${OPTIONS} 