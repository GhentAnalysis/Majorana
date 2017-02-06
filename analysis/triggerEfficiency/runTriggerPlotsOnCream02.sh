#!/bin/bash
cd "/user/tomc/Majorana/CMSSW_8_0_25/src/Majorana/analysis/triggerEfficiency"
source $VO_CMS_SW_DIR/cmsset_default.sh
eval `scram runtime -sh`
eval $command
