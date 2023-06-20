# darkShowerMiniValidation

Source cms installation then do (from slc7):
```bash
cmsrel CMSSW_12_4_9
cd CMSSW_12_4_9/src
git clone git@github.com:prijb/darkShowerMiniValidation.git CheckGen
scramv1 b -j 8
```

Then set correct path names in CheckDS/python/ConfFile_cfg.py and do
```bash
cmsRun CheckDS/python/ConfFile_cfg.py
```

Output tree contains muon branches (like muon_X) and L1 trigger decision branches (l1_name, l1_result, l1_prescale). 
The L1 trigger decisions checked are set in CheckDS/python/ConfFile_cfg.py


NB: may be able to use different CMSSW install for slc7 

