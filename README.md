# Init Majorana repository
```
cmsrel CMSSW_8_0_20
cd CMSSW_8_0_20/src
cmsenv
git cms-init
git clone https://github.com/GhentAnalysis/Majorana
```

# How to produce tuples (tested on T2_BE_IIHE)

Go to Majorana/PatAnalyzer/test
Update productionLabel in submitAll.py
Update datasets.txt with needed crab or pnfs paths
Run ./submitAll.py datasets.txt

Output will be found at
/pnfs/iihe/cms/store/user/$USER/majorana/<datasetname>/crab_<productionLabel>/*/*/trilepton_X.root (for submissions with crab)
/user/$USER/public/majorana/<datasetname>/local_<productionLabel>/trilepton_X.root (for local submissions)
