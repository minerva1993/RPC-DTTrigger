DT + RPC study for improving RPC timing resolution. First, follow the recipes from DTNtuple repo: https://github.com/battibass/DTNtuples

*Setup
```{.Bash}
#First follow here: https://github.com/battibass/DTNtuples
git clone git@github.com:minerva1993/RPC-DTTrigger.git
git cms-addpkg DataFormats/RPCRecHit
cat RPC-DTTrigger/patch.patch | patch -p1
scram b -j5

* Run analyzer
cd RPC-DTTrigger/RPCRecHitDTProducer/
cmsRun RPCRecHitDTProducer_cfg.py #Create sample rootfile from the MiniAOD
cd ../DTRPCTiming/test
cmsRun rpctest_cfg.py #input file path may be different

* Run producer
cd RPC-DTTrigger/DTRPCTimingUpdate/python
cmsRun dtrpcTimeUpdate_cfg.py inputFiles_load="input_PU200.txt"
```
