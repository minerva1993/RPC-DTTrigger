DT + RPC study for improving RPC timing resolution.

*Setup
```{.Bash}
cmsrel CMSSW_10_6_1_pre2
cd CMSSW_10_6_1_pre2/src
git clone git@github.com:minerva1993/RPC-DTTrigger.git
scram b -j5

RPC-DTTrigger/RPCRecHitDTProducer/
cmsRun RPCRecHitDTProducer_cfg.py #Create sample rootfile from the MiniADO in eos
cd ../DTRPCTiming/test
cmsRun rpctest_cfg.py #input file path may be different
```
