DT + RPC study for improving RPC timing resolution.

*Setup
```{.Bash}
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_6_1_patch2
cd CMSSW_10_6_1_patch2/src
cmsenv
git cms-init
git clone git@github.com:minerva1993/RPC-DTTrigger.git
git cms-addpkg DataFormats/RPCRecHit
git cms-addpkg SimMuon/RPCDigitizer
git remote add devel git@github.com:minerva1993/cmssw.git
git remote update
git checkout -t devel/dtrpc_devel
scram b -j5

cd RPC-DTTrigger/RPCRecHitDTProducer/
cmsRun RPCRecHitDTProducer_cfg.py #Create sample rootfile from the MiniADO in eos
cd ../DTRPCTiming/test
cmsRun rpctest_cfg.py #input file path may be different
```
