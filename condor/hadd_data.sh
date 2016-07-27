#!/bin/bash

#hadd /eos/uscms/store/user/clint/Run2016/July20/ljmet_trees/ljmet_Data_ElEl.root /eos/uscms/store/user/clint/Run2016/July20/DoubleEG*2/*.root
hadd -f /eos/uscms/store/user/clint/Run2016/July20/ljmet_trees/ljmet_Data_ElMu.root /eos/uscms/store/user/clint/Run2016/July20/MuonEG*2/*.root
hadd -f /eos/uscms/store/user/clint/Run2016/July20/ljmet_trees/ljmet_Data_MuMu.root /eos/uscms/store/user/clint/Run2016/July20/DoubleMu*2/*.root