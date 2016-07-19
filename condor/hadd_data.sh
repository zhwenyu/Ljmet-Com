#!/bin/bash

#hadd /eos/uscms/store/user/lpctlbsm/clint/Run2016B/July12/ljmet_trees/ljmet_Data_ElEl.root /eos/uscms/store/user/lpctlbsm/clint/Run2016B/July12/DoubleEG*2/*.root
hadd -f /eos/uscms/store/user/lpctlbsm/clint/Run2016B/July12/ljmet_trees/ljmet_Data_ElMu.root /eos/uscms/store/user/lpctlbsm/clint/Run2016B/July12/MuonEG*2/*.root
hadd -f /eos/uscms/store/user/lpctlbsm/clint/Run2016B/July12/ljmet_trees/ljmet_Data_MuMu.root /eos/uscms/store/user/lpctlbsm/clint/Run2016B/July12/DoubleMu*2/*.root