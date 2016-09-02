cd ../
cp -v src/singleLepEventSelector.cc.bk2.beforeTriggerUpdate src/singleLepEventSelector.cc
scram b -j8
cd condor_TTtrilep/
