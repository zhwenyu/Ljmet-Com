#!/bin/bash

rm -iv /uscms_data/d3/rsyarif/Brown2018/TT_BB_SSDL_LJMet_2017data/ljmet.tar 

### Jan 10, 2019 -- More Triggers

python TT_SSDL_Submit_FakeRate_Data_2017dataset.py --prefixFile Prefix_Data_FakeRate_2017dataset.py --inputFile Samples_Data_FakeRate_2017dataset.txt
python TT_SSDL_Submit_Data.py --prefixFile Prefix_Data2017.txt --inputFile Samples_Data2017.txt
python TT_SSDL_Submit_MC.py --prefixFile Prefix_BackgroundMC.txt --inputFile Samples_BackgroundMC.txt
python TT_SSDL_Submit_MC.py --prefixFile Prefix_InclusiveTT.txt --inputFile Samples_InclusiveTT.txt

### Jan 8, 2019

# python TT_SSDL_Submit_FakeRate_Data_2017dataset.py --prefixFile Prefix_Data_FakeRate_2017dataset.py --inputFile Samples_Data_FakeRate_2017dataset.txt


### 2018

# python TT_SSDL_Submit_Data.py --prefixFile Prefix_Data2017.txt --inputFile Samples_Data2017.txt
# 
# python TT_SSDL_Submit_MC.py --prefixFile Prefix_BackgroundMC.txt --inputFile Samples_BackgroundMC.txt
# # python TT_SSDL_Submit_MC_JECUP.py --prefixFile Prefix_BackgroundMC.txt --inputFile Samples_BackgroundMC.txt
# # python TT_SSDL_Submit_MC_JECDOWN.py --prefixFile Prefix_BackgroundMC.txt --inputFile Samples_BackgroundMC.txt
# # python TT_SSDL_Submit_MC_JERUP.py --prefixFile Prefix_BackgroundMC.txt --inputFile Samples_BackgroundMC.txt
# # python TT_SSDL_Submit_MC_JERDOWN.py --prefixFile Prefix_BackgroundMC.txt --inputFile Samples_BackgroundMC.txt
# 
# 
# python TT_SSDL_Submit_MC.py --prefixFile Prefix_InclusiveTT.txt --inputFile Samples_InclusiveTT.txt
# # python TT_SSDL_Submit_MC_JECUP.py --prefixFile Prefix_InclusiveTT.txt --inputFile Samples_InclusiveTT.txt
# # python TT_SSDL_Submit_MC_JECDOWN.py --prefixFile Prefix_InclusiveTT.txt --inputFile Samples_InclusiveTT.txt
# # python TT_SSDL_Submit_MC_JERUP.py --prefixFile Prefix_InclusiveTT.txt --inputFile Samples_InclusiveTT.txt
# # python TT_SSDL_Submit_MC_JERDOWN.py --prefixFile Prefix_InclusiveTT.txt --inputFile Samples_InclusiveTT.txt
# 
# 
#TESTS
# python TT_SSDL_Submit_Data.py --prefixFile Prefix_TESTDATA.txt --inputFile Samples_TESTDATA.txt
# python TT_SSDL_Submit_MC.py --prefixFile Prefix_TEST.txt --inputFile Samples_TEST.txt
