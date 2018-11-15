#!/bin/bash

# input=LJMet80x_TT_SSdilepton_2017_8_23_rizki
# input=LJMet80x_TT_SSdilepton_TpTpCalcOn_2017_9_1_rizki
input=LJMet80x_TT_SSdilepton_BkgMCGamma_2017_11_28_rizki

postStr=/ljmet_trees/
output=$input$postStr

echo "HADDING NOMINAL"
# python -u hadd_sigTT.py $input $output /
python -u hadd_MCbkg.py $input $output /

# echo "HADDING JECUP"
# python -u hadd_sigTT.py $input $output JECUP
# 
# echo "HADDING JECDOWN"
# python -u hadd_sigTT.py $input $output JECDOWN
# 
# echo "HADDING JERUP"
# python -u hadd_sigTT.py $input $output JERUP
# 
# echo "HADDING JERDOWN"
# python -u hadd_sigTT.py $input $output JERDOWN

echo "DONE"

