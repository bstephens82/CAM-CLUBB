#!/usr/bin/env python

amwg_diags=False  #cheyenne only for now
profiles=True
adf_diags=False

test_cases=["b1850.054.f.Riexp0.5"]

cntl_case="b1850.054.f_ztest2"
username="stepheba"
PBS_ACCOUNT="P93300642"
test_first_yr=1
test_nyrs=3
cntl_first_yr=1979
cntl_nyrs=3
m2o_test=True
m2o_cntl=False
m2m=False
queue="regular"
adf_yaml_template="config_amwg_default_plots.yaml"
pbs_template="runder.adf"
adf_compare_obs=False

adf_dir="/glade/work/stepheba/post/ADF/"
profile_dir="/glade/work/stepheba/post/DiagnosticsAndTools/"



# --------------------------------------------------------------------------------------------------
# shouldn't need to modify below

import numpy as np
import sys
import os
import glob

def replace_line(f,str_org,str_replace):
    with open(f, 'r') as file:
        data = file.readlines()
        for line_num,line in enumerate(data):
            if line.startswith(str_org):
                data[line_num]=str_replace+"\n"
    file.close()

    with open(f,'w') as file:
        file.writelines(data)
    file.close()

# if a string occurs more than once, use this to specify
# which instance "n" of the string to replace
def replace_nth_line(f,str_org,str_replace,n):
    instance = 0
    with open(f, 'r') as file:
        data = file.readlines()
        for line_num,line in enumerate(data):
            if line.startswith(str_org):
                instance = instance + 1
                if instance == n:
                    data[line_num]=str_replace+"\n"
    file.close()

    with open(f,'w') as file:
        file.writelines(data)
    file.close()



##############################################################################################
######################################### ADF DIAGS ##########################################
##############################################################################################
if adf_diags:
    
    print("creating ADF batch script...")
    
    test_last_yr = test_first_yr + test_nyrs - 1
    cntl_last_yr = cntl_first_yr + cntl_nyrs - 1

    os.chdir(adf_dir)
 
    for test in test_cases:
        yamlfile="config_"+test+".yaml"
        os.system("cp "+adf_yaml_template+" "+yamlfile)
        replace_line(yamlfile,"    compare_obs:","    compare_obs: "+str(adf_compare_obs))
        replace_nth_line(yamlfile,"    cam_case_name: ","    cam_case_name: "+test,1)
        replace_nth_line(yamlfile,"    start_year:","    start_year: "+str(test_first_yr),1)
        replace_nth_line(yamlfile,"    end_year:","    end_year: "+str(test_last_yr),1)
        if not adf_compare_obs:
            replace_nth_line(yamlfile,"    cam_case_name: ","    cam_case_name: "+cntl_case,2)
            replace_nth_line(yamlfile,"    start_year:","    start_year: "+str(cntl_first_yr),2)
            replace_nth_line(yamlfile,"    end_year:","    end_year: "+str(cntl_last_yr),2)

        pbsfile="runder.adf_"+test
        os.system("cp "+pbs_template+" "+pbsfile)
        replace_line(pbsfile,"./run_adf_diag","./run_adf_diag "+yamlfile)
        os.system("qsub "+pbsfile)


##############################################################################################
#################################### CLUBB PROFILE PLOTS #####################################
##############################################################################################
if profiles:

    print("submitting CLUBB profiles batch script...")

    os.chdir(profile_dir)

    if len(test_cases) == 1:
        profile_file = "CAM_CLUBB_diag_single.py"
        replace_line(profile_file,"case=","case="+'"'+test_cases[0]+'"')
        replace_line(profile_file,"years=","years=["+str(test_first_yr)+"]")
        replace_line(profile_file,"nyear=","nyear=["+str(test_nyrs)+"]")
        os.system("qsub runder.zhun_diags_single")
    else:
        print('multiple profile cases not supported yet.')



#############################################################################
###############################  AMWG DIAGS #################################
#############################################################################
#if $amwg_diags
#then
#
#echo "running AMWG diagnostics..."
#
##test environment
#hash create_postprocess 2>/dev/null || { echo >&2 "The cesm_pp_activate comand has not been issued.  Aborting."; exit 1; }
#
## ----------------------------
##   Model vs. Obs (test)
## ----------------------------
#if $m2o_test
#then
#
#last_yr=$(expr $first_yr + $nyrs - 1)
#
#for test_case in ${test_cases[@]}; do
#
#if [ ! -d "/glade/scratch/$username/archive/$test_case" ]; then
#echo "test_case doesn't exist"
#exit
#fi
#
#if [ -d "./$test_case" ]; then
#echo "removing previous directories..."
#rm -rf $test_case
#rm -rf ../diagnostics-output/atm/climo/$test_case
#rm -rf ../diagnostics-output/atm/diag/$test_case
#rm -rf ../diagnostics-output/atm/diag/$test_case-obs.${first_yr}_${nyrs}
#fi
#
#create_postprocess --caseroot /glade/work/stepheba/post/cesm-postprocess/$test_case
#
#cd /glade/work/stepheba/post/cesm-postprocess/$test_case
#
#./pp_config --set DOUT_S_ROOT=/glade/scratch/$username/archive/$test_case/
#./pp_config --set ATM_GRID=0.9x1.25
#./pp_config --set ATMDIAG_OUTPUT_ROOT_PATH=/glade/scratch/stepheba/diagnostics-output/atm
#./pp_config --set ATMDIAG_test_first_yr=$first_yr
#./pp_config --set ATMDIAG_test_nyrs=$nyrs
#
#vim -E -s atm_averages <<-EOF
#        :%s/None/$PBS_ACCOUNT/
#        :update
#        :quit
#EOF
#
#qcmd -q $queue -- ./atm_averages --wait
#
#vim -E -s atm_diagnostics <<-EOF
#        :%s/None/$PBS_ACCOUNT/
#        :update
#        :quit
#EOF
#
#qcmd -q $queue -- ./atm_diagnostics --wait
#
##zip up the output diagnostics file
#cd /glade/scratch/stepheba/diagnostics-output/atm/diag/
#
#zip -rq $test_case-obs.$first_yr\_$last_yr.zip $test_case-obs.$first_yr\_$last_yr
#
#cd /glade/work/stepheba/post/cesm-postprocess/
#
#done 
#
#else
#
#echo "m2o_test = False"
#
#fi
#
#
## ----------------------------
##   Model vs. Obs (cntl)
## ----------------------------
#if $m2o_cntl
#then
#
## make sure cntl_case exists
#if [ ! -d "/glade/scratch/$username/archive/$cntl_case" ]; then
#echo "cntl_case doesn't exist"
#exit
#fi
#
## remove previous directories
#if [ -d "./$cntl_case" ]; then
#echo "removing previous directories..."
#rm -rf $cntl_case
#rm -rf ../diagnostics-output/atm/climo/$cntl_case
#rm -rf ../diagnostics-output/atm/diag/$cntl_case
#rm -rf ../diagnostics-output/atm/diag/$cntl_case-obs.${first_yr}_${nyrs}
#fi
#
#create_postprocess --caseroot /glade/work/stepheba/post/cesm-postprocess/$cntl_case
#
#cd /glade/work/stepheba/post/cesm-postprocess/$cntl_case
#
#./pp_config --set DOUT_S_ROOT=/glade/scratch/$username/archive/$cntl_case/
#./pp_config --set ATM_GRID=0.9x1.25
#./pp_config --set ATMDIAG_OUTPUT_ROOT_PATH=/glade/scratch/stepheba/diagnostics-output/atm
#./pp_config --set ATMDIAG_test_first_yr=$first_yr
#./pp_config --set ATMDIAG_test_nyrs=$nyrs
#
#vim -E -s atm_averages <<-EOF
#        :%s/None/$PBS_ACCOUNT/
#        :update
#        :quit
#EOF
#
#qcmd -q $queue -- ./atm_averages --wait
#
#vim -E -s atm_diagnostics <<-EOF
#        :%s/None/$PBS_ACCOUNT/
#        :update
#        :quit
#EOF
#
#qcmd -q $queue -- ./atm_diagnostics --wait
#
##zip up the output diagnostics file
##cd /glade/scratch/stepheba/diagnostics-output/atm/diag/
#
##zip -rq $cntl_case-obs.$first_yr\_$last_yr.zip $cntl_case-obs.$first_yr\_$last_yr
#
##DIRECTORY=$1
#
##scp -r /glade/scratch/stepheba/diagnostics-output/atm/diag/$cntl_case-obs.$first_yr\_$last_yr stepheba@moffatt.ucar.edu:/project/diagnostics/external/${DIRECTORY}
#
#else
#
#echo "m2o_cntl = False"
#
#fi
#
#
## ---------------------------------
##   Model (test) vs. Model (cntl)
## ---------------------------------
#if $m2m
#then
#
#cd /glade/work/stepheba/post/cesm-postprocess/$test_case
#
#./pp_config --set ATMDIAG_cntl_casename=$cntl_case
#./pp_config --set ATMDIAG_MODEL_VS_MODEL=True
#./pp_config --set ATMDIAG_cntl_first_yr=$first_yr
#./pp_config --set ATMDIAG_cntl_nyrs=$nyrs
#
## go!
##last_yr=$(expr $first_yr + $nyrs - 1)
#
#qcmd -q $queue -- ./atm_diagnostics --wait
#
#else
#
#echo "m2m = False"
#
#fi; fi
#
