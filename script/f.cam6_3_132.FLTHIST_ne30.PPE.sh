cd /glade/work/pel/cesm_tags/PPE/cime/scripts

./create_newcase  --compset  FLTHIST --res ne30pg3_ne30pg3_mg17 \
--case /glade/p/cesmdata/cseg/runs/cesm2_0/f.cam6_3_132.FLTHIST_ne30.PPE  \
--run-unsupported  --pecount 2160 --project P93300642

## add COSP and setup

cd /glade/p/cesmdata/cseg/runs/cesm2_0/f.cam6_3_132.FLTHIST_ne30.PPE  
./xmlchange --append CAM_CONFIG_OPTS=-cosp
./case.setup

## we want to do an hybrid case from a previous run
./xmlchange RUN_STARTDATE=1995-01-01
./xmlchange STOP_N=2
./xmlchange STOP_OPTION=nyears
./xmlchange RESUBMIT=5
./xmlchange RUN_TYPE=hybrid
./xmlchange RUN_REFCASE=f.cam6_3_107.FLTHIST_v0a.ne30.clm5_1.001
./xmlchange RUN_REFDATE=1994-01-01
./xmlchange GET_REFCASE=TRUE
./xmlchange RUN_REFDIR=cesm2_init

## cam namelist
cd /glade/p/cesmdata/cseg/runs/cesm2_0/f.cam6_3_132.FLTHIST_ne30.PPE

cat <<EOF> user_nl_cam
dust_emis_fact = 1.3
mfilt    =       0,       5,     20,      40,      12,       120,      1,   1
nhtfrq              =       0,     -24,    -24,      -3,       0,       -2,      0,  -8760
ndens               =       2,       2,      2,       2,       2,       1,      2,   1
interpolate_output  =  .true.,  .true., .true., .false., .false., .true.,  .true.
interpolate_nlat    =     192,     192,    192,     192,     192,     192,   192
interpolate_nlon    =     288,     288,    288,     288,     288,     288,   288

empty_htapes = .true.

fincl1 = 'ACTNI', 'ACTNL', 'ACTREI', 'ACTREL', 'AODDUST', 'AODVIS', 'AODVISdn','BURDENBC', 'BURDENDUST', 'BURDENPOM', 'BURDENSEASALT',
'BURDENSO4', 'BURDENSOA', 'CAPE', 'CCN3', 'CDNUMC', 'CH4', 'CLDHGH', 'CLDICE', 'CLDLIQ', 'CLDLOW', 'CLDMED', 'CLDTOT', 'CLOUD', 'CMFMC_DP',
'CT_H2O', 'DCQ', 'DQCORE', 'DTCOND', 'DTCORE', 'DTV', 'EVAPPREC', 'EVAPSNOW', 'FCTI', 'FCTL', 'FICE', 'FLDS', 'FLNS', 'FLNSC', 'FLNT', 'FLNTC', 'FLUT',
'FREQZM', 'FSDS', 'FSDSC', 'FSNS', 'FSNSC', 'FSNT', 'FSNTC', 'FSNTOA', 'ICEFRAC', 'LANDFRAC', 'LHFLX', 'LWCF', 'MPDICE', 'MPDLIQ', 'MPDQ', 'MPDT',
'OCNFRAC', 'OMEGA', 'OMEGA500', 'PBLH', 'PHIS', 'PINT', 'PMID', 'PRECC', 'PRECL', 'PRECSC', 'PRECSL', 'PRECT', 'PS', 'PSL', 'PTEQ', 'PTTEND', 'Q',
'QFLX', 'QRL', 'QRS', 'QTGW', 'RCMTEND_CLUBB', 'RELHUM', 'RVMTEND_CLUBB', 'SHFLX', 'SOLIN', 'SST', 'STEND_CLUBB', 'SWCF',
'T', 'TAUX', 'TAUY', 'TFIX', 'TGCLDIWP', 'TGCLDLWP', 'TMQ', 'TREFHT', 'TS', 'TTGW', 'U', 'U10', 'UBOT', 'UTGWORO', 'UTGW_TOTAL',
'V', 'VBOT', 'VTGWORO', 'VTGW_TOTAL', 'WPRTP_CLUBB', 'WPTHLP_CLUBB', 'Z3', 'ZMDQ', 'ZMDT', 'N2O', 'CO2','CFC11','CFC12',
'CLD_MISR','FISCCP1_COSP','CLD_CAL','CLD_MISR','CLDTOT_CAL','CLDHGH_CAL', 'CLDMED_CAL','CLDLOW_CAL','CLMODIS', 'AODVISdn', 'AODDUSTdn',
'CCN3', 'CDNUMC', 'H2O', 'NUMICE', 'NUMLIQ'

fincl3 = 'PRECT', 'PRECC', 'FLUT', 'U850', 'U200', 'V850', 'V200', 'OMEGA', 'PSL'

fincl4 =  'PRECC','PRECL'

fincl5 = 'Uzm','Vzm','Wzm','THzm', 'VTHzm','WTHzm','UVzm','UWzm'
phys_grid_ctem_nfreq=-6
phys_grid_ctem_zm_nbas=120
phys_grid_ctem_za_nlat=90

fincl7= 'AQSO4_H2O2','AQSO4_O3', 'bc_a1', 'bc_a4', 'dst_a1', 'dst_a2', 'dst_a3', 'ncl_a1',
'ncl_a1', 'ncl_a2', 'ncl_a3', 'pom_a1', 'pom_a4', 'so4_a1', 'so4_a2', 'so4_a3',
'soa_a1', 'num_a1', 'num_a2', 'num_a3', 'num_a4',
'bc_a1SFWET', 'bc_a4SFWET', 'dst_a1SFWET', 'dst_a2SFWET', 'dst_a3SFWET', 'ncl_a1SFWET',
'ncl_a2SFWET', 'ncl_a3SFWET', 'pom_a1SFWET', 'pom_a4SFWET', 'so4_a1SFWET', 'so4_a2SFWET', 'so4_a3SFWET', 'soa_a1SFWET',
'soa_a2SFWET', 'bc_c1SFWET', 'bc_c4SFWET', 'dst_c1SFWET', 'dst_c2SFWET', 'dst_c3SFWET', 'ncl_c1SFWET', 'ncl_c2SFWET',
'ncl_c3SFWET', 'pom_c1SFWET', 'pom_c4SFWET', 'so4_c1SFWET', 'so4_c2SFWET', 'so4_c3SFWET', 'soa_c1SFWET', 'soa_c2SFWET',
'bc_a1DDF', 'bc_a4DDF', 'dst_a1DDF', 'dst_a2DDF', 'dst_a3DDF', 'ncl_a1DDF', 'ncl_a2DDF', 'ncl_a3DDF',
'pom_a1DDF', 'pom_a4DDF', 'so4_a1DDF', 'so4_a2DDF', 'so4_a3DDF', 'soa_a1DDF', 'soa_a2DDF',
'so4_a1_CLXF', 'so4_a2_CLXF', 'SFbc_a4', 'SFpom_a4', 'SFso4_a1', 'SFso4_a2',
'so4_a1_sfgaex1', 'so4_a2_sfgaex1', 'so4_a3_sfgaex1', 'soa_a1_sfgaex1', 'soa_a2_sfgaex1',
'SFdst_a1','SFdst_a2', 'SFdst_a3', 'SFncl_a1', 'SFncl_a2', 'SFncl_a3',
'num_a2_sfnnuc1', 'SFSO2', 'OCN_FLUX_DMS', 'SAD_SULFC', 'SAD_TROP', 'SAD_AERO'

clubb_l_diag_Lscale_from_taus=.true.

EOF

## clm namelists
cat <<EOF> user_nl_clm
use_init_interp = .true.
check_finidat_year_consistency = .false.
fsurdat = '/glade/work/slevis/git/mksurfdata_toolchain/tools/mksurfdata_esmf/surfdata_ne30np4.pg3_SSP5-8.5_78pfts_CMIP6_1850-2100_c230227.nc'
flanduse_timeseries = '/glade/work/slevis/git/mksurfdata_toolchain/tools/mksurfdata_esmf/landuse.timeseries_ne30np4.pg3_SSP5-8.5_78_CMIP6_1850-2100_c230227.nc'
EOF

## copy land mods
cd /glade/p/cesmdata/cseg/runs/cesm2_0/f.cam6_3_132.FLTHIST_ne30.PPE
cp /glade/p/cesmdata/cseg/runs/cesm2_0/f.cam6_3_119.FLTHIST_ne30.r328_gamma0.33_soae.001/SourceMods/src.clm/* SourceMods/src.clm/
#
# Source mods from Ben
#
cp /glade/work/pel/cesm_tags/PPE/SourceMods/src.cam/*.F90 SourceMods/src.cam/

## make cross checking
#cd /glade/p/cesmdata/cseg/runs/cesm2_0/
#xxdiff f.cam6_3_119.FLTHIST_ne30.r328_gamma0.33_soae.001 f.cam6_3_132.FLTHIST_ne30.PPE &


#cd /glade/p/cesmdata/cseg/runs/cesm2_0/f.cam6_3_132.FLTHIST_ne30.PPE
#qcmd -A CESM0023 -- ./case.build

# test

##cd /glade/p/cesmdata/cseg/runs/cesm2_0/f.cam6_3_132.FLTHIST_ne30.PPE
##./xmlchange PROJECT=CESM0023,RESUBMIT=0,STOP_N=1,STOP_OPTION=nmonths
##./xmlchange REST_OPTION=nmonths,REST_N=1
##./case.submit

#cd /glade/p/cesmdata/cseg/runs/cesm2_0/f.cam6_3_132.FLTHIST_ne30.PPE
#./xmlchange PROJECT=CESM0023,RESUBMIT=5,STOP_N=2,STOP_OPTION=nyears
#./xmlchange REST_OPTION=nyears,REST_N=1
#./xmlchange JOB_WALLCLOCK_TIME=12:00:00
#./case.submit
