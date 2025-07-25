#!/bin/bash
source /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.2.2/setup.sh
source ../../InstallArea/setup.sh

input=$1
output=$2

input="root://junoeos01.ihep.ac.cn//eos/juno/groups/Reconstruction/dingxf/Realdata_rec/run_3872/rtraw2rec/soft_trigger/Run3872_rtraw_soft_trigger_001_J25.1.5.root"
# input_correlation="/junofs/production/commissioning/rtraw/2025/0203/RUN.3338.JUNODAQ.Calib-ACU-AmBe-Global-Pos-x0y0z0.ds-4.soft_trigger.20250203112813.001_J25.1.2.rtraw"
output="test.root"

python run_sr.py --evtmax -1 \
    --input $input --rec-file --user-output $output \
    --mutiplicity_cut \
    --muon_veto \
    --save_muon \
    --muon_track_veto \
    --save_track_veto \
    --muon_track_file /datafs/users/wujxy/waterphase_analysis/realdata_analysis/cosmo_water/gl_muon_track/merged_3872.root \
    --bkg_filter
    # --input-correlation $input_correlation \
    # --flasher_ide \
    # --coin_select \

