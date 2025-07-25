#!/bin/bash
source /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/setup.sh
source /datafs/users/dingxf/neutrino-physics/water-phase-analysis/juno-starter/InstallArea/setup.sh

python /datafs/users/dingxf/neutrino-physics/water-phase-analysis/juno-starter/JUNOStarter/share/run_sr.py --input-list $1 --evtmax 1000 --user-output $2
chmod +x $2
