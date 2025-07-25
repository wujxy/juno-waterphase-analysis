# ****************************************************************************/
# Author: Xuefeng Ding <dingxf@ihep.ac.cn> @ IHEP-CAS
#
# Date: 2025 January 27th
# Version: v1.0
# Description: Processor with SniperMT
#   Based on tut_calib2rec.py and
#     $JUNOTOP/Utilities/MtUtilities/share/mt_calib.py
#
# All rights reserved. 2023 copyrighted.
# ****************************************************************************/
# import DBISvc
# import CondDB
# import ParaDB
import argparse
import os.path
import sys

import BufferMemMgr  # noqa: F401
import Geometry  # noqa: F401
import RootIOSvc  # noqa: F401
import RootWriter  # noqa: F401
import Sniper
import SniperProfiling  # noqa: F401

import JUNOStarter  # noqa: F401


def register_options_omilrec(parser):
    parser.add_argument(
        "--loglevel",
        default="Info",
        choices=["Test", "Debug", "Info", "Warn", "Error", "Fatal"],
        help="Set the Log Level",
    )
    parser.add_argument("--np", type=int, default=1, help="number of threads")
    parser.add_argument("--evtmax", type=int, default=10, help="n-event")
    parser.add_argument(
        "--no-profiling",
        dest="profiling",
        action="store_false",
        help="disable profiling",
    )
    parser.set_defaults(profiling=True)
    parser.add_argument(
        "--no-profiling-details",
        dest="profiling_with_details",
        action="store_false",
        help="disable profiling with details saved",
    )
    parser.set_defaults(profiling_with_details=False)
    parser.add_argument(
        "--input",
        default="/junofs/users/dingxf/dingxf_data/b8_uniform_5_calib_EDM.root",
        help="input",
    )
    parser.add_argument("--input-list", default=None, help="input file name")
    parser.add_argument(
        "--user-output", default="sample_ANA_USER.root", help="user output"
    )
    parser.add_argument("--FOO", default="bar", help="set FOO")
    parser.add_argument("--daq-window-delay", default=0, type=int)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="mt_calib2rec", fromfile_prefix_chars="@"
    )
    register_options_omilrec(parser)
    args = parser.parse_args()

    ## load libraries
    # Sniper.loadDll("libMCParamsSvc.so")
    # Sniper.loadDll("libSpmtElecConfigSvc.so")
    # Sniper.loadDll("libPMTSimParamSvc.so")
    Sniper.loadDll("libPMTCalibSvc.so")
    Sniper.loadDll("libCondDB.so")
    Sniper.loadDll("libDBISvc.so")
    Sniper.loadDll("libDeconvolution.so")
    # Sniper.loadDll("libTimeRec.so")
    # Sniper.loadDll("libSPMTCalibAlg.so")

    task = Sniper.Task("TopTask/junotoptask")
    task.setEvtMax(args.evtmax)
    DATA_LOG_MAP = {"Test": 0, "Debug": 2, "Info": 3, "Warn": 4, "Error": 5, "Fatal": 6}
    task.setLogLevel(DATA_LOG_MAP[args.loglevel])

    if args.profiling:
        profiling = task.createSvc("SniperProfiling")
        if args.profiling_with_details:
            profiling.property("SaveDetails").set(True)

    import BufferMemMgr  # noqa

    bufMgr = task.createSvc("BufferMemMgr")
    bufMgr.property("TimeWindow").set([-0.1, 0.1])  # [-0.1 sec, 0.1 sec]

    import RootIOSvc  # noqa

    inputs = []
    if args.input_list:
        if not os.path.exists(args.input_list):
            sys.exit(-1)
        with open(args.input_list) as f:
            for line in f:
                line = line.strip()
                inputs.append(line)
    else:
        inputs.append(args.input)
    task.createSvc("RootInputSvc/InputSvc").property("InputFile").set(inputs)

    task.createSvc("RootWriter").property("Output").set(
        {"USER_OUTPUT": args.user_output}
    )
    task.createSvc("PMTParamSvc")

    alg = task.createAlg("JUNOStarter")
    alg.property("FOO").set(args.FOO)
    # alg.property("evtmax").set(args.evtmax)

    task.show()
    task.run()
