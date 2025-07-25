/*****************************************************************************/
// Author: Xuefeng DING <dingxf@ihep.ac.cn> @ IHEP-CAS
//
// Project: JUNOStarterESD
// Date: 2025 January 27th
// Version: v1.0
// Description:
//   Data processing algorithm
//
// Maintainer:
//   Xuefeng Ding <dingxf@ihep.ac.cn>
//
// All rights reserved. 2024 copyrighted.
/*****************************************************************************/
#include "JUNOStarterESD.h"
#include "PMTCalibSvc/IPMTCalibSvc.hh"
#include <Event/CdLpmtCalibEvt.h>
#include <Event/CdLpmtCalibHeader.h>
#include <Event/CdLpmtElecEvt.h>
#include <Event/CdLpmtElecHeader.h>
#include <Event/CdTriggerEvt.h>
#include <Event/CdTriggerHeader.h>
#include <Event/CdVertexRecEvt.h>
#include <Event/CdVertexRecHeader.h>
#include <Event/WpCalibEvt.h>
#include <Event/WpCalibHeader.h>
#include <EvtNavigator/EvtNavHelper.h>
#include <EvtNavigator/EvtNavigator.h>
#include <EvtNavigator/NavBuffer.h>
#include <Geometry/IPMTParamSvc.h>
#include <Identifier/CdID.h>
#include <Identifier/IDService.h>
#include <Identifier/Identifier.h>
#include <RootWriter/RootWriter.h>
#include <SniperKernel/AlgFactory.h>
#include <SniperKernel/SniperLog.h>
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <cstdlib>
#include <iterator>
#include <limits>
#include <sys/types.h>
#include <tuple>

int64_t safe_subtract(uint64_t a, uint64_t b) {
  constexpr uint64_t mask = (1ULL << 48) - 1;
  a = a & mask + (1ULL << 48);
  b = b & mask;
  auto diff = (a - b) & mask;
  auto flag = (diff >> 47) << 48;
  return diff - flag;
}

bool compare_low48(uint64_t t_ref, uint64_t key) {
  return safe_subtract(t_ref, key) < 0;
}

static std::unordered_map<uint32_t, uint32_t> m_pmt_num_map;
static std::map<uint64_t,
                std::tuple<double, double, double, double, double, double>>
    gl_muon_track_map;
static std::map<
    uint64_t,
    std::tuple<std::tuple<double, double, double, double, double, double>,
               uint64_t>>
    sw_muon_track_map;
static std::vector<std::tuple<int, int, uint64_t>> track_veto_rules;

DECLARE_ALGORITHM(JUNOStarterESD);

JUNOStarterESD::JUNOStarterESD(const std::string &name) : AlgBase(name) {
  i_totPMTs = 0;

  declProp("bias", m_bias);
  declProp("rec_file", m_rec_file);
  declProp("correlate_rtraw", m_correlate_rtraw);
  declProp("muon_veto", m_muon_veto);
  declProp("muon_veto_pre", m_muon_veto_pre);
  declProp("muon_veto_post", m_muon_veto_post);
  declProp("mutiplicity_cut", m_mutiplicity_cut);
  declProp("flasher_ide", m_flasher_ide);
  declProp("coin_select", m_coin_select);
  declProp("bkg_filter", m_bkg_filter);
  declProp("muon_track_veto", m_muon_track_veto);
  declProp("muon_track_file", m_muon_track_file);
  declProp("save_track_veto", m_save_track_veto);
  declProp("save_muon", m_save_muon);
}

bool JUNOStarterESD::initialize() {
  m_cur_entry = -1;
  getParent()->show();
  SniperPtr<RootWriter> svc(getParent(), "RootWriter");
  if (svc.invalid()) {
    LogError << "Can't locate RootWriter" << std::endl;
    return false;
  }

  if (!initPMTPars())
    return false;

  if (m_rec_file) {
    m_edm_path = "/Event/CdLpmtCalib_FPGA";
    ;
  } else {
    m_edm_path = "/Event/CdLpmtCalib";
  }

  m_pmtTable.reserve(i_totPMTs);

  m_tree = svc->bookTree(*getParent(), "USER_OUTPUT/selectevt",
                         "evt after muon veto and multiplicity cut");
  m_tree->Branch("nHits", &m_nhits, "nHits/i");
  m_tree->Branch("nPMTs", &m_npmts, "nPMTs/i");
  m_tree->Branch("nPE", &m_npe, "nPE/F");
  m_tree->Branch("segment", &m_segment, "segment/i");

  m_segment_tree =
      svc->bookTree(*getParent(), "USER_OUTPUT/segment", "muon segment");
  m_segment_tree->Branch("t0", &m_segment_t0, "t0/l");
  m_segment_tree->Branch("t1", &m_segment_t1, "t1/l");
  m_segment_tree->Branch("live_t", &m_live_t, "live_t/l");

  m_pmt_tree =
      svc->bookTree(*getParent(), "USER_OUTPUT/pmt_perform", "pmt performance");
  m_pmt_tree->Branch("pmtid", &m_pmtid, "pmtid/i");
  m_pmt_tree->Branch("firednum", &m_firednum, "firednum/i");
  m_veto_tree = svc->bookTree(*getParent(), "USER_OUTPUT/veto", "muon info");
  m_veto_tree->Branch("n_skip_cd_muon", &m_n_skip_cd_muon, "n_skip_cd_muon/i");
  m_veto_tree->Branch("n_skip_wp_muon", &m_n_skip_wp_muon, "n_skip_wp_muon/i");
  m_veto_tree->Branch("n_skip_empty", &m_n_skip_empty, "n_skip_empty/i");
  m_veto_tree->Branch("n_skip_multi_r", &m_n_skip_multi_r, "n_skip_multi_r/i");
  m_veto_tree->Branch("n_delayed", &m_n_delayed, "n_delayed/i");
  m_veto_tree->Branch("n_skip_muon_track", &m_n_skip_muon_track,
                      "n_skip_muon_track/i");
  m_veto_tree->Branch("n_skip_track_veto", &m_n_skip_track_veto,
                      "n_skip_track_veto/i");
  m_veto_tree->Branch("n_finded_muon_track", &m_n_finded_muon_track,
                      "n_finded_muon_track/i");

  m_muon_tree = svc->bookTree(*getParent(), "USER_OUTPUT/muon", "muon info");
  m_muon_tree->Branch("npe", &m_muon_npe, "npe/F");
  m_muon_tree->Branch("track_tag", &m_muon_track_tag, "track_tag/i");
  m_muon_tree->Branch("npmt", &m_muon_npmt, "npmt/i");
  m_muon_tree->Branch("dt_last_muon", &m_dt_last_muon, "dt_last_muon/l");

  m_selectrec_tree = svc->bookTree(*getParent(), "USER_OUTPUT/selectrec",
                                   "rec after muon veto and multiplicity cut");
  m_selectrec_tree->Branch("x", &m_x, "recx/F");
  m_selectrec_tree->Branch("y", &m_y, "recy/F");
  m_selectrec_tree->Branch("z", &m_z, "recz/F");
  m_selectrec_tree->Branch("px", &m_px, "recpx/F");
  m_selectrec_tree->Branch("py", &m_py, "recpy/F");
  m_selectrec_tree->Branch("pz", &m_pz, "recpz/F");
  m_selectrec_tree->Branch("chisq_alt", &m_chisq_alt, "recchisq_alt/F");
  m_selectrec_tree->Branch("raw_entry", &m_raw_entry, "raw_entry/i");
  // m_selectrec_tree->Branch("Second", &m_evt_sec, "Second/i");
  // m_selectrec_tree->Branch("NanoSec", &m_evt_nsec, "NanoSec/i");
  m_selectrec_tree->Branch("cos_theta", &m_cos_theta, "cos_theta/F");
  m_selectrec_tree->Branch("E_rec_0_d5", &m_E_rec_0_d5, "E_rec/i");
  m_selectrec_tree->Branch("E_rec_m1_md8", &m_E_rec_m1_md8, "E_rec/i");
  m_selectrec_tree->Branch("E_rec_d8_1", &m_E_rec_d8_1, "E_rec/i");
  m_selectrec_tree->Branch("max_q_frac", &m_max_q_frac, "max_q_frac/F");
  m_selectrec_tree->Branch("core_q_frac", &m_core_q_frac, "core_q_frac/F");
  m_selectrec_tree->Branch("q_ratio_1m", &m_q_ratio_1m, "q_ratio_1m/F");
  m_selectrec_tree->Branch("q_ratio_2m", &m_q_ratio_2m, "q_ratio_2m/F");
  m_selectrec_tree->Branch("q_ratio_3m", &m_q_ratio_3m, "q_ratio_3m/F");
  m_selectrec_tree->Branch("oppo_q_frac", &m_oppo_q_frac, "oppo_q_frac/F");
  m_selectrec_tree->Branch("trigger_rate", &m_trigger_rate, "trigger_rate/F");
  m_selectrec_tree->Branch("dt_last_muon_track", &m_dt_last_muon_track,
                           "dt_last_muon_track/i");
  m_selectrec_tree->Branch("dr_last_muon_track", &m_dr_last_muon_track,
                           "dr_last_muon_track/F");

  m_coincidence_tree = svc->bookTree(*getParent(), "USER_OUTPUT/coincidence",
                                     "coincidence events");
  m_coincidence_tree->Branch("x", &m_x, "recx/F");
  m_coincidence_tree->Branch("y", &m_y, "recy/F");
  m_coincidence_tree->Branch("z", &m_z, "recz/F");
  m_coincidence_tree->Branch("px", &m_px, "recpx/F");
  m_coincidence_tree->Branch("py", &m_py, "recpy/F");
  m_coincidence_tree->Branch("pz", &m_pz, "recpz/F");
  m_coincidence_tree->Branch("chisq_alt", &m_chisq_alt, "recchisq_alt/F");
  m_coincidence_tree->Branch("x_coin", &m_x_coin, "x_coin/F");
  m_coincidence_tree->Branch("y_coin", &m_y_coin, "y_coin/F");
  m_coincidence_tree->Branch("z_coin", &m_z_coin, "z_coin/F");
  m_coincidence_tree->Branch("px_coin", &m_px_coin, "px_coin/F");
  m_coincidence_tree->Branch("py_coin", &m_py_coin, "py_coin/F");
  m_coincidence_tree->Branch("pz_coin", &m_pz_coin, "pz_coin/F");
  m_coincidence_tree->Branch("chisq_alt_coin", &m_chisq_alt_coin,
                             "chisq_alt_coin/F");
  m_coincidence_tree->Branch("dt_coin", &m_dt_coin, "dt_coin/i");

  m_cluster_tree =
      svc->bookTree(*getParent(), "USER_OUTPUT/cluster", "cluster events");
  m_cluster_tree->Branch("dt_pre_muon", &m_dt_pre_muon, "dt_pre_muon/i");

  m_muon_track_tree = svc->bookTree(*getParent(), "USER_OUTPUT/trackvetoevt",
                                    "event info under track veto");
  m_muon_track_tree->Branch("dt", &m_dt, "dt/i");
  m_muon_track_tree->Branch("d_track", &m_d_track, "d_track/F");
  m_muon_track_tree->Branch("x", &m_x, "recx/F");
  m_muon_track_tree->Branch("y", &m_y, "recy/F");
  m_muon_track_tree->Branch("z", &m_z, "recz/F");
  m_muon_track_tree->Branch("pz", &m_pz, "recpz/F");
  m_muon_track_tree->Branch("chisq_alt", &m_chisq_alt, "recchisq_alt/F");

  m_t_ref = 0;
  m_t_prev = 0;
  m_t_muon = 0;
  m_t_muon_track = 0;
  m_muon_npe = 0;
  m_muon_npmt = 0;
  m_muon_track_npmt = 0;
  m_muon_track_npe = 0;
  m_segment = std::numeric_limits<uint32_t>::max();
  m_t_muon_prev = std::numeric_limits<uint64_t>::max();
  m_t_muon_track_prev = std::numeric_limits<uint64_t>::max();
  m_dt_last_muon_track = 0;
  m_dr_last_muon_track = 0;
  m_trigger_rate = 0;
  m_n_delayed = 0;
  m_n_skip_cd_muon = 0;
  m_n_skip_wp_muon = 0;
  m_n_skip_empty = 0;
  m_n_skip_multi_r = 0;
  m_n_skip_muon_track = 0;
  m_n_skip_track_veto = 0;
  m_n_finded_muon_track = 0;

  if (m_muon_track_veto) {
    m_muon_veto_post = 50'000ULL; // 50us
  }

  // load muon track map
  if (m_muon_track_veto) {
    LoadMuonTrack();
    track_veto_rules = {{0, 1000, 5'000'000'000ULL},
                        {1000, 3000, 4'000'000'000ULL},
                        {3000, 4000, 2'000'000'000ULL},
                        {4000, 5000, 200'000'000ULL}};
  }

  // get buffer
  SniperDataPtr<JM::NavBuffer> navBuf(getParent(), "/Event");
  if (navBuf.invalid()) {
    LogError << "cannot get the NavBuffer @ /Event" << std::endl;
    return false;
  }
  m_navBuf = navBuf.data();

  LogDebug << "Initialized" << std::endl;
  return true;
}

bool JUNOStarterESD::execute() {
  m_t_pre = m_t_ref;

  LogDebug << "executing" << std::endl;
  LogDebug << "Processing entry " << m_cur_entry << std::endl;

  ++m_cur_entry;
  auto *nav = m_navBuf->curEvt();
  auto ts_ref = nav->TimeStamp();
  uint64_t t_ref = ts_ref.GetSec() * 1'000'000'000ULL + ts_ref.GetNanoSec();
  m_t_ref = t_ref;
  m_t_ref_correct = t_ref;

  // muon track veto
  if (m_muon_track_veto) {
    if (muon_track_veto(nav)) {
      return true;
    }
  }

  // muon veto
  if (m_muon_veto) {
    if (muon_veto(nav)) {
      return true;
    }
  }

  // only for CD events
  if (nav->getDetectorType() != 0) {
    return true;
  }
  auto *calibhdr = JM::getHeaderObject<JM::CdLpmtCalibHeader>(nav, m_edm_path);
  if (calibhdr == nullptr) {
    ++m_n_skip_empty;
    return true;
  }
  auto *calibevt = dynamic_cast<JM::CdLpmtCalibEvt *>(calibhdr->event());
  if (!calibevt) {
    LogFatal << "No event found in the header " << m_cur_entry << std::endl;
    return false;
  }

  // multiplicity cut
  if (m_mutiplicity_cut) {
    if (mutiplicity_cut(nav)) {
      return true;
    }
  }

  // skip periodic trigger type event
  if (m_correlate_rtraw) {
    auto cdtriggerhdr = JM::getHeaderObject<JM::CdTriggerHeader>(nav);
    if (cdtriggerhdr != nullptr) {
      auto cdtriggerevt =
          dynamic_cast<JM::CdTriggerEvt *>(cdtriggerhdr->event());
      for (auto trigger_type : cdtriggerevt->triggerType()) {
        if (trigger_type == "Periodic") {
          m_n_skip_periodic_tri++;
          return true;
        }
      }
    }
  }

  // save the event after muon veto and multiplicity cut
  if (m_rec_file) {
    get_recevt(nav);
  }

  // filter bkg events
  if (m_bkg_filter && check_bkg()) {
    return true;
  }

  // calculate the dr and dt between the last muon track and current event
  if (m_muon_track_veto && m_save_track_veto) {
    m_dt_last_muon_track = std::numeric_limits<int64_t>::max();
    m_dr_last_muon_track = std::numeric_limits<float>::max();
    for (auto it = m_navBuf->current(); it != m_navBuf->begin(); --it) {
      auto *next = it->get();
      auto ts_next = next->TimeStamp();
      uint64_t t_next =
          ts_next.GetSec() * 1'000'000'000ULL + ts_next.GetNanoSec();
      auto t_diff = t_ref - t_next;
      if (t_diff > 5'000'000'000ULL) {
        // 200ms
        break;
      }
      if (sw_muon_track_map.find(t_next) != sw_muon_track_map.end()) {
        auto muon_track = sw_muon_track_map[t_next];
        auto track = std::get<0>(muon_track);
        auto t_track = std::get<1>(muon_track);
        auto ev_xyz = std::array<double, 3>{m_x, m_y, m_z};
        auto dt = safe_subtract(t_ref, t_track);
        auto dr = m_jcs.get_track_distance(ev_xyz, track);
        // if (dr < m_dr_last_muon_track) {
        if (dr < 2'000) {
          m_dr_last_muon_track = dr;
          m_dt_last_muon_track = dt;
          break;
        }
      }
    }
  }

  // check cluster events and get dt_pre_muon
  if (check_cluster()) {
    m_dt_pre_muon = t_ref - m_t_muon;
    m_cluster_tree->Fill();
  }

  // loop over all fired PMTs
  m_npmts = calibevt->calibPMTCol().size();
  m_npe = npe(nav, m_edm_path);
  m_nhits = 0;
  m_E_rec_0_d5 = 0;
  m_E_rec_m1_md8 = 0;
  m_E_rec_d8_1 = 0;
  m_max_q_frac = 0;
  float max_charge = 0;
  float total_charge = 0;
  unsigned int max_q_copyno = 0;
  for (auto *channel : calibevt->calibPMTCol()) {
    m_nhits += channel->size();

    // count pmt fired number
    if (m_npmts < 1000) {
      auto pmtid = channel->pmtId();
      m_pmt_num_map[pmtid]++;
      // LogDebug << "pmtid: " << pmtid << " firednum: " << m_pmt_num_map[pmtid]
      //  << std::endl;
    }

    // process channel data
    if (m_rec_file) {
      unsigned int pmtId = channel->pmtId();
      Identifier id = Identifier(pmtId);
      unsigned int thisid = CdID::module(id);
      auto *pmtprop = m_pmtTable.at(thisid);
      std::array<double, 3> pmt = {pmtprop->pos.X(), pmtprop->pos.Y(),
                                   pmtprop->pos.Z()};
      std::array<double, 3> ev_xyz = {m_x, m_y, m_z};
      std::array<double, 3> ev_pxyz = {m_px, m_py, m_pz};
      // reconstruct energy based on cherenkov angle
      // auto fht = channel->firstHitTime();
      // auto t_corr = fht - m_jcs.get_tof(ev_xyz, pmt);
      auto hit_ev_direct_cosTheta =
          m_jcs.get_hit_ev_direct_cosTheta(ev_xyz, pmt, ev_pxyz);
      if (hit_ev_direct_cosTheta < 1.0 && hit_ev_direct_cosTheta > 0.25) {
        m_E_rec_0_d5++;
      }
      if (hit_ev_direct_cosTheta < -0.6 && hit_ev_direct_cosTheta > -0.8) {
        m_E_rec_m1_md8++;
      }
      if (hit_ev_direct_cosTheta < 1.0 && hit_ev_direct_cosTheta > 0.8) {
        m_E_rec_d8_1++;
      }
      auto charge = channel->sumCharge();
      total_charge += charge;
      if (charge > max_charge) {
        max_charge = charge;
        max_q_copyno = thisid;
      }
    }
  }

  // caculate q fraction in different zones (after get max q pmt)
  m_core_q_frac = 0;
  m_q_ratio_1m = 0;
  m_q_ratio_2m = 0;
  m_q_ratio_3m = 0;
  m_oppo_q_frac = 0;
  if (m_rec_file && m_flasher_ide) {
    m_max_q_frac = max_charge / total_charge;
    // loop again
    auto *max_q_pmtprop = m_pmtTable.at(max_q_copyno);
    std::array<double, 3> max_q_pmt = {
        max_q_pmtprop->pos.X(), max_q_pmtprop->pos.Y(), max_q_pmtprop->pos.Z()};
    for (auto *channel : calibevt->calibPMTCol()) {
      auto *pmtprop = m_pmtTable.at(CdID::module(Identifier(channel->pmtId())));
      std::array<double, 3> pmt = {pmtprop->pos.X(), pmtprop->pos.Y(),
                                   pmtprop->pos.Z()};
      auto dist = m_jcs.get_distance(max_q_pmt, pmt);
      auto charge = channel->sumCharge();
      // core region
      if (dist < 3000) {
        m_core_q_frac += charge / total_charge;
      }
      // radius region
      if (dist < 1000) {
        m_q_ratio_1m += charge / max_charge;
      } else if (dist < 2000) {
        m_q_ratio_2m += charge / max_charge;
      } else if (dist < 3000) {
        m_q_ratio_3m += charge / max_charge;
      }
      // opposite region
      if (dist > 19365 * sqrt(2)) {
        m_oppo_q_frac += charge / total_charge;
      }
    }
    // get trigger rate
    get_trigger_rate(nav);
  }

  // find coincidence events
  if (m_coin_select) {
    find_coincidence(nav);
  }

  m_tree->Fill();

  if (m_rec_file) {
    m_selectrec_tree->Fill();
  }

  return true;
}

bool JUNOStarterESD::finalize() {
  // fill muon tree
  m_veto_tree->Fill();
  // fill segment tree
  m_segment_tree->Fill();
  // fill pmt tree
  for (auto &pmt_num : m_pmt_num_map) {
    m_pmtid = pmt_num.first;
    m_firednum = pmt_num.second;
    m_pmt_tree->Fill();
  }
  LogDebug << "finalizing" << std::endl;

  return true;
}

bool JUNOStarterESD::initPMTPars() {
  // std::lock_guard<std::mutex> lock(QCtrRecAlg_initPMTPars_mutex);

  // get the PMT positions and parameters
  SniperPtr<IPMTParamSvc> pmtsvc(getParent(), "PMTParamSvc");
  if (pmtsvc.invalid()) {
    LogError << "Failed to get service PMTParamSvc" << std::endl;
    return true;
  }
  IPMTParamSvc *m_pmtsvc = pmtsvc.data();

  // SniperPtr<IPMTCalibSvc> calSvc(*getParent(), "PMTCalibSvc");
  // if (calSvc.invalid()) {
  //   LogError << "Can't Locate  PMTCalibSvc." << std::endl;
  //   return false;
  // }

  // Read PMT parameters from CalibSVC
  LogDebug << "Loading PMTPara from SVC(database) ..." << std::endl;
  // LogDebug << "LGain[0]:" << calSvc->getGain().at(0) << std::endl;
  // LogDebug << "LDNR[0]:" << calSvc->getDarkRate().at(0) << std::endl;
  // LogDebug << "LPDE[0]:" << calSvc->getRelativeDE().at(0) << std::endl;

  i_totPMTs = m_pmtsvc->get_NTotal_CD_LPMT();
  // LogInfo << "DarkRate size=" << calSvc->getDarkRate().size()
  //         << ",rpde=" << calSvc->getRelativeDE().size()
  //         << ",getChargeSpec=" << calSvc->getChargeSpec().size() <<
  //         std::endl;
  for (unsigned int ith = 0; ith < i_totPMTs; ith++) {
    PMTProp *thisPMT = new PMTProp(ith);
    PMTTYPE thistype = PMTTYPE::NONE;
    if (m_pmtsvc->isHamamatsu(ith))
      thistype = PMTTYPE::LPMT_DYNODE;
    else if (m_pmtsvc->isNormalNNVT(ith))
      thistype = PMTTYPE::LPMT_MCP;
    else if (m_pmtsvc->isHighQENNVT(ith))
      thistype = PMTTYPE::LPMT_MCP_HIGHQE;
    else if (m_pmtsvc->isHZC(ith))
      thistype = PMTTYPE::SPMT;
    thisPMT->type = thistype;
    thisPMT->pos = 19.18 / 19.434 *
                   TVector3(m_pmtsvc->getPMTX(ith), m_pmtsvc->getPMTY(ith),
                            m_pmtsvc->getPMTZ(ith));
    // thisPMT->dnrate = calSvc->getDarkRate().at(ith);
    // thisPMT->rpde = calSvc->getRelativeDE().at(ith);
    // thisPMT->gain = calSvc->getChargeSpec().at(ith)->GetMean();
    /*
    htemp->Fit("gaus", "RQ0", "", 0.5, 1.5);
    TF1* fitfunc = htemp->GetFunction("gaus");
    thisPMT->gain = fitfunc->GetParameter(1);
    if(!m_pmtsvc->isHamamatsu(ith)) thisPMT->gain = htemp->GetMean();;
    delete fitfunc;
    */
    m_pmtTable.push_back(thisPMT);
  }
  return true;
}

int JUNOStarterESD::npmt(JM::EvtNavigator *nav, std::string edm_path) {
  if (nav->getDetectorType() == 1) {
    auto *wphdr = JM::getHeaderObject<JM::WpCalibHeader>(nav, edm_path);
    if (wphdr != nullptr) {
      auto *wp = dynamic_cast<JM::WpCalibEvt *>(wphdr->event());

      int npmt = 0;
      npmt = wp->calibPMTCol().size();
      return npmt;
    }
  }
  if (nav->getDetectorType() == 0) {
    auto *cdhdr = JM::getHeaderObject<JM::CdLpmtCalibHeader>(nav, edm_path);
    if (cdhdr != nullptr) {
      auto *cd = dynamic_cast<JM::CdLpmtCalibEvt *>(cdhdr->event());

      int npmt = 0;
      npmt = cd->calibPMTCol().size();
      return npmt;
    }
  }
  return 0;
}

float JUNOStarterESD::npe(JM::EvtNavigator *nav, std::string edm_path) {
  auto evt_time = nav->TimeStamp().GetSec() * 1'000'000'000ULL +
                  nav->TimeStamp().GetNanoSec();

  static std::unordered_map<uint64_t, float> m_cd_npe_map;
  static std::unordered_map<uint64_t, float> m_wp_npe_map;

  if (nav->getDetectorType() == 1) {
    if (m_wp_npe_map.find(evt_time) != m_wp_npe_map.end()) {
      return m_wp_npe_map[evt_time];
    }
    auto *wphdr = JM::getHeaderObject<JM::WpCalibHeader>(nav, edm_path);
    if (wphdr != nullptr) {
      auto *wp = dynamic_cast<JM::WpCalibEvt *>(wphdr->event());

      float npe = 0;
      for (auto *channel : wp->calibPMTCol()) {
        for (unsigned int hit_i = 0; hit_i < channel->size(); hit_i++) {
          auto q = channel->charge(hit_i);
          if (q > 100) {
            continue;
          } else {
            npe += q;
          }
        }
      }
      m_cd_npe_map[evt_time] = npe;
      return npe;
    }
  }
  if (nav->getDetectorType() == 0) {
    if (m_cd_npe_map.find(evt_time) != m_cd_npe_map.end()) {
      return m_cd_npe_map[evt_time];
    }

    auto *calibhdr = JM::getHeaderObject<JM::CdLpmtCalibHeader>(nav, edm_path);
    if (calibhdr != nullptr) {
      auto *calib = dynamic_cast<JM::CdLpmtCalibEvt *>(calibhdr->event());

      float npe = 0;
      for (auto *channel : calib->calibPMTCol()) {
        // npe += channel->sumCharge();
        for (unsigned int hit_i = 0; hit_i < channel->size(); hit_i++) {
          auto q = channel->charge(hit_i);
          if (q > 100) {
            continue;
          } else {
            npe += q;
          }
        }
      }
      m_wp_npe_map[evt_time] = npe;
      return npe;
    }
  }
  return 0;
}

bool JUNOStarterESD::check_muon(JM::EvtNavigator *nav) {
  auto npe = JUNOStarterESD::npe(nav, m_edm_path);
  if (nav->getDetectorType() == 1) {
    if (npe > 1000) {
      return true;
    }
  } else if (nav->getDetectorType() == 0) {
    if (npe > 10000) {
      return true;
    }
  }
  return false;
}

bool JUNOStarterESD::muon_veto(JM::EvtNavigator *nav) {
  auto ts_ref = nav->TimeStamp();
  uint64_t t_ref = ts_ref.GetSec() * 1'000'000'000ULL + ts_ref.GetNanoSec();
  bool is_muon = check_muon(nav);
  if (is_muon) {
    // 25 us
    m_t_muon = t_ref;
    if (m_t_muon > m_t_muon_prev + m_muon_veto_post) {
      m_segment_t0 = m_t_muon_prev + m_muon_veto_post;
      m_segment_t1 = m_t_muon;
      m_live_t += m_segment_t1 - m_segment_t0;
      m_segment_tree->Fill();
      ++m_segment;
    } else if (m_t_muon_prev == std::numeric_limits<uint64_t>::max()) {
      // first muon
      m_t_muon_prev = m_t_muon;
      m_segment = 0;
    } else {
      // consecutive muons
      m_t_muon_prev = m_t_muon;
    }
    // save muon
    if (m_save_muon) {
      m_dt_last_muon = safe_subtract(m_t_ref, m_t_muon_track);
      m_muon_npmt = npmt(nav, m_edm_path);
      m_muon_npe = npe(nav, m_edm_path);
      m_muon_track_tag = 0;
      m_muon_tree->Fill();
    }

    // muon tagged
    if (nav->getDetectorType() == 0) {
      ++m_n_skip_cd_muon;
    } else if (nav->getDetectorType() == 1) {
      ++m_n_skip_wp_muon;
    }
    return true;
  }

  if (t_ref < m_t_muon + m_muon_veto_post &&
      t_ref > m_t_muon) { // removed by post muon veto
    ++m_n_delayed;
    return true;
  }

  for (auto it = m_navBuf->current(); it != m_navBuf->end(); ++it) {
    auto *next = it->get();
    auto ts_next = next->TimeStamp();
    uint64_t t_next =
        ts_next.GetSec() * 1'000'000'000ULL + ts_next.GetNanoSec();
    if (t_next > t_ref + m_muon_veto_pre) {
      break;
    }
    if (check_muon(next)) {
      return true; // removed by pre muon veto
    }
  }
  return false;
}

bool JUNOStarterESD::mutiplicity_cut(JM::EvtNavigator *nav) {
  auto ts_ref = nav->TimeStamp();
  uint64_t t_ref = ts_ref.GetSec() * 1'000'000'000ULL + ts_ref.GetNanoSec();

  bool evt_cut_pre = false;
  bool evt_cut_post = false;

  // after search 25 us
  auto it = m_navBuf->current();
  if (it != m_navBuf->end()) {
    ++it;
    for (; it != m_navBuf->end(); ++it) {
      auto *next = it->get();
      if (!next)
        continue;
      if ((check_muon(next)) ||
          (next->getDetectorType() != 0)) { // skip muon events or
                                            // non-CD events
        continue;
      }
      auto ts_next = next->TimeStamp();
      uint64_t t_next =
          ts_next.GetSec() * 1'000'000'000ULL + ts_next.GetNanoSec();
      if (t_next < t_ref + 25'000) {
        evt_cut_post = true;
        break;
      }
      break;
    }
  } else // for end of the buffer
  {
    --it;
    for (; it != m_navBuf->begin(); --it) {
      auto *next = it->get();
      if (!next)
        continue;
      if ((check_muon(next)) ||
          (next->getDetectorType() != 0)) { // skip muon events or
                                            // non-CD events
        continue;
      }
      auto ts_next = next->TimeStamp();
      uint64_t t_next =
          ts_next.GetSec() * 1'000'000'000ULL + ts_next.GetNanoSec();
      if (t_next > t_ref - 25'000) {
        evt_cut_pre = true;
        break;
      }
      break;
    }
  }

  // forwards 25 us
  it = m_navBuf->current();
  if (it != m_navBuf->begin()) {
    --it;
    for (; it != m_navBuf->begin(); --it) {
      auto *pre = it->get();
      if (!pre)
        continue;
      if ((check_muon(pre)) ||
          (pre->getDetectorType() != 0)) { // skip muon events or
                                           // non-CD events
        continue;
      }
      auto ts_pre = pre->TimeStamp();
      uint64_t t_pre = ts_pre.GetSec() * 1'000'000'000ULL + ts_pre.GetNanoSec();
      if (t_pre > t_ref - 25'000) {
        evt_cut_pre = true;
        break;
      }
      break;
    }
  } else // for begin of the buffer
  {
    ++it;
    for (; it != m_navBuf->end(); ++it) {
      auto *pre = it->get();
      if (!pre)
        continue;
      if ((check_muon(pre)) ||
          (pre->getDetectorType() != 0)) { // skip muon events or
                                           // non-CD events
        continue;
      }
      auto ts_pre = pre->TimeStamp();
      uint64_t t_pre = ts_pre.GetSec() * 1'000'000'000ULL + ts_pre.GetNanoSec();
      if (t_pre < t_ref + 25'000) {
        evt_cut_post = true;
        break;
      }
      break;
    }
  }

  if (evt_cut_pre || evt_cut_post) { // mutiplicity cut
    ++m_n_skip_multi_r;
    return true;
  }
  return false;
}

bool JUNOStarterESD::get_trigger_rate(JM::EvtNavigator *nav) {
  auto ts_ref = nav->TimeStamp();
  uint64_t t_ref = ts_ref.GetSec() * 1'000'000'000ULL + ts_ref.GetNanoSec();
  uint64_t t_start = t_ref - 5'000'000;
  uint64_t t_end = t_ref + 5'000'000;
  uint64_t n_evt = 0;
  for (auto it = m_navBuf->begin(); it != m_navBuf->end(); ++it) {
    auto *next = it->get();
    if (!next)
      continue;
    if ((check_muon(next)) ||
        (next->getDetectorType() != 0)) { // skip muon events or
                                          // non-CD events
      continue;
    }
    auto ts_next = next->TimeStamp();
    uint64_t t_next =
        ts_next.GetSec() * 1'000'000'000ULL + ts_next.GetNanoSec();
    if (t_next > t_end) {
      break;
    }
    if (t_next >= t_start) {
      ++n_evt;
    }
  }
  m_trigger_rate = n_evt / 0.01; // 10 ms
  return true;
}

bool JUNOStarterESD::check_isotope() {
  float r = sqrt(m_x * m_x + m_y * m_y);
  if (m_pz > -0.75 && m_z > 12000 && m_z < 20000 && m_chisq_alt > 50 &&
      m_chisq_alt < 150 && r < 2000) {
    return true;
  }
  return false;
}

bool JUNOStarterESD::check_cluster() {
  float r = sqrt(m_x * m_x + m_y * m_y);
  if (m_pz < -0.75 && m_z > -5000 && m_z < 15000 && m_chisq_alt > 50 &&
      m_chisq_alt < 150 && r < 2000) {
    return true;
  }
  return false;
}

bool JUNOStarterESD::find_coincidence(JM::EvtNavigator *nav) {
  if (check_isotope()) {
    LogDebug << "Isotope event found" << std::endl;
    auto ts_ref = nav->TimeStamp();
    uint64_t t_ref = ts_ref.GetSec() * 1'000'000'000ULL + ts_ref.GetNanoSec();
    uint64_t t_coincidence = t_ref + 1'100'000; // 1100 us
    auto it = m_navBuf->current();
    if (it != m_navBuf->end()) {
      ++it;
      for (; it != m_navBuf->end(); ++it) {
        auto *next = it->get();
        if (!next) {
          continue;
        }
        if ((check_muon(next)) ||
            (next->getDetectorType() != 0)) { // skip muon events or
                                              // non-CD events
          continue;
        }
        auto ts_next = next->TimeStamp();
        uint64_t t_next =
            ts_next.GetSec() * 1'000'000'000ULL + ts_next.GetNanoSec();
        if (t_next > t_coincidence) {
          break;
        }
        if (t_next <= t_coincidence) {
          LogDebug << "Coincidence event found" << std::endl;
          auto *rec_coin = JM::getHeaderObject<JM::CdVertexRecHeader>(next);
          if (rec_coin == nullptr) {
            LogInfo << "No rec evt found" << std::endl;
            return true;
          }
          auto *rec_evt_coin =
              dynamic_cast<JM::CdVertexRecEvt *>(rec_coin->event());
          auto vertex = rec_evt_coin->getVertex(0);
          m_dt_coin = t_next - t_ref;
          m_x_coin = vertex->x();
          m_y_coin = vertex->y();
          m_z_coin = vertex->z();
          m_px_coin = vertex->px();
          m_py_coin = vertex->py();
          m_pz_coin = vertex->pz();
          m_chisq_alt_coin = vertex->chisq_alt();
          m_coincidence_tree->Fill();
          break;
        }
      }
    }
  }
  return true;
}

bool JUNOStarterESD::check_bkg() {
  auto r = sqrt(m_x * m_x + m_y * m_y + m_z * m_z);
  auto sig_condition = m_chisq_alt > 70 && m_chisq_alt < 180 && r < 14000 &&
                       m_pz > -0.75 && m_z < 10000;
  if (sig_condition) {
    return false;
  }
  return true;
}

bool JUNOStarterESD::get_recevt(JM::EvtNavigator *nav) {
  auto *rec = JM::getHeaderObject<JM::CdVertexRecHeader>(nav);
  if (rec == nullptr) {
    LogInfo << "No rec evt found" << std::endl;
    return true;
  }
  auto *rec_evt = dynamic_cast<JM::CdVertexRecEvt *>(rec->event());
  auto vertex = rec_evt->getVertex(0);

  // get cos theta between the event and the Sun
  auto solar_t = m_t_ref_correct / 1'000'000'000ULL;
  m_solar_angle.set_t(solar_t);
  m_solar_angle.calculate();
  auto sun_AltAz = std::array<double, 2>{m_solar_angle.get_sun_alt(),
                                         m_solar_angle.get_sun_az()};

  m_x = vertex->x();
  m_y = vertex->y();
  m_z = vertex->z();
  m_px = vertex->px();
  m_py = vertex->py();
  m_pz = vertex->pz();
  // m_evt_sec = ts_ref.GetSec();
  // m_evt_nsec = ts_ref.GetNanoSec();
  m_chisq_alt = vertex->chisq_alt();
  m_t0 = vertex->t0();
  m_raw_entry = m_cur_entry;

  auto ev_p = std::array<double, 3>{m_px, m_py, m_pz};

  m_cos_theta = m_jcs.get_cosTheta(ev_p, sun_AltAz);

  return true;
}

bool JUNOStarterESD::LoadMuonTrack() {
  // openm_muon_track_file
  TFile *track_file = TFile::Open(m_muon_track_file.c_str(), "READ");
  if (!track_file) {
    LogError << "Failed to open muon track file" << std::endl;
    return false;
  } else {
    LogInfo << "Open muon track file" << std::endl;
  }
  TTree *rec_info = (TTree *)track_file->Get("rec_info");
  if (!rec_info) {
    LogError << "Failed to get rec_info tree" << std::endl;
    return false;
  }

  Int_t fSec, fNanoSec;
  double init_x, init_y, init_z, exit_x, exit_y, exit_z;
  rec_info->SetBranchAddress("Trigger_Time.fSec", &fSec);
  rec_info->SetBranchAddress("Trigger_Time.fNanoSec", &fNanoSec);
  rec_info->SetBranchAddress("init_x", &init_x);
  rec_info->SetBranchAddress("init_y", &init_y);
  rec_info->SetBranchAddress("init_z", &init_z);
  rec_info->SetBranchAddress("exit_x", &exit_x);
  rec_info->SetBranchAddress("exit_y", &exit_y);
  rec_info->SetBranchAddress("exit_z", &exit_z);

  LogInfo << "Loading muon track info" << std::endl;

  auto muon_track_num = rec_info->GetEntries();

  LogInfo << "Muon track number: " << muon_track_num << std::endl;

  for (int i = 0; i < muon_track_num; ++i) {
    rec_info->GetEntry(i);
    gl_muon_track_map[fSec * 1'000'000'000ULL + fNanoSec] =
        std::make_tuple(init_x, init_y, init_z, exit_x, exit_y, exit_z);
  }

  return true;
}

bool JUNOStarterESD::find_muon_track(JM::EvtNavigator *nav,
                                     bool update_muon_track) {
  // high 64 bit recover time
  auto t_ref = nav->TimeStamp().GetSec() * 1'000'000'000ULL +
               nav->TimeStamp().GetNanoSec();

  // Find the closest key to t_ref
  auto lower = std::lower_bound(
      gl_muon_track_map.begin(), gl_muon_track_map.end(), t_ref,
      [](const auto &a, uint64_t b) { return compare_low48(a.first, b); });

  if (lower == gl_muon_track_map.end()) {
    // If lower is end(), then the map is empty or all keys are less than t_ref
    if (!gl_muon_track_map.empty()) {
      lower = --gl_muon_track_map.end();
    } else {
      LogError << "Muon track map is empty" << std::endl;
      return false;
    }
  }

  auto upper = lower;
  if (lower != gl_muon_track_map.begin()) {
    --lower;
  }

  // Determine which of the two candidates is closer to t_ref
  int64_t diff_lower = safe_subtract(t_ref, lower->first);
  int64_t diff_upper = safe_subtract(upper->first, t_ref);
  uint64_t closest_key =
      (diff_lower <= diff_upper) ? lower->first : upper->first;

  // Check if the closest key is within the time_range
  int64_t closest_diff = safe_subtract(closest_key, t_ref);
  if (std::abs(closest_diff) <= static_cast<int64_t>(25'000ULL)) {
    if (update_muon_track) {
      if (closest_key != m_t_muon_track) { // find next muon track in the map
        m_dt_last_muon = safe_subtract(closest_key, m_t_muon_track);
        m_t_muon_track = closest_key;
        if (m_muon_track_npe != 0 && m_save_muon) { // save last muon event
          m_muon_npe = m_muon_track_npe;
          m_muon_npmt = m_muon_track_npmt;
          m_muon_track_tag = 1;
          m_muon_tree->Fill();
        }
        m_muon_track_npmt = 0;
        m_muon_track_npe = 0;
        m_n_finded_muon_track++;
      }
      m_muon_track = gl_muon_track_map[closest_key];
      sw_muon_track_map[t_ref] = std::make_tuple(m_muon_track, closest_key);
      if (m_save_muon) {
        m_muon_track_npmt += npmt(nav, m_edm_path);
        m_muon_track_npe += npe(nav, m_edm_path);
      }
    }
    return true;
  }

  return false;
}

bool JUNOStarterESD::muon_track_veto(JM::EvtNavigator *nav) {
  auto ts_ref = nav->TimeStamp();
  uint64_t t_ref = ts_ref.GetSec() * 1'000'000'000ULL + ts_ref.GetNanoSec();
  bool is_muon_track = find_muon_track(nav, true);
  if (is_muon_track) {
    ++m_n_skip_muon_track;
    return true;
  }

  // get event rec info
  get_recevt(nav);

  // muon veto & track veto
  if (m_t_muon_track != 0) {
    // m_dt = safe_subtract(t_ref, m_t_muon_track);
    // auto evt_vertex = std::array<double, 3>{m_x, m_y, m_z};
    // m_d_track = m_jcs.get_track_distance(evt_vertex, m_muon_track);

    // whole detector veto
    if (m_dt < 2'000'000ULL &&                  // 2 ms
        static_cast<int64_t>(m_dt) > -25'000) { // removed by post muon veto
      ++m_n_delayed;
      return true;
    }

    // track veto
    for (auto it = m_navBuf->current(); it != m_navBuf->begin(); --it) {
      auto *next = it->get();
      auto ts_next = next->TimeStamp();
      uint64_t t_next =
          ts_next.GetSec() * 1'000'000'000ULL + ts_next.GetNanoSec();
      auto t_diff = t_ref - t_next;
      if (t_diff > 5'000'000'000ULL) {
        // 5s
        break;
      }
      if (sw_muon_track_map.find(t_next) != sw_muon_track_map.end()) {
        auto muon_track = sw_muon_track_map[t_next];
        auto track = std::get<0>(muon_track);
        auto t_track = std::get<1>(muon_track);
        auto ev_xyz = std::array<double, 3>{m_x, m_y, m_z};
        auto dt = safe_subtract(t_ref, t_track);
        auto dr = m_jcs.get_track_distance(ev_xyz, track);

        for (const auto &[min_d_track, max_d_track, max_dt] :
             track_veto_rules) {
          if (dr > min_d_track && dr < max_d_track &&
              dt < static_cast<int64_t>(max_dt) && dt > -25'000) {
            m_dt = dt;
            m_d_track = dr;
            ++m_n_skip_track_veto;
            if (m_save_track_veto) {
              if (!check_bkg() && !mutiplicity_cut(nav)) {
                m_muon_track_tree->Fill();
              }
            }
            return true;
          }
        }
      }
    }

    // // track veto
    // for (const auto &[min_d_track, max_d_track, max_dt] : track_veto_rules) {
    //   if (m_d_track > min_d_track && m_d_track < max_d_track &&
    //       m_dt<max_dt &&static_cast<int64_t>(m_dt)> - 25'000) {
    //     ++m_n_skip_track_veto;
    //     if (m_save_track_veto) {
    //       if (!check_bkg() && !mutiplicity_cut(nav)) {
    //         m_muon_track_tree->Fill();
    //       }
    //     }
    //     return true;
    //   }
    // }
  }

  return false;
}

// // calculate the dr and dt between the last muon track and current event
// if (m_muon_track_veto && m_save_track_veto) {
//   m_dt_last_muon_track = std::numeric_limits<int64_t>::max();
//   m_dr_last_muon_track = std::numeric_limits<float>::max();
//   for (auto it = m_navBuf->current(); it != m_navBuf->begin(); --it) {
//     auto *next = it->get();
//     auto ts_next = next->TimeStamp();
//     uint64_t t_next =
//         ts_next.GetSec() * 1'000'000'000ULL + ts_next.GetNanoSec();
//     auto t_diff = t_ref - t_next;
//     if (t_diff > 5'000'000'000ULL) {
//       // 200ms
//       break;
//     }
//     if (sw_muon_track_map.find(t_next) != sw_muon_track_map.end()) {
//       auto muon_track = sw_muon_track_map[t_next];
//       auto track = std::get<0>(muon_track);
//       auto t_track = std::get<1>(muon_track);
//       auto ev_xyz = std::array<double, 3>{m_x, m_y, m_z};
//       auto dt = safe_subtract(t_ref, t_track);
//       auto dr = m_jcs.get_track_distance(ev_xyz, track);
//       // if (dr < m_dr_last_muon_track) {
//       if (dr < 2'000) {
//         m_dr_last_muon_track = dr;
//         m_dt_last_muon_track = dt;
//         break;
//       }
//     }
//   }
// }