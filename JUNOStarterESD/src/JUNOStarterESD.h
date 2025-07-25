/*****************************************************************************/
// Author: Xuefeng DING <dingxf@ihep.ac.cn> @ IHEP-CAS
//
// Project: Processor
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
#pragma once
#include "RecTools/PMTTable.h"
#include "SniperKernel/AlgBase.h"
#include <EvtNavigator/EvtNavigator.h>
#include <EvtNavigator/NavBuffer.h>
#include <JunoCoordinateSystem.h>
#include <SolarAngle.h>
#include <cstdint>
#include <vector>
namespace JM {
class CdLpmtCalibEvt;
class WpCalibEvt;
class CdTriggerEvt;
class EvtNavigator;
} // namespace JM
class TTree;
class JUNOStarterESD : public AlgBase {
public:
  JUNOStarterESD(const std::string &name);

  bool initialize();
  bool execute();
  bool finalize();

private:
  static float npe(JM::EvtNavigator *nav, std::string edm_path);
  static int npmt(JM::EvtNavigator *nav, std::string edm_path);
  bool find_muon_track(JM::EvtNavigator *nav, bool update_muon_track);
  bool check_muon(JM::EvtNavigator *nav);
  bool check_isotope();
  bool check_cluster();
  bool get_trigger_rate(JM::EvtNavigator *nav);
  bool get_recevt(JM::EvtNavigator *nav);
  bool muon_track_veto(JM::EvtNavigator *nav);
  bool muon_veto(JM::EvtNavigator *nav);
  bool mutiplicity_cut(JM::EvtNavigator *nav);
  bool find_coincidence(JM::EvtNavigator *nav);
  bool check_bkg();
  bool initPMTPars();
  bool LoadMuonTrack();
  // helper fields
  std::string m_edm_path;
  JM::NavBuffer *m_navBuf = nullptr;
  SolarAngle m_solar_angle;
  JunoCoordinateSystem m_jcs;
  PMTTable m_pmtTable;
  unsigned int i_totPMTs;
  uint32_t m_cur_entry;
  TTree *m_tree;
  TTree *m_pmt_tree;
  TTree *m_muon_tree;
  TTree *m_veto_tree;
  TTree *m_segment_tree;
  TTree *m_selectrec_tree;
  TTree *m_coincidence_tree;
  TTree *m_cluster_tree;
  TTree *m_muon_track_tree;
  bool m_correlate_rtraw;
  // muon event
  bool m_save_muon;
  float m_muon_npe;
  int m_muon_track_tag;
  int64_t m_dt_last_muon;
  int m_muon_npmt;
  // rec event
  int m_nhits;
  int m_npmts;
  float m_npe;
  bool m_rec_file;
  float m_x;
  float m_y;
  float m_z;
  float m_px;
  float m_py;
  float m_pz;
  float m_chisq_alt;
  float m_t0;
  float m_cos_theta;
  int64_t m_dt_last_muon_track;
  float m_dr_last_muon_track;
  uint64_t m_t_ref;
  uint64_t m_t_pre;
  uint64_t m_t_ref_correct;
  uint32_t m_raw_entry;
  uint32_t m_evt_sec;
  uint32_t m_evt_nsec;
  uint32_t m_segment;
  uint64_t m_segment_t0;
  uint64_t m_segment_t1;
  uint32_t m_pmtid;
  uint32_t m_firednum;
  uint64_t m_t_prev;
  uint64_t m_t_muon;
  uint64_t m_t_muon_prev;
  bool m_muon_veto;
  uint64_t m_muon_veto_pre;
  uint64_t m_muon_veto_post;
  bool m_mutiplicity_cut;
  int64_t m_bias;
  uint64_t m_live_t;
  uint32_t m_n_delayed;
  uint32_t m_n_skip_cd_muon;
  uint32_t m_n_skip_wp_muon;
  uint32_t m_n_skip_empty;
  uint32_t m_n_skip_multi_r;
  uint32_t m_n_skip_periodic_tri;
  uint32_t m_E_rec_0_d5;
  uint32_t m_E_rec_m1_md8;
  uint32_t m_E_rec_d8_1;
  // flasher identification
  bool m_flasher_ide;
  float m_max_q_frac;
  float m_trigger_rate;
  float m_core_q_frac;
  float m_q_ratio_3m;
  float m_q_ratio_2m;
  float m_q_ratio_1m;
  float m_oppo_q_frac;
  // find coincidence
  bool m_coin_select;
  float m_x_coin;
  float m_y_coin;
  float m_z_coin;
  float m_px_coin;
  float m_py_coin;
  float m_pz_coin;
  float m_chisq_alt_coin;
  uint64_t m_dt_coin;
  bool m_bkg_filter;
  uint64_t m_dt_pre_muon;
  // muon track veto
  bool m_muon_track_veto;
  std::string m_muon_track_file;
  uint32_t m_n_skip_muon_track;
  uint64_t m_t_muon_track;
  uint64_t m_t_muon_track_prev;
  std::tuple<double, double, double, double, double, double> m_muon_track;
  uint64_t m_gl_fet_high16;
  int m_n_skip_track_veto;
  int m_n_finded_muon_track;
  bool m_save_track_veto;
  uint64_t m_dt;
  float m_d_track;
  float m_muon_track_npe;
  int m_muon_track_npmt;
};
