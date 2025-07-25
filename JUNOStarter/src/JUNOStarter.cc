/*****************************************************************************/
// Author: Xuefeng DING <dingxf@ihep.ac.cn> @ IHEP-CAS
//
// Project: JUNOStarter
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
#include "JUNOStarter.h"
#include <Event/CdLpmtCalibEvt.h>
#include <Event/CdLpmtCalibHeader.h>
#include <Event/CdLpmtElecEvt.h>
#include <Event/CdLpmtElecHeader.h>
#include <Event/CdVertexRecEvt.h>
#include <Event/CdVertexRecHeader.h>
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
#include <TH1D.h>
#include <TTree.h>
#include <chrono>
#include <cstdint>
#include <sys/types.h>

DECLARE_ALGORITHM(JUNOStarter);

JUNOStarter::JUNOStarter(const std::string &name) : AlgBase(name) {
  declProp("FOO", m_foo);
}

bool JUNOStarter::initialize() {
  m_cur_entry = -1;
  getParent()->show();
  SniperPtr<RootWriter> svc(getParent(), "RootWriter");
  if (svc.invalid()) {
    LogError << "Can't locate RootWriter" << std::endl;
    return false;
  }
  m_tree = svc->bookTree(*getParent(), "USER_OUTPUT/bipo", "select pairs");
  m_tree->Branch("p_npe", &m_p_npe, "p_npe/F");
  m_tree->Branch("p_x", &m_p_x, "p_x/F");
  m_tree->Branch("p_y", &m_p_y, "p_y/F");
  m_tree->Branch("p_z", &m_p_z, "p_z/F");
  m_tree->Branch("d_npe", &m_d_npe, "d_npe/F");
  m_tree->Branch("d_x", &m_d_x, "d_x/F");
  m_tree->Branch("d_y", &m_d_y, "d_y/F");
  m_tree->Branch("d_z", &m_d_z, "d_z/F");
  m_tree->Branch("dt", &m_dt, "dt/i");
  m_t_prev = 0;

  LogDebug << "Initialized" << std::endl;
  return true;
}

bool JUNOStarter::execute() {
  auto t0 = std::chrono::high_resolution_clock::now();

  SniperDataPtr<JM::NavBuffer> navBuf(getParent(), "/Event");
  if (navBuf.invalid()) {
    LogError << "Can't locate data: /Event" << std::endl;
    return false;
  }

  ++m_cur_entry;
  auto *nav = navBuf->curEvt();
  if (nav->getDetectorType() != JM::EvtNavigator::DetectorType::CD) {
    return true;
  }
  auto *rechdr = JM::getHeaderObject<JM::CdVertexRecHeader>(nav);
  if (rechdr == nullptr) {
    return true;
  }
  auto *recevt = dynamic_cast<JM::CdVertexRecEvt *>(rechdr->event());
  if (!recevt) {
    LogFatal << "No event found in the header " << m_cur_entry << std::endl;
    return false;
  }

  auto ts_ref = nav->TimeStamp();
  uint64_t t_ref = ts_ref.GetSec() * 1'000'000'000ULL + ts_ref.GetNanoSec();

  if (t_ref < m_t_prev) {
    LogError << "Time order error: " << t_ref << " < " << m_t_prev << " at "
             << m_cur_entry << std::endl;
    return false;
  }

  for (auto it = navBuf->current(); it != navBuf->begin(); --it) {
    auto *prev = (it - 1)->get();
    if (prev->getDetectorType() != JM::EvtNavigator::DetectorType::CD) {
      continue;
    }
    auto *p_rechdr = JM::getHeaderObject<JM::CdVertexRecHeader>(prev);
    if (p_rechdr == nullptr) {
      return true;
    }
    auto *p_recevt = dynamic_cast<JM::CdVertexRecEvt *>(p_rechdr->event());
    if (!p_recevt) {
      LogFatal << "No event found in the header " << m_cur_entry << std::endl;
      return false;
    }
    auto ts_prev = prev->TimeStamp();
    uint64_t t_prev =
        ts_prev.GetSec() * 1'000'000'000ULL + ts_prev.GetNanoSec();
    int64_t dt = static_cast<int64_t>(t_ref - t_prev);
    if (dt <= 1'500'000) { // single event cut: 1.5 ms
                           // we found a coincidence !!
      m_p_npe = p_recevt->getVertex(0)->peSum();
      m_p_x = p_recevt->getVertex(0)->x();
      m_p_y = p_recevt->getVertex(0)->y();
      m_p_z = p_recevt->getVertex(0)->z();
      m_d_npe = recevt->getVertex(0)->peSum();
      m_d_x = recevt->getVertex(0)->x();
      m_d_y = recevt->getVertex(0)->y();
      m_d_z = recevt->getVertex(0)->z();
      m_dt = dt;
      m_tree->Fill();
      return true;
    }
  }

  return true;
}

bool JUNOStarter::finalize() {
  LogDebug << "finalizing" << std::endl;

  return true;
}
