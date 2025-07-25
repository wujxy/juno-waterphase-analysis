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
#include "SniperKernel/AlgBase.h"
class TTree;
class JUNOStarter : public AlgBase {
public:
  JUNOStarter(const std::string &name);

  bool initialize();
  bool execute();
  bool finalize();

private:
  // config fields
  std::string m_foo;

  // helper fields
  uint32_t m_cur_entry;
  TTree *m_tree;
  float m_p_npe;
  float m_p_x;
  float m_p_y;
  float m_p_z;
  float m_d_npe;
  float m_d_x;
  float m_d_y;
  float m_d_z;
  uint32_t m_dt;
  uint64_t m_t_prev;
};
