/*****************************************************************************/
// Author: Xuefeng Ding <dingxf@ihep.ac.cn> @ IHEP-CAS
//
// Date: 2023 November 12
// Version: v1.0
// Description: JunoCoordinateSystem: JCS
//
// All rights reserved. 2023 copyrighted.
/*****************************************************************************/
#pragma once

#include <array>
class JunoCoordinateSystem {
public:
  explicit JunoCoordinateSystem(double north_JCS_phi = 180 - 56.7)
      : m_north_JCS_phi(north_JCS_phi) {}
  [[nodiscard]] std::array<double, 2>
  HCS_AltAz_to_JCS_ThetaPhi(double alt, double az) const;
  [[nodiscard]] double
  get_cosTheta(const std::array<double, 3> &ev_p,
               const std::array<double, 2> &sun_AltAz) const;
  [[nodiscard]] double
  get_hit_ev_direct_cosTheta(const std::array<double, 3> &ev_xyz,
                             const std::array<double, 3> &pmt,
                             const std::array<double, 3> &ev_pxyz) const;
  [[nodiscard]] double get_tof(const std::array<double, 3> &ev_xyz,
                               const std::array<double, 3> &pmt) const;
  [[nodiscard]] double get_distance(const std::array<double, 3> &a,
                                    const std::array<double, 3> &b) const;
  [[nodiscard]] double get_track_distance(
      const std::array<double, 3> &ev_xyz,
      const std::tuple<double, double, double, double, double, double> &track)
      const;
  [[nodiscard]]
  static std::array<double, 3>
  spherical_to_cartesian(const std::array<double, 2> &coord);
  static double dot(const std::array<double, 3> &a,
                    const std::array<double, 3> &b);

private:
  double m_north_JCS_phi;
};