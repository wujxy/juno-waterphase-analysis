/*****************************************************************************/
// Author: Xuefeng Ding <dingxf@ihep.ac.cn> @ IHEP-CAS
//
// Date: 2023 November 12
// Version: v1.0
// Description: JunoCoordinateSystem: JCS
//
// All rights reserved. 2023 copyrighted.
/*****************************************************************************/
#include "JunoCoordinateSystem.h"
#include <array>
#include <cmath>
#include <tuple>

static const double c = 299.792458; // mm/ns

std::array<double, 2>
JunoCoordinateSystem::HCS_AltAz_to_JCS_ThetaPhi(double alt, double az) const {
  auto theta = 90 - alt;           // HCS: 0 -> ground. JCS: 0 -> zenith
  auto phi = m_north_JCS_phi - az; // HCS: clockwise. JCS: anti-clockwise
  return {theta, phi};
}

double JunoCoordinateSystem::get_cosTheta(
    const std::array<double, 3> &ev_p,
    const std::array<double, 2> &sun_AltAz) const {
  auto sun_jcs_theta_phi =
      HCS_AltAz_to_JCS_ThetaPhi(sun_AltAz.at(0), sun_AltAz.at(1));
  auto sun_jcs_xyz = spherical_to_cartesian(sun_jcs_theta_phi);
  auto sun_direction = std::array<double, 3>{
      -sun_jcs_xyz.at(0), -sun_jcs_xyz.at(1), -sun_jcs_xyz.at(2)};
  auto cos_theta =
      dot(ev_p, sun_direction) / (std::sqrt(dot(ev_p, ev_p)) *
                                  std::sqrt(dot(sun_direction, sun_direction)));
  return cos_theta;
}

double JunoCoordinateSystem::get_hit_ev_direct_cosTheta(
    const std::array<double, 3> &ev_xyz, const std::array<double, 3> &pmt,
    const std::array<double, 3> &ev_pxyz) const {
  std::array<double, 3> ev_to_pmt = {pmt.at(0) - ev_xyz.at(0),
                                     pmt.at(1) - ev_xyz.at(1),
                                     pmt.at(2) - ev_xyz.at(2)};
  std::array<double, 3> ev_direct = ev_pxyz;
  auto cos_theta =
      dot(ev_direct, ev_to_pmt) /
      (sqrt(dot(ev_to_pmt, ev_to_pmt)) * sqrt(dot(ev_direct, ev_direct)));

  if (cos_theta > 1.0 || cos_theta < -1.0)
    return NAN;
  else
    return cos_theta;
}

double JunoCoordinateSystem::get_tof(const std::array<double, 3> &ev_xyz,
                                     const std::array<double, 3> &pmt) const {
  std::array<double, 3> ev_to_pmt = {pmt.at(0) - ev_xyz.at(0),
                                     pmt.at(1) - ev_xyz.at(1),
                                     pmt.at(2) - ev_xyz.at(2)};
  auto tof = sqrt(dot(ev_to_pmt, ev_to_pmt)) * 1.33 / c;
  return tof;
}

double
JunoCoordinateSystem::get_distance(const std::array<double, 3> &a,
                                   const std::array<double, 3> &b) const {
  auto dx = a.at(0) - b.at(0);
  auto dy = a.at(1) - b.at(1);
  auto dz = a.at(2) - b.at(2);
  return sqrt(dx * dx + dy * dy + dz * dz);
}

double JunoCoordinateSystem::get_track_distance(
    const std::array<double, 3> &ev_xyz,
    const std::tuple<double, double, double, double, double, double> &track)
    const {
  const double x1 = std::get<0>(track), y1 = std::get<1>(track),
               z1 = std::get<2>(track);
  const double x2 = std::get<3>(track), y2 = std::get<4>(track),
               z2 = std::get<5>(track);

  const double dx = x2 - x1, dy = y2 - y1, dz = z2 - z1;

  const double qx = ev_xyz[0] - x1, qy = ev_xyz[1] - y1, qz = ev_xyz[2] - z1;

  const double cross_x = qy * dz - qz * dy;
  const double cross_y = qz * dx - qx * dz;
  const double cross_z = qx * dy - qy * dx;

  const double cross_norm =
      cross_x * cross_x + cross_y * cross_y + cross_z * cross_z;
  const double p1p2_norm = dx * dx + dy * dy + dz * dz;
  return std::sqrt(cross_norm / p1p2_norm);
}

std::array<double, 3> JunoCoordinateSystem::spherical_to_cartesian(
    const std::array<double, 2> &coord) {
  auto theta = coord.at(0) / 180. * M_PI;
  auto phi = coord.at(1) / 180. * M_PI;
  return {std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi),
          std::cos(theta)};
}

double JunoCoordinateSystem::dot(const std::array<double, 3> &a,
                                 const std::array<double, 3> &b) {
  return a.at(0) * b.at(0) + a.at(1) * b.at(1) + a.at(2) * b.at(2);
}
