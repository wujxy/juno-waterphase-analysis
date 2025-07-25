/*****************************************************************************/
// Author: Xuefeng Ding <dingxf@ihep.ac.cn> @ IHEP-CAS
//
// Date: 2023 November 11
// Version: v1.0
// Description: SolarAngle
//
// All rights reserved. 2023 copyrighted.
/*****************************************************************************/
#pragma once
#include <memory>
struct Time;
struct Location;
struct Position;

class SolarAngle {
public:
  SolarAngle();
  SolarAngle(const SolarAngle &) = delete;
  SolarAngle(SolarAngle &&) = delete;
  SolarAngle &operator=(const SolarAngle &) = delete;
  SolarAngle &operator=(SolarAngle &&) = delete;
  ~SolarAngle();
  // height modifies height/1 AU and is small
  void set_location(double latitude, double longitude);
  void set_t(double unix_tai);
  void calculate();
  [[nodiscard]] int get_year() const;
  [[nodiscard]] int get_month() const;
  [[nodiscard]] int get_day() const;
  [[nodiscard]] int get_hour() const;
  [[nodiscard]] int get_minutes() const;
  [[nodiscard]] double get_seconds() const;
  [[nodiscard]] double get_sun_alt() const;
  [[nodiscard]] double get_sun_az() const;

private:
  std::unique_ptr<struct Time> m_time_ut;
  std::unique_ptr<struct Location> m_location;
  std::unique_ptr<struct Position> m_position;
};