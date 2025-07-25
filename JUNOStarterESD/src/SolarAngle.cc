/*****************************************************************************/
// Author: Xuefeng Ding <dingxf@ihep.ac.cn> @ IHEP-CAS
//
// Date: 2023 November 11
// Version: v1.0
// Description: SolarAngle
//
// All rights reserved. 2023 copyrighted.
/*****************************************************************************/
#include "SolarAngle.h"
#include "SolTrack.h"
#include <ctime>

SolarAngle::SolarAngle() : m_location(std::make_unique<Location>()), m_position(std::make_unique<Position>()) {}

SolarAngle::~SolarAngle() = default;

void SolarAngle::set_location(double latitude, double longitude) {
  // -700 m and ground differs by 1e-7 degree and is neglected
  m_location->longitude = longitude;
  m_location->latitude = latitude;
  m_location->pressure = 101.0325;   // NTP. for refraction correction, not used
  m_location->temperature = 293.15;  // NTP. for refraction correction, not used
}

void SolarAngle::set_t(double unix_tai) {
  // conver TAI time to UT time
  auto unix_t_dbl = unix_tai - 37;                     // good for 2017 -- 2030 within 1 sec.
  auto unix_t = static_cast<std::time_t>(unix_t_dbl);  // good for 2017 -- 2030 probably.
  auto unix_subs = unix_t_dbl - static_cast<double>(unix_t);
  std::tm unix_tm{};
  gmtime_r(&unix_t, &unix_tm);
  struct Time unix_time {
    unix_tm.tm_year + 1900, unix_tm.tm_mon + 1, unix_tm.tm_mday, unix_tm.tm_hour, unix_tm.tm_min,
        unix_tm.tm_sec + unix_subs
  };
  this->m_time_ut = std::make_unique<struct Time>(unix_time);
}

void SolarAngle::calculate() { SolTrack(*m_time_ut, *m_location, m_position.get(), 1, 1, 0, 0); }

double SolarAngle::get_sun_alt() const { return m_position->altitude; }

double SolarAngle::get_sun_az() const { return m_position->azimuthRefract; }

int SolarAngle::get_year() const { return m_time_ut->year; }
int SolarAngle::get_month() const { return m_time_ut->month; }
int SolarAngle::get_day() const { return m_time_ut->day; }
int SolarAngle::get_hour() const { return m_time_ut->hour; }
int SolarAngle::get_minutes() const { return m_time_ut->minute; }
double SolarAngle::get_seconds() const { return m_time_ut->second; }