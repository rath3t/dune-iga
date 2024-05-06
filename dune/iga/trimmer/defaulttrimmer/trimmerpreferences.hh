// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <mutex>

class Preferences
{
public:
  // Get the singleton instance
  static Preferences& getInstance() {
    static Preferences instance; // Guaranteed to be destroyed and initialized on first use
    return instance;
  }

  // Disable copy constructor and assignment operator
  Preferences(const Preferences&)            = delete;
  Preferences& operator=(const Preferences&) = delete;

  int boundaryDivisions() {
    std::lock_guard lock(mtx);
    return boundaryDivisions_;
  }

  void boundaryDivisions(int _boundaryDivisions) {
    std::lock_guard lock(mtx);
    this->boundaryDivisions_ = _boundaryDivisions;
  }

private:
  // Private constructor for Singleton pattern
  Preferences()
      : boundaryDivisions_(7) {}

  std::mutex mtx; // Mutex for thread safety
  int boundaryDivisions_;
};