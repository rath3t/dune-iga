//
// Created by Henri on 03.05.2023.
//

#pragma once

template <typename Name = std::string, typename Period = std::chrono::milliseconds>
class Timer {
  using Clock = std::chrono::high_resolution_clock;
  using Time = std::chrono::time_point<Clock>;

 public:
  void startTimer(Name name) {
    startTimes.emplace(name, Clock::now());
  }

  auto stopTimer(Name name) {
    auto stopTime = Clock::now();
    auto startTime = startTimes.at(name);
    startTimes.erase(name);

    return duration_cast<Period>(stopTime - startTime);
  }

  static std::string makeUniqueName(const std::string&& executableName, const std::string&& type = "") {
    using namespace std::chrono;
    std::string currentTime = fmt::format("_{}", std::chrono::system_clock::now());

    std::ranges::transform(currentTime, currentTime.begin(), [](char ch) {
      return (ch == ' ' or ch == ':') ? '_' : ch;
    });  // replace space and colon with underscore
    auto logFilename = executableName + currentTime + type;

    return logFilename;
  }

 private:
  std::unordered_map<Name, Time> startTimes;
};

