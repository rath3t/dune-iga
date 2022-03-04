#pragma once

#include "../Define.h"

#include "Ref.h"

#include <nlohmann/fifo_map.hpp>
#include <nlohmann/json.hpp>

#include "../Define.h"
#include "../Geometry/Interval.h"

#include <stdexcept>

namespace nlohmann {
  template <>
  struct adl_serializer<Eigen::Ref<const Eigen::VectorXd>>
  {
    template <typename TJson>
    static void to_json(TJson& json, const Eigen::Ref<const Eigen::VectorXd>& vector)
    {
      for (int i = 0; i < vector.size(); i++) {
        json.push_back(vector[i]);
      }
    }

    template <typename TJson>
    static void from_json(const TJson& json, Eigen::Ref<const Eigen::VectorXd>& vector)
    {
      assert(json.size() == vector.size());

      for (int i = 0; i < vector.size(); i++) {
        json.at(i).get_to(vector[i]);
      }
    }
  };
}

namespace nlohmann {
  template <>
  struct adl_serializer<anurbs::linear_algebra::Vector<2>>
  {
    template <typename TJson>
    static void to_json(TJson& json, const anurbs::linear_algebra::Vector<2>& point)
    {
      for (int i = 0; i < point.size(); i++) {
        json.push_back(point[i]);
      }
    }

    template <typename TJson>
    static void from_json(const TJson& json, anurbs::linear_algebra::Vector<2>& point)
    {
      assert(json.size() == point.size());

      for (int i = 0; i < point.size(); i++) {
        json.at(i).get_to(point[i]);
      }
    }
  };

  template <>
  struct adl_serializer<anurbs::linear_algebra::Vector<3>> {
    template <typename TJson>
    static void
    to_json(TJson& json, const anurbs::linear_algebra::Vector<3>& point)
    {
      for (int i = 0; i < point.size(); i++) {
        json.push_back(point[i]);
      }
    }

    template <typename TJson>
    static void from_json(const TJson& json, anurbs::linear_algebra::Vector<3>& point)
    {
      assert(json.size() == point.size());

      for (int i = 0; i < point.size(); i++) {
        json.at(i).get_to(point[i]);
      }
    }
  };
}

namespace nlohmann {
  template <>
  struct adl_serializer<Eigen::Ref<const Eigen::Matrix<double, Eigen::Dynamic, 2>>>
  {
    template <typename TJson>
    static void to_json(TJson& json, const Eigen::Ref<const Eigen::Matrix<double, Eigen::Dynamic, 2>>& points)
    {
      for (int i = 0; i < points.rows(); i++) {
        std::array<double, 2> array;
        for (int j = 0; j < 2; j++) {
          array[j] = points(i, j);
        }
        json.push_back(array);
      }
    }

    template <typename TJson>
    static void from_json(const TJson& json, Eigen::Ref<Eigen::Matrix<double, Eigen::Dynamic, 2>>& points)
    {
      assert(json.size() == points.rows());

      for (int i = 0; i < points.rows(); i++) {
        std::array<double, 2> array;
        json.at(i).get_to(array);
        for (int j = 0; j < 2; j++) {
          points(i, j) = array[j];
        }
      }
    }
  };

  template <>
  struct adl_serializer<Eigen::Ref<const Eigen::Matrix<double, Eigen::Dynamic, 3>>> {
    template <typename TJson>
    static void
    to_json(TJson& json, const Eigen::Ref<const Eigen::Matrix<double, Eigen::Dynamic, 3>>& points)
    {
      for (int i = 0; i < points.rows(); i++) {
        std::array<double, 3> array;
        for (int j = 0; j < 3; j++) {
          array[j] = points(i, j);
        }
        json.push_back(array);
      }
    }

    template <typename TJson>
    static void from_json(const TJson& json, Eigen::Ref<Eigen::Matrix<double, Eigen::Dynamic, 3>>& points)
    {
      assert(json.size() == points.rows());

      for (int i = 0; i < points.rows(); i++) {
        std::array<double, 3> array;
        json.at(i).get_to(array);
        for (int j = 0; j < 3; j++) {
          points(i, j) = array[j];
        }
      }
    }
  };
}

namespace anurbs {

  template <typename TData>
  void to_json(nlohmann::json& json, const Ref<TData>& ref)
  {
    if (ref.key().empty()) {
      throw std::runtime_error("Entity has no key");
    }

    json = ref.key();
  }

  void from_json(const nlohmann::json& json, Interval& interval)
  {
    assert(json.size() == 2);

    const auto t0 = json.at(0).get<double>();
    const auto t1 = json.at(1).get<double>();

    interval = Interval(t0, t1);
  }

  void to_json(nlohmann::json& json, const Interval& interval)
  {
    json.push_back(interval.t0());
    json.push_back(interval.t1());
  }

  template<class TKey, class TValue, class TCompare, class TAllocator>
  using FiFoMap = nlohmann::fifo_map<TKey, TValue,
                                     nlohmann::fifo_map_compare<TKey>, TAllocator>;

  using Json = nlohmann::basic_json<FiFoMap>;

  using ToJson = nlohmann::json;

  static std::string key_from_json(const Json& source)
  {
    const auto it = source.find("key");

    if (it == source.end()) {
      return {};
    }

    return it.value();
  }

  static std::string type_name_from_json(const Json& source)
  {
    const auto it = source.find("type");

    if (it == source.end()) {
      throw std::runtime_error("Type is missing");
    }

    return it.value();
  }

} // namespace anurbs