// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once


#include <mapbox/earcut.hpp>


namespace Dune::IGA {
  template <int dim>
  struct TransformToSpan {
    std::array<double, dim> scaling_{};
    std::array<double, dim> offset_{};

    TransformToSpan() = default;
    explicit TransformToSpan(std::pair<std::array<double, dim>, std::array<double, dim>>& input)
        : scaling_(input.first), offset_(input.second) {};

    // TODO Concept
    template <typename VecType>
    [[nodiscard]] auto toLocal(const VecType& cp) const -> VecType {
      VecType local;
      for (int i = 0; i < dim; ++i) {
        local[i] = (cp[i] - offset_[i]) / scaling_[i];
        local[i] = std::clamp(local[i], 0.0, 1.0);
      }
      return local;
    }

    [[nodiscard]] auto toLocal(const Clipper2Lib::PointD& p) const -> Clipper2Lib::PointD {
      std::array<double, 2> cp{p.x, p.y};
      auto cpLoc = toLocal(cp);
      return {cpLoc[0], cpLoc[1]};
    }
  };

  [[nodiscard]] double calculateBoundaryLoopLength(const std::vector<Boundary>& loop, const TransformToSpan<2>& transformer, unsigned int div = 200) {
    Clipper2Lib::PathD polygon;
    polygon.reserve(div * loop.size());

    for (auto& boundary : loop) {
      auto path = boundary.path(div);
      for (Clipper2Lib::PointD point : path)
        polygon.emplace_back(transformer.toLocal(point));
    }
    return std::fabs(Clipper2Lib::Length(polygon));
  }

  using DomainType = Utilities::Domain<double>;
  struct DomainInformation {
    DomainInformation(const DomainType& d, int i) : domain{d}, localIndex{i} {}
    DomainType domain{};
    int localIndex{};
  };

  auto determineCurvedBoundaries(const std::vector<Boundary>& boundaries) -> std::vector<bool> {
    std::vector<bool> result;
    result.reserve(boundaries.size());
    std::ranges::transform(boundaries, std::back_inserter(result),
                           [](auto& boundary) { return boundary.degree() > 1; });

    return result;
  }

  auto splitBoundariesImpl(const std::vector<Boundary>& boundaries, int maxSubSample, double tolerance, const TransformToSpan<2>& transformer) {
    bool checkTarget = false;

    double targetLength = 0.0;
    if (tolerance not_eq 0.0) {
      checkTarget = true;
      targetLength = calculateBoundaryLoopLength(boundaries, transformer);
    }
    auto refineMap = determineCurvedBoundaries(boundaries);

    // Get all domains
    std::vector<DomainInformation> domains;
    domains.reserve((1 << maxSubSample) * boundaries.size());

    for (auto i : std::views::iota(0u, boundaries.size()))
      domains.emplace_back(boundaries[i].domain, i);

    std::vector<DomainInformation> tempDomainInfos;
    tempDomainInfos.reserve((1 << maxSubSample) * boundaries.size());

    for (auto _ : std::views::iota(0, maxSubSample)) {
      tempDomainInfos.clear();
      tempDomainInfos.insert(tempDomainInfos.end(), domains.begin(), domains.end());
      domains.clear();
      for (auto& domainInfo : tempDomainInfos) {
        bool refine = refineMap[domainInfo.localIndex];

        if (refine) {
          std::array<DomainType, 2> newDomains = Utilities::splitDomainInHalf(domainInfo.domain);
          domains.emplace_back(newDomains[0], domainInfo.localIndex);
          domains.emplace_back(newDomains[1], domainInfo.localIndex);
        } else
          domains.push_back(domainInfo);
      }
      // Check convergence
      if (checkTarget) {
        auto computeCurrentLength = [&](const std::vector<DomainInformation>& domainInfos) -> double {
          Clipper2Lib::PathD path;
          for (const auto& domainInfo : domainInfos) {
            auto u = transformer.toLocal(boundaries[domainInfo.localIndex].nurbsGeometry(domainInfo.domain.left()));
            path.emplace_back(u[0], u[1]);
          }
          return std::fabs(Clipper2Lib::Length(path, true));
        };
        auto actLength = computeCurrentLength(domains);
        if (std::fabs(actLength - targetLength) / targetLength < tolerance)
          break;
      }
    }

    // Create split outerBoundaries from DomainInformation
    std::vector<Boundary> newBoundaries;
    std::ranges::transform(domains, std::back_inserter(newBoundaries), [&boundaries](const auto& domain_info) -> Boundary {
      return {boundaries[domain_info.localIndex].nurbsGeometry, domain_info.domain};
    });
    return newBoundaries;
  }



  template <int dim>
  struct GridBoundarySegment : Dune::BoundarySegment<dim, dim, double> {
    explicit GridBoundarySegment(Boundary& _boundary, const auto& _transformer) : boundary(_boundary), transformer(_transformer) {}

    Dune::FieldVector<double, dim> operator()(const Dune::FieldVector<double, 1>& localI) const override {
      // u has to be mapped on the domain of 0 to 1
      const auto local = std::clamp(localI[0], 0.0, 1.0);
      double u         = Utilities::mapToRange(local, Utilities::Domain<double>{}, boundary.domain);
      return transformer.toLocal(boundary.nurbsGeometry(u));
    }
    TransformToSpan<dim> transformer;
    Boundary boundary;
  };
}


// Add support for Dune::FieldVector in Earcut
namespace mapbox::util {

  template <typename T>
  struct nth<0, Dune::FieldVector<T, 2>> {
    inline static auto get(const Dune::FieldVector<T, 2>& t) { return t[0]; };
  };

  template <typename T>
  struct nth<1, Dune::FieldVector<T, 2>> {
    inline static auto get(const Dune::FieldVector<T, 2>& t) { return t[1]; };
  };
}  // namespace mapbox::util
