// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "ibrageometry.hh"

#include <fstream>
#include <nlohmann/json.hpp>

namespace Dune::IGA {

template <int dimworld, typename InputStringType>
requires(not std::convertible_to<std::string, InputStringType> and
         not std::convertible_to<InputStringType, const char*>)
auto readJson(InputStringType& ibraInputFile) {
  using json = nlohmann::json;

  std::vector<Ibra::Surface<dimworld>> surfaces;
  std::vector<Ibra::Curve2D> curves2D;

  std::vector<Ibra::BrepRepresentation> brepRepresentations;
  std::vector<Ibra::BrepLoopRepresentation> brepLoopRepresentations;
  std::vector<Ibra::BrepTrimRepresentation> brepTrimRepresentations;

  json ibraJson;
  try {
    ibraJson = json::parse(ibraInputFile);

    for (auto& j : ibraJson) {
      auto geo = j.get<Ibra::IbraBase>();

      switch (geo.type) {
        case Ibra::Type::NurbsSurfaceGeometry3D:
          surfaces.push_back(j.get<Ibra::Surface<dimworld>>());
          break;
        case Ibra::Type::NurbsCurveGeometry2D:
          curves2D.push_back(j.get<Ibra::Curve2D>());
          break;
        case Ibra::Type::BrepType:
          brepRepresentations.push_back(j.get<Ibra::BrepRepresentation>());
          break;
        case Ibra::Type::BrepLoopType:
          brepLoopRepresentations.push_back(j.get<Ibra::BrepLoopRepresentation>());
          break;
        case Ibra::Type::BrepTrimType:
          brepTrimRepresentations.push_back(j.get<Ibra::BrepTrimRepresentation>());
          break;
        default:
          // Do nothing -- keep the compiler happy
          break;
      }
    }
  } catch (json::parse_error& ex) {
    DUNE_THROW(Dune::IOError,
               "Error in reading input stream: " << ", parse error at byte " << ex.byte << " What: " << ex.what());
  }

  // Make Connections
  auto brepRepresentation = brepRepresentations[0];

  Ibra::Brep brep{brepRepresentation, curves2D, surfaces, brepLoopRepresentations, brepTrimRepresentations};
  return brep;
}
template <int dimworld>
auto readJson(const std::string& fileName) {
  std::ifstream ibraInputFile;
  ibraInputFile.open(fileName);
  return readJson<dimworld>(ibraInputFile);
}

} // namespace Dune::IGA
