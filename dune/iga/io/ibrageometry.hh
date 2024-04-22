// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <nlohmann/json.hpp>

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/iga/geometrykernel/nurbspatchgeometry.hh>

namespace Dune::IGANEW::Ibra {

enum class Type
{
  NurbsCurveGeometry2D,
  NurbsSurfaceGeometry3D,
  BrepLoopType,
  BrepTrimType,
  BrepType,
  NoType
};

inline Type typeForTypeString(const std::string& typeString) {
  if (typeString == "NurbsCurveGeometry2D")
    return Type::NurbsCurveGeometry2D;
  else if (typeString == "NurbsSurfaceGeometry3D")
    return Type::NurbsSurfaceGeometry3D;
  else if (typeString == "BrepLoop")
    return Type::BrepLoopType;
  else if (typeString == "BrepTrim")
    return Type::BrepTrimType;
  else if (typeString == "Brep")
    return Type::BrepType;
  else
    return Type::NoType;
}

struct IbraBase
{
  std::string typeString;
  std::string key;
  Type type{Type::NoType};

  friend std::ostream& operator<<(std::ostream& stream, const IbraBase& geo) {
    stream << "Type: " << geo.typeString << ",\n Key: " << geo.key << "\n";
    return stream;
  }
};

template <int dim, int worldDim>
requires(dim < 3) && (worldDim >= dim)
struct IbraNURBSData : public IbraBase
{
  std::array<int, dim> degree{};
  std::array<int, dim> n_controlPoints{};
  std::array<std::vector<double>, dim> knots;
  Dune::DynamicMatrix<double> controlPoints;
  std::vector<double> weights;

  template <typename ControlPointType>
  [[nodiscard]] std::vector<std::vector<ControlPointType>> transformControlPoints() const {
    if constexpr (dim == 1) {
      std::vector<ControlPointType> cps;
      for (int i = 0; i < n_controlPoints[0]; ++i) {
        cps.push_back(controlPointAt<ControlPointType>({i}));
      }
      return {cps};
    } else if constexpr (dim == 2) {
      std::vector<std::vector<ControlPointType>> vec;
      std::vector<ControlPointType> cpTemp(n_controlPoints[1]);
      for (int i = 0; i < n_controlPoints[0]; ++i) {
        cpTemp.clear();
        for (int j = 0; j < n_controlPoints[1]; ++j) {
          cpTemp.push_back(controlPointAt<ControlPointType>({i, j}));
        }
        vec.push_back(cpTemp);
      }
      return vec;
    }
    __builtin_unreachable();
  }

  [[nodiscard]] std::array<std::vector<double>, dim> compileKnotVectors() const {
    std::array<std::vector<double>, dim> knotVec;

    // Insert first and last knot repeatedly
    for (int i = 0; i < dim; ++i) {
      knotVec[i].insert(knotVec[i].end(), knots[i][0]);
      knotVec[i].insert(knotVec[i].end(), knots[i].begin(), knots[i].end());
      knotVec[i].insert(knotVec[i].end(), knots[i].back());
    }
    return knotVec;
  }

private:
  template <typename ControlPointType>
  [[nodiscard]] ControlPointType controlPointAt(std::array<int, dim> idx) const {
    int row = (dim == 1) ? idx[0] : idx[0] * n_controlPoints[1] + idx[1];

    Dune::FieldVector<double, worldDim> p;
    for (int k = 0; k < worldDim; ++k)
      p[k] = controlPoints[row][k];

    double w = weights[row];

    return {p, w};
  }
};
template <int worldDim>
using Curve   = IbraNURBSData<1, worldDim>;
using Curve2D = IbraNURBSData<1, 2>;

template <int worldDim>
using Surface = IbraNURBSData<2, worldDim>;

struct BrepTrimRepresentation : IbraBase
{
  std::string geometry;
  Utilities::Domain<double> domain;
};

struct BrepTrim : IbraBase
{
  Curve2D geometry;
  Utilities::Domain<double> domain = {-1, -1};

  BrepTrim(BrepTrimRepresentation& trimRepresentation, std::vector<Curve2D>& allCurves)
      : domain{trimRepresentation.domain},
        IbraBase{static_cast<IbraBase&>(trimRepresentation)} {
    auto it = std::find_if(allCurves.begin(), allCurves.end(),
                           [trimRepresentation](auto x) { return (x.key == trimRepresentation.geometry); });
    if (it != allCurves.end())
      geometry = *it;
    else
      DUNE_THROW(Dune::InvalidStateException, "Couldn't find geometry in BrepTrim: " << key);
  }

  template <typename TrimmingCurve>
  auto asCurve() const {
    using ControlPointType    = typename TrimmingCurve::ControlPointType;
    using ControlPointNetType = typename TrimmingCurve::ControlPointNetType;
    using PatchData           = typename TrimmingCurve::PatchData;

    const auto cp = geometry.transformControlPoints<ControlPointType>()[0];
    std::array dimSize{static_cast<int>(cp.size())};
    std::array knotSpans{geometry.compileKnotVectors()};
    ControlPointNetType controlNet{dimSize, cp};

    return TrimmingCurve(PatchData(knotSpans, controlNet, std::array(geometry.degree)));
  }
};

struct BrepLoopRepresentation : IbraBase
{
  using StringVector = std::vector<std::string>;

  std::string brep;
  std::string face;
  StringVector trims;
};
struct BrepLoop : IbraBase
{
  std::vector<BrepTrim> trims;

  BrepLoop(BrepLoopRepresentation& loopRepresentation, std::vector<BrepTrim>& allTrims)
      : IbraBase{static_cast<IbraBase&>(loopRepresentation)} {
    for (auto& trimRepr : loopRepresentation.trims) {
      auto it = std::find_if(allTrims.begin(), allTrims.end(), [trimRepr](auto x) { return (x.key == trimRepr); });
      if (it != allTrims.end())
        trims.push_back(*it);
      else
        DUNE_THROW(Dune::InvalidStateException, "Couldn't find geometry in BrepTrim: " << key);
    }
  }
};

struct BrepRepresentation : IbraBase
{
  using StringVector = std::vector<std::string>;

  StringVector curve_geometries_2d;
  StringVector surface_geometries_3d;
  StringVector loops;
  StringVector trims;
};

template <int worldDim>
struct Brep : IbraBase
{
  std::vector<Curve2D> curves2D;
  using SurfaceType = Surface<worldDim>;
  std::vector<SurfaceType> surfaces;
  std::vector<BrepLoop> loops;
  std::vector<BrepTrim> trims;

  Brep(const BrepRepresentation& representation, std::vector<Curve2D>& allCurves, std::vector<SurfaceType>& allSurfaces,
       std::vector<Ibra::BrepLoopRepresentation>& allBrepLoopRepresentations,
       std::vector<Ibra::BrepTrimRepresentation>& allBrepTrimRepresentations)
      : IbraBase{static_cast<const BrepRepresentation&>(representation)} {
    // Lambda to check if a given key string is in a list of strings
    auto objectIsInList = [](const std::vector<std::string>& strRepr) {
      return [&](const auto& keyString) -> bool { return std::ranges::find(strRepr, keyString) != strRepr.end(); };
    };

    // Find in all Curves
    std::ranges::copy_if(allCurves, std::back_inserter(curves2D), objectIsInList(representation.curve_geometries_2d),
                         &Curve2D::key);
    // and all Surfaces
    std::ranges::copy_if(allSurfaces, std::back_inserter(surfaces),
                         objectIsInList(representation.surface_geometries_3d), &SurfaceType::key);

    // Find all Trims
    for (auto& trimRepr : allBrepTrimRepresentations) {
      if (objectIsInList(representation.trims)(trimRepr.key))
        trims.emplace_back(trimRepr, curves2D);
    }
    // All Loops
    for (auto& loopRepr : allBrepLoopRepresentations) {
      if (objectIsInList(representation.loops)(loopRepr.key))
        loops.emplace_back(loopRepr, trims);
    }

    // Asserts if all geometries were found
    assert(curves2D.size() == representation.curve_geometries_2d.size());
    assert(surfaces.size() == representation.surface_geometries_3d.size());
    assert(loops.size() == representation.loops.size());
    assert(trims.size() == representation.trims.size());
  }
};

///////
/// Json Read
//////
using json = nlohmann::json;

void getGenerics(const json& j, IbraBase& geometry) {
  j.at("type").get_to(geometry.typeString);
  j.at("key").get_to(geometry.key);

  geometry.type = typeForTypeString(geometry.typeString);
}

void from_json(const json& j, IbraBase& geometry) {
  getGenerics(j, geometry);
}

template <int worldDim>
void from_json(const json& j, Curve<worldDim>& curve) {
  getGenerics(j, curve);

  // Degree
  j.at("degree").get_to(curve.degree[0]);

  // Knots (u)
  std::ranges::copy(j.at("knots"), std::back_inserter(curve.knots[0]));

  // Poles aka Control Points
  j.at("nb_poles").get_to(curve.n_controlPoints[0]);

  curve.controlPoints.resize(curve.n_controlPoints[0], worldDim);

  for (int i = 0; i < curve.n_controlPoints[0]; ++i)
    for (int k = 0; k < worldDim; ++k)
      curve.controlPoints[i][k] = j.at("poles")[i][k];

  if (j.contains("weights"))
    std::ranges::copy(j.at("weights"), std::back_inserter(curve.weights));
  else {
    curve.weights.resize(curve.n_controlPoints[0]);
    std::ranges::fill(curve.weights, 1.0);
  }
}

template <int worldDim>
void from_json(const json& j, Surface<worldDim>& surface) {
  getGenerics(j, surface);

  // Degree
  j.at("degree_u").get_to(surface.degree[0]);
  j.at("degree_v").get_to(surface.degree[1]);

  // Knots (u,v)
  std::ranges::copy(j.at("knots_u"), std::back_inserter(surface.knots[0]));
  std::ranges::copy(j.at("knots_v"), std::back_inserter(surface.knots[1]));

  // Poles aka Control Points
  j.at("nb_poles_u").get_to(surface.n_controlPoints[0]);
  j.at("nb_poles_v").get_to(surface.n_controlPoints[1]);

  int total_controlPoints = surface.n_controlPoints[0] * surface.n_controlPoints[1];
  surface.controlPoints.resize(total_controlPoints, worldDim);

  for (int i = 0; i < total_controlPoints; ++i)
    for (int k = 0; k < worldDim; ++k)
      surface.controlPoints[i][k] = j.at("poles")[i][k];

  if (j.contains("weights"))
    std::ranges::copy(j.at("weights"), std::back_inserter(surface.weights));
  else {
    surface.weights.resize(total_controlPoints);
    std::ranges::fill(surface.weights, 1.0);
  }
}

void from_json(const json& j, BrepRepresentation& brep) {
  getGenerics(j, brep);

  j.at("curve_geometries_2d").get_to(brep.curve_geometries_2d);
  j.at("surface_geometries_3d").get_to(brep.surface_geometries_3d);

  j.at("loops").get_to(brep.loops);
  j.at("trims").get_to(brep.trims);
}

void from_json(const json& j, BrepLoopRepresentation& loop) {
  getGenerics(j, loop);

  j.at("brep").get_to(loop.brep);
  j.at("face").get_to(loop.face);
  j.at("trims").get_to(loop.trims);
}

void from_json(const json& j, BrepTrimRepresentation& trim) {
  getGenerics(j, trim);

  j.at("geometry").get_to(trim.geometry);

  auto dom            = j.at("domain");
  trim.domain.left()  = dom.front();
  trim.domain.right() = dom.back();
};

} // namespace Dune::IGANEW::Ibra

// DUNE_IGA_IBRAGEOMETRY_HH
