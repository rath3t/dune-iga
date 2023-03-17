//
// Created by Henri on 12.03.2023.
//

#ifndef DUNE_IGA_IBRAGEOMETRY_HH
#define DUNE_IGA_IBRAGEOMETRY_HH

#include "nurbspatchdata.hh"

#include <nlohmann/json.hpp>

namespace Dune::IGA::Ibra {

  enum Type {
    NurbsCurveGeometry2D,
    NurbsSurfaceGeometry3D,
    BrepLoopType,
    BrepTrimType,
    BrepType,
    NoType
  };

  Type typeForTypeString(const std::string& typeString) {
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

  class Geometry {
   public:
    std::string typeString;
    std::string key;
    Type type{Type::NoType};

    friend std::ostream& operator<<(std::ostream& stream, const Geometry& geo) {
      stream << "Type: " << geo.typeString << ",\n Key: " << geo.key << "\n";
      return stream;
    }
  };

  // The surface is inherently 3D, you can request ControlPoints with worldDim 2 or 3
  struct Surface : Geometry {
    std::array<int, 2> degree{};
    std::array<int, 2> n_controlPoints{};
    std::array<std::vector<double>, 2> knots;
    Dune::DynamicMatrix<double> controlPoints;

    template <int outputDim> requires(outputDim == 2 || outputDim == 3)
    std::vector<std::vector<typename NURBSPatchData<2, outputDim>::ControlPointType>> compileControlPoints() {

      auto cpAt = [this](int u, int v) -> Dune::FieldVector<double, outputDim> {
        int idx = u * n_controlPoints[1] + v;

        Dune::FieldVector<double, outputDim> p;
        for (int k = 0; k < outputDim; ++k)
          p[k] = controlPoints[idx][k];

        return p;
      };
      using ControlPointType = NURBSPatchData<2, outputDim>::ControlPointType ;
      std::vector<std::vector<ControlPointType>> vec;

      std::vector<ControlPointType> cpTemp(n_controlPoints[1]);
      for (int i = 0; i < n_controlPoints[0]; ++i) {
        cpTemp.clear();
        for (int j = 0; j < n_controlPoints[1]; ++j) {
          cpTemp.push_back({.p = cpAt(i, j), .w = 1});
        }
        vec.push_back(cpTemp);
      }
      return vec;
    }

    std::array<std::vector<double>, 2> compileKnotVectors() {
      std::array<std::vector<double>, 2> knotVec;

      // Insert first and last knot repeatedly
      std::vector<double> subKnotVector;
      for (int i = 0; i < 2; ++i) {
        subKnotVector.clear();
        subKnotVector.insert(subKnotVector.end(), knots[i][0]);
        subKnotVector.insert(subKnotVector.end(), knots[i].begin(), knots[i].end());
        subKnotVector.insert(subKnotVector.end(), knots[i].back());

        knotVec[i] = subKnotVector;
      }
      return knotVec;
    }

  };

  struct Curve : Geometry {
   public:
    static constexpr int dim = 2;
    int degree               = 0;
    int n_controlPoints      = 0;
    std::vector<double> knots;
    DynamicMatrix<double> controlPoints;

    std::vector<double> weights;

    std::vector<typename NURBSPatchData<1, 2>::ControlPointType> compileControlPoints() {
      constexpr int worldDim {2};

      auto cpAt = [this](int j) -> Dune::FieldVector<double, worldDim> {
        Dune::FieldVector<double, worldDim> p;
        for (int k = 0; k < worldDim; ++k)
          p[k] = controlPoints[j][k];
        return p;
      };
      using ControlPointType = NURBSPatchData<1, 2>::ControlPointType;
      std::vector<ControlPointType> vec;
      for (int i = 0; i < n_controlPoints; ++i) {
        vec.push_back({.p = cpAt(i), .w = weights[i]});
      }

      return vec;
    }
    // This might be unused and maybe should go away
    template <typename T> requires std::is_floating_point_v<T>
    std::vector<Dune::FieldVector<T, 2>> compileControlPoints() {
      constexpr int outputDim {2};

      auto cpAt = [this](int u) -> Dune::FieldVector<T, outputDim> {
        Dune::FieldVector<T, outputDim> p;
        for (int k = 0; k < outputDim; ++k)
          p[k] = controlPoints[u][k];

        return p;
      };

      std::vector<Dune::FieldVector<T, outputDim>> vec;
      for (int i = 0; i < n_controlPoints; ++i) {
        vec.push_back(cpAt(i));
      }

      return vec;
    }

    std::vector<double> compileKnotVectors() {
      std::vector<double> knotVec;

      // Insert first and last knot repeatedly
      knotVec.insert(knotVec.end(), knots[0]);
      knotVec.insert(knotVec.end(), knots.begin(), knots.end());
      knotVec.insert(knotVec.end(), knots.back());

      return knotVec;
    }
  };
  using Curve2D = Curve;

  void getGenerics(const Geometry& from, Geometry& to) {
    to.key        = from.key;
    to.typeString = from.typeString;
    to.type       = from.type;
  }

  struct BrepTrimRepresentation : Geometry {
    std::string geometry;
    std::array<double, 2> domain;
  };

  struct BrepTrim : Geometry {
    Curve2D geometry;
    std::array<double, 2> domain = {-1, -1};

    BrepTrim(BrepTrimRepresentation& trimRepresentation, std::vector<Curve2D>& allCurves) {
      getGenerics(trimRepresentation, *this);

      domain = trimRepresentation.domain;

      auto it = std::find_if(allCurves.begin(), allCurves.end(),
                             [trimRepresentation](auto x) { return (x.key == trimRepresentation.geometry); });
      if (it != allCurves.end())
        geometry = *it;
      else
        std::cerr << "Couldn't find geometry in BrepTrim: " << key << std::endl;
    }
  };

  struct BrepLoopRepresentation : Geometry {
    using StringVector = std::vector<std::string>;

    std::string brep;
    std::string face;
    StringVector trims;
  };
  struct BrepLoop : Geometry {
    std::vector<BrepTrim> trims;

    BrepLoop(BrepLoopRepresentation& loopRepresentation, std::vector<BrepTrim>& allTrims) {
      getGenerics(loopRepresentation, *this);

      for (auto& trimRepr : loopRepresentation.trims) {
        auto it = std::find_if(allTrims.begin(), allTrims.end(), [trimRepr](auto x) { return (x.key == trimRepr); });
        if (it != allTrims.end())
          trims.push_back(*it);
        else
          DUNE_THROW(Dune::InvalidStateException, "Couldn't find geometry in BrepTrim: "<< key);
      }
    }
  };

  struct BrepRepresentation : Geometry {
    using StringVector = std::vector<std::string>;

    StringVector curve_geometries_2d;
    StringVector surface_geometries_3d;
    StringVector loops;
    StringVector trims;
  };

  struct Brep : Geometry {
    std::vector<Curve2D> curves2D;
    std::vector<Surface> surfaces;
    std::vector<BrepLoop> loops;
    std::vector<BrepTrim> trims;

    Brep(const BrepRepresentation& representation, std::vector<Curve2D>& allCurves, std::vector<Surface>& allSurfaces,
         std::vector<Ibra::BrepLoopRepresentation>& allBrepLoopRepresentations,
         std::vector<Ibra::BrepTrimRepresentation>& allBrepTrimRepresentations) {
      getGenerics(representation, *this);

      // Lambda Function to check if a given key string is in a list of strings
      auto objectIsInList = [](std::vector<std::string> strRepr, std::string& keyString) -> bool {
        auto it = std::find_if(strRepr.begin(), strRepr.end(), [keyString](auto x) { return (x == keyString); });
        return (it != strRepr.end());
      };

      // Find in all Curves
      for (auto& curve : allCurves) {
        if (objectIsInList(representation.curve_geometries_2d, curve.key)) curves2D.push_back(curve);
      }
      // and all Surfaces
      for (auto& surface : allSurfaces) {
        if (objectIsInList(representation.surface_geometries_3d, surface.key)) surfaces.push_back(surface);
      }
      // Find all Trims
      for (auto& trimRepr : allBrepTrimRepresentations) {
        if (objectIsInList(representation.trims, trimRepr.key)) trims.emplace_back(trimRepr, curves2D);
      }
      // All Loops
      for (auto& loopRepr : allBrepLoopRepresentations) {
        if (objectIsInList(representation.loops, loopRepr.key)) loops.emplace_back(loopRepr, trims);
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

  void getGenerics(const json& j, Geometry& geometry) {
    j.at("type").get_to(geometry.typeString);
    j.at("key").get_to(geometry.key);

    geometry.type = typeForTypeString(geometry.typeString);
  }

  void from_json(const json& j, Geometry& geometry) { getGenerics(j, geometry); }


  void from_json(const json& j, Curve& curve) {
    getGenerics(j, curve);

    // Degree
    j.at("degree").get_to(curve.degree);

    // Knots (u)
    const auto& knots = j.at("knots");

    curve.knots.resize(knots.size());

    for (size_t k = 0; k < knots.size(); ++k)
      curve.knots[k] = knots[k];

    // Poles aka Control Points
    j.at("nb_poles").get_to(curve.n_controlPoints);

    curve.controlPoints.resize(curve.n_controlPoints, 2);

    for (int i = 0; i < curve.n_controlPoints; ++i)
      for (int k = 0; k < 2; ++k)
        curve.controlPoints[i][k] = j.at("poles")[i][k];

    std::vector<double> unit_weights(curve.n_controlPoints);
    std::fill(unit_weights.begin(), unit_weights.end(), 1.0);

    // Access with default value
    auto weights  = j.value("weights", unit_weights);
    curve.weights = weights;
  }

  void from_json(const json& j, Surface& surface) {
    getGenerics(j, surface);

    // Degree
    j.at("degree_u").get_to(surface.degree[0]);
    j.at("degree_v").get_to(surface.degree[1]);

    // Knots (u,v)
    const auto& knots_u = j.at("knots_u");
    const auto& knots_v = j.at("knots_v");

    surface.knots[0].resize(knots_u.size());
    surface.knots[1].resize(knots_v.size());

    for (size_t k = 0; k < knots_u.size(); ++k)
      surface.knots[0][k] = knots_u[k];

    for (size_t k = 0; k < knots_v.size(); ++k)
      surface.knots[1][k] = knots_v[k];

    // Poles aka Control Points
    j.at("nb_poles_u").get_to(surface.n_controlPoints[0]);
    j.at("nb_poles_v").get_to(surface.n_controlPoints[1]);

    surface.controlPoints.resize(surface.n_controlPoints[0] * surface.n_controlPoints[1], 3);

    for (int i = 0; i < surface.n_controlPoints[0] * surface.n_controlPoints[1]; ++i)
      for (int k = 0; k < 3; ++k)
        surface.controlPoints[i][k] = j.at("poles")[i][k];
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
    j.at("domain").get_to(trim.domain);
  }

}  // namespace Dune::IGA::Ibra

#endif  // DUNE_IGA_IBRAGEOMETRY_HH
