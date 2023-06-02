// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "dune/iga/utils/concepts.hh"
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/test/gridcheck.hh>

namespace Dune::IGA {
  template <int dim, int dimworld, LinearAlgebra NurbsGridLinearAlgebraTraits>
  class NURBSGrid;
}
namespace Dune::Capabilities {
  template <std::integral auto dim, std::integral auto dimworld, int codim,
            Dune::IGA::LinearAlgebra NurbsGridLinearAlgebraTraits>
  requires(dim <= 3) struct hasEntity<Dune::IGA::NURBSGrid<dim, dimworld, NurbsGridLinearAlgebraTraits>, codim> {
    static const bool v = true;
  };

  template <Dune::IGA::LinearAlgebra NurbsGridLinearAlgebraTraits>
  struct hasEntity<Dune::IGA::NURBSGrid<2, 2, NurbsGridLinearAlgebraTraits>, 1> {
    static const bool v = true;
  };

  template <Dune::IGA::LinearAlgebra NurbsGridLinearAlgebraTraits>
  struct hasEntity<Dune::IGA::NURBSGrid<2, 2, NurbsGridLinearAlgebraTraits>, 2> {
    static const bool v = true;
  };
}  // namespace Dune::Capabilities

template <std::integral auto dim, std::integral auto dimworld, Dune::IGA::LinearAlgebra NurbsGridLinearAlgebraTraits>
struct Dune::EnableBoundarySegmentIndexCheck<Dune::IGA::NURBSGrid<dim, dimworld, NurbsGridLinearAlgebraTraits>>
    : public std::true_type {};

template <std::integral auto dim, std::integral auto dimworld, Dune::IGA::LinearAlgebra NurbsGridLinearAlgebraTraits>
struct EnableLevelIntersectionIteratorCheck<Dune::IGA::NURBSGrid<dim, dimworld, NurbsGridLinearAlgebraTraits>> {
  static const bool v = true;
};
