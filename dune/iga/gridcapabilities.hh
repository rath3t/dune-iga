// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "dune/iga/utils/concepts.hh"
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/test/gridcheck.hh>

namespace Dune::IGA {
  template <int dim, int dimworld, typename ScalarType>
  class NURBSGrid;
}
namespace Dune::Capabilities {
  template <std::integral auto dim, std::integral auto dimworld, int codim, typename ScalarType>
  requires(dim <= 3) struct hasEntity<Dune::IGA::NURBSGrid<dim, dimworld, ScalarType>, codim> {
    static const bool v = true;
  };

  template <std::integral auto dim, std::integral auto dimworld, int codim, typename ScalarType>
  struct hasEntityIterator<Dune::IGA::NURBSGrid<dim, dimworld, ScalarType>, codim>
      : public hasEntity<Dune::IGA::NURBSGrid<dim, dimworld, ScalarType>, codim> {};

  template <std::integral auto dim, std::integral auto dimworld, int codim, typename ScalarType>
  struct canCommunicate<Dune::IGA::NURBSGrid<dim, dimworld, ScalarType>, codim> {
    static const bool v = false;
  };

  template <std::integral auto dim, std::integral auto dimworld, typename ScalarType>
  struct isLevelwiseConforming<Dune::IGA::NURBSGrid<dim, dimworld, ScalarType>> {
    static const bool v = true;
  };

  template <std::integral auto dim, std::integral auto dimworld, typename ScalarType>
  struct isLeafwiseConforming<Dune::IGA::NURBSGrid<dim, dimworld, ScalarType>> {
    static const bool v = true;
  };

  template <std::integral auto dim, std::integral auto dimworld, typename ScalarType>
  struct hasBackupRestoreFacilities<Dune::IGA::NURBSGrid<dim, dimworld, ScalarType>> {
    static const bool v = false;
  };

  template <std::integral auto dim, std::integral auto dimworld, typename ScalarType>
  struct threadSafe<Dune::IGA::NURBSGrid<dim, dimworld, ScalarType>> {
    static const bool v = false;
  };

  template <std::integral auto dim, std::integral auto dimworld, typename ScalarType>
  struct viewThreadSafe<Dune::IGA::NURBSGrid<dim, dimworld, ScalarType>> {
    static const bool v = false;
  };

  template <std::integral auto dim, std::integral auto dimworld, typename ScalarType>
  struct isCartesian<Dune::IGA::NURBSGrid<dim, dimworld, ScalarType>> {
    static const bool v = false;
  };

}  // namespace Dune::Capabilities

template <std::integral auto dim, std::integral auto dimworld, typename ScalarType>
struct Dune::EnableBoundarySegmentIndexCheck<Dune::IGA::NURBSGrid<dim, dimworld, ScalarType>> : public std::true_type {
};

template <std::integral auto dim, std::integral auto dimworld, typename ScalarType>
struct EnableLevelIntersectionIteratorCheck<Dune::IGA::NURBSGrid<dim, dimworld, ScalarType>> {
  static const bool v = true;
};
