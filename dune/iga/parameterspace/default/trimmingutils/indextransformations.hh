// SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

namespace Dune::IGA::DefaultParameterSpace {

/**
 * \brief Transformations for vertex and edge indices between ParameterSpace Notation and Dune Notation
 */
struct Transformations
{
  template <std::integral T, std::integral U>
  static U mapToDune(T codim, U index) {
    if (codim == 2) /* vertex */
      return vertexIndexToDune(index);

    if (codim == 1) /* edge */
      return edgeIndexToDune(index);

    DUNE_THROW(Dune::IOError, "Mapping only defined for codim 1 and 2");
  }

  template <std::integral T, std::integral U>
  static U mapToParameterSpace(T codim, U index) {
    if (codim == 2) /* vertex */ {
      return vertexIndexToParameterSpace(index);
    }
    if (codim == 1) /* edge */ {
      return edgeIndexToParameterSpace(index);
    }
    DUNE_THROW(Dune::IOError, "Mapping only defined for codim 1 and 2");
  }

private:
  template <std::integral U>
  static U vertexIndexToDune(U index) {
    switch (index) {
      case 0:
      case 1:
        return index;
      case 2:
        return 3;
      case 3:
        return 2;
      default:
        DUNE_THROW(Dune::IOError, "index out of bounds");
    }
  }

  template <std::integral U>
  static U edgeIndexToDune(U index) {
    switch (index) {
      case 0:
        return 2;
      case 1:
        return 1;
      case 2:
        return 3;
      case 3:
        return 0;
      default:
        DUNE_THROW(Dune::IOError, "index out of bounds");
    }
  }

  template <std::integral U>
  static U vertexIndexToParameterSpace(U index) {
    switch (index) {
      case 0:
      case 1:
        return index;
      case 2:
        return 3;
      case 3:
        return 2;
      default:
        DUNE_THROW(Dune::IOError, "index out of bounds");
    }
  }

  template <std::integral U>
  static U edgeIndexToParameterSpace(U index) {
    switch (index) {
      case 0:
        return 3;
      case 1:
        return 1;
      case 2:
        return 0;
      case 3:
        return 2;
      default:
        DUNE_THROW(Dune::IOError, "index out of bounds");
    }
  }
};
} // namespace Dune::IGA::DefaultParameterSpace