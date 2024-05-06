// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once

namespace Dune::IGANEW::DefaultTrim {

/**
 * \brief Transformations for vertex and edge indices between Trimmer Notation and Dune Notation
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
  static U mapToTrimmer(T codim, U index) {
    if (codim == 2) /* vertex */ {
      return vertexIndexToTrimmer(index);
    }
    if (codim == 1) /* edge */ {
      return edgeIndexToTrimmer(index);
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
  static U vertexIndexToTrimmer(U index) {
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
  static U edgeIndexToTrimmer(U index) {
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
} // namespace Dune::IGANEW::DefaultTrim