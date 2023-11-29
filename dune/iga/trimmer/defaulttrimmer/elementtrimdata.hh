
#pragma once

namespace Dune {
  namespace IGANEW {
    namespace DefaultTrim {

      template <int mydim_, typename ScalarType>
      struct ElementTrimData {
        using ctype = ScalarType;

        static constexpr int mydimension = mydim_;
        [[nodiscard]] GeometryType type() const { return GeometryTypes::cube(mydimension); }
        using LocalCoordinate = FieldVector<ctype, mydimension>;

        bool checkInside(const LocalCoordinate& local) const {
          // TODO this functions could be a bit complicated basically we have to make sure the point lies inside the
          // outer boundary loop, but outside the inner loops thus we have to implement something as
          // https://en.wikipedia.org/wiki/Point_in_polygon#:~:text=One%20simple%20way%20of%20finding,an%20even%20number%20of%20times.
          //  maybe what we are searching for is already existing in Clipperlib
          //  https://angusj.com/clipper2/Docs/Units/Clipper/Functions/PointInPolygon.htm looks promising
          return true;
        }

        // outer loops
        // inner loops ....
        // triangulation...
      };

    }  // namespace Trim
  }    // namespace IGANEW
}  // namespace Dune
