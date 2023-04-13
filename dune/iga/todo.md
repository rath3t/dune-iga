nurbsgeometry.hh
- 

- [x] `type()` anpassen, dass `none` ausgeben
- [ ] `referenceElement()` anpassen wg `type`
- [ ] `volume`
- [ ] `integrationElement` ??

nurbsgrid.hh
- 

- [ ] `int numBoundarySegments()`
- [x] `int size` via `Patch`
- [x] `int size(const GeometryType& type)` via `GridView`

nurbsgridentity.hh
- 

- [ ] `unsigned int subEntities`
- [ ] `bool hasBoundaryIntersections()`
- [x] `subEntity(int i)`
- [x] `auto type()`
- [x] `auto referenceElement` <span style="color:orange">Todo: Test if function overloads as expected</span>.

nurbsgridindexsets.hh
- 

- [x] `auto& types(int codim)`
- [ ] `NURBSGridLeafIndexSet(NURBSLeafGridView<GridImpl> const& g)`


nurbsintersection.hh
- 

- [ ] Siehe dune buch Intersection chapter
- [ ] `GlobalCoordinate outerNormal` ??
- [ ] `std::size_t boundarySegmentIndex()` consecutive numbering of boundary intersections


nurbsleafgridview.hh
- 

- [x] `int size(const GeometryType &type) (none case)` <span style="color:orange">Todo: Test</span>.


nurbspatch.hh
- 

- [ ] `int patchBoundaryIndex(const int ocdim1Id)`  simplify with mapToRange<int, int> ??
