nurbsgeometry.hh
- 

- vieles muss angepasst werden, integrationElement  ableitungen / 
- type() und referenceElement() none ausgeben


nurbsgrid.hh
- 

- int numBoundarySegments()
- int size
- int size(const GeometryType& type)

nurbsgridentity.hh
- 

- unsigned int subEntities
- bool hasBoundaryIntersections()
- subEntity(int i)
- auto type()
- auto referenceElement

nurbsgridindexsets.hh
- 

- auto& types(int codim)
- NURBSGridLeafIndexSet(NURBSLeafGridView<GridImpl> const& g)
  //TODO unique index for each geoemtry type
  // vertex 0 ---> n1
  // edges 0--> n2
  // elements of geometry::cube 0 --> n3, full elements
  // elements of geometry::none 0 --> n4, trimmed elements
  // indexset no gap in indices
  // for (ele : elements(gridView))
  //     if(ele.type()==Geometry::none)
  //         ele.IntegrationRule()
  //       else(ele.type()==Geometry::cube)
  //         Usual Gauss
  //    loops over trimmed and full elements
  //    each index should appear once


nurbsintersection.hh
- 

- //TODO dune buch Intersection chapter
-  GlobalCoordinate outerNormal ??
- std::size_t boundarySegmentIndex() consequtive numbering of boundary intersections


nurbsleafgridview.hh
- 

- int size(const GeometryType &type) (none case)


nurbspatch.hh
- 

- int patchBoundaryIndex(const int ocdim1Id)  simplify with map<int, int> ??
- 