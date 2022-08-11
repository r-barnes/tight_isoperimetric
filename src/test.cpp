#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "geometry.hpp"

// #include <CGAL/Exact_predicates_exact_constructions_kernel.h>
// #include <CGAL/Gps_circle_segment_traits_2.h>
// #include <CGAL/General_polygon_set_2.h>
// #include <CGAL/Lazy_exact_nt.h>
#include <cmath>
#include <list>
// typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
// typedef Kernel::Point_2                                   Point_2;
// typedef Kernel::Circle_2                                  Circle_2;
// typedef CGAL::Gps_circle_segment_traits_2<Kernel>         Traits_2;
// typedef CGAL::General_polygon_set_2<Traits_2>             Polygon_set_2;
// typedef Traits_2::General_polygon_2                       Polygon_2;
// typedef Traits_2::General_polygon_with_holes_2            Polygon_with_holes_2;
// typedef Traits_2::Curve_2                                 Curve_2;
// typedef Traits_2::X_monotone_curve_2                      X_monotone_curve_2;
// Construct a polygon from a circle.

TEST_CASE("area of four circles joined by rectangle (1 hole)") {
  // Insert four non-intersecting circles.
  Polygon_set_2 S;
  const auto circ1 = construct_polygon(Point_2(1, 1), 1);
  const auto circ2 = construct_polygon(Point_2(5, 1), 1);
  const auto circ3 = construct_polygon(Point_2(5, 5), 1);
  const auto circ4 = construct_polygon(Point_2(1, 5), 1);
  S.insert(circ1);
  S.insert(circ2);
  S.insert(circ3);
  S.insert(circ4);
  // Compute the union with four rectangles incrementally.
  const auto rect1 = construct_polygon(Point_2(1, 0), Point_2(5, 0),
                                       Point_2(5, 2), Point_2(1, 2));
  const auto rect2 = construct_polygon(Point_2(1, 4), Point_2(5, 4),
                                       Point_2(5, 6), Point_2(1, 6));
  const auto rect3 = construct_polygon(Point_2(0, 1), Point_2(2, 1),
                                       Point_2(2, 5), Point_2(0, 5));
  const auto rect4 = construct_polygon(Point_2(4, 1), Point_2(6, 1),
                                       Point_2(6, 5), Point_2(4, 5));
  S.join(rect1);
  S.join(rect2);
  S.join(rect3);
  S.join(rect4);

  CHECK(S.number_of_polygons_with_holes()==1);
  CHECK(CGAL::to_double(area(S)) == doctest::Approx(4 * 8 - 4 + M_PI));
}

TEST_CASE("point to point distance"){
  const auto distance = point_to_point_distance(Point_2(1,3), Point_2(5,7));
  CHECK(CGAL::to_double(distance) == doctest::Approx(std::pow(1-5,2)+std::pow(3-7,2)));
}

TEST_CASE("Interpoint Interpolator"){
  InterPointInterpolator ipi(Point_2(0, 0), Point_2(10, 5));
  for(int i=0;i<=10;i++){
    const auto interpolated_point = ipi(i/10.0);
    CHECK(CGAL::to_double(interpolated_point.x())==i);
    CHECK(CGAL::to_double(interpolated_point.y())==i/2.0);
  }
}