#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "geometry.hpp"

#include <cmath>
#include <list>

//For two circles of radii R and r and centered at (0,0) and (d,0) intersecting
//in a region shaped like an asymmetric lens.
constexpr double lens_area(const double r, const double R, const double d){
  return r*r*std::acos((d*d+r*r-R*R)/2/d/r) + R*R*std::acos((d*d+R*R-r*r)/2/d/R) - 0.5*std::sqrt((-d+r+R)*(d+r-R)*(d-r+R)*(d+r+R));
}

TEST_CASE("area of four circles joined by rectangle (1 hole)") {
  // Insert four non-intersecting circles.
  Polygon_set_2 S;
  const auto circ1 = construct_polygon(Point_2(1, 1), 1);
  const auto circ2 = construct_polygon(Point_2(5, 1), 1);
  const auto circ3 = construct_polygon(Point_2(5, 5), 1);
  const auto circ4 = construct_polygon(Point_2(1, 5), 1);
  S.join(circ1);
  S.join(circ2);
  S.join(circ3);
  S.join(circ4);
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

TEST_CASE("area of two non-overlapping circles") {
  Polygon_set_2 S;
  const auto circ1 = construct_polygon(Point_2(0, 0), 1);
  const auto circ2 = construct_polygon(Point_2(2, 2), 1);
  S.join(circ1);
  S.join(circ2);
  CHECK(S.number_of_polygons_with_holes()==2);
  CHECK(CGAL::to_double(area(S)) == doctest::Approx(2*M_PI));
}

TEST_CASE("area of two overlapping circles") {
  Polygon_set_2 S;
  const auto circ1 = construct_polygon(Point_2(0, 0), 1);
  const auto circ2 = construct_polygon(Point_2(0.5, 0), 1);
  S.join(circ1);
  S.join(circ2);
  CHECK(S.number_of_polygons_with_holes()==1);
  CHECK(CGAL::to_double(area(S)) == doctest::Approx(2*M_PI - lens_area(1, 1, 0.5)));
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