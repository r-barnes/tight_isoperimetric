// Compile with: clang++ -DBOOST_ALL_NO_LIB -DCGAL_USE_GMPXX=1 -O2 -g -DNDEBUG -Wall -Wextra -pedantic -march=native
// -frounding-math bob.cpp -lgmpxx -lmpfr -lgm

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/Gps_circle_segment_traits_2.h>
#include <CGAL/Lazy_exact_nt.h>

#include <list>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Circle_2 Circle_2;
typedef CGAL::Gps_circle_segment_traits_2<Kernel> Traits_2;
typedef CGAL::General_polygon_set_2<Traits_2> Polygon_set_2;
typedef Traits_2::General_polygon_2 Polygon_2;
typedef Traits_2::General_polygon_with_holes_2 Polygon_with_holes_2;
typedef Traits_2::Curve_2 Curve_2;
typedef Traits_2::X_monotone_curve_2 X_monotone_curve_2;

auto squared_distance(const Traits_2::Point_2& P1, const Traits_2::Point_2& P2) {
  const auto dx = P1.x() - P2.x();
  const auto dy = P1.y() - P2.y();
  return dx * dx + dy * dy;
}

auto area(const Polygon_2& P) {
  double res = 0.0;
  for (auto it = P.curves_begin(); it != P.curves_end(); ++it) {
    if (it->is_linear()) {
      const auto s = it->source();
      const auto t = it->target();
      res += CGAL::to_double((s.x() - t.x()) * (s.y() + t.y()) / 2);
    } else if (it->is_circular()) {
      const auto s = it->source();
      const auto t = it->target();
      res += CGAL::to_double((s.x() - t.x()) * (s.y() + t.y()) / 2);
      const auto ds           = CGAL::to_double(squared_distance(s, t));
      const auto rs           = CGAL::to_double(it->supporting_circle().squared_radius());
      const auto areaSector   = rs * std::asin(std::sqrt(ds) / (std::sqrt(rs) * 2));
      const auto areaTriangle = std::sqrt(ds) * std::sqrt(rs * 4 - ds) / 4;
      res += (areaSector - areaTriangle);
    }
  }
  return res;
}

auto area(const Polygon_with_holes_2& P) {
  auto res = area(P.outer_boundary());
  for (auto it = P.holes_begin(); it != P.holes_end(); ++it) {
    res += area(*it);
  }
  return res;
}

// Construct a polygon from a circle.
Polygon_2 construct_polygon(const Circle_2& circle) {
  // Subdivide the circle into two x-monotone arcs.
  Traits_2 traits;
  Curve_2 curve(circle);
  std::list<CGAL::Object> objects;
  traits.make_x_monotone_2_object()(curve, std::back_inserter(objects));
  CGAL_assertion(objects.size() == 2);
  // Construct the polygon.
  Polygon_2 pgn;
  X_monotone_curve_2 arc;
  std::list<CGAL::Object>::iterator iter;
  for (iter = objects.begin(); iter != objects.end(); ++iter) {
    CGAL::assign(arc, *iter);
    pgn.push_back(arc);
  }
  return pgn;
}

// Construct a polygon from a rectangle.
Polygon_2 construct_polygon(const Point_2& p1, const Point_2& p2, const Point_2& p3, const Point_2& p4) {
  Polygon_2 pgn;
  X_monotone_curve_2 s1(p1, p2);
  pgn.push_back(s1);
  X_monotone_curve_2 s2(p2, p3);
  pgn.push_back(s2);
  X_monotone_curve_2 s3(p3, p4);
  pgn.push_back(s3);
  X_monotone_curve_2 s4(p4, p1);
  pgn.push_back(s4);
  return pgn;
}

// The main program:
int main() {
  // Insert four non-intersecting circles.
  Polygon_set_2 S;
  Polygon_2 circ1, circ2, circ3, circ4;
  circ1 = construct_polygon(Circle_2(Point_2(1, 1), 1));
  S.insert(circ1);
  circ2 = construct_polygon(Circle_2(Point_2(5, 1), 1));
  S.insert(circ2);
  circ3 = construct_polygon(Circle_2(Point_2(5, 5), 1));
  S.insert(circ3);
  circ4 = construct_polygon(Circle_2(Point_2(1, 5), 1));
  S.insert(circ4);
  // Compute the union with four rectangles incrementally.
  Polygon_2 rect1, rect2, rect3, rect4;
  rect1 = construct_polygon(Point_2(1, 0), Point_2(5, 0), Point_2(5, 2), Point_2(1, 2));
  S.join(rect1);
  rect2 = construct_polygon(Point_2(1, 4), Point_2(5, 4), Point_2(5, 6), Point_2(1, 6));
  S.join(rect2);
  rect3 = construct_polygon(Point_2(0, 1), Point_2(2, 1), Point_2(2, 5), Point_2(0, 5));
  S.join(rect3);
  rect4 = construct_polygon(Point_2(4, 1), Point_2(6, 1), Point_2(6, 5), Point_2(4, 5));
  S.join(rect4);

  // Print the output.
  std::list<Polygon_with_holes_2> res;
  S.polygons_with_holes(std::back_inserter(res));
  std::copy(res.begin(), res.end(), std::ostream_iterator<Polygon_with_holes_2>(std::cout, "\n"));
  std::cout << std::endl;

  for (const auto& x : res) {
    std::cout << "Area = " << area(x) << std::endl;
  }

  return 0;
}
