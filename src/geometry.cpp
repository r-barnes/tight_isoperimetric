#include "geometry.hpp"

#include "Arr_circle_segment_traits_2_custom.h"
#include <CGAL/squared_distance_2.h>

#include <iostream> // TODO: Remove
#include <string>

using ACST = CGAL::Arr_circle_segment_traits_2_custom<K>;

auto squared_distance(const Traits_2::Point_2& P1, const Traits_2::Point_2& P2){
  const auto dx = CGAL::to_double(P1.x()) - CGAL::to_double(P2.x());
  const auto dy = CGAL::to_double(P1.y()) - CGAL::to_double(P2.y());
  return dx * dx + dy * dy;
}

// ------ return signed area under the linear segment (P1, P2)
auto area(Traits_2::Point_2 const& P1, Traits_2::Point_2 const& P2)
{
  auto const dx = CGAL::to_double(P1.x()) - CGAL::to_double(P2.x());
  auto const sy = CGAL::to_double(P1.y()) + CGAL::to_double(P2.y());
  return dx * sy / 2;
}

// ------ return signed area under the circular segment (P1, P2, C)
auto area(Traits_2::Point_2 const& P1, Traits_2::Point_2 const& P2, Circle const& C)
{
  auto const dx = CGAL::to_double(P1.x()) - CGAL::to_double(P2.x());
  auto const dy = CGAL::to_double(P1.y()) - CGAL::to_double(P2.y());
  auto const squaredChord = dx * dx + dy * dy;
  auto const chord = std::sqrt(squaredChord);
  auto const squaredRadius = CGAL::to_double(C.squared_radius());
  auto const areaSector = squaredRadius * std::asin(std::min(1.0, chord / (std::sqrt(squaredRadius) * 2)));
  auto const areaTriangle = chord * std::sqrt(std::max(0.0, squaredRadius * 4 - squaredChord)) / 4;
  auto const areaCircularSegment = areaSector - areaTriangle;
  return area(P1, P2) + C.orientation() * areaCircularSegment;
}

// ------ return signed area under the X-monotone curve
auto area(X_monotone_curve_2 const& XCV)
{
  if (XCV.is_linear())
  {
    return area(XCV.source(), XCV.target());
  }
  else if (XCV.is_circular())
  {
    return area(XCV.source(), XCV.target(), XCV.supporting_circle());
  }
  else
  {
    return 0.0;
  }
}

// ------ return area of the simple polygon
auto area(Polygon_2 const& P)
{
  auto res = 0.0;
  for (auto it = P.curves_begin(); it != P.curves_end(); ++it) res += area(*it);
  return res;
}

// ------ return area of the polygon with (optional) holes
auto area(Polygon_with_holes_2 const& P)
{
  auto res = area(P.outer_boundary());
  for (auto it = P.holes_begin(); it != P.holes_end(); ++it) res += area(*it);
  return res;
}

double area(const Polygon_set_2& Pset){
  std::list<Polygon_with_holes_2> res;
  Pset.polygons_with_holes(std::back_inserter(res));
  double ret = 0;
  for(const auto &x : res){
    const auto temp = area(x);
    std::cout<<"Area = "<<temp<<std::endl;
    if(!std::isnan(temp)){
      ret += area(x);
    }
  }
  return ret;
}

/// Constructs a polygon from a circle
Polygon_2 construct_polygon(const Point_2 &pt, const K::FT squared_radius){
  // Subdivide the circle into two x-monotone arcs.
  Traits_2 traits;
  Curve_2 curve(Circle_2(pt, squared_radius));

  std::list<CGAL::Object> objects;
  traits.make_x_monotone_2_object()(curve, std::back_inserter(objects));
  CGAL_assertion(objects.size() == 2);

  // Construct the polygon.
  Polygon_2 pgn;
  X_monotone_curve_2 arc;
  for (auto iter = objects.begin(); iter != objects.end(); ++iter) {
    CGAL::assign(arc, *iter);
    pgn.push_back(arc);
  }

  return pgn;
}

// Construct a polygon from a rectangle.
Polygon_2 construct_polygon (
  const Point_2& p1,
  const Point_2& p2,
  const Point_2& p3,
  const Point_2& p4
){
  Polygon_2 pgn;
  X_monotone_curve_2 s1(p1, p2);    pgn.push_back(s1);
  X_monotone_curve_2 s2(p2, p3);    pgn.push_back(s2);
  X_monotone_curve_2 s3(p3, p4);    pgn.push_back(s3);
  X_monotone_curve_2 s4(p4, p1);    pgn.push_back(s4);
  return pgn;
}


/// Distance between a point and a line segment
K::FT point_to_line_segment_distance(const Point_2 &pt, const Segment_2 &seg){
  return CGAL::squared_distance(pt, seg);
}


/// Distance between two points
K::FT point_to_point_distance(const Point_2 &a, const Point_2 &b){
  return CGAL::squared_distance(a, b);
}


/// Distance between a point and a site
K::FT point_to_site_distance(const Point_2 &pt, const Site_2 &site){
  if (site.is_point()) {
    return point_to_point_distance(pt, site.point());
  } else {
    return point_to_line_segment_distance(pt, site.segment());
  }
}

void print_polygon(const Polygon_with_holes_2 &pwh, std::string filename){
  std::ofstream fout(filename);
  const auto &outer_boundary = pwh.outer_boundary();
  for(auto cit=outer_boundary.curves_begin();cit!=outer_boundary.curves_end();++cit){
    std::list<ACST::Approximate_point_2> output;
    ACST().approximate_2_object()(*cit, 0.1, std::back_inserter(output));
    for(const auto &x: output){
      fout<<x<<std::endl;
    }
    fout<<"#"<<std::endl;
  }
}

void print_polygon(const Polygon_set_2& Pset, std::string filename_prefix){
  std::vector<Polygon_with_holes_2> res;
  res.reserve(Pset.number_of_polygons_with_holes());
  Pset.polygons_with_holes(std::back_inserter(res));
  for(size_t i=0;i<res.size();i++){
    print_polygon(res.at(i), filename_prefix + "_" + std::to_string(i) + ".out");
  }
}

/*
Point_2 CircularArcInterpolator::operator()(const double t) const {
  const Vector_2 x_axis(1,0);
  const auto start_angle = 0; //CGAL::angle(start - circle.center(), x_axis);
  const auto ang_diff = 2*M_PI; //CGAL::angle(start, circle.center(), end);
  const auto cx = circle.center().x();
  const auto cy = circle.center().y();
  const auto radius = std::sqrt(CGAL::to_double(circle.squared_radius()));
  // auto px = start.x();
  // auto py = start.y();
  for(int i=0;i<STEPS;i++){
    const auto t = i / static_cast<double>(STEPS);
    const auto x = cx + radius * std::cos(start_angle + ang_diff*t);
    const auto y = cy + radius * std::sin(start_angle + ang_diff*t);
    // fout<<px<<","<<py<<","<<x<<","<<y<<std::endl;
    fout<<x<<","<<y<<std::endl;
    // px = x;
    // py = y;
  }
}


      if(cit->supporting_circle().orientation()==CGAL::COUNTERCLOCKWISE){
        const Point_2 S(CGAL::to_double(cit->source().x()),CGAL::to_double(cit->source().y()));
        const Point_2 T(CGAL::to_double(cit->target().x()),CGAL::to_double(cit->target().y()));
        // fout<<T.x()<<","<<T.y()<<std::endl;
        draw_arc(cit->supporting_circle(), S, T, fout);
      } else {
        //TODO
        std::cout<<"TODO"<<std::endl;
      }

  }
}*/