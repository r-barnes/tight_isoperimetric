#pragma once

#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/Gps_circle_segment_traits_2.h>
#include <CGAL/IO/WKT.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Segment_Delaunay_graph_2.h>
#include <CGAL/Segment_Delaunay_graph_adaptation_policies_2.h>
#include <CGAL/Segment_Delaunay_graph_adaptation_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_traits_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Voronoi_diagram_2.h>

#include <CGAL/Lazy_exact_nt.h>

#include <stdexcept>

using K = CGAL::Exact_predicates_exact_constructions_kernel;
using Gt = CGAL::Segment_Delaunay_graph_traits_2<K>;
using SDG2 = CGAL::Segment_Delaunay_graph_2<Gt>;
using AT = CGAL::Segment_Delaunay_graph_adaptation_traits_2<SDG2>;
using AP = CGAL::Segment_Delaunay_graph_degeneracy_removal_policy_2<SDG2>;
using VoronoiDiagram = CGAL::Voronoi_diagram_2<SDG2, AT, AP>;
using Site_2 = AT::Site_2;
using Point_2 = AT::Point_2;
using Segment_2 = CGAL::Segment_2<K>;
using Locate_result = VoronoiDiagram::Locate_result;
using Vertex_handle = VoronoiDiagram::Vertex_handle;
using Face_handle = VoronoiDiagram::Face_handle;
using Halfedge_handle = VoronoiDiagram::Halfedge_handle;
using Ccb_halfedge_circulator = VoronoiDiagram::Ccb_halfedge_circulator;
using BHE_Iter = VoronoiDiagram::Bounded_halfedges_iterator;
using Halfedge = VoronoiDiagram::Halfedge;
using Vertex = VoronoiDiagram::Vertex;
using Polygon = CGAL::Polygon_with_holes_2<K>;
using MultiPolygon = std::deque<Polygon>;

using Circle_2 = K::Circle_2;
using Vector_2 = K::Vector_2;
using Traits_2 = CGAL::Gps_circle_segment_traits_2<K>;
using Polygon_set_2 = CGAL::General_polygon_set_2<Traits_2>;
using Polygon_2 = Traits_2::General_polygon_2;
using Polygon_with_holes_2 = Traits_2::General_polygon_with_holes_2;
using Curve_2 = Traits_2::Curve_2;
using X_monotone_curve_2 = Traits_2::X_monotone_curve_2;

class InterPointInterpolator {
 public:
  InterPointInterpolator(const Point_2 &a, const Point_2 &b) : a(a), b(b) {}
  // Returns points interpolated from a at t=0 to b at t=1
  Point_2 operator()(const double t) const {
    const auto m = b - a;
    return a + t * m;
  }
 private:
  const Point_2 a;
  const Point_2 b;
};

class CircularArcInterpolator {
 public:
  CircularArcInterpolator(const X_monotone_curve_2 &curve) : curve(curve) {
    if(!curve.is_circular()){
      throw std::runtime_error("CircularArcInterpolator requires a circular curve!");
    }
  }

  Point_2 operator()(const double t) const;

 private:
  X_monotone_curve_2 curve;
};

Polygon_2 construct_polygon(const Point_2 &pt, const K::FT squared_radius);

Polygon_2 construct_polygon (
  const Point_2& p1,
  const Point_2& p2,
  const Point_2& p3,
  const Point_2& p4
);

double area(const Polygon_2& P);
double area(const Polygon_with_holes_2& P);
double area(const Polygon_set_2& Pset);

K::FT point_to_line_segment_distance(const Point_2 &pt, const Segment_2 &seg);
K::FT point_to_point_distance(const Point_2 &a, const Point_2 &b);
K::FT point_to_site_distance(const Point_2 &pt, const Site_2 &site);

void print_polygon(const Polygon_with_holes_2 &pwh, std::string filename);
void print_polygon(const Polygon_set_2& Pset, std::string filename_prefix);