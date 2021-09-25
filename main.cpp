// Compile with: clang++ -DBOOST_ALL_NO_LIB -DCGAL_USE_GMPXX=1 -O2 -g -DNDEBUG -Wall -Wextra -pedantic -march=native -frounding-math main.cpp -lgmpxx -lmpfr -lgmp
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
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
#include <graph_lite.h>

#include <CGAL/Lazy_exact_nt.h>

#include <algorithm>
#include <cmath>
#include <deque>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <queue>
#include <set>
#include <stdexcept>
#include <unordered_set>

typedef CGAL::Exact_predicates_exact_constructions_kernel              K;
typedef CGAL::Segment_Delaunay_graph_traits_2<K>                       Gt;
typedef CGAL::Segment_Delaunay_graph_2<Gt>                             SDG2;
typedef CGAL::Segment_Delaunay_graph_adaptation_traits_2<SDG2>         AT;
typedef CGAL::Segment_Delaunay_graph_degeneracy_removal_policy_2<SDG2> AP;
typedef CGAL::Voronoi_diagram_2<SDG2, AT, AP>      VoronoiDiagram;
typedef AT::Site_2                                 Site_2;
typedef AT::Point_2                                Point_2;
typedef CGAL::Segment_2<K>                         Segment_2;
typedef VoronoiDiagram::Locate_result              Locate_result;
typedef VoronoiDiagram::Vertex_handle              Vertex_handle;
typedef VoronoiDiagram::Face_handle                Face_handle;
typedef VoronoiDiagram::Halfedge_handle            Halfedge_handle;
typedef VoronoiDiagram::Ccb_halfedge_circulator    Ccb_halfedge_circulator;
typedef VoronoiDiagram::Bounded_halfedges_iterator BHE_Iter;
typedef VoronoiDiagram::Halfedge                   Halfedge;
typedef VoronoiDiagram::Vertex                     Vertex;
typedef CGAL::Polygon_with_holes_2<K>              Polygon;
typedef std::deque<Polygon>                        MultiPolygon;

typedef K::Circle_2                                Circle_2;
typedef CGAL::Gps_circle_segment_traits_2<K>       Traits_2;
typedef CGAL::General_polygon_set_2<Traits_2>      Polygon_set_2;
typedef Traits_2::General_polygon_2                Polygon_2;
typedef Traits_2::General_polygon_with_holes_2     Polygon_with_holes_2;
typedef Traits_2::Curve_2                          Curve_2;
typedef Traits_2::X_monotone_curve_2               X_monotone_curve_2;


/// Creates a hash of a Point_2, used for making O(1) point lookups
// struct Point2Hash {
//   size_t operator()(const Point_2 &pt) const {
//     std::hash<double> hasher;
//     auto seed = hasher(pt.x());
//     // boost::hash_combine from https://stackoverflow.com/q/35985960/752843
//     seed ^= hasher(pt.y()) + 0x9e3779b9 + (seed<<6) + (seed>>2);
//     return seed;
//   }
// };

// typedef std::unordered_set<Point_2, Point2Hash> Point2_Set;
typedef std::set<Point_2> Point2_Set;
typedef std::map<Vertex_handle, int> VH_Int_Map;


/// Holds a more accessible description of the Voronoi diagram
struct MedialData {
  /// Map of vertices comprising the Voronoi diagram
  std::vector<Vertex_handle> vertex_handles;
  /// List of edges in the diagram (pairs of the vertices above)
  std::vector<std::pair<size_t, size_t>> edges;
  /// Medial axis up governor. 1:1 correspondance with edges above.
  std::vector<VoronoiDiagram::Delaunay_graph::Vertex_handle> ups;
  /// Medial axis down governor. 1:1 correspondance with edges above.
  std::vector<VoronoiDiagram::Delaunay_graph::Vertex_handle> downs;
};

struct MedialPoint {
  Point_2 pt;
  /// True if the point is a vertex of the Voronoi diagram
  bool vd_vertex;
  /// Index of MedialData.vertex_handles that starts the Voronoi edge this point
  /// was sampled from
  size_t start_handle;
  /// Index of MedialData.vertex_handles that ends the Voronoi edge this point
  /// was sampled from
  size_t end_handle;
  /// Distance between the point and its support
  K::FT distance;
};


typedef graph_lite::Graph<size_t, MedialPoint> MedialGraph;


class Interpolator {
 public:
  Interpolator(const Point_2 &a, const Point_2 &b) : a(a), b(b) {}
  // Returns points interpolated from a at t=0 to b at t=1
  Point_2 operator()(const double t) const {
    auto m = b - a;
    return a + t * m;
  }
 private:
  Point_2 a;
  Point_2 b;
};



/// Read well-known text from @p filename to obtain shape boundary
MultiPolygon get_wkt_from_file(const std::string& filename){
  std::ifstream fin(filename);
  MultiPolygon mp;
  CGAL::read_multi_polygon_WKT(fin, mp);

  if(mp.empty()){
    throw std::runtime_error("WKT file '" + filename + "' was empty!");
  }
  for(const auto &poly: mp){
    if(poly.outer_boundary().size()==0){
      throw std::runtime_error("WKT file '" + filename + "' contained a polygon without an outer boundary!");
    }
  }

  return mp;
}


/// Converts a MultiPolygon into its corresponding Voronoi diagram
VoronoiDiagram convert_mp_to_voronoi_diagram(const MultiPolygon &mp){
  VoronoiDiagram vd;

  const auto add_segments_to_vd = [&](const auto &poly){
    for(std::size_t i=0;i<poly.size();i++){
      // Modulus to close the loop
      vd.insert(Site_2::construct_site_2(poly[i], poly[(i+1)%poly.size()]));
    }
  };

  for(const auto &poly: mp){                    // For each polygon in MultiPolygon
    std::cout<<poly<<std::endl;                 // Print polygon to screen for debugging
    add_segments_to_vd(poly.outer_boundary());  // Add the outer boundary
    for(const auto &hole : poly.holes()){       // And any holes
      add_segments_to_vd(hole);
    }
  }

  if(!vd.is_valid()){
    throw std::runtime_error("Voronoi Diagram was not valid!");
  }

  return vd;
}


/// Find @p item in collection @p c or add it if not present.
/// Returns the index of `item`'s location
int find_or_add(VH_Int_Map &c, const Vertex_handle &item){
  // Map means we can do this in log(N) time
  if(c.count(item) == 0){
    c.emplace(item, c.size());
    return c.size() - 1;
  }

  return c.at(item);
}


/// Convert a map of <T, int> pairs to a vector of `T` ordered by increasing int
std::vector<Vertex_handle> map_to_ordered_vector(const VH_Int_Map &m){
  std::vector<std::pair<Vertex_handle, int>> to_sort(m.begin(), m.end());
  to_sort.reserve(m.size());
  std::sort(to_sort.begin(), to_sort.end(), [](const auto &a, const auto &b){
    return a.second < b.second;
  });

  std::vector<Vertex_handle> ret;
  ret.reserve(to_sort.size());
  std::transform(begin(to_sort), end(to_sort), std::back_inserter(ret),
    [](auto const& pair){ return pair.first; }
  );

  return ret;
}


/// Find vertex handles which are in the interior of the MultiPolygon
std::set<Vertex_handle> identify_vertex_handles_inside_mp(
  const VoronoiDiagram &vd,
  const MultiPolygon &mp
){
  // Used to accelerate interior lookups by avoiding Point-in-Polygon checks for
  // vertices we've already considered
  std::set<Vertex_handle> considered;
  // The set of interior vertices we are building
  std::set<Vertex_handle> interior;

  for (
      auto edge_iter = vd.bounded_halfedges_begin();
      edge_iter != vd.bounded_halfedges_end();
      edge_iter++
  ) {
    // Determine if an orientation implies an interior vertex
    const auto inside = [](const auto &orientation){
      return orientation == CGAL::ON_ORIENTED_BOUNDARY || orientation == CGAL::POSITIVE;
    };

    // Determine if a vertex is in the interior of the multipolygon and, if so,
    // add it to `interior`
    const auto vertex_in_mp_interior = [&](const Vertex_handle& vh){
      // Skip vertices which have already been considered, since a vertex may
      // be connected to multiple halfedges
      if(considered.count(vh)!=0){
        return;
      }
      // Ensure we don't look at a vertex twice
      considered.insert(vh);
      // Determine if the vertex is inside of any polygon of the MultiPolygon
      const auto inside_of_a_poly = std::any_of(
        mp.begin(), mp.end(), [&](const auto &poly) {
          return inside(CGAL::oriented_side(vh->point(), poly));
        }
      );
      // If the vertex was inside the MultiPolygon make a note of it
      if(inside_of_a_poly){
        interior.insert(vh);
      }
    };

    // Check both vertices of the current halfedge of the Voronoi diagram
    vertex_in_mp_interior(edge_iter->source());
    vertex_in_mp_interior(edge_iter->target());
  }

  return interior;
}


/// The medial axis is formed by building a Voronoi diagram and then removing
/// the edges of the diagram which connect to the concave points of the
/// MultiPolygon. Here, we identify those concave points
Point2_Set identify_concave_points_of_mp(const MultiPolygon &mp){
  Point2_Set concave_points;

  // Determine cross-product, given three points. The sign of the cross-product
  // determines whether the point is concave or convex.
  const auto z_cross_product = [](const Point_2 &pt1, const Point_2 &pt2, const Point_2 &pt3){
    const auto dx1 = pt2.x() - pt1.x();
    const auto dy1 = pt2.y() - pt1.y();
    const auto dx2 = pt3.x() - pt2.x();
    const auto dy2 = pt3.y() - pt2.y();
    return dx1 * dy2 - dy1 * dx2;
  };

  // Loop through all the points in a polygon, get their cross products, and
  // add any concave points to the set we're building.
  // `sense` should be `1` for outer boundaries and `-1` for holes, since holes
  // will have points facing outward.
  const auto consider_polygon = [&](const auto &poly, const double sense){
    for(size_t i=0;i<poly.size()+1;i++){
      const auto zcp = z_cross_product(
        poly[(i + 0) % poly.size()],
        poly[(i + 1) % poly.size()],
        poly[(i + 2) % poly.size()]
      );
      if(sense*zcp < 0){
        concave_points.insert(poly[(i + 1) % poly.size()]);
      }
    }
  };

  // Loop over the polygons of the MultiPolygon, as well as their holes
  for(const auto &poly: mp){
    // Outer boundary has positive sense
    consider_polygon(poly.outer_boundary(), 1);
    for(const auto &hole: poly.holes()){
      // Inner boundaries (holes) have negative (opposite) sense
      consider_polygon(hole, -1);
    }
  }

  return concave_points;
}


/// Print the points which collectively comprise the medial axis
void print_medial_axis_points(const MedialData &md, const std::string &filename){
  std::ofstream fout(filename);
  fout<<"x,y"<<std::endl;
  for (const auto &vh : md.vertex_handles) {
    fout << vh->point().x() << "," << vh->point().y() << std::endl;
  }
}


/// Prints the edges of the medial diagram
void print_medial_axis_edges(const MedialData &md, const std::string &filename){
  std::ofstream fout(filename);
  fout<<"SourceIdx,TargetIdx,UpGovernorIsPoint,DownGovernorIsPoint"<<std::endl;
  for (std::size_t i = 0; i < md.edges.size(); i++) {
    fout << md.edges[i].first        << ","
          << md.edges[i].second      << ","
          << md.ups[i]->is_point()   << "," // Is up-governor a point?
          << md.downs[i]->is_point()        // Is down-governor a point?
          << std::endl;
  }
}


MedialData filter_voronoi_diagram_to_medial_axis(
  const VoronoiDiagram &vd,
  const MultiPolygon &mp
){
  MedialData ret;
  VH_Int_Map handles;

  const auto interior = identify_vertex_handles_inside_mp(vd, mp);
  const auto concave_points = identify_concave_points_of_mp(mp);

  // Returns true if a point is a concave point of the MultiPolygon
  const auto pconcave = [&](const Point_2 &pt){
    return concave_points.count(pt) != 0;
  };

  // The Voronoi diagram is comprised of a number of vertices connected by lines
  // Here, we go through each edge of the Voronoi diagram and determine which
  // vertices it's incident on. We add these vertices to `ret.vertex_handles`
  // so that they will have unique ids.

  // The `up` and `down` refer to the medial axis governors - that which
  // constrains each edge of the Voronoi diagram
  for (
      auto edge_iter = vd.bounded_halfedges_begin();
      edge_iter != vd.bounded_halfedges_end();
      edge_iter++
  ) {
    const Halfedge& halfedge = *edge_iter;
    const Vertex_handle& v1p = halfedge.source();
    const Vertex_handle& v2p = halfedge.target();

    // Filter Voronoi diagram to only the part in the interior of the
    // MultiPolygon
    if(interior.count(v1p)==0 || interior.count(v2p)==0){
      continue;
    }

    // Drop those edges of the diagram which are not part of the medial axis
    if(pconcave(v1p->point()) || pconcave(v2p->point())){
      continue;
    }

    // Get unique ids for edge vertex handle that's part of the medial axis
    const auto id1 = find_or_add(handles, v1p);
    const auto id2 = find_or_add(handles, v2p);
    ret.edges.emplace_back(id1, id2);

    // Keep track of the medial axis governors
    ret.ups.push_back(halfedge.up());
    ret.downs.push_back(halfedge.down());
  }

  ret.vertex_handles = map_to_ordered_vector(handles);

  return ret;
}


/// Distance between a point and a line segment
/// From: https://stackoverflow.com/a/6853926/752843
auto point_to_line_segment_distance(const Point_2 &pt, const Segment_2 &seg){
  return CGAL::squared_distance(pt, seg);
}


/// Distance between two points
auto point_to_point_distance(const Point_2 &a, const Point_2 &b){
  return CGAL::squared_distance(a, b);
}


/// Distance between a point and a site
auto point_to_site_distance(const Point_2 &pt, const Site_2 &site){
  if (site.is_point()) {
    return point_to_point_distance(pt, site.point());
  } else {
    return point_to_line_segment_distance(pt, site.segment());
  }
}


/// Given 2 sites returns the minimum distance between the point and both of them
auto point_to_min_site_distance(const Point_2 &pt, const Site_2 &site1, const Site_2 &site2){
  return std::min(point_to_site_distance(pt, site1), point_to_site_distance(pt, site2));
}


MedialGraph get_medial_axis_graph(const MedialData &md){
  MedialGraph mg;

  size_t vertex_id = 0;

  for(size_t i=0;i<md.vertex_handles.size();i++){
    const auto &vh = md.vertex_handles.at(i);
    mg.add_node_with_prop(vertex_id, MedialPoint{
      .pt = vh->point(),
      .vd_vertex = true,
      .start_handle = vertex_id,
      .end_handle = vertex_id,
      .distance = -1  // We'll fix this below
    });
    vertex_id++;
  }

  // Fix distances to vertex handles
  for(size_t e=0;e<md.edges.size();e++){
    const auto &edge = md.edges.at(e);
    auto &spt = mg.node_prop(mg.find(edge.first));
    spt.distance = point_to_min_site_distance(
      spt.pt, md.ups.at(e)->site(), md.downs.at(e)->site()
    );
    auto &ept = mg.node_prop(mg.find(edge.second));
    ept.distance = point_to_min_site_distance(
      ept.pt, md.ups.at(e)->site(), md.downs.at(e)->site()
    );
  }

  // Get distances to all other points
  constexpr auto STEPS = 100;
  for(size_t e=0;e<md.edges.size();e++){
    const auto &edge = md.edges.at(e);
    const auto &svh = md.vertex_handles.at(edge.first);  // Start vertex
    const auto &evh = md.vertex_handles.at(edge.second); // End vertex
    Interpolator interp(svh->point(), evh->point());
    for(int t=1;t<STEPS;t++){
      // Location of the interpolated point
      const auto interp_point = interp(t/static_cast<double>(STEPS));

      mg.add_node_with_prop(vertex_id, MedialPoint{
        .pt = interp_point,
        .vd_vertex = false,
        .start_handle = edge.first,
        .end_handle = edge.second,
        .distance = point_to_min_site_distance(
          interp_point, md.ups.at(e)->site(), md.downs.at(e)->site()
        )
      });

      if(t==1){
        mg.add_edge(edge.first, vertex_id);
      } else if(t==STEPS-1){
        mg.add_edge(vertex_id, edge.second);
      } else {
        mg.add_edge(vertex_id - 1, vertex_id);
      }
      vertex_id++;
    }
  }

  return mg;
}


/// Prints the points of the medial graph and their distances from supports as CSV
void print_medial_graph_points_and_distances(const MedialGraph &mg, const std::string &filename){
  std::ofstream fout(filename);
  fout<<"x,y,distance"<<std::endl;
  for (const auto&& [node, node_prop, neighbors_view]: mg) {
    fout<<node_prop.pt.x()<<","<<node_prop.pt.y()<<","<<node_prop.distance<<std::endl;
  }
}


/// Constructs a polygon from a circle
Polygon_2 construct_polygon_circle(const Point_2 &pt, const K::FT radius){
  // Subdivide the circle into two x-monotone arcs.
  Traits_2 traits;
  Curve_2 curve (Circle_2(pt, radius));
  std::list<CGAL::Object>  objects;
  traits.make_x_monotone_2_object() (curve, std::back_inserter(objects));
  std::cout<<"asize = "<<objects.size()<<" "<<radius<<std::endl;
  CGAL_assertion (objects.size() == 2);
  // Construct the polygon.
  Polygon_2 pgn;
  X_monotone_curve_2 arc;
  std::list<CGAL::Object>::iterator iter;
  for (iter = objects.begin(); iter != objects.end(); ++iter) {
    CGAL::assign (arc, *iter);
    pgn.push_back (arc);
  }
  return pgn;
}


void iterate_graph(const MedialGraph &mg){
  //Find local maxima
  // std::vector<size_t> local_maxima;
  // for (const auto&& [node, node_prop, neighbors_view]: mg) {
  //   const auto [nbr_begin, nbr_end] = neighbors_view;
  //   bool local_max = true;
  //   for(auto nbr_it=nbr_begin; nbr_it!=nbr_end; ++nbr_it){
  //     if(mg.node_prop(*nbr_it).distance > node_prop.distance){
  //       local_max = false;
  //       break;
  //     }
  //   }
  //   if(local_max){
  //     local_maxima.push_back(node);
  //   }
  // }

  //Find maximum
  const size_t max_node = std::get<0>(*std::max_element(mg.begin(), mg.end(), [&](const auto &a, const auto &b){
    return std::get<1>(a).distance > std::get<1>(b).distance;
  }));

  std::unordered_set<size_t> visited;

  const auto cmp = [&](const size_t left, const size_t right) {
    return mg.node_prop(left).distance > mg.node_prop(right).distance;
  };
  std::priority_queue<size_t, std::vector<size_t>, decltype(cmp)> pq(cmp);
  pq.push(max_node);

  Polygon_set_2 gph;

  while(!pq.empty()){
    const auto c = pq.top();
    pq.pop();

    if(visited.count(c)!=0){
      continue;
    }
    visited.insert(c);

    const auto &props = mg.node_prop(c);

    if(props.distance==0){
      continue;
    }

    gph.join(construct_polygon_circle(props.pt, props.distance));

    const auto [nbr_begin, nbr_end] = mg.neighbors(c);
    for(auto n=nbr_begin;n!=nbr_end;++n){
      if(visited.count(*n)!=0){
        continue;
      }
      pq.push(*n);
    }
  }
}


int main(int argc, char** argv) {
  if(argc!=2){
      std::cerr<<"Syntax: "<<argv[0]<<" <Shape Boundary WKT>"<<std::endl;
      return -1;
  }

  CGAL::set_pretty_mode(std::cout);

  const auto mp = get_wkt_from_file(argv[1]);
  const auto voronoi = convert_mp_to_voronoi_diagram(mp);
  const auto ma_data = filter_voronoi_diagram_to_medial_axis(voronoi, mp);
  const auto medial_graph = get_medial_axis_graph(ma_data);

  iterate_graph(medial_graph);

  print_medial_axis_points(ma_data, "voronoi_points.csv");
  print_medial_axis_edges(ma_data, "voronoi_edges.csv");
  print_medial_graph_points_and_distances(medial_graph, "medial_graph_points.csv");

  // Print out the sources on which those edges rely. A Voronoi edge can be
  // formed by constraints imposed by two edges, two points, or a point and an
  // edge
  {
    std::ofstream fout("voronoi_sources.txt");
    for (std::size_t i = 0; i < ma_data.ups.size(); i++) {
      const auto& s1 = ma_data.ups[i]->site();
      if (s1.is_point()) {
        fout << "Point(" << s1.point().x() << "," << s1.point().y() << ") ";
      } else {
        fout << "Line("
          << s1.segment().source().x() << ","
          << s1.segment().source().y() << ","
          << s1.segment().target().x() << ","
          << s1.segment().target().y() << ") ";
      }

      const auto& s2 = ma_data.downs[i]->site();
      if (s2.is_point()) {
        fout << "Point(" << s2.point().x() << "," << s2.point().y() << ")";
      } else {
        fout << "Line("
          << s2.segment().source().x() << ","
          << s2.segment().source().y() << ","
          << s2.segment().target().x() << ","
          << s2.segment().target().y() << ")";
      }

      fout << std::endl;
    }
  }

  return 0;
}
