#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/WKT.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Segment_Delaunay_graph_2.h>
#include <CGAL/Segment_Delaunay_graph_adaptation_policies_2.h>
#include <CGAL/Segment_Delaunay_graph_adaptation_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_traits_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Voronoi_diagram_2.h>

#include <algorithm>
#include <deque>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <unordered_set>

typedef CGAL::Exact_predicates_inexact_constructions_kernel            K;
typedef CGAL::Segment_Delaunay_graph_traits_2<K>                       Gt;
typedef CGAL::Segment_Delaunay_graph_2<Gt>                             SDG2;
typedef CGAL::Segment_Delaunay_graph_adaptation_traits_2<SDG2>         AT;
typedef CGAL::Segment_Delaunay_graph_degeneracy_removal_policy_2<SDG2> AP;
typedef CGAL::Voronoi_diagram_2<SDG2, AT, AP>      VoronoiDiagram;
typedef AT::Site_2                                 Site_2;
typedef AT::Point_2                                Point_2;
typedef VoronoiDiagram::Locate_result              Locate_result;
typedef VoronoiDiagram::Vertex_handle              Vertex_handle;
typedef VoronoiDiagram::Face_handle                Face_handle;
typedef VoronoiDiagram::Halfedge_handle            Halfedge_handle;
typedef VoronoiDiagram::Ccb_halfedge_circulator    Ccb_halfedge_circulator;
typedef VoronoiDiagram::Bounded_halfedges_iterator BHE_Iter;
typedef VoronoiDiagram::Halfedge                   Halfedge;
typedef VoronoiDiagram::Vertex                     Vertex;
typedef CGAL::Polygon_with_holes_2<K> Polygon;
typedef std::deque<Polygon> MultiPolygon;

struct Point2Hash {
  size_t operator()(const Point_2 &pt) const {
    std::hash<double> hasher;
    auto seed = hasher(pt.x());
    seed ^= hasher(pt.y()) + 0x9e3779b9 + (seed<<6) + (seed>>2);
    return seed;
  }
};

typedef std::unordered_set<Point_2, Point2Hash> Point2Set;


/// Holds a more accessible description of the Voronoi diagram
struct VoronoiData {
  /// Map of vertices comprising the Voronoi diagram
  std::map<Vertex_handle, int> vertex_handles;
  /// List of edges in the diagram (pairs of the vertices above)
  std::vector<std::pair<int, int>> edges;
  /// Medial axis up governor. 1:1 correspondance with edges above.
  std::vector<VoronoiDiagram::Delaunay_graph::Vertex_handle> ups;
  /// Medial axis down governor. 1:1 correspondance with edges above.
  std::vector<VoronoiDiagram::Delaunay_graph::Vertex_handle> downs;
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


/// Converts a list of points (should be a closed loop) into a Voronoi diagram
VoronoiDiagram convert_mp_to_voronoi_diagram(const MultiPolygon &mp){
  VoronoiDiagram vd;

  const auto add_segments_to_vd = [&](const auto &poly){
    for(std::size_t i=0;i<poly.size();i++){
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


/// Find an `item` in `v` or add it if not present.
/// Returns the index of `item`'s location
template<class T, class U>
int find_or_add(std::map<T, int> &c, const U& item){
  // Map means we can do this in log(N) time
  if(c.count(item) == 0){
    c.emplace(item, c.size());
    return c.size() - 1;
  }

  return c.at(item);
}


/// Convert a map of <T, int> pairs to a vector of `T` ordered by increasing int
template<class T>
std::vector<T> map_to_ordered_vector(const std::map<T, int> &m){
  std::vector<std::pair<T, int>> to_sort(m.begin(), m.end());
  std::sort(to_sort.begin(), to_sort.end(), [](const auto &a, const auto &b){
    return a.second < b.second;
  });

  std::vector<T> ret;
  std::transform(begin(to_sort), end(to_sort), std::back_inserter(ret),
    [](auto const& pair){ return pair.first; }
  );

  return ret;
}

std::set<Vertex_handle> identify_vertex_handles_inside_mp(const VoronoiDiagram &vd, const MultiPolygon &mp){
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
      if(considered.count(vh)==0){
        considered.insert(vh);
        const auto inside_of_a_poly = std::any_of(
          mp.begin(), mp.end(), [&](const auto &poly) {
            // return inside(poly.oriented_side(vh->point()));
            return inside(CGAL::oriented_side(vh->point(), poly));
          }
        );
        if(inside_of_a_poly){
          interior.insert(vh);
        }
      }
    };

    vertex_in_mp_interior(edge_iter->source());
    vertex_in_mp_interior(edge_iter->target());
  }

  return interior;
}

Point2Set identify_concave_points_of_mp(const MultiPolygon &mp){
  Point2Set concave_points;

  // Determine cross product, given three points
  const auto z_cross_product = [](const Point_2 &pt1, const Point_2 &pt2, const Point_2 &pt3){
    const auto dx1 = pt2.x() - pt1.x();
    const auto dy1 = pt2.y() - pt1.y();
    const auto dx2 = pt3.x() - pt2.x();
    const auto dy2 = pt3.y() - pt2.y();
    return dx1 * dy2 - dy1 * dx2;
  };

  // Loop through all the points in a polygon and get their cross products
  // Sense should be `1` for outer boundaries and `-1` for holes (since holes)
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

VoronoiData get_voronoi_data_filtered_to_mp(
  const VoronoiDiagram &vd,
  const MultiPolygon &mp
){
  VoronoiData ret;

  const auto interior = identify_vertex_handles_inside_mp(vd, mp);
  const auto concave_points = identify_concave_points_of_mp(mp);

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

    if(pconcave(v1p->point()) || pconcave(v2p->point())){
      continue;
    }

    const auto id1 = find_or_add(ret.vertex_handles, v1p);
    const auto id2 = find_or_add(ret.vertex_handles, v2p);

    ret.edges.emplace_back(id1, id2);

    // Keep track of the medial axis governors
    ret.ups.push_back(halfedge.up());
    ret.downs.push_back(halfedge.down());
  }

  return ret;
}

/*
void reduce_voronoi_data_to_medial_axis_of_mp(const MultiPolygon &mp, VoronoiData &vd){
  // Now, we check each point to see if it falls inside or outside the
  // MultiPolygon
}*/


int main(int argc, char** argv) {
  if(argc!=2){
      std::cerr<<"Syntax: "<<argv[0]<<" <Shape Boundary>"<<std::endl;
      return -1;
  }

  CGAL::set_pretty_mode(std::cout);

  const auto mp = get_wkt_from_file(argv[1]);
  const auto voronoi = convert_mp_to_voronoi_diagram(mp);
  const auto vdata = get_voronoi_data_filtered_to_mp(voronoi, mp);

  // Print the points which collectively comprise the Voronoi diagram
  {
    std::ofstream fout("voronoi_points.csv");
    fout<<"x,y"<<std::endl;
    for (const auto &vh : map_to_ordered_vector(vdata.vertex_handles)) {
      fout << vh->point().x() << "," << vh->point().y() << std::endl;
    }
  }

  // Print out the edges of the Voronoi diagram
  {
    std::ofstream fout("voronoi_edges.csv");
    fout<<"SourceIdx,TargetIdx,UpGovernorIsPoint,DownGovernorIsPoint"<<std::endl;
    for (std::size_t i = 0; i < vdata.edges.size(); i++) {
      fout << vdata.edges[i].first        << ","
            << vdata.edges[i].second      << ","
            << vdata.ups[i]->is_point()   << "," // Is up-governor a point?
            << vdata.downs[i]->is_point()        // Is down-governor a point?
            << std::endl;
    }
  }

  // Print out the sources on which those edges rely. A Voronoi edge can be
  // formed by constraints imposed by two edges, two points, or a point and an
  // edge
  {
    std::ofstream fout("voronoi_sources.txt");
    for (std::size_t i = 0; i < vdata.ups.size(); i++) {
      const auto& s1 = vdata.ups[i]->site();
      if (s1.is_point()) {
        fout << "Point(" << s1.point().x() << "," << s1.point().y() << ") ";
      } else {
        fout << "Line("
          << s1.segment().source().x() << ","
          << s1.segment().source().y() << ","
          << s1.segment().target().x() << ","
          << s1.segment().target().y() << ") ";
      }

      const auto& s2 = vdata.downs[i]->site();
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
