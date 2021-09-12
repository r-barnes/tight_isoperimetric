#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Segment_Delaunay_graph_2.h>
#include <CGAL/Segment_Delaunay_graph_adaptation_policies_2.h>
#include <CGAL/Segment_Delaunay_graph_adaptation_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_traits_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Voronoi_diagram_2.h>

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>

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
typedef CGAL::Polygon_2<K>                         Polygon_2;


/// Holds a more accessible description of the Voronoi diagram
struct VoronoiData {
  /// List of vertices comprising the Voronoi diagram
  std::vector<Vertex_handle> vertex_handles;
  /// List of edges in the diagram (pairs of the vertices above)
  std::vector<std::pair<int, int>> edges;
  /// Medial axis up governor. 1:1 correspondance with edges above.
  std::vector<VoronoiDiagram::Delaunay_graph::Vertex_handle> ups;
  /// Medial axis down governor. 1:1 correspondance with edges above.
  std::vector<VoronoiDiagram::Delaunay_graph::Vertex_handle> downs;
};


/// Read @p filename to obtain shape boundary
std::vector<Point_2> get_boundary_points_from_file(const std::string& filename){
  std::vector<Point_2> points;
  {
    std::ifstream fp(filename);
    if(!fp.good()){
      throw std::runtime_error("Couldn't open file '" + filename + "'!");
    }

    double a;
    double b;
    while (fp >> a >> b) {
      points.emplace_back(a,b);
    }
  }

  if(points.empty()){
    throw std::runtime_error("Points file '" + filename + "' was empty!");
  }

  // Make sure the points describe a closed loop
  if(points.front()!=points.back()){
    points.push_back(points.front());
  }

  return points;
}


/// Converts a list of points (should be a closed loop) into a Voronoi diagram
VoronoiDiagram convert_point_list_to_voronoi_diagram(const std::vector<Point_2> &points){
  VoronoiDiagram vd;

  // Define sites
  for (std::size_t i = 0; i<points.size()-1; i++) {
    vd.insert(Site_2::construct_site_2(points[i], points[i+1]));
  }
  if(!vd.is_valid()){
    throw std::runtime_error("Voronoi Diagram was not valid!");
  }

  return vd;
}


/// Find an `item` in `v` or add it if not present.
/// Returns the index of `item`'s location
template<class T, class U>
int find_or_add(std::vector<T> &v, const U& item){
  auto idx = std::find(v.begin(), v.end(), item);
  if(idx == v.end()){
    v.push_back(item);
    idx = v.end() - 1;
  }

  return std::distance(v.begin(), idx);
}


VoronoiData get_voronoi_data(const VoronoiDiagram &vd){
  VoronoiData ret;

  // The Voronoi diagram is comprised of a number of vertices connected by lines
  // Here, we go through each edge of the Voronoi diagram and determine which
  // vertices it's incident on. We add these vertices to ret.vertex_handles
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

    const auto id1 = find_or_add(ret.vertex_handles, v1p);
    const auto id2 = find_or_add(ret.vertex_handles, v2p);

    ret.edges.emplace_back(id1, id2);

    // Keep track of the medial axis governors
    ret.ups.push_back(halfedge.up());
    ret.downs.push_back(halfedge.down());
  }

  return ret;
}


int main(int argc, char** argv) {
  if(argc!=2){
      std::cerr<<"Syntax: "<<argv[0]<<" <Shape Boundary>"<<std::endl;
      return -1;
  }

  CGAL::set_pretty_mode(std::cout);

  const auto points = get_boundary_points_from_file(argv[1]);
  const auto voronoi = convert_point_list_to_voronoi_diagram(points);
  const auto vdata = get_voronoi_data(voronoi);

  // Print the points which collectively comprise the Voronoi diagram
  {
    std::ofstream fout("voronoi_points.txt");
    for (const auto &vhti : vdata.vertex_handles) {
      fout << vhti->point().x() << " " << vhti->point().y() << std::endl;
    }
  }

  // Print out the edges of the Voronoi diagram
  {
    std::ofstream fout("voronoi_edges.txt");
    for (std::size_t i = 0; i < vdata.edges.size(); i++) {
      fout << vdata.edges[i].first       << " "
            << vdata.edges[i].second      << " "
            << vdata.ups[i]->is_point()   << " "
            << vdata.downs[i]->is_point() << std::endl;
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
