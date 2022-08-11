// Compile with: clang++ -DBOOST_ALL_NO_LIB -DCGAL_USE_GMPXX=1 -O2 -g -DNDEBUG -Wall -Wextra -pedantic -march=native -frounding-math main.cpp -lgmpxx -lmpfr -lgmp
#include "geometry.hpp"

#include <graph_lite.h>

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
  K::FT squared_distance;
};


typedef graph_lite::Graph<size_t, MedialPoint> MedialGraph;


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

    // Avoid adding undirected edge twice
    if(id1>id2){
      continue;
    }

    ret.edges.emplace_back(id1, id2);

    // Keep track of the medial axis governors
    ret.ups.push_back(halfedge.up());
    ret.downs.push_back(halfedge.down());
  }

  ret.vertex_handles = map_to_ordered_vector(handles);

  return ret;
}


/// Given 2 sites returns the minimum distance between the point and both of them
auto point_to_min_site_distance(const Point_2 &pt, const Site_2 &site1, const Site_2 &site2){
  std::cout<<"dist "<<point_to_site_distance(pt, site1)<<" "<<point_to_site_distance(pt, site2)<<std::endl;
  return std::min(point_to_site_distance(pt, site1), point_to_site_distance(pt, site2));
}


MedialGraph get_medial_axis_graph(const MedialData &md){
  MedialGraph mg;

  size_t vertex_id = 0;

  for(const auto &vh: md.vertex_handles){
    mg.add_node_with_prop(vertex_id, MedialPoint{
      .pt = vh->point(),
      .vd_vertex = true,
      .start_handle = vertex_id,
      .end_handle = vertex_id,
      .squared_distance = -1  // We'll fix this below
    });
    vertex_id++;
  }

  // Fix distances to vertex handles
  for(size_t e=0;e<md.edges.size();e++){
    const auto &edge = md.edges.at(e);
    auto &spt = mg.node_prop(mg.find(edge.first));
    spt.squared_distance = point_to_min_site_distance(
      spt.pt, md.ups.at(e)->site(), md.downs.at(e)->site()
    );
    auto &ept = mg.node_prop(mg.find(edge.second));
    ept.squared_distance = point_to_min_site_distance(
      ept.pt, md.ups.at(e)->site(), md.downs.at(e)->site()
    );
    // mg.add_edge(edge.first, edge.second); //TODO
  }

  // NOTE: No squared_distance should be -1 any more

  // Get distances to all other points and build a graph from there
  constexpr auto STEPS = 20;
  for(size_t e=0;e<md.edges.size();e++){
    const auto &edge = md.edges.at(e);
    const auto &svh = md.vertex_handles.at(edge.first);  // Start vertex
    const auto &evh = md.vertex_handles.at(edge.second); // End vertex
    InterPointInterpolator interp(svh->point(), evh->point());
    for(int t=1;t<STEPS;t++){
      // Location of the interpolated point
      const auto interp_point = interp(t/static_cast<double>(STEPS));

      mg.add_node_with_prop(vertex_id, MedialPoint{
        .pt = interp_point,
        .vd_vertex = false,
        .start_handle = edge.first,
        .end_handle = edge.second,
        .squared_distance = point_to_min_site_distance(
          interp_point, md.ups.at(e)->site(), md.downs.at(e)->site()
        )
      });

      if(t==1){
        mg.add_edge(edge.first, vertex_id);
        std::cout<<"graph "<<edge.first<<"-"<<vertex_id<<std::endl;
      } else if(t==STEPS-1){
        mg.add_edge(vertex_id - 1, vertex_id);
        mg.add_edge(vertex_id, edge.second);
        std::cout<<"graph "<<(vertex_id - 1)<<"-"<<vertex_id<<std::endl;
        std::cout<<"graph "<<vertex_id<<"-"<<edge.second<<std::endl;
      } else {
        mg.add_edge(vertex_id - 1, vertex_id);
        std::cout<<"graph "<<(vertex_id - 1)<<"-"<<vertex_id<<std::endl;
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
  for (const auto& node: mg) {
    const auto props = mg.node_prop(node);
    fout<<props.pt.x()<<","<<props.pt.y()<<","<<props.squared_distance<<std::endl;
  }

  std::ofstream fg("graph_connections");
  for (const auto& node: mg) {
    const auto [nbegin, nend] = mg.neighbors(node);
    for(auto n=nbegin;n!=nend;++n){
      fg<<node<<"-"<<(*n)<<std::endl;
    }
  }
}

/*
void draw_arc(const Circle_2 &circle, const Point_2 &start, const Point_2 &end, std::ostream &fout){
  constexpr auto STEPS = 10;

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
*/

void iterate_graph(const MedialGraph &mg){
  // Set up a max-heap priority-queue
  const auto cmp = [&](const size_t left, const size_t right) {
    return mg.node_prop(left).squared_distance > mg.node_prop(right).squared_distance;
  };
  std::priority_queue<size_t, std::vector<size_t>, decltype(cmp)> pq(cmp);

  // Convenience function for getting the square distance of a medial graph node
  const auto sdist = [&](const auto &node){
    return mg.node_prop(node).squared_distance;
  };

  // Find maxima
  //TODO: Identifies a single maximum distance point, which is maybe a bad
  //      assumption!
  {
    std::vector<size_t> max_distance_nodes;
    for(const auto &node: mg){
      if(max_distance_nodes.empty() || sdist(max_distance_nodes.front()) == sdist(node)){
        max_distance_nodes.push_back(node);
      } else if (sdist(node) > sdist(max_distance_nodes.front())) {
        max_distance_nodes.clear();         // Clear the priority queue
        max_distance_nodes.push_back(node); // Add new most distant node
      }
    }

    if(pq.size()>1){
      std::cerr<<"More than one starting node had equal distance!"<<std::endl;
    }

    for(const auto &x: max_distance_nodes){
      pq.push(x);
    }
  }

  std::cerr<<"Maxima found = " << pq.size() << std::endl;

  std::unordered_set<size_t> visited;

  Polygon_set_2 gph;

  int step=0;
  while(!pq.empty()){
    // Visit the next node
    const auto c = pq.top();
    pq.pop();

    // Ensure we don't visit a node twice
    if(visited.count(c)!=0){
      continue;
    }

    // Make a note that we've visited this node
    visited.insert(c);

    // Get node properties
    const auto &props = mg.node_prop(c);

    std::cout<<"Visiting "<<props.pt<<" with sradius "<<props.squared_distance<<std::endl;

    // If the node is just a point then we ignore it.
    // TODO: Maybe we should look at the neighbours anyway.
    if(props.squared_distance==0){
      continue;
    }

    // Add this node's supporting circle to the figure we are building
    gph.join(construct_polygon(props.pt, props.squared_distance));

    print_polygon(gph, "polygon_w_holes_output_step" + std::to_string(step));

    // Add neighbours of the point to the queue
    const auto [nbr_begin, nbr_end] = mg.neighbors(c);
    for(auto n=nbr_begin;n!=nbr_end;++n){
      if(visited.contains(*n)){
        continue;
      }
      pq.push(*n);
    }

    step++;
  }

  std::cout<<"polygons = "<<gph.number_of_polygons_with_holes()<<std::endl;
  std::cout<<"area = "<<area(gph)<<std::endl;
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

  print_medial_axis_points(ma_data, "voronoi_points.csv");
  print_medial_axis_edges(ma_data, "voronoi_edges.csv");
  print_medial_graph_points_and_distances(medial_graph, "medial_graph_points.csv");

  iterate_graph(medial_graph);

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
