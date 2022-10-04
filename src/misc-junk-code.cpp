// TODO(r-barnes): Remove this file eventually

// Find local maxima
//  std::vector<size_t> local_maxima;
//  for (const auto&& [node, node_prop, neighbors_view]: mg) {
//    const auto [nbr_begin, nbr_end] = neighbors_view;
//    bool local_max = true;
//    for(auto nbr_it=nbr_begin; nbr_it!=nbr_end; ++nbr_it){
//      if(mg.node_prop(*nbr_it).distance > node_prop.distance){
//        local_max = false;
//        break;
//      }
//    }
//    if(local_max){
//      local_maxima.push_back(node);
//    }
//  }

void iterate_graph(const MedialGraph& mg) {
  // Find local maxima
  //  std::vector<size_t> local_maxima;
  //  for (const auto&& [node, node_prop, neighbors_view]: mg) {
  //    const auto [nbr_begin, nbr_end] = neighbors_view;
  //    bool local_max = true;
  //    for(auto nbr_it=nbr_begin; nbr_it!=nbr_end; ++nbr_it){
  //      if(mg.node_prop(*nbr_it).distance > node_prop.distance){
  //        local_max = false;
  //        break;
  //      }
  //    }
  //    if(local_max){
  //      local_maxima.push_back(node);
  //    }
  //  }

  // Find maximum
  const size_t max_node = std::get<0>(*std::max_element(mg.begin(), mg.end(), [&](const auto& a, const auto& b) {
    return std::get<1>(a).squared_distance < std::get<1>(b).squared_distance;
  }));

  std::unordered_set<size_t> visited;

  const auto cmp = [&](const size_t left, const size_t right) {
    return mg.node_prop(left).squared_distance > mg.node_prop(right).squared_distance;
  };
  std::priority_queue<size_t, std::vector<size_t>, decltype(cmp)> pq(cmp);
  pq.push(max_node);

  Polygon_set_2 gph;

  while (!pq.empty()) {
    const auto c = pq.top();
    pq.pop();

    if (visited.count(c) != 0) {
      continue;
    }
    visited.insert(c);

    const auto& props = mg.node_prop(c);

    if (props.squared_distance == 0) {
      continue;
    }

    gph.join(construct_polygon_circle(props.pt, props.squared_distance));

    const auto [nbr_begin, nbr_end] = mg.neighbors(c);
    for (auto n = nbr_begin; n != nbr_end; ++n) {
      if (visited.count(*n) != 0) {
        continue;
      }
      pq.push(*n);
    }
  }
  std::cout << std::endl;

  std::cout << "HIIIIIIIIIIIIIIIIII!" << std::endl;
  std::list<Polygon_with_holes_2> res;
  gph.polygons_with_holes(std::back_inserter(res));
  std::cout << "pwh size = " << res.size() << std::endl;
  std::copy(res.begin(), res.end(), std::ostream_iterator<Polygon_with_holes_2>(std::cout, "\n"));
  std::cout << std::endl;
  for (const auto& r : res) {
    std::cout << "r area = " << r.outer_boundary().area() << std::endl;
    print_polygon_with_holes(r, "poly_holes.csv");
  }
}

/// Constructs a polygon from a circle
Polygon_2 construct_polygon_circle(const Point_2& pt, const K::FT squared_radius) {
  // Subdivide the circle into two x-monotone arcs.
  Traits_2 traits;
  Curve_2 curve(Circle_2(pt, squared_radius));

  // This output is correct (TODO)
  // std::cout<<"circ "<<pt.x()<<","<<pt.y()<<","<<std::sqrt(CGAL::to_double(squared_radius))<<std::endl;

  std::list<CGAL::Object> objects;
  traits.make_x_monotone_2_object()(curve, std::back_inserter(objects));
  CGAL_assertion(objects.size() == 2);

  // Construct the polygon.
  Polygon_2 pgn;
  X_monotone_curve_2 arc;
  std::list<CGAL::Object>::iterator iter;
  for (iter = objects.begin(); iter != objects.end(); ++iter) {
    CGAL::assign(arc, *iter);
    // TODO
    // const Point_2 S = arc.source();
    // const Point_2 T = arc.target();
    // const Point_2 S(CGAL::to_double(arc.source().x()),CGAL::to_double(arc.source().y()));
    // const Point_2 T(CGAL::to_double(arc.target().x()),CGAL::to_double(arc.target().y()));
    // std::cout<<"arc "<<T.x()<<","<<T.y()<<std::endl;
    // draw_arc(arc.supporting_circle(), S, T, std::cout);
    pgn.push_back(arc);
  }

  return pgn;
}

void draw_arc(const Circle_2& circle, const Point_2& start, const Point_2& end, std::ostream& fout) {
  constexpr auto STEPS = 10;

  const Vector_2 x_axis(1, 0);
  const auto start_angle = 0;         // CGAL::angle(start - circle.center(), x_axis);
  const auto ang_diff    = 2 * M_PI;  // CGAL::angle(start, circle.center(), end);
  const auto cx          = circle.center().x();
  const auto cy          = circle.center().y();
  const auto radius      = std::sqrt(CGAL::to_double(circle.squared_radius()));
  // auto px = start.x();
  // auto py = start.y();
  for (int i = 0; i < STEPS; i++) {
    const auto t = i / static_cast<double>(STEPS);
    const auto x = cx + radius * std::cos(start_angle + ang_diff * t);
    const auto y = cy + radius * std::sin(start_angle + ang_diff * t);
    // fout<<px<<","<<py<<","<<x<<","<<y<<std::endl;
    fout << x << "," << y << std::endl;
    // px = x;
    // py = y;
  }
}

#mgp = np.loadtxt("/z/circs.csv", dtype = float, skiprows = 1, delimiter = ',')
#ax.scatter(mgp[:, 0], mgp[:, 1]) #, c = mgp[ :, 2], cmap = 'Greens')
#for i in range(len(mgp[:, 0])):
#ax.text(mgp[i, 0], mgp[i, 1], s = "{0:.2f}".format(mgp[i, 2]))
#circle = plt.Circle((mgp[i, 0], mgp[i, 1]), mgp[i, 2])
#ax.add_patch(circle)

#union_lines = np.loadtxt("/z/arcs.csv", dtype = float, delimiter = ',')
#ax.scatter(union_lines[:, 0], union_lines[:, 1])

#union_lines = np.loadtxt("poly_holes.csv", dtype = float, delimiter = ',')
#ax.scatter(union_lines[:, 0], union_lines[:, 1])

#ea = np.loadtxt("v_edgeattributes.txt")
#for eal in ea:
#ax.plot((eal[0], eal[2]), (eal[1], eal[3]))
