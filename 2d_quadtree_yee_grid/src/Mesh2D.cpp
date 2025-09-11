#include <vtkWriter.h>
#include <Mesh2D.h>
#include <TopoDS.hxx>           
#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <TopExp_Explorer.hxx>
#include <BRepBndLib.hxx>
#include <Bnd_Box.hxx>
#include <STEPControl_Reader.hxx>
#include <BRepTools.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <TopAbs_State.hxx>
#include <gp_Pnt.hxx>
#include <vector>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <array>

Mesh2D::Mesh2D(const std::string& stepFile, double cellSize)
  : stepFile_(stepFile), cellSize_(cellSize) {}

// ------------------- Read STEP file ------------------- //
bool Mesh2D::readSTEP() {
  STEPControl_Reader reader;
  if (reader.ReadFile(stepFile_.c_str()) != IFSelect_RetDone) {
    std::cerr << "Cannot read STEP file: " << stepFile_ << std::endl;
    return false;
  }

  reader.TransferRoots();
  TopoDS_Shape shape = reader.OneShape();

  // Assume single face (2D sheet)
  for (TopExp_Explorer faceExp(shape, TopAbs_FACE); faceExp.More(); faceExp.Next()) {
    face_ = TopoDS::Face(faceExp.Current());  // <<--- use TopoDS::Face
    return true;
  }


  std::cerr << "No face found in STEP file." << std::endl;
  return false;
}

// ------------------- Check if point is inside face ------------------- // 
bool Mesh2D::isPointInside(const gp_Pnt& p) const {
  BRepClass_FaceClassifier classifier(face_, p, 1e-7);
  TopAbs_State state = classifier.State();
  return (state == TopAbs_IN || state == TopAbs_ON);
}

// ------------------- Generate uniform grid -------------------
void Mesh2D::generateUniformGrid() {
  double xmin = -2.0, xmax = 2.0;
  double ymin = -2.0, ymax = 2.0;

  int id = 0;
  for (double x = xmin; x <= xmax; x += cellSize_) {
    for (double y = ymin; y <= ymax; y += cellSize_) {
      EzNode n_Ez(x, y, id);   
      HxNode n_Hx(x, y - cellSize_/2, id);
      HyNode n_Hy(x - cellSize_/2, y, id);
      Ez_nodes.push_back(n_Ez);
      Hx_nodes.push_back(n_Hx);
      Hy_nodes.push_back(n_Hy);

      id++;
    }
  }
}

// ------------------- Filter cells outside .STEP geometry -------------------//
void Mesh2D::filterNodesOutsideGeometry() {
  std::vector<EzNode> outside_Ez_Nodes;
  std::vector<HxNode> outside_Hx_Nodes;
  std::vector<HyNode> outside_Hy_Nodes;
  
  int id = 0;
  for (auto &n : Ez_nodes) {
    gp_Pnt p(n.x, n.y, 0);  

    if (isPointInside(p)) {
      n.nodeID = id;        
      n.nodePoint = p;      
      outside_Ez_Nodes.push_back(n);            
      id++;
    }
  }

  id = 0;
  for (auto &n : Hx_nodes) {
    gp_Pnt p(n.x, n.y, 0);  

    if (isPointInside(p)) {
      n.nodeID = id;        
      n.nodePoint = p;      
      outside_Hx_Nodes.push_back(n);            
      id++;
    }
  }

  id = 0;
  for (auto &n : Hy_nodes) {
    gp_Pnt p(n.x, n.y, 0);  

    if (isPointInside(p)) {
      n.nodeID = id;        
      n.nodePoint = p;      
      outside_Hy_Nodes.push_back(n);            
      id++;
    }
  }

  Ez_nodes = outside_Ez_Nodes;
  Hx_nodes = outside_Hx_Nodes;
  Hy_nodes = outside_Hy_Nodes;

}

template <typename NodeType>
void linkNeighbors(std::vector<NodeType> &nodes,
                   std::unordered_map<std::string, NodeType*> &lookup,
                   double cellSize) {
  auto makeKey = [](double x, double y) {
    return std::to_string(x) + "_" + std::to_string(y);
  };


  // build lookup
  for (auto &n : nodes) {
    lookup[makeKey(n.x, n.y)] = &n;
  }

  // assign neighbors
  for (auto &n : nodes) {
    n.left   = lookup.count(makeKey(n.x - cellSize, n.y)) ? lookup[makeKey(n.x - cellSize, n.y)] : nullptr;
    n.right  = lookup.count(makeKey(n.x + cellSize, n.y)) ? lookup[makeKey(n.x + cellSize, n.y)] : nullptr;
    n.top    = lookup.count(makeKey(n.x, n.y + cellSize)) ? lookup[makeKey(n.x, n.y + cellSize)] : nullptr;
    n.bottom = lookup.count(makeKey(n.x, n.y - cellSize)) ? lookup[makeKey(n.x, n.y - cellSize)] : nullptr;
  }
}

void Mesh2D::linkNodeNeighbors() {
  std::unordered_map<std::string, EzNode*> Ez_lookup;
  std::unordered_map<std::string, HxNode*> Hx_lookup;
  std::unordered_map<std::string, HyNode*> Hy_lookup;

  linkNeighbors(Ez_nodes, Ez_lookup, cellSize_);
  linkNeighbors(Hx_nodes, Hx_lookup, cellSize_);
  linkNeighbors(Hy_nodes, Hy_lookup, cellSize_);
}

void Mesh2D::linkCrossNeighbors() {
  auto makeKey = [](double x, double y) {
    return std::to_string(x) + "_" + std::to_string(y);
  };

  // --- lookup maps ---
  std::unordered_map<std::string, EzNode*> Ez_lookup;
  std::unordered_map<std::string, HxNode*> Hx_lookup;
  std::unordered_map<std::string, HyNode*> Hy_lookup;

  for (auto &n : Ez_nodes) Ez_lookup[makeKey(n.x, n.y)] = &n;
  for (auto &n : Hx_nodes) Hx_lookup[makeKey(n.x, n.y)] = &n;
  for (auto &n : Hy_nodes) Hy_lookup[makeKey(n.x, n.y)] = &n;

  // --- Ez -> Hx / Hy ---
  for (auto &n : Ez_nodes) {
    n.hx_top    = Hx_lookup.count(makeKey(n.x, n.y + cellSize_/2)) ? Hx_lookup[makeKey(n.x, n.y + cellSize_/2)] : nullptr;
    n.hx_bottom = Hx_lookup.count(makeKey(n.x, n.y - cellSize_/2)) ? Hx_lookup[makeKey(n.x, n.y - cellSize_/2)] : nullptr;

    n.hy_left   = Hy_lookup.count(makeKey(n.x - cellSize_/2, n.y)) ? Hy_lookup[makeKey(n.x - cellSize_/2, n.y)] : nullptr;
    n.hy_right  = Hy_lookup.count(makeKey(n.x + cellSize_/2, n.y)) ? Hy_lookup[makeKey(n.x + cellSize_/2, n.y)] : nullptr;
  }

  // --- Hx -> Ez ---
  for (auto &n : Hx_nodes) {
    n.ez_top    = Ez_lookup.count(makeKey(n.x, n.y + cellSize_/2)) ? Ez_lookup[makeKey(n.x, n.y + cellSize_/2)] : nullptr;
    n.ez_bottom = Ez_lookup.count(makeKey(n.x, n.y - cellSize_/2)) ? Ez_lookup[makeKey(n.x, n.y - cellSize_/2)] : nullptr;
  }

  // --- Hy -> Ez ---
  for (auto &n : Hy_nodes) {
    n.ez_left   = Ez_lookup.count(makeKey(n.x - cellSize_/2, n.y)) ? Ez_lookup[makeKey(n.x - cellSize_/2, n.y)] : nullptr;
    n.ez_right  = Ez_lookup.count(makeKey(n.x + cellSize_/2, n.y)) ? Ez_lookup[makeKey(n.x + cellSize_/2, n.y)] : nullptr;
  }
}

void Mesh2D::add_ghost_layer() {
  Ez_domain_size = Ez_nodes.size();
  Hx_domain_size = Hx_nodes.size();
  Hy_domain_size = Hy_nodes.size();

  std::cout<<"---------------------------------"<<std::endl;
  std::cout<<"Ez_domain_size = "<<Ez_domain_size<<std::endl;
  std::cout<<"Hx_domain_size = "<<Hx_domain_size<<std::endl;
  std::cout<<"Hy_domain_size = "<<Hy_domain_size<<std::endl;
  std::cout<<"---------------------------------"<<std::endl;

  // Make stable copies of domain nodes (to avoid mixing with newly added ghosts)
  std::vector<HxNode> Hx_copy_nodes = Hx_nodes;
  std::vector<HyNode> Hy_copy_nodes = Hy_nodes;

  int hx_id = 0, hy_id = 0, ez_id = 0;

  // --- Ez nodes: create missing Hx/Hy ghosts ---
  for (auto &n : Ez_nodes) {
    if (n.hx_bottom == nullptr) {
      Hx_nodes.emplace_back(n.x, n.y - cellSize_/2, Hx_domain_size + hx_id);
      n.hx_bottom = &Hx_nodes.back();
      hx_id++;
    }
    if (n.hx_top == nullptr) {
      Hx_nodes.emplace_back(n.x, n.y + cellSize_/2, Hx_domain_size + hx_id);
      n.hx_top = &Hx_nodes.back();
      hx_id++;
    }
    if (n.hy_left == nullptr) {
      Hy_nodes.emplace_back(n.x - cellSize_/2, n.y, Hy_domain_size + hy_id);
      n.hy_left = &Hy_nodes.back();
      hy_id++;
    }
    if (n.hy_right == nullptr) {
      Hy_nodes.emplace_back(n.x + cellSize_/2, n.y, Hy_domain_size + hy_id);
      n.hy_right = &Hy_nodes.back();
      hy_id++;
    }
  }

  // --- Hx nodes: create missing Ez ghosts ---
  for (int i = 0; i < Hx_domain_size; i++) {
    const auto &orig = Hx_copy_nodes[i];
    if (orig.ez_bottom == nullptr) {
      Ez_nodes.emplace_back(orig.x, orig.y - cellSize_/2, Ez_domain_size + ez_id);
      Hx_nodes[i].ez_bottom = &Ez_nodes.back();
      ez_id++;
    }
    if (orig.ez_top == nullptr) {
      Ez_nodes.emplace_back(orig.x, orig.y + cellSize_/2, Ez_domain_size + ez_id);
      Hx_nodes[i].ez_top = &Ez_nodes.back();
      ez_id++;
    }
  }

  // --- Hy nodes: create missing Ez ghosts ---
  for (int i = 0; i < Hy_domain_size; i++) {
    const auto &orig = Hy_copy_nodes[i];
    if (orig.ez_left == nullptr) {
      Ez_nodes.emplace_back(orig.x - cellSize_/2, orig.y, Ez_domain_size + ez_id);
      Hy_nodes[i].ez_left = &Ez_nodes.back();
      ez_id++;
    }
    if (orig.ez_right == nullptr) {
      Ez_nodes.emplace_back(orig.x + cellSize_/2, orig.y, Ez_domain_size + ez_id);
      Hy_nodes[i].ez_right = &Ez_nodes.back();
      ez_id++;
    }
  }

  // record ghost counts
  Ez_ghost_cells = Ez_nodes.size() - Ez_domain_size;
  Hx_ghost_cells = Hx_nodes.size() - Hx_domain_size;
  Hy_ghost_cells = Hy_nodes.size() - Hy_domain_size;
}

void Mesh2D::check_nullptr() {
  for (int i = 0; i < Ez_domain_size; i++) {
    if (Ez_nodes[i].hx_bottom == nullptr) std::cout << "Ez["<<i<<"].hx_bottom == nullptr\n";
    if (Ez_nodes[i].hx_top    == nullptr) std::cout << "Ez["<<i<<"].hx_top == nullptr\n";
    if (Ez_nodes[i].hy_left   == nullptr) std::cout << "Ez["<<i<<"].hy_left == nullptr\n";
    if (Ez_nodes[i].hy_right  == nullptr) std::cout << "Ez["<<i<<"].hy_right == nullptr\n";
  
  }

  for (int i = 0; i < Hx_domain_size; i++) {  
    if (Hx_nodes[i].ez_top    == nullptr) std::cout << "Hx["<<i<<"].ez_top == nullptr\n";
    if (Hx_nodes[i].ez_bottom == nullptr) std::cout << "Hx["<<i<<"].ez_bottom == nullptr\n";
  }

  for (int i = 0; i < Hy_domain_size; i++) {    
    if (Hy_nodes[i].ez_left  == nullptr) std::cout << "Hy["<<i<<"].ez_left == nullptr\n";
    if (Hy_nodes[i].ez_right == nullptr) std::cout << "Hy["<<i<<"].ez_right == nullptr\n";
  }
}

void Mesh2D::check_connects() {

  std::cout<<"Ez_nodes[100].x = "<<Ez_nodes[100].x<<"\n";
  std::cout<<"Ez_nodes[100].y = "<<Ez_nodes[100].y<<"\n";

  std::cout<<"Ez_nodes[100].hx_top->x = "<<Ez_nodes[100].hx_top->x<<"\n";
  std::cout<<"Ez_nodes[100].hx_top->y = "<<Ez_nodes[100].hx_top->y<<"\n";
  std::cout<<"Ez_nodes[100].hx_top->y = "<<Ez_nodes[100].hx_top->nodeID<<"\n";

  std::cout<<"Ez_nodes[100].hx_top->x = "<<Ez_nodes[100].hx_bottom->x<<"\n";
  std::cout<<"Ez_nodes[100].hx_top->y = "<<Ez_nodes[100].hx_bottom->y<<"\n";
  std::cout<<"Ez_nodes[100].hx_top->nodeID = "<<Ez_nodes[100].hx_bottom->nodeID<<"\n";
  std::cout<<"Ez_nodes[100].hx_bottom->ez_top->nodeID = "<<Ez_nodes[100].hx_bottom->ez_top->nodeID<<"\n";


}

void Mesh2D::saveMeshToVTK(const std::string& type, const std::string& filename) const {
  VTKQuadWriter vtk2elements;
  std::vector<double> coords;
  std::vector<int> connectivity;
  std::vector<double> field;
  std::string fieldName;

  if (type == "Ez") {
    fieldName = "Ez";
    // ---- Ez export ----
    for (int i = 0; i < Ez_domain_size; i++) {
      const auto &n = Ez_nodes[i];
      coords.push_back(n.x);
      coords.push_back(n.y);
      coords.push_back(0.0);
      
    }

    for (int i = 0; i < Ez_domain_size; i++) {
      const auto &n = Ez_nodes[i];
      if (n.right && n.top && n.top->right) {
        connectivity.push_back(n.nodeID);
        connectivity.push_back(n.right->nodeID);
        connectivity.push_back(n.top->right->nodeID);
        connectivity.push_back(n.top->nodeID);
        field.push_back(n.fieldValue);
      }
    }
  }
  else if (type == "Hx") {
    fieldName = "Hx";
    // ---- Hx export ----
    for (int i = 0; i < Hx_domain_size; i++) {
      const auto &n = Hx_nodes[i];
      coords.push_back(n.x);
      coords.push_back(n.y);
      coords.push_back(0.0);
      
    }

    for (int i = 0; i < Hx_domain_size; i++) {
      const auto &n = Hx_nodes[i];
      if (n.right && n.top && n.top->right) {
        connectivity.push_back(n.nodeID);
        connectivity.push_back(n.right->nodeID);
        connectivity.push_back(n.top->right->nodeID);
        connectivity.push_back(n.top->nodeID);
        field.push_back(n.fieldValue);
      }
    }
    // often Hx are exported as points, so you can drop the connectivity part if not needed
  }
  else if (type == "Hy") {
    fieldName = "Hy";
    // ---- Hy export ----
    for (int i = 0; i < Hy_domain_size; i++) {
      const auto &n = Hy_nodes[i];
      coords.push_back(n.x);
      coords.push_back(n.y);
      coords.push_back(0.0);
      
    }

    for (int i = 0; i < Hy_domain_size; i++) {
      const auto &n = Hy_nodes[i];
      if (n.right && n.top && n.top->right) {
        connectivity.push_back(n.nodeID);
        connectivity.push_back(n.right->nodeID);
        connectivity.push_back(n.top->right->nodeID);
        connectivity.push_back(n.top->nodeID);
        field.push_back(n.fieldValue);
      }
    }
    // same note as Hx
  }
  else {
    throw std::runtime_error("Unknown type: " + type);
  }

  vtk2elements.set_points(coords);
  vtk2elements.set_cells(connectivity);
  vtk2elements.add_scalar(field, fieldName);
  vtk2elements.write_vtk(filename);

}
