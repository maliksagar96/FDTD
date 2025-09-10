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

  static NodeType ghostField(0,0,0);
  ghostField.fieldValue = 0;

  // build lookup
  for (auto &n : nodes) {
    lookup[makeKey(n.x, n.y)] = &n;
  }

  // assign neighbors
  for (auto &n : nodes) {
    n.left   = lookup.count(makeKey(n.x - cellSize, n.y)) ? lookup[makeKey(n.x - cellSize, n.y)] : &ghostField;
    n.right  = lookup.count(makeKey(n.x + cellSize, n.y)) ? lookup[makeKey(n.x + cellSize, n.y)] : &ghostField;
    n.top    = lookup.count(makeKey(n.x, n.y + cellSize)) ? lookup[makeKey(n.x, n.y + cellSize)] : &ghostField;
    n.bottom = lookup.count(makeKey(n.x, n.y - cellSize)) ? lookup[makeKey(n.x, n.y - cellSize)] : &ghostField;
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

  static EzNode ghostEz(0,0,0);
  static HxNode ghostHx(0,0,0);
  static HyNode ghostHy(0,0,0);

  ghostEz.fieldValue = 0;
  ghostHx.fieldValue = 0;
  ghostHy.fieldValue = 0;

  // --- lookup maps ---
  std::unordered_map<std::string, EzNode*> Ez_lookup;
  std::unordered_map<std::string, HxNode*> Hx_lookup;
  std::unordered_map<std::string, HyNode*> Hy_lookup;

  for (auto &n : Ez_nodes) Ez_lookup[makeKey(n.x, n.y)] = &n;
  for (auto &n : Hx_nodes) Hx_lookup[makeKey(n.x, n.y)] = &n;
  for (auto &n : Hy_nodes) Hy_lookup[makeKey(n.x, n.y)] = &n;

  // --- Ez -> Hx / Hy ---
  for (auto &n : Ez_nodes) {
    n.hx_top    = Hx_lookup.count(makeKey(n.x, n.y + cellSize_/2)) ? Hx_lookup[makeKey(n.x, n.y + cellSize_/2)] : &ghostHx;
    n.hx_bottom = Hx_lookup.count(makeKey(n.x, n.y - cellSize_/2)) ? Hx_lookup[makeKey(n.x, n.y - cellSize_/2)] : &ghostHx;

    n.hy_left   = Hy_lookup.count(makeKey(n.x - cellSize_/2, n.y)) ? Hy_lookup[makeKey(n.x - cellSize_/2, n.y)] : &ghostHy;
    n.hy_right  = Hy_lookup.count(makeKey(n.x + cellSize_/2, n.y)) ? Hy_lookup[makeKey(n.x + cellSize_/2, n.y)] : &ghostHy;
  }

  // --- Hx -> Ez ---
  for (auto &n : Hx_nodes) {
    n.ez_top    = Ez_lookup.count(makeKey(n.x, n.y + cellSize_/2)) ? Ez_lookup[makeKey(n.x, n.y + cellSize_/2)] : &ghostEz;
    n.ez_bottom = Ez_lookup.count(makeKey(n.x, n.y - cellSize_/2)) ? Ez_lookup[makeKey(n.x, n.y - cellSize_/2)] : &ghostEz;
  }

  // --- Hy -> Ez ---
  for (auto &n : Hy_nodes) {
    n.ez_left   = Ez_lookup.count(makeKey(n.x - cellSize_/2, n.y)) ? Ez_lookup[makeKey(n.x - cellSize_/2, n.y)] : &ghostEz;
    n.ez_right  = Ez_lookup.count(makeKey(n.x + cellSize_/2, n.y)) ? Ez_lookup[makeKey(n.x + cellSize_/2, n.y)] : &ghostEz;
  }
}

void Mesh2D::saveMeshToVTK(const std::string& type, const std::string& filename) const {
  VTKQuadWriter vtk2elements;
  std::vector<double> coords;
  std::vector<int> connectivity;
  std::vector<double> field;

  // pick the right container
  const std::vector<EzNode> *Ez = nullptr;
  const std::vector<HxNode> *Hx = nullptr;
  const std::vector<HyNode> *Hy = nullptr;

  if (type == "Ez") Ez = &Ez_nodes;
  else if (type == "Hx") Hx = &Hx_nodes;
  else if (type == "Hy") Hy = &Hy_nodes;
  else {
    throw std::runtime_error("Unknown type: " + type);
  }

  // ---- Ez export ----
  if (Ez) {
    for (auto &n : *Ez) {
      coords.push_back(n.x);
      coords.push_back(n.y);
      coords.push_back(0.0);
      field.push_back(n.fieldValue);
    }

    for (auto &n : *Ez) {
      if (n.right && n.top && n.top->right) {
        connectivity.push_back(n.nodeID);
        connectivity.push_back(n.right->nodeID);
        connectivity.push_back(n.top->right->nodeID);
        connectivity.push_back(n.top->nodeID);
      }
    }
  }

  // ---- Hx export ----
  if (Hx) {
    for (auto &n : *Hx) {
      coords.push_back(n.x);
      coords.push_back(n.y);
      coords.push_back(0.0);
      field.push_back(n.fieldValue);
    }

    for (auto &n : *Hx) {
      if (n.right && n.top && n.top->right) {
        connectivity.push_back(n.nodeID);
        connectivity.push_back(n.right->nodeID);
        connectivity.push_back(n.top->right->nodeID);
        connectivity.push_back(n.top->nodeID);
      }
    }
    // usually Hx are stored as points, so maybe no quads?
  }

  // ---- Hy export ----
  if (Hy) {
    for (auto &n : *Hy) {
      coords.push_back(n.x);
      coords.push_back(n.y);
      coords.push_back(0.0);
      field.push_back(n.fieldValue);
    }

    for (auto &n : *Hy) {
      if (n.right && n.top && n.top->right) {
        connectivity.push_back(n.nodeID);
        connectivity.push_back(n.right->nodeID);
        connectivity.push_back(n.top->right->nodeID);
        connectivity.push_back(n.top->nodeID);
      }
    }
    // same as Hx, no quads unless you want them
  }

  vtk2elements.set_points(coords);
  vtk2elements.set_cells(connectivity);
  vtk2elements.write_vtk(filename);
}



