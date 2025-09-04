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

// ------------------- Read STEP file -------------------
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

// ------------------- Check if point is inside face -------------------
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
      Node n(x, y, id);   // construct node
      nodes_.push_back(n);
      id++;
    }
  }

  // optional: later we can set left/right/top/bottom by indexing
}

// ------------------- Filter cells outside geometry -------------------
void Mesh2D::filterNodesOutsideGeometry() {
  std::vector<Node> insideNodes;
  int id = 0;

  for (auto &n : nodes_) {
    gp_Pnt p(n.x, n.y, 0);  // make a point from node coords
    if (!isPointInside(p)) {
      n.nodeID = id;        // assign new sequential ID
      n.nodePoint = p;      // store gp_Pnt
      insideNodes.push_back(n);
      id++;
    }
  }

  nodes_ = insideNodes;  // keep only valid nodes
}

void Mesh2D::linkNodeNeighbors() {
  // assume uniform grid spacing = cellSize_
  std::unordered_map<std::string, Node*> lookup;

  // build quick lookup by (x,y) string key
  for (auto &n : nodes_) {
    std::string key = std::to_string(n.x) + "_" + std::to_string(n.y);
    lookup[key] = &n;
  }

  for (auto &n : nodes_) {
    // generate neighbor keys
    std::string leftKey   = std::to_string(n.x - cellSize_) + "_" + std::to_string(n.y);
    std::string rightKey  = std::to_string(n.x + cellSize_) + "_" + std::to_string(n.y);
    std::string topKey    = std::to_string(n.x) + "_" + std::to_string(n.y + cellSize_);
    std::string bottomKey = std::to_string(n.x) + "_" + std::to_string(n.y - cellSize_);

    // assign neighbors if they exist
    n.left   = lookup.count(leftKey)   ? lookup[leftKey]   : nullptr;
    n.right  = lookup.count(rightKey)  ? lookup[rightKey]  : nullptr;
    n.top    = lookup.count(topKey)    ? lookup[topKey]    : nullptr;
    n.bottom = lookup.count(bottomKey) ? lookup[bottomKey] : nullptr;
  }
}

void generateEdges() {
  std::unordered_map<long long, Edge*> edgeMap;
  auto makeKey = [](int a, int b) {
    if (a > b) std::swap(a, b);
    return (static_cast<long long>(a) << 32) | b;
  };

  int id_x = 0, id_y = 0;   // separate IDs for Hx and Hy

  for (auto &n : nodes_) {
    // horizontal edge (Hx) - vertical in space
    if (n.top) {
      long long key = makeKey(n.nodeID, n.top->nodeID);
      if (!edgeMap.count(key)) {
        Edge *e = new Edge{id_x, &n, n.top, false, 0.0};
        edges_y.push_back(e);   // vertical edges are Hy
        edgeMap[key] = e;
        id_y++;
      }
      n.Hx_top = edgeMap[key];
      n.top->Hx_bottom = edgeMap[key];
    }

    // vertical edge (Hy) - horizontal in space
    if (n.right) {
      long long key = makeKey(n.nodeID, n.right->nodeID);
      if (!edgeMap.count(key)) {
        Edge *e = new Edge{id_y, &n, n.right, true, 0.0};
        edges_x.push_back(e);   // horizontal edges are Hx
        edgeMap[key] = e;
        id_x++;
      }
      n.Hy_right = edgeMap[key];
      n.right->Hy_left = edgeMap[key];
    }
  }
}

void Mesh2D::saveMeshToVTK(const std::string& filename) const {
  VTKQuadWriter vtk2elements;
  std::vector<double> coords;
  std::vector<int> connectivity;
  std::vector<double> Ez_scalar;

  // Export node coordinates
  for (auto &n : nodes_) {
    coords.push_back(n.x);
    coords.push_back(n.y);
    coords.push_back(0.0);
    Ez_scalar.push_back(n.Ez);
  }

  // Build quad connectivity
  for (auto &n : nodes_) {
     if (n.right && n.top && n.top->right) {
      
      connectivity.push_back(n.nodeID);
      connectivity.push_back(n.right->nodeID);
      connectivity.push_back(n.top->right->nodeID);
      connectivity.push_back(n.top->nodeID);
      
    }
  }

  vtk2elements.set_points(coords);
  vtk2elements.set_cells(connectivity);
  vtk2elements.write_vtk(filename);
}


