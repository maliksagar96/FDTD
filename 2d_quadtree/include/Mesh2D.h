#pragma once

#include <STEPControl_Reader.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Pnt.hxx>
#include <vector>
#include <string>

struct Node {
  double x, y;
  double fieldValue;
  int nodeID;
  gp_Pnt nodePoint;
  Node *left, *right, *top, *bottom;
  Node *Hx_left, *Hx_right, *Hy_top, *Hy_bottom;

  Node(double x_, double y_, int id)
    : x(x_), y(y_), nodeID(id), fieldValue(0), nodePoint(x_, y_, 0.0),
      left(nullptr), right(nullptr), top(nullptr), bottom(nullptr),
      Hx_left(nullptr), Hx_right(nullptr),
      Hy_top(nullptr), Hy_bottom(nullptr) {}
};


class Mesh2D {
  public:
    Mesh2D(const std::string& stepFile, double cellSize);
    bool readSTEP();                              // read STEP file and store TopoDS_Face
    void generateUniformGrid();                   // create uniform grid of nodes
    void filterNodesOutsideGeometry();            // keep only inside nodes
    void linkNodeNeighbors();                     // set left/right/top/bottom
    
    void saveMeshToVTK(const std::string& filename) const;
    void save_Hx_MeshToVTK(const std::string& filename) const;
    const std::vector<Node>& get_Ez_Nodes() const { return Ez_nodes; }
    const std::vector<Node>& get_Hx_Nodes() const { return Hx_nodes; }
    const std::vector<Node>& get_Hy_Nodes() const { return Hy_nodes; }
    
    void generateEdges();

  private:
    std::string stepFile_;
    double cellSize_;
    TopoDS_Face face_;                            // assume single face for 2D
    std::vector<Node> Ez_nodes;
    std::vector<Node> Hx_nodes;
    std::vector<Node> Hy_nodes;

    bool isPointInside(const gp_Pnt& p) const;    // helper for inside/outside check
};
