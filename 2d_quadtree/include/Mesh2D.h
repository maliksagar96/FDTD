#pragma once

#include <STEPControl_Reader.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Pnt.hxx>
#include <vector>
#include <string>

struct Node;

struct Edge {
  int id;
  Node *n1, *n2;   // endpoints
  bool isX;        // true if horizontal (Hx), false if vertical (Hy)
  double H;
};

struct Node {
  double x, y;
  double Ez_fieldValue;
  int nodeID;
  gp_Pnt nodePoint;
  Node *left, *right, *top, *bottom;

  Edge *Hx_top, *Hx_bottom;
  Edge *Hy_left, *Hy_right; 

  

  Node(double x_, double y_, int id)
    : x(x_), y(y_), nodeID(id), Ez_fieldValue(0), nodePoint(x_, y_, 0.0),
      left(nullptr), right(nullptr), top(nullptr), bottom(nullptr),
      Hx_top(nullptr), Hx_bottom(nullptr),
      Hy_left(nullptr), Hy_right(nullptr) {}
};


class Mesh2D {
  public:
    Mesh2D(const std::string& stepFile, double cellSize);
    bool readSTEP();                              // read STEP file and store TopoDS_Face
    void generateUniformGrid();                   // create uniform grid of nodes
    void filterNodesOutsideGeometry();            // keep only inside nodes
    void linkNodeNeighbors();                     // set left/right/top/bottom
    void saveMeshToTXT(const std::string&) const; // save node coords
    void saveMeshToVTK(const std::string& filename) const;
    const std::vector<Node>& getNodes() const { return nodes_; }
    const std::vector<Edge>& getEdges() const { return edges_; }
    void generateEdges();

  private:
    std::string stepFile_;
    double cellSize_;
    TopoDS_Face face_;                            // assume single face for 2D
    std::vector<Node> nodes_;
    std::vector<Edge*> edges_x; // all Hx edges
    std::vector<Edge*> edges_y; // all Hy edges

    bool isPointInside(const gp_Pnt& p) const;    // helper for inside/outside check
};
