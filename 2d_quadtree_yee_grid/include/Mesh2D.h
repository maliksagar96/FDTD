#pragma once

#include <STEPControl_Reader.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Pnt.hxx>
#include <vector>
#include <string>

class HxNode;
class HyNode;
class EzNode;

class BaseNode {
public:
  BaseNode(double x_, double y_, int id)
    : x(x_), y(y_), nodeID(id),
      left(nullptr), right(nullptr), top(nullptr), bottom(nullptr) {}
  virtual ~BaseNode() = default;

  int nodeID;
  double x, y;
  double fieldValue;
  gp_Pnt nodePoint;
  BaseNode *left, *right, *top, *bottom;  
};

class EzNode : public BaseNode {
  public:
    EzNode(double x, double y, int id) : BaseNode(x, y, id) {}
    HxNode *hx_top, *hx_bottom;
    HyNode *hy_left, *hy_right;
};

class HxNode : public BaseNode {
public:
  HxNode(double x, double y, int id) : BaseNode(x, y, id) {}
  EzNode *ez_top, *ez_bottom;
};

class HyNode : public BaseNode {
public:
  HyNode(double x, double y, int id) : BaseNode(x, y, id) {}
  EzNode *ez_left, *ez_right;
};

class Mesh2D {
  public:
    Mesh2D(const std::string& stepFile, double cellSize);
    bool readSTEP();                              // read STEP file and store TopoDS_Face
    void generateUniformGrid();                   // create uniform grid of nodes
    void filterNodesOutsideGeometry();            // keep only inside nodes
    void linkNodeNeighbors();                     // set left/right/top/bottom
    void linkCrossNeighbors();
    void saveMeshToVTK(const std::string&, const std::string&) const;
    void save_Hx_MeshToVTK(const std::string& filename) const;    
    const std::vector<EzNode>& get_Ez_Nodes() const { return Ez_nodes; }
    const std::vector<HxNode>& get_Hx_Nodes() const { return Hx_nodes; }
    const std::vector<HyNode>& get_Hy_Nodes() const { return Hy_nodes; }
    
    void generateEdges();

  private:
    
    bool isPointInside(const gp_Pnt& p) const;    // helper for inside/outside check
    
    std::string stepFile_;
    double cellSize_;
    TopoDS_Face face_;                            
    std::vector<EzNode> Ez_nodes;
    std::vector<HxNode> Hx_nodes;
    std::vector<HyNode> Hy_nodes;    
};
