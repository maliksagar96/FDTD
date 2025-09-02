#pragma once

#include <STEPControl_Reader.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Pnt.hxx>
#include <vector>
#include <string>

struct Cell {
  double x_min, x_max;
  double y_min, y_max;
  gp_Pnt center;
};

class Mesh2D {
  public:
    Mesh2D(const std::string& stepFile, double cellSize);
    bool readSTEP();                           // read STEP file and store TopoDS_Face
    void generateUniformGrid();                // create uniform cells around geometry
    void filterCellsOutsideGeometry();        // remove cells inside geometry
    void saveMeshToTXT(const std::string& filename) const; // save mesh points

    const std::vector<Cell>& getCells() const { return cells_; }

  private:
    std::string stepFile_;
    double cellSize_;
    TopoDS_Face face_;                       // assume single face for 2D
    std::vector<Cell> cells_;
    bool isPointInside(const gp_Pnt& p) const; // helper for inside/outside check
};
