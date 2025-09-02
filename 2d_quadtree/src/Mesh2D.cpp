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

  for (double x = xmin; x < xmax; x += cellSize_) {
    for (double y = ymin; y < ymax; y += cellSize_) {
      Cell c;
      c.x_min = x;
      c.x_max = x + cellSize_;
      c.y_min = y;
      c.y_max = y + cellSize_;
      c.center = gp_Pnt(x + cellSize_/2, y + cellSize_/2, 0);
      cells_.push_back(c);
    }
  }
}

// ------------------- Filter cells outside geometry -------------------
void Mesh2D::filterCellsOutsideGeometry() {
  std::vector<Cell> outsideCells;
  for (const auto& c : cells_) {
    if (!isPointInside(c.center)) {
      outsideCells.push_back(c);
    }
  }
  cells_ = outsideCells;
}

// ------------------- Save mesh to TXT -------------------
void Mesh2D::saveMeshToTXT(const std::string& filename) const {
  std::ofstream fout(filename);
  if (!fout.is_open()) {
    std::cerr << "Cannot open file to write: " << filename << std::endl;
    return;
  }

  for (const auto& c : cells_) {
    fout << c.center.X() << " " << c.center.Y() << " " << c.center.Z() << "\n";
  }

  fout.close();
  std::cout << "Mesh saved to " << filename << " with " << cells_.size() << " points." << std::endl;
}
