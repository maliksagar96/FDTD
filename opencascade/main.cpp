#include <STEPControl_Reader.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS.hxx>
#include <TopExp_Explorer.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <Geom_Curve.hxx>
#include <gp_Pnt.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <iostream>





bool isPointInsideFace(const TopoDS_Face& face, const gp_Pnt& point) {
  // Classify the point relative to the face
  BRepClass_FaceClassifier classifier(face, point, 1e-7); // tolerance
  TopAbs_State state = classifier.State();
  return (state == TopAbs_IN || state == TopAbs_ON);
}


int main() {
  STEPControl_Reader reader;
  if (reader.ReadFile("../square_1m.STEP") != IFSelect_RetDone) {
    std::cerr << "Cannot read STEP file." << std::endl;
    return 1;
  }

  reader.TransferRoots();
  TopoDS_Shape shape = reader.OneShape();

  // Iterate over faces
  for (TopExp_Explorer faceExp(shape, TopAbs_FACE); faceExp.More(); faceExp.Next()) {
    TopoDS_Face face = TopoDS::Face(faceExp.Current());

    // Iterate over edges
    for (TopExp_Explorer edgeExp(face, TopAbs_EDGE); edgeExp.More(); edgeExp.Next()) {
      TopoDS_Edge edge = TopoDS::Edge(edgeExp.Current());

      Standard_Real f, l;
      Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, f, l);
      if (!curve.IsNull()) {
        int n = 40; // number of points per edge
        for (int i = 0; i <= n; ++i) {
          Standard_Real param = f + (l - f) * i / n;
          gp_Pnt p = curve->Value(param);
          // std::cout << p.X() << " " << p.Y() << " " << p.Z() << "\n";
        }
      }
    }
  }

  gp_Pnt p1(0,0,0);   // inside
  gp_Pnt p2(2,2,2);   // outside

  for (TopExp_Explorer faceExp(shape, TopAbs_FACE); faceExp.More(); faceExp.Next()) {
    TopoDS_Face face = TopoDS::Face(faceExp.Current());

    std::cout << "p1 inside? " << isPointInsideFace(face, p1) << "\n";
    std::cout << "p2 inside? " << isPointInsideFace(face, p2) << "\n";
  }


  return 0;
}
