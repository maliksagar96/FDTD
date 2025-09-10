#include <Mesh2D.h>

int main() {
  Mesh2D mesh("../test/square_1m.STEP", 0.5);
  if (!mesh.readSTEP()) return 1;

  mesh.generateUniformGrid();
  mesh.filterNodesOutsideGeometry();
  mesh.linkNodeNeighbors();
  // mesh.generateEdges();
  mesh.saveMeshToVTK("mesh.vtk");  

  return 0;
}
