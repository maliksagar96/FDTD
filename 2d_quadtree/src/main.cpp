#include <Mesh2D.h>

int main() {
  Mesh2D mesh("../test/square_1m.STEP", 0.05);
  if (!mesh.readSTEP()) return 1;

  mesh.generateUniformGrid();
  mesh.filterCellsOutsideGeometry();
  mesh.saveMeshToTXT("mesh.txt");

  return 0;
}
