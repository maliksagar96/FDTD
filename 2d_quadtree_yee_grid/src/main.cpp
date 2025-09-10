#include <Mesh2D.h>
#include <fdtd_2d.h>

using namespace std;

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <control_file.scf>" << std::endl;
    return 1;
  }

  std::ifstream file(argv[1], std::ifstream::binary);
  if (!file) {
    std::cerr << "Error: could not open file " << argv[1] << std::endl;
    return 1;
  }

  Json::Value root;
  file >> root;

  FDTD_2D solver(root);

  Mesh2D mesh("../test/square_1m.STEP", 0.05);
  if (!mesh.readSTEP()) return 1;

  mesh.generateUniformGrid();
  mesh.filterNodesOutsideGeometry();
  mesh.linkNodeNeighbors();
  mesh.linkCrossNeighbors();
  mesh.saveMeshToVTK("Ez", "Ezmesh.vtk");  
  mesh.saveMeshToVTK("Hx", "Hxmesh.vtk");  
  mesh.saveMeshToVTK("Hy", "Hymesh.vtk");  

  vector<EzNode> Ez_nodes = mesh.get_Ez_Nodes();
  vector<HxNode> Hx_nodes = mesh.get_Hx_Nodes();
  vector<HyNode> Hy_nodes = mesh.get_Hy_Nodes(); 

  solver.get_fields(Ez_nodes, Hx_nodes, Hy_nodes);
  solver.TMz_mesh_update();
  return 0;
}
