// Author: Jeff Comer <jcomer2@illinois.edu>
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include "useful.H"
#include "Grid.H"

using namespace std;

///////////////////////////////////////////////////////////////////////
// Driver
int main(int argc, char* argv[]) {
  if ( argc != 6 ) {
    printf("Remove nx ny nz grid points along the respective direction.\n");
    printf("Usage: %s srcGrid nx ny nz outFile\n", argv[0]);
    return 0;
  }

  const char* inFile = argv[1];
  int shaveX = atoi(argv[2]);
  int shaveY = atoi(argv[3]);
  int shaveZ = atoi(argv[4]);
  const char* outFile = argv[argc-1];

  // Load the first grid.
  Grid src(inFile);
  printf("Loaded %s.\n", argv[1]);

  int nx = src.getNx();
  int ny = src.getNy();
  int nz = src.getNz();
  Matrix3 basis = src.getBasis();
  Vector3 origin = src.getOrigin();
  nx -= shaveX;
  ny -= shaveY;
  nz -= shaveZ;

  int hx = shaveX/2;
  int hy = shaveY/2;
  int hz = shaveZ/2;
  Vector3 o = origin + basis.transform(Vector3(hx,hy,hz));

  Grid out(basis, o, nx, ny, nz);
  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {
      for (int iz = 0; iz < nz; iz++) {
	double v = src.getValue(ix+hx, iy+hy, iz+hz);
	out.setValue(ix, iy, iz, v);
      }
    }
  }
  
  char comments[256];
  sprintf(comments, "%s shaved %d %d %d", inFile, shaveX, shaveY, shaveZ);
  out.write(outFile, comments);
  printf("Wrote `%s'.\n", outFile);

  return 0;
}
