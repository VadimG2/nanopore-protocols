// Author: Jeff Comer <jcomer2@illinois.edu>
// Writes the positions of grid points outside the pore.
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
using namespace std;

#include "useful.H"
#include "Grid.H"
#include "Scatter.H"

// Is a point in the material making up the membrane? 
bool inPoreWalls(Vector3 r, double l, double s0, double slope) {
  if (fabs(r.z) > l) return false;
  
  double s = sqrt(r.x*r.x + r.y*r.y);
  if (s < s0 + slope*fabs(r.z)) return false;

  return true;
}

///////////////////////////////////////////////////////////////////////
// Driver
int main(int argc, char* argv[]) {
  if ( argc < 7 ) {
    printf("Usage: %s systemCellFile gridSize poreLength poreDiameter poreAngle outFile\n", argv[0]);
    return 0;
  }
  
  const char* cellFile = argv[1];
  double gridSize = strtod(argv[2], NULL);
  double poreLength = strtod(argv[3], NULL);
  double poreDiameter = strtod(argv[4], NULL);
  double poreAngle = strtod(argv[5], NULL);
  const char* outFile = argv[argc-1];
  
  // Load the system cell.
  Scatter b(cellFile);
  if (b.length() < 3) {
    printf("Error! Invalid systemCellFile `%s'.\n", cellFile);
    printf("File should contain the system's basis vectors.\n");
    printf("Format is:\n");
    printf("ax ay az\n");
    printf("bx by bz\n");
    printf("cx cy cz\n");
  }

  // Make the grid.
  Matrix3 basis = b.topMatrix();
  Grid g(basis, gridSize);

  // Open the file.
  FILE* out = fopen(outFile, "w");
  double l = 0.5*poreLength;
  double s0 = 0.5*poreDiameter;
  double pi = 4.0*atan(1.0);
  double slope = tan(pi*poreAngle/180.0);

  printf("Scanning points.\n");

  int count = 0;
  const int n = g.length();
  for (int i = 0; i < n; i++) {
    Vector3 r(g.getPosition(i));
    
    if (!inPoreWalls(r, l, s0, slope)) {
      // Write this point.
      fprintf(out, "%s\n", r.toString().val());
      count++;
    }
  }
  fclose(out);
  
  double totalVol = g.getVolume();
  double srcVol = (totalVol*count)/n;
  double remainVol = totalVol - srcVol;

  printf("Total points: %d\n", n);
  printf("Source points: %d\n", count);
  printf("Remaining points: %d\n", n-count);
  printf("Source fraction: %.10g\n", double(count)/n);
  printf("Total volume: %.10g\n", totalVol);
  printf("Source volume: %.10g\n", srcVol);
  printf("Remaining volume: %.10g\n", remainVol);
  
  return 0;
}
