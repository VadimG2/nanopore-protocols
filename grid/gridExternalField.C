// Author: Jeff Comer <jcomer2@illinois.edu>
// Apply a uniform external field to a grid.
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include "useful.H"
#include "Grid.H"

using namespace std;

///////////////////////////////////////////////////////////////////////
// Driver
int main(int argc, char* argv[]) {
  if ( argc != 5 ) {
    printf("Usage: %s srcGrid scaleFactor voltageDrop outFile\n", argv[0]);
    return 0;
  }

  const char* inFile = argv[1];
  double scaleFactor = strtod(argv[2], NULL);
  double voltage = strtod(argv[3], NULL);
  const char* outFile = argv[argc-1];

  // Load the first grid.
  Grid src(inFile);
  printf("Loaded %s.\n", argv[1]);

  // Scale.
  printf("Scaling by %g.\n", scaleFactor);
  src.scale(scaleFactor);

  // Write the profile.
  char srcProfile[256];
  sprintf(srcProfile, "%s.profile.dat", inFile);
  printf("Writing a profile along z to `%s'.\n", srcProfile);
  src.profileZ(0.5, 0.5, srcProfile);

  // Add the gradient.
  printf("Adding voltage drop of %g.\n", voltage);
  double lz = src.getExtent().z;
  printf("The length of the system is %g.\n", lz);
  src.addGradient(Vector3(0.0,0.0,-voltage/lz));

  // Write the new grid.
  char comments[256];
  sprintf(comments, "%s with voltage drop of %g", argv[1], voltage);
  printf("Writing the new grid `%s'.\n", outFile);
  src.write(argv[argc-1], comments);
  printf("%s\n", comments);

  char destProfile[256];
  sprintf(destProfile, "%s.profile.dat", outFile);
  printf("Writing a profile along z to `%s'.\n", destProfile);
  src.profileZ(0.5, 0.5, destProfile);

  return 0;
}
