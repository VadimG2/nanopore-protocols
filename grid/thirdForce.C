///////////////////////////////////////////////////////////////////////
// Generate a third force grid.
// Author: Jeff Comer <jcomer2@illinois.edu>

#include <cstdio>
#include <cstdlib>
using namespace std;

#include "useful.H"
#include "Scatter.H"
#include "ThirdForceGrid.H"

///////////////////////////////////////////////////////////////////////
// Driver
int mainThirdForce(int argc, char* argv[])
{
  if ( argc != 7 ) {
    printf(
      "Usage: %s coordFile basisFile dx radius0 sigma outFile\n",
      argv[0]);
    return 0;
  }

  // Extract the parameters.
  double dx = strtod(argv[3],NULL);
  double radius = strtod(argv[4],NULL);
  double sigma = strtod(argv[5],NULL);
  printf("\n********************************\nThird force grid: \n");
  printf("The resolution is %.10g.\n", dx);
  printf("Creating a grid with radius0 = %.10g and sigma = %.10g...\n",
	 radius, sigma);

  // Get the system dimensions.
  Scatter sysVec(argv[2]);
  Matrix3 sys = sysVec.topMatrix();
  printf("System dimensions:\n%s\n", (const char*)sys.toString());
  Vector3 origin = (sys.ex() + sys.ey() + sys.ez())*-0.5;
  printf("System origin:%s\n", (const char*)origin.toString());

  // Get the source points.
  Scatter source(argv[1]);
  printf("\nLoaded %i source points from %s\n", source.length(), argv[1]);

  // Generate the cell decomposition.
  double cutoff = 1.8*(radius+sigma);
  printf("\nGenerating cell decomposition for %s with cutoff = %.10g...\n", argv[1], cutoff);
  CellDecomposition cell(sys, origin, cutoff);
  printf("Cell decomposition of %i cells\n", cell.length());
  cell.decompose(source);
  printf("Cell decomposition contains %i points\n", cell.countPoints());

  // Validate the cell decomposition.
  //const int nCells = cell.length();
  //const int nPoints = source.length();
  
  // Organize the file comments.
  char comments[256];
  sprintf(comments, "ThirdForce grid for %s with dx = %.10g, radius0 = %.10g, and sigma = %.10g", argv[1], dx, radius, sigma);
    
  // Generate the grid.
  ThirdForceGrid third(sys, origin, radius, sigma, dx);
  printf("\nComputing the grid...\n");
  third.compute(source, cell);
  //third.compute(source);
  printf("Writing the grid...\n");
  third.write(argv[6], comments);
  printf("Wrote %s successfully.\n", argv[6]);

  return 0;
}

int main(int argc, char* argv[]) {
  return mainThirdForce(argc, argv);
}
