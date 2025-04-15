///////////////////////////////////////////////////////////////////////
// Generate a third force grid.
// Author: Jeff Comer <jcomer2@illinois.edu>
#ifndef THIRDFORCEGRID_H
#define THIRDFORCEGRID_H

#include "Grid.H"
#include "CellDecomposition.H"
#include "Scatter.H"

class ThirdForceGrid : public Grid {
public:
  ThirdForceGrid(Matrix3 box0, Vector3 origin0, double radius0, double sigma0, double dx) : Grid(box0, origin0, dx) {
    radius = radius0;
    sigma = sigma0;
    
    box = box0;
  }

  // Compute the grid using a cell decomposition, computing the energy between
  // grid points and the source point in the neighborhood of the grid point.
  void compute(const Scatter& source, const CellDecomposition& decomp) {
    int j = 0;

    //IndexList one = decomp.neighborhood(Vector3(0.0));
    for (int ix = 0; ix < nx; ix++) {
      for (int iy = 0; iy < ny; iy++) {
	for (int iz = 0; iz < nz; iz++) {
	  Vector3 r = basis.transform(Vector3(ix, iy, iz)) + origin;
	  val[j] = localEnergy(r, source, decomp.neighborhood(r));
	  j++;
	}
      }
      if (ix % 5 == 0) printf("%2.1f percent complete\n", (100.0*j)/size);
    }
  }

  // Compute the energy at position r due to specified source points.
  double localEnergy(Vector3 r, const Scatter& source, const IndexList& neigh) const {
    const int n = neigh.length();
    //printf("%s\n\n",(const char*)neigh.toString());
   
    double e = 0.0;
    for (int i = 0; i < n; i++) {
      Vector3 r0 = source.get(neigh.get(i));
      // Transform the vector between the grid node and the source point into lattice space.
      Vector3 l = basisInv.transform(r - r0);
      // Wrap it according to the periodic boundaries.
      if (l.x < 0.5*nx) {l.x += nx;}
      if (l.x > 0.5*nx) {l.x -= nx;}
      if (l.y < 0.5*ny) {l.y += ny;}
      if (l.y > 0.5*ny) {l.y -= ny;}
      if (l.z < 0.5*nz) {l.z += nz;}
      if (l.z > 0.5*nz) {l.z -= nz;}
      
      // Transform the wrapped vector back into world space.
      Vector3 d = basis.transform(l);
      
      // Acculumulate the energy.
      e += energy(d.length());
    }

    //printf("%.10g\n", e);
    return e;
  }

  // Compute the grid by computing energy between all grid points and all sources :(.
  void compute(const Scatter& source) {
    int j = 0;

    for (int ix = 0; ix < nx; ix++) {
      for (int iy = 0; iy < ny; iy++) {
	for (int iz = 0; iz < nz; iz++) {
	  Vector3 r = basis.transform(Vector3(ix, iy, iz)) + origin;
	  val[j] = globalEnergy(r, source);

	  j++;
	}
      }
      if (ix % 5 == 0) printf("%2.1f percent complete\n", (100.0*j)/size);
    }
  }
 
  // Compute the energy at position r due to all source points.
  // This can be slow for large systems.
  double globalEnergy(Vector3 r, const Scatter& source) const {
    const int n = source.length();
    
    double e = 0.0;
    for (int i = 0; i < n; i++) {
      Vector3 r0 = source.get(i);
      // Transform the vector between the grid node and the source point into lattice space.
      Vector3 l = basisInv.transform(r - r0);
      // Wrap it according to the periodic boundaries.
      if (l.x < 0.5*nx) {l.x += nx;}
      if (l.x > 0.5*nx) {l.x -= nx;}
      if (l.y < 0.5*ny) {l.y += ny;}
      if (l.y > 0.5*ny) {l.y -= ny;}
      if (l.z < 0.5*nz) {l.z += nz;}
      if (l.z > 0.5*nz) {l.z -= nz;}
      
      // Transform the wrapped vector back into world space.
      Vector3 d = basis.transform(l);
      
      // Acculumulate the energy.
      e += energy(d.length());
    }
    return e;
  }
   
  double energy(double r) const {
    if (r < radius) return radius - r + 0.5*sigma;
    if (r > radius + sigma) return 0.0;
    
    double dr = r - radius;
    return dr*(0.5*dr/sigma - 1.0) + 0.5*sigma;
  }

private:
  double radius, sigma;
  Matrix3 box; // periodic box of the system
  //ThirdForceGrid(const ThirdForceGrid&) {}
  //ThirdForceGrid operator=(const ThirdForceGrid&) {}
};
#endif
