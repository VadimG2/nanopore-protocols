
///////////////////////////////////////////////////////////////////////
// Cell decomposition of points.
// Author: Jeff Comer <jcomer2@illinois.edu>
#ifndef CELLDECOMPOSITION_H
#define CELLDECOMPOSITION_H

#include "Grid.H"
#include "Scatter.H"

struct SpaceCell {
  IndexList point;
  IndexList neigh;
  int ix, iy, iz;
  int index;
};

class CellDecomposition : public Grid {
public:
  // Make an empty cell decomposition with the desired geometry.
  CellDecomposition(Matrix3 box0, Vector3 origin0, double cutoff) : Grid(box0, origin0, cutoff) {

    cut = cutoff;
    cell = new SpaceCell[size];
    prepareCells();
    populateNeighbors();
  }

  CellDecomposition(Vector3 box, double cutoff) : Grid(Matrix3(box.x,box.y,box.z), -0.5*box, cutoff) {
    cut = cutoff;
    cell = new SpaceCell[size];
    prepareCells();
    populateNeighbors();
  }

  virtual ~CellDecomposition() {
    delete[] cell;
  }

  double getCutoff() const {
    return cut;
  }

  // Find the cell index for a given position.
  int getCell(const Vector3& r) const {
    // Transform.
    Vector3 l = basisInv.transform(r-origin);
    // Wrap into the home box.
    int cx = wrap(int(floor(l.x+0.5)), nx);
    int cy = wrap(int(floor(l.y+0.5)), ny);
    int cz = wrap(int(floor(l.z+0.5)), nz);
    
    int j = cz + cy*nz + cx*nz*ny;

    return j;
  }
  
  IndexList getCellContents(int i) {
#ifdef DEBUG 
    if (i < 0 || i >= size) {
      printf("Warning! CellDecomposition::getCellContents out of bounds.\n");
      return IndexList();
    }
#endif
    return cell[i].point;
  }

  // Restart with an empty decomposition.
  void prepareCells() {
    int j = 0;

    for (int ix = 0; ix < nx; ix++) {
      for (int iy = 0; iy < ny; iy++) {
	for (int iz = 0; iz < nz; iz++) {
	  cell[j].ix = ix;
	  cell[j].iy = iy;
	  cell[j].iz = iz;
	  cell[j].index = j;
	  cell[j].point.clear();
	  val[j] = 0.0;

	  j++;
	}
      }
    }
  }

  // Restart with an empty decomposition.
  void clearCells() {
    for (int i = 0; i < size; i++) {
      cell[i].point.clear();
      val[i] = 0.0;
    }
  }

  // Generate the neighbor lists.
  void populateNeighbors() {
    int j = 0;

    for (int ix = 0; ix < nx; ix++) {
      for (int iy = 0; iy < ny; iy++) {
	for (int iz = 0; iz < nz; iz++) {
	  cell[j].neigh = getNeighbors(j);
	  j++;
	}
      }
    }
  }  
  
  void decompose(const Vector3* r, int n) {
    for (int i = 0; i < n; i++) addPoint(r[i], i);
  }

  void decompose(const Scatter& pos) {
    const int n = pos.length();

    for (int i = 0; i < n; i++) addPoint(pos.get(i), i);
  }

  void addPoint(const Vector3& r, int pointIndex) {
    int j = getCell(r);
    cell[j].point.add(pointIndex);
    val[j] += 1.0;
  }

  int countPoints() const {
    int num = 0;
    for (int j = 0; j < size; j++) {
      num += cell[j].point.length();
    }
    return num;
  }

  IndexList neighborhood(Vector3 r) const {
    int j = getCell(r);
    int nNeighs = cell[j].neigh.length();

    IndexList ret;

    // Get the point indices in all neighbors of the home cell (inclusive).
    for (int n = 0; n < nNeighs; n++) {
      int neigh = cell[j].neigh.get(n);
      
      ret.add(cell[neigh].point);
    }

    return ret;
  }

protected:
  double cut; // cutoff length
  SpaceCell* cell; // contains cell attributes

  // Fill the neighbor list with the cell's neighbors (including itself).
  IndexList getNeighbors(int j) {
    const int nynz = ny*nz;
    int ix, iy, iz;
    IndexList ret;

    for (ix = -1; ix <= 1; ix++) {
      for (iy = -1; iy <= 1; iy++) {
	for (iz = -1; iz <= 1; iz++) {
	  int cx = cell[j].ix + ix;
	  if (cx >= nx) cx -= nx;
	  if (cx < 0) cx += nx;

	  int cy = cell[j].iy + iy;
	  if (cy >= ny) cy -= ny;
	  if (cy < 0) cy += ny;

	  int cz = cell[j].iz + iz;
	  if (cz >= nz) cz -= nz;
	  if (cz < 0) cz += nz;

	  int ind = cz + cy*nz + cx*nynz;
	  if (ret.find(ind) == -1) ret.add(ind);
	}
      }
    }

    return ret;
  }

private:  

  //CellDecomposition(const CellDecomposition&) {}
  //CellDecomposition operator=(const CellDecomposition&) {}
};
#endif
