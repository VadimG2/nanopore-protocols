//////////////////////////////////////////////////////////////////////
// Grid base class.
// Author: Jeff Comer <jcomer2@illinois.edu>
#ifndef GRID_H
#define GRID_H

#include "useful.H"
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>

using namespace std;

class Grid {
public:
  Grid(Matrix3 basis0, Vector3 origin0, int nx0, int ny0, int nz0) {
    basis = basis0;
    origin = origin0;
    nx = abs(nx0);
    ny = abs(ny0);
    nz = abs(nz0);
    
    basisInv = basis.inverse();
    size = nx*ny*nz;
    val = new double[size];
    zero();
  }

  Grid(Vector3 box, double dx) {
    dx = fabs(dx);
    box.x = fabs(box.x);
    box.y = fabs(box.y);
    box.z = fabs(box.z);

    // Tile the grid into the system box.
    // The grid spacing is always a bit smaller than dx.
    nx = int(ceil(box.x/dx));
    ny = int(ceil(box.y/dx));
    nz = int(ceil(box.z/dx));
    if (nx <= 0) nx = 1;
    if (ny <= 0) ny = 1;
    if (nz <= 0) nz = 1;
    basis = Matrix3(box.x/nx, box.y/ny, box.z/nz);
    origin = -0.5*box;

    basisInv = basis.inverse();
    size = nx*ny*nz;
    val = new double[size];
    zero();
  }

  Grid(Matrix3 box, int nx0, int ny0, int nz0) {
    nx = nx0;
    ny = ny0;
    nz = nz0;

    // Tile the grid into the system box.
    if (nx <= 0) nx = 1;
    if (ny <= 0) ny = 1;
    if (nz <= 0) nz = 1;
    basis = Matrix3(box.ex()/nx, box.ey()/ny, box.ez()/nz);
    origin = -0.5*(box.ex() + box.ey() + box.ez());

    basisInv = basis.inverse();
    size = nx*ny*nz;
    val = new double[size];
    zero();
  }

  Grid(Matrix3 box, Vector3 origin0, double dx) {
    dx = fabs(dx);
    
    // Tile the grid into the system box.
    // The grid spacing is always a bit larger than dx.
    nx = int(floor(box.ex().length()/dx))-1;
    ny = int(floor(box.ey().length()/dx))-1;
    nz = int(floor(box.ez().length()/dx))-1;
    if (nx <= 0) nx = 1;
    if (ny <= 0) ny = 1;
    if (nz <= 0) nz = 1;

    basis = Matrix3(box.ex()/nx, box.ey()/ny, box.ez()/nz);
    origin = origin0;

    basisInv = basis.inverse();
    size = nx*ny*nz;
    val = new double[size];
    zero();
  }

  Grid(Matrix3 box, double dx) {
    dx = fabs(dx);
    
    // Tile the grid into the system box.
    // The grid spacing is always a bit smaller than dx.
    nx = int(ceil(box.ex().length()/dx));
    ny = int(ceil(box.ey().length()/dx));
    nz = int(ceil(box.ez().length()/dx));
    if (nx <= 0) nx = 1;
    if (ny <= 0) ny = 1;
    if (nz <= 0) nz = 1;

    basis = Matrix3(box.ex()/nx, box.ey()/ny, box.ez()/nz);
    origin = -0.5*(box.ex() + box.ey() + box.ez());

    basisInv = basis.inverse();
    size = nx*ny*nz;
    val = new double[size];
    zero();
  }

  Grid(const Grid& g) {
    nx = g.nx;
    ny = g.ny;
    nz = g.nz;
    basis = g.basis;
    origin = g.origin;
    
    basisInv = g.basisInv;
    size = g.size;
    val = new double[size];
    for (int i = 0; i < size; i++) val[i] = g.val[i];
  }

  Grid& operator=(const Grid& g) {
    delete[] val;

    nx = g.nx;
    ny = g.ny;
    nz = g.nz;
    basis = g.basis;
    origin = g.origin;
    
    basisInv = g.basisInv;
    size = g.size;
    val = new double[size];
    for (int i = 0; i < size; i++) val[i] = g.val[i];

    return *this;
  }

  // Read a grid from a file.
  Grid(const char* fileName) {
    // Open the file.
    FILE* inp = fopen(fileName,"r");
    if (inp == NULL) {
      printf("Grid:Grid Couldn't open file %s.\n",fileName);
      exit(-1);
    }
    //printf("Reading dx file %s...\n", fileName);
    
    size = 0;
    nx = 0;
    ny = 0;
    nz = 0;
    basis = Matrix3(1.0);
    origin = Vector3(0.0);    

    int n = 0;
    double x, y, z;
    char line[256];
    int p, nRead;
    int deltaCount = 0;
    Vector3 base[3];
    while (fgets(line, 256, inp) != NULL) {
      // Ignore comments.
      int len = strlen(line);
      if (line[0] == '#') continue;
      if (len < 2) continue;
    
      if (isInt(line[0]) && n < size) {
	// Read grid values.
	nRead = sscanf(line, "%lf %lf %lf", &x, &y, &z);
	if (size > 0) {
	  switch(nRead) {
	  case 1:
	    val[n] = x;
	    n++;
	    break;
	  case 2:
	    val[n] = x;
	    val[n+1] = y;
	    n += 2;
	    break;
	  case 3:
	    val[n] = x;
	    val[n+1] = y;
	    val[n+2] = z;
	    n += 3;
	    break;
	  }
	}
      } else if (len > 5) {
	// Read the grid parameters.
	char start[6];
	for (int i = 0; i < 5; i++) start[i] = line[i];
	start[5] = '\0';

	if(strcmp("origi", start) == 0) {
	  // Get an origin line.
	  p = firstSpace(line, 256);
	  sscanf(&(line[p+1]), "%lf %lf %lf", &x, &y, &z);
	  origin = Vector3(x, y, z);
	  //printf("Origin: %.10g %.10g %.10g\n", x, y, z);
	} else if(strcmp("delta", start) == 0) {
	  // Get a delta matrix line.
	  p = firstSpace(line, 256);
	  sscanf(&(line[p+1]), "%lf %lf %lf", &x, &y, &z);
	  base[deltaCount] = Vector3(x, y, z);
	  //printf("Delta %d: %.10g %.10g %.10g\n", deltaCount, x, y, z);
	  if (deltaCount < 2) deltaCount = deltaCount + 1;
	} else if(strcmp("objec", start) == 0) {
	  //printf("%s", line);
	  // Get the system dimensions.
	  if (line[7] != '1') continue;
	  int read = sscanf(line, "object 1 class gridpositions counts %d %d %d\n", &nx, &ny, &nz);
	  //printf("Size: %d %d %d\n", nx, ny, nz);
	  if (read == 3) {
	    size = nx*ny*nz;
	    val = new double[size];
	    zero();
	  }
	}
      }
    }
    fclose(inp);

    basis = Matrix3(base[0], base[1], base[2]);
    basisInv = basis.inverse();
    if (size == 0) {
      printf("Error Grid:Grid Improperly formatted dx file %s.\n",fileName);
      exit(-1);
    }
  }
 
  // Write without comments.
  virtual void write(const char* fileName) const {
    write(fileName, "");
  }

  // Writes the grid as a file in the dx format.
  virtual void write(const char* fileName, const char* comments) const {
    // Open the file.
    FILE* out = fopen(fileName,"w");
    if (out == NULL) {
      printf("Couldn't open file %s.\n",fileName);
      exit(-1);
    }

    // Write the header.
    fprintf(out, "# %s\n", comments);
    fprintf(out, "object 1 class gridpositions counts %d %d %d\n", nx, ny, nz);
    fprintf(out, "origin %.10g %.10g %.10g\n", origin.x, origin.y, origin.z);
    fprintf(out, "delta %.10g %.10g %.10g\n", basis.exx, basis.eyx, basis.ezx);
    fprintf(out, "delta %.10g %.10g %.10g\n", basis.exy, basis.eyy, basis.ezy);
    fprintf(out, "delta %.10g %.10g %.10g\n", basis.exz, basis.eyz, basis.ezz);
    fprintf(out, "object 2 class gridconnections counts %d %d %d\n", nx, ny, nz);
    fprintf(out, "object 3 class array type double rank 0 items %d data follows\n", size);
    
    // Write the data.
    int penultima = 3*(size/3);
    int mod = size - penultima;

    int i;
    for (i = 0; i < penultima; i+=3) {
      fprintf(out, "%.10g %.10g %.10g\n", val[i], val[i+1], val[i+2]);
    }
    if (mod == 1) {
      fprintf(out, "%.10g\n", val[size-1]);
    } else if (mod == 2) {
      fprintf(out, "%.10g %.10g\n", val[size-2], val[size-1]);
    }
    fclose(out);
  }

  // Writes the grid data as a single column in the order:
  // nx ny nz ox oy oz dxx dyx dzx dxy dyy dzy dxz dyz dzz val0 val1 val2 ...
  virtual void writeData(const char* fileName) {
    // Open the file.
    FILE* out = fopen(fileName,"w");
    if (out == NULL) {
      printf("Couldn't open file %s.\n",fileName);
      exit(-1);
    }

    fprintf(out, "%d\n%d\n%d\n", nx, ny, nz);
    fprintf(out, "%.10g\n%.10g\n%.10g\n", origin.x, origin.y, origin.z);
    fprintf(out, "%.10g\n%.10g\n%.10g\n", basis.exx, basis.eyx, basis.ezx);
    fprintf(out, "%.10g\n%.10g\n%.10g\n", basis.exx, basis.eyx, basis.ezx);
    fprintf(out, "%.10g\n%.10g\n%.10g\n", basis.exx, basis.eyx, basis.ezx);

    for (int i = 0; i < size; i++) fprintf(out, "%.10g\n", val[i]);
    fclose(out);
  }

  virtual ~Grid() {
    delete[] val;
  }

  virtual void zero() {
    for (int i = 0; i < size; i++) {
      val[i] = 0.0;
    }
  }
  
  bool setValue(int j, double v) {
    if (j < 0 || j >= size) return false;
    val[j] = v;
    return true;
  }

  bool setValue(int ix, int iy, int iz, double v) {
    if (ix < 0 || ix >= nx) return false;
    if (iy < 0 || iy >= ny) return false;
    if (iz < 0 || iz >= nz) return false;
    int j = iz + iy*nz + ix*ny*nz;

    val[j] = v;
    return true;
  }

  double getValue(int j) const {
    if (j < 0 || j >= size) return 0.0;
    return val[j];
  }

  double getValue(int ix, int iy, int iz) const {
    if (ix < 0 || ix >= nx) return 0.0;
    if (iy < 0 || iy >= ny) return 0.0;
    if (iz < 0 || iz >= nz) return 0.0;
    
    int j = iz + iy*nz + ix*ny*nz;
    return val[j];
  }

  Vector3 getPosition(int ix, int iy, int iz) const {
    return basis.transform(Vector3(ix, iy, iz)) + origin;
  }

  Vector3 getPosition(int j) const {
    int iz = j%nz;
    int iy = (j/nz)%ny;
    int ix = j/(nz*ny);

    return basis.transform(Vector3(ix, iy, iz)) + origin;
  }

  Vector3 transformTo(Vector3 r) {
    return basisInv.transform(r-origin);
  }
  Vector3 transformFrom(Vector3 l) {
    return basis.transform(l) + origin;
  }

  int nearestIndex(Vector3 r) {
    Vector3 l = basisInv.transform(r-origin);
    
    int ix = int(floor(l.x + 0.5));
    int iy = int(floor(l.y + 0.5));
    int iz = int(floor(l.z + 0.5));

    ix = wrap(ix, nx);
    iy = wrap(iy, ny);
    iz = wrap(iz, nz);
    
    return iz + iy*nz + ix*ny*nz;
  }

  int length() const {
    return size;
  }
  Vector3 getOrigin() const {return origin;}
  Matrix3 getBasis() const {return basis;}
  Matrix3 getInverseBasis() const {return basisInv;}
  int getNx() const {return nx;}
  int getNy() const {return ny;}
  int getNz() const {return nz;}
  int getSize() const {return size;}
  void setBasis(const Matrix3& b) {
    basis = b;
    basisInv = basis.inverse();
  }
  void setOrigin(const Vector3& o) {
    origin = o;
  }
  Matrix3 getBox() const {
    return Matrix3(nx*basis.ex(), ny*basis.ey(), nz*basis.ez());
  } 
  Vector3 getExtent() const {
    return basis.transform(Vector3(nx,ny,nz));
  } 
  Vector3 getDestination() const {
    return basis.transform(Vector3(nx,ny,nz)) + origin;
  }
  Vector3 getCenter() const {
    return basis.transform(Vector3(nx/2,ny/2,nz/2)) + origin;
  }
  double getCellVolume() const {
    return fabs(basis.det());
  }
  double getVolume() const {
    return getCellVolume()*size;
  }
  Vector3 getCellDiagonal() const {
    return basis.ex() + basis.ey() + basis.ez();
  }

  void shift(double s) {
    for (int i = 0; i < size; i++) val[i] += s;
  }

  void scale(double s) {
    for (int i = 0; i < size; i++) val[i] *= s;
  }

  double max() const {
    double maxValue = val[0];
      for (int i = 1; i < size; i++) {
	if (val[i] > maxValue) maxValue = val[i];
      }
      return maxValue;
  }

  double min() const {
    double minValue = val[0];
      for (int i = 1; i < size; i++) {
	if (val[i] < minValue) minValue = val[i];
      }
      return minValue;
  }

  double mean() const {
    double sum = 0.0;
    for (int i = 0; i < size; i++) sum += val[i];
      return sum/size;
  }

  double meanBoltzmann() const {
    double sum = 0.0;
    double sumWeight = 0.0;
    for (int i = 0; i < size; i++) {
      double weight = exp(-val[i]);

      sumWeight += weight;
      sum += val[i]*weight;
    }
    
    return sum/sumWeight;
  }
  
  void shiftToCenters() {
    origin += basis.transform(Vector3(0.5));
  }

  // Add a constant gradient to the grid.
  void addGradient(Vector3 g) {
    int j = 0;

    for (int ix = 0; ix < nx; ix++) {
      for (int iy = 0; iy < ny; iy++) {
	for (int iz = 0; iz < nz; iz++) {
	  Vector3 dr = basis.transform(Vector3(ix,iy,iz));
	  double pot0 = dr.dot(g);
	  val[j] += pot0;
	  j++;
	}
      }
    }
  }

   // Tile the system cx by cy by cz.
  Grid tile(int cx, int cy, int cz) const {
    const int nynz = ny*nz;
    int nx1 = abs(cx*nx);
    int ny1 = abs(cy*ny);
    int nz1 = abs(cz*nz);
    Grid dest(basis, origin, nx1, ny1, nz1);
    
    const int dy = dest.nz;
    const int dx = dest.nz*dest.ny;
    

    for (int ix = 0; ix < dest.nx; ix++) {
      for (int iy = 0; iy < dest.ny; iy++) {
	for (int iz = 0; iz < dest.nz; iz++) {
	      int dj = iz + iy*dy + ix*dx;
	      int j = wrap(iz,nz) + wrap(iy,ny)*nz + wrap(ix,nx)*nynz;
	      dest.val[dj] = val[j];
	}
      }
    }
    return dest;
  }

  // Get a profile along z, at a position defined by factorX and factorY.
  // factorX = factorY = 0 means that the profile starts from the origin.
  // factorX = factorY = 0.5 means the profile is along the center of the grid.
  void profileZ(double factorX, double factorY, const char* fileName) {
    int ix = int(floor((factorX*nx + 0.5)));
    int iy = int(floor((factorY*ny + 0.5)));
    const int nynz = ny*nz;

    if (ix < 0) ix = 0;
    if (ix >= nx) ix = nx-1;
    if (iy < 0) iy = 0;
    if (iy >= ny) iy = ny-1;
    
    FILE* out = fopen(fileName,"w");
    if (out == NULL) {
      printf("Couldn't open file %s.\n",fileName);
      return;
    }

    for (int iz = 0; iz < nz; iz++) {
      int j = iz + iy*nz + ix*nynz;
      double v = val[j];
      double z = origin.z + iz*basis.ezz;
      fprintf(out, "%0.10g %0.10g\n", z, v);
    }
    fclose(out);
  }

  // Specify the real space x and y.
  void profileZ(double x, double y, double z0, double z1, const char* fileName) {
    int iz0 = int ((z0 - origin.z) / basis.ezz);
    int iz1 = int ((z1 - origin.z) / basis.ezz);

    if (iz0 < 0) iz0 = 0;
    if (iz1 >= nz) iz1 = nz-1;
    
    FILE* out = fopen(fileName,"w");
    if (out == NULL) {
      printf("Couldn't open file %s.\n",fileName);
      return;
    }

    Vector3 position;
    position.x = x;
    position.y = y;

    for (int iz = iz0; iz <= iz1; iz++) {

      double z = origin.z + iz*basis.ezz;
      position.z = z;

      fprintf(out, "%0.10g %0.10g\n", z, interpolatePotential(position));
    }
    fclose(out);
  }


  // Get the average value of a section of the grid.
  double averageSection(int iz) {
    const int nynz = ny*nz;
    double sum = 0.0;

    for (int ix = 0; ix < nx; ix++) {
      for (int iy = 0; iy < ny; iy++) {
	int j = iz + iy*nz + ix*nynz;
	sum += val[j];
      }
    }
      
    double v = sum/(nx*ny);
    return v;
  }

  // Get the average value of a section of the grid.
  double averageZLine(int ix, int iy, int z0, int z1) {
    const int nynz = ny*nz;
    double sum = 0.0;

    for (int iz = z0; iz <= z1; iz++) {
	int j = iz + iy*nz + ix*nynz;
	sum += val[j];
    }
      
    double v = sum/(z1-z0+1);
    return v;
  }

  // Will return the error in a variable
  double averageZLine(int ix, int iy, int z0, int z1, double& error, double& stddev) {
    const int nynz = ny*nz;
    double sum = 0.0;

    for (int iz = z0; iz <= z1; iz++) {
	int j = iz + iy*nz + ix*nynz;
	sum += val[j];
    }
      
    double v = sum/(z1-z0+1);

    // Now calculate the errors
    sum = 0.0;
    for (int iz = z0; iz <= z1; iz++) {
	int j = iz + iy*nz + ix*nynz;
	sum += (v-val[j])*(v-val[j]);
    }

    stddev = sqrt(sum/(z1-z0)); // n-1
    error = 1.96 *sqrt(sum)/(z1-z0+1);

    return v;
  }

  // Get the average value of a section of the grid.
  // Similar to average region, but takes direct grid point inputs
  double averageZBox(int x0, int x1, int y0, int y1, int z0, int z1) {
    const int nynz = ny*nz;
    double sum = 0.0;

    for (int ix = x0; ix <= x1; ix++) {
      for (int iy = y0; iy <= y1; iy++) {
	for (int iz = z0; iz <= z1; iz++) {
	  int j = iz + iy*nz + ix*nynz;
	  sum += val[j];
	}
      }
    }
    
    double v = sum/(z1-z0+1+x1-x0+1+y1-y0+1);
    return v;
  }


  // Compute the average profile along z.
  void averageProfileZ(const char* fileName) {
    const int nynz = ny*nz;
    
    FILE* out = fopen(fileName,"w");
    if (out == NULL) {
      printf("Couldn't open file %s.\n",fileName);
      exit(-1);
    }

   
    for (int iz = 0; iz < nz; iz++) {
      double sum = 0.0;

      for (int ix = 0; ix < nx; ix++) {
	for (int iy = 0; iy < ny; iy++) {
	  int j = iz + iy*nz + ix*nynz;
	  sum += val[j];
	}
      }
      
      double v = sum/(nx*ny);
      double z = origin.z + iz*basis.ezz;
      fprintf(out, "%0.10g %0.10g\n", z, v);
    }

    fclose(out);
  }

  // Compute the average profile along an axis.
  // Assumes that the grid axis with index "axis" is aligned with the world axis of index "axis".
  void averageProfile(const char* fileName, int axis) {
    const int nynz = ny*nz;
    
    FILE* out = fopen(fileName,"w");
    if (out == NULL) {
      printf("Couldn't open file %s.\n",fileName);
      exit(-1);
    }

    int dir0 = wrap(axis, 3);
    int dir1 = (axis+1)%3;
    int dir2 = (axis+2)%3;

    int jump[3];
    jump[0] = nynz;
    jump[1] = nz;
    jump[2] = 1;

    int n[3];
    n[0] = nx;
    n[1] = ny;
    n[2] = nz;
   
    for (int i0 = 0; i0 < n[dir0]; i0++) {
      double sum = 0;

      for (int i1 = 0; i1 < n[dir1]; i1++) {
	for (int i2 = 0; i2 < n[dir2]; i2++) {
	  int j = i0*jump[dir0] + i1*jump[dir1] + i2*jump[dir2];
	  sum += val[j];
	}
      }
      
      double v = sum/(n[dir1]*n[dir2]);
      double x = 0.0;
      switch (dir0) {
      case 0:
	x = origin.x + i0*basis.exx;
	break;
      case 1:
	x = origin.y + i0*basis.eyy;
	break;
      case 2:
	x = origin.z + i0*basis.ezz;
	break;
      }
      fprintf(out, "%0.10g %0.10g\n", x, v);
    }

    fclose(out);
  }

  // Compute the average profile along z.
  void averageProfileZBoltzmann(const char* fileName) {
    const int nynz = ny*nz;
    
    FILE* out = fopen(fileName,"w");
    if (out == NULL) {
      printf("Couldn't open file %s.\n",fileName);
      exit(-1);
    }

   
    for (int iz = 0; iz < nz; iz++) {
      double sumNum = 0.0;
      double sumDen = 0.0;

      for (int ix = 0; ix < nx; ix++) {
	for (int iy = 0; iy < ny; iy++) {
	  int j = iz + iy*nz + ix*nynz;
	  double weight = exp(-val[j]);

	  sumNum += val[j]*weight;
	  sumDen += weight;
	}
      }
      
      double v = sumNum/sumDen;
      double z = origin.z + iz*basis.ezz;
      fprintf(out, "%0.10g %0.10g\n", z, v);
    }

    fclose(out);
  }

  void reflectZ() {
    const int nynz = ny*nz;
    const int cz = nz/2;

    for (int ix = 0; ix < nx; ix++) {
      for (int iy = 0; iy < ny; iy++) {
	for (int iz = 0; iz < cz; iz++) {
	  int j = iz + iy*nz + ix*nynz;
	  int k = (nz-iz-1) + iy*nz + ix*nynz;

	  val[k] = val[j];
	}
      }
    }
  }


  // Create profile along z averaged over a cylinder.
  void cylinderProfileZ(const char* fileName, double radius) {
    const int nynz = ny*nz;
    const double rad2 = radius*radius;

    FILE* out = fopen(fileName,"w");
    if (out == NULL) {
      printf("Couldn't open file %s.\n",fileName);
      exit(-1);
    }

    for (int iz = 0; iz < nz; iz++) {
      double sum = 0.0;

      int n = 0;
      for (int ix = 0; ix < nx; ix++) {
	for (int iy = 0; iy < ny; iy++) {
	  Vector3 r = basis.transform(Vector3(ix,iy,iz)) + origin;

	  if (r.x*r.x + r.y*r.y < rad2) {
	    int j = iz + iy*nz + ix*nynz;
	    sum += val[j];
	    n++;
	  }
	}
      }
      
      double v = sum/n;
      double z = origin.z + iz*basis.ezz;
      fprintf(out, "%0.10g %0.10g\n", z, v);
    }

    fclose(out);
  }

  // Compute the average value over a box.
  double averageRegion(Vector3 r0, Vector3 r1) {
    int i = 0;
    int n = 0;
    double sum = 0.0;

    for (int ix = 0; ix < nx; ix++) {
      for (int iy = 0; iy < ny; iy++) {
	for (int iz = 0; iz < nz; iz++) {
	  Vector3 r = basis.transform(Vector3(ix,iy,iz)) + origin;
	  if (r.x >= r0.x && r.x < r1.x && r.y >= r0.y && r.y < r1.y 
	      && r.z >= r0.z && r.z < r1.z) {
	    sum += val[i];
	    n++;
	  }
	  
	  i++;
	}
      }
    }

    return sum/n;
  }

  // Add two grids. The resulting grids has the dimensions of *this.
  // The added values come from interpolation.
  void addInterpolate(const Grid& g) {
    Vector3 d[27];

    // Make the images of the g in this's basis.
    int j = 0;
    for (int ix = -1; ix <= 1; ix++) {
      for (int iy = -1; iy <= 1; iy++) {
	for (int iz = -1; iz <= 1; iz++) {
	  d[j] = basis.transform(Vector3(ix*nx,iy*ny,iz*nz));
	  j++;
	}
      }
    }

    // Check each point for a value to add.
    int i = 0;
    for (int ix = 0; ix < nx; ix++) {
      for (int iy = 0; iy < ny; iy++) {
	for (int iz = 0; iz < nz; iz++) {
	  Vector3 r = basis.transform(Vector3(ix,iy,iz)) + origin;
	  
	  // Paste only one image.
	  for (int k = 0; k < 27; k++) {
	    Vector3 rp = r + d[k];
	    
	    Vector3 gr = g.basisInv.transform(rp-g.origin);
	    // Paste only one image.
	    if (gr.x < 0.0 || gr.x >= g.nx) continue;
	    if (gr.y < 0.0 || gr.y >= g.ny) continue;
	    if (gr.z < 0.0 || gr.z >= g.nz) continue;
	    
	    // Add the interpolated value.
	    val[i] += g.interpolatePotential(rp);
	    // Add it only once.
	    break;
	  }

	  i++;
	}
      }
    }
  }

  // Add two grids with the same number of nodes.
  bool add(const Grid& g) {
    if (size != g.size) return false;

    for (int i = 0; i < size; i++) val[i] += g.val[i];
    return true;
  }

  // Add two grids. The resulting grids has the dimensions of *this.
  // The added values come from the spatially nearest nodes from g.
  void addGrid(const Grid& g) {
    int i = 0;
    const int nynz = ny*nz;

     for (int ix = 0; ix < g.nx; ix++) {
      for (int iy = 0; iy < g.ny; iy++) {
	for (int iz = 0; iz < g.nz; iz++) {
	  Vector3 r = g.basis.transform(Vector3(ix,iy,iz)) + g.origin;
	  Vector3 l = basisInv.transform(r-origin);
	  
	  // Find the nearest node in *this.
	  int jx = wrap(int(floor(l.x + 0.5)), nx);
	  int jy = wrap(int(floor(l.y + 0.5)), ny);
	  int jz = wrap(int(floor(l.z + 0.5)), nz);
	  int j = jz + jy*nz + jx*nynz;

	  val[j] += g.val[i];
	  i++;
	}
      }
    }
  }

  // Add two grids. The resulting grids has the dimensions of *this.
  // The added values are linearly distributed.
  void addLinear(const Grid& g) {
    const int nynz = ny*nz;

    // Loop through all nodes of the grid g.
    for (int i = 0; i < g.size; i++) {
      Vector3 r = g.getPosition(i);
      Vector3 l = basisInv.transform(r-origin);

      // Find the nearest node in *this.
      int homeX = int(floor(l.x));
      int homeY = int(floor(l.y));
      int homeZ = int(floor(l.z));
      double wx = l.x - homeX;
      double wy = l.y - homeY;
      double wz = l.z - homeZ;

      // Find the eight neighboring nodes and distribute the value among them.
      for (int jx = 0; jx <= 1; jx++) {
	for (int jy = 0; jy <= 1; jy++) {
	  for (int jz = 0; jz <= 1; jz++) {
	    int kx = wrap(homeX + jx, nx);
	    int ky = wrap(homeY + jy, ny);
	    int kz = wrap(homeZ + jz, nz);

	    int ind = kz + ky*nz + kx*nynz;; 
	    double w = ((jx==0)?(1-wx):wx)*((jy==0)?(1-wy):wy)*((jz==0)?(1-wz):wz);
	    val[ind] += w*g.val[i];
	  }
	}
      }
    }
  }

  void crossSectionFactor(const char* fileName, int dir, double factor) const {
    int dir0 = wrap(dir, 3);
    int dim[3];
    dim[0] = nx;
    dim[1] = ny;
    dim[2] = nz;
    
    int section = int(floor(factor*dim[dir0] + 0.5));
    crossSection(fileName, dir, section);
  }

  void crossSection(const char* fileName, int dir, int section) const {
    // Open the potential output file.
    FILE* out = fopen(fileName, "w");
    if (out == NULL) {
      printf("Grid:crossSection Couldn't open file %s.\n", fileName);
      exit(-1);
    }
    

    // Open the coordinate output files.
    char s[256];
    sprintf(s, "%s.rx", fileName);
    FILE* outX = fopen(s, "w");
    sprintf(s, "%s.ry", fileName);
    FILE* outY = fopen(s, "w");
    sprintf(s, "%s.rz", fileName);
    FILE* outZ = fopen(s, "w");
    
    // Choose the directions.
    int dir0 = wrap(dir, 3);
    int dir1 = (dir0 + 1) % 3;
    int dir2 = (dir0 + 2) % 3;
    
    int dim[3];
    dim[0] = nx;
    dim[1] = ny;
    dim[2] = nz;

    int jump[3];
    jump[0] = ny*nz;
    jump[1] = nz;
    jump[2] = 1;

    // Write the cross section.
    int i0 = wrap(section, dim[dir0]);
    for (int i1 = 0; i1 < dim[dir1]; i1++) {
      int k = i0*jump[dir0] + i1*jump[dir1];
      Vector3 r = getPosition(k);
      fprintf(out, "%.10g", val[k]);
      fprintf(outX, "%.10g", r.x);
      fprintf(outY, "%.10g", r.y);
      fprintf(outZ, "%.10g", r.z);

      for (int i2 = 1; i2 < dim[dir2]; i2++) {
	int j = i0*jump[dir0] + i1*jump[dir1] + i2*jump[dir2];
	Vector3 r1 = getPosition(j);
	fprintf(out, " %.10g", val[j]);
	fprintf(outX, " %.10g", r1.x);
	fprintf(outY, " %.10g", r1.y);
	fprintf(outZ, " %.10g", r1.z);
      }
      fprintf(out, "\n");
      fprintf(outX, "\n");
      fprintf(outY, "\n");
      fprintf(outZ, "\n");
    }

    fclose(out);
    fclose(outX);
    fclose(outY);
    fclose(outZ);
  }
  
  bool crop(int ix0, int iy0, int iz0, int ix1, int iy1, int iz1) {
    const int nynz = ny*nz;
    if (ix0 < 0 || ix0 >= nx) ix0 = 0;
    if (ix1 < 0 || ix1 >= nx) ix1 = nx-1;
    if (iy0 < 0 || iy0 >= ny) iy0 = 0;
    if (iy1 < 0 || iy1 >= ny) iy1 = ny-1;
    if (iz0 < 0 || iz0 >= nz) iz0 = 0;
    if (iz1 < 0 || iz1 >= nz) iz1 = nz-1;

    if (ix0 >= ix1) return false;
    if (iy0 >= iy1) return false;
    if (iz0 >= iz1) return false;

    int newNx = ix1-ix0+1;
    int newNy = iy1-iy0+1;
    int newNz = iz1-iz0+1;
    int newSize = newNx*newNy*newNz;
    double* v = new double[newSize];

    // Copy the appropriate data.
    int i = 0;
    for (int ix = ix0; ix <= ix1; ix++) {
      for (int iy = iy0; iy <= iy1; iy++) {
	for (int iz = iz0; iz <= iz1; iz++) {
	  int j = iz + iy*nz + ix*nynz;
	  v[i] = val[j];
	  i++;
	}
      }
    }

    // Determine the new origin and set the new members.
    origin += basis.transform(Vector3(ix0,iy0,iz0));
    nx = newNx;
    ny = newNy;
    nz = newNz;
    size = newSize;

    // Swap the pointers and deallocate the old array.
    double* valOld = val;
    val = v;
    delete[] valOld;
    
    
    return true;
  }

  // Remove parts of the grid outside of the box defined by
  // r0 and r1.
  void crop(Vector3 r0, Vector3 r1) {
    // Find the start and end along x.
    int ix0 = -1;
    int ix1 = -1;
    for (int ix = 0; ix < nx; ix++) {
      Vector3 r = basis.transform(Vector3(ix,0,0)) + origin;
      if (ix0 < 0 && r.x >= r0.x) ix0 = ix;
      if (ix1 < 0 && r.x > r1.x) {
	ix1 = ix-1;
	break;
      }
    }
 
    // Find the start and end along y.
    int iy0 = -1;
    int iy1 = -1;
    for (int iy = 0; iy < ny; iy++) {
      Vector3 r = basis.transform(Vector3(0,iy,0)) + origin;
      if (iy0 < 0 && r.y >= r0.y) iy0 = iy;
      if (iy1 < 0 && r.y > r1.y) {
	iy1 = iy-1;
	break;
      }
    }
 
    // Find the start and end along z.
    int iz0 = -1;
    int iz1 = -1;
    for (int iz = 0; iz < nz; iz++) {
      Vector3 r = basis.transform(Vector3(0,0,iz)) + origin;
      if (iz0 < 0 && r.z >= r0.z) iz0 = iz;
      if (iz1 < 0 && r.z > r1.z) {
	iz1 = iz-1;
	break;
      }
    }

    // Use the index-based cropping function.
    crop(ix0, iy0, iz0, ix1, iy1, iz1);
  }

  // Resample this grid.
  Grid resample(int mx, int my, int mz) const {
    Vector3 ex = (nx*basis.ex())/mx;
    Vector3 ey = (ny*basis.ey())/my;
    Vector3 ez = (nz*basis.ez())/mz;
    Matrix3 newBasis = Matrix3(ex, ey, ez);
    
    Grid ret(newBasis, origin, mx, my, mz);
    int i = 0;
    for (int ix = 0; ix < ret.nx; ix++) {
      for (int iy = 0; iy < ret.ny; iy++) {
	for (int iz = 0; iz < ret.nz; iz++) {
	  Vector3 r = ret.basis.transform(Vector3(ix,iy,iz)) + ret.origin;
	  ret.val[i] = interpolatePotential(r);
	  i++;
	}
      }
    }

    return ret;
  }

  void paste(const Grid& clip) {
    for (int ix = 0; ix < nx; ix++) {
      for (int iy = 0; iy < ny; iy++) {
	for (int iz = 0; iz < nz; iz++) {
	  int j = iz + iy*nz + ix*nz*ny;

	  // Find each node in this grid.
	  Vector3 r = basis.transform(Vector3(ix,iy,iz)) + origin;
	  Vector3 clipR = clip.basisInv.transform(r-clip.origin);

	  // Check that this point is in the clipboard grid.
	  if (clipR.x < 0.0 || clipR.x >= clip.nx) continue;
	  if (clipR.y < 0.0 || clipR.y >= clip.ny) continue;
	  if (clipR.z < 0.0 || clipR.z >= clip.nz) continue;
	  
	  double pot = clip.interpolatePotential(r);
	  val[j] = pot;
	}
      }
    }
  }

  // Make each node the average value of its neighbors.
  void blur() {
    const int nzny = nz*ny;
    double* v = new double[size];

    for (int ix = 0; ix < nx; ix++) {
      for (int iy = 0; iy < ny; iy++) {
	for (int iz = 0; iz < nz; iz++) {
	  int j = iz + iy*nz + ix*nzny;

	  int ix0 = wrap(ix-1,nx);
	  int ix1 = wrap(ix+1,nx);
	  int iy0 = wrap(iy-1,ny);
	  int iy1 = wrap(iy+1,ny);
	  int iz0 = wrap(iz-1,nz);
	  int iz1 = wrap(iz+1,nz);

	  int jx0 = iz + iy*nz + ix0*nzny;
	  int jx1 = iz + iy*nz + ix1*nzny;
	  int jy0 = iz + iy0*nz + ix*nzny;
	  int jy1 = iz + iy1*nz + ix*nzny;
	  int jz0 = iz0 + iy*nz + ix*nzny;
	  int jz1 = iz1 + iy*nz + ix*nzny;

	  v[j] = (val[jx0]+val[jx1]+val[jy0]+val[jy1]+val[jz0]+val[jz1])/6.0;
	}
      }
    }

    // Swap the pointers.
    delete[] val;
    val = v;
  }

  // Make a map of the minimum value along z.
  // Format is "x y z valMin".
  void depthMapZ(const char* fileName) {
    FILE* out = fopen(fileName, "w");
    if (out == NULL) {
      printf("Grid:depthMapZ Couldn't open file %s.\n",fileName);
      exit(-1);
    }
    
    char profileName[256];
    sprintf(profileName, "%s.profile", fileName);
    FILE* outProfile = fopen(profileName, "w");
    if (outProfile == NULL) {
      printf("Grid:depthMapZ Couldn't open file %s.\n",fileName);
      exit(-1);
    }

    int i = 0;
    // Loop over x and y.
    for (int ix = 0; ix < nx; ix++) {
      for (int iy = 0; iy < ny; iy++) {
	
	// Find the minimum value along z.
	int minIz = i;
	double minVal = val[i];
	for (int iz = 0; iz < nz; iz++) {
	  fprintf(outProfile, "%.10g ", val[i]);

	  if (val[i] < minVal) {
	    minVal = val[i];
	    minIz = iz;
	  }
	  i++;
	}

	// Print the map of the minimum values.
	Vector3 r = basis.transform(Vector3(ix,iy,minIz)) + origin;
	fprintf(out, "%.10g %.10g %.10g %.10g\n", r.x, r.y, r.z, minVal);
	fprintf(outProfile, "\n");
      }
    }
    fclose(outProfile);
    fclose(out);
  }

  double getPotential(Vector3 pos) const {
    // Find the nearest node.
    Vector3 l = basisInv.transform(pos - origin);
    int ix = wrap(int(floor(l.x + 0.5)), nx);
    int iy = wrap(int(floor(l.y + 0.5)), ny);
    int iz = wrap(int(floor(l.z + 0.5)), nz);
    
    int j = iz + iy*nz + ix*ny*nz;
    return val[j];
  }

  double interpolatePotential(Vector3 pos) const {
    // Find the home node.
    Vector3 l = basisInv.transform(pos - origin);
    int homeX = int(floor(l.x));
    int homeY = int(floor(l.y));
    int homeZ = int(floor(l.z));
    
    // Get the array jumps.
    int jump[3];
    jump[0] = nz*ny;
    jump[1] = nz;
    jump[2] = 1;

    // Shift the indices in the home array.
    int home[3];
    home[0] = homeX;
    home[1] = homeY;
    home[2] = homeZ;

    // Get the grid dimensions.
    int g[3];
    g[0] = nx;
    g[1] = ny;
    g[2] = nz;

    // Get the interpolation coordinates.
    double w[3];
    w[0] = l.x - homeX;
    w[1] = l.y - homeY;
    w[2] = l.z - homeZ;

    // Find the values at the neighbors.
    double g1[4][4][4];
    for (int ix = 0; ix < 4; ix++) {
      for (int iy = 0; iy < 4; iy++) {
	for (int iz = 0; iz < 4; iz++) {
	  // Wrap around the periodic boundaries. 
	  int jx = ix-1 + home[0];
	  jx = wrap(jx, g[0]);
	  int jy = iy-1 + home[1];
	  jy = wrap(jy, g[1]);
	  int jz = iz-1 + home[2];
	  jz = wrap(jz, g[2]);
	  
	  int ind = jz*jump[2] + jy*jump[1] + jx*jump[0];
	  g1[ix][iy][iz] = val[ind];
	}
      }
    }
    double a0, a1, a2, a3;
  
    // Mix along x.
    double g2[4][4];
    for (int iy = 0; iy < 4; iy++) {
      for (int iz = 0; iz < 4; iz++) {

	a3 = 0.5*(-g1[0][iy][iz] + 3*g1[1][iy][iz] - 3*g1[2][iy][iz] + g1[3][iy][iz]);
	a2 = 0.5*(2*g1[0][iy][iz] - 5*g1[1][iy][iz] + 4*g1[2][iy][iz] - g1[3][iy][iz]);
	a1 = 0.5*(-g1[0][iy][iz] + g1[2][iy][iz]);
	a0 = g1[1][iy][iz];

	g2[iy][iz] = a3*w[0]*w[0]*w[0] + a2*w[0]*w[0] + a1*w[0] + a0;
      }
    }

    // Mix along y.
    double g3[4];
    for (int iz = 0; iz < 4; iz++) {
      a3 = 0.5*(-g2[0][iz] + 3*g2[1][iz] - 3*g2[2][iz] + g2[3][iz]);
      a2 = 0.5*(2*g2[0][iz] - 5*g2[1][iz] + 4*g2[2][iz] - g2[3][iz]);
      a1 = 0.5*(-g2[0][iz] + g2[2][iz]);
      a0 = g2[1][iz];
   
      g3[iz] = a3*w[1]*w[1]*w[1] + a2*w[1]*w[1] + a1*w[1] + a0;
    }

    // Mix along z.
    a3 = 0.5*(-g3[0] + 3*g3[1] - 3*g3[2] + g3[3]);
    a2 = 0.5*(2*g3[0] - 5*g3[1] + 4*g3[2] - g3[3]);
    a1 = 0.5*(-g3[0] + g3[2]);
    a0 = g3[1];

    return a3*w[2]*w[2]*w[2] + a2*w[2]*w[2] + a1*w[2] + a0;
  }

  double interpolateDiffX(Vector3 pos) const {
    // Find the home node.
    Vector3 l = basisInv.transform(pos - origin);
    int homeX = int(floor(l.x));
    int homeY = int(floor(l.y));
    int homeZ = int(floor(l.z));
    
    // Get the array jumps with shifted indices.
    int jump[3];
    jump[0] = nz*ny;
    jump[1] = nz;
    jump[2] = 1;

    // Shift the indices in the home array.
    int home[3];
    home[0] = homeX;
    home[1] = homeY;
    home[2] = homeZ;

    // Shift the indices in the grid dimensions.
    int g[3];
    g[0] = nx;
    g[1] = ny;
    g[2] = nz;

    // Get the interpolation coordinates.
    double w[3];
    w[0] = l.x - homeX;
    w[1] = l.y - homeY;
    w[2] = l.z - homeZ;

    // Find the values at the neighbors.
    double g1[4][4][4];
    for (int ix = 0; ix < 4; ix++) {
      for (int iy = 0; iy < 4; iy++) {
	for (int iz = 0; iz < 4; iz++) {
	  // Wrap around the periodic boundaries. 
	  int jx = ix-1 + home[0];
	  jx = wrap(jx, g[0]);
	  int jy = iy-1 + home[1];
	  jy = wrap(jy, g[1]);
	  int jz = iz-1 + home[2];
	  jz = wrap(jz, g[2]);
	  
	  int ind = jz*jump[2] + jy*jump[1] + jx*jump[0];
	  g1[ix][iy][iz] = val[ind];
	}
      }
    }
    double a0, a1, a2, a3;
  
    // Mix along x, taking the derivative.
    double g2[4][4];
    for (int iy = 0; iy < 4; iy++) {
      for (int iz = 0; iz < 4; iz++) {

	a3 = 0.5*(-g1[0][iy][iz] + 3*g1[1][iy][iz] - 3*g1[2][iy][iz] + g1[3][iy][iz]);
	a2 = 0.5*(2*g1[0][iy][iz] - 5*g1[1][iy][iz] + 4*g1[2][iy][iz] - g1[3][iy][iz]);
	a1 = 0.5*(-g1[0][iy][iz] + g1[2][iy][iz]);
	a0 = g1[1][iy][iz];

	//g2[iy][iz] = a3*w[0]*w[0]*w[0] + a2*w[0]*w[0] + a1*w[0] + a0;
	g2[iy][iz] = 3.0*a3*w[0]*w[0] + 2.0*a2*w[0] + a1;
      }
    }

    // Mix along y.
    double g3[4];
    for (int iz = 0; iz < 4; iz++) {
      a3 = 0.5*(-g2[0][iz] + 3*g2[1][iz] - 3*g2[2][iz] + g2[3][iz]);
      a2 = 0.5*(2*g2[0][iz] - 5*g2[1][iz] + 4*g2[2][iz] - g2[3][iz]);
      a1 = 0.5*(-g2[0][iz] + g2[2][iz]);
      a0 = g2[1][iz];
   
      g3[iz] = a3*w[1]*w[1]*w[1] + a2*w[1]*w[1] + a1*w[1] + a0;
    }

    // Mix along z.
    a3 = 0.5*(-g3[0] + 3*g3[1] - 3*g3[2] + g3[3]);
    a2 = 0.5*(2*g3[0] - 5*g3[1] + 4*g3[2] - g3[3]);
    a1 = 0.5*(-g3[0] + g3[2]);
    a0 = g3[1];

    return -(a3*w[2]*w[2]*w[2] + a2*w[2]*w[2] + a1*w[2] + a0);
  }

  double interpolateDiffY(Vector3 pos) const {
    // Find the home node.
    Vector3 l = basisInv.transform(pos - origin);
    int homeX = int(floor(l.x));
    int homeY = int(floor(l.y));
    int homeZ = int(floor(l.z));
    
    // Get the array jumps with shifted indices.
    int jump[3];
    jump[0] = nz*ny;
    jump[1] = nz;
    jump[2] = 1;

    // Shift the indices in the home array.
    int home[3];
    home[0] = homeX;
    home[1] = homeY;
    home[2] = homeZ;

    // Shift the indices in the grid dimensions.
    int g[3];
    g[0] = nx;
    g[1] = ny;
    g[2] = nz;

    // Get the interpolation coordinates.
    double w[3];
    w[0] = l.x - homeX;
    w[1] = l.y - homeY;
    w[2] = l.z - homeZ;

    // Find the values at the neighbors.
    double g1[4][4][4];
    for (int ix = 0; ix < 4; ix++) {
      for (int iy = 0; iy < 4; iy++) {
	for (int iz = 0; iz < 4; iz++) {
	  // Wrap around the periodic boundaries. 
	  int jx = ix-1 + home[0];
	  jx = wrap(jx, g[0]);
	  int jy = iy-1 + home[1];
	  jy = wrap(jy, g[1]);
	  int jz = iz-1 + home[2];
	  jz = wrap(jz, g[2]);
	  
	  int ind = jz*jump[2] + jy*jump[1] + jx*jump[0];
	  g1[ix][iy][iz] = val[ind];
	}
      }
    }
    double a0, a1, a2, a3;
  
    // Mix along x, taking the derivative.
    double g2[4][4];
    for (int iy = 0; iy < 4; iy++) {
      for (int iz = 0; iz < 4; iz++) {

	a3 = 0.5*(-g1[0][iy][iz] + 3*g1[1][iy][iz] - 3*g1[2][iy][iz] + g1[3][iy][iz]);
	a2 = 0.5*(2*g1[0][iy][iz] - 5*g1[1][iy][iz] + 4*g1[2][iy][iz] - g1[3][iy][iz]);
	a1 = 0.5*(-g1[0][iy][iz] + g1[2][iy][iz]);
	a0 = g1[1][iy][iz];

	g2[iy][iz] = a3*w[0]*w[0]*w[0] + a2*w[0]*w[0] + a1*w[0] + a0;
      }
    }

    // Mix along y.
    double g3[4];
    for (int iz = 0; iz < 4; iz++) {
      a3 = 0.5*(-g2[0][iz] + 3*g2[1][iz] - 3*g2[2][iz] + g2[3][iz]);
      a2 = 0.5*(2*g2[0][iz] - 5*g2[1][iz] + 4*g2[2][iz] - g2[3][iz]);
      a1 = 0.5*(-g2[0][iz] + g2[2][iz]);
      a0 = g2[1][iz];
   
      //g3[iz] = a3*w[1]*w[1]*w[1] + a2*w[1]*w[1] + a1*w[1] + a0;
      g3[iz] = 3.0*a3*w[1]*w[1] + 2.0*a2*w[1] + a1;
    }

    // Mix along z.
    a3 = 0.5*(-g3[0] + 3*g3[1] - 3*g3[2] + g3[3]);
    a2 = 0.5*(2*g3[0] - 5*g3[1] + 4*g3[2] - g3[3]);
    a1 = 0.5*(-g3[0] + g3[2]);
    a0 = g3[1];

    return -(a3*w[2]*w[2]*w[2] + a2*w[2]*w[2] + a1*w[2] + a0);
  }

  double interpolateDiffZ(Vector3 pos) const {
    // Find the home node.
    Vector3 l = basisInv.transform(pos - origin);
    int homeX = int(floor(l.x));
    int homeY = int(floor(l.y));
    int homeZ = int(floor(l.z));
    
    // Get the array jumps with shifted indices.
    int jump[3];
    jump[0] = nz*ny;
    jump[1] = nz;
    jump[2] = 1;

    // Shift the indices in the home array.
    int home[3];
    home[0] = homeX;
    home[1] = homeY;
    home[2] = homeZ;

    // Shift the indices in the grid dimensions.
    int g[3];
    g[0] = nx;
    g[1] = ny;
    g[2] = nz;

    // Get the interpolation coordinates.
    double w[3];
    w[0] = l.x - homeX;
    w[1] = l.y - homeY;
    w[2] = l.z - homeZ;

    // Find the values at the neighbors.
    double g1[4][4][4];
    for (int ix = 0; ix < 4; ix++) {
      for (int iy = 0; iy < 4; iy++) {
	for (int iz = 0; iz < 4; iz++) {
	  // Wrap around the periodic boundaries. 
	  int jx = ix-1 + home[0];
	  jx = wrap(jx, g[0]);
	  int jy = iy-1 + home[1];
	  jy = wrap(jy, g[1]);
	  int jz = iz-1 + home[2];
	  jz = wrap(jz, g[2]);
	  
	  int ind = jz*jump[2] + jy*jump[1] + jx*jump[0];
	  g1[ix][iy][iz] = val[ind];
	}
      }
    }
    double a0, a1, a2, a3;
  
    // Mix along x, taking the derivative.
    double g2[4][4];
    for (int iy = 0; iy < 4; iy++) {
      for (int iz = 0; iz < 4; iz++) {

	a3 = 0.5*(-g1[0][iy][iz] + 3*g1[1][iy][iz] - 3*g1[2][iy][iz] + g1[3][iy][iz]);
	a2 = 0.5*(2*g1[0][iy][iz] - 5*g1[1][iy][iz] + 4*g1[2][iy][iz] - g1[3][iy][iz]);
	a1 = 0.5*(-g1[0][iy][iz] + g1[2][iy][iz]);
	a0 = g1[1][iy][iz];

	g2[iy][iz] = a3*w[0]*w[0]*w[0] + a2*w[0]*w[0] + a1*w[0] + a0;
      }
    }

    // Mix along y.
    double g3[4];
    for (int iz = 0; iz < 4; iz++) {
      a3 = 0.5*(-g2[0][iz] + 3*g2[1][iz] - 3*g2[2][iz] + g2[3][iz]);
      a2 = 0.5*(2*g2[0][iz] - 5*g2[1][iz] + 4*g2[2][iz] - g2[3][iz]);
      a1 = 0.5*(-g2[0][iz] + g2[2][iz]);
      a0 = g2[1][iz];
   
      g3[iz] = a3*w[1]*w[1]*w[1] + a2*w[1]*w[1] + a1*w[1] + a0;
    }

    // Mix along z.
    a3 = 0.5*(-g3[0] + 3*g3[1] - 3*g3[2] + g3[3]);
    a2 = 0.5*(2*g3[0] - 5*g3[1] + 4*g3[2] - g3[3]);
    a1 = 0.5*(-g3[0] + g3[2]);
    a0 = g3[1];

    return -(3.0*a3*w[2]*w[2] + 2.0*a2*w[2] + a1);
  }

  Vector3 interpolateForce(Vector3 pos) const {
    Vector3 f;
    f.x = interpolateDiffX(pos);
    f.y = interpolateDiffY(pos);
    f.z = interpolateDiffZ(pos);
    Vector3 f1 = basisInv.transpose().transform(f);
    return f1;
  }

  static double wrapFloat(double x, double l) {
    int image = int(floor(x/l));
    x -= image*l;
    return x;
  }

  static double wrapDiff(double x, double l) {
    int image = int(floor(x/l));
    x -= image*l;
    if (x >= 0.5*l) x -= l;
    return x;
  }

  static int wrap(int i, int n) {
    if (i < 0) {
      i %= n;
      i += n;
    }
    // The portion above allows i == n, so no else keyword.
    if (i >= n) i %= n;

    return i;
  }

  Vector3 wrap(Vector3 r) const {
    Vector3 l = basisInv.transform(r-origin);
    l.x = wrapFloat(l.x, nx);
    l.y = wrapFloat(l.y, ny);
    l.z = wrapFloat(l.z, nz);
    return basis.transform(l) + origin;
  }

  Vector3 wrapDiff(Vector3 r) const {
    Vector3 l = basisInv.transform(r);
    l.x = wrapDiff(l.x, nx);
    l.y = wrapDiff(l.y, ny);
    l.z = wrapDiff(l.z, nz);
    return basis.transform(l);
  }

protected:
  Vector3 origin;
  Matrix3 basis;
  int nx, ny, nz;
  int size;
  Matrix3 basisInv;
  double* val;
  Grid() {}

private:
  //Grid(const Grid&) {}
  //Grid operator=(const Grid&) {}
};
#endif
