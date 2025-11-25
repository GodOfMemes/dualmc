// Copyright (C) 2017, Dominik Wodniok
// This software may be modified and distributed under the terms
// of the BSD 3-Clause license.
// See the LICENSE.txt file for details.

#ifndef DUALMC_H_INCLUDED
#define DUALMC_H_INCLUDED

/// \file   dualmc.h
/// \author Dominik Wodniok
/// \date   2009

// c includes
#include <cstdint>

// stl includes
#include <unordered_map>
#include <vector>

namespace dualmc {
    

typedef float VertexComponentsType;
typedef uint32_t QuadIndexType;
typedef uint32_t TriIndexType;

/// vertex structure for dual points
struct Vertex {
    /// non-initializing constructor
    Vertex();

    /// initializing constructor
    Vertex(VertexComponentsType x, VertexComponentsType y, VertexComponentsType z);
    
    /// initializing constructor
    Vertex(Vertex const & v);
    
    // components
    VertexComponentsType x,y,z;
};

/// quad indices structure
struct Quad {
    /// non-initializing constructor
    Quad();
    
    /// initializing constructor
    Quad(QuadIndexType i0, QuadIndexType i1,QuadIndexType i2, QuadIndexType i3);
    
    // quad indices
    QuadIndexType i0,i1,i2,i3;
    
};

/// tri indices structure
struct Tri {
  /// non-initializing constructor
  Tri();

  /// initializing constructor
  Tri(TriIndexType i0, TriIndexType i1, TriIndexType i2);

  // tri indices
  TriIndexType i0, i1, i2;

};

/// \class  DualMC
/// \author Dominik Wodniok
/// \date   2009
/// Class which implements the dual marching cubes algorithm from Gregory M. Nielson.
/// Faces and vertices of the standard marching cubes algorithm correspond to
/// vertices and faces in the dual algorithm. As a vertex in standard marching cubes
/// usually is shared by 4 faces, the dual mesh is entirely made from quadrangles.
/// Unfortunately, under rare circumstances the original algorithm can create
/// non-manifold meshes. See the remarks of the original paper on this.
/// The class optionally can guarantee manifold meshes by taking the Manifold
/// Dual Marching Cubes approach from Rephael Wenger as described in
/// chapter 3.3.5 of his book "Isosurfaces: Geometry, Topology, and Algorithms".
template<class T> class DualMC {
public:
    // typedefs
    typedef T VolumeDataType;

    /// Extracts the iso surface for a given volume and iso value.
    /// Output is a list of vertices and a list of indices, which connect
    /// vertices to quads.
    /// The quad mesh either uses shared vertex indices or is a quad soup if
    /// desired.
    void build(
        VolumeDataType const * data,
        int32_t const dimX, int32_t const dimY, int32_t const dimZ,
        VolumeDataType const iso,
        bool const generateManifold,
        bool const generateSoup,
        std::vector<Vertex> & vertices,
        std::vector<Quad> & quads
        );

    void buildTris(
        VolumeDataType const * data,
        int32_t const dimX, int32_t const dimY, int32_t const dimZ,
        VolumeDataType const iso,
        bool const generateManifold,
        bool const generateSoup,
        std::vector<Vertex> & vertices,
        std::vector<Tri> & tris
        );

private:

    /// Extract quad mesh with shared vertex indices.
    void buildSharedVerticesQuads(
        VolumeDataType const iso,
        std::vector<Vertex> & vertices,
        std::vector<Quad> & quads
        );
        
    /// Extract quad soup.
    void buildQuadSoup(
        VolumeDataType const iso,
        std::vector<Vertex> & vertices,
        std::vector<Quad> & quads
        );

    void buildSharedVerticesTris(
        VolumeDataType const iso,
        std::vector<Vertex> & vertices,
        std::vector<Tri> & tris
        );

    /*void buildTriSoup(
      VolumeDataType const iso,
      std::vector<Vertex> & vertices,
      std::vector<Tri> & tris
      );*/


protected:

    /// enum with edge codes for a 12-bit voxel edge mask to indicate
    /// grid edges which intersect the ISO surface of classic marching cubes
    enum DMCEdgeCode {
        EDGE0 = 1,
        EDGE1 = 1 << 1,
        EDGE2 = 1 << 2,
        EDGE3 = 1 << 3,
        EDGE4 = 1 << 4,
        EDGE5 = 1 << 5,
        EDGE6 = 1 << 6,
        EDGE7 = 1 << 7,
        EDGE8 = 1 << 8,
        EDGE9 = 1 << 9,
        EDGE10 = 1 << 10,
        EDGE11 = 1 << 11,
        FORCE_32BIT = 0xffffffff
    };

    /// get the 8-bit in-out mask for the voxel corners of the cell cube at (cx,cy,cz)
    /// and the given iso value
    int getCellCode(int32_t const cx, int32_t const cy, int32_t const cz, VolumeDataType const iso) const;

    /// Get the 12-bit dual point code mask, which encodes the traditional
    /// marching cube vertices of the traditional marching cubes face which
    /// corresponds to the dual point.
    /// This is also where the manifold dual marching cubes algorithm is
    /// implemented.
    int getDualPointCode(int32_t const cx, int32_t const cy, int32_t const cz,
      VolumeDataType const iso, DMCEdgeCode const edge) const;

    /// Given a dual point code and iso value, compute the dual point.
    void calculateDualPoint(int32_t const cx, int32_t const cy, int32_t const cz,
      VolumeDataType const iso, int const pointCode, Vertex &v) const;

    /// Get the shared index of a dual point which is uniquly identified by its
    /// cell cube index and a cube edge. The dual point is computed,
    /// if it has not been computed before.
    QuadIndexType getSharedDualPointIndex(int32_t const cx, int32_t const cy, int32_t const cz,
      VolumeDataType const iso, DMCEdgeCode const edge,
      std::vector<Vertex> & vertices);
    
    /// Compute a linearized cell cube index.
    int32_t gA(int32_t const x, int32_t const y, int32_t const z) const;

protected:
    // static lookup tables needed for (manifold) dual marching cubes

    /// Dual Marching Cubes table
    /// Encodes the edge vertices for the 256 marching cubes cases.
    /// A marching cube case produces up to four faces and ,thus, up to four
    /// dual points.
    static int32_t const dualPointsList[256][4];
    
    /// Table which encodes the ambiguous face of cube configurations, which
    /// can cause non-manifold meshes.
    /// Needed for manifold dual marching cubes.
    static uint8_t const problematicConfigs[256];
    
protected:

    /// convenience volume extent array for x-,y-, and z-dimension
    int32_t dims[3];

    /// convenience volume data point
    VolumeDataType const * data;
    
    /// store whether the manifold dual marching cubes algorithm should be
    /// applied.
    bool generateManifold;
    
    /// Dual point key structure for hashing of shared vertices
    struct DualPointKey {
        // a dual point can be uniquely identified by ite linearized volume cell
        // id and point code
        int32_t linearizedCellID;
        int pointCode;
        /// Equal operator for unordered map
        bool operator==(DualPointKey const & other) const;
    };
    
    /// Functor for dual point key hash generation
    struct DualPointKeyHash {
        size_t operator()(DualPointKey const & k) const {
            return size_t(k.linearizedCellID) | (size_t(k.pointCode) << 32u);
        }
    };
    
    /// Hash map for shared vertex index computations
    std::unordered_map<DualPointKey,QuadIndexType,DualPointKeyHash> pointToIndex;
};

// inline function definitions

//------------------------------------------------------------------------------

inline
Vertex::Vertex(){}

//------------------------------------------------------------------------------

inline
Vertex::Vertex(
    VertexComponentsType x,
    VertexComponentsType y,
    VertexComponentsType z
    ) : x(x), y(y), z(z) {}

//------------------------------------------------------------------------------

inline
Vertex::Vertex(Vertex const & v) : x(v.x), y(v.y), z(v.z) {}

//------------------------------------------------------------------------------

inline
Quad::Quad(){}

//------------------------------------------------------------------------------

inline
Quad::Quad(
    QuadIndexType i0,
    QuadIndexType i1,
    QuadIndexType i2,
    QuadIndexType i3
    ) : i0(i0),i1(i1),i2(i2),i3(i3) {}

//------------------------------------------------------------------------------

inline
Tri::Tri() {}

//------------------------------------------------------------------------------

inline
Tri::Tri(
  TriIndexType i0,
  TriIndexType i1,
  TriIndexType i2
) : i0(i0), i1(i1), i2(i2) {}

//------------------------------------------------------------------------------

template<class T> inline
int32_t DualMC<T>::gA(int32_t const x, int32_t const y, int32_t const z) const {
    return x + dims[0] * (y + dims[1] * z);
}

//------------------------------------------------------------------------------
template<class T> inline
bool DualMC<T>::DualPointKey::operator==(typename DualMC<T>::DualPointKey const & other) const {
    return linearizedCellID == other.linearizedCellID && pointCode == other.pointCode;
}

// Copyright (C) 2017, Dominik Wodniok
// This software may be modified and distributed under the terms
// of the BSD 3-Clause license. See the LICENSE.txt file for details.

/// \file   dualmc.tpp
/// \author Dominik Wodniok
/// \date   2009

//------------------------------------------------------------------------------

template<class T> inline
int DualMC<T>::getCellCode(int32_t const cx, int32_t const cy, int32_t const cz, VolumeDataType const iso) const {
  // determine for each cube corner if it is outside or inside
  int code = 0;
  if (data[gA(cx, cy, cz)] >= iso)
    code |= 1;
  if (data[gA(cx + 1, cy, cz)] >= iso)
    code |= 2;
  if (data[gA(cx, cy + 1, cz)] >= iso)
    code |= 4;
  if (data[gA(cx + 1, cy + 1, cz)] >= iso)
    code |= 8;
  if (data[gA(cx, cy, cz + 1)] >= iso)
    code |= 16;
  if (data[gA(cx + 1, cy, cz + 1)] >= iso)
    code |= 32;
  if (data[gA(cx, cy + 1, cz + 1)] >= iso)
    code |= 64;
  if (data[gA(cx + 1, cy + 1, cz + 1)] >= iso)
    code |= 128;
  return code;
}

//------------------------------------------------------------------------------

template<class T> inline
int DualMC<T>::getDualPointCode(int32_t const cx, int32_t const cy, int32_t const cz, VolumeDataType const iso, DMCEdgeCode const edge) const {
  int cubeCode = getCellCode(cx, cy, cz, iso);

  // is manifold dual marching cubes desired?
  if (generateManifold) {
    // The Manifold Dual Marching Cubes approach from Rephael Wenger as described in
    // chapter 3.3.5 of his book "Isosurfaces: Geometry, Topology, and Algorithms"
    // is implemente here.
    // If a problematic C16 or C19 configuration shares the ambiguous face 
    // with another C16 or C19 configuration we simply invert the cube code
    // before looking up dual points. Doing this for these pairs ensures
    // manifold meshes.
    // But this removes the dualism to marching cubes.

    // check if we have a potentially problematic configuration
    uint8_t const direction = problematicConfigs[uint8_t(cubeCode)];
    // If the direction code is in {0,...,5} we have a C16 or C19 configuration.
    if (direction != 255) {
      // We have to check the neighboring cube, which shares the ambiguous
      // face. For this we decode the direction. This could also be done
      // with another lookup table.
      // copy current cube coordinates into an array.
      int32_t neighborCoords[] = { cx,cy,cz };
      // get the dimension of the non-zero coordinate axis
      unsigned int const component = direction >> 1;
      // get the sign of the direction
      int32_t delta = (direction & 1) == 1 ? 1 : -1;
      // modify the correspong cube coordinate
      neighborCoords[component] += delta;
      // have we left the volume in this direction?
      if (neighborCoords[component] >= 0 && neighborCoords[component] < (dims[component] - 1)) {
        // get the cube configuration of the relevant neighbor
        int neighborCubeCode = getCellCode(neighborCoords[0], neighborCoords[1], neighborCoords[2], iso);
        // Look up the neighbor configuration ambiguous face direction.
        // If the direction is valid we have a C16 or C19 neighbor.
        // As C16 and C19 have exactly one ambiguous face this face is
        // guaranteed to be shared for the pair.
        if (problematicConfigs[uint8_t(neighborCubeCode)] != 255) {
          // replace the cube configuration with its inverse.
          cubeCode ^= 0xff;
        }
      }
    }
  }
  for (int i = 0; i < 4; ++i)
    if (dualPointsList[cubeCode][i] & edge) {
      return dualPointsList[cubeCode][i];
    }
  return 0;
}

//------------------------------------------------------------------------------

template<class T> inline
void DualMC<T>::calculateDualPoint(int32_t const cx, int32_t const cy, int32_t const cz, VolumeDataType const iso, int const pointCode, Vertex & v) const {
  // initialize the point with lower voxel coordinates
  v.x = static_cast<VertexComponentsType>(cx);
  v.y = static_cast<VertexComponentsType>(cy);
  v.z = static_cast<VertexComponentsType>(cz);

  // compute the dual point as the mean of the face vertices belonging to the
  // original marching cubes face
  Vertex p;
  p.x = 0;
  p.y = 0;
  p.z = 0;
  int points = 0;

  // sum edge intersection vertices using the point code
  if (pointCode & EDGE0) {
    p.x += ((float)iso - (float)data[gA(cx, cy, cz)]) / ((float)data[gA(cx + 1, cy, cz)] - (float)data[gA(cx, cy, cz)]);
    points++;
  }

  if (pointCode & EDGE1) {
    p.x += 1.0f;
    p.z += ((float)iso - (float)data[gA(cx + 1, cy, cz)]) / ((float)data[gA(cx + 1, cy, cz + 1)] - (float)data[gA(cx + 1, cy, cz)]);
    points++;
  }

  if (pointCode & EDGE2) {
    p.x += ((float)iso - (float)data[gA(cx, cy, cz + 1)]) / ((float)data[gA(cx + 1, cy, cz + 1)] - (float)data[gA(cx, cy, cz + 1)]);
    p.z += 1.0f;
    points++;
  }

  if (pointCode & EDGE3) {
    p.z += ((float)iso - (float)data[gA(cx, cy, cz)]) / ((float)data[gA(cx, cy, cz + 1)] - (float)data[gA(cx, cy, cz)]);
    points++;
  }

  if (pointCode & EDGE4) {
    p.x += ((float)iso - (float)data[gA(cx, cy + 1, cz)]) / ((float)data[gA(cx + 1, cy + 1, cz)] - (float)data[gA(cx, cy + 1, cz)]);
    p.y += 1.0f;
    points++;
  }

  if (pointCode & EDGE5) {
    p.x += 1.0f;
    p.z += ((float)iso - (float)data[gA(cx + 1, cy + 1, cz)]) / ((float)data[gA(cx + 1, cy + 1, cz + 1)] - (float)data[gA(cx + 1, cy + 1, cz)]);
    p.y += 1.0f;
    points++;
  }

  if (pointCode & EDGE6) {
    p.x += ((float)iso - (float)data[gA(cx, cy + 1, cz + 1)]) / ((float)data[gA(cx + 1, cy + 1, cz + 1)] - (float)data[gA(cx, cy + 1, cz + 1)]);
    p.z += 1.0f;
    p.y += 1.0f;
    points++;
  }

  if (pointCode & EDGE7) {
    p.z += ((float)iso - (float)data[gA(cx, cy + 1, cz)]) / ((float)data[gA(cx, cy + 1, cz + 1)] - (float)data[gA(cx, cy + 1, cz)]);
    p.y += 1.0f;
    points++;
  }

  if (pointCode & EDGE8) {
    p.y += ((float)iso - (float)data[gA(cx, cy, cz)]) / ((float)data[gA(cx, cy + 1, cz)] - (float)data[gA(cx, cy, cz)]);
    points++;
  }

  if (pointCode & EDGE9) {
    p.x += 1.0f;
    p.y += ((float)iso - (float)data[gA(cx + 1, cy, cz)]) / ((float)data[gA(cx + 1, cy + 1, cz)] - (float)data[gA(cx + 1, cy, cz)]);
    points++;
  }

  if (pointCode & EDGE10) {
    p.x += 1.0f;
    p.y += ((float)iso - (float)data[gA(cx + 1, cy, cz + 1)]) / ((float)data[gA(cx + 1, cy + 1, cz + 1)] - (float)data[gA(cx + 1, cy, cz + 1)]);
    p.z += 1.0f;
    points++;
  }

  if (pointCode & EDGE11) {
    p.z += 1.0f;
    p.y += ((float)iso - (float)data[gA(cx, cy, cz + 1)]) / ((float)data[gA(cx, cy + 1, cz + 1)] - (float)data[gA(cx, cy, cz + 1)]);
    points++;
  }

  // divide by number of accumulated points
  float invPoints = 1.0f / (float)points;
  p.x *= invPoints;
  p.y *= invPoints;
  p.z *= invPoints;

  // offset point by voxel coordinates
  v.x += p.x;
  v.y += p.y;
  v.z += p.z;
}

//------------------------------------------------------------------------------

template<class T> inline
QuadIndexType DualMC<T>::getSharedDualPointIndex(
  int32_t const cx, int32_t const cy, int32_t const cz,
  VolumeDataType const iso, DMCEdgeCode const edge,
  std::vector<Vertex> & vertices
) {
  // create a key for the dual point from its linearized cell ID and point code
  DualPointKey key;
  key.linearizedCellID = gA(cx, cy, cz);
  key.pointCode = getDualPointCode(cx, cy, cz, iso, edge);

  // have we already computed the dual point?
  auto iterator = pointToIndex.find(key);
  if (iterator != pointToIndex.end()) {
    // just return the dual point index
    return iterator->second;
  }
  else {
    // create new vertex and vertex id
    QuadIndexType newVertexId = static_cast<QuadIndexType>(vertices.size());
    vertices.emplace_back();
    calculateDualPoint(cx, cy, cz, iso, key.pointCode, vertices.back());
    // insert vertex ID into map and also return it
    pointToIndex[key] = newVertexId;
    return newVertexId;
  }
}

//------------------------------------------------------------------------------

template<class T> inline
void DualMC<T>::build(
  VolumeDataType const * data,
  int32_t const dimX, int32_t const dimY, int32_t const dimZ,
  VolumeDataType const iso,
  bool const generateManifold,
  bool const generateSoup,
  std::vector<Vertex> & vertices,
  std::vector<Quad> & quads
) {

  // set members
  this->dims[0] = dimX;
  this->dims[1] = dimY;
  this->dims[2] = dimZ;
  this->data = data;
  this->generateManifold = generateManifold;

  // clear vertices and quad indices
  vertices.clear();
  quads.clear();

  // generate quad soup or shared vertices quad list
  if (generateSoup) {
    buildQuadSoup(iso, vertices, quads);
  }
  else {
    buildSharedVerticesQuads(iso, vertices, quads);
  }
}

template<class T> inline
void DualMC<T>::buildTris(
  VolumeDataType const * data,
  int32_t const dimX, int32_t const dimY, int32_t const dimZ,
  VolumeDataType const iso,
  bool const generateManifold,
  bool const generateSoup,
  std::vector<Vertex> & vertices,
  std::vector<Tri> & tris
) {

  // set members
  this->dims[0] = dimX;
  this->dims[1] = dimY;
  this->dims[2] = dimZ;
  this->data = data;
  this->generateManifold = generateManifold;

  // clear vertices and quad indices
  vertices.clear();
  tris.clear();

  // generate quad soup or shared vertices quad list
  if (generateSoup) {
    //buildTriSoup(iso, vertices, tris);
  }
  else {
    buildSharedVerticesTris(iso, vertices, tris);
  }
}

//------------------------------------------------------------------------------

template<class T> inline
void DualMC<T>::buildQuadSoup(
  VolumeDataType const iso,
  std::vector<Vertex> & vertices,
  std::vector<Quad> & quads
) {

  int32_t const reducedX = dims[0] - 2;
  int32_t const reducedY = dims[1] - 2;
  int32_t const reducedZ = dims[2] - 2;

  Vertex vertex0;
  Vertex vertex1;
  Vertex vertex2;
  Vertex vertex3;
  int pointCode;

  // iterate voxels
  for (int32_t z = 0; z < reducedZ; ++z)
    for (int32_t y = 0; y < reducedY; ++y)
      for (int32_t x = 0; x < reducedX; ++x) {
        // construct quad for x edge
        if (z > 0 && y > 0) {
          // is edge intersected?
          bool const entering = data[gA(x, y, z)] < iso && data[gA(x + 1, y, z)] >= iso;
          bool const exiting = data[gA(x, y, z)] >= iso && data[gA(x + 1, y, z)] < iso;
          if (entering || exiting) {
            // generate quad
            pointCode = getDualPointCode(x, y, z, iso, EDGE0);
            calculateDualPoint(x, y, z, iso, pointCode, vertex0);

            pointCode = getDualPointCode(x, y, z - 1, iso, EDGE2);
            calculateDualPoint(x, y, z - 1, iso, pointCode, vertex1);

            pointCode = getDualPointCode(x, y - 1, z - 1, iso, EDGE6);
            calculateDualPoint(x, y - 1, z - 1, iso, pointCode, vertex2);

            pointCode = getDualPointCode(x, y - 1, z, iso, EDGE4);
            calculateDualPoint(x, y - 1, z, iso, pointCode, vertex3);

            if (entering) {
              vertices.emplace_back(vertex0);
              vertices.emplace_back(vertex1);
              vertices.emplace_back(vertex2);
              vertices.emplace_back(vertex3);
            }
            else {
              vertices.emplace_back(vertex0);
              vertices.emplace_back(vertex3);
              vertices.emplace_back(vertex2);
              vertices.emplace_back(vertex1);
            }
          }
        }

        // construct quad for y edge
        if (z > 0 && x > 0) {
          // is edge intersected?
          bool const entering = data[gA(x, y, z)] < iso && data[gA(x, y + 1, z)] >= iso;
          bool const exiting = data[gA(x, y, z)] >= iso && data[gA(x, y + 1, z)] < iso;
          if (entering || exiting) {
            // generate quad
            pointCode = getDualPointCode(x, y, z, iso, EDGE8);
            calculateDualPoint(x, y, z, iso, pointCode, vertex0);

            pointCode = getDualPointCode(x, y, z - 1, iso, EDGE11);
            calculateDualPoint(x, y, z - 1, iso, pointCode, vertex1);

            pointCode = getDualPointCode(x - 1, y, z - 1, iso, EDGE10);
            calculateDualPoint(x - 1, y, z - 1, iso, pointCode, vertex2);

            pointCode = getDualPointCode(x - 1, y, z, iso, EDGE9);
            calculateDualPoint(x - 1, y, z, iso, pointCode, vertex3);

            if (exiting) {
              vertices.emplace_back(vertex0);
              vertices.emplace_back(vertex1);
              vertices.emplace_back(vertex2);
              vertices.emplace_back(vertex3);
            }
            else {
              vertices.emplace_back(vertex0);
              vertices.emplace_back(vertex3);
              vertices.emplace_back(vertex2);
              vertices.emplace_back(vertex1);
            }
          }
        }

        // construct quad for z edge
        if (x > 0 && y > 0) {
          // is edge intersected?
          bool const entering = data[gA(x, y, z)] < iso && data[gA(x, y, z + 1)] >= iso;
          bool const exiting = data[gA(x, y, z)] >= iso && data[gA(x, y, z + 1)] < iso;
          if (entering || exiting) {
            // generate quad
            pointCode = getDualPointCode(x, y, z, iso, EDGE3);
            calculateDualPoint(x, y, z, iso, pointCode, vertex0);

            pointCode = getDualPointCode(x - 1, y, z, iso, EDGE1);
            calculateDualPoint(x - 1, y, z, iso, pointCode, vertex1);

            pointCode = getDualPointCode(x - 1, y - 1, z, iso, EDGE5);
            calculateDualPoint(x - 1, y - 1, z, iso, pointCode, vertex2);

            pointCode = getDualPointCode(x, y - 1, z, iso, EDGE7);
            calculateDualPoint(x, y - 1, z, iso, pointCode, vertex3);

            if (exiting) {
              vertices.emplace_back(vertex0);
              vertices.emplace_back(vertex1);
              vertices.emplace_back(vertex2);
              vertices.emplace_back(vertex3);
            }
            else {
              vertices.emplace_back(vertex0);
              vertices.emplace_back(vertex3);
              vertices.emplace_back(vertex2);
              vertices.emplace_back(vertex1);
            }
          }
        }
      }

  // generate triangle soup quads
  size_t const numQuads = vertices.size() / 4;
  quads.reserve(numQuads);
  for (size_t i = 0; i < numQuads; ++i) {
    quads.emplace_back( static_cast<QuadIndexType>(i) * 4
                      , static_cast<QuadIndexType>(i) * 4 + 1
                      , static_cast<QuadIndexType>(i) * 4 + 2
                      , static_cast<QuadIndexType>(i) * 4 + 3
                      );
  }
}

//------------------------------------------------------------------------------

template<class T> inline
void DualMC<T>::buildSharedVerticesQuads(
  VolumeDataType const iso,
  std::vector<Vertex> & vertices,
  std::vector<Quad> & quads
) {


  int32_t const reducedX = dims[0] - 2;
  int32_t const reducedY = dims[1] - 2;
  int32_t const reducedZ = dims[2] - 2;

  QuadIndexType i0, i1, i2, i3;

  pointToIndex.clear();

  // iterate voxels
  for (int32_t z = 0; z < reducedZ; ++z)
    for (int32_t y = 0; y < reducedY; ++y)
      for (int32_t x = 0; x < reducedX; ++x) {
        // construct quads for x edge
        if (z > 0 && y > 0) {
          bool const entering = data[gA(x, y, z)] < iso && data[gA(x + 1, y, z)] >= iso;
          bool const exiting = data[gA(x, y, z)] >= iso && data[gA(x + 1, y, z)] < iso;
          if (entering || exiting) {
            // generate quad
            i0 = getSharedDualPointIndex(x, y, z, iso, EDGE0, vertices);
            i1 = getSharedDualPointIndex(x, y, z - 1, iso, EDGE2, vertices);
            i2 = getSharedDualPointIndex(x, y - 1, z - 1, iso, EDGE6, vertices);
            i3 = getSharedDualPointIndex(x, y - 1, z, iso, EDGE4, vertices);

            if (entering) {
              quads.emplace_back(i0, i1, i2, i3);
            }
            else {
              quads.emplace_back(i0, i3, i2, i1);
            }
          }
        }

        // construct quads for y edge
        if (z > 0 && x > 0) {
          bool const entering = data[gA(x, y, z)] < iso && data[gA(x, y + 1, z)] >= iso;
          bool const exiting = data[gA(x, y, z)] >= iso && data[gA(x, y + 1, z)] < iso;
          if (entering || exiting) {
            // generate quad
            i0 = getSharedDualPointIndex(x, y, z, iso, EDGE8, vertices);
            i1 = getSharedDualPointIndex(x, y, z - 1, iso, EDGE11, vertices);
            i2 = getSharedDualPointIndex(x - 1, y, z - 1, iso, EDGE10, vertices);
            i3 = getSharedDualPointIndex(x - 1, y, z, iso, EDGE9, vertices);

            if (exiting) {
              quads.emplace_back(i0, i1, i2, i3);
            }
            else {
              quads.emplace_back(i0, i3, i2, i1);
            }
          }
        }

        // construct quads for z edge
        if (x > 0 && y > 0) {
          bool const entering = data[gA(x, y, z)] < iso && data[gA(x, y, z + 1)] >= iso;
          bool const exiting = data[gA(x, y, z)] >= iso && data[gA(x, y, z + 1)] < iso;
          if (entering || exiting) {
            // generate quad
            i0 = getSharedDualPointIndex(x, y, z, iso, EDGE3, vertices);
            i1 = getSharedDualPointIndex(x - 1, y, z, iso, EDGE1, vertices);
            i2 = getSharedDualPointIndex(x - 1, y - 1, z, iso, EDGE5, vertices);
            i3 = getSharedDualPointIndex(x, y - 1, z, iso, EDGE7, vertices);

            if (exiting) {
              quads.emplace_back(i0, i1, i2, i3);
            }
            else {
              quads.emplace_back(i0, i3, i2, i1);
            }
          }
        }
      }
}

//template<class T> inline
//void DualMC<T>::buildTriSoup(
//  VolumeDataType const iso,
//  std::vector<Vertex> & vertices,
//  std::vector<Tri> & tris
//)
//{
//  int32_t const reducedX = dims[0] - 2;
//  int32_t const reducedY = dims[1] - 2;
//  int32_t const reducedZ = dims[2] - 2;
//
//  Vertex vertex0;
//  Vertex vertex1;
//  Vertex vertex2;
//  Vertex vertex3;
//  int pointCode;
//
//  // iterate voxels
//  for (int32_t z = 0; z < reducedZ; ++z)
//  {
//    for (int32_t y = 0; y < reducedY; ++y)
//    {
//      for (int32_t x = 0; x < reducedX; ++x) {
//        // construct quad for x edge
//        if (z > 0 && y > 0) {
//          // is edge intersected?
//          bool const entering = data[gA(x, y, z)] < iso && data[gA(x + 1, y, z)] >= iso;
//          bool const exiting = data[gA(x, y, z)] >= iso && data[gA(x + 1, y, z)] < iso;
//          if (entering || exiting) {
//            // generate quad
//            pointCode = getDualPointCode(x, y, z, iso, EDGE0);
//            calculateDualPoint(x, y, z, iso, pointCode, vertex0);
//
//            pointCode = getDualPointCode(x, y, z - 1, iso, EDGE2);
//            calculateDualPoint(x, y, z - 1, iso, pointCode, vertex1);
//
//            pointCode = getDualPointCode(x, y - 1, z - 1, iso, EDGE6);
//            calculateDualPoint(x, y - 1, z - 1, iso, pointCode, vertex2);
//
//            pointCode = getDualPointCode(x, y - 1, z, iso, EDGE4);
//            calculateDualPoint(x, y - 1, z, iso, pointCode, vertex3);
//
//            if (entering) {
//              vertices.emplace_back(vertex0);
//              vertices.emplace_back(vertex1);
//              vertices.emplace_back(vertex2);
//              vertices.emplace_back(vertex3);
//            }
//            else {
//              vertices.emplace_back(vertex0);
//              vertices.emplace_back(vertex3);
//              vertices.emplace_back(vertex2);
//              vertices.emplace_back(vertex1);
//            }
//          }
//        }
//
//        // construct quad for y edge
//        if (z > 0 && x > 0) {
//          // is edge intersected?
//          bool const entering = data[gA(x, y, z)] < iso && data[gA(x, y + 1, z)] >= iso;
//          bool const exiting = data[gA(x, y, z)] >= iso && data[gA(x, y + 1, z)] < iso;
//          if (entering || exiting) {
//            // generate quad
//            pointCode = getDualPointCode(x, y, z, iso, EDGE8);
//            calculateDualPoint(x, y, z, iso, pointCode, vertex0);
//
//            pointCode = getDualPointCode(x, y, z - 1, iso, EDGE11);
//            calculateDualPoint(x, y, z - 1, iso, pointCode, vertex1);
//
//            pointCode = getDualPointCode(x - 1, y, z - 1, iso, EDGE10);
//            calculateDualPoint(x - 1, y, z - 1, iso, pointCode, vertex2);
//
//            pointCode = getDualPointCode(x - 1, y, z, iso, EDGE9);
//            calculateDualPoint(x - 1, y, z, iso, pointCode, vertex3);
//
//            if (exiting) {
//              vertices.emplace_back(vertex0);
//              vertices.emplace_back(vertex1);
//              vertices.emplace_back(vertex2);
//              vertices.emplace_back(vertex3);
//            }
//            else {
//              vertices.emplace_back(vertex0);
//              vertices.emplace_back(vertex3);
//              vertices.emplace_back(vertex2);
//              vertices.emplace_back(vertex1);
//            }
//          }
//        }
//
//        // construct quad for z edge
//        if (x > 0 && y > 0) {
//          // is edge intersected?
//          bool const entering = data[gA(x, y, z)] < iso && data[gA(x, y, z + 1)] >= iso;
//          bool const exiting = data[gA(x, y, z)] >= iso && data[gA(x, y, z + 1)] < iso;
//          if (entering || exiting) {
//            // generate quad
//            pointCode = getDualPointCode(x, y, z, iso, EDGE3);
//            calculateDualPoint(x, y, z, iso, pointCode, vertex0);
//
//            pointCode = getDualPointCode(x - 1, y, z, iso, EDGE1);
//            calculateDualPoint(x - 1, y, z, iso, pointCode, vertex1);
//
//            pointCode = getDualPointCode(x - 1, y - 1, z, iso, EDGE5);
//            calculateDualPoint(x - 1, y - 1, z, iso, pointCode, vertex2);
//
//            pointCode = getDualPointCode(x, y - 1, z, iso, EDGE7);
//            calculateDualPoint(x, y - 1, z, iso, pointCode, vertex3);
//
//            if (exiting) {
//              vertices.emplace_back(vertex0);
//              vertices.emplace_back(vertex1);
//              vertices.emplace_back(vertex2);
//              vertices.emplace_back(vertex3);
//            }
//            else {
//              vertices.emplace_back(vertex0);
//              vertices.emplace_back(vertex3);
//              vertices.emplace_back(vertex2);
//              vertices.emplace_back(vertex1);
//            }
//          }
//        }
//      }
//    }
//  }
//  // generate triangle soup quads
//  size_t const numTris = vertices.size() / 3;
//  tris.reserve(numTris);
//  for (size_t i = 0; i < numQuads; ++i) {
//    tris.emplace_back(static_cast<QuadIndexType>(i) * 4
//      , static_cast<QuadIndexType>(i) * 4 + 1
//      , static_cast<QuadIndexType>(i) * 4 + 2
//      , static_cast<QuadIndexType>(i) * 4 + 3
//    );
//  }
//}

template<class T> inline
void DualMC<T>::buildSharedVerticesTris(
  VolumeDataType const iso,
  std::vector<Vertex> & vertices,
  std::vector<Tri> & tris
) {
  int32_t const reducedX = dims[0] - 2;
  int32_t const reducedY = dims[1] - 2;
  int32_t const reducedZ = dims[2] - 2;

  TriIndexType i0, i1, i2, i3;

  pointToIndex.clear();

  // iterate voxels
  for (int32_t z = 0; z < reducedZ; ++z)
    for (int32_t y = 0; y < reducedY; ++y)
      for (int32_t x = 0; x < reducedX; ++x) {
        // construct quads for x edge
        if (z > 0 && y > 0) {
          bool const entering = data[gA(x, y, z)] < iso && data[gA(x + 1, y, z)] >= iso;
          bool const exiting = data[gA(x, y, z)] >= iso && data[gA(x + 1, y, z)] < iso;
          if (entering || exiting) {
            // get quad
            i0 = getSharedDualPointIndex(x, y, z, iso, EDGE0, vertices);
            i1 = getSharedDualPointIndex(x, y, z - 1, iso, EDGE2, vertices);
            i2 = getSharedDualPointIndex(x, y - 1, z - 1, iso, EDGE6, vertices);
            i3 = getSharedDualPointIndex(x, y - 1, z, iso, EDGE4, vertices);

            if (entering) {
              tris.emplace_back(i0, i1, i2);
              tris.emplace_back(i2, i3, i0);
            }
            else {
              tris.emplace_back(i2, i1, i0);
              tris.emplace_back(i0, i3, i2);
            }
          }
        }

        // construct quads for y edge
        if (z > 0 && x > 0) {
          bool const entering = data[gA(x, y, z)] < iso && data[gA(x, y + 1, z)] >= iso;
          bool const exiting = data[gA(x, y, z)] >= iso && data[gA(x, y + 1, z)] < iso;
          if (entering || exiting) {
            // generate quad
            i0 = getSharedDualPointIndex(x, y, z, iso, EDGE8, vertices);
            i1 = getSharedDualPointIndex(x, y, z - 1, iso, EDGE11, vertices);
            i2 = getSharedDualPointIndex(x - 1, y, z - 1, iso, EDGE10, vertices);
            i3 = getSharedDualPointIndex(x - 1, y, z, iso, EDGE9, vertices);

            if (exiting) {
              tris.emplace_back(i0, i1, i2);
              tris.emplace_back(i2, i3, i0);                            
            }
            else {
              tris.emplace_back(i2, i1, i0);
              tris.emplace_back(i0, i3, i2);
            }
          }
        }

        // construct quads for z edge
        if (x > 0 && y > 0) {
          bool const entering = data[gA(x, y, z)] < iso && data[gA(x, y, z + 1)] >= iso;
          bool const exiting = data[gA(x, y, z)] >= iso && data[gA(x, y, z + 1)] < iso;
          if (entering || exiting) {
            // generate quad
            i0 = getSharedDualPointIndex(x, y, z, iso, EDGE3, vertices);
            i1 = getSharedDualPointIndex(x - 1, y, z, iso, EDGE1, vertices);
            i2 = getSharedDualPointIndex(x - 1, y - 1, z, iso, EDGE5, vertices);
            i3 = getSharedDualPointIndex(x, y - 1, z, iso, EDGE7, vertices);

            if (exiting) {
              tris.emplace_back(i0, i1, i2);
              tris.emplace_back(i2, i3, i0);
            }
            else {
              tris.emplace_back(i2, i1, i0);
              tris.emplace_back(i0, i3, i2);              
            }
          }
        }
      }
}

// Copyright (C) 2018, Dominik Wodniok
// This software may be modified and distributed under the terms
// of the BSD 3-Clause license.
// See the LICENSE.txt file for details.

/// \file   dualmc_table.tpp
/// \author Dominik Wodniok
/// \date   2009


//  Coordinate system
//
//       y
//       |
//       |
//       |
//       0-----x
//      /
//     /
//    z
//

// Cell Corners
// (Corners are voxels. Number correspond to Morton codes of corner coordinates)
//
//       2-------------------3
//      /|                  /|
//     / |                 / |
//    /  |                /  |
//   6-------------------7   |
//   |   |               |   |
//   |   |               |   |
//   |   |               |   |
//   |   |               |   |
//   |   0---------------|---1
//   |  /                |  /
//   | /                 | /
//   |/                  |/
//   4-------------------5
//


//         Cell Edges
//  
//       o--------4----------o
//      /|                  /|
//     7 |                 5 |
//    /  |                /  |
//   o--------6----------o   |
//   |   8               |   9
//   |   |               |   |
//   |   |               |   |
//   11  |               10  |
//   |   o--------0------|---o
//   |  /                |  /
//   | 3                 | 1
//   |/                  |/
//   o--------2----------o
//

// Encodes the edge vertices for the 256 marching cubes cases.
// A marching cube case produces up to four faces and ,thus, up to four
// dual points.
template<class T>
int32_t const DualMC<T>::dualPointsList[256][4] = {
{0, 0, 0, 0}, // 0
{EDGE0 | EDGE3 | EDGE8, 0, 0, 0}, // 1
{EDGE0 | EDGE1 | EDGE9, 0, 0, 0}, // 2
{EDGE1 | EDGE3 | EDGE8 | EDGE9, 0, 0, 0}, // 3
{EDGE4 | EDGE7 | EDGE8, 0, 0, 0}, // 4
{EDGE0 | EDGE3 | EDGE4 | EDGE7, 0, 0, 0}, // 5
{EDGE0 | EDGE1 | EDGE9, EDGE4 | EDGE7 | EDGE8, 0, 0}, // 6
{EDGE1 | EDGE3 | EDGE4 | EDGE7 | EDGE9, 0, 0, 0}, // 7
{EDGE4 | EDGE5 | EDGE9, 0, 0, 0}, // 8
{EDGE0 | EDGE3 | EDGE8, EDGE4 | EDGE5 | EDGE9, 0, 0}, // 9
{EDGE0 | EDGE1 | EDGE4 | EDGE5, 0, 0, 0}, // 10
{EDGE1 | EDGE3 | EDGE4 | EDGE5 | EDGE8, 0, 0, 0}, // 11
{EDGE5 | EDGE7 | EDGE8 | EDGE9, 0, 0, 0}, // 12
{EDGE0 | EDGE3 | EDGE5 | EDGE7 | EDGE9, 0, 0, 0}, // 13
{EDGE0 | EDGE1 | EDGE5 | EDGE7 | EDGE8, 0, 0, 0}, // 14
{EDGE1 | EDGE3 | EDGE5 | EDGE7, 0, 0, 0}, // 15
{EDGE2 | EDGE3 | EDGE11, 0, 0, 0}, // 16
{EDGE0 | EDGE2 | EDGE8 | EDGE11, 0, 0, 0}, // 17
{EDGE0 | EDGE1 | EDGE9, EDGE2 | EDGE3 | EDGE11, 0, 0}, // 18
{EDGE1 | EDGE2 | EDGE8 | EDGE9 | EDGE11, 0, 0, 0}, // 19
{EDGE4 | EDGE7 | EDGE8, EDGE2 | EDGE3 | EDGE11, 0, 0}, // 20
{EDGE0 | EDGE2 | EDGE4 | EDGE7 | EDGE11, 0, 0, 0}, // 21
{EDGE0 | EDGE1 | EDGE9, EDGE4 | EDGE7 | EDGE8, EDGE2 | EDGE3 | EDGE11, 0}, // 22
{EDGE1 | EDGE2 | EDGE4 | EDGE7 | EDGE9 | EDGE11, 0, 0, 0}, // 23
{EDGE4 | EDGE5 | EDGE9, EDGE2 | EDGE3 | EDGE11, 0, 0}, // 24
{EDGE0 | EDGE2 | EDGE8 | EDGE11, EDGE4 | EDGE5 | EDGE9, 0, 0}, // 25
{EDGE0 | EDGE1 | EDGE4 | EDGE5, EDGE2 | EDGE3 | EDGE11, 0, 0}, // 26
{EDGE1 | EDGE2 | EDGE4 | EDGE5 | EDGE8 | EDGE11, 0, 0, 0}, // 27
{EDGE5 | EDGE7 | EDGE8 | EDGE9, EDGE2 | EDGE3 | EDGE11, 0, 0}, // 28
{EDGE0 | EDGE2 | EDGE5 | EDGE7 | EDGE9 | EDGE11, 0, 0, 0}, // 29
{EDGE0 | EDGE1 | EDGE5 | EDGE7 | EDGE8, EDGE2 | EDGE3 | EDGE11, 0, 0}, // 30
{EDGE1 | EDGE2 | EDGE5 | EDGE7 | EDGE11, 0, 0, 0}, // 31
{EDGE1 | EDGE2 | EDGE10, 0, 0, 0}, // 32
{EDGE0 | EDGE3 | EDGE8, EDGE1 | EDGE2 | EDGE10, 0, 0}, // 33
{EDGE0 | EDGE2 | EDGE9 | EDGE10, 0, 0, 0}, // 34
{EDGE2 | EDGE3 | EDGE8 | EDGE9 | EDGE10, 0, 0, 0}, // 35
{EDGE4 | EDGE7 | EDGE8, EDGE1 | EDGE2 | EDGE10, 0, 0}, // 36
{EDGE0 | EDGE3 | EDGE4 | EDGE7, EDGE1 | EDGE2 | EDGE10, 0, 0}, // 37
{EDGE0 | EDGE2 | EDGE9 | EDGE10, EDGE4 | EDGE7 | EDGE8, 0, 0}, // 38
{EDGE2 | EDGE3 | EDGE4 | EDGE7 | EDGE9 | EDGE10, 0, 0, 0}, // 39
{EDGE4 | EDGE5 | EDGE9, EDGE1 | EDGE2 | EDGE10, 0, 0}, // 40
{EDGE0 | EDGE3 | EDGE8, EDGE4 | EDGE5 | EDGE9, EDGE1 | EDGE2 | EDGE10, 0}, // 41
{EDGE0 | EDGE2 | EDGE4 | EDGE5 | EDGE10, 0, 0, 0}, // 42
{EDGE2 | EDGE3 | EDGE4 | EDGE5 | EDGE8 | EDGE10, 0, 0, 0}, // 43
{EDGE5 | EDGE7 | EDGE8 | EDGE9, EDGE1 | EDGE2 | EDGE10, 0, 0}, // 44
{EDGE0 | EDGE3 | EDGE5 | EDGE7 | EDGE9, EDGE1 | EDGE2 | EDGE10, 0, 0}, // 45
{EDGE0 | EDGE2 | EDGE5 | EDGE7 | EDGE8 | EDGE10, 0, 0, 0}, // 46
{EDGE2 | EDGE3 | EDGE5 | EDGE7 | EDGE10, 0, 0, 0}, // 47
{EDGE1 | EDGE3 | EDGE10 | EDGE11, 0, 0, 0}, // 48
{EDGE0 | EDGE1 | EDGE8 | EDGE10 | EDGE11, 0, 0, 0}, // 49
{EDGE0 | EDGE3 | EDGE9 | EDGE10 | EDGE11, 0, 0, 0}, // 50
{EDGE8 | EDGE9 | EDGE10 | EDGE11, 0, 0, 0}, // 51
{EDGE4 | EDGE7 | EDGE8, EDGE1 | EDGE3 | EDGE10 | EDGE11, 0, 0}, // 52
{EDGE0 | EDGE1 | EDGE4 | EDGE7 | EDGE10 | EDGE11, 0, 0, 0}, // 53
{EDGE0 | EDGE3 | EDGE9 | EDGE10 | EDGE11, EDGE4 | EDGE7 | EDGE8, 0, 0}, // 54
{EDGE4 | EDGE7 | EDGE9 | EDGE10 | EDGE11, 0, 0, 0}, // 55
{EDGE4 | EDGE5 | EDGE9, EDGE1 | EDGE3 | EDGE10 | EDGE11, 0, 0}, // 56
{EDGE0 | EDGE1 | EDGE8 | EDGE10 | EDGE11, EDGE4 | EDGE5 | EDGE9, 0, 0}, // 57
{EDGE0 | EDGE3 | EDGE4 | EDGE5 | EDGE10 | EDGE11, 0, 0, 0}, // 58
{EDGE4 | EDGE5 | EDGE8 | EDGE10 | EDGE11, 0, 0, 0}, // 59
{EDGE5 | EDGE7 | EDGE8 | EDGE9, EDGE1 | EDGE3 | EDGE10 | EDGE11, 0, 0}, // 60
{EDGE0 | EDGE1 | EDGE5 | EDGE7 | EDGE9 | EDGE10 | EDGE11, 0, 0, 0}, // 61
{EDGE0 | EDGE3 | EDGE5 | EDGE7 | EDGE8 | EDGE10 | EDGE11, 0, 0, 0}, // 62
{EDGE5 | EDGE7 | EDGE10 | EDGE11, 0, 0, 0}, // 63
{EDGE6 | EDGE7 | EDGE11, 0, 0, 0}, // 64
{EDGE0 | EDGE3 | EDGE8, EDGE6 | EDGE7 | EDGE11, 0, 0}, // 65
{EDGE0 | EDGE1 | EDGE9, EDGE6 | EDGE7 | EDGE11, 0, 0}, // 66
{EDGE1 | EDGE3 | EDGE8 | EDGE9, EDGE6 | EDGE7 | EDGE11, 0, 0}, // 67
{EDGE4 | EDGE6 | EDGE8 | EDGE11, 0, 0, 0}, // 68
{EDGE0 | EDGE3 | EDGE4 | EDGE6 | EDGE11, 0, 0, 0}, // 69
{EDGE0 | EDGE1 | EDGE9, EDGE4 | EDGE6 | EDGE8 | EDGE11, 0, 0}, // 70
{EDGE1 | EDGE3 | EDGE4 | EDGE6 | EDGE9 | EDGE11, 0, 0, 0}, // 71
{EDGE4 | EDGE5 | EDGE9, EDGE6 | EDGE7 | EDGE11, 0, 0}, // 72
{EDGE0 | EDGE3 | EDGE8, EDGE4 | EDGE5 | EDGE9, EDGE6 | EDGE7 | EDGE11, 0}, // 73
{EDGE0 | EDGE1 | EDGE4 | EDGE5, EDGE6 | EDGE7 | EDGE11, 0, 0}, // 74
{EDGE1 | EDGE3 | EDGE4 | EDGE5 | EDGE8, EDGE6 | EDGE7 | EDGE11, 0, 0}, // 75
{EDGE5 | EDGE6 | EDGE8 | EDGE9 | EDGE11, 0, 0, 0}, // 76
{EDGE0 | EDGE3 | EDGE5 | EDGE6 | EDGE9 | EDGE11, 0, 0, 0}, // 77
{EDGE0 | EDGE1 | EDGE5 | EDGE6 | EDGE8 | EDGE11, 0, 0, 0}, // 78
{EDGE1 | EDGE3 | EDGE5 | EDGE6 | EDGE11, 0, 0, 0}, // 79
{EDGE2 | EDGE3 | EDGE6 | EDGE7, 0, 0, 0}, // 80
{EDGE0 | EDGE2 | EDGE6 | EDGE7 | EDGE8, 0, 0, 0}, // 81
{EDGE0 | EDGE1 | EDGE9, EDGE2 | EDGE3 | EDGE6 | EDGE7, 0, 0}, // 82
{EDGE1 | EDGE2 | EDGE6 | EDGE7 | EDGE8 | EDGE9, 0, 0, 0}, // 83
{EDGE2 | EDGE3 | EDGE4 | EDGE6 | EDGE8, 0, 0, 0}, // 84
{EDGE0 | EDGE2 | EDGE4 | EDGE6, 0, 0, 0}, // 85
{EDGE0 | EDGE1 | EDGE9, EDGE2 | EDGE3 | EDGE4 | EDGE6 | EDGE8, 0, 0}, // 86
{EDGE1 | EDGE2 | EDGE4 | EDGE6 | EDGE9, 0, 0, 0}, // 87
{EDGE4 | EDGE5 | EDGE9, EDGE2 | EDGE3 | EDGE6 | EDGE7, 0, 0}, // 88
{EDGE0 | EDGE2 | EDGE6 | EDGE7 | EDGE8, EDGE4 | EDGE5 | EDGE9, 0, 0}, // 89
{EDGE0 | EDGE1 | EDGE4 | EDGE5, EDGE2 | EDGE3 | EDGE6 | EDGE7, 0, 0}, // 90
{EDGE1 | EDGE2 | EDGE4 | EDGE5 | EDGE6 | EDGE7 | EDGE8, 0, 0, 0}, // 91
{EDGE2 | EDGE3 | EDGE5 | EDGE6 | EDGE8 | EDGE9, 0, 0, 0}, // 92
{EDGE0 | EDGE2 | EDGE5 | EDGE6 | EDGE9, 0, 0, 0}, // 93
{EDGE0 | EDGE1 | EDGE2 | EDGE3 | EDGE5 | EDGE6 | EDGE8, 0, 0, 0}, // 94
{EDGE1 | EDGE2 | EDGE5 | EDGE6, 0, 0, 0}, // 95
{EDGE1 | EDGE2 | EDGE10, EDGE6 | EDGE7 | EDGE11, 0, 0}, // 96
{EDGE0 | EDGE3 | EDGE8, EDGE1 | EDGE2 | EDGE10, EDGE6 | EDGE7 | EDGE11, 0}, // 97
{EDGE0 | EDGE2 | EDGE9 | EDGE10, EDGE6 | EDGE7 | EDGE11, 0, 0}, // 98
{EDGE2 | EDGE3 | EDGE8 | EDGE9 | EDGE10, EDGE6 | EDGE7 | EDGE11, 0, 0}, // 99
{EDGE4 | EDGE6 | EDGE8 | EDGE11, EDGE1 | EDGE2 | EDGE10, 0, 0}, // 100
{EDGE0 | EDGE3 | EDGE4 | EDGE6 | EDGE11, EDGE1 | EDGE2 | EDGE10, 0, 0}, // 101
{EDGE0 | EDGE2 | EDGE9 | EDGE10, EDGE4 | EDGE6 | EDGE8 | EDGE11, 0, 0}, // 102
{EDGE2 | EDGE3 | EDGE4 | EDGE6 | EDGE9 | EDGE10 | EDGE11, 0, 0, 0}, // 103
{EDGE4 | EDGE5 | EDGE9, EDGE1 | EDGE2 | EDGE10, EDGE6 | EDGE7 | EDGE11, 0}, // 104
{EDGE0 | EDGE3 | EDGE8, EDGE4 | EDGE5 | EDGE9, EDGE1 | EDGE2 | EDGE10, EDGE6 | EDGE7 | EDGE11}, // 105
{EDGE0 | EDGE2 | EDGE4 | EDGE5 | EDGE10, EDGE6 | EDGE7 | EDGE11, 0, 0}, // 106
{EDGE2 | EDGE3 | EDGE4 | EDGE5 | EDGE8 | EDGE10, EDGE6 | EDGE7 | EDGE11, 0, 0}, // 107
{EDGE5 | EDGE6 | EDGE8 | EDGE9 | EDGE11, EDGE1 | EDGE2 | EDGE10, 0, 0}, // 108
{EDGE0 | EDGE3 | EDGE5 | EDGE6 | EDGE9 | EDGE11, EDGE1 | EDGE2 | EDGE10, 0, 0}, // 109
{EDGE0 | EDGE2 | EDGE5 | EDGE6 | EDGE8 | EDGE10 | EDGE11, 0, 0, 0}, // 110
{EDGE2 | EDGE3 | EDGE5 | EDGE6 | EDGE10 | EDGE11, 0, 0, 0}, // 111
{EDGE1 | EDGE3 | EDGE6 | EDGE7 | EDGE10, 0, 0, 0}, // 112
{EDGE0 | EDGE1 | EDGE6 | EDGE7 | EDGE8 | EDGE10, 0, 0, 0}, // 113
{EDGE0 | EDGE3 | EDGE6 | EDGE7 | EDGE9 | EDGE10, 0, 0, 0}, // 114
{EDGE6 | EDGE7 | EDGE8 | EDGE9 | EDGE10, 0, 0, 0}, // 115
{EDGE1 | EDGE3 | EDGE4 | EDGE6 | EDGE8 | EDGE10, 0, 0, 0}, // 116
{EDGE0 | EDGE1 | EDGE4 | EDGE6 | EDGE10, 0, 0, 0}, // 117
{EDGE0 | EDGE3 | EDGE4 | EDGE6 | EDGE8 | EDGE9 | EDGE10, 0, 0, 0}, // 118
{EDGE4 | EDGE6 | EDGE9 | EDGE10, 0, 0, 0}, // 119
{EDGE4 | EDGE5 | EDGE9, EDGE1 | EDGE3 | EDGE6 | EDGE7 | EDGE10, 0, 0}, // 120
{EDGE0 | EDGE1 | EDGE6 | EDGE7 | EDGE8 | EDGE10, EDGE4 | EDGE5 | EDGE9, 0, 0}, // 121
{EDGE0 | EDGE3 | EDGE4 | EDGE5 | EDGE6 | EDGE7 | EDGE10, 0, 0, 0}, // 122
{EDGE4 | EDGE5 | EDGE6 | EDGE7 | EDGE8 | EDGE10, 0, 0, 0}, // 123
{EDGE1 | EDGE3 | EDGE5 | EDGE6 | EDGE8 | EDGE9 | EDGE10, 0, 0, 0}, // 124
{EDGE0 | EDGE1 | EDGE5 | EDGE6 | EDGE9 | EDGE10, 0, 0, 0}, // 125
{EDGE0 | EDGE3 | EDGE8, EDGE5 | EDGE6 | EDGE10, 0, 0}, // 126
{EDGE5 | EDGE6 | EDGE10, 0, 0, 0}, // 127
{EDGE5 | EDGE6 | EDGE10, 0, 0, 0}, // 128
{EDGE0 | EDGE3 | EDGE8, EDGE5 | EDGE6 | EDGE10, 0, 0}, // 129
{EDGE0 | EDGE1 | EDGE9, EDGE5 | EDGE6 | EDGE10, 0, 0}, // 130
{EDGE1 | EDGE3 | EDGE8 | EDGE9, EDGE5 | EDGE6 | EDGE10, 0, 0}, // 131
{EDGE4 | EDGE7 | EDGE8, EDGE5 | EDGE6 | EDGE10, 0, 0}, // 132
{EDGE0 | EDGE3 | EDGE4 | EDGE7, EDGE5 | EDGE6 | EDGE10, 0, 0}, // 133
{EDGE0 | EDGE1 | EDGE9, EDGE4 | EDGE7 | EDGE8, EDGE5 | EDGE6 | EDGE10, 0}, // 134
{EDGE1 | EDGE3 | EDGE4 | EDGE7 | EDGE9, EDGE5 | EDGE6 | EDGE10, 0, 0}, // 135
{EDGE4 | EDGE6 | EDGE9 | EDGE10, 0, 0, 0}, // 136
{EDGE0 | EDGE3 | EDGE8, EDGE4 | EDGE6 | EDGE9 | EDGE10, 0, 0}, // 137
{EDGE0 | EDGE1 | EDGE4 | EDGE6 | EDGE10, 0, 0, 0}, // 138
{EDGE1 | EDGE3 | EDGE4 | EDGE6 | EDGE8 | EDGE10, 0, 0, 0}, // 139
{EDGE6 | EDGE7 | EDGE8 | EDGE9 | EDGE10, 0, 0, 0}, // 140
{EDGE0 | EDGE3 | EDGE6 | EDGE7 | EDGE9 | EDGE10, 0, 0, 0}, // 141
{EDGE0 | EDGE1 | EDGE6 | EDGE7 | EDGE8 | EDGE10, 0, 0, 0}, // 142
{EDGE1 | EDGE3 | EDGE6 | EDGE7 | EDGE10, 0, 0, 0}, // 143
{EDGE2 | EDGE3 | EDGE11, EDGE5 | EDGE6 | EDGE10, 0, 0}, // 144
{EDGE0 | EDGE2 | EDGE8 | EDGE11, EDGE5 | EDGE6 | EDGE10, 0, 0}, // 145
{EDGE0 | EDGE1 | EDGE9, EDGE2 | EDGE3 | EDGE11, EDGE5 | EDGE6 | EDGE10, 0}, // 146
{EDGE1 | EDGE2 | EDGE8 | EDGE9 | EDGE11, EDGE5 | EDGE6 | EDGE10, 0, 0}, // 147
{EDGE4 | EDGE7 | EDGE8, EDGE2 | EDGE3 | EDGE11, EDGE5 | EDGE6 | EDGE10, 0}, // 148
{EDGE0 | EDGE2 | EDGE4 | EDGE7 | EDGE11, EDGE5 | EDGE6 | EDGE10, 0, 0}, // 149
{EDGE0 | EDGE1 | EDGE9, EDGE4 | EDGE7 | EDGE8, EDGE2 | EDGE3 | EDGE11, EDGE5 | EDGE6 | EDGE10}, // 150
{EDGE1 | EDGE2 | EDGE4 | EDGE7 | EDGE9 | EDGE11, EDGE5 | EDGE6 | EDGE10, 0, 0}, // 151
{EDGE4 | EDGE6 | EDGE9 | EDGE10, EDGE2 | EDGE3 | EDGE11, 0, 0}, // 152
{EDGE0 | EDGE2 | EDGE8 | EDGE11, EDGE4 | EDGE6 | EDGE9 | EDGE10, 0, 0}, // 153
{EDGE0 | EDGE1 | EDGE4 | EDGE6 | EDGE10, EDGE2 | EDGE3 | EDGE11, 0, 0}, // 154
{EDGE1 | EDGE2 | EDGE4 | EDGE6 | EDGE8 | EDGE10 | EDGE11, 0, 0, 0}, // 155
{EDGE6 | EDGE7 | EDGE8 | EDGE9 | EDGE10, EDGE2 | EDGE3 | EDGE11, 0, 0}, // 156
{EDGE0 | EDGE2 | EDGE6 | EDGE7 | EDGE9 | EDGE10 | EDGE11, 0, 0, 0}, // 157
{EDGE0 | EDGE1 | EDGE6 | EDGE7 | EDGE8 | EDGE10, EDGE2 | EDGE3 | EDGE11, 0, 0}, // 158
{EDGE1 | EDGE2 | EDGE6 | EDGE7 | EDGE10 | EDGE11, 0, 0, 0}, // 159
{EDGE1 | EDGE2 | EDGE5 | EDGE6, 0, 0, 0}, // 160
{EDGE0 | EDGE3 | EDGE8, EDGE1 | EDGE2 | EDGE5 | EDGE6, 0, 0}, // 161
{EDGE0 | EDGE2 | EDGE5 | EDGE6 | EDGE9, 0, 0, 0}, // 162
{EDGE2 | EDGE3 | EDGE5 | EDGE6 | EDGE8 | EDGE9, 0, 0, 0}, // 163
{EDGE4 | EDGE7 | EDGE8, EDGE1 | EDGE2 | EDGE5 | EDGE6, 0, 0}, // 164
{EDGE0 | EDGE3 | EDGE4 | EDGE7, EDGE1 | EDGE2 | EDGE5 | EDGE6, 0, 0}, // 165
{EDGE0 | EDGE2 | EDGE5 | EDGE6 | EDGE9, EDGE4 | EDGE7 | EDGE8, 0, 0}, // 166
{EDGE2 | EDGE3 | EDGE4 | EDGE5 | EDGE6 | EDGE7 | EDGE9, 0, 0, 0}, // 167
{EDGE1 | EDGE2 | EDGE4 | EDGE6 | EDGE9, 0, 0, 0}, // 168
{EDGE0 | EDGE3 | EDGE8, EDGE1 | EDGE2 | EDGE4 | EDGE6 | EDGE9, 0, 0}, // 169
{EDGE0 | EDGE2 | EDGE4 | EDGE6, 0, 0, 0}, // 170
{EDGE2 | EDGE3 | EDGE4 | EDGE6 | EDGE8, 0, 0, 0}, // 171
{EDGE1 | EDGE2 | EDGE6 | EDGE7 | EDGE8 | EDGE9, 0, 0, 0}, // 172
{EDGE0 | EDGE1 | EDGE2 | EDGE3 | EDGE6 | EDGE7 | EDGE9, 0, 0, 0}, // 173
{EDGE0 | EDGE2 | EDGE6 | EDGE7 | EDGE8, 0, 0, 0}, // 174
{EDGE2 | EDGE3 | EDGE6 | EDGE7, 0, 0, 0}, // 175
{EDGE1 | EDGE3 | EDGE5 | EDGE6 | EDGE11, 0, 0, 0}, // 176
{EDGE0 | EDGE1 | EDGE5 | EDGE6 | EDGE8 | EDGE11, 0, 0, 0}, // 177
{EDGE0 | EDGE3 | EDGE5 | EDGE6 | EDGE9 | EDGE11, 0, 0, 0}, // 178
{EDGE5 | EDGE6 | EDGE8 | EDGE9 | EDGE11, 0, 0, 0}, // 179
{EDGE4 | EDGE7 | EDGE8, EDGE1 | EDGE3 | EDGE5 | EDGE6 | EDGE11, 0, 0}, // 180
{EDGE0 | EDGE1 | EDGE4 | EDGE5 | EDGE6 | EDGE7 | EDGE11, 0, 0, 0}, // 181
{EDGE0 | EDGE3 | EDGE5 | EDGE6 | EDGE9 | EDGE11, EDGE4 | EDGE7 | EDGE8, 0, 0}, // 182
{EDGE4 | EDGE5 | EDGE6 | EDGE7 | EDGE9 | EDGE11, 0, 0, 0}, // 183
{EDGE1 | EDGE3 | EDGE4 | EDGE6 | EDGE9 | EDGE11, 0, 0, 0}, // 184
{EDGE0 | EDGE1 | EDGE4 | EDGE6 | EDGE8 | EDGE9 | EDGE11, 0, 0, 0}, // 185
{EDGE0 | EDGE3 | EDGE4 | EDGE6 | EDGE11, 0, 0, 0}, // 186
{EDGE4 | EDGE6 | EDGE8 | EDGE11, 0, 0, 0}, // 187
{EDGE1 | EDGE3 | EDGE6 | EDGE7 | EDGE8 | EDGE9 | EDGE11, 0, 0, 0}, // 188
{EDGE0 | EDGE1 | EDGE9, EDGE6 | EDGE7 | EDGE11, 0, 0}, // 189
{EDGE0 | EDGE3 | EDGE6 | EDGE7 | EDGE8 | EDGE11, 0, 0, 0}, // 190
{EDGE6 | EDGE7 | EDGE11, 0, 0, 0}, // 191
{EDGE5 | EDGE7 | EDGE10 | EDGE11, 0, 0, 0}, // 192
{EDGE0 | EDGE3 | EDGE8, EDGE5 | EDGE7 | EDGE10 | EDGE11, 0, 0}, // 193
{EDGE0 | EDGE1 | EDGE9, EDGE5 | EDGE7 | EDGE10 | EDGE11, 0, 0}, // 194
{EDGE1 | EDGE3 | EDGE8 | EDGE9, EDGE5 | EDGE7 | EDGE10 | EDGE11, 0, 0}, // 195
{EDGE4 | EDGE5 | EDGE8 | EDGE10 | EDGE11, 0, 0, 0}, // 196
{EDGE0 | EDGE3 | EDGE4 | EDGE5 | EDGE10 | EDGE11, 0, 0, 0}, // 197
{EDGE0 | EDGE1 | EDGE9, EDGE4 | EDGE5 | EDGE8 | EDGE10 | EDGE11, 0, 0}, // 198
{EDGE1 | EDGE3 | EDGE4 | EDGE5 | EDGE9 | EDGE10 | EDGE11, 0, 0, 0}, // 199
{EDGE4 | EDGE7 | EDGE9 | EDGE10 | EDGE11, 0, 0, 0}, // 200
{EDGE0 | EDGE3 | EDGE8, EDGE4 | EDGE7 | EDGE9 | EDGE10 | EDGE11, 0, 0}, // 201
{EDGE0 | EDGE1 | EDGE4 | EDGE7 | EDGE10 | EDGE11, 0, 0, 0}, // 202
{EDGE1 | EDGE3 | EDGE4 | EDGE7 | EDGE8 | EDGE10 | EDGE11, 0, 0, 0}, // 203
{EDGE8 | EDGE9 | EDGE10 | EDGE11, 0, 0, 0}, // 204
{EDGE0 | EDGE3 | EDGE9 | EDGE10 | EDGE11, 0, 0, 0}, // 205
{EDGE0 | EDGE1 | EDGE8 | EDGE10 | EDGE11, 0, 0, 0}, // 206
{EDGE1 | EDGE3 | EDGE10 | EDGE11, 0, 0, 0}, // 207
{EDGE2 | EDGE3 | EDGE5 | EDGE7 | EDGE10, 0, 0, 0}, // 208
{EDGE0 | EDGE2 | EDGE5 | EDGE7 | EDGE8 | EDGE10, 0, 0, 0}, // 209
{EDGE0 | EDGE1 | EDGE9, EDGE2 | EDGE3 | EDGE5 | EDGE7 | EDGE10, 0, 0}, // 210
{EDGE1 | EDGE2 | EDGE5 | EDGE7 | EDGE8 | EDGE9 | EDGE10, 0, 0, 0}, // 211
{EDGE2 | EDGE3 | EDGE4 | EDGE5 | EDGE8 | EDGE10, 0, 0, 0}, // 212
{EDGE0 | EDGE2 | EDGE4 | EDGE5 | EDGE10, 0, 0, 0}, // 213
{EDGE0 | EDGE1 | EDGE9, EDGE2 | EDGE3 | EDGE4 | EDGE5 | EDGE8 | EDGE10, 0, 0}, // 214
{EDGE1 | EDGE2 | EDGE4 | EDGE5 | EDGE9 | EDGE10, 0, 0, 0}, // 215
{EDGE2 | EDGE3 | EDGE4 | EDGE7 | EDGE9 | EDGE10, 0, 0, 0}, // 216
{EDGE0 | EDGE2 | EDGE4 | EDGE7 | EDGE8 | EDGE9 | EDGE10, 0, 0, 0}, // 217
{EDGE0 | EDGE1 | EDGE2 | EDGE3 | EDGE4 | EDGE7 | EDGE10, 0, 0, 0}, // 218
{EDGE4 | EDGE7 | EDGE8, EDGE1 | EDGE2 | EDGE10, 0, 0}, // 219
{EDGE2 | EDGE3 | EDGE8 | EDGE9 | EDGE10, 0, 0, 0}, // 220
{EDGE0 | EDGE2 | EDGE9 | EDGE10, 0, 0, 0}, // 221
{EDGE0 | EDGE1 | EDGE2 | EDGE3 | EDGE8 | EDGE10, 0, 0, 0}, // 222
{EDGE1 | EDGE2 | EDGE10, 0, 0, 0}, // 223
{EDGE1 | EDGE2 | EDGE5 | EDGE7 | EDGE11, 0, 0, 0}, // 224
{EDGE0 | EDGE3 | EDGE8, EDGE1 | EDGE2 | EDGE5 | EDGE7 | EDGE11, 0, 0}, // 225
{EDGE0 | EDGE2 | EDGE5 | EDGE7 | EDGE9 | EDGE11, 0, 0, 0}, // 226
{EDGE2 | EDGE3 | EDGE5 | EDGE7 | EDGE8 | EDGE9 | EDGE11, 0, 0, 0}, // 227
{EDGE1 | EDGE2 | EDGE4 | EDGE5 | EDGE8 | EDGE11, 0, 0, 0}, // 228
{EDGE0 | EDGE1 | EDGE2 | EDGE3 | EDGE4 | EDGE5 | EDGE11, 0, 0, 0}, // 229
{EDGE0 | EDGE2 | EDGE4 | EDGE5 | EDGE8 | EDGE9 | EDGE11, 0, 0, 0}, // 230
{EDGE4 | EDGE5 | EDGE9, EDGE2 | EDGE3 | EDGE11, 0, 0}, // 231
{EDGE1 | EDGE2 | EDGE4 | EDGE7 | EDGE9 | EDGE11, 0, 0, 0}, // 232
{EDGE0 | EDGE3 | EDGE8, EDGE1 | EDGE2 | EDGE4 | EDGE7 | EDGE9 | EDGE11, 0, 0}, // 233
{EDGE0 | EDGE2 | EDGE4 | EDGE7 | EDGE11, 0, 0, 0}, // 234
{EDGE2 | EDGE3 | EDGE4 | EDGE7 | EDGE8 | EDGE11, 0, 0, 0}, // 235
{EDGE1 | EDGE2 | EDGE8 | EDGE9 | EDGE11, 0, 0, 0}, // 236
{EDGE0 | EDGE1 | EDGE2 | EDGE3 | EDGE9 | EDGE11, 0, 0, 0}, // 237
{EDGE0 | EDGE2 | EDGE8 | EDGE11, 0, 0, 0}, // 238
{EDGE2 | EDGE3 | EDGE11, 0, 0, 0}, // 239
{EDGE1 | EDGE3 | EDGE5 | EDGE7, 0, 0, 0}, // 240
{EDGE0 | EDGE1 | EDGE5 | EDGE7 | EDGE8, 0, 0, 0}, // 241
{EDGE0 | EDGE3 | EDGE5 | EDGE7 | EDGE9, 0, 0, 0}, // 242
{EDGE5 | EDGE7 | EDGE8 | EDGE9, 0, 0, 0}, // 243
{EDGE1 | EDGE3 | EDGE4 | EDGE5 | EDGE8, 0, 0, 0}, // 244
{EDGE0 | EDGE1 | EDGE4 | EDGE5, 0, 0, 0}, // 245
{EDGE0 | EDGE3 | EDGE4 | EDGE5 | EDGE8 | EDGE9, 0, 0, 0}, // 246
{EDGE4 | EDGE5 | EDGE9, 0, 0, 0}, // 247
{EDGE1 | EDGE3 | EDGE4 | EDGE7 | EDGE9, 0, 0, 0}, // 248
{EDGE0 | EDGE1 | EDGE4 | EDGE7 | EDGE8 | EDGE9, 0, 0, 0}, // 249
{EDGE0 | EDGE3 | EDGE4 | EDGE7, 0, 0, 0}, // 250
{EDGE4 | EDGE7 | EDGE8, 0, 0, 0}, // 251
{EDGE1 | EDGE3 | EDGE8 | EDGE9, 0, 0, 0}, // 252
{EDGE0 | EDGE1 | EDGE9, 0, 0, 0}, // 253
{EDGE0 | EDGE3 | EDGE8, 0, 0, 0}, // 254
{0, 0, 0, 0} // 255
};

/// Encodes the ambiguous face of cube configurations, which
/// can cause non-manifold meshes.
/// Non-problematic configurations have a value of 255.
/// The first bit of each value actually encodes a positive or negative
/// direction while the second and third bit enumerate the axis.
template<class T>
uint8_t const DualMC<T>::problematicConfigs[256] = {
255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
255,255,255,255,255,255,255,255,255,255,255,255,255,1,0,255,
255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
255,255,255,255,255,255,255,255,255,255,255,3,255,255,2,255,
255,255,255,255,255,255,255,5,255,255,255,255,255,255,5,5,
255,255,255,255,255,255,4,255,255,255,3,3,1,1,255,255,
255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
255,255,255,255,255,255,255,255,255,255,255,5,255,5,255,5,
255,255,255,255,255,255,255,3,255,255,255,255,255,2,255,255,
255,255,255,255,255,3,255,3,255,4,255,255,0,255,0,255,
255,255,255,255,255,255,255,1,255,255,255,0,255,255,255,255,
255,255,255,1,255,255,255,1,255,4,2,255,255,255,2,255,
255,255,255,0,255,2,4,255,255,255,255,0,255,2,255,255,
255,255,255,255,255,255,4,255,255,4,255,255,255,255,255,255
};


} // END: namespace dualmc
#endif // DUALMC_H_INCLUDED
