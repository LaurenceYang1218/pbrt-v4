// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0

#ifndef PBRT_CPU_AGGREGATES_H
#define PBRT_CPU_AGGREGATES_H

#include <pbrt/pbrt.h>

#include <pbrt/cpu/primitive.h>
#include <pbrt/util/parallel.h>

#include <atomic>
#include <memory>
#include <vector>

namespace pbrt {

Primitive CreateAccelerator(const std::string &name, std::vector<Primitive> prims,
                            const ParameterDictionary &parameters);

struct BVHBuildNode;
struct BVHPrimitive;
struct LinearBVHNode;
struct MortonPrimitive;

// BVHAggregate Definition
class BVHAggregate {
  public:
    // BVHAggregate Public Types
    enum class SplitMethod { SAH, HLBVH, Middle, EqualCounts };

    // BVHAggregate Public Methods
    BVHAggregate(std::vector<Primitive> p, int maxPrimsInNode = 1,
                 SplitMethod splitMethod = SplitMethod::SAH);

    static BVHAggregate *Create(std::vector<Primitive> prims,
                                const ParameterDictionary &parameters);

    Bounds3f Bounds() const;
    pstd::optional<ShapeIntersection> Intersect(const Ray &ray, Float tMax) const;
    bool IntersectP(const Ray &ray, Float tMax) const;

  private:
    // BVHAggregate Private Methods
    BVHBuildNode *buildRecursive(ThreadLocal<Allocator> &threadAllocators,
                                 pstd::span<BVHPrimitive> bvhPrimitives,
                                 std::atomic<int> *totalNodes,
                                 std::atomic<int> *orderedPrimsOffset,
                                 std::vector<Primitive> &orderedPrims);
    BVHBuildNode *buildHLBVH(Allocator alloc,
                             const std::vector<BVHPrimitive> &primitiveInfo,
                             std::atomic<int> *totalNodes,
                             std::vector<Primitive> &orderedPrims);
    BVHBuildNode *emitLBVH(BVHBuildNode *&buildNodes,
                           const std::vector<BVHPrimitive> &primitiveInfo,
                           MortonPrimitive *mortonPrims, int nPrimitives, int *totalNodes,
                           std::vector<Primitive> &orderedPrims,
                           std::atomic<int> *orderedPrimsOffset, int bitIndex);
    BVHBuildNode *buildUpperSAH(Allocator alloc,
                                std::vector<BVHBuildNode *> &treeletRoots, int start,
                                int end, std::atomic<int> *totalNodes) const;
    int flattenBVH(BVHBuildNode *node, int *offset);

    // BVHAggregate Private Members
    int maxPrimsInNode;
    std::vector<Primitive> primitives;
    SplitMethod splitMethod;
    LinearBVHNode *nodes = nullptr;
};

struct KdTreeNode;
struct BoundEdge;

// KdTreeAggregate Definition
class KdTreeAggregate {
  public:
    // KdTreeAggregate Public Methods
    KdTreeAggregate(std::vector<Primitive> p, int isectCost = 5, int traversalCost = 1,
                    Float emptyBonus = 0.5, int maxPrims = 1, int maxDepth = -1);
    static KdTreeAggregate *Create(std::vector<Primitive> prims,
                                   const ParameterDictionary &parameters);
    pstd::optional<ShapeIntersection> Intersect(const Ray &ray, Float tMax) const;

    Bounds3f Bounds() const { return bounds; }

    bool IntersectP(const Ray &ray, Float tMax) const;

  private:
    // KdTreeAggregate Private Methods
    void buildTree(int nodeNum, const Bounds3f &bounds,
                   const std::vector<Bounds3f> &primBounds,
                   pstd::span<const int> primNums, int depth,
                   std::vector<BoundEdge> edges[3], pstd::span<int> prims0,
                   pstd::span<int> prims1, int badRefines);

    // KdTreeAggregate Private Members
    int isectCost, traversalCost, maxPrims;
    Float emptyBonus;
    std::vector<Primitive> primitives;
    std::vector<int> primitiveIndices;
    KdTreeNode *nodes;
    int nAllocedNodes, nextFreeNode;
    Bounds3f bounds;
};

struct UniformGridVoxel;

class UniformGridAggregate {
  public:
    UniformGridAggregate(std::vector<Primitive> p);
    static UniformGridAggregate *Create(std::vector<Primitive> prims, 
                                        const ParameterDictionary &parameters);
    ~UniformGridAggregate();

    pstd::optional<ShapeIntersection> Intersect(const Ray &ray, Float tMax) const;
    
    Bounds3f Bounds() const { return bounds; }

    bool IntersectP(const Ray &ray, Float tMax) const;
  private:
    int posToVoxel(const Point3f &p, int axis) const {
      int v = (int)((p[axis] - bounds.pMin[axis]) * invWidth[axis]);
      return Clamp(v, 0, nVoxels[axis]-1);
    }

    float voxelToPos(int p, int axis) const {
      return bounds.pMin[axis] + p * width[axis];
    }
    
    inline int offset(int x, int y, int z) const {
      return z * nVoxels[0] * nVoxels[1] + y * nVoxels[0] + x;
    }
    
    int nVoxels[3];
    UniformGridVoxel **voxels;
    Vector3f width, invWidth;
    std::vector<Primitive> primitives;
    std::vector<int> primitiveIndexes;
    Bounds3f bounds;
    pthread_rwlock_t rwMutex;
};

}  // namespace pbrt

#endif  // PBRT_CPU_AGGREGATES_H
