/* The copyright in this software is being made available under the BSD
 * Licence, included below.  This software may be subject to other third
 * party and contributor rights, including patent rights, and no such
 * rights are granted under this licence.
 *
 * Copyright (c) 2017-2018, ISO/IEC
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the distribution.
 *
 * * Neither the name of the ISO/IEC nor the names of its contributors
 *   may be used to endorse or promote products derived from this
 *   software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include "geometry.h"

#include "DualLutCoder.h"
#include "OctreeNeighMap.h"
#include "geometry_octree.h"
#include "geometry_intra_pred.h"
#include "io_hls.h"
#include "tables.h"
#include "quantization.h"

#include <set>
#include <random>

namespace pcc {

//============================================================================

enum class DirectMode
{
  kUnavailable,
  kAllPointSame,
  kTwoPoints
};

//============================================================================

class GeometryOctreeEncoder : protected GeometryOctreeContexts {
public:
  GeometryOctreeEncoder(
    const GeometryParameterSet& gps,
    const GeometryBrickHeader& gbh,
    const GeometryOctreeContexts& ctxtMem,
    EntropyEncoder* arithmeticEncoder);

  GeometryOctreeEncoder(const GeometryOctreeEncoder&) = default;
  GeometryOctreeEncoder(GeometryOctreeEncoder&&) = default;
  GeometryOctreeEncoder& operator=(const GeometryOctreeEncoder&) = default;
  GeometryOctreeEncoder& operator=(GeometryOctreeEncoder&&) = default;

  void beginOctreeLevel(const Vec3<int>& planarDepth);

  int encodePositionLeafNumPoints(int count);

  int encodePlanarMode(
    OctreeNodePlanar& planar,
    int plane,
    int dist,
    int neighb,
    int& h,
    int planeId,
    int contextAngle);

  void determinePlanarMode(
    int planeId,
    OctreeNodePlanar& planar,
    OctreePlanarBuffer::Row* planeBuffer,
    int coord1,
    int coord2,
    int coord3,
    uint8_t neighPattern,
    int planarProb[3],
    int planarRate[3],
    int contextAngle);

  void determinePlanarMode(
    int occupancy,
    const bool planarEligible[3],
    PCCOctree3Node& child,
    OctreeNodePlanar& planar,
    uint8_t neighPattern,
    int planarProb[3],
    int contextAngle,
    int contextAnglePhiX,
    int contextAnglePhiY);

  void encodeOccupancyNeighZ(
    int mappedOccupancy,
    int mappedOccIsPredicted,
    int mappedOccPrediction,
    int mappedOccAdjGt0,
    int mappedOccAdjGt1,
    int mappedOccAdjUnocc,
    int mappedPlanarMaskX,
    int mappedFixedMaskX0,
    bool planarPossibleX,
    int mappedPlanarMaskY,
    int mappedFixedMaskY0,
    bool planarPossibleY,
    int mappedPlanarMaskZ,
    int mappedFixedMaskZ0,
    bool planarPossibleZ);

  void encodeOccupancyNeighNZ(
    int neighPattern,
    int mappedOccupancy,
    int mappedOccIsPredicted,
    int mappedOccPrediction,
    int mappedOccAdjGt0,
    int mappedOccAdjGt1,
    int mappedOccAdjUnocc,
    int mappedPlanarMaskX,
    int mappedFixedMaskX0,
    bool planarPossibleX,
    int mappedPlanarMaskY,
    int mappedFixedMaskY0,
    bool planarPossibleY,
    int mappedPlanarMaskZ,
    int mappedFixedMaskZ0,
    bool planarPossibleZ);

  void encodeOccupancyBitwise(
    int neighPattern,
    int mappedOccupancy,
    int mappedOccIsPredicted,
    int mappedOccPrediction,
    int mappedOccAdjGt0,
    int mappedOccAdjGt1,
    int mappedOccAdjUnocc,
    int mappedPlanarMaskX,
    int mappedFixedMaskX0,
    bool planarPossibleX,
    int mappedPlanarMaskY,
    int mappedFixedMaskY0,
    bool planarPossibleY,
    int mappedPlanarMaskZ,
    int mappedFixedMaskZ0,
    bool planarPossibleZ);

  void encodeOccupancyBytewise(int neighPattern, int mappedOccupancy);

  void encodeOccupancy(
    const GeometryNeighPattern& gnp,
    int occupancy,
    int occupancyIsPredicted,
    int occupancyPrediction,
    int planarMaskX,
    int planarMaskY,
    int planarMaskZ,
    bool planarPossibleX,
    bool planarPossibleY,
    bool planarPossibleZ);

  void encodeOrdered2ptPrefix(
    const point_t points[2], Vec3<bool> directIdcm, Vec3<int>& nodeSizeLog2);

  void encodePointPosition(
    const Vec3<int>& nodeSizeLog2AfterPlanar, const Vec3<int32_t>& pos);

  void encodePointPositionAngular(
    const Vec3<int>& nodeSizeLog2,
    const Vec3<int>& nodeSizeLog2AfterPlanar,
    const Vec3<int32_t>& pos,
    const PCCOctree3Node& node,
    const OctreeNodePlanar& planar,
    const Vec3<int>& headPos,
    const int* zLaser,
    const int* thetaLaser,
    int numLasers);

  void encodeQpOffset(int dqp);

  void encodeIsIdcm(DirectMode mode);

  void encodeDirectPosition(
    DirectMode mode,
    bool geom_unique_points_flag,
    bool joint_2pt_idcm_enabled_flag,
    const Vec3<int>& nodeSizeLog2,
    int shiftBits,
    PCCOctree3Node& node,
    OctreeNodePlanar& planar,
    PCCPointSet3& pointCloud,
    bool angularIdcm,
    const Vec3<int>& headPos,
    const int* zLaser,
    const int* thetaLaser,
    int numLasers);

  void encodeThetaRes(int ThetaRes);

  const GeometryOctreeContexts& getCtx() const { return *this; }

public:
  // selects between the bitwise and bytewise occupancy coders
  bool _useBitwiseOccupancyCoder;

  const uint8_t* _neighPattern64toR1;

  EntropyEncoder* _arithmeticEncoder;

  // Planar state
  OctreePlanarState _planar;

  // Azimuthal buffer
  std::vector<int> _phiBuffer;

  // azimuthal elementary shifts
  AzimuthalPhiZi _phiZi;
};

//============================================================================

GeometryOctreeEncoder::GeometryOctreeEncoder(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  const GeometryOctreeContexts& ctxtMem,
  EntropyEncoder* arithmeticEncoder)
  : GeometryOctreeContexts(ctxtMem)
  , _useBitwiseOccupancyCoder(gps.bitwise_occupancy_coding_flag)
  , _neighPattern64toR1(neighPattern64toR1(gps))
  , _arithmeticEncoder(arithmeticEncoder)
  , _planar(gps)
  , _phiBuffer(gps.geom_angular_num_lidar_lasers(), 0x80000000)
  , _phiZi(
      gps.geom_angular_num_lidar_lasers(), gps.geom_angular_num_phi_per_turn)
{
  if (!_useBitwiseOccupancyCoder && !gbh.entropy_continuation_flag) {
    for (int i = 0; i < 10; i++)
      _bytewiseOccupancyCoder[i].init(kDualLutOccupancyCoderInit[i]);
  }
}

//============================================================================

void
GeometryOctreeEncoder::beginOctreeLevel(const Vec3<int>& planarDepth)
{
  for (int i = 0; i < 10; i++) {
    _bytewiseOccupancyCoder[i].resetLut();
  }

  _planar.initPlanes(planarDepth);
}

//============================================================================
// Encode the number of points in a leaf node of the octree.

int
GeometryOctreeEncoder::encodePositionLeafNumPoints(int count)//编码叶子节点的点数
{
  if (count == 1) {
    _arithmeticEncoder->encode(1, _ctxSinglePointPerBlock);//为一则直接二进制编1
  } else {
    _arithmeticEncoder->encode(0, _ctxSinglePointPerBlock);//大于一则编0，然后编点数，解码端解出0即代表非1
    _arithmeticEncoder->encodeExpGolomb(
      uint32_t(count - 2), 0, _ctxPointCountPerBlock);
  }

  return count;
}

//============================================================================

int
GeometryOctreeEncoder::encodePlanarMode(
  OctreeNodePlanar& node,
  int plane,
  int dist,
  int neighb,
  int& h,
  int planeId,
  int contextAngle)
{
  const int mask0 = (1 << planeId);//当前要编码的维度
  const int mask1[3] = {6, 5, 3};//110,101,011,分别指示x,y,z是否为非平面

  bool isPlanar = node.planarMode & mask0;//当前维度是否为平面
  int planeBit = (node.planePosBits & mask0) == 0 ? 0 : 1;//当前维度平面位置

  int discreteDist = (dist <= (2 >> OctreePlanarBuffer::shiftAb) ? 0 : 1);//距离为量化为远、近两个等级
  _arithmeticEncoder->encode(isPlanar, _ctxPlanarMode[planeId]);//编码isPlanar，上下文为维度信息

  if (!isPlanar) {//如果不满足平面，更新平面信息并退出
    node.planarPossible &= mask1[planeId];
    return -1;
  }

  // encode the plane index
  if (contextAngle == -1) {  // angular mode off 如果角度模式没有满足资格，则用原始的平面模式的上下文
    if (plane < 0) {//如果最近已编码节点不是平面，则上下文为维度信息
      _arithmeticEncoder->encode(planeBit, _ctxPlanarPlaneLastIndexZ[planeId]);
      h =
        approxSymbolProbability(planeBit, _ctxPlanarPlaneLastIndexZ[planeId]);
    } else {
      discreteDist += (dist <= (16 >> OctreePlanarBuffer::shiftAb) ? 0 : 1);//将距离变为三档：近、不太远、远
      int lastIndexPlane2d = plane + (discreteDist << 1);//将最近点的距离与平面位置合并为一种上下文，共6个取值
      _arithmeticEncoder->encode(
        planeBit, _ctxPlanarPlaneLastIndex[planeId][neighb][lastIndexPlane2d]);//维度信息、垂直方向邻居占位信息、最近点信息，共3*4*6=72上下文
      h = approxSymbolProbability(
        planeBit, _ctxPlanarPlaneLastIndex[planeId][neighb][lastIndexPlane2d]);
    }
  } else {               // angular mode on
    if (planeId == 2) {  // angular角度模式满足资格，如果是z平面，则采用垂直角
      _arithmeticEncoder->encode(
        planeBit, _ctxPlanarPlaneLastIndexAngular[contextAngle]);
      h = approxSymbolProbability(
        planeBit, _ctxPlanarPlaneLastIndexAngular[contextAngle]);
    } else {  // azimuthal
      _arithmeticEncoder->encode(//如果是x平面或者y平面，则采用水平角
        planeBit, _ctxPlanarPlaneLastIndexAngularPhi[contextAngle]);
      h = approxSymbolProbability(
        planeBit, _ctxPlanarPlaneLastIndexAngularPhi[contextAngle]);
    }
  }
  return planeBit;
}

//============================================================================
// determine Planar mode for one direction

void
GeometryOctreeEncoder::determinePlanarMode(
  int planeId,
  OctreeNodePlanar& planar,
  OctreePlanarBuffer::Row* planeBuffer,
  int coord1,
  int coord2,
  int coord3,
  uint8_t neighPattern,
  int planarProb[3],
  int planarRate[3],
  int contextAngle)
{
  const int kPlanarChildThreshold = 63;
  const int kAdjNeighIdxFromPlanePos[3][2] = {1, 0, 2, 3, 4, 5};
  const int planeSelector = 1 << planeId;

  OctreePlanarBuffer::Elmt* row;//Buffer是用来存储当前层已编码节点的平面信息
  int rowLen = OctreePlanarBuffer::rowSize;//buffer长度
  int closestPlanarFlag;//buffer中存储的最近节点的平面位置，低平面或者高平面
  int closestDist;//buffer中存储的最近节点距离当前节点的距离，被量化为三档：近，不太远，远

  if (!planeBuffer) {//角度模式时不开启buffer
    // angular: buffer disabled
    closestPlanarFlag = 0;
    closestDist = 0;
  } else {
    coord1 =
      (coord1 & OctreePlanarBuffer::maskAb) >> OctreePlanarBuffer::shiftAb;
    coord2 =
      (coord2 & OctreePlanarBuffer::maskAb) >> OctreePlanarBuffer::shiftAb;
    coord3 = coord3 & OctreePlanarBuffer::maskC;//节点应该根据z值存储到相同的buffer

    row = planeBuffer[coord3];//对应维度的缓存区

    int minDist = std::abs(coord1 - int(row[rowLen - 1].a))//最小距离，如果是z平面则根据x、y坐标计算
      + std::abs(coord2 - int(row[rowLen - 1].b));
    int idxMinDist = rowLen - 1;

    for (int idxP = 0; idxP < rowLen - 1; idxP++) {//遍历缓存区中的节点，寻找最小距离
      int dist0 = std::abs(coord1 - int(row[idxP].a))
        + std::abs(coord2 - int(row[idxP].b));
      if (dist0 < minDist) {
        idxMinDist = idxP;
        minDist = dist0;
      }
    }

    // push closest point front
    row[rowLen - 1] = row[idxMinDist];//将寻找到的最小距离的的节点放在最后面

    closestPlanarFlag = row[idxMinDist].planeIdx;//记录获取到的最小距离节点的平面位置，即低平面还是高平面
    closestDist = minDist;//记录最小距离节点的距离

    for (int idxP = 0; idxP < rowLen - 1; idxP++) {
      row[idxP] = row[idxP + 1];//由于最小距离节点在缓存区中的位置发生了变化，将缓存区中的节点重新排序
    }
  }
  const int kAdjNeighIdxFromPlaneMask[3] = {0, 2, 4};//用于获取垂直方向邻居的占位情况
  int adjNeigh = (neighPattern >> kAdjNeighIdxFromPlaneMask[planeId]) & 3;//垂直方向邻居占位情况
  int planeBit = encodePlanarMode(//编码平面位置，低平面还是高平面
    planar, closestPlanarFlag, closestDist, adjNeigh, planarProb[planeId],
    planeId, contextAngle);

  bool isPlanar = (planar.planarMode & planeSelector)//记录当前节点是否为平面，注意，这个只是用来更新概率，并非真正的isPlanar，真正的在encodePlanarMode函数里面
    && planarProb[planeId] > kPlanarChildThreshold;

  planarRate[planeId] =//根据isPlanar的信息来更新概率
    (255 * planarRate[planeId] + (isPlanar ? 256 * 8 : 0) + 128) >> 8;

  if (planeBuffer) {//将当前节点的平面信息、xy坐标存入缓存区中（如果是z平面的话）
    row[rowLen - 1] = {unsigned(coord1), planeBit, unsigned(coord2)};
  }
}

//============================================================================
// determine Planar mode for all directions

void
GeometryOctreeEncoder::determinePlanarMode(
  int occupancy,
  const bool planarEligible[3],
  PCCOctree3Node& node,
  OctreeNodePlanar& planar,
  uint8_t neighPattern,
  int planarProb[3],
  int contextAngle,
  int contextAnglePhiX,
  int contextAnglePhiY)
{
  auto& planeBuffer = _planar._planarBuffer;//buffer

  // determine what planes exist in occupancy
  setPlanesFromOccupancy(occupancy, planar);

  uint8_t planarEligibleMask = 0;
  planarEligibleMask |= planarEligible[2] << 2;
  planarEligibleMask |= planarEligible[1] << 1;
  planarEligibleMask |= planarEligible[0] << 0;
  planar.planarMode &= planarEligibleMask;
  planar.planePosBits &= planarEligibleMask;

  int xx = node.pos[0];
  int yy = node.pos[1];
  int zz = node.pos[2];

  // planar x
  if (planarEligible[0]) {
    determinePlanarMode(
      0, planar, planeBuffer.getBuffer(0), yy, zz, xx, neighPattern,//xyz三个维度上的函数差距为xx yy zz的顺序以及角度
      planarProb, _planar._rate.data(), contextAnglePhiX);//x水平角
  }
  // planar y
  if (planarEligible[1]) {
    determinePlanarMode(
      1, planar, planeBuffer.getBuffer(1), xx, zz, yy, neighPattern,
      planarProb, _planar._rate.data(), contextAnglePhiY);//y水平角
  }
  // planar z
  if (planarEligible[2]) {
    determinePlanarMode(
      2, planar, planeBuffer.getBuffer(2), xx, yy, zz, neighPattern,
      planarProb, _planar._rate.data(), contextAngle);//垂直角
  }
}

//-------------------------------------------------------------------------
// encode occupancy bits (neighPattern10 == 0 case)

void
GeometryOctreeEncoder::encodeOccupancyNeighZ(
  int mappedOccupancy,
  int mappedOccIsPredicted,
  int mappedOccPrediction,
  int mappedOccAdjGt0,
  int mappedOccAdjGt1,
  int mappedOccAdjUnocc,
  int mappedPlanarMaskX,
  int mappedFixedMaskX0,
  bool planarPossibleX,
  int mappedPlanarMaskY,
  int mappedFixedMaskY0,
  bool planarPossibleY,
  int mappedPlanarMaskZ,
  int mappedFixedMaskZ0,
  bool planarPossibleZ)
{
  int numOccupiedAcc = 0;
  //由于NC=0的时候经历了NO=1的判定，因此至少两个子节点被占据
  int maxPerPlaneX = 4 - (mappedPlanarMaskX ? 2 : 1);//如果在x维度上是平面，则该平面上未被占据的最大个数是2，为3的话
  int maxPerPlaneY = 4 - (mappedPlanarMaskY ? 2 : 1);
  int maxPerPlaneZ = 4 - (mappedPlanarMaskZ ? 2 : 1);
  bool sure_planarityX = mappedPlanarMaskX || !planarPossibleX;
  bool sure_planarityY = mappedPlanarMaskY || !planarPossibleY;
  bool sure_planarityZ = mappedPlanarMaskZ || !planarPossibleZ;

  int maskedOccupancy =
    mappedPlanarMaskX | mappedPlanarMaskY | mappedPlanarMaskZ;

  int coded0X[2] = {0, 0};
  int coded0Y[2] = {0, 0};
  int coded0Z[2] = {0, 0};
  if (maskedOccupancy) {
    for (int i = 0; i < 8; i++) {
      if ((maskedOccupancy >> i) & 1) {
        coded0X[(mappedFixedMaskX0 >> i) & 1]++;
        coded0Y[(mappedFixedMaskY0 >> i) & 1]++;
        coded0Z[(mappedFixedMaskZ0 >> i) & 1]++;
      }
    }
  }

  for (int i = 0; i < 8; i++) {
    int bitIdx = kOccBitCodingOrder[i];
    if ((maskedOccupancy >> bitIdx) & 1)
      continue;

    int bitAdjGt0 = (mappedOccAdjGt0 >> bitIdx) & 1;
    int bitAdjGt1 = (mappedOccAdjGt1 >> bitIdx) & 1;
    int bitAdjUnocc = (mappedOccAdjUnocc >> bitIdx) & 1;
    int numAdj = bitAdjGt0 + bitAdjGt1;
    int idxAdj = bitAdjUnocc + 2 * numAdj;
    if (i > 4) {
      static const int8_t kCtxIdxAdjReduc567[6] = {0, 0, 1, 2, 3, 3};
      idxAdj = kCtxIdxAdjReduc567[idxAdj];
    }

    int ctxIdxMapIdx = 3 * idxAdj;
    if (!maskedOccupancy) {
      int bitIsPredicted = (mappedOccIsPredicted >> bitIdx) & 1;
      int bitPrediction = (mappedOccPrediction >> bitIdx) & 1;
      ctxIdxMapIdx = 3 * idxAdj + bitIsPredicted + bitPrediction;
    }

    // NB: There must be at least minOccupied child nodes
    //  -- avoid coding the occupancyBit if it is implied.
    int mask0X = (mappedFixedMaskX0 >> bitIdx) & 1;
    bool bitIsOneX = (sure_planarityX && coded0X[mask0X] >= maxPerPlaneX)
      || (coded0X[0] + coded0X[1] >= 6);

    int mask0Y = (mappedFixedMaskY0 >> bitIdx) & 1;
    bool bitIsOneY = (sure_planarityY && coded0Y[mask0Y] >= maxPerPlaneY)
      || (coded0Y[0] + coded0Y[1] >= 6);

    int mask0Z = (mappedFixedMaskZ0 >> bitIdx) & 1;
    bool bitIsOneZ = (sure_planarityZ && coded0Z[mask0Z] >= maxPerPlaneZ)
      || (coded0Z[0] + coded0Z[1] >= 6);

    // masking for planar is here
    int bit = (mappedOccupancy >> bitIdx) & 1;
    if (!(bitIsOneX || bitIsOneY || bitIsOneZ)) {
      int ctxIdx;
      auto& ctxIdxMap = _ctxIdxMaps[ctxIdxMapIdx];
      ctxIdx = ctxIdxMap.evolve(bit, &ctxIdxMap[i][numOccupiedAcc]);

      _arithmeticEncoder->encode(bit, _ctxOccupancy[ctxIdx]);

      if (!bit) {
        coded0X[mask0X]++;
        coded0Y[mask0Y]++;
        coded0Z[mask0Z]++;
      }
    }

    numOccupiedAcc += bit;
  }
}

//-------------------------------------------------------------------------
// encode occupancy bits (neighPattern10 != 0 case)

void
GeometryOctreeEncoder::encodeOccupancyNeighNZ(
  int neighPattern,
  int mappedOccupancy,
  int mappedOccIsPredicted,
  int mappedOccPrediction,
  int mappedOccAdjGt0,
  int mappedOccAdjGt1,
  int mappedOccAdjUnocc,
  int mappedPlanarMaskX,
  int mappedFixedMaskX0,
  bool planarPossibleX,
  int mappedPlanarMaskY,
  int mappedFixedMaskY0,
  bool planarPossibleY,
  int mappedPlanarMaskZ,
  int mappedFixedMaskZ0,
  bool planarPossibleZ)
{
  // code occupancy using the neighbour configuration context
  // with reduction from 64 states to 9 (or 6).
  int neighPatternR1 = _neighPattern64toR1[neighPattern];//64种状态映射为9种

  //  int neighPattern9 = kNeighPattern64to9[neighPattern];
  int neighPattern5 = kNeighPattern9to5[neighPatternR1];//9种继续映射为5种
  int neighPattern3 = kNeighPattern9to3[neighPatternR1];//5种映射为3种

  //  int neighPattern7 = kNeighPattern10to7[neighPattern10];
  //  int neighPattern5 = kNeighPattern7to5[neighPattern7];

  uint32_t partialOccupancy = 0;
  //planarPossible只能说明当为0时一定不是平面，因此只有0信息有用
  //注意planarMask的取值，在映射之前，在高平面时其取值为0x0f,数值1指示的是对应的子节点是否一定不被占据，而不是是否占据
  //planarPossibleX是根据isplanar判定的，其默认为1，如果为1，说明不了什么，但是如果其为0，说明一定不是平面
  //而planarMask则作用与planarPossible相反，其默认为0，因此如果是0则说明不了什么，但是非0，则说明一定为平面且可以推断出平面的位置
  //统计结果也可以验证这个结论，当mappedPlanarMaskX非0时，planarPossibleX一定为1；当planarPossibleX为0时，mappedPlanarMaskX一定为0
  bool sure_planarityX = mappedPlanarMaskX || !planarPossibleX;
  //这个变量本身没有明确的意义，要配合下面for循环的第一个continue来理解！！！
  //第一个continue对于通过sure_planarityX来判断后续的bitIsOneX至关重要
  //maskedOccupancy由mappedPlanarMaskX、Y、Z组成，以mappedPlanarMaskX为例，
  //如果其为正，说明一定存在平面，而当前节点如果不在对应的平面上，则直接被跳过
  //mappedPlanarMaskX为正，因此留下的节点必定在平面上
  //因此sure_planarityX进入到bitIsOneX判定时，其为true时的情况有两种
  //1，x维度必定不存在平面；2，x维度必定存在平面且当前待编码的子节点就处于这个x平面
  //这两种情况都可以在子节点所处平面未占据的节点数为3时直接确定最后一个为占据


  bool sure_planarityY = mappedPlanarMaskY || !planarPossibleY;
  bool sure_planarityZ = mappedPlanarMaskZ || !planarPossibleZ;
  //由于这个mask指示的非占据情况，因此三者相或即可得出能确定的非占据的所以情况
  int maskedOccupancy =//注意取值为0时说明三个维度上都不是平面
    mappedPlanarMaskX | mappedPlanarMaskY | mappedPlanarMaskZ;

  int coded0X[2] = {0, 0};//存储各维度平面中高平面与低平面不被占据的子节点个数
  int coded0Y[2] = {0, 0};
  int coded0Z[2] = {0, 0};
  if (maskedOccupancy) {//通过平面提供的信息，只要不是全部非平面，可以通过平面信息得知哪些子节点一定非占据
    for (int i = 0; i < 8; i++) {//关心的只是个数，所以不用经过下面的kOccBitCodingOrder[i]的映射
      if ((maskedOccupancy >> i) & 1) {//对应位被确定为非占据，则子节点所在各维度上对应平面中非占据个数加1
        coded0X[(mappedFixedMaskX0 >> i) & 1]++;//mappedFixedMaskX0在映射前为0xf0
        coded0Y[(mappedFixedMaskY0 >> i) & 1]++;//mappedFixedMaskY0在映射前为0xcc
        coded0Z[(mappedFixedMaskZ0 >> i) & 1]++;//mappedFixedMaskZ0在映射前为0xaa
      }
    }
  }

  // NB: it is impossible for pattern to be 0 (handled in Z case).
  for (int i = 0; i < 8; i++) {
    int bitIdx = kOccBitCodingOrder[i];//调整编码顺序以适应后续的旋转、屏蔽
    if ((maskedOccupancy >> bitIdx) & 1)//该子节点不占据，由于解码端同样可以推导，因此不用编码，直接跳过
      continue;
	//另外，这一步对于通过sure_planarityX来判断后续的bitIsOneX至关重要
	//maskedOccupancy由mappedPlanarMaskX、Y、Z组成，以mappedPlanarMaskX为例，
	//如果其为正，说明一定存在平面，而当前节点如果不在对应的平面上，则直接被跳过
	//mappedPlanarMaskX为正时，因此留下的节点必定在平面上
	//因此前面的sure_planarityX为true时的情况有两种
	//1，x维度必定不是个平面；2，x维度是平面且当前子节点处于x平面


    int idx;
    if (i < 4) {
      idx = ((neighPatternR1 - 1) << i) + partialOccupancy + i + 1;//子节点0，1，2，3的六近邻都是9种状态
    } else if (i < 6) {
      idx = ((neighPattern5 - 1) << i) + partialOccupancy + i + 1;//子节点4，5的六近邻都是5种状态
    } else if (i == 6) {
      idx = ((neighPattern3 - 1) << i) + partialOccupancy + i + 1;//子节点6的六近邻状态是3种
    } else if (i == 7) {
      idx = partialOccupancy + i + 1;//子节点7不考虑六近邻
    } else {
      // work around clang -Wsometimes-uninitialized fault
      break;
    }

    int bitAdjGt0 = (mappedOccAdjGt0 >> bitIdx) & 1;
    int bitAdjGt1 = (mappedOccAdjGt1 >> bitIdx) & 1;
    int bitAdjUnocc = (mappedOccAdjUnocc >> bitIdx) & 1;

    int numAdj = bitAdjGt0 + bitAdjGt1;//只有三种情况：0，无相邻子节点；1：只有一个相邻子节点；2，两个及以上子节点
    int idxAdj = bitAdjUnocc + 2 * numAdj;//最大值为5
    if (i > 4) {
      static const int8_t kCtxIdxAdjReduc567[6] = {0, 0, 1, 2, 3, 3};
      idxAdj = kCtxIdxAdjReduc567[idxAdj];
    }

    int ctxIdxMapIdx = 3 * idxAdj;//取值范围为0到17
    if (!maskedOccupancy) {  // 三个平面都判定为非平面
      int bitIsPredicted = (mappedOccIsPredicted >> bitIdx) & 1;
      int bitPrediction = (mappedOccPrediction >> bitIdx) & 1;
      ctxIdxMapIdx = 3 * idxAdj + bitIsPredicted + bitPrediction;//
    }

    // NB: if firt 7 bits are 0, then the last is implicitly 1.
    // masking for planar is here
    int mask0X = (mappedFixedMaskX0 >> bitIdx) & 1;//当前子节点属于x低平面还是高平面，0则说明其位于低平面，1位于高平面
	//如果x是平面，那么x的低平面与高平面一定只有一个平面有点
    bool bitIsOneX = (sure_planarityX && coded0X[mask0X] >= 3)//子节点所在的0-3或者4-7最多1个被占据，且当前子节点位于平面或者当前节点不是平面
      || (coded0X[0] + coded0X[1] >= 7);//七个子节点都不占据的时候

    int mask0Y = (mappedFixedMaskY0 >> bitIdx) & 1;
    bool bitIsOneY = (sure_planarityY && coded0Y[mask0Y] >= 3)
      || (coded0Y[0] + coded0Y[1] >= 7);

    int mask0Z = (mappedFixedMaskZ0 >> bitIdx) & 1;
    bool bitIsOneZ = (sure_planarityZ && coded0Z[mask0Z] >= 3)
      || (coded0Z[0] + coded0Z[1] >= 7);
	//如果bitIsOneX || bitIsOneY || bitIsOneZ，则说明最后一位bit一定为1，
    int bit = (mappedOccupancy >> bitIdx) & 1;
    if (!(bitIsOneX || bitIsOneY || bitIsOneZ)) {
      auto& ctxIdxMap = _ctxIdxMaps[ctxIdxMapIdx];
      int ctxIdx = ctxIdxMap.evolve(bit, &ctxIdxMap[i][idx]);//注意此处与之前所看的上下文编码不同，其得到的ctx的值都是初始化之后再更新的，idx只是个索引，
	  //之前的上下文编码以及平面模式中的上下文编码是固定索引，更新概率；而此处每一个情况存的是索引值，更新的也是这个索引值。这有点类似265中固定编码器，更新索引
      _arithmeticEncoder->encode(bit, _ctxOccupancy[ctxIdx]);

      if (!bit) {//之前是通过平面信息来推断，并不能包含所有的情况，因此还需根据实际情况来继续判断
        coded0X[mask0X]++;
        coded0Y[mask0Y]++;
        coded0Z[mask0Z]++;
      }
    }

    partialOccupancy |= bit << i;
  }
}

//-------------------------------------------------------------------------

void
GeometryOctreeEncoder::encodeOccupancyBitwise(
  int neighPattern,
  int mappedOccupancy,
  int mappedOccIsPredicted,
  int mappedOccPrediction,
  int mappedOccAdjGt0,
  int mappedOccAdjGt1,
  int mappedOccAdjUnocc,
  int mappedPlanarMaskX,
  int mappedFixedMaskX0,
  bool planarPossibleX,
  int mappedPlanarMaskY,
  int mappedFixedMaskY0,
  bool planarPossibleY,
  int mappedPlanarMaskZ,
  int mappedFixedMaskZ0,
  bool planarPossibleZ)

{
  if (neighPattern == 0) {
    encodeOccupancyNeighZ(
      mappedOccupancy, mappedOccIsPredicted, mappedOccPrediction,
      mappedOccAdjGt0, mappedOccAdjGt1, mappedOccAdjUnocc, mappedPlanarMaskX,
      mappedFixedMaskX0, planarPossibleX, mappedPlanarMaskY, mappedFixedMaskY0,
      planarPossibleY, mappedPlanarMaskZ, mappedFixedMaskZ0, planarPossibleZ);
    return;
  }

  encodeOccupancyNeighNZ(
    neighPattern, mappedOccupancy, mappedOccIsPredicted, mappedOccPrediction,
    mappedOccAdjGt0, mappedOccAdjGt1, mappedOccAdjUnocc, mappedPlanarMaskX,
    mappedFixedMaskX0, planarPossibleX, mappedPlanarMaskY, mappedFixedMaskY0,
    planarPossibleY, mappedPlanarMaskZ, mappedFixedMaskZ0, planarPossibleZ);
}

//-------------------------------------------------------------------------

void
GeometryOctreeEncoder::encodeOccupancyBytewise(
  int neighPattern, int mappedOccupancy)
{
  // code occupancy using the neighbour configuration context
  // with reduction from 64 states to 10 (or 6).
  int neighPatternR1 = _neighPattern64toR1[neighPattern];
  auto& bytewiseCoder = _bytewiseOccupancyCoder[neighPatternR1];
  bytewiseCoder.encode(mappedOccupancy, _arithmeticEncoder);
}

//-------------------------------------------------------------------------
// decode node occupancy bits
//

void
GeometryOctreeEncoder::encodeOccupancy(
  const GeometryNeighPattern& gnp,
  int occupancy,
  int occupancyIsPred,
  int occupancyPred,
  int planarMaskX,
  int planarMaskY,
  int planarMaskZ,
  bool planarPossibleX,
  bool planarPossibleY,
  bool planarPossibleZ)
{
  // 3 planars => single child and we know its position
  if (planarMaskX && planarMaskY && planarMaskZ)//说明是三个平面，则可以根据planePos
    return;//需要注意的是，如果为1则一定是三平面，但是如果为0，则不一定不是三平面/NO=1

  if (gnp.neighPattern == 0) {//资格判定，保证singleNode的正确率
    bool singleChild = !popcntGt1(occupancy);//判断NO=1，即是否为真的三平面
    if (planarPossibleX && planarPossibleY && planarPossibleZ) {//这也是一个资格判断,能排除一部分一定不是singlechild的情况，仅此而已
      _arithmeticEncoder->encode(singleChild, _ctxSingleChild);
    }

    if (singleChild) {//是singleNode但是上面没判定出来
      // no siblings => encode index = (z,y,x) not 8bit pattern
      // if mask is not zero, then planar, then child z known from plane index
      if (!planarMaskZ)//依次将planarMask没有判定出来的继续判定
        _arithmeticEncoder->encode(!!(occupancy & 0xaa));//以此判断低平面还是高平面

      if (!planarMaskY)	
        _arithmeticEncoder->encode(!!(occupancy & 0xcc));

      if (!planarMaskX)
        _arithmeticEncoder->encode(!!(occupancy & 0xf0));

      return;
    }
  }

  // at least two child nodes occupied and two planars => we know the occupancy
  if (gnp.neighPattern == 0) {//这个限定是必须的，因为如果不和singleChild保持一致的话，不能确保一定不存在NO=1即实质上为三个平面的情况
    if (planarMaskX && planarMaskY)//判断两个方向的平面是否重合，如果重合，配合planePosition即可确定occupancy
      return;
    if (planarMaskY && planarMaskZ)
      return;
    if (planarMaskX && planarMaskZ)
      return;
  }

  auto neighPattern = gnp.neighPattern;
  //接下来的旋转变换都是基于neighPattern来决定的，以保持一致性
  auto mapOcc = mapGeometryOccupancy(occupancy, neighPattern);//根据六近邻的几何旋转不变性来映射occupancy
  auto mapOccIsP = mapGeometryOccupancy(occupancyIsPred, neighPattern);//指示是否帧内预测
  auto mapOccP = mapGeometryOccupancy(occupancyPred, neighPattern);//指示预测的结果
  auto mapAdjGt0 = mapGeometryOccupancy(gnp.adjacencyGt0, neighPattern);//指示是否有大于0个共面的子节点
  auto mapAdjGt1 = mapGeometryOccupancy(gnp.adjacencyGt1, neighPattern);//指示是否有大于1个共面的子节点
  auto mapAdjUnocc = mapGeometryOccupancy(gnp.adjacencyUnocc, neighPattern);//指示未占据的共面子节点个数

  auto mapPlanarMaskX = mapGeometryOccupancy(planarMaskX, neighPattern);//mask指示是否为非平面以及平面位置，其八位，0代表非平面，0xfc代表高平面
  auto mapPlanarMaskY = mapGeometryOccupancy(planarMaskY, neighPattern);
  auto mapPlanarMaskZ = mapGeometryOccupancy(planarMaskZ, neighPattern);

  auto mapFixedMaskX0 = mapGeometryOccupancy(0xf0, neighPattern);//固定的mask，用来指示各维度的高平面
  auto mapFixedMaskY0 = mapGeometryOccupancy(0xcc, neighPattern);
  auto mapFixedMaskZ0 = mapGeometryOccupancy(0xaa, neighPattern);

  if (_useBitwiseOccupancyCoder)
    encodeOccupancyBitwise(
      neighPattern, mapOcc, mapOccIsP, mapOccP, mapAdjGt0, mapAdjGt1,
      mapAdjUnocc, mapPlanarMaskX, mapFixedMaskX0, planarPossibleX,
      mapPlanarMaskY, mapFixedMaskY0, planarPossibleY, mapPlanarMaskZ,
      mapFixedMaskZ0, planarPossibleZ);//
  else
    encodeOccupancyBytewise(neighPattern, mapOcc);
}

//-------------------------------------------------------------------------
// Encode part of the position of two unordred points  point in a given volume.
void
GeometryOctreeEncoder::encodeOrdered2ptPrefix(
  const point_t points[2], Vec3<bool> directIdcm, Vec3<int>& nodeSizeLog2)
{
  if (nodeSizeLog2[0] >= 1 && directIdcm[0]) {
    bool sameBit = true;
    int ctxIdx = 0;
    while (nodeSizeLog2[0] && sameBit) {
      nodeSizeLog2[0]--;
      int mask = 1 << nodeSizeLog2[0];
      auto bit0 = !!(points[0][0] & mask);
      auto bit1 = !!(points[1][0] & mask);
      sameBit = bit0 == bit1;

      _arithmeticEncoder->encode(sameBit, _ctxSameBitHighx[ctxIdx]);
      ctxIdx = std::min(4, ctxIdx + 1);
      if (sameBit)
        _arithmeticEncoder->encode(bit0);
    }
  }

  if (nodeSizeLog2[1] >= 1 && directIdcm[1]) {
    bool sameX = !directIdcm[0] || points[0][0] == points[1][0];
    bool sameBit = true;
    int ctxIdx = 0;
    while (nodeSizeLog2[1] && sameBit) {
      nodeSizeLog2[1]--;
      int mask = 1 << nodeSizeLog2[1];
      auto bit0 = !!(points[0][1] & mask);
      auto bit1 = !!(points[1][1] & mask);
      sameBit = bit0 == bit1;

      _arithmeticEncoder->encode(sameBit, _ctxSameBitHighy[ctxIdx]);
      ctxIdx = std::min(4, ctxIdx + 1);
      if (!(sameX && !sameBit))
        _arithmeticEncoder->encode(bit0);
    }
  }

  if (nodeSizeLog2[2] >= 1 && directIdcm[2]) {
    bool sameBit = true;
    bool sameXy = (!directIdcm[0] || points[0][0] == points[1][0])
      && (!directIdcm[1] || points[0][1] == points[1][1]);
    int ctxIdx = 0;
    while (nodeSizeLog2[2] && sameBit) {
      nodeSizeLog2[2]--;
      int mask = 1 << nodeSizeLog2[2];
      auto bit0 = !!(points[0][2] & mask);
      auto bit1 = !!(points[1][2] & mask);
      sameBit = bit0 == bit1;

      _arithmeticEncoder->encode(sameBit, _ctxSameBitHighz[ctxIdx]);
      ctxIdx = std::min(4, ctxIdx + 1);
      if (!(sameXy && !sameBit))
        _arithmeticEncoder->encode(bit0);
    }
  }
}

//-------------------------------------------------------------------------
// Encode a position of a point in a given volume.
void
GeometryOctreeEncoder::encodePointPosition(
  const Vec3<int>& nodeSizeLog2AfterPlanar, const Vec3<int32_t>& pos)
{
  for (int k = 0; k < 3; k++) {
    if (nodeSizeLog2AfterPlanar[k] <= 0)
      continue;

    for (int mask = 1 << (nodeSizeLog2AfterPlanar[k] - 1); mask; mask >>= 1) {
      _arithmeticEncoder->encode(!!(pos[k] & mask));
    }
  }
}

//-------------------------------------------------------------------------
// Encode a position of a point in a given volume, using elevation angle prior

void
GeometryOctreeEncoder::encodePointPositionAngular(
  const Vec3<int>& nodeSizeLog2,
  const Vec3<int>& nodeSizeLog2AfterUnordered,
  const Vec3<int32_t>& pos,
  const PCCOctree3Node& child,
  const OctreeNodePlanar& planar,
  const Vec3<int>& headPos,
  const int* zLaser,
  const int* thetaLaser,
  int numLasers)
{
  Vec3<int> posXyz = {(child.pos[0] << nodeSizeLog2[0]) - headPos[0],
                      (child.pos[1] << nodeSizeLog2[1]) - headPos[1],
                      (child.pos[2] << nodeSizeLog2[2]) - headPos[2]};

  // -- PHI --
  // code x or y directly and compute phi of node
  bool codeXorY = std::abs(posXyz[0]) <= std::abs(posXyz[1]);
  if (codeXorY) {  // direct code y
    if (nodeSizeLog2AfterUnordered[1])
      for (int mask = 1 << (nodeSizeLog2AfterUnordered[1] - 1); mask;
           mask >>= 1)
        _arithmeticEncoder->encode(!!(pos[1] & mask));

    posXyz[1] = pos[1] - headPos[1];
    if (planar.planarMode & 1) {
      int mask = 1 << (nodeSizeLog2[0] - 1);
      if (pos[0] & mask)
        posXyz[0] += mask;
    }
  } else {  //direct code x
    if (nodeSizeLog2AfterUnordered[0])
      for (int mask = 1 << (nodeSizeLog2AfterUnordered[0] - 1); mask;
           mask >>= 1)
        _arithmeticEncoder->encode(!!(pos[0] & mask));

    posXyz[0] = pos[0] - headPos[0];
    if (planar.planarMode & 2) {
      int mask = 1 << (nodeSizeLog2[1] - 1);
      if (pos[1] & mask)
        posXyz[1] += mask;
    }
  }

  // Laser
  int laserNode = int(child.laserIndex);

  point_t posPointLidar =
    point_t(pos[0] - headPos[0], pos[1] - headPos[1], pos[2] - headPos[2]);
  int laserIndex = findLaser(posPointLidar, thetaLaser, numLasers);
  encodeThetaRes(laserIndex - laserNode);

  // find predictor
  int phiNode = iatan2(posXyz[1], posXyz[0]);
  int predPhi = _phiBuffer[laserIndex];
  if (predPhi == 0x80000000)
    predPhi = phiNode;

  // elementary shift predictor
  int nShift =
    ((predPhi - phiNode) * _phiZi.invDelta(laserIndex) + 536870912) >> 30;
  predPhi -= _phiZi.delta(laserIndex) * nShift;

  // choose x or y
  int* posXY = codeXorY ? &posXyz[0] : &posXyz[1];
  int idx = codeXorY ? 0 : 1;

  // azimuthal code x or y
  int mask2 = codeXorY ? (nodeSizeLog2AfterUnordered[0] > 0
                            ? 1 << (nodeSizeLog2AfterUnordered[0] - 1)
                            : 0)
                       : (nodeSizeLog2AfterUnordered[1] > 0
                            ? 1 << (nodeSizeLog2AfterUnordered[1] - 1)
                            : 0);
  for (; mask2; mask2 >>= 1) {
    // angles left and right
    int phiR = codeXorY ? iatan2(posXyz[1], posXyz[0] + mask2)
                        : iatan2(posXyz[1] + mask2, posXyz[0]);
    int phiL = phiNode;

    // ctx azimutal
    int angleL = phiL - predPhi;
    int angleR = phiR - predPhi;
    int contextAnglePhi =
      (angleL >= 0 && angleR >= 0) || (angleL < 0 && angleR < 0) ? 2 : 0;
    angleL = std::abs(angleL);
    angleR = std::abs(angleR);
    if (angleL > angleR) {
      contextAnglePhi++;
      int temp = angleL;
      angleL = angleR;
      angleR = temp;
    }
    if (angleR > (angleL << 1))
      contextAnglePhi += 4;

    // entropy coding
    int bit = !!(pos[idx] & mask2);
    _arithmeticEncoder->encode(
      bit, _ctxPlanarPlaneLastIndexAngularPhiIDCM[contextAnglePhi]);
    if (bit) {
      *posXY += mask2;
      phiNode = phiR;
      predPhi = _phiBuffer[laserIndex];
      if (predPhi == 0x80000000)
        predPhi = phiNode;

      // elementary shift predictor
      int nShift =
        ((predPhi - phiNode) * _phiZi.invDelta(laserIndex) + 536870912) >> 30;
      predPhi -= _phiZi.delta(laserIndex) * nShift;
    }
  }

  _phiBuffer[laserIndex] = phiNode;

  // -- THETA --
  int maskz = nodeSizeLog2AfterUnordered[2] > 0
    ? 1 << (nodeSizeLog2AfterUnordered[2] - 1)
    : 0;
  if (!maskz)
    return;

  if (planar.planarMode & 4) {
    int mask = 1 << (nodeSizeLog2[2] - 1);
    if (pos[2] & mask)
      posXyz[2] += mask;
  }

  // Since x and y are known,
  // r is known too and does not depend on the bit for z
  uint64_t xLidar = (int64_t(posXyz[0]) << 8) - 128;
  uint64_t yLidar = (int64_t(posXyz[1]) << 8) - 128;
  uint64_t r2 = xLidar * xLidar + yLidar * yLidar;
  int64_t rInv = irsqrt(r2);

  // code z
  int64_t hr = zLaser[laserIndex] * rInv;
  int fixedThetaLaser =
    thetaLaser[laserIndex] + int(hr >= 0 ? -(hr >> 17) : ((-hr) >> 17));

  int zShift = (rInv << nodeSizeLog2AfterUnordered[2]) >> 18;
  for (; maskz; maskz >>= 1, zShift >>= 1) {
    // determine non-corrected theta
    int64_t zLidar = ((posXyz[2] + maskz) << 1) - 1;
    int64_t theta = zLidar * rInv;
    int theta32 = theta >= 0 ? theta >> 15 : -((-theta) >> 15);
    int thetaLaserDelta = fixedThetaLaser - theta32;

    int thetaLaserDeltaBot = thetaLaserDelta + zShift;
    int thetaLaserDeltaTop = thetaLaserDelta - zShift;
    int contextAngle = thetaLaserDelta >= 0 ? 0 : 1;
    if (thetaLaserDeltaTop >= 0)
      contextAngle += 2;
    else if (thetaLaserDeltaBot < 0)
      contextAngle += 2;

    int bit = !!(pos[2] & maskz);
    _arithmeticEncoder->encode(
      bit, _ctxPlanarPlaneLastIndexAngularIdcm[contextAngle]);
    if (bit)
      posXyz[2] += maskz;
  }
}

//-------------------------------------------------------------------------
// Direct coding of position of points in node (early tree termination).
void
GeometryOctreeEncoder::encodeQpOffset(int dqp)
{
  _arithmeticEncoder->encode(dqp == 0, _ctxQpOffsetIsZero);
  if (dqp == 0) {
    return;
  }
  _arithmeticEncoder->encode(dqp > 0, _ctxQpOffsetSign);
  _arithmeticEncoder->encodeExpGolomb(abs(dqp) - 1, 0, _ctxQpOffsetAbsEgl);
}

//-------------------------------------------------------------------------

template<typename It>
void
setNodeQpsUniform(
  Vec3<int> nodeSizeLog2,
  int qp,
  int geom_qp_multiplier_log2,
  It nodesBegin,
  It nodesEnd)
{
  // Conformance: limit the qp such that it cannot overquantize the node
  qp = std::min(qp, nodeSizeLog2.min() * 8);
  assert(qp % (1 << geom_qp_multiplier_log2) == 0);

  for (auto it = nodesBegin; it != nodesEnd; ++it)
    it->qp = qp;
}

//-------------------------------------------------------------------------
// Sets QP randomly

template<typename It>
void
setNodeQpsRandom(
  Vec3<int> nodeSizeLog2,
  int /* qp */,
  int geom_qp_multiplier_log2,
  It nodesBegin,
  It nodesEnd)
{
  // Conformance: limit the qp such that it cannot overquantize the node
  int maxQp = nodeSizeLog2.min() * 8;

  int seed = getenv("SEED") ? atoi(getenv("SEED")) : 0;
  static std::minstd_rand gen(seed);
  std::uniform_int_distribution<> uniform(0, maxQp);

  // pick a random qp, avoiding unrepresentable values
  for (auto it = nodesBegin; it != nodesEnd; ++it)
    it->qp = uniform(gen) & (~0 << geom_qp_multiplier_log2);
}

//-------------------------------------------------------------------------
// determine delta qp for each node based on the point density

template<typename It>
void
setNodeQpsByDensity(
  Vec3<int> nodeSizeLog2,
  int baseQp,
  int geom_qp_multiplier_log2,
  It nodesBegin,
  It nodesEnd)
{
  // Conformance: limit the qp such that it cannot overquantize the node
  int maxQp = nodeSizeLog2.min() * 8;
  int lowQp = PCCClip(baseQp - 8, 0, maxQp);
  int mediumQp = std::min(baseQp, maxQp);
  int highQp = std::min(baseQp + 8, maxQp);

  // NB: node.qp always uses a step size doubling interval of 8 QPs.
  //     the chosen QPs (and conformance limit) must respect the qp multiplier
  assert(lowQp % (1 << geom_qp_multiplier_log2) == 0);
  assert(mediumQp % (1 << geom_qp_multiplier_log2) == 0);
  assert(highQp % (1 << geom_qp_multiplier_log2) == 0);

  std::vector<int> numPointsInNode;
  std::vector<double> cum_prob;
  int32_t numPointsInLvl = 0;
  for (auto it = nodesBegin; it != nodesEnd; ++it) {
    numPointsInNode.push_back(it->end - it->start);
    numPointsInLvl += it->end - it->start;
  }
  std::sort(numPointsInNode.begin(), numPointsInNode.end());
  double cc = 0;
  for (auto num : numPointsInNode) {
    cc += num;
    cum_prob.push_back(cc / numPointsInLvl);
  }
  int th1 = -1, th2 = -1;
  for (int i = 0; i < cum_prob.size(); i++) {
    if (th1 == -1 && cum_prob[i] > 0.05) {
      th1 = numPointsInNode[i];
    } else if (th2 == -1 && cum_prob[i] > 0.6)
      th2 = numPointsInNode[i];
  }
  for (auto it = nodesBegin; it != nodesEnd; ++it) {
    if (it->end - it->start < th1) {
      it->qp = highQp;
    } else if (it->end - it->start < th2)
      it->qp = mediumQp;
    else
      it->qp = lowQp;
  }
}

//-------------------------------------------------------------------------

template<typename It>
void
calculateNodeQps(
  OctreeEncOpts::QpMethod method,
  Vec3<int> nodeSizeLog2,
  int baseQp,
  int geom_qp_multiplier_log2,
  It nodesBegin,
  It nodesEnd)
{
  auto fn = &setNodeQpsUniform<It>;

  switch (method) {
    using Method = OctreeEncOpts::QpMethod;
  default:
  case Method::kUniform: fn = &setNodeQpsUniform<It>; break;
  case Method::kRandom: fn = &setNodeQpsRandom<It>; break;
  case Method::kByDensity: fn = &setNodeQpsByDensity<It>; break;
  }

  fn(nodeSizeLog2, baseQp, geom_qp_multiplier_log2, nodesBegin, nodesEnd);
}

//-------------------------------------------------------------------------

void
geometryQuantization(
  PCCPointSet3& pointCloud, PCCOctree3Node& node, Vec3<int> nodeSizeLog2)
{
  QuantizerGeom quantizer = QuantizerGeom(node.qp);
  int qpShift = QuantizerGeom::qpShift(node.qp);

  for (int k = 0; k < 3; k++) {
    int quantBitsMask = (1 << nodeSizeLog2[k]) - 1;
    int32_t clipMax = quantBitsMask >> qpShift;

    for (int i = node.start; i < node.end; i++) {
      int32_t pos = int32_t(pointCloud[i][k]);
      int32_t quantPos = quantizer.quantize(pos & quantBitsMask);
      quantPos = PCCClip(quantPos, 0, clipMax);

      // NB: this representation is: |ppppppqqq00|, which, except for
      // the zero padding, is the same as the decoder.
      pointCloud[i][k] = (pos & ~quantBitsMask) | (quantPos << qpShift);
    }
  }
}

//-------------------------------------------------------------------------

void
geometryScale(
  PCCPointSet3& pointCloud, PCCOctree3Node& node, Vec3<int> quantNodeSizeLog2)
{
  QuantizerGeom quantizer = QuantizerGeom(node.qp);
  int qpShift = QuantizerGeom::qpShift(node.qp);

  for (int k = 0; k < 3; k++) {
    int quantBitsMask = (1 << quantNodeSizeLog2[k]) - 1;
    for (int i = node.start; i < node.end; i++) {
      int pos = pointCloud[i][k];
      int lowPart = (pos & quantBitsMask) >> qpShift;
      int lowPartScaled = PCCClip(quantizer.scale(lowPart), 0, quantBitsMask);
      int highPartScaled = pos & ~quantBitsMask;
      pointCloud[i][k] = highPartScaled | lowPartScaled;
    }
  }
}

//-------------------------------------------------------------------------

void
checkDuplicatePoints(
  PCCPointSet3& pointCloud,
  PCCOctree3Node& node,
  std::vector<int>& pointIdxToDmIdx)
{
  auto first = PCCPointSet3::iterator(&pointCloud, node.start);
  auto last = PCCPointSet3::iterator(&pointCloud, node.end);

  std::set<Vec3<int32_t>> uniquePointsSet;
  for (auto i = first; i != last;) {
    if (uniquePointsSet.find(**i) == uniquePointsSet.end()) {
      uniquePointsSet.insert(**i);
      i++;
    } else {
      std::iter_swap(i, last - 1);
      last--;
      pointIdxToDmIdx[--node.end] = -2;  // mark as duplicate
    }
  }
}

//-------------------------------------------------------------------------

DirectMode
canEncodeDirectPosition(
  bool geom_unique_points_flag,
  const PCCOctree3Node& node,
  const PCCPointSet3& pointCloud)
{
  int numPoints = node.end - node.start;
  // Check for duplicated points only if there are less than 10.
  // NB: this limit is rather arbitrary
  if (numPoints > 10)//点数大于10直接不予考虑
    return DirectMode::kUnavailable;

  bool allPointsAreEqual = numPoints > 1 && !geom_unique_points_flag;//不只一个点且都是重复点
  for (auto idx = node.start + 1; allPointsAreEqual && idx < node.end; idx++)
    allPointsAreEqual &= pointCloud[node.start] == pointCloud[idx];

  if (allPointsAreEqual)
    return DirectMode::kAllPointSame;

  if (numPoints > MAX_NUM_DM_LEAF_POINTS)
    return DirectMode::kUnavailable;

  return DirectMode::kTwoPoints;//两个或一个点
}

//-------------------------------------------------------------------------

void
GeometryOctreeEncoder::encodeIsIdcm(DirectMode mode)
{
  bool isIdcm = mode != DirectMode::kUnavailable;
  _arithmeticEncoder->encode(isIdcm, _ctxBlockSkipTh);
}

//-------------------------------------------------------------------------

void
GeometryOctreeEncoder::encodeDirectPosition(
  DirectMode mode,
  bool geom_unique_points_flag,
  bool joint_2pt_idcm_enabled_flag,
  const Vec3<int>& effectiveNodeSizeLog2,
  int shiftBits,
  PCCOctree3Node& node,
  OctreeNodePlanar& planar,
  PCCPointSet3& pointCloud,
  bool angularIdcm,
  const Vec3<int>& headPos,
  const int* zLaser,
  const int* thetaLaser,
  int numLasers)
{
  int numPoints = node.end - node.start;

  switch (mode) {
  case DirectMode::kUnavailable: return;

  case DirectMode::kTwoPoints://有两个点或者一个点
    _arithmeticEncoder->encode(numPoints > 1, _ctxNumIdcmPointsGt1);//两种情况下都是先编码点数
    if (!geom_unique_points_flag && numPoints == 1)
      _arithmeticEncoder->encode(numPoints == 1, _ctxSinglePointPerBlock);
    break;

  case DirectMode::kAllPointSame://重复点
    _arithmeticEncoder->encode(0, _ctxNumIdcmPointsGt1);
    _arithmeticEncoder->encode(0, _ctxSinglePointPerBlock);
    _arithmeticEncoder->encode(numPoints == 2, _ctxSingleIdcmDupPoint);
    if (numPoints > 2)
      _arithmeticEncoder->encodeExpGolomb(
        numPoints - 3, 0, _ctxPointCountPerBlock);

    // only one actual psoition to code
    numPoints = 1;
  }

  // if the points have been quantised, the following representation is used
  // for point cloud positions:
  //          |---| = nodeSizeLog2 (example)
  //   ppppppqqqq00 = cloud[ptidx]
  //          |-|   = effectiveNodeSizeLog2 (example)
  // where p are unquantised bits, qqq are quantised bits, and 0 are zero bits.
  // nodeSizeLog2 is the size of the current node prior to quantisation.
  // effectiveNodeSizeLog2 is the size of the node after quantisation.
  //
  // NB: while nodeSizeLog2 may be used to access the current position bit
  //     in both quantised and unquantised forms, effectiveNodeSizeLog2 cannot
  //     without taking into account the padding.
  //
  // NB: this contrasts with node.pos, which contains the previously coded
  //     position bits ("ppppppq" in the above example) without any padding.
  //
  // When coding the direct mode, the zero padding is removed to permit
  // indexing by the effective node size instead.
  Vec3<int> points[2];
  for (int i = 0; i < numPoints; i++)
    points[i] = pointCloud[node.start + i] >> shiftBits;

  // update node size after planar
  Vec3<int> nodeSizeLog2Rem = effectiveNodeSizeLog2;
  for (int k = 0; k < 3; k++) {
    if (nodeSizeLog2Rem[k] > 0 && (planar.planarMode & (1 << k)))
      nodeSizeLog2Rem[k]--;
  }

  // Indicates which components are directly coded, or coded using angular
  // contextualisation.
  Vec3<bool> directIdcm = !angularIdcm;
  point_t posNodeLidar;

  if (angularIdcm) {
    posNodeLidar = node.pos;
    // todo(df): this should be fixed to take quantisation into account
    // inorder to compare with headPos, shift by shiftBits and changed in
    // the decoder too.
    for (int k = 0; k < 3; k++)
      posNodeLidar[k] <<= effectiveNodeSizeLog2[k];
    posNodeLidar -= headPos;

    bool codeXorY = std::abs(posNodeLidar[0]) <= std::abs(posNodeLidar[1]);
    directIdcm.x() = !codeXorY;
    directIdcm.y() = codeXorY;
  }

  // Jointly code two points
  if (numPoints == 2 && joint_2pt_idcm_enabled_flag) {
    // Apply an implicit ordering to the two points, considering only the
    // directly coded axes
    if (times(points[1], directIdcm) < times(points[0], directIdcm)) {//对两点进行排序
      std::swap(points[0], points[1]);                                 
      pointCloud.swapPoints(node.start, node.start + 1);//点云中两点也进行排序，排序十分重要，决定着能否在坐标不等时节省一个bit
    }

    encodeOrdered2ptPrefix(points, directIdcm, nodeSizeLog2Rem);
  }

  if (angularIdcm) {
    for (int k = 0; k < 3; ++k) {
      int mask = (1 << effectiveNodeSizeLog2[k]);
      for (int i = 0; i < effectiveNodeSizeLog2[k] - nodeSizeLog2Rem[k]; ++i) {
        mask >>= 1;
        if (points[0][k] & mask)
          posNodeLidar[k] += mask;
      }
      mask >>= 1;
      posNodeLidar[k] += mask;
    }
    node.laserIndex = findLaser(posNodeLidar, thetaLaser, numLasers);
  }

  // code points after planar
  for (auto idx = 0; idx < numPoints; idx++) {
    if (angularIdcm) {
      encodePointPositionAngular(
        effectiveNodeSizeLog2, nodeSizeLog2Rem, points[idx], node, planar,
        headPos, zLaser, thetaLaser, numLasers);
    } else
      encodePointPosition(nodeSizeLog2Rem, points[idx]);
  }
}

//-------------------------------------------------------------------------

void
GeometryOctreeEncoder::encodeThetaRes(int thetaRes)
{
  _arithmeticEncoder->encode(thetaRes == 0 ? 1 : 0, _ctxThetaResIsZero);

  if (thetaRes) {
    _arithmeticEncoder->encode(thetaRes > 0 ? 1 : 0, _ctxThetaResSign);
    int absThetaRes = std::abs(thetaRes);
    _arithmeticEncoder->encode(absThetaRes == 1 ? 1 : 0, _ctxThetaResIsOne);
    if (absThetaRes >= 2)
      _arithmeticEncoder->encode(absThetaRes == 2 ? 1 : 0, _ctxThetaResIsTwo);
    if (absThetaRes >= 3)
      _arithmeticEncoder->encodeExpGolomb(absThetaRes - 3, 1, _ctxThetaResExp);
  }
}

//-------------------------------------------------------------------------

void
encodeGeometryOctree(
  const OctreeEncOpts& params,
  const GeometryParameterSet& gps,
  GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMem,
  std::vector<std::unique_ptr<EntropyEncoder>>& arithmeticEncoders,
  pcc::ringbuf<PCCOctree3Node>* nodesRemaining)
{
  auto arithmeticEncoderIt = arithmeticEncoders.begin();
  GeometryOctreeEncoder encoder(gps, gbh, ctxtMem, arithmeticEncoderIt->get());

  // saved state for use with parallel bistream coding.
  // the saved state is restored at the start of each parallel octree level
  std::unique_ptr<GeometryOctreeEncoder> savedState;

  // init main fifo
  //  -- worst case size is the last level containing every input poit
  //     and each point being isolated in the previous level.
  pcc::ringbuf<PCCOctree3Node> fifo(pointCloud.getPointCount() + 1);//用于盛放八叉树节点

  // push the first node
  fifo.emplace_back();
  PCCOctree3Node& node00 = fifo.back();
  node00.start = uint32_t(0);//根节点中所有点中，起始点序号自然为0
  node00.end = uint32_t(pointCloud.getPointCount());//根节点中所有点，末尾点序号为点的个数
  node00.pos = int32_t(0);//根节点坐标为000
  node00.numSiblingsPlus1 = 8;//兄弟节点+1
  node00.siblingOccupancy = 0;//当前节点及兄弟节点占位信息构成的8位occupancy
  node00.qp = 0;//量化步长
  node00.idcmEligible = 0;//IDCM资格，其在节点被刚刚创建时就开始判定，中间受平面编码的影响

  // map of pointCloud idx to DM idx, used to reorder the points
  // after coding.
  std::vector<int> pointIdxToDmIdx(int(pointCloud.getPointCount()), -1);
  int nextDmIdx = 0;

  // generate the list of the node size for each level in the tree
  auto lvlNodeSizeLog2 = mkQtBtNodeSizeList(gps, params.qtbt, gbh);//qtbt表，即根据k、m的值决定每一层八叉树节点的边长

  const int idcmThreshold = gps.geom_planar_mode_enabled_flag//IDCM的阈值，用于判断节点是否进行IDCM的资格
    ? gps.geom_planar_idcm_threshold * 127 * 127
    : 127 * 127 * 127;

  //  Lidar angles for planar prediction
  const int numLasers = gps.geom_angular_num_lidar_lasers();//雷达数目
  const int* thetaLaser = gps.geom_angular_theta_laser.data();//各雷达垂直角
  const int* zLaser = gps.geom_angular_z_laser.data();//雷达相对headposition的z偏移，用于计算修正垂直角
  const int* numPhi = gps.geom_angular_num_phi_per_turn.data();//水平角每一转的角度

  // Lidar position relative to slice origin
  auto headPos = gps.geomAngularOrigin - gbh.geomBoxOrigin;//雷达相对于slic原点的位置

  int deltaAngle = 128 << 18;//各雷达之间垂直角的差异，需大于一定数值
  for (int i = 0; i < numLasers - 1; i++) {
    int d = std::abs(thetaLaser[i] - thetaLaser[i + 1]);
    if (deltaAngle > d)
      deltaAngle = d;
  }

  MortonMap3D occupancyAtlas;//一个缓存区，用于获取邻居信息
  if (gps.neighbour_avail_boundary_log2) {//如果此参数为0，那么只有兄弟节点的占位信息可以被获取
    occupancyAtlas.resize(
      gps.adjacent_child_contextualization_enabled_flag,
      gps.neighbour_avail_boundary_log2);
    occupancyAtlas.clear();
  }

  // the minimum node size is ordinarily 2**0, but may be larger due to
  // early termination for trisoup.
  int minNodeSizeLog2 = gbh.trisoup_node_size_log2;//体素大小，即节点划分至该大小时停止划分

  // prune anything smaller than the minimum node size (these won't be coded)
  // NB: this must result in a cubic node at the end of the list
  // NB: precondition: root node size >= minNodeSizeLog2.
  lvlNodeSizeLog2.erase(//根据体素大小来剔除无需编码的节点
    std::remove_if(
      lvlNodeSizeLog2.begin(), lvlNodeSizeLog2.end(),
      [&](Vec3<int>& size) { return size < minNodeSizeLog2; }),
    lvlNodeSizeLog2.end());
  assert(lvlNodeSizeLog2.back() == minNodeSizeLog2);

  // append a dummy entry to the list so that depth+2 access is always valid
  lvlNodeSizeLog2.emplace_back(lvlNodeSizeLog2.back());

  // the termination depth of the octree phase
  // NB: the tree depth may be greater than the maxNodeSizeLog2 due to
  //     perverse qtbt splitting.
  // NB: by definition, the last two elements are minNodeSizeLog2
  int maxDepth = lvlNodeSizeLog2.size() - 2;

  // generate the qtbt splitting list
  //  - start at the leaf, and work up
  std::vector<int8_t> tree_lvl_partition_list;
  for (int lvl = 0; lvl < maxDepth; lvl++) {
    gbh.tree_lvl_coded_axis_list.push_back(
      ~nonSplitQtBtAxes(lvlNodeSizeLog2[lvl], lvlNodeSizeLog2[lvl + 1]));

    // Conformance: at least one axis must attempt to be coded at each level
    assert(gbh.tree_lvl_coded_axis_list.back() != 0);
  }

  // Determine the desired quantisation depth after qtbt is determined
  if (params.qpOffsetNodeSizeLog2 > 0) {
    // find the first level that matches the scaling node size
    for (int lvl = 0; lvl < maxDepth; lvl++) {
      if (lvlNodeSizeLog2[lvl].min() > params.qpOffsetNodeSizeLog2)
        continue;
      gbh.geom_octree_qp_offset_depth = lvl;
      break;
    }
  }

  // the node size where quantisation is performed
  Vec3<int> quantNodeSizeLog2 = 0;
  int idcmQp = 0;
  int sliceQp = gbh.sliceQp(gps);
  int numLvlsUntilQuantization = 0;//用于代替depth ，该参数为0时就要进行量化
  if (gps.geom_scaling_enabled_flag) {
    // if an invalid depth is set, use tree height instead
    if (gbh.geom_octree_qp_offset_depth < 0)
      gbh.geom_octree_qp_offset_depth = maxDepth;
    numLvlsUntilQuantization = gbh.geom_octree_qp_offset_depth + 1;
  }

  // The number of nodes to wait before updating the planar rate.
  // This is to match the prior behaviour where planar is updated once
  // per coded occupancy.
  int nodesBeforePlanarUpdate = 1;

  if (gps.octree_point_count_list_present_flag)
    gbh.footer.octree_lvl_num_points_minus1.reserve(maxDepth);

  for (int depth = 0; depth < maxDepth; depth++) {//八叉树划分的最大的一个for循环，变量为八叉树层数
    // setyo at the start of each level
    auto fifoCurrLvlEnd = fifo.end();//指向本层节点的最后一个节点
    int numNodesNextLvl = 0;//下一层的节点数
    Vec3<int32_t> occupancyAtlasOrigin = 0xffffffff;

    // derive per-level node size related parameters
    auto nodeSizeLog2 = lvlNodeSizeLog2[depth];//当前层的节点边长
    auto childSizeLog2 = lvlNodeSizeLog2[depth + 1];//当前层节点的子节点边长
    // represents the largest dimension of the current node
    int nodeMaxDimLog2 = nodeSizeLog2.max();//节点三个维度边长中最大的边长

    // if one dimension is not split, atlasShift[k] = 0
    int codedAxesPrevLvl = depth ? gbh.tree_lvl_coded_axis_list[depth - 1] : 7;//上一层被编码的维度，xyz组成三位，最大值为7
    int codedAxesCurLvl = gbh.tree_lvl_coded_axis_list[depth];//当前层被编码的维度，即根据qtbt指示每次划分的方式

    auto pointSortMask = qtBtChildSize(nodeSizeLog2, childSizeLog2);//各维度上子节点相对父节点边长变化量

    // Idcm quantisation applies to child nodes before per node qps
    if (--numLvlsUntilQuantization > 0) {//先平面编码，然后IDCM量化，然后节点量化
      // If planar is enabled, the planar bits are not quantised (since
      // the planar mode is determined before quantisation)
      quantNodeSizeLog2 = nodeSizeLog2;
      if (gps.geom_planar_mode_enabled_flag)
        quantNodeSizeLog2 -= 1;

      for (int k = 0; k < 3; k++)
        quantNodeSizeLog2[k] = std::max(0, quantNodeSizeLog2[k]);

      // limit the idcmQp such that it cannot overquantise the node
      auto minNs = quantNodeSizeLog2.min();//限制idcm量化
      idcmQp = gps.geom_base_qp + gps.geom_idcm_qp_offset;
      idcmQp <<= gps.geom_qp_multiplier_log2;
      idcmQp = std::min(idcmQp, minNs * 8);
    }

    // determing a per node QP at the appropriate level
    if (!numLvlsUntilQuantization) {
      // idcm qps are no longer independent
      idcmQp = 0;
      quantNodeSizeLog2 = nodeSizeLog2;
      calculateNodeQps(
        params.qpMethod,
        nodeSizeLog2, sliceQp, gps.geom_qp_multiplier_log2, fifo.begin(),
        fifoCurrLvlEnd);
    }

    // save context state for parallel coding
    if (depth == maxDepth - 1 - gbh.geom_stream_cnt_minus1)
      if (gbh.geom_stream_cnt_minus1)
        savedState.reset(new GeometryOctreeEncoder(encoder));

    // load context state for parallel coding starting one level later
    if (depth > maxDepth - 1 - gbh.geom_stream_cnt_minus1) {
      encoder = *savedState;
      encoder._arithmeticEncoder = (++arithmeticEncoderIt)->get();
    }

    auto planarDepth = gbh.rootNodeSizeLog2 - nodeSizeLog2;
    encoder.beginOctreeLevel(planarDepth);//初始化平面编码所需的类

    // process all nodes within a single level
    for (; fifo.begin() != fifoCurrLvlEnd; fifo.pop_front()) {//当前层节点编码，遍历当前层的每个节点
      PCCOctree3Node& node0 = fifo.front();//获取当前节点

      // encode delta qp for each octree block
      if (numLvlsUntilQuantization == 0) {
        int qpOffset = (node0.qp - sliceQp) >> gps.geom_qp_multiplier_log2;
        encoder.encodeQpOffset(qpOffset);
      }

      int shiftBits = QuantizerGeom::qpShift(node0.qp);
      auto effectiveNodeSizeLog2 = nodeSizeLog2 - shiftBits;//量化后的当前节点边长
      auto effectiveChildSizeLog2 = childSizeLog2 - shiftBits;//量化后的子节点边长

      // make quantisation work with qtbt and planar.
      int codedAxesCurNode = codedAxesCurLvl;//当前节点被编码的维度，即qtbt的影响
      if (shiftBits != 0) {//考虑量化的影响
        for (int k = 0; k < 3; k++) {
          if (effectiveChildSizeLog2[k] < 0)
            codedAxesCurNode &= ~(4 >> k);
        }
      }

      if (numLvlsUntilQuantization == 0) {
        geometryQuantization(pointCloud, node0, quantNodeSizeLog2);//几何量化
        if (gps.geom_unique_points_flag)
          checkDuplicatePoints(pointCloud, node0, pointIdxToDmIdx);//检查重复点
      }

      GeometryNeighPattern gnp{};
      if (gps.neighbour_avail_boundary_log2) {
        updateGeometryOccupancyAtlas(//更新节点占用信息以便于其他节点查询
          node0.pos, codedAxesPrevLvl, fifo, fifoCurrLvlEnd, &occupancyAtlas,
          &occupancyAtlasOrigin);

        gnp = makeGeometryNeighPattern(//计算NeighPattern
          gps.adjacent_child_contextualization_enabled_flag, node0.pos,
          codedAxesPrevLvl, codedAxesCurLvl, occupancyAtlas);
      } else {//如果缓存区关闭，则换用另一种方法来计算NeighPattern
        // The position of the node in the parent's occupancy map
        int posInParent = 0;
        posInParent |= (node0.pos[0] & 1) << 2;
        posInParent |= (node0.pos[1] & 1) << 1;
        posInParent |= (node0.pos[2] & 1) << 0;
        posInParent &= codedAxesPrevLvl;//计算当前节点位于父节点中的位置

        gnp.neighPattern =//根据occupancy来计算父节点六近邻
          neighPatternFromOccupancy(posInParent, node0.siblingOccupancy);
      }

      // split the current node into 8 children
      //  - perform an 8-way counting sort of the current node's points
      //  - (later) map to child nodes
      std::array<int, 8> childCounts = {};
      countingSort(//利用计数排序算法将当前节点中的所有点根据坐标分给8个子节点，完成子节点划分,
        PCCPointSet3::iterator(&pointCloud, node0.start),
        PCCPointSet3::iterator(&pointCloud, node0.end), childCounts,
        [=](const PCCPointSet3::Proxy& proxy) {
          const auto& point = *proxy;
          return !!(int(point[2]) & pointSortMask[2])//各点坐标与对应节点边长相与即得该点的一个值，然后各点依赖得到的值进行计数排序，从大到小分配给8个子节点
            | (!!(int(point[1]) & pointSortMask[1]) << 1)
            | (!!(int(point[0]) & pointSortMask[0]) << 2);
        });

      // generate the bitmap of child occupancy and count
      // the number of occupied children in node0.
      int occupancy = 0;
      int numSiblings = 0;
      for (int i = 0; i < 8; i++) {
        if (childCounts[i]) {
          occupancy |= 1 << i;//根据八个子节点中是否有点来生成8位occupancy，注意从低位开始
          numSiblings++;
        }
      }

      int contextAngle = -1;//z方向角度编码的context
      int contextAnglePhiX = -1;//x方向角度编码的context
      int contextAnglePhiY = -1;//y方向角度编码的context
      if (gps.geom_angular_mode_enabled_flag) {
        contextAngle = determineContextAngleForPlanar(//获取各个维度上角度编码的context
          node0, headPos, nodeSizeLog2, zLaser, thetaLaser, numLasers,
          deltaAngle, encoder._phiZi, encoder._phiBuffer.data(),
          &contextAnglePhiX, &contextAnglePhiY);
      }

      if (gps.geom_planar_mode_enabled_flag) {
        // update the plane rate depending on the occupancy and local density
        auto occupancy = node0.siblingOccupancy;
        auto numSiblings = node0.numSiblingsPlus1;
        if (!nodesBeforePlanarUpdate--) {//根据occupancy和节点密度来更新平面编码的概率
          encoder._planar.updateRate(occupancy, numSiblings);
          nodesBeforePlanarUpdate = numSiblings - 1;
        }
      }

      OctreeNodePlanar planar;
      if (!isLeafNode(effectiveNodeSizeLog2)) {
        // planar eligibility
        bool planarEligible[3] = {false, false, false};
        if (gps.geom_planar_mode_enabled_flag) {
          encoder._planar.isEligible(planarEligible);//两方面的限制：邻居节点数与概率更新
          if (gps.geom_angular_mode_enabled_flag) {
            if (contextAngle != -1)//如果角度模式资格判断为满足，则平面模式资格直接判为满足
              planarEligible[2] = true;
            planarEligible[0] = (contextAnglePhiX != -1);
            planarEligible[1] = (contextAnglePhiY != -1);
          }

          for (int k = 0; k < 3; k++)
            planarEligible[k] &= (codedAxesCurNode >> (2 - k)) & 1;
        }

        int planarProb[3] = {127, 127, 127};
        // determine planarity if eligible
        if (planarEligible[0] || planarEligible[1] || planarEligible[2])//决定平面模式的信息
          encoder.determinePlanarMode(
            occupancy, planarEligible, node0, planar, gnp.neighPattern,
            planarProb, contextAngle, contextAnglePhiX, contextAnglePhiY);

        node0.idcmEligible &=
          planarProb[0] * planarProb[1] * planarProb[2] <= idcmThreshold;
      }

      // At the scaling depth, it is possible for a node that has previously
      // been marked as being eligible for idcm to be fully quantised due
      // to the choice of QP.  There is therefore nothing to code with idcm.
      if (isLeafNode(effectiveNodeSizeLog2))//叶子节点不能再进行IDCM
        node0.idcmEligible = false;

      if (node0.idcmEligible) {//注意IDCM与正常occupancy的顺序，先进行IDCM的判断与编码，再去正常occupancy
        // todo(df): this is pessimistic in the presence of idcm quantisation,
        // since that is eligible may only meet the point count constraint
        // after quantisation, which is performed after the decision is taken.
        auto mode = canEncodeDirectPosition(
          gps.geom_unique_points_flag, node0, pointCloud);

        encoder.encodeIsIdcm(mode);

        if (mode != DirectMode::kUnavailable) {//当不可能IDCM时继续正常occupancy
          int idcmShiftBits = shiftBits;//似乎没有启用
          auto idcmSize = effectiveNodeSizeLog2;

          if (idcmQp) {
            node0.qp = idcmQp;
            idcmShiftBits = QuantizerGeom::qpShift(idcmQp);
            idcmSize = nodeSizeLog2 - idcmShiftBits;
            geometryQuantization(pointCloud, node0, quantNodeSizeLog2);

            if (gps.geom_unique_points_flag)
              checkDuplicatePoints(pointCloud, node0, pointIdxToDmIdx);
          }

          encoder.encodeDirectPosition(//
            mode, gps.geom_unique_points_flag, gps.joint_2pt_idcm_enabled_flag,
            idcmSize, idcmShiftBits, node0, planar, pointCloud,
            gps.geom_angular_mode_enabled_flag, headPos, zLaser, thetaLaser,
            numLasers);

          // inverse quantise any quantised positions
          geometryScale(pointCloud, node0, quantNodeSizeLog2);

          // point reordering to match decoder's order
          for (auto idx = node0.start; idx < node0.end; idx++)
            pointIdxToDmIdx[idx] = nextDmIdx++;

          // NB: by definition, this is the only child node present
          if (gps.inferred_direct_coding_mode <= 1)
            assert(node0.numSiblingsPlus1 == 1);

          // This node has no children, ensure that future nodes avoid
          // accessing stale child occupancy data.
          if (gps.adjacent_child_contextualization_enabled_flag)
            updateGeometryOccupancyAtlasOccChild(
              node0.pos, 0, &occupancyAtlas);

          continue;
        }
      }

      // when all points are quantized to a single point
      if (!isLeafNode(effectiveNodeSizeLog2)) {
        // encode child occupancy map
        assert(occupancy > 0);

        // planar mode for current node
        // mask to be used for the occupancy coding
        // (bit =1 => occupancy bit not coded due to not belonging to the plane)
        int planarMask[3] = {0, 0, 0};
        maskPlanar(planar, planarMask, codedAxesCurNode);//给出planar以及mask的信息
		//mask的作用是isplanar与planarPosition的融合，把一个1bit符号转换为对应到8位occupancy中的位置
		//注意，如果planarPos=1对应的mask为0x0f，其用法是说明对应位一定不占据，而不是对应位占据
        // generate intra prediction
        bool intraPredUsed = !(planarMask[0] | planarMask[1] | planarMask[2]);//大概率不是平面
        int occupancyIsPredicted = 0;
        int occupancyPrediction = 0;
        if (
          nodeMaxDimLog2 < gps.intra_pred_max_node_size_log2
          && gps.neighbour_avail_boundary_log2 > 0 && intraPredUsed) {
          predictGeometryOccupancyIntra(
            occupancyAtlas, node0.pos, codedAxesPrevLvl, &occupancyIsPredicted,
            &occupancyPrediction);
        }

        encoder.encodeOccupancy(
          gnp, occupancy, occupancyIsPredicted, occupancyPrediction,
          planarMask[0], planarMask[1], planarMask[2],
          planar.planarPossible & 1, planar.planarPossible & 2,
          planar.planarPossible & 4);//planarPossible指示三个方向上经过isPlanar判定的结果
      }

      // update atlas for child neighbours
      // NB: the child occupancy atlas must be updated even if the current
      //     node has no occupancy coded in order to clear any stale state in
      //     the atlas.
      if (gps.adjacent_child_contextualization_enabled_flag)
        updateGeometryOccupancyAtlasOccChild(
          node0.pos, occupancy, &occupancyAtlas);

      // Leaf nodes are immediately coded.  No further splitting occurs.
      if (isLeafNode(effectiveChildSizeLog2)) {
        int childStart = node0.start;

        // inverse quantise any quantised positions
        geometryScale(pointCloud, node0, quantNodeSizeLog2);

        for (int i = 0; i < 8; i++) {
          if (!childCounts[i]) {
            // child is empty: skip
            continue;
          }

          int childEnd = childStart + childCounts[i];
          for (auto idx = childStart; idx < childEnd; idx++)
            pointIdxToDmIdx[idx] = nextDmIdx++;

          childStart = childEnd;

          // if the bitstream is configured to represent unique points,
          // no point count is sent.
          if (gps.geom_unique_points_flag) {
            assert(childCounts[i] == 1);
            continue;
          }

          encoder.encodePositionLeafNumPoints(childCounts[i]);
        }

        // leaf nodes do not get split
        continue;
      }

      // nodeSizeLog2 > 1: for each child:
      //  - determine elegibility for IDCM
      //  - directly code point positions if IDCM allowed and selected
      //  - otherwise, insert split children into fifo while updating neighbour state
      int childPointsStartIdx = node0.start;

      for (int i = 0; i < 8; i++) {//开始创建八个子节点并将八个子节点堆入FIFO队列中,注意节点的顺序不受编码occupancy时旋转变换的影响
        if (!childCounts[i]) {//空节点即跳过
          // child is empty: skip
          continue;
        }

        // create new child
        fifo.emplace_back();
        auto& child = fifo.back();//将子节点堆入FIFO队列中

        int x = !!(i & 4);
        int y = !!(i & 2);
        int z = !!(i & 1);

        child.qp = node0.qp;
        // only shift position if an occupancy bit was coded for the axis
        child.pos[0] = (node0.pos[0] << !!(codedAxesCurLvl & 4)) + x;//计算子节点坐标位置。这种表示方式与avs不同，可以很好的得到其相对于父节点的位置信息
        child.pos[1] = (node0.pos[1] << !!(codedAxesCurLvl & 2)) + y;
        child.pos[2] = (node0.pos[2] << !!(codedAxesCurLvl & 1)) + z;

        child.start = childPointsStartIdx;//新节点中的点的起始索引
        childPointsStartIdx += childCounts[i];
        child.end = childPointsStartIdx;
        child.numSiblingsPlus1 = numSiblings;
        child.siblingOccupancy = occupancy;
        child.laserIndex = node0.laserIndex;

        child.idcmEligible = isDirectModeEligible(//在节点被刚刚创建时，就开始判定其IDCM资格
          gps.inferred_direct_coding_mode, nodeMaxDimLog2, gnp.neighPattern,
          node0, child);

        numNodesNextLvl++;//更新下一层节点数目
      }
    }

    // calculate the number of points that would be decoded if decoding were
    // to stop at this point.
    if (gps.octree_point_count_list_present_flag) {
      int numPtsAtLvl = numNodesNextLvl + nextDmIdx - 1;
      gbh.footer.octree_lvl_num_points_minus1.push_back(numPtsAtLvl);
    }
  }

  // the last element is the number of decoded points
  if (!gbh.footer.octree_lvl_num_points_minus1.empty())
    gbh.footer.octree_lvl_num_points_minus1.pop_back();

  // save the context state for re-use by a future slice if required
  ctxtMem = encoder.getCtx();

  // return partial coding result
  //  - add missing levels to node positions
  //  - inverse quantise the point cloud
  // todo(df): this does not yet support inverse quantisation of node.pos
  if (nodesRemaining) {
    auto nodeSizeLog2 = lvlNodeSizeLog2[maxDepth];
    for (auto& node : fifo) {
      for (int k = 0; k < 3; k++)
        node.pos[k] <<= nodeSizeLog2[k];
      geometryScale(pointCloud, node, quantNodeSizeLog2);
    }
    *nodesRemaining = std::move(fifo);
    return;
  }

  ////
  // The following is to re-order the points according in the decoding
  // order since IDCM causes leaves to be coded earlier than they
  // otherwise would.
  PCCPointSet3 pointCloud2;
  pointCloud2.addRemoveAttributes(
    pointCloud.hasColors(), pointCloud.hasReflectances());
  pointCloud2.resize(pointCloud.getPointCount());

  // copy points with DM points first, the rest second
  int outIdx = nextDmIdx;
  for (int i = 0; i < pointIdxToDmIdx.size(); i++) {
    int dstIdx = pointIdxToDmIdx[i];
    if (dstIdx == -1) {
      dstIdx = outIdx++;
    } else if (dstIdx == -2) {  // ignore duplicated points
      continue;
    }

    pointCloud2[dstIdx] = pointCloud[i];
    if (pointCloud.hasColors())
      pointCloud2.setColor(dstIdx, pointCloud.getColor(i));
    if (pointCloud.hasReflectances())
      pointCloud2.setReflectance(dstIdx, pointCloud.getReflectance(i));
  }
  pointCloud2.resize(outIdx);
  swap(pointCloud, pointCloud2);
}

//============================================================================

void
encodeGeometryOctree(
  const OctreeEncOpts& opt,
  const GeometryParameterSet& gps,
  GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMem,
  std::vector<std::unique_ptr<EntropyEncoder>>& arithmeticEncoders)
{
  encodeGeometryOctree(
    opt, gps, gbh, pointCloud, ctxtMem, arithmeticEncoders, nullptr);
}

//============================================================================

}  // namespace pcc
