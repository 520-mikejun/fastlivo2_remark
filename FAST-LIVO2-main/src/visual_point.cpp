/* 
This file is part of FAST-LIVO2: Fast, Direct LiDAR-Inertial-Visual Odometry.
FAST-LIVO2 是一个快速、直接的激光雷达-惯性-视觉里程计系统，旨在通过融合 LiDAR、IMU 和相机数据实现高效的定位与建图。
本文件定义了 VisualPoint 类，用于管理视觉特征点的 3D 位置和观测信息，是视觉里程计的重要组成部分。

Developer: Chunran Zheng <zhengcr@connect.hku.hk>
开发者：郑春然，香港大学联系邮箱。

For commercial use, please contact me at <zhengcr@connect.hku.hk> or
Prof. Fu Zhang at <fuzhang@hku.hk>.
商业用途请联系开发者或指导教授张福。

This file is subject to the terms and conditions outlined in the 'LICENSE' file,
which is included as part of this source code package.
本文件受源代码包中包含的 'LICENSE' 文件中条款和条件的约束。
*/

#include "visual_point.h"  // 包含 VisualPoint 类的头文件，定义类成员和方法声明
#include "feature.h"       // 包含 Feature 类的头文件，定义特征点的相关属性
#include <stdexcept>       // 标准异常库，用于抛出运行时错误
#include <vikit/math_utils.h> // VIKit 数学工具库，提供向量操作支持

/*
 * VisualPoint 类的构造函数
 * 功能：初始化一个视觉特征点对象，设置其 3D 位置和其他属性。
 * 参数：
 *   - pos: 特征点的 3D 位置 (Vector3d 类型)
 * 原理：创建一个新的视觉点，初始化其位置、表面法向量和状态标志，用于后续的观测管理和优化。
 */
VisualPoint::VisualPoint(const Vector3d &pos)
    : pos_(pos),                   // 设置特征点的 3D 位置
      previous_normal_(Vector3d::Zero()), // 初始化前一次的表面法向量为零向量
      normal_(Vector3d::Zero()),   // 初始化当前表面法向量为零向量
      is_converged_(false),        // 是否收敛标志，默认 false，表示优化尚未完成
      is_normal_initialized_(false), // 法向量是否初始化标志，默认 false
      has_ref_patch_(false)        // 是否有参考图像块标志，默认 false
{
  // 注意：obs_（观测列表）和 ref_patch（参考图像块指针）在构造函数中不显式初始化，
  // 依靠成员变量的默认构造（obs_ 为空列表，ref_patch 为 nullptr）
}

/*
 * VisualPoint 类的析构函数
 * 功能：释放 VisualPoint 对象占用的资源，特别是观测列表中的动态内存。
 * 原理：遍历 obs_ 列表，删除每个 Feature 对象的内存，然后清空列表，确保无内存泄漏。
 */
VisualPoint::~VisualPoint() 
{
  // 遍历观测列表，释放每个 Feature 对象的内存
  for (auto it = obs_.begin(), ite = obs_.end(); it != ite; ++it)
  {
    delete(*it); // 删除动态分配的 Feature 对象
  }
  obs_.clear(); // 清空观测列表
  ref_patch = nullptr; // 参考图像块指针置空（无需删除，由其他机制管理）
}

/*
 * 添加帧引用的函数
 * 功能：将一个特征点观测添加到 VisualPoint 的观测列表中。
 * 参数：
 *   - ftr: 指向 Feature 对象的指针，表示某个帧中的观测
 * 原理：将新的观测插入到列表头部，便于快速访问最近的观测。
 */
void VisualPoint::addFrameRef(Feature *ftr)
{
  obs_.push_front(ftr); // 将特征点插入到观测列表的头部
  // 原理：使用链表（假设 obs_ 是 std::list 或类似结构），头部插入时间复杂度为 O(1)
}

/*
 * 删除特征引用的函数
 * 功能：从观测列表中移除指定的特征点，并释放其内存。
 * 参数：
 *   - ftr: 指向要删除的 Feature 对象的指针
 * 原理：如果删除的是参考图像块，更新相关标志；然后在列表中查找并删除指定特征。
 */
void VisualPoint::deleteFeatureRef(Feature *ftr)
{
  // 检查是否删除的是参考图像块
  if (ref_patch == ftr)
  {
    ref_patch = nullptr;    // 清空参考图像块指针
    has_ref_patch_ = false; // 更新标志，表示当前无参考图像块
  }

  // 遍历观测列表，查找并删除指定特征
  for (auto it = obs_.begin(), ite = obs_.end(); it != ite; ++it)
  {
    if ((*it) == ftr) // 找到匹配的特征指针
    {
      delete((*it)); // 释放 Feature 对象的内存
      obs_.erase(it); // 从列表中移除该项
      return;         // 删除后立即退出，避免继续遍历
    }
  }
}

/*
 * 获取最近视角观测的函数
 * 功能：从观测列表中找到与给定帧位置视角最接近的特征点。
 * 参数：
 *   - framepos: 当前帧的相机位置 (Vector3d 类型)
 *   - ftr: 输出参数，返回找到的 Feature 指针
 *   - cur_px: 当前帧中的像素坐标 (Vector2d 类型，未使用但保留)
 * 返回值：
 *   - true: 找到有效观测
 *   - false: 未找到有效观测
 * 原理：基于视角方向的余弦相似度，找出与当前帧视角最接近的观测，视角差异小于 60°。
 */
bool VisualPoint::getCloseViewObs(const Vector3d &framepos, Feature *&ftr, const Vector2d &cur_px) const
{
  if (obs_.size() <= 0) return false; // 观测列表为空，直接返回 false

  // 计算当前帧到特征点的方向向量并归一化
  Vector3d obs_dir(framepos - pos_);
  obs_dir.normalize();

  auto min_it = obs_.begin(); // 初始化最小夹角的迭代器
  double min_cos_angle = 0;   // 初始化最小夹角余弦值

  // 遍历所有观测，寻找视角最接近的特征
  for (auto it = obs_.begin(), ite = obs_.end(); it != ite; ++it)
  {
    // 计算观测帧到特征点的方向向量并归一化
    Vector3d dir((*it)->T_f_w_.inverse().translation() - pos_);
    dir.normalize();

    // 计算当前帧方向与观测方向的夹角余弦
    double cos_angle = obs_dir.dot(dir);
    if (cos_angle > min_cos_angle) // 如果夹角更小（余弦值更大）
    {
      min_cos_angle = cos_angle; // 更新最小夹角
      min_it = it;               // 更新迭代器
    }
  }
  ftr = *min_it; // 设置输出特征为找到的最优观测

  // 检查夹角是否小于 60°（余弦值 > 0.5）
  // 原理：夹角大于 60° 的观测可能由于视角差异过大而不可靠
  if (min_cos_angle < 0.5)
  {
    return false; // 视角差异过大，返回 false
  }

  return true; // 找到有效观测，返回 true
}

/*
 * 查找最低得分特征的函数
 * 功能：从观测列表中找到得分最低的特征点。
 * 参数：
 *   - framepos: 当前帧的相机位置 (Vector3d 类型，未使用但保留)
 *   - ftr: 输出参数，返回找到的 Feature 指针
 * 原理：遍历所有观测，比较 Feature 的 score_ 属性，找出得分最低的特征。
 */
void VisualPoint::findMinScoreFeature(const Vector3d &framepos, Feature *&ftr) const
{
  auto min_it = obs_.begin(); // 初始化最低得分的迭代器
  float min_score = std::numeric_limits<float>::max(); // 初始化最小得分为浮点数最大值

  // 遍历所有观测，寻找得分最低的特征
  for (auto it = obs_.begin(), ite = obs_.end(); it != ite; ++it)
  {
    if ((*it)->score_ < min_score) // 如果当前得分更低
    {
      min_score = (*it)->score_; // 更新最小得分
      min_it = it;               // 更新迭代器
    }
  }
  ftr = *min_it; // 设置输出特征为最低得分的观测
}

/*
 * 删除非参考图像块特征的函数
 * 功能：从观测列表中删除除参考图像块外的所有特征点。
 * 原理：保留 ref_patch 对应的特征，删除其他观测，释放内存并更新列表。
 */
void VisualPoint::deleteNonRefPatchFeatures()
{
  // 遍历观测列表
  for (auto it = obs_.begin(); it != obs_.end();)
  {
    if (*it != ref_patch) // 如果不是参考图像块
    {
      delete *it;        // 释放 Feature 对象的内存
      it = obs_.erase(it); // 从列表中移除并移动迭代器
    }
    else
    {
      ++it; // 是参考图像块，保留并移动到下一项
    }
  }
}
