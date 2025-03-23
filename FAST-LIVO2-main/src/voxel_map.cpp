/* 
This file is part of FAST-LIVO2: Fast, Direct LiDAR-Inertial-Visual Odometry.
FAST-LIVO2 是一个快速、直接的激光雷达-惯性-视觉里程计系统，旨在通过融合 LiDAR、IMU 和相机数据实现高效的定位与建图。
本文件定义了体素地图管理相关功能，用于构建和维护点云的体素化表示，支持 LiDAR 里程计和地图更新。

Developer: Chunran Zheng <zhengcr@connect.hku.hk>
开发者：郑春然，香港大学联系邮箱。

For commercial use, please contact me at <zhengcr@connect.hku.hk> or
Prof. Fu Zhang at <fuzhang@hku.hk>.
商业用途请联系开发者或指导教授张福。

This file is subject to the terms and conditions outlined in the 'LICENSE' file,
which is included as part of this source code package.
本文件受源代码包中包含的 'LICENSE' 文件中条款和条件的约束。
*/

#include "voxel_map.h" // 包含体素地图相关的头文件，定义类和函数声明

/*
 * 计算点在体坐标系中的协方差
 * 功能：根据点的距离和角度不确定性，计算其在体坐标系中的协方差矩阵。
 * 参数：
 *   - pb: 点在体坐标系中的位置 (Eigen::Vector3d)
 *   - range_inc: 距离不确定性 (float)
 *   - degree_inc: 角度不确定性 (float，单位为度)
 *   - cov: 输出协方差矩阵 (Eigen::Matrix3d)
 * 原理：
 *   - 基于激光雷达的测量模型，考虑距离和角度的噪声影响。
 *   - 使用方向向量和正交基向量分解不确定性，构建协方差矩阵。
 */
void calcBodyCov(Eigen::Vector3d &pb, const float range_inc, const float degree_inc, Eigen::Matrix3d &cov)
{
  if (pb[2] == 0) pb[2] = 0.0001; // 防止 z=0 导致除零，添加微小偏移
  float range = sqrt(pb[0] * pb[0] + pb[1] * pb[1] + pb[2] * pb[2]); // 计算点到原点的距离
  float range_var = range_inc * range_inc; // 距离方差，假设与距离成正比

  // 角度不确定性矩阵（2x2），表示水平和垂直方向的角度方差
  Eigen::Matrix2d direction_var;
  direction_var << pow(sin(DEG2RAD(degree_inc)), 2), 0, 
                   0, pow(sin(DEG2RAD(degree_inc)), 2);

  // 归一化方向向量
  Eigen::Vector3d direction(pb);
  direction.normalize();

  // 方向向量的反对称矩阵，用于计算旋转影响
  Eigen::Matrix3d direction_hat;
  direction_hat << 0, -direction(2), direction(1),
                   direction(2), 0, -direction(0),
                   -direction(1), direction(0), 0;

  // 构造两个正交基向量，形成局部坐标系
  Eigen::Vector3d base_vector1(1, 1, -(direction(0) + direction(1)) / direction(2));
  base_vector1.normalize();
  Eigen::Vector3d base_vector2 = base_vector1.cross(direction); // 叉乘得到第三个正交向量
  base_vector2.normalize();

  // 组合正交基向量到矩阵 N
  Eigen::Matrix<double, 3, 2> N;
  N << base_vector1(0), base_vector2(0),
       base_vector1(1), base_vector2(1),
       base_vector1(2), base_vector2(2);

  // 计算协方差：距离方差沿方向向量分布，角度方差沿正交方向分布
  Eigen::Matrix<double, 3, 2> A = range * direction_hat * N;
  cov = direction * range_var * direction.transpose() + A * direction_var * A.transpose();
}

/*
 * 加载体素地图配置
 * 功能：从 ROS 参数服务器加载体素地图的配置参数。
 * 参数：
 *   - nh: ROS 节点句柄
 *   - voxel_config: 输出配置结构体 (VoxelMapConfig)
 * 原理：通过 ROS 参数接口读取预设值，配置体素地图的层级、体素大小、阈值等。
 */
void loadVoxelConfig(ros::NodeHandle &nh, VoxelMapConfig &voxel_config)
{
  // 是否发布平面地图，默认为 false
  nh.param<bool>("publish/pub_plane_en", voxel_config.is_pub_plane_map_, false);
  
  // LIO 参数
  nh.param<int>("lio/max_layer", voxel_config.max_layer_, 1); // 最大层级，控制八叉树深度
  nh.param<double>("lio/voxel_size", voxel_config.max_voxel_size_, 0.5); // 最大体素大小（单位：米）
  nh.param<double>("lio/min_eigen_value", voxel_config.planner_threshold_, 0.01); // 平面判断的最小特征值阈值
  nh.param<double>("lio/sigma_num", voxel_config.sigma_num_, 3); // 残差筛选的标准差倍数
  nh.param<double>("lio/beam_err", voxel_config.beam_err_, 0.02); // 激光束角度误差（单位：度）
  nh.param<double>("lio/dept_err", voxel_config.dept_err_, 0.05); // 深度误差（单位：米）
  nh.param<vector<int>>("lio/layer_init_num", voxel_config.layer_init_num_, vector<int>{5,5,5,5,5}); // 各层初始点数阈值
  nh.param<int>("lio/max_points_num", voxel_config.max_points_num_, 50); // 每个体素最大点数
  nh.param<int>("lio/max_iterations", voxel_config.max_iterations_, 5); // EKF 最大迭代次数

  // 局部地图参数
  nh.param<bool>("local_map/map_sliding_en", voxel_config.map_sliding_en, false); // 是否启用地图滑动
  nh.param<int>("local_map/half_map_size", voxel_config.half_map_size, 100); // 半地图大小（体素单位）
  nh.param<double>("local_map/sliding_thresh", voxel_config.sliding_thresh, 8); // 滑动触发阈值（单位：米）
}

/*
 * 初始化平面
 * 功能：在体素中根据点云数据初始化平面参数。
 * 参数：
 *   - points: 输入点云数据 (std::vector<pointWithVar>)
 *   - plane: 输出平面参数 (VoxelPlane 指针)
 * 原理：
 *   - 计算点的协方差矩阵和中心点。
 *   - 通过特征分解判断是否为平面，若最小特征值小于阈值，则认为是平面，并计算法向量和不确定性。
 */
void VoxelOctoTree::init_plane(const std::vector<pointWithVar> &points, VoxelPlane *plane)
{
  // 初始化平面参数
  plane->plane_var_ = Eigen::Matrix<double, 6, 6>::Zero(); // 平面参数协方差（6维：3个位置，3个法向量）
  plane->covariance_ = Eigen::Matrix3d::Zero(); // 点云协方差矩阵
  plane->center_ = Eigen::Vector3d::Zero(); // 平面中心
  plane->normal_ = Eigen::Vector3d::Zero(); // 平面法向量
  plane->points_size_ = points.size(); // 点云数量
  plane->radius_ = 0; // 平面半径

  // 计算中心点和协方差矩阵
  for (auto pv : points)
  {
    plane->covariance_ += pv.point_w * pv.point_w.transpose(); // 累加外积
    plane->center_ += pv.point_w; // 累加位置
  }
  plane->center_ = plane->center_ / plane->points_size_; // 计算平均中心
  plane->covariance_ = plane->covariance_ / plane->points_size_ - plane->center_ * plane->center_.transpose(); // 协方差公式

  // 特征分解
  Eigen::EigenSolver<Eigen::Matrix3d> es(plane->covariance_);
  Eigen::Matrix3cd evecs = es.eigenvectors(); // 特征向量
  Eigen::Vector3cd evals = es.eigenvalues(); // 特征值
  Eigen::Vector3d evalsReal = evals.real(); // 取实部

  // 找到最大、最小和中间特征值的索引
  Eigen::Matrix3f::Index evalsMin, evalsMax;
  evalsReal.rowwise().sum().minCoeff(&evalsMin);
  evalsReal.rowwise().sum().maxCoeff(&evalsMax);
  int evalsMid = 3 - evalsMin - evalsMax;

  // 提取特征向量
  Eigen::Vector3d evecMin = evecs.real().col(evalsMin); // 最小特征值对应的向量（法向量）
  Eigen::Vector3d evecMid = evecs.real().col(evalsMid);
  Eigen::Vector3d evecMax = evecs.real().col(evalsMax);

  // 平面参数的雅可比矩阵（用于不确定性传播）
  Eigen::Matrix3d J_Q;
  J_Q << 1.0 / plane->points_size_, 0, 0,
         0, 1.0 / plane->points_size_, 0,
         0, 0, 1.0 / plane->points_size_;

  // 判断是否为平面：最小特征值小于阈值
  if (evalsReal(evalsMin) < planer_threshold_)
  {
    // 计算平面参数的不确定性
    for (int i = 0; i < points.size(); i++)
    {
      Eigen::Matrix<double, 6, 3> J; // 雅可比矩阵
      Eigen::Matrix3d F; // 中间矩阵
      for (int m = 0; m < 3; m++)
      {
        if (m != (int)evalsMin) // 对非最小特征值方向
        {
          Eigen::Matrix<double, 1, 3> F_m =
              (points[i].point_w - plane->center_).transpose() / ((plane->points_size_) * (evalsReal[evalsMin] - evalsReal[m])) *
              (evecs.real().col(m) * evecs.real().col(evalsMin).transpose() + evecs.real().col(evalsMin) * evecs.real().col(m).transpose());
          F.row(m) = F_m;
        }
        else
        {
          Eigen::Matrix<double, 1, 3> F_m = Eigen::Matrix<double, 1, 3>::Zero();
          F.row(m) = F_m;
        }
      }
      J.block<3, 3>(0, 0) = evecs.real() * F; // 法向量部分
      J.block<3, 3>(3, 0) = J_Q; // 中心部分
      plane->plane_var_ += J * points[i].var * J.transpose(); // 累加协方差
    }

    // 设置平面参数
    plane->normal_ << evecs.real()(0, evalsMin), evecs.real()(1, evalsMin), evecs.real()(2, evalsMin); // 法向量
    plane->y_normal_ << evecs.real()(0, evalsMid), evecs.real()(1, evalsMid), evecs.real()(2, evalsMid); // 平面次方向
    plane->x_normal_ << evecs.real()(0, evalsMax), evecs.real()(1, evalsMax), evecs.real()(2, evalsMax); // 平面主方向
    plane->min_eigen_value_ = evalsReal(evalsMin); // 最小特征值
    plane->mid_eigen_value_ = evalsReal(evalsMid); // 中间特征值
    plane->max_eigen_value_ = evalsReal(evalsMax); // 最大特征值
    plane->radius_ = sqrt(evalsReal(evalsMax)); // 平面半径（最大分布范围）
    plane->d_ = -(plane->normal_.dot(plane->center_)); // 平面方程的常数项
    plane->is_plane_ = true; // 标记为平面
    plane->is_update_ = true; // 标记已更新
    if (!plane->is_init_) // 如果未初始化
    {
      plane->id_ = voxel_plane_id; // 分配唯一 ID
      voxel_plane_id++;
      plane->is_init_ = true;
    }
  }
  else
  {
    plane->is_update_ = true; // 标记更新
    plane->is_plane_ = false; // 非平面
  }
}

/*
 * 初始化八叉树
 * 功能：根据临时点云数据初始化体素的八叉树结构。
 * 原理：
 *   - 如果点数超过阈值，尝试初始化平面。
 *   - 若为平面，则停止分裂；否则递归分裂为子体素。
 */
void VoxelOctoTree::init_octo_tree()
{
  if (temp_points_.size() > points_size_threshold_) // 点数超过阈值
  {
    init_plane(temp_points_, plane_ptr_); // 初始化平面
    if (plane_ptr_->is_plane_ == true) // 是平面
    {
      octo_state_ = 0; // 标记为叶节点（平面）
      if (temp_points_.size() > max_points_num_) // 点数过多
      {
        update_enable_ = false; // 禁用更新
        std::vector<pointWithVar>().swap(temp_points_); // 清空临时点云，释放内存
        new_points_ = 0; // 重置新点计数
      }
    }
    else
    {
      octo_state_ = 1; // 标记为非叶节点
      cut_octo_tree(); // 递归分裂
    }
    init_octo_ = true; // 标记已初始化
    new_points_ = 0; // 重置新点计数
  }
}

/*
 * 分裂八叉树
 * 功能：将当前体素分裂为八个子体素，并分配点云。
 * 原理：
 *   - 根据点的位置分配到八个子体素。
 *   - 对每个子体素递归初始化平面或继续分裂。
 */
void VoxelOctoTree::cut_octo_tree()
{
  if (layer_ >= max_layer_) // 达到最大层级
  {
    octo_state_ = 0; // 标记为叶节点
    return;
  }
  for (size_t i = 0; i < temp_points_.size(); i++) // 遍历临时点云
  {
    int xyz[3] = {0, 0, 0}; // 子体素索引
    if (temp_points_[i].point_w[0] > voxel_center_[0]) { xyz[0] = 1; }
    if (temp_points_[i].point_w[1] > voxel_center_[1]) { xyz[1] = 1; }
    if (temp_points_[i].point_w[2] > voxel_center_[2]) { xyz[2] = 1; }
    int leafnum = 4 * xyz[0] + 2 * xyz[1] + xyz[2]; // 计算子体素编号（0-7）

    if (leaves_[leafnum] == nullptr) // 子体素未创建
    {
      leaves_[leafnum] = new VoxelOctoTree(max_layer_, layer_ + 1, layer_init_num_[layer_ + 1], max_points_num_, planer_threshold_);
      leaves_[leafnum]->layer_init_num_ = layer_init_num_; // 继承层级初始化点数
      leaves_[leafnum]->voxel_center_[0] = voxel_center_[0] + (2 * xyz[0] - 1) * quater_length_; // 设置子体素中心
      leaves_[leafnum]->voxel_center_[1] = voxel_center_[1] + (2 * xyz[1] - 1) * quater_length_;
      leaves_[leafnum]->voxel_center_[2] = voxel_center_[2] + (2 * xyz[2] - 1) * quater_length_;
      leaves_[leafnum]->quater_length_ = quater_length_ / 2; // 子体素尺寸减半
    }
    leaves_[leafnum]->temp_points_.push_back(temp_points_[i]); // 分配点到子体素
    leaves_[leafnum]->new_points_++;
  }

  // 处理每个子体素
  for (uint i = 0; i < 8; i++)
  {
    if (leaves_[i] != nullptr && leaves_[i]->temp_points_.size() > leaves_[i]->points_size_threshold_)
    {
      init_plane(leaves_[i]->temp_points_, leaves_[i]->plane_ptr_); // 初始化子体素平面
      if (leaves_[i]->plane_ptr_->is_plane_)
      {
        leaves_[i]->octo_state_ = 0; // 标记为叶节点
        if (leaves_[i]->temp_points_.size() > leaves_[i]->max_points_num_)
        {
          leaves_[i]->update_enable_ = false; // 禁用更新
          std::vector<pointWithVar>().swap(leaves_[i]->temp_points_); // 清空点云
          new_points_ = 0;
        }
      }
      else
      {
        leaves_[i]->octo_state_ = 1; // 标记为非叶节点
        leaves_[i]->cut_octo_tree(); // 递归分裂
      }
      leaves_[i]->init_octo_ = true;
      leaves_[i]->new_points_ = 0;
    }
  }
}

/*
 * 更新八叉树
 * 功能：将新点插入八叉树并更新结构。
 * 参数：
 *   - pv: 新点数据 (pointWithVar)
 * 原理：
 *   - 若未初始化，添加点并检查是否需要初始化。
 *   - 若已初始化，根据平面状态递归插入到子体素或更新当前平面。
 */
void VoxelOctoTree::UpdateOctoTree(const pointWithVar &pv)
{
  if (!init_octo_) // 未初始化
  {
    new_points_++;
    temp_points_.push_back(pv);
    if (temp_points_.size() > points_size_threshold_) { init_octo_tree(); } // 点数超限，初始化
  }
  else
  {
    if (plane_ptr_->is_plane_) // 当前是平面
    {
      if (update_enable_) // 可更新
      {
        new_points_++;
        temp_points_.push_back(pv);
        if (new_points_ > update_size_threshold_) // 新点超限，更新平面
        {
          init_plane(temp_points_, plane_ptr_);
          new_points_ = 0;
        }
        if (temp_points_.size() >= max_points_num_) // 点数达到上限
        {
          update_enable_ = false;
          std::vector<pointWithVar>().swap(temp_points_);
          new_points_ = 0;
        }
      }
    }
    else // 非平面
    {
      if (layer_ < max_layer_) // 未达最大层级
      {
        int xyz[3] = {0, 0, 0};
        if (pv.point_w[0] > voxel_center_[0]) { xyz[0] = 1; }
        if (pv.point_w[1] > voxel_center_[1]) { xyz[1] = 1; }
        if (pv.point_w[2] > voxel_center_[2]) { xyz[2] = 1; }
        int leafnum = 4 * xyz[0] + 2 * xyz[1] + xyz[2];
        if (leaves_[leafnum] != nullptr) // 子体素存在
        {
          leaves_[leafnum]->UpdateOctoTree(pv); // 递归更新
        }
        else // 创建新子体素
        {
          leaves_[leafnum] = new VoxelOctoTree(max_layer_, layer_ + 1, layer_init_num_[layer_ + 1], max_points_num_, planer_threshold_);
          leaves_[leafnum]->layer_init_num_ = layer_init_num_;
          leaves_[leafnum]->voxel_center_[0] = voxel_center_[0] + (2 * xyz[0] - 1) * quater_length_;
          leaves_[leafnum]->voxel_center_[1] = voxel_center_[1] + (2 * xyz[1] - 1) * quater_length_;
          leaves_[leafnum]->voxel_center_[2] = voxel_center_[2] + (2 * xyz[2] - 1) * quater_length_;
          leaves_[leafnum]->quater_length_ = quater_length_ / 2;
          leaves_[leafnum]->UpdateOctoTree(pv);
        }
      }
      else if (update_enable_) // 已达最大层级且可更新
      {
        new_points_++;
        temp_points_.push_back(pv);
        if (new_points_ > update_size_threshold_)
        {
          init_plane(temp_points_, plane_ptr_);
          new_points_ = 0;
        }
        if (temp_points_.size() > max_points_num_)
        {
          update_enable_ = false;
          std::vector<pointWithVar>().swap(temp_points_);
          new_points_ = 0;
        }
      }
    }
  }
}

/*
 * 查找对应体素
 * 功能：根据点的位置找到对应的八叉树节点。
 * 参数：
 *   - pw: 点的位置 (Eigen::Vector3d)
 * 返回值：对应的 VoxelOctoTree 指针
 * 原理：递归遍历八叉树，找到包含该点的叶节点。
 */
VoxelOctoTree *VoxelOctoTree::find_correspond(Eigen::Vector3d pw)
{
  if (!init_octo_ || plane_ptr_->is_plane_ || (layer_ >= max_layer_)) return this; // 未初始化、是平面或达最大层级，返回当前节点

  int xyz[3] = {0, 0, 0};
  xyz[0] = pw[0] > voxel_center_[0] ? 1 : 0;
  xyz[1] = pw[1] > voxel_center_[1] ? 1 : 0;
  xyz[2] = pw[2] > voxel_center_[2] ? 1 : 0;
  int leafnum = 4 * xyz[0] + 2 * xyz[1] + xyz[2];

  return (leaves_[leafnum] != nullptr) ? leaves_[leafnum]->find_correspond(pw) : this; // 递归查找或返回当前节点
}

/*
 * 插入新点
 * 功能：将新点插入八叉树并返回插入位置的节点。
 * 参数：
 *   - pv: 新点数据 (pointWithVar)
 * 返回值：插入位置的 VoxelOctoTree 指针
 * 原理：类似 UpdateOctoTree，但不触发初始化，仅插入点。
 */
VoxelOctoTree *VoxelOctoTree::Insert(const pointWithVar &pv)
{
  if ((!init_octo_) || (init_octo_ && plane_ptr_->is_plane_) || (init_octo_ && (!plane_ptr_->is_plane_) && (layer_ >= max_layer_)))
  {
    new_points_++;
    temp_points_.push_back(pv);
    return this; // 未初始化、是平面或达最大层级，直接插入当前节点
  }

  if (init_octo_ && (!plane_ptr_->is_plane_) && (layer_ < max_layer_))
  {
    int xyz[3] = {0, 0, 0};
    xyz[0] = pv.point_w[0] > voxel_center_[0] ? 1 : 0;
    xyz[1] = pv.point_w[1] > voxel_center_[1] ? 1 : 0;
    xyz[2] = pv.point_w[2] > voxel_center_[2] ? 1 : 0;
    int leafnum = 4 * xyz[0] + 2 * xyz[1] + xyz[2];
    if (leaves_[leafnum] != nullptr) { return leaves_[leafnum]->Insert(pv); } // 子体素存在，递归插入
    else // 创建新子体素
    {
      leaves_[leafnum] = new VoxelOctoTree(max_layer_, layer_ + 1, layer_init_num_[layer_ + 1], max_points_num_, planer_threshold_);
      leaves_[leafnum]->layer_init_num_ = layer_init_num_;
      leaves_[leafnum]->voxel_center_[0] = voxel_center_[0] + (2 * xyz[0] - 1) * quater_length_;
      leaves_[leafnum]->voxel_center_[1] = voxel_center_[1] + (2 * xyz[1] - 1) * quater_length_;
      leaves_[leafnum]->voxel_center_[2] = voxel_center_[2] + (2 * xyz[2] - 1) * quater_length_;
      leaves_[leafnum]->quater_length_ = quater_length_ / 2;
      return leaves_[leafnum]->Insert(pv);
    }
  }
  return nullptr; // 未处理情况，返回空指针
}

/*
 * 状态估计
 * 功能：使用扩展卡尔曼滤波 (EKF) 更新系统状态（位置和姿态）。
 * 参数：
 *   - state_propagat: 预测状态 (StatesGroup)
 * 原理：
 *   - 将体坐标系点变换到世界坐标系。
 *   - 构建点到平面的残差，迭代优化状态。
 *   - 更新协方差矩阵，判断收敛。
 */
void VoxelMapManager::StateEstimation(StatesGroup &state_propagat)
{
  cross_mat_list_.clear(); // 清空交叉矩阵列表
  cross_mat_list_.reserve(feats_down_size_); // 预分配空间
  body_cov_list_.clear(); // 清空体坐标系协方差列表
  body_cov_list_.reserve(feats_down_size_);

  // 计算每个点的协方差和交叉矩阵
  for (size_t i = 0; i < feats_down_body_->size(); i++)
  {
    V3D point_this(feats_down_body_->points[i].x, feats_down_body_->points[i].y, feats_down_body_->points[i].z);
    if (point_this[2] == 0) { point_this[2] = 0.001; } // 防止 z=0
    M3D var;
    calcBodyCov(point_this, config_setting_.dept_err_, config_setting_.beam_err_, var); // 计算协方差
    body_cov_list_.push_back(var);
    point_this = extR_ * point_this + extT_; // 变换到外部坐标系
    M3D point_crossmat;
    point_crossmat << SKEW_SYM_MATRX(point_this); // 反对称矩阵
    cross_mat_list_.push_back(point_crossmat);
  }

  vector<pointWithVar>().swap(pv_list_); // 清空点列表
  pv_list_.resize(feats_down_size_);

  int rematch_num = 0; // 重匹配次数
  MD(DIM_STATE, DIM_STATE) G, H_T_H, I_STATE; // EKF 矩阵
  G.setZero();
  H_T_H.setZero();
  I_STATE.setIdentity();

  bool flg_EKF_inited, flg_EKF_converged, EKF_stop_flg = 0;
  for (int iterCount = 0; iterCount < config_setting_.max_iterations_; iterCount++) // 迭代优化
  {
    double total_residual = 0.0;
    pcl::PointCloud<pcl::PointXYZI>::Ptr world_lidar(new pcl::PointCloud<pcl::PointXYZI>);
    TransformLidar(state_.rot_end, state_.pos_end, feats_down_body_, world_lidar); // 变换点云到世界坐标系
    M3D rot_var = state_.cov.block<3, 3>(0, 0); // 旋转协方差
    M3D t_var = state_.cov.block<3, 3>(3, 3); // 平移协方差

    // 更新点云协方差
    for (size_t i = 0; i < feats_down_body_->size(); i++)
    {
      pointWithVar &pv = pv_list_[i];
      pv.point_b << feats_down_body_->points[i].x, feats_down_body_->points[i].y, feats_down_body_->points[i].z;
      pv.point_w << world_lidar->points[i].x, world_lidar->points[i].y, world_lidar->points[i].z;
      M3D cov = body_cov_list_[i];
      M3D point_crossmat = cross_mat_list_[i];
      cov = state_.rot_end * cov * state_.rot_end.transpose() + (-point_crossmat) * rot_var * (-point_crossmat.transpose()) + t_var;
      pv.var = cov; // 世界坐标系协方差
      pv.body_var = body_cov_list_[i]; // 体坐标系协方差
    }
    ptpl_list_.clear(); // 清空点到平面残差列表

    BuildResidualListOMP(pv_list_, ptpl_list_); // 构建残差列表

    // 计算总残差
    for (int i = 0; i < ptpl_list_.size(); i++)
    {
      total_residual += fabs(ptpl_list_[i].dis_to_plane_);
    }
    effct_feat_num_ = ptpl_list_.size(); // 有效特征点数
    cout << "[ LIO ] Raw feature num: " << feats_undistort_->size() << ", downsampled feature num:" << feats_down_size_ 
         << " effective feature num: " << effct_feat_num_ << " average residual: " << total_residual / effct_feat_num_ << endl;

    /*** 计算测量雅可比矩阵 H 和测量协方差 ***/
    MatrixXd Hsub(effct_feat_num_, 6); // 测量雅可比矩阵
    MatrixXd Hsub_T_R_inv(6, effct_feat_num_); // H^T * R^-1
    VectorXd R_inv(effct_feat_num_); // 测量噪声协方差逆
    VectorXd meas_vec(effct_feat_num_); // 测量向量
    meas_vec.setZero();

    for (int i = 0; i < effct_feat_num_; i++)
    {
      auto &ptpl = ptpl_list_[i];
      V3D point_this(ptpl.point_b_);
      point_this = extR_ * point_this + extT_;
      M3D point_crossmat;
      point_crossmat << SKEW_SYM_MATRX(point_this);

      V3D point_world = state_propagat.rot_end * point_this + state_propagat.pos_end;
      Eigen::Matrix<double, 1, 6> J_nq; // 平面参数雅可比
      J_nq.block<1, 3>(0, 0) = point_world - ptpl_list_[i].center_;
      J_nq.block<1, 3>(0, 3) = -ptpl_list_[i].normal_;

      M3D var = state_propagat.rot_end * extR_ * ptpl_list_[i].body_cov_ * (state_propagat.rot_end * extR_).transpose(); // 点协方差
      double sigma_l = J_nq * ptpl_list_[i].plane_var_ * J_nq.transpose(); // 平面不确定性

      R_inv(i) = 1.0 / (0.001 + sigma_l + ptpl_list_[i].normal_.transpose() * var * ptpl_list_[i].normal_); // 测量噪声协方差逆

      /*** 计算测量雅可比矩阵 H ***/
      V3D A(point_crossmat * state_.rot_end.transpose() * ptpl_list_[i].normal_);
      Hsub.row(i) << VEC_FROM_ARRAY(A), ptpl_list_[i].normal_[0], ptpl_list_[i].normal_[1], ptpl_list_[i].normal_[2];
      Hsub_T_R_inv.col(i) << A[0] * R_inv(i), A[1] * R_inv(i), A[2] * R_inv(i), 
                             ptpl_list_[i].normal_[0] * R_inv(i), ptpl_list_[i].normal_[1] * R_inv(i), ptpl_list_[i].normal_[2] * R_inv(i);
      meas_vec(i) = -ptpl_list_[i].dis_to_plane_; // 测量值（负距离）
    }

    EKF_stop_flg = false;
    flg_EKF_converged = false;

    /*** 迭代卡尔曼滤波更新 ***/
    MatrixXd K(DIM_STATE, effct_feat_num_); // 卡尔曼增益
    auto &&HTz = Hsub_T_R_inv * meas_vec; // H^T * R^-1 * z
    H_T_H.block<6, 6>(0, 0) = Hsub_T_R_inv * Hsub; // H^T * R^-1 * H
    MD(DIM_STATE, DIM_STATE) &&K_1 = (H_T_H.block<DIM_STATE, DIM_STATE>(0, 0) + state_.cov.block<DIM_STATE, DIM_STATE>(0, 0).inverse()).inverse();
    G.block<DIM_STATE, 6>(0, 0) = K_1.block<DIM_STATE, 6>(0, 0) * H_T_H.block<6, 6>(0, 0);
    auto vec = state_propagat - state_; // 状态偏差
    VD(DIM_STATE) solution = K_1.block<DIM_STATE, 6>(0, 0) * HTz + vec.block<DIM_STATE, 1>(0, 0) - G.block<DIM_STATE, 6>(0, 0) * vec.block<6, 1>(0, 0);

    state_ += solution; // 更新状态
    auto rot_add = solution.block<3, 1>(0, 0); // 旋转增量
    auto t_add = solution.block<3, 1>(3, 0); // 平移增量
    if ((rot_add.norm() * 57.3 < 0.01) && (t_add.norm() * 100 < 0.015)) { flg_EKF_converged = true; } // 判断收敛（旋转<0.01°，平移<0.015m）
    V3D euler_cur = state_.rot_end.eulerAngles(2, 1, 0); // 当前欧拉角

    /*** 重匹配判断 ***/
    if (flg_EKF_converged || ((rematch_num == 0) && (iterCount == (config_setting_.max_iterations_ - 2)))) { rematch_num++; }

    /*** 收敛判断和协方差更新 ***/
    if (!EKF_stop_flg && (rematch_num >= 2 || (iterCount == config_setting_.max_iterations_ - 1)))
    {
      state_.cov.block<DIM_STATE, DIM_STATE>(0, 0) = 
          (I_STATE.block<DIM_STATE, DIM_STATE>(0, 0) - G.block<DIM_STATE, DIM_STATE>(0, 0)) * state_.cov.block<DIM_STATE, DIM_STATE>(0, 0); // 更新协方差
      position_last_ = state_.pos_end; // 更新最后位置
      geoQuat_ = tf::createQuaternionMsgFromRollPitchYaw(euler_cur(0), euler_cur(1), euler_cur(2)); // 更新四元数
      EKF_stop_flg = true;
    }
    if (EKF_stop_flg) break;
  }
}

/*
 * 变换 LiDAR 点云
 * 功能：将体坐标系点云变换到世界坐标系。
 * 参数：
 *   - rot: 旋转矩阵 (Eigen::Matrix3d)
 *   - t: 平移向量 (Eigen::Vector3d)
 *   - input_cloud: 输入点云
 *   - trans_cloud: 输出变换后的点云
 */
void VoxelMapManager::TransformLidar(const Eigen::Matrix3d rot, const Eigen::Vector3d t, const PointCloudXYZI::Ptr &input_cloud,
                                     pcl::PointCloud<pcl::PointXYZI>::Ptr &trans_cloud)
{
  pcl::PointCloud<pcl::PointXYZI>().swap(*trans_cloud); // 清空输出点云
  trans_cloud->reserve(input_cloud->size()); // 预分配空间
  for (size_t i = 0; i < input_cloud->size(); i++)
  {
    pcl::PointXYZINormal p_c = input_cloud->points[i];
    Eigen::Vector3d p(p_c.x, p_c.y, p_c.z);
    p = (rot * (extR_ * p + extT_) + t); // 变换公式：先外参变换，再世界变换
    pcl::PointXYZI pi;
    pi.x = p(0);
    pi.y = p(1);
    pi.z = p(2);
    pi.intensity = p_c.intensity;
    trans_cloud->points.push_back(pi);
  }
}

/*
 * 构建体素地图
 * 功能：根据降采样后的世界坐标系点云构建体素地图。
 * 原理：
 *   - 计算每个点的协方差并分配到体素。
 *   - 初始化每个体素的八叉树结构。
 */
void VoxelMapManager::BuildVoxelMap()
{
  float voxel_size = config_setting_.max_voxel_size_;
  float planer_threshold = config_setting_.planner_threshold_;
  int max_layer = config_setting_.max_layer_;
  int max_points_num = config_setting_.max_points_num_;
  std::vector<int> layer_init_num = config_setting_.layer_init_num_;

  std::vector<pointWithVar> input_points;
  for (size_t i = 0; i < feats_down_world_->size(); i++)
  {
    pointWithVar pv;
    pv.point_w << feats_down_world_->points[i].x, feats_down_world_->points[i].y, feats_down_world_->points[i].z;
    V3D point_this(feats_down_body_->points[i].x, feats_down_body_->points[i].y, feats_down_body_->points[i].z);
    M3D var;
    calcBodyCov(point_this, config_setting_.dept_err_, config_setting_.beam_err_, var); // 计算体坐标系协方差
    M3D point_crossmat;
    point_crossmat << SKEW_SYM_MATRX(point_this);
    var = (state_.rot_end * extR_) * var * (state_.rot_end * extR_).transpose() + 
          (-point_crossmat) * state_.cov.block<3, 3>(0, 0) * (-point_crossmat).transpose() + state_.cov.block<3, 3>(3, 3); // 世界坐标系协方差
    pv.var = var;
    input_points.push_back(pv);
  }

  uint plsize = input_points.size();
  for (uint i = 0; i < plsize; i++)
  {
    const pointWithVar p_v = input_points[i];
    float loc_xyz[3];
    for (int j = 0; j < 3; j++)
    {
      loc_xyz[j] = p_v.point_w[j] / voxel_size; // 转换为体素坐标
      if (loc_xyz[j] < 0) { loc_xyz[j] -= 1.0; } // 处理负坐标
    }
    VOXEL_LOCATION position((int64_t)loc_xyz[0], (int64_t)loc_xyz[1], (int64_t)loc_xyz[2]); // 体素位置
    auto iter = voxel_map_.find(position);
    if (iter != voxel_map_.end()) // 体素已存在
    {
      voxel_map_[position]->temp_points_.push_back(p_v);
      voxel_map_[position]->new_points_++;
    }
    else // 创建新体素
    {
      VoxelOctoTree *octo_tree = new VoxelOctoTree(max_layer, 0, layer_init_num[0], max_points_num, planer_threshold);
      voxel_map_[position] = octo_tree;
      voxel_map_[position]->quater_length_ = voxel_size / 4;
      voxel_map_[position]->voxel_center_[0] = (0.5 + position.x) * voxel_size;
      voxel_map_[position]->voxel_center_[1] = (0.5 + position.y) * voxel_size;
      voxel_map_[position]->voxel_center_[2] = (0.5 + position.z) * voxel_size;
      voxel_map_[position]->temp_points_.push_back(p_v);
      voxel_map_[position]->new_points_++;
      voxel_map_[position]->layer_init_num_ = layer_init_num;
    }
  }
  for (auto iter = voxel_map_.begin(); iter != voxel_map_.end(); ++iter)
  {
    iter->second->init_octo_tree(); // 初始化所有体素的八叉树
  }
}

/*
 * 从体素位置生成 RGB 颜色
 * 功能：根据点的位置生成可视化颜色。
 * 参数：
 *   - input_point: 输入点位置 (V3D)
 * 返回值：RGB 颜色向量 (V3F)
 * 原理：根据体素坐标的和生成周期性颜色。
 */
V3F VoxelMapManager::RGBFromVoxel(const V3D &input_point)
{
  int64_t loc_xyz[3];
  for (int j = 0; j < 3; j++)
  {
    loc_xyz[j] = floor(input_point[j] / config_setting_.max_voxel_size_); // 体素坐标
  }

  VOXEL_LOCATION position((int64_t)loc_xyz[0], (int64_t)loc_xyz[1], (int64_t)loc_xyz[2]);
  int64_t ind = loc_xyz[0] + loc_xyz[1] + loc_xyz[2];
  uint k((ind + 100000) % 3); // 周期性颜色索引
  V3F RGB((k == 0) * 255.0, (k == 1) * 255.0, (k == 2) * 255.0); // R, G, B
  return RGB;
}

/*
 * 更新体素地图
 * 功能：将新点云数据更新到体素地图中。
 * 参数：
 *   - input_points: 输入点云数据
 * 原理：类似 BuildVoxelMap，但仅更新现有地图。
 */
void VoxelMapManager::UpdateVoxelMap(const std::vector<pointWithVar> &input_points)
{
  float voxel_size = config_setting_.max_voxel_size_;
  float planer_threshold = config_setting_.planner_threshold_;
  int max_layer = config_setting_.max_layer_;
  int max_points_num = config_setting_.max_points_num_;
  std::vector<int> layer_init_num = config_setting_.layer_init_num_;
  uint plsize = input_points.size();
  for (uint i = 0; i < plsize; i++)
  {
    const pointWithVar p_v = input_points[i];
    float loc_xyz[3];
    for (int j = 0; j < 3; j++)
    {
      loc_xyz[j] = p_v.point_w[j] / voxel_size;
      if (loc_xyz[j] < 0) { loc_xyz[j] -= 1.0; }
    }
    VOXEL_LOCATION position((int64_t)loc_xyz[0], (int64_t)loc_xyz[1], (int64_t)loc_xyz[2]);
    auto iter = voxel_map_.find(position);
    if (iter != voxel_map_.end()) { voxel_map_[position]->UpdateOctoTree(p_v); }
    else
    {
      VoxelOctoTree *octo_tree = new VoxelOctoTree(max_layer, 0, layer_init_num[0], max_points_num, planer_threshold);
      voxel_map_[position] = octo_tree;
      voxel_map_[position]->layer_init_num_ = layer_init_num;
      voxel_map_[position]->quater_length_ = voxel_size / 4;
      voxel_map_[position]->voxel_center_[0] = (0.5 + position.x) * voxel_size;
      voxel_map_[position]->voxel_center_[1] = (0.5 + position.y) * voxel_size;
      voxel_map_[position]->voxel_center_[2] = (0.5 + position.z) * voxel_size;
      voxel_map_[position]->UpdateOctoTree(p_v);
    }
  }
}

/*
 * 并行构建残差列表
 * 功能：使用 OpenMP 并行计算点到平面的残差。
 * 参数：
 *   - pv_list: 输入点云数据
 *   - ptpl_list: 输出点到平面残差列表
 * 原理：
 *   - 为每个点找到对应体素，计算残差。
 *   - 使用互斥锁确保线程安全。
 */
void VoxelMapManager::BuildResidualListOMP(std::vector<pointWithVar> &pv_list, std::vector<PointToPlane> &ptpl_list)
{
  int max_layer = config_setting_.max_layer_;
  double voxel_size = config_setting_.max_voxel_size_;
  double sigma_num = config_setting_.sigma_num_;
  std::mutex mylock;
  ptpl_list.clear();
  std::vector<PointToPlane> all_ptpl_list(pv_list.size()); // 所有残差
  std::vector<bool> useful_ptpl(pv_list.size()); // 是否有效标志
  std::vector<size_t> index(pv_list.size());
  for (size_t i = 0; i < index.size(); ++i)
  {
    index[i] = i;
    useful_ptpl[i] = false;
  }
  #ifdef MP_EN
    omp_set_num_threads(MP_PROC_NUM); // 设置线程数
    #pragma omp parallel for // 并行处理
  #endif
  for (int i = 0; i < index.size(); i++)
  {
    pointWithVar &pv = pv_list[i];
    float loc_xyz[3];
    for (int j = 0; j < 3; j++)
    {
      loc_xyz[j] = pv.point_w[j] / voxel_size;
      if (loc_xyz[j] < 0) { loc_xyz[j] -= 1.0; }
    }
    VOXEL_LOCATION position((int64_t)loc_xyz[0], (int64_t)loc_xyz[1], (int64_t)loc_xyz[2]);
    auto iter = voxel_map_.find(position);
    if (iter != voxel_map_.end())
    {
      VoxelOctoTree *current_octo = iter->second;
      PointToPlane single_ptpl;
      bool is_sucess = false;
      double prob = 0;
      build_single_residual(pv, current_octo, 0, is_sucess, prob, single_ptpl); // 计算残差
      if (!is_sucess) // 如果失败，尝试邻近体素
      {
        VOXEL_LOCATION near_position = position;
        if (loc_xyz[0] > (current_octo->voxel_center_[0] + current_octo->quater_length_)) { near_position.x = near_position.x + 1; }
        else if (loc_xyz[0] < (current_octo->voxel_center_[0] - current_octo->quater_length_)) { near_position.x = near_position.x - 1; }
        if (loc_xyz[1] > (current_octo->voxel_center_[1] + current_octo->quater_length_)) { near_position.y = near_position.y + 1; }
        else if (loc_xyz[1] < (current_octo->voxel_center_[1] - current_octo->quater_length_)) { near_position.y = near_position.y - 1; }
        if (loc_xyz[2] > (current_octo->voxel_center_[2] + current_octo->quater_length_)) { near_position.z = near_position.z + 1; }
        else if (loc_xyz[2] < (current_octo->voxel_center_[2] - current_octo->quater_length_)) { near_position.z = near_position.z - 1; }
        auto iter_near = voxel_map_.find(near_position);
        if (iter_near != voxel_map_.end()) { build_single_residual(pv, iter_near->second, 0, is_sucess, prob, single_ptpl); }
      }
      if (is_sucess)
      {
        mylock.lock();
        useful_ptpl[i] = true;
        all_ptpl_list[i] = single_ptpl;
        mylock.unlock();
      }
      else
      {
        mylock.lock();
        useful_ptpl[i] = false;
        mylock.unlock();
      }
    }
  }
  for (size_t i = 0; i < useful_ptpl.size(); i++)
  {
    if (useful_ptpl[i]) { ptpl_list.push_back(all_ptpl_list[i]); } // 收集有效残差
  }
}

/*
 * 计算单个残差
 * 功能：计算点到平面的残差和概率。
 * 参数：
 *   - pv: 输入点数据
 *   - current_octo: 当前体素
 *   - current_layer: 当前层级
 *   - is_sucess: 是否成功标志
 *   - prob: 输出概率
 *   - single_ptpl: 输出点到平面数据
 * 原理：
 *   - 检查点是否在平面范围内，计算距离和不确定性。
 *   - 使用高斯分布评估匹配概率。
 */
void VoxelMapManager::build_single_residual(pointWithVar &pv, const VoxelOctoTree *current_octo, const int current_layer, bool &is_sucess,
                                            double &prob, PointToPlane &single_ptpl)
{
  int max_layer = config_setting_.max_layer_;
  double sigma_num = config_setting_.sigma_num_;
  double radius_k = 3; // 半径倍数
  Eigen::Vector3d p_w = pv.point_w;

  if (current_octo->plane_ptr_->is_plane_) // 是平面
  {
    VoxelPlane &plane = *current_octo->plane_ptr_;
    float dis_to_plane = fabs(plane.normal_.dot(p_w) + plane.d_); // 点到平面距离
    float dis_to_center = (plane.center_ - p_w).squaredNorm(); // 到中心的平方距离
    float range_dis = sqrt(dis_to_center - dis_to_plane * dis_to_plane); // 平面内距离

    if (range_dis <= radius_k * plane.radius_) // 在平面范围内
    {
      Eigen::Matrix<double, 1, 6> J_nq;
      J_nq.block<1, 3>(0, 0) = p_w - plane.center_;
      J_nq.block<1, 3>(0, 3) = -plane.normal_;
      double sigma_l = J_nq * plane.plane_var_ * J_nq.transpose() + plane.normal_.transpose() * pv.var * plane.normal_; // 总不确定性
      if (dis_to_plane < sigma_num * sqrt(sigma_l)) // 在标准差范围内
      {
        is_sucess = true;
        double this_prob = 1.0 / (sqrt(sigma_l)) * exp(-0.5 * dis_to_plane * dis_to_plane / sigma_l); // 高斯概率
        if (this_prob > prob) // 更新最佳匹配
        {
          prob = this_prob;
          pv.normal = plane.normal_;
          single_ptpl.body_cov_ = pv.body_var;
          single_ptpl.point_b_ = pv.point_b;
          single_ptpl.point_w_ = pv.point_w;
          single_ptpl.plane_var_ = plane.plane_var_;
          single_ptpl.normal_ = plane.normal_;
          single_ptpl.center_ = plane.center_;
          single_ptpl.d_ = plane.d_;
          single_ptpl.layer_ = current_layer;
          single_ptpl.dis_to_plane_ = plane.normal_.dot(p_w) + plane.d_;
        }
      }
    }
  }
  else if (current_layer < max_layer) // 非平面且可分裂
  {
    for (size_t leafnum = 0; leafnum < 8; leafnum++)
    {
      if (current_octo->leaves_[leafnum] != nullptr)
      {
        build_single_residual(pv, current_octo->leaves_[leafnum], current_layer + 1, is_sucess, prob, single_ptpl); // 递归处理子体素
      }
    }
  }
}

/*
 * 发布体素地图
 * 功能：将体素地图中的平面可视化为 ROS 消息。
 * 原理：
 *   - 遍历所有体素，收集更新过的平面。
 *   - 根据协方差生成颜色和大小，发布为圆柱体标记。
 */
void VoxelMapManager::pubVoxelMap()
{
  double max_trace = 0.25; // 最大协方差迹
  double pow_num = 0.2; // 幂函数调整
  ros::Rate loop(500); // 发布频率
  float use_alpha = 0.8; // 透明度
  visualization_msgs::MarkerArray voxel_plane;
  voxel_plane.markers.reserve(1000000); // 预分配空间
  std::vector<VoxelPlane> pub_plane_list;

  for (auto iter = voxel_map_.begin(); iter != voxel_map_.end(); iter++)
  {
    GetUpdatePlane(iter->second, config_setting_.max_layer_, pub_plane_list); // 收集更新平面
  }
  for (size_t i = 0; i < pub_plane_list.size(); i++)
  {
    V3D plane_cov = pub_plane_list[i].plane_var_.block<3, 3>(0, 0).diagonal();
    double trace = plane_cov.sum();
    if (trace >= max_trace) { trace = max_trace; }
    trace = trace * (1.0 / max_trace);
    trace = pow(trace, pow_num); // 归一化和幂调整
    uint8_t r, g, b;
    mapJet(trace, 0, 1, r, g, b); // 生成颜色
    Eigen::Vector3d plane_rgb(r / 256.0, g / 256.0, b / 256.0);
    double alpha = pub_plane_list[i].is_plane_ ? use_alpha : 0; // 平面透明度
    pubSinglePlane(voxel_plane, "plane", pub_plane_list[i], alpha, plane_rgb); // 发布单个平面
  }
  voxel_map_pub_.publish(voxel_plane);
  loop.sleep();
}

/*
 * 获取更新平面
 * 功能：递归收集八叉树中更新过的平面。
 * 参数：
 *   - current_octo: 当前体素
 *   - pub_max_voxel_layer: 最大发布层级
 *   - plane_list: 输出平面列表
 */
void VoxelMapManager::GetUpdatePlane(const VoxelOctoTree *current_octo, const int pub_max_voxel_layer, std::vector<VoxelPlane> &plane_list)
{
  if (current_octo->layer_ > pub_max_voxel_layer) { return; } // 超过最大层级
  if (current_octo->plane_ptr_->is_update_) { plane_list.push_back(*current_octo->plane_ptr_); } // 添加更新平面
  if (current_octo->layer_ < current_octo->max_layer_ && !current_octo->plane_ptr_->is_plane_)
  {
    for (size_t i = 0; i < 8; i++)
    {
      if (current_octo->leaves_[i] != nullptr) { GetUpdatePlane(current_octo->leaves_[i], pub_max_voxel_layer, plane_list); } // 递归子体素
    }
  }
}

/*
 * 发布单个平面
 * 功能：将平面参数转换为 ROS 可视化标记。
 * 参数：
 *   - plane_pub: 输出标记数组
 *   - plane_ns: 命名空间
 *   - single_plane: 平面参数
 *   - alpha: 透明度
 *   - rgb: 颜色
 */
void VoxelMapManager::pubSinglePlane(visualization_msgs::MarkerArray &plane_pub, const std::string plane_ns, const VoxelPlane &single_plane,
                                     const float alpha, const Eigen::Vector3d rgb)
{
  visualization_msgs::Marker plane;
  plane.header.frame_id = "camera_init";
  plane.header.stamp = ros::Time();
  plane.ns = plane_ns;
  plane.id = single_plane.id_;
  plane.type = visualization_msgs::Marker::CYLINDER; // 使用圆柱体表示平面
  plane.action = visualization_msgs::Marker::ADD;
  plane.pose.position.x = single_plane.center_[0];
  plane.pose.position.y = single_plane.center_[1];
  plane.pose.position.z = single_plane.center_[2];
  geometry_msgs::Quaternion q;
  CalcVectQuation(single_plane.x_normal_, single_plane.y_normal_, single_plane.normal_, q); // 计算四元数
  plane.pose.orientation = q;
  plane.scale.x = 3 * sqrt(single_plane.max_eigen_value_); // x 方向尺寸
  plane.scale.y = 3 * sqrt(single_plane.mid_eigen_value_); // y 方向尺寸
  plane.scale.z = 2 * sqrt(single_plane.min_eigen_value_); // z 方向厚度
  plane.color.a = alpha;
  plane.color.r = rgb(0);
  plane.color.g = rgb(1);
  plane.color.b = rgb(2);
  plane.lifetime = ros::Duration();
  plane_pub.markers.push_back(plane);
}

/*
 * 计算向量四元数
 * 功能：根据三个正交向量计算旋转四元数。
 * 参数：
 *   - x_vec, y_vec, z_vec: 平面的三个正交方向
 *   - q: 输出四元数
 */
void VoxelMapManager::CalcVectQuation(const Eigen::Vector3d &x_vec, const Eigen::Vector3d &y_vec, const Eigen::Vector3d &z_vec,
                                      geometry_msgs::Quaternion &q)
{
  Eigen::Matrix3d rot;
  rot << x_vec(0), x_vec(1), x_vec(2),
         y_vec(0), y_vec(1), y_vec(2),
         z_vec(0), z_vec(1), z_vec(2);
  Eigen::Matrix3d rotation = rot.transpose();
  Eigen::Quaterniond eq(rotation);
  q.w = eq.w();
  q.x = eq.x();
  q.y = eq.y();
  q.z = eq.z();
}

/*
 * 颜色映射
 * 功能：将标量值映射到 RGB 颜色（Jet 配色方案）。
 * 参数：
 *   - v: 输入值
 *   - vmin, vmax: 值范围
 *   - r, g, b: 输出颜色分量
 */
void VoxelMapManager::mapJet(double v, double vmin, double vmax, uint8_t &r, uint8_t &g, uint8_t &b)
{
  r = 255;
  g = 255;
  b = 255;

  if (v < vmin) { v = vmin; }
  if (v > vmax) { v = vmax; }

  double dr, dg, db;

  if (v < 0.1242) // 蓝色
  {
    db = 0.504 + ((1. - 0.504) / 0.1242) * v;
    dg = dr = 0.;
  }
  else if (v < 0.3747) // 青色
  {
    db = 1.;
    dr = 0.;
    dg = (v - 0.1242) * (1. / (0.3747 - 0.1242));
  }
  else if (v < 0.6253) // 绿色到黄色
  {
    db = (0.6253 - v) * (1. / (0.6253 - 0.3747));
    dg = 1.;
    dr = (v - 0.3747) * (1. / (0.6253 - 0.3747));
  }
  else if (v < 0.8758) // 黄色到红色
  {
    db = 0.;
    dr = 1.;
    dg = (0.8758 - v) * (1. / (0.8758 - 0.6253));
  }
  else // 红色
  {
    db = 0.;
    dg = 0.;
    dr = 1. - (v - 0.8758) * ((1. - 0.504) / (1. - 0.8758));
  }

  r = (uint8_t)(255 * dr);
  g = (uint8_t)(255 * dg);
  b = (uint8_t)(255 * db);
}

/*
 * 地图滑动
 * 功能：根据当前位置滑动地图，移除超出范围的体素。
 * 原理：
 *   - 检查移动距离是否超过阈值。
 *   - 若超过，清理超出指定范围的体素。
 */
void VoxelMapManager::mapSliding()
{
  if ((position_last_ - last_slide_position).norm() < config_setting_.sliding_thresh) // 未达滑动阈值
  {
    std::cout << RED << "[DEBUG]: Last sliding length " << (position_last_ - last_slide_position).norm() << RESET << "\n";
    return;
  }

  last_slide_position = position_last_; // 更新滑动位置
  double t_sliding_start = omp_get_wtime();
  float loc_xyz[3];
  for (int j = 0; j < 3; j++)
  {
    loc_xyz[j] = position_last_[j] / config_setting_.max_voxel_size_;
    if (loc_xyz[j] < 0) { loc_xyz[j] -= 1.0; }
  }
  clearMemOutOfMap((int64_t)loc_xyz[0] + config_setting_.half_map_size, (int64_t)loc_xyz[0] - config_setting_.half_map_size,
                   (int64_t)loc_xyz[1] + config_setting_.half_map_size, (int64_t)loc_xyz[1] - config_setting_.half_map_size,
                   (int64_t)loc_xyz[2] + config_setting_.half_map_size, (int64_t)loc_xyz[2] - config_setting_.half_map_size); // 清理地图
  double t_sliding_end = omp_get_wtime();
  std::cout << RED << "[DEBUG]: Map sliding using " << t_sliding_end - t_sliding_start << " secs" << RESET << "\n";
}

/*
 * 清理超出地图范围的内存
 * 功能：删除超出指定范围的体素，释放内存。
 * 参数：
 *   - x_max, x_min, y_max, y_min, z_max, z_min: 地图范围边界
 */
void VoxelMapManager::clearMemOutOfMap(const int& x_max, const int& x_min, const int& y_max, const int& y_min, const int& z_max, const int& z_min)
{
  int delete_voxel_cout = 0;
  for (auto it = voxel_map_.begin(); it != voxel_map_.end(); )
  {
    const VOXEL_LOCATION& loc = it->first;
    bool should_remove = loc.x > x_max || loc.x < x_min || loc.y > y_max || loc.y < y_min || loc.z > z_max || loc.z < z_min;
    if (should_remove)
    {
      delete it->second; // 释放体素内存
      it = voxel_map_.erase(it); // 从地图中移除
      delete_voxel_cout++;
    }
    else
    {
      ++it;
    }
  }
  std::cout << RED << "[DEBUG]: Delete " << delete_voxel_cout << " root voxels" << RESET << "\n";
}
