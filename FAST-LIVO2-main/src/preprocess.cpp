#include "preprocess.h"  // 包含预处理类的头文件，定义了类成员和方法声明

// 定义 LiDAR 返回类型常量，用于解析 Livox LiDAR 数据中的反射标签
#define RETURN0 0x00      // 单次回波
#define RETURN0AND1 0x10  // 多次回波（仅包含第一次和最后一次）

/*
 * Preprocess 类的默认构造函数
 * 功能：初始化预处理类的成员变量，设置默认参数。
 * 原理：为不同 LiDAR 类型的点云处理提供初始配置，参数主要用于特征提取和点云过滤。
 */
Preprocess::Preprocess() 
    : feature_enabled(0),  // 是否启用特征提取，默认关闭 (0)
      lidar_type(AVIA),    // 默认 LiDAR 类型为 Livox Avia
      blind(0.01),         // 盲区距离（单位：米），默认 0.01m，低于此距离的点将被过滤
      point_filter_num(1)  // 点云降采样间隔，默认每 1 个点取 1 个（无降采样）
{
  // 初始化特征提取和点云处理的参数
  inf_bound = 10;            // 无限远边界（单位：米），用于判断点是否过远
  N_SCANS = 6;               // 扫描线数，默认 6，适用于 Livox Avia
  group_size = 8;            // 平面判断的点组大小，默认 8 个点一组
  disA = 0.01;               // 平面判断距离参数 A，默认 0.01
  disA = 0.1;                // 注意：此处重复赋值，可能是笔误，应为 disB = 0.1，用于平面判断距离参数 B
  p2l_ratio = 225;           // 平面到线的距离比例阈值，默认 225
  limit_maxmid = 6.25;       // 最大-中间距离比例阈值，用于 Livox Avia 的平面判断
  limit_midmin = 6.25;       // 中间-最小距离比例阈值，用于 Livox Avia 的平面判断
  limit_maxmin = 3.24;       // 最大-最小距离比例阈值，用于其他 LiDAR 的平面判断
  jump_up_limit = 170.0;     // 跳跃上升角度阈值（单位：度），默认 170°
  jump_down_limit = 8.0;     // 跳跃下降角度阈值（单位：度），默认 8°
  cos160 = 160.0;            // 边缘跳跃判断的夹角阈值（单位：度），默认 160°
  edgea = 2;                 // 边缘跳跃判断的距离比例参数 a，默认 2
  edgeb = 0.1;               // 边缘跳跃判断的距离差参数 b，默认 0.1
  smallp_intersect = 172.5;  // 小平面交叉角度阈值（单位：度），默认 172.5°
  smallp_ratio = 1.2;        // 小平面距离比例阈值，默认 1.2
  given_offset_time = false; // 是否提供时间偏移，默认 false

  // 将角度转换为余弦值，便于后续计算
  // 原理：提前计算余弦值，避免在循环中重复计算三角函数，提高效率
  jump_up_limit = cos(jump_up_limit / 180 * M_PI);      // 将 170° 转换为余弦值 (~-0.9848)
  jump_down_limit = cos(jump_down_limit / 180 * M_PI);  // 将 8° 转换为余弦值 (~0.9903)
  cos160 = cos(cos160 / 180 * M_PI);                    // 将 160° 转换为余弦值 (~-0.9397)
  smallp_intersect = cos(smallp_intersect / 180 * M_PI); // 将 172.5° 转换为余弦值 (~-0.9914)
}

/*
 * Preprocess 类的析构函数
 * 功能：释放 Preprocess 对象的资源。
 * 原理：当前为空，因为没有动态分配的资源需要手动释放。
 */
Preprocess::~Preprocess() {}

/*
 * 设置预处理参数的函数
 * 功能：允许外部配置预处理参数，覆盖默认值。
 * 参数：
 *   - feat_en: 是否启用特征提取
 *   - lid_type: LiDAR 类型（如 AVIA, VELO16 等）
 *   - bld: 盲区距离
 *   - pfilt_num: 点云降采样间隔
 * 原理：提供灵活性，允许在运行时调整预处理行为。
 */
void Preprocess::set(bool feat_en, int lid_type, double bld, int pfilt_num)
{
  feature_enabled = feat_en;   // 设置特征提取开关
  lidar_type = lid_type;       // 设置 LiDAR 类型
  blind = bld;                 // 设置盲区距离
  point_filter_num = pfilt_num; // 设置降采样间隔
}

/*
 * 处理 Livox 自定义格式点云的函数
 * 功能：将 Livox ROS 驱动的 CustomMsg 转换为统一的点云格式，并可选提取特征。
 * 参数：
 *   - msg: Livox 点云消息指针
 *   - pcl_out: 输出点云（表面点）
 * 原理：调用 avia_handler 处理 Livox 数据，然后将结果赋值给输出点云。
 */
void Preprocess::process(const livox_ros_driver::CustomMsg::ConstPtr &msg, PointCloudXYZI::Ptr &pcl_out)
{
  avia_handler(msg);  // 处理 Livox Avia 点云数据
  *pcl_out = pl_surf; // 将处理后的表面点云赋值给输出
}

/*
 * 处理标准 ROS 点云消息的函数
 * 功能：根据 LiDAR 类型处理不同格式的点云数据。
 * 参数：
 *   - msg: ROS PointCloud2 消息指针
 *   - pcl_out: 输出点云（表面点）
 * 原理：通过 switch 分支调用对应 LiDAR 的处理函数，支持多种 LiDAR 类型。
 */
void Preprocess::process(const sensor_msgs::PointCloud2::ConstPtr &msg, PointCloudXYZI::Ptr &pcl_out)
{
  switch (lidar_type)
  {
  case OUST64:       // Ouster OS1-64
    oust64_handler(msg);
    break;
  case VELO16:       // Velodyne VLP-16
    velodyne_handler(msg);
    break;
  case L515:         // Intel RealSense L515
    l515_handler(msg);
    break;
  case XT32:         // Hesai XT32
    xt32_handler(msg);
    break;
  case PANDAR128:    // Hesai Pandar128
    Pandar128_handler(msg);
    break;
  default:
    printf("Error LiDAR Type: %d \n", lidar_type); // 未知类型，输出错误信息
    break;
  }
  *pcl_out = pl_surf; // 将处理后的表面点云赋值给输出
}

/*
 * 处理 Livox Avia 点云的函数
 * 功能：解析 Livox CustomMsg，提取点云并可选进行特征提取。
 * 参数：
 *   - msg: Livox 点云消息指针
 * 原理：
 *   - 将点云分为表面点 (pl_surf) 和角点 (pl_corn)，并存储完整点云 (pl_full)。
 *   - 如果启用特征提取，基于几何特性识别平面和边缘。
 *   - 否则仅进行基本过滤和降采样。
 */
void Preprocess::avia_handler(const livox_ros_driver::CustomMsg::ConstPtr &msg)
{
  pl_surf.clear();  // 清空表面点云容器
  pl_corn.clear();  // 清空角点云容器
  pl_full.clear();  // 清空完整点云容器
  double t1 = omp_get_wtime(); // 记录开始时间，用于性能分析
  int plsize = msg->point_num; // 获取点云数量
  printf("[ Preprocess ] Input point number: %d \n", plsize); // 输出输入点数

  // 预分配内存，提高效率
  pl_corn.reserve(plsize);  // 为角点云预留空间
  pl_surf.reserve(plsize);  // 为表面点云预留空间
  pl_full.resize(plsize);   // 为完整点云分配空间

  // 清空并预分配每个扫描线的缓冲区
  for (int i = 0; i < N_SCANS; i++)
  {
    pl_buff[i].clear();
    pl_buff[i].reserve(plsize);
  }
  uint valid_num = 0; // 有效点计数器

  if (feature_enabled) // 如果启用特征提取
  {
    // 遍历点云，提取有效点并按扫描线分组
    for (uint i = 1; i < plsize; i++)
    {
      // 检查点是否有效：扫描线编号合法且标签为多次回波 (RETURN0AND1)
      if ((msg->points[i].line < N_SCANS) && ((msg->points[i].tag & 0x30) == 0x10))
      {
        // 填充完整点云数据
        pl_full[i].x = msg->points[i].x;
        pl_full[i].y = msg->points[i].y;
        pl_full[i].z = msg->points[i].z;
        pl_full[i].intensity = msg->points[i].reflectivity; // 反射强度
        pl_full[i].curvature = msg->points[i].offset_time / float(1000000); // 偏移时间（单位：ms）

        // 判断点是否为新点（与前一点有明显差异）
        if ((abs(pl_full[i].x - pl_full[i - 1].x) > 1e-7) || 
            (abs(pl_full[i].y - pl_full[i - 1].y) > 1e-7) ||
            (abs(pl_full[i].z - pl_full[i - 1].z) > 1e-7))
        {
          pl_buff[msg->points[i].line].push_back(pl_full[i]); // 按扫描线存储
        }
      }
    }

    // 统计特征提取时间
    static int count = 0;
    static double time = 0.0;
    count++;
    double t0 = omp_get_wtime();

    // 对每个扫描线提取特征
    for (int j = 0; j < N_SCANS; j++)
    {
      if (pl_buff[j].size() <= 5) continue; // 点数太少，跳过
      pcl::PointCloud<PointType> &pl = pl_buff[j];
      plsize = pl.size();
      vector<orgtype> &types = typess[j]; // 存储点的类型信息
      types.clear();
      types.resize(plsize);
      plsize--;

      // 计算每个点的距离信息
      for (uint i = 0; i < plsize; i++)
      {
        types[i].range = pl[i].x * pl[i].x + pl[i].y * pl[i].y; // 到原点的水平距离平方
        vx = pl[i].x - pl[i + 1].x; // 与下一点的 x 差
        vy = pl[i].y - pl[i + 1].y; // 与下一点的 y 差
        vz = pl[i].z - pl[i + 1].z; // 与下一点的 z 差
        types[i].dista = vx * vx + vy * vy + vz * vz; // 两点间的距离平方
      }
      types[plsize].range = pl[plsize].x * pl[plsize].x + pl[plsize].y * pl[plsize].y;

      give_feature(pl, types); // 提取特征（平面、边缘等）
    }
    time += omp_get_wtime() - t0;
    printf("Feature extraction time: %lf \n", time / count); // 输出平均特征提取时间
  }
  else // 不启用特征提取，仅进行基本处理
  {
    for (uint i = 0; i < plsize; i++)
    {
      if (msg->points[i].line < N_SCANS) // 检查扫描线编号是否有效
      {
        valid_num++;
        // 填充点云数据
        pl_full[i].x = msg->points[i].x;
        pl_full[i].y = msg->points[i].y;
        pl_full[i].z = msg->points[i].z;
        pl_full[i].intensity = msg->points[i].reflectivity;
        pl_full[i].curvature = msg->points[i].offset_time / float(1000000); // 偏移时间（ms）

        // 处理时间偏移，确保连续性
        if (i == 0)
          pl_full[i].curvature = fabs(pl_full[i].curvature) < 1.0 ? pl_full[i].curvature : 0.0;
        else
        {
          // 如果时间跳跃过大，使用前一点时间加上固定增量 (假设 24000Hz)
          pl_full[i].curvature = fabs(pl_full[i].curvature - pl_full[i - 1].curvature) < 1.0
                                     ? pl_full[i].curvature
                                     : pl_full[i - 1].curvature + 0.004166667f; // 100ms / 24000
        }

        // 降采样并过滤盲区内的点
        if (valid_num % point_filter_num == 0)
        {
          if (pl_full[i].x * pl_full[i].x + pl_full[i].y * pl_full[i].y + pl_full[i].z * pl_full[i].z >= blind_sqr)
          {
            pl_surf.push_back(pl_full[i]); // 添加到表面点云
          }
        }
      }
    }
  }
  printf("[ Preprocess ] Output point number: %zu \n", pl_surf.points.size()); // 输出处理后的点数
}

/*
 * 处理 Intel RealSense L515 点云的函数
 * 功能：解析 L515 的 RGB 点云数据，转换为统一的格式。
 * 参数：
 *   - msg: ROS PointCloud2 消息指针
 * 原理：提取 XYZ 和 RGB 数据，过滤盲区内的点，忽略特征提取。
 */
void Preprocess::l515_handler(const sensor_msgs::PointCloud2::ConstPtr &msg)
{
  pl_surf.clear();  // 清空表面点云
  pl_corn.clear();  // 清空角点云
  pl_full.clear();  // 清空完整点云
  pcl::PointCloud<pcl::PointXYZRGB> pl_orig; // 原始点云（带 RGB）
  pcl::fromROSMsg(*msg, pl_orig); // 从 ROS 消息转换为 PCL 格式
  int plsize = pl_orig.size();    // 获取点云数量
  pl_corn.reserve(plsize);        // 预留角点云空间
  pl_surf.reserve(plsize);        // 预留表面点云空间

  double time_stamp = msg->header.stamp.toSec(); // 获取时间戳
  for (int i = 0; i < pl_orig.points.size(); i++)
  {
    if (i % point_filter_num != 0) continue; // 降采样

    // 计算点到原点的距离平方
    double range = pl_orig.points[i].x * pl_orig.points[i].x + 
                   pl_orig.points[i].y * pl_orig.points[i].y + 
                   pl_orig.points[i].z * pl_orig.points[i].z;
    if (range < blind_sqr) continue; // 过滤盲区内的点

    PointType added_pt; // 新点
    added_pt.x = pl_orig.points[i].x;
    added_pt.y = pl_orig.points[i].y;
    added_pt.z = pl_orig.points[i].z;
    added_pt.normal_x = pl_orig.points[i].r; // RGB 值存储在 normal 字段
    added_pt.normal_y = pl_orig.points[i].g;
    added_pt.normal_z = pl_orig.points[i].b;
    added_pt.curvature = 0.0; // L515 不提供时间偏移，设为 0

    pl_surf.points.push_back(added_pt); // 添加到表面点云
  }
  cout << "pl size:: " << pl_orig.points.size() << endl; // 输出原始点数
}

/*
 * 处理 Ouster OS1-64 点云的函数
 * 功能：解析 Ouster 点云数据，支持特征提取。
 * 参数：
 *   - msg: ROS PointCloud2 消息指针
 * 原理：根据 feature_enabled 决定是否提取特征，否则仅过滤和转换数据。
 */
void Preprocess::oust64_handler(const sensor_msgs::PointCloud2::ConstPtr &msg)
{
  pl_surf.clear();
  pl_corn.clear();
  pl_full.clear();
  pcl::PointCloud<ouster_ros::Point> pl_orig; // 原始点云
  pcl::fromROSMsg(*msg, pl_orig);
  int plsize = pl_orig.size();
  pl_corn.reserve(plsize);
  pl_surf.reserve(plsize);

  if (feature_enabled) // 特征提取模式
  {
    for (int i = 0; i < N_SCANS; i++)
    {
      pl_buff[i].clear();
      pl_buff[i].reserve(plsize); // 预分配扫描线缓冲区
    }

    // 遍历点云，按扫描线分组
    for (uint i = 0; i < plsize; i++)
    {
      double range = pl_orig.points[i].x * pl_orig.points[i].x + 
                     pl_orig.points[i].y * pl_orig.points[i].y + 
                     pl_orig.points[i].z * pl_orig.points[i].z;
      if (range < blind_sqr) continue;

      PointType added_pt;
      added_pt.x = pl_orig.points[i].x;
      added_pt.y = pl_orig.points[i].y;
      added_pt.z = pl_orig.points[i].z;
      added_pt.intensity = pl_orig.points[i].intensity;
      added_pt.normal_x = 0;
      added_pt.normal_y = 0;
      added_pt.normal_z = 0;
      double yaw_angle = atan2(added_pt.y, added_pt.x) * 57.3; // 偏航角（度）
      if (yaw_angle >= 180.0) yaw_angle -= 360.0; // 归一化到 [-180, 180]
      if (yaw_angle <= -180.0) yaw_angle += 360.0;

      added_pt.curvature = pl_orig.points[i].t / 1e6; // 时间偏移（ms）

      if (pl_orig.points[i].ring < N_SCANS) // 按扫描线存储
      {
        pl_buff[pl_orig.points[i].ring].push_back(added_pt);
      }
    }

    // 对每个扫描线提取特征
    for (int j = 0; j < N_SCANS; j++)
    {
      PointCloudXYZI &pl = pl_buff[j];
      int linesize = pl.size();
      vector<orgtype> &types = typess[j];
      types.clear();
      types.resize(linesize);
      linesize--;
      for (uint i = 0; i < linesize; i++)
      {
        types[i].range = sqrt(pl[i].x * pl[i].x + pl[i].y * pl[i].y);
        vx = pl[i].x - pl[i + 1].x;
        vy = pl[i].y - pl[i + 1].y;
        vz = pl[i].z - pl[i + 1].z;
        types[i].dista = vx * vx + vy * vy + vz * vz;
      }
      types[linesize].range = sqrt(pl[linesize].x * pl[linesize].x + pl[linesize].y * pl[linesize].y);
      give_feature(pl, types);
    }
  }
  else // 非特征提取模式
  {
    double time_stamp = msg->header.stamp.toSec();
    for (int i = 0; i < pl_orig.points.size(); i++)
    {
      if (i % point_filter_num != 0) continue;

      double range = pl_orig.points[i].x * pl_orig.points[i].x + 
                     pl_orig.points[i].y * pl_orig.points[i].y + 
                     pl_orig.points[i].z * pl_orig.points[i].z;
      if (range < blind_sqr) continue;

      PointType added_pt;
      added_pt.x = pl_orig.points[i].x;
      added_pt.y = pl_orig.points[i].y;
      added_pt.z = pl_orig.points[i].z;
      added_pt.intensity = pl_orig.points[i].intensity;
      added_pt.normal_x = 0;
      added_pt.normal_y = 0;
      added_pt.normal_z = 0;
      double yaw_angle = atan2(added_pt.y, added_pt.x) * 57.3;
      if (yaw_angle >= 180.0) yaw_angle -= 360.0;
      if (yaw_angle <= -180.0) yaw_angle += 360.0;

      added_pt.curvature = pl_orig.points[i].t / 1e6; // 时间偏移（ms）

      pl_surf.points.push_back(added_pt);
    }
  }
}

/*
 * 处理 Velodyne VLP-16 点云的函数
 * 功能：解析 Velodyne 点云，支持特征提取和时间偏移计算。
 * 参数：
 *   - msg: ROS PointCloud2 消息指针
 * 原理：
 *   - 如果点云提供时间偏移，直接使用；否则根据偏航角计算。
 *   - 支持特征提取，分组处理扫描线。
 */
#define MAX_LINE_NUM 64
void Preprocess::velodyne_handler(const sensor_msgs::PointCloud2::ConstPtr &msg)
{
  pl_surf.clear();
  pl_corn.clear();
  pl_full.clear();

  pcl::PointCloud<velodyne_ros::Point> pl_orig;
  pcl::fromROSMsg(*msg, pl_orig);
  int plsize = pl_orig.points.size();
  pl_surf.reserve(plsize);

  bool is_first[MAX_LINE_NUM];           // 每条扫描线是否为首点
  double yaw_fp[MAX_LINE_NUM] = {0};     // 首点的偏航角
  double omega_l = 3.61;                 // 扫描角速度（度/ms）
  float yaw_last[MAX_LINE_NUM] = {0.0};  // 末点的偏航角
  float time_last[MAX_LINE_NUM] = {0.0}; // 上一点的时间偏移

  // 检查是否提供时间偏移
  if (pl_orig.points[plsize - 1].t > 0) { given_offset_time = true; }
  else
  {
    given_offset_time = false;
    memset(is_first, true, sizeof(is_first));
    double yaw_first = atan2(pl_orig.points[0].y, pl_orig.points[0].x) * 57.29578; // 首点偏航角
    double yaw_end = yaw_first;
    int layer_first = pl_orig.points[0].ring;
    for (uint i = plsize - 1; i > 0; i--)
    {
      if (pl_orig.points[i].ring == layer_first)
      {
        yaw_end = atan2(pl_orig.points[i].y, pl_orig.points[i].x) * 57.29578; // 末点偏航角
        break;
      }
    }
  }

  if (feature_enabled)
  {
    for (int i = 0; i < N_SCANS; i++)
    {
      pl_buff[i].clear();
      pl_buff[i].reserve(plsize);
    }

    for (int i = 0; i < plsize; i++)
    {
      PointType added_pt;
      added_pt.normal_x = 0;
      added_pt.normal_y = 0;
      added_pt.normal_z = 0;
      int layer = pl_orig.points[i].ring;
      if (layer >= N_SCANS) continue;
      added_pt.x = pl_orig.points[i].x;
      added_pt.y = pl_orig.points[i].y;
      added_pt.z = pl_orig.points[i].z;
      added_pt.intensity = pl_orig.points[i].intensity;
      added_pt.curvature = pl_orig.points[i].t / 1000.0; // 时间偏移（ms）

      if (!given_offset_time)
      {
        double yaw_angle = atan2(added_pt.y, added_pt.x) * 57.2957;
        if (is_first[layer])
        {
          yaw_fp[layer] = yaw_angle;
          is_first[layer] = false;
          added_pt.curvature = 0.0;
          yaw_last[layer] = yaw_angle;
          time_last[layer] = added_pt.curvature;
          continue;
        }

        // 计算时间偏移：角度差除以角速度
        if (yaw_angle <= yaw_fp[layer]) 
          added_pt.curvature = (yaw_fp[layer] - yaw_angle) / omega_l;
        else 
          added_pt.curvature = (yaw_fp[layer] - yaw_angle + 360.0) / omega_l;

        if (added_pt.curvature < time_last[layer]) 
          added_pt.curvature += 360.0 / omega_l;

        yaw_last[layer] = yaw_angle;
        time_last[layer] = added_pt.curvature;
      }

      pl_buff[layer].points.push_back(added_pt);
    }

    for (int j = 0; j < N_SCANS; j++)
    {
      PointCloudXYZI &pl = pl_buff[j];
      int linesize = pl.size();
      if (linesize < 2) continue;
      vector<orgtype> &types = typess[j];
      types.clear();
      types.resize(linesize);
      linesize--;
      for (uint i = 0; i < linesize; i++)
      {
        types[i].range = sqrt(pl[i].x * pl[i].x + pl[i].y * pl[i].y);
        vx = pl[i].x - pl[i + 1].x;
        vy = pl[i].y - pl[i + 1].y;
        vz = pl[i].z - pl[i + 1].z;
        types[i].dista = vx * vx + vy * vy + vz * vz;
      }
      types[linesize].range = sqrt(pl[linesize].x * pl[linesize].x + pl[linesize].y * pl[linesize].y);
      give_feature(pl, types);
    }
  }
  else
  {
    for (int i = 0; i < plsize; i++)
    {
      PointType added_pt;
      added_pt.normal_x = 0;
      added_pt.normal_y = 0;
      added_pt.normal_z = 0;
      added_pt.x = pl_orig.points[i].x;
      added_pt.y = pl_orig.points[i].y;
      added_pt.z = pl_orig.points[i].z;
      added_pt.intensity = pl_orig.points[i].intensity;
      added_pt.curvature = pl_orig.points[i].t / 1000.0;

      if (!given_offset_time)
      {
        int layer = pl_orig.points[i].ring;
        double yaw_angle = atan2(added_pt.y, added_pt.x) * 57.2957;

        if (is_first[layer])
        {
          yaw_fp[layer] = yaw_angle;
          is'Brien[layer] = false;
          added_pt.curvature = 0.0;
          yaw_last[layer] = yaw_angle;
          time_last[layer] = added_pt.curvature;
          continue;
        }

        if (yaw_angle <= yaw_fp[layer]) 
          added_pt.curvature = (yaw_fp[layer] - yaw_angle) / omega_l;
        else 
          added_pt.curvature = (yaw_fp[layer] - yaw_angle + 360.0) / omega_l;

        if (added_pt.curvature < time_last[layer]) 
          added_pt.curvature += 360.0 / omega_l;

        yaw_last[layer] = yaw_angle;
        time_last[layer] = added_pt.curvature;
      }

      if (i % point_filter_num == 0)
      {
        if (added_pt.x * added_pt.x + added_pt.y * added_pt.y + added_pt.z * added_pt.z > blind_sqr)
        {
          pl_surf.points.push_back(added_pt);
        }
      }
    }
  }
}

/*
 * 处理 Hesai Pandar128 点云的函数
 * 功能：解析 Pandar128 点云，转换并排序。
 * 参数：
 *   - msg: ROS PointCloud2 消息指针
 * 原理：提取 XYZ 和时间戳，按时间排序以便后续处理。
 */
void Preprocess::Pandar128_handler(const sensor_msgs::PointCloud2::ConstPtr &msg)
{
  pl_surf.clear();
  pcl::PointCloud<Pandar128_ros::Point> pl_orig;
  pcl::fromROSMsg(*msg, pl_orig);
  int plsize = pl_orig.points.size();
  pl_surf.reserve(plsize);

  for (int i = 0; i < plsize; i++)
  {
    PointType added_pt;
    added_pt.normal_x = 0;
    added_pt.normal_y = 0;
    added_pt.normal_z = 0;
    added_pt.x = pl_orig.points[i].x;
    added_pt.y = pl_orig.points[i].y;
    added_pt.z = pl_orig.points[i].z;
    added_pt.curvature = pl_orig.points[i].timestamp * 1000.f; // 时间戳（ms）

    if (i % point_filter_num == 0)
    {
      if (added_pt.x * added_pt.x + added_pt.y * added_pt.y + added_pt.z * added_pt.z > blind_sqr)
      {
        pl_surf.points.push_back(added_pt);
      }
    }
  }

  // 定义比较函数，按时间排序
  auto comparePoints = [](const PointType& a, const PointType& b) -> bool
  {
    return a.curvature < b.curvature;
  };
  std::sort(pl_surf.points.begin(), pl_surf.points.end(), comparePoints); // 按时间升序排序
}

/*
 * 处理 Hesai XT32 点云的函数
 * 功能：解析 XT32 点云，支持特征提取和时间偏移计算。
 * 参数：
 *   - msg: ROS PointCloud2 消息指针
 * 原理：类似 Velodyne，基于偏航角计算时间偏移，或直接使用提供的时间戳。
 */
void Preprocess::xt32_handler(const sensor_msgs::PointCloud2::ConstPtr &msg)
{
  pl_surf.clear();
  pl_corn.clear();
  pl_full.clear();
  pcl::PointCloud<xt32_ros::Point> pl_orig;
  pcl::fromROSMsg(*msg, pl_orig);
  int plsize = pl_orig.points.size();
  pl_surf.reserve(plsize);

  bool is_first[MAX_LINE_NUM];
  double yaw_fp[MAX_LINE_NUM] = {0};
  double omega_l = 3.61;
  float yaw_last[MAX_LINE_NUM] = {0.0};
  float time_last[MAX_LINE_NUM] = {0.0};

  if (pl_orig.points[plsize - 1].timestamp > 0) { given_offset_time = true; }
  else
  {
    given_offset_time = false;
    memset(is_first, true, sizeof(is_first));
    double yaw_first = atan2(pl_orig.points[0].y, pl_orig.points[0].x) * 57.29578;
    double yaw_end = yaw_first;
    int layer_first = pl_orig.points[0].ring;
    for (uint i = plsize - 1; i > 0; i--)
    {
      if (pl_orig.points[i].ring == layer_first)
      {
        yaw_end = atan2(pl_orig.points[i].y, pl_orig.points[i].x) * 57.29578;
        break;
      }
    }
  }

  double time_head = pl_orig.points[0].timestamp;

  if (feature_enabled)
  {
    for (int i = 0; i < N_SCANS; i++)
    {
      pl_buff[i].clear();
      pl_buff[i].reserve(plsize);
    }

    for (int i = 0; i < plsize; i++)
    {
      PointType added_pt;
      added_pt.normal_x = 0;
      added_pt.normal_y = 0;
      added_pt.normal_z = 0;
      int layer = pl_orig.points[i].ring;
      if (layer >= N_SCANS) continue;
      added_pt.x = pl_orig.points[i].x;
      added_pt.y = pl_orig.points[i].y;
      added_pt.z = pl_orig.points[i].z;
      added_pt.intensity = pl_orig.points[i].intensity;
      added_pt.curvature = pl_orig.points[i].timestamp / 1000.0;

      if (!given_offset_time)
      {
        double yaw_angle = atan2(added_pt.y, added_pt.x) * 57.2957;
        if (is_first[layer])
        {
          yaw_fp[layer] = yaw_angle;
          is_first[layer] = false;
          added_pt.curvature = 0.0;
          yaw_last[layer] = yaw_angle;
          time_last[layer] = added_pt.curvature;
          continue;
        }

        if (yaw_angle <= yaw_fp[layer]) 
          added_pt.curvature = (yaw_fp[layer] - yaw_angle) / omega_l;
        else 
          added_pt.curvature = (yaw_fp[layer] - yaw_angle + 360.0) / omega_l;

        if (added_pt.curvature < time_last[layer]) 
          added_pt.curvature += 360.0 / omega_l;

        yaw_last[layer] = yaw_angle;
        time_last[layer] = added_pt.curvature;
      }

      pl_buff[layer].points.push_back(added_pt);
    }

    for (int j = 0; j < N_SCANS; j++)
    {
      PointCloudXYZI &pl = pl_buff[j];
      int linesize = pl.size();
      if (linesize < 2) continue;
      vector<orgtype> &types = typess[j];
      types.clear();
      types.resize(linesize);
      linesize--;
      for (uint i = 0; i < linesize; i++)
      {
        types[i].range = sqrt(pl[i].x * pl[i].x + pl[i].y * pl[i].y);
        vx = pl[i].x - pl[i + 1].x;
        vy = pl[i].y - pl[i + 1].y;
        vz = pl[i].z - pl[i + 1].z;
        types[i].dista = vx * vx + vy * vy + vz * vz;
      }
      types[linesize].range = sqrt(pl[linesize].x * pl[linesize].x + pl[linesize].y * pl[linesize].y);
      give_feature(pl, types);
    }
  }
  else
  {
    for (int i = 0; i < plsize; i++)
    {
      PointType added_pt;
      added_pt.normal_x = 0;
      added_pt.normal_y = 0;
      added_pt.normal_z = 0;
      added_pt.x = pl_orig.points[i].x;
      added_pt.y = pl_orig.points[i].y;
      added_pt.z = pl_orig.points[i].z;
      added_pt.intensity = pl_orig.points[i].intensity;
      added_pt.curvature = (pl_orig.points[i].timestamp - time_head) * 1000.f;

      if (i % point_filter_num == 0)
      {
        if (added_pt.x * added_pt.x + added_pt.y * added_pt.y + added_pt.z * added_pt.z > blind_sqr)
        {
          pl_surf.points.push_back(added_pt);
        }
      }
    }
  }
}



/*
 * 特征提取函数
 * 功能：从点云中提取几何特征（如平面、边缘），为后续 SLAM 或建图提供关键点。
 * 参数：
 *   - pl: 输入点云 (PointCloudXYZI 类型)，包含点的 XYZ 坐标和时间偏移 (curvature)
 *   - types: 点类型容器 (vector<orgtype>)，存储每个点的特征信息（如距离、类型等）
 * 原理：
 *   - 分为三个阶段：平面提取、边缘提取和平面点降采样。
 *   - 平面提取：基于点组的几何连续性判断平面。
 *   - 边缘提取：检测相邻点间的距离和角度跳跃，识别边缘特征。
 *   - 输出：将点分类为平面 (Real_Plane)、边缘 (Edge_Jump) 等，并存储到 pl_surf 和 pl_corn。
 */
void Preprocess::give_feature(pcl::PointCloud<PointType> &pl, vector<orgtype> &types)
{
  int plsize = pl.size(); // 获取点云大小
  int plsize2; // 用于限制循环范围的变量
  if (plsize == 0)
  {
    printf("something wrong\n"); // 点云为空，输出错误提示并退出
    return;
  }
  uint head = 0; // 起始点索引，跳过盲区内的点

  // 跳过盲区内的点（距离小于 blind_sqr）
  // 原理：盲区内的点通常是噪声或无效点，需过滤以提高特征提取的可靠性
  while (types[head].range < blind_sqr)
  {
    head++;
  }

  // 第一阶段：表面（平面）特征提取
  plsize2 = (plsize > group_size) ? (plsize - group_size) : 0; // 避免越界，留出 group_size 个点

  Eigen::Vector3d curr_direct(Eigen::Vector3d::Zero()); // 当前点组的方向向量
  Eigen::Vector3d last_direct(Eigen::Vector3d::Zero()); // 上一个点组的方向向量

  uint i_nex = 0, i2; // i_nex: 下一组点的结束索引；i2: 当前处理的点索引副本
  uint last_i = 0;    // 上次处理的起始点索引
  uint last_i_nex = 0;// 上次处理的结束点索引
  int last_state = 0; // 上次状态（1 表示平面，0 表示非平面）
  int plane_type;     // 平面判断结果

  // 遍历点云，提取平面特征
  for (uint i = head; i < plsize2; i++)
  {
    if (types[i].range < blind_sqr) { continue; } // 跳过盲区内的点

    i2 = i; // 记录当前点索引

    // 判断当前点是否属于平面，返回类型并更新 i_nex 和 curr_direct
    plane_type = plane_judge(pl, types, i, i_nex, curr_direct);

    if (plane_type == 1) // 是平面
    {
      // 标记点组内的点为平面类型
      // 原理：边界点标记为 Poss_Plane（可能平面），中间点标记为 Real_Plane（真实平面）
      for (uint j = i; j <= i_nex; j++)
      {
        if (j != i && j != i_nex) { types[j].ftype = Real_Plane; } // 中间点
        else { types[j].ftype = Poss_Plane; } // 边界点
      }

      // 检查与上一个平面的方向一致性
      // 原理：如果方向差异过大（夹角接近 90°），可能是边缘平面
      if (last_state == 1 && last_direct.norm() > 0.1)
      {
        double mod = last_direct.transpose() * curr_direct; // 计算方向向量夹角余弦
        if (mod > -0.707 && mod < 0.707) { types[i].ftype = Edge_Plane; } // 夹角在 45°~135°，标记为边缘平面
        else { types[i].ftype = Real_Plane; } // 否则为真实平面
      }

      i = i_nex - 1; // 跳到点组末尾，准备处理下一组
      last_state = 1; // 更新状态为平面
    }
    else // 非平面（plane_type == 2 或其他）
    {
      i = i_nex; // 跳到下一组起始点
      last_state = 0; // 更新状态为非平面
    }

    // 更新上一次处理的索引和方向
    last_i = i2;
    last_i_nex = i_nex;
    last_direct = curr_direct;
  }

  // 第二阶段：边缘特征提取
  plsize2 = plsize > 3 ? plsize - 3 : 0; // 留出前后点进行比较
  for (uint i = head + 3; i < plsize2; i++)
  {
    // 跳过盲区内的点或已标记为平面的点
    if (types[i].range < blind_sqr || types[i].ftype >= Real_Plane) { continue; }

    // 跳过距离过小的点（避免除零或噪声干扰）
    if (types[i - 1].dista < 1e-16 || types[i].dista < 1e-16) { continue; }

    Eigen::Vector3d vec_a(pl[i].x, pl[i].y, pl[i].z); // 当前点的向量
    Eigen::Vector3d vecs[2]; // 前后点的相对向量

    // 计算与前后点的夹角和距离特性
    for (int j = 0; j < 2; j++)
    {
      int m = -1; // 默认检查前点
      if (j == 1) { m = 1; } // j=1 时检查后点

      if (types[i + m].range < blind_sqr) // 前后点在盲区内
      {
        if (types[i].range > inf_bound) { types[i].edj[j] = Nr_inf; } // 当前点过远
        else { types[i].edj[j] = Nr_blind; } // 当前点在盲区附近
        continue;
      }

      vecs[j] = Eigen::Vector3d(pl[i + m].x, pl[i + m].y, pl[i + m].z); // 前后点的向量
      vecs[j] = vecs[j] - vec_a; // 相对向量

      // 计算夹角余弦值
      // 原理：余弦值用于判断方向变化，接近 1 表示同向，接近 -1 表示反向
      types[i].angle[j] = vec_a.dot(vecs[j]) / vec_a.norm() / vecs[j].norm();
      if (types[i].angle[j] < jump_up_limit) { types[i].edj[j] = Nr_180; } // 接近 180°
      else if (types[i].angle[j] > jump_down_limit) { types[i].edj[j] = Nr_zero; } // 接近 0°
    }

    // 计算前后向量的夹角余弦
    types[i].intersect = vecs[Prev].dot(vecs[Next]) / vecs[Prev].norm() / vecs[Next].norm();

    // 边缘跳跃判断条件
    // 原理：结合夹角、距离跳跃和几何特性识别边缘点
    if (types[i].edj[Prev] == Nr_nor && types[i].edj[Next] == Nr_zero && 
        types[i].dista > 0.0225 && types[i].dista > 4 * types[i - 1].dista)
    {
      if (types[i].intersect > cos160) // 前后向量夹角小于 160°
      {
        if (edge_jump_judge(pl, types, i, Prev)) { types[i].ftype = Edge_Jump; } // 确认边缘跳跃
      }
    }
    else if (types[i].edj[Prev] == Nr_zero && types[i].edj[Next] == Nr_nor && 
             types[i - 1].dista > 0.0225 && types[i - 1].dista > 4 * types[i].dista)
    {
      if (types[i].intersect > cos160)
      {
        if (edge_jump_judge(pl, types, i, Next)) { types[i].ftype = Edge_Jump; }
      }
    }
    else if (types[i].edj[Prev] == Nr_nor && types[i].edj[Next] == Nr_inf)
    {
      if (edge_jump_judge(pl, types, i, Prev)) { types[i].ftype = Edge_Jump; }
    }
    else if (types[i].edj[Prev] == Nr_inf && types[i].edj[Next] == Nr_nor)
    {
      if (edge_jump_judge(pl, types, i, Next)) { types[i].ftype = Edge_Jump; }
    }
    else if (types[i].edj[Prev] > Nr_nor && types[i].edj[Next] > Nr_nor)
    {
      if (types[i].ftype == Nor) { types[i].ftype = Wire; } // 细线特征（如电线）
    }
  }

  // 第三阶段：小平面调整和点云输出
  plsize2 = plsize - 1;
  double ratio; // 相邻点距离比例
  for (uint i = head + 1; i < plsize2; i++)
  {
    // 跳过盲区内的点或距离过小的点
    if (types[i].range < blind_sqr || types[i - 1].range < blind_sqr || types[i + 1].range < blind_sqr) { continue; }
    if (types[i - 1].dista < 1e-8 || types[i].dista < 1e-8) { continue; }

    if (types[i].ftype == Nor) // 未分类的普通点
    {
      // 计算前后距离比例
      if (types[i - 1].dista > types[i].dista) { ratio = types[i - 1].dista / types[i].dista; }
      else { ratio = types[i].dista / types[i - 1].dista; }

      // 判断小平面：夹角接近 180° 且距离比例较小
      if (types[i].intersect < smallp_intersect && ratio < smallp_ratio)
      {
        if (types[i - 1].ftype == Nor) { types[i - 1].ftype = Real_Plane; }
        if (types[i + 1].ftype == Nor) { types[i + 1].ftype = Real_Plane; }
        types[i].ftype = Real_Plane;
      }
    }
  }

  // 将特征点存储到表面点云和角点云
  int last_surface = -1; // 上一个平面点的起始索引
  for (uint j = head; j < plsize; j++)
  {
    if (types[j].ftype == Poss_Plane || types[j].ftype == Real_Plane) // 平面点
    {
      if (last_surface == -1) { last_surface = j; } // 记录平面起始点

      // 每隔 point_filter_num 个点取一个
      if (j == uint(last_surface + point_filter_num - 1))
      {
        PointType ap;
        ap.x = pl[j].x;
        ap.y = pl[j].y;
        ap.z = pl[j].z;
        ap.curvature = pl[j].curvature;
        pl_surf.push_back(ap); // 添加到表面点云
        last_surface = -1; // 重置
      }
    }
    else // 非平面点
    {
      if (types[j].ftype == Edge_Jump || types[j].ftype == Edge_Plane) 
      { 
        pl_corn.push_back(pl[j]); // 边缘点添加到角点云
      }
      if (last_surface != -1) // 处理未完成降采样的平面点
      {
        PointType ap;
        for (uint k = last_surface; k < j; k++) // 计算平均值
        {
          ap.x += pl[k].x;
          ap.y += pl[k].y;
          ap.z += pl[k].z;
          ap.curvature += pl[k].curvature;
        }
        ap.x /= (j - last_surface);
        ap.y /= (j - last_surface);
        ap.z /= (j - last_surface);
        ap.curvature /= (j - last_surface);
        pl_surf.push_back(ap); // 添加平均点到表面点云
        last_surface = -1;
      }
    }
  }
}

/*
 * 发布点云的函数
 * 功能：将处理后的点云转换为 ROS 消息并发布。
 * 参数：
 *   - pl: 输入点云
 *   - ct: 当前时间戳
 * 原理：将点云转换为 sensor_msgs::PointCloud2 格式，设置坐标系和时间戳。
 */
void Preprocess::pub_func(PointCloudXYZI &pl, const ros::Time &ct)
{
  pl.height = 1; // 设置点云高度为 1（无结构化）
  pl.width = pl.size(); // 设置点云宽度为点数
  sensor_msgs::PointCloud2 output;
  pcl::toROSMsg(pl, output); // 转换为 ROS 消息
  output.header.frame_id = "livox"; // 设置坐标系为 "livox"
  output.header.stamp = ct; // 设置时间戳
}

/*
 * 平面判断函数
 * 功能：判断一组点是否构成平面，返回类型并计算方向向量。
 * 参数：
 *   - pl: 输入点云
 *   - types: 点类型容器
 *   - i_cur: 当前点索引
 *   - i_nex: 输出下一组点的结束索引
 *   - curr_direct: 输出当前方向向量
 * 返回值：
 *   - 0: 非平面
 *   - 1: 平面
 *   - 2: 包含盲区点
 * 原理：
 *   - 使用距离和几何特性判断点组是否平坦。
 *   - 计算点间距离分布和方向一致性。
 */
int Preprocess::plane_judge(const PointCloudXYZI &pl, vector<orgtype> &types, uint i_cur, uint &i_nex, Eigen::Vector3d &curr_direct)
{
  // 计算点组的最大距离阈值
  // 原理：距离阈值随当前点到原点的距离线性增加 (disA * range + disB)
  double group_dis = disA * types[i_cur].range + disB;
  group_dis = group_dis * group_dis; // 平方形式便于比较

  double two_dis; // 两点间的距离平方
  vector<double> disarr; // 存储点间距离
  disarr.reserve(20); // 预分配空间

  // 检查初始点组（group_size 个点）
  for (i_nex = i_cur; i_nex < i_cur + group_size; i_nex++)
  {
    if (types[i_nex].range < blind_sqr)
    {
      curr_direct.setZero(); // 包含盲区点，返回无效
      return 2;
    }
    disarr.push_back(types[i_nex].dista); // 存储点间距离
  }

  // 扩展点组直到超出距离阈值
  for (;;)
  {
    if ((i_cur >= pl.size()) || (i_nex >= pl.size())) break; // 防止越界

    if (types[i_nex].range < blind_sqr)
    {
      curr_direct.setZero();
      return 2;
    }
    vx = pl[i_nex].x - pl[i_cur].x; // 计算与当前点的距离
    vy = pl[i_nex].y - pl[i_cur].y;
    vz = pl[i_nex].z - pl[i_cur].z;
    two_dis = vx * vx + vy * vy + vz * vz;
    if (two_dis >= group_dis) { break; } // 超出阈值，结束扩展
    disarr.push_back(types[i_nex].dista);
    i_nex++;
  }

  // 计算点组的宽度（垂直于方向的最大距离）
  double leng_wid = 0;
  double v1[3], v2[3]; // v1: 点到起始点的向量；v2: 垂直向量
  for (uint j = i_cur + 1; j < i_nex; j++)
  {
    if ((j >= pl.size()) || (i_cur >= pl.size())) break;
    v1[0] = pl[j].x - pl[i_cur].x;
    v1[1] = pl[j].y - pl[i_cur].y;
    v1[2] = pl[j].z - pl[i_cur].z;

    // 计算垂直向量（叉积）
    v2[0] = v1[1] * vz - vy * v1[2];
    v2[1] = v1[2] * vx - v1[0] * vz;
    v2[2] = v1[0] * vy - vx * v1[1];

    double lw = v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2];
    if (lw > leng_wid) { leng_wid = lw; } // 更新最大宽度
  }

  // 判断是否为平面
  // 原理：如果长度平方与宽度的比值小于阈值，表示点分布过于细长，不是平面
  if ((two_dis * two_dis / leng_wid) < p2l_ratio)
  {
    curr_direct.setZero();
    return 0;
  }

  // 对距离数组排序（冒泡排序）
  uint disarrsize = disarr.size();
  for (uint j = 0; j < disarrsize - 1; j++)
  {
    for (uint k = j + 1; k < disarrsize; k++)
    {
      if (disarr[j] < disarr[k])
      {
        leng_wid = disarr[j];
        disarr[j] = disarr[k];
        disarr[k] = leng_wid;
      }
    }
  }

  // 检查距离分布是否过于均匀
  if (disarr[disarr.size() - 2] < 1e-16)
  {
    curr_direct.setZero();
    return 0;
  }

  // 根据 LiDAR 类型进一步验证
  if (lidar_type == AVIA)
  {
    double dismax_mid = disarr[0] / disarr[disarrsize / 2]; // 最大与中间距离比例
    double dismid_min = disarr[disarrsize / 2] / disarr[disarrsize - 2]; // 中间与最小距离比例

    if (dismax_mid >= limit_maxmid || dismid_min >= limit_midmin)
    {
      curr_direct.setZero(); // 距离分布不均匀，不是平面
      return 0;
    }
  }
  else
  {
    double dismax_min = disarr[0] / disarr[disarrsize - 2]; // 最大与最小距离比例
    if (dismax_min >= limit_maxmin)
    {
      curr_direct.setZero();
      return 0;
    }
  }

  // 是平面，设置方向向量并归一化
  curr_direct << vx, vy, vz;
  curr_direct.normalize();
  return 1;
}

/*
 * 边缘跳跃判断函数
 * 功能：确认某点是否为边缘跳跃点。
 * 参数：
 *   - pl: 输入点云
 *   - types: 点类型容器
 *   - i: 当前点索引
 *   - nor_dir: 检查方向（0: 前，1: 后）
 * 返回值：
 *   - true: 是边缘跳跃点
 *   - false: 不是边缘跳跃点
 * 原理：基于前后点的距离变化判断跳跃特性。
 */
bool Preprocess::edge_jump_judge(const PointCloudXYZI &pl, vector<orgtype> &types, uint i, Surround nor_dir)
{
  // 检查前后点是否有效（不在盲区内）
  if (nor_dir == 0) // 前向
  {
    if (types[i - 1].range < blind_sqr || types[i - 2].range < blind_sqr) { return false; }
  }
  else if (nor_dir == 1) // 后向
  {
    if (types[i + 1].range < blind_sqr || types[i + 2].range < blind_sqr) { return false; }
  }

  // 获取前后点的距离
  double d1 = types[i + nor_dir - 1].dista; // 近点距离
  double d2 = types[i + 3 * nor_dir - 2].dista; // 远点距离
  double d;

  // 确保 d1 是较大值
  if (d1 < d2)
  {
    d = d1;
    d1 = d2;
    d2 = d;
  }

  // 开方计算实际距离
  d1 = sqrt(d1);
  d2 = sqrt(d2);

  // 判断跳跃条件
  // 原理：如果距离比例过大 (d1 > edgea * d2) 或距离差过大 (d1 - d2 > edgeb)，不是边缘
  if (d1 > edgea * d2 || (d1 - d2) > edgeb) { return false; }

  return true; // 满足条件，是边缘跳跃点
}
