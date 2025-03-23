/* 
LIVMapper 类的主要任务是实现多传感器融合的实时定位与建图（SLAM，Simultaneous Localization and Mapping），结合激光雷达（LiDAR）、惯性测量单元（IMU）和相机（视觉）数据，构建环境地图并估计机器人在其中的位姿。
其核心功能包括：
传感器数据处理与同步：接收并同步 LiDAR、IMU 和图像数据，处理时间戳和数据格式。
状态估计：通过激光-惯性（LIO）或视觉-惯性（VIO）里程计，实时估计机器人的位置、姿态、速度等状态。
地图构建：基于体素地图（Voxel Map）管理点云数据，构建并更新环境的三维地图。
数据发布与保存：将估计的位姿、点云地图和路径发布到 ROS 话题，并支持保存点云数据（PCD 格式）。
*/

#include "LIVMapper.h"

// LIVMapper 类的构造函数，初始化成员变量
LIVMapper::LIVMapper(ros::NodeHandle &nh)
    : extT(0, 0, 0),           // 初始化 LiDAR 到 IMU 的平移向量为零向量
      extR(M3D::Identity())    // 初始化 LiDAR 到 IMU 的旋转矩阵为单位矩阵
{
    // 初始化外参向量，长度为 3 或 9，初始值为 0
    extrinT.assign(3, 0.0);         // LiDAR 到 IMU 的平移向量
    extrinR.assign(9, 0.0);         // LiDAR 到 IMU 的旋转矩阵（展平为向量）
    cameraextrinT.assign(3, 0.0);   // 相机到 LiDAR 的平移向量
    cameraextrinR.assign(9, 0.0);   // 相机到 LiDAR 的旋转矩阵（展平为向量）

    // 初始化预处理和 IMU 处理的智能指针
    p_pre.reset(new Preprocess());  // 点云预处理对象
    p_imu.reset(new ImuProcess());  // IMU 数据处理对象

    // 从 ROS 参数服务器读取参数
    readParameters(nh);

    // 初始化体素地图配置
    VoxelMapConfig voxel_config;    // 定义体素地图配置对象
    loadVoxelConfig(nh, voxel_config); // 从 ROS 加载体素地图配置

    // 初始化各种点云数据的智能指针
    visual_sub_map.reset(new PointCloudXYZI());          // 视觉子地图点云
    feats_undistort.reset(new PointCloudXYZI());         // 去畸变后的特征点云
    feats_down_body.reset(new PointCloudXYZI());         // 降采样后的机体坐标系点云
    feats_down_world.reset(new PointCloudXYZI());        // 降采样后的世界坐标系点云
    pcl_w_wait_pub.reset(new PointCloudXYZI());          // 等待发布的点云（世界坐标系）
    pcl_wait_pub.reset(new PointCloudXYZI());            // 等待发布的点云（通用）
    pcl_wait_save.reset(new PointCloudXYZRGB());         // 用于保存的 RGB 点云
    pcl_wait_save_intensity.reset(new PointCloudXYZI()); // 用于保存的带强度点云

    // 初始化体素地图管理和 VIO 管理器
    voxelmap_manager.reset(new VoxelMapManager(voxel_config, voxel_map)); // 体素地图管理器
    vio_manager.reset(new VIOManager());                                 // 视觉-惯性里程计管理器

    // 设置根目录并初始化文件和组件
    root_dir = ROOT_DIR;             // 项目根目录
    initializeFiles();               // 初始化文件输出流
    initializeComponents();          // 初始化系统组件

    // 初始化路径消息的头部
    path.header.stamp = ros::Time::now();  // 设置当前时间戳
    path.header.frame_id = "camera_init";  // 设置参考坐标系为 "camera_init"
}

// LIVMapper 类的析构函数，当前为空
LIVMapper::~LIVMapper() {}

// 从 ROS 参数服务器读取参数
void LIVMapper::readParameters(ros::NodeHandle &nh)
{
    // 读取通用参数
    nh.param<string>("common/lid_topic", lid_topic, "/livox/lidar");       // LiDAR 数据话题
    nh.param<string>("common/imu_topic", imu_topic, "/livox/imu");         // IMU 数据话题
    nh.param<bool>("common/ros_driver_bug_fix", ros_driver_fix_en, false); // 是否启用 ROS 驱动修复
    nh.param<int>("common/img_en", img_en, 1);                             // 是否启用图像处理（1 为启用）
    nh.param<int>("common/lidar_en", lidar_en, 1);                         // 是否启用 LiDAR 处理（1 为启用）
    nh.param<string>("common/img_topic", img_topic, "/left_camera/image"); // 图像数据话题

    // 读取视觉-惯性里程计（VIO）参数
    nh.param<bool>("vio/normal_en", normal_en, true);                     // 是否启用法向量计算
    nh.param<bool>("vio/inverse_composition_en", inverse_composition_en, false); // 是否启用逆合成方法
    nh.param<int>("vio/max_iterations", max_iterations, 5);               // 最大迭代次数
    nh.param<double>("vio/img_point_cov", IMG_POINT_COV, 100);            // 图像点协方差
    nh.param<bool>("vio/raycast_en", raycast_en, false);                  // 是否启用光线投射
    nh.param<bool>("vio/exposure_estimate_en", exposure_estimate_en, true); // 是否启用曝光估计
    nh.param<double>("vio/inv_expo_cov", inv_expo_cov, 0.2);              // 逆曝光协方差
    nh.param<int>("vio/grid_size", grid_size, 5);                         // 网格大小
    nh.param<int>("vio/grid_n_height", grid_n_height, 17);                // 网格高度数
    nh.param<int>("vio/patch_pyrimid_level", patch_pyrimid_level, 3);     // 补丁金字塔层数
    nh.param<int>("vio/patch_size", patch_size, 8);                       // 补丁大小
    nh.param<double>("vio/outlier_threshold", outlier_threshold, 1000);   // 离群点阈值

    // 读取时间偏移和无人机相关参数
    nh.param<double>("time_offset/exposure_time_init", exposure_time_init, 0.0); // 初始曝光时间偏移
    nh.param<double>("time_offset/img_time_offset", img_time_offset, 0.0);       // 图像时间偏移
    nh.param<double>("time_offset/imu_time_offset", imu_time_offset, 0.0);       // IMU 时间偏移
    nh.param<bool>("uav/imu_rate_odom", imu_prop_enable, false);                 // 是否启用 IMU 传播里程计
    nh.param<bool>("uav/gravity_align_en", gravity_align_en, false);             // 是否启用重力对齐

    // 读取评估和 IMU 参数
    nh.param<string>("evo/seq_name", seq_name, "01");                    // 序列名称
    nh.param<bool>("evo/pose_output_en", pose_output_en, false);         // 是否输出姿态
    nh.param<double>("imu/gyr_cov", gyr_cov, 1.0);                       // 陀螺仪协方差
    nh.param<double>("imu/acc_cov", acc_cov, 1.0);                       // 加速度计协方差
    nh.param<int>("imu/imu_int_frame", imu_int_frame, 3);                // IMU 初始帧数
    nh.param<bool>("imu/imu_en", imu_en, false);                         // 是否启用 IMU
    nh.param<bool>("imu/gravity_est_en", gravity_est_en, true);          // 是否估计重力
    nh.param<bool>("imu/ba_bg_est_en", ba_bg_est_en, true);              // 是否估计加速度和陀螺仪偏置

    // 读取预处理参数
    nh.param<double>("preprocess/blind", p_pre->blind, 0.01);            // 盲区距离
    nh.param<double>("preprocess/filter_size_surf", filter_size_surf_min, 0.5); // 表面点滤波器大小
    nh.param<int>("preprocess/lidar_type", p_pre->lidar_type, AVIA);     // LiDAR 类型
    nh.param<int>("preprocess/scan_line", p_pre->N_SCANS, 6);            // 扫描线数
    nh.param<int>("preprocess/point_filter_num", p_pre->point_filter_num, 3); // 点云滤波数量
    nh.param<bool>("preprocess/feature_extract_enabled", p_pre->feature_enabled, false); // 是否启用特征提取

    // 读取 PCD 保存和外参校准参数
    nh.param<int>("pcd_save/interval", pcd_save_interval, -1);           // PCD 保存间隔
    nh.param<bool>("pcd_save/pcd_save_en", pcd_save_en, false);          // 是否启用 PCD 保存
    nh.param<bool>("pcd_save/colmap_output_en", colmap_output_en, false); // 是否启用 Colmap 输出
    nh.param<double>("pcd_save/filter_size_pcd", filter_size_pcd, 0.5);  // PCD 滤波器大小
    nh.param<vector<double>>("extrin_calib/extrinsic_T", extrinT, vector<double>()); // LiDAR 到 IMU 平移向量
    nh.param<vector<double>>("extrin_calib/extrinsic_R", extrinR, vector<double>()); // LiDAR 到 IMU 旋转矩阵
    nh.param<vector<double>>("extrin_calib/Pcl", cameraextrinT, vector<double>());   // 相机到 LiDAR 平移向量
    nh.param<vector<double>>("extrin_calib/Rcl", cameraextrinR, vector<double>());   // 相机到 LiDAR 旋转矩阵
    nh.param<double>("debug/plot_time", plot_time, -10);                 // 调试绘图时间
    nh.param<int>("debug/frame_cnt", frame_cnt, 6);                      // 调试帧数

    // 读取发布相关参数
    nh.param<double>("publish/blind_rgb_points", blind_rgb_points, 0.01); // RGB 点盲区距离
    nh.param<int>("publish/pub_scan_num", pub_scan_num, 1);              // 发布扫描次数
    nh.param<bool>("publish/pub_effect_point_en", pub_effect_point_en, false); // 是否发布有效点
    nh.param<bool>("publish/dense_map_en", dense_map_en, false);         // 是否启用密集地图

    // 计算盲区距离的平方
    p_pre->blind_sqr = p_pre->blind * p_pre->blind;
}

// 初始化系统组件
void LIVMapper::initializeComponents() 
{
    // 设置表面点降采样滤波器的叶大小
    downSizeFilterSurf.setLeafSize(filter_size_surf_min, filter_size_surf_min, filter_size_surf_min);

    // 将外参从向量形式转换为矩阵形式
    extT << VEC_FROM_ARRAY(extrinT);  // LiDAR 到 IMU 平移向量
    extR << MAT_FROM_ARRAY(extrinR);  // LiDAR 到 IMU 旋转矩阵

    // 将外参赋值给体素地图管理器
    voxelmap_manager->extT_ << VEC_FROM_ARRAY(extrinT);
    voxelmap_manager->extR_ << MAT_FROM_ARRAY(extrinR);

    // 从 ROS 命名空间加载相机模型，若失败抛出异常
    if (!vk::camera_loader::loadFromRosNs("laserMapping", vio_manager->cam)) throw std::runtime_error("Camera model not correctly specified.");

    // 配置 VIO 管理器参数
    vio_manager->grid_size = grid_size;                  // 网格大小
    vio_manager->patch_size = patch_size;                // 补丁大小
    vio_manager->outlier_threshold = outlier_threshold;  // 离群点阈值
    vio_manager->setImuToLidarExtrinsic(extT, extR);     // 设置 IMU 到 LiDAR 外参
    vio_manager->setLidarToCameraExtrinsic(cameraextrinR, cameraextrinT); // 设置 LiDAR 到相机外参
    vio_manager->state = &_state;                        // 设置状态指针
    vio_manager->state_propagat = &state_propagat;       // 设置传播状态指针
    vio_manager->max_iterations = max_iterations;        // 最大迭代次数
    vio_manager->img_point_cov = IMG_POINT_COV;          // 图像点协方差
    vio_manager->normal_en = normal_en;                  // 是否启用法向量
    vio_manager->inverse_composition_en = inverse_composition_en; // 是否启用逆合成
    vio_manager->raycast_en = raycast_en;                // 是否启用光线投射
    vio_manager->grid_n_width = grid_n_width;            // 网格宽度数
    vio_manager->grid_n_height = grid_n_height;          // 网格高度数
    vio_manager->patch_pyrimid_level = patch_pyrimid_level; // 补丁金字塔层数
    vio_manager->exposure_estimate_en = exposure_estimate_en; // 是否启用曝光估计
    vio_manager->colmap_output_en = colmap_output_en;    // 是否启用 Colmap 输出
    vio_manager->initializeVIO();                        // 初始化 VIO

    // 配置 IMU 处理模块
    p_imu->set_extrinsic(extT, extR);                    // 设置外参
    p_imu->set_gyr_cov_scale(V3D(gyr_cov, gyr_cov, gyr_cov)); // 设置陀螺仪协方差
    p_imu->set_acc_cov_scale(V3D(acc_cov, acc_cov, acc_cov)); // 设置加速度计协方差
    p_imu->set_inv_expo_cov(inv_expo_cov);               // 设置逆曝光协方差
    p_imu->set_gyr_bias_cov(V3D(0.0001, 0.0001, 0.0001)); // 设置陀螺仪偏置协方差
    p_imu->set_acc_bias_cov(V3D(0.0001, 0.0001, 0.0001)); // 设置加速度计偏置协方差
    p_imu->set_imu_init_frame_num(imu_int_frame);        // 设置 IMU 初始帧数

    // 根据参数禁用 IMU 相关功能
    if (!imu_en) p_imu->disable_imu();                   // 禁用 IMU
    if (!gravity_est_en) p_imu->disable_gravity_est();   // 禁用重力估计
    if (!ba_bg_est_en) p_imu->disable_bias_est();        // 禁用偏置估计
    if (!exposure_estimate_en) p_imu->disable_exposure_est(); // 禁用曝光估计

    // 根据启用状态设置 SLAM 模式
    slam_mode_ = (img_en && lidar_en) ? LIVO : imu_en ? ONLY_LIO : ONLY_LO; // LIVO, ONLY_LIO 或 ONLY_LO
}

// 初始化文件输出流
void LIVMapper::initializeFiles() 
{
    // 如果启用 PCD 保存和 Colmap 输出，执行脚本
    if (pcd_save_en && colmap_output_en)
    {
        // 设置 Colmap 输出脚本路径
        const std::string folderPath = std::string(ROOT_DIR) + "/scripts/colmap_output.sh";
        
        // 构造赋予脚本执行权限的命令
        std::string chmodCommand = "chmod +x " + folderPath;
        
        // 执行 chmod 命令，若失败则输出错误
        int chmodRet = system(chmodCommand.c_str());  
        if (chmodRet != 0) {
            std::cerr << "Failed to set execute permissions for the script." << std::endl;
            return;
        }

        // 执行脚本，若失败则输出错误
        int executionRet = system(folderPath.c_str());
        if (executionRet != 0) {
            std::cerr << "Failed to execute the script." << std::endl;
            return;
        }
    }

    // 打开 Colmap 点云输出文件
    if(colmap_output_en) fout_points.open(std::string(ROOT_DIR) + "Log/Colmap/sparse/0/points3D.txt", std::ios::out);
    // 打开 PCD 位置文件（JSON 格式）
    if(pcd_save_interval > 0) fout_pcd_pos.open(std::string(ROOT_DIR) + "Log/PCD/scans_pos.json", std::ios::out);
    // 打开调试文件（预处理矩阵）
    fout_pre.open(DEBUG_FILE_DIR("mat_pre.txt"), std::ios::out);
    // 打开调试文件（输出矩阵）
    fout_out.open(DEBUG_FILE_DIR("mat_out.txt"), std::ios::out);
}

// 初始化 ROS 订阅者和发布者
void LIVMapper::initializeSubscribersAndPublishers(ros::NodeHandle &nh, image_transport::ImageTransport &it) 
{
    // 订阅 LiDAR 数据，根据类型选择回调函数
    sub_pcl = p_pre->lidar_type == AVIA ? 
              nh.subscribe(lid_topic, 200000, &LIVMapper::livox_pcl_cbk, this): 
              nh.subscribe(lid_topic, 200000, &LIVMapper::standard_pcl_cbk, this);
    // 订阅 IMU 数据
    sub_imu = nh.subscribe(imu_topic, 200000, &LIVMapper::imu_cbk, this);
    // 订阅图像数据
    sub_img = nh.subscribe(img_topic, 200000, &LIVMapper::img_cbk, this);
    
    // 创建点云发布者
    pubLaserCloudFullRes = nh.advertise<sensor_msgs::PointCloud2>("/cloud_registered", 100); // 注册后的完整点云
    pubNormal = nh.advertise<visualization_msgs::MarkerArray>("visualization_marker", 100);   // 法向量标记
    pubSubVisualMap = nh.advertise<sensor_msgs::PointCloud2>("/cloud_visual_sub_map_before", 100); // 视觉子地图
    pubLaserCloudEffect = nh.advertise<sensor_msgs::PointCloud2>("/cloud_effected", 100);    // 有效点云
    pubLaserCloudMap = nh.advertise<sensor_msgs::PointCloud2>("/Laser_map", 100);            // 激光地图
    pubOdomAftMapped = nh.advertise<nav_msgs::Odometry>("/aft_mapped_to_init", 10);          // 优化后里程计
    pubPath = nh.advertise<nav_msgs::Path>("/path", 10);                                     // 路径
    plane_pub = nh.advertise<visualization_msgs::Marker>("/planner_normal", 1);              // 平面标记
    voxel_pub = nh.advertise<visualization_msgs::MarkerArray>("/voxels", 1);                 // 体素标记
    pubLaserCloudDyn = nh.advertise<sensor_msgs::PointCloud2>("/dyn_obj", 100);              // 动态对象点云
    pubLaserCloudDynRmed = nh.advertise<sensor_msgs::PointCloud2>("/dyn_obj_removed", 100);  // 移除动态对象后的点云
    pubLaserCloudDynDbg = nh.advertise<sensor_msgs::PointCloud2>("/dyn_obj_dbg_hist", 100);  // 动态对象调试历史
    mavros_pose_publisher = nh.advertise<geometry_msgs::PoseStamped>("/mavros/vision_pose/pose", 10); // MAVROS 姿态
    pubImage = it.advertise("/rgb_img", 1);                                                  // RGB 图像
    pubImuPropOdom = nh.advertise<nav_msgs::Odometry>("/LIVO2/imu_propagate", 10000);        // IMU 传播里程计

    // 创建定时器，每 0.004 秒调用 IMU 传播回调
    imu_prop_timer = nh.createTimer(ros::Duration(0.004), &LIVMapper::imu_prop_callback, this);
    // 为体素地图管理器设置平面发布者
    voxelmap_manager->voxel_map_pub_= nh.advertise<visualization_msgs::MarkerArray>("/planes", 10000);
}

// 处理首帧 LiDAR 数据
void LIVMapper::handleFirstFrame() 
{
    // 如果不是首帧
    if (!is_first_frame)
    {
        // 记录首帧 LiDAR 时间
        _first_lidar_time = LidarMeasures.last_lio_update_time;
        // 设置 IMU 处理模块的首帧时间（仅用于日志）
        p_imu->first_lidar_time = _first_lidar_time;
        // 标记为首帧
        is_first_frame = true;
        // 输出提示信息
        cout << "FIRST LIDAR FRAME!" << endl;
    }
}

// 重力对齐函数
void LIVMapper::gravityAlignment() 
{
    // 如果 IMU 初始化完成且重力对齐未完成
    if (!p_imu->imu_need_init && !gravity_align_finished) 
    {
        // 输出开始信息
        std::cout << "Gravity Alignment Starts" << std::endl;
        // 定义世界坐标系重力方向 (0, 0, -1)
        V3D ez(0, 0, -1), gz(_state.gravity);
        // 计算从当前重力到世界重力的四元数
        Quaterniond G_q_I0 = Quaterniond::FromTwoVectors(gz, ez);
        // 转换为旋转矩阵
        M3D G_R_I0 = G_q_I0.toRotationMatrix();

        // 应用旋转矩阵对齐状态
        _state.pos_end = G_R_I0 * _state.pos_end;    // 对齐位置
        _state.rot_end = G_R_I0 * _state.rot_end;    // 对齐旋转
        _state.vel_end = G_R_I0 * _state.vel_end;    // 对齐速度
        _state.gravity = G_R_I0 * _state.gravity;    // 对齐重力
        // 标记重力对齐完成
        gravity_align_finished = true;
        // 输出完成信息
        std::cout << "Gravity Alignment Finished" << std::endl;
    }
}

// 处理 IMU 数据
void LIVMapper::processImu() 
{
    // double t0 = omp_get_wtime(); // 记录开始时间（注释掉）

    // 调用 IMU 处理模块更新状态并去畸变点云
    p_imu->Process2(LidarMeasures, _state, feats_undistort);

    // 如果启用重力对齐，执行对齐
    if (gravity_align_en) gravityAlignment();

    // 更新传播状态
    state_propagat = _state;
    // 更新体素地图管理器的状态
    voxelmap_manager->state_ = _state;
    // 更新体素地图管理器的去畸变点云
    voxelmap_manager->feats_undistort_ = feats_undistort;

    // double t_prop = omp_get_wtime(); // 记录传播结束时间（注释掉）
    // std::cout << "[ Mapping ] feats_undistort: " << feats_undistort->size() << std::endl; // 输出去畸变点云数量（注释掉）
    // std::cout << "[ Mapping ] predict cov: " << _state.cov.diagonal().transpose() << std::endl; // 输出协方差（注释掉）
    // std::cout << "[ Mapping ] predict sta: " << state_propagat.pos_end.transpose() << state_propagat.vel_end.transpose() << std::endl; // 输出状态（注释掉）
}

// 状态估计和地图构建
void LIVMapper::stateEstimationAndMapping() 
{
    // 根据标志选择处理模式
    switch (LidarMeasures.lio_vio_flg) 
    {
        case VIO:
            // 处理视觉-惯性里程计
            handleVIO();
            break;
        case LIO:
        case LO:
            // 处理激光-惯性或仅激光里程计
            handleLIO();
            break;
    }
}

// 处理 VIO（视觉-惯性里程计）
void LIVMapper::handleVIO() 
{
    // 将旋转矩阵转换为欧拉角
    euler_cur = RotMtoEuler(_state.rot_end);
    // 记录当前状态到预处理文件
    fout_pre << std::setw(20) << LidarMeasures.last_lio_update_time - _first_lidar_time << " " << euler_cur.transpose() * 57.3 << " "
             << _state.pos_end.transpose() << " " << _state.vel_end.transpose() << " " << _state.bias_g.transpose() << " "
             << _state.bias_a.transpose() << " " << V3D(_state.inv_expo_time, 0, 0).transpose() << std::endl;
    
    // 检查点云是否为空
    if (pcl_w_wait_pub->empty() || (pcl_w_wait_pub == nullptr)) 
    {
        // 输出提示信息并返回
        std::cout << "[ VIO ] No point!!!" << std::endl;
        return;
    }
    
    // 输出原始特征点数量
    std::cout << "[ VIO ] Raw feature num: " << pcl_w_wait_pub->points.size() << std::endl;

    // 根据时间差设置是否绘制调试图
    if (fabs((LidarMeasures.last_lio_update_time - _first_lidar_time) - plot_time) < (frame_cnt / 2 * 0.1)) 
    {
        vio_manager->plot_flag = true; // 启用绘图
    } 
    else 
    {
        vio_manager->plot_flag = false; // 禁用绘图
    }

    // 处理当前帧
    vio_manager->processFrame(LidarMeasures.measures.back().img, _pv_list, voxelmap_manager->voxel_map_, LidarMeasures.last_lio_update_time - _first_lidar_time);

    // 如果启用 IMU 传播，更新 EKF 状态
    if (imu_prop_enable) 
    {
        ekf_finish_once = true;          // 标记 EKF 完成一次
        latest_ekf_state = _state;       // 更新最新 EKF 状态
        latest_ekf_time = LidarMeasures.last_lio_update_time; // 更新最新时间
        state_update_flg = true;         // 标记状态已更新
    }

    // 发布世界坐标系点云和 RGB 图像
    publish_frame_world(pubLaserCloudFullRes, vio_manager);
    publish_img_rgb(pubImage, vio_manager);

    // 将更新后的状态记录到输出文件
    euler_cur = RotMtoEuler(_state.rot_end);
    fout_out << std::setw(20) << LidarMeasures.last_lio_update_time - _first_lidar_time << " " << euler_cur.transpose() * 57.3 << " "
             << _state.pos_end.transpose() << " " << _state.vel_end.transpose() << " " << _state.bias_g.transpose() << " "
             << _state.bias_a.transpose() << " " << V3D(_state.inv_expo_time, 0, 0).transpose() << " " << feats_undistort->points.size() << std::endl;
}

// 处理 LIO（激光-惯性里程计）
void LIVMapper::handleLIO() 
{    
    // 将旋转矩阵转换为欧拉角
    euler_cur = RotMtoEuler(_state.rot_end);
    // 记录当前状态到预处理文件
    fout_pre << setw(20) << LidarMeasures.last_lio_update_time - _first_lidar_time << " " << euler_cur.transpose() * 57.3 << " "
             << _state.pos_end.transpose() << " " << _state.vel_end.transpose() << " " << _state.bias_g.transpose() << " "
             << _state.bias_a.transpose() << " " << V3D(_state.inv_expo_time, 0, 0).transpose() << endl;
           
    // 检查去畸变点云是否为空
    if (feats_undistort->empty() || (feats_undistort == nullptr)) 
    {
        // 输出提示信息并返回
        std::cout << "[ LIO ]: No point!!!" << std::endl;
        return;
    }

    // 记录开始时间
    double t0 = omp_get_wtime();

    // 对去畸变点云进行降采样
    downSizeFilterSurf.setInputCloud(feats_undistort);
    downSizeFilterSurf.filter(*feats_down_body);
    
    // 记录降采样结束时间
    double t_down = omp_get_wtime();

    // 更新体素地图管理器的点云数据
    feats_down_size = feats_down_body->points.size();        // 记录降采样点云大小
    voxelmap_manager->feats_down_body_ = feats_down_body;    // 设置机体坐标系点云
    transformLidar(_state.rot_end, _state.pos_end, feats_down_body, feats_down_world); // 转换为世界坐标系
    voxelmap_manager->feats_down_world_ = feats_down_world;  // 设置世界坐标系点云
    voxelmap_manager->feats_down_size_ = feats_down_size;    // 设置点云大小
    
    // 如果地图未初始化，构建体素地图
    if (!lidar_map_inited) 
    {
        lidar_map_inited = true;
        voxelmap_manager->BuildVoxelMap();
    }

    // 记录地图构建结束时间
    double t1 = omp_get_wtime();

    // 执行状态估计并更新状态
    voxelmap_manager->StateEstimation(state_propagat);
    _state = voxelmap_manager->state_;        // 更新状态
    _pv_list = voxelmap_manager->pv_list_;    // 更新点平面列表

    // 记录状态估计结束时间
    double t2 = omp_get_wtime();

    // 如果启用 IMU 传播，更新 EKF 状态
    if (imu_prop_enable) 
    {
        ekf_finish_once = true;          // 标记 EKF 完成一次
        latest_ekf_state = _state;       // 更新最新 EKF 状态
        latest_ekf_time = LidarMeasures.last_lio_update_time; // 更新最新时间
        state_update_flg = true;         // 标记状态已更新
    }

    // 如果启用姿态输出，保存到文件
    if (pose_output_en) 
    {
        static bool pos_opend = false;    // 文件是否已打开
        static int ocount = 0;            // 输出计数器
        std::ofstream outFile, evoFile;   // 文件流
        if (!pos_opend) 
        {
            // 首次打开文件
            evoFile.open(std::string(ROOT_DIR) + "Log/result/" + seq_name + ".txt", std::ios::out);
            pos_opend = true;
            if (!evoFile.is_open()) ROS_ERROR("open fail\n");
        } 
        else 
        {
            // 以追加模式打开文件
            evoFile.open(std::string(ROOT_DIR) + "Log/result/" + seq_name + ".txt", std::ios::app);
            if (!evoFile.is_open()) ROS_ERROR("open fail\n");
        }
        // 定义输出矩阵和四元数
        Eigen::Matrix4d outT;
        Eigen::Quaterniond q(_state.rot_end);
        evoFile << std::fixed;  // 设置固定格式输出
        // 输出时间、位置和四元数
        evoFile << LidarMeasures.last_lio_update_time << " " << _state.pos_end[0] << " " << _state.pos_end[1] << " " << _state.pos_end[2] << " "
                << q.x() << " " << q.y() << " " << q.z() << " " << q.w() << std::endl;
    }
    
    // 更新欧拉角并发布里程计
    euler_cur = RotMtoEuler(_state.rot_end);
    geoQuat = tf::createQuaternionMsgFromRollPitchYaw(euler_cur(0), euler_cur(1), euler_cur(2));
    publish_odometry(pubOdomAftMapped);

    // 记录里程计发布结束时间
    double t3 = omp_get_wtime();

    // 将点云转换到世界坐标系并更新体素地图
    PointCloudXYZI::Ptr world_lidar(new PointCloudXYZI());
    transformLidar(_state.rot_end, _state.pos_end, feats_down_body, world_lidar);
    for (size_t i = 0; i < world_lidar->points.size(); i++) 
    {
        // 更新点平面列表中的世界坐标
        voxelmap_manager->pv_list_[i].point_w << world_lidar->points[i].x, world_lidar->points[i].y, world_lidar->points[i].z;
        M3D point_crossmat = voxelmap_manager->cross_mat_list_[i];  // 交叉矩阵
        M3D var = voxelmap_manager->body_cov_list_[i];              // 协方差
        // 计算更新后的协方差
        var = (_state.rot_end * extR) * var * (_state.rot_end * extR).transpose() +
              (-point_crossmat) * _state.cov.block<3, 3>(0, 0) * (-point_crossmat).transpose() + _state.cov.block<3, 3>(3, 3);
        voxelmap_manager->pv_list_[i].var = var;  // 更新协方差
    }
    // 更新体素地图
    voxelmap_manager->UpdateVoxelMap(voxelmap_manager->pv_list_);
    std::cout << "[ LIO ] Update Voxel Map" << std::endl;
    _pv_list = voxelmap_manager->pv_list_;  // 更新点平面列表
    
    // 记录体素地图更新结束时间
    double t4 = omp_get_wtime();

    // 如果启用地图滑动，执行滑动操作
    if(voxelmap_manager->config_setting_.map_sliding_en)
    {
        voxelmap_manager->mapSliding();
    }
    
    // 选择点云（密集或降采样）
    PointCloudXYZI::Ptr laserCloudFullRes(dense_map_en ? feats_undistort : feats_down_body);
    int size = laserCloudFullRes->points.size();
    PointCloudXYZI::Ptr laserCloudWorld(new PointCloudXYZI(size, 1));

    // 将点云转换到世界坐标系
    for (int i = 0; i < size; i++) 
    {
        RGBpointBodyToWorld(&laserCloudFullRes->points[i], &laserCloudWorld->points[i]);
    }
    *pcl_w_wait_pub = *laserCloudWorld;  // 更新等待发布的点云

    // 发布相关数据
    if (!img_en) publish_frame_world(pubLaserCloudFullRes, vio_manager);  // 发布世界坐标系点云
    if (pub_effect_point_en) publish_effect_world(pubLaserCloudEffect, voxelmap_manager->ptpl_list_); // 发布有效点云
    if (voxelmap_manager->config_setting_.is_pub_plane_map_) voxelmap_manager->pubVoxelMap(); // 发布体素地图
    publish_path(pubPath);          // 发布路径
    publish_mavros(mavros_pose_publisher); // 发布 MAVROS 姿态

    // 更新帧数和平均耗时
    frame_num++;
    aver_time_consu = aver_time_consu * (frame_num - 1) / frame_num + (t4 - t0) / frame_num;

    // 输出时间统计信息
    printf("\033[1;34m+-------------------------------------------------------------+\033[0m\n");
    printf("\033[1;34m|                         LIO Mapping Time                    |\033[0m\n");
    printf("\033[1;34m+-------------------------------------------------------------+\033[0m\n");
    printf("\033[1;34m| %-29s | %-27s |\033[0m\n", "Algorithm Stage", "Time (secs)");
    printf("\033[1;34m+-------------------------------------------------------------+\033[0m\n");
    printf("\033[1;36m| %-29s | %-27f |\033[0m\n", "DownSample", t_down - t0);
    printf("\033[1;36m| %-29s | %-27f |\033[0m\n", "ICP", t2 - t1);
    printf("\033[1;36m| %-29s | %-27f |\033[0m\n", "updateVoxelMap", t4 - t3);
    printf("\033[1;34m+-------------------------------------------------------------+\033[0m\n");
    printf("\033[1;36m| %-29s | %-27f |\033[0m\n", "Current Total Time", t4 - t0);
    printf("\033[1;36m| %-29s | %-27f |\033[0m\n", "Average Total Time", aver_time_consu);
    printf("\033[1;34m+-------------------------------------------------------------+\033[0m\n");

    // 记录更新后的状态到输出文件
    euler_cur = RotMtoEuler(_state.rot_end);
    fout_out << std::setw(20) << LidarMeasures.last_lio_update_time - _first_lidar_time << " " << euler_cur.transpose() * 57.3 << " "
             << _state.pos_end.transpose() << " " << _state.vel_end.transpose() << " " << _state.bias_g.transpose() << " "
             << _state.bias_a.transpose() << " " << V3D(_state.inv_expo_time, 0, 0).transpose() << " " << feats_undistort->points.size() << std::endl;
}

// 保存 PCD 文件
void LIVMapper::savePCD() 
{
    // 如果启用 PCD 保存且有点云数据且不定期保存
    if (pcd_save_en && (pcl_wait_save->points.size() > 0 || pcl_wait_save_intensity->points.size() > 0) && pcd_save_interval < 0) 
    {
        // 定义保存路径
        std::string raw_points_dir = std::string(ROOT_DIR) + "Log/PCD/all_raw_points.pcd";
        std::string downsampled_points_dir = std::string(ROOT_DIR) + "Log/PCD/all_downsampled_points.pcd";
        pcl::PCDWriter pcd_writer;  // PCD 文件写入器

        // 如果启用图像处理
        if (img_en)
        {
            // 定义降采样点云
            pcl::PointCloud<pcl::PointXYZRGB>::Ptr downsampled_cloud(new pcl::PointCloud<pcl::PointXYZRGB>);
            pcl::VoxelGrid<pcl::PointXYZRGB> voxel_filter;  // 体素滤波器
            voxel_filter.setInputCloud(pcl_wait_save);      // 设置输入点云
            voxel_filter.setLeafSize(filter_size_pcd, filter_size_pcd, filter_size_pcd); // 设置滤波叶大小
            voxel_filter.filter(*downsampled_cloud);        // 执行滤波
            
            // 保存原始点云
            pcd_writer.writeBinary(raw_points_dir, *pcl_wait_save);
            std::cout << GREEN << "Raw point cloud data saved to: " << raw_points_dir 
                      << " with point count: " << pcl_wait_save->points.size() << RESET << std::endl;
            
            // 保存降采样点云
            pcd_writer.writeBinary(downsampled_points_dir, *downsampled_cloud);
            std::cout << GREEN << "Downsampled point cloud data saved to: " << downsampled_points_dir 
                      << " with point count after filtering: " << downsampled_cloud->points.size() << RESET << std::endl;

            // 如果启用 Colmap 输出
            if(colmap_output_en)
            {
                // 输出 Colmap 文件头
                fout_points << "# 3D point list with one line of data per point\n";
                fout_points << "#  POINT_ID, X, Y, Z, R, G, B, ERROR\n";
                // 遍历降采样点云并输出
                for (size_t i = 0; i < downsampled_cloud->size(); ++i) 
                {
                    const auto& point = downsampled_cloud->points[i];
                    fout_points << i << " "
                                << std::fixed << std::setprecision(6)
                                << point.x << " " << point.y << " " << point.z << " "
                                << static_cast<int>(point.r) << " "
                                << static_cast<int>(point.g) << " "
                                << static_cast<int>(point.b) << " "
                                << 0 << std::endl;
                }
            }
        }
        else
        {      
            // 保存带强度的原始点云
            pcd_writer.writeBinary(raw_points_dir, *pcl_wait_save_intensity);
            std::cout << GREEN << "Raw point cloud data saved to: " << raw_points_dir 
                      << " with point count: " << pcl_wait_save_intensity->points.size() << RESET << std::endl;
        }
    }
}

// 主运行函数
void LIVMapper::run() 
{
    // 设置 ROS 循环频率为 5000Hz
    ros::Rate rate(5000);
    while (ros::ok()) 
    {
        // 处理 ROS 回调
        ros::spinOnce();
        // 如果数据包未同步，等待
        if (!sync_packages(LidarMeasures)) 
        {
            rate.sleep();
            continue;
        }
        // 处理首帧
        handleFirstFrame();

        // 处理 IMU 数据
        processImu();

        // 执行状态估计和地图构建
        stateEstimationAndMapping();
    }
    // 保存 PCD 文件
    savePCD();
}

// 单次 IMU 传播
void LIVMapper::prop_imu_once(StatesGroup &imu_prop_state, const double dt, V3D acc_avr, V3D angvel_avr)
{
    // 获取 IMU 平均加速度范数
    double mean_acc_norm = p_imu->IMU_mean_acc_norm;
    // 校正加速度和角速度
    acc_avr = acc_avr * G_m_s2 / mean_acc_norm - imu_prop_state.bias_a;
    angvel_avr -= imu_prop_state.bias_g;

    // 计算旋转增量
    M3D Exp_f = Exp(angvel_avr, dt);
    // 传播 IMU 姿态
    imu_prop_state.rot_end = imu_prop_state.rot_end * Exp_f;

    // 计算全局坐标系下的加速度
    V3D acc_imu = imu_prop_state.rot_end * acc_avr + V3D(imu_prop_state.gravity[0], imu_prop_state.gravity[1], imu_prop_state.gravity[2]);

    // 传播 IMU 位置
    imu_prop_state.pos_end = imu_prop_state.pos_end + imu_prop_state.vel_end * dt + 0.5 * acc_imu * dt * dt;

    // 更新速度
    imu_prop_state.vel_end = imu_prop_state.vel_end + acc_imu * dt;
}

// IMU 传播回调函数
void LIVMapper::imu_prop_callback(const ros::TimerEvent &e)
{
    // 如果 IMU 未初始化或无新数据或 EKF 未完成，返回
    if (p_imu->imu_need_init || !new_imu || !ekf_finish_once) { return; }
    mtx_buffer_imu_prop.lock();  // 加锁
    new_imu = false;  // 重置新 IMU 数据标志
    // 如果启用 IMU 传播且缓冲区不为空
    if (imu_prop_enable && !prop_imu_buffer.empty())
    {
        static double last_t_from_lidar_end_time = 0;  // 上次时间差
        if (state_update_flg)
        {
            // 初始化传播状态
            imu_propagate = latest_ekf_state;
            // 丢弃无用的 IMU 数据包
            while ((!prop_imu_buffer.empty() && prop_imu_buffer.front().header.stamp.toSec() < latest_ekf_time))
            {
                prop_imu_buffer.pop_front();
            }
            last_t_from_lidar_end_time = 0;
            // 遍历 IMU 缓冲区进行传播
            for (int i = 0; i < prop_imu_buffer.size(); i++)
            {
                double t_from_lidar_end_time = prop_imu_buffer[i].header.stamp.toSec() - latest_ekf_time;
                double dt = t_from_lidar_end_time - last_t_from_lidar_end_time;
                V3D acc_imu(prop_imu_buffer[i].linear_acceleration.x, prop_imu_buffer[i].linear_acceleration.y, prop_imu_buffer[i].linear_acceleration.z);
                V3D omg_imu(prop_imu_buffer[i].angular_velocity.x, prop_imu_buffer[i].angular_velocity.y, prop_imu_buffer[i].angular_velocity.z);
                prop_imu_once(imu_propagate, dt, acc_imu, omg_imu);  // 单次传播
                last_t_from_lidar_end_time = t_from_lidar_end_time;
            }
            state_update_flg = false;  // 重置状态更新标志
        }
        else
        {
            // 使用最新 IMU 数据进行传播
            V3D acc_imu(newest_imu.linear_acceleration.x, newest_imu.linear_acceleration.y, newest_imu.linear_acceleration.z);
            V3D omg_imu(newest_imu.angular_velocity.x, newest_imu.angular_velocity.y, newest_imu.angular_velocity.z);
            double t_from_lidar_end_time = newest_imu.header.stamp.toSec() - latest_ekf_time;
            double dt = t_from_lidar_end_time - last_t_from_lidar_end_time;
            prop_imu_once(imu_propagate, dt, acc_imu, omg_imu);
            last_t_from_lidar_end_time = t_from_lidar_end_time;
        }

        // 发布传播后的里程计
        V3D posi, vel_i;
        Eigen::Quaterniond q;
        posi = imu_propagate.pos_end;
        vel_i = imu_propagate.vel_end;
        q = Eigen::Quaterniond(imu_propagate.rot_end);
        imu_prop_odom.header.frame_id = "world";
        imu_prop_odom.header.stamp = newest_imu.header.stamp;
        imu_prop_odom.pose.pose.position.x = posi.x();
        imu_prop_odom.pose.pose.position.y = posi.y();
        imu_prop_odom.pose.pose.position.z = posi.z();
        imu_prop_odom.pose.pose.orientation.w = q.w();
        imu_prop_odom.pose.pose.orientation.x = q.x();
        imu_prop_odom.pose.pose.orientation.y = q.y();
        imu_prop_odom.pose.pose.orientation.z = q.z();
        imu_prop_odom.twist.twist.linear.x = vel_i.x();
        imu_prop_odom.twist.twist.linear.y = vel_i.y();
        imu_prop_odom.twist.twist.linear.z = vel_i.z();
        pubImuPropOdom.publish(imu_prop_odom);
    }
    mtx_buffer_imu_prop.unlock();  // 解锁
}

// 将点云从机体坐标系转换到世界坐标系
void LIVMapper::transformLidar(const Eigen::Matrix3d rot, const Eigen::Vector3d t, const PointCloudXYZI::Ptr &input_cloud, PointCloudXYZI::Ptr &trans_cloud)
{
    // 清空目标点云并预分配空间
    PointCloudXYZI().swap(*trans_cloud);
    trans_cloud->reserve(input_cloud->size());
    // 遍历输入点云
    for (size_t i = 0; i < input_cloud->size(); i++)
    {
        pcl::PointXYZINormal p_c = input_cloud->points[i];
        Eigen::Vector3d p(p_c.x, p_c.y, p_c.z);
        // 应用外参和状态变换
        p = (rot * (extR * p + extT) + t);
        PointType pi;
        pi.x = p(0);
        pi.y = p(1);
        pi.z = p(2);
        pi.intensity = p_c.intensity;
        trans_cloud->points.push_back(pi);  // 添加转换后的点
    }
}

// 将单点从机体坐标系转换到世界坐标系
void LIVMapper::pointBodyToWorld(const PointType &pi, PointType &po)
{
    V3D p_body(pi.x, pi.y, pi.z);
    V3D p_global(_state.rot_end * (extR * p_body + extT) + _state.pos_end);
    po.x = p_global(0);
    po.y = p_global(1);
    po.z = p_global(2);
    po.intensity = pi.intensity;
}

// 模板函数：将点从机体坐标系转换到世界坐标系
template <typename T> void LIVMapper::pointBodyToWorld(const Matrix<T, 3, 1> &pi, Matrix<T, 3, 1> &po)
{
    V3D p_body(pi[0], pi[1], pi[2]);
    V3D p_global(_state.rot_end * (extR * p_body + extT) + _state.pos_end);
    po[0] = p_global(0);
    po[1] = p_global(1);
    po[2] = p_global(2);
}

// 模板函数：返回转换后的点
template <typename T> Matrix<T, 3, 1> LIVMapper::pointBodyToWorld(const Matrix<T, 3, 1> &pi)
{
    V3D p(pi[0], pi[1], pi[2]);
    p = (_state.rot_end * (extR * p + extT) + _state.pos_end);
    Matrix<T, 3, 1> po(p[0], p[1], p[2]);
    return po;
}

// 将 RGB 点从机体坐标系转换到世界坐标系
void LIVMapper::RGBpointBodyToWorld(PointType const *const pi, PointType *const po)
{
    V3D p_body(pi->x, pi->y, pi->z);
    V3D p_global(_state.rot_end * (extR * p_body + extT) + _state.pos_end);
    po->x = p_global(0);
    po->y = p_global(1);
    po->z = p_global(2);
    po->intensity = pi->intensity;
}

// 标准点云回调函数
void LIVMapper::standard_pcl_cbk(const sensor_msgs::PointCloud2::ConstPtr &msg)
{
    // 如果禁用 LiDAR，返回
    if (!lidar_en) return;
    mtx_buffer.lock();  // 加锁
    // 检查时间戳是否回退
    if (msg->header.stamp.toSec() < last_timestamp_lidar)
    {
        ROS_ERROR("lidar loop back, clear buffer");
        lid_raw_data_buffer.clear();
    }
    // 创建新的点云指针
    PointCloudXYZI::Ptr ptr(new PointCloudXYZI());
    p_pre->process(msg, ptr);  // 处理点云数据
    // 添加到缓冲区
    lid_raw_data_buffer.push_back(ptr);
    lid_header_time_buffer.push_back(msg->header.stamp.toSec());
    last_timestamp_lidar = msg->header.stamp.toSec();  // 更新最后时间戳

    mtx_buffer.unlock();  // 解锁
    sig_buffer.notify_all();  // 通知等待线程
}

// Livox 点云回调函数
void LIVMapper::livox_pcl_cbk(const livox_ros_driver::CustomMsg::ConstPtr &msg_in)
{
    // 如果禁用 LiDAR，返回
    if (!lidar_en) return;
    mtx_buffer.lock();  // 加锁
    // 复制消息
    livox_ros_driver::CustomMsg::Ptr msg(new livox_ros_driver::CustomMsg(*msg_in));
    // 检查 IMU 和 LiDAR 时间同步
    if (abs(last_timestamp_imu - msg->header.stamp.toSec()) > 1.0 && !imu_buffer.empty())
    {
        double timediff_imu_wrt_lidar = last_timestamp_imu - msg->header.stamp.toSec();
        printf("\033[95mSelf sync IMU and LiDAR, HARD time lag is %.10lf \n\033[0m", timediff_imu_wrt_lidar - 0.100);
    }

    // 获取当前时间戳
    double cur_head_time = msg->header.stamp.toSec();
    ROS_INFO("Get LiDAR, its header time: %.6f", cur_head_time);
    // 检查时间戳是否回退
    if (cur_head_time < last_timestamp_lidar)
    {
        ROS_ERROR("lidar loop back, clear buffer");
        lid_raw_data_buffer.clear();
    }
    // 处理点云数据
    PointCloudXYZI::Ptr ptr(new PointCloudXYZI());
    p_pre->process(msg, ptr);

    // 检查点云是否有效
    if (!ptr || ptr->empty()) {
        ROS_ERROR("Received an empty point cloud");
        mtx_buffer.unlock();
        return;
    }

    // 添加到缓冲区
    lid_raw_data_buffer.push_back(ptr);
    lid_header_time_buffer.push_back(cur_head_time);
    last_timestamp_lidar = cur_head_time;

    mtx_buffer.unlock();  // 解锁
    sig_buffer.notify_all();  // 通知等待线程
}

// IMU 回调函数
void LIVMapper::imu_cbk(const sensor_msgs::Imu::ConstPtr &msg_in)
{
    // 如果禁用 IMU，返回
    if (!imu_en) return;

    // 如果未收到 LiDAR 数据，返回
    if (last_timestamp_lidar < 0.0) return;
    // 复制消息并应用时间偏移
    sensor_msgs::Imu::Ptr msg(new sensor_msgs::Imu(*msg_in));
    msg->header.stamp = ros::Time().fromSec(msg->header.stamp.toSec() - imu_time_offset);
    double timestamp = msg->header.stamp.toSec();

    // 检查 IMU 和 LiDAR 是否同步
    if (fabs(last_timestamp_lidar - timestamp) > 0.5 && (!ros_driver_fix_en))
    {
        ROS_WARN("IMU and LiDAR not synced! delta time: %lf .\n", last_timestamp_lidar - timestamp);
    }

    // 如果启用 ROS 驱动修复，调整时间戳
    if (ros_driver_fix_en) timestamp += std::round(last_timestamp_lidar - timestamp);
    msg->header.stamp = ros::Time().fromSec(timestamp);

    mtx_buffer.lock();  // 加锁

    // 检查时间戳是否回退
    if (last_timestamp_imu > 0.0 && timestamp < last_timestamp_imu)
    {
        mtx_buffer.unlock();
        sig_buffer.notify_all();
        ROS_ERROR("imu loop back, offset: %lf \n", last_timestamp_imu - timestamp);
        return;
    }

    last_timestamp_imu = timestamp;  // 更新最后时间戳

    // 添加到 IMU 缓冲区
    imu_buffer.push_back(msg);
    mtx_buffer.unlock();  // 解锁
    // 如果启用 IMU 传播，添加到传播缓冲区
    if (imu_prop_enable)
    {
        mtx_buffer_imu_prop.lock();  // 加锁
        if (imu_prop_enable && !p_imu->imu_need_init) { prop_imu_buffer.push_back(*msg); }
        newest_imu = *msg;  // 更新最新 IMU 数据
        new_imu = true;     // 标记有新数据
        mtx_buffer_imu_prop.unlock();  // 解锁
    }
    sig_buffer.notify_all();  // 通知等待线程
}

// 从 ROS 消息获取图像
cv::Mat LIVMapper::getImageFromMsg(const sensor_msgs::ImageConstPtr &img_msg)
{
    // 转换为 OpenCV 格式（BGR8）
    cv::Mat img;
    img = cv_bridge::toCvCopy(img_msg, "bgr8")->image;
    return img;
}

// 图像回调函数
void LIVMapper::img_cbk(const sensor_msgs::ImageConstPtr &msg_in)
{
    // 如果禁用图像处理，返回
    if (!img_en) return;
    // 复制消息
    sensor_msgs::Image::Ptr msg(new sensor_msgs::Image(*msg_in));
    
    // 获取时间戳并应用偏移
    double msg_header_time = msg->header.stamp.toSec() + img_time_offset;
    // 检查时间戳重复
    if (abs(msg_header_time - last_timestamp_img) < 0.001) return;
    ROS_INFO("Get image, its header time: %.6f", msg_header_time);
    // 如果未收到 LiDAR 数据，返回
    if (last_timestamp_lidar < 0) return;

    // 检查时间戳是否回退
    if (msg_header_time < last_timestamp_img)
    {
        ROS_ERROR("image loop back. \n");
        return;
    }

    mtx_buffer.lock();  // 加锁

    // 校正图像时间
    double img_time_correct = msg_header_time;

    // 检查时间间隔是否过小
    if (img_time_correct - last_timestamp_img < 0.02)
    {
        ROS_WARN("Image need Jumps: %.6f", img_time_correct);
        mtx_buffer.unlock();
        sig_buffer.notify_all();
        return;
    }

    // 获取图像并添加到缓冲区
    cv::Mat img_cur = getImageFromMsg(msg);
    img_buffer.push_back(img_cur);
    img_time_buffer.push_back(img_time_correct);

    last_timestamp_img = img_time_correct;  // 更新最后时间戳
    mtx_buffer.unlock();  // 解锁
    sig_buffer.notify_all();  // 通知等待线程
}

// 数据包同步函数
bool LIVMapper::sync_packages(LidarMeasureGroup &meas)
{
    // 检查数据缓冲区是否为空
    if (lid_raw_data_buffer.empty() && lidar_en) return false;
    if (img_buffer.empty() && img_en) return false;
    if (imu_buffer.empty() && imu_en) return false;

    // 根据 SLAM 模式处理
    switch (slam_mode_)
    {
    case ONLY_LIO:
    {
        // 初始化最后更新时间
        if (meas.last_lio_update_time < 0.0) meas.last_lio_update_time = lid_header_time_buffer.front();
        if (!lidar_pushed)
        {
            // 将 LiDAR 数据推入测量缓冲区
            meas.lidar = lid_raw_data_buffer.front();
            if (meas.lidar->points.size() <= 1) return false;

            meas.lidar_frame_beg_time = lid_header_time_buffer.front(); // 设置开始时间
            meas.lidar_frame_end_time = meas.lidar_frame_beg_time + meas.lidar->points.back().curvature / double(1000); // 计算结束时间
            meas.pcl_proc_cur = meas.lidar;  // 设置当前处理点云
            lidar_pushed = true;             // 标记已推送
        }

        // 检查 IMU 数据是否足够
        if (imu_en && last_timestamp_imu < meas.lidar_frame_end_time)
        {
            return false;
        }

        // 创建测量组
        struct MeasureGroup m;
        m.imu.clear();
        m.lio_time = meas.lidar_frame_end_time;
        mtx_buffer.lock();  // 加锁
        // 添加 IMU 数据
        while (!imu_buffer.empty())
        {
            if (imu_buffer.front()->header.stamp.toSec() > meas.lidar_frame_end_time) break;
            m.imu.push_back(imu_buffer.front());
            imu_buffer.pop_front();
        }
        lid_raw_data_buffer.pop_front();  // 移除已处理数据
        lid_header_time_buffer.pop_front();
        mtx_buffer.unlock();  // 解锁
        sig_buffer.notify_all();  // 通知等待线程

        meas.lio_vio_flg = LIO;  // 设置为 LIO 模式
        meas.measures.push_back(m); // 添加测量组
        lidar_pushed = false;    // 重置推送标志
        return true;

        break;
    }

    case LIVO:
    {
        // LIVO 模式下，LIO 和 VIO 时间同步
        EKF_STATE last_lio_vio_flg = meas.lio_vio_flg;
        switch (last_lio_vio_flg)
        {
        case WAIT:
        case VIO:
        {
            // 计算图像捕获时间
            double img_capture_time = img_time_buffer.front() + exposure_time_init;
            if (meas.last_lio_update_time < 0.0) meas.last_lio_update_time = lid_header_time_buffer.front();

            // 获取最新时间
            double lid_newest_time = lid_header_time_buffer.back() + lid_raw_data_buffer.back()->points.back().curvature / double(1000);
            double imu_newest_time = imu_buffer.back()->header.stamp.toSec();

            // 检查图像时间是否有效
            if (img_capture_time < meas.last_lio_update_time + 0.00001)
            {
                img_buffer.pop_front();
                img_time_buffer.pop_front();
                ROS_ERROR("[ Data Cut ] Throw one image frame! \n");
                return false;
            }

            // 检查数据是否足够
            if (img_capture_time > lid_newest_time || img_capture_time > imu_newest_time)
            {
                return false;
            }

            // 创建测量组
            struct MeasureGroup m;
            m.imu.clear();
            m.lio_time = img_capture_time;
            mtx_buffer.lock();  // 加锁
            // 添加 IMU 数据
            while (!imu_buffer.empty())
            {
                if (imu_buffer.front()->header.stamp.toSec() > m.lio_time) break;
                if (imu_buffer.front()->header.stamp.toSec() > meas.last_lio_update_time) m.imu.push_back(imu_buffer.front());
                imu_buffer.pop_front();
            }
            mtx_buffer.unlock();  // 解锁
            sig_buffer.notify_all();  // 通知等待线程

            // 更新当前和下一帧点云
            *(meas.pcl_proc_cur) = *(meas.pcl_proc_next);
            PointCloudXYZI().swap(*meas.pcl_proc_next);

            int lid_frame_num = lid_raw_data_buffer.size();
            int max_size = meas.pcl_proc_cur->size() + 24000 * lid_frame_num;
            meas.pcl_proc_cur->reserve(max_size);
            meas.pcl_proc_next->reserve(max_size);

            // 分割点云数据
            while (!lid_raw_data_buffer.empty())
            {
                if (lid_header_time_buffer.front() > img_capture_time) break;
                auto pcl(lid_raw_data_buffer.front()->points);
                double frame_header_time(lid_header_time_buffer.front());
                float max_offs_time_ms = (m.lio_time - frame_header_time) * 1000.0f;

                for (int i = 0; i < pcl.size(); i++)
                {
                    auto pt = pcl[i];
                    if (pcl[i].curvature < max_offs_time_ms)
                    {
                        pt.curvature += (frame_header_time - meas.last_lio_update_time) * 1000.0f;
                        meas.pcl_proc_cur->points.push_back(pt);
                    }
                    else
                    {
                        pt.curvature += (frame_header_time - m.lio_time) * 1000.0f;
                        meas.pcl_proc_next->points.push_back(pt);
                    }
                }
                lid_raw_data_buffer.pop_front();
                lid_header_time_buffer.pop_front();
            }

            meas.measures.push_back(m);  // 添加测量组
            meas.lio_vio_flg = LIO;      // 设置为 LIO 模式
            return true;
        }

        case LIO:
        {
            // 处理 VIO
            double img_capture_time = img_time_buffer.front() + exposure_time_init;
            meas.lio_vio_flg = VIO;  // 设置为 VIO 模式
            meas.measures.clear();
            double imu_time = imu_buffer.front()->header.stamp.toSec();

            // 创建测量组
            struct MeasureGroup m;
            m.vio_time = img_capture_time;
            m.lio_time = meas.last_lio_update_time;
            m.img = img_buffer.front();
            mtx_buffer.lock();  // 加锁
            img_buffer.pop_front();
            img_time_buffer.pop_front();
            mtx_buffer.unlock();  // 解锁
            sig_buffer.notify_all();  // 通知等待线程
            meas.measures.push_back(m); // 添加测量组
            lidar_pushed = false;       // 重置推送标志
            return true;
        }

        default:
        {
            return false;
        }
        }
        break;
    }

    case ONLY_LO:
    {
        // 处理仅激光里程计
        if (!lidar_pushed) 
        { 
            if (lid_raw_data_buffer.empty()) return false;
            meas.lidar = lid_raw_data_buffer.front();
            meas.lidar_frame_beg_time = lid_header_time_buffer.front();
            meas.lidar_frame_end_time = meas.lidar_frame_beg_time + meas.lidar->points.back().curvature / double(1000);
            lidar_pushed = true;             
        }
        struct MeasureGroup m;
        m.lio_time = meas.lidar_frame_end_time;
        mtx_buffer.lock();  // 加锁
        lid_raw_data_buffer.pop_front();
        lid_header_time_buffer.pop_front();
        mtx_buffer.unlock();  // 解锁
        sig_buffer.notify_all();  // 通知等待线程
        lidar_pushed = false;     // 重置推送标志
        meas.lio_vio_flg = LO;    // 设置为 LO 模式
        meas.measures.push_back(m); // 添加测量组
        return true;
        break;
    }

    default:
    {
        printf("!! WRONG SLAM TYPE !!");
        return false;
    }
    }
    ROS_ERROR("out sync");
}

/**
 * @brief 发布 RGB 图像
 * 
 * 该函数用于从 VIOManager 获取当前的 RGB 图像，并以 ROS 兼容格式发布。
 * 适用于视觉 SLAM（VIO/LIVO）系统，用于可视化或后续处理。
 * 
 * @param pubImage     ROS 图像发布器（image_transport），用于发布图像数据
 * @param vio_manager  视觉-惯性-里程计（VIO）管理器，提供当前 RGB 图像
 */
void LIVMapper::publish_img_rgb(const image_transport::Publisher &pubImage, VIOManagerPtr vio_manager)
{
    // 🟢 1️⃣ 获取最新的 RGB 图像
    cv::Mat img_rgb = vio_manager->img_cp; 

    // 🟢 2️⃣ 创建 ROS 兼容的图像消息
    cv_bridge::CvImage out_msg;
    out_msg.header.stamp = ros::Time::now();  // 设置时间戳
    out_msg.encoding = sensor_msgs::image_encodings::BGR8; // 设置编码格式
    out_msg.image = img_rgb; // 赋值图像数据

    // 🟢 3️⃣ 发布 ROS 图像消息
    pubImage.publish(out_msg.toImageMsg());
}


// 发布世界坐标系点云
/**
 * @brief 发布世界坐标系下的点云帧
 * 
 * 该函数负责发布 LiDAR 扫描的点云数据，支持将 LiDAR 点云与相机图像融合，以生成带颜色信息的 RGB 点云。
 * 并且，若启用了 PCD 存储功能，则可以保存点云到磁盘。
 * 
 * @param pubLaserCloudFullRes  ROS 发布器，用于发布点云
 * @param vio_manager           视觉-惯性-里程计（VIO）管理器，用于获取 RGB 图像数据
 */
void LIVMapper::publish_frame_world(const ros::Publisher &pubLaserCloudFullRes, VIOManagerPtr vio_manager)
{
    // 🟢 1️⃣ 检查待发布点云是否为空
    if (pcl_w_wait_pub->empty()) return;

    // 创建一个 RGB 点云对象
    PointCloudXYZRGB::Ptr laserCloudWorldRGB(new PointCloudXYZRGB());

    // 🟢 2️⃣ 处理图像信息（如果启用了视觉 VIO）
    if (img_en)
    {
        static int pub_num = 1;  // 计数器，控制发布频率
        *pcl_wait_pub += *pcl_w_wait_pub; // 累加点云数据

        if (pub_num == pub_scan_num) // 达到设定的发布间隔
        {
            pub_num = 1;  // 重置计数器
            size_t size = pcl_wait_pub->points.size();
            laserCloudWorldRGB->reserve(size);

            // 获取最新的 RGB 图像
            cv::Mat img_rgb = vio_manager->img_rgb;

            // 遍历所有点云，将 LiDAR 点投影到图像上，并获取颜色信息
            for (size_t i = 0; i < size; i++)
            {
                PointTypeRGB pointRGB;
                pointRGB.x = pcl_wait_pub->points[i].x;
                pointRGB.y = pcl_wait_pub->points[i].y;
                pointRGB.z = pcl_wait_pub->points[i].z;

                // 🟡 3️⃣ 计算点云在相机坐标系中的位置
                V3D p_w(pcl_wait_pub->points[i].x, pcl_wait_pub->points[i].y, pcl_wait_pub->points[i].z);
                V3D pf(vio_manager->new_frame_->w2f(p_w)); 
                if (pf[2] < 0) continue; // 丢弃相机后方的点

                // 🟡 4️⃣ 计算像素坐标
                V2D pc(vio_manager->new_frame_->w2c(p_w));

                // 确保点落在相机图像范围内
                if (vio_manager->new_frame_->cam_->isInFrame(pc.cast<int>(), 3))
                {
                    // 获取插值后的像素值（RGB）
                    V3F pixel = vio_manager->getInterpolatedPixel(img_rgb, pc);
                    pointRGB.r = pixel[2];
                    pointRGB.g = pixel[1];
                    pointRGB.b = pixel[0];

                    // 只保留距离大于 blind_rgb_points 的点
                    if (pf.norm() > blind_rgb_points) laserCloudWorldRGB->push_back(pointRGB);
                }
            }
        }
        else
        {
            pub_num++;  // 增加计数器
        }
    }

    // 🟢 5️⃣ 发布 ROS 点云消息
    sensor_msgs::PointCloud2 laserCloudmsg;
    if (img_en)
    {
        pcl::toROSMsg(*laserCloudWorldRGB, laserCloudmsg); // 转换 RGB 点云
    }
    else 
    { 
        pcl::toROSMsg(*pcl_w_wait_pub, laserCloudmsg); // 转换无颜色的点云
    }
    laserCloudmsg.header.stamp = ros::Time::now(); // 设置时间戳
    laserCloudmsg.header.frame_id = "camera_init"; // 设置坐标系
    pubLaserCloudFullRes.publish(laserCloudmsg);  // 发布点云

    // 🟢 6️⃣ 处理 PCD 存储
    if (pcd_save_en)
    {
        int size = feats_undistort->points.size();
        PointCloudXYZI::Ptr laserCloudWorld(new PointCloudXYZI(size, 1));
        static int scan_wait_num = 0;  // 计数器

        if (img_en)
        {
            *pcl_wait_save += *laserCloudWorldRGB; // 保存 RGB 点云
        }
        else
        {
            *pcl_wait_save_intensity += *pcl_w_wait_pub; // 保存无颜色点云
        }
        scan_wait_num++;

        // 如果达到设定的保存间隔，则将点云存入文件
        if ((pcl_wait_save->size() > 0 || pcl_wait_save_intensity->size() > 0) && pcd_save_interval > 0 && scan_wait_num >= pcd_save_interval)
        {
            pcd_index++;  // 增加索引
            string all_points_dir = string(ROOT_DIR) + "Log/PCD/" + to_string(pcd_index) + ".pcd";
            pcl::PCDWriter pcd_writer;

            if (pcd_save_en)
            {
                cout << "current scan saved to /PCD/" << all_points_dir << endl;
                if (img_en)
                {
                    pcd_writer.writeBinary(all_points_dir, *pcl_wait_save); // 保存 RGB 点云
                    PointCloudXYZRGB().swap(*pcl_wait_save); // 清空缓存
                }
                else
                {
                    pcd_writer.writeBinary(all_points_dir, *pcl_wait_save_intensity); // 保存带强度信息的点云
                    PointCloudXYZI().swap(*pcl_wait_save_intensity); // 清空缓存
                }        

                // 保存位姿数据
                Eigen::Quaterniond q(_state.rot_end);
                fout_pcd_pos << _state.pos_end[0] << " " << _state.pos_end[1] << " " << _state.pos_end[2] << " "
                             << q.w() << " " << q.x() << " " << q.y() << " " << q.z() << endl;
                scan_wait_num = 0;  // 重置计数器
            }
        }
    }

    // 🟢 7️⃣ 清空缓存，防止点云重复发布
    if (laserCloudWorldRGB->size() > 0) 
        PointCloudXYZI().swap(*pcl_wait_pub); 

    PointCloudXYZI().swap(*pcl_w_wait_pub);
}


/**
 * @brief 发布视觉子地图（Visual Submap）
 * 
 * 该函数用于发布当前的视觉子地图点云。视觉子地图通常用于
 * 视觉-惯性-激光 SLAM（LIVO）或纯视觉 SLAM（VIO），
 * 作为局部地图的一部分，用于特征匹配、回环检测等。
 * 
 * @param pubSubVisualMap  ROS 发布器，用于发布视觉子地图点云
 */
void LIVMapper::publish_visual_sub_map(const ros::Publisher &pubSubVisualMap)
{
    // 获取视觉子地图点云
    PointCloudXYZI::Ptr laserCloudFullRes(visual_sub_map);

    // 获取点云的大小，如果为空，则直接返回
    int size = laserCloudFullRes->points.size();
    if (size == 0) return; 

    // 创建新的点云指针，用于存储需要发布的子地图点云
    PointCloudXYZI::Ptr sub_pcl_visual_map_pub(new PointCloudXYZI());

    // 复制视觉子地图数据
    *sub_pcl_visual_map_pub = *laserCloudFullRes;

    /*** 发布子地图 ***/
    if (1)  // 这里的 `if (1)` 实际没有意义，可能用于调试
    {
        sensor_msgs::PointCloud2 laserCloudmsg;

        // 将 PCL 点云转换为 ROS 点云消息
        pcl::toROSMsg(*sub_pcl_visual_map_pub, laserCloudmsg);

        // 设置时间戳，确保数据的时序同步
        laserCloudmsg.header.stamp = ros::Time::now();

        // 设置坐标系，"camera_init" 代表点云的参考框架
        laserCloudmsg.header.frame_id = "camera_init";

        // 通过 ROS 话题发布视觉子地图点云
        pubSubVisualMap.publish(laserCloudmsg);
    }
}


/**
 * @brief 发布影响点云（Effect Feature Point Cloud）
 * 
 * 该函数用于发布某些特定的点云数据，例如关键点、特征点等。通常在 LiDAR-IMU 里程计（LIO）或
 * 视觉-惯性-激光 SLAM（LIVO）系统中，这些点用于优化、匹配或地图构建。
 * 
 * @param pubLaserCloudEffect  ROS 发布器，用于发布点云
 * @param ptpl_list            存储点到平面的特征点信息（PointToPlane 类型）
 */
void LIVMapper::publish_effect_world(const ros::Publisher &pubLaserCloudEffect, const std::vector<PointToPlane> &ptpl_list)
{
    // 获取特征点云的数量
    int effect_feat_num = ptpl_list.size();

    // 创建 PCL 点云指针，并为点云分配空间
    PointCloudXYZI::Ptr laserCloudWorld(new PointCloudXYZI(effect_feat_num, 1));

    // 遍历 ptpl_list，将点云数据转换为 PCL 格式
    for (int i = 0; i < effect_feat_num; i++)
    {
        laserCloudWorld->points[i].x = ptpl_list[i].point_w_[0];  // X 坐标
        laserCloudWorld->points[i].y = ptpl_list[i].point_w_[1];  // Y 坐标
        laserCloudWorld->points[i].z = ptpl_list[i].point_w_[2];  // Z 坐标
        // 注意：此处没有设置点云的强度 (intensity)，如果需要可以修改
    }

    /*** 将 PCL 格式的点云转换为 ROS 兼容的消息 ***/
    sensor_msgs::PointCloud2 laserCloudFullRes3;
    pcl::toROSMsg(*laserCloudWorld, laserCloudFullRes3);  // 转换 PCL -> ROS

    // 设置时间戳，确保数据的时序同步
    laserCloudFullRes3.header.stamp = ros::Time::now();

    // 设置点云坐标系，"camera_init" 代表点云的参考框架
    laserCloudFullRes3.header.frame_id = "camera_init";

    // 发布点云到 ROS 话题
    pubLaserCloudEffect.publish(laserCloudFullRes3);
}


/**
 * @brief 设置位姿戳（PoseStamped 数据填充）
 * 
 * 该函数用于填充位置和姿态信息，将 `_state` 的最新位置信息
 * 以及 `geoQuat` 四元数旋转赋值给输出的 `out` 变量。
 * 
 * @tparam T  目标数据类型（通常是 geometry_msgs::Pose 或 PoseStamped 类型）
 * @param out 输出的位姿对象
 */
template <typename T> 
void LIVMapper::set_posestamp(T &out)
{
    // 设置位置坐标（来自状态估计 _state）
    out.position.x = _state.pos_end(0);
    out.position.y = _state.pos_end(1);
    out.position.z = _state.pos_end(2);

    // 设置四元数旋转（来自 geoQuat）
    out.orientation.x = geoQuat.x;
    out.orientation.y = geoQuat.y;
    out.orientation.z = geoQuat.z;
    out.orientation.w = geoQuat.w;
}


/**
 * @brief 发布里程计（Odometry）数据
 * 
 * 该函数负责发布当前的位姿信息，并在 ROS 坐标系中广播 TF 变换，以支持机器人或无人机的导航。
 * 
 * @param pubOdomAftMapped  ROS 发布器，发布里程计消息
 */
void LIVMapper::publish_odometry(const ros::Publisher &pubOdomAftMapped)
{
    // 设置里程计的全局坐标系（父坐标系）
    odomAftMapped.header.frame_id = "camera_init";

    // 设置子坐标系（用于描述里程计估计的坐标系）
    odomAftMapped.child_frame_id = "aft_mapped";

    // 设置时间戳，确保数据是最新的
    odomAftMapped.header.stamp = ros::Time::now();

    // 填充当前位姿信息
    set_posestamp(odomAftMapped.pose.pose);

    /*** 发布 TF 变换 ***/
    static tf::TransformBroadcaster br;  // TF 变换广播器

    tf::Transform transform;  // 变换矩阵
    tf::Quaternion q;         // 四元数表示的旋转

    // 设置位置信息（来自状态估计 _state）
    transform.setOrigin(tf::Vector3(
        _state.pos_end(0), 
        _state.pos_end(1), 
        _state.pos_end(2)
    ));

    // 设置旋转信息（使用四元数表示）
    q.setW(geoQuat.w);
    q.setX(geoQuat.x);
    q.setY(geoQuat.y);
    q.setZ(geoQuat.z);
    transform.setRotation(q);

    // 发送 TF 变换，将 "camera_init" 作为参考坐标系，将 "aft_mapped" 作为当前估计的位姿
    br.sendTransform(tf::StampedTransform(transform, odomAftMapped.header.stamp, "camera_init", "aft_mapped"));

    // 发布 Odometry 消息，以便其他 ROS 节点接收并使用
    pubOdomAftMapped.publish(odomAftMapped);
}


// 发布 MAVROS 姿态
/**
 * @brief 发布 MAVROS 兼容的位姿信息
 * 
 * 该函数用于向 MAVROS 发送位姿 (PoseStamped) 数据，使无人机能够使用
 * 视觉里程计数据进行导航或控制。
 * 
 * @param mavros_pose_publisher  ROS 发布器，用于发布位姿信息
 */
void LIVMapper::publish_mavros(const ros::Publisher &mavros_pose_publisher)
{
    // 设置消息的时间戳，使用当前 ROS 时间
    msg_body_pose.header.stamp = ros::Time::now();

    // 设置坐标系的参考框架，这里使用 "camera_init" 作为起始坐标系
    msg_body_pose.header.frame_id = "camera_init";

    // 设置位姿信息，将当前状态的位姿赋值给 msg_body_pose.pose
    set_posestamp(msg_body_pose.pose);

    // 通过 ROS 话题发布该位姿信息，以便 MAVROS 或其他组件接收和使用
    mavros_pose_publisher.publish(msg_body_pose);
}


//

// 发布路径的函数
void LIVMapper::publish_path(const ros::Publisher pubPath)
{
    // 设置当前姿态（位置和方向）到消息体中
    set_posestamp(msg_body_pose.pose);  // 调用模板函数将当前状态 (_state.pos_end 和 geoQuat) 设置到 msg_body_pose 的 pose 字段

    // 更新消息的时间戳为当前时间
    msg_body_pose.header.stamp = ros::Time::now();  // 设置消息头的时间戳为当前 ROS 时间

    // 设置消息的参考坐标系为 "camera_init"
    msg_body_pose.header.frame_id = "camera_init";  // 指定路径的坐标系为初始相机坐标系

    // 将当前姿态消息添加到路径的姿态列表中
    path.poses.push_back(msg_body_pose);  // 将当前的 msg_body_pose 添加到 path 的 poses 向量中，构建路径

    // 通过 ROS 发布路径消息
    pubPath.publish(path);  // 使用传入的发布者 pubPath 发布整个路径消息 (nav_msgs::Path 类型)
}
