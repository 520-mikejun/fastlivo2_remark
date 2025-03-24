#include "vio.h"

// VIOManager 类的构造函数
VIOManager::VIOManager()
{
    // 可选：取消注释以初始化点云降采样滤波器，体素大小为 x、y、z 方向上的 0.2 米
    // downSizeFilter.setLeafSize(0.2, 0.2, 0.2);
}

// VIOManager 类的析构函数
VIOManager::~VIOManager()
{
    // 清理动态分配的内存
    delete visual_submap;  // 删除视觉子地图对象
    for (auto& pair : warp_map) delete pair.second;  // 删除 warp_map 中的所有 warp 对象
    warp_map.clear();  // 清空 warp_map 容器
    for (auto& pair : feat_map) delete pair.second;  // 删除 feat_map 中的所有特征对象
    feat_map.clear();  // 清空 feat_map 容器
}

// 设置 IMU 到激光雷达的外参变换
// 输入：transl (平移向量), rot (旋转矩阵)
void VIOManager::setImuToLidarExtrinsic(const V3D &transl, const M3D &rot)
{
    Pli = -rot.transpose() * transl;  // 计算从激光雷达到 IMU 坐标系的平移
    Rli = rot.transpose();  // 存储从激光雷达到 IMU 坐标系的逆旋转
}

// 设置激光雷达到相机的外参变换
// 输入：R (旋转向量), P (平移向量)
void VIOManager::setLidarToCameraExtrinsic(vector<double> &R, vector<double> &P)
{
    Rcl << MAT_FROM_ARRAY(R);  // 将旋转向量转换为矩阵
    Pcl << VEC_FROM_ARRAY(P);  // 将平移向量转换为 Eigen 向量
}

// 初始化视觉-惯性里程计系统
void VIOManager::initializeVIO()
{
    visual_submap = new SubSparseMap;  // 分配一个新的稀疏视觉子地图

    // 获取相机内参
    fx = cam->fx();  // x 方向焦距
    fy = cam->fy();  // y 方向焦距
    cx = cam->cx();  // 主点 x 坐标
    cy = cam->cy();  // 主点 y 坐标
    image_resize_factor = cam->scale();  // 图像缩放因子

    printf("内参: %.6lf, %.6lf, %.6lf, %.6lf\n", fx, fy, cx, cy);  // 记录内参信息

    width = cam->width();  // 图像宽度
    height = cam->height();  // 图像高度

    printf("宽度: %d, 高度: %d, 缩放因子: %f\n", width, height, image_resize_factor);  // 记录图像尺寸和缩放因子

    Rci = Rcl * Rli;  // 计算相机到 IMU 的旋转矩阵
    Pci = Rcl * Pli + Pcl;  // 计算相机到 IMU 的平移向量

    V3D Pic;  // 相机在 IMU 坐标系中的位置
    M3D tmp;  // 临时矩阵用于计算雅可比
    Jdphi_dR = Rci;  // 初始化旋转对角度的雅可比矩阵
    Pic = -Rci.transpose() * Pci;  // 计算相机相对于 IMU 的位置
    tmp << SKEW_SYM_MATRX(Pic);  // 生成 Pic 的反对称矩阵
    Jdp_dR = -Rci * tmp;  // 计算位置对旋转的雅可比矩阵

    // 根据网格大小初始化网格维度
    if (grid_size > 10)
    {
        grid_n_width = ceil(static_cast<double>(width / grid_size));  // 计算网格宽度
        grid_n_height = ceil(static_cast<double>(height / grid_size));  // 计算网格高度
    }
    else
    {
        grid_size = static_cast<int>(height / grid_n_height);  // 根据高度调整网格大小
        grid_n_height = ceil(static_cast<double>(height / grid_size));  // 重新计算网格高度
        grid_n_width = ceil(static_cast<double>(width / grid_size));  // 重新计算网格宽度
    }
    length = grid_n_width * grid_n_height;  // 计算总网格数量

    // 如果启用了光线投射功能
    if (raycast_en)
    {
        // cv::Mat img_test = cv::Mat::zeros(height, width, CV_8UC1);  // 可选：创建测试图像
        // uchar* it = (uchar*)img_test.data;  // 指向测试图像数据的指针

        border_flag.resize(length, 0);  // 初始化边界标志向量，所有值为 0

        std::vector<std::vector<V3D>>().swap(rays_with_sample_points);  // 清空并重新分配光线采样点容器
        rays_with_sample_points.reserve(length);  // 为每个网格预留空间
        printf("网格大小: %d, 网格高度: %d, 网格宽度: %d, 总网格数: %d\n", grid_size, grid_n_height, grid_n_width, length);

        float d_min = 0.1;  // 最小采样深度
        float d_max = 3.0;  // 最大采样深度
        float step = 0.2;  // 采样步长
        for (int grid_row = 1; grid_row <= grid_n_height; grid_row++)  // 遍历所有网格行
        {
            for (int grid_col = 1; grid_col <= grid_n_width; grid_col++)  // 遍历所有网格列
            {
                std::vector<V3D> SamplePointsEachGrid;  // 存储每个网格的采样点
                int index = (grid_row - 1) * grid_n_width + grid_col - 1;  // 计算网格索引

                // 如果是边界网格，标记为 1
                if (grid_row == 1 || grid_col == 1 || grid_row == grid_n_height || grid_col == grid_n_width) 
                    border_flag[index] = 1;

                int u = grid_size / 2 + (grid_col - 1) * grid_size;  // 计算网格中心的像素 u 坐标
                int v = grid_size / 2 + (grid_row - 1) * grid_size;  // 计算网格中心的像素 v 坐标
                // it[ u + v * width ] = 255;  // 可选：在测试图像中标记网格中心

                // 沿深度方向采样
                for (float d_temp = d_min; d_temp <= d_max; d_temp += step)
                {
                    V3D xyz;  // 采样点在相机坐标系中的 3D 坐标
                    xyz = cam->cam2world(u, v);  // 将像素坐标转换为相机坐标系下的单位向量
                    xyz *= d_temp / xyz[2];  // 根据深度缩放坐标
                    SamplePointsEachGrid.push_back(xyz);  // 添加采样点到当前网格
                }
                rays_with_sample_points.push_back(SamplePointsEachGrid);  // 将采样点添加到光线容器
            }
        }
        // 可选：调试输出采样点信息并显示测试图像
        // printf("光线采样点数量: %d, 容量: %d, 第一网格容量: %d, 第一网格采样点数: %d\n",
        //        rays_with_sample_points.size(), rays_with_sample_points.capacity(),
        //        rays_with_sample_points[0].capacity(), rays_with_sample_points[0].size());
        // cv::imshow("img_test", img_test);
        // cv::waitKey(1);
    }

    // 如果启用了 Colmap 输出功能
    if (colmap_output_en)
    {
        pinhole_cam = dynamic_cast<vk::PinholeCamera*>(cam);  // 将相机转换为针孔相机模型
        fout_colmap.open(DEBUG_FILE_DIR("Colmap/sparse/0/images.txt"), ios::out);  // 打开 Colmap 图像输出文件
        fout_colmap << "# 图像列表，每张图像包含两行数据:\n";  // 写入文件头注释
        fout_colmap << "#   IMAGE_ID, QW, QX, QY, QZ, TX, TY, TZ, CAMERA_ID, NAME\n";  // 描述图像数据格式
        fout_colmap << "#   POINTS2D[] as (X, Y, POINT3D_ID)\n";  // 描述二维点数据格式

        fout_camera.open(DEBUG_FILE_DIR("Colmap/sparse/0/cameras.txt"), ios::out);  // 打开 Colmap 相机输出文件
        fout_camera << "# 相机列表，每相机一行数据:\n";  // 写入文件头注释
        fout_camera << "#   CAMERA_ID, MODEL, WIDTH, HEIGHT, PARAMS[]\n";  // 描述相机数据格式
        fout_camera << "1 PINHOLE " << width << " " << height << " "
                    << std::fixed << std::setprecision(6)  // 控制浮点数精度为 6 位
                    << fx << " " << fy << " "
                    << cx << " " << cy << std::endl;  // 写入相机参数
        fout_camera.close();  // 关闭相机文件
    }

    // 初始化各种容器
    grid_num.resize(length);  // 网格类型向量
    map_index.resize(length);  // 地图索引向量
    map_dist.resize(length);  // 地图距离向量
    update_flag.resize(length);  // 更新标志向量
    scan_value.resize(length);  // 扫描值向量

    patch_size_total = patch_size * patch_size;  // 计算补丁的总像素数
    patch_size_half = static_cast<int>(patch_size / 2);  // 计算补丁半宽
    patch_buffer.resize(patch_size_total);  // 初始化补丁缓冲区
    warp_len = patch_size_total * patch_pyrimid_level;  // 计算变换补丁的总长度
    border = (patch_size_half + 1) * (1 << patch_pyrimid_level);  // 计算边界宽度

    retrieve_voxel_points.reserve(length);  // 为检索体素点预留空间
    append_voxel_points.reserve(length);  // 为追加体素点预留空间

    sub_feat_map.clear();  // 清空子特征地图
}

// 重置网格状态
void VIOManager::resetGrid()
{
    fill(grid_num.begin(), grid_num.end(), TYPE_UNKNOWN);  // 将所有网格类型设置为未知
    fill(map_index.begin(), map_index.end(), 0);  // 将所有地图索引清零
    fill(map_dist.begin(), map_dist.end(), 10000.0f);  // 将所有地图距离设置为初始值 10000
    fill(update_flag.begin(), update_flag.end(), 0);  // 将所有更新标志清零
    fill(scan_value.begin(), scan_value.end(), 0.0f);  // 将所有扫描值清零

    retrieve_voxel_points.clear();  // 清空检索体素点容器
    retrieve_voxel_points.resize(length);  // 调整大小以匹配网格数量

    append_voxel_points.clear();  // 清空追加体素点容器
    append_voxel_points.resize(length);  // 调整大小以匹配网格数量

    total_points = 0;  // 重置总点数为 0
}

// 计算投影雅可比矩阵
// 输入：p (3D 点坐标)
// 输出：J (2x3 雅可比矩阵)
void VIOManager::computeProjectionJacobian(V3D p, MD(2, 3) & J)
{
    const double x = p[0];  // 点的 x 坐标
    const double y = p[1];  // 点的 y 坐标
    const double z_inv = 1. / p[2];  // 深度倒数
    const double z_inv_2 = z_inv * z_inv;  // 深度倒数的平方
    J(0, 0) = fx * z_inv;  // u 对 x 的偏导数
    J(0, 1) = 0.0;  // u 对 y 的偏导数
    J(0, 2) = -fx * x * z_inv_2;  // u 对 z 的偏导数
    J(1, 0) = 0.0;  // v 对 x 的偏导数
    J(1, 1) = fy * z_inv;  // v 对 y 的偏导数
    J(1, 2) = -fy * y * z_inv_2;  // v 对 z 的偏导数
}

// 从图像中提取补丁
// 输入：img (输入图像), pc (中心像素坐标), patch_tmp (输出补丁), level (金字塔层级)
void VIOManager::getImagePatch(cv::Mat img, V2D pc, float *patch_tmp, int level)
{
    const float u_ref = pc[0];  // 补丁中心的 u 坐标
    const float v_ref = pc[1];  // 补丁中心的 v 坐标
    const int scale = (1 << level);  // 当前层级的缩放因子
    const int u_ref_i = floorf(pc[0] / scale) * scale;  // 调整后的整数 u 坐标
    const int v_ref_i = floorf(pc[1] / scale) * scale;  // 调整后的整数 v 坐标
    const float subpix_u_ref = (u_ref - u_ref_i) / scale;  // u 方向的亚像素偏移
    const float subpix_v_ref = (v_ref - v_ref_i) / scale;  // v 方向的亚像素偏移
    const float w_ref_tl = (1.0 - subpix_u_ref) * (1.0 - subpix_v_ref);  // 左上角权重
    const float w_ref_tr = subpix_u_ref * (1.0 - subpix_v_ref);  // 右上角权重
    const float w_ref_bl = (1.0 - subpix_u_ref) * subpix_v_ref;  // 左下角权重
    const float w_ref_br = subpix_u_ref * subpix_v_ref;  // 右下角权重
    
    // 遍历补丁中的每个像素
    for (int x = 0; x < patch_size; x++)
    {
        uint8_t *img_ptr = (uint8_t *)img.data + (v_ref_i - patch_size_half * scale + x * scale) * width + (u_ref_i - patch_size_half * scale);  // 指向图像数据的指针
        for (int y = 0; y < patch_size; y++, img_ptr += scale)
        {
            // 使用双线性插值计算补丁像素值
            patch_tmp[patch_size_total * level + x * patch_size + y] =
                w_ref_tl * img_ptr[0] + w_ref_tr * img_ptr[scale] + w_ref_bl * img_ptr[scale * width] + w_ref_br * img_ptr[scale * width + scale];
        }
    }
}

// 将视觉点插入体素地图
// 输入：pt_new (新的视觉点)
void VIOManager::insertPointIntoVoxelMap(VisualPoint *pt_new)
{
    V3D pt_w(pt_new->pos_[0], pt_new->pos_[1], pt_new->pos_[2]);  // 提取点的世界坐标
    double voxel_size = 0.5;  // 体素大小设置为 0.5 米
    float loc_xyz[3];  // 体素坐标
    for (int j = 0; j < 3; j++)
    {
        loc_xyz[j] = pt_w[j] / voxel_size;  // 将世界坐标转换为体素坐标
        if (loc_xyz[j] < 0) { loc_xyz[j] -= 1.0; }  // 处理负坐标的情况
    }
    VOXEL_LOCATION position((int64_t)loc_xyz[0], (int64_t)loc_xyz[1], (int64_t)loc_xyz[2]);  // 创建体素位置对象
    
    auto iter = feat_map.find(position);  // 在特征地图中查找该体素
    if (iter != feat_map.end())
    {
        iter->second->voxel_points.push_back(pt_new);  // 如果存在，添加点到体素中
        iter->second->count++;  // 增加体素中的点计数
    }
    else
    {
        VOXEL_POINTS *ot = new VOXEL_POINTS(0);  // 如果不存在，创建新的体素点对象
        ot->voxel_points.push_back(pt_new);  // 添加点到新体素中
        feat_map[position] = ot;  // 将新体素插入特征地图
    }
}

// 计算仿射变换矩阵（基于单应性）
// 输入：cam (相机模型), px_ref (参考像素坐标), xyz_ref (参考 3D 坐标), normal_ref (参考法向量), T_cur_ref (当前到参考的变换), level_ref (参考层级), A_cur_ref (输出仿射矩阵)
void VIOManager::getWarpMatrixAffineHomography(const vk::AbstractCamera &cam, const V2D &px_ref, const V3D &xyz_ref, const V3D &normal_ref,
                                               const SE3 &T_cur_ref, const int level_ref, Matrix2d &A_cur_ref)
{
    // 创建单应矩阵
    const V3D t = T_cur_ref.inverse().translation();  // 获取逆变换的平移部分
    const Eigen::Matrix3d H_cur_ref =
        T_cur_ref.rotation_matrix() * (normal_ref.dot(xyz_ref) * Eigen::Matrix3d::Identity() - t * normal_ref.transpose());  // 计算单应矩阵
    
    // 使用单应投影计算仿射变换矩阵 A_ref_cur
    const int kHalfPatchSize = 4;  // 补丁半宽
    V3D f_du_ref(cam.cam2world(px_ref + Eigen::Vector2d(kHalfPatchSize, 0) * (1 << level_ref)));  // 计算 u 方向偏移的单位向量
    V3D f_dv_ref(cam.cam2world(px_ref + Eigen::Vector2d(0, kHalfPatchSize) * (1 << level_ref)));  // 计算 v 方向偏移的单位向量
    const V3D f_cur(H_cur_ref * xyz_ref);  // 将参考点投影到当前帧
    const V3D f_du_cur = H_cur_ref * f_du_ref;  // 将 u 偏移点投影到当前帧
    const V3D f_dv_cur = H_cur_ref * f_dv_ref;  // 将 v 偏移点投影到当前帧
    V2D px_cur(cam.world2cam(f_cur));  // 将当前点投影到像素坐标
    V2D px_du_cur(cam.world2cam(f_du_cur));  // 将 u 偏移点投影到像素坐标
    V2D px_dv_cur(cam.world2cam(f_dv_cur));  // 将 v 偏移点投影到像素坐标
    A_cur_ref.col(0) = (px_du_cur - px_cur) / kHalfPatchSize;  // 计算仿射矩阵的 u 方向列
    A_cur_ref.col(1) = (px_dv_cur - px_cur) / kHalfPatchSize;  // 计算仿射矩阵的 v 方向列
}

// 计算仿射变换矩阵（基于深度）
// 输入：cam (相机模型), px_ref (参考像素坐标), f_ref (参考单位向量), depth_ref (参考深度), T_cur_ref (当前到参考的变换), level_ref (参考层级), pyramid_level (金字塔层级), halfpatch_size (补丁半宽), A_cur_ref (输出仿射矩阵)
void VIOManager::getWarpMatrixAffine(const vk::AbstractCamera &cam, const Vector2d &px_ref, const Vector3d &f_ref, const double depth_ref,
                                     const SE3 &T_cur_ref, const int level_ref, const int pyramid_level, const int halfpatch_size,
                                     Matrix2d &A_cur_ref)
{
    // 计算仿射变换矩阵 A_ref_cur
    const Vector3d xyz_ref(f_ref * depth_ref);  // 根据深度计算参考点的 3D 坐标
    Vector3d xyz_du_ref(cam.cam2world(px_ref + Vector2d(halfpatch_size, 0) * (1 << level_ref) * (1 << pyramid_level)));  // 计算 u 方向偏移点的 3D 坐标
    Vector3d xyz_dv_ref(cam.cam2world(px_ref + Vector2d(0, halfpatch_size) * (1 << level_ref) * (1 << pyramid_level)));  // 计算 v 方向偏移点的 3D 坐标
    xyz_du_ref *= xyz_ref[2] / xyz_du_ref[2];  // 按深度比例调整 u 偏移点
    xyz_dv_ref *= xyz_ref[2] / xyz_dv_ref[2];  // 按深度比例调整 v 偏移点
    const Vector2d px_cur(cam.world2cam(T_cur_ref * (xyz_ref)));  // 将参考点投影到当前帧像素坐标
    const Vector2d px_du(cam.world2cam(T_cur_ref * (xyz_du_ref)));  // 将 u 偏移点投影到当前帧像素坐标
    const Vector2d px_dv(cam.world2cam(T_cur_ref * (xyz_dv_ref)));  // 将 v 偏移点投影到当前帧像素坐标
    A_cur_ref.col(0) = (px_du - px_cur) / halfpatch_size;  // 计算仿射矩阵的 u 方向列
    A_cur_ref.col(1) = (px_dv - px_cur) / halfpatch_size;  // 计算仿射矩阵的 v 方向列
}

// 执行仿射变换
// 输入：A_cur_ref (仿射矩阵), img_ref (参考图像), px_ref (参考像素坐标), level_ref (参考层级), search_level (搜索层级), pyramid_level (金字塔层级), halfpatch_size (补丁半宽), patch (输出补丁)
void VIOManager::warpAffine(const Matrix2d &A_cur_ref, const cv::Mat &img_ref, const Vector2d &px_ref, const int level_ref, const int search_level,
                            const int pyramid_level, const int halfpatch_size, float *patch)
{
    const int patch_size = halfpatch_size * 2;  // 计算补丁总大小
    const Matrix2f A_ref_cur = A_cur_ref.inverse().cast<float>();  // 计算逆仿射矩阵并转换为浮点型
    if (isnan(A_ref_cur(0, 0)))
    {
        printf("仿射变换矩阵包含 NaN，可能相机没有平移\n");  // 检查矩阵是否有效
        return;
    }

    float *patch_ptr = patch;  // 指向输出补丁的指针
    for (int y = 0; y < patch_size; ++y)
    {
        for (int x = 0; x < patch_size; ++x)  // 遍历补丁中的每个像素
        {
            Vector2f px_patch(x - halfpatch_size, y - halfpatch_size);  // 计算补丁内的相对像素坐标
            px_patch *= (1 << search_level);  // 根据搜索层级缩放
            px_patch *= (1 << pyramid_level);  // 根据金字塔层级缩放
            const Vector2f px(A_ref_cur * px_patch + px_ref.cast<float>());  // 变换到参考图像坐标
            if (px[0] < 0 || px[1] < 0 || px[0] >= img_ref.cols - 1 || px[1] >= img_ref.rows - 1)
                patch_ptr[patch_size_total * pyramid_level + y * patch_size + x] = 0;  // 如果超出图像范围，设为 0
            else
                patch_ptr[patch_size_total * pyramid_level + y * patch_size + x] = (float)vk::interpolateMat_8u(img_ref, px[0], px[1]);  // 插值计算像素值
        }
    }
}

// 获取最佳搜索层级
// 输入：A_cur_ref (仿射矩阵), max_level (最大层级)
// 输出：最佳搜索层级
int VIOManager::getBestSearchLevel(const Matrix2d &A_cur_ref, const int max_level)
{
    int search_level = 0;  // 初始化搜索层级
    double D = A_cur_ref.determinant();  // 计算仿射矩阵的行列式
    while (D > 3.0 && search_level < max_level)  // 如果行列式大于 3 且未达到最大层级
    {
        search_level += 1;  // 增加搜索层级
        D *= 0.25;  // 行列式按 1/4 缩减
    }
    return search_level;  // 返回最佳搜索层级
}

// 计算归一化互相关（NCC）
// 输入：ref_patch (参考补丁), cur_patch (当前补丁), patch_size (补丁大小)
// 输出：NCC 值
double VIOManager::calculateNCC(float *ref_patch, float *cur_patch, int patch_size)
{
    double sum_ref = std::accumulate(ref_patch, ref_patch + patch_size, 0.0);  // 计算参考补丁像素值总和
    double mean_ref = sum_ref / patch_size;  // 计算参考补丁均值

    double sum_cur = std::accumulate(cur_patch, cur_patch + patch_size, 0.0);  // 计算当前补丁像素值总和
    double mean_curr = sum_cur / patch_size;  // 计算当前补丁均值

    double numerator = 0, demoniator1 = 0, demoniator2 = 0;  // 初始化分子和分母
    for (int i = 0; i < patch_size; i++)
    {
        double n = (ref_patch[i] - mean_ref) * (cur_patch[i] - mean_curr);  // 计算协方差项
        numerator += n;  // 累加分子
        demoniator1 += (ref_patch[i] - mean_ref) * (ref_patch[i] - mean_ref);  // 计算参考补丁方差
        demoniator2 += (cur_patch[i] - mean_curr) * (cur_patch[i] - mean_curr);  // 计算当前补丁方差
    }
    return numerator / sqrt(demoniator1 * demoniator2 + 1e-10);  // 返回 NCC 值，添加小值避免除零
}

// 从视觉稀疏地图中检索点
// 输入：img (当前图像), pg (点云数据), plane_map (平面地图)
void VIOManager::retrieveFromVisualSparseMap(cv::Mat img, vector<pointWithVar> &pg, const unordered_map<VOXEL_LOCATION, VoxelOctoTree *> &plane_map)
{
    if (feat_map.size() <= 0) return;  // 如果特征地图为空，直接返回
    double ts0 = omp_get_wtime();  // 记录开始时间

    visual_submap->reset();  // 重置视觉子地图
    sub_feat_map.clear();  // 清空子特征地图

    float voxel_size = 0.5;  // 体素大小

    if (!normal_en) warp_map.clear();  // 如果未启用法向量，清空变换地图

    cv::Mat depth_img = cv::Mat::zeros(height, width, CV_32FC1);  // 创建深度图像
    float *it = (float *)depth_img.data;  // 指向深度图像数据的指针

    int loc_xyz[3];  // 体素坐标

    // 遍历点云数据，生成深度图和子特征地图
    for (int i = 0; i < pg.size(); i++)
    {
        V3D pt_w = pg[i].point_w;  // 获取世界坐标点
        for (int j = 0; j < 3; j++)
        {
            loc_xyz[j] = floor(pt_w[j] / voxel_size);  // 计算体素坐标
            if (loc_xyz[j] < 0) { loc_xyz[j] -= 1.0; }  // 处理负坐标
        }
        VOXEL_LOCATION position(loc_xyz[0], loc_xyz[1], loc_xyz[2]);  // 创建体素位置

        auto iter = sub_feat_map.find(position);  // 在子特征地图中查找
        if (iter == sub_feat_map.end()) { sub_feat_map[position] = 0; }  // 如果不存在，插入新位置
        else { iter->second = 0; }  // 如果存在，标记为 0

        V3D pt_c(new_frame_->w2f(pt_w));  // 将世界坐标转换为当前帧坐标
        if (pt_c[2] > 0)  // 如果深度大于 0
        {
            V2D px = new_frame_->cam_->world2cam(pt_c);  // 投影到像素坐标
            if (new_frame_->cam_->isInFrame(px.cast<int>(), border))  // 如果在图像范围内
            {
                float depth = pt_c[2];  // 获取深度值
                int col = int(px[0]);  // 像素列
                int row = int(px[1]);  // 像素行
                it[width * row + col] = depth;  // 记录深度值到深度图
            }
        }
    }

    vector<VOXEL_LOCATION> DeleteKeyList;  // 需要删除的体素位置列表

    // 遍历子特征地图，检查体素点是否在视野内
    for (auto &iter : sub_feat_map)
    {
        VOXEL_LOCATION position = iter.first;  // 获取体素位置
        auto corre_voxel = feat_map.find(position);  // 在特征地图中查找对应体素

        if (corre_voxel != feat_map.end())
        {
            bool voxel_in_fov = false;  // 标记体素是否在视野内
            std::vector<VisualPoint *> &voxel_points = corre_voxel->second->voxel_points;  // 获取体素中的点
            int voxel_num = voxel_points.size();  // 体素中的点数

            for (int i = 0; i < voxel_num; i++)
            {
                VisualPoint *pt = voxel_points[i];  // 获取当前点
                if (pt == nullptr || pt->obs_.size() == 0) continue;  // 如果点为空或无观测，跳过

                V3D norm_vec(new_frame_->T_f_w_.rotation_matrix() * pt->normal_);  // 计算当前帧下的法向量
                V3D dir(new_frame_->T_f_w_ * pt->pos_);  // 计算当前帧下的方向向量
                if (dir[2] < 0) continue;  // 如果深度为负，跳过

                V2D pc(new_frame_->w2c(pt->pos_));  // 投影到当前帧像素坐标
                if (new_frame_->cam_->isInFrame(pc.cast<int>(), border))  // 如果在图像范围内
                {
                    voxel_in_fov = true;  // 标记体素在视野内
                    int index = static_cast<int>(pc[1] / grid_size) * grid_n_width + static_cast<int>(pc[0] / grid_size);  // 计算网格索引
                    grid_num[index] = TYPE_MAP;  // 标记网格类型为地图点
                    Vector3d obs_vec(new_frame_->pos() - pt->pos_);  // 计算观测向量
                    float cur_dist = obs_vec.norm();  // 计算当前距离
                    if (cur_dist <= map_dist[index])  // 如果距离更近
                    {
                        map_dist[index] = cur_dist;  // 更新最小距离
                        retrieve_voxel_points[index] = pt;  // 更新检索点
                    }
                }
            }
            if (!voxel_in_fov) { DeleteKeyList.push_back(position); }  // 如果体素不在视野内，加入删除列表
        }
    }

    // 如果启用了光线投射模块
    if (raycast_en)
    {
        for (int i = 0; i < length; i++)  // 遍历所有网格
        {
            if (grid_num[i] == TYPE_MAP || border_flag[i] == 1) continue;  // 如果是地图点或边界网格，跳过

            for (const auto &it : rays_with_sample_points[i])  // 遍历当前网格的光线采样点
            {
                V3D sample_point_w = new_frame_->f2w(it);  // 将采样点转换到世界坐标
                for (int j = 0; j < 3; j++)
                {
                    loc_xyz[j] = floor(sample_point_w[j] / voxel_size);  // 计算体素坐标
                    if (loc_xyz[j] < 0) { loc_xyz[j] -= 1.0; }  // 处理负坐标
                }

                VOXEL_LOCATION sample_pos(loc_xyz[0], loc_xyz[1], loc_xyz[2]);  // 创建采样点位置

                auto corre_sub_feat_map = sub_feat_map.find(sample_pos);  // 在子特征地图中查找
                if (corre_sub_feat_map != sub_feat_map.end()) break;  // 如果存在，跳出循环

                auto corre_feat_map = feat_map.find(sample_pos);  // 在特征地图中查找
                if (corre_feat_map != feat_map.end())
                {
                    bool voxel_in_fov = false;  // 标记体素是否在视野内
                    std::vector<VisualPoint *> &voxel_points = corre_feat_map->second->voxel_points;  // 获取体素中的点
                    int voxel_num = voxel_points.size();  // 体素中的点数
                    if (voxel_num == 0) continue;  // 如果体素为空，跳过

                    for (int j = 0; j < voxel_num; j++)
                    {
                        VisualPoint *pt = voxel_points[j];  // 获取当前点
                        if (pt == nullptr || pt->obs_.size() == 0) continue;  // 如果点为空或无观测，跳过

                        V3D norm_vec(new_frame_->T_f_w_.rotation_matrix() * pt->normal_);  // 计算当前帧下的法向量
                        V3D dir(new_frame_->T_f_w_ * pt->pos_);  // 计算当前帧下的方向向量
                        if (dir[2] < 0) continue;  // 如果深度为负，跳过
                        dir.normalize();  // 归一化方向向量

                        V2D pc(new_frame_->w2c(pt->pos_));  // 投影到当前帧像素坐标
                        if (new_frame_->cam_->isInFrame(pc.cast<int>(), border))  // 如果在图像范围内
                        {
                            voxel_in_fov = true;  // 标记体素在视野内
                            int index = static_cast<int>(pc[1] / grid_size) * grid_n_width + static_cast<int>(pc[0] / grid_size);  // 计算网格索引
                            grid_num[index] = TYPE_MAP;  // 标记网格类型为地图点
                            Vector3d obs_vec(new_frame_->pos() - pt->pos_);  // 计算观测向量
                            float cur_dist = obs_vec.norm();  // 计算当前距离
                            if (cur_dist <= map_dist[index])  // 如果距离更近
                            {
                                map_dist[index] = cur_dist;  // 更新最小距离
                                retrieve_voxel_points[index] = pt;  // 更新检索点
                            }
                        }
                    }
                    if (voxel_in_fov) sub_feat_map[sample_pos] = 0;  // 如果体素在视野内，添加到子特征地图
                    break;  // 跳出采样点循环
                }
                else
                {
                    auto iter = plane_map.find(sample_pos);  // 在平面地图中查找
                    if (iter != plane_map.end())
                    {
                        VoxelOctoTree *current_octo = iter->second->find_correspond(sample_point_w);  // 查找对应的八叉树节点
                        if (current_octo->plane_ptr_->is_plane_)  // 如果是平面
                        {
                            pointWithVar plane_center;  // 平面中心点
                            VoxelPlane &plane = *current_octo->plane_ptr_;  // 获取平面数据
                            plane_center.point_w = plane.center_;  // 设置平面中心
                            plane_center.normal = plane.normal_;  // 设置平面法向量
                            visual_submap->add_from_voxel_map.push_back(plane_center);  // 添加到视觉子地图
                            break;  // 跳出采样点循环
                        }
                    }
                }
            }
        }
    }

    // 删除不在视野内的体素
    for (auto &key : DeleteKeyList)
    {
        sub_feat_map.erase(key);  // 从子特征地图中移除
    }

    // 处理检索到的点，计算变换和误差
    for (int i = 0; i < length; i++)
    {
        if (grid_num[i] == TYPE_MAP)  // 如果是地图点
        {
            VisualPoint *pt = retrieve_voxel_points[i];  // 获取检索点
            V2D pc(new_frame_->w2c(pt->pos_));  // 投影到当前帧像素坐标

            V3D pt_cam(new_frame_->w2f(pt->pos_));  // 转换到当前帧坐标
            bool depth_continous = false;  // 标记深度是否连续
            for (int u = -patch_size_half; u <= patch_size_half; u++)  // 检查补丁内的深度连续性
            {
                for (int v = -patch_size_half; v <= patch_size_half; v++)
                {
                    if (u == 0 && v == 0) continue;  // 跳过中心点
                    float depth = it[width * (v + int(pc[1])) + u + int(pc[0])];  // 获取深度值
                    if (depth == 0.) continue;  // 如果深度为 0，跳过
                    double delta_dist = abs(pt_cam[2] - depth);  // 计算深度差
                    if (delta_dist > 0.5)  // 如果深度差大于 0.5
                    {
                        depth_continous = true;  // 标记为不连续
                        break;
                    }
                }
                if (depth_continous) break;  // 如果不连续，跳出循环
            }
            if (depth_continous) continue;  // 如果深度不连续，跳过该点

            Feature *ref_ftr;  // 参考特征
            std::vector<float> patch_wrap(warp_len);  // 变换补丁
            int search_level;  // 搜索层级
            Matrix2d A_cur_ref_zero;  // 仿射变换矩阵

            if (!pt->is_normal_initialized_) continue;  // 如果法向量未初始化，跳过

            // 如果启用了法向量处理
            if (normal_en)
            {
                float phtometric_errors_min = std::numeric_limits<float>::max();  // 初始化最小光度误差

                if (pt->obs_.size() == 1)  // 如果只有一次观测
                {
                    ref_ftr = *pt->obs_.begin();  // 获取唯一的参考特征
                    pt->ref_patch = ref_ftr;  // 设置参考补丁
                    pt->has_ref_patch_ = true;  // 标记已设置参考补丁
                }
                else if (!pt->has_ref_patch_)  // 如果尚未设置参考补丁
                {
                    for (auto it = pt->obs_.begin(), ite = pt->obs_.end(); it != ite; ++it)  // 遍历所有观测
                    {
                        Feature *ref_patch_temp = *it;  // 当前参考补丁
                        float *patch_temp = ref_patch_temp->patch_;  // 补丁数据
                        float phtometric_errors = 0.0;  // 光度误差
                        int count = 0;  // 计数器
                        for (auto itm = pt->obs_.begin(), itme = pt->obs_.end(); itm != itme; ++itm)  // 比较所有其他观测
                        {
                            if ((*itm)->id_ == ref_patch_temp->id_) continue;  // 跳过自身
                            float *patch_cache = (*itm)->patch_;  // 其他补丁数据
                            for (int ind = 0; ind < patch_size_total; ind++)  // 计算光度误差
                            {
                                phtometric_errors += (patch_temp[ind] - patch_cache[ind]) * (patch_temp[ind] - patch_cache[ind]);
                            }
                            count++;
                        }
                        phtometric_errors = phtometric_errors / count;  // 计算平均光度误差
                        if (phtometric_errors < phtometric_errors_min)  // 如果误差更小
                        {
                            phtometric_errors_min = phtometric_errors;  // 更新最小误差
                            ref_ftr = ref_patch_temp;  // 更新参考特征
                        }
                    }
                    pt->ref_patch = ref_ftr;  // 设置参考补丁
                    pt->has_ref_patch_ = true;  // 标记已设置参考补丁
                }
                else { ref_ftr = pt->ref_patch; }  // 如果已设置，直接使用参考补丁
            }
            else  // 如果未启用法向量处理
            {
                if (!pt->getCloseViewObs(new_frame_->pos(), ref_ftr, pc)) continue;  // 获取最近的观测特征，若失败则跳过
            }

            // 计算仿射变换
            if (normal_en)
            {
                V3D norm_vec = (ref_ftr->T_f_w_.rotation_matrix() * pt->normal_).normalized();  // 计算归一化的法向量
                V3D pf(ref_ftr->T_f_w_ * pt->pos_);  // 将点转换到参考帧坐标
                SE3 T_cur_ref = new_frame_->T_f_w_ * ref_ftr->T_f_w_.inverse();  // 计算当前到参考的变换

                getWarpMatrixAffineHomography(*cam, ref_ftr->px_, pf, norm_vec, T_cur_ref, 0, A_cur_ref_zero);  // 获取仿射变换矩阵
                search_level = getBestSearchLevel(A_cur_ref_zero, 2);  // 获取最佳搜索层级
            }
            else
            {
                auto iter_warp = warp_map.find(ref_ftr->id_);  // 在变换地图中查找
                if (iter_warp != warp_map.end())
                {
                    search_level = iter_warp->second->search_level;  // 使用缓存的搜索层级
                    A_cur_ref_zero = iter_warp->second->A_cur_ref;  // 使用缓存的仿射矩阵
                }
                else
                {
                    getWarpMatrixAffine(*cam, ref_ftr->px_, ref_ftr->f_, (ref_ftr->pos() - pt->pos_).norm(), new_frame_->T_f_w_ * ref_ftr->T_f_w_.inverse(),
                                        ref_ftr->level_, 0, patch_size_half, A_cur_ref_zero);  // 计算新的仿射变换矩阵
                    search_level = getBestSearchLevel(A_cur_ref_zero, 2);  // 获取最佳搜索层级
                    Warp *ot = new Warp(search_level, A_cur_ref_zero);  // 创建新的变换对象
                    warp_map[ref_ftr->id_] = ot;  // 插入变换地图
                }
            }

            // 生成变换补丁
            for (int pyramid_level = 0; pyramid_level <= patch_pyrimid_level - 1; pyramid_level++)
            {
                warpAffine(A_cur_ref_zero, ref_ftr->img_, ref_ftr->px_, ref_ftr->level_, search_level, pyramid_level, patch_size_half, patch_wrap.data());
            }

            getImagePatch(img, pc, patch_buffer.data(), 0);  // 从当前图像提取补丁

            float error = 0.0;  // 初始化误差
            for (int ind = 0; ind < patch_size_total; ind++)  // 计算光度误差
            {
                error += (ref_ftr->inv_expo_time_ * patch_wrap[ind] - state->inv_expo_time * patch_buffer[ind]) *
                         (ref_ftr->inv_expo_time_ * patch_wrap[ind] - state->inv_expo_time * patch_buffer[ind]);
            }

            if (ncc_en)  // 如果启用了 NCC 检查
            {
                double ncc = calculateNCC(patch_wrap.data(), patch_buffer.data(), patch_size_total);  // 计算 NCC
                if (ncc < ncc_thre) continue;  // 如果 NCC 小于阈值，跳过
            }

            if (error > outlier_threshold * patch_size_total) continue;  // 如果误差大于阈值，跳过

            // 将点添加到视觉子地图
            visual_submap->voxel_points.push_back(pt);  // 添加体素点
            visual_submap->propa_errors.push_back(error); Swapped(error);  // 添加传播误差
            visual_submap->search_levels.push_back(search_level);  // 添加搜索层级
            visual_submap->errors.push_back(error);  // 添加误差
            visual_submap->warp_patch.push_back(patch_wrap);  // 添加变换补丁
            visual_submap->inv_expo_list.push_back(ref_ftr->inv_expo_time_);  // 添加逆曝光时间
        }
    }
    total_points = visual_submap->voxel_points.size();  // 更新总点数

    printf("[ VIO ] 从视觉稀疏地图检索到 %d 个点\n", total_points);  // 记录检索到的点数
}

// 计算雅可比矩阵并更新扩展卡尔曼滤波器
// 输入：img (当前图像)
void VIOManager::computeJacobianAndUpdateEKF(cv::Mat img)
{
    if (total_points == 0) return;  // 如果没有点，直接返回

    compute_jacobian_time = update_ekf_time = 0.0;  // 初始化时间记录

    // 从高到低遍历金字塔层级
    for (int level = patch_pyrimid_level - 1; level >= 0; level--)
    {
        if (inverse_composition_en)  // 如果启用了逆向组合
        {
            has_ref_patch_cache = false;  // 重置参考补丁缓存标志
            updateStateInverse(img, level);  // 使用逆向方法更新状态
        }
        else
            updateState(img, level);  // 使用正向方法更新状态
    }
    state->cov -= G * state->cov;  // 更新状态协方差
    updateFrameState(*state);  // 更新帧状态
}

// 生成视觉地图点
// 输入：img (当前图像), pg (点云数据)
void VIOManager::generateVisualMapPoints(cv::Mat img, vector<pointWithVar> &pg)
{
    if (pg.size() <= 10) return;  // 如果点云数量少于 10，直接返回

    // double t0 = omp_get_wtime();  // 记录开始时间

    // 遍历点云数据，生成新的视觉地图点
    for (int i = 0; i < pg.size(); i++)
    {
        if (pg[i].normal == V3D(0, 0, 0)) continue;  // 如果法向量为零，跳过

        V3D pt = pg[i].point_w;  // 获取世界坐标点
        V2D pc(new_frame_->w2c(pt));  // 投影到当前帧像素坐标

        if (new_frame_->cam_->isInFrame(pc.cast<int>(), border))  // 如果在图像范围内
        {
            int index = static_cast<int>(pc[1] / grid_size) * grid_n_width + static_cast<int>(pc[0] / grid_size);  // 计算网格索引

            if (grid_num[index] != TYPE_MAP)  // 如果不是地图点
            {
                float cur_value = vk::shiTomasiScore(img, pc[0], pc[1]);  // 计算 Shi-Tomasi 角点分数
                if (cur_value > scan_value[index])  // 如果分数更高
                {
                    scan_value[index] = cur_value;  // 更新扫描值
                    append_voxel_points[index] = pg[i];  // 更新追加点
                    grid_num[index] = TYPE_POINTCLOUD;  // 标记为点云类型
                }
            }
        }
    }

    // 处理从体素地图添加的点
    for (int j = 0; j < visual_submap->add_from_voxel_map.size(); j++)
    {
        V3D pt = visual_submap->add_from_voxel_map[j].point_w;  // 获取世界坐标点
        V2D pc(new_frame_->w2c(pt));  // 投影到当前帧像素坐标

        if (new_frame_->cam_->isInFrame(pc.cast<int>(), border))  // 如果在图像范围内
        {
            int index = static_cast<int>(pc[1] / grid_size) * grid_n_width + static_cast<int>(pc[0] / grid_size);  // 计算网格索引

            if (grid_num[index] != TYPE_MAP)  // 如果不是地图点
            {
                float cur_value = vk::shiTomasiScore(img, pc[0], pc[1]);  // 计算 Shi-Tomasi 角点分数
                if (cur_value > scan_value[index])  // 如果分数更高
                {
                    scan_value[index] = cur_value;  // 更新扫描值
                    append_voxel_points[index] = visual_submap->add_from_voxel_map[j];  // 更新追加点
                    grid_num[index] = TYPE_POINTCLOUD;  // 标记为点云类型
                }
            }
        }
    }

    int add = 0;  // 新增点计数器
    for (int i = 0; i < length; i++)  // 遍历所有网格
    {
        if (grid_num[i] == TYPE_POINTCLOUD)  // 如果是点云点
        {
            pointWithVar pt_var = append_voxel_points[i];  // 获取追加点
            V3D pt = pt_var.point_w;  // 获取世界坐标

            V3D norm_vec(new_frame_->T_f_w_.rotation_matrix() * pt_var.normal);  // 计算当前帧下的法向量
            V3D dir(new_frame_->T_f_w_ * pt);  // 计算当前帧下的方向向量
            dir.normalize();  // 归一化方向向量
            double cos_theta = dir.dot(norm_vec);  // 计算视角余弦值
            V2D pc(new_frame_->w2c(pt));  // 投影到当前帧像素坐标

            float *patch = new float[patch_size_total];  // 分配补丁内存
            getImagePatch(img, pc, patch, 0);  // 从当前图像提取补丁

            VisualPoint *pt_new = new VisualPoint(pt);  // 创建新的视觉点

            Vector3d f = cam->cam2world(pc);  // 将像素坐标转换为相机坐标系单位向量
            Feature *ftr_new = new Feature(pt_new, patch, pc, f, new_frame_->T_f_w_, 0);  // 创建新的特征
            ftr_new->img_ = img;  // 设置特征图像
            ftr_new->id_ = new_frame_->id_;  // 设置特征 ID
            ftr_new->inv_expo_time_ = state->inv_expo_time;  // 设置逆曝光时间

            pt_new->addFrameRef(ftr_new);  // 添加特征引用
            pt_new->covariance_ = pt_var.var;  // 设置协方差
            pt_new->is_normal_initialized_ = true;  // 标记法向量已初始化

            if (cos_theta < 0) { pt_new->normal_ = -pt_var.normal; }  // 如果视角反向，调整法向量
            else { pt_new->normal_ = pt_var.normal; }  // 否则使用原始法向量
            
            pt_new->previous_normal_ = pt_new->normal_;  // 设置前一法向量

            insertPointIntoVoxelMap(pt_new);  // 将点插入体素地图
            add += 1;  // 增加新增点计数
        }
    }

    printf("[ VIO ] 新增 %d 个视觉地图点\n", add);  // 记录新增点数
}

// 更新视觉地图点
// 输入：img (当前图像)
void VIOManager::updateVisualMapPoints(cv::Mat img)
{
    if (total_points == 0) return;  // 如果没有点，直接返回

    int update_num = 0;  // 更新点计数器
    SE3 pose_cur = new_frame_->T_f_w_;  // 获取当前帧位姿
    for (int i = 0; i < total_points; i++)  // 遍历所有视觉子地图点
    {
        VisualPoint *pt = visual_submap->voxel_points[i];  // 获取当前点
        if (pt == nullptr) continue;  // 如果点为空，跳过
        if (pt->is_converged_)  // 如果点已收敛
        { 
            pt->deleteNonRefPatchFeatures();  // 删除非参考补丁特征
            continue;  // 跳过
        }

        V2D pc(new_frame_->w2c(pt->pos_));  // 投影到当前帧像素坐标
        bool add_flag = false;  // 是否添加新特征的标志
        
        float *patch_temp = new float[patch_size_total];  // 分配临时补丁内存
        getImagePatch(img, pc, patch_temp, 0);  // 从当前图像提取补丁

        Feature *last_feature = pt->obs_.back();  // 获取最后一次观测特征

        // 检查是否需要添加新特征的条件
        SE3 pose_ref = last_feature->T_f_w_;  // 获取参考帧位姿
        SE3 delta_pose = pose_ref * pose_cur.inverse();  // 计算位姿差
        double delta_p = delta_pose.translation().norm();  // 计算平移距离
        double delta_theta = (delta_pose.rotation_matrix().trace() > 3.0 - 1e-6) ? 0.0 : std::acos(0.5 * (delta_pose.rotation_matrix().trace() - 1));  // 计算旋转角度
        if (delta_p > 0.5 || delta_theta > 0.3) add_flag = true;  // 如果平移大于 0.5 或旋转大于 0.3 度，标记为需要添加

        Vector2d last_px = last_feature->px_;  // 获取最后像素坐标
        double pixel_dist = (pc - last_px).norm();  // 计算像素距离
        if (pixel_dist > 40) add_flag = true;  // 如果像素距离大于 40，标记为需要添加

        // 控制观测特征的数量
        if (pt->obs_.size() >= 30)
        {
            Feature *ref_ftr;  // 参考特征
            pt->findMinScoreFeature(new_frame_->pos(), ref_ftr);  // 找到得分最小的特征
            pt->deleteFeatureRef(ref_ftr);  // 删除该特征
        }
        if (add_flag)  // 如果需要添加新特征
        {
            update_num += 1;  // 增加更新计数
            update_flag[i] = 1;  // 标记为已更新
            Vector3d f = cam->cam2world(pc);  // 将像素坐标转换为相机坐标系单位向量
            Feature *ftr_new = new Feature(pt, patch_temp, pc, f, new_frame_->T_f_w_, visual_submap->search_levels[i]);  // 创建新特征
            ftr_new->img_ = img;  // 设置特征图像
            ftr_new->id_ = new_frame_->id_;  // 设置特征 ID
            ftr_new->inv_expo_time_ = state->inv_expo_time;  // 设置逆曝光时间
            pt->addFrameRef(ftr_new);  // 添加特征引用
        }
    }
    printf("[ VIO ] 在视觉子地图中更新 %d 个点\n", update_num);  // 记录更新点数
}

// 函数：更新参考Patch
void VIOManager::updateReferencePatch(const unordered_map<VOXEL_LOCATION, VoxelOctoTree *> &plane_map)
{
  // 如果总点数为0，则直接返回，不执行后续操作
  // 原理：避免在没有数据的情况下进行无意义的计算
  if (total_points == 0) return;

  // 遍历视觉子地图中的所有体素点
  // 原理：逐点处理以更新每个点的参考Patch信息
  for (int i = 0; i < visual_submap->voxel_points.size(); i++)
  {
    // 获取当前体素点指针
    // 原理：通过指针操作直接访问和管理点数据，提高效率
    VisualPoint *pt = visual_submap->voxel_points[i];

    // 如果点的法向量未初始化，则跳过此点
    // 原理：未初始化的点无法参与后续计算，缺少必要信息
    if (!pt->is_normal_initialized_) continue;
    // 如果点已收敛，则跳过此点
    // 原理：已收敛的点无需再次更新，其状态已稳定
    if (pt->is_converged_) continue;
    // 如果点的观测次数少于5次，则跳过此点
    // 原理：观测次数太少的数据可靠性不足，不足以更新参考Patch
    if (pt->obs_.size() <= 5) continue;
    // 如果更新标志为0，则跳过此点
    // 原理：通过标志控制哪些点需要更新，避免不必要的计算
    if (update_flag[i] == 0) continue;

    // 获取点的世界坐标
    // 原理：后续计算需要基于世界坐标系中的位置进行空间分析
    const V3D &p_w = pt->pos_;
    // 定义局部坐标数组，用于存储转换后的网格坐标
    float loc_xyz[3];
    // 将世界坐标转换为体素网格坐标
    for (int j = 0; j < 3; j++)
    {
      // 将世界坐标除以体素分辨率（0.5）转换为网格坐标
      // 原理：将连续坐标离散化为整数网格坐标，便于体素索引
      loc_xyz[j] = p_w[j] / 0.5;
      // 如果坐标为负数，则向下取整（减1）
      // 原理：确保负数坐标正确映射到网格中，避免索引错误
      if (loc_xyz[j] < 0) { loc_xyz[j] -= 1.0; }
    }
    // 创建体素位置对象，使用整数坐标表示
    // 原理：使用整数坐标便于在哈希表中查找对应体素
    VOXEL_LOCATION position((int64_t)loc_xyz[0], (int64_t)loc_xyz[1], (int64_t)loc_xyz[2]);
    // 在平面映射中查找对应的体素
    // 原理：利用哈希表快速定位点所在的体素位置
    auto iter = plane_map.find(position);
    // 如果找到对应体素
    if (iter != plane_map.end())
    {
      // 定义当前八叉树体素指针
      VoxelOctoTree *current_octo;
      // 查找与点位置对应的八叉树节点
      // 原理：八叉树用于高效组织和管理空间数据，快速找到最近平面
      current_octo = iter->second->find_correspond(p_w);
      // 如果当前节点包含平面
      if (current_octo->plane_ptr_->is_plane_)
      {
        // 获取平面引用
        // 原理：直接操作平面数据以计算点与平面的关系
        VoxelPlane &plane = *current_octo->plane_ptr_;
        // 计算点到平面的有符号距离
        // 原理：使用平面方程 ax + by + cz + d = 0 计算点到平面的距离
        float dis_to_plane = plane.normal_(0) * p_w(0) + plane.normal_(1) * p_w(1) + plane.normal_(2) * p_w(2) + plane.d_;
        // 计算距离的绝对值
        // 原理：绝对值用于后续阈值判断，忽略方向
        float dis_to_plane_abs = fabs(dis_to_plane);
        // 计算点到平面中心的平方距离
        // 原理：欧几里得距离平方，用于判断点与平面中心的接近程度
        float dis_to_center = (plane.center_(0) - p_w(0)) * (plane.center_(0) - p_w(0)) +
                              (plane.center_(1) - p_w(1)) * (plane.center_(1) - p_w(1)) + 
                              (plane.center_(2) - p_w(2)) * (plane.center_(2) - p_w(2));
        // 计算点在平面上的投影距离
        // 原理：使用勾股定理，range_dis = sqrt(dis_to_center - dis_to_plane^2)
        float range_dis = sqrt(dis_to_center - dis_to_plane * dis_to_plane);
        // 如果投影距离小于平面的3倍半径
        // 原理：确保点在平面范围内，3倍半径作为经验阈值
        if (range_dis <= 3 * plane.radius_)
        {
          // 定义雅可比矩阵，用于不确定性计算
          // 原理：雅可比矩阵用于将平面参数的不确定性传播到距离上
          Eigen::Matrix<double, 1, 6> J_nq;
          // 前三列为点到中心的向量
          J_nq.block<1, 3>(0, 0) = p_w - plane.center_;
          // 后三列为负法向量
          J_nq.block<1, 3>(0, 3) = -plane.normal_;
          // 计算不确定性（方差）
          // 原理：sigma_l = J * Cov * J^T，计算点到平面的距离方差
          double sigma_l = J_nq * plane.plane_var_ * J_nq.transpose();
          // 叠加点的法向量协方差
          // 原理：考虑点本身的协方差对不确定性的贡献
          sigma_l += plane.normal_.transpose() * pt->covariance_ * plane.normal_;

          // 如果点到平面的距离小于3倍标准差
          // 原理：3倍标准差是统计上的显著性阈值，判断点是否属于平面
          if (dis_to_plane_abs < 3 * sqrt(sigma_l))
          {
            // 如果前一法向量与当前平面法向量方向相反，则取反
            // 原理：确保法向量方向一致，避免符号问题
            if (pt->previous_normal_.dot(plane.normal_) < 0) { pt->normal_ = -plane.normal_; }
            else { pt->normal_ = plane.normal_; }

            // 计算法向量更新的幅度
            // 原理：范数表示法向量变化的大小，用于判断收敛
            double normal_update = (pt->normal_ - pt->previous_normal_).norm();

            // 更新前一法向量为当前法向量
            // 原理：保存当前状态以便下次比较
            pt->previous_normal_ = pt->normal_;

            // 如果法向量更新幅度小于0.0001且观测次数大于10，则认为点已收敛
            // 原理：小更新幅度和高观测次数表明点状态稳定
            if (normal_update < 0.0001 && pt->obs_.size() > 10)
            {
              pt->is_converged_ = true;
              // 可选择将收敛点加入列表（已注释）
              // 原理：记录收敛点用于后续分析或可视化
              // visual_converged_point.push_back(pt);
            }
          }
        }
      }
    }

    // 初始化最大得分，用于选择最佳参考Patch
    float score_max = -1000.;
    // 遍历点的所有观测
    // 原理：通过比较所有观测的得分，选出最佳参考Patch
    for (auto it = pt->obs_.begin(), ite = pt->obs_.end(); it != ite; ++it)
    {
      // 获取当前观测的参考Patch
      Feature *ref_patch_temp = *it;
      // 获取Patch数据（像素值数组）
      float *patch_temp = ref_patch_temp->patch_;
      // 初始化NCC（归一化互相关）计算变量
      float NCC_up = 0.0;    // 分子
      float NCC_down1 = 0.0; // 分母第一部分
      float NCC_down2 = 0.0; // 分母第二部分
      float NCC = 0.0;       // NCC值
      float score = 0.0;     // 综合得分
      int count = 0;         // 计数器

      // 将点投影到参考帧坐标系
      // 原理：计算点在参考帧中的位置，用于后续方向性判断
      V3D pf = ref_patch_temp->T_f_w_ * pt->pos_;
      // 计算参考帧中的法向量
      // 原理：法向量用于判断点的可见性
      V3D norm_vec = ref_patch_temp->T_f_w_.rotation_matrix() * pt->normal_;
      // 归一化投影点
      // 原理：归一化便于计算方向余弦
      pf.normalize();
      // 计算投影点与法向量的夹角余弦
      // 原理：余弦值反映点与Patch的视角一致性
      double cos_angle = pf.dot(norm_vec);
      // 如果夹角大于20度（cos值小于0.86），则跳过（已注释）
      // 原理：限制视角差异以保证Patch匹配质量
      // if(fabs(cos_angle) < 0.86) continue;

      // 初始化参考Patch均值
      float ref_mean;
      // 如果均值未计算或接近0，则计算均值
      // 原理：避免重复计算均值，提高效率
      if (abs(ref_patch_temp->mean_) < 1e-6)
      {
        // 计算Patch像素值的和
        float ref_sum = std::accumulate(patch_temp, patch_temp + patch_size_total, 0.0);
        // 计算均值
        ref_mean = ref_sum / patch_size_total;
        // 保存均值到Patch对象
        ref_patch_temp->mean_ = ref_mean;
      }

      // 遍历其他观测，计算NCC
      // 原理：通过与其他观测比较，评估当前Patch的代表性
      for (auto itm = pt->obs_.begin(), itme = pt->obs_.end(); itm != itme; ++itm)
      {
        // 跳过与当前参考Patch相同的观测
        // 原理：避免自比较，减少冗余计算
        if ((*itm)->id_ == ref_patch_temp->id_) continue;
        // 获取其他Patch数据
        float *patch_cache = (*itm)->patch_;

        // 初始化其他Patch均值
        float other_mean;
        // 如果均值未计算或接近0，则计算均值
        if (abs((*itm)->mean_) < 1e-6)
        {
          // 计算Patch像素值的和
          float other_sum = std::accumulate(patch_cache, patch_cache + patch_size_total, 0.0);
          // 计算均值
          other_mean = other_sum / patch_size_total;
          // 保存均值到Patch对象
          (*itm)->mean_ = other_mean;
        }

        // 计算NCC的分子和分母
        // 原理：NCC = Σ((x - μx)(y - μy)) / sqrt(Σ(x - μx)^2 * Σ(y - μy)^2)
        for (int ind = 0; ind < patch_size_total; ind++)
        {
          // 分子：两Patch去均值后的像素乘积和
          NCC_up += (patch_temp[ind] - ref_mean) * (patch_cache[ind] - other_mean);
          // 分母：参考Patch去均值后的平方和
          NCC_down1 += (patch_temp[ind] - ref_mean) * (patch_temp[ind] - ref_mean);
          // 分母：其他Patch去均值后的平方和
          NCC_down2 += (patch_cache[ind] - other_mean) * (patch_cache[ind] - other_mean);
        }
        // 计算NCC并取绝对值，累加
        // 原理：绝对值确保正相关性，累加用于平均
        NCC += fabs(NCC_up / sqrt(NCC_down1 * NCC_down2));
        count++;
      }

      // 计算平均NCC
      // 原理：平均值反映Patch与其他观测的总体相似性
      NCC = NCC / count;

      // 综合得分：NCC加上夹角余弦
      // 原理：结合光度一致性（NCC）和几何一致性（cos_angle）
      score = NCC + cos_angle;

      // 保存得分到Patch对象
      ref_patch_temp->score_ = score;

      // 如果当前得分高于最大得分，则更新最佳参考Patch
      // 原理：选择得分最高的Patch作为参考，提升匹配质量
      if (score > score_max)
      {
        score_max = score;
        pt->ref_patch = ref_patch_temp;
        pt->has_ref_patch_ = true;
      }
    }
  }
}

// 函数：将参考Patch投影到当前帧
void VIOManager::projectPatchFromRefToCur(const unordered_map<VOXEL_LOCATION, VoxelOctoTree *> &plane_map)
{
  // 如果总点数为0，则直接返回
  // 原理：避免在无数据时执行投影操作
  if (total_points == 0) return;
  // 如果不是特定帧（已注释），可跳过
  // 原理：调试用，可限制特定帧执行
  // if(new_frame_->id_ != 2) return;

  // 定义Patch大小为25
  // 原理：固定Patch大小，便于一致性处理
  int patch_size = 25;
  // 定义输出目录，用于保存投影结果图像
  string dir = string(ROOT_DIR) + "Log/ref_cur_combine/";

  // 初始化结果图像（全黑），分别存储投影、带法向量的投影和密集投影
  cv::Mat result = cv::Mat::zeros(height, width, CV_8UC1);
  cv::Mat result_normal = cv::Mat::zeros(height, width, CV_8UC1);
  cv::Mat result_dense = cv::Mat::zeros(height, width, CV_8UC1);

  // 克隆当前帧图像，用于计算光度误差
  // 原理：保留原始图像以便后续比较
  cv::Mat img_photometric_error = new_frame_->img_.clone();

  // 获取图像数据指针，用于直接操作像素值
  uchar *it = (uchar *)result.data;
  uchar *it_normal = (uchar *)result_normal.data;
  uchar *it_dense = (uchar *)result_dense.data;

  // 定义像素结构体，包含位置和值
  // 原理：便于存储和操作投影后的像素数据
  struct pixel_member
  {
    Vector2f pixel_pos;   // 像素位置
    uint8_t pixel_value;  // 像素值
  };

  // 初始化计数器，用于命名输出文件
  int num = 0;
  // 遍历所有体素点
  for (int i = 0; i < visual_submap->voxel_points.size(); i++)
  {
    // 获取当前体素点
    VisualPoint *pt = visual_submap->voxel_points[i];

    // 如果法向量已初始化
    // 原理：只有法向量初始化的点才能进行投影
    if (pt->is_normal_initialized_)
    {
      // 定义参考特征指针
      Feature *ref_ftr;
      // 获取参考Patch
      ref_ftr = pt->ref_patch;
      // 将点投影到当前帧相机坐标系
      // 原理：通过相机模型将世界坐标转换为像素坐标
      V2D pc(new_frame_->w2c(pt->pos_));
      // 使用先验位姿投影，用于比较
      V2D pc_prior(new_frame_->w2c_prior(pt->pos_));

      // 计算参考帧中的法向量
      // 原理：法向量用于判断点的可见性方向
      V3D norm_vec(ref_ftr->T_f_w_.rotation_matrix() * pt->normal_);
      // 计算参考帧中的点位置
      V3D pf(ref_ftr->T_f_w_ * pt->pos_);

      // 如果点与法向量方向相反，则反转法向量
      // 原理：确保法向量指向相机方向，提升投影一致性
      if (pf.dot(norm_vec) < 0) norm_vec = -norm_vec;

      // 获取当前帧和参考帧图像
      cv::Mat img_cur = new_frame_->img_;
      cv::Mat img_ref = ref_ftr->img_;

      // 计算当前帧到参考帧的变换
      // 原理：通过位姿变换计算两帧之间的相对关系
      SE3 T_cur_ref = new_frame_->T_f_w_ * ref_ftr->T_f_w_.inverse();
      // 定义仿射变换矩阵
      Matrix2d A_cur_ref;
      // 计算仿射变换矩阵（基于单应性）
      // 原理：通过点、法向量和位姿计算图像间的仿射变换
      getWarpMatrixAffineHomography(*cam, ref_ftr->px_, pf, norm_vec, T_cur_ref, 0, A_cur_ref);

      // 计算最佳搜索层级
      // 原理：根据变换矩阵选择合适的金字塔层级，提升效率
      int search_level = getBestSearchLevel(A_cur_ref.inverse(), 2);

      // 计算变换行列式
      // 原理：行列式反映变换的缩放程度，过大表示畸变严重
      double D = A_cur_ref.determinant();
      // 如果行列式过大（变换过于剧烈），跳过
      if (D > 3) continue;

      // 计数器加1，用于唯一标识投影结果
      num++;

      // 拼接当前帧和参考帧图像，用于可视化
      cv::Mat ref_cur_combine_temp;
      int radius = 20; // 定义矩形框半径
      cv::hconcat(img_cur, img_ref, ref_cur_combine_temp);
      // 转换为彩色图像，便于绘制标记
      cv::cvtColor(ref_cur_combine_temp, ref_cur_combine_temp, CV_GRAY2BGR);

      // 获取当前帧Patch数据
      // 原理：提取当前帧对应位置的像素值，用于误差计算
      getImagePatch(img_cur, pc, patch_buffer.data(), 0);

      // 初始化误差变量
      float error_est = 0.0; // 估计误差
      float error_gt = 0.0;  // 真实误差（未使用）

      // 计算光度误差
      // 原理：比较参考Patch和当前Patch的光度一致性
      for (int ind = 0; ind < patch_size_total; ind++)
      {
        // 计算参考和当前帧Patch的光度差平方和，考虑曝光时间
        error_est += (ref_ftr->inv_expo_time_ * visual_submap->warp_patch[i][ind] - state->inv_expo_time * patch_buffer[ind]) *
                     (ref_ftr->inv_expo_time_ * visual_submap->warp_patch[i][ind] - state->inv_expo_time * patch_buffer[ind]);
      }
      // 定义调试信息字符串，显示曝光时间和误差
      std::string ref_est = "ref_est " + std::to_string(1.0 / ref_ftr->inv_expo_time_);
      std::string cur_est = "cur_est " + std::to_string(1.0 / state->inv_expo_time);
      std::string cur_propa = "cur_gt " + std::to_string(error_gt);
      std::string cur_optimize = "cur_est " + std::to_string(error_est);

      // 在图像上绘制调试信息
      cv::putText(ref_cur_combine_temp, ref_est, cv::Point2f(ref_ftr->px_[0] + img_cur.cols - 40, ref_ftr->px_[1] + 40), cv::FONT_HERSHEY_COMPLEX, 0.4,
                  cv::Scalar(0, 255, 0), 1, 8, 0);
      cv::putText(ref_cur_combine_temp, cur_est, cv::Point2f(pc[0] - 40, pc[1] + 40), cv::FONT_HERSHEY_COMPLEX, 0.4, cv::Scalar(0, 255, 0), 1, 8, 0);
      cv::putText(ref_cur_combine_temp, cur_propa, cv::Point2f(pc[0] - 40, pc[1] + 60), cv::FONT_HERSHEY_COMPLEX, 0.4, cv::Scalar(0, 0, 255), 1, 8, 0);
      cv::putText(ref_cur_combine_temp, cur_optimize, cv::Point2f(pc[0] - 40, pc[1] + 80), cv::FONT_HERSHEY_COMPLEX, 0.4, cv::Scalar(0, 255, 0), 1, 8, 0);

      // 绘制参考和当前帧的矩形框，用于可视化
      cv::rectangle(ref_cur_combine_temp, cv::Point2f(ref_ftr->px_[0] + img_cur.cols - radius, ref_ftr->px_[1] - radius),
                    cv::Point2f(ref_ftr->px_[0] + img_cur.cols + radius, ref_ftr->px_[1] + radius), cv::Scalar(0, 0, 255), 1);
      cv::rectangle(ref_cur_combine_temp, cv::Point2f(pc[0] - radius, pc[1] - radius), cv::Point2f(pc[0] + radius, pc[1] + radius),
                    cv::Scalar(0, 255, 0), 1);
      cv::rectangle(ref_cur_combine_temp, cv::Point2f(pc_prior[0] - radius, pc_prior[1] - radius),
                    cv::Point2f(pc_prior[0] + radius, pc_prior[1] + radius), cv::Scalar(255, 255, 255), 1);
      // 绘制中心点
      cv::circle(ref_cur_combine_temp, cv::Point2f(ref_ftr->px_[0] + img_cur.cols, ref_ftr->px_[1]), 1, cv::Scalar(0, 0, 255), -1, 8);
      cv::circle(ref_cur_combine_temp, cv::Point2f(pc[0], pc[1]), 1, cv::Scalar(0, 255, 0), -1, 8);
      cv::circle(ref_cur_combine_temp, cv::Point2f(pc_prior[0], pc_prior[1]), 1, cv::Scalar(255, 255, 255), -1, 8);
      // 保存投影结果图像
      cv::imwrite(dir + std::to_string(new_frame_->id_) + "_" + std::to_string(ref_ftr->id_) + "_" + std::to_string(num) + ".png",
                  ref_cur_combine_temp);

      // 定义像素变换矩阵，存储投影结果
      std::vector<std::vector<pixel_member>> pixel_warp_matrix;

      // 遍历Patch区域，计算投影像素
      for (int y = 0; y < patch_size; ++y)
      {
        vector<pixel_member> pixel_warp_vec;
        for (int x = 0; x < patch_size; ++x)
        {
          // 计算Patch内的像素偏移
          Vector2f px_patch(x - patch_size / 2, y - patch_size / 2);
          // 根据搜索层级缩放
          px_patch *= (1 << search_level);
          // 计算参考帧中的像素位置
          const Vector2f px_ref(px_patch + ref_ftr->px_.cast<float>());
          // 插值获取像素值
          // 原理：双线性插值获取参考帧中的像素值
          uint8_t pixel_value = (uint8_t)vk::interpolateMat_8u(img_ref, px_ref[0], px_ref[1]);

          // 通过仿射变换计算当前帧中的像素位置
          const Vector2f px(A_cur_ref.cast<float>() * px_patch + pc.cast<float>());
          // 如果超出图像边界，则跳过
          if (px[0] < 0 || px[1] < 0 || px[0] >= img_cur.cols - 1 || px[1] >= img_cur.rows - 1)
            continue;
          else
          {
            // 记录像素位置和值
            pixel_member pixel_warp;
            pixel_warp.pixel_pos << px[0], px[1];
            pixel_warp.pixel_value = pixel_value;
            pixel_warp_vec.push_back(pixel_warp);
          }
        }
        pixel_warp_matrix.push_back(pixel_warp_vec);
      }

      // 初始化边界值，用于确定投影区域
      float x_min = 1000;
      float y_min = 1000;
      float x_max = 0;
      float y_max = 0;

      // 计算变换区域的边界
      for (int i = 0; i < pixel_warp_matrix.size(); i++)
      {
        vector<pixel_member> pixel_warp_row = pixel_warp_matrix[i];
        for (int j = 0; j < pixel_warp_row.size(); j++)
        {
          float x_temp = pixel_warp_row[j].pixel_pos[0];
          float y_temp = pixel_warp_row[j].pixel_pos[1];
          if (x_temp < x_min) x_min = x_temp;
          if (y_temp < y_min) y_min = y_temp;
          if (x_temp > x_max) x_max = x_temp;
          if (y_temp > y_max) y_max = y_temp;
        }
      }
      // 取整计算边界
      int x_min_i = floor(x_min);
      int y_min_i = floor(y_min);
      int x_max_i = ceil(x_max);
      int y_max_i = ceil(y_max);
      // 计算逆变换矩阵
      // 原理：逆变换用于从当前帧映射回参考帧
      Matrix2f A_cur_ref_Inv = A_cur_ref.inverse().cast<float>();
      // 遍历边界内的像素，生成密集投影
      for (int i = x_min_i; i < x_max_i; i++)
      {
        for (int j = y_min_i; j < y_max_i; j++)
        {
          Eigen::Vector2f pc_temp(i, j);
          // 逆变换回参考帧
          Vector2f px_patch = A_cur_ref_Inv * (pc_temp - pc.cast<float>());
          // 检查是否在Patch范围内
          if (px_patch[0] > (-patch_size / 2 * (1 << search_level)) && px_patch[0] < (patch_size / 2 * (1 << search_level)) &&
              px_patch[1] > (-patch_size / 2 * (1 << search_level)) && px_patch[1] < (patch_size / 2 * (1 << search_level)))
          {
            // 计算参考帧像素位置并插值
            const Vector2f px_ref(px_patch + ref_ftr->px_.cast<float>());
            uint8_t pixel_value = (uint8_t)vk::interpolateMat_8u(img_ref, px_ref[0], px_ref[1]);
            // 更新结果图像（带法向量投影）
            it_normal[width * j + i] = pixel_value;
          }
        }
      }
    }
  }
  // 再次遍历所有体素点，生成密集投影
  for (int i = 0; i < visual_submap->voxel_points.size(); i++)
  {
    VisualPoint *pt = visual_submap->voxel_points[i];

    // 如果法向量未初始化，则跳过
    if (!pt->is_normal_initialized_) continue;

    // 获取参考特征
    Feature *ref_ftr;
    // 投影到当前帧
    V2D pc(new_frame_->w2c(pt->pos_));
    ref_ftr = pt->ref_patch;

    // 定义仿射变换矩阵
    Matrix2d A_cur_ref;
    // 计算仿射变换矩阵（基于深度和变换）
    // 原理：使用深度和相对位姿计算仿射变换
    getWarpMatrixAffine(*cam, ref_ftr->px_, ref_ftr->f_, (ref_ftr->pos() - pt->pos_).norm(), new_frame_->T_f_w_ * ref_ftr->T_f_w_.inverse(), 0, 0,
                        patch_size_half, A_cur_ref);
    // 计算最佳搜索层级
    int search_level = getBestSearchLevel(A_cur_ref.inverse(), 2);
    // 检查变换行列式
    double D = A_cur_ref.determinant();
    if (D > 3) continue;

    // 获取当前帧和参考帧图像
    cv::Mat img_cur = new_frame_->img_;
    cv::Mat img_ref = ref_ftr->img_;
    // 遍历Patch区域
    for (int y = 0; y < patch_size; ++y)
    {
      for (int x = 0; x < patch_size; ++x)
      {
        // 计算Patch内的像素偏移
        Vector2f px_patch(x - patch_size / 2, y - patch_size / 2);
        px_patch *= (1 << search_level);
        // 计算参考帧像素位置
        const Vector2f px_ref(px_patch + ref_ftr->px_.cast<float>());
        // 插值获取像素值
        uint8_t pixel_value = (uint8_t)vk::interpolateMat_8u(img_ref, px_ref[0], px_ref[1]);

        // 通过仿射变换计算当前帧像素位置
        const Vector2f px(A_cur_ref.cast<float>() * px_patch + pc.cast<float>());
        // 如果超出边界，则跳过
        if (px[0] < 0 || px[1] < 0 || px[0] >= img_cur.cols - 1 || px[1] >= img_cur.rows - 1)
          continue;
        else
        {
          // 更新结果图像（密集投影）
          int col = int(px[0]);
          int row = int(px[1]);
          it[width * row + col] = pixel_value;
        }
      }
    }
  }
  // 定义组合图像，用于保存最终结果
  cv::Mat ref_cur_combine;
  cv::Mat ref_cur_combine_normal;
  cv::Mat ref_cur_combine_error;

  // 拼接投影结果和当前帧图像
  cv::hconcat(result, new_frame_->img_, ref_cur_combine);
  cv::hconcat(result_normal, new_frame_->img_, ref_cur_combine_normal);

  // 转换为彩色图像，便于可视化
  cv::cvtColor(ref_cur_combine, ref_cur_combine, CV_GRAY2BGR);
  cv::cvtColor(ref_cur_combine_normal, ref_cur_combine_normal, CV_GRAY2BGR);
  // 计算光度误差
  // 原理：比较投影结果与当前帧的差异
  cv::absdiff(img_photometric_error, result_normal, img_photometric_error);
  cv::hconcat(img_photometric_error, new_frame_->img_, ref_cur_combine_error);

  // 保存结果图像
  cv::imwrite(dir + std::to_string(new_frame_->id_) + "_0_" + ".png", ref_cur_combine);
  cv::imwrite(dir + std::to_string(new_frame_->id_) + "_0_" + "photometric" + ".png", ref_cur_combine_error);
  cv::imwrite(dir + std::to_string(new_frame_->id_) + "_0_" + "normal" + ".png", ref_cur_combine_normal);
}

// 函数：预计算参考Patch的雅可比矩阵
void VIOManager::precomputeReferencePatches(int level)
{
  // 获取当前时间，用于性能测量
  double t1 = omp_get_wtime();
  // 如果总点数为0，则直接返回
  if (total_points == 0) return;

  // 定义雅可比矩阵变量
  MD(1, 2) Jimg;   // 图像梯度雅可比
  MD(2, 3) Jdpi;   // 投影雅可比
  MD(1, 3) Jdphi, Jdp, JdR, Jdt; // 状态相关雅可比

  // 定义H矩阵维度，总点数乘以Patch大小
  const int H_DIM = total_points * patch_size_total;

  // 初始化H_sub_inv矩阵，用于存储逆向雅可比
  H_sub_inv.resize(H_DIM, 6);
  H_sub_inv.setZero();
  // 定义点的反对称矩阵
  M3D p_w_hat;

  // 遍历所有点
  for (int i = 0; i < total_points; i++)
  {
    // 计算当前层级的缩放因子
    const int scale = (1 << level);

    // 获取当前体素点
    VisualPoint *pt = visual_submap->voxel_points[i];
    // 获取参考帧图像
    cv::Mat img = pt->ref_patch->img_;

    // 如果点为空，则跳过
    if (pt == nullptr) continue;

    // 计算点到参考帧相机的深度
    double depth((pt->pos_ - pt->ref_patch->pos()).norm());
    // 计算参考帧中的点坐标
    V3D pf = pt->ref_patch->f_ * depth;
    // 获取参考帧像素坐标
    V2D pc = pt->ref_patch->px_;
    // 获取参考帧旋转矩阵
    M3D R_ref_w = pt->ref_patch->T_f_w_.rotation_matrix();

    // 计算投影雅可比矩阵
    // 原理：描述三维点到二维像素的投影关系
    computeProjectionJacobian(pf, Jdpi);
    // 计算点的反对称矩阵
    p_w_hat << SKEW_SYM_MATRX(pt->pos_);

    // 获取参考帧像素坐标的浮点值
    const float u_ref = pc[0];
    const float v_ref = pc[1];
    // 计算整数像素坐标
    const int u_ref_i = floorf(pc[0] / scale) * scale;
    const int v_ref_i = floorf(pc[1] / scale) * scale;
    // 计算子像素偏移
    const float subpix_u_ref = (u_ref - u_ref_i) / scale;
    const float subpix_v_ref = (v_ref - v_ref_i) / scale;
    // 计算双线性插值权重
    const float w_ref_tl = (1.0 - subpix_u_ref) * (1.0 - subpix_v_ref);
    const float w_ref_tr = subpix_u_ref * (1.0 - subpix_v_ref);
    const float w_ref_bl = (1.0 - subpix_u_ref) * subpix_v_ref;
    const float w_ref_br = subpix_u_ref * subpix_v_ref;

    // 遍历Patch区域，计算雅可比
    for (int x = 0; x < patch_size; x++)
    {
      // 计算图像数据指针位置
      uint8_t *img_ptr = (uint8_t *)img.data + (v_ref_i + x * scale - patch_size_half * scale) * width + u_ref_i - patch_size_half * scale;
      for (int y = 0; y < patch_size; ++y, img_ptr += scale)
      {
        // 计算水平方向梯度（双线性插值）
        float du = 0.5f * (
            (w_ref_tl * img_ptr[scale] + w_ref_tr * img_ptr[scale * 2] + w_ref_bl * img_ptr[scale * width + scale] +
             w_ref_br * img_ptr[scale * width + scale * 2]) -
            (w_ref_tl * img_ptr[-scale] + w_ref_tr * img_ptr[0] + w_ref_bl * img_ptr[scale * width - scale] + 
             w_ref_br * img_ptr[scale * width]));
        // 计算垂直方向梯度（双线性插值）
        float dv = 0.5f * (
            (w_ref_tl * img_ptr[scale * width] + w_ref_tr * img_ptr[scale + scale * width] + w_ref_bl * img_ptr[width * scale * 2] +
             w_ref_br * img_ptr[width * scale * 2 + scale]) -
            (w_ref_tl * img_ptr[-scale * width] + w_ref_tr * img_ptr[-scale * width + scale] + w_ref_bl * img_ptr[0] + 
             w_ref_br * img_ptr[scale]));

        // 组合图像梯度雅可比
        Jimg << du, dv;
        // 缩放梯度，考虑层级
        Jimg = Jimg * (1.0 / scale);

        // 计算旋转相关的雅可比
        JdR = Jimg * Jdpi * R_ref_w * p_w_hat;
        // 计算平移相关的雅可比
        Jdt = -Jimg * Jdpi * R_ref_w;

        // 将雅可比存入H_sub_inv矩阵
        H_sub_inv.block<1, 6>(i * patch_size_total + x * patch_size + y, 0) << JdR, Jdt;
      }
    }
  }
  // 标记参考Patch缓存已生成
  has_ref_patch_cache = true;
}

// 函数：使用逆向方法更新状态
void VIOManager::updateStateInverse(cv::Mat img, int level)
{
  // 如果总点数为0，则直接返回
  if (total_points == 0) return;

  // 保存旧状态，用于回滚
  StatesGroup old_state = (*state);
  // 定义像素坐标变量
  V2D pc;
  // 定义雅可比矩阵
  MD(1, 2) Jimg;
  MD(2, 3) Jdpi;
  MD(1, 3) Jdphi, Jdp, JdR, Jdt;
  // 定义测量向量和H矩阵
  VectorXd z;
  MatrixXd H_sub;
  // 标记EKF是否结束
  bool EKF_end = false;
  // 初始化上一次误差为最大值
  float last_error = std::numeric_limits<float>::max();
  // 初始化时间记录变量
  compute_jacobian_time = update_ekf_time = 0.0;
  // 定义位置的反对称矩阵
  M3D P_wi_hat;
  // 标记测量向量是否初始化
  bool z_init = true;
  // 定义H矩阵维度
  const int H_DIM = total_points * patch_size_total;

  // 初始化测量向量
  z.resize(H_DIM);
  z.setZero();

  // 初始化H矩阵
  H_sub.resize(H_DIM, 6);
  H_sub.setZero();

  // 迭代更新状态
  for (int iteration = 0; iteration < max_iterations; iteration++)
  {
    // 记录开始时间
    double t1 = omp_get_wtime();
    // 初始化外点计数器
    double count_outlier = 0;
    // 如果未缓存参考Patch，则预计算
    if (has_ref_patch_cache == false) precomputeReferencePatches(level);
    // 初始化测量计数和误差
    int n_meas = 0;
    float error = 0.0;
    // 获取当前旋转和平移
    M3D Rwi(state->rot_end);
    V3D Pwi(state->pos_end);
    // 计算相机到世界的变换
    P_wi_hat << SKEW_SYM_MATRX(Pwi);
    Rcw = Rci * Rwi.transpose();
    Pcw = -Rci * Rwi.transpose() * Pwi + Pci;

    // 定义点的反对称矩阵
    M3D p_hat;

    // 遍历所有点
    for (int i = 0; i < total_points; i++)
    {
      // 初始化Patch误差
      float patch_error = 0.0;

      // 计算当前层级缩放因子
      const int scale = (1 << level);

      // 获取当前体素点
      VisualPoint *pt = visual_submap->voxel_points[i];

      // 如果点为空，则跳过
      if (pt == nullptr) continue;

      // 将点投影到相机坐标系
      V3D pf = Rcw * pt->pos_ + Pcw;
      // 转换为像素坐标
      pc = cam->world2cam(pf);

      // 获取像素坐标的浮点值
      const float u_ref = pc[0];
      const float v_ref = pc[1];
      // 计算整数像素坐标
      const int u_ref_i = floorf(pc[0] / scale) * scale;
      const int v_ref_i = floorf(pc[1] / scale) * scale;
      // 计算子像素偏移
      const float subpix_u_ref = (u_ref - u_ref_i) / scale;
      const float subpix_v_ref = (v_ref - v_ref_i) / scale;
      // 计算双线性插值权重
      const float w_ref_tl = (1.0 - subpix_u_ref) * (1.0 - subpix_v_ref);
      const float w_ref_tr = subpix_u_ref * (1.0 - subpix_v_ref);
      const float w_ref_bl = (1.0 - subpix_u_ref) * subpix_v_ref;
      const float w_ref_br = subpix_u_ref * subpix_v_ref;

      // 获取参考Patch数据
      vector<float> P = visual_submap->warp_patch[i];
      // 遍历Patch区域
      for (int x = 0; x < patch_size; x++)
      {
        // 计算图像数据指针位置
        uint8_t *img_ptr = (uint8_t *)img.data + (v_ref_i + x * scale - patch_size_half * scale) * width + u_ref_i - patch_size_half * scale;
        for (int y = 0; y < patch_size; ++y, img_ptr += scale)
        {
          // 计算当前像素值（双线性插值）
          double res = w_ref_tl * img_ptr[0] + w_ref_tr * img_ptr[scale] + w_ref_bl * img_ptr[scale * width] +
                       w_ref_br * img_ptr[scale * width + scale] - P[patch_size_total * level + x * patch_size + y];
          // 保存残差
          z(i * patch_size_total + x * patch_size + y) = res;
          // 累加Patch误差
          patch_error += res * res;
          // 获取预计算的雅可比
          MD(1, 3) J_dR = H_sub_inv.block<1, 3>(i * patch_size_total + x * patch_size + y, 0);
          MD(1, 3) J_dt = H_sub_inv.block<1, 3>(i * patch_size_total + x * patch_size + y, 3);
          // 更新雅可比
          JdR = J_dR * Rwi + J_dt * P_wi_hat * Rwi;
          Jdt = J_dt * Rwi;
          // 保存到H矩阵
          H_sub.block<1, 6>(i * patch_size_total + x * patch_size + y, 0) << JdR, Jdt;
          n_meas++;
        }
      }
      // 保存Patch误差
      visual_submap->errors[i] = patch_error;
      // 累加总误差
      error += patch_error;
    }

    // 计算平均误差
    error = error / n_meas;

    // 记录雅可比计算时间
    compute_jacobian_time += omp_get_wtime() - t1;

    // 记录EKF开始时间
    double t3 = omp_get_wtime();

    // 如果误差减小，则更新状态
    if (error <= last_error)
    {
      old_state = (*state);
      last_error = error;

      // 计算H矩阵的转置
      auto &&H_sub_T = H_sub.transpose();
      // 初始化H^T*H和G矩阵
      H_T_H.setZero();
      G.setZero();
      // 计算H^T*H
      H_T_H.block<6, 6>(0, 0) = H_sub_T * H_sub;
      // 计算卡尔曼增益的逆
      MD(DIM_STATE, DIM_STATE) &&K_1 = (H_T_H + (state->cov / img_point_cov).inverse()).inverse();
      // 计算H^T*z
      auto &&HTz = H_sub_T * z;
      // 计算状态差
      auto vec = (*state_propagat) - (*state);
      // 计算增益矩阵G
      G.block<DIM_STATE, 6>(0, 0) = K_1.block<DIM_STATE, 6>(0, 0) * H_T_H.block<6, 6>(0, 0);
      // 计算状态更新量
      auto solution = -K_1.block<DIM_STATE, 6>(0, 0) * HTz + vec - G.block<DIM_STATE, 6>(0, 0) * vec.block<6, 1>(0, 0);
      // 更新状态
      (*state) += solution;
      // 获取旋转和平移更新量
      auto &&rot_add = solution.block<3, 1>(0, 0);
      auto &&t_add = solution.block<3, 1>(3, 0);

      // 如果更新量足够小，则认为EKF收敛
      if ((rot_add.norm() * 57.3f < 0.001f) && (t_add.norm() * 100.0f < 0.001f)) { EKF_end = true; }
    }
    else
    {
      // 如果误差增大，则回滚状态
      (*state) = old_state;
      EKF_end = true;
    }

    // 记录EKF更新时间
    update_ekf_time += omp_get_wtime() - t3;

    // 如果达到最大迭代次数或EKF结束，则退出循环
    if (iteration == max_iterations || EKF_end) break;
  }
}

// 函数：更新状态（正向方法）
void VIOManager::updateState(cv::Mat img, int level)
{
  // 如果总点数为0，则直接返回
  if (total_points == 0) return;

  // 保存旧状态，用于回滚
  StatesGroup old_state = (*state);

  // 定义测量向量和H矩阵
  VectorXd z;
  MatrixXd H_sub;
  // 标记EKF是否结束
  bool EKF_end = false;
  // 初始化上一次误差为最大值
  float last_error = std::numeric_limits<float>::max();

  // 定义H矩阵维度
  const int H_DIM = total_points * patch_size_total;
  // 初始化测量向量
  z.resize(H_DIM);
  z.setZero();
  // 初始化H矩阵（考虑曝光时间则为7维，否则为6维）
  H_sub.resize(H_DIM, 7);
  H_sub.setZero();

  // 迭代更新状态
  for (int iteration = 0; iteration < max_iterations; iteration++)
  {
    // 记录开始时间
    double t1 = omp_get_wtime();

    // 获取当前旋转和平移
    M3D Rwi(state->rot_end);
    V3D Pwi(state->pos_end);
    // 计算相机到世界的变换
    Rcw = Rci * Rwi.transpose();
    Pcw = -Rci * Rwi.transpose() * Pwi + Pci;
    // 计算平移雅可比
    Jdp_dt = Rci * Rwi.transpose();

    // 初始化误差和测量计数
    float error = 0.0;
    int n_meas = 0;

    // 并行遍历所有点（如果启用多线程）
    #ifdef MP_EN
      omp_set_num_threads(MP_PROC_NUM);
      #pragma omp parallel for reduction(+:error, n_meas)
    #endif
    for (int i = 0; i < total_points; i++)
    {
      // 定义雅可比矩阵
      MD(1, 2) Jimg;
      MD(2, 3) Jdpi;
      MD(1, 3) Jdphi, Jdp, JdR, Jdt;

      // 初始化Patch误差
      float patch_error = 0.0;
      // 获取搜索层级和金字塔层级
      int search_level = visual_submap->search_levels[i];
      int pyramid_level = level + search_level;
      int scale = (1 << pyramid_level);
      float inv_scale = 1.0f / scale;

      // 获取当前体素点
      VisualPoint *pt = visual_submap->voxel_points[i];

      // 如果点为空，则跳过
      if (pt == nullptr) continue;

      // 将点投影到相机坐标系
      V3D pf = Rcw * pt->pos_ + Pcw;
      // 转换为像素坐标
      V2D pc = cam->world2cam(pf);

      // 计算投影雅可比
      computeProjectionJacobian(pf, Jdpi);
      // 计算点的反对称矩阵
      M3D p_hat;
      p_hat << SKEW_SYM_MATRX(pf);

      // 获取像素坐标的浮点值
      float u_ref = pc[0];
      float v_ref = pc[1];
      // 计算整数像素坐标
      int u_ref_i = floorf(pc[0] / scale) * scale;
      int v_ref_i = floorf(pc[1] / scale) * scale;
      // 计算子像素偏移
      float subpix_u_ref = (u_ref - u_ref_i) / scale;
      float subpix_v_ref = (v_ref - v_ref_i) / scale;
      // 计算双线性插值权重
      float w_ref_tl = (1.0 - subpix_u_ref) * (1.0 - subpix_v_ref);
      float w_ref_tr = subpix_u_ref * (1.0 - subpix_v_ref);
      float w_ref_bl = (1.0 - subpix_u_ref) * subpix_v_ref;
      float w_ref_br = subpix_u_ref * subpix_v_ref;

      // 获取参考Patch数据和曝光时间
      vector<float> P = visual_submap->warp_patch[i];
      double inv_ref_expo = visual_submap->inv_expo_list[i];

      // 遍历Patch区域
      for (int x = 0; x < patch_size; x++)
      {
        // 计算图像数据指针位置
        uint8_t *img_ptr = (uint8_t *)img.data + (v_ref_i + x * scale - patch_size_half * scale) * width + u_ref_i - patch_size_half * scale;
        for (int y = 0; y < patch_size; ++y, img_ptr += scale)
        {
          // 计算水平方向梯度
          float du = 0.5f * (
              (w_ref_tl * img_ptr[scale] + w_ref_tr * img_ptr[scale * 2] + w_ref_bl * img_ptr[scale * width + scale] +
               w_ref_br * img_ptr[scale * width + scale * 2]) -
              (w_ref_tl * img_ptr[-scale] + w_ref_tr * img_ptr[0] + w_ref_bl * img_ptr[scale * width - scale] + 
               w_ref_br * img_ptr[scale * width]));
          // 计算垂直方向梯度
          float dv = 0.5f * (
              (w_ref_tl * img_ptr[scale * width] + w_ref_tr * img_ptr[scale + scale * width] + w_ref_bl * img_ptr[width * scale * 2] +
               w_ref_br * img_ptr[width * scale * 2 + scale]) -
              (w_ref_tl * img_ptr[-scale * width] + w_ref_tr * img_ptr[-scale * width + scale] + w_ref_bl * img_ptr[0] + 
               w_ref_br * img_ptr[scale]));

          // 组合图像梯度雅可比，考虑曝光时间
          Jimg << du, dv;
          Jimg = Jimg * state->inv_expo_time;
          Jimg = Jimg * inv_scale;
          // 计算状态相关的雅可比
          Jdphi = Jimg * Jdpi * p_hat;
          Jdp = -Jimg * Jdpi;
          JdR = Jdphi * Jdphi_dR + Jdp * Jdp_dR;
          Jdt = Jdp * Jdp_dt;

          // 计算当前像素值和残差
          double cur_value = w_ref_tl * img_ptr[0] + w_ref_tr * img_ptr[scale] + w_ref_bl * img_ptr[scale * width] + 
                            w_ref_br * img_ptr[scale * width + scale];
          double res = state->inv_expo_time * cur_value - inv_ref_expo * P[patch_size_total * level + x * patch_size + y];

          // 保存残差
          z(i * patch_size_total + x * patch_size + y) = res;

          // 累加Patch误差
          patch_error += res * res;
          n_meas += 1;

          // 根据是否估计曝光时间，保存雅可比
          if (exposure_estimate_en) { H_sub.block<1, 7>(i * patch_size_total + x * patch_size + y, 0) << JdR, Jdt, cur_value; }
          else { H_sub.block<1, 6>(i * patch_size_total + x * patch_size + y, 0) << JdR, Jdt; }
        }
      }
      // 保存Patch误差
      visual_submap->errors[i] = patch_error;
      // 累加总误差
      error += patch_error;
    }

    // 计算平均误差
    error = error / n_meas;

    // 记录雅可比计算时间
    compute_jacobian_time += omp_get_wtime() - t1;

    // 记录EKF开始时间
    double t3 = omp_get_wtime();

    // 如果误差减小，则更新状态
    if (error <= last_error)
    {
      old_state = (*state);
      last_error = error;

      // 计算H矩阵的转置
      auto &&H_sub_T = H_sub.transpose();
      // 初始化H^T*H和G矩阵
      H_T_H.setZero();
      G.setZero();
      // 计算H^T*H
      H_T_H.block<7, 7>(0, 0) = H_sub_T * H_sub;
      // 计算卡尔曼增益的逆
      MD(DIM_STATE, DIM_STATE) &&K_1 = (H_T_H + (state->cov / img_point_cov).inverse()).inverse();
      // 计算H^T*z
      auto &&HTz = H_sub_T * z;
      // 计算状态差
      auto vec = (*state_propagat) - (*state);
      // 计算增益矩阵G
      G.block<DIM_STATE, 7>(0, 0) = K_1.block<DIM_STATE, 7>(0, 0) * H_T_H.block<7, 7>(0, 0);
      // 计算状态更新量
      MD(DIM_STATE, 1) solution = -K_1.block<DIM_STATE, 7>(0, 0) * HTz + vec - G.block<DIM_STATE, 7>(0, 0) * vec.block<7, 1>(0, 0);

      // 更新状态
      (*state) += solution;
      // 获取旋转、平移和曝光更新量
      auto &&rot_add = solution.block<3, 1>(0, 0);
      auto &&t_add = solution.block<3, 1>(3, 0);
      auto &&expo_add = solution.block<1, 1>(6, 0);

      // 如果更新量足够小，则认为EKF收敛
      if ((rot_add.norm() * 57.3f < 0.001f) && (t_add.norm() * 100.0f < 0.001f)) EKF_end = true;
    }
    else
    {
      // 如果误差增大，则回滚状态
      (*state) = old_state;
      EKF_end = true;
    }

    // 记录EKF更新时间
    update_ekf_time += omp_get_wtime() - t3;

    // 如果达到最大迭代次数或EKF结束，则退出循环
    if (iteration == max_iterations || EKF_end) break;
  }
}

// 函数：更新帧状态
void VIOManager::updateFrameState(StatesGroup state)
{
  // 获取当前旋转和平移
  M3D Rwi(state.rot_end);
  V3D Pwi(state.pos_end);
  // 计算相机到世界的变换
  Rcw = Rci * Rwi.transpose();
  Pcw = -Rci * Rwi.transpose() * Pwi + Pci;
  // 更新当前帧的位姿
  new_frame_->T_f_w_ = SE3(Rcw, Pcw);
}

// 函数：绘制跟踪点
void VIOManager::plotTrackedPoints()
{
  // 获取总点数
  int total_points = visual_submap->voxel_points.size();
  // 如果总点数为0，则直接返回
  if (total_points == 0) return;

  // 遍历所有点
  for (int i = 0; i < total_points; i++)
  {
    // 获取当前体素点
    VisualPoint *pt = visual_submap->voxel_points[i];
    // 投影到当前帧
    V2D pc(new_frame_->w2c(pt->pos_));

    // 如果当前误差小于先验误差，则为内点
    if (visual_submap->errors[i] <= visual_submap->propa_errors[i])
    {
      // 绘制绿色圆点表示内点
      cv::circle(img_cp, cv::Point2f(pc[0], pc[1]), 7, cv::Scalar(0, 255, 0), -1, 8);
    }
    else
    {
      // 绘制红色圆点表示外点
      cv::circle(img_cp, cv::Point2f(pc[0], pc[1]), 7, cv::Scalar(255, 0, 0), -1, 8);
    }
  }
}

// 函数：获取插值像素值（彩色图像）
V3F VIOManager::getInterpolatedPixel(cv::Mat img, V2D pc)
{
  // 获取像素坐标的浮点值
  const float u_ref = pc[0];
  const float v_ref = pc[1];
  // 计算整数像素坐标
  const int u_ref_i = floorf(pc[0]);
  const int v_ref_i = floorf(pc[1]);
  // 计算子像素偏移
  const float subpix_u_ref = (u_ref - u_ref_i);
  const float subpix_v_ref = (v_ref - v_ref_i);
  // 计算双线性插值权重
  const float w_ref_tl = (1.0 - subpix_u_ref) * (1.0 - subpix_v_ref);
  const float w_ref_tr = subpix_u_ref * (1.0 - subpix_v_ref);
  const float w_ref_bl = (1.0 - subpix_u_ref) * subpix_v_ref;
  const float w_ref_br = subpix_u_ref * subpix_v_ref;
  // 计算图像数据指针位置（彩色图像，3通道）
  uint8_t *img_ptr = (uint8_t *)img.data + ((v_ref_i) * width + (u_ref_i)) * 3;
  // 计算B、G、R通道的插值像素值
  float B = w_ref_tl * img_ptr[0] + w_ref_tr * img_ptr[0 + 3] + w_ref_bl * img_ptr[width * 3] + w_ref_br * img_ptr[width * 3 + 0 + 3];
  float G = w_ref_tl * img_ptr[1] + w_ref_tr * img_ptr[1 + 3] + w_ref_bl * img_ptr[1 + width * 3] + w_ref_br * img_ptr[width * 3 + 1 + 3];
  float R = w_ref_tl * img_ptr[2] + w_ref_tr * img_ptr[2 + 3] + w_ref_bl * img_ptr[2 + width * 3] + w_ref_br * img_ptr[width * 3 + 2 + 3];
  // 返回插值像素值
  V3F pixel(B, G, R);
  return pixel;
}

// 函数：为Colmap导出数据
void VIOManager::dumpDataForColmap()
{
  // 定义静态计数器，用于生成文件名
  static int cnt = 1;
  // 格式化计数器为5位数字
  std::ostringstream ss;
  ss << std::setw(5) << std::setfill('0') << cnt;
  std::string cnt_str = ss.str();
  // 定义图像保存路径
  std::string image_path = std::string(ROOT_DIR) + "Log/Colmap/images/" + cnt_str + ".png";

  // 去畸变并保存图像
  cv::Mat img_rgb_undistort;
  pinhole_cam->undistortImage(img_rgb, img_rgb_undistort);
  cv::imwrite(image_path, img_rgb_undistort);

  // 获取当前帧的四元数和平移
  Eigen::Quaterniond q(new_frame_->T_f_w_.rotation_matrix());
  Eigen::Vector3d t = new_frame_->T_f_w_.translation();
  // 写入Colmap格式数据
  fout_colmap << cnt << " "
              << std::fixed << std::setprecision(6)  // 保证浮点数精度为6位
              << q.w() << " " << q.x() << " " << q.y() << " " << q.z() << " "
              << t.x() << " " << t.y() << " " << t.z() << " "
              << 1 << " "  // 相机ID（假设为1）
              << cnt_str << ".png" << std::endl;
  fout_colmap << "0.0 0.0 -1" << std::endl;
  // 计数器加1
  cnt++;
}

// 函数：处理一帧图像
void VIOManager::processFrame(cv::Mat &img, vector<pointWithVar> &pg, const unordered_map<VOXEL_LOCATION, VoxelOctoTree *> &feat_map, double img_time)
{
  // 检查图像尺寸是否匹配，若不匹配则调整大小
  if (width != img.cols || height != img.rows)
  {
    // 如果图像为空，则打印警告
    if (img.empty()) printf("[ VIO ] Empty Image!\n");
    // 调整图像大小
    cv::resize(img, img, cv::Size(img.cols * image_resize_factor, img.rows * image_resize_factor), 0, 0, CV_INTER_LINEAR);
  }
  // 克隆原始彩色图像
  img_rgb = img.clone();
  // 克隆图像用于绘制
  img_cp = img.clone();

  // 如果是彩色图像，则转换为灰度图
  if (img.channels() == 3) cv::cvtColor(img, img, CV_BGR2GRAY);

  // 创建新帧对象
  new_frame_.reset(new Frame(cam, img));
  // 更新帧状态
  updateFrameState(*state);

  // 重置网格（未定义，可能用于特征管理）
  resetGrid();

  // 记录开始时间
  double t1 = omp_get_wtime();

  // 从视觉稀疏地图中检索数据
  retrieveFromVisualSparseMap(img, pg, feat_map);

  // 记录检索结束时间
  double t2 = omp_get_wtime();

  // 计算雅可比并更新EKF
  computeJacobianAndUpdateEKF(img);

  // 记录EKF结束时间
  double t3 = omp_get_wtime();

  // 生成视觉地图点
  generateVisualMapPoints(img, pg);

  // 记录生成结束时间
  double t4 = omp_get_wtime();

  // 绘制跟踪点
  plotTrackedPoints();

  // 如果启用投影，则执行Patch投影
  if (plot_flag) projectPatchFromRefToCur(feat_map);

  // 记录投影结束时间
  double t5 = omp_get_wtime();

  // 更新视觉地图点
  updateVisualMapPoints(img);

  // 记录更新结束时间
  double t6 = omp_get_wtime();

  // 更新参考Patch
  updateReferencePatch(feat_map);

  // 记录参考Patch更新结束时间
  double t7 = omp_get_wtime();

  // 如果启用Colmap输出，则导出数据
  if (colmap_output_en) dumpDataForColmap();

  // 帧计数器加1
  frame_count++;
  // 计算平均总时间（排除投影时间）
  ave_total = ave_total * (frame_count - 1) / frame_count + (t7 - t1 - (t5 - t4)) / frame_count;

  // 打印时间统计信息
  printf("\033[1;34m+-------------------------------------------------------------+\033[0m\n");
  printf("\033[1;34m|                         VIO Time                            |\033[0m\n");
  printf("\033[1;34m+-------------------------------------------------------------+\033[0m\n");
  printf("\033[1;34m| %-29s | %-27zu |\033[0m\n", "Sparse Map Size", feat_map.size());
  printf("\033[1;34m+-------------------------------------------------------------+\033[0m\n");
  printf("\033[1;34m| %-29s | %-27s |\033[0m\n", "Algorithm Stage", "Time (secs)");
  printf("\033[1;34m+-------------------------------------------------------------+\033[0m\n");
  printf("\033[1;32m| %-29s | %-27lf |\033[0m\n", "retrieveFromVisualSparseMap", t2 - t1);
  printf("\033[1;32m| %-29s | %-27lf |\033[0m\n", "computeJacobianAndUpdateEKF", t3 - t2);
  printf("\033[1;32m| %-27s   | %-27lf |\033[0m\n", "-> computeJacobian", compute_jacobian_time);
  printf("\033[1;32m| %-27s   | %-27lf |\033[0m\n", "-> updateEKF", update_ekf_time);
  printf("\033[1;32m| %-29s | %-27lf |\033[0m\n", "generateVisualMapPoints", t4 - t3);
  printf("\033[1;32m| %-29s | %-27lf |\033[0m\n", "updateVisualMapPoints", t6 - t5);
  printf("\033[1;32m| %-29s | %-27lf |\033[0m\n", "updateReferencePatch", t7 - t6);
  printf("\033[1;34m+-------------------------------------------------------------+\033[0m\n");
  printf("\033[1;32m| %-29s | %-27lf |\033[0m\n", "Current Total Time", t7 - t1 - (t5 - t4));
  printf("\033[1;32m| %-29s | %-27lf |\033[0m\n", "Average Total Time", ave_total);
  printf("\033[1;34m+-------------------------------------------------------------+\033[0m\n");
}
