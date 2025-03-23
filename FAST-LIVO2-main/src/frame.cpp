// 包含必要的头文件
#include <boost/bind.hpp>          // Boost 库的 bind 功能，用于绑定回调函数
#include "feature.h"               // 自定义头文件，定义特征点类 Feature
#include "frame.h"                 // 自定义头文件，声明 Frame 类
#include "visual_point.h"          // 自定义头文件，定义视觉点类
#include <stdexcept>               // 标准异常库，用于抛出运行时错误
#include <vikit/math_utils.h>      // VIKit 数学工具库，提供数学计算支持
#include <vikit/performance_monitor.h> // VIKit 性能监控工具，用于性能分析
#include <vikit/vision.h>          // VIKit 视觉工具库，提供图像处理函数

// 初始化静态成员变量 frame_counter_，用于为每个 Frame 对象分配唯一 ID
int Frame::frame_counter_ = 0;

/*
 * Frame 类的构造函数
 * 功能：初始化一个图像帧对象，关联相机模型并处理输入图像。
 * 参数：
 *   - cam: 指向抽象相机模型的指针，用于投影和图像尺寸校验。
 *   - img: 输入的图像数据 (cv::Mat 类型)，通常为灰度图像。
 * 原理：通过自增的 frame_counter_ 分配唯一 ID，并调用 initFrame 初始化帧。
 */
Frame::Frame(vk::AbstractCamera *cam, const cv::Mat &img)
    : id_(frame_counter_++),  // 为当前帧分配唯一 ID，frame_counter_ 自增
      cam_(cam)               // 保存相机模型指针，用于后续投影和校验
{
  initFrame(img);  // 调用初始化函数处理输入图像
}

/*
 * Frame 类的析构函数
 * 功能：释放 Frame 对象占用的资源，特别是特征点的内存。
 * 原理：遍历特征点容器 fts_，逐一删除每个特征点的内存，避免内存泄漏。
 */
Frame::~Frame()
{
  // 使用 std::for_each 遍历特征点向量 fts_，并删除每个特征点的动态内存
  std::for_each(fts_.begin(), fts_.end(), [&](Feature *i) { delete i; });
}

/*
 * 初始化帧的函数
 * 功能：验证输入图像并将其赋值给成员变量 img_，为后续特征提取做准备。
 * 参数：
 *   - img: 输入的图像数据 (cv::Mat 类型)。
 * 原理：通过一系列检查确保图像有效（非空、尺寸匹配、灰度格式），然后保存图像。
 * 异常处理：若图像不符合要求，抛出运行时错误。
 */
void Frame::initFrame(const cv::Mat &img)
{
  // 检查图像是否为空，若为空则抛出异常
  if (img.empty()) { throw std::runtime_error("Frame: provided image is empty"); }

  // 检查图像尺寸是否与相机模型匹配，若不匹配则抛出异常
  // 原理：相机模型定义了期望的图像宽度和高度，输入图像必须一致以确保投影正确
  if (img.cols != cam_->width() || img.rows != cam_->height())
  {
    throw std::runtime_error("Frame: provided image has not the same size as the camera model");
  }

  // 检查图像是否为灰度格式 (CV_8UC1)，若不是则抛出异常
  // 原理：视觉里程计通常使用灰度图像以减少计算量并提高鲁棒性
  if (img.type() != CV_8UC1) { throw std::runtime_error("Frame: provided image is not grayscale"); }

  // 将输入图像赋值给成员变量 img_，后续用于特征提取或金字塔构建
  img_ = img;
}

/// Frame 类的工具函数命名空间
namespace frame_utils
{

/*
 * 创建图像金字塔的函数
 * 功能：从原始图像生成多层图像金字塔，用于多尺度特征检测和匹配。
 * 参数：
 *   - img_level_0: 第 0 层图像（原始图像，灰度格式）。
 *   - n_levels: 金字塔层数。
 *   - pyr: 输出图像金字塔容器 (ImgPyr 类型，通常为 std::vector<cv::Mat>)。
 * 原理：
 *   - 图像金字塔通过逐层降采样生成，每层分辨率减半。
 *   - 多尺度图像有助于在不同分辨率下检测特征，提高鲁棒性和精度。
 *   - 使用 vk::halfSample 函数实现高效降采样。
 */
void createImgPyramid(const cv::Mat &img_level_0, int n_levels, ImgPyr &pyr)
{
  // 调整金字塔容器大小以容纳指定层数
  pyr.resize(n_levels);

  // 将第 0 层设置为原始图像
  pyr[0] = img_level_0;

  // 循环生成后续层级的图像
  for (int i = 1; i < n_levels; ++i)
  {
    // 为当前层分配内存，分辨率为上一层的一半
    // 原理：每层宽度和高度减半，保持灰度格式 (CV_8U)
    pyr[i] = cv::Mat(pyr[i - 1].rows / 2, pyr[i - 1].cols / 2, CV_8U);

    // 使用 VIKit 的 halfSample 函数从上一层降采样到当前层
    // 原理：halfSample 通常通过平均或插值方式减少像素，保留主要特征
    vk::halfSample(pyr[i - 1], pyr[i]);
  }
}

} // namespace frame_utils
