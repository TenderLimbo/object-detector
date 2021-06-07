#include "Manager.h"

Manager::Manager(std::string window_name_) {
    window_name = std::move(window_name_);
    cv::namedWindow(window_name, cv::WINDOW_AUTOSIZE);
}

int Manager::read(const std::string &name) {
  src = imread(name, cv::IMREAD_COLOR);
  if (src.empty()) {
    std::cout << "Could not open or find the image" << std::endl;
    return -1;
  }
  return 0;
}

std::vector<std::vector<Point3> > Manager::getMatrix() const {
  std::vector<std::vector<Point3>> matrix(src.rows, std::vector<Point3>(src.cols));
  for (int i = 0; i < src.rows; i++) {
    for (int j = 0; j < src.cols; j++) {
      cv::Vec3b vec = src.at<cv::Vec3b>(i, j);
      matrix[i][j] = {vec[0], vec[1], vec[2]};
    }
  }
  return matrix;
}

void Manager::show(const std::vector<std::vector<uint8_t>> &matrix) {
  dst = cv::Mat(src.rows,src.cols,CV_8UC1);
  for (int i = 0; i < dst.rows; i++) {
    for (int j = 0; j < dst.cols; j++) {
      dst.at<uchar>(i,j) = matrix[i][j];
    }
  }
  cv::imshow(window_name, dst);
}

void Manager::showRes(const std::vector<std::vector<Point3>> &matrix) {
  dst = cv::Mat(matrix.size(),matrix[0].size(), CV_8UC3);
  for (int i = 0; i < dst.rows; i++) {
    for (int j = 0; j < dst.cols; j++) {
      dst.at<cv::Vec3b>(i, j)[0] = matrix[i][j].blue;
      dst.at<cv::Vec3b>(i, j)[1] = matrix[i][j].green;
      dst.at<cv::Vec3b>(i, j)[2] = matrix[i][j].red;
    }
  }
  cv::imshow(window_name, dst);
}

Manager::~Manager() {
  cv::waitKey(0);
  cv::destroyWindow(window_name);
}