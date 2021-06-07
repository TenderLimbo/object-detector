#pragma once

#include <opencv2/opencv.hpp>
#include <iostream>
#include <vector>

struct Point3 {
    uint8_t blue = 0;
    uint8_t green = 0;
    uint8_t red = 0;
};

struct Point {
    int x = 0;
    int y = 0;
};

struct Rect {
    Point down_left;
    Point up_right;
};

class Manager {
    cv::Mat src;
    cv::Mat dst;
    std::string window_name;
public:
    Manager(std::string window_name_);
    int read(const std::string &name);
    std::vector<std::vector<Point3>> getMatrix() const;
    void show(const std::vector<std::vector<uint8_t>> &matrix);
    void showRes(const std::vector<std::vector<Point3>> &matrix);
    ~Manager();
};

