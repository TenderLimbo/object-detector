#include <cmath>
#include <climits>
#include "Manager.h"

using namespace std;

template<typename T>
using Matrix = vector<vector<T>>;

double GaussFunc(int x, int y, double sigma) {
  return exp(-(x * x + y * y) / (2 * sigma * sigma)) / (2 * M_PI * sigma * sigma);
}

template<typename T, typename U>
Matrix<T> ElementWiseMultiplication(const Matrix<T> &m, const Matrix<U> &n) {
  Matrix<T> res(m.size(), vector<T>(m[0].size()));
  for (int i = 0; i < m.size(); i++) {
    for (int j = 0; j < m[i].size(); j++) {
      res[i][j] = m[i][j] * static_cast<T>(n[i][j]);
    }
  }
  return res;
}

template<typename T>
T SumMatrixElements(const Matrix<T> &matrix) {
  T sum = 0;
  for (auto &line : matrix)
    for (auto &el : line)
      sum += el;
  return sum;
}

template<typename T>
Matrix<T> GetSubmatrix(const Matrix<T> &matrix, int x, int y, int height, int width) {
  Matrix<T> submatrix(height, vector<T>(width));
  for (int i = x; i < x + height; ++i) {
    for (int j = y; j < y + width; ++j) {
      submatrix[i - x][j - y] = matrix[i][j];
    }
  }
  return submatrix;
}

Matrix<double> CreateKern(int size, double sigma) {
  int shift = (int) ((size - 1) / 2);
  Matrix<double> kernel(size, vector<double>(size));
  double sum = 0;
  for (int i = -shift; i <= shift; ++i) {
    for (int j = -shift; j <= shift; ++j) {
      kernel[i + shift][j + shift] = GaussFunc(i, j, sigma);
      sum += kernel[i + shift][j + shift];
    }
  }
  for (auto &line : kernel) {
    for (auto &el : line) {
      el /= sum;
    }
  }
  return kernel;
}

Matrix<uint8_t> GaussFilter(const Matrix<uint8_t> &src, int size, double sigma) {
  Matrix<double> kernel = CreateKern(size, sigma);
  int shift = (int) ((size - 1) / 2);

  Matrix<double> extended_matrix(shift * 2 + src.size(), vector<double>(shift * 2 + src[0].size()));
  for (int i = shift; i < src.size(); i++) {
    for (int j = shift; j < src[i].size(); j++) {
      extended_matrix[i][j] = static_cast<double>(src[i - shift][j - shift]);
    }
  }
  Matrix<uint8_t> out(src);
  for (int i = 0; i < out.size(); ++i) {
    for (int j = 0; j < out[i].size(); ++j) {
      double val = SumMatrixElements(ElementWiseMultiplication(kernel,
                                                               GetSubmatrix(extended_matrix, i, j, size, size)));
      out[i][j] = static_cast<uint8_t>(val);
    }
  }
  return out;
}

Matrix<uint8_t> SobelOperator(const Matrix<uint8_t> &src, const Matrix<int8_t> &Gx,
                              const Matrix<int8_t> &Gy, Matrix<uint8_t> &degrees) {
  degrees.assign(src.size(), vector<uint8_t>(src[0].size()));
  Matrix<uint8_t> out(src);
  for (int i = 1; i < src.size() - 1; ++i) {
    for (int j = 1; j < src[i].size() - 1; ++j) {
      Matrix<uint8_t> submatrix = GetSubmatrix(src, i - 1, j - 1, 3, 3);
      double x_pixel = SumMatrixElements(ElementWiseMultiplication(Gx, submatrix));
      double y_pixel = SumMatrixElements(ElementWiseMultiplication(Gy, submatrix));
      out[i][j] = static_cast<uint8_t>((sqrt(x_pixel * x_pixel + y_pixel * y_pixel)));
      int theta = static_cast<int>(atan(y_pixel / x_pixel) * 180 / M_PI);
      theta %= 180;
      if (theta <= 22.5 || theta > 157.5)
        degrees[i][j] = 0;
      else if (theta <= 67.5)
        degrees[i][j] = 1;
      else if (theta <= 112.5)
        degrees[i][j] = 2;
      else if (theta <= 157.5)
        degrees[i][j] = 3;
    }
  }
  return move(out);
}

template<typename T>
bool Check(const Matrix<T> &m, int i, int j, int val_i, int val_j) {
  if (j < 0 || j > m[0].size() - 1 || i < 0 || i > m.size() - 1 ||
      val_j < 0 || val_j > m[0].size() - 1 || val_i < 0 || val_i > m.size() - 1)
    return false;
  else {
    return m[i][j] <= m[val_i][val_j];
  }
}

template<typename T>
Matrix<T> NonMaximumSuppression(const Matrix<T> &src, const Matrix<uint8_t> &degrees) {
  Matrix<T> out(src);
  for (int i = 0; i < out.size(); i++) {
    for (int j = 0; j < out[i].size(); j++) {
      switch (degrees[i][j]) {
        case 0:
          if (Check(out, i, j, i, j - 1) || Check(out, i, j, i, j + 1))
            out[i][j] = 0;
          break;
        case 1:
          if (Check(out, i, j, i - 1, j - 1) || Check(out, i, j, i + 1, j + 1))
            out[i][j] = 0;
          break;
        case 2:
          if (Check(out, i, j, i - 1, j) || Check(out, i, j, i + 1, j))
            out[i][j] = 0;
          break;
        case 3:
          if (Check(out, i, j, i - 1, j + 1) || Check(out, i, j, i + 1, j - 1))
            out[i][j] = 0;
          break;
      }
    }
  }
  return out;
}

template<typename T>
Matrix<T> DoubleThresholding(const Matrix<T> &src, double lower_border, double upped_border) {
  Matrix<T> out(src.size(), vector<T>(src[0].size()));
  int down = (int) (lower_border * 255);
  int up = (int) (upped_border * 255);
  for (int i = 0; i < src.size(); ++i) {
    for (int j = 0; j < src[i].size(); ++j) {
      if (src[i][j] >= up)
        out[i][j] = 255;
      else {
        if (src[i][j] <= down)
          out[i][j] = 0;
        else {
          out[i][j] = 127;
        }
      }
    }
  }
  return out;
}

template<typename T>
Matrix<T> BlobExtraction(const Matrix<T> &src, vector<Point> &white_points) {
  Matrix<T> out(src);
  Matrix<int8_t> moves({{-1, -1, -1, 0, 0, 1, 1, 1}, {-1, 0, 1, -1, 1, -1, 0, 1}});
  int shift = 5;
  for (int i = 0; i < src.size(); ++i) {
    for (int j = 0; j < src[i].size(); ++j) {
      if (out[i][j] == 255 && i > shift && i < src.size() - shift && j > shift && j < src[0].size() - shift)
        white_points.push_back({i, j});
      if (out[i][j] == 127) {
        int up_count = 0;
        for (int k = 0; k < moves[0].size(); ++k) {
          if (i + moves[0][k] < 0 || i + moves[0][k] > out.size() - 1 ||
              j + moves[1][k] < 0 || j + moves[1][k] > out[0].size()) {
            out[i][j] = 0;
            continue;
          }
          if (out[i + moves[0][k]][j + moves[1][k]] == 255) {
            ++up_count;
            break;
          }
        }
        if (up_count > 0 && i > shift && i < src.size() - shift && j > shift && j < src[0].size() - shift) {
          out[i][j] = 255;
          white_points.push_back({i, j});
        } else
          out[i][j] = 0;
      }
    }
  }
  return out;
}

Matrix<uint8_t> YUVfromRGB(const Matrix<Point3> &matrix) {
  Matrix<uint8_t> m(matrix.size(), vector<uint8_t>(matrix[0].size()));
  for (int i = 0; i < matrix.size(); i++) {
    for (int j = 0; j < matrix[i].size(); j++) {
      m[i][j] = 0.257 * matrix[i][j].red + 0.504 * matrix[i][j].green + 0.098 * matrix[i][j].blue + 16;
    }
  }
  return m;
}

Matrix<Point> Clustering(const vector<Point> &white_points) {
  vector<int> pixel_to_clusters(white_points.size(), -1);
  pixel_to_clusters[0] = 0;
  Matrix<Point> clusters;
  clusters.push_back({});
  bool new_cluster = true;
  int space = 30;
  for (int i = 0; i < white_points.size(); ++i) {
    if (pixel_to_clusters[i] == -1) {
      clusters.push_back({white_points[i]});
      pixel_to_clusters[i] = clusters.size() - 1;
    }
    for (int j = 0; j < white_points.size(); ++j) {
      if (pixel_to_clusters[j] == -1) {
        if (static_cast<int>(sqrt(
            (white_points[i].x - white_points[j].x) * (white_points[i].x - white_points[j].x) +
                (white_points[i].y - white_points[j].y) * (white_points[i].y - white_points[j].y))) < space) {
          new_cluster = false;
          pixel_to_clusters[j] = pixel_to_clusters[i];
          clusters[pixel_to_clusters[j]].push_back(white_points[j]);
        }
      }
    }
    if (new_cluster) {
      clusters.push_back({white_points[i]});
      pixel_to_clusters[i] = clusters.size() - 1;
    }
  }
  return clusters;
}

vector<Rect> FindingEdges(const Matrix<Point> &clusters) {
  vector<Rect> rects;
  for (const auto &cluster : clusters) {
    if (cluster.size() < 50)
      continue;
    int left = INT_MAX, right = 0, down = INT_MAX, up = 0;
    for (const auto &el : cluster) {
      if (el.x < down)
        down = el.x;
      if (el.x > up)
        up = el.x;
      if (el.y < left)
        left = el.y;
      if (el.y > right)
        right = el.y;
    }
    rects.push_back({{down, left}, {up, right}});
  }
  return rects;
}

void DrawRects(const vector<Rect> &rects, Matrix<Point3> &color_matrix) {
  for (const auto &rect : rects) {
    for (int i = rect.down_left.y; i < rect.up_right.y; i++) {
      color_matrix[rect.down_left.x][i] = {0, 0, 255};
      color_matrix[rect.up_right.x][i] = {0, 0, 255};
    }
    for (int i = rect.down_left.x; i < rect.up_right.x; i++) {
      color_matrix[i][rect.down_left.y] = {0, 0, 255};
      color_matrix[i][rect.up_right.y] = {0, 0, 255};
    }
  }
}

int main() {
  Manager manager("Detector");
  manager.read("DroneImage.jpg");
  Matrix<Point3> matrix(manager.getMatrix());
  Matrix<uint8_t> grey_matrix(YUVfromRGB(matrix));
  // manager.show(grey_matrix, "Grey");
  const int size = 7;
  const double sigma = 1;
  Matrix<uint8_t> after_gauss_matrix(GaussFilter(grey_matrix, size, sigma));
  //manager.show(after_gauss_matrix, "AfterGauss");

  Matrix<int8_t> Gx({{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}});
  Matrix<int8_t> Gy({{-1, -2, -1}, {0, 0, 0}, {1, 2, 1}});

  Matrix<uint8_t> degrees;
  Matrix<uint8_t> after_sobel(SobelOperator(after_gauss_matrix, Gx, Gy, degrees));
  // manager.show(after_sobel, "AfterSobel");

  Matrix<uint8_t> after_NonMaximumSuppression(NonMaximumSuppression(after_sobel, degrees));
  // manager.show(after_NonMaximumSuppression, "AfterNonMaximumSuppression");

  const double low_border = 0.3, upper_border = 0.45;
  Matrix<uint8_t> after_thresholding(DoubleThresholding(after_NonMaximumSuppression, low_border, upper_border));
  //manager.show(after_thresholding, "AfterThresholding");

  vector<Point> white_points;
  Matrix<uint8_t> after_extraction(BlobExtraction(after_thresholding, white_points));
 // manager.show(after_extraction, "After Extraction");

  Matrix<Point> clusters(Clustering(white_points));
  vector<Rect> rects(FindingEdges(clusters));

  Matrix<Point3> color_matrix(manager.getMatrix());
  DrawRects(rects, color_matrix);

  manager.showRes(color_matrix);
  return 0;
}