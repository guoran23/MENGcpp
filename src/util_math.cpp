#include "util_math.h"

void UtilMath::lgwt(int n, double a, double b, std::vector<double> &x1d,
                    std::vector<double> &w1d) {
  if (n <= 1) {
    throw std::invalid_argument("lgwt Error: n must be greater than 1");
  }
  if (static_cast<int>(x1d.size()) < n) {
    throw std::invalid_argument("lgwt Error: x1d must have a size larger than n.");
  }
  if (static_cast<int>(w1d.size()) < n) {
    throw std::invalid_argument("lgwt Error: w1d must have a size larger than n.");
  }
  const double eps = std::numeric_limits<double>::epsilon();

  int n0 = n - 1;
  int n1 = n0 + 1;
  int n2 = n0 + 2;

  // Allocate memory for dynamic arrays
  std::vector<double> xu(n1);
  std::vector<double> y(n1), Lp(n1);
  // n1 X n2, initial L=0.0
  std::vector<std::vector<double>> L(n1, std::vector<double>(n2, 0.0));
  std::vector<double> y0(n1, 2.0);

  // Initial guess
  for (int i = 0; i < n1; ++i) {
    xu[i] = -1.0 + 2.0 * i / (n1 - 1.0);
    y[i] = std::cos((2.0 * i + 1.0) * M_PI / (2.0 * n0 + 2.0)) +
           (0.27 / n1) * std::sin(M_PI * xu[i] * n0 / n2);
  }
  //  Legendre-Gauss Vandermonde Matrix, L(n1,n2) and Lp(n1)
  //   Compute the zeros of the N+1 Legendre Polynomial
  // using the recursion relation and the Newton-Raphson method
  //  Iterate until new points are uniformly within epsilon of old points
  while (true) {
    // Compute the maximum absolute difference between y and y0
    double max_diff = 0.0;
    for (int i = 0; i < n1; ++i) {
      max_diff = std::max(max_diff, std::abs(y[i] - y0[i]));
    }

    // Exit the loop if the difference is smaller than epsilon
    if (max_diff <= eps) {
      break;
    }

    // Initialize Legendre-Gauss Vandermonde Matrix
    for (int i = 0; i < n1; ++i) {
      L[i][0] = 1.0;
      L[i][1] = y[i];
    }

    for (int k = 1; k < n1; ++k) {
      for (int i = 0; i < n1; ++i) {
        L[i][k + 1] =
            ((2.0 * k + 1.0) * y[i] * L[i][k] - k * L[i][k - 1]) / (k + 1);
      }
    }

    // Derivative of Legendre polynomial
    for (int i = 0; i < n1; ++i) {
      Lp[i] = n2 * (L[i][n1 - 1] - y[i] * L[i][n2 - 1]) / (1.0 - y[i] * y[i]);
    }

    // Update y0 and y
    y0 = y;
    for (int i = 0; i < n1; ++i) {
      y[i] = y0[i] - L[i][n2 - 1] / Lp[i];
    }
  }

  // Compute x1d and w1d
  double n2_o_n1 = static_cast<double>(n2) / n1;
  for (int i = 0; i < n; ++i) {
    x1d[i] = (a * (1.0 - y[i]) + b * (1.0 + y[i])) / 2.0;
    w1d[i] =
        (b - a) / ((1.0 - y[i] * y[i]) * Lp[i] * Lp[i]) * n2_o_n1 * n2_o_n1;
  }

  // Reverse x1d and w1d, then x1d is in ascending order
  std::reverse(x1d.begin(), x1d.end());
  std::reverse(w1d.begin(), w1d.end());
}
