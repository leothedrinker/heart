#define _USE_MATH_DEFINES
#include <string.h>
#include <array>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <chrono>
#include <thread>
#include <string>

constexpr std::size_t dim_x = 80;
constexpr std::size_t dim_y = 22;

template<std::size_t N>
struct StringLiteral {

    static constexpr std::size_t size = N - 1;

    // thanks to https://ctrpeach.io/posts/cpp20-string-literal-template-parameters/
    constexpr StringLiteral(const char (&str)[N]) {
        std::copy_n(str, N - 1, value);
    }
    
    char value[N - 1];
};


template <std::size_t size_x_, std::size_t size_y_, StringLiteral lum_lit>
class Canvas {

public:
    static constexpr std::size_t dim = size_x_ * size_y_;
    static constexpr std::size_t size_x = size_x_;
    static constexpr std::size_t size_y = size_y_;
    // static constexpr StringLiteral

    Canvas(std::string_view title) : title_(title) {
        std::cout << "\x1b[2J";
        reset();
    }

    Canvas(const Canvas&) = delete;
    Canvas(Canvas&&) = default;
    Canvas& operator=(const Canvas&) = delete;
    Canvas& operator=(Canvas&&) = default;

    void reset() {
        lum_.fill(-1);
        zbuff_.fill(0);
    }

    void render() {
        std::cout << "\x1b[H";
        size_t l = size_x > title_.size() ? (size_x - title_.size()) / 2 : 0;
        std::cout << std::setfill('=') << std::setw(l) << ' ' << title_ << ' ' << std::setfill('=') << std::setw(l) << "\n";
        
        for (int idx = 0; idx <= dim; ++idx) {
            std::cout << (idx % size_x ? render_lum(lum_[idx]) : '\n');
        }
    }

    void draw_point(int x, int y, double z_inv, double s) {
        if (x >= 0 && x < size_x && y >= 0 && y < size_y) {
            int idx = x + size_x * y;
            if (z_inv > zbuff_[idx]) {
                zbuff_[idx] = z_inv;
                lum_[idx] = s;
            }
        }
    }

private:

    std::array<double, dim> lum_;
    std::array<double, dim> zbuff_;
    const std::string title_;

    inline char render_lum(double s) {
        if (s >= 0) {
            return lum_lit.value[static_cast<int>(s * lum_lit.size)];
        } else {
            return ' ';
        }
    }
    
};

template<std::size_t t, std::size_t p>
struct Heart {
    static constexpr std::size_t theta_n = t, phi_n = p;
    double scale;
    double z; // distance to viewer
};

struct Rotation {
    double w1 = 0, w2 = 0, w3 = 0; 
    mutable double cos_w1, cos_w2, cos_w3;
    mutable double sin_w1, sin_w2, sin_w3;

    void compute_ratio() const {
      cos_w1 = cos(w1);
      cos_w2 = cos(w2);
      cos_w3 = cos(w3);
      sin_w1 = sin(w1);
      sin_w2 = sin(w2);
      sin_w3 = sin(w3);
    }

    inline std::tuple<double, double, double> rotate(const std::tuple<double, double, double>& vec) const {
      const auto &[x, y, z] = vec;
      return {
        x * cos_w2 * cos_w3 + y * cos_w2 * sin_w3 + z * sin_w2,
        x * sin_w1 * sin_w2 * (-cos_w3 - cos_w1 * sin_w3) + y * (cos_w1 * cos_w3 - sin_w1 * sin_w2 * sin_w3) + z * sin_w1 * cos_w2,
        x * (sin_w1 * sin_w3 - cos_w1 * sin_w2 * cos_w3) + y * (-sin_w1 * cos_w3 - cos_w1 * sin_w2 * sin_w3) + z * cos_w1 * cos_w2
      };
    }
};

using SurfaceT = std::tuple<double, double, double>;  // <x, y, z>
using NormalT = std::tuple<double, double, double>;  // <Nx, Ny, Nz>

int main() {

    // canvas set up
    Canvas<160, 44, R"(.,-~:;=!*#$@)"> canvas("Merry Christmas");


    // heart discretization
    constexpr std::size_t heart_theta_n = 300;
    constexpr std::size_t heart_phi_n = 100;
    constexpr std::size_t heart_n = heart_theta_n * heart_phi_n;
    constexpr double d_theta{2 * M_PI / heart_theta_n}, d_phi{M_PI / heart_phi_n};
    
    Rotation rot{};
    
    // precompute all heart discretizations
    std::array<std::pair<SurfaceT, NormalT>, heart_n> heart;
    
    for (int ti = 0; ti < heart_theta_n; ++ti) {
      double t = d_theta * ti;
      double ct{cos(t)}, c2t{cos(2*t)}, c3t{cos(3*t)};
      double st{sin(t)}, s2t{sin(2*t)}, s3t{sin(3*t)};

      for (int pi = 0; pi < heart_phi_n; ++pi) {
        double p = d_phi * pi;
        double cp{cos(p)}, sp{sin(p)};

        int idx = ti * heart_phi_n + pi;

        double theta_x{1.875 * st - 0.5 * s3t}; 
        double theta_y{1.875 * ct - 0.5 * c2t - 0.25 * c3t};
        double theta_xp{1.875 * ct - 1.5 * c3t};  // theta_x's derivative wrt to theta
        double theta_yp{s2t + 0.75 * s3t - 1.875 * st};  // theta_y's derivative wrt theta

        heart[idx] = {
          {sp * theta_x, sp * theta_y, cp},
          {sp * sp * theta_yp, sp * sp * theta_xp, sp * cp * (theta_x * theta_yp - theta_y * theta_xp)}
        };
      }
    }
    
    for(;;) {

        canvas.reset();

        rot.compute_ratio();

        for (auto& [point, norm]: heart) {
          auto& [x, y, z] = point;
          auto& [nx, ny, nz] = norm;

          // rotation
          auto [xr, yr, zr] = rot.rotate(point);
          auto [nxr, nyr, nzr] = rot.rotate(norm);

          // projection onto screen
          double d = 1 / (zr + 5);
          int xt = canvas.size_x / 2 + canvas.size_x / 2 * d * xr;
          int yt = canvas.size_y / 2 - canvas.size_y * d * yr;

          // finding luminance
          double lum = std::max(0., (nyr + nzr - nxr) / sqrt(nxr * nxr + nyr * nyr + nzr * nzr) / sqrt(3));

          canvas.draw_point(xt, yt, d, lum);
        }

        canvas.render();
        
        rot.w2 += 0.05;
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
}