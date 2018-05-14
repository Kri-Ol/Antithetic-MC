// https://en.wikipedia.org/wiki/Antithetic_variates

#include <cmath>
#include <iostream>
#include <random>
#include <tuple>
#include <string>

#define func auto


static inline func to_integrate(const double x) -> double {
    return 1.0 / (1.0 + x);
}


static func normal_MC(const int           N,
                            std::mt19937& rng) -> std::tuple<double, double> {

    auto u01{ std::uniform_real_distribution<double>{0.0, 1.0 } };

    auto avg{ 0.0 };
    auto var{ 0.0 };

    for (int k = 0; k != N; ++k) {
        auto u = u01(rng);

        auto v = to_integrate(u);

        avg += v;
        var += v * v;
    }

    avg /= double(N);
    var /= double(N);

    var = std::max(0.0, var - avg * avg);
    return { avg, std::sqrt(var / (double(N - 1))) };
}

static func antithetic_MC(const int           N,
                                std::mt19937& rng) -> std::tuple<double, double> {

    auto u01{ std::uniform_real_distribution<double>{0.0, 1.0 } };

    auto avg{ 0.0 };
    auto var{ 0.0 };

    for (int k = 0; k != N; ++k) {
        auto u = u01(rng);

        auto v = 0.5*(to_integrate(u) + to_integrate(1.0 - u));

        avg += v;
        var += v * v;
    }

    avg /= double(N);
    var /= double(N);

    var = std::max(0.0, var - avg * avg);
    return { avg, std::sqrt(var / (double(N - 1))) };
}


func main() -> int {
    constexpr int N{ 3000 };

    auto rng{ std::mt19937{123456789UL} };

    auto [avg, stddev] = normal_MC(N, rng);
    std::cout << avg << "   " << stddev << "\n";

    std::tie(avg, stddev) = antithetic_MC(N / 2, rng);
    std::cout << avg << "   " << stddev << "\n";

    {
        double q;
        std::cin >> q;
        std::cout << q;
    }

    return 0;
}