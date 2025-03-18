#ifndef MULTIFRACTALFORMALISM_H
#define MULTIFRACTALFORMALISM_H

#include <vector>
#include <functional>
#include <cstdint>

class MultifractalFormalism {
public:
    // Основной метод для расчёта гельдеровских экспонент
    static std::vector<std::vector<float>> CountHeroldExpIso(
            const std::vector<std::vector<uint8_t>>& im,
            int rmax = 20,
            std::function<void(float)> reportProgress = nullptr
    );

private:
    // Вспомогательные методы
    static int getSumHelper(
            const std::vector<std::vector<int>>& ar,
            int lx, int rx, int ly, int ry
    );

    static int getSum(
            const std::vector<std::vector<int>>& ar,
            int lx, int rx, int ly, int ry
    );

    // Метод для расчёта коэффициента наклона линейной регрессии
    static float calculateSlope(const std::vector<double>& x, const std::vector<double>& y);
};

#endif // MULTIFRACTALFORMALISM_H
#ifndef MULTIFRACTALFORMALISM_H
#define MULTIFRACTALFORMALISM_H

#include <vector>
#include <functional>
#include <cstdint>

class MultifractalFormalism {
public:
    // Основной метод для расчёта гельдеровских экспонент
    static std::vector<std::vector<float>> CountHeroldExpIso(
        const std::vector<std::vector<uint8_t>>& im,
        int rmax = 20,
        std::function<void(float)> reportProgress = nullptr
    );

private:
    // Вспомогательные методы
    static int getSumHelper(
        const std::vector<std::vector<int>>& ar,
        int lx, int rx, int ly, int ry
    );

    static int getSum(
        const std::vector<std::vector<int>>& ar,
        int lx, int rx, int ly, int ry
    );

    // Метод для расчёта коэффициента наклона линейной регрессии
    static float calculateSlope(const std::vector<double>& x, const std::vector<double>& y);
};

#endif // MULTIFRACTALFORMALISM_H
