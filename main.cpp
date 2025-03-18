#include <iostream>
#include <opencv2/opencv.hpp>
#include "MultifractalFormalism.h"
#include "FasterBoxCounter.h"

int main() {
    // Путь к изображению
    std::string imagePath = "/Users/dudberoll/CLionProjects/cw/Texture Image Sample.jpg"; // Замените на путь к вашему изображению

    // Чтение изображения в градациях серого
    cv::Mat image = cv::imread(imagePath, cv::IMREAD_GRAYSCALE);
    if (image.empty()) {
        std::cerr << "Ошибка: не удалось загрузить изображение!" << std::endl;
        return -1;
    }

    std::cout << "Изображение успешно загружено: "
              << image.cols << "x" << image.rows << std::endl;

    // Конвертация изображения в std::vector для обработки
    std::vector<std::vector<uint8_t>> imageVector(image.rows, std::vector<uint8_t>(image.cols));
    for (int i = 0; i < image.rows; ++i) {
        for (int j = 0; j < image.cols; ++j) {
            imageVector[i][j] = image.at<uint8_t>(i, j);
        }
    }

    // Вызов функции расчёта гельдеровских экспонент
    try {
        auto result = MultifractalFormalism::CountHeroldExpIso(
                imageVector, 20, [](float progress) {
                    std::cout << "Прогресс: " << progress * 100 << "%" << std::endl;
                });

        // Вывод результата (например, только первые 5x5 элементов)
        std::cout << "Результат (первая часть):" << std::endl;
        for (int i = 0; i < std::min(5, static_cast<int>(result.size())); ++i) {
            for (int j = 0; j < std::min(5, static_cast<int>(result[i].size())); ++j) {
                std::cout << result[i][j] << " ";
            }
            std::cout << std::endl;
        }
        
        // Пример использования FasterBoxCounter для расчета локальных Box размерностей
        // Создаем бинарное изображение с порогом 128
        std::vector<std::vector<bool>> binaryImage(image.rows, std::vector<bool>(image.cols));
        for (int i = 0; i < image.rows; ++i) {
            for (int j = 0; j < image.cols; ++j) {
                binaryImage[i][j] = (image.at<uint8_t>(i, j) > 128);
            }
        }
        
        // Инициализируем FasterBoxCounter с размером окна 2^5 = 32
        FasterBoxCounter boxCounter(5);
        boxCounter.Prepare(binaryImage);
        
        // Рассчитываем фрактальные размерности для нескольких точек
        std::cout << "\nЛокальные Box размерности (примеры):" << std::endl;
        for (int i = 0; i < std::min(5, image.rows); i += image.rows/5) {
            for (int j = 0; j < std::min(5, image.cols); j += image.cols/5) {
                double fractalDim = boxCounter.GetFractalDim(i, j);
                std::cout << "Точка (" << i << "," << j << "): " << fractalDim << std::endl;
            }
        }
        
    } catch (const std::exception& ex) {
        std::cerr << "Ошибка при расчёте: " << ex.what() << std::endl;
        return -1;
    }

    return 0;
}
