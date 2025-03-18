# Минимальная версия CMake
cmake_minimum_required(VERSION 3.10)

# Название проекта
project(cw)

# Устанавливаем стандарт C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Указываем исходные файлы проекта
set(SOURCES
        main.cpp   # Главный файл программы
        MultifractalFormalism.cpp # Файл с классом
)

# Добавляем исполняемый файл
add_executable(${PROJECT_NAME} ${SOURCES})

# Найти OpenCV
find_package(OpenCV REQUIRED)

# Указать путь к заголовочным файлам OpenCV
include_directories(${OpenCV_INCLUDE_DIRS})

# Линковка с библиотеками OpenCV
target_link_libraries(${PROJECT_NAME} ${OpenCV_LIBS})

# Опционально включаем предупреждения компилятора
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -Wpedantic)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(${PROJECT_NAME} PRIVATE /W4)
endif()

# Установка директории для бинарного файла
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
