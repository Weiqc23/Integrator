cmake_minimum_required(VERSION 3.10)
project(JFNK)
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()
# 添加调试信息
set(CMAKE_BUILD_TYPE Debug)
# 设置 C++ 编译器标志，启用 C++11 标准
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# 检查编译器是否为 GNU C 或 GNU C++ 编译器
if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  # 如果是 GNU 编译器，更新 C++ 编译器标志
  # -Wall 开启所有常见的编译警告
  # -Wno-long-long 不显示关于 long long 类型的警告
  # -pedantic 严格遵循标准，显示更多的标准兼容性警告
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic")
endif()

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# 设置策略 CMP0115
if(POLICY CMP0115)
  cmake_policy(SET CMP0115 NEW)
endif()

# 指定vcpkg工具链文件路径
set(CMAKE_TOOLCHAIN_FILE "D:/vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE STRING "Vcpkg toolchain file")

find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED COMPONENTS)
find_package(Python3 COMPONENTS Development NumPy REQUIRED)

# target_include_directories(test PRIVATE ${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS})
include_directories(${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS})
include_directories(${EIGEN3_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

set(SOURCE_FILES   
    # integrator.cpp
    integrator_liner.cpp
    points.cpp
    preprocess.cpp
)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -DNDEBUG")
add_executable(test main_test.cpp ${SOURCE_FILES})

target_link_libraries(test PRIVATE Eigen3::Eigen)
target_link_libraries(test PRIVATE Eigen3::Eigen ${Boost_LIBRARIES})
target_link_libraries(test PRIVATE ${Python3_LIBRARIES})