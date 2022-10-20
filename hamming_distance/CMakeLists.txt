cmake_minimum_required(VERSION 3.6)

project(hamming CXX)

find_package(OpenMP)

add_executable (hamming new_hamming.cpp)
add_executable (fragment_hamming fragment_hamming.cpp)
add_executable (fragment_tree_hamming fragment_tree_hamming.cpp)
add_executable (fragment_tree_hamming_new fragment_tree_hamming_new.cpp)

target_link_libraries(fragment_hamming OpenMP::OpenMP_CXX)
target_link_libraries(hamming OpenMP::OpenMP_CXX)
target_link_libraries(fragment_tree_hamming OpenMP::OpenMP_CXX)
target_link_libraries(fragment_tree_hamming_new OpenMP::OpenMP_CXX)
