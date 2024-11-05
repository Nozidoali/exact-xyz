add_executable(main main.cpp)
target_link_libraries(main cxyz)

add_executable(backend backend.cpp)
target_link_libraries(backend cxyz)

add_executable(evaluate evaluate.cpp)
target_link_libraries(evaluate cxyz)