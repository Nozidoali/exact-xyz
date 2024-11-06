add_executable(prepare_dicke prepare_dicke.cpp)
target_link_libraries(prepare_dicke cxyz)

add_executable(backend backend.cpp)
target_link_libraries(backend cxyz)

add_executable(evaluate evaluate.cpp)
target_link_libraries(evaluate cxyz)