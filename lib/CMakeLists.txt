add_library(
    ${PACKAGE_NAME}
    SHARED
    prepare_general.cpp
    prepare_sparse.cpp
    prepare_structured.cpp
    qcircuit.cpp
    qgate.cpp
    qstate.cpp
    resyn.cpp
)
target_link_libraries(${PACKAGE_NAME} PUBLIC xyz_headers)

