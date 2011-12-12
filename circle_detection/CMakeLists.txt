add_library( circdet SHARED circdet.cpp )

install( TARGETS circdet LIBRARY DESTINATION lib )
install( FILES circdet.hpp DESTINATION include )