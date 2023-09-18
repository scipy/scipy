# This spec file is read by gcc when linking.  It is used to specify the
# standard libraries we need in order to link with various sanitizer libs.

*link_libasan: -lrt -ldl -lrt -lpthread -lm

*link_libhwasan: -ldl -lrt -lpthread -lm

*link_libtsan: -lrt -ldl -lrt -lpthread -lm

*link_libubsan: -ldl -lrt -lpthread -lm

*link_liblsan: -ldl -lrt -lpthread -lm

