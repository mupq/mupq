HOST_SYMCRYPTO_SRC = \
	$(SYMCRYPTO_SRC_COMMON) \
	mupq/common/keccakf1600.c \
	mupq/pqclean/common/aes.c \
	mupq/pqclean/common/sha2.c

obj-host/libsymcrypto.a: $(call hostobjs,$(HOST_SYMCRYPTO_SRC))

HOST_LDLIBS += -lsymcrypto
HOST_LIBDEPS += obj-host/libsymcrypto.a
