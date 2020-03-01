# Copyright (c) 2019. Lloyd T. Elliott

CFLAGS += -DBUILD_DATE='"$(shell date "+%B %Y")"' -DGIT_VERSION='"$(shell git describe --always --abbrev=10)"' $(COMMON) $(foreach i,$(INCLUDE),-I$(i))

all:
	gcc kgen.c -o kgen                        `# executable name        ` \
        $(CFLAGS)                                 `# system flags           ` \
	-std=c17 -O3                              `# optimization           ` \
	-DNDEBUG -DMKL_ILP64                      `# optimization           ` \
	-fno-trapping-math                        `# optimization           ` \
	-funsafe-math-optimizations               `# optimization           ` \
	-fno-rounding-math                        `# optimization           ` \
	-fcx-limited-range                        `# optimization           ` \
	-fno-signed-zeros                         `# optimization           ` \
	-floop-nest-optimize                      `# optimization           ` \
	-static-libgcc -static-libstdc++ -static  `# static build           ` \
	-flto -fuse-linker-plugin                 `# static build           ` \
	-march=native -mtune=generic -m64         `# target architecture    ` \
	-I${MKLROOT}/include                      `# MKL specific flags     ` \
	-L${MKLROOT}/lib/intel64                  `# MKL specific flags     ` \
	-lmkl_intel_ilp64 -lmkl_sequential        `# MKL specific flags     ` \
	-lmkl_core -lpthread                      `# MKL specific flags     ` \
	-Wl,--no-as-needed                        `# MKL specific flags     ` \
	-L. -lm -ldl -lqrupdate -lgsl             `# required libraries     `
