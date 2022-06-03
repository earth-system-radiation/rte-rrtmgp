#
# Top-level Makefile
#
.PHONY: libs tests check docs
all: libs tests check docs

libs:
	make -C build -j
	make -C tests -j 1
	make -C examples/all-sky -j
	make -C examples/rfmip-clear-sky -j

tests:
	make -C examples/rfmip-clear-sky tests
	make -C examples/all-sky         tests
	make -C tests                    tests

check:
	make -C examples/rfmip-clear-sky check
	make -C examples/all-sky         check
	make -C tests                    check

docs:
	@cd doc; ./build_documentation.sh

clean:
	make -C build clean
	make -C examples/rfmip-clear-sky clean
	make -C examples/all-sky         clean
	make -C tests                    clean
	rm -rf public
