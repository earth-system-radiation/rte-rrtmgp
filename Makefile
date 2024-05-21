#
# Top-level Makefile
#
.PHONY: libs tests check
all:    libs tests check

libs:
	$(MAKE) -C build                    $@

tests:
	$(MAKE) -C tests                    $@
	$(MAKE) -C examples/rfmip-clear-sky $@
	$(MAKE) -C examples/all-sky         $@

check:
	$(MAKE) -C tests                    $@
	$(MAKE) -C examples/rfmip-clear-sky $@
	$(MAKE) -C examples/all-sky         $@

docs:
	@cd doc; ./build_documentation.sh

clean:
	$(MAKE) -C build                    $@
	$(MAKE) -C tests                    $@
	$(MAKE) -C examples/rfmip-clear-sky $@
	$(MAKE) -C examples/all-sky         $@
	rm -rf public
