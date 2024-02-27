#
# Top-level Makefile
#
.PHONY: libs tests check docs

libs:
	$(MAKE) -C build
	$(MAKE) -C tests
	$(MAKE) -C examples/all-sky
	$(MAKE) -C examples/rfmip-clear-sky

tests:
	$(MAKE) -C examples/rfmip-clear-sky $@
	$(MAKE) -C examples/all-sky         $@
	$(MAKE) -C tests                    $@

check:
	$(MAKE) -C examples/rfmip-clear-sky $@
	$(MAKE) -C examples/all-sky         $@
	$(MAKE) -C tests                    $@

docs:
	@cd doc; ./build_documentation.sh

clean:
	$(MAKE) -C build                    $@
	$(MAKE) -C examples/rfmip-clear-sky $@
	$(MAKE) -C examples/all-sky         $@
	$(MAKE) -C tests                    $@
	rm -rf public
