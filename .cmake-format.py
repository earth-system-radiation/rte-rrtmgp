with section("parse"):
    additional_commands = {"check_python3_package": {"pargs": 1, "kwargs": {"CODE": 1}}}

with section("format"):
    dangle_parens = True
    max_lines_hwrap = 0
    keyword_case = "upper"
    autosort = True

with section("lint"):
    # The formatter sometimes fails to fit the code into the line limit (C0301) and can
    # disagree with the linter regarding the indentation (C0307):
    disabled_codes = ["C0301", "C0307"]
    # Names of local variables must be in lowercase but sometimes we need to
    # override standard CMake variables:
    local_var_pattern = "CMAKE_[0-9A-Z_]+|[a-z][0-9a-z_]+"
    # The standard names of the languages in CMake are C and Fortran. Names of
    # private variables must be in lowercase but can have substings "C" and
    # "Fortran":
    private_var_pattern = (
        "([a-z_][0-9a-z_]*_)?(C|Fortran)(_[a-z_][0-9a-z_]*)?|[a-z_][0-9a-z_]+"
    )
    # The standard name of the language in CMake is Fortran. Names of public
    # variables must be in uppercase but can have substring "Fortran":
    public_var_pattern = "([A-Z][0-9A-Z_]*_)?Fortran(_[A-Z][0-9A-Z_]*)?|[A-Z][0-9A-Z_]+"
