import re

import os
import argparse

def fortran_to_c(fortran_code):
    subroutine_pattern = re.compile(
        r'subroutine\s+(\w+)\s*'
        r'\((.*?)\)\s*&?\s*'
        r'bind\(C,\s*name="(.*?)"\)',
        re.DOTALL | re.IGNORECASE
    )
    var_pattern = re.compile(
        r'^\s*((?:integer|logical\(\w+\)|real\(\w+\)))\s*,?\s*' # Match variable type
        r'(?:dimension\s*\(([^)]*)\)\s*,?\s*)?'                 # Match dimensions (optional)
        r'(?:target\s*,?\s*&?\s*)?'                             # Match target (optional)
        r'\s*intent\(\s*(\w+)\s*\)\s*::\s*'                     # Match intent
        r'([\w, ]+(?:\s*&\s*\n[\w, ]+)*)',                      # Capture variable names
        re.MULTILINE                                            # Allow multiline processing
    )

    type_map = {
        'integer': 'int',
        'real(wp)': 'Float',
        'logical(wl)': 'Bool'
    }

    match = subroutine_pattern.search(fortran_code)
    if not match:
        return "", -1

    subroutineName = match.group(1)
    funcParams     = match.group(2)
    cFunctionName  = match.group(3)

    funcParams    = funcParams.replace(' ', '')
    # Remove comments
    funcParams    = re.sub(r'&(![^&\n]*)\n', '', funcParams)
    funcParams    = funcParams.replace('\n', '')
    funcParams    = funcParams.replace('&', '')
    funcParams    = funcParams.split(',')
    matchedParams = {}

    matches = var_pattern.findall(fortran_code)

    for match in matches:
        fortran_type = match[0]
        dimensions   = match[1] if match[1] else None
        intent       = match[2]
        variables    = match[3]

        variables = variables.replace('&', '')
        variables = variables.replace(' ', '')
        variables = variables.replace('\n', '')
        variables = variables.split(',')

        for var in variables:
            type           = type_map.get(fortran_type, 'void*')
            qualifier      = 'const' if intent == 'in' else ''
            pointer_or_ref = '*' if dimensions else '&'

            matchedParams[var] = {
                'type': type,
                'qualifier': qualifier,
                'pointer_or_ref': pointer_or_ref,
                'dimensions': dimensions
            }

    args = ""

    for idx, param in enumerate(funcParams):
        matchedParam   = matchedParams.get(param, {})
        qualifier      = matchedParam.get("qualifier", None)
        type           = matchedParam.get("type", "void")
        pointer_or_ref = matchedParam.get("pointer_or_ref", "*")
        dimensions     = matchedParam.get("dimensions", None)

        isLast = idx == len(funcParams) - 1

        strQualifier  = f"{qualifier} " if qualifier else ""
        strDimensions = "\n" if isLast else ",\n"

        if dimensions:
            dimensions = dimensions.replace('&', '').replace('\n', '').replace(' ', '')
            strDimensions = f" // Dims: ({dimensions})\n" if isLast else f", // Dims: ({dimensions})\n"

        args += f"\t{strQualifier}{type}{pointer_or_ref} {param}{strDimensions}"

    return f"void {cFunctionName}(\n" + args + ");", 0

def extract_between(text, start_str, end_str):
    start = text.find(start_str)
    if start == -1:
        return None
    end = text.find(end_str, start + len(start_str))
    if end == -1:
        return None
    return text[start:end]

def extract_subroutine_names(fortran_code):
    subroutine_pattern = re.compile(r'end subroutine\s+(\w+)', re.IGNORECASE)
    return subroutine_pattern.findall(fortran_code)

def extract_cbinds(fortran_code):
    cbinds = []
    passed = 0

    subroutines = extract_subroutine_names(fortran_code)

    for subroutine_name in subroutines:
        subroutine_code = extract_between(fortran_code, f"subroutine {subroutine_name}", f"end subroutine {subroutine_name}")

        if subroutine_code:
            res = subroutine_code.find("bind(C, name=")

            # Skip subroutine if not binded
            if res == -1:
                continue
        else:
            print("Invalid subroutine declaration! exiting...")
            return cbinds

        print(f"Processing: {subroutine_name} ...", end=" ")
        cbind, err = fortran_to_c(subroutine_code)

        if err:
            cbind = f"// void {subroutine_name}(...); Fix me!"
            print(f"FAILED")
        else:
            passed += 1
            print(f"PASSED")

        cbinds.append(cbind)

    print(f"Finished parsing: {passed}/{len(cbinds)}")

    return cbinds

def main():
    parser = argparse.ArgumentParser(description="Fortran subroutine definition to C definition parser")

    parser.add_argument("-if", "--in-file", dest="input_files", type=str, nargs="+", required=True,
                        help="Path(s) to the Fortran file(s) to be parsed")
    parser.add_argument("-of", "--out-file", dest="output_file", type=str, required=True,
                        help="Path where the output C definition file will be saved")

    args = parser.parse_args()

    inputFiles = args.input_files
    outputFile = args.output_file

    for file in inputFiles:
        if not os.path.isfile(file):
            print(f"Error: Input file '{file}' does not exist.")
            return

    copyright ="""/* This code is part of RRTMGP

Contacts: Robert Pincus and Eli Mlawer
email:  rrtmgp@aer.com

Copyright 2024-
Trustees of Columbia University in the City of New York
All right reserved.

Use and duplication is permitted under the terms of the
    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause

This header files defines the C bindings for the kernels used in RRTMGP
Adapted from code written by Chiel van Heerwaarden at Wageningen University and Research
*/"""

    subroutines = []

    for inputFile in inputFiles:
        print(f"Parsing {inputFile} to {outputFile}:")

        with open(inputFile, "r") as f:
            fortran_code = f.read()

        currentSubroutines = extract_cbinds(fortran_code)
        subroutines.extend(currentSubroutines)

        print()

    with open(outputFile, "w+") as f:
        f.write(copyright)

        f.write("\n#pragma once")
        f.write("\n")
        f.write("\n#include \"rte_types.h\"")
        f.write("\n")
        f.write("\nextern \"C\" ")
        f.write("{\n")

        for idx, subroutine in enumerate(subroutines):
            # Add one more tab to every line
            f.write("\t")
            subroutine = subroutine.replace("\t", "\t\t").replace(");", "\t);")
            f.write(subroutine)

            if idx == len(subroutines) - 1:
                f.write("\n")
            else:
                f.write("\n\n")

        f.write("}\n")


if __name__ == "__main__":
    main()
