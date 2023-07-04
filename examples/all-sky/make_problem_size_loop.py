#! /usr/bin/env python
#
# This script...
#
import argparse
import csv 
import os

# template: exe ncol nlay nloops output_file k_dist clouds aeorsols 

# specify: kdist, optional clouds, aerosols. specify nloops 
# Toggle clouds and aerosols?  
# Loop over sets of ncol, nlay, 
# output name 


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Description here ")
    # Argument parseing described at 
    # https://stackoverflow.com/questions/15753701/how-can-i-pass-a-list-as-a-command-line-argument-with-argparse
    parser.add_argument("-x", "--executable", 
                        type=str, 
                        default="./rrtmgp_allsky",
                        help="Path to exectuable")
    parser.add_argument("-c", "--ncol", 
                        type=lambda items: [int(i) for i in list(csv.reader([items]))[0]], 
                        default="2,4,8,16",
                        help="Number of columns  (multiple e.g. 2,4,8,16)")
    parser.add_argument("-l", "--nlay", 
                        type=lambda items: [int(i) for i in list(csv.reader([items]))[0]], 
                        default="32, 64, 96",
                        help="Number of layers (multiple e.g. 32,64.96)")
    parser.add_argument("-i", "--nloops", 
                        type=int, default=1,
                        help="Number of loops (same for all ncol)")
    parser.add_argument("-o", "--output_file", 
                        type=str, 
                        default="rrtmgp-allsky-output.nc",
                        help="Path to output file")
    parser.add_argument("-k", "--k_distribution", 
                        type=str, 
                        required = True, 
                        help="Path to gas optics file [required]")
    parser.add_argument("-cl", "--cloud_optics", 
                        type=str, default="",
                        help="Path to cloud optics file")
    parser.add_argument("-a", "--aerosol_optics", 
                        type=str, default="",
                        help="Path to aerosol optics file")
    args = parser.parse_args()

    # Can't supply aerosols without clouds 
    if(args.cloud_optics == "" and args.aerosol_optics != ""): 
        raise AssertionError("Need to supply cloud optics if providing aerosol optics")

    # Every combo of ncol, nlay 
    for l in args.nlay: 
        for i in args.ncol: 
          print(f"{args.executable} {i:6d} {l:4d} {args.nloops:3d} " + \
                f"{args.output_file} {args.k_distribution}"         + \
                f"{args.cloud_optics} {args.aerosol_optics}") 
