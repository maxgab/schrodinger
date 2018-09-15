#!/usr/bin/env python3
"""loads custom parameters"""
import os
import sys
import numpy as np
def file_input(*files):
    """checks for input files

    Args:
        files (str): input files
    Returns:
        list: input files"""
    return_files = []
    for item in files:
        while os.path.isfile(item) == 0:
            item = input("Cannot find file " + item + ". Please input file location: ")
        return_files.append(item)
    if len(return_files) == 1:
        return_files = str(return_files[0])
    return return_files

def set_options(*params):
    """set custom options if -o called for float

    Args:
        params (list): parameters and default values [description, default_parameter]
    Returns:
        list: retrieved parameters"""
    params = np.array(params)
    if len(sys.argv) == 2:
        if sys.argv[1] == "-o":
            for index, item in enumerate(params):
                custom_inp = input(str(item[0]) + " [default: " + item[1] + "]: ")
                if custom_inp != "":
                    params[index][1] = float(custom_inp)
    return list(map(float, params[:, 1]))
