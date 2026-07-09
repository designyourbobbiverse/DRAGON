'''
GAZE1D.py
DRAGONGAZE
Created by Bobbie Markwick on 07/07/2026.
'''

from Config import *
import os


def h5FileName(n):
    if n < 0: return hdf_dir + "/" + h5_base_filename
    return h5FileName(-1) + "_" + str(n).zfill(5) + ".h5"
    
def imgFileName(n):
    return img_dir + "/" + img_base_filename + "_" + str(n).zfill(5) + ".png"

def imgFileName(n, key):
    return img_dir + "/" + fileprefix[key] + "_" + img_base_filename + "_" + str(n).zfill(5) + ".png"

def imgFileName(n, key, axis):
    return img_dir + "/" + fileprefix[key] + "_" + axis + "_" + img_base_filename + "_" + str(n).zfill(5) + ".png"

def fileExists(n):
    return os.path.exists(h5FileName(n))
