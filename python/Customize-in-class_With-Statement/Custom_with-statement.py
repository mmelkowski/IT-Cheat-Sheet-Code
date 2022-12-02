# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 12:17:16 2019

@author: mmelkowski
"""


class file():
    def __init__(self, fname):
        self.fname = fname

    def __enter__(self):
        # Opening of file
        self.fd = open(self.fname, "r")
        return self.fd

    def __exit__(self, type, value, traceback):
        # Exception handling here
        self.fd.close()


if __name__ == "__main__":
    obj_file = file("test.txt")
    with obj_file as f:
        print(f.readlines())
