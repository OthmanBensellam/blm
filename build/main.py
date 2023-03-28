# ctypes_test.py
import ctypes
import pathlib

libname = pathlib.Path().absolute() / "libtools.so"
print("1")
c_lib = ctypes.CDLL(libname)
print("2")
#x, y = 6, 2.3
#answer = c_lib.cmult(x, y)
