from ctypes import CDLL, Structure, POINTER, c_double, c_int
import numpy as np
import os

class CoordArray(Structure):
    _fields_ = [
        ("data", POINTER(c_double)),
        ("size", c_int)
    ]

# Load library
lib_path = os.path.join(os.path.dirname(__file__), 'build/lib/libminimal_test.dylib')
print(f"Loading library from: {lib_path}")
lib = CDLL(lib_path)

# Set up function signatures - without asterisks
lib.create_test_coords.restype = CoordArray
lib.free_coords.argtypes = [POINTER(c_double)]

try:
    # Get coordinates
    print("Calling create_test_coords...")
    result = lib.create_test_coords()
    print(f"Got array of size: {result.size}")

    # Convert to numpy array
    coords = np.ctypeslib.as_array(result.data, shape=(result.size // 3, 3))
    print("\nCoordinates:")
    print(coords)

    # Clean up
    print("\nCleaning up...")
    lib.free_coords(result.data)
    print("Done!")

except Exception as e:
    print(f"Error: {e}")