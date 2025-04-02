from ctypes import CDLL, Structure, POINTER, c_double, c_int, c_bool, c_char_p
import numpy as np
import os
import atexit

class CoordArray(Structure):
    _fields_ = [
        ("data", POINTER(c_double)),
        ("size", c_int)
    ]

# Load library
lib_path = os.path.join(os.path.dirname(__file__), '../..', 'build/lib/molecule_interface.dylib')
lib = CDLL(lib_path)

# Set up function signatures
lib.init_state.argtypes = [c_char_p, c_char_p]
lib.init_state.restype = c_bool

lib.try_update.argtypes = [c_int]
lib.try_update.restype = CoordArray

lib.accept_update.restype = c_bool

lib.get_current_coords.restype = CoordArray

lib.free_coords.argtypes = [POINTER(c_double)]

# Register cleanup
atexit.register(lib.cleanup_state)

def init_molecule(seq_file: str, coord_file: str) -> bool:
    """Initialize the molecule with sequence and coordinate files."""
    return lib.init_state(
        seq_file.encode('utf-8'),
        coord_file.encode('utf-8')
    )

def try_update(change_idx: int) -> np.ndarray | None:
    """Try updating the molecule structure and return coordinates if successful."""
    result = lib.try_update(change_idx)
    
    if result.size == 0:
        # print('result size: 'result.size)
        return None
        
    # Convert to numpy array and reshape to Nx3
    coords = np.ctypeslib.as_array(result.data, shape=(result.size // 3, 3))
    
    # Make a copy of the data so we can free the original
    coords_copy = coords.copy()
    
    # Free the original data
    lib.free_coords(result.data)
    
    return coords_copy

def accept_update() -> bool:
    """Accept the last successful update."""
    return lib.accept_update()

def get_current_coords() -> np.ndarray | None:
    """Get coordinates of current molecule state."""
    result = lib.get_current_coords()
    
    if result.size == 0:
        return None
        
    # Convert to numpy array and reshape to Nx3
    coords = np.ctypeslib.as_array(result.data, shape=(result.size // 3, 3))
    
    # Make a copy of the data so we can free the original
    coords_copy = coords.copy()
    
    # Free the original data
    lib.free_coords(result.data)
    
    return coords_copy

def clear_state():
    atexit.register(lib.cleanup_state)