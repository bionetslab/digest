#!/usr/bin/python3

import sys
import os
import inspect

# =============================================================================
# Load all paths of required folders with scripts
# =============================================================================
UTILS_DIR = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
BASE_DIR = os.path.dirname(UTILS_DIR)
sys.path.insert(0, BASE_DIR+'/setup')
sys.path.insert(0, BASE_DIR+'/d_utils')
sys.path.insert(0, BASE_DIR+'/mappers')

# =============================================================================
# Set directories
# ============================================================================
FILES_DIR = BASE_DIR+'mapping_files/'

