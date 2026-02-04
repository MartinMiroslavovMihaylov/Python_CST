# CST_Constructor/__init__.py

from importlib.metadata import PackageNotFoundError, version

# --- Version ---
try:
    __version__ = version(__name__)
except PackageNotFoundError:
    __version__ = "0+unknown"

# --- Import classes ---
from .CST_Constructor import CST_Commands, Curves, PostProcess

__all__ = ["__version__", "CST_Commands", "Curves", "PostProcess"]