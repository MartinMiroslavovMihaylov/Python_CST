# Python_CST/__init__.py

from importlib.metadata import PackageNotFoundError, version

# --- Version ---
try:
    __version__ = version(__name__)
except PackageNotFoundError:
    __version__ = "0+unknown"  # fallback for source installs

# --- Import classes/functions from Constructor module ---
from CST_Constructor.CST_Constructor import (
    CST_Commands,
    Curves,
)

# --- Exported names ---
__all__ = [
    "__version__",
    "CST_Commands",
    "Curves"
]