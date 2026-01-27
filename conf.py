# docs-src/conf.py
# Configuration file for the Sphinx documentation builder.

import os
import sys
import shutil
import subprocess

# --- Determine repo root (CODE_DIR) ---

# First try environment variable set by GitHub Actions
CODE_DIR = os.environ.get("CODE_DIR")
if not CODE_DIR or not os.path.isdir(CODE_DIR):
    # Fallback: assume repo is one folder up from docs-src
    CODE_DIR = os.path.abspath("../repo")

# Add repo root to sys.path for autodoc
sys.path.insert(0, CODE_DIR)
print(f"[conf.py] Added CODE_DIR to sys.path: {CODE_DIR}")

# --- Expose code examples under 'examples-src' ---

def _mount_examples():
    src = os.path.join(CODE_DIR, "Examples")
    dst = os.path.abspath(os.path.join(os.path.dirname(__file__), "examples-src"))

    if not os.path.isdir(src):
        print(f"[conf.py] No code examples at {src}; skipping.", file=sys.stderr)
        return

    try:
        if os.path.exists(dst):
            if os.path.islink(dst):
                os.unlink(dst)
            else:
                shutil.rmtree(dst)

        os.symlink(src, dst, target_is_directory=True)
        print(f"[conf.py] Linked examples-src -> {src}")

    except Exception:
        # Windows fallback to copy if symlink fails
        try:
            subprocess.run(
                ["cmd", "/c", "mklink", "/J", dst, src], check=True, shell=True
            )
            print(f"[conf.py] Junction examples-src -> {src}")
        except Exception:
            shutil.copytree(src, dst)
            print(f"[conf.py] Copied examples to '{dst}' (symlink/junction unavailable)")

_mount_examples()

# --- Project information ---
project = "Python-CST Constructor Script"
author = "Martin Mihaylov"
release = "27.01.2026"

# --- General configuration ---
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_rtd_theme",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

autodoc_typehints = "description"

# --- Mock heavy imports ---
autodoc_mock_imports = [
    "numpy",
    "pandas",
    "matplotlib",
    "matlab",
]

# --- HTML output options ---
html_theme = "sphinx_rtd_theme"
html_theme_options = {
    "collapse_navigation": False,
    "navigation_depth": 4,
    "titles_only": False,
    "includehidden": True,
    "sticky_navigation": True,
}
html_static_path = []
