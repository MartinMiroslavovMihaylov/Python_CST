# docs-src/conf.py
# Configuration file for Sphinx documentation builder.

import os
import sys
import shutil
import subprocess

# --- Add main repo to sys.path ---
CODE_DIR = os.environ.get("CODE_DIR")
if not CODE_DIR or not os.path.isdir(CODE_DIR):
    CODE_DIR = os.path.abspath("../repo")  # fallback if CODE_DIR not set

sys.path.insert(0, CODE_DIR)
print(f"[conf.py] Added CODE_DIR to sys.path: {CODE_DIR}")

def _mount_examples():
    src = os.path.join(CODE_DIR, "Examples")
    dst = os.path.join(os.path.dirname(__file__), "examples-src")

    if not os.path.isdir(src):
        print(f"[conf.py] No code examples at {src}; skipping.", file=sys.stderr)
        return

    # Remove existing examples-src
    if os.path.exists(dst):
        try:
            if os.path.islink(dst):
                os.unlink(dst)
            elif os.path.isdir(dst):
                shutil.rmtree(dst)
            else:
                os.remove(dst)
        except Exception as e:
            print(f"[conf.py] Failed to remove old examples-src: {e}", file=sys.stderr)
            return

    # Try symlink first (works on Linux / GitHub Actions)
    try:
        os.symlink(src, dst, target_is_directory=True)
        print(f"[conf.py] Linked examples-src -> {src}")
        return
    except Exception:
        pass

    # Try Windows junction
    if os.name == "nt":
        try:
            subprocess.run(["cmd", "/c", "mklink", "/J", dst, src], check=True, shell=True)
            print(f"[conf.py] Junction examples-src -> {src}")
            return
        except Exception:
            pass

    # Fallback: copy examples
    try:
        shutil.copytree(src, dst)
        print(f"[conf.py] Copied examples to '{dst}' (symlink/junction unavailable)")
    except Exception as e:
        print(f"[conf.py] Failed to copy examples: {e}", file=sys.stderr)

_mount_examples()

# --- Project info ---
project = "Python-CST Constructor Script"
author = "Martin Mihaylov"
release = "27.01.2026"

# --- Sphinx extensions ---
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

# --- Mock heavy imports only on CI ---
autodoc_mock_imports = [
    "numpy",
    "pandas",
    "matplotlib",
    "matlab",
    "scipy"
]

# Detect CI environment
if os.environ.get("GITHUB_ACTIONS") == "true":
    autodoc_mock_imports.append("cst")
    print("[conf.py] Running on GitHub Actions: 'cst' will be mocked.")

# --- HTML options ---
html_theme = "sphinx_rtd_theme"
html_theme_options = {
    "collapse_navigation": False,
    "navigation_depth": 4,
    "titles_only": False,
    "includehidden": True,
    "sticky_navigation": True,
}
html_static_path = []
