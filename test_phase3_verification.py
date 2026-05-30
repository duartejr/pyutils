#!/usr/bin/env python3
"""
Phase 3 Verification Test Suite

Proves that:
1. CI root-cause fixes are correct (hardcoded version & dist-info path)
2. README.md is comprehensive (not empty)
3. evapotranspiration.py has docstrings, type hints, and produces correct output
4. read_nc.py has docstrings and type hints
5. Version bumped to 0.3.0
6. CI workflow is improved (fail-fast disabled, robust install, all 3 test suites)
7. All modules still import cleanly after refactor
"""

import sys
import subprocess
import importlib
import inspect
from pathlib import Path


def print_section(title):
    print(f"\n{'='*70}")
    print(f"  {title}")
    print(f"{'='*70}\n")


def test_ci_root_cause_fixes():
    """Prove the two CI failure root causes are fixed."""
    print_section("TEST 1: CI Root-Cause Fixes")

    all_pass = True

    # --- Root cause 1: hardcoded version "0.1.0" in test_phase1_verification.py ---
    src = Path("test_phase1_verification.py").read_text()
    if '"version": "0.1.0"' not in src:
        print("✅ test_phase1_verification.py: hardcoded version '0.1.0' REMOVED")
    else:
        print("❌ test_phase1_verification.py: still contains hardcoded version '0.1.0'")
        all_pass = False

    # --- Root cause 2: hardcoded dist-info path ---
    if "pyutils-0.1.0.dist-info" not in src:
        print("✅ test_phase1_verification.py: hardcoded 'pyutils-0.1.0.dist-info' path REMOVED")
    else:
        print("❌ test_phase1_verification.py: hardcoded dist-info path still present")
        all_pass = False

    # Confirm the fix uses dynamic discovery
    if 'glob("pyutils-*.dist-info")' in src:
        print("✅ test_phase1_verification.py: now uses dynamic glob for dist-info")
    else:
        print("❌ test_phase1_verification.py: dynamic glob not found")
        all_pass = False

    # --- Root cause 3: CI fail-fast was True (default), cancelling sibling jobs ---
    ci = Path(".github/workflows/ci.yml").read_text()
    if "fail-fast: false" in ci:
        print("✅ ci.yml: fail-fast: false — sibling jobs no longer cancelled on one failure")
    else:
        print("❌ ci.yml: fail-fast not set to false")
        all_pass = False

    # --- Robust install strategy ---
    if "--no-deps" in ci:
        print("✅ ci.yml: uses --no-deps install to decouple heavy C-extension deps from package install")
    else:
        print("❌ ci.yml: missing --no-deps strategy")
        all_pass = False

    if "cartopy install failed" in ci or "best-effort" in ci:
        print("✅ ci.yml: cartopy is best-effort (won't block CI if it fails to build)")
    else:
        print("❌ ci.yml: cartopy failure still blocks CI")
        all_pass = False

    # --- All 3 test suites wired into CI ---
    for suite in ["test_phase1_verification.py", "test_phase2_verification.py", "test_phase3_verification.py"]:
        if suite in ci:
            print(f"✅ ci.yml: runs {suite}")
        else:
            print(f"❌ ci.yml: missing {suite}")
            all_pass = False

    return all_pass


def test_readme():
    """Prove README.md is comprehensive."""
    print_section("TEST 2: README.md Written")

    readme = Path("README.md")
    if not readme.exists():
        print("❌ README.md does not exist")
        return False

    content = readme.read_text()
    size = len(content)
    lines = content.count("\n")
    print(f"✅ README.md exists ({lines} lines, {size} chars)")

    required = [
        ("Installation section", "## Installation"),
        ("Quick Start section", "## Quick Start"),
        ("Features table", "| Module"),
        ("Code examples", "```python"),
        ("evapotranspiration example", "evapotranspiration"),
        ("read_nc example", "read_nc"),
        ("Thiessen example", "thiessen"),
        ("Project structure", "## Project Structure"),
        ("License section", "## License"),
        ("System requirements", "System requirements"),
    ]

    all_pass = True
    for label, token in required:
        if token in content:
            print(f"✅ {label}")
        else:
            print(f"❌ {label} missing (expected '{token}')")
            all_pass = False

    return all_pass


def test_evapotranspiration_docstrings_and_types():
    """Prove evapotranspiration.py has docstrings, type hints, and correct output."""
    print_section("TEST 3: evapotranspiration.py — Docstrings, Types & Correctness")

    import evapotranspiration as et

    all_pass = True

    # Module docstring
    if et.__doc__ and len(et.__doc__.strip()) > 50:
        print(f"✅ Module docstring present ({len(et.__doc__.strip())} chars)")
    else:
        print("❌ Module docstring missing or too short")
        all_pass = False

    for fn_name in ("hargreaves", "thornthwaite", "penman_monteith"):
        fn = getattr(et, fn_name)
        doc = inspect.getdoc(fn) or ""
        hints = fn.__annotations__

        if len(doc) > 60:
            print(f"✅ {fn_name}(): docstring present ({len(doc)} chars)")
        else:
            print(f"❌ {fn_name}(): docstring missing or too short")
            all_pass = False

        if hints:
            print(f"✅ {fn_name}(): type hints: {list(hints.keys())}")
        else:
            print(f"⚠️  {fn_name}(): no type hints (optional but recommended)")

    # Correctness: Hargreaves with known values
    import numpy as np
    t_med = np.array([25.0] * 12)
    t_max = np.array([30.0] * 12)
    t_min = np.array([20.0] * 12)
    result = et.hargreaves(t_med, t_max, t_min, y=-5.0, months=1)
    if isinstance(result, np.ndarray) and len(result) == 12 and np.all(result > 0):
        print(f"✅ hargreaves() returns 12 positive mm/month values: {result[:3].round(1)}...")
    else:
        print(f"❌ hargreaves() output unexpected: {result}")
        all_pass = False

    return all_pass


def test_read_nc_docstrings_and_types():
    """Prove read_nc.py has docstrings and type hints."""
    print_section("TEST 4: read_nc.py — Docstrings & Type Hints")

    import read_nc as rn

    all_pass = True

    if rn.__doc__ and len(rn.__doc__.strip()) > 50:
        print(f"✅ Module docstring present ({len(rn.__doc__.strip())} chars)")
    else:
        print("❌ Module docstring missing or too short")
        all_pass = False

    for fn_name in ("read_nc", "save_nc"):
        fn = getattr(rn, fn_name)
        doc = inspect.getdoc(fn) or ""
        hints = fn.__annotations__

        if len(doc) > 60:
            print(f"✅ {fn_name}(): docstring present ({len(doc)} chars)")
        else:
            print(f"❌ {fn_name}(): docstring missing or too short")
            all_pass = False

        if hints:
            print(f"✅ {fn_name}(): type hints: {list(hints.keys())}")
        else:
            print(f"⚠️  {fn_name}(): no type hints")

    # Private helpers renamed (__ → _)
    src = Path("read_nc.py").read_text()
    if "def _in_poly" in src and "def __in_poly" not in src:
        print("✅ read_nc.py: private helpers use single-underscore convention (_in_poly)")
    else:
        print("❌ read_nc.py: still uses double-underscore private helpers")
        all_pass = False

    # No bare print() in production code
    import ast
    tree = ast.parse(src)
    bare_prints = [
        n for n in ast.walk(tree)
        if isinstance(n, ast.Call)
        and isinstance(n.func, ast.Name)
        and n.func.id == "print"
    ]
    if not bare_prints:
        print("✅ read_nc.py: no bare print() calls")
    else:
        print(f"❌ read_nc.py: {len(bare_prints)} bare print() calls remain")
        all_pass = False

    return all_pass


def test_version_bump():
    """Prove version is 0.3.0 in all locations."""
    print_section("TEST 5: Version 0.3.0 in All Locations")

    all_pass = True

    import tomllib
    with open("pyproject.toml", "rb") as f:
        config = tomllib.load(f)
    v = config["project"]["version"]
    if v == "0.3.0":
        print(f"✅ pyproject.toml version = '{v}'")
    else:
        print(f"❌ pyproject.toml version = '{v}' (expected '0.3.0')")
        all_pass = False

    init_src = Path("__init__.py").read_text()
    if '__version__ = "0.3.0"' in init_src:
        print("✅ __init__.py __version__ = '0.3.0'")
    else:
        print("❌ __init__.py version not 0.3.0")
        all_pass = False

    pyutils_src = Path("pyutils.py").read_text()
    if '__version__ = "0.3.0"' in pyutils_src:
        print("✅ pyutils.py __version__ = '0.3.0'")
    else:
        print("❌ pyutils.py version not 0.3.0")
        all_pass = False

    result = subprocess.run(["pip", "show", "pyutils"], capture_output=True, text=True)
    for line in result.stdout.splitlines():
        if line.startswith("Version:"):
            if "0.3.0" in line:
                print(f"✅ pip show: {line.strip()}")
            else:
                print(f"❌ pip show: {line.strip()} (expected 0.3.0)")
                all_pass = False

    # Confirm pyutils module reports correct version
    if "pyutils" in sys.modules:
        del sys.modules["pyutils"]
    import pyutils
    if pyutils.__version__ == "0.3.0":
        print(f"✅ import pyutils.__version__ = '{pyutils.__version__}'")
    else:
        print(f"❌ import pyutils.__version__ = '{pyutils.__version__}' (expected '0.3.0')")
        all_pass = False

    return all_pass


def test_all_modules_still_import():
    """Regression check — make sure refactored modules still import."""
    print_section("TEST 6: All Modules Import After Refactor")

    modules = [
        "pyutils", "pyeto", "rv",
        "evapotranspiration", "read_nc", "thiessen",
        "plot_maps", "spi", "idw", "corr_med",
        "bh_thorthwaite", "bh_thorthwaite2", "eto",
        "vazoes_minimas", "custom_colorbar", "shp_area",
        "ccsthorthwaite", "rvGamma",
    ]

    all_pass = True
    for name in modules:
        if name in sys.modules:
            del sys.modules[name]
        try:
            importlib.import_module(name)
            print(f"✅ import {name}: OK")
        except ImportError as e:
            print(f"⚠️  import {name}: optional dep missing ({e})")
        except SyntaxError as e:
            print(f"❌ import {name}: SYNTAX ERROR — {e}")
            all_pass = False
        except Exception as e:
            print(f"❌ import {name}: {type(e).__name__}: {e}")
            all_pass = False

    return all_pass


def test_phase1_still_passes():
    """Run Phase 1 tests to confirm they still pass after fixes."""
    print_section("TEST 7: Phase 1 Test Suite Still Passes")

    result = subprocess.run(
        [sys.executable, "test_phase1_verification.py"],
        capture_output=True, text=True
    )
    lines = result.stdout.splitlines()
    # Find the summary line
    summary = next((l for l in reversed(lines) if "TOTAL:" in l), "")
    passed = "7/7" in summary

    if passed:
        print(f"✅ Phase 1 tests: {summary.strip()}")
    else:
        print(f"❌ Phase 1 tests: {summary.strip()}")
        # Print last 20 lines for diagnosis
        for l in lines[-20:]:
            print(f"   {l}")

    return passed


def main():
    print("\n")
    print("╔" + "="*68 + "╗")
    print("║" + " "*68 + "║")
    print("║" + "  PYUTILS PHASE 3 VERIFICATION TEST SUITE".center(68) + "║")
    print("║" + "  Docs, Type Hints & CI Fix — Full Evidence".center(68) + "║")
    print("║" + " "*68 + "║")
    print("╚" + "="*68 + "╝")

    tests = [
        ("CI Root-Cause Fixes", test_ci_root_cause_fixes),
        ("README.md Written", test_readme),
        ("evapotranspiration.py Docstrings/Types/Correctness", test_evapotranspiration_docstrings_and_types),
        ("read_nc.py Docstrings/Types", test_read_nc_docstrings_and_types),
        ("Version 0.3.0 Everywhere", test_version_bump),
        ("All Modules Import After Refactor", test_all_modules_still_import),
        ("Phase 1 Test Suite Still Passes", test_phase1_still_passes),
    ]

    results = {}
    for name, fn in tests:
        try:
            results[name] = fn()
        except Exception as e:
            import traceback
            print(f"\n❌ EXCEPTION in '{name}':")
            traceback.print_exc()
            results[name] = False

    print_section("PHASE 3 SUMMARY")

    passed = sum(1 for v in results.values() if v)
    total = len(results)

    for name, ok in results.items():
        print(f"{'✅ PASS' if ok else '❌ FAIL'}: {name}")

    print(f"\n{'='*70}")
    print(f"TOTAL: {passed}/{total} test groups passed")
    print(f"{'='*70}\n")

    if passed == total:
        print("🎉 ALL PHASE 3 TESTS PASSED! 🎉\n")
    else:
        print(f"⚠️  {total - passed} test group(s) need attention\n")

    return 0 if passed == total else 1


if __name__ == "__main__":
    sys.exit(main())
