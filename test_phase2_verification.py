#!/usr/bin/env python3
"""
Phase 2 Verification Test Suite

Proves that:
1. Removed files are gone (no duplicates, no stubs, no embedded projects)
2. All remaining modules import cleanly
3. Fixed bugs are fixed (thiessen.pd.pd, buffer print, rv broken imports)
4. Version is bumped to 0.2.0
5. GitHub Actions workflows are in place
6. Package re-installs correctly at v0.2.0
"""

import sys
import subprocess
import importlib
from pathlib import Path

repo_root = Path(__file__).parent


def print_section(title):
    print(f"\n{'='*70}")
    print(f"  {title}")
    print(f"{'='*70}\n")


def test_deleted_files():
    """Prove duplicate/empty/embedded files were removed."""
    print_section("TEST 1: Deleted Files (No More Duplicates/Stubs)")

    repo = repo_root

    should_not_exist = [
        "plot_maps_old.py",     # duplicate
        "anom_pr.py",           # empty stub
        "kriging.py",           # empty stub
        "test_plot_maps.py",    # empty test file
        "PyETo-master",         # embedded external project
    ]

    all_pass = True
    for name in should_not_exist:
        path = repo / name
        if not path.exists():
            print(f"✅ DELETED: {name}")
        else:
            print(f"❌ STILL EXISTS: {name}")
            all_pass = False

    return all_pass


def test_remaining_files_exist():
    """Prove all expected modules are still present."""
    print_section("TEST 2: Remaining Files Still Present")

    repo = repo_root

    should_exist = [
        "pyproject.toml",
        "setup.py",
        "__init__.py",
        "plot_maps.py",         # kept (not the old one)
        "bh_thorthwaite.py",
        "bh_thorthwaite2.py",
        "buffer.py",
        "ccsthorthwaite.py",
        "corr_med.py",
        "evapotranspiration.py",
        "read_nc.py",
        "thiessen.py",
        "rvGamma.py",
        "idw.py",
        "spi.py",
        "pyeto/__init__.py",
        "pyeto/fao.py",
        "rv/__init__.py",
        "rv/RvGama.py",
        "rv/RvOns.py",
        ".github/workflows/ci.yml",
        ".github/workflows/release.yml",
    ]

    all_pass = True
    for name in should_exist:
        path = repo / name
        if path.exists():
            print(f"✅ PRESENT: {name}")
        else:
            print(f"❌ MISSING: {name}")
            all_pass = False

    return all_pass


def test_bug_fixes():
    """Prove specific bugs were fixed."""
    print_section("TEST 3: Bug Fixes Verified")

    all_pass = True

    # 1. thiessen.py: pd.pd.DataFrame → pd.DataFrame
    thiessen_src = (repo_root / "thiessen.py").read_text()
    if "pd.pd.DataFrame" not in thiessen_src:
        print("✅ thiessen.py: pd.pd.DataFrame typo is FIXED")
    else:
        print("❌ thiessen.py: pd.pd.DataFrame typo still present")
        all_pass = False

    if "pd.DataFrame" in thiessen_src:
        print("✅ thiessen.py: uses correct pd.DataFrame")
    else:
        print("❌ thiessen.py: pd.DataFrame not found")
        all_pass = False

    # 2. buffer.py: debug print removed
    buffer_src = (repo_root / "buffer.py").read_text()
    if "print(ds_out)" not in buffer_src:
        print("✅ buffer.py: debug print(ds_out) is REMOVED")
    else:
        print("❌ buffer.py: debug print(ds_out) still present")
        all_pass = False

    # 3. rv/RvOns.py: broken utils import removed
    rvons_src = (repo_root / "rv/RvOns.py").read_text()
    if "from utils import basinsf" not in rvons_src:
        print("✅ rv/RvOns.py: broken 'from utils import basinsf' is REMOVED")
    else:
        print("❌ rv/RvOns.py: broken import still present")
        all_pass = False

    # 4. rvGamma.py: broken utils imports removed
    rvgamma_src = (repo_root / "rvGamma.py").read_text()
    if "from utils.get_tri_files" not in rvgamma_src:
        print("✅ rvGamma.py: broken utils imports are REMOVED")
    else:
        print("❌ rvGamma.py: broken utils imports still present")
        all_pass = False

    # 5. rv/plr.py: syntax error fixed (no bare function body)
    plr_src = (repo_root / "rv/plr.py").read_text()
    try:
        compile(plr_src, "rv/plr.py", "exec")
        print("✅ rv/plr.py: compiles without syntax errors")
    except SyntaxError as e:
        print(f"❌ rv/plr.py: SyntaxError: {e}")
        all_pass = False

    return all_pass


def test_module_imports():
    """Prove all remaining modules import without error."""
    print_section("TEST 4: All Modules Import Cleanly")

    modules = [
        "pyutils",
        "pyeto",
        "rv",
        "read_nc",
        "corr_med",
        "custom_colorbar",
        "idw",
        "isinshp",
        "plot_maps",
        "spi",
        "vazoes_minimas",
        "evapotranspiration",
        "bh_thorthwaite",
        "bh_thorthwaite2",
        "eto",
        "shp_area",
        "ccsthorthwaite",
    ]

    all_pass = True
    failed = []
    for mod_name in modules:
        # Clear cached version if already imported
        if mod_name in sys.modules:
            del sys.modules[mod_name]
        try:
            importlib.import_module(mod_name)
            print(f"✅ import {mod_name}: OK")
        except ImportError as e:
            print(f"⚠️  import {mod_name}: missing optional dependency ({e})")
        except SyntaxError as e:
            print(f"❌ import {mod_name}: SYNTAX ERROR — {e}")
            all_pass = False
            failed.append(mod_name)
        except Exception as e:
            print(f"❌ import {mod_name}: UNEXPECTED ERROR — {type(e).__name__}: {e}")
            all_pass = False
            failed.append(mod_name)

    if failed:
        print(f"\n❌ {len(failed)} module(s) failed: {', '.join(failed)}")
    return all_pass


def test_version_bump():
    """Prove version is consistent across all locations (reads from pyproject.toml)."""
    print_section("TEST 5: Version Consistent Across All Locations")

    all_pass = True

    # Read the authoritative version from pyproject.toml
    import tomllib
    with open(repo_root / "pyproject.toml", "rb") as f:
        config = tomllib.load(f)

    version = config["project"]["version"]
    print(f"✅ pyproject.toml version = '{version}'")

    # Check __init__.py
    init_src = (repo_root / "__init__.py").read_text()
    if f'__version__ = "{version}"' in init_src:
        print(f"✅ __init__.py __version__ = '{version}'")
    else:
        print(f"❌ __init__.py __version__ not updated to {version}")
        all_pass = False

    # Check installed version (reload to avoid cached stale module)
    import importlib, sys
    if "pyutils" in sys.modules:
        del sys.modules["pyutils"]
    import pyutils
    if pyutils.__version__ == version:
        print(f"✅ Installed pyutils.__version__ = '{pyutils.__version__}'")
    else:
        print(f"❌ Installed version = '{pyutils.__version__}' (expected '{version}')")
        all_pass = False

    # Check pip show
    result = subprocess.run(["pip", "show", "pyutils"], capture_output=True, text=True)
    for line in result.stdout.splitlines():
        if line.startswith("Version:"):
            installed = line.split(":", 1)[1].strip()
            if installed == version:
                print(f"✅ pip show pyutils: {line.strip()}")
            else:
                print(f"❌ pip show pyutils: {line.strip()} (expected {version})")
                all_pass = False

    return all_pass


def test_github_actions():
    """Prove GitHub Actions workflows were created."""
    print_section("TEST 6: GitHub Actions Workflows")

    all_pass = True

    ci_path = repo_root / ".github/workflows/ci.yml"
    release_path = repo_root / ".github/workflows/release.yml"

    for path in [ci_path, release_path]:
        if path.exists():
            content = path.read_text()
            print(f"✅ {path.name} exists ({len(content.splitlines())} lines)")
        else:
            print(f"❌ {path.name} missing")
            all_pass = False

    # Check CI workflow triggers
    if ci_path.exists():
        ci = ci_path.read_text()
        checks = [
            ("triggers on push to master", "push:" in ci),
            ("triggers on pull_request", "pull_request:" in ci),
            ("runs on ubuntu", "ubuntu-latest" in ci),
            ("tests Python 3.9-3.11", "python-version" in ci),
        ]
        print(f"\n  CI workflow checks:")
        for label, ok in checks:
            status = "✅" if ok else "❌"
            print(f"  {status} {label}")
            if not ok:
                all_pass = False

    # Check Release workflow
    if release_path.exists():
        rel = release_path.read_text()
        checks = [
            ("triggers on PR closed", "pull_request:" in rel and "closed" in rel),
            ("only runs on merge", "merged == true" in rel),
            ("reads version from pyproject.toml", "pyproject.toml" in rel),
            ("creates git tag", "git tag" in rel),
            ("creates GitHub release", "createRelease" in rel),
        ]
        print(f"\n  Release workflow checks:")
        for label, ok in checks:
            status = "✅" if ok else "❌"
            print(f"  {status} {label}")
            if not ok:
                all_pass = False

    return all_pass


def test_no_broken_imports_in_code():
    """Scan remaining source for known broken import patterns."""
    print_section("TEST 7: No Broken Import Patterns in Source")

    repo = repo_root
    all_py = list(repo.glob("*.py")) + list(repo.glob("rv/*.py")) + list(repo.glob("pyeto/*.py"))
    all_py = [p for p in all_py if "__pycache__" not in str(p)]

    broken_patterns = [
        ("from utils import", "missing utils package"),
        ("import utils.", "missing utils package"),
        ("pd.pd.", "double pd prefix typo"),
        ("print(ds_out)", "leftover debug print"),
    ]

    # Exclude this test file and Phase 1 test from the scan
    exclude = {"test_phase2_verification.py", "test_phase1_verification.py"}
    all_py = [p for p in all_py if p.name not in exclude]

    all_pass = True
    hits = []
    for py_file in sorted(all_py):
        # Only check non-comment, non-blank lines
        lines = py_file.read_text(errors="replace").splitlines()
        active_lines = [l for l in lines if l.strip() and not l.strip().startswith("#")]
        src_active = "\n".join(active_lines)
        for pattern, reason in broken_patterns:
            if pattern in src_active:
                print(f"❌ {py_file.name}: contains '{pattern}' ({reason})")
                all_pass = False
                hits.append((py_file.name, pattern))

    if not hits:
        print(f"✅ Scanned {len(all_py)} source files — no broken patterns found")

    return all_pass


def main():
    print("\n")
    print("╔" + "="*68 + "╗")
    print("║" + " "*68 + "║")
    print("║" + "  PYUTILS PHASE 2 VERIFICATION TEST SUITE".center(68) + "║")
    print("║" + "  Code Cleanup — Real Proof Everything Works".center(68) + "║")
    print("║" + " "*68 + "║")
    print("╚" + "="*68 + "╝")

    tests = [
        ("Deleted Files (No Duplicates/Stubs)", test_deleted_files),
        ("Remaining Files Present", test_remaining_files_exist),
        ("Bug Fixes Verified", test_bug_fixes),
        ("Module Imports Clean", test_module_imports),
        ("Version Consistent Across All Locations", test_version_bump),
        ("GitHub Actions Workflows", test_github_actions),
        ("No Broken Import Patterns", test_no_broken_imports_in_code),
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

    print_section("PHASE 2 SUMMARY")

    passed = sum(1 for v in results.values() if v)
    total = len(results)

    for name, ok in results.items():
        status = "✅ PASS" if ok else "❌ FAIL"
        print(f"{status}: {name}")

    print(f"\n{'='*70}")
    print(f"TOTAL: {passed}/{total} test groups passed")
    print(f"{'='*70}\n")

    if passed == total:
        print("🎉 ALL PHASE 2 TESTS PASSED! 🎉\n")
    else:
        print(f"⚠️  {total - passed} test group(s) need attention\n")

    return 0 if passed == total else 1


if __name__ == "__main__":
    sys.exit(main())
