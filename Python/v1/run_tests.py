"""Run all unit tests without pytest:  python3 run_tests.py"""
import sys, glob, importlib.util
sys.path.insert(0, "src")
passed = failed = 0
for f in sorted(glob.glob("tests/unit/test_*.py")):
    spec = importlib.util.spec_from_file_location("m", f); m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    for n in sorted(x for x in dir(m) if x.startswith("test_")):
        try: getattr(m, n)(); passed += 1; print("PASS", f, n)
        except Exception as e: failed += 1; print("FAIL", f, n, repr(e))
print(f"{passed} passed, {failed} failed"); sys.exit(1 if failed else 0)
