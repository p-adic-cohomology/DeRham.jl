from pathlib import Path
from fractions import Fraction
from math import prod
import json, tomllib, hashlib

root = Path(__file__).resolve().parent
native = tomllib.loads((root / "native-results.toml").read_text())
assert native["VERSION"] == "1.13.2"
assert native["basis_lengths"] == [0, 0, 1, 0, 0]
assert native["exact_g_coordinate"] == "1//48"
assert native["exact_Tg"] == "1//2"
exponents = native["exponent_vectors"]
assert len(exponents) == 126
assert len({tuple(ev) for ev in exponents}) == 126
assert all(len(ev) == 6 and sum(ev) == 4 and min(ev) >= 0 for ev in exponents)

def expected(ev):
    if any(e % 2 for e in ev):
        return Fraction(0)
    return Fraction(prod(prod(range(1, e, 2)) for e in ev), 2)

def modular(value, modulus):
    return (value.numerator * pow(value.denominator, -1, modulus)) % modulus

exact = [expected(ev) for ev in exponents]
assert exact.count(Fraction(0)) == 105
assert exact.count(Fraction(1, 2)) == 15
assert exact.count(Fraction(3, 2)) == 6
for probe in native["probes"]:
    modulus = int(probe["modulus"])
    assert modulus == 3 ** probe["M"]
    observed = list(map(int, probe["actual_row"]))
    target = [modular(c, modulus) for c in exact]
    assert target == list(map(int, probe["expected_row"]))
    assert observed == [modular(c / 2, modulus) for c in exact]
    assert sum(a != b for a,b in zip(observed,target)) == 21
    index = probe["g_column_one_based"] - 1
    assert exponents[index] == [2, 2, 0, 0, 0, 0]
    assert observed[index] == modular(Fraction(1,4),modulus)
    assert target[index] == modular(Fraction(1,2),modulus)
    assert all((2*a-b) % modulus == 0 for a,b in zip(observed,target))
for name, expected_hash in json.loads((root / "artifact-hashes.json").read_text()).items():
    assert hashlib.sha256((root/name).read_bytes()).hexdigest() == expected_hash, name
print(json.dumps({"status":"pass","independent_closed_form_columns":126,
    "native_entry_comparisons":378,"zero_columns":105,"nonzero_columns":21,
    "precisions":[2,3,4],"missing_factor":2},indent=2))
