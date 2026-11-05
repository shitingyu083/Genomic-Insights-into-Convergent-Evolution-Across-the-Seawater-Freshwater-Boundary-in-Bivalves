"""Compare text outputs of the scratch run (<run>) against the pristine repository copy (<ref>).
Numeric cells are equal if |a-b| <= 1e-6 * max(1, |a|, |b|). Also maps _test_out/ files to their shipped twins.
usage: python compare_outputs.py <ref_repo> <run_repo> > report.tsv"""
import os, sys, re, hashlib, json
REF, RUN = sys.argv[1], sys.argv[2]
TXT = (".tsv", ".csv", ".txt", ".json", ".nwk", ".md")
NUM = re.compile(r"^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$")
TWINS = {  # scratch output -> shipped file (relative to repo)
    "_test_out/A_convergence_v4": "R1_reanalysis/_R1_work/reanalysis/A_convergence_v4",
    "_test_out/pgls_span_nov1/ref_pgls_span_results.tsv": "R1_reanalysis/_R1_work/analysis/pgls_span/pgls_span_results_nov1.tsv",
    "_test_out/pgls_span_nov1/ref_pgls_span_results_rawMb.tsv": "R1_reanalysis/_R1_work/analysis/pgls_span/pgls_span_results_rawMb_nov1.tsv",
    "_test_out/pgls_span_nov1/pgls_span_results.tsv": "R1_reanalysis/_R1_work/analysis/pgls_span/pgls_span_results_nov1_iqtree.tsv",
    "_test_out/busco_tmp": None,
}


def cells(text):
    return [re.split(r"[\t,]", l) for l in text.splitlines()]


def same_num(a, b):
    try:
        x, y = float(a), float(b)
    except ValueError:
        return a == b
    if x != x and y != y:
        return True
    return abs(x - y) <= 1e-6 * max(1.0, abs(x), abs(y))


def compare(pa, pb):
    a = open(pa, encoding="utf-8", errors="replace").read()
    b = open(pb, encoding="utf-8", errors="replace").read()
    if a == b:
        return "identical", 0, ""
    ca, cb = cells(a), cells(b)
    if len(ca) != len(cb):
        return "rows_differ", abs(len(ca) - len(cb)), f"{len(ca)} vs {len(cb)} lines"
    nd, mx, ex = 0, 0.0, ""
    for i, (ra, rb) in enumerate(zip(ca, cb)):
        if len(ra) != len(rb):
            nd += 1; ex = ex or f"line {i+1} width"; continue
        for u, v in zip(ra, rb):
            if u == v:
                continue
            if NUM.match(u.strip()) and NUM.match(v.strip()):
                if not same_num(u, v):
                    nd += 1; d = abs(float(u) - float(v)); mx = max(mx, d)
                    ex = ex or f"line {i+1}: {u} vs {v}"
            else:
                nd += 1; ex = ex or f"line {i+1}: {u[:40]!r} vs {v[:40]!r}"
    return ("numeric_equal" if nd == 0 else "differs"), nd, (ex + (f"; max abs diff {mx:.3g}" if mx else ""))


rows = []
for dp, dn, fn in os.walk(RUN):
    if "_test_logs" in dp:
        continue
    for f in fn:
        if not f.endswith(TXT):
            continue
        p = os.path.join(dp, f)
        rel = os.path.relpath(p, RUN).replace("\\", "/")
        ref = os.path.join(REF, rel)
        for k, v in TWINS.items():
            if rel.startswith(k):
                if v is None:
                    ref = None
                else:
                    ref = os.path.join(REF, v + rel[len(k):]) if rel != k and not v.endswith(".tsv") else os.path.join(REF, v)
        if ref is None:
            continue
        if not os.path.exists(ref):
            if os.path.getmtime(p) > os.path.getmtime(os.path.join(RUN, "README.md")):
                rows.append(("new_output", rel, 0, ""))
            continue
        if open(p, "rb").read() == open(ref, "rb").read():
            continue
        st, n, ex = compare(ref, p)
        rows.append((st, rel, n, ex))
for r in sorted(rows):
    print("\t".join(map(str, r)))
