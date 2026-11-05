"""Cell-value comparison of every .xlsx in two repository copies: compare_xlsx.py <ref_repo> <run_repo>."""
import os, sys, openpyxl

REF, RUN = sys.argv[1], sys.argv[2]


def same(a, b):
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return abs(a - b) <= 1e-6 * max(1.0, abs(a))
    return a == b


for root, dirs, names in os.walk(RUN):
    for n in names:
        if not n.endswith(".xlsx") or n.startswith("~$"):
            continue
        pb = os.path.join(root, n); rel = os.path.relpath(pb, RUN); pa = os.path.join(REF, rel)
        if not os.path.exists(pa):
            print("new_output", rel, sep="\t"); continue
        if os.path.getmtime(pb) <= os.path.getmtime(pa) + 1:
            print("not_rewritten", rel, sep="\t"); continue
        wa, wb = openpyxl.load_workbook(pa, data_only=True), openpyxl.load_workbook(pb, data_only=True)
        nd, first = 0, ""
        if wa.sheetnames != wb.sheetnames:
            print("differs", rel, "sheet names", sep="\t"); continue
        for s in wa.sheetnames:
            ra = list(wa[s].iter_rows(values_only=True)); rb = list(wb[s].iter_rows(values_only=True))
            if len(ra) != len(rb):
                nd += 1; first = first or f"{s}: {len(ra)} vs {len(rb)} rows"
            for i, (x, y) in enumerate(zip(ra, rb)):
                for j, (a, b) in enumerate(zip(x, y)):
                    if not same(a, b):
                        nd += 1; first = first or f"{s}!R{i + 1}C{j + 1}: {a!r} vs {b!r}"
        print("cells_equal" if nd == 0 else "differs", rel, nd, first, sep="\t")
