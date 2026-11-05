"""Fig. 2f replacement panel drawn from the archived ancestral reconstructions (R0 layout).

Source: Results_Dec31/MultipleSequenceAlignments/Final_Unique_Alignments_dedup_species and
Ancestral_Results_Final_Species (RAxML-NG 1.2.2 --ancestral, LG+G+F, species_tree_fixed.nwk);
the archived ATG18 (OG0004469) run had failed and was rerun on the pruned tree (fig2f_asr.py).
Columns: stem-branch changes shared by >= 2 freshwater units, all ancestral/derived PP >= 0.8.
"""
import pathlib
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

H = pathlib.Path(__file__).resolve().parent
PKG = H.parents[2] / 'SUBMIT_R1_revised_20261024'
OUT = PKG / 'Figures_R1' / 'panels_for_Fig2'
mpl.rcParams.update({'font.family': 'Arial', 'font.weight': 'bold', 'pdf.fonttype': 42, 'ps.fonttype': 42})
GENES = [('ATG7', 'OG0005680', 709), ('ATG18', 'OG0004469', 516), ('CNX', 'OG0004287', 623)]
UNITS = {'Unionida': ['Sol', 'Pst', 'Awo', 'Ude'], 'Dreissenidae': ['Dpo', 'Dro'], 'Corbicula': ['Cfl'],
         'Limnoperna': ['Lfo']}
UNIT_OF = {sp: u for u, l in UNITS.items() for sp in l}
R0_ORDER = ('Sol Pst Awo Ude Cbi Sco Tna Mar Dro Dpo Mve Sso Cfl Ama Csi Rph Mme Sbr Pmag Mye Cfa Aja Yba Air '
            'Apu Pmax Apl Lfo Bpl Mvi Mco Mga Mca Pfu Pim Sgl Oed Cvi Cni Cho Can Cgi').split()
RED, BLUE, DARK, GREY = '#ca2a28', '#2371a9', '#333333', '#888888'
AA = {**dict.fromkeys('AILMFWV', '#4a7fc1'), **dict.fromkeys('KR', '#d0473f'), **dict.fromkeys('DE', '#b05cc6'),
      **dict.fromkeys('NQST', '#3fa34d'), 'G': '#e08a2c', 'P': '#c9b52b', **dict.fromkeys('HY', '#2aa6a6'),
      'C': '#e57fa0'}


def fasta(p):
    seq, k = {}, None
    for line in open(p):
        line = line.strip()
        if line.startswith('>'):
            k = line[1:].replace('.pep', '')
            seq[k] = []
        elif k:
            seq[k].append(line)
    return {k: ''.join(v) for k, v in seq.items()}


def parse(s):
    s, i = s.strip().rstrip(';'), [0]

    def node():
        ch = []
        if s[i[0]] == '(':
            i[0] += 1
            while True:
                ch.append(node())
                c = s[i[0]]
                i[0] += 1
                if c == ')':
                    break
        j = i[0]
        while i[0] < len(s) and s[i[0]] not in ',();':
            i[0] += 1
        return [s[j:i[0]].split(':')[0].replace('.pep', ''), ch]
    return node()


def tips(n):
    return [n[0]] if not n[1] else [t for c in n[1] for t in tips(c)]


def order(n, rank):
    if n[1]:
        n[1].sort(key=lambda c: min(rank.get(t, 1e9) for t in tips(c)))
        for c in n[1]:
            order(c, rank)


def height(n):
    return 0 if not n[1] else 1 + max(height(c) for c in n[1])


def main():
    tree = parse(open(H / 'species_tree_fixed.nwk').read())
    tree[1] = [c for c in tree[1] if c[0] != 'Lgi']  # outgroup not drawn, as in R0
    order(tree, {sp: k for k, sp in enumerate(R0_ORDER)})
    rows = tips(tree)
    assert sorted(rows) == sorted(R0_ORDER), set(rows) ^ set(R0_ORDER)
    if rows != R0_ORDER:
        print('tree order differs from R0 display order:', [(a, b) for a, b in zip(rows, R0_ORDER) if a != b][:6])
    y = {sp: k for k, sp in enumerate(rows)}
    s1 = pd.read_csv(PKG / 'Supplementary_Table_S1_retained_bivalve_genomes_R1.tsv', sep='\t')
    names = dict(s1[['ID', 'Species']].values)
    d = pd.read_csv(H / 'fig2f_shared_changes.tsv', sep='\t')
    d = d[d['all_pp_ge_0.8'].astype(bool)].copy()
    d['units'] = d.changes.map(lambda s: [c.split(':')[0] for c in s.split(';')])

    TW, NX0, NX1, GAP = 7.0, 7.0, 18.6, 0.8
    cols, blocks, x = [], [], NX1 + 0.4
    for g, og, ncol in GENES:
        aln = fasta(H / f'{og}_species.fasta')
        assert {len(v) for v in aln.values()} == {ncol}, og
        x0 = x
        for r in d[d.gene == g].sort_values('aln_col').itertuples():
            cols.append((g, og, x, r, aln))
            x += 1
        blocks.append((g, x0, x))
        x += GAP
    nrow, W = len(rows), x
    fig = plt.figure(figsize=(W * 0.152, (nrow + 8.6) * 0.152))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(-0.4, W)
    ax.set_ylim(nrow + 5.4, -3.2)
    ax.axis('off')

    n = tree
    while 'Mye' in tips(n):  # largest clade with Sol but without Mye = Heteroconchia
        n = next(c for c in n[1] if 'Sol' in tips(c))
    het = set(tips(n))
    print('rows:', ' '.join(rows), '| Heteroconchia', len(het))
    for lab, members, fc in (('Heteroconchia', het, '#dceaf6'), ('Pteriomorphia', set(rows) - het, '#f7e0e0')):
        ks = [y[s] for s in members]
        a, b = min(ks), max(ks) + 1
        assert b - a == len(ks), lab  # each clade is one contiguous block of rows
        ax.add_patch(Rectangle((NX0, a), NX1 - NX0, b - a, fc=fc, ec='none', zorder=0))
        ax.text(NX0 + 0.55, (a + b) / 2, lab, rotation=90, ha='center', va='center', fontsize=10, zorder=2)
    step, fw = TW / height(tree), set(UNIT_OF)

    def draw(n):
        xn = TW - height(n) * step
        if not n[1]:
            return y[n[0]] + 0.5
        ys = [draw(c) for c in n[1]]
        for yc, c in zip(ys, n[1]):
            ax.plot([xn, TW - height(c) * step], [yc, yc], color=RED if set(tips(c)) <= fw else 'black',
                    lw=0.9, solid_capstyle='butt', zorder=3)
        ax.plot([xn, xn], [min(ys), max(ys)], color=RED if set(tips(n)) <= fw else 'black', lw=0.9, zorder=3)
        return (min(ys) + max(ys)) / 2
    yr = draw(tree)
    ax.plot([-0.3, 0], [yr, yr], color='black', lw=0.9)
    for sp in rows:
        ax.text(NX0 + 1.15, y[sp] + 0.5, names[sp], style='italic', fontsize=7, va='center', ha='left', zorder=2)
    yd, yg = nrow + 0.9, nrow + 2.25
    for g, og, xc, r, aln in cols:
        c, dv = r.aln_col - 1, r.derived
        for u in r.units:  # every changed unit must carry the derived residue in at least one tip
            assert any(sp in aln and aln[sp][c] == dv for sp in UNITS[u]), (g, r.aln_col, u)
        hit = {sp for u in r.units for sp in UNITS[u]}
        for sp in rows:
            if sp not in aln:
                ax.add_patch(Rectangle((xc, y[sp]), 1, 1, fc='#d9d9d9', ec='white', lw=0.4))
                continue
            ch = aln[sp][c]
            der = ch == dv
            ax.add_patch(Rectangle((xc, y[sp]), 1, 1, fc=AA.get(ch, '#f2f2f2'), alpha=1 if der else 0.22,
                                   ec='white', lw=0.4))
            ax.text(xc + 0.5, y[sp] + 0.54, ch, ha='center', va='center', fontsize=6.3,
                    color='black' if der else GREY)
            if der and sp in hit:
                ax.add_patch(Rectangle((xc + 0.07, y[sp] + 0.07), 0.86, 0.86, fc='none', ec=RED, lw=1.1, zorder=4))
        ax.text(xc + 0.5, -0.25, str(r.aln_col), rotation=90, ha='center', va='bottom', fontsize=6.3)
        ax.scatter(xc + 0.5, yd, s=15, c=BLUE if r.kind == 'parallel' else DARK, edgecolors='none', zorder=4)
    for g, x0, x1 in blocks:
        ax.plot([x0 + 0.1, x1 - 0.1], [yd + 0.75, yd + 0.75], color='black', lw=0.8)
        ax.text((x0 + x1) / 2, yg, g, ha='center', va='center', fontsize=10)
    ly = nrow + 3.75
    for k, (col, lab) in enumerate(((BLUE, 'Parallel change'), (DARK, 'Convergent change'))):
        ax.scatter(0.3, ly + 1.25 * k, s=15, c=col, edgecolors='none')
        ax.text(0.9, ly + 1.25 * k, lab, va='center', fontsize=7.5)
    ax.add_patch(Rectangle((10.5, ly - 0.43), 0.86, 0.86, fc='none', ec=RED, lw=1.1))
    ax.text(11.9, ly, 'Derived residue on a freshwater stem branch', va='center', fontsize=7.5)
    ax.plot([10.4, 11.4], [ly + 1.25, ly + 1.25], color=RED, lw=1.2)
    ax.text(11.9, ly + 1.25, 'Freshwater lineage', va='center', fontsize=7.5)
    ax.add_patch(Rectangle((33.0, ly - 0.5), 1, 1, fc='#d9d9d9', ec='none'))
    ax.text(34.4, ly, 'Taxon absent from orthogroup', va='center', fontsize=7.5)
    fig.savefig(OUT / 'Fig2f_R1.png', dpi=300)
    fig.savefig(OUT / 'Fig2f_R1.pdf')
    hdr = ['gene', 'orthogroup', 'aln_col', 'Mye_pos', 'derived', 'kind', 'stem_changes_anc(PP)->der(PP)',
           'n_comparators_with_derived', 'n_comparators'] + rows
    src = [[g, og, r.aln_col, '' if pd.isna(r.Mye_pos) else int(r.Mye_pos), r.derived, r.kind, r.changes,
            r.n_comparators_with_derived, r.n_comparators] + [aln[sp][r.aln_col - 1] if sp in aln else 'absent'
                                                             for sp in rows] for g, og, xc, r, aln in cols]
    pd.DataFrame(src, columns=hdr).to_csv(OUT / 'Fig2f_source_data.tsv', sep='\t', index=False)
    print(len(cols), 'columns', {g: round(x1 - x0) for g, x0, x1 in blocks}, 'inches', fig.get_size_inches().round(2))


if __name__ == '__main__':
    main()
