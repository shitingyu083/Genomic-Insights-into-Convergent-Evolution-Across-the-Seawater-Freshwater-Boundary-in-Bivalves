#!/usr/bin/env python3
"""R1 Table S8: Fig. 2f candidate sites (ATG7, ATG18/WIPI2, CNX), run on the server.
Reads the archived OrthoFinder species alignments and RAxML-NG 1.2.2 marginal ancestral states
(--ancestral, LG+G+F, fixed species tree; Results_Dec31/Ancestral_Results_Final_Species).
The archived OG0004469 (WIPI2) run failed (Apu, Cgi, Pim absent from the alignment), so it is rerun
here with the same command on the species tree pruned to the 40 taxa present.
Outputs: ~/R1_revision_20261024/fig2f_asr/{fig2f_shared_changes.tsv, fig2f_figcols.tsv, fig2f_summary.json}"""
import os, re, json, subprocess, collections
H = os.path.expanduser('~')
O = os.environ.get("R0_DATA", H + '/###qutongjinhua_0710cxcf') + '/ortho_clean/orthofinder_output/Results_Dec31'
AL = O + '/MultipleSequenceAlignments/Final_Unique_Alignments_dedup_species'
AR = O + '/Ancestral_Results_Final_Species'
TREE = O + '/PhylipAlignments_fixed_strict/species_tree_fixed.nwk'
W = os.environ.get("R1_WORK", H + '/R1_revision_20261024') + '/fig2f_asr'
UNITS = {'Unionida': ['Sol', 'Pst', 'Awo', 'Ude'], 'Dreissenidae': ['Dpo', 'Dro'], 'Corbicula': ['Cfl'], 'Limnoperna': ['Lfo']}
FW = {t for v in UNITS.values() for t in v}
# Fig. 2f (R0): column labels, comparator (background) residues, derived residues shown for Ude/Awo/Dpo/Dro/Cfl
FIG = {'ATG7': ('OG0005680', [19, 27, 22, 28, 65, 70, 71, 82, 109, 140, 145, 152, 181, 194, 217, 221, 226, 227],
                'LDRPIHKRVMNSSETNLQ', 'MEKINYRKILYTKDSPAP'),
       'ATG18': ('OG0004469', [8, 9, 18, 25, 127, 139, 168, 198, 232], 'DEVANPDYL', 'EQLGKHEHI'),
       'CNX': ('OG0004287', [47, 54, 195, 333, 401, 423], 'FTAFAV', 'LSSLSI')}
PPMIN = 0.8


def cl(l): return l[:-4] if l.endswith('.pep') else l


def fasta(p):
    s, k = {}, None
    for l in open(p):
        l = l.strip()
        if l.startswith('>'): k = cl(l[1:].split()[0]); s[k] = []
        elif k: s[k].append(l)
    return {k: ''.join(v).upper() for k, v in s.items()}


def parse(t):
    t = t.strip().rstrip(';'); i = [0]
    def node():
        ch = []
        if t[i[0]] == '(':
            i[0] += 1
            while True:
                ch.append(node()); c = t[i[0]]; i[0] += 1
                if c == ')': break
        m = re.match(r"([^,():;]*)(:[-0-9.eE+]+)?", t[i[0]:]); i[0] += m.end()
        return [m.group(1), float(m.group(2)[1:]) if m.group(2) else 0.0, ch]
    return node()


def prune(n, keep):
    lab, bl, ch = n
    if not ch: return n if cl(lab) in keep else None
    ch = [c for c in (prune(c, keep) for c in ch) if c]
    if not ch: return None
    if len(ch) == 1: return [ch[0][0], ch[0][1] + bl, ch[0][2]]
    return [lab, bl, ch]


def nwk(n, top=True):
    s = ('(' + ','.join(nwk(c, False) for c in n[2]) + ')' if n[2] else '') + n[0]
    return s + ';' if top else s + ':%.6g' % n[1]


def walk(n, par, tips, up=None):
    lab = cl(n[0]) or f'_u{len(par)}'; par[lab] = up
    tips[lab] = {lab} if not n[2] else set().union(*[walk(c, par, tips, lab) for c in n[2]])
    return tips[lab]


def rerun(og):
    aln = fasta(f'{AL}/{og}_species.fasta')
    tp = f'{W}/{og}_tree_pruned.nwk'; open(tp, 'w').write(nwk(prune(parse(open(TREE).read()), set(aln))) + '\n')
    log = open(f'{AR}/{og}_anc.raxml.log').read()
    cmd = next(l for l in log.splitlines() if l.startswith('raxml-ng --ancestral')).split()
    for i, a in enumerate(cmd[:-1]):
        if a == '--msa': cmd[i + 1] = f'{AL}/{og}_species.fasta'
        if a == '--tree': cmd[i + 1] = tp
        if a in ('--pref', '--prefix'): cmd[i + 1] = f'{W}/{og}_anc'
    cmd = ['nice', '-n', '19', H + '/mambaforge/bin/raxml-ng'] + cmd[1:] + ['--redo']
    print('rerun:', ' '.join(c.replace(H, '~') for c in cmd[3:]))
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    return f'{W}/{og}_anc.raxml'


def load(og):
    pre = f'{AR}/{og}_anc.raxml' if os.path.exists(f'{AR}/{og}_anc.raxml.ancestralStates') else rerun(og)
    aln = fasta(f'{AL}/{og}_species.fasta')
    par, tips = {}, {}; walk(parse(open(pre + '.ancestralTree').read()), par, tips)
    st = dict(l.rstrip('\n').split('\t')[:2] for l in open(pre + '.ancestralStates') if '\t' in l)
    pp = {}
    with open(pre + '.ancestralProbs') as f:
        next(f)
        for l in f:
            x = l.rstrip('\n').split('\t'); pp[(x[0], int(x[1]))] = max(map(float, x[3:]))
    return pre, aln, par, tips, st, pp


def state(n, c, aln, st, pp):
    """residue and marginal posterior at alignment column c (0-based); tips are observed (PP 1)."""
    if n in aln: return aln[n][c], 1.0
    return st[n][c], pp.get((n, c + 1), float('nan'))


def stems(par, tips):
    out = {}
    for u, ts in UNITS.items():
        ts = set(ts) & set(tips)
        m = min((n for n, s in tips.items() if ts <= s), key=lambda n: len(tips[n]))
        out[u] = (par[m], m, len(tips[m]) == len(ts))   # (parent node, stem node, monophyletic)
    return out


def upos(seq, c): return None if seq[c] in '-.' else len(seq[:c + 1].replace('-', '').replace('.', ''))


def colmap(seq, sites):
    pos = [i for i, ch in enumerate(seq) if ch not in '-.']
    return [pos[s - 1] if 0 < s <= len(pos) else None for s in sites]


def scan(aln, st, pp, S, comp):
    """stem-branch changes shared (same derived residue) by >= 2 freshwater units (Zhang & Kumar 1997 logic)."""
    rows, L = [], len(next(iter(aln.values())))
    for c in range(L):
        ch = {}
        for u, (a, d, mono) in S.items():
            if all(aln[t][c] in '-X?' for t in UNITS[u] if t in aln): continue
            ra, pa = state(a, c, aln, st, pp); rd, pd = state(d, c, aln, st, pp)
            if ra != rd and rd not in '-X?' and ra not in '-X?': ch[u] = (ra, pa, rd, pd)
        by = collections.defaultdict(list)
        for u, v in ch.items(): by[v[2]].append(u)
        for res, us in by.items():
            if len(us) < 2: continue
            na = len({ch[u][0] for u in us})
            kind = 'parallel' if na == 1 else ('convergent' if na == len(us) else 'mixed')
            rows.append([c + 1, res, len(us), kind,
                         ';'.join(f'{u}:{ch[u][0]}({ch[u][1]:.3f})->{ch[u][2]}({ch[u][3]:.3f})' for u in us),
                         int(all(ch[u][1] >= PPMIN and ch[u][3] >= PPMIN for u in us)),
                         sum(aln[t][c] == res for t in comp), len(comp)])
    return rows


def main():
    os.makedirs(W, exist_ok=True)
    allrows, figrows, summ = [], [], {}
    for g, (og, sites, bg, dv) in FIG.items():
        pre, aln, par, tips, st, pp = load(og)
        L = len(next(iter(aln.values())))
        assert all(len(s) == L for s in aln.values()) and all(len(s) == L for s in st.values()), og
        top = [cl(c[0]) for c in parse(open(pre + '.ancestralTree').read())[2]]
        S = stems(par, tips)
        comp = sorted(t for t in aln if t not in FW and t != 'Lgi')
        maj = []
        for c in range(L):
            cnt = collections.Counter(aln[t][c] for t in comp if aln[t][c] not in '-X?')
            maj.append(cnt.most_common(1)[0] if cnt else ('-', 0))
        def score(cols):
            sb = sum(c is not None and maj[c][0] == b for c, b in zip(cols, bg))
            sd = sum(c is not None and aln[t][c] == d for c, d in zip(cols, dv) for t in ('Dpo', 'Dro', 'Cfl') if t in aln)
            return sb, sd
        systems = {'alignment column': [s - 1 if s <= L else None for s in sites]}
        for t in sorted(aln): systems[f'{t} ungapped'] = colmap(aln[t], sites)
        ranked = sorted(systems, key=lambda k: -(score(systems[k])[0] + score(systems[k])[1] / 3))
        best = ranked[0]
        # independent of numbering: columns anywhere with comparator majority = figure background and Dpo/Dro/Cfl = figure derived
        anyhit = [sum(maj[c][0] == b and all(aln.get(t, '-')[c] == d for t in ('Dpo', 'Dro', 'Cfl')) for c in range(L))
                  for b, d in zip(bg, dv)]
        rows = scan(aln, st, pp, S, comp)
        for r in rows: allrows.append([g, og, r[0], upos(aln['Mye'], r[0] - 1) if 'Mye' in aln else ''] + r[1:])
        for s, c, b, d in zip(sites, systems[best], bg, dv):
            if c is None: figrows.append([g, s, '', b, d] + [''] * 4); continue
            sh = [r for r in rows if r[0] == c + 1]
            figrows.append([g, s, c + 1, b, d, f'{maj[c][0]}({maj[c][1]}/{len(comp)})',
                            ''.join(aln.get(t, '?')[c] for t in ('Sol', 'Pst', 'Awo', 'Ude', 'Dpo', 'Dro', 'Cfl', 'Lfo')),
                            ';'.join(f'{u}:{state(a, c, aln, st, pp)[0]}->{state(dd, c, aln, st, pp)[0]}' for u, (a, dd, m) in S.items()),
                            '|'.join(f'{x[1]}:{x[2]}units:{x[3]}' for x in sh) or 'none'])
        summ[g] = dict(og=og, prefix=pre.replace(H, '~'), ntaxa=len(aln), ncol=L, root_children=top,
                       stems={u: [a, d, m] for u, (a, d, m) in S.items()}, best_system=best,
                       top_systems=[(k, score(systems[k])) for k in ranked[:4]], nsites=len(sites),
                       anyhit_per_site=anyhit, shared_rows=len(rows), shared_pp08=sum(r[5] for r in rows),
                       shared_absent_comp=sum(r[6] == 0 for r in rows))
    with open(W + '/fig2f_shared_changes.tsv', 'w') as f:
        f.write('gene\torthogroup\taln_col\tMye_pos\tderived\tn_units\tkind\tchanges\tall_pp_ge_0.8\tn_comparators_with_derived\tn_comparators\n')
        for r in allrows: f.write('\t'.join(map(str, r)) + '\n')
    with open(W + '/fig2f_figcols.tsv', 'w') as f:
        f.write('gene\tfig_site\taln_col\tfig_background\tfig_derived\tcomparator_majority\tFW_Sol_Pst_Awo_Ude_Dpo_Dro_Cfl_Lfo\tstem_states\tshared_change_here\n')
        for r in figrows: f.write('\t'.join(map(str, r)) + '\n')
    json.dump(summ, open(W + '/fig2f_summary.json', 'w'), indent=1)
    for g, s in summ.items():
        print(g, s['og'], 'taxa', s['ntaxa'], 'cols', s['ncol'], 'root', s['root_children'][:3])
        print('  stems', s['stems'])
        print('  top systems', s['top_systems'][:3], 'max', s['nsites'], 3 * s['nsites'])
        print('  anyhit per fig site', s['anyhit_per_site'])
        print('  shared rows', s['shared_rows'], 'pp08', s['shared_pp08'], 'absent_comp', s['shared_absent_comp'])


if __name__ == '__main__':
    main()
