#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
V58 — gera todas as figuras do manuscrito V48 numa rodada.

    cd /home/juan && python3 V58_gera_figuras.py

Nao ha nada para editar.  O script descobre sozinho os arquivos de pegada, o
arquivo de SASA, a montagem do 3VBS, a cadeia de VP1 e o residuo 22, e imprime
cada escolha que fez para voce conferir.

NAO altera nem sobrescreve nenhum original.  NAO executa simulacao.  NAO
recalcula o descritor: le os valores que os seus scripts ja gravaram.  Escreve
apenas em /home/juan/V58_figuras/.

Uma figura por arquivo, SVG vetorial + PNG 600 dpi, rotulos em ingles, sem
bloco de texto dentro da figura.  Se um insumo faltar, a figura e pulada com
aviso: nenhuma figura sai desenhada com dado parcial em silencio.
"""

import os, sys, csv, glob, random, subprocess, shutil

try:
    import numpy as np
except ImportError:
    sys.exit("numpy ausente: pip install --user numpy")
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    sys.exit("matplotlib ausente: pip install --user matplotlib")
try:
    from scipy.stats import spearmanr
except ImportError:
    spearmanr = None

BASE    = "/home/juan"
RAIZES  = [f"{BASE}/V49_pegadas", f"{BASE}/abertos", f"{BASE}/paper_artifacts",
           f"{BASE}/figuras_V28", f"{BASE}/estruturas", BASE]
SAIDA   = f"{BASE}/V58_figuras"
PISO    = 106
SEMENTE = 20261116
os.makedirs(SAIDA, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "legend.frameon": False, "figure.dpi": 600,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
    "svg.fonttype": "none",
})
CINZA, ESCURO, DESTAQUE = "#9a9a9a", "#2b2b2b", "#b2182b"
LARG = 3.35

def salva(fig, nome):
    for ext in ("svg", "png"):
        fig.savefig(os.path.join(SAIDA, f"{nome}.{ext}"), format=ext)
    plt.close(fig)
    print(f"    -> {nome}.svg + .png")

def le(p):
    with open(p, newline="", encoding="utf-8", errors="replace") as f:
        d = "\t" if p.endswith(".tsv") else ","
        r = list(csv.DictReader(f, delimiter=d))
    return ([c.strip() for c in r[0].keys()] if r else []), r

def col(cs, alvos):
    for c in cs:
        if c.strip().lower() in alvos:
            return c
    return None

def num(v):
    try:
        return float(str(v).replace(",", "."))
    except (TypeError, ValueError):
        return None

def verdade(v):
    return str(v).strip().lower() in {"1", "true", "sim", "yes", "t", "s", "verdadeiro"}

print("=" * 74)
print("V58 — figuras do manuscrito V48")
print("=" * 74)

# ========================================================== 1. PEGADAS
print("\n[1] localizando pegadas")
peg, origem = {}, {}
for r in RAIZES:
    for p in sorted(glob.glob(os.path.join(r, "**", "*_pegada*.tsv"), recursive=True)):
        base = os.path.basename(p)
        pdb = base.split("_")[0].upper()
        if pdb in peg and "_V49" not in base:
            continue
        cs, rs = le(p)
        cb = col(cs, {"blocking", "bloqueio", "blocking_atoms"})
        ck = col(cs, {"shell_pct", "casca_pct", "percentil_radial"})
        cf = col(cs, {"in_footprint", "pegada", "contato"})
        cn = col(cs, {"resseq", "residuo", "resi", "num", "residue"})
        if not (cb and cf):
            continue
        b, k, f_, n_ = [], [], [], []
        for d in rs:
            vb = num(d.get(cb))
            if vb is None:
                continue
            b.append(vb)
            k.append(num(d.get(ck)) if ck else float("nan"))
            f_.append(verdade(d.get(cf)))
            n_.append(num(d.get(cn)) if cn else None)
        if b:
            peg[pdb] = (np.array(b), np.array(k, float), np.array(f_), n_)
            origem[pdb] = p

if not peg:
    print("    nenhum *_pegada*.tsv encontrado em:")
    for r in RAIZES:
        print("      ", r)
    sys.exit(1)

print(f"    {len(peg)} complexos")
for k_ in sorted(peg):
    nv = "V49" if "_V49" in origem[k_] else "   "
    print(f"      {k_} {nv}  {len(peg[k_][0]):5d} residuos   {origem[k_]}")

codigos = sorted(peg)
emb = codigos[:]
random.Random(SEMENTE).shuffle(emb)
CAL, RET = set(emb[:len(emb)//2]), set(emb[len(emb)//2:])
print(f"\n    divisao semente {SEMENTE} (ordem alfabetica embaralhada)")
print(f"      calibracao: {', '.join(sorted(CAL))}")
print(f"      retido:     {', '.join(sorted(RET))}")
print("      CONFIRA contra a divisao do V49 antes de publicar as figuras 1 e 2")

def junta(cods):
    return (np.concatenate([peg[c][0] for c in cods]),
            np.concatenate([peg[c][1] for c in cods]),
            np.concatenate([peg[c][2] for c in cods]))

# ========================================================== FIG 1
print("\n[FIG 1] distribuicao da ocupacao, metade de calibracao")
try:
    b, k, f = junta(sorted(CAL))
    bf, bn = b[f], b[~f]
    print(f"    contato {len(bf)}  nao-contato {len(bn)}")
    print(f"    contato: mediana {np.median(bf):.0f}  p95 {np.percentile(bf,95):.0f}")
    topo = float(np.percentile(b, 99.5))
    bins = np.linspace(0, topo, 60)
    fig, ax = plt.subplots(figsize=(LARG, 2.3))
    ax.hist(bn, bins=bins, color=CINZA, alpha=.85, density=True, linewidth=0,
            label="Non-contact residues")
    ax.hist(bf, bins=bins, histtype="step", color=ESCURO, linewidth=1.0,
            density=True, label="Observed contact residues")
    ax.axvline(PISO, color=DESTAQUE, lw=.9, ls="--")
    ax.set_xlabel("Outward-cone atomic occupancy (atoms)")
    ax.set_ylabel("Density")
    ax.set_xlim(0, topo); ax.legend(loc="upper right")
    salva(fig, "Fig1_occupancy_distribution_calibration")
except Exception as e:
    print(f"    PULADA: {e}")

# ========================================================== FIG 2
print("\n[FIG 2] varredura do limiar, metade retida")
try:
    b, k, f = junta(sorted(RET))
    print(f"    contato {int(f.sum())}  nao-contato {int((~f).sum())}")
    ts = np.arange(0, int(np.percentile(b, 99.5)) + 1)
    sens = np.array([(b[f] <= t).mean() for t in ts]) * 100
    espe = np.array([(b[~f] > t).mean() for t in ts]) * 100
    i = int(np.argmin(abs(ts - PISO)))
    print(f"    em {PISO}: sensibilidade {sens[i]:.1f}%  especificidade {espe[i]:.1f}%")
    fig, ax = plt.subplots(figsize=(LARG, 2.3))
    ax.plot(ts, sens, color=ESCURO, lw=1.0, label="Contact retention")
    ax.plot(ts, espe, color=CINZA,  lw=1.0, label="Non-contact rejection")
    ax.axvline(PISO, color=DESTAQUE, lw=.9, ls="--")
    ax.set_xlabel("Occupancy ceiling (atoms)")
    ax.set_ylabel("Percentage of residues")
    ax.set_xlim(ts.min(), ts.max()); ax.set_ylim(0, 100)
    ax.legend(loc="center right")
    salva(fig, "Fig2_threshold_sweep_heldout")
except Exception as e:
    print(f"    PULADA: {e}")

# ========================================================== FIG 3
print("\n[FIG 3] ocupacao x percentil radial")
try:
    b, k, f = junta(codigos)
    m = np.isfinite(k) & np.isfinite(b)
    b, k, f = b[m], k[m], f[m]
    print(f"    residuos com os dois valores: {len(b)}")
    if spearmanr and len(b) > 2:
        print(f"    Spearman: {spearmanr(k, b).statistic:+.3f}")
    fig, ax = plt.subplots(figsize=(LARG, 2.5))
    ax.scatter(k[~f], b[~f], s=1.2, c=CINZA, lw=0, alpha=.5, label="Non-contact residues")
    ax.scatter(k[f],  b[f],  s=1.6, c=ESCURO, lw=0, alpha=.8, label="Observed contact residues")
    ax.axhline(PISO, color=DESTAQUE, lw=.9, ls="--")
    ax.set_xlabel("Radial Cα percentile within the shell")
    ax.set_ylabel("Outward-cone atomic occupancy (atoms)")
    ax.set_xlim(0, 100); ax.legend(loc="upper right", markerscale=4)
    salva(fig, "Fig3_occupancy_vs_radial_percentile")
except Exception as e:
    print(f"    PULADA: {e}")

# ========================================================== FIG 4
print("\n[FIG 4] retencao de SASA x percentil radial")
p_sasa = None
for r in RAIZES:
    for p in sorted(glob.glob(os.path.join(r, "**", "*.tsv"), recursive=True)
                  + glob.glob(os.path.join(r, "**", "*.csv"), recursive=True)):
        bn_ = os.path.basename(p).lower()
        if not any(t in bn_ for t in ("sasa", "retenc", "v28_15", "acess")):
            continue
        try:
            with open(p, encoding="utf-8", errors="replace") as f_:
                cab = f_.readline().lower()
        except OSError:
            continue
        if any(t in cab for t in ("retenc", "retain", "fracao", "fraction")):
            p_sasa = p
            break
    if p_sasa:
        break

if not p_sasa:
    print("    PULADA: nao achei tabela de SASA por residuo com coluna de retencao")
else:
    try:
        print(f"    arquivo: {p_sasa}")
        cs, rs = le(p_sasa)
        print(f"    colunas: {cs}")
        cr = col(cs, {"retencao", "retained", "fracao", "fraction", "ret"})
        ck = col(cs, {"shell_pct", "casca_pct", "percentil_radial", "percentil"})
        cn = col(cs, {"resseq", "residuo", "resi", "num", "residue"})
        if not (cr and ck):
            raise RuntimeError(f"coluna faltando: retencao={cr} percentil={ck}")
        R, K, N = [], [], []
        for d in rs:
            a, c = num(d.get(cr)), num(d.get(ck))
            if a is None or c is None:
                continue
            R.append(a); K.append(c); N.append(num(d.get(cn)) if cn else None)
        R, K = np.array(R), np.array(K)
        print(f"    residuos: {len(R)}")
        if spearmanr and len(R) > 2:
            print(f"    Spearman: {spearmanr(K, R).statistic:+.3f}")
        fig, ax = plt.subplots(figsize=(LARG, 2.5))
        ax.scatter(K, R, s=3, c=CINZA, lw=0, alpha=.7)
        j = next((i for i, v in enumerate(N) if v == 22), None)
        if j is not None:
            ax.scatter([K[j]], [R[j]], s=22, facecolors="none",
                       edgecolors=DESTAQUE, linewidths=1.0)
            print(f"    residuo 22: percentil {K[j]:.2f}  retencao {R[j]:.4f}")
        else:
            print("    residuo 22 nao localizado; nao foi marcado")
        ax.set_xlabel("Radial Cα percentile within the shell")
        ax.set_ylabel("Retained solvent-accessible area (fraction)")
        ax.set_xlim(0, 100); ax.set_ylim(0, 1.05)
        salva(fig, "Fig4_sasa_retention_vs_radial")
    except Exception as e:
        print(f"    PULADA: {e}")

# ========================================================== FIG 5
print("\n[FIG 5] ocupacao do sitio 22 nos depositos EV-A71")
EV = {"3VBF", "3VBH", "3VBO", "3VBR", "3VBS", "3VBU", "6I2K"}
alvo = {}
for pdb in sorted(set(peg) & EV):
    b, k, f, n_ = peg[pdb]
    j = next((i for i, v in enumerate(n_) if v == 22), None)
    if j is not None:
        alvo[pdb] = (b[j], k[j])
if not alvo:
    print("    PULADA: residuo 22 nao achado nas pegadas dos depositos EV-A71")
else:
    try:
        for p_, (vb, vk) in sorted(alvo.items()):
            print(f"    {p_}: ocupacao {vb:.0f}  percentil {vk:.2f}" if np.isfinite(vk)
                  else f"    {p_}: ocupacao {vb:.0f}")
        nomes = sorted(alvo); vals = [alvo[n][0] for n in nomes]
        fig, ax = plt.subplots(figsize=(LARG, 2.1))
        ax.bar(nomes, vals, color=CINZA, width=.55, linewidth=0)
        ax.axhline(PISO, color=DESTAQUE, lw=.9, ls="--")
        ax.set_ylabel("Outward-cone atomic occupancy (atoms)")
        ax.set_xlabel("Deposited structure")
        ax.set_ylim(0, max(vals) * 1.12)
        salva(fig, "Fig5_site22_occupancy_across_deposits")
    except Exception as e:
        print(f"    PULADA: {e}")

# ========================================================== 6. CAPSIDEO
print("\n[2] figuras de capsideo")
pymol_bin = shutil.which("pymol")
if not pymol_bin:
    print("    PyMOL nao encontrado no PATH — parte de capsideo pulada.")
    print("    os cinco graficos acima nao dependem dele.")
else:
    print(f"    PyMOL: {pymol_bin}")
    # localiza a montagem do 3VBS: prefere arquivo com 'montagem'/'assembly'
    cand = []
    for r in RAIZES:
        for ext in ("cif", "pdb", "cif.gz", "pdb.gz", "ent"):
            cand += glob.glob(os.path.join(r, "**", f"*3[Vv][Bb][Ss]*.{ext}"), recursive=True)
    def pontua(p):
        b = os.path.basename(p).lower()
        s = 0
        for t, v in (("montagem", 40), ("assembly", 40), ("240", 20),
                     ("expand", 15), ("bio", 10)):
            if t in b:
                s += v
        return (s, os.path.getsize(p))
    cand = sorted(set(cand), key=pontua, reverse=True)
    if not cand:
        print("    nenhum arquivo 3VBS encontrado — parte de capsideo pulada")
    else:
        print(f"    candidatos ({len(cand)}), em ordem de preferencia:")
        for p in cand[:6]:
            print(f"      {os.path.getsize(p)/1e6:8.1f} MB  {p}")
        estrut = cand[0]
        print(f"    usando: {estrut}")

        # identifica a cadeia de VP1 e confirma o residuo 22, sem chutar
        cadeia = None
        try:
            import gemmi
            st = gemmi.read_structure(estrut)
            st.setup_entities()
            mod = st[0]
            inv = []
            for ch in mod:
                aas = [r for r in ch if r.name not in ("HOH", "WAT")]
                if not aas:
                    continue
                r22 = next((r for r in aas if r.seqid.num == 22), None)
                inv.append((ch.name, len(aas), r22.name if r22 else None))
            print(f"    cadeias no modelo: {len(inv)}")
            for nm, n, r22 in inv[:8]:
                print(f"      cadeia {nm:<3} {n:4d} residuos   residuo 22: {r22 or '-'}")
            if len(inv) > 8:
                print(f"      ... e mais {len(inv)-8}")
            # regra declarada: VP1 e a cadeia com 270-310 residuos que tem o 22
            alvos = [c for c in inv if c[2] and 270 <= c[1] <= 310]
            if alvos:
                cadeia = alvos[0][0]
                print(f"    VP1 por regra (270-310 residuos e tem o 22): cadeia "
                      f"{cadeia}, residuo 22 = {alvos[0][2]}")
            else:
                comr22 = [c for c in inv if c[2]]
                if comr22:
                    cadeia = max(comr22, key=lambda c: c[1])[0]
                    print(f"    nenhuma cadeia em 270-310; usando a mais longa com "
                          f"residuo 22: cadeia {cadeia}")
                    print("    CONFIRA se essa e mesmo VP1 antes de publicar")
        except ImportError:
            print("    gemmi ausente; nao consegui identificar a cadeia")
        except Exception as e:
            print(f"    falha ao inspecionar a estrutura: {e}")

        if not cadeia:
            print("    cadeia de VP1 indeterminada — capsideo pulado")
        else:
            pml = os.path.join(SAIDA, "_capsideo_gerado.pml")
            with open(pml, "w") as f:
                f.write(f"""
set assembly, 1
set orthoscopic, 1
set antialias, 2
set ray_shadows, 0
set depth_cue, 0
set specular, 0.15
set ambient, 0.18
set direct, 0.55
set ray_opaque_background, 0
bg_color white
load {estrut}, capsid
select vp1,    capsid and polymer and chain {cadeia}
select site22, capsid and polymer and chain {cadeia} and resi 22
python
from pymol import cmd
print("  atomos capsid =", cmd.count_atoms("capsid"))
print("  atomos vp1    =", cmd.count_atoms("vp1"))
print("  atomos site22 =", cmd.count_atoms("site22"))
python end
hide everything
show spheres, capsid
color grey80, capsid
color firebrick, site22
set sphere_scale, 1.15, site22
orient capsid
zoom capsid, 4
png {SAIDA}/Fig6_capsid_outer_surface.png, width=3000, height=3000, dpi=600, ray=1
clip slab, 60
png {SAIDA}/Fig7_capsid_equatorial_section.png, width=3000, height=3000, dpi=600, ray=1
clip slab, 1000
hide everything
show cartoon, vp1
color grey70, vp1
show spheres, site22
color firebrick, site22
orient vp1
zoom vp1, 3
png {SAIDA}/Fig8_protomer_site22.png, width=3000, height=2400, dpi=600, ray=1
""")
            print(f"    pml gerado: {pml}")
            print("    renderizando (pode levar varios minutos na particula inteira)")
            try:
                r = subprocess.run([pymol_bin, "-cq", pml], capture_output=True,
                                   text=True, timeout=5400)
                for ln in (r.stdout or "").splitlines():
                    if ln.strip():
                        print("      " + ln)
                if r.returncode != 0:
                    print(f"      PyMOL saiu com codigo {r.returncode}")
                    for ln in (r.stderr or "").splitlines()[-15:]:
                        print("      ! " + ln)
            except subprocess.TimeoutExpired:
                print("      estourou 90 min; rode o pml a mao se precisar")

print("\n" + "=" * 74)
for p in sorted(glob.glob(os.path.join(SAIDA, "*.svg"))
              + glob.glob(os.path.join(SAIDA, "*.png"))):
    print(f"  {os.path.getsize(p)/1024:8.0f} kB  {p}")
print(f"\nfiguras em {SAIDA}")
print("nenhum arquivo original foi alterado")
print("=" * 74)
