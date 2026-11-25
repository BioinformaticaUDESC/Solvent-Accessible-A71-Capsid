#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
V63 — SASA do VP1 de 3VBS sem passar pelo formato PDB.

    cd /home/juan && python3 V63_sasa_3vbs.py

O V61 falhou porque convertia a estrutura para texto PDB, que so aceita
cadeia de UM caractere: a montagem tem A1 ate A60, D1 ate D60, e `A10`
estourou.  A cadeia isolada deu zero pelo mesmo motivo — `A1` nao sobrevive
a conversao e o filtro nao casava.

Aqui o freesasa recebe coordenadas e raios direto (calcCoord), sem texto
PDB no meio.  Os raios vem do classificador padrao do proprio freesasa
(ProtOr), com tabela por elemento como reserva.

Para o estado montado so entram os atomos que podem ocluir a cadeia de
referencia: a propria cadeia mais tudo dentro de 12 A de qualquer atomo
dela.  Com sonda de 1,4 A nenhum atomo alem disso altera a area da
referencia, entao o recorte e exato e nao aproximacao.

Grava /home/juan/V63_SASA_vs_RADIAL_3VBS.tsv — nome novo, nada sobrescrito.
Ao final compara com o que o manuscrito afirma e avisa se divergir.
"""

import os, sys, glob
import numpy as np

PUBLICADO = {"ret22": 0.9118, "rho": 0.329, "n": 293,
             "perc22": 1.63, "iso22": 147.67, "mont22": 134.64}
CORTE_VIZ = 12.0

BASE   = "/home/juan"
RAIZES = [f"{BASE}/figuras_V28", f"{BASE}/abertos", f"{BASE}/paper_artifacts",
          f"{BASE}/estruturas", BASE]
SAIDA  = f"{BASE}/V58_figuras"
TABELA = f"{BASE}/V63_SASA_vs_RADIAL_3VBS.tsv"
os.makedirs(SAIDA, exist_ok=True)

try:
    import gemmi
except ImportError:
    sys.exit("gemmi ausente: pip install --user gemmi")
try:
    import freesasa
except ImportError:
    sys.exit("freesasa ausente: pip install --user freesasa")
try:
    from scipy.spatial import cKDTree
except ImportError:
    sys.exit("scipy ausente: pip install --user scipy")

print("=" * 74); print("V63 — SASA do VP1 de 3VBS"); print("=" * 74)

# ------------------------------------------------------ estrutura
cand = []
for r in RAIZES:
    for e in ("cif","pdb","cif.gz","pdb.gz","ent"):
        cand += glob.glob(os.path.join(r,"**",f"*3[Vv][Bb][Ss]*.{e}"), recursive=True)
def pontua(p):
    b = os.path.basename(p).lower()
    return (sum(v for t,v in (("assembly",40),("montagem",40),("240",20),
                              ("expand",15)) if t in b), os.path.getsize(p))
cand = sorted(set(cand), key=pontua, reverse=True)
if not cand: sys.exit("nenhum 3VBS encontrado")
estrut = cand[0]
print(f"\nestrutura: {estrut}")

st = gemmi.read_structure(estrut)
st.setup_entities(); st.remove_ligands_and_waters(); st.remove_hydrogens()
mod = st[0]

_aa = {}
def eh_aa(nome):
    if nome not in _aa:
        i = gemmi.find_tabulated_residue(nome)
        _aa[nome] = bool(i and i.is_amino_acid())
    return _aa[nome]

inv = []
for ch in mod:
    aas = [r for r in ch if eh_aa(r.name)]
    if aas:
        inv.append((ch.name, len(aas), any(r.seqid.num == 22 for r in aas)))
vp1 = [c for c in inv if 270 <= c[1] <= 310 and c[2]]
print(f"cadeias: {len(inv)}   VP1: {len(vp1)}")
if not vp1: sys.exit("nenhuma cadeia VP1 identificada")
REF = vp1[0][0]
print(f"cadeia de referencia: {REF}  ({vp1[0][1]} residuos)")

# ------------------------------------------------------ atomos
RAIO_EL = {"C":1.88,"N":1.64,"O":1.46,"S":1.77,"P":1.80,"SE":1.90}
clf = freesasa.Classifier()
def raio(resname, atname, elem):
    try:
        v = clf.radius(resname, atname)
        if v and v > 0: return float(v)
    except Exception:
        pass
    return RAIO_EL.get((elem or "").upper().strip(), 1.80)

xyz, rad, dono = [], [], []           # dono = (cadeia, resseq) ou None
for ch in mod:
    marca = (ch.name == REF)
    for r in ch:
        if not eh_aa(r.name): continue
        chave = (ch.name, r.seqid.num) if marca else None
        for a in r:
            xyz.append((a.pos.x, a.pos.y, a.pos.z))
            rad.append(raio(r.name, a.name, a.element.name))
            dono.append(chave)
xyz = np.asarray(xyz, dtype=float)
rad = np.asarray(rad, dtype=float)
eh_ref = np.array([d is not None for d in dono])
print(f"atomos de aminoacido: {len(xyz)}   na cadeia {REF}: {int(eh_ref.sum())}")

# ------------------------------------------------------ centroide / percentil
cas = []
for ch in mod:
    for r in ch:
        if not eh_aa(r.name): continue
        a = r.find_atom("CA", "*")
        if a: cas.append((a.pos.x, a.pos.y, a.pos.z))
cas = np.asarray(cas, dtype=float)
centro = cas.mean(axis=0)
ordenados = np.sort(np.linalg.norm(cas - centro, axis=1))
print(f"centroide: ({centro[0]:.1f}, {centro[1]:.1f}, {centro[2]:.1f})   "
      f"Ca {len(cas)}   raio {ordenados[0]:.1f} a {ordenados[-1]:.1f} A")
def percentil(x):
    return 100.0 * np.searchsorted(ordenados, x) / len(ordenados)

# ------------------------------------------------------ SASA
def sasa(coords, radios):
    r = freesasa.calcCoord(coords.ravel().tolist(), radios.tolist())
    return np.array([r.atomArea(i) for i in range(len(radios))])

print("\n[1] cadeia isolada")
A_iso_at = sasa(xyz[eh_ref], rad[eh_ref])
print(f"    atomos {int(eh_ref.sum())}   area total {A_iso_at.sum():.1f} A^2")

print(f"[2] cadeia na particula (vizinhanca de {CORTE_VIZ:.0f} A)")
viz = set(np.flatnonzero(eh_ref).tolist())
arv = cKDTree(xyz)
for j in arv.query_ball_point(xyz[eh_ref], CORTE_VIZ):
    viz.update(j)
idx = np.array(sorted(viz))
print(f"    atomos no recorte: {len(idx)} de {len(xyz)}")
A_mont_at_todos = sasa(xyz[idx], rad[idx])
pos_de = {g: k for k, g in enumerate(idx)}
A_mont_at = np.array([A_mont_at_todos[pos_de[g]] for g in np.flatnonzero(eh_ref)])
print(f"    area total da cadeia na particula: {A_mont_at.sum():.1f} A^2")
print(f"    enterrado pela montagem: "
      f"{100*(1 - A_mont_at.sum()/A_iso_at.sum()):.1f}%")

# ------------------------------------------------------ por residuo
donos_ref = [d for d in dono if d is not None]
iso_res, mont_res = {}, {}
for k, g in enumerate(donos_ref):
    iso_res[g]  = iso_res.get(g, 0.0)  + A_iso_at[k]
    mont_res[g] = mont_res.get(g, 0.0) + A_mont_at[k]

cad = next(ch for ch in mod if ch.name == REF)
linhas = []
for r in cad:
    if not eh_aa(r.name): continue
    a = r.find_atom("CA", "*")
    if not a: continue
    g = (REF, r.seqid.num)
    if g not in iso_res: continue
    iso, mont = iso_res[g], mont_res[g]
    rr = float(np.linalg.norm(np.array([a.pos.x,a.pos.y,a.pos.z]) - centro))
    linhas.append({"residue": r.seqid.num, "resname": r.name, "chain": REF,
                   "sasa_isolada": round(iso,3), "sasa_montada": round(mont,3),
                   "retencao": round(mont/iso,6) if iso > 0 else "",
                   "raio_ca": round(rr,3), "shell_pct": round(percentil(rr),4)})

cols = ["residue","resname","chain","sasa_isolada","sasa_montada",
        "retencao","raio_ca","shell_pct"]
with open(TABELA,"w",encoding="utf-8") as f:
    f.write("\t".join(cols)+"\n")
    for d in linhas:
        f.write("\t".join(str(d[c]) for c in cols)+"\n")
print(f"\nresiduos na tabela: {len(linhas)}")
print(f"-> {TABELA}")

# ------------------------------------------------------ verificacao
print("\n" + "-"*74)
print("VERIFICACAO contra o manuscrito")
print("-"*74)
val = [d for d in linhas if d["retencao"] != ""]
R = np.array([d["retencao"]  for d in val], float)
K = np.array([d["shell_pct"] for d in val], float)
N = [d["residue"] for d in val]
nz = R > 0
print(f"  residuos com area isolada > 0   aqui {int(nz.sum()):>6}   texto {PUBLICADO['n']:>6}")
from scipy.stats import spearmanr
rho = spearmanr(K, R).statistic
print(f"  Spearman retencao x percentil   aqui {rho:+6.3f}   texto {PUBLICADO['rho']:+6.3f}")
j = next((i for i,v in enumerate(N) if v == 22), None)
if j is None:
    print("  residuo 22 ausente da tabela")
else:
    d22 = next(d for d in val if d["residue"] == 22)
    print(f"  residuo 22  {d22['resname']}")
    print(f"    area isolada   aqui {d22['sasa_isolada']:8.2f}   texto {PUBLICADO['iso22']:8.2f}")
    print(f"    area montada   aqui {d22['sasa_montada']:8.2f}   texto {PUBLICADO['mont22']:8.2f}")
    print(f"    retencao       aqui {R[j]:8.4f}   texto {PUBLICADO['ret22']:8.4f}")
    print(f"    percentil      aqui {K[j]:8.2f}   texto {PUBLICADO['perc22']:8.2f}")

aviso = []
if j is not None and abs(R[j]-PUBLICADO["ret22"]) > 0.01: aviso.append("retencao do 22")
if abs(rho-PUBLICADO["rho"]) > 0.03: aviso.append("Spearman")
if abs(int(nz.sum())-PUBLICADO["n"]) > 3: aviso.append("n de residuos")
if aviso:
    print(f"\n  >>> DIVERGE em: {', '.join(aviso)}")
    print("      Nao mexi no manuscrito.  Me manda esta saida.")
else:
    print("\n  tudo dentro da tolerancia do que o texto afirma")

# decis extremos, como o texto relata
lo, hi = K <= 10, K >= 90
if lo.any() and hi.any():
    print(f"\n  decil interno (percentil <=10): n={int(lo.sum())}  "
          f"mediana da retencao {np.median(R[lo]):.3f}")
    print(f"  decil externo (percentil >=90): n={int(hi.sum())}  "
          f"mediana da retencao {np.median(R[hi]):.3f}")

# ------------------------------------------------------ Fig4
print("\n[FIG 4]")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({
    "font.family":"sans-serif","font.sans-serif":["DejaVu Sans","Arial"],
    "font.size":8,"axes.labelsize":8,"xtick.labelsize":7,"ytick.labelsize":7,
    "axes.spines.top":False,"axes.spines.right":False,"axes.linewidth":.6,
    "xtick.major.width":.6,"ytick.major.width":.6,"legend.frameon":False,
    "figure.dpi":600,"svg.fonttype":"none"})
fig, ax = plt.subplots(figsize=(3.35,2.6))
ax.scatter(K, R, s=3.5, c="#9a9a9a", lw=0, alpha=.75)
if j is not None:
    ax.scatter([K[j]],[R[j]], s=28, facecolors="none",
               edgecolors="#b2182b", linewidths=1.2)
ax.set_xlabel("Radial Cα percentile within the shell")
ax.set_ylabel("Retained accessible area (fraction)")
ax.set_xlim(0,100); ax.set_ylim(0,1.05)
for e in ("svg","png"):
    fig.savefig(os.path.join(SAIDA,f"Fig4_sasa_retention_vs_radial.{e}"),
                format=e, bbox_inches="tight", pad_inches=0.04)
plt.close(fig)
print(f"  -> {SAIDA}/Fig4_sasa_retention_vs_radial.svg + .png   pontos {len(R)}")

print("\n" + "="*74)
print("nenhum arquivo original foi alterado")
print("="*74)
