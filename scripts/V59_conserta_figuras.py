#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
V59 — conserta as figuras que sairam erradas no V58 e gera a que faltou.

    cd /home/juan && python3 V59_conserta_figuras.py

O que muda em relacao ao V58:

  Fig3  rotulo do eixo y estava cortado -> rotulo mais curto e bbox explicito
  Fig4  nao foi gerada -> busca pelo nome exato que o V54 usa
        (V28_SASA_vs_RADIAL*) e, se nao achar, LISTA os candidatos que
        considerou e por que rejeitou cada um
  Fig5  so apareceu o 6I2K -> o residuo 22 dos depositos 3VB* nao esta nas
        pegadas; esta nos *perfil_radial*, na coluna `residue`, que e de onde
        o V49 leu.  Agora le de la.
  Fig6  e Fig7: o sitio nao aparecia porque eu marcava UMA cadeia num
        particulado de 240.  Agora marca o residuo 22 em TODAS as copias de
        VP1, que e o que a figura antiga do V13 mostrava.

Nao altera nem sobrescreve original.  Escreve so em /home/juan/V58_figuras/.
Nao recalcula descritor: le o que os seus scripts gravaram.
"""

import os, sys, csv, glob, subprocess, shutil

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
try:
    from scipy.stats import spearmanr
except ImportError:
    spearmanr = None

BASE   = "/home/juan"
RAIZES = [f"{BASE}/V49_pegadas", f"{BASE}/abertos", f"{BASE}/paper_artifacts",
          f"{BASE}/figuras_V28", f"{BASE}/estruturas", BASE]
SAIDA  = f"{BASE}/V58_figuras"
PISO   = 106
os.makedirs(SAIDA, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 8, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": .6, "xtick.major.width": .6, "ytick.major.width": .6,
    "legend.frameon": False, "figure.dpi": 600, "svg.fonttype": "none",
})
CINZA, ESCURO, DESTAQUE = "#9a9a9a", "#2b2b2b", "#b2182b"
LARG = 3.35

def salva(fig, nome):
    for ext in ("svg", "png"):
        fig.savefig(os.path.join(SAIDA, f"{nome}.{ext}"), format=ext,
                    bbox_inches="tight", pad_inches=0.03)
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
    return str(v).strip().lower() in {"1","true","sim","yes","t","s","verdadeiro"}

print("=" * 74)
print("V59 — conserto das figuras")
print("=" * 74)

# ============================================================ FIG 3 (rotulo)
print("\n[FIG 3] regerando com rotulo que cabe")
peg = {}
for r in RAIZES:
    for p in sorted(glob.glob(os.path.join(r, "**", "*_pegada*.tsv"), recursive=True)):
        b_ = os.path.basename(p); pdb = b_.split("_")[0].upper()
        if pdb in peg and "_V49" not in b_:
            continue
        cs, rs = le(p)
        cb = col(cs, {"blocking","bloqueio","blocking_atoms"})
        ck = col(cs, {"shell_pct","casca_pct","percentil_radial"})
        cf = col(cs, {"in_footprint","pegada","contato"})
        if not (cb and cf):
            continue
        B,K,F = [],[],[]
        for d in rs:
            vb = num(d.get(cb))
            if vb is None: continue
            B.append(vb); K.append(num(d.get(ck)) if ck else float("nan"))
            F.append(verdade(d.get(cf)))
        if B: peg[pdb] = (np.array(B), np.array(K,float), np.array(F))
print(f"    {len(peg)} complexos")
try:
    b = np.concatenate([peg[c][0] for c in peg])
    k = np.concatenate([peg[c][1] for c in peg])
    f = np.concatenate([peg[c][2] for c in peg])
    m = np.isfinite(k) & np.isfinite(b); b,k,f = b[m],k[m],f[m]
    print(f"    residuos: {len(b)}")
    if spearmanr: print(f"    Spearman: {spearmanr(k,b).statistic:+.3f}")
    fig, ax = plt.subplots(figsize=(LARG, 2.6))
    ax.scatter(k[~f], b[~f], s=1.2, c=CINZA, lw=0, alpha=.5, label="Non-contact residues")
    ax.scatter(k[f],  b[f],  s=1.6, c=ESCURO, lw=0, alpha=.8, label="Observed contact residues")
    ax.axhline(PISO, color=DESTAQUE, lw=.9, ls="--")
    ax.set_xlabel("Radial Cα percentile within the shell")
    ax.set_ylabel("Outward-cone occupancy (atoms)")
    ax.set_xlim(0, 100)
    ax.legend(loc="upper right", markerscale=4)
    salva(fig, "Fig3_occupancy_vs_radial_percentile")
except Exception as e:
    print(f"    PULADA: {e}")

# ============================================================ FIG 4 (SASA)
print("\n[FIG 4] retencao de SASA x percentil radial")
alvo, considerados = None, []
# 1) nome exato que o V54 usa
for r in RAIZES:
    g = glob.glob(os.path.join(r, "**", "V28_SASA_vs_RADIAL*"), recursive=True)
    g = [x for x in g if x.endswith((".tsv", ".csv"))]
    if g:
        alvo = sorted(g, key=len)[0]
        break
# 2) qualquer tabela com retencao e percentil
if not alvo:
    for r in RAIZES:
        for p in sorted(glob.glob(os.path.join(r, "**", "*.tsv"), recursive=True)
                      + glob.glob(os.path.join(r, "**", "*.csv"), recursive=True)):
            bn = os.path.basename(p).lower()
            if not any(t in bn for t in ("sasa","retenc","acess","radial","v28")):
                continue
            try:
                cs, _ = le(p)
            except Exception:
                continue
            cr = col(cs, {"retencao","retained","fracao","fraction","ret"})
            ck = col(cs, {"shell_pct","casca","casca_pct","percentil","percentil_radial"})
            considerados.append((p, cs, cr, ck))
            if cr and ck:
                alvo = p
                break
        if alvo:
            break

if not alvo:
    print("    NAO LOCALIZADA. candidatos que examinei e o que faltou em cada:")
    if not considerados:
        print("      nenhum arquivo com 'sasa', 'retenc', 'acess', 'radial' ou 'v28' no nome")
        print("      rode:  find /home/juan -iname '*sasa*' -o -iname '*V28*' | head -40")
    for p, cs, cr, ck in considerados[:20]:
        falta = []
        if not cr: falta.append("coluna de retencao")
        if not ck: falta.append("coluna de percentil")
        print(f"      {p}")
        print(f"         colunas: {cs}")
        print(f"         faltou: {', '.join(falta)}")
    print("    me manda essa lista que eu ajusto os nomes de coluna")
else:
    try:
        print(f"    arquivo: {alvo}")
        cs, rs = le(alvo)
        print(f"    colunas: {cs}")
        cr = col(cs, {"retencao","retained","fracao","fraction","ret"})
        ck = col(cs, {"shell_pct","casca","casca_pct","percentil","percentil_radial"})
        cn = col(cs, {"residue","resseq","residuo","resi","num"})
        ci = col(cs, {"sasa_isolada","sasa_isolated","isolada","iso_a2"})
        print(f"    usando: retencao={cr}  percentil={ck}  residuo={cn}  isolada={ci}")
        R,K,N = [],[],[]
        for d in rs:
            a, c = num(d.get(cr)), num(d.get(ck))
            if a is None or c is None: continue
            R.append(a); K.append(c); N.append(num(d.get(cn)) if cn else None)
        R,K = np.array(R), np.array(K)
        print(f"    residuos: {len(R)}")
        if spearmanr and len(R) > 2:
            print(f"    Spearman retencao x percentil: {spearmanr(K,R).statistic:+.3f}")
        fig, ax = plt.subplots(figsize=(LARG, 2.6))
        ax.scatter(K, R, s=3, c=CINZA, lw=0, alpha=.7)
        j = next((i for i,v in enumerate(N) if v == 22), None)
        if j is not None:
            ax.scatter([K[j]],[R[j]], s=24, facecolors="none",
                       edgecolors=DESTAQUE, linewidths=1.1)
            print(f"    residuo 22: percentil {K[j]:.2f}  retencao {R[j]:.4f}")
        else:
            print("    residuo 22 nao achado na coluna de numeracao; nao marcado")
        ax.set_xlabel("Radial Cα percentile within the shell")
        ax.set_ylabel("Retained accessible area (fraction)")
        ax.set_xlim(0, 100); ax.set_ylim(0, 1.05)
        salva(fig, "Fig4_sasa_retention_vs_radial")
    except Exception as e:
        print(f"    PULADA: {e}")

# ============================================================ FIG 5 (perfis)
print("\n[FIG 5] sitio 22 nos depositos — lendo dos perfis radiais")
perfis = set()
for r in RAIZES:
    perfis |= {os.path.realpath(p)
               for p in glob.glob(os.path.join(r, "**", "*perfil_radial*"), recursive=True)
               if p.endswith((".tsv", ".csv"))}
print(f"    {len(perfis)} arquivos de perfil radial")
alvo22 = {}
for p in sorted(perfis):
    cs, rs = le(p)
    cb = col(cs, {"blocking_atoms","blocking","bloqueio"})
    ck = col(cs, {"shell_pct","casca_pct"})
    cn = col(cs, {"residue","resseq","residuo","resi"})
    if not (cb and cn):
        continue
    for d in rs:
        if str(d.get(cn, "")).strip() != "22":
            continue
        vb = num(d.get(cb))
        if vb is None:
            break
        nome = os.path.basename(p)
        rot = None
        for tag in ("3VBF","3VBH","3VBO","3VBR","3VBS","3VBU","6I2K"):
            if tag.lower() in nome.lower():
                rot = tag
                break
        if rot is None:
            rot = nome.split("_")[0].upper()
        alvo22[rot] = (vb, num(d.get(ck)) if ck else None,
                       str(d.get("resname","")).strip(), nome)
        break

if not alvo22:
    print("    NAO LOCALIZADO residuo 22 em nenhum perfil radial")
    print("    perfis examinados:")
    for p in sorted(perfis)[:20]:
        print("      ", p)
else:
    for k_ in sorted(alvo22):
        vb, vk, rn, nome = alvo22[k_]
        print(f"    {k_:<6} {rn:<4} ocupacao {vb:6.0f}  "
              f"percentil {('%.2f'%vk) if vk is not None else 'n/d'}  "
              f"{'passa' if vb <= PISO else 'REPROVA'}   [{nome}]")
    try:
        nomes = sorted(alvo22); vals = [alvo22[n][0] for n in nomes]
        fig, ax = plt.subplots(figsize=(max(LARG, .62*len(nomes)+1.1), 2.3))
        ax.bar(nomes, vals, color=CINZA, width=.6, linewidth=0)
        ax.axhline(PISO, color=DESTAQUE, lw=.9, ls="--")
        ax.set_ylabel("Outward-cone occupancy (atoms)")
        ax.set_xlabel("Deposited structure")
        ax.set_ylim(0, max(vals)*1.12)
        salva(fig, "Fig5_site22_occupancy_across_deposits")
    except Exception as e:
        print(f"    PULADA: {e}")

# ============================================================ FIG 6 e 7
print("\n[CAPSIDEO] marcando o residuo 22 em TODAS as copias de VP1")
pm = shutil.which("pymol")
if not pm:
    print("    PyMOL ausente; pulado")
    sys.exit(0)

cand = []
for r in RAIZES:
    for ext in ("cif","pdb","cif.gz","pdb.gz","ent"):
        cand += glob.glob(os.path.join(r, "**", f"*3[Vv][Bb][Ss]*.{ext}"), recursive=True)
def pontua(p):
    b = os.path.basename(p).lower()
    return (sum(v for t,v in (("montagem",40),("assembly",40),("240",20),
                              ("expand",15),("bio",10)) if t in b),
            os.path.getsize(p))
cand = sorted(set(cand), key=pontua, reverse=True)
if not cand:
    sys.exit("    nenhum arquivo 3VBS encontrado")
estrut = cand[0]
print(f"    estrutura: {estrut}")

try:
    import gemmi
except ImportError:
    sys.exit("    gemmi ausente; nao consigo listar as cadeias de VP1")

st = gemmi.read_structure(estrut)
st.setup_entities()
mod = st[0]
vp1, outras = [], 0
for ch in mod:
    aas = [r for r in ch if r.name not in ("HOH","WAT")]
    if not aas:
        continue
    tem22 = any(r.seqid.num == 22 for r in aas)
    if 270 <= len(aas) <= 310 and tem22:
        vp1.append(ch.name)
    else:
        outras += 1
print(f"    cadeias no modelo: {len(vp1)+outras}   VP1 (270-310 res e tem o 22): {len(vp1)}")
if not vp1:
    sys.exit("    nenhuma cadeia VP1 identificada; nao desenho a figura")
print(f"    primeiras: {', '.join(vp1[:12])}{' ...' if len(vp1)>12 else ''}")
if len(vp1) != 60:
    print(f"    ATENCAO: esperado 60 copias de VP1 num capsideo T=1 pseudo-T=3; "
          f"achei {len(vp1)}.  Confira antes de publicar.")

sel = "+".join(vp1)
pml = os.path.join(SAIDA, "_capsideo_v59.pml")
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
select vp1all,  capsid and polymer and chain {sel}
select site22,  vp1all and resi 22
python
from pymol import cmd
print("  atomos capsid =", cmd.count_atoms("capsid"))
print("  atomos vp1all =", cmd.count_atoms("vp1all"))
print("  atomos site22 =", cmd.count_atoms("site22"))
print("  copias do sitio =", len(cmd.get_model("site22 and name CA").atom))
python end
hide everything
show spheres, capsid
color grey80, capsid
color firebrick, site22
set sphere_scale, 1.4, site22
orient capsid
zoom capsid, 4
png {SAIDA}/Fig6_capsid_outer_surface.png, width=3000, height=3000, dpi=600, ray=1
clip slab, 60
png {SAIDA}/Fig7_capsid_equatorial_section.png, width=3000, height=3000, dpi=600, ray=1
clip slab, 1000
""")
print(f"    pml: {pml}")
print("    renderizando")
try:
    r = subprocess.run([pm, "-cq", pml], capture_output=True, text=True, timeout=5400)
    for ln in (r.stdout or "").splitlines():
        if ln.strip():
            print("      " + ln)
    if r.returncode != 0:
        for ln in (r.stderr or "").splitlines()[-15:]:
            print("      ! " + ln)
except subprocess.TimeoutExpired:
    print("      estourou 90 min; rode o pml a mao")

print("\n" + "=" * 74)
print(f"figuras em {SAIDA}")
print("=" * 74)
