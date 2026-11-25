#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
V64 — acabamento final das figuras 2D.

    cd /home/juan && python3 V63_sasa_3vbs.py     # primeiro, cria a tabela
    cd /home/juan && python3 V64_acabamento.py    # depois, o acabamento

Fig0  rotulo do eixo y em duas linhas, pad maior
Fig1  rotulo "106" junto da tracejada
Fig2  idem
Fig3  idem
Fig4  da tabela do V63; residuo 22 em vermelho com rotulo "site 22"
Fig5  valor no topo de cada barra, residuo no eixo x (H22 / Q22), teto rotulado

Tudo em SVG vetorial e PNG 600 dpi, fundo branco solido.

Os valores das barras e os nomes dos residuos sao LIDOS dos seus arquivos,
nao digitados.  Se sairem diferentes de 863/857/840/846, o script mostra o
que leu — o numero do arquivo manda, nao o do briefing.

Nao altera nem sobrescreve original.  Escreve so em /home/juan/V58_figuras/.
"""

import os, sys, csv, glob, random
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE   = "/home/juan"
RAIZES = [f"{BASE}/V49_pegadas", f"{BASE}/abertos", f"{BASE}/paper_artifacts",
          f"{BASE}/figuras_V28", f"{BASE}/estruturas", BASE]
SAIDA  = f"{BASE}/V58_figuras"
TAB63  = f"{BASE}/V63_SASA_vs_RADIAL_3VBS.tsv"
PISO, SEMENTE = 106, 20261116
os.makedirs(SAIDA, exist_ok=True)

plt.rcParams.update({
    "font.family":"sans-serif","font.sans-serif":["DejaVu Sans","Arial"],
    "font.size":8,"axes.labelsize":8,"xtick.labelsize":7,"ytick.labelsize":7,
    "legend.fontsize":7,"axes.spines.top":False,"axes.spines.right":False,
    "axes.linewidth":.6,"xtick.major.width":.6,"ytick.major.width":.6,
    "legend.frameon":False,"figure.dpi":600,"svg.fonttype":"none",
    "figure.facecolor":"white","axes.facecolor":"white"})
CL, ESC, DEST = "#c8c8c8", "#6b6b6b", "#b2182b"
CINZA = "#9a9a9a"
LARG = 3.35

def salva(fig, nome, pad=0.06):
    for e in ("svg","png"):
        fig.savefig(os.path.join(SAIDA,f"{nome}.{e}"), format=e, dpi=600,
                    bbox_inches="tight", pad_inches=pad, facecolor="white")
    plt.close(fig); print(f"    -> {nome}.svg + .png (600 dpi)")

def marca_piso(ax, vertical=True, texto="106"):
    """rotulo discreto junto da tracejada"""
    if vertical:
        ax.axvline(PISO, color=DEST, lw=.9, ls="--")
        y = ax.get_ylim()[1]
        ax.annotate(texto, xy=(PISO, y), xytext=(3, -2),
                    textcoords="offset points", color=DEST,
                    fontsize=7, va="top", ha="left")
    else:
        ax.axhline(PISO, color=DEST, lw=.9, ls="--")
        x = ax.get_xlim()[1]
        ax.annotate(texto, xy=(x, PISO), xytext=(-3, 3),
                    textcoords="offset points", color=DEST,
                    fontsize=7, va="bottom", ha="right",
                    bbox=dict(boxstyle="square,pad=0.15", facecolor="white",
                              edgecolor="none", alpha=0.85))

def le(p):
    with open(p, newline="", encoding="utf-8", errors="replace") as f:
        d = "\t" if p.endswith(".tsv") else ","
        r = list(csv.DictReader(f, delimiter=d))
    return ([c.strip() for c in r[0].keys()] if r else []), r

def col(cs, alvos):
    for c in cs:
        if c.strip().lower() in alvos: return c
    return None

def num(v):
    if v is None: return None
    s = str(v).strip().replace(" ","")
    if not s or s.lower() in {"na","nan","none","-"}: return None
    pc = s.endswith("%");  s = s[:-1] if pc else s
    try: x = float(s.replace(",","."))
    except ValueError: return None
    return x/100.0 if pc else x

def verdade(v):
    return str(v).strip().lower() in {"1","true","sim","yes","t","s","verdadeiro"}

print("="*74); print("V64 — acabamento"); print("="*74)

# ======================================================= FIG 4
print("\n[FIG 4] da tabela do V63")
if not os.path.exists(TAB63):
    print(f"    {TAB63} NAO EXISTE.")
    print("    Rode antes:  python3 V63_sasa_3vbs.py")
    print("    Nao ha arquivo de SASA do capsideo em /home/juan — o unico")
    print("    caminho e recalcular.  Fig4 pulada.")
else:
    cs, rs = le(TAB63)
    R,K,N,RN = [],[],[],[]
    for d in rs:
        a, c = num(d.get("retencao")), num(d.get("shell_pct"))
        if a is None or c is None: continue
        R.append(a); K.append(c); N.append(num(d.get("residue")))
        RN.append(str(d.get("resname","")).strip())
    R,K = np.array(R), np.array(K)
    print(f"    residuos: {len(R)}")
    try:
        from scipy.stats import spearmanr
        print(f"    Spearman: {spearmanr(K,R).statistic:+.3f}")
    except ImportError: pass
    j = next((i for i,v in enumerate(N) if v == 22), None)
    fig, ax = plt.subplots(figsize=(LARG,2.7))
    ax.scatter(K, R, s=3.5, c=CINZA, lw=0, alpha=.75)
    if j is not None:
        ax.scatter([K[j]],[R[j]], s=30, facecolors="none",
                   edgecolors=DEST, linewidths=1.3, zorder=5)
        ax.annotate(f"{RN[j].capitalize()}22" if RN[j] else "site 22",
                    xy=(K[j],R[j]), xytext=(12,-10), textcoords="offset points",
                    color=DEST, fontsize=7,
                    arrowprops=dict(arrowstyle="-", lw=.6, color=DEST))
        print(f"    residuo 22 ({RN[j]}): percentil {K[j]:.2f}  retencao {R[j]:.4f}")
    else:
        print("    residuo 22 ausente da tabela; nao marcado")
    ax.set_xlabel("Radial Cα percentile within the shell")
    ax.set_ylabel("Retained accessible area (fraction)")
    ax.set_xlim(0,100); ax.set_ylim(0,1.05)
    salva(fig, "Fig4_sasa_retention_vs_radial", pad=0.10)
    print(f"    CONFERIR: {len(R)} pontos plotados")

# ======================================================= pegadas
print("\n[carga] pegadas")
peg = {}
for r in RAIZES:
    for p in sorted(glob.glob(os.path.join(r,"**","*_pegada*.tsv"), recursive=True)):
        b_ = os.path.basename(p); pdb = b_.split("_")[0].upper()
        if pdb in peg and "_V49" not in b_: continue
        cs, rs = le(p)
        cb = col(cs,{"blocking","bloqueio","blocking_atoms"})
        ck = col(cs,{"shell_pct","casca_pct","percentil_radial"})
        cf = col(cs,{"in_footprint","pegada","contato"})
        cn = col(cs,{"resseq","residuo","resi","num","residue"})
        if not (cb and cf): continue
        B,K2,F,NN,RR = [],[],[],[],[]
        for d in rs:
            vb = num(d.get(cb))
            if vb is None: continue
            B.append(vb); K2.append(num(d.get(ck)) if ck else float("nan"))
            F.append(verdade(d.get(cf))); NN.append(num(d.get(cn)) if cn else None)
            RR.append(str(d.get("resname","")).strip())
        if B: peg[pdb] = (np.array(B), np.array(K2,float), np.array(F), NN, RR)
print(f"    {len(peg)} complexos")

ordem = []
for r in RAIZES:
    g = glob.glob(os.path.join(r,"**","validacao_filtro.tsv"), recursive=True)
    if g:
        _, rs = le(sorted(g, key=len)[0])
        for d in rs:
            c = str(d.get("pdb","")).strip().upper()
            if c in peg and c not in ordem: ordem.append(c)
        break
base = ordem if ordem else sorted(peg)
e = base[:]; random.Random(SEMENTE).shuffle(e)
CAL, RET = set(e[:len(e)//2]), set(e[len(e)//2:])
print(f"    divisao {'validacao_filtro.tsv' if ordem else 'alfabetica'}: "
      f"cal {len(CAL)}  ret {len(RET)}")

def junta(cods):
    return (np.concatenate([peg[c][0] for c in cods]),
            np.concatenate([peg[c][1] for c in cods]),
            np.concatenate([peg[c][2] for c in cods]))

# ======================================================= FIG 1
print("\n[FIG 1]")
b,k,f = junta(sorted(CAL))
print(f"    contato {int(f.sum())}  nao-contato {int((~f).sum())}")
topo = float(np.percentile(b,99.5)); bins = np.linspace(0,topo,60)
fig, ax = plt.subplots(figsize=(LARG,2.4))
ax.hist(b[~f], bins=bins, color=CINZA, alpha=.85, density=True, lw=0,
        label="Non-contact residues")
ax.hist(b[f],  bins=bins, histtype="step", color="#2b2b2b", lw=1.0,
        density=True, label="Observed contact residues")
ax.set_xlabel("Outward-cone occupancy (atoms)")
ax.set_ylabel("Density"); ax.set_xlim(0,topo)
ax.legend(loc="upper right")
marca_piso(ax, vertical=True)
salva(fig, "Fig1_occupancy_distribution_calibration", pad=0.10)

# ======================================================= FIG 2
print("\n[FIG 2]")
b,k,f = junta(sorted(RET))
ts = np.arange(0, int(np.percentile(b,99.5))+1)
sens = np.array([(b[f] <= t).mean() for t in ts])*100
espe = np.array([(b[~f] > t).mean() for t in ts])*100
i = int(np.argmin(abs(ts-PISO)))
print(f"    contato {int(f.sum())}  nao-contato {int((~f).sum())}")
print(f"    no {PISO}: sens {sens[i]:.1f}%  espec {espe[i]:.1f}%")
fig, ax = plt.subplots(figsize=(LARG,2.4))
ax.plot(ts, sens, color="#2b2b2b", lw=1.0, label="Contact retention")
ax.plot(ts, espe, color=CINZA,     lw=1.0, label="Non-contact rejection")
ax.set_xlabel("Occupancy ceiling (atoms)")
ax.set_ylabel("Percentage of residues")
ax.set_xlim(ts.min(), ts.max()); ax.set_ylim(0,100)
ax.legend(loc="center right")
marca_piso(ax, vertical=True)
salva(fig, "Fig2_threshold_sweep_heldout", pad=0.10)

# ======================================================= FIG 3
print("\n[FIG 3]")
b,k,f = junta(sorted(peg))
m = np.isfinite(k) & np.isfinite(b); b,k,f = b[m],k[m],f[m]
print(f"    residuos: {len(b)}")
fig, ax = plt.subplots(figsize=(LARG,2.7))
ax.scatter(k[~f], b[~f], s=1.2, c=CL,    lw=0, alpha=.55, label="Non-contact residues")
ax.scatter(k[f],  b[f],  s=1.6, c="#2b2b2b", lw=0, alpha=.85, label="Observed contact residues")
ax.set_xlabel("Radial Cα percentile within the shell")
ax.set_ylabel("Outward-cone occupancy (atoms)")
ax.set_xlim(0,100)
ax.legend(loc="upper right", markerscale=4)
marca_piso(ax, vertical=False)
salva(fig, "Fig3_occupancy_vs_radial_percentile", pad=0.10)

# ======================================================= FIG 5
print("\n[FIG 5]")
EV = ("3VBF","3VBH","3VBO","3VBR","3VBS","3VBU","6I2K")
alvo = {}
perfis = set()
for r in RAIZES:
    perfis |= {os.path.realpath(p)
               for p in glob.glob(os.path.join(r,"**","*perfil_radial*"), recursive=True)
               if p.endswith((".tsv",".csv"))}
for p in sorted(perfis):
    cs, rs = le(p)
    cb = col(cs,{"blocking_atoms","blocking","bloqueio"})
    cn = col(cs,{"residue","resseq","residuo","resi"})
    if not (cb and cn): continue
    for d in rs:
        if str(d.get(cn,"")).strip() != "22": continue
        vb = num(d.get(cb))
        if vb is None: break
        nm = os.path.basename(p)
        rot = next((t for t in EV if t.lower() in nm.lower()), nm.split("_")[0].upper())
        alvo[rot] = (vb, str(d.get("resname","")).strip().upper())
        break
for pdb in [x for x in peg if x in EV]:
    if pdb in alvo: continue
    B,K2,F,NN,RR = peg[pdb]
    j2 = next((i for i,v in enumerate(NN) if v == 22), None)
    if j2 is not None: alvo[pdb] = (B[j2], RR[j2].upper())

if not alvo:
    print("    nada encontrado; Fig5 pulada")
else:
    UMA = {"HIS":"H","GLN":"Q","ASN":"N","LYS":"K","ARG":"R","GLU":"E",
           "ASP":"D","SER":"S","THR":"T","TYR":"Y","ALA":"A","GLY":"G"}
    nomes = sorted(alvo)
    vals  = [alvo[n][0] for n in nomes]
    rots  = []
    for n in nomes:
        rn = alvo[n][1]
        letra = UMA.get(rn, "")
        rots.append(f"{n}\n{letra}22" if letra else n)
        print(f"    {n:<6} {rn or '?':<4} ocupacao {alvo[n][0]:6.0f}")
    if len(nomes) < 4:
        print("    ATENCAO: o texto reporta quatro depositos")
    fig, ax = plt.subplots(figsize=(max(LARG,.80*len(nomes)+1.3), 2.5))
    barras = ax.bar(rots, vals, color=CINZA, width=.6, linewidth=0)
    for bb, v in zip(barras, vals):
        ax.annotate(f"{v:.0f}", xy=(bb.get_x()+bb.get_width()/2, v),
                    xytext=(0,2), textcoords="offset points",
                    ha="center", va="bottom", fontsize=7, color="#2b2b2b")
    ax.set_ylabel("Outward-cone occupancy (atoms)")
    ax.set_xlabel("Deposited structure")
    ax.set_ylim(0, max(vals)*1.18)
    marca_piso(ax, vertical=False, texto="ceiling = 106")
    salva(fig, "Fig5_site22_occupancy_across_deposits", pad=0.10)

# ======================================================= FIG 0
print("\n[FIG 0] so o rotulo; a geometria vem do V62")
f0 = os.path.join(SAIDA, "Fig0_cone_definition.svg")
print(f"    o V62 ja grava a figura; aqui so confirmo que existe:")
print(f"    {'OK' if os.path.exists(f0) else 'AUSENTE — rode V62'}  {f0}")
print("    o rotulo em duas linhas entra na proxima execucao do V62")

print("\n" + "="*74)
for p in sorted(glob.glob(os.path.join(SAIDA,"*.png"))):
    print(f"  {os.path.getsize(p)/1024:8.0f} kB  {os.path.basename(p)}")
print(f"\nfiguras em {SAIDA}")
print("="*74)
