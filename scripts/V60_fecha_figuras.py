#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
V60 — fecha as tres pendencias das figuras.

    cd /home/juan && python3 V60_fecha_figuras.py

(A) Fig4 saiu VAZIA.  O arquivo foi achado mas nenhuma linha passou na leitura
    numerica.  Aqui eu imprimo o cabecalho e as primeiras linhas antes de
    qualquer coisa, aceito virgula decimal e porcentagem, e, se nao houver
    coluna de retencao utilizavel, calculo retencao = montada / isolada.

(B) Fig5 perdeu o 6I2K.  O 22 do 6I2K esta nas pegadas e o dos 3VB* esta nos
    perfis radiais.  Agora junta as duas fontes e diz de onde veio cada barra.

(C) Fig1 e Fig2 nao batem com o texto (96,5% / 65,3%).  Monto as DUAS
    convencoes de divisao a partir da semente 20261116 — a ordem do
    validacao_filtro.tsv e a ordem alfabetica — calculo sens/espec no 106 em
    cada uma, digo qual reproduz o texto e redesenho com essa.

Nao altera nem sobrescreve original.  Nao recalcula descritor.
Escreve so em /home/juan/V58_figuras/.
"""

import os, sys, csv, glob, random
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
try:
    from scipy.stats import spearmanr
except ImportError:
    spearmanr = None

BASE   = "/home/juan"
RAIZES = [f"{BASE}/V49_pegadas", f"{BASE}/abertos", f"{BASE}/paper_artifacts",
          f"{BASE}/figuras_V28", f"{BASE}/estruturas", BASE]
SAIDA  = f"{BASE}/V58_figuras"
PISO, SEMENTE = 106, 20261116
os.makedirs(SAIDA, exist_ok=True)

plt.rcParams.update({
    "font.family":"sans-serif","font.sans-serif":["DejaVu Sans","Arial"],
    "font.size":8,"axes.labelsize":8,"xtick.labelsize":7,"ytick.labelsize":7,
    "legend.fontsize":7,"axes.spines.top":False,"axes.spines.right":False,
    "axes.linewidth":.6,"xtick.major.width":.6,"ytick.major.width":.6,
    "legend.frameon":False,"figure.dpi":600,"svg.fonttype":"none"})
CINZA, ESCURO, DESTAQUE = "#9a9a9a", "#2b2b2b", "#b2182b"
LARG = 3.35

def salva(fig, nome):
    for e in ("svg","png"):
        fig.savefig(os.path.join(SAIDA,f"{nome}.{e}"), format=e,
                    bbox_inches="tight", pad_inches=0.04)
    plt.close(fig); print(f"    -> {nome}.svg + .png")

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
    """aceita virgula decimal, porcentagem, espacos"""
    if v is None: return None
    s = str(v).strip().replace(" ", "")
    if not s or s.lower() in {"na","nan","none","-",""}: return None
    pc = s.endswith("%")
    if pc: s = s[:-1]
    s = s.replace(",", ".")
    try:
        x = float(s)
    except ValueError:
        return None
    return x/100.0 if pc else x

def verdade(v):
    return str(v).strip().lower() in {"1","true","sim","yes","t","s","verdadeiro"}

print("=" * 74); print("V60 — fechamento das figuras"); print("=" * 74)

# ==================================================== (A) Fig4
print("\n[A] Fig4 — retencao de SASA")
alvo = None
for r in RAIZES:
    g = [x for x in glob.glob(os.path.join(r,"**","V28_SASA_vs_RADIAL*"), recursive=True)
         if x.endswith((".tsv",".csv"))]
    if g: alvo = sorted(g, key=len)[0]; break
if not alvo:
    for r in RAIZES:
        for p in sorted(glob.glob(os.path.join(r,"**","*.tsv"), recursive=True)
                      + glob.glob(os.path.join(r,"**","*.csv"), recursive=True)):
            if any(t in os.path.basename(p).lower() for t in ("sasa","retenc","acess")):
                alvo = p; break
        if alvo: break

if not alvo:
    print("    arquivo nao localizado")
else:
    print(f"    arquivo: {alvo}")
    cs, rs = le(alvo)
    print(f"    colunas ({len(cs)}): {cs}")
    print(f"    linhas: {len(rs)}")
    print("    tres primeiras linhas:")
    for d in rs[:3]:
        print("      " + " | ".join(f"{k}={d[k]}" for k in cs))

    cr  = col(cs, {"retencao","retained","fracao","fraction","ret","retention"})
    ck  = col(cs, {"shell_pct","casca","casca_pct","percentil","percentil_radial",
                   "radial_pct","percentile"})
    ci  = col(cs, {"sasa_isolada","sasa_isolated","isolada","iso_a2","sasa_iso","iso"})
    cm  = col(cs, {"sasa_montada","sasa_assembled","montada","assembled",
                   "sasa_240","sasa_capsid","mont"})
    cn  = col(cs, {"residue","resseq","residuo","resi","num"})
    print(f"    mapeado: retencao={cr}  percentil={ck}  isolada={ci}  "
          f"montada={cm}  residuo={cn}")

    R,K,N = [],[],[]
    sem_ret = sem_k = 0
    for d in rs:
        c = num(d.get(ck)) if ck else None
        if c is None: sem_k += 1; continue
        a = num(d.get(cr)) if cr else None
        if a is None and ci and cm:                 # calcula a retencao
            vi, vm = num(d.get(ci)), num(d.get(cm))
            if vi and vi > 0 and vm is not None:
                a = vm / vi
        if a is None: sem_ret += 1; continue
        if a > 1.5: a = a/100.0                     # veio em porcentagem
        R.append(a); K.append(c); N.append(num(d.get(cn)) if cn else None)

    print(f"    aproveitadas {len(R)}   sem percentil {sem_k}   sem retencao {sem_ret}")
    if not R:
        print("    AINDA VAZIA. me manda as tres linhas impressas acima que eu")
        print("    acerto o mapeamento de coluna.")
    else:
        R, K = np.array(R), np.array(K)
        if spearmanr and len(R) > 2:
            print(f"    Spearman retencao x percentil: {spearmanr(K,R).statistic:+.3f}  n={len(R)}")
        print(f"    retencao: min {R.min():.3f}  mediana {np.median(R):.3f}  max {R.max():.3f}")
        fig, ax = plt.subplots(figsize=(LARG,2.6))
        ax.scatter(K, R, s=3, c=CINZA, lw=0, alpha=.7)
        j = next((i for i,v in enumerate(N) if v == 22), None)
        if j is not None:
            ax.scatter([K[j]],[R[j]], s=26, facecolors="none",
                       edgecolors=DESTAQUE, linewidths=1.1)
            print(f"    residuo 22: percentil {K[j]:.2f}  retencao {R[j]:.4f}")
        else:
            print("    residuo 22 nao marcado (sem coluna de numeracao ou ausente)")
        ax.set_xlabel("Radial Cα percentile within the shell")
        ax.set_ylabel("Retained accessible area (fraction)")
        ax.set_xlim(0,100); ax.set_ylim(0,1.05)
        salva(fig, "Fig4_sasa_retention_vs_radial")

# ==================================================== carrega pegadas
print("\n[carga] pegadas")
peg, orig = {}, {}
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
        B,K,F,NN = [],[],[],[]
        for d in rs:
            vb = num(d.get(cb))
            if vb is None: continue
            B.append(vb); K.append(num(d.get(ck)) if ck else float("nan"))
            F.append(verdade(d.get(cf))); NN.append(num(d.get(cn)) if cn else None)
        if B:
            peg[pdb] = (np.array(B), np.array(K,float), np.array(F), NN)
            orig[pdb] = p
print(f"    {len(peg)} complexos")

# ==================================================== (B) Fig5
print("\n[B] Fig5 — sitio 22, juntando perfis radiais e pegadas")
EV = ("3VBF","3VBH","3VBO","3VBR","3VBS","3VBU","6I2K")
alvo22 = {}
perfis = set()
for r in RAIZES:
    perfis |= {os.path.realpath(p)
               for p in glob.glob(os.path.join(r,"**","*perfil_radial*"), recursive=True)
               if p.endswith((".tsv",".csv"))}
for p in sorted(perfis):
    cs, rs = le(p)
    cb = col(cs,{"blocking_atoms","blocking","bloqueio"})
    cn = col(cs,{"residue","resseq","residuo","resi"})
    ck = col(cs,{"shell_pct","casca_pct"})
    if not (cb and cn): continue
    for d in rs:
        if str(d.get(cn,"")).strip() != "22": continue
        vb = num(d.get(cb))
        if vb is None: break
        nome = os.path.basename(p)
        rot = next((t for t in EV if t.lower() in nome.lower()), nome.split("_")[0].upper())
        alvo22[rot] = (vb, num(d.get(ck)) if ck else None,
                       str(d.get("resname","")).strip(), "perfil radial")
        break
for pdb in [k for k in peg if k in EV]:
    if pdb in alvo22: continue
    B,K,F,NN = peg[pdb]
    j = next((i for i,v in enumerate(NN) if v == 22), None)
    if j is not None:
        alvo22[pdb] = (B[j], K[j] if np.isfinite(K[j]) else None, "", "pegada")

if not alvo22:
    print("    nada encontrado")
else:
    for k_ in sorted(alvo22):
        vb, vk, rn, fonte = alvo22[k_]
        print(f"    {k_:<6} {rn:<4} ocupacao {vb:6.0f}  "
              f"percentil {('%.2f'%vk) if vk is not None else ' n/d'}  "
              f"{'passa' if vb <= PISO else 'REPROVA'}   [{fonte}]")
    nomes = sorted(alvo22); vals = [alvo22[n][0] for n in nomes]
    print(f"    {len(nomes)} depositos: {', '.join(nomes)}")
    if len(nomes) < 4:
        print("    ATENCAO: o texto reporta quatro depositos (3VBF, 3VBH, 3VBS, 6I2K).")
        print("    Faltando algum, a figura nao cobre o que o texto afirma.")
    fig, ax = plt.subplots(figsize=(max(LARG,.70*len(nomes)+1.3), 2.35))
    ax.bar(nomes, vals, color=CINZA, width=.6, linewidth=0)
    ax.axhline(PISO, color=DESTAQUE, lw=.9, ls="--")
    ax.set_ylabel("Outward-cone occupancy (atoms)")
    ax.set_xlabel("Deposited structure")
    ax.set_ylim(0, max(vals)*1.14)
    salva(fig, "Fig5_site22_occupancy_across_deposits")

# ==================================================== (C) divisao
print("\n[C] Fig1 e Fig2 — qual convencao de divisao reproduz o texto")
ordem_tsv = []
for r in RAIZES:
    g = glob.glob(os.path.join(r,"**","validacao_filtro.tsv"), recursive=True)
    if g:
        p = sorted(g, key=len)[0]
        print(f"    validacao_filtro.tsv: {p}")
        _, rs = le(p)
        for d in rs:
            c = str(d.get("pdb","")).strip().upper()
            if c and c in peg and c not in ordem_tsv:
                ordem_tsv.append(c)
        break
if not ordem_tsv:
    print("    validacao_filtro.tsv nao localizado; so consigo a convencao alfabetica")

def avalia(cal, ret):
    b = np.concatenate([peg[c][0] for c in ret])
    f = np.concatenate([peg[c][2] for c in ret])
    return (b[f] <= PISO).mean()*100, (b[~f] > PISO).mean()*100

convs = {}
if ordem_tsv:
    e = ordem_tsv[:]; random.Random(SEMENTE).shuffle(e)
    convs["A (ordem do validacao_filtro.tsv)"] = (set(e[:len(e)//2]), set(e[len(e)//2:]))
alf = sorted(peg); e = alf[:]; random.Random(SEMENTE).shuffle(e)
convs["B (ordem alfabetica)"] = (set(e[:len(e)//2]), set(e[len(e)//2:]))

melhor, menor = None, 1e9
print(f"\n    {'convencao':<34}{'sens':>8}{'espec':>8}   distancia de 96,5/65,3")
for rot,(cal,ret) in convs.items():
    s, p_ = avalia(cal, ret)
    d = abs(s-96.5) + abs(p_-65.3)
    print(f"    {rot:<34}{s:7.1f}%{p_:7.1f}%   {d:6.2f}")
    if d < menor: menor, melhor = d, (rot, cal, ret)

rot, CAL, RET = melhor
print(f"\n    adotada: {rot}")
print(f"      calibracao: {', '.join(sorted(CAL))}")
print(f"      retido:     {', '.join(sorted(RET))}")
if menor > 2.0:
    print("    ATENCAO: nenhuma das duas reproduz 96,5/65,3 de perto.")
    print("    Entao a diferenca nao e so convencao de divisao — pode ser o")
    print("    conjunto de complexos que entrou aqui.  Me manda esta saida.")

for rotulo, cods, nome, tipo in (
        ("calibracao", sorted(CAL), "Fig1_occupancy_distribution_calibration", "hist"),
        ("retido",     sorted(RET), "Fig2_threshold_sweep_heldout",            "sweep")):
    b = np.concatenate([peg[c][0] for c in cods])
    f = np.concatenate([peg[c][2] for c in cods])
    print(f"\n    {nome}: {rotulo}, contato {int(f.sum())}, nao-contato {int((~f).sum())}")
    if tipo == "hist":
        topo = float(np.percentile(b, 99.5)); bins = np.linspace(0, topo, 60)
        fig, ax = plt.subplots(figsize=(LARG,2.3))
        ax.hist(b[~f], bins=bins, color=CINZA, alpha=.85, density=True, lw=0,
                label="Non-contact residues")
        ax.hist(b[f],  bins=bins, histtype="step", color=ESCURO, lw=1.0,
                density=True, label="Observed contact residues")
        ax.axvline(PISO, color=DESTAQUE, lw=.9, ls="--")
        ax.set_xlabel("Outward-cone occupancy (atoms)")
        ax.set_ylabel("Density"); ax.set_xlim(0, topo)
        ax.legend(loc="upper right")
    else:
        ts = np.arange(0, int(np.percentile(b,99.5))+1)
        sens = np.array([(b[f] <= t).mean() for t in ts])*100
        espe = np.array([(b[~f] > t).mean() for t in ts])*100
        i = int(np.argmin(abs(ts-PISO)))
        print(f"      no {PISO}: sens {sens[i]:.1f}%  espec {espe[i]:.1f}%")
        fig, ax = plt.subplots(figsize=(LARG,2.3))
        ax.plot(ts, sens, color=ESCURO, lw=1.0, label="Contact retention")
        ax.plot(ts, espe, color=CINZA,  lw=1.0, label="Non-contact rejection")
        ax.axvline(PISO, color=DESTAQUE, lw=.9, ls="--")
        ax.set_xlabel("Occupancy ceiling (atoms)")
        ax.set_ylabel("Percentage of residues")
        ax.set_xlim(ts.min(), ts.max()); ax.set_ylim(0,100)
        ax.legend(loc="center right")
    salva(fig, nome)

print("\n" + "=" * 74)
print(f"figuras em {SAIDA}")
print("=" * 74)
