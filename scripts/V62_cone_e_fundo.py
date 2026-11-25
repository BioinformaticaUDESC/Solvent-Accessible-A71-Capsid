#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
V62 — a figura do cone, e a correcao de fundo das renderizacoes.

    cd /home/juan && python3 V62_cone_e_fundo.py

(1) FIGURA DO CONE.  O manuscrito define o descritor em palavras e nunca o
    mostra.  Esta figura NAO e um desenho esquematico: e um corte real da sua
    montagem 3VBS.  Toma o Ca do residuo 22, o centroide da particula, o eixo
    radial externo entre os dois, e desenha o cone de semiangulo 45 graus e
    30 A sobre os atomos de verdade.  Os atomos que o descritor conta saem
    marcados; os de fora, em cinza claro.  O numero impresso e a contagem
    real, e tem de bater com a ocupacao da tabela.

    Isso mostra uma coisa que so o texto dizia: o cone conta atomos da
    propria cadeia lateral e dos vizinhos sequenciais, e nao e teste de linha
    de visada.

(2) FUNDO DAS RENDERIZACOES.  Os .pml usavam ray_opaque_background, 0, que
    deixa o fundo TRANSPARENTE — aparece preto em visualizador escuro e pode
    virar bloco preto no PDF da revista.  Aqui sai 1, fundo branco solido.

Nao altera nem sobrescreve original.  Escreve so em /home/juan/V58_figuras/.
"""

import os, sys, glob, math, shutil, subprocess
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE   = "/home/juan"
RAIZES = [f"{BASE}/figuras_V28", f"{BASE}/abertos", f"{BASE}/paper_artifacts",
          f"{BASE}/estruturas", BASE]
SAIDA  = f"{BASE}/V58_figuras"
RAIO, SEMIANG, RES = 30.0, 45.0, 22
os.makedirs(SAIDA, exist_ok=True)

plt.rcParams.update({
    "font.family":"sans-serif","font.sans-serif":["DejaVu Sans","Arial"],
    "font.size":8,"axes.labelsize":8,"xtick.labelsize":7,"ytick.labelsize":7,
    "legend.fontsize":7,"axes.spines.top":False,"axes.spines.right":False,
    "axes.linewidth":.6,"xtick.major.width":.6,"ytick.major.width":.6,
    "legend.frameon":False,"figure.dpi":600,"svg.fonttype":"none"})
CINZA, ESCURO, DESTAQUE = "#c8c8c8", "#6b6b6b", "#b2182b"

try:
    import gemmi
except ImportError:
    sys.exit("gemmi ausente: pip install --user gemmi")

print("=" * 74); print("V62 — figura do cone e correcao de fundo"); print("=" * 74)

# ---------------------------------------------------------- estrutura
cand = []
for r in RAIZES:
    for ext in ("cif","pdb","cif.gz","pdb.gz","ent"):
        cand += glob.glob(os.path.join(r,"**",f"*3[Vv][Bb][Ss]*.{ext}"), recursive=True)
def pontua(p):
    b = os.path.basename(p).lower()
    return (sum(v for t,v in (("assembly",40),("montagem",40),("240",20),
                              ("expand",15)) if t in b), os.path.getsize(p))
cand = sorted(set(cand), key=pontua, reverse=True)
if not cand:
    sys.exit("nenhum arquivo 3VBS encontrado")
estrut = cand[0]
print(f"\nestrutura: {estrut}")

st = gemmi.read_structure(estrut)
st.setup_entities(); st.remove_ligands_and_waters(); st.remove_hydrogens()
mod = st[0]

inv = [(ch.name, sum(1 for _ in ch), any(r.seqid.num == RES for r in ch)) for ch in mod]
vp1 = [c for c in inv if 270 <= c[1] <= 310 and c[2]]
if not vp1:
    print("  inventario:")
    for n,L,t in inv[:10]: print(f"    {n:<5} {L:4d}  tem o {RES}: {t}")
    sys.exit("nenhuma cadeia VP1 identificada")
REF = vp1[0][0]
print(f"cadeia de referencia: {REF}   VP1 no modelo: {len(vp1)}")

# ---------------------------------------------------------- coordenadas
pos, cas = [], []
for ch in mod:
    for r in ch:
        for a in r:
            pos.append([a.pos.x, a.pos.y, a.pos.z])
        ca = r.find_atom("CA", "*")
        if ca: cas.append([ca.pos.x, ca.pos.y, ca.pos.z])
pos = np.array(pos); cas = np.array(cas)
centro = cas.mean(axis=0)
print(f"atomos: {len(pos)}   Ca: {len(cas)}")
print(f"centroide: ({centro[0]:.1f}, {centro[1]:.1f}, {centro[2]:.1f})")

cad = next(ch for ch in mod if ch.name == REF)
r22 = next((r for r in cad if r.seqid.num == RES), None)
if r22 is None: sys.exit(f"residuo {RES} ausente na cadeia {REF}")
ca22 = r22.find_atom("CA", "*")
if ca22 is None: sys.exit(f"residuo {RES} sem CA")
p = np.array([ca22.pos.x, ca22.pos.y, ca22.pos.z])
print(f"residuo {RES}: {r22.name}   raio do Ca {np.linalg.norm(p-centro):.1f} A")

# ---------------------------------------------------------- o descritor
u = (p - centro); u /= np.linalg.norm(u)          # eixo radial externo
d = pos - p
dist = np.linalg.norm(d, axis=1)
viz = dist > 1e-6                                  # exclui o proprio vertice
cosang = np.zeros(len(d))
cosang[viz] = (d[viz] @ u) / dist[viz]
dentro = viz & (dist <= RAIO) & (cosang >= math.cos(math.radians(SEMIANG)))
n_dentro = int(dentro.sum())
print(f"\nocupacao do cone (45 graus, {RAIO:.0f} A): {n_dentro} atomos")
print("  CONFIRA: tem de bater com a ocupacao do residuo 22 na sua tabela")
print("  (o texto publica 840 para o 3VBS)")

# ---------------------------------------------------------- corte 2D
e1 = u
tmp = np.array([0.,0.,1.])
if abs(tmp @ e1) > 0.9: tmp = np.array([1.,0.,0.])
e2 = tmp - (tmp @ e1)*e1; e2 /= np.linalg.norm(e2)
e3 = np.cross(e1, e2)

loc = pos - p
X, Y, Z = loc @ e2, loc @ e1, loc @ e3
MARGEM, ESPESSURA = 38.0, 4.0
fatia = (np.abs(Z) <= ESPESSURA) & (np.abs(X) <= MARGEM) & (Y >= -14) & (Y <= MARGEM)
print(f"atomos na fatia de +-{ESPESSURA:.0f} A: {int(fatia.sum())}")

f_in  = fatia & dentro
f_out = fatia & ~dentro

fig, ax = plt.subplots(figsize=(3.35, 3.0))
ax.scatter(X[f_out], Y[f_out], s=5, c=CINZA,  lw=0, alpha=.85,
           label="Particle atoms outside the cone")
ax.scatter(X[f_in],  Y[f_in],  s=5, c=ESCURO, lw=0, alpha=.95,
           label="Atoms counted by the descriptor")

t = math.radians(SEMIANG)
for s in (-1, 1):
    ax.plot([0, s*RAIO*math.sin(t)], [0, RAIO*math.cos(t)],
            color=DESTAQUE, lw=.9)
arco = np.linspace(-t, t, 200)
ax.plot(RAIO*np.sin(arco), RAIO*np.cos(arco), color=DESTAQUE, lw=.9)
ax.annotate("", xy=(0, MARGEM*0.97), xytext=(0, 0),
            arrowprops=dict(arrowstyle="->", lw=.8, color=DESTAQUE))
ax.scatter([0], [0], s=30, facecolors="none", edgecolors=DESTAQUE, linewidths=1.2)

ax.set_xlabel("Tangential distance from the Cα (Å)")
ax.set_ylabel("Distance along outward\nradial axis (Å)")
ax.set_xlim(-MARGEM, MARGEM); ax.set_ylim(-14, MARGEM)
ax.set_aspect("equal")
ax.legend(loc="lower center", ncol=1)
for ext in ("svg","png"):
    fig.savefig(os.path.join(SAIDA, f"Fig0_cone_definition.{ext}"),
                format=ext, dpi=600, bbox_inches="tight", pad_inches=0.10,
                facecolor="white")
plt.close(fig)
print(f"  -> {SAIDA}/Fig0_cone_definition.svg + .png")

# ---------------------------------------------------------- fundo opaco
print("\n[fundo] re-renderizando Fig6, Fig7 e Fig8 com fundo branco solido")
pm = shutil.which("pymol")
if not pm:
    print("  PyMOL ausente; pulado")
else:
    sel = "+".join(c[0] for c in vp1)
    pml = os.path.join(SAIDA, "_capsideo_v62.pml")
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
set ray_opaque_background, 1
bg_color white
load {estrut}, capsid
select vp1all, capsid and polymer and chain {sel}
select site22, vp1all and resi {RES}
select vp1ref, capsid and polymer and chain {REF}
select ref22,  vp1ref and resi {RES}
python
from pymol import cmd
print("  capsid =", cmd.count_atoms("capsid"))
print("  site22 =", cmd.count_atoms("site22"),
      " copias =", len(cmd.get_model("site22 and name CA").atom))
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
hide everything
show cartoon, vp1ref
color grey70, vp1ref
show spheres, ref22
color firebrick, ref22
orient vp1ref
zoom vp1ref, 3
png {SAIDA}/Fig8_protomer_site22.png, width=3000, height=2400, dpi=600, ray=1
""")
    print(f"  pml: {pml}")
    try:
        r = subprocess.run([pm, "-cq", pml], capture_output=True, text=True, timeout=5400)
        for ln in (r.stdout or "").splitlines():
            if ln.strip(): print("    " + ln)
        if r.returncode != 0:
            for ln in (r.stderr or "").splitlines()[-12:]: print("    ! " + ln)
    except subprocess.TimeoutExpired:
        print("    estourou 90 min; rode o pml a mao")

print("\n" + "=" * 74)
print(f"figuras em {SAIDA}")
print("nenhum arquivo original foi alterado")
print("=" * 74)
