#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V53 - os dois valores que faltam: resolucao do 3VBS e versao do LightDock.

O comando que eu mandei antes quebrou no zsh por erro meu: colei um fragmento
de marcacao no fim da linha. Isto e a mesma coisa como arquivo, sem citacao
aninhada e sem nada para o shell interpretar.

Somente leitura.

  python3 V53_dois_valores.py
"""
import glob, os, sys, subprocess

print("=" * 70)
print("resolucao e contagem de atomos do 3VBS")
print("=" * 70)
try:
    import gemmi
except ImportError:
    sys.exit("falta gemmi")

cands = []
for raiz in ("/home/juan/figuras_V28", "/home/juan/abertos",
             "/home/juan/cif_proveniencia", "/home/juan/cif_criterio", "/home/juan"):
    cands += [p for p in glob.glob(os.path.join(raiz, "**", "3VBS*.cif"), recursive=True)
              if "assembly" not in os.path.basename(p).lower()]
cands = sorted(set(os.path.realpath(p) for p in cands), key=len)
if not cands:
    print("nenhum 3VBS*.cif localizado")
else:
    p = cands[0]
    print(f"arquivo: {p}")
    if len(cands) > 1:
        print(f"  (ha {len(cands)} copias; usando a de caminho mais curto)")
    bl = gemmi.cif.read(p).sole_block()
    chaves = ("_refine.ls_d_res_high", "_em_3d_reconstruction.resolution",
              "_reflns.d_resolution_high", "_exptl.method",
              "_em_3d_reconstruction.method", "_struct.title")
    for k in chaves:
        v = bl.find_value(k)
        if v and v.strip() not in (".", "?"):
            print(f"  {k:<38}{v.strip().strip(chr(39))}")
    st = gemmi.read_structure(p)
    st.setup_entities()
    print(f"  cadeias no deposito                   {len(st[0])}")
    print(f"  atomos no deposito                    "
          f"{sum(1 for c in st[0] for r in c for a in r)}")
    nome = max(st.assemblies,
               key=lambda a: sum(len(g.operators) for g in a.generators)).name
    st.transform_to_assembly(nome, gemmi.HowToNameCopiedChain.AddNumber)
    print(f"  expandido, assembleia '{nome}'              {len(st[0])} cadeias, "
          f"{sum(1 for c in st[0] for r in c for a in r)} atomos")

print()
print("=" * 70)
print("versao do LightDock")
print("=" * 70)
achou = False
try:
    import lightdock
    v = getattr(lightdock, "__version__", None)
    print(f"  import lightdock -> {v if v else 'modulo sem __version__'}")
    print(f"  caminho: {os.path.dirname(lightdock.__file__)}")
    achou = bool(v)
except ImportError:
    print("  lightdock nao importavel neste ambiente")
for cmd in (["pip", "show", "lightdock"], ["lightdock", "--version"],
            ["lightdock3", "--version"]):
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        saida = (r.stdout or r.stderr).strip()
        if saida:
            print(f"  $ {' '.join(cmd)}")
            for ln in saida.splitlines()[:4]:
                print(f"      {ln}")
            achou = True
    except Exception:
        continue
if not achou:
    print("\n  Nao localizada. Se o docking rodou em outro ambiente ou em outra")
    print("  maquina, a versao tem de vir do registro daquela execucao; nao")
    print("  inferir da versao instalada hoje.")
