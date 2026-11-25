#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
V66 — monta o repositorio do projeto para subir ao GitHub.

    cd /home/juan && python3 V66_monta_repositorio.py

Monta /home/juan/repo_a71 com dados, scripts, figuras e manifesto de hashes.
COPIA, nunca move: nada em /home/juan e alterado ou removido.

NAO faz push.  Ele prepara tudo, faz o commit local e IMPRIME os comandos de
push para voce rodar depois de revisar.  O repositorio e PUBLICO; vale olhar a
lista de arquivos antes de publicar.

Repositorio: https://github.com/BioinformaticaUDESC/Solven-Acessiblre-A71-Capsid
"""

import os, sys, glob, shutil, hashlib, subprocess, datetime

BASE   = "/home/juan"
DESTINO = f"{BASE}/repo_a71"
REMOTO = "https://github.com/BioinformaticaUDESC/Solven-Acessiblre-A71-Capsid.git"
RAIZES = [f"{BASE}/V49_pegadas", f"{BASE}/abertos", f"{BASE}/paper_artifacts",
          f"{BASE}/figuras_V28", f"{BASE}/V58_figuras", f"{BASE}/V65_submissao",
          BASE]

if os.path.exists(DESTINO):
    sys.exit(f"{DESTINO} ja existe. Apague ou renomeie antes de rodar de novo.")

print("=" * 74)
print("V66 — montagem do repositorio")
print("=" * 74)

for sub in ("data/per_residue", "data/radial_profiles", "scripts", "figures",
            "manuscript"):
    os.makedirs(os.path.join(DESTINO, sub), exist_ok=True)

copiados = []

def copia(origem, sub, novo_nome=None):
    if not os.path.isfile(origem):
        return False
    alvo = os.path.join(DESTINO, sub, novo_nome or os.path.basename(origem))
    if os.path.exists(alvo):
        return False
    shutil.copy2(origem, alvo)
    copiados.append(os.path.relpath(alvo, DESTINO))
    return True

def busca(padrao, so_raiz=None):
    achados = []
    for r in ([so_raiz] if so_raiz else RAIZES):
        if not os.path.isdir(r):
            continue
        achados += glob.glob(os.path.join(r, "**", padrao), recursive=True)
    return sorted({os.path.realpath(a) for a in achados})

print("\n[1] tabelas por residuo")
for p in busca("*_pegada*.tsv"):
    copia(p, "data/per_residue")
for p in busca("*perfil_radial*.tsv") + busca("*perfil_radial*.csv"):
    copia(p, "data/radial_profiles")
for nome in ("validacao_filtro.tsv", "V63_SASA_vs_RADIAL_3VBS.tsv",
             "V61_SASA_vs_RADIAL_3VBS.tsv"):
    for p in busca(nome):
        copia(p, "data")
print(f"    {len([c for c in copiados if c.startswith('data')])} arquivos")

print("\n[2] scripts")
SCRIPTS = ["filtro_externo.py", "piso_um_eixo_26.py", "piso_calibra_testa.py",
           "perfil_radial_vp1.py",
           "V49_filtro_canonico.py", "V50_assembleia_oficial.py",
           "V51_completude_nos_39.py", "V52_rastreio_e_robustez.py",
           "V53_dois_valores.py", "V54_selecao_e_denominador.py",
           "V58_gera_figuras.py", "V59_conserta_figuras.py",
           "V60_fecha_figuras.py", "V62_cone_e_fundo.py",
           "V63_sasa_3vbs.py", "V64_acabamento.py", "V65_monta_paineis.py"]
for nome in SCRIPTS:
    for p in busca(nome):
        copia(p, "scripts")
for p in busca("V49_METODO.txt"):
    copia(p, "scripts")
print(f"    {len([c for c in copiados if c.startswith('scripts')])} arquivos")

print("\n[3] figuras")
for p in busca("Fig[0-8]_*.png") + busca("Fig[0-8]_*.svg"):
    copia(p, "figures")
for p in busca("Figure_[1-4]_*.png") + busca("Figure_[1-4]_*.tif"):
    copia(p, "figures")
print(f"    {len([c for c in copiados if c.startswith('figures')])} arquivos")

print("\n[4] manuscrito")
for p in (busca("Parkinsonismo_*Submissao*.docx") +
          busca("COVER_LETTER_*.docx") +
          busca("SUPPLEMENTARY_METHODS*.md")):
    copia(p, "manuscript")
n_ms = len([c for c in copiados if c.startswith('manuscript')])
print(f"    {n_ms} arquivos")
if n_ms == 0:
    print("    nenhum encontrado — copie o .docx e a cover letter a mao para")
    print(f"    {DESTINO}/manuscript/ antes do push")

# ------------------------------------------------------------ manifesto
print("\n[5] manifesto SHA-256")
linhas = []
for rel in sorted(copiados):
    caminho = os.path.join(DESTINO, rel)
    h = hashlib.sha256()
    with open(caminho, "rb") as f:
        for bloco in iter(lambda: f.read(1 << 20), b""):
            h.update(bloco)
    linhas.append(f"{h.hexdigest()}  {rel}")
with open(os.path.join(DESTINO, "MANIFEST.sha256"), "w") as f:
    f.write("\n".join(linhas) + "\n")
print(f"    {len(linhas)} arquivos com hash")

# ------------------------------------------------------------ README
hoje = datetime.date.today().isoformat()
n_peg = len(glob.glob(os.path.join(DESTINO, "data/per_residue/*")))
n_perf = len(glob.glob(os.path.join(DESTINO, "data/radial_profiles/*")))
n_fig = len(glob.glob(os.path.join(DESTINO, "figures/*")))
n_scr = len(glob.glob(os.path.join(DESTINO, "scripts/*")))

with open(os.path.join(DESTINO, "README.md"), "w") as f:
    f.write(f"""# Solvent-accessible yet externally occluded

Data and analysis code for the manuscript *Solvent-Accessible Yet Externally
Occluded: Directional Analysis of a Docking-Selected Site in the Mature
Enterovirus A71 Capsid*, submitted to *Proteins: Structure, Function, and
Bioinformatics*.

## What this is

A protein-protein docking run against an isolated chain can select a site that
the assembled particle makes unreachable. Solvent-accessible surface area does
not detect this, because it is not directional. This repository holds the data
and code for **outward-cone atomic occupancy**: the number of atoms of the
assembled particle inside a cone of 45 degrees half-angle, vertex at a
residue's CA, axis along the outward radial vector from the particle centroid,
searched to 30 A.

## Layout

| path | contents |
|---|---|
| `data/per_residue/` | per-residue occupancy and contact labels, {n_peg} files |
| `data/radial_profiles/` | radial profiles of the EV-A71 deposits, {n_perf} files |
| `data/` | validation table and the SASA-versus-radial table for 3VBS |
| `scripts/` | calibration, filter and figure code, {n_scr} files |
| `figures/` | individual panels and the four composite figures, {n_fig} files |
| `manuscript/` | submitted manuscript and cover letter |
| `MANIFEST.sha256` | SHA-256 of every file above |

## Reproducing the figures

The scripts read the tables in `data/` and write to a new directory; they do not
modify their inputs. They expect the assembled coordinate files, which are not
redistributed here: download them from the PDB (3VBF, 3VBH, 3VBS, 6I2K and the
reference set) and point the scripts at them.

```
python3 scripts/V63_sasa_3vbs.py      # recomputes the SASA table
python3 scripts/V62_cone_e_fundo.py   # cone figure and capsid renders
python3 scripts/V64_acabamento.py     # the 2D figures
python3 scripts/V65_monta_paineis.py  # the four composites
```

Requires Python 3 with numpy, scipy, matplotlib, gemmi, freesasa and Pillow;
PyMOL for the capsid renders.

## What is not here

The exploratory docking layer that first nominated the site is **not
reproducible** and is not included. Its run tree is no longer available and the
surviving summary tables are not sufficient to reconstruct those executions.
The manuscript states this. Nothing in the calibration or evaluation depends on
that layer.

## Structures

Coordinates come from the Protein Data Bank and are not redistributed here.
The case study uses 3VBF, 3VBH, 3VBS (Wang et al., *Nat Struct Mol Biol* 2012,
doi:10.1038/nsmb.2255) and 6I2K (Zhou et al., *Nat Microbiol* 2019,
doi:10.1038/s41564-018-0319-z).

## Citing

See `CITATION.cff`. Assembled {hoje}.
""")
print(f"\n[6] README.md escrito")

with open(os.path.join(DESTINO, "CITATION.cff"), "w") as f:
    f.write(f"""cff-version: 1.2.0
message: "If you use this code or these data, please cite the manuscript."
title: "Solvent-accessible yet externally occluded: data and analysis code"
authors:
  - family-names: Teixeira
    given-names: Juan Philippe
    orcid: "https://orcid.org/0000-0001-7323-2919"
date-released: "{hoje}"
repository-code: "https://github.com/BioinformaticaUDESC/Solven-Acessiblre-A71-Capsid"
license: MIT
""")

with open(os.path.join(DESTINO, ".gitignore"), "w") as f:
    f.write("__pycache__/\n*.pyc\n.ipynb_checkpoints/\n*.cif\n*.pdb\n*.ent\n"
            "*.cif.gz\n*.pdb.gz\n.DS_Store\n")

print("[7] CITATION.cff e .gitignore escritos")
print("\n    LICENSE nao foi escrito: escolha a licenca voce.")
print("    Usual: MIT para o codigo, CC-BY-4.0 para os dados.")

# ------------------------------------------------------------ git
print("\n[8] git")
def git(*a):
    r = subprocess.run(["git", "-C", DESTINO, *a], capture_output=True, text=True)
    if r.stdout.strip(): print("   ", r.stdout.strip()[:400])
    if r.returncode != 0 and r.stderr.strip():
        print("   !", r.stderr.strip()[:300])
    return r.returncode == 0

git("init", "-q")
git("add", "-A")
git("-c", "user.name=Juan Philippe Teixeira",
    "-c", "user.email=08638028938@edu.udesc.br",
    "commit", "-q", "-m",
    "Data and analysis code for the outward-cone occupancy manuscript")
git("branch", "-M", "main")
git("remote", "add", "origin", REMOTO)

tam = sum(os.path.getsize(os.path.join(d, f))
          for d, _, fs in os.walk(DESTINO) for f in fs
          if ".git/" not in os.path.join(d, f))

print("\n" + "=" * 74)
print(f"ARQUIVOS NO REPOSITORIO ({len(copiados)} copiados, {tam/1e6:.1f} MB)")
print("=" * 74)
for rel in sorted(copiados):
    print(f"  {os.path.getsize(os.path.join(DESTINO, rel))/1024:9.0f} kB  {rel}")

print("\n" + "=" * 74)
print("REVISE A LISTA ACIMA — o repositorio e PUBLICO")
print("=" * 74)
print(f"""
Se estiver certa, suba com:

    cd {DESTINO}
    git push -u origin main

Se o repositorio remoto ja tiver commits, use antes:

    git pull --rebase origin main

Depois do push, conecte o repositorio ao Zenodo para emitir o DOI:
  zenodo.org -> GitHub -> ligar o repositorio -> criar um release no GitHub
O DOI sai do release, e e ele que entra no Data availability do manuscrito.
""")
