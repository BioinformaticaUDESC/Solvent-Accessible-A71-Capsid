#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V51 - o vies de completude atomica DENTRO dos 39 retidos.

POR QUE. O descritor conta atomos. Um modelo com ~1 atomo por residuo nao e
equiparavel a um com ~7,7 a 8,0, inclusive por nao ter cadeias laterais que
ocupariam o cone. O diagnostico que eu reportei foi feito sobre o inventario
bruto de 45, que inclui complexos depois rejeitados pelo controle de montagem e
portanto NAO e a populacao das estatisticas de interesse (39 complexos, 583
residuos). A extensao do vies nos 39 e desconhecida e nao deve ser declarada
pequena por analogia.

O QUE FAZ:
  1. monta os 39 a partir das saidas da auditoria, nao de lista fixa
  2. identifica os modelos de cadeia principal pela completude ja medida
  3. compara 39 contra os de representacao mais completa: n, mediana, Q75, Q90,
     Q95, maximo, media, e contribuicao por complexo
  4. repete na populacao Picornaviridae de avaliacao externa
  5. mede o impacto na classificacao em 106 atomos, sem recalibrar nada
  6. confirma que os modelos incompletos estao fora das duas metades dos 26
  7. grava TSV com resultados, denominadores e caminho de cada arquivo

NAO recalibra o limiar. NAO roda simulacao. NAO altera nada.

  python3 V51_completude_nos_39.py
"""
import os, re, sys, csv, glob, math
import numpy as np

BASE = "/home/juan"
RAIZES = [f"{BASE}/abertos", f"{BASE}/paper_artifacts", f"{BASE}/figuras_V28", BASE]
CRITERIO = f"{BASE}/V48_CRITERIO_MONTAGEM.tsv"
TITULOS = f"{BASE}/V48_TITULOS_DEPOSITOS.tsv"
COMPLET = f"{BASE}/V48_COMPLETUDE_MODELO.tsv"
METODO = f"{BASE}/V49_METODO.txt"
PISO = 106.0          # fixo, lido do metodo de registro; nunca recalibrado aqui
MIN_APR = 5.0         # atomos por residuo abaixo disso = representacao incompleta
RES_BAIXA = 8.0       # resolucao a partir da qual a entrada e externa


def novo(p):
    if not os.path.exists(p):
        return p
    r, e = os.path.splitext(p)
    n = 2
    while os.path.exists(f"{r}_{n:02d}{e}"):
        n += 1
    return f"{r}_{n:02d}{e}"


def ler(p):
    with open(p, newline="", encoding="utf8", errors="replace") as f:
        h = f.read(4096)
        f.seek(0)
        return list(csv.DictReader(f, delimiter="\t" if h.count("\t") >= h.count(",") else ","))


def num(v):
    try:
        return float(str(v).replace(",", ".").strip())
    except Exception:
        return float("nan")


def col(cs, ch):
    for c in cs:
        if c and c.strip().lower() in ch:
            return c
    for c in cs:
        if c and any(k in c.strip().lower() for k in ch):
            return c
    return None


def verdadeiro(v):
    return str(v).strip().lower() in ("1", "true", "sim", "yes", "t", "y")


for f_ in (CRITERIO, TITULOS, COMPLET):
    if not os.path.exists(f_):
        sys.exit(f"falta {f_}; rode V48f, V48h e V48j antes")

# ---------------------------------------------- quem sai, e por qual clausula
crit = {r["pdb"].strip().upper(): r.get("decisao", "").strip() for r in ler(CRITERIO)}
sem_ops = {k for k, v in crit.items() if v != "MONTAGEM ICOSAEDRICA"}
loc = set()
reso = {}
for r in ler(TITULOS):
    pdb = r.get("pdb", "").strip().upper()
    reso[pdb] = num(r.get("resolucao", ""))
    if "REVISAR" in str(r.get("triagem", "")).upper() and pdb not in sem_ops:
        loc.add(pdb)
fora = sem_ops | loc
retidos = sorted(set(crit) - fora)
print(f"clausula 1, sem operadores ({len(sem_ops)}): {sorted(sem_ops)}")
print(f"clausula 2, reconstrucao localizada ({len(loc)}): {sorted(loc)}")
print(f"retidos: {len(retidos)}\n")

# ---------------------------------------------- completude ja medida
apr = {}
for r in ler(COMPLET):
    apr[r.get("pdb", "").strip().upper()] = num(r.get("atomos_por_residuo_todas", ""))
incompletos = sorted(p for p in retidos
                     if not math.isnan(apr.get(p, float("nan"))) and apr[p] < MIN_APR)
completos = [p for p in retidos if p not in incompletos]
print(f"representacao incompleta (< {MIN_APR:.0f} atomos/residuo): {incompletos}")
for p in incompletos:
    print(f"   {p}  {apr[p]:.2f} atomos/residuo   resolucao {reso.get(p, float('nan'))} A")
print(f"representacao completa: {len(completos)} complexos, "
      f"{min(apr[p] for p in completos if not math.isnan(apr.get(p,float('nan')))):.2f}"
      f" a {max(apr[p] for p in completos if not math.isnan(apr.get(p,float('nan')))):.2f}"
      f" atomos/residuo\n")

# ---------------------------------------------- as duas metades dos 26
cal = tes = []
if os.path.exists(METODO):
    txt = open(METODO, encoding="utf8", errors="replace").read()
    m = re.search(r"calibracao\s+([0-9A-Z \t]+)\n", txt)
    n_ = re.search(r"teste\s+([0-9A-Z \t]+)\n", txt)
    cal = m.group(1).split() if m else []
    tes = n_.group(1).split() if n_ else []
print(f"metades lidas de {os.path.basename(METODO)}: "
      f"calibracao {len(cal)}, teste {len(tes)}")
dentro = [p for p in incompletos if p in cal or p in tes]
print(f"modelos incompletos dentro das metades: {dentro or 'nenhum'}")
if not dentro and cal and tes:
    print("  -> o desempenho medido nos 26 nao e afetado por eles. Documentado.\n")
else:
    print("  -> ATENCAO: ha modelo incompleto nas metades; o desempenho muda.\n")

# ---------------------------------------------- carrega as pegadas
peg = {}
vistos = set()
for raiz in RAIZES:
    for p in glob.glob(os.path.join(raiz, "**", "*_pegada.tsv"), recursive=True):
        rp = os.path.realpath(p)
        if rp in vistos:
            continue
        vistos.add(rp)
        pdb = os.path.basename(p).split("_")[0].upper()
        if pdb not in retidos or pdb in peg:
            continue
        rows = ler(rp)
        if not rows:
            continue
        cs = list(rows[0].keys())
        cb = col(cs, {"blocking", "bloqueio", "blocking_atoms"})
        cf = col(cs, {"in_footprint", "pegada", "contato"})
        if not (cb and cf):
            continue
        cont, nao = [], []
        for r in rows:
            b = num(r.get(cb))
            if math.isnan(b):
                continue
            (cont if verdadeiro(r.get(cf)) else nao).append(b)
        peg[pdb] = {"arq": rp, "contato": cont, "nao": nao}
falta = [p for p in retidos if p not in peg]
print(f"pegadas lidas: {len(peg)} de {len(retidos)}"
      + (f"   SEM PEGADA: {falta}" if falta else ""))


def resumo(pdbs, rot):
    v = np.array([b for p in pdbs if p in peg for b in peg[p]["contato"]])
    if not len(v):
        return None
    return {"grupo": rot, "complexos": len([p for p in pdbs if p in peg]),
            "n_contato": len(v), "media": round(float(v.mean()), 1),
            "mediana": round(float(np.median(v)), 1),
            "Q75": round(float(np.percentile(v, 75)), 1),
            "Q90": round(float(np.percentile(v, 90)), 1),
            "Q95": round(float(np.percentile(v, 95)), 1),
            "maximo": int(v.max()),
            "acima_do_piso": int((v > PISO).sum()),
            "pct_acima": round(100 * float((v > PISO).sum()) / len(v), 1)}


externos = [p for p in retidos if not math.isnan(reso.get(p, float("nan")))
            and reso[p] >= RES_BAIXA]
ext_completos = [p for p in externos if p not in incompletos]

linhas = []
print("\n" + "=" * 78)
print("(1) distribuicao da ocupacao em residuos de contato, estratificada")
print("=" * 78)
print(f"{'grupo':<36}{'cx':>4}{'n':>6}{'med':>6}{'Q75':>6}{'Q90':>6}{'Q95':>6}"
      f"{'max':>6}{'>106':>7}")
for pdbs, rot in ((retidos, f"todos os {len(retidos)} retidos"),
                  (completos, f"{len(completos)} de representacao completa"),
                  (incompletos, f"{len(incompletos)} de cadeia principal"),
                  (externos, f"externos >= {RES_BAIXA:.0f} A ({len(externos)})"),
                  (ext_completos, f"externos completos ({len(ext_completos)})"),
                  ([p for p in retidos if p in cal or p in tes],
                   f"os 26 da calibracao e teste")):
    r = resumo(pdbs, rot)
    if not r:
        continue
    linhas.append(r)
    print(f"{rot:<36}{r['complexos']:>4}{r['n_contato']:>6}{r['mediana']:>6.0f}"
          f"{r['Q75']:>6.0f}{r['Q90']:>6.0f}{r['Q95']:>6.0f}{r['maximo']:>6}"
          f"{r['acima_do_piso']:>7}")

print("\n" + "=" * 78)
print(f"(2) impacto na classificacao em {PISO:.0f} atomos, sem recalibrar")
print("=" * 78)
for pdbs, rot in ((retidos, "todos os retidos"), (completos, "so completos")):
    c = np.array([b for p in pdbs if p in peg for b in peg[p]["contato"]])
    n = np.array([b for p in pdbs if p in peg for b in peg[p]["nao"]])
    if not len(c) or not len(n):
        continue
    sens = float((c <= PISO).mean())
    espe = float((n > PISO).mean())
    print(f"  {rot:<22} sensibilidade {sens:6.1%} ({int((c<=PISO).sum())}/{len(c)})"
          f"   especificidade {espe:6.1%} ({int((n>PISO).sum())}/{len(n)})")

print("\n" + "=" * 78)
print("(3) contribuicao por complexo nos retidos")
print("=" * 78)
print(f"{'pdb':<7}{'res A':>7}{'at/res':>8}{'n contato':>11}{'med':>6}{'max':>7}"
      f"{'>106':>6}  grupo")
porcx = []
for p in retidos:
    if p not in peg:
        continue
    c = peg[p]["contato"]
    g = ("cadeia principal" if p in incompletos else
         ("externo" if p in externos else "calibracao/teste" if p in cal + tes
          else "retido, fora das metades"))
    med = float(np.median(c)) if c else float("nan")
    print(f"{p:<7}{reso.get(p, float('nan')):>7.2f}{apr.get(p, float('nan')):>8.2f}"
          f"{len(c):>11}{med:>6.0f}{(max(c) if c else 0):>7.0f}"
          f"{sum(1 for x in c if x > PISO):>6}  {g}")
    porcx.append({"pdb": p, "resolucao": reso.get(p, ""),
                  "atomos_por_residuo": apr.get(p, ""), "grupo": g,
                  "n_contato": len(c), "mediana": round(med, 1) if c else "",
                  "maximo": int(max(c)) if c else "",
                  "acima_do_piso": sum(1 for x in c if x > PISO),
                  "n_nao_contato": len(peg[p]["nao"]), "arquivo": peg[p]["arq"]})

s1 = novo(f"{BASE}/V51_ESTRATIFICADO.tsv")
with open(s1, "w", newline="", encoding="utf8") as f:
    w = csv.DictWriter(f, fieldnames=list(linhas[0].keys()), delimiter="\t",
                       extrasaction="ignore")
    w.writeheader()
    for r in linhas:
        w.writerow(r)
s2 = novo(f"{BASE}/V51_POR_COMPLEXO.tsv")
with open(s2, "w", newline="", encoding="utf8") as f:
    w = csv.DictWriter(f, fieldnames=list(porcx[0].keys()), delimiter="\t",
                       extrasaction="ignore")
    w.writeheader()
    for r in porcx:
        w.writerow(r)
print(f"\n-> {s1}\n-> {s2}")
print(f"\nO piso de {PISO:.0f} foi usado como constante. Nenhuma recalibracao.")
print("Compare a linha dos retidos com a dos de representacao completa: a")
print("diferenca entre elas e o tamanho do vies NESTA populacao, que e a que")
print("vale. O diagnostico anterior, sobre os 45 brutos, nao serve para isso.")
