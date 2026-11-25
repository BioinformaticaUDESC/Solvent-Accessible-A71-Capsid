#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V54 - estimativa ciente da selecao, e o denominador 285 contra 293.

PARTE A - o vies de selecao da regra. O limiar de 106 saiu so da metade de
calibracao, mas a ESCOLHA de qual regra adotar (so radial, so cone, as duas)
foi orientada pelo desempenho na metade retida. Independencia do limiar e
independencia do modelo sao coisas diferentes e aqui so a primeira vale. Isto
mede a segunda por reamostragem aninhada POR COMPLEXO, usando so os dados que
ja existem:

  laco externo   deixa um dos 26 complexos de fora
  laco interno   nos 25 restantes: calibra o limiar de cada eixo no quantil de
                 sensibilidade declarado E escolhe a regra por desempenho
                 interno, sem olhar o complexo de fora
  avaliacao      aplica a regra e o limiar escolhidos ao complexo de fora

A regra de escolha interna tem de ser declarada, e esta e: entre as regras que
atingem a sensibilidade alvo no treino interno, a de maior especificidade; se
nenhuma atingir, a de maior sensibilidade. Formalizar isso agora e em si uma
limitacao, menor do que escolher a olho na metade retida, e esta dito.

PARTE B - o denominador da correlacao de SASA. O V28_14 reporta 285 residuos
com area isolada mensuravel; o V28_15 correlaciona retencao em 293. Para area
isolada nula ou quase nula a razao e indefinida ou instavel. Aqui se conta
quantos dos 293 tem area isolada pequena e se recalcula a correlacao sem eles,
para ver se o valor se move.

NAO roda simulacao. NAO recalibra o limiar publicado. Somente leitura.

  python3 V54_selecao_e_denominador.py
"""
import os, re, csv, glob, math, sys
import numpy as np

BASE = "/home/juan"
RAIZES = [f"{BASE}/abertos", f"{BASE}/paper_artifacts", f"{BASE}/figuras_V28", BASE]
METODO = f"{BASE}/V49_METODO.txt"
SENS = 0.95
PISO_PUB = 106.0
CORTES = (0.0, 1.0, 5.0)      # area isolada minima a testar, em A^2


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


def vdd(v):
    return str(v).strip().lower() in ("1", "true", "sim", "yes", "t", "y")


def spearman(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    m = ~(np.isnan(a) | np.isnan(b))
    a, b = a[m], b[m]
    if len(a) < 3:
        return float("nan"), 0

    def rk(x):
        o = np.argsort(x, kind="mergesort")
        r = np.empty(len(x), float)
        r[o] = np.arange(1, len(x) + 1)
        for v in np.unique(x):
            i = x == v
            r[i] = r[i].mean()
        return r
    return float(np.corrcoef(rk(a), rk(b))[0, 1]), len(a)


# ===================================================== PARTE A
print("=" * 78)
print("(A) estimativa ciente da selecao, reamostragem aninhada por complexo")
print("=" * 78)
if not os.path.exists(METODO):
    sys.exit(f"falta {METODO}; rode o V49 antes")
txt = open(METODO, encoding="utf8", errors="replace").read()
cal = (re.search(r"calibracao\s+([0-9A-Z \t]+)\n", txt) or [None, ""])[1].split()
tes = (re.search(r"teste\s+([0-9A-Z \t]+)\n", txt) or [None, ""])[1].split()
os26 = cal + tes
print(f"os 26 lidos de {os.path.basename(METODO)}: {len(os26)}")

dados, vistos = {}, set()
for r in RAIZES:
    for p in glob.glob(os.path.join(r, "**", "*_pegada.tsv"), recursive=True):
        rp = os.path.realpath(p)
        if rp in vistos:
            continue
        vistos.add(rp)
        pdb = os.path.basename(p).split("_")[0].upper()
        if pdb not in os26 or pdb in dados:
            continue
        rows = ler(rp)
        if not rows:
            continue
        cs = list(rows[0].keys())
        cb = col(cs, {"blocking", "bloqueio", "blocking_atoms"})
        ck = col(cs, {"shell_pct", "casca_pct"})
        cf = col(cs, {"in_footprint", "pegada", "contato"})
        if not (cb and ck and cf):
            continue
        v = []
        for row in rows:
            b, k = num(row.get(cb)), num(row.get(ck))
            if math.isnan(b) or math.isnan(k):
                continue
            v.append((k, b, vdd(row.get(cf))))
        dados[pdb] = v
print(f"pegadas com casca e ocupacao: {len(dados)} de {len(os26)}")
falta = [p for p in os26 if p not in dados]
if falta:
    print(f"  sem dados: {falta}")


def q_np(v, q):
    return float(np.quantile(np.asarray(v, float), q))


def pisos(treino):
    """limiares dos dois eixos, pelo quantil de sensibilidade, so no treino"""
    peg = [(k, b) for p in treino for (k, b, f) in dados[p] if f]
    if not peg:
        return None, None
    return q_np([k for k, _ in peg], 1 - SENS), q_np([b for _, b in peg], SENS)


REGRAS = {
    "so cone": lambda k, b, ts, tb: b <= tb,
    "so radial": lambda k, b, ts, tb: k >= ts,
    "as duas": lambda k, b, ts, tb: (k >= ts) and (b <= tb),
}


def avalia(pdbs, regra, ts, tb):
    f = REGRAS[regra]
    c = [f(k, b, ts, tb) for p in pdbs for (k, b, fp) in dados[p] if fp]
    n = [f(k, b, ts, tb) for p in pdbs for (k, b, fp) in dados[p] if not fp]
    if not c or not n:
        return None
    return (sum(c) / len(c), 1 - sum(n) / len(n), len(c), len(n))


def escolhe(treino, ts, tb):
    """regra de escolha DECLARADA: maior especificidade entre as que atingem o
    alvo de sensibilidade; se nenhuma atingir, maior sensibilidade"""
    cand = []
    for rot in REGRAS:
        m = avalia(treino, rot, ts, tb)
        if m:
            cand.append((rot, m[0], m[1]))
    if not cand:
        return None
    ok = [c for c in cand if c[1] >= SENS]
    return max(ok, key=lambda c: c[2])[0] if ok else max(cand, key=lambda c: c[1])[0]


pres = [p for p in os26 if p in dados]
escolhas, acertos_c, acertos_n = {}, [], []
print(f"\n{'fora':<8}{'regra escolhida':<14}{'piso cone':>11}{'piso radial':>13}"
      f"{'sens':>8}{'espec':>8}")
linhas = []
for fora in pres:
    treino = [p for p in pres if p != fora]
    ts, tb = pisos(treino)
    if ts is None:
        continue
    rot = escolhe(treino, ts, tb)
    m = avalia([fora], rot, ts, tb)
    escolhas[rot] = escolhas.get(rot, 0) + 1
    f = REGRAS[rot]
    acertos_c += [f(k, b, ts, tb) for (k, b, fp) in dados[fora] if fp]
    acertos_n += [f(k, b, ts, tb) for (k, b, fp) in dados[fora] if not fp]
    s = f"{m[0]:>8.1%}{m[1]:>8.1%}" if m else f"{'-':>8}{'-':>8}"
    print(f"{fora:<8}{rot:<14}{tb:>11.0f}{ts:>13.1f}{s}")
    linhas.append({"deixado_fora": fora, "regra_escolhida": rot,
                   "piso_cone": round(tb, 1), "piso_radial": round(ts, 1),
                   "sens_no_complexo": round(m[0]*100, 1) if m else "",
                   "espec_no_complexo": round(m[1]*100, 1) if m else "",
                   "n_contato": sum(1 for x in dados[fora] if x[2]),
                   "n_nao_contato": sum(1 for x in dados[fora] if not x[2])})

if acertos_c:
    S = sum(acertos_c) / len(acertos_c)
    E = 1 - sum(acertos_n) / len(acertos_n)
    print(f"\nagregado sobre as {len(linhas)} particoes externas:")
    print(f"  sensibilidade  {S:6.1%}  ({sum(acertos_c)}/{len(acertos_c)})")
    print(f"  especificidade {E:6.1%}  ({len(acertos_n)-sum(acertos_n)}/{len(acertos_n)})")
    print(f"  regra escolhida internamente: "
          + ", ".join(f"{k} em {v} de {len(linhas)}" for k, v in sorted(escolhas.items())))
    print(f"\n  comparacao: a estimativa publicada, com a regra escolhida olhando a")
    print(f"  metade retida, e 96,5 % / 65,3 %. A diferenca entre as duas e o")
    print(f"  tamanho do vies de selecao do modelo nesta amostra.")

# ===================================================== PARTE B
print("\n" + "=" * 78)
print("(B) denominador da correlacao de SASA: 285 contra 293")
print("=" * 78)
alvo = None
for r in RAIZES:
    g = glob.glob(os.path.join(r, "**", "V28_SASA_vs_RADIAL*"), recursive=True)
    if g:
        alvo = sorted(g, key=len)[0]
        break
sasa = None
for r in RAIZES:
    g = [x for x in glob.glob(os.path.join(r, "**", "V28_SASA_3VBS*"), recursive=True)
         if x.endswith((".tsv", ".csv"))]
    if g:
        sasa = sorted(g, key=len)[0]
        break
for rot, p in (("cruzado", alvo), ("SASA", sasa)):
    print(f"  {rot:<9}{p if p else 'nao localizado'}")
if not alvo:
    print("\n  sem o arquivo cruzado nao da para refazer a correlacao")
else:
    rows = ler(alvo)
    cs = list(rows[0].keys()) if rows else []
    print(f"  colunas: {cs}")
    ciso = col(cs, {"sasa_isolada", "sasa_isolated", "isolada", "iso_a2"})
    cret = col(cs, {"retencao", "retained", "fracao", "fraction", "ret"})
    ccas = col(cs, {"shell_pct", "casca", "casca_pct", "percentil"})
    print(f"  usando: isolada={ciso}  retencao={cret}  casca={ccas}")
    if not (cret and ccas):
        print("  colunas insuficientes; abra o arquivo e me diga os nomes")
    else:
        tot = len(rows)
        for corte in CORTES:
            sel = []
            for r_ in rows:
                iso = num(r_.get(ciso)) if ciso else float("nan")
                ret, cas = num(r_.get(cret)), num(r_.get(ccas))
                if math.isnan(ret) or math.isnan(cas):
                    continue
                if ciso and not math.isnan(iso) and iso <= corte:
                    continue
                sel.append((ret, cas))
            rho, n = spearman([x[0] for x in sel], [x[1] for x in sel])
            print(f"  area isolada > {corte:>4.1f} A^2 : n = {n:>4}   Spearman {rho:+.3f}")
        print(f"  linhas no arquivo: {tot}")
        print("\n  Se o n cair de 293 para perto de 285 e o rho nao se mover, a")
        print("  correlacao publicada nao depende dos residuos de area quase nula,")
        print("  e o texto deve declarar o corte usado. Se mover, o valor muda.")

if linhas:
    s = novo(f"{BASE}/V54_SELECAO_ANINHADA.tsv")
    with open(s, "w", newline="", encoding="utf8") as f:
        w = csv.DictWriter(f, fieldnames=list(linhas[0].keys()), delimiter="\t")
        w.writeheader()
        for r_ in linhas:
            w.writerow(r_)
    print(f"\n-> {s}")
