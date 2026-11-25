#!/usr/bin/env python3
"""Piso do filtro: calibrado em metade dos complexos, medido na outra.

O piso atual (58,6% de casca, 182 atomos bloqueadores) foi determinado
nos mesmos 6I2K e 3J8F em que foi medido. Aqui a divisao e sorteada com
semente declarada, o piso sai so da metade de calibracao, e o numero
que vale e o da metade de teste, que o piso nunca viu.

O filtro EXCLUI sitio impossivel. Entao calibra-se para SENSIBILIDADE
alta (nao rejeitar pegada real) e mede-se a ESPECIFICIDADE obtida --
quanto da superficie ele consegue descartar a esse custo.
"""
import random
import sys
from pathlib import Path

try:
    import numpy as np
except ImportError:
    sys.exit("falta numpy")

SEMENTE = 20261116
SENS_ALVO = 0.95
PISO_ATUAL = (58.6, 182)
FORA = {"23BG", "23PS"}       # sem montagem; casca de 53 A


def le(p):
    ln = Path(p).read_text().splitlines()
    cab = ln[0].split("\t")
    return [dict(zip(cab, l.split("\t"))) for l in ln[1:] if l.strip()]


complexos = []
for d in le("validacao_filtro.tsv"):
    if d["pdb"] in FORA or not d.get("auc_casca"):
        continue
    f = Path(f"{d['pdb']}_{d['cadeia']}_pegada.tsv")
    if f.exists():
        complexos.append((d["pdb"], f))
    else:
        print(f"  sem arquivo por residuo: {d['pdb']}")

print(f"{len(complexos)} complexos com dados por residuo")
random.Random(SEMENTE).shuffle(complexos)
meio = len(complexos) // 2
cal, tes = complexos[:meio], complexos[meio:]
print(f"semente {SEMENTE}")
print(f"  calibracao ({len(cal)}): {' '.join(p for p, _ in cal)}")
print(f"  teste      ({len(tes)}): {' '.join(p for p, _ in tes)}\n")


def junta(grupo):
    cas, blo, peg = [], [], []
    for _, f in grupo:
        for r in le(f):
            try:
                cas.append(float(r["shell_pct"]))
                blo.append(float(r["blocking"]))
            except ValueError:
                continue
            peg.append(r["in_footprint"].strip().lower() == "yes")
    return np.array(cas), np.array(blo), np.array(peg)


cc, cb, cp = junta(cal)
tc, tb, tp = junta(tes)
print(f"calibracao: {cp.sum()} residuos de pegada em {len(cp)} totais")
print(f"teste:      {tp.sum()} residuos de pegada em {len(tp)} totais\n")

t_casca = float(np.quantile(cc[cp], 1 - SENS_ALVO))
t_bloq = float(np.quantile(cb[cp], SENS_ALVO))
print(f"piso calibrado para {SENS_ALVO:.0%} de sensibilidade:")
print(f"  casca   >= {t_casca:.1f}%   (atual: {PISO_ATUAL[0]}%)")
print(f"  bloqueio <= {t_bloq:.0f}     (atual: {PISO_ATUAL[1]})\n")


def mede(rot, c, b, p, tc_, tb_):
    passa = (c >= tc_) & (b <= tb_)
    sens = passa[p].mean() if p.any() else float("nan")
    espec = (~passa[~p]).mean()
    print(f"  {rot:<34} sensibilidade {sens:6.1%}   "
          f"especificidade {espec:6.1%}")
    return sens, espec


print("desempenho:")
mede("calibracao, piso calibrado", cc, cb, cp, t_casca, t_bloq)
mede("TESTE, piso calibrado", tc, tb, tp, t_casca, t_bloq)
mede("TESTE, piso atual do paper", tc, tb, tp, *PISO_ATUAL)
mede("calibracao, piso atual", cc, cb, cp, *PISO_ATUAL)

print("\nsitio 22 do EV-A71 contra o piso calibrado:")
for alvo in ("3VBS", "3VBF", "3VBH"):
    for f in Path(".").glob(f"{alvo}_*_pegada.tsv"):
        for r in le(f):
            if r["residue"] == "22":
                c22, b22 = float(r["shell_pct"]), float(r["blocking"])
                print(f"  {f.name:<24} casca {c22:5.2f}% (piso {t_casca:.1f}) "
                      f"bloqueio {b22:5.0f} (piso {t_bloq:.0f})  "
                      f"{'EXCLUIDO' if (c22 < t_casca or b22 > t_bloq) else 'passa'}")
