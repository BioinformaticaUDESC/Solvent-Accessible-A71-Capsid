#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
V98 - Tabela 4 nos 22 incluidos e nos 26 todos, lado a lado, para
      responder ao pedido de mostrar que o ranking nao muda.

A separacao entre incluido e excluido NAO vem de memoria nem de lista
escrita: e recalculada aqui, residuo a residuo, comparando a ocupacao
recalculada (cone_45_30) com a depositada (blocking_pub), pela mesma
regra fixada antes de qualquer complexo ser testado -- no maximo um
atomo de diferenca em todo residuo e nenhuma travessia do limiar 106.

Imprime as duas tabelas, a ordem dos descritores em cada uma, e diz
explicitamente se alguma posicao trocou e se algum intervalo pareado
passou a incluir zero.

Nao altera nada.  Grava /home/juan/V83_benchmark/V98_tabela4_26.tsv.

Uso:
    python3 V98_tabela4_com26.py
"""

import os
import sys
import glob
import numpy as np

try:
    from scipy.stats import rankdata
except ImportError:
    sys.exit("scipy necessario.")

BENCH = "/home/juan/V83_benchmark"
SEMENTE = 20261125
NBOOT = 5000
LIMIAR = 106
TOLERANCIA = 1
REF = "cone_45_30"

BASES = [("sasa_1.4", "Accessible area, 1.4 A probe", +1),
         ("sasa_5.0", "Accessible area, 5.0 A probe", +1),
         ("sasa_10.0", "Accessible area, 10.0 A probe", +1),
         ("prof_media", "Residue depth", -1),
         ("shell_pct", "Radial percentile", +1),
         ("hse_fora", "Half-sphere exposure", -1)]
DIRECAO = {c: d for c, _, d in BASES}
DIRECAO[REF] = -1
DIRECAO["blocking_pub"] = -1


def barra(t):
    print("\n" + "=" * 78)
    print(t)
    print("=" * 78)


def auc(score, y):
    y = np.asarray(y, dtype=bool)
    n1, n0 = int(y.sum()), int((~y).sum())
    if n1 == 0 or n0 == 0:
        return np.nan
    r = rankdata(np.asarray(score, dtype=float))
    return (r[y].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0)


def le_tsv(caminho):
    with open(caminho) as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        br = {c: [] for c in cols}
        for linha in fh:
            p = linha.rstrip("\n").split("\t")
            if len(p) != len(cols):
                continue
            for c, v in zip(cols, p):
                br[c].append(v)
    out = {}
    for c, vs in br.items():
        if c == "contato":
            out[c] = np.array([v.strip().lower() == "yes" for v in vs], dtype=bool)
        else:
            out[c] = np.array([np.nan if v == "NA" else float(v) for v in vs])
    return out


def carrega():
    casos = {}
    for a in sorted(glob.glob(os.path.join(BENCH, "*_descritores.tsv"))):
        pid = os.path.basename(a).split("_")[0].upper()
        d = le_tsv(a)
        p = os.path.join(BENCH, "%s_profundidade.tsv" % pid)
        if os.path.exists(p):
            dp = le_tsv(p)
            mapa = {int(r): i for i, r in enumerate(dp["residuo"])}
            d["prof_media"] = np.array(
                [dp["prof_media"][mapa[int(r)]] if int(r) in mapa else np.nan
                 for r in d["residuo"]])
        else:
            d["prof_media"] = np.full(len(d["residuo"]), np.nan)
        casos[pid] = d
    return casos


def veredito(d):
    a, b = d[REF], d["blocking_pub"]
    m = np.isfinite(a) & np.isfinite(b)
    dif = (a[m] - b[m]).astype(int)
    flips = int(np.sum((b[m] <= LIMIAR) != (a[m] <= LIMIAR)))
    mx = int(np.abs(dif).max()) if dif.size else 99
    return (mx <= TOLERANCIA and flips == 0), mx, flips


def ic(v):
    v = np.asarray(v)
    v = v[np.isfinite(v)]
    return (np.nan, np.nan) if v.size == 0 else (float(np.percentile(v, 2.5)),
                                                 float(np.percentile(v, 97.5)))


def tabela(casos, pids, rotulo):
    rng = np.random.default_rng(SEMENTE)
    idx = rng.integers(0, len(pids), size=(NBOOT, len(pids)))

    def pooled(col):
        v = np.concatenate([casos[p][col] for p in pids])
        y = np.concatenate([casos[p]["contato"] for p in pids])
        m = np.isfinite(v)
        base = auc(DIRECAO[col] * v[m], y[m])
        am = []
        for linha in idx:
            vv = np.concatenate([casos[pids[j]][col] for j in linha])
            yy = np.concatenate([casos[pids[j]]["contato"] for j in linha])
            mm = np.isfinite(vv)
            if yy[mm].sum() and (~yy[mm]).sum():
                am.append(auc(DIRECAO[col] * vv[mm], yy[mm]))
        return base, ic(am)

    def por_cplx(col):
        o = {}
        for p in pids:
            v, y = casos[p][col], casos[p]["contato"]
            m = np.isfinite(v)
            o[p] = (auc(DIRECAO[col] * v[m], y[m])
                    if m.sum() > 20 and y[m].sum() else np.nan)
        return o

    ref_arr = np.array([por_cplx(REF)[p] for p in pids])
    b_ref, i_ref = pooled(REF)
    print("\n  %s  (%d complexos, %d residuos, %d contatos)"
          % (rotulo, len(pids),
             sum(len(casos[p]["residuo"]) for p in pids),
             int(sum(casos[p]["contato"].sum() for p in pids))))
    print("  %-30s %8s %18s %22s %8s"
          % ("descritor", "AUC", "IC95", "dif pareada (IC95)", "vence"))
    print("  %-30s %8.4f  [%.4f, %.4f] %22s %8s"
          % ("Outward-cone occupancy 45/30", b_ref, i_ref[0], i_ref[1], "-", "-"))
    linhas = [("Outward-cone occupancy (45 deg, 30 A)", b_ref, i_ref, None, None, None)]
    for col, rot, _ in BASES:
        if col not in casos[pids[0]]:
            continue
        b, i95 = pooled(col)
        ac = por_cplx(col)
        dif = ref_arr - np.array([ac[p] for p in pids])
        m = np.isfinite(dif)
        boot = np.array([dif[m][rng.integers(0, m.sum(), m.sum())].mean()
                         for _ in range(NBOOT)])
        d_ic = ic(boot)
        vence = int((dif[m] > 0).sum())
        cruza = d_ic[0] <= 0 <= d_ic[1]
        print("  %-30s %8.4f  [%.4f, %.4f]   %.4f [%.4f, %.4f] %4d/%-3d%s"
              % (rot[:30], b, i95[0], i95[1], dif[m].mean(), d_ic[0], d_ic[1],
                 vence, int(m.sum()), "  <- inclui zero" if cruza else ""))
        linhas.append((rot, b, i95, dif[m].mean(), d_ic, (vence, int(m.sum()))))
    return linhas


def main():
    barra("V98 - Tabela 4 com 22 e com 26")
    casos = carrega()
    if not casos:
        sys.exit("nenhum TSV em %s." % BENCH)

    incl, excl = [], []
    print("\n  %-8s %8s %8s %10s" % ("pdb", "max dif", "trocas", "veredito"))
    for pid in sorted(casos):
        ok, mx, fl = veredito(casos[pid])
        print("  %-8s %8d %8d %10s" % (pid, mx, fl, "inclui" if ok else "EXCLUI"))
        (incl if ok else excl).append(pid)
    print("\n  incluidos %d   fora da tolerancia %d  (%s)"
          % (len(incl), len(excl), ", ".join(excl) if excl else "-"))

    if not excl:
        print("\n  Os quatro ainda nao foram reconstruidos e medidos.")
        print("  Rode o V97, depois V83 --so 6ADR,6EIT,7C9W,8B8R, depois V86.")
        return

    t22 = tabela(casos, incl, "A.  ANALISE PRINCIPAL, so os incluidos")
    t26 = tabela(casos, sorted(casos), "B.  ROBUSTEZ, todos os complexos")

    barra("C.  A ORDEM MUDA?")
    o22 = [r[0] for r in sorted(t22[1:], key=lambda r: -r[1])]
    o26 = [r[0] for r in sorted(t26[1:], key=lambda r: -r[1])]
    print("  ordem por AUC, %d complexos:" % len(incl))
    for i, r in enumerate(o22, 1):
        print("    %d. %s" % (i, r))
    print("\n  ordem por AUC, %d complexos:" % len(casos))
    for i, r in enumerate(o26, 1):
        print("    %d. %s" % (i, r))
    print("\n  >>> ordem identica: %s" % ("SIM" if o22 == o26 else "NAO"))

    perdeu = [r[0] for r in t26[1:] if r[4] and r[4][0] <= 0 <= r[4][1]]
    print("  >>> intervalos pareados que passam a incluir zero com %d: %s"
          % (len(casos), ", ".join(perdeu) if perdeu else "nenhum"))

    p = os.path.join(BENCH, "V98_tabela4_26.tsv")
    with open(p, "w") as fh:
        fh.write("conjunto\tdescritor\tauc\tic_lo\tic_hi\tdif\tdif_lo\tdif_hi\tvence\tde\n")
        for rot_conj, tt in (("incluidos_%d" % len(incl), t22),
                             ("todos_%d" % len(casos), t26)):
            for r in tt:
                if r[3] is None:
                    fh.write("%s\t%s\t%.4f\t%.4f\t%.4f\t\t\t\t\t\n"
                             % (rot_conj, r[0], r[1], r[2][0], r[2][1]))
                else:
                    fh.write("%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%d\t%d\n"
                             % (rot_conj, r[0], r[1], r[2][0], r[2][1],
                                r[3], r[4][0], r[4][1], r[5][0], r[5][1]))
    barra("saida")
    print("  %s" % p)
    print("\n  Para o texto: a analise principal continua sendo a dos %d" % len(incl))
    print("  incluidos.  A tabela com %d entra como checagem de robustez," % len(casos))
    print("  nao como resultado alternativo -- os quatro estao fora porque a")
    print("  reconstrucao nao os reproduziu, e isso nao muda por conveniencia.")


if __name__ == "__main__":
    main()
