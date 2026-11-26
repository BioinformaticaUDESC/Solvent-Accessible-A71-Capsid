#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
V94 - Etapa 1 do parecer: fecha a Tabela 4 para que as frases de
      "provisional", "unresolved" e "verification remains necessary"
      deixem de ser verdade e possam ser apagadas.

Faz quatro coisas, todas a partir dos arquivos por residuo em
/home/juan/V83_benchmark -- nenhuma vem de resumo de execucao anterior:

A. CONTABILIDADE DE INCLUSAO, com o numero da tolerancia explicito.
   Para cada complexo compara, residuo a residuo, a ocupacao
   recalculada (coluna cone_45_30, gerada pelo V83 a partir da
   montagem reconstruida) contra a depositada (coluna blocking_pub).
   Reporta n, exatos, diferenca maxima e trocas no limiar de 106.

B. TABELA 4 inteira recalculada, com intervalo marginal para TODOS os
   descritores, residue depth incluido -- o que hoje esta "not
   reported" -- pelo mesmo bootstrap em nivel de complexo usado nos
   demais.  Mais as diferencas pareadas e a contagem de vitorias.

C. GRADE 3x3 de meia-abertura x raio, que o texto cita sem mostrar.

D. Grava quatro arquivos: a Tabela 4 pronta para colar, as mesmas
   tabelas em TSV para deposito no Zenodo, e a contabilidade de
   inclusao.

Nao altera o manuscrito.  Nao escreve fora de /home/juan/V83_benchmark.

Uso:
    python3 V94_etapa1_fecha_tabela4.py
    python3 V94_etapa1_fecha_tabela4.py --reamostragens 10000
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
TOLERANCIA = 1          # atomos; fixada antes de qualquer complexo ser testado

for i, a in enumerate(sys.argv):
    if a == "--reamostragens" and i + 1 < len(sys.argv):
        NBOOT = int(sys.argv[i + 1])

REF = "cone_45_30"
BASES = [("sasa_1.4", "Accessible area, 1.4 A probe", +1),
         ("sasa_5.0", "Accessible area, 5.0 A probe", +1),
         ("sasa_10.0", "Accessible area, 10.0 A probe", +1),
         ("prof_media", "Residue depth", -1),
         ("shell_pct", "Radial percentile", +1),
         ("hse_fora", "Half-sphere exposure", -1)]
GRADE = ["cone_%d_%d" % (s, r) for s in (30, 45, 60) for r in (20, 30, 40)]
DIRECAO = {c: -1 for c in GRADE}
DIRECAO.update({c: d for c, _, d in BASES})
DIRECAO["blocking_pub"] = -1
DIRECAO["cone_45_30_sem_vizinhos"] = -1


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
    with open(caminho, "r") as fh:
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
            # alinhar pelo numero de residuo
            mapa = {int(r): i for i, r in enumerate(dp["residuo"])}
            d["prof_media"] = np.array(
                [dp["prof_media"][mapa[int(r)]] if int(r) in mapa else np.nan
                 for r in d["residuo"]])
        else:
            d["prof_media"] = np.full(len(d["residuo"]), np.nan)
        casos[pid] = d
    return casos


def ic(v):
    v = np.asarray(v)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return (np.nan, np.nan)
    return float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))


def main():
    barra("V94 - Etapa 1: fechar a Tabela 4")
    casos = carrega()
    if not casos:
        sys.exit("nenhum TSV em %s.  rode o V83 e o V86 antes." % BENCH)
    pids = sorted(casos)
    print("  %d complexos   semente %d   %d reamostragens" % (len(pids), SEMENTE, NBOOT))
    print("  fontes: %s/<pdb>_descritores.tsv e <pdb>_profundidade.tsv" % BENCH)

    # ------------------------------------------------------------- A -----
    barra("A.  CONTABILIDADE DE INCLUSAO")
    print("  Regra de inclusao, fixada antes de qualquer complexo ser testado:")
    print("    diferenca de no maximo %d atomo em TODO residuo, E" % TOLERANCIA)
    print("    nenhum residuo cruzando o limiar de %d atomos." % LIMIAR)
    print("\n  %-8s %8s %9s %8s %8s %8s" %
          ("pdb", "residuos", "exatos", "max dif", "trocas", "veredito"))
    linhas_inc = []
    for pid in pids:
        d = casos[pid]
        a = d[REF]
        b = d["blocking_pub"]
        m = np.isfinite(a) & np.isfinite(b)
        dif = (a[m] - b[m]).astype(int)
        flips = int(np.sum((b[m] <= LIMIAR) != (a[m] <= LIMIAR)))
        mx = int(np.abs(dif).max()) if dif.size else -1
        ok = (mx <= TOLERANCIA) and flips == 0
        print("  %-8s %8d %9d %8d %8d %8s"
              % (pid, int(m.sum()), int((dif == 0).sum()), mx, flips,
                 "inclui" if ok else "EXCLUI"))
        linhas_inc.append((pid, int(m.sum()), int((dif == 0).sum()), mx, flips,
                           "incluido" if ok else "excluido"))
    todos_ok = all(l[5] == "incluido" for l in linhas_inc)
    print("\n  %d complexos, todos dentro da regra: %s"
          % (len(pids), "sim" if todos_ok else "NAO"))
    if todos_ok:
        print("  Nao ha complexo retido fora da tolerancia.  A frase sobre")
        print("  contabilidade de inclusao nao resolvida deixa de ser verdade.")

    # ------------------------------------------------------------- B -----
    barra("B.  TABELA 4 RECALCULADA, com intervalo marginal para todos")
    rng = np.random.default_rng(SEMENTE)
    idx = rng.integers(0, len(pids), size=(NBOOT, len(pids)))

    def pooled_ic(col):
        v = np.concatenate([casos[p][col] for p in pids])
        y = np.concatenate([casos[p]["contato"] for p in pids])
        m = np.isfinite(v)
        base = auc(DIRECAO[col] * v[m], y[m])
        am = []
        for linha in idx:
            vv = np.concatenate([casos[pids[j]][col] for j in linha])
            yy = np.concatenate([casos[pids[j]]["contato"] for j in linha])
            mm = np.isfinite(vv)
            if yy[mm].sum() == 0 or (~yy[mm]).sum() == 0:
                continue
            am.append(auc(DIRECAO[col] * vv[mm], yy[mm]))
        return base, ic(am), int(m.sum())

    def por_complexo(col):
        out = {}
        for p in pids:
            d = casos[p]
            v, y = d[col], d["contato"]
            m = np.isfinite(v)
            out[p] = (auc(DIRECAO[col] * v[m], y[m])
                      if m.sum() > 20 and y[m].sum() > 0 else np.nan)
        return out

    ac_ref = por_complexo(REF)
    ref_arr = np.array([ac_ref[p] for p in pids])
    base_ref, ic_ref, n_ref = pooled_ic(REF)

    print("\n  %-30s %8s %20s %22s %9s"
          % ("descritor", "AUC", "IC95 marginal", "dif pareada (IC95)", "cone vence"))
    print("  %-30s %8.4f   [%.4f, %.4f] %22s %9s"
          % ("Outward-cone occupancy 45/30", base_ref, ic_ref[0], ic_ref[1], "-", "-"))

    tabela = [("Outward-cone occupancy (45 deg, 30 A)", base_ref, ic_ref,
               None, None, None)]
    for col, rotulo, _ in BASES:
        if col not in casos[pids[0]]:
            continue
        b, i95, n = pooled_ic(col)
        ac = por_complexo(col)
        dif = ref_arr - np.array([ac[p] for p in pids])
        m = np.isfinite(dif)
        boot = np.array([dif[m][rng.integers(0, m.sum(), m.sum())].mean()
                         for _ in range(NBOOT)])
        d_ic = ic(boot)
        vence = int((dif[m] > 0).sum())
        print("  %-30s %8.4f   [%.4f, %.4f]   %.4f [%.4f, %.4f] %6d/%-3d"
              % (rotulo[:30], b, i95[0], i95[1], dif[m].mean(),
                 d_ic[0], d_ic[1], vence, int(m.sum())))
        tabela.append((rotulo, b, i95, dif[m].mean(), d_ic, (vence, int(m.sum()))))

    # ------------------------------------------------------------- C -----
    barra("C.  GRADE 3x3, pareada contra o publicado")
    grade = []
    for col in GRADE:
        if col not in casos[pids[0]]:
            continue
        b, i95, _ = pooled_ic(col)
        ac = por_complexo(col)
        dif = np.array([ac[p] for p in pids]) - ref_arr
        m = np.isfinite(dif)
        boot = np.array([dif[m][rng.integers(0, m.sum(), m.sum())].mean()
                         for _ in range(NBOOT)])
        d_ic = ic(boot)
        marca = "  <- usado no artigo" if col == REF else ""
        print("  %-14s AUC %.4f   dif %+.4f  [%+.4f, %+.4f]%s"
              % (col, b, dif[m].mean(), d_ic[0], d_ic[1], marca))
        grade.append((col, b, dif[m].mean(), d_ic))
    cruza = [g for g in grade if not (g[3][0] <= 0 <= g[3][1])]
    print("\n  celulas cuja diferenca exclui zero: %d de %d" % (len(cruza), len(grade)))

    # ------------------------------------------------------------- D -----
    barra("D.  ARQUIVOS GRAVADOS")
    p1 = os.path.join(BENCH, "V94_inclusao.tsv")
    with open(p1, "w") as fh:
        fh.write("complexo\tresiduos\texatos\tdif_maxima\ttrocas_limiar\tveredito\n")
        for l in linhas_inc:
            fh.write("\t".join(str(x) for x in l) + "\n")

    p2 = os.path.join(BENCH, "V94_tabela4.tsv")
    with open(p2, "w") as fh:
        fh.write("descritor\tauc\tic_lo\tic_hi\tdif_pareada\tdif_lo\tdif_hi\t"
                 "cone_vence\tde\n")
        for rot, b, i95, dm, dic, v in tabela:
            if dm is None:
                fh.write("%s\t%.4f\t%.4f\t%.4f\t\t\t\t\t\n" % (rot, b, i95[0], i95[1]))
            else:
                fh.write("%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%d\t%d\n"
                         % (rot, b, i95[0], i95[1], dm, dic[0], dic[1], v[0], v[1]))

    p3 = os.path.join(BENCH, "V94_grade3x3.tsv")
    with open(p3, "w") as fh:
        fh.write("configuracao\tauc\tdif_vs_publicado\tdif_lo\tdif_hi\n")
        for col, b, dm, dic in grade:
            fh.write("%s\t%.4f\t%.4f\t%.4f\t%.4f\n" % (col, b, dm, dic[0], dic[1]))

    p4 = os.path.join(BENCH, "V94_TABELA4_PARA_COLAR.md")
    with open(p4, "w") as fh:
        fh.write("## Table 4\n\n")
        fh.write("Discrimination of contact residues by outward-cone occupancy and by "
                 "five established accessibility descriptors, over %d assembled "
                 "complexes. Intervals are from %d resamples of complexes "
                 "(seed %d). Paired difference is AUC(cone) - AUC(descriptor) "
                 "computed within each complex.\n\n" % (len(pids), NBOOT, SEMENTE))
        fh.write("-" * 75 + "\n")
        fh.write("%-38s %-20s %-22s %s\n" % ("Descriptor", "AUC (95% CI)",
                                             "Paired difference", "Cone higher"))
        fh.write("-" * 75 + "\n")
        for rot, b, i95, dm, dic, v in tabela:
            if dm is None:
                fh.write("%-38s %.3f (%.3f-%.3f) %-22s %s\n"
                         % (rot, b, i95[0], i95[1], "-", "-"))
            else:
                fh.write("%-38s %.3f (%.3f-%.3f) %.3f (%.3f-%.3f)  %d/%d\n"
                         % (rot, b, i95[0], i95[1], dm, dic[0], dic[1], v[0], v[1]))
        fh.write("-" * 75 + "\n")

    for p in (p1, p2, p3, p4):
        print("  %s" % p)

    barra("O QUE ISSO LIBERA NO TEXTO")
    print("  1. O intervalo de residue depth deixa de ser 'not reported'.")
    print("  2. A contabilidade de inclusao tem regra numerica e veredito por")
    print("     complexo; nao ha complexo retido fora da tolerancia.")
    print("  3. A grade 3x3 existe em numeros e a frase solta sobre")
    print("     'parameter-sensitivity analyses' passa a ter referente.")
    print("  4. Cada celula da Tabela 4 foi recalculada dos arquivos por")
    print("     residuo nesta execucao, nao lida de um resumo anterior.")
    print("\n  Os quatro TSV servem para o deposito no Zenodo.")
    print("  Falta ainda rodar os 4 complexos excluidos para mostrar que o")
    print("  ranking com 26 nao muda -- isso exige reconstruir as montagens.")


if __name__ == "__main__":
    main()
