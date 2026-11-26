#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
V101 - Recupera a convencao da divisao e do quantil, e so entao roda a
       grade do ponto de operacao.

Por que: o V100 derivou a divisao com np.random.default_rng (PCG64) e
obteve metades que nao batem com as documentadas -- 8 dos 13 em comum,
nem igual nem complementar.  Mesma semente, gerador diferente, outra
permutacao.  Com a divisao errada, o limiar sai 109 em vez de 106 e a
grade inteira perde o sentido.

Isto NAO e escolher parametro para melhorar resultado.  O alvo e
publicado e fixo: a lista nominal dos 13 complexos retidos, que esta
nos registros do projeto, e o limiar 106 com 96.5% de retencao e 65.3%
de rejeicao.  O script procura qual combinacao de convencoes reproduz
esse alvo, e se nenhuma reproduzir, ele diz isso e usa a lista nominal
diretamente, declarando que a permutacao nao foi recuperada.

Convencoes testadas:
  ordenacao   -- ordem da validacao_filtro.tsv com as exclusoes antes da
                 permutacao, com as exclusoes depois, e ordem alfabetica
  gerador     -- default_rng (PCG64), RandomState (MT19937), e o shuffle
                 do modulo random
  metade      -- primeira metade como calibracao ou como retida
  quantil     -- indice piso ou interpolacao linear

Grava /home/juan/V83_benchmark/V101_grade_operacao.tsv.

Uso:
    python3 V101_recupera_convencao.py
"""

import os
import sys
import glob
import random
import itertools
import numpy as np

BENCH = "/home/juan/V83_benchmark"
SEMENTE = 20261116
QUANTIL = 95.0
LIMIAR_PUB, RETEN_PUB, REJEI_PUB = 106, 96.5, 65.3
EXCLUIDOS = {"23BG", "23PS"}

RETIDO_REGISTRO = {"1V9U", "6CX1", "6I2K", "6PPO", "6PSF", "6SNW", "7VXZ",
                   "7VY5", "7VY6", "7VYL", "7VYM", "7XXA", "8B8R"}

SEMIANGS = [30, 45, 60]
RAIOS = [20, 30, 40]


def barra(t):
    print("\n" + "=" * 78)
    print(t)
    print("=" * 78)


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


def ordem_validacao():
    cam = None
    for p in ("/home/juan/abertos/validacao_filtro.tsv",
              "/home/juan/validacao_filtro.tsv"):
        if os.path.exists(p):
            cam = p
            break
    if cam is None:
        ach = glob.glob("/home/juan/**/validacao_filtro.tsv", recursive=True)
        cam = ach[0] if ach else None
    if cam is None:
        return None, None
    ordem = []
    with open(cam, errors="replace") as fh:
        cols = [c.strip().lower() for c in fh.readline().rstrip("\n").split("\t")]
        ic = next((i for i, c in enumerate(cols)
                   if c in ("complexo", "complex", "pdb", "pdb_id", "id")), 0)
        for linha in fh:
            p = linha.rstrip("\n").split("\t")
            if len(p) > ic:
                pid = p[ic].strip().upper()
                if pid and pid not in ordem:
                    ordem.append(pid)
    return cam, ordem


def permutacoes(nomes):
    """(rotulo, lista permutada) para cada convencao de gerador."""
    n = len(nomes)
    out = []

    p = np.random.default_rng(SEMENTE).permutation(n)
    out.append(("default_rng.permutation(n)", [nomes[i] for i in p]))

    p = np.random.RandomState(SEMENTE).permutation(n)
    out.append(("RandomState.permutation(n)", [nomes[i] for i in p]))

    a = list(nomes)
    np.random.RandomState(SEMENTE).shuffle(a)
    out.append(("RandomState.shuffle(lista)", a))

    a = list(nomes)
    random.Random(SEMENTE).shuffle(a)
    out.append(("random.Random.shuffle", a))

    a = list(np.random.default_rng(SEMENTE).permutation(np.array(nomes, dtype=object)))
    out.append(("default_rng.permutation(lista)", [str(x) for x in a]))

    return out


def main():
    barra("V101 - recuperar a convencao, depois rodar a grade")
    casos = {}
    for a in sorted(glob.glob(os.path.join(BENCH, "*_descritores.tsv"))):
        casos[os.path.basename(a).split("_")[0].upper()] = le_tsv(a)
    if not casos:
        sys.exit("nenhum TSV em %s." % BENCH)
    print("  %d complexos medidos" % len(casos))

    cam, ordem = ordem_validacao()
    if ordem is None:
        print("  validacao_filtro.tsv nao encontrada.")
        ordem = []
    else:
        print("  ordem de %s: %d identificadores" % (cam, len(ordem)))

    disp = set(casos)
    ordens = []
    if ordem:
        ordens.append(("tabela, exclusoes ANTES",
                       [c for c in ordem if c in disp and c not in EXCLUIDOS]))
        ordens.append(("tabela, exclusoes DEPOIS",
                       [c for c in ordem if c in disp or c in EXCLUIDOS]))
    ordens.append(("alfabetica", sorted(disp)))

    barra("A.  QUAL CONVENCAO REPRODUZ A METADE RETIDA DOCUMENTADA?")
    print("  alvo: %s\n" % ", ".join(sorted(RETIDO_REGISTRO)))
    achou = None
    for rot_o, base in ordens:
        base = [c for c in base if c in disp]
        if len(base) % 2 or len(base) < 4:
            print("  %-26s  %d complexos, impar ou curto demais" % (rot_o, len(base)))
            continue
        for rot_g, perm in permutacoes(base):
            metade = len(perm) // 2
            for rot_m, ret in (("1a metade = retida", set(perm[:metade])),
                               ("2a metade = retida", set(perm[metade:]))):
                comum = len(ret & RETIDO_REGISTRO)
                bate = (ret == RETIDO_REGISTRO)
                if bate or comum >= 11:
                    print("  %-26s %-30s %-20s  %2d/13%s"
                          % (rot_o, rot_g, rot_m, comum, "   <== BATE" if bate else ""))
                if bate and achou is None:
                    achou = (rot_o, rot_g, rot_m, sorted(set(perm) - ret), sorted(ret))
    if achou is None:
        print("\n  Nenhuma das %d combinacoes reproduz a lista documentada."
              % (len(ordens) * 5 * 2))
        print("  A permutacao original nao e recuperavel do que esta no disco.")
        print("  Sigo com a lista nominal dos registros, que e o dado")
        print("  documentado, e o texto deve dizer que a divisao e fornecida")
        print("  por enumeracao e nao por re-execucao da permutacao.")
        ret = sorted(RETIDO_REGISTRO & disp)
        cal = sorted(disp - RETIDO_REGISTRO - EXCLUIDOS)
        # a calibracao sao os outros 13 do split de 26
        cal = [c for c in cal if c in disp][:13] if len(cal) > 13 else cal
        origem = "lista nominal dos registros"
    else:
        rot_o, rot_g, rot_m, cal, ret = achou
        origem = "%s | %s | %s" % (rot_o, rot_g, rot_m)
        print("\n  Convencao recuperada: %s" % origem)

    print("\n  calibracao (%d): %s" % (len(cal), ", ".join(cal)))
    print("  retida     (%d): %s" % (len(ret), ", ".join(ret)))
    cal = [c for c in cal if c in casos]
    ret = [c for c in ret if c in casos]

    # ------------------------------------------------- quantil + alvo ---
    def opera(col, metodo):
        oc = np.concatenate([casos[p][col][casos[p]["contato"]] for p in cal])
        oc = oc[np.isfinite(oc)]
        if oc.size == 0:
            return None
        lim = float(np.percentile(oc, QUANTIL, method=metodo))
        v = np.concatenate([casos[p][col] for p in ret])
        y = np.concatenate([casos[p]["contato"] for p in ret])
        m = np.isfinite(v)
        v, y = v[m], y[m]
        return (lim,
                100.0 * np.mean(v[y] <= lim) if y.sum() else np.nan,
                100.0 * np.mean(v[~y] > lim) if (~y).sum() else np.nan,
                int(y.sum()), int((~y).sum()))

    barra("B.  QUAL CONVENCAO DE QUANTIL DA 106 / 96.5 / 65.3 ?")
    melhor_q = None
    for metodo in ("lower", "linear", "higher", "nearest", "midpoint"):
        r = opera("cone_45_30", metodo)
        if r is None:
            continue
        lim, rt, rj, nc, nn = r
        bate = (abs(lim - LIMIAR_PUB) < 0.51 and abs(rt - RETEN_PUB) <= 0.3
                and abs(rj - REJEI_PUB) <= 0.3)
        print("  %-10s limiar %7.1f   retencao %5.1f%%   rejeicao %5.1f%%%s"
              % (metodo, lim, rt, rj, "   <== BATE" if bate else ""))
        if bate and melhor_q is None:
            melhor_q = metodo
    if melhor_q is None:
        print("\n  Nenhuma convencao de quantil reproduz o alvo com esta divisao.")
        print("  Sem reproduzir o ponto de operacao publicado, a grade abaixo")
        print("  nao sustenta afirmacao sobre a robustez DELE.  Ela so mostra")
        print("  como o ponto de operacao se move nesta re-execucao.")
        melhor_q = "linear"
    else:
        print("\n  Convencao de quantil: %s" % melhor_q)

    # ------------------------------------------------------------ grade -
    barra("C.  GRADE DO PONTO DE OPERACAO")
    base = opera("cone_45_30", melhor_q)
    print("  %-10s %-7s %9s %11s %11s %9s %9s"
          % ("abertura", "raio", "limiar", "retencao", "rejeicao", "d reten", "d rejei"))
    linhas = []
    for sa, ra in itertools.product(SEMIANGS, RAIOS):
        col = "cone_%d_%d" % (sa, ra)
        if col not in casos[next(iter(casos))]:
            continue
        r = opera(col, melhor_q)
        if r is None:
            continue
        lim, rt, rj, _, _ = r
        marca = "  <- usado no artigo" if (sa, ra) == (45, 30) else ""
        print("  %-10s %-7s %9.1f %10.1f%% %10.1f%% %+9.1f %+9.1f%s"
              % ("%d graus" % sa, "%d A" % ra, lim, rt, rj,
                 rt - base[1], rj - base[2], marca))
        linhas.append((sa, ra, lim, rt, rj, rt - base[1], rj - base[2]))

    p = os.path.join(BENCH, "V101_grade_operacao.tsv")
    with open(p, "w") as fh:
        fh.write("# divisao: %s\n# quantil: %s\n" % (origem, melhor_q))
        fh.write("semiangulo\traio\tlimiar\tretencao\trejeicao\td_reten\td_rejei\n")
        for l in linhas:
            fh.write("%d\t%d\t%.1f\t%.2f\t%.2f\t%.2f\t%.2f\n" % l)
    print("\n  %s" % p)

    barra("LEITURA")
    if linhas:
        print("  limiar de %.0f a %.0f atomos" % (min(l[2] for l in linhas),
                                                  max(l[2] for l in linhas)))
        print("  retencao de %.1f%% a %.1f%%" % (min(l[3] for l in linhas),
                                                 max(l[3] for l in linhas)))
        print("  rejeicao de %.1f%% a %.1f%%" % (min(l[4] for l in linhas),
                                                 max(l[4] for l in linhas)))
    print("\n  O limiar muda muito com a abertura, e tem de mudar: o cone de")
    print("  60 graus varre varias vezes mais volume que o de 30.  O que o")
    print("  parecerista quer saber nao e se o limiar e estavel, e sim se o")
    print("  PONTO DE OPERACAO -- retencao e rejeicao apos recalibrar -- se")
    print("  mantem.  E isso que as duas ultimas colunas mostram.")


if __name__ == "__main__":
    main()
