#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
V86 - Profundidade de residuo na particula montada, o quarto descritor
      que faltava no benchmark.  Mais a busca pela origem do 34.612.

DEFINICAO, para poder ser nomeada corretamente no texto:
  profundidade de um residuo = media, sobre os atomos do residuo, da
  distancia ao ponto mais proximo de solvente CONTIGUO AO EXTERIOR da
  particula.  E a definicao de Chakravarty e Varshney (distancia a agua
  do seio do solvente), calculada em grade de 1 A em vez de MSMS.

A armadilha que isso evita: como a vizinhanca e recortada em 25 A ao
redor da cadeia de referencia, a face interna da caixa cai em espaco
vazio ARTIFICIAL -- vazio por truncamento, nao por ser solvente.
Preencher a partir de toda a borda inventaria solvente dentro do
capsideo.  Entao o preenchimento parte so dos pontos cujo raio ao
centroide da particula excede o raio maximo dos atomos, que sao
exterior inequivoco, e avanca para dentro pelo espaco livre.

Direcao fixada antes de rodar: contato = exposto = profundidade BAIXA,
entao a AUC e calculada sobre -profundidade.

Le as montagens de /home/juan/V78_montagens e os rotulos dos TSV de
/home/juan/V83_benchmark.  Escreve arquivos NOVOS; nao sobrescreve os
do V83.

Uso:
    python3 V86_profundidade.py
    python3 V86_profundidade.py --so-busca-34612
    python3 V86_profundidade.py --espaco 1.5      # grade mais grossa
"""

import os
import re
import sys
import glob
import time
import numpy as np

try:
    from scipy.spatial import cKDTree
    from scipy.stats import rankdata
    from scipy import ndimage
except ImportError:
    sys.exit("scipy necessario.")

MONTAGENS = "/home/juan/V78_montagens"
BENCH     = "/home/juan/V83_benchmark"
SONDA     = 1.4
MARGEM    = 25.0
SEMENTE   = 20261125
NBOOT     = 5000

ESPACO = 1.0
for i, a in enumerate(sys.argv):
    if a == "--espaco" and i + 1 < len(sys.argv):
        ESPACO = float(sys.argv[i + 1])
SO_BUSCA = "--so-busca-34612" in sys.argv

RAIO_PROTOR = {"C": 1.88, "N": 1.64, "O": 1.46, "S": 1.77, "P": 1.80}
RAIO_PADRAO = 1.80


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


def raios(elem):
    r = np.full(len(elem), RAIO_PADRAO)
    for e, v in RAIO_PROTOR.items():
        r[elem == e] = v
    return r


def le_tsv(caminho):
    with open(caminho, "r") as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        dados = {c: [] for c in cols}
        for linha in fh:
            p = linha.rstrip("\n").split("\t")
            if len(p) != len(cols):
                continue
            for c, v in zip(cols, p):
                dados[c].append(v)
    out = {}
    for c, vs in dados.items():
        if c == "contato":
            out[c] = np.array([v.strip().lower() == "yes" for v in vs], dtype=bool)
        else:
            out[c] = np.array([np.nan if v == "NA" else float(v) for v in vs])
    return out


# --------------------------------------------------------- profundidade ---

def profundidade(xyz_sub, rad_sub, centro_part, rmax_part, residuos_atomos,
                 espaco=ESPACO, sonda=SONDA):
    """Devolve {residuo: (profundidade media dos atomos, profundidade do CA)}."""
    mn = xyz_sub.min(axis=0) - 6.0
    mx = xyz_sub.max(axis=0) + 6.0
    eixos = [np.arange(mn[k], mx[k] + espaco, espaco) for k in range(3)]
    forma = [len(e) for e in eixos]
    if np.prod(forma) > 40e6:
        return None, "grade grande demais (%s)" % "x".join(map(str, forma))

    malha = np.stack(np.meshgrid(*eixos, indexing="ij"), axis=-1)
    pts = malha.reshape(-1, 3)

    arv = cKDTree(xyz_sub)
    d, i = arv.query(pts, k=1, workers=-1)
    livre = d > (rad_sub[i] + sonda)
    vol = livre.reshape(forma)

    # componentes conexas do espaco livre
    lab, n = ndimage.label(vol)
    if n == 0:
        return None, "nenhum espaco livre na grade"

    # semente: so pontos inequivocamente FORA da particula
    raio_pt = np.linalg.norm(pts - centro_part, axis=1).reshape(forma)
    semente = vol & (raio_pt > rmax_part + 2.0)
    if not semente.any():
        return None, "a caixa nao alcanca o exterior da particula"
    rot_ext = np.unique(lab[semente])
    rot_ext = rot_ext[rot_ext > 0]
    if rot_ext.size == 0:
        return None, "nenhuma componente livre toca o exterior"

    bulk = np.isin(lab, rot_ext)
    pts_bulk = pts[bulk.reshape(-1)]
    if len(pts_bulk) < 100:
        return None, "solvente externo com %d pontos so" % len(pts_bulk)

    arv_bulk = cKDTree(pts_bulk)
    out = {}
    for r, (idx_atomos, idx_ca) in residuos_atomos.items():
        if not idx_atomos:
            continue
        dd, _ = arv_bulk.query(xyz_sub[idx_atomos], k=1, workers=-1)
        dca = float(dd[idx_atomos.index(idx_ca)]) if idx_ca in idx_atomos else np.nan
        out[r] = (float(np.mean(dd)), dca)
    return out, None


def processa(pid):
    npz = os.path.join(MONTAGENS, "%s_montagem.npz" % pid)
    tsv = os.path.join(BENCH, "%s_descritores.tsv" % pid)
    if not (os.path.exists(npz) and os.path.exists(tsv)):
        return None

    z = np.load(npz, allow_pickle=False)
    xyz = z["xyz"].astype(np.float64)
    masc = z["mascara"]
    cad = np.char.decode(z["cad"], "ascii")
    nome = np.char.decode(z["nome"], "ascii")
    elem = np.char.decode(z["elem"], "ascii")
    seq = z["seq"]
    op = z["op"]
    cad_ref = bytes(z["cad_ref"]).decode()
    ref_op = max(int(z["id_op"]), 0)

    nuvem = xyz[masc]
    centro = nuvem.mean(axis=0)
    rmax = float(np.linalg.norm(nuvem - centro, axis=1).max())

    m_ref_global = (op == ref_op) & (cad == cad_ref)
    xyz_ref = xyz[m_ref_global]

    # recorte de 25 A, igual ao do V83
    arv_ref = cKDTree(xyz_ref)
    dmin, _ = arv_ref.query(nuvem, k=1, distance_upper_bound=MARGEM, workers=-1)
    sub = np.flatnonzero(np.isfinite(dmin) & (dmin <= MARGEM))
    xyz_sub = nuvem[sub]
    rad_sub = raios(elem[masc][sub])
    seq_sub = seq[masc][sub]
    cad_sub = cad[masc][sub]
    op_sub = op[masc][sub]
    nome_sub = nome[masc][sub]

    e_ref = (cad_sub == cad_ref) & (op_sub == ref_op)
    residuos_atomos = {}
    for j in np.flatnonzero(e_ref):
        r = int(seq_sub[j])
        a, ca = residuos_atomos.setdefault(r, ([], -1))
        a.append(int(j))
        if nome_sub[j] == "CA":
            residuos_atomos[r] = (a, int(j))

    prof, err = profundidade(xyz_sub, rad_sub, centro, rmax, residuos_atomos)
    if prof is None:
        print("  %-6s %s" % (pid, err))
        return None

    d = le_tsv(tsv)
    res = d["residuo"].astype(int)
    pm = np.array([prof.get(int(r), (np.nan, np.nan))[0] for r in res])
    pc = np.array([prof.get(int(r), (np.nan, np.nan))[1] for r in res])

    saida = os.path.join(BENCH, "%s_profundidade.tsv" % pid)
    with open(saida, "w") as fh:
        fh.write("residuo\tprof_media\tprof_CA\tcontato\tcone_45_30\n")
        for k in range(len(res)):
            fh.write("%d\t%s\t%s\t%s\t%s\n"
                     % (res[k],
                        "NA" if not np.isfinite(pm[k]) else "%.3f" % pm[k],
                        "NA" if not np.isfinite(pc[k]) else "%.3f" % pc[k],
                        "yes" if d["contato"][k] else "no",
                        "%.0f" % d["cone_45_30"][k]))
    return {"pid": pid, "residuo": res, "prof_media": pm, "prof_CA": pc,
            "contato": d["contato"], "cone_45_30": d["cone_45_30"]}


# ------------------------------------------------------------- 34612 -----

def busca_34612():
    barra("BUSCA  -  de onde vem o 34.612")
    alvos = ("34612", "34,612", "34 612", "34.612")
    ign = ("/.git/", "site-packages", "/envs/", "miniconda", "anaconda",
           "/.cache/", "/node_modules/")
    exts = (".txt", ".md", ".tsv", ".csv", ".log", ".py", ".json", ".yaml", ".yml")

    print("  1. ocorrencias literais do numero em arquivos de texto")
    achou = 0
    for caminho in glob.iglob("/home/juan/**/*", recursive=True):
        if not caminho.endswith(exts) or any(x in caminho for x in ign):
            continue
        try:
            if os.path.getsize(caminho) > 80e6:
                continue
            with open(caminho, "r", errors="replace") as fh:
                for n, linha in enumerate(fh, 1):
                    if any(a in linha for a in alvos):
                        print("     %s:%d" % (caminho, n))
                        print("        %s" % linha.strip()[:160])
                        achou += 1
                        if achou > 40:
                            print("     (mais de 40, parando)")
                            return
        except (OSError, UnicodeDecodeError):
            continue
    if achou == 0:
        print("     nenhuma.  O numero nao aparece escrito em nenhum arquivo")
        print("     de texto sob /home/juan.")

    print("\n  2. tabelas com contagem de linhas entre 34.000 e 35.200")
    perto = []
    for ext in ("tsv", "csv"):
        for caminho in glob.iglob("/home/juan/**/*.%s" % ext, recursive=True):
            if any(x in caminho for x in ign):
                continue
            try:
                if not (0.5e6 < os.path.getsize(caminho) < 400e6):
                    continue
                with open(caminho, "rb") as fh:
                    n = sum(1 for _ in fh) - 1
            except OSError:
                continue
            if 34000 <= n <= 35200:
                perto.append((n, caminho))
    if perto:
        for n, c in sorted(perto):
            print("     %7d linhas  %s" % (n, c))
    else:
        print("     nenhuma.")

    print("\n  Se as duas buscas vierem vazias, o numero nao e recuperavel do")
    print("  disco e a frase precisa ser reescrita em torno da populacao que")
    print("  a figura mostra.  Isso e edicao do manuscrito, nao faco sozinho.")


def main():
    barra("V86 - profundidade de residuo, e a origem do 34.612")
    if SO_BUSCA:
        busca_34612()
        return

    print("  grade de %.1f A   sonda %.1f A   recorte %.0f A" % (ESPACO, SONDA, MARGEM))
    print("  preenchimento semeado so no exterior da particula\n")

    arqs = sorted(glob.glob(os.path.join(MONTAGENS, "*_montagem.npz")))
    casos = []
    for k, a in enumerate(arqs, 1):
        pid = os.path.basename(a).split("_")[0].upper()
        t0 = time.time()
        c = processa(pid)
        if c is None:
            continue
        casos.append(c)
        y = c["contato"]
        m = np.isfinite(c["prof_media"])
        a_prof = auc(-c["prof_media"][m], y[m])
        a_cone = auc(-c["cone_45_30"][m], y[m])
        print("  %-6s %4d residuos  %3d contatos   AUC prof %.3f   cone %.3f   %4.0f s"
              % (pid, int(m.sum()), int(y.sum()), a_prof, a_cone, time.time() - t0))

    if not casos:
        sys.exit("nada processado.")

    barra("AUC agrupada e diferenca pareada")
    pm = np.concatenate([c["prof_media"] for c in casos])
    pc = np.concatenate([c["prof_CA"] for c in casos])
    cn = np.concatenate([c["cone_45_30"] for c in casos])
    yy = np.concatenate([c["contato"] for c in casos])
    for rot, v in (("profundidade media dos atomos", pm),
                   ("profundidade do CA", pc),
                   ("cone 45/30 (referencia)", cn)):
        m = np.isfinite(v)
        print("  %-32s AUC %.4f   (n=%d)" % (rot, auc(-v[m], yy[m]), int(m.sum())))

    rng = np.random.default_rng(SEMENTE)
    difs, vence = [], 0
    for c in casos:
        y = c["contato"]
        m = np.isfinite(c["prof_media"])
        if y[m].sum() == 0 or (~y[m]).sum() == 0:
            continue
        a_cone = auc(-c["cone_45_30"][m], y[m])
        a_prof = auc(-c["prof_media"][m], y[m])
        difs.append(a_cone - a_prof)
        vence += int(a_cone > a_prof)
    difs = np.asarray(difs)
    boot = np.array([difs[rng.integers(0, len(difs), len(difs))].mean()
                     for _ in range(NBOOT)])
    print("\n  cone menos profundidade, pareado por complexo")
    print("    diferenca media %+.4f   IC95 [%+.4f, %+.4f]   vence em %d/%d"
          % (difs.mean(), np.percentile(boot, 2.5), np.percentile(boot, 97.5),
             vence, len(difs)))
    print("\n  semente %d, %d reamostragens, escolhida para esta analise" % (SEMENTE, NBOOT))

    busca_34612()
    barra("saida")
    print("  um %s/<pdb>_profundidade.tsv por complexo." % BENCH)
    print("  os TSV do V83 nao foram tocados.")


if __name__ == "__main__":
    main()
