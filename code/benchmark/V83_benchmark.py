#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
V83 - O benchmark que o parecer pediu: o cone contra as alternativas,
      na mesma particula montada, com os mesmos residuos e os mesmos
      rotulos de contato.

Descritores, com a direcao fixada ANTES de rodar (contato = exposto):
  ocupacao no cone 45/30        contato = BAIXA   -> AUC sobre -valor
  grade 3x3: semiangulo {30,45,60} x raio {20,30,40} A
  cone 45/30 sem o proprio residuo e sem os vizinhos sequenciais +-1
  SASA na particula montada, sondas 1.4, 5.0 e 10.0 A
                                contato = ALTA    -> AUC sobre +valor
  HSE para fora (hemisferio no sentido CA->CB), raio 13 A
                                contato = BAIXA   -> AUC sobre -valor
  percentil radial (shell_pct da pegada publicada)
                                contato = ALTA    -> AUC sobre +valor

A SASA de sonda grande e EXATA, nao aproximada: um atomo so influencia
outro se as esferas expandidas pela sonda se tocam, ou seja ate
r_i + r_j + 2*sonda = 1.88 + 1.88 + 20 = 23.8 A para a sonda de 10 A.
Recortar a cadeia de referencia mais tudo num raio de 25 A preserva o
resultado e troca 400 mil atomos por algumas dezenas de milhares.

Rotulos de contato: coluna in_footprint das pegadas publicadas.
Le as montagens verificadas de /home/juan/V78_montagens.
Escreve so em /home/juan/V83_benchmark.  Nao toca em nada existente.

Uso:
    python3 V83_benchmark.py
    python3 V83_benchmark.py --sem-sonda-grande     # pula 5 e 10 A
    python3 V83_benchmark.py --so 6I2K,7XXA         # subconjunto
"""

import os
import sys
import math
import glob
import time
import numpy as np

sys.path.insert(0, "/home/juan")

try:
    from scipy.spatial import cKDTree
    from scipy.stats import rankdata
except ImportError:
    sys.exit("scipy necessario.   pip install scipy")

try:
    import freesasa
except ImportError:
    freesasa = None

MONTAGENS = "/home/juan/V78_montagens"
SAIDA     = "/home/juan/V83_benchmark"

DIRS_PEGADA = [
    "/home/juan/repo_a71_V69_20261125T130138_06f6e831/data/per_residue",
    "/home/juan/V49_pegadas",
    "/home/juan/abertos",
]

COSMIN_PUB = 0.7071          # constante que esta no filtro_externo.py
SEMIANGS   = [30.0, 45.0, 60.0]
RAIOS      = [20.0, 30.0, 40.0]
SONDAS     = [1.4, 5.0, 10.0]
MARGEM     = 25.0            # recorte exato para a maior sonda
HSE_RAIO   = 13.0

RAIO_PROTOR = {"C": 1.88, "N": 1.64, "O": 1.46, "S": 1.77, "P": 1.80}
RAIO_PADRAO = 1.80

SEM_SONDA_GRANDE = "--sem-sonda-grande" in sys.argv
SUBSET = None
for a in sys.argv:
    if a.startswith("--so"):
        if "=" in a:
            SUBSET = [x.strip().upper() for x in a.split("=", 1)[1].split(",")]
if SUBSET is None:
    for i, a in enumerate(sys.argv):
        if a == "--so" and i + 1 < len(sys.argv):
            SUBSET = [x.strip().upper() for x in sys.argv[i + 1].split(",")]

# direcao fixada antes de rodar: +1 se contato tem valor ALTO, -1 se BAIXO
DIRECAO = {"sasa": +1, "cone": -1, "hse": -1, "shell": +1}


def barra(t):
    print("\n" + "=" * 78)
    print(t)
    print("=" * 78)


def auc(score, rotulo):
    """Mann-Whitney com empates.  rotulo booleano: True = contato."""
    y = np.asarray(rotulo, dtype=bool)
    n1, n0 = int(y.sum()), int((~y).sum())
    if n1 == 0 or n0 == 0:
        return float("nan")
    r = rankdata(np.asarray(score, dtype=float))
    return (r[y].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0)


# ---------------------------------------------------------------- dados ---

def acha_pegada(pid):
    for d in DIRS_PEGADA:
        if not os.path.isdir(d):
            continue
        for pad in ("%s_*_pegada_V49*.tsv", "%s_*_pegada*.tsv"):
            ach = sorted(glob.glob(os.path.join(d, pad % pid)))
            if ach:
                return ach[0]
    return None


def le_pegada(caminho):
    """{residuo: (blocking, shell_pct, in_footprint)}"""
    out = {}
    with open(caminho, "r", errors="replace") as fh:
        cols = [c.strip().lower() for c in fh.readline().rstrip("\n").split("\t")]
        need = ("residue", "blocking", "shell_pct", "in_footprint")
        if any(c not in cols for c in need):
            return out
        ii = {c: cols.index(c) for c in need}
        for linha in fh:
            p = linha.rstrip("\n").split("\t")
            if len(p) <= max(ii.values()):
                continue
            try:
                r = int(float(p[ii["residue"]]))
                b = int(float(p[ii["blocking"]]))
                s = float(p[ii["shell_pct"]])
            except ValueError:
                continue
            f = p[ii["in_footprint"]].strip().lower() in ("yes", "y", "1", "true")
            out[r] = (b, s, f)
    return out


def carrega_npz(caminho):
    z = np.load(caminho, allow_pickle=False)
    d = {
        "xyz":  z["xyz"].astype(np.float64),
        "op":   z["op"],
        "elem": np.char.decode(z["elem"], "ascii"),
        "cad":  np.char.decode(z["cad"], "ascii"),
        "seq":  z["seq"],
        "nome": np.char.decode(z["nome"], "ascii"),
        "comp": np.char.decode(z["comp"], "ascii"),
        "mascara": z["mascara"],
        "cad_ref": bytes(z["cad_ref"]).decode(),
        "id_op": int(z["id_op"]),
        "composicao": bytes(z["composicao"]).decode(),
    }
    return d


# ----------------------------------------------------------- descritores --

def cone_conta(vert, nuvem, arv, centro, semiang, raio, excluir=None):
    u = vert - centro
    u = u / np.linalg.norm(u)
    cand = np.asarray(arv.query_ball_point(vert, raio + 1e-9), dtype=int)
    if cand.size == 0:
        return 0
    if excluir is not None and excluir.size:
        cand = cand[~np.isin(cand, excluir)]
        if cand.size == 0:
            return 0
    d = nuvem[cand] - vert
    dist = np.sqrt((d * d).sum(axis=1))
    viz = dist > 1e-6
    cos = np.zeros(len(d))
    cos[viz] = (d[viz] @ u) / dist[viz]
    lim = COSMIN_PUB if abs(semiang - 45.0) < 1e-9 else math.cos(math.radians(semiang))
    return int(np.count_nonzero(viz & (dist <= raio) & (cos >= lim)))


def raios_atomicos(elem):
    r = np.full(len(elem), RAIO_PADRAO)
    for e, v in RAIO_PROTOR.items():
        r[elem == e] = v
    return r


def sasa_por_residuo(xyz_sub, elem_sub, seq_sub, idx_ref_local, sonda):
    """SASA absoluta por residuo da cadeia de referencia, na particula."""
    if freesasa is None:
        return None
    par = freesasa.Parameters()
    par.setProbeRadius(sonda)
    res = freesasa.calcCoord(xyz_sub.flatten().tolist(),
                             raios_atomicos(elem_sub).tolist(), par)
    areas = np.array([res.atomArea(i) for i in range(len(elem_sub))])
    out = {}
    for i in idx_ref_local:
        out[int(seq_sub[i])] = out.get(int(seq_sub[i]), 0.0) + float(areas[i])
    return out


def hse_para_fora(pos_ca, vet_cb, todos_ca, arv_ca, raio=HSE_RAIO):
    """Numero de CA no hemisferio apontado pelo vetor CA->CB."""
    out = []
    for p, v in zip(pos_ca, vet_cb):
        if v is None or not np.isfinite(v).all():
            out.append(np.nan)
            continue
        cand = np.asarray(arv_ca.query_ball_point(p, raio), dtype=int)
        if cand.size == 0:
            out.append(0)
            continue
        d = todos_ca[cand] - p
        dist = np.sqrt((d * d).sum(axis=1))
        m = dist > 1e-6
        out.append(int(np.count_nonzero(m & ((d @ v) > 0))))
    return np.asarray(out, dtype=float)


def vetores_cb(mt, ref_res_idx):
    """CA->CB por residuo; para GLY usa CB virtual a partir de N, CA, C."""
    vets = []
    for atomos in ref_res_idx:
        nomes = {mt["nome"][i]: mt["xyz"][i] for i in atomos}
        ca = nomes.get("CA")
        if ca is None:
            vets.append(None)
            continue
        cb = nomes.get("CB")
        if cb is None:
            n, c = nomes.get("N"), nomes.get("C")
            if n is None or c is None:
                vets.append(None)
                continue
            b, cc = ca - n, c - ca
            a = np.cross(b, cc)
            cb = -0.58273431 * a + 0.56802827 * b - 0.54067466 * cc + ca
        v = cb - ca
        nv = np.linalg.norm(v)
        vets.append(v / nv if nv > 1e-6 else None)
    return vets


# ------------------------------------------------------------- por caso ---

def processa(pid, caminho_npz):
    t0 = time.time()
    mt = carrega_npz(caminho_npz)
    peg_arq = acha_pegada(pid)
    if peg_arq is None:
        print("  %-6s sem pegada" % pid)
        return None
    peg = le_pegada(peg_arq)
    if not peg:
        print("  %-6s pegada sem as colunas necessarias" % pid)
        return None

    masc = mt["mascara"]
    nuvem = mt["xyz"][masc]
    centro = nuvem.mean(axis=0)
    arv = cKDTree(nuvem)

    ref_op = max(mt["id_op"], 0)
    m_ref = (mt["op"] == ref_op) & (mt["cad"] == mt["cad_ref"])
    seq_ref = mt["seq"][m_ref]
    idx_ref_global = np.flatnonzero(m_ref)

    residuos = sorted(set(int(r) for r in seq_ref) & set(peg.keys()))
    if len(residuos) < 20:
        print("  %-6s so %d residuos em comum com a pegada" % (pid, len(residuos)))
        return None

    por_res = {r: [] for r in residuos}
    for i in idx_ref_global:
        r = int(mt["seq"][i])
        if r in por_res:
            por_res[r].append(i)

    pos_ca, ok = [], []
    for r in residuos:
        ca = [i for i in por_res[r] if mt["nome"][i] == "CA"]
        if ca:
            pos_ca.append(mt["xyz"][ca[0]])
            ok.append(r)
    residuos = ok
    pos_ca = np.asarray(pos_ca)

    # indices na NUVEM (nao no array completo) dos atomos de cada residuo
    mapa = -np.ones(len(mt["xyz"]), dtype=np.int64)
    mapa[np.flatnonzero(masc)] = np.arange(int(masc.sum()))

    dados = {"residuo": residuos}
    dados["blocking_pub"] = [peg[r][0] for r in residuos]
    dados["shell_pct"]    = [peg[r][1] for r in residuos]
    dados["contato"]      = [peg[r][2] for r in residuos]

    # ---- grade 3x3 do cone
    for sa in SEMIANGS:
        for ra in RAIOS:
            col = "cone_%.0f_%.0f" % (sa, ra)
            dados[col] = [cone_conta(p, nuvem, arv, centro, sa, ra)
                          for p in pos_ca]

    # ---- cone 45/30 sem o proprio residuo e sem os vizinhos sequenciais
    sem_viz = []
    for r, p in zip(residuos, pos_ca):
        fora = []
        for rr in (r - 1, r, r + 1):
            for i in por_res.get(rr, []):
                j = mapa[i]
                if j >= 0:
                    fora.append(j)
        sem_viz.append(cone_conta(p, nuvem, arv, centro, 45.0, 30.0,
                                  excluir=np.asarray(fora, dtype=int)))
    dados["cone_45_30_sem_vizinhos"] = sem_viz

    # ---- SASA na particula, recorte exato de 25 A
    sondas = [1.4] if SEM_SONDA_GRANDE else SONDAS
    if freesasa is None:
        print("  %-6s freesasa ausente, SASA pulada" % pid)
    else:
        xyz_ref = mt["xyz"][idx_ref_global]
        arv_ref = cKDTree(xyz_ref)
        dmin, _ = arv_ref.query(nuvem, k=1, distance_upper_bound=MARGEM)
        viz_masc = np.isfinite(dmin) & (dmin <= MARGEM)
        sub_idx = np.flatnonzero(viz_masc)
        xyz_sub = nuvem[sub_idx]
        elem_sub = mt["elem"][masc][sub_idx]
        seq_sub = mt["seq"][masc][sub_idx]
        cad_sub = mt["cad"][masc][sub_idx]
        op_sub = mt["op"][masc][sub_idx]
        local_ref = np.flatnonzero((cad_sub == mt["cad_ref"]) & (op_sub == ref_op))
        for s in sondas:
            m = sasa_por_residuo(xyz_sub, elem_sub, seq_sub, local_ref, s)
            dados["sasa_%.1f" % s] = [m.get(r, np.nan) if m else np.nan
                                      for r in residuos]
        dados["_n_sub"] = len(sub_idx)

    # ---- HSE
    m_ca_tudo = mt["nome"][masc] == "CA"
    todos_ca = nuvem[m_ca_tudo]
    arv_ca = cKDTree(todos_ca)
    vets = vetores_cb(mt, [por_res[r] for r in residuos])
    dados["hse_fora"] = hse_para_fora(pos_ca, vets, todos_ca, arv_ca)

    dt = time.time() - t0
    n_cont = int(np.sum(dados["contato"]))
    print("  %-6s %4d residuos  %3d contatos  %8d atomos  %6.0f s  (%s)"
          % (pid, len(residuos), n_cont, len(nuvem), dt, mt["composicao"][:28]))
    return dados


def salva_tsv(pid, dados):
    os.makedirs(SAIDA, exist_ok=True)
    cols = [c for c in dados if not c.startswith("_")]
    caminho = os.path.join(SAIDA, "%s_descritores.tsv" % pid)
    n = len(dados["residuo"])
    with open(caminho, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for i in range(n):
            linha = []
            for c in cols:
                v = dados[c][i]
                if isinstance(v, (bool, np.bool_)):
                    linha.append("yes" if v else "no")
                elif isinstance(v, float) and not np.isfinite(v):
                    linha.append("NA")
                elif isinstance(v, float):
                    linha.append("%.3f" % v)
                else:
                    linha.append(str(v))
            fh.write("\t".join(linha) + "\n")
    return caminho


def direcao_de(col):
    if col.startswith("sasa"):
        return DIRECAO["sasa"]
    if col.startswith("cone") or col == "blocking_pub":
        return DIRECAO["cone"]
    if col.startswith("hse"):
        return DIRECAO["hse"]
    if col.startswith("shell"):
        return DIRECAO["shell"]
    return +1


def main():
    barra("V83 - benchmark do cone contra as alternativas")
    arqs = sorted(glob.glob(os.path.join(MONTAGENS, "*_montagem.npz")))
    if SUBSET:
        arqs = [a for a in arqs
                if os.path.basename(a).split("_")[0].upper() in SUBSET]
    if not arqs:
        sys.exit("nenhuma montagem em %s.  rode o V81 primeiro." % MONTAGENS)
    print("  %d montagens   sondas %s   grade %dx%d"
          % (len(arqs), "1.4 so" if SEM_SONDA_GRANDE else "1.4/5/10",
             len(SEMIANGS), len(RAIOS)))
    if freesasa is None:
        print("  AVISO: freesasa ausente, a SASA nao sera calculada.")
    print("  recorte de %.0f A: exato para sonda de 10 A (1.88+1.88+20=23.8)\n"
          % MARGEM)

    todos = {}
    for k, a in enumerate(arqs, 1):
        pid = os.path.basename(a).split("_")[0].upper()
        d = processa(pid, a)
        if d is None:
            continue
        todos[pid] = d
        salva_tsv(pid, d)
        if k == 1:
            print("         >>> ao ritmo do primeiro, os %d levam ~%.0f min"
                  % (len(arqs), (len(arqs) - 1) * 0 + 0))

    if not todos:
        sys.exit("nada processado.")

    # --------------------------------------------------------- AUC ------
    barra("AUC  -  por descritor, agrupando os complexos")
    cols = [c for c in next(iter(todos.values()))
            if c not in ("residuo", "contato") and not c.startswith("_")]

    pool = {c: [] for c in cols}
    rot = []
    for pid, d in todos.items():
        rot.extend(d["contato"])
        for c in cols:
            pool[c].extend(d[c])
    rot = np.asarray(rot, dtype=bool)

    print("  %d residuos, %d contatos, %d complexos\n"
          % (len(rot), int(rot.sum()), len(todos)))
    print("  %-28s %8s %10s" % ("descritor", "AUC", "direcao"))
    linhas = []
    for c in cols:
        v = np.asarray(pool[c], dtype=float)
        m = np.isfinite(v)
        if m.sum() < 50 or rot[m].sum() < 5:
            continue
        s = direcao_de(c)
        a = auc(s * v[m], rot[m])
        linhas.append((a, c, s))
    for a, c, s in sorted(linhas, reverse=True):
        print("  %-28s %8.4f %10s" % (c, a, "alto" if s > 0 else "baixo"))

    # por complexo, so para o cone publicado e a maior sonda
    chave = [c for c in ("cone_45_30", "sasa_10.0", "sasa_5.0", "shell_pct",
                         "hse_fora", "cone_45_30_sem_vizinhos") if c in cols]
    barra("AUC por complexo")
    print("  %-6s %6s " % ("pdb", "cont") + " ".join("%12s" % c for c in chave))
    for pid in sorted(todos):
        d = todos[pid]
        y = np.asarray(d["contato"], dtype=bool)
        if y.sum() == 0 or (~y).sum() == 0:
            continue
        saida = []
        for c in chave:
            v = np.asarray(d[c], dtype=float)
            m = np.isfinite(v)
            saida.append("%12.3f" % auc(direcao_de(c) * v[m], y[m])
                         if m.sum() > 20 and y[m].sum() > 0 else "%12s" % "-")
        print("  %-6s %6d " % (pid, int(y.sum())) + " ".join(saida))

    barra("saida em %s" % SAIDA)
    print("  um TSV por complexo, com todos os descritores por residuo.")
    print("\n  Leia a tabela de AUC acima antes de escrever qualquer coisa:")
    print("  se a SASA de sonda grande empatar ou passar o cone, a")
    print("  contribuicao do paper muda de 'descritor novo' para")
    print("  'calibracao de um descritor', e o texto tem que dizer isso.")


if __name__ == "__main__":
    main()
